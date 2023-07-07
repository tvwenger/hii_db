"""
wise.py

Utilities for adding WISE Catalog information to the database.

Copyright(C) 2020-2021 by
Trey V. Wenger; tvwenger@gmail.com

GNU General Public License v3 (GNU GPLv3)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

2020-04-01 Trey V. Wenger
2021-09-30 Trey V. Wenger reorganization
"""

import sqlite3
import numpy as np
from astropy.coordinates import SkyCoord


def gen_catalog(db, wisefile):
    """
    Read WISE Catalog CSV file and populate WISE Catalog table.

    Inputs:
        db :: string
            Database filename
        wisefile :: string
            WISE Catalog CSV filename

    Returns: Nothing
    """
    print("Generating WISE Catalog table...")

    # Read the WISE catalog
    wise = np.genfromtxt(
        wisefile,
        delimiter=",",
        dtype=None,
        names=True,
        autostrip=True,
        encoding="utf-8",
    )
    wise_coords = SkyCoord(wise["GLong"], wise["GLat"], frame="galactic", unit="deg")
    wise_ra = wise_coords.fk5.ra.deg
    wise_dec = wise_coords.fk5.dec.deg

    # Populate Catalog
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        rows = [
            (
                w["GName"],
                w["Name"] if w["Name"] != "" else None,
                w["HII_Name"] if w["HII_Name"] != "" else None,
                w["Catalog"],
                ra,
                dec,
                w["GLong"],
                w["GLat"],
                w["Size"],
                w["KDAR"].replace("(", "").replace(")", "")
                if w["KDAR"] != ""
                else None,
                w["Distance_Method"] if w["Distance_Method"] != "" else None,
                w["Distance_Author"] if w["Distance_Author"] != "" else None,
            )
            for w, ra, dec in zip(wise, wise_ra, wise_dec)
        ]
        cur.executemany(
            """
        INSERT INTO Catalog
        (gname, alias, hii_name, catalog, ra, dec, glong, glat,
        radius, kdar, dist_method, dist_author)
        VALUES
        (?,?,?,?,?,?,?,?,?,?,?,?)
        """,
            rows,
        )
    print("Done!")
    print()


def gen_groups(db, wisefile):
    """
    Read WISE Catalog CSV file and populate Groups table, then match
    WISE Catalog sources to Groups.

    Inputs:
        db :: string
            Database filename
        wisefile :: string
            WISE Catalog CSV filename

    Returns: Nothing
    """
    print("Generating Groups table...")

    # Read the WISE catalog
    wise = np.genfromtxt(
        wisefile,
        delimiter=",",
        dtype=None,
        names=True,
        autostrip=True,
        encoding="utf-8",
    )

    # Get group names/info
    all_groups = []
    group_info = []
    group_members = []
    for w in wise:
        if w["Group"] == "":
            continue
        group = w["Group"].strip()

        # check if group is already found
        if group in all_groups:
            # check that group data is the same
            idx = all_groups.index(group)
            if group_info[idx][0] == "":
                group_info[idx][0] = w["Group_VLSR"]
            elif w["Group_VLSR"] != "" and group_info[idx][0] != w["Group_VLSR"]:
                print("Group {0} different VLSR".format(group))
                continue

            if group_info[idx][1] == "":
                group_info[idx][1] = w["Group_e_VLSR"]
            elif w["Group_e_VLSR"] != "" and group_info[idx][1] != w["Group_e_VLSR"]:
                print("Group {0} different e_VLSR".format(group))
                continue

            if group_info[idx][2] == "":
                group_info[idx][2] = w["Group_KDAR"]
            elif w["Group_KDAR"] != "" and group_info[idx][2] != w["Group_KDAR"]:
                print("Group {0} different KDAR".format(group))
                continue

            # add member
            group_members[idx].append(w["GName"])
        else:
            # add group
            all_groups.append(group)
            group_info.append([w["Group_VLSR"], w["Group_e_VLSR"], w["Group_KDAR"]])
            group_members.append([w["GName"]])

    # Populate Groups
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        rows = [
            (
                group,
                data[0] if data[0] != "" else None,
                data[1] if data[1] != "" else None,
                data[2].replace("(", "").replace(")", "") if data[2] != "" else None,
            )
            for group, data in zip(all_groups, group_info)
        ]
        cur.executemany(
            """
        INSERT INTO Groups
        (name, vlsr, e_vlsr, kdar)
        VALUES
        (?,?,?,?)
        """,
            rows,
        )
    print("Done!")
    print()

    # Match Catalog to Groups
    print("Matching Catalog to Groups...")
    rows = []
    for group, members in zip(all_groups, group_members):
        # get Group ID
        with sqlite3.connect(db) as conn:
            cur = conn.cursor()
            cur.execute("PRAGMA foreign_keys = ON")
            cur.execute("SELECT id FROM Groups WHERE name=?", [group])
            group_id = cur.fetchall()[0][0]
        for member in members:
            with sqlite3.connect(db) as conn:
                # Get Catalog ID
                cur = conn.cursor()
                cur.execute("PRAGMA foreign_keys = ON")
                cur.execute("SELECT id FROM Catalog WHERE gname=?", [member])
                catalog_id = cur.fetchall()[0][0]
            rows.append((catalog_id, group_id))

    # Insert into CatalogGroups
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        cur.executemany(
            """
        INSERT INTO CatalogGroups
        (catalog_id, group_id)
        VALUES
        (?,?)
        """,
            rows,
        )
    print("Done!")
    print()


def add_detections(db, wisefile):
    """
    Read WISE Catalog and add detections to database.

    Inputs:
        db :: string
            Database filename
        wisefile :: string
            WISE Catalog CSV filename

    Returns: Nothing
    """
    print("Adding WISE Detections...")

    # Read the WISE Catalog
    wise = np.genfromtxt(
        wisefile,
        delimiter=",",
        dtype=None,
        names=True,
        autostrip=True,
        encoding="utf-8",
    )

    # Add WISE Detections to table, split multi-component sources
    # into two rows in Detections table
    rows = []
    for wiseidx, w in enumerate(wise):
        # skip non-detections
        if w["VLSR"] == "":
            continue
        if w["VLSR"] == "0.0" and w["FWHM"] == "":
            continue

        # Get coordinate
        coord = SkyCoord(
            w["GLong_Observed"], w["GLat_Observed"], frame="galactic", unit="deg"
        )
        ra = coord.fk5.ra.deg
        dec = coord.fk5.dec.deg

        # Split multi-component sources, handle missing data
        ncomp = len(w["VLSR"].split(";"))
        vlsrs = np.array([float(f) for f in w["VLSR"].split(";")])
        if len(w["e_VLSR"]) > 0:
            e_vlsrs = np.array([float(f) for f in w["e_VLSR"].split(";")])
        else:
            e_vlsrs = np.array([np.nan] * ncomp)

        fwhms = np.array([float(f) for f in w["FWHM"].split(";")])
        if len(w["e_FWHM"]) > 0:
            e_fwhms = np.array([float(f) for f in w["e_FWHM"].split(";")])
        else:
            e_fwhms = np.array([np.nan] * ncomp)

        if len(w["T_L"]) > 0:
            tls = np.array([float(f) for f in w["T_L"].split(";")])
        else:
            tls = np.array([np.nan] * ncomp)

        if len(w["e_T_L"]) > 0:
            e_tls = np.array([float(f) for f in w["e_T_L"].split(";")])
        else:
            e_tls = np.array([np.nan] * ncomp)

        # Fix Caswell & Haynes sources
        if "Caswell" in w["Author"]:
            tls = tls * 1000.0  # K to mK
            e_tls = np.array([10.0 for _ in e_tls])  # mK
            e_vlsrs = np.array([2.5 for _ in e_vlsrs])  # km/s
            e_fwhms = np.array([1.5 for _ in e_fwhms])  # km/s

        # Confirm that each parameter has the same number of components
        if np.any(
            np.array([len(e_vlsrs), len(fwhms), len(e_fwhms), len(tls), len(e_tls)])
            != ncomp
        ):
            print("PROBLEM WITH {0} COMPONENTS".format(w["GName"]))
            print(vlsrs)
            print(e_vlsrs)
            print(fwhms)
            print(e_fwhms)
            print(tls)
            print(e_tls)
            continue

        # Set other properties
        if "alpha" in w["Lines"]:
            line_unit = "optical"
            cont_unit = "optical"
            line_freq = 456700000.0
        else:
            line_unit = "mK"
            cont_unit = "mK"
            line_freq = 30000.0/w["Wavelength"]
        resolution = np.nan
        if w["Resolution"] != "":
            resolution = float(w["Resolution"].replace("~", ""))  # arcmin
        beam_area = (
            np.pi * (resolution * 60.0) ** 2.0 / (4.0 * np.log(2.0))
        )  # sq. arcsec
        area = np.pi * (w["HRDS_Size"]) ** 2.0 / (4.0 * np.log(2.0))  # sq. arcsec
        area_unit = "arcsec2"
        if area == 0.0:
            area = None
            area_unit = None
        cont_freq = np.nan
        cont = np.nan
        e_cont = np.nan
        cont_unit = np.nan
        if w["HRDS_Flux"] > 0.0:
            cont_freq = 9000.0
            cont = w["HRDS_Flux"]
            e_cont = w["HRDS_e_Flux"]
            cont_unit = "mJy"

        # Loop over components, sorted by tl, and populate table
        comps = ["a", "b", "c", "d", "e"]
        sortind = np.argsort(tls)
        for vlsr, e_vlsr, fwhm, e_fwhm, tl, e_tl, comp in zip(
            vlsrs[sortind],
            e_vlsrs[sortind],
            fwhms[sortind],
            e_fwhms[sortind],
            tls[sortind],
            e_tls[sortind],
            comps,
        ):
            if ncomp == 1:
                comp = None
            row = (
                w["GName"],
                ra,
                dec,
                w["GLong_Observed"],
                w["GLat_Observed"],
                line_freq,
                comp,
                tl,
                e_tl,
                line_unit,
                vlsr,
                e_vlsr,
                fwhm,
                e_fwhm,
                cont_freq,
                cont,
                e_cont,
                cont_unit,
                area,
                area_unit,
                w["T_e"],
                w["Lines"],
                beam_area,
                w["Telescope"],
                w["Author"],
                "WISE Catalog",
            )
            rows.append(row)

    # Populate detections table
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        cur.executemany(
            """
        INSERT INTO Detections
        (name, ra, dec, glong, glat, line_freq, component,
        line, e_line, line_unit, vlsr, e_vlsr, fwhm, e_fwhm,
        cont_freq, cont, e_cont, cont_unit, area, area_unit, te,
        lines, beam_area, telescope, author, source) VALUES
        (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
            rows,
        )
    print("Done!")
    print()

    # Match WISE detections to Catalog
    print("Matching WISE Detections to WISE Catalog...")
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")

        # Get WISE catalog gnames and positions
        cur.execute("SELECT id, gname, ra, dec, catalog FROM Catalog")
        wisecat = np.array(
            cur.fetchall(),
            dtype=[
                ("id", "i"),
                ("gname", "U100"),
                ("ra", "f8"),
                ("dec", "f8"),
                ("catalog", "U1"),
            ],
        )
        cat_coords = SkyCoord(wisecat["ra"], wisecat["dec"], frame="fk5", unit="deg")

        # Get the WISE detection gnames
        cur.execute(
            """
        SELECT id, name, ra, dec FROM Detections
        WHERE source="WISE Catalog"
        """
        )
        wisedet = np.array(
            cur.fetchall(),
            dtype=[("id", "i"), ("name", "U100"), ("ra", "f8"), ("dec", "f8")],
        )
        det_coords = SkyCoord(wisedet["ra"], wisedet["dec"], frame="fk5", unit="deg")

        # Match and calculate separations
        matches = np.array(
            [np.where(wisecat["gname"] == gname)[0][0] for gname in wisedet["name"]]
        )
        seps = [
            coord.separation(cat_coords[match]).arcsec
            for coord, match in zip(det_coords, matches)
        ]

        # Check that Catalog entires have 'catalog' == 'K'
        bad = wisecat[matches]["catalog"] != "K"
        if np.sum(bad) > 0:
            print("The following sources have catalog != K:")
            print(wisecat[matches][bad]["gname"])

        # Populate Catalog-Detections
        rows = [
            (int(wisecat["id"][match]), int(det["id"]), float(sep))
            for det, match, sep in zip(wisedet, matches, seps)
        ]
        cur.executemany(
            """
        INSERT INTO CatalogDetections
        (catalog_id, detection_id, separation) VALUES (?, ?, ?)
        """,
            rows,
        )
    print("Done!")
    print()
    
    print("Update with GDIGS sources")

    # Read the GDIGS data
    gdigsfile = 'data/gdigs_newregions7.tab'
    gdigs_new = np.genfromtxt(
        gdigsfile,
        delimiter="&",
        comments="%",
        dtype=None,
        autostrip=True,
        encoding="utf-8",
        names = ['gname', 'glong', 'glat', 'radius', 't_l', 'e_t_l', 'vlsr', 'e_vlsr', 'fwhm', 'e_fwhm', 'rms', 'snr']
        )

    # Fix names
    gdigs_new['gname'] = [s.replace('$' , '') for s in gdigs_new['gname']]
    gdigs_new['gname'] = [s.replace('*' , '') for s in gdigs_new['gname']]

    # get column names   
    query = "SELECT * FROM Detections WHERE name='G012.576+00.221'"
    c = conn.execute(query)
    columns = [column[0] for column in c.description]
    
    # add detections
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        for i in range(len(gdigs_new)):
            query = "UPDATE Detections SET" +\
                    " line_freq=5.76" +\
                    ", telescope='GBT'" +\
                    ", lines='H95-H117'" +\
                    ", author='Linville et al. (2023)'" +\
                    ", line=" + str(gdigs_new['t_l'][i]) +\
                    ", e_line=" + str(gdigs_new['e_t_l'][i]) +\
                    ", line_unit='mK'" +\
                    ", vlsr=" + str(gdigs_new['vlsr'][i]) +\
                    ", e_vlsr=" + str(gdigs_new['e_vlsr'][i]) +\
                    ", fwhm=" + str(gdigs_new['fwhm'][i]) +\
                    ", e_fwhm=" + str(gdigs_new['e_fwhm'][i]) +\
                    " WHERE name='" + gdigs_new['gname'][i] + "'"
            cur.execute(query)
        else:
            # copy previous row
            query = "INSERT INTO Detections (" + ", ".join(columns[1:]) + ") SELECT " +\
             ", ".join(columns[1:]) + " FROM Detections WHERE name='" +\
              gdigs_new['gname'][i][0:15] + "'"
            cur.execute(query)

            # insert new info
            query = "UPDATE Detections SET" +\
                    " line=" + str(gdigs_new['t_l'][i]) +\
                    ", e_line=" + str(gdigs_new['e_t_l'][i]) +\
                    ", line_unit='mK'" +\
                    ", vlsr=" + str(gdigs_new['vlsr'][i]) +\
                    ", e_vlsr=" + str(gdigs_new['e_vlsr'][i]) +\
                    ", fwhm=" + str(gdigs_new['fwhm'][i]) +\
                    ", e_fwhm=" + str(gdigs_new['e_fwhm'][i]) +\
                    " WHERE id=" + str(cur.lastrowid)
            cur.execute(query)
            
        print("Done!")