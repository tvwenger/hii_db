"""
wise.py

Utilities for adding WISE Catalog information to the database.

Copyright(C) 2020-2025 by
Trey V. Wenger; tvwenger@gmail.com
L. D. Anderson;
This code is licensed under MIT license (see LICENSE for details)
"""

import sqlite3
import numpy as np
import pandas as pd
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
    wise = pd.read_csv(wisefile, low_memory=False)
    print(len(wise))
    print(len(wise["GName"].unique()))

    # Populate Catalog
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        rows = [
            (
                w["GName"],
                w["Alias"],
                w["Name"],
                w["Catalog"],
                w["RA (J2000)"],
                w["Dec (J2000)"],
                w["GLong"],
                w["GLat"],
                w["Radius"],
                w["KDAR"],
                w["DMethod"],
                w["Author_KDAR"],
            )
            for i, w in wise.iterrows()
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
    wise = pd.read_csv(wisefile, low_memory=False, dtype={"Group": object})

    # Get group names/info
    all_groups = []
    group_info = []
    group_members = []
    for i, w in wise.iterrows():
        if not isinstance(w["Group"], str):
            continue
        group = w["Group"].strip()

        # check if group is already found
        if group in all_groups:
            # check that group data is the same
            idx = all_groups.index(group)
            if group_info[idx][0] == "":
                group_info[idx][0] = w["Group_VLSR"]
            elif isinstance(w["Group"], str) and group_info[idx][0] != w["Group_VLSR"]:
                print(
                    "Group {0} different VLSR".format(group)
                    + f" {w['Group_VLSR']} vs {group_info[idx][0]}"
                )
                continue

            if group_info[idx][1] == "":
                group_info[idx][1] = w["Group_e_VLSR"]
            elif (
                isinstance(w["Group_e_VLSR"], str)
                and group_info[idx][1] != w["Group_e_VLSR"]
            ):
                print("Group {0} different e_VLSR".format(group))
                continue

            if group_info[idx][2] == "":
                group_info[idx][2] = w["Group_KDAR"]
            elif (
                isinstance(w["Group_KDAR"], str)
                and group_info[idx][2] != w["Group_KDAR"]
            ):
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
                data[2],
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
    wise = pd.read_csv(wisefile, low_memory=False)

    # Add WISE Detections to table, split multi-component sources
    # into two rows in Detections table
    rows = []
    for wiseidx, w in wise.iterrows():
        # skip non-detections
        if np.isnan(w["VLSR"]):
            continue
        if w["VLSR"] == 0.0 and np.isnan(w["FWHM"]):
            continue

        # Get coordinate
        coord = SkyCoord(
            w["GLong_Observed"], w["GLat_Observed"], frame="galactic", unit="deg"
        )
        ra = coord.fk5.ra.deg
        dec = coord.fk5.dec.deg

        # Split multi-component sources, handle missing data
        vlsr_cols = ["VLSR"] + [f"VLSR{i}" for i in range(2, 5)]
        vlsrs = np.array([w[col] for col in vlsr_cols if not np.isnan(w[col])])
        e_vlsrs = np.array([w["e_" + col] for col in vlsr_cols if not np.isnan(w[col])])
        ncomp = len(vlsrs)

        fwhm_cols = ["FWHM"] + [f"FWHM{i}" for i in range(2, 5)]
        fwhms = np.array(
            [w[col] for col, vcol in zip(fwhm_cols, vlsr_cols) if not np.isnan(w[vcol])]
        )
        e_fwhms = np.array(
            [
                w["e_" + col]
                for col, vcol in zip(fwhm_cols, vlsr_cols)
                if not np.isnan(w[vcol])
            ]
        )

        tl_cols = ["TL"] + [f"TL{i}" for i in range(2, 5)]
        tls = np.array(
            [w[col] for col, vcol in zip(tl_cols, vlsr_cols) if not np.isnan(w[vcol])]
        )
        e_tls = np.array(
            [
                w["e_" + col]
                for col, vcol in zip(tl_cols, vlsr_cols)
                if not np.isnan(w[vcol])
            ]
        )

        # Fix Caswell & Haynes sources
        if "Caswell" in w["Author"]:
            tls = tls * 1000.0  # K to mK
            e_tls = np.array([10.0 for _ in e_tls])  # mK
            e_vlsrs = np.array([2.5 for _ in e_vlsrs])  # km/s
            e_fwhms = np.array([1.5 for _ in e_fwhms])  # km/s

        # Confirm that each parameter has the same number of components
        if np.any(
            np.array(
                [
                    np.isnan(vlsrs).sum(),
                    np.isnan(e_vlsrs).sum(),
                    np.isnan(fwhms).sum(),
                    np.isnan(e_fwhms).sum(),
                    np.isnan(tls).sum(),
                    np.isnan(e_tls).sum(),
                ]
            )
            != 0
        ):
            print("PROBLEM WITH {0} COMPONENTS".format(w["GName"]))
            print(vlsrs)
            print(e_vlsrs)
            print(fwhms)
            print(e_fwhms)
            print(tls)
            print(e_tls)
            print("================================")

        # Set other properties
        obs_type = None
        if w["Wavelength"] > 1.0e6:
            line_unit = "optical"
            cont_unit = "optical"
            line_freq = 456700000.0
        else:
            line_unit = "mK"
            cont_unit = "mK"
            obs_type = "peak"
            line_freq = 30000.0 / w["Wavelength"]  # cm -> MHz
        resolution = np.nan
        if not np.isnan(w["Resolution"]):
            resolution = float(w["Resolution"])  # arcmin
        beam_area = (
            np.pi * (resolution * 60.0) ** 2.0 / (4.0 * np.log(2.0))
        )  # sq. arcsec
        # area = np.pi * (w["HRDS_Size"]) ** 2.0 / (4.0 * np.log(2.0))  # sq. arcsec
        # area_unit = "arcsec2"
        area = 0.0
        if area == 0.0:
            area = None
            area_unit = None
        cont_freq = np.nan
        cont = np.nan
        e_cont = np.nan
        cont_unit = np.nan
        # if w["HRDS_Flux"] > 0.0:
        #     cont_freq = 9000.0
        #     cont = w["HRDS_Flux"]
        #     e_cont = w["HRDS_e_Flux"]
        #     cont_unit = "mJy"

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
                w["Te"],
                w["e_Te"],
                # w["Lines"],
                beam_area,
                obs_type,
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
        cont_freq, cont, e_cont, cont_unit, area, area_unit, te, e_te,
        beam_area, type, telescope, author, source) VALUES
        (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
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
