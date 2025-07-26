"""
wenger_shrds_2021.py

Utilities for adding Wenger+2021 SHRDS Full Catalog data to the
database.

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

import os
import numpy as np
import glob
import sqlite3
from astropy.coordinates import SkyCoord
from astropy.io import fits
from .utils import parse_region_coord

# These are the SHRDS calibrators we're ignoring
calibrators = [
    "0823-500",
    "1253-055",
    "1421-490",
    "1934-638",
    "0537-441",
    "0723-008",
    "0727-115",
    "1129-58",
    "1148-671",
    "0906-47",
    "1036-52",
    "1613-586",
    "1714-336",
    "1714-397",
    "j1322-6532",
    "orion-kl",
    "1829-207",
    "shrds073",
]


def add_detections(db):
    """
    Read SHRDS detections and populate detections table.
    Also populate Catalog->Detections, Fields, and Fields->Detections.

    Inputs:
        db :: string
            Database filename

    Returns: Nothing
    """
    print("Adding Wenger+2021 SHRDS Fields...")

    # Read quality factor data
    qf_data = np.genfromtxt(
        "data/wenger_shrds_2021/shrds_qfs.txt", dtype=None, names=True, encoding="utf-8"
    )

    # Loop over Fields
    data = []
    fields = glob.glob("data/wenger_shrds_2021/*")
    fields.sort()
    for field in fields:
        if not os.path.isdir(field):
            continue
        myfield = field.split("/")[-1]
        if myfield in calibrators:
            # skip calibrators
            continue

        # Get field center position from continuum image
        img = glob.glob(os.path.join(field, "*.fits"))
        if len(img) != 1:
            raise ValueError("Issue with FITS images: {0}".format(field))
        hdu = fits.open(img[0])[0]
        field_coord = SkyCoord(
            hdu.header["CRVAL1"], hdu.header["CRVAL2"], unit="deg", frame="fk5"
        )
        # HPBW of ATCA at 7 GHz is ~500 arcsec
        hpbw = 500.0
        # mosaics are undefined
        if "mos" in myfield:
            hpbw = np.nan
        row = (
            myfield,
            field_coord.fk5.ra.to("deg").value,
            field_coord.fk5.dec.to("deg").value,
            field_coord.galactic.l.to("deg").value,
            field_coord.galactic.b.to("deg").value,
            hpbw,
        )

        # Add row and get insert ID
        lastid = -1
        with sqlite3.connect(db) as conn:
            cur = conn.cursor()
            cur.execute("PRAGMA foreign_keys = ON")
            cur.execute(
                """
            INSERT INTO Fields
            (name, ra, dec, glong, glat, hpbw) VALUES (?,?,?,?,?,?)
            """,
                row,
            )
            lastid = cur.lastrowid

        # Find region files
        regs = glob.glob(os.path.join(field, "*.rgn"))
        regs.sort()

        # Loop over regions, read data, and save field/region
        for reg in regs:
            myreg = reg.split("/")[-1]
            data.append([myfield, lastid, myreg])
    print("Done!")
    print()

    # Loop over regions, populate rows
    line_rows = []
    cont_rows = []
    print("Adding Wenger+2021 SHRDS data to Detections and FieldsDetections...")
    for field, field_id, reg in data:
        gname = reg
        gname = gname.replace(".notaper.rgn", "")
        gname = gname.replace(".notaper.imsmooth.rgn", "")

        # Get position from region file
        fname = os.path.join("data", "wenger_shrds_2021", field, reg)
        coord = parse_region_coord(fname)
        ra = coord.fk5.ra.deg
        dec = coord.fk5.dec.deg
        glong = coord.galactic.l.deg
        glat = coord.galactic.b.deg

        # Loop over peak/total
        for datatype in ["peak", "total"]:
            # Read continuum info
            fname = os.path.join(
                "data",
                "wenger_shrds_2021",
                field,
                "{0}.clean.{1}.continfo.txt".format(reg, datatype),
            )
            cont_data = np.genfromtxt(fname, dtype=None, names=True, encoding="UTF-8")

            # Read spectrum info
            fname = os.path.join(
                "data",
                "wenger_shrds_2021",
                field,
                "{0}.clean.{1}.wt.specinfo.txt".format(reg, datatype),
            )
            spec_data = np.genfromtxt(fname, dtype=None, names=True, encoding="UTF-8")

            # Get QFs
            ind = np.where((qf_data["field"] == field) * (qf_data["source"] == gname))[
                0
            ]
            if len(ind) != 1:
                print(ind)
                raise ValueError(
                    "No match for SHRDS Full Catalog {0} {1} in QF file".format(
                        field, gname
                    )
                )
            ind = ind[0]
            cont_qf = qf_data[ind]["imqf"]
            if "imsmooth" in reg:
                cont_qf = qf_data[ind]["smoimqf"]
            line_qf = "A"
            if qf_data[ind]["blend"] in ["Y", "X"]:
                line_qf = "F"

            # Change QFs to numbers
            qf_letter = ["A", "B", "C", "D", "F"]
            qf_number = [1, 2, 3, 4, 5]
            cont_qf = qf_number[qf_letter.index(cont_qf)]
            line_qf = qf_number[qf_letter.index(line_qf)]

            # Add line rows
            taper = "notaper"
            if "notaper.imsmooth" in reg:
                taper = "notapsm"
            unit = "mJy/beam"
            if datatype == "total":
                unit = "mJy"
            for spec in spec_data:
                # Get lineid and component
                lineid = spec["lineid"]
                comp = None
                if "(" in spec["lineid"]:
                    comp = spec["lineid"][-2]
                    lineid = spec["lineid"][:-3]
                if lineid == "all":
                    lineid = "H88-H112"
                lineid = lineid.replace("a", "")

                # match closest frequency in cont_data to estimate beam size and region size
                match = np.argmin(np.abs(spec["frequency"] - cont_data["frequency"]))
                area = np.nan
                if datatype == "total":
                    area = cont_data["area_arcsec"][match]
                beam_area = cont_data["beam_arcsec"][match]

                # Fix small errors
                e_velo = np.max([spec["e_velo"], 0.1])
                e_line = np.max([spec["e_line"], 0.01])
                e_fwhm = np.max([spec["e_fwhm"], 0.01])
                rms = np.max([spec["rms"], 0.01])
                e_line2cont = np.max([spec["e_line2cont"], 0.0001])
                e_elec_temp = np.max([spec["e_elec_temp"], 0.1])

                # Catch bad data
                my_line_qf = line_qf
                if np.isnan(spec["line"]):
                    my_line_qf = 5
                line = spec["line"]
                velo = spec["velo"]
                fwhm = spec["fwhm"]
                elec_temp = spec["elec_temp"]
                line2cont = spec["line2cont"]
                linesnr = spec["linesnr"]
                if (
                    e_line > line
                    or e_fwhm > fwhm
                    or e_elec_temp > elec_temp
                    or e_line2cont > line2cont
                    or elec_temp < 1000.0
                    or line2cont > 0.5
                ):
                    line = None
                    e_line = None
                    velo = None
                    e_velo = None
                    fwhm = None
                    e_fwhm = None
                    elec_temp = None
                    e_elec_temp = None
                    line2cont = None
                    e_line2cont = None
                    linesnr = None
                    line_qf = 5
                row = (
                    gname,
                    ra,
                    dec,
                    glong,
                    glat,
                    spec["frequency"],
                    comp,
                    line,
                    e_line,
                    unit,
                    velo,
                    e_velo,
                    fwhm,
                    e_fwhm,
                    rms,
                    linesnr,
                    my_line_qf,
                    spec["frequency"],
                    spec["cont"],
                    spec["rms"],
                    unit,
                    cont_qf,
                    area,
                    "arcsec2",
                    beam_area,
                    line2cont,
                    e_line2cont,
                    elec_temp,
                    e_elec_temp,
                    lineid,
                    "ATCA",
                    "Wenger et al. (2021)",
                    "SHRDS Full Catalog",
                    datatype,
                    taper,
                    int(field_id),
                )
                line_rows.append(row)

            # Add cont rows
            taper = "notaper"
            if "notaper.imsmooth" in reg:
                taper = "notapsm"
            unit = "mJy/beam"
            if datatype == "total":
                unit = "mJy"
            for cont in cont_data:
                # Get lineid
                lineid = "cont"
                if cont["spw"] != "cont":
                    lineid = "spw{0}".format(cont["spw"])

                # Area is undefinied for peak
                area = np.nan
                if datatype == "total":
                    area = cont["area_arcsec"]

                # Save data
                row = (
                    gname,
                    ra,
                    dec,
                    glong,
                    glat,
                    cont["frequency"],
                    cont["cont"],
                    cont["e_cont_A"],
                    unit,
                    cont_qf,
                    area,
                    "arcsec2",
                    cont["beam_arcsec"],
                    lineid,
                    "ATCA",
                    "Wenger et al. (2021)",
                    "SHRDS Full Catalog",
                    datatype,
                    taper,
                    int(field_id),
                )
                cont_rows.append(row)

    # Populate Detections table
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        cur.executemany(
            """
        INSERT INTO Detections
        (name, ra, dec, glong, glat, line_freq, component,
        line, e_line, line_unit, vlsr, e_vlsr, fwhm, e_fwhm, spec_rms,
        line_snr, line_qf, cont_freq, cont, e_cont, cont_unit,
        cont_qf, area, area_unit, beam_area, linetocont, e_linetocont, te, e_te,
        lines, telescope, author, source, type, taper, field_id) VALUES
        (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
            line_rows,
        )

        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        cur.executemany(
            """
        INSERT INTO Detections
        (name, ra, dec, glong, glat, cont_freq, cont, e_cont, cont_unit,
        cont_qf, area, area_unit, beam_area,
        lines, telescope, author, source, type, taper, field_id) VALUES
        (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
            cont_rows,
        )
    print("Done!")
    print()

    # Match SHRDS detections to Catalog
    print("Matching Wenger+2021 SHRDS Detections to WISE Catalog...")
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

        # Get the SHRDS detection gnames
        cur.execute(
            """
        SELECT id, name, ra, dec FROM Detections
        WHERE source="SHRDS Full Catalog"
        """
        )
        dets = np.array(
            cur.fetchall(),
            dtype=[("id", "i"), ("name", "U100"), ("ra", "f8"), ("dec", "f8")],
        )
        det_coords = SkyCoord(dets["ra"], dets["dec"], frame="fk5", unit="deg")

        # Match and calculate separations
        nomatch = []
        rows = []
        for det, coord in zip(dets, det_coords):
            # Fix a few with bad region gnames
            check_gname = det["name"]
            match = np.where(wisecat["gname"] == check_gname)[0]
            if len(match) == 0:
                nomatch.append(check_gname)
                continue
            if len(match) > 1:
                print("Multiple matches for SHRDS Full Catalog {0}".format(det["name"]))
                continue
            match = match[0]
            sep = coord.separation(cat_coords[match]).arcsec
            rows.append((int(wisecat["id"][match]), int(det["id"]), sep))
        if len(nomatch) > 0:
            print("Could not find SHRDS Full Catalog matches for:")
            print(np.unique(nomatch))

        # Populate Catalog-Detections
        cur.executemany(
            """
        INSERT INTO CatalogDetections
        (catalog_id, detection_id, separation) VALUES (?, ?, ?)
        """,
            rows,
        )
    print("Done!")
    print()
