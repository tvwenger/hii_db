"""
wenger_shrds_2019.py

Utilities for adding Wenger+2019 SHRDS Bright Catalog data to the
database.

Copyright(C) 2020-2025 by
Trey V. Wenger; tvwenger@gmail.com
L. D. Anderson;
This code is licensed under MIT license (see LICENSE for details)
"""

import os
import numpy as np
import glob
import sqlite3
from astropy.coordinates import SkyCoord
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
]


def add_detections(db, data_dir="data"):
    """
    Read SHRDS detections and populate detections table.
    Also populate Catalog->Detections.

    Inputs:
        db :: string
            Database filename.
        data_dir :: string
            Path to data directory

    Returns: Nothing
    """
    print("Adding Wenger+2019 SHRDS data to Detections...")

    # Read quality factor data
    qf_data_notap = np.genfromtxt(
        os.path.join(
            data_dir, "rrl_surveys", "wenger_shrds_2019", "shrds_qfs_notap.txt"
        ),
        dtype=None,
        names=True,
        encoding="utf-8",
    )
    qf_data_uvtap = np.genfromtxt(
        os.path.join(
            data_dir, "rrl_surveys", "wenger_shrds_2019", "shrds_qfs_uvtap.txt"
        ),
        dtype=None,
        names=True,
        encoding="utf-8",
    )

    # Loop over epochs
    data = []
    fields = glob.glob(os.path.join(data_dir, "rrl_surveys", "wenger_shrds_2019", "*"))
    fields.sort()
    for field in fields:
        if not os.path.isdir(field):
            continue
        myfield = field.split("/")[-1]
        if myfield in calibrators:
            # skip calibrators
            continue

        # Find region files
        regs = glob.glob(os.path.join(field, "*.rgn"))
        regs.sort()

        # Loop over regions, read data, and save sb/field/region
        for reg in regs:
            myreg = reg.split("/")[-1]
            data.append([myfield, myreg])

    # Loop over regions, populate rows
    rows = []
    for field, reg in data:
        gname = reg
        gname = gname.replace(".notaper.rgn", "")
        gname = gname.replace(".uvtaper.rgn", "")

        # Get position from region file
        fname = os.path.join(data_dir, "rrl_surveys", "wenger_shrds_2019", field, reg)
        coord = parse_region_coord(fname)
        ra = coord.fk5.ra.deg
        dec = coord.fk5.dec.deg
        glong = coord.galactic.l.deg
        glat = coord.galactic.b.deg

        # Read spectrum info
        fname = os.path.join(
            data_dir,
            "rrl_surveys",
            "wenger_shrds_2019",
            field,
            "{0}.clean.wt.specinfo.txt".format(reg),
        )
        qf_data = qf_data_notap
        if "uvtaper" in reg:
            fname = os.path.join(
                data_dir,
                "rrl_surveys",
                "wenger_shrds_2019",
                field,
                "{0}.clean.uvtaper.wt.specinfo.txt".format(reg),
            )
            qf_data = qf_data_uvtap
        spec_data = np.genfromtxt(fname, dtype=None, names=True, encoding="UTF-8")

        # Add rows
        taper = "notaper"
        if "uvtaper" in reg:
            taper = "uvtaper"
        unit = "mJy/beam"
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

            # Get QFs
            check_comp = "a" if comp is None else comp
            ind = np.where(
                (qf_data["field"] == field)
                * (qf_data["gname"] == gname)
                * (qf_data["comp"] == check_comp)
            )[0]
            if len(ind) != 1:
                print(ind)
                raise ValueError(
                    "No match for Bright Catalog {0} {1} {2} in QF file".format(
                        field, reg, comp
                    )
                )
            ind = ind[0]
            cont_qf = qf_data[ind]["cont_qf"]
            line_qf = qf_data[ind]["rrl_qf"]

            # Change QFs to numbers
            qf_letter = ["A", "B", "C", "D", "F"]
            qf_number = [1, 2, 3, 4, 5]
            cont_qf = qf_number[qf_letter.index(cont_qf)]
            line_qf = qf_number[qf_letter.index(line_qf)]
            if np.isnan(spec["line"]):
                line_qf = 5

            # Fix small errors
            e_velo = np.max([spec["e_velo"], 0.1])
            e_line = np.max([spec["e_line"], 0.1])
            e_fwhm = np.max([spec["e_fwhm"], 0.1])
            rms = np.max([spec["rms"], 0.1])
            e_line2cont = np.max([spec["e_line2cont"], 0.001])
            e_elec_temp = np.max([spec["e_elec_temp"], 0.1])

            # Catch bad data
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
                "mJy/beam",
                velo,
                e_velo,
                fwhm,
                e_fwhm,
                rms,
                linesnr,
                line_qf,
                spec["frequency"],
                spec["cont"],
                spec["rms"],
                unit,
                cont_qf,
                line2cont,
                e_line2cont,
                elec_temp,
                e_elec_temp,
                lineid,
                "ATCA",
                "Wenger et al. (2019)",
                "SHRDS Bright Catalog",
                taper,
            )
            rows.append(row)

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
        cont_qf, linetocont, e_linetocont, te, e_te,
        lines, telescope, author, source, taper) VALUES
        (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
            rows,
        )
    print("Done!")
    print()

    # Match SHRDS detections to Catalog
    print("Matching Wenger+2019 SHRDS Detections to WISE Catalog...")
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
        WHERE source="SHRDS Bright Catalog"
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
            if det["name"] == "G352.313-00.440":
                check_gname = "G352.313-00.442"
            if det["name"] == "G359.467-00.172":
                check_gname = "G359.467-00.173"
            if det["name"] == "G348.891-00.179":
                check_gname = "G348.892-00.179"
            if det["name"] == "G312.037+00.084":
                check_gname = "G312.091+00.069"
            if det["name"] == "G337.617-00.064":
                check_gname = "G337.617-00.065"
            if det["name"] == "G288.233-01.117":
                check_gname = "G288.233-01.118"
            match = np.where(wisecat["gname"] == check_gname)[0]
            if len(match) == 0:
                nomatch.append(check_gname)
                continue
            if len(match) > 1:
                print(
                    "Multiple matches for SHRDS Bright Catalog {0}".format(det["name"])
                )
                continue
            match = match[0]
            sep = coord.separation(cat_coords[match]).arcsec
            rows.append((int(wisecat["id"][match]), int(det["id"]), sep))

        if len(nomatch) > 0:
            print("Could not find SHRDS Bright Catalog matches for:")
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
