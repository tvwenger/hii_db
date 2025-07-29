"""
wenger_te_2019.py

Utilities for adding Wenger+2019 VLA data to the database.

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


def add_detections(db, data_dir="data"):
    """
    Read VLA detections and populate detections table.
    Also populate Catalog->Detections.

    Inputs:
        db :: string
            Database filename
        data_dir :: string
            Path to data directory

    Returns: Nothing
    """
    print("Adding Wenger+2019 VLA data to Detections...")

    # Read quality factor data
    qf_data = np.genfromtxt(
        os.path.join(data_dir, "te", "wenger_te_2019", "te_qfs.txt"),
        dtype=None,
        names=True,
        usecols=(0, 1, 2, 3, 4, 5, 6),
        encoding="utf-8",
    )

    # Loop over SBs
    data = []
    sbs = glob.glob(os.path.join(data_dir, "te", "wenger_te_2019", "*"))
    sbs.sort()
    for sb in sbs:
        if not os.path.isdir(sb):
            continue
        mysb = sb.split("/")[-1]

        # Loop over fields
        fields = glob.glob(os.path.join(sb, "*"))
        fields.sort()
        for field in fields:
            if not os.path.isdir(field):
                continue
            myfield = field.split("/")[-1]
            if myfield[0] != "G":
                # skip calibrators
                continue

            # Find region files
            regs = glob.glob(os.path.join(field, "*.rgn"))
            regs.sort()

            # Loop over regions, read data, and save sb/field/region
            for reg in regs:
                myreg = reg.split("/")[-1]
                data.append([mysb, myfield, myreg])

    # Loop over regions, populate rows
    rows = []
    for sb, field, reg in data:
        gname = reg
        gname = gname.replace(".notaper.rgn", "")
        gname = gname.replace(".uvtaper.rgn", "")
        gname = gname.replace(".uvtaper.imsmooth.rgn", "")

        # Get position from region file
        fname = os.path.join(data_dir, "te", "wenger_te_2019", sb, field, reg)
        coord = parse_region_coord(fname)
        ra = coord.fk5.ra.deg
        dec = coord.fk5.dec.deg
        glong = coord.galactic.l.deg
        glat = coord.galactic.b.deg

        # Loop over peak/total
        for datatype in ["peak", "total"]:
            # Read continuum info
            fname = os.path.join(
                data_dir,
                "te",
                "wenger_te_2019",
                sb,
                field,
                "{0}.clean.{1}.continfo.txt".format(reg, datatype),
            )
            cont_data = np.genfromtxt(fname, dtype=None, names=True, encoding="UTF-8")

            # Read spectrum info
            fname = os.path.join(
                data_dir,
                "te",
                "wenger_te_2019",
                sb,
                field,
                "{0}.clean.{1}.wt.specinfo.txt".format(reg, datatype),
            )
            spec_data = np.genfromtxt(fname, dtype=None, names=True, encoding="UTF-8")

            # Get QFs
            ind = np.where(
                (qf_data["SB"] == sb)
                * (qf_data["field"] == field)
                * (qf_data["source"] == gname)
            )[0]
            if len(ind) != 1:
                print(ind)
                raise ValueError(
                    "No match for VLA {0} {1} {2} in QF file".format(sb, field, gname)
                )
            ind = ind[0]
            cont_qf = qf_data[ind]["cqf_notap"]
            line_qf = qf_data[ind]["lqf_notap"]
            if "uvtaper" in reg:
                cont_qf = qf_data[ind]["cqf_uvtap"]
                line_qf = qf_data[ind]["lqf_uvtap"]

            # Change QFs to numbers
            qf_letter = ["A", "B", "C", "D", "F"]
            qf_number = [1, 2, 3, 4, 5]
            cont_qf = qf_number[qf_letter.index(cont_qf)]
            line_qf = qf_number[qf_letter.index(line_qf)]

            # Add rows
            taper = "notaper"
            if "uvtaper" in reg:
                taper = "uvtaper"
            if "uvtaper.imsmooth" in reg:
                taper = "uvtapsm"
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
                    lineid = "H87-H93"
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
                row = (
                    gname,
                    ra,
                    dec,
                    glong,
                    glat,
                    spec["frequency"],
                    comp,
                    spec["line"],
                    e_line,
                    unit,
                    spec["velo"],
                    e_velo,
                    spec["fwhm"],
                    e_fwhm,
                    rms,
                    spec["linesnr"],
                    line_qf,
                    spec["frequency"],
                    spec["cont"],
                    spec["rms"],
                    unit,
                    cont_qf,
                    area,
                    "arcsec2",
                    beam_area,
                    spec["line2cont"],
                    e_line2cont,
                    spec["elec_temp"],
                    e_elec_temp,
                    lineid,
                    "JVLA",
                    "Wenger et al. (2019)",
                    "VLA Te",
                    datatype,
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
        cont_qf, area, area_unit, beam_area, linetocont, e_linetocont, te, e_te,
        lines, telescope, author, source, type, taper) VALUES
        (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
            rows,
        )
    print("Done!")
    print()

    # Match VLA detections to Catalog
    print("Matching VLA Detections to WISE Catalog...")
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

        # Get the VLA detection gnames
        cur.execute(
            """
        SELECT id, name, ra, dec FROM Detections
        WHERE source="VLA Te"
        """
        )
        det = np.array(
            cur.fetchall(),
            dtype=[("id", "i"), ("name", "U100"), ("ra", "f8"), ("dec", "f8")],
        )
        det_coords = SkyCoord(det["ra"], det["dec"], frame="fk5", unit="deg")

        # Match and calculate separations
        bad = np.zeros(len(det), dtype=bool)
        for i, gname in enumerate(det["name"]):
            if gname not in wisecat["gname"]:
                bad[i] = True
                print(f"No match for {gname}")
        det = det[~bad]
        det_coords = det_coords[~bad]

        matches = np.array(
            [np.where(wisecat["gname"] == gname)[0][0] for gname in det["name"]]
        )
        seps = [
            coord.separation(cat_coords[match]).arcsec
            for coord, match in zip(det_coords, matches)
        ]

        # Populate Catalog-Detections
        rows = [
            (int(wisecat["id"][match]), int(det["id"]), float(sep))
            for det, match, sep in zip(det, matches, seps)
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
