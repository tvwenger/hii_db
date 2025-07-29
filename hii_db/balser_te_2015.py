"""
balser_te_2015.py

Utilities for adding Balser+2015 data to the database.

Copyright(C) 2020-2025 by
Trey V. Wenger; tvwenger@gmail.com
L. D. Anderson;
This code is licensed under MIT license (see LICENSE for details)
"""

import os
import numpy as np
import sqlite3
from astropy.coordinates import SkyCoord


def add_detections(db, datatype, data_dir="data"):
    """
    Add Balser+2015 GBT and 140 Foot data to Detections table.
    Also populate Catalog->Detections.

    Inputs:
        db :: string
            database filename
        datatype :: string
            Data to add. Either "GBT" or "140 Foot"
        data_dir :: string
            Path to data directory

    Returns: Nothing
    """
    print("Adding Balser+2015 {0} data to Detections...".format(datatype))
    if datatype == "GBT":
        fname = os.path.join(data_dir, "te", "balser_2015", "gbtRRL_balser2015.csv")
    elif datatype == "140 Foot":
        fname = os.path.join(data_dir, "te", "balser_2015", "140foot_balser2015.csv")
    else:
        raise ValueError("Invalid data type: {0}".format(datatype))
    data = np.genfromtxt(
        fname, dtype=None, names=True, delimiter=",", autostrip=True, encoding="UTF-8"
    )
    coords = SkyCoord(data["glong"], data["glat"], frame="galactic", unit="deg")

    # Loop over rows in data
    rows = []
    for dat, coord in zip(data, coords):
        ra = coord.fk5.ra.deg
        dec = coord.fk5.dec.deg
        glong = coord.galactic.l.deg
        glat = coord.galactic.b.deg

        # Convert peak fluxes to mJy assuming 2 K/Jy for GBT and
        # 0.27 K/Jy for 140 Foot
        if datatype == "GBT":
            line_peak = dat["tl"] / 2.0
            e_line_peak = dat["e_tl"] / 2.0
            cont_peak = dat["tc"] / 2.0
            e_cont_peak = dat["e_tc"] / 2.0
        elif datatype == "140 Foot":
            line_peak = dat["tl"] / 0.27
            e_line_peak = dat["e_tl"] / 0.27
            cont_peak = dat["tc"] / 0.27
            e_cont_peak = dat["e_tc"] / 0.27
        if cont_peak > 0.0:
            linetocont_peak = line_peak / cont_peak
            e_linetocont_peak = np.sqrt(
                linetocont_peak**2.0
                * (
                    e_line_peak**2.0 / line_peak**2.0
                    + e_cont_peak**2.0 / cont_peak**2.0
                )
            )
        else:
            linetocont_peak = None
            e_linetocont_peak = None

        # Convert total fluxes using Balser+1995 and assuming
        # eta_b = 0.92 for X-band GBT
        # eta_b = 0.70 for X-band 140 Foot
        # frequency = 8556 MHz
        # wavelength = 3.5 cm
        if dat["cont_width"] > 0.0:
            if datatype == "GBT":
                eta_b = 0.92
            elif datatype == "140 Foot":
                eta_b = 0.70
            line_total = (
                2.647 * dat["tl"] / eta_b * (dat["cont_width"] / 60.0 / 3.5) ** 2.0
            )
            e_line_total = np.sqrt(
                line_total**2.0
                * (
                    (dat["e_tl"] / dat["tl"]) ** 2.0
                    + 4.0 * (dat["e_cont_width"] / dat["cont_width"]) ** 2.0
                )
            )
            cont_total = (
                2.647 * dat["tc"] / eta_b * (dat["cont_width"] / 60.0 / 3.5) ** 2.0
            )
            e_cont_total = np.sqrt(
                cont_total**2.0
                * (
                    (dat["e_tc"] / dat["tc"]) ** 2.0
                    + 4.0 * (dat["e_cont_width"] / dat["cont_width"]) ** 2.0
                )
            )
            linetocont_total = line_total / cont_total
            e_linetocont_total = np.sqrt(
                linetocont_total**2.0
                * (
                    e_line_total**2.0 / line_total**2.0
                    + e_cont_total**2.0 / cont_total**2.0
                )
            )

            # Compute area
            area = np.pi * dat["cont_width"] ** 2.0 / (4 * np.log(2.0))
            e_area = 2.0 * area * (dat["e_cont_width"] / dat["cont_width"])
        else:
            cont_peak = None
            e_cont_peak = None
            cont_total = None
            e_cont_total = None
            line_total = None
            e_line_total = None
            linetocont_total = None
            e_linetocont_total = None
            area = None
            e_area = None

        # Average line frequency and other telescope-dependent info
        cont_freq = 8665.0
        if datatype == "GBT":
            line_freq = 8902.0
            line_qf = 1
            beam_area = 8180.0
            lines = "H87-H93"
            telescope = "GBT"
            author = "Balser et al. (2011, 2015)"
        elif datatype == "140 Foot":
            line_freq = 8447.0
            line_qf = dat["lqf"]
            beam_area = 41770.0
            lines = "H91;H92"
            telescope = "NRAO 140 Foot"
            author = "Quireza et al. (2006a), Balser et al. (2015)"

        # Create peak row
        row = (
            dat["gname"],
            ra,
            dec,
            glong,
            glat,
            line_freq,
            line_peak,
            e_line_peak,
            "mJy/beam",
            dat["vlsr"],
            dat["e_vlsr"],
            dat["fwhm"],
            dat["e_fwhm"],
            int(line_qf),
            cont_freq,
            cont_peak,
            e_cont_peak,
            "mJy/beam",
            int(dat["cqf"]),
            linetocont_peak,
            e_linetocont_peak,
            dat["te"],
            dat["e_te"],
            None,
            None,
            None,
            beam_area,
            lines,
            telescope,
            author,
            "Balser+2015",
            "peak",
        )
        rows.append(row)

        # Create total row
        row = (
            dat["gname"],
            ra,
            dec,
            glong,
            glat,
            line_freq,
            line_total,
            e_line_total,
            "mJy",
            dat["vlsr"],
            dat["e_vlsr"],
            dat["fwhm"],
            dat["e_fwhm"],
            int(line_qf),
            cont_freq,
            cont_total,
            e_cont_total,
            "mJy",
            int(dat["cqf"]),
            linetocont_total,
            e_linetocont_total,
            dat["te"],
            dat["e_te"],
            area,
            e_area,
            "arcsec2",
            beam_area,
            lines,
            telescope,
            author,
            "Balser+2015",
            "peak",
        )
        rows.append(row)

    # Populate Detections table
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        cur.executemany(
            """
        INSERT INTO Detections
        (name, ra, dec, glong, glat, line_freq,
        line, e_line, line_unit, vlsr, e_vlsr, fwhm, e_fwhm,
        line_qf, cont_freq, cont, e_cont, cont_unit,
        cont_qf, linetocont, e_linetocont, te, e_te,
        area, e_area, area_unit, beam_area, lines, telescope,
        author, source, type) VALUES
        (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
            rows,
        )
    print("Done!")
    print()

    # Match Balser+2015 detections to Catalog
    print("Matching {0} Detections to WISE Catalog...".format(datatype))
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")

        # Get WISE catalog gnames and positions
        cur.execute("SELECT id, gname, alias, ra, dec, catalog, radius FROM Catalog")
        wise = np.array(
            cur.fetchall(),
            dtype=[
                ("id", "i"),
                ("gname", "U100"),
                ("alias", "U100"),
                ("ra", "f8"),
                ("dec", "f8"),
                ("catalog", "U1"),
                ("radius", "f8"),
            ],
        )
        wise_coords = SkyCoord(wise["ra"], wise["dec"], frame="fk5", unit="deg")

        # Get the Balser+2015 detection gnames
        if datatype == "GBT":
            telescope = "GBT"
        elif datatype == "140 Foot":
            telescope = "NRAO 140 Foot"
        cur.execute(
            """
        SELECT id, name, ra, dec FROM Detections
        WHERE source="Balser+2015" AND telescope=?
        """,
            [telescope],
        )
        dets = np.array(
            cur.fetchall(),
            dtype=[("id", "i"), ("name", "U100"), ("ra", "f8"), ("dec", "f8")],
        )
        det_coords = SkyCoord(dets["ra"], dets["dec"], frame="fk5", unit="deg")

    # Loop over detections and find matches
    rows = []
    for det, coord in zip(dets, det_coords):

        # Match based on GName
        catalog_id = None
        match = np.where(wise["gname"] == det["name"])[0]
        if len(match) == 1:
            match = match[0]
            catalog_id = wise["id"][match]

        # If that doesn't work, matched based on alias
        catalog_id = None
        match = np.where(wise["alias"] == det["name"])[0]
        if len(match) == 1:
            match = match[0]
            catalog_id = wise["id"][match]

        # If that doesn't work, match based on position
        # Within smallest WISE region
        if catalog_id is None:
            seps = coord.separation(wise_coords).arcsec
            # within WISE region
            inside = np.where(seps - wise["radius"] < 0.0)[0]
            if len(inside) > 0:
                # pick smallest region containing the source
                minsizeind = np.argmin(wise["radius"][inside])
                match = inside[minsizeind]
                catalog_id = wise["id"][match]

        # If that doesn't work, match based on position
        # Closest within 1.5x WISE region radius
        if catalog_id is None:
            seps = coord.separation(wise_coords).arcsec
            # within WISE region
            inside = np.where(seps - 1.5 * wise["radius"] < 0.0)[0]
            if len(inside) > 0:
                # pick smallest region containing the source
                minsizeind = np.argmin(wise["radius"][inside])
                match = inside[minsizeind]
                catalog_id = wise["id"][match]

        # No match :(
        if catalog_id is None:
            print("Found no Catalog match for {0} {1}".format(datatype, det["name"]))
            print(
                "Min sep: {0} {1:.2f} arcsec".format(
                    wise["gname"][np.argmin(seps)], np.min(seps)
                )
            )
            continue

        # Compute separation to Catalog match
        sep = coord.separation(wise_coords[match]).arcsec
        rows.append((int(catalog_id), int(det["id"]), sep))

    # Populate Catalog-Detections
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        cur.executemany(
            """
        INSERT INTO CatalogDetections
        (catalog_id, detection_id, separation) VALUES (?, ?, ?)
        """,
            rows,
        )
    print("Done!")
    print()
