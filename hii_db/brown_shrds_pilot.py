"""
wenger_brown_pilot.py

Utilities for adding Brown+2017 SHRDS pilot survey data to the
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
import sqlite3
from astropy.coordinates import SkyCoord


def add_detections(db):
    """
    Read SHRDS Pilot data and populate detections table.
    Also populate Catalog->Detections.

    Inputs:
        db :: string
            Database filename
    """
    print("Adding Brown+2017 SHRDS Pilot data to Detections...")

    # Read data
    data = np.genfromtxt(
        os.path.join("data", "brown_shrds_pilot", "shrds_pilot.txt"),
        dtype=[
            ("gname", "U15"),
            ("lineid", "U3"),
            ("vlsr", "f8"),
            ("e_vlsr", "f8"),
            ("fwhm", "f8"),
            ("e_fwhm", "f8"),
            ("cont_flux", "f8"),
            ("rms", "f8"),
            ("rrl_flux", "f8"),
            ("e_rrl_flux", "f8"),
            ("te", "f8"),
            ("e_te", "f8"),
            ("rrl_snr", "f8"),
        ],
    )
    glong = np.array([float(gname[1:8]) for gname in data["gname"]])
    glat = np.array([float(gname[8:]) for gname in data["gname"]])
    coords = SkyCoord(glong, glat, frame="galactic", unit="deg")

    # Populate rows
    rows = []
    for dat, coord in zip(data, coords):
        ra = coord.fk5.ra.to("deg").value
        dec = coord.fk5.dec.to("deg").value
        glong = coord.galactic.l.to("deg").value
        glat = coord.galactic.b.to("deg").value
        row = (
            dat["gname"],
            ra,
            dec,
            glong,
            glat,
            dat["rrl_flux"],
            dat["e_rrl_flux"],
            "mJy/beam",
            dat["vlsr"],
            dat["e_vlsr"],
            dat["fwhm"],
            dat["e_fwhm"],
            dat["rms"],
            dat["rrl_snr"],
            int(1),
            int(5),
            dat["te"],
            dat["e_te"],
            "H" + dat["lineid"],
            "ATCA",
            "Brown et al. (2017)",
            "SHRDS Pilot",
        )
        rows.append(row)

    # Populate Detections table
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        cur.executemany(
            """
        INSERT INTO Detections
        (name, ra, dec, glong, glat, line, e_line, line_unit,
        vlsr, e_vlsr, fwhm, e_fwhm, spec_rms,
        line_snr, line_qf, cont_qf, te, e_te,
        lines, telescope, author, source) VALUES
        (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
            rows,
        )
    print("Done!")
    print()

    # Match SHRDS Pilot detections to Catalog
    print("Matching Brown+2017 SHRDS Pilot Detections to WISE Catalog...")
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
        WHERE source="SHRDS Pilot"
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
            # Fix bad gnames
            check_gname = det["name"]
            if det["name"] == "G324.662-00.331":
                check_gname = "G324.642-00.321"
            if det["name"] == "G327.313-00.536":
                check_gname = "G327.300-00.548"

            # Match based on gname
            catalog_id = None
            match = np.where(wisecat["gname"] == check_gname)[0]
            if len(match) == 1:
                match = match[0]
                catalog_id = wisecat["id"][match]
            else:
                # Match based on separation
                seps = coord.separation(cat_coords).to("arcsec").value
                if np.min(seps) < 60.0:
                    match = np.argmin(seps)
                    catalog_id = wisecat["id"][match]

            if catalog_id is None:
                nomatch.append(det["name"])
                continue
            sep = coord.separation(cat_coords[match]).arcsec
            rows.append((int(catalog_id), int(det["id"]), sep))

        if len(nomatch) > 0:
            print("Could not find SHRDS Pilot matches for:")
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
