"""
add_parallax.py

Add parallax table to database.

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

import argparse
import sqlite3
import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u


def reset(db):
    """
    Reset database tables.

    Inputs:
      db :: string
        The database filename

    Returns: Nothing
    """
    print("Resetting Parallax and CatalogParallax tables...")
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")

        # Delete table if exists
        cur.execute("DROP TABLE IF EXISTS CatalogParallax")
        cur.execute("DROP TABLE IF EXISTS Parallax")

        # Create table
        cur.execute(
            """
        CREATE TABLE Parallax
        (id integer primary key autoincrement,
        gname string,
        alias string,
        ra real,
        dec real,
        glong real,
        glat real,
        plx real,
        e_plx real,
        mux real,
        e_mux real,
        muy real,
        e_muy real,
        vlsr real,
        e_vlsr real,
        author text)
        """
        )
        cur.execute(
            """
        CREATE TABLE CatalogParallax
        (catalog_id int,
        parallax_id int,
        separation float,
        FOREIGN KEY(catalog_id) REFERENCES Catalog(id),
        FOREIGN KEY(parallax_id) REFERENCES Parallax(id))
        """
        )
    print("Done.")
    print()


def add_sources(db, datafile, reffile):
    """
    Add sources to table.

    Inputs:
        db :: string
            Database filename
        datafile :: string
            Parallax data filename
        reffile :: string
            Parallax references filename

    Returns: Nothing
    """
    # Reid Reid+2019 data, get glong and glat
    print("Adding Reid+2019 Parallax sources...")
    data = np.genfromtxt(datafile, dtype=None, names=True, encoding="UTF-8")
    coords = SkyCoord(
        data["ra"],
        data["dec"],
        frame="fk5",
        unit=(u.hourangle, u.deg),
        pm_ra_cosdec=data["mux"] * u.mas / u.yr,
        pm_dec=data["muy"] * u.mas / u.yr,
    )

    # Read references and match with data
    refs = np.genfromtxt(
        reffile, dtype=None, names=True, delimiter=";", autostrip=True, encoding="UTF-8"
    )
    match_refs = []
    for dat in data:
        dat_refs = dat["ref"].split(",")
        match_refs.append(
            [
                refs[refs["id"] == int(dat_ref)]["author"][0]
                if dat_ref != "private"
                else "Private Communication"
                for dat_ref in dat_refs
            ]
        )

    # Populate rows
    rows = [
        [
            dat["gname"],
            dat["alias"] if dat["alias"] != "None" else None,
            coord.fk5.ra.deg,
            coord.fk5.dec.deg,
            coord.galactic.l.deg,
            coord.galactic.b.deg,
            dat["plx"],
            dat["e_plx"],
            dat["mux"],
            dat["e_mux"],
            dat["muy"],
            dat["e_muy"],
            float(dat["vlsr"]),
            float(dat["e_vlsr"]),
            ";".join(match_ref),
        ]
        for dat, coord, match_ref in zip(data, coords, match_refs)
    ]

    # Add data to database
    print("Populating Parallax table...")
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        cur.executemany(
            """
        INSERT INTO Parallax
        (gname, alias, ra, dec, glong, glat, plx, e_plx, mux, e_mux,
        muy, e_muy, vlsr, e_vlsr, author)
        VALUES
        (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """,
            rows,
        )
    print("Done!")
    print()

    # Match Parallax sources to WISE Catalog
    print("Matching Parallax sources to WISE Catalog...")
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

        # Get the parallax data
        cur.execute("SELECT id, gname, ra, dec FROM Parallax")
        dets = np.array(
            cur.fetchall(),
            dtype=[("id", "i"), ("gname", "U100"), ("ra", "f8"), ("dec", "f8")],
        )
        det_coords = SkyCoord(dets["ra"], dets["dec"], frame="fk5", unit="deg")

    # Loop over catalog and find matches
    rows = []
    for cat, coord in zip(wise, wise_coords):
        # Match based on GName
        parallax_id = None
        match = np.where(dets["gname"] == cat["gname"])[0]
        if len(match) == 1:
            match = match[0]
            parallax_id = dets["id"][match]
            print(
                "Matching WISE {0} to parallax {1} based on GName".format(
                    cat["gname"], dets["gname"][match]
                )
            )

        # If that doesn't work, matched based on alias
        if parallax_id is None:
            match = np.where(cat["alias"] == dets["gname"])[0]
            if len(match) == 1:
                match = match[0]
                parallax_id = dets["id"][match]
                print(
                    "Matching WISE {0} to parallax {1} based on alias".format(
                        cat["gname"], dets["gname"][match]
                    )
                )

        # If that doesn't work, match based on position within WISE region radius
        if parallax_id is None:
            seps = coord.separation(det_coords).arcsec
            # within WISE region
            min_sepind = np.min(seps - cat["radius"])
            if min_sepind < 0.0:
                match = np.argmin(seps - cat["radius"])
                parallax_id = dets["id"][match]
                print(
                    "Matching WISE {0} to parallax {1} based on WISE size".format(
                        cat["gname"], dets["gname"][match]
                    )
                )

        # If that doesn't work, match based on position within 90 arcsec
        if parallax_id is None:
            seps = coord.separation(det_coords).arcsec
            # within WISE region
            min_sepind = np.min(seps)
            if min_sepind < 90.0:
                match = np.argmin(seps)
                parallax_id = dets["id"][match]
                print(
                    "Matching WISE {0} to parallax {1} based on "
                    "separation ({2:.1f} arcsec)".format(
                        cat["gname"], dets["gname"][match], np.min(seps)
                    )
                )

        # No match :(
        if parallax_id is None:
            continue

        # Compute separation to Catalog match
        sep = coord.separation(det_coords[match]).arcsec
        rows.append((int(cat["id"]), int(dets["id"][match]), sep))

    # Populate CatalogParallax
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        cur.executemany(
            """
        INSERT INTO CatalogParallax
        (catalog_id, parallax_id, separation) VALUES (?, ?, ?)
        """,
            rows,
        )
    print("Done!")
    print()


def main(
    db,
    datafile="data/reid_2019/reid2019_merge.txt",
    reffile="data/reid_2019/reid2019_refs.txt",
):
    """
    Reset and populate Parallax table.

    Inputs:
        db :: string
            Database filename
        datafile :: string
            Parallax data filename
        reffile :: string
            Parallax references filename

    Returns: Nothing
    """
    # Reset database
    reset(db)

    # Add sources and match to WISE Catalog
    add_sources(db, datafile, reffile)


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(
        description="HII Region Database Parallax Table Generator",
        prog="add_parallax.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    PARSER.add_argument("db", type=str, help="The database filename")
    PARSER.add_argument(
        "--data",
        type=str,
        default="data/reid_2019/reid2019_merge.txt",
        help="The parallax data filename",
    )
    PARSER.add_argument(
        "--refs",
        type=str,
        default="data/reid_2019/reid2019_refs.txt",
        help="The parallax data references filename",
    )
    ARGS = vars(PARSER.parse_args())
    main(ARGS["db"], datafile=ARGS["data"], reffile=ARGS["refs"])
