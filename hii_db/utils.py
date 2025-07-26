"""
utils.py

General utilities.

Copyright(C) 2020-2023 by
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
2021-09-30 Trey V. Wenger - reorganization
2023-05-24 Trey V. Wenger - Adjust schema to v5
"""

import os
import sqlite3
from astropy.coordinates import SkyCoord


def reset(db, wise_only=False):
    """
    Drop tables if they already exist, create new table schema.

    Inputs:
        db :: string
            Database filename
        wise_only :: boolean
            If True, skip Fields and FieldsDetections tables

    Returns: Nothing
    """
    print("Resetting database...")
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")

        # Delete tables if present (need to delete link tables first)
        cur.execute("DROP TABLE IF EXISTS CatalogDetections")
        cur.execute("DROP TABLE IF EXISTS CatalogGroups")
        cur.execute("DROP TABLE IF EXISTS Catalog")
        cur.execute("DROP TABLE IF EXISTS Groups")
        cur.execute("DROP TABLE IF EXISTS Detections")
        cur.execute("DROP TABLE IF EXISTS Fields")

        # WISE Catalog table
        cur.execute(
            """
        CREATE TABLE Catalog
        (id integer primary key autoincrement,
        gname text,
        alias text,
        hii_name text,
        catalog text,
        ra real,
        dec real,
        glong real,
        glat real,
        radius real,
        kdar text,
        dist_method text,
        dist_author text)
        """
        )

        # WISE Catalog Groups table
        cur.execute(
            """
        CREATE TABLE Groups
        (id integer primary key autoincrement,
        name text,
        vlsr real,
        e_vlsr real,
        kdar text)
        """
        )

        # Fields table
        cur.execute(
            """
        CREATE TABLE Fields
        (id integer primary key autoincrement,
        name text,
        ra real,
        dec real,
        glong real,
        glat real,
        hpbw real)
        """
        )

        # Detections table
        cur.execute(
            """
        CREATE TABLE Detections
        (id integer primary key autoincrement,
        name text,
        ra real,
        dec real,
        glong real,
        glat real,
        line_freq real,
        component text,
        line real,
        e_line real,
        line_unit text,
        vlsr real,
        e_vlsr real,
        fwhm real,
        e_fwhm real,
        spec_rms real,
        line_qf integer,
        line_snr real,
        cont_freq real,
        cont real,
        e_cont real,
        cont_unit text,
        area real,
        e_area real,
        area_unit text,
        cont_qf integer,
        pb_level real,
        linetocont real,
        e_linetocont real,
        te real,
        e_te real,
        lines text,
        beam_area real,
        telescope text,
        author text,
        source text,
        type text,
        taper text,
        field_id int,
        FOREIGN KEY(field_id) REFERENCES Fields(id))
        """
        )

        # WISE Catalog <=> Groups
        cur.execute(
            """
        CREATE TABLE CatalogGroups
        (catalog_id int,
        group_id int,
        FOREIGN KEY(catalog_id) REFERENCES Catalog(id),
        FOREIGN KEY(group_id) REFERENCES Groups(id))
        """
        )

        # WISE Catalog <=> Detections
        cur.execute(
            """
        CREATE TABLE CatalogDetections
        (catalog_id int,
        detection_id int,
        separation float,
        FOREIGN KEY(catalog_id) REFERENCES Catalog(id),
        FOREIGN KEY(detection_id) REFERENCES Detections(id))
        """
        )

    print("Done!")
    print()


def parse_region_coord(regionfile):
    """
    Parse a CASA region file to get the RA, Dec position. Return
    the position as an astropy coordinate object.

    Inputs:
      regionfile: string
        the filename of the region file

    Returns:
      coord: astropy.SkyCoord
        the coordinate
    """
    # Parse file
    if not os.path.exists(regionfile):
        raise ValueError("{0} does not exist!".format(regionfile))
    with open(regionfile, "r") as f:
        line = f.readline()
        line = f.readline()
        mypart = line.split("[[")[1].split("]")[0]
        rapart = mypart.split(",")[0]
        decpart = mypart.split(",")[1][1:]
        ra_h, ra_m, ra_s = rapart.split(":")
        dec_sn = decpart[0]
        dec_d, dec_m, dec_s, dec_ss = decpart[1:].split(".")

    # Create astropy coordinate
    mycoord = "{0}h{1}m{2}s {3}{4}d{5}m{6}.{7}s".format(
        ra_h, ra_m, ra_s, dec_sn, dec_d, dec_m, dec_s, dec_ss
    )
    coord = SkyCoord(mycoord, frame="fk5")
    return coord
