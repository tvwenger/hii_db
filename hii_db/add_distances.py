"""
add_distances.py

Create kinematic distance table to database.

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
from kd import pdf_kd
import pickle


def get_data(db):
    """
    Get all detections with VLSR

    Inputs:
        db :: string
            Database filename

    Returns:
        newdata :: numpy record array
            Relevant data with median VLSR and e_VLSR for every
            WISE Catalog source with a detection.
    """
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        cur.execute(
            """
        SELECT cat.id, cat.gname, cat.glong, cat.glat,
        det.component, det.vlsr, det.e_vlsr
        FROM Catalog cat
        INNER JOIN CatalogDetections catdet ON catdet.catalog_id = cat.id
        INNER JOIN Detections det on catdet.detection_id = det.id
        WHERE det.vlsr IS NOT NULL AND
        (det.source = "WISE Catalog" OR
        (det.source = "SHRDS Full Catalog" AND det.lines = "H88-H112" AND
         det.line_qf < 5) OR
        (det.source = "VLA Te" AND det.lines = "H87-H93"))
        """
        )
        data = np.array(
            cur.fetchall(),
            dtype=[
                ("id", "i4"),
                ("gname", "U100"),
                ("glong", "f8"),
                ("glat", "f8"),
                ("component", "U4"),
                ("vlsr", "f8"),
                ("e_vlsr", "f8"),
            ],
        )
        # fix nan e_vlsr
        data["e_vlsr"][np.isnan(data["e_vlsr"])] = 0.1

    # Remove any source that has one or more multiple velocity components
    bad_sources = [row["gname"] for row in data if row["component"] != "None"]
    bad_sources = np.unique(bad_sources)
    bad_rows = [i for i, row in enumerate(data) if row["gname"] in bad_sources]
    data = np.delete(data, bad_rows)

    # Get median VLSR per unqiue WISE Catalog source
    gnames, unique_idx = np.unique(data["gname"], return_index=True)
    newdata = data[unique_idx]
    for i, gname in enumerate(gnames):
        # update VLSR and e_VLSR to median
        match = data["gname"] == gname
        newdata[i]["vlsr"] = np.median(data[match]["vlsr"])
        newdata[i]["e_vlsr"] = np.median(data[match]["e_vlsr"])
    return newdata


def compute_distances(
    data,
    db,
    num_samples=5000,
    batchsize=100,
    rotcurve="reid2019_rotcurve",
    tablename="Distances_Reid2019",
):
    """
    Compute the Monte Carlo Reid+2019 kinematic distances.
    Reset database table and populate new table.

    Inputs:
        data :: numpy record array
            VLSR data returned by get_data()
        db :: string
            Database filename
        num_samples :: integer
            Number of Monte Carlo samples to generate
        batchsize :: integer
            Batch size
        rotcurve :: string
            Rotation curve
        tablename :: string
            Name for new table

    Returns: Nothing
    """
    print("Resetting {0} table...".format(tablename))
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")

        # Delete table if exists
        cur.execute("DROP TABLE IF EXISTS {0}".format(tablename))

        # Create table
        cur.execute(
            """
        CREATE TABLE {0}
        (catalog_id integer,
        Rgal real,
        Rgal_err_neg real,
        Rgal_err_pos real,
        Rgal_kde blob,
        Rtan real,
        Rtan_err_neg real,
        Rtan_err_pos real,
        Rtan_kde blob,
        near real,
        near_err_neg real,
        near_err_pos real,
        near_kde blob,
        far real,
        far_err_neg real,
        far_err_pos real,
        far_kde blob,
        distance_kde blob,
        tangent real,
        tangent_err_neg real,
        tangent_err_pos real,
        tangent_kde blob,
        vlsr_tangent real,
        vlsr_tangent_err_neg real,
        vlsr_tangent_err_pos real,
        vlsr_tangent_kde blob,
        FOREIGN KEY(catalog_id) REFERENCES Catalog(id))
        """.format(
                tablename
            )
        )
    print("Done.")
    print()

    # Compute kinematic distances in groups of 500
    rows = []
    for i in range(len(data) // batchsize + 1):
        start = i * batchsize
        end = (i + 1) * batchsize
        mydata = data[start:end]
        print("Computing kinematic distances (group {0})...".format(i))
        kdout = pdf_kd.pdf_kd(
            mydata["glong"],
            mydata["glat"],
            mydata["vlsr"],
            velo_err=mydata["e_vlsr"],
            rotcurve=rotcurve,
            num_samples=num_samples,
        )
        print("Done.")
        print()
        rows += [
            [
                int(mydata["id"][i]),
                kdout["Rgal"][i],
                kdout["Rgal_err_neg"][i],
                kdout["Rgal_err_pos"][i],
                sqlite3.Binary(pickle.dumps(kdout["Rgal_kde"][i])),
                kdout["Rtan"][i],
                kdout["Rtan_err_neg"][i],
                kdout["Rtan_err_pos"][i],
                sqlite3.Binary(pickle.dumps(kdout["Rtan_kde"][i])),
                kdout["near"][i],
                kdout["near_err_neg"][i],
                kdout["near_err_pos"][i],
                sqlite3.Binary(pickle.dumps(kdout["near_kde"][i])),
                kdout["far"][i],
                kdout["far_err_neg"][i],
                kdout["far_err_pos"][i],
                sqlite3.Binary(pickle.dumps(kdout["far_kde"][i])),
                sqlite3.Binary(pickle.dumps(kdout["distance_kde"][i])),
                kdout["tangent"][i],
                kdout["tangent_err_neg"][i],
                kdout["tangent_err_pos"][i],
                sqlite3.Binary(pickle.dumps(kdout["tangent_kde"][i])),
                kdout["vlsr_tangent"][i],
                kdout["vlsr_tangent_err_neg"][i],
                kdout["vlsr_tangent_err_pos"][i],
                sqlite3.Binary(pickle.dumps(kdout["vlsr_tangent_kde"][i])),
            ]
            for i in range(len(mydata))
        ]

    # Save to database
    print("Populating {0} table...".format(tablename))
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        cur.executemany(
            """
        INSERT INTO {0} (catalog_id,
        Rgal, Rgal_err_neg, Rgal_err_pos, Rgal_kde,
        Rtan, Rtan_err_neg, Rtan_err_pos, Rtan_kde,
        near, near_err_neg, near_err_pos, near_kde,
        far, far_err_neg, far_err_pos, far_kde,
        distance_kde,
        tangent, tangent_err_neg, tangent_err_pos, tangent_kde,
        vlsr_tangent, vlsr_tangent_err_neg, vlsr_tangent_err_pos, vlsr_tangent_kde)
        VALUES
        (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        """.format(
                tablename
            ),
            rows,
        )
    print("Done!")
    print()


def main(
    db,
    num_samples=5000,
    batchsize=100,
    rotcurve="reid19_rotcurve",
    tablename="Distances_Reid2019",
):
    """
    Compute the Monte Carlo Reid+2019 kinematic distances.
    Reset database table and populate new table.

    Inputs:
        db :: string
            Database filename
        num_samples :: integer
            Number of Monte Carlo samples to generate
        batchsize :: integer
            Batch size
        rotcurve :: string
            Rotation curve
        tablename :: string
            Name for new table

    Returns: Nothing
    """
    # Get the data
    data = get_data(db)
    print("Found {0} unique Catalog sources with VLSR".format(len(data)))
    print()

    # Generate and add kinematic distances
    compute_distances(
        data,
        db,
        num_samples=num_samples,
        batchsize=batchsize,
        rotcurve=rotcurve,
        tablename=tablename,
    )


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(
        description="HII Region Database Distances Table Generator",
        prog="add_distances.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    PARSER.add_argument("db", type=str, help="Database filename")
    PARSER.add_argument(
        "-n",
        "--num_samples",
        type=int,
        default=5000,
        help="Number of Monte Carlo samples",
    )
    PARSER.add_argument(
        "-b",
        "--batchsize",
        type=int,
        default=100,
        help="Batch size",
    )
    PARSER.add_argument(
        "-r", "--rotcurve", type=str, default="reid19_rotcurve", help="Rotation curve"
    )
    PARSER.add_argument(
        "-t", "--tablename", type=str, default="Distances_Reid2019", help="Table name"
    )
    ARGS = vars(PARSER.parse_args())
    main(
        ARGS["db"],
        num_samples=ARGS["num_samples"],
        batchsize=ARGS["batchsize"],
        rotcurve=ARGS["rotcurve"],
        tablename=ARGS["tablename"],
    )
