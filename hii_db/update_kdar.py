"""
update_kdar.py

Update Catalog table kinematic distance ambiguity resolutions with
preliminary IGPS and ATCA data.

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

2023-02-09 Trey V. Wenger
"""

import argparse
import sqlite3
import pandas as pd


def main(db, fnames):
    """
    Update Catalog table in database with kinematic distance ambiguities from
    a set of formatted datasets.

    Inputs:
        db :: string
            Database filename
        fnames :: list of strings
            Files containing new KDAR data. Format must be space-delimited
            with, at minimum, headers 'gname' and 'KDAR'

    Returns: Nothing
    """
    # Load data
    merged_data = None
    for fname in fnames:
        # read data
        data = pd.read_csv(
            fname, delim_whitespace=True, comment="#", usecols=["gname", "KDAR", "QF"]
        )
        data = data.assign(dist_method=fname, dist_author="Wenger")

        # drop missing data
        bad = data["KDAR"] == "?"
        data = data.drop(data[bad].index)

        # numerate kdar QF
        data.loc[data["QF"] == "OG", "QF"] = "1"
        data.loc[data["QF"] == "A", "QF"] = "1"
        data.loc[data["QF"] == "B", "QF"] = "2"
        data.loc[data["QF"] == "C", "QF"] = "3"
        data.loc[data["QF"] == "T", "QF"] = "2"
        data["QF"] = pd.to_numeric(data["QF"])

        if merged_data is None:
            merged_data = data.copy()
        else:
            merged_data = pd.concat([merged_data, data], ignore_index=True)

    # Check conflicts
    print("KDAR Conflicts:")
    groups = merged_data.groupby("gname")
    for group in groups:
        if (group[1]["KDAR"].values != group[1]["KDAR"].iloc[0]).any():
            print(group[0])
            print(group[1])
            print()

    # Keep best KDAR
    merged_data = merged_data.loc[
        merged_data.groupby("gname")["QF"].idxmin()
    ].reset_index(drop=True)

    # Update Catalog
    rows = [
        [dat[1]["KDAR"], dat[1]["dist_method"], dat[1]["dist_author"], dat[1]["gname"]]
        for dat in merged_data.iterrows()
    ]

    # Add data to database
    print("Updating Catalog table with new KDARs...")
    with sqlite3.connect(db) as conn:
        cur = conn.cursor()
        cur.execute("PRAGMA foreign_keys = ON")
        cur.executemany(
            """
            UPDATE Catalog
            SET kdar=?, dist_method=?, dist_author=?
            WHERE gname=?
            """,
            rows,
        )
    print("Done!")
    print()


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(
        description="HII Region Database Update KDARs",
        prog="update_kdar.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    PARSER.add_argument("db", type=str, help="Database filename")
    PARSER.add_argument(
        "--files",
        nargs="+",
        help="KDAR data files",
    )
    ARGS = vars(PARSER.parse_args())
    main(ARGS["db"], ARGS["files"])
