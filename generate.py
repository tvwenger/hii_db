"""
generate.py

Generate the HII region database, following the schema described
in schema.txt.

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

import argparse
from hii_db import utils, wise
from hii_db import balser_te_2015, wenger_te_2019
from hii_db import brown_shrds_pilot, wenger_shrds_2019, wenger_shrds_2021


def main(db, wisefile, wise_only=False):
    # Reset the database
    utils.reset(db, wise_only=wise_only)

    # Add WISE data to Catalog, Groups, and Detections tables
    wise.gen_catalog(db, wisefile)
    wise.gen_groups(db, wisefile)
    wise.add_detections(db, wisefile)

    if not wise_only:
        # Add Balser+2015 Detections
        balser_te_2015.add_detections(db, "140 Foot")
        balser_te_2015.add_detections(db, "GBT")

        # Add VLA Te Detections
        wenger_te_2019.add_detections(db)

        # Add SHRDS Detections
        brown_shrds_pilot.add_detections(db)
        wenger_shrds_2019.add_detections(db)
        wenger_shrds_2021.add_detections(db)


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(
        description="HII Region Database Generator",
        prog="generate.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    PARSER.add_argument("db", type=str, help="The database filename")
    PARSER.add_argument(
        "--wise",
        type=str,
        default="wise/wise_hii_V3.0_hrds.csv",
        help="WISE Catalog CSV filename",
    )
    PARSER.add_argument(
        "--wiseonly",
        action="store_true",
        default=False,
        help="Generate database with WISE Catalog only",
    )
    ARGS = vars(PARSER.parse_args())
    main(ARGS["db"], ARGS["wise"], wise_only=ARGS["wiseonly"])
