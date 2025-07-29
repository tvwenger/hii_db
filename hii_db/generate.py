"""
generate.py

Generate the HII region database, following the schema described
in schema.txt.

Copyright(C) 2020-2025 by
Trey V. Wenger; tvwenger@gmail.com
L. D. Anderson;
This code is licensed under MIT license (see LICENSE for details)
"""

import argparse
from hii_db import utils, wise
from hii_db import balser_te_2015, wenger_te_2019
from hii_db import brown_shrds_pilot, wenger_shrds_2019, wenger_shrds_2021


def main(db, wisefile, wise_only=False, data_dir="data"):
    # Reset the database
    utils.reset(db, wise_only=wise_only)

    # Add WISE data to Catalog, Groups, and Detections tables
    wise.gen_catalog(db, wisefile)
    wise.gen_groups(db, wisefile)
    wise.add_detections(db, wisefile)

    if not wise_only:
        # Add Balser+2015 Detections
        balser_te_2015.add_detections(db, "140 Foot", data_dir=data_dir)
        balser_te_2015.add_detections(db, "GBT", data_dir=data_dir)

        # Add VLA Te Detections
        wenger_te_2019.add_detections(db, data_dir=data_dir)

        # Add SHRDS Detections
        brown_shrds_pilot.add_detections(db, data_dir=data_dir)
        wenger_shrds_2019.add_detections(db, data_dir=data_dir)
        wenger_shrds_2021.add_detections(db, data_dir=data_dir)


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
    PARSER.add_argument(
        "--data_dir",
        type=str,
        default="data/",
        help="Path to data directory",
    )
    ARGS = vars(PARSER.parse_args())
    main(
        ARGS["db"], ARGS["wise"], wise_only=ARGS["wiseonly"], data_dir=ARGS["data_dir"]
    )
