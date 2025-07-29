"""
utils.py

General utilities.

Copyright(C) 2020-2025 by
Trey V. Wenger; tvwenger@gmail.com
L. D. Anderson;
This code is licensed under MIT license (see LICENSE for details)
"""

import os
import sqlite3
import pandas as pd
import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u


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
        cur.execute("DROP TABLE IF EXISTS CatalogParallax")
        cur.execute("DROP TABLE IF EXISTS Parallax")
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


def match_by_name(
    df1_0, df2_0, name1="Name", name2="Name", extension="2", order_by="Year"
):
    """
    Matches rows from df1_0 to df2_0 based on name1 and name2, adding columns with an extension

    Parameters:
        df1 (pd.DataFrame): The main dataframe.
        df2 (pd.DataFrame): The dataframe containing the rows to match.
        name1 (str): The name of df1 on which to match
        name2 (str): The name of df2 on which to match
        extension (str): The string to append to new column names in matched_rows_df.
        order_by (str): The column by which to order the matches (e.g., "Year"). If None, no ordering is done.

    Returns:
        matched_rows_df (pd.DataFrame): A dataframe with matched row data.
        unmatched_rows_df (pd.DataFrame): A dataframe with rows from df2 that did not match.
    """
    # Initialize lists for matched and unmatched rows
    matched_rows = []
    unmatched_rows = []

    # copy dfs so the originals aren't modified
    df1 = df1_0.copy()
    df2 = df2_0.copy()

    # Split 'Name' column by semicolons and replace NaN/None with empty list
    df1[f"{name1}_Split"] = (
        df1[name1].str.split(";").apply(lambda x: x if isinstance(x, list) else [])
    )
    df2[f"{name2}_Split"] = (
        df2[name2].str.split(";").apply(lambda x: x if isinstance(x, list) else [])
    )

    # rename
    common_columns = df1.columns.intersection(df2.columns)
    df2 = df2.rename(
        columns={
            col: f"{col}_{extension}"
            for col in list(set(list(common_columns) + [f"{name2}_Split"]))
        }
    )

    # Iterate over each row in df1
    for i, row in df1.iterrows():
        df1_names_set = set(row[f"{name1}_Split"])

        # Filter df2 rows where Name_Split matches df1 names
        matching_rows = df2[
            df2[f"{name2}_Split" + "_" + extension].apply(
                lambda x: bool(df1_names_set & set(x))
            )
        ]

        # If there are matching rows, process the first match
        if not matching_rows.empty:
            if order_by and order_by in matching_rows.columns:
                matching_rows = matching_rows.sort_values(by=order_by, ascending=False)
            matched_row = matching_rows.iloc[0]

            # Use pd.concat() to merge row and matched_row (add suffix to columns of matched_row)
            matched_row_info = pd.concat([row, matched_row])
            matched_rows.append(matched_row_info)
        else:
            # If no match, keep the original row
            matched_rows.append(row)
            unmatched_rows.append(row)

    # Return the matched and unmatched rows as a DataFrame
    matched_rows_df = pd.DataFrame(matched_rows).reset_index(drop=True)
    unmatched_rows_df = pd.DataFrame(unmatched_rows).reset_index(drop=True)
    return (matched_rows_df, unmatched_rows_df)


def match_by_distance(df1, df2, size=0.0, extension="2", order_by="Year"):
    """
    Matches rows from df1 to df2 based on distance criteria and returns a list of matched rows.

    Parameters:
        df1 (pd.DataFrame): The main dataframe.  Assumes GLong and GLat columns.
        df2 (pd.DataFrame): The dataframe containing the rows to match.  Assumes GLong and GLat columns.
        size (float): The distance threshold in degrees to consider a match.
        extension (str): The string to append to new column names in matched_rows_df.
        order_by (str): The column by which to order the matches (e.g., "Year"). If None, no ordering is done.

    Returns:
        matched_rows_df (pd.DataFrame): A dataframe with matched row data.
        unmatched_rows_df (pd.DataFrame): A dataframe with rows from df2 that did not match.
    """
    # Initialize lists for matched and unmatched rows
    matched_rows = []
    unmatched_rows = []

    # Compute distances between df1 and df2
    distances = np.sqrt(
        (
            np.asarray(df1["GLong"].values)[:, np.newaxis]
            - np.asarray(df2["GLong"].values)[np.newaxis, :]
        )
        ** 2
        + (
            np.asarray(df1["GLat"].values)[:, np.newaxis]
            - np.asarray(df2["GLat"].values)[np.newaxis, :]
        )
        ** 2
    )

    # Determine where distances are within the input size
    wh = distances < size

    # Iterate over each row in df1
    for i, row in df1.iterrows():
        # If there are matches within the size threshold
        if np.sum(wh[i, :]):
            matched_indices = np.where(wh[i, :])[0]  # Get the indices of the matches
            matching_rows = df2.iloc[matched_indices]

            # If order_by column exists, sort the matched rows by that column
            if order_by and order_by in df2.columns:
                matching_rows = matching_rows.sort_values(by=order_by, ascending=False)
            matched_row = matching_rows.iloc[0]  # Pick the first match

            # Calculate the distance to the first match
            match_distance = distances[i, matched_indices[0]]

            # Use pd.concat() to merge row and matched_row (add suffix to columns of matched_row)
            matched_row_info = pd.concat([row, matched_row.add_suffix(f"_{extension}")])
            matched_rows.append(matched_row_info)

            # Append the distance of the match to matched_row_info
            matched_row_info["Match_Distance"] = match_distance
        else:
            # If no match, keep the original row
            matched_rows.append(row)
            unmatched_rows.append(row)

    # Return the matched and unmatched rows as a DataFrame
    matched_rows_df = pd.DataFrame(matched_rows).reset_index(drop=True)
    unmatched_rows_df = pd.DataFrame(unmatched_rows).reset_index(drop=True)
    return (matched_rows_df, unmatched_rows_df)


def make_df(dir, columns, dtypes):
    """
    Makes large df from all .pkl files in "dir", with columns "columns"

    Parameters:
        dir (str): The directory to pull from.
        columns (list): The columns to include in the final datafram.
        dtypes (dict): Dictionary mapping columns to data types
    Returns:
        df (pd.DataFrame): A dataframe containing all the information from the directory.
    """

    df = pd.DataFrame()

    # Load each .pkl file into the dataframe
    for filename in os.listdir(dir):
        if filename.endswith(".pkl") and not filename.endswith("_2.pkl"):
            print(filename)
            file_path = os.path.join(dir, filename)

            # Load the .pkl file into a pandas DataFrame and include only specified columns
            this_df = pd.read_pickle(file_path)[columns]

            # Replace empty with NaN
            this_df = this_df.replace(r"^\s*$", np.nan, regex=True).astype(dtypes)

            # Strip units by extracting the 'value' attribute from each column
            this_df = this_df.map(lambda x: x.value if isinstance(x, u.Quantity) else x)
            df = pd.concat([df, this_df], ignore_index=True)
    return df.reset_index(drop=True)


def add_region(df, **kwargs):
    """
    Adds a new region to the given df using provided information
    in `kwargs`.  Computes RA, Dec, and GName.

    Paramaters:
    df (pd.DataFrame): The dataframe to which the new region data will be added.
    kwargs (dict): A dictionary containing the necessary data for the new region.
    """

    # Add the new row to the DataFrame using the provided kwargs
    df.loc[len(df)] = kwargs

    # ICRS
    galactic_coord = SkyCoord(
        l=kwargs["GLong"], b=kwargs["GLat"], frame="galactic", unit=(u.deg, u.deg)
    )
    icrs_coord = galactic_coord.transform_to("icrs")
    df["RA (J2000)"] = icrs_coord.ra.value
    df["Dec (J2000)"] = icrs_coord.dec.value
    df["GName"] = generate_gname(df["GLong"], df["GLat"])


def generate_gname(glong, glat):
    """
    Creates a GName of the form Glll.lll+/-bb.bbb

    Parameters:
    glong (float): The Galactic longitude
    glat (float): The Galactic latitude

    Returns:
    gname: GName of the form Glll.lll+/-bb.bbb
    """
    glong_str = glong.apply(lambda x: f"{np.floor(x * 1000) / 1000.:07.3f}")
    glat_str = glat.apply(lambda x: f"{np.floor(x * 1000) / 1000.:+07.3f}")

    # Combine them into the desired format
    return "G" + glong_str + glat_str
