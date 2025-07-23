import os
import pandas as pd
import sqlite3

# Directory containing the .pkl files
pkl_directory = '/Users/loren/papers/wise/python/rrl_surveys/'

# SQLite database file (change to the appropriate path for your database)
db_file = 'rrl_surveys.db'

# Column names
columns = ["Name", "GName", "GLong", "GLat", "RA (J2000)", "Dec (J2000)", 
           "TL", "e_TL", "VLSR", "e_VLSR", "FWHM", "e_FWHM", "RMS"]

# Connect to the SQLite database (it will create the database if it doesn't exist)
conn = sqlite3.connect(db_file)
cursor = conn.cursor()

# Create the table with the specified columns (if it doesn't exist already)
create_table_query = f"""
CREATE TABLE IF NOT EXISTS catalog_data (
    Name TEXT,
    GName TEXT,
    GLong REAL,
    GLat REAL,
    "RA (J2000)" REAL,
    "Dec (J2000)" REAL,
    TL REAL,
    e_TL REAL,
    VLSR REAL,
    e_VLSR REAL,
    FWHM REAL,
    e_FWHM REAL,
    RMS REAL
);
"""

cursor.execute(create_table_query)
conn.commit()

print("Table created successfully (if it didn't exist already).")



# Function to insert data into the database
def insert_data(df):
    # Insert the DataFrame into the SQL database
    for row in df.itertuples(index=False, name=None):
        cursor.execute(f"""
        INSERT INTO catalog_data 
        (Name, GName, GLong, GLat, "RA (J2000)", "Dec (J2000)", TL, e_TL, VLSR, e_VLSR, FWHM, e_FWHM, RMS)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""", row)
    conn.commit()

# Loop through each .pkl file in the directory
for filename in os.listdir(pkl_directory):
    if filename.endswith('.pkl') and not filename.endswith('_2.pkl'):
        file_path = os.path.join(pkl_directory, filename)
        
        # Load the .pkl file into a pandas DataFrame
        df = pd.read_pickle(file_path)
        
        # Ensure the DataFrame contains the correct columns
        if all(col in df.columns for col in columns):
            print(f"Processing {filename}...")
            # Insert the DataFrame data into the SQL database
            insert_data(df)
        else:
            print(f"Skipping {filename}: Missing one or more required columns")

# Close the database connection
conn.close()
print("Database populated successfully!")
