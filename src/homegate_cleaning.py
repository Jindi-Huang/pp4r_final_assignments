import pandas as pd
import os
import numpy as np
import sys

# Get the grandparent directory of the script's directory
script_path = os.path.abspath(sys.argv[0])
grandparent_dir = os.path.dirname(os.path.dirname(os.path.dirname(script_path)))
os.chdir(grandparent_dir) # change this to your path

def data_cleaning(input_path, output_path):
    """Clean the data and export
    
    Args:
        input_path (str): Path to the raw data (Homegate_data.csv)
        output_path (str): Path to the processed data
    """
    df = pd.read_csv(input_path)

    df.rename(columns=lambda x: x.rstrip(':'), inplace=True)

    new_column_names = [
        'Net Rent (CHF)',
        'Additional Expenses (CHF)',
        'Gross Rent (CHF)',
        'Available Date',
        'Property Type',
        'Number of Rooms',
        'Number of Toilets',
        'Floor Level',
        'Living Area (sq. m.)',
        'Year Built',
        'Listing ID',
        'Object Reference',
        'Address',
        'URL Link',
        'Room Height (m.)',
        'Last Refurbishment Year',
        'Number of Floors',
        'Floor Space (sq. m.)',
        'Land Area (sq. m.)',
        'Number of Apartments',
        'Property Volume (cubic m.)',
    ]

    df.columns = new_column_names

    pattern = r'^(?P<Street>.*?),\s*(?P<Postcode>\d+)\s+(?P<City>.*)$'
    address_df = df['Address'].str.extract(pattern)
    df = pd.concat([df, address_df], axis=1)
    df = df.drop(columns=['Address'])
    df['City'] = df['City'].replace("Zürich, 8050 Zürich", "Zürich")

    columns_to_convert = ['Net Rent (CHF)',
        'Additional Expenses (CHF)',
        'Gross Rent (CHF)',]
    for column in columns_to_convert:
        df[column] = df[column].replace('On request', np.nan) 
        df[column] = df[column].str.replace('CHF', '')  
        df[column] = df[column].str.replace(',', '')   
        df[column] = df[column].str.replace('.–', '')  
    df[columns_to_convert] = df[columns_to_convert].apply(pd.to_numeric, errors='coerce').astype('Int64')
    
    df["Living Area (sq. m.)"] = df["Living Area (sq. m.)"].str.replace(' m2', '')  # Remove " m2"
    df['Floor Space (sq. m.)'] = df['Floor Space (sq. m.)'].str.replace(' m2', '')  
    df['Land Area (sq. m.)'] = df['Land Area (sq. m.)'].str.replace(' m2', '')
    df['Property Volume (cubic m.)'] = df['Property Volume (cubic m.)'].str.replace(' m3', '')
    df["Room Height (m.)"] = df["Room Height (m.)"].str.replace(' m', '') # Remove " m"

    columns_to_convert = ['Living Area (sq. m.)',
        'Room Height (m.)',
        'Number of Rooms',
        'Number of Toilets',
        'Floor Level',
        'Year Built']
    for column in columns_to_convert:
        df[column] = pd.to_numeric(df[column], errors='coerce')

    columns_to_convert = ['Year Built', 'Last Refurbishment Year', 'Postcode']
    for column in columns_to_convert:
        df[column] = pd.to_numeric(df[column], errors='coerce').astype('Int64')

    df['City'] = df['City'].replace("Zurich", "Zürich")
    df['City'] = df['City'].replace("Zürich Wiedikon", "Zürich")
    df['City'] = df['City'].replace("ZÃ¼rich", "Zürich")
    df['City'] = df['City'].replace("Zürich-Oerlikon", "Zürich")
    df['City'] = df['City'].replace("Zurigo", "Zürich")
    df['City'] = df['City'].replace("Glattpark (Opfikon)", "Opfikon")
    df['City'] = df['City'].replace("Glattpark", "Opfikon")
    df['City'] = df['City'].replace("Kilchberg ZH", "Kilchberg")

    df.to_csv(output_path, index=False)

if __name__ == "__main__":
    print
    data_cleaning(
        snakemake.input[0],
        snakemake.output[0]
    )