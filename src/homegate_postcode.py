import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import folium as fl

import json
import matplotlib as mpl
from osgeo import gdal
import geopandas as gpd





def create_map(input_path, output_path):
    """Analyze the data and export the results

    Args:
        input_path (str): Path to the processed data (Homegate_data_cleaned.csv)
        output_path (str): Path to the output directory
    """

    df = pd.read_csv(input_path)

    # getting postcode boundaries from source: https://github.com/mikpan/ch-maps/blob/master/geo/ch-plz.geojson
    zh_2 = fl.Map(location=[47.3769, 8.5417], zoom_start=12, tiles="Stamen Terrain")

    with open('data/ch-plz.geojson', 'r') as geojson_file:
        zurich_postcodes_geojson_2 = json.load(geojson_file)
    fl.GeoJson(zurich_postcodes_geojson_2).add_to(zh_2)

    gdal.SetConfigOption('SHAPE_RESTORE_SHX', 'YES')
    my_dir = os.getcwd()
    in_vector = my_dir + "/data/PLZO_SHP_LV95" 
    gdf = gpd.read_file(in_vector)
    os.chdir(in_vector)

    # creating aggregated data for postcodes - here average prices
    aggr_df = df.groupby('Postcode')['Gross Rent (CHF)'].mean().reset_index()
    aggr_df['Postcode'] = aggr_df['Postcode'].astype(int)
    post_code_list = aggr_df["Postcode"].values
   
    # selecting the relevant indices
    indices = []
    for i in range(len(zurich_postcodes_geojson_2["features"])):
        if zurich_postcodes_geojson_2["features"][i]["id"] in post_code_list:
            indices.append(i)

    indices = np.array(indices)

    # combining them into a new geojson file
    zurich_boundaries = zurich_postcodes_geojson_2.copy()
    zurich_boundaries["features"] = [zurich_postcodes_geojson_2["features"][i] for i in indices]

    zh_boundaries = fl.Map(location=[47.3769, 8.5417], zoom_start=11, tiles="Stamen Terrain")
    fl.GeoJson(zurich_boundaries).add_to(zh_boundaries)

    # merging the geojson file with our own aggr_df
    for i in range(len(zurich_boundaries["features"])):
        zurich_boundaries["features"][i]["average"] = int(aggr_df.loc[aggr_df['Postcode'] == zurich_boundaries["features"][i]["id"], "Gross Rent (CHF)"].values)
    

    zh_boundaries = fl.Map(location=[47.3769, 8.5417], zoom_start=11, tiles="Stamen Terrain")

    for feature in zurich_boundaries["features"]:
        # for one area:
        my_opacity = (feature["average"]-1500) / 3500

        geojson_layer = fl.GeoJson(feature,
                style_function = lambda x, opacity=my_opacity: {
                'fillColor': "red",
                'fillOpacity': opacity,
                'color': 'black',
                'weight': 1
            }).add_to(zh_boundaries)
        
        geojson_layer.add_to(zh_boundaries)
    # saving the files
    os.chdir(my_dir)
    zh_boundaries.save(output_path)



    

# my_dir = '/Users/jindi/Documents/GitHub/pp4r_final_assignments' # change this to your path
# os.chdir(my_dir) 
    
# create_map("data/Homegate_data_cleaned.csv", "/Users/jindi/Documents/GitHub/pp4r_final_assignments/output/htmls/zurich_map.html")

if __name__ == "__main__":   
    create_map(
        snakemake.input[0], 
        snakemake.output[0]
    )