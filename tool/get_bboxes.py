# check a parquet generated has all 8000 subgrid boxes
import pandas as pd
import geopandas as gpd
filename = 'subbox_temps_anomalies.parquet'
# Import the grid8k function from eqarea.py located in the ../steps directory
import sys
sys.path.append('./steps')
from eqarea import grid8k

# Use the grid8k() function to get the bounding boxes (lat_min, lat_max, lon_min, lon_max)
bboxes = grid8k()

# Convert the list of bboxes to a DataFrame
gdf = pd.DataFrame(bboxes, columns=['lat_min', 'lat_max', 'lon_min', 'lon_max'])

# Round the coordinates to 2 decimal places
rounded_bboxes = gdf.round(2)

# Load the parquet file
df_parquet = pd.read_parquet(filename)

# Extract unique bounding boxes from the parquet DataFrame
parquet_bboxes = df_parquet[['lat_min', 'lat_max', 'lon_min', 'lon_max']].drop_duplicates()

# Round these coordinates to 2 decimal places
parquet_bboxes = parquet_bboxes.round(2)

# Convert the rounded parquet bboxes to tuples for easy comparison
parquet_bboxes_tuples = set(parquet_bboxes.apply(tuple, axis=1))

# Convert the rounded generated grid boxes to tuples
generated_bboxes_tuples = set(rounded_bboxes.apply(tuple, axis=1))

# Find the intersection of both sets to get matching bboxes
matching_bboxes = parquet_bboxes_tuples.intersection(generated_bboxes_tuples)

# Count the number of matches and the total number of unique bboxes in the parquet
num_matching_bboxes = len(matching_bboxes)
num_unique_parquet_bboxes = len(parquet_bboxes_tuples)

# Print the results
print(f"Number of matching bboxes: {num_matching_bboxes}")
print(f"Total number of unique bboxes in the parquet: {num_unique_parquet_bboxes}")
