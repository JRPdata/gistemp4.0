# takes a SBBX txt output and uses an entire year to generate an ocean map of what grid boxes is covered
# Hansen2010 (Fig13, bottom right) discusses an all year ice free ocean map, but the 8000 subgrid boxes in the SBBX mask changes from month to month
# this outputs inclusively what is covered for a whole year for each grid box
# run to get the txt file first: python3 tool/SBBX_to_txt.py tmp/input/SBBX.ERSSTv5
filename = 'tmp/input/SBBX.ERSSTv5.txt'
import geojson
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, box

import sys

# Get the year from the command-line argument
if len(sys.argv) != 2:
    print("Usage: python program.py <year>")
    sys.exit(1)
year = int(sys.argv[1])

# Create a Plate Carree map
fig = plt.figure(figsize=(20, 10))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax.set_global()

# Add map features
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')
ax.add_feature(cfeature.COASTLINE, edgecolor='black')
ax.add_feature(cfeature.BORDERS, edgecolor='black')

# Read the original text file and create a list of polygons
polygons = []
with open(filename, 'r') as f:
    # Skip the header lines
    next(f)
    next(f)

    # Iterate over the lines
    for line in f:
        # Split the line into columns
        cols = line.split()

        # Check if the year matches
        if int(cols[4]) == year:
            # Extract latitude and longitude values
            # Extract latitude and longitude values
            lat1 = float(cols[0])
            lat2 = float(cols[1])
            lon1 = float(cols[2])
            lon2 = float(cols[3])

            # Create a polygon from the rectangle
            poly = box(lon1, lat1, lon2, lat2)
            polygons.append(poly)

# Combine the polygons
combined_poly = polygons[0]
for i, poly in enumerate(polygons[1:]):
    combined_poly = combined_poly.union(poly)
    print(f"Combined {i+1} of {len(polygons)-1} polygons")

# Add the combined polygon to the map
ax.add_geometries([combined_poly], ccrs.PlateCarree(), facecolor='lightblue', edgecolor='black')

# Save to file
plt.savefig(f'combined_grid_boxes_map_{year}.png', dpi=300)

import geojson

# ... (rest of the code remains the same)

# Save the combined polygon to a GeoJSON file
with open("ersst_sea_ice_mask.geojson", "w") as f:
    geojson.dump(geojson.MultiPolygon([[[list(zip(poly.exterior.coords.xy[0], poly.exterior.coords.xy[1]))] for poly in combined_poly.geoms]]), f)

import fiona

# ...

# Create a schema for the shapefile
schema = {
    'geometry': 'Polygon',
    'properties': {
        'id': 'int'
    }
}
import fiona

with fiona.open("ersst_sea_ice_mask.shp", "w", driver="ESRI Shapefile", schema=schema, crs="EPSG:4326") as writer:
    for i, poly in enumerate(combined_poly.geoms):
        writer.write({
            'geometry': poly.__geo_interface__,
            'properties': {'id': i}
        })
