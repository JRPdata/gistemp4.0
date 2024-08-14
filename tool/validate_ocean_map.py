# verify ocean_map.py works properly
import fiona
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import shape

# Read the shapefile
with fiona.open("ersst_sea_ice_mask.shp", "r") as src:
    geoms = [shape(poly["geometry"]) for poly in src]

# Create a new map to validate
fig = plt.figure(figsize=(20, 10))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

# Set the extent of the map
ax.set_global()

# Add map features
ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')
ax.add_feature(cfeature.COASTLINE, edgecolor='black')
ax.add_feature(cfeature.BORDERS, edgecolor='black')

# Add the mask to the map
ax.add_geometries(geoms, ccrs.PlateCarree(), facecolor='lightblue', edgecolor='black')

# Save the validation map to a file
plt.savefig("ersst_sea_ice_mask_validation.png", dpi=300)
