# quick tool that outputs lat,lon range for filename SBBX txt output
filename = 'tmp/input/SBBX.ERSSTv5.orig.txt'

with open(filename, 'r') as f:
    # Skip the header lines
    next(f)
    next(f)

    # Initialize min/max variables
    min_lat = float('inf')
    max_lat = float('-inf')
    min_lon = float('inf')
    max_lon = float('-inf')

    # Iterate over the lines
    for line in f:
        # Split the line into columns
        cols = line.split()

        # Extract latitude and longitude values
        lat = float(cols[0])
        lon = float(cols[1])

        # Update min/max variables
        min_lat = min(min_lat, lat)
        max_lat = max(max_lat, lat)
        min_lon = min(min_lon, lon)
        max_lon = max(max_lon, lon)

# Print the extents
print(f'Latitude range: {min_lat} to {max_lat}')
print(f'Longitude range: {min_lon} to {max_lon}')
