# gets ERSSTv5 data (PSL/NOAA sst.mnmean.nc file) into subboxes (multi-process regridding); outputs parquets (for reference) and SBBX used by GISTEMP.
# Run with pypy3 and latest numpy for significant speedup.
# Only global anomaly has been checked: 1 month is 0.01 C off
# Options: 3 different area methods possible, see process_batch() (lines 396-398)
# Requirement 1: place 'good' txt output of SBBX into orig_sst_txt_path (run tool/SBBX_to_txt.py on tmp/input/SBBX.ERSSTv5) (we require since we can't generate the subgrid boxes per month mask),
# Requirement 2: place (consolidated) ERSST5 (sst.mnmean.nc file) in ersst_monthly_means

import traceback
# delete subbox_temps.parquet if you have fresh ersst data
# performance testing # about 25-30% faster with pyp3 (conda)
import time
# not a package, but the fort.py file in gistemp's directory (tool/)
from fort import File

import struct
import pandas as pd
import glob
import shutil
import os

from datetime import datetime, timedelta
import netCDF4 as nc
import numpy as np
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.geometry import box
from shapely.strtree import STRtree
from tqdm import tqdm  # For progress bar
from pyproj import Proj, transform
from pyproj import Transformer
# Import the grid8k function from eqarea.py located in the ../steps directory
import sys
sys.path.append('steps')
from eqarea import grid8k

# performance
import concurrent
from concurrent.futures import ProcessPoolExecutor
# for tqdm multi-process
import multiprocessing

start_time = time.time()

# generate only when monthly means have updated (takes > hour with pypy3)
do_generate_subbox_parquets = True

# baseline for averaging subboxes
baseline_start_year = 1951
baseline_end_year = 1980
ersst_monthly_means = 'tmp/input/sst.mnmean.nc'
intermediate_path = 'intermediate_chunks'
# SBBX_to_txt.f output for the reference subboxes to include
# i.e.: run gistemp step 0, rename reference file to .orig after extracting, run SBBX_to_txt (fortran program)
orig_sst_txt_path = 'tmp/input/SBBX.ERSSTv5.orig.txt'
subbox_temps_parquet_path = 'subbox_temps.parquet'
subbox_temps_anomalies_parquet_path = 'subbox_temps_anomalies.parquet'
baseline_parquet_path = 'subbox_baseline_temps.parquet'
output_sbbx = 'tmp/input/SBBX.ERSSTv5'
# year to extract from sst orig.txt (from SBBX_to_txt.f)
reference_subbox_year = 2023
# Clip the SBBX data starting from more reliably records to match gistemp
start_clip_year = 1880
# particular to SBBX data (hardcoded in program)
missing_flag = 9999
# doesn't seem to generate more accurate data so we don't use
sea_surface_cutoff_temp = -1.77
# step 3 original infills all subboxes so long as there is any overlapping grid data not masked or meeting the surface cutoff temp


def find_masked_subboxes_in_orig_txt():
    global orig_sst_txt_path
    global baseline_start_year
    global missing_flag
    # Define the column names manually if the file lacks headers
    column_names = ['lat_min', 'lat_max', 'lon_min', 'lon_max', 'year', 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul',
                    'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    # Read the data, skipping the first two rows and using the defined column names
    df = pd.read_csv(orig_sst_txt_path, skiprows=2, sep=r'\s+', header=None, names=column_names)

    # Define the magic number
    magic_number = np.float32(missing_flag)

    # This includes temperature columns and possibly lat/lon columns
    df = df.astype({
        'lat_min': np.float32,
        'lat_max': np.float32,
        'lon_min': np.float32,
        'lon_max': np.float32,
        'year': np.int16,
        'Jan': np.float32,
        'Feb': np.float32,
        'Mar': np.float32,
        'Apr': np.float32,
        'May': np.float32,
        'Jun': np.float32,
        'Jul': np.float32,
        'Aug': np.float32,
        'Sep': np.float32,
        'Oct': np.float32,
        'Nov': np.float32,
        'Dec': np.float32
    })

    # Filter for year 1951
    df = df[df['year'] == np.int16(baseline_start_year)]

    # Drop the unnecessary 'year' column
    df = df.drop(columns=['year'])

    # Reset the index after filtering and dropping the column
    df.reset_index(drop=True, inplace=True)

    # Create a list of month columns
    months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    # Initialize an empty list to hold processed rows
    processed_rows = []
    processed_rows_unmasked = []

    # Iterate over each row to create the multi-index DataFrame
    for index, row in df.iterrows():
        for i, month in enumerate(months):
            # Determine if the value is equal to the magic number
            is_magic = row[month] == magic_number

            # Append the processed row only if is_magic is True
            if is_magic:
                processed_rows.append([
                    np.int16(round(row['lat_min'] * 100)),
                    np.int16(round(row['lat_max'] * 100)),
                    np.int16(round(row['lon_min'] * 100)),
                    np.int16(round(row['lon_max'] * 100)),
                    np.int8(i + 1)
                ])
            else:
                processed_rows_unmasked.append([
                    np.int16(round(row['lat_min'] * 100)),
                    np.int16(round(row['lat_max'] * 100)),
                    np.int16(round(row['lon_min'] * 100)),
                    np.int16(round(row['lon_max'] * 100)),
                    np.int8(i + 1)
                ])

    # Create the final DataFrame
    processed_df = pd.DataFrame(
        processed_rows,
        columns=['lat_min_100', 'lat_max_100', 'lon_min_100', 'lon_max_100', 'month']
    )
    processed_df_unmasked = pd.DataFrame(
        processed_rows_unmasked,
        columns=['lat_min_100', 'lat_max_100', 'lon_min_100', 'lon_max_100', 'month']
    )

    # Convert the latitude, longitude, and month columns to float32 and bool respectively
    processed_df = processed_df.astype({
        'lat_min_100': np.int16,
        'lat_max_100': np.int16,
        'lon_min_100': np.int16,
        'lon_max_100': np.int16,
        'month': np.int8
    })

    # Convert the latitude, longitude, and month columns to float32 and bool respectively
    processed_df_unmasked = processed_df_unmasked.astype({
        'lat_min_100': np.int16,
        'lat_max_100': np.int16,
        'lon_min_100': np.int16,
        'lon_max_100': np.int16,
        'month': np.int8
    })

    return processed_df, processed_df_unmasked


def generate_grid8k_polygons():
    # Generate the bounding boxes using grid8k()

    # Convert generator to list
    bboxes = list(grid8k())

    # Convert list to NumPy array for further processing
    bboxes_array = np.array(bboxes, dtype=np.float32)

    # Create polygons using the original float32 coordinates
    polygons = [box(lon_min, lat_min, lon_max, lat_max) for lat_min, lat_max, lon_min, lon_max in bboxes_array]

    # Multiply by 100, round, and convert to int16
    rounded_int100_bboxes = np.round(bboxes_array * 100).astype(np.int16)

    return polygons, rounded_int100_bboxes


def extract_time_tuples_from_netcdf(ds):
    """
    Extract (year, month) tuples from the time dimension of a NetCDF file using netCDF4.

    Parameters:
    - nc_file: Path to the NetCDF file.

    Returns:
    - A list of (year, month) tuples corresponding to the time dimension.
    """
    # Check if 'time' variable is present
    if 'time' not in ds.variables:
        raise ValueError("Time variable 'time' not found in the NetCDF file.")

    # Extract the time variable
    time_data = ds.variables['time']

    # Extract the time units and base date from the attributes
    units = time_data.units
    if 'days since' in units:
        base_date_str = units.split('days since ')[1]
        base_date = datetime.strptime(base_date_str, "%Y-%m-%d %H:%M:%S")
    else:
        raise ValueError("Unsupported time units format: {}".format(units))

    # Convert time values to datetime objects
    time_values = time_data[:]
    dates = [base_date + timedelta(days=int(day)) for day in time_values]

    # Extract year and month tuples
    time_tuples = [(np.int16(date.year), np.int8(date.month)) for date in dates]

    return time_tuples


def read_ersst_monthly_means_absolute(file_path):
    # Open the nc file
    ds = nc.Dataset(file_path, 'r')

    # Extract the temperature variable
    temp_var = ds.variables['sst']  # assuming the variable name is 'sst'

    # Get the temperature data
    temp_data = temp_var[:]

    # Get the latitude and longitude coordinates (unmasked)
    lats = ds.variables['lat'][:].data
    lons = ds.variables['lon'][:].data

    time_tuples = extract_time_tuples_from_netcdf(ds)

    # Close the nc file
    ds.close()

    return temp_data, lats, lons, time_tuples


def is_within_polygon(lon, lat, polygon):
    point = Point(lon, lat)
    return polygon.contains(point)


def create_spatial_index(lats, lons, resolution=2):
    """
    Create a spatial index for the regular grid cells.

    Parameters:
    - lats: 2D array of latitudes for the grid
    - lons: 2D array of longitudes for the grid
    - resolution: Grid resolution in degrees (default: 2)

    Returns:
    - str_tree: Spatial index (STRtree) for grid cells
    - index_list: List of tuples (lat_index, lon_index, Polygon)
    """
    grid_cells = []
    index_list = []

    res_div_2 = np.float32(resolution / 2)
    # Create bounding boxes for each grid cell, using signed lons for the polygon
    # Original data has lons from 0 .. 358

    # Remember: the spatial index uses polygons with signed lons, and the
    signed_offset = np.float32(-360.0)
    lats = np.array(lats).astype(np.float32)
    lons = np.array(lons).astype(np.float32)
    signed_lons = [(lon + signed_offset) if lon > 180 else lon for lon in lons]
    for lat_idx, lat in enumerate(lats):
        for lon_idx, lon in enumerate(signed_lons):
            min_lat = lat - res_div_2
            max_lat = lat + res_div_2
            min_lon = lon - res_div_2
            max_lon = lon + res_div_2

            lat_idx16 = np.int16(lat_idx)
            lon_idx16 = np.int16(lon_idx)

            # Handle longitude wrapping at the International Date Line
            if min_lon < -180.0:
                # Wraps around the western hemisphere
                bbox1 = Polygon([(-180.0, min_lat), (-180.0, max_lat),
                                 (max_lon, max_lat), (max_lon, min_lat)])
                bbox2 = Polygon([(min_lon + 360.0, min_lat),
                                 (min_lon + 360.0, max_lat),
                                 (180.0, max_lat), (180.0, min_lat)])
                grid_cells.extend([bbox1, bbox2])
                index_list.extend([(lat_idx16, lon_idx16, bbox1), (lat_idx16, lon_idx16, bbox2)])
            elif max_lon > 180.0:
                # Wraps around the eastern hemisphere
                bbox1 = Polygon([(min_lon, min_lat), (min_lon, max_lat),
                                 (180.0, max_lat), (180.0, min_lat)])
                bbox2 = Polygon([(-180.0, min_lat),
                                 (-180.0, max_lat),
                                 (max_lon - 360.0, max_lat),
                                 (max_lon - 360.0, min_lat)])
                grid_cells.extend([bbox1, bbox2])
                index_list.extend([(lat_idx16, lon_idx16, bbox1), (lat_idx16, lon_idx16, bbox2)])
            else:
                bbox = Polygon([(min_lon, min_lat), (min_lon, max_lat),
                                (max_lon, max_lat), (max_lon, min_lat)])
                grid_cells.append(bbox)
                index_list.append((lat_idx16, lon_idx16, bbox))

    # Create the spatial index
    str_tree = STRtree(grid_cells)

    return str_tree, index_list


# Use 6933 to preserve area
def latlon_to_projected(lats, lons, proj_epsg=6933):
    proj = Proj(init=f'epsg:{proj_epsg}')
    lons, lats = np.meshgrid(lons, lats)
    x, y = transform(Proj(init='epsg:6933'), proj, lons, lats)
    return x, y


# Function to transform polygon coordinates and calculate area
def calculate_area_in_projection(polygon, eq_transformer):
    # Transform coordinates
    if isinstance(polygon, Polygon) and polygon.area > 0:
        transformed_coords = [eq_transformer.transform(x, y) for x, y in polygon.exterior.coords]
        transformed_polygon = Polygon(transformed_coords)
        # Calculate area in the projected coordinate system
        return transformed_polygon.area
    else:
        return np.nan

# Don't do any transformation
def calculate_box_poly_area(polygon):
    if isinstance(polygon, Polygon) and polygon.area > 0:
        return polygon.area
    else:
        return np.nan


def calculate_spherical_polygon_area(polygon):
    """
    Calculate the area of a polygon using latitudes projected onto a sphere.

    Parameters:
        polygon (Polygon): A shapely Polygon object with latitude and longitude coordinates.

    Returns:
        float: The area of the polygon on the sphere.
    """
    if not isinstance(polygon, Polygon) or polygon.is_empty:
        return np.nan

    # Extract the bounding box of the polygon (min_lon, min_lat, max_lon, max_lat)
    min_lon, min_lat, max_lon, max_lat = polygon.bounds

    # Convert degrees to radians
    min_lon_rad = np.radians(min_lon)
    max_lon_rad = np.radians(max_lon)
    min_lat_rad = np.radians(min_lat)
    max_lat_rad = np.radians(max_lat)

    # Calculate the delta longitude
    delta_lon = max_lon_rad - min_lon_rad

    # Calculate the area using the formula: delta_lon * (sin(lat_max) - sin(lat_min))
    area = delta_lon * (np.sin(max_lat_rad) - np.sin(min_lat_rad))

    return np.abs(area)  # Return absolute area


def process_batch(batch_idx, polygon_batch, rounded_bbox_batch, temp_data, str_tree, index_list, time_tuples,
                    chunk_folder, masked_constant_type, cutoff_temp, start_chunk_num, progress_counter, lock):
    chunk_paths = []
    chunk_num = start_chunk_num + 1
    eq_transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    for poly_idx in range(len(polygon_batch)):
        poly = polygon_batch[poly_idx]
        subbox_id = chunk_num

        # Find overlapping grid cells
        overlapping_indices = str_tree.query(poly, predicate='intersects')
        poly_overlap_areas_and_indices = []

        for idx in overlapping_indices:
            overlap_poly = poly.intersection(index_list[idx][2])
            #area = calculate_box_poly_area(overlap_poly)
            #area = calculate_area_in_projection(overlap_poly, eq_transformer)
            area = calculate_spherical_polygon_area(overlap_poly)
            poly_overlap_areas_and_indices.append((area, index_list[idx][0], index_list[idx][1]))

        # Normalize weights
        poly_overlap_areas = [area for area, _, _ in poly_overlap_areas_and_indices]
        total_area = sum(area for area in poly_overlap_areas if not np.isnan(area))
        weights = [(area / total_area) if total_area > 0 and not np.isnan(total_area) else np.nan for area in
                   poly_overlap_areas]

        chunk_data_rows = []

        mask_this_subbox = False
        # Processing each time step
        for time_idx in range(temp_data.shape[0]):
            year, month = time_tuples[time_idx]

            weighted_temp = 0
            valid_weights_sum = 0
            valid_weight_count = 0
            corrected_weights = []

            weights_array = np.array(weights, dtype=np.float32)
            valid_weights_mask = ~np.isnan(weights_array) & (weights_array > 0)

            temps = []
            for (_, lat_idx, lon_idx), weight, is_valid_weight in zip(poly_overlap_areas_and_indices, weights,
                                                                      valid_weights_mask):
                if is_valid_weight:
                    temp = temp_data[time_idx, lat_idx, lon_idx]
                    if not isinstance(temp, masked_constant_type):
                        #if temp >= cutoff_temp:
                        corrected_weights.append(weight)
                        valid_weights_sum += weight
                        valid_weight_count += 1
                        temps.append(temp)
                        #else:
                        #    corrected_weights.append(None)
                        #    temps.append(None)
                    else:
                        corrected_weights.append(None)
                        temps.append(None)
                else:
                    corrected_weights.append(None)
                    temps.append(None)

            if mask_this_subbox:
                break

            if valid_weights_sum > 0:
                normed_weights = [w / valid_weights_sum if w is not None else None for w in corrected_weights]
            else:
                normed_weights = [None] * len(corrected_weights)

            for weight, temp in zip(normed_weights, temps):
                if weight is not None:
                    weighted_temp += (temp * weight)

            if valid_weight_count == 0:
                weighted_temp = np.nan

            year_mo = f"{year}_{month}"
            rounded_bbox = rounded_bbox_batch[poly_idx]
            #lon_min_100, lat_min_100, lon_max_100, lat_max_100 = rounded_bbox.bounds
            lat_min_100, lat_max_100, lon_min_100, lon_max_100 = rounded_bbox
            abs_temp_float32 = np.float32(weighted_temp)
            chunk_data_rows.append({
                'subbox_id': subbox_id,
                'year_month': year_mo,
                'year': year,
                'month': month,
                'lat_min_100': lat_min_100,
                'lat_max_100': lat_max_100,
                'lon_min_100': lon_min_100,
                'lon_max_100': lon_max_100,
                'temperature': abs_temp_float32
            })

        if not mask_this_subbox:
            chunk_df = pd.DataFrame(chunk_data_rows)
            chunk_df.set_index(['subbox_id', 'year_month'], inplace=True)
            chunk_name = f'chunk_{chunk_num:04}.parquet'
            chunk_path = os.path.join(chunk_folder, f'{chunk_name}')
            chunk_df.to_parquet(chunk_path)

            chunk_paths.append(chunk_path)

        chunk_num += 1

        # Increment the progress counter
        with lock:
            progress_counter.value += 1

    return chunk_paths


def process_gridboxes_parallel(polygons, rounded_bboxs, temp_data, str_tree, index_list, time_tuples):
    """
    Process gridboxes in parallel and concatenate results.

    Parameters:
    - polygons: List of shapely.geometry.Polygon objects
    - temp_data: 3D array of temperature data (time, lat, unsigned lon)
    - str_tree: Spatial index (STRtree) for grid cells (signed lons)
    - index_list: List of tuples (lat, unsigned lon, Polygon)
    - time_tuples: List of (year, month) tuples

    Returns:
    - full_df: Concatenated DataFrame from all chunks
    """
    global intermediate_path
    # from step 4 of gistemp (don't include sea ice)
    global sea_surface_cutoff_temp
    # Apply the condition to mask values less than the cutoff temperature
    # mask = temp_data < sea_surface_cutoff_temp
    # Update the existing mask
    # temp_data = np.ma.masked_where(mask, temp_data)
    # free memory
    # del mask

    chunk_folder = intermediate_path
    os.makedirs(chunk_folder, exist_ok=True)

    # speed up determining whether value is masked
    masked_constant_type = np.ma.core.MaskedConstant

    # Initialize progress bar
    num_poly = len(polygons)
    chunk_paths = []

    chunk_completion_filename = 'chunks_complete.txt'
    chunk_completion_file_path = os.path.join(intermediate_path, chunk_completion_filename)
    chunks_complete_fid = open(chunk_completion_file_path, 'w')
    # Process polygons in parallel

    num_workers = os.cpu_count()

    num_polygons = len(polygons)
    batch_size = np.ceil(num_polygons / num_workers).astype(np.int16)  # Calculate batch size based on the number of workers

    # Split polygons into batches
    polygon_batches = [polygons[i:i + batch_size] for i in range(0, num_polygons, batch_size)]
    rounded_bbox_batches = [rounded_bboxs[i:i + batch_size] for i in range(0, num_polygons, batch_size)]
    start_chunk_nums = np.array([i for i in range(0, num_polygons, batch_size)]).astype(np.int16)

    cutoff_temp = sea_surface_cutoff_temp

    with multiprocessing.Manager() as manager:
        # Shared progress counter
        progress_counter = manager.Value('i', 0)
        lock = manager.Lock()

        # Initialize tqdm progress bar
        with tqdm(total=num_polygons, desc="Processing gridboxes") as pbar:
            def update_progress():
                # Update the progress bar based on the shared counter
                pbar.n = progress_counter.value
                pbar.last_print_n = pbar.n
                pbar.update(0)

            # Periodically update the progress bar
            def periodic_update():
                try:
                    while True:
                        update_progress()
                        time.sleep(1)  # Update every second (adjust as needed)
                except Exception as exc:
                    # in case of finished
                    return

            # Start progress update thread
            import threading
            progress_thread = threading.Thread(target=periodic_update, daemon=True)
            progress_thread.start()

            with ProcessPoolExecutor(max_workers=num_workers) as executor:
                futures = {
                    executor.submit(process_batch, batch_idx, polygon_batches[batch_idx], rounded_bbox_batches[batch_idx],
                                    temp_data, str_tree, index_list, time_tuples, chunk_folder, masked_constant_type, cutoff_temp, start_chunk_nums[batch_idx], progress_counter, lock): batch_idx
                    for batch_idx in range(len(polygon_batches))
                }

                # Monitor progress
                while any(future.running() for future in futures):
                    time.sleep(1)  # Adjust sleep duration as needed
                    update_progress()
                # Ensure the progress bar is fully updated
                update_progress()

                for future in futures:
                    try:
                        result = future.result()  # Wait for each batch to complete (1 per thread)
                        for chunk_path in result:
                            chunk_name = os.path.basename(chunk_path)
                            chunks_complete_fid.write(f'{chunk_name}\n')
                        chunks_complete_fid.flush()
                        chunk_paths.extend(result)
                    except Exception as exc:
                        print(traceback.format_exc())
                        print(f'Polygon processing generated an exception: {exc}')
                        chunks_complete_fid.close()
                        print("Exiting.")
                        # Exit everything
                        executor.shutdown(wait=False, cancel_futures=True)
                        sys.exit(1)  # Or raise the exception to exit

    chunks_complete_fid.write("FINISHED CHUNKS")
    chunks_complete_fid.close()

    # Concatenate all chunks
    print("Concatenating chunks...")
    full_df = concatenate_chunks(sorted(chunk_paths))

    print("Finished generating initial df (absolute temps)...")
    return full_df


def read_parquet(chunk_path):
    return pd.read_parquet(chunk_path)


def concatenate_chunks(chunk_paths):
    with ProcessPoolExecutor() as executor:
        df_list = list(executor.map(read_parquet, chunk_paths))
    full_df = pd.concat(df_list)

    # Check if the index is already set (it should be)
    if not {'subbox_id', 'year_month'}.issubset(full_df.index.names):
        print("Setting index....")
        full_df.set_index(['subbox_id', 'year_month'], inplace=True)
    return full_df


# returns only a complete chunk set otherwise []
def get_completed_chunk_set():
    global intermediate_path
    # Define the path to the completion file
    chunk_completion_filename = 'chunks_complete.txt'
    chunk_completion_file_path = os.path.join(intermediate_path, chunk_completion_filename)

    # Check if intermediate_path exists
    if not os.path.exists(intermediate_path):
        print(f"Intermediate path '{intermediate_path}' does not exist.")
        return []

    # Check if the completion file exists
    if not os.path.isfile(chunk_completion_file_path):
        print(f"Completion file '{chunk_completion_filename}' does not exist.")
        return []

    # Read the last line of the completion file
    with open(chunk_completion_file_path, 'r') as f:
        lines = f.readlines()
        if not lines:
            print("Completion file is empty.")
            return []

        last_line = lines[-1].strip()
        if last_line != "FINISHED CHUNKS":
            print(f"Last line in completion file is not the magic string: {last_line}")
            return []

    # Check if all chunk files listed in the completion file exist
    listed_files = {line.strip() for line in lines[:-1]}  # Exclude the last line which is the magic string
    missing_files = [f for f in listed_files if not os.path.isfile(os.path.join(intermediate_path, f))]

    if missing_files:
        return []

    chunk_paths = [os.path.join(intermediate_path, f) for f in sorted(listed_files)]
    return chunk_paths


def cleanup_chunks():
    global intermediate_path
    # Define the pattern for chunk files and the completion file
    chunk_file_pattern = os.path.join(intermediate_path, 'chunk_*.parquet')
    chunk_completion_filename = 'chunks_complete.txt'
    chunk_completion_file_path = os.path.join(intermediate_path, chunk_completion_filename)

    # Delete all chunk files
    for chunk_file in glob.glob(chunk_file_pattern):
        try:
            os.remove(chunk_file)
            #print(f"Deleted file: {chunk_file}")
        except:
            pass

    # Delete the completion file if it exists
    if os.path.isfile(chunk_completion_file_path):
        try:
            os.remove(chunk_completion_file_path)
            #print(f"Deleted file: {chunk_completion_file_path}")
        except:
            pass

    # Remove the intermediate path directory if it exists
    if os.path.isdir(intermediate_path):
        try:
            shutil.rmtree(intermediate_path)
            #print(f"Deleted directory: {intermediate_path}")
        except:
            pass


# have a mean that returns nan if any values are nan
def custom_mean(series):
    if series.isna().any():
        return np.float32(np.nan)
    return series.mean()


def calculate_and_add_monthly_subbox_anomalies(temp_data_df, start_year, end_year):
    global baseline_parquet_path
    print("Getting masked subboxes from original SBBX...")
    masked_subboxes_df, unmasked_subboxes_df = find_masked_subboxes_in_orig_txt()

    if not os.path.exists(baseline_parquet_path):
        print("Calculating baseline temps...")
        # Calculate subbox baselines by month
        subbox_baseline_by_month_df = temp_data_df[
            (temp_data_df['year'] >= start_year) &
            (temp_data_df['year'] <= end_year)
            ].groupby(['subbox_id', 'month'])['temperature'].agg(custom_mean).reset_index()
        subbox_baseline_by_month_df.set_index(['subbox_id', 'month'], inplace=True)

        print("Writing baseline parquet...")
        subbox_baseline_by_month_df.to_parquet(baseline_parquet_path)
    else:
        print("Loading baseline temps from parquet (already existing)...")
        subbox_baseline_by_month_df = pd.read_parquet(baseline_parquet_path)

    # Add anomaly column to original df
    anomaly_col = f'anomaly_{start_year}_{end_year}'
    print("Calculating Anomalies ...")

    print("Adding anomaly column...")
    temp_data_df[anomaly_col] = np.float32(np.nan)

    print("Resetting index for merge...")
    # Reset index to perform merge
    temp_data_df_reset = temp_data_df.reset_index()

    print("Merging baseline temps...")
    # Merge temp_data_df with subbox_baseline_by_month_df on subbox_id and month
    merged_df = pd.merge(temp_data_df_reset, subbox_baseline_by_month_df.reset_index(),
                         on=['subbox_id', 'month'],
                         suffixes=('', '_baseline'))

    print("Calculating the anomaly temps in merged_df...")
    # Calculate the anomaly
    merged_df[anomaly_col] = merged_df['temperature'] - merged_df['temperature_baseline']

    print("Setting merged_df index ")
    merged_df.set_index(['subbox_id', 'year_month'], inplace=True)

    print("Updating anomaly temps column in temp_data_df...")
    # Update the original DataFrame with the anomaly column
    temp_data_df[anomaly_col] = merged_df[anomaly_col]

    print("Masking months per subbox per original SBBX (row)")
    # Next we need to mask the months for the subboxes that need masking per GISS (~2000 or so data points)
    # Create a unique key for the subboxes for efficient matching so we can merge the nans on matching subboxes
    print("Creating masked_subboxes_df keys")
    masked_subboxes_df['key'] = (
            masked_subboxes_df['lat_min_100'].astype(str) + '_' +
            masked_subboxes_df['lat_max_100'].astype(str) + '_' +
            masked_subboxes_df['lon_min_100'].astype(str) + '_' +
            masked_subboxes_df['lon_max_100'].astype(str) + '_' +
            masked_subboxes_df['month'].astype(str)
    )
    print("Creating temp_data_df keys")
    temp_data_df['key'] = (
            temp_data_df['lat_min_100'].astype(str) + '_' +
            temp_data_df['lat_max_100'].astype(str) + '_' +
            temp_data_df['lon_min_100'].astype(str) + '_' +
            temp_data_df['lon_max_100'].astype(str) + '_' +
            temp_data_df['month'].astype(str)
    )
    # Merge DataFrames on the 'key' column
    print("Merging dfs...")
    merged_df = pd.merge(temp_data_df, masked_subboxes_df[['key']], on='key', how='left', indicator=True)
    # Update the 'anomaly_col' where there is a match
    #temp_data_df.loc[merged_df['_merge'] == 'both', anomaly_col] = np.float32(np.nan)

    # Update the 'anomaly_col' where there is a match
    print("Updating on merged both...")
    temp_data_df.loc[
        temp_data_df['key'].isin(merged_df.loc[merged_df['_merge'] == 'both', 'key']), anomaly_col] = np.float32(np.nan)

    print("Masking anomaly_col where no match is found in unmasked_subboxes_df")

    # Create unique keys for both temp_data_df and unmasked_subboxes_df
    print("Creating unmasked_subboxes_df keys")
    unmasked_subboxes_df['key'] = (
            unmasked_subboxes_df['lat_min_100'].astype(str) + '_' +
            unmasked_subboxes_df['lat_max_100'].astype(str) + '_' +
            unmasked_subboxes_df['lon_min_100'].astype(str) + '_' +
            unmasked_subboxes_df['lon_max_100'].astype(str) + '_' +
            unmasked_subboxes_df['month'].astype(str)
    )

    # Mask values in anomaly_col where there is no corresponding match in unmasked_subboxes_df
    print("Masking unmatched keys in temp_data_df...")
    temp_data_df.loc[
        ~temp_data_df['key'].isin(unmasked_subboxes_df['key']), anomaly_col] = np.float32(np.nan)


    print("Dropping key column...")
    # Drop the temporary 'key' column
    temp_data_df.drop(columns='key', inplace=True)

    print("Done calculating anomalies...")

    return temp_data_df


def generate_subbox_parquets():
    global baseline_start_year
    global baseline_end_year
    global ersst_monthly_means
    global subbox_temps_parquet_path
    global subbox_temps_anomalies_parquet_path
    global reference_subbox_year

    # right now the function is not to overwrite subbox_temps.parquet (manually deletion is required)
    if not os.path.exists(subbox_temps_parquet_path):
        completed_chunk_set = get_completed_chunk_set()
        if not completed_chunk_set:
            # new calculation, or, restart calculation since it failed part way through
            cleanup_chunks()

            temp_data, lats, lons, time_tuples = read_ersst_monthly_means_absolute(ersst_monthly_means)

            # Assuming `filename` is your file and `year` is the target year
            polygons, rounded_bboxs = generate_grid8k_polygons()

            # Create spatial index and process polygons
            str_tree, index_list = create_spatial_index(lats, lons)
            print("Process polygons (regrid data into 8k subboxes)...")
            subbox_temps_df = process_gridboxes_parallel(polygons, rounded_bboxs, temp_data, str_tree, index_list, time_tuples)
        else:
            print("Concatenating chunks (reusing last intermediate chunked data from crash?)")
            subbox_temps_df = concatenate_chunks(completed_chunk_set)

        # store the subbox df
        print("Save subbox df (absolute only)...")
        subbox_temps_df.to_parquet(subbox_temps_parquet_path, compression='snappy', index=True)
        print("Done writing absolute temperatures parquet.")
    else:
        print(f"Loading existing parquet for absolute temps: {subbox_temps_parquet_path}")
        subbox_temps_df = pd.read_parquet(subbox_temps_parquet_path)
        print("Done loading absolute temps parquet.")

    print("Calculate and add column for monthly subbox anomalies")
    # calculate the anomalies for the subboxes, modifying the df
    calculate_and_add_monthly_subbox_anomalies(subbox_temps_df, baseline_start_year, baseline_end_year)

    print("Done calculating subbox anomalies parquet")

    print("Saving subbox df for temps and anomalies")
    subbox_temps_df.to_parquet(subbox_temps_anomalies_parquet_path, compression='snappy', index=True)
    print("Done writing temps, anomalies parquet")

    print("Cleaning up intermediate chunks")
    cleanup_chunks()
    print("Finished with all parquets.\n")

    print("subbox_temps_df:")
    print(subbox_temps_df)

    return subbox_temps_df


def get_parquet():
    global do_generate_subbox_parquets
    global subbox_temps_anomalies_parquet_path
    print(f"Getting parquets (this will take a while)...")
    if do_generate_subbox_parquets:
        subbox_temps_df = generate_subbox_parquets()
    else:
        try:
            # load the parquets
            print(f"Loading parquet from {subbox_temps_anomalies_parquet_path}")
            subbox_temps_df = pd.read_parquet(subbox_temps_anomalies_parquet_path)
            print("Finished loading parquet.")
        except:
            raise ValueError("Error loading parquet for subbox temp anomalies.\nSet do_generate_subbox_parquets = True to generate, or specify correct path.\nExiting")
            exit(1)
    print("Finished getting parquet.")

    return subbox_temps_df


def df_to_sbbx(df, output_file):
    global baseline_start_year
    global baseline_end_year
    anomaly_col = f'anomaly_{baseline_start_year}_{baseline_end_year}'
    # Ensure the year column is a float
    df['year'] = df['year'].astype(float)

    # Filter the DataFrame based on the start_clip_year
    df_filtered = df[df['year'] >= start_clip_year]

    # Extract unique year-month periods
    year_months = set(df_filtered.index.get_level_values('year_month'))
    sorted_year_months = sorted(year_months, key=lambda x: (int(x.split('_')[0]), int(x.split('_')[1])))
    first_year_month = sorted_year_months[0]
    last_year_month = sorted_year_months[-1]
    # Convert year_month to integers
    first_year, first_month = map(int, first_year_month.split('_'))
    last_year, last_month = map(int, last_year_month.split('_'))

    # Adjust title and record calculations
    title = f"Monthly Sea Surface Temperature anom (C) ERSSTv5 {first_month:02d}/{first_year} - {last_month:02d}/{last_year}".ljust(
        80)

    global missing_flag
    with open(output_file, 'wb') as f:
        s = File(f, bos='>')

        # Header fields: (missing_flag: 9999)
        num_months = (last_year - first_year) * 12 + (last_month - first_month + 1)
        header_record = [num_months, 0, 6, num_months, 0, first_year, missing_flag, 0, title.encode('ascii')]
        header = struct.pack('>8i80s', *header_record)
        s.writeline(header)

        subbox_ids = df_filtered.index.get_level_values('subbox_id').unique()
        for subbox_id in tqdm(subbox_ids, desc="Writing subboxes to SBBX"):
            subbox_df = df_filtered.loc[subbox_id]
            time_series = np.full(num_months, 9999.0, dtype=np.float32)
            for _, row in subbox_df.iterrows():
                # must force int
                year = int(row['year'])
                month = int(row['month'])
                anomaly = row[anomaly_col]
                index = (year - first_year) * 12 + month - first_month
                if 0 <= index < num_months and not np.isnan(anomaly):
                    time_series[index] = anomaly

            # Record format: 7 integers + num_months floats
            zero_int = int(0)
            zero_float = float(0)

            # Extract and scale latitude and longitude values
            # Must round to get identical subboxes
            lat_lon_values = [int(row['lat_min_100']), int(row['lat_max_100']), int(row['lon_min_100']), int(row['lon_max_100'])]

            # Prepare the record list
            record = [num_months] + lat_lon_values + [zero_int, zero_int, zero_float] + time_series.tolist()

            # Pack the record
            record_format = f'>7if{num_months}f'
            record_bytes = struct.pack(record_format, *record)

            # Write the record
            s.writeline(record_bytes)


if __name__ == '__main__':
    # Get the subbox df with the anomalies (generate or from disk)
    temps_and_anomalies_subbox_temps_df = get_parquet()

    # Convert to SBBX format (clipping the data for years when there is not useful data)
    df_to_sbbx(temps_and_anomalies_subbox_temps_df, output_sbbx)

    print(f"Finished!")
    end_time = time.time()
    print(f"Total running time: {end_time - start_time} seconds")
