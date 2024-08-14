# hacky tool to compare to SBBX txt output files (average RMSE)
import pandas as pd
import numpy as np

# Load reference file
ref_file = pd.read_csv('tmp/input/SBBX.ERSSTv5.orig.txt', skiprows=2, sep=r'\s+')

# Load other files
files = ['tmp/input/SBBX.ERSSTv5.cutoff_project_area.txt', 
         'tmp/input/SBBX.ERSSTv5.no_cutoff_geo_area.txt', 
         'tmp/input/SBBX.ERSSTv5.no_cutoff_project_area.txt']

for file in files:
    # Load file
    df = pd.read_csv(file, skiprows=2, sep=r'\s+')
    
    # Create unique key
    df['key'] = df.iloc[:, :5].apply(lambda x: ' '.join(map(str, x)), axis=1)
    ref_file['key'] = ref_file.iloc[:, :5].apply(lambda x: ' '.join(map(str, x)), axis=1)
    
    # Merge files on unique key
    merged = pd.merge(df, ref_file, on='key', suffixes=('_file', '_ref'))


    # Compute RMSE for each area and month
    rmse_values = []
    for _, group in merged.groupby(merged.columns[:4].tolist()):
        for col in range(6, 17):
            # Filter out rows with NaN or 9999.00
            group_filtered = group[(group.iloc[:, col] != 9999.00) & (group.iloc[:, col] != np.nan) & 
                                (group.iloc[:, col + 17] != 9999.00) & (group.iloc[:, col + 17] != np.nan)]

            # Compute RMSE
            values = group_filtered.iloc[:, col] - group_filtered.iloc[:, col + 17]
            rmse = np.sqrt((values ** 2).mean())
            rmse_values.append(rmse)

    # Convert rmse_values to a pandas Series
    rmse_series = pd.Series(rmse_values)

     # Drop NaN values
    rmse_series = rmse_series.dropna()

    # Average RMSE for the entire file
    avg_rmse = rmse_series.mean()
    print(f'Average RMSE for {file}: {avg_rmse}')
