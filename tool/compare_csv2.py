import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np
import os

# Check if two file arguments are provided
if len(sys.argv) < 3:
    print("Usage: python script.py file1.csv file2.csv subtitlestring")
    sys.exit(1)
subtitle = ""
if len(sys.argv) == 4:
    subtitle = "\n"+sys.argv[3]

# Load the files, skipping the first two rows and the first column
df1 = pd.read_csv(sys.argv[1], skiprows=2, usecols=range(1, 13))
df2 = pd.read_csv(sys.argv[2], skiprows=2, usecols=range(1, 13))

# Select only numeric columns
df1 = df1.apply(pd.to_numeric, errors='coerce')
df2 = df2.apply(pd.to_numeric, errors='coerce')

diff_values = []
for (index1, row1), (index2, row2) in zip(df1.iterrows(), df2.iterrows()):
    for col in range(len(row1)):
        if index1 == index2:  # Ensure corresponding rows
            # make sure 2 decimal precision (as gistemp is only 2 decimals)
            diff = round(abs(row1[col] - row2[col]), 2)
            if not np.isnan(diff):
                if diff > 0.01:
                    print(diff)
                diff_values.append(diff)

diff_values = pd.Series(diff_values)
total_diff = diff_values[diff_values > 0].count()


# Create a histogram with a bin size of 0.01
plt.hist(diff_values[diff_values > 0], bins=np.arange(0, diff_values.max() + 0.01, 0.001), edgecolor='black')
plt.xlabel('Absolute Difference')
plt.ylabel('Count')
plt.title(f'Histogram of Absolute Differences (Total Differences: {total_diff}){subtitle}')

# Save the figure to a PNG file
plt.savefig(f'compare_csv_histogram.png', bbox_inches='tight')
