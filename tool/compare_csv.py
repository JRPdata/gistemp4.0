# generates a comparison between csv files (for global monthly means csv)
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

# Calculate the absolute differences between the two files
diff_values = (df1 - df2).abs().stack()

total_diff = diff_values[diff_values > 0].count()

# Create a histogram with a bin size of 0.01
plt.hist(diff_values[diff_values > 0], bins=np.arange(0, diff_values.max() + 0.01, 0.01), edgecolor='black')
plt.xlabel('Absolute Difference')
plt.ylabel('Count')
plt.title(f'Histogram of Absolute Differences (Total Differences: {total_diff}){subtitle}')


# Save the figure to a PNG file
plt.savefig(f'compare_csv_histogram.png', bbox_inches='tight')

print(f'Figure saved to compare_csv_histogram.png')
