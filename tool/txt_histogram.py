# geenrates a histogram of anomaly tempeartures for ERSST SBBX txt output
# specify optional fields year, month for more specificity
import numpy as np
import matplotlib.pyplot as plt
import sys

# Define a function to filter data based on year and month
def filter_data(data, year=None, month=None):
    if year is not None:
        data = data[data[:, 4] == year]
    if month is not None:
        data = data[:, 5+month-1]  # Assuming 9999.00 is the fill value
    else:
        data = data[:, 5:16]
    data = data[data != 9999.00]
    return data


# Check if input file path is provided
if len(sys.argv) < 2:
    print("Usage: ersst_txt_hist.py --year=[year] --month_num=[1-12] <path_to_ersst5_txt_file.txt>")
    sys.exit(1)

input_file_path = sys.argv[1]


# Parse command-line arguments
year = None
month = None
for arg in sys.argv[1:]:
    if arg.startswith('--year='):
        year = int(arg.split('=')[1])
    elif arg.startswith('--month_num='):
        month = int(arg.split('=')[1])
    else:
        input_file_path = arg

# Load data from file
data = np.loadtxt(input_file_path, skiprows=2)

# Filter data based on year and month
data = filter_data(data, year, month)

# Calculate statistics and create histogram
mean = np.mean(data)
median = np.median(data)
p67 = np.percentile(data, 67)
p95 = np.percentile(data, 95)

plt.hist(data, bins=np.arange(min(data), max(data) + 0.5, 0.5), edgecolor='black')

# Add annotations
plt.annotate(f'Mean: {mean:.2f}', xy=(0.05, 0.9), xycoords='axes fraction')
plt.annotate(f'Median: {median:.2f}', xy=(0.05, 0.85), xycoords='axes fraction')
plt.annotate(f'67% interval: [{p67:.2f}, {2*p67-p95:.2f}]', xy=(0.05, 0.8), xycoords='axes fraction')
plt.annotate(f'95% interval: [{p95:.2f}, {2*p95-p67:.2f}]', xy=(0.05, 0.75), xycoords='axes fraction')

# Set title and labels
month_str = year_str = comma_str = ''
if month:
    month_str = f'Month num: {month}'
    print(month_str)
if year:
    year_str = f'Year: {year}'
if year and month:
    comma_str = ', '
title_str = f'ERSST5 SBBX Temperatures {year_str}{comma_str}{month_str}'
plt.title(title_str)
plt.xlabel('SBBX Temperatures (C)')
plt.ylabel('# Months')

# Show plot
plt.show()
