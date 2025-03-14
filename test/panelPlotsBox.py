import matplotlib.pyplot as plt
import numpy as np
import glob

# Lists to store extracted data
ngrid_values = []
r_C_values = []
E_values = []
epair_values = []
gap_values = []

# Read multiple files
file_list = glob.glob("collect_cat_box_Paper/collect_cat_SVbas_146Gd_2.0d0_hfb_*.dat")  # Adjust path to your data files

for filename in file_list:
    with open(filename, "r") as file:
        lines = file.readlines()
        for i in range(0, len(lines), 2):  # Each data block spans two lines
            # Extract ngrid from the first line
            ngrid = int(lines[i].split(",")[0].split("=")[1])  

            # Extract r_C, E, epair, gap from the second line
            values = lines[i + 1].split()[3:9]  # Skip Z and N
            r_C, E, epair,epairn, gap, gapn = map(float, values)

            # Store values
            ngrid_values.append(ngrid)
            r_C_values.append(r_C)
            E_values.append(E)
            epair_values.append(epair)
            gap_values.append(gap)

# Convert lists to numpy arrays and sort by ngrid
ngrid_values, r_C_values, E_values, epair_values, gap_values = zip(*sorted(zip(ngrid_values, r_C_values, E_values, epair_values, gap_values)))

# Create a 4-panel plot
fig, axes = plt.subplots(4, 1, figsize=(6, 10), sharex=True)
fig.subplots_adjust(hspace=0.3)

# Define data and labels
data_series = [r_C_values, E_values, epair_values, gap_values]
y_labels = [r'$r_C$', r'$E$', r'$e_{pair}$', r'$gap$']

for ax, y, label in zip(axes, data_series, y_labels):
    ax.plot(np.array(ngrid_values)*0.3, y, marker='o', linestyle='-', color='b')
    ax.set_ylabel(label, fontsize=12)
    ax.grid(True)

axes[-1].set_xlabel(r'$ngrid$', fontsize=12)

plt.show()
