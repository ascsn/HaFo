import matplotlib.pyplot as plt
import numpy as np
import glob
import matplotlib as mpl
import matplotlib.ticker as ticker
import scienceplots

# Lists to store extracted data
ngrid_values = []
r_C_values = []
E_values = []
epair_values = []
gap_values = []

mpl.use("TkAgg")
plt.style.use("science")
mpl.rcParams["text.usetex"] = True
mpl.rcParams["axes.linewidth"] = 1.5
mpl.rcParams["axes.formatter.useoffset"] = False

# Initialize lists for cut25 and cut100
box_sizes_25, r_C_25, E_25, epair_25, gap_25 = [], [], [], [], []
box_sizes_100, r_C_100, E_100, epair_100, gap_100 = [], [], [], [], []

fontDict = {"family": "helvetica", "fontsize": 22, "fontweight": "bold"}
fontDictTitle = {"family": "helvetica", "fontsize": 22, "fontweight": "bold"}


# Read multiple files
file_list = glob.glob("collect_cat_box_Paper/collect_cat_SVbas_120Sn_0.0d0_hfb_*.dat")  # Adjust path to your data files

print(file_list)

for filename in file_list:
    with open(filename, "r") as file:
        lines = file.readlines()
        
        # Extract cutoff type:
        if "Cut25" in filename:
            cut_type = "cut25"
        elif "Cut100" in filename:
            cut_type = "cut100"
        else:
            continue

        for i in range(0, len(lines), 2):  # Each data block spans two lines
            # Extract ngrid from the first line
            ngrid = int(lines[i].split(",")[0].split("=")[1])  

            # Extract r_C, E, epair, gap from the second line
            values = lines[i + 1].split()[3:9]  # Skip Z and o
            r_C, E, epair,epairn, gap, gapn = map(float, values)
            # print(values)
            # Store values
            if cut_type == "cut25":   
                box_sizes_25.append(ngrid*0.3)
                r_C_25.append(r_C)
                E_25.append(E)
                epair_25.append(epair)
                gap_25.append(gap)
            elif cut_type == "cut100":
                box_sizes_100.append(ngrid*0.3)
                r_C_100.append(r_C)
                E_100.append(E)
                epair_100.append(epair)
                gap_100.append(gap)


# Convert lists to numpy arrays and sort by ngrid
box_sizes_25, r_C_25, E_25, epair_25, gap_25 = zip(*sorted(zip(box_sizes_25, r_C_25, E_25, epair_25, gap_25)))
box_sizes_100, r_C_100, E_100, epair_100, gap_100 = zip(*sorted(zip(box_sizes_100, r_C_100, E_100, epair_100, gap_100)))
                                                   
# Convert to NumPy arrays and sort by box size
def sort_data(box_sizes, values):
    box_sizes = np.array(box_sizes)
    sorted_indices = np.argsort(box_sizes)
    return box_sizes[sorted_indices], np.array(values)[sorted_indices]

box_sizes_25, r_C_25 = sort_data(box_sizes_25, r_C_25)
box_sizes_25, E_25 = sort_data(box_sizes_25, E_25)
box_sizes_25, epair_25 = sort_data(box_sizes_25, epair_25)
box_sizes_25, gap_25 = sort_data(box_sizes_25, gap_25)

box_sizes_100, r_C_100 = sort_data(box_sizes_100, r_C_100)
box_sizes_100, E_100 = sort_data(box_sizes_100, E_100)
box_sizes_100, epair_100 = sort_data(box_sizes_100, epair_100)
box_sizes_100, gap_100 = sort_data(box_sizes_100, gap_100)

# Define consistent colors
cut25_color = "blue"
cut100_color = "red"

# Define label padding
labelpadding = [None, None, None, None]
ylabelList = ["Proton rms radius [fm]", "Total energy [MeV]", "Pairing energy [MeV]", "Proton gap [MeV]"]
xvals = [box_sizes_25, box_sizes_100]
yvals_25 = [r_C_25 - r_C_25[-1], E_25-E_25[-1], epair_25-epair_25[-1], gap_25 - gap_25[-1]]
yvals_100 = [r_C_100 - r_C_100[-1], E_100-E_100[-1], epair_100-epair_100[-1], gap_100 - gap_100[-1]]
# yvals_25 = [r_C_25, E_25, epair_25, gap_25]
# yvals_100 = [r_C_100, E_100, epair_100, gap_100]

label_coords = [ (-0.,.5), (-0.5,.5), (-0.5,.5), (-0.5,.5)]
# yranges = [[-0.00019,0.00006], [-0.12,0.42], [-1.49,0.99], [-0.12, 0.14]]
# Define custom formatter for first panel
def scientific_formatter(x, pos):
    return rf"{x:.1f}$\times 10^{-4}$"

# Create a 4-panel plot without extra space
fig, axs = plt.subplots(4, 1, figsize=(6, 12), sharex=True, gridspec_kw={'hspace': 0})

for i in range(4):
    if i ==3:
        labels = ["cut25", "cut100"]
    else:
        labels = [None, None]
    axs[i].plot(box_sizes_25, yvals_25[i], marker="o", linestyle="-",label=labels[0], color=cut25_color)
    axs[i].plot(box_sizes_100, yvals_100[i], marker="o", linestyle="--",label=labels[1], color=cut100_color)
    axs[i].tick_params(
        axis="both",
        which="major",
        labelsize=16,
        labelright=False,
        direction="in",
        right=True,
        top=True,
        size=8,
    )
    axs[i].tick_params(
        axis="both",
        which="minor",
        labelsize=16,
        labelright=False,
        direction="in",
        right=True,
        top=True,
        size=4,
    )
    # axs[i].set_ylim(yranges[i])
    axs[i].set_ylabel(ylabelList[i], fontdict=fontDict)  
    axs[i].label_coords = label_coords[i]

# Set y-axis tick labels to have 2 decimal places
for ax in axs[1:]:
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))

axs[0].set_title(r"$^{120}$Sn, SV-bas, spherical 1D, vary box size, width=0.0 fm", fontdict= fontDict)
axs[3].set_xlabel("Box size [fm]", labelpad = labelpadding[3], fontdict=fontDictTitle)

# # Add legend only in the bottom panel
axs[3].legend()


# Improve layout
plt.tight_layout()

# Save the plot
plt.savefig("SVbas_146Gd_hfb_props_cut25_100_w00.pdf")

# plt.show()
