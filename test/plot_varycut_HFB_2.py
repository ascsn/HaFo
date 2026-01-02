import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from pathlib import Path
import re # Import regex for parsing

# Read the content of the gnuplot file
gnuplot_file_path = "pairing/varycut_HFB_2.gnu.txt"
file_content = Path(gnuplot_file_path).read_text()

# --- 1. Parse Data from the gnuplot file ---
# Find the start of the data section based on the gnuplot 'plot' command structure
# The actual data is usually after the plot commands, or interspersed.
# In this specific file, the data lines appear at the very end.
data_lines = []
in_data_section = False
for line in file_content.splitlines():
    if line.strip().startswith('p  "varycut_HFB_2.dat"'):
        # This line usually comes before data in a combined gnuplot file
        # In this specific file, the data is after a series of comments/settings
        pass
    if re.match(r'^-?\d+\s+\d+\.\d+', line.strip()): # Lines that start with an optional minus, digits, then space, then digits.decimal
        data_lines.append(line.strip())

if not data_lines:
    # If no data lines found this way, it might be in an unexpected format, or truly missing
    print(f"Warning: No data lines found in {gnuplot_file_path}. Plot will be empty.")
    data_np = np.empty((0, 3)) # Empty array
else:
    parsed_data = []
    for line in data_lines:
        parts = line.split()
        # Gnuplot 'u 1:($2+...)' means use the 1st data column (index 1 in parts), and 2nd data column (index 2 in parts)
        # Gnuplot 'u 1:($3-...)' means use the 1st data column (index 1 in parts), and 3rd data column (index 3 in parts)
        # However, observing the data:
        # -80   4.645548  -1016.36      -0.000 -14.016   0.000   1.590   0.000   8.816
        # Col 1: 4.645548 (x-axis, cutoff energy) -> parts[1]
        # Col 2: -1016.36 (total energy) -> parts[2]
        # Col 9: 8.816 (charge radius) -> parts[8] (assuming 0-indexed after split)
        try:
            x_val = float(parts[1])
            energy_val = float(parts[2])
            charge_radius_val = float(parts[8]) # This assumes the 9th column (index 8) is charge radius
            parsed_data.append([x_val, energy_val, charge_radius_val])
        except (ValueError, IndexError) as e:
            print(f"Skipping malformed data line: '{line}' due to error: {e}")
            continue
    data_np = np.array(parsed_data)

# Define offsets from gnuplot commands
energy_offset = 1016.3533
charge_radius_offset = 8.8180

# --- 2. Matplotlib Setup ---
mpl.use("Agg") # Set backend for headless environment
mpl.rcParams["text.usetex"] = False # Disable LaTeX as it requires external setup
mpl.rcParams["axes.linewidth"] = 1.5
mpl.rcParams["axes.formatter.useoffset"] = False
mpl.rcParams["font.family"] = "DejaVu Sans" # Default font
mpl.rcParams["font.size"] = 18 # Base font size, will adjust specific elements later

# Font settings derived from gnuplot "TimesNewRoman,30" (adjusting for non-latex)
fontDict = {"fontsize": 28, "fontweight": "bold"} # For labels
fontDictTitle = {"fontsize": 30, "fontweight": "bold"} # For titles

# Colors and markers from gnuplot 'set style line'
# Note: gnuplot pt (point type) mapping to matplotlib markers
# pt 7: circle (o), pt 5: filled square (s), pt 4: filled circle (.), pt 6: triangle up (^), pt 9: triangle down (v)
line_styles = {
    1: {'color': '#FF0000', 'marker': 'o', 'linestyle': '-', 'linewidth': 4, 'markersize': 10},
    2: {'color': '#00AA00', 'marker': 's', 'linestyle': '-', 'linewidth': 4, 'markersize': 10},
    3: {'color': '#0000FF', 'marker': '.', 'linestyle': '-', 'linewidth': 4, 'markersize': 10}, # Filled circle
    4: {'color': '#AA00AA', 'marker': '^', 'linestyle': '-', 'linewidth': 4, 'markersize': 10},
    5: {'color': '#00CCAA', 'marker': 'v', 'linestyle': '-', 'linewidth': 4, 'markersize': 10},
}

# Gnuplot plot command labels
g_values = ["g=0.0", "g=0.01", "g=0.02", "g=0.03", "g=0.04"]

# --- 3. Figure and Subplots ---
# Gnuplot 'set multiplot' with 'set size' and 'set origin' suggests a 2x1 layout
# Top plot (charge radius): size 0.8,0.3, origin 0.2,0.6
# Bottom plot (energy): size 0.8,0.6, origin 0.2,0.2
# This implies a ratio of heights 0.3 : 0.6 = 1:2. Shared x-axis.
fig, axs = plt.subplots(2, 1, figsize=(10, 15), sharex=True, gridspec_kw={'height_ratios': [0.3, 0.6]})
plt.subplots_adjust(hspace=0) # gnuplot 'set tmargin 0', 'set bmargin 0' etc.

# --- 4. Bottom Subplot (Total Energy) ---
ax_bottom = axs[1]
ax_bottom.set_xlabel("cutoff energy [MeV]", fontdict=fontDict, labelpad=5)
ax_bottom.set_ylabel("total energy [MeV]", fontdict=fontDict, labelpad=10)
ax_bottom.set_xlim(14, 101)
ax_bottom.set_ylim(0, 3.3)

# Plotting
if data_np.shape[0] > 0:
    x_energy = data_np[:, 0]
    y_energy = data_np[:, 1] + energy_offset

    # Only plot the first series (i 0) as only its data is available in the provided file
    style = line_styles[1] # Use style 1 for g=0.0
    ax_bottom.plot(
        x_energy, y_energy,
        color=style['color'],
        marker=style['marker'],
        linestyle=style['linestyle'],
        linewidth=style['linewidth'],
        markersize=style['markersize'],
        label=g_values[0] # Label for g=0.0
    )
else:
    print("No energy data to plot.")

ax_bottom.legend(loc='upper right', prop={'size': 16}, bbox_to_anchor=(1, 0.9)) # Adjusted bbox_to_anchor slightly
ax_bottom.tick_params(
    axis="both",
    which="major",
    labelsize=22, # Keep labelsize as set in previous interaction, or adjust
    direction="in",
    right=True,
    top=True,
    size=8,
)
ax_bottom.tick_params(
    axis="both",
    which="minor",
    labelsize=22, # Keep labelsize as set in previous interaction, or adjust
    direction="in",
    right=True,
    top=True,
    size=4,
)


# --- 5. Top Subplot (Charge Radius) ---
ax_top = axs[0]
plt.setp(ax_top.get_xticklabels(), visible=False) # Gnuplot 'unset xlabel'
ax_top.set_ylabel("charge radius [fm]", fontdict=fontDict, labelpad=10)
ax_top.set_xlim(14, 101)
ax_top.set_ylim(-0.012, 0.007)

# Plotting
if data_np.shape[0] > 0:
    x_radius = data_np[:, 0]
    y_radius = data_np[:, 2] - charge_radius_offset # Using 9th column (index 8) from raw data

    # Only plot the first series (i 0) as only its data is available
    style = line_styles[1] # Use style 1 for g=0.0
    ax_top.plot(
        x_radius, y_radius,
        color=style['color'],
        marker=style['marker'],
        linestyle=style['linestyle'],
        linewidth=style['linewidth'],
        markersize=style['markersize'],
        label=g_values[0] # Use label for g=0.0, gnuplot uses t "" so we add it only for reference
    )
else:
    print("No charge radius data to plot.")

ax_top.legend(loc='upper right', prop={'size': 16}, bbox_to_anchor=(1, 0.9)) # Adjusted bbox_to_anchor slightly
ax_top.tick_params(
    axis="both",
    which="major",
    labelsize=22, # Keep labelsize as set in previous interaction, or adjust
    direction="in",
    right=True,
    top=True,
    size=8,
)
ax_top.tick_params(
    axis="both",
    which="minor",
    labelsize=22, # Keep labelsize as set in previous interaction, or adjust
    direction="in",
    right=True,
    top=True,
    size=4,
)

# --- Final Touches ---
fig.align_ylabels(axs)
plt.tight_layout() # Handles margins/size equivalent to gnuplot settings

# Save output
plt.savefig("varycut_HFB_2.pdf", dpi=300, bbox_inches="tight")