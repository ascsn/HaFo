import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d, PchipInterpolator # Import for interpolation
import scienceplots

# Data extracted from the gnuplot script
data_blocks = {
    "SVmin_HFB_cut30": np.array([
        [0.0, 55.92], [0.3, 54.943], [0.5, 54.214], [0.6, 54.02], [0.7, 53.91],
        [0.8, 54.107], [0.9, 54.441], [1.0, 55.084], [1.2, 56.89]
    ]),
    "SVmin_BCS_cut30": np.array([
        [0.0, 55.905], [0.2, 55.011], [0.3, 54.892], [0.5, 54.657], [0.7, 54.801],
        [0.8, 55.155], [1.0, 56.305], [1.2, 58.000], [1.5, 63.344]
    ]),
    "FyDrHFB_cut30": np.array([
        [0.0, 78.391], [0.5, 69.62], [0.7, 64.161], [0.8, 58.603], [0.9, 55.434],
        [1.0, 57.505], [1.1, 60.89], [1.2, 68.86]
    ]),
    "FyDrHFB_cut15": np.array([
        [0.0, 65.633], [0.4, 60.87], [0.5, 57.93], [0.6, 55.291], [0.7, 52.002],
        [0.8, 50.46], [0.9, 53.797], [1.0, 56.67], [1.1, 60.166], [1.2, 67.292]
    ]),
}

# Matplotlib setup (based on your previous preferences and gnuplot commands)
mpl.use("Agg") # Changed backend to Agg for headless environment
plt.style.use("science") # Ensure 'scienceplots' is installed for this to work
mpl.rcParams["text.usetex"] = True
mpl.rcParams["axes.linewidth"] = 1.5
mpl.rcParams["axes.formatter.useoffset"] = False

# Font settings for labels and titles (from previous turns, using numeric fontweight for max boldness)
fontDict = {"family": "serif", "fontsize": 28, "fontweight": 900}
fontDictTitle = {"family": "serif", "fontsize": 28, "fontweight": 900}

# Figure and Subplots (2 rows, 1 column)
fig, axs = plt.subplots(2, 1, figsize=(8, 14), constrained_layout=True)

# --- Top Subplot (FyDr plot, now axs[0]) ---
ax_top = axs[0]
ax_top.set_xlabel("", fontdict=fontDict) # Keep xlabel empty for top plot
ax_top.set_xlim(-0.02, 1.22)
ax_top.set_ylim(49, 79.5)
ax_top.set_yticks(np.arange(50, 80, 5))
ax_top.set_xticks(np.arange(0.0, 1.3, 0.2)) # Set x-ticks at 0.2 intervals
ax_top.set_ylabel(r"quality measure $\chi^2$", fontdict=fontDict, labelpad=10)

# Plotting Data Block 2 (FyDrHFB_cut30) with interpolation - Red line
x_raw_2 = data_blocks["FyDrHFB_cut30"][:, 0]
y_raw_2 = data_blocks["FyDrHFB_cut30"][:, 1]
f_interp_2 = interp1d(x_raw_2, y_raw_2, kind='cubic')
x_smooth_2 = np.linspace(x_raw_2.min(), x_raw_2.max(), 500)
y_smooth_2 = f_interp_2(x_smooth_2)
ax_top.plot(
    x_smooth_2, y_smooth_2,
    label=r"Fy($\Delta r$,HFB) like, $\epsilon_{\text{cut}}$=30MeV",
    color='red', # Changed to red
    linestyle='-',
    linewidth=3
)

# Plotting Data Block 3 (FyDrHFB_cut15) with PchipInterpolator - Green line
x_raw_3 = data_blocks["FyDrHFB_cut15"][:, 0]
y_raw_3 = data_blocks["FyDrHFB_cut15"][:, 1]
# Changed to PchipInterpolator for monotonic interpolation
f_interp_3 = PchipInterpolator(x_raw_3, y_raw_3)
x_smooth_3 = np.linspace(x_raw_3.min(), x_raw_3.max(), 500)
y_smooth_3 = f_interp_3(x_smooth_3)
ax_top.plot(
    x_smooth_3, y_smooth_3,
    label=r"Fy($\Delta r$,HFB) like, $\epsilon_{\text{cut}}$=15MeV",
    color='green', # Changed to green
    linestyle='--',
    linewidth=3
)


# Tick parameters for top plot
ax_top.tick_params(
    axis="both",
    which="major",
    labelsize=22,
    direction="in",
    right=True,
    top=True,
    size=8,
)
ax_top.tick_params(
    axis="both",
    which="minor",
    labelsize=22,
    direction="in",
    right=True,
    top=True,
    size=4,
)
# Adjusted bbox_to_anchor to bring legend inside the plot
ax_top.legend(loc='upper left', bbox_to_anchor=(0.4, 0.99), prop={'size': 16})

# --- Bottom Subplot (SVMin plot, now axs[1]) ---
ax_bottom = axs[1]
ax_bottom.set_xlabel(r"$\Gamma_{\text{fold}}$ (fm)", fontdict=fontDict) 
ax_bottom.set_xlim(-0.02, 1.22)
ax_bottom.set_ylim(53.7, 58.6)
ax_bottom.set_yticks(np.arange(54, 59, 1))
ax_bottom.set_xticks(np.arange(0.0, 1.3, 0.2)) # Set x-ticks at 0.2 intervals
ax_bottom.set_ylabel(r"quality measure $\chi^2$", fontdict=fontDict, labelpad=10)

# Plotting Data Block 0 (SVmin_HFB_cut30) with interpolation - Red line
x_raw_0 = data_blocks["SVmin_HFB_cut30"][:, 0]
y_raw_0 = data_blocks["SVmin_HFB_cut30"][:, 1]
f_interp_0 = interp1d(x_raw_0, y_raw_0, kind='cubic')
x_smooth_0 = np.linspace(x_raw_0.min(), x_raw_0.max(), 500)
y_smooth_0 = f_interp_0(x_smooth_0)
ax_bottom.plot(
    x_smooth_0, y_smooth_0,
    label=r"SV-min like, HFB, $\epsilon_{\text{cut}}$=30MeV",
    color='red', # Changed to red
    linestyle='-',
    linewidth=3
)

# Plotting Data Block 1 (SVmin_BCS_cut30) with interpolation - Green line
x_raw_1 = data_blocks["SVmin_BCS_cut30"][:, 0]
y_raw_1 = data_blocks["SVmin_BCS_cut30"][:, 1]
f_interp_1 = interp1d(x_raw_1, y_raw_1, kind='cubic')
x_smooth_1 = np.linspace(x_raw_1.min(), x_raw_1.max(), 500)
y_smooth_1 = f_interp_1(x_smooth_1)
ax_bottom.plot(
    x_smooth_1, y_smooth_1,
    label=r"SV-min like, BCS, $\epsilon_{\text{cut}}$=30MeV",
    color='green', # Changed to green
    linestyle='--',
    linewidth=3
)

# Tick parameters for bottom plot
ax_bottom.tick_params(
    axis="both",
    which="major",
    labelsize=22,
    direction="in",
    right=True,
    top=True,
    size=8,
)
ax_bottom.tick_params(
    axis="both",
    which="minor",
    labelsize=22,
    direction="in",
    right=True,
    top=True,
    size=4,
)
ax_bottom.legend(loc='upper right', prop={'size': 16})

# Align Y-labels
fig.align_ylabels(axs)

plt.tight_layout()

plt.savefig("fit-width_python.pdf", dpi=300, bbox_inches="tight")
plt.show()