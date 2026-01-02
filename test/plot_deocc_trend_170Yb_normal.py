import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from pathlib import Path
import scienceplots

# --- Data Parsing from the provided snippet ---
# Gnuplot indices mapping to cuts:
# i 0 -> cut = 40
# i 1 -> cut = 60
# i 3 -> cut = 20

# Manually extracted and structured data based on cut values and columns
# Columns: [Delta_shell (col 2), E_B (col 8), r_C (col 13)]
data = {
    20: np.array([ # This is gnuplot index 3
        [2, -1374.529, 5.240],
        [4, -1374.537, 5.240],
        [6, -1374.535, 5.240],
        [8, -1374.535, 5.240],
        [10, -1374.535, 5.240],
        [12, -1374.523, 5.240],
        [14, -1374.530, 5.240]
    ]),
    40: np.array([ # This is gnuplot index 0
        [2, -1374.2818, 5.2902],
        [4, -1374.5733, 5.2916],
        [6, -1374.5795, 5.2916],
        [8, -1374.5838, 5.2915],
        [10, -1374.5778, 5.2915],
        [12, -1374.5721, 5.2914],
        [14, -1374.5608, 5.2914]
    ]),
    60: np.array([ # This is gnuplot index 1
        [2, -1374.2187, 5.2890],
        [4, -1374.5150, 5.2917],
        [6, -1374.6366, 5.2921],
        [8, -1374.6449, 5.2920],
        [10, -1374.6171, 5.2919],
        [12, -1374.5585, 5.2918],
        [14, -1374.5528, 5.2918]
    ])
}

# --- Matplotlib Setup ---
mpl.use("Agg") # Use Agg backend for non-interactive plotting
plt.style.use("science") # Removed as it caused an error in the previous turn
mpl.rcParams["text.usetex"] = True
mpl.rcParams["axes.linewidth"] = 1.5
mpl.rcParams["axes.formatter.useoffset"] = False # Matches 'set form x' for default formatting

# Font settings (Gnuplot "TimesNewRoman,30" -> Matplotlib fontdict)
fontDict = {"family": "Times New Roman", "fontsize": 30, "fontweight": "normal"}
fontDictTitle = {"family": "Times New Roman", "fontsize": 36, "fontweight": "normal"}

# Gnuplot line styles to Matplotlib properties
# pt 7 corresponds to 'o' (circle) marker in Matplotlib
# ps 0.8 is marker size
plot_styles = {
    "ls 1": {'color': '#FF0000', 'lw': 4, 'marker': '', 'ms': 8, 'ls': '-'}, # red
    "ls 2": {'color': '#00AA00', 'lw': 4, 'marker': '', 'ms': 8, 'ls': '-'}, # green
    "ls 3": {'color': '#0000FF', 'lw': 4, 'marker': '', 'ms': 8, 'ls': '-'}, # blue
}

# Figure and Subplots (2 rows, 1 column)
# Gnuplot 'set size 1.03,1.55' suggests a tall plot.
# A ratio of 1.55/1.03 approx 1.5. For width 8, height ~12.
# Using (8,14) as a good previous example for two vertical plots
fig, axs = plt.subplots(2, 1, figsize=(10, 15), constrained_layout=True, sharex=True)

# --- Top Subplot: Charge Radius ---
# Corresponds to 'set origin 0.2,0.8' in gnuplot
ax_top = axs[0]
ax_top.set_xlim(1, 15)
ax_top.set_ylim(-2.95, 0.5)
ax_top.set_ylabel(r"charge radius $\left(10^{-3}\mathrm{fm}\right)$", fontdict=fontDict, labelpad=10) # offset 0.5
ax_top.set_xlabel("", fontdict=fontDict) # Keep xlabel empty for top plot
ax_top.set_title(r"$^{170}\mathrm{Yb}$, SV-bas, HFB", fontdict=fontDictTitle, y=1.0 - (0.3/ax_top.get_figure().get_size_inches()[1]*14)) # Adjust y for offset -0.3

# Plotting data for top subplot (charge radius)
# i 3 -> cut=20 (ls 1)
x_20_rc = data[20][:, 0]
y_20_rc = (data[20][:, 2] - 5.240) * 1e3
ax_top.plot(x_20_rc, y_20_rc, label=r"$\epsilon_{\text{cut}}$ = 20 MeV", **plot_styles["ls 1"])

# i 0 -> cut=40 (ls 2)
x_40_rc = data[40][:, 0]
y_40_rc = (data[40][:, 2] - 5.2914) * 1e3
ax_top.plot(x_40_rc, y_40_rc, label=r"$\epsilon_{\text{cut}}$ = 40 MeV", **plot_styles["ls 2"])

# i 1 -> cut=60 (ls 3)
x_60_rc = data[60][:, 0]
y_60_rc = (data[60][:, 2] - 5.2918) * 1e3
ax_top.plot(x_60_rc, y_60_rc, label=r"$\epsilon_{\text{cut}}$ = 60 MeV", **plot_styles["ls 3"])

ax_top.legend(loc='lower right', prop={'size': 24}) # set key bottom right

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

# --- Bottom Subplot: Total Energy ---
# Corresponds to 'set origin 0.2,0.2' in gnuplot
ax_bottom = axs[1]
ax_bottom.set_xlim(1, 15)
ax_bottom.set_ylim(-0.1, 0.37)
ax_bottom.set_xlabel(r"$\Delta_{\text{shell}}$", fontdict=fontDict, labelpad=0) # offset 0
ax_bottom.set_ylabel("total energy (MeV)", fontdict=fontDict, labelpad=10) # offset 1

# Plotting data for bottom subplot (total energy)
# i 3 -> cut=20 (ls 1)
x_20_eb = data[20][:, 0]
y_20_eb = data[20][:, 1] + 1374.530
ax_bottom.plot(x_20_eb, y_20_eb, label=r"$\epsilon_{\text{cut}}$ = 20 MeV", **plot_styles["ls 1"])

# i 0 -> cut=40 (ls 2)
x_40_eb = data[40][:, 0]
y_40_eb = data[40][:, 1] + 1374.5608
ax_bottom.plot(x_40_eb, y_40_eb, label=r"$\epsilon_{\text{cut}}$ = 40 MeV", **plot_styles["ls 2"])

# i 1 -> cut=60 (ls 3)
x_60_eb = data[60][:, 0]
y_60_eb = data[60][:, 1] + 1374.5528
ax_bottom.plot(x_60_eb, y_60_eb, label=r"$\epsilon_{\text{cut}}$ = 60 MeV", **plot_styles["ls 3"])

ax_bottom.legend(loc='upper right', prop={'size': 24}) # set key top right

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

# Align Y-labels (helpful for multi-panel plots)
fig.align_ylabels(axs)

# Save the plot
plt.savefig("trend-deocc-170Yb_3_python.pdf", dpi=300, bbox_inches="tight")
plt.show()