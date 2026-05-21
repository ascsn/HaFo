from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scienceplots
from matplotlib import gridspec

# --- 1. Parse Data ---
data = {
    20: np.array([
    [20, 2, 408, -1374.529, 5.240], [20, 4, 776, -1374.537, 5.240],
    [20, 6, 1322, -1374.535, 5.240], [20, 8, 2068, -1374.535, 5.240],
    [20, 10, 3050, -1374.535, 5.240], [20, 12, 4306, -1374.523, 5.240],
    [20, 14, 5844, -1374.530, 5.240]
]),
    40: np.array([
    [40, 2, 408, -1374.2818, 5.2902], [40, 4, 776, -1374.5733, 5.2916],
    [40, 6, 1322, -1374.5795, 5.2916], [40, 8, 2068, -1374.5838, 5.2915],
    [40, 10, 3050, -1374.5778, 5.2915], [40, 12, 4306, -1374.5721, 5.2914],
    [40, 14, 5844, -1374.5608, 5.2914]
]),
    60: np.array([
    [60, 2, 408, -1374.2187, 5.2890], [60, 4, 776, -1374.5150, 5.2917],
    [60, 6, 1322, -1374.6366, 5.2921], [60, 8, 2068, -1374.6449, 5.2920],
    [60, 10, 3050, -1374.6171, 5.2919], [60, 12, 4306, -1374.5585, 5.2918],
    [60, 14, 5844, -1374.5528, 5.2918]
])
}

# --- 2. Matplotlib Setup ---
mpl.use("Agg") # Use Agg backend for headless environment
plt.style.use("science") # Ensure 'scienceplots' is installed for this to work
mpl.rcParams["text.usetex"] = True # Explicitly set to False for this environment
mpl.rcParams["axes.linewidth"] = 1.5
mpl.rcParams["axes.formatter.useoffset"] = False
# mpl.rcParams["font.family"] = "DejaVu Sans" # Set a common sans-serif font
mpl.rcParams['axes.labelsize'] = 25    # font size of the x and y labels
mpl.rcParams['xtick.labelsize'] = 28   # font size of the x tick labels
mpl.rcParams['ytick.labelsize'] = 28   # font size of the y tick labels
mpl.rcParams['legend.fontsize'] = 24   # legend font size
mpl.rcParams['axes.titlesize'] = 22    # font size of the axes title
mpl.rcParams['lines.linewidth'] = 4.0


# Font settings
fontDict = {"family": "helvetica", "fontsize": 32, "fontweight": 900}
fontDictTitle = {"family": "helvetica", "fontsize": 32, "fontweight": 900}

# Colors matching gnuplot styles
colors = {
    1: "#FF0000", # red
    2: "#00AA00", # green
    3: "#0000FF", # blue
}

# --- 3. Figure and Subplots ---
fig = plt.figure(figsize=(10, 15))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 1], figure=fig, hspace=0)
# plt.subplots_adjust(hspace=0)

# --- 4. Bottom Subplot (Total Energy) ---
ax_bottom = fig.add_subplot(gs[1, 0])
ax_bottom.set_xlabel(r"$N_{\text{st}}^{1/3}$", fontdict=fontDict, labelpad=0)
ax_bottom.set_ylabel(r"$\Delta E$ (MeV)", fontdict=fontDict, labelpad=10)
ax_bottom.set_xlim(7, 19)
ax_bottom.set_ylim(-0.1, 0.37)

# Plotting with interpolation and markers
series_data_bottom = {
    20: {'raw': data[20], 'offset_E': 1374.530, 'color_key': 1},
    40: {'raw': data[40], 'offset_E': 1374.5608, 'color_key': 2},
    60: {'raw': data[60], 'offset_E': 1374.5528, 'color_key': 3},
}

for cut_val, s_data in series_data_bottom.items():
    x_raw = s_data['raw'][:, 2] ** (1/3)
    y_raw = s_data['raw'][:, 3] + s_data['offset_E']
    # f_interp = interp1d(x_raw, y_raw, kind='cubic')
    # x_smooth = np.linspace(x_raw.min(), x_raw.max(), 500)

    ax_bottom.plot(
        x_raw, y_raw,
        label=f"{cut_val} MeV",
        color=colors[s_data['color_key']], linestyle='-', linewidth=4
    )
    # ax_bottom.plot(
    #     x_raw, y_raw,
    #     marker='o', markersize=8, linestyle='None',
    #     color=colors[s_data['color_key']]
    # )

leg = ax_bottom.legend(
    title=r"$\varepsilon_{\rm cut}$",
    title_fontsize=52,
    loc='upper right',
    prop={'size': 24},
    bbox_to_anchor=(0.8, 1.0),  
    # frameon=True,
    # edgecolor='black',
)
leg.get_title().set_fontweight('bold')
ax_bottom.text(0.88, 0.1, r'(b)', transform=ax_bottom.transAxes, fontsize=36, verticalalignment='top')
ax_bottom.tick_params(
    axis="both",
    which="major",
    labelsize=28, # Increased tick-label size
    direction="in",
    right=True,
    top=True,
    size=8,
)
ax_bottom.tick_params(
    axis="both",
    which="minor",
    labelsize=28, # Increased tick-label size
    direction="in",
    right=True,
    top=True,
    size=4,
)

# --- 5. Top Subplot (Charge Radius) ---
ax_top = fig.add_subplot(gs[0, 0], sharex=ax_bottom)
ax_top.set_xlabel("", fontdict=fontDict, labelpad=0) # X-label removed for top plot due to sharex
ax_top.set_ylabel(r"$\Delta R_{\rm{ch}}$ ($10^{-3}$ fm)", fontdict=fontDict, labelpad=10)
ax_top.set_xlim(7, 19)
ax_top.set_ylim(-2.95, 0.5)
# ax_top.set_title(r"$^{170}$Yb, SV-bas, HFB", fontdict=fontDictTitle, pad=10)
ax_top.text(0.25, 0.35, r'$^{170}$Yb', transform=ax_top.transAxes, fontsize=72, verticalalignment='top', fontweight='bold')
ax_top.text(0.2, 0.2, r'SV-bas, HFB', transform=ax_top.transAxes, fontsize=54, verticalalignment='top', fontweight='bold')
ax_top.text(0.88, 0.1, r'(a)', transform=ax_top.transAxes, fontsize=36, verticalalignment='top', fontweight='bold')
# Make x-axis tick labels invisible for the top plot
plt.setp(ax_top.get_xticklabels(), visible=False)

# Plotting with interpolation and markers
series_data_top = {
    20: {'raw': data[20], 'offset_RC': 5.240, 'color_key': 1},
    40: {'raw': data[40], 'offset_RC': 5.2914, 'color_key': 2},
    60: {'raw': data[60], 'offset_RC': 5.2918, 'color_key': 3},
}

for cut_val, s_data in series_data_top.items():
    x_raw = s_data['raw'][:, 2] ** (1/3)
    y_raw = (s_data['raw'][:, 4] - s_data['offset_RC']) * 1e3
    # f_interp = interp1d(x_raw, y_raw, kind='cubic')
    # x_smooth = np.linspace(x_raw.min(), x_raw.max(), 500)

    ax_top.plot(
        x_raw, y_raw,
        label=rf"$\epsilon_{{\text{{cut}}}}$ = {cut_val} MeV",
        color=colors[s_data['color_key']], linestyle='-', linewidth=4
    )
    # ax_top.plot(
    #     x_raw, y_raw,
    #     marker='o', markersize=8, linestyle='None',
    #     color=colors[s_data['color_key']]
    # )

# ax_top.legend(loc='lower right', prop={'size': 20})
ax_top.tick_params(
    axis="both",
    which="major",
    labelsize=28, # Increased tick-label size
    direction="in",
    right=True,
    top=True,
    size=8,
)
ax_top.tick_params(
    axis="both",
    which="minor",
    labelsize=28, # Increased tick-label size
    direction="in",
    right=True,
    top=True,
    size=4,
)

# --- Final Touches ---
fig.align_ylabels([ax_bottom, ax_top]) # Align y-labels
plt.tight_layout()  # Reintroduce tight_layout

# Save output
plt.savefig("trend-deocc-170Yb_3_python.pdf", dpi=300, bbox_inches="tight")