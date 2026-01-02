import re
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
from matplotlib import gridspec
import scienceplots
import numpy as np
from pathlib import Path

# Apply the 'science' style for the plot design
plt.style.use('science')

# --- Set global font sizes for all text elements ---
mpl.rcParams['axes.labelsize'] = 25    # font size of the x and y labels
mpl.rcParams['xtick.labelsize'] = 28   # font size of the x tick labels
mpl.rcParams['ytick.labelsize'] = 28   # font size of the y tick labels
mpl.rcParams['legend.fontsize'] = 24   # legend font size
mpl.rcParams['axes.titlesize'] = 22    # font size of the axes title
mpl.rcParams['lines.linewidth'] = 4.0
plt.style.use("science") # Apply scienceplots style first
mpl.rcParams["text.usetex"] = True
mpl.rcParams["axes.linewidth"] = 1.5
mpl.rcParams["axes.formatter.useoffset"] = False

# Increased font sizes and using numeric fontweight for maximum boldness
# 900 typically corresponds to 'black' or 'heavy'
fontDict = {"family": "helvetica", "fontsize": 32, "fontweight": 900}
fontDictTitle = {"family": "helvetica", "fontsize": 32, "fontweight": 900}

# Parameters
widths = ["0.0d0", "1.0d0", "2.0d0", "3.0d0"]
cutoffs = ["Cut25"]
methods = ["hfb"]
colors = ['#FF0000', '#00AA00', '#0000FF', '#AA00AA', '#00CCAA']


# Store data
data = {}

# Helper function for linestyle
def get_linestyle(method):
    return '--' if method == 'bcs' else '-'


# Load data
for width in widths:
    for cutoff in cutoffs:
        for method in methods:
            filename = f"collect_cat_box_paper_new/collect_cat_SVbas_unknown_{width}_{method}_{cutoff}.dat"
            path = Path(filename)
            if not path.exists():
                print(f"File not found: {path}")
                continue

            key = (method, width, cutoff)
            data[key] = {"ngrid": [], "E": [], "r_C": []}

            with open(path) as f:
                lines = f.readlines()
                for i in range(0, len(lines) - 1, 2):
                    header = lines[i]
                    match = re.search(r'ngrid=(\d+)', header)
                    if not match:
                        continue
                    ngrid = int(match.group(1))

                    values = lines[i + 1].strip().split()#008000
                    if len(values) >= 4:
                        r_C = float(values[3])
                        E = float(values[4])
                        data[key]["ngrid"].append(ngrid)
                        data[key]["r_C"].append(r_C)
                        data[key]["E"].append(E)

# Plotting section using gridspec
with plt.style.context(['science']):
    # Use `GridSpec` to define the subplot layout
    # fig = plt.figure(figsize=(10, 15))
    # gs = gridspec.GridSpec(2, 1, figure=fig, hspace=0)
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10, 15), sharex=True, gridspec_kw={'hspace': 0})

    # Top plot: E vs ngrid for various widths at Cut25
    # ax0 = fig.add_subplot(gs[0, 0])
    lines0 = {}
    for i, width in enumerate(widths):
        for method in ["hfb", "bcs"]:
            key = (method, width, "Cut25")
            if key in data:
                label = rf'$\Gamma_{{\text{{fold}}}} = {width[:-2]}\ \mathrm{{fm}}$ ({method.upper()})'
                line, = ax0.plot(np.array(data[key]["ngrid"])*0.3, data[key]["r_C"], colors[i], linewidth=3.0, linestyle=get_linestyle(method), label=label)
                # line, = ax0.plot(np.array(data[key]["ngrid"])*0.3, data[key]["E"], colors[i], linewidth=3.0, linestyle=get_linestyle(method),  label=label)
                lines0[(method, width)] = line

    #ax0.set_title(r"$\varepsilon_{\text{cut}}$ = 25 MeV", fontdict=fontDictTitle)
    ax0.set_xlim(14, 40)
    ax0.set_ylabel("charge radius (fm)", fontdict=fontDict)
    ax0.text(0.1, 0.05, r"(a)", fontsize=40, ha='center', va='bottom', transform=ax0.transAxes)

    # ax0.set_ylim(-1021, -1019)
    
    # Bottom plot: r_C vs ngrid for various widths at Cut25
    # ax1 = fig.add_subplot(gs[1, 0])
    lines1 = {}
    for i, width in enumerate(widths):
        for method in ["hfb", "bcs"]:
            key = (method, width, "Cut25")
            if key in data:
                label = rf'$\Gamma_{{\text{{fold}}}} = {width[:-2]}\ \mathrm{{fm}}$ ({method.upper()})'
                line, = ax1.plot(np.array(data[key]["ngrid"])*0.3, data[key]["E"], colors[i], linewidth=3.0, linestyle=get_linestyle(method),  label=label)
                #line, = ax1.plot(np.array(data[key]["ngrid"])*0.3, data[key]["r_C"], colors[i], linewidth=3.0, linestyle=get_linestyle(method), label=label)
                lines1[(method, width)] = line

    ax1.set_xlim(14, 40)
    # ax1.set_ylim(4.640, 4.64399)

    ax1.set_xlabel("box size (fm)", fontdict=fontDict)
    ax1.set_ylabel("total energy (MeV)", fontdict=fontDict)
    ax1.text(0.1, 0.85, r"(b)", fontsize=40, ha='center', va='bottom', transform=ax1.transAxes)

    # Apply tick and label settings to both axes
    for ax in [ax0, ax1]:
        ax.tick_params(
            axis="both",
            which="major",
            labelsize=28,
            labelright=False,
            direction="in",
            right=True,
            top=True,
            size=8,
        )
        ax.tick_params(
            axis="both",
            which="minor",
            labelsize=28,
            labelright=False,
            direction="in",
            right=True,
            top=True,
            size=4,
        )

    # Remove x-tick labels from the top panel
    # ax0.set_xticklabels([])

    # Define the desired order for the legend
    desired_order = [
        ("bcs", "0.0d0"),
        ("hfb", "0.0d0"),
        ("bcs", "1.0d0"),
        ("hfb", "1.0d0"),
        #("bcs", "1.5d0"),
        #("hfb", "1.5d0"),
        ("bcs", "2.0d0"),
        ("hfb", "2.0d0"),
        ("bcs", "3.0d0"),
        ("hfb", "3.0d0"),
        ("bcs", "3.5d0"),
        ("hfb", "3.5d0"),
    ]

    # Create handles and labels for the legend in the desired order
    handles = []
    labels = []
    for method, width in desired_order:
        key = (method, width)
        if key in lines1:
            handles.append(lines1[key])
            label = rf'$\tilde{{R}}_{{\mathcal{{F}}}} = {width[:-2]}$'
            labels.append(label)

    # Add the legend to the bottom plot (ax1)
    ax1.legend(handles=handles, labels=labels, prop={'size': 24})

    # Align the y-labels
    fig.align_ylabels([ax0, ax1])
    plt.tight_layout(pad=1.0, h_pad=0.1)

    plt.savefig("panel_plots_box_Sn.pdf", dpi=300, bbox_inches="tight")
    # plt.show()