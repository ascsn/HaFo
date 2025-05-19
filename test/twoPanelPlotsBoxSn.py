import re
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import scienceplots
import numpy as np
from pathlib import Path


mpl.use("TkAgg")
plt.style.use("science")
mpl.rcParams["text.usetex"] = True
mpl.rcParams["axes.linewidth"] = 1.5
mpl.rcParams["axes.formatter.useoffset"] = False

fontDict = {"family": "helvetica", "fontsize": 22, "fontweight": "bold"}
fontDictTitle = {"family": "helvetica", "fontsize": 18, "fontweight": "bold"}

# Parameters
widths = ["0.0d0", "1.5d0", "3.0d0"]
cutoffs = ["Cut25", "Cut100"]
methods = ["hfb", "bcs"]
colors = ['red', 'blue', 'black', 'purple']


# Store data
data = {}

# Helper function for linestyle
def get_linestyle(method):
    return '--' if method == 'bcs' else '-'


# Load data
for width in widths:
    for cutoff in cutoffs:
        for method in methods:
            filename = f"collect_cat_box_Paper/collect_cat_SVbas_120Sn_{width}_{method}_{cutoff}.dat"
            path = Path(filename)
            if not path.exists():
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

                    values = lines[i + 1].strip().split()
                    if len(values) >= 4:
                        r_C = float(values[3])
                        E = float(values[4])
                        data[key]["ngrid"].append(ngrid)
                        data[key]["r_C"].append(r_C)
                        data[key]["E"].append(E)

# Plot
fig, axs = plt.subplots(2, 1, figsize=(8, 12), sharex=True, constrained_layout=True)
# fig.subplots_adjust(hspace=0, wspace=0)

# Top-left: E vs ngrid for various widths at Cut25
ax = axs[0]
for i, width in enumerate(widths):
    for method in ["hfb", "bcs"]:
        key = (method, width, "Cut25")
        if key in data:
            label = rf'$\Gamma_{{\text{{fold}}}} = {width[:-2]}\ \mathrm{{fm}}$ ({method.upper()})'
            ax.plot(np.array(data[key]["ngrid"])*0.3, data[key]["E"], colors[i], linewidth=2.0, linestyle=get_linestyle(method),  label=label)
ax.tick_params(
    axis="both",
    which="major",
    labelsize=16,
    labelright=False,
    direction="in",
    right=True,
    top=True,
    size=8,
)
ax.tick_params(
    axis="both",
    which="minor",
    labelsize=16,
    labelright=False,
    direction="in",
    right=True,
    top=True,
    size=4,
)
ax.set_title("Cutoff = 25 MeV", fontdict=fontDictTitle)
ax.text(0.05, 0.75, "(a)", transform=ax.transAxes, fontsize=16)
ax.set_xlim(14, 40)
ax.set_ylim(-1022,-1018.5)
ax.set_ylabel("Energy [MeV]", fontdict=fontDict)
# ax.legend(prop={'size': 12})

# Bottom-left: r_C vs ngrid for various widths at Cut25
ax = axs[1]
for i, width in enumerate(widths):
    for method in ["hfb", "bcs"]:
        key = (method, width, "Cut25")
        if key in data:
            label = rf'$\Gamma_{{\text{{fold}}}} = {width[:-2]}\ \mathrm{{fm}}$ ({method.upper()})'
            ax.plot(np.array(data[key]["ngrid"])*0.3, data[key]["r_C"], colors[i], linewidth=2.0, linestyle=get_linestyle(method), label=label)
ax.set_xlabel("Box Size [fm]", fontdict=fontDict)
ax.set_ylabel("Charge Radius [fm]", fontdict=fontDict, labelpad=10)
ax.set_xlim(14, 40)
ax.set_ylim(4.640, 4.644)
ax.tick_params(
    axis="both",
    which="major",
    labelsize=22,
    labelright=False,
    direction="in",
    right=True,
    top=True,
    size=8,
)
ax.tick_params(
    axis="both",
    which="minor",
    labelsize=22,
    labelright=False,
    direction="in",
    right=True,
    top=True,
    size=4,
)
ax.text(0.05, 0.75, "(b)", transform=ax.transAxes, fontsize=16)
ax.legend(prop={'size': 12})

plt.tight_layout()

plt.savefig("panel_plots_box_Sn.pdf", dpi=600, bbox_inches="tight")
# plt.show()
