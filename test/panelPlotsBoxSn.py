import re
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import scienceplots
import numpy as np


mpl.use("TkAgg")
plt.style.use("science")
mpl.rcParams["text.usetex"] = True
mpl.rcParams["axes.linewidth"] = 1.5
mpl.rcParams["axes.formatter.useoffset"] = False

fontDict = {"family": "helvetica", "fontsize": 18, "fontweight": "bold"}
fontDictTitle = {"family": "helvetica", "fontsize": 18, "fontweight": "bold"}


# Parameters
widths = ["0.0d0", "1.5d0", "3.0d0"]
cutoffs = ["Cut25", "Cut100"]
data = {(w, c): {"ngrid": [], "E": [], "r_C": []} for w in widths for c in cutoffs}


def split_line(line):
    match = re.match(r"([a-z]+)([0-9]+)", line, re.I)
    if match:
        items = match.groups()
    return items[-1]


# Load data
for width in widths:
    for cutoff in cutoffs:
        filename = f"collect_cat_box_paper_new/collect_cat_SVbas_unknown_{width}_hfb_{cutoff}.dat"
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()
                for i in range(0, len(lines), 2):
                    # First line: extract ngrid
                    header = lines[i].strip()
                    if "ngrid=" not in header:
                        continue
                    match = re.search(r'ngrid=(\d+)', header)
                    if not match:
                        continue
                    ngrid = int(match.group(1))

                    # Second line: extract r_C and E
                    fields = lines[i+1].strip().split()
                    if len(fields) >= 4:
                        r_C = float(fields[3])
                        E = float(fields[4])
                        data[(width, cutoff)]["ngrid"].append(ngrid)
                        data[(width, cutoff)]["r_C"].append(r_C)
                        data[(width, cutoff)]["E"].append(E)
        except FileNotFoundError:
            print(f"File not found: {filename}")

# Plot
fig, axs = plt.subplots(2, 2, figsize=(12, 10), sharex=True)

# Remove space between plots
plt.subplots_adjust(hspace=0, wspace=0)

# Panel 1: E vs ngrid for widths at Cut25
ax = axs[0, 0]
for width in widths:
    key = (width, "Cut25")
    if data[key]["ngrid"]:
        ax.plot(np.array(data[key]["ngrid"])*0.3, data[key]["E"], label=f'Width = {width[:-2]}', marker='o')
# ax.set_title("Energy vs ngrid (Cutoff = 25)")
# ax.set_xlabel("ngrid")
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
ax.set_ylabel("Energy [MeV]", fontdict=fontDict)
ax.legend()

# Panel 2: r_C vs ngrid for widths at Cut25
ax = axs[1, 0]
for width in widths:
    key = (width, "Cut25")
    if data[key]["ngrid"]:
        ax.plot(np.array(data[key]["ngrid"])*0.3, data[key]["r_C"], label=rf'$\Gamma_{{\text{{fold}}}} = ${width[:-2]}', marker='o')
# ax.set_title("r_C vs ngrid (Cutoff = 25)")
ax.set_xlabel("Box Size [fm]", fontdict=fontDict)
ax.set_ylabel("Charge Radius [fm]", fontdict=fontDict, labelpad=20)
ax.set_xlim(14, 40)
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
ax.text(0.05, 0.75, "(c)", transform=ax.transAxes, fontsize=16)
ax.legend()

# Panel 3: E vs ngrid for Cut25 and Cut100 at width = 1.5d0
ax = axs[0, 1]
for cutoff in cutoffs:
    key = ("1.5d0", cutoff)
    if cutoff == "Cut25":
        cutoff = "25"
    elif cutoff == "Cut100":
        cutoff = "100"
    else:
        continue
    if data[key]["ngrid"]:
        ax.plot(np.array(data[key]["ngrid"])*0.3, data[key]["E"], label=f'Cutoff = {cutoff}', marker='o')
ax.set_title(r"$\Gamma_{\text{fold}} = 1.5$ fm", fontdict=fontDictTitle)
# ax.set_xlabel("ngrid")
ax.set_xlim(14, 40)
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
ax.text(0.05, 0.75, "(b)", transform=ax.transAxes, fontsize=16)
ax.legend()

# Panel 4: r_C vs ngrid for Cut25 and Cut100 at width = 1.5d0
ax = axs[1, 1]
for cutoff in cutoffs:
    key = ("1.5d0", cutoff)
    if cutoff == "Cut25":
        cutoff = "25"
    elif cutoff == "Cut100":
        cutoff = "100"
    else:
        continue
    if data[key]["ngrid"]:
        ax.plot(np.array(data[key]["ngrid"])*0.3, data[key]["r_C"], label=f'Cutoff = {cutoff}', marker='o')
# ax.set_title("r_C vs ngrid (Width = 1.5)")
ax.set_xlabel("Box Size [fm]", fontdict=fontDict)
ax.set_xlim(14, 40)
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
ax.text(0.05, 0.75, "(d)", transform=ax.transAxes, fontsize=16)
# ax.set_ylabel("Charge Radius (r_C)")
ax.legend()

plt.tight_layout()
plt.savefig("panel_plots_box_Sn.pdf", dpi=600)
# plt.show()
