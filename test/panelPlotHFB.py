"""
This code is used to generate panel plots by reading data from the collect_cat_ files.
The code reads data using read_data function, and then plots the data using make_plots
function. The preprocess_line function is used to preprocess data that is not separated 
by space.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, re
import scienceplots


dataPath = "./collect_cat_box"
# cutoff = np.arange(15, 16)


plt.style.use("science")
mpl.rcParams["text.usetex"] = True
mpl.rcParams["axes.linewidth"] = 1.5
mpl.rcParams["axes.formatter.useoffset"] = False

colors = [
    "r",
    "b",
    "g",
    "k",
    "purple",
    "m",
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "orange",
    "y",
]


def preprocess_line(line):
    # Add a space before negative signs that are not part of scientific notation
    line = re.sub(r"(\d)(-)", r"\1 -", line)
    # Ensure that six decimal places are not split by adding spaces correctly
    line = re.sub(r"(-?\d\.\d{6})(?=-?\d)", r"\1 ", line)
    return line


def read_data(pairingType=None, functional=None):
    # Read data from collect_cat files
    fileList = []
    hfbwidths = []
    totalEnergyArray = []
    chargeRadiusArray = []
    pairingEnergyArray = []
    neutronGapArray = []
    for file in os.listdir(dataPath):
        fileList.append(file)
    # Sort the files in increasing order of either hfbwidth or box size.
    fileList = sorted(fileList)
    for file in fileList:
        totalEnergy = []
        chargeRadius = []
        pairingEnergy = []
        neutronGap = []
        # Checks the name of the file to contain pairing type, functional name,
        # element name and whether a cutoff is set
        if (
            pairingType in file
            and functional in file
            and "120Sn" in file
            and "Cut100" in file
        ):
            print("Pair type", pairingType)
            # Create an array of elements separated by "_" in a filename.
            fileName = file.rstrip().split("_")
            print(fileName)
            # Split the element in the fileName array separated by "d"
            # hfbwidths.append(fileName[4].split("d")[0])
            # Split the element in the fileName array separated by "." and round
            # to the nearest integer
            hfbwidths.append(round(float(fileName[-1].split(".")[0]) * 0.3))
            with open(dataPath + "/" + file) as f:
                data = f.read()
            for line in data.splitlines():
                if line.startswith(" Z,N,r_C,E,epair,gap,sumuv="):
                    dataLine = preprocess_line(line).split()
                else:
                    dataLine = line.rstrip().split()
                # Check if size of array dataLine > 1
                if len(dataLine) > 1:
                    totalEnergy.append(float(dataLine[4]))
                    chargeRadius.append(float(dataLine[3]))
                    pairingEnergy.append(float(dataLine[6]))
                    neutronGap.append(float(dataLine[8]))
            totalEnergy = np.array(totalEnergy)
            chargeRadius = np.array(chargeRadius)
            pairingEnergy = np.array(pairingEnergy)
            neutronGap = np.array(neutronGap)
        # Collect all arrays into one big array
        if len(totalEnergy) > 0:
            totalEnergyArray.append(totalEnergy)
            chargeRadiusArray.append(chargeRadius)
            pairingEnergyArray.append(pairingEnergy)
            neutronGapArray.append(neutronGap)
    totalEnergyArray = np.array(totalEnergyArray)
    chargeRadiusArray = np.array(chargeRadiusArray)
    pairingEnergyArray = np.array(pairingEnergyArray)
    neutronGapArray = np.array(neutronGapArray)
    print(totalEnergyArray.shape)
    # totalEnergyArray = totalEnergyArray[totalEnergyArray.astype(bool)]
    return (
        hfbwidths,
        totalEnergyArray,
        chargeRadiusArray,
        pairingEnergyArray,
        neutronGapArray,
    )


def make_plots(data, labels):
    (
        hfbwidth,
        totalEnergyArray,
        chargeRadiusArray,
        pairingEnergyArray,
        neutronGapArray,
    ) = data
    print(hfbwidth)
    labelpadding = [0, 20, 0, 15]
    rows = [0, 1, 0, 1]
    cols = [0, 0, 1, 1]
    # yrange = [[-1037, -1013], [4.618, 4.66], [-99, 10], [-0.2, 2.8]]
    fig, ax = plt.subplots(nrows=2, ncols=2, sharex=False, figsize=(16, 12), dpi=600)
    for j in range(4):
        dict1 = labels[j]
        xl = dict1["xlabel"]
        yl = dict1["ylabel"]
        fl = dict1["filename"]
        ttl = dict1["title"]
        newData = []
        for i, val in enumerate(data[j + 1]):
            newData.append(val)
        ax[rows[j], cols[j]].plot(
            hfbwidth,
            newData,
            linestyle="-",
            marker="o",
            label=f"box size = {hfbwidth[i]} fm",
            linewidth=3,
            color="b",
        )

        # ax[j].set_ylim(yrange[j])
        ax[rows[j], cols[j]].set_xlabel(xl, fontsize=18)
        ax[rows[j], cols[j]].set_ylabel(yl, labelpad=labelpadding[j], fontsize=18)
        ax[rows[j], cols[j]].tick_params(
            axis="both",
            which="major",
            labelsize=16,
            labelright=False,
            direction="in",
            right=True,
            top=True,
            size=8,
        )
        ax[rows[j], cols[j]].tick_params(
            axis="both",
            which="minor",
            labelsize=16,
            labelright=False,
            direction="in",
            right=True,
            top=True,
            size=4,
        )
    # fig.text(0.5, 0.07, r"$\epsilon_{{cut}}$ [MeV]", ha="center", fontsize=16)
    plt.subplots_adjust(hspace=0)
    fig.text(
        0.5,
        0.9,
        r"$^{120}$Sn, SV-bas, spherical 1D, vary box size",
        ha="center",
        fontsize=18,
    )
    # ax[0].set_title(r"$^{120}$Sn", fontsize=18)
    # plt.legend(loc="lower right", prop={"size": 8})
    fileName = "SVbas_120Sn_hfb_props_Diff_box.pdf"
    plt.savefig(fileName)
    plt.close()


def main():
    pairingType = ["hfb"]
    eltName = r"$^{120}$Sn"
    functional = "SVbas"
    for pairtype in pairingType:
        data = read_data(pairtype, functional)
        # Total Energy plot
        labels1 = {
            "xlabel": r"box size [fm]",
            "ylabel": r"Total Energy [MeV]",
            "filename": f"TotalEnergy_{pairtype}_{eltName}_{functional}_Cut100.pdf",
            "title": f"{eltName}",
        }

        # Radii plot
        labels2 = {
            "xlabel": r"box size [fm]",
            "ylabel": r"Charge Radius [fm]",
            "filename": rf"Radius_{pairtype}_{eltName}_{functional}_Cut100.pdf",
            "title": f"{eltName}",
        }

        # Neutron pairing energy plot
        labels3 = {
            "xlabel": r"box size [fm]",
            "ylabel": r"Pairing Energy [MeV]",
            "filename": rf"PairingEnergy_{pairtype}_{eltName}_{functional}_Cut100.pdf",
            "title": eltName,
        }

        # Average neutron pairing gap plot
        labels4 = {
            "xlabel": r"box size [fm]",
            "ylabel": r"Pairing Gap [MeV]",
            "filename": rf"PairingGap_{pairtype}_{eltName}_{functional}_Cut100.pdf",
            "title": eltName,
        }

        labels = [labels1, labels2, labels3, labels4]
        make_plots(data, labels)

        os.system(rf"mv *_{pairtype}_*.pdf plots_paper/")


if __name__ == "__main__":
    main()
