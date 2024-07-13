import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import scienceplots


dataPath = "./collect_cat"
cutoff = np.arange(15, 121)


plt.style.use("science")
mpl.rcParams["text.usetex"] = True
mpl.rcParams["axes.linewidth"] = 1.5

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


def make_plots(data, legends, **kwargs):
    xl = kwargs["xlabel"]
    yl = kwargs["ylabel"]
    fl = kwargs["filename"]
    ttl = kwargs["title"]
    fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
    for i, val in enumerate(data):
        ax.plot(
            cutoff, val, label=f"hfbwidth={legends[i]}", linewidth=2, color=colors[i]
        )

    ax.set_xlabel(r"$\epsilon_{{cut}}$", fontsize=22)
    ax.set_ylabel(yl, fontsize=22)
    ax.set_title(ttl, fontsize=18)
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=18,
        labelright=False,
        direction="in",
        right=True,
        top=True,
        size=8,
    )
    ax.tick_params(
        axis="both",
        which="minor",
        labelsize=18,
        labelright=False,
        direction="in",
        right=True,
        top=True,
        size=4,
    )
    plt.legend(loc="lower left", prop={"size": 10})
    plt.savefig(fl)
    plt.close()


def read_data(pairingType=None):
    # Read data from collect_cat files
    fileList = []
    hfbwidths = []
    totalEnergyArray = []
    chargeRadiusArray = []
    pairingEnergyArray = []
    neutronGapArray = []
    for file in os.listdir(dataPath):
        fileList.append(file)
    fileList = sorted(fileList)
    for file in fileList:
        totalEnergy = []
        chargeRadius = []
        pairingEnergy = []
        neutronGap = []
        if pairingType in file:
            fileName = file.rstrip().split("_")
            hfbwidths.append(fileName[-2])
            with open(dataPath + "/" + file) as f:
                for line in f:
                    dataLine = line.rstrip().split()
                    if len(dataLine) > 1:
                        totalEnergy.append(float(dataLine[4]))
                        chargeRadius.append(float(dataLine[3]))
                        pairingEnergy.append(float(dataLine[6]))
                        neutronGap.append(float(dataLine[8]))
                totalEnergy = np.array(totalEnergy)
                chargeRadius = np.array(chargeRadius)
                pairingEnergy = np.array(pairingEnergy)
                neutronGap = np.array(neutronGap)
        if len(totalEnergy) > 0:
            totalEnergyArray.append(totalEnergy)
            chargeRadiusArray.append(chargeRadius)
            pairingEnergyArray.append(pairingEnergy)
            neutronGapArray.append(neutronGap)
    totalEnergyArray = np.array(totalEnergyArray)
    chargeRadiusArray = np.array(chargeRadiusArray)
    pairingEnergyArray = np.array(pairingEnergyArray)
    neutronGapArray = np.array(neutronGapArray)
    return (
        hfbwidths,
        totalEnergyArray,
        chargeRadiusArray,
        pairingEnergyArray,
        neutronGapArray,
    )


def main():
    pairingType = ["bcs"]
    for pairtype in pairingType:
        (
            hfbwidth,
            totalEnergyArray,
            chargeRadiusArray,
            pairingEnergyArray,
            neutronGapArray,
        ) = read_data(pairtype)
        # Total Energy plot
        labels = {
            "xlabel": r"$\epsilon_{{cut}}$",
            "ylabel": r"Total Energy (MeV)",
            "filename": f"TotalEnergy_{pairtype}.pdf",
            "title": r"Total Energy",
        }
        make_plots(totalEnergyArray, hfbwidth, **labels)

        # Total Energy plot
        labels = {
            "xlabel": r"$\epsilon_{{cut}}$",
            "ylabel": r"Total Energy (MeV)",
            "filename": f"TotalEnergy_{pairtype}.pdf",
            "title": r"$^{120}$Sn",
        }
        make_plots(totalEnergyArray, hfbwidth, **labels)

        # Radii plot
        labels = {
            "xlabel": r"$\epsilon_{{cut}}$",
            "ylabel": r"Charge Radius (fm)",
            "filename": f"Radius_{pairtype}.pdf",
            "title": r"$^{120}$Sn",
        }
        make_plots(chargeRadiusArray, hfbwidth, **labels)

        # Neutron pairing energy plot
        labels = {
            "xlabel": r"$\epsilon_{{cut}}$",
            "ylabel": r"Pairing Energy (MeV)",
            "filename": f"PairingEnergy_{pairtype}.pdf",
            "title": r"$^{120}$Sn",
        }
        make_plots(pairingEnergyArray, hfbwidth, **labels)

        # Average neutron pairing gap plot
        labels = {
            "xlabel": r"$\epsilon_{{cut}}$",
            "ylabel": r"Pairing Gap (MeV)",
            "filename": f"PairingGap_{pairtype}.pdf",
            "title": r"$^{120}$Sn",
        }
        make_plots(neutronGapArray, hfbwidth, **labels)

    os.system(f"mv *.pdf plots/")


if __name__ == "__main__":
    main()
