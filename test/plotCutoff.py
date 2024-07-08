import numpy as np
import matplotlib.pyplot as plt
import os

dataPath = "./collect_cat"
cutoff = np.arange(18, 121)


def make_plots(energy, legends):
    fig, ax = plt.subplots(dpi=150)
    for i, val in enumerate(energy):
        ax.plot(cutoff, energy[i], label=f"{legends[i]}")

    ax.set_xlabel(rf"$\epsilon_{{cut}}$")
    ax.set_ylabel("Total Energy (MeV)")
    plt.legend()
    plt.savefig("TotalEnergy.pdf")
    plt.close()


def main():
    # Read data from collect_cat files
    hfbwidths = ["0.0", "1.0", "2.0"]
    totalEnergyArray = []
    for file in os.listdir(dataPath):
        totalEnergy = []
        if "collect" in file:
            with open(dataPath + "/" + file) as f:
                for line in f:
                    dataLine = line.rstrip().split()
                    if len(dataLine) > 1:
                        totalEnergy.append(float(dataLine[4]))
                totalEnergy = np.array(totalEnergy)
        if len(totalEnergy) > 0:
            totalEnergyArray.append(totalEnergy)
    totalEnergyArray = np.array(totalEnergyArray)

    make_plots(totalEnergyArray, hfbwidths)


if __name__ == "__main__":
    main()
