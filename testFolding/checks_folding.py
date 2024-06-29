import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

testResultsPath = "./results_data/"
plotResultsPath = "./plots_data/"

colors = [
    "r",
    "b",
    "g",
    "k",
    "orange",
    "purple",
    "y",
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
]


def plot():
    jsList = []
    wpList = []
    fwList = []
    for file in os.listdir(testResultsPath):
        fileSplit = file.strip("test.dat").split("_")
        jsList.append(fileSplit[2])
        wpList.append(fileSplit[-1])
        fwList.append(fileSplit[4])

    jsList = np.asarray(sorted(list(set(jsList))))
    wpList = np.asarray(sorted(list(set(wpList))), dtype=np.float64)
    fwList = np.asarray(sorted(list(set(fwList))), dtype=np.float64)
    print("wplist", wpList)
    for js in jsList:
        for fw in fwList:
            fig, ax = plt.subplots(dpi=150)
            cmap = plt.cm.viridis
            norm = plt.Normalize(wpList[0], wpList[-1])
            colors = cmap(norm(np.asarray(wpList, dtype=np.float64)))
            for i, wp in enumerate(wpList):
                fileName = f"test_js_{js}_fw_{fw}_wp_{wp}.dat"
                data = np.loadtxt(f"{testResultsPath}{fileName}")
                ax.plot(
                    data.T[0],
                    data.T[2],
                    color=colors[i],
                    linestyle="-",
                    label=rf"Folded: {wpList[i]}",
                )
                ax.plot(
                    data.T[0],
                    data.T[3],
                    color=colors[i],
                    linestyle="--",
                    label=rf"Exact: {wpList[i]}",
                )
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax)
            cbar.set_label("wp width")
            # ax.set_yscale("log")
            ax.set_xlabel("r (fm)")
            ax.set_ylabel(r"folded $\varphi$")
            ax.set_title(f"js = {js}, foldwidth = {fw}")
            # plt.legend(prop={"size": 4}, loc="upper right")
            plotFileName = f"Fig_Js_{js}_fw_{fw}.pdf"
            plt.savefig(plotFileName)
            os.system(f"mv {plotFileName} {plotResultsPath}")
            plt.close()
            plt.cla()


if __name__ == "__main__":
    os.system("gfortran -c -w ../source_1D/givens.f")
    os.system("gfortran -c Gauss_folding.f90")
    os.system("gfortran -c k2psi.f")
    os.system(
        "gfortran -o testFolding testFolding.f90 Gauss_folding.o k2psi.o givens.o"
    )

    os.system("./testFolding")
    plot()
