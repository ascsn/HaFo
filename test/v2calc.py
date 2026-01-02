# Code to do a constrained calculation using HaFo code for 120Sn
# Written: Sudhanva Lalit, August 2025
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess

# Pairing strengths for pairgaps
# foldWidth = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
foldWidth = [3.0]  # in fm
# pairStrengths = [316.5, 4340, 26720]
# pairStrengths = [408.8428, 613.4474, 1567.3926, 4641.2138,11318.6759] # pairgap = 1.866
# pairStrengths = [386.17, 576.945, 1470.08, 4398.35, 10975.81, 20073.04]  # pairgap = 1.610
# pairStrengths = [308, 482.84, 1343.467, 4137.54, 9932.994, 17806.51 ]  # pairgap = 1.31
pairStrengths = [30235.012020]  
# PairTypes
pairTypes = ["bcs", "hfb"]


def main():
    nucleus = "120Sn"
    functional = "SV-bas"
    pairType = pairTypes[1]  # Choose "hfb" for HFB calculation

    v2scan = np.arange(0.0, 0.12, 1e-3)
    for fold, pairStrength in zip(foldWidth, pairStrengths):
        print(
            f"Running calculation for fold width: {fold} fm, pairStrength: {pairStrength} MeV"
        )
        epair_n_vals = []
        gap_n_vals = []
        energyvals = []
        rms_vals = []
        dynamic_gap_vals = []
        for v2 in v2scan:
            print(f"Running scan for v2 = {v2}")
            process = subprocess.run(
                [
                    "sh",
                    "scanv2.bat",
                    f"{fold}",
                    f"{pairStrength}",
                    f"{v2*v2:1.6e}",
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            if process.stderr:
                print("Error:", process.stderr)
            # Scan for011 for output
            with open("for011", "r") as f:
                lines = f.readlines()
                for line in lines:
                    values = line.strip().split()
                    if "Z,N,r_C,E,epair,gap,sumuv=" in values:
                        # print("Values:", values)
                        (
                            Z,
                            N,
                            r_C,
                            E,
                            epair_p,
                            epair_n,
                            gap_p,
                            gap_n,
                            sumuv_p,
                            sumuv_n,
                        ) = map(float, values[1:])
                        print(
                            f"Z: {Z}, N: {N}, r_C: {r_C}, E: {E}, epair_p: {epair_p}, epair_n: {epair_n}, gap_p: {gap_p}, gap_n: {gap_n}, sumuv_p: {sumuv_p}, sumuv_n: {sumuv_n}"
                        )

                    if values and values[-1] == "43":
                        (
                            j,
                            l,
                            energy,
                            weight,
                            evariance,
                            rms,
                            multiplicity,
                            dynamic_gap,
                            wgtherm,
                            x,  # unknown variable
                            y,  # unknown variable
                            z,  # unknown variable
                            number,
                        ) = map(float, values)
                        print(
                            f"2j: {j}, l: {l}, energy: {energy}, weight: {weight}, evariance: {evariance}, rms: {rms}, multiplicity: {multiplicity}, dynamic_gap: {dynamic_gap}, wgtherm: {wgtherm}, x: {x}, y: {y}, z: {z}, number: {number}"
                        )

                        # Collect data for all v2 values
                        epair_n_vals.append(epair_n)
                        gap_n_vals.append(gap_n)
                        energyvals.append(E)
                        rms_vals.append(r_C)
                        dynamic_gap_vals.append(dynamic_gap)
        # Save results to a file
        results_file = f"uvscan_{nucleus}_{functional}.txt"
        with open(results_file, "a") as f:
            f.write(f"#{pairType}, Foldwidth = {fold}, pairStrength = {pairStrength}\n")
            np.savetxt(
                f,
                np.column_stack(
                    (
                        v2scan,
                        epair_n_vals,
                        gap_n_vals,
                        energyvals,
                        rms_vals,
                        dynamic_gap_vals,
                    )
                ),
                header="v2, epair_n, gap_n, energy, rms, dynamic_gap",
                fmt="%1.6e",
                # delimiter="\t",
            )
            f.write("\n")

    # Plot the results


if __name__ == "__main__":
    main()
