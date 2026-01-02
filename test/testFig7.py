# Code to do a constrained calculation using HaFo code for 120Sn
# Written: Sudhanva Lalit, October 2, 2025
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess

# Pairing strengths for pairgaps
# foldWidth = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]
foldWidth = [3.0]  # in fm
pairStrengths = [30235.012020]  # for Rf = 3 SV-bas and 120Sn
# pairStrengths = [1529.696198]  # for Rf= 1 SV-bas and 120Sn
# PairTypes
pairTypes = ["bcs", "hfb"]

# Define grid sizes for the new study (e.g., number of grid points or box size)
grid_sizes = np.arange(51,129, 2) # Example grid sizes (e.g., number of points)

# --- Function to safely convert HaFo output to float, treating '-127' as 'E-127' ---
def safe_float_conversion(value_str):
    """Converts the string value, replacing non-standard exponent formats."""
    # Check for the specific non-standard format 'XXX-YYY' where XXX is float part
    if value_str.count('-') == 1 and not value_str.startswith('-') and 'E' not in value_str and 'D' not in value_str:
        # Example: '1.832375-127' -> '1.832375E-127'
        float_part, exponent_part = value_str.rsplit('-', 1)
        if len(exponent_part) == 3: # Assuming exponent part is always 3 digits (e.g., 127)
             # Reconstruct the string with 'E-' inserted
             value_str = float_part + 'E-' + exponent_part
             
    # Specifically catch the case we want to zero out: 1.832375E-127 (or similar very small numbers)
    # Python can handle E-127 fine, but if you strictly want E-127 to be zero:
    if 'E-' in value_str or 'D-' in value_str:
        if int(value_str.split('E-')[-1] if 'E-' in value_str else value_str.split('D-')[-1]) >= 100:
            return 0.0

    return float(value_str)

# ----------------------------------------------------------------------------------

def main():
    nucleus = "120Sn"
    functional = "SV-bas"
    pairType = pairTypes[1]

    v2scan = np.arange(0.0, 0.12, 1e-3)
    
    optimal_v2_vs_grid = []
    grid_size_for_plot_optimal = []
    all_scan_data = [] # Stores (v2, Energy, grid_size) tuples

    colors = plt.cm.viridis(np.linspace(0, 1, len(grid_sizes)))

    for i, grid_size in enumerate(grid_sizes):
        print(f"\n--- Starting calculations for Grid Size: {grid_size} (Color Index: {i}) ---")
        
        for fold, pairStrength in zip(foldWidth, pairStrengths):
            energyvals = []
            v2_for_energy = []
            
            for v2 in v2scan:
                process = subprocess.run(
                    [
                        "sh",
                        "scanv2.bat",
                        f"{fold}",
                        f"{pairStrength}",
                        f"{v2*v2:1.6e}",
                        f"{grid_size}"
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                )

                if process.stderr:
                    print("Error during HaFo run:", process.stderr)
                
                # Scan for011 for output
                try:
                    with open("for011", "r") as f:
                        lines = f.readlines()
                        current_E = None
                        for line in lines:
                            values = line.strip().split()
                            if "Z,N,r_C,E,epair,gap,sumuv=" in values:
                                # APPLY THE FIX HERE: Use the new conversion function
                                E = safe_float_conversion(values[4]) 
                                current_E = E
                                break 
                        
                        if current_E is not None:
                            energyvals.append(current_E)
                            v2_for_energy.append(v2)
                            all_scan_data.append((v2, current_E, grid_size)) 
                    process = subprocess.run(['cp', 'for011', f"collect_for011_36/for011_{grid_size*0.3}.txt"], check=True)        
                except FileNotFoundError:
                    print(f"Warning: for011 not found for v2={v2}, grid_size={grid_size}. Skipping data collection.")
                except ValueError as e:
                    print(f"Parsing Error for v2={v2}, grid_size={grid_size}. Error: {e}")
                
            # --- Analysis for the current grid_size (for optimal v^2 plot) ---
            if energyvals:
                min_E_index = np.argmin(energyvals)
                optimal_v2 = v2_for_energy[min_E_index]
                min_E = energyvals[min_E_index]

                print(f"Box Size: {grid_size*0.3}, Min Energy: {min_E:1.4f} MeV, Optimal v^2: {optimal_v2:1.6e}")

                grid_size_for_plot_optimal.append(grid_size)
                optimal_v2_vs_grid.append(optimal_v2)
                
            else:
                print(f"Could not find energy data for Grid Size: {grid_size}. Check HaFo logs.")
                
    # ===============================
    # Plot 1: Optimal v^2 vs Grid Size 
    # ===============================
    if grid_size_for_plot_optimal:
        plt.figure(figsize=(10, 6))
        plt.plot(np.array(grid_size_for_plot_optimal)*0.3, optimal_v2_vs_grid, marker='o', linestyle='-')
        
        plt.xlabel("Box Size (fm)")
        plt.ylabel(r"Optimal $v^2$ (Corresponds to Min $E$)")
        plt.title(f"Optimal $v^2$ vs. Box Size for {nucleus} ({functional})")
        plt.grid(True)
        plt.savefig("optimal_v2_vs_box.png")
        print("\nPlot 1 (Optimal v^2 vs Box Size) saved as optimal_v2_vs_box.png")

    # ----------------------------------------------------------------------------------

    # ===============================
    # Plot 2: v^2 vs Energy (Gradient Plot)
    # ===============================
    if all_scan_data:
        # Group scan data by grid size for plotting
        data_by_grid = {}
        for v2, E, grid_size in all_scan_data:
            if grid_size not in data_by_grid:
                data_by_grid[grid_size] = {'v2': [], 'E': []}
            data_by_grid[grid_size]['v2'].append(v2)
            data_by_grid[grid_size]['E'].append(E)
        
        # --- FIX: Explicitly create Figure (fig) and Axes (ax) ---
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Plot each grid size scan using the assigned color, plotting ON the 'ax' object
        for i, grid_size in enumerate(grid_sizes):
            v2_data = data_by_grid[grid_size]['v2']
            E_data = data_by_grid[grid_size]['E']
            
            # Ensure the data is sorted by v2 for a smooth line plot
            sorted_indices = np.argsort(v2_data)
            v2_sorted = np.array(v2_data)[sorted_indices]
            E_sorted = np.array(E_data)[sorted_indices]
            
            ax.plot( # <--- Use ax.plot() instead of plt.plot()
                v2_sorted, 
                E_sorted, 
                marker='.', 
                linestyle='-', 
                color=colors[i], 
                label=f'Box Size: {grid_size*0.3} fm'
            )

        # Create a custom color bar proxy
        sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=grid_sizes[0], vmax=grid_sizes[-1]))
        sm.set_array([]) 
        
        # --- FIX: Pass the Axes object 'ax' to plt.colorbar() ---
        cbar = plt.colorbar(sm, ax=ax, label='Box Size (fm)')
        cbar.set_ticks(grid_sizes*0.3)

        ax.set_xlabel(r"$v^2$ parameter (constrained value)") # <--- Use ax.set_xlabel()
        ax.set_ylabel("Energy $E$ (MeV)") # <--- Use ax.set_ylabel()
        ax.set_title(f"Constrained Energy $E(v^2)$ Scan for ${nucleus}$ ({functional})") # <--- Use ax.set_title()
        ax.grid(True) # <--- Use ax.grid()
        fig.savefig("v2_vs_energy_gradient.png") # <--- Use fig.savefig()
        plt.show()

if __name__ == "__main__":
    main()