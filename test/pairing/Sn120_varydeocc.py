import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib import gridspec
import scienceplots

# Apply the 'science' style for the plot design
plt.style.use('science')

# --- Set global font sizes for all text elements ---
mpl.rcParams['axes.labelsize'] = 25    # font size of the x and y labels
mpl.rcParams['xtick.labelsize'] = 28   # font size of the x tick labels
mpl.rcParams['ytick.labelsize'] = 28   # font size of the y tick labels
mpl.rcParams['legend.fontsize'] = 20   # legend font size
mpl.rcParams['axes.titlesize'] = 22    # font size of the axes title
mpl.rcParams['lines.linewidth'] = 4.0

# Data provided in the prompt, separated by gnuplot block markers.
data_str = """
#0:collect_cat_FyDr_120Sn_cut15_wid0_0.0d0_bcs.dat
 0.0  2   322  50  70  4.659498   -1019.91       -0.573 -19.218    0.382   1.739   1.499  11.051
 0.0  3   464  50  70  4.659498   -1019.91       -0.573 -19.219    0.382   1.738   1.499  11.055
 0.0  4   636  50  70  4.659498   -1019.91       -0.573 -19.219    0.382   1.738   1.499  11.055
 0.0  5   866  50  70  4.659498   -1019.91       -0.573 -19.219    0.382   1.738   1.499  11.055
 0.0  6  1122  50  70  4.659498   -1019.91       -0.573 -19.219    0.382   1.738   1.499  11.055
 0.0  7  1420  50  70  4.659498   -1019.91       -0.573 -19.219    0.382   1.738   1.499  11.055
 0.0  8  1792  50  70  4.659498   -1019.91       -0.573 -19.219    0.382   1.738   1.499  11.055
 0.0  9  2200  50  70  4.659498   -1019.91       -0.573 -19.219    0.382   1.738   1.499  11.055
 0.0  10 2658   50  70  4.659498   -1019.91       -0.573 -19.219    0.382   1.738   1.499  11.055
#1:collect_cat_FyDr_120Sn_cut15_wid0_0.0d0_hfb.dat
 0.0  2   322  50  70  4.668905   -1021.05       -0.662 -26.828    0.422   2.351   1.567  11.412
 0.0  3   464  50  70  4.668917   -1021.05       -0.668 -26.830    0.421   2.351   1.588  11.413
 0.0  4   636  50  70  4.668905   -1021.05       -0.662 -26.828    0.422   2.351   1.567  11.413
 0.0  5   866  50  70  4.668917   -1021.05       -0.668 -26.830    0.421   2.351   1.588  11.413
 0.0  6  1122  50  70  4.668917   -1021.05       -0.668 -26.830    0.421   2.351   1.588  11.413
 0.0  7  1420  50  70  4.668917   -1021.05       -0.668 -26.830    0.421   2.351   1.588  11.413
 0.0  8  1792  50  70  4.668906   -1021.05       -0.662 -26.828    0.422   2.351   1.567  11.413
 0.0  9  2200  50  70  4.668906   -1021.05       -0.662 -26.828    0.422   2.351   1.567  11.413
 0.0  10 2658   50  70  4.668917   -1021.05       -0.668 -26.830    0.421   2.351   1.588  11.413
#2:collect_cat_FyDr_120Sn_cut15_wid1_1.0d0_bcs.dat
 1.0  2   322  50  70  4.652495   -1020.38       -0.939 -18.124    0.464   1.665   2.024  10.882
 1.0  3   464  50  70  4.652495   -1020.38       -0.939 -18.124    0.464   1.665   2.024  10.886
 1.0  4   636  50  70  4.652495   -1020.38       -0.939 -18.124    0.464   1.665   2.024  10.886
 1.0  5   866  50  70  4.652495   -1020.38       -0.939 -18.124    0.464   1.665   2.024  10.886
 1.0  6  1122  50  70  4.652495   -1020.38       -0.939 -18.124    0.464   1.665   2.024  10.886
 1.0  7  1420  50  70  4.652495   -1020.38       -0.939 -18.124    0.464   1.665   2.024  10.886
 1.0  8  1792  50  70  4.652495   -1020.38       -0.939 -18.124    0.464   1.665   2.024  10.886
 1.0  9  2200  50  70  4.652495   -1020.38       -0.939 -18.124    0.464   1.665   2.024  10.886
 1.0  10 2658   50  70  4.652495   -1020.38       -0.939 -18.124    0.464   1.665   2.024  10.886
#3:collect_cat_FyDr_120Sn_cut15_wid1_1.0d0_hfb.dat
 1.0  2   322  50  70  4.661447   -1020.45       -1.132 -19.977    0.502   1.831   2.253  10.910
 1.0  3   464  50  70  4.661447   -1020.45       -1.132 -19.977    0.502   1.831   2.253  10.910
 1.0  4   636  50  70  4.661508   -1020.46       -1.133 -20.038    0.503   1.830   2.255  10.952
 1.0  5   866  50  70  4.661508   -1020.46       -1.133 -20.038    0.503   1.830   2.255  10.952
 1.0  6  1122  50  70  4.661508   -1020.46       -1.133 -20.038    0.503   1.830   2.255  10.952
 1.0  7  1420  50  70  4.661445   -1020.45       -1.131 -19.971    0.502   1.830   2.251  10.910
 1.0  8  1792  50  70  4.661444   -1020.45       -1.131 -19.974    0.502   1.831   2.252  10.911
 1.0  9  2200  50  70  4.661444   -1020.45       -1.131 -19.974    0.502   1.831   2.252  10.911
 1.0  10 2658   50  70  4.661443   -1020.45       -1.131 -19.973    0.502   1.831   2.252  10.911
#4:collect_cat_FyDr_120Sn_cut45_wid0_0.0d0_bcs.dat
 0.0  2   322  50  70  4.655275   -1019.28       -0.670 -18.709    0.361   1.559   1.855  12.003
 0.0  3   464  50  70  4.658941   -1019.76       -0.734 -21.827    0.357   1.623   2.059  13.445
 0.0  4   636  50  70  4.665090   -1020.57       -0.848 -27.538    0.362   1.765   2.342  15.598
 0.0  5   866  50  70  4.672698   -1021.57       -1.028 -35.353    0.378   1.943   2.718  18.194
 0.0  6  1122  50  70  4.675092   -1021.82       -1.082 -37.682    0.383   1.969   2.828  19.134
 0.0  7  1420  50  70  4.675425   -1021.86       -1.089 -37.993    0.383   1.964   2.845  19.342
 0.0  8  1792  50  70  4.675436   -1021.86       -1.089 -38.004    0.383   1.962   2.846  19.370
 0.0  9  2200  50  70  4.675436   -1021.86       -1.089 -38.004    0.383   1.962   2.846  19.373
 0.0  10 2658   50  70  4.675436   -1021.86       -1.089 -38.004    0.383   1.962   2.846  19.374
#5:collect_cat_FyDr_120Sn_cut45_wid0_0.0d0_hfb.dat
 0.0  2   322  50  70  4.669001   -1019.81       -0.723 -26.489    0.414   2.085   1.745  12.702
 0.0  3   464  50  70  4.669291   -1019.89       -0.744 -27.092    0.413   2.098   1.803  12.914
 0.0  4   636  50  70  4.669266   -1019.89       -0.727 -27.087    0.411   2.097   1.768  12.916
 0.0  5   866  50  70  4.669246   -1019.87       -0.758 -26.893    0.414   2.090   1.830  12.870
 0.0  6  1122  50  70  4.669249   -1019.86       -0.758 -26.861    0.414   2.088   1.832  12.862
 0.0  7  1420  50  70  4.669317   -1019.89       -0.765 -27.100    0.414   2.097   1.847  12.921
 0.0  8  1792  50  70  4.669242   -1019.88       -0.758 -26.972    0.413   2.092   1.833  12.894
 0.0  9  2200  50  70  4.668705   -1019.77       -0.748 -26.187    0.411   2.057   1.820  12.729
 0.0  10 2658   50  70  4.669307   -1019.89       -0.759 -27.097    0.414   2.097   1.835  12.922
#6:collect_cat_FyDr_120Sn_cut45_wid1_1.0d0_bcs.dat
 1.0  2   322  50  70  4.657260   -1019.59       -0.967 -19.381    0.416   1.603   2.322  12.093
 1.0  3   464  50  70  4.658725   -1019.86       -1.094 -21.074    0.417   1.620   2.620  13.006
 1.0  4   636  50  70  4.660283   -1020.11       -1.252 -22.757    0.425   1.641   2.943  13.872
 1.0  5   866  50  70  4.661357   -1020.29       -1.362 -23.943    0.431   1.648   3.163  14.529
 1.0  6  1122  50  70  4.661538   -1020.32       -1.379 -24.144    0.430   1.641   3.210  14.712
 1.0  7  1420  50  70  4.661550   -1020.32       -1.380 -24.157    0.429   1.638   3.216  14.744
 1.0  8  1792  50  70  4.661550   -1020.32       -1.380 -24.157    0.429   1.638   3.217  14.749
 1.0  9  2200  50  70  4.661550   -1020.32       -1.380 -24.157    0.429   1.638   3.217  14.749
 1.0  10 2658   50  70  4.661550   -1020.32       -1.380 -24.157    0.429   1.638   3.217  14.749
#7:collect_cat_FyDr_120Sn_cut45_wid1_1.0d0_hfb.dat
 1.0  2   322  50  70  4.659234   -1019.66       -1.344 -22.410    0.485   1.749   2.772  12.810
 1.0  3   464  50  70  4.661668   -1019.92       -1.573 -24.463    0.509   1.787   3.094  13.687
 1.0  4   636  50  70  4.662140   -1019.97       -1.611 -24.882    0.511   1.783   3.153  13.956
 1.0  5   866  50  70  4.662180   -1019.97       -1.614 -24.920    0.511   1.779   3.159  14.007
 1.0  6  1122  50  70  4.662184   -1019.97       -1.616 -24.921    0.511   1.778   3.165  14.013
 1.0  7  1420  50  70  4.662185   -1019.97       -1.616 -24.921    0.510   1.778   3.167  14.016
 1.0  8  1792  50  70  4.662186   -1019.97       -1.616 -24.922    0.510   1.778   3.169  14.020
 1.0  9  2200  50  70  4.662187   -1019.97       -1.617 -24.922    0.510   1.777   3.170  14.022
 1.0  10 2658   50  70  4.662187   -1019.97       -1.617 -24.922    0.510   1.777   3.170  14.023
"""

# Use a context manager for plotting to ensure resources are handled correctly
with plt.style.context(['science']):
    # Manually split the data into 8 blocks based on the gnuplot 'index'
    data_sections = []
    block_indices = [i for i, line in enumerate(data_str.splitlines()) if line.strip().startswith('#')]
    
    for i in range(len(block_indices)):
        start = block_indices[i] + 1
        end = block_indices[i+1] if i+1 < len(block_indices) else len(data_str.splitlines())
        section_lines = [line for line in data_str.splitlines()[start:end] if line.strip()]
        if section_lines:
            data_sections.append(np.genfromtxt(section_lines))

    # Use `GridSpec` to define the subplot layout and control the spacing.
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig)

    # The gnuplot script uses columns 3, 6, and 9. Python uses 0-based indexing.
    col_x = 2
    col_y1 = 8
    col_y2 = 5

    # Define common label text
    label_x = r'$N_\mathrm{st}^{1/3}$'
    label_y2 = r'neutron pairing energy (MeV)'
    label_y1 = r'charge radius ($10^{-3}$ fm)'
    label1=r'$\tilde{R}_\mathcal{F}$=0, HFB'
    label2=r'$\tilde{R}_\mathcal{F}$=0, BCS'
    label3=r'$\tilde{R}_\mathcal{F}$=1, HFB'
    label4=r'$\tilde{R}_\mathcal{F}$=1, BCS'

    # --- Top-Left Plot (Cut 15 MeV, Pairing Energy) ---
    ax_tl = fig.add_subplot(gs[0, 0])
    ax_tl.set_ylabel(label_y1)
    ax_tl.set_xlim(6.8, 13.95)
    ax_tl.set_ylim(-1, 19.9)

    # ax_tl.set_title(r'Fy($\Delta r$), $^{120}$Sn, $\varepsilon_\mathrm{cut}$ = 15 MeV', loc='center')

    # Plotting the data with the specified transformations for each block
    ax_tl.plot(data_sections[1][:, col_x]**(1/3), (data_sections[1][:, col_y2] - 4.668917) * 1e3, color='red', linewidth=3.0)
    ax_tl.plot(data_sections[0][:, col_x]**(1/3), (data_sections[0][:, col_y2] - 4.659498) * 1e3, color='#00aa00', linewidth=3.0)
    ax_tl.plot(data_sections[3][:, col_x]**(1/3), (data_sections[3][:, col_y2] - 4.661443) * 1e3, color='blue', linewidth=3.0)
    ax_tl.plot(data_sections[2][:, col_x]**(1/3), (data_sections[2][:, col_y2] - 4.652495) * 1e3, color='purple', linewidth=3.0)
    ax_tl.text(0.35, 0.95, r'$\varepsilon_\mathrm{cut}$ = 15 MeV', transform=ax_tl.transAxes, fontsize=30, verticalalignment='top')

    # --- Bottom-Left Plot (Cut 15 MeV, Charge Radius) ---
    ax_bl = fig.add_subplot(gs[1, 0])
    ax_bl.set_ylabel(label_y2)
    ax_bl.set_xlabel(label_x)
    ax_bl.set_xlim(6.8, 13.95)
    ax_bl.set_ylim(-1, 19.9)


    ax_bl.plot(data_sections[1][:, col_x]**(1/3), data_sections[1][:, col_y1] + 26.830, label=label1, color = 'red', linewidth = 3.0)
    ax_bl.plot(data_sections[0][:, col_x]**(1/3), data_sections[0][:, col_y1] + 19.219, label=label2, color='#00aa00', linewidth = 3.0)
    ax_bl.plot(data_sections[3][:, col_x]**(1/3), data_sections[3][:, col_y1] + 19.973, label=label3, color = 'blue', linewidth = 3.0)
    ax_bl.plot(data_sections[2][:, col_x]**(1/3), data_sections[2][:, col_y1] + 18.124, label=label4, color = 'purple', linewidth = 3.0)


    # --- Top-Right Plot (Cut 45 MeV, Pairing Energy) ---
    ax_tr = fig.add_subplot(gs[0, 1])
    ax_tr.set_xlim(6.8, 13.95)
    ax_tr.set_ylim(-22, 2)
    ax_tr.text(0.35, 0.15, r'$\varepsilon_\mathrm{cut}$ = 45 MeV', transform=ax_tr.transAxes, fontsize=30, verticalalignment='top')

    # ax_tr.set_title(r'Fy($\Delta r$), $^{120}$Sn, $\varepsilon_\mathrm{cut}$ = 45 MeV', loc='center')

    ax_tr.plot(data_sections[5][:, col_x]**(1/3), (data_sections[5][:, col_y2] - 4.669307) * 1e3, label=label1, color='red', linewidth=3.0)
    ax_tr.plot(data_sections[4][:, col_x]**(1/3), (data_sections[4][:, col_y2] - 4.675436) * 1e3, label=label2, color='#00aa00', linewidth=3.0)
    ax_tr.plot(data_sections[7][:, col_x]**(1/3), (data_sections[7][:, col_y2] - 4.662187) * 1e3, label=label3, color='blue', linewidth=3.0)
    ax_tr.plot(data_sections[6][:, col_x]**(1/3), (data_sections[6][:, col_y2] - 4.661550) * 1e3, label=label4, color='purple', linewidth=3.0)

    # --- Bottom-Right Plot (Cut 45 MeV, Charge Radius) ---
    ax_br = fig.add_subplot(gs[1, 1])
    ax_br.set_xlabel(label_x)
    ax_br.set_xlim(6.8, 13.95)
    ax_br.set_ylim(-1, 19.9)

    ax_br.plot(data_sections[5][:, col_x]**(1/3), data_sections[5][:, col_y1] + 27.097, label=label1, color='red', linewidth=3.0)
    ax_br.plot(data_sections[4][:, col_x]**(1/3), data_sections[4][:, col_y1] + 38.004, label=label2, color='#00aa00', linewidth=3.0)
    ax_br.plot(data_sections[7][:, col_x]**(1/3), data_sections[7][:, col_y1] + 24.922, label=label3, color='blue', linewidth=3.0)
    ax_br.plot(data_sections[6][:, col_x]**(1/3), data_sections[6][:, col_y1] + 24.157, label=label4, color='purple', linewidth=3.0)
    ax_br.legend()

    # --- Final Layout Adjustments ---
    # Align y-labels for the left-side subplots
    fig.align_ylabels([ax_tl, ax_bl])

    # Increase the tick length for all axes
    for ax in [ax_tl, ax_bl, ax_tr, ax_br]:
        ax.tick_params(axis='both', which='major', length=8)
        ax.tick_params(axis='y', which='minor', length=4)
        ax.tick_params(axis='x', which='minor', length=0)
    
    # Set x-ticks to be separated by 1 for the bottom panels
    ax_bl.set_xticks(np.arange(7, 14, 1))
    ax_br.set_xticks(np.arange(7, 14, 1))
    ax_tl.set_xticks(np.arange(7, 14, 1))
    ax_tr.set_xticks(np.arange(7, 14, 1))
    
    # Remove x-tick labels from the upper panels
    ax_tl.set_xticklabels([])
    ax_tr.set_xticklabels([])
    
    # Remove y-tick labels and y-labels from the right panels
    ax_tr.set_yticklabels([])
    ax_br.set_yticklabels([])
    ax_tr.set_ylabel('') # Explicitly remove y-label from right panel
    ax_br.set_ylabel('') # Explicitly remove y-label from right panel

    # Set the horizontal and vertical spacing to zero to remove gaps between plots.
    fig.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.1)
    plt.subplots_adjust(wspace=0, hspace=0)
    
    # Save the figure to a file and show it
    plt.savefig('Sn120-varydeocc_corrected_cut15_python.pdf', dpi=600, bbox_inches='tight')
    plt.show()