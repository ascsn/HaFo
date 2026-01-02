set terminal postscript eps enhanced color "TimesNewRoman,30"
set colorsequence classic
set output "trend-deocc-170Yb_3.eps"
set border 15 lw 1
set tmargin 0
set bmargin 0
set rmargin 0
set lmargin 0
set size 1.03,1.55
set style line 1 lt 1 lc rgb "#FF0000" lw 4 pt 7 ps 0.8
set style line 2 lt 1 lc rgb "#00AA00" lw 4 pt 7 ps 0.8
set style line 3 lt 1 lc rgb "#0000FF" lw 4 pt 7 ps 0.8
set style line 4 lt 1 lc rgb "#AA00AA" lw 4 pt 7 ps 0.8
set style line 5 lt 1 lc rgb "#009977" lw 4 pt 7 ps 0.8
set style line 6 lt 1 lc rgb "#000000" lw 4
set style line 7 lt 0 lc rgb "#FF0000" lw 3
set style line 8 lt 0 lc rgb "#000000" lw 3
set multiplot
#
set size 0.8,0.6
set origin 0.2,0.2
set xrange [1:15]
set xlabel "{/Symbol= D}_{shell}" offset 0
set yrange [-0.1:0.37]
set ylabel "total energy  [MeV]" offset 1
set form x
#set key at graph 0.4,0.66
set key top right
# IMPORTANT: Change the data file name here
p  "trend-deocc-170Yb_3_data.txt"  i 3 u 2:($8+1374.530) t "{/Symbol e}_{cut} = 20 MeV" w l ls 1,\
   "trend-deocc-170Yb_3_data.txt"  i 0 u 2:($8+1374.5608) t "{/Symbol e}_{cut} = 40 MeV" w l ls 2,\
   "trend-deocc-170Yb_3_data.txt"  i 1 u 2:($8+1374.5528) t "{/Symbol e}_{cut} = 60 MeV" w l ls 3
#
set size 0.8,0.6
set origin 0.2,0.8
set yrange [-2.95:0.5]
set ylabel "charge radius  [10^{-3}fm]" offset 0.5
set xlabel "" offset 0
set form x ""
set key bottom right
set title "^{170}Yb, SV-bas, HFB" offset 0,-0.3
# IMPORTANT: Change the data file name here
p  "trend-deocc-170Yb_3_data.txt"  i 3 u 2:(($13-5.240)*1e3) t "{/Symbol e}_{cut} = 20 MeV" w l ls 1,\
   "trend-deocc-170Yb_3_data.txt"  i 0 u 2:(($13-5.2914)*1e3) t "{/Symbol e}_{cut} = 40 MeV" w l ls 2,\
   "trend-deocc-170Yb_3_data.txt"  i 1 u 2:(($13-5.2918)*1e3) t "{/Symbol e}_{cut} = 60 MeV" w l ls 3
#
!epstopdf trend-deocc-170Yb_3.eps
!rm -f trend-deocc-170Yb_3.eps
!evince trend-deocc-170Yb_3.pdf
#!open -a preview trend-deocc-170Yb_3.pdf
q