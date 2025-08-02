set terminal postscript eps enhanced color "TimesNewRoman,30"
set colorsequence classic
set output "varycutn_BCS.eps"
set border 15 lw 1
set tmargin 0
set bmargin 0
set rmargin 0
set lmargin 0
set size 1.43,1.03
set style line 1 lt 1 lc rgb "#FF0000" lw 2 pt 7 ps 1.5
set style line 2 lt 1 lc rgb "#00AA00" lw 2 pt 5 ps 1.5
set style line 3 lt 1 lc rgb "#0000FF" lw 2 pt 4 ps 1.5
set style line 4 lt 1 lc rgb "#AA00AA" lw 2 pt 6 ps 1.5
set style line 5 lt 1 lc rgb "#00CCAA" lw 2 pt 9 ps 1.5
set style line 6 lt 1 lc rgb "#000000" lw 2
set style line 7 lt 0 lc rgb "#FF0000" lw 3
set style line 8 lt 0 lc rgb "#000000" lw 3
set multiplot
#
set size 0.8,0.8
set origin 0.2,0.2
set xrange [14:121]
set yrange [0:1]
#set ytics 0.001
set xlabel "cutoff energy  [MeV]"
set ylabel "(E_{pair}(cut)-E_{pair}(15)/E_{pair}(120)" offset 0
#set yrange [0.0029:-0.0019]
set key at graph 1.5,0.9
set label "^{120}Sn, UNEDF1, BCS" at graph 0.33,0.09
p  "varycut_BCS.gnu" i 0 u 1:(-$7/(9.6000)) t "width=0.0" w l ls 1,\
  "varycut_BCS.gnu" i 1 u 1:(-$7/(4.2700)) t "width=0.5" w l ls 2,\
  "varycut_BCS.gnu" i 2 u 1:(-$7/(1.3500)) t "width=1.0" w l ls 3,\
  "varycut_BCS.gnu" i 3 u 1:(-$7/(0.73004)) t "width=1.5" w l ls 4,\
  "varycut_BCS.gnu" i 4 u 1:(-$7/(0.73999)) t "width=1.6" w l ls 5,\
  "varycut_BCS.gnu" i 6 u 1:(-$7/(0.97998)) t "width=1.8" w l ls 6,\
  "varycut_BCS.gnu" i 8 u 1:(-$7/(1.8900)) t "width=2.0" w l ls 7
#
!eps2pdf -f varycutn_BCS.eps
!rm -f varycutn_BCS.eps
!evince varycutn_BCS.pdf
q
