set terminal postscript eps enhanced color "TimesNewRoman,30"
set colorsequence classic
set output "varycutp_BCS.eps"
set border 15 lw 1
set tmargin 0
set bmargin 0
set rmargin 0
set lmargin 0
set size 2.43,3.43
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
set xrange [14:101]
#set yrange [0:1]
#set ytics 0.001
set xlabel "cutoff energy  [MeV]"
set ylabel "E_B  [MeV]" offset 0
#set yrange [0.0029:-0.0019]
unset key 
#set label "^{120}Sn, UNEDF1, BCS" at graph 0.33,0.09
p  "collect_146Gd/values_146Gd_0.0d0_bcs.dat" u 1:($3) t "width=0.0" w l ls 1,\
  "collect_146Gd/values_146Gd_0.5d0_bcs.dat" u 1:($3) t "width=0.5" w l ls 2,\
  "collect_146Gd/values_146Gd_1.0d0_bcs.dat" u 1:($3) t "width=1.0" w l ls 3,\
  "collect_146Gd/values_146Gd_1.5d0_bcs.dat" u 1:($3) t "width=1.5" w l ls 4,\
  "collect_146Gd/values_146Gd_2.0d0_bcs.dat" u 1:($3) t "width=2.0" w l ls 5,\
  "collect_146Gd/values_146Gd_3.0d0_bcs.dat" u 1:($3) t "width=3.0" w l ls 6,\
  "collect_146Gd/values_146Gd_4.0d0_bcs.dat" u 1:($3) t "width=4.0" w l ls 7,\
  "collect_146Gd/values_146Gd_5.0d0_bcs.dat" u 1:($3) t "width=5.0" w l ls 8
#
set size 0.8,0.8
set origin 0.2,1.0
set xrange [14:101]
set form x ""
#set yrange [0:1]
#set ytics 0.001
set xlabel ""
set ylabel "r_C  [fm]" offset 0
#set yrange [0.0029:-0.0019]
unset key 
#set label "^{120}Sn, UNEDF1, BCS" at graph 0.33,0.09
p  "collect_146Gd/values_146Gd_0.0d0_bcs.dat" u 1:($2) t "width=0.0" w l ls 1,\
  "collect_146Gd/values_146Gd_0.5d0_bcs.dat" u 1:($2) t "width=0.5" w l ls 2,\
  "collect_146Gd/values_146Gd_1.0d0_bcs.dat" u 1:($2) t "width=1.0" w l ls 3,\
  "collect_146Gd/values_146Gd_1.5d0_bcs.dat" u 1:($2) t "width=1.5" w l ls 4,\
  "collect_146Gd/values_146Gd_2.0d0_bcs.dat" u 1:($2) t "width=2.0" w l ls 5,\
  "collect_146Gd/values_146Gd_3.0d0_bcs.dat" u 1:($2) t "width=3.0" w l ls 6,\
  "collect_146Gd/values_146Gd_4.0d0_bcs.dat" u 1:($2) t "width=4.0" w l ls 7,\
  "collect_146Gd/values_146Gd_5.0d0_bcs.dat" u 1:($2) t "width=5.0" w l ls 8
#
set size 0.8,0.8
set origin 0.2,1.8
set xrange [14:101]
set form x ""
#set yrange [0:1]
#set ytics 0.001
set xlabel ""
set ylabel "E_{pair}  [MeV]" offset 0
#set yrange [0.0029:-0.0019]
unset key 
#set label "^{120}Sn, UNEDF1, BCS" at graph 0.33,0.09
p  "collect_146Gd/values_146Gd_0.0d0_bcs.dat" u 1:($4) t "width=0.0" w l ls 1,\
  "collect_146Gd/values_146Gd_0.5d0_bcs.dat" u 1:($4) t "width=0.5" w l ls 2,\
  "collect_146Gd/values_146Gd_1.0d0_bcs.dat" u 1:($4) t "width=1.0" w l ls 3,\
  "collect_146Gd/values_146Gd_1.5d0_bcs.dat" u 1:($4) t "width=1.5" w l ls 4,\
  "collect_146Gd/values_146Gd_2.0d0_bcs.dat" u 1:($4) t "width=2.0" w l ls 5,\
  "collect_146Gd/values_146Gd_3.0d0_bcs.dat" u 1:($4) t "width=3.0" w l ls 6,\
  "collect_146Gd/values_146Gd_4.0d0_bcs.dat" u 1:($4) t "width=4.0" w l ls 7,\
  "collect_146Gd/values_146Gd_5.0d0_bcs.dat" u 1:($4) t "width=5.0" w l ls 8
#
set size 0.8,0.8
set origin 0.2,2.6
set xrange [14:101]
set form x ""
#set yrange [0:1]
#set ytics 0.001
set xlabel ""
set ylabel "gap  [MeV]" offset 0
#set yrange [0.0029:-0.0019]
unset key 
#set label "^{120}Sn, UNEDF1, BCS" at graph 0.33,0.09
p  "collect_146Gd/values_146Gd_0.0d0_bcs.dat" u 1:($5) t "width=0.0" w l ls 1,\
  "collect_146Gd/values_146Gd_0.5d0_bcs.dat" u 1:($5) t "width=0.5" w l ls 2,\
  "collect_146Gd/values_146Gd_1.0d0_bcs.dat" u 1:($5) t "width=1.0" w l ls 3,\
  "collect_146Gd/values_146Gd_1.5d0_bcs.dat" u 1:($5) t "width=1.5" w l ls 4,\
  "collect_146Gd/values_146Gd_2.0d0_bcs.dat" u 1:($5) t "width=2.0" w l ls 5,\
  "collect_146Gd/values_146Gd_3.0d0_bcs.dat" u 1:($5) t "width=3.0" w l ls 6,\
  "collect_146Gd/values_146Gd_4.0d0_bcs.dat" u 1:($5) t "width=4.0" w l ls 7,\
  "collect_146Gd/values_146Gd_5.0d0_bcs.dat" u 1:($5) t "width=5.0" w l ls 8
#
set size 0.8,0.8
set origin 1.2,0.2
set xrange [14:101]
set form x
set yrange [0:1.1]
#set ytics 0.001
set xlabel "cutoff energy  [MeV]"
set ylabel "" offset 0
#set yrange [0.0029:-0.0019]
set key at graph 1.5,0.9
#set label "^{120}Sn, UNEDF1, BCS" at graph 0.33,0.9
p  "collect_146Gd/values_146Gd_0.0d0_bcs.dat" u 1:($7) t "width=0.0" w l ls 1,\
  "collect_146Gd/values_146Gd_0.5d0_bcs.dat" u 1:($7) t "width=0.5" w l ls 2,\
  "collect_146Gd/values_146Gd_1.0d0_bcs.dat" u 1:($7) t "width=1.0" w l ls 3,\
  "collect_146Gd/values_146Gd_1.5d0_bcs.dat" u 1:($7) t "width=1.5" w l ls 4,\
  "collect_146Gd/values_146Gd_2.0d0_bcs.dat" u 1:($7) t "width=2.0" w l ls 5,\
  "collect_146Gd/values_146Gd_3.0d0_bcs.dat" u 1:($7) t "width=3.0" w l ls 6,\
  "collect_146Gd/values_146Gd_4.0d0_bcs.dat" u 1:($7) t "width=4.0" w l ls 7,\
  "collect_146Gd/values_146Gd_5.0d0_bcs.dat" u 1:($7) t "width=5.0" w l ls 8
#
set size 0.8,0.8
set origin 1.2,1.0
set xrange [14:101]
set form x ""
set yrange [-0.09:0]
#set ytics 0.001
set xlabel ""
set ylabel "" offset 0
#set yrange [0.0029:-0.0019]
unset key 
#set label "^{120}Sn, UNEDF1, BCS" at graph 0.33,0.09
p  "collect_146Gd/values_146Gd_0.0d0_bcs.dat" u 1:($6) t "width=0.0" w l ls 1,\
  "collect_146Gd/values_146Gd_0.5d0_bcs.dat" u 1:($6) t "width=0.5" w l ls 2,\
  "collect_146Gd/values_146Gd_1.0d0_bcs.dat" u 1:($6) t "width=1.0" w l ls 3,\
  "collect_146Gd/values_146Gd_1.5d0_bcs.dat" u 1:($6) t "width=1.5" w l ls 4,\
  "collect_146Gd/values_146Gd_2.0d0_bcs.dat" u 1:($6) t "width=2.0" w l ls 5,\
  "collect_146Gd/values_146Gd_3.0d0_bcs.dat" u 1:($6) t "width=3.0" w l ls 6,\
  "collect_146Gd/values_146Gd_4.0d0_bcs.dat" u 1:($6) t "width=4.0" w l ls 7,\
  "collect_146Gd/values_146Gd_5.0d0_bcs.dat" u 1:($6) t "width=5.0" w l ls 8
#
set size 0.8,0.8
set origin 1.2,1.8
set xrange [14:101]
set form x ""
set yrange [0:2]
#set ytics 0.001
set xlabel ""
set ylabel "" offset 0
#set yrange [0.0029:-0.0019]
unset key
#set label "^{120}Sn, UNEDF1, BCS" at graph 0.33,0.09
p  "collect_146Gd/values_146Gd_0.0d0_bcs.dat" u 1:($8) t "width=0.0" w l ls 1,\
  "collect_146Gd/values_146Gd_0.5d0_bcs.dat" u 1:($8) t "width=0.5" w l ls 2,\
  "collect_146Gd/values_146Gd_1.0d0_bcs.dat" u 1:($8) t "width=1.0" w l ls 3,\
  "collect_146Gd/values_146Gd_1.5d0_bcs.dat" u 1:($8) t "width=1.5" w l ls 4,\
  "collect_146Gd/values_146Gd_2.0d0_bcs.dat" u 1:($8) t "width=2.0" w l ls 5,\
  "collect_146Gd/values_146Gd_3.0d0_bcs.dat" u 1:($8) t "width=3.0" w l ls 6,\
  "collect_146Gd/values_146Gd_4.0d0_bcs.dat" u 1:($8) t "width=4.0" w l ls 7,\
  "collect_146Gd/values_146Gd_5.0d0_bcs.dat" u 1:($8) t "width=5.0" w l ls 8
#
set size 0.8,0.8
set origin 1.2,2.6
set xrange [14:101]
set form x ""
set yrange [-0.49:0.5]
#set ytics 0.001
set xlabel ""
set ylabel "" offset 0
#set yrange [0.0029:-0.0019]
unset key 
set label "^{120}Sn" at graph 1.1,0.8
set label "UNEDF1" at graph 1.1,0.7
set label "BCS" at graph 1.1,0.6
p  "collect_146Gd/values_146Gd_0.0d0_bcs.dat" u 1:($9) t "width=0.0" w l ls 1,\
  "collect_146Gd/values_146Gd_0.5d0_bcs.dat" u 1:($9) t "width=0.5" w l ls 2,\
  "collect_146Gd/values_146Gd_1.0d0_bcs.dat" u 1:($9) t "width=1.0" w l ls 3,\
  "collect_146Gd/values_146Gd_1.5d0_bcs.dat" u 1:($9) t "width=1.5" w l ls 4,\
  "collect_146Gd/values_146Gd_2.0d0_bcs.dat" u 1:($9) t "width=2.0" w l ls 5,\
  "collect_146Gd/values_146Gd_3.0d0_bcs.dat" u 1:($9) t "width=3.0" w l ls 6,\
  "collect_146Gd/values_146Gd_4.0d0_bcs.dat" u 1:($9) t "width=4.0" w l ls 7,\
  "collect_146Gd/values_146Gd_5.0d0_bcs.dat" u 1:($9) t "width=5.0" w l ls 8
#
!eps2pdf -f varycutp_BCS.eps
!rm -f varycutp_BCS.eps
!evince varycutp_BCS.pdf
q
