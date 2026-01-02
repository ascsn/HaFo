#set terminal postscript eps enhanced color "TimesNewRoman,34"
set terminal epslatex color colortext font 17 solid round standalone
set colorsequence classic
#set output "varyocc_widths_3color.eps"
set output "varyocc_widths_3color.tex"
set border 15 lw 1
set tmargin 0
set bmargin 0
set rmargin 0
set lmargin 0
set size 2.23,2.88
set origin 0.0,0.0
set style line 1 lt 1 lc rgb "#FF0000" lw 5 pt 7 ps 1.5
set style line 2 lt 1 lc rgb "#00AA00" lw 5 pt 7 ps 1.5
set style line 3 lt 1 lc rgb "#0000FF" lw 5 pt 7 ps 1.5
set style line 4 lt 1 lc rgb "#AA00AA" lw 5 pt 7 ps 1.5
set style line 5 lt 1 lc rgb "#009977" lw 5 pt 7 ps 1.5
set style line 6 lt 1 lc rgb "#000000" lw 2
set style line 7 lt 0 lc rgb "#FF0000" lw 5
set style line 8 lt 0 lc rgb "#000000" lw 5
set palette defined (0 "black", 0.2 "blue", 0.4 "green", 0.8 "orange", 1.0 "red")
unset colorbox
set cbrange [0:20]
set multiplot
#
unset arrow
set size 0.8,0.6
set origin 0.2,2.2
set xrange [-47:89]
set yrange [0:18]
set ylabel "s.p. radius (fm)" offset 1
set xlabel "" offset 0,0.2
set xtics 20
set ytics 5
unset key
#set key bottom right 
set arrow 1 from 0,graph 0 to 0,graph 1 nohead lt 0 lw 2
set arrow 2 from -47,4.643 to 0,4.643 nohead lt 0 lw 2
set label 1 "BCS, $\\tilde{R}_\\mathcal{F}$ = 0 fm" at graph 0.55,0.10
set label 2 "(a)" at graph 0.05,0.9
#set label "$^{120}$Sn, SV-bas, BCS" at graph 0.03,0.0
p  "varyocc_widths_3.gnu" i 0 u 3:6:2 t "$\\tilde{R}_\\mathcal{F}$=0" w p ls 1 palette
#
unset arrow
set size 0.8,0.6
set origin 0.2,0.2
set xrange [-47:89]
set yrange [0:18]
set ylabel "" offset 1
set xtics 20
set ytics 5
set xlabel "s.p. energy  (MeV)" offset 0,0.2
#set label "$^{120}$Sn, SV-bas, HFB" at graph 0.03,0.09
set arrow 1 from 0,graph 0 to 0,graph 1 nohead lt 0 lw 2
set arrow 2 from -47,4.643 to 0,4.643 nohead lt 0 lw 2
set label 1 "HFB, $\\tilde{R}_\\mathcal{F}$ = 0 fm" at graph 0.55,0.10
set label 2 "(g)" at graph 0.05,0.9
p  "varyocc_widths_3.gnu" i 1 u 3:6:2 t "$\\tilde{R}_\\mathcal{F}$=0" w p ls 1 palette
#
unset arrow
set size 0.8,0.6
set origin 0.2,0.8
set xrange [-47:89]
set form x ""
set yrange [0:18]
set ylabel "s.p. radius (fm)" offset 1
set ytics 5
set xtics 20
set xlabel "" offset 0
#set label "$^{120}$Sn, SV-bas, HFB" at graph 0.03,0.09
set arrow 1 from 0,graph 0 to 0,graph 1 nohead lt 0 lw 2
set arrow 2 from -47,4.643 to 0,4.643 nohead lt 0 lw 2
set label 1 "HFB, $\\tilde{R}_\\mathcal{F}$ = 1 fm" at graph 0.55,0.10
set label 2 "(e)" at graph 0.05,0.9
#set log cb
#set form cb "$10^{%T}$"
p  "varyocc_widths_3.gnu" i 2  u 3:6:2 t "$\\tilde{R}_\\mathcal{F}$=1" w p ls 3 palette
#p  "varyocc_widths_3.gnu" i 2  u ($4>3e-8 ? $3:1/0):6:2 t "$\\tilde{R}_\\mathcal{F}$=1" w p ls 3  palette
#
unset arrow
unset label
set size 0.8,0.6
set origin 0.2,1.4
set xrange [-47:89]
set yrange [0:18]
set ylabel "" offset 1
set xlabel "" offset 0
set xtics 10
set ytics 5
set arrow 1 from 0,graph 0 to 0,graph 1 nohead lt 0 lw 2
set arrow 2 from -47,4.636 to 0,4.636 nohead lt 0 lw 2
set label 1 "HFB, $\\tilde{R}_\\mathcal{F}$ = 2 fm" at graph 0.55,0.10
set label 2 "(c)" at graph 0.05,0.9
p  "varyocc_widths_3.gnu" i 3 u 3:6:2 t "$\\tilde{R}_\\mathcal{F}$=2" w p ls 5 palette
#
unset arrow
set size 0.8,0.6
set form x
set origin 1.2,2.2
set xrange [-47:89]
set yrange [0:4.45]
set ytics 1
set xtics 20
set ylabel "pairing gap $\\Delta_n$  (MeV)" offset 0
set xlabel "" offset 0
unset key
set arrow 1 from 0,graph 0 to 0,graph 1 nohead lt 0 lw 2
#set label "$^{120}$Sn, SV-bas, BCS" at graph 0.03,0.09
set label 1 "BCS, $\\tilde{R}_\\mathcal{F}$ = 0 fm" at graph 0.55,0.90
set label 2 "(b)" at graph 0.05,0.9
#p  "varywid_comb.gnu" i 0 u 3:8:2 t "$\\tilde{R}_\\mathcal{F}$=0" w p ls 1
p  "varyocc_widths_3.gnu" i 0 u 3:8:2 t "$\\tilde{R}_\\mathcal{F}$=0" w p ls 1 palette
#
set form x 
unset label
set size 0.8,0.6
set origin 1.2,0.2
set xrange [-47:89]
set yrange [0:4.45]
set ytics 1
set ylabel "" offset 0
set xlabel "" offset 0
set xlabel "s.p. energy  (MeV)" offset 0,0.2
#set key top right box
set arrow 1 from 0,graph 0 to 0,graph 1 nohead lt 0 lw 2
#set label "$^{120}$Sn, SV-bas, HFB" at graph 0.03,0.09
set label 1 "HFB, $\\tilde{R}_\\mathcal{F}$ = 0 fm" at graph 0.55,0.90
set label 2 "(h)" at graph 0.05,0.9
p  "varyocc_widths_3.gnu" i 1 u 3:8:2 t "$\\tilde{R}_\\mathcal{F}$=0" w p ls 1 palette
#
unset label
set size 0.8,0.6
set origin 1.2,0.8
set form x ""
set xrange [-47:89]
set yrange [0:4.45]
set ytics 1
set ylabel "pairing gap $\\Delta_n$ (MeV)" offset 0
set xlabel "" offset 0
#set key bottom right
#set label "$^{120}$Sn, SV-bas, HFB" at graph 0.03,0.09
set label 1 "HFB, $\\tilde{R}_\\mathcal{F}$ = 1 fm" at graph 0.55,0.90
set label 2 "(f)" at graph 0.05,0.9
#p  "varyocc_widths_3.gnu" i 2 u 3:8:2 t "$\\tilde{R}_\\mathcal{F}$=1" w p ls 3 palette
p  "varyocc_widths_3.gnu" i 2  u ($4>3e-8 ? $3:1/0):8:2 t "$\\tilde{R}_\\mathcal{F}$=1" w p ls 3  palette
#
unset colorbox
set size 0.8,0.6
set origin 1.2,1.4
set form x ""
set xrange [-47:89]
set yrange [0:4.45]
set ytics 1
set ylabel "" offset 0
set xlabel "" offset 0
set colorbox user origin  2.05,0.75 size 0.05,1.5
set cblabel "angular momentum $\\ell$" offset 0.2
set label 1 "HFB, $\\tilde{R}_\\mathcal{F}$ = 2 fm" at graph 0.55,0.90
set label 2 "(d)" at graph 0.05,0.9
p  "varyocc_widths_3.gnu" i 3 u 3:8:2 t "$\\tilde{R}_\\mathcal{F}$=2" w p ls 5 palette
#
#!eps2pdf -f varyocc_widths_3color.eps
#!rm -f varyocc_widths_3color.eps
#!evince varyocc_widths_3color.pdf
#!open -a preview  varyocc_widths_3color.pdf
q
