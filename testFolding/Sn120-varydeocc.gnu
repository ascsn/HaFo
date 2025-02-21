set terminal po eps enhanced color "TimesNewRoman,30"
set colorsequence classic
set output "Sn120-varydeocc.eps"
set border 15 lw 1
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set size 2.03,2.18
set multiplot
#
set origin 0.2,0.2
set size 0.8,0.8
set xlabel "DEOCC"
set ylabel "energy [MeV]" offset 1.5
unset key 
set xrange [1:19]
#set yrange [-1020.7:-1018.1]
#set ytics 1
p "Summary_120Sn_SVbas.txt"  i 1 u 2:6 t "cut05" w lp pt 7,\
 "Summary_120Sn_SVbas.txt"  i 3 u 2:6 t "cut10" w lp pt 7,\
 "Summary_120Sn_SVbas.txt"  i 5 u 2:6 t "cut20" w lp pt 7,\
 "Summary_120Sn_SVbas.txt"  i 7 u 2:6 t "cut40" w lp pt 7,\
 "Summary_120Sn_SVbas.txt"  i 9 u 2:6 t "cut60" w lp pt 7,\
 "Summary_120Sn_SVbas.txt"  i 0 u 2:6 not w lp lt 0 lc 1 pt 6,\
 "Summary_120Sn_SVbas.txt"  i 2 u 2:6 not w lp lt 0 lc 2 pt 6,\
 "Summary_120Sn_SVbas.txt"  i 4 u 2:6 not w lp lt 0 lc 3 pt 6,\
 "Summary_120Sn_SVbas.txt"  i 6 u 2:6 not w lp lt 0 lc 4 pt 6,\
 "Summary_120Sn_SVbas.txt"  i 8 u 2:6 not w lp lt 0 lc 5 pt 6
#
set origin 0.2,1.0
set size 0.8,0.8
set xlabel ""
set form x ""
set ylabel "charge radius [fm]" offset 1.5
set key at graph 2.05,1.42 maxcolumns 4 maxrows 3
set xrange [1:19]
#set yrange [-1020.7:-1018.1]
#set ytics 1
set title "SV-bas, width=1 fm, pairing tuned to cutoff"
p "Summary_120Sn_SVbas.txt"  i 1 u 2:14 t "cut05 HFB" w lp pt 7,\
 "Summary_120Sn_SVbas.txt"  i 3 u 2:14 t "cut10 HFB" w lp pt 7,\
 "Summary_120Sn_SVbas.txt"  i 5 u 2:14 t "cut20 HFB" w lp pt 7,\
 "Summary_120Sn_SVbas.txt"  i 7 u 2:14 t "cut40 HFB" w lp pt 7,\
 "Summary_120Sn_SVbas.txt"  i 9 u 2:14 t "cut60 HFB" w lp pt 7,\
 "Summary_120Sn_SVbas.txt"  i 0 u 2:14 t "cut05 BCS" w lp lt 0 lc 1 pt 6,\
 "Summary_120Sn_SVbas.txt"  i 2 u 2:14 t "cut10 BCS" w lp lt 0 lc 2 pt 6,\
 "Summary_120Sn_SVbas.txt"  i 4 u 2:14 t "cut20 BCS" w lp lt 0 lc 3 pt 6,\
 "Summary_120Sn_SVbas.txt"  i 6 u 2:14 t "cut40 BCS" w lp lt 0 lc 4 pt 6,\
 "Summary_120Sn_SVbas.txt"  i 8 u 2:14 t "cut60 BCS" w lp lt 0 lc 5 pt 6
#
unset title
set origin 1.2,0.2
set size 0.8,0.8
set xlabel "DEOCC"
set form x 
set ylabel "" offset 1.5
unset key 
set xrange [1:19]
#set yrange [-1020.7:-1018.1]
#set ytics 1
p "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 1 u 2:5 t "cut05" w lp pt 7,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 3 u 2:5 t "cut10" w lp pt 7,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 5 u 2:5 t "cut20" w lp pt 7,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 7 u 2:5 t "cut40" w lp pt 7,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 9 u 2:5 t "cut60" w lp pt 7,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 0 u 2:5 not w lp lt 0 lc 1 pt 6,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 2 u 2:5 not w lp lt 0 lc 2 pt 6,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 4 u 2:5 not w lp lt 0 lc 3 pt 6,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 6 u 2:5 not w lp lt 0 lc 4 pt 6,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 8 u 2:5 not w lp lt 0 lc 5 pt 6
#
set origin 1.2,1.0
set size 0.8,0.8
set xlabel ""
set form x ""
set ylabel "" offset 1.5
#unset key 
set xrange [1:19]
#set yrange [-1020.7:-1018.1]
#set ytics 1
set title "Fy({/Symbol= D}r,HFB), width=0 fm, pairing tuned to cutoff"
p "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 1 u 2:13 t "cut05" w lp pt 7,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 3 u 2:13 t "cut10" w lp pt 7,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 5 u 2:13 t "cut20" w lp pt 7,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 7 u 2:13 t "cut40" w lp pt 7,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 9 u 2:13 t "cut60" w lp pt 7,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 0 u 2:13 not w lp lt 0 lc 1 pt 6,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 2 u 2:13 not w lp lt 0 lc 2 pt 6,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 4 u 2:13 not w lp lt 0 lc 3 pt 6,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 6 u 2:13 not w lp lt 0 lc 4 pt 6,\
 "Summary_120Sn_FyDr_wid0_tsym_T.txt"  i 8 u 2:13 not w lp lt 0 lc 5 pt 6
#
!eps2pdf -f  Sn120-varydeocc.eps
!rm -f  Sn120-varydeocc.eps
#!open -a preview Sn120-varydeocc.pdf
!evince Sn120-varydeocc.pdf
q
