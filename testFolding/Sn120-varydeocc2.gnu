set terminal po eps enhanced color "TimesNewRoman,30"
set colorsequence classic
set output "Sn120-varydeocc2.eps"
set border 15 lw 1
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set size 2.03,2.09
set multiplot
#
set origin 0.2,0.2
set size 0.8,0.8
set xlabel "DEOCC"
set ylabel "neutr. pairing energy [MeV]" offset 1.5
unset key 
set xrange [1.5:10.5]
set yrange [-0.2:1.75]
#set ytics 1
p "collect_cat/collect_cat_FyDr_120Sn_cut30_wid0_0.0d0_hfb_shell.dat"  i 0  u 1:($5+25.144) t "wid=0,HFB" w lp pt 7,\
 "collect_cat/collect_cat_FyDr_120Sn_cut30_wid0_0.0d0_bcs_shell.dat"  i 0  u 1:($5+13.008) t "wid=0,BCS" w lp pt 7,\
 "collect_cat/collect_cat_FyDr_120Sn_cut30_wid0_1.0d0_hfb.dat"  i 0  u 1:($5+22.006) t "wid=1,HFB" w lp pt 7,\
 "collect_cat/collect_cat_FyDr_120Sn_cut30_wid0_1.0d0_bcs.dat"  i 0  u 1:($5+17.244) t "wid=1,BCS" w lp pt 7
#
set origin 0.2,1.0
set size 0.8,0.8
set xlabel ""
set form x ""
set ylabel "charge radius [10^{-3}fm]" offset 1
set key at graph 2.15,1.28 maxcolumns 4 maxrows 1
set xrange [1.5:10.5]
set yrange [-2.2:0.2]
#set ytics 1
set title "Fy({/Symbol= D}r), ^{120}Sn"
p "collect_cat/collect_cat_FyDr_120Sn_cut30_wid0_0.0d0_hfb_shell.dat"  i 0  u 1:(($2-4.664499)*1e3) t "width=0,HFB" w lp pt 7,\
 "collect_cat/collect_cat_FyDr_120Sn_cut30_wid0_0.0d0_bcs_shell.dat"  i 0  u 1:(($2-4.654200)*1e3) t "width=0,BCS" w lp pt 7,\
 "collect_cat/collect_cat_FyDr_120Sn_cut30_wid0_1.0d0_hfb.dat"  i 0  u 1:(($2-4.657514)*1e3) t "width=1,HFB" w lp pt 7,\
 "collect_cat/collect_cat_FyDr_120Sn_cut30_wid0_1.0d0_bcs.dat"  i 0  u 1:(($2-4.646520)*1e3) t "width=1,BCS" w lp pt 7
#
unset title
set origin 1.2,0.2
set size 0.8,0.8
set xlabel "DEOCC"
set form x 
set ylabel "" offset 1.5
unset key 
set xrange [1.5:10.5]
set yrange [-0.2:1.3]
#set ytics 1
p "collect_cat/collect_cat_SVMin_120Sn_cut30_HFB_wid0_0.0d0_hfb.dat"  i 0  u 1:($5+15.481) t "cut05" w lp pt 7,\
 "collect_cat/collect_cat_SVMin_120Sn_cut30_BCS_wid0_0.0d0_bcs.dat"  i 0  u 1:($5+17.609) t "cut05" w lp pt 7,\
 "collect_cat/collect_cat_SVMin_120Sn_cut30_HFB_wid1_1.0d0_hfb.dat"  i 0  u 1:($5+16.102) t "cut05" w lp pt 7,\
 "collect_cat/collect_cat_SVMin_120Sn_cut30_BCS_wid1_1.0d0_bcs.dat"  i 0  u 1:($5+17.214) t "cut05" w lp pt 7
#
set origin 1.2,1.0
set size 0.8,0.8
set xlabel ""
set form x ""
set ylabel "" offset 1.5
#unset key 
set xrange [1.5:10.5]
set yrange [-0.28:0.04]
#set ytics 1
set title "SV-min, ^{120}Sn"
p "collect_cat/collect_cat_SVMin_120Sn_cut30_HFB_wid0_0.0d0_hfb.dat"  i 0  u 1:(($2-4.645424)*1e3) t "cut05" w lp pt 7,\
 "collect_cat/collect_cat_SVMin_120Sn_cut30_BCS_wid0_0.0d0_bcs.dat"  i 0  u 1:(($2-4.644230)*1e3) t "cut05" w lp pt 7,\
 "collect_cat/collect_cat_SVMin_120Sn_cut30_HFB_wid1_1.0d0_hfb.dat"  i 0  u 1:(($2-4.643109)*1e3) t "cut05" w lp pt 7,\
 "collect_cat/collect_cat_SVMin_120Sn_cut30_BCS_wid1_1.0d0_bcs.dat"  i 0  u 1:(($2-4.641120)*1e3) t "cut05" w lp pt 7
#
!eps2pdf -f  Sn120-varydeocc2.eps
!rm -f  Sn120-varydeocc2.eps
#!open -a preview Sn120-varydeocc2.pdf
!evince Sn120-varydeocc2.pdf
q
