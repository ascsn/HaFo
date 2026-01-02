set terminal po eps enhanced color "TimesNewRoman,33"
set colorsequence classic
set output "Sn120-varydeocc_cut15.eps"
set border 15 lw 1
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set size 1.43,1.96
set multiplot
#
set origin 0.2,0.2
set size 0.6,0.8
set xlabel "{/Symbol= D}_{shell}"
set ylabel "neutron pairing energy [MeV]" offset 1.5
set key at graph 0.99,0.98 spacing 1.1
set xrange [1.5:10.5]
#set yrange [-0.07:0.005]
set yrange [-1:19.9]
#set ytics 1
p "../collect_cat/collect_cat_FyDr_120Sn_cut15_wid0_0.0d0_hfb.dat"  i 0  u 2:($8+26.830) t "{/Symbol= G}_{fold}=0, HFB" w lp lw 4 pt 7,\
 "../collect_cat/collect_cat_FyDr_120Sn_cut15_wid0_0.0d0_bcs.dat"  i 0  u 2:($8+19.219) t "{/Symbol= G}_{fold}=0, BCS" w lp lw 4 lc rgb '#00aa00' pt 7,\
 "../collect_cat/collect_cat_FyDr_120Sn_cut15_wid1_1.0d0_hfb.dat"  i 0  u 2:($8+19.973) t "{/Symbol= G}_{fold}=1 fm, HFB" w lp lw 4 pt 7,\
 "../collect_cat/collect_cat_FyDr_120Sn_cut15_wid1_1.0d0_bcs.dat"  i 0  u 2:($8+18.124) t "{/Symbol= G}_{fold}=1 fm, BCS" w lp lw 4 pt 7
#
set origin 0.2,1.0
set size 0.6,0.8
set xlabel ""
set form x ""
set ylabel "charge radius [10^{-3}fm]" offset 1
#set key at graph 2.15,1.28 maxcolumns 4 maxrows 1
set xrange [1.5:10.5]
#set yrange [-0.015:0.07]
set yrange [-22:2]
#set ytics 1
set key at graph 0.99,0.36
set title "Fy({/Symbol= D}r), ^{120}Sn, {/Symbol= e}_{cut} = 15 MeV"
p "../collect_cat/collect_cat_FyDr_120Sn_cut15_wid0_0.0d0_hfb.dat"  i 0  u 2:(($5-4.668917)*1e3) t "{/Symbol= G}_{fold}=0, HFB" w lp lw 4 pt 7,\
 "../collect_cat/collect_cat_FyDr_120Sn_cut15_wid0_0.0d0_bcs.dat"  i 0  u 2:(($5-4.659498)*1e3) t "{/Symbol= G}_{fold}=0, BCS" w lp lw 4 lc rgb '#00aa00' pt 7,\
 "../collect_cat/collect_cat_FyDr_120Sn_cut15_wid1_1.0d0_hfb.dat"  i 0  u 2:(($5-4.661443)*1e3) t "{/Symbol= G}_{fold}=1 fm, HFB" w lp lw 4 pt 7,\
 "../collect_cat/collect_cat_FyDr_120Sn_cut15_wid1_1.0d0_bcs.dat"  i 0  u 2:(($5-4.652495)*1e3) t "{/Symbol= G}_{fold}=1 fm, BCS" w lp lw 4 pt 7
#
unset title
set origin 0.8,0.2
set size 0.6,0.8
set xlabel "{/Symbol= D}_{shell}"
set form x 
set ylabel "" offset 1.5
set form y ""
set key at graph 0.99,0.98 invert
set xrange [1.5:10.5]
set yrange [-1:19.9]
#set ytics 1
p "../collect_cat/collect_cat_FyDr_120Sn_cut45_wid0_0.0d0_hfb.dat"  i 0  u 2:($8+27.097) t "{/Symbol= G}_{fold}=0, HFB" w lp lw 4 pt 7,\
 "../collect_cat/collect_cat_FyDr_120Sn_cut45_wid0_0.0d0_bcs.dat"  i 0  u 2:($8+38.004) t "{/Symbol= G}_{fold}=0, BCS" w lp lc rgb '#00aa00' lw 4 pt 7,\
 "../collect_cat/collect_cat_FyDr_120Sn_cut45_wid1_1.0d0_hfb.dat"  i 0  u 2:($8+24.922) t "{/Symbol= G}_{fold}=1 fm, HFB" w lp lw 4 pt 7,\
 "../collect_cat/collect_cat_FyDr_120Sn_cut45_wid1_1.0d0_bcs.dat"  i 0  u 2:($8+24.157) t "{/Symbol= G}_{fold}=1 fm, BCS" w lp lw 4 pt 7
#
set origin 0.8,1.0
set size 0.6,0.8
set xlabel ""
set form x ""
set ylabel "" offset 1.5
#unset key 
set key at graph 0.99,0.36 noinvert
set xrange [1.5:10.5]
set yrange [-22:2]
#set ytics 1
set title "Fy({/Symbol= D}r), ^{120}Sn, {/Symbol= e}_{cut} = 45 MeV"
p "../collect_cat/collect_cat_FyDr_120Sn_cut45_wid0_0.0d0_hfb.dat"  i 0  u 2:(($5-4.669307)*1e3) t "{/Symbol= G}_{fold}=0, HFB" w lp lw 4 pt 7,\
 "../collect_cat/collect_cat_FyDr_120Sn_cut45_wid0_0.0d0_bcs.dat"  i 0  u 2:(($5-4.675436)*1e3) t "{/Symbol= G}_{fold}=0, BCS" w lp lc rgb '#00aa00' lw 4 pt 7,\
 "../collect_cat/collect_cat_FyDr_120Sn_cut45_wid1_1.0d0_hfb.dat"  i 0  u 2:(($5-4.662187)*1e3) t "{/Symbol= G}_{fold}=1 fm, HFB" w lp lw 4 pt 7,\
 "../collect_cat/collect_cat_FyDr_120Sn_cut45_wid1_1.0d0_bcs.dat"  i 0  u 2:(($5-4.661550)*1e3) t "{/Symbol= G}_{fold}=1 fm, BCS" w lp lw 4 pt 7
#
!eps2pdf -f  Sn120-varydeocc_cut15.eps
!rm -f  Sn120-varydeocc_cut15.eps
#!open -a preview Sn120-varydeocc_cut15.pdf
!evince Sn120-varydeocc_cut15.pdf
q
