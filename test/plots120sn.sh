########################  Energy ###############################################
reset
set terminal postscript enhanced color font "Bold-Helvetica,22" portrait
#set size ratio 1.61
set lmargin 0;set rmargin 0;set tmargin 0;set bmargin 0

set output '120SnProperties.eps'
set border 31 lw 2.0

x = 1.0
y = 1.0
xml = 0.20
xmr = 0.06
ymb = 0.14
ymt = 0.04

set multiplot


#-----1-----------

#---Ye = 0.5

set size x-xml-xmr,(y-ymb-ymt)/2
set origin xml,ymb
set xlabel "/Symbol{e} (MeV)" 
set label "E (MeV)" at screen xml/5,y/2.2 rotate left

unset key
unset grid

set xrange [10:101]
#set yrange [-14.0:70]

set mxtics 10
set xtics scale 2
set mytics 4
set ytics scale 2

#set logscale x
#set xtics format '10^{%T}'

#set label '(b)' at 1.4e-4,47.5
#set label 'T = 4 MeV' at 1.4e-4,62.5 #font 'bold-times-roman, 22'
#set label 'Ye = 0.4' at 1.4e-4,55.5 #font 'bold-times-roman, 22'

#set label 'E_o' at 8e-4,5.5 font 'bold-times-roman, 22' textcolor rgb 'red'
#set label 'E_d' at 8e-3,8.0 font 'bold-times-roman, 22' textcolor rgb 'brown'
#set label 'E_{/Symbol t}' at 5e-3,13.0 font 'bold-times-roman, 22' textcolor rgb 'blue'
#set label 'E_{3}' at 3e-3,-10.0 font 'bold-times-roman, 22' textcolor rgb 'dark-green'
#set label 'E_{/Symbol a}' at 0.008,-8.0 font 'bold-times-roman, 22' textcolor rgb 'olive-green'
#set label 'E_e' at 0.01,40.0 font 'bold-times-roman, 22' textcolor rgb 'brown4'
#set label 'E_{tot}' at 2e-2,27.5 font 'bold-times-roman, 22' textcolor rgb 'black'
#set arrow head filled from 5e-3,12.0 to 4e-3,-1.0 lc rgb 'black'

set style line 1 dt 1 lc rgb "red" lw 4
set style line 2 dt 1 lc rgb "brown" lw 4
set style line 3 dt 1 lc rgb "blue" lw 4
set style line 4 dt 1 lc rgb "dark-green" lw 4
set style line 5 dt 1 lc rgb "olive-green" lw 4
set style line 6 dt 1 lc rgb "brown4" lw 4
set style line 7 dt 2 lc rgb "blue" lw 4
set style line 8 dt 1 lc rgb "black" lw 4

p "collect_cat/collect_cat_Svbas_120Sn_0.0d0_bcs_Cut100" u 1:2 w l ls 1,"EnergyT4Ye0.4.dat" u 1:($3*100) w l ls 2,"EnergyT4Ye0.4.dat" u 1:($4*20) w l ls 3,"EnergyT4Ye0.4.dat" u 1:($5*1000) w l ls 4,"EnergyT4Ye0.4.dat" u 1:6 w l ls 5,"EnergyT4Ye0.4.dat" u 1:7 w l ls 6,"EnergyT4Ye0.4.dat" u 1:9 w l ls 8


unset ytics
unset label
unset xlabel
unset arrow





#----------------------2----------------

#---Ye = 0.1

set origin xml,(y+ymb-ymt)/2
unset ylabel
set xtics format ''
#set ytics format '10^{%T}'
set xrange [1e-4:0.1]
set yrange [-1.9:13.9]

set mxtics 10
set xtics scale 2
set mytics 4
set ytics scale 2

set logscale x

#set label 'APR' at 0.00015,12.0 font 'bold-times-roman, 22'
#set arrow from 0.0003,12.0 to 0.0007,12.0 nohead lw 2.0 dt 1
#set label 'Virial' at 0.00015,11.0 font 'bold-times-roman, 22'
#set arrow from 0.0003,11.0 to 0.0007,11.0 nohead lw 2.0 dt 2

set label '(a)' at 1.4e-4,9.5

set label 'T = 4 MeV' at 1.4e-4,12.5 #font 'bold-times-roman, 22'
set label 'Ye = 0.1' at 1.4e-4,11. #font 'bold-times-roman, 22'

set label 'E_o' at 4e-2,2.5 font 'bold-times-roman, 22' textcolor rgb 'red'
set label 'E_d' at 9e-3,1.8 font 'bold-times-roman, 22' textcolor rgb 'brown'
set label 'E_{/Symbol t}' at 3e-2,-1. font 'bold-times-roman, 22' textcolor rgb 'blue'
set label 'E_{3}' at 1.2e-3,0.9 font 'bold-times-roman, 22' textcolor rgb 'dark-green'
set label 'E_{/Symbol a}' at 3e-4,-1.0 font 'bold-times-roman, 22' textcolor rgb 'olive-green'
set label 'E_e' at 2e-2,8.0 font 'bold-times-roman, 22' textcolor rgb 'brown4'
set label 'E_{tot}' at 8e-4,7.5 font 'bold-times-roman, 22' textcolor rgb 'black'

set style line 1 dt 1 lc rgb "red" lw 4
set style line 2 dt 1 lc rgb "brown" lw 4
set style line 3 dt 1 lc rgb "blue" lw 4
set style line 4 dt 1 lc rgb "dark-green" lw 4
set style line 5 dt 1 lc rgb "olive-green" lw 4
set style line 6 dt 1 lc rgb "brown4" lw 4
set style line 7 dt 1 lc rgb "blue" lw 4
set style line 8 dt 1 lc rgb "black" lw 4

plot "EnergyT4Ye0.1.dat" u 1:2 w l ls 1,"EnergyT4Ye0.1.dat" u 1:($3*100) w l ls 2,"EnergyT4Ye0.1.dat" u 1:($4*5) w l ls 3,"EnergyT4Ye0.1.dat" u 1:5 w l ls 4,"EnergyT4Ye0.1.dat" u 1:6 w l ls 5,"EnergyT4Ye0.1.dat" u 1:7 w l ls 6,"EnergyT4Ye0.1.dat" u 1:9 w l ls 8




unset multiplot




