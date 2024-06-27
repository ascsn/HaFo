#!/usr/bin/gnuplot
reset

set terminal postscript eps color enhanced "Helvetica" 22
set output "plotTestCase.eps"
#set format y "10^{%L}"
set mxtics 9
set mytics 5

#set yrange [1e-80:1e1]
#set xrange [1e8:5e11]
#set logscale y
set colorsequence classic
#set ytics nomirror scale 2
set xtics auto
#unset key
#set key center at screen 0.4,0.86

#unset key
set origin 0,0
#set size 0.5,1.
#set rmargin 4
#set bmargin 4
set xlabel "r "
set ylabel "folded {/Symbol j}"
set style line 1 dt 1 lc rgb "black" lw 4 ps 1.5
set style line 2 dt 3 lc rgb "black" lw 4 ps 1.5
set style line 3 dt 1 lc rgb "red" pt 3 lw 4 ps 1.5
set style line 4 dt 3 lc rgb "red" pt 3 lw 4 ps 1.5
set style line 5 dt 1 lc rgb "web-blue" pt 4 lw 4 ps 1.5
set style line 6 dt 3 lc rgb "web-blue" pt 4 lw 4 ps 1.5
set style line 7 dt 1 lc rgb "dark-green" pt 2 lw 4 ps 1.5
set style line 8 dt 3 lc rgb "dark-green" pt 2 lw 4 ps 1.5
set style line 9 dt 1 lc rgb "dark-magenta" lw 4
set style line 10 dt 3 lc rgb "dark-magenta" lw 4
set style line 11 dt 1 lc rgb "medium-blue" lw 4
set style line 12 dt 3 lc rgb "medium-blue" lw 4
set style line 13 dt 3 lc rgb "dark-spring-green" pt 3 lw 4 ps 1.5
set style line 14 dt 3 lc rgb "dark-olivegreen" pt 2 lw 4 ps 1.5
set style line 15 dt 3 lc rgb "royalblue" pt 2 lw 4 ps 1.5
set style line 16 dt 3 lc rgb "dark-grey" pt 2 lw 4 ps 1.5
#set style line 

set key top right font 'Helvetica, 14'

plot "testFoldingData.dat" using 1:3 w lines ls 1 ti "Folded Wavefunction",\
    "testFoldingData.dat" using 1:4 w lines ls 15 ti "Exact Result (eqn 6)"

unset label
unset ylabel
unset ytics
