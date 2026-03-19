set terminal po eps enhanced color "TimesNewRoman,30"
set colorsequence classic
set output "MCC_isovector_surface.eps"
set border 15 lw 1
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set size 1.63,1.38
set style fill solid 0.2
set boxwidth 0.8
set multiplot
#
set origin 0.2,0.4
set size 0.4,0.8
set ylabel "MCC with asm" offset 00
set xtics ("surface" 1,"l*s" 2,"pairing" 3,"bulk T=1" 4) offset 0,-5 rotate by 90
unset key 
set xrange [0.5:4.5]
set yrange [0:1]
#set ytics 1
#set label "Fy(IVP3), HFBwidth=0.9 fm" at graph 0.2,0.92
p "MCC_isovector_surface.gnu" u 1:3  not w boxes lw 3
#
set origin 0.8,0.4
set size 0.8,0.8
set boxwidth 0.2
set ylabel "MCC with E^*(^{236}U)" offset 1.2
set xtics ("surface" 1,"l*s" 2,"pairing" 3,"bulk T=1" 4) offset 0,-5 rotate by 90
#set xtics (""  4,""  3,""  2,""  1) offset 0,-5 rotate by 90
set key top left spacing 1.2
set xrange [0.5:4.5]
set yrange [0:1]
#set ytics 1
#set label "Fy(IVP3), HFBwidth=0.9 fm" at graph 0.2,0.92
set title "Fy(IVP3) HFBwidth=0.9 fm with free asm=a@^{(s)}_-" offset -12
p "MCC_isovector_surface.gnu"  u ($1-0.3):4  t "^{236}U" w boxes lw 3,\
 "MCC_isovector_surface.gnu"  u ($1-0.1):5   t "^{238}U" w boxes lw 3,\
 "MCC_isovector_surface.gnu"  u ($1+0.1):6  t "^{240}Pt" w boxes lw 3,\
 "MCC_isovector_surface.gnu"  u ($1+0.3):7  t "^{242}Cm" w boxes lw 3
#
!eps2pdf -f  MCC_isovector_surface.eps
!rm -f  MCC_isovector_surface.eps
#!open -a preview MCC_isovector_surface.pdf
!evince MCC_isovector_surface.pdf
q
# MCC 
# nr  group   asm   E*(236U)  DE(238U)   DE(240Pt)  DE(242Cm)
1 surface    0.696  0.385    0.415   0.308   0.232 
#1 surface'  0.696  0.390    0.416   0.309   0.235
2 ls         0.380  0.227    0.270   0.204   0.128
3 pairing    0.347  0.489    0.451   0.588   0.639
4 J/L        0.316  0.003    0.029   0.006   0.018


