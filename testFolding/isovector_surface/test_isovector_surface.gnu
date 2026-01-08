set terminal po eps enhanced color "TimesNewRoman,30"
set colorsequence classic
set output "test_isovector_surface.eps"
set border 15 lw 1
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set size 1.03,1.03
set multiplot
#
set origin 0.18,0.2
set size 0.8,0.8
set xlabel "isovector surface parameter 'asm'"
set ylabel "{/Symbol= c}^2" offset 00
unset key 
set xrange [-0.03:0.04]
#set yrange [-1020.7:-1018.1]
#set ytics 1
set label "Fy(IVP3), HFBwidth=0.9 fm" at graph 0.2,0.92
p "test_isovector_surface.gnu" u 1:2 smooth csplines not w l lw 3,\
  38.3 not w l lt 0 lw 2,  39.3 not w l lt 0 lw 2
#
!eps2pdf -f  test_isovector_surface.eps
!rm -f  test_isovector_surface.eps
#!open -a preview test_isovector_surface.pdf
!evince test_isovector_surface.pdf
q
# HFBwidtgh=0.9, set Fy(IVP3)
# asm      chi**2
-0.03      43.74
-0.02      41.018  41.064
-0.01      39.250
 0.0       38.413
 0.01      38.407
 0.02      39.199
 0.04      42.739




