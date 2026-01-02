#set terminal po eps enhanced color "TimesNewRoman,30"
set terminal epslatex color colortext font 17 solid round standalone
set colorsequence classic
#set output "fit-width_FyIVP3.eps"
set output "fit-width_FyIVP3.tex"
set border 15 lw 1
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set size 1.03,1.74
set multiplot
#
set origin 0.2,0.21
set size 0.8,1.2
set xlabel "folding radius $\\tilde{R}_\\mathcal{F}$ (fm)"
set xrange [-0.02:1.1]
#set xrange [-0.02:1.22]
#set yrange [50.2:57]
set yrange [38.3:45.5]
set ytics 1
set ylabel "quality measure $\\chi^2$" offset 0.0,3
set key at graph 0.43,0.12
#set title "SV-min like, HFB"
set label 1 "Fy(IVP3)" at graph 0.45,0.9
set label 2 "(b)" at graph 0.05, 0.05
p "fit-width_FyIVP3.gnu" i 2 u 1:($2) smooth cspline not w l lt 1 lw 4,\
#
set origin 0.2,1.41
set size 0.8,0.3
set xlabel ""
set form x ""
set xrange [-0.02:1.1]
#set xrange [0.48:1.12]
#set yrange [53.7:55.2]
set yrange [53.8:55.45]
set ytics 1
set ylabel "" offset 0.4
unset key 
set label 1 "SV-min like, HFB" at graph 0.33,0.75
set label 2 "(a)" at graph 0.05,0.15
p "fit-width_FyIVP3.gnu" i 1 u 1:($2) smooth csplines not  w l lw 4
#
#!eps2pdf -f  fit-width_FyIVP3.eps
#!rm -f  fit-width_FyIVP3.eps
#!open -a preview fit-width_FyIVP3.pdf
#!evince fit-width_FyIVP3.pdf
q
#0: SVmin v**2 cutoff 30
0.0   55.07
0.3   54.646    
0.5   54.189    
0.6   53.982    
0.7   53.848    
0.8   54.071    
0.9   54.262    
1.0   54.591    
1.2   56.308    


#1: SVmin v**2 cutoff 15
0.0   55.32
0.3   54.801    
0.5   54.472    
0.6   54.236    
0.7   54.148    
0.8   54.303    
0.9   54.769    
1.0   55.130    
1.2   56    *


#2: FyIVP3
0.0  45.162    
0.1  43.193    
0.4  40.4   40.914    
0.6  39.299    
0.8  38.643    
0.9  38.413    
1.0  39.019    
1.1  40.888    







