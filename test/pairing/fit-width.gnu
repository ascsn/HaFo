set terminal po eps enhanced color "TimesNewRoman,30"
set colorsequence classic
set output "fit-width.eps"
set border 15 lw 1
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set size 1.03,2.36
set multiplot
#
set origin 0.2,0.2
set size 0.8,0.8
set xlabel "HFB folding width  [fm]"
set xrange [-0.02:1.22]
#set xtics 0.1
set ylabel "quality measure {/Symbol= c}^2" offset 0.4
set key top left
set title "SV-min like fit, varying {/Symbol= G}_{fold} for {/Symbol= e}_{cut} = 30 MeV"
p "fit-width.gnu" i 0 u 1:($2) t "SV-min like, HFB" w lp lw 3 pt 7 ps 1.5,\
 "fit-width.gnu" i 1 u 1:($2) t "SV-min like, BCS" w lp lt 2 lc rgb '#00aa00' lw 3 pt 7 ps 1.5
#
set origin 0.2,1.4
set size 0.8,0.8
set xrange [-0.02:1.22]
#set xtics 0.1
set ylabel "quality measure {/Symbol= c}^2" offset 0.4
set key top right
set title "Fy({/Symbol= D}r,HFB) like fit, HFB, for two {/Symbol= e}_{cut}"
p "fit-width.gnu" i 2 u 1:($2) t "Fy({/Symbol= D}r,HFB) like, {/Symbol= e}_{cut}=30" w lp lw 3 pt 7 ps 1.5,\
 "fit-width.gnu" i 3 u 1:($2) t "Fy({/Symbol= D}r,HFB) like, {/Symbol= e}_{cut}=15" w lp lt 2 lc rgb '#00aa00' lw 3 pt 7 ps 1.5
# "fit-width.gnu" i 4 u 1:($2) not w lp lt 0 lc 1 lw 3 pt 6 ps 1.5,\
#
!epstopdf fit-width.eps
!rm -f  fit-width.eps
#!open -a preview fit-width.pdf
!evince fit-width.pdf
q
# SVmin fit, HFB cut30, old cut
#wid   chi**2
0.0   55.92
0.3   54.943
0.5   54.214
0.6   54.02
0.7   53.91
0.8   54.107
0.9   54.441
1.0   55.084
1.2   56.89


# SVmin fit, BCS cut30, old cut
#wid   chi**2
0.0   55.905    
0.2   55.011    
0.3   54.892    
0.5   54.657    
0.7   54.801    
0.8   55.155    
1.0   56.305    
1.2   58.000
1.5   63.344    


# FyDrHFB cut30 fit varying width
#wid   chi**2
0.0   78.391
0.5   69.62
0.7   64.161
0.8   58.603    
0.9   55.434    
1.0   57.505    
1.1   60.89
1.2   68.86


# FyDrHFB cut15 fit varying width
#wid   chi**2
0.0  65.633    
0.4  60.87
0.5  57.93
0.6  55.291
0.7  52.002
0.8  50.46
0.9  53.797
1.0  56.67
1.1  60.166
1.2  67.292


# FyDrHFB cut30 fit varying width, new cutoff
#wid  4sh         6sh
 0   66.032     69.926    
0.4  64.341    
0.5  63.101    74.303?
0.8  57.57
0.9  54.352    
1.0  55.007    
1.1  59.029    
1.2  65.3
