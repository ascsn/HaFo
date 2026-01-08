set terminal po eps enhanced color "TimesNewRoman,30"
set colorsequence classic
set output "CoD_isovector_surface.eps"
set border 15 lw 1
set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0
set size 1.03,2.23
set style fill solid 0.2
set boxwidth 0.8
set multiplot
#
set origin 0.2,0.4
set size 0.8,0.8
set ylabel "CoD with asm" offset 00
set xtics ("DE(236U)  " 19,"E(236U*)  " 18,"E(236U)   " 17,"asm" 16,\
"delhgradxi" 15,"delh1xi   " 14,"delfxi    " 13,"hgradxi   " 12,\
"h1xi      " 11,"fxi       " 10,"hgradp    "  9,"asp       "  8,\
"c1nabj    "  7,"c0nabj    "  6,"L         "  5,"J         "  4,\
"incompr. K"  3,"rho0-sym  "  2,"E/A-sym   "  1) offset 0,-5 rotate by 90
unset key 
set xrange [0.5:19.5]
set yrange [0:1]
set arrow from 5.5,0 to 5.5,1 nohead lt 0 lw 2
set arrow from 7.5,0 to 7.5,1 nohead lt 0 lw 2
set arrow from 9.5,0 to 9.5,1 nohead lt 0 lw 2
set arrow from 15.5,0 to 15.5,1 nohead lt 0 lw 2
set arrow from 16.5,0 to 16.5,1 nohead lt -1 lw 2
set label "bulk" at 3,0.8 center rotate by 0
set label "l*s" at 6.5,0.8 center
set label "T=0 surface" at 8.5,0.77 center rotate by 90
set label "pairing" at 12.5,0.8 center rotate by 0
set label "T=1 surface" at 16,0.77 center rotate by 90
#set ytics 1
#set label "Fy(IVP3), HFBwidth=0.9 fm" at graph 0.2,0.92
p "CoD_isovector_surface.gnu" i 0 u 1:($3**2)  not w boxes lw 3
#
set origin 0.2,1.25
set size 0.8,0.8
set ylabel "CoD with E^*(^{236}U)" offset 00
set xtics ("" 19,"" 18,"" 17,"" 16,"" 15,"" 14,"" 13,"" 12,\
"" 11,"" 10,""  9,""  8,""  7,""  6,""  5,""  4,\
""  3,""  2,""  1) offset 0,-5 rotate by 90
unset key 
set xrange [0.5:19.5]
set yrange [0:1]
#set ytics 1
#set label "Fy(IVP3), HFBwidth=0.9 fm" at graph 0.2,0.92
set title "Fy(IVP3) HFBwidth=0.9 fm with free asm=a@^{(s)}_-"
p "CoD_isovector_surface.gnu" i 1 u 1:($3**2)  not w boxes lw 3
#
!eps2pdf -f  CoD_isovector_surface.eps
!rm -f  CoD_isovector_surface.eps
#!open -a preview CoD_isovector_surface.pdf
!evince CoD_isovector_surface.pdf
q
# align asm
#alignement for observable: 16 asm                          
# average,variance=  5.00000E-03  1.66232E-02
   19   19 -0.47743DE(236U)                      
   18   18 -0.09682E(236U*)                      
   17   17  0.21444E(236U)                       
   16   16  1       asm
   15   15  0.21298 delhgradxi                   
   14   14  0.21947 delh1xi                      
   13   13 -0.28005 delfxi                       
   12   12  0.40411 hgradxi                      
   11   11 -0.29159 h1xi                         
   10   10  0.19255 fxi                          
    9    9 -0.45121 hgradp                       
    8    8 -0.71303 asp                          
    7    7  0.20424 c1nabj                       
    6    6  0.57354 c0nabj                       
    5    5  0.12476 L                            
    4    4  0.38566 J                            
    3    3 -0.57754 incompr. K                   
    2    2  0.28208 rho0-sym                     
    1    1  0.48685 E/A-sym                      


# align isomer
#alignement for observable: 19DE(236U)                      
# average,variance=   2.3302      0.37466    
   19   19  1
   18   18  0.55584E(236U*)                      
   17   17 -0.02457E(236U)                       
   16   16 -0.47743 asm                          
   15   15 -0.66599 delhgradxi                   
   14   14  0.29976 delh1xi                      
   13   13 -0.05474 delfxi                       
   12   12 -0.32774 hgradxi                      
   11   11  0.48322 h1xi                         
   10   10 -0.37307 fxi                          
    9    9  0.32623 hgradp                       
    8    8  0.52455 asp                          
    7    7  0.29883 c1nabj                       
    6    6 -0.39072 c0nabj                       
    5    5  0.10633 L                            
    4    4  0.13704 J                            
    3    3  0.21820 incompr. K                   
    2    2 -0.09820 rho0-sym                     
    1    1 -0.74151 E/A-sym                      
