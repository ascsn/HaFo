set colorsequence classic
set ylabel "V_pair"
p "varyfolding.gnu" i 0 u 1:2 t "HFB" w lp pt 7,\
 "varyfolding.gnu" i 1 u 1:2 t "BCS" w lp pt 7
pause -1
set ylabel "radius"
p "varyfolding.gnu" i 0 u 1:3 t "HFB" w lp pt 7,\
 "varyfolding.gnu" i 1 u 1:3 t "BCS" w lp pt 7
pause -1
set ylabel "energy"
p "varyfolding.gnu" i 0 u 1:4 t "HFB" w lp pt 7,\
 "varyfolding.gnu" i 1 u 1:4 t "BCS" w lp pt 7
pause -1
#set ylabel "pairing energy"
#p "varyfolding.gnu" i 0 u 1:5 t "HFB" w lp pt 7,\
# "varyfolding.gnu" i 1 u 1:5 t "BCS" w lp pt 7
#pause -1
#set ylabel "gap"
#p "varyfolding.gnu" i 0 u 1:6 t "HFB" w lp pt 7,\
# "varyfolding.gnu" i 1 u 1:6 t "BCS" w lp pt 7
#pause -1
set ylabel "pair.pot"
p "varyfolding.gnu" i 0 u 1:7 t "HFB" w lp pt 7,\
 "varyfolding.gnu" i 1 u 1:7 t "BCS" w lp pt 7
pause -1
set ylabel "pair.pot/strength"
p "varyfolding.gnu" i 0 u 1:($7/$2) t "HFB" w lp pt 7,\
 "varyfolding.gnu" i 1 u 1:($7/$2) t "BCS" w lp pt 7
pause -1
set ylabel "pair.dens."
p "varyfolding.gnu" i 0 u 1:(sqrt($8)) t "HFB" w lp pt 7,\
 "varyfolding.gnu" i 1 u 1:(sqrt($8)) t "BCS" w lp pt 7
pause -1
q
#                        HFB                                                    BCS
# 50/70    r   E   E_pair   gap_n   pairpot   pairdens
0.0  260  4.637807  -1024.22   -18.338   1.859   -3.6992   2.053E-04   4.646508  -1023.53  -15.557   1.666   -3.3049    1.639E-04
0.5  370  4.636513  -1021.40   -18.714   1.883   -4.4518   1.446E-04   4.646447  -1023.78  -16.518   1.735   -3.9455    1.137E-04
1.0  950  4.633581  -1009.91   -18.257   1.853   -6.3506   4.468E-05   4.646826  -1023.71  -17.621   1.811   -5.7446    3.656E-05
1.5 3500  4.628548  -978.786   -18.317   1.859   -9.9656   8.106E-06   4.648108  -1023.64  -20.036   1.974   -9.3090    7.074E-06


# BCS
0.0  275  4.646059  -1024.61  -18.125  1.855   -3.5902     1.704E-04
0.5  385  4.646185  -1024.48  -18.217  1.859   -4.1612     1.168E-04
1.0  965  4.646735  -1023.99  -18.379  1.865   -5.8808     3.713E-05
1.5 3410  4.648214  -1023.14  -18.418  1.860   -8.8415     6.722E-06
