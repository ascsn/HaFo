#!/bin/bash
# Batch file for systematic variation of cutoff and folding width.
# Test case is 120Sn with 51 radial grid points and 17 shell above
# the Fermi energy. Proton pairing is practically switched off.
# The input parameters are:
#  $1   = cutoff energy 'cutfac'
#  $2   = folding width 'hfbwidth'
#  $3   = pairing strength '-fxi'
# A note on the pairing strength: The Fayans functional has three (four)
# pairing parameters. Variation of pairing strength has thus different
# direction. We vary here only the first pairing parameter 'fxi'.  This
# parameter is negative. To have a üpositive entry, we enter the absolute
# value and augment the minus sign in the procedure.

# Set HFB widths
#hfbwidths=("0.0d0" "1.0d0" "1.5d0" "2.0d0" "3.0d0")

#ipair=(7)    # for BCS
ipair=(7 8)    # for HFB

# Set hfbwidth
hfbwidth="0.0d0"

pairStrengthsBCS=("3.1" "7.08" "17.35" "39.85" "117.8") # 120Sn, cut100
pairStrengthsHFB=("2.9383" "6.181" "13.651" "35.142" "104.24") # Sn120 Cut100


for j in "${ipair[@]}"; do
if [ $j -eq 7 ]; then
    cond="bcs"
    #pairStrength="119.4"
    #combined=$(paste -d "," <(printf "%s\n" "${hfbwidths[@]}") <(printf "%s\n" "${pairStrengthsBCS[@]}"))
else
    cond="hfb"
    #pairStrength="107"
    #combined=$(paste -d "," <(printf "%s\n" "${hfbwidths[@]}") <(printf "%s\n" "${pairStrengthsHFB[@]}"))
fi
#while IFS=, read -r hfbwidth pairStrength; do
echo "cond = ${j}, hfbwidth = ${hfbwidth}, pairStrength = ${pairStrength}" 
for i in $(seq 3 1 20); do
echo "1,,    IWUNIT" > for005
echo "'1D FyDrHFB_cutoff15_wid0'" >> for005
echo "   0   0.300     NGRID,RSTEP" >> for005
echo "  -1  1.00000E-10  1.00000E-10    200   9820   3 -1    iwr,dabr,dabe,itpr,itmax,ifplot" >> for005
echo "     0.20   110.00     0.15    -1.00  x0damp,e0damp,addnew,wall" >> for005
echo " 1, 0, 0,-1, 1, 3, 1, 0, 0,13, if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk" >> for005
echo "1,,,                       ifform" >> for005
echo "0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12" >> for005
echo "0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12" >> for005
echo "  -9.83607919134      0.616867867133      0.219426189580     avp,h1vp,h2vp" >> for005
echo "   23.7133523672     -0.294456682777       23.2671270000    " >> for005
echo "  0.333333330000      8.000000000000E-02" >> for005
echo "  0.564427600000       0.00000000000      0.448923060000    " >> for005
echo "   0.00000000000       23.2671270000       0.00000000000    " >> for005
echo "   0.00000000000      0.333333333333    " >> for005
echo "  -4.45430330000       4.24374220000       3.07895230000      0.666666670000      0.160000000000    " >> for005
echo "   0.00000000000      0.300000000000       0.00000000000       0.00000000000     2" >> for005
echo "   0.00000000000       0.00000000000    " >> for005
echo "  0.197210090000      5.319830200000E-02   0.00000000000       0.00000000000    " >> for005
echo "    0.00000000000    ,1.0, ${j},   0.00000000000    ,   0.00000000000    ,0,,  -15.0D0,${hfbwidth},   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
echo "   50 70    0    0    ${i}  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
echo "0, ifgsc" >> for005
echo "-1, to stop" >> for005


../hafom_huge.exe < for005 > /dev/null

echo "cutoff=-${i},hfbwidth=${hfbwidth},shells=${i}" >> collect_cat_FyDr_120Sn_cut15_wid0_${hfbwidth}_${cond}_shell_Cut100.dat
grep Z,N for011 >>  collect_cat_FyDr_120Sn_cut15_wid0_${hfbwidth}_${cond}_shell_Cut100.dat

done
echo "${cond} Done!"
mv collect_cat_FyDr_120Sn_cut15_wid0_${hfbwidth}_${cond}_shell_Cut100.dat FyDrTests/
#done <<< "${combined}"
done