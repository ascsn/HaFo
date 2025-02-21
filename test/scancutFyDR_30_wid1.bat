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
hfbwidth="1.0d0"

#pairStrengthsBCS=("3.29" "7.85" "18.27" "42.81" "130.07") # 120Sn, cut100
#pairStrengthsHFB=("3.0304" "7.137" "14.765" "37.38" "116") # Sn120 Cut100


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
echo "'1D: FyDrHFB_cutoff30_wid1'" >> for005
echo "0,0.3,              NGRID,RSTEP" >> for005
echo "-1,1E-8,1E-9,278,5780,3,-1, IWR,DABR/E,ITPR/MAX,IFPL" >> for005
echo "0.2,90.,0.2,-1,,               X0DAMP,E0DAMP,ADDNEW,WALL" >> for005
echo " 1, 0, 0,-1, 1, 3, 1, 0, 0,13, if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk" >> for005
echo "1,,,                       ifform" >> for005
echo "0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12" >> for005
echo "0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12" >> for005
echo "  -10.33326005      0.6043403239      0.3238756404     avp,h1vp,h2vp" >> for005
echo "   154.2319398     -0.8784338824       227.9718659    " >> for005
echo "  0.3333333300      8.0000000000E-02" >> for005
echo "  0.5660748429       0.000000000      0.4514299287    " >> for005
echo "   0.000000000       227.9718659       0.000000000    " >> for005
echo "   0.000000000      0.3333333333    " >> for005
echo "  -7.271470384       4.862385245       8.528948242      0.6666666700      0.1600000000    " >> for005
echo "   0.000000000      0.3000000000    0 0 2" >> for005
echo "   0.000000000       0.000000000    " >> for005
echo "  0.1977863261     -1.6548496597E-02   0.000000000       0.000000000    " >> for005
echo "  0.00000000000,1.0,${j},   0.0000000000       0.0000000000    ,0,,  -30.00,${hfbwidth},, POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
echo "  50 70  0 0 ${i}" >> for005
echo "0         IFGSC" >> for005
echo "-1,,    IWUNIT" >> for005


../hafom_huge.exe < for005 > /dev/null

echo "cutoff=-${i},hfbwidth=${hfbwidth},pair.strength=${pairStrength}" >> collect_cat_FyDr_120Sn_cut30_wid1_${hfbwidth}_${cond}_shell_Cut100.dat
grep Z,N for011 >>  collect_cat_FyDr_120Sn_cut30_wid1_${hfbwidth}_${cond}_shell_Cut100.dat

done
echo "${cond} Done!"
mv collect_cat_FyDr_120Sn_cut30_wid1_${hfbwidth}_${cond}_shell_Cut100.dat FyDrTests/
#done <<< "${combined}"
done