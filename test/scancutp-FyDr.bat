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

nucleus="146Gd"
edf="FyDr"

# Set HFB widths
hfbwidths=("0.0d0" "1.0d0" "1.5d0" "2.0d0" "3.0d0")

#ipair=(7)    # for BCS
ipair=(8)    # for HFB

# Set hfbwidth
#hfbwidth="3.0d0"

#pairStrengthsBCS=("4.485" "6.548" "14.1" "43.4" "137.9")
pairStrengthsBCS=("3.072" "6.615" "13.55" "42.65" "150.05") # Cut100
#pairStrengthsHFB=("4.358" "6.375" "13.62" "37.3" "114.4")
pairStrengthsHFB=("3.349" "6.21" "12.83" "36.2" "143.35") # Cut100

for j in "${ipair[@]}"; do
if [ $j -eq 7 ]; then
    cond="bcs"
#    pairStrength="3.072"
    combined=$(paste -d "," <(printf "%s\n" "${hfbwidths[@]}") <(printf "%s\n" "${pairStrengthsBCS[@]}"))
else
    cond="hfb"
#    pairStrength="143.35"
    combined=$(paste -d "," <(printf "%s\n" "${hfbwidths[@]}") <(printf "%s\n" "${pairStrengthsHFB[@]}"))
fi
while IFS=, read -r hfbwidth pairStrength; do
echo "cond = ${j}, pairStrength = ${pairStrength}" 
for i in $(seq 15 1 100); do
echo "1,,    IWUNIT" > for005
echo " 'Fy(Dr)' " >> for005
echo " 71,0.3,,,,                  NGRID,RSTEP" >> for005
echo "-1,1E-8,1E-9,278,9000,3,-1,,, IWR,DABR/E,ITPR/MAX,IFPL" >> for005
echo "0.3,110.,0.2,-1,,,               X0DAMP,E0DAMP,ADDNEW,WALL" >> for005
echo " 1, 0, 0, 1, 1, 3, 0, 0, 0,13,,, if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk" >> for005
echo "1,,,                       ifform" >> for005
echo "0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12" >> for005
echo "0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12" >> for005
echo "  -9.542989   0.6323225   0.1343786  avp,h1vp,h2vp     " >> for005
echo "   4.182360   0.2537760   1.215020  avm,h1vm,h2vm" >> for005
echo "  0.3333333   8.0000000E-02  sigmavol,rho0fy" >> for005
echo "  0.6047266    0.000000  0.6656000  asp,hsp,hgradp " >> for005
echo "   0.000000   -2000.0     0.000000  asm,hsm,hgradm " >> for005
echo "   0.000000    0.3333333  hcoul,sigmacoul" >> for005
#echo "  -4.472000    4.229000    3.227000  0.6666667   0.1600000  fxi,h1xi,hgradxi,gammaxi,rho0xi" >> for005
echo "  -${pairStrength}    4.229000    3.227000  0.6666667   0.1600000  fxi,h1xi,hgradxi,gammaxi,rho0xi" >> for005
echo "   0.0  0.0  0 0 1,,,     delfxi ecstab isovectormodel" >> for005
echo "   0.000000    0.000000     c0tau,c1tau " >> for005
echo "   0.1879220  -2.3700000E-02 0.000000  0.000000,,,    c0nabj,c1nabj,c0tens,c1tens" >> for005
echo "  0,1.0,${j},700,700,0,0,-${i},${hfbwidth},,,       power,acoul,ipair,pairfp,pairfn,ilnrsp,temprt,cutfac" >> for005
echo "64,82,,,17,,,,,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
echo " 0         IFGSC" >> for005
echo "-1,,    IWUNIT" >> for005


../hafom_huge.exe < for005 > /dev/null

fileName="collect_cat_${edf}_${nucleus}_${hfbwidth}_${cond}_Cut100.dat"
echo "cutoff=-${i},hfbwidth=${hfbwidth},pair.strength=${pairStrength}" >> ${fileName}
grep Z,N for011 >> ${fileName}

done
echo "${cond} Done!"
mv ${fileName} collect_cat/
done <<< "${combined}"
done
