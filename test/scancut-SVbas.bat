#!/bin/bash
# Batch file for systematic variation of cutoff and folding width.
# Test case is 120Sn with 51 radial grid points and 17 shell above
# the Fermi energy. Proton pairing is practically switched off.
# The input parameters are:
#  $1   = cutoff energy 'cutfac'
#  $2   = folding width 'hfbwidth'
#  $3   = neutron pairing strength

# Set HFB widths
hfbwidths=("0.0d0" "1.0d0" "2.0d0" "3.0d0")

# Set BCS pairing strengths corresponding to those hfbwidths
# fixing lower cutoff pairing gap
pairStrengthsBCS=("577" "1884" "18700" "43510")
# fixing higher cutoff pairing gap


# Set HFB pairing strengths corresponding to those hfbwidths
pairStrengthsHFB=("529" "1725" "13460" "35064")

#ipair=(7)    # for BCS
ipair=(7 8)    # for HFB

# Set hfbwidth
hfbwidth="0.0d0"

for j in "${ipair[@]}"; do
if [ $j -eq 7 ]; then
    cond="bcs"
    #pairStrength="13000"
    combined=$(paste -d "," <(printf "%s\n" "${hfbwidths[@]}") <(printf "%s\n" "${pairStrengthsBCS[@]}"))
else
    cond="hfb"
#    pairStrength="529"
    combined=$(paste -d "," <(printf "%s\n" "${hfbwidths[@]}") <(printf "%s\n" "${pairStrengthsHFB[@]}"))
fi
while IFS=, read -r hfbwidth pairStrength; do
echo "cond = ${cond}, hfbwidth=${hfbwidth}, pairStrength = ${pairStrength}" 
for i in $(seq 15 1 120); do
echo "1,,    IWUNIT" > for005
echo "'SV-bas HFB'" >> for005
echo "71,0.3,              NGRID,RSTEP" >> for005
#echo "51,0.3,              NGRID,RSTEP" >> for005
echo "-1,1E-8,1E-9,100,3780,3,-1, IWR,DABR/E,ITPR/MAX,IFPL" >> for005
echo "0.2,100.,0.2,-1,,,,,               X0DAMP,E0DAMP,ADDNEW,WALL" >> for005
echo "1,0,0,1,-1,3,0,0,0,10,,, if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk" >> for005
echo "1,,,                       ifform" >> for005
echo "0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12" >> for005
echo "0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12" >> for005
echo " -1879.639,   313.749,   112.677, 12527.376,   124.634,      t0,t1,t2,t3,t4" >> for005
echo "0.25854594    ,-.38168900    ,-2.8236379    ,0.12322898    ,   34.11167,       x0,x1,x2,x3,b4p" >> for005
echo "   0.00000   0.20113   0.00000,,,      riftm,rho0pr,powerp" >> for005
echo "     0.000,  0.000000,  0.000000,           tyuk,xyuk,ayuk" >> for005
echo "  0.300000,1.0,${j},67.4618,${pairStrength},0,0.0,-${i},${hfbwidth},,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
echo "50,70,,,17,,,,,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
echo " 0         IFGSC" >> for005
echo "-1,,    IWUNIT" >> for005


../hafom_huge.exe < for005 > /dev/null

echo "cutoff=-${i},hfbwidth=${hfbwidth},pair.strength=${pairStrength}" >> collect_cat_SVbas_${hfbwidth}_${cond}_Cut100.dat
grep Z,N for011 >>  collect_cat_SVbas_${hfbwidth}_${cond}_Cut100.dat

done
echo "${cond} Done!"
mv collect_cat_SVbas_${hfbwidth}_${cond}_Cut100.dat collect_cat/
done <<< "${combined}"
done