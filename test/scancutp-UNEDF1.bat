#!/bin/bash
# Batch file for systematic variation of cutoff and folding width.
# Test case is 120Sn with 51 radial grid points and 17 shell above
# the Fermi energy. Proton pairing is practically switched off.
# The input parameters are:
#  $1   = cutoff energy 'cutfac'
#  $2   = folding width 'hfbwidth'
#  $3   = neutron pairing strength

nucleus="146Gd"
edf="UNEDF1"

ipair=(5) # For BCS
#ipair=(6) # For HFB

# Set hfbwidths
#hfbwidths=("0.0d0" "1.0d0" "1.5d0" "2.0d0" "3.0d0")

hfbwidth="3.0d0"

pairStrengthsBCS=("253" "858" "2981" "11150" "81650")
pairStrengthsHFB=("215" "860" "3070" "10000" "25000")

for j in "${ipair[@]}"; do
if [ $j -eq 5 ]; then
    cond="bcs"
    pairStrength="81650"
#    combined=$(paste -d "," <(printf "%s\n" "${hfbwidths[@]}") <(printf "%s\n" "${pairStrengthsBCS[@]}"))
else
    cond="hfb"
    pairStrength="1835"
#    combined=$(paste -d "," <(printf "%s\n" "${hfbwidths[@]}") <(printf "%s\n" "${pairStrengthsHFB[@]}"))
fi
#while IFS=, read -r hfbwidth pairStrength; do
echo "cond = ${j}, hfbwidth=${hfbwidth}, pairStrength = ${pairStrength}" 
for i in $(seq 15 1 120); do
echo " 1,                   IWUNIT " > for005
echo " 'UNEDF1' " >> for005
echo " 71,0.3,,,,                  NGRID,RSTEP" >> for005
echo "-1,1.E-8,1.E-9,278,9000,3,-1,,   IWR,DABR,DABE,ITPR,ITMAX,IRAD,IFPLOT,IDISC " >> for005
echo " 0.4,100.,0.2,-1.0,,,              XODAMP,E0DAMP,ADDNEW,WALL" >> for005
echo " 1,0,0,1,-2,2,0,0,,,,,      IF:EX,TM,ENRG,RSO,ZPE,DFMS,RHOC,ENAV" >> for005
echo " 1,,                      IFFORM" >> for005
echo " 0.64,0.0,0.0,0.0,0.0,    APRO0,APR011,APRO12,BPRO11,BPRO12" >> for005
echo " 0.64,0.0,0.0,0.0,0.0,    ANEU0,ANEU11,ANEU12,BNEU11,BNEU12" >> for005
echo "-2078.32802325563625, 239.400812041522045, 1574.24304337273043, 14263.6462470775750, 76.7361441233364018" >> for005
echo " 0.543043709322549584E-01, -5.07814768633786962, -1.36657048130840808,-0.161168850358621302,  71.3165222295833985," >> for005
echo " 0.000,0.000,0.000,0.000,0.000,    t5,t6,t7,t8,t9" >> for005
echo " 0.000,0.000,0.000,0.000000,       x5,x6,x7,x8" >> for005
echo " 0.000,0.000,0.000,                tyuk,xyuk,ayuk" >> for005
echo " 0.2700180115,1.0,${j},${pairStrength},258.11805,0,0.0,-${i},${hfbwidth},,,,,power,acoul,ipair,pairfp,pairfn,ilnrsp,temprt,cutfac,hfbwidth" >> for005
echo " 64,82,,,17,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
echo " 0,     ifgsc" >> for005
echo " -1,                   IWUNIT" >> for005

../hafom_huge.exe < for005 > /dev/null

fileName="collect_cat_${edf}_${nucleus}_${hfbwidth}_${cond}.dat"
echo "cutoff=-${i},hfbwidth=${hfbwidth},pair.strength=${pairStrength}" >> ${fileName}
grep Z,N for011 >>  ${fileName}

done
echo "${cond} Done!"
mv ${fileName} collect_cat/
#done <<< "${combined}"
done