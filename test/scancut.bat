#!/bin/bash
# Batch file for systematic variation of cutoff and folding width.
# Test case is 120Sn with 51 radial grid points and 17 shell above
# the Fermi energy. Proton pairing is practically switched off.
# The input parameters are:
#  $1   = cutoff energy 'cutfac'
#  $2   = folding width 'hfbwidth'
#  $3   = neutron pairing strength
# hfbwidth=(0d0 0.5d0 1d0 1.5d0 2d0)

ipair=(5)

# Set hfbwidth
hfbwidth="1.9d0"

# for j in "${hfbwidth[@]}"; do
for j in "${ipair[@]}"; do
if [ $j -eq 5 ]; then
    cond="bcs"
    pairStrength="9261"
else
    cond="hfb"
    pairStrength="1835"
fi
echo "cond = ${j}, pairStrength = ${pairStrength}" 
for i in $(seq 15.0 1. 120.0); do
echo " 1,                   IWUNIT " > for005
echo " 'UNEDF1' " >> for005
echo " 51,0.3,,,,                  NGRID,RSTEP" >> for005
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
echo " 0.2700180115,1.0,${j},25.0,${pairStrength},0,0.0,-${i},${hfbwidth},,,,,power,acoul,ipair,pairfp,pairfn,ilnrsp,temprt,cutfac,hfbwidth" >> for005
echo " 50,70,,,17,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
echo " 0,     ifgsc" >> for005
echo " -1,                   IWUNIT" >> for005

../hafom_huge.exe < for005 > /dev/null

echo "cutoff=-${i},hfbwidth=${hfbwidth},pair.strength=${pairStrength}" >> collect_cat_${hfbwidth}_${cond}.dat
grep Z,N for011 >>  collect_cat_${hfbwidth}_${cond}.dat

done
echo "${cond} Done!"
done