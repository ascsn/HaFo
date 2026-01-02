#!/bin/bash
# Batch file for systematic variation of cutoff and folding width.
# Test case is 120Sn with 51 radial grid points and 17 shell above
# the Fermi energy. Proton pairing is practically switched off.
# The input parameters are:
#  $1   = cutoff energy 'cutfac'
#  $2   = folding width 'hfbwidth'
#  $3   = neutron pairing strength

nucleus="120Sn"
edf="SVbas"

# Test box size, width = 0.0 for pairgap = 1.61 at cutoff = 25 MeV
#pairStrengthsBCS=("537")  

# Test box size, width = 0.0 for pairgap = 1.61 at cutoff = 25 MeV
#pairStrengthsHFB=("456.9")  

# Test box size, width = 1.0 for pairgap = 1.61 at cutoff = 25 MeV
#pairStrengthsHFB=("385.3")  

# Test box size, width = 1.5 for pairgap = 1.61 at cutoff = 25 MeV
#pairStrengthsHFB=("4531")  

# Test box size, width = 2.0 for pairgap = 1.61 at cutoff = 25 MeV
#pairStrengthsHFB=("10993")  

# Test box size, width = 3.0 for pairgap = 1.61 at cutoff = 25 MeV
#pairStrengthsHFB=("30710")  

# Test box size, width = 0.0 for pairgap = 1.61 at cutoff = 100 MeV
#pairStrengthsHFB=("456.9")  

# Test box size, width = 1.0 for pairgap = 1.61 at cutoff = 100 MeV
#pairStrengthsHFB=("1500.2")  

# Test box size, width = 1.5 for pairgap = 1.61 at cutoff = 100 MeV
#pairStrengthsHFB=("4539")  

# Test box size, width = 2.0 for pairgap = 1.61 at cutoff = 100 MeV
#pairStrengthsHFB=("10993")  

# Test box size, width = 3.0 for pairgap = 1.61 at cutoff = 100 MeV
#pairStrengthsHFB=("30300")  

#ipair=(7)    # for BCS
ipair=(8)    # for HFB

# Set hfbwidth
hfbwidth="3.0d0"

# Set grid size
ngrid="71"

for j in "${ipair[@]}"; do
if [ $j -eq 7 ]; then
    cond="bcs"
    pairStrength="537"
else
    cond="hfb"
    pairStrength="30235.012020"
fi
i=25
#while IFS=, read -r ngrid pairStrength; do
echo "cond = ${cond}, hfbwidth=${hfbwidth}, pairStrength = ${pairStrength}" 
for ngrid in $(seq 51 2 129); do
echo "1,,    IWUNIT" > for005
echo "'SV-bas HFB'" >> for005
echo "${ngrid},0.3,              NGRID,RSTEP" >> for005
#echo "51,0.3,              NGRID,RSTEP" >> for005
echo "-1,1E-8,1E-9,200,20000,3,-1, IWR,DABR/E,ITPR/MAX,IFPL" >> for005
echo "0.3,100.,0.2,-1,,,,,               X0DAMP,E0DAMP,ADDNEW,WALL" >> for005
echo "1,0,0,-1,1,3,0,0,0,10,,, if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk" >> for005
echo "1,,,                       ifform" >> for005
echo "0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12" >> for005
echo "0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12" >> for005
echo " -1879.639,   313.749,   112.677, 12527.376,   124.634,      t0,t1,t2,t3,t4" >> for005
echo "0.25854594    ,-.38168900    ,-2.8236379    ,0.12322898    ,   34.11167,       x0,x1,x2,x3,b4p" >> for005
echo "   0.00000   0.20113   0.00000,,,      riftm,rho0pr,powerp" >> for005
echo "     0.000,  0.000000,  0.000000,           tyuk,xyuk,ayuk" >> for005
echo "  0.300000,1.0,${j},10,${pairStrength},0,0.0,-${i},${hfbwidth},,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
echo "50,70,,,7,,,,,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
echo " 0         IFGSC" >> for005
echo "-1,,    IWUNIT" >> for005


../hafom_huge.exe < for005 > /dev/null

#fileName="collect_cat_${edf}_${nucleus}_${hfbwidth}_${cond}_Cut${i}_${ngrid}.dat"
fileName="collect_cat_${edf}_${nucleus}_${hfbwidth}_${cond}_Cut${i}_sh.dat"
echo "ngrid=${ngrid},hfbwidth=${hfbwidth},pair.strength=${pairStrength}" >> ${fileName}
grep Z,N for011 >> ${fileName}

done
echo "${fileName} Done!"
echo "${cond} Done!"
mv ${fileName} collect_cat_box_paper_new/
#done <<< "${combined}"
done