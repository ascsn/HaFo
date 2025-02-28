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

nucleus="120Sn"
edf="FyDr"

ipair=(7 8)    # 7 for BCS, 8 for HFB

# Set hfbwidth
hfbwidth="0.0d0"

for j in "${ipair[@]}"; do
if [ $j -eq 7 ]; then
    cond="bcs"
else
    cond="hfb"
fi
echo "cond = ${j}, hfbwidth = ${hfbwidth}, pairStrength = ${pairStrength}" 
for i in $(seq 2 1 10); do
echo "1,,    IWUNIT" > for005
echo "'1D: FyDrHFB_cutoff15_wid0'" >> for005
echo "0,0.3,              NGRID,RSTEP" >> for005
echo "-1,1E-8,1E-9,278,5780,3,-1, IWR,DABR/E,ITPR/MAX,IFPL" >> for005
echo "0.2,90.,0.2,-1,,               X0DAMP,E0DAMP,ADDNEW,WALL" >> for005
echo " 1, 0, 0,-1, 1, 3, 1, 0, 0,13, if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk" >> for005
echo "1,,,                       ifform" >> for005
echo "0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12" >> for005
echo "0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12" >> for005
echo "  -9.786370533      0.6170974498      0.2123324350     avp,h1vp,h2vp" >> for005
echo "   25.25410957     -0.3578875136       25.94021760    " >> for005
echo "  0.3333333300      8.0000000000E-02" >> for005
echo "  0.5613396367       0.000000000      0.4355504811    " >> for005
echo "   0.000000000       25.94021760       0.000000000    " >> for005
echo "   0.000000000      0.3333333333    " >> for005
echo "  -4.452112339       4.241248563       3.034027663      0.6666666700      0.1600000000    " >> for005
echo "   0.000000000      0.3000000000    0 0 2" >> for005
echo "   0.000000000       0.000000000    " >> for005
echo "  0.1967729788      5.0471294254E-02   0.000000000       0.000000000    " >> for005
echo "  0.00000000000,1.0,${j},   0.0000000000       0.0000000000    ,0,,  -15.00,${hfbwidth},, POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
echo "  50 70  0 0 ${i}" >> for005
echo "0         IFGSC" >> for005
echo "-1,,    IWUNIT" >> for005


../../hafom_huge.exe < for005 > /dev/null

fileName="collect_cat_${edf}_${nucleus}_cut15_wid0_${hfbwidth}_${cond}.dat"
echo "cutoff=-${i},hfbwidth=${hfbwidth},shells=${i}" >> ${fileName}
grep Z,N for011 >> ${fileName}

done
echo "${cond} Done!"
mv ${fileName} collect_cat/
done