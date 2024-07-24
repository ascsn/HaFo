#!/bin/bash
# Batch file for systematic variation of cutoff and folding width.
# Test case is 146Gd radial grid points and 17 shell above
# the Fermi energy. Neutron pairing is practically switched off.
# 
# The oproton pairing strength is entered in the variable 'pairStrength'.
# The following values should be taken:
#  width   strength for BCS     strength for HFB
#  0.0          253                    253
#  0.5          350                    355
#  1.0          855                    893
#  1.5         2980                   3080
#  2.0        11300                  10850

#  0.0          237.8                 186      -- cutoff=100 pairing gap=1.835/1.838
#  0.5          332.25                 270.7
#  1.0          829                    798
#  1.5         2935                   2882
#  2.0        11105                  10015 

ipair=(5)

# Set hfbwidth
hfbwidth="1.5d0"
#hfbwidth="3.0d0"
#hfbwidth="4.0d0"
hfbwidth="5.0d0"

# for j in "${hfbwidth[@]}"; do
for j in "${ipair[@]}"; do
if [ $j -eq 5 ]; then
    cond="bcs"
#    pairStrength="46950"
#    pairStrength="124800"
    pairStrength="306500"
else
    cond="hfb"
    pairStrength="2882"
fi
echo "cond = ${j}, pairStrength = ${pairStrength}" 
for i in $(seq 99 1 100); do
echo " 1,                   IWUNIT " >> for005
echo " 'Fy(std), BCS mixed pairing' " >> for005
echo " 71,0.3,,,,                  NGRID,RSTEP" >> for005
echo "-1,1E-9,1E-9,300,2000,3,-1,,, IWR,DABR/E,ITPR/MAX,IFPL" >> for005
echo "0.2,90.,0.2,-1,,,               X0DAMP,E0DAMP,ADDNEW,WALL" >> for005
echo " 1, 0, 0, 1, 1, 3, 0, 0, 0,13,,, if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk" >> for005
echo "1,,,                       ifform" >> for005
echo "0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12" >> for005
echo "0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12" >> for005
echo "  -9.495922     0.6271056      0.1452334   avp,h1vp,h2vp     " >> for005
echo "   12.47412     -1.387537       19.37946   avm,h1vm,h2vm" >> for005
echo "  0.3333333     8.0000000E-02              sigmavol,rho0fy" >> for005
echo "  0.5241615      0.000000  9.9219924E-02   asp,hsp,hgradp " >> for005
echo "   0.000000     -2000.000      0.000000    asm,hsm,hgradm " >> for005
echo "  0.000000      0.3333333                  hcoul,sigmacoul" >> for005
echo " -1.6398493       1.1311550      1.02627800E-02  0.66666670  0.16 fxi,h1xi,hgradxi,gammaxi,rh0xi" >> for005
echo "   0.0  0.0,5,,,,  1     delfxi ecstab isovectormodel" >> for005
echo "   0.000000     0.000000                   c0tau,c1tau " >> for005
echo "  0.1896399   2.4966800E-02  0.0  0.0      c0nabj,c1nabj,c0tens,c1tens" >> for005
echo "  0,1.0,${j},${pairStrength},28,0,0,-${i},${hfbwidth},,,,,       power,acoul,ipair,pairfp,pairfn,ilnrsp,temprt,cutfac" >> for005
echo " 50,70,,,17,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
echo " 0,     ifgsc" >> for005
echo " -1,                   IWUNIT" >> for005

../hafom_huge.exe < for005 > /dev/null

echo "cutoff=-${i},hfbwidth=${hfbwidth},pair.strength=${pairStrength}" >> collect_146Gd_${hfbwidth}_${cond}_Fy.dat
grep Z,N for011 >>  collect_146Gd_${hfbwidth}_${cond}_Fy.dat
echo "Cutoff = ${i}"
done
echo "${cond} Done!"
done
