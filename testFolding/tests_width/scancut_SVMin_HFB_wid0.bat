#!/bin/bash
# Batch file for systematic variation of cutoff and folding width.
# Test case is 120Sn with 51 radial grid points and 17 shell above
# the Fermi energy. Proton pairing is practically switched off.

nucleus="120Sn"
edf="SVMin"

#ipair=(7)    # for BCS
ipair=(7 8)    # for HFB

# Set hfbwidth
hfbwidth="0.0d0"

for j in "${ipair[@]}"; do
if [ $j -eq 7 ]; then
    cond="bcs"
    #pairStrength="5984"
    #combined=$(paste -d "," <(printf "%s\n" "${hfbwidths[@]}") <(printf "%s\n" "${pairStrengthsBCS[@]}"))
else
    cond="hfb"
#    pairStrength="30290."
    #combined=$(paste -d "," <(printf "%s\n" "${hfbwidths[@]}") <(printf "%s\n" "${pairStrengthsHFB[@]}"))
fi
#while IFS=, read -r hfbwidth pairStrength; do
echo "cond = ${cond}, hfbwidth=${hfbwidth}, pairStrength = ${pairStrength}" 
for i in $(seq 2 1 20); do
echo "1,,    IWUNIT" > for005
echo "'SV-MIN HFB'" >> for005
echo "0,0.3,              NGRID,RSTEP" >> for005
#echo "51,0.3,              NGRID,RSTEP" >> for005
echo "-1,1E-8,1E-9,278,5780,3,-1, IWR,DABR/E,ITPR/MAX,IFPL" >> for005
echo "0.2,90.,0.2,-1,,,,,               X0DAMP,E0DAMP,ADDNEW,WALL" >> for005
echo "1,0,0,-1,1,3,0,0,0,11,,, if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk" >> for005
echo "1,,,                       ifform" >> for005
echo "0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12" >> for005
echo "0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12" >> for005
echo "  -1948.8999502,   308.02577052,  -69.158031953,   12958.468481,   120.06090398,      t0,t1,t2,t3,t4" >> for005
echo "  0.34467230040,  -1.0483321855,   1.3722367298,  0.48060206968,   37.146007709,       x0,x1,x2,x3,b4p" >> for005
echo "   0.0000000, 0.0000000,                      c1, c4 " >> for005
echo "   0.0000000,   0.0000000,                    c1p,c4p " >> for005
echo "   0.0000000,,,                      b3tilde " >> for005
echo "   0.0000000,   0.0000000,                     b12,b12p " >> for005
echo "   0.0000000,   0.0000000,               b7,b7p " >> for005
echo "   0.0000000,               pow7" >> for005
echo "   0.0000000000, 0.19652400010, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000,      riftm,rho0pr,powerp,ecstab,alpinv" >> for005
echo "     0.000,  0.000000,  0.000000,           tyuk,xyuk,ayuk" >> for005
echo "   0.28535398169,1.0,${j}, 488.57136963, 450.93449381,0,,  -30.00,${hfbwidth},, POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
echo "50,70,,,${i},,,,,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
echo " 0         IFGSC" >> for005
echo "-1,,    IWUNIT" >> for005


../../hafom_huge.exe < for005 > /dev/null

fileName="collect_cat_${edf}_${nucleus}_cut30_HFB_wid0_${hfbwidth}_${cond}.dat"
echo "cutoff=-${i},hfbwidth=${hfbwidth},shells=${i}" >> ${fileName}
grep Z,N for011 >> ${fileName}

done
echo "${cond} Done!"
mv ${fileName} collect_cat/
#done <<< "${combined}"
done