#!/bin/bash

echo "#SV-bas HFB cut25" >> collect_box.dat

#for i in $(seq 51 2 121)
for i in 51 53 55 57 59  61 63 65 67 69  71 73 75 77 79  81 83 85 87 89  91 93 95 97 99  101 103 105 107 109  111 113 115 117 119  121 123 125 127 129 
do


echo "$i" >> collect_box.dat

echo "1,,    IWUNIT" > for005
echo "'SV-bas HFB'" >> for005
echo "$i,0.3,              NGRID,RSTEP" >> for005
echo "-1,1E-8,1E-9,200,20000,3,-1, IWR,DABR/E,ITPR/MAX,IFPL" >> for005
echo "0.4,100.,0.2,-1,,,,,               X0DAMP,E0DAMP,ADDNEW,WALL" >> for005
echo "1,0,0,-1,1,3,0,0,0,10,,, if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk" >> for005
echo "1,,,                       ifform" >> for005
echo "0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12" >> for005
echo "0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12" >> for005
echo " -1879.639,   313.749,   112.677, 12527.376,   124.634,      t0,t1,t2,t3,t4" >> for005
echo "0.25854594    ,-.38168900    ,-2.8236379    ,0.12322898    ,   34.11167,       x0,x1,x2,x3,b4p" >> for005
echo "   0.00000   0.20113   0.00000,,,      riftm,rho0pr,powerp" >> for005
echo "     0.000,  0.000000,  0.000000,           tyuk,xyuk,ayuk" >> for005
echo "  0.300000,1.0,8,67.4618,4553,0,0.0,-25,1.5d0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
echo "50,70,,,17,,,,,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
echo " 0         IFGSC" >> for005
echo "-1,,    IWUNIT" >> for005


../hafom_huge.exe < for005 > /dev/null 

grep Z,N  for011 >> collect_box.dat

done
