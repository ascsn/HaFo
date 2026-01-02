echo "1,,    IWUNIT" > for005.v2
echo "'SV-bas HFB'" >> for005.v2
echo "$4,0.3,              NGRID,RSTEP" >> for005.v2
echo "-1,1E-8,1E-9,200,20000,3,-1, IWR,DABR/E,ITPR/MAX,IFPL" >> for005.v2
echo "0.3,100.,0.2,-1,,,,,               X0DAMP,E0DAMP,ADDNEW,WALL" >> for005.v2
#echo "0.05,100.,0.1,-1,,,,,               X0DAMP,E0DAMP,ADDNEW,WALL" >> for005.v2
echo "1,0,0,1,-1,3,0,0,0,10,,, if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk" >> for005.v2
echo "1,,,                       ifform" >> for005.v2
echo "0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12" >> for005.v2
echo "0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12" >> for005.v2
echo " -1879.639,   313.749,   112.677, 12527.376,   124.634,      t0,t1,t2,t3,t4" >> for005.v2
echo "0.25854594    ,-.38168900    ,-2.8236379    ,0.12322898    ,   34.11167,       x0,x1,x2,x3,b4p" >> for005.v2
echo "   0.00000   0.20113   0.00000,,,      riftm,rho0pr,powerp" >> for005.v2
echo "     0.000,  0.000000,  0.000000,           tyuk,xyuk,ayuk" >> for005.v2
#echo "  0.300000,1.0,8,67.4618,1500,0,0.0,-15,1D0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005.v2
#echo "  0.300000,1.0,8,67.4618,323.1,0,0.0,-100,0D0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005.v2
echo "  0.300000,1.0,8,67.4618,$2,0,0.0,-25,$1,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005.v2
#echo "  0.300000,1.0,7,67.4618,470,0,0.0,-15,0D0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005.v2
#echo "  0.300000,1.0,8,67.4618,3800,0,0.0,-100,1.5D0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005.v2
#echo "  0.300000,1.0,8,67.4618,34500,0,0.0,-100,3D0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005.v2
#echo "  0.300000,1.0,8,67.4618,32000,0,0.0,-100,3D0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005.v2
echo "50,70,,,2,,,,,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005.v2
#echo "50,70,,,3,,,,,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005.v2
#echo "50,70,,,6,,,,,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005.v2
echo " 0         IFGSC" >> for005.v2
#echo " 77,$1" >> for005.v2
#echo " 43,$3" >> for005.v2
echo " 36,$3" >> for005.v2
#echo " 50,$1" >> for005.v2
echo "-1,,    IWUNIT" >> for005.v2

../hafotest.exe < for005.v2 > /dev/null
#../hafotest.new < for005.v2 > /dev/null
#../hafotest.new_nowall < for005.v2 > /dev/null
#cat fort.654 >> collect.v2-space2
#cat fort.654 >> addv2
cat fort.654 >> collect.v2-space2


