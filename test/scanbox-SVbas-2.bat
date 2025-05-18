rm -rf $1
mkdir $1
cd $1

echo "$1" >> ../collect_box_2.dat
cp ../hfbparams .

echo "1,,    IWUNIT" > for005
echo "'SV-bas HFB'" >> for005
echo "$1,0.3,              NGRID,RSTEP" >> for005
echo "-1,1E-8,1E-9,200,20000,3,-1, IWR,DABR/E,ITPR/MAX,IFPL" >> for005
echo "0.05,120.,0.1,-1,,,,,               X0DAMP,E0DAMP,ADDNEW,WALL" >> for005
#echo "0.2,100.,0.1,-1,,,,,               X0DAMP,E0DAMP,ADDNEW,WALL" >> for005
echo "1,0,0,-1,1,3,0,0,0,10,,, if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk" >> for005
echo "1,,,                       ifform" >> for005
echo "0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12" >> for005
echo "0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12" >> for005
echo " -1879.639,   313.749,   112.677, 12527.376,   124.634,      t0,t1,t2,t3,t4" >> for005
echo "0.25854594    ,-.38168900    ,-2.8236379    ,0.12322898    ,   34.11167,       x0,x1,x2,x3,b4p" >> for005
echo "   0.00000   0.20113   0.00000,,,      riftm,rho0pr,powerp" >> for005
echo "     0.000,  0.000000,  0.000000,           tyuk,xyuk,ayuk" >> for005
#echo "  0.300000,1.0,8,67.4618,4553,0,0.0,-25,1.5d0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
echo "  0.300000,1.0,8,67.4618,30710,0,0.0,-25,3d0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
#echo "  0.300000,1.0,8,67.4618,34000,0,0.0,-25,3d0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
#echo "  0.300000,1.0,7,10.0,40000,0,0.0,-25,3d0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
#echo "  0.300000,1.0,7,10.0,12000,0,0.0,-25,2d0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
#echo "  0.300000,1.0,7,10.0,24000,0,0.0,-25,2.5d0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
#echo "  0.300000,1.0,7,10.0,1600,0,0.0,-25,1d0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
#echo "50,70,,,23.7,,,,,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
#echo "50,70,,,21,,,,,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
echo "50,70,,,9,,,,,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
echo " 0         IFGSC" >> for005
echo "-1,,    IWUNIT" >> for005


../../hafom_huge.exe < for005 > /dev/null 

#grep Z,N  for011 >> ../collect_box.$1
grep Z,N  for011 >> ../collect_box_2.dat


cd ..
#rm -rf $1
