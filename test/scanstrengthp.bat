#echo "1,,    IWUNIT" > for005
#echo "'SV-bas HFB'" >> for005
#echo "71,0.3,              NGRID,RSTEP" >> for005
#echo "-1,1E-8,1E-9,100,3780,3,-1, IWR,DABR/E,ITPR/MAX,IFPL" >> for005
#echo "0.2,100.,0.2,-1,,,,,               X0DAMP,E0DAMP,ADDNEW,WALL" >> for005
#echo "1,0,0,1,-1,3,0,0,0,10,,, if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk" >> for005
#echo "1,,,                       ifform" >> for005
#echo "0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12" >> for005
#echo "0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12" >> for005
#echo " -1879.639,   313.749,   112.677, 12527.376,   124.634,      t0,t1,t2,t3,t4" >> for005
#echo "0.25854594    ,-.38168900    ,-2.8236379    ,0.12322898    ,   34.11167,       x0,x1,x2,x3,b4p" >> for005
#echo "   0.00000   0.20113   0.00000,,,      riftm,rho0pr,powerp" >> for005
#echo "     0.000,  0.000000,  0.000000,           tyuk,xyuk,ayuk" >> for005
#echo "  0.300000,1.0,8,67.4618,$1,0,0.0,-15,1D0,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,..." >> for005
#echo "50,74,,,3,,,,,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
#echo " 0         IFGSC" >> for005
#echo "-1,,    IWUNIT" >> for005

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
echo " 0.270018011502708211,1.0,$3,$1,25.8,0,0.0,-15,$2,,,,,power,acoul,ipair,pairfp,pairfn,ilnrsp,temprt" >> for005
echo " 64,82,,,3,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
echo " 0,     ifgsc" >> for005
#echo " 1,                   IWUNIT " >> for005
#echo " 'UNEDF1' " >> for005
#echo " 71,0.3,,,,                  NGRID,RSTEP" >> for005
#echo "-1,1.E-8,1.E-9,278,9000,3,-1,,   IWR,DABR,DABE,ITPR,ITMAX,IRAD,IFPLOT,IDISC " >> for005
#echo " 0.4,100.,0.2,-1.0,,,              XODAMP,E0DAMP,ADDNEW,WALL" >> for005
#echo " 1,0,0,1,-2,2,0,0,,,,,      IF:EX,TM,ENRG,RSO,ZPE,DFMS,RHOC,ENAV" >> for005
#echo " 1,,                      IFFORM" >> for005
#echo " 0.64,0.0,0.0,0.0,0.0,    APRO0,APR011,APRO12,BPRO11,BPRO12" >> for005
#echo " 0.64,0.0,0.0,0.0,0.0,    ANEU0,ANEU11,ANEU12,BNEU11,BNEU12" >> for005
#echo "-2078.32802325563625, 239.400812041522045, 1574.24304337273043, 14263.6462470775750, 76.7361441233364018" >> for005
#echo " 0.543043709322549584E-01, -5.07814768633786962, -1.36657048130840808,-0.161168850358621302,  71.3165222295833985," >> for005
#echo " 0.000,0.000,0.000,0.000,0.000,    t5,t6,t7,t8,t9" >> for005
#echo " 0.000,0.000,0.000,0.000000,       x5,x6,x7,x8" >> for005
#echo " 0.000,0.000,0.000,                tyuk,xyuk,ayuk" >> for005
#echo " 0.270018011502708211,1.0,5,$1,25.8,0,0.0,-15,$2,,,,,power,acoul,ipair,pairfp,pairfn,ilnrsp,temprt" >> for005
#echo " 68,82,,,3,,  nprot,nneut,nmaxp,nmaxn,iparsh" >> for005
#echo " 0,     ifgsc" >> for005
echo " -1,                   IWUNIT" >> for005


../haforpa.exe < for005 > /dev/null
grep Z,N for011
