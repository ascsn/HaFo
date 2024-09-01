PROGRAM prepNMPsets
IMPLICIT NONE
!
! Routine to produce input for the 1D SHF code from NMP.
! NMP are read from standard input.
!

REAL(8), PARAMETER :: h2m= (20.74982D0+20.72126D0)/2D0
REAL(8), PARAMETER :: rhonm_trans = 0.160

REAL(8) :: EA,rhoeq,mstar,K
REAL(8) :: b0,b1,b3,power
REAL(8) :: c0,c3,powfac
INTEGER :: ik

  READ(5,*) EA,rhoeq,mstar,K
  CALL params_from_NMP(power,B0,B3,B1,rhoeq,EA,h2m,mstar,K)
  B3=6D0*B3

! recouple to C parameters
  powfac = rhonm_trans**power/6.0D0
  c0     = b0+powfac*b3
  c3     = b0-powfac*b3
  ik=NINT(K)

  OPEN(1,file='for001')

  WRITE(1,'(a,i3,a,f4.2,a)')  "011,'SV-K",ik,"-m",mstar,"'"
  WRITE(1,'(a)')&
  "10,0,,,,      NWDR",&
  "  1,    2,   3,   4,  5,  6,  11,  9, 10,  12,   8,0,0,0,0,0,  MESSGROESSEN",&
  "1.0,0.040,0.04,0.02,20.,20.,0.36,0.12,0.12,0.01,90.1,0.,0.,0., FEHLER",&
  "1.0,  1.0, 1.0, 1.0,0.0,0.0, 1.0,1.0,1.0, 1.0, 1.0,0.,0.,0.,0,0,0, % O. ABSOLUT",&
  "00.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., KORRREKTUR",&
  "016,1000,0,         IAUNT,IAOB,",&
  "0,1,0,1,1,1,0,3,0,,,,  NMINV,IF:EX,TM,ENRG,RSO,ZPE,GSC,TOND,MFOLD",&
  "1.0,0.2,.1E-8,.1E-7,7,0,2,,,,   ACOUL,ADDNEW,DABE,DABR,ipair,ilnrsp,iparsh",&
  "0.2,20,0.2,20,            RSYMAX,NSYMAX,RNTMAX,NNTMAX",&
  "0,0,4,2820,0,,                IWRTE,IST1,IFCURV,ITMAX,IEXTAN",&
  "38,09,0,2.0,0.0001,2,,11,, NPARIN,NUMAX,IVAR,FREEP,ENDCON,IFTINP"
  WRITE(1,'(1pg14.6,a)')&
  c0,"  , 0.40247E-01, 1,    c0",&
  c3,"  ,  1.4145    , 4,    c3",&
  100D0,"  , 0.20374    , 3,    b2",&
  -300D0,"  ,  1.5433    , 6,   c0p",&
  -30D0,"   , 0.44751    , 7,   b1p",&
  150D0,"   , 0.87043    , 5,    t4",&
  0D0,"   ,  15.223    , 8,   b2p",&
  -6000D0,"  ,  57.451    , 9,   c3p",&
  30D0,"  ,  22.151    ,15,   b4p",&
  600, "  ,  26.662    ,16,pairfp",&
  600, "  ,  40.103    ,17,pairfn",&
  0.2D0,"  , 0.48553E-02,18,rho0pr",&
  power,"  ,  0.0010    ,10, power",&
  b1,"  , 0.95915    , 2,    b1"
  WRITE(1,'(a)')&
  "    1D-4    ,   0.7    ,38,str_sq",&
  "  0.00000    , 0.20000    ,14, riftm",&
  "  0.00000    ,  0.0000    ,11,  byuk",&
  "  0.00000    ,  0.0000    ,12, byukp",&
  "  0.00000    ,  0.0000    ,13,  ayuk",&
  "  0.00000    ,  0.0000    ,19,powerp",&
  "  0.00000    ,  0.0000    ,20,b3tild",&
  "  0.00000    ,  0.0000    ,21,   b12",&
  "  0.00000    ,  0.0000    ,22,  b12p",&
  "  0.00000    ,  0.0000    ,23,    b7",&
  "  0.00000    ,  0.0000    ,24,   b7p",&
  "  0.00000    ,  0.0000    ,25,  pow7",&
  "  0.00000    ,  0.0000    ,26,    c1",&
  "  0.00000    ,  0.0000    ,27,   c1p",&
  "  0.00000    ,  0.0000    ,28,    c4",&
  "  0.00000    ,  0.0000    ,29,   c4p",&
  "  0.00000    ,  0.0000    ,30,ecstab",&
  "  0.00000    ,  0.0000    ,31,alpinv",&
  "  0.00000    ,  0.0000    ,32,betinv",&
  "  0.00000    ,  0.0000    ,33,alpin2",&
  "  0.00000    ,  0.0000    ,34,betin2",&
  "  0.00000    ,  0.0000    ,35,riftm1",&
  "  0.00000    ,  0.0000    ,36,   bl2",&
  "  0.00000    ,  0.0000    ,37,  bl2p",&
  "-11, 'stop'"

CLOSE(1)
STOP

CONTAINS
SUBROUTINE params_from_NMP(power_from_K,B0_deduced,B3_deduced,&
                           B1_deduced,&
                           rhotot,eovera,hbar2m,effmas,compr)
IMPLICIT NONE
!
! computes isoscalar Skyrme parameters from NMP.
! Input parameters (NMP) are:
!  rhotot = equilibrium density of symmetric matter,
!  eovera = binding energy per nucleon at equilibrium,
!  hbar2m = hbar**2/(2*nucleon_mass),
!  effmas = effective mass m*/m
!  compr  = incompressibility at equilibrium.
! Output parameters are:
!  power_from_K = power of density dependence,
!  B0_deduced   = B0 parameter,
!  B1_deduced   = B1 parameter,
!  B3_deduced   = B3 parameter.
!
!  Note: The B3 parameter deduced here is factor 6 smaller than
!        the equivalent used in the 1D code.
!
  REAL(8), INTENT(IN) :: rhotot,eovera,hbar2m,effmas,compr
  REAL(8), INTENT(OUT) :: power_from_K,B0_deduced,B3_deduced,B1_deduced

! (3/5)*(3*pi^2)^*(2/3)  factor for tau_n/p
  REAL(8), PARAMETER :: f2d3= 2D0/3D0, f5d3=5D0/3D0 
  REAL(8), PARAMETER :: betab=5.742467999D0 
  REAL(8), PARAMETER :: betat=betab*0.5D0**f2d3
  REAL(8) :: counter,denom,b1,power
  REAL(8) :: hbmterm      ! kinetic energy term
  REAL(8) :: b1term       ! kinetic interaction term
  REAL(8) :: rhopowp1     ! rho**(1+power)/2

  b1       = hbar2m*(1D0/effmas-1D0)/rhotot
  B1_deduced = b1

  hbmterm  = hbar2m*betat*rhotot**f2d3
  b1term   = b1*betat*rhotot**f5d3

  counter = 10D0*b1term-2D0*hbmterm-compr
  denom   = 9D0*eovera+6D0*b1term-3D0*hbmterm
  power   = counter/denom-1D0
  power_from_K = power


  rhopowp1 = rhotot**(1D0+power)/2D0
  B3_deduced = (hbmterm/3D0-f2d3*b1term-eovera)/(power*rhopowp1)
  B0_deduced = 2D0*(eovera-hbmterm-b1term-B3_deduced*rhopowp1)/rhotot

  RETURN
END SUBROUTINE params_from_NMP

END PROGRAM prepNMPsets
