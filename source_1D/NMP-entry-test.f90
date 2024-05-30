MODULE NMP_map
! Routine to determine the volume parameters of the Fayans
! functional from given nuclear matter parameters (NMP).
!
! Presently without kinetic parameters.
!
IMPLICIT NONE


REAL(8),PRIVATE, PARAMETER :: PI=3.14159265358979D0
REAL(8),PRIVATE, PARAMETER :: betap=0.6D0*(1.5D0*pi**2)**(2D0/3D0)
REAL(8),PRIVATE, PARAMETER :: beta=0.6D0*(3D0*pi**2)**(2D0/3D0)
REAL(8),PRIVATE, PARAMETER :: r23=2D0/3D0

REAL(8),PRIVATE :: bmass0,bmass1,r0fy,r0fy2,efermfy,c0tau,c1tau,h2vm
REAL(8),PRIVATE :: ap,am,hp,hm                        ! auxiliary Fayans parameters
REAL(8),PRIVATE :: xp,xps,h2var,dp,efacfy,rho0fy2
REAL(8),PRIVATE :: rhonm,eovera,compr,effmas,asymme,dasymm,skappa,&
           bmassp,bmassn,rho0fy,sigmavol
REAL(8),PRIVATE :: rhs,rhsp,asy(0:3)
LOGICAL,PRIVATE,PARAMETER :: ttest=.TRUE.

CONTAINS
SUBROUTINE NMP_Fayans(avp,h1vp,h2vp,avm,h1vm,h2vmin,sigmavolin,&
                      c0t,c1t,rho0fyin,&
                      rhonmin,eoverain,comprin,effmasin,asymmein,dasymmin,skappain,&
                      bmasspin,bmassnin)

! NMP as input, also input the fixed Fayans parameter 'rho0fy'
REAL(8), INTENT(IN) :: rhonmin,eoverain,comprin,effmasin,asymmein,dasymmin,&
                       skappain,bmasspin,bmassnin,rho0fyin,sigmavolin,h2vmin
! Resulting Fayans volume parameters
REAL(8), INTENT(OUT) :: avp,h1vp,h2vp,avm,h1vm,c0t,c1t

!REAL(8) :: compfunct
REAL(8) :: r0fy,r0fy2,efermfy,ftest
REAL(8):: energ0,eslope0,energ1, eslope1, energ2,eslope2,ekin0,ekin1,ekin2
INTEGER :: i

! Map NMP to internal variables
  rhonm=rhonmin
  eovera=eoverain
  compr=comprin
  effmas=effmasin
  asymme=asymmein
  dasymm=dasymmin
  skappa=skappain
  bmassp=bmasspin
  bmassn=bmassnin
  rho0fy=rho0fyin
  sigmavol=sigmavolin
  h2vm=h2vmin

! Deactivate kinetic term
  c0tau=0D0
  c1tau=0D0

! Isoscalar and isovector inverse masses
  bmass0 = 0.5D0*(bmassn+bmassp)
  bmass1 = 0.5D0*(bmassn-bmassp)

! A couple of Fayans scaling parameters
  rho0fy2 = 2D0*rho0fy
  r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
  r0fy2 = r0fy*r0fy
  efermfy = (1.125D0*pi)**(2D0/3D0)*bmass0/r0fy2
  efacfy = efermfy*rho0fy*(2D0/3D0)
  IF(ttest) WRITE(6,'(a,3(1pg13.5))') 'eferm,efac,r0fy=',efermfy,efacfy,r0fy
!  pfacfy = 2D0*efermfy/(3D0*rho0fy)
!  rho0inv = 0.5D0/rho0fy

!  tests
   IF(ttest) THEN
     DO i=1,50
      CALL testfunct(i*0.004D0,energ0,eslope0,ftest)
      WRITE(6,'(5(1pg16.9))') i*0.004D0,energ0,ftest
     END DO
     CALL testfunct(rhonm-1D-5,energ2,eslope2,ekin2)
     CALL testfunct(rhonm,energ0,eslope0,ekin0)
     CALL testfunct(rhonm+1D-5,energ1,eslope1,ekin1)
     WRITE(6,'(a,5(1pg14.7))') 'TEST:',energ0,eslope0,(energ1-energ2)/2D-5,&
           (energ1+energ2-2D0*energ0)/(1D-5)**2*rhonm**2*9D0,&
           (ekin1+ekin2-2D0*ekin0)/(1D-5)**2*rhonm**2*9D0
!  avp= -9.49577
!  h1vp= -5.95543/avp
!  h2vp= -1.37930/avp 
     h2var=-1.37930D0/(-9.49577D0)
     ftest = compfunct(h2var)
     WRITE(6,'(a,2(1pg13.5))') 'test compfunct:',h2var,ftest
   END IF
   

! Linear eqs. for isoscalar parameters
  xp = rhonm/(2D0*rho0fy)
  xps = xp**sigmavol
  IF(ttest) WRITE(6,'(a,3(1pg13.5))') 'rhonm,xp,xps=',rhonm,xp,xps

!  dp = 1D0+h2vp*xps
!  ap = ((sigmavol+dp)*xps*rhs-xps*xp*rhsp)/(sigmavol*dp*xp*xps)
!  hp = (sigmavol*dp-sigmavol-dp)*rhs-xp*rhsp)/(sigmavol*dp*xp*xps)
            
  IF(.NOT.ttest) THEN
   DO i=-1,10
     h2var=i*1D0
     ftest=compfunct(h2var)
     WRITE(6,'(2(1pg13.5))') h2var,ftest
   END DO
  END IF

  h2var=1.37930/9.49577
     ftest=compfunct(h2var)
     WRITE(6,'(2(1pg13.5))') h2var,ftest
  CALL bisect(h2var,0D0,1D1)

  IF(ttest) WRITE(6,'(a,4(1pg13.5))') 'ap,hp=',ap,hp,efermfy,efacfy

  avp=ap/efacfy
  h1vp=hp/ap
  h2vp=h2var

  c0t=c0tau
  c1t=c1tau  

  CALL timeoddparams(h2vm)
  avm=am/efacfy
  h1vm=hm/am

END SUBROUTINE NMP_Fayans

SUBROUTINE bisect(xout,xlow,xup)
! Finds roots of 'compfunct(x)=0'by bisection.
! 'xlow'and 'xup' define intervall for search.

REAL(8), INTENT(IN) :: xlow,xup
REAL(8), INTENT(OUT) :: xout

INTEGER, PARAMETER :: maxiter=50

INTEGER :: iter
REAL(8) :: f1,f3,f2,x1,x2,x3

  x1=xlow
  f1=compfunct(x1)
  x2=xup
  f2=compfunct(x2)
  IF(ttest) WRITE(6,'(a,6(1pg13.5))') 'bisection:',x1,x2,f1,f2
  IF(f1*f2>0D0) THEN
    STOP "bisection did no embrace zero"
  ELSE IF(f1==0D0) THEN
    xout=xlow
    RETURN
  ELSE IF(f2==0D0) THEN
    xout=xup
    RETURN
  END IF
  x3=0.5D0*(x1+x2)

  IF(ttest) WRITE(6,'(a)') 'bisection: iter,x1,x2,x3,f1,f2,f3'
  DO iter=1,maxiter
    f3=compfunct(x3)
!    IF(ttest) WRITE(6,'(i4,6(1pg13.5))') iter,x1,x2,x3,f1,f2,f3
    IF(f1*f3<0D0) THEN
      x2 = x3
      f2 = f3
    ELSE IF(f2*f3<0D0) THEN
      x1 = x3
      f1 = f3
    ELSE IF(f3==0D0) THEN
      xout=x3
      RETURN
    ELSE
      STOP "undecided in BISECTION"
    END IF
    x3 = 0.5D0*(x1+x2)
    f3 = compfunct(x3)
  END DO

  xout = x3
  IF(ttest) WRITE(6,'(a,2(1pg13.5))') 'BISECT:',ap,hp

  RETURN

END SUBROUTINE bisect

REAL(8) FUNCTION compfunct(h2var)
REAL(8), INTENT(IN OUT) :: h2var

REAL(8) :: det


! solution for 'ap' and 'hp' coefficients
  dp = 1D0+h2var*xps
  rhs = (rho0fy2*dp)*(eovera-betap*rhonm**r23*(bmass0+c0tau*rhonm))/xp
  rhsp = -(rho0fy2*dp)**2*&
         (r23*betap*rhonm**r23*(bmass0+0.6D0*c0tau*rhonm))/rhonm
  det = sigmavol*dp
  ap = ((sigmavol+dp)*rhs-rhsp)/det
  hp = ((sigmavol+dp-sigmavol*dp)*rhs-rhsp)/det/xps
!  IF(ttest) WRITE(6,'(a,4(1pg13.5))') 'xps,rhs,rhsp=',xps,rhs,rhsp
!  IF(ttest) WRITE(6,'(a,6(1pg13.5))') 'COMPFUNCT:',h2var,dp,ap/efacfy,hp/ap,rhs,rhsp

! final incompressibility matching
  compfunct = betap*rhonm**r23*(-2D0*bmass0+10D0*c0tau*rhonm) &
      +(-9D0*sigmavol*(1D0+sigmavol)*xps*xp*hp/dp &
      -18D0*(ap-(sigmavol+1D0)*hp*xps)*sigmavol*h2var*xps*xp/dp**2 &
      +9D0*(ap-hp*xps)*xp*(2D0*(sigmavol*h2var*xps)**2/dp**3 - &
              sigmavol*(sigmavol-1D0)*h2var*xps/dp**2    ))/rho0fy2 &
      -compr

  IF(ttest) WRITE(6,'(a,7(1pg13.5))') 'compr:',&
      betap*rhonm**r23*(-2D0*bmass0+10D0*c0tau*rhonm),&
      -9D0*sigmavol*(1D0+sigmavol)*xps*xp*hp/dp/rho0fy2, &
      -18D0*(ap-sigmavol*hp*xps)*sigmavol*h2var*xps*xp/dp**2/rho0fy2, &
      +9D0*(ap-hp*xps)*xp*(2D0*(sigmavol*h2var*xps)**2/dp**3 - &
              sigmavol*(sigmavol-1D0)*h2var*xps/dp**2    )/rho0fy2, &
      compr

END FUNCTION compfunct

SUBROUTINE timeoddparams(h2varm)
REAL(8), INTENT(IN) :: h2varm

REAL(8) :: det,c1contr,dm,xm


! solution for 'am' and 'hm' coefficients
  dm = 1D0+h2varm*xp
  c1contr=c1tau*(5D0/3D0)*betap*rhonm**(-0.3333333333333D0)
  rhs = rho0fy2*dm*(asymme-bmass0*5D0/9D0*betap*rhonm**0.66666666D0&
                    +c1contr)/xp
  rhsp = rho0fy2*dm**2*(dasymm-&
              bmass0*10D0/9D0*betap*rhonm**0.66666666D0+c1contr/2D0)/3D0/xp
  am = -(-(2D0+h2varm*xp)*rhs+rhsp)/dm
  hm = -(-rhs+rhsp)/dm/xp
  IF(ttest) WRITE(6,'(a,6(1pg13.5))') 'TIMEODD:',h2varm,dm,xp,rhs,rhsp
  IF(ttest) WRITE(6,'(a,6(1pg13.5))') 'TIMEODD:',asymme,dasymm
  IF(ttest) WRITE(6,'(a,6(1pg15.7))') 'TIMEODD:',&
    bmass0*5D0/9D0*betap*rhonm**0.66666666D0+xp*(am-hm*xp)/dm/rho0fy2,&
    bmass0*10D0/9D0*betap*rhonm**0.66666666D0+&
    3D0*xp*(am-hm*xp*(2D0+h2varm*xp))/dm**2/rho0fy2,&
    bmass0*10D0/9D0*betap*rhonm**0.66666666D0,&
    3D0*xp*(am-hm*xp*(2D0+h2varm*xp))/dm**2/rho0fy2,&
    3D0*xp*((am-2D0*hm*xp)/dm-(am-hm*xp)*h2varm*xp/dm**2)/rho0fy2

!   asy(0) = bmass0*5D0/9D0*betap*rhonm**0.66666666D0 +xp*(am-hm*xp)/dm/rho0fy2
  asy(0) = xp*(am-hm*xp)/dm/rho0fy2
  IF(ttest) WRITE(6,'(a,6(1pg15.7))') 'TIMEODD:',rhonm,asy(0)

  rhonm=rhonm+1D-4
  xp = rhonm/(2D0*rho0fy)
  xps = xp**sigmavol
  dm = 1D0+h2varm*xp
  asy(3) =  bmass0*5D0/9D0*betap*rhonm**0.66666666D0+xp*(am-hm*xp)/dm/rho0fy2
  rhonm=rhonm-2D-4
  xp = rhonm/(2D0*rho0fy)
  xps = xp**sigmavol
  dm = 1D0+h2varm*xp
  asy(1) =  bmass0*5D0/9D0*betap*rhonm**0.66666666D0+xp*(am-hm*xp)/dm/rho0fy2
  rhonm=rhonm+1D-4
  xp = rhonm/(2D0*rho0fy)
  xps = xp**sigmavol
  dm = 1D0+h2varm*xp
  WRITE(6,'(a,1pg13.5)') 'slope:',3D0*rhonm*(asy(3)-asy(1))/2D-4

!  am=efacfy*12.4287
!  hm=(-17.177052/12.4287)*am
  xm = 2D-3/rho0fy2
  asy(3) = (am-hm*xp)/(1D0+h2varm*xp)*xm**2/xp/rho0fy2
  asy(2) = 0D0
  asy(1) = (am-hm*xp)/(1D0+h2varm*xp)*xm**2/xp/rho0fy2
  WRITE(6,'(a,7(1pg13.5))') '2.deriv:',rhonm,&
       rhonm**2*(asy(3)-2D0*asy(2)+asy(1))/8D-6,&
       am/efacfy,hm/am,h2varm

!  am=efacfy*12.4287
!  hm=(-17.177052/12.4287)*am
  dm = 1D0+h2varm*xp
  IF(ttest) WRITE(6,'(a,6(1pg15.7))') 'direct:',&
    bmass0*5D0/9D0*betap*rhonm**0.66666666D0+xp*(am-hm*xp)/dm/rho0fy2,&
    bmass0*10D0/9D0*betap*rhonm**0.66666666D0+&
    3D0*xp*(am-hm*xp*(2D0+h2varm*xp))/dm**2/rho0fy2,&
    +xp*(am-hm*xp)/dm/rho0fy2



END SUBROUTINE timeoddparams


SUBROUTINE testfunct(rhoact,energ,eslope,ekin)

REAL(8), INTENT(IN) :: rhoact

REAL(8), INTENT(OUT) :: energ,eslope,ekin

REAL(8) :: avp=-9.49577,h1vp=-5.95543,h2var=-1.37930, &
           sigmavol=0.333333
REAL(8) :: avm=12.4287,h1vm= -17.1771/12.4287,dm
REAL(8) :: rho0fya2=0.16D0

  ap=avp*efacfy
  hp=0D0 !h1vp*ap/avp
  h2var=0D0 !-1.37930/avp
!  WRITE(6,'(a,5(1pg16.9))') 'avp,h1vp,h2vp=',avp,h1vp,h2var,ap

  xp = rhoact/rho0fya2
  xps = xp**sigmavol
  energ = bmass0*betap*rhoact**r23 + &
          (ap-hp*xps)/(1D0+h2var*xps) * xp/rho0fya2
  ekin = bmass0*betap*rhoact**r23 
  eslope = rho0fy2**2*r23*bmass0*betap*rhoact**r23/rhoact + &
              (ap*(1D0-(1D0-sigmavol)*h2var*xps)  -  &
               hp*xps*(1D0+sigmavol-sigmavol*h2var*xps))/ &
              (1D0+h2var*xps)**2
  dp=1D0+h2var*xps
!  eslope = (rho0fy2*dp)**2*r23*bmass0*betap*rhoact**r23/rhoact + &
!           (ap-(1D0+sigmavol)*hp*xps)*dp -&
!           (ap-hp*xps)*sigmavol*h2var*xps
  eslope = (rho0fy2*dp)**2*r23*bmass0*betap*rhoact**r23/rhoact + &
           (dp+sigmavol-sigmavol*dp)*ap -&
           (dp+sigmavol)*hp*xps
  WRITE(6,'(a,5(1pg16.7))') 'in TEST: ',rhoact,energ,ekin,&
          (ap-hp*xps)/(1D0+h2var*xps) * xp/rho0fya2
!  WRITE(6,'(a,5(1pg16.7))') 'EKIN:',rhoact,bmass0,betap,rhoact**r23 
!  WRITE(6,'(a,5(1pg15.7))') 'in TEST: ',rhoact,energ,&
!     bmass0*betap*rhoact**r23,&
!     (ap-hp*xps)/(1D0+h2var*xps) * xp/rho0fya2

  h2vm= 240.973/avm
  dm=1D0+h2vm*xp
  am=avm*efacfy
  hm=h1vm*am
  WRITE(6,'(a,7(1pg13.5))') 'TEST asymm:',&
    bmass0*5D0/9D0*betap*rhoact**0.66666666D0+ xp*(am-hm*xp)/dm/rho0fya2,&
   bmass0*10D0/9D0*betap*rhoact**0.66666666D0&
   +3D0*(am-hm*xp*(2D0+h2vm*xp))/(dm*dm)/rho0fya2,&
   bmass0*10D0/9D0*betap*rhoact**0.66666666D0,&
   xp*(am-hm*xp)/dm/rho0fya2,&
   +3D0*(am-hm*xp*(2D0+h2vm*xp))/(dm*dm)/rho0fya2


   

END SUBROUTINE testfunct

!  avp= -9.49577
!  h1vp= -5.95543/avp
!  h2vp= -1.37930/avp 
!  sigmavol= 0.333333
!  rho0fy= 0.800000E-01


END MODULE NMP_map
