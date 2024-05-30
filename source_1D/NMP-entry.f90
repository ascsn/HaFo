
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
REAL(8),PRIVATE :: xp,xps,h2var,dp,efacfy,rho0fy2     ! auxiliary Fayans parameters
REAL(8),PRIVATE :: c0delrho,c1delrho                  ! axuiliary Skyrme parameter 
REAL(8),PRIVATE :: c0rho,c1rho,c0alpha,c1alpha        ! auxiliary Skyrme parameters
REAL(8),PRIVATE :: rho23,rhoalpha,compralp   ! auxiliary Skyrme parameters
REAL(8),PRIVATE :: rhonm,eovera,compr,effmas,asymme,dasymm,d2asymm,skappa,&
           bmassp,bmassn,rho0fy,sigmavol
REAL(8),PRIVATE :: rhs,rhsp,asy(0:3)
LOGICAL,PRIVATE,PARAMETER :: ttest=.FALSE.
LOGICAL,PRIVATE,PARAMETER :: ttest2=.FALSE.

CONTAINS
SUBROUTINE NMP_Skyrme(t0,t1,t2,t3,x0,x1,x2,x3,alpha,&
                      rhonmin,eoverain,comprin,effmasin,asymmein,dasymmin,skappain,&
                      bmasspin,bmassnin,c0delrhoin,c1delrhoin)

! NMP as input, also input the fixed Fayans parameter 'rho0fy'
REAL(8), INTENT(IN) :: rhonmin,eoverain,comprin,effmasin,asymmein,dasymmin,&
                       skappain,bmasspin,bmassnin,c0delrhoin,c1delrhoin
! Resulting Fayans volume parameters
REAL(8), INTENT(OUT) :: t0,t1,t2,t3,x0,x1,x2,x3,alpha

REAL(8) :: t0ot,t1ot,t2ot,t3ot,x0ot,x1ot,x2ot,x3ot,eta
!REAL(8) :: compfunct
!REAL(8) :: rhoalpha,compralp
REAL(8):: comptest
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
  c0delrho=c0delrhoin
  c1delrho=c1delrhoin


!
  IF(ttest) WRITE(6,'(a,12(1pg14.6))') &
     'NMP_Skrme called with: rhonm,eovera,compr,asymme'//&
     'dasymm,bmassp,bmassn,effmas,skappa=',&
      rhonm,eovera,compr,asymme,dasymm,bmassp,bmassn,effmas,skappa

! Isoscalar and isovector inverse masses
  bmass0 = 0.5D0*(bmassn+bmassp)
  bmass1 = 0.5D0*(bmassn-bmassp)

! kinetic term
  c0tau=bmass0/rhonm*(1D0/effmas-1D0)
  c1tau=bmass0/rhonm*(1D0/effmas-1D0-skappa)
!  c1tau=0D0
  IF(ttest) WRITE(6,'(a,2(1pg13.5))') 'c0tau,c1tau=',c0tau,c1tau

! isoscalar density terms
  rho23=rhonm**0.666666667D0
  rhs=EoverA/rhonm-bmass0*betap*rho23/rhonm-c0tau*betap*rho23
  rhsp=-(2D0*bmass0*betap*rho23/rhonm+5D0*c0tau*betap*rho23)/3D0
  IF(ttest) WRITE(6,'(a,3(1pg13.5))') 'rhs,rhsp,rho23=',rhs,rhsp,rho23

  IF(ttest2) THEN
    DO i=1,100
      comptest=compfunctSk(i*0.01D0)
      WRITE(6,'(a,1pg13.5)') 'COMPFUNCT:',comptest
    END DO
  END IF

  CALL bisect(compfunctSk,alpha,1D-2,8D0)
! ==> result: c0rho,c0alpha,alpha

! isovector density terms
!  DO i=-70,-60
!    eta=i/100D0  
   eta=-1D0   !   problem with sign ???
  rhoalpha=rhonm**alpha
  rhs=asymme/rhonm-bmass0*betap*rho23*5D0/9D0/rhonm&
       -c0tau*betap*rho23*5D0/9D0 &
      +eta*c1tau*betap*rho23*5D0/3D0
  rhsp=dasymm/(3D0*rhonm)-bmass0*betap*rho23*10D0/27D0/rhonm&
       -c0tau*betap*rho23*25D0/27D0 &
       +eta*c1tau*betap*rho23*25D0/9D0
  c1alpha=(rhsp-rhs)/(alpha*rhoalpha)
  c1rho=(rhs-c1alpha*rhoalpha)

  IF(ttest) WRITE(6,'(a,f5.2,4(1pg13.5))') &
    ' eta,c0rho,c0alpha,c1rho,c1alpha=',eta,c0rho,c0alpha,c1rho,c1alpha
!  END DO

  IF(ttest) WRITE(6,'(a,4(1pg13.5))') &
    ' c0tau,c0delrho,c1tau,c1delrho=',c0tau,c0delrho,c1tau,c1delrho

  CALL mapSk_CtoT(c0rho,c0alpha,c1rho,c1alpha,c0tau,c1tau,c0delrho,c1delrho,&
                  t0ot,t1ot,t2ot,t3ot,x0ot,x1ot,x2ot,x3ot)

  t0=t0ot
  t1=t1ot
  t2=t2ot
  t3=t3ot
  x0=x0ot
  x1=x1ot
  x2=x2ot
  x3=x3ot

  RETURN

END SUBROUTINE NMP_Skyrme

REAL(8) FUNCTION compfunctSk(alpha)
REAL(8), INTENT(IN) :: alpha

REAL(8) :: det

  rhoalpha=rhonm**alpha
  c0rho=((1D0+alpha)*rhs-rhsp)/(alpha)
  c0alpha=(-rhs+rhsp)/(alpha*rhoalpha)
  compralp=-2D0*bmass0*betap*rho23+10D0*c0tau*betap*rhonm*rho23&
           +9D0*alpha*(1D0+alpha)*c0alpha*rhonm*rhoalpha
  IF(ttest) WRITE(6,'(a,7(1pg13.5))') &
     'COMP: alpha,rhoalpha,rhs,rhsp,c0rho,c0alpha,compralp=',&
     alpha,rhoalpha,rhs,rhsp,c0rho,c0alpha,compralp
  compfunctSk=compralp-compr


END FUNCTION compfunctSk


SUBROUTINE mapSk_CtoT(c0rhoin,c0alphain,c1rhoin,c1alphain,c0tauin,c1tauin,&
                      c0delrhoin,c1delrhoin,&
                      t0,t1,t2,t3,x0,x1,x2,x3)
REAL(8), INTENT(IN) :: c0rhoin,c0alphain,c1rhoin,c1alphain,c0tauin,c1tauin,&
                       c0delrhoin,c1delrhoin
REAL(8), INTENT(OUT) :: t0,t1,t2,t3,x0,x1,x2,x3

  t0=c0rhoin*8D0/3D0
  x0=-0.5D0*(c0rhoin+3D0*c1rhoin)/c0rhoin
  t3=16D0*c0alphain
  x3=-0.5D0*(c0alphain+3D0*c1alphain)/c0alphain
  t1=(c0tauin-4D0*c0delrhoin)*4D0/3D0
  x1=(-3D0*c1tauin+12D0*c1delrhoin-c0tauin+4D0*c0delrhoin)*2D0/(3D0*t1)
  t2=(3D0*c0tauin+4D0*c0delrhoin-6D0*c1tauin-8D0*c1delrhoin)*4D0/3D0
  x2=(6D0*c1tauin+8D0*c1delrhoin-0.5D0*t2)/t2

END SUBROUTINE mapSk_CtoT

SUBROUTINE mapSk_TtoC(c0rhoot,c0alphaot,c1rhoot,c1alphaot,c0tauot,c1tauot,&
                      c0delrhoot,c1delrhoot,&
                      t0,t1,t2,t3,x0,x1,x2,x3)
REAL(8), INTENT(OUT) :: c0rhoot,c0alphaot,c1rhoot,c1alphaot,c0tauot,c1tauot,&
                       c0delrhoot,c1delrhoot
REAL(8), INTENT(IN) :: t0,t1,t2,t3,x0,x1,x2,x3

  c0rhoot=t0*3D0/8D0
  c0alphaot=t3/16D0
  c0tauot=t1*3D0/16D0+t2*5D0/16D0+t2*x2/4D0
  c0delrhoot=-t1*9D0/64D0+t2*5D0/64D0+t2*x2/16D0

  c1rhoot=-t0*(0.5D0+x0)/4D0
  c1alphaot=-t3*(0.5D0+x3)/24D0
  c1tauot=-(t1*(0.5D0+x1)-t2*(0.5D0+x2))/8D0
  c1delrhoot=(3D0*t1*(0.5D0+x1)+t2*(0.5D0+x2))/32D0

END SUBROUTINE mapSk_TtoC


SUBROUTINE NMP_Fayans(avp,h1vp,h2vp,avm,h1vm,h2vmout,h2vmin,sigmavolin,&
                      c0t,c1t,rho0fyin,&
                      rhonmin,eoverain,comprin,effmasin,asymmein,dasymmin,d2asymmin,skappain,&
                      bmasspin,bmassnin,tksyminput)

! NMP as input, also input the fixed Fayans parameter 'rho0fy'
! The option 'tksyminpout' switches to new version where K_sym='d2asymmin' is the input
! parameter and 'h2vm' output; the old option uses 'h2vmin' as input.
!
REAL(8), INTENT(IN) :: rhonmin,eoverain,comprin,effmasin,asymmein,dasymmin,d2asymmin,&
                       skappain,bmasspin,bmassnin,rho0fyin,sigmavolin,h2vmin
! Resulting Fayans volume parameters
REAL(8), INTENT(OUT) :: avp,h1vp,h2vp,avm,h1vm,c0t,c1t,h2vmout
LOGICAL, INTENT(IN) :: tksyminput    ! switches to input 'd2asymmin' instead of 'h2vmin'

!REAL(8) :: compfunct
REAL(8) :: r0fy,r0fy2,efermfy,ftest,Ksym,Jkin,Lkin,Ksymkin
REAL(8):: energ0,eslope0,energ1, eslope1, energ2,eslope2,ekin0,ekin1,ekin2,dm
REAL(8),PARAMETER :: ksymfac=6D0
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
  IF(tksyminput) d2asymm=d2asymmin
!
  IF(ttest) WRITE(6,'(a,12(1pg14.6))') &
     'NMP_Fayans called with: rhonm,eovera,compr,asymme'//&
     'dasymm,bmassp,bmassn,rho0fy,sigmavol,h2vm,effmas,skappa=',&
      rhonm,eovera,compr,asymme,dasymm,bmassp,bmassn,rho0fy,sigmavol,&
      h2vm,effmas,skappa


! Isoscalar and isovector inverse masses
  bmass0 = 0.5D0*(bmassn+bmassp)
  bmass1 = 0.5D0*(bmassn-bmassp)

! kinetic term
  c0tau=bmass0/rhonm*(1D0/effmas-1D0)
  c1tau=bmass0/rhonm*(1D0/effmas-1D0-skappa)
  c0t=c0tau
  c1t=c1tau
!  c1tau=0D0
  IF(ttest) WRITE(6,'(a,2(1pg13.5))') 'c0tau,c1tau=',c0tau,c1tau
! Deactivate kinetic term
!  IF(effmas.NE.1D0 .OR. skappa.NE.0D0) STOP "NMP entry not yet for m*/m"
!  c0tau=0D0
!  c1tau=0D0

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
   IF(ttest2) THEN
     DO i=1,50
      CALL testfunctFy(i*0.004D0,energ0,eslope0,ftest)
      WRITE(6,'(5(1pg16.9))') i*0.004D0,energ0,ftest
     END DO
     CALL testfunctFy(rhonm-1D-5,energ2,eslope2,ekin2)
     CALL testfunctFy(rhonm,energ0,eslope0,ekin0)
     CALL testfunctFy(rhonm+1D-5,energ1,eslope1,ekin1)
     WRITE(6,'(a,5(1pg14.7))') 'TEST:',energ0,eslope0,(energ1-energ2)/2D-5,&
           (energ1+energ2-2D0*energ0)/(1D-5)**2*rhonm**2*9D0,&
           (ekin1+ekin2-2D0*ekin0)/(1D-5)**2*rhonm**2*9D0
!  avp= -9.49577
!  h1vp= -5.95543/avp
!  h2vp= -1.37930/avp 
!old     h2var=-1.37930D0/(-9.49577D0)
!avp= -10.436376
!h1vp= 0.54178836    
!h2vp= 0.42032398
     h2var=0.42032398  !/(-10.436376)
     ftest = compfunctFy(h2var)
     WRITE(6,'(a,2(1pg13.5))') 'test compfunct:',h2var,ftest
   END IF
   

! Linear eqs. for isoscalar parameters
  xp = rhonm/(2D0*rho0fy)
  xps = xp**sigmavol
  IF(ttest) WRITE(6,'(a,3(1pg13.5))') 'rhonm,xp,xps=',rhonm,xp,xps

!  dp = 1D0+h2vp*xps
!  ap = ((sigmavol+dp)*xps*rhs-xps*xp*rhsp)/(sigmavol*dp*xp*xps)
!  hp = (sigmavol*dp-sigmavol-dp)*rhs-xp*rhsp)/(sigmavol*dp*xp*xps)
            
!  IF(.NOT.ttest) THEN
!   DO i=-1,10
!     h2var=i*1D0
!     ftest=compfunct(h2var)
!     WRITE(6,'(2(1pg13.5))') h2var,ftest
!   END DO
!  END IF

  h2var=1.37930/9.49577      ! starting value
  IF(ttest) THEN
     ftest=compfunctFy(h2var)
     WRITE(6,'(2(1pg13.5))') h2var,ftest
  END IF
  CALL bisect(compfunctFy,h2var,0D0,1D1)

  IF(ttest) WRITE(6,'(a,4(1pg13.5))') 'ap,hp=',ap,hp,efermfy,efacfy

  avp=ap/efacfy
  h1vp=hp/ap
  h2vp=h2var

  c0t=c0tau
  c1t=c1tau  

  IF(tksyminput) THEN
    IF(ttest) THEN 
      WRITE(6,'(a,6(1pg13.5))') 'enter BISECT-1:', &
       rhonm,eovera,compr,asymme,dasymm,d2asymm
!  
      DO i=0,-1 !16
        h2vm=EXP(i*0.5D0)
        CALL timeoddparams(h2vm)
        dm = 1D0+h2vm*xp
        WRITE(6,'(a,3(1pg13.5))') 'test K_sym:',h2vm,    &
         10D0/9D0*betap*rhonm**0.66666666D0 * &
         (-bmass0+(c0tau*5D0+c1tau*15D0)*rhonm) &
        -(9D0*xp**2)/rho0fy  * (hm+h2vm*am)/dm**3,&
        d2asymm
      END DO
    END IF
!
    Jkin= 5D0/9D0*betap*rhonm**0.66666666D0 * &
         (bmass0+(c0tau+c1tau*3D0)*rhonm) 
    Lkin =   10D0/9D0*betap*rhonm**0.66666666D0 * &
       (bmass0+(c0tau*2.5D0+c1tau*7.5D0)*rhonm) 
    Ksymkin =   10D0/9D0*betap*rhonm**0.66666666D0 * &
       (-bmass0+(c0tau*5D0+c1tau*15D0)*rhonm) 
    h2vm=(ksymfac*(dasymm-Lkin-3D0*(asymme-Jkin))/(d2asymm-Ksymkin)-1D0)/xp
        CALL timeoddparams(h2vm)
!    CALL bisect(comp1functFy,h2var,0D0,1D2)
!    h2var=EXP(h2var)
    Ksym=Ksymkin+ksymfac*(dasymm-Lkin-3D0*(asymme-Jkin))/(1D0+h2vm*xp) 
    h2vmout=h2vm
    dm = 1D0+h2vm*xp
    IF(ttest) WRITE(2,'(a,5(1pg13.5))') 'test Ksym:',&
       h2vm,h2var,Ksym,Ksymkin+ksymfac*(dasymm-Lkin-3D0*(asymme-Jkin))/(1D0+h2var*xp),&
         10D0/9D0*betap*rhonm**0.66666666D0 * &
         (-bmass0+(c0tau*5D0+c1tau*15D0)*rhonm) &
        -(9D0*xp**2)/rho0fy  * (hm+h2vm*am)/dm**3
!    h2vmout=h2var
  ELSE
    CALL timeoddparams(h2vm)
  END IF
  avm=am/efacfy
  h1vm=hm/am
  IF(ttest) THEN
    d2asymm=Ksymkin+ksymfac*(dasymm-Lkin-3D0*(asymme-Jkin))/(1D0+h2vm*xp)
    WRITE(2,'(a,5(1pg13.5))') 'test parm:',am,hm,avm,h1vm,h2vm
    WRITE(2,'(a,6(1pg13.5))') 'NMP:',asymme,dasymm,d2asymm,Jkin,Lkin,Ksymkin
    WRITE(2,'(a,5(1pg13.5))') ' Ksym_new=',&
    Ksymkin+ksymfac*(dasymm-Lkin-3D0*(asymme-Jkin))/(1D0+h2vm*xp),&
    (ksymfac*(dasymm-Lkin-3D0*(asymme-Jkin))/(d2asymm-Ksymkin)-1D0)/xp,d2asymm
  END IF
  IF(ttest) WRITE(6,'(a,4(1pg13.5))') 'avp,h1vp,h2vp=',avp,h1vp,h2vp
  IF(ttest) WRITE(6,'(a,4(1pg13.5))') 'avm,h1vm,h2vm=',avm,h1vm,h2vmout

!  test Ksymm
  IF(1==2 .AND. .NOT.tksyminput) THEN
    dm = 1D0+h2vm*xp
    Jkin= 5D0/9D0*betap*rhonm**0.66666666D0 * &
       (bmass0+(c0tau+c1tau*3D0)*rhonm) 
    d2asymm = 10D0/9D0*betap*rhonm**0.66666666D0 * &
       (-bmass0+(c0tau*5D0+c1tau*15D0)*rhonm) &
      -(9D0*xp**2)/rho0fy  * (hm+h2vm*am)/dm**3 
    Lkin =   10D0/9D0*betap*rhonm**0.66666666D0 * &
       (bmass0+(c0tau*2.5D0+c1tau*7.5D0)*rhonm) 
    Ksymkin =   10D0/9D0*betap*rhonm**0.66666666D0 * &
       (-bmass0+(c0tau*5D0+c1tau*15D0)*rhonm) 
    Ksym=Ksymkin+ksymfac*(Jkin+dasymm-3D0*asymme)/(1D0+h2vm*xp) 
    CALL bisect(comp1functFy,h2var,0D0,1D2)
    WRITE(6,'(a,5(1pg13.5))') 'K_symm,H2VM=',d2asymm,EXP(h2var),h2var,h2vm,-10D0/9D0*betap*rhonm**0.66666666D0*bmass0
    WRITE(2,'(a,5(1pg13.5))') 'K_symm,H2VM=',d2asymm,h2vm,&
      -(9D0*xp**2)/rho0fy  * (hm+h2vm*am)/dm**3,Ksym,5D0/9D0*betap*rhonm**0.66666666D0*bmass0
    WRITE(2,'(a,5(1pg13.5))') 'K_symm-shrt=',&
      Ksymkin+ksymfac*(dasymm-Lkin-3D0*(asymme-Jkin))/(1+h2vm*xp),am,hm,avm,h1vm
  END IF

END SUBROUTINE NMP_Fayans

SUBROUTINE bisect(compfunct,xout,xlow,xup)
! Finds roots of 'compfunct(x)=0'by bisection.
! 'xlow'and 'xup' define intervall for search.

REAL(8) :: compfunct
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
    IF(ttest) WRITE(6,'(i4,6(1pg13.5))') iter,x1,x2,x3,f1,f2,f3
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
!      WRITE(6,'(i4,6(1pg13.5))') iter,x1,x2,x3,f1,f2,f3
      STOP "undecided in BISECTION"
    END IF
    x3 = 0.5D0*(x1+x2)
    f3 = compfunct(x3)
  END DO

  xout = x3
  IF(ttest) WRITE(6,'(a,2(1pg13.5))') 'BISECT:',ap,hp

  RETURN

END SUBROUTINE bisect

REAL(8) FUNCTION compfunctFy(h2var)
REAL(8), INTENT(IN OUT) :: h2var

REAL(8) :: det


! solution for 'ap' and 'hp' coefficients
  dp = 1D0+h2var*xps
  rhs = (rho0fy2*dp)*(eovera-betap*rhonm**r23*(bmass0+c0tau*rhonm))/xp
  rhsp = -(rho0fy2*dp)**2*&
         (r23*betap*rhonm**r23*(bmass0+2.5D0*c0tau*rhonm))/rhonm
  det = sigmavol*dp
  ap = ((sigmavol+dp)*rhs-rhsp)/det
  hp = ((sigmavol+dp-sigmavol*dp)*rhs-rhsp)/det/xps
!  IF(ttest) WRITE(6,'(a,4(1pg13.5))') 'xps,rhs,rhsp=',xps,rhs,rhsp
!  IF(ttest) WRITE(6,'(a,6(1pg13.5))') 'COMPFUNCT:',h2var,dp,ap/efacfy,hp/ap,rhs,rhsp

! final incompressibility matching
  compfunctFy = betap*rhonm**r23*(-2D0*bmass0+10D0*c0tau*rhonm) &
      +(-9D0*sigmavol*(1D0+sigmavol)*xps*xp*hp/dp &
      -18D0*(ap-(sigmavol+1D0)*hp*xps)*sigmavol*h2var*xps*xp/dp**2 &
      +9D0*(ap-hp*xps)*xp*(2D0*(sigmavol*h2var*xps)**2/dp**3 - &
              sigmavol*(sigmavol-1D0)*h2var*xps/dp**2    ))/rho0fy2 &
      -compr

  IF(ttest2) WRITE(6,'(a,7(1pg13.5))') 'compr:',&
      betap*rhonm**r23*(-2D0*bmass0+10D0*c0tau*rhonm),&
      -9D0*sigmavol*(1D0+sigmavol)*xps*xp*hp/dp/rho0fy2, &
      -18D0*(ap-sigmavol*hp*xps)*sigmavol*h2var*xps*xp/dp**2/rho0fy2, &
      +9D0*(ap-hp*xps)*xp*(2D0*(sigmavol*h2var*xps)**2/dp**3 - &
              sigmavol*(sigmavol-1D0)*h2var*xps/dp**2    )/rho0fy2, &
      compr

END FUNCTION compfunctFy


SUBROUTINE timeoddparams(h2varm)
REAL(8), INTENT(IN) :: h2varm

REAL(8) :: det,c0contr,c1contr,dm,xm,rho23,rho53,eta


! solution for 'am' and 'hm' coefficients
  rho23=rhonm**(2D0/3D0)
  rho53=rho23*rhonm
  dm = 1D0+h2varm*xp
  eta=-1D0
  c0contr=c0tau*(5D0/9D0)*betap*rho53
  c1contr=c1tau*(5D0/3D0)*betap*rho53
  rhs = rho0fy2*dm*(asymme-bmass0*5D0/9D0*betap*rho23&
                    -c0contr+eta*c1contr)/xp
  rhsp = rho0fy2*dm**2*(dasymm-&
           bmass0*10D0/9D0*betap*rho23-5D0*c0contr+eta*5D0*c1contr)/3D0/xp
  am = -(-(2D0+h2varm*xp)*rhs+rhsp)/dm
  hm = -(-rhs+rhsp)/dm/xp
  IF(ttest) WRITE(6,'(a,6(1pg13.5))') 'TIMEODD: h2varm...',h2varm,dm,xp,rhs,rhsp
  IF(ttest) WRITE(6,'(a,6(1pg13.5))') 'TIMEODD: asym...',asymme,dasymm
  IF(ttest) WRITE(6,'(a,6(1pg15.7))') 'TIMEODD:',&
    bmass0*5D0/9D0*betap*rhonm**0.66666666D0+xp*(am-hm*xp)/dm/rho0fy2,&
    bmass0*10D0/9D0*betap*rhonm**0.66666666D0+&
    3D0*xp*(am-hm*xp*(2D0+h2varm*xp))/dm**2/rho0fy2,&
    bmass0*10D0/9D0*betap*rhonm**0.66666666D0,&
    3D0*xp*(am-hm*xp*(2D0+h2varm*xp))/dm**2/rho0fy2,&
    3D0*xp*((am-2D0*hm*xp)/dm-(am-hm*xp)*h2varm*xp/dm**2)/rho0fy2

!   asy(0) = bmass0*5D0/9D0*betap*rhonm**0.66666666D0 +xp*(am-hm*xp)/dm/rho0fy2
!  asy(0) = xp*(am-hm*xp)/dm/rho0fy2
  IF(ttest) THEN

    WRITE(6,'(a,6(1pg15.7))') 'TIMEODD:',rhonm,asy(0)

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

!    am=efacfy*12.4287
!    hm=(-17.177052/12.4287)*am
    xm = 2D-3/rho0fy2
    asy(3) = (am-hm*xp)/(1D0+h2varm*xp)*xm**2/xp/rho0fy2
    asy(2) = 0D0
    asy(1) = (am-hm*xp)/(1D0+h2varm*xp)*xm**2/xp/rho0fy2
    WRITE(6,'(a,7(1pg13.5))') '2.deriv:',rhonm,&
       rhonm**2*(asy(3)-2D0*asy(2)+asy(1))/8D-6,&
       am/efacfy,hm/am,h2varm

!    am=efacfy*12.4287
!    hm=(-17.177052/12.4287)*am
    dm = 1D0+h2varm*xp
    WRITE(6,'(a,6(1pg15.7))') 'direct:',&
      bmass0*5D0/9D0*betap*rhonm**0.66666666D0+xp*(am-hm*xp)/dm/rho0fy2,&
      bmass0*10D0/9D0*betap*rhonm**0.66666666D0+&
      3D0*xp*(am-hm*xp*(2D0+h2varm*xp))/dm**2/rho0fy2,&
      +xp*(am-hm*xp)/dm/rho0fy2

  END IF

END SUBROUTINE timeoddparams


REAL(8) FUNCTION comp1functFy(h2varml)
REAL(8), INTENT(IN OUT) :: h2varml     ! log of h2vm

REAL(8) :: det,c0contr,c1contr,dm,xm,rho23,rho53,eta
REAL(8) :: h2varm

! isovector conmpressibility K_sym for given 'h2varm' and
! effective isovector parameters. The two other isovector
! parameters  'am' and 'hm' are also computed here by a 
! call to 'timeoddparams'.

  h2varm=EXP(h2varml)
!  h2varm=h2varml/(1D0-h2varml)
  CALL timeoddparams(h2varm)
  dm = 1D0+h2varm*xp

!  comp1functFy =
!      -bmass0*10D0/9D0*betap*rhonm**0.66666666D0+&
!      +c0tau*50D0/9D0*betap*rhonm**1.66666666D0+&
!      +c1tau*50D0/3D0*betap*rhonm**1.66666666D0+&
!      +(9D0*xp)**2/rho0fy2  * &
!       (hm+h2varm*am-2D0*hm*h2varm*xp)/dm*+3
  comp1functFy = 10D0/9D0*betap*rhonm**0.66666666D0 * &
       (-bmass0+(c0tau*5D0+c1tau*15D0)*rhonm) &
      -(9D0*xp**2)/rho0fy  * (hm+h2varm*am)/dm**3 &
      - d2asymm


END FUNCTION comp1functFy


SUBROUTINE testfunctFy(rhoact,energ,eslope,ekin)

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


   

END SUBROUTINE testfunctFy



END MODULE NMP_map
