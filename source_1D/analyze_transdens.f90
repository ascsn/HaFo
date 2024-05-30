PROGRAM analyze_transdens

IMPLICIT REAL(8) (a-h,o-z)

INTEGER, PARAMETER :: lrpa=1, kmode=200, kdat=6
REAL(8),PARAMETER :: pi=3.1415926D0
REAL(8),PARAMETER :: emin=8D0,emax=33D0,dele=0.1D0
REAL(8),PARAMETER :: ewidth=dele*4D0

REAL(8),ALLOCATABLE :: dens(:,:,:),formfa(:,:,:),ffl(:),fft(:)
CHARACTER(30) :: atext
CHARACTER(90) :: atextl
REAL(8) :: energy(3,kmode),accff(3),slope(kmode),average(kmode),bel(3,4,kmode)
REAL(8) :: ffdat(kdat,1:2,kmode),ffdat0(kdat,1:2),strength(4)
REAL(8) :: ein(0:kdat),angle(0:kdat), sigmalog(0:kdat)   ! data points
ein(0)=50.4; angle(0)= 15.0
ein(1)=49.4; angle(1)= 92.9
ein(2)=50.4; angle(2)= 92.9
ein(3)=49.9; angle(3)=116.9
ein(4)=49.4; angle(4)=140.9
ein(5)=49.9; angle(5)=164.9
ein(6)=50.4; angle(6)=164.9

! 
! read data and store
!

DO icase=1,3
  IF(icase==1) OPEN(1,file='transdens.res')
  IF(icase==2) OPEN(1,file='transcurm.res')
  IF(icase==3) OPEN(1,file='transcurp.res')

  READ(1,'(a)') atext
  WRITE(6,'(a)') atext
  READ(1,'(8x,i5)') ngrid
  WRITE(6,'(a,i5)') '# running with ngrid=',ngrid
  READ(1,'(1x)')
  
  IF(icase==1) ALLOCATE (dens(ngrid,1:3,kmode),formfa(ngrid,1:3,kmode),ffl(ngrid),fft(ngrid))
  
  DO iset=1,kmode

!    READ(1,'(a)')  atextl
!    WRITE(6,'(a)') atextl
!    BACKSPACE(1)
!    READ(1,'(10x,i5,9x,e12.5,7x,4e12.5)',ERR=99,END=99) mode,energy(icase,iset),bel(icase,:,iset)
    READ(1,'(10x,i5,9x,e12.5,7x,4e13.5)',END=99) mode,energy(icase,iset),bel(icase,:,iset)
    IF(icase>1 .AND. energy(1,iset).NE.energy(icase,iset)) STOP 'energies do not match'
    IF(icase==3) WRITE(6,'(a,i5,6(1pg13.5))') '# mode,energy:',mode,energy(icase,iset),bel(icase,:,iset)
!    WRITE(6,'(a,i5,1pg13.5)') '# mode,energy:',mode,energy(icase,iset)
  
    DO i=1,ngrid
      READ(1,*) x,dens(i,icase,iset)
    END DO
    READ(1,'(1x)')
  
    IF(iset==1) THEN
      rstep=x/(ngrid-1)
      IF(icase>1 .AND. rstep.NE.rstepold) STOP 'grid spacings do not match'
      frpirs = 4D0*pi*rstep
      rleng  = (ngrid-1)*rstep
      delq   = pi/rleng
      delql  = 0.5D0*delq
      rstepold=rstep
    END IF

    IF(icase==3) THEN  
      DO iq=1,ngrid
        q1     = iq*delql-delql
        IF(q1 == 0D0) q1 = delql/10D0
        r      = 0D0
        accff   = 0D0
        DO i=2,ngrid
          r      = rstep+r
          vol    = r*r  !
          accff(1)  = vol*bessel(r*q1,lrpa)*dens(i,1,iset)+accff(1)
          IF(lrpa>0) accff(2)  = vol*bessel(r*q1,lrpa-1)*dens(i,2,iset)+accff(2)
          accff(3)  = vol*bessel(r*q1,lrpa+1)*dens(i,3,iset)+accff(3)
        END DO
        accff   = frpirs*accff
        formfa(iq,:,iset) = accff/q1
        formfa(iq,1,iset) = SQRT(2D0*lrpa+1D0)*formfa(iq,1,iset)
!        formfa(iq,2,iset) = (SQRT(lrpa+1D0)*formfa(iq,2,iset)-SQRT(1D0*lrpa)*formfa(iq,3,iset))&
        formfa(iq,2,iset) = (SQRT(lrpa+1D0)*formfa(iq,2,iset)+SQRT(1D0*lrpa)*formfa(iq,3,iset))&
                            *197.32D0/938D0
        WRITE(6,'(1x,f8.3,f8.4,5(1pg13.5))') energy(1,iset),q1,formfa(iq,1:2,iset)
      END DO
      WRITE(6,'(1x)')
  !    WRITE(6,'(1x)')
    END IF
    maxset=iset

  END DO
99 CONTINUE

  CLOSE(1)
END DO
WRITE(6,'(a,i5)') '# max. modes=',maxset
!
! compute cross sections and slopes
!
OPEN(8,file='slopes')
WRITE(8,'(1x/1x/a)') '#mode energy slope'
DO iset=1,maxset
  IF(energy(1,iset)>ein(1)-10D0) EXIT
!  ffl(:)=formfa(:,1,iset)
!  fft(:)=formfa(:,2,iset)
  CALL slopes(energy(1,iset),slope(iset),average(iset),.FALSE.)
  WRITE(8,'(i5,f8.3,6(1pg13.5))') iset,energy(1,iset),slope(iset),average(iset),sigmalog(6)-sigmalog(0)
END DO
!
! Average to smooth distibution over energy and print
!
WRITE(8,'(1x/1x/a)')'# slopes from smoothed cross section'
nbin=(emax-emin)/dele
DO ibin=0,nbin
  ebin=emin+ibin*dele
  IF(ebin>ein(1)-10D0) EXIT
  ffdat0=0D0
  DO idat=1,kdat
    qeff = SQRT((1D0-COS(angle(idat)*pi/180D0))*2D0*ein(idat)*(ein(idat)-ebin)/192.32**2)
    iq=qeff/delql+1
    q1=iq*delql-delql
    q2=q1*delql
    strength=0D0
    DO iset=1,maxset
      weight=bel(3,1,iset)*EXP(-((energy(1,iset)-ebin)/ewidth)**2)
      ffdat(idat,:,iset)=((qeff-q1)*formfa(iq+1,:,iset)+(q2-qeff)*formfa(iq,:,iset))/delql
      ffdat0(idat,:)=ffdat0(idat,:)+weight*ffdat(idat,:,iset)**2
      strength=strength+bel(3,:,iset)*EXP(-((energy(1,iset)-ebin)/ewidth)**2)
!  WRITE(6,'(a,5(1pg13.5))') '#FFDAT0:',ebin,ffdat0(idat,:),ffdat(idat,:,iset)
    END DO
  END DO
  CALL  slopes(ebin,slope(1),average(1),.TRUE.)
  WRITE(8,'(9(1pg13.5))') ebin,slope(1),average(1),strength
END DO

STOP
CONTAINS

SUBROUTINE slopes(energact,slopeout,averout,tbin)
REAL(8),INTENT(IN) :: energact
REAL(8),INTENT(OUT) :: slopeout,averout 
LOGICAL,INTENT(IN) :: tbin
REAL(8),PARAMETER :: faclong=1D0, factrans=1D0

IF(.NOT.tbin) THEN
  WRITE(7,'(a,i4,1pg13.5)') '# mode:',iset,energy(1,iset)
  DO idat=0,kdat
    qeff = SQRT((1D0-COS(angle(idat)*pi/180D0))*2D0*ein(idat)*(ein(idat)-energact)/192.32**2)
    WRITE(9,'(i5,1pg13.5)') idat,qeff
    iq=qeff/delql+1
    IF(iq>ngrid-1) THEN
      WRITE(6,'(a,2i5,3(1pg13.5))') 'idat,iq,qeff,ein,energy=',idat,iq,qeff,ein(idat),energy(1,iset)
      STOP 'too large effective q'
    END IF
    q1=iq*delql-delql
    q2=q1*delql
    ffdat(idat,:,iset)=((qeff-q1)*formfa(iq+1,:,iset)+(q2-qeff)*formfa(iq,:,iset))/delql
    ffdat0(idat,:) = ffdat(idat,:,iset)**2
    sigma = faclong*ffdat0(idat,1)+factrans*ffdat0(idat,2)*(1D0+TAN(angle(idat)*pi/360D0)**2)
    sigmalog(idat) = LOG(sigma)
    WRITE(7,'(5(1pg13.5))') angle(idat),ein(idat),sigma
  END DO
  WRITE(7,'(1x/1x)')
ELSE
  DO idat=1,kdat
    sigma = faclong*ffdat0(idat,1)+factrans*ffdat0(idat,2)*(1D0+TAN(angle(idat)*pi/360D0)**2)
    sigmalog(idat) = LOG(sigma)
  END DO
END IF

  aver0=kdat
  averx=SUM(angle(1:kdat))/aver0
  averx2=SUM(angle(1:kdat)**2)/aver0
  averf=SUM(sigmalog(1:kdat))/aver0
  averfx=SUM(sigmalog(1:kdat)*angle(1:kdat))/aver0
  slopeout=(averfx-averf*averx)/(averx2-averx**2)
  averout=(averf*averx2-averfx*averx)/(averx2-averx**2)

END SUBROUTINE slopes

END PROGRAM analyze_transdens

!-----bessel------------------------------------------------------

FUNCTION bessel(x,l)
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'params.inc'

!     spherical bessel function of angular momentum 'l' at
!     argument 'x'.

!     uses series expansion for x<l=1 and asymptotic expression
!     for x>l-1.

!     precision is tested only up to 5 digits.

!INCLUDE 'constn.inc'



REAL(8), INTENT(IN)                      :: x
INTEGER, INTENT(IN)                      :: l
REAL(8), PARAMETER :: xargmx=1D10


!------------------------------------------------------------------------------

!c       if(l.gt.7) stop ' l>7 in "bessel"'

!      case of zero argument

IF(x == 0D0) THEN
  IF(l == 0) THEN
    b      = 1D0
  ELSE
    b      = 0D0
  END IF
  
!      very large argument  --> set result to 0D0
  
ELSE IF(ABS(x) > xargmx) THEN
  b      = 0D0
  
!      small arguments  --> use series expansion
  
ELSE IF(ABS(x) <= l-1D0) THEN
  add   = 1D0
  DO i=1,l
    add   = x/(2*i+1)*add
  END DO
  arg2  = -0.25*x*x
  b     = add
  DO i=1,l+1
    add   = arg2*add/(i*(l+i+0.5D0))
    b     = add + b
  END DO
  
!      large arguments  --> use expansion in sin and cos
  
ELSE
  s      = SIN(x)/x
  c      = COS(x)/x
  IF(l == 0) THEN
    b      = s
  ELSE IF(l == 1) THEN
    b      = s/x-c
  ELSE IF(l >= 2) THEN
    xinv   = 1D0/x
    bjm    = s
    b      = s/x-c
    DO lrec=2,l
      bjp    = ((lrec+lrec)*xinv-xinv)*b-bjm
      bjm    = b
      b      = bjp
    END DO
    b      = b
  END IF
END IF

bessel = b
RETURN
END FUNCTION bessel
