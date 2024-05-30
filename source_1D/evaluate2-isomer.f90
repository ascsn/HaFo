!-----evaluate_barriers------------------------------------------------
 

PROGRAM evaluate_barriers

!   reads data from 'ensurf' for symmetric path
!   computes g.s., first barrier and isomer

!IMPLICIT DOUBLE PRECISION (a-h,o-z)
IMPLICIT NONE

INTEGER,PARAMETER :: kparmx=14
!INCLUDE 'comfit.inc'
REAL(8) :: ZARRAY(kparmx,kparmx),DPARMO(kparmx)
REAL(8) :: par(kparmx),workdv(0:kparmx,4)
CHARACTER (LEN=3) :: ident,id2

INTEGER :: nparmx,i,j
REAL(8) :: barr1,engs,barr1sav,engssav,eniso,enisosav,barr2,barr2sav
REAL(8) :: errengs,erreniso,errbarr1,errbarr2,parmnatfac

!---------------------------------------------------------------------

!     read 'dparmo' and 'zarray' from data file

OPEN(31,FILE='zarray.dat')
READ(31,'(i5)') nparmx
CLOSE (31)
READ(5,'(a3)') id2
READ(id2,'(i3)') i
parmnatfac = i/1D2


!nparmx = 0
WRITE(6,*) ' NPARMX=',nparmx

!     read barriers from 'ensurf.*'
DO i=0,nparmx
  
!       prepare identifier
  
  OPEN(31,STATUS='scratch')
  WRITE(31,'(i3)') i+100
  REWIND(31)
  READ(31,'(a3)') ident
  CLOSE(31)
  ident(1:1) = 's'
  
!       read and store
  
  OPEN(31,FILE='ensurf.'//id2//'.'//ident)
  WRITE(6,*) ' OPENED:','ensurf.'//id2//'.'//ident
!test        write(6,'(a,i5)')
!test     &   'ensurf.'//ident,i
  IF(i == 0) THEN
    CALL read_ensurf(engssav,barr1sav,enisosav,31)
    WRITE(6,*) ' i,barrdirsav,ehssav=',i,engssav,barr1sav,enisosav
    workdv(0,1) = engssav
    workdv(0,2) = barr1sav
    workdv(0,3) = enisosav
  ELSE
    CALL read_ensurf(engs,barr1,eniso,31)
    workdv(i,1) = (engs-engssav)
    workdv(i,2) = (barr1-barr1sav)
    workdv(i,3) = (eniso-enisosav)
    WRITE(6,*) ' i,barrdir,ehs=',i,engs,barr1,eniso
  END IF
  CLOSE(31)
  
END DO
!  output for erroranalysis
OPEN(7,file='Pu-isomer.cmp.'//id2)
WRITE(7,'(a)') ' E_B 240Pu'
WRITE(7,'(1pg18.10)') workdv(0,1)
WRITE(7,'(50(1pg18.10))') workdv(1:nparmx,1)
WRITE(7,'(a)') ' 1.barrier 240Pu'
WRITE(7,'(1pg18.10)') workdv(0,2)
WRITE(7,'(50(1pg18.10))') workdv(1:nparmx,2)
WRITE(7,'(a)') ' isomer 240Pu'
WRITE(7,'(1pg18.10)') workdv(0,3)
WRITE(7,'(50(1pg18.10))') workdv(1:nparmx,3)
CLOSE(7)

!     compute variances

errengs   = 0.0
errbarr1   = 0.0
erreniso   = 0.0

IF(nparmx > 0) THEN
  DO i=1,nparmx
    errengs  = workdv(i,1)*workdv(i,1) + errengs
    errbarr1  = workdv(i,2)*workdv(i,2) + errbarr1
    erreniso  = workdv(i,3)*workdv(i,3) + erreniso
  END DO
  errengs  = SQRT(errengs)/parmnatfac
  errbarr1  = sqrt(errbarr1)/parmnatfac
  erreniso  = SQRT(erreniso)/parmnatfac
END IF

!     print

WRITE(6,'(a,2(1pg13.5)/a,2(1pg13.5)/a,2(1pg13.5)/a,2(1pg13.5))')  &
    ' binding energy:       ',engssav,errengs,    &
    ' 1. fission barrier    ',barr1sav,errbarr1,  &
    ' isomer energy:        ',enisosav,erreniso  

STOP
END PROGRAM evaluate_barriers

!-----read_ensurf----------------------------------------------------

SUBROUTINE read_ensurf(engs,barr1,eniso,iunit)

! Reads information from file and extracts fission barrier.
!       iunit      = unit number for file with input data
!       barrdir    = barrier from maximum of listed values
!       barrint    = barrier from interpolated maximum
! Input file is ssumed to be sorted in order of decreasing 'beta'.


!IMPLICIT REAL(8) (a-h,o-z)
IMPLICIT NONE
REAL(8), INTENT(OUT)    :: barr1
REAL(8), INTENT(OUT)    :: engs,eniso
INTEGER, INTENT(IN OUT) :: iunit

INTEGER, PARAMETER :: kmaxen=200

REAL(8) :: e(kmaxen)
LOGICAL :: tsmin
LOGICAL,PARAMETER :: tinterpol=.true.
INTEGER :: i,imax(2),imin(2),nmax,nmin
REAL(8) :: emax(2),emin(2),bmin,bmax,d1e,d2e,x
REAL(8) :: beta(6),bsav
REAL(8) :: beta2max=1.5   ! search is done only for smaller deformations


!-------------------------------------------------------------------

DO i =1,10
  READ(iunit,'(1x)')
!        read(iunit,'(a)') a
!        write(6,'(a)') a
END DO
emax  = -1.0D33
emin  = 1.0D33
imax  = 0
imin  = 0
tsmin = .true.
nmax=1
nmin=1
DO i=1,kmaxen
  READ(iunit,*,ERR=19,END=19) beta,e(i)
!  WRITE(6,'(6f9.3,1pg13.5)') beta,e(i)
  IF(i.GT.1 .AND. ABS(beta(1)-bsav).GT.0.1D0) CYCLE
  IF(i.GT.1 .AND. beta(1).LT.bsav) CYCLE
  bsav = beta(1)
  IF(beta(1)>beta2max) CYCLE
  IF(tsmin) THEN
    IF(e(i) < emin(nmin)) THEN
      emin(nmin) = e(i)
      imin(nmin) = i
      bmin = beta(1)
    ELSE
      tsmin = .false.
      nmin=1+nmin
      WRITE(*,*) ' switch from minimum: i=',i
!      GO TO 19
    END IF
  ELSE
!test        write(6,'(i,1pg13.5)') i,e(i)
    IF(e(i) > emax(nmax)) THEN
      emax(nmax) = e(i)
      imax(nmax) = i
      bmax = beta(1)
    ELSE
      tsmin = .true.
      nmax=1+nmax
      IF(nmax==3) GOTO 19
      WRITE(*,*) ' switch from maximum: i=',i
    END IF
  END IF
END DO
STOP ' KMAXEN not large enough'

19   CONTINUE
engs = emin(1)
eniso = emin(2)-engs
barr1 = emax(1)-engs


IF(tinterpol) THEN

  ! interpolate to minimum energy
  IF(imin(1)<2) STOP " IMIN(1)<2" 
  d1e     = (e(imin(1)+1)-e(imin(1)-1))/2.0D0
  d2e     = e(imin(1)+1)-2.0D0*e(imin(1))+e(imin(1)-1)
  x       = -d1e/d2e
  engs    = e(imin(1))+x*d1e+x*x*d2e/2.0D0

  IF(imin(2)<2) STOP " IMIN(2)<2" 
  d1e     = (e(imin(2)+1)-e(imin(2)-1))/2.0D0
  d2e     = e(imin(2)+1)-2.0D0*e(imin(2))+e(imin(2)-1)
  x       = -d1e/d2e
  eniso   = e(imin(2))+x*d1e+x*x*d2e/2.0D0-engs

  ! interpolate to maximum energy
  IF(imax(1)>=kmaxen) STOP "IMAX(1) too large"
  d1e     = (e(imax(1)+1)-e(imax(1)-1))/2.0D0
  d2e     = e(imax(1)+1)-2.0D0*e(imax(1))+e(imax(1)-1)
  x       = -d1e/d2e
  barr1   = e(imax(1))+x*d1e+x*x*d2e/2.0D0 -engs

!WRITE(6,'(a,1pg13.5)') &
!  ' g.s.: ',engs-e(imin(1)), &
!  ' barr1:',barr1-e(imax(1))+e(imin(1)), &
!  ' iso:  ',eniso-e(imin(2)), &
!  ' barr2:',barr2-e(imax(2))+e(imin(1))

END IF

RETURN

! 99   stop ' premature end of input file'

END SUBROUTINE read_ensurf
