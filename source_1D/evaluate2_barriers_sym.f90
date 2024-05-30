!-----evaluate_barriers------------------------------------------------
 

PROGRAM evaluate_barriers

!   version with extracts data from one CHF path
!   uses natural parameter variation

!IMPLICIT DOUBLE PRECISION (a-h,o-z)
IMPLICIT NONE

INTEGER,PARAMETER :: kparmx=14
!INCLUDE 'comfit.inc'

REAL(8),PARAMETER :: parmnatfac=0.5D0  ! step factor for natural params
REAL(8) :: ZARRAY(kparmx,kparmx),DPARMO(kparmx)
REAL(8) :: par(kparmx),workdv(0:kparmx,3)
CHARACTER (LEN=3) :: ident,id2

INTEGER :: nparmx,i,j
REAL(8) :: barrint,barrdir,ehs,ehsdir,barrdirsav
REAL(8) :: barrintsav,ehsdirsav,ehssav,errdir,errehs

!---------------------------------------------------------------------


OPEN(31,FILE='zarray.dat')
READ(31,'(i5)') nparmx
CLOSE (31)
!nparmx = 0
WRITE(6,*) ' NPARMX=',nparmx
!WRITE(6,*) ' DPARMO=',(dparmo(j),j=1,nparmx)
READ(5,'(a3)') id2

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
    CALL read_ensurf(barrdirsav,ehsdirsav,barrintsav,ehssav,31)
    workdv(0,1) = barrintsav
    workdv(0,3) = ehssav
    WRITE(6,*) ' i,barrdirsav,ehsdirsav=',i,barrdirsav,ehsdirsav
    WRITE(6,*) ' i,barrintsav,ehssav=',i,barrintsav,ehssav
  ELSE
    CALL read_ensurf(barrdir,ehsdir,barrint,ehs,31)
    workdv(i,1) = (barrint-barrintsav)
!          workdv(i,2) = (barrint-barrintsav)/dparmo(i)
    workdv(i,3) = (ehs-ehssav)
    WRITE(6,*) ' i,barrdir,ehsdir=',i,barrdir,ehsdir
    WRITE(6,*) ' i,barrint,ehs=',i,barrint,ehs
  END IF
  CLOSE(31)
  
END DO

!  output for erroranalysis
OPEN(7,file='HS-barrier.cmp.'//id2)
WRITE(7,'(a)') ' E_B 266Hs'
WRITE(7,'(1pg18.10)') workdv(0,3)
WRITE(7,'(50(1pg18.10))') workdv(1:nparmx,3)
WRITE(7,'(a)') ' barrier 266Hs'
WRITE(7,'(1pg18.10)') workdv(0,1)
WRITE(7,'(50(1pg18.10))') workdv(1:nparmx,1)
CLOSE(7)

!     compute variances

errdir   = 0.0
errehs   = 0.0

IF(nparmx > 0) THEN
  DO i=1,nparmx
      errdir  = workdv(i,1)*workdv(i,1) + errdir
      errehs  = workdv(i,3)*workdv(i,3) + errehs
  END DO
  errdir  = SQRT(errdir)/parmnatfac
  errehs  = SQRT(errehs)/parmnatfac
END IF

!     print

WRITE(6,'(a,2(1pg13.5)/a,2(1pg13.5)/a,2(1pg13.5))')  &
    ' fission barrier in 266Hs - direct:',barrdirsav,errdir,  &
    ' binding energy in 266Hs:            ',ehssav,errehs

STOP
END PROGRAM evaluate_barriers

!-----read_ensurf----------------------------------------------------

SUBROUTINE read_ensurf(barrdir,energdir,barrint,energ0,iunit)

!     reads information from file and
!     extracts fission barrier.
!       iunit      = unit number for file with input data
!       barrdir    = barrier from maximum of listed values
!       barrint    = barrier from interpolated maximum


!IMPLICIT REAL(8) (a-h,o-z)
IMPLICIT NONE
REAL(8), INTENT(OUT)    :: barrdir,barrint
REAL(8), INTENT(OUT)    :: energ0,energdir
INTEGER, INTENT(IN OUT) :: iunit

INTEGER, PARAMETER :: kmaxen=100

REAL(8) :: e(kmaxen)
LOGICAL :: tsmin
INTEGER :: i,imax,imin,maxen,minen
REAL(8) :: emax,emin,d1e,d2e,x
REAL(8) :: dummy(6)


!-------------------------------------------------------------------

DO i =1,12
  READ(iunit,'(1x)')
!        read(iunit,'(a)') a
!        write(6,'(a)') a
END DO
maxen = 0
emax  = -1.0D33
emin  = 1.0D33
imax  = 0
imin  = 0
tsmin = .true.
DO i=1,kmaxen
!old        read(iunit,*,err=19,end=19) x1,x2,x3,x4,x5,x6,x7,e(i)
  READ(iunit,*,ERR=19,END=19) dummy,e(i)
  maxen = i
  IF(tsmin) THEN
    minen = i
    IF(e(i) < emin) THEN
      emin = e(i)
      imin = i
    ELSE
      tsmin = .false.
    END IF
  ELSE
!test        write(6,'(i,1pg13.5)') i,e(i)
    IF(e(i) > emax) THEN
      emax = e(i)
      imax = i
    END IF
  END IF
END DO
STOP ' KMAXEN not large enough'

19   CONTINUE
energdir = emin
barrdir = emax-energdir

! interpolate to minimum energy
IF(imin<2) STOP " IMIN<2" 
d1e     = (e(imin+1)-e(imin-1))/2.0D0
d2e     = e(imin+1)-2.0D0*e(imin)+e(imin-1)
x       = -d1e/d2e
energ0  = e(imin)+x*d1e+x*x*d2e/2.0D0

! interpolate to maximum energy
IF(imax>=maxen) STOP "IMAX too large"
d1e     = (e(imax+1)-e(imax-1))/2.0D0
d2e     = e(imax+1)-2.0D0*e(imax)+e(imax-1)
x       = -d1e/d2e
barrint = e(imax)+x*d1e+x*x*d2e/2.0D0-energ0




RETURN

! 99   stop ' premature end of input file'

END SUBROUTINE read_ensurf
