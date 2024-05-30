!-----evaluate_barriers------------------------------------------------
 

PROGRAM evaluate_barriers

!     version with extracts data from one CHF path

!IMPLICIT DOUBLE PRECISION (a-h,o-z)
IMPLICIT NONE

INTEGER,PARAMETER :: kparmx=14
!INCLUDE 'comfit.inc'
REAL(8) :: ZARRAY(kparmx,kparmx),DPARMO(kparmx)
REAL(8) :: par(kparmx),workdv(kparmx,3)
CHARACTER (LEN=3) :: ident

INTEGER :: nparmx,i,j
REAL(8) :: barrdir,ehs,barrdirsav,ehssav,errdir,errehs

!---------------------------------------------------------------------

!     read 'dparmo' and 'zarray' from data file

OPEN(31,FILE='zarray.dat')
READ(31,'(i5)') nparmx
READ(31,'(17g24.16)') (dparmo(j),j=1,nparmx)
DO i=1,nparmx
  READ(31,'(17g24.16)') (zarray(j,i),j=1,nparmx)
END DO
CLOSE (31)

!nparmx = 0
WRITE(6,*) ' NPARMX=',nparmx
WRITE(6,*) ' DPARMO=',(dparmo(j),j=1,nparmx)

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
  
  OPEN(31,FILE='ensurf.'//ident)
  WRITE(6,*) ' OPENED:','ensurf.'//ident
!test        write(6,'(a,i5)')
!test     &   'ensurf.'//ident,i
  IF(i == 0) THEN
    CALL read_ensurf(barrdirsav,ehssav,31)
    WRITE(6,*) ' i,barrdirsav,ehssav=',i,barrdirsav,ehssav
  ELSE
    CALL read_ensurf(barrdir,ehs,31)
    workdv(i,1) = (barrdir-barrdirsav)/dparmo(i)
!          workdv(i,2) = (barrint-barrintsav)/dparmo(i)
    workdv(i,3) = (ehs-ehssav)/dparmo(i)
    WRITE(6,*) ' i,barrdir,ehs=',i,barrdir,ehs
  END IF
  CLOSE(31)
  
END DO

!     compute variances

errdir   = 0.0
!      errint   = 0.0
errehs   = 0.0

IF(nparmx > 0) THEN
  DO i=1,nparmx
    DO j=1,nparmx
      errdir  = workdv(i,1)*workdv(j,1)*zarray(i,j) + errdir
!          errint  = workdv(i,2)*workdv(j,2)*zarray(i,j)
!     &              + errint
      errehs  = workdv(i,3)*workdv(j,3)*zarray(i,j) + errehs
    END DO
  END DO
  errdir  = SQRT(errdir)
!        errint  = sqrt(errint)
  errehs  = SQRT(errehs)
END IF

!     print

WRITE(6,'(a,2(1pg13.5)/a,2(1pg13.5)/a,2(1pg13.5))')  &
    ' fission barrier - direct:',barrdirsav,errdir,  &
    ' binding energy:          ',ehssav,errehs

STOP
END PROGRAM evaluate_barriers

!-----read_ensurf----------------------------------------------------

SUBROUTINE read_ensurf(barrdir,energ0,iunit)

! Reads information from file and extracts fission barrier.
!       iunit      = unit number for file with input data
!       barrdir    = barrier from maximum of listed values
!       barrint    = barrier from interpolated maximum
! Input file is ssumed to be sorted in order of decreasing 'beta'.


!IMPLICIT REAL(8) (a-h,o-z)
IMPLICIT NONE
REAL(8), INTENT(OUT)    :: barrdir
REAL(8), INTENT(OUT)    :: energ0
INTEGER, INTENT(IN OUT) :: iunit

INTEGER, PARAMETER :: kmaxen=200

REAL(8) :: e(kmaxen)
LOGICAL :: tsmin
INTEGER :: i,imax,imin
REAL(8) :: emax,emin,bmin,bmax
REAL(8) :: beta(12)
REAL(8) :: beta2max=0.9   ! search is done only for smaller deformations


!-------------------------------------------------------------------

DO i =1,9
  READ(iunit,'(1x)')
!        read(iunit,'(a)') a
!        write(6,'(a)') a
END DO
emax  = -1.0D33
emin  = 1.0D33
imax  = 0
imin  = 0
tsmin = .false.
DO i=1,kmaxen
!old        read(iunit,*,err=19,end=19) x1,x2,x3,x4,x5,x6,x7,e(i)
  READ(iunit,*,ERR=19,END=19) beta,e(i)
  IF(beta(4)>beta2max) CYCLE
  IF(tsmin) THEN
    IF(e(i) < emin) THEN
      emin = e(i)
      imin = i
      bmin = beta(4)
    ELSE
!      tsmin = .false.
      GO TO 19
    END IF
  ELSE
!test        write(6,'(i,1pg13.5)') i,e(i)
    IF(e(i) > emax) THEN
      emax = e(i)
      imax = i
      bmax = beta(4)
    ELSE
      tsmin = .true.
!      GO TO 19
    END IF
  END IF
END DO
STOP ' KMAXEN not large enough'

19   CONTINUE
energ0 = emin
barrdir = emax!-energ0

!WRITE(*,*) ' emin,emax,bmin,bmax,imin,imax=',emin,emax,bmin,bmax,imin,imax

RETURN

! 99   stop ' premature end of input file'

END SUBROUTINE read_ensurf
