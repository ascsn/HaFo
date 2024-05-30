!-----evaluate_surface------------------------------------------------
 

PROGRAM evaluate_surface

!     version with extracts data from one CHF path

!IMPLICIT DOUBLE PRECISION (a-h,o-z)
IMPLICIT NONE

INTEGER,PARAMETER :: kparmx=14
!INCLUDE 'comfit.inc'
REAL(8) :: parmnatfac
REAL(8) :: ZARRAY(kparmx,kparmx),DPARMO(kparmx),as(kparmx)
REAL(8) :: par(kparmx),workdv(0:kparmx,3)
CHARACTER (LEN=3) :: ident

INTEGER :: nparmx,i,j
REAL(8) :: barrint,barrdir,ehs,ehsdir,barrdirsav,asurf,asurfi
REAL(8) :: barrintsav,ehsdirsav,ehssav,errdir,errehs,dum(3)

!---------------------------------------------------------------------

!     read 'dparmo' and 'zarray' from data file

READ(5,*) parmnatfac
parmnatfac = parmnatfac/1D2
OPEN(31,FILE='zarray.dat')
READ(31,'(i5)') nparmx
WRITE(6,*) ' NPARMX=',nparmx
zarray = 0D0
DO i=1,nparmx
  zarray(i,i) = 1D0/parmnatfac
END DO
CLOSE (31)


!nparmx = 0
WRITE(6,*) ' NPARMX,PARMNATFAC=',nparmx,parmnatfac
!WRITE(6,*) ' DPARMO=',(dparmo(j),j=1,nparmx)


OPEN(31,FILE='asurf.collect')
OPEN(32,FILE='asurf.cmp')
WRITE(32,'(a)') ' surface energy'
READ(31,'(1x)')
READ(31,*) dum,asurf
WRITE(32,'(1pg18.10)') asurf
DO i=1,nparmx
 WRITE(6,*) ' i=',i
 READ(31,'(1x)')
READ(31,'(1x)')
 READ(31,*) dum,asurfi
 as(i) = (asurfi-asurf)
END DO
 WRITE(32,'(100(1pg18.10))') as(1:nparmx)
CLOSE(31)

STOP

END PROGRAM evaluate_surface
