!-----evaluate_zarray------------------------------------------------
 

PROGRAM evaluate_zarray

! analyze zarray.dat

!IMPLICIT DOUBLE PRECISION (a-h,o-z)
IMPLICIT NONE

INTEGER,PARAMETER :: kparmx=20
!INCLUDE 'comfit.inc'
REAL(8) :: ZARRAY(kparmx,kparmx),Znew(kparmx,kparmx),DPARMO(kparmx)
REAL(8) :: par(kparmx),workdv(0:kparmx,3)
CHARACTER (LEN=3) :: ident

INTEGER :: nparmx,i,j,nerr
REAL(8) :: det
REAL(8) :: barrintsav,ehsdirsav,ehssav,errdir,errehs

!---------------------------------------------------------------------

!     read 'dparmo' and 'zarray' from data file

OPEN(31,FILE='zarray.dat')
READ(31,'(i5)') nparmx
!WRITE(6,*) ' NPARMX=',nparmx
READ(31,'(17g24.16)') (dparmo(j),j=1,nparmx)
!WRITE(6,*) 'DPARMO:',(dparmo(j),j=1,nparmx)
DO i=1,nparmx
  READ(31,'(17g24.16)') (zarray(j,i),j=1,nparmx)
END DO
CLOSE (31)

DO i=1,nparmx; DO j=1,nparmx; 
  znew(j,i) = zarray(j,i)/SQRT(zarray(j,j)*zarray(i,i))
END DO; END DO;

call matin2(znew,nparmx,kparmx,nerr,det)

WRITE(6,*) ' Determinant=',det

STOP
END PROGRAM evaluate_zarray
