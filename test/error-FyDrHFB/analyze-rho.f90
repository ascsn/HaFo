IMPLICIT REAL(8) (a-h,o-z)
INTEGER,PARAMETER :: ndatmx=13,ngrid=61
REAL(8) :: dat(1:ngrid,0:ndatmx,1:4)  ! last index stands for type of density
REAL(8) :: r(1:ngrid),acc(1:4)
CHARACTER(20) :: acase

OPEN(1,file='collect-error-rho.dat')
READ(1,'(a)') acase
WRITE(6,'(a)') '#'//acase
DO i=0,13
  READ(1,'(1x)')
  DO j=1,ngrid
    READ(1,*) jdum,r(j),dat(j,i,:)
  END DO
END DO

WRITE(6,'(a)') '#   r   rho-coul   error     rho-prot   error      rho-neut   error    rho-mass  error'
DO j=1,ngrid
  acc=0D0
  DO i=1,ndatmx
    acc = acc + (dat(j,i,:)-dat(j,0,:))**2
  END DO
  WRITE(6,'(f7.1,8(1pg13.5))') r(j),(dat(j,0,icase),SQRT(acc(icase))*10D0,icase=1,4)
END DO

STOP
END
