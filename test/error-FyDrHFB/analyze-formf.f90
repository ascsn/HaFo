IMPLICIT REAL(8) (a-h,o-z)
INTEGER,PARAMETER :: ndatmx=13,ngrid=61
REAL(8) :: dat(1:ngrid,0:ndatmx,1:4)  ! last index stands for type of density
REAL(8) :: q(1:ngrid),acc(1:4),covar(1:ngrid,1:ngrid,1:4),var(1:ngrid,1:4)
CHARACTER(20) :: acase

OPEN(1,file='collect-error-formf.dat')
READ(1,'(a)') acase
WRITE(6,'(a)') '#'//acase
DO i=0,13
  READ(1,'(1x)')
  DO j=1,ngrid
    READ(1,*) jdum,q(j),dat(j,i,:)
  END DO
END DO


! evaluate and print uncertainities point by point
WRITE(6,'(a)') '#   q   formf-coul   error     formf-prot   error     formf-neut   error   formf-mass  error'
DO j=1,ngrid
  acc=0D0
  DO i=1,ndatmx
    dat(j,i,:) = dat(j,i,:)-dat(j,0,:)
    acc = acc + dat(j,i,:)**2
  END DO
  WRITE(6,'(f7.1,8(1pg13.5))') q(j),(dat(j,0,icase),SQRT(acc(icase))*10D0,icase=1,4)
  var(j,:) = acc
END DO
  WRITE(6,'(1x/1x)')

! evaluate and print uncertainities point by point
var=SQRT(var)
WRITE(6,'(a)') '#   correlation matrix', &
    '#  r_1  r_2    covariances coul-prot-neut-mass'
DO j1=1,ngrid
  DO j2=1,ngrid;
    covar(j2,j1,:)=0D0
    DO i=1,ndatmx
      covar(j2,j1,:) = covar(j2,j1,:) + dat(j1,i,:)*dat(j2,i,:)
    END DO
    WRITE(6,'(2f5.1,4(1pg13.5))') q(j1),q(j2),covar(j2,j1,:)/var(j2,:)/var(j1,:)
  END DO
  WRITE(6,'(1x)')
END DO


STOP
END
