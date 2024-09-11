IMPLICIT REAL(8) (a-h,o-z)
INTEGER,PARAMETER :: ndatmx=13
REAL(8) :: dat(0:ndatmx)
OPEN(1,file='collect-error.dat')
READ(1,'(1x)')
DO i=0,13
  READ(1,*) dat(i)
END DO
acc=0D0
!WRITE(6,'(i3,4(1pg13.5))') 0,dat(0)
DO i=1,ndatmx
  acc = acc + (dat(i)-dat(0))**2
!  WRITE(6,'(i3,4(1pg13.5))') i,dat(i),(dat(i)-dat(0))*10D0,(dat(i)-dat(0))**2*100D0,acc*100D0
END DO
WRITE(6,'(a,1pg13.5)') 'extrapolation error=',SQRT(acc)*10D0
!WRITE(6,'(a,1pg13.5)') 'extrapolation error=',SQRT(SUM((dat(1:ndatmx)-dat(0))**2))*10D0
STOP
END
