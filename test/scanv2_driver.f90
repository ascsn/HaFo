CHARACTER(10) :: xx
!DO i=0,300
!  WRITE(xx,'(f6.4)') i*0.0001
!DO i=0,300
!DO i=301,500
DO i=0,400
  WRITE(xx,'(f10.8)') (i*0.0002)**2
!  WRITE(xx,'(f10.8)') (i*0.0005)**2
  CALL SYSTEM('./scanv2.bat '//xx)
END DO
CALL SYSTEM('echo "  " >> collect.v2-space2')
STOP
END

