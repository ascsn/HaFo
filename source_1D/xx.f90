SUBROUTINE find_surfe(formfqs,q,rdms,xmin,xmax,x,helmpower)
  REAL(8), INTENT(IN) :: xmin,xmax,helmpower
  REAL(8), INTENT(OUT) :: x
  REAL(8), INTENT(IN) :: formfqs,q,rdms
  
  REAL(8) :: f1,f2,f3,x1,x2,x3
  INTEGER :: i

  INTEGER, PARAMETER :: nbracket=40
  REAL(8), PARAMETER :: epsil=1D-11
  LOGICAL,PARAMETER :: tprint=.FALSE.

  x1=xmin
  x2=xmax
  x3=x1
  IF(tprint) THEN
    x1=SQRT(0.6D0)*q*x3
    WRITE(311,'(6(1pg13.5))') x3,formfqs,formf_cond(x3),&
     3D0*(SIN(q*rdms)-q*rdms*COS(q*rdms))*(q*rdms)**(-3),&
     SINH(x1)/(x1*COSH(x1)**2)
    x3=x2
    x1=SQRT(0.6D0)*q*x3
    WRITE(311,'(6(1pg13.5))') x3,formfqs,formf_cond(x3),&
     3D0*(SIN(q*rdms)-q*rdms*COS(q*rdms))*(q*rdms)**(-3),&
     SINH(x1)/(x1*COSH(x1)**2)
    x1=xmin
  END IF
  f1=formf_cond(x1)
  f2=formf_cond(x2)
  IF(f1*f2>0D0) STOP 'initially not brackted'
  DO i=1,nbracket
    IF(ABS(x2-x1)<epsil) EXIT
    x3 = 0.5D0*(x1+x2)
    f3 = formf_cond(x3)
    IF(tprint) WRITE(311,'(i3,6(1pg13.5))') i,x1,x2,x3,f1,f2,f3
    IF(f3==0D0) EXIT
    IF(f1*f3<0D0) THEN
      x2=x3 
    ELSE IF(f2*f3<0D0) THEN
      x1=x3
    END IF
    f3=formf_cond(x3)
  END DO
  
  x=x3  

  IF(tprint) THEN
    x1=SQRT(0.6D0)*q*x3
    WRITE(311,'(6(1pg13.5))') x,formfqs,f3,&
     3D0*(SIN(q*rdms)-q*rdms*COS(q*rdms))*(q*rdms)**(-3),&
     SINH(x1)/(x1*COSH(x1)**2)
  END IF

CONTAINS
REAL(8) FUNCTION formf_cond(s)
  REAL(8), INTENT(IN) :: s
  REAL(8) :: qr,qs
  REAL(8), PARAMETER :: sq35=SQRT(0.6D0)

  qr = q*rdms
  IF(helmpower==0D0) THEN
    qs = q*s*sq35
    formf_cond = 3D0*(SIN(qr)-qr*COS(qr))*qr**(-3) * &
                 SINH(qs)/(qs*COSH(qs)**2)  -  formfqs
  ELSE
    qs = q*s/SQRT(helmpower)  
    formf_cond = 3D0*(SIN(qr)-qr*COS(qr))*qr**(-3) * &
                 COSH(qs)**(-helmpower)  -  formfqs
  END IF
  
  RETURN

END FUNCTION formf_cond

END SUBROUTINE find_surfe
