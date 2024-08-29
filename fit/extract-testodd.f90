REAL :: expiso(999),thiso(999),expgap(999),thgap(999),exprms(999),thrms(999)
INTEGER :: naiso1(999),nziso1(999),naiso2(999),nziso2(999), &
           nagap(999),nzgap(999),narms(999),nzrms(999)
CHARACTER(90) :: line

niso = 0
ngap = 0
nrms = 0

OPEN(1,file='for002',action='READ')

DO i=1,99999

  READ(1,'(a)',err=99,end=29) line

  IF(index(line,"isotope sh")>0) THEN
    niso = niso+1
!    WRITE(6,'(a)') line
    BACKSPACE(1)
    READ(1,'(17x,i4,1x,i4,3x,i4,1x,i4,14x,2e12.4)') &
      naiso1(niso),nziso1(niso),naiso2(niso),nziso2(niso),expiso(niso),thiso(niso)
  END IF

  IF(index(line,"3gap")>0) THEN
    ngap = ngap+1
    BACKSPACE(1)
    READ(1,'(6x,i4,1x,i4,14x,2e12.4)') &
      nagap(ngap),nzgap(ngap),expgap(ngap),thgap(ngap)
  END IF

  IF(index(line,"rms radius")>0 .AND. nrms==0) THEN
    DO j=1,9999
      READ(1,'(a)') line
!      WRITE(6,'(a)') line
      IF(index(line,"rms error")>0) GOTO 19
      BACKSPACE(1)
      READ(1,*,err=19,end=99) narms(j),nzrms(j),exprms(j),thrms(j)
      nrms=j
    END DO
    STOP 'rms loop exhausted'
  19 CONTINUE
  END IF

END DO

29 CONTINUE

!WRITE(6,'(a,3i4)') ' niso,ngap,nrms=',niso,ngap,nrms

WRITE(6,'(a)') '#isotopes shifts:'
DO j=1,niso
  IF(j>1 .AND. ABS(nzsav-nziso1(j))>2) WRITE(6,'(2x)')
  WRITE(6,'(4i4,2(1pg12.4))') naiso1(j),nziso1(j),naiso2(j),nziso2(j),expiso(j),thiso(j)
  nzsav= nziso1(j)
END DO

WRITE(6,'(/1x/a)') '#3gaps:'
DO j=1,ngap
  IF(j>1 .AND. ABS(nasav-nagap(j))>2) WRITE(6,'(2x)')
  WRITE(6,'(2i4,2(1pg12.4))') nagap(j),nzgap(j),expgap(j),thgap(j)
  IF(nagap(j)==145) WRITE(6,'(2x)')
  nasav = nagap(j)
END DO


WRITE(6,'(/1x/a)') '#rms radii:'
DO j=1,nrms
  IF(j>1 .AND. ABS(nasav-narms(j))>2) WRITE(6,'(2x)')
  WRITE(6,'(2i4,2(1pg13.5))') narms(j),nzrms(j),exprms(j),thrms(j)
  nasav = narms(j)
END DO

STOP 'regular end'
99 STOP 'error or premature end'
END
