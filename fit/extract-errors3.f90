PROGRAM extract_errors
IMPLICIT NONE
!
! To extract key values from fit output 'for002', all for E,R,sigma,r:
!   average errors (absolute)
!   r.m.s. errors (absolute)
!   sub-chi**2
!   errors for 208Pb
!
  CHARACTER(88) :: aline,aold1,aold2
  INTEGER, PARAMETER :: maxline=99999

  REAL(8) :: error(4,4) ! storage for errors in fixed order as in 'for002'
!  REAL(8) :: energ(39:52),rms(39:52)  ! storage for energy and rm radius
  REAL(8) :: dummy(8),chi2=0D0,chi2b=0D0,chi2c=0D0,datin
  INTEGER :: iread,icount=0,ib,idum(2),na,nz

  OPEN(2,file='for002.SVmin-softcut')
  OPEN(3,file='errors.SVmin-softcut')
!  energ=0D0
!  rms=0D0
   DO iread=1,maxline
     aold2 = aold1
     aold1 = aline
     READ(2,'(a)',ERR=19,END=19) aline
     IF (aline(1:11)=='   binding-') THEN
       DO ib=1,99
         READ(2,*,ERR=29) na,nz,dummy(1:2),datin
         IF(nz==82 .AND. na==208) THEN
           error(4,1) = datin
         ELSE IF(nz==92 .AND. na==218) THEN
           backspace(2)
           READ(2,*) na,nz,dummy(1:6),error(3,1)
         END IF
       END DO
29     CONTINUE
     ELSE IF (aline(1:8)=='   diffr') THEN
       DO ib=1,99
         READ(2,*,ERR=39) na,nz,dummy(1:2),datin
         IF(nz==82 .AND. na==208) THEN
           error(4,2) = dummy(1)-dummy(2)
           print*, "Error diff:", nz, na, error(4, 2)
         ELSE IF(nz==60 .AND. na==142) THEN
           backspace(2)
           READ(2,*) na,nz,dummy(1:6),error(3,2)
           print*, "diff:", nz, na, error(3, 2)
         END IF
       END DO
39     CONTINUE
     ELSE IF (aline(1:8)=='   surfa') THEN
       DO ib=1,99
         READ(2,*,ERR=49) na,nz,dummy(1:2),datin
         IF(nz==82 .AND. na==208) THEN
           error(4,3) =  dummy(1)-dummy(2)
         ELSE IF(nz==60 .AND. na==142) THEN
           backspace(2)
           READ(2,*) na,nz,dummy(1:6),error(3,3)
         END IF
       END DO
49     CONTINUE
     ELSE IF (aline(1:8)=='   rms r') THEN
       DO ib=1,99
         READ(2,*,ERR=59) na,nz,dummy(1:2),datin
         IF(nz==82 .AND. na==208) THEN
           error(4,4) =  dummy(1)-dummy(2)
         ELSE IF(nz==88 .AND. na==214) THEN
           backspace(2)
           READ(2,*) na,nz,dummy(1:6),error(3,4)
         END IF
       END DO
59     CONTINUE
     ELSE IF(aline(23:34)=='r.m.s. error') THEN
       icount=1+icount
       IF(icount<5) THEN
          READ(aline,'(41x,1e12.4,21x,1e12.4)') error(1,icount),error(2,icount)
!          READ(aold2,*) dummy,error(3,icount)
       END IF
     END IF
     IF(icount>10) EXIT
   END DO
19 CONTINUE
!   WRITE(3,'(1x/1x)') 
   WRITE(3,'(a)') '# nr.obs rms-error average-error  subchi error-208Pb observable'   
   WRITE(3,'(a,4f9.4,a)') '1 ',error(1:4,1),' energy E [MeV]'
   WRITE(3,'(a,4f9.4,a)') '2 ',error(1:4,2),' diffr.radius [fm]'
   WRITE(3,'(a,4f9.4,a)') '3 ',error(1:4,3),' surface thick. [fm]'
   WRITE(3,'(a,4f9.4,a)') '4 ',error(1:4,4),' rms radius [fm]'
  CLOSE(2)
  CLOSE(3)
!  OPEN(3,file='error-oneline')
  WRITE(6,'(16f9.4)') error
!  CLOSE(3)
  STOP  'regular end'
99 STOP 'premature error or EoF'
END PROGRAM extract_errors
