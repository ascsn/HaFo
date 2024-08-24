PROGRAM extract_obs

! To extract basic observables from 'for002'.


IMPLICIT NONE
INTEGER, PARAMETER :: numnuc=199,nmaxp=126,nmaxn=192

CHARACTER*80 :: inline
CHARACTER*8 :: force

INTEGER :: na,nz,nn
INTEGER :: iread,i
REAL :: exp,fit
LOGICAL :: ten=.FALSE.



exp=0.0
fit=0.0

DO iread=1,99999
  READ(5,'(a)',ERR=99,END=19) inline
!  IF(ten) WRITE(6,'(a)') inline
  IF(inline(22:28) == 'comment') THEN
    READ(inline,'(29x,a8)',ERR=98,END=98) force
  END IF
  IF(inline(4:17) == 'binding-energy') THEN
    ten=.TRUE.
    WRITE(6,'(a)') 'binding-energy'
    DO i=1,numnuc
      fit=0.0
      READ(5,*,ERR=18,END=98) na,nz,exp,fit
      IF(fit.NE.0.0)  WRITE(6,'(2i4,2f10.3,4(f8.3))') &
       na,nz,exp,fit,fit-exp
    END DO
 18 CONTINUE  
  END IF
  IF(inline(4:16) == 'diffr.-radius') THEN
    WRITE(6,'(a)') 'diffr.-radius'
    DO i=1,numnuc
      fit=0.0
      READ(5,*,ERR=17,END=98) na,nz,exp,fit
      IF(fit.NE.0.0)  WRITE(6,'(2i4,2f10.3,4(f8.3))') &
       na,nz,exp,fit,fit-exp
    END DO
 17 CONTINUE  
  END IF
  IF(inline(4:20) == 'surface thickness') THEN
    WRITE(6,'(a)') 'surface thickness'
    DO i=1,numnuc
      fit=0.0
      READ(5,*,ERR=16,END=98) na,nz,exp,fit
      IF(fit.NE.0.0)  WRITE(6,'(2i4,2f10.3,4(f8.3))') &
       na,nz,exp,fit,fit-exp
    END DO
 16 CONTINUE  
  END IF
  IF(inline(4:13) == 'rms radius') THEN
    WRITE(6,'(a)') 'rms radius'
    DO i=1,numnuc
      fit=0.0
      READ(5,*,ERR=19,END=98) na,nz,exp,fit
      IF(fit.NE.0.0)  WRITE(6,'(2i4,2f10.3,4(f8.3))') &
       na,nz,exp,fit,fit-exp
    END DO
  END IF
END DO



STOP 'file count exhausted'

19 CONTINUE

STOP 'regular end'

98 STOP 'error or EoF in reading nucleus'
99 STOP 'error in reading text'

END PROGRAM extract_obs


