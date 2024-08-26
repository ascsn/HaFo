PROGRAM extract_obs

! To extract basic observables from 'for002'.


IMPLICIT NONE
INTEGER, PARAMETER :: numnuc=99,nmaxp=126,nmaxn=192

CHARACTER*80 :: inline
CHARACTER*8 :: force

INTEGER :: na,nz,nn
INTEGER :: iread,i
REAL :: expbind,fitbind,expr,fitr,gap3(2)
REAL :: ebind(2,nmaxp,nmaxn),rms(2,nmaxp,nmaxn)




ebind=0.0
rms=0.0

DO iread=1,99999
  READ(5,'(a)',ERR=99,END=19) inline
  IF(inline(22:28) == 'comment') THEN
    READ(inline,'(29x,a8)',ERR=98,END=98) force
  END IF
  IF(inline(4:17) == 'binding-energy') THEN
    DO i=1,numnuc
      READ(5,*,ERR=18,END=98) na,nz,expbind,fitbind
      ebind(1,nz,na-nz)=expbind
      ebind(2,nz,na-nz)=fitbind
    END DO
 18 CONTINUE  
  END IF
  IF(inline(4:13) == 'rms radius') THEN
    DO i=1,numnuc
      READ(5,*,ERR=19,END=98) na,nz,expr,fitr
      rms(1,nz,na-nz)=expr
      rms(2,nz,na-nz)=fitr
    END DO
  END IF
END DO



STOP 'file count exhausted'

19 CONTINUE

WRITE(6,'(a)') '# Summary of results for the parameterization:'//force
WRITE(6,'(a)') '#0:   A   Z     E_B(exp)  E_B(fit) difference  '
DO nz=2,nmaxp; DO nn=2,nmaxn
  IF(ebind(2,nz,nn).NE.0.0)  WRITE(6,'(2i4,2f10.3,4(f8.3))') &
       nn+nz,nz,ebind(1,nz,nn),ebind(2,nz,nn),ebind(1,nz,nn)-ebind(2,nz,nn)
END DO;END DO
WRITE(6,'(//a)') '#1:   A   Z   r_rms(exp)  r_rms(fit)  difference '
DO nz=2,nmaxp; DO nn=2,nmaxn
  IF(rms(2,nz,nn).NE.0.0) &
    WRITE(6,'(2i4,2f10.4,4(f8.4))') &
       nn+nz,nz,rms(1,nz,nn),rms(2,nz,nn),rms(1,nz,nn)-rms(2,nz,nn)
END DO;END DO
WRITE(6,'(//a)') '#2: 3gap-2neutrons energies  '
WRITE(6,'(a)') '#   A   Z   exp    fit     difference '
DO nz=4,nmaxp-2; DO nn=4,nmaxn-2
  IF(ebind(2,nz,nn).NE.0.0 .AND. ebind(2,nz,nn+2).NE.0.0 .AND. &
     ebind(2,nz,nn-2).NE.0.0) THEN
    gap3(:)=0.5*(ebind(:,nz,nn+2)-2.0*ebind(:,nz,nn)+ebind(:,nz,nn-2))
    WRITE(6,'(2i4,2f10.4,4(f8.4))') &
       nn+nz,nz,gap3(:),gap3(1)-gap3(2)
  ENDIF
END DO;END DO
WRITE(6,'(//a)') '#3: 3gap-2protons energies  '
WRITE(6,'(a)') '#   A   Z   exp    fit     difference '
DO nz=4,nmaxp-2; DO nn=4,nmaxn-2
  IF(ebind(2,nz,nn).NE.0.0 .AND. ebind(2,nz+2,nn).NE.0.0 .AND. &
     ebind(2,nz-2,nn).NE.0.0) THEN
    gap3(:)=0.5*(ebind(:,nz+2,nn)-2.0*ebind(:,nz,nn)+ebind(:,nz-2,nn))
    WRITE(6,'(2i4,2f10.4,4(f8.4))') &
       nn+nz,nz,gap3(:),gap3(1)-gap3(2)
  ENDIF
END DO;END DO
WRITE(6,'(//a)') '#4: 3gap-2neutrons rms radii  '
WRITE(6,'(a)') '#   A   Z   exp    fit     difference '
DO nz=4,nmaxp-2; DO nn=4,nmaxn-2
  IF(rms(2,nz,nn).NE.0.0 .AND. rms(2,nz,nn+2).NE.0.0 .AND. &
     rms(2,nz,nn-2).NE.0.0) THEN
    gap3(:)=0.5*(rms(:,nz,nn+2)-2.0*rms(:,nz,nn)+rms(:,nz,nn-2))
    WRITE(6,'(2i4,2f10.3,4(f8.4))') &
       nn+nz,nz,gap3(:),gap3(1)-gap3(2)
  ENDIF
END DO;END DO
WRITE(6,'(//a)') '#5: 3gap-2protons rmsradii  '
WRITE(6,'(a)') '#   A   Z   exp    fit     difference '
DO nz=4,nmaxp-2; DO nn=4,nmaxn-2
  IF(rms(2,nz,nn).NE.0.0 .AND. rms(2,nz+2,nn).NE.0.0 .AND. &
     rms(2,nz-2,nn).NE.0.0) THEN
    gap3(:)=0.5*(rms(:,nz+2,nn)-2.0*rms(:,nz,nn)+rms(:,nz-2,nn))
    WRITE(6,'(2i4,2f10.4,4(f8.4))') &
       nn+nz,nz,gap3(:),gap3(1)-gap3(2)
  ENDIF
END DO;END DO
WRITE(6,'(//a)') '#6: 3gap-1neutrons energies  '
WRITE(6,'(a)') '#   A   Z   exp    fit     difference '
DO nz=4,nmaxp-1; DO nn=4,nmaxn-1
  IF(ebind(2,nz,nn).NE.0.0 .AND. ebind(2,nz,nn+1).NE.0.0 .AND. &
     ebind(2,nz,nn-1).NE.0.0) THEN
    gap3(:)=0.5*(ebind(:,nz,nn+1)-2.0*ebind(:,nz,nn)+ebind(:,nz,nn-1))
    WRITE(6,'(2i4,2f10.4,4(f8.4))') &
       nn+nz,nz,gap3(:),gap3(1)-gap3(2)
  ENDIF
END DO;END DO
WRITE(6,'(//a)') '#7: 3gap-1protons energies  '
WRITE(6,'(a)') '#   A   Z   exp    fit     difference '
DO nz=4,nmaxp-1; DO nn=4,nmaxn-1
  IF(ebind(2,nz,nn).NE.0.0 .AND. ebind(2,nz+1,nn).NE.0.0 .AND. &
     ebind(2,nz-1,nn).NE.0.0) THEN
    gap3(:)=0.5*(ebind(:,nz+1,nn)-2.0*ebind(:,nz,nn)+ebind(:,nz-1,nn))
    WRITE(6,'(2i4,2f10.4,4(f8.4))') &
       nn+nz,nz,gap3(:),gap3(1)-gap3(2)
  ENDIF
END DO;END DO
WRITE(6,'(//a)') '#8: 3gap-1neutrons rms radii  '
WRITE(6,'(a)') '#   A   Z   exp    fit     difference '
DO nz=4,nmaxp-1; DO nn=4,nmaxn-1
  IF(rms(2,nz,nn).NE.0.0 .AND. rms(2,nz,nn+1).NE.0.0 .AND. &
     rms(2,nz,nn-1).NE.0.0) THEN
    gap3(:)=0.5*(rms(:,nz,nn+1)-2.0*rms(:,nz,nn)+rms(:,nz,nn-1))
    WRITE(6,'(2i4,2f10.4,4(f8.4))') &
       nn+nz,nz,gap3(:),gap3(1)-gap3(2)
  ENDIF
END DO;END DO
WRITE(6,'(//a)') '#9: 3gap-1protons rmsradii  '
WRITE(6,'(a)') '#   A   Z   exp    fit     difference '
DO nz=4,nmaxp-1; DO nn=4,nmaxn-1
  IF(rms(2,nz,nn).NE.0.0 .AND. rms(2,nz+1,nn).NE.0.0 .AND. &
     rms(2,nz-1,nn).NE.0.0) THEN
    gap3(:)=0.5*(rms(:,nz+1,nn)-2.0*rms(:,nz,nn)+rms(:,nz-1,nn))
    WRITE(6,'(2i4,2f10.4,4(f8.4))') &
       nn+nz,nz,gap3(:),gap3(1)-gap3(2)
  ENDIF
END DO;END DO

STOP 'regular end'

98 STOP 'error or EoF in reading nucleus'
99 STOP 'error in reading text'

END PROGRAM extract_obs


