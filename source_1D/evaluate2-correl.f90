PROGRAM evaluate_correl

IMPLICIT REAL(8) (a-h,o-z)

!     reads variance information from file and
!     computes the variance matrix of two selected
!     observables.


!INCLUDE 'comfit.inc'

CHARACTER (LEN=3) :: ident

!     'workdv' stores the gradients for observables
!     as indicated in the heading line for each set.

LOGICAL,PARAMETER :: tvarellips=.false.,talignmat=.true.,tplane=.true.

REAL(8) :: parmnatfac
INTEGER,PARAMETER :: kobs=59
INTEGER,PARAMETER :: kparmx=20
REAL(8) :: ZARRAY(kparmx,kparmx),DPARMO(kparmx)
REAL(8) :: par(kparmx)
REAL(8) :: workdv(0:kparmx,50),alig(kobs*kobs)
!REAL(8) :: workdv(0:kparmx,kobs),alig(kobs*kobs)
REAL(8) :: almat(kobs,kobs)
INTEGER :: ipo(2,kobs*kobs)
INTEGER :: irevpo(kobs)
INTEGER :: j,nobsmax

!     group types, with pygmy
INTEGER :: iplpoi(kobs)

!     group types, without pygmy
!      data iplpoi/03,10,04,12,07,11,05,06,14,15,
!     &            27,09,02,08,01,13,16,17,18,19,
!     &            20,21,22,23,24,25,26,28,29,30/

!     old sorting
!      data iplpoi/01,02,03,10,04,12,07,11,13,05,
!     &            06,08,09,30,14,20,22,21,15,16,
!     &            17,18,24,19,26,27,28,29,23,25/

CHARACTER (LEN=20) :: obs(kobs)
CHARACTER (LEN=2) :: aobs(59) = &
  (/'01','02','03','04','05','06','07','08','09',  &
    '10','11','12','13','14','15','16','17','18','19',  &
    '20','21','22','23','24','25','26','27','28','29',  &
    '30','31','32','33','34','35','36','37','38','39',  &
    '40','41','42','43','44','45','46','47','48','49',  &
    '50','51','52','53','54','55','56','57','58','59'/)
!DATA obs/'  E/A   ','  rho   ','   K    ','  m*/m  ',  &
!    '  asym  ','  dasym ','  kappa ',' E/Aneut',  &
!    'slopneut','  GMR   ','  GDR   ','  GQR   ',  &
!    '  GDR16 ',' polariz','  skin  ','skinsn74',  &
!    'skinsn82','skinsn90','skinca48','pol132sn',  &
!    ' pol48ca',' pol68ni',' pygmyPb','skinni68',  &
!    ' pygmySn',' dskin  ',' rmsn   ',' rdsmn  ',  &
!    ' formfn ','E_B148Sn','mpo208Pb','mpo132Sn', 'mpol68Ni', &
!    'mpol48Ca','qpo208Pb'/

LOGICAL,PARAMETER :: ttest=.true.

!---------------------------------------------------------------------

! define pointer
DO i=1,kobs
  iplpoi(i) = i
END DO
!iplpoi(1:35)= &
!  (/03,10,04,12,07,11,05,06,14,15, &
!    27,09,23,02,08,01,13,16,17,18,  &
!    19,20,21,22,24,25,26,28,29,30, 31,32,33,34,35/)


!     compute reverse pointer

IF(ttest) THEN
  WRITE(6,'(a)') 'irevpoi before:'
  WRITE(6,'(10i5)') irevpo
  WRITE(6,'(a)') 'iplpoi before:'
  WRITE(6,'(10i5)') iplpoi
END IF
DO i=1,kobs
  DO n=1,kobs
    itest = iplpoi(n)
    IF(i == itest) THEN
      irevpo(i) = n
      IF(ttest) WRITE(*,*) ' i,n,itest=',i,n,itest,irevpo(i)
    END IF
  END DO
END DO
IF(ttest) THEN
  WRITE(6,'(a)') 'irevpoi after:'
  WRITE(6,'(10i5)') irevpo
END IF

!     check pointers

DO i=1,kobs
  j = iplpoi(irevpo(i))
  IF(j /= i) THEN
    WRITE(6,'(a,2i5)') ' wrong pointer: i,j=',i,j
    STOP
  END IF
END DO

!     read nr. of parameters

READ(5,'(a3)') ident
OPEN(31,FILE='erroranalysis.'//ident)
READ(31,'(27x,i5,f8.2)') nparmx,parmnatfac
IF(ttest) WRITE(6,'(a,i5)') ' nparmx=',nparmx
IF(nparmx>kparmx) STOP ' NPARMX too large. Enhance KPARMX'

! override parameter variation

DO i=1,nparmx
  READ(31,'(1x)') 
END DO

! define simple diagonal 'zarray'

zarray = 0D0
DO i=1,nparmx
  zarray(i,i) = 1D0/parmnatfac**2
END DO

! read variations of observables

DO iobs=1,kobs
  READ(31,'(a)',err=99,end=99) obs(iobs)
  READ(31,*) workdv(0,iobs)
  READ(31,*) workdv(1:nparmx,iobs)
  IF(ttest) WRITE(*,*) ' observable read: nr=',iobs
  nobsmax=iobs
END DO
99 CONTINUE

CLOSE (31)

!     select, compute and print

  OPEN(31,FILE='varellips.dat')
  
  WRITE(31,'(a)') ' rho vs. E/A:'
  CALL varellips(workdv(0,1),workdv(0,2),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' rho vs. E/A: align=',align
  CLOSE(31)
IF(tvarellips) THEN
  OPEN(31,FILE='varellips.dat')
  
  WRITE(31,'(a)') ' polariz vs. skin:'
  CALL varellips(workdv(0,14),workdv(0,15),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' polariz vs. skin: align=',align
  
  WRITE(31,'(a)') ' asym vs. skin:'
  CALL varellips(workdv(0,5),workdv(0,15),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' asym vs. skin:    align=',align
  
  WRITE(31,'(a)') ' dasym vs. skin:'
  CALL varellips(workdv(0,6),workdv(0,15),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' dasym vs. skin:   align=',align
  
  WRITE(31,'(a)') ' kappa vs. skin:'
  CALL varellips(workdv(0,7),workdv(0,15),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' kappa vs. skin:   align=',align
  
  WRITE(31,'(a)') ' K vs. skin:'
  CALL varellips(workdv(0,3),workdv(0,15),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' K vs. skin:       align=',align
  
  WRITE(31,'(a)') ' m*/m vs. skin:'
  CALL varellips(workdv(0,4),workdv(0,15),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' m*/m vs. skin:    align=',align
  
  WRITE(31,'(a)') ' slope_neut vs. skin:'
  CALL varellips(workdv(0,9),workdv(0,15),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' slope_nt vs. skin:align=',align
  
  WRITE(31,'(a)') ' GDR vs. skin:'
  CALL varellips(workdv(0,11),workdv(0,15),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' GDR vs. skin:    align=',align
  
  WRITE(31,'(a)') ' asym vs. daysm:'
  CALL varellips(workdv(0,5),workdv(0,6),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' asym vs. daysm:  align=',align
  
  WRITE(31,'(a)') ' polariz  vs. GDR:'
  CALL varellips(workdv(0,14),workdv(0,11),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' polariz vs. GDR: align=',align
  
  WRITE(31,'(a)') ' polariz  vs. asym:'
  CALL varellips(workdv(0,14),workdv(0,5),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' polariz vs. asym:align=',align
  
  WRITE(31,'(a)') ' polariz  vs. dasym:'
  CALL varellips(workdv(0,14),workdv(0,6),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' polariz vs.dasym:align=',align
  
  WRITE(31,'(a)') ' polariz  vs. slope_neut:'
  CALL varellips(workdv(0,14),workdv(0,9),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' polariz vs.slope:align=',align
  
  WRITE(31,'(a)') ' skinsn74 vs. skinsn82:'
  CALL varellips(workdv(0,16),workdv(0,17),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' skinsn74 vs. skinsn82: align=',align
  
  WRITE(31,'(a)') ' skinsn74 vs. skinsn90:'
  CALL varellips(workdv(0,16),workdv(0,18),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' skinsn74 vs. skinsn90: align=',align
  
  WRITE(31,'(a)') ' skinsn82 vs. skinsn90:'
  CALL varellips(workdv(0,17),workdv(0,18),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' skinsn82 vs. skinsn90: align=',align
  
  WRITE(31,'(a)') ' skinpb vs. skinsn82:'
  CALL varellips(workdv(0,15),workdv(0,17),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' skinpb vs. skinsn82: align=',align
  
  WRITE(31,'(a)') ' skinSn vs. pygmySn:'
  CALL varellips(workdv(0,17),workdv(0,25),nparmx,31,align)
  WRITE(31,'(5x)')
  WRITE(31,'(5x)')
  WRITE(6,'(a,f8.5)') ' skinSn vs. pygmySn: align=',align
  
  CLOSE(31)
END IF

!     systematics of correlations

WRITE(7,'(a)') ' alignment matrix:'
is = 0
DO i=1,nobsmax ! kobs
  DO j=1,i-1
    is = 1+is
    CALL varellips(workdv(0,i),workdv(0,j),nparmx,31,alig(is))
    ipo(1,is) = i
    ipo(2,is) = j
    IF(ttest) WRITE(6,'(a,2i5,1pg13.5)') ' i,j,alig=',i,j,alig(is)
  END DO
  WRITE(7,'(100(f7.4))') (alig(k),k=is-i+2,is)
END DO
IF(ttest) WRITE(*,*) ' ipo(1):',ipo(1,1:is)
IF(ttest) WRITE(*,*) ' ipo(2):',ipo(2,1:is)

!     output of correlation matrix for color scale plot

IF(talignmat) THEN
  IF(ttest) WRITE(*,*) ' now print alignment matrix'
  is = 0
  DO i=1,nobsmax   ! kobs
    DO j=1,i-1
      is = 1+is
      almat(i,j) = alig(is)
      almat(j,i) = alig(is)
    END DO
    almat(i,i) = 1D0
  END DO
  
  OPEN(11,FILE='alignmatrix.gnu')
  WRITE(11,'(a)') 'set term post eps enhanced color "bold" 16'
  WRITE(11,'(a)') 'set output "alignmatrix.eps"'
  WRITE(11,'(a)') 'set border 15 lw 3'
  WRITE(11,'(a)') 'set tmargin 0'
  WRITE(11,'(a)') 'set bmargin 0'
  WRITE(11,'(a)') 'set lmargin 0'
  WRITE(11,'(a)') 'set rmargin 0'
  WRITE(11,'(a)') 'set size 1.05,1.15'
  WRITE(11,'(a)') 'set multiplot'
  WRITE(11,'(a)') '#'
  WRITE(11,'(a)') 'set origin 0.07,-0.05'
  WRITE(11,'(a)') 'set size 1,1.4'
  WRITE(11,'(a)') 'set pm3d map'
  WRITE(11,'(a)') 'set xrange [0:36]'
  WRITE(11,'(a)') 'set yrange [0:36]'
  WRITE(11,'(a)') 'set cbrange [0:1]'
  WRITE(11,'(a,90(a,i2,a,i3,a),a)')  &
      'set ytics (',('"',i,': '//obs(iplpoi(i))//'"',i,',' ,i=1,kobs-1),  &
      '"',kobs,':'//obs(iplpoi(kobs))//'"',kobs,')'
  WRITE(11,'(a)') 'sp "alignmatrix.gnu" u 1:2:3 not'
!        write(11,'(a)') 'pause -1'
  WRITE(11,'(a)') 'q'
  CALL flush(11)
  DO i=1,nobsmax ! kobs
    xi = i-0.5D0
    ia = iplpoi(i)
    DO j=1,nobsmax ! kobs
      xj = j-0.5D0
      ja = iplpoi(j)
      WRITE(11,'(2f8.1,1pg13.5)') xi,xj,almat(ia,ja)
      WRITE(11,'(2f8.1,1pg13.5)') xi,xj+1D0,almat(ia,ja)
    END DO
    WRITE(11,'(5x)')
    DO j=1,nobsmax ! kobs
      xj = j-0.5D0
      ja = iplpoi(j)
      WRITE(11,'(2f8.1,1pg13.5)') xi+1D0,xj,almat(ia,ja)
      WRITE(11,'(2f8.1,1pg13.5)') xi+1D0,xj+1D0,almat(ia,ja)
    END DO
    WRITE(11,'(5x)')
    WRITE(11,'(5x)')
  END DO
  CLOSE(11)
  
END IF

!     output sorted for each observable separately

DO i=1,is
  aligmin = alig(i)
  DO j=i+1,is
    IF(alig(j) < aligmin) THEN
      a = alig(j)
      alig(j) = aligmin
      aligmin = a
      alig(i) = a
      isav = ipo(1,j)
      ipo(1,j) = ipo(1,i)
      ipo(1,i) = isav
      isav = ipo(2,j)
      ipo(2,j) = ipo(2,i)
      ipo(2,i) = isav
    END IF
  END DO
END DO

DO i=1,is
  WRITE(7,'(a,1x,a,f8.5)') obs(ipo(1,i)), obs(ipo(2,i)),alig(i)
END DO

DO j=1,nobsmax ! kobs
  ja = iplpoi(j)
  jb = irevpo(j)
  OPEN(11,FILE='align.'//aobs(j))
  WRITE(11,'(a,i3,a)') '#alignement for observable:',j,obs(ja)
  WRITE(11,'(a,2(1pg13.5))') '# average,variance=', &
     workdv(0,j),variance(workdv(0,j),nparmx)
  DO i=is,1,-1
    IF(irevpo(ipo(1,i)) == j) THEN
      WRITE(11,'(2i5,1x,a,f8.5)')  &
          irevpo(ipo(2,i)),ipo(2,i),obs(ipo(2,i)),alig(i)
    ELSE IF(irevpo(ipo(2,i)) == j) THEN
      WRITE(11,'(2i5,1x,a,f8.5)')  &
          irevpo(ipo(1,i)),ipo(1,i),obs(ipo(1,i)),alig(i)
    END IF
  END DO
  CLOSE(11)
END DO

IF(tplane) THEN
  WRITE(*,*) ' alignment in plane'
  DO j2=1,nobsmax ! kobs
    OPEN(13,FILE='cosplane.'//aobs(j2))
    WRITE(6,'(a,i5)') ' cosines with '//obs(j2)//'  j=',j2
    WRITE(13,'(a/5a)') ' cosines with '//obs(j2)//' in order of:',  &
        obs(3),obs(5),obs(6),obs(4),obs(7)
    CALL cosplane(workdv(0,3),workdv(0,5),workdv(0,6),  &
                  workdv(0,4),workdv(0,7),workdv(0,43),  &
                  workdv(0,j2),nparmx)
    CLOSE(13)
  END DO
END IF

STOP

CONTAINS
!-----variance -----------------------------------------------------

REAL(8) FUNCTION variance(grad1,nparmx)
IMPLICIT REAL(8) (a-h,o-z)


REAL(8), INTENT(IN)             :: grad1(0:nparmx)
INTEGER, INTENT(IN)                      :: nparmx

!     computes variance for observable on 'grad1'

INTEGER :: i,j
REAL(8) :: err11

!---------------------------------------------------------------
err11   = 0.0

IF(nparmx > 0) THEN
  DO i=1,nparmx
    DO j=1,nparmx
      err11  = grad1(i)*grad1(j)*zarray(i,j) + err11
    END DO
  END DO
END IF

variance = SQRT(err11)

RETURN
END FUNCTION variance

!-----varellips -----------------------------------------------------

SUBROUTINE varellips(grad1,grad2,nparmx,iu,alignact)
IMPLICIT REAL(8) (a-h,o-z)


REAL(8), INTENT(IN)             :: grad1(0:nparmx)
REAL(8), INTENT(IN)             :: grad2(0:nparmx)
INTEGER, INTENT(IN)                      :: nparmx
INTEGER, INTENT(IN)                  :: iu
REAL(8), INTENT(OUT)            :: alignact

!INCLUDE 'comfit.inc'

!     computes 2x2 variance matrix for observables on 'grad1'
!     and 'grad2'.



INTEGER, PARAMETER :: kplot=50
INTEGER :: i,j
REAL(8) :: x,err11,err22,err12,eig1,eig2,angle

!---------------------------------------------------------------


!     compute variances

err11   = 0.0
err22   = 0.0
err12   = 0.0

IF(nparmx > 0) THEN
  DO i=1,nparmx
    DO j=1,nparmx
      err11  = grad1(i)*grad1(j)*zarray(i,j) + err11
      err22  = grad2(i)*grad2(j)*zarray(i,j) + err22
      err12  = grad1(i)*grad2(j)*zarray(i,j) + err12
    END DO
  END DO
END IF
CALL  diag2(err11,err12,err22,eig1,eig2,angle)

!     print eigenvalues

WRITE(6,'(a,3(1pg13.5)/a,3(1pg13.5))')  &
    ' variance  matrix:',err11,err22,err12,  &
    ' eigenvalues etc :',eig1,eig2,angle/3.1415926*180D0

!     compute error ellipsoid and print

alignact = err12**2/(err11*err22)
det = err11*err22-err12**2
erri11 = err22/det
erri22 = err11/det
erri12 = -err12/det
xlim = SQRT(ABS(erri22))/SQRT(erri11*erri22-erri12**2)
WRITE(6,'(a,4(1pg12.4))') 'sub-matrix:',erri11,erri22,erri12,xlim
DO j=-kplot,kplot
  x = j*xlim/kplot
  diff = SQRT(MAX(1D0/erri22 +x**2*(erri12**2-erri11*erri22)/erri22**2,  &
      1D-20))
  rsum  = -erri12*x/erri22
  WRITE(iu,'(3(1pg15.7))') x+grad1(0),rsum+diff+grad2(0)
END DO
DO j=-kplot+1,kplot
  x = -j*xlim/kplot
  diff = SQRT(MAX(1D0/erri22 +x**2*(erri12**2-erri11*erri22)/erri22**2,  &
      1D-20))
  rsum  = -erri12*x/erri22
  WRITE(iu,'(3(1pg15.7))') x+grad1(0),rsum-diff+grad2(0)
END DO
WRITE(iu,'(a,1pg13.5)') '# alignement: ',alignact

RETURN
END SUBROUTINE varellips

!-----diag2 -----------------------------------------------------

SUBROUTINE diag2(a11,a12,a22,eig1,eig2,angle)
IMPLICIT REAL(8) (a-h,o-z)



REAL(8), INTENT(IN OUT)         :: a11
REAL(8), INTENT(IN)             :: a12
REAL(8), INTENT(IN OUT)         :: a22
REAL(8), INTENT(OUT)            :: eig1
REAL(8), INTENT(OUT)            :: eig2
REAL(8), INTENT(OUT)            :: angle

!     diagonalizes 2x2 matrix,
!     returns eigenvectors and angle of rotation

!---------------------------------------------------------------

deig = SQRT((0.5D0*(a11-a22))**2+a12**2)
eig1 = 0.5D0*(a11+a22)+deig
eig2 = 0.5D0*(a11+a22)-deig
angle = ATAN((a11-eig1)/a12)

RETURN
END SUBROUTINE diag2
!test c
!test       implicit double precision (a-h,o-z)
!test c
!test       angin = 3.1415926D0/6.0
!test       c = cos(angin)
!test       s = sin(-angin)
!test       call diag2(c*c,s*c,s*s,eig1,eig2,angle)
!test       write(6,'(5(1pg13.5))') angin,eig1,eig2,angle,
!test      &  angle+3.1415926D0/4D0
!test
!test
!test       stop
!test       end
!-----cosplane--------------------------------------------------------

SUBROUTINE cosplane(g1,g2,g3,g4,g5,g6,gobs,nparmx)
IMPLICIT REAL(8) (a-h,o-z)


REAL(8), INTENT(IN)             :: g1(0:kparmx)
REAL(8), INTENT(IN)             :: g2(0:kparmx)
REAL(8), INTENT(IN)             :: g3(0:kparmx)
REAL(8), INTENT(IN)             :: g4(0:kparmx)
REAL(8), INTENT(IN OUT)         :: g5(0:kparmx)
REAL(8), INTENT(IN OUT)         :: g6(0:kparmx)
REAL(8), INTENT(IN)             :: gobs(0:kparmx)
INTEGER, INTENT(IN OUT)                  :: nparmx

!     Computes cosine of observable 'gobs' with the plane
!     of observables 'g1'--'g6'.
!     The variance matrix 'zarray' comes through common (comfit.inc).


!INCLUDE 'comfit.inc'
INTEGER, PARAMETER :: ksub=6
REAL(8) :: g(0:kparmx,ksub)
REAL(8) :: zsub(ksub,ksub),zinv(ksub,ksub),gobs_sub(ksub)
INTEGER :: i,j,nsub
REAL(8) :: determ,xngsub

!---------------------------------------------------------------------

g(:,1) = g1
g(:,2) = g2
g(:,3) = g3
g(:,4) = g4
g(:,5) = g5
g(:,6) = g6

!     prepare matrix and vector in subspace

DO i=1,ksub
  DO j=1,ksub
    zsub(i,j) = zoverlap(g(0,i),g(0,j),nparmx)
  END DO
END DO
DO i=1,ksub
  gobs_sub(i) = zoverlap(g(0,i),gobs,nparmx)
END DO

!     evaluate cosines for increasing number vectors

xngobs = zoverlap(gobs,gobs,nparmx)
DO nsub=1,ksub
  zinv  = zsub
  CALL matin2(zinv,nsub,ksub,nerror,determ)
  xngsub = 0D0
  DO i=1,nsub
    DO j=1,nsub
      xngsub = xngsub + gobs_sub(i)*zinv(i,j)*gobs_sub(j)
    END DO
  END DO
  WRITE(13,'(a,i4,a,3(1pg13.5))')  &
      ' nr. of obs =',nsub,': norm,subnorm,cos=', xngobs,xngsub,(xngsub/xngobs)
END DO

RETURN
END SUBROUTINE cosplane

!-----zoverlap--------------------------------------------------------

FUNCTION zoverlap(g1,g2,nparmx)
IMPLICIT REAL(8) (a-h,o-z)

!     computes the overlap of 'g1' and 'g2' using norm matrix 'zarray'.


REAL(8), INTENT(IN)             :: g1(0:kparmx)
REAL(8), INTENT(IN)             :: g2(0:kparmx)
INTEGER, INTENT(IN)                      :: nparmx

!INCLUDE 'comfit.inc'


!---------------------------------------------------------------------

xno = 0D0
DO i=1,nparmx
  DO j=1,nparmx
    xno = xno + g1(i)*zarray(i,j)*g2(j)
  END DO
END DO
zoverlap = xno
RETURN
END FUNCTION zoverlap

!-----matin2----------------------------------------------------------

SUBROUTINE matin2 (a,nmax,ndim,nerror,determ)
IMPLICIT REAL(8) (a-h,o-z)


REAL(8), INTENT(IN OUT)         :: a(ndim,ndim)
INTEGER, INTENT(IN)                      :: nmax
INTEGER, INTENT(IN)                  :: ndim
INTEGER, INTENT(IN OUT)                  :: nerror
REAL(8), INTENT(IN OUT)         :: determ
!      include 'params.inc'
INTEGER, PARAMETER :: kindex=500

LOGICAL, PARAMETER :: tunclean=.true.

!     inverts the matrix 'a'.

!     the parameters on the i/o list are:

!       a       = at input:  matrix to be inverted
!                 at output: inverse matrix
!       nmax    = number of actual rows and columns
!       ndim    = dimension of matrix in the calling routine
!       nerror  = error flag
!       determ  = determinant of input matrix 'a'


DIMENSION  work(ndim,ndim)
DIMENSION asav(ndim,ndim),alud(ndim,ndim),rhs(ndim)
DIMENSION INDEX(ndim)

!-----------------------------------------------------------------------------

DO i=1,nmax
  DO j=1,nmax
    work(i,j) = 0D0
    asav(j,i) = a(j,i)
  END DO
  work(i,i) = 1D0
END DO

CALL ludcmp(a,nmax,ndim,INDEX,det)

DO j=1,nmax
  CALL lubksb(a,nmax,ndim,INDEX,work(1,j))
END DO

DO i=1,nmax
  DO j=1,nmax
    alud(i,j) = a(i,j)
    a(i,j) =  work(i,j)
  END DO
END DO

IF(tunclean) RETURN


!     clean solution iteratively

DO i=1,nmax
  DO j=1,nmax
    rhs(j) = 0D0
    work(j,i) = a(j,i)
  END DO
  rhs(i) = 1D0
  CALL mprove(asav,alud,nmax,ndim,INDEX,rhs,work(1,i))
END DO


DO i=1,nmax
  DO j=1,nmax
    a(i,j) = work(i,j)
  END DO
END DO


RETURN
END SUBROUTINE matin2


!-----from numerical recipes

SUBROUTINE mprove(a,alud,n,np,indx,b,x)
IMPLICIT REAL(8) (a-h,o-z)

REAL(8), INTENT(IN)             :: a(np,np)
REAL(8), INTENT(IN OUT)         :: alud(np,np)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                  :: np
INTEGER, INTENT(IN OUT)                  :: indx(n)
REAL(8), INTENT(IN)             :: b(n)
REAL(8), INTENT(IN OUT)         :: x(n)
INTEGER, PARAMETER :: nmax=100
DIMENSION  r(nmax)
REAL*8 sdp

DO i=1,n
  sdp =-b(i)
  DO j=1,n
    sdp = sdp + a(i,j)*x(j)
  END DO
  r(i) = sdp
END DO
CALL lubksb(alud,n,np,indx,r)
DO i=1,n
  x(i) = x(i) - r(i)
END DO

RETURN
END SUBROUTINE mprove


SUBROUTINE ludcmp(a,n,np,indx,d)
IMPLICIT REAL(8) (a-h,o-z)

REAL(8), INTENT(IN OUT)         :: a(np,np)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                  :: np
INTEGER, INTENT(OUT)                     :: indx(n)
REAL(8), INTENT(OUT)            :: d


INTEGER, PARAMETER :: nmax=500
REAL(8), PARAMETER :: tiny=1.0D-20
INTEGER :: i,imax,j,k
REAL(8) :: aamax,dum,rsum,vv(nmax)

d=1.d0
DO  i=1,n
  aamax=0.d0
  DO  j=1,n
    IF (ABS(a(i,j)) > aamax) aamax=ABS(a(i,j))
  END DO
  IF (aamax == 0.d0) PAUSE 'singular matrix in ludcmp'
  vv(i)=1.d0/aamax
END DO
DO  j=1,n
  DO  i=1,j-1
    rsum=a(i,j)
    DO  k=1,i-1
      rsum=rsum-a(i,k)*a(k,j)
    END DO
    a(i,j)=rsum
  END DO
  aamax=0.d0
  DO  i=j,n
    rsum=a(i,j)
    DO  k=1,j-1
      rsum=rsum-a(i,k)*a(k,j)
    END DO
    a(i,j)=rsum
    dum=vv(i)*ABS(rsum)
    IF (dum >= aamax) THEN
      imax=i
      aamax=dum
    END IF
  END DO
  IF (j /= imax)THEN
    DO  k=1,n
      dum=a(imax,k)
      a(imax,k)=a(j,k)
      a(j,k)=dum
    END DO
    d=-d
    vv(imax)=vv(j)
  END IF
  indx(j)=imax
  IF(a(j,j) == 0.d0)a(j,j)=tiny
  IF(j /= n)THEN
    dum=1.d0/a(j,j)
    DO  i=j+1,n
      a(i,j)=a(i,j)*dum
    END DO
  END IF
END DO
RETURN
END SUBROUTINE ludcmp

SUBROUTINE lubksb(a,n,np,indx,b)
IMPLICIT REAL(8) (a-h,o-z)

REAL(8), INTENT(IN)             :: a(np,np)
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN)                  :: np
INTEGER, INTENT(IN)                      :: indx(n)
REAL(8), INTENT(IN OUT)         :: b(n)


INTEGER :: i,ii,j,ll
REAL(8) :: rsum

ii=0
DO  i=1,n
  ll=indx(i)
  rsum=b(ll)
  b(ll)=b(i)
  IF (ii /= 0)THEN
    DO  j=ii,i-1
      rsum=rsum-a(i,j)*b(j)
    END DO
  ELSE IF (rsum /= 0.d0) THEN
    ii=i
  END IF
  b(i)=rsum
END DO
DO  i=n,1,-1
  rsum=b(i)
  DO  j=i+1,n
    rsum=rsum-a(i,j)*b(j)
  END DO
  b(i)=rsum/a(i,i)
END DO
RETURN
END SUBROUTINE lubksb


END PROGRAM evaluate_correl
