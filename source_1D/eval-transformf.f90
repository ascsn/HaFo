PROGRAM eval_transformf
!
! Reads transition formfactors for the RPA states and
! builds smooth frequency map of transition formfactors.
!
IMPLICIT NONE
!
REAL(8), PARAMETER :: sqpi=1.772454D0
REAL(8), PARAMETER :: widlow=0.2D0,widmin=8D0,widslop=6D0    !  208Pb
!REAL(8), PARAMETER ::  widlow=1D-1,widmin=6D0,widslop=6D0    !  208Pb
!REAL(8), PARAMETER ::  widlow=1D-1,widmin=8D0,widslop=5D0    !  208Pb
!REAL(8), PARAMETER ::  widlow=1D-1,widmin=80D0,widslop=5D0    !  test
!REAL(8), PARAMETER ::  widlow=1D-1,widmin=12D0,widslop=8D0   !  16O
!REAL(8), PARAMETER ::  widlow=1D-1,widmin=9D0,widslop=6D0    !  140Ce

REAL(8),PARAMETER :: enmin=0D0,enmax=37D0,enstep=widlow/5D0
!REAL(8),PARAMETER :: enmin=5D0,enmax=21D0,enstep=widlow/5D0
INTEGER :: nenerg,nenergmax,ityp

INTEGER,PARAMETER :: kgrid=199,kstates=399
LOGICAL,PARAMETER :: ttest=.true.

REAL(8) :: rad(kgrid),formf(kgrid,4,kstates),omega(kstates),bel(4)
INTEGER :: ngrid,nstates,nst,nstin,n
REAL(8) :: r,energ,dynwid,facwid,prhoact(4),prho0(4),transelement(4),fac
REAL(8) :: trans3  ! r**3 transition density for isovector (ityp=4)
REAL(8) :: xlimit

CHARACTER(LEN=12) :: filename

rad= 0D0
formf = 0D0
omega= 0D0
bel = 0D0

! reading data

WRITE(6,*) ' enter file name from which data are read (without extension "res"):'
!READ(5,'(a)') filename
filename = 'transform'
!OPEN(24,file='transdens.res')
OPEN(24,file=filename(1:LEN_TRIM(filename))//'.res')
OPEN(25,file=filename(1:LEN_TRIM(filename))//'.map')
READ(24,'(1x)') 
READ(24,'(8x,i5)') ngrid
IF(ttest) WRITE(*,*) ' ngrid=',ngrid
IF(ngrid > kgrid) STOP ' too large grid'
READ(24,'(1x)') 

DO nst=1,kstates

  READ(24,'(10x,i5,9x,e13.5,5x,4e13.5)',ERR=99,END=98) &
    nstin,omega(nst),bel
  IF(ttest) WRITE(*,*) ' nst,nstin,omega=',nst,nstin,omega(nst)
  IF(nst /= nstin) STOP ' states do not match'

  DO n=1,ngrid
!    IF(ttest) WRITE(*,*) 'n=',n
    READ(24,*,ERR=99,END=98) rad(n),formf(n,:,nst)
  END DO
!  READ(24,'(1x)',ERR=99,END=98) 
!  READ(24,'(1x)',ERR=99,END=98) 

  nstates = nst


  WRITE(7,'(i5,f10.3,4(1pg13.5))')  nst,omega(nst),formf(1,:,nst)

END DO
STOP ' too many states'

98 CONTINUE
CLOSE(24)
WRITE(7,'(2x/2x)')

!  folding to smooth distributions and printing

IF(ttest) WRITE(*,*) ' reading completed - now print folded densties'

nenergmax = NINT((enmax-enmin)/enstep)
DO nenerg=0,nenergmax
  energ = nenerg*enstep+enmin
  transelement = 0D0
  trans3 = 0D0
  DO n=1,ngrid
    prhoact = 0D0
    DO nst=1,nstates
      dynwid = max(widlow,(omega(nst)-widmin)/widslop)
      prhoact = formf(n,:,nst)**2 &
               *exp(-((omega(nst)-energ)/dynwid)**2)/dynwid + prhoact
    END DO
!    IF(n==1) prho0=prhoact/rad(n)**2
    WRITE(25,'(f12.5,f8.2,4(1pg13.5))') energ,rad(n), &
     (xlimit(prhoact(ityp)),ityp=1,4)
!     (xlimit(prhoact(ityp)/prho0(ityp)/rad(n)**2),ityp=1,4)
  END DO
  WRITE(25,'(1x/a)') '#  '
END DO

STOP 'regular end'
99 STOP ' error in input'
END PROGRAM eval_transformf


REAL(8) FUNCTION xlimit(x)
REAL(8) :: x
IF(ABS(x) > 1D-80) THEN
   xlimit = x
ELSE
   xlimit = 0D0
END IF
RETURN
END
