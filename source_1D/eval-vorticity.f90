PROGRAM eval_vorticita
!
! Reads vorticities for the RPA states and
! builds smooth frequency map of them.
!
IMPLICIT NONE
!
REAL(8), PARAMETER :: sqpi=1.772454D0
REAL(8), PARAMETER :: widlow=1D0,widmin=60D0,widslop=6D0    !  208Pb
!REAL(8), PARAMETER ::  widlow=1D-1,widmin=6D0,widslop=6D0    !  208Pb
!REAL(8), PARAMETER ::  widlow=1D-1,widmin=8D0,widslop=5D0    !  208Pb
!REAL(8), PARAMETER ::  widlow=1D-1,widmin=80D0,widslop=5D0    !  test
!REAL(8), PARAMETER ::  widlow=1D-1,widmin=12D0,widslop=8D0   !  16O
!REAL(8), PARAMETER ::  widlow=1D-1,widmin=9D0,widslop=6D0    !  140Ce

REAL(8),PARAMETER :: enmin=0D0,enmax=37D0,enstep=widlow/5D0
INTEGER :: nenerg,nenergmax,ityp

INTEGER,PARAMETER :: kgrid=99,kstates=199
LOGICAL,PARAMETER :: ttest=.false.

REAL(8) :: rad(kgrid),vort(kgrid,4,kstates),omega(kstates),bel(4)
INTEGER :: ngrid,nstates,nst,nstin,n
REAL(8) :: r,energ,dynwid,facwid,vortact(4),transelement(4),fac
REAL(8) :: xlimit

rad= 0D0
vort = 0D0
omega= 0D0
bel = 0D0

! reading data

OPEN(24,file='vorticity.res')
READ(24,'(1x)') 
READ(24,'(8x,i5)') ngrid
IF(ttest) WRITE(*,*) ' ngrid=',ngrid
IF(ngrid > kgrid) STOP ' too large grid'
READ(24,'(1x)') 

DO nst=1,kstates

  IF(ttest) WRITE(*,*) ' now reading state nr.',nst
  READ(24,'(10x,i5,9x,e13.5)',ERR=99,END=98) nstin,omega(nst)
  IF(ttest) WRITE(*,*) ' nst,nstin,omega=',nst,nstin,omega(nst)
  IF(nst /= nstin) STOP ' states do not match'

  DO n=1,ngrid
    READ(24,*,ERR=99,END=98) rad(n),vort(n,:,nst)
  END DO
!  READ(24,'(1x)',ERR=99,END=98) 
!  READ(24,'(1x)',ERR=99,END=98) 

  nstates = nst

! correct
!  DO ityp=1,4
!    fac = SQRT(bel(ityp))/SUM(vort(:,ityp,nst)*rad(:)**3)
    fac = SUM(vort(:,4,nst)*rad(:)**3)
    fac = fac/sqrt(fac**2)
    vort(:,:,nst) = fac*vort(:,:,nst)
!  END DO

  WRITE(7,'(i5,f10.3,1pg13.5)')  nst,omega(nst),SUM(vort(:,4,nst)*rad(:)**3)

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
  DO n=1,ngrid
    vortact = 0D0
    DO nst=1,nstates
      dynwid = max(widlow,(omega(nst)-widmin)/widslop)
      vortact = vort(n,:,nst)*exp(-((omega(nst)-energ)/dynwid)**2) + vortact
    END DO
    WRITE(6,'(f12.5,f8.2,4(1pg13.5))') energ,rad(n), &
     (xlimit(vortact(ityp)),ityp=1,4)
    transelement = transelement + rad(n)**5*vortact
  END DO
  WRITE(6,'(1x/a)') '#  '
  WRITE(7,'(f12.5,5(1pg13.5))') energ, &
    (xlimit((transelement(ityp)*rad(2)*4*3.1415926D0)**2),ityp=1,4)
END DO

STOP 'regular end'
99 STOP ' error in input'
END

REAL(8) FUNCTION xlimit(x)
REAL(8) :: x
IF(ABS(x) > 1D-80) THEN
   xlimit = x
ELSE
   xlimit = 0D0
END IF
RETURN
END
