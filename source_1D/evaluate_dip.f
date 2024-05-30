c-----evaluate_dip-----------------------------------------------------
c
      SUBROUTINE evaluate_dip(rhoin,fin,rhomin,rhomax,ngrid,rstep,delq)
c
c     Computes minimal and maximal density for given distribution.
c     List parameters:
c       rhoin    = array of local density 
c       fin      = array of formfactor of rhoin
c       rhomin   = minial density
c       rhomax   = maximal density
c       ngrid    = number of grid points
c       rstep    = local grid spacing
c       qstep    = spacing of momentum grid
c
      IMPLICIT REAL(8)  (a-h,o-z)      
c
      REAL(8), INTENT(IN) :: rhoin(ngrid),fin(ngrid),rstep,delq
      REAL(8), INTENT(OUT) :: rhomin,rhomax
      INTEGER, INTENT(IN) :: ngrid
c
      REAL(8), PARAMETER :: precis=1D-15
      INTEGER, PARAMETER :: itqmax=30
      LOGICAL, PARAMETER :: ttest=.TRUE.
c
c---------------------------------------------------------------------
c
      WRITE(6,*) 'FIN:'
      WRITE(6,'(5(1pg13.5))') fin
c
c     rough localization of minima and maxima
c
      rhomax0=0D0
      irlow = NINT(fin(1))/20+2
      DO ir=9,ngrid
        IF(rhoin(ir)>rhomax0) THEN
          rhomax0 = rhoin(ir)
          irmax=ir
          WRITE(6,'(a,i5,1pg13.5)') 'IR,rhomax0=',IR,rhomax0
        END IF
      END DO
      rhomin0=1D3
      DO ir=1,irmax
        IF(rhoin(ir)<rhomin0) THEN
          rhomin0 = rhoin(ir)
          irmin=ir
        END IF
      END DO


c
c     find exact minimum by searching for zero of derivative
c     with bracketing
c
      IF(irmin==1) THEN
        rhomin = rhoin(1)
      ELSE
        rlow=(irmin-2)*rstep
        rup=irmin*rstep
        CALL FourBesDerv(fin,rlow,rholow,ngrid,rstep,delq)
        CALL FourBesDerv(fin,rup,rhoup,ngrid,rstep,delq)
        IF(ttest) WRITE(6,'(a,2i4,6(1pg13.5))') 'CENTRDEP:',
     &        irmin,irmax,rlow,rup,rholow,rhoup
        IF(rholow*rhoup>0D0) THEN
          WRITE(*,*) 'EVALUATE_DIP: wrong start minimum'
          rhomin =0D0
          rhomax=rhomax0
          RETURN
        END IF
        DO it=1,itqmax
          rmid=(rup+rlow)/2D0
          CALL FourBesDerv(fin,rmid,rhomid,ngrid,rstep,delq)
          IF(ttest) WRITE(6,'(a,i3,6(1pg13.5))') 'CENTRDEP:',it,
     &        rlow,rup,rmid,rholow,rhoup,rhomid
          IF(ABS(rhomid)<precis) EXIT
          IF(rhomid*rholow<0D0) THEN
            rup=rmid
            rhoup=rhomid
          ELSE
            rlow=rmid
            rholow=rhomid
          END IF
        END DO
        CALL FourBesBack(fin,rmid,rhomin,ngrid,rstep,delq)
        IF(ttest) WRITE(6,'(2(1pg13.5))') rmid,rhomin
      END IF

c
c     find maximum
c
      rlow=(irmax-2)*rstep
      rup=(irmax)*rstep
      CALL FourBesDerv(fin,rlow,rholow,ngrid,rstep,delq)
      CALL FourBesDerv(fin,rup,rhoup,ngrid,rstep,delq)
      IF(ttest) THEN
        CALL FourBesBack(fin,rlow,rhol,ngrid,rstep,delq)
        CALL FourBesBack(fin,rup,rhou,ngrid,rstep,delq)
        WRITE(6,'(a,2i4,6(1pg13.5))') 'CENTRDEP:',
     &        irmin,irmax,rlow,rup,rholow,rhoup,rhol,rhou
      END IF
      IF(rholow*rhoup>0D0) THEN
        WRITE(*,*) 'EVALUATE_DIP: wrong start maximum'
        rhomin=0D0
        rhomax=rhomax0
        RETURN
      END IF
      DO it=1,itqmax
        rmid=(rup+rlow)/2D0
        CALL FourBesDerv(fin,rmid,rhomid,ngrid,rstep,delq)
        IF(ABS(rhomid)<precis) EXIT
        IF(rhomid*rholow<0D0) THEN
          rup=rmid
          rhoup=rhomid
        ELSE
          rlow=rmid
          rholow=rhomid
        END IF
      END DO
      CALL FourBesBack(fin,rmid,rhomax,ngrid,rstep,delq)
c
      WRITE(6,'(a,2i5,4(1pg13.5))') 
     &   'EVALUATE_DIP: irmin,irmax,rhomin,rhomax=',
     &   irmin,irmax,rhoin(irmin),rhomin,rhoin(irmax),rhomax
c
      RETURN


      END SUBROUTINE evaluate_dip

c-----FourBesBack-----------------------------------------------------
c
      SUBROUTINE FourBesBack(fin,ract,rhoout,ngrid,rstep,delq)
c
c     Backward Fourier-Bessel transformation.
c     List parameters:
c       fin      = array of given formfactor
c       ract     = actual radial point
c       rhoout   = density at oint 'ract'
c       ngrid    = number of grid points
c       rstep    = local grid spacing
c       qstep    = spacing of momentum grid
c
      IMPLICIT REAL(8)  (a-h,o-z)      
c
      REAL(8), INTENT(IN) :: fin(ngrid),ract,rstep,delq
      REAL(8), INTENT(OUT) :: rhoout
      INTEGER, INTENT(IN) :: ngrid
c
      REAL(8), PARAMETER :: PI=3.141592654 
c
c     charge density by fourier back-transformation
c
      dpi3   = 0.5D0/(pi*pi)
      IF(ract==0D0) THEN
c                                                         r=0.0 separately
        acc    = 0D0
        DO iq=2,ngrid
          qa     = delq*iq-delq
          acc    = qa*qa*fin(iq)+acc
        ENDDO
        rhoout = ABS(acc*delq*dpi3)
      ELSE
c                                                         now all other r
        dc     = COS(ract*delq)
        ds     = SIN(ract*delq)
        c      = 1D0
        s      = 0D0
        acc    = 0D0
        qa     = 0D0
        DO iq=2,ngrid
          qa     = delq+qa
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
          acc    = qa*s*fin(iq)+acc
        ENDDO
!        WRITE(6,'(a,i4,4(1pg13.5))') 'FourBesBack:',
!     &          ngrid,acc,delq,dpi3,ract
        rhoout = (acc*delq*dpi3/ract)
      END IF

      RETURN

      END SUBROUTINE FourBesBack
c-----FourBesDerv-----------------------------------------------------
c
      SUBROUTINE FourBesDerv(fin,ract,rhoout,ngrid,rstep,delq)
c
c     First derivative via Backward Fourier-Bessel transformation.
c     List parameters:
c       fin      = array of given formfactor
c       ract     = actual radial point
c       rhoout   = density at oint 'ract'
c       ngrid    = number of grid points
c       rstep    = local grid spacing
c       qstep    = spacing of momentum grid
c
      IMPLICIT REAL(8)  (a-h,o-z)      
c
      REAL(8), INTENT(IN) :: fin(ngrid),ract,rstep,delq
      REAL(8), INTENT(OUT) :: rhoout
      INTEGER, INTENT(IN) :: ngrid
c
      REAL(8), PARAMETER :: PI=3.141592654 
c
c     charge density by fourier back-transformation
c
      dpi3   = 0.5D0/(pi*pi)
      IF(ract==0D0) THEN
c                                                         r=0.0 separately
        rhoout = 0D0
      ELSE
c                                                         now all other r
        dc     = COS(ract*delq)
        ds     = SIN(ract*delq)
        c      = 1D0
        s      = 0D0
        acc    = 0D0
        qa     = 0D0
        DO iq=2,ngrid
          qa     = delq+qa
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
          acc    = qa*(qa*c-s/ract)*fin(iq)+acc
        ENDDO
!        WRITE(6,'(a,i4,4(1pg13.5))') 'FourBesDerv:',
!     &          ngrid,acc,delq,dpi3,ract
        rhoout = (acc*delq*dpi3/ract)
      END IF

      RETURN

      END SUBROUTINE FourBesDerv
