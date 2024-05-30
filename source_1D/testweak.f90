IMPLICIT NONE
INTEGER,PARAMETER :: ngrid=50
REAL(8),PARAMETER :: rstep=0.2D0
REAL(8),PARAMETER :: rdiffrin=6D0,sigmain=1D0,totalin=100D0
REAL(8),PARAMETER :: pi=3.14159265358979D0
REAL(8) :: rleng,delq,frpirs,q
REAL(8) :: dpi3,acc,c,s,qs,qa,x,r,ds,dc,qr

REAL(8) :: rhop(ngrid),rhon(ngrid),sjsp(ngrid),sjsn(ngrid),rhochr(ngrid),formf(ngrid)
INTEGER :: iq,ir

REAl(8) :: rmsw,rmsc,rmsp,rmsn,rms,  &
    rdiffrw,rdiffrc,rdiffrp,rdiffrn,rdiffr, surfw,surfc,surfp,surfn,surf


!----------------------------------------------------------------

!     prepare parameters

frpirs = 4D0*pi*rstep
rleng  = (ngrid-1)*rstep
delq   = pi/rleng

! prepare Helm model

formf(1) = totalin
DO iq=2,ngrid
  q = (iq-1)*delq
  qr = q*rdiffrin
  formf(iq)   = totalin*(SIN(qr)/qr-COS(qr))*3D0/(qr*qr) &
               *exp(-0.5D0*(sigmain*q)**2)
END DO

!     charge density by fourier back-transformation

dpi3   = 0.5D0/(pi*pi)
!                                                         r=0.0 separately
acc    = 0D0
DO iq=2,ngrid
  qa     = delq*iq-delq
  acc    = qa*qa*formf(iq)+acc
END DO
rhop(1) = (acc*delq*dpi3)
!                                                         now all other r
r      = 0D0
DO ir=2,ngrid
  r      = rstep+r
  
  dc     = COS(r*delq)
  ds     = SIN(r*delq)
  c      = 1D0
  s      = 0D0
  acc    = 0D0
  qa     = 0D0
  DO iq=2,ngrid
    qa     = delq+qa
    x      = c*dc-s*ds
    s      = s*dc+c*ds
    c      = x
    acc    = qa*s*formf(iq)+acc
  END DO
  rhop(ir) = (acc*delq*dpi3/r)
END DO

rhon = rhop


CALL radii(rhop,rhon,sjsp,sjsn,rstep,ngrid,0,0D0,  &
    rhochr, rmsw,rmsc,rmsp,rmsn,rms,  &
    rdiffrw,rdiffrc,rdiffrp,rdiffrn,rdiffr, surfw,surfc,surfp,surfn,surf)


STOP
END


