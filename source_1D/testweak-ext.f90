IMPLICIT NONE
INTEGER,PARAMETER :: ngrid=201
REAL(8),PARAMETER :: rstep=0.1D0
REAL(8),PARAMETER :: rdiffrin=6D0,sigmain=1D0,totalin=100D0
REAL(8),PARAMETER :: pi=3.14159265358979D0
REAL(8) :: rleng,delq,frpirs,q
REAL(8) :: dpi3,acc,c,s,qs,qa,x,r,ds,dc,qr,cmwid,anumb

REAL(8) :: rhop(ngrid),rhon(ngrid),sjsp(ngrid),sjsn(ngrid),rhochr(ngrid),formf(ngrid)
INTEGER :: iq,ir

REAl(8) :: rmsw,rmsc,rmsp,rmsn,rms,  &
    rdiffrw,rdiffrc,rdiffrp,rdiffrn,rdiffr, surfw,surfc,surfp,surfn,surf


!----------------------------------------------------------------

!     prepare parameters

frpirs = 4D0*pi*rstep
rleng  = (ngrid-1)*rstep
delq   = pi/rleng

READ(5,'(1x)')
DO ir=1,ngrid
   READ(5,*) r,rhop(ir),rhon(ir)
END DO
sjsp=0D0
sjsn=0D0

anumb = 208D0
cmwid  = 1.58D0/(6D0*anumb**0.66666667D0)

CALL radii(rhop,rhon,sjsp,sjsn,rstep,ngrid,0,cmwid,  &
    rhochr, rmsw,rmsc,rmsp,rmsn,rms,  &
    rdiffrw,rdiffrc,rdiffrp,rdiffrn,rdiffr, surfw,surfc,surfp,surfn,surf)


!write(6,'(14x,a,4(/a,4(1pg13.4)))') &
!   '     rms       rdiffr     surface       halo',&
!   '   mass radii:',rms,rdiffr,surf,&
!    rms-sqrt(0.6*rdiffr*rdiffr+3.0*surf*surf),&
!   ' charge radii:',rmsc,rdiffrc,surfc,&
!    rmsc-sqrt(0.6*rdiffrc*rdiffrc+3.0*surfc*surfc)
!write(6,'(2(a,4(1pg13.4)/),2(a,2(1pg13.4)/))') &
!   ' proton radii:',rmsp,rdiffrp,surfp,&
!    rmsp-sqrt(0.6*rdiffrp*rdiffrp+3.0*surfp*surfp),&
!   ' neutr. radii:',rmsn,rdiffrn,surfn,&
!    rmsn-sqrt(0.6*rdiffrn*rdiffrn+3.0*surfn*surfn),&
!   ' diffe. radii:',rmsn-rmsp,rdiffrn-rdiffrp


STOP
END


