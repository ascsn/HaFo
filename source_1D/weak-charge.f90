!-----radii    ----------------------------------------------------------
 

SUBROUTINE radii(rhop,rhon,sjsp,sjsn,rstep,ngrid,ifrso,ifrmf,cmwid,  &
    rhochr, rmsw,rmsc,rmsp,rmsn,rms,  &
    rdiffrw,rdiffrc,rdiffrp,rdiffrn,rdiffr, surfw,surfc,surfp,surfn,surf,&
    fow,str_squark)

!     Computes rms-radius, diffraction-radius and surface-thickness
!     for proton-, neutron-, total- and charge-density.
!     Densities are to be given on an equidistant 1D radial grid.

!     The charge-density includes folding with intrinsic electric and
!     magnetic formfactors of proton and neutron.
!     The contribution of the anomalous magnetic moment is included
!     only if ifrso=1.

!     Input variables are:
!      rhop    = array of radial proton density (rhop(1) --> at r=0)
!      rhon    = array of radial neutron density
!      sjsp    = Nabla of spin-orbit current for protons (used for ifrso=1)
!      sjsn    = Nabla of spin-orbit current for neutrons (used for ifrso=1)
!      rstep   = spacing of the (equidistant) radial grid
!      ngrid   = number of radial grid points
!      ifrso   = switch to l*s (or tensor) contribution to formfactor
!      ifrmf   = switch to RMF instead of SkyrmeHF
!      cmwid   = width for c.m. correction
!                 self-consistent c.m. correction uses
!                    cmwid  = 3/(8*<P^2>) = 0.375*(hbar^2/2M)/(ekin)
!                 rough estimate uses
!                    cmwid  = 1.58/(6.0*(nprot+nneut)**0.66666667)

!     Output variables are:
!      rhochr  = charge density inluding all nucleon formfactors
!      rmsw    = weak-charge r.m.s. radius
!      rmsc    = charge r.m.s. radius
!      rmsp    = proton r.m.s. radius
!      rmsn    = neutron r.m.s. radius
!      rms     = total r.m.s. radius
!      rdiffrc = weak-charge diffraction radius (sharp- or box-radius)
!      rdiffrc = charge diffraction radius (sharp- or box-radius)
!      rdiffrp = proton diffraction radius
!      rdiffrn = neutron diffraction radius
!      rdiffr  = total diffraction radius
!      surfw   = weak-charge surface thickness
!      surfc   = charge surface thickness
!      surfp   = proton surface thickness
!      surfn   = neutron surface thickness
!      surf    = total surface thickness

!     the following fields are workspace for formfactors:
!      foc    = complete charge formfactor, including also the folding with
!               the anomalous magnetic moment
!      fop    = formfactor of the proton distribution
!      fon    = formfactor of the neutorn distribution
!      fo     = formfactor of the total nucleon distribution


IMPLICIT NONE

REAL(8), INTENT(IN)         :: rhop(ngrid)
REAL(8), INTENT(IN)         :: rhon(ngrid)
REAL(8), INTENT(IN)         :: sjsp(ngrid)
REAL(8), INTENT(IN)         :: sjsn(ngrid)
REAL(8), INTENT(IN)         :: rstep
REAL(8), INTENT(IN)         :: str_squark
INTEGER, INTENT(IN)         :: ngrid
INTEGER, INTENT(IN)         :: ifrso
INTEGER, INTENT(IN)         :: ifrmf
REAL(8), INTENT(IN)         :: cmwid
REAL(8), INTENT(OUT)         :: rhochr(ngrid)
REAL(8), INTENT(OUT)         :: rmsw
REAL(8), INTENT(OUT)         :: rmsc
REAL(8), INTENT(OUT)         :: rmsp
REAL(8), INTENT(OUT)         :: rmsn
REAL(8), INTENT(OUT)         :: rms
REAL(8), INTENT(OUT)         :: rdiffrw
REAL(8), INTENT(OUT)         :: rdiffrc
REAL(8), INTENT(OUT)         :: rdiffrp
REAL(8), INTENT(OUT)         :: rdiffrn
REAL(8), INTENT(OUT)         :: rdiffr
REAL(8), INTENT(OUT)         :: surfw
REAL(8), INTENT(OUT)         :: surfc
REAL(8), INTENT(OUT)         :: surfp
REAL(8), INTENT(OUT)         :: surfn
REAL(8), INTENT(OUT)         :: surf



REAL(8) :: foc(ngrid),fop(ngrid),fon(ngrid),fo(ngrid),fow(ngrid)
REAL(8) :: fsop(ngrid),fson(ngrid)

REAL(8) :: delq,obs
REAL(8) :: fbinterpol
INTEGER :: iqmax,ir
INTEGER :: ifrsoa=0
INTEGER,PARAMETER :: ifplot=0

LOGICAL, PARAMETER :: tapv0=.TRUE. ! to compute the weak-charge asymmetry P_AV
LOGICAL :: tapv
!DATA  zero,half,one,pi/0D0,0.5D0,1D0,3.14159265358979D0/

!--------------------------------------------------------------------------

  WRITE(6,*) 'before FBTRANS'
  CALL FLUSH(6)

!     prepare basic formfactors

CALL  fbtrans(rhop,rhon,rstep,ngrid,cmwid,fop,fon,delq,iqmax)
IF(ifrso == 1) THEN
  CALL  fbtrans(sjsp,sjsn,rstep,ngrid,cmwid,fsop,fson,delq,iqmax)
ELSE
  fsop = 0D0
  fson = 0D0
END IF

!     compose electro-magnetic formfactor
ifrsoa=ifrso
WRITE(*,*) 'IFRSO...=',ifrso,ifrsoa,ifrmf,str_squark,cmwid
CALL formf_elmag(fop,fon,fsop,fson,delq,iqmax,ifrsoa,ifrmf,foc)

!     compose weak-charge formfactor

CALL formf_weak(fop,fon,fsop,fson,str_squark,delq,iqmax,ifrsoa,ifrmf,fow)

WRITE(*,*) '# test FOC,FOW:'
DO ir=1,iqmax
  WRITE(*,'(i5,f8.3,2(1pg13.5))') ir,(ir-1)*delq,foc(ir),fow(ir)
  CALL FLUSH(6)
END DO

!nprot = NINT(fop(1))
!nneut = NINT(fon(1))
tapv = tapv0 .AND. (NINT(fop(1))==82 .AND. NINT(fon(1))==126 &
                    .OR. NINT(fop(1))==20 .AND. NINT(fon(1))==28 )
WRITE(*,*) 'TAPV=',tapv,NINT(fop(1)),NINT(fon(1))
!     evaluate charge density and all form parameters
!     (charge density is used as work space in the first three calls)

WRITE(*,*) 'FOW:',fow
CALL FLUSH(6)
CALL form_parameters(fow,delq,iqmax,rstep,ngrid, rhochr,rmsw,rdiffrw,surfw)
IF(TAPV) THEN
 OPEN(131,file='dwc.dat')
 WRITE(131,'(0pf7.3,1pg20.12)') ((ir-1)*rstep,rhochr(ir),ir=1,ngrid)
 CLOSE(131)
END IF
WRITE(*,*) 'FOP:',fop
CALL FLUSH(6)
CALL form_parameters(fop,delq,iqmax,rstep,ngrid, rhochr,rmsp,rdiffrp,surfp)
WRITE(*,*) 'FON:',fon
CALL FLUSH(6)
CALL form_parameters(fon,delq,iqmax,rstep,ngrid, rhochr,rmsn,rdiffrn,surfn)
fo(:) = fop(:)+fon(:)
WRITE(*,*) 'FO:',fo
CALL FLUSH(6)
CALL form_parameters(fo,delq,iqmax,rstep,ngrid, rhochr,rms,rdiffr,surf)
WRITE(*,*) 'FOC:',foc
CALL FLUSH(6)
CALL form_parameters(foc,delq,iqmax,rstep,ngrid, rhochr,rmsc,rdiffrc,surfc)
IF(TAPV) THEN
 OPEN(131,file='dec.dat')
 WRITE(131,'(0pf7.3,1pg20.12)') ((ir-1)*rstep,rhochr(ir),ir=1,ngrid)
 CLOSE(131)
END IF




!     write protocol

WRITE(6,'(a)') ' form parameters from new subroutines:'
WRITE(06,'(16x,a,6(/1x,a,4(1pg13.5)))')  &
    ' diff.radius    surface     rms-radius   total nr.',  &
    ' charge density:',rdiffrc,surfc,rmsc,foc(1),  &
    '   mass density:',rdiffr,surf,rms,fo(1),  &
    ' proton density:',rdiffrp,surfp,rmsp,fop(1),  &
    'neutron density:',rdiffrn,surfn,rmsn,fon(1), &
    'weak-ch.density:',rdiffrw,surfw,rmsw,fow(1), &
    'diff. np,wc,wc0:',rmsn-rmsp,rmsw-rmsc,rmsw-5.503D0

obs = fbinterpol(0.475D0,fow,iqmax,delq,(ngrid-1)*rstep)
WRITE(6,'(a,2(1pg13.5))') ' F_W(0.457)=',obs,obs/fow(1)

IF(ifplot>0) THEN
  WRITE(11,'(a)') ' formfactors:','  q   F_p    F_n    F_W  '
  WRITE(11,'(4(1pg13.5))') ((ir-1)*delq,fop(ir),fon(ir),fow(ir),ir=1,ngrid)
END IF


RETURN
END SUBROUTINE radii

!-----formf_elmag----------------------------------------------------------

SUBROUTINE formf_elmag(fop,fon,fsop,fson,delq,iqmax,ifrso,ifrmf,foc)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: fop(iqmax)
REAL(8), INTENT(IN)             :: fon(iqmax)
REAL(8), INTENT(IN)             :: fsop(iqmax)
REAL(8), INTENT(IN)             :: fson(iqmax)
REAL(8), INTENT(IN)             :: delq
INTEGER, INTENT(IN)             :: iqmax
INTEGER, INTENT(IN)             :: ifrso
INTEGER, INTENT(IN)             :: ifrmf
REAL(8), INTENT(OUT)            :: foc(iqmax)

!     Composes the charge formfactor from given proton and neutron
!     formfactors.

!     Input variables are:
!      fop     = 1D array with proton formfactor
!      fon     = 1D array with neutron formfactor
!      fsop    = 1D array with proton spin-orbit formfactor
!      fson    = 1D array with neutron spin-orbit formfactor
!      delq    = spacing in momentum-space grid
!      iqmax   = maximum nr. of points in momentum-space grid
!      ifrso   = 1 --> charge density shall contain contribution from
!                      spin-orbit current (Skyrme) or tensor current (RMF)
!      ifrmf   = 1 --> for case of RMF, identify 'fop','fon' with vector
!                      density and 'fsop','fson' with tensor current.
!                0 --> SkyrmeHF
!     Output::
!      foc     = 1D array with charge formfactor




REAL(8) :: q1,GEpC,GEnC,GMpC,GMnC
INTEGER :: iq

!----------------------------------------------------------------

foc(1) = fop(1)
DO iq=2,iqmax
  q1 = delq*(iq-1)
  CALL formf_nucleon(q1,ifrso,ifrmf,GEpC,GEnC,GMpC,GMnC)
  foc(iq) = (fop(iq)*GEpC+fon(iq)*GEnC + GMpC*fsop(iq)+GMnC*fson(iq))
END DO

RETURN
END SUBROUTINE formf_elmag


!-----formf_weak----------------------------------------------------------

SUBROUTINE formf_weak(fop,fon,fsop,fson,str_squark,delq,iqmax,ifrso,ifrmf,fow)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: fop(iqmax)
REAL(8), INTENT(IN)             :: fon(iqmax)
REAL(8), INTENT(IN)             :: fsop(iqmax)
REAL(8), INTENT(IN)             :: fson(iqmax)
REAL(8), INTENT(IN)             :: str_squark
REAL(8), INTENT(IN)             :: delq
INTEGER, INTENT(IN)             :: iqmax
INTEGER, INTENT(IN)             :: ifrso
INTEGER, INTENT(IN)             :: ifrmf
REAL(8), INTENT(OUT)            :: fow(iqmax)

!     Composes the weak-charge formfactor from given proton and neutron
!     formfactors. Magnetic contributions are not included prresently.

!     Input variables are:
!      fop        = 1D array with proton formfactor
!      fon        = 1D array with neutron formfactor
!      str_squark = strength of strange-quark contribution
!      delq       = spacing in momentum-space grid
!      iqmax      = maximum nr. of points in momentum-space grid
!      ifrso   = 1 --> charge density shall contain contribution from
!                      spin-orbit current (Skyrme) or tensor current (RMF)
!      ifrmf   = 1 --> for case of RMF, identify 'fop','fon' with vector
!                      density and 'fsop','fson' with tensor current.
!                0 --> SkyrmeHF
!     Output::
!      fow        = 1D array with weak-charge formfactor


REAL(8),PARAMETER :: hbar2m=20.75D0,mc2=938.8D0  ! nucleon mass (roughly)
!REAL(8),PARAMETER :: sinthetaw2=0.23D0
!REAL(8),PARAMETER :: rad_corrn=0.9878            ! radiat.corr. neutrons
!REAL(8),PARAMETER :: rad_corrp=0.0721D0/(1D0-4D0*sinthetaw2)  ! "" protons

REAL(8) :: q1,GEpC,GEnC,GMpC,GMnC,darfac,gmz,gmn,gpz,gpn,gs,qfac,gsm
INTEGER :: iq
REAL(8) ::  Qweakp=0.0713D0        ! proton weak charge            
REAL(8) ::  Qweakn=-0.9888D0       ! neutron weak charge                
REAL(8) ::  xmgq=-0.017D0          ! magnetic moment of s-quark    

!----------------------------------------------------------------

!fow(1) = fop(1)
!WRITE(*,*) 'test FOW: str_squark=',str_squark
DO iq=1,iqmax
  q1 = delq*(iq-1)
  CALL formf_nucleon(q1,ifrso,ifrmf,GEpC,GEnC,GMpC,GMnC)
  qfac = q1*q1*0.5D0*hbar2m/mc2
  gs = str_squark*qfac/(1D0+4.97D0*qfac)
  gsm = xmgq*0.5D0*hbar2m/mc2
  gpz = Qweakp*GEpC+Qweakn*GEnC+Qweakn*gs
  gpn = Qweakn*GEpC+Qweakp*GEnC+Qweakn*gs
  gmz = Qweakp*GMpC+Qweakn*GMnC+Qweakn*gsm
  gmn = Qweakn*GMpC+Qweakp*GMnC+Qweakn*gsm
  fow(iq) = fop(iq)*gpz+fon(iq)*gpn + fsop(iq)*gmz+fson(iq)*gmn
!  gpz = rad_corrp*(0.25D0*(GEpC-GEnC)-sinthetaw2*GEpC-0.25D0*gs)
!  gpn = rad_corrn*(0.25D0*(GEnC-GEpC)-sinthetaw2*GEnC-0.25D0*gs)
!  gmz = rad_corrp*(0.25D0*(GMpC-GMnC)-sinthetaw2*GMpC)
!  gmn = rad_corrn*(0.25D0*(GMnC-GMpC)-sinthetaw2*GMnC)
!  fow(iq) = 4D0*(fop(iq)*gpz+fon(iq)*gpn + fsop(iq)*gmz+fson(iq)*gmn)
!WRITE(*,'(i5,5(1pg13.5))') iq,fsop(iq),fson(iq),schsmp,schsmn
END DO

RETURN
END SUBROUTINE formf_weak



!-----form_parameters-------------------------------------------------------

SUBROUTINE form_parameters(formf,delq,iqmax,rstep,ngrid,  &
    rhoout,rms,rdiffr,surf)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: formf(iqmax)
REAL(8), INTENT(IN)             :: delq
INTEGER, INTENT(IN)             :: iqmax
REAL(8), INTENT(IN)             :: rstep
INTEGER, INTENT(IN)             :: ngrid
REAL(8), INTENT(OUT)             :: rhoout(ngrid)
REAL(8), INTENT(OUT)             :: rms
REAL(8), INTENT(OUT)             :: rdiffr
REAL(8), INTENT(OUT)             :: surf

!     Computes the three form parameters for a given formfactor.
!     Input variables are:
!      formf   = 1D array with formfactor in steps of 'delq'
!      delq    = step size in momentum space
!      iqmax   = maximum nr. of grid points in momentum space
!      rstep   = spacing of the (equidistant) radial grid
!      ngrid   = number of radial grid points
!     output paramaters:
!      rhoout  = local density distribution from Fourier-Bessel back transf.
!      rms     = r.m.s. radius
!      rdiffr  = diffraction radius
!      surf    = surface thickness


REAL(8) :: findqzero,fbinterpol    ! function rooutines

REAL(8) :: q1,qr,finter,dpi3,qa,acc,r,dc,ds,c,s,x,rmom2
INTEGER :: iq,ir

REAL(8),PARAMETER :: pi=3.14159265358979D0


!------------------------------------------------------------------------------

!     diffraction radius

rdiffr  = 4.493D0/findqzero(formf,1,iqmax,delq,rstep*(ngrid-1))

!     surface thickness

q1     = 5.600D0/rdiffr
qr     = q1*rdiffr
finter = fbinterpol(q1,formf,iqmax,delq,rstep*(ngrid-1))
surf   = SQRT(ABS(2D0/(q1*q1)*LOG(ABS((SIN(qr)/qr-COS(qr))  &
    *3D0*formf(1)/(qr*qr*finter)))))

!     r.m.s. radius from density (derivative definition numerically unstable)
!     density from Fourier-Bessel transformation


!     charge density by fourier back-transformation

dpi3   = 0.5D0/(pi*pi)
!                                                         r=0.0 separately
acc    = 0D0
DO iq=2,iqmax
  qa     = delq*iq-delq
  acc    = qa*qa*formf(iq)+acc
END DO
rhoout(1) = (acc*delq*dpi3)
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
  DO iq=2,iqmax
    qa     = delq+qa
    x      = c*dc-s*ds
    s      = s*dc+c*ds
    c      = x
    acc    = qa*s*formf(iq)+acc
  END DO
  rhoout(ir) = (acc*delq*dpi3/r)
  
END DO

!     now r.m.s. radius by straightforward r-space integration

rmom2  = 0D0
r      = 0D0
DO ir=2,ngrid
  r     = rstep+r
  rmom2 = (r*r)*(r*r)*rhoout(ir) + rmom2
END DO
rmom2  = 4D0*pi*rstep*rmom2/formf(1)
rms    = SQRT(rmom2)

RETURN
END SUBROUTINE form_parameters



!-----formf_nucleon----------------------------------------------------------

SUBROUTINE formf_nucleon(q1,ifrso,ifrmf,GEpC,GEnC,GMpC,GMnC)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: q1
INTEGER, INTENT(IN)             :: ifrso
INTEGER, INTENT(IN)             :: ifrmf
REAL(8), INTENT(OUT)            :: GEpC
REAL(8), INTENT(OUT)            :: GEnC
REAL(8), INTENT(OUT)            :: GMpC
REAL(8), INTENT(OUT)            :: GMnC

!     Wrapper routine to call nucleon formfactor of choice 'imodel'.
!     The choice is set in this routine by a parameters statement.


!     Input variables are:
!      q1      = momentum at which fomrfactor is desired
!      ifrso   = 1 --> charge density shall contain contribtion from
!                      spin-orbit current
!                0 --> raw Sachs formfactors are used and no magnetic
!                      contribution; this option is the same in RMF
!                      and in SkyrmeHF.
!      ifrmf   = 1 --> prepares for use in RMF, variables are
!                      used in place  GEpC <--> G1pC, GMpC <--> G2pC
!                      and similarly for neutrons
!     Output::
!      GEpC  = electric proton formfactor (charge distribution)
!      GEnC  = electric neutron formfactor (charge distribution)
!      GMpC  = magnetic proton formfactor (charge distribution)
!      GMnC  = magnetic neutron formfactor (charge distribution)


!REAL(8) :: xmgne,xmgpr,d4m
!DATA xmgne,xmgpr,d4m/.04197,-.06131,.010989/
REAL(8),PARAMETER :: xmgpr = 2.793D0, xmgne=-1.913D0
!REAL(8),PARAMETER :: anmagp=xmgpr-1D0,anmagn=xmgne
REAL(8),PARAMETER :: d4m=0.011043D0,d2m=SQRT(d4m)

REAL(8) :: schsep,schsen,schsmp,schsmn,schses,schsev,schsm,darfac,tau

! switches formfactor model:
!  1  =  simple dipole
!  2  =  for protons Bernauer, for neutrons Simon&Walther
!  3  =  all Simon&Walther 
! 
INTEGER,PARAMETER :: imodel=3

!---------------------------------------------------------------------------

IF(imodel==1) THEN
  CALL formf_nucleon_dip(q1,schsep,schsen,schsmp,schsmn)
ELSE IF(imodel==2 ) THEN
  CALL formf_nucleon_SW(q1,schsep,schsen,schsmp,schsmn)
  CALL formf_nucleon_polyxdip(q1,schsep,schsmp)
ELSE IF(imodel==3) THEN
  CALL formf_nucleon_SW(q1,schsep,schsen,schsmp,schsmn)
ELSE
  STOP ' FORMF_NUCLEON: this IMODEL not provided'
END IF

IF(ifrso == 1) THEN
  tau = (q1*q1)*d4m
  IF(ifrmf==1) THEN
    GEpC = (schsep+xmgpr*tau*schsmp)/(1+tau)
    GEnC = (schsen+xmgne*tau*schsmn)/(1+tau)
    GMpC = (xmgpr*schsmp-schsep)/(1+tau)*d2m
    GMnC = (xmgne*schsmn-schsen)/(1+tau)*d2m
  ELSE
    darfac = 1D0/SQRT(1D0+tau)
    GMpC = -schsmp*xmgpr*(d4m+d4m) - schsep*d4m
    GMnC = -schsmn*xmgne*(d4m+d4m) - schsen*d4m
    GEpC = schsep*darfac
    GEnC = schsen*darfac
  END IF
ELSE
  GMpC = 0D0
  GMnC = 0D0
  GEpC = schsep
  GEnC = schsen
END IF


RETURN

END SUBROUTINE formf_nucleon

!-----formf_nucleon_dip--------------------------------------------------------

SUBROUTINE formf_nucleon_dip(q1,schsep,schsen,schsmp,schsmn)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: q1
REAL(8), INTENT(OUT)             :: schsep
REAL(8), INTENT(OUT)             :: schsen
REAL(8), INTENT(OUT)             :: schsmp
REAL(8), INTENT(OUT)             :: schsmn

!     Computes intrinsic proton and neutron formfactors in
!     single dipole form

REAL(8),PARAMETER :: rp2=0.769D0, rn2=-0.116D0
!REAL(8),PARAMETER :: xmgpr = 2.793D0, xmgne=-1.913D0

REAL(8) :: dipform,q2

!---------------------------------------------------------------------------

  q2 = q1**2
  dipform = (1D0+q2*rp2/12D0)**(-2)
  schsep = dipform
  schsen = -dipform*(q2*rn2/6D0)/(1D0+q2*0.04416D0)
  schsmp = dipform
  schsmn = dipform

RETURN

END SUBROUTINE formf_nucleon_dip


!-----formf_nucleon_SW--------------------------------------------------------

SUBROUTINE formf_nucleon_SW(q1,schsep,schsen,schsmp,schsmn)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: q1
REAL(8), INTENT(OUT)             :: schsep
REAL(8), INTENT(OUT)             :: schsen
REAL(8), INTENT(OUT)             :: schsmp
REAL(8), INTENT(OUT)             :: schsmn

!     Computes intrinsic proton and neutron formfactors in
!     Simon & Walther parameterization

!     Input variables are:
!      q1      = momentum at which fomrfactor is desired
!      ifrso   = 1 --> charge density shall contain contribtion from
!                      spin-orbit current
!     Output::
!      schsep  = electric proton formfactor
!      schsen  = electric neutron formfactor
!      schsmp  = magnetic proton formfactor
!      schsmn  = magnetic neutron formfactor
!      darfac  = Darwin correction (only for 'ifrso=1')


!      the parameters in data are for the sachs-formfactors:
!       sa..,sm.. are for the isoscalar electric formfactor
!       va..,vm.. are for the isovector electric formfactor
!       gp..,gm.. are for the magnetic formfactors of proton and neutron

!       isoscalar and isovector formfactors are recoupled to
!       the electric formfactors of proton and neutron.

!     the sachs-formfactors are then recoupled to the nucleon
!     vector formfactor (form1..) and the nucleon tensor
!     formfactor (form2..).


!     parameters for scalar and vector electric formfactor
!     (walther, priv.comm., sept. 86)

REAL(8) ::  sa1,sa2,sa3,sa4,sm1,sm2,sm3,sm4,  &
            va1,va2,va3,va4,vm1,vm2,vm3,vm4
DATA    sa1,    sa2,    sa3,   sa4,   sm1,  sm2,  sm3,  sm4  &
    /2.2907,-0.6777,-0.7923,0.1793, 15.75,26.68,41.01,134.2/
DATA    va1,    va2,    va3,   va4,   vm1,  vm2,  vm3,  vm4  &
    /0.3681, 1.2263,-0.6316,0.0372,  5.00,15.02,44.08,154.2/

!     magnetic formfactor , simon et al n.p. a333 (1980) 381

REAL(8) :: gp1,gp2,gp3,gp4,gm1,gm2,gm3,gm4  
DATA    gp1,    gp2,    gp3,   gp4,   gm1,  gm2,  gm3,  gm4  &
    / 0.694,  0.719,-0.418 ,0.005 ,  8.50,15.02,44.08,355.4/

!     old version of the
!     parameters for the magnetic formfactor (arenhoevel 1980)

!      data gp1,gp2,gp3,gp4,gm1,gm2,gm3,gm4/.794,.594,-.393,.005,
!     *                      8.73,15.02,44.08,355.40/

!     the anomanlous magnetic moments and the darwin factor 'd4m'
REAL(8) :: xmgne,xmgpr,d4m
DATA xmgne,xmgpr,d4m/.04197,-.06131,.010989/

REAL(8) :: schses,schsev,schsm,qa

!---------------------------------------------------------------------------

  
!      isoscalar and isovector electrical formactors on 'schses'
!      and 'schsev', recoupled then to proton and neutron electric
!      formfactors  'schsep' and 'schsen'. proton and neutron
!      magnetic formfactor on 'schsmp' and 'schsmn'.
  
qa     = q1*q1
schses = sa1/(1D0+qa/sm1)+sa2/(1D0+qa/sm2)+sa3/(1D0+qa/sm3)  &
      +sa4/(1D0+qa/sm4)
schsev = va1/(1D0+qa/vm1)+va2/(1D0+qa/vm2)+va3/(1D0+qa/vm3)  &
      +va4/(1D0+qa/vm4)
schsep = 0.5D0*(schses+schsev)
schsen = 0.5D0*(schses-schsev)
schsm  = ( gp1/(1D0+qa/gm1)+gp2/(1D0+qa/gm2)+gp3/(1D0+qa/gm3)  &
         +gp4/(1D0+qa/gm4) )
schsmp = schsm
schsmn = schsm


RETURN

END SUBROUTINE formf_nucleon_SW

!-----formf_nucleon_polyxdip---------------------------------------------------

SUBROUTINE formf_nucleon_polyxdip(q1,schsep,schsmp)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: q1
REAL(8), INTENT(OUT)             :: schsep
REAL(8), INTENT(OUT)             :: schsmp

!     Computes intrinsic proton and neutron formfactors in
!     polynomial x dipole form

!REAL(8),PARAMETER :: xmgpr = 2.793D0, xmgne=-1.913D0
REAL(8),PARAMETER :: ae(1:8)=(/-0.4980,+5.4592,-34.7281,+124.3173,-262.9808,+329.1395,-227.3306,+66.6980/)

REAL(8),PARAMETER :: am(1:8)=(/+0.2472,-4.9123,+29.7509,-84.0430,+129.3256,-111.1068,+49.9753,-9.1659/)


REAL(8) :: acc,q2
INTEGER :: i

!---------------------------------------------------------------------------

  q2=q1**2*(0.19731D0**2)
  acc=0D0
  DO i=8,1,-1
    acc = (acc+ae(i))*q2
  END DO
  schsep = (1D0+acc)/(1D0+q2/0.71D0)**2


  acc=0D0
  DO i=8,1,-1
    acc = (acc+am(i))*q2
  END DO
  schsmp =  (1D0+acc)/(1D0+q2/0.71D0)**2

RETURN

END SUBROUTINE formf_nucleon_polyxdip


!-----fbtrans ------------ part of the hartree-fock package ----------------------

SUBROUTINE fbtrans(rhop,rhon,rstep,ngrid,cmwid,fop,fon,delq,iqmax)

IMPLICIT NONE

REAL*8, INTENT(IN)              :: rhop(ngrid)
REAL*8, INTENT(IN)              :: rhon(ngrid)
REAL(8), INTENT(IN)             :: rstep
INTEGER, INTENT(IN)             :: ngrid
REAL(8), INTENT(IN)             :: cmwid
REAL(8), INTENT(OUT)             :: fop(ngrid)
REAL(8), INTENT(OUT)             :: fon(ngrid)
REAL(8), INTENT(OUT)             :: delq
INTEGER, INTENT(OUT)             :: iqmax

!     Input variables are:
!      rhop    = 1D array of radial proton density (rhop(1) --> at r=0)
!      rhon    = 1D array of radial neutron density
!      rstep   = spacing of the (equidistant) radial grid
!      ngrid   = number of radial grid points
!      cmwid   = width parameter for center of mass correction
!     output paramaters:
!      fop     = proton formfactor
!      fon     = neutron    "
!      delq    = spacing of the momentum-space grid
!      iqmax   = maximum number of grid points in momentum space


REAL(8) :: frpirs,rleng,q1,sump,sumn,r,dc,ds,s,c,x,exphf
INTEGER :: ir,iq

REAL(8),PARAMETER :: pi=3.14159265358979D0
!DATA  zero,half,one,pi/0D0,0.5D0,1D0,3.14159265358979D0/


!----------------------------------------------------------------

!     prepare parameters

frpirs = 4D0*pi*rstep
rleng  = (ngrid-1)*rstep
delq   = pi/rleng
iqmax  = ngrid

!     Fourier-Bessel transformation at q=0
sump   = 0D0
sumn   = 0D0
r      = 0D0
DO ir=2,ngrid
  r     = rstep+r
  sump   = r*r*rhop(ir)+sump
  sumn   = r*r*rhon(ir)+sumn
END DO
fop(1) = frpirs*sump
fon(1) = frpirs*sumn


!     Fourier-Bessel transformation for q>0

DO iq=2,iqmax
  q1     = iq*delq-delq
  dc     = COS(q1*rstep)
  ds     = SIN(q1*rstep)
  r      = 0D0
  sump   = 0D0
  sumn   = 0D0
  c      = 1D0
  s      = 0D0
  DO ir=2,ngrid
    r      = rstep+r
    x      = c*dc-s*ds
    s      = s*dc+c*ds
    c      = x
    sump   = r*s*rhop(ir)+sump
    sumn   = r*s*rhon(ir)+sumn
  END DO
  exphf  = EXP(cmwid*q1*q1)
  fop(iq) = exphf*frpirs*sump/q1
  fon(iq) = exphf*frpirs*sumn/q1
END DO

RETURN
END SUBROUTINE fbtrans
!-----fbinterpol -------------------------------------------------------

REAL(8) FUNCTION fbinterpol(qin,formfa,iqmax,delq,rlengbox)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: qin
REAL(8), INTENT(IN)             :: formfa(iqmax)
INTEGER, INTENT(IN)             :: iqmax
REAL(8), INTENT(IN)             :: delq
REAL(8), INTENT(IN)             :: rlengbox

!     Evaluates and returns formfactor at the momentum 'qin' by
!     fourier-bessel interpolation
!     from formfactors at fourier-grid points given on array 'formfa'
!     in the span '1..iqmax' and in steps of 'delq'.
!     Input variables:
!      qin     = momentum at which formfactor is wanted
!      formfa  = 1D array carrying formfactor at grid point in momentum space
!      iqmax   = maximum nr. of grid points in momentum space
!      delq   = grid spacing in momentum space
!      rlengbox = length of spatial numerical box



!     'precis' defines the neighbourhood of a pole-term for which the
!     routine switches to an explicit laurent expansion.
REAL(8),PARAMETER :: precis=1D-9

REAL(8),PARAMETER :: pi=3.14159265358979D0

REAL(8) :: qeff,xfx,xfx2,signum,acc,acc2,qdif,sinq
INTEGER :: ixfx,ipole,i

!-----------------------------------------------------------------------

qeff   = qin*rlengbox
xfx    = (qeff/pi)
xfx2   = xfx*xfx
ixfx   = (xfx+precis)

!     the fourier-bessel interpolation runs into problems with round-off
!     errors if 'qeff' comes too close to a multiple of 'pi'. we switch
!     to an explicit treatment of the pole-term in that case (i.e. first
!     case in the if-closure).

IF(ABS(xfx-ixfx) < precis) THEN
!                                             pole term separately
  ipole  = ixfx+1
  signum = -1D0
  acc    = 0D0
  DO i=2,iqmax
    signum = -signum
    IF(i /= ipole)  &
        acc    = formfa(i)*signum/(1D0-xfx2/FLOAT((i-1)*(i-1))) + acc
  END DO
  signum = (MOD(ipole,2)*2-1)
  qdif   = (xfx-ixfx)
  acc2   = signum*formfa(ipole)*(ixfx*ixfx)/((ixfx+ixfx)+qdif)
  qdif   = qdif*pi
  sinq   = (1D0-qdif*qdif*0.166666667)*signum
  fbinterpol  = 2.0*(acc*qdif+acc2*pi)*sinq/qeff
ELSE
!                                             normal case
  signum = -1D0
  acc    = 0D0
  DO i=2,iqmax
    signum = -signum
    acc    = formfa(i)*signum/(1D0-xfx2/FLOAT((i-1)*(i-1))) + acc
  END DO
  fbinterpol  = (acc+acc)*SIN(qeff)/qeff
END IF

RETURN
END FUNCTION fbinterpol
!-----findqzero ------------------------------------------------------

REAL(8) FUNCTION findqzero(foin,iqstrt,iqmax,delq,rlengbox)

IMPLICIT NONE

REAL(8), INTENT(IN)         :: foin(iqmax)
INTEGER, INTENT(IN)         :: iqstrt
INTEGER, INTENT(IN)         :: iqmax
REAL(8), INTENT(IN)         :: delq
REAL(8), INTENT(IN)         :: rlengbox

!     evaluates and returns the momentum at which the formfactor
!     given on array 'foin' has a zero.
!     scan for zero starts at the field position 'iqstrt'.
!     'iqmax' is the actual length of the array 'foin'.
!     'delq' is the stepsize in the momentum.
!     Input variables:
!      iqstrt  = grid point at which search starts
!      foin    = 1D array carrying formfactor at grid point in momentum space
!      iqmax   = maximum nr. of grid points in momentum space
!      delq    = grid spacing in momentum space
!      rlengbox = length of spatial numerical box



!     some numerical parameters for termination of the search
INTEGER,PARAMETER :: itq0mx=3
REAL(8),PARAMETER :: endq0=1D-4
!DATA itq0mx,endq0/3,0.0001/

REAL(8) :: fbinterpol        ! external function

INTEGER :: iqmaxm,iscan,iq1,iq2,itq0
REAL(8) :: q,q1,q2,q3,f1,f2,f3,q3old,qa1,qa2,c11,c2,c3,ca

!---------------------------------------------------------------------

!     scan changing sign in 'foin'.
!      'iq1' is position before the zero and 'iq1' after the zero.

iqmaxm = iqmax-1
DO iscan=iqstrt,iqmaxm
  IF(foin(iscan)*foin(iscan+1) <= 0D0) EXIT
END DO

iq1   = iscan
iq2   = iq1+1
IF(foin(iq1) <= 0D0) iq2   = iq1-1
q1    = iq1*delq-delq
q2    = iq2*delq-delq
f1    = foin(iq1)
f2    = foin(iq2)

!     come closer to 'foin'=0.0 by repeated linear, interpolation.

q3     = q1-f1*(q2-q1)/(f2-f1)
DO itq0=1,itq0mx
  f3     = fbinterpol(q3,foin,iqmax,delq,rlengbox)
  IF(f1*f3 > 0D0) THEN
    q1     = q3
    f1     = f3
  ELSE
    q2     = q3
    f2     = f3
  END IF
  q3old  = q3
  q3     = q1-f1*(q2-q1)/(f2-f1)
  IF(ABS(q3-q3old) < endq0) EXIT
END DO

!     improve position of the zero by quadratic interpolation:

!                         the coefficients of the interpolating binomial

f1     = fbinterpol(q1,foin,iqmax,delq,rlengbox)
f2     = fbinterpol(q2,foin,iqmax,delq,rlengbox)
f3     = fbinterpol(q3,foin,iqmax,delq,rlengbox)
ca     = 1D0/((q1-q2)*(q2-q3)*(q3-q1))
c11     = ca*(q3*(q3-q2)*q2*f1+q1*(q1-q3)*q3*f2+q2*(q2-q1)*q1*f3)
c2     = ca*((q2*q2-q3*q3)*f1+(q3*q3-q1*q1)*f2+(q1*q1-q2*q2)*f3)
c3     = ((q3-q2)*f1+(q1-q3)*f2+(q2-q1)*f3)*ca

!                         the two zeroes of the interpolating binomial

qa1    = ( SQRT(ABS(c2*c2-4D0*c11*c3))-c2)/(c3+c3)
qa2    = (-SQRT(ABS(c2*c2-4D0*c11*c3))-c2)/(c3+c3)

!                         select the zero which is closer to q3

q      = qa1
IF(ABS(qa1-q3) > ABS(qa2-q3)) q=qa2

!     return zero thus found.

findqzero  = q

RETURN
END FUNCTION findqzero
