!-----ini_fyfunc------------ part of the hartree-fock package -----------------
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-28  Time: 21:54:07

SUBROUTINE ini_fyfunc()

INCLUDE 'params.inc'

!     Initialize auxiliary parameters for Fayans functional.

!     Parameters are input per common /skypar/,

INCLUDE 'skypar.inc'
INCLUDE 'constn.inc'

PARAMETER (onedtwo=1.0D0/2.0D0)
PARAMETER (onedsix=1.0D0/6.0D0)
PARAMETER (hc=197.3285879)
PARAMETER (dsmall=1D-30)

!------------------------------------------------------------------------------

!     switch to Fayans model

!      IF(.NOT.ioptsk.EQ.13) STOP "FYFUNC needs IOPTSK=13"

!     auxiliary variables

!          general
r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
r0fy2 = r0fy*r0fy
efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
efacfy = efermfy*rho0fy*(2D0/3D0)
ufacfy = efermfy/3D0
pfacfy = 2D0*efermfy/(3D0*rho0fy)
rho0inv = 0.5D0/rho0fy
!         spin-orbit
c0fac = 2D0*efermfy/(3D0*rho0fy)*r0fy2
clsp = c0fac*c0nabj
clsm = c0fac*c1nabj
ctensp = c0fac*c0tens
ctensm = c0fac*c1tens
!         pairing
!      ufacfypr = pfacfy/(2D0*rho0xi)
ufacfypr = pfacfy/(rho0xi)     ! ???
!      ufacfypr = efermfy/(3D0*rho0fy*rho0xi)
pgradfac = -2D0*hgradxi*r0fy2*pfacfy/rho0xi
!      pgradfac = 2D0*ufacfypr*hgradp*r0fy2
ipair = 7
tpair = .true.

!      WRITE(2,'(a,6(1pg16.7))') 'rho0fy,efacfy,efermfy,r0fy=',
!     & rho0fy,efacfy,efermfy

RETURN
END SUBROUTINE ini_fyfunc

!-----fyfunc------------ part of the hartree-fock package ---------------------

SUBROUTINE fyfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln, sp,sn,ssp,ssn,  &
    frp,frn,potp,potn,bmasp,bmasn,bmaslp,bmasln, wlsp,wlsn,energd)

INCLUDE 'params.inc'

!     the basic Fayans energy density and derived quantities.

!     Input per list:
!      rp,rn        = densities
!      rc           = charge density (used in Coulomb exchange)
!      d1p,d1n      = derivative of densities
!      d2p,d2n      = Laplacian of densities
!      tp,tn        = kinetic energy densities
!      tlp,tln      = L**2 densities
!      sp,sn        = spin orbit densities
!      ssp,ssn      = derivative of spin-orbit densities
!      frp,frn      = folded densities (for Yukawa term)

!     Forces are input per common /skypar/ und /tswitc/.
!     Ouput per list:
!      energd        = energy density
!      potp,potn     = local potentials
!      bmasp,bmasn   = inverse mass
!      bmaslp,bmasln = inverse L**2 mass
!      wlsp,wlsn     = spin-orbit potential

INCLUDE 'skypar.inc'
INCLUDE 'numpar.inc'
INCLUDE 'wavfun.inc'
INCLUDE 'constn.inc'


REAL, INTENT(IN)                         :: rp
REAL, INTENT(IN)                         :: rn
REAL, INTENT(IN)                         :: rc
REAL, INTENT(IN)                         :: d1p
REAL, INTENT(IN)                         :: d1n
REAL, INTENT(IN)                         :: d2p
REAL, INTENT(IN)                         :: d2n
REAL, INTENT(IN)                         :: tp
REAL, INTENT(IN)                         :: tn
REAL, INTENT(IN OUT)                     :: tlp
REAL, INTENT(IN OUT)                     :: tln
REAL, INTENT(IN)                         :: sp
REAL, INTENT(IN)                         :: sn
REAL, INTENT(IN)                         :: ssp
REAL, INTENT(IN)                         :: ssn
REAL, INTENT(IN OUT)                     :: frp
REAL, INTENT(IN OUT)                     :: frn
REAL, INTENT(OUT)                        :: potp
REAL, INTENT(OUT)                        :: potn
REAL, INTENT(OUT)                        :: bmasp
REAL, INTENT(OUT)                        :: bmasn
REAL, INTENT(IN OUT)                     :: bmaslp
REAL, INTENT(IN OUT)                     :: bmasln
REAL, INTENT(OUT)                        :: wlsp
REAL, INTENT(OUT)                        :: wlsn
REAL, INTENT(OUT)                        :: energd
REAL, PARAMETER :: onedtwo=1.0D0/2.0D0
REAL, PARAMETER :: onedsix=1.0D0/6.0D0
REAL, PARAMETER :: hc=197.3285879
REAL, PARAMETER :: dsmall=1D-30

LOGICAL,:: SAVE :: tfirst=.true.

!------------------------------------------------------------------------------

!     switch to Fayans model

IF(.NOT.ioptsk == 13) STOP "FYFUNC needs IOPTSK=13"

!     auxiliary variables

!      r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
!      r0fy2 = r0fy*r0fy
!      efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
!      efacfy = efermfy*rho0fy*(2D0/3D0)
!      ufacfy = efermfy/3D0
!      rho0inv = 0.5D0/rho0fy
!      clsp = 2D0*efermfy/(3D0*rho0fy)*r0fy2*c0nabj
!      clsm = 2D0*efermfy/(3D0*rho0fy)*r0fy2*c1nabj

!     total densities and other auxiliary quantities

rs  = rp+rn
ts  = tp+tn
ss  = sp+sn
sss = ssp+ssn
d2s = d2p+d2n
d1s = d1p+d1n
xplus = rs*rho0inv
xminus = (rp-rn)*rho0inv
xpsigma = MAX(xplus,dsmall)**sigmavol
xmsigma = xminus  !**sigmavol
d1sr2 = d1s*d1s*r0fy2*rho0inv**2
d1mr2 = d1sr2 ! (d1p-d1n)**2*r0fy2*rho0inv**2
!      energ1 =  efacfy*(avp*(xplus*xplus)*
!     &   (1D0-h1vp*xpsigma)/(1D0+h2vp*xpsigma))

!     energy  (direct Coulomb separately, pairing part in HFBpot)

!      IF(rc==zero) THEN
!        WRITE(911,'(a,8(1pg13.5))') 'FYFUNC:',
!     &    rs,1D0/rho0inv,rp-rn,avm,h1vm,h2vm,efacfy,
!     &    efacfy*avm*(xminus*xminus)*(1D0-h1vm*xplus)/
!     &    (1D0+h2vm*xplus)/rs
!      END IF
energd = efacfy*(avp*(xplus*xplus)* (1D0-h1vp*xpsigma)/(1D0+h2vp*xpsigma)  &
    +avm*(xminus*xminus)* (1D0-h1vm*xplus)/(1D0+h2vm*xplus)  &
    +asp*d1sr2/(1D0+hsp*xpsigma+hgradp*d1sr2)  &
    +asm*d1mr2/(1D0+hsm*xplus+hgradm*d1mr2) )  &
    +c0tau*rs*ts+c1tau*(rp-rn)*(tp-tn)
!      energ2=energd-energ1
!     &         +ctensp*ss*ss+ctensm*(sp-sn)**2
!     &         +clsp*rs*sss+clsm*(rp-rn)*(ssp-ssn)
!old     &                 (1D0-h1vm*xmsigma)/(1D0+h2vm*xmsigma)

!     mean field potentials (surface term needs extra loop later)
!                           (direct Coulomb separately, pairing part in HFBpot)

potp = ufacfy*( 2D0*avp*xplus*(1D0-h1vp*xpsigma)/(1D0+h2vp*xpsigma)  &
    -avp*sigmavol*xpsigma*xplus*(h1vp+h2vp)/ (1D0+h2vp*xpsigma)**2  &
    +2D0*avm*xminus*(1D0-h1vm*xplus)/(1D0+h2vm*xplus)  &
    -avm*xminus*xminus*(h1vm+h2vm)/(1D0+h2vm*xplus)**2  &
    -asp*hsp*sigmavol*d1sr2*xpsigma/MAX(xplus,dsmall)/  &
    (1D0+hsp*xpsigma+hgradp*d1sr2)**2  &
    -asm*hsm*d1mr2/(1D0+hsm*xplus+hgradm*d1mr2)**2 )  &
    +c0tau*ts+c1tau*(tp-tn)
!     &       +c0nabj*sss+c1nabj*(ssp-ssn)
potn = ufacfy*( 2D0*avp*xplus*(1D0-h1vp*xpsigma)/(1D0+h2vp*xpsigma)  &
    -avp*sigmavol*xpsigma*xplus*(h1vp+h2vp)/ (1D0+h2vp*xpsigma)**2  &
    -2D0*avm*xminus*(1D0-h1vm*xplus)/(1D0+h2vm*xplus)  &
    -avm*xminus*xminus*(h1vm+h2vm)/(1D0+h2vm*xplus)**2  &
    -asp*hsp*sigmavol*d1sr2*xpsigma/MAX(xplus,dsmall)/  &
    (1D0+hsp*xpsigma+hgradp*d1sr2)**2  &
    -asm*hsm*d1mr2/(1D0+hsm*xplus+hgradm*d1mr2)**2 )  &
    +c0tau*ts-c1tau*(tp-tn)
!     &       +c0nabj*sss-c1nabj*(ssp-ssn)
bmasp = c0tau*rs+c1tau*(rp-rn)
bmasn = c0tau*rs-c1tau*(rp-rn)

!     Coulomb exchange

IF(texcha) THEN
  excfac = 1.4179979
  hxpsigmac = hcoul*xplus**sigmacoul
  add0 = -excfac*rc**0.333333333333D0
  add = add0*(1D0-hxpsigmac)
  energd = 0.75D0*add*rc + energd
  potp = add - 0.75D0*sigmacoul*add0*hxpsigmac*rc*rho0inv/xplus + potp
  potn =  - 0.75D0*sigmacoul*add0*hxpsigmac*rc*rho0inv/xplus + potn
END IF
!      energ3=energd-energ2

!     spin-orbit

potp   = -clsp*sss - cslm*(ssp-ssn)   + potp
!new      potp   = -clsp*sss - clsm*(ssp-ssn)   + potp
potn   = -clsp*sss + clsm*(ssp-ssn)   + potn
wlsp   = (clsp*d1s+clsm*(d1p-d1n)) +ctensp*ss+ctensm*(sp-sn)
wlsn   = (clsp*d1s-clsm*(d1p-d1n)) +ctensp*ss-ctensm*(sp-sn)
energd = -clsp*rs*sss - clsm*(rp-rn)*(ssp-ssn)  &
    +ctensp*ss*ss+ctensm*(sp-sn)**2 + energd
!      energ4=energd-energ3

!
IF(tfirst) THEN
  WRITE(6,'(a,6(1pg13.5))') 'r0fy,efermfy,C_0,ufacfy,rhoinv=',  &
      r0fy,efermfy,2D0*efermfy/(3D0*rho0fy),ufacfy,rho0inv,asp
  WRITE(6,'(a,4(1pg13.5))') 'potp etc=', rp,rn,xpsigma,potp,potn,energd
  WRITE(6,*) 'clsp,c0nabj,clsm,c1nabj=',clsp,c0nabj,clsm,c1nabj
  tfirst=.false.
END IF
!      WRITE(911,'(a,8(1pg16.7))')'rs,xplus,rho0inv,E=',rs,xplus,rho0inv,
!     &   energ1/rs,energ2/rs,energ3/rs,energ4/rs

RETURN
END SUBROUTINE fyfunc
!-----fyfunc_grad------------ part of the hartree-fock package -----------------

SUBROUTINE fyfunc_grad()

INCLUDE 'params.inc'

!     The derivative contribution to the gradient term of the Fayans functional-

!     Parameters are input per common /skypar/,
!     potentials and densities per /wavfun/.

!     The result is accumulated onto the array 'pot'.

INCLUDE 'skypar.inc'
INCLUDE 'numpar.inc'
INCLUDE 'wavfun.inc'
INCLUDE 'grid.inc'
INCLUDE 'constn.inc'

REAL(:: 8), allocatable :: dpotaux(:),gpotadd(:)
PARAMETER (onedtwo=1.0D0/2.0D0)
PARAMETER (onedsix=1.0D0/6.0D0)
PARAMETER (hc=197.3285879)
PARAMETER (dsmall=1D-30)

!------------------------------------------------------------------------------

allocate(dpotaux(1:kgr2),gpotadd(kgr2))
!      r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
!      r0fy2 = r0fy*r0fy
!      efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
!      efacfy = efermfy*rho0fy*(2D0/3D0)
!      ufacfy = efermfy/3D0
!      rho0inv = 0.5D0/rho0fy
!      WRITE(6,*) 'FYFUNC_GRAD:',r0fy,r0fy2,efermfy,efacfy,
!     &         efacfy*asp*r0fy2
DO i=1,ngrid
  j   = i+ngrid
  rp  = rho(i)
  rn  = rho(j)
  xminus = (rp-rn)*rho0inv
  xmsigma = xminus !**sigmavol
  rs  = rp+rn
  xplus = rs*rho0inv
  xpsigma = xplus**sigmavol
  d1p = d1r(i)
  d1n = d1r(j)
  d1mr2 = (d1p-d1n)**2*r0fy2*rho0inv**2
  d1s = d1p+d1n
  d1sr2 = d1s*d1s*r0fy2*rho0inv**2
  dpotaux(i) = ufacfy*asp*r0fy2*(1D0+hsp*xpsigma)/  &
      (1D0+hsp*xpsigma+hgradp*d1sr2)**2
  dpotaux(j) = ufacfy*asm*r0fy2*(1D0+hsm*xplus)/  &
      (1D0+hsm*xplus+hgradm*d1sr2)**2
!
!          --> pairing contribution still missing
!
!        WRITE(6,*) i,dpotaux(i),dpotaux(j)
END DO
CALL derv15(dpotaux(1),gpotadd(1), one,ngrid,rstep)
CALL derv15(dpotaux(ngridp),gpotadd(ngridp), one,ngrid,rstep)
!      WRITE(6,*) 'surface pots:'
DO i=1,ngrid
  j   = i+ngrid
  d1p = d1r(i)
  d1n = d1r(j)
  d2p = d2r(i)
  d2n = d2r(j)
  d1mr = (d1p-d1n)*rho0inv
  d1sr = (d1p+d1n)*rho0inv
  d2mr = (d2p-d2n)*rho0inv
  d2sr = (d2p+d2n)*rho0inv
  pot(i) = pot(i) -2D0*(gpotadd(i)*d1sr+dpotaux(i)*d2sr)  &
      -2D0*(gpotadd(j)*d1mr+dpotaux(j)*d2mr)
  pot(j) = pot(j) -2D0*(gpotadd(i)*d1sr+dpotaux(i)*d2sr)  &
      +2D0*(gpotadd(j)*d1mr+dpotaux(j)*d2mr)
!        WRITE(6,*) (i-1)*rstep,
!     &       dpotaux(i),-2D0*(gpotadd(i)*d1sr+dpotaux(i)*d2sr)
END DO

deallocate(dpotaux,gpotadd)
RETURN
END SUBROUTINE fyfunc_grad
!-----fy_hfbpot------------ part of the hartree-fock package -----------------

SUBROUTINE fy_hfbpot()

INCLUDE 'params.inc'

!     The pairing potential of the Fayans functional.
!     Also used to compute pairing energy.

!     Parameters are input per common /skypar/,
!     potentials and densities per /wavfun/.

!     The result is accumulated onto the array 'pot'.

INCLUDE 'skypar.inc'
INCLUDE 'numpar.inc'
INCLUDE 'wavfun.inc'
INCLUDE 'result.inc'
INCLUDE 'grid.inc'
INCLUDE 'constn.inc'

!------------------------------------------------------------------------------

!      WRITE(6,*) 'FY_HFBPOT'
!      r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
!      r0fy2 = r0fy*r0fy
!      efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
!      pfacfy = 2D0*efermfy/(3D0*rho0fy)
!      ufacfypr = pfacfy/(2D0*rho0xi)
epair = 0D0
DO i=1,ngrid
  ract = (i-1)*rstep
  in = i+ngrid
  d1s = d1r(i)+d1r(in)
  d1sr2 = (d1s*r0fy/rho0xi)**2
  xplus = (rho(i)+rho(in))/rho0xi
  xplusgamma = xplus**gammaxi
  rhfac = pfacfy*(fxi+h1xi*xplusgamma+hgradxi*d1sr2)
  pairpt(i)  = paircr(i)  * rhfac
  pairpt(in) = paircr(in) * rhfac
  paircrtot2 = paircr(i)**2+paircr(in)**2
!       epair = epair + 0.5D0*rhfac*paircrtot2*(ract*ract)
  epair = epair + rhfac*paircrtot2*(ract*ract)
END DO
epair = frpirs*epair
!      WRITE(6,'(a,4(1pg13.5))') 'FYpair:',pfacfy*fxi,pfacfy*h1xi,rho0xi,
!     &     epair

RETURN
END SUBROUTINE fy_hfbpot
!-----fy_pairpot------------ part of the hartree-fock package -----------------

SUBROUTINE fy_pairpot()

INCLUDE 'params.inc'

!     Non-gradient part of the  pairing contribnution
!     to potential from the Fayans functional.

!     Parameters are input per common /skypar/,
!     potentials and densities per /wavfun/.

!     The result is accumulated onto the array 'pot'.

INCLUDE 'skypar.inc'
INCLUDE 'numpar.inc'
INCLUDE 'wavfun.inc'
INCLUDE 'grid.inc'
INCLUDE 'constn.inc'
REAL(:: 8), allocatable :: dpotaux(:),gpotadd(:)

!------------------------------------------------------------------------------

allocate(dpotaux(1:kgr2),gpotadd(kgr2))
!      r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
!      r0fy2 = r0fy*r0fy
!      efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
!      ufacfypr = efermfy/(3D0*rho0fy*rho0xi)
!      pgradfac = -2D0*ufacfypr*hgradp*r0fy2
DO i=1,ngrid
  dpotaux(i) = (paircr(i)**2+paircr(i+ngrid)**2)
END DO
!      WRITE(6,*) 'paircr**2:',SUM(dpotaux(1:ngrid))
CALL derv15(dpotaux(1),gpotadd(1), one,ngrid,rstep)
!      WRITE(6,*) 'FY_PAIRPOT test pairing potential:',
!     &   ufacfypr*h1xi*gammaxi,pgradfac
DO i=1,ngrid
  j = i+ngrid
  d1sr = (d1r(i)+d1r(j))/rho0xi
  d1sr2 = (d1sr*r0fy)**2
  d2sr = (d2r(i)+d2r(j))/rho0xi
  xplus = (rho(i)+rho(j))/rho0xi
  IF(gammaxi==1D0) THEN
    xplusgamma = 1D0
  ELSE
    xplusgamma = xplus**(gammaxi-1D0)
  END IF
  add  = ufacfypr*dpotaux(i)*h1xi*gammaxi*xplusgamma  &
      +pgradfac*(dpotaux(i)*d2sr+gpotadd(i)*d1sr)
!        WRITE(6,*) i,add,dpotaux(i),gpotadd(i)
!        WRITE(6,*) i,add,pot(i),pot(j)
  pot(i) = pot(i) + add
  pot(j) = pot(j) + add
!        ppadd(i) = add
!        ppadd(j) = xplusgamma
END DO
!      WRITE(6,'(5(1pg13.5))') ppadd(1:ngrid)

deallocate(dpotaux,gpotadd)
RETURN
END SUBROUTINE fy_pairpot
!-----attach_fy--------------------------------------------------------------------

SUBROUTINE attach_fy(par,jpact)

use nmp_map

!     the explicit Fayans model parameters
!     are attached according to the values on the field 'par'
!     and the pointers on the field 'jp' .

INCLUDE 'params.inc'

INCLUDE 'skypar.inc'
INCLUDE 'numpar.inc'
INCLUDE 'comfit.inc'


REAL, INTENT(IN OUT)                     :: par(kparmx)
INTEGER, INTENT(IN OUT)                  :: jpact(kparmx)


DATA zero/0D0/

LOGICAL :: :: tfirst=.true.,tmapsurf=.false.

!--------------------------------------------------------------------

IF(ioptsk /= 13) STOP "ATTACH_FY requires IOPTSK=13"

!     first deal with the NMP entry

IF(iftinp==4) THEN
  
  eoain   = par(jpact(1))
  rhoin   = par(jpact(4))
  compin = par(jpact(11))
  asymin = par(jpact(6))
  dasyin = par(jpact(9))
  h2vm   = par(jpact(12))
  efmain = par(jpact(2))
  skapin = par(jpact(7))
  sigmavol = par(jpact(10))
  rho0fy  = par(jpact(22))
  IF(ifdfms == 2) THEN
    amfc0(1) = 20.73553 ! 20.7339829978 !
    amfc0(2) = 20.73553 ! 20.7339829978 !
  ELSE IF(ifdfms == 3) THEN
    amfc0(1) = 20.749811 ! 20.74982   ! 20.7505   !20.747
    amfc0(2) = 20.721249 ! 20.72126   ! 20.7219   !20.719
  ELSE
    STOP "IFDFMS must be 2 or 3 for Fayans"
  END IF
!        WRITE(2,'(a,8(1pg13.5))') 'rho,eoa,comp,efma,asym,dasy,skap=',
!     &         rhoin,eoain,compin,efmain,asymin,dasyin,skapin,rho0fy
  CALL  nmp_fayans(avp,h1vp,h2vp,avm,h1vm,h2vm,sigmavol,  &
      c0tau,c1tau,rho0fy, rhoin,eoain,compin,efmain,asymin,  &
      dasyin,skapin,amfc0(1),amfc0(2))
!        WRITE(2,'(a,7(1pg13.5))') 'rho,eoa,comp,efma,asym,dasy,skap=',
!     &         rhoin,eoain,compin,efmain,asymin,dasyin,skapin
!        WRITE(2,'(a,7(1pg13.5))') 'avp,h1vp,h2vp,avm,h1vm,h2vm=',
!     &    avp,h1vp,h2vp,avm,h1vm,h2vm
  IF(tmapsurf) h2vm = h2vm/avm
  
  asp    = par(jpact(3))
  asm    = par(jpact(8))
  hsp = par(jpact(23))
  hsm = par(jpact(24))
  hgradp = par(jpact(13))
  hgradm = par(jpact(14))
  IF(hsp<-1D3) hsp=h2vp
  IF(hsm<-1D3) hsm=h2vm    ! ??
  c0nabj = par(jpact(5))
  c1nabj = par(jpact(15))
  c0tens = par(jpact(36))       ! T=0 tensor coupling
  c1tens = par(jpact(37))       ! T=1 tensor coupling
  hcoul  = par(jpact(21))
  sigmacoul = 0.33333333333333D0
  fxi = par(jpact(16))
  h1xi = par(jpact(17))
  hgradxi = par(jpact(19))
  rho0xi = par(jpact(18))       ! default 0.16
  gammaxi = par(jpact(20))
  ecstab = par(jpact(30))
  IF(ABS(asp)>1D-10 .AND. tmapsurf) THEN
    hgradp = par(jpact(13))/asp
    IF(hsp>-1D3) THEN
      hsp = hsp/asp
    ELSE
      hsp=h2vp
    END IF
  END IF
  IF(ABS(asm)>1D-10 .AND. tmapsurf) THEN
    hgradm = par(jpact(14))/asm
    IF(hsm>-1D3) THEN
      hsm = hsm/asm
    ELSE
      hsm=h2vm    ! ??
    END IF
  END IF
  RETURN
END IF

!     now all the rest

avp    = par(jpact(1))
avm    = par(jpact(6))
sigmavol = par(jpact(10))
rho0fy  = par(jpact(22))
asp    = par(jpact(3))
asm    = par(jpact(8))
hsp = par(jpact(23))
hsm = par(jpact(24))
IF(iftinp /= 2) THEN
  h1vp   = par(jpact(4))
  h2vp   = par(jpact(11))
  h1vm   = par(jpact(9))
  h2vm   = par(jpact(12))
  hgradp = par(jpact(13))
  hgradm = par(jpact(14))
  IF(hsp<-1D3) hsp=h2vp
  IF(hsm<-1D3) hsm=h2vm    ! ??
  IF(tfirst) THEN
    WRITE(111,'(a,1pg16.8,a,1pg12.4)') 'h1vp=',h1vp*avp,',',  &
        dparmo(jpact(4))*ABS(avp)
    WRITE(111,'(a,1pg16.8,a,1pg12.4)') 'h2vp=',h2vp*avp,',',  &
        dparmo(jpact(11))*ABS(avp)
    WRITE(111,'(a,1pg16.8,a,1pg12.4)') 'h1vm=',h1vm*avm,',',  &
        dparmo(jpact(9))*ABS(avm)
    WRITE(111,'(a,1pg16.8,a,1pg12.4)') 'h2vm=',h2vm*avm,',',  &
        dparmo(jpact(12))*ABS(avm)
    WRITE(111,'(a,1pg16.8,a,1pg12.4)') 'hsp=',hsp*asp,',',  &
        dparmo(jpact(13))*ABS(asp)
    WRITE(111,'(a,1pg16.8,a,1pg12.4)') 'hgradp=',hgradp*asp,',',  &
        dparmo(jpact(14))*ABS(asp)
    tfirst = .false.
  END IF
ELSE
  h1vp   = par(jpact(4))/avp
  h2vp   = par(jpact(11))/avp
  h1vm   = par(jpact(9))/avm
  h2vm   = par(jpact(12))/avm
  IF(ABS(asp)>1D-10) THEN
    hgradp = par(jpact(13))/asp
    IF(hsp>-1D3) THEN
      hsp = hsp/asp
    ELSE
      hsp=h2vp
    END IF
  END IF
  IF(ABS(asm)>1D-10) THEN
    hgradm = par(jpact(14))/asm
    IF(hsm>-1D3) THEN
      hsm = hsm/asm
    ELSE
      hsm=h2vm    ! ??
    END IF
  END IF
END IF
c0nabj = par(jpact(5))
c1nabj = par(jpact(15))
c0tau  = par(jpact(2))
c1tau  = par(jpact(7))
c0tens = par(jpact(36))       ! T=0 tensor coupling
c1tens = par(jpact(37))       ! T=1 tensor coupling
hcoul  = par(jpact(21))
sigmacoul = 0.33333333333333D0
fxi = par(jpact(16))
h1xi = par(jpact(17))
hgradxi = par(jpact(19))
rho0xi = par(jpact(18))       ! default 0.16
gammaxi = par(jpact(20))
ecstab = par(jpact(30))

RETURN
END SUBROUTINE attach_fy
!-----naming_fy----------------------------------------------------------

SUBROUTINE naming_fy(zuord)

!     genererates array 'zuord' which carries the names of the
!     associated Fayans parameters

!      zuord    = output array of names

INCLUDE 'params.inc'
INCLUDE 'constn.inc'

INCLUDE 'comfit.inc'



CHARACTER (LEN=6), INTENT(OUT)           :: zuord(kparmx)


zuord(:) = '      '

zuord(1)='avp   '
zuord(4)='h1vp  '
zuord(11)='h2vp  '
zuord(6)='avm   '
zuord(9)='h1vm  '
zuord(12)='h2vm  '
zuord(10)='sigmav'
zuord(22)='rho0fy'
zuord(3)='asp   '
zuord(23)='hsp   '
zuord(13)='hgradp'
zuord(8)='asm   '
zuord(24)='hsm   '
zuord(14)='hgradm'
zuord(5)='c0nabj'
zuord(15)='c1nabj'
zuord(2)='c0tau '
zuord(7)='c1tau '
zuord(36)='c0tens'
zuord(37)='c1tens'
zuord(21)='hcoul '
zuord(16)='fxi   '
zuord(17)='h1xi  '
zuord(19)='hgradx'
zuord(18)='rho0xi'
zuord(20)='gammax'
zuord(30)='ecstab'
IF(iftinp==4) THEN
  zuord(1)='EoverA'
  zuord(4)='rho_nm'
  zuord(11)='compr '
  zuord(6)='asymm '
  zuord(9)='dasym '
  zuord(2)='effmas'
  zuord(7)=' kappa'
END IF

RETURN
END SUBROUTINE naming_fy
