c-----ini_fyfunc------------ part of the hartree-fock package -----------------
c
      subroutine ini_fyfunc()
c
      include 'params.inc'
c
c     Initialize auxiliary parameters for Fayans functional.
c
c     Parameters are input per common /skypar/, 
c
      include 'skypar.inc'
      include 'constn.inc'
c
      parameter (onedtwo=1.0D0/2.0D0)
      parameter (onedsix=1.0D0/6.0D0)
      parameter (hc=197.3285879)
      parameter (dsmall=1D-30)
c
c------------------------------------------------------------------------------
c
c     switch to Fayans model
c
!      IF(.NOT.ioptsk.EQ.13) STOP "FYFUNC needs IOPTSK=13"
c
c     auxiliary variables
c
!          general
      r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
      r0fy2 = r0fy*r0fy
      efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
      efacfy = efermfy*rho0fy*(2D0/3D0)
      ufacfy = efermfy/3D0
      pfacfy = 2D0*efermfy/(3D0*rho0fy)
      rho0inv = 0.5D0/rho0fy
!         spin-orbit
      C0fac = 2D0*efermfy/(3D0*rho0fy)*r0fy2
      clsp = C0fac*c0nabj
      clsm = C0fac*c1nabj
      ctensp = C0fac*c0tens
      ctensm = C0fac*c1tens
!         pairing
!      ufacfypr = pfacfy/(2D0*rho0xi)
      ufacfypr = pfacfy/(rho0xi)     ! ???
!      ufacfypr = efermfy/(3D0*rho0fy*rho0xi)
      pgradfac = -2D0*r0fy2*pfacfy/rho0xi
!      IF(ipair.NE.7 .AND. ipair.NE.8) 
!     &   STOP 'Fy functional requires IPAIR=7 or 8'
      tpair = .TRUE.

!      WRITE(2,'(a,6(1pg16.7))') 'rho0fy,efacfy,efermfy,pfacfy=',
!     &    rho0fy,efacfy,efermfy,pfacfy
!      CALL flush(02)

      return
      end

c-----fyfunc------------ part of the hartree-fock package ---------------------
c
      subroutine fyfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,
     &                  sp,sn,ssp,ssn,
     &                  frp,frn,potp,potn,bmasp,bmasn,bmaslp,bmasln,
     &                  wlsp,wlsn,energd)
c
      include 'params.inc'
c
c     the basic Fayans energy density and derived quantities.
c
c     Input per list:
c      rp,rn        = densities
c      rc           = charge density (used in Coulomb exchange)
c      d1p,d1n      = derivative of densities
c      d2p,d2n      = Laplacian of densities
c      tp,tn        = kinetic energy densities
c      tlp,tln      = L**2 densities
c      sp,sn        = spin orbit densities
c      ssp,ssn      = derivative of spin-orbit densities
c      frp,frn      = folded densities (for Yukawa term)
c                     
c     Forces are input per common /skypar/ und /tswitc/.
c     Ouput per list:
c      energd        = energy density
c      potp,potn     = local potentials
c      bmasp,bmasn   = inverse mass
c      bmaslp,bmasln = inverse L**2 mass
c      wlsp,wlsn     = spin-orbit potential
c
      include 'skypar.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'constn.inc'
c
      parameter (onedtwo=1.0D0/2.0D0)
      parameter (onedsix=1.0D0/6.0D0)
      parameter (hc=197.3285879)
      parameter (dsmall=1D-30)
      LOGICAL,SAVE :: tfirst=.true.
c
c------------------------------------------------------------------------------
c
c     switch to Fayans model
c
      IF(.NOT.ioptsk.EQ.13) STOP "FYFUNC needs IOPTSK=13"
c
c     auxiliary variables
c
!      r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
!      r0fy2 = r0fy*r0fy
!      efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
!      efacfy = efermfy*rho0fy*(2D0/3D0)
!      ufacfy = efermfy/3D0
!      rho0inv = 0.5D0/rho0fy
!      clsp = 2D0*efermfy/(3D0*rho0fy)*r0fy2*c0nabj
!      clsm = 2D0*efermfy/(3D0*rho0fy)*r0fy2*c1nabj
c
c     total densities and other auxiliary quantities
c
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
c
c     energy  (direct Coulomb separately, pairing part in HFBpot)
c
!      IF(rc==zero) THEN
!        WRITE(911,'(a,8(1pg13.5))') 'FYFUNC:',
!     &    rs,1D0/rho0inv,rp-rn,avm,h1vm,h2vm,efacfy,
!     &    efacfy*avm*(xminus*xminus)*(1D0-h1vm*xplus)/
!     &    (1D0+h2vm*xplus)/rs
!      END IF
      energd = efacfy*(avp*(xplus*xplus)*
     &                 (1D0-h1vp*xpsigma)/(1D0+h2vp*xpsigma)
     &                +avm*(xminus*xminus)*
     &                 (1D0-h1vm*xplus)/(1D0+h2vm*xplus)
     &                +asp*d1sr2/(1D0+hsp*xpsigma+hgradp*d1sr2)
     &                +asm*d1mr2/(1D0+hsm*xplus+hgradm*d1mr2)
     &         )
     &         +c0tau*rs*ts+c1tau*(rp-rn)*(tp-tn)
!      energ2=energd-energ1
!     &         +ctensp*ss*ss+ctensm*(sp-sn)**2
!     &         +clsp*rs*sss+clsm*(rp-rn)*(ssp-ssn)
cold     &                 (1D0-h1vm*xmsigma)/(1D0+h2vm*xmsigma)
c
c     mean field potentials (surface term needs extra loop later)
c                           (direct Coulomb separately, pairing part in HFBpot)
c
      potp = ufacfy*(
     &           2D0*avp*xplus*(1D0-h1vp*xpsigma)/(1D0+h2vp*xpsigma)
     &          -avp*sigmavol*xpsigma*xplus*(h1vp+h2vp)/
     &             (1D0+h2vp*xpsigma)**2
     &          +2D0*avm*xminus*(1D0-h1vm*xplus)/(1D0+h2vm*xplus)
     &          -avm*xminus*xminus*(h1vm+h2vm)/(1D0+h2vm*xplus)**2
     &          -asp*hsp*sigmavol*d1sr2*xpsigma/MAX(xplus,dsmall)/
     &               (1D0+hsp*xpsigma+hgradp*d1sr2)**2
     &          -asm*hsm*d1mr2/(1D0+hsm*xplus+hgradm*d1mr2)**2
     &       ) 
     &       +c0tau*ts+c1tau*(tp-tn)
!     &       +c0nabj*sss+c1nabj*(ssp-ssn)
      potn = ufacfy*(
     &           2D0*avp*xplus*(1D0-h1vp*xpsigma)/(1D0+h2vp*xpsigma)
     &          -avp*sigmavol*xpsigma*xplus*(h1vp+h2vp)/
     &             (1D0+h2vp*xpsigma)**2
     &          -2D0*avm*xminus*(1D0-h1vm*xplus)/(1D0+h2vm*xplus)
     &          -avm*xminus*xminus*(h1vm+h2vm)/(1D0+h2vm*xplus)**2
     &          -asp*hsp*sigmavol*d1sr2*xpsigma/MAX(xplus,dsmall)/
     &               (1D0+hsp*xpsigma+hgradp*d1sr2)**2
     &          -asm*hsm*d1mr2/(1D0+hsm*xplus+hgradm*d1mr2)**2
     &       ) 
     &       +c0tau*ts-c1tau*(tp-tn)
!     &       +c0nabj*sss-c1nabj*(ssp-ssn)
      bmasp = c0tau*rs+c1tau*(rp-rn)
      bmasn = c0tau*rs-c1tau*(rp-rn)
c
c     Coulomb exchange
c
      if(texcha) then
        excfac = 1.4179979
        hxpsigmac = hcoul*xplus**sigmacoul
        add0 = -excfac*rc**0.333333333333D0
        add = add0*(1D0-hxpsigmac)
        energd = 0.75D0*add*rc + energd
        potp = add - 0.75D0*sigmacoul*add0*hxpsigmac*rc*rho0inv/xplus
     &        + potp
        potn =  - 0.75D0*sigmacoul*add0*hxpsigmac*rc*rho0inv/xplus
     &        + potn
      endif
!      energ3=energd-energ2
c
c     spin-orbit
c
cold      potp   = -clsp*sss - cslm*(ssp-ssn)   + potp
      potp   = -clsp*sss - clsm*(ssp-ssn)   + potp
      potn   = -clsp*sss + clsm*(ssp-ssn)   + potn
      wlsp   = (clsp*d1s+clsm*(d1p-d1n)) 
     &        +ctensp*ss+ctensm*(sp-sn)
      wlsn   = (clsp*d1s-clsm*(d1p-d1n)) 
     &        +ctensp*ss-ctensm*(sp-sn)
      energd = -clsp*rs*sss - clsm*(rp-rn)*(ssp-ssn)
     &         +ctensp*ss*ss+ctensm*(sp-sn)**2
     &         + energd
!      energ4=energd-energ3

!
      IF(tfirst) THEN
        WRITE(6,'(a,6(1pg13.5))') 'r0fy,efermfy,C_0,ufacfy,rhoinv=',
     &        r0fy,efermfy,2D0*efermfy/(3D0*rho0fy),ufacfy,rho0inv,asp
        WRITE(6,'(a,4(1pg13.5))') 'potp etc=',
     &        rp,rn,xpsigma,potp,potn,energd
        WRITE(6,*) 'clsp,c0nabj,clsm,c1nabj=',clsp,c0nabj,clsm,c1nabj
        tfirst=.false.
      END IF
!      WRITE(911,'(a,8(1pg16.7))')'rs,xplus,rho0inv,E=',rs,xplus,rho0inv,
!     &   energ1/rs,energ2/rs,energ3/rs,energ4/rs
c
      return
      end
c-----fyfunc_grad------------ part of the hartree-fock package -----------------
c
      subroutine fyfunc_grad()
c
      include 'params.inc'
c
c     The derivative contribution to the gradient term of the Fayans functional-
c
c     Parameters are input per common /skypar/, 
c     potentials and densities per /wavfun/.
c
c     The result is accumulated onto the array 'pot'.
c
      include 'skypar.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'grid.inc'
      include 'constn.inc'
c
      REAL(8), ALLOCATABLE :: dpotaux(:),gpotadd(:)
      parameter (onedtwo=1.0D0/2.0D0)
      parameter (onedsix=1.0D0/6.0D0)
      parameter (hc=197.3285879)
      parameter (dsmall=1D-30)
c
c------------------------------------------------------------------------------
c
      ALLOCATE(dpotaux(1:KGR2),gpotadd(KGR2))
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
        dpotaux(i) = ufacfy*asp*r0fy2*(1D0+hsp*xpsigma)/
     &                   (1D0+hsp*xpsigma+hgradp*d1sr2)**2
        dpotaux(j) = ufacfy*asm*r0fy2*(1D0+hsm*xplus)/
     &                   (1D0+hsm*xplus+hgradm*d1sr2)**2
!
!          --> pairing contribution still missing
!
!        WRITE(6,*) i,dpotaux(i),dpotaux(j)
      END DO     
      CALL derv15(dpotaux(1),gpotadd(1), one,ngrid,rstep)
      CALL derv15(dpotaux(ngridp),gpotadd(ngridp), one,ngrid,rstep)
!      WRITE(6,*) 'surface pots:'
!      dpotaux(1:ngrid) = ufacfy*asp*r0fy2
!      dpotaux(ngrid+1:ngrid+ngrid) = 0D0
!      gpotadd=0D0
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
        pot(i) = pot(i)
     &          -2D0*(gpotadd(i)*d1sr+dpotaux(i)*d2sr)
     &          -2D0*(gpotadd(j)*d1mr+dpotaux(j)*d2mr)
        pot(j) = pot(j)
     &          -2D0*(gpotadd(i)*d1sr+dpotaux(i)*d2sr)
     &          +2D0*(gpotadd(j)*d1mr+dpotaux(j)*d2mr)
!        WRITE(6,*) (i-1)*rstep,
!     &       dpotaux(i),-2D0*(gpotadd(i)*d1sr+dpotaux(i)*d2sr)
      END DO     
c
      DEALLOCATE(dpotaux,gpotadd)
      return
      end
c-----fy_hfbpot------------ part of the hartree-fock package -----------------
c
      subroutine fy_hfbpot()
c
      include 'params.inc'
c
c     The pairing potential of the Fayans functional.
c     Also used to compute pairing energy.
c
c     Parameters are input per common /skypar/, 
c     potentials and densities per /wavfun/.
c
c     The result is accumulated onto the array 'pot'.
c
      include 'skypar.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'result.inc'
      include 'grid.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
!      WRITE(6,*) 'FY_HFBPOT'
!      r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
!      r0fy2 = r0fy*r0fy
!      efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
!      pfacfy = 2D0*efermfy/(3D0*rho0fy)
!      ufacfypr = pfacfy/(2D0*rho0xi)
      epair = 0D0
      epairn=0D0
      accgrad=0D0
      do i=1,ngrid
        ract = (i-1)*rstep
        in = i+ngrid
        d1s = d1r(i)+d1r(in)
        d1sr2 = (d1s*r0fy/rho0xi)**2
        xplus = (rho(i)+rho(in))/rho0xi
        xplusgamma = xplus**gammaxi
!        rhfac = pfacfy*(fxi+h1xi*xplusgamma+hgradxi*d1sr2)
        rhfacp = pfacfy*(fxip+h1xip*xplusgamma+hgradxip*d1sr2)
        rhfacn = pfacfy*(fxin+h1xin*xplusgamma+hgradxin*d1sr2)
        pairpt(i)  = paircr(i)  * rhfacp
        pairpt(in) = paircr(in) * rhfacn
!        paircrtot2 = paircr(i)**2+paircr(in)**2
!       epair = epair + 0.5D0*rhfac*paircrtot2*(ract*ract)
        epair = epair + (ract*ract)*   
     &          (rhfacp*paircr(i)**2+rhfacn*paircr(in)**2)
        epairn = epairn + (ract*ract)*rhfacn*paircr(in)**2
        accgrad=accgrad+(ract*ract)*d1sr2
      enddo
      epair = frpirs*epair
!      IF(nprot==20 .AND. nneut==20) THEN
!        WRITE(6,'(a,4(1pg13.5))') 'FYpair:',
!     &      pfacfy*fxip,pfacfy*h1xip,pfacfy*hgradxip,epair
!      END IF
        WRITE(6,'(a,4(1pg13.5))') 'FYpair:',epair,frpirs*epairn,
     &     frpirs*accgrad
!        write(6,'(a)') ' Fy pair potential for neutrons:'
!        write(6,'(5(1pg12.4))') (pairpt(i),i=ngrid+1,2*ngrid)

      return
      end
c-----fy_pairpot------------ part of the hartree-fock package -----------------
c
      subroutine fy_pairpot()
c
      include 'params.inc'
c
c     Non-gradient part of the  pairing contribnution 
c     to potential from the Fayans functional.
c
c     Parameters are input per common /skypar/, 
c     potentials and densities per /wavfun/.
c
c     The result is accumulated onto the array 'pot'.
c
      include 'skypar.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'grid.inc'
      include 'constn.inc'
      REAL(8), ALLOCATABLE :: dpotauxp(:),gpotaddp(:)
      REAL(8), ALLOCATABLE :: dpotauxn(:),gpotaddn(:)
c
c------------------------------------------------------------------------------
c
      ALLOCATE(dpotauxp(1:KGR2),gpotaddp(KGR2))
      ALLOCATE(dpotauxn(1:KGR2),gpotaddn(KGR2))
!      r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
!      r0fy2 = r0fy*r0fy
!      efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
!      ufacfypr = efermfy/(3D0*rho0fy*rho0xi)
      DO i=1,ngrid
        dpotauxp(i) = paircr(i)**2
        dpotauxn(i) = paircr(i+ngrid)**2
      END DO
!      WRITE(6,*) 'paircr**2:',SUM(dpotaux(1:ngrid))
      CALL derv15(dpotauxp(1),gpotaddp(1), one,ngrid,rstep)
      CALL derv15(dpotauxn(1),gpotaddn(1), one,ngrid,rstep)
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
        add  = ufacfypr*gammaxi*xplusgamma*
     &         (dpotauxp(i)*h1xip+dpotauxn(i)*h1xin)
     &        +pgradfac*(
     &           hgradxip*(dpotauxp(i)*d2sr+gpotaddp(i)*d1sr)
     &          +hgradxin*(dpotauxn(i)*d2sr+gpotaddn(i)*d1sr) )
!        WRITE(6,*) i,add,dpotaux(i),gpotadd(i)
!        WRITE(6,*) i,add,pot(i),pot(j)
!        add=0D0
        pot(i) = pot(i) + add
        pot(j) = pot(j) + add
!        ppadd(i) = add
!        ppadd(j) = xplusgamma
      END DO
!      WRITE(6,'(5(1pg13.5))') ppadd(1:ngrid)

      DEALLOCATE(dpotauxp,gpotaddp,dpotauxn,gpotaddn)
      return
      end
c-----attach_fy--------------------------------------------------------------------
c
      subroutine attach_fy(par,jpact)
c
      USE NMP_map
c
c     the explicit Fayans model parameters 
c     are attached according to the values on the field 'par'
c     and the pointers on the field 'jp' .
c
      include 'params.inc'
c
      include 'skypar.inc'
      include 'numpar.inc'
      include 'comfit.inc'
c
      dimension par(kparmx),jpact(kparmx)
c
      data zero/0D0/
      LOGICAL :: tfirst=.true.,tmapsurf=.false.
!     isovectormodel  switches option for isovector pairing
!                     0 -->  no isovector pairing
!                     1 -->  same factor on all terms
!                     2 -->  only on 'fxi' 
!      INTEGER,PARAMETER :: isovectormodel=1
c
c--------------------------------------------------------------------
c
      IF(ioptsk.NE.13) STOP "ATTACH_FY requires IOPTSK=13"
c
c     defaults
c
      delfxi=0D0
c
c     attach pairing parameters and recoujpled to p/n parameters
c
      fxi = par(jpact(16))
      delfxi = par(jpact(25))
      h1xi = par(jpact(17))
      hgradxi = par(jpact(19))
      rho0xi = par(jpact(18))       ! default 0.16
      gammaxi = par(jpact(20))
!
!               map pairing parameters according to pairing model
!
      IF(isovectormodel==0) THEN
        fxip=fxi
        h1xip=h1xi
        hgradxip=hgradxi
        fxin=fxi
        h1xin=h1xi
        hgradxin=hgradxi
      ELSE IF(isovectormodel==2) THEN
        fxip=fxi+delfxi
        h1xip=h1xi
        hgradxip=hgradxi
        fxin=fxi-delfxi
        h1xin=h1xi
        hgradxin=hgradxi
      ELSE IF(isovectormodel==1) THEN
        xipfac=(fxi+delfxi)/fxi
        fxip=fxi*xipfac
        h1xip=h1xi*xipfac
        hgradxip=hgradxi*xipfac
        xinfac=(fxi-delfxi)/fxi
        fxin=fxi*xinfac
        h1xin=h1xi*xinfac
        hgradxin=hgradxi*xinfac
      END IF
c
c     first deal with the NMP entry
c
      IF(iftinp==4) THEN
c
        eoain   = par(jpact(1))
        rhoin   = par(jpact(4))
        compin = par(jpact(11))
        asymin = par(jpact(6))
        dasyin = par(jpact(9))
        h2vm   = par(jpact(12))
        h2vm=MAX(h2vm,0D0)              ! avoid zeroes
        efmain = par(jpact(2))
        skapin = par(jpact(7))
        sigmavol = par(jpact(10))
        rho0fy  = par(jpact(22))
        IF(ifdfms.eq.2) THEN
          amfc0(1) = 20.73553 ! 20.7339829978 !
          amfc0(2) = 20.73553 ! 20.7339829978 !
        ELSE IF(ifdfms.eq.3) THEN
          amfc0(1) = 20.749811 ! 20.74982   ! 20.7505   !20.747
          amfc0(2) = 20.721249 ! 20.72126   ! 20.7219   !20.719
        ELSE
          STOP "IFDFMS must be 2 or 3 for Fayans"
        END IF
!        WRITE(2,'(a,8(1pg13.5))') 'rho,eoa,comp,efma,asym,dasy,skap=',
!     &         rhoin,eoain,compin,efmain,asymin,dasyin,skapin,rho0fy
        call  NMP_Fayans(avp,h1vp,h2vp,avm,h1vm,h2vm,sigmavol,
     &                   c0tau,c1tau,rho0fy,
     &                   rhoin,eoain,compin,efmain,asymin,
     &                   dasyin,skapin,amfc0(1),amfc0(2))
!        WRITE(2,'(a,7(1pg13.5))') 'rho,eoa,comp,efma,asym,dasy,skap=',
!     &         rhoin,eoain,compin,efmain,asymin,dasyin,skapin
!        WRITE(2,'(a,7(1pg13.5))') 'avp,h1vp,h2vp,avm,h1vm,h2vm=',
!     &    avp,h1vp,h2vp,avm,h1vm,h2vm
        IF(tmapsurf) h2vm = h2vm/avm
c
        asp    = par(jpact(3))
        asm    = par(jpact(8))
        hsp = par(jpact(23))
        hsm = par(jpact(24))
        hgradp = par(jpact(13))
!        hgradp=MAX(hgradp,0D0)            ! avoid zeroes
        hgradm = par(jpact(14))
        IF(hsp<-1D3) hsp=h2vp
        IF(hsm<-1D3) hsm=h2vm    ! ??
        c0nabj = par(jpact(5))
        c1nabj = par(jpact(15))
        c0tens = par(jpact(36))       ! T=0 tensor coupling
        c1tens = par(jpact(37))       ! T=1 tensor coupling
        hcoul  = par(jpact(21))
        sigmacoul = 0.33333333333333D0
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
c
c     now all the rest
c    
      avp    = par(jpact(1))
      avm    = par(jpact(6))
      sigmavol = par(jpact(10))
      rho0fy  = par(jpact(22))
      asp    = par(jpact(3))
      asm    = par(jpact(8))
      hsp = par(jpact(23))
      hsm = par(jpact(24))
      IF(iftinp.NE.2) THEN
        h1vp   = par(jpact(4))
        h2vp   = par(jpact(11))
        h1vm   = par(jpact(9))
        h2vm   = par(jpact(12))
        h2vm=MAX(h2vm,0D0)              ! avoid zeroes
        hgradp = par(jpact(13))
!        hgradp=MAX(hgradp,0D0)            ! avoid zeroes
        hgradm = par(jpact(14))
        IF(hsp<-1D3) hsp=h2vp
        IF(hsm<-1D3) hsm=h2vm    ! ??
        IF(tfirst) THEN
          WRITE(111,'(a,1pg16.8,a,1pg12.4)') 'h1vp=',h1vp*avp,',',
     &             dparmo(jpact(4))*ABS(avp)
          WRITE(111,'(a,1pg16.8,a,1pg12.4)') 'h2vp=',h2vp*avp,',',
     &             dparmo(jpact(11))*ABS(avp)
          WRITE(111,'(a,1pg16.8,a,1pg12.4)') 'h1vm=',h1vm*avm,',',
     &             dparmo(jpact(9))*ABS(avm)
          WRITE(111,'(a,1pg16.8,a,1pg12.4)') 'h2vm=',h2vm*avm,',',
     &             dparmo(jpact(12))*ABS(avm)
          WRITE(111,'(a,1pg16.8,a,1pg12.4)') 'hsp=',hsp*asp,',',
     &             dparmo(jpact(13))*ABS(asp)
          WRITE(111,'(a,1pg16.8,a,1pg12.4)') 'hgradp=',hgradp*asp,',',
     &             dparmo(jpact(14))*ABS(asp)
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
      ecstab = par(jpact(30))
c
      RETURN
      END
c-----naming_fy----------------------------------------------------------
c
      subroutine naming_fy(zuord)
c
c     genererates array 'zuord' which carries the names of the
c     associated Fayans parameters
c
c      zuord    = output array of names
c
      include 'params.inc'
      include 'constn.inc'
c
      include 'comfit.inc'
c
      character*6 zuord(kparmx)
c
      zuord(:) = '      '
c
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
       zuord(25)='delfxi'
       IF(iftinp==4) THEN
         zuord(1)='EoverA'
         zuord(4)='rho_nm'
         zuord(11)='compr '
         zuord(6)='asymm '
         zuord(9)='dasym '
         zuord(2)='effmas'
         zuord(7)=' kappa'
       END IF
c
      return
      end
c-----fyfunc_grad_linres------------ part of the hartree-fock package ---------
c
      subroutine fyfunc_grad_linres(potpart,facrhop,facrhon,
     &                              facd1rp,facd1rn,facd2rp,facd2rn)
c
      include 'params.inc'
c
c     The derivative contribution to the gradient term of the Fayans functional.
c     Version to deduce linear response with respect to gradient terms.
c
c     Parameters are input per common /skypar/, 
c     potentials and densities per /wavfun/.
c
c     List parameters:
c       potpart  = contribution from gradoent terms to local potential
c       facrho   = factor on density
c       facd1r   = factor on 1. derivative of density
c       facd2r   = factor on 2. derivative of density
c
c     The result is accumulated onto the array 'pot'.
c
      include 'skypar.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'grid.inc'
      include 'constn.inc'
c
      REAL(8),INTENT(IN) :: facrhop,facd1rp,facd2rp
      REAL(8),INTENT(IN) :: facrhon,facd1rn,facd2rn
      REAL(8),INTENT(OUT) :: potpart(kgr2)

      REAL(8), ALLOCATABLE :: dpotaux(:),gpotadd(:)
      parameter (onedtwo=1.0D0/2.0D0)
      parameter (onedsix=1.0D0/6.0D0)
      parameter (hc=197.3285879)
      parameter (dsmall=1D-30)
c
c------------------------------------------------------------------------------
c
      ALLOCATE(dpotaux(1:KGR2),gpotadd(KGR2))
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
        rp  = rho(i)*facrhop
        rn  = rho(j)*facrhon
        xminus = (rp-rn)*rho0inv
        xmsigma = xminus !**sigmavol
        rs  = rp+rn
        xplus = rs*rho0inv
        xpsigma = xplus**sigmavol
        d1p = d1r(i)*facd1rp
        d1n = d1r(j)*facd1rn
        d1mr2 = (d1p-d1n)**2*r0fy2*rho0inv**2
        d1s = d1p+d1n
        d1sr2 = d1s*d1s*r0fy2*rho0inv**2
        dpotaux(i) = ufacfy*asp*r0fy2*(1D0+hsp*xpsigma)/
     &                   (1D0+hsp*xpsigma+hgradp*d1sr2)**2
        dpotaux(j) = ufacfy*asm*r0fy2*(1D0+hsm*xplus)/
     &                   (1D0+hsm*xplus+hgradm*d1sr2)**2
!
!          --> pairing contribution still missing
!
!        WRITE(6,*) i,dpotaux(i),dpotaux(j)
      END DO     
      CALL derv15(dpotaux(1),gpotadd(1), one,ngrid,rstep)
      CALL derv15(dpotaux(ngridp),gpotadd(ngridp), one,ngrid,rstep)
!      WRITE(6,*) 'surface pots:'
!      dpotaux(1:ngrid) = ufacfy*asp*r0fy2
!      dpotaux(ngrid+1:ngrid+ngrid) = 0D0
!      gpotadd=0D0
      DO i=1,ngrid
        j   = i+ngrid
        d1p = d1r(i)*facd1rp
        d1n = d1r(j)*facd1rn
        d2p = d2r(i)*facd2rp
        d2n = d2r(j)*facd2rn
        d1mr = (d1p-d1n)*rho0inv
        d1sr = (d1p+d1n)*rho0inv
        d2mr = (d2p-d2n)*rho0inv
        d2sr = (d2p+d2n)*rho0inv
        potpart(i) = 
     &          -2D0*(gpotadd(i)*d1sr+dpotaux(i)*d2sr)
     &          -2D0*(gpotadd(j)*d1mr+dpotaux(j)*d2mr)
        potpart(j) = 
     &          -2D0*(gpotadd(i)*d1sr+dpotaux(i)*d2sr)
     &          +2D0*(gpotadd(j)*d1mr+dpotaux(j)*d2mr)
!        WRITE(6,*) (i-1)*rstep,
!     &       dpotaux(i),-2D0*(gpotadd(i)*d1sr+dpotaux(i)*d2sr)
      END DO     
c
      DEALLOCATE(dpotaux,gpotadd)
      return
      end
c-----fyfunc_grad_linres------------ part of the hartree-fock package ---------
c
      subroutine fyfunc_grad_linres2(potpart)
c
      include 'params.inc'
c
c     The derivative contribution to the gradient term of the Fayans functional.
c     Version to deduce linear response with respect to gradient terms.
c
c     List parameters:
c       potpart  = contribution from gradient terms to local potential
c
c     Present version not for the isovector gradient term!!!
c
      include 'skypar.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'grid.inc'
      include 'constn.inc'
      include 'linres.inc'
c
      REAL(8),INTENT(OUT) :: potpart(kgr2)

      REAL(8),PARAMETER :: reldiff=1.0D-4    ! relative step for derivatives 
      REAL(8), ALLOCATABLE :: dpotaux(:,:),gpotadd(:,:)
      parameter (onedtwo=1.0D0/2.0D0)
      parameter (onedsix=1.0D0/6.0D0)
      parameter (hc=197.3285879)
      parameter (dsmall=1D-30)
      LOGICAL,PARAMETER :: ttest=.FALSE.
c
c------------------------------------------------------------------------------
c
      IF(asm.NE.0D0) 
     &   STOP 'isovector gradient residual intercation not yet active'
!
      ALLOCATE(dpotaux(1:ngrid,0:1),gpotadd(1:ngrid,0:1))
c
c     prepare kernel fields
c
      alrpa=lrpa*1D0
      afac   = alrpa+alrpa*alrpa
      parity = 1-2*mod(lrpa,2)
      facrel=1D0+reldiff
      IF(ttest) WRITE(111,'(a)') 'kernel fields'
      DO i=2,ngrid
        j   = i+ngrid
        rs  = rho(i)+rho(j)
        d1s = d1r(i)+d1r(j)
        CALL fyfunc_gradpotential(gradmf0,rs,d1s)
        CALL fyfunc_gradpotential(gradmfr,rs*facrel,d1s)
        gradmfr=(gradmfr-gradmf0)/(rs*reldiff)
        CALL fyfunc_gradpotential(gradmfd,rs,d1s*facrel)
        gradmfd=(gradmfd-gradmf0)/(d1s*reldiff)
        dpotaux(i,0) = -2D0*gradmf0
        dpotaux(i,1) = -2D0*gradmfr*(prho(i)+prho(i+ngrid))
     &                 -2D0*gradmfd*(dprho(i)+dprho(i+ngrid))
        IF(ttest) WRITE(111,'(f6.2,4(1pg13.5))') (i-1)*rstep,
     &     dpotaux(i,0:1),gradmfr,gradmfd
      END DO     
c
c     derivative of kernel fields
c
      call extrap(dpotaux(1,0),one)
      call extrap(dpotaux(1,1),parity)
      CALL derv15(dpotaux(1,0),gpotadd(1,0), one,ngrid,rstep)
      CALL derv15(dpotaux(1,1),gpotadd(1,1), parity,ngrid,rstep)
c
c     compose to gradient contribution
c
      IF(ttest) WRITE(111,'(a)') 'contributions after derivative'
      DO i=2,ngrid
        j   = i+ngrid
        rinv   = rinvs(i)
        d1p = d1r(i)
        d1n = d1r(j)
        d2p = d2r(i)
        d2n = d2r(j)
!        d1mr = (d1p-d1n)*rho0inv
        d1sr = (d1p+d1n)*rho0inv
!        d2mr = (d2p-d2n)*rho0inv
        d2sr = (d2p+d2n)*rho0inv
        d1prs = dprho(i)+dprho(j)
        d2pr  = d2prho(i)+(2.0*dprho(i)-afac*prho(i)*rinv)*rinv
        d2pn  = d2prho(j)+(2.0*dprho(j)-afac*prho(j)*rinv)*rinv
        potpart(i) = 
     &          (gpotadd(i,1)*d1sr+dpotaux(i,1)*d2sr
     &          +gpotadd(i,0)*d1prs+dpotaux(i,0)*(d2pr+d2pn))
        potpart(j) = potpart(i)
        IF(ttest) WRITE(111,'(f6.2,4(1pg13.5))') (i-1)*rstep,
     &          gpotadd(i,1)*d1sr,dpotaux(i,1)*d2sr,
     &          gpotadd(i,0)*d1prs,dpotaux(i,0)*(d2pr+d2pn)
          dpotaux(i,0) = gpotadd(i,0)*d1sr+dpotaux(i,0)*d2sr
        IF(ttest)WRITE(111,'(f6.2,4(1pg13.5))') (i-1)*rstep,dpotaux(i,0)
      END DO     
      call extrap(dpotaux(1,0),one)
      CALL derv15(dpotaux(1,0),gpotadd(1,0), one,ngrid,rstep)
      DO i=2,ngrid
        WRITE(111,'(f6.2,4(1pg13.5))') (i-1)*rstep,gpotadd(i,0),
     &    potpart(i)
      END DO     
      WRITE(111,'(1x/1x)')
c
      DEALLOCATE(dpotaux,gpotadd)
      return
      end
c-----fyfunc-gradpotential------------ part of the hartree-fock package ------
c
      subroutine fyfunc_gradpotential(gradmf,rs,d1s)
c
      include 'params.inc'
c
c     The gradien-in-gradient contribution to the mean-field potential
c     of the Fayans functional.
c
c     Input list parameters:
c       rs     = total density
c       d1s    = 1. derivative of the total density
c     Output list parameters:
c       gradmf  = kernel of the isoscalar gradient mean field
c     Fayans parameters are communcated via 'skypar.inc'.
c
      include 'skypar.inc'
      include 'constn.inc'
      include 'linres.inc'
c
      REAL(8),INTENT(OUT) :: gradmf
      REAL(8),INTENT(IN) :: rs,d1s
c
      REAL(8) :: xpsigma,d1sr2
c
c-------------------------------------------------------------------------
c
        xpsigma = (rs*rho0inv)**sigmavol
        d1sr2 = d1s*d1s*r0fy2*rho0inv**2
c     
        gradmf = ufacfy*asp*r0fy2*(1D0+hsp*xpsigma)/
     &           (1D0+hsp*xpsigma+hgradp*d1sr2)**2
      return
      end
