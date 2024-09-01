c-----ppotnew------------------------------------------------------
c
      subroutine ppotnew()
c
      include 'params.inc'
      parameter (reldiff=1.0D-4)    ! relative step for derivatives 
c      parameter (reldiff=1.0D-7)    ! relative step for derivatives 
      LOGICAL ttest
      PARAMETER (ttest=.false.)
c
c     the response potentials 'rpot', 'rpotls' and 'rbmass'
c     computed here by numerical derivatives of the mean.field 
c     given in 'skfunc'.
c
      include 'wavfun.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'potinv.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      REAL(8) :: frho(KGR2)        ! folded density in case of Yukawa
      REAL(8) :: gprho(KGR2,2)     ! gradient of repsonse density 'prho'
      REAL(8) :: unitr(KGR2,2)     ! unit vector of radial direction
      REAL(8),ALLOCATABLE :: potadd(:,:)
c
      data xminlw/1.0e-20/
c
c--------------------------------------------------------------------
c
c     preparatory steps
c
      parity = 1-2*mod(lrpa,2)
      fc     =1.0D0+reldiff           ! factor for variation
      afac   = alrpa+alrpa*alrpa
      do ir=1,ngrid
        work1(ir) = (ir-1)*rstep
      enddo
      call grad(work1,unitr,1)
      call grad(prho,gprho,2)
      call derv15(prho(1),dprho(1),parity,ngrid,rstep)
      call derv15(prho(ngrid+1),dprho(ngrid+1),parity,ngrid,rstep)
      call derv25(prho(1),d2prho(1),parity,ngrid,rstep)
      call derv25(prho(ngrid+1),d2prho(ngrid+1),parity,ngrid,rstep)
      if(tyukaw) then
        do i=1,ngrid
          frho(i)       = rho(i)                       ! use 'd1t' as workspace
          frho(i+ngrid) = rho(i+ngrid)
        enddo
        call lapfil(frho,work1,ayms2,-1)
        call lapfil(frho(ngridp),work1,ayms2,-1)
c        write(06,'(a,20(/1x,10g13.5))')
c     &        'testprint folded density:',(frho(i),i=1,ngrid)
c        write(06,'(20(/1x,10g13.5))')
c     &        (frho(i),i=ngridp,ngrid2)
      else
        do i=1,ngrid
          frho(i)       = 0.0D0
          frho(i+ngrid) = 0.0D0
        enddo
      endif
      t12co  = -0.125*(3.0*(t1+t1*0.5D0*x1)-t2-t2*0.5D0*x2)
      t12q   = t12co+0.125*(3.0*(0.5D0*t1+t1*x1)+0.5D0*t2+t2*x2)

c
c     loop over radial grid points
c

      do ir=1,ngrid
c
c      base point
c
       ic   = ir+ngrid
       rinv   = rinvs(ir)
       rp   = rho(ir)
       rn   = rho(ic)
       rc   = rp              !    not always true
       d1p  = d1r(ir)
       d1n  = d1r(ic)
       d2p  = d2r(ir)
       d2n  = d2r(ic)
       tp   = tau(ir)
       tn   = tau(ic)
       tlp  = taul(ir)
       tln  = taul(ic)
       sp   = sj(ir)
       sn   = sj(ic)
       ssp  = sjs(ir)
       ssn  = sjs(ic)
       frp  = frho(ir)
       frn  = frho(ic)
!       WRITE(9,'(a,i3)') 'IR:',ir
       call alfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,sp,sn,ssp,ssn,
     &             frp,frn,potp0,potn0,bmasp0,bmasn0,bmaslp0,bmasln0,
     &             wlsp0,wlsn0,energd)
c
c      variations, one after another
c
       call alfunc(rp*fc,rn,rc*fc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,sp,sn,
     &             ssp,ssn,frp*fc,frn,
     &             potp,potn,bmasp,bmasn,bmaslp,bmasln,wlsp,wlsn,energd)
       IF(rp<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = prho(ir)/(rp*reldiff)
       END IF
       rpot(ir) = (potp-potp0)*diffinv + rpot(ir)
       rpot(ic) = (potn-potn0)*diffinv + rpot(ic)
       rbmass(ir) = (bmasp-bmasp0)*diffinv + rbmass(ir)
       rbmass(ic) = (bmasn-bmasn0)*diffinv + rbmass(ic)
       rpotls(ir,1) = ((wlsp-wlsp0)*diffinv)*unitr(ir,1) + rpotls(ir,1)
       rpotls(ir,2) = ((wlsp-wlsp0)*diffinv)*unitr(ir,2) + rpotls(ir,2)
       rpotls(ic,1) = ((wlsn-wlsn0)*diffinv)*unitr(ir,1) + rpotls(ic,1)
       rpotls(ic,2) = ((wlsn-wlsn0)*diffinv)*unitr(ir,2) + rpotls(ic,2)
!       ap=(potp-potp0)*diffinv
!       an=(potn-potn0)*diffinv
c       
       call alfunc(rp,rn*fc,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,sp,sn,
     &             ssp,ssn,frp,frn*fc,
     &             potp,potn,bmasp,bmasn,bmaslp,bmasln,wlsp,wlsn,energd)
       IF(rn<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = prho(ic)/(rn*reldiff)
       END IF
       rpot(ir) = (potp-potp0)*diffinv + rpot(ir)
       rpot(ic) = (potn-potn0)*diffinv + rpot(ic)
       rbmass(ir) = (bmasp-bmasp0)*diffinv + rbmass(ir)
       rbmass(ic) = (bmasn-bmasn0)*diffinv + rbmass(ic)
       rpotls(ir,1) = ((wlsp-wlsp0)*diffinv)*unitr(ir,1) + rpotls(ir,1)
       rpotls(ir,2) = ((wlsp-wlsp0)*diffinv)*unitr(ir,2) + rpotls(ir,2)
       rpotls(ic,1) = ((wlsn-wlsn0)*diffinv)*unitr(ir,1) + rpotls(ic,1)
       rpotls(ic,2) = ((wlsn-wlsn0)*diffinv)*unitr(ir,2) + rpotls(ic,2)
!       WRITE(9,'(4(1pg13.5))')
!     &    ap+(potp-potp0)*diffinv,an+(potn-potn0)*diffinv
c       
       call alfunc(rp,rn,rc,d1p*fc,d1n,d2p,d2n,tp,tn,tlp,tln,sp,sn,
     &             ssp,ssn,frp,frn,
     &             potp,potn,bmasp,bmasn,bmaslp,bmasln,wlsp,wlsn,energd)
       IF(ABS(d1p)<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = dprho(ir)/(d1p*reldiff)
       END IF
       rpot(ir) = (potp-potp0)*diffinv + rpot(ir)
       rpot(ic) = (potn-potn0)*diffinv + rpot(ic)
       rbmass(ir) = (bmasp-bmasp0)*diffinv + rbmass(ir)
       rbmass(ic) = (bmasn-bmasn0)*diffinv + rbmass(ic)
       IF(ABS(d1p)<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = 1.0D0/(d1p*reldiff)
       END IF
       rpotls(ir,1) = (wlsp-wlsp0)*diffinv*gprho(ir,1) + rpotls(ir,1)
       rpotls(ir,2) = (wlsp-wlsp0)*diffinv*gprho(ir,2) + rpotls(ir,2)
       rpotls(ic,1) = (wlsn-wlsn0)*diffinv*gprho(ir,1) + rpotls(ic,1)
       rpotls(ic,2) = (wlsn-wlsn0)*diffinv*gprho(ir,2) + rpotls(ic,2)
!       IF(ioptsk==130) THEN
!         rpot(ir) = (potadd(ir,3)-potadd(ir,0))*diffinv + rpot(ir)
!         rpot(ic) = (potadd(ic,3)-potadd(ic,0))*diffinv + rpot(ic)
!       END IF
!       WRITE(9,'(4(1pg13.5))') potp-potp0,potn-potn0,
!     &      bmasp-bmasp0,bmasn-bmasn0
c       
       call alfunc(rp,rn,rc,d1p,d1n*fc,d2p,d2n,tp,tn,tlp,tln,sp,sn,
     &             ssp,ssn,frp,frn,
     &             potp,potn,bmasp,bmasn,bmaslp,bmasln,wlsp,wlsn,energd)
       IF(ABS(d1n)<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = dprho(ic)/(d1n*reldiff)
       END IF
       rpot(ir) = (potp-potp0)*diffinv + rpot(ir)
       rpot(ic) = (potn-potn0)*diffinv + rpot(ic)
       rbmass(ir) = (bmasp-bmasp0)*diffinv + rbmass(ir)
       rbmass(ic) = (bmasn-bmasn0)*diffinv + rbmass(ic)
       IF(ABS(d1n)<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = 1.0D0/(d1n*reldiff)
       END IF
       rpotls(ir,1) = (wlsp-wlsp0)*diffinv*gprho(ic,1) + rpotls(ir,1)
       rpotls(ir,2) = (wlsp-wlsp0)*diffinv*gprho(ic,2) + rpotls(ir,2)
       rpotls(ic,1) = (wlsn-wlsn0)*diffinv*gprho(ic,1) + rpotls(ic,1)
       rpotls(ic,2) = (wlsn-wlsn0)*diffinv*gprho(ic,2) + rpotls(ic,2)
c       
       call alfunc(rp,rn,rc,d1p,d1n,d2p*fc,d2n,tp,tn,tlp,tln,sp,sn,
     &             ssp,ssn,frp,frn,
     &             potp,potn,bmasp,bmasn,bmaslp,bmasln,wlsp,wlsn,energd)
       d2pr  = d2prho(ir)+(2.0*dprho(ir)-afac*prho(ir)*rinv)*rinv
       IF(ABS(d2p)<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = d2pr/(d2p*reldiff)
       END IF
       rpot(ir) = (potp-potp0)*diffinv + rpot(ir)
       rpot(ic) = (potn-potn0)*diffinv + rpot(ic)
       rbmass(ir) = (bmasp-bmasp0)*diffinv + rbmass(ir)
       rbmass(ic) = (bmasn-bmasn0)*diffinv + rbmass(ic)
       rpotls(ir,1) = ((wlsp-wlsp0)*diffinv)*unitr(ir,1) + rpotls(ir,1)
       rpotls(ir,2) = ((wlsp-wlsp0)*diffinv)*unitr(ir,2) + rpotls(ir,2)
       rpotls(ic,1) = ((wlsn-wlsn0)*diffinv)*unitr(ir,1) + rpotls(ic,1)
       rpotls(ic,2) = ((wlsn-wlsn0)*diffinv)*unitr(ir,2) + rpotls(ic,2)
       IF(ttest) WRITE(9,'(i3,4(1pg13.5))') ir,
     &    (potp-potp0)/(d2p*reldiff),(potn-potn0)/(d2p*reldiff),
     &   t12q,t12co
c       
       call alfunc(rp,rn,rc,d1p,d1n,d2p,d2n*fc,tp,tn,tlp,tln,sp,sn,
     &             ssp,ssn,frp,frn,
     &             potp,potn,bmasp,bmasn,bmaslp,bmasln,wlsp,wlsn,energd)
       d2pn    = d2prho(ic)+(2.0*dprho(ic)-afac*prho(ic)*rinv)*rinv
       IF(ABS(d2n)<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = d2pn/(d2n*reldiff)
       END IF
       rpot(ir) = (potp-potp0)*diffinv + rpot(ir)
       rpot(ic) = (potn-potn0)*diffinv + rpot(ic)
       rbmass(ir) = (bmasp-bmasp0)*diffinv + rbmass(ir)
       rbmass(ic) = (bmasn-bmasn0)*diffinv + rbmass(ic)
       rpotls(ir,1) = ((wlsp-wlsp0)*diffinv)*unitr(ir,1) + rpotls(ir,1)
       rpotls(ir,2) = ((wlsp-wlsp0)*diffinv)*unitr(ir,2) + rpotls(ir,2)
       rpotls(ic,1) = ((wlsn-wlsn0)*diffinv)*unitr(ir,1) + rpotls(ic,1)
       rpotls(ic,2) = ((wlsn-wlsn0)*diffinv)*unitr(ir,2) + rpotls(ic,2)
       IF(ttest) WRITE(9,'(i3,4(1pg13.5))') ir,
     &    (potp-potp0)/(d2n*reldiff),(potn-potn0)/(d2n*reldiff),
     &   t12co,t12q
c       
       call alfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp*fc,tn,tlp,tln,sp,sn,
     &             ssp,ssn,frp,frn,
     &             potp,potn,bmasp,bmasn,bmaslp,bmasln,wlsp,wlsn,energd)
       IF(tp<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = ptau(ir)/(tp*reldiff) 
       END IF
       rpot(ir) = (potp-potp0)*diffinv + rpot(ir)
       rpot(ic) = (potn-potn0)*diffinv + rpot(ic)
       rbmass(ir) = (bmasp-bmasp0)*diffinv + rbmass(ir)
       rbmass(ic) = (bmasn-bmasn0)*diffinv + rbmass(ic)
       rpotls(ir,1) = ((wlsp-wlsp0)*diffinv)*unitr(ir,1) + rpotls(ir,1)
       rpotls(ir,2) = ((wlsp-wlsp0)*diffinv)*unitr(ir,2) + rpotls(ir,2)
       rpotls(ic,1) = ((wlsn-wlsn0)*diffinv)*unitr(ir,1) + rpotls(ic,1)
       rpotls(ic,2) = ((wlsn-wlsn0)*diffinv)*unitr(ir,2) + rpotls(ic,2)
!       WRITE(9,'(4(1pg13.5))') potp-potp0,potn-potn0,
!     &      bmasp-bmasp0,bmasn-bmasn0
c       
       call alfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn*fc,tlp,tln,sp,sn,
     &             ssp,ssn,frp,frn,
     &             potp,potn,bmasp,bmasn,bmaslp,bmasln,wlsp,wlsn,energd)
       IF(tn<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = ptau(ic)/(tn*reldiff) 
       END IF
       rpot(ir) = (potp-potp0)*diffinv + rpot(ir)
       rpot(ic) = (potn-potn0)*diffinv + rpot(ic)
       rbmass(ir) = (bmasp-bmasp0)*diffinv + rbmass(ir)
       rbmass(ic) = (bmasn-bmasn0)*diffinv + rbmass(ic)
       rpotls(ir,1) = ((wlsp-wlsp0)*diffinv)*unitr(ir,1) + rpotls(ir,1)
       rpotls(ir,2) = ((wlsp-wlsp0)*diffinv)*unitr(ir,2) + rpotls(ir,2)
       rpotls(ic,1) = ((wlsn-wlsn0)*diffinv)*unitr(ir,1) + rpotls(ic,1)
       rpotls(ic,2) = ((wlsn-wlsn0)*diffinv)*unitr(ir,2) + rpotls(ic,2)
!       WRITE(9,'(4(1pg13.5))') potp-potp0,potn-potn0,
!     &      bmasp-bmasp0,bmasn-bmasn0
c       
       call alfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,sp*fc,sn,
     &             ssp,ssn,frp,frn,
     &             potp,potn,bmasp,bmasn,bmaslp,bmasln,wlsp,wlsn,energd)
       IF(ABS(sp)<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = 1.0D0/(sp*reldiff)
       END IF
cassumezero       rpot(ir) = (potp-potp0)*diffinv*ptau(ir) + rpot(ir)
cassumezero       rpot(ic) = (potn-potn0)*diffinv*ptau(ir) + rpot(ic)
cassumezero       rbmass(ir) = (bmasp-bmasp0)*diffinv*ptau(ir) + rbmass(ir)
cassumezero       rbmass(ic) = (bmasn-bmasn0)*diffinv*ptau(ir) + rbmass(ic)
       rpotls(ir,1) = (wlsp-wlsp0)*diffinv*psj(ir,1) + rpotls(ir,1)
       rpotls(ir,2) = (wlsp-wlsp0)*diffinv*psj(ir,2) + rpotls(ir,2)
       rpotls(ic,1) = (wlsn-wlsn0)*diffinv*psj(ir,1) + rpotls(ic,1)
       rpotls(ic,2) = (wlsn-wlsn0)*diffinv*psj(ir,2) + rpotls(ic,2)
c       
       call alfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,sp,sn*fc,
     &             ssp,ssn,frp,frn,
     &             potp,potn,bmasp,bmasn,bmaslp,bmasln,wlsp,wlsn,energd)
       IF(ABS(sn)<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = 1.0D0/(sn*reldiff)
       END IF
cassumezero       rpot(ir) = (potp-potp0)*diffinv*ptau(ic) + rpot(ir)
cassumezero       rpot(ic) = (potn-potn0)*diffinv*ptau(ic) + rpot(ic)
cassumezero       rbmass(ir) = (bmasp-bmasp0)*diffinv*ptau(ic) + rbmass(ir)
cassumezero       rbmass(ic) = (bmasn-bmasn0)*diffinv*ptau(ic) + rbmass(ic)
       rpotls(ir,1) = (wlsp-wlsp0)*diffinv*psj(ic,1) + rpotls(ir,1)
       rpotls(ir,2) = (wlsp-wlsp0)*diffinv*psj(ic,2) + rpotls(ir,2)
       rpotls(ic,1) = (wlsn-wlsn0)*diffinv*psj(ic,1) + rpotls(ic,1)
       rpotls(ic,2) = (wlsn-wlsn0)*diffinv*psj(ic,2) + rpotls(ic,2)
c
       call alfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,sp,sn,
     &             ssp*fc,ssn,frp,frn,
     &             potp,potn,bmasp,bmasn,bmaslp,bmasln,wlsp,wlsn,energd)
       IF(ABS(ssp)<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = psjs(ir)/(ssp*reldiff)
       END IF
       rpot(ir) = (potp-potp0)*diffinv + rpot(ir)
       rpot(ic) = (potn-potn0)*diffinv + rpot(ic)
       rbmass(ir) = (bmasp-bmasp0)*diffinv + rbmass(ir)
       rbmass(ic) = (bmasn-bmasn0)*diffinv + rbmass(ic)
       rpotls(ir,1) = ((wlsp-wlsp0)*diffinv)*unitr(ir,1) + rpotls(ir,1)
       rpotls(ir,2) = ((wlsp-wlsp0)*diffinv)*unitr(ir,2) + rpotls(ir,2)
       rpotls(ic,1) = ((wlsn-wlsn0)*diffinv)*unitr(ir,1) + rpotls(ic,1)
       rpotls(ic,2) = ((wlsn-wlsn0)*diffinv)*unitr(ir,2) + rpotls(ic,2)
!       WRITE(9,'(4(1pg13.5))') potp-potp0,potn-potn0,
!     &      bmasp-bmasp0,bmasn-bmasn0
c
       call alfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,sp,sn,
     &             ssp,ssn*fc,frp,frn,
     &             potp,potn,bmasp,bmasn,bmaslp,bmasln,wlsp,wlsn,energd)
       IF(ABS(ssn)<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = psjs(ic)/(ssn*reldiff)
       END IF
       rpot(ir) = (potp-potp0)*diffinv + rpot(ir)
       rpot(ic) = (potn-potn0)*diffinv + rpot(ic)
       rbmass(ir) = (bmasp-bmasp0)*diffinv + rbmass(ir)
       rbmass(ic) = (bmasn-bmasn0)*diffinv + rbmass(ic)
       rpotls(ir,1) = ((wlsp-wlsp0)*diffinv)*unitr(ir,1) + rpotls(ir,1)
       rpotls(ir,2) = ((wlsp-wlsp0)*diffinv)*unitr(ir,2) + rpotls(ir,2)
       rpotls(ic,1) = ((wlsn-wlsn0)*diffinv)*unitr(ir,1) + rpotls(ic,1)
       rpotls(ic,2) = ((wlsn-wlsn0)*diffinv)*unitr(ir,2) + rpotls(ic,2)
!       WRITE(9,'(4(1pg13.5))') potp-potp0,potn-potn0,
!     &      bmasp-bmasp0,bmasn-bmasn0
c
c       write(6,'(a,i4,8(1pg12.4))')
c     &   ' ir,rpot,rbmass=',ir,rpot(ir),rpot(ic),rbmass(ir),rbmass(ic),
c     &    rpotls(ir,1),rpotls(ir,2),rpotls(ic,1),rpotls(ic,2)
c
      enddo
c
c     contributions from gradient-in-gradient term
c
      IF(ioptsk==13) THEN
        ALLOCATE(potadd(ngrid2,0:0))
        CALL fyfunc_grad_linres2(potadd(1,0))
        rpot(1:ngrid2) = rpot(1:ngrid2) + potadd(1:ngrid2,0)*rho0inv
        DEALLOCATE(potadd)
      END IF
c
c
      parrho = (1-mod(lrpa,2)*2)
      call extrap(rpot(1),parrho)
      call extrap(rpot(ngrid+1),parrho)
      call grad(rbmass,drbmas,2)
      call extrap(rbmass(1),parrho)
      call extrap(rbmass(ngrid+1),parrho)
c
      return
      end
c-----qpotnew------------------------------------------------------
c
      subroutine qpotnew()
c
      include 'params.inc'
      parameter (reldiff=1.0D-7)    ! relative step for derivatives 
c
c     The response potentials 'rkin' and 'rspin' 
c     computed here by numerical derivatives of the mean.field 
c     given in 'alfunc'.
c     The time-odd parts are deduced from the time-even parts
c     assuming that the energy functional is LINEAR in the 
c     kinetic density 'tau' and spin-orbit density 'sjs'.
c
c     The restoration of the c.m. mode is not yet perfect.
c     That is partially due to (lacking) precision of numerical
c     derivatives. One may have to improve on that.
c
      include 'wavfun.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'potinv.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      REAL(8) :: fqcur(KGR2,2)     ! workspace for spin potential
c
      data xminlw/1.0e-20/
c
c--------------------------------------------------------------------
c
c     preparatory steps
c
      parity = 1-2*mod(lrpa,2)
      fc     =1.0D0+reldiff           ! factor for variation
      fm     =1.0D0-reldiff           ! factor for variation
      afac   = alrpa+alrpa*alrpa
      call rotsig(sigma(1),rsigma(1,1),1)
      call rotsig(sigma(1),rsigma(1,2),2)
      call rotsig(sigma(ngrid+1),rsigma(ngrid+1,1),1)
      call rotsig(sigma(ngrid+1),rsigma(ngrid+1,2),2)
c
c     loop over radial grid points
c
      do ir=1,ngrid
c
c      base point
c
       ic   = ir+ngrid
       rinv   = rinvs(ir)
       rp   = rho(ir)
       rn   = rho(ic)
       rs   = rp+rn
       rc   = rp              !    not always true
       d1p  = d1r(ir)
       d1n  = d1r(ic)
       d2p  = d2r(ir)
       d2n  = d2r(ic)
       tp   = tau(ir)
       tn   = tau(ic)
       ts   = tp+tn
       sp   = sj(ir)
       sn   = sj(ic)
       ssp  = sjs(ir)
       ssn  = sjs(ic)
       frp  = 0.0D0
       frn  = 0.0D0
       tlp  = 0D0
       tln  = 0D0
       call alfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,sp,sn,ssp,ssn,
     &             frp,frn,potp0,potn0,bmasp0,bmasn0,bmaslp0,bmasln0,
     &             wlsp0,wlsn0,energ0)
c
c      variation of tp+tn
c
       dt = reldiff*ts*0.5D0
       call alfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp+dt,tn+dt,tlp,tln,
     &             sp,sn,ssp,ssn,
     &             frp,frn,potp,potn,bmasp,bmasn,bmaslp,bmasln,
     &             wlsp,wlsn,energd)
       IF(rs<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = 2.0D0*(energd-energ0)/((dt+dt)*rs)
       END IF
       add = diffinv*(qcur(ir,1)+qcur(ic,1))
       rkin(ir,1) = add + rkin(ir,1)
       rkin(ic,1) = add + rkin(ic,1)
       add = diffinv*(qcur(ir,2)+qcur(ic,2))
       rkin(ir,2) = add + rkin(ir,2)
       rkin(ic,2) = add + rkin(ic,2)
c
c      variation of tp-tn
c
       call alfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp+dt,tn-dt,tlp,tln,
     &             sp,sn,ssp,ssn,
     &             frp,frn,potp,potn,bmasp,bmasn,bmaslp,bmasln,
     &             wlsp,wlsn,energd)
       if(rp.ne.rn) then
         diffinv = 2.0D0*(energd-energ0)/((dt+dt)*(rp-rn))
       else
         diffinv = 0.0D0
       endif
       add = diffinv*(qcur(ir,1)-qcur(ic,1))
       rkin(ir,1) = add + rkin(ir,1)
       rkin(ic,1) = -add + rkin(ic,1)
       add = diffinv*(qcur(ir,2)-qcur(ic,2))
       rkin(ir,2) = add + rkin(ir,2)
       rkin(ic,2) = -add + rkin(ic,2)
c
c      variation of ssp+ssn
c
       ds = reldiff*(ssp+ssn)*0.5D0
       call alfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,
     &             sp,sn,ssp+ds,ssn+ds,
     &             frp,frn,potp,potn,bmasp,bmasn,bmaslp,bmasln,
     &             wlsp,wlsn,energd)
       IF(rs<1D-100)  THEN
         diffinv = 0D0
       ELSE
         diffinv = (energd-energ0)/((ds+ds)*rs)
       END IF
       add = diffinv*(qcur(ir,1)+qcur(ic,1))
       fqcur(ir,1) = add
       fqcur(ic,1) = add
       add = diffinv*(qcur(ir,2)+qcur(ic,2))
       fqcur(ir,2) = add
       fqcur(ic,2) = add
c       add = diffinv*(rqcur(ir,1)+rqcur(ic,1))
c       rspin(ir) = add + rspin(ir)
c       rspin(ic) = add + rspin(ic)
       add = diffinv*(rsigma(ir,1)+rsigma(ic,1))
       rkin(ir,1) = add + rkin(ir,1)
       rkin(ic,1) = add + rkin(ic,1)
       add = diffinv*(rsigma(ir,2)+rsigma(ic,2))
       rkin(ir,2) = add + rkin(ir,2)
       rkin(ic,2) = add + rkin(ic,2)
c
c      variation of ssp-ssn
c
       call alfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,
     &             sp,sn,ssp+ds,ssn-ds,
     &             frp,frn,potp,potn,bmasp,bmasn,bmaslp,bmasln,
     &             wlsp,wlsn,energd)
       if(rp.ne.rn) then
         diffinv = (energd-energ0)/((ds+ds)*(rp-rn)) 
       else
         diffinv = 0.0D0
       endif
       add = diffinv*(qcur(ir,1)-qcur(ic,1))
       fqcur(ir,1) = add + fqcur(ir,1)
       fqcur(ic,1) = -add + fqcur(ic,1)
       add = diffinv*(qcur(ir,2)-qcur(ic,2))
       fqcur(ir,2) = add + fqcur(ir,2)
       fqcur(ic,2) = -add + fqcur(ic,2)
c       rspin(ir) = add + rspin(ir)
c       rspin(ic) = -add + rspin(ic)
       addc = add
       add = diffinv*(rsigma(ir,1)-rsigma(ic,1))
       rkin(ir,1) = add + rkin(ir,1)
       rkin(ic,1) = -add + rkin(ic,1)
       add1 = add
       add = diffinv*(rsigma(ir,2)-rsigma(ic,2))
       rkin(ir,2) = add + rkin(ir,2)
       rkin(ic,2) = -add + rkin(ic,2)
c       write(14,'(i5,4g12.4)') ir,addc,add1,add,diffinv
c
c       write(14,'(a,i4,8(1pg12.4))')
c     &   ' ir,rkin,fqcur=',ir,
c     &    rkin(ir,1),rkin(ir,2),rkin(ic,1),rkin(ic,2),
c     &    fqcur(ir,1),fqcur(ir,2),fqcur(ic,1),fqcur(ic,2)
      enddo
c
c     spin potential "afterburn"
c
!      fqcur=2D0*fqcur; rkin=2D0*rkin
      parrho = (1-mod(lrpa,2)*2)
      call extrap(fqcur(1,1),-parrho)
      call extrap(fqcur(ngrid+1,1),-parrho)
      call extrap(fqcur(1,2),-parrho)
      call extrap(fqcur(ngrid+1,2),-parrho)
      call rot(fqcur,rspin,2)
c
      call extrap(rkin(1,1),-parrho)
      call extrap(rkin(ngrid+1,1),-parrho)
      call extrap(rkin(1,2),-parrho)
      call extrap(rkin(ngrid+1,2),-parrho)
      call div(rkin,drkin,2)
c
      return
      end
c
      subroutine alfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,sp,sn,
     &             ssp,ssn,frp,frn,potp0,potn0,
     &             bmasp0,bmasn0,bmaslp0,bmasln0,wlsp0,wlsn0,energd)

      include 'params.inc'
      include 'skypar.inc'

      IF(ioptsk==13) THEN
       call fyfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,sp,sn,ssp,ssn,
     &             frp,frn,potp0,potn0,bmasp0,bmasn0,bmaslp0,bmasln0,
     &             wlsp0,wlsn0,energd)
!       CALL fyfunc_grad()
      ELSE
       call skfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,sp,sn,ssp,ssn,
     &             frp,frn,potp0,potn0,bmasp0,bmasn0,bmaslp0,bmasln0,
     &             wlsp0,wlsn0,energd)
      END IF

      return
      end
    
