c-----hpsi  ------------ part of the hartree-fock package ----------------------
c
      subroutine test_hpsi()
c
      include 'params.inc'
c
c     test action of the mean-field hamiltonian term by term
c
c     input:
c       wfin    = input wave function field ,
c       hwfout  = output wave function field ,
c       d1wfin  = first derivative of 'wfin'.
c       d2wfin  = second derivative of 'wfin'.
c       ifd1wf  = decides on the computation of 'd1wfin':
c                 ifd1wf>0 --> the 'd1wfin' is interpreted as
c                              workspace and the first derivative
c                              is computed anew;
c                 ifd1wf=0 --> the 'd1wfin' contains already the
c                              first derivative if 'wfin'.
c       jsa     = combined angular momentum quantum number
c                 of wfin and hwfout,
c       iso     = isospin of the wave functions,
c       subfac  = constant energy subtracted from hamiltonian.
c
c     the fields for the mean-field hamiltonian are transferred
c     via common /wavfun/:
c       pot     = potential
c       potls   = spin-orbit potential
c       bmass   = inverse effective mass
c       d1b     = derivative of inverse effective mass
c
c     grid properties are transferred via common /grid/
c
      dimension wfin(kg1),hwfout(kg1),d1wfin(kg1),d2wfin(kg1)
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'constn.inc'
      include 'result.inc'
      include 'skypar.inc'
      include 'linres.inc'
c
c-------------------------------------------------------------------------------
c
c     the quantum numbers of the wave function
c      la      = orbital angular momentum of state  wfin>
c      parity  = parity for the radial wave function (note that it
c                contains an extra factor 'r' such that the 1s1/2
c                state has parity=-1 etc)
c      isoadd  = stack-pointer for potential fields pot, potls
c                and isoadd
c      xl0     = spin-orbit splitting for state  wfin>
c
      WRITE(6,'(a)') 'testprint partial s.p. energies'
      WRITE(6,'(a)') ' iso  jsa   n   norm    pot    ls     mass    '//
     &      '  kin    kinetic  s.p.energy   s.p.kin'
      rmstot=0D0
      totnum=0D0
      do 1000 iso=1,2
      isoadd = (iso-1)*ngrid
      isnpad = (iso-1)*kjmx
      do 40 jsa=1,jsmax(iso)
      jsnpl = jsa+isnpad
      nmax  = npla(1,jsnpl)
      if(nmax.le.0) goto 40
      la     = jsa/2
      xlp    = (la*la+la)
      parity = (mod(la,2)*2-1)
      isoadd = (iso-1)*ngrid
      xl0    = wls(jsa)
      phasjs = float(jsa-la)
      do 70 n=1,nmax
      npl    = npla(n+1,jsnpl)
      nplm   = npl-1
      nplshl = npl/ngrid+1
      wfin(1:ngrid)=wf(npl:nplm+ngrid)
c
c     the derivatives of the wave function
c
!      if(ifd1wf.gt.0)
      call derv17(wfin(1),d1wfin(1),parity,ngrid,rstep)
      call derv25(wfin(1),d2wfin(1),parity,ngrid,rstep)
c
c     norms and overlaps
c
!      xnorm1=0D0
!      xnorm2=0D0
!      ovlp12=0D0
!      hwfout(1:ngrid)=savlr(npl:npl+ngrid,2)
!      do ir=2,ngrid
!        irs    = ir+isoadd
!        xnorm1 = xnorm1 + wfin(ir)*wfin(ir)
!        xnorm2 = xnorm2 + hwfout(ir)*hwfout(ir)
!        ovlp12 = ovlp12 + hwfout(ir)*wfin(ir)
!      enddo
!      xnorm1=rstep*xnorm1
!      xnorm2=rstep*xnorm2
!      ovlp12=rstep*ovlp12
c
c     the action of the potential
c
      spepot = 0D0
      xrms = 0D0
      do ir=2,ngrid
        irs    = ir+isoadd
        potact = pot(irs)
        spepot = spepot + potact*wfin(ir)*wfin(ir)
        xrms = xrms + (ir-1)**2*wfin(ir)*wfin(ir)
      enddo
      xrms=rstep**3*xrms
      rmstot=rmstot+xrms*weight(nplshl)*phasjs
      xrms=sqrt(xrms)
      totnum=totnum+weight(nplshl)*phasjs
      spepot=rstep*spepot
c
c     the action of the spin-orbit potential
c
      spels = 0D0
      do ir=2,ngrid
        irs    = ir+isoadd
        potact = xl0*potls(irs)
        spels = spels + potact*wfin(ir)*wfin(ir)
      enddo
      spels=rstep*spels
c
c     the action of the effective mass
c
      spems=0D0
      do ir=2,ngrid
        irs    = ir+isoadd
        potact =(bmass(irs)*zentrf(ir)+bmassl(irs))*xlp
     &          +d1b(irs)*rinvs(ir)
        hwfout(ir) = potact*wfin(ir)
     &              -d1b(irs)*d1wfin(ir)-bmass(irs)*d2wfin(ir)
        spems = spems + hwfout(ir)*wfin(ir)
      enddo
      spems=rstep*spems
c
c     the action of the mere kinetic energy
c
      spekin=0D0
      amfcac = amfc(iso)
      hwfout(1) = amfcac*d2wfin(1)
      do ir=2,ngrid
c      r      = (ir-1)*rstep
        potact = amfcac*xlp*zentrf(ir)
        hwfout(ir) = potact*wfin(ir) - amfcac*d2wfin(ir)
        spekin = spekin + hwfout(ir)*wfin(ir)
      enddo
      spekin=rstep*spekin
c
c     printout
c
!      WRITE(6,'(3i5,8(1pg13.5)/15x,8(1pg13.5))') iso,jsa,n,xrms,
      WRITE(6,'(3i5,8(1pg13.5))') iso,jsa,n,xrms,
     &   spepot,spels,spems,spekin,spepot+spels+spems,spenrg(nplshl),
     &   spkin(nplshl)
c
 70   continue
 40   continue
 1000 continue
      WRITE(6,'(a,1pg13.5)') 'test rms_tot=',sqrt(rmstot/totnum)
c
      return
      end
