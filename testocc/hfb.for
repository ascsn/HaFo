c-----hfbpot------------------------------------------------------
c
      subroutine hfbpot
c
      include 'params.inc'
c
c     computes the pair density
c
c       sum_n  u_n*v_n*|phi_n(r)|^2  ,
c
c     subsequently the pair potential on 'pairpt' from /wavfun/,
c     and finally the dynamic gaps on 'delatf' from /shells/
c     as expectation value of 'pairpt'.
c
c     'ipair' is used here to indicate the switch to surface pairing
c             5,6 --> volume pairing
c             7,8 --> surface pairing
c
c
      include 'wavfun.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c--------------------------------------------------------------------
c
c
c     zero field before accumulation
c
      do i=1,kgr2
        paircr(i) = zero
      enddo
c
c     accumulate:
c
c     iso-loop, initialise "shift constants" 'im,ia,ie' for the
c     densities, 'ngw' for the weights, and 'nadd' for the
c     wave functions.
c
      do 20 iso=1,2
      ie     = ngrid*iso
      im     = ie-ngrid+1
      ia     = im+1
      if(iso.eq.1) then
       nmax   = nmaxp
      else
       nmax   = nmaxn
      endif
      ngw    = (iso-1)*nmaxp
      nadd   = ngw*ngrid-ie
c
c     loop over states for given isospin.
c      'nshl' counts the states in harmonic oscillator ordering.
c      'js' is the according composed ang.mom.,
c      'll' is the orbital ang.mom.,
c      'wl' is the spin-orbit splitting,
c      'om' is the actual phase-space factor
c
      do 30 nshl=1,nmax
      js     = jsp(nshl)
      ll     = js/2
c      xlp    = float(ll*ll+ll+1)
c      wl     = wls(js)
      v2     = weight(nshl+ngw)
      uv     = sqrt(max(v2-v2*v2,xsmall))
      om     = zwj(js)*uv*0.5      ! factor 1/2 because "state>0"
      nadd   = nadd+ngrid
      if(om.lt.xsmall) goto 30
c                                  the case r=0.0 needs special treatment
      if(ll.lt.1) then
       x      = wfs(im+nadd)
       paircr(im) = (x*x)*om + paircr(im)
      endif
c                                  now all r>0.0
      r      = zero
      do 40  ir=ia,ie
      r      = r+rstep
      wfa    = wf(nadd+ir)/r
      paircr(ir) = om*(wfa*wfa) + paircr(ir)
   40 continue
   30 continue
   20 continue
c
c     now multiply with strength to obtain local pair-potential
c
      if(ipair.le.6) then
        v0prot = -pairfp
        v0neut = -pairfn
        do i=1,ngrid
          pairpt(i)  = v0prot * paircr(i)
          in = i+ngrid
          pairpt(in) = v0neut * paircr(in)
        enddo
      else
        v0prot = -pairfp
        v0neut = -pairfn
        do i=1,ngrid
          in = i+ngrid
          rhofac = one-(rho(i)+rho(in))/rho0pr
          pairpt(i)  = v0prot * paircr(i)  * rhofac
          pairpt(in) = v0neut * paircr(in) * rhofac
        enddo
      endif
c
c     finally compute the actual gaps as s.p. expectation values with
c     the pair potential
c
      nadd   = -ngrid
      do 60 istate=1,nmaxto
      if(istate.le.nmaxp) then
        isoadd = 0
        iso    = 1
      else
        isoadd = ngrid
        iso    = 2
      endif
      nadd   = ngrid+nadd
      sum    = zero
      pptsum = zero
      do ir=1,ngrid
        wfa    = wf(nadd+ir)
        pptsum = max(abs(pairpt(isoadd+ir)),pptsum)
        sum    = (wfa*wfa)*pairpt(isoadd+ir) + sum
      enddo
      pptmax(iso) = pptsum
      deltaf(istate) = -rstep*sum
   60 continue
c
      if(iwr.gt.-1) then
        write(6,'(a)') ' dynamics gaps for protons, directly:'
        write(6,'(5(1pg12.4))') (deltaf(i),i=1,nmaxp)
        write(6,'(a)') ' dynamics gaps for neutrons, directly:'
        write(6,'(5(1pg12.4))') (deltaf(i+nmaxp),i=1,nmaxn)
      endif
        write(6,'(a)') ' pair current:'
        write(6,'(5(1pg12.4))') (paircr(i+nmaxp),i=1,nmaxn)
        write(6,'(a)') ' pair potential:'
        write(6,'(5(1pg12.4))') (pairpt(i+nmaxp),i=1,nmaxn)
c
      return
      end
c-----pairdn------------ part of the hartree-fock package ----------------------
c                        changed to fit into the axial H.F.code
c
      subroutine pairdn(e,delta,gw,ph,partnm,eferm,epair,avdelt,avg,
     &                  nmax,ndim,iter,eps,iab)
c
      include params.inc
c
c     the include file 'params.inc' allows
c     to make the whole subroutine double precision by
c     transfering a line 'implicit double precision (a-h,o-z)'.
c     Furthermore it carries the common block /CONSTN/ from which
c     the pseudo-constants ONE=1.0 and ZERO=0.0 are taken here.
c
c
c     * * * * * * * * * *   pairing with full interaction   * * * * * * * * * *
c
c     determines the fermi-energy 'eferm'
c     such that the "particle number condition"
c     is fulfilled. this is achieved by iteration.
c
c     input via list:
c      e       = array(1..nmax) of single particle energies.
c      delta   = array(1..nmax) of gap parameters.
c      gw      = array(1..nmax) of occupation probabilities,
c                input is used as initila guess for the iteration.
c      ph      = array(1..nmax) of the multiplicities of the single
c                particle levels.
c                This multiplicities excludes (!) the spin degeneracy,
c                e.g. a 1s1/2 state has ph = 1 .
c      partnm  = desired particle number
c      eferm   = fermi-energy, input is used as initial guess;
c                if eferm=0.0 is input then a suitable initial guess
c                is computed from input 'gw'.
c      epair   = pairing energy
c      avdelt  = average gap
c      avg     = averag pairing force
c      nmax    = maximum number of single particle states
c                (needs to be smaller than the parameter 'ndim')
c      ndim    = dimensioning of the arrays, particularly of 'vpairf',
c                in the calling routine
c      iter    = number of pairing iterations,
c                abs(iter) is the number of iterations
c      eps     = end condition for pairing iteration;
c                iteration terminates if error in particle number
c                is smaller than 'eps' and if 'eferm' and/or 'delta'
c                are smaller than 'eps'.
c      iab     = print level; significant values are -2, -1, 0 and 1;
c                larger value prints more.
c
c     output via list:
c      gw      = array(1..nmax) of the occupation probabilities
c                finally achieved.
c      eferm   = fermi-energy.
c
      dimension e(ndim),gw(ndim),ph(ndim),delta(ndim)
c
c     several numerical parameters, here for nuclei
c
c      data dltini/1.0/          ! initial gaps
c      data stepd /0.4/          ! stepsize for iteration of gap
      data stepl /0.4/          ! stepsize for iteration of Fermi energy
      data delefe/5.0/          ! difference of Fermi energy for derivative
c
      data one,zero,half,xsmall/1.0,0.0,0.5,1.0E-20/
c
c-------------------------------------------------------------------------------
c
      if(iab.gt.-1) then
        write(6,'(/a,i4,1(a,f9.4)/a,10(/1x,10f8.3))')
     &   ' *** pairing with dynamic gap entered with: iter=',iter,
     &   '  eferm=',eferm,
     &   '  single particle energies:',
     &   (e(i),i=1,nmax)
        write(6,'(a,10(/1x,10f8.4))') '  delta:',(delta(i),i=1,nmax)
c        write(6,'(a,10(/1x,10f8.4))') '    ph:',(ph(i),i=1,nmax)
      endif
c
c      prepare fermi-energy 'eferm', if not yet initialised.
c
      if(eferm.eq.zero) then
       partsm = zero
       do 120 it=1,nmax
       itsum  = it
       partsm = 2.0*ph(it) + partsm
       if(partsm-partnm.ge.zero) goto 110
  120  continue
       stop ' not enough states too reach particle number'
  110  continue
       eferm  = 0.5*(e(itsum)+e(itsum+1))
       if(iab.ge.1) write(06,'(1x,a,2g11.4,a,g11.4)')
     &   'initialization: partsm,partnm=',partsm,partnm,' eferm=',eferm
       do 130 n=1,nmax
         if(n.le.itsum) then
           gw(n) = one
         else
           gw(n) = zero
         endif
  130  continue
c
      endif
c
c     the coupled gradient iteration,
c     misusing 'gw' as workspace
c
      do 200 it=1,iter
c
c
        elam   = eferm
        elam2  = eferm+delefe             ! varied to compute derivative
        parnum = zero
        parnm2 = zero
        do 230 na=1,nmax
          edif   = e(na)-elam
          delta2 = delta(na)*delta(na)
          gw(na) = half-half*edif/sqrt(edif*edif+delta2)
          parnum = 2.0*gw(na)*ph(na)    + parnum
          edif2  = e(na)-elam2
          parnm2 = 2.0*ph(na)*(half-half*edif2/sqrt(edif2*edif2+delta2))
     &            + parnm2
  230   continue
        slope  = (parnm2-parnum)/delefe
        slope  = sign(max(abs(slope),one),slope)
        error  = parnum-partnm
        eferm  = eferm - stepl*error/slope
c
        if(iab.ge.0)
     &   write(6,'(a,i5,a,5g12.4)')
     &    ' iter=',it,' eferm,parnum,parnm2,slope,error_partnm=',
     &                  eferm,parnum,parnm2,slope,error
        if(iab.ge.1)
     &   write(6,'(a,20(/1x,5g12.4))')
     &    ' occupation v**2:',(gw(na),na=1,nmax)
        if(abs(error).lt.eps) goto 299
c
  200 continue
      if(iab.ge.-1)
     & write(6,'(2(a,g11.3),1(a,i3))')
     &  ' ---> no convergence in PAIRDN!  err(part.nr)=',error,
     &  '  nz=',partnm,'  nmzx=',nmax
      return
c
  299 continue
c
c     finally compute average gap and pairing energy
c
      if(iab.ge.-1)
     &   write(6,'(a,i4,a,3g13.5)')
     &  ' PAIRDN: iteration=',it,' part.nr.=',parnum
      sumuv  = 0.0
      sumduv = 0.0
      do 330 na=1,nmax
        deltaa = delta(na)
        edif   = e(na)-eferm
        equasi = sqrt(edif*edif+deltaa*deltaa)
        v2     = half-half*edif/equasi
        vol    = ph(na)*sqrt(v2-v2*v2)
        sumuv  = vol        + sumuv
        sumduv = vol*deltaa + sumduv
        sumepa = ph(na)*half*(deltaa*deltaa)/equasi
     &          + sumepa
 330  continue
      sumuv  = max(sumuv,xsmall)
      avdelt = sumduv/sumuv
      epair  = -sumduv
      avg    = epair/sumuv**2
      if(iab.gt.0) write(6,'(5x,a,3g12.4)')
     &   ' av.Delta,av.G,pair.E=',avdelt,avg,epair
c
      return
      end
c-----lagpar------------ part of the hartree-fock package ----------------------
c
      subroutine lagpar
c
      include 'params.inc'
c
c     determines the matrix of Lagrange parameters for HFB iteration
c     and stores result on PAIRMT on common /shells/.
c     the action of the HFB Hamiltonian and the wavefunction is
c     stored on the workspace SAVLR from common /linres/
c
      dimension savpri(5)
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      data xminv/1.0E-05/
c
c-------------------------------------------------------------------------------
c
c     reset matrix of Lagrange parameters
c
      do na=1,nmaxp
        do nb=1,nmaxp
          pairmp(na,nb) = zero
        enddo
      enddo
      do na=1,nmaxn
        do nb=1,nmaxn
          pairmn(na,nb) = zero
        enddo
      enddo
c     iso-loop
c      'isoadd' is a pointer for the potential and mass fields;
c      'isnpad' is a pointer for the wave function pointer 'npla';
c      on 'egx' is accumulated lateron the shell-energy.
c
      do 1000 iso=1,2
      if(iso.eq.1) then
        isoadd = 0
        isnpad = 0
        eferm  = efep
      else
        isoadd = ngrid
        isnpad = kjmx
        eferm  = efen
      endif
c
c     square root of inverse free mass on 'sqb';
c
      call derv15(bmass(isoadd+1),d1b(isoadd+1),one,ngrid,rstep)
c
c  FIRST BIG LOOP: compute HFB-Hamiltonian on wavefunction and
c                  matrix of Lagrange parameters
c
c     loop over composed ang.mom. 'js'
c
      do 40 js=1,jsmax(iso)
      jsnpl = js+isnpad
      nmax  = npla(1,jsnpl)
      if(nmax.gt.0) then
c
c      loop over radial quantum number 'n'.
c       accumulating matrix elements with HFB Hamiltonian (mix of
c       mean field and gap potential)
c       'npl' is wave function pointer,
c       'nplshl' is number of state in harmonic oscillator ordering.
c
       do 70 n=1,nmax
       npl    = npla(n+1,jsnpl)
       nplm   = npl-1
       nplshl = npl/ngrid+1
c
c      action of the mean-field hamiltonian
c
       call hpsi(wf(npl),aux3,wfs(npl),aux2,0,js,iso,zero)
       if(iwr.ge.0) then
         testsp = zero
         do ir=2,ngrid
           testsp   = aux3(ir)*wf(nplm+ir)    + testsp
         enddo
         delcm(nplshl) = testsp*rstep
       endif
c
c      complement to actual HFB Hamiltonian
c
       v2     = weight(nplshl)
       uv     = sqrt(max(v2-v2*v2,xsmall))
       do ir=1,ngrid
         aux3(ir) = v2*aux3(ir)+uv*pairpt(ir+isoadd)*wf(nplm+ir)
         savlr(nplm+ir,1) = aux3(ir)                    ! save action of HFB
       enddo
       if(iwr.ge.0) then
         testsp = zero
         do ir=1,ngrid
           testsp   = aux3(ir)*wf(nplm+ir)    + testsp
         enddo
         drkin(nplshl) = rstep*testsp
       endif
c
c      accumulate to matrix of Lagrange parameters, thereby symmetrising
c
       do 71 nb=1,nmax
       nplb   = npla(nb+1,jsnpl)
       nplmb  = nplb-1
       nplshb = nplb/ngrid+1
       sum    = zero
       do ir=2,ngrid
         sum    = wf(nplmb+ir)*aux3(ir) + sum
       enddo
       sum    = 0.5*rstep*sum
       if(iso.eq.1) then
         pairmp(nplshb,nplshl) = sum + pairmp(nplshb,nplshl)
         pairmp(nplshl,nplshb) = sum + pairmp(nplshl,nplshb)
       else
         nxa = nplshl-nmaxp
         nxb = nplshb-nmaxp
         pairmn(nxb,nxa) = sum + pairmn(nxb,nxa)
         pairmn(nxa,nxb) = sum + pairmn(nxa,nxb)
       endif
   71  continue
c
   70  continue
      endif
   40 continue
c
c  END OF FIRST BIG LOOP
c
c
c
c  SECONG BIG LOOP: subract Lagrange parameters from stored wavefunctions
c
c
c     loop over composed ang.mom. 'js'
c
      do 140 js=1,jsmax(iso)
      jsnpl = js+isnpad
      nmax  = npla(1,jsnpl)
      if(nmax.gt.0) then
c
c      double loop over radial quantum number 'n', 'nb'.
c
       if(iwr.gt.0)  write(6,'(a,2i5)')
     &   ' Lagrange matrix for: iso,js=',iso,js
       do 170 n=1,nmax
       npl    = npla(n+1,jsnpl)
       nplm   = npl-1
       nplshl = npl/ngrid+1
c
       do 171 nb=1,nmax
       nplb   = npla(nb+1,jsnpl)
       nplmb  = nplb-1
       nplshb = nplb/ngrid+1
       if(iso.eq.1) then
         fac    = pairmp(nplshb,nplshl)
       else
         fac    = pairmn(nplshb-nmaxp,nplshl-nmaxp)
       endif
       savpri(nb) = fac
       do ir=2,ngrid
         savlr(nplm+ir,1) = savlr(nplm+ir,1) - fac*wf(nplmb+ir)
       enddo
  171  continue
       if(iwr.gt.0) write(6,'(1x,5(1pg12.4))') (savpri(nb),nb=1,nmax)
       vpact = 1.0/max(sqrt(weight(nplshl)),xminv)
       do ir=2,ngrid
         savlr(nplm+ir,1) = vpact*savlr(nplm+ir,1)
       enddo
c
  170  continue
      endif
  140 continue
c
c  END OF SECOND BIG LOOP
c
c    test prints
c
      if(iwr.ge.0) then
        if(iso.eq.1) then
         npr1 = 1
         npr2 = nmaxp
        else
         npr1 = nmaxp+1
         npr2 = nmaxn+nmaxp
        endif
        write(6,'(a,10(/1x,10f8.3))')
     &        ' spe direct:',(delcm(i),i=npr1,npr2)
        write(6,'(a,10(/1x,10f8.3))') ' spe HFB:',(drkin(i),i=npr1,npr2)
      endif
c
 1000 continue
c
      return
      end

