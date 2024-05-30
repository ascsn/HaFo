c-----parmat------------ part of the hartree-fock package ----------------------c                        changed to fit into the axial H.F.code
c
      subroutine parmat(e,vpairf,v0,eferm,iso,ndim,iab)
c
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   !
      include params.inc
c
c     the include file 'params.inc' allows
c     to make the whole subroutine double precision by
c     transfering a line 'implicit double precision (a-h,o-z)'.
c     Furthermore it carries the common block /CONSTN/ from which
c     the pseudo-constants ONE=1.0 and ZERO=0.0 are taken here.
c
c
c     * * * * * * * * * *   matrix of pairing force   * * * * * * * * * *
c
c     computes the matrix of the pairing force
c
c     input via list:
c      e       = array(1..nmax) of single particle energies.
c                (only used for testprints)
c      v0      = forefactor of the delta-force for the pairing
c      eferm   = Fermi energy (only used for testprints)
c      iso     = isospin of particle sort (1=proton, 2=neutron);
c                the maximum number of single particle states
c                if deduced from 'iso', it needs to be smaller
c                than the parameter 'ndim')
c      ndim    = dimensioning of the arrays, particularly of 'vpairf',
c                in the calling routine
c      iab     = print level; significant values are -2, -1, 0 and 1;
c                larger value prints more.
c
c     output via list:
c      vpairf  = array(1..nmax,1..nmax) of pairing matrix elements
c
      dimension e(ndim)
      dimension vpairf(ndim,ndim)
c
      include 'grid.inc'
      include 'wavfun.inc'
      include 'constn.inc'
c
c
c-------------------------------------------------------------------------------
c
c     the basic pointers to the wavefunctions
c
      nadd   = (iso-1)*nmaxp*ngrid-ngrid
      if(iso.eq.1) then
       nmax   = nmaxp
      else
       nmax   = nmaxn
      endif
c
c     big loop over states
c
      npl1 = nadd
      do 10 n1=1,nmax
        npl1 = ngrid + npl1
        npl2 = nadd
        do 20 n2=1,n1
          npl2 = ngrid + npl2
          vsum = 0.0
          do 30 ir=2,ngrid
            xwf  = wf(ir+npl1)*wf(ir+npl2)
            vsum = zentrf(ir)*(xwf*xwf) + vsum
   30     continue
          vsum = v0*rstep*vsum/fourpi
          vpairf(n1,n2) = vsum
          vpairf(n2,n1) = vsum
   20   continue
   10 continue
c
c     testprints
c
      if(iab.ge.0) then
        write(6,'(a,i2,a)')
     &  ' the pairing matrix for particle sort ISO=',iso,' :'
        do 50 n1=1,nmax
          write(6,'(10f8.3)')
     &    (vpairf(n2,n1),n2=1,n1)
   50   continue
        write(6,'(a)')
     &  ' the pairing matrix elements versus energy differences:',
     &  '  |e(1)-e(2)|  sum |e(i)-eferm|    vpairf '
        do 60 n1=1,nmax
          do 61 n2=1,n1
            write(6,'(3f10.3)')
     &      abs(e(n1)-e(n2)),abs(e(n1)-eferm)+abs(e(n2)-eferm),
     &      vpairf(n2,n1)
   61     continue
   60   continue
      endif
c
      return
      end
c-----flpair------------ part of the hartree-fock package ----------------------c                        changed to fit into the axial H.F.code
c
      subroutine flpair(e,delta,gw,ph,vpairf,partnm,eferm,epair,
     &                  nmax,ndim,iter,eps,iab)
c
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   !
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
c     determines the fermi-energy 'eferm' and the gaps 'delta(n)'
c     such that the "gap-equation" and the "particle number condition"
c     are fulfilled simultanously. this is achieved by iteration.
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
c      vpairf  = array(1..nmax,1..nmax) of pairing matrix elements
c      partnm  = desired particle number
c      eferm   = fermi-energy, input is used as initial guess;
c                if eferm=0.0 is input then a suitable initial guess
c                is computed from input 'gw'.
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
      dimension vpairf(ndim,ndim)
c
c     several numerical parameters, here for nuclei
c
      data dltini/1.0/          ! initial gaps
      data stepd /0.4/          ! stepsize for iteration of gap
      data stepl /0.4/          ! stepsize for iteration of Fermi energy
      data delefe/5.0/          ! difference of Fermi energy for derivative
c
      data one,zero,xsmall/1.0,0.0,1.0E-20/
c
c-------------------------------------------------------------------------------
c
      if(iab.ge.1)
     & write(6,'(/a,i4,1(a,f9.4)/a,10(/1x,10f8.3))')
     &   ' *** full V pairing entered with: iter=',iter,
     &   '  eferm=',eferm,
     &   '                        single particle energies:',
     &   (e(i),i=1,nmax)
c
c      prepare fermi-energy 'eferm' and gaps 'delta', if not yet initialised.
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
         delta(n) = dltini
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
        sumder = 0.0
        sumepa = 0.0
        do 210 na=1,nmax
          sumgap = 0.0
          do 220 nb=1,nmax
            equasi = e(nb)-eferm
            sumgap = vpairf(na,nb)*ph(nb)
     &              *delta(nb)/sqrt(equasi*equasi+delta(nb)*delta(nb))
     &              + sumgap
  220     continue
          error  = delta(na)-sumgap*0.5
          sumder = gw(na)*ph(na)*error*error + sumder
          gw(na) = delta(na) - stepd*error
          edif   = e(na)-eferm
          sumepa = ph(na)*0.25*sumgap
     &            *delta(na)/sqrt(edif*edif+delta(na)*delta(na))
     &            + sumepa
  210   continue
        errgap = sqrt(2.0*sumder/partnm)
        epair  = -sumepa
c
        sumuv  = 0.0
        sumpr1 = 0.0
        sumpr2 = 0.0
        do 230 na=1,nmax
          delta(na) = gw(na)
          deltaa = delta(na)
          edif   = e(na)-eferm
          v2doub = one-edif/sqrt(edif*edif+deltaa*deltaa)
          gw(na) = v2doub*0.5
          sumpr1 = ph(na)*v2doub + sumpr1
          edif   = e(na)-(eferm+delefe)
          v2doub = one-edif/sqrt(edif*edif+deltaa*deltaa)
          sumpr2 = ph(na)*v2doub + sumpr2
  230   continue
        slope  = (sumpr2-sumpr1)/delefe
        slope  = sign(max(abs(slope),one),slope)
        error  = sumpr1-partnm
        eferm  = eferm - stepl*error/slope
c
        if(iab.ge.0)
     &   write(6,'(a,i5,a,2g12.4)')
     &    ' iter=',it,'  errgap,error_partnm=',errgap,error
        if(iab.ge.1)
     &   write(6,'(a,20(/1x,5g12.4))')
     &    ' occupation v**2:',(gw(na),na=1,nmax)
        if(errgap.lt.eps .and. abs(error).lt.eps) goto 299
c
  200 continue
      if(iab.ge.-1)
     & write(6,'(3(a,g11.3),1(a,i3))')
     &  ' ---> no convergence in flpair! err(gap-eq)=',errgap,
     &  '  err(part.nr)=',error,'  nz=',partnm,'  nmzx=',nmax
      return
c
  299 continue
c
      if(iab.ge.-1) then
        write(6,'(a,i4,a,3g13.5)')
     &  ' flpair: iteration=',it,
     &  '  gap-eq.,part.nr.,epair=',errgap,sumpr1,epair
        sumuv  = 0.0
        sumduv = 0.0
        do 330 na=1,nmax
          deltaa = delta(na)
          edif   = e(na)-eferm
          v2     = 0.5-0.5*edif/sqrt(edif*edif+deltaa*deltaa)
          vol    = ph(na)*sqrt(v2-v2*v2)
          sumuv  = vol        + sumuv
          sumduv = vol*deltaa + sumduv
  330   continue
        sumuv  = max(sumuv,xsmall)
        avdelt = sumduv/sumuv
        avg    = epair/sumuv**2
        write(6,'(20x,a,2g12.4)')
     &   ' average delta, average G=',avdelt,avg
        write(6,'(20x,a,2g12.4)')
     &   ' pairing energy,-sumduv=',epair,sumduv
      endif
c
      return
      end

