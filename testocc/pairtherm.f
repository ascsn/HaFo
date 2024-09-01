c-----pairdn------------ part of the hartree-fock package ----------------------
c                        changed to fit into the axial H.F.code
c
      subroutine pairdn(e,delta,gw,ph,pairgw,partnm,eferm,epair,avdelt,
     &                  avg,efadd,rtmp,nmax,ndim,iter,eps,iab)
c
      include 'params.inc'
c
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
c      pairgw  = space for dynamical cut-off factor
c      partnm  = desired particle number
c      eferm   = fermi-energy, input is used as initial guess;
c                if eferm=0.0 is input then a suitable initial guess
c                is computed from input 'gw'.
c      epair   = pairing energy
c      avdelt  = average gap
c      avg     = averag pairing force
c      efadd   = pairing band above 'eferm', in case of dynamical pair space
c                see routine  CUTOFF
c      rtmp    = temperature (in cooperation with BCS)
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
      dimension e(ndim),gw(ndim),ph(ndim),delta(ndim),pairgw(ndim)
c
c     several numerical parameters, here for nuclei
c
c      data dltini/1.0/          ! initial gaps
c      data stepd /0.4/          ! stepsize for iteration of gap
      data stepl /0.2/          ! stepsize for iteration of Fermi energy
      data delefe/3.0/          ! difference of Fermi energy for derivative
      data ifcut/1/             ! switch for flexibl cutoff
      data slpmin/2.0/          ! minimum slope in dN/dE
c
c      data one,zero,xsmall/1.0,0.0,1.0E-20/
      data one,zero,half,xsmall/1.0,0.0,0.5,1.0E-20/
      data cutfac/1.65/         ! active pairing space , default 1.65
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
        write(6,'(a,10(/1x,10f8.4))') '    ph:',(ph(i),i=1,nmax)
      endif
c
c                               preset pairing weights
c
      if(eferm.ne.zero .and. ifcut.eq.1) then
       actnum = partnm + cutfac*partnm**0.6666667 
       delcut = 0.1
       efadd  = 20.0/partnm**0.33333333
       call cutoff(e,pairgw,ph,actnum,nmax,eferm,efadd,delcut,
     &             tstnum,iter,eps,iab)
      else
       do n=1,nmax
        pairgw(n) = 1.0
       enddo
      endif
c
c      prepare fermi-energy 'eferm', if not yet initialised.
c
      if(eferm.eq.zero) then
       partsm = zero
       do it=1,nmax
        itsum  = it
        partsm = 2.0*ph(it) + partsm
        if(partsm-partnm.ge.zero) goto 110
       enddo
       stop ' not enough states too reach particle number'
  110  continue
       eferm  = 0.5*(e(itsum)+e(itsum+1))
       if(iab.ge.1) write(06,'(1x,a,2g11.4,a,g11.4)')
     &   'initialization: partsm,partnm=',partsm,partnm,' eferm=',eferm
       do n=1,nmax
         if(n.le.itsum) then
           gw(n) = one
         else
           gw(n) = zero
         endif
       enddo
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
C  
        If (RTMP.EQ.ZERO) THEN
           DO na=1,nmax
              edif   = e(na)-elam
              delta2 = delta(na)*pairgw(na)
              delta2 = delta2*delta2
              gw(na) = half-half*edif/sqrt(edif*edif+delta2)
              parnum = 2.0*gw(na)*ph(na)    + parnum
              edif2  = e(na)-elam2
              parnm2 = 2.0*ph(na)
     &             *(half-half*edif2/sqrt(edif2*edif2+delta2))
     &             + parnm2
           ENDDO
        ELSE
           DO na=1,nmax
C
              edif   = e(na)-elam
              delta2 = delta(na)*pairgw(na)
              delta2 = delta2*delta2
              EQUASI = sqrt(edif*edif+delta2)
              V2     = half-half*edif/EQUASI
              TEXP   = EXP(EQUASI/RTMP)+1
              gw(na) = 1 / TEXP + (1 - 2 / TEXP) * V2
C
              parnum = 2.0*gw(na)*ph(na)    + parnum
C
              edif2  = e(na)-elam2
              EQUASI = sqrt(edif2*edif2+delta2)
              V2     = half-half*edif2/EQUASI
              TEXP   = EXP(EQUASI/RTMP)+1
              parnm2 = 2.0*ph(na)
     &             *(1 / TEXP  + (1 - 2/ TEXP)  * V2)
     &             + parnm2
           ENDDO
        ENDIF
C
        slope  = (parnm2-parnum)/delefe
        slope  = sign(max(abs(slope),slpmin),slope)
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
      sumuv  = 0.0
      sumduv = 0.0
C
      IF (RTMP.EQ.ZERO) THEN
         DO  na=1,nmax
            deltaa = delta(na)*pairgw(na)
            edif   = e(na)-eferm
            equasi = sqrt(edif*edif+deltaa*deltaa)
            v2     = half-half*edif/equasi
            vol    = ph(na)*sqrt(max(v2-v2*v2,xsmall))
            sumuv  = vol        + sumuv
            sumduv = vol*deltaa + sumduv
            sumepa = ph(na)*half*(deltaa*deltaa)/equasi + sumepa
         ENDDO
      ELSE
         DO  na=1,nmax
            deltaa = delta(na)*pairgw(na)
            edif   = e(na)-eferm
            equasi = sqrt(edif*edif+deltaa*deltaa)
            v2     = half-half*edif/equasi
            TEXP   = EXP(EQUASI/RTMP)+1
            TTAN   = (TEXP-2)/TEXP
            vol    = ph(na)*sqrt(max(v2-v2*v2,xsmall))*TTAN
            sumuv  = vol        + sumuv
            sumduv = vol*deltaa + sumduv
            sumepa = ph(na)*half*(deltaa*deltaa*TTAN*TTAN)/equasi 
     &           + sumepa
         ENDDO
      ENDIF
C
      sumuv  = max(sumuv,xsmall)
      avdelt = sumduv/sumuv
      epair  = -sumduv
      avg    = epair/sumuv**2
      if(iab.ge.-1)
     &   write(6,'(a,i4,a,g13.5,a,2g13.5))')
     &  ' PAIRDN: iteration=',it,' part.nr.=',parnum,
     &  ' eferm,epair=',eferm,epair
      if(iab.gt.0) write(6,'(5x,a,3g12.4)')
     &   ' av.Delta,av.G,pair.E=',avdelt,avg,epair
c
      return
      end

