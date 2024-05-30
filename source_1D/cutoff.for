c-----cutoff------------ part of the hartree-fock package ----------------------
c
      subroutine cutoff(e,gw,ph,anz,nmax,eferm,efadd,delta,
     &                  partnm,iter,eps,iab)
c
c     this routine determines the effective "Fermi-energy" and
c     cut-off distribution for the pairing space.
c
c     the include file 'params.inc' communicates the dimensioning of
c     the fields 'e' and 'gw' with the parameter 'ksl2'. it allows
c     furthermore to make the whole subroutine double precision by
c     transfering a line 'implicit double precision (a-h,o-z)'.
c
      include 'params.inc'
c
c     input via list:
c      e       = array(1..nmax) of single particle energies.
c      gw      = array(1..nmax) of cut-off weights,
c                input is used as initial guess for the iteration.
c      ph      = array(1..nmax) of the multiplicities of the single
c                particle levels, to be provided by the calling routine.
c      anz     = desired effective particle number
c      nmax    = maximum number of single particle states
c                (needs to be smaller than the parameter 'ksl2'
c                 and smaller than 46, for dimensioning reasons).
c      eferm   = true fermi-energy as given by the pairing routine.
c      efadd   = additional "fermi-energy" for the pairing-space weights
c      delta   = relative width of the distribution
c      iter    = number of pairing iterations
c      eps     = end condition for iteration;
c                iteration terminates if error in effective particle number
c                is smaller than 'eps' and if 'eferm' and/or 'delta'
c                are smaller than 'eps'.
c      iab     = print level; significant values are -2, -1, 0 and 1;
c                larger value prints more.
c
c     output via list:
c      gw      = array(1..nmax) of the occupation probabilities
c                finally achieved.
c      efadd   = additional "fermi-energy" for the pairing-space weights
c      partnm  = particle number, from adding 'gw' with degeneracy 'ph'.
c
c     some auxiliary fields are:
c      dparn   = array(1..3) of variation in the particle number as
c                consequence of 'deferm' and 'ddelta'.
c      dgap    = array(1..3) of variation in the gap equation as
c                consequence of 'deferm' and 'ddelta'.
c      gr      = double array of variation of particle number and of
c                gap equation, used for construction of the tangantial
c                step
c
      dimension e(ksl2),gw(ksl2),ph(ksl2)
      dimension dparn(2),gr(1,1)
c
cold      data ph/1.,2.,1.,3.,2.,1.,4.,3.,2.,1.,5.,4.,3.,2.,1.,6.,5.,
cold     &        4.,3.,2.,1.,7.,6.,5.,4.,3.,2.,1.,8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        9.,8.,7.,6.,5.,4.,3.,2.,1.,10./
      data gr,dparn/3*0./
c
c     further numerical variables are:
c      deferm  = step in the fermi-energy, to evaluate
c                numerically the derivatives.
c      ddelta  = step in the gap, to evaluate numerically the
c                derivatives; used only in the "constant force" approach.
c      dampin  = damping factor for the tangential step
c      efmax,dlmax,detmin,defmin,ddlmin  = limitations for the step
c
c     note: these parameters are adjusted for an application with
c           the nuclear shell model. it is assumed that the single
c           particle energies are given in units of [mev] and display
c           a reasonable nuclear spectrum.
c
      data dampin,deferm/0.9,1.0/
      data efmax,detmin,defmin/3.,.1e-3,.1e-3/
c
c      data xmaxlw,x10lw/1.0e-30,1.0e-10/
      data one,zero/1.0,0.0/
c
c-------------------------------------------------------------------------------
c
c      set switch variables
c
      iwarng = 0
c
      if(iab.ge.1)
     & write(6,'(/a,i4,3(a,f9.4)/a,10(/1x,10f8.3))')
     &   ' *** cutoff entered with: iter=',iter,'  delta=',delta,
     &   '  efadd=',efadd,'  deferm=',deferm,
     &   '                        single particle energies:',
     &   (e(i),i=1,nmax)
c
c      prepare fermi-energy 'eferm', if not yet initialised.
c
      if(efadd.eq.zero) then
       isum   = 0
       do it=1,nmax
        itsum  = it
        isum   = int(ph(it)+.000001) +isum
        if((one+one)*isum-anz.ge.0) goto 110
       enddo
       stop ' not enough states too reach particle number'
  110  continue
       efadd  = e(itsum)+deferm-eferm
       if(iab.ge.1) write(06,'(5x,a,i5,a,g11.4)')
     &   'initialistion: isum=',isum,'  efadd=',efadd
c
      endif
c
c
c-------------------------------------------------------------------------------
c
c     iterate Fermi energy and the distribution
c
c-------------------------------------------------------------------------------
c
      do it=1,iter
c
c      particle number at 'eferm' (k=1) and 'eferm' varied (k=2).
c
       do ka=1,2
c                                        set actual fermi-energy 'elam'
        k      = 3-ka
        if(k.eq.2) then
         elam   = eferm+efadd+deferm
        else
         elam   = eferm+efadd
        endif
        deltac  = efadd*delta
c                                         compute actual particle number
        parnum = zero
        do i=1,nmax
         gw(i)  = 1.0/(1.0+exp((e(i)-elam)/deltac))
         parnum = gw(i)*ph(i)    + parnum
        enddo
c                                         store particle number
        dparn(k) = parnum+parnum-anz
       enddo
c
c      gradient step to next 'eferm'   (including step limiters)
c
       gr(1,1) = (dparn(2)-dparn(1))/deferm
       det    = sign(max(abs(gr(1,1)),detmin),gr(1,1))
       efold  = efadd
       efstep = dampin*dparn(1)/det
       efadd  = efadd-sign(min(abs(efstep),efmax),efstep)
       deferm = efold-efadd
       deferm = sign(max(abs(deferm),defmin),deferm)
       if(iab.eq.0)
     &  write(6,'(a,i4,2(a,g11.3))')
     &   ' pairing-weight iteration it=',it,'  :  efadd=',efadd,
     &   '  dparn(1)=',dparn(1)
       if(iab.ge.1)
     &  write(6,'(a,i4,2(a,g11.3)/a,2g11.3,a,g11.3/
     &            a,10f7.4,10(/6x,10f7.4))')
     &   ' iteration at it=',it,' :   efadd=',efadd,
     &   '  det=',det,
     &   ' dparn=',dparn(1),dparn(2),' ;   gr=',gr(1,1),
     &   '   gw=',gw
       if(abs(dparn(1)).lt.eps)   goto 350
      enddo
      iwarng = 1
  350 continue
c
c-------------------------------------------------------------------------------
c
c     finally the particle-number is checked and printed
c
c-------------------------------------------------------------------------------
c
      partnm = zero
      do i=1,nmax
       partnm = gw(i)*ph(i) + partnm
      enddo
c
      if(iab.ge.-1)
     &  write(6,'(a,i4,2(a,g13.5))')
     &   ' pair finished with ',it,' iterations: part.nr.=',partnm
c
c     in case of incomplete convergence print a warning
c
      if(iab.ge.-1 .and. iwarng.eq.1)
     & write(6,'(2(a,g11.3),2(a,i3))')
     &  ' ---> no convergence in cutoff! err(part.nr)=',dparn(1)
c
      return
      end

