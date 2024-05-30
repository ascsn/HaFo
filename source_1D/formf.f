c-----fourf ------------ part of the hartree-fock package ----------------------
c
      subroutine fourf(q1,
     &                 rhop,rhon,sjsp,sjsn,
     &                 rstep,ngrid,ifrso,ifform,cmwid,
     &                 foc0,fop,fon,fo,foc)
c
      implicit double precision (a-h,o-z)
c
c     Input variables are:
c      q1    : radial momentum at which the fourier-transformation
c              is to be evaluated
c      rhop    = array of radial proton density (rhop(1) --> at r=0)
c      rhon    = array of radial neutron density
c      sjsp    = Nabla of spin-orbit current for protons (used for ifrso=1)
c      sjsn    = Nabla of spin-orbit current for neutrons (used for ifrso=1)
c      rstep   = spacing of the (equidistant) radial grid
c      ngrid   = number of radial grid points
c      ifrso   = 1 --> chrarge density contains contribtion from
c                      spin-orbit current
c      ifform  = 1 --> nucleon formfactors are folded in
c                0 --> mere proton density is analyzed
c      cmwid : parameter for center of mass correction
c     output paramaters:
c      foc0 :   charge formfactor
c      fop  :   proton "
c      fon  :   neutron    "
c      fo   :   total      "
c      foc  :   charge     "   including spin-orbit current
c
      dimension rhop(ngrid),rhon(ngrid)         ! input information
      dimension sjsp(ngrid),sjsn(ngrid)         ! input information
c
c      data  zero,one,pi/0.0,1.0,3.1415926/
      data  zero,half,one,pi/0.0,0.5,1.0,3.1415926/
c
c      the parameters in data are for the sachs-formfactors:
c       sa..,sm.. are for the isoscalar electric formfactor
c       va..,vm.. are for the isovector electric formfactor
c       gp..,gm.. are for the magnetic formfactors of proton and neutron
c
c       isoscalar and isovector formfactors are recoupled to
c       the electric formfactors of proton and neutron.
c
c     the sachs-formfactors are then recoupled to the nucleon
c     vector formfactor (form1..) and the nucleon tensor
c     formfactor (form2..).
c
c
c     parameters for scalar and vector electric formfactor
c     (walther, priv.comm., sept. 86)
c
      data    sa1,    sa2,    sa3,   sa4,   sm1,  sm2,  sm3,  sm4
     *    /2.2907,-0.6777,-0.7923,0.1793, 15.75,26.68,41.01,134.2/
      data    va1,    va2,    va3,   va4,   vm1,  vm2,  vm3,  vm4
     *    /0.3681, 1.2263,-0.6316,0.0372,  5.00,15.02,44.08,154.2/
c
c     magnetic formfactor , simon et al n.p. a333 (1980) 381
c
      data    gp1,    gp2,    gp3,   gp4,   gm1,  gm2,  gm3,  gm4
     *    / 0.694,  0.719,-0.418 ,0.005 ,  8.50,15.02,44.08,355.4/
c
c     old version of the
c     parameters for the magnetic formfactor (arenhoevel 1980)
c
c      data gp1,gp2,gp3,gp4,gm1,gm2,gm3,gm4/.794,.594,-.393,.005,
c     *                      8.73,15.02,44.08,355.40/
c
c     the anomanlous magnetic moments and the darwin factor 'd4m'
c
      data apro0,apro11,apro12,bpro11
     *    / 0.74, 0.0  ,  0.0 ,  0.0  /
      data aneu0,aneu11,aneu12,bneu11
     *    / 0.74, 0.0  ,  0.0 ,  0.0  /
      data xmgne,xmgpr,d4m/.04197,-.06131,.010989/
c
c----------------------------------------------------------------
c
      frpirs = 4*pi*rstep
c
c     fourier-bessel transformation for proton- and neutron-density
c     given on field 'rho' on common /wavfun/.
c
      dc     = cos(q1*rstep)
      ds     = sin(q1*rstep)
      r      = zero
      sump   = zero
      sumn   = zero
      c      = one
      s      = zero
       do 10 i=2,ngrid
       r      = rstep+r
       x      = c*dc-s*ds
       s      = s*dc+c*ds
       c      = x
       sump   = r*s*rhop(i)+sump
       sumn   = r*s*rhon(i)+sumn
   10  continue
      sump   = frpirs*sump/q1
      sumn   = frpirs*sumn/q1
c
c     fourier-bessel transformation for spin-orbit density of
c     protons and neutrons given on field 'sjs'.
c
      if(ifrso.eq.1) then
       sumsp  = zero
       sumsn  = zero
       r      = zero
       c      = one
       s      = zero
        do 20 i=2,ngrid
        r      = rstep+r
        x      = c*dc-s*ds
        s      = s*dc+c*ds
        c      = x
        sumsp  = r*s*sjsp(i)+sumsp
        sumsn  = r*s*sjsn(i)+sumsn
   20   continue
       sumsp=frpirs*sumsp/q1
       sumsn=frpirs*sumsn/q1
      else
       sumsp  = zero
       sumsn  = zero
      endif
c
c     three cases for the nucleon folding:
c      ifform = 0  --> no folding is done
c      ifform = 1  --> folding with free nucleon formfactor parametrised
c                      by sum of lorentzians, see data on sa.., va..
c                      and gp..
c      ifform =-1  --> density dependent folding for charge density
c
      if(ifform.eq.0) then
       schsmp = 0.0
c       schcmn = 0.0
       sumc   = sump
      elseif(ifform.gt.0) then
c
c      isoscalar and isovector electrical formactors on 'schses'
c      and 'schsev', recoupled then to proton and neutron electric
c      formfactors  'schsep' and 'schsen'. proton and neutron
c      magnetic formfactor on 'schsmp' and 'schsmn'.
c
       qa     = q1*q1
       schses = sa1/(1.+qa/sm1)+sa2/(1.+qa/sm2)+sa3/(1.+qa/sm3)
     &         +sa4/(1.+qa/sm4)
       schsev = va1/(1.+qa/vm1)+va2/(1.+qa/vm2)+va3/(1.+qa/vm3)
     &         +va4/(1.+qa/vm4)
       schsep = half*(schses+schsev)
       schsen = half*(schses-schsev)
       schsm  = ( gp1/(one+qa/gm1)+gp2/(one+qa/gm2)+gp3/(one+qa/gm3)
     &           +gp4/(one+qa/gm4) )
       schsmp = schsm*xmgpr - schsep*d4m
       schsmn = schsm*xmgne - schsen*d4m
c
       sumc   = sump*schsep+sumn*schsen
c
      else
c
c      charge fromfactor with density-dependent intrinsic nucleon radii,
c      density dependence is given by the parameters 'apro..','bpro..',
c      'aneu..' and 'bneu..' transferred via common /rmsnuc/.
c
c       iradd = ngrid
c       if(ifform.eq.-2) iradd=0
       fq2    = q1*q1
       dc     = cos(q1*rstep)
       ds     = sin(q1*rstep)
       r      = zero
       sumc   = zero
       c      = one
       s      = zero
       do 30 i=2,ngrid
       r      = rstep+r
       x      = c*dc-s*ds
       s      = s*dc+c*ds
       c      = x
c                                   proton-radius as fkt. of density:
       rhotot = rhop(i)+rhon(i)
       rmspi2 = ( apro0+apro11*rhotot
     &           +apro12*rhotot*rhotot
     &           +bpro11*(rhop(i)-rhon(i)) )/6.0
c                                   neutron-radius as fkt. of density:
       rmsni2 = aneu0+aneu11*rhotot
     &         +aneu12*rhotot*rhotot
     &         +bneu11*(rhop(i)-rhon(i))
c                                   accumlate folding:
       fintp  = exp(-rmspi2*fq2)
       fintn  = 12./(12.+(0.64+rmsni2)*fq2)-12./(12.+(0.64-rmsni2)*fq2)
       sumc   = r*s*(rhop(i)*fintp + rhon(i)*fintn) + sumc
   30  continue
       sumc   = frpirs*sumc/q1
c
      endif
c
c      factor for darwin correction 'darfac'
c      and the factor for the c.m. correction 'exphf'.
c
       qa     = q1*q1
       if(ifrso.eq.1) then
        darfac = one/sqrt(one+qa*d4m)
       else
        darfac = one
       endif
       exphf  = exp(cmwid*qa)
c
c     finally assembling the formfactors
c
      foc0   = exphf*sumc
      fop    = exphf*sump
      fon    = exphf*sumn
      fo     = exphf*(sump+sumn)
      foc    = exphf*(sumc*darfac + (schsmp*sumsp+schsmn*sumsn))
c
      return
      end
