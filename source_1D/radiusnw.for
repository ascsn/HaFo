c-----radius-----------------------------------------------------------
c
      subroutine radius(rhop,rhon,sjsp,sjsn,
     &                  nprot,nneut,rstep,ngrid,amfpr,amfne,zpe,ifzpe,
     &                  iwr,ifrso,ifform,
     &                  formf,rhochr,
     &                  rmsc,rmsp,rmsn,rms,
     &                  rdmsc,rdmsp,rdmsn,rdms,
     &                  surfc,surfp,surfn,surf)
c
      implicit double precision (a-h,o-z)
c
c     * * * * *    rms-radius, dms-radius and surface-thickness    * * * * *
c     * * * *    for proton-, neutron-, total- and charge-density    * * * *
c     * * * * * * *    given on an equidistant radial grid     * * * * * * *
c
c     the charge-density includes folding with intrinsic electric and
c     magentic formfactors of proton and neutron.
c     the contribution of the anomalous magnetic moment is included
c     only if ifrso=1.
c
c     Input variables are:
c      rhop    = array of radial proton density (rhop(1) --> at r=0)
c      rhon    = array of radial neutron density
c      sjsp    = Nabla of spin-orbit current for protons (used for ifrso=1)
c      sjsn    = Nabla of spin-orbit current for neutrons (used for ifrso=1)
c      nprot   = number of protons
c      nneut   = number of neutrons
c      rstep   = spacing of the (equidistant) radial grid
c      ngrid   = number of radial grid points
c      amfpr   = hbar^2/2m_proton
c      amfne   = hbar^2/2m_neutron
c      zpe     = c.m. corection energy (>0), in case of <P^2>/2mA recipe
c      ifzpe   = 1 -->  <P^2>/2mA recipe, requires input of 'zpe'
c                0 -->  m-m/A recipe
c      iwr     = print level
c                -2 suppresses any printing
c                -1 gives reasonable amount
c                > -1 only for testing
c      ifrso   = 1 --> chrarge density contains contribtion from
c                      spin-orbit current
c      ifform  = 1 --> nucleon formfactors are folded in
c                0 --> mere proton density is analyzed
c
c     Output variables are:
c      formf   = charge formfactor
c      rhochr  = charge density inluding all nucleon formfactors
c      rmsc    = charge r.m.s. radius
c      rmsp    = proton r.m.s. radius
c      rmsn    = neutron r.m.s. radius
c      rms     = total r.m.s. radius
c      rdmsc   = charge diffraction radius (sharp- or box-radius)
c      rdmsp   = proton diffraction radius
c      rdmsn   = neutron diffraction radius
c      rdms    = total diffraction radius
c      surfc   = charge surface thickness
c      surfp   = proton surface thickness
c      surfn   = neutron surface thickness
c      surf    = total surface thickness
c
c     the following fields 'foc' etc are workspace for formfactors:
c      foc0   = charge formfactor, including folding only with electric
c               nucleon distribution
c      foc    = complete charge formfactor, including also the folding with
c               the anomalous magnetic moment
c      fop    = formfactor of the proton distribution
c      fon    = formfactor of the neutorn distribution
c      fo     = formfactor of the total nucleon distribution
c
      parameter (kfrm=50)        ! maximum dimension of formfactor fields
c
      dimension foc(kfrm),fop(kfrm),fon(kfrm),fo(kfrm),foc0(kfrm)
c
      dimension rhop(ngrid),rhon(ngrid)         ! input information
      dimension sjsp(ngrid),sjsn(ngrid)         ! input information
      dimension rhochr(ngrid),formf(kfrm)       ! results
c
      common /radpar/delq,rleng,iqmax           ! to communicate to subr.
c
      data  zero,half,one,pi/0.0,0.5,1.0,3.1415926/
c
c-------------------------------------------------------------------------------
c
c     preparation of auxiliary parameters
c
      anumb  = one*(nprot+nneut)
      kqmax  = min(kfrm,ngrid*3/4)       ! ???
      iqmax  = kqmax
c      iqmax  = ngrid-1
      rleng  = (ngrid-1)*rstep
      delq   = 3.1415926/rleng
      frpirs = 4*pi*rstep
c
c     compute 'cmwid' as width for for c.m.correction on the formfactor
c
      if(ifzpe.eq.0) then
       cmwid  = 1.58/(6.0*anumb**0.66666667)
      else
       cmwid  = 0.375*0.5*(amfpr+amfne)/(anumb*zpe)
      endif
c
c      formfactor by fourier-bessel transformation
c
c                                                 q=0.0 separately
      sump   = zero
      sumn   = zero
c      sumsp  = zero
      r      = zero
       do 500 ir=2,ngrid
       r     = rstep+r
       sump   = r*r*rhop(ir)+sump
       sumn   = r*r*rhon(ir)+sumn
  500  continue
      foc0(1) = frpirs*sump
      fop(1)  = foc0(1)
      fon(1)  = frpirs*sumn
      fo(1)   = fop(1)+fon(1)
      foc(1)  = foc0(1)
c                                                      all other q
      do 530 iq=2,iqmax
      q1     = iq*delq-delq
      call fourf(q1,
     &           rhop,rhon,sjsp,sjsn,
     &           rstep,ngrid,ifrso,ifform,cmwid,
     &           foc0(iq),fop(iq),fon(iq),fo(iq),foc(iq))
  530 continue
c
c     finding first diffraction radii 'rdms..' at first zero
c     of the formfactor,
c     and second diffraction radii 'rsec..' at second zero.
c
c     the endings indicate the type of radius:
c      ....c  = radius of the complete charge distribution
c      ....c0 = radius of the charge distribution without spin-orbit effects
c      ....p  = radius of the proton distribution
c      ....n  = radius of the neutron distribution
c      ....   = radius of the total mass distribution
c
      rdmsc  = 4.493/qzero(foc,1)
      rsecc  = 7.7253/qzero(foc,int(4.493/(delq*rdmsc))+2)
      if(iwr.ge.-1) then
       rdmsc0 = 4.493/qzero(foc0,1)
       rsecc0 = 7.7253/qzero(foc0,int(4.493/(delq*rdmsc0))+2)
       rdms   = 4.493/qzero(fo,1)
       rsec   = 7.7253/qzero(fo,int(4.493/(delq*rdms))+2)
       if(nprot.gt.0) then
        rdmsp  = 4.493/qzero(fop,1)
        rsecp  = 7.7253/qzero(fop,int(4.493/(delq*rdmsp))+2)
       else
        rdmsp  = zero
        rsecp  = zero
       endif
       if(nneut.gt.0) then
        rdmsn  = 4.493/qzero(fon,1)
        rsecn  = 7.7253/qzero(fon,int(4.493/(delq*rdmsn))+2)
       else
        rdmsn  = zero
        rsecn  = zero
       endif
      endif
c      rdmsca=rdmsc/anumb**.3333333
c
c      evaluating surface width from height of first maximum.
c      the 'q1'="momentum_of_first_maximum" is estimated from
c      'rdmsc' as in a square-well density-distribution.
c
      q1     = 5.600/rdmsc
      qr     = q1*rdmsc
      surfc  = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &              *3.*nprot/(qr*qr*fbint(q1,foc))))))
      if(iwr.ge.-1) then
       qr     = q1*rdms
       surf   = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &           *3.*(nprot+nneut)/(qr*qr*fbint(q1,fo))))))
       if(nprot.gt.0) then
        qr     = q1*rdmsp
        surfp  = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &                *3.*nprot/(qr*qr*fbint(q1,fop))))))
       else
        surfp  = zero
       endif
       if(nneut.gt.0) then
        qr     = q1*rdmsn
        surfn  = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &             *3.*nneut/(qr*qr*fbint(q1,fon))))))
       else
        surfn  = zero
       endif
       qr     = q1*rdmsc
       surfc0 = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &            *3.*nprot/(qr*qr*fbint(q1,foc0))))))
      endif
c
c     rms-radii from curvature of formfactor at q=0
c
      dqt    = delq/4.0
      dqtd2  = 6.0/(dqt*dqt)
      rmsc   = sqrt(abs(dqtd2*(1.25*foc(1)
     &                  -1.33333333*fbint(dqt,foc)
     &                  +fbint((dqt+dqt),foc)/12.0)/foc(1)))
      if(iwr.ge.-1) then
       rmsc0  = sqrt(abs(dqtd2*(1.25*foc0(1)
     &                   -1.33333333*fbint(dqt,foc0)
     &                +fbint((dqt+dqt),foc0)/12.0)/foc0(1)))
       if(nprot.gt.0) then
        rmsp   = sqrt(abs(dqtd2*(1.25*fop(1)
     &                   -1.33333333*fbint(dqt,fop)
     &                +fbint((dqt+dqt),fop)/12.0)/fop(1)))
       else
        rmsp   = zero
       endif
       if(nneut.gt.0) then
        rmsn   = sqrt(abs(dqtd2*(1.25*fon(1)
     &                   -1.33333333*fbint(dqt,fon)
     &                +fbint((dqt+dqt),fon)/12.0)/fon(1)))
       else
        rmsn   = zero
       endif
        rms    = sqrt(abs(dqtd2*(1.25*fo(1)
     &                   -1.33333333*fbint(dqt,fo)
     &                +fbint((dqt+dqt),fo)/12.0)/fo(1)))
      endif
c
c     save complete charge formfactor on 'formf'.
c
      do 540 iq=1,iqmax
  540 formf(iq)=foc(iq)
c
c     testprints of spin-orbit current
c
      if(ifrso.eq.1 .and. iwr.ge.0) then
       write(06,'(a,20(/1x,10g12.4))') ' foc0:',
     &                 (foc0(iq),iq=1,iqmax)
       write(06,'(a,20(/1x,10g12.4))') ' foc:',
     &                 (foc(iq),iq=1,iqmax)
      endif
c
c     charge density by fourier back-transformation
c
      dpi3   = half/(pi*pi)
c                                                         r=0.0 separately
      sum    = zero
      do 550 iq=2,iqmax
      qa     = delq*iq-delq
  550 sum    = qa*qa*formf(iq)+sum
      rhochr(1) = abs(sum*delq*dpi3)
c                                                         now all other r
      r      = zero
      do 590 ir=2,ngrid
      r      = rstep+r
c
      dc     = cos(r*delq)
      ds     = sin(r*delq)
      c      = one
      s      = zero
      sum    = zero
      qa     = zero
      do 595 iq=2,iqmax
      qa     = delq+qa
      x      = c*dc-s*ds
      s      = s*dc+c*ds
      c      = x
      sum    = qa*s*foc(iq)+sum
  595 continue
      rhochr(ir) = abs(sum*delq*dpi3/r)
c
  590 continue
c
c     print results
c
      if(iwr.ge.-1)
     & write(06,'(16x,a,5(/1x,a,4g13.4))')
     &   ' diff.radius    surface     sec.radius   rms-radius',
     &   ' charge density:',rdmsc,surfc,rsecc,rmsc,
     &   '   mass density:',rdms,surf,rsec,rms,
     &   ' proton density:',rdmsp,surfp,rsecp,rmsp,
     &   'neutron density:',rdmsn,surfn,rsecn,rmsn,
     &   ' charge no s.o.:',rdmsc0,surfc0,rsecc0,rmsc0
c
      if(iwr.ge.-1) write(06,'(1h0,a,50(/1x,10g13.5))')
     &            'testprint rhochr:', (rhochr(i),i=1,ngrid)
c
      return
      end
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
c-----fbint ------------ part of the hartree-fock package ----------------------
c
      function fbint(qin,formfa)
c
      implicit double precision (a-h,o-z)
c
c     evaluates and returns formfactor at the momentum 'qeff=q*rleng' by
c     fourier-bessel interpolation
c     from formfactors at fourier-grid points given on array 'formfa'
c     in the span '1..iqmax'.
c
      parameter (kfrm=50)        ! maximum dimension of formfactor fields
      dimension formfa(kfrm)
c
      common /radpar/delq,rleng,iqmax           ! to communicate to subr.
c
c     'precis' defines the neighbourhood of a pole-term for which the
c     routine switches to an explicit laurent expansion.
c
      data precis/1.e-2/
      data  zero,one,pi/0.0,1.0,3.1415926/
c
c-----------------------------------------------------------------------
c
      qeff   = qin*rleng
      xfx    = (qeff/pi)
      xfx2   = xfx*xfx
      ixfx   = (xfx+precis)
c
c     the fourier-bessel interpolation runs into problems with round-off
c     errors if 'qeff' comes too close to a multiple of 'pi'. we switch
c     to an explicit treatment of the pole-term in that case (i.e. first
c     case in the if-closure).
c
      if(abs(xfx-ixfx) .lt. precis) then
c                                             pole term separately
       ipole  = ixfx+1
       signum = -one
       sum    =  zero
       do 20 i=2,iqmax
       signum = -signum
   20  if(i.ne.ipole)
     *  sum    = formfa(i)*signum/(one-xfx2/float((i-1)*(i-1))) + sum
       signum = (mod(ipole,2)*2-1)
       qdif   = (xfx-ixfx)
       sum2   = signum*formfa(ipole)*(ixfx*ixfx)/((ixfx+ixfx)+qdif)
       qdif   = qdif*pi
       sinq   = (one-qdif*qdif*0.166666667)*signum
       fbint  = 2.0*(sum*qdif+sum2*pi)*sinq/qeff
      else
c                                             normal case
       signum = -one
       sum    = zero
       do 10 i=2,iqmax
       signum = -signum
   10  sum    = formfa(i)*signum/(one-xfx2/float((i-1)*(i-1))) + sum
       fbint  = (sum+sum)*sin(qeff)/qeff
      endif
c
      return
      end
c-----qzero ------------ part of the hartree-fock package ----------------------
c
      function qzero(foin,iqstrt)
c
      implicit double precision (a-h,o-z)
c
c     evaluates and returns the momentum at which the formfactor
c     given on array 'foin' has a zero.
c     scan for zero starts at the field position 'iqstrt'.
c     'iqmax' is the actual length of the array 'foin'.
c     'delq' is the stepsize in the momentum.
c
      parameter (kfrm=50)        ! maximum dimension of formfactor fields
      dimension foin(kfrm)
c
      common /radpar/delq,rleng,iqmax           ! to communicate to subr.
c
c     some numerical parameters for termination of the search
c
      data itq0mx,endq0/3,0.0001/
      data  zero,one/0.0,1.0/
c
c---------------------------------------------------------------------
c
c     scan changing sign in 'foin'.
c      'iq1' is position before the zero and 'iq1' after the zero.
c
      iqmaxm = iqmax-1
      do 10 iscan=iqstrt,iqmaxm
      if(foin(iscan)*foin(iscan+1).le.zero) goto 19
   10 continue
   19 continue
      iq1   = iscan
      iq2   = iq1+1
      if(foin(iq1).le.zero) iq2   = iq1-1
      q1    = iq1*delq-delq
      q2    = iq2*delq-delq
      f1    = foin(iq1)
      f2    = foin(iq2)
c
c     come closer to 'foin'=0.0 by repeated linear, interpolation.
c
      q3     = q1-f1*(q2-q1)/(f2-f1)
      do 20 itq0=1,itq0mx
      f3     = fbint(q3,foin)
      if(f1*f3.gt.zero) then
       q1     = q3
       f1     = f3
      else
       q2     = q3
       f2     = f3
      endif
      q3old  = q3
      q3     = q1-f1*(q2-q1)/(f2-f1)
      if(abs(q3-q3old).lt.endq0) goto 29
   20 continue
   29 continue
c
c     improve position of the zero by quadratic interpolation:
c
c                         the coefficients of the interpolating binomial
c
      f1     = fbint(q1,foin)
      f2     = fbint(q2,foin)
      f3     = fbint(q3,foin)
      ca     = one/((q1-q2)*(q2-q3)*(q3-q1))
      c11     = ca*(q3*(q3-q2)*q2*f1+q1*(q1-q3)*q3*f2+q2*(q2-q1)*q1*f3)
      c2     = ca*((q2*q2-q3*q3)*f1+(q3*q3-q1*q1)*f2+(q1*q1-q2*q2)*f3)
      c3     = ((q3-q2)*f1+(q1-q3)*f2+(q2-q1)*f3)*ca
c
c                         the two zeroes of the interpolating binomial
c
      qa1    = ( sqrt(abs(c2*c2-4.*c11*c3))-c2)/(c3+c3)
      qa2    = (-sqrt(abs(c2*c2-4.*c11*c3))-c2)/(c3+c3)
c
c                         select the zero which is closer to q3
c
      q      = qa1
      if(abs(qa1-q3).gt.abs(qa2-q3)) q=qa2
c
c     return zero thus found.
c
      qzero  = q
c
      return
      end

