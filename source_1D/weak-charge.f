c-----radii    ----------------------------------------------------------
c
      subroutine radii(rhop,rhon,sjsp,sjsn,nprot,nneut,rstep,ngrid,
     &                 ifrso,cmwid,rhochr,
     &                  rmsc,rmsp,rmsn,rms,
     &                  rdiffrc,rdiffrp,rdiffrn,rdiffr,
     &                  surfc,surfp,surfn,surf)

c
      implicit double precision (a-h,o-z)
c
c     Compute rms-radius, diffraction-radius and surface-thickness 
c     for proton-, neutron-, total- and charge-density.
c     Densities are to be given on an equidistant 1D radial grid.
c
c     The charge-density includes folding with intrinsic electric and
c     magnetic formfactors of proton and neutron.
c     The contribution of the anomalous magnetic moment is included
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
c      cmwid   = width for c.m. correction
c                 self-consistent c.m. correction uses
c                    cmwid  = 3/(8*<P^2>) = 0.375*(hbar^2/2M)/(ekin)
c                 rough estimate uses
c                    cmwid  = 1.58/(6.0*(nprot+nneut)**0.66666667)
c
c     Output variables are:
c      rhochr  = charge density inluding all nucleon formfactors
c      rmsc    = charge r.m.s. radius
c      rmsp    = proton r.m.s. radius
c      rmsn    = neutron r.m.s. radius
c      rms     = total r.m.s. radius
c      rdiffrc = charge diffraction radius (sharp- or box-radius)
c      rdiffrp = proton diffraction radius
c      rdiffrn = neutron diffraction radius
c      rdiffr  = total diffraction radius
c      surfc   = charge surface thickness
c      surfp   = proton surface thickness
c      surfn   = neutron surface thickness
c      surf    = total surface thickness
c
c     the following fields 'foc' etc are workspace for formfactors:
c      foc    = complete charge formfactor, including also the folding with
c               the anomalous magnetic moment
c      fop    = formfactor of the proton distribution
c      fon    = formfactor of the neutorn distribution
c      fo     = formfactor of the total nucleon distribution
c
      dimension foc(ngrid),fop(ngrid),fon(ngrid),fo(ngrid)
      dimension fsop(ngrid),fson(ngrid)
c
      dimension rhop(ngrid),rhon(ngrid)
      dimension sjsp(ngrid),sjsn(ngrid)
      dimension rhochr(ngrid)
c
      data  zero,half,one,pi/0D0,0.5D0,1D0,3.14159265358979D0/
c
c--------------------------------------------------------------------------
c
c     prepare basic formfactors
c
      call  fbtrans(rhop,rhon,rstep,ngrid,cmwid,fop,fon,delq,iqmax)
      if(ifrso.eq.1) then
        call  fbtrans(sjsp,sjsn,rstep,ngrid,cmwid,fsop,fson,delq,iqmax)
      else
        fsop = 0D0
        fson = 0D0
      endif
c
c     compose electro-magnetic formfactor
c
      call formf_elmag(fop,fon,fsop,fson,delq,iqmax,ifrso,foc)
c
c     evaluate charge density and all form parameters
c     (charge density is used as work space in the first three calls)
c
      call form_parameters(fop,delq,iqmax,rstep,ngrid,
     &                     rhochr,rmsp,rdiffrp,surfp)
      call form_parameters(fon,delq,iqmax,rstep,ngrid,
     &                     rhochr,rmsn,rdiffrn,surfn)
      fo(:) = fop(:)+fon(:)
      call form_parameters(fo,delq,iqmax,rstep,ngrid,
     &                     rhochr,rms,rdiffr,surf)
      call form_parameters(foc,delq,iqmax,rstep,ngrid,
     &                     rhochr,rmsc,rdiffrc,surfc)
c
c     write protocol
c
      write(6,'(a)') ' form parameters from new subroutines:'
        write(06,'(16x,a,4(/1x,a,4(1pg13.5)))')
     &   ' diff.radius    surface     rms-radius   total nr.',
     &   ' charge density:',rdiffrc,surfc,rmsc,foc(1),
     &   '   mass density:',rdiffr,surf,rms,fo(1),
     &   ' proton density:',rdiffrp,surfp,rmsp,fop(1),
     &   'neutron density:',rdiffrn,surfn,rmsn,fon(1)
c
      return
      end
c
c-----formf_elmag----------------------------------------------------------
c
      subroutine formf_elmag(fop,fon,fsop,fson,delq,iqmax,ifrso,foc)
c
      implicit double precision (a-h,o-z)
c
c     Composes the charge formfactor from given proton and neutron 
c     formfactors.
c
c     Input variables are:
c      fop     = 1D array with proton formfactor
c      fon     = 1D array with neutron formfactor
c      fsop    = 1D array with proton spin-orbit formfactor
c      fson    = 1D array with neutron spin-orbit formfactor
c      delq    = spacing in momentum-space grid
c      iqmax   = maximum nr. of points in momentum-space grid
c      ifrso   = 1 --> chrarge density shall contain contribtion from
c                      spin-orbit current
c     Output::
c      foc     = 1D array with charge formfactor
c
      dimension fop(iqmax),fon(iqmax) 
      dimension fsop(iqmax),fson(iqmax)
      dimension foc(iqmax)
c
      data  zero,half,one,pi/0D0,0.5D0,1D0,3.14159265358979D0/
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
      foc(1) = fop(1)          
      do iq=2,iqmax
       q1 = delq*(iq-1)
c
c
c      isoscalar and isovector electrical formactors on 'schses'
c      and 'schsev', recoupled then to proton and neutron electric
c      formfactors  'schsep' and 'schsen'. proton and neutron
c      magnetic formfactor on 'schsmp' and 'schsmn'.
c
       qa     = q1*q1
       schses = sa1/(1D0+qa/sm1)+sa2/(1D0+qa/sm2)+sa3/(1D0+qa/sm3)
     &         +sa4/(1D0+qa/sm4)
       schsev = va1/(1D0+qa/vm1)+va2/(1D0+qa/vm2)+va3/(1D0+qa/vm3)
     &         +va4/(1D0+qa/vm4)
       schsep = 0.5D0*(schses+schsev)
       schsen = 0.5D0*(schses-schsev)
       schsm  = ( gp1/(1D0+qa/gm1)+gp2/(1D0+qa/gm2)+gp3/(1D0+qa/gm3)
     &           +gp4/(1D0+qa/gm4) )
       schsmp = schsm*xmgpr - schsep*d4m
       schsmn = schsm*xmgne - schsen*d4m
c
       sumc   = fop(iq)*schsep+fon(iq)*schsen
c
c      factor for darwin correction 'darfac'
c      and the factor for the c.m. correction 'exphf'.
c
       if(ifrso.eq.1) then
        darfac = 1D0/sqrt(1D0+qa*d4m)
       else
        darfac = 1D0
       endif
c
c      factor for c.m. correction
c
       foc(iq) = (sumc*darfac + (schsmp*fsop(iq)+schsmn*fson(iq)))
      enddo
c
      return
      end




c-----form_parameters-------------------------------------------------------
c
      subroutine form_parameters(formf,delq,iqmax,rstep,ngrid,
     &                           rhoout,rms,rdiffr,surf)
c
      implicit double precision (a-h,o-z)
c
c     Computes the three form parameters for a given formfactor.
c     Input variables are:
c      formf   = 1D array with formfactor in steps of 'delq'
c      delq    = step size in momentum space
c      iqmax   = maximum nr. of grid points in momentum space
c      rstep   = spacing of the (equidistant) radial grid
c      ngrid   = number of radial grid points
c     output paramaters:
c      rhoout  = local density distribution from Fourier-Bessel back transf.
c      rms     = r.m.s. radius
c      rdiffr  = diffraction radius
c      surf    = surface thickness
c
      real*8 formf(iqmax),rhoout(ngrid)
c
      data  zero,half,one,pi/0D0,0.5D0,1D0,3.14159265358979D0/
c
c
c------------------------------------------------------------------------------
c  
      write(6,'(a)') ' form_parameters entered'
c
c     diffraction radius
c
      rdiffr  = 4.493D0/findqzero(formf,1,iqmax,delq,rstep*(ngrid-1))
c
c     surface thickness
c
      q1     = 5.600D0/rdiffr
      qr     = q1*rdiffr
      finter = fbinterpol(q1,formf,iqmax,delq,rstep*(ngrid-1))
      surf   = sqrt(abs(2D0/(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &              *3D0*formf(1)/(qr*qr*finter)))))
c
c     r.m.s. radius from density (derivative definition numerically unstable)
c     density from Fourier-Bessel transformation
c
c
c     charge density by fourier back-transformation
c
      dpi3   = 0.5D0/(pi*pi)
c                                                         r=0.0 separately
      sum    = zero
      do iq=2,iqmax
        qa     = delq*iq-delq
        sum    = qa*qa*formf(iq)+sum
      enddo
      rhoout(1) = abs(sum*delq*dpi3)
c                                                         now all other r
      r      = zero
      do ir=2,ngrid
        r      = rstep+r
c
        dc     = cos(r*delq)
        ds     = sin(r*delq)
        c      = one
        s      = zero
        sum    = zero
        qa     = zero
        do iq=2,iqmax
          qa     = delq+qa
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
          sum    = qa*s*formf(iq)+sum
        enddo
        rhoout(ir) = abs(sum*delq*dpi3/r)
c
      enddo
c
c     now r.m.s. radius by straightforward r-space integration
c
      rmom0  = zero
      rmom2  = zero
      r      = zero
      do ir=2,ngrid
        r     = rstep+r
        rmom0 = (r*r)*rhoout(ir) + rmom0
        rmom2 = (r*r)*(r*r)*rhoout(ir) + rmom2
      enddo
      rmom2  = 4D0*pi*rstep*rmom2/formf(1)
      rms    = sqrt(rmom2)
      write(6,'(a,2(1pg13.5))')'FORM PARAMS: rms=',rms,sqrt(rmom2/rmom0)
c
      return
      end
c
c-----fbtrans ------------ part of the hartree-fock package ----------------------
c
      subroutine fbtrans(rhop,rhon,rstep,ngrid,cmwid,fop,fon,delq,iqmax)
c
      implicit double precision (a-h,o-z)
c
c     Input variables are:
c      rhop    = 1D array of radial proton density (rhop(1) --> at r=0)
c      rhon    = 1D array of radial neutron density
c      rstep   = spacing of the (equidistant) radial grid
c      ngrid   = number of radial grid points
c      cmwid   = width parameter for center of mass correction
c     output paramaters:
c      fop     = proton formfactor
c      fon     = neutron    "
c      delq    = spacing of the momentum-space grid
c      iqmax   = maximum number of grid points in momentum space
c
      real*8 rhop(ngrid),rhon(ngrid) 
      real*8 fop(ngrid),fon(ngrid)   
c
      data  zero,half,one,pi/0D0,0.5D0,1D0,3.14159265358979D0/
c
c
c----------------------------------------------------------------
c
c     prepare parameters
c
      frpirs = 4D0*pi*rstep
      rleng  = (ngrid-1)*rstep
      delq   = pi/rleng
      iqmax  = ngrid
c
c     Fourier-Bessel transformation at q=0
c
      sump   = zero
      sumn   = zero
      r      = zero
      do ir=2,ngrid
        r     = rstep+r
        sump   = r*r*rhop(ir)+sump
        sumn   = r*r*rhon(ir)+sumn
      enddo
      fop(1) = frpirs*sump
      fon(1) = frpirs*sumn

c
c     Fourier-Bessel transformation for q>0
c
      do iq=2,iqmax
        q1     = iq*delq-delq
        dc     = cos(q1*rstep)
        ds     = sin(q1*rstep)
        r      = zero
        sump   = zero
        sumn   = zero
        c      = one
        s      = zero
         do i=2,ngrid
          r      = rstep+r
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
          sump   = r*s*rhop(i)+sump
          sumn   = r*s*rhon(i)+sumn
         enddo
        exphf  = exp(cmwid*q1*q1)
        fop(iq) = exphf*frpirs*sump/q1
        fon(iq) = exphf*frpirs*sumn/q1
      enddo
c
      return
      end
c-----fbinterpol -------------------------------------------------------
c
      function fbinterpol(qin,formfa,iqmax,delq,rlengbox)
c
      implicit double precision (a-h,o-z)
c
c     Evaluates and returns formfactor at the momentum 'qeff=q*rleng' by
c     fourier-bessel interpolation
c     from formfactors at fourier-grid points given on array 'formfa'
c     in the span '1..iqmax' and in steps of 'delq'.
c
      dimension formfa(iqmax)
c
c     'precis' defines the neighbourhood of a pole-term for which the
c     routine switches to an explicit laurent expansion.
c
      data precis/1.d-9/          ! /1.e-2/
      data  zero,one,pi/0D0,1D0,3.14159265358979D0/
c
c-----------------------------------------------------------------------
c
      qeff   = qin*rlengbox
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
       do i=2,iqmax
         signum = -signum
         if(i.ne.ipole)
     *    sum    = formfa(i)*signum/(one-xfx2/float((i-1)*(i-1))) + sum
       enddo
       signum = (mod(ipole,2)*2-1)
       qdif   = (xfx-ixfx)
       sum2   = signum*formfa(ipole)*(ixfx*ixfx)/((ixfx+ixfx)+qdif)
       qdif   = qdif*pi
       sinq   = (one-qdif*qdif*0.166666667)*signum
       fbinterpol  = 2.0*(sum*qdif+sum2*pi)*sinq/qeff
      else
c                                             normal case
       signum = -one
       sum    = zero
       do i=2,iqmax
         signum = -signum
         sum    = formfa(i)*signum/(one-xfx2/float((i-1)*(i-1))) + sum
       enddo
       fbinterpol  = (sum+sum)*sin(qeff)/qeff
      endif
c
      return
      end
c-----findqzero ------------ part of the hartree-fock package ----------------------
c
      function findqzero(foin,iqstrt,iqmax,delq,rlengbox)
c
      implicit double precision (a-h,o-z)
c
c     evaluates and returns the momentum at which the formfactor
c     given on array 'foin' has a zero.
c     scan for zero starts at the field position 'iqstrt'.
c     'iqmax' is the actual length of the array 'foin'.
c     'delq' is the stepsize in the momentum.
c
      dimension foin(iqmax)
c
c     some numerical parameters for termination of the search
c
      data itq0mx,endq0/3,0.0001/
      data  zero,one/0D0,1D0/
c
c---------------------------------------------------------------------
c
c     scan changing sign in 'foin'.
c      'iq1' is position before the zero and 'iq1' after the zero.
c
      iqmaxm = iqmax-1
      do iscan=iqstrt,iqmaxm
        if(foin(iscan)*foin(iscan+1).le.zero) goto 19
      enddo
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
      do itq0=1,itq0mx
        f3     = fbinterpol(q3,foin,iqmax,delq,rlengbox)
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
      enddo
   29 continue
c
c     improve position of the zero by quadratic interpolation:
c
c                         the coefficients of the interpolating binomial
c
      f1     = fbinterpol(q1,foin,iqmax,delq,rlengbox)
      f2     = fbinterpol(q2,foin,iqmax,delq,rlengbox)
      f3     = fbinterpol(q3,foin,iqmax,delq,rlengbox)
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
      findqzero  = q
c
      return
      end
