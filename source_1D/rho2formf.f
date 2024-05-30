c
      program rho2formf
c
      implicit double precision (a-h,o-z)
c
      parameter (kfrm=299)        ! maximum dimension of formfactor fields
      parameter (kgrid=kfrm)
c
c     computes radial formfactor for density read in from unit 5
c
      dimension rin(kgrid),rho(kgrid),formf(kgrid)
c
      data  zero,half,one,pi/0.0,0.5,1.0,3.1415926/
c
      common /radpar/delq,rleng,iqmax           ! to communicate to subr.
c
c-----------------------------------------------------------------------------
c
c     reading density
c
      read(5,'(1x)') 
      read(5,'(1x)') 
      do ir=1,kgrid
        read(5,*,err=199,end=199) rin(ir),rho(ir)
        ngrid = ir
      enddo
      stop ' insufficient dimensions '
 199  continue
      rstep = rin(2)-rin(1)
      write(6,*) '# test grid:',rin(ngrid)-rstep*(ngrid-1)
c
c     compute formfactor 
c
      iqmax  = ngrid
      rleng  = (ngrid-1)*rstep
      delq   = 3.1415926D0/rleng
      frpirs = 4*pi*rstep
c                                                 q=0.0 separately
      sum    = 0.0D0
      r      = 0.0D0
      do ir=2,ngrid
        r     = rstep+r
        rho(ir) = rho(ir)/100.0
        sum   = r*r*rho(ir)+sum
      enddo
      formf(1)  = frpirs*sum
c                                                      all other q
      do iq=2,iqmax
        q1     = iq*delq-delq
        call fourf0(q1,rho,rstep,ngrid,formf(iq))
      enddo
c
c     print
c
      write(6,'(a)') '# q formfactor '
      do iq=1,iqmax
        q1     = iq*delq-delq
        write(6,'(1x,2(1pg13.5))') q1,formf(iq)
      enddo
c
c
      write(6,'(a)') '     ','     ','# q formfactor '
      write(6,'(1x,2(1pg13.5))') 0.0D0,formf(1)
      do iq=2,10*iqmax
        q1     = (iq-1)*delq*0.1D0
        write(6,'(1x,2(1pg13.5))') q1,fbint(q1,formf)
      enddo
c
      stop
      end

c-----fourf0------------ part of the hartree-fock package ----------------------
c
      subroutine fourf0(q1,rho,rstep,ngrid,foc)
c
      implicit double precision (a-h,o-z)
c
c     Input variables are:
c      q1    : radial momentum at which the fourier-transformation
c              is to be evaluated
c      rhoc    = array of radial charge density (rhop(1) --> at r=0)
c      rstep   = spacing of the (equidistant) radial grid
c      ngrid   = number of radial grid points
c     output paramaters:
c      foc  :   charge    
c
      dimension rho(ngrid)
c
c      data  zero,one,pi/0.0,1.0,3.1415926/
      data  zero,half,one,pi/0.0,0.5,1.0,3.1415926/
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
      sum    = zero
      c      = one
      s      = zero
      do i=2,ngrid
        r      = rstep+r
        x      = c*dc-s*ds
        s      = s*dc+c*ds
        c      = x
        sum    = r*s*rho(i)+sum
      enddo
      foc    = frpirs*sum/q1
c
      return
      end
c-----fbint ------------ part of the hartree-fock package ----------------------
c
      double precision function fbint(qin,formfa)
c
      implicit double precision (a-h,o-z)
c
c     evaluates and returns formfactor at the momentum 'qeff=q*rleng' by
c     fourier-bessel interpolation
c     from formfactors at fourier-grid points given on array 'formfa'
c     in the span '1..iqmax'.
c
      parameter (kfrm=299)        ! maximum dimension of formfactor fields
      dimension formfa(kfrm)
c
      common /radpar/delq,rleng,iqmax           ! to communicate to subr.
c
c     'precis' defines the neighbourhood of a pole-term for which the
c     routine switches to an explicit laurent expansion.
c
      data precis/1.d-9/          ! /1.e-2/
      data  zero,one,pi/0.0d0,1.0d0,3.14159265358979d0/
c
c-----------------------------------------------------------------------
c
      if(iqmax.gt.kfrm) stop ' too large grid in FBINT '
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
       fbint  = 2.0*(sum*qdif+sum2*pi)*sinq/qeff
      else
c                                             normal case
       signum = -one
       sum    = zero
       do i=2,iqmax
         signum = -signum
         sum    = formfa(i)*signum/(one-xfx2/float((i-1)*(i-1))) + sum
       enddo
       fbint  = (sum+sum)*sin(qeff)/qeff
      endif
c
      return
      end
