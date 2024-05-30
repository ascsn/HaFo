c========================================================================
c========================================================================
c    
c     Here follows the Pairing and Lipkin-Nogami package of
c      -->   G.A.Leander and W.Nazarewicz, march 1985   <--
c     revised by P.-G. Reinhard, june 1994
c
c     BCS is standard, but includes correction 'gpair*v^2' to
c     the s.p. energies.
c
c     The method of Lipkin and Nogami and coworkers is explained,
c     e.g. in Nucl. Phys. A201 (1973) 357 .
c
c========================================================================
c-----bcslip-------------------------------------------------------------------
c
      subroutine bcslip(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,
     &                  vv2,vln,epair,ecorr,entot,erspln,suv2rs,
     &                  ilipng,iftest)
c
ctest      implicit double precision (a-h,o-z)
c
ctest      parameter (kmaxln=46)
      include 'params.inc'                       ! activate that for
      parameter (kmaxln=ksl2)                    ! HaFo package
c
c------------------------------------------------------------------------------
c
c     Program to solve the Lipkin-Nogami, or BCS, equations
c       gap equation:      0 = f1 = delta-gpair*sum(u*v)
c       particle number:   0 = f2 = sum(2v**2)-n
c       p.n. fluctuation:  0 = f3 = xlam2-gpair/4*(sum(u**3*v)sum(v**3*u)-
c                                   sum(u**4*v**4))/((sum(u**2*v**2))**2 -
c                                   sum(u**4*v**4))
c     the f3 is omitted for BCS.
c
c     The 3 equations are solved by Newton's method.
c     The matrix of derivatives is computed by finite differences.
c     The internal numerical parameters are explained below.
c
c------------------------------------------------------------------------------
c
c     The input parameters are:
c
c       spen   = array of s.p. energies
c       ph     = array of multiplicities, e.g. a j=7/2 state has ph = 4 .
c       nspd   = lowest level in pairing
c       nspg   = highest level in pairing
c       npart  = particle number
c       gpair  = strength of pairing force
c       delta  = initial value for pairing gap
c       xlam   = initial value for Fermi energy
c       xlam2  = initial value for fluctuation constraint
c       vv2    = array of internal occupation probabilities (initial)
c       erspln = Skyrme response energy for 'lambda_2', needed for
c                Skyrme H.F. codes, to be provided in a separte routine.
c                In all other cases set to 0.0 .
c       suv2rs = u^2*v^2 from Skyrme response routine, for normalization.
c                In all other cases set to 1.0 .
c       ilipng = switch between BCS and Lipkin-Nogami
c                 1 --> Lipkin-Nogami
c                 0 --> BCS
c       iftest = switch for testprints
c
c     The output parameters are:
c       delta  = result for pairing gap
c       xlam   = result for Fermi energy
c       xlam2  = result for fluctuation constraint
c       vv2    = array of internal occupation probabilities
c       vlb    = array of occupation probabilities for external observables
c       epair  = naive pairing energy , -delta^2/gpair
c       ecorr  = consistent pairing energy as correction for the hf energy:
c                =-delta**2/gpair
c                 -gpair*sum{v(i)**4-v(i)**2}
c                 -4*xlam2*sum_(u(i)*v(i))**2
c                = -(4*xlam2-gpair)*sum_(u(i)*v(i))**2
c       entot  = total s.q.p. energy (= s.p. sum - g * npart/2 for normal state)
c
c------------------------------------------------------------------------------
c
c
      dimension spen(kmaxln),vv2(kmaxln),vln(kmaxln),ph(kmaxln)
c
ctest      dimension derv(3,3)
c
c------------------------------------------------------------------------------
c
c     some useful constants, allow easy switch single <--> double precision
c
      data zero,one/0.0,1.0/
c
c------------------------------------------------------------------------------
c
c     iteration parameters:
c      itbcsl  = max.nr. iterations in 3 pairing parameters
c      iterv2  = max.nr. iterations in 'vv2'
c      epsf1   = precision on pairing equation 'f1'
c      epsf2   = precision on pairing equation 'f2'
c      epsf3   = precision on pairing equation 'f3'
c      epsv2   = precision on 'vv2'
c      epsdv2  = limit on change of energy in determining 'vv2'
c      detmin  = limiting minimum value for the determinant in linear eq.
c      ddllim  = limit on stepsize in 'delta'
c      dxllim  = limit on stepsize in 'xlam'
c      dx2lim  = limit on stepsize in 'xlam2'
c      dampin  = general reduction of the step, top stabilize
c      delstp  = step in 'delta' for first derivative in finite difference
c      xlstp   = step in 'xlam' for first derivative in finite difference
c      xl2stp  = step in 'xlam2' for first derivative in finite difference
c
      data itbcsl,iterv2,epsf1,epsf2,epsf3,epsv2,epsdv2
     &    /    50,    50,     5*1.0D-8                 /
      data detmin,ddllim,dxllim,dx2lim,dampin
     &    /1.0D-1,  0.2,    3.0,   1.0, 0.8/
      data delstp,xlstp,xl2stp
     &    /  0.01,  0.01,  0.01/
c
c---------------------------------------------------------------------
c
c     preparation and tests
c
      if(nspg.gt.kmaxln) stop ' last level > "KMAXLN" in BCSLIP'
      if(ilipng.ne.1) then
       xlam2 = zero
       dxlam2 = zero
      endif
c
c    iteration of the three pairing parameters 'delta', 'xlam', and 'xlam2' .
c
      do 1000 itbcs=1,itbcsl
      if(iftest.gt.0)
     &   write(6,'(a,i5,a)')  ' iteration nr.',itbcs,' :'
c
c     initial call to pairing equations
c
      call eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,
     &           vv2,vv2,erspln,suv2rs,iftest,ilipng,
     &           iterv2, epsv2,epsdv2,dampin,
     &           f1,f2,f3)
      if(iftest.gt.0) then
        write(6,'(a,3g12.3/a,3g12.3)')
     &   ' delta,xlam,xlam2=',delta,xlam,xlam2,
     &   '         f1,f2,f3=',f1,f2,f3
      endif
c
c     test termination
c
      if(itbcs.gt.1 .and. abs(f1).lt.epsf1 .and. abs(f2).lt.epsf2
     &   .and. abs(f3).lt.epsf3) goto 1099
c
c     derivatives by finite differences
c
      call eqdef(spen,ph,nspd,nspg,npart,gpair,delta*(one+delstp),
     &            xlam,xlam2,vv2,vln,erspln,suv2rs,iftest,ilipng,  ! vln = aux
     &            iterv2, epsv2,epsdv2,dampin,
     &            f1del,f2del,f3del)
      call eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam+xlstp,
     &            xlam2,vv2,vln,erspln,suv2rs,iftest,ilipng,       ! vln = aux
     &            iterv2, epsv2,epsdv2,dampin,
     &            f1xl,f2xl,f3xl)
      if(ilipng.eq.1)
     & call eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2+xl2stp,
     &            vv2,vln,erspln,suv2rs,iftest,ilipng,             ! vln = aux
     &            iterv2, epsv2,epsdv2,dampin,
     &            f1xl2,f2xl2,f3xl2)
      ddelac = delta*delstp
      ddf1   = (f1del-f1)/ddelac
      ddf2   = (f2del-f2)/ddelac
      dlf1   = (f1xl-f1)/xlstp
      dlf2   = (f2xl-f2)/xlstp
      if(ilipng.eq.1) then
        dl2f1  = (f1xl2-f1)/xl2stp
        dl2f2  = (f2xl2-f2)/xl2stp
        ddf3   = (f3del-f3)/ddelac
        dlf3   = (f3xl-f3)/xlstp
        dl2f3  = (f3xl2-f3)/xl2stp
c
c     solve the linear equations for the case of 3 steps
c
c   define matrix in orders
c      derv(1,1) = ddf1
c      derv(2,1) = ddf2
c      derv(3,1) = ddf3
c      derv(1,2) = dlf1
c      derv(2,2) = dlf2
c      derv(3,2) = dlf3
c      derv(1,3) = dl2f1
c      derv(2,3) = dl2f2
c      derv(3,3) = dl2f3
c
c    matrix of inverse is computed explicitely and stored on bij
c
        det    =(ddf1*dlf2*dl2f3-ddf1*dl2f2*dlf3-ddf2*dlf1*dl2f3
     #          +ddf2*dl2f1*dlf3+ddf3*dlf1*dl2f2-ddf3*dl2f1*dlf2)
        det    = sign(max(abs(det),detmin),det)        ! lower limit for 'det'
        detinv = one/det
        b11 = (dlf2*dl2f3-dl2f2*dlf3)*detinv
        b12 = -(dlf1*dl2f3-dl2f1*dlf3)*detinv
        b13 = (dlf1*dl2f2-dl2f1*dlf2)*detinv
        b21 = -(ddf2*dl2f3-dl2f2*ddf3)*detinv
        b22 = (ddf1*dl2f3-dl2f1*ddf3)*detinv
        b23 = -(ddf1*dl2f2-dl2f1*ddf2)*detinv
        b31 = (ddf2*dlf3-dlf2*ddf3)*detinv
        b32 = -(ddf1*dlf3-dlf1*ddf3)*detinv
        b33 = (ddf1*dlf2-dlf1*ddf2)*detinv
ctest      if(iftest.ge.1) write(6,'(a,3g12.4,2(/19x,3g12.4))')
ctest     &   ' derivative matrix:',derv
ctest      call  MATIN2 (derv,3,3,NERROR,DET)
ctest      if(iftest.ge.1)
ctest     &  write(6,'(a,9g12.4)')
ctest     &   'd:',derv,
ctest     &   'b:',b11,b21,b31,b12,b22,b32,b13,b23,b33
ctest      ddelta = -derv(1,1)*f1-derv(1,2)*f2-derv(1,3)*f3
ctest      dxlam  = -derv(2,1)*f1-derv(2,2)*f2-derv(2,3)*f3
ctest      dxlam2 = -derv(3,1)*f1-derv(3,2)*f2-derv(3,3)*f3
        ddelta = -b11*f1-b12*f2-b13*f3
        dxlam  = -b21*f1-b22*f2-b23*f3
        dxlam2 = -b31*f1-b32*f2-b33*f3
        if(iftest.gt.0)
     &   write(6,'(a,4g12.3)')
     &    ' det,ddelta,dxlam,dxlam2=',det,ddelta,dxlam,dxlam2
      else
c
c     case of BCS, step in two quantities
c
        det    = -(dlf1*ddf2-dlf2*ddf1)
        det    = sign(max(abs(det),detmin),det)      ! lower limit for 'det'
        ddelta = (f2*dlf1-f1*dlf2)/det
        dxlam  = (f1*ddf2-f2*ddf1)/det
      endif
c
c     limit stepsize if necessary
c
      ddelta = dampin*sign(min(abs(ddelta/delta),ddllim),ddelta)
      dxlam  = dampin*sign(min(abs(dxlam ),dxllim),dxlam )
      if(ilipng.eq.1)
     &   dxlam2 = dampin*sign(min(abs(dxlam2),dx2lim),dxlam2)
c
c     step forward
c
      delta  = delta*exp(ddelta)
      xlam   = xlam  + dxlam
      if(ilipng.eq.1)
     &    xlam2  = xlam2 + dxlam2
c
 1000 continue
c
      if(iftest.ge.0) write(6,'(a)')
     & ' no convergence of pairing parameters in BCSLIP'
c
 1099 continue
      if(iftest.ge.0) write(6,'(a,i4,6(1pg11.3))')
     & ' BCSLIP: iter,delta,xlam,xlam2,f1,f2,f3=',
     &   itbcs,delta,xlam,xlam2,f1,f2,f3
c
c     compute some observables after convergence
c
      suv    = zero
      sv2    = zero
      sv4    = zero
      suv2   = zero
      sspe   = zero
      do 200 k=nspd,nspg
      v2     = vv2(k)
      uv     = sqrt(v2-v2*v2)
      vln(k) = uv                                     ! auxiliary storage
      xmult  = ph(k)
      suv    = suv    + xmult*uv
      sv2    = sv2    + xmult*v2
      sv4    = sv4    + xmult*v2*v2
      suv2   = suv2   + xmult*(v2 - v2*v2)
      sspe   = sspe   + xmult*spen(k)*v2
  200 continue
cold      ecorr  = -delta*delta/gpair-gpair*sv4-4.0*xlam2*suv2+gpair*sv2
      ecorr  = -delta*delta/gpair-(4.0*xlam2-gpair)*suv2
      entot  =  sspe - gpair * npart/2
      epair1 = -gpair*sv2
      efluc  = -4.*xlam2*suv2
      epair  = -gpair*suv
cold      etot   = 2.0*sspe+epair1+efluc+epair
cold      epairr = epair
      epair1 = epair1+efluc+epair
c
c           corrected coefficients, vln(i)=<a^+_i a_i>_ln
c
      if(ilipng.eq.1) then
        sum2=0.
        sum4=0.
        sumc=0.
        sumc2=0.
c
        do 120 k=nspd,nspg
        v2    = vv2(k)
        xmult = ph(k)
cold        uv    = sqrt(v2-v2*v2)
        uv    = vln(k)
        s     = uv+uv
        c     = v2+v2-one
        sum2  = sum2  + xmult*s**2
        sum4  = sum4  + xmult*s**4
        sumc  = sumc  + xmult*s*s*c
        sumc2 = sumc2 + xmult*(s*c)**2
  120   continue
c
        od=sum2*0.25/(sumc2*sum2-sumc*sumc+0.25*sum2*(sum2*sum2-sum4))
c
        do 130 k=nspd,nspg
        v2    = vv2(k)
        uv    = vln(k)
cold        uv    = sqrt(v2-v2*v2)
        s     = uv+uv
        c     = v2+v2-one
        vln(k)= v2+od*s*s*(c*sum2-sumc)
  130   continue
      else
        do 135 k=nspd,nspg
        vln(k)= vv2(k)
  135   continue
      endif
c
      return
      end
c-----eqdef-------------------------------------------------------------------
c
      subroutine eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,
     &                 vv2in,vv2out,erspln,suv2rs,iftest,ilipng,
     &                 iterv2, epsv2,epsdv2,dampin,
     &                 f1,f2,f3)
c
c---------------------------------------------------------------------
c
ctest      implicit double precision (a-h,o-z)
c
ctest      parameter (kmaxln=46)
      include 'params.inc'
      parameter (kmaxln=ksl2)
c
c---------------------------------------------------------------------
c
c     routine to compute the three pairing equations, as explained
c     in the header of 'BCSLIP'.
c
c     input:
c       spen   = array of s.p. energies
c       ph     = array of multiplicities, e.g. a j=7/2 state has ph = 4 .
c       nspd   = lowest level in pairing
c       nspg   = highest level in pairing
c       npart  = particle number
c       gpair  = strength of pairing force
c       delta  = pairing gap
c       xlam   = initial value for Fermi energy
c       xlam2  = initial value for fluctuation constraint
c       vv2in  = array of initial internal occupation probabilities
c       erspln - Skyrme response energy for 'lambda_2'
c       suv2rs -   u^2*v^2 from response routine, for normalization
c       iftest = switch for testprints
c       ilipng = switch between BCS and Lipkin-Nogami,
c                 1 --> f3=0
c       iterv2 = max.nr. iterations in 'vv2'
c       epsv2  = precision on 'vv2'
c       epsdv2 = limit on change of energy in determining 'vv2'
c       dampin = general reduction of the step, top stabilize
c     output:
c       f1     = delta-gpair*sum(u*v) =  gap equation
c       f2     = sum(2v**2)-n         =  particle number condition
c       f3     = xlam2-gpair/4*(sum(u**3*v)sum(v**3*u)-
c                               sum(u**4*v**4))/((sum(u**2*v**2))**2 -
c                               sum(u**4*v**4))
c              = p.n. fluctuation
c       vv2out = new occupation weights
c
c---------------------------------------------------------------------
c
c
      dimension spen(kmaxln),vv2in(kmaxln),vv2out(kmaxln),ph(kmaxln)
c
      data zero,half,one/0.0,0.5,1.0/
c
c---------------------------------------------------------------------
c
c     preparation
c
      fnn    = npart
      if(nspg.gt.kmaxln) stop ' last level > "kmaxln"'
      d2     = delta*delta
      alpha  = 4.*xlam2-gpair
      alphd2 = alpha*half*d2
c
c     reset accumulators for combinations of 'u' and 'v'.
c
      sv2   = zero                                  ! sum v^2
      suv   = zero                                  ! sum u*v
      if(ilipng.eq.1) then
       suv4  = zero                                  ! sum u^4 v^4
       suv2  = zero                                  ! sum u^2 v^2
       su3v  = zero                                  ! sum u*v u^2
       suv3  = zero                                  ! sum u*v v^2
      endif
c
      do 500 k=nspd,nspg                            ! each state separately
c
c     iteration loop to compute occupation 'vv2'
c
      e     = spen(k)-xlam
      v2    = vv2in(k)                                ! start from previous value
      x     = e+alpha*v2
c
      do 6 it=1,iterv2
      eq    = sqrt(d2+x*x)
      v2    = half-half*x/eq
      f     = e+alpha*v2-x
      derf  = -alphd2/(eq*eq*eq)-one
ctest      dx    = dampin*f/derf         ! opionally switch on if unstable
      dx    = f/derf
      if(abs(f).lt.epsv2 .and. abs(dx).lt.epsdv2) go to 7
      x     = x-dx
    6 continue
c
      write(6,'(a)')
     &  ' maximal number of iterations exceeded in "EQDEF"'
c
    7 continue
      if(iftest.gt.1)
     & write(6,'(a,i5,a,i5,g12.3)')
     &  ' state nr.',k,' :  iterv,v2',it,v2
      vv2out(k)=v2
c
c     accumulate sums with new 'v2'
c
      xmult = ph(k)
      u2    = one-v2
      uv    = sqrt(u2*v2)
      sv2   = sv2  + xmult*v2
      suv   = suv  + xmult*uv
      if(ilipng.eq.1) then
       suv4  = suv4 + xmult*u2*u2*v2*v2
       suv2  = suv2 + xmult*u2*v2
       su3v  = su3v + xmult*uv*u2
       suv3  = suv3 + xmult*uv*v2
      endif
c
  500 continue
c
c     final assembling of the three pairing equations
c
      f1     = delta-gpair*suv                    ! gap equation
      f2     = sv2+sv2-fnn                        ! particle number
c     
c     the numerator from Skyrme response is a bit special.
c     it was computed with an older "v2". The corresponding
c     "sum V2*U2" is stored in "suv2rs". In order to avoid
c     all to frequent recomputation of that expensive piece,
c     we rescale it here just by the new "suv2".
c
      if(ilipng.eq.1) then
       denom  = (suv2*suv2-suv4)                   !  denominator/32
       enume1 = -gpair/4.0*(su3v*suv3-suv4)        !  1. numerator/32
       enume2 = suv2*suv2*erspln/(suv2rs*suv2rs*32.0)   !  2. numerator/32
       f3     = xlam2+(enume1+enume2)/denom        ! p.n. fluctuations
cold   f3=xlam2-gpair/4.*(su3v*suv3-suv4)/(suv2*suv2-suv4)
       endif
c
      return
      end
c================================================================
      function gavrg(ibcs,iz,in,iprot,rhosc,n1,n2,xmeff,a0)
c
c         gavrg calculates the strength of the monopole pairing
c         interaction using the average gap method. the average
c         gap model is based on moller and nix, nucl. phys.a536
c         (1992) 20. the average pair density is computed using
c         the semi-classical approximation  for   the  harmonic
c         oscillator model modified  by  the  isospin-dependent
c         factor derived based on the wigner-kirkwood expansion
c         for the woods-saxon potential (+ s.o. + coulomb).
c
c         input:    iz:    proton  number
c                   in:    neutron number
c                   iprot: 0 (neutrons) or 1 (protons)
c                   n1:    number of first h.f. level used
c                          in the bcs (ln) calculations.
c                          (preferably n1=1)
c                   n2:    number of last h.f. level used
c                          in the bcs (ln) calculations.
c                          (preferably n2=npart, i.e., iz or in)
c                   xmeff: effective mass (=1 in the woods
c                          saxon model)
c                   ibcs:  1 (bsc) or .ne.1 (ln)
c                   a0:    surface parameter in a Saxon Woods shape,
c                          is 0.58 * surface width from formfactor
c
c        output:    gavrg: pairing strength
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c      written by w. nazarewicz, 1994                               c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      implicit real*8 (a-h,o-z)
      parameter (kmaxln=2024)            ! dimensions of pairing fields
c
      dimension alphap(2),bethap(2)

      data alphap /5.13,-2.35/,
     :     bethap /-1.28e-2,-2.38e-4/
c
c          pair density at the fermi level
c
ctest      write(6,'(a/4i4,1f8.3,2i4,2f8.3)')
ctest     &   ' ibcs,iz,in,iprot,rhosc,n1,n2,xmeff,a0',
ctest     &     ibcs,iz,in,iprot,rhosc,n1,n2,xmeff,a0
      a=iz+in
      ntot=iz
      if (iprot.eq.0) ntot=in
      alpha=alphap(iprot+1)
      betha=bethap(iprot+1)
      fac=(2.*ntot/a)**(1.d0/3.)
      part3=(3.*ntot)**(2.0/3.)
      hom=41/a**(1.0/3)*fac
      tz=(in-iz)/a
      tz2=tz*tz
      rhosc=part3/hom/(1.-0.25/part3)/2.*(1.+alpha*tz2+betha*a*tz2)
c
c                here: correction for the effective mass
c
       r=1.14*a**(1.d0/3)
       eta=a0/r
c
       xm=log(xmeff)
       pieta=(3.1415926*eta)**2
       xmm=xmeff*(1.-eta*xm)**3*(1.+2*eta*xm*pieta/(1.+pieta))
c
      rhosc=rhosc*xmm**1.5
      robar=rhosc
c

      if (ibcs.eq.1) then
c
c              bcs variant
c
         if (iprot.eq.0) delavr=13.3/sqrt(a)
         if (iprot.eq.1) delavr=13.9/sqrt(a)

      else
c
c              ln variant
c
         if (iprot.eq.0) delavr=8.6/sqrt(a)
         if (iprot.eq.1) delavr=9.9/sqrt(a)
c
      endif
c
      np=ntot/2
c
      y1=(-np+n1-1)/robar
      y2=(-np+n2  )/robar
c
      gavrg=2./robar/(log(sqrt(y2**2+delavr**2)+y2)-
     *                log(sqrt(y1**2+delavr**2)+y1))
c
      return
      end

