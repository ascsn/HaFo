      program tstbcs
c
      implicit double precision (a-h,o-z)
c
      parameter (kmaxln=46)
c
      dimension spen(kmaxln),vv2(kmaxln),vln(kmaxln),ph(kmaxln)
c
      data zero,half,one/0.0,0.5,1.0/
      data ph/1.,2.,1.,3.,2.,1.,4.,3.,2.,1.,5.,4.,3.,2.,1.,6.,5.,
     &        4.,3.,2.,1.,7.,6.,5.,4.,3.,2.,1.,8.,7.,6.,5.,4.,3.,2.,1.,
     &        9.,8.,7.,6.,5.,4.,3.,2.,1.,10./
ctest      data ph/kmaxln*1.0/
c
      ifold  = 0
      iftang = 0
      ifextd = 0
      ilipng = 1
      gpair = 0.05
      dele   = 1.0
      nstate = 10
      npart  = 20 
      sump   = 0.0
      do 10 n=1,nstate
      spen(n) = n*dele 
      if(n.le.6) then
       vv2(n)  = 1.0
       spen(n) = -dele
      else
       vv2(n)  = 0.0
       spen(n) = dele
      endif
   10 continue
c
      delta = 3.0
      xlam  = 0.0
      xlam2 = 0.5
      if(ifold.eq.1) then
c
c        the old pairing routine , no gpair*v^2
c
       call pairol(spen,vv2,npart,nstate,gpair,xlam,epair,delta,1,
     &                  partnm,0054,1.0D-8,iab)
       write(6,'(2(/a,2(1pg12.4)))')
     &  ' pairing parameters: delta,xlam=',delta,xlam,
     &  ' results:                  epair=',epair
       write(6,'(a,40(/1x,5g12.3))')
     & ' occupation v2:',(vv2(k),k=1,nstate)
      endif
c
      if(ifextd.ne.0) then
c
c        Lipkin-Nogami with separate external 'delta' iteration,
c         compact version of the old code
c
       call bcslpd(spen,ph,1,nstate,npart,gpair,delta,xlam,xlam2,
     &            vv2,vln,epair,ecorr,entot,zero,one)
       write(6,'(2(/a,3(1pg12.4)))')
     &  ' pairing parameters: delta,xlam,xlam2=',delta,xlam,xlam2,
     &  ' results:            epair,ecorr,entot=',epair,ecorr,entot
       write(6,'(a,40(/1x,5g12.3))')
     & ' occupation v2:',(vv2(k),k=1,nstate)
      endif
      if(iftang.eq.0) then
c
c        Lipkin-Nogami with Newton iteration of all 3 parameters
c          'delta', 'xlam', 'xlam2'.
c        ! This is the actual version
c        contains also switch to BCS with gpair*v^2
c
       call bcslip(spen,ph,1,nstate,npart,gpair,delta,xlam,xlam2,
     &                  vv2,vln,epair,ecorr,entot,zero,one,ilipng,1)
      else
c
c        Lipkin-Nogami with Newton iteration of all 3 parameters
c          'delta', 'xlam', 'xlam2'.
c        This version has errors in d/dxlam2,
c        but it tries a direct computation of the derivatives
c
       call bcslol(spen,ph,1,nstate,npart,gpair,delta,xlam,xlam2,
     &            vv2,vln,epair,ecorr,entot,zero,one)
      endif
c
      write(6,'(2(/a,3(1pg12.4)))')
     &  ' pairing parameters: delta,xlam,xlam2=',delta,xlam,xlam2,
     &  ' results:            epair,ecorr,entot=',epair,ecorr,entot
      write(6,'(a,40(/1x,5g12.3))')
     & ' occupation v2:',(vv2(k),k=1,nstate)
c
      stop
      end
c-----bcslpd-------------------------------------------------------------------
c
      subroutine bcslpd(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,
     &                  vv2,vln,epair,ecorr,entot,erspln,suv2rs)
c
      implicit double precision (a-h,o-z)
c
      parameter (kmaxln=46)
c
c     variant with interlaced iteration:
c     outer iteration on 'delta' and inner on 'xlam' & 'xlam2'.
c
c       gpair  = strength of pairing force
c       delta  = pairing gap
c       erespo = response energy , for the fluctuation part
c       npart  = particle number
c       nspd   = lowest level in pairing
c       nspg   = highest level in pairing
c       spen   = array of s.p. energies
c       ph     = array of multiplicities, e.g. a j=7/2 state has ph = 4 .
c       vv2    = array of internal occupation probabilities
c     erspln - Skyrme response energy for 'lambda_2'             !mod
c              this information is carried further through on 'erespo'
c              via common /bcslop/
c     suv2rs -   u^2*v^2 from response routine, for normalization  !mod
c              this information is carried further through on 'suv2ro'
c              via common /bcslop/
c     epair  - pairing energy.
c     ecorr  - correction for the hf energy: 
c                 ecorr=-d**2/g-g*sum_v(i)**4-4*xl2*sum_(u(i)*v(i))**2+g*nsp
c     entot  - total s.q.p. energy (= s.p. sum - g * npart/2 for normal state)
c
c    the routine solves the three equations
c     gap equation:      0 = f1 = delta-gpair*sum(u*v)
c     particle number:   0 = f2 = sum(2v**2)-n
c
c     p.n. fluctuation:  0 = f3 = xlam2-gpair/4*(sum(u**3*v)sum(v**3*u)-
c                                 sum(u**4*v**4))/((sum(u**2*v**2))**2 -
c                                  sum(u**4*v**4))
c
c
      dimension spen(kmaxln),vv2(kmaxln),vln(kmaxln),ph(kmaxln)
      dimension v2conv(3,50)  ! for testprints
c
      data zero,half,one/0.0,0.5,1.0/
c
c     iteration parameters
c      itbcsl  = max.nr. iterations in 3 pairing parameters
c      iterv2  = max.nr. iterations in 'vv2'
c      epsf1   = precision on pairing equation 'f1'
c      epsf2   = precision on pairing equation 'f2'
c      epsf3   = precision on pairing equation 'f3'
c      epsv2   = precision on 'vv2'
c      epsdv2  = limit on change of energy in determining 'vv2'
c      detmin  = limiting minimum value for the determinant in linear eq.
c      deltmn  = minimum value for gap 'delta'
c      dampin  = general reduction of the step, top stabilize
c      iftest  = switch for testprints
c
      data itbcsl,iterv2,epsf1,epsf2,epsf3,epsv2,epsdv2
     &    /    50,    50,     5*1.0D-6                 /
      data detmin,ddllim,dxllim,dx2lim,dv2lim,deltmn,dampin
     &    /1.0D-3,  1.0,    0.5,   1.0,   1.0,1.0E-4, 0.5/
      data iftest
     &    /     1/
      data xsmall
     &    /1.D-20/
      data delstp
     &    / 0.001/
c
c---------------------------------------------------------------------
c
      fnn    = npart
      if(nspg.gt.kmaxln) stop ' last level > "kmaxln"'
      dv2lma = dv2lim*gpair                ! limit v2-step with 'gpair'
c
c  first iteration loop: iterating 'delta'
c
      do 1001 itdel=1,itbcsl
      delsav = delta
      if(iftest.ne.0)
     &   write(6,'(a,i5,a)')  ' ITERATION NR.',itdel,' :'
      do 1002 itvard=0,1
      delta  = delsav*(one+itvard*delstp)
c
c  second iteration loop: iteration of the two pairing parameters
c                            'xlam' ,  and 'xlam2' .
c
      do 1000 itbcs=1,itbcsl
      if(iftest.gt.1)
     &   write(6,'(a,i5,a)')  ' iteration nr.',itbcs,' :'
      d2     = delta*delta
      alpha  = 4.*xlam2-gpair
ctest      alpha  = zero        ! test
c
c---------------------------------------------------------------
c
c   First block:
c
c   compute the three criteria for given pairing parameters
c
c---------------------------------------------------------------
c
c
c     reset accumulators for combinations of 'u' and 'v'.
c
      sv2   = zero                                  ! sum v^2
      suv   = zero                                  ! sum u*v
      suv4  = zero                                  ! sum u^4 v^4
      suv2  = zero                                  ! sum u^2 v^2
      su3v  = zero                                  ! sum u*v u^2
      suv3  = zero                                  ! sum u*v v^2
c
      do 500 k=nspd,nspg                            ! each state separately
c
c     iteration loop to compute occupation 'vv2'
c
      e     = spen(k)-xlam
      v2    = vv2(k)                                ! start from previous value
cold      x     = (half-v2)*sqrt(d2/(v2-v2*v2))         ! restore from initial v2
      x     = e+alpha*v2
c
      do 6 it=1,iterv2
      eq    = sqrt(d2+x*x)
      v2    = 0.5*(1.-x/eq)
      f     = e+alpha*v2-x
      v2conv(1,it) = x
      v2conv(2,it) = v2
      v2conv(3,it) = f
      derf  = -alpha*0.5*d2/eq**3-1.
      dx    = f/derf
      dx    = dampin*sign(min(abs(dx),dv2lma),dx)
      if(abs(f).lt.epsv2 .and. abs(dx).lt.epsdv2) go to 7
      x     = x-dx
    6 continue
c
      write(6,'(a)')
     &  ' maximal number of iterations exceeded in "EQDEF"'
      if(iftest.gt.2) then
         write(6,'(a)') '  iter   x    v2    f '
         do 90 it=1,iterv2
   90    write(6,'(1x,i5,3g12.3)')
     &   it,v2conv(1,it),v2conv(2,it),v2conv(3,it)
      endif
c
    7 continue
      if(iftest.gt.1)
     & write(6,'(a,i5,a,i5,4(1pg11.4))')
     &  ' state ',k,' :  iterv,v2,x,f,eq=',it,v2,x,f,eq
      vv2(k)=v2
c
c     accumulate sums with new 'v2'
c
      xmult = ph(k)
      u2    = one-v2
      uv    = sqrt(u2*v2)
      sv2   = sv2  + xmult*v2
      suv   = suv  + xmult*uv
      suv4  = suv4 + xmult*u2*u2*v2*v2
      suv2  = suv2 + xmult*u2*v2
      su3v  = su3v + xmult*uv*u2
      suv3  = suv3 + xmult*uv*v2
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
      denom  = max((suv2*suv2-suv4),xsmall)       !  denominator/32
      enume1 = -gpair/4.0*(su3v*suv3-suv4)        !  1. numerator/32
      enume2 = suv2*suv2*erspln/(suv2rs*suv2rs*32.0)   !  2. numerator/32
      if(iftest.ge.2) write(6,'(a,3g12.4)')
     & ' in ETOT: enume1,enume2,denom ',
     &            enume1,enume2,denom
      f3     = xlam2+(enume1+enume2)/denom        ! p.n. fluctuations
cold  f3=xlam2-gpair/4.*(su3v*suv3-suv4)/(suv2*suv2-suv4)
c
      if(iftest.gt.1) then
        write(6,'(a,4g12.3)')
     &   ' xlam,xlam2,f2,f3==',xlam,xlam2,f2,f3
      endif
      if(abs(f2).lt.epsf2 .and. abs(f3).lt.epsf3) goto 1099
c
c---------------------------------------------------------------
c
c   Second block:
c
c   compute the matrix of gradients explicitely
c
c---------------------------------------------------------------
c
c     accumulate the various sums of gradients on 'v2' with u and v weights
c
      slm    = zero          ! sum d/dl uv
      slv    = zero          ! sum d/dl v^2
      slum   = zero          ! sum uv d/dl u^2
      slmu   = zero          ! sum u^2 d/dl uv
      slvm   = zero          ! sum uv d/dl v^2
      slmv   = zero          ! sum v^2 d/dl uv
      sluuv2 = zero          ! sum v^4 u^2 d/dl u^2
      slvvu2 = zero          ! sum u^4 v^2 d/dl v^2
      sluv   = zero          ! sum v^2 d/dl u^2
      slvu   = zero          ! sum u^2 d/dl v^2
      sdm    = zero          ! sum d/dd uv
      sdv    = zero          ! sum d/dd v^2
      sdum   = zero          ! sum uv d/dd u^2
      sdmu   = zero          ! sum u^2 d/dd uv
      sdvm   = zero          ! sum uv d/dd v^2
      sdmv   = zero          ! sum v^2 d/dd uv
      sduuv2 = zero          ! sum v^4 u^2 d/dd u^2
      sdvvu2 = zero          ! sum u^4 v^2 d/dd v^2
      sduv   = zero          ! sum v^2 d/dd u^2
      sdvu   = zero          ! sum u^2 d/dd v^2
c
      do 100 k=nspd,nspg
c
      e      = spen(k)-xlam
      v2     = vv2(k)
cold      x      = (half-v2)*sqrt(d2/(v2-v2*v2))
      x      = e+alpha*v2
      u2     = one-v2
      uv     = sqrt(u2*v2)
      eq     = delta/(uv+uv)
      d2e    = d2/(eq*eq)
      denom  = ph(k)/(eq+eq-alpha*d2e)             ! denominator in derivative * multiplicity
      dlv    = denom*d2e                           ! d/dlambda v^2
      ddv    = denom*delta*x/(eq*eq)               ! d/ddelta v^2
      dlu    = -dlv                                ! d/dlambda u^2
      ddu    = -ddv                                ! d/ddelta u^2
      dlm    = dlv*(u2-v2)/(uv+uv)                 ! d/dlambda uv
      ddm    = ddv*(u2-v2)/(uv+uv)                 ! d/ddelta uv
c
      slm    = slm    + dlm
      slv    = slv    + dlv
      slum   = slum   + dlu*uv
      slmu   = slmu   + dlm*u2
      slvm   = slvm   + dlv*uv
      slmv   = slmv   + dlm*v2
      sluuv2 = sluuv2 + dlu*u2*v2*v2
      slvvu2 = slvvu2 + dlv*v2*u2*u2
      sluv   = sluv   + dlu*v2
      slvu   = slvu   + dlv*u2
      sdm    = sdm    + ddm
      sdv    = sdv    + ddv
      sdum   = sdum   + ddu*uv
      sdmu   = sdmu   + ddm*u2
      sdvm   = sdvm   + ddv*uv
      sdmv   = sdmv   + ddm*v2
      sduuv2 = sduuv2 + ddu*u2*v2*v2
      sdvvu2 = sdvvu2 + ddv*v2*u2*u2
      sduv   = sduv   + ddu*v2
      sdvu   = sdvu   + ddv*u2
  100 continue
c
c     compute the gradients , reuse old sums from loop '500' .
c       note: d/dl2 f1 = 0  ,  d/dl2 f2 = 0 , d/dl2 f3 = 1 .
c       --> linear response presently not implemented !
c
      ddf1   = one-gpair*sdm
      dlf1   = -gpair*slm
      ddf2   = sdv+sdv
      dlf2   = slv+slv
      denom  = one/(suv2*suv2-suv4)
      enume1 = -gpair/4.0*(su3v*suv3-suv4)
      dlu4v4 = sluuv2+slvvu2
      dlu4v4 = dlu4v4+dlu4v4
      ddu4v4 = sduuv2+sdvvu2
      ddu4v4 = ddu4v4+ddu4v4
      ddf3   = -0.25*gpair*denom*
     &          ( (sdum+sdmu)*suv3 + (sdvm+sdmv)*su3v - ddu4v4
     &           -enume1*denom*(2.0*(sduv+sdvu)-ddu4v4)  )
      dlf3   = -0.25*gpair*denom*
     &          ( (slum+slmu)*suv3 + (slvm+slmv)*su3v - dlu4v4
     &           -enume1*denom*(2.0*(sluv+slvu)-dlu4v4)  )
c
c     solve the linear equations for the three steps
c
cold      det    = -(dlf1*ddf2-dlf2*ddf1)
cold      det    = sign(max(abs(det),detmin),det)        ! lower limit for 'det'
cold      ddelta = (f2*dlf1-f1*dlf2)/det
cold      dxlam  = (f1*ddf2-f2*ddf1)/det
          dxlam  = -f2/dlf2
cwrong      dxlam2 = -f3+ddelta*ddf3+dxlam*dlf3
cold      dxlam2 = -f3-ddelta*ddf3-dxlam*dlf3
      dxlam2 = -f3-dxlam*dlf3
      dxlam2 = -f3
      if(iftest.gt.1)
     & write(6,'(a,4g12.3)')
     &  ' det,ddelta,dxlam,dxlam2=',det,ddelta,dxlam,dxlam2
c
c     limit stepsize if necessary
c
cold      ddelta = dampin*sign(min(abs(ddelta),ddllim),ddelta)
      dxlam  = dampin*sign(min(abs(dxlam ),dxllim),dxlam )
      dxlam2 = dampin*sign(min(abs(dxlam2),dx2lim),dxlam2)
c
c     step forward
c
cold      delta  = abs(delta + ddelta)
      xlam   = xlam  + dxlam
      xlam2  = xlam2 + dxlam2
ctest      xlam2  = zero
c
 1000 continue
      write(6,'(a)')
     & ' no convergence of "xlam" & "xlam2" in BCSLIP'
 1099 continue
c
c     save results for 'delta' step
c
      if(itvard.eq.0) then
        f1bas  = f1
        xlsav  = xlam
        xl2sav = xlam2
      else
        f1var  = f1
      endif
 1002 continue
c
c     compute next 'delta' by gradient step
c
      delta  = delsav
      if(abs(f1bas).lt.epsf1) goto 1098
      ddf1   = (f1var-f1bas)/(delsav*delstp)
      ddelta = -f1bas/ddf1
      ddelta = dampin*sign(min(abs(ddelta/delta),ddllim),ddelta)
      delta  = delta*exp(ddelta)
      if(iftest.ge.1) write(6,'(a,i4,3(1pg12.4))')
     & ' Delta variation: itdel,f1,ddelta,delta=',
     &   itdel,f1bas,ddelta,delta
 1001 continue
      write(6,'(a)')
     & ' no convergence of "delta" in BCSLIP'
 1098 continue
      xlam   = xlsav
      xlam2  = xl2sav
c
c     compute some observables after convergence
c
      sv4    = zero
      suv2   = zero
      sstate = zero
      sspe   = zero
      do 200 k=nspd,nspg
      v2     = vv2(k)
      xmult  = ph(k)
      sv4    = sv4    + xmult*v2*v2
      suv2   = suv2   + xmult*(v2 - v2*v2)
      sstate = sstate + xmult*one
      sspe   = sspe   + xmult*spen(k)*v2
  200 continue
      ecorr  = -delta*delta/gpair-gpair*sv4-4.0*xlam2*suv2+gpair*sstate
      entot  =  sspe - gpair * npart/2
      epair1 = -gpair*sv2
      efluc  = -4.*xlam2*suv2
      epair  = -gpair*suv
      etot   = 2.0*sspe+epair1+efluc+epair
      epairr = epair
      epair1 = epair1+efluc+epair
c
c           corrected coefficients, vln(i)=<a^+_i a_i>_ln
c
      sum2=0.
      sum4=0.
      sumc=0.
      sumc2=0.
c
      do 120 k=nspd,nspg
      v2    = vv2(k)
      xmult = ph(k)
      uv    = sqrt(v2-v2*v2)
      s     = uv+uv
      c     = v2+v2-one
      sum2  = sum2  + xmult*s**2
      sum4  = sum4  + xmult*s**4
      sumc  = sumc  + xmult*s*s*c
      sumc2 = sumc2 + xmult*(s*c)**2
  120 continue
c
      od=sum2*0.25/(sumc2*sum2-sumc*sumc+0.25*sum2*(sum2*sum2-sum4))
c
      do 130 k=nspd,nspg
      v2    = vv2(k)
      uv    = sqrt(v2-v2*v2)
      s     = uv+uv
      c     = v2+v2-one
      vln(k)= v2+od*s*s*(c*sum2-sumc)
  130 continue
c
      return
      end
c-----bcslol-------------------------------------------------------------------
c
      subroutine bcslol(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,
     &                  vv2,vln,epair,ecorr,entot,erspln,suv2rs)
c
      implicit double precision (a-h,o-z)
c
      parameter (kmaxln=46)
c
c       gpair  = strength of pairing force
c       delta  = pairing gap
c       erespo = response energy , for the fluctuation part
c       npart  = particle number
c       nspd   = lowest level in pairing
c       nspg   = highest level in pairing
c       spen   = array of s.p. energies
c       ph     = array of multiplicities, e.g. a j=7/2 state has ph = 4 .
c       vv2    = array of internal occupation probabilities
c     erspln - Skyrme response energy for 'lambda_2'             !mod
c              this information is carried further through on 'erespo'
c              via common /bcslop/
c     suv2rs -   u^2*v^2 from response routine, for normalization  !mod
c              this information is carried further through on 'suv2ro'
c              via common /bcslop/
c     epair  - pairing energy.
c     ecorr  - correction for the hf energy: 
c                 ecorr=-d**2/g-g*sum_v(i)**4-4*xl2*sum_(u(i)*v(i))**2+g*nsp
c     entot  - total s.q.p. energy (= s.p. sum - g * npart/2 for normal state)
c
c    the routine solves the three equations
c     gap equation:      0 = f1 = delta-gpair*sum(u*v)
c     particle number:   0 = f2 = sum(2v**2)-n
c
c     p.n. fluctuation:  0 = f3 = xlam2-gpair/4*(sum(u**3*v)sum(v**3*u)-
c                                 sum(u**4*v**4))/((sum(u**2*v**2))**2 -
c                                  sum(u**4*v**4))
c
c
      dimension spen(kmaxln),vv2(kmaxln),vln(kmaxln),ph(kmaxln)
      dimension v2conv(3,50)  ! for testprints
c
      data zero,half,one/0.0,0.5,1.0/
c
c     iteration parameters
c      itbcsl  = max.nr. iterations in 3 pairing parameters
c      iterv2  = max.nr. iterations in 'vv2'
c      epsf1   = precision on pairing equation 'f1'
c      epsf2   = precision on pairing equation 'f2'
c      epsf3   = precision on pairing equation 'f3'
c      epsv2   = precision on 'vv2'
c      epsdv2  = limit on change of energy in determining 'vv2'
c      detmin  = limiting minimum value for the determinant in linear eq.
c      deltmn  = minimum value for gap 'delta'
c      dampin  = general reduction of the step, top stabilize
c      iftest  = switch for testprints
c
      data itbcsl,iterv2,epsf1,epsf2,epsf3,epsv2,epsdv2
     &    /   150,    50,     5*1.0D-8                 /
      data detmin,ddllim,dxllim,dx2lim,dv2lim,deltmn,dampin
     &    /1.0D-3,  1.0,    0.5,   1.0,   1.0,1.0E-4, 0.5/
      data iftest
     &    /     1/
      data xsmall
     &    /1.D-20/
c
c---------------------------------------------------------------------
c
      write(6,'(//a/)') ' simultanous & tangential step '
      fnn    = npart
      if(nspg.gt.kmaxln) stop ' last level > "kmaxln"'
      dv2lma = dv2lim*gpair                ! limit v2-step with 'gpair'
c
c  Leading iteration loop: iteration of the three pairing parameters
c                            'delta' , 'xlam' ,  and 'xlam2' .
c
      do 1000 itbcs=1,itbcsl
      if(iftest.ne.0)
     &   write(6,'(a,i5,a)')  ' iteration nr.',itbcs,' :'
      d2     = delta*delta
      alpha  = 4.*xlam2-gpair
ctest      alpha  = zero        ! test
c
c---------------------------------------------------------------
c
c   First block:
c
c   compute the three criteria for given pairing parameters
c
c---------------------------------------------------------------
c
c
c     reset accumulators for combinations of 'u' and 'v'.
c
      sv2   = zero                                  ! sum v^2
      suv   = zero                                  ! sum u*v
      suv4  = zero                                  ! sum u^4 v^4
      suv2  = zero                                  ! sum u^2 v^2
      su3v  = zero                                  ! sum u*v u^2
      suv3  = zero                                  ! sum u*v v^2
c
      do 500 k=nspd,nspg                            ! each state separately
c
c     iteration loop to compute occupation 'vv2'
c
      e     = spen(k)-xlam
      v2    = vv2(k)                                ! start from previous value
cold      x     = (half-v2)*sqrt(d2/(v2-v2*v2))         ! restore from initial v2
      x     = e+alpha*v2
c
      do 6 it=1,iterv2
      eq    = sqrt(d2+x*x)
      v2    = 0.5*(1.-x/eq)
      f     = e+alpha*v2-x
      v2conv(1,it) = x
      v2conv(2,it) = v2
      v2conv(3,it) = f
      derf  = -alpha*0.5*d2/eq**3-1.
      dx    = f/derf
      dx    = dampin*sign(min(abs(dx),dv2lma),dx)
      if(abs(f).lt.epsv2 .and. abs(dx).lt.epsdv2) go to 7
      x     = x-dx
    6 continue
c
      write(6,'(a)')
     &  ' maximal number of iterations exceeded in "EQDEF"'
      if(iftest.gt.2) then
         write(6,'(a)') '  iter   x    v2    f '
         do 90 it=1,iterv2
   90    write(6,'(1x,i5,3g12.3)')
     &   it,v2conv(1,it),v2conv(2,it),v2conv(3,it)
      endif
c
    7 continue
      if(iftest.gt.1)
     & write(6,'(a,i5,a,i5,4(1pg11.4))')
     &  ' state ',k,' :  iterv,v2,x,f,eq=',it,v2,x,f,eq
      vv2(k)=v2
c
c     accumulate sums with new 'v2'
c
      xmult = ph(k)
      u2    = one-v2
      uv    = sqrt(u2*v2)
      sv2   = sv2  + xmult*v2
      suv   = suv  + xmult*uv
      suv4  = suv4 + xmult*u2*u2*v2*v2
      suv2  = suv2 + xmult*u2*v2
      su3v  = su3v + xmult*uv*u2
      suv3  = suv3 + xmult*uv*v2
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
      denom  = max((suv2*suv2-suv4),xsmall)       !  denominator/32
      enume1 = -gpair/4.0*(su3v*suv3-suv4)        !  1. numerator/32
      enume2 = suv2*suv2*erspln/(suv2rs*suv2rs*32.0)   !  2. numerator/32
      write(6,'(a,3g12.4)')
     & ' in ETOT: enume1,enume2,denom ',
     &            enume1,enume2,denom
      f3     = xlam2+(enume1+enume2)/denom        ! p.n. fluctuations
cold  f3=xlam2-gpair/4.*(su3v*suv3-suv4)/(suv2*suv2-suv4)
c
      if(iftest.ne.0) then
        write(6,'(a,4g12.3/a,3g12.3)')
     &   ' delta,xlam,xlam2,alpha=',delta,xlam,xlam2,alpha,
     &   '         f1,f2,f3=',f1,f2,f3
      endif
      if(abs(f1).lt.epsf1 .and. abs(f2).lt.epsf2
     &   .and. abs(f3).lt.epsf3) goto 1099
c
c---------------------------------------------------------------
c
c   Second block:
c
c   compute the matrix of gradients explicitely
c
c---------------------------------------------------------------
c
c     accumulate the various sums of gradients on 'v2' with u and v weights
c
      slm    = zero          ! sum d/dl uv
      slv    = zero          ! sum d/dl v^2
      slum   = zero          ! sum uv d/dl u^2
      slmu   = zero          ! sum u^2 d/dl uv
      slvm   = zero          ! sum uv d/dl v^2
      slmv   = zero          ! sum v^2 d/dl uv
      sluuv2 = zero          ! sum v^4 u^2 d/dl u^2
      slvvu2 = zero          ! sum u^4 v^2 d/dl v^2
      sluv   = zero          ! sum v^2 d/dl u^2
      slvu   = zero          ! sum u^2 d/dl v^2
      sdm    = zero          ! sum d/dd uv
      sdv    = zero          ! sum d/dd v^2
      sdum   = zero          ! sum uv d/dd u^2
      sdmu   = zero          ! sum u^2 d/dd uv
      sdvm   = zero          ! sum uv d/dd v^2
      sdmv   = zero          ! sum v^2 d/dd uv
      sduuv2 = zero          ! sum v^4 u^2 d/dd u^2
      sdvvu2 = zero          ! sum u^4 v^2 d/dd v^2
      sduv   = zero          ! sum v^2 d/dd u^2
      sdvu   = zero          ! sum u^2 d/dd v^2
c
      do 100 k=nspd,nspg
c
      e      = spen(k)-xlam
      v2     = vv2(k)
cold      x      = (half-v2)*sqrt(d2/(v2-v2*v2))
      x      = e+alpha*v2
      u2     = one-v2
      uv     = sqrt(u2*v2)
      eq     = delta/(uv+uv)
      d2e    = d2/(eq*eq)
      denom  = ph(k)/(eq+eq-alpha*d2e)             ! denominator in derivative * multiplicity
      dlv    = denom*d2e                           ! d/dlambda v^2
      ddv    = denom*delta*x/(eq*eq)               ! d/ddelta v^2
      dlu    = -dlv                                ! d/dlambda u^2
      ddu    = -ddv                                ! d/ddelta u^2
      dlm    = dlv*(u2-v2)/(uv+uv)                 ! d/dlambda uv
      ddm    = ddv*(u2-v2)/(uv+uv)                 ! d/ddelta uv
c
      slm    = slm    + dlm
      slv    = slv    + dlv
      slum   = slum   + dlu*uv
      slmu   = slmu   + dlm*u2
      slvm   = slvm   + dlv*uv
      slmv   = slmv   + dlm*v2
      sluuv2 = sluuv2 + dlu*u2*v2*v2
      slvvu2 = slvvu2 + dlv*v2*u2*u2
      sluv   = sluv   + dlu*v2
      slvu   = slvu   + dlv*u2
      sdm    = sdm    + ddm
      sdv    = sdv    + ddv
      sdum   = sdum   + ddu*uv
      sdmu   = sdmu   + ddm*u2
      sdvm   = sdvm   + ddv*uv
      sdmv   = sdmv   + ddm*v2
      sduuv2 = sduuv2 + ddu*u2*v2*v2
      sdvvu2 = sdvvu2 + ddv*v2*u2*u2
      sduv   = sduv   + ddu*v2
      sdvu   = sdvu   + ddv*u2
  100 continue
c
c     compute the gradients , reuse old sums from loop '500' .
c       note: d/dl2 f1 = 0  ,  d/dl2 f2 = 0 , d/dl2 f3 = 1 .
c       --> linear response presently not implemented !
c
      ddf1   = one-gpair*sdm
      dlf1   = -gpair*slm
      ddf2   = sdv+sdv
      dlf2   = slv+slv
      denom  = one/(suv2*suv2-suv4)
      enume1 = -gpair/4.0*(su3v*suv3-suv4)
      dlu4v4 = sluuv2+slvvu2
      dlu4v4 = dlu4v4+dlu4v4
      ddu4v4 = sduuv2+sdvvu2
      ddu4v4 = ddu4v4+ddu4v4
      ddf3   = -0.25*gpair*denom*
     &          ( (sdum+sdmu)*suv3 + (sdvm+sdmv)*su3v - ddu4v4
     &           -enume1*denom*(2.0*(sduv+sdvu)-ddu4v4)  )
      dlf3   = -0.25*gpair*denom*
     &          ( (slum+slmu)*suv3 + (slvm+slmv)*su3v - dlu4v4
     &           -enume1*denom*(2.0*(sluv+slvu)-dlu4v4)  )
c
c     solve the linear equations for the three steps
c
      det    = -(dlf1*ddf2-dlf2*ddf1)
      det    = sign(max(abs(det),detmin),det)        ! lower limit for 'det'
      ddelta = (f2*dlf1-f1*dlf2)/det
      dxlam  = (f1*ddf2-f2*ddf1)/det
cwrong      dxlam2 = -f3+ddelta*ddf3+dxlam*dlf3
      dxlam2 = -f3-ddelta*ddf3-dxlam*dlf3
ctest      dxlam2 = -f3
      if(iftest.gt.0)
     & write(6,'(a,4g12.3/a,4g12.3)')
     &  ' det,ddelta,dxlam,dxlam2=',det,ddelta,dxlam,dxlam2,
     &  ' dlf1,dlf2,ddf1,ddf2=',dlf1,dlf2,ddf1,ddf2
c
c     limit stepsize if necessary
c
cold      ddelta = dampin*sign(min(abs(ddelta),ddllim),ddelta)
      ddelta = dampin*sign(min(abs(ddelta/delta),ddllim),ddelta)
      dxlam  = dampin*sign(min(abs(dxlam ),dxllim),dxlam )
      dxlam2 = dampin*sign(min(abs(dxlam2),dx2lim),dxlam2)
c
c     step forward
c
cold      delta  = abs(delta + ddelta)
      delta  = delta*exp(ddelta)
      xlam   = xlam  + dxlam
      xlam2  = xlam2 + dxlam2
ctest      xlam2  = zero
c
 1000 continue
      write(6,'(a)')
     & ' no convergence of pairing parameters in BCSLIP'
 1099 continue
c
c     compute some observables after convergence
c
      sv4    = zero
      suv2   = zero
      sstate = zero
      sspe   = zero
      do 200 k=nspd,nspg
      v2     = vv2(k)
      xmult  = ph(k)
      sv4    = sv4    + xmult*v2*v2
      suv2   = suv2   + xmult*(v2 - v2*v2)
      sstate = sstate + xmult*one
      sspe   = sspe   + xmult*spen(k)*v2
  200 continue
      ecorr  = -delta*delta/gpair-gpair*sv4-4.0*xlam2*suv2+gpair*sstate
      entot  =  sspe - gpair * npart/2
      epair1 = -gpair*sv2
      efluc  = -4.*xlam2*suv2
      epair  = -gpair*suv
      etot   = 2.0*sspe+epair1+efluc+epair
      epairr = epair
      epair1 = epair1+efluc+epair
c
c           corrected coefficients, vln(i)=<a^+_i a_i>_ln
c
      sum2=0.
      sum4=0.
      sumc=0.
      sumc2=0.
c
      do 120 k=nspd,nspg
      v2    = vv2(k)
      xmult = ph(k)
      uv    = sqrt(v2-v2*v2)
      s     = uv+uv
      c     = v2+v2-one
      sum2  = sum2  + xmult*s**2
      sum4  = sum4  + xmult*s**4
      sumc  = sumc  + xmult*s*s*c
      sumc2 = sumc2 + xmult*(s*c)**2
  120 continue
c
      od=sum2*0.25/(sumc2*sum2-sumc*sumc+0.25*sum2*(sum2*sum2-sum4))
c
      do 130 k=nspd,nspg
      v2    = vv2(k)
      uv    = sqrt(v2-v2*v2)
      s     = uv+uv
      c     = v2+v2-one
      vln(k)= v2+od*s*s*(c*sum2-sumc)
  130 continue
c
      return
      end
c-----bcslip-------------------------------------------------------------------
c
      subroutine bcslip(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,
     &                  vv2,vln,epair,ecorr,entot,erspln,suv2rs,
     &                  ilipng,iftest)
c
      implicit double precision (a-h,o-z)
c
      parameter (kmaxln=46)
chafo      include 'params.inc'                       ! activate that for
chafo      parameter (kmaxln=ksl2)                    ! HaFo package
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
      sstate = zero
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
      sstate = sstate + xmult*one
      sspe   = sspe   + xmult*spen(k)*v2
  200 continue
      ecorr  = -delta*delta/gpair-gpair*sv4-4.0*xlam2*suv2+gpair*sstate
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
      implicit double precision (a-h,o-z)
c
      parameter (kmaxln=46)
chafo      include 'params.inc'
chafo      parameter (kmaxln=ksl2)
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
c-----pairol------------ part of the hartree-fock package ----------------------
c
      subroutine pairol(e,gw,nz,nmax,gp,eferm,epair,delta,ipair,
     &                  partnm,iter,eps,iab)
c
      implicit double precision (a-h,o-z)
c
c     this is a pairing pairing routine. but it is presently used
c     as a routine to determine the thermal occupation of the
c     single particle states. this is achieved in the "constant gap"
c     option, using 'iter<0' as a switch.
c
c     the include file 'params.inc' communicates the dimensioning of
c     the fields 'e' and 'gw' with the parameter 'ksl2'. it allows
c     furthermore to make the whole subroutine double precision by
c     transfering a line 'implicit double precision (a-h,o-z)'.
c
      parameter (kmaxln=46)
c
c     include 'params.inc'
c
c     * * * * * * * * * *     pairing iteration      * * * * * * * * * * * *
c
c     determines the fermi-energy 'eferm' and the gap 'delta'
c     such that the "gap-equation" and the "particle number condition"
c     are fulfilled simultanously. this is achieved by iteration
c     with a gradient scheme (newton's tangential formula)
c     where the gradient is evaluated by finite differences.
c
c     there are two schemes possible:
c     1. the "constant force" approach, where the pairing force is given
c                                       in 'gp'.
c     2. the "constant gap" approach,   where the gap is kept fixed and
c                                       is given in 'delta'.
c     there is only 'eferm' to be determined in the "constant gap"
c     approach. the gap-equation is used at the end to determine an
c     according force 'gp'.
c
c     input via list:
c      e       = array(1..nmax) of single particle energies.
c      gw      = array(1..nmax) of occupation probabilities,
c                input is used as initila guess for the iteration.
c      nz      = desired particle number
c      nmax    = maximum number of single particle states
c                (needs to be smaller than the parameter 'ksl2'
c                 and smaller than 46, for dimensioning reasons).
c      gp      = size of pairing-force, is input in case of
c                "constant force" approach.
c      eferm   = fermi-energy, input is used as initial guess;
c                if eferm=0.0 is input then a suitable initial guess
c                is computed from input 'gw'.
c      epair   = pairing energy
c      delta   = pairing gap, is input in case of "constant gap" approach.
c      ipair   = switch for type of pairing
c                 0 --> constant Delta
c                 1 --> constant Force
c                 3 --> Fermi distribuion using T=Delta
c      iter    = number of pairing iterations
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
c      gp      = size of pairing-force, is output in case of
c                "constant gap" approach.
c      delta   = pairing gap, is output in case of
c                "constant force" approach.
c      eferm   = fermi-energy.
c      sum     = particle number, from adding 'gw' with degeneracy 'ph'.
c
c     some auxiliary fields are:
c      ph      = array(1..46) of the multiplicities of the single
c                particle levels. array is ordered according to
c                the harmonic oscillator scheme, with rising shell
c                number and lowering angular momentum within each shell.
c      dparn   = array(1..3) of variation in the particle number as
c                consequence of 'deferm' and 'ddelta'.
c      dgap    = array(1..3) of variation in the gap equation as
c                consequence of 'deferm' and 'ddelta'.
c      gr      = double array of variation of particle number and of
c                gap equation, used for construction of the tangantial
c                step
c
      dimension e(kmaxln),gw(kmaxln),ph(kmaxln)
      dimension dparn(3),dgap(3),gr(2,2)
c
      data ph/1.,2.,1.,3.,2.,1.,4.,3.,2.,1.,5.,4.,3.,2.,1.,6.,5.,
     &        4.,3.,2.,1.,7.,6.,5.,4.,3.,2.,1.,8.,7.,6.,5.,4.,3.,2.,1.,
     &        9.,8.,7.,6.,5.,4.,3.,2.,1.,10./
ctest      data ph/kmaxln*1.0/
      data gr,dparn,dgap/10*0./
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
      data dampin,deferm,ddelta/0.9,1.0,0.1/
      data efmax,dlmax,detmin,defmin,ddlmin/3.,.6,.1e-3,.1e-3,.1e-3/
c
      data xmaxlw,x10lw/1.0e-30,1.0e-10/
      data one,zero/1.0,0.0/
c
c-------------------------------------------------------------------------------
c
c      set switch variables
c
      iwarng = 0
c
      if(iab.ge.1)
     & write(6,'(/a,i4,4(a,f9.4)/a,10(/1x,10f8.3))')
     &   ' *** pair entered with: iter=',iter,'  delta=',delta,
     &   '  eferm=',eferm,'  gp=',gp,'  deferm=',deferm,
     &   '                        single particle energies:',
     &   (e(i),i=1,nmax)
c
c      prepare fermi-energy 'eferm', if not yet initialised.
c
      if(eferm.eq.zero) then
       isum   = 0
       do 120 it=1,nmax
       itsum  = it
       isum   = int(ph(it)+.000001) +isum
       if(isum+isum-nz.ge.0) goto 110
  120  continue
       stop ' not enough states too reach particle number'
  110  continue
       eferm  = e(itsum)+deferm
       if(iab.ge.1) write(06,'(5x,a,i5,a,g11.4)')
     &   'initialistion: isum=',isum,'  eferm=',eferm
c
      endif
c
      if(ipair.eq.1) then
c
c-------------------------------------------------------------------------------
c
c     case of constant pairing-force, given by 'gp'.
c      gap 'delta', fermi-energy 'eferm', and occupations 'gw'
c      are computed.
c
c-------------------------------------------------------------------------------
c
       if(ipair.ne.0 .and. ipair.ne.1 .and. ipair.ne.3)
     &   stop ' invalid option "ipair" in subroutine "pairol"'
       delti  = delta
       if(abs(delti).lt.x10lw) delti = one
       do 200 it=1,abs(iter)
c
c       gap-equation & particle-number at 'eferm' & 'delti' (for k=1)
c       and at 'eferm' varied (k=2) and at 'delta' varied (k=3).
c
        do 210 ka=1,3
c                                          actual fermi-energy and gap
        k      = 4-ka
        if(k.eq.1) then
         elam   = eferm
         delact = delti
        elseif(k.eq.2) then
         elam   = eferm+deferm
         delact = delti
        else
         elam   = eferm
         delact = delti+ddelta
        endif
c                                         compute actual quasiparticle energies
c                                         and occupation-weights,
c                                         accumulate gap-equation and
c                                         particle number
        parnum = zero
        gapeq  = zero
         do 220 i=1,nmax
         equasi = sqrt((e(i)-elam)*(e(i)-elam)+delact*delact)
         gw(i)  = 0.5-0.5*(e(i)-elam)/equasi
         gapeq  = ph(i)/equasi   + gapeq
         parnum = gw(i)*ph(i)    + parnum
  220    continue
c                                         store actual part.number and
c                                         gap-equation
        dgap(k)  = 0.5*gp*gapeq-one
        dparn(k) = parnum+parnum-nz
  210   continue
c
c       construction of the gradient-matrix and its inverse
c        ( determinant 'det' is limited in value to avoid singularities )
c
        gr(1,1) = ( dgap(2)- dgap(1))/deferm
        gr(1,2) = (dparn(2)-dparn(1))/deferm
        gr(2,1) = ( dgap(3)- dgap(1))/ddelta
        gr(2,2) = (dparn(3)-dparn(1))/ddelta
        det     = gr(1,1)*gr(2,2)-gr(1,2)*gr(2,1)
        det     = sign(max(abs(det),detmin),det)
        if(abs(det).lt.xmaxlw)  det    = sign(xmaxlw,det)
c
        if(iab.ge.1) then
         write(6,'(a,i4,3(a,g11.3)/a,3g11.3,a,3g11.3/a,4g11.3)')
     &    ' iteration at it=',it,' :  delti=',delti,'  eferm=',eferm,
     &    '  det=',det,
     &    ' dgap=',dgap,' ,  dparn=',dparn,
     &    '   gr=',gr
         write(6,'(a,10f7.4,10(/6x,10f7.4))')
     &    '   gw=',gw
        endif
c
        grsav   =  gr(1,1)
        gr(1,1) =  gr(2,2)/det
        gr(1,2) = -gr(1,2)/det
        gr(2,1) = -gr(2,1)/det
        gr(2,2) = grsav/det
c
c       tangential step  (including some step limiters)
c
        efold  = eferm
        delold = delti
        efstep = dampin*(gr(1,1)*dgap(1)+gr(2,1)*dparn(1))
        dlstep = dampin*(gr(1,2)*dgap(1)+gr(2,2)*dparn(1))
        eferm  = eferm-sign(min(abs(efstep),efmax),efstep)
        delti  = abs(delti-sign(min(abs(dlstep),dlmax),dlstep) )
        deferm = efold-eferm
        ddelta = delold-delti
        if(abs(deferm).lt.defmin) deferm = defmin
        if(abs(ddelta).lt.ddlmin) ddelta = ddlmin
c
        if(iab.eq.0)
     &   write(6,'(a,i4,2(a,g11.3))')
     &    ' pairing iteration it=',it,'  :  delti=',delti,
     &    '  eferm=',eferm
c
        if(abs(dgap(1)).lt.eps .and. abs(dparn(1)).lt.eps)   goto 299
  200   continue
        iwarng = 1
c
  299   continue
        delta  = delti
c
      else
c
c-------------------------------------------------------------------------------
c
c      case of constant paring-gap, given by 'delta'.   
c      or case of Temperature, given by Delta
c
c-------------------------------------------------------------------------------
c
       do 310 it=1,abs(iter)
c
c      particle number at 'eferm' (k=1) and 'eferm' varied (k=2).
c
        do 320 ka=1,2
c                                         set actual fermi-energy 'elam'
        k      = 3-ka
        if(k.eq.2) then
         elam   = eferm+deferm
        else
         elam   = eferm
        endif
c                                         compute actual particle number
        parnum = zero
	if(ipair.eq.0) then                          ! pairing
         do 330 i=1,nmax
         equasi = sqrt( (e(i)-elam)*(e(i)-elam) + delta*delta )
         gw(i)  = 0.5-0.5*(e(i)-elam)/equasi
         parnum = gw(i)*ph(i)    + parnum
  330    continue
        elseif(ipair.eq.3) then                      ! thermal
         do 331 i=1,nmax
         equasi = sqrt( (e(i)-elam)*(e(i)-elam) + delta*delta )
         gw(i)  = 1.0/(1.0+exp((e(i)-elam)/delta))
         parnum = gw(i)*ph(i)    + parnum
  331    continue
	endif
c                                         store particle number
        dparn(k) = parnum+parnum-nz
  320   continue
c
c      gradient step to next 'eferm'   (including step limiters)
c
        gr(1,1) = (dparn(2)-dparn(1))/deferm
        det    = sign(max(abs(gr(1,1)),detmin),gr(1,1))
        efold  = eferm
        efstep = dampin*dparn(1)/det
        eferm  = eferm-sign(min(abs(efstep),efmax),efstep)
        deferm = efold-eferm
        deferm = sign(max(abs(deferm),defmin),deferm)
        if(iab.eq.0)
     &   write(6,'(a,i4,2(a,g11.3))')
     &    ' pairing iteration it=',it,'  :  eferm=',eferm,
     &    '  dparn(1)=',dparn(1)
        if(iab.ge.1)
     &   write(6,'(a,i4,2(a,g11.3)/a,2g11.3,a,g11.3/
     &             a,10f7.4,10(/6x,10f7.4))')
     &    ' iteration at it=',it,' :   eferm=',eferm,
     &    '  det=',det,
     &    ' dparn=',dparn(1),dparn(2),' ;   gr=',gr(1,1),
     &    '   gw=',gw
        if(abs(dparn(1)).lt.eps)   goto 350
  310   continue
        iwarng = 1
  350  continue
c
c      determine force 'gp' according to given 'delta' and the
c      equilibrium reached.
c
       gapeq  = zero
        do 360 i=1,nmax
  360   gapeq = ph(i)/sqrt((e(i)-eferm)*(e(i)-eferm)+delta*delta)
     *           + gapeq
       gp     = 2.0/gapeq
c
      endif
c
c-------------------------------------------------------------------------------
c
c     finally the gap-equation and particle-number are checked and printed
c
c-------------------------------------------------------------------------------
c
      if(ipair.eq.0 .or. ipair.eq.1) then
       gphalf = gp*0.5
       del2   = max(xmaxlw,delta*delta)
       gapeq  = zero
       partnm = zero
        do 410 i=1,nmax
        partnm = gw(i)*ph(i) + partnm
  410   gapeq  = gphalf*ph(i)/sqrt(del2+(e(i)-eferm)**2) + gapeq
c
        if(iab.ge.-1)
     &   write(6,'(a,i4,2(a,g13.5))')
     &   ' pair finished with ',it,' iterations: gap-eq.=',gapeq,
     &   '  part.nr=',partnm
c                                  the pairing energy
        epair    = -delta*delta/(gp+xsmall)
      elseif(ipair.eq.3) then
        entrop = zero
        partnm = zero
        do 411 i=1,nmax
        partnm = gw(i)*ph(i) + partnm
	occ    = gw(i)
	occd   = 1.0-occ
	entrop = ph(i)*(occ*log(occ)+occd*log(occd)) + entrop
  411   continue 
c
        if(iab.ge.-1)
     &   write(6,'(a,i4,2(a,g13.5))')
     &   ' pair finished with ',it,' iterations: entropy=',entrop,
     &   '  part.nr=',partnm
c                                 the T*S
	epair   = Delta*entrop
      endif
c
c     in case of incomplete convergence print a warning
c
      if(iwr.ge.-1 .and. iwarng.eq.1)
     & write(6,'(2(a,g11.3),2(a,i3))')
     &  ' ---> no convergence in pair! err(gap-eq)=',dgap(1),
     &  '  err(part.nr)=',dparn(1),'  nz=',nz,'  nmzx=',nmax
c
      return
      end
C-----MATIN2----------------------------------------------------------
C
      SUBROUTINE  MATIN2 (A,NMAX,NDIM,NERROR,DETERM)
C
c      INCLUDE 'PARAMS.INC'
      implicit double precision (a-h,o-z)
      parameter (KPARMX=3)
      PARAMETER (KINDEX=KPARMX*KPARMX)
C
C     INVERTS THE MATRIX 'C'.
C
C     THE PARAMETERS ON THE I/O LIST ARE:
C
C       A       = AT INPUT:  MATRIX TO BE INVERTED
C                 AT OUTPUT: INVERSE MATRIX
C       NMAX    = NUMBER OF ACTUAL ROWS AND COLUMNS
C       NDIM    = DIMENSION OF MATRIX IN THE CALLING ROUTINE
C       NERROR  = ERROR FLAG
C       DETERM  = DETERMINANT OF INPUT MATRIX 'A'
C
C       THE ROUTINE IS TAKEN FROM 'CERNLIB' AND MODIFIED FOR THIS
C       PACKAGE.
C
      DIMENSION A(KINDEX),INDEX(KINDEX)
cold      data r0,r1/0.0,1.0/
C
C-----------------------------------------------------------------------------
C
      IF(NMAX.GT.KINDEX) STOP ' TOO LARGE DIMENSION IN MATIN2'
C
      DETER=1.00
      N=NMAX
C     IEMAT=N+N2
      IEMAT=N+0
      IDIM=NMAX
      NMIN1=N-1
C
C   UMSETZEN DES NMAX*NMAX-ANTEILS DER MATRIX A IN GEPACKTE
C   EIN-DIMENSIONALE FORM.
C
      DO 400 I=2,IDIM
      DO 400 J=1,IDIM
      IPACK  = (I-1)*NMAX + J
      IMATR  = (I-1)*NDIM + J
      A(IPACK) = A(IMATR)
  400 CONTINUE
C        THE ROUTINE DOES ITS OWN EVALUATION FOR DOUBLE SUBSCRIPTING OF
C        ARRAY A.
      IPIVC=1-IDIM
C        MAIN LOOP TO INVERT THE MATRIX
      DO 11 MAIN=1,N
      PIVOT=0.00
      IPIVC=IPIVC+IDIM
C        SEARCH FOR NEXT PIVOT IN COLUMN MAIN.
      IPIVC1=IPIVC+MAIN-1
      IPIVC2=IPIVC +NMIN1
      DO 2 I1=IPIVC1,IPIVC2
      IF(ABS(A(I1))-ABS(PIVOT)) 2,2,1
    1 PIVOT=A(I1)
      LPIV=I1
    2 CONTINUE
C        IS PIVOT DIFFERENT FROM ZERO
      IF(PIVOT) 3,15,3
C        GET THE PIVOT-LINE INDICATOR AND SWAP LINES IF NECESSARY
    3 ICOL=LPIV-IPIVC+1
      INDEX(MAIN)=ICOL
      IF(ICOL-MAIN) 6,6,4
C        COMPLEMENT THE DETERMINANT
    4 DETER=-DETER
C        POINTER TO LINE PIVOT FOUND
      ICOL=ICOL-IDIM
C        POINTER TO EXACT PIVOT LINE
      I3=MAIN-IDIM
      DO 5 I=1,IEMAT
      ICOL=ICOL+IDIM
      I3=I3+IDIM
      SWAP=A(I3)
      A(I3)=A(ICOL)
    5 A(ICOL)=SWAP
C        COMPUTE DETERMINANT
    6 DETER=DETER*PIVOT
      PIVOT=1./PIVOT
C        TRANSFORM PIVOT COLUMN
      I3=IPIVC+NMIN1
      DO 7 I=IPIVC,I3
    7 A(I)=-A(I)*PIVOT
      A(IPIVC1)=PIVOT
C        PIVOT ELEMENT TRANSFORMED
C
C        NOW CONVERT REST OF THE MATRIX
      I1=MAIN-IDIM
C        POINTER TO PIVOT LINE ELEMENTS
      ICOL=1-IDIM
C        GENERAL COLUMN POINTER
      DO 10 I=1,IEMAT
      ICOL=ICOL+IDIM
      I1=I1+IDIM
C        POINTERS MOVED
      IF(I-MAIN) 8,10,8
C        PIVOT COLUMN EXCLUDED
    8 JCOL=ICOL+NMIN1
      SWAP=A(I1)
      I3=IPIVC-1
      DO 9 I2=ICOL,JCOL
      I3=I3+1
    9 A(I2)=A(I2)+SWAP*A(I3)
      A(I1)=SWAP*PIVOT
   10 CONTINUE
   11 CONTINUE
C        NOW REARRANGE THE MATRIX TO GET RIGHT INVERS
      DO 14 I1=1,N
      MAIN=N+1-I1
      LPIV=INDEX(MAIN)
      IF(LPIV-MAIN) 12,14,12
   12 ICOL=(LPIV-1)*IDIM+1
      JCOL=ICOL+NMIN1
      IPIVC=(MAIN-1)*IDIM+1-ICOL
      DO 13 I2=ICOL,JCOL
      I3=I2+IPIVC
      SWAP=A(I2)
      A(I2)=A(I3)
   13 A(I3)=SWAP
   14 CONTINUE
      DETERM=DETER
      NERROR=0
C   RUECKSETZEN DER GEPACKTEN EINDIMENSIONALEN FORM IN DEN
C   NMAX*NMAX-ANTEIL DER MATRIX C
C
      DO 500 I=IDIM,2,-1
      DO 500 J=IDIM,1,-1
      IPACK  = (I-1)*NMAX + J
      IMATR  = (I-1)*NDIM + J
  500 A(IMATR) = A(IPACK)
C
      RETURN
   15 NERROR=MAIN
      DETERM=DETER
      RETURN
      END

