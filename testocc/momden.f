c-----momden------------ part of the hartree-fock package ---------------------
c
      subroutine momden
c
      include 'params.inc'
c
c     * * * * * * * * *  density in momentum space   * * * * * * * * * *
c
c     accumulates the momentum density and prints it.
c
c     for this task the routine needs some information on the shell
c     structure; this is provided by the following fields:
c      jsp   = total and orbital angular momentum for the state 'nshl'.
c              (for a more detailed explanation see 'initgr')
c      zwj   = "nr. of states in shell with combined ang.mom. 'js'",
c              i.e. (2*j+1)/4*pi where j is total angular momentum.
c      wls   = strength of spin-orbit splitting for ang.mom. 'js'
c     (note that 'js' and also 'jsp' is a combined quantum number
c      for total and orbital ang.mom.)
c
      include 'skypar.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'constn.inc'
c
c     provide space for accumulation
c
      dimension auxmom(3*kg1)
      equivalence (aux1,auxmom)
c
c      dimension zwj(19)
c      data zwj/.15915494,.15915494,.31830989,.31830989,.47746483,
c     &         .47746483,.63661977,.63661977,.79577472,.79577472,
c     &         .95492966,.95492966,1.1140846,1.1140846,1.2732395,
c     &         1.2732395,1.4232945,1.4232945,1.5915494/
c
c------------------------------------------------------------------------------
c
      open(88,file='momden.res')
      delqmd = 0.02
      nkqmax = 250
      if(nkqmax.gt.3*kg1) stop ' too large nlqmax in momden'
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
      nadd   = ngw*ngrid-ngrid
c
c     zero fields before accumulation
c
      do i=1,nkqmax
        auxmom(i) = zero
      enddo
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
      om     = zwj(js)*weight(nshl+ngw)/17.4532
      nadd   = nadd+ngrid
      if(om.lt.xsmall) goto 30
      if(2*ll.eq.js) js=js-1
      write(88,'(a,4i5/a)')
     & ' momentum density for state iso,nr.state,l,2*j,=',
     &   iso,nshl,ll,js,
     & '     q       density   '
c      sumtst = zero
cc      iqmax  = kqmax
c                                                 q=0.0 separately
      if(ll.eq.0) then
       sum    = zero
       r      = zero
       do ir=2,ngrid
         r     = rstep+r
         sum   = r*wf(nadd+ir)+sum
       enddo
c       sumtst  = sum*sum + sumtst
       write(88,'(1x,f8.3,1pg12.4)') 0.0,sum*sum
       auxmom(1) = om*sum*sum + auxmom(1)
      endif
c                                                      all other q
      do 530 iq=2,nkqmax
      q1     = iq*delqmd-delqmd
c      dc     = cos(q1*rstep)
c      ds     = sin(q1*rstep)
      r      = zero
      sum    = zero
c      c      = one
c      s      = zero
       do 510 i=2,ngrid
       r      = rstep+r
c       x      = c*dc-s*ds
c       s      = s*dc+c*ds
c       c      = x
       sum    = r*bessel(r*q1,ll)*wf(nadd+i)+sum
  510  continue
c      sumtst   = (sum*q1)**2 + sumtst
      write(88,'(1x,f8.3,1pg12.4)') q1,sum*sum
      auxmom(iq) = om*sum*sum + auxmom(iq)
  530 continue
c      write(88,'(a,2i4,a,g12.5)')
c     & ' state iso,na=',iso,nshl,': norm=',sumtst
c
      write(88,'(1x/1x)') 
   30 continue
c
c     now print
c
      write(88,'(/a,i2,a/a,999(/f9.3,g12.4))')
     & ' momentum density for iso=',iso,':',
     & '     q     density',
     & (delqmd*(iq-1),auxmom(iq),iq=1,nkqmax)
c
c     test total particle munber
c
      summom = zero
      do 590 iq=1,nkqmax
      q1     = iq*delqmd-delqmd
      summom = auxmom(iq)*q1*q1 + summom
  590 continue
      write(88,'(a,g12.5)')
     & ' test total density: ',fourpi*delqmd*summom
c
   20 continue
c
c     test bessel
c
c      do 310 ll=0,2
c      do 310 ir=1,ngrid,2
c      r     = rstep*(ir-1)
c      sum   = zero
c      do 320 iq=1,2*nkqmax
c      q1     = iq*delqmd-delqmd
c      sum   = (bessel(r*q1,ll)*q1)**2 + sum
c  320 continue
c      write(88,'(a,i3,2g12.5)')
c     & ' test bessel: ll,r,bessel=',ll,r,fourpi*delqmd*sum*r*r
c  310 continue
c
      return
      end
c-----bessel------------------------------------------------------
c
      function bessel(x,l)
c
      include 'params.inc'
c
c     spherical bessel function of angular momentum 'l' at
c     argument 'x'.
c
c     uses series expansion for x<l=1 and asymptotic expression
c     for x>l-1.
c
c     precision is tested only up to 5 digits.
c
      include 'constn.inc'
c
      data xargmx/1.0D10/
c
c------------------------------------------------------------------------------
c
cc       if(l.gt.7) stop ' l>7 in "bessel"'
c
c      case of zero argument
c
       if(x.eq.zero) then
        if(l.eq.0) then
         b      = one
        else
         b      = zero
        endif
c
c      very large argument  --> set result to zero
c
       elseif(abs(x).gt.xargmx) then
         b      = zero
c
c      small arguments  --> use series expansion
c
       elseif(abs(x).le.l-one) then
        add   = one
        do i=1,l
          add   = x/(2*i+1)*add
        enddo
        arg2  = -0.25*x*x
        b     = add
        do i=1,l+1
          add   = arg2*add/(i*(l+i+0.5D0))
          b     = add + b
        enddo
c
c      large arguments  --> use expansion in sin and cos
c
       else
        s      = sin(x)/x
        c      = cos(x)/x
        if(l.eq.0) then
         b      = s
        elseif(l.eq.1) then
         b      = s/x-c
        elseif(l.ge.2) then
         xinv   = one/x
         bjm    = s
         b      = s/x-c
         do lrec=2,l
           bjp    = ((lrec+lrec)*xinv-xinv)*b-bjm
           bjm    = b
           b      = bjp
         enddo
         b      = b
        endif
       endif
c
       bessel = b
       return
       end
