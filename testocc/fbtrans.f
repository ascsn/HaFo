c-----fbtpot--------------------------------------------------------------
c
      subroutine fbtpot
c
c     computes and prints Fourier-Bessel transformed local potentials
c
      include 'params.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'result.inc'
      include 'potinv.inc'
      include 'constn.inc'
c
c     parameters of the Woods-Saxon potential
c
      data wsrad,wsv0,wskap,wswid/1.24,53.3,0.63,0.63/
c
c------------------------------------------------------------------------
c
      radact = wsrad*anumb**0.3333333
c
c     prepare Coulomb potential of homogenesouly charged sphere
c
      chgtot = 0.0
      do i=1,ngrid
        ract = (i-1)*rstep
        if(ract.le.radact) chgtot = rsquar(i) + chgtot
      enddo
      chgtot = frpirs*chgtot
      do i=1,ngrid
        ract = (i-1)*rstep
        if(ract.le.radact) then
          rhochr(i) = (aprot-1.0)/chgtot
        else
          rhochr(i) = 0.0
        endif
      enddo
      call coulmb(rhochr(1),binv,1.0D0)
c
c     compose WS potentials on 'aux1'
c
      wsvpr = -wsv0*(1.0+wskap*(aneut-aprot)/anumb)
      wsvne = -wsv0*(1.0-wskap*(aneut-aprot)/anumb)
      do i=1,ngrid
        ract = (i-1)*rstep
        wsfac = 1.0/(1.0+exp((ract-radact)/wswid))
        aux1(i) = wsvpr*wsfac+binv(i)
        aux1(i+ngrid) = wsvne*wsfac
      enddo
c
c     print in coordinate space
c
      write(11,'(1x/a/a)')
     &  ' Print potentials:',
     &  '   r     U_prot   U_neut     WS_prot    WS_neut  '
      do i=2,ngrid
        ract = (i-1)*rstep
        write(11,'(1x,f8.4,5(1pg12.4))') 
     &   ract,pot(i),pot(i+ngrid),aux1(i),aux1(i+ngrid),binv(i)
      enddo
c
c     loop over 'q', compute Fourier-Bessel transformations, print
c
      write(11,'(1x/a/a)')
     &  ' Print Fourier-Bessel transformed potentials:',
     &  '   q     U_prot   U_neut     WS_prot    WS_neut  '
      sumup   = zero
      sumun   = zero
      sumwsp   = zero
      sumwsn   = zero
      do i=1,ngrid
        ract = (i-1)*rstep
        sumup  = rsquar(i)*pot(i) + sumup
        sumun  = rsquar(i)*pot(i+ngrid) + sumun
        sumwsp = rsquar(i)*aux1(i) + sumwsp
        sumwsn = rsquar(i)*aux1(i+ngrid) + sumwsn
      enddo
      sumup   = frpirs*sumup
      sumun   = frpirs*sumun
      sumwsp  = frpirs*sumwsp
      sumwsn  = frpirs*sumwsn
      write(11,'(1x,f8.4,4(1pg12.4))') 
     &   0.0D0,abs(sumup),abs(sumun),abs(sumwsp),abs(sumwsn)
c
      do iq=2,iqmax
        q1     = iq*delq-delq
c
        dc     = cos(q1*rstep)
        ds     = sin(q1*rstep)
        r      = zero
        sumup   = zero
        sumun   = zero
        sumwsp   = zero
        sumwsn   = zero
        c      = one
        s      = zero
        do i=2,ngrid
          r      = rstep+r
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
          sumup  = r*s*pot(i)+sumup
          sumun  = r*s*pot(i+ngrid)+sumun
          sumwsp = r*s*aux1(i)+sumwsp
          sumwsn = r*s*aux1(i+ngrid)+sumwsn
        enddo
        sumup   = frpirs*sumup/q1
        sumun   = frpirs*sumun/q1
        sumwsp  = frpirs*sumwsp/q1
        sumwsn  = frpirs*sumwsn/q1
c
        write(11,'(1x,f8.4,4(1pg12.4))') 
     &   q1,abs(sumup),abs(sumun),abs(sumwsp),abs(sumwsn)
c
      enddo
c
      return
      end
