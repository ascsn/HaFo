c-----enrwig--------------------------------------------------------
c
      subroutine enrwig
c
      include 'params.inc'
c
c     Wigner energy in LDA:
c     uses densities given in common /wavfun/
c
c     Comment: The LDA model seems to have a general problem.
c              It produces energies which grow with A
c              while the Wigner energy should shrink as 47MeV/A.
c
      include 'numpar.inc'
      include 'grid.inc'
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'constn.inc'
c
c
c-----------------------------------------------------------------------
c
c     the surface energy by comparison with LDA
c
      sum    = 0.0
      do ir=2,ngrid
        rhop   = rho(ir)
        rhon   = rho(ir+ngrid)
        ract   = (ir-1)*rstep
        sum    = (ract*ract)*funwig(rhop,rhon) + sum
      enddo
      sum     = sum*frpirs
c
      write(11,'(a,1pg12.4)') 
     &  ' the Wigner energy in LDA=',sum
      write(6,'(a,1pg12.4)') 
     &  ' the Wigner energy in LDA=',sum
c
      return
      end
c-----funwig--------------------------------------------------------
c
      function funwig(rhop,rhon)
c
      include 'params.inc'
c
c     the model functional for the Wigner energy
c
      data wig_coupl/1.0D0/     ! the overall coupling strength
      data rhoswitch/0.08D0/    ! switching density
c
c-----------------------------------------------------------------------
c
      rhot = rhop+rhon
      xi   = (rhon-rhop)/rhot
      funwig = wig_coupl*0.5D0
     &        *(1.0D0-abs(xi))**4
     &        *(rhot)*exp(-(rhot/rhoswitch)**2)
c
      return
      end







