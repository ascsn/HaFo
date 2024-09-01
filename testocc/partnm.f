c-----nmpart--------------------------------------------------------
c
      subroutine nmpart(rhonm)
c
      include 'params.inc'
c
c     separate pieces of E/A in symmetric nuclear matter
c      rhonm  = input: total density
c     The force paxrameters are entered via common /skypar/
c
      include 'skypar.inc'
      include 'wavfun.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
      data betab /5.742467999/     ! (3/5)*(3*pi^2)^*(2/3)  factor for tau_n/p
      data f2d3  /0.666666667/
c
c-----------------------------------------------------------------------
c
      rhopow = rhonm**power *half
      rhop   = procon*rhonm
      taup   = betab*rhop**f2d3*rhop
      taunm  = taup+taup
      if(rhoswm.gt.zero) then
        tanfac = tanh(rhonm/rhoswm)
      else
        tanfac = 1.0
      endif
      ekinA  = (amfc(1)+amfc(2))*taup/rhonm
      epotA  = half*(t0a-half*t0b)*rhonm
      epowA  = rhopow*(t3a-half*t3b)*rhonm
      emassA = (t1p2-half*t1m2)*taunm*tanfac
      if(exmass.ne.zero) then
         epotpA = 2.0d0*exfunc(rhop)/rhonm
      else
         epotpA = half*(-half*t0b)*rhonm
      endif
      if(ioptsk.eq.2 .or. ioptsk.eq.3)  epowA = dirfun(rhonm)/rhonm
c
      write(02,'(a/a/5f8.3)')
     &   ' parts of E/A:',
     &   '    ekin     epot     epower    emasseff      epotiso     ',
     &   ekinA,epotA,epowA,emassA,epotpA

c
      return
      end




