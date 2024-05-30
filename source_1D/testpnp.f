c-----testpnp-----------------------------------------------------------
c
      subroutine testpnp
c
c     to investigate particle-number projection
c
c
      include 'params.inc'
      include 'grid.inc'
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
      logical twrit17                     ! flag to write details to unit 17
      parameter (twrit17=.true.)
c
      complex*16 rotrho(KGR2),rotchi(KGR2)
      complex*16 gaugephas,accrho,accchi,accnorm,addnorm,rotnorm
      complex*16 addpn,addrho,addchi,addrho3,accrho3,projphas
      complex*16 dcrit,rotneut,projnum,cphasrot
      external rotnorm
c
      data maxphi/71/
      data minshift,maxshift/-60,300/
      data ioffneut/0/                    ! offset presently 0
c
c-----------------------------------------------------------------------
c
      ioffset= ngrid*ioffneut
      delphi = 2.0D0*PI/maxphi
      frpirs = 4.0D0*PI*rstep
c      power  = 7.0D0/3.0D0
      sphas1 = 2.0D0*PI/3.0D0
      sphas2 = 4.0D0*PI/3.0D0
      cphasrot =(1.0D0+EXP(DCMPLX(0.0D0,sphas1))
     &           +EXP(DCMPLX(0.0D0,sphas2))      )/3.0D0
      write(16,'(/a,i5)') 
     &   ' PNP projection for several shifts: nr. points=',maxphi
      write(16,'(/a)') 
     &   ' ibeta beta  aver.N    rho**2  rho**7/3  chi**2 '
      wcrit = weight(nmaxp+7)
      bcrit = -log(wcrit/(1.0D0-wcrit))/2.0D0
c
c     write header for protocol
c     
      if(twrit17) then
        open(17,form='unformatted',file='pnp.data')
        write(17) ngrid,frpirs
        write(17) ((ir-1)*rstep,ir=1,ngrid)
        write(17) (delphi*iphi-PI,iphi=0,maxphi-1)
        write(17) (bcrit+0.5D-1*ishift,ishift=minshift,maxshift)
      endif
c      write(16,'(a,2(1pg14.7),i5)')  ' wcrit etc:',wcrit,bcrit,nmaxp+7
      do ishift=minshift,maxshift
        realshift = bcrit+0.5D-1*ishift
        accnorm = DCMPLX(0.0d0,0.0d0)
        accrho  = DCMPLX(0.0d0,0.0d0)
        accrho3 = DCMPLX(0.0d0,0.0d0)
        accchi  = DCMPLX(0.0d0,0.0d0)
        accpn   = DCMPLX(0.0d0,0.0d0)
        projnum = DCMPLX(0.0d0,0.0d0)
        testnum = DCMPLX(0.0d0,0.0d0)
        do iphi=0,maxphi-1
          phi       = delphi*iphi-PI
          gaugephas = EXP(DCMPLX(2.0d0*realshift,2.0d0*phi))
          projphas = EXP(DCMPLX(-aneut*realshift,-aneut*phi))
          addnorm = rotnorm(gaugephas)*projphas
c          write(6,'(5x,4(1pg12.4))') gaugephas,addnorm
          call rotdens(rotrho,rotchi,gaugephas)
          call rotnumb(rotneut,gaugephas)
          if(twrit17) 
     &       write(17) addnorm,(rotrho(ir+ioffset),ir=1,ngrid)
          if(iphi.eq.0) actnum = rotneut
c          write(6,'(a,2(1pg12.4))') ' gaugephas=',gaugephas
c          write(6,'(100(/i5,4(1pg12.4)))') 
c     &       (ir,rotrho(ir),rotchi(ir),ir=2,ngrid)
          addrho  = DCMPLX(0.0d0,0.0d0)
          addrho3 = DCMPLX(0.0d0,0.0d0)
          addchi  = DCMPLX(0.0d0,0.0d0)
          addpn   = DCMPLX(0.0d0,0.0d0)
          do ir=2,ngrid
            irc     = ir+ioffset
            ract2   = (ir-1)*rstep
            ract2   = ract2*ract2
            addpn   = addpn  + ract2*rotrho(irc)
            addrho  = addrho + ract2*rotrho(irc)**2
            if(abs(rotrho(ir)).gt.1.0D-99) then
c              power = 2.333333D0
              rr = rotrho(ir)
              ri = dimag(rotrho(irc))
              rphase = atan2(ri,rr)
              rabsv  = rr*rr+ri*ri
              rphase = rphase*power
              rabsv  = rabsv**(power*0.5D0)
              addrho3 = addrho3 + ract2*rabsv*
     &          (exp(DCMPLX(0.0D0,rphase))    )
c     &           +exp(DCMPLX(0.0D0,rphase+sphas1))
c     &           +exp(DCMPLX(0.0D0,rphase+sphas2)))
c              addrho3 = addrho3 + ract2*exp(log(rotrho(ir))*power)
c     &       addrho3 = addrho3 + ract2*rotrho(ir)**2.333333d0
            endif
            addchi  = addchi + ract2*rotchi(irc)**2
c            if(iphi.eq.0) write(6,'(i5,4(1pg12.4))')
c     &       ir,rotrho(ir),rotchi(ir)
          enddo
          projnum = projnum + addnorm*rotneut
          testnum = testnum + addnorm*(rotneut)**power*cphasrot
          accpn   = accpn   + addnorm*addpn*frpirs
          accnorm = accnorm + addnorm
          accrho  = accrho  + addnorm*addrho*frpirs
          accrho3 = accrho3 + addnorm*addrho3*frpirs
          accchi  = accchi  + addnorm*addchi*frpirs/gaugephas
c          if(iphi.eq.0) 
c     &      write(6,'(5x,4(1pg12.4))') addrho*frpirs,addchi*frpirs
        enddo
        gaugephas = EXP(DCMPLX(2.0d0*realshift,2.0d0*delphi*(maxphi/2)))
        dcrit    = 1.0D0-wcrit+wcrit*(gaugephas)
c        dcrit    = 1.0D0-wcrit-wcrit*exp(2.0d0*realshift)
        write(16,'(1x,i6,f11.7,8(1pg18.10))')
     &    ishift,realshift,actnum,real(testnum/accnorm),
     &    real(accrho/accnorm),
     &    real(accchi/accnorm),real(accrho3/accnorm)
c     &    ishift,realshift,cdabs(dcrit),real(accrho/accnorm),
c     &    real(accchi/accnorm),real(accrho3/accnorm)
c     &    realshift,real(accpn/accnorm),real(accrho/accnorm),
c     &    real(accchi/accnorm),real(accrho3/accnorm)
      enddo
c
      return
      end


c-----rotnumb---------- part of the hartree-fock package ----------------------
c
      subroutine rotnumb(rotneut,gaugephas)
c
      include 'params.inc'
c
c     the gauge-angle 'rotgauge' rotated neutron number 'rotneut'
c
      include 'wavfun.inc'
      include 'result.inc'
      include 'grid.inc'
c
      complex*16 rotneut
      complex*16 denom,gaugephas,rhofac,chifac
c
c------------------------------------------------------------------------------
c
c
c     loop over isospin
c
      rotneut = DCMPLX(0.0D0,0.0D0)
      phasej = zwj(1)*0.5
      do iso=2,2
        isadd  = ngrid*nmaxp*(iso-1)
        if(iso.eq.1) then
         nend   = nmaxp
         isoadd = 0
        else
         nend   = nmaxn
         isoadd = nmaxp
        endif
c
c       loop over states in oscillator ordering
c
        do n=1,nend
          js     = jsp(n)
          ll     = js/2
          j2     = ((js+1)/2)*2-1
          ie     = n*ngrid+isadd
          ia     = ie-ngrid+1
          occfac = zwj(js)/phasej
          vfac   = weight(isoadd+n)
          ufac   = sqrt(max(1.0d0-vfac,xsmall))
          vfac   = sqrt(max(vfac,xsmall))
          denom  = 1.0d0/(ufac*ufac+vfac*vfac*gaugephas)
          rhofac = occfac*vfac*vfac*gaugephas*denom
          rotneut = rhofac + rotneut
        enddo
c
      enddo
c
      return
      end

c-----rotdens---------- part of the hartree-fock package ----------------------
c
      subroutine rotdens(rotrho,rotchi,gaugephas)
c
      include 'params.inc'
c
c     the gauge-angle 'rotgauge' rotated density 'rotrho'
c     and pair density 'rotchi', only for neutrons .
c
      include 'wavfun.inc'
      include 'result.inc'
      include 'grid.inc'
c
      complex*16 rotrho(KGR2),rotchi(KGR2)
      complex*16 denom,gaugephas,rhofac,chifac
c
c------------------------------------------------------------------------------
c
      do ir=1,ngrid
        rotrho(ir) = DCMPLX(0.0d0,0.0d0)
        rotchi(ir) = DCMPLX(0.0d0,0.0d0)
      enddo
c
c     loop over isospin
c
      phasej = zwj(1)*0.5
      do iso=2,2
        isadd  = ngrid*nmaxp*(iso-1)
        if(iso.eq.1) then
         nend   = nmaxp
         isoadd = 0
        else
         nend   = nmaxn
         isoadd = nmaxp
        endif
c
c       loop over states in oscillator ordering
c
        do n=1,nend
          js     = jsp(n)
          ll     = js/2
          j2     = ((js+1)/2)*2-1
          ie     = n*ngrid+isadd
          ia     = ie-ngrid+1
          occfac = zwj(js)/phasej
          vfac   = weight(isoadd+n)
          ufac   = sqrt(max(1.0d0-vfac,xsmall))
          vfac   = sqrt(max(vfac,xsmall))
          denom  = 1.0d0/(ufac*ufac+vfac*vfac*gaugephas)
          rhofac = occfac*vfac*vfac*gaugephas*denom
          chifac = occfac*ufac*vfac*gaugephas*denom
c          write(6,'(1x,a,4(1pg12.4))')
c     &     ' phasej,occfac,rhofac=',phasej,occfac,rhofac
c
          if(ll.eq.0) then
            rhosp  = wfs(ia)**2*phasej
          else
            rhosp  = 0.0
          endif
          rotrho(1) = rotrho(1) + rhosp*rhofac
          rotchi(1) = rotchi(1) + rhosp*chifac
c
          do i=ia+1,ie
            r     = (i-ia)*rstep
            rhosp = (wf(i)/r)**2*phasej
            irho  = i-ia
            rotrho(irho) = rotrho(irho) + rhosp*rhofac
            rotchi(irho) = rotchi(irho) + rhosp*chifac
          enddo
c
        enddo
c
      enddo
c
      return
      end


c-----rotnorm---------- part of the hartree-fock package ----------------------
c
      complex*16 function rotnorm(gaugephas)
c
      include 'params.inc'
c
c     the gauge-angle 'rotgauge' rotated norm overlap, only for neutrons .
c
      include 'wavfun.inc'
      include 'result.inc'
      include 'grid.inc'
c
      complex*16 prod,gaugephas
c
c------------------------------------------------------------------------------
c
c     loop over isospin
c
      prod   = DCMPLX(1.0d0,0.0d0)
      phasej = zwj(1)*0.5
      do iso=2,2
        isadd  = ngrid*nmaxp*(iso-1)
        if(iso.eq.1) then
         nend   = nmaxp
         isoadd = 0
        else
         nend   = nmaxn
         isoadd = nmaxp
        endif
c
c       loop over states in oscillator ordering
c
        do n=1,nend
          js     = jsp(n)
          ll     = js/2
          if(2*ll.eq.js) js=js-1
          mmult = (js+1)/2                   ! half the multiplicity
          j2     = ((js+1)/2)*2-1
          ie     = n*ngrid+isadd
          ia     = ie-ngrid+1
          occfac = weight(isoadd+n)*zwj(js)/phasej
          v2     = weight(isoadd+n)
          u2     = 1.0d0-v2
          prod   = prod*(u2+v2*gaugephas)**mmult        
c
        enddo
c
      enddo
      rotnorm = prod
c
      return
      end




