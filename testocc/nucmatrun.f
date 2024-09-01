c-----nucmutrun---------------------------------------------------------------
c
      program nucmatrun
c
      include 'params.inc'
c
c     routine to compute nuclear matter properties for Skyrme forces
c     force parameters are set in blockdata
c
      include 'skypar.inc'
      include 'wavfun.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
c
      call initfc
ctest      write(*,*) ' amfc=',amfc
      do irho=1,20
        rhonm = irho*0.01
        ratiopn = 0.5D0
        call enernm(rhonm,ratiopn,0.0D0,energm,efmasp,efmasn)
        write(*,*) ' rhonm,ratiopn,energm=',rhonm,ratiopn,energm
      enddo
      stop
      end
c-----block ------------ part of the hartree-fock package ---------------------
c
      blockdata  block
c
      include 'params.inc'
c
c     settings for SLy6
c
      include 'skypar.inc'
      include 'wavfun.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
      data  T0,T1,T2,T3,T4
     &     /-2479.500,   462.180,  -448.610, 13673.000,   122.000/
      data  T5,T6,T7,T8,T9,TYUK/6*0.0/
      data  X0,X1,X2,X3,b4p
     &     / 0.82500,  -0.46500,  -1.00000,   1.35500,  61.00000/
      data  X4,X5,X6,X7,X8,XYUK, AYUK,rhoswm,exmass/9*0.0D0/ 
      data  POWER/0.166666667/
      data  ioptsk,ifdfms,ifzpe,ifex
     &      /0,0,1,0/
c
      end
c-----enernm--------------------------------------------------------
c
      subroutine enernm(rhonm,procon,temper,energm,efmasp,efmasn)
c
      include 'params.inc'
c
c     E/A for nuclear matter
c      rhonm  = input: total density
c      procon = input: proton content (0.5 for symm.matter)
c      temper = temperature
c      energm = output: E/A in MeV
c      effmas = output: effective mass m*/m
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
      rhop   = procon*rhonm
      rhon   = rhonm-rhop
      if(temper.eq.zero) then
        taup   = betab*rhop**f2d3*rhop
        taun   = betab*rhon**f2d3*rhon
      elseif(procon.eq.half) then
        call hotnm(rhop,taup,temper)
        rhon   = rhop
        taun   = taup
      else
        stop ' non zero T only for symmetric nuclear matter '
      endif
      zpefac = one
      call skfunc(rhop,rhon,zero,zero,zero,zero,zero,
     &              taup,taun,zero,zero,zero,zero,
     &              zero,zero,
     &              potp,potn,bmasp,bmasn,potlsp,potlsn,
     &              energd)
      enrg   = zpefac*(amfc(1)*taup+amfc(2)*taun) + energd
c        if(iwr.gt.-2)
c     &   write(6,'(a,g12.4)') ' ekin=',(amfc(1)*taup+amfc(2)*taun)/rhonm
ctest      write(*,*) ' energd,zpefac,amfc,rhonm=',
ctest     &             energd,zpefac,amfc,rhonm
      energm = enrg/rhonm
      efmasp = zpefac*amfc(1)/(zpefac*amfc(1)+bmasp)
      efmasn = zpefac*amfc(2)/(zpefac*amfc(1)+bmasn)
c
      return
      end

c-----hotnm --------------------------------------------------------
c
      subroutine hotnm(rhonmp,taunmp,temper)
c
      include 'params.inc'
c
c     solves Fermi equilibrium for hot nuclear matter:
c      rhonmp = input: total proton density
c      taunmp = kinetic energy density
c      temper = temperature
c
      include 'skypar.inc'
      include 'constn.inc'
c
      data f1d3  /0.333333333/
      data itmax/50/
      data delfp/0.01/
      data epsfp/1.0E-5/
      data ipgrid/200/
      data delp/0.01/
c
c-----------------------------------------------------------------------
c
      fermp  = (3.*pi*pi*rhonmp)**f1d3            ! first guess Fermi momentum
      amfca0 = 0.5*(amfc(1)+amfc(2))              ! average hbar^2/2m
      do iter=1,itmax
        sum0   = 0.0
        sum1   = 0.0
        fermp2 = fermp*fermp
        fermpp = (fermp+delfp)*(fermp+delfp)
        do ip=1,ipgrid
          pact2  = ip*delp
          pact2  = pact2*pact2
c          write(6,'(a,i5,4(1pg12.4))')
c     &    ' ip,pact,amfca0,frmp2,sum0=',ip,pact,amfca0,frmp2,sum0
          sum0   = pact2/(1.0+exp(amfca0*(pact2-fermp2)/temper)) + sum0
          sum1   = pact2/(1.0+exp(amfca0*(pact2-fermpp)/temper)) + sum1
        enddo
        sum0   = delp*sum0/(pi*pi)
        sum1   = delp*sum1/(pi*pi)
        drhodf = (sum1-sum0)/delfp
        fpstep = -(sum0-rhonmp)/sign(max(abs(drhodf),xsmall),drhodf)
        fpstep = sign(min(abs(fpstep),10.0*delfp),fpstep)
        fermp  = fermp+fpstep
c        write(6,'(a,i5,3(1pg12.4))')
c     &   ' HOTNM: iter,fermp,fpstep,rhonmp=',iter,fermp,fpstep,sum0
        if(abs(fpstep).lt.epsfp) goto 99
      enddo
c
c     Fermi momentum is determined. Now compute kinetic energy.
   99 continue
      sumk   = 0.0
      sum0   = 0.0
      fermp2 = fermp*fermp
      do ip=1,ipgrid
        pact2  = ip*delp
        pact2  = pact2*pact2
        vol    = pact2/(1.0+exp(amfca0*(pact2-fermp2)/temper))
        sum0   = vol       + sum0
        sumk   = vol*pact2 + sumk
      enddo
      sum0   = delp*sum0/(pi*pi)
      sumk   = delp*sumk/(pi*pi)
      taunmp = sumk
      averen = sumk/sum0*amfca0
      fermp  = (3.*pi*pi*rhonmp)**f1d3            ! first guess Fermi momentum
      avere0  = 0.6*fermp*fermp*amfca0
      write(6,'(a,5(1pg12.4))')
     &   ' HOTNM: rhonmp(1..2),taunmp,averen,avere0=',
     &    rhonmp,sum0,sumk,averen,avere0
      write(11,'(a,5(1pg12.4))')
     &   ' HOTNM: rhonmp(1..2),taunmp,averen,avere0=',
     &    rhonmp,sum0,sumk,averen,avere0
c
      return
      end
c-----skfunc------------ part of the hartree-fock package ---------------------
c
      subroutine skfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,sp,sn,ssp,ssn,
     &                  frp,frn,potp,potn,bmasp,bmasn,wlsp,wlsn,energd)
c
      include 'params.inc'
c
c     the basic Skyrme energy density and derived quantities.
c
c     Input per list:
c      rp,rn        = densities
c      rc           = charge density
c      d1p,d1n      = derivative of densities
c      d2p,d2n      = Laplacian of densities
c      tp,tn        = kinetic energy densities
c      sp,sn        = spin orbit densities
c      ssp,ssn      = derivative of spin-orbit densities
c      frp,frn      = folded densities (for Yukawa term)
c     Forces are input per common /skypar/ und /tswitc/.
c     Ouput per list:
c      energd       = energy density
c      potp,potn    = local potentials
c      bmasp,bmasn  = inverse mass
c      wlsp,wlsn    = spin-orbit potential
c
c     !!! Note: the option 'texten' needs yet careful checking
c
      include 'skypar.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
c     total densities and other auxiliary quantities
c
      rs  = rp+rn
      ts  = tp+tn
      ss  = sp+sn
      sss = ssp+ssn
      d2s = d2p+d2n
      d1s = d1p+d1n
      rspow  = rs**power
      if(rhoswm.gt.zero) then
        tanfac = tanh(rs/rhoswm)
      else
        tanfac = 1.0
      endif
      tyab   = tya-tyb
      t8zwh  = half*t8zw
      riftm  = one*iftm                                          !one*iftm
      excfac = -1.418033*acoul
c
c     basic (density dependent) part of Skyrme functional 
c
      potp = ( t0a*rs-t0b*rp  
     &       +rspow*(t3a*rs-t3b*rp +t3as*rs
     &               -t3bs*(rp*rp+rn*rn)/max(rs,xsmall) ) )
      potn = ( t0a*rs-t0b*rn  
     &       +rspow*(t3a*rs-t3b*rn +t3as*rs
     &               -t3bs*(rp*rp+rn*rn)/max(rs,xsmall) ) )
      bmasp = zero
      bmasn = zero
      wlsp  = zero
      wlsn  = zero
      energd =  +half*(t0a*(rs*rs)-t0b*(rp*rp+rn*rn))
     &          +half*rspow*(t3a*(rs*rs)-t3b*(rp*rp+rn*rn))
c      if(rs.eq.0.14d0) then
c        write(6,'(a,4g12.4)') ' skfunc:',t0a,t0b,t3a,t3b
c      endif
c
c     variants of pure density dependence
c
      if(texcha) then
        add    = excfac*rc**0.3333333
        potp   = add + potp
        energd = 0.75*add*rc + energd
      endif
      if(texmas) then
        potp  = expot(rp) + potp
        potn  = expot(rn) + potn
        energd = exfunc(rp) + exfunc(rn) + energd
      elseif(tdirpa) then
        add  = dirpot(rs)
        potp = add + potp
        potn = add + potn
        energd = dirfun(rs) + energd
      elseif(b5p.ne.zero) then
        potp   = 2.66666667*b5p*rp**1.6666667 + potp
        potn   = 2.66666667*b5p*rn**1.6666667 + potn
        energd = b5p*(rp**2.6666667+rn**2.6666667) + energd
      elseif(texten) then
        potp = ((t7a+t7a)*rs*rs*rs-t7b*rs*(rs*rp+rp*rp+rn*rn)) + potp
        potn = ((t7a+t7a)*rs*rs*rs-t7b*rs*(rs*rn+rp*rp+rn*rn)) + potn
        if(t7.ne.zero) stop ' option T7.NE.0 not active '
      endif
      if(tyukaw) then
       frs    = frp+frn
       potp   = (tya*frs-tyb*frp) + potp
       potn   = (tya*frs-tyb*frn) + potn
       energd = half*(rs*tya*frs+rp*tyb*frp+rn*tyb*frn) + energd
      endif

c
c     standard kinetic and gradient terms
c
      potp   = (t1p2*ts-t1m2*tp)*tanfac-t13m2*d2s+t13p2*d2p + potp
      potn   = (t1p2*ts-t1m2*tn)*tanfac-t13m2*d2s+t13p2*d2n + potn
c      write(6,'(7g12.4)') 
c     &  d2s,d2p,d2n,t1p2,t1m2,t13m2,t13p2
      bmasp  = (t1p2*rs-t1m2*rp)*tanfac  + bmasp
      bmasn  = (t1p2*rs-t1m2*rn)*tanfac  + bmasn
      wlsp   = riftm*(t1m20*sp-t1p2x*ss) + wlsp
      wlsn   = riftm*(t1m20*sn-t1p2x*ss) + wlsn
      energd = (t1p2*rs*ts-t1m2*(rp*tp+rn*tn))*tanfac 
     &        + half*riftm*(t1m20*(sp*sp+sn*sn)-t1p2x*ss*ss)
     &        + half*(-t13m2*d2s*rs + t13p2*(d2p*rp+d2n*rn))
     &        + energd
c
c
c     extended kinetic terms
c
      if(texten) then
        potp  =  +t5p6*(rs+rs)*ts-t5m6*(rs*tp+rp*tp+rn*tn)
     &           -t53m6d*((rs+rs)*d2s+d1s*d1s)
     &           +t53p6*(d2s*rp+(rs+rs)*d2p+rp*d2p+rn*d2n
     &                                      +(d1s+d1s)*d1p)
     &           +t5m60h*(sp*sp+sn*sn)-t5p6xh*ss*ss        + potp
        potn  =  +t5p6*(rs+rs)*ts-t5m6*(rs*tp+rp*tp+rn*tn)
     &           -t53m6d*((rs+rs)*d2s+d1s*d1s)
     &           +t53p6*(d2s*rn+(rs+rs)*d2n+rp*d2p+rn*d2n
     &                                      +(d1s+d1s)*d1n)
     &           +t5m60h*(sp*sp+sn*sn)-t5p6xh*ss*ss        + potn
        bmasp  = rs*(t5p6*rs-t5m6*rp) + bmasp
        bmasn  = rs*(t5p6*rs-t5m6*rn) + bmasn
        wlsp   = riftm*rs*(t5m60*sp-t5p6x*ss) + wlsp
        wlsn   = riftm*rs*(t5m60*sp-t5p6x*ss) + wlsn
        energd = +tp*rs*(t5p6*rs-t5m6*rp)+tn*rs*(t5p6*rs-t5m6*rn)
     &           +rs*t5m60h*(sp*sp+sn*sn)-t5p6xh*ss*ss
     &           +energd
      elseif(rhoswm.gt.zero) then
       add    = (t1p2*rs*ts-t1m2*(rp*tp+rn*tn))
     &          *(1.0-tanfac*tanfac)/rhoswm     
       potp   = add + potp
       potn   = add + potn
c              energy is already contained in 'tanfac'-term above       
      elseif(b6.ne.zero) then
        potp    = half*(b6*ts*ts+b6p*tp*tp) + potp
        potn    = half*(b6*ts*ts+b6p*tn*tn) + potn
        bmasp   = b6*ts*rs+b6p*tp*rp + bmasp
        bmasn   = b6*ts*rs+b6p*tn*rn + bmasn
        energd  = half*b6*rs*ts*ts
     &           +half*b6p*(rp*tp*tp+rn*tn*tn)   + energd
      endif
      if (tssmas) then    
        stop ' option TSSMAS not ective '
c          potp  = -t9*d2ts + potp        ! need d2ts
c          potn  = -t9*d2ts + potn
c          bmasp = -t9*d2s +  bmasp
c          bmasn = -t9*d2s +  bmasn
      endif                 

c
c     standard spin-orbit
c      
      potp   = -b4*sss - b4p*ssp   + potp
      potn   = -b4*sss - b4p*ssn   + potn
      wlsp   = (b4*d1s+b4p*d1p+riftm*(t1m20*sp-t1p2x*ss)) + wlsp
      wlsn   = (b4*d1s+b4p*d1n+riftm*(t1m20*sn-t1p2x*ss)) + wlsn
      energd = -b4*rs*sss - b4p*(rp*ssp+rn*ssn) + energd
c
c     extended spin-orbit terms
c
      if (texten) then
        potp   = -t8zw*(rs*(sss+sss+ssp)+rp*ssp+rn*ssn)  + potp
        potn   = -t8zw*(rs*(sss+sss+ssn)+rp*ssp+rn*ssn)  + potn
        wlsp   =  t8zw*(rs*(d1s+d1s+d1p)+rp*d1s)  + wlsp
        wlsn   =  t8zw*(rs*(d1s+d1s+d1n)+rn*d1s)  + wlsn
        energd = -t8zwh*( rp*(rs*(sss+sss+ssp)+rp*ssp+rn*ssn)
     &                   +rn*(rs*(sss+sss+ssn)+rp*ssp+rn*ssn))
     &           + energd
      endif

c
c     ioptsk.eq.5 --> rho-dependent kinetic and l*s, direct version
c

      if (ioptsk.eq.5) then
        potp = potp + ( -2*b4*c4*rs*sss -b4p*c4p*(rp*ssp+rn*ssn 
     &              + rs*ssp))
     &              + (2*t1p2*c1*rs*ts-b1p*c1p*(rp*tp+rn*tn +c1p*rs*tp))
        potn = potn + ( -2*b4*c4*rs*sss -b4p*c4p*(rp*ssp+rn*ssn 
     &              + rs*ssn))
     &              + (2*t1p2*c1*rs*ts-b1p*c1p*(rp*tp+rn*tn +c1p*rs*tn))
       bmasp = bmasp + t1p2*c1*rs*rs -t1m2*c1p*rs*rp
       bmasn = bmasn + t1p2*c1*rs*rs -t1m2*c1p*rs*rn
       wlsp   = wlsp + (b4*c4*(rs+rs)*d1s+b4p*c4p*(rs*d1p+d1s*rp))
       wlsn   = wlsn + (b4*c4*(rs+rs)*d1s+b4p*c4p*(rs*d1n+d1s*rn))
      endif
c
c     ioptsk.eq.6 --> rho-dependent kinetic and l*s, inverse version
c


c     b1-, b1p-, b4- and b4p-terms that appeared before must be
c     substracted to not add them twice; same procedure with wlsp, wlsn, bmasp and bmasn

      if (ioptsk.eq.6) then

        potp = potp + b4*sss + b4p*ssp - b4*(1/(1-c4*rs)
     &              + c4*rs/(1-c4*rs)*(1-c4*rs))*sss - b4p*(1/(1-c4*rs)
     &              + c4p*(rp+rn)/(1-c4*rs)*(1-c4*rs))*ssp
     &              - t1p2*ts + t1m2*tp +t1p2*(ts/(1-c1*rs)
     &              + c1*rs*ts/(1-c1*rs)*(1-c1*rs))
     &              - t1m2*(c1p*(rp*tp+rn*tn)/(1-c1p*rs)*(1-c1p*rs)
     &              + tp/(1-c1p*rs))
        potn = potn + b4*sss + b4p*ssn - b4*(1/(1-c4*rs)
     &              + c4*rs/(1-c4*rs)*(1-c4*rs))*sss - b4p*(1/(1-c4*rs)
     &              + c4p*(rp+rn)/(1-c4*rs)*(1-c4*rs))*ssn
     &              - t1p2*ts + t1m2*tn +t1p2*(ts/(1-c1*rs)
     &              + c1*rs*ts/(1-c1*rs)*(1-c1*rs))
     &              - t1m2*(c1p*(rp*tp+rn*tn)/(1-c1p*rs)*(1-c1p*rs)
     &              + tn/(1-c1p*rs))
        bmasp = bmasp - (t1p2*rs-t1m2*rp) + t1p2*rs/(1-c1*rs) 
     &                - t1m2*rp/(1-c1p*rs)
        bmasn = bmasn - (t1p2*rs-t1m2*rn) + t1p2*rs/(1-c1*rs)
     &                - t1m2*rn/(1-c1p*rs)

        wlsp   = wlsp - (b4*d1s+b4p*d1p) + b4*(c4*rs*d1s/((1-c4*rs)
     &                *(1-c4*rs)) + d1s/(1-c4*rs))
     &                + b4p*(c4p*rp*d1s/((1-c4p*rs)*(1-c4p*rs)) 
     &                + d1p/(1-c4p*rs))
        wlsn   = wlsn - (b4*d1s+b4p*d1n) + b4*(c4*rs*d1s/((1-c4*rs)
     &                *(1-c4*rs)) + d1s/(1-c4*rs))
     &                + b4p*(c4p*rn*d1s/((1-c4p*rs)*(1-c4p*rs)) 
     &                + d1n/(1-c4p*rs))
       
      endif
c
ctest      write(*,*) ' energd,potp,potn,bmasp,bmasn,wlsp,wlsn=',
ctest     &             energd,potp,potn,bmasp,bmasn,wlsp,wlsn

      return
      end

c-----initfc------------ part of the hartree-fock package ----------------------
c
      subroutine initfc
c
      include 'params.inc'
c
c     * * * * * * *      initialise force parameters      * * * * * * * *
c
c     initialises auxiliary force parameters from the force parameters
c     given via common /skypar/.
c     the given parameters are the skyrme parameters
c     't0-t8,tyuk', 'x0-x8,xyuk', 'ayuk' and 'power';
c     furthermore the pairing parameters 'pairfp' and 'pairfn',
c     and the coulomb strength 'acoul'
c     all the other parameters in common /skypar/ are derived quantities.
c
      include 'skypar.inc'
      include 'wavfun.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
      data iprib1/0/
c
cold      data x10up/1.0e10/
c
c-------------------------------------------------------------------------------
c
c     determine effective masses according to options 'ifdfms'
c     and 'ifzpe'. in fact, 'amfc' is half the inverse mass in
c     units mev/fm**2.
c
      if(ifdfms.eq.0) then
       amfc(1) = 20.7525
       amfc(2) = 20.7525
      elseif(ifdfms.eq.1) then
       amfc(1) = 20.750
       amfc(2) = 20.722
      elseif(ifdfms.eq.2) then
       amfc(1) = 20.73553
       amfc(2) = 20.73553
cskp       amfc(1) = 20.73
cskp       amfc(2) = 20.73
      elseif(ifdfms.eq.3) then
c       Mail from K.Rutz 15.4.1998 gives NIST-values as most precise:
c         hbar**2/(2m_n) = 20.721 259(6) MeV*fm**2 (NIST)
c         hbar**2/(2m_p) = 20.749 821(6) MeV*fm**2 (NIST)
c
       amfc(1) = 20.74982   ! 20.7505   !20.747
       amfc(2) = 20.72126   ! 20.7219   !20.719
c       amfc(1) = 20.73621
c       amfc(2) = 20.73621
      else
       stop ' IFDFMS must be out of [0..3]'
      endif
      if(ifzpe.ne.1) then
        if(ifdfms.eq.3) then
          shellN  = (1.5*anumb)**0.33333333
          effmfc  = 2.0/(shellN+0.33333333/shellN)      ! N.P. A422 (1984) 157
          amfc(1) = amfc(1)-amfc(1)*effmfc/anumb
          amfc(2) = amfc(2)-amfc(2)*effmfc/anumb
        else
          amfc(1) = amfc(1)-amfc(1)/anumb
          amfc(2) = amfc(2)-amfc(2)/anumb
        endif
      endif
c
c     combinations of skyrme-parameters which occur frequently:
c     the simple assigments are (hopefully) self-explaining.
c     note: in case of yukawa folding, i.e. tyuk.ne.0, the
c           zero-range equivalent strength of the yukawa
c           potential is subtracted from the zer0-range two-
c           body force. thus the t0 still describes the overall
c           strength.
c
cold      aysm3  = one/(fourpi*ayuk*ayuk*ayuk)
c
      t0b    = t0*(0.5+x0)-tyuk*(0.5+xyuk)
      t3b    = t3*(0.5+x3)/6.0
      t0a    = t0*(1.0+0.5*x0)-tyuk*(1.0+0.5*xyuk)
      if(ioptsk.eq.2 .or. ioptsk.eq.3) then
        t3a    = 0.0
        b3dir  = t3*(1.0+0.5*x3)/12.0
c        write(2,'(a,2g12.4)') ' INIT: b0dir,gamdir=',b0dir,gamdir
      else
        t3a    = t3*(1.0+0.5*x3)/6.0
        b3dir  = 0.0
      endif
      if(exmass.ne.zero) then
        coefx  = t3b
        t3b    = 0.0                   ! to disable old isovector part
      endif
      t7a    = t7*(1.0+0.5*x7)/24.0
      t7b    = t7*(0.5+x7)/24.0
cold      t3v    = 0.25*t3
      t3as   = power*0.5*t3a
      t3bs   = power*0.5*t3b
      t8zw   = t8/12.0
      t1p2   = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)
      t1m2   = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)
      if(ioptsk.ne.4) then
        b4     = 0.5*t4
      else
        b4     = -0.5*b4p
        t4     = -b4p
        t1p2   = 0.5*t1m2
      endif
      if(b4p.eq.zero) b4p = b4
      if(ioptsk.eq.4 .and. iprib1.le.10) then
        write(6,'(a,1pg13.5)') ' B1=',t1p2-0.5*t1m2
        iprib1 = 1 + iprib1
      endif
      if(tyuk.eq.zero) then
       t13p2  = 0.125*(1.5*t1+t1*x1*3.0+0.5*t2+t2*x2)
       t13m2  = 0.125*(3.0*t1+t1*x1*1.5-t2-t2*x2*0.5)
      else
       t13p2  = 0.0
       t13m2  = 0.0
      endif
      t1m20  = (t1-t2)/8.0
      t1p2x  = (x1*t1+x2*t2)/8.0
      t5p6   = (t5+t5*x5*0.5+t6+t6*x6*0.5)/24.0
      t5m6   = (0.5*t5+t5*x5-0.5*t6-t6*x6)/24.0
      t53p6  = (1.5*t5+t5*x5*3.0+0.5*t6+t6*x6)/96.0
      t53m6  = (3.0*t5+t5*x5*1.5-t6-t6*x6*0.5)/96.0
      t5m60  = (t5-t6)/48.0
      t5p6x  = (x5*t5+x6*t6)/48.0
      t5m60h = 0.5*t5m60
      t5p6xh = 0.5*t5p6x
      t53m6d = 2.0*t53m6
      if(ayuk.ne.zero) then
       ayms2  = 1.0/(ayuk*ayuk)
      else
       ayms2  = zero
      endif
      tya    = tyuk*(1.0+0.5*xyuk)
      tyb    = tyuk*(0.5+xyuk)
c
c
c     switch force-option: normal parametrisation --> isky=8
c                          extend. parametrisation --> isky=10
c                          folded yukawa param. --> isky=9
c
      if(abs(t5)+abs(t6)+abs(t7)+abs(t8).gt.zero) then
       if(tyuk.eq.zero) then
        isky   = 10
       else
        write(06,'(a,a)') ' incompatible parameters: ',
     *            ' tyuk and t5-t8 not simultanously treated.'
        stop
       endif
      elseif(tyuk.ne.zero) then
       isky   = 9
      else
       isky   = 8
      endif
      if(isky.gt.8 .and. ipair.eq.3)
     & stop ' extended Skyrme not compatible with Lipkin-Nogami'
c
c     option for the force translated to logical variables
c     for faster switching inside the do-loop.
c      'texcha' = exchange contribution to coulomg force
c      'tyukaw' = finite range two body force with yukawa folding
c      'texten' = extended parametrisation ( using t5-t8)
c
      texcha= (ifex.ne.0)
      tyukaw= (isky.eq.9)
      texten= (isky.gt.9)
      tssmas= (t9.ne.zero)
      tpair = (iparac.ge.7)
      texmas= (exmass.ne.zero)
      tdirpa= (ioptsk.eq.2) .or. (ioptsk.eq.3)
      trmfki= ioptsk.eq.4
c
      return
      end
c-----exfunc------------ part of the hartree-fock package --------------------
c
      double precision function exfunc(rhoact)
c
c     the explicit exchange functional. Input is
c      'rhoact'  via list
c      'exmass'  via  common /skypar/
c      'coefx'   via  common /skypar/
c
      include 'params.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------
c
c
c     reduced Fermi-momentum 'eta' employs 29.6088132=3*pi**2
c
      eta    = (29.6088132*rhoact)**0.333333333d0*2.0/exmass
      eta2   = eta*eta
      if(ioptsk.ne.3) then
        exfunc = coefx/24.0 * (
     &           1.5d0*eta2*eta2-eta2
     &           -4.0*eta2*eta*atan(eta)
     &           +(3.0d0*eta2+1.0d0)*log(1.0d0+eta2)      )
      else
        exfunc = coefx*eta**6/(one+0.11111111111*eta2)
      endif
c
      return
      end
c-----expot ------------ part of the hartree-fock package --------------------
c
      double precision function expot(rhoact)
c
c     the explicit mean field potental from the exchange functional. Input is
c      'rhoact'  via list
c      'exmass'  via  common /skypar/
c      'coefx'   via  common /skypar/
c
      include 'params.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------
c
c
c     reduced Fermi-momentum 'eta' employs 29.6088132=3*pi**2
c
      eta    = (29.6088132*rhoact)**0.333333333d0*2.0/exmass
      eta2   = eta*eta
      if(ioptsk.ne.3) then
        expot  = coefx*eta/(3.0d0*rhoact) * (
     &           0.5d0*eta2*eta
     &           -eta2*atan(eta)
     &           +0.5d0*eta*log(1.0d0+eta2)      )     * 0.5
      else
        expot  = coefx*eta/(3.0d0*rhoact)
     &          *eta**5*(6.0+0.444444444444*eta2)
     &          /(one+0.11111111111*eta2)**2
      endif
c
      return
      end
c-----exenre------------ part of the hartree-fock package --------------------
c
      double precision function exenre(rhoact)
c
c     the rearrangement energy from the exchange functional. Input is
c      'rhoact'  via list
c      'exmass'  via  common /skypar/
c      'coefx'   via  common /skypar/
c
      include 'params.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------
c
c
c     reduced Fermi-momentum 'eta' employs 29.6088132=3*pi**2
c
      eta    = (29.6088132*rhoact)**0.333333333d0*2.0/exmass
      eta2   = eta*eta
      atane  = eta2*atan(eta)
      xloge  = log(1.0d0+eta2)
      if(ioptsk.ne.3) then
        expot  = coefx*eta/(3.0d0*rhoact) * (
     &           0.5d0*eta2*eta-atane+0.5d0*eta*xloge      )    *   0.5
        exfunc = coefx/24.0 * (
     &           1.5d0*eta2*eta2-eta2-4.0*eta*atane
     &           +(3.0d0*eta2+1.0d0)*xloge      )
        exenre = exfunc - 0.5d0*rhoact*expot
      else
        exenre = coefx/27.0 * eta**8 /(one+0.11111111111*eta2)**2
      endif
c
      return
      end
c-----dirfun------------ part of the hartree-fock package --------------------
c
      double precision function dirfun(rhoact)
c
c     the direct-part functional in Pade approximation. Input is
c      'rhoact'  via list
c      'b3dir'  via  common /skypar/
c      'power'   via  common /skypar/
c
      include 'params.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------
c
c
cold      dirfun = 0.5*b0dir*rhoact*rhoact/(one+gamdir*rhoact)
      dirfun = b3dir*rhoact**3/(one+power*rhoact)
c
      return
      end
c-----dirpot------------ part of the hartree-fock package --------------------
c
      double precision function dirpot(rhoact)
c
c     the m.f. potential for direct-part functional in Pade approximation. 
c     Input is
c      'rhoact'  via list
c      'b3dir'  via  common /skypar/
c      'power'   via  common /skypar/
c
      include 'params.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------
c
c
cold      dirpot = 0.5*b0dir*rhoact*(2.0+gamdir*rhoact)
cold     &         /(one+gamdir*rhoact)**2
      dirpot = b3dir*rhoact*rhoact*(3.0+2.0*power*rhoact)
     &         /(one+power*rhoact)**2
c      write(2,'(a,2g12.4)') ' DIRPOT: b0dir,gamdir=',b0dir,gamdir
c
      return
      end
c-----direnr------------ part of the hartree-fock package --------------------
c
      double precision function direnr(rhoact)
c
c     the rearr. energy for direct-part functional in Pade approximation. 
c     Input is
c      'rhoact'  via list
c      'b3dir'  via  common /skypar/
c      'power'   via  common /skypar/
c
      include 'params.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------
c
c
      direnr = -0.5*b3dir*rhoact**3/(one+power*rhoact)**2
c
      return
      end












