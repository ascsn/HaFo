c-----skfunc------------ part of the hartree-fock package ---------------------
c
      subroutine skfunc(rp,rn,rc,d1p,d1n,d2p,d2n,tp,tn,tlp,tln,
     &                  sp,sn,ssp,ssn,
     &                  frp,frn,potp,potn,bmasp,bmasn,bmaslp,bmasln,
     &                  wlsp,wlsn,energd)
c
      include 'params.inc'
c
      parameter (diffrho=5.0D-5)    ! for finite derivative
      parameter (diffi=0.5D0/diffrho)    ! for finite derivative
c
c     the basic Skyrme energy density and derived quantities.
c
c     Input per list:
c      rp,rn        = densities
c      rc           = charge density
c      d1p,d1n      = derivative of densities
c      d2p,d2n      = Laplacian of densities
c      tp,tn        = kinetic energy densities
c      tlp,tln      = L**2 densities
c      sp,sn        = spin orbit densities
c      ssp,ssn      = derivative of spin-orbit densities
c      frp,frn      = folded densities (for Yukawa term)
c                     
c     Forces are input per common /skypar/ und /tswitc/.
c     Ouput per list:
c      energd        = energy density
c      potp,potn     = local potentials
c      bmasp,bmasn   = inverse mass
c      bmaslp,bmasln = inverse L**2 mass
c      wlsp,wlsn     = spin-orbit potential
c
c     !!! Note: the option 'texten' needs yet careful checking
c
      include 'skypar.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'constn.inc'
c
      parameter (onedtwo=1.0D0/2.0D0)
      parameter (onedsix=1.0D0/6.0D0)
      parameter (hc=197.3285879)
c
c------------------------------------------------------------------------------
c
c     switch to QMC model
c
      IF(ioptsk==12) THEN
        CALL SetParameters
        potp = hc*UPotProton(rp,rn,tp,tn,d2p,d2n,ssp,ssn,d1p,d1n,
     &                       0D0,0D0)
        potn = hc*UPotNeutron(rp,rn,tp,tn,d2p,d2n,ssp,ssn,d1p,d1n,
     &                        0D0,0D0)
        wlsp = hc*WPotProton(rp,rn,tp,tn,d2p,d2n,ssp,ssn,d1p,d1n,
     &                       0D0,0D0)
        wlsn = hc*WPotNeutron(rp,rn,tp,tn,d2p,d2n,ssp,ssn,d1p,d1n,
     &                        0D0,0D0)
        bmasp = hc*OneOver2MeffProton(rp,rn,tp,tn,d2p,d2n,ssp,ssn,
     &                                0D0,0D0,0D0,0D0)
        bmasn = hc*OneOver2MeffNeutron(rp,rn,tp,tn,d2p,d2n,ssp,ssn,
     &                                 0D0,0D0,0D0,0D0)
        energd = hc*EnergyDensity(rp,rn,tp,tn,d2p,d2n,ssp,ssn,d1p,d1n,
     &                            0D0,0D0)
!        WRITE(6,'(10(1pg13.5)/7x,8(1pg13.5))')
!     &   'SKFUNC:',rp,rn,tp,tn,d2p,d2n,ssp,ssn,d1p,d1n,
!     &    potp,potn,wlsp,wlsn,bmasp,bmasn,energd
!        STOP
        IF(texcha) THEN
          excfac = -1.418033*acoul
          add    = excfac*rc**0.33333333333333D0
          potp   = add + potp
          energd = 0.75D0*add*rc + energd
!          WRITE(2,'(a,1pg13.5)') 'excfac=',excfac
        END IF
        RETURN    
      END IF
c
c     deactivate l-mass
c
      bmaslp=0D0
      bmasln=0D0
c
c     total densities and other auxiliary quantities
c
      rs  = rp+rn
      ts  = tp+tn
      ss  = sp+sn
      sss = ssp+ssn
      d2s = d2p+d2n
      d1s = d1p+d1n
c      if(talpinv) then
c        rspow  = (rs+betinv*(rs*rs))/(1D0+alpinv*rs)
c      else
        rspow  = rs**power
c      endif
      if(rhoswm.gt.zero) then
        tanfac = tanh(rs/rhoswm)
      else
        tanfac = 1.0
      endif
      tyab   = tya-tyb
      t8zwh  = half*t8zw
      riftm1 = riftm
      if(.not.(ioptsk.eq.10 .or. ioptsk.eq.11)) then
        riftm  = one*iftm         !one*iftm
        riftm1 = riftm
      endif
      if(ioptsk.eq.7) then
        riftm = 0D0
        riftm1 = 0D0
      endif
      IF(TADSIC) THEN
        excfac = -1.418033*acoul*(1D0-(1D0/nprot)**0.33333333333D0)
      ELSE
        excfac = -1.418033*acoul
      END IF
      t1m20r = t1m20*riftm1
      t1p2xr = t1p2x*riftm
!      WRITE(2,*) ' SKFUNC:riftm,riftm1=',riftm,riftm1
c
c     basic (density dependent) part of Skyrme functional 
c
      if(ioptsk.ne.7) then
        if(talpinv) then
          alpeff = alpinv+betinv
          rsinv = 1D0/(1D0+alpeff*rs)
          rsfac = (rs+betinv*(rs*rs))*rsinv*alpeff
          rsfac2 = (1D0+betinv*(rs+rs)+(alpeff*betinv)*rs*rs)*
     &             (rsinv*rsinv)*alpeff
          potp =  t0a*rs-t0b*rp  
     &           +rsfac*(t3a*rs-t3b*rp)
     &           +half*rsfac2*(t3a*(rs*rs)-t3b*(rp*rp+rn*rn))
          potn =  t0a*rs-t0b*rn  
     &           +rsfac*(t3a*rs-t3b*rn)
     &           +half*rsfac2*(t3a*(rs*rs)-t3b*(rp*rp+rn*rn))
          energd =  +half*(t0a*(rs*rs)-t0b*(rp*rp+rn*rn))
     &              +half*rsfac*(t3a*(rs*rs)-t3b*(rp*rp+rn*rn))
        elseif(tpowerp) then
          rtilde = rn-rp
          rspowp = rs**powerp
          potp = ( t0a*rs-t0b*rp  
     &            + b3eff*rspow*rs
     &            - b3peff*rspowp*(rtilde*rtilde)/max(rs,xsmall) 
     &            + b3peff2*rspowp*rtilde                      )
          potn = ( t0a*rs-t0b*rn  
     &            + b3eff*rspow*rs
     &            - b3peff*rspowp*(rtilde*rtilde)/max(rs,xsmall) 
     &            - b3peff2*rspowp*rtilde                      )
          energd =  +half*(t0a*(rs*rs)-t0b*(rp*rp+rn*rn))
     &              +b3scal*rspow*(rs*rs)
     &              -b3vect*rspowp*(rtilde*rtilde)
        else
          potp = ( t0a*rs-t0b*rp  
     &           +rspow*(t3a*rs-t3b*rp +t3as*rs
     &                   -t3bs*(rp*rp+rn*rn)/max(rs,xsmall) ) )
          potn = ( t0a*rs-t0b*rn  
     &           +rspow*(t3a*rs-t3b*rn +t3as*rs
     &                   -t3bs*(rp*rp+rn*rn)/max(rs,xsmall) ) )
          energd =  +half*(t0a*(rs*rs)-t0b*(rp*rp+rn*rn))
     &              +half*rspow*(t3a*(rs*rs)-t3b*(rp*rp+rn*rn))
        endif
        if(tb3til) then
          potp   = (b3tilde*(2.0D0+power)/3.0D0)*rp**(1.0D0+power)
     &            + potp
          potn   = (b3tilde*(2.0D0+power)/3.0D0)*rn**(1.0D0+power)
     &            + potn
          energd = b3tilde/3.0D0*
     &             (rp**(2.0D0+power) + rn**(2.0D0+power))
     &            + energd
        endif
        if(tb12) then
          rspowm1 = rspow/rs
          potp   = -b12*rspow*d2s
     &             -b12*(power*0.5D0)*rspowm1*(d1s*d1s)
     &             -b12p*rspow*d2p
     &             -b12p*power*rspowm1*d1s*d1p
     &             +b12p*(power*0.5D0)*rspowm1*(d1p*d1p+d1n*d1n)
     &            + potp
          potn   = -b12*rspow*d2s
     &             -power*0.5D0*b12*rspowm1*(d1s*d1s)
     &             -b12p*rspow*d2n
     &             -b12p*power*rspowm1*d1s*d1n
     &             +b12p*(power*0.5D0)*rspowm1*(d1p*d1p+d1n*d1n)
     &            + potn
          energd = 0.5D0*b12*rspow*(d1s*d1s) 
     &            +0.5D0*b12p*rspow*(d1p*d1p+d1n*d1n)
     &            + energd
        endif
        if(tdens2) then
          rspow7 = rs**pow7
          t7b  = onedsix*b7p
          t7a  = onedsix*b7+onedtwo*t7b
          b7a  = t7a*(pow7+2.0D0)*onedtwo
          b7pa = t7b
          b7pb = t7b*pow7*onedtwo
          potp =  rspow7*(b7a*rs - b7pa*rp
     &              -b7pb*(rp*rp+rn*rn)/max(rs,xsmall) ) 
     &          + potp
          potn =  rspow7*(b7a*rs - b7pa*rn
     &              -b7pb*(rp*rp+rn*rn)/max(rs,xsmall) ) 
     &          + potn
          energd = rspow7*onedtwo*(t7a*(rs*rs)-t7b*(rp*rp+rn*rn))
     &           + energd
        endif
        if(talpin2) then
          alpeff = alpin2+betin2
          rsinv = 1D0/(1D0+alpeff*rs)
          rsfac = (rs+betin2*(rs*rs))*rsinv*alpeff
          rsfac2 = (1D0+betin2*(rs+rs)+(alpeff*betin2)*rs*rs)*
     &             (rsinv*rsinv)*alpeff
          potp = potp
     &           +rsfac*(b7*rs-b7p*rp)
     &           +half*rsfac2*(b7*(rs*rs)-b7p*(rp*rp+rn*rn))
          potn = potn 
     &           +rsfac*(b7*rs-b7p*rn)
     &           +half*rsfac2*(b7*(rs*rs)-b7p*(rp*rp+rn*rn))
          energd = energd
     &             +half*rsfac*(b7*(rs*rs)-b7p*(rp*rp+rn*rn))
        endif
        if(tdenslow) then
          rspow7 = rs**pow7
          rspow3 = rs**0.3333333333333333D0
          rsinv = 1D0/(1D0+betin2*rspow3)
          rsinv2 = rsinv*rsinv
          rpn2 = (rp*rp+rn*rn)/max(rs,xsmall)
          betin3 = betin2/3D0
!          ignore factor 1/2
          potp = potp
     &          +b7*rs*rspow7*rsinv*((2+pow7)-betin3*rspow3*rsinv)
     &          -b7p*rspow7*rsinv*(2D0*rp+pow7*rpn2
     &                             -betin3*rspow3*rsinv*rpn2)
          potn = potn
     &          +b7*rs*rspow7*rsinv*((2+pow7)-betin3*rspow3*rsinv)
     &          -b7p*rspow7*rsinv*(2D0*rn+pow7*rpn2
     &                             -betin3*rspow3*rsinv*rpn2)
          energd = energd
     &            +(b7*rs*rs-b7p*(rp*rp+rn*rn))*rspow7*rsinv
        endif
      else
        energd = bhf_poten(rp,rn)
        potp   = (bhf_poten(rp+diffrho,rn)-bhf_poten(rp-diffrho,rn))
     &          *diffi
        potn   = (bhf_poten(rp,rn+diffrho)-bhf_poten(rp,rn-diffrho))
     &          *diffi
      endif
      bmasp = zero
      bmasn = zero
      wlsp  = zero
      wlsn  = zero
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
      if(ioptsk.ne.7) then
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
       energd = half*(rs*tya*frs-rp*tyb*frp-rn*tyb*frn) + energd
      endif
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
!      wlsp   = riftm*(t1m20*sp-t1p2x*ss) + wlsp
!      wlsn   = riftm*(t1m20*sn-t1p2x*ss) + wlsn
      energd = (t1p2*rs*ts-t1m2*(rp*tp+rn*tn))*tanfac 
     &        + half*(t1m20r*(sp*sp+sn*sn)-t1p2xr*ss*ss)
     &        + half*(-t13m2*d2s*rs + t13p2*(d2p*rp+d2n*rn))
     &        + energd
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
      if(tssmas) then    
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
      wlsp   = (b4*d1s+b4p*d1p+(t1m20r*sp-t1p2xr*ss)) + wlsp
      wlsn   = (b4*d1s+b4p*d1n+(t1m20r*sn-t1p2xr*ss)) + wlsn
      energd = -b4*rs*sss - b4p*(rp*ssp+rn*ssn) + energd
c
c     extended spin-orbit terms
c
      if(texten) then
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
c     note:  t1p2 = b1 , t1m2 = b1p
c

ctest      write(2,*) 
ctest     &  ' TKINDIR: tkindir,c1,c1p,c4,c4p=',tkindir,c1,c1p,c4,c4p
      if(tkindir) then
ctest        write(2,'(4(1pg16.8))')
ctest     &       t1p2*c1,t1m2*c1p,b4*c4,b4p*c4p
ctest        stop
        potp = potp + ( -2*c4*rs*sss 
     &                  -c4p*(rp*ssp+rn*ssn+rs*ssp) )
     &              + (2*c1*rs*ts-c1p*(rp*tp+rn*tn+rs*tp))
        potn = potn + ( -2*c4*rs*sss 
     &                  -c4p*(rp*ssp+rn*ssn+rs*ssn) )
     &              + (2*c1*rs*ts-c1p*(rp*tp+rn*tn+rs*tn))
       bmasp = bmasp + c1*rs*rs - c1p*rs*rp
       bmasn = bmasn + c1*rs*rs - c1p*rs*rn
       wlsp   = wlsp + (c4*(rs+rs)*d1s+c4p*(rs*d1p+d1s*rp))
       wlsn   = wlsn + (c4*(rs+rs)*d1s+c4p*(rs*d1n+d1s*rn))
       energd =  c1*rs*rs*ts - c1p*rs*(rp*tp+rn*tn)
     &         - c4*rs*rs*sss - c4p*rs*(rp*ssp+rn*ssn)
     &         + energd
cbaumann: the following lines contain the old version which
cbaumann: complies with the diploma thesis of I. Baumann
cbaumann        potp = potp + ( -2*b4*c4*rs*sss 
cbaumann     &                  -b4p*c4p*(rp*ssp+rn*ssn+rs*ssp) )
cbaumann     &              + (2*t1p2*c1*rs*ts-t1m2*c1p*(rp*tp+rn*tn+rs*tp))
cbaumann        potn = potn + ( -2*b4*c4*rs*sss 
cbaumann     &                  -b4p*c4p*(rp*ssp+rn*ssn+rs*ssn) )
cbaumann     &              + (2*t1p2*c1*rs*ts-t1m2*c1p*(rp*tp+rn*tn+rs*tn))
cbaumann       bmasp = bmasp + t1p2*c1*rs*rs - t1m2*c1p*rs*rp
cbaumann       bmasn = bmasn + t1p2*c1*rs*rs - t1m2*c1p*rs*rn
cbaumann       wlsp   = wlsp + (b4*c4*(rs+rs)*d1s+b4p*c4p*(rs*d1p+d1s*rp))
cbaumann       wlsn   = wlsn + (b4*c4*(rs+rs)*d1s+b4p*c4p*(rs*d1n+d1s*rn))
cbaumann       energd =  t1p2*c1*rs*rs*ts - t1m2*c1p*rs*(rp*tp+rn*tn)
cbaumann     &         - b4*c4*rs*rs*sss - b4p*c4p*rs*(rp*ssp+rn*ssn)
cbaumann     &         + energd
      endif
c
c     ioptsk.eq.6 --> rho-dependent kinetic and l*s, inverse version
c


c     b1-, b1p-, b4- and b4p-terms that appeared before must be
c     substracted to not add them twice; same procedure with wlsp, wlsn, bmasp and bmasn

      if(tkininv) then
        stop ' IOPTSK=6 presently disabled '
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

      return
      end

c-----bhf_poten-----------------------------------------------------
c
      double precision function bhf_poten(rhop,rhon)
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   !
c
c     BHF energy functional from Baldo et al, NPA 736 (2004) 241
c
c
c     coefficients symmetric matter low rho
c
      parameter (b1=-116.734211421291D0)
      parameter (b2= 310.796396961038D0)
      parameter (b3=-584.213036193276D0)
      parameter (b4= 587.659458730963D0)
      parameter (b5=-290.107854880778D0)
      parameter (b6=  55.540461239575D0)
c
c     coefficients symmetric matter large rho
c
      parameter (esym0=1.0D0*(b1+1.0D0*(b2+1.0D0*(b3+1.0D0*(b4
     &                +1.0D0*(b5+1.0D0*b6)))))   )
      parameter (a1=-14.4391981216314D0)
      parameter (a2= 16.1424105528328D0)
c
c     critical density netween 'low' and 'large'
c
      parameter (rho0=0.155)
c
c     coefficients asymmetric matter 
c
      parameter (c1=-40.1046417029240D0)
      parameter (c2= 41.3865460399823D0)
      parameter (c3=-32.1381502519477D0)
      parameter (c4= 15.1308678586138D0)
      parameter (c5= -2.71394288779458D0)
c
c-------------------------------------------------------------------------
c
c     the symmetric EoS
c
      rhotot = rhop+rhon
      x      = rhotot/rho0
      if(x.lt.1.0D0) then
        esym  = x*(b1+x*(b2+x*(b3+x*(b4+x*(b5+x*b6)))))
      else
        y     = x-1.0D0
        esym  = y*(a1+a2*y)+esym0
      endif
c
c     the asymmetric EoS
c
      easym = x*(c1+x*(c2+x*(c3+x*(c4+x*c5))))
c
c     compose finally
c
      if(rhop.eq.rhon) then
        zeta2 = 0.0D0
      else
        zeta2  = (rhon-rhop)/rhotot
        zeta2  = zeta2*zeta2
      endif
      bhf_poten = rhotot*(easym*zeta2+(1.0D0-zeta2)*esym)
c      write(6,'(1x,4(1pg13.5))') rhotot,zeta2,esym,easym
c
      return
      end

!      include 'Special_QMC_HFEquations_sub.f'
      include 'Special_QMC_With_Pion_HFEquations_sub.f'


