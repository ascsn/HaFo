... in nucmat ...
      if(ioptsk.eq.4) then
          bmsp = 0.0
c          call rmfkin(rhonm0*0.5,rhonm0*0.5,zero,zero,zero,zero,
c     &                zero,zero,1,energk,enreak,
c     &                potp,potn,bmsp,bmsn,plsp,pldn)        
          effmiv = bmsp/hqqd2m + effmiv
c          write(6,'(a,4g12.4)')
c     &     ' rhonm0,bmsp,effmiv,hqqd2m=',rhonm0,bmsp,effmiv,hqqd2m
      endif

...  in lanpar ...
      if(ioptsk.eq.4) then
          bmsp = 0.0
          call rmfkin(rhonm0*0.5,rhonm0*0.5,zero,zero,zero,zero,
     &                zero,zero,1,energk,enreak,
     &                potp,potn,bmsp,bmsn,plsp,pldn)        
          effmiv = bmsp/hqqd2m + effmiv
      endif

... in xkappa ...
      if(ioptsk.eq.4) then
        bmsp = 0.0
        bmsn = 0.0
        call rmfkin(rhop,rhon,zero,zero,zero,zero,zero,zero,1,
     &              energk,enreak,
     &              potp,potn,bmsp,bmsn,plsp,pldn)        
        amfcpp = bmsp + amfcpp
        amfcpn = bmsn + amfcpn
      endif

... in enernm ...
        if(ioptsk.eq.4) then
          call rmfkin(rhop,rhop,taup,taup,zero,zero,zero,zero,0,
     &                energk,enreak,
     &                potp,potn,bmsp,bmsn,plsp,pldn)        
          enrg   = energk + enrg
        endif

        if(ioptsk.eq.4) then
          call rmfkin(rhop,rhon,taup,taun,zero,zero,zero,zero,0,
     &                energk,enreak,
     &                potp,potn,bmsp,bmsn,plsp,pldn)        
          enrg   = energk + enrg
        endif


... in skypot ...
      else
c                                     finally add RMF-kinetic terms, all here
        if(trmfki) then
          if(i.eq.ia) then
            d1pa = d2p/3.0
            d1na = d2n/3.0
          else
            d1pa = d1p/x
            d1na = d1n/x
          endif
          call rmfkin(rp,rm,tp,tm,ssp,ssm,d1pa,d1na,1,
     &                energk,enreak,
     &                pot(i),pot(j),bmass(i),bmass(j),potls(i),potls(j))      
          enrear = -2.0*(x*x)*enreak + enrear
        endif



c-----rmfkin------------ part of the hartree-fock package --------------------
c
      subroutine rmfkin(rp,rn,tp,tn,sjsp,sjsn,drp,drn,ifpot,
     &                  energk,enreak,
     &                  potp,potn,bmsp,bmsn,plsp,plsn)
c
c     the kinetic+l*s functional as motivated by the RMF.
c     input:
c       rp,rn      = proton- and neutron-density
c       tp,tn      = kinetic densities
c       sjsp,sjsn  = nabla of spin-orbit densities
c       drp,drn    = gradient of densities
c       ifpot      = switch to accumulate potential
c                    (0 --> only energy is returned)
c     input via common /forces/
c       b1kin      = overall force parameter
c       b5kin      = force parameter in denominator (set b1kin*m/hbar**2?)
c     in- and output:
c       potp,potn  = contribution to the local potential (added up)
c       bmsp,bmsn  = contribution to the inv. effective mass (added up)
c       plsp,plpn  = contribution to the spin-orbit potential (added up)
c     output:
c       energk     = contribution to the integrand of total energy
c       enreak     = contribution to the integrand of rearrangement energy
c
      include 'params.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
      data ipri/0/
c
c------------------------------------------------------------------------
c
c     preparatory combinations
c
      rs     = rp+rn
      tsjs   = tp+tn - (sjsp+sjsn)
      drs    = drp+drn
      denom  = one/(one-b5kin*rs)
      denom2 = 0.5*b1kin*(denom*denom)
      denom3 = denom2*denom*(one+b5kin*rs)
c
c     the terms
c
      energk = denom2*rs*tsjs
      if(ifpot.ne.1) return
c
      uloc   = denom3*tsjs
      enreak = 0.5*(energk-rs*uloc)      
      potp = uloc + potp
      potn = uloc + potn
      bmsa = denom2*rs
      bmsp = bmsa + bmsp
      bmsn = bmsa + bmsn
      plsa = denom3*drs
      plsp = plsa + plsp
      plsn = plsa + plsn
c
      if(ipri.lt.0) then
        ipri = 1 + ipri
        write(6,'(a,9g12.4)')
     &   ' b1kin,b5kin,rs,ts,energk,enreak,uloc,bmsa,plsa=',
     &     b1kin,b5kin,rs,tp+tn,energk,enreak,uloc,bmsa,plsa
      endif
c
      return
      end


