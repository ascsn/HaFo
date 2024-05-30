c-----sepenr------------------------------------------------------------
c
      subroutine sepenr()
c
      include 'params.inc'
c
c     computes and prints the separate energy contributions
c     for a finite nucleus:
c      energm = output: E/A in MeV
c      enr0sc = b_0 part
c      enr0vc = b'_0 part
c      enr1sc = b_1 part
c      enr1vc = b'_1 part
c      enr3sc = b_3 part
c      enr3vc = b'_3 part
c      ekinA  = kinetic energy
c     The force parameters are entered via common /skypar/
c
c     This routine is to be called in 'subroutine hafosu'
c     after the end of the iteration (line '1009 continue')
c     and before 'call radius'.
c
      include 'skypar.inc'
      include 'grid.inc'
      include 'wavfun.inc'
      include 'result.inc'
      include 'numpar.inc'
      include 'potinv.inc'
      include 'constn.inc'
c
c-----------------------------------------------------------------------
c
c     prepare Coulomb potential
c
      if(acoul.ne.zero) then
       call coulmb(rhochr(1),binv,acoul)
      else
       do ir=1,ngrid
         binv(ir) = 0.0
       enddo
      endif
c
      enr0sc = 0.0
      enr0vc = 0.0
      enr1sc = 0.0
      enr1vc = 0.0
      enr2sc = 0.0
      enr2vc = 0.0
      enr3sc = 0.0
      enr3vc = 0.0
      enr4   = 0.0
      enrlst = 0.0
      ecoul  = 0.0
      ecolex = 0.0
      if(ifex.eq.1) then
        excfac = -1.418033*acoul
      else
        excfac = 0.0
      endif
      do ir=2,ngrid
        ract   = (ir-1)*rstep
        r2     = ract*ract
        rhop   = rho(ir)
        rhon   = rho(ir+ngrid)
        rhonm  = rhop+rhon
        rhopow = max(rhonm,xsmall)**power
        taup   = tau(ir)
        taun   = tau(ir+ngrid)
        taunm  = taup+taun
        d2p    = d2r(ir)
        d2n    = d2r(ir+ngrid)
        d2nm   = d2p+d2n
        sp     = sj(ir)
        sn     = sj(ir+ngrid)
        snm    = sp+sm
        ssp    = sjs(ir)
        ssn    = sjs(ir+ngrid)
        ssnm   = ssp+ssn
        enr0sc = r2  * t0a*rhonm*rhonm                  + enr0sc
        enr0vc = -r2 * t0b*(rhop*rhop+rhon*rhon)        + enr0vc
        enr1sc = r2  * t1p2*rhonm*taunm                 + enr1sc
        enr1vc = -r2 * t1m2*(rhop*taup+rhon*taun)       + enr1vc
        enr2sc = -r2 * t13m2*rhonm*d2nm                 + enr2sc
        enr2vc = r2  * t13p2*(rhop*d2p+rhon*d2n)        + enr2vc
        enr3sc = r2  * rhopow*t3a*rhonm*rhonm           + enr3sc
        enr3vc = -r2 * rhopow*t3b*(rhop*rhop+rhon*rhon) + enr3vc
        enr4   = -r2 * (b4*rhonm*ssnm
     &                  +b4p*(rhop*ssp+rhon*ssn)  )     + enr4
        enrlst = r2  * (t1m20*(sp*sp+sn*sn)-t1p2x*snm*snm) + enrlst
        ecoul  = r2  * rhop*binv(ir)                    + ecoul
        ecolex = r2  * excfac*rhochr(ir)**1.3333333     + ecolex
      enddo
c
      enr0sc = enr0sc*frpirs*half
      enr0vc = enr0vc*frpirs*half    !  0.25
      enr1sc = enr1sc*frpirs               
      enr1vc = enr1vc*frpirs
      enr2sc = enr2sc*frpirs*half               
      enr2vc = enr2vc*frpirs*half
      enr3sc = enr3sc*frpirs*half
      enr3vc = enr3vc*frpirs*half    !  0.25
      enr4   = enr4  *frpirs
      enrlst = enrlst*frpirs*iftm*half
      ecoul  = ecoul *frpirs*half
      ecolex = ecolex*frpirs*0.75
c
      write(6,'(/a/a,5(/a,3(1pg15.7)),2(/a,26x,1pg15.7)/a,
     &         2(/a,1pg15.7)/a,3(1pg15.7)/a,26x,1pg15.7/
     &           a,3(1pg15.7)/)')
     &  ' Print separate contributions to SHF energy:',
     &  '           isoscalar    isovector      total   ',
     &  ' t0-terms:',enr0sc,enr0vc,enr0sc+enr0vc,
     &  ' t1-terms:',enr1sc,enr1vc,enr1sc+enr1vc,
     &  ' t2-terms:',enr2sc,enr2vc,enr2sc+enr2vc,
     &  ' t3-terms:',enr3sc,enr3vc,enr3sc+enr3vc,
     &  ' t0+t3   :',enr3sc+enr0sc,enr3vc+enr0vc,
     &               enr3sc+enr3vc+enr0sc+enr0vc,
     &  ' ls-term :',enr4,
     &  ' ls-t1_t2:',enrlst,
     &  '             proton      neutron       total   ', 
     &  ' Coulomb :',ecoul*acoul,
     &  ' Coul.exc:',ecolex,
     &  ' kinetic :',ekinpr,ekinne,(ekinpr+ekinne),
     &  ' TOTAL   :',(ekinpr+ekinne)+enr0sc+enr0vc+enr1sc+enr1vc
     &               +enr2sc+enr2vc+enr3sc+enr3vc+enr4+enrlst
     &               +acoul*(ecoul+ecolex)-zpe,
     &  ' t0&t3   :',enr0sc+enr3sc,enr0vc+enr3vc,
     &               enr0sc+enr0vc+enr3sc+enr3vc
c
      return
      end







