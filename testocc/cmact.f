c-----cmact ------------ part of the hartree-fock package ----------------------
c
      subroutine cmact(wfin,cwfout,d1wfin,d2wfin,jsa,iso,weighta,
     &                 addspe,adddel)
c
      include 'params.inc'
c
c     action of the center-of-mass correction
c
c                     2
c         |hwfout> = P  |wfin> - sum  P|wf > (fac) <wf |P|wfin>
c                     cm            b     b           b
c
c     where '(fac)' is a weight composed of pairing factors and
c     angular momentum couplings.
c
c     List input parameters are:
c       wfin    = input wave function field ,
c       cwfout  = output wave function field ,
c       d2wfin  = second derivative if 'wfin'.
c       jsa     = combined angular momentum quantum number
c                 of wfin and hwfout,
c       iso     = isospin of the wave functions,
c       weighta = occupation weight of input state
c       addspe  = c.m. correction to s.p.e
c       adddel  = c.m. correction to gap
c
c     the fields for the mean-field hamiltonian are transferred
c     via common /wavfun/:
c       pot     = potential
c       potls   = spin-orbit potential
c       bmass   = inverse effective mass
c       d1b     = derivative of inverse effective mass
c       wf      = all wavefunctions
c       wfs     = first derivatives of all wavefunctions
c
c     grid properties are transferred via common /grid/
c
      dimension wfin(kg1),cwfout(kg1),d1wfin(kg1),d2wfin(kg1)
c
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'grid.inc'
      include 'constn.inc'
c
      data zpetot,zpedir,zpeex/3*0.0D0/
      logical ttest
      data ttest/.false./
      data v2lim/1.0D-3/                ! limiting v2a in uva*uvb term
c
c------------------------------------------------------------------------------
c
c     the quantum numbers of the wave function
c      la      = orbital angular momentum of state  wfin>
c      parity  = parity for the radial wave function (note that it
c                contains an extra factor 'r' such that the 1s1/2
c                state has parity=-1 etc)
c      isoadd  = stack-pointer for potential fields pot, potls
c                and isoadd
c      xl0     = spin-orbit splitting for state  wfin>
c
      isnpad = (iso-1)*kjmx
      xja    = float((jsa-1)/2)+0.5
      la     = jsa/2
      xla    = one*la
      xla2   = xla*xla+xla       
      xlp    = (la*la+la)
      parity = (mod(la,2)*2-1)
      isoadd = (iso-1)*ngrid
      xl0    = wls(jsa)
      v2a    = weighta
      uva    = sqrt(max(v2a-v2a*v2a,xsmall))
      h2mfac = half*(amfc(1)+amfc(2))/anumb
      if(temprt.gt.zero) 
     &  stop ' CMACT not yet adapted for nonzero temperature'
c      call derv25(wfin(1),d2wfin(1),parity,ngrid,rstep)
c
c     the direct term
c
      if(ttest) then
        sumdir = 0.0D0
        do ir=2,ngrid
          sumdir = wfin(ir)*cwfout(ir) + sumdir
        enddo
      endif
c
      cwfout(1) = d2wfin(1)
      do ir=2,ngrid
c        cwfout(ir) = -d2wfin(ir)+xla2*zentrf(ir)*wfin(ir)
c     &              + cwfout(ir)
        cwfout(ir) = (d2wfin(ir)-xla2*zentrf(ir)*wfin(ir))*h2mfac
     &              + cwfout(ir)
      enddo
c
      if(ttest) then
        do ir=2,ngrid
          sumdir = -wfin(ir)*cwfout(ir) + sumdir
        enddo
        sumdir = rstep*sumdir
      endif
      sumspe = sumdir                 ! c.m. correction to s.p.e.
      sumdel = sumdir                 ! c.m. correction to gap
c      write(6,'(a,1pg13.5,a)')  ' sumdir=',sumdir
c      write(6,'(5(1pg13.5))') (d2wfin(i),i=1,ngrid)
c
c     the exchange term - requiring summation over all states
c
c
c      the loop over the 'b'-angular momenta 'jsb'.
c        the derived quantities are:
c         'xjb'   = total angular momentum
c         'xlb'   = orbital angular momentum
c         'lb'    =    "       "       "      as integer
c
       do 40 jsb=1,jsmax(iso)
       nmaxb  = npla(1,jsb+isnpad)
       if(nmaxb.le.0) goto 40
       xjb    = float((jsb-1)/2)+0.5
       lb     = jsb/2
       xlb    = one*lb
c
c      check triangular relation for angular momentum selection rules
c
       if(mod(la+lb,2).ne.1) goto 40
       if(iabs(la-lb).gt.1) goto 40
       if(abs(xja-xjb).gt.one) goto 40
c
c      compute factors for reduced matrix element.
c
       if(iabs(la-lb).eq.1) then
        if(xja-xla.eq.xjb-xlb) then
         sixfcu = sqrt((xla+xjb+2.5)*(xla+xjb-0.5)
     *                 /((xla+xla+one)*(xlb+xlb+one)))
        else
         sixfcu = sqrt((xjb-xla+1.5)*(xla-xjb+1.5)
     *                 /((xla+xla+one)*(xlb+xlb+one)))
     *            *2.0*(xjb-xlb)
        endif
c        sixtst = sixjot(xlb,xjb,0.5D0,xja,xla,1.0D0)
c     &          *sqrt((xja+xja+1.0D0)*(xjb+xjb+1.0D0))*1.0D0
c     &          *(1.0D0-2.0D0*mod(xlb+xja+1.5D0,2.0D0))
c        write(6,'(a,4f5.1,2(1pg13.5))') ' test sixjot: l,j,lj=',
c     &       xla,xja,xlb,xjb,sixfcu,sixtst
        if(xla.gt.xlb) then
         facdru = sqrt(xla)
         facctu = -xla
         facctb = xlb
        else
         facdru = -sqrt(xlb)
         facctu = xlb
         facctb = -xla
        endif
       else
        sixfcu = zero
        facdru = zero
        facctu = zero
       endif
c
c      the transition matrix element
c
       do nshlb=1,nmaxb
         nptb   = npla(nshlb+1,jsb+isnpad)
         nptmb  = nptb-1
         v2b    = weight(nptb/ngrid+1)
         uvb    = sqrt(max(v2b-v2b*v2b,xsmall))
c  
         radmtu = zero
         do i=2,ngrid
           ra     = (i-1)*rstep
           radmtu = wf(i+nptmb)
     *                   *(-d1wfin(i)+facctu*wfin(i)/ra)
     *              + radmtu
         enddo
         radmtu = rstep*radmtu
         if(v2a.gt.xsmall) then
c         phfac= (1.0D0-2.0D0*mod(xlb+xja+1.5D0,2.0D0))
c     &         *(1.0D0-2.0D0*mod(xla+xjb+1.5D0,2.0D0))
c         wfac = -2.0D0*radmtu*(facdru*sixfcu)**2*wfac
           wfac = 2.0D0*radmtu*(facdru*sixfcu)**2*h2mfac
           wface = -wfac*v2b/(xja+xja+one)
           wfacd = -wfac*uvb/(xja+xja+one)
           wfac = wfac*(v2b+uva*uvb/max(v2a,v2lim))/(xja+xja+one)
           sume = 0.0D0
           sumd = 0.0D0
           do i=2,ngrid
             ra     = (i-1)*rstep
c           cwfout(i) = (-wfs(i+nptmb)+facctb*wf(i+nptmb)/ra)*wfac
c     &                 + cwfout(i)
             cwfb = (wfs(i+nptmb)+facctu*wf(i+nptmb)/ra)
             wfa  = wfin(i)
             cwfout(i) = cwfb*wfac  + cwfout(i)
             sume = (wfa*cwfb) + sume
             sumd = (wfa*cwfb) + sumd
           enddo
           sumspe = sumspe+rstep*wface*sume
           sumdel = sumdel+rstep*wfacd*sumd
         endif
c
       enddo
   40  continue
c
c      store c.m. correcrtions to s.p.e  and gap
c
       addspe = -sumspe
       adddel = -sumdel
c
c      finally compute overlap
c
       if(ttest) then
         sumex = 0.0D0
         do i=2,ngrid
           sumex = wfin(i)*cwfout(i) + sumex
         enddo
         sumex = rstep*sumex
         sumex  = sumex-sumdir
         phasjs = (xja+xja+one)
         zpedir = v2a*phasjs*sumdir + zpedir
         zpeex  = v2a*phasjs*sumex*0.5D0  + zpeex
         zpetot = zpedir+zpeex
c
         write(6,'(a,i5,4(1pg13.5))')
     &    ' CMACT: jsa,zpe,zpe_dir,zpe_ex,spe_dir=',jsa,
     &      zpetot,
     &      zpedir,
     &      zpeex,
     &      sumdir
       endif
c
       return
       end

c-----rvariance--------- part of the hartree-fock package ----------------------
c
      subroutine rvariance
c
      include 'params.inc'
c
c     * * * * * *    variance of r    * * * * * *
c
c     all necessary information is handled via common's.
c     input is:
c      wf       are the wave functions, eneterd via 'wavfun'
c     and further information on shell ordering, grid etc, comes
c     via the according common's (see there for further explanations).
c
c     output is printed directly within this routine
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c-----------------------------------------------------------------------
c
c     check whether 'zpe' is wanted
c
      if(ifzpe.lt.1 .and. ifzpe.gt.-3) then
       zpe    = zero
      else
c
c-------------------------------------------------------------------------------
c
c      full evaluation of momentum width <p**2> and subsequent classical
c      zpe correction.
c
c-------------------------------------------------------------------------------
c
c      isospin loop. the same isospin applies for the 'a'-state and
c      the 'b'-state, because the p-operator does not mix isospin.
c
       diag   = 0.0
       offdig = 0.0
       do 10 iso=1,2
       isnpad = (iso-1)*kjmx
c
c      'a'-loop over the angular momenta 'jsa'.
c        the derived quantities are:
c         'xja'   = total angular momentum
c         'xla'   = orbital angular momentum
c         'la'    =    "       "       "      as integer
c
       do 20 jsa=1,jsmax(iso)
       nmaxa  = npla(1,jsa+isnpad)
       if(nmaxa.le.0) goto 20
       xja    = float((jsa-1)/2)+0.5
       la     = jsa/2
       xla    = one*la
       xla2   = xla*xla+xla
c
c      the diagonal contribution from p**2.
c
       do 30 nshla=1,nmaxa
       npta   = npla(nshla+1,jsa+isnpad)
       nptma  = npta-1
       sum    = 0.0
       do i=2,ngrid
         wfa = wf(i+nptma)*((i-1)*rstep)
c        sum    = ( (wfs(i+nptma)*wfs(i+nptma))
         sum    = wfa*wfa  + sum
       enddo
       sum    = rstep*sum
       v2a    = weight(npta/ngrid+1)
       if(temprt.gt.zero) then
        ttana  = ttthrm(npta/ngrid+1)
        v2a    = (v2a-0.5*(one-ttana))/ttana
        v2a    = (one+ttan)*0.5*v2a
c         write(6,'(a,3(1pg12.4))')
c     &    ' v2a,newfac,thfac=',weight(npta/ngrid+1),v2a,
c     &    (one+ttthrm(npta/ngrid+1))*0.5         
       endif
       diag   = v2a*(xja+xja+one)*sum + diag
   30  continue
c
c      the loop over the 'b'-angular momenta 'jsb'.
c        the derived quantities are:
c         'xjb'   = total angular momentum
c         'xlb'   = orbital angular momentum
c         'lb'    =    "       "       "      as integer
c
       do 40 jsb=1,jsa
       nmaxb  = npla(1,jsb+isnpad)
       if(nmaxb.le.0) goto 40
       xjb    = float((jsb-1)/2)+0.5
       lb     = jsb/2
       xlb    = one*lb
c
c      check triangular relation for angular momentum selection rules
c
       if(mod(la+lb,2).ne.1) goto 40
       if(iabs(la-lb).gt.1) goto 40
       if(abs(xja-xjb).gt.one) goto 40
c
c      compute factors for reduced matrix element.
c
       if(iabs(la-lb).eq.1) then
        if(xja-xla.eq.xjb-xlb) then
         sixfcu = sqrt((xla+xjb+2.5)*(xla+xjb-0.5)
     *                 /((xla+xla+one)*(xlb+xlb+one)))
        else
         sixfcu = sqrt((xjb-xla+1.5)*(xla-xjb+1.5)
     *                 /((xla+xla+one)*(xlb+xlb+one)))
     *            *2.0*(xjb-xlb)
        endif
        if(xla.gt.xlb) then
         facdru = sqrt(xla)
         facctu = -xla
        else
         facdru = -sqrt(xlb)
         facctu = xlb
        endif
       else
        sixfcu = zero
        facdru = zero
        facctu = zero
       endif
c
c      the radial matrix elements for given jsa and jsb, sums running
c      over all nshla (to jsa) and nshlb (to jsb).
c
       do 50 nshla=1,nmaxa
       npta   = npla(nshla+1,jsa+isnpad)
       nptma  = npta-1
       v2a    = weight(npta/ngrid+1)
       if(temprt.gt.zero) then
        ttana  = ttthrm(npta/ngrid+1)
        v2a    = (v2a-0.5*(one-ttana))/ttana
c        v2a    = (one+ttan)*0.5*v2a
       endif
       uva    = sqrt(max(v2a-v2a*v2a,xsmall))
c
       do 50 nshlb=1,nmaxb
       nptb   = npla(nshlb+1,jsb+isnpad)
       nptmb  = nptb-1
       v2b    = weight(nptb/ngrid+1)
       if(temprt.gt.zero) then
        ttanb  = ttthrm(nptb/ngrid+1)
        v2b    = (v2b-0.5*(one-ttanb))/ttanb
c        v2a    = (one+ttan)*0.5*v2a
       endif
       uvb    = sqrt(max(v2b-v2b*v2b,xsmall))
c
       radmtu = zero
       do i=2,ngrid
         ra     = (i-1)*rstep
         radmtu = wf(i+nptmb)*wf(i+nptma)*ra + radmtu
       enddo
       radmtu = rstep*radmtu
       if(temprt.gt.zero) then
         wtha  = (one-ttthrm(npta/ngrid+1))*0.5
         wthb  = (one-ttthrm(nptb/ngrid+1))*0.5
         wthca = one-wtha
         wthcb = one-wthb
         wfac = wtha*wthb*(-0.5*(v2a+v2b)+(v2a*v2b-uva*uvb))
     &          +wthb*wthca*v2a
     &          +wthca*wthcb*(v2a*v2b-uva*uvb)
       else
         wfac = (v2a*v2b-uva*uvb)
       endif
       offdig = (radmtu*facdru*sixfcu)**2*wfac
     *         + offdig
c
   50  continue
   40  continue
   20  continue
   10  continue
c
c      finally assembling the zpe-correction  (in MeV).
c
       diag = diag/(anumb*anumb)
       offdig = offdig/(anumb*anumb)
       rwidth = (diag-offdig-offdig)
       peffec = 2.25D0/rwidth
       zper   = peffec*half*(amfc(1)+amfc(2))/anumb
c
c      optionally comparision with simple estimates for 'zpe'
c
       if(iwr.ge.-1) then
        write(06,'(5(a,1pg12.5))') 
     &                ' zpe_r=',zper,
     &                ' rwidth=',rwidth,
     &                ' diag=',diag,
     &                ' offdig=',offdig,
     &                ' redun=',4.0D0*pwidth*rwidth/9.0D0
c
       endif
c
      endif
c
      return
      end

c-----structure_factor--------- part of the hartree-fock package ----------------------
c
      subroutine structure_factor
c
c     this routine is not yet ready for use
c
      include 'params.inc'
c
c     * * * * * *    structure factor and correlation function    * * * * * *
c
c     all necessary information is handled via common's.
c     input is:
c      wf       are the wave functions, eneterd via 'wavfun'
c     and further information on shell ordering, grid etc, comes
c     via the according common's (see there for further explanations).
c
c     output is printed directly within this routine
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'constn.inc'
      REAL(8) :: strucfac(kqmax,2),correlfunc(kqmax,2)
c
c-----------------------------------------------------------------------
c
c
c-------------------------------------------------------------------------------
c
c      full evaluation of momentum width <p**2> and subsequent classical
c      zpe correction.
c
c-------------------------------------------------------------------------------
c
c     q-space parameters
c
      iqmax  = min(kfrm,ngrid*3/4)
      rleng  = (ngrid-1)*rstep
      delq   = 3.1415926/rleng
      frpirs = 4*pi*rstep
      strucfac = 0D0
c
c      isospin loop. the same isospin applies for the 'a'-state and
c      the 'b'-state, because the p-operator does not mix isospin.
c
       do 10 iso=1,2
       isnpad = (iso-1)*kjmx
c
c      'a'-loop over the angular momenta 'jsa'.
c        the derived quantities are:
c         'xja'   = total angular momentum
c         'xla'   = orbital angular momentum
c         'la'    =    "       "       "      as integer
c
       do 20 jsa=1,jsmax(iso)
       nmaxa  = npla(1,jsa+isnpad)
       if(nmaxa.le.0) goto 20
       xja    = float((jsa-1)/2)+0.5
       la     = jsa/2
       xla    = one*la
       xla2   = xla*xla+xla
c
c      the loop over the 'b'-angular momenta 'jsb'.
c        the derived quantities are:
c         'xjb'   = total angular momentum
c         'xlb'   = orbital angular momentum
c         'lb'    =    "       "       "      as integer
c
       do 40 jsb=1,jsa
       nmaxb  = npla(1,jsb+isnpad)
       if(nmaxb.le.0) goto 40
       xjb    = float((jsb-1)/2)+0.5
       lb     = jsb/2
       xlb    = one*lb
c
c      check triangular relation for angular momentum selection rules
c
       if(mod(la+lb,2).ne.1) goto 40
       if(iabs(la-lb).gt.1) goto 40
       if(abs(xja-xjb).gt.one) goto 40
c
c      compute factors for reduced matrix element.
c
       if(iabs(la-lb).eq.1) then
        if(xja-xla.eq.xjb-xlb) then
         sixfcu = sqrt((xla+xjb+2.5)*(xla+xjb-0.5)
     *                 /((xla+xla+one)*(xlb+xlb+one)))
        else
         sixfcu = sqrt((xjb-xla+1.5)*(xla-xjb+1.5)
     *                 /((xla+xla+one)*(xlb+xlb+one)))
     *            *2.0*(xjb-xlb)
        endif
        if(xla.gt.xlb) then
         facdru = sqrt(xla)
         facctu = -xla
        else
         facdru = -sqrt(xlb)
         facctu = xlb
        endif
       else
        sixfcu = zero
        facdru = zero
        facctu = zero
       endif
c
c      the radial matrix elements for given jsa and jsb, sums running
c      over all nshla (to jsa) and nshlb (to jsb).
c
       do 50 nshla=1,nmaxa
       npta   = npla(nshla+1,jsa+isnpad)
       nptma  = npta-1
       v2a    = weight(npta/ngrid+1)
       if(temprt.gt.zero) then
        ttana  = ttthrm(npta/ngrid+1)
        v2a    = (v2a-0.5*(one-ttana))/ttana
c        v2a    = (one+ttan)*0.5*v2a
       endif
       uva    = sqrt(max(v2a-v2a*v2a,xsmall))
c
       do 50 nshlb=1,nmaxb
       nptb   = npla(nshlb+1,jsb+isnpad)
       nptmb  = nptb-1
       v2b    = weight(nptb/ngrid+1)
       if(temprt.gt.zero) then
        ttanb  = ttthrm(nptb/ngrid+1)
        v2b    = (v2b-0.5*(one-ttanb))/ttanb
c        v2a    = (one+ttan)*0.5*v2a
       endif
       uvb    = sqrt(max(v2b-v2b*v2b,xsmall))
c
       do iq=1,iqmax
         q1     = iq*delq-delq
         radmtu = zero
         do i=2,ngrid
           ra     = (i-1)*rstep
           radmtu = wf(i+nptmb)*wf(i+nptma)*bessel(q1*ra,0) + radmtu
         enddo
         radmtu = rstep*radmtu
         if(temprt.gt.zero) then
           wtha  = (one-ttthrm(npta/ngrid+1))*0.5
           wthb  = (one-ttthrm(nptb/ngrid+1))*0.5
           wthca = one-wtha
           wthcb = one-wthb
           wfac = wtha*wthb*(-0.5*(v2a+v2b)+(v2a*v2b-uva*uvb))
     &            +wthb*wthca*v2a
     &            +wthca*wthcb*(v2a*v2b-uva*uvb)
         else
           wfac = (v2a*v2b-uva*uvb)
         endif
         strucfac(iq,iso) = -(radmtu*facdru*sixfcu)**2*wfac
     *                     + strucfac(iq,iso)
       enddo
c
   50  continue
   40  continue
   20  continue
   10  continue
       strucfac(1:iqmax,1) = strucfac(1:iqmax,1)/aprot
       strucfac(1:iqmax,2) = strucfac(1:iqmax,2)/aneut
c
c      printing structure factor
c
!       strucfac(1,:) = 0D0
       OPEN(922,file='structure_factor.res')
       WRITE(922,'(a/a)') '# structure factor',
     &                    '#  q     proton   neutron '
       do iq=1,iqmax
         q1     = iq*delq-delq
         WRITE(922,'(f8.4,2(1pg13.5))') q1,strucfac(iq,:)
       enddo
c
c     correlation function
      dpi3   = half/(pi*pi)
      DO iso=1,2
c                                                         r=0.0 separately
      acc    = zero
      do iq=2,iqmax
        qa     = delq*iq-delq
        acc    = qa*qa*strucfac(iq,iso)+acc
      enddo
      correlfunc(1,iso) = abs(acc*delq*dpi3)
c                                                         now all other r
      r      = zero
      do ir=2,ngrid
        r      = rstep+r
        acc    = zero
        qa     = zero
        do iq=2,iqmax
          qa     = delq+qa
          acc    = qa*qa*bessel(qa,0)*strucfac(iq,iso)+acc
        enddo
        correlfunc(ir,iso) = abs(acc*delq*dpi3/r)
c
      enddo
      END DO
c
c      printing correlation function
c
       WRITE(922,'(1x/1x/a/a)') '# correlation function',
     &                    '#  r     proton   neutron '
       do i=1,iqmax
         WRITE(922,'(f8.4,2(1pg13.5))') rstep*(i-1),correlfunc(i,:)
       enddo
c
      return
      end
