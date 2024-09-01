c-----gscat --------------------------------------------------------------------
c
      subroutine gscat(dgsrho,dhhrho,qqrpa,qqhf,pprpa,pphf,
     &                 maxmod,nameq,namep,ichr,resfac,cmfac,ifppsc)
c
      include 'params.inc'
      parameter(krpmat=krpmax*(krpmax+1)/2)
c
c     computes contribution of the scattering diagrams to the gsc.
c     uses
c
c               (sc)
c      del^2-obs     = (2<p^2>-2<p^2>_0) <robs[v_q,q]>
c
c                     +(2<q^2>-2<q^2>_0) <robs[v_p,p]>
c
c     where
c
c       [h_0-epsilon_alpha,robs] = obs.
c
c     this version computes the effect on the full spatial density.
c     the densities are accumulated on tye workspaces 'dgsrho' and
c     'dhhrho'.
c
c     redundance is usually guaranteed because both, qqrpa and
c     qqrpa0=qqhf, are generated from an rpa-recipe.
c
c
      character*6 nameq,namep
      character*2 ichr(kichr)
      dimension   qqrpa(krpmat),qqhf(krpmat),
     &            pprpa(krpmat),pphf(krpmat)
      dimension   dgsrho(kgr2),dhhrho(kgr2)
c
      include 'wavfun.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'result.inc'
      include 'gscres.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
c     set counters
c
      do 1010 ir=1,ngrid2
      dgsrho(ir) = zero
      dhhrho(ir) = zero
 1010 continue
c
c     first part: effect of q-response
c
c     loop over the initial modes:
c      reading response wave functions from buffer file 18
c      preparing response hamiltonian
c
      do 210 imode=1,maxmod
      write(6,'(a,i5)')
     &  ' q-response in gscat for mode nr.',imode
      nameq(4:4) = ichr(imode)
      call buffql(qlrout,'i',nameq)
c
c     the response hamiltonian for 'imode'
c
      call qcurre(resfac)
      call qpot(resfac)
c
c     now the second mode 'imod2'
c
      do 220 imod2=1,maxmod
      if(imode.le.imod2) then
       ipoint = (imod2*(imod2-1))/2+imode
      else
       ipoint = (imode*(imode-1))/2+imod2
      endif
      widfac = 2.0*(pprpa(ipoint)-pphf(ipoint))
cc      if(maxmod.ge.1)
cc     &   write(14,'(a,3g12.4)') ' widfac: pprpa,pphf,diff=',
cc     &   pprpa(ipoint),pphf(ipoint),widfac
c
c     read second response momentum
c
      nameq(4:4) = ichr(imod2)
      call buffql(qlrout,'i',nameq)
c
c     accumulate
c
c                                                      hh-parts
      call qqcomm(resfac,widfac,2)
      call mfresp(dhhrho)
c                                                      pp-parts
      call qqcomm(resfac,widfac,1)
      call mfresp(dgsrho)
c
  220 continue
  210 continue
c
c     second part: effect of p-response
c
c     loop over the initial modes:
c      reading response wave functions from buffer file 18
c      preparing response hamiltonian
c
      do 310 imode=1,maxmod
      namep(2:2) = ichr(imode)
      call buffql(plrout,'i',namep)
      write(6,'(a,i5)')
     &  ' p-response in gscat for mode nr.',imode
c
c     the response hamiltonian for 'imode'
c
      call pdensi(resfac)
      if(iwr.ge.00) write(6,'(/a,30(/1x,5g12.4))')
     & ' prho in gscat:',(prho(ir),ir=1,ngrid2)
      call ppot(resfac,cmfac)
c
c     now the second mode 'imod2'
c
      do 320 imod2=1,maxmod
      if(imode.le.imod2) then
       ipoint = (imod2*(imod2-1))/2+imode
      else
       ipoint = (imode*(imode-1))/2+imod2
      endif
      widfac = 2.0*(qqrpa(ipoint)-qqhf(ipoint))
cc      if(maxmod.eq.1)
cc     &   write(14,'(a,3g12.4)') ' widfac: qqrpa,qqhf,diff=',
cc     &   qqrpa(ipoint),qqhf(ipoint),widfac
c
c     read second response momentum
c
      namep(2:2) = ichr(imod2)
      call buffql(qlrout,'i',namep)
c
c     accumulate
c
c                                                      hh-parts
      call ppcomm(resfac,widfac,2)
      call mfresp(dhhrho)
c                                                      pp-parts
      call ppcomm(resfac,widfac,1)
      call mfresp(dgsrho)
c
  320 continue
  310 continue
      parll  = one
      call extrap(dgsrho(1),parll)
      call extrap(dgsrho(ngrid+1),parll)
      call extrap(dhhrho(1),parll)
      call extrap(dhhrho(ngrid+1),parll)
c
c     change sign of the correlation densities
c      (this if an ad hoc measure to obtain agreement with
c       estimates in the schematic model.)
c
       do 330 ir=1,ngrid
       ira = ir+ngrid
       dgsrho(ir) = -dgsrho(ir)
       dgsrho(ira) = -dgsrho(ira)
       dhhrho(ir) = -dhhrho(ir)
       dhhrho(ira) = -dhhrho(ira)
  330  continue
c
c     evaluate global observables and print result
c
      scpppr  = zero
      scppne  = zero
      schhpr  = zero
      schhne  = zero
      ifprin = 0
      if(ifprin.ne.0) write(14,'(/a/a,a/3a)')
     &  ' scattering part of correlation density:',
     &  '     pp-part                hh-part    ',
     &  '    total        ',
     &  (' proton     neutron  ',ir=1,3)
      do 110 ir=1,ngrid
      r   = (ir-1)*rstep
      r4  = r*r
      r4  = r4*r4
      ira = ngrid+ir
      work1(ir) = dgsrho(ir)+dhhrho(ir)
      work2(ir) = dgsrho(ira)+dhhrho(ira)
      scpppr  = dgsrho(ir)*r4 + scpppr
      scppne  = dgsrho(ira)*r4 + scppne
      schhpr  = dhhrho(ir)*r4 + schhpr
      schhne  = dhhrho(ira)*r4 + schhne
  110 if(ifprin.ne.0) write(14,'(1x,f5.2,8g12.3)')
     &  r,dgsrho(ir),dgsrho(ira),dhhrho(ir),dhhrho(ira),
     &  work1(ir),work2(ir)
c
c     formfactor
c
      facp = one/aprot
      if(aneut.ne.0) facn = one/aneut
      qsurf  = 5.600/rdms
      qrdms   = 4.493/rdms
      if(ifprin.ne.0) then
       write(14,'(/2(a,f8.4)/2a)')
     &  ' formfactors:  qsurf=',qsurf,'  qrdms=',qrdms,
     &  '  iq   q      proton     neutron'
c                                                      all  q
       delqh = 0.5*delq
       do 530 iq=1,ngrid
       q1     = iq*delqh-delqh
       dformp = facp*fbform(q1,work1)
       dformn = facn*fbform(q1,work2)
       if(iq.le.ngrid .and. ifprin.ne.0)
     &   write(14,'(1x,i3,f6.2,2(2x,6g11.3))')
     &   iq,q1,dformp,dformn
  530  continue
      endif
c
c     print r.m.s. radii
c
      scpppr  = frpirs*scpppr/(aprot*3.0)
      if(aneut.ne.0) scppne  = frpirs*scppne/(aneut*3.0)
      schhpr  = frpirs*schhpr/(aprot*3.0)
      if(aneut.ne.0) schhne  = frpirs*schhne/(aneut*3.0)
      write(14,'(/a,3(/5x,a,2g12.4))')
     & ' delta sigma for scat.diagrams(p/n):','pp-part=',scpppr,scppne,
     & 'hh-part=',schhpr,schhne,' total=',scpppr+schhpr,scppne+schhne
      delsig(1,1,2,irunc) = scpppr
      delsig(2,1,2,irunc) = scppne
      delsig(1,2,2,irunc) = schhpr
      delsig(2,2,2,irunc) = schhne
c
      call formcr(dgsrho(1),aprot,rdmsp,surfp,
     &            rdmscr(1,1,2,irunc),surfcr(1,1,2,irunc))
      call formcr(dhhrho(1),aprot,rdmsp,surfp,
     &            rdmscr(1,2,2,irunc),surfcr(1,2,2,irunc))
      if(aneut.ne.0) then
       call formcr(dgsrho(ngrid+1),aneut,rdmsn,surfn,
     &             rdmscr(2,1,2,irunc),surfcr(2,1,2,irunc))
       call formcr(dhhrho(ngrid+1),aneut,rdmsn,surfn,
     &             rdmscr(2,2,2,irunc),surfcr(2,2,2,irunc))
      endif
      write(14,'(a,2g12.4)')
     &   ' delta_sigma from formf.: p/n=',
     &  surfcr(1,1,2,irunc)+surfcr(1,2,2,irunc),
     &  surfcr(2,1,2,irunc)+surfcr(2,2,2,irunc)
      write(14,'(a,2g12.4)')
     &   '       delta_diff.radius: p/n=',
     &  rdmscr(1,1,2,irunc)+rdmscr(1,2,2,irunc),
     &  rdmscr(2,1,2,irunc)+rdmscr(2,2,2,irunc)
c
      return
      end
c-----qqcomm-------------------------------------------------------------------
c
      subroutine qqcomm(resfac,widfac,ifppsc)
c
      include 'params.inc'
c
c     repsonse hamiltonian of q-type acting on another q-type
c     wavefunction. result is accumulted on 'savlr' as workspace.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'result.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
c      reset accumulator field 'savlr'
c
      do 1010 ir=1,kx11
 1010 savlr(ir,klrpa) = zero
c
      include 'initdo.inc'
c
c     action of linear response hamiltonian
c
      iph    = xj-xja
      fac = resfac*widfac*(1-2*mod(abs(iph),2))
      jdifb = xja-xj+lrpa+1
c
c     accumulate to response density
c
       if(ifppsc.ne.2) then
        call multr(qlrout(nadd+1,jdif),work1(1),rstep,ngrid)
        sig    = float(mod(la,2)*2-1)
        call derv15(work1(1),aux1(1),sig,ngrid,rstep)
        if(resfac.ne.zero)
     &   call hrespq(iso,ll,xj,nadd,xja,xla,jdifb,
     &               work1(1),aux1(1),hlinrs)
cpro        call projec(hlinrs(1),jsa,nstack,nmax)
cc        fac = resfac*widfac
        do 41 ir=2,ngrid
        savlr(nadd+ir,klrpa) =  fac*hlinrs(ir)
     &                        + savlr(nadd+ir,klrpa)
   41   continue
       endif
c
c     subtract commuted part
c
      if(ifppsc.ne.1 .and. resfac.ne.zero) then
        call hrespq(iso,la,xja,nadd,xj,xl,jdif,
     &              wf(nadd+1),wfs(nadd+1),hlinrs)
       do 100 nb=1,nmax
       jsb   = jsp(nb)
       wgtb  = weight(nb+nstack/ngrid)
       if(js.eq.jsb .and. wgtb.gt.half) then
        lb    = jsb/2
        xlb   = lb
        xjb   = (jsb-1)/2+half
        naddb = nstack + (nb-1)*ngrid
c
        sum   = zero
        do 110 ir=2,ngrid
        r     = (ir-1)*rstep
  110   sum   = r*hlinrs(ir)*wf(naddb+ir) + sum
        sum   = -rstep*sum*fac
cc        write(14,'(a,2i4,a,g12.3)')
cc     &  ' <wf|v|wf> in qqcomm for na,nb=',na,nb,' is:',sum
c
        do 120 ir=2,ngrid
  120   savlr(nadd+ir,klrpa) = sum*qlrout(naddb+ir,jdifb)
     &                        + savlr(nadd+ir,klrpa)
c
       endif
c
  100  continue
c
       endif
c
      include 'enddo.inc'
c
      return
      end
c-----ppcomm-------------------------------------------------------------------
c
      subroutine ppcomm(resfac,widfac,ifppsc)
c
      include 'params.inc'
c
c     repsonse hamiltonian of p-type acting on another p-type
c     wavefunction. result is accumulted on 'savlr' as workspace.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'result.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
c      reset accumulator field 'savlr'
c
      do 1010 ir=1,kx11
 1010 savlr(ir,klrpa) = zero
c
      include 'initdo.inc'
c
c     action of linear response hamiltonian
c
      jdifb = xja-xj+lrpa+1
      iph    = xj-xja
      fac = resfac*widfac*(1-2*mod(abs(iph),2))
c
c     accumulate
c
       if(ifppsc.ne.2) then
        call multr(qlrout(nadd+1,jdif),work1(1),rstep,ngrid)
        sig    = float(mod(la,2)*2-1)
        call derv15(work1(1),aux1(1),sig,ngrid,rstep)
        if(resfac.ne.zero)
     &     call hrespp(iso,ll,xj,nadd,xja,xla,jdifb,
     &               work1(1),aux1(1))
cpro        call projec(hlinrs(1),jsa,nstack,nmax)
cc         fac = resfac*widfac
        do 41 ir=2,ngrid
        savlr(nadd+ir,klrpa) =  fac*hlinrs(ir)
     &                        + savlr(nadd+ir,klrpa)
   41   continue
       endif
c
c     subtract commuted part
c
      if(ifppsc.ne.1 .and. resfac.ne.zero) then
         call hrespp(iso,la,xja,nadd,xj,xl,jdif,
     &               wf(nadd+1),wfs(nadd+1))
       do 100 nb=1,nmax
       jsb   = jsp(nb)
       wgtb  = weight(nb+nstack/ngrid)
       if(js.eq.jsb .and. wgtb.gt.half) then
        lb    = jsb/2
        xlb   = lb
        xjb   = (jsb-1)/2+half
        naddb = nstack + (nb-1)*ngrid
c
        sum   = zero
        do 110 ir=2,ngrid
        r     = (ir-1)*rstep
  110   sum   = r*hlinrs(ir)*wf(naddb+ir) + sum
        sum   = -rstep*sum*fac
cc        write(14,'(a,2i4,a,g12.3)')
cc     &  ' <wf|v|wf> in ppcomm for na,nb=',na,nb,' is:',sum
c
        do 120 ir=2,ngrid
  120   savlr(nadd+ir,klrpa) = sum*qlrout(naddb+ir,jdifb)
     &                        + savlr(nadd+ir,klrpa)
       endif
c
  100  continue
      endif
c
      include 'enddo.inc'
c
       return
       end
c-----mfresp-------------------------------------------------------------------
c
      subroutine mfresp(rsprho)
c
      include 'params.inc'
c
c     mean-field repsonse on the wavefunctions which are stored
c     on the workspace 'savlr(*,klrpa)'. result is accumulated
c     to a density on 'rsprho'.
c
      dimension rsprho(kgr2)
c
      include 'wavfun.inc'
      include 'potinv.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'result.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
c     set accumulator for response density to zero
c
c      do 1020 ir=1,ngrid2
c 1020 rsprho(ir) = zero
c
c    set multiplicity and normalization factor
c
      accfac = mlrpa*d4pi
c
c     iso-loop
c      'isoadd' is a pointer for the potential and mass fields;
c      'isnpad' is a pointer for the wave function pointer 'npla'.
c
      do 1000 iso=1,2
      isoadd = (iso-1)*ngrid
      isnpad = (iso-1)*kjmx
      nstack = (iso-1)*nmaxp*ngrid
      if(iso.eq.1) then
       nmaxal = nmaxp
      else
       nmaxal = nmaxn
      endif
      if(nmaxal.eq.0) goto 1000
c
c     square root of inverse mass 'bmass' on 'sqb';
c     d/dr**2 'sqb' / 'sqb' on 'd2b';
c     d/dr('bmass')/r on 'd1b'.
c
      do 10 i=1,ngrid
   10 sqb(i) = sqrt(bmass(i+isoadd))
c
      call derv25(sqb(1),d2b(1),one,ngrid,rstep)
      do 20 i=1,ngrid
   20 d2b(i) = d2b(i)/sqb(i)
c
      call derv15(bmass(isoadd+1),d1b(1),one,ngrid,rstep)
      d1b(1)=(16.0*bmass(2+isoadd)-15.0*bmass(1+isoadd)
     *        -bmass(3+isoadd))/(6.0*rstep**rstep)
      do 30 i=2,ngrid
   30 d1b(i) = d1b(i)*rinvs(i)
c
c     loop over composed ang.mom. 'js'
c      'll' is the orbital ang.mom.;
c      'sig' is the parity of the wave function 'wf' (it is -parity
c      of the state because 'wf' is multiplied by r);.
c      'phasjs' is the degeneracy of the shell.
c      'jsnpl' is pointer to to nucleon shell in pointer field 'npla'.
c
      do 40 js=1,jsmax(iso)
      jsnpl = js+isnpad
      nmax  = npla(1,jsnpl)
      if(nmax.gt.0) then
       ll     = js/2
       sig    = float(mod(ll,2)*2-1)
       xlp    = float(ll*ll+ll)
       xl0    = wls(js)
       phasjs = float(js-ll)
c
c      the effective potential for given 'js' and 'iso' on 'potl';
c      its minimum on 'umin';
c      the zentrifugal force on 'aux3'.
c
       potl(1) = pot(1+isoadd)+xl0*potls(1+isoadd)+d1b(1)
       do 60 i=2,ngrid
       n      = i+isoadd
       y      = xlp*zentrf(i)
       y       = pot(n)+xl0*potls(n)+bmass(n)*y+d1b(i)
       potl(i) = y
   60  continue
c
c      loop over radial quantum number 'n'.
c       accumulating energies and performing one w.f. iteration.
c       'npl' is wave function pointer,
c       'nplshl' is the number of the state in harmonic oscillator ordering.
c
       do 70 n=1,nmax
       npl    = npla(n+1,jsnpl)
       nplm   = npl-1
       nplshl = npl/ngrid+1
c
c      project and prepare wavefunction
c
       call projec(savlr(npl,klrpa),js,nstack,nmaxal)
       call multr(savlr(npl,klrpa),aux1(1),rstep,ngrid)
c
c      inversion to obtain response
c
       irep   = -1
       umin   = spenrg(nplshl)
       call invert(aux1,aux1,irep)
cc       call projec(aux1(1),js,nstack,nmaxal)
c
c      accumulate to density properly subtracting spurios components
c
       aux1(1) = zero
       do 310 ir=2,ngrid
  310  aux1(ir) = rinvs(ir)*aux1(ir)
       call projec(aux1(1),js,nstack,nmaxal)
cc       ovproj = zero
cc       do 310 ir=1,ngrid
cc       ovproj = wf(nplm+ir)*aux1(ir) + ovproj
cc  310  continue
cc       ovproj = rstep*ovproj
cc       sumtst = zero
       do 320 ir=2,ngrid
cc       aux1(ir) = aux1(ir) - ovproj*wf(nplm+ir)
cc       sumtst = wf(nplm+ir)*aux1(ir) + sumtst
       ri = rinvs(ir)
  320  rsprho(isoadd+ir) = accfac*ri*wf(nplm+ir)*aux1(ir)
     &                    + rsprho(isoadd+ir)
cc       write(14,'(a,i4,a,g12.4)')
cc     &  ' test orton. for na=',nplshl,': error=',sumtst*frpirs
c
   70  continue
      endif
c
   40 continue
 1000 continue
c
c     restore 'd1b'
c
      call derv17(bmass(1),d1b(1),one,ngrid,rstep)
      call derv17(bmass(ngrid+1),d1b(ngrid+1),one,ngrid,rstep)
c
       return
       end
ctestc-----gscat2--------------------------------------------------------------------
ctestc
ctest      subroutine gscat2(qqrpa,qqhf,pprpa,pphf,
ctest     &                  maxmod,nameq,namep,ichr,resfac,cmfac)
ctestc
ctest      include 'params.inc'
ctest      parameter(krpmat=krpmax*(krpmax+1)/2)
ctestc
ctestc     computes contribution of the scattering diagrams to the gsc.
ctestc     uses
ctestc
ctestc               (sc)
ctestc      del^2-obs     = (2<p^2>-2<p^2>_0) <q[v_q,robs]>
ctestc
ctestc                     +(2<q^2>-2<q^2>_0) <p[v_p,robs]>
ctestc
ctestc     where
ctestc
ctestc       [h_0-epsilon_alpha,robs] = obs.
ctestc
ctestc     presently we run obs=rms-radius^2.
ctestc
ctestc     this program serves as a countercheck of the routine 'gscat'.
ctestc
ctestc     redundance is usually guaranteed because both, qqrpa and
ctestc     qqrpa0=qqhf, are generated from an rpa-recipe.
ctestc
ctestc
ctest      character*6 nameq,namep
ctest      character*1 ichr(62)
ctest      dimension   qqrpa(krpmat),qqhf(krpmat),
ctest     &            pprpa(krpmat),pphf(krpmat)
ctestc
ctest      include 'wavfun.inc'
ctest      include 'numpar.inc'
ctest      include 'grid.inc'
ctest      include 'linres.inc'
ctest      include 'result.inc'
ctest      include 'constn.inc'
ctestc
ctestc------------------------------------------------------------------------------
ctestc
ctestc     prepare 'robs' auxiliary storage 'savlr'.
ctestc
ctest      call rspobs
ctestc
ctestc     set counters
ctestc
ctest      rmsctq = zero
ctest      rmsctp = zero
ctestc
ctestc     first part: effect of q-response
ctestc
ctestc     loop over the initial modes:
ctestc      reading response wave functions from buffer file 18
ctestc      preparing response hamiltonian
ctestc
ctest      do 210 imode=1,maxmod
ctest      nameq(4:4) = ichr(imode)
ctest      call buffql(qlrout,'i',nameq)
ctestc
ctestc     the response hamiltonian for 'imode'
ctestc
ctest      call qcurre(resfac)
ctest      call qpot(resfac)
ctest      call qrobs(resfac)
ctestc
ctestc     now the second mode 'imod2'
ctestc
ctest      do 220 imod2=1,maxmod
ctest      if(imode.le.imod2) then
ctest       ipoint = (imod2*(imod2-1))/2+imode
ctest      else
ctest       ipoint = (imode*(imode-1))/2+imod2
ctest      endif
ctest      widfac = 2.0*(pprpa(ipoint)-pphf(ipoint))
ctest      if(maxmod.eq.1)
ctest     &   write(14,'(a,3g12.4)') ' widfac: pprpa,pphf,diff=',
ctest     &   pprpa(ipoint),pphf(ipoint),widfac
ctestc
ctestc     read second response momentum
ctestc
ctest      nameq(4:4) = ichr(imod2)
ctest      call buffql(qlrout,'i',nameq)
ctestc
ctestc     accumulate
ctestc
ctest      rmsctq = widfac*qpover(qlrout,plrout) + rmsctq
ctestc
ctest  220 continue
ctest  210 continue
ctestc
ctest      write(6,'(a,g12.4,i5)') ' resfac,ifsepi=',resfac,ifsepi
ctestc
ctestc     second part: effect of p-response
ctestc
ctestc     loop over the initial modes:
ctestc      reading response wave functions from buffer file 18
ctestc      preparing response hamiltonian
ctestc
ctest      do 310 imode=1,maxmod
ctest      namep(2:2) = ichr(imode)
ctest      call buffql(plrout,'i',namep)
ctestc
ctestc     the response hamiltonian for 'imode'
ctestc
ctest      call pdensi(resfac)
ctest      if(iwr.ge.00) write(6,'(/a,30(/1x,5g12.4))')
ctest     & ' prho in gscat:',(prho(ir),ir=1,ngrid2)
ctest      call ppot(resfac,cmfac)
ctest      call probs(resfac)
ctestc
ctestc     now the second mode 'imod2'
ctestc
ctest      do 320 imod2=1,maxmod
ctest      if(imode.le.imod2) then
ctest       ipoint = (imod2*(imod2-1))/2+imode
ctest      else
ctest       ipoint = (imode*(imode-1))/2+imod2
ctest      endif
ctest      widfac = 2.0*(qqrpa(ipoint)-qqhf(ipoint))
ctest      if(maxmod.eq.1)
ctest     &   write(14,'(a,3g12.4)') ' widfac: qqrpa,qqhf,diff=',
ctest     &   qqrpa(ipoint),qqhf(ipoint),widfac
ctestc
ctestc     read second response momentum
ctestc
ctest      namep(2:2) = ichr(imod2)
ctest      call buffql(qlrout,'i',namep)
ctestc
ctestc     accumulate
ctestc
ctest      rmsctp = widfac*qpover(qlrout,plrout) + rmsctp
ctestc
ctest  320 continue
ctest  310 continue
ctestc
ctestc     print result
ctestc
ctest      rmsctp = mlrpa*rmsctp/(3.0*anumb)
ctest      rmsctq = mlrpa*rmsctq/(3.0*anumb)
ctest      write(14,'(/a/1x,3(a,g12.3))')
ctest     & ' delta sigma for scat.diagrams (total):',' q-part=',rmsctq,
ctest     & ' p-part=',rmsctp,' total=',rmsctp+rmsctq
ctestc
ctest      return
ctest      end
ctestc-----rspobs-------------------------------------------------------------------
ctestc
ctest      subroutine rspobs
ctestc
ctest      include 'params.inc'
ctestc
ctestc     mean-field repsonse on the observable 'rms-radius'
ctestc     stored on 'savlr'.
ctestc
ctest      include 'wavfun.inc'
ctest      include 'potinv.inc'
ctest      include 'grid.inc'
ctest      include 'linres.inc'
ctest      include 'result.inc'
ctest      include 'constn.inc'
ctestc
ctestc------------------------------------------------------------------------------
ctestc
ctestc     iso-loop
ctestc      'isoadd' is a pointer for the potential and mass fields;
ctestc      'isnpad' is a pointer for the wave function pointer 'npla'.
ctestc
ctest      do 1000 iso=1,2
ctest      isoadd = (iso-1)*ngrid
ctest      isnpad = (iso-1)*kjmx
ctest      nstack = (iso-1)*nmaxp*ngrid
ctest      if(iso.eq.1) then
ctest       nmaxal = nmaxp
ctest      else
ctest       nmaxal = nmaxn
ctest      endif
ctestc
ctestc     square root of inverse mass 'bmass' on 'sqb';
ctestc     d/dr**2 'sqb' / 'sqb' on 'd2b';
ctestc     d/dr('bmass')/r on 'd1b'.
ctestc
ctest      do 10 i=1,ngrid
ctest   10 sqb(i) = sqrt(bmass(i+isoadd))
ctestc
ctest      call derv25(sqb(1),d2b(1),one,ngrid,rstep)
ctest      do 20 i=1,ngrid
ctest   20 d2b(i) = d2b(i)/sqb(i)
ctestc
ctest      call derv15(bmass(isoadd+1),d1b(1),one,ngrid,rstep)
ctest      d1b(1)=(16.0*bmass(2+isoadd)-15.0*bmass(1+isoadd)
ctest     *        -bmass(3+isoadd))/(6.0*rstep**rstep)
ctest      do 30 i=2,ngrid
ctest   30 d1b(i) = d1b(i)*rinvs(i)
ctestc
ctestc     loop over composed ang.mom. 'js'
ctestc      'll' is the orbital ang.mom.;
ctestc      'sig' is the parity of the wave function 'wf' (it is -parity
ctestc      of the state because 'wf' is multiplied by r);.
ctestc      'phasjs' is the degeneracy of the shell.
ctestc      'jsnpl' is pointer to to nucleon shell in pointer field 'npla'.
ctestc
ctest      do 40 js=1,jsmax(iso)
ctest      jsnpl = js+isnpad
ctest      nmax  = npla(1,jsnpl)
ctest      if(nmax.gt.0) then
ctest       ll     = js/2
ctest       sig    = float(mod(ll,2)*2-1)
ctest       xlp    = float(ll*ll+ll)
ctest       xl0    = wls(js)
ctest       phasjs = float(js-ll)
ctestc
ctestc      the effective potential for given 'js' and 'iso' on 'potl';
ctestc      its minimum on 'umin';
ctestc      the zentrifugal force on 'aux3'.
ctestc
ctest       potl(1) = pot(1+isoadd)+xl0*potls(1+isoadd)+d1b(1)
ctest       do 60 i=2,ngrid
ctest       n      = i+isoadd
ctest       y      = xlp*zentrf(i)
ctest       y       = pot(n)+xl0*potls(n)+bmass(n)*y+d1b(i)
ctest       potl(i) = y
ctest   60  continue
ctestc
ctestc      loop over radial quantum number 'n'.
ctestc       accumulating energies and performing one w.f. iteration.
ctestc       'npl' is wave function pointer,
ctestc       'nplshl' is the number of the state in harmonic oscillator ordering.
ctestc
ctest       do 70 n=1,nmax
ctest       npl    = npla(n+1,jsnpl)
ctest       nplm   = npl-1
ctest       nplshl = npl/ngrid+1
ctestc
ctestc      observable on wavefunction
ctestc
ctest       aux1(1) = zero
ctest       do 90 ir=2,ngrid
ctest       r     = (ir-1)*rstep
ctest       aux1(ir) = r*r*wf(nplm+ir)*rinvs(ir)
ctest   90  continue
ctest       call projec(aux1(1),js,nstack,nmaxal)
ctest       call multr(aux1(1),aux1(1),rstep,ngrid)
ctestc
ctestc      inversion to obtain response
ctestc
ctest       irep   = -1
ctest       umin   = spenrg(nplshl)
ctest       call invert(aux1,aux1,irep)
ctestc
ctestc      store on upper corner of 'savlr', wavefunction still in 'wf'-mode.
ctestc
ctest       do 310 ir=1,ngrid
ctest  310  aux1(ir) = rinvs(ir)*aux1(ir)
ctest       call projec(aux1(1),js,nstack,nmaxal)
ctest       do 320 ir=1,ngrid
ctest  320  savlr(nplm+ir,klrpa) = (ir-1)*rstep*aux1(ir)
ctestc
ctest      endif
ctestc
ctest   70  continue
ctest   40  continue
ctest 1000  continue
ctestc
ctestc
ctestc     restore 'd1b'
ctestc
ctest      call derv17(bmass(1),d1b(1),one,ngrid,rstep)
ctest      call derv17(bmass(ngrid+1),d1b(ngrid+1),one,ngrid,rstep)
ctestc
ctest       return
ctest       end
ctestc-----qrobs -------------------------------------------------------------------
ctestc
ctest      subroutine qrobs(resfac)
ctestc
ctest      include 'params.inc'
ctestc
ctestc     repsonse hamiltonian of q-type acting on the mf-repsonse of
ctestc     the observable stored on 'savlr'. result is transferred via
ctestc     'plrout'.
ctestc
ctest      include 'wavfun.inc'
ctest      include 'grid.inc'
ctest      include 'linres.inc'
ctest      include 'result.inc'
ctest      include 'constn.inc'
ctestc
ctestc------------------------------------------------------------------------------
ctestc
ctestc
ctest      include 'initdo.inc'
ctestc
ctestc     action of linear response hamiltonian
ctestc
ctest      sig    = float(mod(la,2)*2-1)
ctest      call derv15(savlr(nadd+1,klrpa),aux1,sig,ngrid,rstep)
ctest      if(resfac.ne.zero)
ctest     &  call hrespq(iso,la,xja,nadd,xj,xl,jdif,
ctest     &              savlr(nadd+1,klrpa),aux1,hlinrs)
ctestc
ctestc     store result and divide 1/r tocome to 'qp'-mode
ctestc
ctest       do 40 ir=1,ngrid
ctest       plrout(nadd+ir,jdif) =  resfac * hlinrs(ir)
ctest   40  continue
ctestc
ctestc     subtract commuted part
ctestc
ctest      if(resfac.ne.zero) then
ctest       call hrespq(iso,la,xja,nadd,xj,xl,jdif,
ctest     &              wf(nadd+1),wfs(nadd+1),hlinrs)
ctest       do 100 nb=1,nmax
ctest       jsb   = jsp(nb)
ctest       wgtb  = weight(nb+nstack/ngrid)
ctest       if(js.eq.jsb .and. wgtb.gt.half) then
ctest        lb    = jsb/2
ctest        xlb   = lb
ctest        xjb   = (jsb-1)/2+half
ctest        naddb = nstack + (nb-1)*ngrid
ctestc
ctest        sum   = zero
ctest        do 110 ir=2,ngrid
ctest        r     = (ir-1)*rstep
ctest  110   sum   = r*hlinrs(ir)*wf(naddb+ir) + sum
ctest        sum   = resfac * rstep*sum
ctestcc        write(14,'(a,2i4,a,g12.3)')
ctestcc     &  ' <wf|v|wf> in qrobs for na,nb=',na,nb,' is:',sum
ctestc
ctest        do 120 ir=2,ngrid
ctest  120   plrout(nadd+ir,jdif) = plrout(nadd+ir,jdif)
ctest     &                       - sum*savlr(naddb+ir,klrpa)*rinvs(ir)
ctestc
ctest       endif
ctestc
ctest  100  continue
ctest      endif
ctestc
ctest      parll  = (1-mod(ll,2)*2)
ctest      call extrap(plrout(nadd+1,jdif),parll)
ctest      call projec(plrout(nadd+1,jdif),js,nstack,nmax)
ctestc
ctest      include 'enddo.inc'
ctestc
ctest      return
ctest      end
ctestc-----probs -------------------------------------------------------------------
ctestc
ctest      subroutine probs(resfac)
ctestc
ctest      include 'params.inc'
ctestc
ctestc     repsonse hamiltonian of p-type acting on the mf-repsonse of
ctestc     the observable stored on 'savlr'. result is transferred via
ctestc     'plrout'.
ctestc
ctest      include 'wavfun.inc'
ctest      include 'grid.inc'
ctest      include 'linres.inc'
ctest      include 'result.inc'
ctest      include 'constn.inc'
ctestc
ctestc------------------------------------------------------------------------------
ctestc
ctest      include 'initdo.inc'
ctestc
ctestc     action of linear response hamiltonian
ctestc
ctest      sig    = float(mod(la,2)*2-1)
ctest      call derv15(savlr(nadd+1,klrpa),aux1,sig,ngrid,rstep)
ctest      if(resfac.ne.zero)
ctest     &   call hrespp(iso,la,xja,nadd,xj,xl,jdif,
ctest     &               savlr(nadd+1,klrpa),aux1)
ctestc
ctestc     store result and divide 1/r to come to 'qp'-mode
ctestc
ctest       do 40 ir=1,ngrid
ctest       plrout(nadd+ir,jdif) =  resfac * hlinrs(ir)
ctest   40  continue
ctestc
ctestc     subtract commuted part
ctestc
ctest      if(resfac.ne.zero) then
ctest        call hrespp(iso,la,xja,nadd,xj,xl,jdif,
ctest     &              wf(nadd+1),wfs(nadd+1))
ctest       do 100 nb=1,nmax
ctest       jsb   = jsp(nb)
ctest       wgtb  = weight(nb+nstack/ngrid)
ctest       if(js.eq.jsb .and. wgtb.gt.half) then
ctest        lb    = jsb/2
ctest        xlb   = lb
ctest        xjb   = (jsb-1)/2+half
ctest        naddb = nstack + (nb-1)*ngrid
ctestc
ctest        sum   = zero
ctest        do 110 ir=2,ngrid
ctest        r     = (ir-1)*rstep
ctest  110   sum   = r*hlinrs(ir)*wf(naddb+ir) + sum
ctest        sum   = resfac * rstep*sum
ctestcc        write(14,'(a,2i4,a,g12.3)')
ctestcc     &  ' <wf|v|wf> in probs for na,nb=',na,nb,' is:',sum
ctestc
ctest        do 120 ir=2,ngrid
ctest  120   plrout(nadd+ir,jdif) = plrout(nadd+ir,jdif)
ctest     &                       - sum*savlr(naddb+ir,klrpa)*rinvs(ir)
ctestc
ctest       endif
ctestc
ctest  100  continue
ctest      endif
ctestc
ctest      parll  = (1-mod(ll,2)*2)
ctest      call extrap(plrout(nadd+1,jdif),parll)
ctest      call projec(plrout(nadd+1,jdif),js,nstack,nmax)
ctestc
ctest       include 'enddo.inc'
ctestc
ctest       return
ctest       end
c-----d2qrho--------------------------------------------------------------------
c
      subroutine d2qrho(varrpr,var0pr,varrsc,var0sc,varrvc,var0vc,
     &                  maxmod,ifgrho)
c
c     computes the second collective derivative of the density.
c     misuses 'wf', 'plrout' and 'd2prho' as workspaces. therefore
c     this routine must be the very last which is called.
c
      include 'params.inc'
c
      character*10 text(2)
c
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      data xmaxlw/1.0e-20/
      data text/' <p|rho|p>','d/dq**2rho'/
      data icont/1/
c
c----------------------------------------------------------------------
c
c
      do 210 ir=1,ngrid
      prho(ir) = zero
      prho(ir+ngrid) = zero
      dprho(ir) = zero
      dprho(ir+ngrid) = zero
      d2prho(ir) = zero
      d2prho(ir+ngrid) = zero
      sigma(ir) = zero
      sigma(ir+ngrid) = zero
      qcur(ir,1) = zero
      qcur(ir+ngrid,1) = zero
      qcur(ir,2) = zero
      qcur(ir+ngrid,2) = zero
  210 continue
      ovpp = zero
      ovpq = zero
      ovqq = zero
      ovpp2 = zero
      ovpq2 = zero
      ovqq2 = zero
c      ovunpr = zero
c      ovpqun = zero
c
c  first part: computation of <p_{ph}|rho|p_{ph}>
c
c   accumulation of [q,q]_ph and [p,p]_ph is done on 'ptau' misused
c   as workspace.
c   trigger is decreasing 'jdif'.
c
      do 377 ir=1,ngrid
      ptau(ir) = zero
      ptau(ir+ngrid) = zero
  377 continue
      jsav = -1
c
      include 'initdo.inc'
c
c     momentum from [t,q]
c
      redmat = redume(xj,alrpa,xja)
      do 80 ir=1,ngrid
      work1(ir) = redmat*qrad(isoadd+ir)*wf(nadd+ir)
      aux1(ir) = work1(ir)*rinvs(ir)
      qlrout(nadd+ir,jdif) = aux1(ir)
   80 continue
c      ovunpr = overlp(qlrout(nadd+1,jdif),qlrout(nadd+1,jdif),
c     &                rstep,ngrid) + ovunpr
      call tpsi(work1(1),work3(1),aux2(1),js,iso,zero)
      call tpsi(wf(nadd+1),work1(1),aux2(1),jsa,iso,zero)
      do 90 ir=2,ngrid
      work3(ir) =
     &         (work3(ir)-redmat*qrad(isoadd+ir)*work1(ir))
     &        *rinvs(ir)
      plrout(nadd+ir,jdif) = work3(ir)
   90 continue
c      ovpqun = overlp(qlrout(nadd+1,jdif),plrout(nadd+1,jdif),
c     &                rstep,ngrid) + ovpqun
      parll  = (1-mod(ll,2)*2)
      call extrap(plrout(nadd+1,jdif),parll)
      call extrap(work3(1),parll)
      call projec(plrout(nadd+1,jdif),js,nstack,nmax)
      call projc2(work3(1),js,nstack,nmax)
      call extrap(aux1(1),parll)
      call projc2(aux1(1),js,nstack,nmax)
      call extrap(qlrout(nadd+1,jdif),parll)
      call projec(qlrout(nadd+1,jdif),js,nstack,nmax)
      do 85 ir=1,ngrid
      work3(ir) = (ir-1)*rstep*work3(ir)
   85 work1(ir) = redmat*qrad(isoadd+ir)*work3(ir)
      call tpsi(work1(1),work2(1),aux2(1),jsa,iso,zero)
      call tpsi(work3(1),work1(1),aux2(1),js,iso,zero)
      do 95 ir=2,ngrid
      work3(ir) =
     &         (work2(ir)-redmat*qrad(isoadd+ir)*work1(ir))
      savlr(nadd+ir,jdif) = work3(ir)
     &        *rinvs(ir)
      work1(ir) = redmat*qrad(isoadd+ir)*aux1(ir)
   95 continue
      parll  = (1-mod(ll,2)*2)
      call extrap(savlr(nadd+1,jdif),parll)
      call projec(savlr(nadd+1,jdif),jsa,nstack,nmax)
      call extrap(work1(1),parll)
      call projec(work1(1),jsa,nstack,nmax)
c
c     accumulate curvatures and reset counter
c
      if(jdif.le.jsav) then
       ovpp2  = overlp(ptau(ngrid+1),ptau(ngrid+1),rstep,ngrid)+ovpp2
       ovpq2  = overlp(ptau(ngrid+1),ptau(1),rstep,ngrid)+ovpq2
       ovqq2  = overlp(ptau(1),ptau(1),rstep,ngrid)+ovqq2
       do 381 ir=1,ngrid
       ptau(ir) = zero
       ptau(ir+ngrid) = zero
  381  continue
      endif
      jsav = jdif
c
      do 75 ir=2,ngrid
      ic    = isoadd+ir
      vol   = wf(nadd+ir)*rinvs(ir)
      qcur(ic,1) = work1(ir)*vol + qcur(ic,1)
      prho(ic) = savlr(nadd+ir,jdif)*vol+prho(ic)
      ptau(ir) = work1(ir) + ptau(ir)
      ptau(ngrid+ir) = savlr(nadd+ir,jdif) + ptau(ngrid+ir)
   75 continue
c
      ovpp  = overlp(plrout(nadd+1,jdif),plrout(nadd+1,jdif),
     &                rstep,ngrid)+ovpp
      ovpq  = overlp(qlrout(nadd+1,jdif),plrout(nadd+1,jdif),
     &                rstep,ngrid)+ovpq
      ovqq  = overlp(qlrout(nadd+1,jdif),qlrout(nadd+1,jdif),
     &                rstep,ngrid)+ovqq
c
      include 'enddo.inc'
c
      call buffql(plrout,'o','pcoll ')
      call buffql(qlrout,'o','qcoll ')
c
      ovpp2  = overlp(ptau(ngrid+1),ptau(ngrid+1),rstep,ngrid)+ovpp2
      ovpq2  = overlp(ptau(ngrid+1),ptau(1),rstep,ngrid)+ovpq2
      ovqq2  = overlp(ptau(1),ptau(1),rstep,ngrid)+ovqq2
c
      call extrap(prho(1),one)
      call extrap(prho(ngrid+1),one)
      call extrap(qcur(1,1),one)
      call extrap(qcur(ngrid+1,1),one)
      widfac = one
      call gscacc(dprho,aux2,plrout,plrout,widfac)
      call extrap(dprho(1),one)
      call extrap(dprho(ngrid+1),one)
      call gscacc(qcur(1,2),aux2,qlrout,qlrout,widfac)
      call extrap(qcur(1,2),one)
      call extrap(qcur(ngrid+1,2),one)
c
c     redundance with q-definition
c
c      write(14,'(a,g12.3)')   ' unprojected: ovqq=',ovunpr
      cosa = ovpq*ovpq/(ovqq*ovpp)
      write(14,'(a,4g12.3)')
     & ' redundance: ovpp,ovpq,ovqq,cos=',ovpp,ovpq,ovqq,cosa
c      write(14,'(a,3g12.3)')
c     & ' ovpq proj - unproj:',ovpq,ovpqun,ovpqun/ovpq
      fac  = (ovpp/ovpq)**2
      cosa = ovpq2*ovpq2/(ovqq2*ovpp2)
      write(14,'(a,4g12.3)')
     & ' red. qq,pp: ovpp,ovpq,ovqq,cos=',ovpp2,ovpq2,ovqq2,cosa
c
      fac4pl = mlrpa*d4pi
      ovqq1 = zero
      ovqp1 = zero
      ovpp1 = zero
      ovqq2 = zero
      ovqp2 = zero
      ovpp2 = zero
      do 450 ir=1,ngrid
      r     = (ir-1)*rstep
      r2    = r*r
      qcur(ir,1) = fac4pl*fac*qcur(ir,1)
      qcur(ir+ngrid,1) = fac4pl*fac*qcur(ir+ngrid,1)
      prho(ir) = fac4pl*prho(ir)
      prho(ir+ngrid) = fac4pl*prho(ir+ngrid)
      qcur(ir,2) = fac*qcur(ir,2)
      qcur(ir+ngrid,2) = fac*qcur(ir+ngrid,2)
      ovqq1 = qcur(ir,1)*qcur(ir,1)*r2 + ovqq1
      ovqp1 = qcur(ir,1)*prho(ir)*r2 + ovqp1
      ovpp1 = prho(ir)*prho(ir)*r2 + ovpp1
      ovqq2 = qcur(ir,2)*qcur(ir,2)*r2 + ovqq2
      ovqp2 = qcur(ir,2)*dprho(ir)*r2 + ovqp2
      ovpp2 = dprho(ir)*dprho(ir)*r2 + ovpp2
  450 continue
      cosa = ovqp2*ovqp2/(ovqq2*ovpp2)
      write(14,'(a,4g12.3)')
     & ' <a|b>: ovpp,ovpq,ovqq,cos=',ovpp2,ovqp2,ovqq2,cosa
      cosa = ovqp1*ovqp1/(ovqq1*ovpp1)
      write(14,'(a,4g12.3)')
     & ' <|ab>: ovpp,ovpq,ovqq,cos=',ovpp1,ovqp1,ovqq1,cosa
c
c  second part: computation of d/dq**2 rho
c
c
c     provide necessary quantum numbers and pointers
c
      do 310 iso=1,2
      if(iso.eq.1) then
       isoadd = 0
      else
       isoadd = ngrid
      endif
c
      jsa   = 1
      js    = 1+lrpa+lrpa
c
c      compute [t,q]rho
c
      do 180 ir=1,ngrid
      work2(ir) = rho(isoadd+ir)*rstep*(ir-1)
  180 work1(ir) = qrad(isoadd+ir)*work2(ir)
      call tpsi(work1(1),work3(1),aux2(1),js,iso,zero)
      call tpsi(work2(1),work1(1),aux2(1),jsa,iso,zero)
      do 190 ir=2,ngrid
      work3(ir) =
     &         (work3(ir)-qrad(isoadd+ir)*work1(ir))
  190 continue
      parll = 2*mod(lrpa,2)-1
      call extrap(work3(1),parll)
c
c      compute [t,q][t,q]rho
c
      do 380 ir=1,ngrid
  380 work1(ir) = qrad(isoadd+ir)*work3(ir)
      call tpsi(work1(1),work2(1),aux2(1),jsa,iso,zero)
      call tpsi(work3(1),work1(1),aux2(1),js,iso,zero)
      do 390 ir=2,ngrid
      work3(ir) =
     &         (work2(ir)-qrad(isoadd+ir)*work1(ir))
     &        *rinvs(ir)*d4pi
      d2prho(isoadd+ir) = work3(ir)
  390 continue
      call extrap(d2prho(isoadd+1),one)
c
  310 continue
c
      if(lrpa.eq.-1) then
c
c     the case lrpa=1 with direct second derivative
c
       do 610 ir=1,ngrid
       r = (ir-1)*rstep
       work1(ir) = r*rho(ir)
       work2(ir) = r*rho(ir+ngrid)
  610  continue
       parity = -one
       call derv25(work1(1),ptau(1),parity,ngrid,rstep)
       call derv25(work2(1),ptau(ngrid+1),parity,ngrid,rstep)
       fac1 = 4.0*amfc(1)**2
       fac2 = 4.0*amfc(2)**2
       do 620 ir=2,ngrid
       ptau(ir) = rinvs(ir)*ptau(ir)*fac1
       ptau(ir+ngrid) = rinvs(ir)*ptau(ir+ngrid)*fac2
  620  continue
c
c       density from p**2
c
       do 110 iso=1,2
       if(iso.eq.1) then
        nmax   = nmaxp
        nstack = 0
        isoadd = 0
       else
       nmax   = nmaxn
        nstack = nmaxp*ngrid
        isoadd = ngrid
       endif
       do 120 na=1,nmax
       jsa   = jsp(na)
       la    = jsa/2
       xla   = la
       xja   = (jsa-1)/2+0.5
       wgt   = weight(na+nstack/ngrid)
       if(wgt.lt.0.5) goto 120
       nadd  = nstack + (na-1)*ngrid
       parity = 2*mod(la,2)-1
c
       call tpsi(wf(nadd+1),work1(1),aux2(1),jsa,iso,zero)
       fac   = sqrt(xja+xja+one)/sqfrpi
       do 130 ir=2,ngrid
       work1(ir) = fac*rinvs(ir)*work1(ir)/amfc(iso)
       work2(ir) = fac*rinvs(ir)*wf(nadd+ir)
  130  continue
       parity  = -(mod(la,2)*2-1)
       call extrap(work1(1),parity)
       call extrap(work2(1),parity)
       call projec(work1(1),jsa,nstack,nmax)
       fac = 4.0*amfc(iso)**2
       do 140 ir=1,ngrid
       ic    = isoadd+ir
       sigma(ic) = fac*work1(ir)*work2(ir)+sigma(ic)
  140  continue
c
  120  continue
  110  continue
c
       call extrap(ptau(1),one)
       call extrap(ptau(ngrid+1),one)
c
      write(14,'(/a/2a)')
     &  ' direct densities:',
     &  '  ir   r    <p|p>    <|p**2>  p**2<>      ',
     &  '  <p|p>    <|p**2>  p**2<>      '
      do 260 ir=1,ngrid
      ic = ngrid+ir
      r = rstep*(ir-1)
      write(14,'(1x,i3,f6.2,2(2x,3g11.3))')
     & ir,r,zero,sigma(ir),ptau(ir),zero,sigma(ic),ptau(ic)
  260 continue
c
       do 410 iso=1,2
       isoadd = (iso-1)*ngrid
       ovtt = overlp(prho(isoadd+1),prho(isoadd+1),rstep,ngrid)
       ovtd = overlp(prho(isoadd+1),sigma(isoadd+1),rstep,ngrid)
       ovdd = overlp(sigma(isoadd+1),sigma(isoadd+1),rstep,ngrid)
       cosa = ovtd*ovtd/(ovtt*ovdd)
       write(14,'(a,4g12.4)')
     &  ' test <|p**2>: <tt>,<td>,<dd>,cos=',ovtt,ovtd,ovdd,cosa
  410  continue
c
       do 420 iso=1,2
       isoadd = (iso-1)*ngrid
       ovtt = overlp(d2prho(isoadd+1),d2prho(isoadd+1),rstep,ngrid)
       ovtd = overlp(d2prho(isoadd+1),ptau(isoadd+1),rstep,ngrid)
       ovdd = overlp(ptau(isoadd+1),ptau(isoadd+1),rstep,ngrid)
       cosa = ovtd*ovtd/(ovtt*ovdd)
       write(14,'(a,4g12.4)')
     &  ' test p**2rho: <tt>,<td>,<dd>,cos=',ovtt,ovtd,ovdd,cosa
  420  continue
      endif
c
c     printings
c
      if(maxmod.eq.1 .and. ifgrho.gt.1) then
       facsum = (one+one)/mlrpa
c
       do 230 iso=1,2
       isoadd = (iso-1)*ngrid
       write(14,'(/a,i3/2a)')
     &   ' densities: iso=',iso,
     &   '  ir   r      <p|p>     <|p**2>    sum         p**2<>      ',
     &   '   <q|q>    <|qq>'
       do 230 ir=1,ngrid
       ic = isoadd+ir
       r = rstep*(ir-1)
       write(14,'(1x,i3,f6.2,2(2x,6g11.3))')
     &  ir,r,dprho(ic),prho(ic),(dprho(ic)+prho(ic))*facsum,
     &  mlrpa*d2prho(ic),qcur(ic,2),qcur(ic,1)
  230  continue
c
c     fourier-bessel transform
c
       qsurf  = 5.600/rdmsc
       if(iso.eq.1) then
        fac   = one/aprot
       else
        fac   = one/aneut
       endif
       do 530 iso=1,2
       isoadd = (iso-1)*ngrid
       write(14,'(/a,i5,a,f8.4/2a)')
     &   ' formfactors: iso=',iso,'  qsurf=',qsurf,
     &   '  iq   q      <p|p>     <|p**2>    sum         p**2<>      ',
     &   '  <q|q>    <|qq>'
c
       delqh = 0.5*delq
       do 530 iq=1,ngrid
       q1     = iq*delqh-delqh
        sppp   = fac*fbform(q1,dprho(1+isoadd))
        sp2p   = fac*fbform(q1,prho(1+isoadd))
        sd2p   = fac*fbform(q1,d2prho(1+isoadd))
        sqqp   = fac*fbform(q1,qcur(1+isoadd,2))
        sq2p   = fac*fbform(q1,qcur(1+isoadd,1))
       write(14,'(1x,i3,f6.2,2(2x,6g11.3))')
     &  iq,q1,sppp,sp2p,facsum*(sppp+sp2p),sd2p,sqqp,sq2p
  530  continue
c
      endif
c
c     r.m.s. corrections
c
      sppp2  = zero
      sppn2  = zero
      sp2p2  = zero
      sp2n2  = zero
      sd2p2  = zero
      sd2n2  = zero
      sqqp2  = zero
      sqqn2  = zero
      sq2p2  = zero
      sq2n2  = zero
      r      = zero
       do 500 ir=2,ngrid
       r     = rstep+r
       r4    = r*r
       r4    = r4*r4
       sppp2  = r4*dprho(ir)+sppp2
       sppn2  = r4*dprho(ir+ngrid)+sppn2
       sp2p2  = r4*prho(ir)+sp2p2
       sp2n2  = r4*prho(ir+ngrid)+sp2n2
       sd2p2  = r4*d2prho(ir)+sd2p2
       sd2n2  = r4*d2prho(ir+ngrid)+sd2n2
       sqqp2  = r4*qcur(ir,2)+sqqp2
       sqqn2  = r4*qcur(ir+ngrid,2)+sqqn2
       sq2p2  = r4*qcur(ir,1)+sq2p2
       sq2n2  = r4*qcur(ir+ngrid,1)+sq2n2
  500  continue
      facp = frpirs/aprot
      facn = frpirs/aneut
      sppp2  = facp*sppp2
      sppn2  = facn*sppn2
      sp2p2  = facp*sp2p2
      sp2n2  = facn*sp2n2
      sqqp2  = facp*sqqp2
      sqqn2  = facn*sqqn2
      write(14,'(/a/2a)')
     &  ' r.m.s. corrections: ',
     &  '    <p|p>     <|p**2>    sum         p**2<>      ',
     &  '  <q|q>    <|qq>'
      write(14,'(2x,6g11.3/2x,6g11.3)')
     &     sppp2,sp2p2,facsum*(sppp2+sp2p2),
     &     facp*sd2p2,sqqp2,facp*sq2p2,
     &     sppn2,sp2n2,facsum*(sppn2+sp2n2),
     &     facn*sd2n2,sqqn2,facn*sq2n2
c
c     combine to new 'recipe'
c
      qfac   = (ovpq/ovpp)**2/ovqq
      pfac   = ovqq/ovpq**2
      rpoldp = 0.125*(varrpr-var0pr)**2/(varrpr*var0pr)/(3.0)
      rpfacp = 0.125*((varrsc-var0sc)**2/(varrsc*var0sc)
     &                +(varrvc-var0vc)**2/(varrvc*var0vc))/(3.0)
      rpfacs = 0.125*(varrsc-var0sc)**2/(varrsc*var0sc)/(3.0*anumb)
      write(14,'(/a/9x,a)')
     &   ' new recipe from full <qoq> and <pop>:',
     &   '   <qoq>       <pop>       mixed'
      qpart  = rpfacp*qfac*sqqp2
      ppart  = rpfacp*pfac*sppp2
      write(14,'(a,3g12.3)')
     &   ' proton:',2.0*qpart,2.0*ppart,qpart+ppart
      qpart  = rpfacs*qfac*(aprot*sqqp2+aneut*sqqn2)
      ppart  = rpfacs*pfac*(aprot*sppp2+aneut*sppn2)
      write(14,'(a,3g12.3)')
     &   ' total:',2.0*qpart,2.0*ppart,qpart+ppart
      simfcp = (qfac*sqqp2/pfac+sppp2)/(sppp2+sp2p2)
      simfcs = (qfac*(aprot*sqqp2+aneut*sqqn2)/pfac
     &          +aprot*sppp2+aneut*sppn2)/
     &         (aprot*(sppp2+sp2p2)+aneut*(sppn2+sp2n2))
      write(14,'(2(a,g12.3))')
     &   ' factor new/old: proton=',simfcp,'  total=',simfcs,
     &   ' factor sum{salar+vector}/proton:',rpfacp/rpoldp
c
c     new densities and formfactors
c
      rpfacp = 3.0*rpfacp
      rpfacs = 3.0*anumb*rpfacs
      if(ifgrho.gt.1) write(14,'(3(/a))')
     &  ' densities:',
     &  '              new recipe              simple recipe',
     &  ' ir   r     proton     total       proton    total'
      do 510 ir=1,ngrid
      fulrhp = rpfacp*(qfac*qcur(ir,2)+pfac*dprho(ir))
      fulrhs = rpfacs*(qfac*(qcur(ir,2)+qcur(ir+ngrid,2))
     &                +pfac*(dprho(ir)+dprho(ir+ngrid)))
      simrhp = rpfacp*pfac*(dprho(ir)+prho(ir))
      simrhs = rpfacs*pfac*(dprho(ir)+prho(ir)+
     &                      dprho(ir+ngrid)+prho(ir+ngrid))
      if(ifgrho.gt.1) write(14,'(i5,f6.2,4g12.3)')
     &  ir,(ir-1)*rstep,fulrhp,fulrhs,simrhp,simrhs
      dprho(ir) = fulrhp
      dprho(ir+ngrid) = fulrhs
      prho(ir) = simrhp
      prho(ir+ngrid) = simrhs
  510 continue
c
      if(ifgrho.gt.1) then
       write(14,'(3(/a))')
     &  ' formfactors:',
     &  '              new recipe              simple recipe',
     &  ' iq   q     proton     total       proton    total'
       delqh = 0.5*delq
       do 515 iq=1,ngrid
       q1     = iq*delqh-delqh
       write(14,'(i5,f7.3,4g12.3)')
     &   iq,q1,fbform(q1,dprho(1))/aprot,
     &         fbform(q1,dprho(1+ngrid))/anumb,
     &         fbform(q1,prho(1))/aprot,
     &         fbform(q1,prho(1+ngrid))/anumb
  515  continue
      endif
c
      call formcr(dprho(1),aprot,rdmsp,surfp,delrp,delsp)
      call formcr(dprho(1+ngrid),anumb,rdms,surf,delrt,delst)
      call formcr(prho(1),aprot,rdmsp,surfp,delr0p,dels0p)
      call formcr(prho(1+ngrid),anumb,rdms,surf,delr0t,dels0t)
      write(14,'(2(/a),2(/a,4g12.3))')
     &  '                  new recipe              simple recipe',
     &  '                 proton     total       proton    total',
     &  ' delta_sigma:',delsp,delst,dels0p,dels0t,
     &  ' diff.radius:',delrp,delrt,delr0p,delr0t
c
      return
      end
c-----projc2------------------------------------------------------
c
      subroutine projc2(awf,js,nstack,nmax)
c
      include 'params.inc'
c
c     projects out all wave function 'wf' (given in 'wavfun.inc')
c     from the wave function 'awf'.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'constn.inc'
c
      dimension awf(kg1)
c
c---------------------------------------------------------------------
c
c     check bounds
c
      if(nmax.gt.max(nmaxp,nmaxn)) then
       write(6,'(a,i8,a,2i8)')
     &   ' in "projec" was nmax=',nmax,' > nmaxp,nmaxn=',nmaxp,nmaxn
       stop ' "nmax" too large in "projec" '
      endif
c
c     loop over the a-states
c
      do 10 na=1,nmax
      jsa   = jsp(na)
      wgt   = weight(na+nstack/ngrid)
      if(js.eq.jsa .and. wgt.gt.half) then
       la    = jsa/2
       xl    = la
       xja   = (jsa-1)/2+half
       nadd  = nstack + (na-1)*ngrid
c
c      overlap
c
       sum   = zero
       do 110 ir=2,ngrid
  110  sum   = (ir-1)*rstep*awf(ir)*wf(nadd+ir) + sum
       sum   = rstep*sum
c
c      project 'wf' twice such that [a,a]_{ph} is produced
c
       sum   = 2.0*sum
       do 120 ir=2,ngrid
  120  awf(ir) = awf(ir) - sum*wf(nadd+ir)*rinvs(ir)
c
      endif
   10 continue
c
      ll     = js/2
      if(ll.eq.0) then
       call extrap(awf,one)
      else
       awf(1) = zero
      endif
c
      return
      end
c-----fbform------------------------------------------------------
c
      function fbform(q1,rhoin)
c
      include 'params.inc'
c
c     fourier-bessel transformation of 'rhoin' at momentum 'q1'.
c     resulting formfactor on 'fbform'.
c
      include 'grid.inc'
      include 'constn.inc'
c
      dimension rhoin(kg1)
c
c---------------------------------------------------------------------
c
      r      = zero
      sum    = zero
      if(q1.ne.zero) then
       dc     = cos(q1*rstep)
       ds     = sin(q1*rstep)
       c      = one
       s      = zero
       do 510 ir=2,ngrid
       r      = rstep+r
       x      = c*dc-s*ds
       s      = s*dc+c*ds
       c      = x
       vol    = r*s
       sum    = vol*rhoin(ir)+sum
  510  continue
       fbform = frpirs*sum/q1
      else
       do 520 ir=2,ngrid
       r       = rstep+r
       vol     = r*r
       sum    = vol*rhoin(ir)+sum
  520  continue
       fbform = frpirs*sum
      endif
c
       return
       end
c-----formcr-----------------------------------------------------------
c
      subroutine formcr(rhoin,partn,rdmsin,surfin,delrad,delsur)
c
      include 'params.inc'
c
c     correction to diffraction-radius and surface from
c     correlation density on 'rhoin' with particle number 'partn',
c     radius 'rdmsin', and surface 'surfin'.
c
      include 'grid.inc'
      include 'constn.inc'
      include 'result.inc'
      include 'wavfun.inc'
c
      dimension rhoin(kg1)
c
c---------------------------------------------------------------------
c
c    surface correction
c
      qsurf  = 5.600/rdmsin
      q1     = qsurf
      if(partn.ne.zero) dform  = fbform(q1,rhoin)/partn
      qr     = q1*rdmsin
      formbx = ((sin(qr)/qr-cos(qr))*3.0/(qr*qr))
      formf0 = exp(-0.5*(surfin*q1)**2)*formbx
      delsur = -dform/(surfin*q1*q1*formf0)
c
c    radius correction
c
      qrdms   = 4.493/rdmsin
      q1     = qrdms
      r      = zero
      sump   = zero
      dc     = cos(q1*rstep)
      ds     = sin(q1*rstep)
      c      = one
      s      = zero
      do 510 ir=2,ngrid
      r      = rstep+r
      x      = c*dc-s*ds
      s      = s*dc+c*ds
      c      = x
      vol    = r*r*(c-s/(q1*r))
      sump   = vol*rho(ir)+sump
      sumn   = vol*rho(ir)+sumn
  510 continue
      if(partn.eq.aprot) then
       forms = frpirs*sump/q1
      elseif(partn.eq.aneut) then
       forms = frpirs*sumn/q1
      elseif(partn.eq.(aprot+aneut)) then
       forms = frpirs*(sump+sumn)/q1
      else
       write(6,'(a,g12.3)') ' invalid particle number:',partn
       forms = 1.0e20
      endif
c
      delrad = rdmsin*fbform(q1,rhoin)/(q1*forms)
c
      return
      end
c-----d2qenr--------------------------------------------------------------------
c
      subroutine d2qenr(varrpr,var0pr,varrsc,var0sc,varrvc,var0vc,
     &                  maxmod)
c
c     computes the second collective derivative of the energy
c     and according collectively filtered gsc.
c
      include 'params.inc'
c
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      data xmaxlw/1.0e-20/
      data icont/1/
c
c----------------------------------------------------------------------
c
c     first part: isoscalar contribution, uses 'qlrout' and 'plrout'
c                 from previous call to 'd2qrho'.
c
      call buffql(qlrout,'i','qcoll ')
      ovqq  = qpover(qlrout,qlrout)
      call commhq(ovqhq,ovqh0q,ovhqhq,ovh0qq,1,one,1)
      call buffql(plrout,'i','pcoll ')
      ovpq  = qpover(qlrout,plrout)
      call commhp(ovpqb,ovpqb0,ovqbqb,ovqbb0,one,zero,1)
      rpfacs = 0.125*(varrsc-var0sc)**2/(varrsc*var0sc)/(3.0)
      corens = rpfacs*(ovqhq/ovqq+ovpqb*ovqq/ovpq**2)
c
c     second part: isovector contribution, modifies previous 'qlrout'
c                  and 'plrout'.
c
      call buffql(qlrout,'i','qcoll ')
      call buffql(plrout,'i','pcoll ')
      facpr  = aneut/anumb
      facne  = -aprot/anumb
      include 'initdo.inc'
      if(iso.eq.1) then
       fac    = facpr
      else
       fac    = facne
      endif
      do 110 ir=1,ngrid
      qlrout(nadd+ir,jdif) = fac*qlrout(nadd+ir,jdif)
      plrout(nadd+ir,jdif) = fac*plrout(nadd+ir,jdif)
  110 continue
      include 'enddo.inc'
      call buffql(qlrout,'o','qcoll ')
      call buffql(plrout,'o','pcoll ')
c
      ovqq  = qpover(qlrout,qlrout)
      call commhq(ovqhq,ovqh0q,ovhqhq,ovh0qq,1,one,1)
      call buffql(plrout,'i','pcoll ')
      ovpq  = qpover(qlrout,plrout)
      call commhp(ovpqb,ovpqb0,ovqbqb,ovqbb0,one,zero,1)
      rpfacv = 0.125*(varrvc-var0vc)**2/(varrvc*var0vc)/(3.0)
      corenv = rpfacv*(ovqhq/ovqq+ovpqp*ovqq/ovpq**2)
c
c     printings
c
      write(14,'(/3(a,g12.3))')
     &  ' correlation energy: scalar=',corens,'  vector=',corenv,
     &  '  total=',corens+corenv
      write(14,'(a,4g12.3)')
     &  '   ovqhq,ovqq,ovphp,ovpq:',ovqhq,ovqq,ovpqb,ovpq
c
      return
      end
