c-----sumrul-----------------------------------------------------
c
      subroutine sumrul
c
      include 'params.inc'
      parameter(krpmat=krpmax*(krpmax+1)/2)   ! auxiliary for RPA matrices
      parameter(kreczw=kgr2*kvarln,           ! auxiliary for buffer files
     &          krec17=5*kreczw)
      logical tnols,tnomas,tnopot,tnocoul     ! switch for tests
      parameter (tnols=.false.,tnomas=.false.,tnopot=.false.)
      parameter (tnocoul=.false.)
c      ---> see also test switches at the bottom of 'linres.inc'.
c
c     rpa & fluid dynamics routine for states with natural parity:
c     is called from 'hafosu' if 'ifgsc'=1.
c
c     Note: The normalization of the multipole moments is not always
c           correct. For example, for the L=21 channel multiply all
c           strengths by a factor 4.19 (= 4*Pi/3).
c
c     the options are read-in at the beginning of this routine.
c     they are:
c      lrpa    = angular momentum of the mode
c      ifprim  = computation of ph-matrix elements for one mode only
c      ifsepi  = switch to separable interaction
c                0 --> full residual interaction
c                1 --> single separable term Q*Q
c                2 --> Skyrme-SepRPA in VPM style
c                      (take selection of modes from standard
c                       RPA input for modes, only local modes allowed)
c      cutmod  = cut-off for redundant modes. eigenvalues of mass matrix
c                below cutmod are discarded.
c      resfac  = factor on residual interaction (standard 1.0).
c                < 0 switches to damping of the residual interaction
c                with abs(resfac).
c      cmfac   = strength of r_cm*r_cm residual interaction which
c                serves allow inversion of the hamiltonian in case
c                of the isoscalar dipole mode (cm motion).
c      iforth  = level of pre-orthogonalization of the states
c                -1 --> no pre-orthogonalization
c                0 --> pre-orthogonalization but no second run
c                1 --> pre-orthogonalization with one extra run using
c                      resfac=0.0
c                2 --> pre-orthogonalization with one extra run using
c                      full resfac
c      ifirpa  = switch for inverse rpa-moments
c      itrpa   = number of iterations in inverse h-action
c      epsrpa  = stepsize in inverse h-action
c      dabrpa  = required precision in inverse h-action
c      ifgrho  = switch for computation of ground-state correlations
c                1 --> collective recipe
c                2 --> full corr. energy (incl. collective recipe)
c                3 --> full corr. density (incl. collective recipe)
c                4 --> transition densities are collected and printed.
c      ifppsc  = switch for pp- or hh- part in the scattering diagram:
c                0 --> all parts (standard)
c                1 --> only pp-parts
c                2 --> only hh-parts
c      ifanlr  = level of analysis; lowest is -1, final modes are
c                composed explicitely and widths are computes for 0,
c                even higher values give more detailed printout.
c
c     the selection of modes is read in this routine line by line and
c     the information on the modes is evaluated in 'qinit'.
c     the mode parameters are
c      typmod  = local form
c                 'mom=m' --> initialize r**(lrpa+m)
c                 'mot=m' --> initialize r**(lrpa+m)*Woods_Saxon
c                 'bes=q' --> initialize j_lrpa(q*r/r)
c                 'ph=(p,h)' --> explicit particle-hole state (p,h)
c                 'end'   --> terminates input of modes
c                the number of modes is determined dynamically by reading
c                until the 'end' option
c      isorpa  = isospin of the mode (1,2)
c      jselec  = option for picking selective (j,ja,iso)-transitions from
c                a local operator
c                '(js,jsa,iso)' --> only transition jsa -> js is taken
c                'all+n'    --> all transitions separately up to js=jsa+n
c                'nojj'     --> no selective picking
c      isplit  = splitting of upper and lower momentum component:
c                 0 --> no splitting = strictly local operator
c                 1 --> select lower component of corresponding momentum
c                 2 --> select upper component of corresponding momentum
c      ihwgt   = hamiltonian weight h**ihwgt on q
c      ih0wgt  = mean-field hamiltonian weight h0**ihwgt on q
c      ipwgt   = hamiltonian weights on p, p = h**ipwgt q (use 0,1,-1);
c                standard is ipwgt=1, for 'ph' use ipwgt=0, low modes
c                may prefer ipwgt=-1
c      iprshl  = switches projection on-off valence shell
c                0 --> no projection
c                1 --> projects q on valence shell
c                2 --> projects q off valence shell
c                the valence shell is defined by the unoccupied levels
c                in the selection provided by nmaxp and nmaxn.
c
!      resulting B(EL) values are normalized to the following sum rules:
!      isovector L=1: <[rY_10,[H.rY_10]]> = (3NZ/4piA)*hbar**2/2m
!      isoscalar L=2: <r^2Y_20,[H,r^2Y_20]]> = (5Ar^2/2pi)*hbar**/m
!                     but should by (5Ar^2/3pi)*hbar**2/m
!                     ==> correction factor 3/2
c
      include 'skypar.inc'
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'numpar.inc'
      include 'gscres.inc'
      include 'constn.inc'
      include 'linres.inc'
c
      dimension brpa(krpmat),crpa(krpmat),brpa0(krpmat),crpa0(krpmat),
     &          qqhf(krpmat),qqrpa(krpmat),qqrpa0(krpmat),
     &          pphf(krpmat),pprpa(krpmat),pprpa0(krpmat),
     &          vect(krpmax,krpmax),vec2(krpmax,krpmax),
     &          alpha(krpmax,krpmax),
     &          qpinv(krpmax,krpmax),
     &          broot(krpmax),croot(krpmax),
     &          broot0(krpmax),croot0(krpmax),
     &          qtrapr(krpmax),qtrane(krpmax)
      double precision kappa(krpmax,krpmax),eta(krpmax,krpmax)
      REAL(8) :: workrpa(krpmax)
      dimension dgsrho(kgr2),dhhrho(kgr2),qrpwid(krpmax)
      dimension nparsv(krpmax),nholsv(krpmax),isophs(krpmax)
      equivalence (qqhf(1),brpa(1)),(qqrpa(1),crpa(1)),
     &            (qqrpa0(1),crpa0(1)),(vect(1,3),isophs(1)),
     &            (vect(1,1),nparsv(1)),(vect(1,2),nholsv(1))
c
      logical secrun
c
      character*2 ichr(kichr)
      character*16 typmod
      character*12 jselec
      character*50 toper
      character*192 atext
      character*(7) qop,pop,qbop,qrop,prop,xop,yop,callop
c      data ichr/'1','2','3','4','5','6','7','8','9','0',
c     &          'a','b','c','d','e','f','g','h','i','j',
c     &          'k','l','m','n','o','p','q','r','s','t',
c     &          'u','v','w','x','y','z',
c     &          'a','b','c','d','e','f','g','h','i','j',
c     &          'k','l','m','n','o','p','q','r','s','t',
c     &          'u','v','w','x','y','z',
c     &          '!','@','#','$','%','^','&','*','(',')','-','_'/
      data     qop,      pop,    qbop,    qrop,    prop
     &    /'qin   ','p     ','[h,p ]','qres  ','pres  '/
      data     xop,      yop
     &    /'xsep  ','ysep   '/
      data qqhf/krpmat*0.0D0/
      data ifprim,cutmod,resfac,cmfac,ecut,ifppsc
     &    /     1,1.0D-5,   1.0D0,  0.0D0,1.0D+20,0/
c
      data ifop14,ifop17/0,0/
      LOGICAL,PARAMETER :: testrot=.FALSE.
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
c
c      IF(ioptsk==13) THEN
c        WRITE(6,*) ' RPA not yet applicable with IOPTSK=13'
c        WRITE(14,*) ' RPA not yet applicable with IOPTSK=13'
c        RETURN
c      END IF
c      if(iftm.ne.0) stop ' RPA only for IFTM=0'
      ifsepi = 0
      resdmp = 0.0d0      
      tspinorbit = t4.ne.zero .OR. (ioptsk==13 .AND. clsp.NE.0D0)
c
c      set pieces of force to zero for testing parts of residual interaction
c
       if(tnopot) then
         t0 = 0.0D0
         x0 = 0.0D0
         t3 = 0.0D0
         x3 = 0.0D0
         t7  = 0.0D0
         t8  = 0.0D0
      endif
      if(tnomas) then
         t1 = 0.0D0
         x1 = 0.0D0
         t2 = 0.0D0
         x2 = 0.0D0
         t5  = 0.0D0
         t6  = 0.0D0
      endif
      if(tnols) then
         t4 = 0.0D0
         b4 = 0.0D0
         b4p =0.0D0
      endif
      if(tnocoul) then
         acoul = 0.0D0
      endif
       write(14,*) ' t0,x0,t3,x3=',t0,x0,t3,x3
       write(14,*) ' t4,b4,b4p=',t4,b4,b4p
c
c      check bounds
c
       do i=1,kichr
         ichr(i) = char(mod(i,256))//char(i/256)
c         write(6,'(1x,a,i5,a2)') ' i,key=',i,ichr(i)
       enddo
       if(krpmax.gt.kichr) stop ' krpmax too large'
c
c      read parameters for 'sumrul'
c
       iforth  = 0
       qobserv = 0.0d0
       read(5,*) lrpa,resfac,cmfac,cutmod,ifsepi,resdmp,iforth,ecut
       read(5,*) itrpa,epsrpa,dabrpa,ifanlr,ifirpa,ifgrho,ifppsc,qobserv
       if(2*lrpa+1.gt.klrpa) 
     &  stop ' code was not compiled for that large LRPA '
c
c     open file 14 for concentrated output of matrix elements
c
      if(ifgrho.eq.4) then
        open(unit=17,status='scratch',form='unformatted',
     &       access='direct',recl=5*kvarln*ngrid2)
      elseif(ifgrho.eq.2 .or. ifgrho.eq.3) then
        open(unit=17,status='scratch',form='unformatted',
     &       access='direct',recl=kvarln*krpmat)
      endif
      if(ifop14.eq.0) then
       CLOSE(14)
       open(unit=14,status='REPLACE',form='formatted',file='for014')
       ifop14 = 1
      endif
c
c     write a header
c
      write(14,'(a/1x,2a/1x,2(a,i4),5x,2(a,i4))')
     &  '1',
     &  '  title:',title,
     &  'nucleus:  nprot=',nprot,'  nneut=',nneut,
     &  ' shells:  nmaxp=',nmaxp,'  namxn=',nmaxn
      write(14,'(a,f5.2,a,i5,a,f5.2)')
     &  '   grid:  rstep:',rstep,'  ngrid=',ngrid,
     &  '   wall=',wall
      write(14,'(1x,a,5g13.6/7x,a,5g13.6/7x,a,g13.6)')
     &  'force: t0-t4=',t0,t1,t2,t3,t4,
     &        ' x0-x3=',x0,x1,x2,x3,b4p,
     &        ' power=',power
      write(14,'(7x,3(a,g13.6)/)')
     &        ' tyuk=',tyuk,' xyuk=',xyuk,' ayuk=',ayuk
      if(iforth.ne.0) write(14,'(/a)')
     &       ' modes are preorthonormalized'
      write(14,'(9x,a,2(/a,3g12.5,5x,g12.5))')
     &  '   r_rms      r_diffr     surface          energy',
     &  ' mass:   ', rms, rdms, surf,etotal,
     &  ' charge: ',rmsc,rdmsc,surfc
      spemxp = -9999.0
      do n=1,nmaxp
        if(weight(n).gt.half) spemxp=max(spenrg(n),spemxp)
      enddo
      spemxn = -9999.0
      do n=nmaxp+1,nmaxp+nmaxn
        if(weight(n).gt.half) spemxn=max(spenrg(n),spemxn)
      enddo
      write(14,'(a,2(1pg12.4))')
     &  ' continuum threshold protons/neutron:',spemxp,spemxn
c
cold      if(KSCDIP.eq.1 .and. lrpa.eq.1) call iniobs()
      if(KSCDIP.eq.1) call iniobs()
c
c      sumrul-part
c
c      first supply derivative of 'bmass' on 'd1b' for use in 'hpsi'
c
      call derv17(bmass(1),d1b(1),one,ngrid,rstep)
      call derv17(bmass(ngrid+1),d1b(ngrid+1),one,ngrid,rstep)
c
c     provide derivative of "rho"
c
      call derv17(rho(1),d1r(1), one,ngrid,rstep)
      call derv17(rho(ngrid+1),d1r(ngrid+1), one,ngrid,rstep)
c
c     write options and parameters for rpa part
c
      write(6,'(/1x,a,i3,3(a,g10.3)/10x,a,i5,a,g10.3)')
     &   ' rpa with: lrpa=',lrpa,'  res.int.fac=',resfac,
     &   '  cminv.fac=',cmfac,'  cutmode=',cutmod,
     &   ' ifsepi=',ifsepi,'   resdmp=',resdmp
      write(6,'(21x,a,i5,a,f7.3,a,g12.4)')
     &    '  itrpa=',itrpa,' epsrpa=',epsrpa,' daprpa=',dabrpa
      write(14,'(/1x,a,i3,3(a,g10.3)/10x,a,i5,a,g10.3,a,i5)')
     &   ' rpa with: lrpa=',lrpa,'  res.int.fac=',resfac,
     &   '  cminv.fac=',cmfac,'  cutmode=',cutmod,
     &   ' ifsepi=',ifsepi,'   resdmp=',resdmp,' ifgrho=',ifgrho
      if(KSCDIP.eq.1) write(14,'(a)') ' now running r**3 for dipole '

c
c      set parameters inherited from 'extrsp'
c
       ifcmpl = 0
       enext  = zero
       enwid  = zero
c
c     set the other angular momentum parameters
c
      alrpa=lrpa
      mlrpa=2*lrpa+1
      call setred
c      write(14,'(a,i4)') ' IFGRHO=',ifgrho
c
c     initialise buffering for response wave functions
c
      call buffql(qlrout,'c','      ')
c
c     here follows the main part, the computation of q-p-rpa.
c
c     compute response operators and commutators, store on scratch.
c     the type of the modes are read line by line. the loop over
c     the modes is programmed explicitely because the loop index
c     'imode' is changed dynamically depending on the initilization
c     in 'qinit'.
c     'iflagr' sets a flag for reading input from file 05;
c     'iflagr>0' indicates that a given string of initializations
c     is yet to be worked up.
c
      iflagr = 0
c
c     save some options and set counters for pre-orthonormalization.
c
      resfcs = resfac
      if(iforth.gt.0) then
       if(iforth.eq.1)  resfac = 0.001
       secrun = .false.
      else
       secrun = .true.
      endif
!
!     quick test of "div"  ==> "grad" and "div" have wrong sign
!
!      DO ir=1,ngrid
!        qrad(ir) = ((ir-1)*rstep)**2
!      END DO
!      CALL grad(qrad,rpotls,1)
!      CALL div(rpotls,qrad,1)
!      WRITE(14,'(a)') '# test grad and div'
!      DO ir=2,ngrid
!        ri = 1D0/((ir-1)*rstep)
!        WRITE(14,'(i5,f8.2,3(1pg13.5))') 
!     &    ir,(ir-1)*rstep,qrad(ir),rpotls(ir,1)*ri,rpotls(ir,2)*ri
!      END DO
!      STOP
!
!     quick test of "rot"
!
      IF(testrot) THEN
      DO ir=1,ngrid
        work1(ir) = ((ir-1)*rstep)*rho(ir)
      END DO
      CALL grad(work1,rpotls,1)
      WRITE(14,'(a,3(1pg13.5))') 'after GRAD:',
     &  SUM(rpotls(1:ngrid,1)*rsquar(1:ngrid)),
     &  SUM(rpotls(1:ngrid,2)*rsquar(1:ngrid)),
     &  SUM(rho(1:ngrid)*rsquar(1:ngrid))
      DO ir=1,ngrid
!        rpotls(ir,2) = ((ir-1)*rstep)*rho(ir)
        rpotls(ir,1) = 0D0
      END DO
      WRITE(14,'(a,2(1pg13.5))') 'before ROT:',
     &  SUM(rpotls(1:ngrid,1)*rsquar(1:ngrid)),
     &  SUM(rpotls(1:ngrid,2)*rsquar(1:ngrid))
      CALL rot(rpotls,qrad,1)
      WRITE(14,'(a)') '# test rot upper'
      DO ir=2,ngrid
        WRITE(14,'(i5,f8.2,3(1pg13.5))') 
     &    ir,(ir-1)*rstep,qrad(ir)
      END DO
      DO ir=1,ngrid
        rpotls(ir,1) = ((ir-1)*rstep)*rho(ir)*SQRT(alrpa/MLRPA)
        rpotls(ir,2) = -((ir-1)*rstep)*rho(ir)*SQRT((alrpa+1D0)/MLRPA)
!        rpotls(ir,1) = ((ir-1)*rstep)*SQRT(alrpa/MLRPA)
!        rpotls(ir,2) = -((ir-1)*rstep)*SQRT((alrpa+1D0)/MLRPA)
      END DO
!      rpotls=0D0
      CALL grad(work1,rpotls,1)
      WRITE(14,'(a,2(1pg13.5))') 'after GRAD:',
     &  SUM(rpotls(1:ngrid,1)*rsquar(1:ngrid)),
     &  SUM(rpotls(1:ngrid,2)*rsquar(1:ngrid))
      CALL rot(rpotls,qrad,1)
      WRITE(14,'(a)') '# test rot full'
      DO ir=2,ngrid
        WRITE(14,'(i5,f8.2,3(1pg13.5))') 
     &    ir,(ir-1)*rstep,qrad(ir)
      END DO
      CALL grad(work1,rpotls,1)
      DO ir=1,ngrid
!        rpotls(ir,1) = (((ir-1)*rstep))*rho(ir)
        rpotls(ir,2) = 0D0
      END DO
      CALL rot(rpotls,qrad,1)
      WRITE(14,'(a)') '# test rot lower'
      DO ir=2,ngrid
        WRITE(14,'(i5,f8.2,3(1pg13.5))') 
     &    ir,(ir-1)*rstep,qrad(ir)
      END DO
      STOP
      END IF
c
c     the gsc runs twice, the second time with resfac = 0.1*resfcs.
c     this uses pseudo-loop on 340.
c
      irunc = 0
  340 continue
      irunc = 1 + irunc
c      write(14,'(a,2i4)') ' IFGRHO,IRUNC=',ifgrho,irunc
c
c     the pre-orthonormalization uses pseudo-loop on 330.
c
  330 continue
      imode  = 0
   30 continue
      if(.not.secrun .or. iforth.le.0) then
       if(iflagr.eq.0)  then
        write(6,'(a)') ' now reading typmod'
        read(5,*) typmod,isorpa,jselec,isplit,ihwgt,ih0wgt,ipwgt,iprshl
        write(6,'(1x,a,i5,a,5i5)')
     &    typmod,isorpa,jselec,isplit,ihwgt,ih0wgt,ipwgt,iprshl
       endif
c                                                     setup the Q operator
       call qinit(typmod,jselec,isplit,ihwgt,ih0wgt,ipwgt,
     &            iprshl,imode,ifprim,iflagr,xkqrd,resfac,cmfac,
     &            nparsv,nholsv,isophs,toper)
c       write(14,'(a,g12.4)') ' resfac=',resfac
       if(iflagr.le.-999) then
        maxmod = imode
        goto 39
       endif
      else
c                                            pick operators from 'qres'
       imode = 1+imode
       if(imode.gt.maxmod) goto 39
       write(14,'(a,i3)') ' recoupled mode nr.',imode
c       write(14,'(a,g12.4)') ' resfac=',resfac
       qrop(5:6) = ichr(imode)
       prop(5:6) = ichr(imode)
       call buffql(qlrout,'i',qrop)
       call buffql(plrout,'i',prop)
       pop(2:3 ) = ichr(imode)
       call buffql(plrout,'o',pop)
       ipwgt = -1000
      endif
c                                            store preliminary q-operator
      qop(4:5 ) = ichr(imode)
      call buffql(qlrout,'o',qop)
c                                            orthogonalise on previous q
c      write(14,'(a,i4)') ' after IORTH: IFGRHO=',ifgrho
      if(imode.gt.1 .and. iforth.ge.0) then
       if(krpmax.le.4) stop ' krpmax<5 !!! '
       do imod2=1,imode-1
         pop(2:3 ) = ichr(imod2)
         call buffql(plrout,'i',pop)
         ov     = qpover(qlrout,plrout)
         ipoint = (imod2*(imod2+1))/2
         vect(imod2,krpmax) = ov/brpa(ipoint)
         if(ifsepi.ne.0) 
     &       write(14,*) ' imod2,ov,brpa=',imod2,ov,brpa(ipoint)
       enddo
       vect(imode,krpmax) = one
       call addwfs(vect(1,krpmax),qop(1:3),ichr,imode)
c                                            normalise new operator
       ov    = qpover(qlrout,qlrout)
       call multwf(qlrout,one/sqrt(ov))
c
       write(6,'(a,i3,a)')
     &    ' q-operator nr.',imode,' has been orthogonalised.'
      endif
c      write(14,'(a,i4)') ' after IORTH: IFGRHO=',ifgrho
c                                            store orthogon. q-operator
      qop(4:5 ) = ichr(imode)
      call buffql(qlrout,'o',qop)
c                                            compute pointer names
      qop(4:5 ) = ichr(imode)
      pop(2:6 ) = ichr(imode)
      qbop(5:6) = ichr(imode)
      if(ifsepi.ne.2) then
        callop = qbop
      else
ctest        call prjshl(1,qlrout)
        call commvq(ovpqb,ovpqb0,ovqbqb,ovqbb0,1,resfac)
        yop(5:6) = ichr(imode)
        call buffql(plrout,'o',yop)
        xop(5:6) = ichr(imode)
        callop   = xop
      endif
c                                            set level of analysis in 'moment'
      ifanly = ifanlr
c                                                  compute & store p, [h,p]
      call moment(ovqq,ovqhq,ovqp,ovpp,ovpqb,ovqh0q,qqinv,widrpa,
     &            qrlpr,qrlne,resfac,cmfac,
     &            eshifq,eshift,eshifi,
     &            qop,pop,callop,ipwgt,ifanly)
c      write(14,'(a,i4)') ' after MOMENT: IFGRHO=',ifgrho
c                                            save transition density
ctest      if(ifgrho.eq.4 .and. secrun) then
      if(ifgrho.eq.4) then
        write(17,rec=imode) (prho(i1),i1=1,ngrid2),
     &          (sigma(i2),i2=1,ngrid2),
     &          (qcur(ir,1),ir=1,ngrid2),
     &          (qcur(ir,2),ir=1,ngrid2)
      endif
c                                            save multipole transition
      qtrapr(imode) = qrlpr
      qtrane(imode) = qrlne
c                                            save diagonal elements
c     ipoint = (imode*(imode+1))/2
c     brpa(ipoint) = ovqhq
c     vec2(imode,imode) = brpa(ipoint)
c     crpa(ipoint) = ovpqb
c     qpinv(imode,imode) = ovqp+ovqp
c                                            separable kappa_inv matrix
      if(ifsepi.eq.2) then
        xop(5:6 ) = ichr(imode)
        yop(5:6 ) = ichr(imode)
        write(6,'(2a)') ' XOP=',xop
        call buffql(qlrout,'i',xop)
        call buffql(plrout,'i',yop)
        do imod2=1,imode
          pop(2:3 ) = ichr(imod2)
          call buffql(savlr,'i',pop)
          kappa(imode,imod2) = 2.0*qpover(qlrout,savlr)
          kappa(imod2,imode) = kappa(imode,imod2)
          qop(4:5) = ichr(imod2)
          call buffql(savlr,'i',qop)
          eta(imode,imod2) = 2.0*qpover(plrout,savlr)
          eta(imod2,imode) = eta(imode,imod2)
        enddo
        call buffql(qlrout,'i',qop)
        call buffql(plrout,'i',pop)
        ipoint = (imode*(imode+1))/2
        brpa(ipoint) = -qpover(qlrout,plrout)
      else
c                                            now all stuff for std. RPA
c                                            mass overlaps
        call buffql(plrout,'i','[h,q] ')
        call buffql(savlr,'i','[h0,q]')
        ipoint = (imode*(imode-1))/2
        do imod2=1,imode
          ipoint = 1+ipoint
          qop(4:5 ) = ichr(imod2)
          call buffql(qlrout,'i',qop)
          brpa(ipoint) = -qpover(qlrout,plrout)
          if(secrun) brpa0(ipoint) = -qpover(qlrout,savlr)
        enddo
c                                            spring overlaps
        call buffql(qlrout,'i','[h,p] ')
        call buffql(savlr,'i','[h0,p]')
        ipoint = (imode*(imode-1))/2
        do imod2=1,imode
          ipoint = 1+ipoint
          pop(2:3 ) = ichr(imod2)
          call buffql(plrout,'i',pop)
          crpa(ipoint) = qpover(qlrout,plrout)
          if(secrun) crpa0(ipoint) = qpover(plrout,savlr)
        enddo
      endif
c      write(14,'(a,i4)') ' after overlaps: IFGRHO=',ifgrho
c
c     end of explicit loop '30', terminates by jump to '39'.
c
      goto 30
   39 continue
c
c     read frequency window for plot of strength
c     and prepare output file for strength
c
      if(ifsepi.ne.2 .and. ifgrho.eq.0) then
        widlow = 0.2
        widmin = 8D0  ! 10D0
        widslop = 5D0
        read(5,*) omega1,omega2,delomg,widomg,isosep
cold        open(61,file='strength_full.res')
        open(61,file='lsmooth.inp')
c     &  "'"//title//"'",
        write(61,'(a/6g13.5,2x,a)')
     &  title,
     &  omega1,omega2,delomg,widlow,widmin,widslop,' windor, widths'
c     &  omega1,omega2,delomg,widomg, 'e_init,e_end,e_step,width'
c        if(isosep.eq.1) then
c          write(61,'(a)')
c     &      '  N   X             N             N    Y  '
c        else
          write(61,'(a)')
     &      '  N   X             N       N   Y  '
c        endif
      endif
c
c     optionally separable RPA
c
      if(ifsepi.eq.2) then
        call seprpa(kappa,eta,maxmod,xop,yop,ichr,resfac)
        return
      endif 
c      write(14,'(a,i4)') ' after LSMOOTH: IFGRHO=',ifgrho
c
c     optionally save the mass and spring overlaps
c
      if((ifgrho.eq.2 .or. ifgrho.eq.3) .and. secrun) then
       maxmat = (maxmod*maxmod+maxmod)/2
       write(17,rec=1) (brpa(i),i=1,maxmat)
       write(17,rec=2) (crpa(i),i=1,maxmat)
c       write(17,rec=1) (brpa0(i),i=1,maxmat)
c       write(17,rec=2) (crpa0(i),i=1,maxmat)
      endif
c
c     the norm overlap
c
      do imode=1,maxmod
        pop(2:3 ) = ichr(imode)
        call buffql(plrout,'i',pop)
        do imod2=1,maxmod
          qop(4:5 ) = ichr(imod2)
          call buffql(qlrout,'i',qop)
          qpinv(imod2,imode) = 2.0*qpover(qlrout,plrout)
        enddo
      enddo
      if(maxmod.eq.1) qpisav = qpinv(1,1)
c                                            inverse of norm matrix
      call  matin2 (qpinv,maxmod,krpmax,nerror,determ)
      write(6,'(a,i4,a,g12.4)')
     &    ' norm inversion finished with: nerror=',nerror,
     &    '  determ=',determ
c
c     print the matrices
c
      call primat(brpa,maxmod,' the mass matrix "brpa":')
      call primat(crpa,maxmod,' the spring matrix "crpa":')
      if(secrun) then
       call primat(brpa0,maxmod,' the mass matrix "brpa0":')
       call primat(crpa0,maxmod,' the spring matrix "crpa0":')
      endif
      write(6,'(/a)') ' the inverse norm matrix :'
      do imode=1,maxmod
        write(6,'(1x,5g12.4)') (qpinv(imod2,imode),imod2=1,maxmod)
      enddo
      if(maxmod.eq.1) then
       write(14,'(a,4g12.3)')
     &  ' <q[h,q]>,<p[h,p]>,<[q,p]>,<[q,p]>**-1:',
     &  brpa(1),crpa(1),qpisav,qpinv(1,1)
       write(14,'(a,2g12.3)') ' qtra:',qtrapr(1),qtrane(1)
      endif
c
c     diagonalize rpa-matrices from full hamiltonian
c
      call rpdiag(brpa,crpa,broot,croot,
     &            vect,vec2,alpha,qpinv,workrpa,
     &            cutmod,maxmod,maxout)
      call prirpa(vec2,qtrapr,qtrane,croot,maxmod,maxout)
ctest
      if(ifgrho.eq.4) then
        if(iforth.eq.0) then
          call rtrans(vec2,croot,alpha,maxout,maxmod)
        else
          stop "IFGRHO=4 requires IFORTH=0"
        endif
      endif
c
c     use transformation on alpha and vec2 to compute new diagonal modes
c     explicitely
c
      if(ifanlr.ge.00 .or. .not.secrun) then
        do imode=1,maxout
          if(croot(imode).gt.ecut) goto 699
          call addwfs(alpha(1,imode),qop(1:3),ichr,maxmod)
          qrop(5:6) = ichr(imode)
          call buffql(qlrout,'o',qrop)
          call addwfs(vec2(1,imode),pop(1:1),ichr,maxmod)
          prop(5:6) = ichr(imode)
          call buffql(qlrout,'o',prop)
        enddo
  699   continue
        maxout = imode-1
      endif
c      write(14,'(a,i4)') ' after RTRANS: IFGRHO=',ifgrho
c
c     end of explicit loop '330'. reset counter.
c
      if(.not.secrun) then
       secrun = .true.
       resfac = resfcs
       maxmod = maxout
       if(maxmod.le.0) stop ' no mode left over'
       goto 330
      endif
c
c     finish file for strength distribution
c     (print isoscalar and isovector using temporary files as storage)
c
      if(ifsepi.ne.2 .and. ifgrho.eq.0) then
        close(61)
        open(unit=2,status='unknown',file='strength_full.res')
        call lsmoot(widomg)
        open(61,file='lsmooth.inp')
        open(616,status='scratch',form='unformatted')
        do ispool=1,9999
          read(61,'(a)',err=9199,end=9199) atext
          write(616) atext
          ispoolmax = ispool
        enddo
        stop ' file count exhausted while storing "lsmooth.inp"'
 9199   continue
        rewind(61)
        rewind(616)
        widlow = 0.2
        widmin = 8D0 ! 10D0
        widslop = 5D0
        write(61,'(a/6g13.5,2x,a)')
     &  title,
     &  omega1,omega2,delomg,widlow,widmin,widslop,' windor, widths'
          write(61,'(a)')
     &      '  N   X             N   N N  N         Y  '
        read(616) atext
        read(616) atext
        read(616) atext
        do ispool=4,ispoolmax
          read(616) atext
          write(61,'(a)') atext
        enddo
        close(616)
        close(61)
        write(2,'(a,5x)') '     '
        write(2,'(a,5x)') '     '
        call lsmoot(widomg)
        close(2)
      endif
c      write(14,'(a,i4)') ' after writing LSMOOTH: IFGRHO=',ifgrho
c
c     compute variance matrices from full hamiltonian
c
      if(maxmod.le.10)
     & call cormhf(qqhf,pphf,qqrpa,alpha,vec2,qpinv,
     &             maxmod,maxout,qop,pop,ichr)
      call cormrp(qqrpa,pprpa,alpha,vec2,qpinv,croot,
     &            maxmod,maxout)
c      write(14,'(a,i4)') ' after CORMRP: IFGRHO=',ifgrho
c
c     diagonalize rpa-matrices from mean-field hamiltonian
c     and compute variances
c
      if(ifgrho.ne.0) then
       write(14,*) ' start GSC branch'
       call rpdiag(brpa0,crpa0,broot0,croot0,
     &             vect,vec2,alpha,qpinv,workrpa,
     &             cutmod,maxmod,maxout)
       if(ifanly.ge.0) then
        call prirpa(vec2,qtrapr,qtrane,croot0,maxmod,maxout)
       else
        do imode=1,maxout
          energy = sqrt(abs(croot0(imode)))
          croot0(imode) = energy
        enddo
       endif
       call cormrp(qqrpa0,pprpa0,alpha,vec2,qpinv,croot0,
     &             maxmod,maxout)
c
c     testprint redundance
c
       if(maxmod.eq.1) then
        redun = ovqp*ovqp/(ovqq*ovpp)
        write(14,'(a,4g12.3)')
     &   ' redundance: <qq>,<qp>,<pp>,cosang=',ovqq,ovqp,ovpp,redun
       endif
c
c     compute ground-state correlations to radii etc
c
       call gsc(qqrpa,qqrpa0,qqhf,pprpa,pprpa0,pphf,
     &          qtrapr,qtrane,maxmod,maxout,
     &          varrpr,var0pr,varrsc,var0sc,varrvc,var0vc,
     &          delsgp,delsgs)
      endif
      if(ifgrho.eq.1)
     & call gscrho(dgsrho,dhhrho,qqrpa,qqrpa0,pprpa,pprpa0,
     &             maxmod,qop,pop,ichr)
      if(ifgrho.eq.2)
     & call gscenr(qqrpa,qqrpa0,pprpa,pprpa0,vect,maxmod,maxout)
c
c     new form of gsc. reformulate widths.
c
      if(ifgrho.ne.0) then
       if(maxmod.eq.1) then
        write(14,'(2(/1x,a,4g12.3))')
     &  ' hf: <pp>,<qq>,<qp>,red.=',pprpa0(1),qqrpa0(1),one/qpinv(1,1),
     &    4.0*pprpa0(1)*qqrpa0(1)/qpinv(1,1)**2,
     &  ' rp: <pp>,<qq>,<qp>,red.=',pprpa(1),qqrpa(1),one/qpinv(1,1),
     &    4.0*pprpa(1)*qqrpa(1)/qpinv(1,1)**2
       endif
       if(ifgrho.eq.3) then
ctest        call gscat2(qqrpa,qqrpa0,pprpa,pprpa0,
ctest     &             maxmod,qop,pop,ichr,resfac,cmfac)
        call gscat(dgsrho,dhhrho,qqrpa,qqrpa0,pprpa,pprpa0,
     &              maxmod,qop,pop,ichr,resfac,cmfac,ifppsc)
       endif
       call quadwd(qqrpa,qqrpa0,maxmod,maxout)
       call quadwd(pprpa,pprpa0,maxmod,maxout)
       if(ifgrho.eq.2.or.ifgrho.eq.3)
     &  call gscenr(qqrpa,qqrpa0,pprpa,pprpa0,vect,maxmod,maxout)
      endif
      if(ifgrho.eq.3)
     & call gscrho(dgsrho,dhhrho,qqrpa,qqrpa0,pprpa,pprpa0,
     &             maxmod,qop,pop,ichr)
ctest      if(maxmod.gt.1) then
ctest       call radius
ctest       write(14,'(/a/9x,a,2(/a,3g12.5,5x,g12.5))')
ctest     &  ' radii after correlation:',
ctest     &  '   r_rms      r_diffr     surface          energy',
ctest     &  ' mass:   ', rms, rdms, surf,etotal,
ctest     &  ' charge: ',rmsc,rdmsc,surfc
ctest      endif
c
c     finish pseudo-loop 340 for gsc and optionally restart
c     with reduced 'resfac'
c
      if(ifgrho.eq.3 .and. resfac.eq.resfcs) then
       resfac = 0.1*resfcs
       goto 340
      endif
      if(ifgrho.eq.3) call prigsc
c
c     print superposition weights
c
      if(ifanly.ge.0) then
        do iw=6,14,8
          write(iw,'(a)')
     &      ' superposition weights for q-operator:'
          do imode=1,maxout
            write(iw,'(10(1pg8.1))') (alpha(i1,imode)**2,i1=1,maxmod)
          enddo
          write(iw,'(a)')
     &      ' superposition weights for p-operator:'
          do imode=1,maxout
            write(iw,'(10(1pg8.1))') (vec2(i1,imode)**2,i1=1,maxmod)
          enddo
        enddo
      endif
c
c     compose transition densities
c
c      if(ifgrho.eq.4) call rtrans(vec2,croot,maxout,maxmod)
c
c      analyse properties of new q-operators
c
c
c      test orthogonality of results
c
      if(ifanlr.ge.1 .and. maxmod.lt.15) then
        do imode=1,maxout
          ipoint = (imode*(imode-1))/2
          do imod2=1,imode
            ipoint = 1+ipoint
            qrop(5:6 ) = ichr(imode)
            prop(5:6 ) = ichr(imod2)
            call buffql(qlrout,'i',qrop)
            call buffql(plrout,'i',prop)
            brpa(ipoint) = qpover(qlrout,plrout)
            vec2(imode,imod2) = brpa(ipoint)
            vec2(imod2,imode) = brpa(ipoint)
          enddo
        enddo
        do iw=6,14,8
          write(iw,'(/a)')
     &       ' check orthogonality of resulting vetcors:'
          do imode=1,maxout
            write(iw,'(a,i2,a,5g12.4,20(/8x,5g12.4))')
     &        ' mode',imode,':',(vec2(imod2,imode),imod2=1,maxmod)
          enddo
        enddo
      endif
c
c      compute energies etc for the new modes
c
      if(ifanlr.ge.0) then
       do 91 imode=1,maxout
       do iwunit=6,14,8
         write(iwunit,'(/a,i3,a)')
     &    ' coupled mode nr.',imode,': '
         write(iwunit,'(a,10(/5x,5(1pg12.4)))')
     &    ' composition weights (modes normalised to unit mass):',
     &    ( (alpha(ip,imode))**2 ,ip=1,maxmod)
       enddo
       qrop(5:5) = ichr(imode)
       prop(5:6) = ichr(imode)
       call buffql(qlrout,'i',qrop)
       call moment(ovqq,ovqhq,ovqp,ovpp,ovphp,ovqh0q,qqinv,widrpa,
     &             qrlpr,qrlne,resfac,cmfac,
     &             eshifq,eshift,eshifi,
     &             qrop,prop,'[h,p] ',-999,ifanlr+1)
       broot(imode)  = ovqp
       qrpwid(imode) = half/widrpa
ctest3
      if(ifgrho.eq.4)
     & write(17,rec=imode) (prho(i1),i1=1,ngrid2),
     &          (sigma(i2),i2=1,ngrid2),
     &          (qcur(ir,1),ir=1,ngrid2),
     &          (qcur(ir,2),ir=1,ngrid2)
c
       if(imode.le.2) then
        call buffql(plrout,'i','[h,p] ')
        call buffql(savlr,'i','[h0,p]')
        ipoint = 2*(imode-1)
        do imod2=1,2
          ipoint = 1+ipoint
          prop(5:6) = ichr(imod2)
          call buffql(qlrout,'i',prop)
          crpa(ipoint) = qpover(qlrout,savlr)
          brpa(ipoint) = qpover(qlrout,plrout)
        enddo
ctest3
        if(imode.eq.2 .and. ifgrho.eq.4) then
         write(14,'(/a/4x,a,2(/1x,a,4g12.4))')
     &    ' q*q-matrix of first two modes:',
     &    '   <p_1 h p_1>   <p_2 h p_1>     <p_1 h p_2>   <p_2 h p_2>',
     &    '  h:',crpa(1),crpa(2),crpa(3),crpa(4),
     &    ' h0:',brpa(1),brpa(2),brpa(3),brpa(4)
c                         compute third functional derivatives
        read(17,rec=1) (prho(i1),i1=1,ngrid2),
     &          (sigma(i2),i2=1,ngrid2),
     &          (qcur(ir,1),ir=1,ngrid2),
     &          (qcur(ir,2),ir=1,ngrid2)
        read(17,rec=2) (aux1(i1),i1=1,ngrid),(aux2(i1),i1=1,ngrid),
     &          (sigma(i2),i2=1,ngrid2),
     &          (qcur(ir,1),ir=1,ngrid2),
     &          (qcur(ir,2),ir=1,ngrid2)
         sum112 = zero
         sum122 = zero
         sum11  = zero
         sum12  = zero
         t3d16  = t3/16.0*(power+2.0)*(power+one)*power*frpirs
         powm1  = power-one
         do ir=2,ngrid
           r2     = (ir-1)*rstep
           r2     = r2*r2
           rs     = (rho(ir)+rho(ir+ngrid))
           rhopow = r2*rs**powm1
           rho1   = prho(ir)+prho(ir+ngrid)
           rho2   = aux1(ir)+aux2(ir)
           sum112 = rhopow*rho1*rho1*rho2 + sum112
           sum122 = rhopow*rho1*rho2*rho2 + sum122
           rhopow = rhopow*rs
           sum11  = rhopow*rho1*rho1 + sum11
           sum12  = rhopow*rho1*rho2 + sum12
         enddo
         write(14,'(4(a,g12.3))')
     &    ' <[h,p1],p1],p2]>=',t3d16*sum112,
     &    ' <[h,p1],p2],p2]>=',t3d16*sum122,
     &    ' <[h,p1]p1]>=',t3d16*sum11,
     &    ' <[h,p1]p2]>=',t3d16*sum12
         endif
       endif
c                                 <q h**-1 q>
       call buffql(qlrout,'i',qrop)
       call invmrq(ovqhmq,ovpqb,ovpqb0,resfac,cmfac)
c                                 <p h**-1 p>
       call buffql(plrout,'i',prop)
       call invmrp(ovphmp,resfac)
c                                 print
       write(6,'(5(a,g12.4))')
     &   ' check: <qhq>=',ovqhq,' <php>=',ovphp,' <qh-1q>=',ovqhmq,
     &   ' <ph-1p>=',ovphmp,' <qp>=',ovqp
       write(14,'(5(a,g11.4))')
     &   ' qhq=',ovqhq,' php=',ovphp,' qh-1q=',ovqhmq,
     &   ' ph-1p=',ovphmp,' qp=',ovqp
c
       enrpa  = sqrt(ovqhq*ovphp)/ovqp
       dele2p = half*enrpa*(ovphmp-one)
       dele2q = half*enrpa*(enrpa*enrpa*ovqhmq-one)
       dele   = sqrt(abs(dele2p+dele2q))
       write(6,'(a,g10.3,a,2g10.3,a,g10.3)')
     &   ' h-1 weight: enrpa=',enrpa,' del-2 q,p=',dele2p,dele2q,
     &   ' del-e=',dele
       write(14,'(a,g10.3,a,2g10.3,a,g10.3)')
     &   ' h-1 weight: enrpa=',enrpa,' del-2 q,p=',dele2p,dele2q,
     &   ' del-e=',dele
       e1sq   = ovqhq/max(ovqhmq,xsmall)
       e3sq   = ovphp/max(ovqhq,xsmall)
       del3   = 0.5D0*sqrt(abs(e3sq-e1sq))
       write(14,'(a,g11.4,a,2g11.4)')
     &   '   del_e3 =',del3,'  e1,e3=',sqrt(abs(e1sq)),sqrt(abs(e3sq))
c
c       call pritra(imode)
   91  continue
c       if(ifgrho.ne.0)
c     &  call gscrh2(dgsrho,dhhrho,qrpwid,broot,maxmod,qrop,prop,
c     &              ichr)
c
      endif
c
   20 continue
   10 continue
c
      if(ifgrho.eq.4 .or. ifgrho.eq.2) close(unit=17)
c
      if(ifgrho.eq.1) then
       isorpa = 1
       call initqr('mom',zero,one*lrpa)
       call d2qrho(varrpr,var0pr,varrsc,var0sc,varrvc,var0vc,
     &         maxmod,ifgrho)
       call d2qenr(varrpr,var0pr,varrsc,var0sc,varrvc,var0vc,
     &         maxmod)
      endif
c
      return
      end
c-----rpdiag--------------------------------------------------------------------
c
      subroutine rpdiag(brpa,crpa,broot,croot,
     &                  vect,beta,alpha,qpinv,work1,
     &                  cutmod,maxmod,maxout)
c
      include 'params.inc'
      parameter(krpmat=krpmax*(krpmax+1)/2)
c
c     diagonalizes the rpa-matrices brpa and crpa with relative
c     norm-overlap given on qpinv.
c
c     the list parameters are:
c      brpa    = mass matrix (input)
c      crpa    = spring matrix (input)
c      broot   = diagonal elements of mass matrix
c      croot   = diagonal elements of spring matrix
c                and square of the final rpa-energies
c      vect    = workspace for diagonalization of the mass matrix
c                and storage for part of the composition of 'p'-operators.
c      beta    = workspace for diagonalization of the spring mass matrix
c                and other intermediate matrix operations;
c                and finally transformation of momenta 'p'
c      alpha   = transformation of coordinates 'q'
c      qpinv   = inverse conjugation matrix
c      work1   = linear workspace
c      maxmod  = maximum number of initial ph modes
c      maxout  = maximum number of final rpa modes
c
      dimension   brpa(krpmat),crpa(krpmat),
     &            broot(krpmax),croot(krpmax),
     &            vect(krpmax,krpmax),
     &            beta(krpmax,krpmax),
     &            alpha(krpmax,krpmax),
     &            qpinv(krpmax,krpmax),
     &            work1(krpmax)
c
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
c     save mass matrix on 'vec2' for tests
c
      do imode=1,maxmod
        ipoint = (imode*(imode-1))/2
        do imod2=1,imode
          ipoint = 1+ipoint
          beta(imode,imod2) = -brpa(ipoint)
          beta(imod2,imode) = -brpa(ipoint)
        enddo
      enddo
c
c     diagonalise mass matrix
c
      call givens(brpa,broot,vect,maxmod,maxmod,krpmax)
      write(6,'(a,10(/1x,5g12.4))')
     &    ' roots of brpa:',(-broot(i),i=1,maxmod)
      imomax = 0
      do imode=1,maxmod
        broot(imode) = -broot(imode)
        if(broot(imode)/broot(1).lt.cutmod .and. imomax.eq.0)
     &     imomax = imode-1
        write(6,'(a,i2,a,20(/1x,5g12.4))')
     &    ' vect(*,',imode,') :',(vect(ip,imode),ip=1,maxmod)
      enddo
      if(imomax.ne.0) then
       maxout = imomax
       write(6,'(a,i5)') ' cutoff at mode nr.',maxout
       write(14,'(a,i5)') ' cutoff at mode nr.',maxout
      else
       maxout = maxmod
      endif
c
c     scale the transformation
c
      do imode=1,maxout
        xnorm  = sqrt(broot(imode))
        xnormi = one/xnorm
        do imod2=1,maxmod
          alpha(imod2,imode) = xnormi*vect(imod2,imode)
          vect(imod2,imode)  = xnorm*vect(imod2,imode)
        enddo
      enddo
c
c     test the transformation
c
      write(6,'(/a)') ' test the mass matrix after transformation:'
      xnorm  = zero
      do imode=1,maxout
        do imod2=1,maxout
          sum    = zero
          do i1=1,maxmod
          do i2=1,maxmod
            ipoint = (i1*(i1-1)/2) + i2
            sum    = sum+alpha(i1,imode)*alpha(i2,imod2)*beta(i1,i2)
          enddo
          enddo
          if(imode.eq.imod2) then
           xnorm  = (one-sum)**2 + xnorm
          else
           xnorm  = sum*sum +xnorm
          endif
          work1(imod2) = sum
        enddo
        write(6,'(1x,10f8.5)')  (work1(imod2),imod2=1,maxout)
      enddo
      xnorm  = sqrt( xnorm/(maxout*maxout) )
      write(14,'(a,g12.4)')
     &  ' average deviation of diag. mass matrix=',xnorm
      write(6,'(a,g12.4)')
     &  ' average deviation of diag. mass matrix=',xnorm
c
c     combine transformations for 'p'
c
      do imode=1,maxout
        do imod2=1,maxmod
          work1(imod2) = vect(imod2,imode)
        enddo
        do imod2=1,maxmod
          sum    = zero
          do i1=1,maxmod
            sum    = qpinv(imod2,i1)*work1(i1) + sum
          enddo
          vect(imod2,imode) = sum+sum
        enddo
      enddo
c
c     transform spring matrix
c
      ipoint = 0
      do imode=1,maxout
      do imod2=1,imode
        ipoint = 1+ipoint
        sum    = zero
        ip     = 0
        do i1=1,maxmod
        do i2=1,i1
          ip     = 1+ip
          sum    = sum+vect(i1,imode)*vect(i2,imod2)*crpa(ip)
          if(i1.ne.i2) sum = sum+vect(i2,imode)*vect(i1,imod2)*crpa(ip)
        enddo
        enddo
        brpa(ipoint)  = sum
      enddo
      enddo
      ipmax  = (maxout*(maxout+1))/2
      do ipoint=1,ipmax
        crpa(ipoint) = brpa(ipoint)
      enddo
c
c     diagonalise spring matrix
c
      call givens(crpa,croot,beta,maxout,maxout,krpmax)
      write(6,'(a,10(/1x,5g12.4))')
     &    ' roots of crpa:',( croot(i),i=1,maxmod)
      write(6,'(a,5(/1x,5g12.4))')
     &   ' rpa-energies:',
     &   (sqrt(abs(croot(imode))),imode=1,maxout)
c
c     compose final transformation for 'q'
c
      do imod2=1,maxmod
        do imode=1,maxout
          sum    = zero
          do i1=1,maxout
            sum    = beta(i1,imode)*alpha(imod2,i1) + sum
          enddo
          work1(imode) = sum
        enddo
        do imode=1,maxout
          alpha(imod2,imode) = work1(imode)
        enddo
      enddo
      do imode=1,maxout
        write(6,'(a,i2,a,20(/1x,5g12.4))')
     &    ' alpha(*,',imode,') :',(alpha(ip,imode),ip=1,maxmod)
      enddo
c
c     compose final transformation for 'p'
c
      do imode=1,maxout
        do imod2=1,maxmod
          sum    = zero
          do i1=1,maxout
            sum    = beta(i1,imode)*vect(imod2,i1) + sum
          enddo
          work1(imod2) = sum
        enddo
        do imod2=1,maxmod
          beta(imod2,imode) = work1(imod2)
        enddo
      enddo
      do imode=1,maxout
        sum = 0D0
        do ip=1,maxmod
          sum = beta(ip,imode)**2 + sum
        enddo
        write(6,'(a,i3,a)') ' beta(*,',imode,'):'
        do ip=1,maxmod
          write(6,'(i5,5(1pg12.4))')
     &     ip,beta(ip,imode),beta(ip,imode)**2/sum
        enddo
        write(6,'(1x)')
      enddo
c
      return
      end
c     include 'givens.f'
c-----specor--------------------------------------------------------------------
c
      subroutine specor
c
      include 'params.inc'
c
c     polarisation correction to single particle energies
c
      dimension speadd(ksl2),spedir(ksl2),specr(ksl2)
c
      include 'skypar.inc'
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'numpar.inc'
      include 'constn.inc'
      include 'linres.inc'
c
      data xminlw/1.e-20/
      data ifop14/0/
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
c     read general options
c
      cmfac  = zero
      read(05,*)   lmax
      read(05,*)   itrpa,epsrpa,dabrpa
      mlmax  = lmax+lmax+1
c
c     open file 14 for concentrated output of matrix elements
c
      if(ifop14.eq.0) then
       open(unit=14,status='unknown',form='formatted',file='for014')
       ifop14 = 1
      endif
c
c     write a header
c
      write(14,'(a/1x,2a/1x,2(a,i4),5x,2(a,i4))')
     &  '1',
     &  '  title:',title,
     &  'nucleus:  nprot=',nprot,'  nneut=',nneut,
     &  ' shells:  nmaxp=',nmaxp,'  namxn=',nmaxn
      write(14,'(1x,a,5g13.6/7x,a,5g13.6/7x,a,g13.6/)')
     &  'force: t0-t4=',t0,t1,t2,t3,t4,
     &        ' x0-x3=',x0,x1,x2,x3,b4p,
     &        ' power=',power
      write(14,'(1x,2(a,f10.3)/18x,2(a,f10.3))')
     &  'some properties: rms-radius=',rms,'    diffr.radius=',rdms,
     &  '   surface=',surf,'          energy=',etotal
c
c     header
      write(6,'(a)') ' running specor:'
      write(14,'(a)') ' running specor:'
      write(6,'(21x,a,i5,a,f7.3,a,g12.4)')
     &    '  itrpa=',itrpa,' epsrpa=',epsrpa,' daprpa=',dabrpa
      resfac = one
c
c     loop over the states.
c     states are read in line by line.
c     loop is terminated if state-nr.<0 is encountered
c
      do 10 istate=1,nmaxto+1
      read(05,*)   na,iso
      if(na.le.0) return
      write(6,'(/a,i5,a,i5)') ' now: na=',na,'  iso=',iso
c
      if(iso.eq.1) then
       nmax   = nmaxp
       nstack = 0
       isoadd = 0
      else
       nmax   = nmaxn
       nstack = nmaxp*ngrid
       isoadd = ngrid
      endif
c
c     check bounds
c
      if(na.gt.nmax) then
       write(14,'(a,i5,a)') ' state nr.',na,' too large.'
       goto 10
      endif
c
      jsa   = jsp(na)
      la    = jsa/2
      xla   = la
      xja   = (jsa-1)/2+0.5D0
      wgt   = weight(na+nstack/ngrid)
      nadd  = nstack + (na-1)*ngrid
      subfac = spenrg(na+nstack/ngrid)
c
c     sub-header for compressed output
c
      write(14,'(/2(a,i3),a,f4.1,i3,a/a,a)')
     &    ' for iso=',iso,' & na=',na,' with ja,la=',xja,la,' :',
     &    '  lrpa   direct   response    core      total      spe',
     &    '    spe0'
c
c     loop over polarisation channels
c
      coradd = zero
      do nb=1,nmax
        specr(nb)  = zero
      enddo
      ja2    = xja+xja
      mlrpmx = min(4*la+1,ja2+ja2+1)
      if(mlrpmx.gt.klrpa .or. mlrpmx.gt.mlmax) then
       mlrpmx = min(klrpa,mlmax)
       write(6,'(/a,i3,a,a,i3/a)')
     &    ' the state nr.',na,' requires more polarisation channels',
     &    ' than provided by klrpa=',klrpa,'.',
     &    ' check convergence!'
      endif
      do 30 mlrpa=1,mlrpmx,4
      lrpa   = mlrpa/2
      alrpa  = lrpa
      call setred
      isorpa = 1
      write(6,'(a,i3)') ' run specor with lrpa=',lrpa
c
c     compute density for state (iso,na) on 'prho','ptau',...
c     and the according hartree-fock potentials on 'rpot',...
c
      call densna(iso,na,prho,ptau,psj,psjs,dprho,d2prho)
c
      redmat = redume(xja,alrpa,xja)
      wigfac = wigner(xja,alrpa,xja,xja,zero,-xja)
      wigfa2 = wigner(xja,alrpa,xja,xja-one,zero,one-xja)
c      write(6,'(2(a,g12.4),g12.4)') ' testprint: redmat=',redmat,
c     &    '  wigner=',wigfac,wigfa2
      rfac   = redmat*wigfac
c      write(6,'(a,20(/1x,5g12.4))') ' testprint prho:',
c     &    (rfac*prho(ir),ir=1,ngrid)
c      write(6,'(4x,5g12.4))') (rfac*prho(ir),ir=ngridp,ngrid2)
c      write(6,'(a,20(/1x,5g12.4))') ' testprint psjs:',
c     &    (rfac*psjs(ir),ir=1,ngrid)
c      write(6,'(4x,5g12.4))') (rfac*psjs(ir),ir=ngridp,ngrid2)
      call ppot(resfac,cmfac)
c      write(6,'(a,20(/1x,5g12.4))') ' testprint ppot:',
c     &    (rfac*rpot(ir),ir=1,ngrid)
c      write(6,'(4x,5g12.4))') (rfac*rpot(ir),ir=ngridp,ngrid2)
c
c     initialise bare perturbed wave functions on 'qlrout'
c
      call delmf0
      write(6,'(a,g12.4)') ' testprint |qlrout|**2=',qnorm(qlrout)
c
c     compute direct contribution to s.p.energy
c
      do nb=1,nmax
        jsb   = jsp(nb)
        lb    = jsb/2
        xlb   = lb
        xjb   = (jsb-1)/2+0.5D0
        wgt   = weight(nb+nstack/ngrid)
        naddb = nstack + (nb-1)*ngrid
        call hrespp(iso,lb,xjb,naddb,xjb,xlb,lrpa+1,
     &              wf(naddb+1),wfs(naddb+1))
        wigfab = wigner(xjb,alrpa,xjb,xjb,zero,-xjb)
        spesum = zero
        r      = zero
        do ir=2,ngrid
          r      = rstep+r
          spesum = wf(naddb+ir)*hlinrs(ir)*r + spesum
        enddo
        spedir(nb) = rstep*wigfab*spesum
      enddo
c
c     compute response to perturbed mean-field
c
      ovqp   = zero
      call invmrq(ovqpm1,ovpqb,ovpqb0,resfac,cmfac)
      call pdensi(resfac)
c      write(6,'(a,20(/1x,5g12.4))') ' testprint prho:',
c     &    (rfac*prho(ir),ir=1,ngrid)
c      write(6,'(4x,5g12.4))') (rfac*prho(ir),ir=ngridp,ngrid2)
c      write(6,'(a,20(/1x,5g12.4))') ' testprint psjs:',
c     &    (rfac*psjs(ir),ir=1,ngrid)
c      write(6,'(4x,5g12.4))') (rfac*psjs(ir),ir=ngridp,ngrid2)
      call ppot(resfac,cmfac)
c      write(6,'(a,20(/1x,5g12.4))') ' testprint ppot:',
c     &    (rfac*rpot(ir),ir=1,ngrid)
c      write(6,'(4x,5g12.4))') (rfac*rpot(ir),ir=ngridp,ngrid2)
c
c     compute expectation value of polarisation hamiltonian
c     with s.p. state 'nb'
c
      do nb=1,nmax
        jsb   = jsp(nb)
        lb    = jsb/2
        xlb   = lb
        xjb   = (jsb-1)/2+0.5D0
        wgt   = weight(nb+nstack/ngrid)
        naddb = nstack + (nb-1)*ngrid
        wigfab = wigner(xjb,alrpa,xjb,xjb,zero,-xjb)
        call hrespp(iso,lb,xjb,naddb,xjb,xlb,lrpa+1,
     &              wf(naddb+1),wfs(naddb+1))
        spesum = zero
        r      = zero
        do ir=2,ngrid
          r      = rstep+r
          spesum = wf(naddb+ir)*hlinrs(ir)*r + spesum
        enddo
        spesum = rstep*wigfab*spesum
        speadd(nb) = -(spesum+spesum)
        specr(nb)  = speadd(nb)+spedir(nb)+specr(nb)
      enddo
c
      write(6,'(a,i3,a,f4.1,i3,a,2i3/a,10(/1x,5g12.4))')
     &    ' for na=',na,' with ja,la=',xja,la,' and for mode iso,l=',
     &    isorpa,lrpa,' response pol. spe:',(speadd(nb),nb=1,nmax)
      write(6,'(a,10(/1x,5g12.4))')
     &    ' direct pol. spe:',(spedir(nb),nb=1,nmax)
      coradd = ovqpm1 + coradd
      write(6,'(a,g12.4)')
     &    ' core-response  :',ovqpm1
      write(14,'(1x,i5,3f10.3)') lrpa,spedir(na),speadd(na),ovqpm1
c
   30 continue
      write(6,'(a,a,i3,a,f4.1,i3)')
     &    ' total polarisation energy',
     &    ' for na=',na,' with ja,la=',xja,la
      write(6,'(5g12.4)')  (specr(nb),nb=1,nmax)
      cor1   = specr(na)+coradd
      cor2   = specr(na)+coradd+subfac
      write(6,'(a,i4/10x,3(a,g12.4))')
     &    ' incl. core response for na=',na,' coradd=',coradd,
     &    ' cor.spe=',cor1,'  spe=',cor2
      write(14,'(a,5x,f10.3,5x,4f10.3)')
     &  ' total',specr(na),coradd,cor1,cor2,subfac
c
   20 continue
   10 continue
c
      return
      end
c-----densna---------------------------------------------------------------
c
      subroutine densna(iso,na,arho,atau,asj,asjs,ad1r,ad2r)
c
      include 'params.inc'
c
c     computes for a nucleon in state (iso,na)
c     the density 'arho', the kinetic density 'atau' and
c     the spin-orbit density 'asj' from the radial single particle
c     wave functions given on 'wf' and their derivatives given on
c     'wfs'.
c
c     status: 20. dec. 1988, p.-g.r.
c
      dimension arho(kgr2),atau(kgr2),asj(kgr2,2),asjs(kgr2),
     &          ad1r(kgr2),ad2r(kgr2)
c
      include 'skypar.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      data omlow/1.0D-6/
c
c     zero fields before accumulation
c
      do i=1,ngrid+ngrid
        arho(i)  = zero
        atau(i)  = zero
        asjs(i)  = zero
        asj(i,1) = zero
        asj(i,2) = zero
      enddo
c
c     derived pointers for given 'iso'
c
      ie     = ngrid*iso
      im     = ie-ngrid+1
      ia     = im+1
      if(iso.eq.1) then
       nmax   = nmaxp
       nstack = 0
       ngw    = 0
       isoadd = 0
      else
       nmax   = nmaxn
       nstack = nmaxp*ngrid
       ngw    = nmaxp
       isoadd = ngrid
      endif
c
c     derived quantum numbers for given 'na'
c     loop over states for given isospin.
c      'na' counts the states in harmonic oscillator ordering.
c      'js' is the according composed ang.mom.,
c      'll' is the orbital ang.mom.,
c      'wl' is the spin-orbit splitting,
c      'om' is the actual phase-space factor
c
      js     = jsp(na)
      xja    = (js-1)/2+0.5D0
      ll     = js/2
      xla    = ll
      xlp    = float(ll*ll+ll+1)
      wl     = wls(js)
      om     = redume(xja,alrpa,xja)
     &        * wigner(xja,alrpa,xja,xja,zero,-xja)
      if(weight(na+ngw).gt.half) om    = -om
      wigfac = wigner(xja,alrpa,xja,xja,zero,-xja)
      if(weight(na+ngw).gt.half) wigfac = -wigfac
      write(6,'(a,g12.4)')
     &   ' in densna:   redmat=',om
c                                     to be divided by 'mlrpa' ???
      nadd  = nstack + (na-1)*ngrid
c                                     adapt 'nadd' to loops from 'ia'-'ie'
      nadd  = nadd - isoadd
      if(abs(om).lt.omlow) goto 30
c                         the case r=0.0 needs special treatment
      if(ll.lt.1) then
       x      = wfs(im+nadd)
       arho(im) = (x*x)*om + arho(im)
      elseif(ll.eq.1) then
       x      = (16.0*wf(ia+nadd)-wf(ia+nadd+1))/(6.0*rstep*rstep)
       atau(im) =.75*om*(x*x) + atau(im)
      endif
c                         now all r>0.0
      r      = 0.0D0
      do  i=ia,ie
        r      = r+rstep
        ri     = one/r
        phas   = (ri*ri)*om
c
        wfa    = wf(nadd+i)
        wfa2   = wfa*wfa
        arho(i) = (phas*wfa2) + arho(i)
c
        wfas   = wfs(nadd+i)
        atau(i) = (wfas*wfas-((wfa*wfas)+(wfa*wfas)-xlp*wfa2*ri)*ri)
     &        *phas + atau(i)
      enddo
c     write(6,'(a,20(/1x,5g12.4))')
c    & ' in densna: prho=',(arho(i),i=1,ngrid)
c
   30 continue
c
c     some derivatives: d/dr 'arho' --> 'ad1r'
c                       d/dr**2 'arho' + 2/r d/dr 'rho' --> 'ad2r'
c                       d/dr 'asj' + 2/r 'asj' --> 'asjs'
c     'im' is the initial position for given 'iso'.
c
      zwh = 2.0/rstep
      call derv15(arho(im),ad1r(im), one, ngrid,rstep)
      call derv25(arho(im),ad2r(im), one, ngrid,rstep)
      ad2r(im) = ad2r(im) + (8.0*ad1r(im+1)-ad1r(im+2))/(3.0*rstep)
      do i=im+1,ie
        ad2r(i) = ad2r(i) + zwh*ad1r(i)/(i-im)
      enddo
c
c     accumulate the spin-orbit density 'psj'
c
      if(t4.eq.zero) return
      do ilamp=1,2
        lamp   = lrpa-3+ilamp+ilamp
        if(lamp.ge.0) then
          alamp  = lamp
          call redspo(xja,xla,alrpa,alamp,xja,xla,facder,facrad)
c          facder = facder
c          facrad = facrad
          facder = -facder  !/sq2   !*sq2 
          facrad = -facrad  !/sq2   !*sq2 
          do ir=1,ngrid
            rinv   = rinvs(ir)
            iro    = isoadd+ir
            ira    = nadd+iro
            asj(iro,ilamp) = wigfac*wf(ira)*rinv*rinv
     &                      * (facder*wfs(ira) + facrad*wf(ira)*rinv)
     &                      + asj(iro,ilamp)
          enddo
        endif
      enddo
      parrho = (1-mod(lrpa,2)*2)
      call extrap(asj(isoadd+1,1),-parrho)
      call extrap(asj(isoadd+1,2),-parrho)
c
c     compute 'asjs' from 'asj' as divergence
c
      do ilamp=1,2
        lamp   = lrpa-3+ilamp+ilamp
        if(lamp.ge.0) then
          alamp  = lamp
          if(ilamp.eq.1) then
           facder = sqrt(alrpa/mlrpa)
           facrad = -facder*(alrpa-one)
          else
           facder = -sqrt(alamp/mlrpa)
           facrad = facder*(alamp+one)
          endif
          call derv15(asj(isoadd+1,ilamp),work1(1),-parrho,ngrid,rstep)
          do ir=2,ngrid
            irn    = isoadd+ir
            asjs(irn) = facder*work1(ir)+facrad*asj(irn,ilamp)*rinvs(ir)
     &                 + asjs(irn)
         enddo
        endif
      enddo
      call extrap(asjs(isoadd+1),parrho)
c
c  20 continue
c
      write(6,'(a,50(/1a,5g12.4))')
     &  ' Nabla*J_ls protons:',(asjs(ir),ir=1,ngrid)
      write(6,'(a,50(/1a,5g12.4))')
     &  ' Nabla*J_ls neutrons:',(asjs(ir+ngrid),ir=1,ngrid)
c
      return
      end
c-----delmf0---------------------------------------------------------------
c
      subroutine delmf0
c
      include 'params.inc'
c
c     action of the mean-field hamiltonian on 'pot','bmass',...
c     q(lm) on-the s.p wave function
c
      include 'wavfun.inc'
      include 'result.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
c
c--------------------------------------------------------------------
c
      ovqq   = zero
c
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     now all quantum numbers are provided store qlrout.
c
      call hrespp(iso,la,xja,nadd,xj,xl,jdif,
     &            wf(nadd+1),wfs(nadd+1))
      redmat = redume(xj,alrpa,xja)
      do ir=2,ngrid
        qlrout(nadd+ir,jdif) = redmat*hlinrs(ir)
      enddo
      parity = 1-2*mod(la+lrpa,2)
      call extrap(qlrout(nadd+1,jdif),parity)
c
c     project out all wave function with same 'js'
c
      call projec(qlrout(nadd+1,jdif),js,nstack,nmax)
c
      ovqq   = overlp(qlrout(nadd+1,jdif),qlrout(nadd+1,jdif),
     &                rstep,ngrid)   +  ovqq
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
      return
      end
c-----qinit ------------------------------------------------------------------
c
      subroutine qinit(typmod,jselec,isplit,ihwgt,ih0wgt,ipwgt,
     &                 iprshl,imode,ifprim,iflagr,xkqrd,resfac,cmfac,
     &                 nparsv,nholsv,isophs,toper)
c
      include 'params.inc'
c
c     initialization of the q-operator according to options given in
c     the list; these are the parameters read in from file 5; for
c     explanation see haeder of 'sumrul'.
c     there are furthermore the list-parameters:
c      'iflagr' = read-flag, the routine returns the number of 'automatic'
c                 modes (from 'ph=all' or 'jj=all+n') if this option is
c                 encountered. reading from file 5 is then omitted until
c                 all 'automatic' modes are worked up.
c      'nparsv' = field for storing particle states in 'automatic'
c                 selection of ph-modes:
c      'nholsv' = field for storing hole states in 'automatic'
c                 selection of ph-modes:
c
      character*7  afix
      character*16 typmod
      character*12 jselec
      character*50 toper
      dimension nparsv(krpmax),nholsv(krpmax),isophs(krpmax)
c
      include 'wavfun.inc'
      include 'result.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
      include 'numpar.inc'
c
c--------------------------------------------------------------------
c
c     trace end condition;
c     count number of mode imode and count down read-flag iflagr
c
      if(typmod(1:3).eq.'end')  then
       iflagr = -999
       return
      else
       imode = 1+imode
      endif
c
      open(unit=3,status='scratch',form='formatted')
      rewind 3
c
c     analyse typmod
c
      lentyp = len(typmod)
      iposeq = index(typmod,'=')
      write(6,'(1x,a)') typmod
      if(iposeq.le.0) stop ' invalid mode input typmod'
      if(typmod(1:3).eq.'mom' .or. typmod(1:3).eq.'bes'
     &   .or. typmod(1:3).eq.'loc' .or. typmod(1:3).eq.'mod'
     &   .or. typmod(1:3).eq.'bed' .or. typmod(1:3).eq.'moj'
     &   .or. typmod(1:3).eq.'mot' .or. typmod(1:3).eq.'bej'
     &   .or. typmod(1:3).eq.'msc') then
       write(3,'(a)') typmod(iposeq+1:lentyp)
       rewind 3
       read(3,*) xkqrd
       rewind 3
      elseif(typmod(1:3).eq.'ph=') then
       write(3,'(2a)') typmod(4:lentyp),','
       rewind 3
       read(3,*) npart,nhole,isoini
       rewind 3
      elseif(typmod(1:3).eq.'pha') then
       if(iflagr.eq.0) then
        lencom = index(typmod,',')
        if(lencom.le.0) stop ' no comma in ph-option of typmod'
        write(3,'(2a)') typmod(5:lentyp),','
        write(14,'(2a)') typmod(5:lentyp),','
        rewind 3
        read(3,*) ecutlw,ecutup
        rewind 3
        call setph(nend,ecutlw,ecutup,nparsv,nholsv,isophs)
        iflagr = nend
        if(imode+nend-1.gt.krpmax) stop ' too many ph-states'
       endif
       npart  = nparsv(iflagr)
       nhole  = nholsv(iflagr)
       isoini = isophs(iflagr)
      else
       stop ' unknown typmod'
      endif
c
c     analyze 'jselec'
c
      if(jselec(1:3).eq.'jj=') then
       write(3,'(a)') jselec(4:lentyp)
       rewind 3
       read(3,*) jspart,nhole,isoini
       rewind 3
      elseif(jselec(1:5).eq.'jjall') then
       if(iflagr.eq.0) then
        call setjj(nend,nparsv,nholsv,isophs)
        iflagr = nend
        if(imode+nend-1.gt.krpmax) stop ' too many ph-states'
       endif
       jspart = nparsv(iflagr)
       nhole  = nholsv(iflagr)
       isoini = isophs(iflagr)
      elseif(jselec(1:4).eq.'nojj') then
       jspart = 0
      else
       stop ' invalid option jselec'
      endif
c
c     count down flag
c
      iflagr = max(0,iflagr-1)
c      write(14,'(a,2i5)')
c     &   ' QINIT: imode,iflagr=',imode,iflagr
      if(imode+iflagr.gt.krpmax)
     &  stop ' too many modes!'
c
c     setup printing of the mode and print
c
      lentop = 1
      if(iprshl.eq.1) then
       toper(lentop:lentop+4) = 'pr-on'
       lentop = 5
      elseif(iprshl.eq.2) then
       toper(lentop:lentop+5) = 'pr-off'
       lentop = 6
      endif
      if(ihwgt.gt.0) then
       if(ihwgt.gt.9) stop ' ihwgt > 9 is too large'
       toper(lentop:lentop+4) = ' h**'//char(ihwgt+48)
       lentop = lentop+5
      elseif(ihwgt.lt.0) then
c       if(ihwgt.lt.-9 .and. itrpa.gt.10) stop ' ihwgt < -9 is too small'
       toper(lentop:lentop+5) = ' h**-'//char(abs(ihwgt)+48)
       lentop = lentop+6
      endif
      if(ih0wgt.gt.0) then
       if(ih0wgt.gt.9) stop ' ih0wgt > 9 is too large'
       toper(lentop:lentop+5) = ' h0**'//char(ih0wgt+48)
       lentop = lentop+6
      elseif(ih0wgt.lt.0) then
c       if(ih0wgt.lt.-9) stop ' ih0wgt < -9 is too small'
       toper(lentop:lentop+6) = ' h0**-'//char(abs(ih0wgt)+48)
       lentop = lentop+7
      endif
      if(typmod(1:3).eq.'mom'.or.typmod(1:3).eq.'mod'
     &   .or.typmod(1:3).eq.'moj' .or.typmod(1:3).eq.'mot'
     &   .or.typmod(1:3).eq.'msc') then
       write(3,'(f5.2)') xkqrd
       rewind 3
       read(3,'(a5)') afix(1:5)
       rewind 3
       toper(lentop:lentop+10) = ' r**'//char(lrpa+48)//'+'//afix(1:5)
       lentop = lentop+11
      elseif(typmod(1:3).eq.'bes'.or.typmod(1:3).eq.'bed'
     &       .or.typmod(1:3).eq.'bej' ) then
       write(3,'(f5.2)') xkqrd
       rewind 3
       read(3,'(a)') afix(1:5)
       rewind 3
       toper(lentop:lentop+12) =
     &   'j_'//char(lrpa+48)//'('//afix(1:5)//'r/r)'
       lentop = lentop+13
      elseif(typmod(1:3).eq.'loc') then
       write(3,'(f3.0)') xkqrd
       rewind 3
       read(3,'(a3)') afix(1:3)
       rewind 3
       toper(lentop:lentop+10) = ' delta(r-'//afix(1:3)//')'
       lentop = lentop+13
      elseif(typmod(1:2).eq.'ph') then
       write(3,'(i2,a,i2,a,i1)') npart,',',nhole,',',isoini
       rewind 3
       read(3,'(a)') afix(1:7)
       rewind 3
       toper(lentop:lentop+11) = ' ph('//afix(1:7)//')'
       lentop = lentop+12
      endif
      if(jselec(1:2).eq.'jj') then
       write(3,'(i2,a,i2,a,i1)') jspart,',',nhole,',',isoini
       rewind 3
       read(3,'(a)') afix(1:7)
       rewind 3
       toper(lentop:lentop+11) = ' jj('//afix(1:7)//')'
       lentop = lentop+12
      endif
c                                      here comes jselec if ready
      if(ipwgt.eq.-1) then
       toper(lentop:lentop+11) = ' , p=h**-1 q'
       lentop = lentop+11
      elseif(ipwgt.eq.0) then
       toper(lentop:lentop+5) = ' , p=q'
       lentop = lentop+5
      elseif(ipwgt.eq.1) then
       toper(lentop:lentop+7) = ' , p=h q'
       lentop = lentop+7
      endif
c
      write(6,'(/a,i3,a,i1,a,a)')
     & ' mode nr.',imode,': iso=',isorpa,', ',toper(1:lentop)
      write(14,'(/a,i3,a,i1,a,a)')
     & ' mode nr.',imode,': iso=',isorpa,', ',toper(1:lentop)
      CALL FLUSH(14)
c
c     set up the local part of operator 'q' for given 'lrpa' and 'isorpa'.
c     store the radial dependence of 'q' on 'qrad'.
c
      call initqr(typmod(1:3),xkqrd,powqr)
      if(iwr.ge.0 .or. ifanlr.gt.0) then
       write(6,'(/1x,a/1x,a,5g12.4,10(/1x,5g12.4))')
     &  'testprint of the radial operator:',
     &  ' protons:',(qrad(ir),ir=1,ngrid)
       write(6,'(1x,a,5g12.4,10(/1x,5g12.4))')
     &  'neutrons:',(qrad(ir+ngrid),ir=1,ngrid)
      endif
c
c     optionally derivatives for split-momentum initialisation
c
      if(isplit.ne.0) then
       parity = 1-mod(lrpa,2)*2
       do 110 iso=1,2
       isoadd = (iso-1)*ngrid
       do ir=1,ngrid
         drkin(isoadd+ir) = zero
       enddo
c
       if(isplit.eq.1) ilamp   = 1
       if(isplit.eq.2) ilamp   = 2
       lamp   = lrpa-3+ilamp+ilamp
       if(lamp.ge.0) then
        if(lamp.lt.lrpa) then
         facder = sqrt(alrpa/mlrpa)
         facrad = (alrpa+one)*facder
        else
         facder = -sqrt((alrpa+one)/mlrpa)
         facrad = -alrpa*facder
        endif
        call derv17(qrad(isoadd+1),work1(1),parity,ngrid,rstep)
        rkin(isoadd+1,ilamp) = (facder+facrad)*work1(1)
        do ir=2,ngrid
          rkin(isoadd+ir,ilamp) = facder*work1(ir)
     &                           +facrad*qrad(isoadd+ir)*rinvs(ir)
        enddo
        parmas = (mod(lrpa,2)*2-1)
        call extrap(rkin(isoadd+1,ilamp),parmas)
        write(6,'(a,i5,20(/1x,10g12.4))')
     &    ' testprint rkin for ilamp=',ilamp,
     &   (rkin(isoadd+ir,ilamp),ir=1,ngrid)
       else
        do ir=1,ngrid
          rkin(isoadd+ir,ilamp) = zero
        enddo
       endif
c                                              the derivative potential
       if(ilamp.eq.1) then
        facder =  sqrt(alrpa/(alrpa+alrpa+one))
        facrad = -facder*(alrpa-one)
       else
        facder = -sqrt((alrpa+one)/(alrpa+alrpa+one))
        facrad = facder*(alrpa+one+one)
       endif
       call derv17(rkin(isoadd+1,ilamp),work1,-parity,ngrid,rstep)
       do ir=2,ngrid
         rinv   = rinvs(ir)
         drkin(isoadd+ir) = ( facder*work1(ir)
     &                     + facrad*rkin(isoadd+ir,ilamp)*rinv )
     &                     + drkin(isoadd+ir)
       enddo
c                                             zero the rest
       ilampc = 3-ilamp
       do ir=1,ngrid
         rkin(isoadd+ir,ilampc) = zero
         rqcur(isoadd+ir,1)     = zero
       enddo
c
  110  continue
      endif
c
c     optionally prepare potential fields for spin-orbit initialization
c     (potentials are set to 0 except 'rpotls' for later call of
c     'hlinrs').
c
      if(typmod(3:3).eq.'j') then
       do ir=1,ngrid2
         rpot(ir)     = zero
         rbmass(ir)   = zero
         drbmas(ir,1) = zero
         drbmas(ir,2) = zero
       enddo
       call grad(qrad,rpotls,2)
      endif
c
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     now all quantum numbers are provided compute qlrout.
c
      if(typmod(1:2).ne.'ph') then
       if(typmod(3:3).eq.'j') then
        call hrespp(iso,la,xja,nadd,xj,xl,jdif,wf(nadd+1),wfs(nadd+1))
        do ir=1,ngrid
          qlrout(nadd+ir,jdif) = hlinrs(ir)
        enddo
       elseif(jspart.eq.0 .or.
     &    (js.eq.jspart .and. na.eq.nhole .and. iso.eq.isoini) ) then
        redmat = redume(xj,alrpa,xja)
        if(typmod(3:3).eq.'t') then
          adamp = 1.0D0/e0damp
c          adamp = 100.0D0
          write(7,'(a,1pg12.4)')  ' damped Q: e0damp=',e0damp
          do ir=2,ngrid
            qlrout(nadd+ir,jdif) = wf(nadd+ir)*rinvs(ir)
          enddo
          parity = 1-2*mod(la,2)
          call extrap(qlrout(nadd+1,jdif),parity)
          call t0inv(qlrout(nadd+1,jdif),qlrout(nadd+1,jdif),
     &               work1(1),adamp,la)
          do ir=2,ngrid
            qlrout(nadd+ir,jdif) = redmat*qrad(isoadd+ir)
     &                            *qlrout(nadd+ir,jdif)
          enddo
          call t0inv(qlrout(nadd+1,jdif),qlrout(nadd+1,jdif),
     &               work1(1),adamp,ll)
        else
          do ir=2,ngrid
            qlrout(nadd+ir,jdif) = redmat*qrad(isoadd+ir)
     &                            *wf(nadd+ir)*rinvs(ir)
          enddo
        endif
        parity = 1-2*mod(la+lrpa,2)
        call extrap(qlrout(nadd+1,jdif),parity)
        call projec(qlrout(nadd+1,jdif),js,nstack,nmax)
       else
        do ir=1,ngrid
          qlrout(nadd+ir,jdif) = zero
        enddo
       endif
c
c     initialise pure ph-state. states are given via 'xkqrd'.
c
      else
       if(isoini.eq.2) then
        nparts = npart+nmaxp
        nholes = nhole+nmaxp
       else
        nparts = npart
        nholes = nhole
       endif
       if(na.eq.nhole .and. js.eq.jsp(npart) .and. iso.eq.isoini) then
        if(npart.gt.nmax)
     &    stop ' too large particle state called in init. of q'
        if(nhole.gt.nmax)
     &    stop ' too large hole state called in init. of q'
        if(wgt.lt.half) stop ' hole state was not occupied'
        wgtprt = weight(nparts)
        if(wgtprt.gt.half) stop ' particle state was not empty'
        naddpr = (nparts-1)*ngrid
        do ir=2,ngrid
          qlrout(nadd+ir,jdif) = wf(naddpr+ir) * rinvs(ir)
        enddo
        parll  = 1-mod(ll,2)*2
        call extrap(qlrout(nadd+1,jdif),parll)
        write(6,'(4(a,i4))')
     &    ' 1p1h initialised for npart=',npart,' nhole=',nhole,
     &    ' js=',js,'  isoini=',isoini
       else
        do ir=1,ngrid
          qlrout(nadd+ir,jdif) = zero
        enddo
       endif
      endif
c
c     optionally compute momentum for specific 'ilamp' and the
c     according 'qlrout' by linear response.
c
      if(isplit.ne.0) then
       call hrespq(iso,la,xja,nadd,xj,xl,jdif,
     &             wf(nadd+1),wfs(nadd+1),plrout(nadd+1,jdif))
       call projec(plrout(nadd+1,jdif),js,nstack,nmax)
      endif
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
c     restore 'qlrout' according to 'plrout' with specific 'ilamp'
c
      if(isplit.ne.0) call invmrp(ovqpm1,resfac)
c
c     optionally project windows
c
      if(iprshl .eq. 1) then
       call prjshl(1,qlrout)
      elseif(iprshl .eq. 2) then
       call prjshl(0,qlrout)
      endif
c
c     optionally the h0- or hrpa-damping
c
      if(ih0wgt.lt.0) then
c                                             inverse h0**(-n)
       inqpow = abs(ih0wgt)
       if(iwr.ge.-1)
     &  write(6,'(a,i1,a)')
     &   ' initialise "qlrout" with h0**(-',inqpow,')'
       if(inqpow.gt.1) then
        do itini=1,inqpow-1
          call h0inv(qlrout,qlrout,avinv,1)
        enddo
       endif
       call h0inv(qlrout,qlrout,avinv,1)
      elseif(ih0wgt.gt.0) then
c                                             h0**(+n)
       if(ih0wgt.gt.1) then
        do itini=1,ih0wgt-1
          call h0act(qlrout,qlrout,avinv,0)
        enddo
       endif
       call h0act(qlrout,qlrout,avh0,1)
      endif
c
      if(ihwgt.gt.0) then
c                                             h**n q
       do ith=1,ihwgt
         if(mod(ith,2).eq.1)
     &     call commhq(ovqp,ovqp0,ovpp,ovp0p0,1,resfac,0)
         if(mod(ith,2).eq.0)
     &     call commhp(ovpqb,ovpqb0,ovqbqb,ovqbb0,resfac,zero,1)
       enddo
       if(mod(ihwgt,2).eq.1) call copywf(plrout,qlrout)
      elseif(ihwgt.lt.0) then
c                                             inverse h**(-n)
       itinvm = abs(ihwgt)
       if(iwr.ge.-1)
     &  write(6,'(a,i2,a)')
     &   ' initialise "qlrout" with hrpa**(-',itinvm,')'
       do itinv=1,itinvm
         if(mod(itinv,2).eq.1) 
     &      call invmrq(ovqpm1,ovpq,ovpq0,resfac,cmfac)
         if(mod(itinv,2).eq.0) call invmrp(ovqpm1,resfac)
       enddo
       if(mod(itinvm,2).eq.1) call copywf(plrout,qlrout)
      endif
c     optionally project windows
c
      if(iprshl .eq. 1) then
       call prjshl(1,qlrout)
      elseif(iprshl .eq. 2) then
       call prjshl(0,qlrout)
      endif
c
c     normalise 'qlrout'
c
      if(ifprim.ne.1) then
ccc      if(ifprim.eq.-9999) then
       qnorm  = qpover(qlrout,qlrout)
c       write(14,'(a,g12.3)') ' norm in qinit=',qnorm
       qnorm  = one/sqrt(qnorm)
       call multwf(qlrout,qnorm)
c       do 770 ir=1,ngrid
c       qrad(ir+ngrid) = qnorm*qrad(ir+ngrid)
c  770  qrad(ir) = qnorm*qrad(ir)
      endif
c
      return
      end
c-----prjshl------------------------------------------------------
c
      subroutine prjshl(ionshl,alrout)
c
      include 'params.inc'
c
c     projects from 'alrout' all states from the given shell
c     above the fermi-surface.
c     the action of the projection is determined by
c      ionshl = 1 --> projects 'qlrout' "on" the shell above,
c      ionshl = 0 --> the shell above is projected away from 'qlrout'
c
      include 'wavfun.inc'
      include 'linres.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'constn.inc'
c
      dimension alrout(KX11,KLRPA)
c      dimension nprinc(56),nmin(2)
c
c     'nprinc' translates book-keeping label 'na' or 'nb' into
c     the principle quantum number of the harmonic oscillator.
c
c      data nprinc/1,2,2,3,3,3,4,4,4,4,5,5,5,5,5,6,6,6,6,6,6,
c     &            7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,
c     &            10,10,10,10,10,10,10,10,10,10,11/
c     data nminfl/1,2, 1,6, 1,8, 2,20, 3,28,  5,50,  8,82, 12,126,
c    &            17,184, 21,258, 30,350/
c
c-------------------------------------------------------------------------
c
c     compute projection shell below fermi niveau
c
c      if(nmaxp.gt.56 .or. nmaxn.gt.56)
c     &    stop ' too large nr. of shells in "prjshl"'
c
      include 'initdo.inc'
c
      do ir=1,ngrid
        work1(ir) = zero
      enddo
c
      do nb=1,nmax
        jsb   = jsp(nb)
        lb    = jsb/2
        wgtb  = weight(nb+nstack/ngrid)
cold        if(jsb.eq.js .and. wgtb.lt.half
cold       &   .and. abs(nprinc(na)-nprinc(nb)).le.0) then
        if(jsb.eq.js .and. wgtb.lt.half) then
          write(06,'(a,5i5)')
     &      ' project: iso,na,nb,jsa,jsb=',iso,na,nb,jsa,jsb
          xlb   = lb
          xjb   = (jsb-1)/2+0.5D0
          naddb = nstack + (nb-1)*ngrid
c
c         overlap
c 
          sum   = zero
          do ir=2,ngrid
            sum = rfact(ir)*alrout(nadd+ir,jdif)*wf(naddb+ir) + sum
          enddo
          sum   = rstep*sum
c
c         accumulate parts of 'qlrout' on 'wf'-shell
c
          do ir=1,ngrid
            work1(ir) = sum*wf(naddb+ir)*rinvs(ir) + work1(ir)
          enddo
c
       endif
      enddo
c
c     prepare 'qlrout' according to 'ionshl'
c
      if(ionshl.eq.1) then
        do ir=1,ngrid
          alrout(nadd+ir,jdif) = work1(ir)
        enddo
      else
        do ir=1,ngrid
          alrout(nadd+ir,jdif) = alrout(nadd+ir,jdif) - work1(ir)
        enddo
      endif
c
      include 'enddo.inc'
c
      return
      end
c-----initqr-----------------------------------------------------------------
c
      subroutine initqr(typloc,xkqrd,powqr)
c
      include 'params.inc'
c
c     initialises the local collective operator 'qrad' which is the
c     local part of the operator 'q' for given 'lrpa' and 'isorpa'.
c
      character*3 typloc
      dimension beszer(0:11)
c      double precision s,c,x,x2,x4,bessel
c
      include 'wavfun.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'grid.inc'
      include 'constn.inc'
c
      data beszer/3.14,4.49,5.76,6.99,8.19,9.36,10.51,11.66,12.79,
     &            13.92,15.03,16.14/
c
c--------------------------------------------------------------------
c
c     set some scaling factors depending on angular momentum
c     and isospin
c
      if(lrpa.eq.1) then
       rscale = 2.0/isorpa
      else
       rscale = 2.0*anumb*rms*rms/isorpa
      endif
      if(lrpa.eq.0) then
       powqr  = 2.0
       opnorm = one
      else
       powqr  = alrpa
       opnorm = one
      endif
      if(isorpa.eq.2) then
       facpr  = aneut/anumb
       facne  = -aprot/anumb
      elseif(isorpa.eq.3) then
       facpr  = 1.0D0
       facne  = 0.0D0
      else
       facpr  = 1.0D0
       facne  = 1.0D0
      endif
      if(typloc.eq.'bes'.or.typloc.eq.'bed'.or.typloc.eq.'bej') then
c                     spherical bessel of order lrpa
       if(lrpa.gt.11) stop ' lrpa>11 not implemented in "initqr"'
       qrad(1+ngrid) = zero
       qrad(1) = zero
       do i=2,ngrid
         r      = (i-1)*rstep
         x      = r*xkqrd*beszer(lrpa)/rdms
         bessl  = bessel(x,lrpa)
         qrad(i+ngrid) = facne*bessl
         qrad(i) = facpr*bessl
       enddo
      elseif(typloc.eq.'mom'.or.typloc.eq.'mod'
     &       .or.typloc.eq.'moj' .or.typloc.eq.'mot') then
c                     power initialization
       do ir=2,ngrid
         r      = (ir-1)*rstep
         rform  = r**powqr*(r/rdms)**xkqrd
         qrad(ir) = opnorm*facpr*rform
         qrad(ngrid+ir) = opnorm*facne*rform
       enddo
      elseif(typloc.eq.'msc') then
c                     power initialization with c.m. correction (L=0)
       rfac2  = 1.6666666666667D0*rms*rms  * KCMCORR
       do ir=2,ngrid
         r      = (ir-1)*rstep
         rform  = r**3-rfac2*r
         qrad(ir) = opnorm*facpr*rform
         qrad(ngrid+ir) = opnorm*facne*rform
       enddo
      elseif(typloc(1:2).eq.'ph') then
c                     zero 'qrad' in case of ph-excitation
       do ir=1,ngrid
         qrad(ir) = zero
         qrad(ngrid+ir) = zero
       enddo
      elseif(typloc.eq.'loc') then
c                     local initialization
       iloc  = xkqrd
       do ir=2,ngrid
         if(ir.ne.iloc) then
          qrad(ir) = zero
          qrad(ngrid+ir) = zero
         else
          qrad(ir) = opnorm*facpr
          qrad(ngrid+ir) = opnorm*facne
         endif
       enddo
      else
       stop ' inappropriate option typmod in initqr!'
      endif
      if(typloc.eq.'mod' .or. typloc.eq.'bed') then
c       write(6,'(a)') ' testprint saxon-woods:'
       do ir=1,ngrid
         r      = (ir-1)*rstep
         qrad(ir) = one/(one+exp((r-rdmsp)/(0.5*surfp)))*qrad(ir)
         if(nmaxn.gt.0) then
          irn = ir+ngrid
          qrad(irn) = one/(one+exp((r-rdmsn)/(0.5*surfn)))*qrad(irn)
         endif
       enddo
      endif
ctest
      if(ifanlr.gt.0) then
       write(6,'(a,20(/1x,5g12.4))')
     &  ' testprint qrad:',(qrad(i),i=1,ngrid)
       write(6,'(1x,5g12.4)')
     &  (qrad(i),i=ngridp,ngrid2)
      endif
ctest
      qrad(1) = zero
      qrad(ngridp) = zero
c
      return
      end
c-----setph -------------------------------------------------------------------
c
      subroutine setph(nend,ecutlw,ecutup,nparsv,nholsv,isophs)
c
      include 'params.inc'
c
c     determines the ph states in the valence space
c
      dimension nparsv(krpmax),nholsv(krpmax),isophs(krpmax)
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'result.inc'
      include 'skypar.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      data xminrd/1.0D-5/
c
c--------------------------------------------------------------------
c     set counter
c
      write(14,'(2a,i5)')  ' single particle spectra from setph: ',
     &         '  krpmax=',krpmax
      nme    = 0
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     now all a.m. quantum numbers are provided, run the loop
c     over 'nb'
c
      if(wgt.le.half) goto 30
      redmat = redume(xj,alrpa,xja)
      if(abs(redmat).le.xminrd) goto 30
c
c     loop over the b-states
c
      do 40 nb=1,nmax
      jsb   = jsp(nb)
      if(jsb.ne.js) goto 40
      lb    = jsb/2
      xlb   = lb
      xjb   = (jsb-1)/2+0.5D0
      wgtb  = weight(nb+nstack/ngrid)
      if(wgtb.gt.half) goto 40
      naddb = nstack + (nb-1)*ngrid
      subfcb = spenrg(nb+nstack/ngrid)
c
c     assign quantum numbers in pointer fields 'iniqr', 'iniiso' and
c     'xkqr'.
c
      npart=nb
      nhole=na
      edif = abs(subfac-subfcb)
      if(edif.lt.ecutup .and. edif.gt.ecutlw
     &   .and. subfac.ne.zero .and. subfcb.ne.zero) then
       nme=1+nme
       if(nme.le.krpmax) then
        nparsv(nme) = npart
        nholsv(nme) = nhole
        isophs(nme) = iso
        iwarn = 0
        nend  = nme
       else
        iwarn = 1
       endif
       if(iwr.ge.-1)
     &  write(6,'(i4,2x,2i3,4(1pe13.4))') nme,npart,nhole,edif
       if(iwr.ge.-1)
     &  write(14,'(i4,2x,2i3,(1pe13.4),i3,2(i4,0pf6.1))') 
     &    nme,npart,nhole,edif,iso,la,xja,lb,xjb
      endif
c
   40 continue
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
c
      if(iwarn.eq.1) stop ' too many ph modes in setph'
c      nend  = nme
c
c     print modes
c
c      write(6,'(a)') ' the explicit ph-modes from setph:'
c      write(6,'(a,20(/1x,10(i2,1x,i2,5x)))') ' npart,nhole:',
c     &    (nparsv(imode),nholsv(imode),imode=1,nend)
c      write(14,'(a)') ' the explicit ph-modes from setph:'
c      write(14,'(a,20(/1x,5(i2,1x,i2,5x)))') ' npart,nhole:',
c     &    (nparsv(imode),nholsv(imode),imode=1,nend)
c
      return
      end
c-----setjj -------------------------------------------------------------------
c
      subroutine setjj(nend,jssav,nholsv,isophs)
c
      include 'params.inc'
c
c     determines the (js,na,iso) combinations in the valence space
c
      dimension jssav(krpmax),nholsv(krpmax),isophs(krpmax)
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'result.inc'
      include 'skypar.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
c
c     shells and js-transitions tabulated
c
      dimension numsav(7),nminsv(7),nmaxsv(7),nextsv(7),
     &          jsmsav(7),jsblsv(7)
      data numsav/  6, 8, 20, 28, 50, 82,126/
      data nminsv/  1, 2,  4,  4,  7, 10, 15/
      data nmaxsv/  2, 3,  6,  7, 11, 16, 22/
      data nextsv/  2, 0,  0,  7, 11, 16, 22/
      data jsmsav/  5, 7,  9,  9, 11, 13, 15/
      data jsblsv/  3, 0,  0,  7,  9, 11, 13/
c
      data xminrd/1.0D-5/
c
c--------------------------------------------------------------------
c     set counter
c
      nme    = 0
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     now all a.m. quantum numbers are provided, run the loop
c     over 'nb'
c
c
c     set actual triggers
c
      if(iso.eq.1) then
       npnumb = nprot
      else
       npnumb = nneut
      endif
      do icase=1,7
        if(numsav(icase).eq.npnumb) then
         nminac = nminsv(icase)
         nmaxac = nmaxsv(icase)
         nextra = nextsv(icase)
         jsmaxa = jsmsav(icase)
         jsbloc = jsblsv(icase)
         goto 219
        endif
      enddo
      stop ' particle number not appropriate for rpa program'
  219 continue
c
      if(na.lt.nminac .or. na.gt.nmaxac) goto 30
      if(na.eq.nextra) then
       if(js.gt.jsmaxa+2) goto 30
      else
       if(js.gt.jsmaxa .or. js.eq.jsbloc) goto 30
      endif
      redmat = redume(xj,alrpa,xja)
      if(abs(redmat).le.xminrd) goto 30
c
c     assign quantum numbers in pointer fields
c
      nhole=na
      nme=1+nme
      jssav(nme)  = js
      nholsv(nme) = nhole
      isophs(nme) = iso
      if(iwr.ge.0)
     &  write(6,'(3i4)') js,nhole,iso
      if(iwr.ge.1)
     &  write(14,'(3i4)') js,nhole,iso
c
   40 continue
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
      nend  = nme
c
c     print modes
c
      write(6,'(a)') ' the explicit ph-modes from setjj:'
      write(6,'(a,20(/1x,10(3i3,5x)))') ' js,nhole,iso:',
     &    (jssav(imode),nholsv(imode),isophs(imode),imode=1,nend)
      write(14,'(a)') ' the explicit ph-modes from setjj:'
      write(14,'(a,20(/1x,5(3i3,5x)))') ' js,nhole,iso:',
     &    (jssav(imode),nholsv(imode),isophs(imode),imode=1,nend)
c
c
      return
      end
c-----moment-------------------------------------------------------------------
c
      subroutine moment(ovqq,ovqhq,ovqp,ovpp,ovpqb,ovqh0q,qqinv,widrpa,
     &                  qrlpr,qrlne,resfac,cmfac,
     &                  eshifq,eshift,eshifi,
     &                  nameq,namep,namehp,
     &                  ipinit,ifanly)
c
      include 'params.inc'
c
c     for given q-operator, the [h,q] and [h,[h,q]] is computed.
c     q is input on 'qlrout',
c     the '[h,q]' is output on pointer '[h,q] ' on file for018.
c     p is computed according to 'ipinit':
c       ipinit = 1  -->  p = [h,q]
c       ipinit = 0  -->  p ~ q
c       ipinit = -1 -->  p = h**(-1)[q]
c     it is output on 'plrout' and on pointer 'namep' on file for018.
c     the '[h,p]' is output on pointer '[h,p] ' on file for018.
c     the fields 'qlrout' and 'plrout' are transferred via
c     common /linres/.
c
c     there are furthermore some prints of the results, and
c     some analysing steps, switched by 'ifanly'.
c
c     the list parameters are:
c      ovqq    = overlap <q*q>
c      ovqp    = overlap <q [h,q]>
c      ovpqb   = overlap <[q,h] [h,[h,q]]>
c      nameq   = key_name of q operator for storage in 'buffql'
c      namep   = key_name of p operator for storage in 'buffql'
c      namehp  = key_name of [h,p] operator for storage in 'buffql'
c      ipinit  = switch for computation of 'p' (see above).
c      ifanly  = switch for final analysis of gsc and tda
c
      character*(7) nameq,namep,namehp
c
      include 'skypar.inc'
      include 'wavfun.inc'
      include 'linres.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'result.inc'
      include 'constn.inc'
c
c----------------------------------------------------------------------------
c
c      check input
c
       ovqq    = qnorm(qlrout)
       if(ovqq.eq.zero) then
        write(6,'(a)') ' field qlrout is empty. --> quit.'
        stop
       endif
       call h0inv(qlrout,savlr,qqinv,0)
c
c      compute [h,q]
c
       call commhq(ovqhq,ovqh0q,ovhqhq,ovh0qq,1,resfac,1)
       call buffql(plrout,'o','[h,q] ')
       call buffql(savlr,'o','[h0,q]')
       call h0inv(plrout,savlr,ppinv,0)
c
c      compute p
c
       if(ipinit.eq.-1) then
        if(ifsepi.eq.2) stop ' IFSEPI=2 incompatible with IPINIT=-1'
        call invmrq(ovqpm1,ovpqbm,ovpq0m,resfac,cmfac)
        call buffql(qlrout,'i',nameq)
       elseif(ipinit.eq.0) then
        if(ifsepi.eq.2) stop ' IFSEPI=2 incompatible with IPINIT=0'
        call copywf(qlrout,plrout)
       elseif(ipinit.le.-999) then
        call buffql(plrout,'i',namep)
       endif
       ovqp   = qpover(qlrout,plrout)
       ovpp   = qpover(plrout,plrout)
       if(ipinit.gt.-999) then
        call buffql(plrout,'o',namep)
       else
        call buffql(qlrout,'i','[h,q] ')
        ovphq  = qpover(qlrout,plrout)
       endif
c
c      in case of VPM compute separable X operator
c
       if(ifsepi.eq.2) then
         call commvp(ovpqb,ovpqb0,ovqbqb,ovqbb0,resfac,cmfac,1)
         call buffql(qlrout,'o',namehp)
         write(14,'(a,2g13.5)')
     &    ' <[P,[V,P]]>,<[V,P]**2>=',
     &    qpover(qlrout,plrout),qpover(qlrout,qlrout)
         return
       endif
c
c      compute [h,p]
c
c      call copywf(qlrout,savlr)
       call multwf(savlr,zero)
cc       call commhp(ovpqb,ovpqb0,ovqbqb,ovqbb0,resfac,zero,1)
       call commhp(ovpqb,ovpqb0,ovqbqb,ovqbb0,resfac,cmfac,1)
       call buffql(qlrout,'o','[h,p] ')
       call buffql(savlr,'o','[h0,p]')
       call buffql(savlr,'i',nameq)
       ovqqb  = qpover(savlr,qlrout)
c
c      now compute some energy-averages and rpa-properties
c
       qwid   = ovqq
       redund = ovqhq*ovqhq/(ovqq*ovhqhq)
       redun0 = ovqh0q*ovqh0q/(ovqq*ovh0qq)
       enqhq  = ovqhq/ovqq
       enqhq0 = ovqh0q/ovqq
       enqh2q = sqrt(ovhqhq/ovqq)
       enqh20 = sqrt(ovh0qq/ovqq)
       enwidq = sqrt(abs(ovhqhq/ovqq-enqhq*enqhq))
       enwd0q = sqrt(abs(ovh0qq/ovqq-enqhq0*enqhq0))
c
       enrpa  = sign(sqrt(abs(ovpqb*ovqhq/(ovqp*ovqp))),ovpqb)
       widrpa = enrpa/ovqhq*half
       enrpa0 = sign(sqrt(abs(ovpqb0*ovqhq/(ovqp*ovqp))),ovpqb0/ovqp)
       widrp0 = enrpa0/ovqhq*half
       resinq = mlrpa*half*(ovpqb-ovpqb0)/(ovqp*ovqp)
       eshifq = resinq*ovqq/mlrpa
       resinp = mlrpa*half*(ovqhq-ovqh0q)
       eshifp = resinp*ovpp/(mlrpa*ovqp*ovqp)
c
       enph2p = sqrt(ovqbqb/ovpp)
       enph20 = sqrt(ovqbb0/ovpp)
       enphp  = ovpqb/ovpp
       enphp0 = ovpqb0/ovpp
       enwidp = sqrt(abs(enph2p*enph2p-enphp*enphp))
       enwd0p = sqrt(abs(enph20*enph20-enphp0*enphp0))
c
       enqhqi = ovqq/max(xsmall,qqinv)
       enphpi = ovpp/max(xsmall,ppinv)
       enqhqa = sqrt(enqhq0*enqhqi)
       enphpa = sqrt(enphp0*enphpi)
c
       if(ifanly.ge.1)
     &  write(6,'(a,20(/5x,a,g12.5,a,g12.5))')
     &    '   raw overlaps:',
     &    '             (q**2)=',  ovqq,
     &    '           (q[h,q])=',  ovqhq,
     &    '       ([q,h][h,q])=',  ovhqhq,
     &    '           (p[h,p])=', ovpqb,
     &    '          (q[h0,q])=',ovqh0q,
     &    '     ([q,h0][h0,q])=',ovh0qq,
     &    '          (p[h0,p])=',ovpqb0
       write(6,'(a,4(/5x,2(a,g12.5)))')
     &    '   rpa properties:',
     &    '     energy=', enrpa,'     width=',widrpa,
     &    '  h0-energy=',enrpa0,'  h0-width=',widrp0,
     &    '    kappa-q=',resinq,'   e-shift=',eshifq,
     &    '    kappa-p=',resinp,'   e-shift=',eshifp
       write(6,'(a,g12.5,a,g12.5/a,2(/5x,3(a,g12.5)))')
     &    '   redundance:       q--p =',redund,
     &    '      q--p0 =',redun0,
     &    '   q-averaged energies:',
     &    '    e1--q =', enqhq,'    e2--q =',enqh2q,
     &    '  delta-e =',enwidq,
     &    '   e01--q =',enqhq0,'   e02--q =',enqh20,
     &    ' delta-e0 =',enwd0q
       write(6,'(a,2(/5x,3(a,g12.5)))')
     &    '   p-averaged energies:',
     &    '    e1--p =', enphp,'    e2--p =',enph2p,
     &    '  delta-e =',enwidp,
     &    '   e01--p =',enphp0,'   e02--p =',enph20,
     &    ' delta-e0 =',enwd0p
c
      if(ifanly.ge.1) then
        write(14,'(a,2(/5x,3(a,g12.5)))')
     &    '   q-averaged energies:',
     &    '    e1--q =', enqhq,'    e2--q =',enqh2q,
     &    '  delta-e =',enwidq,
     &    '   e01--q =',enqhq0,'   e02--q =',enqh20,
     &    ' delta-e0 =',enwd0q
        write(14,'(a,2(/5x,3(a,g12.5)))')
     &    '   p-averaged energies:',
     &    '    e1--p =', enphp,'    e2--p =',enph2p,
     &    '  delta-e =',enwidp,
     &    '   e01--p =',enphp0,'   e02--p =',enph20,
     &    ' delta-e0 =',enwd0p
      endif
      if(ifanly.ge.0) then
       write(14,'(a,4(/5x,2(a,g12.5)))')
     &    '   rpa properties:',
     &    '     energy=', enrpa,'     width=',widrpa,
     &    '  h0-energy=',enrpa0,'  h0-width=',widrp0,
     &    '    kappa-q=',resinq,'   e-shift=',eshifq,
     &    '    kappa-p=',resinp,'   e-shift=',eshifp
      else
       write(14,'(4(a,f10.3))')
     &    '  e01-q=',enqhq0,'  h0-en.=',enrpa0,'  rpa-en.=',enrpa,
     &    '  shift=',eshifq
      endif
c
       write(6,'(a,2(/5x,2(a,g12.5)))')
     &    '   inverse h0 energies:',
     &    '  e0-1--q =',enqhqi,'  average =',enqhqa,
     &    '  e0-1--p =',enphpi,'  average =',enphpa
       if(ipinit.lt.-1000) then
        dele2p = half*enrpa*(ovhqhq-2.0*ovphq+ovpp)
        enrpa2 = enrpa*enrpa
        dele2q = half/enrpa*(ovqbqb-2.0*enrpa2*ovqqb+enrpa2**2*ovqq)
        dele   = sqrt(abs(dele2p+dele2q))
        write(6,'(a,3(/3(a,g11.4)))')
     &   ' h**0 weighted error:',
     &   '  qhhq=',ovhqhq,' phq=',ovphq,' pp=',ovpp,
     &   '  phhp=',ovqbqb,' qhp=',ovqqb,' qq=',ovqq,
     &   '  dele2-p=',dele2p,' dele2-q=',dele2q,'  dele=',dele
        write(14,'(a,3(/3(a,g11.4)))')
     &   ' h**0 weighted error:',
     &   '  qhhq=',ovhqhq,' phq=',ovphq,' pp=',ovpp,
     &   '  phhp=',ovqbqb,' qhp=',ovqqb,' qq=',ovqq,
     &   '  dele2-p=',dele2p,' dele2-q=',dele2q,'  dele=',dele
       endif
c
c      compute the energy-moments
c
      if(ifanly.ge.2) then
       emohf0 = ovqq
       emohf1 = ovqp0
       emohf2 = ovp0p0
       emohf3 = ovpqb0
       emorp1 = emohf1
       emorp3 = emorp1*enrpa*enrpa
       reseff = mlrpa*half*(enrpa*enrpa-enqhq*enqhq)/(enqhq*ovqq)
       write(6,'(a,3(/5x,2(a,g12.5))/5x,a,g12.5)')
     &     '   energy moments:',
     &     '   hf mom. 0 =',emohf0,'   hf mom. 1 =',emohf1,
     &     '   hf mom. 2 =',emohf2,'   hf mom. 3 =',emohf3,
     &     '  rpa mom. 1 =',emorp1,'  rpa mom. 3 =',emorp3,
     &     '  eff. kappa =',reseff
      endif
c
c     transition density and fraction of sumrule
c
      fac    = one/sqrt(abs(2.0*enrpa))
      if(isorpa.eq.1) then
       facpr  = fac
       facne  = fac
      else
       facpr  = fac*aneut/anumb
       facne  = -fac*aprot/anumb
      endif
      if(ifgrho.ne.0 .and. ifanly.ge.1) then
       write(6,'(/a,20(/1x,10g12.4))')
     &     ' transition density for protons:',
     &     (prho(ir),ir=1,ngrid)
       write(6,'(/a,20(/1x,10g12.4))')
     &     ' transition density for neutrons:',
     &     (prho(ir),ir=ngridp,ngrid2)
      endif
      if(lrpa.eq.0) then
       apow   = alrpa+2.0
      else
       apow   = alrpa
      endif
c      if(KSCDIP.eq.1 .and. lrpa.eq.1) then
      if(KSCDIP.eq.1) then
c        rcorfc = rms*rms*5.0/3.0
        write(14,'(a,f8.4)') ' now running QRADOB'
        sump   = zero
        sumn   = zero
        do ir=2,ngrid
c          rfac   = ((ir-1)*rstep)
c          rfac   = rfac*(rfac*rfac-rcorfc)
          rfac   = qradob(ir)*((ir-1)*rstep)**2
          sump   = rfac*prho(ir) + sump
          sumn   = rfac*prho(ir+ngrid) + sumn
        enddo
      else
        sump   = zero
        sumn   = zero
        do ir=2,ngrid
          rfac   = ((ir-1)*rstep)**(apow+2.0)
          sump   = rfac*prho(ir) + sump
          sumn   = rfac*prho(ir+ngrid) + sumn
        enddo
      endif
      qrlpr  = frpirs*sump
      qrlne  = frpirs*sumn
c                                         correct for double countings
      qrlpr  = qrlpr*D4PI*sqrt(2.0D0*alrpa+1.0D0)
      qrlne  = qrlne*D4PI*sqrt(2.0D0*alrpa+1.0D0)
c
      write(14,'(a,2g12.4)') ' qrlpr,qrlne=', qrlpr,qrlne
      qrlto  = qrlpr+qrlne
      qrlvc  = facpr*qrlpr+facne*qrlne
      widfac = 0.125/(widrpa*ovqp*ovqp)
      belp   = qrlpr*qrlpr*widfac
      belt   = qrlto*qrlto*widfac
      belv   = qrlvc*qrlvc*widfac
      ewsrp  = belp*enrpa
      ewsrt  = belt*enrpa
      ewsrv  = belv*enrpa
      write(6,'(3(a,g12.4))')
     &     ' r**lrpa  total:  <q**2>=',belt,
     &     '  ewsr=',ewsrt
      write(6,'(3(a,g12.4))')
     &     ' r**lrpa charge:  <q**2>=',belp,
     &     '  ewsr=',ewsrp
      write(6,'(3(a,g12.4))')
     &     ' r**lrpa vector:  <q**2>=',belv,
     &     '  ewsr=',ewsrv,' <r**l p>=',qrlvc
      if(ifanly.ge.0) then
       write(14,'(3(a,g12.4))')
     &     ' r**lrpa  total:  <q**2>=',belt,
     &     '  ewsr=',ewsrt
       write(14,'(3(a,g12.4))')
     &     ' r**lrpa charge:  <q**2>=',belp,
     &     '  ewsr=',ewsrp
      else
c       write(14,'(2(a,g11.3))')
c     &   ' delta(r**lrpa): scalar=',belt,'   proton=',belp
      endif
c
c     print b(el)-values for 1ph-transitions
c
      if(ipinit.eq.0) then
       qrlt0  = qrlpr+qrlne
       qrlt0  = qrlt0*qrlt0
       qrlt1  = (aneut*qrlpr-aprot*qrlne)/anumb
       qrlt1  = qrlt1*qrlt1
       write(14,'(a,4g12.5)')
     &   ' e,b(el)"s: ',enrpa0,qrlpr*qrlpr,qrlt0,qrlt1
      endif
c
      if(ifanly.le.2) return
c
c      some analysis of the single-mode response
c
c
      if(mlrpa.ne.3 .or. isorpa.ne.1) then
c
c     correction to energy, using direct energy moments
c
       ecor   = half*(enrpa
     &              -enrpa*(widrpa*qwid+pwid/widrpa))
       resinq = (ovqp-ovqp0)/ovqq
       resinp = (ovpqb-ovpqb0)/ovpp
       epertq = -0.25*resinq*resinq/enqhq0
       epertp = -0.25*resinp*resinp/enphp0
       epert  = half*(epertq+epertp)
       ecord  = ecor-epert
       write(6,'(1x,2(a,f9.4)/20x,3(a,f9.4))')
     &     'correlation energies: ecor=',ecor,
     &     '  ecor-reduced=',ecord,
     &     ' pert.=',epert,' pert-q=',epertq,' pert-p=',epertp
       write(14,'(1x,2(a,f9.4)/20x,3(a,f9.4))')
     &     'correlation energies: ecor=',ecor,
     &     '  ecor-reduced=',ecord,
     &     ' pert.=',epert,' pert-q=',epertq,' pert-p=',epertp
      endif
      ecored = ecord
c
c     inverse energy rpa
c
      if(ifirpa.eq.1) then
       call buffql(qlrout,'i',nameq)
       call h0inv(qlrout,savlr,qqinv,0)
       call invmrq(ovqpm1,ovpqb,ovpqb0,resfac,cmfac)
c
c      final analysis
c
       ovppm1 = qnorm(plrout)
       pwidm1 = ovppm1/(4.0*ovqpm1*ovqpm1)
       erpam1 = sqrt(ovqp/ovqpm1)
       qwidm1 = sqrt(ovqp*ovqpm1)
       resinv = mlrpa*half*(ovpqb-ovpqb0)/(ovqpm1*ovqpm1)
       eshifi = resinv*ovqq/mlrpa
       write(6,'(/1x,4(a,g11.3))')
     &     ' invmrq: erpa=',erpam1,' qwidm1=',qwidm1,
     &     ' pwidm1=',pwidm1,' kappa=',resinv
       write(14,'(/1x,a,f9.3,2(a,0pg11.3)/8x,2(a,g12.4))')
     &     ' invmom: erpa=',erpam1,
     &    ' pwid=',pwidm1,' qwid=',qwidm1,
     &    ' inv-kappa=',resinv,'  eshift=',eshifi
c
       wdrp0m = sqrt(abs(ovpqb0/ovqpm1))/ovqpm1
       enphp0 = ovpqb0/ovppm1
c
       call h0inv(plrout,savlr,ppinv,0)
       enqhqi = ovqq/max(xsmall,qqinv)
       enphpi = ovpp/max(xsmall,ppinv)
       enqhqa = sqrt(enqhq0*enqhqi)
       enphpa = sqrt(enphp0*enphpi)
c
       widrpm = half/qwidm1
       ecor   = half*(erpam1
     &              -erpam1*(widrpm*qwid+pwidm1/widrpm))
       epertq = -0.25*resinq*resinq/enqhqa
       epertp = -0.25*resinp*resinp/enphpa
       epert  = half*(epertq+epertp)
       ecord  = ecor-epert
       write(6,'(1x,2(a,f9.4)/20x,3(a,f9.4))')
     &     'inv. corr.  energies: ecor=',ecor,
     &     '  ecor-reduced=',ecord,
     &     '  pert.=',epert,' pert-q=',epertq,' pert-p=',epertp
       write(14,'(1x,2(a,f9.4)/20x,3(a,f9.4))')
     &     'inv. corr.  energies: ecor=',ecor,
     &     '  ecor-reduced=',ecord,
     &     '  pert.=',epert,' pert-q=',epertq,' pert-p=',epertp
      else
       eshifi = zero
      endif
c
      if(mlrpa.ne.3 .or. isorpa.ne.1) then
c
c      test tda
c
       call buffql(plrout,'i',nameq)
       call commhp(ovpqb,ovpqb0,ovqbqb,ovqbb0,one,zero,0)
       pwid   = qnorm(plrout)
       enrtd  = ovpqb/pwid
       en0td  = ovpqb0/pwid
       resint = mlrpa*half*(enrtd-en0td)/pwid
       eshift = resint*ovqq/mlrpa
       enrtd2 = sqrt(ovqbqb/pwid)
       en0td2 = sqrt(ovqbb0/pwid)
       widtda = sqrt(abs(enrtd2*enrtd2-enrtd*enrtd))
       widtd0 = sqrt(abs(en0td2*en0td2-en0td*en0td))
       erpatd = sqrt(abs(en0td*en0td+en0td*resint*pwid*2.0/mlrpa))
       write(6,'(a,3(/5x,3(a,g12.4)))')
     &     '  the tamm-dancoff aproach:',
     &     ' tda-energy=',enrtd,'   2.moment=',enrtd2,
     &     '   width=',widtda,
     &     ' tda-energ0=',en0td,'   2.moment=',en0td2,
     &     '   width=',widtd0,
     &     ' eff. kappa=',resint,' rpa-energy=',erpatd
       write(14,'(a,3(/5x,3(a,g12.4)))')
     &     '  the tamm-dancoff aproach:',
     &     ' tda-energy=',enrtd,'   2.moment=',enrtd2,
     &     '   width=',widtda,
     &     ' tda-energ0=',en0td,'   2.moment=',en0td2,
     &     '   width=',widtd0,
     &     ' tda-kappa=',resint,' eshift=',eshift,
     &     ' rpa-energy=',erpatd
      else
c
c      test c.m.motion
c
       pwid   = ovpp
       call cminit
       call commhp(ovpqb,ovpqb0,ovqbqb,ovqbb0,one,cmfac,0)
       exwid  = ovpp
       exwid  = max(exwid,xsmall)
       cosang = ovqp*ovqp/(pwid*exwid)
       residu = ovqbqb/exwid
       resid0 = ovqbb0/exwid
       enrtd  = ovpqb
       en0td  = ovpqb0
       write(6,'(a/4(a,g11.3)/4(a,g11.3))')
     &     ' third the c.m.test:',
     &     ' p-width=',  pwid,' ex.wid=', exwid,' cos.angle=',cosang,
     &     '    ovqp=',ovqp,
     &     ' res.no.=',residu,' res.n0=',resid0,' res.ener.=',enrtd,
     &     ' res.en0=',en0td
c
      endif
c
      return
      end
c-----prirpa-----------------------------------------------------
c
      subroutine prirpa(beta,qtrapr,qtrane,omega,maxmod,maxout)
c
      include 'params.inc'
      parameter(krpmat=krpmax*(krpmax+1)/2)
c
c     computes b(el)-values and prints spectrum.
c
      include 'wavfun.inc'
      include 'result.inc'
      include 'constn.inc'
      include 'linres.inc'
      common /bels/belsavpr(krpmax),belsavne(krpmax),
     &             belsavsc(krpmax),belsavvc(krpmax)
c
      dimension omega(krpmax),
     &          beta(krpmax,krpmax),
     &          qtrapr(krpmax),qtrane(krpmax)
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
c     compute b(el)-values and energy-weighted sumrule
c
      write(6,'(22x,a/a,a)')
     &   'proton               isoscalar             isovector',
     &   ' nr  energy      b(el)     e*b(el)     b(el)     e*b(el)  ',
     &   '   b(el)     e*b(el)'
      write(14,'(22x,a/a,a)')
     &   'proton               isoscalar             isovector',
     &   ' nr  energy      b(el)     e*b(el)     b(el)     e*b(el)  ',
     &   '   b(el)     e*b(el)'
      e3srpr = zero
      e3srsc = zero
      e3srvc = zero
      e2srpr = zero
      e2srsc = zero
      e2srvc = zero
      ewsrpr = zero
      ewsrsc = zero
      ewsrvc = zero
      einrpr = zero
      einrsc = zero
      einrvc = zero
      sumblp = zero
      sumbls = zero
      sumblv = zero
      do imode=1,maxout
        energy = sqrt(abs(omega(imode)))
        omega(imode) = energy
c
        facpr  = aneut/anumb
        facne  = -aprot/anumb
        sumpr  = zero
        sumne  = zero
        sumsc  = zero
        sumvc  = zero
        do imod2=1,maxmod
          sumpr  = qtrapr(imod2)*beta(imod2,imode) + sumpr
          sumne  = qtrane(imod2)*beta(imod2,imode) + sumne
          sumsc  = (qtrapr(imod2)+qtrane(imod2))*beta(imod2,imode) 
     &            + sumsc
          sumvc  = (facpr*qtrapr(imod2)+facne*qtrane(imod2))
     &        *beta(imod2,imode) + sumvc
        enddo
c
        bl2pr  = sumpr*sumpr/(4.0*omega(imode))  
        bl2ne  = sumne*sumne/(4.0*omega(imode))  
        bl2sc  = sumsc*sumsc/(4.0*omega(imode))
        bl2vc  = sumvc*sumvc/(4.0*omega(imode))
!        IF(lrpa==2) THEN
!          bl2pr  = bl2pr*0.6666667D0
!          bl2ne  = bl2ne*0.6666667D0
!          bl2sc  = bl2sc*0.6666667D0
!          bl2vc  = bl2vc*0.6666667D0
!        END IF
        ewpr   = bl2pr*energy
        ewsc   = bl2sc*energy
        ewvc   = bl2vc*energy
        ewsrpr = ewpr+ewsrpr
        ewsrsc = ewsc+ewsrsc
        ewsrvc = ewvc+ewsrvc
        energ2 = energy*energy
        e2srpr = bl2pr*energ2+e2srpr
        e2srsc = bl2sc*energ2+e2srsc
        e2srvc = bl2vc*energ2+e2srvc
        energ3 = energ2*energy
        e3srpr = bl2pr*energ3+e3srpr
        e3srsc = bl2sc*energ3+e3srsc
        e3srvc = bl2vc*energ3+e3srvc
        einrpr = bl2pr/energy+einrpr
        einrsc = bl2sc/energy+einrsc
        einrvc = bl2vc/energy+einrvc
        sumblp = bl2pr+sumblp
        sumbls = bl2sc+sumbls
        sumblv = bl2vc+sumblv
c
        write(6,'(i3,7g11.3)')
     &   imode,energy,bl2pr,ewpr,bl2sc,ewsc,bl2vc,ewvc
        write(14,'(i3,7g11.3)')
     &   imode,energy,bl2pr,ewpr,bl2sc,ewsc,bl2vc,ewvc
        write(61,'(i3,7(1pg20.12))')
     &   imode,energy,bl2pr,ewpr,bl2sc,ewsc,bl2vc,ewvc
        IF(imode==1) THEN
          OPEN(661,file='lowmode')
          WRITE(661,'(1pg20.12)') energy
          CLOSE(661)
        END IF
        belsavpr(imode) = bl2pr
        belsavne(imode) = bl2ne
        belsavsc(imode) = bl2sc
        belsavvc(imode) = bl2vc
      enddo
      write(61,'(1x/1x)')
      sumblv = max(sumblv,xsmall)
      ewsrvc = max(ewsrvc,xsmall)
      einrvc = max(einrvc,xsmall)
      eav2pr = sqrt(abs(e2srpr/sumblp-(ewsrpr/sumblp)**2))
      eav2sc = sqrt(abs(e2srsc/sumbls-(ewsrsc/sumbls)**2))
      eav2vc = sqrt(abs(e2srvc/sumblv-(ewsrvc/sumblv)**2))
      eav3pr = 0.5D0*sqrt(abs(e3srpr/ewsrpr-(ewsrpr/einrpr)))
      eav3sc = 0.5D0*sqrt(abs(e3srsc/ewsrsc-(ewsrsc/einrsc)))
      eav3vc = 0.5D0*sqrt(abs(e3srvc/ewsrvc-(ewsrvc/einrvc)))
      eav31vc = sqrt(abs(e3srvc/ewsrvc))
      write(6,'(a,5x,6(1pg11.3)/)')
     &   ' total:  ',sumblp,ewsrpr,sumbls,ewsrsc,sumblv,ewsrvc
      write(14,'(a,5x,6(1pg11.3))')
     &   ' total:  ',sumblp,ewsrpr,sumbls,ewsrsc,sumblv,ewsrvc
      write(14,'(50x,a)') (factor 4*pi to TRK)
      write(14,'(2(a,5x,3(9x,1pg13.5)))')
     &   ' einv :  ',einrpr,einrsc,einrvc
      CALL centroid_window(omega,maxout)

      write(14,'(2(a,5x,3(9x,1pg13.5))/)')
     &   ' e0m1-av:',sumblp/einrpr,sumbls/einrsc,sumblv/einrvc
      write(14,'(2(a,5x,3(9x,1pg13.5))/)')
     &   ' e10-av: ',ewsrpr/sumblp,ewsrsc/sumbls,ewsrvc/sumblv
      write(14,'(2(a,5x,3(9x,1pg13.5))/)')
     &   ' e21-av: ',e2srpr/ewsrpr,e2srsc/ewsrsc,e2srvc/ewsrvc
      write(14,'(2(a,5x,3(9x,1pg13.5))/)')
     &   ' e1m1-av:',SQRT(ewsrpr/einrpr),SQRT(ewsrsc/einrsc),
     &               SQRT(ewsrvc/einrvc)
      write(14,'(2(a,5x,3(9x,1pg13.5))/)')
     &   ' e20-av :',SQRT(e2srpr/sumblp),SQRT(e2srsc/sumbls),
     &               SQRT(e2srvc/sumblv)
      write(14,'(2(a,5x,3(9x,1pg13.5))/)')
     &   ' e31-av :',SQRT(e3srpr/ewsrpr),SQRT(e3srsc/ewsrsc),
     &               SQRT(e3srvc/ewsrvc)
      write(14,'((a,5x,3(9x,1pg13.5))/(a,5x,4(9x,1pg13.5)))')
     &   ' e2-var: ',eav2pr,eav2sc,eav2vc,
     &   ' e3-var: ',eav3pr,eav3sc,eav3vc
      write(14,'(1x)')
! print the same in different ordering
      WRITE(14,'(2a)')  'averages:      e0m1         e10         e21  ',
     &      '         e1m1         e20          e31    '
      WRITE(14,'(a,6(1pg13.5))') 'isoscalar:',
     &   sumbls/einrsc,ewsrsc/sumbls,e2srsc/ewsrsc,SQRT(ewsrsc/einrsc),
     &   SQRT(e2srsc/sumbls),SQRT(e3srsc/ewsrsc)
      WRITE(14,'(a,6(1pg13.5))') 'isovector:',
     &   sumblv/einrvc,ewsrvc/sumblv,e2srvc/ewsrvc,SQRT(ewsrvc/einrvc),
     &   SQRT(e2srvc/sumblv),SQRT(e3srvc/ewsrvc)
      WRITE(14,'(a)') '# the same averages ordered for plot:'
      WRITE(14,'(a)') '#nr  isoscalar    isovector    moment'
      WRITE(14,'(i3,2(1pg13.5),2x,a)') 
     &          0,sumbls/einrsc,sumblv/einrvc,'e0m1'
      WRITE(14,'(i3,2(1pg13.5),2x,a)') 
     &          1,SQRT(ewsrsc/einrsc),SQRT(ewsrvc/einrvc),'e1m1'
      WRITE(14,'(i3,2(1pg13.5),2x,a)') 
     &          2,ewsrsc/sumbls,ewsrvc/sumblv,'e10'
      WRITE(14,'(i3,2(1pg13.5),2x,a)') 
     &          3,SQRT(e2srsc/sumbls),SQRT(e2srvc/sumblv),'e20'
      WRITE(14,'(i3,2(1pg13.5),2x,a)') 
     &          4,e2srsc/ewsrsc,e2srvc/ewsrvc,'e21'
      WRITE(14,'(i3,2(1pg13.5),2x,a)') 
     &          5,SQRT(e3srsc/ewsrsc),SQRT(e3srvc/ewsrvc),'e31'

      open(33,file='xxpolar')
!      write(33,'(2i5,2f9.4,4(1pg24.16))')
!     &   nprot,nneut,rms,rmsn-rmsp,einrvc,omega(1),sumbls,
!     &   einrvc*(8D0*PI/3D0)*1.44D0
      write(33,'(2i5,f9.4,5(1pg20.12))')
     &   nprot,nneut,rms,rmsn-rmsp,rmsc,omega(1),sumbls,
     &   einrvc*(8D0*PI/3D0)*1.44D0
      close(33)
c
      return
      end
c-----primat-----------------------------------------------------
c
      subroutine primat(upmat,maxmod,text)
c
      include 'params.inc'
      parameter(krpmat=krpmax*(krpmax+1)/2)
c
c     prints matrix upmat given in upper triangular form.
c
      character*(*) text
      dimension upmat(krpmat)
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
      write(6,'(/a)') text
      do imode=1,maxmod
        ipoini = (imode*(imode-1))/2+1
        ipoinf = (imode*(imode+1))/2
        write(6,'(1x,5g12.4)') (upmat(ip),ip=ipoini,ipoinf)
      enddo
c
      return
      end
c-----gsc   --------------------------------------------------------------------
c
      subroutine gsc(qqrpa,qqrpa0,qqhf,pprpa,pprpa0,pphf,
     &               qtrapr,qtrane,maxmod,maxout,
     &               varrpr,var0pr,varrsc,var0sc,varrvc,var0vc,
     &               delsgp,delsgs)
c
      include 'params.inc'
      parameter(krpmat=krpmax*(krpmax+1)/2)
c
c     computes ground state correlations on variance of multipoles
c     and on densities.
c
c     the list parameters are:
c      qqrpa   = q*q expectation over rpa ground state
c                in upper triangular form
c      qqrpa0  = q*q expectation over rpa ground state with h0
c                in upper triangular form
c      qqhf    = q*q expectation over hf ground state
c                in upper triangular form
c      pprpa   = p*p expectation over rpa ground state
c                in upper triangular form
c      pprpa0  = p*p expectation over rpa ground state with h0
c                in upper triangular form
c      pphf    = q*q expectation over hf ground state
c                in upper triangular form
c      qtrapr  = transition moments of initial momenta for protons
c      qtrane  = transition moments of initial momenta for neutrons
c      maxmod  = maximum number of initial ph modes
c      maxout  = maximum number of final rpa modes
c
      dimension   qqrpa(krpmat),qqrpa0(krpmat),qqhf(krpmat),
     &            pprpa(krpmat),pprpa0(krpmat),pphf(krpmat),
     &            qtrapr(krpmax),qtrane(krpmax)
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'result.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
      write(6,'(a)') ' gsc entered'
c
c     compute variance of multipole operators
c
      call varmom(qqrpa,qtrapr,qtrane,maxmod,
     &            maxout,varrsc,varrvc,varrpr)
      call varmom(qqrpa0,qtrapr,qtrane,maxmod,
     &            maxout,var0sc,var0vc,var0pr)
      call varmom(qqhf,qtrapr,qtrane,maxmod,
     &            maxout,varhsc,varhvc,varhpr)
c
c     correlation effects
c
      corblp = varrpr-var0pr
      corbls = varrsc-var0sc
      corblv = varrvc-var0vc
      do iw=6,14,8
        write(iw,'(a/11x,a)')
     &    ' the variance of the multipole operator:',
     &    '     rpa        rpa0         hf      correlation'
        write(iw,'(a,4g12.4)')
     &    ' isoscalar:',varrsc,var0sc,varhsc,corbls,
     &    ' isovector:',varrvc,var0vc,varhvc,corblv,
     &    '    proton:',varrpr,var0pr,varhpr,corblp
      enddo
c
      call sclcor(smsclp,smscln,smsclt,r2l2p,r2l2n,r2l2t)
      if(nprot.gt.0) then
       scalp  = aprot*smsclp/fourpi
       delsgp = corblp/(scalp*scalp) * 0.25*(5.0*rmsp*rmsp/(6.0*pi))
      else
       scalp  = zero
       delsgn = zero
      endif
      scalt  = anumb*smsclt/fourpi
      delsgs = corbls/(scalt*scalt) * 0.25*(5.0*rms*rms/(6.0*pi))
c      write(14,'(2(a,g12.4))')
c     &    ' old delta-sigma**2:  proton=',delsgp,
c     &    '  total=',delsgs
c      write(6,'(2(a,g12.4))')
c     &    ' old delta-sigma**2:  proton=',delsgp,
c     &    '  total=',delsgs
      delsgs = (corbls/varrsc)*delsgs
      delsgp = (corblp/varrpr)*delsgp
      write(14,'(2(a,g12.4))')
     &    ' new delta-sigma**2:  proton=',delsgp,
     &    '  total=',delsgs
      write(6,'(2(a,g12.4))')
     &    ' new delta-sigma**2:  proton=',delsgp,
     &    '  total=',delsgs
c
      delsgs = corbls**2/varrsc *(pi+pi)
     &        / (3.0*mlrpa*anumb*anumb*r2l2t)
      delsgp = corblp**2/varrpr *(pi+pi)
     &        / (3.0*mlrpa*aprot*aprot*r2l2p)
      write(14,'(2(a,g12.4))')
     &    ' delta-sigma**2 from q-scaling: proton=',delsgp,
     &    '  total=',delsgs
c
      return
      end

c-----sclcor --------------------------------------------------------------------
c
      subroutine sclcor(smsclp,smscln,smsclt,r2l2p,r2l2n,r2l2t)
c
      include 'params.inc'
c
c     compute moments of the r**l for the gsc on the rms-radius
c     by the collective scaling model   ( according to f. barranco ).
c
c     input via common's
c     output via list:
c      smsclp = scaling factor for protons
c      smscln = scaling factor for neutrons
c      smsclt = scaling factor for total density
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'result.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
      if(lrpa.ne.0) then
       powmul = alrpa-one
       derfac = 3.0+alrpa
      else
       powmul = one
       derfac = 5.0
      endif
      smsclp = zero
      smscln = zero
      r2l2p  = zero
      r2l2n  = zero
      do ir=2,ngrid
        r      = (ir-1)*rstep
        rpow   = r**powmul
        vol    = r*r*rho(ir)
        smsclp = vol*derfac*r*rpow              + smsclp
        r2l2p  = vol*rpow*rpow                  + r2l2p
        vol    = r*r*rho(ngrid+ir)
        smscln = vol*derfac*r*rpow              + smscln
        r2l2n  = vol*rpow*rpow                  + r2l2n
      enddo
      smsclp = frpirs*smsclp
      smscln = frpirs*smscln
      r2l2p  = frpirs*r2l2p
      r2l2n  = frpirs*r2l2n
      smsclt = (smsclp+smscln)/anumb
      r2l2t  = (r2l2p+r2l2n)/anumb
      if(nprot.gt.0) then
       smsclp = smsclp/aprot
       r2l2p  = r2l2p/aprot
      else
       smsclp = zero
       r2l2p  = zero
      endif
      if(nneut.gt.0) then
       smscln = smscln/aneut
       r2l2n  = r2l2n/aneut
      else
       smscln = zero
       r2l2n  = zero
      endif
      write(6,'(3(a,g12.4))')
     &   ' multipole scaling (3+l)<r**l>: prot=',smsclp,
     &   '  neut=',smscln,'   total=',smsclt
      write(14,'(3(a,g12.4))')
     &   ' multipole scaling (3+l)<r**l>: prot=',smsclp,
     &   '  neut=',smscln,'   total=',smsclt
      write(6,'(3(a,g12.4))')
     &   ' <r**(2l-2)>: prot=',r2l2p,' neut=',r2l2n,'  total=',r2l2t
      write(14,'(3(a,g12.4))')
     &   ' <r**(2l-2)>: prot=',r2l2p,' neut=',r2l2n,'  total=',r2l2t
c
      return
      end
c-----varmom--------------------------------------------------------------------
c
      subroutine varmom(qqvar,qtrapr,qtrane,maxmod,
     &                  maxout,varysc,varyvc,varypr)
c
      include 'params.inc'
      parameter(krpmat=krpmax*(krpmax+1)/2)
c
c     computes variance of multipoles for qq-overlap given on 'qqvar'
c
c     the list parameters are:
c      qqvar   = q*q overlap (input)
c      qtrapr  = transition moments of initial momenta for protons
c      qtrane  = transition moments of initial momenta for neutrons
c      maxmod  = maximum number of initial ph modes
c      maxout  = maximum number of final rpa modes
c      varysc  = variance of scalar multipole operator
c      varyvc  = variance of vector multipole operator
c      varypr  = variance of proton multipole operator
c
      dimension   qqvar(krpmat),
     &            qtrapr(krpmax),qtrane(krpmax)
c
      include 'wavfun.inc'
      include 'constn.inc'
      include 'linres.inc'
c
c------------------------------------------------------------------------------
c
c     compute variance of multipole operators
c
      facpr  = aneut/anumb
      facne  = -aprot/anumb
      varypr = zero
      varyvc = zero
      varysc = zero
      ipoint = 0
      do imode=1,maxmod
        do imod2=1,imode-1
          ipoint = 1 + ipoint
          varysc = qqvar(ipoint)*(qtrapr(imode)+qtrane(imode))
     &              *(qtrapr(imod2)+qtrane(imod2))*2.0 + varysc
          varyvc = qqvar(ipoint)*(facpr*qtrapr(imode)
     &            +facne*qtrane(imode))
     &            *(facpr*qtrapr(imod2)+facne*qtrane(imod2))*2.0
     &            + varyvc
          varypr = qqvar(ipoint)*qtrapr(imode)*qtrapr(imod2)*2.0 
     &            + varypr
        enddo
        ipoint = 1 + ipoint
        varysc = qqvar(ipoint)*(qtrapr(imode)+qtrane(imode))
     &          *(qtrapr(imode)+qtrane(imode)) + varysc
        varyvc = qqvar(ipoint)*(facpr*qtrapr(imode)+facne*qtrane(imode))
     &          *(facpr*qtrapr(imode)+facne*qtrane(imode)) + varyvc
        varypr = qqvar(ipoint)*qtrapr(imode)*qtrapr(imode) + varypr
      enddo
c
c     renormalize the variance to account for different normalization
c     in 'qtrapr' and 'qtrane'.
c
      renorm = (mlrpa*d4pi)**2
      varypr = renorm*varypr
      varysc = renorm*varysc
      varyvc = renorm*varyvc
c
      return
      end
c-----cormhf--------------------------------------------------------------------
c
      subroutine cormhf(qqhf,pphf,workmt,alpha,beta,qpinv,
     &                  maxmod,maxout,nameq,namep,ichr)
c
      include 'params.inc'
      parameter(krpmat=krpmax*(krpmax+1)/2)
c
c     computes ground state correlations on variance of multipoles
c     and on densities.
c
c     the list parameters are:
c      qqhf    = q*q expectation over hf ground state
c                in upper triangular form
c      pphf    = q*q expectation over hf ground state
c                in upper triangular form
c      workmt  = workspace for matrix in upper triangular form
c      alpha   = transformation of coordinates 'q'
c      beta    = transformation of momenta 'p'
c      qpinv   = inverse conjugation matrix
c      omega   = rpa energies
c      maxmod  = maximum number of initial ph modes
c      maxout  = maximum number of final rpa modes
c      nameq   = basis name for initial q-operator
c      namep   = basis name for initial p-operator
c      ichr    = pointer names for the various modes
c
      character*(7) nameq,namep
      character*2 ichr(kichr)
      dimension   workmt(krpmat),qqhf(krpmat),
     &            pphf(krpmat),
     &            alpha(krpmax,krpmax),
     &            beta(krpmax,krpmax),qpinv(krpmax,krpmax),
     &            omega(krpmax)
c
      include 'linres.inc'
      include 'constn.inc'
c
      data ifmeso/1/
c
c------------------------------------------------------------------------------
c
      write(6,'(a)') ' cormhf entered'
c
c    compute q*q correlation matrix for hf state
c         (use 'workmt' workspace)
c
c                                              the overlap of the initial q-op.
      ipoint = 0
      do imode=1,maxmod
        nameq(4:5 ) = ichr(imode)
        call buffql(qlrout,'i',nameq)
        do imod2=1,imode
          ipoint = 1 + ipoint
          nameq(4:5 ) = ichr(imod2)
          call buffql(plrout,'i',nameq)
          workmt(ipoint) = qpover(qlrout,plrout)
        enddo
      enddo
c                                              transform with inverse conj.
      ipoint = 0
      do imode=1,maxmod
      do imod2=1,imode
        ipoint = 1 + ipoint
        sum    = zero
        ip     = 0
        do i1=1,maxmod
          do i2=1,i1-1
            ip     = 1+ip
            sum    = (qpinv(imode,i1)*qpinv(imod2,i2)
     &              + qpinv(imode,i2)*qpinv(imod2,i1))*workmt(ip) + sum
          enddo
          ip     = 1+ip
c         sum    = qpinv(i1,imode)*qpinv(i1,imod2)*workmt(ip) + sum
          sum    = qpinv(imode,i1)*qpinv(imod2,i1)*workmt(ip) + sum
        enddo
        qqhf(ipoint) = sum
      enddo
      enddo
      if(maxmod.eq.1)
     &  write(14,'(a,3g12.3)')
     &   ' qq: <qq>,qpinv,qqhf=',workmt(1),qpinv(1,1),qqhf(1)
c
c    compute p*p correlation matrix for hf state
c
c                                              the overlap of the initial q-op.
      ipoint = 0
      do imode=1,maxmod
        namep(2:3 ) = ichr(imode)
        call buffql(qlrout,'i',namep)
        do imod2=1,imode
          ipoint = 1 + ipoint
          namep(2:3 ) = ichr(imod2)
          call buffql(plrout,'i',namep)
          workmt(ipoint) = qpover(qlrout,plrout)
        enddo
      enddo
c                                              transform with inverse conj.
      ipoint = 0
      do imode=1,maxmod
      do imod2=1,imode
        ipoint = 1 + ipoint
        sum    = zero
        ip     = 0
        do i1=1,maxmod
          do i2=1,i1-1
            ip     = 1+ip
            sum    = (qpinv(i1,imode)*qpinv(i2,imod2)
     &              + qpinv(i2,imode)*qpinv(i1,imod2))*workmt(ip) + sum
          enddo
          ip     = 1+ip
          sum    = qpinv(i1,imode)*qpinv(i1,imod2)*workmt(ip) + sum
        enddo
        pphf(ipoint) = sum
      enddo
      enddo
      if(maxmod.eq.1)
     &  write(14,'(a,3g12.3)')
     &   ' pp: <pp>,qpinv,pphf=',workmt(1),qpinv(1,1),pphf(1)
c
      return
      end
c-----cormrp--------------------------------------------------------------------
c
      subroutine cormrp(qqrpa,pprpa,alpha,beta,qpinv,omega,
     &                  maxmod,maxout)
c
      include 'params.inc'
      parameter(krpmat=krpmax*(krpmax+1)/2)
c
c     computes ground state correlations on variance of multipoles
c     and on densities.
c
c     the list parameters are:
c      qqrpa   = q*q expectation over rpa ground state
c                in upper triangular form
c      pprpa   = q*q expectation over rpa ground state
c                in upper triangular form
c      alpha   = transformation of coordinates 'q'
c      beta    = transformation of momenta 'p'
c      qpinv   = inverse conjugation matrix
c      omega   = rpa energies
c      maxmod  = maximum number of initial ph modes
c      maxout  = maximum number of final rpa modes
c
      dimension   qqrpa(krpmat),
     &            pprpa(krpmat),
     &            alpha(krpmax,krpmax),
     &            beta(krpmax,krpmax),qpinv(krpmax,krpmax),
     &            omega(krpmax)
c
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
      write(6,'(a)') ' cormrp entered'
c
c    compute q*q correlation matrix for rpa state
c
      ipoint = 0
      do imode=1,maxmod
      do imod2=1,imode
        ipoint = 1 + ipoint
        sum    = zero
        do i1=1,maxout
          sum    = beta(imode,i1)*beta(imod2,i1)/(2.0*omega(i1)) + sum
        enddo
c       qqrpa(ipoint) = sum
        qqrpa(ipoint) = sum*half
      enddo
      enddo
c
c    compute p*p correlation matrix for rpa state
c
      ipoint = 0
      do imode=1,maxmod
      do imod2=1,imode
        ipoint = 1 + ipoint
        sum    = zero
        do i1=1,maxout
          sum    = alpha(imode,i1)*alpha(imod2,i1)*(0.5*omega(i1)) + sum
        enddo
c       pprpa(ipoint) = sum
        pprpa(ipoint) = sum*half
      enddo
      enddo
c
      return
      end
c-----quadwd--------------------------------------------------------------------
c
      subroutine quadwd(qqrpa,qqrpa0,maxmod,maxout)
c
      include 'params.inc'
      parameter(krpmat=krpmax*(krpmax+1)/2)
c
c     rewrites widths 'qqrpa' and 'qqrpa0' to the form quadratic
c     in the difference. result is returned such that the
c     old routines 'gscrho' and 'gscenr' can be used as before.
c
c     'qqrpa' and 'qqrpa0' are at entry the inverse widths,
c                              at exit the forefactors of the correlations.
c
      dimension   qqrpa(krpmat),qqrpa0(krpmat),
     &            qinv(krpmax,krpmax),difinv(krpmax,krpmax),
     &            savear(krpmax)
c
      include 'linres.inc'
      include 'wavfun.inc'
      include 'constn.inc'
c
      equivalence (savlr(1,1),qinv(1,1)),(plrout(1,1),difinv(1,1)),
     &            (aux1(1),savear(1))
c
c------------------------------------------------------------------------------
c
c
c     print at entry
c
cc      write(14,'(a,2g12.3)')
cc     &   ' entry quadwd: qqrpa,qqrpa0=',qqrpa(1),qqrpa0(1)
c
c     expand q-matrices for inversion
c
      ipoint = 0
      do imode=1,maxmod
      do imod2=1,imode
        ipoint = 1 + ipoint
        qinv(imod2,imode) = qqrpa(ipoint)
        qinv(imode,imod2) = qqrpa(ipoint)
        difinv(imod2,imode) = qqrpa(ipoint)-qqrpa0(ipoint)
        difinv(imode,imod2) =       difinv(imod2,imode)
      enddo
      enddo
c
c     invert inverse q-width
c
      call  matin2 (qinv,maxmod,krpmax,nerror,determ)
      write(6,'(a,i5,g12.3)')
     &  ' inversion of <qq> with nerreor,determ=',nerror,determ
c
cc      write(14,'(a,g12.3)') ' diff qqrpa-qqrpa0=',difinv(1,1)
cc      write(14,'(a,g12.3)') ' inverse qqrpa=',qinv(1,1)
c
c     product of q-width and difference, result on 'qinv'
c
      do imode=1,maxmod
        do imod2=1,maxmod
          sum = zero
          do i=1,maxmod
            sum = qinv(i,imode)*difinv(imod2,i) + sum
          enddo
          savear(imod2) = sum
        enddo
        do imod2=1,maxmod
          qinv(imod2,imode) = savear(imod2)
        enddo
      enddo
c
cc      write(14,'(a,g12.3)')
cc    &  ' product q*diff=',qinv(1,1)
c
c     product of 'difinv'  and 'qinv' on 'difinv'
c
      do imode=1,maxmod
        do imod2=1,maxmod
          sum = zero
          do i=1,maxmod
            sum = difinv(i,imode)*qinv(imod2,i) + sum
          enddo
          savear(imod2) = sum
        enddo
        do imod2=1,maxmod
          difinv(imod2,imode) = savear(imod2)
        enddo
      enddo
c
cc      write(14,'(a,g12.3)')
cc     &  ' product diff*q*diff=',difinv(1,1)
c
c     compress to triangular form, set 'qqrpa0' to zero
c
      ipoint = 0
      do imode=1,maxmod
      do imod2=1,imode
        ipoint = 1 + ipoint
        qqrpa(ipoint) = difinv(imod2,imode)
        qqrpa0(ipoint) = zero
      enddo
      enddo
c
      return
      end
c-----gscenr--------------------------------------------------------------------
c
      subroutine gscenr(qqrpa,qqrpa0,pprpa,pprpa0,xrpa,
     &                  maxmod,maxout)
c
      include 'params.inc'
      parameter(krpmat=krpmax*(krpmax+1)/2)
c
c     computes ground state correlations on energy
c
c     the list parameters are:
c      qqrpa   = q*q expectation over rpa ground state
c                in upper triangular form
c      qqrpa0  = q*q expectation over rpa ground state with h0
c                in upper triangular form
c      pprpa   = p*p expectation over rpa ground state
c                in upper triangular form
c      pprpa0  = p*p expectation over rpa ground state with h0
c                in upper triangular form
c      xrpa    = workspace for mass- or spring-overlap stored in
c                upper triangular form
c      maxmod  = maximum number of initial ph modes
c      maxout  = maximum number of final rpa modes
c
      dimension   qqrpa(krpmat),qqrpa0(krpmat),
     &            pprpa(krpmat),pprpa0(krpmat),
     &            xrpa(krpmat)
c
      include 'linres.inc'
      include 'gscres.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
      write(6,'(a)') ' gscenr entered'
      maxmat = (maxmod*maxmod+maxmod)/2
c
c     fold the mass matrix with the correlation in p, factor 1/2 included.
c
      read(17,rec=1) (xrpa(i),i=1,maxmat)
      varepp = zero
      enrpp  = zero
      ipoint = 0
      do 110 imode=1,maxmod
      do 120 imod2=1,imode-1
      ipoint = 1 + ipoint
      varepp = -(pprpa(ipoint)-pprpa0(ipoint))*xrpa(ipoint) + varepp
      enrpp  = -pprpa(ipoint)*xrpa(ipoint) + enrpp
  120 continue
      ipoint = 1 + ipoint
      varepp = -half*(pprpa(ipoint)-pprpa0(ipoint))*xrpa(ipoint)
     &        + varepp
      enrpp  = -half*pprpa(ipoint)*xrpa(ipoint) + enrpp
  110 continue
c
c     fold the spring matrix with the correlation in p, factor 1/2 included.
c
      read(17,rec=2) (xrpa(i),i=1,maxmat)
      vareqq = zero
      enrqq  = zero
      ipoint = 0
      do 210 imode=1,maxmod
      do 220 imod2=1,imode-1
      ipoint = 1 + ipoint
      vareqq = (qqrpa(ipoint)-qqrpa0(ipoint))*xrpa(ipoint) + vareqq
      enrqq  = qqrpa(ipoint)*xrpa(ipoint) + enrqq
  220 continue
      ipoint = 1 + ipoint
      vareqq = half*(qqrpa(ipoint)-qqrpa0(ipoint))*xrpa(ipoint)
     &        + vareqq
      enrqq  = half*qqrpa(ipoint)*xrpa(ipoint) + enrqq
  210 continue
c
c     print results
c
      varepp = varepp+varepp
      vareqq = vareqq+vareqq
      enrpp  = enrpp+enrpp
      enrqq  = enrqq+enrqq
      energc(1,irunc) = mlrpa*vareqq
      energc(2,irunc) = mlrpa*varepp
      write(6,'(2(/3(a,g11.4)))')
     &  ' correlation energies:  varepp=',varepp,'  vareqq=',vareqq,
     &  '  total=',mlrpa*(varepp+vareqq),
     &  '    test rpa-energies:  varepp=',enrpp,'  vareqq=',enrqq,
     &  '  total=',enrpp+enrqq
      write(14,'(2(/3(a,g11.4)))')
     &  ' correlation energies:  varepp=',varepp,'  vareqq=',vareqq,
     &  '  total=',mlrpa*(varepp+vareqq)
ctest     &  '    test rpa-energies:  varepp=',enrpp,'  vareqq=',enrqq,
ctest     &  '  total=',enrpp+enrqq
c      write(14,'(a,g11.4)')
c     &  ' correlation energy (without subtraction):',
c     &  mlrpa*(varepp+vareqq)
c
      return
      end
c-----gscrho--------------------------------------------------------------------
c
      subroutine gscrho(dgsrho,dhhrho,qqrpa,qqhf,pprpa,pphf,
     &                  maxmod,nameq,namep,ichr)
c
      include 'params.inc'
      parameter(krpmat=krpmax*(krpmax+1)/2)
c
c     computes coupling of q- and p-response to density.
c     uses
c
c               1               1              1                 1
c     del_rho = - tr{ del_q^2 ----- del_q^2 ------- <p rho p> ------- }
c               8             <q^2>         <q^2>_0           <q^2>_0
c
c               1               1              1                 1
c             + - tr{ del_p^2 ----- del_p^2 ------- <q rho q> ------- }
c               8             <p^2>         <p^2>_0           <p^2>_0
c
c     the <p rho p> and <q rho q> is replaced by assuming redundance
c     thus actually using
c
c
c               1               1
c     del_rho = - tr{ del_q^2 ----- del_q^2  <q rho q> }
c               2             <q^2>
c
c               1               1
c             + - tr{ del_p^2 ----- del_p^2  <p rho p> }
c               2             <p^2>
c
c     redundance is usually guaranteed because both, qqrpa and
c     qqrpa0=qqhf, are generated from an rpa-recipe.
c
c
      character*(7) nameq,namep
      character*2 ichr(kichr)
      dimension   qqrpa(krpmat),qqhf(krpmat),
     &            pprpa(krpmat),pphf(krpmat)
      dimension   dgsrho(kgr2),dhhrho(kgr2)
c     dimension   qrpwid(krpmax)
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'result.inc'
      include 'gscres.inc'
      include 'constn.inc'
c
      logical toccwg
      data toccwg/.true./
c
c------------------------------------------------------------------------------
c
c     clear values before accumulation
c
      do ir=1,ngrid+ngrid
        dhhrho(ir) = zero
        dgsrho(ir) = zero
      enddo
      do n=1,nmaxto
        wgtc(n,irunc) = zero
      enddo
c
c     loop over the initial modes:
c      reading response wave functions from buffer file 18
c      preparing response hamiltonian
c
      do 210 imode=1,maxmod
      nameq(4:5) = ichr(imode)
      namep(2:3) = ichr(imode)
      call buffql(qlrout,'i',nameq)
      call buffql(plrout,'i',namep)
      write(6,'(a,i5)')
     &  ' response in gscrho for mode nr.',imode
c
c     fields for the response hamiltonian
c
c      call qcurre
c      call qpot
c
      do 220 imod2=1,maxmod
      if(imode.le.imod2) then
       ipoint = (imod2*(imod2-1))/2+imode
      else
       ipoint = (imode*(imode-1))/2+imod2
      endif
      widfac = 0.5D0*(qqrpa(ipoint)-qqhf(ipoint))
      if(maxmod.eq.1)
     &   write(14,'(a,3g12.4)') ' widfac: qqrpa,qqhf,diff=',
     &   qqrpa(ipoint),qqhf(ipoint),widfac
c
c     read second response momentum
c
      namep(2:3) = ichr(imod2)
      call buffql(savlr,'i',namep)
      call gscacc(dgsrho,dhhrho,plrout,savlr,widfac)
c
c     now come the contributions from the fluctuation of the p*p width
c
      widfac = 0.5D0*(pprpa(ipoint)-pphf(ipoint))
      if(maxmod.eq.1)
     &   write(14,'(a,3g12.4)') ' widfac: pprpa,pphf,diff=',
     &   pprpa(ipoint),pphf(ipoint),widfac
c
c     read second response momentum
c
      nameq(4:5) = ichr(imod2)
      call buffql(savlr,'i',nameq)
      call gscacc(dgsrho,dhhrho,qlrout,savlr,widfac)
c
c     close loop over both modes
c
  220 continue
  210 continue
c
c     define density at zero
c
      call extrap(dgsrho(1),one)
      call extrap(dgsrho(ngrid+1),one)
      call extrap(dhhrho(1),one)
      call extrap(dhhrho(ngrid+1),one)
c
c     print the results
c
c      sumpp = zero
c      sumqq = zero
c      sumpq = zero
      do ir=1,ngrid
        work1(ir) = dgsrho(ir)+dhhrho(ir)
        ic = ir+ngrid
        work2(ir) = dgsrho(ic)+dhhrho(ic)
      enddo
c      sumpp = sumpp/sumqq
c      sumpq = sumpq/sumqq
c      sumqq = sumqq/sumqq
c      redun = sumpq*sumpq/sumpp
c      write(14,'(a,4g11.3)')
c     & ' gsc: <qq,qq>,<qq,pp>,<pp,pp>,cos=',
c     &   sumqq,sumpq,sumpp,redun
c
c      maxfil = 14
      maxfil = 6
      do iw=6,maxfil,8
        write(iw,'(/a/11x,a,15x,a,15x,a/1x,4a)')
     &    ' correlation densities:',
     &    '  pp-part  ','  hh-part  ','total',
     &    '  r  ',('    proton      neutron   ',i=1,3)
        do ir=1,ngrid
          write(iw,'(1x,f6.2,6g12.3)')   (ir-1)*rstep,
     &    dgsrho(ir),dgsrho(ir+ngrid),dhhrho(ir),dhhrho(ir+ngrid),
     &    work1(ir),work2(ir)
        enddo
      enddo
c
c     formfactor
c
      facp = one/aprot
      if(aneut.ne.zero) facn = one/aneut
      qsurf  = 5.600/rdms
      qrdms   = 4.493/rdms
      if(maxfil.gt.6) then
       write(14,'(/2(a,f8.4)/2a)')
     &  ' formfactors:  qsurf=',qsurf,'  qrdms=',qrdms,
     &  '  iq   q      proton     neutron'
c                                                      all  q
       delqh = 0.5D0*delq
       do iq=1,ngrid
         q1     = iq*delqh-delqh
         dformp = facp*fbform(q1,work1)
         dformn = facn*fbform(q1,work2)
         if(iq.le.ngrid) write(14,'(1x,i3,f6.2,2(2x,6g11.3))')
     &     iq,q1,dformp,dformn
       enddo
      endif
c
      call formcr(dgsrho(1),aprot,rdmsp,surfp,
     &            rdmscr(1,1,1,irunc),surfcr(1,1,1,irunc))
      call formcr(dhhrho(1),aprot,rdmsp,surfp,
     &            rdmscr(1,2,1,irunc),surfcr(1,2,1,irunc))
      if(aneut.ne.0) then
       call formcr(dgsrho(ngrid+1),aneut,rdmsn,surfn,
     &             rdmscr(2,1,1,irunc),surfcr(2,1,1,irunc))
       call formcr(dhhrho(ngrid+1),aneut,rdmsn,surfn,
     &             rdmscr(2,2,1,irunc),surfcr(2,2,1,irunc))
      endif
      write(14,'(a,2g12.4)')
     &   ' delta_sigma from formf.: p/n=',
     &  surfcr(1,1,1,irunc)+surfcr(1,2,1,irunc),
     &  surfcr(2,1,1,irunc)+surfcr(2,2,1,irunc)
      write(14,'(a,2g12.4)')
     &   '       delta_diff.radius: p/n=',
     &  rdmscr(1,1,1,irunc)+rdmscr(1,2,1,irunc),
     &  rdmscr(2,1,1,irunc)+rdmscr(2,2,1,irunc)
c
c     difference in r.m.s.radii, total contribution
c
      dnorpr = zero
      dnorne = zero
      drmspr = zero
      drmsne = zero
      drmsp1 = zero
      drmsp2 = zero
      drmsn1 = zero
      drmsn2 = zero
      do ir=1,ngrid
        r2     = rsquar(ir)
        r4     = r2*r2
        dnorpr = work1(ir)*r2 + dnorpr
        drmspr = work1(ir)*r4 + drmspr
        dnorne = work2(ir)*r2 + dnorne
        drmsne = work2(ir)*r4 + drmsne
        ic     = ir+ngrid
        drmsp1 = dgsrho(ir)*r4 + drmsp1
        drmsp2 = dhhrho(ir)*r4 + drmsp2
        drmsn1 = dgsrho(ic)*r4 + drmsn1
        drmsn2 = dhhrho(ic)*r4 + drmsn2
      enddo
      dnorpr = frpirs*dnorpr
      dnorne = frpirs*dnorne
      drmspr = frpirs*drmspr
      drmsne = frpirs*drmsne
      drmsto = (drmspr+drmsne)/anumb
      drmspr = drmspr/aprot
      drmsp1 = frpirs*drmsp1/aprot
      drmsp2 = frpirs*drmsp2/aprot
      if(aneut.ne.zero) then
       drmsne = drmsne/aneut
       drmsn1 = frpirs*drmsn1/aneut
       drmsn2 = frpirs*drmsn2/aneut
      endif
      delsig(1,1,1,irunc) = drmsp1/3.0
      delsig(2,1,1,irunc) = drmsn1/3.0
      delsig(1,2,1,irunc) = drmsp2/3.0
      delsig(2,2,1,irunc) = drmsn2/3.0
c
      write(6,'(a)') '  now the total corrections:'
      write(14,'(a)') '  now the total corrections:'
      write(6,'(2(a,g12.4))')
     &  ' test normalisation: dnorpr=',dnorpr,' dnorne=',dnrone
      write(6,'(1x,3(a,f9.4)/1x,2(a,f9.5))')
     &    'correlations: drmspr=',drmspr,
     &    ' drmsne=',drmsne,' drmsto=',drmsto,
     &    '   cdelqop=',delqop,'  delq0=',delq0
c
ctest      rnew   = sqrt(abs(rmsp*rmsp+drmspr))
c      drmspr = rnew-rmsp
ctest      crmsp  = rnew
ctest      rnew   = sqrt(abs(rmsn*rmsn+drmsne))
c      drmsne = rnew-rmsn
ctest      crmsn  = rnew
ctest      rnew   = sqrt(abs(rms*rms+drmsto))
c      drmsto = rnew-rms
ctest      crms   = rnew
ctest      write(6,'(1x,3(a,f9.4))') 'new radii:      rmsp=',crmsp,
ctest     &    '   rmsn=', crmsn,'    rms=',  crms
      write(6,'(1x,4(a,g11.3))') ' delta sigma**2 :  proton=',drmspr/3,
     &    ' neutron=',drmsne/3,' total=',drmsto/3
ctest      write(14,'(1x,3(a,f9.4))') 'new radii:      rmsp=',crmsp,
ctest     &    '   rmsn=', crmsn,'    rms=',  crms
      write(14,'(1x,4(a,f8.4))')
     &    ' delta sigma**2 rms:  proton=',drmspr/3,
     &    ' neutron=',drmsne/3,' total=',drmsto/3
c      write(14,'(2(a,2g12.3))')
c     &    ' d_s**2 from <pobsp>_{p/n}:',drmsp1/3,drmsn1/3,
c     &    ' , from <qobsq>_{p/n}:',drmsp2/3,drmsn2/3
c
      write(14,'(/a/a,20(/1x,5g12.4))')
     &   ' correlation on occupation weights:',
     &   ' protons:',(wgtc(na,irunc),na=1,nmaxp)
      sum = 0.0D0
      do na=1,nmaxp
        sum = ph(na)*wgtc(na,irunc) + sum
c        write(14,'(a,i4,4g12.4)')
c     &     ' na,ph,wgtc,prod=',na,ph(na),wgtc(na,irunc),
c     &     ph(na)*wgtc(na,irunc),sum
      enddo
      write(14,'(a,g12.4)') ' sum: ',sum
      write(14,'(a,20(/1x,5g12.4))')
     &   ' neutrons:',(wgtc(na,irunc),na=nmaxp+1,nmaxto)
      sum = 0.0D0
      do na=nmaxp+1,nmaxto
        sum = ph(na)*wgtc(na,irunc) + sum
c        write(14,'(a,i4,4g12.4)')
c     &     ' na,ph,wgtc,prod=',na,ph(na),wgtc(na,irunc),
c     &     ph(na)*wgtc(na,irunc),sum
      enddo
      write(14,'(a,g12.4)') ' sum: ',sum
c
c     test densities
c
      if(maxmod.eq.1 .and. .not.toccwg) then
c
      do ir=1,ngrid+ngrid
        dhhrho(ir) = zero
        dgsrho(ir) = zero
      enddo
      nameq(4:5) = ichr(1)
      namep(2:3) = ichr(1)
      call buffql(qlrout,'i',nameq)
      call buffql(plrout,'i',namep)
      widfac = one
      call gscacc(dgsrho,dhhrho,plrout,plrout,widfac)
      call gscacc(dgsrho,dhhrho,qlrout,qlrout,widfac)
      call extrap(dgsrho(1),one)
      call extrap(dgsrho(ngrid+1),one)
      call extrap(dhhrho(1),one)
      call extrap(dhhrho(ngrid+1),one)
c
c     print the results
c
      sumpp = zero
      sumqq = zero
      sumpq = zero
      do 410 ir=1,ngrid
      r2    = (ir-1)*rstep
      r2    = r2*r2
      ic    = ir+ngrid
      sumpp = (dgsrho(ir)*dgsrho(ir)+dgsrho(ic)*dgsrho(ic))*r2 + sumpp
      sumpq = (dgsrho(ir)*dhhrho(ir)+dgsrho(ic)*dhhrho(ic))*r2 + sumpq
      sumqq = (dhhrho(ir)*dhhrho(ir)+dhhrho(ic)*dhhrho(ic))*r2 + sumqq
  410 continue
      sumpp = sumpp/sumqq
      sumpq = sumpq/sumqq
      sumqq = sumqq/sumqq
      redun = sumpq*sumpq/sumpp
      write(14,'(a,4g11.3)')
     & ' unnorm: <qq,qq>,<qq,pp>,<pp,pp>,cos=',
     &   sumqq,sumpq,sumpp,redun
ctest      else
c
c
c     analyze full density
c
ctest       write(14,'(a/a)')
ctest     &  ' uncorrelated density',
ctest     &  '  r       proton      neutron'
ctest       write(14,'(100(1x,f6.2,2g12.3/))')
ctest     &  ((ir-1)*rstep,rho(ir),rho(ir+ngrid),ir=1,ngrid)
ctest       do 710 ir=1,ngrid
ctest       rho(ir) = work1(ir) + rho(ir)
ctest       rho(ir+ngrid) = work2(ir) + rho(ir+ngrid)
ctest  710  continue
ctest       write(14,'(a/a)')
ctest     &  ' correlated density',
ctest     &  '  r       proton      neutron'
ctest       write(14,'(100(1x,f6.2,2g12.3/))')
ctest     &  ((ir-1)*rstep,rho(ir),rho(ir+ngrid),ir=1,ngrid)
      endif
c
      return
      end
c-----gscacc--------------------------------------------------------------------
c
      subroutine gscacc(dgsrho,dhhrho,awf,bwf,widfac)
c
      include 'params.inc'
c
c     accumulates correlations density
c
      dimension   dgsrho(kgr2),dhhrho(kgr2),
     &            awf(kx11,klrpa),bwf(kx11,klrpa)
      dimension   addwgi(ksl2)
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'result.inc'
      include 'gscres.inc'
      include 'constn.inc'
c
      logical toccpp
      data toccpp/.true./
c
c------------------------------------------------------------------------------
c
      accfac = mlrpa*widfac*d4pi
      do n=1,nmaxto
        addwgi(n) = 0.0D0
      enddo
c
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     accumulate ph-density, contribution from rho_pp
c
      phact = ph(na+nstack/ngrid)
      accact = accfac*phact
      do ir=2,ngrid
        dgsrho(isoadd+ir) = awf(nadd+ir,jdif)*bwf(nadd+ir,jdif)
     &                     *accact    ! accfac
     &                     + dgsrho(isoadd+ir)
      enddo
c
c     prepare p_hh acting on wf, for later use in diagonal-density
c
c      call hrespq(iso,la,xja,nadd,xj,xl,jdif,
c     &            wf(nadd+1),wfs(nadd+1),hlinrs)
c
c     loop over the b-states
c       select only states with the same angular momenta, thus
c       the jdif remains also the same.
c
      do 120 nb=1,na
      jsb   = jsp(nb)
      if(jsb.ne.jsa) goto 120
      wgtb  = weight(nb+nstack/ngrid)
      naddb = nstack + (nb-1)*ngrid
      subfab = spenrg(nb+nstack/ngrid)
c
c     accumulate ph-density, contribution from rho_hh
c
c                                                             momentum-overlap
      sump   = zero
      r      = zero
      do ir=2,ngrid
        r      = rstep + r
        sump   = awf(nadd+ir,jdif)*bwf(naddb+ir,jdif)*(r*r) + sump
      enddo
      sump   = rstep*accfac*sump
      if(na.ne.nb) then
       sump   = sump+sump
      else
       wgtc(na+nstack/ngrid,irunc) = -sump  !*phact
     &               + wgtc(na+nstack/ngrid,irunc)
      endif
c                                                           wocc upper
      if(toccpp .and. na.eq.nb) then
        do nd=1,nmax
          jsd   = jsp(nd)
          wgtd  = weight(nd+nstack/ngrid)
          if(jsd.eq.js .and. wgtd.lt.0.5D0)then
            naddd = nstack + (nd-1)*ngrid
            suma   = zero
            sumb   = zero
            r      = zero
            do ir=2,ngrid
              r      = rstep + r
              suma   = awf(nadd+ir,jdif)*(wf(naddd+ir)*r) + suma
              sumb   = bwf(nadd+ir,jdif)*(wf(naddd+ir)*r) + sumb
            enddo
            suma   = rstep*suma
            sumb   = rstep*sumb
            addwgi(nd+nstack/ngrid) = suma*sumb*ph(na)
     &                        +addwgi(nd+nstack/ngrid)
c            write(14,'(a,2i5,2g12.4)')
c     &       ' na,nd,addwgi=',na,nd,addwgi(nd+nstack/ngrid),suma*sumb
          endif
        enddo
      endif
c                                                             now subtract
      do ir=2,ngrid
        dhhrho(isoadd+ir) =- wf(nadd+ir)*wf(naddb+ir)
     &                     * zentrf(ir)*sump + dhhrho(isoadd+ir)
      enddo
c
c     accumulate diagonal-density, contribution from p_hh
c
c                                                             response-overlap
c      sump   = zero
c      r      = zero
c      do 135 ir=2,ngrid
c      r      = rstep + r
c     sump   = wf(naddb+ir)*hlinrs(ir)*r + sump
c 135 continue
c     sump   = rstep*widfac*sump
c     if(na.ne.nb) sump   = sump+sump
c                                                             now accumulate
c      do 45 ir=2,ngrid
c   45 dhhrho(isoadd+ir) = wf(nadd+ir)*plrout(naddb+ir,jdif)
c     &                    *rinvs(ir)*sump
c     &                   + dhhrho(isoadd+ir)
c
  120 continue
c
c     accumulate diagonal-density, contribution from p_pp
c
c                                                            action of response
c                                                            on second momentum
c      call multr(savlr(nadd+1,jdif),work1(1),rstep,ngrid)
c      parity = (mod(ll,2)*2-1)
c      call derv17(work1(1),work2(1),parity,ngrid,rstep)
c      call hrespq(iso,ll,xj,nadd,xja,xla,jdif,
c     &            work1(1),work2(1),hlinrs)
c      call projec(hlinrs,jsa,nstack,nmax)
c                                                            accumulate
c      do 140 ir=2,ngrid
c  140 dhhrho(isoadd+ir) = -wf(nadd+ir)*hlinrs(ir)*rinvs(ir)*widfac
c     &                   + dhhrho(isoadd+ir)
c
c
      include 'enddo.inc'
c
c     finally set the correlations on occupation weights in the
c     particle domain
c
      if(.not.toccpp) return
C
C     ISOSPIN LOOP
C
      DO ISO=1,2
        IF(ISO.EQ.1) THEN
         NMAX   = NMAXP
         NSTACK = 0
         ISOADD = 0
        ELSE
         NMAX   = NMAXN
         NSTACK = NMAXP*NGRID
         ISOADD = NGRID
        ENDIF
C
C       LOOP OVER THE A-STATES
C
        DO NA=1,NMAX
          npl = NA+NSTACK/NGRID 
          WGT = WEIGHT(npl)
          IF(WGT.LT.0.5D0) then
            wgtc(npl,irunc) = addwgi(npl)*accfac/ph(npl) 
     &                       + wgtc(npl,irunc)
          ENDIF
c          write(14,'(a,2g12.4)')
c     &   ' na,npl,addwgi=',addwgi(npl),wgtc(npl,irunc)
        ENDDO
      ENDDO
c
      return
      end
c-----prigsc--------------------------------------------------------------------
c
      subroutine prigsc
c
      include 'params.inc'
c
c     printout of results from gsc
c
      dimension rsf(2)
c
      include 'gscres.inc'
      include 'wavfun.inc'
      include 'result.inc'
      include 'constn.inc'
c
      data rsf/1.0D0,0.1/
c
c------------------------------------------------------------------------------
c
      open(15,file='corroccup.res')
      write(15,'(a,2i5)')
     &  '#   for NPROT,NNEUT=',nprot,nneut
      do 10 ifile=6,14,8
c
c     first, the full details
c
      write(ifile,'(/a/20x,a,29x,a/11x,a,10x,a,11x,a,10x,a/
     &               1x,a,4x,a,7x,a,8x,a,7x,a)')
     &  ' detailed results on gsc:',
     &  'protons','neutrons',
     &  'direct','scattering','direct','scattering',
     &  'rsf',('pp       hh',i=1,4)
      write(ifile,'(1x,a,2(/f4.1,1x,4f9.5,1x,4f9.5))')
     &  'delta_sigma_rms:',
     &  (rsf(irun),(((delsig(iso,ibox,iscat,irun),ibox=1,2),iscat=1,2)
     &   ,iso=1,2),irun=1,2)
      write(ifile,'(1x,a,2(/f4.1,1x,4f9.5,1x,4f9.5))')
     &  'delta diffraction radius:',
     &  (rsf(irun),(((rdmscr(iso,ibox,iscat,irun),ibox=1,2),iscat=1,2)
     &   ,iso=1,2),irun=1,2)
      write(ifile,'(1x,a,2(/f4.1,1x,4f9.5,1x,4f9.5))')
     &  'delta surface:',
     &  (rsf(irun),(((surfcr(iso,ibox,iscat,irun),ibox=1,2),iscat=1,2)
     &   ,iso=1,2),irun=1,2)
c
c     second, uncorrected correlations
c
      write(ifile,'(1x,a/13x,a,8x,a,5x,a,5x,a,2(/1x,a,2x,3f9.5),
     &             /1x,a,2x,3f9.5,1pg13.5)')
     & 'g.s.correlations for protons, uncorrected:',
     & 'r','R','surface','e',
     & 'direct',delsig(1,1,1,1)+delsig(1,2,1,1),
     &          rdmscr(1,1,1,1)+rdmscr(1,2,1,1),
     &          surfcr(1,1,1,1)+surfcr(1,2,1,1),
     & 'scatt.',delsig(1,1,2,1)+delsig(1,2,2,1),
     &          rdmscr(1,1,2,1)+rdmscr(1,2,2,1),
     &          surfcr(1,1,2,1)+surfcr(1,2,2,1),
     & ' total',delsig(1,1,1,1)+delsig(1,2,1,1)
     &          +delsig(1,1,2,1)+delsig(1,2,2,1),
     & rdmscr(1,1,1,1)+rdmscr(1,2,1,1)+rdmscr(1,1,2,1)+rdmscr(1,2,2,1),
     & surfcr(1,1,1,1)+surfcr(1,2,1,1)+surfcr(1,1,2,1)+surfcr(1,2,2,1),
     & energc(1,1)+energc(2,1)
c
c     third, correlations reduced by pp-diagrams
c
      write(ifile,'(1x,a/13x,a,8x,a,5x,a,5x,a,2(/1x,a,2x,3f9.5),
     &             /1x,a,2x,3f9.5,1pg13.5)')
     & 'g.s.correlations for protons, pp-corrected:',
     & 'r','r','surface','e',
     & 'direct',delsig(1,1,1,1)+delsig(1,2,1,1)-100.0*delsig(1,1,1,2),
     &          rdmscr(1,1,1,1)+rdmscr(1,2,1,1)-100.0*rdmscr(1,1,1,2),
     &          surfcr(1,1,1,1)+surfcr(1,2,1,1)-100.0*surfcr(1,1,1,2),
     & 'scatt.',delsig(1,1,2,1)+delsig(1,2,2,1)-100.0*delsig(1,1,2,2),
     &          rdmscr(1,1,2,1)+rdmscr(1,2,2,1)-100.0*rdmscr(1,1,2,2),
     &          surfcr(1,1,2,1)+surfcr(1,2,2,1)-100.0*surfcr(1,1,2,2),
     & ' total',delsig(1,1,1,1)+delsig(1,2,1,1)-100.0*delsig(1,1,1,2)
     &          +delsig(1,1,2,1)+delsig(1,2,2,1)-100.0*delsig(1,1,2,2),
     & rdmscr(1,1,1,1)+rdmscr(1,2,1,1)+rdmscr(1,1,2,1)+rdmscr(1,2,2,1)
     & -100.0*rdmscr(1,1,1,2)-100.0*rdmscr(1,1,2,2),
     & surfcr(1,1,1,1)+surfcr(1,2,1,1)+surfcr(1,1,2,1)+surfcr(1,2,2,1)
     & -100.0*surfcr(1,1,1,2)-100.0*surfcr(1,1,2,2),
     & energc(1,1)+energc(2,1)-100.0*energc(1,2)-100.0*energc(2,2)
c
c     fourth, correlations reduced by 2.order-diagrams
c
      write(ifile,'(1x,a/13x,a,8x,a,5x,a,5x,a,2(/1x,a,2x,3f9.5))')
     & 'g.s.correlations for protons, 2.order corrected:',
     & 'r','r','surface','e',
     & 'direct',delsig(1,1,1,1)+delsig(1,2,1,1)-100.0*delsig(1,1,1,2)
     &          -100.0*delsig(1,2,1,2),
     &          rdmscr(1,1,1,1)+rdmscr(1,2,1,1)-100.0*rdmscr(1,1,1,2)
     &          -100.0*rdmscr(1,2,1,2),
     &          surfcr(1,1,1,1)+surfcr(1,2,1,1)-100.0*surfcr(1,1,1,2)
     &          -100.0*surfcr(1,2,1,2),
     & 'scatt.',delsig(1,1,2,1)+delsig(1,2,2,1)-100.0*delsig(1,1,2,2)
     &          -100.0*delsig(1,2,2,2),
     &          rdmscr(1,1,2,1)+rdmscr(1,2,2,1)-100.0*rdmscr(1,1,2,2)
     &          -100.0*rdmscr(1,2,2,2),
     &          surfcr(1,1,2,1)+surfcr(1,2,2,1)-100.0*surfcr(1,1,2,2)
     &          -100.0*surfcr(1,2,2,2)
      write(ifile,'(1x,a,2x,3f9.5,1pg13.5)')
     & ' total',delsig(1,1,1,1)+delsig(1,2,1,1)-100.0*delsig(1,1,1,2)
     &          -100.0*delsig(1,2,1,2)
     &         +delsig(1,1,2,1)+delsig(1,2,2,1)-100.0*delsig(1,1,2,2)
     &          -100.0*delsig(1,2,2,2),
     & rdmscr(1,1,1,1)+rdmscr(1,2,1,1)+rdmscr(1,1,2,1)+rdmscr(1,2,2,1)
     & -100.0*rdmscr(1,1,1,2)-100.0*rdmscr(1,1,2,2)
     & -100.0*rdmscr(1,2,1,2)-100.0*rdmscr(1,2,2,2),
     & surfcr(1,1,1,1)+surfcr(1,2,1,1)+surfcr(1,1,2,1)+surfcr(1,2,2,1)
     & -100.0*surfcr(1,1,1,2)-100.0*surfcr(1,1,2,2)
     & -100.0*surfcr(1,2,1,2)-100.0*surfcr(1,2,2,2),
     & energc(1,1)+energc(2,1)-100.0*energc(1,2)-100.0*energc(2,2)
c
c     fifth, correlations reduced by scattering part of pp-diagrams
c
      write(ifile,'(1x,a/13x,a,8x,a,5x,a,5x,a,2(/1x,a,2x,3f9.5),
     &             /1x,a,2x,3f9.5,1pg13.5)')
     & 'g.s.correlations for protons, pp-corrected for scat.part:',
     & 'r','r','surface','e',
     & 'direct',delsig(1,1,1,1)+delsig(1,2,1,1),
     &          rdmscr(1,1,1,1)+rdmscr(1,2,1,1),
     &          surfcr(1,1,1,1)+surfcr(1,2,1,1),
     & 'scatt.',delsig(1,1,2,1)+delsig(1,2,2,1)-100.0*delsig(1,1,2,2),
     &          rdmscr(1,1,2,1)+rdmscr(1,2,2,1)-100.0*rdmscr(1,1,2,2),
     &          surfcr(1,1,2,1)+surfcr(1,2,2,1)-100.0*surfcr(1,1,2,2),
     & ' total',delsig(1,1,1,1)+delsig(1,2,1,1)-100.0*delsig(1,1,1,2)
     &          +delsig(1,1,2,1)+delsig(1,2,2,1)-100.0*delsig(1,1,2,2),
     & rdmscr(1,1,1,1)+rdmscr(1,2,1,1)+rdmscr(1,1,2,1)+rdmscr(1,2,2,1)
     & -100.0*rdmscr(1,1,2,2),
     & surfcr(1,1,1,1)+surfcr(1,2,1,1)+surfcr(1,1,2,1)+surfcr(1,2,2,1)
     & -100.0*surfcr(1,1,2,2),
     & energc(1,1)+energc(2,1)-100.0*energc(1,2)-100.0*energc(2,2)
c
c     sixth, correlations reduced by scat. part of 2.order-diagrams
c
      write(ifile,'(1x,a/13x,a,8x,a,5x,a,5x,a,2(/1x,a,2x,3f9.5))')
     & 'g.s.correlations for protons, scat. part 2.order corrected:',
     & 'r','r','surface','e',
     & 'direct',delsig(1,1,1,1)+delsig(1,2,1,1),
     &          rdmscr(1,1,1,1)+rdmscr(1,2,1,1),
     &          surfcr(1,1,1,1)+surfcr(1,2,1,1),
     & 'scatt.',delsig(1,1,2,1)+delsig(1,2,2,1)-100.0*delsig(1,1,2,2)
     &          -100.0*delsig(1,2,2,2),
     &          rdmscr(1,1,2,1)+rdmscr(1,2,2,1)-100.0*rdmscr(1,1,2,2)
     &          -100.0*rdmscr(1,2,2,2),
     &          surfcr(1,1,2,1)+surfcr(1,2,2,1)-100.0*surfcr(1,1,2,2)
     &          -100.0*surfcr(1,2,2,2)
      write(ifile,'(1x,a,2x,3f9.5,1pg13.5)')
     & ' total',delsig(1,1,1,1)+delsig(1,2,1,1)
     &         +delsig(1,1,2,1)+delsig(1,2,2,1)-100.0*delsig(1,1,2,2)
     &          -100.0*delsig(1,2,2,2),
     & rdmscr(1,1,1,1)+rdmscr(1,2,1,1)+rdmscr(1,1,2,1)+rdmscr(1,2,2,1)
     & -100.0*rdmscr(1,1,2,2)-100.0*rdmscr(1,2,2,2),
     & surfcr(1,1,1,1)+surfcr(1,2,1,1)+surfcr(1,1,2,1)+surfcr(1,2,2,1)
     & -100.0*surfcr(1,1,2,2)-100.0*surfcr(1,2,2,2),
     & energc(1,1)+energc(2,1)-100.0*energc(1,2)-100.0*energc(2,2)
c
c     seventh, print the weights
c
      if(ifile.eq.6) then
        ifilo = 6
      else
        ifilo = 15
      endif
      write(ifilo,'(a/a,11x,a/5x,a,8x,a)')
     &  '# correlations on proton occupation weights:',
     &  '#    na  spenerg   uncorrected     corrected'
      sum = 0.0D0
      summ = 0.0D0
      sump = 0.0D0
      do na=1,nmaxp
        sum = ph(na)*wgtc(na,1) + sum
        if(wgtc(na,1).LT.0.0D0) then
          summ = ph(na)*wgtc(na,1) + summ
        else
          sump = ph(na)*wgtc(na,1) + sump
        endif
        write(ifilo,'(1x,i4,4g12.4)')
     &    na,spenrg(na),
     &    wgtc(na,1),
     &    wgtc(na,1)-100.0*wgtc(na,2)
      enddo
      write(ifilo,'(a,3(1pg12.4)/a/a)') '# sum=',sum,summ,sump,
     &      '     ','     '
      write(ifilo,'(a/4x,a,11x,a/a,8x,a)')
     &  '# correlations on neutron occupation weights:',
     &  '#   na  spenerg    uncorrected     corrected'
      sum = 0.0D0
      summ = 0.0D0
      sump = 0.0D0
      do na=nmaxp+1,nmaxto
        sum = ph(na)*wgtc(na,1) + sum
        if(wgtc(na,1).LT.0.0D0) then
          summ = ph(na)*wgtc(na,1) + summ
        else
          sump = ph(na)*wgtc(na,1) + sump
        endif
        write(ifilo,'(1x,i4,4g12.4)')
     &    na,spenrg(na),
     &    wgtc(na,1),
     &    wgtc(na,1)-100.0*wgtc(na,2)
      enddo
      write(ifilo,'(a,3(1pg12.4)/a/a)') '# sum=',sum,summ,sump,
     &      '     ','     '
c
   10 continue
c
      return
      end
cc-----gscrh2--------------------------------------------------------------------
cc
c      subroutine gscrh2(dgsrho,dhhrho,qrpwid,qpnor,maxmod,nameq,namep,
c     &                  ichr)
cc
c      include 'params.inc'
cc
cc     computes coupling of p-response to rms radii
cc
c      character*(7) nameq,namep
c      character*2 ichr(kichr)
c      dimension   dgsrho(kgr2),dhhrho(kgr2),
c     &            qrpwid(krpmax),qpnor(krpmax)
cc
c      include 'wavfun.inc'
c      include 'grid.inc'
c      include 'linres.inc'
c      include 'result.inc'
c      include 'constn.inc'
cc
cc------------------------------------------------------------------------------
cc
cc     clear values before accumulation
cc
c      do 90 ir=1,ngrid+ngrid
c      dhhrho(ir) = zero
c   90 dgsrho(ir) = zero
cc
cc     loop over the rpa-modes:
cc      reading response wave functions from buffer file 18
cc      preparing response hamiltonian
cc
c      do 210 imode=1,maxmod
c      nameq(5:5) = ichr(imode)
c      namep(5:5) = ichr(imode)
c      call buffql(qlrout,'i',nameq)
c      call buffql(plrout,'i',namep)
cc
cc     fields for the response hamiltonian
cc
c      call qcurre(one)
c      call qpot(one)
cc
c      do 220 imod2=1,maxmod
cc
cc     compute hartree-fock width matrix
cc
c      nameq(5:5) = ichr(imod2)
c      call buffql(savlr,'i',nameq)
c      widfac = -0.25*qpover(qlrout,savlr)
c      if(imod2.eq.imode) widfac = 0.25*qrpwid(imode) + widfac
c      widfac = widfac/(qpnor(imode)*qpnor(imod2))
cc
cc     read second response momentum
cc
c      namep(5:5) = ichr(imod2)
c      call buffql(savlr,'i',namep)
cc
cc     provide necessary quantum numbers and pointers with
cc     do loops 10, 20 and 30
cc
c      include 'initdo.inc'
cc
cc     accumulate ph-density, contribution from rho_pp
cc
c      do 50 ir=2,ngrid
c   50 dgsrho(isoadd+ir) = plrout(nadd+ir,jdif)*savlr(nadd+ir,jdif)
c     &                    *widfac
c     &                   + dgsrho(isoadd+ir)
cc
cc     prepare p_hh acting on wf, for later use in diagonal-density
cc
c      call hrespq(iso,la,xja,nadd,xj,xl,jdif,
c     &            wf(nadd+1),wfs(nadd+1),hlinrs)
cc
cc     loop over the b-states
cc       select only states with the same angular momenta, thus
cc       the jdif remains also the same.
cc
c      do 120 nb=1,na
c      jsb   = jsp(nb)
c      if(jsb.ne.jsa) goto 120
c      wgtb  = weight(nb+nstack/ngrid)
c      naddb = nstack + (nb-1)*ngrid
c      subfab = spenrg(nb+nstack/ngrid)
cc
cc     accumulate ph-density, contribution from rho_hh
cc
cc                                                             momentum-overlap
c      sump   = zero
c      r      = zero
c      do 130 ir=2,ngrid
c      r      = rstep + r
c      sump   = plrout(nadd+ir,jdif)*savlr(naddb+ir,jdif)*r*r + sump
c  130 continue
c      sump   = rstep*widfac*sump
c      if(na.ne.nb) sump   = sump+sump
cc                                                             now subtract
c      do 40 ir=2,ngrid
c   40 dgsrho(isoadd+ir) =- wf(nadd+ir)*wf(naddb+ir)
c     &                   * zentrf(ir)*sump + dgsrho(isoadd+ir)
cc
cc     accumulate diagonal-density, contribution from p_hh
cc
cc                                                             response-overlap
c      sump   = zero
c      r      = zero
c      do 135 ir=2,ngrid
c      r      = rstep + r
c      sump   = wf(naddb+ir)*hlinrs(ir)*r + sump
c  135 continue
c      sump   = rstep*widfac*sump
c      if(na.ne.nb) sump   = sump+sump
cc                                                             now accumulate
c      do 45 ir=2,ngrid
c   45 dhhrho(isoadd+ir) = wf(nadd+ir)*plrout(naddb+ir,jdif)
c     &                    *rinvs(ir)*sump
c     &                   + dhhrho(isoadd+ir)
cc
c  120 continue
cc
cc     accumulate diagonal-density, contribution from p_pp
cc
cc                                                            action of response
cc                                                            on second momentum
c      call multr(savlr(nadd+1,jdif),work1(1),rstep,ngrid)
c      parity = (mod(ll,2)*2-1)
c      call derv17(work1(1),work2(1),parity,ngrid,rstep)
c      call hrespq(iso,ll,xj,nadd,xja,xla,jdif,
c     &            work1(1),work2(1),hlinrs)
c      call projec(hlinrs,jsa,nstack,nmax)
cc                                                            accumulate
c      do 140 ir=2,ngrid
c  140 dhhrho(isoadd+ir) = -wf(nadd+ir)*hlinrs(ir)*rinvs(ir)*widfac
c     &                   + dhhrho(isoadd+ir)
c
cc
c      include 'enddo.inc'
cc
cc     close loop over both modes
cc
c  220 continue
c  210 continue
cc
cc     define density at zero
cc
c      call extrap(dgsrho(1),one)
c      call extrap(dgsrho(ngrid+1),one)
c      call extrap(dhhrho(1),one)
c      call extrap(dhhrho(ngrid+1),one)
cc
cc     print the results
cc
c      do 310 ir=1,ngrid
c      work1(ir) = dgsrho(ir)+dhhrho(ir)
c      ic = ir+ngrid
c      work2(ir) = dgsrho(ic)+dhhrho(ic)
c  310 continue
c      do 320 iw=6,14,8
c      write(iw,'(/a/11x,a,9x,8x,a,8x,11x,a/1x,3a)')
c     &    ' correlation densities:',
c     &    'ph-part','pp-hh-part','total',
c     &    ('    proton      neutron   ',i=1,3)
c      do 330 ir=1,ngrid
c      write(iw,'(1x,f6.2,6g12.3)')   (ir-1)*rstep,
c     &    dgsrho(ir),dgsrho(ir+ngrid),dhhrho(ir),dhhrho(ir+ngrid),
c     &    work1(ir),work2(ir)
c  330 continue
c  320 continue
cc
cc     difference in r.m.s.radii
cc
c      dnorpr = zero
c      dnorne = zero
c      drmspr = zero
c      drmsne = zero
c      do 410 ir=1,ngrid
c      r2     = rsquar(ir)
c      dnorpr = dgsrho(ir)*r2 + dnorpr
c      drmspr = dgsrho(ir)*r2*r2 + drmspr
c      ic     = ir+ngrid
c      dnorne = dgsrho(ic)*r2 + dnorne
c      drmsne = dgsrho(ic)*r2*r2 + drmsne
c  410 continue
c      dnorpr = frpirs*dnorpr
c      dnorne = frpirs*dnorne
c      drmspr = frpirs*drmspr
c      drmsne = frpirs*drmsne
c      drmsto = (drmspr+drmsne)/anumb
c      drmspr = drmspr/aprot
c      drmsne = drmsne/aneut
cc
c      write(6,'(2(a,g12.4))')
c     &  ' test normalisation: dnorpr=',dnorpr,' dnorne=',dnrone
c      write(6,'(1x,3(a,f9.4)/1x,2(a,f9.5))')
c     &    'correlations: drmspr=',drmspr,
c     &    ' drmsne=',drmsne,' drmsto=',drmsto,
c     &    '   cdelqop=',delqop,'  delq0=',delq0
cc
c      rnew   = sqrt(abs(rmsp*rmsp+drmspr))
c      drmspr = rnew-rmsp
c      crmsp  = rnew
c      rnew   = sqrt(abs(rmsn*rmsn+drmsne))
c      drmsne = rnew-rmsn
c      crmsn  = rnew
c      rnew   = sqrt(abs(rms*rms+drmsto))
c      drmsto = rnew-rms
c      crms   = rnew
c      write(6,'(1x,3(a,f9.4))') 'new radii:      rmsp=',crmsp,
c     &    '   rmsn=', crmsn,'    rms=',  crms
c      write(6,'(1x,4(a,f9.4))') 'differences:  drmspr=',drmspr,
c     &    ' drmsne=',drmsne,' drmsto=',drmsto
c      write(14,'(1x,3(a,f9.4))') 'new radii:      rmsp=',crmsp,
c     &    '   rmsn=', crmsn,'    rms=',  crms
c      write(14,'(1x,4(a,f9.4))') 'differences:  drmspr=',drmspr,
c     &    ' drmsne=',drmsne,' drmsto=',drmsto
cc
c      return
c      end
c-----rtrans-------------------------------------------------------------------
c
      subroutine rtrans(vec2,croot,alpha,maxout,maxmod)
c
      include 'params.inc'
c
c     compute transition densities for the 'maxout' coupled modes.
c     use the coupling weights on 'vec2' with 'croot'.
c     the 'maxmod' uncoupled transition denisties are read from file 17.
c
      include 'wavfun.inc'
      include 'linres.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'constn.inc'
      common /bels/belsavpr(krpmax),belsavne(krpmax),
     &             belsavsc(krpmax),belsavvc(krpmax)
c
      dimension croot(krpmax),vec2(krpmax,krpmax),alpha(krpmax,krpmax)
c
c-------------------------------------------------------------------
c
      facpr  = aneut/anumb
      facne  = -aprot/anumb
      do 10 imode=1,maxout
c
      do ir=1,ngrid
        work1(ir) = zero
        work2(ir) = zero
        work3(ir) = zero
        aux1(ir) = zero
        aux2(ir) = zero
        aux3(ir) = zero
      enddo
      do imod2=1,maxmod
        read(17,rec=imod2) (prho(i1),i1=1,ngrid2),
     &         (sigma(i2),i2=1,ngrid2),
     &         (qcur(ir,1),ir=1,ngrid2),
     &         (qcur(ir,2),ir=1,ngrid2)
        reno = D4PI*1D0*mlrpa
        fac    = vec2(imod2,imode)*reno/sqrt(croot(imode)*4.0)
        facc   = alpha(imod2,imode)*reno*sqrt(croot(imode)/4.0)
        do ir=1,ngrid
          work1(ir) = fac*prho(ir) + work1(ir)
          work2(ir) = fac*prho(ir+ngrid) + work2(ir)
          aux1(ir)  = facc*qcur(ir,1) + aux1(ir)
          aux2(ir)  = facc*qcur(ir,2) + aux2(ir)
          aux3(ir)  = facc*qcur(ir+ngrid,1) + aux3(ir)
          work3(ir) = facc*qcur(ir+ngrid,2) + work3(ir)
        enddo
      enddo
c
c     restore accumulated densities on according fields
c
      do ir=1,ngrid
        irn    = ir+ngrid
        prho(ir)    = work1(ir)
        prho(irn)   = work2(ir)
        qcur(ir,1)  = aux1(ir)
        qcur(irn,1) = aux3(ir)
        qcur(ir,2)  = aux2(ir)
        qcur(irn,2) = work3(ir)
      enddo
!
!     normalize transition density to B(EL) according to
!     B(EL) = |int dr r**(l+2) prho|**2
!
      accp = 0D0
      accb = 0D0
      do ir=1,ngrid
        irn    = ir+ngrid
        r4 = ((ir-1)*rstep)**4
        accp = accp + r4*prho(ir)
        accn = accn + r4*prho(ir)
      enddo     
      accp = rstep*accp
      accn = rstep*accn
      facp = SQRT(belsavpr(imode))/accp
      facn = SQRT(belsavne(imode))/accn
      prho(1:ngrid) = facp*prho(1:ngrid)
      prho(ngridp:ngrid2) = facn*prho(ngridp:ngrid2)
      qcur(1:ngrid,1) = facp*qcur(1:ngrid,1)
      qcur(ngridp:ngrid2,1) = facn*qcur(ngridp:ngrid2,1)
      qcur(1:ngrid,2) = facp*qcur(1:ngrid,2)
      qcur(ngridp:ngrid2,2) = facn*qcur(ngridp:ngrid2,2)
c
c     accumulate spin density
c
      if(tspinorbit) then
       do ir=1,ngrid
         work1(ir) = zero
         work2(ir) = zero
       enddo
       do imod2=1,maxmod
         read(17,rec=imod2) (aux1(i1),i1=1,ngrid),
     &           (aux2(i2),i2=1,ngrid),
     &           (sigma(i3),i3=1,ngrid2)
!         fac    = vec2(imod2,imode)/sqrt(croot(imode)*2.0)
         facc   = alpha(imod2,imode)
         do ir=1,ngrid
           work1(ir) = facc*sigma(ir) + work1(ir)
           work2(ir) = facc*sigma(ir+ngrid) + work2(ir)
         enddo
       enddo
c
       do ir=1,ngrid
         irn    = ir+ngrid
         sigma(ir)  = work1(ir)
         sigma(irn) = work2(ir)
       enddo
c
      endif
c
      call pritra(imode,croot(imode))
c
   10 continue
      write(62,'(1x/1x/1x)')
      write(63,'(1x/1x/1x)')
c
      return
      end
c-----pritra-------------------------------------------------------------------
c
      subroutine pritra(imode,omega)
c
      include 'params.inc'
c
c     print transition densities etc from data given on 'prho' and
c     'qcur' in /linres/
c
      include 'wavfun.inc'
      include 'linres.inc'
      include 'result.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
      common /bels/belsavpr(krpmax),belsavne(krpmax),
     &             belsavsc(krpmax),belsavvc(krpmax)
      logical top62
      data top62/.false./
      logical top24,top25
      data top24,top25/.false.,.false./
      data ifmeso/1/
      data gsp,gsn,quenchp,quenchn/5.58D0,-3.82D0,0.68D0,0.64D0/
      parameter (ifmagcurr=1)      ! switch to magnetic current as well
c
c-------------------------------------------------------------------
c
c     effective g-factors
c
      gspe = gsp*quenchp/2D0
      gsne = gsn*quenchn/2D0
      gsT0 = (gsp+gsn)*0.5D0*0.7D0
      gsT1 = (gsn-gsp)*0.5D0*0.7D0
c
c     print the header
c
      if(.not.top62) then
        open(62,file='tor_strength.res')
        write(62,'(2a/a/a)')
     &   '# omega     ISGDR_0   ISGDR_1  toroid_0    toroid_1   ',
     &    ' vortex_0    vortex_1   compres_0 compres_1 vort_01  ',
     &   '   6  40      0.20000     1 80 3  ',
     &   '  X N N  N N  N N  N N  N Y'
        if(ifmagcurr.eq.1) then
          open(63,file='torm_strength.res')
          write(63,'(2a/a/a)')
     &    '# omega     ISGDR_0   ISGDR_1  toroid_0    toroid_1   ',
     &    ' vortex_0    vortex_1   compres_0 compres_1 vort_01  ',
     &    '   6  40      0.20000     1 80 3  ',
     &    '  X N N  N N  N N  N N  N Y'
        endif
        top62 = .true.
      endif
      if(.not.top24) then
        open(24,file='transdens.res')
        open(29,file='transdensj.res')
        open(291,file='transcurp.res')
        open(292,file='transcurrot.res')
        open(293,file='transform.res')
        open(294,file='transcurm.res')
        open(295,file='velocitym.res')
        open(296,file='velocityp.res')
        top24 = .true.
      endif
      if(.not.top25) then
        open(25,file='vorticity.res')
        open(251,file='rhodivvel.res')
        if(ifmagcurr.eq.1) open(26,file='vorticitym.res')
        top25 = .true.
      endif
      if(imode .eq. 1) then
        rewind(24)
        write(24,'(a/a,i5/a)')
     &   '# transition densities for RPA states:',
     &   '# ngrid=',ngrid,
     &   '# r  prho_prot    prho_neut     prho_total    prho_isovect'
        rewind(293)
        write(293,'(a/a,i5/a)')
     &   '# transition formfactor for RPA states:',
     &   '# ngrid=',ngrid,
     &   '# r  form_prot    form_neut     form_total    form_isovect'
        rewind(29)
        write(29,'(a/a,i5/a)')
     &   '# div transition current for RPA states:',
     &   '# ngrid=',ngrid,
     &   '# r  divj_prot    divj_neut     divj_total    divj_isovect'
        rewind(291)
        write(291,'(a/a,i5/a)')
     &   '# upper comp. transition current for RPA states:',
     &   '# ngrid=',ngrid,
     &   '# r    j+_prot      j+_neut       j+_total      j+_isovect'
        rewind(294)
        write(294,'(a/a,i5/a)')
     &   '# lower comp. transition current for RPA states:',
     &   '# ngrid=',ngrid,
     &   '# r    j-_prot      j-_neut       j-_total      j-_isovect'
        rewind(292)
        write(292,'(a/a,i5/a)')
     &   '# rot transition current for RPA states:',
     &   '# ngrid=',ngrid,
     &   '# r  rotj_prot    rotj_neut     rotj_total    rotj_isovect'
        rewind(25)
        write(25,'(a/a,i5/a)')
     &   '# vorticity densities for RPA states:',
     &   '# ngrid=',ngrid,
     &   '# r  vort_prot    vort_neut     vort_total    vort_isovect'
        rewind(251)
        write(251,'(a/a,i5/a)')
     &   '# rho*div(velocity) for RPA states:',
     &   '# ngrid=',ngrid,
     &   '# r rhodivv_prot rhodivv_neut rhodivv_total rhodivv_isovect'
        if(ifmagcurr.eq.1) then
          rewind(26)
          write(26,'(a/a,i5/a)')
     &     '# vorticity densities for RPA states:',
     &     '# ngrid=',ngrid,
     &     '# r  vort_prot    vort_neut     vort_total    vort_isovect'
        endif
        rewind(295)
        write(295,'(a/a,i5/a)')
     &   '# lower comp. transition velocity for RPA states:',
     &   '# ngrid=',ngrid,
     &   '# r    v-_prot      v-_neut       v-_total      v-_isovect'
        rewind(296)
        write(296,'(a/a,i5/a)')
     &   '# upper comp. transition velocity for RPA states:',
     &   '# ngrid=',ngrid,
     &   '# r    v+_prot      v+_neut       v+_total      v+_isovect'
      endif
c
c     obsolete switch
c
      iftest=-1
      if(iftest.eq.1) then
       write(6,'(a)') ' rtrans presently disenabled'
       return
      endif
c
c     test rotation
c
c      call grad(prho,drbmas,1)
c      call rot(drbmas,work2,1)
c      call div(drbmas,work3,1)
c      write(91,'(a)') '#   r   prho   rot*grad div*grad  '
c      do ir=1,ngrid
c        write(91,'(1x,f7.1,5(1pg12.4))')
c     &    (ir-1)*rstep,prho(ir),work2(ir),work3(ir),work3(ir)/prho(ir)
c      enddo
c
c     compute localisation of the mode
c
      rtdens = zero
      rtrms  = zero
      rtrho  = zero
      do ir=2,ngrid
        vol    = (ir-1)*rstep
        vol    = vol*vol
        irn    = ir+ngrid
        add    = vol*(prho(ir)*prho(ir)+prho(irn)*prho(irn))
        rtdens = add + rtdens
        rtrms  = add*vol + rtrms
        rtrho  = add*(rho(ir)+rho(irn)) + rtrho
      enddo
      write(14,'(3(a,g11.3))')
     &  ' localization: rms=',sqrt(rtrms/rtdens),
     &  ' av_rho=',rtrho/rtdens,'  r_tr**2=',rtdens
      if(iftest.eq.0) return
c
c
      write(6,'(a,i3,a/2a)')
     & ' transition density and vorticity for mode nr.',imode,' :',
     & '    r     proton   isoscalar  isovector',
     &         '     div(j)     rot(j)     vorticity'
c     write(14,'(a,i3,a/2a)')
c    & ' transition density and vorticity for mode nr.',imode,' :',
c    & '    r     proton   isoscalar  isovector',
c    &         '     div(j)     rot(j)     vorticity'
      write(14,'(a,i3,a/2a)')
     & ' proton transition density and vorticity for mode nr.',
     & imode,' :',
     & '   r   prho    rho/d(j)     j-        j+    ',
     &         '  div(j)    rot(j)  vorticity'
c
c     note: the continuity equation is modified here by a factor sqrt(omega),
c           there is probably a problem with the normalization of "j",
c           we redefine the current as follows:
c
!      qcur = qcur*sqrt(omega/2D0)       ! move normalization to 'rtrans'?
c
c     prepare rotation of density*spin-density for spin-orbit
c     contribution to effective current.
c
      if(tspinorbit) then
       do ir=1,ngrid
         rp     = rho(ir)
         rn     = rho(ir+ngrid)
         sp     = sigma(ir)
         sn     = sigma(ir+ngrid)
         work2(ir) = (rp+rn)*(sp+sn)
         work3(ir) = rp*sp-rn*sn
       enddo
       if(ifmeso.eq.1) then
         call rotsig(sigma(1),qcurmag(1,1),1)
         call rotsig(sigma(1),qcurmag(1,2),2)
         call rotsig(sigma(ngridp),qcurmag(ngridp,1),1)
         call rotsig(sigma(ngridp),qcurmag(ngridp,2),2)
         parity = 1-mod(lrpa,2)*2
         call extrap(qcurmag(1,1),-parity)
         call extrap(qcurmag(1,2),-parity)
         call extrap(qcurmag(ngridp,1),-parity)
         call extrap(qcurmag(ngridp,2),-parity)
       endif
      endif
c
c     meson exchange correction to isovector current,
c     unchanged isoscalar current,
c     compose proton and neutron current from that.
c
      if(ifmeso.eq.1) then
        facpr  = aneut/anumb
        facne  = -aprot/anumb
        confac = amfc(1)+amfc(2)
        confpr = confac*aneut/anumb
        confne = -confac*aprot/anumb
        b1   = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)     !    =  b1
        rhofac = 2*b1/confac
        confac = 1D0
        do ilamp=1,2
        do ir=1,ngrid
          irn    = ir+ngrid
          curtp      = confac*qcur(ir,ilamp)
     &                +rhofac*(rho(irn)*qcur(ir,ilamp)
     &                     -rho(ir)*qcur(irn,ilamp))
          curtn      = confac*qcur(irn,ilamp)
     &                +rhofac*(rho(ir)*qcur(irn,ilamp)
     &                    -rho(irn)*qcur(ir,ilamp))
          qcur(ir,ilamp)  = curtp
          qcur(irn,ilamp) = curtn
        enddo
        enddo
        if(tspinorbit) then
          rhofac = t4/(amfc(1)+amfc(2))
          do ilamp=1,2
          do ir=1,ngrid
            qcur(ir,ilamp)  = qcur(ir,ilamp) - rhofac*(
     &              rho(irn)*qcurmag(ir,ilamp)
     &              -rho(ir)*qcurmag(irn,ilamp)  )
            qcur(irn,ilamp)  = qcur(irn,ilamp) + rhofac*(
     &              rho(irn)*qcurmag(ir,ilamp)
     &              -rho(ir)*qcurmag(irn,ilamp)  )
          enddo
          enddo
        endif
      endif
c
      parity = 1-mod(lrpa,2)*2
      call extrap(qcur(1,1),-parity)
      call extrap(qcur(1,2),-parity)
      call extrap(qcur(ngridp,1),-parity)
      call extrap(qcur(ngridp,2),-parity)
      if(tspinorbit) then
       do ir=1,ngrid
         x      = aux1(ir)
         aux1(ir) = aux3(ir)
         aux3(ir) = x
       enddo
       call div(aux1,aux2,1)
c       write(14,'(a,10(/1x,5g12.4))')
c     &    ' div(rot(rho*sigma)):',(aux2(ir),ir=1,ngrid)
      endif
c
c     current with magnetic contribution to current on 'qcurmag'.
c
      if(ifmagcurr.eq.1) then
        call rotsig(sigma(1),qcurmag(1,1),1)
        call rotsig(sigma(1),qcurmag(1,2),2)
        call rotsig(sigma(ngridp),qcurmag(ngridp,1),1)
        call rotsig(sigma(ngridp),qcurmag(ngridp,2),2)
        parity = 1-mod(lrpa,2)*2
        call extrap(qcurmag(1,1),-parity)
        call extrap(qcurmag(1,2),-parity)
        call extrap(qcurmag(ngridp,1),-parity)
        call extrap(qcurmag(ngridp,2),-parity)
c        write(14,'(a,10(1pg12.4))')
c     &   'sigma:',sigma(1:10)
c        write(14,'(a,10(1pg12.4))')
c     &   'curmag-1:',qcurmag(1:10,1)
c        write(14,'(a,10(1pg12.4))')
c     &   'curmag-2:',qcurmag(1:10,2)
c        qcurmag(1:ngrid,1:2) = qcur(1:ngrid,1:2)
c     &                        +gspe*qcurmag(1:ngrid,1:2)
c        qcurmag(ngridp:ngrid2,1:2) = qcur(ngridp:ngrid2,1:2)
c     &                            +gsne*qcurmag(ngridp:ngrid2,1:2)
        facpr  = aneut/anumb
        facne  = -aprot/anumb
        do lcomp=1,2
        do ir=1,ngrid
          qcm0 = qcurmag(ir,lcomp)+qcurmag(ir+ngrid,lcomp)
          qcm1 = facpr*qcurmag(ir,lcomp)+facne*qcurmag(ir+ngrid,lcomp)
          qcm0 = gsT0*qcm0*0.5D0
          qcm1 = gsT1*qcm1*0.5D0
          qcurmag(ir,lcomp) = (qcm1-facne*qcm0)
     &                       + qcur(ir,lcomp)
          qcurmag(ir+ngrid,lcomp) = (qcm1-facpr*qcm0)
     &                       + qcur(ir+ngrid,lcomp)
        enddo
        enddo
      endif
c
c     longitudinal, transverse currents and vorticity
c
!      qcur(:,1) = 0D0          ! skip upper/lower componant temporarily
      call div(qcur,dqcur,2)
      call rot(qcur,rqcur,2)
      rsigma = qcur
      do ir=2,ngrid
        qcur(ir,1) = qcur(ir,1)/rho(ir)
        qcur(ir,2) = qcur(ir,2)/rho(ir)
        work3(ir)  = (omega)*prho(ir)/(dqcur(ir)*41.5D0)
        irn    = ir+ngrid
        qcur(irn,1) = qcur(irn,1)/rho(irn)
        qcur(irn,2) = qcur(irn,2)/rho(irn)
      enddo
!      call div(qcur,dqcur,2)
!      call rot(qcur,rqcur,2)
      work3(1)   = zero
      call extrap(qcur(1,1),-parity)
      call extrap(qcur(1,2),-parity)
      call extrap(qcur(ngridp,1),-parity)
      call extrap(qcur(ngridp,2),-parity)
      do iw=295,296
        write(iw,'(a,i5,a,1pg13.5,a,4(1pg13.5))')
     &  '# mode nr.',imode,'  energy=',omega,' bel=',
     &   belsavpr(imode),belsavne(imode),
     &   belsavsc(imode),belsavvc(imode)
      enddo
      facpr  = aneut/anumb
      facne  = -aprot/anumb
      do ir=1,ngrid
        irn    = ir+ngrid
        r      = (ir-1)*rstep
        write(295,'(1x,f6.2,6(1pg13.5))')  r,qcur(ir,1),
     &    qcur(irn,1),qcur(ir,1)+qcur(irn,1),
     &    facpr*qcur(ir,1)+facne*qcur(irn,1)
        write(296,'(1x,f6.2,6(1pg13.5))')  r,qcur(ir,2),
     &    qcur(irn,2),qcur(ir,2)+qcur(irn,2),
     &    facpr*qcur(ir,2)+facne*qcur(irn,2)
      enddo
c
c     divergence of velocity field
c
      write(251,'(a,i5,a,1pg13.5,a,4(1pg13.5))')
     &  '# mode nr.',imode,'  energy=',omega,' bel=',
     &   belsavpr(imode),belsavne(imode),
     &   belsavsc(imode),belsavvc(imode)
      call div(qcur,aux3,2)
      do ir=1,ngrid+ngrid
        aux3(ir) = rho(ir)*aux3(ir)
      enddo
      do ir=1,ngrid
        irn    = ir+ngrid
        r      = (ir-1)*rstep
        write(251,'(1x,f6.2,6(1pg13.5))')  r,aux3(ir),
     &    aux3(irn),aux3(ir)+aux3(irn),
     &    facpr*aux3(ir)+facne*aux3(irn)
      enddo
c
c     rotation of velocity field multiplied by density
c
      call rot(qcur,aux3,2)
      do ir=1,ngrid+ngrid
        aux3(ir) = rho(ir)*aux3(ir)
      enddo
      qcur = rsigma
c
c     compute global moments (toroidal and vortex)
c
      facpr  = aneut/anumb
      facne  = -aprot/anumb
      rfac2  = 1.6666666666667D0*rms*rms  * KCMCORR
      dipt0  = 0D0
      dipt1  = 0D0
      dipj0  = 0D0
      dipj1  = 0D0
      sisg0  = 0D0
      sisg1  = 0D0
      sisd0  = 0D0
      sisd1  = 0D0
      toro0  = 0D0
      toro1  = 0D0
      vort0  = 0D0
      vort1  = 0D0
      vorl0  = 0D0        ! lowest order vorticity
      vorl1  = 0D0        ! lowest order vorticity
c
      diptp = 0D0      
      diptn = 0D0      
      dipjp = 0D0      
      dipjn = 0D0      
c
      do ir=1,ngrid
        r      = (ir-1)*rstep
        r2     = r*r
        r3     = r2*r
        irn    = ir+ngrid
        vorf   = r2*r3
        torf   = r2*(r3-rfac2*r)
        diptp  = r3*prho(ir) + diptp
        diptn  = r3*prho(irn) + diptn
        dipjp  = r3*dqcur(ir) + dipjp
        dipjn  = r3*dqcur(irn) + dipjn
        dipt0  = r3*(prho(ir)+prho(irn)) + dipt0
        dipt1  = r3*(facpr*prho(ir)+facne*prho(irn)) + dipt1
        dipj0  = r3*(dqcur(ir)+dqcur(irn)) + dipj0
        dipj1  = r3*(facpr*dqcur(ir)+facne*dqcur(irn)) + dipj1
        sisg0  = torf*(prho(ir)+prho(irn)) + sisg0
        sisg1  = vorf*(facpr*prho(ir)+facne*prho(irn)) + sisg1
        sisd0  = torf*(dqcur(ir)+dqcur(irn)) + sisd0
        sisd1  = vorf*(facpr*dqcur(ir)+facne*dqcur(irn)) + sisd1
        toro0  = torf*(rqcur(ir,1)+rqcur(irn,1)) + toro0
        toro1  = vorf*(facpr*rqcur(ir,1)+facne*rqcur(irn,1)) + toro1
c        vorf   = r3
c        vort0  = vorf*(aux3(ir)+aux3(irn)) + vort0
c        vort1  = vorf*(facpr*aux3(ir)+facne*aux3(irn)) + vort1
        vort0  = torf*(aux3(ir)+aux3(irn)) + vort0
        vort1  = vorf*(facpr*aux3(ir)+facne*aux3(irn)) + vort1
        vorl0  = r3*(aux3(ir)+aux3(irn)) + vorl0
        vorl1  = r3*(facpr*aux3(ir)+facne*aux3(irn)) + vorl1
      enddo      
      currfac = 0.2103    ! transforms Skyrme current to Schrödinger current
      diptp  = frpirs*diptp/41.5
      diptn  = frpirs*diptn/41.5
      dipjp  = frpirs*dipjp
      dipjn  = frpirs*dipjn
      dipt0  = frpirs*dipt0/41.5
      dipt1  = frpirs*dipt1/41.5
      dipj0  = frpirs*dipj0
      dipj1  = frpirs*dipj1
      sisg0  = frpirs*sisg0/(2*lrpa+3)/2D0
      sisg1  = frpirs*sisg1/(2*lrpa+3)/2D0
      sisd0  = frpirs*sisd0/(2*lrpa+3)/2D0*currfac 
      sisd1  = frpirs*sisd1/(2*lrpa+3)/2D0*currfac
      toro0  = frpirs*toro0*sqrt(lrpa/(1D0+lrpa))/(2*lrpa+3)/2D0*currfac
      toro1  = frpirs*toro1*sqrt(lrpa/(1D0+lrpa))/(2*lrpa+3)/2D0*currfac
      vort0  = frpirs*vort0/(2*lrpa+3)/2D0*currfac
      vort1  = frpirs*vort1/(2*lrpa+3)/2D0*currfac
      write(62,'(1x,f8.2,16(1pg12.4))')
     &  omega,sisg0**2,sisg1**2,(toro0)**2,(toro1)**2,
     &  vort0**2,vort1**2,
     &   sisd0**2,sisd1**2,(toro0-sisd0)**2,(toro1-sisd1)**2,
     &   vorl0**2,vorl1**2,dipt0**2,dipt1**2,dipj0**2,dipj1**2
!      write(65,'(1x,f8.2,16(1pg12.4))') omega,sisd0,sisd1,toro0,toro1
      WRITE(662,'(f8.3,8(1pg13.5))')
     &   omega,
     &   diptp,diptn,dipjp,dipjn,omega*diptp/dipjp,omega*diptn/dipjn,
     &   omega*dipt0/dipj0,omega*dipt1/dipj1

c
c     print arrays
c
      write(24,'(a,i5,a,1pg13.5,a,4(1pg13.5))')
     &  '# mode nr.',imode,'  energy=',omega,' bel=',
     &   belsavpr(imode),belsavne(imode),
     &   belsavsc(imode),belsavvc(imode)
      write(25,'(a,i5,a,1pg13.5,a,4(1pg13.5))')
     &  '# mode nr.',imode,'  energy=',omega,' bel=',
     &   belsavpr(imode),belsavne(imode),
     &   belsavsc(imode),belsavvc(imode)
      write(29,'(a,i5,a,1pg13.5,a,4(1pg13.5))')
     &  '# mode nr.',imode,'  energy=',omega,' bel=',
     &   belsavpr(imode),belsavne(imode),
     &   belsavsc(imode),belsavvc(imode)
      do iw=291,294
        write(iw,'(a,i5,a,1pg13.5,a,4(1pg13.5))')
     &  '# mode nr.',imode,'  energy=',omega,' bel=',
     &   belsavpr(imode),belsavne(imode),
     &   belsavsc(imode),belsavvc(imode)
      enddo

      do ir=1,ngrid
        irn    = ir+ngrid
        r      = (ir-1)*rstep
        write(6,'(1x,f6.2,6g11.3)')  r,
     &    prho(ir),prho(ir)+prho(irn),facpr*prho(ir)+facne*prho(irn),
     &    dqcur(ir),rqcur(ir,1),aux3(ir)
        write(24,'(1x,f6.2,6(1pg13.5))')  r,prho(ir),
     &    prho(irn),prho(ir)+prho(irn),facpr*prho(ir)+facne*prho(irn)
        write(29,'(1x,f6.2,6(1pg13.5))')  r,dqcur(ir),
     &    dqcur(irn),dqcur(ir)+dqcur(irn),
     &    facpr*dqcur(ir)+facne*dqcur(irn)
        write(291,'(1x,f6.2,6(1pg13.5))')  r,qcur(ir,2),
     &    qcur(irn,2),qcur(ir,2)+qcur(irn,2),
     &    facpr*qcur(ir,2)+facne*qcur(irn,2)
        write(294,'(1x,f6.2,6(1pg13.5))')  r,qcur(ir,1),
     &    qcur(irn,1),qcur(ir,1)+qcur(irn,1),
     &    facpr*qcur(ir,1)+facne*qcur(irn,1)
        write(292,'(1x,f6.2,6(1pg13.5))')  r,rqcur(ir,1),
     &    rqcur(irn,1),rqcur(ir,1)+rqcur(irn,1),
     &    facpr*rqcur(ir,1)+facne*rqcur(irn,1)
        write(25,'(1x,f6.2,6(1pg13.5))')  r,aux3(ir),
     &    aux3(irn),aux3(ir)+aux3(irn),facpr*aux3(ir)+facne*aux3(irn)
c        write(14,'(1x,f4.1,1x,1pg9.2,6(1x,g9.2))') r,
c     &    prho(ir),r*r*prho(ir)*prho(ir),rho(ir)
        write(14,'(1x,f4.1,1x,1pg9.2,6(1x,g9.2))') r,
     &    prho(ir),work3(ir),qcur(ir,1),qcur(ir,2),
     &    dqcur(ir),rqcur(ir,1),aux3(ir)*rho(ir)
      enddo
      write(14,'(a,6(1pg12.4))')
     &   ' moments: toroidal (T=0/1), vortex (T=0/1)=',
     &    toro0,toro1,vort0,vort1,sisd0,sisd1
      write(24,'(1x)')
      write(25,'(1x)')
      write(29,'(1x)')
      write(291,'(1x)')
      write(294,'(1x)')
      write(292,'(1x)')
c
      iqmax  = kqmax
      write(6,'(a/a,a)')
     &  ' formfactor:',
     &  '  iq   q        prho       div(j)       rot(j)   ',
     &  ' ff-long   ff-trans'
      write(14,'(a/a,a)')
     &  ' formfactor:',
     &  '  iq   q        prho       div(j)       rot(j)   ',
     &  ' ff-long   ff-trans'
      delql  = 0.5D0*delq
      do iq=1,ngrid
        q1     = iq*delql-delql
        if(q1.eq.zero) q1 = delql/10D0
        r      = zero
        sumd   = zero
        sumr   = zero
        sump   = zero
        sumn   = zero
        do i=2,ngrid
          r      = rstep+r
          vol    = r*r*bessel(r*q1,lrpa)
          sump   = vol*prho(i)+sump
          sumn   = vol*prho(i+ngrid)+sumn
          sumd   = vol*dqcur(i)+sumd
          sumr   = vol*rqcur(i,1)+sumr
        enddo
        sump   = frpirs*sump
        sumn   = frpirs*sumn
        sumd   = frpirs*sumd
        sumr   = frpirs*sumr
        forlon = sumd/q1
        fortra = sumr/q1
        write(6,'(1x,i3,f7.4,5g13.4)')
     &       iq-1,q1,sump,sumd,sumr,forlon,fortra
        write(14,'(1x,i3,f7.4,5(1pg20.12))')
     &       iq-1,q1,sump,sumn,sump+sumn,facpr*sump+facne*sumn
        write(293,'(1x,f7.4,5(1pg13.4))')
     &       q1,sump,sumn,sump+sumn,facpr*sump+facne*sumn
      enddo
c
c     the whole stuff again for the current with magnetic contribution
c
      if(ifmagcurr.ne.1) return
c
c     longitudinal, transverse currents and vorticity
c
!      qcurmag(:,2) = 0D0          ! skip upper/lower componant temporarily
      call div(qcurmag,dqcur,2)
      call rot(qcurmag,rqcur,2)
      do ir=2,ngrid
        qcurmag(ir,1) = qcurmag(ir,1)/rho(ir)
        qcurmag(ir,2) = qcurmag(ir,2)/rho(ir)
        work3(ir)  = prho(ir)/dqcur(ir)
        irn    = ir+ngrid
        qcurmag(irn,1) = qcurmag(irn,1)/rho(irn)
        qcurmag(irn,2) = qcurmag(irn,2)/rho(irn)
      enddo
      work3(1)   = zero
      call extrap(qcurmag(1,1),-parity)
      call extrap(qcurmag(1,2),-parity)
      call extrap(qcurmag(ngridp,1),-parity)
      call extrap(qcurmag(ngridp,2),-parity)
c
c     rotation of velocity field multiplied by density
c
      call rot(qcurmag,aux3,2)
      do ir=1,ngrid+ngrid
        aux3(ir) = rho(ir)*aux3(ir)
      enddo
c
c     compute global moments (toroidal and vortex)
c
      facpr  = aneut/anumb
      facne  = -aprot/anumb
      rfac2  = 1.6666666666667D0*rms*rms  * KCMCORR
      dipt0  = 0D0
      dipt1  = 0D0
      sisg0  = 0D0
      sisg1  = 0D0
      sisd0  = 0D0
      sisd1  = 0D0
      toro0  = 0D0
      toro1  = 0D0
      vort0  = 0D0
      vort1  = 0D0
      vorl0  = 0D0        ! lowest order vorticity
      vorl1  = 0D0        ! lowest order vorticity
      do ir=1,ngrid
        r      = (ir-1)*rstep
        r2     = r*r
        r3     = r2*r
        irn    = ir+ngrid
        torf   = r2*(r3-rfac2*r)
        dipt0  = r3*(prho(ir)+prho(irn)) + dipt0
        dipt1  = r3*(facpr*prho(ir)+facne*prho(irn)) + dipt1
        sisg0  = torf*(prho(ir)+prho(irn)) + sisg0
        sisg1  = torf*(facpr*prho(ir)+facne*prho(irn)) + sisg1
        sisd0  = torf*(dqcur(ir)+dqcur(irn)) + sisd0
        sisd1  = torf*(facpr*dqcur(ir)+facne*dqcur(irn)) + sisd1
        toro0  = torf*(rqcur(ir,1)+rqcur(irn,1)) + toro0
        toro1  = torf*(facpr*rqcur(ir,1)+facne*rqcur(irn,1)) + toro1
        vorf   = r2*r3
c        vorf   = r3
c        vort0  = vorf*(aux3(ir)+aux3(irn)) + vort0
c        vort1  = vorf*(facpr*aux3(ir)+facne*aux3(irn)) + vort1
        vort0  = torf*(aux3(ir)+aux3(irn)) + vort0
        vort1  = torf*(facpr*aux3(ir)+facne*aux3(irn)) + vort1
        vorl0  = r3*(aux3(ir)+aux3(irn)) + vorl0
        vorl1  = r3*(facpr*aux3(ir)+facne*aux3(irn)) + vorl1
      enddo      
      dipt0  = frpirs*dipt0
      dipt1  = frpirs*dipt1
      sisg0  = frpirs*sisg0
      sisg1  = frpirs*sisg1
      sisd0  = frpirs*sisd0/(2*lrpa+3)/2D0
      sisd1  = frpirs*sisd1/(2*lrpa+3)/2D0
      toro0  = frpirs*toro0*sqrt(lrpa/(1D0+lrpa))/(2*lrpa+3)/2D0
      toro1  = frpirs*toro1*sqrt(lrpa/(1D0+lrpa))/(2*lrpa+3)/2D0
      vort0  = frpirs*vort0/(2*lrpa+3)/2D0
      vort1  = frpirs*vort1/(2*lrpa+3)/2D0
      write(14,'(a,6(1pg12.4))')
     &   ' moments (+magnetic): toroidal (T=0/1), vortex (T=0/1)=',
     &    toro0,toro1,vort0,vort1,sisd0,sisd1
      write(63,'(1x,f8.2,12(1pg12.4))')
     &  omega,sisg0**2,sisg1**2,toro0**2,toro1**2,vort0**2,vort1**2,
     &   sisd0**2,sisd1**2,(toro0-sisd0)**2,(toro1-sisd1)**2,
     &   vorl0**2,vorl1**2
      write(26,'(a,i5,a,1pg13.5,a,4(1pg13.5))')
     &  '# mode nr.',imode,'  energy=',omega,' bel=',
     &   belsavpr(imode),belsavne(imode),
     &   belsavsc(imode),belsavvc(imode)
      do ir=1,ngrid
        irn    = ir+ngrid
        r      = (ir-1)*rstep
        write(26,'(1x,f6.2,6(1pg13.5))')  r,aux3(ir),
     &    aux3(irn),aux3(ir)+aux3(irn),facpr*aux3(ir)+facne*aux3(irn)
      enddo
      write(26,'(1x)')
c
      return
      end
c-----invmom------------------------------------------------------
c
      subroutine invmom
c
      return
      end
c-----invmrq------------------------------------------------------
c
      subroutine invmrq(ovqpm1,ovpqb,ovpqb0,resfac,cmfac)
c
      include 'params.inc'
c
c     determines inverse energy-weighted sumrul by itertation
c     of
c            [h,p] = q
c
c     iteration of p is done in place at 'plrout'.
c     the [h,p] is computed using routine 'commhp'.
c     the q is input on 'qlrout'; but the field 'qlrout' is
c     destroyed by the routine.
c     the field 'savlr' is used as workspace.
c
      include 'wavfun.inc'
      include 'potinv.inc'
      include 'skypar.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
c-------------------------------------------------------------------------------
c
c     check option and set parameters
c
      if(itrpa.le.0) return
c
c     prepare mass-field for kinetic-energy inversion
c
      sqfree = sqrt(0.5D0*(amfc(1)+amfc(2)))
      do i=1,ngrid
        sqb(i) = sqfree
      enddo
c
c     save q-operator on 'savlr' because the field 'qlrout' is used
c     as workspace for the result of 'commhp'
c
      qwid   = qnorm(qlrout)
      call copywf(qlrout,savlr)
c
c     initialise plrout
c
      call h0inv(qlrout,plrout,avinv,1)
c     call copywf(qlrout,plrout)
c
c     iterate
c
      do 1000 iter=1,itrpa
c                                      action of [h,p]
      call copywf(savlr,qlrout)
      call commhp(ovpqb,ovpqb0,ovqbqb,ovqbb0,resfac,cmfac,0)
c                                      compute normalisation
      pfac   = qwid/qpover(qlrout,savlr)
c                                      propagate to improved guess for 'plrout'
c     step with h0-inversion
      if(epsrpa.gt.zero) then
       call admlwf(qlrout,savlr,pfac,-one)
       ovdevi = qnorm(qlrout)
       call h0inv(qlrout,qlrout,avinv,1)
       ovdinv = qnorm(qlrout)
       call admlwf(plrout,qlrout,pfac,-abs(epsrpa))
       ovqpm1 = qpover(savlr,plrout)
      else
c     step with kinetic energy inversion
      ovqpm1 = zero
      ovdevi = zero
      ovdinv = zero
c
      include 'initdo.inc'
c
c     compute mismatch on 'aux3'
c     imaginary mismatch on 'aux2' in case of 'ifcmpl.eq.1'
c
      do ir=1,ngrid
        aux3(ir) = (pfac*qlrout(ir+nadd,jdif)-savlr(ir+nadd,jdif))
      enddo
      ovdevi = overlp(aux3(1),aux3(1),rstep,ngrid) + ovdevi
c
c     potential field for kinetic-energy inversion
c
      xlp    = float(ll*ll+ll)
      xlfac  = sqfree**2*xlp
      do i=1,ngrid
        potl(i) = e0damp + xlfac*zentrf(i)
      enddo
c                                   w.f. iteration by h-inversion
      irep  = -1
      call invert(aux3,aux3,irep)
      ovdinv = overlp(aux3(1),aux3(1),rstep,ngrid) + ovdinv
      do ir=1,ngrid
        plrout(nadd+ir,jdif) = pfac*plrout(nadd+ir,jdif)
     &                        -abs(epsrpa)*aux3(ir)
      enddo
c
      ovqpm1 = overlp(savlr(nadd+1,jdif),plrout(nadd+1,jdif),
     &                rstep,ngrid)    +  ovqpm1
c
      include 'enddo.inc'
c
      endif
c                                      prints
      precis = log10(ovdevi/qwid)
      precin = log10(ovdinv/qwid)
      em1    = qwid/ovqpm1
      write(6,'(1x,a,i3,a,2f5.1,a,f6.2,a,g12.4)')
     &   'inv-it=',iter,' prec.=',precis,precin,' pfac=',pfac,
     &   ' e(-1)=',em1
      if(ovdevi/qwid.lt.dabrpa) goto 1009
c
 1000 continue
      write(6,'(a,a,i5,a)') ' inverse momentum not fully converged',
     &   ' within',i4,' iterations!'
      write(14,'(a,g12.4)')
     &  ' invmrq not fully converged. final precision=',precis
 1009 continue
c
      write(6,'(2(a,g12.4))')
     &  ' after invmrq: ovpqb=',ovpqb,'  ovqpm1=',ovqpm1
c
      return
      end
c-----invmrp------------------------------------------------------
c
      subroutine invmrp(ovqpm1,resfac)
c
      include 'params.inc'
c
c     determines inverse of p by itertation
c     of
c            [h,q] = p
c
c     iteration of q is done in place at 'qlrout'.
c     the [h,q] is computed using routine 'commhq'.
c     the p is input on 'plrout'; but the field 'plrout' is
c     destroyed by the routine.
c     the field 'savlr' is used as workspace.
c
      include 'wavfun.inc'
      include 'potinv.inc'
      include 'skypar.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
c-------------------------------------------------------------------------------
c
c     check option and set parameters
c
      if(itrpa.le.0) return
c
c     prepare mass-field for kinetic-energy inversion
c
      sqfree = sqrt(0.5D0*(amfc(1)+amfc(2)))
      do i=1,ngrid
        sqb(i) = sqfree
      enddo
c
c     save q-operator on 'savlr' because the field 'qlrout' is used
c     as workspace for the result of 'commhp'
c
      call copywf(plrout,savlr)
      pwid   = qnorm(plrout)
      write(6,'(a,g12.4)')
     &   ' start invmrp with pwid=',pwid
c
c     initialise qlrout
c
c      if(qnorm(qlrout).eq.zero)
c     & call h0inv(plrout,qlrout,avinv,1)
      call h0inv(plrout,qlrout,avinv,1)
c     call copywf(plrout,qlrout)
c
c     iterate
c
      do 1000 iter=1,itrpa
c                                      action of [h,q]
      call commhq(ovqp,ovqp0,ovpp,ovp0p0,1,resfac,0)
c                                      compute normalisation
      pfac   = pwid/qpover(plrout,savlr)
c                                      propagate to improved guess for 'plrout'
c     step with h0-inversion
      if(epsrpa.gt.zero) then
       call admlwf(plrout,savlr,pfac,-one)
       ovdevi = qnorm(plrout)
       call h0inv(plrout,plrout,avinv,1)
       ovdinv = qnorm(plrout)
       call admlwf(qlrout,plrout,pfac,-abs(epsrpa))
       ovqpm1 = qpover(savlr,qlrout)
      else
c     step with kinetic energy inversion
      ovqpm1 = zero
      ovdevi = zero
      ovdinv = zero
c
      include 'initdo.inc'
c
c     compute mismatch on 'aux3'
c
      do ir=1,ngrid
        aux3(ir) = (pfac*plrout(ir+nadd,jdif)-savlr(ir+nadd,jdif))
      enddo
      ovdevi = overlp(aux3(1),aux3(1),rstep,ngrid) + ovdevi
c
c     potential field for kinetic-energy inversion
c
      xlp    = float(ll*ll+ll)
      xlfac  = sqfree**2*xlp
      do i=1,ngrid
        potl(i) = e0damp + xlfac*zentrf(i)
      enddo
c                                   w.f. iteration by h-inversion
      irep  = -1
      call invert(aux3,aux3,irep)
      ovdinv = overlp(aux3(1),aux3(1),rstep,ngrid) + ovdinv
c
      do ir=1,ngrid
        qlrout(nadd+ir,jdif) = pfac*qlrout(nadd+ir,jdif)
     &                        -abs(epsrpa)*aux3(ir)
      enddo
c
      ovqpm1 = overlp(savlr(nadd+1,jdif),qlrout(nadd+1,jdif),
     &                rstep,ngrid)    +  ovqpm1
c
      include 'enddo.inc'
c
      endif
c                                      prints
      precis = log10(ovdevi/pwid)
      precin = log10(ovdinv/pwid)
      em1    = pwid/ovqpm1
      write(6,'(1x,a,i3,a,2f5.1,a,f6.2,a,g11.4)')
     &   'inv-it=',iter,' prec.=',precis,precin,' pfac=',pfac,
     &   ' e(-2)=',em1
      if(ovdevi/pwid.lt.dabrpa) goto 1009
c
 1000 continue
      write(6,'(a,a,i5,a)') ' inverse momentum not fully converged',
     &   ' within',i4,' iterations!'
      write(14,'(a,g12.4)')
     &  ' invmrp not fully converged. final precision=',precis
 1009 continue
c
c     final analysis
c
      ovqqm2 = qnorm(qlrout)
      pwidm1 = ovqqm2/(4.0*ovqpm1*ovqpm1)
      erpam1 = sign(sqrt(abs(ovqp/ovqpm1)),ovqp*ovqpm1)
      qwidm1 = sign(sqrt(abs(ovqp*ovqpm1)),ovqp*ovqpm1)
      write(6,'(/1x,4(a,g11.3))')
     &    ' invmrp: erpa=',erpam1,' qwidm1=',qwidm1,
     &    ' pwidm1=',pwidm1
c
      return
      end
c-----h0inv ------------------------------------------------------
c
      subroutine h0inv(rwfin,rwfout,avinv,ifprj)
c
      include 'params.inc'
c
c     action of the inverse mean-field hamiltonian for the set
c     of s.p. states in 'wfin'. result is on 'wfout'. the
c     average inverse energy is on 'avinv'.
c
c     some auxiliary fields are needed: 'npa' for intermediate storage
c     of the wave function pointers 'npl'; 'wls' = spin-orbit splitting
c     for given 'js'; 'pkin' = the zentrifugal force for given orbital
c     ang.mom.; 'd1b' = d/dr 'bmass' / r  .
c
      dimension rwfin(kx11,klrpa),rwfout(kx11,klrpa)
c
      include 'wavfun.inc'
      include 'linres.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'numpar.inc'
      include 'potinv.inc'
      include 'constn.inc'
c
c-------------------------------------------------------------------------------
c
      avinv  = zero
      isosav = -1
      if(ifprj .eq. 1)       ovqq   = zero
c
      include 'initdo.inc'
c
      if(iso.ne.isosav) then
c
c      square root of inverse mass 'bmass' on 'sqb';
c      d/dr**2 'sqb' / 'sqb' on 'd2b';
c      d/dr 'bmass' /r on 'd1b'.
c
       do i=1,ngrid
         sqb(i) = sqrt(bmass(i+isoadd))
       enddo
c
       call derv25(sqb(1),d2b(1),one,ngrid,rstep)
       do i=1,ngrid
         d2b(i) = d2b(i)/sqb(i)
       enddo
c
       isosav = iso
      endif
c
c     the effective potential for given 'js' and 'iso' on 'potl';
c     its minimum on 'umin';
c     the zentrifugal force on 'aux3'.
c
      xlp    = (ll*ll+ll)
      xl0    = wls(js)
      d1bdr1 =(16.0*bmass(2+isoadd)-15.0*bmass(1+isoadd)
     *         -bmass(3+isoadd))/(6.0*rstep*rstep)
      potl(1) = pot(1+isoadd)+xl0*potls(1+isoadd)+d1bdr1
      work1(1) = zero
      r      = zero
      do ir=2,ngrid
        r      = rstep+r
        irs    = ir+isoadd
        potl(ir) = pot(irs)+xl0*potls(irs)
     &          +bmass(irs)*xlp*zentrf(ir)+d1b(irs)*rinvs(ir)
        work1(ir) = r*rwfin(nadd+ir,jdif)
      enddo
c
      umin   = subfac
      call invert(work1,work1,-1)
c
c     test inversion
c
c      ifinvt = 0
c      if(ifinvt.eq.1) then
c       write(6,'(3(a,i3),a,g12.4)')
c     &   ' test inversion at iso=',iso,' jsa=',jsa,'  js=',js
c       call hpsi(work1(1),work2(1),aux1(1),aux2(1),
c     &           1,js,iso,subfac)
c       sumin  = zero
c       sumout = zero
c       sumio  = zero
c       r      = zero
c       do 180 ir=2,ngrid
c       r      = rstep+r
c       sumin  = rwfin(nadd+ir,jdif)*rwfin(nadd+ir,jdif)*r*r + sumin
c       sumio  = rwfin(nadd+ir,jdif)*r*work2(ir)  + sumio
c       sumout = work2(ir)*work2(ir)   + sumout
c  180  continue
c       angle  = one - sumio*sumio/max(sumin*sumout,xsmall)
c       write(6,'(a,3g12.4,a,g12.4)')
c     &   '  in,out,in-out=',sumin,sumout,sumio,
c     &   '  angle**2=',angle
c      endif
c      ifwfts = 0
c      if(ifwfts.eq.1) then
c       write(6,'(a,20(/1x,5g12.4))')
c     &  ' rwfin:',(rwfin(nadd+i,jdif)*(i-1)*rstep,i=1,ngrid)
c       write(6,'(a,20(/1x,5g12.4))')
c     &  ' invert:',(work1(i),i=1,ngrid)
c       write(6,'(a,20(/1x,5g12.4))')
c     &  ' hpsi*invert:',(work2(i),i=1,ngrid)
c      endif
c
c     average inverse energy
c
      sum    = zero
      r      = zero
      do i=2,ngrid
        r      = rstep+r
        sum    = r*rwfin(nadd+i,jdif)*work1(i) + sum
        rwfout(nadd+i,jdif) = rinvs(i)*work1(i)
      enddo
      avinv  = rstep*sum + avinv
c
      if(ifprj .eq. 1) then
       parity = 1-2*mod(la+lrpa,2)
       call extrap(rwfout(nadd+1,jdif),parity)
       call projec(rwfout(nadd+1,jdif),js,nstack,nmax)
       ovqq   = overlp(rwfout(nadd+1,jdif),rwfout(nadd+1,jdif),
     &                 rstep,ngrid)   +  ovqq
      endif
c
c     testprints
c
c     write(6,'(3(a,i3),a,g12.4)')
c    &  ' at iso=',iso,' jsa=',jsa,'  js=',js,' :  avinv=',rstep*sum
c     sum    = zero
c     r      = zero
c     do 190 ir=2,ngrid
c     r      = rstep+r
c 190 sum    = r*wf(nadd+ir)*rwfin(nadd+ir,jdif) + sum
c     write(6,'(21x,2(a,g12.4))')
c    &   '  subfac=',subfac,'   overlap=',rstep*sum
c
      include 'enddo.inc'
c
      return
      end
c-----h0act ------------------------------------------------------
c
      subroutine h0act(rwfin,rwfout,avh0,ifprj)
c
      include 'params.inc'
c
c     action of the mean-field hamiltonian for the set
c     of s.p. states in 'wfin'. result is on 'wfout'. the
c     average h0-energy is on 'avh0'.
c
c     some auxiliary fields are needed: 'npa' for intermediate storage
c     of the wave function pointers 'npl'; 'wls' = spin-orbit splitting
c     for given 'js'; 'pkin' = the zentrifugal force for given orbital
c     ang.mom.; 'd1b' = d/dr 'bmass' / r  .
c
      dimension rwfin(kx11,klrpa),rwfout(kx11,klrpa)
c
      include 'wavfun.inc'
      include 'linres.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'numpar.inc'
      include 'potinv.inc'
      include 'constn.inc'
c
c-------------------------------------------------------------------------------
c
      avh0   = zero
      isosav = -1
      if(ifprj .eq. 1)       ovqq   = zero
c
      include 'initdo.inc'
c
c     action of mean-field hamiltonian on actual wave function
c
      call multr(rwfin(nadd+1,jdif),work1,rstep,ngrid)
      call hpsi(work1(1),work2(1),aux1(1),aux2(1),
     &          1,js,iso,subfac)
c
c     average h0-energy and 1/r scaling of resulting wave function
c
      sum    = zero
      r      = zero
      do i=2,ngrid
        r      = rstep+r
        sum    = r*rwfin(nadd+i,jdif)*work2(i) + sum
        rwfout(nadd+i,jdif) = rinvs(i)*work2(i)
      enddo
      avh0   = rstep*sum + avh0
c
      if(ifprj .eq. 1) then
       parity = 1-2*mod(la+lrpa,2)
       call extrap(rwfout(nadd+1,jdif),parity)
       call projec(rwfout(nadd+1,jdif),js,nstack,nmax)
       ovqq   = overlp(rwfout(nadd+1,jdif),rwfout(nadd+1,jdif),
     &                 rstep,ngrid)   +  ovqq
      endif
c
      include 'enddo.inc'
c
      return
      end
c-----cminit------------------------------------------------------
c
      subroutine cminit
c
      include 'params.inc'
c
c     action of the c.m. momentum on-the s.p wave function
c      produces result on 'plraj' overwriting old value.
c      overlap with old 'plrout' is accummulated on 'ovpq' and
c      norm of new 'plrout' on 'ovpp'.
c
      include 'wavfun.inc'
      include 'result.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
c
c-------------------------------------------------------------------
c
      ovqp   = zero
      ovpp   = zero
      ovoo   = zero
c
      if(lrpa.ne.1 .and. isorpa.ne.1) stop 'wrong call to "cmlraj"'
c
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     momentum from [t,q]
c
       redmat = redume(xj,alrpa,xja)
       do ir=1,ngrid
         work1(ir) = redmat*qrad(isoadd+ir)*wf(nadd+ir)
       enddo
       call tpsi(work1(1),work3(1),aux2(1),js,iso,zero)
       call tpsi(wf(nadd+1),work1(1),aux2(1),jsa,iso,zero)
       do ir=2,ngrid
         work3(ir) =
     &            (work3(ir)-redmat*qrad(isoadd+ir)*work1(ir))
     &            *rinvs(ir)
         savlr(nadd+ir,jdif) = work3(ir)
       enddo
      parll  = (1-mod(ll,2)*2)
      call extrap(savlr(nadd+1,jdif),parll)
c
c     action of the kinetic energy on 'q->wf', store on 'hlinrs'
c
      iphase = xl+xja-half
      redfac = (1-2*mod(iphase,2))*sqrt((xj+xj+one)*(xja+xja+one))
     *        * sixjot(xl,xj,half,xja,xla,one)
      if(xl.gt.xla) then
       facder = redfac*sqrt(xla+one)
       facrad = -facder*(xla+one)
      else
       facder = -redfac*sqrt(xla)
       facrad = facder*xla
      endif
      do ir=2,ngrid
        rinv   = rinvs(ir)
        work1(ir) = plrout(nadd+ir,jdif)
        plrout(nadd+ir,jdif) =
     &        (facder*wfs(nadd+ir)+facrad*wf(nadd+ir)*rinv)*rinv
      enddo
      parll  = (1-mod(ll,2)*2)
      call extrap(plrout(nadd+1,jdif),parll)
c
c     project out all wave function with same 'js'
c
      call projec(plrout(nadd+1,jdif),js,nstack,nmax)
      call projec(savlr(nadd+1,jdif),js,nstack,nmax)
c
c     accumulate overlaps
c
      ovoo   = overlp(savlr(nadd+1,jdif),savlr(nadd+1,jdif),
     &                rstep,ngrid)+ovoo
      ovqp   = overlp(savlr(nadd+1,jdif),plrout(nadd+1,jdif),
     &                rstep,ngrid) + ovqp
      ovpp   = overlp(plrout(nadd+1,jdif),plrout(nadd+1,jdif),
     &                rstep,ngrid)    +  ovpp
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
      cosang = ovqp*ovqp/(ovoo*ovpp)
      write(14,'(a,4g12.3)')
     &  ' <n|n>,<n|o>,<o,o>,cos=',ovpp,ovqp,ovoo,cosang
c
      return
      end
c-----commhq--------------------------------------------------------------------
c
      subroutine commhq(ovqp,ovqp0,ovpp,ovp0p0,ifprj,resfac,ifsav0)
c
      include 'params.inc'
c
c     computes p = [h,q] for given q
c     with p and q given by its action on wave functions:
c         q |wf> <--> 'qlrout'
c         p |wf> <--> 'plrout'
c
c     the response wave functions are handled via common /linres/.
c     they mean:
c         'qlrout' = given q = input to routine
c         'plrout' = resulting [h,q] = output of routine
c
c     during evaluation of [h,q] some overlaps are accumulated:
c      ovqp    = <wf|q [h,q]|wf>
c      ovqp0   = <wf|q [h0,q]|wf>
c      ovpp    = <wf|p p|wf> = <wf|[q,h] [h,q]|wf>
c      ovp0p0  = <wf|[q,h0] [h0,q]|wf>
c     where h is the full hamiltonian and h0 the mean-field part
c
c     further list parameters are:
c      ifprj   = switch for ph-projection of the final [h,q]
c      resfac  = weight factor for the residual interaction;
c                should be 1.0
c      ifsav0  = switch for saving the h0-commutator on 'savlr'
c
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
      data xmaxlw/1.0D-20/
c
c----------------------------------------------------------------------
c
c      write(6,'(a,i5,g12.4,i5)')
c     &  ' ifprj,resfac,ifsav0=',ifprj,resfac,ifsav0
c
c     prepare the response current 'qcur' and its derivative 'dqcur'
c
      call qcurre(resfac)
c      write(6,'(a)') ' after QCURRE'
      if(.not.tyukaw .and. ayuk.ne.zero) then
c       write(6,'(a)') ' damping switch'
       dampay = ayuk*ayuk
c       do 915 idamp=1,8
       call t0inv(qcur(1,1),qcur(1,1),work1,dampay,lrpa-1)
       call t0inv(qcur(ngrid+1,1),qcur(ngrid+1,1),work1,dampay,lrpa-1)
       call t0inv(qcur(1,2),qcur(1,2),work1,dampay,lrpa+1)
       call t0inv(qcur(ngrid+1,2),qcur(ngrid+1,2),work1,dampay,lrpa+1)
       if(imfold.eq.1) then
        call t0inv(qcur(1,1),qcur(1,1),work1,dampay,lrpa-1)
        call t0inv(qcur(ngrid+1,1),qcur(ngrid+1,1),work1,dampay,lrpa-1)
        call t0inv(qcur(1,2),qcur(1,2),work1,dampay,lrpa+1)
        call t0inv(qcur(ngrid+1,2),qcur(ngrid+1,2),work1,dampay,lrpa+1)
       endif
       call rot(qcur,rqcur,2)
c  915  continue
      elseif(resdmp.ne.zero) then
       call filter(sigma(1),sigma(1),lrpa)
       call filter(sigma(ngrid+1),sigma(ngrid+1),lrpa)
       call filter(qcur(1,1),qcur(1,1),lrpa-1)
       call filter(qcur(ngrid+1,1),qcur(ngrid+1,1),lrpa-1)
       call filter(qcur(1,2),qcur(1,2),lrpa+1)
       call filter(qcur(ngrid+1,2),qcur(ngrid+1,2),lrpa+1)
       call rot(qcur,rqcur,2)
      endif
c
c     prepare the kinetic response potential for given current
c
      call qpot(resfac)
!      write(14,'(a)') ' QPOT'
      if(.not.tyukaw .and. ayuk.ne.zero) then
       call t0inv(rkin(1,1),rkin(1,1),work1,dampay,lrpa-1)
       call t0inv(rkin(ngrid+1,1),rkin(ngrid+1,1),work1,dampay,lrpa-1)
       call t0inv(rkin(1,2),rkin(1,2),work1,dampay,lrpa+1)
       call t0inv(rkin(ngrid+1,2),rkin(ngrid+1,2),work1,dampay,lrpa+1)
       call div(rkin,drkin,2)
      endif
c
c     prepare action of the q-response hamiltonian on the wave
c     function. result on 'hlinrs'.
c
      ovhh   = zero
      ovtt   = zero
      ovht   = zero
c
      ovqp0  = zero
      ovp0p0 = zero
      ovqp   = zero
      ovpp   = zero
c
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     commutator of hf-hamiltonian on q-wavefunction
c
c     the computation via h is overriden. there seems to
c     be a problem with the consistence in nabla -> qcur.
c     use commutator with t instead (which is o.k. for
c     local q).
c
c
c     set test switches.
c     !! attention: test switches are deactivated by 'c' columns.
c
      ift    = 0
      ifttst = 0
c
      parll  = (1-mod(ll,2)*2)
c
      call multr(qlrout(nadd+1,jdif),work1(1),rstep,ngrid)
      call hpsi(work1(1),work2(1),aux1(1),aux2(1),1,js,iso,subfac)
      do ir=2,ngrid
        work2(ir) = work2(ir)*rinvs(ir)
      enddo
      call extrap(work2(1),parll)
c
c     action of linear response hamiltonian
c
ctest      ifnocu = 1
ctest      if(resfac.ne.zero .and. ifnocu.ne.1) then
      if(resfac.ne.zero) then
        call hrespq(iso,la,xja,nadd,xj,xl,jdif,
     &              wf(nadd+1),wfs(nadd+1),hlinrs)
      else
        do ir=1,ngrid
          hlinrs(ir) = 0D0
        enddo
      endif
cc
cc      compute p using [t,q].
cc
c      if(ift.eq.1 .or. ifttst.eq.1) then
c       redmat = redume(xj,alrpa,xja)
c       do ir=1,ngrid
c         work1(ir) = redmat*qrad(isoadd+ir)*wf(nadd+ir)
c        enddo
c       call tpsi(work1(1),work3(1),aux2(1),js,iso,zero)
c       call tpsi(wf(nadd+1),work1(1),aux2(1),jsa,iso,zero)
c       do 90 ir=2,ngrid
c       work3(ir) =
c     &          (work3(ir)-redmat*qrad(isoadd+ir)*work1(ir))
c     &         *rinvs(ir)
c       plrout(nadd+ir,jdif) = work3(ir)
c   90  continue
c      endif
c
c     compute p using [h,q] and hlinrs
c
      if(ift.ne.1) then
       do ir=2,ngrid
         plrout(nadd+ir,jdif) = work2(ir)+ resfac * hlinrs(ir)
       enddo
       parll  = (1-mod(ll,2)*2)
       if(ifsav0.eq.1) then
        do ir=2,ngrid
          savlr(nadd+ir,jdif) = work2(ir)
        enddo
        call extrap(savlr(nadd+1,jdif),parll)
        if(ifprj.eq.1)
     &   call projec(savlr(nadd+1,jdif),js,nstack,nmax)
       endif
      endif
      call extrap(plrout(nadd+1,jdif),parll)
cc
cc     test work2-work3 versus hlinrs
cc
c      if(ifttst.eq.1) then
c       call hrespq(iso,la,xja,nadd,xj,xl,jdif,
c     &             wf(nadd+1),wfs(nadd+1),hlinrs)
cc
cc      check overlap of [h,q] versus [t,q]
cc
c       if(ifprj.eq.1) call projec(work3(1),js,nstack,nmax)
c       if(ifprj.eq.1) call projec(work2(1),js,nstack,nmax)
c       if(ifprj.eq.1) call projec(hlinrs(1),js,nstack,nmax)
cc      write(6,'(3(a,i2)/a,a)')
cc    &    ' testprint of [h,q] and hlinrs for iso=',
cc    &    iso,'  na=',na,'  jdif=',jdif,
cc    &      '   r    [h-t,q]     hlinrs      quotient    [h,q]  [t,q]'
cc      do 1010 ir=2,ngrid/2
cc      write(6,'(1x,f4.1,5g12.5)') (ir-1)*rstep,work2(ir)-work3(ir),
cc    &       hlinrs(ir),
cc    &   (work2(ir)-work3(ir))
cc    &    /sign(max(abs(hlinrs(ir)),xmaxlw),hlinrs(ir)),
cc    &    work2(ir),work3(ir)
cc1010  continue
cc
cc      test overlap of [t,q] versus [h,q]+hlinrs
cc
c       do 1020 ir=1,ngrid
c 1020  work1(ir) = work2(ir)+hlinrs(ir)
c       xno2   = overlp(work1(1),work1(1),rstep,ngrid)
c       xno3   = overlp(work3(1),work3(1),rstep,ngrid)
c       xno23  = overlp(work1(1),work3(1),rstep,ngrid)
c       cosang = (1.0-xno23*xno23/max(xno2*xno3,xmaxlw))
c       write(6,'(3(a,g11.4))')
c     &    ' compare p: norm-h=',xno2,'  norm-t=',xno3,
c     &    '  angle=',cosang
c       ovhh   = xno2+ovhh
c       ovtt   = xno3+ovtt
c       ovht   = xno23+ovht
c      endif
c
      if(ifprj.eq.1)
     & call projec(plrout(nadd+1,jdif),js,nstack,nmax)
      if(ifprj.eq.1) call projec(work2(1),js,nstack,nmax)
c
      ovqp0  = overlp(qlrout(nadd+1,jdif),work2(1),rstep,ngrid)
     &        + ovqp0
      ovp0p0 = overlp(work2(1),work2(1),rstep,ngrid)  + ovp0p0
      ovqp   = overlp(qlrout(nadd+1,jdif),plrout(nadd+1,jdif),
     &                rstep,ngrid)    +  ovqp
      ovpp   = overlp(plrout(nadd+1,jdif),plrout(nadd+1,jdif),
     &                rstep,ngrid)    +  ovpp
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
cc
cc     testprint t <--> h
cc
c      if(ifttst.eq.1) then
c       cosang = (1.0-ovht*ovht/max(ovtt*ovhh,xmaxlw))
c       write(6,'(3(a,g11.4))')
c     &     ' total p: norm-h=',ovhh,'  norm-t=',ovtt,
c     &     '  angle=',cosang
c      endif
c
      return
      end
c-----hrespq-----------------------------------------------------------------
c
      subroutine hrespq(iso,la,xja,nadd,xj,xl,jdif,
     &                  wfin,wfsin,hwfout)
c
      include 'params.inc'
c
c     action of linear response hamiltonian h_q (for q-like perturbations)
c     on the single particle h.f. state from 'wf' in common /wavfun/
c     according to the quantum numbers
c      iso    = isospin
c      la     = orbital angular momentum of 'wf'
c      xja    = total angular momentum of 'wf'
c      xl     = orbital angular momentum of result on 'hwfout'
c      xlj    = total angular momentum of result on 'hwfout'
c      jdif   = difference quantum number between angular momenta
c
c     the response potentials 'rkin', 'drkin' have been computed
c     in 'qpot' and are transferred via common /linres/.
c     the resonse current 'rqcur' has been computed in 'qcurre'
c     and is transferred via common /linres/.
c
c     the fields for the mean-field hamiltonian are transferred
c     via common /wavfun/.
c
c     the input fields are:
c      wfin    = wave function on which the response hamiltonian acts
c      wfsin   = derivative of the wave function on which
c                the response hamiltonian acts
c     these input fields are the r-weighted wave functions.
c
c     the result is:
c      hwfout  = response hamiltonian on 'wfin'
c     the output field is not r-weighted (as all response wave functions).
c
      dimension wfin(kg1),wfsin(kg1),hwfout(kg1)
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      data ifnodd/0/
c
c-----------------------------------------------------------------------
c
      xla   = la
      ll    = xl
      isoadd = ngrid*(iso-1)
c
c     clear accumulator field
c
      do ir=1,ngrid
        hwfout(ir) = zero
      enddo
      if(ifsepi.eq.1 .or. ifnodd.eq.1) return
c
c     loop over  lambda' = l+-1
c
      do 20 ilamp=1,2
      lamp   = lrpa-3+2*ilamp
      if(lamp.ge.0) then
       alamp  = lamp
c
c      derivative on wave function
c
       call redcur(xj,xl,alrpa,alamp,xja,xla,facder,facrad)
       do ir=2,ngrid
         rinv   = rinvs(ir)
         hwfout(ir) = rkin(isoadd+ir,ilamp)
     &               *(facder*wfsin(ir)+facrad*wfin(ir)*rinv)*rinv
     &               + hwfout(ir)
       enddo
      endif
c
   20 continue
c
c      derivative on 'rkin'
c
      redmat = 0.5D0*redume(xj,alrpa,xja)
      isoadd = ngrid*(iso-1)
      do ir=2,ngrid
        rinv   = rinvs(ir)
        hwfout(ir) = redmat*drkin(isoadd+ir) * wfin(ir)*rinv
     &              + hwfout(ir)
      enddo
c
c     spin contribution
c      (the factors 'sq2' are determined empirically such that
c       galileian invariance is restored)
c
      if(tspinorbit) then
       redmat = redsig(xj,xl,alrpa,alrpa,xja,xla)
c      rfacd  = -t4*redmat   * sq2
c      rfacco = -t4h*redmat  * sq2
cold       rfacd  = -t4*redmat
cold       rfacco = -t4h*redmat
       rfacd  = -(b4+b4p)*redmat  * (1.0D0)  ! * (-2.0) !  sq2
       rfacco = -b4*redmat        * (1.0D0)  ! * (-2.0) !  sq2
       isoco  = ngrid-isoadd
c       write(14,'(a,5i5)') ' iso,la,isoadd=',iso,la,isoadd
       do ir=2,ngrid
         hwfout(ir) = redmat*rspin(isoadd+ir) * wfin(ir)*rinvs(ir)
     &               + hwfout(ir)
c         write(14,'(i5,2g12.4)')
c     &     ir,
c     &  (rfacd*rqcur(isoadd+ir,1)+rfacco*rqcur(isoco+ir,1))/redmat,
c     &     rspin(isoadd+ir)
       enddo
      endif
c
      parll  = (1-mod(ll,2)*2)
      call extrap(hwfout(1),parll)
c
      return
      end
c-----qcurre------------------------------------------------------
c
      subroutine qcurre(resfac)
c
      include 'params.inc'
c
c     \\jq(lam)\\ ,dr.\\jq(lam)\\
c     jq(l)  =sum(lam=l+1,l-1) \\jq(lam)\\   t(lam,m)
c     dqcur is the radial derivative of \\jq\\
c
c       qcur    = output  \\jq(lam)\\
c       dqcur   = output of radial derivative of \\qcur(lam)\\
c
c     the fields for the mean-field hamiltonian are transferred
c     via common /wavfun/.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      data ireadi/0/
c
c------------------------------------------------------------------------------
c
c      write(6,'(a,g12.4,i3)')
c     &  ' QCURRE: resfac,ifresi=',resfac,ifresi
c
c     clear fields
c
      do i=1,ngrid2
        qcur(i,1)  = zero
        qcur(i,2)  = zero
      enddo
c
      if(resfac.eq.zero .or. ifresi.eq.1) then
       do i=1,ngrid2
         dqcur(i) = zero
         rqcur(i,1) = zero
         sigma(i) = zero
       enddo
       return
      endif
c
      parcur = (1-mod(lrpa,2)*2)
c
c     isospin loop
c
      do 10 iso=1,2
      if(iso.eq.1) then
       nmax   = nmaxp
       iradd  = 0
       nstack = 0
      else
       nmax   = nmaxn
       iradd  = ngrid
       nstack = nmaxp*ngrid
      endif
c
c     reset accumulator for spin-current
c
      do ir=1,ngrid
        work3(ir) = zero
      enddo
c
c     'lamp' is lambda-prime
c
      do 100 ilamp=1,2
      lamp   = lrpa-3+2*ilamp
      alamp  = lamp
      if(lamp.lt.0) go to 100
c
c     clear accumulator field
c
      do ir=1,ngrid
        work1(ir) = zero
      enddo
c
c     loop over the a-states
c
      do 20 na=1,nmax
      jsa   = jsp(na)
      la    = jsa/2
      xla   = la
      xja   = (jsa-1)/2+0.5D0
      wgt   = weight(na+nstack/ngrid)
      if(wgt.lt.0.5D0) goto 20
      nadd  = nstack + (na-1)*ngrid
c
c     loop over the j,l to which they couple
c
      do 30 jdif=1,mlrpa
      xj     = xja+jdif-lrpa-1
      if(xj.lt.0D0) goto 30
      ll     = (xj-0.5D0)
      if(mod(lrpa+la+ll,2).ne.0) then
       ll     = (xj+0.5D0)
      endif
      xl     = ll
      if(ll.gt.lrpa+la .or. ll.lt.abs(la-lrpa)) go to 30
      iphase = abs(xja-xj)
      phase  = 1-2*mod(iphase,2)
c
c     first term: derivative acting on 'wf'
c
c      write(91,'(/a,4i5)')
c     &  ' iso,ilamp,na,jdif=',iso,ilamp,na,jdif
      call redcur(xj,xl,alrpa,alamp,xja,xla,facder,facrad)
c      write(91,'(a,2(1pg12.4))') ' for WF: facder,facrad=',facder,facrad
      facder = -facder
      facrad = -facrad
c     facder = -facder*phase
c     facrad = -facrad*phase
      do ir=2,ngrid
        rinv   = rinvs(ir)
        work1(ir) = (facder*wfs(nadd+ir)+facrad*wf(nadd+ir)*rinv)*rinv
     &             * qlrout(nadd+ir,jdif) + work1(ir)
      enddo
c      write(91,'(a,1pg12.4)') ' norm=',overlp(work1,work1,rdum,idum)
c      do ir=1,ngrid
c        write(91,'(1x,f7.1,1pg12.4)') (ir-1)*rstep,work1(ir)
c      enddo
c
c     second term: derivative acting on 'qlrout'
c
      parity = 1-2*mod(la+lrpa,2)
c     call derv15(qlrout(nadd+1,jdif),work2(1),parity,ngrid,rstep)
      call derv17(qlrout(nadd+1,jdif),work2(1),parity,ngrid,rstep)
c
      call redcur(xja,xla,alrpa,alamp,xj,xl,facder,facrad)
c      write(91,'(a,2(1pg12.4))') 
c     &     ' for QLROUT: facder,facrad=',facder,facrad
      facder = facder*phase
      facrad = phase*facrad+facder
c     facder = facder
c     facrad = facrad+facder
      do ir=2,ngrid
        rinv   = rinvs(ir)
        work1(ir) = wf(nadd+ir)*rinv
     &             * (facder*work2(ir)
     &                  +facrad*qlrout(nadd+ir,jdif)*rinv)
     &             + work1(ir)
      enddo
c      write(91,'(a,1pg12.4)') ' norm=',overlp(work1,work1,rdum,idum)
c      do ir=1,ngrid
c        write(91,'(1x,f7.1,1pg12.4)') (ir-1)*rstep,work1(ir)
c      enddo
c
c     accumulate spin-density on 'work3'
c
      if(ilamp.eq.2) then
c       isphas = abs(xja-xj+alrpa)
       isphas = abs(xja-xj)
       sphase = 1-2*mod(isphas,2)
c???                                         sq2 ???
c       redfac = sq2/mlrpa
       redfac = -1.0D0/mlrpa
     &         * ( sphase*redsig(xja,xla,alrpa,alrpa,xj,xl)
     &            -redsig(xj,xl,alrpa,alrpa,xja,xla) )
c      write(6,'(a,g12.4,a,5f4.1)')
c    &  ' testprint redfac=',redfac,' for xj,xl,xja,xla,sphase=',
c    &  xj,xl,xja,xla,sphase
c      write(6,'(3(a,g12.4))')
c    &   ' <ja,la,l,l,j,l>=',redsig(xja,xla,alrpa,alrpa,xj,xl),
c    &   ' <j,l,l,l,ja,la>=',redsig(xj,xl,alrpa,alrpa,xja,xla),
c    &   ' sq2=',sq2
       do ir=2,ngrid
         rinv   = rinvs(ir)
         ira    = nadd+ir
         work3(ir) = redfac*qlrout(ira,jdif)*wf(ira)*rinv + work3(ir)
       enddo
      endif
c
   30 continue
   20 continue
c
c     now store on 'qcur'
c
c     fac    = d4pi/mlrpa
      fac    = one/mlrpa
      do ir=1,ngrid
        qcur(iradd+ir,ilamp) = fac*work1(ir)
      enddo
      call extrap(qcur(iradd+1,ilamp),-parcur)
c
  100 continue
c
c     store spin-current from 'work3' on 'sigma'
c     for given 'iso'
c
      do ir=2,ngrid
        irn    = iradd+ir
        sigma(irn) = work3(ir)
      enddo
      call extrap(sigma(iradd+1),parcur)
c
   10 continue
c
c     compute 'rqcur' from 'qcur' as rotation
c
      call rot(qcur,rqcur,2)
c
c     write(6,'(a,20(/1x,10g12.4))')
c    &    ' testprint of rqcur:',(rqcur(ir,1),ir=1,ngrid)
c     write(6,'(10(/1x,10g12.4))') (rqcur(ir+ngrid,1),ir=1,ngrid)
c     write(14,'(a,20(/1x,5g12.4))')
c    &    ' testprint of qcur(*,1):',(qcur(ir,1),ir=1,ngrid)
c     write(14,'(10(/1x,5g12.4))') (qcur(ir+ngrid,1),ir=1,ngrid)
c     write(14,'(a,20(/1x,5g12.4))')
c    &    ' testprint of qcur(*,2):',(qcur(ir,2),ir=1,ngrid)
c     write(14,'(10(/1x,5g12.4))') (qcur(ir+ngrid,2),ir=1,ngrid)
c     write(14,'(a,20(/1x,5g12.4))')
c    &    ' testprint of sigma:',(sigma(ir),ir=1,ngrid)
c     write(6,'(20(/1x,10g12.4))') (sigma(ir+ngrid),ir=1,ngrid)
c     write(6,'(20(/1x,10g12.4))') (rsigma(ir,2),ir=1,ngrid)
c     write(6,'(20(/1x,10g12.4))') (rsigma(ir+ngrid,2),ir=1,ngrid)
c
      return
      end
c-----qpot  ------------------------------------------------------
c
      subroutine qpot(resfac)
c
      include 'params.inc'
      logical tpotnew
      parameter (tpotnew=.TRUE.)
c
c     the kinetic potential 'rkin' for given response current 'qcur'.
c     the derivative of 'rkin' is accumulated on 'drkin' using
c     'hlinrs' as workspace.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'linres.inc'
      include 'constn.inc'
c
c--------------------------------------------------------------------
c
      if(tpotnew .and. ifgsc.eq.7)
     &   stop " TDHF requires TPOTNEW=.FALSE."
      parity = 1-2*mod(lrpa,2)
c
       do i=1,ngrid2
         drkin(i) = zero
         rkin(i,1) = zero
         rkin(i,2) = zero
         rspin(i) = zero
       enddo
      if(resfac.eq.zero .or. ifsepi.eq.1 .or. oddoff  
     &   .or. tsemicl) then
        drkin(1:ngrid2) = 0D0
        rkin(1:ngrid2,1:2) = 0D0
        rspin(1:ngrid2) = 0D0
        return
      endif
c
c     overrule rkin for testing separable interaction
c
      if(ifsepi.eq.1) return
      if(tpotnew) then
        call qpotnew()
        return
      endif
c
      t12co2 = 0.5D0*(t1+t1*0.5D0*x1+t2+t2*0.5D0*x2)
      t12q2  = t12co2-0.5D0*(0.5D0*t1+t1*x1-0.5D0*t2-t2*x2)
      if(tnob1) then
        t12co2 = 0D0
        t12q2  = 0D0
      endif
ctest      write(14,'(a,4(1pg12.4))') ' QPOT: t1...x2=',t1,t2,x1,x2
      do 10 ilamp=1,2
      do ir=1,ngrid
        rkin(ir,ilamp)       = t12q2*qcur(ir,ilamp)
     &                        + t12co2*qcur(ngrid+ir,ilamp)
        rkin(ngrid+ir,ilamp) = t12q2*qcur(ngrid+ir,ilamp)
     &                        + t12co2*qcur(ir,ilamp)
      enddo
c
c     contribution from extended parametrisation
c
      if(isky.eq.10) then
       fac56p = 2.0*t5p6
       fac56m = fac56p-2.0*t5m6
       do i=1,ngrid
         ic     = i+ngrid
         rs     = rho(i)+rho(ic)
         rkin(i,ilamp) = rs*(fac56m*qcur(i,ilamp)+fac56p*qcur(ic,ilamp))
     &                   + rkin(i,ilamp)
         rkin(ic,ilamp)= rs*(fac56m*qcur(ic,ilamp)+fac56p*qcur(i,ilamp))
     &                 + rkin(ic,ilamp)
       enddo
      endif
c
c     @@@ contribution from surface sensitive effective mass                !new
c                                                                             .
      if(t9.ne.zero) then                                                    .
         do ir=1,ngrid                                                     .
            irc    = ngrid + ir                                               .
            t9add  = -2.0 * t9 *                                              .
     $               (d2r(ir) + d2r(irc))                                     .
     $              /(rho(ir) + rho(irc) + xsmall)                            .
            rkin(ir,ilamp)       = t9add * qcur(ir,ilamp)                     .
     $                           + rkin(ir,ilamp)                             .
            rkin(ngrid+ir,ilamp) = t9add * qcur(irc,ilamp)                    .
     $                           + rkin(irc,ilamp)                            .
         enddo
      endif                                                                !new
c
c     contribution from spin current
c
      if(tspinorbit) then
c       write(6,'(a,50(/1x,5g12.4))')
c     &  ' sigma protons:',(sigma(ir),ir=1,ngrid)
c       write(6,'(a,50(/1x,5g12.4))')
c     &  ' sigma neutrons:',(sigma(ir+ngrid),ir=1,ngrid)
       call rotsig(sigma(1),aux1,ilamp)
       call rotsig(sigma(ngrid+1),aux2,ilamp)
       b4sum = b4+b4p
       do ir=1,ngrid
         irn    = ngrid+ir
         rkin(ir,ilamp)  = rkin(ir,ilamp)
     &                    - b4sum*aux1(ir) - b4*aux2(ir)
         rkin(irn,ilamp) = rkin(irn,ilamp)
     &                    - b4sum*aux2(ir) - b4*aux1(ir)
c         write(14,'(2i5,3g12.4)')   ir,ilamp,
c     &   -0.5D0*(b4sum*aux1(ir)+b4*aux2(ir)-b4sum*aux2(ir)-b4*aux1(ir))
       enddo
      endif
c
   10 continue
c
c     the derivative kinetic potential
c
      call div(rkin,drkin,2)
c
c     the spin potential
c      
      if(tspinorbit) then
c        call rot(qcur,rqcur,2)
        b4sum = b4+b4p
        do ir=1,ngrid
          irn    = ngrid+ir
          rspin(ir)  = -b4sum*rqcur(ir,1)-b4*rqcur(irn,1)
          rspin(irn) = -b4sum*rqcur(irn,1)-b4*rqcur(ir,1)
c         write(14,'(i5,3g12.4)')   ir,(rspin(ir)-rspin(irn))*0.5D0
        enddo
c        write(14,'(a,g12.4)') ' rspin done. b4sum =',b4sum
      endif
c
      return
      end
c-----commhp--------------------------------------------------------------------
c
      subroutine commhp(ovpq,ovpq0,ovqq,ovq0q0,resfac,cmfac,ifsav0)
c
      include 'params.inc'
c
c     computes q = [h,p] for given p
c     with p and q given by its action on wave functions:
c         p |wf> <--> 'plrout'
c         q |wf> <--> 'qlrout'
c
c     the response wave functions are handled via common /linres/.
c     they mean:
c         'plrout' = given p = input to routine
c         'qlrout' = resulting [h,p] = output of routine
c
c     during evaluation of [h,p] some overlaps are accumulated:
c      ovpq    = <wf|p [h,p]|wf>
c      ovpq0   = <wf|p [h0,p]|wf>
c      ovqq    = <wf|q q|wf> = <wf|[p,h] [h,p]|wf>
c      ovq0q0  = <wf|[p,h0] [h0,p]|wf>
c     where h is the full hamiltonian and h0 the mean-field part
c     action of the p(lm) on-the s.p wave function
c
c     further list parameters are:
c      resfac  = weight factor for the residual interaction;
c                should be 1.0
c      ifsav0  = switch for saving the h0-commutator on 'savlr'
c
c     result is restored on 'qlrout' overwriting old value
c
c     where:
c
c       qlrout  = output qblraj field ,
c
c     the fields for the mean-field hamiltonian are transferred
c     via common /wavfun/.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c----------------------------------------------------------------------
c
c     prepare the response current 'qcur' and its derivative 'dqcur'
c
      call pdensi(resfac)
      if(.not.tyukaw .and. ayuk.ne.zero) then
c       write(6,'(a)') ' HCOMMP: damping switch '
       dampay = ayuk*ayuk
c       do 915 idamp=1,8
       call t0inv(prho(1),prho(1),work1,dampay,lrpa)
       call t0inv(prho(ngrid+1),prho(ngrid+1),work1,dampay,lrpa)
       call t0inv(ptau(1),ptau(1),work1,dampay,lrpa)
       call t0inv(ptau(ngrid+1),ptau(ngrid+1),work1,dampay,lrpa)
       if(imfold.eq.1) then
        call t0inv(ptau(1),ptau(1),work1,dampay,lrpa)
        call t0inv(ptau(ngrid+1),ptau(ngrid+1),work1,dampay,lrpa)
       endif
       if(tspinorbit) then
        if(lrpa.gt.0) then
         call t0inv(psj(1,1),psj(1,1),work1,dampay,lrpa-1)
         call t0inv(psj(ngrid+1,1),psj(ngrid+1,1),work1,dampay,lrpa-1)
         if(imfold.eq.1) then
          call t0inv(psj(1,1),psj(1,1),work1,dampay,lrpa-1)
          call t0inv(psj(ngrid+1,1),psj(ngrid+1,1),work1,dampay,lrpa-1)
         endif
        endif
        call t0inv(psj(1,2),psj(1,2),work1,dampay,lrpa+1)
        call t0inv(psj(ngrid+1,2),psj(ngrid+1,2),work1,dampay,lrpa+1)
        if(imfold.eq.1) then
         call t0inv(psj(1,2),psj(1,2),work1,dampay,lrpa+1)
         call t0inv(psj(ngrid+1,2),psj(ngrid+1,2),work1,dampay,lrpa+1)
        endif
        call div(psj,psjs,2)
       endif
c  915  continue
      elseif(resdmp.ne.zero) then
c       write(6,'(a,2(1pg13.5))') 
c     &  ' HCOMMP: filter switch, resdmp,zero= ',resdmp,zero
       call filter(prho(1),prho(1),lrpa)
cc       write(6,'(a)') ' repeat filter:'
cc       call filter(prho(1),prho(1),lrpa)
       call filter(prho(ngrid+1),prho(ngrid+1),lrpa)
       call filter(ptau(1),ptau(1),lrpa)
       call filter(ptau(ngrid+1),ptau(ngrid+1),lrpa)
       if(tspinorbit) then
        if(lrpa.gt.0) then
         call filter(psj(1,1),psj(1,1),lrpa-1)
         call filter(psj(ngrid+1,1),psj(ngrid+1,1),lrpa-1)
        endif
        call filter(psj(1,2),psj(1,2),lrpa+1)
        call filter(psj(ngrid+1,2),psj(ngrid+1,2),lrpa+1)
        call div(psj,psjs,2)
       endif
      endif
c
c     prepare the kinetic response potential for given current
c
      call ppot(resfac,cmfac)
      if(.not.tyukaw .and. ayuk.ne.zero) then
c       write(6,'(a)') ' HCOMMP: damping switch '
       call t0inv(rpot(1),rpot(1),work1,dampay,lrpa)
       call t0inv(rpot(ngrid+1),rpot(ngrid+1),work1,dampay,lrpa)
       call t0inv(rbmass(1),rbmass(1),work1,dampay,lrpa)
       call t0inv(rbmass(ngrid+1),rbmass(ngrid+1),work1,dampay,lrpa)
       if(imfold.eq.1) then
        call t0inv(rbmass(1),rbmass(1),work1,dampay,lrpa)
        call t0inv(rbmass(ngrid+1),rbmass(ngrid+1),work1,dampay,lrpa)
       endif
       call grad(rbmass,drbmas,2)
       if(tspinorbit) then
        if(lrpa.gt.0) then
         call t0inv(rpotls(1,1),rpotls(1,1),work1,dampay,lrpa-1)
         call t0inv(rpotls(ngrid+1,1),rpotls(ngrid+1,1),work1,
     &              dampay,lrpa-1)
        endif
        call t0inv(rpotls(1,2),rpotls(1,2),work1,dampay,lrpa+1)
        call t0inv(rpotls(ngrid+1,2),rpotls(ngrid+1,2),work1,
     &             dampay,lrpa+1)
       endif
      endif
c
c     prepare action of the p-response hamiltonian on the wave
c     function. result on 'hlinrs'.
c
      ovpq   = zero
      ovpq0  = zero
      ovqq   = zero
      ovq0q0 = zero
c
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     action of linear response hamiltonian
c
      if(resfac.ne.zero)
     &  call hrespp(iso,la,xja,nadd,xj,xl,jdif,wf(nadd+1),wfs(nadd+1))
c
c     commutator of hf-hamiltonian on p-wavefunction
c
      call multr(plrout(nadd+1,jdif),work1(1),rstep,ngrid)
      call hpsi(work1(1),work2(1),aux1(1),aux2(1),1,js,iso,subfac)
c
c     now  store plrout.
c
      if(resfac.eq.0D0) then
        do ir=2,ngrid
          work2(ir) = work2(ir)*rinvs(ir)
          qlrout(nadd+ir,jdif) = work2(ir) 
        enddo
      else
        do ir=2,ngrid
          work2(ir) = work2(ir)*rinvs(ir)
          qlrout(nadd+ir,jdif) = work2(ir) + resfac * hlinrs(ir)
        enddo
      endif
      parll  = (1-mod(ll,2)*2)
      if(ifsav0.eq.1) then
       do ir=2,ngrid
         savlr(nadd+ir,jdif) = work2(ir)
       enddo
       call extrap(savlr(nadd+1,jdif),parll)
       call projec(savlr(nadd+1,jdif),js,nstack,nmax)
      endif
      call extrap(qlrout(nadd+1,jdif),parll)
      call extrap(work2(1),parll)
      call projec(qlrout(nadd+1,jdif),js,nstack,nmax)
      call projec(work2(1),js,nstack,nmax)
c
      ovpq   = overlp(plrout(nadd+1,jdif),qlrout(nadd+1,jdif),
     &                rstep,ngrid)  +  ovpq
      ovpq0  = overlp(plrout(nadd+1,jdif),work2(1),
     &                rstep,ngrid) + ovpq0
      ovqq   = overlp(qlrout(nadd+1,jdif),qlrout(nadd+1,jdif),
     &                rstep,ngrid)  +  ovqq
      ovq0q0 = overlp(work2(1),work2(1),rstep,ngrid) + ovq0q0
c
c     commutator of hf-hamiltonian on imaginary p-wavefunction
c     and subtraction of 'enext'-i*'enwid'
c     in case of complex response
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
      return
      end
c-----hrespp-----------------------------------------------------------------
c
      subroutine hrespp(iso,la,xja,nadd,xj,xl,jdif,wfin,wfsin)
c
      include 'params.inc'
c
c     action of linear response hamiltonian h_p (for p-like perturbations)
c     on the single particle h.f. state from 'wf' in common /wavfun/
c     according to the quantum numbers
c      iso    = isospin
c      la     = orbital angular momentum of 'wf'
c      xja    = total angular momentum of 'wf'
c      xl     = orbital angular momentum of result on'hlinrs'
c      xlj    = total angular momentum of result on'hlinrs'
c      jdif   = difference quantum number between angular momenta
c      nadd   = pointer to the wave function in 'wf'
c
c     the response potentials 'rpot', 'rpotls', 'rbmass' and 'drbmas'
c     have been computed in 'ppot' and are transferred via
c     common /linres/.
c
c     the result is:
c      hlinrs  = output via common /linres/
c
c     the fields for the mean-field hamiltonian are transferred
c     via common /wavfun/.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      dimension wfin(kg1),wfsin(kg1)
      dimension facr(2)
c      dimension help(kg1),save1(kg1)
c      dimension aux4(kg1)
c
c--------------------------------------------------------------------
c
      xla    = la
      ll     = xl
      isoadd = (iso-1)*ngrid
c
c     clear accumulator field
c
      do ir=1,ngrid
        hlinrs(ir) = zero
      enddo
c
c     first: the complicated term with the derivative of the inverse
c            mass
c
c     loop over  lambda' = l+-1
c
      do 20 ilamp=1,2
      lamp   = lrpa-3+2*ilamp
      if(lamp.ge.0) then
       alamp  = lamp
c
c      derivative on wave function, subcase: lambda = la-1
c
       call redcur(xj,xl,alrpa,alamp,xja,xla,facder,facrad)
       do ir=2,ngrid
         rinv   = rinvs(ir)
         hlinrs(ir) = -drbmas(isoadd+ir,ilamp)
     &               * (facder*wfsin(ir)+facrad*wfin(ir)*rinv)*rinv
     &               + hlinrs(ir)
       enddo
c
      endif
   20 continue
c
c     second: simple mass term and potential terms
c
      parity = 1-2*mod(la,2)
      call derv25(wfin(1),aux2(1),-parity,ngrid,rstep)
      redmat = redume(xj,alrpa,xja)
      amfac  = la+la*la
      do 210 ir=2,ngrid
      rinv   = rinvs(ir)
      irs    = isoadd+ir
      hlinrs(ir) = redmat*rinv
     &            * (-rbmass(irs)*(aux2(ir)-amfac*(rinv*rinv)*wfin(ir))
     &               +rpot(irs)*wfin(ir) )
     &            + hlinrs(ir)
  210 continue
c
c     third: the spin-orbit term
c
      if(tspinorbit) then
       do 310 ilamp=1,2
       lamp   = lrpa-3+2*ilamp
       if(lamp.ge.0) then
        alamp  = lamp
        call redspo(xj,xl,alrpa,alamp,xja,xla,facder,facrad)
        facder = -sq2*facder
        facrad = -sq2*facrad
        facr(ilamp) = facrad
        do ir=2,ngrid
          rinv   = rinvs(ir)
          hlinrs(ir) = rpotls(isoadd+ir,ilamp)
     &                * (facder*wfsin(ir)+facrad*wfin(ir)*rinv)*rinv
     &                + hlinrs(ir)
        enddo
       endif
  310  continue
c
c       test section
c
c        iftest = 0
c        if(alrpa.eq.1 .and. iso.eq.1 .and. iftest.eq.1) then
c         do 840 ir=2,ngrid
c         work1(ir) = zero
c  840    work2(ir) = zero
c         do 850 ilamp=1,2
c         lamp   = lrpa-3+2*ilamp
c         if(lamp.ge.0) then
c          alamp = lamp
c          call redspo(xj,xl,alrpa,alamp,xja,xla,facder,facrad)
c          facder = -sq2*facder
c          facrad = -sq2*facrad
c          do 860 ir=2,ngrid
c          work1(ir) = rpotls(isoadd+ir,ilamp)
c     &               * (facder*wfsin(ir))*rinvs(ir)
c     &               + work1(ir)
c          work2(ir) = rpotls(isoadd+ir,ilamp)
c     &               * (facrad*wfin(ir)*rinvs(ir))*rinvs(ir)
c     &               + work2(ir)
c  860     continue
c         endif
c  850    continue
c         write(6,'(a,6f4.1/a)')
c     &    ' hrespp: j,l,l,lambda,ja,la=',xj,xl,alrpa,alamp,xja,xla,
c     &   ' hlinrs    model    hl.der  mod.der  hl.rad   mod.r1  mod.r2'
cc
cc        test spin-orbit response
cc
c         isoadd = (iso-1)*ngrid
c         iscadd = (2-iso)*ngrid
c         do 810 ir=1,ngrid
c  810    aux3(ir) =-t4h*(d1r(ir+isoadd)*2.0+d1r(iscadd+ir))*rinvs(ir)
c     &              * 41.505
c         call redspo(xj,xl,one,zero,xja,xla,facder,facrad)
c         call derv15(aux3(1),help(1),one,ngrid,rstep)
c         redmat = redume(xj,alrpa,xja)*(xja+xja*xja-xla-xla*xla-0.75)
c         do 820 ir=1,ngrid
c         rinv   = rinvs(ir)
c         aux2(ir) = 4.24264*aux3(ir)
c     &           * (facrad*wfin(ir)*rinv)*rinv
c         aux1(ir) = 4.24264*aux3(ir)
c     &           * (facder*wfsin(ir))*rinv
c  820    aux3(ir) = (help(ir)*redmat*wfin(ir)*rinv )
c         do 830 ir=1,ngrid
c         write(6,'(7f9.3)') work1(ir)+work2(ir),
c     &                      aux1(ir)+aux2(ir)+aux3(ir),
c     &                      work1(ir),aux1(ir),
c     &                      work2(ir),aux2(ir),aux3(ir)
c         work1(ir) = work1(ir)+work2(ir)
c  830    continue
cc
cc        hlinrs by commutator with h0
cc
c         iphas  = xl+xja-half
c         redmat = (1-mod(iphas,2)*2)*sixjot(xl,xj,half,xja,xla,one)
c     &           *sqrt((xj+xj+one)*(xja+xja+one))
c         if(xl.eq.xla+one) then
c          facder = -sqrt(xla+one)*redmat
c          facrad = -(xla+one)*facder
c         elseif(xl.eq.xla-one) then
c          facder = sqrt(xla)*redmat
c          facrad = xla*facder
c         else
c          facder = zero
c          facrad = zero
c         endif
c         parity = (mod(la,2)*2-1)
c         call derv17(wfin(1),aux1(1),parity,ngrid,rstep)
c         do 710 ir=1,ngrid
c  710    aux1(ir) = facder*aux1(ir)+facrad*wfin(ir)*rinvs(ir)
c         if(xj.gt.xl) then
c          js     = xj+xj
c         else
c          js     = xj+xj+one
c         endif
c         call hlspsi(aux1(1),aux2(1),js,iso)
c         if(xja.gt.xla) then
c          jsa    = xja+xja
c         else
c          jsa    = xja+xja+one
c         endif
c         call hlspsi(wfin(1),aux1(1),jsa,iso)
c         call derv17(aux1(1),aux3(1),parity,ngrid,rstep)
c         do 720 ir=2,ngrid
c         aux1(ir) = (facder*aux3(ir)+facrad*aux1(ir)*rinvs(ir))
c         aux2(ir) = -41.505*(aux2(ir)-aux1(ir))*rinvs(ir)
c  720    aux1(ir) = aux2(ir)/work1(ir)
c         write(6,'(a)') ' hlinrs versus [hls,wf]  quotient'
c         do 730 ir=2,ngrid
c  730   write(6,'(1x,3f9.3)') work1(ir),aux2(ir),aux1(ir)
c        endif
c
c       end test section
c
      endif
c
      parll  = (1-mod(ll,2)*2)
      call extrap(hlinrs,parll)
c
      return
      end
c-----pdensi------------------------------------------------------
c
      subroutine pdensi(resfac)
c
      include 'params.inc'
c
c     computes p-response denities 'prho','ptau','psj'.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'result.inc'
      include 'constn.inc'
c
c--------------------------------------------------------------------
c
c     clear fields
c
      do i=1,ngrid2
        prho(i) = zero
        ptau(i) = zero
        psj(i,1) = zero
        psj(i,2) = zero
      enddo
      if(resfac.eq.zero .or. ifresi.eq.1) then
       do i=1,ngrid2
         psjs(i) = zero
       enddo
       return
      endif
c
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     accumulate the densities 'pden', 'ptau' and 'psjs'
c
      parity = 1-2*mod(la+lrpa,2)
      call derv17(plrout(nadd+1,jdif),work1(1),parity,ngrid,rstep)
      iph    = xj-xja
      redfac = 2.0
     &        * redume(xja,alrpa,xj)*(1-2*mod(abs(iph),2))/mlrpa
      taufac = 0.5D0*(xla*xla+xla+xl*xl+xl-alrpa*alrpa-alrpa)
      do 40 ir=2,ngrid
      rinv   = rinvs(ir)
      ira    = nadd+ir
      iro    = isoadd+ir
      prho(iro) = redfac*plrout(ira,jdif)*wf(ira)*rinv + prho(iro)
      ptau(iro) = redfac
     &           * (work1(ir)*(wfs(ira)-wf(ira)*rinv)
     &              +taufac*plrout(ira,jdif)*wf(ira)*rinv*rinv)*rinv
     &           + ptau(iro)
   40 continue
c
c     test plrout
c
c      iftest = 0
c      if(lrpa.eq.1 .and. isorpa.eq.1 .and. jsa.eq.1 .and. iso.eq.1
c     &   .and. iftest.eq.1) then
c       redmat = 2.0 * redume(xj,alrpa,xja)/mlrpa
c       write(6,'(a,i2,a,4f4.1/2(a,f9.4,a,f9.4)/a)')
c     &     ' jdif=',jdif,' j,l,ja,la=',xj,xl,xja,xla,
c     &     '  redmat=',redmat,'  redfac=',redfac
c       work2(1) = zero
c       do 910 ir=2,ngrid
c       compar = -41.505*(wfs(nadd+ir)-wf(nadd+ir)*rinvs(ir))*rinvs(ir)
c       work2(ir) = plrout(isoadd+ir,jdif)/compar
c  910  continue
c       write(6,'(1x,5g12.4)') (work2(ir),ir=1,ngrid)
c      endif
c
c     accumulate the spin-orbit density 'psj'
c
      do 100 ilamp=1,2
      lamp   = lrpa-3+ilamp+ilamp
      if(lamp.ge.0) then
c
c      set plrout to wf/r for testing psjs in case of l=0
c
c       iftest = 0
c       if(mlrpa.eq.1 .and. isorpa.eq.1 .and. iftest.eq.1) then
c        do 920 ir=2,ngrid
c  920   plrout(nadd+ir,jdif) = wf(nadd+ir)*rinvs(ir)
c        parplr = 1-2*mod(la,2)
c        call extrap(plrout(nadd+1,jdif1),parplr)
c       endif
c
       testfac =  sq2
       alamp  = lamp
       call redspo(xj,xl,alrpa,alamp,xja,xla,facder,facrad)
       facder = facder/mlrpa
       facrad = facrad/mlrpa
       facder = -facder  * testfac
       facrad = -facrad  * testfac
       do ir=1,ngrid
         rinv   = rinvs(ir)
         ira    = nadd+ir
         iro    = isoadd+ir
         psj(iro,ilamp) = plrout(ira,jdif)
     &                   * (facder*wfs(ira) + facrad*wf(ira)*rinv)*rinv
     &                   + psj(iro,ilamp)
       enddo
       iphase = abs(xja-xj)
       phase  = 1-2*mod(iphase,2)
       call redspo(xja,xla,alrpa,alamp,xj,xl,facder,facrad)
       facder = -phase*facder/mlrpa   * testfac
       facrad = -phase*facrad/mlrpa   * testfac  + facder
       facder =  facder
       facrad =  facrad
c      write(6,'(44x,2(a,f9.4))')
c    &    '   facder=',facder,'  facrad=',facrad
       do ir=1,ngrid
         rinv   = rinvs(ir)
         ira    = nadd+ir
         iro    = isoadd+ir
         psj(iro,ilamp) =(facder*work1(ir)+facrad*plrout(ira,jdif)*rinv)
     &                   * wf(ira)*rinv       + psj(iro,ilamp)
       enddo
      endif
  100 continue
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
      parrho = (1-mod(lrpa,2)*2)
      call extrap(prho(1),parrho)
      call extrap(ptau(1),parrho)
      call extrap(prho(ngrid+1),parrho)
      call extrap(ptau(ngrid+1),parrho)
      call extrap(psj(1,1),-parrho)
      call extrap(psj(ngrid+1,1),-parrho)
      call extrap(psj(1,2),-parrho)
      call extrap(psj(ngrid+1,2),-parrho)
c
c     compute 'psjs' from 'psj' as divergence
c
      call div(psj,psjs,2)
c
c      iftest = 0
c      if(isorpa.eq.1 .and. iftest.eq.1) then
c        write(6,'(a,i3,10(/1x,5g12.4))')
c     &   ' psjs: for isorpa=',isorpa,
c     &   (psjs(ir),ir=1,ngrid)
c        if(mlrpa.eq.1) write(6,'(a,10(/1x,5g12.4))')
c     &   ' sjs:',(sjs(ir),ir=1,ngrid)
c      endif
c      if(lrpa.eq.1 .and. isorpa.eq.1 .and. iftest.eq.1) then
c        call derv17(sjs(1),work1(1),one,ngrid,rstep)
c        write(6,'(a,10(/1x,5g12.4))') ' test psjs by explicit deriv.:',
c     &   (-41.505*work1(ir),ir=1,ngrid)
c      endif
c
      return
      end
c-----ppot  ------------------------------------------------------
c
      subroutine ppot(resfac,cmfac)

      include 'params.inc'
      logical tpotnew
      parameter (tpotnew=.TRUE.)
c
c     the response potentials 'rpot', 'rpotls' and 'rbmass'.
c
      include 'wavfun.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'potinv.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      data xminlw/1.0D-20/
c
c--------------------------------------------------------------------
c
      if(tpotnew .and. tsemicl) 
     &   STOP "semi-classical tau needs old response potentials"
      if(tpotnew .and. ifgsc.eq.7)
     &   stop " TDHF requires TPOTNEW=.FALSE."
      parity = 1-2*mod(lrpa,2)
c
c      if(resfac.eq.zero .or. ifsepi.eq.1 .or. evenoff) then
       do  i=1,ngrid+ngrid
         rpot(i) = zero
         rpotls(i,1) = zero
         rpotls(i,2) = zero
         drbmas(i,1) = zero
         drbmas(i,2) = zero
         rbmass(i) = zero
       enddo
       if(resfac.eq.zero .or. evenoff) return
c      endif
c
c     prepare strength of additional cm-force
c
      if(lrpa.eq.1 .and. cmfac.ne.zero) then
       raver  = zero
       do ir=2,ngrid
         r      = (ir-1)*rstep
         raver  = r*r*r*(prho(ir)+prho(ir+ngrid)) + raver
       enddo
       raver  = frpirs*raver
       cmrfac = cmfac*raver
c       write(6,'(a,g12.4)') ' cmrfac=',cmrfac
      endif
c
c     overrule rpot in case of the separable interaction
c     of the type q_lm * q_lm with radial form r^l.
c
      if(ifsepi.eq.1) then
c
c     expectation value of 'qrad'-operator
c
       raver  = zero
       do ir=2,ngrid
         r      = (ir-1)*rstep
         r      = r**(2+lrpa)
         raver  = r*(prho(ir)+prho(ir+ngrid)) + raver
       enddo
       raver  = frpirs*raver
       qsepfc = raver
c
c      setup response potential from separable r^lrpa
c
       rpot(1)       = zero
       rpot(1+ngrid) = zero
       do ir=2,ngrid
         r      = (ir-1)*rstep
         r      = r**lrpa*qsepfc
         rpot(ir)       = r
         rpot(ir+ngrid) = r
       enddo
c
       if(iwr.ge.00)
     &  write(6,'(a,2g12.3,a,20(/1x,5g12.3))')
     &   ' resfac,qsepfc=',resfac,qsepfc,' and the rpot:',
     &   (rpot(ir),ir=1,ngrid2)
c
       return
      endif
c
c
c     optionally additional cm-potential
c
      if(lrpa.eq.1 .and. cmfac.ne.zero) then
       do ir=2,ngrid
         r      = (ir-1)*rstep
         rpot(ir)       = cmrfac*r   + rpot(ir)
         rpot(ir+ngrid) = cmrfac*r   + rpot(ir+ngrid)
       enddo
      endif
c
c     coulomb  contribution
c
      if(acoul.ne.zero) then
         call coulml(prho(1),aux1,lrpa)
         do ir=1,ngrid
           rpot(ir) = aux1(ir)+rpot(ir)
         enddo
         if(ifex.eq.1) then
            excfac = -1.41803*abs(acoul)/3.0
            twdthr = -0.66666666667
         endif
      endif
c
c     new version of Skyrme residual interactionc
c
      if(tpotnew) then
!        call testgrad()
        call ppotnew()
        return
      endif
c
c     terms depending on 'rho'
c
      t0x    = t0+0.5D0*t0*x0
      t0y    = 0.5D0*t0+t0*x0
      t3x    = (t3+t3*0.5D0*x3)*(power+one)*(power+2.0)/12.0
      t3y1   = (0.5D0*t3+t3*x3)/6.0
      t3y2   = 0.5D0*t3y1*power*(power-one)
      b5     = (t5*(one+x5*0.5D0)+t6*(one+x6*0.5D0))/24.0
      do iso=1,2
        isoadd = (iso-1)*ngrid
        isoco  = (2-iso)*ngrid
        do ir=1,ngrid
          irq    = isoadd+ir
          irco   = isoco +ir
          rq     = rho(irq)
          rc     = rho(irco)
          rs     = rq+rc
          if(rs.gt.1.0D-99) then
            rpowm2 = rs**(power-2.0)
          else
            rpowm2 = 0D0
          endif
          prq    = prho(irq)
          prc    = prho(irco)
          prs    = prq+prc
          rpot(irq) = t0x*prs-t0y*prq
     &           +t3x*rpowm2*(rs*rs)*prs
     &           -t3y1*rpowm2*rs*(power*(rq*prs+rq*prq+rc*prc)+rs*prq)  ! ??
     &           -t3y2*rpowm2*prs*(rq*rq+rc*rc)
     &           + rpot(irq)
          WRITE(9,'(2i5,1pg13.5)') iso,ir,t0x*prs-t0y*prq
     &           +t3x*rpowm2*(rs*rs)*prs
     &           -t3y1*rpowm2*rs*(power*(rq*prs+rq*prq+rc*prc)+rs*prq)  ! ??
     &           -t3y2*rpowm2*prs*(rq*rq+rc*rc)
     &           + (2-iso)*excfac*rq**twdthr*prq
        enddo
      enddo
c
c     terms depending on 'tau'
c
ctest      write(14,'(a,4(1pg12.4))') ' PPOT: t1...x2=',t1,t2,x1,x2
      t12co  = 0.25*(t1+t1*0.5D0*x1+t2+t2*0.5D0*x2)
      t12q   = t12co-0.25*(0.5D0*t1+t1*x1-0.5D0*t2-t2*x2)
      if(tnob1) then
        t12co = 0D0
        t12q  = 0D0
      endif
      if(tsemicl) then
        t12cod = 0.5D0*t12co
        t12qd  = 0.5D0*t12q
        do ir=1,ngrid
          irc    = ngrid+ir
          pr     = prho(ir)
          prc    = prho(irc)
          pt     = (3.*pi*pi*rho(ir))**(0.66666667D0)*0.6D0*pr
          ptc    = (3.*pi*pi*rho(irc))**(0.66666667D0)*0.6D0*prc
          rpot(ir)  = 2D0*(t12q*pt+t12co*ptc) + rpot(ir)
          rpot(irc) = 2D0*(t12q*ptc+t12co*pt) + rpot(irc)
          rbmass(ir)  = 0D0
          rbmass(irc) = 0D0
        enddo
      else
        do ir=1,ngrid
          irc    = ngrid+ir
          pt     = ptau(ir)
          ptc    = ptau(irc)
          rpot(ir)  = t12q*pt+t12co*ptc + rpot(ir)
          rpot(irc) = t12q*ptc+t12co*pt + rpot(irc)
          pr     = prho(ir)
          prc    = prho(irc)
          rbmass(ir)  = t12q*pr+t12co*prc
          rbmass(irc) = t12q*prc+t12co*pr
          if(b5.ne.zero) then
           add       = 2.0*b5*((rho(ir)+rho(irc))*(pt+ptc)
     &                         +(tau(ir)+tau(irc))*(pr+prc))
           rpot(ir)  = add + rpot(ir)
           rpot(irc) = add + rpot(irc)
           add       = 2.0*b5*(rho(ir)+rho(irc))*(pr+prc)
           rbmass(ir)  = add + rbmass(ir)
           rbmass(irc) = add + rbmass(irc)
          endif
        enddo
      endif
c
c      coulomb exchange term in local density approximation
c
      if(acoul.ne.zero .and. ifex.eq.1) then
        do ir=1,ngrid
          rp     = rho(ir)
          if(rp.gt.xminlw)
     &        rpot(ir) = excfac*rp**twdthr*prho(ir) + rpot(ir)
        enddo
      endif
c
c     terms depending on laplace->'rho'
c
      if(imfold.eq.1) then
       dampay = ayuk*ayuk
       call t0inv(prho(1),prho(1),work1,dampay,lrpa)
       call t0inv(prho(ngrid+1),prho(ngrid+1),work1,dampay,lrpa)
      endif
      t12co  = -0.125*(3.0*(t1+t1*0.5D0*x1)-t2-t2*0.5D0*x2)
      t12q   = t12co+0.125*(3.0*(0.5D0*t1+t1*x1)+0.5D0*t2+t2*x2)
      if(tnosurf) then
        t12co  = 0D0
        t12q   = 0D0
      endif
      call derv15(prho(1),dprho(1),parity,ngrid,rstep)
      call derv15(prho(ngrid+1),dprho(ngrid+1),parity,ngrid,rstep)
      call derv25(prho(1),d2prho(1),parity,ngrid,rstep)
      call derv25(prho(ngrid+1),d2prho(ngrid+1),parity,ngrid,rstep)
      afac   = alrpa+alrpa*alrpa
      do 30 ir=2,ngrid
      irc    = ngrid+ir
      rinv   = rinvs(ir)
      dpr     = d2prho(ir)+(2.0*dprho(ir)-afac*prho(ir)*rinv)*rinv
      dprc    = d2prho(irc)+(2.0*dprho(irc)-afac*prho(irc)*rinv)*rinv
      rpot(ir)  = t12q*dpr+t12co*dprc + rpot(ir)
      rpot(irc) = t12q*dprc+t12co*dpr + rpot(irc)
      if(tsemicl) then
        rpot(ir)  = t12qd*dpr+t12cod*dprc + rpot(ir)
        rpot(irc) = t12qd*dprc+t12cod*dpr + rpot(irc)
      endif
c        @@@ added to include surface sens. m*:
      if(t9.ne.zero) then
         t9dpra      = t9*(dpr + dprc)
         rbmass(ir)  = -t9dpra + rbmass(ir)
         rbmass(irc) = -t9dpra + rbmass(irc)
      endif
   30 continue
c
c     terms depending on laplace->'tau'                                     !new
c                                                                              .
      if (t9.ne.zero) then                                                    .
         call derv15(ptau(1),dptau(1),parity,ngrid,rstep)                     .
         call derv15(ptau(ngrid+1),dptau(ngrid+1),parity,ngrid,rstep)         .
         call derv25(ptau(1),d2ptau(1),parity,ngrid,rstep)                    .
         call derv25(ptau(ngrid+1),d2ptau(ngrid+1),parity,ngrid,rstep)        .
         if(imfold.eq.1) then
          call t0inv(dptau(1),dptau(1),work1,dampay,lrpa)
          call t0inv(dptau(ngrid+1),dptau(ngrid+1),work1,dampay,lrpa)
          call t0inv(d2ptau(1),d2ptau(1),work1,dampay,lrpa)
          call t0inv(d2ptau(ngrid+1),d2ptau(ngrid+1),work1,dampay,lrpa)
         endif
         do ir=2,ngrid                                                     .
            afac   = alrpa+alrpa*alrpa                                        .
            irc    = ngrid+ir                                                 .
            rinv   = rinvs(ir)                                                .
            dpt    = d2ptau(ir)+(2.0*dptau(ir)-afac*ptau(ir)*rinv)*rinv       .
            dptc   = d2ptau(irc)                                              .
     $               +(2.0*dptau(irc)-afac*ptau(irc)*rinv)*rinv               .
            t9dpta  = t9*(dpt + dptc)                                         .
            rpot(ir)  = - t9dpta + rpot(ir)                                   .
            rpot(irc) = - t9dpta + rpot(irc)                                  .
         enddo
      end if                                                              !new
c      write(14,'(a)') ' response potentials and masses ',
c     &   '  r   rpot_p   rpot_n   bmass_p    bmass_n   '
c      do ir=1,ngrid
c        write(14,'(1x,f7.2,4(1pg12.4))')
c     &  rstep*(ir-1),rpot(ir),rpot(ir+ngrid),
c     &  rbmass(ir),rbmass(ir+ngrid)
c      enddo
c
c     derivative tensor of 'rbmass' on 'drbmas'
c
      call grad(rbmass,drbmas,2)
c
      if(t4.eq.zero) return
c
c     terms depending on 'sjs'
c
      t4co   = -b4
      t4q    = -b4-b4p
      do ir=2,ngrid
        irc    = ngrid+ir
        ps     = psjs(ir)
        psc    = psjs(irc)
        rpot(ir)  = t4q*ps+t4co*psc + rpot(ir)
        rpot(irc) = t4q*psc+t4co*ps + rpot(irc)
      enddo
c
c     the spin-orbit potential (iftm=1 not yet implemented)
c
      call derv17(prho(1),work1(1),parity,ngrid,rstep)
      call derv17(prho(ngrid+1),work2(1),parity,ngrid,rstep)
      do 310 ilamp=1,2
      lamp   = lrpa-3+ilamp+ilamp
      if(lamp.ge.0) then
       if(lamp.lt.lrpa) then
        facder = sqrt(alrpa/mlrpa)
        facrad = (alrpa+one)*facder
       else
        facder = -sqrt((alrpa+one)/mlrpa)
        facrad = -alrpa*facder
       endif
       do 320 ir=2,ngrid
       irn    = ngrid+ir
cold       rlspr  = ( facder*work1(ir)+facrad*prho(ir)*rinvs(ir) )*(b4+b4p)
cold       rlsne  = ( facder*work2(ir)+facrad*prho(irn)*rinvs(ir) )*b4
cold       rpotls(ir,ilamp)  = rlspr+rlspr+rlsne
cold       rpotls(irn,ilamp) = rlspr+rlsne+rlsne
       rlspr  = ( facder*work1(ir)+facrad*prho(ir)*rinvs(ir) )
       rlsne  = ( facder*work2(ir)+facrad*prho(irn)*rinvs(ir) )
       rpotls(ir,ilamp)  = rlspr*(b4+b4p)+rlsne*b4
       rpotls(irn,ilamp) = rlspr*b4+rlsne*(b4+b4p)
       if(iftm.eq.1) then
        rlspr  = psj(ir,ilamp)
        rlsne  = psj(irn,ilamp)
        rlssum = t1p2x*(rlspr+rlsne)
        rpotls(ir,ilamp)  = t1m20*rlspr-rlssum + rpotls(ir,ilamp)
        rpotls(irn,ilamp) = t1m20*rlsne-rlssum + rpotls(irn,ilamp)
       endif
  320  continue
       parls  = (mod(lrpa,2)*2-1)
       call extrap(rpotls(1,ilamp),parls)
       call extrap(rpotls(ngrid+1,ilamp),parls)
      else
       do ir=2,ngrid
         rpotls(ir,ilamp) = zero
         rpotls(ngrid+ir,ilamp) = zero
       enddo
      endif
  310 continue
c      iftest = 0
c      if(alrpa.eq.1 .and. iftest.eq.1) then
c        write(6,'(a,20(/1x,5g12.4))')
c     &       ' test prho:',
c     &       (prho(ir),ir=1,ngrid)
c        call derv17(rho(1),work1(1),one,ngrid,rstep)
c        write(6,'(a,10(/1x,5g12.4))') ' test d/dr rho:',
c     &   (-41.505*work1(ir),ir=1,ngrid)
c       endif
c
      return
      end
c-----rputil------------------------------------------------------
c
      subroutine rputil
      return
      end
c-----div   ------------------------------------------------------
c
      subroutine div(vecfld,divfld,maxfld)
c
      include 'params.inc'
c
c     computes the divergence of the vector-field 'vecfld' on
c     the scalar field 'divfld'.
c     all components are handled for 'maxfld'=2 whereas only
c     the proton part is handled for 'maxfld'=1.
c     the according collective angular momentum 'alrpa' is
c     transferred via 'linres.inc' .
c     uses 'hlinrs' as workspace.
c
      dimension vecfld(kg1*maxfld,2),divfld(kg1*maxfld)
!      dimension vecfld(kgr2,2),divfld(kgr2)
!      REAL(8), INTENT(IN OUT) :: vecfld(*,2),divfld(*)
c
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
c
c--------------------------------------------------------------------
c
      parity = 1-2*mod(lrpa,2)
c
      do ir=1,ngrid*maxfld
        divfld(ir) = zero
      enddo
c
      do 10 ilamp=1,2
c
c     the derivative potential, contracted over the 'ilamp'.
c
      if(ilamp.eq.1) then
       facder =  sqrt(alrpa/(alrpa+alrpa+one))
       facrad = -facder*(alrpa-one)
      else
       facder = -sqrt((alrpa+one)/(alrpa+alrpa+one))
       facrad = facder*(alrpa+one+one)
      endif
c     call derv15(vecfld(1,ilamp),hlinrs,-parity,ngrid,rstep)
      call derv17(vecfld(1,ilamp),hlinrs,-parity,ngrid,rstep)
      do ir=2,ngrid
        rinv   = rinvs(ir)
        divfld(ir) = ( facder*hlinrs(ir)
     &                + facrad*vecfld(ir,ilamp)*rinv )
     &              + divfld(ir)
      enddo
      if(maxfld.gt.1) then
c      call derv15(vecfld(ngrid+1,ilamp),hlinrs,-parity,ngrid,rstep)
       call derv17(vecfld(ngrid+1,ilamp),hlinrs,-parity,ngrid,rstep)
       do ir=2,ngrid
         rinv   = rinvs(ir)
         divfld(ir+ngrid) = ( facder*hlinrs(ir)
     &                 + facrad*vecfld(ngrid+ir,ilamp)*rinv )
     &               + divfld(ir+ngrid)
      enddo
      endif
c
   10 continue
c
      call extrap(divfld(1),parity)
      if(maxfld.gt.1) call extrap(divfld(ngrid+1),parity)
c
      return
      end
c-----rot   ------------------------------------------------------
c
      subroutine rot(vecfld,rotfld,maxfld)
c
      include 'params.inc'
c
c     computes the rotation of the vector-field 'vecfld' on
c     the scalar field 'rotfld'.
c     (this "scalar" field represents the one vector spherical
c      harmonics "y_{llm}").
c     all components are handled for 'maxfld'=2 whereas only
c     the proton part is handled for 'maxfld'=1.
c     the according collective angular momentum 'alrpa' is
c     transferred via 'linres.inc' .
c     uses 'hlinrs' as workspace.
c
      dimension vecfld(kgr2,2),rotfld(kgr2)
c
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
c
c--------------------------------------------------------------------
c
      parity = 1-2*mod(lrpa,2)
c
      do ir=1,ngrid*maxfld
        rotfld(ir) = zero
      enddo
!      WRITE(14,'(a,2(1pg13.5))') 'ROT:',
!     &  SUM(vecfld(1:ngrid,1)*rsquar(1:ngrid)),
!     &  SUM(vecfld(1:ngrid,2)*rsquar(1:ngrid))
c
      do 10 ilamp=1,2
      lamp   = lrpa-3+2*ilamp
      alamp  = lamp
      if(lamp.lt.0) go to 10
c
c     compute 'rotfld' from 'vecfld' as rotation
c
cold  if(ilamp.eq.1) then
cold   facder = sqrt(alrpa/mlrpa)
cold   facrad = -facder*(alrpa-one)
cold  else
cold   facder = sqrt(alamp/mlrpa)
cold   facrad = facder*(alamp+one)
cold  endif
      if(ilamp.eq.1) then
       facder = sqrt((alrpa+one)/mlrpa)
       facrad = -facder*(alrpa-one)
      else
       facder = sqrt(alrpa/mlrpa)
       facrad = facder*(alrpa+one+one)
      endif
c
      do 20 iso=1,maxfld
      if(iso.eq.1) then
       iradd  = 0
      else
       iradd  = ngrid
      endif
c
      call derv17(vecfld(iradd+1,ilamp),hlinrs,-parity,ngrid,rstep)
!      call derv17(vecfld(iradd+1,ilamp),hlinrs,parity,ngrid,rstep)
      do ir=2,ngrid
        irn    = iradd+ir
        rotfld(irn) = facder*hlinrs(ir)
     &             +facrad*vecfld(irn,ilamp)*rinvs(ir) + rotfld(irn)
      enddo
!      WRITE(14,'(a,i3,5(1pg13.5))') 'ROT2:',ilamp,
!     &   SUM(rotfld(1:ngrid)*rsquar(1:ngrid)),
!     &   SUM(hlinrs(1:ngrid)*rsquar(1:ngrid))
c
   20 continue
   10 continue
c
      call extrap(rotfld(1),parity)
      if(maxfld.gt.1) call extrap(rotfld(ngrid+1),parity)
c
      return
      end
c-----grad  ------------------------------------------------------
c
      subroutine grad(sclfld,gradfl,maxfld)
c
      include 'params.inc'
c
c     computes the gradient of the scalar field 'sclfld' on
c     the vector field 'gradfl'.
c     (this "scalar" field represents the one vector spherical
c      harmonics "y_{llm}").
c     all components are handled for 'maxfld'=2 whereas only
c     the proton part is handled for 'maxfld'=1.
c     the according collective angular momentum 'alrpa' is
c     transferred via 'linres.inc' .
c     uses 'hlinrs' as workspace.
c
      dimension gradfl(kgr2,2),sclfld(kg1*maxfld)
c
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
      dimension workgrad(kg1)
c
c--------------------------------------------------------------------
c
      gradfl(1:kg1,1:2)=0D0
!      WRITE(14,'(a,2(1pg13.5))') 'GRAD:',
!     &  SUM(gradfl(1:ngrid,2)*rsquar(1:ngrid))
      workgrad=0D0
      parity = 1-2*mod(lrpa,2)
c
c     derivative tensor of 'sclfld'
c
      do 110 iso=1,maxfld
      isoadd = (iso-1)*ngrid
!     call derv15(sclfld(isoadd+1),gradfl(isoadd+1,2),
      call derv15(sclfld(isoadd+1),workgrad,
     &            parity,ngrid,rstep)
      do 110 ilamp=1,2
      lamp   = lrpa-3+ilamp+ilamp
      if(lamp.ge.0) then
       if(lamp.lt.lrpa) then
        facder = sqrt(alrpa/mlrpa)
        facrad = (alrpa+one)*facder
       else
        facder = -sqrt((alrpa+one)/mlrpa)
        facrad = -alrpa*facder
       endif
!      gradfl(isoadd+1,ilamp) = (facder+facrad)*gradfl(isoadd+1,2)
       gradfl(isoadd+1,ilamp) = (facder+facrad)*workgrad(isoadd+1)
       do ir=2,ngrid
!         gradfl(isoadd+ir,ilamp) = facder*gradfl(isoadd+ir,2)
         gradfl(isoadd+ir,ilamp) = facder*workgrad(isoadd+ir)
     &                     +facrad*sclfld(isoadd+ir)*rinvs(ir)
       enddo
       parmas = (mod(lrpa,2)*2-1)
       call extrap(gradfl(isoadd+1,ilamp),parmas)
      else
       do ir=2,ngrid
         gradfl(isoadd+ir,ilamp) = zero
       enddo
      endif
  110 continue
c
!      WRITE(14,'(a,2(1pg13.5))') 'GRAD:',
!     &  SUM(gradfl(1:ngrid,1)*rsquar(1:ngrid)),
!     &  SUM(gradfl(1:ngrid,2)*rsquar(1:ngrid))
      return
      end
c-----rotsig------------------------------------------------------
c
      subroutine rotsig(sigcur,rsigm,ilamp)
c
      include 'params.inc'
c
c     computes the rotation of the spin-vector field 'sigcur'
c     for component 'ilamp'. result on rsigm.
c     the spin-vector field has one component of the type "y_{llm}"
c     the result is a current-vector with components
c     "y_{ll+1m}" for 'ilamp'=2 and "y_{ll-1m}" for 'ilamp'=1.
c
c     the rotation is build only for one isospin at a time.
c
c     the routine uses workspace 'work1' from common /linres/.
c
      dimension sigcur(kg1),rsigm(kg1)
c
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
c
c--------------------------------------------------------------------
c
      parcur = (1-mod(lrpa,2)*2)
c
      lamp   = lrpa-3+2*ilamp
      alamp  = lamp
      if(lamp.ge.0) then
       if(ilamp.eq.1) then
        facder = sqrt((alrpa+one)/mlrpa)
        facrad = facder*(alrpa+one)
       else
        facder = sqrt(alrpa/mlrpa)
        facrad = -facder*alrpa
       endif
       call derv17(sigcur,work1(1),parcur,ngrid,rstep)
       do ir=2,ngrid
         rsigm(ir) = facder*work1(ir)+facrad*sigcur(ir)*rinvs(ir)
       enddo
       call extrap(rsigm(1),-parcur)
      endif
c
      return
      end
c-----hlspsi------------------------------------------------------
c
      subroutine hlspsi(wfin,hwfout,jsa,iso)
c
      include 'params.inc'
c
c     action of the l*s part of the mean-field hamiltonian
c
c     where:
c       wfin    = input wave function field ,
c       hwfout  = output wave function field ,
c       jsa     = combined angular momentum of wfin and hwfout,
c       iso     = isospin of the wave functions,
c       subfac  = factor for subtraction of wfin.
c
c     the fields for the mean-field hamiltonian are transferred
c     via common /wavfun/.
c
      dimension wfin(kg1),hwfout(kg1)
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
c
c     the quantum numbers of the wave function
c      la      = orbital angular momentum of state  wfin>
c      parity  = parity for the radial wave function (note that it
c                contains an extra factor 'r' such that the 1s1/2
c                state has parity=-1 etc)
c      isoadd  = stack-pointer for potential fields pot, potls
c                and isoadd
c      xl0     = spin-robit splitting for state  wfin>
c
      la     = jsa/2
      xlp    = (la*la+la)
      parity = (mod(la,2)*2-1)
      isoadd = (iso-1)*ngrid
      xl0    = wls(jsa)
c
      hwfout(1) = zero
      do ir=2,ngrid
        irs    = ir+isoadd
        potact = xl0*potls(irs)
        hwfout(ir) = potact*wfin(ir)
      enddo
c
      return
      end
c-----coulml------------------------------------------------------------
c
      subroutine coulml(crho,coulpt,langmo)
c
      include 'params.inc'
c
c     computes coulomb potential 'coulpt'
c     for given charge density 'crho' of multipolarity 'langmo'.
c     the routine is build analogous to the laplace-inversion 'invert'
c     and specialised to the case of the coulomb field.
c
      dimension crho(kg1),coulpt(kg1)
c
      include 'potinv.inc'
      include 'skypar.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'constn.inc'
c
       data fsx,ezw/.83333333,.83333333e-1/
      data vpeq3/6.0317561/
cmc   'vpeq3' is 'e**2*4*pi/3'. value changed to use e**2=2
ccluster       data  vpeq3/8.3775804/
c
c-----------------------------------------------------------------------
c
c      square-root of inverse eff. mass on 'sqb',
c      second derivative of 'sqb' on 'd2b'.
c      finally assembling  the effective potential on 'poteff'.
c
       xlp    = langmo+langmo*langmo
       do i=2,ngrid
          poteff(i) = xlp/rsquar(i)
       enddo
c
c      transform to effective r.h.s., store result on r.
c
       qlmom  = zero
       vor    = 3.0*vpeq3*abs(acoul)
       do i=1,ngrid
         ract   = (i-1)*rstep
         coulpt(i)   = crho(i)*ract*vor
         qlmom  = crho(i)*ract**(2.0+langmo) + qlmom
       enddo
       qlmom  = qlmom*rstep*fourpi/(2.0*langmo+one)
c
c     prepare parameters for kinetic energy and loop-ends.
c
       offdig = 1./(rstep*rstep)
       offdg2 = offdig+offdig
       nmaxm  = ngrid-1
c
c      first pass: eliminating the lower off-diagonal element.
c
       poteff(1) = zero
       binv(1) = zero
       y      = one
       if(langmo.eq.1) then
        ainv(1) = zero
       else
        ainv(1) = one
       endif
       x      = zero
       rim1   = coulpt(1)
c
       do i=2,nmaxm
c                                         forward elim. for r.h.s.
         x      = ezw*(coulpt(i+1)+rim1)+fsx*coulpt(i)- x*(binv(i-1)/y)
         rim1   = coulpt(i)
         coulpt(i)   = x
c                                         forward elim. for matrix
         ui     = poteff(i)
         z      = -offdig+ezw*ui
         binv(i) = z
         y      = offdg2+fsx*ui-z*(binv(i-1)/y)
         ainv(i) = y
       enddo
c                                         last elim. step explicitely
c                                         because coulpt(i+1) is not availabl
       coulpt(ngrid) = ezw*rim1+fsx*coulpt(ngrid)-x*(binv(nmaxm)/y)
       ui     = poteff(ngrid)
       z      = -offdig+ezw*ui
       binv(ngrid) = z
       y      = offdg2+fsx*ui-z*(binv(nmaxm)/y)
       ainv(ngrid) = y
c
c      in the third pass resolve the remaining upper triangular equation
c      backward.
c
       x      = 1.44000*qlmom/(rstep*ngrid-rstep)**langmo
       coulpt(ngrid) = x*rinvs(ngrid)
       do i=nmaxm,2,(-1)
         x      = (coulpt(i)-x*binv(i+1))/ainv(i)
         coulpt(i)   = x*rinvs(i)
       enddo
       parity = 1-2*mod(langmo,2)
       call extrap(coulpt,parity)
c
       return
       end
ctestc-----laplact------------------------------------------------------------------
ctestc
ctest      subroutine laplact(fin,fout,d2work,ll,subfac)
ctestc
ctest      include 'params.inc'
ctestc
ctestc     action of the Laplacian
ctestc
ctestc     where:
ctestc       fin    = input field ,
ctestc       fout   = output field ,
ctestc       d2work = workspace for second derivative
ctestc       ll     = angular momentum  of fin and fout,
ctestc       subfac = constant subtracted from Laplacian operator
ctestc
ctestc     grid properties are transferred via common /grid/
ctestc
ctest      dimension fin(kg1),fout(kg1),d2work(kg1)
ctestc
ctest      include 'grid.inc'
ctest      include 'skypar.inc'
ctest      include 'constn.inc'
ctestc
ctestc------------------------------------------------------------------------------
ctestc
ctest      xlp    = (ll*ll+ll)
ctest      parity = (1-mod(ll,2)*2)
ctestc
ctestc     augment by r-factor
ctestc
ctest      do ir=1,ngrid
ctest        fout(ir)  = (ir-1)*rstep*fin(ir)
ctest      enddo
ctestc
ctestc     the derivatives of the field
ctestc
ctest      call derv25(fout,d2work,-parity,ngrid,rstep)
ctestc
ctestc     the action of the Laplacian operator
ctestc
ctest      do ir=2,ngrid
ctest        fout(ir) = ((xlp*zentrf(ir)-subfac)*fout(ir) - d2work(ir))
ctest     &            *rinvs(ir)
ctest      enddo
ctest      call extrap(fout,parity)
ctestc
ctest      return
ctest      end

c-----projec------------------------------------------------------
c
      subroutine projec(awf,js,nstack,nmax)
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
       do ir=2,ngrid
         sum   = rfact(ir)*awf(ir)*wf(nadd+ir) + sum
       enddo
       sum   = rstep*sum
c
c      project 'wf'
c
       do ir=2,ngrid
         awf(ir) = awf(ir) - sum*wf(nadd+ir)*rinvs(ir)
       enddo
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
c-----overlp------------------------------------------------------
c
      function overlp(awf,bwf,rdummy,idummy)
c
      include 'params.inc'
c
c     the overlap of 'awf' with 'bwf' where both wave functions
c     are from the linear response class which does not contain
c     the volume factor 'r' but the angular factor 'sqfrpi'.
c
      dimension awf(kg1),bwf(kg1)
c
      include 'grid.inc'
c
      nsw    = mod(ngrid,4)
      if(nsw.eq.3) then
       sum    = rsquar(1)*awf(1)*bwf(1)
     &         + rsquar(2)*awf(2)*bwf(2)
     &         + rsquar(3)*awf(3)*bwf(3)
      elseif(nsw.eq.2) then
       sum    = rsquar(1)*awf(1)*bwf(1)
     &         + rsquar(2)*awf(2)*bwf(2)
      elseif(nsw.eq.1) then
       sum    = rsquar(1)*awf(1)*bwf(1)
      else
       sum    = 0D0
      endif
      nsw    = 1+nsw
      do ir=nsw,ngrid,4
        sum    = rsquar(ir)*awf(ir)*bwf(ir)
     &          + rsquar(ir+1)*awf(ir+1)*bwf(ir+1)
     &          + rsquar(ir+2)*awf(ir+2)*bwf(ir+2)
     &          + rsquar(ir+3)*awf(ir+3)*bwf(ir+3)
     &          + sum
      enddo
      overlp = rstep*sum
c
      return
      end
c-----multr ------------------------------------------------------
c
      subroutine multr(fldin,fldout,rdummy,idummy)
c
      include 'params.inc'
c
      dimension fldin(kg1),fldout(kg1)
c
      include 'grid.inc'
c
      do ir=1,ngrid
        fldout(ir) = rfact(ir)*fldin(ir)
      enddo
c
      return
      end
c-----qnorm ------------------------------------------------------
c
      function qnorm(qwf)
c
      include 'params.inc'
c
c     the expectation value of the q-operator
c
c          <  qwf*qwf  >
c
c     where qwf is a hermitian operator handled by wave
c     functions.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      dimension qwf(kx11,klrpa)
c
c--------------------------------------------------------------------
c
      sum    = zero
c
      include 'initdo.inc'
c
c     now all quantum numbers are provided. do the space integral.
c
      do ir=2,ngrid
        factor = qwf(nadd+ir,jdif)
        sum    = rsquar(ir)*factor*factor + sum
      enddo
c
      include 'enddo.inc'
c
      qnorm  = rstep*sum
c
      return
      end
c-----qpover-----------------------------------------------------------------
c
      function qpover(qwf,pwf)
c
      include 'params.inc'
c
c     the overlap between a q- and a p-operator
c
c          <  qwf*pwf  >
c
c     where qwf and pwf are hermitian operators handled by wave
c     functions.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      dimension qwf(kx11,klrpa),pwf(kx11,klrpa)
c
c--------------------------------------------------------------------
c
      sum    = zero
c
      include 'initdo.inc'
c
c     now all quantum numbers are provided. do the space integral.
c
      sum=overlp(qwf(nadd+1,jdif),pwf(nadd+1,jdif),rdummy,idummy)+sum
c
      include 'enddo.inc'
c
      qpover  = sum
cold      qpover  = rstep*sum
c
      return
      end
c-----multwf------------------------------------------------------
c
      subroutine multwf(qwf,factor)
c
      include 'params.inc'
c
c     multiplies response wavefunction on 'qwf' with constant
c     factor 'factor'.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      dimension qwf(kx11,klrpa)
c
c--------------------------------------------------------------------
c
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     now all quantum numbers are provided. do the copying
c
      do ir=1,ngrid
        qwf(nadd+ir,jdif) = factor*qwf(nadd+ir,jdif)
      enddo
c
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
      return
      end
c-----admlwf------------------------------------------------------
c
      subroutine admlwf(awf,bwf,afac,bfac)
c
      include 'params.inc'
c
c     multiplies 'awf' with 'afac' and 'bwf' with 'bfac'. adds
c     both terms and returns result on 'awf'.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      dimension awf(kx11,klrpa),bwf(kx11,klrpa)
c
c--------------------------------------------------------------------
c
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     now all quantum numbers are provided. do the copying
c
      do ir=1,ngrid
        awf(nadd+ir,jdif) = afac*awf(nadd+ir,jdif)
     &                     +bfac*bwf(nadd+ir,jdif)
      enddo
c
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
      return
      end
c-----copywf------------------------------------------------------
c
      subroutine copywf(qwf,pwf)
c
      include 'params.inc'
c
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      dimension pwf(kx11,klrpa),qwf(kx11,klrpa)
c
c--------------------------------------------------------------------
c
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     now all quantum numbers are provided. do the copying
c
      do ir=1,ngrid
        pwf(nadd+ir,jdif) = qwf(nadd+ir,jdif)
      enddo
c
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
      return
      end
c-----addwfs------------------------------------------------------
c
      subroutine addwfs(vec,wfnam,ichr,maxmod)
c
      include 'params.inc'
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      dimension vec(krpmax)
c
      character*(*) wfnam
      character*(7) qop
      character*2 ichr(krpmax)
c
c------------------------------------------------------------------------------
c
      iposnm = len(wfnam)+1
      qop    = '      '
      qop(1:iposnm-1) = wfnam
      do 510 imod2=1,maxmod
      qop(iposnm:iposnm+1) = ichr(imod2)
      call buffql(savlr,'i',qop)
c
      include 'initdo.inc'
c
      fac    = vec(imod2)
      if(imod2.eq.1) then
       do ir=1,ngrid
         qlrout(nadd+ir,jdif) = fac*savlr(nadd+ir,jdif)
       enddo
      else
       do ir=1,ngrid
         qlrout(nadd+ir,jdif) = fac*savlr(nadd+ir,jdif)
     &                       +qlrout(nadd+ir,jdif)
       enddo
      endif
c
      include 'enddo.inc'
c
  510 continue
c
      return
      end
c-----suppos------------------------------------------------------
c
      subroutine suppos(fac1,fac2,name1,name2,nameot)
c
      include 'params.inc'
c
c     reads fields with key_name 'name1' on 'qlrout' and
c     with key_name 'name2' on 'plrout'.
c     superposes
c       'savlr' = 'fac1'*'qlrout' + 'fac2'*plrout' ,
c     and writes result to key_name 'nameot' in file 18.
c
      character*(7) name1,name2,nameot
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
      call buffql(qlrout,'i',name1)
      call buffql(plrout,'i',name2)
c
      include 'initdo.inc'
c
      do ir=1,ngrid
        savlr(nadd+ir,jdif) = fac1*qlrout(nadd+ir,jdif)
     &                     + fac2*plrout(nadd+ir,jdif)
      enddo
c
      include 'enddo.inc'
c
      call buffql(savlr,'o',nameot)
c
      return
      end
c-----filter-------------------------------------------------------------------
c
      subroutine filter(fldin,fldout,lfilt)
c
      include 'params.inc'
c
c     cutoff in momentum space at 'resdmp' (from /linres/).
c     input on 'fldin'. output on 'fldout'. ('fldin' and 'fldout'
c     can be the same fields).
c     the workspaces 'aux1' and 'aux2' are needed.
c
      dimension fldin(kg1),fldout(kg1)
      dimension auxfou(kgr2)
c
      include 'wavfun.inc'
      include 'linres.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
      equivalence (aux1(1),auxfou(1))
c
      data ressav/1.0D10/
c
c-------------------------------------------------------------------
c
      if(resdmp.gt.ressav) then
       write(6,'(a,f8.4)') ' cutoff in momentumspace=',resdmp
       write(14,'(a,f8.4)') ' cutoff in momentumspace=',resdmp
       write(6,'(a,30(/1x,5g12.4))')
     &  ' fldin:',(fldin(ir),ir=1,ngrid)
      endif
c
c     prepare constants
c
      delql  = 0.2*delq
c      ressav = resdmp
      iqfilt = resdmp/delql+1
      iqfilt = min(iqfilt,kgr2)
      delql  = resdmp/(iqfilt-1)
c      resdmp = (iqfilt-1)*delql*1.000001
c      if(resdmp.ne.ressav) then
c       write(6,'(a,f8.4)') ' cutoff in momentumspace=',resdmp
c       write(14,'(a,f8.4)') ' cutoff in momentumspace=',resdmp
c      endif
c
c     formfactor of 'fldin' stored on workspace 'aux1'
c
      do iq=1,iqfilt
        q1     = iq*delql-delql
        sum    = zero
        do ir=2,ngrid
          r      = rfact(ir)
          sum    = r*r*bessel(r*q1,lfilt)*fldin(ir) + sum
        enddo
        auxfou(iq) = rstep*sum *(one-((iq-one)/iqfilt)**2)
c      auxfou(iq) = rstep*sum
      enddo
c
c     fourier-bessel backtransformation
c
      do ir=1,ngrid
        r      = rfact(ir)
        sum    = zero
        do iq=2,iqfilt
          q1     = iq*delql-delql
          sum    = q1*q1*bessel(r*q1,lfilt)*auxfou(iq) + sum
        enddo
        fldout(ir) = 2.0*delql*sum/pi
      enddo
c
      if(resdmp.gt.ressav) then
       write(6,'(a,30(/1x,5g12.4))')
     &  ' fldout:',(fldout(ir),ir=1,ngrid)
       ressav = resdmp
      endif
c
      return
      end
c-----buffql-------------------------------------------------------------------
c
      subroutine buffql(awf,direc,name)
c
      include 'params.inc'
      parameter (knames=5*krpmax+1 , kspsav=512,
     &           kblock=kspsav/kvarln,
     &           kspfil=kvarln*kblock)
c
c     buffering response wave function 'awf' on file 18
c
c     'direc' defines the direction of the action:
c             'c'  = clear pointer fields (to initialise or reset)
c             'i'  = input, i.e. reading 'awf' from pointer 'name'
c             'o'  = output, i.e. writing reading 'awf' to pointer 'name'
c
c     'name'  is a key name for the wave function on 'awf':
c             the routine looks up the table of key names,
c             if it is already defined then the according pointer
c             is accessed,
c             if it is new then 'name' is appended to the tableau
c             and an according pointer is assigned (in the order
c             of definitions)
c
      character*1 direc
      character*(7) names(knames)
      character*(*) name
      dimension awf(kx11,klrpa)
c
      include 'wavfun.inc'
      include 'linres.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'constn.inc'
c
      data names/knames*'xxxxxxx'/
      data lastnm,nblock,nblwf,nwflen/0,0,0,0/
      data ifop18/0/
c
c-----------------------------------------------------------------------------
c
      if(direc.eq.'c') then
c
c     open file 18 as scratch file for 'qlrout'
c
      if(ifop18.ne.0) then
       close(unit=18)
      endif
      open(unit=18,access='direct',form='unformatted',recl=kspfil,
     &      status='scratch')
      ifop18 = 1
c
c     clear pointer fields to initialise new set of wave functions
c
c                                    compute blocks of the wave functions
       nwflen = nmaxto*ngrid
       nblwf  = (nwflen-1)/kblock+1
       nblock = mlrpa*nblwf
c      write(6,'(3(a,i4))')
c    &  ' set in buffql: nwflen=',nwflen,'  nblwf=',nblwf,
c    &  '  nblock=',nblock
c                                    reset pointer fields
       do ipoint=1,knames
         names(ipoint) = 'xxxxxxx'
       enddo
       lastnm = 0
c
      elseif(direc.eq.'o') then
c
c     output to file 18
c
c                                    check old pointer, compute new pointer
       do ipoint=1,lastnm
         ipact  = ipoint
         if(name.eq.names(ipoint)) goto 29
       enddo
       ipact  = lastnm+1
       lastnm = ipact
       if(lastnm.gt.knames) stop ' too much pointers on file 18'
       names(lastnm) = name
       write(6,'(3a,i3,a)')
     &  ' name ',name,' has been stored. pointer nr.',lastnm,
     &  ' is assigned to it.'
   29  continue
c                                    write output
       nrecor = (ipact-1)*nblock
       nrec1  = nrecor
       do jdif=1,mlrpa
       do ip=1,nblwf
         nrecor = 1+nrecor
         iend   = ip*kblock
         istart = iend-kblock+1
         iend   = min(nwflen,iend)
         write(18,rec=nrecor) (awf(i,jdif),i=istart,iend)
       enddo
       enddo
c      write(6,'(3a,2i5)')
c    &  ' file with name=',name,' written from/to record=',
c    &  nrec1,nrecor
c
      elseif(direc.eq.'i') then
c
c     input from file 18
c
c                                    check pointer
       do ipoint=1,lastnm
         ipact  = ipoint
         if(name.eq.names(ipoint)) goto 59
       enddo
       write(6,'(a,a)')
     &   ' trying to read from undefined pointer on file 18. name=',
     &   name
       stop
   59  continue
c
       nrecor = (ipact-1)*nblock
       nrec1  = nrecor
       do jdif=1,mlrpa
       do ip=1,nblwf
         nrecor = 1+nrecor
         iend   = ip*kblock
         istart = iend-kblock+1
         iend   = min(nwflen,iend)
         read(18,rec=nrecor) (awf(i,jdif),i=istart,iend)
       enddo
       enddo
c      write(6,'(3a,2i5)')
c    &  ' file with name=',name,' read from/to record=',
c    &  nrec1,nrecor
c
      else
c      dummy = sqrt(zero)
c      write(6,'(a,2g12.4,a,g12.4)')
c    &    ' zero,one=',zero,one,' dummy=',dummy
       stop ' wrong argument "direc" in buffql'
      endif
c
      return
      end
c-----resmte-----------------------------------------------------
c
      subroutine resmte
c
      include 'params.inc'
      parameter(krpmat=krpmax*(krpmax+1)/2)
      parameter(kreczw=kgr2*kvarln,
     &          krec17=5*kreczw)
c
c     computes matrix elements of a particle-hole operator and
c     the optimal coupling strength if this operator is used
c     for a separable force.
c
c     is called from 'hafosu' if 'ifgsc'=4.
c
c     the options are read-in at the beginning of this routine.
c     they are:
c      lrpa    = angular momentum of the mode
c      cmfac   = strength of r_cm*r_cm residual interaction which
c                serves allow inversion of the hamiltonian in case
c                of the isoscalar dipole mode (cm motion).
c      ifirpa  = switch for inverse rpa-moments
c      itrpa   = number of iterations in inverse h-action
c      epsrpa  = stepsize in inverse h-action
c      dabrpa  = required precision in inverse h-action
c      ifgrho  = switch for computation of ground-state correlarions
c      ifanlr  = level of analysis; lowest is -1, final modes are
c                composed explicitely and widths are computes for 0,
c                even higher values give more detailed printout.
c
c     the selection of the operator is read in this routine
c     and evaluated in 'qinit'. the mode parameters are
c      typmod  = local form
c                 'mom=m' --> initialize r**(lrpa+m)
c                 'bes=q' --> initialize j_lrpa(q*r/r)
c                 'ph=(p,h)' --> explicit particle-hole state (p,h)
c                 'end'   --> terminates input of modes
c                the number of modes is determined dynamically by reading
c                until the 'end' option
c      isorpa  = isospin of the mode (1,2)
c      jselec  = option for picking selective (j,ja,iso)-transitions from
c                a local operator
c                '(js,jsa,iso)' --> only transition jsa -> js is taken
c                'all+n'    --> all transitions separately up to js=jsa+n
c                'nojj'     --> no selective picking
c      isplit  = splitting of upper and lower momentum component:
c                 0 --> no splitting = strictly local operator
c                 1 --> select lower component of corresponding momentum
c                 2 --> select upper component of corresponding momentum
c      ihwgt   = hamiltonian weight h**ihwgt on q
c      ih0wgt  = mean-field hamiltonian weight h0**ihwgt on q
c      ipwgt   = hamiltonian weights on p, p = h**ipwgt q (use 0,1,-1);
c                standard is ipwgt=1, for 'ph' use ipwgt=0, low modes
c                may prefer ipwgt=-1
c      iprshl  = switches projection on-off valence shell
c                0 --> no projection
c                1 --> projects q on valence shell
c                2 --> projects q off valence shell
c                the valence shell is defined by the unoccupied levels
c                in the selection provided by nmaxp and nmaxn.
c
c
      include 'skypar.inc'
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'numpar.inc'
      include 'constn.inc'
      include 'linres.inc'
c
      dimension nparsv(krpmax),nholsv(krpmax),isophs(krpmax)
c
      character*2 ichr(1)
      character*16 typmod
      character*12 jselec
      character*50 toper
      character*(7) qop,pop,qbop,qrop,prop
      data ichr/'1'/
      data     qop,      pop,    qbop,    qrop,    prop
     &    /'qin   ','p     ','[h,p ]','qres  ','pres  '/
      data resfac,cmfac
     &    /   1.0D0,  0D0/
c
      data ifop14,ifop17/0,0/
c
c  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
c
c      check bounds
c
       if(krpmax.gt.62) stop ' krpmax too large'
c
c      read parameters for 'sumrul'
c
       read(5,*)     lrpa,resfac,cmfac,ifsepi,resdmp
       read(5,*)     itrpa,epsrpa,dabrpa,ifanlr,ifirpa,ifgrho
c
c     open file 14 for concentrated output of matrix elements
c
      if(ifgrho.ne.0)
     &  open(unit=17,status='scratch',form='unformatted',
     &      access='direct',recl=5*kvarln*ngrid2)
      if(ifop14.eq.0) then
       open(unit=14,status='unknown',form='formatted',file='for014')
       ifop14 = 1
      endif
c
c     write a header
c
      write(14,'(a/1x,2a/1x,2(a,i4),5x,2(a,i4))')
     &  '1',
     &  '  title:',title,
     &  'nucleus:  nprot=',nprot,'  nneut=',nneut,
     &  ' shells:  nmaxp=',nmaxp,'  namxn=',nmaxn
      write(14,'(1x,a,5g13.6/7x,a,5g13.6/7x,a,g13.6)')
     &  'force: t0-t4=',t0,t1,t2,t3,t4,
     &        ' x0-x3=',x0,x1,x2,x3,b4p,
     &        ' power=',power
      write(14,'(7x,3(a,g13.6)/)')
     &        ' tyuk=',tyuk,' xyuk=',xyuk,' ayuk=',ayuk
      write(14,'(1x,2(a,f10.3)/18x,2(a,f10.3))')
     &  'some properties: rms-radius=',rms,'    diffr.radius=',rdms,
     &  '   surface=',surf,'          energy=',etotal
c
c      sumrul-part
c
c      first supply derivative of 'bmass' on 'd1b' for use in 'hpsi'
c
      call derv17(bmass(1),d1b(1),one,ngrid,rstep)
      call derv17(bmass(ngrid+1),d1b(ngrid+1),one,ngrid,rstep)
c
c     provide derivative of "rho"
c
      call derv17(rho(1),d1r(1), one,ngrid,rstep)
      call derv17(rho(ngrid+1),d1r(ngrid+1), one,ngrid,rstep)
c
c     write options and parameters for rpa part
c
      write(6,'(/1x,a,i3,2(a,g10.3)/10x,a,i5,a,g10.3)')
     &   ' matrix elements for: lrpa=',lrpa,'  res.int.fac=',resfac,
     &   '  cminv.fac=',cmfac,
     &   ' ifsepi=',ifsepi,'   resdmp=',resdmp
      write(6,'(21x,a,i5,a,f7.3,a,g12.4)')
     &    '  itrpa=',itrpa,' epsrpa=',epsrpa,' daprpa=',dabrpa
      write(14,'(/1x,a,i3,2(a,g10.3)/10x,a,i5,a,g10.3)')
     &   ' matrix elements for: lrpa=',lrpa,'  res.int.fac=',resfac,
     &   '  cminv.fac=',cmfac,
     &   ' ifsepi=',ifsepi,'   resdmp=',resdmp
c
c      set parameters inherited from 'extrsp'
c
       ifcmpl = 0
       enext  = zero
       enwid  = zero
c
c     set the other angular momentum parameters
c
      alrpa=lrpa
      mlrpa=2*lrpa+1
      call setred
c
c     initialise buffering for response wave functions
c
      call buffql(qlrout,'c','      ')
c
c     compute only the matrix elements of 'qrad'
c
       if((weight(nmaxp).gt.xsmall .and. weight(nmaxto).gt.xsmall))
     &  stop
     &   ' not enough unoccupied states to compute matrix elements'
        read(5,*) typmod,isorpa,jselec,isplit,ihwgt,ih0wgt,ipwgt,iprshl
c                                                          prepare Q-operators
        call qinit(typmod,jselec,isplit,ihwgt,ih0wgt,ipwgt,
     &             iprshl,imode,ifprim,iflagr,xkqrd,one,cmfac,
     &             nparsv,nholsv,isophs,toper)
        qop(4:5 ) = ichr(1)
c                                                          save in buffer
        call buffql(qlrout,'o',qop)
        pop(2:3 ) = ichr(1)
        qbop(5:6) = ichr(1)
        ifanly = 1+ifanlr
c                                                          prepare P-operators
        call moment(ovqq,ovqhq,ovqp,ovpp,ovpqb,ovqh0q,qqinv,widrpa,
     &              qrlpr,qrlne,resfac,cmfac,
     &              eshifq,eshift,eshifi,
     &              qop,pop,qbop,ipwgt,ifanly)
c
        resinq = mlrpa*eshifq/ovqq
        resint = mlrpa*eshift/ovqq
        resini = mlrpa*eshifi/ovqq
        write(14,'(a,a/1x,a,4(1pe12.4))')
     &    '  typmod     <q**2>  rpa-kappa ',
     &    '   tda-kappa inv.-kappa',
     &    typmod,ovqq,resinq,resint,resini
c                                                          save in buffer
        call buffql(qlrout,'i',qop)
        call primte(typmod,xkqrd,toper,
     &              ovqq,ovqh0q,qqinv,resinq,resint,resini)
        return
c
      end
c-----primte------------------------------------------------------
c
      subroutine primte(typmod,xkqrd,toper,
     &                  qwidfl,ovqh0q,qqinv,resinq,resint,resini)
c
      include 'params.inc'
c     print the matrix elements of a given multipole-operator.
c     uses 'qrad' as workspace.
c
c     the fields for the mean-field hamiltonian are transferred
c     via common /wavfun/.
c
      character*16 typmod
      character*50 toper
c
      include 'wavfun.inc'
      include 'result.inc'
      include 'skypar.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      data ifopen,xminrd/0,1.0D-10/
c
c--------------------------------------------------------------------
c
c     open protocol file nr. 10 and write a header
c
      open(unit=10,status='unknown',file='for010',form='formatted')
      write(10,'(3f10.2,a50)') anumb,aprot,aneut,title
      parity = 1-mod(lrpa,2)*2
      write(10,'(i5,f10.1,i5)') lrpa,parity,1-isorpa
      if(isorpa.eq.1) then
       resscl = -resinq
       resvct = zero
      else
       resscl = zero
       resvct = -resinq
      endif
      write(10,'(2g13.5,2f4.1,3g13.5)')
     &  resscl,resvct,zero,one,-resint,-resini,qwidfl
      lev    = nmaxp+nmaxn
      write(10,'(i5)') lev
c
c     print the s.p. energies in the rpa input format
c
      if(iwr.ge.-1) then
       write(14,'(/a/a)') ' print s.p.energies in rpa input format:',
     &    ' n  l 2j itz         s.p.energy'
       do 710 iso=1,2
       if(iso.eq.1) then
        nstack = 0
        nmax   = nmaxp
        isoadd = 0
       else
        nstack = nmaxp
        nmax   = nmaxn
        isoadd = kjmx
       endif
       do 710 na=1,nmax
       jsa    = jsp(na)
       la=(jsa/2)
       xja=(jsa-1)/2+0.5D0
       ja2=xja+xja
       neff   = na+nstack
       itz    = 2-iso
       write(14,'(4i3,8x,f10.5)') neff,la,ja2,itz,spenrg(neff)
       write(10,'(4i3,8x,f10.5)') neff,la,ja2,itz,spenrg(neff)
  710  continue
      endif
c
c     set up the local operator 'q'
c
      call initqr(typmod(1:3),xkqrd,powqr)
      if(iwr.ge.0) then
       write(14,'(/1x,a/1x,a,5g12.4,10(/1x,5g12.4))')
     &  'testprint of the radial operator:',
     &  ' protons:',(qrad(ir),ir=1,ngrid)
       write(14,'(1x,a,5g12.4,10(/1x,5g12.4))')
     &  'neutrons:',(qrad(ir+ngrid),ir=1,ngrid)
      endif
c
c     some settings before the loops start
c
      if(iwr.ge.-1)
     & write(14,'(/a,i3,a,a/a/a)')
     &  ' radial transition elements for multipole=',lrpa,
     &  ' and operator=',toper,
     &  ' given in the form:',
     &  '  nh   np    <operator>    <r**l>    red.mat.el.'
c
      nme    = 0
      ovqq   = zero
      ovqp0  = zero
      ovqpm1 = zero
c
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     now all a.m. quantum numbers are provided, run the loop
c     over 'nb'
c
      if(wgt.le.half) goto 30
      redmat = redume(xj,alrpa,xja)
      if(abs(redmat).le.xminrd) goto 30
c
c     set isospin factor to remove these from the matrix elements
c
      if(isorpa.eq.1) then
       unfac  = one
       facrla = one
      else
       if(iso.eq.1) then
        facrla  = aneut/anumb
        unfac  = one
       else
        facrla  = aprot/anumb
        unfac  = -one
       endif
      endif
c
c     loop over the b-states
c
      do 40 nb=1,nmax
      jsb   = jsp(nb)
      if(jsb.ne.js) goto 40
      lb    = jsb/2
      xlb   = lb
      xjb   = (jsb-1)/2+0.5D0
      wgtb  = weight(nb+nstack/ngrid)
      if(wgtb.gt.half) goto 40
      naddb = nstack + (nb-1)*ngrid
      subfcb = spenrg(nb+nstack/ngrid)
c
c     radial matrix element
c
      sumrl  = zero
      overlp = zero
      r      = zero
      do ir=2,ngrid
        r      = rstep+r
        sumrl  = wf(naddb+ir)*r**powqr*wf(nadd+ir) + sumrl
        overlp = wf(naddb+ir)*qlrout(nadd+ir,jdif)*r + overlp
      enddo
      sumrl  = facrla*rstep*sumrl
      overlp = rstep*overlp*unfac/redmat
c
c     print the result
c
      isonp=(iso-1)*nmaxp
      npart=nb+isonp
      nhole=na+isonp
      if(npart.gt.0 .and. iwrit.ge.-1) then
       nme=1+nme
       addend = (overlp*redmat)**2*abs(wgt-wgtb)
       ovqq   = addend + ovqq
       ovqp0  = addend*abs(subfac-subfcb) + ovqp0
       ovqpm1 = addend/abs(subfac-subfcb) + ovqpm1
       write(14,'(2i4,4(1pe13.4))') nhole,npart,overlp,sumrl,redmat,
     &       ovtest
       write(10,'(2i4,3(1pe13.4))') nhole,npart,overlp,sumrl,redmat
      endif
c
   40 continue
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
c     print size of phase space
c
      write(10,'(2i4,3(1pe13.4))') -01,0,zero,zero,zero
      write(14,'(/a,i5,2(a,g11.4)/)')
     &    ' explicit sum: phase space nme=',
     &    nme,' q-width=',ovqq,' e(1)-average=',ovqp0/max(ovqq,xminrd),
     &    '               e(-1)-average=',ovqq/max(ovqpm1,xminrd),
     &    ' sqrt(e(-1)*e(1))=',sqrt(abs(ovqp0*ovqpm1))
c
c     print wavefunctions
c
      call priwf(14)
c
c     compute moments for the gsc on the rms-radius by the collective
c     model   ( according to f. barranco ).
c
      parity = 1-mod(lrpa,2)*2
      call derv15(qrad(1),work1(1),parity,ngrid,rstep)
      call derv15(qrad(ngrid+1),work2(1),parity,ngrid,rstep)
      powmul = alrpa-one
      sfscal = zero
      sfdiff = zero
      smscal = zero
      smdiff = zero
      do ir=2,ngrid
        r      = (ir-1)*rstep
        rpow   = r**powmul
        vol    = r*r*rho(ir)
        sfscal = vol*(3.0*qrad(ir)+r*work1(ir))      + sfscal
        smscal = vol*(3.0+alrpa)*r*rpow              + smscal
        sfdiff = vol*rdms*(2.0*qrad(ir)/r+work1(ir)) + sfdiff
        smdiff = vol*rdms*(2.0+alrpa)*rpow           + smdiff
      enddo
      fac    = frpirs/aprot
      write(14,'(/2a/a,4g12.4)')
     &  ' delta-sigma**2 factors:    f-scal      f-diff   ',
     &  ' mult.-scal  mult.-diff',
     &  ' charge:                ',
     &  fac*sfscal,fac*sfdiff,fac*smscal,fac*smdiff
      write(10,'(4g13.5)')
     &  fac*sfscal,fac*sfdiff,fac*smscal,fac*smdiff
      do 320 ir=2,ngrid
      r      = (ir-1)*rstep
      rpow   = r**powmul
      vol    = r*r*rho(ngrid+ir)
      sfscal = vol*(3.0*qrad(ngrid+ir)+r*work2(ir))      + sfscal
      smscal = vol*(3.0+alrpa)*r*rpow                    + smscal
      sfdiff = vol*rdms*(2.0*qrad(ngrid+ir)/r+work2(ir)) + sfdiff
      smdiff = vol*rdms*(2.0+alrpa)*rpow                 + smdiff
  320 continue
      fac    = frpirs/anumb
      write(14,'(a,4g12.4)')
     &  ' total:                 ',
     &  fac*sfscal,fac*sfdiff,fac*smscal,fac*smdiff
      write(10,'(4g13.5)')
     &  fac*sfscal,fac*sfdiff,fac*smscal,fac*smdiff
c
c     compute missing strength and reevaluate kappa accordingly
c
      fracqq = ovqq/qwidfl
      fracqe = ovqp0/ovqh0q
      fracqi = ovqq/max(xsmall,qqinv)
      write(14,'(3(a,g12.4))')
     &   ' exhaustion of widths: qq=',fracqq,'  e*qq=',fracqe,
     &   ' qq/e=',fracqi
      write(10,'(3(a,g12.4))')
     &   ' exhaustion of widths: qq=',fracqq,'  e*qq=',fracqe,
     &   ' qq/e=',fracqi
      write(10,'(3(a,g12.4))')
     &   '  renormalised kappas: qq=',-resinq/fracqq,
     &   '  e*qq=',-resinq/fracqe,
     &   ' qq/e=',-resinq/fracqi
c
      return
      end
c-----redume------------------------------------------------------
c
      function redume(bj,al,aj)
c
      include 'params.inc'
c
c     reduced m.e. (y(bj)11y(al)11y(aj) edmonds notation
c
c     =(-1)**(ja+0.5+al)*sqrt(((2*aj+1)*(2*al+1)*(2*bj+1))
c     /(4*pi)) * wigner(aj,al,bj,0.5,0.0,-0.5)
c
c     wigner(j1,j2,j,m1,m2,m) is a three-j symbol
c
      include 'redmts.inc'
      include 'linres.inc'
c
      data half,zero,pi4,xlarg/0.5D0,0D0,12.5663704,1.0D20/
c
c     check table
c
      jsa    = aj+aj
      jdif   = bj-aj+al+1
      ilamp  = al-alrpa
      ttable = jsa.le.kjrpmx .and. jdif.le.klrpa .and. abs(ilamp).le.1
      ilamp  = 2+ilamp
      if(ttable) then
       x     = redmt(jsa,jdif,ilamp)
       if(x.lt.xlarg) then
        redume = x
        return
       endif
      endif
c
c     compute from scratch
c
      iphase = abs(bj-aj+al)
      phase  = 1-2*mod(iphase,2)
      iph    = (aj+half+al)
      iph    = 1-2*mod(iph,2)
      wig    =  wigner(bj,al,aj,half,zero,-half)
      if(wig.ne.zero) then
       redume = iph*sqrt( 8.0*(aj+half)*(al+half)*(bj+half)/pi4 )
     &         * wig*phase
      else
       redume = zero
      endif
c
      if(ttable) then
       redmt(jsa,jdif,ilamp) = redume
      else
       write(6,'(a)') ' non-standard call to "redume":'
      endif
c
      return
      end
c-----redcur------------------------------------------------------
c
      subroutine redcur(xja,xla,al,alamp,xj,xl,facder,facrad)
c
      include 'params.inc'
c
c     reduced matrix element for the current part
c      ( ja,la | t(lrpa,lamp)*nabla | j,l )
c     in edmonds notation ; t(l,ld) is a vector harmonics
c
c     attention: the quantum numbers are not checked for consistency.
c                this needs to be done by the calling routine.
c
c
c                                          /a b c\
c     sixjot(a,b,c,d,e,f)  is six-j symbol \d e f/
c
c
c     the radial part 'facrad' is computed for the derivatives on 'wf'
c     which contains the radial factor 'r'.
c
      include 'redmts.inc'
      include 'linres.inc'
c
      data zero,half,one,xlarg/0D0,0.5D0,1.0D0,1.0D20/
c
c     check table
c
      jsa   = xja+xla+half
      jdif  = xj-xja+lrpa+1
      if(alamp.eq.al+one) then
       ilamp  = 2
       ttable = .true.
      elseif(alamp.eq.al-one) then
       ilamp  = 1
       ttable = .true.
      else
       ttable = .false.
      endif
      ttable = jsa.le.kjrpmx .and. jdif.le.klrpa .and. al.eq.alrpa
     &         .and. ttable
      if(ttable) then
       x     = dercr(jsa,jdif,ilamp)
       if(x.lt.xlarg) then
        facder = x
        facrad =  radcr(jsa,jdif,ilamp)
        return
       endif
      endif
c
c     no table entry found. compute from scratch
c
      sumrad = zero
      sumder = zero
      consix = sixjot(xla,xja,half,xj,xl,al)
      do 10 lpcase=1,2
      xlp    = xl-3+lpcase+lpcase
      if(xlp.ge.zero .and. abs(xlp-xla).le.alamp
     &   .and. alamp.le.xlp+xla) then
        add    =  wigner(xla,alamp,xlp,zero,zero,zero)
     &           *sixjot(alamp,one,al,xl,xla,xlp)*consix
     &           *sqrt((xl+(lpcase-1))*(xlp+xlp+one))
        sumder = sumder+(1-2*mod(lpcase,2))*add
        sumrad = sumrad-add*(xl+(lpcase-1))
      endif
   10 continue
c
      iphase = (xj+half+xla+xl)
      fac    = sqrt((xj+half)*(al+half)*(alamp+half)*(xja+half)
     &              *(xla+half)*8.0/3.141592654)
     &        * (1-2*mod(iphase,2))
c
c     add "experimental" phase
c
      facder = fac*sumder
      facrad = fac*sumrad
c
c     save on table
c
      if(ttable) then
       dercr(jsa,jdif,ilamp) = facder
       radcr(jsa,jdif,ilamp) = facrad
      else
       write(6,'(a)') ' non-standard call to "redcur"'
      endif
c
      return
      end
c-----redspo------------------------------------------------------
c
c     subroutine redspo(xja,xla,al,alamp,xj,xl,facder,facrad)
      subroutine redspo(xj,xl,al,alamp,xja,xla,facder,facrad)
c
      include 'params.inc'
c
c     reduced matrix element for the current part
c      ( ja,la | t(lrpa,lamp)*nabla | j,l )
c     in edmonds notation ; t(l,ld) is a vector harmonics
c
c     attention: the quantum numbers are not checked for consistency.
c                this needs to be done by the calling routine.
c
c
c                                          /a b c\
c     sixjot(a,b,c,d,e,f)  is six-j symbol \d e f/
c
c
c     the radial part 'facrad' is computed for the derivatives on 'wf'
c     which contains the radial factor 'r'.
c
      include 'redmts.inc'
      include 'linres.inc'
c
      data zero,half,one,xlarg/0D0,0.5D0,1.0D0,1.0D20/
c
c     check table
c
      jsa   = xja+xla+half
      jdif  = xj-xja+lrpa+1
      if(alamp.eq.al+one) then
       ilamp  = 2
       ttable = .true.
      elseif(alamp.eq.al-one) then
       ilamp  = 1
       ttable = .true.
      else
       ttable = .false.
      endif
      ttable = jsa.le.kjrpmx .and. jdif.le.klrpa .and. al.eq.alrpa
     &         .and. ttable
      if(ttable) then
       x     = dersp(jsa,jdif,ilamp)
       if(x.lt.xlarg) then
        facder = x
        facrad = radsp(jsa,jdif,ilamp)
        return
       endif
      endif
c
c     no table entry found. compute from scratch
c
      sumrad = zero
      sumder = zero
      do 10 lpcase=1,2
      xlp    = xla-3+lpcase+lpcase
      if(xlp.ge.zero .and. abs(xlp-xl).le.alamp
     &   .and. alamp.le.xlp+xl) then
       do 20 jpcase=1,2
       xjp    = xlp+(jpcase-1.5)
       if(xjp.gt.zero .and. abs(xjp-xj).le.alamp
     &   .and. alamp.le.xjp+xj .and. abs(xjp-xja).le.one) then
        iphase = (xj+xja+al)
        add    = (1-2*mod(iphase,2))
     &          * sqrt((xjp+half)*(xla+(lpcase-1)))
     &          * sixjot(alamp,one,al,xja,xj,xjp)
     &          * xinjot(half,half,one,xjp,xja,one,xlp,xla,one)
     &          * redume(xj,alamp,xjp)
        sumder = sumder+(1-2*mod(lpcase,2))*add
        sumrad = sumrad-add*(xla+(lpcase-1))
       endif
   20  continue
      endif
   10 continue
c
      fac    = sqrt((al+half)*(xja+half)*72.0  *2.0  )
      facder = fac*sumder
      facrad = fac*sumrad
c
c     save on table
c
      if(ttable) then
       dersp(jsa,jdif,ilamp) = facder
       radsp(jsa,jdif,ilamp) = facrad
      else
       write(6,'(a)') ' non-standard call to "redspo"'
      endif
c
      return
      end
c-----redsig------------------------------------------------------
c
      function redsig(xj,xl,al,alamp,xja,xla)
c
      include 'params.inc'
c
c     reduced matrix element for the spin part
c      ( ja,la | t(lrpa,lamp)*sigma | j,l )
c     in edmonds notation ; t(l,ld) is a vector harmonics
c
c     attention: the quantum numbers are not checked for consistency.
c                this needs to be done by the calling routine.
c
      include 'redmts.inc'
      include 'linres.inc'
c
      data zero,half,one,xlarg/0D0,0.5D0,1.0D0,1.0D20/
c
c     check table
c
      jsa   = xja+xla+half
      jdif  = xj-xja+lrpa+1
      ttable = jsa.le.kjrpmx .and. jdif.le.klrpa .and. al.eq.alrpa
     &         .and. alamp.eq.al
      if(ttable) then
       x     = redsg(jsa,jdif)
       if(x.lt.xlarg) then
        redsig = x
        return
       endif
      endif
c
c     no table entry found. compute from scratch
c
c     alamp  = al
      iphase = xl
      iphase = (1-mod(iphase,2)*2)
      redsig = iphase*xinjot(half,half,one,xj,xja,al,xl,xla,alamp)
     &        * wigner(xl,alamp,xla,zero,zero,zero)
     &        * sqrt( (xj+half)*(xja+half)*(al+half)
     &               *(xl+half)*(xla+half)*(alamp+half)*96.0/3.1415926 )
c
c     save on table
c
      if(ttable) then
       redsg(jsa,jdif) = redsig
      else
       write(6,'(a)') ' non-standard call to "redsig"'
      endif
c
      return
      end
c-----setred------------------------------------------------------
c
      subroutine setred
c
      include 'params.inc'
c
c     sets the storage fields for the reduced matrix elements
c     to a very large number in order to indicate that they
c     are not yet initialized properly.
c
      include 'redmts.inc'
c
      data xlarg/1.0D20/
c
      xlargp = xlarg*(1.0001)
      do jdif=1,klrpa
      do js=1,kjrpmx
        redmt(js,jdif,3) = xlargp
        redsg(js,jdif) = xlargp
        do ilamp=1,2
          redmt(js,jdif,ilamp) = xlargp
          dercr(js,jdif,ilamp) = xlargp
          radcr(js,jdif,ilamp) = xlargp
          dersp(js,jdif,ilamp) = xlargp
          radsp(js,jdif,ilamp) = xlargp
      enddo
      enddo
      enddo
c
      return
      end
c-----iniobs------------------------------------------------------
c
      subroutine iniobs()
c
      include 'params.inc'
c
c     initialises Bessel function as observable 
c     presently only for LRPA=1
c     subtraction term is stored in last entry of array for observable
c
      include 'grid.inc'
      include 'wavfun.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
c     check
c
      if(lrpa.ne.1 .and. qobserv.le.xsmall)
     &   stop  ' QOBSERV=0.0 requires LRPA=1 '
c
c     compute Bessel operator
c
      qobs2     = 30D0/qobserv**3                     ! to obtain r**3
      qradob(1) = zero
      if(qobserv.gt.xsmall) then
        do i=2,ngrid
          rq      = (i-1)*rstep*qobserv
ctest          qradob(i) = (bessel(rq,lrpa)-rq/3.0d0)*qobs2
          qradob(i) = bessel(rq,lrpa)
        enddo
        if(lrpa.eq.1) call derv17(qradob(1),aux1(1),-one,ngrid,rstep)
        write(14,'(a,50(/1x,5g12.4))')
     &     ' QRADOB:',(qradob(i),i=1,ngrid)
      else
        do i=2,ngrid
          r      = (i-1)*rstep
          qradob(i) = r*(r*r)
          aux1(i)   = 3.0D0*(r*r)
        enddo
      endif
c
      if(lrpa.ne.1 .OR. KCMCORR.NE.1) return
c
c     accumulate ground state expectation value for c.m. subtraction
c
      sumnor = 0D0
      sumbes = 0D0
      do i=2,ngrid
        r      = (i-1)*rstep
        vol    = r*r*(rho(i)+rho(i+ngrid))
        sumnor = vol + sumnor
        sumbes = vol*(2.0*qradob(i)/r+aux1(i)) + sumbes
      enddo
      qradob(ngridp) = sumbes/sumnor/3.0
c
      do ifile=6,14,8
        write(ifile,'(2(a,g12.4),g12.4)') 
     &    ' Bessel observable initialised for q=',qobserv,
     &    ' , subtractive term is ',qradob(ngridp),rms*rms*5.0/3.0
      enddo
      do i=2,ngrid
        r      = (i-1)*rstep
        qradob(i) = qradob(i)-r*qradob(ngridp)
      enddo
c
      return
      end


      include 'gsc.for'
      include 'lsmooth_sub.f'
      include 'seprpa.f'
      include 'ppotnew.f'
      include 'lintdhf.f'
      include 'centroids.f'
