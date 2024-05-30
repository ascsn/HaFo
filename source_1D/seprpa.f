c-----seprpa------------------------------------------------------------------
c
      subroutine seprpa(kappa,eta,maxmod,xop,yop,ichr,resfac)
c
      include 'params.inc'
c
c     separable RPA in VPM style
c      'kappa'  = at input: inverse of X*X interaction matrix
c      'xop'    = key word for pointer on X operators 
c      'yop'    = key word for pointer on Y operators
c      'ichr'   = keys for counting states in X and Y 
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'skypar.inc'
      include 'potinv.inc'
      include 'constn.inc'
c
      double precision kappa(krpmax,krpmax),eta(krpmax,krpmax)
      double precision kapsav(krpmax,krpmax)
      dimension dxover(krpmax)
      dimension nparsv(krpmax),nholsv(krpmax),isophs(krpmax)
      character*50 toper
      character*6 xop,yop
      character*2 ichr(kichr)
      complex*16  comega,cdet
      complex*16  resp(kg1)
      complex*16 grmatr(krpmax,krpmax)
      logical    qlike
c
      logical tnoodd,ttest,tresph
c      data tnoodd,ttest,tresph/.false.,.false.,.true./
      data tnoodd,ttest,tresph/.false.,.false.,.false./
c
c----------------------------------------------------------------------
c
      write(14,'(a,1pg12.4)') ' SEPRPA: resfac=',resfac
      if(tresph) then
        ifrsph = 1
      else
        ifrsph = 0
      endif
      strfac = (1.666667*rms*rms)*PI
c
c     read bounds
c
      read(5,*) omega1,omega2,delomg,widomg,isosep
      write(6,*) 'omega1,omega2,delomg,widomg,isosep=',
     &            omega1,omega2,delomg,widomg,isosep
      skapfc  = -1.0                                  ! fetch factor
      isorpa = isosep
c      write(14,*) 'omega1,omega2,delomg,widomg,isorpa=',
c     &            omega1,omega2,delomg,widomg,isorpa
c
c      write(6,'(20a)') ichr
      open(61,file='strength_full.res')
      open(62,file='strength_ph.res')
c
c     invert to obtain 'kappa'   (matrix inversion is done in place)
c
      write(14,'(a)')  ' KAPPA**-1:'
      do imode=1,maxmod
        write(14,'(1x,20g12.4)') (kappa(imod2,imode),imod2=1,maxmod)
        do imod2=1,maxmod
          kapsav(imod2,imode) = kappa(imod2,imode)
        enddo
      enddo
      call  matin2 (kappa,maxmod,krpmax,nerror,determ)
      write(14,'(a,g12.4)')
     &   ' determinant of KAPPA:',determ
      write(14,'(a)')  ' test inversion:'
      do imode=1,maxmod
        do imod2=1,maxmod
          dxover(imod2) = 0.0
          do imod3=1,maxmod
            dxover(imod2) = dxover(imod2) 
     &               + kappa(imode,imod3)*kapsav(imod3,imod2)
          enddo      
        enddo
        write(14,'(1x,20g12.4)') (dxover(imod2),imod2=1,maxmod)
      enddo
c
      write(14,'(a)')  ' KAPPA:'
      do imode=1,maxmod
        write(14,'(1x,20g12.4)') (kappa(imod2,imode),imod2=1,maxmod)
      enddo
      if(.not.tnoodd) then
        write(14,'(a)')  ' ETA**-1:'
        do imode=1,maxmod
          write(14,'(1x,20g12.4)') (eta(imod2,imode),imod2=1,maxmod)
        enddo
        call  matin2 (eta,maxmod,krpmax,nerror,determ)
        write(14,'(a,g12.4)')
     &   ' determinant of ETA:',determ
        write(14,'(a)')  ' ETA:'
        do imode=1,maxmod
          write(14,'(1x,20g12.4)') (eta(imod2,imode),imod2=1,maxmod)
        enddo
c                                        extend to super-matrix on kappa
        if(2*maxmod.gt.krpmax) stop ' SEPRPA: too many terms'
        do imode=1,maxmod
        do imod2=1,maxmod
          kappa(imod2+maxmod,imode+maxmod) = eta(imod2,imode)
ctest          kappa(imod2+maxmod,imode+maxmod) = 0.0
          kappa(imod2+maxmod,imode) = 0.0
          kappa(imod2,imode+maxmod) = 0.0
        enddo
        enddo
        maxmod = 2*maxmod
ctest
ctest        write(14,'(a)')  ' full KAPPA:'
ctest        do imode=1,maxmod
ctest          write(14,'(1x,20g12.4)') (kappa(imod2,imode),imod2=1,maxmod)
ctest        enddo
ctest
      endif
c      do imode=1,maxmod
c      do imod2=1,maxmod
c        kappa(imod2,imode) = 0.0
c      enddo
c      enddo
c
c     construct operator for the observable 
c     = multipole operator for protons
c
c      maxmod = 1 + maxmod
      iflagr = 0
c      call qinit('msc=2           ','nojj        ',0,0,0,1,
c     &            0,maxmod,0,iflagr,xkqrd,resfac,cmfac,
c     &            nparsv,nholsv,isophs,toper)
      call qinit('mom=0           ','nojj        ',0,0,0,1,
     &            0,maxmod,0,iflagr,xkqrd,resfac,cmfac,
     &            nparsv,nholsv,isophs,toper)
c
c     reduce to pure proton (neutron)
c
!      include 'initdo.inc'
!        IF(iso==1) THEN
!          qlrout(nadd+1:nadd+ngrid,jdif)=0D0
!          plrout(nadd+1:nadd+ngrid,jdif)=0D0
!        END IF
!      include 'enddo.inc'
c
      write(14,'(a,i5)') ' nr. of modes now:',maxmod
      if(maxmod.gt.krpmax) stop ' SEPRPA: too many modes'
c      write(14,'(2a)') ' XOP=',xop
      xop(5:5) = ichr(maxmod)
c      write(14,'(2a)') ' XOP=',xop
      call buffql(qlrout,'o',xop)
      do imode=1,maxmod
        kappa(maxmod,imode) = 0.0
        kappa(imode,maxmod) = 0.0
      enddo
      kappa(maxmod,maxmod) = 0.0
c      write(14,'(a,i5)') ' now MAXMOD=',maxmod
c
      if(ttest) then
        call moment(ovqq,ovqhq,ovqp,ovpp,ovpqb,ovqh0q,qqinv,widrpa,
     &              qrlpr,qrlne,1.0D0,0.000001D0,
     &              eshifq,eshift,eshifi,
     &              xop,'dummy ','dummy3',1,-1)
        sum_m3 = qpover(qlrout,plrout)                   ! M3 sumrule exact
        call buffql(qlrout,'i',xop)
        call commvq(ovpqb,ovpqb0,ovqbqb,ovqbb0,1,1.0D0)
        sum_m1 = qpover(qlrout,plrout)                   ! M1 sumrule exact
c        write(14,'(a/a/7g12.4)')
c     &  ' test observable: ',
c     &  ' ovqq,ovqhq,ovqo,ovpp,ovpqb,ovqh0q,qqinv=',
c     &    ovqq,ovqhq,ovqo,ovpp,ovpqb,ovqh0q,qqinv
c      write(14,'(a,2g12.4)')
c     &   ' <Q|Q>,<P|P>=',qpover(qlrout,qlrout),qpover(plrout,plrout)
        call buffql(qlrout,'i',xop)
        call invmrq(ovqpm1,ovpqb,ovpqb0,1.0D0,0.00000001D0)
        call buffql(qlrout,'i',xop)
        write(14,'(a,g12.4)')
     &   ' 2*<D H**-1 D>=',2.0*qpover(qlrout,plrout)
c
        call buffql(plrout,'i','dummy ')
        if(.not.tnoodd) call buffql(savlr,'i',xop)
        do imode=1,maxmod
c                                                   read X operator
          if(tnoodd .or. imode.le.maxmod/2) then
            xop(5:5) = ichr(imode)
            call buffql(qlrout,'i',xop)
            dxover(imode) = 2.0*qpover(qlrout,plrout)
          elseif(.not.tnoodd .and. imode.gt.maxmod/2 
     &           .and. imode.lt.maxmod) then
            yop(5:5) = ichr(imode-maxmod/2)
            call buffql(qlrout,'i',yop)
            dxover(imode) = 2.0*qpover(qlrout,savlr)
          endif
        enddo
        write(14,'(a,20(/1x,5g12.4))') 
     &   ' <[D,H],X>=',(dxover(imode),imode=1,maxmod)

        sum = 0.0
        do imode=1,maxmod/2
        do imod2=1,maxmod/2
          sum = dxover(imod2)*kappa(imod2,imode)*dxover(imode) 
     &         + sum
        enddo
        enddo
        write(14,'(a,2g12.4)') 
     &    ' M_3 sum rule: from kappa, exact=',sum,sum_m3
        if(.not.tnoodd) then
          sum = 0.0
          do imode=maxmod/2+1,maxmod-1
          do imod2=maxmod/2+1,maxmod-1
            sum = dxover(imod2)*kappa(imod2,imode)*dxover(imode) 
     &           + sum
          enddo
          enddo
          write(14,'(a,2g12.4)') 
     &    ' M_1 sum rule: from kappa, exact=',sum,sum_m1
        endif
        write(14,'(a)')  ' new KAPPA:'
        do imode=1,maxmod
          write(14,'(1x,20g12.4)') (kappa(imod2,imode),imod2=1,maxmod)
        enddo
      endif
c
c     loop over frequencies at which the strengths are computed
c
c     print header
c
      write(14,'(a)') 
     &  ' the strength function from separable RPA:',
     &  '  omega   real(R)   imag(R) '
      write(61,'(a)') 
     &  '# the strength function from separable RPA:',
     &  '#  omega   real(R)   imag(R) ',
     &  '6,24,0.005,0.01,40,100,',
     &  ' X N Y'
      write(62,'(a)') 
     &  '# the strength function for ph states:',
     &  '#  omega   real(R)   imag(R) ',
     &  '6,24,0.005,0.01,40,100,',
     &  ' X N Y'
      write(6,'(a)') 
     &  ' the strength function from separable RPA:',
     &  '  omega   real(R)   imag(R) '
c
c     big loop over frequencies
c
      trk = 0D0
      trkm3 = 0D0
      nomega = (omega2-omega1)/delomg
      do iomega=0,nomega
        omega  = omega1 + delomg*iomega
        comega = cmplx(omega,widomg)
c
        if(tnoodd) then
          isigmn = +1
        else
          isigmn = -1
        endif
        do imode=1,maxmod
        do isig=isigmn,1,2                              ! like ops. or crossover
c                                                   read X operator
          if(tnoodd .or. imode.le.maxmod/2 .or. 
     &       imode.eq.maxmod) then
            xop(5:5) = ichr(imode)
c            write(14,'(2a)') ' XOP=',xop
            call buffql(qlrout,'i',xop)
            qlike = .true.
          else
            yop(5:5) = ichr(imode-maxmod/2)
            call buffql(qlrout,'i',yop)
            qlike = .false.
          endif
          if(isig.lt.0) qlike = .not.qlike
c                                                   response on (PLROUT,SAVLR)
          if(tresph) then
            call respph(comega,isig)
          else
            call respwf(comega,isig)
          endif
c                                                   compute overlap with
c                                                   all other X operators
          do imod2=1,maxmod
            if(tnoodd .or. imod2.le.maxmod/2 .or. 
     &         imod2.eq.maxmod) then
              xop(5:5) = ichr(imod2)
              call buffql(qlrout,'i',xop)
              if(qlike) then
                grmatr(imod2,imode) =
     &            cmplx(qpover(qlrout,plrout),qpover(qlrout,savlr))
              endif            
            else
              yop(5:5) = ichr(imod2-maxmod/2)
              call buffql(qlrout,'i',yop)
              if(.not.qlike) then
                grmatr(imod2,imode) =
     &            cmplx(qpover(qlrout,plrout),qpover(qlrout,savlr))
              endif            
            endif
ctest            grmatr(imod2,imode) =
ctest     &        cmplx(qpover(qlrout,plrout),qpover(qlrout,savlr))
          enddo
        enddo  
        enddo
        write(62,'(1x,1pg14.6,2x,2(1pg13.5))') 
     &     omega,grmatr(maxmod,maxmod)*strfac
        CALL FLUSH(62)
c
c       now the inversion in separable space takes place
c 
ctest        write(14,'(a)') ' R_00 before inversion '
ctest        do imode=1,maxmod
ctest          write(14,'(1x,10(1pg12.4))') 
ctest     &      (grmatr(imod2,imode),imod2=1,imode)
ctest        enddo
        call cmatin2(grmatr,maxmod,krpmax,nerror, cdet)
ctest        write(14,'(a)') ' R_00 after inversion '
ctest        do imode=1,maxmod
ctest          write(14,'(1x,10(1pg12.4))') 
ctest     &      (grmatr(imod2,imode),imod2=1,imode)
ctest        enddo
        do imode=1,maxmod
        do imod2=1,maxmod
          grmatr(imod2,imode) = 
     &       grmatr(imod2,imode) - skapfc*kappa(imod2,imode)
c     &       grmatr(imod2,imode) - kappa(imod2,imode)
        enddo
        enddo
ctest        write(14,'(a)') ' R before inversion '
ctest        do imode=1,maxmod
ctest          write(14,'(1x,10(1pg12.4))') 
ctest     &      (grmatr(imod2,imode),imod2=1,imode)
ctest        enddo
        call cmatin2(grmatr,maxmod,krpmax,nerror, cdet)
ctest        write(14,'(a)') ' R after inversion '
ctest        do imode=1,maxmod
ctest          write(14,'(1x,10(1pg12.4))') 
ctest     &      (grmatr(imod2,imode),imod2=1,imode)
ctest        enddo
c                                                           print result
        write(14,'(1x,1pg13.5,2x,2(1pg13.5))') 
     &     omega,grmatr(maxmod,maxmod)*strfac
        write(61,'(1x,1pg13.5,2x,2(1pg13.5))') 
     &     omega,grmatr(maxmod,maxmod)*strfac
        write(6,'(1x,1pg13.5,2x,2(1pg13.5))') 
     &     omega,grmatr(maxmod,maxmod)*strfac
        CALL FLUSH(61)
        CALL FLUSH(14)
        CALL FLUSH(6)
        trk = trk+delomg*AIMAG(grmatr(maxmod,maxmod)*strfac)*omega
        trkm2 = trkm2
     &         +delomg*AIMAG(grmatr(maxmod,maxmod)*strfac)*omega**3
      enddo
      write(61,'(a,2(1pg13.5))') '# TRK=',trk,sqrt(trkm2/trk)
      write(62,'(a)') '#'
      write(14,'(a,2(1pg13.5))') '# TRK=',trk,sqrt(trkm2/trk)
      write(6,'(a,2(1pg13.5))') '# TRK=',trk,sqrt(trkm2/trk)
c
c
      return
      end
c-----respwf------------------------------------------------------------------
c
      subroutine respwf(comega,isig)
c
      include 'params.inc'
c
c     wavefunction of free response for complex freqency 'comega'.
c     input is operator wavefunction on 'qlrout'.
c     result is: real part on 'plrout' and imaginary part on 'savlr'
c     all via 'linres.inc'
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'skypar.inc'
      include 'potinv.inc'
      include 'constn.inc'
      complex*16  comega
      complex*16  resp(kg1)
c
c-----------------------------------------------------------------------
c
c                                                   loop over states
c                                                   compute response
          include 'initdo.inc'
c          redmat = redume(xj,alrpa,xja)
          do i=1,ngrid
            sqb(i) = sqrt(bmass(i+isoadd))
          enddo
          call derv25(sqb(1),d2b(1),one,ngrid,rstep)
          do i=1,ngrid
            d2b(i) = d2b(i)/sqb(i)
          enddo
          xlp    = (ll*ll+ll)
          xl0    = wls(js)
          do ir=1,ngrid
            irs    = ir+isoadd
            potl(ir) = pot(irs)+xl0*potls(irs)
     &          +bmass(irs)*xlp*zentrf(ir)+d1b(irs)*rinvs(ir)
c            aux1(ir) = (ir*rstep-rstep)*qlrout(nadd+ir,jdif)
            aux1(ir) = rfact(ir)*qlrout(nadd+ir,jdif)
          enddo
c                                                   invert and
c                                                   save real on PLROUT
c                                                   and imaginary on SAVLR
          call cinver(resp,aux1,comega+cmplx(subfac,zero))
          do ir=1,ngrid
            plrout(nadd+ir,jdif) = resp(ir)
            savlr(nadd+ir,jdif) = dimag(resp(ir))
          enddo
          call cinver(resp,aux1,-comega+cmplx(subfac,zero))
          signum = isig*ONE
          fac  = wgt
          do ir=2,ngrid
            plrout(nadd+ir,jdif) = 
     &           fac*rinvs(ir)*(signum*resp(ir) + plrout(nadd+ir,jdif))
            savlr(nadd+ir,jdif) = 
     &           fac*rinvs(ir)*(signum*dimag(resp(ir))
     &                            +savlr(nadd+ir,jdif))
          enddo
          parll  = (1-mod(ll,2)*2)
          call extrap(plrout(nadd+1,jdif),parll)
          call extrap( savlr(nadd+1,jdif),parll)
          include 'enddo.inc'
c
      return
      end

c-----commvp--------------------------------------------------------------------
c
      subroutine commvp(ovpq,ovpq0,ovqq,ovq0q0,resfac,cmfac,ifsav0)
c
      include 'params.inc'
c
c     computes q = [v_res,p] for given p
c     with p and q given by its action on wave functions:
c         p |wf> <--> 'plrout'
c         q |wf> <--> 'qlrout'
c
c     the response wave functions are handled via common /linres/.
c     they mean:
c         'plrout' = given p = input to routine
c         'qlrout' = resulting [h,p] = output of routine
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
      write(14,'(a,1pg12.4)') ' COMMVP: resfac=',resfac
      call pdensi(resfac)
c                                     folding, if damping is desired
      if(ayuk.ne.zero) then
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
       if(t4.ne.zero) then
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
       call filter(prho(1),prho(1),lrpa)
cc       write(6,'(a)') ' repeat filter:'
cc       call filter(prho(1),prho(1),lrpa)
       call filter(prho(ngrid+1),prho(ngrid+1),lrpa)
       call filter(ptau(1),ptau(1),lrpa)
       call filter(ptau(ngrid+1),ptau(ngrid+1),lrpa)
      endif
c
c     prepare the kinetic response potential for given current
c
      call ppot(resfac,cmfac)
      if(ayuk.ne.zero) then
       call t0inv(rpot(1),rpot(1),work1,dampay,lrpa)
       call t0inv(rpot(ngrid+1),rpot(ngrid+1),work1,dampay,lrpa)
       call t0inv(rbmass(1),rbmass(1),work1,dampay,lrpa)
       call t0inv(rbmass(ngrid+1),rbmass(ngrid+1),work1,dampay,lrpa)
       if(imfold.eq.1) then
        call t0inv(rbmass(1),rbmass(1),work1,dampay,lrpa)
        call t0inv(rbmass(ngrid+1),rbmass(ngrid+1),work1,dampay,lrpa)
       endif
       call grad(rbmass,drbmas,2)
       if(t4.ne.zero) then
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
c     now  store qlrout.
c
      do ir=2,ngrid
        qlrout(nadd+ir,jdif) =  resfac * hlinrs(ir)
      enddo
      parll  = (1-mod(ll,2)*2)
      call extrap(qlrout(nadd+1,jdif),parll)
      call projec(qlrout(nadd+1,jdif),js,nstack,nmax)
c
c      write(14,'(a,2i5)')
c     &  ' XOP for na,jdif=',na,jdif
c      write(14,'(1x,5(1pg12.4))') (qlrout(nadd+ir,jdif),ir=1,ngrid)
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
      if(ifrsph.eq.1) call prjshl(1,qlrout)
c
      return
      end
c-----commvq--------------------------------------------------------------------
c
      subroutine commvq(ovqp,ovqp0,ovpp,ovp0p0,ifprj,resfac)
c
      include 'params.inc'
c
c     computes p = [v,q] for given q
c     with p and q given by its action on wave functions:
c         q |wf> <--> 'qlrout'
c         p |wf> <--> 'plrout'
c
c     the response wave functions are handled via common /linres/.
c     they mean:
c         'qlrout' = given q = input to routine
c         'plrout' = resulting [v,q] = output of routine
c
c
c     further list parameters are:
c      ifprj   = switch for ph-projection of the final [h,q]
c      resfac  = weight factor for the residual interaction;
c                should be 1.0
c
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
      data xmaxlw/1.0e-20/
c
c----------------------------------------------------------------------
c
c     prepare the response current 'qcur' and its derivative 'dqcur'
c
      write(14,'(a,1pg12.4)') ' COMMVQ: resfac=',resfac
      call qcurre(resfac)
      if(ayuk.ne.zero) then
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
      if(ayuk.ne.zero) then
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
c     action of linear response hamiltonian
c
      if(resfac.ne.zero)
     &  call hrespq(iso,la,xja,nadd,xj,xl,jdif,
     &              wf(nadd+1),wfs(nadd+1),hlinrs)
cc
c
c     compute [v,q] from hlinrs
c
      do ir=2,ngrid
        plrout(nadd+ir,jdif) = resfac * hlinrs(ir)
      enddo
      parll  = (1-mod(ll,2)*2)
      call extrap(plrout(nadd+1,jdif),parll)
cc
      if(ifprj.eq.1)
     & call projec(plrout(nadd+1,jdif),js,nstack,nmax)
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
      if(ifrsph.eq.1) call prjshl(1,plrout)
c
      return
      end
c-----respph-------------------------------------------------------------------
c
      subroutine respph(comega,isig)
c
      include 'params.inc'
c
c     free response function at 'comega' 
c     from explicit ph states in window 'ecutlw'...'ecutup'.
c     the input wave function comes on 'wfin',
c     the result is returned in 'resp'.
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'result.inc'
      include 'skypar.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
c
c      complex*16  resp(kg1)
      complex*16  comega,respwf,cfac,csumph

      data xminrd/1.0e-5/
      data ecutlw,ecutup/0,-30/                         ! bounds for ph states
      data iftspr/1/
c
c--------------------------------------------------------------------
c
      if(ecutup.lt.ZERO) then
        read(5,*) ecutlw,ecutup
        write(14,'(a,2g12.4)') 
     &   ' PH response with ecutlw,ecutup=',ecutlw,ecutup
      endif
      signum = isig*ONE
c
      do jdif=1,klrpa
      do ir=1,kx11
        plrout(ir,jdif) = 0.0
        savlr(ir,jdif)  = 0.0
      enddo
      enddo
      csumph = 0.0
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
c      write(6,'(a,1i5,2g12.4)')
c     &    ' na,wgt,redmat=',na,wgt,redmat
c
c     loop over the b-states
c
      do 40 nb=1,nmax
      jsb   = jsp(nb)
      if(jsb.ne.js) goto 40
      lb    = jsb/2
      xlb   = lb
      xjb   = (jsb-1)/2+0.5
      wgtb  = weight(nb+nstack/ngrid)
c      write(6,'(a,2i5,2g12.4)')
c     &    ' na,nb,wgtb=',na,nb,wgtb
      if(wgtb.gt.half) goto 40
      naddb = nstack + (nb-1)*ngrid
      subfcb = spenrg(nb+nstack/ngrid)
c
c     assign quantum numbers in pointer fields 'iniqr', 'iniiso' and
c     'xkqr'.
c
      edif = abs(subfac-subfcb)
c      write(6,'(a,2i5,2g12.4)')
c     &    ' na,nb,edif=',na,nb,edif
      if(edif.lt.ecutup .and. edif.gt.ecutlw
     &   .and. subfac.ne.zero .and. subfcb.ne.zero) then
        sum = 0.0
        do ir=2,ngrid
          sum = sum + rfact(ir)*wf(naddb+ir)*qlrout(nadd+ir,jdif)
        enddo
        sum  = rstep*sum
        if(iftspr.eq.1) write(14,'(1x,1g12.4)')   edif,sum
        cfac = -(sum/(comega-cmplx(edif,zero))
     &          +signum*sum/(-comega-cmplx(edif,zero)) )
        do ir=2,ngrid
          respwf = wf(naddb+ir)*rinvs(ir)*cfac
          plrout(nadd+ir,jdif) = respwf + plrout(nadd+ir,jdif) 
          savlr(nadd+ir,jdif)  = imag(respwf) + savlr(nadd+ir,jdif) 
        enddo
        csumph = sum*cfac + csumph
        csumsb = sum*cfac + csumsb
        parll  = (1-mod(ll,2)*2)
        call extrap(plrout(nadd+1,jdif),parll)
        call extrap( savlr(nadd+1,jdif),parll)
c        write(6,'(a,3i5,4g12.4)')
c     &    ' na,nb,jdif,edif,sum=',na,nb,jdif,edif,sum,
c     &    wfnorm(wf(naddb+1)),wfnorm(qlrout(nadd+1,jdif))
      endif
c
   40 continue
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
      iftspr = 0
c
      return
      end

c-----cinver------------ part of the hartree-fock package ----------------------
c
      subroutine cinver(r,rhs,cumin)
c
      include 'params.inc'
c
      complex*16 r(kg1)
      complex*16 cainv(kg1),cbinv(kg1),cpotef(kg1)
      complex*16 cumin,x,y,rim1,ui,z
      dimension  rhs(kg1)
c
      include 'potinv.inc'
      include 'grid.inc'
      include 'constn.inc'
c
       data fsx,ezw/.83333333333333,.8333333333333e-1/
c
c
c      determines 'r' as solution of the inhomogenous linear diff. eq.:
c       " (-d/dx (sqb(x)**2) d/dx + (potl(x) - cumin ) ) r = rhs  "
c      for given 'rhs' where 'sqb**2' is x-dependent inverse mass,
c      'potl' a potential and 'cumin' a shift constant.
c
c      the arrays in the calling list mean:
c       r      = the solution of the equation
c       rhs    = the right hand side (inhomogeneity) of the equation.
c                r and rhs may be the same arrays in the calling routine
c                rhs is then overwritten.
c      the parameter in the calling list:
c       cumin   = constant energy shift
c
c      information on the hamiltonian is transferred via common /potinv/:
c       sqb    = inverse square-root of x-dependent effective mass
c       potl   = effective potential before mass transformation,
c                to be prepared in subroutine 'welfu'.
c
c      information on the grid is transferred via common /grid/.
c
c      note: the program is written for solving the schroedinger equatio
c            in radially symmetric problems. thus it implies the rhs and
c            cpotef to be such that r(x) vanishes at the origin.
c
c      the solution proceeds as follows:
c
c      by introducing scaled wave functions
c       " reff = r/sqb "
c      and
c       " rhseff = rhs/sqb ".
c      the equations is transformed to a laplacian equation
c       " (-d/dx**2 + cpotef(x) ) reff = rhseff "
c      with constant mass and an effective potential
c       " cpotef = (potl(x)-cumin)/(sqb(x)*sqb(x)) + d2b(x) "
c      where
c       " d2b = (d/dx**2 sqb)/sqb " .
c      after solving the laplacian equation, the solution is transformed
c      back by
c       " r   = reff*sqb" .
c
c      the d/dx**2 is treated effectively in 5-point precision.
c      in addition to a standard 3-point formula for d/dx**2
c      a complementing 3-point averaging of ueff*reff and rhseff
c      is used which both together make up the 5-point precision.
c      (see r.zurmuehl, "praktische mathematik", ch. vii, sect. 3,
c       springer verlag, berlin 1963)
c
c
c-----------------------------------------------------------------------
c
c      square-root of inverse eff. mass on 'sqb',
c      second derivative of 'sqb' on 'd2b'.
c      finally assembling  the effective potential on 'cpotef'.
c       this part can be overridden if the routine is called for
c       the same 'potl' and 'sqb' again and if the workspaces
c       'sqm' and 'cpotef' are not destroyed in the meantime;
c       this situation should be indicated by 'irep>0' in the call.
c
c      if(irep.le.0) then
       sqm(1) = one/sqb(1)
        do 80 i=2,ngrid
        sqm(i)    = one/sqb(i)
        cpotef(i) = d2b(i)+(potl(i)-cumin)*(sqm(i)*sqm(i))
   80   continue
c      endif
c
c      transform to effective r.h.s., store result on r.
c
       do 17 i=1,ngrid
   17  r(i)   = rhs(i)*sqm(i)
c
c     prepare parameters for kinetic energy and loop-ends.
c
        offdig = one/(rstep*rstep)
        offdg2 = offdig+offdig
        nmaxm  = ngrid-1
c
c      check repetition switch. if new potl,sqm are indicated be irep=0
c      then begin first pass eliminating the lower off-diagonal element.
c      new diagonal and upper off-diagonal elemnts are saved on a and b.
c      the second pass, forward elimination of the r.h.s., is done withi
c      the same loop; result is stored on field t.
c
c      if(irep.le.0) then
c
        cbinv(1) = zero
        cainv(1) = one
        x      = zero
        y      = one
        rim1   = r(1)
c
        do 13 i=2,nmaxm
c                                         forward elim. for r.h.s.
        x      = ezw*(r(i+1)+rim1)+fsx*r(i)- x*(cbinv(i-1)/y)
        rim1   = r(i)
        r(i)   = x
c                                         forward elim. for matrix
        ui     = cpotef(i)
        z      = -offdig+ezw*ui
        cbinv(i) = z
        y      = offdg2+fsx*ui-z*(cbinv(i-1)/y)
        cainv(i) = y
   13   continue
c                                         last elim. step explicitely
c                                         because r(i+1) is not availabl
        r(ngrid) = ezw*rim1+fsx*r(ngrid)-x*(cbinv(nmaxm)/y)
        ui     = cpotef(ngrid)
        z      = -offdig+ezw*ui
        cbinv(ngrid) = z
        y      = offdg2+fsx*ui-z*(cbinv(nmaxm)/y)
        cainv(ngrid) = y
c
c      if unchanged potl,sqm are indicated by irep=1 then perform
c      in the second pass the forward elimination for the r.h.s.
c      using the diagonal and off-diagonal elements of the previous
c      call.
c
c       else
c        x      = zero
c        rim1   = r(1)
c        do 14 i=2,nmaxm
c        x      = ezw*(r(i+1)+rim1)+fsx*r(i)- x*cbinv(i-1)/cainv(i-1)
c        rim1   = r(i)
c        r(i)   = x
c   14   continue
c        r(ngrid)  = (ezw*rim1+fsx*r(ngrid)-x*cbinv(nmaxm)/cainv(nmaxm))
c       endif
c
c      in the third pass resolve the remaining upper triangular equation
c      backward.
c
       x      = r(ngrid)/cainv(ngrid)
       r(ngrid) = x*sqm(ngrid)
       do 15 i=nmaxm,2,(-1)
       x      = (r(i)-x*cbinv(i+1))/cainv(i)
       r(i)   = x*sqm(i)
   15  continue
       r(1)   = zero
c
       return
       end
c-----cmatin2-----------------------------------------------------------
c
      subroutine  cmatin2 (a,nmax,ndim,nerror,determ)
c
      include 'params.inc'
      parameter (kindex=120)
c
c     inverts the matrix 'a'.
c
c     the parameters on the i/o list are:
c
c       a       = at input:  matrix to be inverted
c                 at output: inverse matrix
c       nmax    = number of actual rows and columns
c       ndim    = dimension of matrix in the calling routine
c       nerror  = error flag
c       determ  = determinant of input matrix 'a'
c
c       the routine is taken from 'cernlib' and modified for this
c       package.
c     complex version by Jan O. Babst 1996/07/23
c
      complex*16 a(*), determ, deter, pivot, swap
      dimension index(kindex)
c
c-----------------------------------------------------------------------------
c
      if(nmax.gt.kindex) stop ' too large dimension in matin2'
c
      deter=cmplx(1.00, 0.0)
      n=nmax
c     iemat=n+n2
      iemat=n+0
      idim=nmax
      nmin1=n-1
c
c   umsetzen des nmax*nmax-anteils der matrix a in gepackte
c   ein-dimensionale form.
c
      do 400 i=2,idim
      do 400 j=1,idim
      ipack  = (i-1)*nmax + j
      imatr  = (i-1)*ndim + j
      a(ipack) = a(imatr)
  400 continue
c        the routine does its own evaluation for double subscripting of
c        array a.
      ipivc=1-idim
c        main loop to invert the matrix
      do 11 main=1,n
      pivot=cmplx(0.00, 0.00)
      ipivc=ipivc+idim
c        search for next pivot in column main.
      ipivc1=ipivc+main-1
      ipivc2=ipivc +nmin1
      do 2 i1=ipivc1,ipivc2
      if(abs(a(i1))-abs(pivot)) 2,2,1
    1 pivot=a(i1)
      lpiv=i1
    2 continue
c        is pivot different from zero
      if(abs(pivot)) 3,15,3
c        get the pivot-line indicator and swap lines if necessary
    3 icol=lpiv-ipivc+1
      index(main)=icol
      if(icol-main) 6,6,4
c        complement the determinant
    4 deter=-deter
c        pointer to line pivot found
      icol=icol-idim
c        pointer to exact pivot line
      i3=main-idim
      do 5 i=1,iemat
      icol=icol+idim
      i3=i3+idim
      swap=a(i3)
      a(i3)=a(icol)
    5 a(icol)=swap
c        compute determinant
    6 deter=deter*pivot
      pivot=1./pivot
c        transform pivot column
      i3=ipivc+nmin1
      do 7 i=ipivc,i3
    7 a(i)=-a(i)*pivot
      a(ipivc1)=pivot
c        pivot element transformed
c
c        now convert rest of the matrix
      i1=main-idim
c        pointer to pivot line elements
      icol=1-idim
c        general column pointer
      do 10 i=1,iemat
      icol=icol+idim
      i1=i1+idim
c        pointers moved
      if(i-main) 8,10,8
c        pivot column excluded
    8 jcol=icol+nmin1
      swap=a(i1)
      i3=ipivc-1
      do 9 i2=icol,jcol
      i3=i3+1
    9 a(i2)=a(i2)+swap*a(i3)
      a(i1)=swap*pivot
   10 continue
   11 continue
c        now rearrange the matrix to get right invers
      do 14 i1=1,n
      main=n+1-i1
      lpiv=index(main)
      if(lpiv-main) 12,14,12
   12 icol=(lpiv-1)*idim+1
      jcol=icol+nmin1
      ipivc=(main-1)*idim+1-icol
      do 13 i2=icol,jcol
      i3=i2+ipivc
      swap=a(i2)
      a(i2)=a(i3)
   13 a(i3)=swap
   14 continue
      determ=deter
      nerror=0
c   ruecksetzen der gepackten eindimensionalen form in den
c   nmax*nmax-anteil der matrix c
c
      do 500 i=idim,2,-1
      do 500 j=idim,1,-1
      ipack  = (i-1)*nmax + j
      imatr  = (i-1)*ndim + j
  500 a(imatr) = a(ipack)
c
      return
   15 nerror=main
      determ=deter
      return
      end
