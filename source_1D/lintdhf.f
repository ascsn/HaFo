c-----lintdhf----------------------------------------------------
c
      subroutine lintdhf
c
      include 'params.inc'
c
c     Linearized TDHF using the routines of the RPA package.
c     The Q-like and P-like excitation operators are thought
c     to be combined to a generator like operator C^+ = Q + i P. 
c     The Q and P are handled on the arrays QLROUT and PLROUT
c     from LINRES.INC, as usual.
c     The C^+ is propagated by a predictor-corrector of second
c     order.
c     The arrays PROP0 and PROP1 serve as storage for the
c     stepping operators.
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
      dimension nparsv(krpmax),nholsv(krpmax),isophs(krpmax)
      dimension prop0(KTIMWF,KLRPA,2),prop1(KTIMWF,KLRPA,2)
      dimension state(KTIMWF,KLRPA,2),observ(KTIMWF,KLRPA)
      character*16 typmod
      character*12 jselec
      character*50 toper
      logical   top14,tdamp
c
c---------------------------------------------------------------------
c
      if(KTIMWF.ne.KX11) stop ' KTIMWF=KX11 needed for LINTDHF'
c
      inquire(14,opened=top14)
      if(top14) close(14)
      open(unit=14,status='unknown',form='formatted',file='for014')
      resdmp = 0.0d0      
      dampenrg = 0.0D0
      nabsorb = 0
      write(6,'(a)') ' now reading basic rpa:'
      read(5,*) lrpa,resfac,cmfac,resdmp,ifprj,dampenrg
      write(6,*) lrpa,resfac,cmfac,resdmp,ifprj,dampenrg
      write(6,'(a)') ' now reading timestep:'
      read(5,*) timestep,numsteps,modpri,nabsorb,powabsorb
      write(6,*) timestep,numsteps, modpri,nabsorb,powabsorb
      write(6,'(a)') ' now reading typmod'
      read(5,*) typmod,isorpa,jselec,isplit,ihwgt,ih0wgt,ipwgt,iprshl
      write(6,'(1x,a,i5,a,5i5)')
     &    typmod,isorpa,jselec,isplit,ihwgt,ih0wgt,ipwgt,iprshl
      tdamp = dampenrg.ne.0.0D0
      if(tdamp) dampfac = 0.5D0*(amfc(1)+amfc(2))/dampenrg
c
c     overrule the expensive 'filter' in q-space by simple
c     kinetic-energy daming in q-space
c
      if(resdmp.ne.0.0) then
        ayuk = 1.0D0/resdmp 
        resdmp = 0.0
      endif
 
c
c     supply derivative of 'bmass' on 'd1b' for use in 'hpsi'
c
      call derv17(bmass(1),d1b(1),one,ngrid,rstep)
      call derv17(bmass(ngrid+1),d1b(ngrid+1),one,ngrid,rstep)
c
c     provide derivative of "rho"
c
      call derv17(rho(1),d1r(1), one,ngrid,rstep)
      call derv17(rho(ngrid+1),d1r(ngrid+1), one,ngrid,rstep)
c
c     initialize QLROUT and PLROUT
c     save QLROUT on SAVLR as observable
c
      alrpa = lrpa*1.0D0
      mlrpa=2*lrpa+1
      call setred
      ifanlr = 1
      iwr    = 0
      do i=1,KTIMWF
      do j=1,KLRPA
        state(i,j,2) = 0.0D0
      enddo
      enddo
      call qinit(typmod,jselec,isplit,ihwgt,ih0wgt,ipwgt,
     &           iprshl,imode,ifprim,iflagr,xkqrd,resfac,cmfac,
     &           nparsv,nholsv,isophs,toper)
      write(6,'(a,1pg13.5)')
     &  ' test initial Q:  <Q|Q>=',
     &  qpover(qlrout,qlrout)
      call  copywf(qlrout,state(1,1,1))
      call  copywf(qlrout,observ)
c
c     time stepping and prints
c
      if(modpri.eq.1) then
        call copywf(state(1,1,1),qlrout)
        call commhq(ovqp,ovqp0,ovpp,ovp0p0,ifprj,resfac,0)
        write(6,'(3(a,1pg13.5))') 
     &     ' <Q|Q>=',qpover(state(1,1,1),state(1,1,1)),
     &     ' <P|HQ>=',qpover(state(1,1,1),plrout),
     &     ' <P|P>=',qpover(plrout,plrout)
        call copywf(state(1,1,1),plrout)
        call commhp(ovqp,ovqp0,ovpp,ovp0p0,resfac,cmfac,0)
        write(6,'(3(a,1pg13.5))') 
     &     ' <P|P>=',qpover(state(1,1,1),state(1,1,1)),
     &     ' <Q|HP>=',qpover(state(1,1,1),qlrout),
     &     ' <Q|Q>=',qpover(qlrout,qlrout)
      endif
      write(14,'(a)')
     &  '#  time   multipole-moment  1-norm   energy',
     &  '  0.0   0.0  0.0   0.0'
      do itime=1,numsteps
        time = itime*timestep
c              predictor
        call copywf(state(1,1,1),qlrout)
c        if(modpri.eq.1) write(6,'(a,1pg13.5)') 
c     &    ' 0: <Q|Q>=',qpover(qlrout,qlrout)
        call commhq(ovqp,ovqp0,ovpp,ovp0p0,ifprj,resfac,0)
        if(tdamp) call dampresp(plrout,dampfac)
        call copywf(plrout,prop0(1,1,2))
        qhq = qpover(qlrout,plrout)
c        if(modpri.eq.1)  write(6,'(2(a,1pg13.5))') 
c     &     ' <P|P>=',qpover(plrout,plrout)
c        if(modpri.eq.1) then
c          call commhp(ovqp,ovqp0,ovpp,ovp0p0,resfac,cmfac,0)
c        endif
        call copywf(state(1,1,2),plrout)
        pwid = qpover(plrout,plrout)
c        if(modpri.eq.1) write(6,'(a,1pg13.5)') ' <P|P>=',pwid
        if(pwid.lt.1.0D-20) then
          call copywf(plrout,qlrout)
        else
          call commhp(ovqp,ovqp0,ovpp,ovp0p0,resfac,cmfac,0)
        if(tdamp) call dampresp(qlrout,dampfac)
        endif
        php = qpover(plrout,qlrout)
c        if(modpri.eq.1)  write(6,'(2(a,1pg13.5))') 
c     &    ' <Q|Q>=',qpover(qlrout,qlrout)
c        if(mod(itime,modpri).eq.0) write(6,'(a,3(1pg12.4))')
c     &   ' qhq,php,energy2=',qhq,php,qhq+php
        call copywf(qlrout,prop0(1,1,1))
        call admlwf(state(1,1,1),prop0(1,1,1),1.0D0,timestep)
        call admlwf(state(1,1,2),prop0(1,1,2),1.0D0,-timestep)
        if(modpri.eq.1) write(6,'(a,1pg13.5)') 
     &    ' new=',qpover(state(1,1,2),state(1,1,2))
        signal = qpover(state(1,1,2),observ)
        if(modpri.eq.1) 
     &   write(6,'(a,f10.4,1pg13.5)') '0: time,signal=',time,signal
c              corrector
        call copywf(state(1,1,1),qlrout)
        if(modpri.eq.1) write(6,'(a,1pg13.5)') 
     &    ' 1: <Q|Q>=',qpover(qlrout,qlrout)
        call commhq(ovqp,ovqp0,ovpp,ovp0p0,ifprj,resfac,0)
        if(tdamp) call dampresp(plrout,dampfac)
        if(modpri.eq.1) write(6,'(a,1pg13.5)') 
     &    '  <P|P>=',qpover(plrout,plrout)
        call copywf(plrout,prop1(1,1,2))
        call copywf(state(1,1,2),plrout)
        if(modpri.eq.1) write(6,'(a,1pg13.5)') 
     &    '  <P|P>=',qpover(plrout,plrout)
        call commhp(ovqp,ovqp0,ovpp,ovp0p0,resfac,cmfac,0)
        if(tdamp) call dampresp(qlrout,dampfac)
        if(modpri.eq.1) write(6,'(a,1pg13.5)') 
     &    ' <Q|Q>=',qpover(qlrout,qlrout)
        call copywf(qlrout,prop1(1,1,1))
        call add3wf(state(1,1,1),prop0(1,1,1),prop1(1,1,1),
     &              -timestep*0.5D0,timestep*0.5D0)
        call add3wf(state(1,1,2),prop0(1,1,2),prop1(1,1,2),
     &              timestep*0.5D0,-timestep*0.5D0)
c              optionally absprbing bounds
        if(nabsorb.gt.0) then
          call absoresp(state(1,1,1),nabsorb,powabsorb)
          call absoresp(state(1,1,2),nabsorb,powabsorb)
        endif
c              measure and print        
        if(mod(itime,modpri).eq.0) then
          signal = qpover(state(1,1,2),observ)
          xnorm  = qpover(state(1,1,1),state(1,1,1))
     &            + qpover(state(1,1,2),state(1,1,2))
          write(6,'(a,f10.4,5(1pg13.5))') 
     &        'time,signal,norm,qhq,php,energy2=',
     &         time,signal,xnorm-1.0D0,qhq,php,qhq+php
          write(14,'(f10.4,3(1pg13.5))') 
     &         time,signal,xnorm-1.0D0,qhq+php
        endif
      enddo      
c
      return
      end
c-----add3wf------------------------------------------------------
c
      subroutine add3wf(awf,bwf,cwf,bfac,cfac)
c
      include 'params.inc'
c
c     add and multiply operation
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      dimension awf(kx11,klrpa),bwf(kx11,klrpa),cwf(kx11,klrpa)
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
        awf(nadd+ir,jdif) = awf(nadd+ir,jdif)
     &                     +bfac*bwf(nadd+ir,jdif)
     &                     +cfac*cwf(nadd+ir,jdif)
      enddo
c
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
      return
      end
c-----dampresp------------------------------------------------------
c
      subroutine dampresp(awf,dampfac)
c
      include 'params.inc'
c
c     kinetic energy damping of response wavefunction 'awf'
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      dimension awf(kx11,klrpa)
c
c--------------------------------------------------------------------
c
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     now all quantum numbers are provided
c
      work1(1) = 0.0D0
      do ir=2,ngrid
        work1(ir) = awf(nadd+ir,jdif)*(ir-1)*rstep
      enddo
      call t0inv(work1,work1,work2,dampfac,la)
      do ir=2,ngrid
        awf(nadd+ir,jdif) = work1(ir)*rinvs(ir)
      enddo
      parity = 1-2*mod(la,2)
      call extrap(awf(nadd+1,jdif),parity)
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
      return
      end
c-----absoresp------------------------------------------------------
c
      subroutine absoresp(awf,nabsorb,powabsorb)
c
      include 'params.inc'
c
c     absorbing bounds for response wavefunction 'awf'
c      nabsorb    = number of absorbing points (within grid)
c      powabsorb  = power of absorbing term
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      dimension awf(kx11,klrpa)
c
c--------------------------------------------------------------------
c
c     prepare mask function
c
      nabs1 = ngrid-nabsorb
      do ir=nabs1,ngrid-1
        work1(ir) = cos(pi/2*(ir-nabs1)/(nabsorb))**powabsorb
      enddo
      work1(ngrid) = 0.0D0
c
c     provide necessary quantum numbers and pointers with
c     do loops 10, 20 and 30
c
      include 'initdo.inc'
c
c     now all quantum numbers are provided
c
      do ir=nabs1,ngrid
        awf(nadd+ir,jdif) = awf(nadd+ir,jdif)*work1(ir)
      enddo
c
c     close the loops 10, 20 and 30
c
      include 'enddo.inc'
c
      return
      end
