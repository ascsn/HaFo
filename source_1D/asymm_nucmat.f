c-----asymm_nucmat-----------------------------------------------------
c
      subroutine asymm_nucmat(rhonm,iwunit)
c
      include 'params.inc'
c
c     computes asymmetric nuclear matter for given equiblibrium
c     density 'rhonm'. Results are printed on output unit 'iwunit'.
c
      real*8 energm(-50:50),d2energ(-49:49)
      data delprocon /0.01D0/
c
c----------------------------------------------------------------------
c
      WRITE(39,'(1x/1x/a,1pg13.5)') '# asymm. matter for rhonm=',rhonm
c
c     compute energies
c
      DO i=-50,50
        procon=i*delprocon
        call enernm(rhonm,0.5D0+procon,zero,energm(i),efmasp,efmasn)
      END DO
c
c     compute 2. derivatives
c
      delp2 = (2D0*delprocon)**2
      DO i=-49,49
        d2energ(i) = (energm(i-1)-2D0*energm(i)+energm(i+1))/delp2
      END DO
c
c     print
c
      DO i=-48,48
        procon=i*delprocon*2D0
        d4energ = (d2energ(i-1)-2D0*d2energ(i)+d2energ(i+1))/delp2
        IF(i==0) THEN
          WRITE(39,'(i5,f6.2,3(1pg13.5),a)') 
     &       i,procon,energm(i),d2energ(i)/2D0,d4energ/24D0,'    .  '
        ELSE
          WRITE(39,'(i5,f6.2,5(1pg13.5),a)') 
     &       i,procon,energm(i),d2energ(i)/2D0,d4energ/24D0,
     &       (energm(i)-energm(0))/procon**2
        END IF
      END DO
c
      RETURN
      END


c-----ASYM_EQUI------- part of the Hartree-Fock package ------------------
c
      subroutine asym_equi()
c
c     Evaluates the stationary state for asymmetric nuclear matter from
c     minimising the energy as function of the density.
c     This is done for a series of proton contents 'procon'
c     and printed to file unit=2
c
c     The force parameters are communicated via common/skypar/.
c
c     The list parameters are:
c       procon  = relative proton content
c       rhonm   = density
c       eovera  = energy per particle
c
      include 'params.inc'
      include 'skypar.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'wavfun.inc'
      include 'constn.inc'
c
      data rhoini,drhoit,itrmax
     *    / 0.16 , 0.0001 ,    200/
      data deltpc/0.01/            ! variation of proton content
      data ifprnt/0/
c
c-----------------------------------------------------------------------
c
      WRITE(2,'(a)') 'Equilibrium asymmetric matter:',
     &               ' procon   rho     E/A  '
      DO ipro=50,1,-2
        procon = ipro/100D0
c
c     iteration to find the minimum energy
c
      eold   = 1.e10
      rhonm0 = rhoini
      do iter=1,itrmax
        rhonmp = rhonm0+drhoit
        call enernm(rhonmp,procon,zero,eplus,efmasp,efmasn)
        call enernm(rhonm0,procon,zero,energy,efmasp,efmasn)
        rhonmm = rhonm0-drhoit
        call enernm(rhonmm,procon,zero,eminus,efmasp,efmasn)
        rhonm0 = -0.5*drhoit*(eplus-eminus)/(eplus-energy-energy+eminus)
     *           +rhonm0
        if(ifprnt.ge.0)
     &   write(06,'(a,i3,a,2g11.4,a,3g11.4)') ' it=',iter,' rho=',
     &      rhonmm+drhoit,rhonm0,' e+0-=',eplus,energy,eminus
        if(abs(eold-energy).lt.dabe) goto 29
        eold   = energy
      enddo
      write(06,'(a,i5,a)') ' no convergency after ',itrmax,
     *            ' steps in subroutine nucmat !'
      write(02,'(a,i5,a)') ' no convergency after ',itrmax,
     *            ' steps in subroutine nucmat !'
c
   29 continue
c
      WRITE(2,'(f7.2,f8.3,1pg13.5)') procon,rhonm0,energy
      END DO
c
      RETURN
      END
