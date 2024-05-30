      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   !
      data betab /5.742467999/     ! (3/5)*(3*pi^2)^*(2/3)  factor for tau_n/p
      data f2d3  /0.666666667/
      do i=1,200
        rho = 0.001*i
c
        rhop   = 0.5D0*rho
        rhon   = 0.5D0*rho
        taup   = betab*rhop**f2d3*rhop
        taun   = betab*rhon**f2d3*rhon
        esymm = (bhf_poten(rho*0.5D0,rho*0.5D0) 
     &         +  20.73553*(taup+taun))/rho
c
        tau = betab*rho**f2d3
        easymm = bhf_poten(rho,0.0D0)/rho  + 20.73553*tau
c
        write(6,'(1x,3(3pg14.6))') rho,esymm,easymm
        write(6,'(1x,3(3pg14.6))') 
     &    rho,bhf_toten(rho*0.5D0,rho*0.5D0),bhf_toten(rho,0.0D0)
      enddo
      rho0  = 0.155
      rho02 = rho0/2.0D0
      drho  = 1.0D-4
      drho2 = 2.0D0*drho
      compr = 9.0D0*rho0**2*
     &        (bhf_toten(rho02+drho,rho02+drho)
     &         -2.0D0*bhf_toten(rho02,rho02)
     &         +bhf_toten(rho02-drho,rho02-drho)             )
     &        /drho2**2
      write(6,'(a,1pg14.5)') ' compr.=',compr
      asymme = (bhf_toten(rho02+drho,rho02-drho)
     &          -bhf_toten(rho02,rho02))/drho2**2*rho0**2
      write(6,'(a,1pg14.5)') ' asymme=',asymme
      stop
      end

c-----bhf_toten----------------------------------------------------
c
      double precision function bhf_toten(rhop,rhon)
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   !
c
      data betab /5.742467999/     ! (3/5)*(3*pi^2)^*(2/3)  factor for tau_n/p
      data f2d3  /0.666666667/
c
c-----------------------------------------------------------------
c
      rho = rhop+rhon
      taup   = betab*rhop**f2d3*rhop
      taun   = betab*rhon**f2d3*rhon
      bhf_toten = (bhf_poten(rhop,rhon)+20.73553*(taup+taun))/rho
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




