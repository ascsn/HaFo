c
c    in this version the " bigJ " term   have been suppressed 
c    and the EnergyDensity has been added
c
c    moreover: the trivial part of kinetic energy is dropped here
c
c    Use of the code:
c    There are 6 functions:
c                          UPotProton(arg),UPotNeutron(arg),
c                          WPotProton(arg),WPotNeutron(arg),
c                          OneOver2MeffProton(arg),
c                          OneOver2MeffNeutron(arg)
c                          EnergyDensity(arg)
c     which follow exactely the notation of Vautherin Brink
c
c    The arguments are :
c         arg=rhop,rhon,taup,taun,nabla2rhop,nabla2rhon,
c             nablaDotBigJp,nablaDotBigJn,gradrhop,gradrhon,BigJp,BigJn
c
c       they are the NUMERICAL values of the quantities  {rhop,...} for  
c       the considered point and direction
c
c       nabla2rhop means laplacian acting on rhop
c       nablaDotBigJp means nabla.Jp
c       gradrhop means gradient of rhop
c       BigJp is the Jp vector density defined bu WB
c
c       The working unit is the fermi**power for any quantity
c
c       The parameters are set by SetParameters which must be called once
c
!      call SetParameters
c***********************************************************************
c  Here are the relevant routines. 
c***********************************************************************
      subroutine attach_QMC(par,jpact)
c      implicit none
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   !
      include 'skypar.inc'
      include 'comfit.inc'
      real*8 par(kparmx)
      integer jpact(kparmx)

      b0     = par(jpact(1))        !  Gsigma
      b1     = par(jpact(2))        !  Gomega
      b2     = par(jpact(3))        !  Grho
      b3     = par(jpact(4))        !  isoscalarMu
      b0p    = par(jpact(6))        !  Msigma
      b1p    = par(jpact(7))        !  Momega
      b2p    = par(jpact(8))        !  Mrho
      b3p    = par(jpact(9))        !  isovectorMu
      power  = par(jpact(10))       !  dpol
      pairfp = par(jpact(16))
      pairfn = par(jpact(17))
      rho0pr = par(jpact(18))       ! default 0.16

!      WRITE(6,'(a)') 'ATTACH_QMC'
!      WRITE(6,'(5(1pg13.5))') par(1:18)
!      WRITE(6,'(20(i5))') jpact(1:18)
!      WRITE(6,'(6(1pg13.5))') b0,b1,b2,b3,b0p,b1p,b2p,b3p,power
!      WRITE(6,'(a,3(1pg13.5))') 'PAIR: ',pairfp,pairfn,rho0pr

      return
      end

c****************************************************************************
      subroutine SetParameters
!      implicit none
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   !
      include 'skypar.inc'
      real*8 NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,Grho,dpol
      real*8 isoscalarMu,isovectorMu
      common/Parameters/NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,
     &                  Grho,dpol,isoscalarMu,isovectorMu
      real*8 hc    !  hbar*c
      parameter (hc=197.3285879)
!      logical tfirst
!      data tfirst/.true./
c
!      WRITE(6,'(6(1pg13.5))') b0,b1,b2,b3,b0p,b1p,b2p,b3p,power
      Gsigma = b0 
      Gomega = b1
      Grho   = b2
      isoscalarMu = b3
      Msigma = b0p
      Momega = b1p
      Mrho = b2p
      isovectorMu = b3p
      dpol = power

c      Msigma=3.547375
c      Momega=3.952769472457305
c      Mrho=3.9020929407591343
c      Gsigma=11.517164500018843
c      Gomega=8.327580982701214
c      Grho=8.348185363640384
c      isovectorMu=4.7
c      isoscalarMu=0.88
c      dpol=0.17969968018830365
c

c
      NucleonMass=4.753458673288399D0
      amfc0(1) = hc/(2D0*NucleonMass)
      amfc0(2) = hc/(2D0*NucleonMass)

c
c
c      if(tfirst) then
!      WRITE(6,'(a,6(1pg13.5)/4x,6(1pg13.5))')
!     &   'QMC:',Gsigma,Gomega,Grho,Msigma,Momega,Mrho,
!     &          dpol,isoscalarMu,isovectorMu,amfc0
c        print*,"sigma mass=",Msigma
c        print*,"Gsigma=",Gsigma
c        print*,"Gomega=",Gomega
c        print*,"Grho=",Grho
c        print*,"scalar polarisability dpol=",dpol
c        tfirst=.false.
c      endif

      return
      end

c****************************************************************************
      doubleprecision function UPotProton(rhop,rhon,taup,taun,
     &           nabla2rhop,nabla2rhon,nablaDotBigJp,nablaDotBigJn,
     &           gradrhop,gradrhon,BigJp,BigJn)
      implicit none
      real*8 rhop,rhon,taup,taun,nabla2rhop,nabla2rhon,nablaDotBigJp
      real*8 nablaDotBigJn,gradrhop,gradrhon,BigJp,BigJn
      real*8 NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,Grho,dpol
      real*8 isoscalarMu,isovectorMu
      common/Parameters/NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,
     &                  Grho,dpol,isoscalarMu,isovectorMu
c
      UPotProton= nabla2rhon*(Gomega/Momega**2 - (3*Grho)/(8.*Mrho**2) - 
     -     Gsigma/Msigma**2 - Gsigma/(2.*NucleonMass**2)) +
     -  nabla2rhop*((3*Gomega)/(4.*Momega**2) + (3*Grho)/(16.*Mrho**2) - 
     -     (3*Gsigma)/(4.*Msigma**2) - Gsigma/(4.*NucleonMass**2)) +
     -  nablaDotBigJp*((-3*(Gsigma + Gomega*(-1 + 2*isoscalarMu)))/
     -      (8.*NucleonMass**2) -
     -     (3*Grho*(-1 + 2*isovectorMu))/(32.*NucleonMass**2)) -
     -  ((Gsigma + Gomega*(-1 + 2*isoscalarMu))*nablaDotBigJn)/
     -   (4.*NucleonMass**2) + (Gomega - Grho/2.)*rhon -
     -  (3*dpol*Gsigma**2*rhon**2)/
     &  (4.*(1 + dpol*Gsigma*(rhon + rhop))**4) -
     -  (3*dpol*Gsigma**2*rhop**2)/
     &  (4.*(1 + dpol*Gsigma*(rhon + rhop))**4) +
     -  (dpol*Gsigma**2*(rhon + rhop)**2)/
     &            (2.*(1 + dpol*Gsigma*(rhon + rhop))**2) -
     -  (Gsigma*(rhon + rhop))/(1 + dpol*Gsigma*(rhon + rhop)) +
     -  rhop*(Gomega/2. + Grho/8. +
     -     Gsigma/(2.*(1 + dpol*Gsigma*(rhon + rhop))**3)) +
     -  (Grho/(4.*Mrho**2) + Gsigma/(2.*NucleonMass**2))*taun +
     -  (Gomega/(2.*Momega**2)
     &      + Grho/(8.*Mrho**2) - Gsigma/(2.*Msigma**2) +
     -     Gsigma/(4.*NucleonMass**2))*taup
      return
      end
c****************************************************************************
      doubleprecision function UPotNeutron(rhop,rhon,taup,taun,
     &             nabla2rhop,nabla2rhon,nablaDotBigJp,nablaDotBigJn,
     &             gradrhop,gradrhon,BigJp,BigJn)
      implicit none
      real*8 rhop,rhon,taup,taun,nabla2rhop,nabla2rhon,nablaDotBigJp,
     1 nablaDotBigJn,gradrhop,gradrhon,BigJp,BigJn
      real*8 NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,Grho,dpol
      real*8 isoscalarMu,isovectorMu
      common/Parameters/NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,
     &                  Grho,dpol,isoscalarMu,isovectorMu
c
      UPotNeutron= nabla2rhop*(Gomega/Momega**2 - (3*Grho)/(8.*Mrho**2)- 
     -     Gsigma/Msigma**2 - Gsigma/(2.*NucleonMass**2)) +
     -  nabla2rhon*((3*Gomega)/(4.*Momega**2) + (3*Grho)/(16.*Mrho**2) - 
     -     (3*Gsigma)/(4.*Msigma**2) - Gsigma/(4.*NucleonMass**2)) +
     -  nablaDotBigJn*((-3*(Gsigma + Gomega*(-1 + 2*isoscalarMu)))/
     -      (8.*NucleonMass**2) -
     -     (3*Grho*(-1 + 2*isovectorMu))/(32.*NucleonMass**2)) -
     -  ((Gsigma + Gomega*(-1 + 2*isoscalarMu))*nablaDotBigJp)/
     -   (4.*NucleonMass**2) + (Gomega - Grho/2.)*rhop -
     -  (3*dpol*Gsigma**2*rhon**2)/
     &  (4.*(1 + dpol*Gsigma*(rhon + rhop))**4) -
     -  (3*dpol*Gsigma**2*rhop**2)/
     &  (4.*(1 + dpol*Gsigma*(rhon + rhop))**4) +
     -  (dpol*Gsigma**2*(rhon+rhop)**2)/
     &  (2.*(1+dpol*Gsigma*(rhon+rhop))**2) -
     -  (Gsigma*(rhon + rhop))/(1 + dpol*Gsigma*(rhon + rhop)) +
     -  rhon*(Gomega/2. + Grho/8. +
     -     Gsigma/(2.*(1 + dpol*Gsigma*(rhon + rhop))**3)) +
     -  (Gomega/(2.*Momega**2)+Grho/(8.*Mrho**2)-Gsigma/(2.*Msigma**2) + 
     -     Gsigma/(4.*NucleonMass**2))*taun +
     -  (Grho/(4.*Mrho**2) + Gsigma/(2.*NucleonMass**2))*taup
      return
      end
c****************************************************************************
      doubleprecision function WPotProton(rhop,rhon,taup,taun,
     &             nabla2rhop,nabla2rhon,nablaDotBigJp,nablaDotBigJn,
     1             gradrhop,gradrhon,BigJp,BigJn)
      implicit none
      real*8 rhop,rhon,taup,taun,nabla2rhop,nabla2rhon,nablaDotBigJp,
     1 nablaDotBigJn,gradrhop,gradrhon,BigJp,BigJn
      real*8 NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,Grho,dpol,
     &      isoscalarMu,isovectorMu
      common/Parameters/NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,
     &         Grho,dpol,isoscalarMu,isovectorMu
      WPotProton=  gradrhon*(-Gomega/(4.*NucleonMass**2) + 
     &     Gsigma/(4.*NucleonMass**2) + 
     -     (Gomega*isoscalarMu)/(2.*NucleonMass**2)) + 
     -  gradrhop*((-3*Gomega)/(8.*NucleonMass**2) - 
     -  (3*Grho)/(32.*NucleonMass**2) + (3*Gsigma)/(8.*NucleonMass**2) + 
     -   (3*Gomega*isoscalarMu)/(4.*NucleonMass**2) + 
     -   (3*Grho*isovectorMu)/(16.*NucleonMass**2))
      return
      end
c****************************************************************************
      doubleprecision function WPotNeutron(rhop,rhon,taup,taun,
     &          nabla2rhop,nabla2rhon,nablaDotBigJp,nablaDotBigJn,
     1          gradrhop,gradrhon,BigJp,BigJn)
      implicit none
      real*8 rhop,rhon,taup,taun,nabla2rhop,nabla2rhon,nablaDotBigJp,
     1 nablaDotBigJn,gradrhop,gradrhon,BigJp,BigJn
      real*8 NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,Grho,dpol,
     &  isoscalarMu,isovectorMu
      common/Parameters/NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,
     &                  Grho,dpol,isoscalarMu,isovectorMu
      WPotNeutron=  gradrhop*(-Gomega/(4.*NucleonMass**2) + 
     &     Gsigma/(4.*NucleonMass**2) + 
     -     (Gomega*isoscalarMu)/(2.*NucleonMass**2)) + 
     -  gradrhon*((-3*Gomega)/(8.*NucleonMass**2) - 
     -  (3*Grho)/(32.*NucleonMass**2) + (3*Gsigma)/(8.*NucleonMass**2) + 
     -   (3*Gomega*isoscalarMu)/(4.*NucleonMass**2) + 
     -   (3*Grho*isovectorMu)/(16.*NucleonMass**2))
      return
      end
c****************************************************************************
      doubleprecision function OneOver2MeffProton(rhop,rhon,taup,taun,
     &            nabla2rhop,nabla2rhon,nablaDotBigJp,nablaDotBigJn,
     1            gradrhop,gradrhon,BigJp,BigJn)
      implicit none
      real*8 rhop,rhon,taup,taun,nabla2rhop,nabla2rhon,nablaDotBigJp,
     1 nablaDotBigJn,gradrhop,gradrhon,BigJp,BigJn
      real*8 NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,Grho,dpol,
     &   isoscalarMu,isovectorMu
      common/Parameters/NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,
     &                  Grho,dpol,isoscalarMu,isovectorMu
c      OneOver2MeffProton=   1/(2.*NucleonMass) + 
      OneOver2MeffProton=  
     &   (Grho/(4.*Mrho**2) + 
     &   Gsigma/(2.*NucleonMass**2))*
     -   rhon + (Gomega/(2.*Momega**2) + Grho/(8.*Mrho**2) - 
     -     Gsigma/(2.*Msigma**2) + Gsigma/(4.*NucleonMass**2))*rhop
      return
      end
c****************************************************************************
      doubleprecision function OneOver2MeffNeutron(rhop,rhon,taup,taun,
     &       nabla2rhop,nabla2rhon,nablaDotBigJp,nablaDotBigJn,
     1       gradrhop,gradrhon,BigJp,BigJn)
      implicit none
      real*8 rhop,rhon,taup,taun,nabla2rhop,nabla2rhon,nablaDotBigJp,
     1 nablaDotBigJn,gradrhop,gradrhon,BigJp,BigJn
      real*8 NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,Grho,dpol,
     & isoscalarMu,isovectorMu
      common/Parameters/NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,
     &                  Grho,dpol,isoscalarMu,isovectorMu
c      OneOver2MeffNeutron= 1/(2.*NucleonMass) +
      OneOver2MeffNeutron= 
     & (Gomega/(2.*Momega**2) + 
     &   Grho/(8.*Mrho**2) - 
     -     Gsigma/(2.*Msigma**2) + Gsigma/(4.*NucleonMass**2))*rhon + 
     -  (Grho/(4.*Mrho**2) + Gsigma/(2.*NucleonMass**2))*rhop
      return
      end
c****************************************************************************
      doubleprecision function EnergyDensity(rhop,rhon,taup,taun,
     &          nabla2rhop,nabla2rhon,nablaDotBigJp,nablaDotBigJn,
     1          gradrhop,gradrhon,BigJp,BigJn)
      implicit none
      real*8 rhop,rhon,taup,taun,nabla2rhop,nabla2rhon,nablaDotBigJp,
     1 nablaDotBigJn,gradrhop,gradrhon,BigJp,BigJn
      real*8 NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,Grho,dpol,
     & isoscalarMu,isovectorMu
      common/Parameters/NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,
     &                  Grho,dpol,isoscalarMu,isovectorMu
        EnergyDensity=    (nabla2rhon*((12*Gomega)/Momega**2 + 
     &  (3*Grho)/Mrho**2 - 
     -    (12*Gsigma)/Msigma**2 - (4*Gsigma)/NucleonMass**2)*rhon)/32. + 
     -  (nabla2rhop*((8*Gomega)/Momega**2 - (3*Grho)/Mrho**2 - 
     -     (8*Gsigma)/Msigma**2 - (4*Gsigma)/NucleonMass**2)*rhon)/16. - 
     -  (3*(4*Gsigma + Gomega*(-4 + 8*isoscalarMu) + 
     -       Grho*(-1 + 2*isovectorMu))*nablaDotBigJn*rhon)/
     -   (32.*NucleonMass**2) + 
     -  ((Gomega - Gsigma - 2*Gomega*isoscalarMu)*nablaDotBigJp*rhon)/
     -   (4.*NucleonMass**2) + (nabla2rhop*
     -     ((12*Gomega)/Momega**2 + (3*Grho)/Mrho**2 - 
     -    (12*Gsigma)/Msigma**2 - (4*Gsigma)/NucleonMass**2)*rhop)/32. + 
     -  (nabla2rhon*((8*Gomega)/Momega**2 - (3*Grho)/Mrho**2 - 
     -     (8*Gsigma)/Msigma**2 - (4*Gsigma)/NucleonMass**2)*rhop)/16. + 
     -  ((Gomega - Gsigma - 2*Gomega*isoscalarMu)*nablaDotBigJn*rhop)/
     -   (4.*NucleonMass**2) - (3*
     -     (4*Gsigma + Gomega*(-4 + 8*isoscalarMu) + 
     -       Grho*(-1 + 2*isovectorMu))*nablaDotBigJp*rhop)/
     -   (32.*NucleonMass**2) + 
     -  (Grho*(rhon**2 - 8*rhon*rhop + rhop**2) + 
     -     4*Gomega*(rhon**2 + 4*rhon*rhop + rhop**2) + 
     -     (4*Gsigma*(rhon**2 + rhop**2 - 
     -          2*rhon**2*(1 + dpol*Gsigma*(rhon + rhop))**2 - 
     -          4*rhon*rhop*(1 + dpol*Gsigma*(rhon + rhop))**2 - 
     -          2*rhop**2*(1 + dpol*Gsigma*(rhon + rhop))**2))/
     -      (1 + dpol*Gsigma*(rhon + rhop))**3)/16. + 
     -  ((            (4*Gomega*rhon)/Momega**2+(Grho*rhon)/Mrho**2 - 
     -       (4*Gsigma*rhon)/Msigma**2 + (2*Grho*rhop)/Mrho**2 + 
     -       (2*Gsigma*(rhon + 2*rhop))/NucleonMass**2)*taun)/8. + 
     -  ((            (4*Gomega*rhop)/Momega**2 - 
     -     (4*Gsigma*rhop)/Msigma**2 + (Grho*(2*rhon + rhop))/Mrho**2 + 
     -       (2*Gsigma*(2*rhon + rhop))/NucleonMass**2)*taup)/8.
!     -  ((4/NucleonMass +(4*Gomega*rhon)/Momega**2+(Grho*rhon)/Mrho**2 - 
!     -       (4*Gsigma*rhon)/Msigma**2 + (2*Grho*rhop)/Mrho**2 + 
!     -       (2*Gsigma*(rhon + 2*rhop))/NucleonMass**2)*taun)/8. + 
!     -  ((4/NucleonMass + (4*Gomega*rhop)/Momega**2 - 
!     -     (4*Gsigma*rhop)/Msigma**2 + (Grho*(2*rhon + rhop))/Mrho**2 + 
!     -       (2*Gsigma*(2*rhon + rhop))/NucleonMass**2)*taup)/8.
      return
      end

