c
c    in this version the " bigJ " term   have been suppressed and the EnergyDensity has been added
c
c               compile with      gfortran -ffixed-line-length-132
c
c    Use of the code:
c    There are 6 functions:
c                                               UPotProton(arg),UPotNeutron(arg),
c                                               WPotProton(arg),WPotNeutron(arg),
c                                               OneOver2MeffProton(arg),OneOver2MeffNeutron(arg)
c                                               EnergyDensity(arg)
c     which follow exactely the notation of Vautherin Brink
c
c    The arguments are :
c               arg=rhop,rhon,taup,taun,nabla2rhop,nabla2rhon,nablaDotBigJp,nablaDotBigJn,
c                       gradrhop,gradrhon,BigJp,BigJn
c
c       they are the NUMERICAL values of the quantities  {rhop,...} for  the considered point and direction
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

      WRITE(6,'(a)') 'ATTACH_QMC'
      WRITE(6,'(5(1pg13.5))') par(1:18)
      WRITE(6,'(20(i5))') jpact(1:18)
      WRITE(6,'(6(1pg13.5))') b0,b1,b2,b3,b0p,b1p,b2p,b3p,power
      WRITE(6,'(a,3(1pg13.5))') 'PAIR: ',pairfp,pairfn,rho0pr

      return
      end

      subroutine SetParameters()
!      implicit none
      IMPLICIT REAL(8) (A-H,O-Z)
      include 'skypar.inc'
      real*8 freeR
      real*8 h2mp,h2mn
      real*8 NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,Grho,dpol,
     &       isoscalarMu,isovectorMu
      real*8 mpi,fpi,ga,pionfactor,ActualPionFactor
      common/Parameters/NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,
     &                  Grho,dpol,isoscalarMu,isovectorMu
      common/Pion/mpi,fpi,ga,pionfactor,ActualPionFactor
      real*8 hc    !  hbar*c
      parameter (hc=197.3285879d0)
      integer, save :: ncalled=0
!
      NucleonMass=4.753458673288399D0
      h2mp = hc/(2D0*NucleonMass)
      h2mn = hc/(2D0*NucleonMass)
      amfc0(1) = hc/(2D0*NucleonMass)
      amfc0(2) = hc/(2D0*NucleonMass)

      Msigma=3.54735721887194
      Momega=3.952769472457305
!      Mrho=3.9020929407591343
      mpi=0.709471443774388      
      fpi=0.4712917447929863
      ga=1.26
c   **************************************************************
c**** here choose the case with or without pion (pionfactor=0 or 1)****
c   **************************************************************
      pionfactor=1D0                   ! setting at compile time
      if(pionfactor.eq.0D0) then
        Gsigma = 11.84669
        Gomega = 8.268321
        Grho = 7.681668
        Msigma = 3.660819
c      Gsigma=11.365282651688423
c      Gomega=7.959633529464553
c      Grho=7.556370766685059
      dpol=0.17969968018830365
      end if
      if(pionfactor.eq.1)then
        Gsigma = 11.84669
        Gomega = 8.268321
        Grho = 7.681668
        Msigma = 3.660819
c      Gsigma=10.69248392612409
c      Gomega=7.762652316091267
c      Grho=6.700101658782705
      dpol=0.17969968018830365
      end if
c
      ActualPionFactor=(-9*ga**2*mpi**2)/(32.*fpi**2)
      isovectorMu=4.7
      isoscalarMu=0.88
c
c     overwrite default settings with parameters from 'for001'
c
      ibinit=1
      IF(ibinit==1) THEN
      Gsigma = b0 
      Gomega = b1
      Grho   = b2
      isoscalarMu = b3
      Msigma = b0p
      Momega = b1p
      Mrho = b2p
      isovectorMu = b3p
      dpol = power
      END IF
c
      freeR=1
c
      IF(ncalled<3) THEN
        WRITE(6,*) ' b:',b0,b1,b2,b3,b0p,b1p,b2p,b3p
        print*,"pion factor=",pionfactor
        print*,"actual pion factor=",ActualPionFactor
        print*,"sigma mass=",Msigma
        print*,"omege mass=",Momega
        print*,"rho mass=",Mrho
        print*,"Gsigma=",Gsigma
        print*,"Gomega=",Gomega
        print*,"Grho=",Grho
        print*,"Isovect=",isovectorMu
        print*,"Isoscal=",isoscalarMu
        print*,"scalar polarisability dpol=",dpol
        print*,"free nucleon radius :",freeR
        ncalled = ncalled +1
      END IF

        return
        end
c****************************************************************************
        doubleprecision function kfderho(rho)
        implicit none
        real*8 rho,pi
        data pi/3.141592/
        kfderho=(3*pi*pi*rho)**(1D0/3D0)
        return
        end
c****************************************************************************
        doubleprecision function UPotProton(rhop,rhon,taup,taun,
     &             nabla2rhop,nabla2rhon,nablaDotBigJp,nablaDotBigJn,
     &             gradrhop,gradrhon,BigJp,BigJn)
        implicit none
      real*8 rhop,rhon,taup,taun,nabla2rhop,nabla2rhon,nablaDotBigJp
      real*8 nablaDotBigJn,gradrhop,gradrhon,BigJp,BigJn
      real*8 NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,Grho,dpol
      real*8 isoscalarMu,isovectorMu
      real*8 mpi,fpi,ga,pionfactor,ActualPionFactor
      real*8 kfderho,kfp,kfn,rhop23
      real*8 pi
      data pi/3.141592/   
      common/Parameters/NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,
     &                  Grho,dpol,isoscalarMu,isovectorMu
      common/Pion/mpi,fpi,ga,pionfactor,ActualPionFactor      
      logical, save :: tfirst=.true.
c
      rhop=MAX(rhop,1D-80)
      rhon=MAX(rhon,1D-80)
      kfp=kfderho(rhop)
      kfn=kfderho(rhon)
      UPotProton= nabla2rhon*(Gomega/Momega**2 - 
     &            (3*Grho)/(8.*Mrho**2) - 
     -     Gsigma/Msigma**2 - Gsigma/(2.*NucleonMass**2)) + 
     -  nabla2rhop*((3*Gomega)/(4.*Momega**2) +
     &    (3*Grho)/(16.*Mrho**2) - 
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
     &  (2.*(1 + dpol*Gsigma*(rhon + rhop))**2) - 
     -  (Gsigma*(rhon + rhop))/(1 + dpol*Gsigma*(rhon + rhop)) + 
     -  rhop*(Gomega/2. + Grho/8. + 
     -     Gsigma/(2.*(1 + dpol*Gsigma*(rhon + rhop))**3)) + 
     -  (Grho/(4.*Mrho**2) + Gsigma/(2.*NucleonMass**2))*taun + 
     -  (Gomega/(2.*Momega**2) 
     &      + Grho/(8.*Mrho**2) - Gsigma/(2.*Msigma**2) +
     -     Gsigma/(4.*NucleonMass**2))*taup
        IF(pionfactor.NE.0D0) UPotProton=UPotProton+
     &   pionfactor*(ActualPionFactor*kFp*(8*kFn*kFp + 4*kFp**2 - 
     -      4*kFp*mpi*atan((2*kFp)/mpi) + 
     &      8*kFp*mpi*atan((-kFn + kFp)/mpi) - 
     -      8*kFp*mpi*atan((kFn + kFp)/mpi) + 
     -      mpi**2*Log(1 + (4*kFp**2)/mpi**2) + 
     -      2*kFn**2*Log(((kFn + kFp)**2 + mpi**2)/
     -         ((-kFn + kFp)**2 + mpi**2)) - 
     -      2*kFp**2*Log(((kFn + kFp)**2 + mpi**2)/
     -         ((-kFn + kFp)**2 + mpi**2)) + 
     -      2*mpi**2*
     &      Log(((kFn + kFp)**2 + mpi**2)/((-kFn + kFp)**2 + mpi**2))
     -      ))/(18.*3D0**0.6666666666666666D0*Pi**3.3333333333333335D0*
     -    rhop**0.6666666666666666D0)
        IF(tfirst) THEN
          WRITE(*,*) 'Prot:',rhop,rhon,UPotProton
          WRITE(*,*) 'kf:',kfp,kfn
          WRITE(*,*) nabla2rhon,nabla2rhop,nablaDotBigJn,nablaDotBigJp
          WRITE(*,*) taup,taun
          tfirst=.false.
        END IF
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
      real*8 mpi,fpi,ga,pionfactor,ActualPionFactor
      real*8 kfderho,kfp,kfn
      real*8 pi
      data pi/3.141592/         
      common/Parameters/NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,
     &                  Grho,dpol,isoscalarMu,isovectorMu
      common/Pion/mpi,fpi,ga,pionfactor,ActualPionFactor
      logical, save :: tfirst=.true.
c      
      rhop=MAX(rhop,1D-80)
      rhon=MAX(rhon,1D-80)
      kfn=kfderho(rhon)
      kfp=kfderho(rhop)
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
     -  (dpol*Gsigma**2*(rhon + rhop)**2)/
     &  (2.*(1 + dpol*Gsigma*(rhon + rhop))**2) - 
     -  (Gsigma*(rhon + rhop))/(1 + dpol*Gsigma*(rhon + rhop)) + 
     -  rhon*(Gomega/2. + Grho/8. + 
     -     Gsigma/(2.*(1 + dpol*Gsigma*(rhon + rhop))**3)) + 
     -  (Gomega/(2.*Momega**2) + Grho/(8.*Mrho**2) -
     &   Gsigma/(2.*Msigma**2) + 
     -     Gsigma/(4.*NucleonMass**2))*taun + 
     -  (Grho/(4.*Mrho**2) + 
     &  Gsigma/(2.*NucleonMass**2))*taup
        IF(pionfactor.NE.0D0) UPotNeutron=UPotNeutron+
     &  pionfactor*(ActualPionFactor*kFn*(4*kFn**2 + 8*kFn*kFp - 
     -   4*kFn*mpi*atan((2*kFn)/mpi) - 8*kFn*mpi*atan((-kFn+kFp)/mpi) - 
     -   8*kFn*mpi*atan((kFn + kFp)/mpi) + 
     -   mpi**2*Log(1 + (4*kFn**2)/mpi**2) - 
     -   2*kFn**2*Log(((kFn + kFp)**2 + mpi**2)/
     -      ((-kFn + kFp)**2 + mpi**2)) + 
     -   2*kFp**2*Log(((kFn + kFp)**2 + mpi**2)/
     -      ((-kFn + kFp)**2 + mpi**2)) + 
     -   2*mpi**2*
     &   Log(((kFn + kFp)**2 + mpi**2)/((-kFn + kFp)**2 + mpi**2))
     -      ))/(18.*3**0.6666666666666666D0*Pi**3.3333333333333335D0*
     -    rhon**0.6666666666666666D0)
        IF(tfirst) THEN
          WRITE(*,*) 'Neut.:',rhop,rhon,UPotNeutron
          tfirst=.false.
        END IF
        return
        end
c****************************************************************************
      doubleprecision function WPotProton(rhop,rhon,taup,taun,
     &             nabla2rhop,nabla2rhon,nablaDotBigJp,nablaDotBigJn,
     1             gradrhop,gradrhon,BigJp,BigJn)
      implicit none
      real*8 rhop,rhon,taup,taun,nabla2rhop,nabla2rhon,nablaDotBigJp,
     1 nablaDotBigJn,gradrhop,gradrhon,BigJp,BigJn
      real*8 NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,Grho,dpol
      real*8 isoscalarMu,isovectorMu
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
      real*8 NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,Grho,dpol
      real*8 isoscalarMu,isovectorMu
      common/Parameters/NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,
     &                  Grho,dpol,isoscalarMu,isovectorMu
      WPotNeutron= gradrhop*(-Gomega/(4.*NucleonMass**2) +
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
!      OneOver2MeffProton=    1/(2.*NucleonMass) +
      OneOver2MeffProton=
     &   (Grho/(4.*Mrho**2) + Gsigma/(2.*NucleonMass**2))*
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
!      OneOver2MeffNeutron= 1/(2.*NucleonMass) +
      OneOver2MeffNeutron=
     & (Gomega/(2.*Momega**2) + Grho/(8.*Mrho**2) - 
     -     Gsigma/(2.*Msigma**2) + Gsigma/(4.*NucleonMass**2))*rhon + 
     -  (Grho/(4.*Mrho**2) + Gsigma/(2.*NucleonMass**2))*rhop
        return
        end
c****************************************************************************
c  we now call TheEnergyDensity["QMC","no bigJ"]   because of the change in mathematica
c  TheEnergyDensity["QMC"] is now the full energy density
c
      doubleprecision function EnergyDensity(rhop,rhon,taup,taun,
     &          nabla2rhop,nabla2rhon,nablaDotBigJp,nablaDotBigJn,
     1          gradrhop,gradrhon,BigJp,BigJn)
      implicit none
      real*8 rhop,rhon,taup,taun,nabla2rhop,nabla2rhon,nablaDotBigJp,
     1 nablaDotBigJn,gradrhop,gradrhon,BigJp,BigJn
      real*8 NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,Grho,dpol,
     & isoscalarMu,isovectorMu
      real*8 mpi,fpi,ga,pionfactor,ActualPionFactor
      real*8 kfderho,kfp,kfn     
      real*8 pi
      data pi/3.141592/     
      common/Parameters/NucleonMass,Msigma,Momega,Mrho,Gsigma,Gomega,
     &                  Grho,dpol,isoscalarMu,isovectorMu
      common/Pion/mpi,fpi,ga,pionfactor,ActualPionFactor     
      logical, save :: tfirst=.true.
      rhop=MAX(rhop,1D-80)
      rhon=MAX(rhon,1D-80)
      kfn=kfderho(rhon)
      kfp=kfderho(rhop)
        EnergyDensity=        (nabla2rhon*((12*Gomega)/Momega**2 +
     & (3*Grho)/Mrho**2 - 
     -   (12*Gsigma)/Msigma**2 - (4*Gsigma)/NucleonMass**2)*rhon)/32. + 
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
     -  (( (4*Gomega*rhon)/Momega**2 + (Grho*rhon)/Mrho**2 - 
     -       (4*Gsigma*rhon)/Msigma**2 + (2*Grho*rhop)/Mrho**2 + 
     -       (2*Gsigma*(rhon + 2*rhop))/NucleonMass**2)*taun)/8. + 
     -  (( (4*Gomega*rhop)/Momega**2 - 
     -      (4*Gsigma*rhop)/Msigma**2 + (Grho*(2*rhon + rhop))/Mrho**2 + 
     -       (2*Gsigma*(2*rhon + rhop))/NucleonMass**2)*taup)/8.
        IF(pionfactor.NE.0D0)
     &   EnergyDensity=EnergyDensity+pionfactor*(ActualPionFactor*
     &       ((rhon**2*
     -        (24*kFn**4 - 4*kFn**2*mpi**2 - 
     -         32*kFn**3*mpi*atan((2*kFn)/mpi) + 
     -         (12*kFn**2*mpi**2+mpi**4)*Log(1D0+(4*kFn**2)/mpi**2)))/
     -       kFn**6 + (rhop**2*(24*kFp**4 - 4*kFp**2*mpi**2 - 
     -         32*kFp**3*mpi*atan((2*kFp)/mpi) + 
     -         (12*kFp**2*mpi**2 + mpi**4)*Log(1 + (4*kFp**2)/mpi**2)))/
     -       kFp**6 + (4*rhon*rhop*
     -         (4*kFn*kFp*(3*kFn**2 + 3*kFp**2 - mpi**2) + 
     -           16*(kFn**3 - kFp**3)*mpi*atan((kFn - kFp)/mpi) - 
     -           16*(kFn**3 + kFp**3)*mpi*atan((kFn + kFp)/mpi) + 
     -           (-3*kFn**4 - 3*kFp**4 + 6*kFp**2*mpi**2 + mpi**4 + 
     -              6*kFn**2*(kFp**2 + mpi**2))*
     -            Log((kFn**2 + 2*kFn*kFp + kFp**2 + mpi**2)/
     -              (kFn**2 - 2*kFn*kFp + kFp**2 + mpi**2))))/
     &            (kFn**3*kFp**3))
     -    )/48.
        IF(tfirst) THEN
          WRITE(*,*) 'Energy:',rhop,rhon,EnergyDensity
          tfirst=.false.
        END IF
        return
        end


        
