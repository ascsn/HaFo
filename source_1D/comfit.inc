C INCLUDE BLOCK: COMFIT
c
      parameter (kparmx=38)    ! max. number fit parameters
c      parameter (kparmx=17)    ! max. number fit parameters
      parameter (kdat=550)      ! max. number of fit data
      parameter (kderij=12*kdat)
      INTEGER, PARAMETER :: kdifmax=20,k3gapmax=90
!      parameter (parmnatfac=0.5D0)  ! step factor for natural params
c
c     table if fit nuclei:
c       IAY    =  mass of nucleus
c       IPROT  = charge of nucleus
c       ZPROT  = IPROT but as real number (to cope with neutron matter)
c       BETA   = beta20 from B(E2)
c       BETACT = temporary values
c       ERRDIF = error on difference of observable, see NUCDIF
c       NUCDIF = pointer on isotopic difference
c                NUCDIF(1,IDIF) = 1. nucleus as listed in IAY,IPROT
c                NUCDIF(2,IDIF) = 2. nucleus as listed in IAY,IPROT
c                NUCDIF(3,IDIF) = Nr. of observable
c       MAXDIF = max. nr. of isotopic differences
c       parmnat = matrix of natural parameter vectors
c       n3gap  = book-keeping for 3gap
c                  entry 1-3 = the 3 nuclei forming the 3gap
c                  entry 4   = proton gap or neutron gap
c                  entry 5   = type of observable (typically 1 or 4)
c 
      logical tsimple
      character*60  commen
      COMMON   /HAFO/ ERRDIF(kdifmax),err3gap(k3gapmax),zprot(kdat),
     &                rsymax,rntmax,nntmax,nsymax,
     &                IAY(kdat),IPROT(kdat),NUCDIF(3,kdifmax),MAXDIF,
     &                n3gap(5,k3gapmax),max3gap,
     &                commen
C
      COMMON /FITPAR/ ZARRAY(kparmx,kparmx),DPARMO(kparmx),CHISU0(12),
     &                curvature(kparmx,kparmx),
     &                parmnat(kparmx,kparmx),
     *                yexp(12,kdat),sigmay(12,kdat),yfit(12,kdat),
     *                sumchi(14),chi(14,kdat),chi0sav(14,kdat),              ! 14 = 12+2
     &                varextrap(12,kdat),                  ! estimated error
     *                CHISQR,CHIA,chimin,
     *                SOPROT,SONEUT,WIDPR,WIDNE,FLUC,freep,endcon,
     &                parmnatfac,rstep_explicit,
     *                IPARVA(kparmx),JP(kparmx),LIFIT(12),IFIT1(12),
     *                NPTS,IWRTE,IST1,IFITR,NPTSG,NDATS,NPARIN,
     *                IFIT,IDAT1,IFALT2,ISTEUR, IFTINP,IB2FIX,nparfit,
     &                ngrid_explicit,
     &                tsimple
C
      common /fitwrk/ chimom(kparmx),fderij(kderij,kparmx),
     &                fjacobian(kderij,0:kparmx)
!
      LOGICAL, PARAMETER :: tdiffobs=.TRUE.
!      REAL(8), ALLOCATABLE :: fjacobian(:,:)               ! Jacobian incl. difference obs.
C END INCLUDE BLOCK
