!-----sample_reasonable------------------------------------------
!
      subroutine sample_reasonable(par,dpar,nparmx,nactive)
!
!     Takes optimal parametrization with corresponding covariance
!     matrix and computes a Gaussian distribution of "reasonable"
!     parametrizations by Monte-Carlo sampling.
!
!     par   = optimal parameters
!     dpar  = uncertainty on optimal parameters
!     --> curvature is communicated via common 'comfit'
!
      include 'params.inc'
      include 'comfit.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'extobs.inc'
      include 'constn.inc'
!
      REAL(8) :: par(kparmx)             ! optimum parametrization
      REAL(8) :: dpar(kparmx)            ! uncertainty on parameters
      INTEGER :: nparmx                  ! max. number of used parameters
      INTEGER :: nactive                 ! nr. of varied parameters
!
      REAL(8) :: a(kparmx),b(kparmx)     ! trial parametrizations
      REAL(8) :: binprob(1:21)           ! bins for sampling
!      REAL(8) :: bintest(-10:10)           ! bins for sampling
      INTEGER, PARAMETER :: nsamples=4000
      LOGICAL, PARAMETER :: tsample=.TRUE.
      REAL :: ran2
!
!--------------------------------------------------------------------
!
      IF(.NOT.tsample) RETURN
!
      binprob=0D0
!      CALL RANDOM_SEED
!      bintest=0D0
!      OPEN(702,file='ensemble.res')
      OPEN(702,file='ensemble2.res')
      chi0 = ffverg(par)
      chi2=0D0
      WRITE(702,'(a,1pg13.5)') '# reasonable ensemble, chi**2=',chi0
      WRITE(702,'(a)') '# isample chitry chiexact nucmat rms...'
      iseed = -1 !13739                      ! seed for random generator
      isample=0                          ! counter for sample
      DO itry=1,199999999
!        IF(MOD(itra,1000)==0) iseed=13739+itry/1000
        a = 0D0
        DO ipar=1,nactive
!          CALL RANDOM_NUMBER(ran2)
!          ranumb = 2D0*ran2-1D0
          ranumb = 2D0*ran1(iseed)-1D0
          a(1:nparmx)=1.1D0*ranumb*parmnat(1:nparmx,ipar)+a(1:nparmx)
!          ibin = NINT(ranumb/0.1D0)
!          ibin=MIN(ibin,10)
!          ibin=MAX(ibin,-10)
!          bintest(ibin)=1+bintest(ibin)
        END DO
        acc = 0D0
        DO ipar=1,nparmx
          acc = SUM(a(1:nparmx)*curvature(1:nparmx,ipar))*a(ipar) + acc
        END DO
        probab = EXP(-acc)
!        WRITE(703,'(i9,2(1pg13.5))') itry,acc,probab
!        IF(probab>0.5D0 .AND. ran1(iseed)< probab) THEN
        IF(probab>0.02D0) THEN
          b(1:nparmx) = par(1:nparmx)+a(1:nparmx)          
          call attach(b,jp)
          call initfc
          CALL nucmat(rhonm0,eovra0,compr0,efmas0,f0,f0p,g0,g0p,
     &                symen0,dasym0,d2asy0,skapp0,0,-2)
!          IF((compr0<270D0 .AND. rhonm0<0.153D0) .OR.
!     &       (compr0<270D0 .AND. rhonm0>0.1575D0) .OR.
!     &       (compr0>275D0 .AND. rhonm0<0.153D0) .OR.
!     &       (compr0>275D0 .AND. rhonm0>0.1575D0) ) THEN
          chi2=ffverg(b)
          WRITE(2,'(a,2i7,2(1pg13.5))') 
     &        ' try: ',itry,isample+1,acc,chi2-chi0
!          IF(chi2-chi0<2D0*acc) THEN
          IF(chi2-chi0<4D0) THEN
            isample = 1+isample
            ibin = INT(acc/0.1D0) + 1
            ibin = MIN(ibin,21)
            binprob(ibin) = 1+binprob(ibin)
!            call attach(b,jp)
!            call initfc
!            CALL nucmat(rhonm0,eovra0,compr0,efmas0,f0,f0p,g0,g0p,
!     &                symen0,dasym0,d2asy0,skapp0,0,-2)
            WRITE(2,'(2i7,5(1pg13.5))') isample,itry,
     &                  compr0,efmas0,symen0,dasym0,skapp0
            WRITE(2,'(20(1pg13.5))') b(1:nparmx)
            nprot  = 20
            nneut  = 28
            nnumb  = nprot+nneut
            anumb  = nnumb
            nmaxp  = 6
            nmaxn  = 7
            call hafosu
            rmsc48=rmsc
            rmsp48=rmsp
            rmsn48=rmsn
            nprot  = 114
            nneut  = 184
            nnumb  = nprot+nneut
            anumb  = nnumb
            nmaxp  = 0
            nmaxn  = 0
            call hafosu
            rmscSHE=rmsc
            rmspSHE=rmsp
            rmsnSHE=rmsn
            nprot  = 82
            nneut  = 126
            nnumb  = nprot+nneut
            anumb  = nnumb
            nmaxp  = 16
            nmaxn  = 22
            call hafosu
            rmsc208=rmsc
            rmsp208=rmsp
            rmsn208=rmsn
            if(abs(b4p).gt.1.0E-7) T4 = (t4+b4p)*0.666666667
            call ASYMT(T0,T1,T2,T3,t4,X0,X1,X2,X3,power,
     &               (amfc0(1)+amfc0(2))/2D0,
     +               YRHO0,eovera,xKOMPr,a4,YEFFM,YXKF,YELLO,
     +               surfe,surfas,YH,YP,YG,YQT1,YQT2)
            WRITE(702,'(i5,i9,42(1pg13.5))')
     &           isample,itry,acc,chi2-chi0,rhonm0,eovra0,compr0,efmas0,
     &           symen0,dasym0,skapp0,
     &           surfe,surfas,b4,b4p,power,pairfp,pairfn,rho0pr,
     &           rmsc48,rmsp48,rmsn48,rms,rmsc208,rmsp208,rmsn208,
     &           rmscSHE,rmspSHE,rmsnSHE,surf,surfp,surfn
            CALL FLUSH(702)
!          END IF
          END IF
        END IF
        IF(isample==nsamples) EXIT
      END DO
      IF(isample<nsamples) WRITE(702,*) "ensemble not completely full"
      WRITE(702,'(1x/1x)')
      WRITE(702,'(a)') '# probability distribution:'
      DO i=1,21
        WRITE(702,'(f5.2,1pg13.5)') i*0.1D0-0.05D0,binprob(i)/nsamples
      END DO
!      WRITE(702,'(1x/1x)')
!      WRITE(702,'(a)') '# RAN1:'
!      DO i=-10,10
!        WRITE(702,'(f5.2,1pg13.5)') i*0.1D0,bintest(i)
!      END DO
      CLOSE(702)
      RETURN
      END


