c-----sensivity---------------------------------------------------------
c
      subroutine sensitivity(alpha,vectz,rootz,nterms,ndatsact)
c
!     Sensitivity analysis      
!     Switches to include difference observables with 'tdiffobs' from
!     'comfit.inc'.
c
      include 'params.inc'
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'comfit.inc'
c
      REAL(8),INTENT(IN) :: alpha(kparmx,kparmx)
      REAL(8),INTENT(IN) :: vectz(kparmx,kparmx),rootz(kparmx)
      INTEGER,INTENT(IN) :: nterms,ndatsact

      REAL(8),ALLOCATABLE :: Sdet(:,:)   ! detailed sensitivity 
      REAL(8) :: rootinv(kparmx)
      REAL(8),ALLOCATABLE :: Saux(:,:),Ssav(:,:),covariance(:,:)
      REAL(8),ALLOCATABLE :: sensgroup(:)
      REAL(8),PARAMETER :: xlow=1D-20
      REAL(8) :: acc
      LOGICAL,PARAMETER :: ttest=.FALSE.
      INTEGER :: j1,j2,icount
c
c
c---------------------------------------------------------------
c
      ALLOCATE(Sdet(kderij,kparmx),Saux(kderij,kparmx))
      ALLOCATE(sensgroup(1:kdat),covariance(kparmx,kparmx))
      IF(ttest) THEN
        WRITE(2,'(a,3i5)') 'SENSITIVITY: ndats,npts,nterms=',
     &                  ndatsact,npts,nterms
      END IF
!     raw sensitivity
      IF(ttest) WRITE(2,'(a)') 'FDERIJ,SDET:'
      IF(tdiffobs) THEN
        ilauf=0
        do idat=1,ndatsact
        do inuc=1,npts
          ilauf=ilauf+1
          do jpar=1,nterms
             Sdet(ilauf,jpar)=fjacobian(ilauf,jpar)
          enddo
        enddo
        enddo
      ELSE
        ilauf=0
        do idat=1,ndatsact
        do inuc=1,npts
          ilauf=ilauf+1
          do jpar=1,nterms
             Sdet(ilauf,jpar)=fderij(ilauf,jpar)/sigmay(idat,inuc)
 !            IF(ttest) WRITE(2,'(4i5,2(1pg13.5))') 
 !    &        idat,inuc,ilauf,jpar,Sdet(ilauf,jpar),fderij(ilauf,jpar)
          enddo
        enddo
        enddo
      END IF
!     scale to deal with conditioned covarince matrix
        ilauf=0
        do idat=1,ndatsact
        do inuc=1,npts
          ilauf=ilauf+1
          do jpar=1,nterms
             Sdet(ilauf,jpar)=Sdet(ilauf,jpar)/SQRT(alpha(jpar,jpar))
          enddo
        enddo
        enddo
!     transform to covarianc*sensitivity
      rootinv(1:nterms)=1D0/rootz(1:nterms)
      WRITE(2,'(a/200(1pg13.5))') 'inverse roots:',rootinv(1:nterms)
      Ssav=Sdet
      DO jcut=1,4
      Sdet=Ssav

!     covariance matrix  (misuse Saux)
      DO j1=1,nterms; DO j2=1,nterms
        covariance(j1,j2)= 
     &     SUM(vectz(j1,1:nterms)*rootinv(1:nterms)*vectz(j2,1:nterms))/
     &       SQRT(alpha(j1,j1)*alpha(j2,j2))
      END DO; END DO
!     sensitivity matrix
      ilauf=0
      do idat=1,ndatsact
      do inuc=1,npts
        ilauf=ilauf+1
!        IF(ttest) WRITE(2,'(a,3i5)') 'SENSITIVITY:',idat,inuc,ilauf
        do jroot=1,nterms
           Saux(ilauf,jroot)= 
     &        SUM(vectz(1:nterms,jroot)*Sdet(ilauf,1:nterms))*
     &        rootinv(jroot)
        enddo
        do jpar=1,nterms
           Sdet(ilauf,jpar)= 
     &        SUM(vectz(jpar,1:nterms)*Saux(ilauf,1:nterms))/
     &           SQRT(alpha(jpar,jpar))
        enddo
      enddo
      enddo
!     sum within groups
      WRITE(2,'(a,i3)') 
     &   '#accumulated sensitivities per data group: cut last',jcut-1
!     &   ' param.   sensitivities'
      do jpar=1,nterms
        ilauf=0
        sensgroup=0D0
        do idat=1,ndatsact
          do inuc=1,npts
            ilauf=ilauf+1
            sensgroup(idat)=sensgroup(idat)+Sdet(ilauf,jpar)**2
          enddo
        enddo
        sensall=SUM(sensgroup)
        sensgroup=sensgroup/sensall
        DO idat=ndatsact,1,-1
          sensgroup(idat)=SUM(sensgroup(1:idat))
        END DO
        WRITE(2,'(i3,200(1pg12.4))') jpar,sensgroup(1:ndatsact)
      enddo
      WRITE(2,'(1x/1x)')
      rootinv(nterms-jcut+1)=0D0

!     compute and print accumulated sensitivty for each data point
      WRITE(2,'(a/200(1pg13.5))') 'covariance:',
     &    (SQRT(covariance(j1,j1)),j1=1,nterms)
      WRITE(2,'(a)') 'parameter averaged sensivity for data point'
      ilauf=0
      icount=0
      do idat=1,ndatsact
        do inuc=1,npts
          ilauf=ilauf+1
          DO jpar=1,nterms
            Sdet(ilauf,jpar)=Sdet(ilauf,jpar)**2/covariance(jpar,jpar)
          END DO
          IF(SUM(Sdet(ilauf,:))>0D0) THEN
            icount=icount+1
            WRITE(2,'(3i4,1pg13.5)') icount,idat,inuc,
     &            SUM(Sdet(ilauf,1:nterms))
          END IF
      enddo
      enddo

!     compute and print relative sensitivty for each data point
      WRITE(2,'(a)') 'relative sensivity for data point'
      ilauf=0
      icount=0
      do idat=1,ndatsact
        do inuc=1,npts
          ilauf=ilauf+1
          IF(SUM(Sdet(ilauf,:))>0D0) THEN
            icount=icount+1
            acc = SUM(Sdet(ilauf,1:nterms))
            WRITE(2,'(3i4,200(1pg12.4))') icount,idat,inuc,
     &            (SUM(Sdet(ilauf,1:jpar))/acc,jpar=1,nterms)
          END IF
      enddo
      enddo

      END DO     ! jcut loop
!
      DEALLOCATE(Sdet,Saux,sensgroup,covariance)

!
      RETURN
      END

