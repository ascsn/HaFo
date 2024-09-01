      subroutine correl_data(nterms)

      include 'params.inc'
      include 'numpar.inc'
      include 'comfit.inc'

      INTEGER,INTENT(IN) :: nterms

      REAL(8),ALLOCATABLE :: corr_data(:,:),rootc(:)
      REAL(8),ALLOCATABLE :: vect(:,:)
      INTEGER,ALLOCATABLE :: ipoi(:)
      LOGICAL,PARAMETER :: tcond=.true.

      nall=ndats*npts
      ALLOCATE(ipoi(nall))
      ipoi=0

      iact=0
      ilauf=0
      do idat=1,ndats
      do ipar=1,npts
        ilauf=ilauf+1
        IF(chi(idat,ipar) > 1D-2) THEN
          iact = iact+1
          ipoi(iact) = ilauf
        END IF
      END DO
      END DO
      nall = iact

      WRITE(2,*) 'CORREL_DATA: nterms,nall=',nterms,nall,ndats*npts
      ALLOCATE(corr_data(nall,nall),vect(nall,nall),rootc(nall))

      DO i=1,nall; DO j=1,nall
        acc = 0D0
        DO n=1,nterms
          acc = acc + fderij(ipoi(i),n)*
     &                SUM(zarray(n,:)*fderij(ipoi(j),:))
        END DO
        corr_data(i,j) = acc
      END DO; END DO

      WRITE(2,*)  'CORR_DATA computed, trace=',SUM(corr_data(:,:))

      IF(tcond) THEN
        DO i=1,nall; DO j=1,nall
          IF(corr_data(i,i)<1D-20 .OR. corr_data(j,j)<1D-20) THEN
             corr_data(i,j) = 0D0
          ELSE
             vect(i,j) = corr_data(i,j)/
     &                   SQRT(corr_data(i,i)*corr_data(j,j))
          END IF
        END DO; END DO
        corr_data=vect
        WRITE(2,*)  'CORR_DATA conditioned, trace=',SUM(corr_data(:,:))
      ELSE
        WRITE(2,*) 'CORR_DATA diagonal elements'
        DO i=1,nall
          WRITE(2,'(i4,1pg13.5)') i,corr_data(i,i)
        END DO      
      END IF

      call diag_real_symm(corr_data,rootc,vect,
     &                    nall,nall,(nall*nall+nall)/2)

      WRITE(2,'(a)') 'roots of data correlations:'
      WRITE(2,'(i4,1pg13.5)') (n,rootc(nall-n+1),n=1,nall)

      accroot = 0D0
      DO n=nall-1,1,-1
        IF(rootc(n)<1D-20) EXIT
        rootc(n) = rootc(n)+rootc(n+1)
        accroot=rootc(n)
      END DO
      WRITE(2,*) 'cumulated roots:'
      DO n=nall,1,-1
        IF(rootc(n)<1D-20) EXIT
        WRITE(2,'(i4,1pg13.5)') n,rootc(n)/accroot
      END DO
      DEALLOCATE(corr_data,vect,rootc)

      return
      end
