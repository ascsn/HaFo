      include 'params.inc'
c
      include 'comfit.inc'
c
      logical tdiagz
      parameter (tdiagz=.true.)
      dimension a(kparmx),ares(kparmx),sigmaa(kparmx),sigma0(kparmx)
      real*8 alpha(kparmx,kparmx)
      dimension b(kparmx),beta(kparmx)
      dimension vectz(kparmx,kparmx),rootz(kparmx),diagal(KPARMX)
      dimension alphdr(kparmx),betadr(kparmx)
      real*16 sum


      open(34,form='unformatted',file='alpha.dat')
      read(34) nterms
      read(34) ((alpha(j,k),j=1,nterms),k=1,nterms)
      close(34)

       if(tdiagz) then
         call diag_real_symm(alpha,rootz,vectz,
     &                       nterms,kparmx,(nterms*nterms+nterms)/2)
         write(6,'(a)') ' roots of ZARRAY before conditioning:'
         write(6,'(5(1pg13.5))') (rootz(i),i=1,nterms)
       endif

      do j=1,nterms
        diagal(j) = alpha(j,j)
      enddo
      

      do j=1,nterms
      do k=1,nterms
c       zarray(j,k)=alpha(j,k)/sqrt(alpha(j,j)*alpha(k,k))
       zarray(j,k)=alpha(j,k)/sqrt(diagal(j)*diagal(k))
c       zarray(j,k)=alpha(j,k)
      enddo
      enddo

      write(6,'(a)') ' conditioned matrix ZARRY before inversion'
      offmax = 0D0
      do j=1,nterms
        do k=1,nterms
          if(j.ne.k) offmax = max(offmax,abs(zarray(j,k)))
        enddo
        write(6,'(5(1pg14.6))') (zarray(k,j),k=1,nterms)
      enddo
      write(6,'(a,1pg13.5)') ' maximum off diagonal=',offmax
      write(6,'(1x)')


      call matin2(zarray,nterms,kparmx,nerr,det)

      do j=1,nterms
      do k=1,nterms
        zarray(j,k)=zarray(j,k)/sqrt(diagal(j)*diagal(k))
      enddo
      enddo

      if(tdiagz) then
         call diag_real_symm(zarray,beta,vectz,
     &                       nterms,kparmx,(nterms*nterms+nterms)/2)
         write(6,'(a)') ' roots of ZARRAY after conditioning:'
         write(6,'(5(1pg13.5))') (beta(i),i=1,nterms)
         write(6,'(a)') ' test product of roots:'
         write(6,'(5(1pg13.5))') 
     &    (beta(i)*rootz(nterms-i+1)-1D0,i=1,nterms)
      endif

      write(6,'(a)') ' test inversion by multiplication:'
      errmax = 0D0
      do j=1,nterms
        do k=1,nterms
          sum = 0
          do l=1,nterms
            sum = sum + zarray(j,l)*alpha(l,k)
          enddo
          rootz(k) = sum
          if(j.eq.k) sum = sum - 1
          errmax = max(errmax,abs(sum))
        enddo
        write(6,'(5(1pg14.6))') (rootz(k),k=1,nterms)
      enddo
      write(6,'(a,1pg13.5)') ' maximum deviation=',errmax

      stop
      end


      include 'matin2old.f'
      include 'givens.f'
