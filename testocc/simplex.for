c-----simplex---------------------------------------------------------
c
      subroutine simplex(a,nterms,nwdr)
c
c     chi**2 minimization by simplex method
c     call the 'amoeba' from numerical recipes
c
c     dimension statement valid for nterms up to 'kparmx'
c
c
      include 'params.inc'
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'comfit.inc'
c
      external ffverg
c
      dimension a(kparmx)
c
c     workspaces
c
      dimension chisav(0:kparmx),parsav(0:kparmx,kparmx)
      dimension parmom(kparmx)
      equivalence (chisav,zarray),(parsav,fderij),(chimom,parmom)
c
      data  eps/0.00001/
c
c----------------------------------------------------------------------
c
c     initialize first guess from paramaters on 'a' and
c     typical variations on 'dparmo'
c
      write(2,'(a,i5,17(1pg12.4))')
     &   ' Init SIMPLEX:  jpar   params=',
     &   0,(a(ipar),ipar=1,nterms)
      chisav(0) = ffverg(a)
      write(6,'(a,1pg12.4)')  ' chi**2=',chisav(0)
      write(2,'(a,1pg12.4)')  ' chi**2=',chisav(0)
      do i=1,kparmx
       parsav(0,i) = a(i)
       parmom(i)   = a(i)
      enddo
      do jpar=1,nterms
        parmom(jpar)=a(jpar)+dparmo(jpar)
        write(2,'(a,i5,17(1pg12.4))')
     &   ' SIMPLEX:  jpar   params=',
     &   jpar,(parmom(ipar),ipar=1,nterms)
        chisav(jpar) = ffverg(parmom)
        write(2,'(a,1pg12.4)') ' chi**2=',chisav(jpar)
        do i=1,kparmx
         parsav(jpar,i)=parmom(i)
        enddo
        parmom(jpar)=a(jpar)-dparmo(jpar)
      enddo
      write(2,'(a/a,30(/1x,i5,1pg12.4))')
     &   ' Initialize SIMPLEX with:','  jpar   chi**2 ',
     &   (jpar,chisav(jpar),jpar=0,nterms)
c
c     simplex iteration by modified 'amoeba'
c
      call amoeba(parsav,chisav,kparmx+1,kparmx,nterms,eps,ffverg,iter,
     &             abs(nwdr))
c
c     pick optimal solution
c
      chimin = +1.0E10
      do jpar=0,nterms
        if(chisav(jpar).lt.chimin) then
          jmin = jpar
        endif
      enddo
      do i=1,kparmx
        a(i) = parsav(jmin,i)
      enddo
c
      return
c
      end
c-----amoeba-----------------------------------------------------------
c
      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter,itmax)
c
c     from numerical recipes, slightly modifed,
c     version for double precision
c
      INTEGER iter,mp,ndim,np,NMAX,ITMAX
      DOUBLE PRECISION ftol,p(mp,np),y(mp),funk
c      PARAMETER (NMAX=20,ITMAX=5000)
      PARAMETER (NMAX=40)        ! must be >= KPARMX
      EXTERNAL funk
CU    USES amotry,funk
      INTEGER i,ihi,ilo,inhi,j,m,n
      DOUBLE PRECISION rtol,sum,swap,ysave,ytry,psum(NMAX),amotry
      REAL(8) :: chisav,parwrit(NMAX)
      chisav=1D99
      if(np.gt.nmax) stop ' insufficient dimension in AMOEBA'
      do j=1,np
        psum(j) = p(1,j)
      enddo
      iter=0
1     continue
      do 12 n=1,ndim
        sum=0.d0
        do 11 m=1,ndim+1
          sum=sum+p(m,n)
11      continue
        psum(n)=sum
12    continue
2     ilo=1
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 13 i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        endif
13    continue
      rtol=2.d0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))
      if (rtol.lt.ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do 14 n=1,ndim
          swap=p(1,n)
          p(1,n)=p(ilo,n)
          p(ilo,n)=swap
14      continue
        return
      endif
      if (iter.ge.ITMAX) return
      iter=iter+2
      write(2,'(a,i5)') ' ITER=',iter
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0d0)
        if(ytry.lt.chisav) then
          write(6,'(a,2g12.4)') 
     &     ' update FOR003. chi old/new=',ytry,chisav
          chisav = ytry
          DO j=1,np
            parwrit(j) = p(1,j)
          END DO
          CALL write_params(parwrit,ytry,iter)
!          rewind(3)
!          do j=1,ndim
!           write(3,'(g15.8,a,e12.5,a,i2,a,a)') p(1,j)
!          enddo
!          write(3,'(a,1pg13.5)') ' chi**2=',ytry
!          call flush(3)
        endif

      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0d0)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5d0)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5d0*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              write(2,'(a,i5,17(1pg14.6))')
     &        ' in AMOEBA: i,params=',i,(psum(ipar),ipar=1,ndim)
              write(6,'(a,i5,17(1pg12.4))')
     &        ' in AMOEBA: i,params=',i,(psum(ipar),ipar=1,ndim)
              y(i)=funk(psum)
              write(6,'(a,1pg12.4)')
     &        ' chi**2=',y(i)
              write(2,'(a,1pg12.4)')
     &        ' chi**2=',y(i)
            endif
16        continue
          iter=iter+ndim
          goto 1
        endif
      else
        iter=iter-1
      endif
      goto 2
      END
      FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
      INTEGER ihi,mp,ndim,np,NMAX
      DOUBLE PRECISION amotry,fac,p(mp,np),psum(np),y(mp),funk
      PARAMETER (NMAX=20)
      EXTERNAL funk
CU    USES funk
      INTEGER j
cold      DOUBLE PRECISION fac1,fac2,ytry,ptry(NMAX)
      DOUBLE PRECISION fac1,ytry,ptry(np)
!      WRITE(2,'(a,3i5)') 'AMOTRY: mp,np,ndim=',mp,np,ndim
      do j=1,np
        ptry(j) = p(1,j)
      enddo
      fac1=(1.d0-fac)/(ndim+1)
c      fac2=fac1-fac
      do 11 j=1,ndim
c        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
        ptry(j)=psum(j)*fac1+p(ihi,j)*fac
        write(6,'(a,2i5,3(1pg12.4))')
     &   ' j,fac,fac1,p,psum,ptry=',j,fac,fac1,p(ihi,j),psum(j),ptry(j)
!        write(2,'(a,2i5,3(1pg12.4))')
!     &   ' j,fac,fac1,p,psum,ptry=',j,fac,fac1,p(ihi,j),psum(j),ptry(j)
11    continue
      write(6,'(a,17(1pg12.4))')
     &   ' in AMOTRY: params=',(ptry(ipar),ipar=1,ndim)
      write(2,'(a,17(1pg14.6))')
     &   ' in AMOTRY: params=',(ptry(ipar),ipar=1,ndim)
      ytry=funk(ptry)
      write(6,'(a,1(1pg12.4))')
     &   ' chi**2=',ytry
      write(2,'(a,1(1pg12.4))')
     &   ' chi**2=',ytry
      if (ytry.lt.y(ihi)) then
        y(ihi)=ytry
        do 12 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)
12      continue
      endif
      amotry=ytry
!      WRITE(2,*) 'AMOTRY finished'
      return
      END
c-----write_params-------------------------------------------------------------
c
      SUBROUTINE write_params(paract,chiact,iter)
c
      INCLUDE 'params.inc'
      INCLUDE 'comfit.inc'
      INCLUDE 'skypar.inc'

      REAL(8), INTENT(IN) :: paract(kparmx),chiact
      INTEGER, INTENT(IN) :: iter

      CHARACTER(6) :: zuord(kparmx)

      CALL naming(zuord,ioptsk)
      REWIND(3)
      DO j=1,nparin
        WRITE(3,'(1pg20.12,a,0pg12.5,a,i2,a,a)')  
     &        paract(j),',',dparmo(j),',',
     &        iparva(j),',',zuord(iparva(j))
      END DO
      WRITE(3,'(a/a,i4,1pg13.5)')
     &     '-11, "to stop"',
     &     ' iteration,chi**2=',iter,chiact
      CALL flush(3)

      RETURN
      END
