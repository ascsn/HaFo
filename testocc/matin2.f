c-----matin2----------------------------------------------------------
c
      subroutine matin2 (a,nmax,ndim,nerror,determ)
c
      include 'params.inc'
      parameter (kindex=500)
      logical tunclean
      parameter (tunclean=.true.)
c
c     inverts the matrix 'a'.
c
c     the parameters on the i/o list are:
c
c       a       = at input:  matrix to be inverted
c                 at output: inverse matrix
c       nmax    = number of actual rows and columns
c       ndim    = dimension of matrix in the calling routine
c       nerror  = error flag
c       determ  = determinant of input matrix 'a'
c
c
      dimension a(ndim,ndim),work(ndim,ndim)
      dimension asav(ndim,ndim),alud(ndim,ndim),rhs(ndim)
      dimension index(ndim)
c
c-----------------------------------------------------------------------------
c
      do i=1,nmax
        do j=1,nmax
          work(i,j) = 0D0
          asav(j,i) = a(j,i)
        enddo
        work(i,i) = 1D0
      enddo

      call ludcmp(a,nmax,ndim,index,determ)

      do j=1,nmax
        call lubksb(a,nmax,ndim,index,work(1,j))
      enddo

      do i=1,nmax
      do j=1,nmax
        alud(i,j) = a(i,j)
        a(i,j) =  work(i,j)
      enddo
      enddo

      if(tunclean) return

c
c     clean solution iteratively
c
      do i=1,nmax
        do j=1,nmax
          rhs(j) = 0D0
          work(j,i) = a(j,i)
        enddo
        rhs(i) = 1D0
        call mprove(asav,alud,nmax,ndim,index,rhs,work(1,i))
      enddo

   
      do i=1,nmax
      do j=1,nmax
        a(i,j) = work(i,j)
      enddo
      enddo


      return
      end


c-----from numerical recipes
c
      subroutine mprove(a,alud,n,np,indx,b,x)
      parameter (nmax=100)
      implicit double precision (a-h,o-z)
      dimension a(np,np),alud(np,np),indx(n),b(n),x(n),r(nmax)
      real*8 sdp
      
      do i=1,n
        sdp =-b(i)
        do j=1,n
          sdp = sdp + a(i,j)*x(j)
        enddo
        r(i) = sdp
      enddo
      call lubksb(alud,n,np,indx,r)
      do i=1,n
        x(i) = x(i) - r(i)
      enddo

      return
      end

c
      SUBROUTINE ludcmp(a,n,np,indx,d)
      implicit double precision (a-h,o-z)
      INTEGER n,np,indx(n),NMAX
      DOUBLE PRECISION d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      DOUBLE PRECISION aamax,dum,sum,vv(NMAX)
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
      SUBROUTINE lubksb(a,n,np,indx,b)
      implicit double precision (a-h,o-z)
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
