c-----matin2----------------------------------------------------------
c
      subroutine matin2 (a2,nmax,ndim,nerror,determ)
c
      include 'params.inc'
      parameter (kindex=500)
c
c     inverts the matrix 'c'.
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
c       the routine is taken from 'cernlib' and modified for this
c       package.
c
      dimension a2(ndim,ndim),index(kindex)
      real*16 a(ndim,ndim)

      do i=1,nmax
      do j=1,nmax
        a(i,j) = a2(i,j)
      enddo
      enddo
 
      call  matin4 (a,nmax,ndim,nerror,determ)

      do i=1,nmax
      do j=1,nmax
        a2(i,j) = a(i,j)
      enddo
      enddo

      return
      end

c
c-----matin4----------------------------------------------------------
c
      subroutine matin4 (a,nmax,ndim,nerror,determ)
c
      implicit real*16 (a-h,o-z)
      parameter (kindex=500)
c
c     inverts the matrix 'c'.
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
c       the routine is taken from 'cernlib' and modified for this
c       package.
c
      real*16 a(*)
      integer index(kindex)
      real*8 determ
c
c-----------------------------------------------------------------------------
c
      if(nmax.gt.kindex) stop ' too large dimension in matin2'
c
      deter=1D0
      n=nmax
c     iemat=n+n2
      iemat=n+0
      idim=nmax
      nmin1=n-1
c
c   umsetzen des nmax*nmax-anteils der matrix a in gepackte
c   ein-dimensionale form.
c
      do 400 i=2,idim
      do 400 j=1,idim
      ipack  = (i-1)*nmax + j
      imatr  = (i-1)*ndim + j
      a(ipack) = a(imatr)
  400 continue
c        the routine does its own evaluation for double subscripting of
c        array a.
      ipivc=1-idim
c        main loop to invert the matrix
      do 11 main=1,n
      pivot=0
      ipivc=ipivc+idim
c        search for next pivot in column main.
      ipivc1=ipivc+main-1
      ipivc2=ipivc +nmin1
      do 2 i1=ipivc1,ipivc2
      if(abs(a(i1))-abs(pivot)) 2,2,1
    1 pivot=a(i1)
      lpiv=i1
    2 continue
c        is pivot different from zero
      if(pivot) 3,15,3
c        get the pivot-line indicator and swap lines if necessary
    3 icol=lpiv-ipivc+1
      index(main)=icol
      if(icol-main) 6,6,4
c        complement the determinant
    4 deter=-deter
c        pointer to line pivot found
      icol=icol-idim
c        pointer to exact pivot line
      i3=main-idim
      do 5 i=1,iemat
      icol=icol+idim
      i3=i3+idim
      swap=a(i3)
      a(i3)=a(icol)
    5 a(icol)=swap
c        compute determinant
    6 deter=deter*pivot
      pivot=1/pivot
c        transform pivot column
      i3=ipivc+nmin1
      do 7 i=ipivc,i3
    7 a(i)=-a(i)*pivot
      a(ipivc1)=pivot
c        pivot element transformed
c
c        now convert rest of the matrix
      i1=main-idim
c        pointer to pivot line elements
      icol=1-idim
c        general column pointer
      do 10 i=1,iemat
      icol=icol+idim
      i1=i1+idim
c        pointers moved
      if(i-main) 8,10,8
c        pivot column excluded
    8 jcol=icol+nmin1
      swap=a(i1)
      i3=ipivc-1
      do 9 i2=icol,jcol
      i3=i3+1
    9 a(i2)=a(i2)+swap*a(i3)
      a(i1)=swap*pivot
   10 continue
   11 continue
c        now rearrange the matrix to get right invers
      do 14 i1=1,n
      main=n+1-i1
      lpiv=index(main)
      if(lpiv-main) 12,14,12
   12 icol=(lpiv-1)*idim+1
      jcol=icol+nmin1
      ipivc=(main-1)*idim+1-icol
      do 13 i2=icol,jcol
      i3=i2+ipivc
      swap=a(i2)
      a(i2)=a(i3)
   13 a(i3)=swap
   14 continue
      determ=deter
      nerror=0
c   ruecksetzen der gepackten eindimensionalen form in den
c   nmax*nmax-anteil der matrix c
c
      do 500 i=idim,2,-1
      do 500 j=idim,1,-1
      ipack  = (i-1)*nmax + j
      imatr  = (i-1)*ndim + j
  500 a(imatr) = a(ipack)
c
      return
   15 nerror=main
      determ=deter
      return
      end
