c-----evaluate_barriers------------------------------------------------
c
      program evaluate_bind
c
      implicit double precision (a-h,o-z)
c
      include 'comfit.inc'
c
      character*3 ident
c
      dimension par(kparmx),workdv(kparmx,3)
c
c---------------------------------------------------------------------
c
c     read 'dparmo' and 'zarray' from data file
c
      open(31,file='zarray.dat')
      read(31,'(i5)') nparmx
      read(31,'(17g14.6)') (dparmo(j),j=1,nparmx)
      do i=1,nparmx
        read(31,'(17g14.6)') (zarray(j,i),j=1,nparmx)
      enddo
      close (31)
c
c      nparmx = 0
c
c     read barriers from 'ensurf.*'
      do i=0,nparmx
c
c       prepare identifier
c
        open(31,status='scratch')
        write(31,'(i3)') i+100
        rewind(31)
        read(31,'(a3)') ident
        close(31)
        ident(1:1) = 's' 
c
c       read and store
c
        open(31,file='ensurf.'//ident)
ctest        write(6,'(a,i5)')
ctest     &   'ensurf.'//ident,i
        if(i.eq.0) then
          call read_ensurf(ehssav,31)
        else
          call read_ensurf(ehs,31)
          workdv(i,3) = (ehs-ehssav)/dparmo(i)
        endif
        close(31)
c
      enddo
c
c     compute variances
c
      errdir   = 0.0
      errehs   = 0.0
c
      if(nparmx.gt.0) then
        do i=1,nparmx
        do j=1,nparmx
          errehs  = workdv(i,3)*workdv(j,3)*zarray(i,j)
     &              + errehs
        enddo
        enddo
        errehs  = sqrt(errehs)
      endif
c
c     print
c
      write(6,'(a,2(1pg13.5)/a,2(1pg13.5)/a,2(1pg13.5))')
     & ' binding energy in 252No:',ehssav,errehs
c
      stop
      end

c-----read_ensurf----------------------------------------------------
c
      subroutine read_ensurf(ehs,iunit)
c
c     reads information from file and 
c     extracts fission barrier.
c       iunit      = unit number for file with input data
c
      implicit double precision (a-h,o-z)
c
      parameter (KMAXEN=10)
      parameter (EBalpha=-28.2)
c
      dimension e(KMAXEN)
c      character*80 a
c
c-------------------------------------------------------------------
c
      do i =1,12
        read(iunit,'(1x)')
c        read(iunit,'(a)') a
c        write(6,'(a)') a
      enddo
      read(iunit,*,err=99,end=99) x1,x2,x3,x4,x5,x6,energ0 
      ehs = energ0
      write(6,'(a,3(1pg13.5))') 
     &  ' E(266),E(262),Qalpha=',energ0,energ1,qalpha
c
      return
c 
 99   stop ' premature end of input file'
c
      end
