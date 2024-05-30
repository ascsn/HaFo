c-----evaluate_barriers------------------------------------------------
c
      program evaluate_barriers
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
      nparmx = 0
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
          call read_ensurf(barrdirsav,barrintsav,ehssav,31)
        else
          call read_ensurf(barrdir,barrint,ehs,31)
          workdv(i,1) = (barrdir-barrdirsav)/dparmo(i)
          workdv(i,2) = (barrint-barrintsav)/dparmo(i)
          workdv(i,3) = (ehs-ehssav)/dparmo(i)
        endif
        close(31)
c
      enddo
c
c     compute variances
c
      errdir   = 0.0
      errint   = 0.0
      errehs   = 0.0
c
      if(nparmx.gt.0) then
        do i=1,nparmx
        do j=1,nparmx
          errdir  = workdv(i,1)*workdv(j,1)*zarray(i,j)
     &              + errdir
          errint  = workdv(i,2)*workdv(j,2)*zarray(i,j)
     &              + errint
          errehs  = workdv(i,3)*workdv(j,3)*zarray(i,j)
     &              + errehs
        enddo
        enddo
        errdir  = sqrt(errdir)
        errint  = sqrt(errint)
        errehs  = sqrt(errehs)
      endif
c
c     print
c
      write(6,'(a,2(1pg13.5)/a,2(1pg13.5)/a,2(1pg13.5))')
     & ' fission barrier in 120/182 - direct:',barrdirsav,errdir,
     & ' fission barrier in 120/182 - inter.:',barrintsav,errint,
     & ' binding energy in 264Hs:            ',ehssav,errehs
c
      stop
      end

c-----read_ensurf----------------------------------------------------
c
      subroutine read_ensurf(barrdir,barrint,ehs,iunit)
c
c     reads information from file and 
c     extracts fission barrier.
c       iunit      = unit number for file with input data
c       barrdir    = barrier from maximum of listed values
c       barrint    = barrier from interpolated maximum
c
      implicit double precision (a-h,o-z)
c
      parameter (KMAXEN=5)
c
      dimension e(KMAXEN)
c      character*80 a
c
c-------------------------------------------------------------------
c
      do i =1,8
        read(iunit,'(1x)')
c        read(iunit,'(a)') a
c        write(6,'(a)') a
      enddo
      read(iunit,*,err=99,end=99) x1,x2,x3,x4,x5,x6,x7,energ0 
      do i =1,8
        read(iunit,'(1x)')
      enddo
      maxen = 0
      emax  = -1.0D33
      imax  = 0
      do i=1,KMAXEN
        read(iunit,*,err=19,end=19) x1,x2,x3,x4,x5,x6,x7,e(i) 
        maxen = i
ctest        write(6,'(i,1pg13.5)') i,e(i)
        if(e(i).gt.emax) then
          emax = e(i)
          imax = i
        endif
      enddo
      do i =1,8
        read(iunit,'(1x)')
      enddo
      read(iunit,*,err=19,end=19) x1,x2,x3,x4,x5,x6,x7,ehs 
c
 19   continue
      barrdir = emax-energ0
      imax    = min(max(imax,2),maxen-1)  
      d1e     = (e(imax+1)-e(imax-1))/2.0D0
      d2e     = e(imax+1)-2.0D0*e(imax)+e(imax-1)
      x       = -d1e/d2e
      barrint = e(imax)+x*d1e+x*x*d2e/2.0D0-energ0
ctest      write(6,'(i5,4(1pg13.5))')
ctest     &  imax,energ0,d1e,d2e,x
c
      return
c 
 99   stop ' premature end of input file'
c
      end
