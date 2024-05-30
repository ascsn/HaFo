c-----evaluate_deform------------------------------------------------
c
      program evaluate_deform
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
        open(31,file='crankn.'//ident)
ctest        write(6,'(a,i5)')
ctest     &   'ensurf.'//ident,i
        
        if(i.eq.0) then
          call read_ensurf(emfsav,eprojsav,31)
        else
          call read_ensurf(emf,eproj,31)
          workdv(i,1) = (emf-emfsav)/dparmo(i)
          workdv(i,2) = (eproj-eprojsav)/dparmo(i)
        endif
        close(31)
c
      enddo
c
c     compute variances
c
      errmf    = 0.0
      errprj   = 0.0
c
      if(nparmx.gt.0) then
        do i=1,nparmx
        do j=1,nparmx
          errmf   = workdv(i,1)*workdv(j,1)*zarray(i,j)
     &              + errmf 
          errprj  = workdv(i,2)*workdv(j,2)*zarray(i,j)
     &              + errprj
          errehs  = workdv(i,3)*workdv(j,3)*zarray(i,j)
     &              + errehs
        enddo
        enddo
        errmf   = sqrt(errmf )
        errprj  = sqrt(errprj)
        errehs  = sqrt(errehs)
      endif
c
c     print
c
      write(6,'(a,2(1pg13.5)/a,2(1pg13.5))')
     & ' binding energy 168Er - mean-field:',emf,errmf ,
     & ' binding energy 168Er - projected: ',eproj,errprj
c
      stop
      end

c-----read_ensurf----------------------------------------------------
c
      subroutine read_ensurf(emf,eproj,iunit)
c
c     reads information from file and 
c     extracts fission barrier.
c       iunit      = unit number for file with input data
c       barrdir    = barrier from maximum of listed values
c       barrint    = barrier from interpolated maximum
c
      implicit double precision (a-h,o-z)
c
      parameter (KMAXEN=10)
c
      dimension e(KMAXEN)
c      character*80 a
c
c-------------------------------------------------------------------
c
      do i =1,9999
        read(iunit,'(1x)',err=19,end=19)
      enddo
 19   continue
      backspace iunit
      read(iunit,*,err=99,end=99) 
     &   x1,x2,x3,x4,ammass,amwid,xx7,energ,ecm
c
      emf = energ-ecm
      eproj = energy-0.5*amwid/ammass
c
      return
c 
 99   stop ' premature end of input file'
c
      end
