c-----jdata -----------------------------------------------------------
c
      program jdata
c
c     extracts the key observables from the FB data file of Joerg Friedrich
c
      parameter (kmaxsv=100)
      character*194 asave(0:kmaxsv)
      character*100 a
      character*8   akey,aperce,apercr,apercs,aprms
      dimension     isort(kmaxsv),ipoint(kmaxsv)
c
c----------------------------------------------------------------------
c
      write(6,'(a/a)')
     & ' compressed output of experimental evaluation of radii:',
     & '    A   Z   r_msg   r_cmp  r_diff   sigma  halo_m  halo_c'
      open(19,form='formatted',status='scratch')
      iset   = 0
      isetmx = 0
      do 10 iread=1,9999
      read(5,'(a)',err=999,end=999) a
      if(a(1:6).eq.'   A =') then
c                                                        read nucleus
        iset  = 1 + iset
        write(19,'(a)') a
        backspace 19
        read(19,'(8x,f8.2,16x,f8.2)')  anumb,znumb
        backspace 19
c                                                        find rms
        do 20 ir2=1,999
          read(5,'(a)',err=999,end=999) a
          if(a(1:13).eq.'   rms-Radius') then
           write(19,'(a)') a
           backspace 19
           read(19,'(23x,f8.3,15x,f8.3)') rmsmsg,rmscmp
           backspace 19
           goto 29
          endif
   20   continue
   29   continue
c                                                        find R_diffr
        do 30 ir2=1,999
          read(5,'(a)',err=999,end=999) a
          if(a(1:12).eq.' NR.   Q MIN') then
           read(5,'(1x/1x)')
           read(5,*) nr,qmin,fq,rdiffr
           goto 39
          endif
   30   continue
   39   continue
c                                                        find sigma
        do 40 ir2=1,999
          read(5,'(a)',err=999,end=999) a
          if(a(29:36).eq.'SIGMA(1)') then
           read(5,'(1x,4(/1x))')
           read(5,*) nr,qmax,fq,sigma
           goto 49
          endif
   40   continue
   49   continue
c                                                        print
        ia   = anumb+0.000001
        iz   = znumb+0.000001
        haloms = rmsmsg-sqrt(0.6*rdiffr**2+3.0*sigma**2)
        halocm = rmscmp-sqrt(0.6*rdiffr**2+3.0*sigma**2)
        write(6,'(1x,2i4,6f8.3)')
     &   ia,iz,rmsmsg,rmscmp,rdiffr,sigma,haloms,halocm
      endif
   10 continue
  999 continue
c
      end

