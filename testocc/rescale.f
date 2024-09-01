      character*80 a
      read(5,*)
     &   lrpa,resfac,cmfac,cutmod,ifsepi,resdmp,iforth,ecut
      resdmp = 2.7
      ecut   = 200.0
      cmfac  = 0.005
      write(6,'(i2,a,f3.1,a,f6.4,a,g8.1,a,i1,a,f3.1,a,i1,
     &        a,f4.0,a,a)')
     &   lrpa,',',resfac,',',cmfac,',',cutmod,',',ifsepi,',',
     &   resdmp,',',iforth,',',ecut,',',
     &  'lrpa,resfac,cmfac,cutmod,ifsepi,resdmp,iforth,ecut'
      do iread=1,999
      read(5,'(a)',end=99) a
      write(6,'(a)') a
      enddo
   99 continue
      stop
      end

