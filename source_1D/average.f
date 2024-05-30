c-----average----------------------------------------------------------
c
      program average
c
c     extract spectral information from file 'for014'
c     and compute peak-weighted average
c
      parameter (kdim=99)
      character*80 a
      dimension  e(kdim),be1p(kdim),ebe1p(kdim)
      dimension  be10(kdim),ebe10(kdim),be11(kdim),ebe11(kdim)
      data       icount/0/
c
c----------------------------------------------------------------------
c
      open(unit=19,status='scratch')
      do 1000 iread=1,9999
      do 1000 iread2=1,999
        read(5,'(a)',err=99,end=99) a
        if(a(4:9).eq.'title:') then
          if(icount.ne.0) then
            aver = evalue(e,ebe11,icount)
            write(6,'(a,1pg12.4)')
     &       ' average energy=',aver
            write(6,'(a)') ' '
          endif
          write(6,'(a)') a(10:len(a))
          do 20 i=1,kdim
            e(i)     = 0.0
            be1p(i)   = 0.0
            ebe1p(i)  = 0.0
            be10(i)   = 0.0
            ebe10(i)  = 0.0
            be11(i)   = 0.0
            ebe11(i)  = 0.0
   20     continue
          icount = 0
        elseif(a(1:4).eq.' nr ') then
          icount = 0
          do 10 i=1,99
            read(5,*,err=19,end=99)
     &        nr,e(i),be1p(i),ebe1p(i),be10(i),ebe10(i),be11(i),ebe11(i)
            icount = 1+icount
   10     continue
   19     continue
        endif
 1000 continue
   99 continue
      if(icount.ne.0) then
            aver = evalue(e,ebe11,icount)
            write(6,'(a,1pg12.4)')
     &       ' average energy=',aver
      endif
      stop
      end
c-----evalue----------------------------------------------------------
c
      function evalue(e,be,ndim)
c
c     iterate weighted average
c
      parameter (kdim=99)
      dimension  e(kdim),be(kdim)
      data       width/5.0/
c
c----------------------------------------------------------------------
c
c     first guess
c
      sum0 = 0.0
      sum1 = 0.0
      do i=1,ndim
       sum0  = be(i)      + sum0
       sum1  = be(i)*e(i) + sum1
      enddo      
      avere = sum1/sum0
c
c     now iterate
c
      widinv = 0.5/width**2
c      avere = 10.0
      do iter=1,10
        averol = avere
        sum0 = 0.0
        sum1 = 0.0
        do i=1,ndim
         wfac  = exp(-widinv*(e(i)-avere)**2)
         sum0  = wfac*be(i)      + sum0
         sum1  = wfac*be(i)*e(i) + sum1
        enddo      
        avere = sum1/sum0
c        write(6,'(a,i5,1pg13.4)')   ' iter,avere=',iter,avere
        if(abs(avere-averol).lt.01) goto 99
      enddo
 99   continue
      evalue = avere
      return
      end
