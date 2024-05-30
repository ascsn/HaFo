c-----extract----------------------------------------------------------
c
      program extrac
c
c     extract spectral information from file 'for014'
c
      parameter (kdim=99)
      character*80 a
      dimension  e(kdim),be1p(kdim),ebe1p(kdim),width(kdim)
      dimension  be10(kdim),ebe10(kdim),be11(kdim),ebe11(kdim)
      data       icount/0/
c
c----------------------------------------------------------------------
c
      open(unit=19,status='scratch')
      do 1000 iread=1,9999
        read(5,'(a)',err=99,end=99) a
        if(a(4:9).eq.'TITLE:') then
          if(icount.ne.0) then
            do 30 i=1,icount
            write(6,'(i5,8(1pg12.4))')
     &       i,e(i),be1p(i),ebe1p(i),be10(i),ebe10(i),be11(i),ebe11(i),
     &         width(i)
   30       continue
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
            width(i) = 0.0
   20     continue
          icount = 0
        elseif(a(1:4).eq.' NR ') then
          write(6,'(3x,a,a)') a(1:80),'        width'
          do 10 i=1,99
            read(5,*,err=19,end=99)
     &        nr,e(i),be1p(i),ebe1p(i),be10(i),ebe10(i),be11(i),ebe11(i)
   10     continue
   19     continue
        elseif(a(4:9).eq.'Del_e3') then
               icount = 1+icount
               rewind 19
               write(19,'(a)') a(12:23)
               rewind 19
               read(19,*,err=99) width(icount)
        endif
 1000 continue
   99 continue
      if(icount.ne.0) then
        do 40 i=1,icount
            write(6,'(i5,8(1pg12.4))')
     &       i,e(i),be1p(i),ebe1p(i),be10(i),ebe10(i),be11(i),ebe11(i),
     &         width(i)
   40   continue
      endif
      stop
      end
