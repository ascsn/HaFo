      PROGRAM SPLIT
      CHARACTER*142 A
      do 10 i1=1,9
        if(i1.eq.1) then
          open(11,form='formatted',status='new',file='x1')
        elseif(i1.eq.2) then
          open(11,form='formatted',status='new',file='x2')
        elseif(i1.eq.3) then
          open(11,form='formatted',status='new',file='x3')
        elseif(i1.eq.4) then
          open(11,form='formatted',status='new',file='x4')
        elseif(i1.eq.5) then
          open(11,form='formatted',status='new',file='x5')
        elseif(i1.eq.6) then
          open(11,form='formatted',status='new',file='x6')
        elseif(i1.eq.7) then
          open(11,form='formatted',status='new',file='x7')
        elseif(i1.eq.8) then
          open(11,form='formatted',status='new',file='x8')
        elseif(i1.eq.9) then
          open(11,form='formatted',status='new',file='x9')
        endif
        do 20 x2=1,5000
          read(5,'(a)',err=900,end=900) a
          lena = length(a)
          if(lena.gt.0) then
           write(11,'(a)') a(1:lena)
          else
           write(11,'(1x)')
          endif
   20   continue
        close(11)
   10 continue
  900 STOP
      END
c-----length-----------------------------------------------------------
c
      function length(text)
c
      character*142 text
c
      do 10 l=142,1,-1
      ileng = l
      if(text(l:l).ne.' ') goto 19
   10 continue
      ileng = 0
   19 continue
      length = ileng
c
      return
      end

