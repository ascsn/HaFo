      character*2 a
      logical tjob
      data tjob/.true./
c
      if(tjob) open(19,status='scratch')
      do iz=40,60,2
c      do iz=142,200,2
c        ins = 202+(iz-110)*2
c        ine = 240-202+ins
        if(tjob) then
          write(19,'(i2)') iz
          rewind 19
          read(19,'(a2)') a
          rewind 19
          open(2,file='sh'//a//'.job')
          write(2,'(a)') 'cd /tmp/mpt218'
        endif
c        ins = 130+(iz-110)*2
c        ine = 240-130+ins
        ins = 48+(iz-50)
        ine = 100+(iz-50)
        do in=ins,ine,2
          if(tjob) then
            write(2,'(a,2i4)') '~/hafo/levden_skp.bat',iz,in
          else
            write(6,'(a,2i4)') '~/hafo/levden_skp.bat',iz,in
          endif
        enddo
        if(tjob) then
          write(2,'(a)')
     &     'rm -f for027.*',
     &     'rm -f *.out',
     &     'cd /home/mpt218/hafo/res',
     &     'cat /tmp/mpt218/shg*.res >> sn_skp.res',
     &     'gzip -f /tmp/mpt218/shg*.res'
          close(2)
        endif
      enddo
      stop
      end


