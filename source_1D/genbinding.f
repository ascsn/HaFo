      character*2 a
      logical tjob
      data tjob/.false./
c
      if(tjob) open(19,status='scratch')
      do iz=182,200,2
c      do iz=140,200,2
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
        ins = 130+(iz-110)*2
        ine = 240-130+ins
c        ins = 48+(iz-50)
c        ine = 100+(iz-50)
        do in=ins,ine,2
          if(tjob) then
            write(2,'(a,2i4)') '~/hafo/binding_ski3.bat',iz,in
          else
            write(6,'(a,2i4)') '~/hafo/binding_ski3.bat',iz,in
          endif
        enddo
      enddo
      stop
      end


