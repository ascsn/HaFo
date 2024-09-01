      character*20 a
      do ir=1,999
        read(5,'(a)',err=999,end=999) a
        if(a(1:8).eq.' mode nr') then
          read(5,'(8x,g12.4,7x,g12.4,8x,g12.4,6x,g12.4)') 
     &                 E0,ep,erpa,VRes
          write(6,'(1x,4(1pg13.4))') e0,ep,erpa,vres 
        endif
      enddo
  999 stop
      end

