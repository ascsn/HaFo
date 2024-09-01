      character*80 a
      open(11,form='formatted',status='old',file='for011')
      write(6,'(a/a)')
     &   ' neutron radii:',
     &   '  Z   N   A   rms,rdms,surf,smom,halo '
      do iread=1,29999
        read(11,'(a)',err=999,end=999) a
        if(a(2:4).eq.'Z,N') then
          backspace 11
          read(11,'(19x,2i5)') IZ,IN
          backspace 11
          backspace 11
          read(11,'(14x,5g13.5)') rms,rdms,surf,smom,halo
          read(11,'(1x)')
          write(6,'(1x,2i5,5(1pg13.5))')
     &     Iz,IN,rms,rdms,surf,smom,halo
        endif
      enddo
  999 continue
      stop
      end

