      do i=1,9999
        read(5,*,end=99) nz,nn,Elow
        read(5,*,end=99) nz,nn,Emid
        read(5,*,end=99) nz2,nn2,Eup
        write(6,'(1x,2i4,3(1pg13.5))') 
     &    nz,nn,Elow-2*Emid+Eup,Elow-Emid,Emid-Eup
      enddo
 99   stop
      end
