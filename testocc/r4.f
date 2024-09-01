      do ir=1,999
      read(5,*,err=99,end=99) ir,r,x1,x2,x3,rho
      write(6,'(1x,i5,3g12.4)') ir,r,rho,rho*r**4
      enddo
   99 stop
      end
