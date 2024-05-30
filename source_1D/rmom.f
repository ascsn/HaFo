c    reads density and prints r^n*density for plotting
      dimension x(0:8),rr(0:8)
      data x/9*0.0/
      do ir=1,9999
      read(5,*,err=99,end=99) iact,r,rho
      do im=0,8
      rr(im) = rho*r**im
      x(im) = max(x(im),rr(im))
      enddo
      write(6,'(1x,i4,f7.2,10g12.4)')
     &  iact,r,(rr(im),im=0,8)
      enddo
   99 continue
      write(6,'(a,10g12.4)') 'maxima=',x
      stop
      end

