      dimension e(150,200)
      parameter (ealpha=-28.3)
c
      do i=1,150
      do j=1,200
        e(i,j) = 0.0
      enddo
      enddo
c
      open(1,file='energ.ski3')
c
      do in=1,9999
        read(1,*,err=99,end=99) nz,nn,energ
        e(nz,nn) = energ
      enddo
 99   continue
c
      do nz=104,130,2
      do nn=142,190,2
        if(e(nz,nn).ne.0.0 .and. e(nz-2,nn-2).ne.0.0) then
          qalpha = e(nz,nn) - e(nz-2,nn-2) - ealpha
          write(6,'(2i5,f8.1)') nz-1,nn-1,qalpha
          write(6,'(2i5,f8.1)') nz-1,nn+1,qalpha
        endif
      enddo
      write(6,'(2x)')
      do nn=142,190,2
        if(e(nz,nn).ne.0.0 .and. e(nz-2,nn-2).ne.0.0) then
          qalpha = e(nz,nn) - e(nz-2,nn-2) - ealpha
          write(6,'(2i5,f8.1)') nz+1,nn-1,qalpha
          write(6,'(2i5,f8.1)') nz+1,nn+1,qalpha
        endif
      enddo
      write(6,'(2x)')
      enddo
c
      stop
      end
