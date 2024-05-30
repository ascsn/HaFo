c      data  enucm/15.902/    ! SLy6  
c      data  enucm/15.96/    ! SkI3  
c      data  enucm/15.752/    ! skms
c      data enucm/15.925/     !   SkI4   
c      data enucm/15.835/     !  SkO4   
      data enucm/15.750/     !   SkO4'  

      do i=1,32
        read(5,*) nz,esh1,esh2,etotal,epair,zpe
        eldm = (etotal-esh1-esh2-epair)
        na = 2*nz
        a  = na
        a13 = a**(-0.3333333)
        etota = etotal/a
        eldma = eldm/a
        eldza = (eldm+zpe)/a
        write(6,'(1x,i4,f8.5,7f8.3)')
     &    na,a13,etota,eldma,eldza,(etota-epair/a+enucm)/a13,
     &     (eldma+enucm)/a13,(eldza+enucm)/a13
      enddo
      stop
      end
