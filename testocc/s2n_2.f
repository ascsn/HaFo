c      
c     reads the results for 'binding energy' from 'fort.23'
c     and deduces S_2n and S_p from it
c
      dimension energ(3,18)
      character*20 atext
c
      do ir1=1,99
        do izl=98,132,2
        do inl=182,186,2
         read(5,*,err=98,end=98) iz,in,energy
c         write(6,'(2i4,2f10.1)') ia,iz,dummy,energy
         izs = (iz-96)/2
         ins = (in-180)/2
         energ(ins,izs) = energy
        enddo
        enddo
        write(6,'(a)') 
     &     '  iz  in  energy  S_2n   S_2p  delta_2n delta_2p '
        do izs=2,17
          iz = 2*izs+96
          in = 184
          write(6,'(1x,2i4,f8.1,2f6.1,2f6.1)')
     &       iz,in,energ(2,izs),
     &       energ(3,izs)-energ(2,izs),
     &       energ(2,izs+1)-energ(2,izs),
     &       energ(1,izs)-2.0*energ(2,izs)+energ(3,izs),
     &       energ(2,izs+1)-2.0*energ(2,izs)+energ(2,izs-1)
        enddo
      enddo
c
 98   stop ' error or EoF in data reading'
      end
