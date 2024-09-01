c      
c     reads the results for 'binding energy' from 'for002'
c     and deduces S_2n and S_p from it if 'for002' is properly arranged
c
      dimension energ(3,18)
      character*20 atext
c
      do ir1=1,29999
        read(5,'(a)',err=99,end=99) atext
c        write(6,'(a)') atext
        if(atext(1:10).eq.'   binding') then
c          do izl=98,132,2
c          do inl=182,186,2
          do izl=82,82
          do inl=118,134,2
           read(5,*,err=98,end=98) ia,iz,dummy,energy
           write(6,'(2i4,2f10.1)') ia,iz,dummy,energy
c           izs = (iz-96)/2
           izs = (iz-80)/2
           in  = ia-iz
c           ins = (in-180)/2
           ins = (in-116)/2
c           energ(ins,izs) = energy
           energ(izs,ins) = energy
          enddo
          enddo
          write(6,'(a)') 
     &     '  iz  in  energy  S_2n   S_2p  delta_2n delta_2p '
c          do izs=2,17
c            iz = 2*iz+96
c            in = 184
          do izs=2,9
            iz = 82
            in = 116+2*izs
            write(6,'(1x,2i4,f8.1,2f6.1,2f6.1)')
     &       iz,in,energ(1,izs),
     &       energ(1,izs+1)-energ(1,izs),
     &       energ(1,izs+1)-2.0*energ(1,izs)+energ(1,izs-1)
c            write(6,'(1x,2i4,f8.1,2f6.1,2f6.1)')
c     &       iz,in,energ(2,izs),
c     &       energ(2,izs)-energ(2,izs),
c     &       energ(2,izs+1)-energ(2,izs),
c     &       energ(1,izs)-2.0*energ(2,izs)+energ(3,izs),
c     &       energ(2,izs+1)-2.0*energ(2,izs)+energ(2,izs-1)
          enddo
        endif
      enddo
c
      

 99   stop ' error or EoF in main reading '
 98   stop ' error or EoF in data reading'
      end
