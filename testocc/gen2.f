c-----geninp------------------------------------------------------------
c
      program geninp
c
c     to generate 'for005' file for computing shell gaps.
c     requires basic input (forces etc) from 'for005.gen'.
c
      character*100 a
      dimension temprt(5)
c
      data nzfix/28/              ! fixed proton number
      data nnmin,nnmax/26,52/    ! limits neutron number
      data temprt/0.01,0.1,0.2,0.4,0.6/
c
c-----------------------------------------------------------------------
c
      open(1,file='for005.gen')
c
      if(mod(nnmin,2).ne.0) stop ' proton number must be even'
c
      do itemp=1,5
      do nn=nnmin,nnmax,2
      do nz=nzfix-2,nzfix+2,2
        rewind 1
        do i2=1,14
          read(1,'(a)') a
          write(6,'(a)') a(1:lenstr(a))        
        enddo
        read(1,*) power,acoul,ipair,pairfp,pairfn,ilnrsp,rtemp
        write(6,'(1x,f8.6,f4.1,i2,2(1pg13.5),i2,0pf6.2)')
     &          power,acoul,ipair,pairfp,pairfn,0,temprt(itemp)
        write(6,'(5i4)') nz,nn,0,0,2
        write(6,'(i2)') 0
      enddo
      enddo
      enddo
      write(6,'(a)') ' -1, "to stop"'
c
      stop
      end
c-----lenstr-----------------------------------------------------------------
c
      function lenstr(a)
c
c     length of string 'a' until last nonempty character
c
      parameter (kmax=20000,kfmax=100,krchar=134,kwchar=20)
c
      character*(*) a
      character*(krchar) b
      b = a
c
      do 10 i=1,krchar
      istr   = krchar+1-i
      if(b(istr:istr).ne.' ') goto 19
   10 continue
      lenstr = 0
      return
c
   19 continue
      lenstr = istr
      return
      end


