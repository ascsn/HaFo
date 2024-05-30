c-----levels-------------------------------------------------------------
c
      program levels
c
c     to extract single particle levels from file 'for027' 
c     and to write a selected window to an output file
c
      character*70 atext
      character*3  typ(16)
      character*8  force
      dimension spenrg(255,16)
c
      data typ/'100','102','104','106','108','110','112','114','116',
     &         '118','120','122','124','126','128','130'/
c
c------------------------------------------------------------------------
c
c     collect all data
c
      open(11,file='for027')
c      do nprotl=1,16
      do iouter=1,29
      do nprotl=1,5
c        open(11,file='for027.Z'//typ(nprotl))
        nprot = 0
        do i1=1,999
        do i2=1,999
          read(11,'(a)',end=99) atext
c          write(6,'(a)') atext
cold          if(atext(1:3).eq.'skp') then
          if(atext(1:4).eq.' * *') then
            read(11,'(a)') force
c            write(6,'(a)') atext
            read(11,'(a)') atext
c            write(6,'(a)') atext
            read(11,'(a17,i5)') atext,nprot
c            write(6,'(a,i5)') ' nprot=',nprot
          elseif(nprot.ne.0 .and. atext(10:14).eq.'spect') then
            read(11,'(1x)') 
            do isp=1,999
              read(11,*,err=19,end=19)  j2,l,energy,weight,rms,multip
c              write(6,*)  j2,l,energy,weight,rms,multip
              ispmax = isp
              spenrg(isp,nprotl) = energy
              multip = multip/2
              do i=2,multip
                read(11,'(1x)') 
              enddo
            enddo
 19         continue
            goto 18
          endif
        enddo
        enddo
 18     continue
c        close(11)
      enddo
c
c     print all s.p. energies sorted in columns for plotting
c
      write(6,'(2a)') ' force:',force
c      do ipriln=0,50,10
      do ipriln=0,20,10
        write(6,'(a,i3,a,i3)') 
     &   ' proton spectra for states',ipriln+1,' until ',ipriln+10
c        do nprotl=1,16
        do nprotl=1,5
c          nprot = 98+2*nprotl
          nprot = 118+2*nprotl
          write(6,'(1x,i4,10f7.2)') 
     &       nprot,(spenrg(isp,nprotl),isp=ipriln+1,ipriln+10)
        enddo
        write(6,'(1x)')
      enddo
      enddo
c
 99   stop
      end





