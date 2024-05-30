      parameter (KBIN=1500)
      parameter (KSTATE=4000)
      parameter (WBIN=0.2)      ! width in MeV
      parameter (WIDTH=5.0)
      parameter (EMINBIN=-100.0)
c
      dimension denslevel(KBIN,2),ener(KSTATE,2)
      character*20 aline
      character*12 title
c
c----------------------------------------------------------------------
c
c     spool to first block proton states
c
      nr = 5
      do iread=1,9999
        read(nr,'(a)',err=98,end=98) aline
        if(aline(1:4).eq.'  2j') goto 999
        if(aline(1:4).eq.' * *') then
          read(nr,'(a)',err=98,end=98) title
        endif
      enddo
 999  continue
c
c     read first block of proton states
c
      ind=0
cold 9    read(nr,101)j2,l,energ
      do iread=1,KSTATE
       read(nr,*) j2,l,energ
c      write(6,101)j2,l,energ
       ind=ind+1
       ener(ind,1)=energ
       if (j2.lt.0) goto 19
      enddo
      stop ' Nr. of states exhausted. enhance KSTATE'
 19   continue
      ilelev=ind
c
c     spool to second block proton states
c
      do iread=1,9999
        read(nr,'(a)',err=98,end=98) aline
        if(aline(1:4).eq.' * *') goto 997
      enddo
 997  continue
      do iread=1,9999
        read(nr,'(a)',err=98,end=98) aline
        if(aline(1:4).eq.'  2j') goto 998
      enddo
 998  continue
c
c     read second block of proton states
c
      ind=0
cold 9    read(nr,101)j2,l,energ
      do iread=1,KSTATE
       read(nr,*) j2,l,energ
c      write(6,101)j2,l,energ
       ind=ind+1
       ener(ind,2)=energ
       if (j2.lt.0) goto 29
      enddo
      stop ' Nr. of states exhausted. enhance KSTATE'
 29   continue
      if(ilelev.ne.ind) stop ' nr. of free states does not match'
      write(*,*)'ind ',ind

c
c     compose level densities
c
      do j=1,2
      do n=1,KBIN
        denslevel(n,j) = 0.0
      enddo      
      enddo      

      do n=1,KBIN
        enbin = n*WBIN+EMINBIN
        do j=1,ilelev       
          denslevel(n,1) =  denslevel(n,1) 
     &                  + exp(-((ener(j,1)-enbin)/WIDTH)**2)
          denslevel(n,2) =  denslevel(n,2) 
     &                  + exp(-((ener(j,2)-enbin)/WIDTH)**2)
        enddo
        write(6,'(2x,f7.2,3g13.5)')
     &    enbin,denslevel(n,1),denslevel(n,2),
     &           denslevel(n,1)-denslevel(n,2)
      enddo
c
 98   continue
c
c
      stop
      end
