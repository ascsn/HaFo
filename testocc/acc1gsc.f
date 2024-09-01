
c-----acc1GSc----------------------------------------------------
c 
      program accgsc
c
c     accumulates GSC results from 'for016' and adds it up for
c     a given nucleus. This version imots LRPA=1. 
c
      character*12 title
      dimension dfull(1:4,0:11),dlocal(1:4,0:11),dum(6),acc(1:4,3) 
c
c------------------------------------------------------------------
c
      open(16,form='formatted',status='old',file='for016')
      nprot = 0
      nneut = 0
      do iouter=1,9999
        read(16,'(a12)',err=999,end=998) title
        if(title(1:4).eq.'FRPA') then
          do icase=1,3
          do iacc=1,4
            acc(iacc,icase) = 0.0
          enddo
          enddo
          backspace(16)
          read(16,'(12x,3i4)') lrpa,nprot,nneut
          backspace(16)
          nprots = nprot
          nneuts = nneut 
          write(6,'(3a,2i4/4a)') 
     &     ' Now analyzing data for:',title,
     &     ' nprot,nneut=',nprot,nneut  ,
     &     'lrpa   r: full   local   diff   ',
     &     '       R: full   local   diff   ',
     &     '       s: full   local   diff   ',
     &     '       E: full   local   diff   '
          do ilrpa=0,11
            read(16,'(12x,3i4,10g13.5)') 
     &        LRp2,Nprot,nneut,dum,(dfull(iacc,lrpa),iacc=1,4)
            read(16,'(12x,3i4,10g13.5)',err=999,end=997) 
     &        LRp,nprot,nneut,dum,(dlocal(iacc,lrpa),iacc=1,4)
            if(Lrp.ne.lrpa) stop ' inconsistent data LRPA'
            if(lrpa.ne.1) then 
            do iacc=1,4
              acc(iacc,1) = dfull(iacc,lrpa) + acc(iacc,1) 
              acc(iacc,2) = dlocal(iacc,lrpa) + acc(iacc,2) 
              acc(iacc,3) = dfull(iacc,lrpa)-dlocal(iacc,lrpa)
     &                     + acc(iacc,3) 
            enddo 
            endif 
            write(6,'(i3,2x,3(3f8.4,2x),3f8.2)')
     &       lrpa,(dfull(iacc,lrpa),dlocal(iacc,lrpa),
     &            (dfull(iacc,lrpa)-dlocal(iacc,lrpa)),iacc=1,4)
            read(16,'(12x,3i4)',err=99,end=99) lrpa,nprot,nneut
            backspace(16)
            if(nprot.ne.nprots) goto 99
            if(nneut.ne.nneuts) goto 99
          enddo 
   99     continue
          write(6,'(6x,4(a,2x))')
     &      '------------------------',
     &      '------------------------',
     &      '------------------------',
     &      '------------------------'
          write(6,'(a,3(3f8.4,2x),3f8.2)')
     &    'total',(acc(iacc,1),acc(iacc,2),acc(iacc,3),iacc=1,4)
          
        endif 
      enddo
c
  998 stop
  999 stop
  997 stop ' missing LRPA'
      end
 
