c-----tabul -------------------------------------------------------------------
c
      program tabul
c
      parameter (KTABMX=50,KLMAX=11)
      parameter (KSUBT=KTABMX*(KLMAX+1))
      parameter (KTOTAL=14*KSUBT)
c
c     evaluates the tables of correlation results
c
      character*132 a
      character*40 tabnam(KTABMX),key,case(2)
      dimension    lmax(KTABMX)
      dimension    datfld(7,0:KLMAX,2,KTABMX)
c
      data case/' 1. case: full correlations              ',
     &          ' 2. case: perturbative diagram subtracted'/
      data datfld/KTOTAL*0.0/
c
      open(unit=2,status='UNKNOWN',form='FORMATTED',file='TABUL.OUT')
      open(unit=3,status='SCRATCH',form='FORMATTED')
c
      itabmx = 0
      do 1000 iread=1,9999
      read(5,'(a)',err=910,end=920) a
      write(2,'(a)') a
c                                        branch to # commands
      if(a(1:1).eq.'#') then
c                                        ':' signals evaluation and print
       if(a(2:2).eq.':') then
        i1     = index(a,'"')
        i2     = index(a(i1+1:132),'"')+i1
        if(i2.eq.i1) stop ' no matching " found in'//a
c                                        find table entry for 1. keyword
        i1     = i1+1
        i2     = i2-1
        ileng  = i2-i1+1
        do 10 i=1,itabmx
        key    = tabnam(i)
        if(key(1:ileng).eq.a(i1:i2)) then
         itab1 = i
         goto 19
        endif
   10   continue
        stop ' keyword '//key//'not found'
   19   continue
        inext  = i2+2
        if(a(inext:inext).eq.'-') then
c                                        branch to difference of tables
         i3     = index(a(inext+1:132),'"')+inext
         if(i3.eq.inext) stop ' missing " after - '
         i4     = index(a(i3+1:132),'"')+i3
         if(i4.eq.i3) stop ' no matching " found in 2. scan'
c                                        find table entry for 2. keyword
         i3     = 1+i3
         i4     = i4-1
         ileng  = i4-i3+1
         do 20 i=1,itabmx
         key    = tabnam(i)
         if(key(1:ileng).eq.a(i3:i4)) then
          itab2 = i
         goto 29
         endif
   20    continue
         stop ' keyword '//key//'not found'
   29    continue
c                                        produce table
         write(2,'(/a)')
     &     ' Difference:'//a(i1:i2)//'-'//a(i3:i4)
         do 30 icase=1,2
         write(2,'(a/a/a)')
     &     case(icase),
     &     '     L      r         R      surface      E  ',
     &     '----------------------------------------------------------'
         srms   = 0.0
         srdiff = 0.0
         ssurf  = 0.0
         senerg = 0.0
         do 40 l=0,lmax(itab1)
         rms    = datfld(1,l,icase,itab1)+datfld(4,l,icase,itab1)
     &           -datfld(1,l,icase,itab2)-datfld(4,l,icase,itab2)
         rdiff  = datfld(2,l,icase,itab1)+datfld(5,l,icase,itab1)
     &           -datfld(2,l,icase,itab2)-datfld(5,l,icase,itab2)
         surf   = datfld(3,l,icase,itab1)+datfld(6,l,icase,itab1)
     &           -datfld(3,l,icase,itab2)-datfld(6,l,icase,itab2)
         energ  = datfld(7,l,icase,itab1)
     &           -datfld(7,l,icase,itab2)
         write(2,'(i6,3f10.4,f10.2)')
     &     l,rms,rdiff,surf,energ
         srms   = rms+srms
         srdiff = rdiff+srdiff
         ssurf  = surf+ssurf
         senerg = energ+senerg
   40    continue
         write(2,'(a)')
     &     '----------------------------------------------------------'
         write(2,'(a,3f10.4,f10.2)')
     &     ' total',srms,srdiff,ssurf,senerg
         write(2,'(a/)')
     &     '----------------------------------------------------------'
   30    continue
        else
c                                        branch to direct table
c                                        produce table
         write(2,'(/a)')
     &     ' Table of '//a(i1:i2)
         do 50 icase=1,2
         write(2,'(a/a/a)')
     &     case(icase),
     &     '     L      r         R      surface      E  ',
     &     '----------------------------------------------------------'
         srms   = 0.0
         srdiff = 0.0
         ssurf  = 0.0
         senerg = 0.0
         do 60 l=0,lmax(itab1)
         rms    = datfld(1,l,icase,itab1)+datfld(4,l,icase,itab1)
         rdiff  = datfld(2,l,icase,itab1)+datfld(5,l,icase,itab1)
         surf   = datfld(3,l,icase,itab1)+datfld(6,l,icase,itab1)
         energ  = datfld(7,l,icase,itab1)
         write(2,'(i6,3f10.4,f10.2)')
     &     l,rms,rdiff,surf,energ
         srms   = rms+srms
         srdiff = rdiff+srdiff
         ssurf  = surf+ssurf
         senerg = energ+senerg
   60    continue
         write(2,'(a)')
     &     '----------------------------------------------------------'
         write(2,'(a,3f10.4,f10.2)')
     &     ' total',srms,srdiff,ssurf,senerg
         write(2,'(a/)')
     &     '----------------------------------------------------------'
   50    continue
        endif
       else
c                                        no ':' in a(2:2) means new
c                                        table entry
        ileng = length(a)-1
        if(ileng.gt.40) stop 'too long keyword'
        itabmx = 1+itabmx
        if(itabmx.gt.KTABMX) stop ' table space KTABMX exhausted'
        key(1:ileng) = a(2:ileng+1)
        tabnam(itabmx) = key
        do 100 iscanl=1,999
        do 110 ir2=1,99
        read(5,'(a)',err=930,end=940) a
        write(2,'(a)') a
        if(a(1:2).eq.'L=') goto 119
        if(a(1:4).eq.'#end') goto 199
  110   continue
        stop ' count exhausted while scanning "L="'
  119   continue
        ileng  = length(a)
        write(3,'(a)') a(3:ileng)
        rewind 3
        read(3,*) L
        rewind 3
        if(L.gt.KLMAX) stop ' too large L encountered in line'//a
        do 120 icase=1,2
        read(5,'(a)',err=950,end=950) a
        write(2,'(a)') a
        i1    = index(a,':')+1
        i2    = length(a)
        write(3,'(a)') a(i1:i2)
        rewind 3
        read(3,*) (datfld(i,L,icase,itabmx),i=1,7)
        rewind 3
  120   continue
        do 130 i=1,7
        if(i.le.3 .or. i.eq.7) then
         datfld(i,L,2,itabmx) = datfld(i,L,1,itabmx)
     &                         -100.0*datfld(i,L,2,itabmx)
        else
         datfld(i,L,2,itabmx) = datfld(i,L,1,itabmx)
        endif
  130   continue
  100   continue
        stop ' no #end mark found'
  199   continue
        lmax(itabmx) = L
       endif
      endif
 1000 continue
      stop ' file count exhausted'
c
  910 stop ' error in reading input'
c
  920 stop ' normal EoF'
c
  930 stop ' error while scanning "L="'
c
  940 stop ' EoF while scanning "L="'
c
  950 stop ' Error or EoF while scanning data for DATFLD'
c
      end
C-----LENGTH-------------------------------------------------------------------
C
      FUNCTION LENGTH(A)
      PARAMETER (KMAXLN=132)
      CHARACTER*(KMAXLN) A
C
C     RETURNS THE LENGTH OF 'A', I.E. THE LAST NON-BLANK PLACE
C     MINIMUM VALUE IS '1'
C
      L = KMAXLN
      DO 10 I=KMAXLN,2,-1
      IF(A(I:I) .NE. ' ') GOTO 19
      L = L-1
   10 CONTINUE
   19 LENGTH = L
      RETURN
      END
