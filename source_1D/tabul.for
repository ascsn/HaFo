c-----tabul -------------------------------------------------------------------
c
      program tabul
c
      parameter (KTABMX=50,KLMAX=11,KMAXLN=80,KOBSER=4,KCASE=5)
      parameter (KSUBT=KTABMX*(KLMAX+1))
      parameter (KTOTAL=14*KSUBT)
c
c     evaluates the tables of correlation results
c
      character*(KMAXLN) a
      character*40 tabnam(KTABMX),key
      character*30 case(KCASE)
      dimension    lmax(KTABMX)
      dimension    datfld(KOBSER,0:KLMAX,KCASE,KTABMX),
     &             sumfld(KOBSER,KCASE)
      data case/' UNCORRECTED',' PP-CORRECTED',
     &          ' 2.ORDER CORRECTED',' PP-CORRECTED FOR SCAT.PART',
     &          ' SCAT. PART 2.ORDER CORRECTED'/
c
      do 310 itab=1,KTABMX
      do 310 icase=1,KCASE
      do 310 l=1,KLMAX
      do 310 iobser=1,KOBSER
        datfld(iobser,l,icase,itab) = 0.0
  310 continue
c
      open(unit=5,status='UNKNOWN',form='FORMATTED',file='TABUL.INP')
      open(unit=2,status='UNKNOWN',form='FORMATTED',file='TABUL.OUT')
      open(unit=3,status='SCRATCH',form='FORMATTED')
c
      itabmx = 0
      do 1000 iread=1,9999
      read(5,'(a)',err=910,end=920) a
      write(2,'(a)') a
      write(6,'(a)') a
c                                        branch to # commands
      if(a(1:1).eq.'#') then
c                                        ':' signals evaluation and print
       if(a(2:2).eq.':') then
        i1     = index(a,'"')
        i2     = index(a(i1+1:KMAXLN),'"')+i1
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
         i3     = index(a(inext+1:KMAXLN),'"')+inext
         if(i3.eq.inext) stop ' missing " after - '
         i4     = index(a(i3+1:KMAXLN),'"')+i3
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
         do 30 icase=1,KCASE
         write(2,'(a/a/a)')
     &     case(icase),
     &     '     L      r         R      surface      E  ',
     &     '----------------------------------------------------------'
         srms   = 0.0
         srdiff = 0.0
         ssurf  = 0.0
         senerg = 0.0
         do 40 l=0,lmax(itab1)
         rms    = datfld(1,l,icase,itab1)
     &           -datfld(1,l,icase,itab2)
         rdiff  = datfld(2,l,icase,itab1)
     &           -datfld(2,l,icase,itab2)
         surf   = datfld(3,l,icase,itab1)
     &           -datfld(3,l,icase,itab2)
         energ  = datfld(4,l,icase,itab1)
     &           -datfld(4,l,icase,itab2)
         write(2,'(i6,3f10.4,f10.2)')
     &     l,rms,rdiff,surf,energ
         srms   = rms+srms
         srdiff = rdiff+srdiff
         ssurf  = surf+ssurf
         senerg = energ+senerg
   40    continue
         sumfld(1,icase) = srms
         sumfld(2,icase) = srdiff
         sumfld(3,icase) = ssurf
         sumfld(4,icase) = senerg
         write(2,'(a)')
     &     '----------------------------------------------------------'
         write(2,'(a,3f10.4,f10.2)')
     &     ' total',srms,srdiff,ssurf,senerg
         write(2,'(a/)')
     &     '----------------------------------------------------------'
   30    continue
         write(2,'(/a,2(/30x,a)/)')
     &    ' SUMMARY OF TOTAL CONTRIBUTIONS:',
     &    '     r         R       surface       E  ',
     &    '------------------------------------------'
         do 31 icase=1,KCASE
         write(2,'(a30,3f10.4,f10.2)')
     &     case(icase),(sumfld(iobs,icase),iobs=1,KOBSER)
   31    continue
        else
c                                        branch to direct table
c                                        produce table
         write(2,'(/a)')
     &     ' Table of '//a(i1:i2)
         do 50 icase=1,KCASE
         write(2,'(a/a/a)')
     &     case(icase),
     &     '     L      r         R      surface      E  ',
     &     '----------------------------------------------------------'
         srms   = 0.0
         srdiff = 0.0
         ssurf  = 0.0
         senerg = 0.0
         do 60 l=0,lmax(itab1)
         rms    = datfld(1,l,icase,itab1)
         rdiff  = datfld(2,l,icase,itab1)
         surf   = datfld(3,l,icase,itab1)
         energ  = datfld(4,l,icase,itab1)
         write(2,'(i6,3f10.4,f10.2)')
     &     l,rms,rdiff,surf,energ
         srms   = rms+srms
         srdiff = rdiff+srdiff
         ssurf  = surf+ssurf
         senerg = energ+senerg
   60    continue
         sumfld(1,icase) = srms
         sumfld(2,icase) = srdiff
         sumfld(3,icase) = ssurf
         sumfld(4,icase) = senerg
         write(2,'(a)')
     &     '----------------------------------------------------------'
         write(2,'(a,3f10.4,f10.2)')
     &     ' total',srms,srdiff,ssurf,senerg
         write(2,'(a/)')
     &     '----------------------------------------------------------'
   50    continue
         write(2,'(/a,2(/30x,a)/)')
     &    ' SUMMARY OF TOTAL CONTRIBUTIONS:',
     &    '     r         R       surface       E  ',
     &    '-----------------------------------------'
         do 51 icase=1,KCASE
         write(2,'(a30,3f10.4,f10.2)')
     &     case(icase),(sumfld(iobs,icase),iobs=1,KOBSER)
   51    continue
        endif
       else
c
c                                        no ':' in a(2:2) means new
c                                        table entry
c
        ileng = length(a)-1
        if(ileng.gt.40) stop 'too long keyword'
        itabmx = 1+itabmx
        if(itabmx.gt.KTABMX) stop ' table space KTABMX exhausted'
        key(1:ileng) = a(2:ileng+1)
        tabnam(itabmx) = key
c                                        scan LRPA and G.S.C.
        lmaxa = -99
        do 100 ir2=1,999
        read(5,'(a)',err=930,end=940) a
        write(2,'(a)') a
        if(index(a,'LRPA=').ne.0) then
         inl = index(a,'LRPA=')+5
         write(3,'(a)') a(inl:inl+2)
         rewind 3
         read(3,*) l
         rewind 3
         lmaxa = max(l,lmaxa)
         if(l.gt.KLMAX) stop ' too large LRPA'
         do 110 icase=1,KCASE
         read(5,'(a)',err=930,end=940) a
         write(2,'(a)') a
         incol  = index(a,':')
         incomm = index(a,',')
         lencol = incol-incomm-1
ctest         write(6,'(2a/a,4i5)')
ctest     &    ' a=',a,' l,incol,incomm,lencol=',l,incol,incomm,lencol
         if(a(1:6).ne.' G.S.C' .or.
     &       a(incomm+1:incol-1).ne.case(icase)(1:lencol) )  then
          write(6,'(a/a)')
     &     ' inconsistent data in:',a
          stop
         endif
         ileng = length(a)
         write(3,'(a)') a(39:ileng)
         rewind 3
         read(3,*) (datfld(iobs,l,icase,itabmx),iobs=1,KOBSER)
         rewind 3
  110    continue
        elseif(a(1:4).eq.'#end') then
         goto 119
        endif
  100   continue
        stop ' no #end mark found'
  119   continue
        lmax(itabmx) = lmaxa
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
      PARAMETER (KMAXLN=80)
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
