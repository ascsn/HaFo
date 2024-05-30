c-----lsmoot--------------------------------------------------------------------
c
      program lsmoot
c
      parameter (klinmx=200,kcolmx=12,komgmx=9000)
c
c     transforms disrete spectra with transition strengths into a
c     smooth distrubution using Lorentzian folding with width 1 in units
c     of the input spectra.
c
c     The input file is directly the 'for014' as produced by 'haforpa'.
c     The smoothing parameters are given at compile time in data.
c     Input string is interrupted by the first occuring blank line.
c     Output is produced. Then reading continues at the first nonblank
c     line with the sequence as explained above.
c     The EoF mark stops the program.
c
      character*80 comment,a
      dimension xin(klinmx),yin(klinmx),dlinin(kcolmx),widin(klinmx)
      dimension stren(komgmx),omega(komgmx)
c
      logical tlower
      data tlower/.true./
c
      data sqpi/1.772454/
      data xmin,xmax,delx,width,isospi/0,20,0.05,0.5,2/
c
c-------------------------------------------------------------------------------
      do i=1,80
        comment(i:i) = ' '
      enddo
c
c     open scratch file
c
      open(unit=19,status='scratch')
c
      do iread=1,9999
       read(5,'(a)',err=999,end=999) a
       if(a(1:3).eq.' nr') goto 996
      enddo
      stop 'file count exhausted'
 996  continue
      maxcol = 8
      nxcol  = 2
      nycol  = 5 + (isospi-1)*2
      nwcol  = 0
      ncount = 0
c
c                                        reading data
c                                        and accumulating b(el)-values in window
       sumbel = 0.0
       sumews = 0.0
       sume2  = 0.0
       maxdat = 0
       do 20 idat=1,klinmx
       read(5,'(a)',err=993,end=999) a
       if(lenstr(a).eq.0 .or. a(1:6).eq.' total') then
        goto 29
       else
        maxdat = 1+maxdat
        write(19,'(a)') a
        rewind 19
        read(19,*) (dlinin(j),j=1,maxcol)
        rewind 19
        xin(maxdat) = dlinin(nxcol)
        yin(maxdat) = dlinin(nycol)
        if(nwcol.gt.0) then
         widin(maxdat) = dlinin(nwcol)
        else
         widin(maxdat) = width
        endif
        if(xin(maxdat).gt.xmin .and. xin(maxdat).lt.xmax) then
         sumbel = yin(maxdat) + sumbel
         sumews = yin(maxdat)*xin(maxdat) + sumews
         sume2  = yin(maxdat)*xin(maxdat)**2 + sume2
        endif
       endif
   20  continue
       read(5,'(a)') a
       if(lenstr(a).ne.0) stop ' too much data'
c
   29  continue
       if(ncount.gt.1) write(6,'(a)') 'NEXT:'
       write(6,'(a,a)')  comment(1:30)
       write(6,'(a/a/)')
     &   'X:E','Y:S(E) arb. units'
       write(6,'(a,a)')  'N:',comment(1:78)
       avener = sumews/sumbel
       avenr2 = sqrt(abs(sume2/sumbel-avener*avener))
       write(6,'(2(A,F6.2),2(a,g12.4)/2(a,g12.4))')
     &   'N: IN WINDOW XMIN=',xmin,' XMAX=',xmax,
     &   ' : B(EL)=',sumbel,'  EWSR=',sumews,
     &   'N:                 AVERAGE ENERGY=',avener,
     &   '  ENERGY-WIDTH=',avenr2
       if(nwcol.eq.0) then
        write(6,'(a,g11.4)')
     &   'N: SCAN WITH FIXED WIDTH=',width
       else
        write(6,'(a)')
     &   'N: SCAN WITH VARIABLE WIDTH.'
       endif
c                                        testprint
c      write(6,'(a,i5)')
c    &     ' maxdat=',maxdat
c      write(6,'(a,10(/1x,5g12.4))')
c    &     ' xin:',(xin(idat),idat=1,maxdat)
c      write(6,'(a,10(/1x,5g12.4))')
c    &     ' yin:',(yin(idat),idat=1,maxdat)
c                                        computing smooth curves
       write(6,'(a,i1,a)') 'H:     X          Y',ncount,',LS'
       numpoi = (xmax-xmin)/delx + 1
       delx   = (xmax-xmin)/(numpoi-1)
       ymax   = zero
       if(numpoi.gt.komgmx) stop ' too small OMEGA grid'
       strmax = 0.0
       istrmx = 0
       do 30 ix=1,numpoi
       x      = xmin + (ix-1)*delx
       y      = 0.0
       do 40 idat=1,maxdat
       xdiff  = (x-xin(idat))
       widac  = widin(idat)
c       widac  = 0.775*widac                 ! make it comparable to Gaussian width
       y      = yin(idat)*widac/(xdiff*xdiff+widac*widac) + y
       if(tlower) then
         xdiff  = (x+xin(idat))
         widac  = widin(idat)
         y      = -yin(idat)*widac/(xdiff*xdiff+widac*widac) + y
       endif
   40  continue
       omega(ix) = x
       stren(ix) = y
       if(y.gt.strmax) then
         strmax = y
         istrmx = ix
       endif
       write(6,'(1x,2g12.4)') x,y
       ymax   = max(y,ymax)
   30  continue
       write(6,'(/a,1pg12.5/)') ' Maximum Y:',ymax
       write(6,'(a,f5.2,a,f5.2,a,a)')
     &  'S: DIST(2.0) XMIN(',xmin,') XMAX(',xmax,') YMIN(0.0)',
     &  ' SCAL(0.8) THICKN(0.7)'
c
c      analyse peak and FWHM
c
       do ix=istrmx,numpoi
         ixup = ix
         if(stren(ix).lt.strmax*0.5) goto 19
       enddo
 19    continue
       do ix=istrmx,1,-1
         ixdw = ix
         if(stren(ix).lt.strmax*0.5) goto 18
       enddo
 18    continue
       write(6,'(a,4(1pg12.4))')
     &  '# peak,FWHM,upper,lower=',
     &    omega(istrmx),omega(ixup)-omega(ixdw),
     &    omega(ixup),omega(ixdw)
c      endif
 1000 continue
      stop ' normal end '
c
  991 stop ' error in one of the heading lines'
c
  992 stop ' error in column specifier line'
c
  993 stop ' error in data line'
c
  999 stop ' EoF on input string'
c
      end
c-----lenstr-----------------------------------------------------------------
c
      function lenstr(a)
c
c     length of string 'a' until last nonempty character
c
      parameter (kmax=200,kfmax=20,krchar=80,kwchar=20)
c
      character*80 a
c
      do 10 i=1,krchar
      istr   = krchar+1-i
      if(a(istr:istr).ne.' ') goto 19
   10 continue
      lenstr = 0
      return
c
   19 continue
      lenstr = istr
      return
      end
