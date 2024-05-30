c-----lsmoot-------------------------------------------------------------------
c
      subroutine lsmoot(width)
c
      implicit double precision (a-h,o-z)
c
      parameter (klinmx=900,kcolmx=12)
c
c     transforms disrete spectra with transition strengths into a
c     smooth distrubution using Lorentzian folding with width 1 in units
c     of the input spectra.
c
c     The soubroutine is completely independent of the calling
c     code. Communication is done via the file 'strength_full.res'.
c     The input file contains the discrete spectra in the form:
c      1. line      = comment which is printed in output
c      2. line      = X-window: XMIN,XMAX,DELX for plot
c      3. line      = specification of columns (X, Y, W or N)
c      4. ... lines = spectral point; X-column contains energies,
c                     Y-column contains strengths; N-columns are
c                     omitted
c     Input string is interrupted by the first occuring blank line.
c     Output is produced. Then reading continues at the first nonblank
c     line with the sequence as explained above.
c     The EoF mark stops the program.
c
      character*82 comment
      character*256 a
      dimension xin(klinmx),yin(klinmx),dlinin(kcolmx),widin(klinmx)
c
      logical tlower
      data tlower/.true./
c
!      data width/0.5/
      data sqpi/1.772454/
      data PI/3.1415926D0/
      parameter (widsr=1.0)   ! width for low-energy sum rule
      parameter (cutsr0=4.0)  ! lower end of low-energy sum rule
      parameter (cutsr1=9.5) ! upper end of low-energy sum rule
c
c-------------------------------------------------------------------------------
c
c     open scratch file
c
      open(unit=19,status='scratch')
      open(unit=1,status='old',file='lsmooth.inp')
c
      ncount = 0
c
      do 1000 iread=1,9999
c      read(1,'(a)',err=991,end=999) a(1:50)
      read(1,'(a)',err=999,end=999) a(1:50)
      write(6,'(a)') a
      if(lenstr(a).ne.0) then
       ncount  = 1+ncount
       comment = a
       read(1,*,err=992,end=999) xmin,xmax,delx,width
c       write(6,'(1x,4f8.3,)') xmin,xmax,delx,width
       read(1,'(a)',err=992,end=999) a
c       write(6,'(a)') a
c                                        analysing column structure
       maxcol = 0
       nxcol  = 0
       nycol  = 0
       nwcol  = 0
       do i=1,80
         if(a(i:i) .ne. ' ') then
           maxcol = 1+maxcol
           if(a(i:i) .eq. 'X') then
            nxcol  = maxcol
           elseif(a(i:i) .eq. 'Y') then
            nycol  = maxcol
           elseif(a(i:i) .eq. 'W') then
            nwcol  = maxcol
           endif
         endif
       enddo
       if(nxcol.eq.0) stop ' missing X-column'
       if(nycol.eq.0) stop ' missing Y-column'
c                                        reading data
c                                        and accumulating b(el)-values in window
       sumbel = 0.0
       sumews = 0.0
       sume2  = 0.0
       maxdat = 0
       do idat=1,klinmx
         read(1,'(a)',err=993,end=999) a
         write(6,'(a)') a
         if(lenstr(a).eq.0) then
           goto 29
         else
           write(19,'(a)') a
           write(6,'(a)') a
           rewind 19
           read(19,*) (dlinin(j),j=1,maxcol)
           rewind 19
c           if(idat.ne.2 .or. dlinin(nxcol).gt.6.0) maxdat = 1+maxdat 
           maxdat = 1+maxdat   
           xin(maxdat) = dlinin(nxcol)
           yin(maxdat) = dlinin(nycol)
           if(nwcol.gt.0) then
            widin(maxdat) = dlinin(nwcol)
           else
c            widin(maxdat) = width + max(0.1*(xin(maxdat)-10.0),0.0)
            widin(maxdat) = width 
           endif
           if(xin(maxdat).gt.xmin.and.xin(maxdat).lt.xmax) then
            sumbel = yin(maxdat) + sumbel
            sumews = yin(maxdat)*xin(maxdat) + sumews
            sume2  = yin(maxdat)*xin(maxdat)**2 + sume2
           endif
         endif
       enddo
       read(1,'(a)') a
       if(lenstr(a).ne.0) stop ' too much data'
c
   29  continue
       close(1)
c       open(unit=2,status='unknown',file='strength_full.res')
c       if(ncount.gt.1) write(2,'(a)') 'NEXT:'
       write(2,'(a,a)')  '#'//comment(1:30)
c       write(2,'(a/a/)')
c     &   'X:E','Y:S(E) arb. units'
c       write(2,'(a,a)')  'N:',comment(1:78)
       avener = sumews/sumbel
       avenr2 = sqrt(abs(sume2/sumbel-avener*avener))
       write(2,'(2(A,F6.2),2(a,g12.4)/2(a,g12.4))')
     &   '#: IN WINDOW XMIN=',xmin,' XMAX=',xmax,
     &   ' : B(EL)=',sumbel,'  EWSR=',sumews,
     &   '#:                 AVERAGE ENERGY=',avener,
     &   '  ENERGY-WIDTH=',avenr2
       if(nwcol.eq.0) then
         write(2,'(a,g11.4)')
     &    '#: SCAN WITH FIXED WIDTH=',width
       else
         write(2,'(a)')
     &    '#: SCAN WITH VARIABLE WIDTH.'
       endif
c                                        testprint
c      write(2,'(a,i5)')
c    &     ' maxdat=',maxdat
c      write(2,'(a,10(/1x,5g12.4))')
c    &     ' xin:',(xin(idat),idat=1,maxdat)
c      write(2,'(a,10(/1x,5g12.4))')
c    &     ' yin:',(yin(idat),idat=1,maxdat)
c                                        computing smooth curves
       write(2,'(a)') 'H:     X          Y,L'
       numpoi = (xmax-xmin)/delx + 1
       delx   = (xmax-xmin)/(numpoi-1)
       ymax   = zero
       omax   = zero
       ewsr   = zero            ! accumulator for energy-weighted sum rule
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
       write(2,'(1x,2g12.4)') x,y/PI
       ewsr = y*x + ewsr
       if(y.gt.ymax) then
         ymax = y
         omax = x
       endif
   30  continue
       write(2,'(a,1pg12.5)') '# Maximum Y  :',ymax
       write(2,'(a,1pg12.5)') '# omega there:',omax
       write(2,'(a,1pg12.5)') '# EWSR       :',ewsr*delx/3.1415926
c       write(2,'(a,f5.2,a,f5.2,a,a)')
c     &  'S: DIST(2.0) XMIN(',xmin,') XMAX(',xmax,') YMIN(0.0)',
c     &  ' SCAL(0.8) THICKN(0.7)'
c
c     low energy sum rule
c
       ewsrlow = zero           ! accumulator for low energy sum rule
       do idat=1,maxdat
         x = xin(idat)
         y = yin(idat)
         if(x.gt.cutsr0) then
           ewsrlow = ewsrlow 
     &              + yin(idat)*0.5D0*(erf((cutsr1-x)/widsr)+1D0)
     &                /xin(idat)
         endif
       enddo
       write(2,'(a,f4.1,a,1pg13.5)')
     &  '# EWSR(',cutsr1,')    :',ewsrlow
       open(66,file='ewsrlow.res')
       write(66,'(1pg24.16)') ewsrlow
       close(66)
c      
      endif
 1000 continue
cc       open(unit=2,status='unknown',file='strength_full.res')
      close(2)
c
      return
c
  991 stop ' error in one of the heading lines'
c
  992 stop ' error in column specifier line'
c
  993 stop ' error in data line'
c
  999 return
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
