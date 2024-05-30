c-----evaluate_correl------------------------------------------------
c
      program evaluate_correl
c
      implicit double precision (a-h,o-z)
c
c     reads variance information from file and
c     computes the variance matrix of two selected
c     observables.
c
c
      include 'comfit.inc'
c
      character*3 ident
c
c     'workdv' stores the gradients for observables:
c       1  E/A
c       2  rho
c       3  K
c       4  m*/m
c       5  asym
c       6  dasym
c       7  kappa
c       8  E/A_neut at 0.1
c       9  slope of E/A_neut at 0.1
c      10  GMR
c      11  GDR
c      12  GQR
c      13  GDR in 16O
c      14  dipole polarizability in 208Pb
c      15  skin in 208Pb
c      16  skin in 124Sn
c      17  skin in 132Sn
c      18  skin in 140Sn
c      19  skin in 48Ca
c      20  polariz. in 132Sn
c      21  polariz. in 48Ca
c      22  polariz. in 68Ni
c      23  pygmy resonance in 208Pb
c      24  skin in 68Ni
c      25  pygmy resonance in 132Sn
c      26  diffraction skin in 208Pb
c      27  neutron rms  in 208Pb
c      28  neutron diffraction radius  in 208Pb
c      29  neutron formfactor at fixed Q  in 208Pb
c      30  binding energy in very neutron rich 148Sn
c
c     pointers sort the observables in closer groups
c       1 <-  1  E/A
c       2 <-  2  rho
c       3 <-  3  K
c       4 <- 10  GMR
c       5 <-  4  m*/m
c       6 <- 12  GQR
c       7 <-  7  kappa
c       8 <- 11  GDR
c       9 <- 13  GDR in 16O
c      10 <-  5  asym
c      11 <-  6  dasym
c      12 <-  8  E/A_neut at 0.1
c      13 <-  9  slope of E/A_neut at 0.1
c      14 <- 30  binding energy in very neutron rich 148Sn
c      15 <- 14  dipole polarizability in 208Pb
c      16 <- 20  polariz. in 132Sn
c      17 <- 22  polariz. in 68Ni
c      18 <- 21  polariz. in 48Ca
c      19 <- 15  skin in 208Pb
c      20 <- 16  skin in 124Sn
c      21 <- 17  skin in 132Sn
c      22 <- 18  skin in 140Sn
c      23 <- 24  skin in 68Ni
c      24 <- 19  skin in 48Ca
c      25 <- 26  diffraction skin in 208Pb
c      26 <- 27  neutron rms  in 208Pb
c      27 <- 28  neutron diffraction radius  in 208Pb
c      28 <- 29  neutron formfactor at fixed Q  in 208Pb
c      29 <- 23  pygmy resonance in 208Pb
c      30 <- 25  pygmy resonance in 132Sn
c
      logical tvarellips,talignmat
      parameter (tvarellips=.true.,talignmat=.true.)
      parameter (kobs=30)
      dimension workdv(0:kparmx,kobs),alig(kobs*kobs)
      dimension almat(kobs,kobs)
      dimension ipo(2,kobs*kobs)
      dimension iplpoi(kobs),irevpo(kobs)
      data iplpoi/01,02,03,10,04,12,07,11,13,05,
     &            06,08,09,30,14,20,22,21,15,16,
     &            17,18,24,19,26,27,28,29,23,25/
      data irevpo/01,02,03,05,10,11,07,12,13,04,
     &            08,06,09,15,19,20,21,22,24,16,
     &            18,17,29,23,30,25,26,27,28,14/
c
      character*8 obs(kobs)
      character*2 aobs(kobs)
      data obs/'  E/A   ','  rho   ','   K    ','  m*/m  ',
     &         '  asym  ','  dasym ','  kappa ',' E/Aneut',
     &         'slopneut','  GMR   ','  GDR   ','  GQR   ',
     &         '  GDR16 ',' polariz','  skin  ','skinsn74',
     &         'skinsn82','skinsn90','skinca48','pol132sn',
     &         ' pol48ca',' pol68ni',' pygmyPb','skinni68',
     &         ' pygmySn',' dskin  ',' rmsn   ',' rdsmn  ',
     &         ' formfn ','E_B148Sn'/
      data aobs
     &        /'01','02','03','04','05','06','07','08','09',
     &    '10','11','12','13','14','15','16','17','18','19',
     &    '20','21','22','23','24','25','26','27','28','29',
     &    '30'/
c
      logical ttest
      parameter (ttest=.true.)
c
c---------------------------------------------------------------------
c
c     check pointers
c
      do i=1,kobs
        j = iplpoi(irevpo(i))
        if(j.ne.i) then
          write(6,'(a,2i5)') ' wrong pointer: i,j=',i,j
          stop
        endif
      enddo     
c
c     read 'dparmo' and 'zarray' from data file
c
      open(31,file='erroranalysis')
      read(31,'(15x,i5)') nparmx
      if(ttest) write(6,'(a,i5)') ' nparmx=',nparmx
      do i=1,nparmx
        read(31,'(100g24.16)') (zarray(j,i),j=1,nparmx)
      enddo
c
      if(ttest) write(*,*) ' after zarray'
c                                               NMP
      read(31,'(1x)')
      do i=0,nparmx
        read(31,*) (workdv(i,j),j=1,9)
      enddo
c                                               GMR
      read(31,'(1x)')
      read(31,*) workdv(0,10)
      read(31,*) (workdv(i,10),i=1,nparmx)
c                                               GDR
      read(31,'(1x)')
      read(31,*) workdv(0,11)
      read(31,*) (workdv(i,11),i=1,nparmx)
c                                               GQR
      read(31,'(1x)')
      read(31,*) workdv(0,12)
      read(31,*) (workdv(i,12),i=1,nparmx)
c                                               GDR16
      read(31,'(1x)')
      read(31,*) workdv(0,13)
      read(31,*) (workdv(i,13),i=1,nparmx)
c                                               polariz
      read(31,'(1x)')
      read(31,*) workdv(0,14)
      read(31,*) (workdv(i,14),i=1,nparmx)
c                                               skinpb
      read(31,'(1x)')
      read(31,*) workdv(0,15)
      read(31,*) (workdv(i,15),i=1,nparmx)
c                                               skinsn74
      read(31,'(1x)')
      read(31,*) workdv(0,16)
      read(31,*) (workdv(i,16),i=1,nparmx)
c                                               skinsn82
      read(31,'(1x)')
      read(31,*) workdv(0,17)
      read(31,*) (workdv(i,17),i=1,nparmx)
c                                               skinsn90
      read(31,'(1x)')
      read(31,*) workdv(0,18)
      read(31,*) (workdv(i,18),i=1,nparmx)
c                                               skini48ca
      read(31,'(1x)')
      read(31,*) workdv(0,19)
      read(31,*) (workdv(i,19),i=1,nparmx)
c                                               pol132sn
      read(31,'(1x)')
      read(31,*) workdv(0,20)
      read(31,*) (workdv(i,20),i=1,nparmx)
c                                                pol48ca
      read(31,'(1x)')
      read(31,*) workdv(0,21)
      read(31,*) (workdv(i,21),i=1,nparmx)
c                                                pol68Ni
      read(31,'(1x)')
      read(31,*) workdv(0,22)
      read(31,*) (workdv(i,22),i=1,nparmx)
c                                                pygmy
      read(31,'(1x)')
      read(31,*) workdv(0,23)
      read(31,*) (workdv(i,23),i=1,nparmx)
c                                               skinni68
      read(31,'(1x)')
      read(31,*) workdv(0,24)
      read(31,*) (workdv(i,24),i=1,nparmx)
c                                               pygmy132Sn
      read(31,'(1x)')
      read(31,*) workdv(0,25)
      read(31,*) (workdv(i,25),i=1,nparmx)
c                                               dskin
      read(31,'(1x)')
      read(31,*) workdv(0,26)
      read(31,*) (workdv(i,26),i=1,nparmx)
c                                               rmsn
      read(31,'(1x)')
      read(31,*) workdv(0,27)
      read(31,*) (workdv(i,27),i=1,nparmx)
c                                               rdmsn
      read(31,'(1x)')
      read(31,*) workdv(0,28)
      read(31,*) (workdv(i,28),i=1,nparmx)
c                                               formfn
      read(31,'(1x)')
      read(31,*) workdv(0,29)
      read(31,*) (workdv(i,29),i=1,nparmx)
c                                               E_B 148Sn
      read(31,'(1x)')
      read(31,*) workdv(0,30)
      read(31,*) (workdv(i,30),i=1,nparmx)
c                                               
      close (31)
c
c     select, compute and print
c
      if(tvarellips) then
      open(31,file='varellips.dat')
c
      write(31,'(a)') ' polariz vs. skin:'
      call varellips(workdv(0,14),workdv(0,15),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' polariz vs. skin: align=',align
c
      write(31,'(a)') ' asym vs. skin:'
      call varellips(workdv(0,5),workdv(0,15),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' asym vs. skin:    align=',align
c
      write(31,'(a)') ' dasym vs. skin:'
      call varellips(workdv(0,6),workdv(0,15),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' dasym vs. skin:   align=',align
c
      write(31,'(a)') ' kappa vs. skin:'
      call varellips(workdv(0,7),workdv(0,15),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' kappa vs. skin:   align=',align
c
      write(31,'(a)') ' K vs. skin:'
      call varellips(workdv(0,3),workdv(0,15),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' K vs. skin:       align=',align
c
      write(31,'(a)') ' m*/m vs. skin:'
      call varellips(workdv(0,4),workdv(0,15),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' m*/m vs. skin:    align=',align
c
      write(31,'(a)') ' slope_neut vs. skin:'
      call varellips(workdv(0,9),workdv(0,15),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' slope_nt vs. skin:align=',align
c
      write(31,'(a)') ' GDR vs. skin:'
      call varellips(workdv(0,11),workdv(0,15),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' GDR vs. skin:    align=',align
c
      write(31,'(a)') ' asym vs. daysm:'
      call varellips(workdv(0,5),workdv(0,6),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' asym vs. daysm:  align=',align
c
      write(31,'(a)') ' polariz  vs. GDR:'
      call varellips(workdv(0,14),workdv(0,11),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' polariz vs. GDR: align=',align
c
      write(31,'(a)') ' polariz  vs. asym:'
      call varellips(workdv(0,14),workdv(0,5),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' polariz vs. asym:align=',align
c
      write(31,'(a)') ' polariz  vs. dasym:'
      call varellips(workdv(0,14),workdv(0,6),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' polariz vs.dasym:align=',align
c
      write(31,'(a)') ' polariz  vs. slope_neut:'
      call varellips(workdv(0,14),workdv(0,9),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' polariz vs.slope:align=',align
c
      write(31,'(a)') ' skinsn74 vs. skinsn82:'
      call varellips(workdv(0,16),workdv(0,17),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' skinsn74 vs. skinsn82: align=',align
c
      write(31,'(a)') ' skinsn74 vs. skinsn90:'
      call varellips(workdv(0,16),workdv(0,18),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' skinsn74 vs. skinsn90: align=',align
c
      write(31,'(a)') ' skinsn82 vs. skinsn90:'
      call varellips(workdv(0,17),workdv(0,18),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' skinsn82 vs. skinsn90: align=',align
c
      write(31,'(a)') ' skinpb vs. skinsn82:'
      call varellips(workdv(0,15),workdv(0,17),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' skinpb vs. skinsn82: align=',align
c
      write(31,'(a)') ' skinSn vs. pygmySn:'
      call varellips(workdv(0,17),workdv(0,25),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') ' skinSn vs. pygmySn: align=',align
c
      write(31,'(a)') ' skinCa vs. polarizCa:'
      call varellips(workdv(0,19),workdv(0,21),nparmx,31,align)
      write(31,'(5x)')
      write(31,'(5x)')
      write(6,'(a,f8.5)') '  skinCa vs. polarizCa: align=',align
c
      close(31)
      endif
c
c     systematics of correlations
c
      write(7,'(a)') ' alignement matrix:'
      is = 0
      do i=1,kobs
      do j=1,i-1
        is = 1+is
        call varellips(workdv(0,i),workdv(0,j),nparmx,31,alig(is))
        ipo(1,is) = i
        ipo(2,is) = j
      enddo
      write(7,'(100(f7.4))') (alig(k),k=is-i+2,is)
      enddo
c
c     output of correlation matrix for color scale plot
c
      if(talignmat) then

        is = 0
        do i=1,kobs
          do j=1,i-1
            is = 1+is
            almat(i,j) = alig(is)
            almat(j,i) = alig(is)
          enddo
          almat(i,i) = 1D0
        enddo

        open(11,file='alignmatrix.gnu')
        write(11,'(a)') 
     &    'set term post eps enhanced color "roman bold" 16'
        write(11,'(a)') 'set output "alignmatrix.eps"'
        write(11,'(a)') 'set size 1,1.4'
        write(11,'(a)') 'set pm3d map'
        write(11,'(a)') 'set xrange [0:31]'
        write(11,'(a)') 'set yrange [0:31]'
        write(11,'(a)') 'set cbrange [0:1]'
        write(11,'(a,30(a,i2,a,i3,a))') 
     &  'set ytics (',('"',i,': '//obs(iplpoi(i))//'"',i,','
     &              ,i=1,kobs-1),
     &         '"'//obs(iplpoi(kobs))//'"',kobs,')'
        write(11,'(a)') 'sp "alignmatrix.gnu" u 1:2:3 not'
c        write(11,'(a)') 'pause -1'
        write(11,'(a)') 'q'
        do i=1,kobs
          xi = i-0.5D0
          ia = iplpoi(i)
          do j=1,kobs
            xj = j-0.5D0
            ja = iplpoi(j)
            write(11,'(2f8.1,1pg13.5)') xi,xj,almat(ia,ja)
            write(11,'(2f8.1,1pg13.5)') xi,xj+1D0,almat(ia,ja)
          enddo
          write(11,'(5x)')
          do j=1,kobs
            xj = j-0.5D0
            ja = iplpoi(j)
            write(11,'(2f8.1,1pg13.5)') xi+1D0,xj,almat(ia,ja)
            write(11,'(2f8.1,1pg13.5)') xi+1D0,xj+1D0,almat(ia,ja)
          enddo
          write(11,'(5x)')
          write(11,'(5x)')
        enddo
        close(11)

      endif
c
c     output sorted for each observable separately
c
      do i=1,is
        aligmin = alig(i)
        do j=i+1,is
          if(alig(j).lt.aligmin) then
            a = alig(j)
            alig(j) = aligmin
            aligmin = a
            alig(i) = a
            isav = ipo(1,j)
            ipo(1,j) = ipo(1,i)
            ipo(1,i) = isav
            isav = ipo(2,j)
            ipo(2,j) = ipo(2,i)
            ipo(2,i) = isav
          endif
        enddo
      enddo
c
      do i=1,is
        write(7,'(a,1x,a,f8.5)')
     &    obs(ipo(1,i)), obs(ipo(2,i)),alig(i)
      enddo
c
      do j=1,kobs
        ja = iplpoi(j)
        jb = irevpo(j)
        open(11,file='align.'//aobs(j))
        write(11,'(a,i3,a)')
     &    '#alignement for observable:',j,obs(ja)
        do i=is,1,-1
          if(irevpo(ipo(1,i)).eq.j) then
            write(11,'(2i3,1x,a,f8.5)')
     &         irevpo(ipo(2,i)),ipo(2,i),obs(ipo(2,i)),alig(i)
          elseif(irevpo(ipo(2,i)).eq.j) then
            write(11,'(2i3,1x,a,f8.5)')
     &         irevpo(ipo(1,i)),ipo(1,i),obs(ipo(1,i)),alig(i)
          endif
        enddo
        close(11)
      enddo
c
      stop
      end

c-----varellips -----------------------------------------------------
c
      subroutine varellips(grad1,grad2,nparmx,iu,align)
c
      implicit double precision (a-h,o-z)
c
      include 'comfit.inc'
c
c     computes 2x2 variance matrix for observables on 'grad1'
c     and 'grad2'. The correlation matrix 'zarray' comes 
c     through 'comfit.inc'.
c
      dimension grad1(0:nparmx),grad2(0:nparmx)
c 
      parameter (kplot=50)
c
c---------------------------------------------------------------
c
c
c     compute variances
c
      err11   = 0.0
      err22   = 0.0
      err12   = 0.0
c
      if(nparmx.gt.0) then
        do i=1,nparmx
        do j=1,nparmx
          err11  = grad1(i)*grad1(j)*zarray(i,j)
     &              + err11
          err22  = grad2(i)*grad2(j)*zarray(i,j)
     &              + err22
          err12  = grad1(i)*grad2(j)*zarray(i,j)
     &              + err12
        enddo
        enddo
      endif
      call  diag2(err11,err12,err22,eig1,eig2,angle)
c
c     print eigenvalues
c
      write(6,'(a,3(1pg13.5)/a,3(1pg13.5))')
     & ' variance  matrix:',err11,err22,err12,
     & ' eigenvalues etc :',eig1,eig2,angle/3.1415926*180D0
c
c     compute error ellipsoid and print
c
      align = err12**2/(err11*err22)
      det = err11*err22-err12**2
      erri11 = err22/det
      erri22 = err11/det
      erri12 = -err12/det
      xlim = sqrt(abs(erri22))/sqrt(erri11*erri22-erri12**2)
      write(6,*) erri11,erri22,erri12,xlim
      do j=-kplot,kplot
        x = j*xlim/kplot
        diff = sqrt(max(1D0/erri22
     &              +x**2*(erri12**2-erri11*erri22)/erri22**2,
     &              1D-20))
        sum  = -erri12*x/erri22
        write(iu,'(3(1pg13.5))')
     &      x+grad1(0),sum+diff+grad2(0)
      enddo      
      do j=-kplot+1,kplot
        x = -j*xlim/kplot
        diff = sqrt(max(1D0/erri22
     &              +x**2*(erri12**2-erri11*erri22)/erri22**2,
     &              1D-20))
        sum  = -erri12*x/erri22
        write(iu,'(3(1pg13.5))')
     &      x+grad1(0),sum-diff+grad2(0)
      enddo      
      write(iu,'(a,1pg13.5)') '# alignement: ',align
c
      return
      end

c-----diag2 -----------------------------------------------------
c
      subroutine diag2(a11,a12,a22,eig1,eig2,angle)
c
      implicit double precision (a-h,o-z)
c
c     diagonalizes 2x2 matrix, 
c     returns eigenvectors and angle of rotation
c
c---------------------------------------------------------------
c
      deig = sqrt((0.5D0*(a11-a22))**2+a12**2)
      eig1 = 0.5D0*(a11+a22)+deig
      eig2 = 0.5D0*(a11+a22)-deig
      angle = atan((a11-eig1)/a12)
c
      return
      end
ctest c
ctest       implicit double precision (a-h,o-z)
ctest c
ctest       angin = 3.1415926D0/6.0
ctest       c = cos(angin)
ctest       s = sin(-angin)
ctest       call diag2(c*c,s*c,s*s,eig1,eig2,angle)
ctest       write(6,'(5(1pg13.5))') angin,eig1,eig2,angle,
ctest      &  angle+3.1415926D0/4D0
ctest 
ctest     
ctest       stop
ctest       end
