      parameter (kobs=30)
      character*40 atext
      character*22 obs(kobs)
      data obs/'  E/A                 ',
     &         '  rho                 ',
     &         '   K                  ',
     &         '  m*/m                ',
     &         '  asym                ',
     &         '  d/drho asym         ',
     &         '  kappa               ',
     &         ' E/A neutrons         ',
     &         ' slope E/A neutrons   ',
     &         '  GMR 208Pb           ',
     &         '  GDR 208Pb           ',
     &         '  GQR 208Pb           ',
     &         '  GDR 16O             ',
     &         ' polarizability 208Pb ',
     &         ' rms skin 208Pb       ',
     &         ' rms skin 74Sn        ',
     &         ' rms skin 82Sn        ',
     &         ' rms skin 90Sn        ',
     &         ' rms skin 48Ca        ',
     &         ' polarizability 132Sn ',
     &         ' polarizability 48Ca  ',
     &         ' polarizability 68Ni  ',
     &         ' pygmy resonance 208Pb',
     &         ' rms skin 68Ni        ',
     &         ' pygmy resonance 132Sn',
     &         ' diffr.skin 208Pb     ',
     &         ' rms. radius neutron  ',
     &         ' diffr. radius neutron',
     &         ' formf.neutron(q=0.45)',
     &         ' binding energy 148Sn '/
      dimension aligns(kobs),iobs(kobs)
c
c----------------------------------------------------
c
      read(5,'(a)') atext
      write(6,'(a)') atext
      do i=1,kobs-1
        read(5,'(i3,9x,f8.5)') iobs(i),aligns(i)
      enddo
c
      write(6,'(2a)')
     & 'set terminal postscript eps enhanced color solid 30'
c     &' "roman bold" 30'
      write(6,'(a)')
     &'set output "align-skin.eps"',
     &'set border 15 lw 3',
     &'set tmargin 0',
     &'set bmargin 0',
     &'set rmargin 0',
     &'set lmargin 0',
     &'set size 1.03,1.03',
     &'set form x ""',
     &'set yrange [80:100]',
     &'set xrange [0.3:15.7]',
     &'set multiplot',
     &'set size 0.8,0.9',
     &'set origin 0.2,0.1',
     &'set ylabel "{/bold= alignement with     [%]}"'
      do i=1,kobs-1
       write(6,'(3a,i2,a)')
     &  'set label "',obs(iobs(i)),'" at ',i,',81 rotate by 90'        
      enddo
      write(6,'(a)')
     &'p "align-skin.gnu" u 1:3 not w boxes lw 9',
     &'q'  
      write(6,'(a)') atext
      do i=1,kobs-1
       write(6,'(2i3,f8.3,a)')
     &     i,iobs(i),aligns(i)*100,obs(iobs(i))
      enddo
c
      stop
      end

