c-----wavelet-----------------------------------------------------------
c
      program wavelet
c
c     Spectral analysis with cos^2 wavelets
c     from input file which contains time in column one and
c     signal in column 2.
c
      parameter (KDIM=30000)
      parameter (KWIDMAX=20)
      parameter (KOMGFAC=1)
      parameter (OMEGAMAX=40.0)
      parameter (PI=3.1415926)
      logical   TMTV
      parameter (TMTV=.true.)
      logical   TEXTEND
      parameter (TEXTEND=.false.)
c
      dimension signal(0:KDIM),time(0:KDIM),filter(0:KDIM)
c
c-----------------------------------------------------------------------
c
      read(5,'(1x)')
      read(5,'(1x)')
      read(5,'(1x)')
      do i=0,KDIM
        read(5,*,err=99,end=99) time(i),signal(i)
        ntimax = i
      enddo
      stop ' file count exhausted'
c
 99   continue
      if(TEXTEND .and. ntimax.lt.KDIM/2) then
        do i=ntimax+1,ntimax*2
          signal(i) = 0.0
        enddo
        timax  = 2.0*time(ntimax)
        ntimax = 2*ntimax
      else
        timax  = time(ntimax)
      endif
c      write(6,*) ' ntimax,timax=',ntimax,timax
      deltim = timax/ntimax
      if(PI/deltim.lt.OMEGAMAX)
     &  stop ' OMEGAMAX cannot represented for given DELTIM'
      delomg = PI/timax*KOMGFAC
      do iomeg=0,ntimax
        omega = iomeg*delomg
        if(omega.gt.OMEGAMAX) goto 19
        nomegmax = iomeg
      enddo
 19   continue
c
      if(TMTV) then
        widmax = PI*KWIDMAX/timax
        widmin = widmax/KWIDMAX
        write(6,'(a)') '$ data=contour'
        write(6,'(a,i4,a,f8.4,a,f8.4,a,f8.4)')
     &   '% ny =',KWIDMAX,' ymin=',widmin,' ymax=',
     &       widmax,' vymin=0.0 vymax=',widmax
        write(6,'(a,i4,a,f8.4,a,f8.4)')
     &   '% nx =',nomegmax+1,' xmin= 0.0 xmax=',
     &       delomg*nomegmax,' vxmin=0.0 vxmax=',delomg*nomegmax
        write(6,'(a)')
     &    '% contstyle=2',
     &    '% logsteps=2',
     &    '% interp=2',
     &    '% xlabel="energy [MeV]"',
     &    '% ylabel="width [MeV]"',
     &    '% toplabel ="wavelet analysis of spectra"',
     &    '% grid=true',
     &    '% equalscale=true',
     &    '% fitpage=false',
     &    '% xyratio=0.8',
     &    '% cmin=1E-3',
     &    '% cmax=10'
      endif
c
      pih    = PI*0.5
      do nwid=1,KWIDMAX
        timaxact = timax/nwid
        widomg   = 1.0/timaxact
c
        ilim = timaxact/deltim
c        write(6,*) ' widomg=',widomg
        xnorm = 0.0
        do i=0,ntimax
          t = i*deltim
          if(i.le.ilim) then          
            filter(i) = signal(i)*cos(pih*t/timaxact)**2 !/sqrt(timaxact)
          else
            filter(i) = 0.0
          endif
          xnorm = filter(i)**2 + xnorm
        enddo
        xnorm = 1.0/sqrt(deltim*xnorm)
c
        do iomeg=0,nomegmax
          omega = iomeg*delomg
          sum = 0.0
          do i=0,ntimax
            sum = filter(i)*sin(omega*deltim*i) + sum
          enddo
          sum = sum*xnorm/sqrt(1.0*ntimax)
          if(TMTV) then
            write(6,'(4(1pg12.4))') sum*sum
          else
            write(6,'(4(1pg12.4))') widomg,omega,sum,sum*sum
          endif
        enddo
        write(6,'(a5)') '     '
        write(6,'(a5)') '     '
c
      enddo
c
      stop
      end

      

