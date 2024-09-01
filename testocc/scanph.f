c-----scanph-------------------------------------------------------------------
c
      subroutine scanph(ecutlw,ecutup,lmode,nparmode)
c
      include 'params.inc'
c
c     determines the ph states in the valence space for
c     angular momentum 'lmode' and parity 'nparmode'.
c     'ecutlw' and 'ecutup' set the energy window for the
c     search.
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'result.inc'
      include 'skypar.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      INTEGER,PARAMETER :: kmaxph=50
      REAL(8) :: xsav(3,kmaxph),x2sav(3),nl2(0:4)
      INTEGER :: nlsav(0:4,kmaxph)
c
c--------------------------------------------------------------------
c
c     write header
c
      write(11,'(//a,2i5)') ' 1ph spectra for L,parity=',
     &   lmode,nparmode
      write(11,'(a)') 'iso na   ja la  nb   jb lb     eph'
c
c     provide necessary quantum numbers and pointers with
c
      DO iso=1,2
        IF(iso.EQ.1) THEN
         nmax   = nmaxp
         nstack = 0
         isoadd = 0
        ELSE
         nmax   = nmaxn
         nstack = nmaxp*ngrid
         isoadd = ngrid
        ENDIF
c
        nmaxph=0
        DO na=1,nmax
          jsa   = jsp(na)
          la    = jsa/2
          xla   = la
          xja   = ((jsa-1)/2)+0.5
          wgt   = weight(na+nstack/ngrid)
          IF(wgt.LT.0.5D0) CYCLE
          nadd  = nstack + (na-1)*ngrid
          subfac = spenrg(na+nstack/ngrid)
c          write(11,*) ' na,ja,la,xj=',na,ja,la,xj
          if(wgt.le.half) CYCLE
c
          DO nb=1,nmax
            jsb   = jsp(nb)
            lb    = jsb/2
            xlb   = lb
            xjb   = (jsb-1)/2+0.5D0
            wgtb  = weight(nb+nstack/ngrid)
c
c           selection rules
c
            IF(wgtb.gt.half) CYCLE
            IF(lmode.LT.NINT(ABS(xja-xjb))) CYCLE
            IF(lmode.GT.NINT(xja+xjb)) CYCLE
            ndifpar = 1-2*MOD(la+lb,2)
            IF(nparmode.NE.ndifpar) CYCLE
c
c           values and print
c
            subfcb = spenrg(nb+nstack/ngrid)  
            npart=nb
            nhole=na
            edif = ABS(subfac-subfcb)
            if(edif.lt.ecutup .and. edif.gt.ecutlw
     &         .and. subfac.ne.zero .and. subfcb.ne.zero) THEN
              nmaxph=1+nmaxph
              xsav(1,nmaxph) = xja
              xsav(2,nmaxph) = xjb
              xsav(3,nmaxph) = edif
              nlsav(0,nmaxph) = iso
              nlsav(1,nmaxph) = na
              nlsav(2,nmaxph) = la
              nlsav(3,nmaxph) = nb
              nlsav(4,nmaxph) = lb
            END IF
c
          END DO
        END DO
c
!        WRITE(11,'(a,2i5,200(1pg13.5))')
!     &   ' iso,nmaxph,edifs=',iso,nmaxph,xsav(3,1:nmaxph)
        DO na=1,nmaxph
          emin=xsav(3,na)
          DO nb=na+1,nmaxph
           IF(xsav(3,nb).LT.emin) THEN
             x2sav = xsav(:,na)
             nl2 = nlsav(:,na)
             xsav(:,na) = xsav(:,nb)
             nlsav(:,na) = nlsav(:,nb)
             xsav(:,nb) = x2sav
             nlsav(:,nb) = nl2   
             emin = xsav(3,na)
           END IF
           END DO
        END DO
        DO na=1,nmaxph
          write(11,'(i2,2(i4,f5.1,i3),1pg13.5)') nlsav(0,na), 
     &         nlsav(1,na),xsav(1,na),nlsav(2,na),
     &         nlsav(3,na),xsav(2,na),nlsav(4,na),xsav(3,na)
        END DO
c
c       uncommented output of lowest state for later analysis
c 
        WRITE(36,'(i3,i2,2(i4,f5.1,i3),1pg20.12)') lmode,nlsav(0,1), 
     &         nlsav(1,1),xsav(1,1),nlsav(2,1),
     &         nlsav(3,1),xsav(2,1),nlsav(4,1),xsav(3,1)
c
      END DO
c
c
      return
      end
c-----dos2ph-------------------------------------------------------------------
c
      subroutine dos2ph(ecutlw,ecutup,l2ph)
c
      include 'params.inc'
c
c     Determines the density of 2ph states in the valence space
c     for angular momentum 'l2ph'. Partity is natural parity
c     determined from angular momentum.
c     'ecutlw' and 'ecutup' set the energy window for the
c     search.
c     Builds forst a table of 1ph states from which the
c     2ph states for angular momentum 'l2ph' are selected.
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'result.inc'
      include 'skypar.inc'
      include 'grid.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      INTEGER,PARAMETER :: kmaxph=200,klmaxph=15
      LOGICAL,PARAMETER :: ttest=.FALSE.
!      REAL(8),PARAMETER :: ebin=0.1D0,emargin=2D1*ebin,wid=4D0*ebin
      REAL(8),PARAMETER :: ebin=0.1D0,emargin=2D1*ebin,wid=ebin
      REAL(8) :: xsav(4,kmaxph,0:klmaxph),x2sav(4),nl2(0:4)
      INTEGER :: nlsav(0:4,kmaxph,0:klmaxph),nmax1ph(0:klmaxph)
      REAL(8),ALLOCATABLE :: dos2phsav(:)
c
c--------------------------------------------------------------------
c
c
c     First step: generate tables of 1ph states.
c
c     provide necessary quantum numbers and pointers with
c
      DO l1ph=0,klmaxph
        npar1ph=1-2*MOD(l1ph,2)
        DO iso=1,2
          IF(iso.EQ.1) THEN
           nmax   = nmaxp
           nstack = 0
           isoadd = 0
          ELSE
           nmax   = nmaxn
           nstack = nmaxp*ngrid
           isoadd = ngrid
          ENDIF
c  
          nmaxph=0
          DO na=1,nmax
            jsa   = jsp(na)
            la    = jsa/2
            xla   = la
            xja   = ((jsa-1)/2)+0.5
            wgt   = weight(na+nstack/ngrid)
!            IF(wgt.LT.0.5D0) CYCLE
            nadd  = nstack + (na-1)*ngrid
            subfac = spenrg(na+nstack/ngrid)
c            write(11,*) ' na,ja,la,xj=',na,ja,la,xj
!            if(wgt.le.half) CYCLE
c
            DO nb=1,nmax
              jsb   = jsp(nb)
              lb    = jsb/2
              xlb   = lb
              xjb   = (jsb-1)/2+0.5D0
              wgtb  = weight(nb+nstack/ngrid)
c
c             selection rules
c
              IF(wgtb.LE.wgt) CYCLE
              IF(l1ph.LT.NINT(ABS(xja-xjb))) CYCLE
              IF(l1ph.GT.NINT(xja+xjb)) CYCLE
              ndifpar = 1-2*MOD(la+lb,2)
              IF(npar1ph.NE.ndifpar) CYCLE
c
c           values and print
c
              subfcb = spenrg(nb+nstack/ngrid)  
              npart=nb
              nhole=na
              edif = ABS(subfac-subfcb)
              if(edif.lt.ecutup .and. edif.gt.0D0
     &           .and. subfac.ne.zero .and. subfcb.ne.zero) THEN
                nmaxph=1+nmaxph
                xsav(1,nmaxph,l1ph) = xja
                xsav(2,nmaxph,l1ph) = xjb
                xsav(3,nmaxph,l1ph) = edif
                xsav(4,nmaxph,l1ph) = ABS(wgt-wgtb)
                nlsav(0,nmaxph,l1ph) = iso
                nlsav(1,nmaxph,l1ph) = na
                nlsav(2,nmaxph,l1ph) = la
                nlsav(3,nmaxph,l1ph) = nb
                nlsav(4,nmaxph,l1ph) = lb
              END IF
c  
            END DO
          END DO
          nmax1ph(l1ph) = nmaxph
c
c       optionally write protocol of 1ph states
c
          IF(ttest) THEN
            WRITE(11,'(a,2i4/a)') '1ph states for iso,L=',iso,l1ph,
     &      ' iso   na   la  nb   lb  epsilona  epsilonb edif'
            DO na=1,nmaxph
              emin=xsav(3,na,l1ph)
              DO nb=na+1,nmaxph
               IF(xsav(3,nb,l1ph).LT.emin) THEN
                 x2sav = xsav(:,na,l1ph)
                 nl2 = nlsav(:,na,l1ph)
                 xsav(:,na,l1ph) = xsav(:,nb,l1ph)
                 nlsav(:,na,l1ph) = nlsav(:,nb,l1ph)
                 xsav(:,nb,l1ph) = x2sav
                 nlsav(:,nb,l1ph) = nl2   
                 emin = xsav(3,na,l1ph)
               END IF
               END DO
            END DO
            DO na=1,nmaxph
              write(11,'(i2,2(i4,f5.1,i3),2(1pg13.5))')nlsav(0,na,l1ph), 
     &             nlsav(1,na,l1ph),xsav(1,na,l1ph),nlsav(2,na,l1ph),
     &             nlsav(3,na,l1ph),xsav(2,na,l1ph),nlsav(4,na,l1ph),
     &             xsav(3,na,l1ph),xsav(4,na,l1ph)
            END DO
          END IF
        END DO
      END DO
c
c     now compute DoS for 2ph states
c 
      OPEN(111,file='DoS_2ph',position='APPEND')
      espan=ecutup-ecutlw-2D0*emargin
      nbinmax=NINT(espan/ebin)
!      ebin=espan/nbinmax
      ALLOCATE(dos2phsav(nbinmax))
      dos2phsav=0D0
      nbinmarg=NINT(emargin/ebin)
      IF(ttest) WRITE(6,'(a,2f8.2,2i4)') 'start DoS:',
     &    espan,emargin,nbinmax,nbinmarg
      WRITE(111,'(a)') '# DoS of 2ph states'
      CALL flush(111)
      DO l1=0,klmaxph-1
        IF(ttest) WRITE(6,'(a,i4)') 'start l1=',l1
        l2=l1+1
        DO n1=1,nmax1ph(l1); DO n2=1,nmax1ph(l2)
          n1a=nlsav(1,n1,l1)
          n1b=nlsav(3,n1,l1)
          l1a=nlsav(2,n1,l1)
          l1b=nlsav(4,n1,l1)
          x1a=xsav(1,n1,l1)
          x1b=xsav(2,n1,l1)
          iso1=nlsav(0,n1,l1)
          n2a=nlsav(1,n2,l2)
          n2b=nlsav(3,n2,l2)
          l2a=nlsav(2,n2,l2)
          l2b=nlsav(4,n2,l2)
          x2a=xsav(1,n2,l2)
          x2b=xsav(2,n2,l2)
          iso2=nlsav(0,n2,l2)
          IF(n1a.EQ.n2a .AND. l1a.EQ.l2a .AND. x1a.EQ.x2a .AND.
     &       iso1.EQ.iso2) CYCLE
          IF(n1a.EQ.n2b .AND. l1a.EQ.l2b .AND. x1a.EQ.x2b .AND.
     &       iso1.EQ.iso2) CYCLE
          IF(n1b.EQ.n2a .AND. l1b.EQ.l2a .AND. x1b.EQ.x2a .AND.
     &       iso1.EQ.iso2) CYCLE
          IF(n1b.EQ.n2b .AND. l1b.EQ.l2b .AND. x1b.EQ.x2b .AND.
     &       iso1.EQ.iso2) CYCLE
          e2ph=xsav(3,n1,l1)+xsav(3,n2,l2)
          nbin2ph=NINT(e2ph/ebin)
          IF(ttest) WRITE(6,'(a,3i4,2f8.2)')
     &       'in loop:',n1,n2,nbin2ph,e2ph,nbin2ph*ebin+ecutlw
          DO nbin=nbin2ph-nbinmarg,nbin2ph+nbinmarg
            IF(nbin>nbinmax .OR. nbin<1) CYCLE
            eact=nbin*ebin+ecutlw
            dos2phsav(nbin)=dos2phsav(nbin)+
     &       xsav(4,n1,l1)*xsav(4,n2,l2)*EXP(-((eact-e2ph)/wid)**2)
          END DO
        END DO; END DO        
        IF(ttest) WRITE(6,'(a,i4)') 'done l1=',l1
!        CALL flush(6)
      END DO
      IF(ttest) WRITE(6,'(a)') 'end of accumulation'
!      CALL flush(6)
c
      DO nbin=1,nbinmax
        eact=nbin*ebin+ecutlw+emargin
        WRITE(111,'(i3,2(1pg13.5))') l2ph,eact,dos2phsav(nbin)
      END DO
      CLOSE(111)
c
      DEALLOCATE(dos2phsav)
c
      return
      end
