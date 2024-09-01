c-----hfb_efermi------------------------------------------------------
c
      subroutine hfb_iterat()
c
      include 'params.inc'
c
c     Full HFB solution iterating the HFB equations together
c     with updates of self-consistent fields. Requires that
c     the switch tnatural=.TRUE. is set in subroutine HFB_eqs.
c
c     Test version 2. January 2018
c
      include 'wavfun.inc'
      include 'potinv.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'constn.inc'
      include 'result.inc'
c
      REAL(8) :: partnumgoal,efermiso
      REAL(8), PARAMETER :: pnslope=1D0,defermmin=0.1D0
      REAL(8), PARAMETER :: pnepsil=1D-10
      INTEGER, PARAMETER :: itermax=10
      LOGICAL, PARAMETER :: tprint=.TRUE.
c
c---------------------------------------------------------------------
c
      IF(nblockp.NE.0 .OR. nblockn.NE.0) STOP 'HFB not for odd nuclei'
c
!      IF(tprint) THEN
!        WRITE(6,'(a,100(/5(1pg13.5)))') 'HFB entry:'
!        WRITE(6,'(a,100(/5(1pg13.5)))') 'proton s.p.energies:',
!     &    spenrg(1:nmaxp)
!        WRITE(6,'(a,100(/5(1pg13.5)))') 'neutron s.p.energies:',
!     &    spenrg(nmaxp+1:nmaxp+nmaxn)
!        WRITE(6,'(a,100(/5(1pg13.5)))') 'U/V:',
!     &            SQRT((1D0-weight(nmaxp+1:nmaxp+nmaxn))/
!     &                 weight(nmaxp+1:nmaxp+nmaxn))
!        WRITE(6,'(a,100(/5(1pg13.5)))') 'V:',
!     &            SQRT(weight(nmaxp+1:nmaxp+nmaxn))
!        WRITE(6,'(a,100(/5(1pg13.5)))') 'U:',
!     &            SQRT(1D0-weight(nmaxp+1:nmaxp+nmaxn))
!      END IF
      DO iso=1,2
        ie=ngrid*iso
        im=ie-ngrid+1
        WRITE(6,'(a,4i3)') 'old densit: iso=',iso,im,ia,ie
        WRITE(6,'(a)') '# r  rho  tau  sj paircr  d2r  sjs'
        DO i=im,ie
          WRITE(6,'(f6.1,6(1pg13.5))') 
     &      rstep*(i-im),rho(i),tau(i),sj(i),paircr(i),d2r(i),sjs(i)
        END DO
      END DO
!      WRITE(6,'(a,100(/3(1pg13.5)))') '#old  rho_p  rho_n chi_n:',
!     &  (rho(i),rho(i+ngrid),paircr(i+ngrid),i=1,ngrid)
      acc=0D0
      accp=0D0
      r=0D0
      ia=ngrid+2
      ie=2*ngrid
      DO i=ia,ie
        r       = rstep+r
        r2      = r*r
        acc  = r2*rho(i) +acc
        accp = r2*paircr(i) +accp
      END DO
      acc = frpirs*acc
      accp = frpirs*accp
      WRITE(6,'(a,i2,2(1pg13.5))') 'HFB init: iso,norm=',
     &   iso,acc,accp
c
!        WRITE(6,'(a,100(/5(1pg13.5)))') 'paircr_n:',
!     &          (paircr(i),i=ngrid+1,2*ngrid)
      epair = epp+epn
      etotal = ekin+epotsk+epair
      IF(tprint) WRITE(6,'(a,4(1pg13.5))') 
     &   'HFB iteration before: etotal,epair=',etotal,epair,ekin,epotsk
      CALL hfbpot
!        WRITE(6,'(a,100(/5(1pg13.5)))') 'paircr_n:',
!     &          (paircr(i),i=ngrid+1,2*ngrid)
      iwr = -1
      DO iter=1,itermax
        DO iso=1,2
          IF(iso==1) THEN
            partnumgoal = aprot
            efermiso = efep
          ELSE
            partnumgoal = aneut
            efermiso = efen
          END IF
!          CALL  hfb_efermi(iso,efermiso,partnumgoal)
          CALL  hfb_natural(iso,efermiso,partnumgoal,100)
        END DO
!        call densit(weight)
        call skypot(iter)
        call hfbpot
        ekinacc  = zero
        r       = zero
        do i=2,ngrid
          r       = rstep+r
          r2      = r*r
          ekinacc  = r2*(amfc(1)*tau(i)+amfc(2)*tau(ngrid+i)) + ekinacc
        enddo
        write(6,'(2(1pg13.5))') frpirs,ekinacc
        ekin  = frpirs*ekinacc
        epair = epp+epn
        etotal = ekin+epotsk+epair
        IF(tprint) WRITE(6,'(a,i4,4(1pg13.5))') 
     &   'HFB iteration: iter,etotal,epair=',
     &    iter,etotal,epair,ekin,epotsk
      END DO

      RETURN
      END
c-----hfb_efermi------------------------------------------------------
c
      subroutine hfb_efermi(iso,efermiso,partnumgoal)
c
      include 'params.inc'
c
c     Iterates the Fermi energy in the HFB equations
c     to adjust wanted particle number.
c
      include 'wavfun.inc'
      include 'potinv.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
      REAL(8), INTENT(IN OUT) :: efermiso
      REAL(8), INTENT(IN) :: partnumgoal
      INTEGER, INTENT(IN) :: iso
c
      REAL(8) :: partnumout,pndiff
      REAL(8), PARAMETER :: pnslope=0.5D0,defermmin=0.1D0
      REAL(8), PARAMETER :: pnepsil=1D-10
      INTEGER, PARAMETER :: itermax=50
      LOGICAL, PARAMETER :: tprint=.FALSE.
c
c---------------------------------------------------------------------
c
c     trial steps to encircle particle number
c
      IF(tprint) WRITE(6,'(a,i3,3(1pg13.5))') 
     &   'HFB_Fermi for iso,eferm,part.goal=',iso,efermiso,partnumgoal
      CALL  hfb_eqs(iso,efermiso,partnumout)
      pndiff = partnumout-partnumgoal
      IF(tprint) WRITE(6,'(a,3(1pg13.5))')
     &  'HFB_Fermi init: part.nums=',partnumout,partnumgoal,pndiff
!      DO i=-20,20
!        eferm1=efermiso+i*0.0D0
!        CALL hfb_eqs(iso,eferm1,pn1)
!        WRITE(6,'(a,i3,2(1pg13.5))') 'test: i,eferm1,pn1=',
!     &     i,eferm1,pn1
!      END DO
      IF(partnumout>partnumgoal) THEN
        eferm2=efermiso
        pn2=partnumout
        eferm1=eferm2
        pnslopa=pnslope
        DO iter=1,itermax
          eferm1=eferm1-MAX(pndiff*pnslopa,iter*defermmin)
!          pnslopa=2D0*pnslopa
          CALL hfb_eqs(iso,eferm1,pn1)
          IF(tprint) WRITE(6,'(a,i3,4(1pg13.5))')
     &       'HFB_Fermi: iter,part.nums=',iter,pn1,partnumgoal
     &            ,pn1-partnumgoal,eferm1
          IF(pn1.LE.partnumgoal) GOTO 19
        END DO
        WRITE(6,'(a,3(1pg13.5))') 'HFB Fermi: eferm1,pn1,pnslopa=',
     &     eferm1,pn1,pnslopa
        STOP 'HFB_Fermi: not bracketed'
      ELSE IF(partnumout<partnumgoal) THEN
        eferm1=efermiso
        eferm2=eferm1
        pn1=partnumout
        DO iter=1,itermax
          eferm2=eferm2-MIN(pndiff*pnslope,-iter*defermmin)
          CALL hfb_eqs(iso,eferm2,pn2)
          IF(tprint) WRITE(6,'(a,i3,3(1pg13.5))')
     &       'HFB_Fermi: iter,part.nums=',iter,eferm2,pn2,partnumgoal
          IF(pn2.GE.partnumgoal) GOTO 19
        END DO
        STOP 'HFB_Fermi: not bracketed'
      ELSE
        WRITE(6,'(a)') 'particle number fits immediately'
        RETURN
      END IF
 19   CONTINUE
c
c     bracketing
c
      IF(tprint) WRITE(6,'(a,4(1pg13.5))') 'bracket:',eferm1,eferm2,
     &           pn1,pn2
      iterbrack=NINT(LOG((pn2-pn1)/pnepsil)/LOG(2D0))
      DO iter=1,iterbrack
        eferm3=0.5D0*(eferm1+eferm2)
        CALL hfb_eqs(iso,eferm3,pn3)
        pndiff=pn3-partnumgoal
        IF(tprint) WRITE(6,'(a,i3,3(1pg13.5))') 
     &    'HFB_Fermi bracket: iter,part.nums.=',iter,pn3,pndiff,eferm3
        IF(ABS(pndiff)<pnepsil) EXIT
        IF(pndiff<0D0) THEN
          eferm1=eferm3
          pn1=pn3
        ELSE
          eferm2=eferm3
          pn2=pn3
        END IF
      END DO
      efermiso=eferm3
      partnumout=pn3
      IF(ABS(pn2-partnumgoal)<ABS(pn3-partnumgoal)) THEN
        efermiso=eferm2
        partnumout=pn2
      ELSE IF(ABS(pn1-partnumgoal)<ABS(pn3-partnumgoal)) THEN
        efermiso=eferm1
        partnumout=pn1
      END IF
      WRITE(6,'(a,3(1pg13.5))')
     &  'HFB_Fermi finished: efermiso,part.nums.=',efermiso,
     &    partnumout,partnumgoal,partnumout-partnumgoal
c
      RETURN
      END
c-----hfb_eqs---------------------------------------------------------
c
      subroutine hfb_eqs(iso,efermiso,partnumout)
c
      include 'params.inc'
c
c     Solves the HFB equations as matrix equations in the
c     space of given s.p. states for nucleon type 'iso'.
c
c
      include 'wavfun.inc'
      include 'potinv.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'constn.inc'
c
      REAL(8), INTENT(IN OUT) :: efermiso
      REAL(8), INTENT(OUT) :: partnumout
      INTEGER, INTENT(IN) :: iso
c
!      REAL(8) :: efermhfb(2)
      REAL(8), ALLOCATABLE :: hmatrix(:,:),gapmatrix(:,:)
      REAL(8), ALLOCATABLE :: fullh(:,:),uvm(:,:),equasi(:),spenat(:)
      REAL(8), ALLOCATABLE :: um(:,:),vm(:,:),phasfac(:),paircut(:)
      INTEGER, ALLOCATABLE :: jssav(:)
      LOGICAL,PARAMETER :: tprint=.TRUE.,tprint0=.TRUE.
      LOGICAL,PARAMETER :: tnatural=.true.,toffdiag=.TRUE.
      INTEGER,PARAMETER :: ifcut=1
c
c--------------------------------------------------------------------
c
c     settings and offsets for given 'iso'
c
      IF(ifcut==2 .AND. .NOT.tphpair) 
     &  STOP 'HFB_EQS: ifcut=2 requires option TPHPAIR=.TRUE.'
      IF(tprint0) WRITE(6,'(a,i3,a,1pg13.5)') 'HFB eqs for iso=',iso,
     &          ', eferm=',efermiso
!      IF(iso==1)
!     & WRITE(6,'(a,5(1pg13.5))') 'iso,pots:',iso,pot(2),potls(2),
!     &      bmass(2),d1b(2)
!      efermhfb(iso) = efermiso
      isoadd = (iso-1)*ngrid
      isnpad = (iso-1)*kjmx
      IF(iso==1) THEN
        nmaxa=nmaxp
        isonpl=0
      ELSE IF(iso==2) THEN
        nmaxa=nmaxn
        isonpl=nmaxp
      ELSE
        STOP 'HFB: ISO must be 1 or 2'
      END IF
      nmaxa2 = 2*nmaxa
      IF(tprint) WRITE(6,'(a,4i5)') 'HFB: iso,nmax=',iso,nmaxa
c
      ALLOCATE(hmatrix(nmaxa,nmaxa),gapmatrix(nmaxa,nmaxa))
      ALLOCATE(um(nmaxa,nmaxa),vm(nmaxa,nmaxa),paircut(nmaxa))
      ALLOCATE(fullh(nmaxa2,nmaxa2),uvm(nmaxa2,nmaxa2))
      ALLOCATE(equasi(nmaxa2),phasfac(nmaxa2),jssav(nmaxa2))
      ALLOCATE(spenat(nmaxa2))
      hmatrix=0D0
      gapmatrix=0D0
      fullh=0D0
c
c     set up pairing cutoff (presently only for fixed energy band)
c
      efadd=-cutfac
      efcut=efermiso+efadd
      IF(ifcut.eq.1 .AND. cutfac<0D0) then
        DO n=1,nmaxa
          paircut(n) = SQRT(1D0/(1D0+
     &                   EXP((spenrg(n+isonpl)-efcut)/(0.1D0*efadd))))
        END DO
        IF(tprint) WRITE(6,'(a,1pg13.5)') 'HFB_matrix: pairing cutoff=',
     &             efadd
      END IF
      IF(tprint) WRITE(6,'(a,100(/2(1pg13.5)))') 's.p.energies,cut:',
     &     (spenrg(n+isonpl),paircut(n),n=1,nmaxa)
c
c     set up Hamiltonian matrix and gap matrix
c
!      do i=1,ngrid
!        sqb(i) = sqrt(bmass(i+isoadd))
!      enddo
!c
!      call derv25(sqb(1),d2b(1),one,ngrid,rstep)
!      do i=1,ngrid
!        d2b(i) = d2b(i)/sqb(i)
!      enddo
c
      call derv15(bmass(isoadd+1),d1b(isoadd+1),one,ngrid,rstep)
!      d1b(isoadd+1)=(16.0*bmass(2+isoadd)-15.0*bmass(1+isoadd)
!     *        -bmass(3+isoadd))/(6.0*rstep**rstep)
!      do i=2,ngrid
!        d1b(isoadd+i) = d1b(isoadd+i)*rinvs(i)
!      enddo

      DO js=1,jsmax(iso)
        jsnpl = js+isnpad
        nmax  = npla(1,jsnpl)
        DO n=1,nmax
          npl    = npla(n+1,jsnpl)
          nplm = npl-1
          nplshl = npl/ngrid+1-isonpl
          jssav(nplshl) = js
          jssav(nplshl+nmaxa) = js
          phasfac(nplshl) = 2D0*zwj(js)/0.15915494D0
          phasfac(nplshl+nmaxa) = 2D0*zwj(js)/0.15915494D0
          call hpsi(wf(npl),aux3,wfs(npl),aux2,1,js,iso,0D0)
          DO n2=1,nmax
            npl2    = npla(n2+1,jsnpl)
            nplm2 = npl2-1
            nplshl2 = npl2/ngrid+1-isonpl
!              WRITE(6,'(6i5)') iso,js,n,n2,nplshl,nplshl2
            acch = 0D0
            accgap = 0D0
!            xnorm=0D0
            DO i=2,ngrid
!              xnorm  = wf(nplm+i)*wf(nplm+i) + xnorm
              acch = wf(nplm2+i)*aux3(i) + acch
              accgap = wf(nplm2+i)*wf(nplm+i)*pairpt(i+isoadd)+accgap
            END DO
            acch=rstep*acch
            accgap=rstep*accgap
            IF(n2==n) THEN
              acch2=SQRT(rstep*SUM(aux3(1:ngrid)**2)-acch**2)
              IF(tprint) WRITE(6,'(a,3i4,4(1pg18.10))') 
     &                   'js,n,nplshl,spe,spgap,var=',
     &              js,n,nplshl,acch,accgap,acch2
              acch=acch-efermiso
            END IF
            IF(toffdiag .OR. n2==n) THEN
              hmatrix(nplshl2,nplshl) = acch
              gapmatrix(nplshl2,nplshl) = accgap
              fullh(nplshl2,nplshl) = acch       
              fullh(nplshl2+nmaxa,nplshl+nmaxa) = -acch       
              fullh(nplshl2,nplshl+nmaxa) = accgap
              fullh(nplshl2+nmaxa,nplshl) = accgap
            END IF
          END DO
        END DO
      END DO
      IF(ifcut.eq.2 .AND. cutfac<0D0) then
        CALL diag_real_symm(hmatrix(:,:),equasi(:),vm(:,:),
     &                      nmaxa,nmaxa,(nmaxa+nmaxa*nmaxa)/2)
            WRITE(6,'(a,100(/5(1pg13.5)))') 'h energies:',
     &           equasi(1:nmaxa)
            paircut(1:nmaxa) = SQRT(1D0/(1D0+
     &            EXP((equasi(1:nmaxa)+efermiso-efcut)/(0.1D0*efadd))))
            WRITE(6,'(a,100(/5(1pg13.5)))') 'cutoff:',
     &           paircut(1:nmaxa)
            DO n1=1,nmaxa;  DO n2=1,nmaxa;
              uvm(n1,n2) = SUM(vm(n1,:)*paircut(1:nmaxa)*vm(n2,:))
            END DO; END DO
            DO n1=1,nmaxa;  DO n2=1,nmaxa;
              vm(n1,n2) = SUM(gapmatrix(n1,:)*uvm(1:nmaxa,n2))
            END DO; END DO
            DO n1=1,nmaxa;  DO n2=1,nmaxa;
              gapmatrix(n1,n2) = SUM(uvm(1:nmaxa,n1)*vm(:,n2))
              fullh(n1,n2+nmaxa) = gapmatrix(n1,n2)
              fullh(n1+nmaxa,n2) = gapmatrix(n1,n2)
            END DO; END DO
      ELSE IF(ifcut.eq.1 .AND. cutfac<0D0 .AND. TPHPAIR) then
        DO js=1,jsmax(iso)
          jsnpl = js+isnpad
          nmax  = npla(1,jsnpl)
          DO n=1,nmax
            npl    = npla(n+1,jsnpl)
            nplshl = npl/ngrid+1-isonpl
            spe = fullh(nplshl,nplshl)+efermiso
            paircut(nplshl) = SQRT(1D0/(1D0+
     &                   EXP((spe-efcut)/(0.1D0*efadd))))
          END DO
          DO n=1,nmax
            npl    = npla(n+1,jsnpl)
            nplshl = npl/ngrid+1-isonpl
            DO n2=1,nmax
              npl2    = npla(n2+1,jsnpl)
              nplshl2 = npl2/ngrid+1-isonpl
              IF(nplshl>nmaxa .OR. nplshl2>nmaxa) THEN
                WRITE(6,'(5i4)') js,n,n2,nplshl,nplshl2
                STOP 'out of range'
              END IF
              accgap = gapmatrix(nplshl2,nplshl)
              IF(TPHPAIR)
     &         accgap = accgap*paircut(nplshl)*paircut(nplshl2)
              gapmatrix(nplshl2,nplshl) = accgap
              fullh(nplshl2,nplshl+nmaxa) = accgap
              fullh(nplshl2+nmaxa,nplshl) = accgap
            END DO
          END DO
        END DO
      END IF
c
c     diagonalize
c
      IF(tprint) THEN
        WRITE(6,'(a,i3,a,i3/a)') 'HFB for iso=',iso,
     &   ', nmax=',nmaxa,'h-matrix:'
        DO n1=1,nmaxa
          WRITE(6,'(30(1pg13.5))') (fullh(n2,n1),n2=1,nmaxa)
        END DO
        WRITE(6,'(a)') 'gap-matrix:'
        DO n1=nmaxa+1,nmaxa2
          WRITE(6,'(30(1pg13.5))') (fullh(n2,n1),n2=1,nmaxa)
        END DO
!        WRITE(6,'(a)') 'gap2-matrix:'
!        DO n1=1,nmaxa
!          WRITE(6,'(30(1pg13.5))') (fullh(n2,n1),n2=nmaxa+1,nmaxa2)
!        END DO
!        WRITE(6,'(a)') 'h2-matrix:'
!        DO n1=nmaxa+1,nmaxa2
!          WRITE(6,'(30(1pg13.5))') (fullh(n2,n1),n2=nmaxa+1,nmaxa2)
!        END DO
      END IF
      CALL diag_real_symm(fullh(:,:),equasi(:),uvm(:,:),
     &                    nmaxa2,nmaxa2,(nmaxa2+nmaxa2*nmaxa2)/2)
!      WRITE(6,'(a)') 'V,U matrix along diagonal:'
!      offdigv=0D0
!      offdigu=0D0
!      DO nq=1,nmaxa2
!        IF(equasi(nq)<0D0) THEN
!          umax=0D0
!          vmax=0D0
!          DO n2=1,nmaxa
!            vact = uvm(n2,nq)
!            uact = uvm(n2+nmaxa,nq)
!            IF(ABS(vact)>ABS(vmax)) THEN
!              vmax=vact
!              umax=uact
!              nmax=n2
!            END IF
!          END DO 
!          DO n2=1,nmaxa
!            IF(n2.NE.nmax) THEN
!              offdigv = MAX(offdigv,ABS(uvm(n2,nq)))
!              offdigu = MAX(offdigu,ABS(uvm(n2+nmaxa,nq)))
!            END IF
!          END DO
!          WRITE(6,'(a,i3,4(1pg13.5))') 'nq,v,u=',n1,vmax,umax
!        END IF
!      END DO
!      WRITE(6,'(a,2(1pg13.5))') 'UVM off diagonality=',offdigv,offfdigu
!      WRITE(6,'(a)') 'U matrix:'
!      DO nq=1,nmaxa2
!        IF(equasi(nq)<0D0)
!     &   WRITE(6,'(i3,30(1pg13.5))') nq,(uvm(n2,nq),n2=nmaxa+1,nmaxa2)
!      END DO
      accn=0D0
      DO nq=1,nmaxa2
        IF(ABS(equasi(nq))<1D-30) equasi(nq)=0D0
        accv2=0D0
        DO na=1,nmaxa
          IF(equasi(nq)<0D0) THEN
            accn = uvm(na,nq)**2*phasfac(na) + accn
            accv2 = uvm(na,nq)**2 + accv2
          END IF
        END DO
        DO n1=1,nmaxa
          IF(nq.LE.nmaxa) 
     &         vm(n1,nq) = SUM(fullh(n1,1:nmaxa)*uvm(1:nmaxa,nq))
        END DO
        IF(tprint .AND. nq.LE.nmaxa) 
     &         WRITE(6,'(a,i5,3(1pg13.5))') 'nq,v2,e,spe=',
     &          nq,accv2,equasi(nq),
     &          SUM(vm(:,nq)*uvm(1:nmaxa,nq))/accv2
      END DO
      IF(tprint0) THEN
        WRITE(6,'(a,1pg13.5/a)') ' particle number=',accn,
     &    'quasiparticle energies:'
        IF(tprint) WRITE(6,'(10(5pg13.5))') equasi(1:nmaxa2)
      END IF
      partnumout = accn
c
c     accumulate new densities, use 'aux1' and 'aux2' for new w.f. 
c     and its derivative
c
      IF(.NOT.tnatural) 
     &  CALL  hfb_densit(iso,nmaxa2,uvm,equasi,jssav,phasfac)
c
c     compute natural orbitals, use 'gapmatrix' for density matrix
c
      IF(tnatural) THEN
        gapmatrix=0D0
        DO nq=1,nmaxa2
          IF(equasi(nq)<0D0) THEN
            DO na=1,nmaxa; DO nb=1,nmaxa
              IF(jssav(na)==jssav(nb)) 
     &            gapmatrix(na,nb) = gapmatrix(na,nb) +
     &                             uvm(na,nq)*uvm(nb,nq) !*phasfac(na)
            END DO; END DO
          END IF
        END DO
        IF(tprint0) THEN
          WRITE(6,'(a)') 'rho matrix:'
          accv2=0D0
          DO na=1,nmaxa
            WRITE(6,'(10(1pg13.5))') (gapmatrix(nb,na),nb=1,nmaxa)
            accv2=accv2+gapmatrix(na,na)*phasfac(na)
          END DO
          CALL off_diagonality(gapmatrix,nmaxa,nmaxa,'rho')
          WRITE(6,'(a,1pg13.5)') 'part.numb from diagonal:',accv2
        END IF
        CALL diag_real_symm(gapmatrix(:,:),equasi(:),vm(:,:),
     &                      nmaxa,nmaxa,(nmaxa+nmaxa*nmaxa)/2)
        IF(tprint) THEN
          WRITE(6,'(a)') 'nat.orb. v-matrix:'
          DO n1=1,nmaxa
            WRITE(6,'(30(1pg13.5))') (vm(n2,n1),n2=1,nmaxa)
          END DO
        END IF       
        IF(tprint0) WRITE(6,'(a/a)') 'natural orbitals:','v2:'
        accv2=0D0
        DO na=1,nmaxa
!          weight(na+isonpl) = equasi(n)
          v2max=0D0
          DO n1=1,nmaxa
           v2act = vm(n1,na)**2
           IF(v2act>v2max) THEN
             v2max=v2act
             nplnq=n1
           END IF
          END DO
          accv2 = accv2 + equasi(na)*phasfac(nplnq)
          acceps=efermiso
          DO n1=1,nmaxa; DO n2=1,nmaxa
             acceps = acceps + vm(n1,na)*hmatrix(n1,n2)*vm(n2,na)
          END DO; END DO
          spenat(na) = acceps
          IF(tprint0 .AND. tphpair)
     &       WRITE(6,'(i3,2(1pg13.5))') na,equasi(na),acceps !phasfac(na)
          IF(tprint0 .AND. .NOT.tphpair)
     &       WRITE(6,'(i3,2(1pg13.5))') na,equasi(na)*paircut(na),
     &          acceps*paircut(na) 
        END DO
        IF(tprint0) WRITE(6,'(a,1pg13.5)') 'final part.num.=',accv2
c
c     recover rho and modify occupations
c     use hmatrix and gapmatrx as workspace
c
!        WRITE(6,'(a)') 'test back trafo rho:'
!        WRITE(6,'(a)') 'occups:'
!        WRITE(6,'(30(1pg13.5))') (equasi(n2),n2=1,nmaxa)
!        hmatrix(1:nmaxa,1:nmaxa) = uvm(1:nmaxa,1:nmaxa)
!        CALL  matin2 (hmatrix,nmaxa,nmaxa,nerror,determ)
!        WRITE(6,'(a,i5,1pg13.5)') 'MATIN2: nerror,determ=',
!     &   nerror,determ
!        WRITE(6,'(a)') 'test invert:'
!        DO n1=1,nmaxa; DO n2=1,nmaxa
!          gapmatrix(n1,n2) =SUM(hmatrix(n1,1:nmaxa)*uvm(1:nmaxa,n2))
!        END DO; END DO
!        CALL  off_diagonality(gapmatrix,nmaxa,nmaxa,
!     &                      'invert')
!        WRITE(6,'(a)') 'diagonal elements invert unit:'
!        WRITE(6,'(30(1pg13.5))') (gapmatrix(n2,n2),n2=1,nmaxa)
!        DO nq=1,nmaxa; DO ng=1,nmaxa
!           gapmatrix(nq,ng) =SUM(hmatrix(nq,1:nmaxa)*vm(1:nmaxa,ng))
!        END DO; END DO
!        WRITE(6,'(a)') 'unitary combination:'
!        DO n1=1,nmaxa
!          WRITE(6,'(30(1pg13.5))') (gapmatrix(n2,n1),n2=1,nmaxa)
!        END DO
!        DO nq1=1,nmaxa; DO nq2=1,nmaxa
!           hmatrix(nq1,nq2) =SUM(gapmatrix(nq1,1:nmaxa)*equasi(1:nmaxa)
!     *                         *gapmatrix(nq2,1:nmaxa))
!        END DO; END DO
!        CALL  off_diagonality(hmatrix,nmaxa,nmaxa,'off.diag unit matrix')
!        WRITE(6,'(a)') 'off.diag unit matrix:'
!        DO n1=1,nmaxa
!          WRITE(6,'(30(1pg13.5))') (hmatrix(n2,n1),n2=1,nmaxa)
!        END DO
!      WRITE(6,'(a)') 'diagonal elements unit matrix:'
!      WRITE(6,'(30(1pg13.5))') (hmatrix(n2,n2),n2=1,nmaxa)
!      WRITE(6,'(a)') 'test forward uvm:'
!      DO n1=1,nmaxa; DO n2=1,nmaxa
!         gapmatrix(n1,n2) =SUM(uvm(1:nmaxa2,n1)*uvm(1:nmaxa2,n2))
!      END DO; END DO
!      CALL  off_diagonality(gapmatrix,nmaxa,nmaxa,
!     &                      'off.diag forward unit')
!      WRITE(6,'(a)') 'diagonal elements forward unit:'
!      WRITE(6,'(30(1pg13.5))') (gapmatrix(n2,n2),n2=1,nmaxa)
c
c     define new 'uvm' from natural orbitals
c
        DO n=1,nmaxa
          spe = spenat(n)+efermiso
          paircut(n) = SQRT(1D0/(1D0+EXP((spe-efcut)/(0.1D0*efadd))))
          WRITE(6,'(a,i4,5(1pg13.5))') '# test nat: n,spe,paircut=',
     &       n,spe,paircut(n),equasi(n),SQRT((1D0-equasi(n))*paircut(n))
          uvm(1:nmaxa,n) = vm(1:nmaxa,n)*
     &        SQRT(MAX(equasi(n),1D-40)*paircut(n))
          uvm(nmaxa+1:nmaxa2,n) = vm(1:nmaxa,n)*
     &         SQRT(MAX(1D0-equasi(n),1D-40)*paircut(n))
        END DO
        spenat(1:nmaxa)=-1D0;spenat(nmaxa+1:nmaxa2)=1D0;
        CALL  hfb_densit(iso,nmaxa2,uvm,spenat,jssav,phasfac)
c
      END IF
c
      DEALLOCATE(um,vm,hmatrix,gapmatrix,fullh,uvm,equasi,phasfac)
      DEALLOCATE(spenat)
c
      return
      end
c-----hfb_densit---------------------------------------------------------
c
      subroutine hfb_densit(iso,nmaxa2,uvm,equasi,jssav,phasfac)
c
      include 'params.inc'
c
c     Accumulates densities. Version for the HFB package.
c
      include 'wavfun.inc'
      include 'potinv.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'constn.inc'
c
      INTEGER, INTENT(IN) :: nmaxa2,iso,jssav(nmaxa2)
      REAL(8), INTENT(IN) :: uvm(nmaxa2,nmaxa2),equasi(nmaxa2)
      REAL(8), INTENT(IN) :: phasfac(nmaxa2)

      LOGICAL,PARAMETER :: tprint0=.TRUE.,tprint=.TRUE.
c
c      dimension zwj(21)
c
c-------------------------------------------------------------------------------
c
c
!      nmaxa = nmaxa2
      isoadd = (iso-1)*ngrid
      isnpad = (iso-1)*kjmx
      ie     = ngrid*iso
      im     = ie-ngrid+1
      ia     = im+1
      rho(im:ie)=0D0 
      tau(im:ie)=0D0 
      taul(im:ie)=0D0
      sj(im:ie)=0D0
      paircr(im:ie)=0D0
      ngw    = (iso-1)*nmaxp
      nadd   = ngw*ngrid-ie
      IF(iso==1) THEN
        nmaxa=nmaxp
        isonpl=0
      ELSE IF(iso==2) THEN
        nmaxa=nmaxn
        isonpl=nmaxp
      ELSE
        STOP 'HFB: ISO must be 1 or 2'
      END IF
c
      acctot=0D0
      DO nq=1,nmaxa2
        IF(equasi(nq)<0D0) THEN
          aux1=0D0
          aux2=0D0
          accno=0D0
          v2max=0D0
          DO js=1,jsmax(iso)
            jsnpl = js+isnpad
            ll     = js/2
            xlp    = ll*ll+ll+1D0
            wl     = wls(js)
            nmax  = npla(1,jsnpl)
            sig    = 2D0*MOD(ll,2)-1D0
c
c     new wave functions from HFB coefficients
c
            DO n=1,nmax
              npl    = npla(n+1,jsnpl)
              nplm = npl-1
              nplshl = npl/ngrid+1-isonpl
              IF(jssav(nplshl).NE.js) STOP 'HFB_deriv: js not matching'
              aux1(1:ngrid) = uvm(nplshl,nq)*wf(npl:nplm+ngrid)
     &                         +aux1(1:ngrid)
              aux2(1:ngrid) = (uvm(nplshl+nmaxa,nq)
     &                            *wf(npl:nplm+ngrid))
     &                         +aux2(1:ngrid)
              v2add=uvm(nplshl,nq)**2
!              WRITE(6,'(a,3i3,1pg13.5)') 'nq,js,n,v2add=',nq,js,n,v2add
              accno = v2add + accno
              acctot = v2add*phasfac(nplshl) + acctot
              IF(v2add>v2max) THEN
                jsnq=js
                v2max=v2add
              END IF
            END DO
          END DO
!          WRITE(6,'(a,(/5(1pg13.5)))') 'ratio U/V:',
!     &        aux2(1:ngrid)/aux1(1:ngrid)
          IF(tprint0)
     &     WRITE(6,'(a,2i4,3(1pg13.5))') 'nq,js,norms=',nq,jsnq,
     &          rstep*2D0*zwj(js)/0.15915494D0*SUM(aux1*aux1),
     &          rstep*2D0*zwj(js)/0.15915494D0*SUM(aux1*aux2),
     &          accno
c
          js=jsnq
          ll     = js/2
          xlp    = ll*ll+ll+1D0
          wl     = wls(js)
          sig    = 2D0*MOD(ll,2)-1D0
          CALL derv15(aux1,aux3,sig,ngrid,rstep)
          om     = zwj(js)
          ompair = om/2D0
!              nadd   = nadd+ngrid
c                                    the case r=0.0 needs special treatment
          if(ll.lt.1) then
            x      = aux1(1)
            rho(im) = (x*x)*om + rho(im)
            paircr(im) = (x*aux2(1))*om + paircr(im)
          elseif(ll.eq.1) then
            x      = (16.0*aux1(2)-aux1(3))/(6.0*rstep*rstep)
            tau(im) =0.75*om*(x*x) + tau(im)
          endif
c                                    now all r>0.0
          r      = 0D0      
          DO i=ia,ie
            r      = r+rstep
            ri     = one/r
            phas   = (ri*ri)*om
c  
            wfa    = aux1(i-isoadd)
            wfa2   = wfa*wfa
c         if(wfa2.gt.xsmall) then
            rho(i) = (phas*wfa2) + rho(i)
            taul(i) = (phas*wfa2)*xlp + taul(i)
            sj(i)  = (phas*wfa2)*wl*ri + sj(i)
c  
            wfas   = aux3(i-isoadd)
            tau(i)=(wfas*wfas-((wfa*wfas)+(wfa*wfas)-xlp*wfa2*ri)*ri)
     &               *phas  + tau(i)
            paircr(i) = ompair*(ri*ri)*wfa*aux2(i-isoadd) + paircr(i)
          END DO
        END IF
      END DO
      IF(tprint)
     & WRITE(6,'(a,i4,1pg13.5)') 'HFB_densit test: iso,part.num',
     &    iso,acctot
c
c     some derivatives: d/dr 'rho' --> 'd1r'
c                       d/dr**2 'rho' + 2/r d/dr 'rho' --> 'd2r'
c                       d/dr 'sj' + 2/r 'sj' --> 'sjs'
c     'im' is the initial position for given 'iso'.
c
        zwh    = 2D0/rstep
        call derv15(rho(im),d1r(im), one,ngrid,rstep)
        call derv25(rho(im),d2r(im), one,ngrid,rstep)
        d2r(im) = d2r(im) + (8.0*d1r(im+1)-d1r(im+2))/(3.0*rstep)
        do i=im+1,ie
          d2r(i) = d2r(i) + zwh*d1r(i)/(i-im)
        enddo
        if(t9.ne.zero) then
         call derv15(tau(im),d1t(im), one,ngrid,rstep)
         call derv25(tau(im),d2t(im), one,ngrid,rstep)
         d2t(im) = d2t(im) + (8.0*d1t(im+1)-d1t(im+2))/(3.0*rstep)
         do i=im+1,ie
           d2t(i) = d2t(i) + zwh*d1t(i)/(i-im)
         enddo
        endif
c  
        call derv15(sj(im),sjs(im),-one,ngrid,rstep)
        sjs(im) = 3.0*sjs(im)
        do i=im+1,ie
          sjs(i) = sjs(i) + zwh*sj(i)/(i-im)
        enddo
c
c     test density
c
      IF(tprint)  THEN
        WRITE(6,'(a,4i3)') 'HFB densit: iso=',iso,im,ia,ie
        WRITE(6,'(a)') '# r  rho  tau  sj paircr  d2r  sjs'
        DO i=im,ie
          WRITE(6,'(f6.1,6(1pg13.5))') 
     &      (i-im)*rstep,rho(i),tau(i),sj(i),paircr(i),d2r(i),sjs(i)
        END DO
      END IF
!      WRITE(6,'(a,100(/3(1pg13.5)))') '#new  rho_p  rho_n chi_n:',
!     &  (rho(i),rho(i+ngrid),paircr(i+ngrid),i=1,ngrid)
!      IF(iso==2)  WRITE(6,'(a,100(/5(1pg13.5)))') '#new paircr_n:',
!     &          (paircr(i),i=ngrid+1,2*ngrid)
c
      acc=0D0
      accp=0D0
      r=0D0
      DO i=ia,ie
        r       = rstep+r
        r2      = r*r
        acc  = r2*rho(i) +acc
        accp = r2*paircr(i) +accp
      END DO
      acc = frpirs*acc
      accp = frpirs*accp
      IF(tprint0)
     & WRITE(6,'(a,i2,2(1pg13.5))') 'HFB densit: iso,norm=',
     &   iso,acc,accp
c

      RETURN
      END
c-----off_diagonality-------------------------------------------------
c
      SUBROUTINE off_diagonality(rmatrix,nmax,ndim,atext)
c
      INTEGER,INTENT(IN) :: nmax,ndim
      REAL(8),INTENT(IN) :: rmatrix(ndim,ndim)
      CHARACTER(*) :: atext
c
c---------------------------------------------------------------------
c
!      WRITE(6,'(a)') 'state-wise off-diagonality for matrix '//atext
      DO n=1,nmax
        add = SUM(rmatrix(n,:)**2)-rmatrix(n,n)**2
        acc = acc + add
!        WRITE(6,'(i3,1pg13.5)') n,SQRT(add)
      END DO
      WRITE(6,'(a,1pg13.5)') 'average off-diagonality for: '//atext,
     &     SQRT(acc/nmax)
      RETURN
      END
c-----hfb_matrix------------------------------------------------------
c
      subroutine hfb_matrix2()
c
      include 'params.inc'
c
c     Solves the two-components HFB equations for quasi-particles
c     as matrix equations in the space of given s.p. states. 
c
c
      include 'wavfun.inc'
      include 'potinv.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'constn.inc'
!      include 'result.inc'
c
      REAL(8), ALLOCATABLE :: hmatrix(:,:,:),gapmatrix(:,:,:)
      REAL(8), ALLOCATABLE :: fullh(:,:,:),uvm(:,:,:),equasi(:,:)
      REAL(8), ALLOCATABLE :: um(:,:,:),vm(:,:,:),phasfac(:,:)
      INTEGER, ALLOCATABLE :: ncount(:)
      REAL(8) :: efermhfb(2)
c
c--------------------------------------------------------------------
c
      nmaxm=MAX(nmaxp,nmaxn)
      ALLOCATE(hmatrix(nmaxm,nmaxm,2),gapmatrix(nmaxm,nmaxm,2))
      ALLOCATE(um(nmaxm,nmaxm,2),vm(nmaxm,nmaxm,2),ncount(nmaxm))
      nmaxm2=2*nmaxm
      ALLOCATE(fullh(nmaxm2,nmaxm2,2),uvm(nmaxm2,nmaxm2,2))
      ALLOCATE(equasi(nmaxm2,2),phasfac(nmaxm2,2))
      hmatrix=0D0
      gapmatrix=0D0
      ncount=0
      efermhfb(1) = efep
      efermhfb(2) = efen
c
c     set up Hamiltonian matrix and gap matrix
c
      DO iso=1,2
        efermiso=efermhfb(iso)
        WRITE(6,'(a,i3,a,1pg13.5,a,i3)') 'HFB for iso=',iso,
     &          ', eferm=',efermiso,', jsmax=',jsmax(iso)
        isoadd = (iso-1)*ngrid
        isnpad = (iso-1)*kjmx
        IF(iso==1) THEN
          nmaxa=nmaxp
          isonpl=0
        ELSE IF(iso==2) THEN
          nmaxa=nmaxn
          isonpl=nmaxp
        END IF
        nmaxa2 = 2*nmaxa
        IF(nmaxa2>nmaxm2) STOP 'HFB: not enough space'
      do i=1,ngrid
        sqb(i) = sqrt(bmass(i+isoadd))
      enddo
c
      call derv25(sqb(1),d2b(1),one,ngrid,rstep)
      do i=1,ngrid
        d2b(i) = d2b(i)/sqb(i)
      enddo
c
      call derv15(bmass(isoadd+1),d1b(1),one,ngrid,rstep)
      d1b(1)=(16.0*bmass(2+isoadd)-15.0*bmass(1+isoadd)
     *        -bmass(3+isoadd))/(6.0*rstep**rstep)
      do i=2,ngrid
        d1b(i) = d1b(i)*rinvs(i)
      enddo

        DO js=1,jsmax(iso)
          jsnpl = js+isnpad
          nmax  = npla(1,jsnpl)
          DO n=1,nmax
            npl    = npla(n+1,jsnpl)
            nplm = npl-1
            nplshl = npl/ngrid+1 -isonpl
!            ncount(nplshl) = nplshl
            WRITE(6,'(6i5)') iso,jsmax(iso),js,nmax,n,nplshl
            phasfac(nplshl,iso) = 2D0*zwj(js)/0.15915494D0
            phasfac(nplshl+nmaxa,iso) = 2D0*zwj(js)/0.15915494D0
            call hpsi(wf(npl),aux3,wfs(npl),aux2,1,js,iso,0D0)
            DO n2=1,nmax
              npl2    = npla(n2+1,jsnpl)
              nplm2 = npl2-1
              nplshl2 = npl2/ngrid+1 -isonpl
!              WRITE(6,'(6i5)') iso,js,n,n2,nplshl,nplshl2
              acch = 0D0
              accgap = 0D0
!              xnorm=0D0
              DO i=2,ngrid
!                xnorm  = wf(nplm+i)*wf(nplm+i) + xnorm
                acch = wf(nplm2+i)*aux3(i) + acch
                accgap = wf(nplm2+i)*wf(nplm+i)*pairpt(i+isoadd)+accgap
              END DO
              acch=0.3D0*acch
              accgap=0.3D0*accgap
              IF(n2==n) THEN
                acch=acch-efermiso
                WRITE(6,'(a,2i4,3(1pg18.10))') 'js,n,spe=',js,n,acch,
     &                       accgap
              END IF
              hmatrix(nplshl2,nplshl,iso) = acch
              gapmatrix(nplshl2,nplshl,iso) = accgap
              fullh(nplshl2,nplshl,iso) = acch       
              fullh(nplshl2+nmaxa,nplshl+nmaxa,iso) = -acch       
              fullh(nplshl2,nplshl+nmaxa,iso) = accgap
              fullh(nplshl2+nmaxa,nplshl,iso) = accgap
            END DO
          END DO
        END DO
!        WRITE(6,*) ncount
c
c     diagonalize
c
        WRITE(6,'(a,i3,a,i3/a)') 'HFB for iso=',iso,', nmax=',nmaxa,
     &   'h-matrix:'
        DO n1=1,nmaxa
          WRITE(6,'(30(1pg13.5))') (fullh(n2,n1,iso),n2=1,nmaxa)
        END DO
        WRITE(6,'(a)') 'h-matrix_2:'
        DO n1=nmaxa+1,nmaxa2
          WRITE(6,'(30(1pg13.5))') (fullh(n2,n1,iso),n2=nmaxa+1,nmaxa2)
        END DO
        WRITE(6,'(a)') 'gap-matrix:'
        DO n1=nmaxa+1,nmaxa2
          WRITE(6,'(30(1pg13.5))') (fullh(n2,n1,iso),n2=1,nmaxa)
        END DO
        WRITE(6,'(a)') 'gap-matrix_2:'
        DO n1=1,nmaxa
          WRITE(6,'(30(1pg13.5))') (fullh(n2,n1,iso),n2=nmaxa+1,nmaxa2)
        END DO
        CALL diag_real_symm(fullh(:,:,iso),equasi(:,iso),uvm(:,:,iso),
     &                    nmaxa2,nmaxm2,(nmaxm2+nmaxm2*nmaxm2)/2)
        DO nq=1,nmaxa2
          IF(ABS(equasi(nq,iso))<1D-30) equasi(nq,iso)=0D0
          WRITE(6,'(a,i5,1pg13.5)') ' nq,equasi=',nq,equasi(nq,iso)
        END DO
        accn=0D0
        DO na=1,nmaxa    ! canonical
          accv2=0D0
          accuv=0D0
          DO nq=1,nmaxa   ! qp states
            nq2=nq+nmaxa
!            IF(equasi(nq,iso)<0D0) THEN
              accn = uvm(na,nq,iso)**2*phasfac(na,iso) + accn
              accuv = uvm(na,nq2,iso)**2 + accuv
              accv2 = uvm(na,nq,iso)**2 + accv2
!            END IF
          END DO
          WRITE(6,'(a,i5,3(1pg13.5))')'na,v2=',na,accv2,accuv,
     &      fullh(na,na,iso) !,spenrg(na)
        END DO
        WRITE(6,'(a,1pg13.5/a)') ' particle number=',accn
!        DO nq=1,nmaxa2
!          accv2=0D0
!          DO na=1,nmaxa
!            IF(ABS(equasi(nq,iso))<1D-30) equasi(nq,iso)=0D0
!            IF(equasi(nq,iso)<0D0) THEN
!              accn = uvm(na,nq,iso)**2*phasfac(na,iso) + accn
!              accv2 = uvm(na,nq,iso)**2 + accv2
!            END IF
!          END DO
!          WRITE(6,'(a,i5,2(1pg13.5))') 'nq,v2=',nq,accv2,equasi(nq,iso)
!        END DO
!        WRITE(6,'(a,1pg13.5/a)') ' particle number=',accn,
!     &    'quasiparticle energies:'
!        WRITE(6,'(10(5pg13.5))') equasi(1:nmaxa2,iso)
c
      END DO
c
!      STOP  "HFB"
      DEALLOCATE(um,vm,hmatrix,gapmatrix,fullh,uvm,equasi,phasfac,
     &           ncount)
      return
      end

      include 'hfb-natural.f'
