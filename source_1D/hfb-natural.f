c-----hfb_eqs---------------------------------------------------------
c
      subroutine hfb_natural(iso,efermiso,partnum,iternata)
c
      include 'params.inc'
c
c     Solves the HFB equations for natural orbitals as matrix equations 
c     in the space of given s.p. states for nucleon type 'iso'.
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
      REAL(8), INTENT(IN OUT) :: partnum
      INTEGER, INTENT(IN) :: iso,iternata
c
!      REAL(8) :: efermhfb(2)
      REAL(8), ALLOCATABLE :: cm(:,:),cmnew(:,:),symcond(:,:)
      REAL(8), ALLOCATABLE :: hcm(:,:)
      REAL(8), ALLOCATABLE :: hmatrix(:,:),gapmatrix(:,:)
      REAL(8), ALLOCATABLE :: phasfac(:),paircut(:),spenat(:),gapnat(:)
      REAL(8), ALLOCATABLE :: gwa(:),v2a(:),uva(:),lambda(:,:),dummy(:)
      INTEGER, ALLOCATABLE :: jssav(:)
      LOGICAL,PARAMETER :: tprint=.FALSE.,tprint0=.FALSE.
      LOGICAL,PARAMETER :: tprint1=.FALSE.
      LOGICAL,PARAMETER :: toffdiag=.TRUE.,tnewv2=.FALSE.
      INTEGER,PARAMETER :: ifcut=1
!      INTEGER,PARAMETER :: iternata=0!100
      REAL(8),PARAMETER :: delstepbas=3.0D0,e0bas=10D0
c
c--------------------------------------------------------------------
c
      IF(iternata<0) RETURN
c
      isoadd = (iso-1)*ngrid
      isnpad = (iso-1)*kjmx
      IF(iso==1) THEN
        nmaxa=nmaxp
        isonpl=0
        efadd=efaddp
        efermiso=efep
      ELSE IF(iso==2) THEN
        nmaxa=nmaxn
        isonpl=nmaxp
        efadd=efaddn
        efermiso=efen
      ELSE
        STOP 'HFB: ISO must be 1 or 2'
      END IF
      IF(tprint0) WRITE(6,'(a,2i4,4(1pg13.5))') 
     &  'HFB natural for iso,nmaxa,eferm,partnum=',iso,nmaxa,efermiso,
     &   partnum
      ALLOCATE(lambda(nmaxa,nmaxa),cmnew(nmaxa,nmaxa),gapnat(nmaxa))
      ALLOCATE(hmatrix(nmaxa,nmaxa),gapmatrix(nmaxa,nmaxa))
      ALLOCATE(gwa(nmaxa),hcm(nmaxa,nmaxa))
      ALLOCATE(paircut(nmaxa),phasfac(nmaxa),spenat(nmaxa))
      ALLOCATE(jssav(nmaxa),v2a(nmaxa),uva(nmaxa),cm(nmaxa,nmaxa))
      ALLOCATE(dummy(nmaxa),symcond(nmaxa,nmaxa))
c
c     copy old weights to new workspace 
c      
      DO js=1,jsmax(iso)
        jsnpl = js+isnpad
        nmax  = npla(1,jsnpl)
        DO n2=1,nmax
          nplq    = npla(n2+1,jsnpl)
          nq = nplq/ngrid+1-isonpl
          nqi=nq+isonpl
!          WRITE(6,'(3i5,4(1pg13.5))') 
!     &     js,n2,nq,weight(nqi),v2a(nq),phpair(nqi),paircut(nq)
          v2a(nq)=weight(nqi)
          paircut(nq)=phpair(nqi)
        END DO
      END DO
      gwa=v2a
c
c     set up Hamiltonian matrix, gap matrix, and initial value for
c     unitary transformation to natural orbits
c
      call derv15(bmass(isoadd+1),d1b(isoadd+1),one,ngrid,rstep)

      hmatrix=0D0
      gapmatrix=0D0
      cm=0D0
      DO js=1,jsmax(iso)
        jsnpl = js+isnpad
        nmax  = npla(1,jsnpl)
        DO n=1,nmax
          npl    = npla(n+1,jsnpl)
          nplm = npl-1
          nplshl = npl/ngrid+1-isonpl
          jssav(nplshl) = js
          phasfac(nplshl) = zwj(js)/0.15915494D0  ! without spin degeneracy
          call hpsi(wf(npl),aux3,wfs(npl),aux2,1,js,iso,0D0)
          DO n2=1,nmax
            npl2    = npla(n2+1,jsnpl)
            nplm2 = npl2-1
            nplshl2 = npl2/ngrid+1-isonpl
!              WRITE(6,'(6i5)') iso,js,n,n2,nplshl,nplshl2
            acch = 0D0
            accgap = 0D0
!            xnorm=0D0
            IF(thfbfold) THEN
              DO i=2,ngrid
                acch = wf(nplm2+i)*aux3(i) + acch
                accgap = wffold(nplm2+i)*wffold(nplm+i)*pairpt(i+isoadd)
     &                  +accgap
              END DO
            ELSE
              DO i=2,ngrid
                acch = wf(nplm2+i)*aux3(i) + acch
                accgap = wf(nplm2+i)*wf(nplm+i)*pairpt(i+isoadd)+accgap
              END DO
            END IF
            acch=rstep*acch
            accgap=rstep*accgap
            IF(n2==n) THEN
              acch2=SQRT(rstep*SUM(aux3(1:ngrid)**2)-acch**2)
              IF(tprint) WRITE(6,'(a,3i4,4(1pg18.10))') 
     &                   'js,n,nplshl,spe,spgap,var=',
     &              js,n,nplshl,acch,accgap,acch2
              acch=acch  !-efermiso
            END IF
            IF(toffdiag .OR. n2==n) THEN
              hmatrix(nplshl2,nplshl) = acch
              gapmatrix(nplshl2,nplshl) = accgap
            END IF
          END DO
          cm(nplshl,nplshl) = 1D0
          spenat(nplshl) = hmatrix(nplshl,nplshl)  !+efermiso
          gapnat(nplshl) = -gapmatrix(nplshl,nplshl)
!          paircut(nplshl) = (1D0/(1D0+
!     &                     EXP((spenat(nplshl)-efcut)/(0.1D0*efadd))))
        END DO
      END DO
      IF(tnewv2) THEN
!        IF(iso==1) WRITE(6,'(a,200(1pg13.5))') 'v2a:   ',v2a(1:nmaxa)
!        CALL pairdn(spenat,gapnat,gwa,phasfac,paircut,dummy,
!     &              v2a,
!     &              cutfac,partnum,efermiso,epairact,deltpact,deltp2,
!     &              eqminp,gpp,efadd,0D0,nmaxa,nmaxa,1,0,
!     &              50,epspr,-1)
!        IF(ecstab.NE.0)
!         CALL pair_stab(spenat,gapnat,gwa,phasfac,paircut,v2a,
!     &                   epairstab,0D0,nmaxa,nmaxa,1)
!        IF(tprint0) 
!     &    write(6,'(a,2(1pg13.4))') 'epair after STAB=',epairstab
        IF(tprint1 .AND. iso==1) THEN
          WRITE(6,'(a,200(1pg13.5))') 'v2a:   ',v2a(1:nmaxa)
          WRITE(6,'(a,200(1pg13.5))') 'gapnat:',gapnat(1:nmaxa)
        END IF
        IF(ecstab>0)
     &    CALL pair_stab(spenat,gapnat,gwa,phasfac,paircut,v2a,
     &                   epairstab,ecstab,nmaxa,nmaxa,1)
        IF(tprint0) 
     &        write(6,'(a,1pg13.4)') 'epair=',epairstab
        IF(tprint1 .AND. iso==1) THEN
          WRITE(6,'(a,200(1pg13.5))') 'v2a:   ',v2a(1:nmaxa)
          WRITE(6,'(a,200(1pg13.5))') 'gapnat:',gapnat(1:nmaxa)
        END IF
        CALL pairdn(spenat,gapnat,gwa,phasfac,paircut,dummy,
     &              v2a,
     &              cutfac,partnum,efermiso,epairact,deltpact,deltp2,
     &              eqminp,gpp,efadd,0D0,nmaxa,nmaxa,1,0,
     &              50,epspr,-1,tlipnog)
        IF(ecstab>0)
     &    CALL pair_stab(spenat,gapnat,gwa,phasfac,paircut,v2a,
     &                   epairstab,0D0,nmaxa,nmaxa,1)
        IF(tprint0) WRITE(6,'(a,2(1pg13.5))') 'after PAIRDN: epair=',
     &              epairact,epairstab
!     &              10,epspr,-1)
      END IF
!      RETURN
      IF(tphpair) THEN
        uva(:)=SQRT((1D0-v2a(:))*v2a(:))*paircut(:)
      ELSE
        uva(:)=SQRT((1D0-v2a(:))*v2a(:))
      END IF
        avsymcond=0D0
        symcondmax=0D0
        symcond=0D0
        DO na=1,nmaxa; DO nb=1,nmaxa
          IF(jsp(na).NE.jsp(nb)) CYCLE
          symcond(na,nb) = 
     &         +(v2a(nb)*hmatrix(na,nb)-v2a(na)*hmatrix(nb,na))
     &         +(uva(nb)*gapmatrix(na,nb)-uva(na)*gapmatrix(nb,na))
          symcondmax=MAX(symcondmax,ABS(symcond(na,nb)))
          avsymcond=avsymcond+symcond(na,nb)**2
        END DO; END DO
        symcondi(iso)=SQRT(avsymcond/partnum)
        IF(tprint0) THEN      
          WRITE(6,'(a,i5)') 'HFB natural: iter=',0
          WRITE(6,'(a,7(1pg13.5))') ' partnum,efermi,epair,cutfac=',
     &          partnum,SUM(v2a*phasfac)*2D0,efermiso,epair,cutfac
          IF(tprint1) THEN
            WRITE(6,'(a)')'# n  js  spe  v2  uv  cut   phasfac'
            DO na=1,nmaxa
            WRITE(6,'(2i3,5(1pg13.5))') na,jssav(na),spenat(na),v2a(na),
     &            gapnat(na),paircut(na),phasfac(na)
              nb=na+isonpl
              WRITE(6,'(6x,5(1pg13.5))') spenrg(nb),weight(nb),
     &            deltaf(nb),phpair(nb),ph(nb)
            END DO
          END IF
!          WRITE(6,'(a)') 'asymmetry matrix:'
!          DO na=1,nmaxa
!            WRITE(6,'(i3,50(1pg13.5))') na,symcond(:,na)
!          END DO
          WRITE(6,'(a,i2,4(1pg13.5))')
     &       'HFB-natural init.symm. condition: iso,max,average=',
     &        iso,symcondmax,symcondi(iso)
          IF(tprint) THEN
            WRITE(6,'(a)') 'h matrix:'
            DO na=1,nmaxa
              WRITE(6,'(i3,a,500(1pg13.5))') na,': ',hmatrix(:,na)
            END DO
            WRITE(6,'(a)') 'gap matrix:'
            DO na=1,nmaxa
              WRITE(6,'(i3,a,500(1pg13.5))') na,': ',gapmatrix(:,na)
            END DO
          END IF
        END IF

      IF(iternata.LE.0) RETURN
c
c     iterate
c
      iab=-2
      IF(tprint1) iab=-1
      DO iter=1,iternata

c     h->psi in configuration space
        DO na=1,nmaxa 
          v2act=v2a(na)
          uvact=uva(na)
          DO n1=1,nmaxa
            hcm(n1,na) = v2act*SUM(hmatrix(n1,:)*cm(:,na))+
     &                   uvact*SUM(gapmatrix(n1,:)*cm(:,na))
          END DO
        END DO
c     matrix of Lagrangian parameters and symmetry condition
!        lambda=0D0
!        uva(:)=uva(:)*paircut(:)
        symcond=0D0
        lambda=0D0
        DO na=1,nmaxa; DO nb=1,na
          IF(jsp(na).NE.jsp(nb)) CYCLE
          cbhca = SUM(cm(:,nb)*hcm(:,na))
          cahcb = SUM(cm(:,na)*hcm(:,nb))
          lambda(na,nb)= 0.5D0*(cbhca+cahcb)
          lambda(nb,na)=lambda(na,nb)
          symcond(na,nb)= (cbhca-cahcb)
          symcond(nb,na)= -symcond(na,nb)
!          DO n1=1,nmaxa; DO n2=1,nmaxa
!            lambda(na,nb)= cm(n1,na)*cm(n2,nb)*0.5D0*(
!     &            (v2a(na)+v2a(nb))*hmatrix(n1,n2)
!     &           +(uva(na)+uva(nb))*gapmatrix(n1,n2)  )
!     &           + lambda(na,nb)
!            symcond(na,nb)= cm(n1,na)*cm(n2,nb)*(
!     &            (v2a(na)-v2a(nb))*hmatrix(n1,n2)
!     &           +(uva(na)-uva(nb))*gapmatrix(n1,n2)  )
!     &           + symcond(na,nb)
!          END DO; END DO
        END DO; END DO

!        nitsub=1
!        DO itsub=1,nitsub
c     gradient step
        spemin=MINVAL(spenat,nmaxa)
        spemax=MAXVAL(spenat,nmaxa)
        delstep=delstepbas/(spemax-spemin)
        cmnew=0D0
        DO na=1,nmaxa; DO n1=1,nmaxa
          IF(jsp(na).NE.jsp(n1)) CYCLE      
          IF(e0bas==0D0) THEN
            dampstep = delstep
          ELSE
            dampstep = delstepbas/(hmatrix(n1,n1)-hmatrix(1,1)+e0bas)
 !    &                  /nitsub
          END IF
          cmnew(n1,na) = cm(n1,na) - dampstep*(
     &        v2a(na)*SUM(hmatrix(n1,:)*cm(:,na))
     &       +uva(na)*SUM(gapmatrix(n1,:)*cm(:,na))
     &       -SUM(cm(n1,:)*lambda(:,na))     )
        END DO; END DO

c     ortho-normalize
        cm=cmnew
        DO na=1,nmaxa
          DO nb=1,na-1
            IF(jsp(na).NE.jsp(nb)) CYCLE
            ovlab=SUM(cm(:,nb)*cm(:,na))
            cm(:,na)=cm(:,na)-cm(:,nb)*ovlab
          END DO
          ovlaa=1D0/SQRT(SUM(cm(:,na)*cm(:,na)))
          cm(:,na)=cm(:,na)*ovlaa
        END DO

!        END DO          ! end of sub-iterations
c     new energies and occupations
        DO na=1,nmaxa
          spenat(na) = 0D0
          gapnat(na) = 0D0
          DO n1=1,nmaxa
            IF(jsp(na).NE.jsp(n1)) CYCLE
            spenat(na) =spenat(na)+cm(n1,na)*SUM(hmatrix(n1,:)*cm(:,na))
            gapnat(na) = gapnat(na)+
     &                   cm(n1,na)*SUM(gapmatrix(n1,:)*cm(:,na))
          END DO
        END DO
        IF(ecstab.NE.0)
     &   CALL pair_stab(spenat,gapnat,gwa,phasfac,paircut,v2a,
     &                  epairstab,ecstab,nmaxa,nmaxa,1)
        call pairdn(spenat,gapnat,gwa,phasfac,paircut,dummy,v2a,
     &              cutfac,partnum,efermi,epairact,deltpact,deltp2,
     &              eqminp,gpp,efadd,0D0,nmaxa,nmaxa,1,0,
     &              50,epspr,iab,tlipnog)
        IF(tphpair) THEN
          uva(:)=SQRT((1D0-v2a(:))*v2a(:))*paircut(:)
        ELSE
          uva(:)=SQRT((1D0-v2a(:))*v2a(:))
        END IF
        avsymcond=0D0
        symcondmax=0D0
        DO na=1,nmaxa
          avsymcond=avsymcond+SUM(symcond(:,na)**2)
          symcondmax=MAX(symcondmax,MAXVAL(ABS(symcond(:,na))))
        END DO
        IF(tprint0) THEN      
          WRITE(6,'(a,i5)') 'HFB natural: iter=',iter
          WRITE(6,'(a,6(1pg13.5))') ' partnum,efermi,epair,delstep=',
     &          partnum,efermi,epairact,delstep
          WRITE(6,'(a,i4,4(1pg13.5))')
     &       'symmetry condition: iter,max,average=',
     &        iter,symcondmax,SQRT(avsymcond/partnum)
         END IF
         IF(tprint1) THEN      
          WRITE(6,'(a)')'# n  js  spe  v2  uv  cut   phasfac'
          DO na=1,nmaxa
            WRITE(6,'(2i3,5(1pg13.5))') na,jssav(na),spenat(na),v2a(na),
     &            uva(na),paircut(na),phasfac(na)
          END DO
          IF(tprint) THEN
            WRITE(6,'(a)') 'symmetry condition in detail:'
            DO na=1,nmaxa
              WRITE(6,'(i3,a,500(1pg13.5))') na,': ',symcond(:,na)
            END DO
          END IF
        END IF

      END DO

      symcondf(iso)=SQRT(avsymcond/partnum)
      WRITE(6,'(a,i2,4(1pg13.5))')
     &       'HFB final symmetry condition: iso,max,average=',
     &        iso,symcondmax,symcondf(iso)
!      CALL  hfb_nat_densit(iso,nmaxa,cm,v2a,uva,jssav,phasfac)
      CALL  hfb_nat_mapwf(iso,nmaxa,cm,v2a,paircut,jssav,phasfac)

      DEALLOCATE(lambda,cmnew,gapnat,hmatrix,gapmatrix)
      DEALLOCATE(paircut,phasfac,spenat,jssav,v2a,uva,cm)
      DEALLOCATE(dummy,symcond)

      IF(iso==1) THEN
        efaddp=efadd
        efep=efermiso
      ELSE
        efaddn=efadd
        efen=efermiso
      END IF

!      IF(iso==2) STOP

      RETURN
      END
c-----hfb_nat_densit---------------------------------------------------------
c
      subroutine hfb_nat_densit(iso,nmaxa,cm,v2a,uva,jssav,phasfac)
c
      include 'params.inc'
c
c     Accumulates densities. Version for the HFB-natural package.
c
      include 'wavfun.inc'
      include 'potinv.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'constn.inc'
c
      INTEGER, INTENT(IN) :: nmaxa,iso,jssav(nmaxa)
      REAL(8), INTENT(IN) :: cm(nmaxa,nmaxa),v2a(nmaxa),uva(nmaxa)
      REAL(8), INTENT(IN) :: phasfac(nmaxa)

      LOGICAL,PARAMETER :: tprint0=.FALSE.,tprint=.FALSE.
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
        IF(nmaxa.NE.nmaxp) STOP 'HFB_nat_densit: wrong # of states'
        isonpl=0
      ELSE IF(iso==2) THEN
        IF(nmaxa.NE.nmaxn) STOP 'HFB_nat_densit: wrong # of states'
        isonpl=nmaxp
      ELSE
        STOP 'HFB: ISO must be 1 or 2'
      END IF
c
      acctot=0D0
      DO nq=1,nmaxa
        aux1=0D0
        aux2=0D0
        aux3=0D0
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
            aux1(1:ngrid)=cm(nplshl,nq)*wf(npl:nplm+ngrid)
     &                    +aux1(1:ngrid)
            aux3(1:ngrid)=cm(nplshl,nq)*wfs(npl:nplm+ngrid)
     &                    +aux3(1:ngrid)
!            aux2(1:ngrid)=SQRT(MAX((1D0-v2a(nplshl))*v2a(nplshl),1D-40))
!     &                     *cm(nplshl,nq)*wf(npl:nplm+ngrid)
!     &                    +aux2(1:ngrid)
            v2add=v2a(nplshl)*cm(nplshl,nq)**2
!            WRITE(6,'(a,3i3,1pg13.5)') 'nq,js,n,v2add=',nq,js,n,v2add
            accno = v2add + accno
            acctot = v2add*phasfac(nplshl) + acctot
            IF(v2add>v2max) THEN
              jsnq=js
              v2max=v2add
            END IF
          END DO
        END DO
!        WRITE(6,'(a,(/5(1pg13.5)))') 'ratio U/V:',
!     &        aux2(1:ngrid)/aux1(1:ngrid)
        IF(tprint0)
     &   WRITE(6,'(a,2i4,3(1pg13.5))') 'nq,js,norms=',nq,jsnq,
     &        rstep*2D0*zwj(js)/0.15915494D0*SUM(aux1*aux1),
     &        rstep*2D0*zwj(js)/0.15915494D0*SUM(aux1*aux2),
     &        accno
c
        js=jssav(nq) ! jsnq
        ll     = js/2
        xlp    = ll*ll+ll+1D0
        wl     = wls(js)
        sig    = 2D0*MOD(ll,2)-1D0
!        CALL derv15(aux1,aux3,sig,ngrid,rstep)
        om     = zwj(js)*v2a(nq)
        ompair = zwj(js)*uva(nq)/2D0
        IF(tprint) WRITE(6,'(a,2i5,5(1pg13.5))') 'js etc:',
     &         js,ll,wl,sig,om,aux1(1)
!              nadd   = nadd+ngrid
c                                    the case r=0.0 needs special treatment
        if(ll.lt.1) then
          x      = aux3(1)
          rho(im) = (x*x)*om + rho(im)
          paircr(im) = (x*aux1(1))*ompair + paircr(im)
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
          wfa    = aux1(i-isoadd)
          wfa2   = wfa*wfa
c         if(wfa2.gt.xsmall) then
          rho(i) = (phas*wfa2) + rho(i)
          taul(i) = (phas*wfa2)*xlp + taul(i)
          sj(i)  = (phas*wfa2)*wl*ri + sj(i)
c  
          wfas   = aux3(i-isoadd)
          tau(i)=(wfas*wfas-((wfa*wfas)+(wfa*wfas)-xlp*wfa2*ri)*ri)
     &             *phas  + tau(i)
          paircr(i) = ompair*(ri*ri)*wfa2 + paircr(i)
        END DO
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
        WRITE(6,'(a,4i3)') 'HFB-natural densit: iso=',iso,im,ia,ie
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
c-----hfb_nat_densit---------------------------------------------------------
c
      subroutine hfb_nat_mapwf(iso,nmaxa,cm,v2a,paircut,jssav,phasfac)
c
      include 'params.inc'
c
c     Map new wavefunctions from configuration space to coordinate space.
c
      include 'wavfun.inc'
      include 'potinv.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      INTEGER, INTENT(IN) :: nmaxa,iso,jssav(nmaxa)
      REAL(8), INTENT(IN) :: cm(nmaxa,nmaxa),v2a(nmaxa),paircut(nmaxa)
      REAL(8), INTENT(IN) :: phasfac(nmaxa)

      LOGICAL,PARAMETER :: tprint0=.FALSE.,tprint=.FALSE.
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
      ngw    = (iso-1)*nmaxp
      nadd   = ngw*ngrid-ie
      IF(iso==1) THEN
        IF(nmaxa.NE.nmaxp) STOP 'HFB_nat_densit: wrong # of states'
        isonpl=0
      ELSE IF(iso==2) THEN
        IF(nmaxa.NE.nmaxn) STOP 'HFB_nat_densit: wrong # of states'
        isonpl=nmaxp
      ELSE
        STOP 'HFB: ISO must be 1 or 2'
      END IF
c
c     new wave functions from HFB coefficients
c
      acctot=0D0
      DO js=1,jsmax(iso)
        jsnpl = js+isnpad
        nmax  = npla(1,jsnpl)
        DO n2=1,nmax
          nplq    = npla(n2+1,jsnpl)
          nplmq = nplq-1
          nq = nplq/ngrid+1-isonpl
          aux1=0D0
          DO n=1,nmax
            npl    = npla(n+1,jsnpl)
            nplm = npl-1
            nplshl = npl/ngrid+1-isonpl
            IF(jssav(nplshl).NE.js) STOP 'HFB_deriv: js not matching'
            aux1(1:ngrid)=cm(nplshl,nq)*wf(npl:nplm+ngrid)
     &                    +aux1(1:ngrid)
          END DO
          savlr(nplq:nplmq+ngrid,1)=aux1(1:ngrid)
          nqi=nq+isonpl
!          WRITE(6,'(3i5,4(1pg13.5))') 
!     &     js,n2,nq,weight(nqi),v2a(nq),phpair(nqi),paircut(nq)
          weight(nqi)=v2a(nq)
          phpair(nqi)=paircut(nq)
        END DO
      END DO
!      IF(ISO==2) STOP
c
c     move to correct place
c
      DO js=1,jsmax(iso)
        jsnpl = js+isnpad
        nmax  = npla(1,jsnpl)
        ll     = js/2
        sig    = 2D0*MOD(ll,2)-1D0
        DO n=1,nmax
          npl    = npla(n+1,jsnpl)
          nplm = npl-1
          nplshl = npl/ngrid+1-isonpl
          wf(npl:nplm+ngrid) = savlr(npl:nplm+ngrid,1)
          call derv17(wf(npl),wfs(npl),sig,ngrid,rstep)
          IF(thfbfold) CALL Gauss_folding_m(wf(npl),wffold(npl),rstep,
     &                         ngrid,js,hfbwidth,MAXVAL(jsmax(1:2)))
        END DO
      END DO
c

      RETURN
      END
