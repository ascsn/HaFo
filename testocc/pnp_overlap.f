c-----pnp_overlap------------ part of the hartree-fock package ----------------
c
      subroutine pnp_overlap(v2occ,ph,pairgw,spe,nmax,ndim,iso,
     &                       nprot,nneut,rho,chi)
c
      include 'params.inc'
c
c     the include file 'params.inc' allows
c     to make the whole subroutine double precision by
c     transfering a line 'implicit double precision (a-h,o-z)'.
C
c     the routine evaluates and prints a couple of overlaps (norm,
c      one-body H, two-body H, ...) for particle-number projection.
c      v2occ   = array(1..nmax) of the occupation probabilities
c      nmax    = maximum number of single particle states
c                (needs to be smaller than the parameter 'ndim')
c      ndim    = dimensioning of the arrays, particularly of 'vpairf',
c                in the calling routine
c      iso     = isospin
c
      include 'grid.inc'
      REAL(8) :: v2occ(ndim),ph(ndim),pairgw(ndim),spe(ndim)
      REAL(8) :: rho(ngrid),chi(ngrid)

      INTEGER,PARAMETER :: nphi=100
      INTEGER,PARAMETER :: nphi2=2*nphi
      REAL(8) :: accphase(-2:nphi2),delphase(-2:nphi2)
      COMPLEX(8) :: overln(-2:nphi2),overlh(-2:nphi2),overlh2(-2:nphi2),
     &              overlp(-2:nphi2)
      COMPLEX(8) :: accn,acch,acch2,accp,apph,apph2,appp
      COMPLEX(8) :: accnln,acchln,acch2ln,accpln
      COMPLEX(8) :: accng,acchg,acch2g,accpg
      COMPLEX(8) :: curvn,curvh,curvh2,curvp, slopn,sloph,sloph2,slopp
      COMPLEX(8) :: skewn,curv4n
      COMPLEX(8) :: overold,overnew,epair2
      COMPLEX(8) :: crho(KGR2),cchi(KGR2),cchib(KGR2),cchi2(KGR2)
      COMPLEX(8) :: cv2wgt(KSL2),cuvwgt(KSL2),cuvwgtb(KSL2)
      COMPLEX(8) :: cv2acc(KSL2),cuvacc(KSL2),cuvaccb(KSL2)  ! projected v**2 , u*v
      COMPLEX(8) :: ci=CMPLX(0D0,1D0)
      REAL(8),PARAMETER :: pi=3.141592654D0
c
c-----------------------------------------------------------------------
c
      IF(ISO==1) OPEN(877,file='normoverlap.res')
      WRITE(877,'(a,i4/a)') '# overlaps for iso=',iso,
     &     '#  phi  norm_abs  norm_phase h_abs h_phase h2_abs h2_phase'
       cchi2=CMPLX(0D0,0D0)
      cv2acc=CMPLX(0D0,0D0)
      cuvacc=CMPLX(0D0,0D0)
      cuvaccb=CMPLX(0D0,0D0)
      accn=CMPLX(0D0,0D0)
      acch=CMPLX(0D0,0D0)
      acch2=CMPLX(0D0,0D0)
      accp=CMPLX(0D0,0D0)
      accphase=0D0
      oldphase=0D0
      delphi=1D0/nphi*pi
      DO n=-2,nphi*2
        phi=n*delphi
        overln(n)=CMPLX(1D0,0D0)
        overlh(n)=CMPLX(0D0,0D0)
        partnm=0D0
        DO i=1,nmax
          overnew = v2occ(i)*(EXP(CMPLX(0D0,phi))-1D0)
          cv2wgt(i)=(overnew+v2occ(i))/(1D0+overnew)
          cuvwgt(i)=SQRT(v2occ(i)-v2occ(i)**2)*EXP(CMPLX(0D0,phi))/
     &              (1D0+overnew)
          cuvwgtb(i)=SQRT(v2occ(i)-v2occ(i)**2)/(1D0+overnew)
          overlh(n)=overlh(n) + spe(i)*ph(i)*cv2wgt(i)
          overnew=(1+overnew)**ph(i)
!          addphase=ATAN2(AIMAG(overnew),REAL(overnew))
          overln(n)=overln(n)*overnew
!     &           (1+v2occ(i)*(EXP(CMPLX(0D0,phi))-1D0))**ph(i)
!     &                   ph(i)*pairgw(i))
          partnm=partnm+v2occ(i)*ph(i)
        END DO
        CALL cdensit(cchib,cuvwgtb,iso)
        CALL cdensit(cchi,cuvwgt,iso)
        cchib(1:ngrid)=cchib(1:ngrid)/2D0
        cchi(1:ngrid)=cchi(1:ngrid)/2D0
        crho(1:ngrid)=rfact(1:ngrid)**2*cchi(1:ngrid)*cchib(1:ngrid)
        overlp(n)=SUM(crho(1:ngrid))*frpirs
        CALL cdensit(crho,cv2wgt,iso)
        crho(1:ngrid)=rfact(1:ngrid)**2*crho(1:ngrid)**2
        overlh2(n)=SUM(crho(1:ngrid))*frpirs
        IF(n==0) THEN
          accspe=REAL(overlh(n))
          acc2b =REAL(overlh2(n))
          accp0 =REAL(overlp(n))
        END IF
        overln(n)=overln(n)*EXP(CMPLX(0D0,-partnm*phi))
!        overlh(n)=overlh(n)*overln(n)
!        overlh2(n)=overlh2(n)*overln(n)
!        overlp(n)=overlp(n)*overln(n)
!        overln=overlh*overln
        IF(n>0) THEN
          cchi2(1:ngrid) =cchi2(1:ngrid)+
     &                    overln(n)*cchi(1:ngrid)*cchib(1:ngrid)
          cv2acc(1:nmax) = cv2acc(1:nmax) +  overln(n)*cv2wgt(1:nmax)
          cuvacc(1:nmax) = cuvacc(1:nmax) +  overln(n)*cuvwgt(1:nmax)
          cuvaccb(1:nmax) = cuvaccb(1:nmax) +  overln(n)*cuvwgtb(1:nmax)
          accn = accn + overln(n)
          acch = acch + overlh(n)*overln(n)
          acch2 = acch2 + overlh2(n)*overln(n)
          accp = accp + overlp(n)*overln(n)
          overold=overln(n)/overold
          addphase=ATAN2(AIMAG(overold),REAL(overold))
          accphase(n)=0.5D0*(addphase+oldphase)+accphase(n-1)
          delphase(n)=(addphase-oldphase)/delphi
          oldphase=addphase
        END IF
        overold=overln(n)
!        overlh(n)=-overlh(n)

        WRITE(877,'(10(1pg13.5))') phi,ABS(overln(n)),accphase(n)/pi,
     &     ABS(overlh(n)),ATAN2(AIMAG(-overlh(n)),REAL(-overlh(n)))/pi,
     &     ABS(overlh2(n)),ATAN2(AIMAG(overlh2(n)),REAL(overlh2(n)))/pi,
     &     ABS(overlp(n)),ATAN2(AIMAG(overlp(n)),REAL(overlp(n)))/pi,
     &     delphase(n)

      END DO
      exacth=REAL(acch/accn)
      exacth2=REAL(acch2/accn)
      exactp=REAL(accp/accn)
      cchi2(1:ngrid) =cchi2(1:ngrid)/accn

! print projected occupation numbers
      cv2acc(1:nmax)=cv2acc(1:nmax)/accn
      cuvacc(1:nmax)=cuvacc(1:nmax)/accn
      cuvaccb(1:nmax)=cuvaccb(1:nmax)/accn
      WRITE(73,'(a)') '# iso istate  v2(mf)   uv(mf)  v2(pnp)   uv(pnp)'
      DO i=1,nmax
        WRITE(73,'(2i4,8(1pg13.5))') iso,i,v2occ(i),
     &         SQRT(v2occ(i)-v2occ(i)**2),REAL(cv2acc(i)),
     &         REAL(cuvacc(i)),REAL(cuvaccb(i)),
     &         REAL(cv2acc(i))+REAL(cuvacc(i)**2)/REAL(cv2acc(i)),
     &         REAL(cv2acc(i))+REAL(cuvaccb(i)**2)/REAL(cv2acc(i))
      END DO
! compute pairing energy from projected chi
      CALL cdensit(cchi,cuvacc,iso)
      CALL cdensit(cchib,cuvaccb,iso)
      cchib(1:ngrid)=cchib(1:ngrid)/2D0
      cchi(1:ngrid)=cchi(1:ngrid)/2D0
      crho(1:ngrid)=rfact(1:ngrid)**2*cchi(1:ngrid)*cchib(1:ngrid)
      epair2=SUM(crho(1:ngrid))*frpirs
      CALL cdensit(crho,cv2acc,iso)
      WRITE(73,'(a,3i4,6(1pg13.5))') '#epair:',iso,
     &    nprot,nneut,exactp,REAL(epair2),
     &    REAL(epair2)/(exactp),
     &     SUM(chi(1:ngrid)**2*rfact(1:ngrid)**2)*frpirs
      WRITE(73,'(a,3i4,6(1pg13.5))') '#E1,E2:',iso,
     &     nprot,nneut,exacth,
     &     SUM(spe(1:nmax)*ph(1:nmax)*REAL(cv2acc(1:nmax))),
     &     SUM(spe(1:nmax)*ph(1:nmax)*v2occ(1:nmax)),
     &     exacth2,
     &     SUM(REAL(crho(1:ngrid))**2*rfact(1:ngrid)**2)*frpirs,
     &     SUM(rho(1:ngrid)**2*rfact(1:ngrid)**2)*frpirs
      WRITE(73,'(a,3i4)') '# density for iso,nprot,nneut=',
     &  iso,nprot,nneut
      DO i=1,ngrid
         WRITE(73,'(i4,0pf8.2,6(1pg13.5))') 
     &        i,(i-1)*rstep,REAL(crho(i)),rho(i),
     &       REAL(cchi(i)),REAL(cchib(i)),chi(i),
     &        SQRT(REAL(cchi2(i)))
      END DO
!        WRITE(877,'(6(1pg13.5))') phi,ABS(overln),
!     &    ATAN2(AIMAG(overln),REAL(overln))/pi,
!     &    REAL(overln),AIMAG(overln),accphase/pi
!     &    accphase/pi,REAL(overln),AIMAG(overln)

!      WRITE(877,'(a,9(1pg13.5))') '# energies:',accspe,REAL(acch/accn),
!     &    acc2b,REAL(acch2/accn),accp0,REAL(accp/accn)


!
!  built topGOA, integrate to energy correction, and print
!

      WRITE(877,'(1x/1x)')

!      curvh=REAL(overlh(1)-overlh(0))/delphi**2
!      curvh2=REAL(overlh2(1)-overlh2(0))/delphi**2
!      curvp=REAL(overlp(1)-overlp(0))/delphi**2
!      sloph=AIMAg(overlh(1)-overlh(0))/delphi*ci
!      sloph2=AIMAG(overlh2(1)-overlh2(0))/delphi*ci
!      slopp=AIMAG(overlp(1)-overlp(0))/delphi*ci
      curvn=REAL(overln(1)-overln(0)*2D0+overln(-1))/delphi**2
      curvh=REAL(overlh(1)-overlh(0)*2D0+overlh(-1))/delphi**2/2D0
      curvh2=REAL(overlh2(1)-overlh2(0)*2D0+overlh2(-1))/delphi**2/2D0
      curvp=REAL(overlp(1)-overlp(0)*2D0+overlp(-1))/delphi**2/2D0
      slopn=AIMAg(overln(1)-overln(-1))/delphi*ci/2D0
      sloph=AIMAg(overlh(1)-overlh(-1))/delphi*ci/2D0
      sloph2=AIMAG(overlh2(1)-overlh2(-1))/delphi*ci/2D0
      slopp=AIMAG(overlp(1)-overlp(-1))/delphi*ci/2D0
      WRITE(877,'(a,8(1pg13.5))') '#slopes:',sloph,sloph2,slopp
      WRITE(877,'(a,8(1pg13.5))') '#curva.:',curvh,curvh2,curvp
      acch=CMPLX(0D0,0D0)
      acch2=CMPLX(0D0,0D0)
      accp=CMPLX(0D0,0D0)
      accng=CMPLX(0D0,0D0)
      acchg=CMPLX(0D0,0D0)
      acch2g=CMPLX(0D0,0D0)
      accpg=CMPLX(0D0,0D0)
      DO n=0,nphi*2
        phi=n*delphi
        phih=phi/2D0
        overlnGOA=EXP(2D0*curvn*SIN(phih)**2)
        apph=accspe+curvh*4D0*SIN(phih)**2+sloph*SIN(phi)
        apph2=acc2b+curvh2*4D0*SIN(phih)**2+sloph2*SIN(phi)
        appp=accp0+curvp*4D0*SIN(phih)**2+slopp*SIN(phi)
        WRITE(877,'(f8.5,16(1pg13.5))') phi,
     &     (overln(n)),overlnGOA,
     &     (overlh(n)),(apph),
     &     (overlh2(n)),(apph2),
     &     (overlp(n)),(appp)
        IF(n>0) THEN
          acch = acch + apph*overln(n)
          acch2 = acch2 + apph2*overln(n)
          accp = accp + appp*overln(n)
!        apph=accspe*EXP((curvh*4D0*SIN(phih)**2+sloph*SIN(phi))/accspe)
!        apph2=acc2b*EXP((curvh2*4D0*SIN(phih)**2
!     &                    +sloph2*SIN(phi))/acc2b)
!        appp=accp0*EXP((curvp*4D0*SIN(phih)**2+slopp*SIN(phi))/accp0)
          accng = accng + overlnGOA
          acchg = acchg + apph*overlnGOA
          acch2g = acch2g + apph2*overlnGOA
          accpg = accpg + appp*overlnGOA
        END IF
      END DO

!
! LN energy   (miraculous factor 4 yet to be found)
!
      WRITE(877,'(1x/1x)')

      overlh=overlh*overln
      overlh2=overlh2*overln
      overlp=overlp*overln
      curvn=REAL(overln(1)-overln(0)*2D0+overln(-1))/delphi**2
      curvh=REAL(overlh(1)-overlh(0)*2D0+overlh(-1))/delphi**2
      curvh2=REAL(overlh2(1)-overlh2(0)*2D0+overlh2(-1))/delphi**2
      curvp=REAL(overlp(1)-overlp(0)*2D0+overlp(-1))/delphi**2
      slopn=AIMAg(overln(1)-overln(-1))/delphi*ci/2D0
      sloph=AIMAg(overlh(1)-overlh(-1))/delphi*ci/2D0
      sloph2=AIMAG(overlh2(1)-overlh2(-1))/delphi*ci/2D0
      slopp=AIMAG(overlp(1)-overlp(-1))/delphi*ci/2D0
      skewn=(overln(2)-2D0*overln(1)+2D0*overln(-1)-overln(-2))/
     &      (2D0*delphi**3)
      curv4n=(overln(2)-4D0*overln(1)+6D0*overln(0)-4D0*overln(-1)
     &        +overln(-2))/delphi**4
      WRITE(877,'(a,8(1pg13.5))') '#slopes:',slopn,sloph,sloph2,slopp
      WRITE(877,'(a,8(1pg13.5))') '#curva.:',curvn,curvh,curvh2,curvp

      curvn=-curvn
      curvh=-curvh
      curvh2=-curvh2
      curvp=-curvp
      skewn=ci*skewn

      accnln = curv4n-curvn**2-skewn**2/curvn  ! denominator
      acchln = curvh-overlh(0)*curvn-AIMAG(sloph)*skewn/curvn
      acch2ln = curvh2-overlh2(0)*curvn-AIMAG(sloph2)*skewn/curvn
      accpln = curvp-overlp(0)*curvn-AIMAG(slopp)*skewn/curvn

! print all stuff
!      fac=SQRT(REAL(curv4n-curvn**2)/2D0)  !4D0
      fac= curvn*4D0*(curv4n-curvn**2)/(curv4n+5D0*curvn**2)
      WRITE(877,'(a,6(1pg13.5))') '# curv etc:',
     &    REAL(curvn),AIMAG(slopn),REAL(curv4n),fac
      WRITE(877,'(a,a)') '#type iso nprot nneut     mf     exact  ',
     &     '   topGOA-exno   topGOA-full     Lipkin-Nogami  GOA'
      WRITE(877,'(a,i2,2i4,7(1pg15.7))') 'energ1:',iso,nprot,nneut,
     &  accspe,exacth-accspe,REAL(acch/accn)-accspe,
     &  REAL(acchg/accng)-accspe,-REAL(curvn*acchln/accnln),
     &  -REAL(fac*acchln/accnln),-REAL(curvh/curvn-overlh(0))
      WRITE(877,'(a,i2,2i4,7(1pg15.7))') 'energ2:',iso,nprot,nneut,
     &  acc2b,exacth2-acc2b,REAL(acch2/accn)-acc2b,
     &  REAL(acch2g/accng)-acc2b,-REAL(curvn*acch2ln/accnln),
     &  -REAL(fac*acch2ln/accnln),-REAL(curvh2/curvn-overlh2(0))
      WRITE(877,'(a,i2,2i4,7(1pg15.7))') 'energp:',iso,nprot,nneut,
     &  accp0,exactp-accp0,REAL(accp/accn)-accp0,
     &  REAL(accpg/accng)-accp0,-REAL(curvn*accpln/accnln),
     &  -REAL(fac*accpln/accnln),-REAL(curvp/curvn-overlp(0))

      IF(iso==1) WRITE(877,'(1x/1x)')
      IF(ISO==2) CLOSE(877)
c
      RETURN
      END
c-----cdensit------------ part of the hartree-fock package ---------------------
c
      subroutine cdensit(cfield,caddwgt,iso)
c
      include 'params.inc'
c
c     accumulates the complex density field 'cfield' with
c     accumulation weight 'caddwgt' and wavefunctions from
c     COMMON WAVFUN.
c
      include 'grid.inc'
      include 'wavfun.inc'
      include 'constn.inc'

      COMPLEX(8) :: cfield(*),caddwgt(*)
      COMPLEX(8) :: om
c
c      dimension zwj(21)
c
c-------------------------------------------------------------------------------
c
      cfield(1:ngrid)=CMPLX(0D0,0D0)
c
      ie     = ngrid !*iso
      im     = ie-ngrid+1
      ia     = im+1
      if(iso.eq.1) then
       nmax   = nmaxp
      else
       nmax   = nmaxn
      endif
      ngw    = (iso-1)*nmaxp
      nadd   = ngw*ngrid-ie
c
c     loop over states for given isospin.
c      'nshl' counts the states in harmonic oscillator ordering.
c      'js' is the according composed ang.mom.,
c      'll' is the orbital ang.mom.,
c      'wl' is the spin-orbit splitting,
c      'om' is the actual phase-space factor
c
      DO nshl=1,nmax
        js     = jsp(nshl)
        ll     = js/2
        xlp    = float(ll*ll+ll+1)
        wl     = wls(js)
        IF(tphpair) THEN
!          om     = zwj(js)*caddwgt(nshl+ngw) 
          om     = zwj(js)*caddwgt(nshl) 
        ELSE
!          om     = zwj(js)*caddwgt(nshl+ngw)*phpair(nshl+ngw) 
          om     = zwj(js)*caddwgt(nshl)*phpair(nshl) 
        END IF
        nadd   = nadd+ngrid
        if(ABS(om)<xsmall) CYCLE
c                                  the case r=0.0 needs special treatment
        if(ll.lt.1) then
         x      = wfs(im+nadd)
         cfield(im) = (x*x)*om + cfield(im)
        endif
c                                  now all r>0.0
        r      = zero
        do i=ia,ie
          r      = r+rstep
          ri     = one/r
          wfa    = wf(nadd+i)
          cfield(i) = (ri*ri)*om*(wfa*wfa) + cfield(i)
        enddo
c
      END DO
c
      return
      end


c-----pnshift_overlap------------ part of the hartree-fock package ----------------
c
      subroutine pnshift_overlap(v2occ,ph,pairgw,spe,dgap,nmax,ndim,iso,
     &                       nprot,nneut)
c
c     scan continuous range of neutron numbers around 'nneut' and
c     compute norm overlap
c
      include 'params.inc'
      REAL(8) :: v2occ(ndim),ph(ndim),pairgw(ndim),spe(ndim),dgap(ndim)

      INTEGER,PARAMETER :: nshift=50
      REAL(8),PARAMETER :: shift=0.2D0
      REAL(8) :: v2s(ndim,-nshift:nshift),ttthrm(ndim)

      WRITE(611,'(a,5i4)') '#PNSHIFT with:',nmax,ndim,iso,nprot,nneut
      DO i=-nshift,nshift
        partnm=nneut+i*shift
        WRITE(6,'(a,i5,f8.3)') 'CALL PAIRDN with:',i,partnm
        call pairdn(spe,dgap,v2s(1,i),ph,pairgw,ttthrm,v2occ,
     &               1.65D0,partnm,efen,epn,deltn,deltn2,eqminn,gpn,
     &               efaddn,0D0,nmax,ndim,1,0,
     &               50,1D-6,-1,.FALSE.)
!        WRITE(6,'(a)') 'from PAIRDN:'
!        WRITE(611,'(a,i4,f8.3)') 'CALL PAIRDN with:',i,partnm 
!        WRITE(611,'(100(1pg13.5))') v2s(:,i)
      END DO
      DO i=-nshift,nshift
      DO j=-nshift,nshift
        prod=1D0
        DO n=1,nmax
          prod=prod*(SQRT(1D0-v2s(n,i)**2)*SQRT(1D0-v2s(n,j)**2)
     &               +v2s(n,j)*v2s(n,i))**ph(n)
        END DO
        WRITE(611,'(2f8.3,1pg13.5)') nneut+i*shift,nneut+j*shift,prod
      END DO
      WRITE(611,'(1x)')
      END DO

      return
      end
