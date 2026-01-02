c-----extobs_all--------------------------------------------------------
c
      subroutine extobs_all(par,workdv,nparmx,iextobs)
c
c     This routine is the general call to external observables
c     for error analysis. It occupies the switching field 'nextobs'
c     according to the given observable 'iextobs' and calls the
c     appropriate subroutiune.
c     Information on the wanted observable is communicated via
c     common block 'extobsc'.
c
      include 'params.inc'
c
      include 'extobs.inc'
      include 'numpar.inc'
      include 'comfit.inc'
      include 'skypar.inc'
      include 'constn.inc'
      dimension par(kparmx),workdv(kparmx,1)

c
c-------------------------------------------------------------------
c
      nextobs = 0            
      nextobs(iextobs) = 1
      ifplot = 0                      ! suppress extensive prints
!      WRITE(2,'(a,i3)') ' extobs: iextobs=',iextobs

      IF(iextobs <= 20) THEN
        call extobs(par,workdv,nparmx)
      ELSE
        call extobsmf(par,workdv,nparmx)
      END IF
      CALL flush(2)

      RETURN
      END

c-----extobs--------------------------------------------------------
c
      subroutine extobs(par,workdv,nparmx)
c
c     Extrapolate to selected observables with extrapolation errors.
c     The error matrix 'zarray' is transferred via common /fitpar/ .
c     Observables are selected with setting options to 1:
c      igmrpb     = GMR in 208Pb
c      igdrpb     = GDR in 208Pb
c      iqgrpb     = GQR in 208Pb
c      igdro      = GDR in 16O
c      idipol     = dipole polarizability in 208Pb
c      idipsn     = dipole polarizability in 132Sn
c      idipca     = dipole polarizability in 48Ca
c      idipni     = dipole polarizability in 68Ni
c      ipygmy     = low energy sum rule in 208Pb (pygmy strength)
c      ipygsn     = low energy sum rule in 132Sn (pygmy strength)
c      imonop     = monopole polarizability in 208Pb
c      imonsn     = monopole polarizability in 132Sn
c      imonni     = monopole polarizability in 68Ni
c      imonca     = monopole polarizability in 48Ca
c      iquadp     = quadrupole polarizability in 208Pb
c     All resonances are estimated in local RPA.
c
      include 'params.inc'
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'comfit.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'extobs.inc'
      include 'constn.inc'
c
c     'workdv' is field for derivatives with respect to parameter-variation.
c
      dimension   par(kparmx),workdv(kparmx,1)
c
      logical ttest
      parameter (ttest=.true.)
c
c-----------------------------------------------------------------------
c
c     RPA code not yet suited for Fayans functional --> return
c
!      IF(ioptsk==13) RETURN
c
      itest =  igmrpb-igmrpb*igmrpb
     &       + igdrpb-igdrpb*igdrpb
     &       + igqrpb-igqrpb*igqrpb
     &       + igdro-igdro*igdro
     &       + idipol-idipol*idipol
     &       + idipsn-idipsn*idipsn
     &       + idipca-idipca*idipca
     &       + idipni-idipni*idipni
     &       + idipzr-idipzr*idipzr
     &       + ipygmy-ipygmy*ipygmy
     &       + ipygsn-ipygsn*ipygsn
     &       + imonop-imonop*imonop
      if(itest.ne.0) stop 'EXTOBS: options not 0 or 1'
      if(igmrpb+igdrpb+igqrpb+igdro+idipol+idipsn+idipca+idipzr
     &   +idipni+ipygmy+ipygsn+imonop+imonsn+imonni+imonca+iquadp
     &   .ne.1) then
        WRITE(2,'(43i2)') nextobs
        CALL flush(2)
        stop "EXTOBS: only one option can be active"
      endif
ctest      ifplot=3
ctest      write(2,'(a)') 'EXTOBS entered'
ctest      write(6,'(a)') 'EXTOBS entered'
c
c     loop over variation of parameters
c
      do i=0,nparmx
        if(i.gt.0) par(i) = par(i) + dparmo(i)
        call attach(par,jp)
        call writ05()
        call system('./haforpa.exe < for005 > for006')
        if(i.eq.0) then
         resensav = resenergy()
        else
         resenact = resenergy()
         workdv(i,1) = (resenact-resensav)/dparmo(i)
         endif
c
         if(i.gt.0) par(i) = par(i) - dparmo(i)
      enddo
c
c     accumulate extrapolation errors by folding with 'zarray'.
c
      erresen   = 0.0
c  
      if(nparmx.gt.0) then
       do i=1,nparmx
       do j=1,nparmx
         erresen   = workdv(i,1)*workdv(j,1)*zarray(i,j) 
     &              + erresen
       enddo
       enddo
       erresen   = sqrt(erresen)
      endif

      if(ttest) then
        write(6,'(a,17g12.4)') 
     &   ' workdv1:',(workdv(i,1),i=1,nparmx)
        write(6,'(a)') ' zarray:'
        do i=1,nparmx
          write(6,'(17g12.4)') (zarray(j,i),j=1,nparmx)
        enddo
      endif
c
      if(igmrpb.eq.1) then
        write(02,'(a,f8.3,1pg16.8,1pg12.4)')
     &    ' GMR in 208PB:  goal,energy,uncertainty=',
     &    13.7,resensav,erresen
        if(ttest) write(06,'(a,f8.3,1pg16.8,1pg12.4)')
     &    ' GMR in 208PB:  goal,energy,uncertainty=',
     &    13.7,resensav,erresen
        write(81,'(a)') ' GMR 208Pb'
      endif
      if(igdrpb.eq.1) then
        write(02,'(a,f8.3,1pg16.8,1pg12.4)')
     &    ' GDR in 208PB:  goal,energy,uncertainty=',
     &    13.6,resensav,erresen
        if(ttest) write(06,'(a,f8.3,1pg16.8,1pg12.4)')
     &    ' GDR in 208PB:  goal,energy,uncertainty=',
     &    13.6,resensav,erresen
        write(81,'(a)') ' GDR 208Pb'
      endif
      if(igqrpb.eq.1) then
        write(02,'(a,f8.3,1pg16.8,1pg12.4)')
     &    ' GQR in 208PB:  goal,energy,uncertainty=',
     &    10.9,resensav,erresen
        if(ttest) write(06,'(a,f8.3,1pg16.8,1pg12.4)')
     &    ' GQR in 208PB:  goal,energy,uncertainty=',
     &    10.9,resensav,erresen
        write(81,'(a)') ' GQR 208Pb'
      endif
      if(igdro.eq.1) then
        write(02,'(a,f8.3,1pg16.8,1pg12.4)')
     &    ' GDR in 16O:  goal,energy,uncertainty=',
     &    23.0,resensav,erresen
        if(ttest) write(06,'(a,f8.3,1pg16.8,1pg12.4)')
     &    ' GDR in 16O:  goal,energy,uncertainty=',
     &    23.0,resensav,erresen
        write(81,'(a)') ' GDR 16O'
      endif
      if(idipol.eq.1) then
        write(02,'(a,1pg16.8,1pg12.4)')
     &    ' polariz. in 208Pb: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,1pg16.8,1pg12.4)')
     &    ' polariz. in 208Pb: value,uncertainty=',
     &       resensav,erresen
        write(81,'(a)') ' polariz. 208Pb'
      elseif(idipca.eq.1) then
        write(02,'(a,1pg16.8,1pg12.4)')
     &    ' polariz. in 48Ca: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,1pg16.8,1pg12.4)')
     &    ' polariz. in 38Ca: value,uncertainty=',
     &       resensav,erresen
        write(81,'(a)') ' polariz. 48Ca'
      elseif(idipsn.eq.1) then
        write(02,'(a,1pg16.8,1pg12.4)')
     &    ' polariz. in 132Sn: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,1pg16.8,1pg12.4)')
     &    ' polariz. in 132Sn: value,uncertainty=',
     &       resensav,erresen
        write(81,'(a)') ' polariz. 132Sn'
      elseif(idipni.eq.1) then
        write(02,'(a,1pg16.8,1pg12.4)')
     &    ' polariz. in 68Ni: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,1pg16.8,1pg12.4)')
     &    ' polariz. in 68Ni: value,uncertainty=',
     &       resensav,erresen
      elseif(idipzr.eq.1) then
        write(02,'(a,1pg16.8,1pg12.4)')
     &    ' polariz. in 90Zr: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,1pg16.8,1pg12.4)')
     &    ' polariz. in 90Zr: value,uncertainty=',
     &       resensav,erresen
        write(81,'(a)') ' polariz. 90Zr'
      elseif(ipygmy.eq.1) then
        write(02,'(a,1pg16.8,1pg12.4)')
     &    ' pygmy res. 208Pb: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,1pg16.8,1pg12.4)')
     &    ' pygmy res. 208Pb: value,uncertainty=',
     &       resensav,erresen
        write(81,'(a)') ' pygmy res. 208Pb'
      elseif(ipygsn.eq.1) then
        write(02,'(a,1pg16.8,1pg12.4)')
     &    ' pygmy res. 132Sn: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,1pg16.8,1pg12.4)')
     &    ' pygmy res. 132Sn: value,uncertainty=',
     &       resensav,erresen
        write(81,'(a)') ' pygmy res. 132Sn'
      elseif(imonop.eq.1) then
        write(02,'(a,1pg16.8,1pg12.4)')
     &    ' monop. polariz. in 208Pb: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,1pg16.8,1pg12.4)')
     &    ' monop. polariz. in 208Pb: value,uncertainty=',
     &       resensav,erresen
        write(81,'(a)') ' monop.pol. 208Pb'
      elseif(imonsn.eq.1) then
        write(02,'(a,1pg16.8,1pg12.4)')
     &    ' monop. polariz. in 132Sn: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,1pg16.8,1pg12.4)')
     &    ' monop. polariz. in 132Sn: value,uncertainty=',
     &       resensav,erresen
        write(81,'(a)') ' monop.pol. 132Sn'
      elseif(imonni.eq.1) then
        write(02,'(a,1pg16.8,1pg12.4)')
     &    ' monop. polariz. in 68Ni: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,1pg16.8,1pg12.4)')
     &    ' monop. polariz. in 68Ni: value,uncertainty=',
     &       resensav,erresen
        write(81,'(a)') ' monop.pol. 68Ni'
      elseif(imonca.eq.1) then
        write(02,'(a,1pg16.8,1pg12.4)')
     &    ' monop. polariz. in 48Ca: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,1pg16.8,1pg12.4)')
     &    ' monop. polariz. in 48Ca: value,uncertainty=',
     &       resensav,erresen
        write(81,'(a)') ' monop.pol. 48Ca'
      elseif(iquadp.eq.1) then
        write(02,'(a,1pg16.8,1pg12.4)')
     &    ' quadr. polariz. in 208Pb: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,1pg16.8,1pg12.4)')
     &    ' quadr. polariz. in 208Pb: value,uncertainty=',
     &       resensav,erresen
        write(81,'(a)') ' quadr.pol. 208Pb'
      endif
c
c     collect
c
      write(81,'(1pg24.16)') 
     &    resensav
      write(81,'(100(1pg24.16))') 
     &    (workdv(i,1),i=1,nparmx)

c
      return
      end
c-----writ05--------------------------------------------------------
c
      subroutine writ05()
c
c     Write input file 'for005' for calculation of giant resonances.
c      igdrpb     = GDR in 208Pb
c      igqrpb     = GQR in 208Pb
c      igdro      = GDR in 16O
c      idipol     = dipole polarizability in 208Pb
c      idipsn     = dipole polarizability in 132Sn
c      idipca     = dipole polarizability in 48Ca
c
      include 'params.inc'
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'comfit.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'extobs.inc'
      include 'constn.inc'
c
      data iu/31/               ! unit for output
c
c-----------------------------------------------------------------------
c
      if(igmrpb+igdrpb+igqrpb+igdro+idipol+idipsn+idipca+idipzr
     &   +idipni+ipygmy+ipygsn+imonop+imonsn+imonni+imonca+iquadp
     &   .ne.1)  stop "WRIT05: no valid option "
c
!      WRITE(2,*) ' WRIT05:',b0,b1,b2,b3,b0p,b1p,b2p,b3p,power
      open(iu,file='for005')
c
c     the mean field part
c
      write(iu,'(a)')
     &   '1,,    IWUNIT',
     &   "' resonances ion 208 Pb'",
     &   '0,0.3,              NGRID,RSTEP',
     &   '-1,1E-10,1E-11,278,5780,3,-1, IWR,DABR/E,ITPR/MAX,IFPL',
     &   '0.2,90.,0.2,-1,,               X0DAMP,E0DAMP,ADDNEW,WALL'
      write(iu,'(9(i2,a),i2,a,5x,a)') 
     &    ifex,',',iftm,',',ifenrg,',',ifrso,',',ifzpe,',',
     &    ifdfms,',',ifrhoc,',',ifenav,',',imfold,',',ioptsk,
     &    ', if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk'
      write(iu,'(a)')
     &   '1,,,                       ifform',
     &   '0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12',
     &   '0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12'
      IF(ioptsk==13) THEN            ! Fayans functional
          WRITE(iu,'(3(1pg18.10),a)') avp,h1vp,h2vp,' avp,h1vp,h2vp'
          WRITE(iu,'(5(1pg18.10))') avm,h1vm,h2vm
          WRITE(iu,'(5(1pg18.10))') sigmavol,rho0fy
          WRITE(iu,'(5(1pg18.10))') asp,hsp,hgradp 
          WRITE(iu,'(5(1pg18.10))') asm,hsm,hgradm 
          WRITE(iu,'(5(1pg18.10))') hcoul,sigmacoul
          WRITE(iu,'(5(1pg18.10))') fxi,h1xi,hgradxi,gammaxi,rho0xi
          WRITE(iu,'(4(1pg18.10),i2)') delfxi,ecstab,delh1xi,delhgradxi,
     &            isovectormodel
          WRITE(iu,'(5(1pg18.10))') c0tau,c1tau 
          WRITE(iu,'(5(1pg18.10))') c0nabj,c1nabj,c0tens,c1tens
      ELSE IF(ioptsk==12) THEN
        write(iu,'(5(1pg20.12,a)/5(1pg20.12,a),a)')
     &  b0,',',b1,',',b2,',',b3,',',t4,',      t0,t1,t2,t3,t4',
     &  b0p,',',b1p,',',b2p,',',b3p,',',b4p,',  x0,x1,x2,x3,b4p'
      ELSE
        write(iu,'(5(1pg20.12,a)/5(1pg20.12,a),a)')
     &  t0,',',t1,',',t2,',',t3,',',t4,',      t0,t1,t2,t3,t4',
     &  x0,',',x1,',',x2,',',x3,',',b4p,',       x0,x1,x2,x3,b4p'
      END IF
      if(ioptsk.eq.5 .or. ioptsk.eq.6 .or. ioptsk.eq.11) then
         write(iu,'(2(1pg20.12,a)/2(1pg20.12,a))')
     &      c1,',',c4,',                      c1 c4',
     &      c1p,',',c4p,',                    c1p,c4p'
      endif
      if(ioptsk.eq.8 .or. ioptsk.eq.11) then
         write(iu,'(3(1pg20.12),a/2(1pg20.12,a))')
     &      b3tilde,tte,tto,',       b3tilde,tte,tto',
     &      b12,',',b12p,',                     b12,b12p'
      endif
      if(ioptsk.eq.9 .or. ioptsk.eq.11) then
         write(iu,'(2(1pg20.12,a)/1pg20.12,a)')
     &      b7,',',b7p,',               b7,b7p',
     &      pow7,',               pow7'
      endif
      if(ioptsk.eq.10 .or. ioptsk.eq.11) then
         write(iu,'(8(1pg20.12),a)')
     &      riftm,rho0pr,powerp,ecstab,alpinv,betinv,
     &      alpin2,betin2,
     &    ',      riftm,rho0pr,powerp,ecstab,alpinv'
      endif
      if(ioptsk.le.4 .or. ioptsk.eq.7) then
         write(iu,'(5(f10.3,a)/4(f10.6,a))')
     &    t5,',',t6,',',t7,',',t8,',',t9,',    t5,t6,t7,t8,t9',
     &    x5,',',x6,',',x7,',',x8,',           x5,x6,x7,x8'
      endif
       if(ioptsk.NE.13) write(iu,'((f10.3,a),2(f10.6,a))')
     &    tyuk,',',xyuk,',',ayuk,',           tyuk,xyuk,ayuk'
       write(iu,'(f10.6,a,a)')
     &     power,',1.0,-1,10.0,10.0',
     &  ',0,,,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,...'
        if(igmrpb.eq.1 .or. igdrpb.eq.1 .or. igqrpb.eq.1 
     &     .or. idipol.eq.1 .or. imonop.eq.1 .or. iquadp.eq.1) 
     &   write(iu,'(a)') 
     &    "82,126,16,22,,,  nprot,nneut,nmaxp,nmaxn,iparsh",
     &    "1         IFGSC",
     &    "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,,,,,,,,,,,,",
     &    "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,,,,,,"
        if(ipygmy.eq.1) 
     &   write(iu,'(a/a/2a/2a)') 
     &    "82,126,37,46,,,  nprot,nneut,nmaxp,nmaxn,iparsh",
     &    "1         IFGSC",
     &    "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,,,,,,,,,,,,,,,,,",
     &    ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,",
     &    "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,,,,,,",
     &    ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"
        if(igdro.eq.1) 
     &   write(iu,'(a)') 
     &    "8,8,3,3,,,  nprot,nneut,nmaxp,nmaxn,iparsh",
     &    "1         IFGSC",
     &    "1,1,1,,,,,,,,,,,,,,",
     &    "1,1,1,,,,,,,,,,,,,,"
        if(idipsn.eq.1 .or. imonsn.eq.1) 
!     &   write(iu,'(a)') 
!     &    "20,20,6,6,,,,  nprot,nneut,nmaxp,nmaxn,iparsh",
     &   write(iu,'(a)') 
     &    "50,82,11,16,,,  nprot,nneut,nmaxp,nmaxn,iparsh",
     &    "1         IFGSC",
     &    "1,1,1,1,1,1,1,1,1,1,1,,,,,,,,,,,,,",
     &    "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,,,,,,,,,,,,,"
        if(ipygsn.eq.1) 
     &   write(iu,'(a/a/a/2a)') 
     &    "50,82,29,37,,,  nprot,nneut,nmaxp,nmaxn,iparsh",
     &    "1         IFGSC",
     &    "1,1,1,1,1,1,1,1,1,1,1,,,,,,,,,,,,,,,,,,,,,,,,,",
     &    "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,,,,,,,,,,,,,,,",
     &    ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,"
        if(idipca.eq.1 .or. imonca.eq.1) 
     &   write(iu,'(a)') 
     &    "20,28,6,7,,,  nprot,nneut,nmaxp,nmaxn,iparsh",
     &    "1         IFGSC",
     &    "1,1,1,1,1,1,,,,,,,,,,,,,,",
     &    "1,1,1,1,1,1,1,,,,,,,,,,,,,,"
        if(idipni.eq.1 .or. imonni.eq.1) 
     &   write(iu,'(a)') 
     &    "28,40,7,10,,,  nprot,nneut,nmaxp,nmaxn,iparsh",
     &    "1         IFGSC",
     &    "1,1,1,1,1,1,1,,,,,,,,,,,,,",
     &    "1,1,1,1,1,1,1,1,1,1,,,,,,,"
        if(idipzr.eq.1) 
     &   write(iu,'(a)') 
     &    "40,50,10,11,,,  nprot,nneut,nmaxp,nmaxn,iparsh",
     &    "1         IFGSC",
     &    "1,1,1,1,1,1,1,1,1,1,,,,,,,,,,,,",
     &    "1,1,1,1,1,1,1,1,1,1,1,,,,,,"
      if(igdrpb.eq.1) call system('cp for005 for005.GRhead')
c
c     GMR specific part
c
      if(igmrpb.eq.1) then
        write(iu,'(a)') 
     &    "0,1.0,0.0000001,1.0E-5,0,,0,200.0,,,",
     &    "100,0.1,0.0001,-1,0,0,0,,,,        ",
     &    " 'mom=0',1,'nojj',0,0,0,1,0,,,",
     &    " 'mom=0',2,'nojj',0,0,0,1,0,,",
     &    " 'mom=2',2,'nojj',0,0,0,1,0,,",
     &    " 'mom=2',1,'nojj',0,0,0,1,0,,",
     &    " 'mom=4',2,'nojj',0,0,0,1,0,,",
     &    " 'mom=4',1,'nojj',0,0,0,1,0,,",
     &    " 'mom=6',2,'nojj',0,0,0,1,0,,",
     &    " 'mom=6',1,'nojj',0,0,0,1,0,,",
     &    " 'bes=1.0',2,'nojj',0,0,0,1,0,,",
     &    " 'bes=1.0',1,'nojj',0,0,0,1,0,,",
     &    " 'end',,,,,,,,,,",
     &    "  7,20,0.1,0.5,2, ",
     &    " -1,                   IWUNIT"
      endif
c
c     GDR specific part
c
      if(igdrpb.eq.1 .or. igdro.eq.1) then
        write(iu,'(a)') 
     &    "1,1.0,0.0000001,1.0E-5,0,,0,200.0,,,",
     &    "100,0.1,0.0001,-1,0,0,0,,,,        ",
     &    " 'mom=0',1,'nojj',0,0,0,1,0,,,",
     &    " 'mom=0',2,'nojj',0,0,0,1,0,,",
     &    " 'mom=2',2,'nojj',0,0,0,1,0,,",
     &    " 'mom=2',1,'nojj',0,0,0,1,0,,",
     &    " 'mom=4',2,'nojj',0,0,0,1,0,,",
     &    " 'mom=4',1,'nojj',0,0,0,1,0,,",
     &    " 'mom=6',2,'nojj',0,0,0,1,0,,",
     &    " 'mom=6',1,'nojj',0,0,0,1,0,,",
     &    " 'bes=1.0',2,'nojj',0,0,0,1,0,,",
     &    " 'bes=1.0',1,'nojj',0,0,0,1,0,,",
     &    " 'end',,,,,,,,,,",
     &    "  7,20,0.1,0.5,2, ",
     &    " -1,                   IWUNIT"
        if(igdrpb.eq.1) call system('cp for005 for005.GQR')
      endif
c
c     pygmy resonance specific part
c
      if(ipygmy.eq.1 .or. ipygsn.eq.1) then
        write(iu,'(a)') 
     &    "1,1.0,0.0000001,1.0E-5,0,,0,200.0,,,",
     &    "100,0.1,0.0001,-1,0,0,0,,,,        ",
     &    " 'mom=0',1,'nojj',0,0,0,1,0,,,",
     &    " 'mom=0',2,'nojj',0,0,0,1,0,,",
     &    " 'pha=0,20',2,'nojj',0,0,0,0,0,,",
     &    " 'mom=4',2,'nojj',0,0,0,1,0,,",
     &    " 'mom=4',1,'nojj',0,0,0,1,0,,",
     &    " 'mom=6',2,'nojj',0,0,0,1,0,,",
     &    " 'mom=6',1,'nojj',0,0,0,1,0,,",
     &    " 'bes=1.0',2,'nojj',0,0,0,1,0,,",
     &    " 'bes=1.0',1,'nojj',0,0,0,1,0,,",
     &    " 'bes=1.4',2,'nojj',0,0,0,1,0,,",
     &    " 'bes=1.4',1,'nojj',0,0,0,1,0,,",
     &    " 'end',,,,,,,,,,",
     &    "  7,20,0.1,0.5,2, ",
     &    " -1,                   IWUNIT"
      endif
c
c     GQR specific part
c
      if(igqrpb.eq.1) then
        write(iu,'(a)') 
     &    "2,1.0,0.0000001,1.0E-5,0,,0,200.0,,,",
     &    "600,0.03,0.00001,-1,0,0,0,,,,        ",
     &    " 'mom=0',1,'nojj',0,-6,0,-1,0,,,",
     &    " 'mom=0',1,'nojj',0,0,0,1,0,,,",
     &    " 'mom=0',2,'nojj',0,0,0,1,0,,",
     &    " 'mom=2',2,'nojj',0,0,0,1,0,,",
     &    " 'mom=2',1,'nojj',0,0,0,1,0,,",
     &    " 'mom=4',2,'nojj',0,0,0,1,0,,",
     &    " 'mom=4',1,'nojj',0,0,0,1,0,,",
     &    " 'mom=6',2,'nojj',0,0,0,1,0,,",
     &    " 'mom=6',1,'nojj',0,0,0,1,0,,",
     &    " 'bes=1.0',2,'nojj',0,0,0,1,0,,",
     &    " 'bes=1.0',1,'nojj',0,0,0,1,0,,",
     &    " 'end',,,,,,,,,,",
     &    "  7,20,0.1,0.5,2, ",
     &    " -1,                   IWUNIT"
        call system('cp for005 for005.GQR')
      endif
c
c     part for dipole polarizability
c
!      if(idipol.eq.1 .or. idipsn.eq.1 .or. idipca.eq.1
!     &    .or. idipni.eq.1 .or. idipzr.eq.1) then
      if(idipol.eq.1 .or. idipsn.eq.1
     &    .or. idipni.eq.1 .or. idipzr.eq.1) then
        write(iu,'(a)') 
     &    "1,1.0,1D-12,1.0E-5,0,,0,200.0,,,",
     &    "5000,0.02,1D-10,-1,0,0,0,,,,        ",
     &    " 'mom=0',2,'nojj',0,0,0,-1,0,,",
     &    " 'end',,,,,,,,,,",
     &    "  7,20,0.1,0.5,2, ",
     &    " -1,                   IWUNIT"
      endif
      if(idipca.eq.1) then
        write(iu,'(a)') 
     &    "1,1.0,1D-12,1.0E-5,0,,0,200.0,,,",
     &    "5000,0.02,1D-10,-1,0,0,0,,,,        ",
     &    " 'mom=0',2,'nojj',0,0,0,-1,0,,",
     &    " 'end',,,,,,,,,,",
     &    "  7,20,0.1,0.5,2, ",
     &    " -1,                   IWUNIT"
      endif
      if(imonop.eq.1 .or. imonsn.eq.1 .or. imonni.eq.1
     &   .or. imonca.eq.1) then
        write(iu,'(a)') 
     &    "10,1.0,0.0000001,1.0E-5,0,,0,200.0,,,",              
     &    "20000,0.02,0.000001,-1,0,0,0,,,,        ",
     &    " 'mom=0',1,'nojj',0,-12,0,-1,0,,",
     &    " 'end',,,,,,,,,,",
     &    "  7,20,0.1,0.5,2, ",
     &    " -1,                   IWUNIT"
!        write(iu,'(a)') 
!     &    "0,1.0,0.0000001,1.0E-5,0,,0,200.0,,,",              
!     &    "2000,0.02,0.000001,-1,0,0,0,,,,        ",
!     &    " 'mom=0',2,'nojj',0,0,0,-1,0,,",
!     &    " 'end',,,,,,,,,,",
!     &    "  7,20,0.1,0.5,2, ",
!     &    " -1,                   IWUNIT"
      endif
      if(iquadp.eq.1) then
!        write(iu,'(a)') 
!     &    "2,1.0,0.0000001,1.0E-5,0,,0,200.0,,,",
!        write(iu,'(a)') 
!     &    "2,1.0,0.0000001,1.0E-5,0,,0,200.0,,,",
!     &    "2000,0.02,0.000001,-1,0,0,0,,,,        ",
!     &    " 'mom=0',1,'nojj',0,0,-8,-1,0,,",
!     &    " 'end',,,,,,,,,,",
!     &    "  7,20,0.1,0.5,2, ",
!     &    " -1,                   IWUNIT"
       write(iu,'(a)') 
     &    "2,1.0,0.0000001,1.0E-5,0,,0,200.0,,,",
     &    "2000,0.02,0.000001,-1,0,0,0,,,,        ",
     &    " 'mom=0',2,'nojj',0,0,0,-1,0,,",
     &    " 'end',,,,,,,,,,",
     &    "  7,20,0.1,0.5,2, ",
     &    " -1,                   IWUNIT"
       endif
c
      close(iu)
c
      return
      end
c-----resenergy--------------------------------------------------------
c
      real*8 function resenergy()
c
c     Reads result file 'lsmooth.inp' and extracts resonance energy.
c      igdrpb     = GDR in 208Pb
c      igqrpb     = GQR in 208Pb
c      igdro      = GDR in 16O
c      idipol     = dipole polarizability in 208Pb
c      idipsn     = dipole polarizability in 132Sn
c      idipca     = dipole polarizability in 48Ca
c      imonop     = monopole polarizability in 208Pb
c
      include 'params.inc'
      include 'constn.inc'
      include 'extobs.inc'
c
      data iu/31/               ! unit for input
c
c-----------------------------------------------------------------------
c
ctest      write(6,'(a,i4)') 'RESENERGY called: igdrpb=',igdrpb
ctest      write(2,'(a,i4)') 'RESENERGY called: igdrpb=',igdrpb
c
      if(igmrpb+igdrpb+igqrpb+igdro+idipol+idipsn+idipca+idipzr
     &   +idipni+ipygmy+ipygsn+imonop+imonsn+imonni+imonca+iquadp
     &   .ne.1)  stop "RESENERG: no valid option "
c
c     polarizability
c
      if(idipol.eq.1 .or. idipsn.eq.1 .or. idipca.eq.1 .or. idipzr.eq.1
     &    .or. idipni.eq.1 .or. imonop.eq.1 .or. imonsn.eq.1 
     &    .or. imonni.eq.1 .or. imonca.eq.1 .or. iquadp.eq.1) then
c        call system('tail -n 4 for014 | head -n 1 > xxpolar')
        open(iu,file='xxpolar')
c        read(iu,'(14x,3(9x,g13.5))') y1,y2,respolariz
        read(iu,*) i1,i2,y1,y2,respolariz,resomeg,resbel,polar
        close(iu)
        resenergy = polar !resomeg ! respolariz*(8D0*PI/3D0)
!        WRITE(2,'(a,2i5,5(1pg13.5))') 'XXPOLAR:',
!     &           i1,i2,y1,y2,respolariz,resomeg,resbel
!        STOP
        IF(imonop==1) resenergy= resomeg   ! resbel
        return
      endif

c
c     all GR
c
      open(iu,file='lsmooth.inp')
c
      if(igdrpb.eq.1 .or. igqrpb.eq.1 .or. igdro.eq.1) then
        do i=1,4
          read(iu,'(1x)',err=99,end=99)
        enddo
        read(iu,*,err=99,end=99) idum,resin
        resenergy = resin
      elseif(igmrpb.eq.1) then
        do i=1,3
          read(iu,'(1x)',err=99,end=99)
        enddo
        read(iu,*,err=99,end=99) idum,resin
        resenergy = resin
      endif
      close(iu)
c
c     pygmy resonance
c
      if(ipygmy.eq.1 .or. ipygsn.eq.1) then
        open(iu,file='ewsrlow.res')
        read(iu,*) resenergy
        close(iu)
      endif
c
ctest      write(6,'(a)') 'RESENERGY finished'
ctest      write(2,'(a)') 'RESENERGY finished'
c
      return
c
99    stop 'RESENERGY: premature end of LSMOOTH.INP'
      end


c-----extobsmf--------------------------------------------------------
c
      subroutine extobsmf(par,workdv,nparmx)
c
c     Extrapolate to special observables in standard nuclei 
c     including extrapolation error using
c     matrix 'zarray' which is transferred via common /fitpar/ .
c
c     The following observables are computed (alternatively):
c      iskinpb     = neutron rsm skin in 208Pb
c      ilevelsn    = neutron level sequence in 132Sn
c      ishe        = ??
c      iqalpha     = Q_alpha for 120/182
c      is2n        = S_2n for 120/182
c      is2p        = S_2p for 120/182
c      iisotop     = isotope shift 214Pb-208Pb
c      iskinsn74   = neutron rms skin in 124Sn
c      iskinsn82   = neutron rms skin in 132Sn
c      iskinsn90   = neutron rms skin in 140Sn
c      iskinca     = neutron rms skin in 48Ca
c      iskinni     = neutron rms skin in 68Ni
c      iskindf     = neutron diffraction skin in 208Pb
c      irmsneut    = neutron rms radius in 208Pb
c      irdfneut    = neutron diffraction radius in 208Pb
c      ifrmneut    = weak-charge formfactor at given Q in 208Pb
c                    (formerly neutron formfactor at given Q in 208Pb)
c      ibindexo    = binding energy of very neutron rich 148Sn
c      irdfcoul    = charge diffraction radius in 208Pb
c      isigcoul    = charge surface thickness in 208Pb
c      irmscoul    = charge r.m.s. radius
c
c
      include 'params.inc'
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'comfit.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'extobs.inc'
      include 'constn.inc'
c
c     'workdv' is field for derivatives with respect to parameter-variation.
c
      dimension   par(kparmx),workdv(kparmx,1)
      parameter (qtestneut=0.3978D0)   ! Q [fm^-1] for testing neutron formfactor 208Pb
!      parameter (qtestneut=0.475D0)   ! Q [fm^-1] for testing neutron formfactor
!      parameter (qtestneut=0.778D0)   ! Q [fm^-1] for testing neutron formfactor
c
c
c-----------------------------------------------------------------------
c
!      WRITE(2,'(a)') ' extobsmf: '
      itest = iskinpb-iskinpb*iskinpb
     &       + ilevelsn-ilevelsn*ilevelsn
     &       + ishe-ishe*ishe
     &       + iqalpha-iqalpha*iqalpha
     &       + is2n-is2n*is2n
     &       + is2p-is2p*is2p
     &       + iskinsn74-iskinsn74*iskinsn74
     &       + iskinsn82-iskinsn82*iskinsn82
     &       + iskinsn90-iskinsn90*iskinsn90
     &       + iskinca-iskinca*iskinca
     &       + iskinni-iskinni*iskinni
     &       + iskindf-iskindf*iskindf
     &       + irmsneut-irmsneut*irmsneut
     &       + irdfneut-irdfneut*irdfneut
     &       + ifrmneut-ifrmneut*ifrmneut
     &       + ibindexo-ibindexo*ibindexo
     &       + irdfcoul-irdfcoul*irdfcoul
     &       + isigcoul-isigcoul*isigcoul
     &       + irmscoul-irmscoul*irmscoul
      if(itest.ne.0) stop 'EXTOBS: options not 0 or 1'
      if(iskinpb+ilevelsn+ishe+iqalpha+is2n+is2p+iisotop
     &   +iskinsn74+iskinsn82+iskinsn90+iskinca+iskinni
     &   +iskindf+irmsneut+irdfneut+ifrmneut+ibindexo
     &   +irdfcoul+isigcoul+irmscoul
     &   .ne.1) then
        WRITE(2,'(43i2)') nextobs
        CALL flush(2)
        stop "EXTOBSMF: only one option can be active"
      endif
c
c     loop over variation of parameters
c
      do i=0,nparmx
        if(i.gt.0) par(i) = par(i) + dparmo(i)
        call attach(par,jp)
c
        if(iskinpb.eq.1 .or. iskindf.eq.1 .or. irmsneut.eq.1
     &     .or. irdfneut.eq.1 .or. ifrmneut.eq.1 
     &     .or. irdfcoul.eq.1 .or. isigcoul.eq.1
     &     .or. irmscoul.eq.1) then
          nprot  = 82
          nneut  = 126
!          if(ifrmneut==1) then
!            nprot=20
!            nneut=48
!          endif
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          if(iskinpb.eq.1) observ = rmsn-rmsp
          if(iskindf.eq.1) observ = rdmsn-rdmsp
          if(irmsneut.eq.1) observ = rmsn
          if(irdfneut.eq.1) observ = rdmsn
          if(ifrmneut.eq.1) then
            observ = fbint(qtestneut,fow)
          endif
          if(irdfcoul.eq.1) observ = rdmsc
          if(isigcoul.eq.1) observ = surfc
          if(irmscoul.eq.1) observ = rmsc
        elseif(iskinca.eq.1) then
          nprot  = 20
          nneut  = 28
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = rmsn -rmsp
        elseif(iskinni.eq.1) then
          nprot  = 28
          nneut  = 40
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = rmsn-rmsp
        elseif(iskinsn74.eq.1) then
          nprot  = 50
          nneut  = 74
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = rmsn-rmsp
        elseif(iskinsn82.eq.1) then
          nprot  = 50
          nneut  = 82
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = rmsn-rmsp
        elseif(iskinsn90.eq.1) then
          nprot  = 50
          nneut  = 90
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = rmsn-rmsp
        elseif(ibindexo.eq.1) then
          nprot  = 50
          nneut  = 98
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = etotal
        elseif(ilevelsn.eq.1) then
          nprot  = 50
          nneut  = 82
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = spenrg(nmaxp+16)-spenrg(nmaxp+14)
        elseif(ishe.eq.1) then
          nprot  = 120
          nneut  = 182
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = etotal
        elseif(iqalpha.eq.1) then
          nprot  = 120
          nneut  = 182
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = etotal
c
          nprot  = 118
          nneut  = 180
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = 28.3-etotal+observ
        elseif(is2n.eq.1) then
          nprot  = 120
          nneut  = 182
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = etotal
c
          nprot  = 120
          nneut  = 180
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = etotal-observ
        elseif(is2p.eq.1) then
          nprot  = 120
          nneut  = 182
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = etotal
c
          nprot  = 118
          nneut  = 182
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = etotal-observ
        elseif(iisotop.eq.1) then
          nprot  = 82
          nneut  = 132
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = rmsc**2
          nprot  = 82
          nneut  = 126
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = rmsc**2-observ
        endif
cold      ngrid  = 0
c
c       collect array of deviations
c
        if(i.eq.0) then
          observsav = observ
        else
          workdv(i,1) = (observ-observsav)/dparmo(i)
        endif
c
        if(i.gt.0) par(i) = par(i) - dparmo(i)
      enddo
c
c     accumulate extrapolation errors by folding with 'zarray'.
c
      errobs   = 0.0
c
      if(nparmx.gt.0) then
        do i=1,nparmx
        do j=1,nparmx
          errobs  = workdv(i,1)*workdv(j,1)*zarray(i,j)
     &              + errobs
        enddo
        enddo
        errobs  = sqrt(errobs)
      endif
c
c     final print
c
        if(iskinpb.eq.1) then
          write(02,'(a,2f8.3)')
     &     'rms neutron skin 208Pb: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') ' skin 208Pb'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(iskinca.eq.1) then
          write(02,'(a,2f8.3)')
     &     'rms neutron skin 48Ca: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') ' skin 48Ca'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(iskinni.eq.1) then
          write(02,'(a,2f8.3)')
     &     'rms neutron skin 68Ni: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') ' skin 68Ni'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(iskinsn74.eq.1) then
          write(02,'(a,2f8.3)')
     &     'rms neutron skin 124Sn: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') ' skin 124Sn'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(iskinsn82.eq.1) then
          write(02,'(a,2f8.3)')
     &     'rms neutron skin 132Sn: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') ' skin 132Sn'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(iskinsn90.eq.1) then
          write(02,'(a,2f8.3)')
     &     'rms neutron skin 140Sn: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') ' skin 140Sn'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(ilevelsn.eq.1) then
          write(02,'(a,2f8.3)')
     &     'neutron level diff. 132Sn: value,uncertainty=',
     &     observsav,errobs
        elseif(ishe.eq.1) then
          write(02,'(a,f10.3,f8.3)')
     &     'binding energy 120/182: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') 'EB(120/182)'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(iqalpha.eq.1) then
          write(02,'(a,2f8.3)')
     &     'Q_alpha 120/182: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') 'Qalpha(120/182)'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(is2n.eq.1) then
          write(02,'(a,2f8.3)')
     &     'S_2n 120/182: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') 'S2n(120/182)'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(is2p.eq.1) then
          write(02,'(a,2f8.3)')
     &     'S_2p 120/182: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') 'S2p(120/182)'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(iisotop.eq.1) then
          write(02,'(a,2f8.3)')
     &     'iso-shift 208Pb: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') 'r(214)-r(208)'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(iskindf.eq.1) then
          write(02,'(a,2f8.3)')
     &     'diff.skin 208Pb: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') 'dskin 208Pb'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(irmsneut.eq.1) then
          write(02,'(a,2f8.3)')
     &     'neutron rms: value,uncertainty=',
     &     observsav,errobs
!          write(2,'(a)') ' rmsn 208Pb'
!          write(2,'(1pg24.16)')  observsav
!          write(2,'(100(1pg24.16))') 
!     &      (workdv(i,1),i=1,nparmx)
          write(81,'(a)') ' rmsn 208Pb'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(irdfneut.eq.1) then
          write(02,'(a,2f8.3)')
     &     'neutron diff.rad.: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') 'rdmsn 208Pb'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(ifrmneut.eq.1) then
          write(02,'(a,2f8.3)')
     &     'weak-charge F_W(0.475): value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') 'formw 208Pb'
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(ibindexo.eq.1) then
          write(02,'(a,f12.3,f8.3)')
     &     'binding energy 148Sn: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') ' E_B 148Sn '
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(irdfcoul.eq.1) then
          write(02,'(a,f12.3,f8.3)')
     &     ' charge R_diffr 208Pb: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') ' ch.R_diffr 208Pb '
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(isigcoul.eq.1) then
          write(02,'(a,f12.3,f8.3)')
     &     ' charge surfac.th. 208Pb: value,uncertainty=',
     &     observsav,errobs
          write(81,'(a)') ' ch.surf.th 208Pb '
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        elseif(irmscoul.eq.1) then
          write(02,'(a,f12.3,f8.3)')
     &     ' charge rms radius 208Pb: value,uncertainty=',
     &     observsav,errobs
!          write(2,'(a)') ' charge r.m.s. 208Pb '
!          write(2,'(1pg24.16)')  observsav
!          write(2,'(100(1pg24.16))') 
!     &      (workdv(i,1),i=1,nparmx)
          write(81,'(a)') ' charge r.m.s. 208Pb '
          write(81,'(1pg24.16)')  observsav
          write(81,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
        endif
c
      return
      end

c-----extobsFq--------------------------------------------------------
c
      subroutine extobsFq(par,workdv,nparmx)
c
c     Extrapolate to weak-charge formfactor
c     including extrapolation error using
c     matrix 'zarray' which is transferred via common /fitpar/ .
c
c     The q-values at which the errors are computed are given
c     here by data statements.
c
c
      include 'params.inc'
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'comfit.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'extobs.inc'
      include 'constn.inc'
c
c     'workdv' is field for derivatives with respect to parameter-variation.
c
      dimension   par(kparmx),workdv(kparmx,1)
      parameter (qwmax=2D0,nqwmax=40)
      double precision formfw(0:kparmx,nqwmax)
c
c-----------------------------------------------------------------------
c
!      WRITE(2,'(a)') ' extobsFw: '
c
c     loop over variation of parameters
c
      open(281,file='error-Fw')
      qwstep = qwmax/nqwmax
      do i=0,nparmx
        if(i.gt.0) par(i) = par(i) + dparmo(i)
        call attach(par,jp)
c
!        nprot  =  20
!        nneut  =  28
!        nnumb  = nprot+nneut
!        anumb  = nnumb
!        nmaxp  = 0
!        nmaxn  = 0
!        call hafosu
!        write(2,'(a)') '# F_W for 48Ca','#  i  q  F_W'
!        do iqw=1,nqwmax
!          write(2,'(i5,f10.5,1pg13.5)') 
!     &          iqw,qwstep*iqw,fbint(qwstep*iqw,fow)
!        enddo
        nprot  =  82
        nneut  = 126
        nnumb  = nprot+nneut
        anumb  = nnumb
        nmaxp  = 0
        nmaxn  = 0
        call hafosu
!        write(2,'(a)') '# F_W for 208Pb','#  i  q  F_W'
        do iqw=1,nqwmax
          formfw(i,iqw) = fbint(qwstep*iqw,fow)
!          write(2,'(i5,f10.5,1pg13.5)') 
!     &          iqw,qwstep*iqw,fbint(qwstep*iqw,fow)
        enddo
        if(i.gt.0) par(i) = par(i) - dparmo(i)
      enddo
cold      ngrid  = 0
c
c     loop over all momenta and compute error relations
c
      do iqw=1,nqwmax
c
c       collect array of deviations
c
        do i=1,nparmx
            workdv(i,1) = (formfw(i,iqw)-formfw(0,iqw))/dparmo(i)
        enddo
c
c       accumulate extrapolation errors by folding with 'zarray'.
c
        observsav = formfw(0,iqw)
        errobs   = 0.0
        if(nparmx.gt.0) then
          do i=1,nparmx
          do j=1,nparmx
            errobs  = workdv(i,1)*workdv(j,1)*zarray(i,j)
     &                + errobs
          enddo
          enddo
          errobs  = sqrt(errobs)
        endif
c
c     final print
c
        write(02,'(a,f5.2,a,2f8.3)')
     &     'F_W(q=',iqw*qwstep,') value,uncertainty=',
     &     observsav,errobs
        write(281,'(a,f4.2,a)') 'F_W(',iqw*qwstep,')'
        write(281,'(1pg24.16)')  observsav
        write(281,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
c
      enddo
      close(281)
c
      return
      end

c-----extbarr--------------------------------------------------------
c
      subroutine extbarr(par,nparmx)
c
c     Prepare input files for extrapolation to fission barrier
c     in Z=120,N=182,
c     Matrix 'zarray' which is transferred via common /fitpar/ .
c
      include 'params.inc'
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'comfit.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'constn.inc'
c
c     'workdv' is field for derivatives with respect to parameter-variation.
c
      dimension   par(kparmx)
c
      character*3 ident
c
c-----------------------------------------------------------------------
c
c     'dparmo' and 'zarray' in data file for later use.
c
      open(31,file='zarray.dat')
      write(31,'(i5)') nparmx
      write(31,'(17g24.16)') (dparmo(j),j=1,nparmx)
      do i=1,nparmx
        write(31,'(17g24.16)') (zarray(j,i),j=1,nparmx)
      enddo
      close (31)
c
c     prepare input for 'skyax.exe'
c
      open(32,file='allsky.bat')
      do i=0,nparmx
        if(i.gt.0) par(i) = par(i) + dparmo(i)
        call attach(par,jp)
c
c       prepare identifier
c
        open(31,status='scratch')
        write(31,'(i3)') i+100
        rewind(31)
        read(31,'(a3)') ident
        close(31)
        ident(1:1) = 's' 
c
c       write for005.*
c
        nprot = 120
        nneut = 182
c                                         for the ground state
        open(31,file='skyinpgs.'//ident)
        write(31,'(a/a/2i5,a,2(/4(1pg17.9),a),(/5(1pg17.9),a)
     &             /i4,0pf12.8,i4,f6.1,i4,a)')
     &   " 100000          ITMXT   max. total iterations",
     &   "'"//ident//":SHE '",
     &   nprot,nneut,"  nprot,nneut",
     &   t0,t1,t2,t3,"t0,t1,t2,t3",
     &   x0,x1,x2,x3,"  x0,x1,x2,x3",
     &   t4/2,b4p,power,amfc(1),amfc(2),
     &         " b4,b4p,alpha,amfc",
     &   ifex,riftm,ifzpe,acoul,ioptsk," if:ex,ls,cm,CHG"
        if(ioptsk.eq.11) then
          write(31,'(3(1pg17.9),a/3(1pg17.9),a)')
     &      alpinv,betinv,riftm1,'    alpinv,betinv,riftm1',
     &      b7,b7p,pow7, '    b7,b7p,pow7'
        endif
        write(31,'(i4,4(1pg17.9),a)')
     &   ipair,pairfp,pairfn,ecstab,rho0pr,
     &      "  ipair,pairfp,pairfn"
        write(31,'(a)')
     &   "    1 0 0  1.0  1.0    0    0    0    "//
     &   "       grid: SYM,NR,NZ,DR,DZ,ICUTR,ICUTZ,IDRTYP",
     &   "  0.0 0.0 1.9 1.9   31   31 0   "//
     &   "       init: BETAIN, OCTINI,  DEOCC(2),MMAX(2)",
     &   "  0.0 0.00 0.35  0.025   2  0 0  "//
     &   "       CHF: swicon,betmin,betmax,betaw,lqcon,ifresh",
     &   "  0.0  0.0  1.0  "//
     &   "       thermal: tempmn,tempmx,tempw",
     &   "  4000  4000  800    0   "//
     &   "       print opt.: IPRINT,IPLOT,ITMAX,ITEST",
     &   "   100.00    0.40  0.20  1.0  0.30  "//
     &   "       step iter:E0INV,X0INV,ADDrho,adpot,ADDOCC ",
     &   " 1.E-7,  1.0,                  "//
     &   "       limit iter:ERRENR,ERRVAR,",
     &   " 9999, 1,1.5,                    ",
     &   " 1.E-6,1.E-6,3000,3000,100,5, ",
     &   " 0.5,1.0, -2.1,               ",
     &   " 1,1,1,1,                   ",
     &   " 1.D-5,20,6.0D-1,           ",
     &   " 1.D-7 ,-50.0,   20,   1.0,",
     &   " 0.1E-3, 0.8, 0.0,  0.3,  999, 1,     "//
     &   "      CHF iter:QERR,C0,D0,QEPS,ITCONS(999),idampc",
     &   " 0,  10.0,  15.0,  0.03,   "//
     &   "      stab.CHF:ifstab,stwid,stslop,rhocut",
     &   " 3.5,0.30,    0,   0.05,     "//
     &   "      damp.CHF:DALP,DGAM,IFLIN,SWIFAC,",
     &   " 0,   ",
     &   "  -1,,,,,,,,,,,,,,,   ",
     &   "-1, "
        close(31)
c                                         for the barrier
        open(31,file='skyinpbar.'//ident)
        write(31,'(a/a/2i5,a,2(/4(1pg17.9),a),(/5(1pg17.9),a)
     &             /i4,0pf12.8,i4,f6.1,i4,a)')
     &   " 100000          ITMXT   max. total iterations",
     &   "'"//ident//":SHE '",
     &   nprot,nneut,"  nprot,nneut",
     &   t0,t1,t2,t3,"t0,t1,t2,t3",
     &   x0,x1,x2,x3,"  x0,x1,x2,x3",
     &   t4/2,b4p,power,amfc(1),amfc(2),
     &         " b4,b4p,alpha,amfc",
     &   ifex,riftm,ifzpe,acoul,ioptsk," if:ex,ls,cm,CHG"
        if(ioptsk.eq.11) then
          write(31,'(3(1pg17.9),a/3(1pg17.9),a)')
     &      alpinv,betinv,riftm1,'    alpinv,betinv,riftm1',
     &      b7,b7p,pow7, '    b7,b7p,pow7'
        endif
        write(31,'(i4,4(1pg17.9),a)')
     &   ipair,pairfp,pairfn,ecstab,rho0pr,
     &      "  ipair,pairfp,pairfn"
        write(31,'(a)')
     &   "    1 0 0  1.0  1.0    0    0    0    "//
     &   "       grid: SYM,NR,NZ,DR,DZ,ICUTR,ICUTZ,IDRTYP",
     &   "  0.0 0.0 1.9 1.9   31   31  0  "//
     &   "       init: BETAIN, OCTINI,  DEOCC(2),MMAX(2)",
     &   "  1.0 0.225 0.35  0.025   2  0 0  "//
     &   "       CHF: swicon,betmin,betmax,betaw,lqcon,ifresh",
     &   "  0.0  0.0  1.0  "//
     &   "       thermal: tempmn,tempmx,tempw",
     &   "  4000  4000  800    0   "//
     &   "       print opt.: IPRINT,IPLOT,ITMAX,ITEST",
     &   "   100.00    0.40  0.20  1.0  0.30  "//
     &   "       step iter:E0INV,X0INV,ADDrho,adpot,ADDOCC ",
     &   " 1.E-7,  1.0,                  "//
     &   "       limit iter:ERRENR,ERRVAR,",
     &   " 9999, 1,1.5,                    ",
     &   " 1.E-6,1.E-6,3000,3000,100,5, ",
     &   " 0.5,1.0, -2.1,               ",
     &   " 1,1,1,1,                   ",
     &   " 1.D-5,20,6.0D-1,           ",
     &   " 1.D-7 ,-50.0,   20,   1.0,",
     &   " 0.1E-3, 0.8, 0.0,  0.3,  999, 1,     "//
     &   "      CHF iter:QERR,C0,D0,QEPS,ITCONS(999),idampc",
     &   " 0,  10.0,  15.0,  0.03,   "//
     &   "      stab.CHF:ifstab,stwid,stslop,rhocut",
     &   " 3.5,0.30,    0,   0.05,     "//
     &   "      damp.CHF:DALP,DGAM,IFLIN,SWIFAC,",
     &   " 0,   ",
     &   "  -1,,,,,,,,,,,,,,,   ",
     &   "-1, "
        close(31)
c                                             for g.s. of 264Hs
        open(31,file='skyinphs.'//ident)
        write(31,'(a/a/2i5,a,2(/4(1pg17.9),a),(/5(1pg17.9),a)
     &             /i4,0pf12.8,i4,f6.1,i4,a)')
     &   " 100000          ITMXT   max. total iterations",
     &   "'"//ident//":SHE '",
     &   108,156,"  nprot,nneut",
     &   t0,t1,t2,t3,"t0,t1,t2,t3",
     &   x0,x1,x2,x3,"  x0,x1,x2,x3",
     &   t4/2,b4p,power,amfc(1),amfc(2),
     &         " b4,b4p,alpha,amfc",
     &   ifex,riftm,ifzpe,acoul,ioptsk," if:ex,ls,cm,CHG"
        if(ioptsk.eq.11) then
          write(31,'(3(1pg17.9),a/3(1pg17.9),a)')
     &      alpinv,betinv,riftm1,'    alpinv,betinv,riftm1',
     &      b7,b7p,pow7, '    b7,b7p,pow7'
        endif
        write(31,'(i4,4(1pg17.9),a)')
     &   ipair,pairfp,pairfn,ecstab,rho0pr,
     &      "  ipair,pairfp,pairfn"
        write(31,'(a)')
     &   "    1 0 0  1.0  1.0    0    0    0    "//
     &   "       grid: SYM,NR,NZ,DR,DZ,ICUTR,ICUTZ,IDRTYP",
     &   "  0.22 0.0 1.9 1.9   31   31  0   "//
     &   "       init: BETAIN, OCTINI,  DEOCC(2),MMAX(2)",
     &   "  0.0 0.0 0.35  0.025   2  0 0  "//
     &   "       CHF: swicon,betmin,betmax,betaw,lqcon,ifresh",
     &   "  0.0  0.0  1.0  "//
     &   "       thermal: tempmn,tempmx,tempw",
     &   "  4000  4000  800    0   "//
     &   "       print opt.: IPRINT,IPLOT,ITMAX,ITEST",
     &   "   100.00    0.40  0.20  1.0  0.30  "//
     &   "       step iter:E0INV,X0INV,ADDrho,adpot,ADDOCC ",
     &   " 1.E-7,  1.0,                  "//
     &   "       limit iter:ERRENR,ERRVAR,",
     &   " 9999, 1,1.5,                    ",
     &   " 1.E-6,1.E-6,3000,3000,100,5, ",
     &   " 0.5,1.0, -2.1,               ",
     &   " 1,1,1,1,                   ",
     &   " 1.D-5,20,6.0D-1,           ",
     &   " 1.D-7 ,-50.0,   20,   1.0,",
     &   " 0.1E-3, 0.8, 0.0,  0.3,  999, 1,     "//
     &   "      CHF iter:QERR,C0,D0,QEPS,ITCONS(999),idampc",
     &   " 0,  10.0,  15.0,  0.03,   "//
     &   "      stab.CHF:ifstab,stwid,stslop,rhocut",
     &   " 3.5,0.30,    0,   0.05,     "//
     &   "      damp.CHF:DALP,DGAM,IFLIN,SWIFAC,",
     &   " 0,   ",
     &   "  -1,,,,,,,,,,,,,,,   ",
     &   "-1, "
        close(31)
        if(i.gt.0) par(i) = par(i) - dparmo(i)
        write(32,'(a)')
     &   'skyax.exe < skyinpgs.'//ident//' > /dev/null',
     &   'skyax.exe < skyinpbar.'//ident//' > /dev/null',
     &   'skyax.exe < skyinphs.'//ident//' > /dev/null'
      enddo
      write(32,'(a)')
     &  'rm -f twoQP*',
     &  'rm -f for007*',
     &  'rm -f local*',
     &  'rm -f densi*',
     &  'rm -f conver*',
     &  'rm -f quadru*',
     &  'rm -f fields*'
      close(32)
c
c     hafo input for computation of surface energy
c
      do i=0,nparmx
        if(i.gt.0) par(i) = par(i) + dparmo(i)
        call attach(par,jp)
        open(31,status='scratch')
        write(31,'(i3)') i+100
        rewind(31)
        read(31,'(a3)') ident
        close(31)
        ident(1:1) = 's' 
        open(31,file='for005.'//ident)
c
c     the mean field part
c
        write(31,'(a)')
     &   '1,,    IWUNIT',
     &   "' 1D SHF input'",
     &   '0,0.3,              NGRID,RSTEP',
     &   '-1,1E-8,1E-9,278,5780,3,-1, IWR,DABR/E,ITPR/MAX,IFPL',
     &   '0.2,90.,0.2,-1,,               X0DAMP,E0DAMP,ADDNEW,WALL'
        write(31,'(9(i2,a),i2,a,5x,a)') 
     &    ifex,',',iftm,',',ifenrg,',',ifrso,',',ifzpe,',',
     &    ifdfms,',',ifrhoc,',',ifenav,',',imfold,',',ioptsk,
     &    ', if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk'
        write(31,'(a)')
     &   '1,,,                       ifform',
     &   '0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12',
     &   '0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12'
        IF(ioptsk==13) THEN            ! Fayans functional
          WRITE(31,'(3(1pg18.10),a)') avp,h1vp,h2vp,' avp,h1vp,h2vp'
          WRITE(31,'(5(1pg18.10))') avm,h1vm,h2vm
          WRITE(31,'(5(1pg18.10))') sigmavol,rho0fy
          WRITE(31,'(5(1pg18.10))') asp,hsp,hgradp 
          WRITE(31,'(5(1pg18.10))') asm,hsm,hgradm 
          WRITE(31,'(5(1pg18.10))') hcoul,sigmacoul
          WRITE(31,'(5(1pg18.10))') fxi,h1xi,hgradxi,gammaxi,rho0xi
          WRITE(31,'(2(1pg18.10),i2)') delfxi,ecstab,isovectormodel
          WRITE(31,'(5(1pg18.10))') c0tau,c1tau 
          WRITE(31,'(5(1pg18.10))') c0nabj,c1nabj,c0tens,c1tens
        ELSE 
          write(31,'(5(1pg20.12,a)/5(g20.12,a))')
     &    t0,',',t1,',',t2,',',t3,',',t4,',      t0,t1,t2,t3,t4',
     &    x0,',',x1,',',x2,',',x3,',',b4p,',       x0,x1,x2,x3,b4p'
        END IF
        if(ioptsk.eq.5 .or. ioptsk.eq.6 .or. ioptsk.eq.11) then
          write(31,'(2(1pg20.12,a)/2(1pg20.12,a))')
     &      c1,',',c4,',                      c1 c4',
     &      c1p,',',c4p,',                    c1p,c4p'
        endif
        if(ioptsk.eq.8 .or. ioptsk.eq.11) then
          write(31,'(3(1pg20.12),a/2(1pg20.12,a))')
     &      b3tilde,tte,tto,',   b3tilde,tte,tto',
     &      b12,',',b12p,',                     b12,b12p'
        endif
        if(ioptsk.eq.9 .or. ioptsk.eq.11) then
          write(31,'(2(1pg20.12,a)/1pg20.12,a)')
     &      b7,',',b7p,',               b7,b7p',
     &      pow7,',               pow7'
        endif
        if(ioptsk.eq.10 .or. ioptsk.eq.11) then
          write(31,'(8(1pg20.12),a)')
     &      riftm,rho0pr,powerp,ecstab,alpinv,betinv,
     &      alpin2,betin2,
     &    ',      riftm,rho0pr,powerp,ecstab,alpinv'
        endif
        if(ioptsk.le.4 .or. ioptsk.eq.7) then
          write(31,'(5(f10.3,a)/4(f10.6,a))')
     &    t5,',',t6,',',t7,',',t8,',',t9,',    t5,t6,t7,t8,t9',
     &    x5,',',x6,',',x7,',',x8,',           x5,x6,x7,x8'
        endif
        IF(ioptsk.NE.13) 
     &     write(31,'((f10.3,a),2(f10.6,a))')
     &      tyuk,',',xyuk,',',ayuk,',           tyuk,xyuk,ayuk'
        write(31,'(f15.11,a,i2,a,f9.3,a,f9.3,a,f9.3,a)')
     &     power,',1.0,',ipair,',',pairfp,',',pairfn,
     &  ',0,,',cutfac,',,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,...'
        if(i.gt.0) par(i) = par(i) - dparmo(i)
        close(31)
      enddo
c
      return
      end


c-----extobschain--------------------------------------------------------
c
      subroutine extobschain(par,workdv,nparmx)
c
c     Extrapolate to chain of very exotic Sn isotopes
c     including extrapolation error using
c     matrix 'zarray' which is transferred via common /fitpar/ .
c
c
c
      include 'params.inc'
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'comfit.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'constn.inc'
c
c     'workdv' is field for derivatives with respect to parameter-variation.
c
      dimension   par(kparmx),workdv(0:kparmx,3)
c
      logical ts2n
      parameter (ts2n=.true.)
c
c
c-----------------------------------------------------------------------
c
c     loop over chain
c
      if(ts2n) then
        write(02,'(a)') 'Chain of exotic Sn isotopes:',
     &                  '  Z   N     S_2n   uncertainty'
      else
        write(02,'(a)') 'Chain of exotic Sn isotopes:',
     &                  '  Z   N   energy   uncertainty'
      endif
      do nn=48,104,2
c
c     loop over variation of parameters
c
        do i=0,nparmx
          if(i.gt.0) par(i) = par(i) + dparmo(i)
          call attach(par,jp)
c
          nmaxp  = 0
          nmaxn  = 0
          nprot = 50
          nneut = nn
          nnumb  = nprot+nneut
          anumb  = nnumb
          call hafosu
          observ = etotal
          obsrad = rmsc
          if(ts2n) then
            call attach(par,jp)
            nmaxp  = 0
            nmaxn  = 0
            nprot = 50
            nneut  = nn+2
            nnumb  = nprot+nneut
            anumb  = nnumb
            call hafosu
            observs = observ-etotal
          endif
cold      ngrid  = 0
c
c       collect array of deviations
c
          if(i.eq.0) then
            workdv(0,1) = observ
            workdv(0,2) = obsrad
            IF(ts2n) workdv(0,3) = observs
          else
            workdv(i,1) = (observ-workdv(0,1))/dparmo(i)
            workdv(i,2) = (obsrad-workdv(0,2))/dparmo(i)
            IF(ts2n) workdv(i,3) = (observs-workdv(0,3))/dparmo(i)
          endif
c
          if(i.gt.0) par(i) = par(i) - dparmo(i)
        enddo
c
c       accumulate extrapolation errors by folding with 'zarray'.
c
        errobs   = 0.0
        errrad   = 0.0
        errobss   = 0.0
c
        if(nparmx.gt.0) then
          do i=1,nparmx
          do j=1,nparmx
            errobs  = workdv(i,1)*workdv(j,1)*zarray(i,j)
     &                + errobs
            errrad  = workdv(i,2)*workdv(j,2)*zarray(i,j)
     &                + errrad
            IF(ts2n) errobss  = workdv(i,3)*workdv(j,3)*zarray(i,j)
     &                         + errobss
          enddo
          enddo
          errobs  = sqrt(errobs)
          errrad  = sqrt(errrad)
          IF(ts2n) errobss  = sqrt(errobss)
        endif
c
c       final print
c
        write(02,'(2i5,f10.3,f6.2,f6.3,f8.5,f7.3,f6.3)') nprot,nneut,
     &    workdv(0,1),errobs,workdv(0,2),errrad,workdv(0,3),errobss
c
      enddo
c
      return
      end
c-----energy_deform--------------------------------------------------------
c
      function energy_deform(nprotin,nneutin)
c
c     Prepare input files for computation of g.s. for 264Hs,
c     compute, and return ZPE-corrected energy on  'ehs'.
c
      include 'params.inc'
c
      include 'numpar.inc'
c      include 'wavfun.inc'
      include 'comfit.inc'
      include 'skypar.inc'
c      include 'result.inc'
      include 'constn.inc'
c
      character*40 zin
      logical tcheap                 ! to replace ZPE_rot by simple estimate
      parameter (tcheap=.true.)   
c
c     'workdv' is field for derivatives with respect to parameter-variation.
c
c-----------------------------------------------------------------------
c
      nprot = nprotin
      nneut = nneutin
c                                         for the ground state
      open(31,file='skyinphs.tst')
      write(31,'(a/a/2i5,a,2(/4(1pg17.9),a),(/5(1pg17.9),a)
     &           /i4,0pf12.8,i4,f6.1,i4,a)')
     & " 100000          ITMXT   max. total iterations",
     & 'tst: 264Hs ',
     & nprot,nneut,"  nprot,nneut",
     & t0,t1,t2,t3,"t0,t1,t2,t3",
     & x0,x1,x2,x3,"  x0,x1,x2,x3",
     & t4/2,b4p,power,amfc(1),amfc(2),
     &       " b4,b4p,alpha,amfc",
     & ifex,riftm,ifzpe,acoul,ioptsk," if:ex,ls,cm,CHG,ioptsk"
      if(ioptsk.eq.11) then                            !new
        write(31,'(3(1pg18.10))')  alpinv,betinv,riftm1
        if(abs(b7p)+abs(b7)>0D0 .AND. b7 == 0D0) then
          b7print = 1D-7
        else
          b7print = b7
        endif
        write(31,'(3(1pg18.10))')  b7print,b7p,pow7
      endif 
       write(31,'(i4,4(1pg17.9),a)')
     & ipair,pairfp,pairfn,0,rho0pr,"  ipair,pairfp,pairfn"
      write(31,'(a)')
     & "    1 0 0  1.0  1.0    0    0    0    "//
     & "       grid: SYM,NR,NZ,DR,DZ,ICUTR,ICUTZ,IDRTYP"
      if(tcheap) then
        write(31,'(a)')
     & "  0.2 0.0 1.9 1.9   31   31  0  "//
     & "       init: BETAIN, OCTINI,  DEOCC(2),MMAX(2)",
     &   "  0.0 0.00 0.35  0.025   2  0 1  "//
     & "       CHF: swicon,betmin,betmax,betaw,lqcon,ifresh"
      else
        write(31,'(a)')
     & "  0.2 0.0 2.9 2.9   31   31    "//
     & "       init: BETAIN, OCTINI,  DEOCC(2),MMAX(2)",
     &   "  0.0 0.00 0.35  0.025   2  0 1  "//
     & "       CHF: swicon,betmin,betmax,betaw,lqcon,ifresh"
      endif
      write(31,'(a)')
     & "  0.0  0.0  1.0  "//
     & "       thermal: tempmn,tempmx,tempw",
     & "  4000  4000  600    0   "//
     & "       print opt.: IPRINT,IPLOT,ITMAX,ITEST",
     & "   100.00    0.40  0.20  1.0  0.30  "//
     & "       step iter:E0INV,X0INV,ADDrho,adpot,ADDOCC ",
     & " 4.E-7,  1.0,                  "//
     & "       limit iter:ERRENR,ERRVAR,",
     & " 9999, 1,1.5,                    ",
     & " 1.E-6,1.E-6,3000,3000,100,5, ",
     & " 0.5,1.0, -2.1,               ",
     & " 1,1,1,1,                   ",
     & " 1.D-5,20,6.0D-1,           ",
     & " 1.D-7 ,-50.0,   20,   1.0,",
     & " 0.1E-3, 0.8, 0.0,  0.3,  999, 1,     "//
     & "      CHF iter:QERR,C0,D0,QEPS,ITCONS(999),idampc",
     & " 0,  10.0,  15.0,  0.03,   "//
     & "      stab.CHF:ifstab,stwid,stslop,rhocut",
     & " 3.5,0.30,    0,   0.05,     "//
     & "      damp.CHF:DALP,DGAM,IFLIN,SWIFAC,",
     & " 0,   ",
     & "  -1,,,,,,,,,,,,,,,   ",
     & "-1, "
      close(31)
c
      if(tcheap) then
        call system('rm ensurf.tst')
      else
        call system('rm zpecor.tst')
      endif
      call system('./skyax.exe < skyinphs.tst > /dev/null')
      if(tcheap) then
        open(unit=31,file='ensurf.tst')
        do i=1,8
           read(31,'(1x)',err=98,end=98)
        enddo
        read(31,*) x1,x2,x3,x4,x5,x6,eraw
        close(31)
        zperot = 1.9
      else
        open(unit=31,file='zpecor.tst')
        read(31,'(1x)',err=99,end=99) zin
        read(31,'(1x)',err=99,end=99) zin
        read(31,'(1x)',err=99,end=99) zin
        read(31,'(43x,a)',err=99,end=99) zin
        close(31)
        open(unit=31,status='scratch')
        write(31,'(a)') zin
        write(6,*) zin
        rewind(31)
        read(31,*) zperot,x1,eraw
      close(31)
      endif
      energy_deform = eraw-zperot
c
      return
 98   stop 'FIT-DEFORM: error in reading ENSURF'
 99   stop 'FIT-DEFORM: error in reading ZPECOR'
      end      
