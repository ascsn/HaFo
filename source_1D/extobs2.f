c-----nucvarnat----------------------------------------------------------
c
      subroutine nucvarnat(par,nparmx)
c
c     computes extrapolated errors on nuclear matter properties
c     using the variation in terms of natural parameters 'parmnat'
c
c      par    = array of model parameters
c      dparmo = array of variations of model parameters
c      nparmx = number of varied parameters
c               0 --> no evaluation of extrapolated errors
c
      include 'params.inc'
      include 'constn.inc'
c      parameter (zero=0.0)
c
      include 'comfit.inc'
      include 'skypar.inc'
c
      REAL(8) :: DEOVRA(kparmx),DRHONM(kparmx),DCOMPR(kparmx),
     *          DSYMEN(kparmx),dskapp(kparmx),DDSYME(kparmx),
     &          DD2SYME(kparmx)
      REAL(8) :: DEFMAS(kparmx),DDENEUT(kparmx),DENEUT0(kparmx)
      REAL(8) :: dg0(kparmx),dg0p(kparmx)
      REAL(8) :: par(kparmx),DENEUT12(kparmx)
      REAL(8) :: pars(kparmx),dtto(kparmx),dtte(kparmx)
      REAL(8) :: dgammx(kparmx)    ! misuse for delta_fxi
      REAL(8) :: db4(kparmx),db4p(kparmx),dvpp(kparmx),dvpn(kparmx)
      REAL(8) :: drho0pr(kparmx),db2(kparmx),db2p(kparmx),dhsp(kparmx)
      character*3 id2
      logical tsymderiv
      parameter (tsymderiv=.true.)
      logical tsecord,tstsecord
      parameter (tsecord=.false.,tstsecord=.false.)
      character*1 asec
      REAL(8) ::  work2(kparmx,kparmx,22),work2dum(kparmx,kparmx) 
      REAL(8) ::  work0(22),dwork(kparmx,22)       !  workspace for variation
      REAL(8), ALLOCATABLE :: eneutsav(:),deneutsav(:,:)
      REAL(8), ALLOCATABLE :: esymsav(:),desymsav(:,:)
c
c----------------------------------------------------------------------
c
      IF(nntmax>0) ALLOCATE(eneutsav(nntmax),deneutsav(kparmx,nntmax))
      IF(nsymax>0) ALLOCATE(esymsav(nsymax),desymsav(kparmx,nsymax))
!      IF(nntmax*nsymax==0) RETURN
c
      pars(1:nparmx) = par(1:nparmx)
      IF(tsymderiv) THEN
        nparlow = -nparmx
      ELSE
        nparlow = 0
      END IF
      pars = par
      do i=nparlow,nparmx
      IF(tsecord) THEN
        nparm2=abs(i)
        nparlw2=-abs(i)
      ELSE
        nparm2=0
        nparlw2=0
      END IF
      do j=nparlw2,nparm2
        if(j.ne.0 .and. i.eq.0) CYCLE
        if(j.ne.0 .and. j.eq.-i) CYCLE
        if(j.eq.0 .and. i.eq.0) par(1:nparmx) = pars(1:nparmx)
        if(j.eq.0 .and. i.ne.0) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        if(j.ne.0 .and. j.ne.i) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))/2D0
     &             + SIGN(1,j)*parmnatfac*parmnat(1:nparmx,abs(j))/2D0
        if(j.ne.0 .and. j.eq.i) par(1:nparmx) = pars(1:nparmx)
     &             - SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
c        WRITE(2,'(a,i3,20(1pg13.5))') ' i,params:',i,par(1:nparmx)
        CALL ATTACH(PAR,JP)
        call initfc
        call nucmat(rhonm0,eovra0,compr0,efmas0,f0,f0p,g0,g0p,
     &              symen0,dasym0,d2asy0,skapp0,1,-2)
        d2asy0=h2vm       ! for test
        if(ioptsk.eq.4)  call xkappa(skapp0)
c
c     neutron matter at rho=0.01 fm**-3 and slope thereof
c
        call enernm(0.10d0,zero,zero,eneut0,efmasp,efmasn)
        call enernm(0.1d0+0.0005,zero,zero,eneutp,efmasp,efmasn)
        call enernm(0.1d0-0.0005,zero,zero,eneutm,efmasp,efmasn)
        deneut = (eneutp-eneutm)/(0.001)
        call enernm(0.12d0,zero,zero,eneut12,efmasp,efmasn)
!        eneut12 = compr0-6D0*dasym0             ! override for test
c      call enernm(0.16d0,zero,zero,eneut1,efmasp,efmasn)
c      write(2,'(/a,2g13.5)')
c     &    ' E/A in neutron matter at rho=0.10,slop:',eneut0,deneut
c
c     now store the result onto the proper place
C
        CALL collectw(i,j,eovra0,eovera,deovra,work2(:,:,1),tsymderiv)
        CALL collectw(i,j,rhonm0,rhonm,drhonm,work2(:,:,2),tsymderiv)
        CALL collectw(i,j,compr0,compr,dcompr,work2(:,:,3),tsymderiv)
        CALL collectw(i,j,efmas0,effmas,defmas,work2(:,:,4),tsymderiv)
        CALL collectw(i,j,symen0,symenr,dsymen,work2(:,:,5),tsymderiv)
        CALL collectw(i,j,skapp0,skappa,dskapp,work2(:,:,6),tsymderiv)
        CALL collectw(i,j,dasym0,dsymme,ddsyme,work2(:,:,7),tsymderiv)
        CALL collectw(i,j,d2asy0,d2symme,dd2syme,work2(:,:,22),
     &                tsymderiv)
        CALL collectw(i,j,eneut0,eneut00,deneut0,work2(:,:,8),tsymderiv)
        CALL collectw(i,j,deneut,deneut00,ddeneut,work2(:,:,9),
     &                tsymderiv)
        CALL collectw(i,j,eneut12,eneut120,deneut12,work2(:,:,10),
     &                tsymderiv)
        IF(ioptsk.NE.13) THEN 
!          CALL collectw(i,j,t4/2D0,b40,db4,work2(:,:,11),tsymderiv)
          CALL collectw(i,j,(t4+b4p)/2D0,b40,db4,work2(:,:,11),
     &                  tsymderiv)
          CALL collectw(i,j,b4p,b4p0,db4p,work2(:,:,12),tsymderiv)
          CALL collectw(i,j,pairfp,vpp0,dvpp,work2(:,:,13),tsymderiv)
          CALL collectw(i,j,pairfn,vpn0,dvpn,work2(:,:,14),tsymderiv)
          CALL collectw(i,j,g0,g00,dg0,work2(:,:,20),tsymderiv)
          CALL collectw(i,j,g0p,g0p0,dg0p,work2(:,:,21),tsymderiv)
          IF(ABS(tte)>1D-5) THEN
            CALL collectw(i,j,tte,tte0,dtte,work2(:,:,18),tsymderiv)
            CALL collectw(i,j,tto,tto0,dtto,work2(:,:,19),tsymderiv)
          END IF
          IF(hgradxi.NE.0D0) THEN
            CALL collectw(i,j,hgradxi,rho0pr0,drho0pr,work2(:,:,15),
     &                  tsymderiv)
          ELSE
            CALL collectw(i,j,rho0pr,rho0pr0,drho0pr,work2(:,:,15),
     &                  tsymderiv)
          END IF
!          CALL collectw(i,j,b2,b20,db2,work2(:,:,16),tsymderiv)
          CALL collectw(i,j,b2+0.5*b2p,b20,db2,work2(:,:,16),tsymderiv)
          CALL collectw(i,j,b2p,b2p0,db2p,work2(:,:,17),tsymderiv)
        ELSE
          CALL collectw(i,j,c0nabj,b40,db4,work2(:,:,11),tsymderiv)
          CALL collectw(i,j,c1nabj,b4p0,db4p,work2(:,:,12),tsymderiv)
          CALL collectw(i,j,fxi,vpp0,dvpp,work2(:,:,13),tsymderiv)
          CALL collectw(i,j,h1xi,vpn0,dvpn,work2(:,:,14),tsymderiv)
          CALL collectw(i,j,hgradxi,rho0pr0,drho0pr,work2(:,:,15),
     &                  tsymderiv)
          CALL collectw(i,j,asp,b20,db2,work2(:,:,16),tsymderiv)
          CALL collectw(i,j,hgradp,b2p0,db2p,work2(:,:,17),tsymderiv)
          CALL collectw(i,j,hsp,hsp0,dhsp,work2(:,:,18),tsymderiv)
          CALL collectw(i,j,delfxi,gammx0,dgammx,work2(:,:,19),
     &                  tsymderiv)
          CALL collectw(i,j,delh1xi,work0(20),dwork(:,20),work2(:,:,20),
     &                  tsymderiv)
          CALL collectw(i,j,delhgradxi,work0(21),dwork(:,21),
     &                  work2(:,:,21),tsymderiv)
          CALL collectw(i,j,asm,work0(22),dwork(:,22),work2(:,:,22),
     &                  tsymderiv)
        END IF
c
c       neutron and symmetric EoS
c
        IF(j.EQ.0 .AND. nntmax>0) THEN
          rhostep = rntmax/nntmax
          DO n=1,nntmax
            call enernm(rhostep*n,zero,zero,eneut0,efmasp,efmasn)
            CALL collectw(i,j,eneut0,eneutsav(n),deneutsav(:,n),
     &                    work2dum,tsymderiv)
          END DO
        END IF
        IF(j.EQ.0 .AND. nsymax>0) THEN
          rhostep = rsymax/nsymax
          DO n=1,nsymax
!            call enernm(rhostep*n,half,zero,easym,efmasp,efmasn)
!         insert to compute symmetry energy
             delrho=1D-4
             rhoact=rhostep*n
             call enernm(rhoact,half-delrho,zero,energ2,efmasp,efmasn)
             call enernm(rhoact,half,zero,energ0,efmasp,efmasn)
             call enernm(rhoact,half+delrho,zero,energ1,efmasp,efmasn)
             asymme = (energ1-2.0*energ0+energ2)/(delrho*delrho)/8D0
             easym=asymme
!         end of insert
            CALL collectw(i,j,easym,esymsav(n),desymsav(:,n),
     &                    work2dum,tsymderiv)
          END DO
        END IF
c
        par(1:nparmx) = pars(1:nparmx)
c        WRITE(2,'(20(1pg13.5))') DEOVRA(I) ,DRHONM(I) , 
c     &    DCOMPR(I) , DEFMAS(I) , DSYMEN(I) , DDSYME(I) ,
c     &    dskapp(i) , DENEUT0(i), DDENEUT(i)
C
!        IF(I.GT.0) PAR(1:nparmx) = pars(1:nparmx)
      END DO
      END DO
!     to restore center values
!      call nucmat(rhonm0,eovra0,compr0,efmas0,f0,f0p,g0,g0p,
!     &              symen0,dasym0,d2asy0,skapp0,1,-2)
c
c     write protokol of derivatives (CONTINUED IN EXTOBS)
c
      open(39,status='scratch')
      write(39,'(i4)') nint(parmnatfac*100+1000)
      rewind(39)
      read(39,'(1x,a3)') id2
      close(39)
      open(82,file='erroranalysis.'//id2)
      write(82,'(a,i5,f8.2)') 'natural parameters: nparmx=',
     &     nparmx,parmnatfac
      do n=1,nparmx
        write(82,'(100(1pg24.16))') parmnat(1:nparmx,n)
      enddo
!      write(82,'(a)') 
!     &  ' E/A rhonm  K  m*/m  asym  dasym kappa  E/Aneut slope'
!        write(82,'(100(1pg24.16))') 
!     &    EOVRA0,RHONM0,COMPR0,EFMAS0,SYMEN0,dasym0,skapp0,
!     &    eneut00,deneut00
!      do i=1,nparmx
!        write(82,'(100(1pg24.16))') 
!     &    DEOVRA(I),DRHONM(I),DCOMPR(I),DEFMAS(I),DSYMEN(I),
!     &    DDSYME(I),Dskapp(i),DENEUT0(I),DDENEUT(I)
!      enddo
      IF(tsecord) THEN
        asec='#'
      ELSE
        asec=' '
      END IF
      WRITE(82,'(a)') asec//'E/A-sym'
      WRITE(82,'(1pg24.16)') eovra0
      WRITE(82,'(100(1pg24.16))') deovra(1:nparmx)
      IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,1)
      WRITE(82,'(a)') asec//'rho0-sym'
      WRITE(82,'(1pg24.16)') rhonm0
      WRITE(82,'(100(1pg24.16))') drhonm(1:nparmx)
      IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,2)
      WRITE(82,'(a)') asec//'incompr. K'
      WRITE(82,'(1pg24.16)') compr0
      WRITE(82,'(100(1pg24.16))') dcompr(1:nparmx)
      IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,3)
      WRITE(82,'(a)') asec//'m*/m'
      WRITE(82,'(1pg24.16)') efmas0
      WRITE(82,'(100(1pg24.16))') defmas(1:nparmx)
      IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,4)
      WRITE(82,'(a)') asec//'kappa'
      WRITE(82,'(1pg24.16)') skapp0
      WRITE(82,'(100(1pg24.16))') dskapp(1:nparmx)
      IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,6)
      WRITE(82,'(a)') asec//'J'
      WRITE(82,'(1pg24.16)') symen0
      WRITE(82,'(100(1pg24.16))') dsymen(1:nparmx)
      IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,5)
      WRITE(82,'(a)') asec//'L'
      WRITE(82,'(1pg24.16)') dasym0
      WRITE(82,'(100(1pg24.16))') ddsyme(1:nparmx)
      IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,7)
      WRITE(82,'(a)') asec//'K_sym'
      WRITE(82,'(1pg24.16)') d2asy0
      WRITE(82,'(100(1pg24.16))') dd2syme(1:nparmx)
      IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,22)
      WRITE(82,'(a)') asec//'K_sym-6*L'
      WRITE(82,'(1pg24.16)') d2asy0-6D0*dasym0
      WRITE(82,'(100(1pg24.16))') dd2syme(1:nparmx)-6D0*ddsyme(1:nparmx)
      IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,22)-
     &    6D0*work2(1:nparmx,1:nparmx,7)
      WRITE(82,'(a)') asec//'G0'
      WRITE(82,'(1pg24.16)') g00
      WRITE(82,'(100(1pg24.16))') dg0(1:nparmx)
      WRITE(82,'(a)') asec//'G0P'
      WRITE(82,'(1pg24.16)') g0p0
      WRITE(82,'(100(1pg24.16))') dg0p(1:nparmx)
      IF(ioptsk.EQ.13) THEN
        WRITE(82,'(a)') asec//'c0nabj'
        WRITE(82,'(1pg24.16)') b40
        WRITE(82,'(100(1pg24.16))') db4(1:nparmx)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,11)
        WRITE(82,'(a)') asec//'c1nabj'
        WRITE(82,'(1pg24.16)') b4p0
        WRITE(82,'(100(1pg24.16))') db4p(1:nparmx)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,12)
        WRITE(82,'(a)') asec//'asp'
        WRITE(82,'(1pg24.16)') b20
        WRITE(82,'(100(1pg24.16))') db2(1:nparmx)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,16)
        WRITE(82,'(a)') asec//'hgradp'
        WRITE(82,'(1pg24.16)') b2p0
        WRITE(82,'(100(1pg24.16))') db2p(1:nparmx)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,17)
        WRITE(82,'(a)') asec//'hsp'
        WRITE(82,'(1pg24.16)') hsp0
        WRITE(82,'(100(1pg24.16))') dhsp(1:nparmx)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,18)
        WRITE(82,'(a)') asec//'fxi'
        WRITE(82,'(1pg24.16)') vpp0
        WRITE(82,'(100(1pg24.16))') dvpp(1:nparmx)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,13)
        WRITE(82,'(a)') asec//'h1xi'
        WRITE(82,'(1pg24.16)') vpn0
        WRITE(82,'(100(1pg24.16))') dvpn(1:nparmx)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,14)
        WRITE(82,'(a)') asec//'hgradxi'
        WRITE(82,'(1pg24.16)') rho0pr0
        WRITE(82,'(100(1pg24.16))') drho0pr(1:nparmx)
        IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,15)
        WRITE(82,'(a)') asec//'delfxi'
        WRITE(82,'(1pg24.16)') gammx0
        WRITE(82,'(100(1pg24.16))') dgammx(1:nparmx)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,19)
        WRITE(82,'(a)') asec//'delh1xi'
        WRITE(82,'(1pg24.16)') work0(20)
        WRITE(82,'(100(1pg24.16))') dwork(1:nparmx,20)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,20)
        WRITE(82,'(a)') asec//'delhgradxi'
        WRITE(82,'(1pg24.16)') work0(21)
        WRITE(82,'(100(1pg24.16))') dwork(1:nparmx,21)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,21)
        WRITE(82,'(a)') asec//'asm'
        WRITE(82,'(1pg24.16)') work0(22)
        WRITE(82,'(100(1pg24.16))') dwork(1:nparmx,22)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,22)
      ELSE
!        WRITE(82,'(a)') asec//'b4'
        WRITE(82,'(a)') asec//'B4_0'
        WRITE(82,'(1pg24.16)') b40
        WRITE(82,'(100(1pg24.16))') db4(1:nparmx)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,11)
        WRITE(82,'(a)') asec//'b4p'
        WRITE(82,'(1pg24.16)') b4p0
        WRITE(82,'(100(1pg24.16))') db4p(1:nparmx)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,12)
        WRITE(82,'(a)') asec//'b2'
        WRITE(82,'(1pg24.16)') b20
        WRITE(82,'(100(1pg24.16))') db2(1:nparmx)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,16)
        WRITE(82,'(a)') asec//'b2p'
        WRITE(82,'(1pg24.16)') b2p0
        WRITE(82,'(100(1pg24.16))') db2p(1:nparmx)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,17)
        WRITE(82,'(a)') asec//'pairfp'
        WRITE(82,'(1pg24.16)') vpp0
        WRITE(82,'(100(1pg24.16))') dvpp(1:nparmx)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,13)
        WRITE(82,'(a)') asec//'pairfn'
        WRITE(82,'(1pg24.16)') vpn0
        WRITE(82,'(100(1pg24.16))') dvpn(1:nparmx)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,14)
        WRITE(82,'(a)') asec//'rho0pr'
        WRITE(82,'(1pg24.16)') rho0pr0
        WRITE(82,'(100(1pg24.16))') drho0pr(1:nparmx)
        IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,15)
        WRITE(82,'(a)') asec//'tte'
        WRITE(82,'(1pg24.16)') tte0
        WRITE(82,'(100(1pg24.16))') dtte(1:nparmx)
        IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,18)
        WRITE(82,'(a)') asec//'tto'
        WRITE(82,'(1pg24.16)') tto0
        WRITE(82,'(100(1pg24.16))') dtto(1:nparmx)
        IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,19)
      END IF
      WRITE(82,'(a)') asec//'d_rho E/A-neut'
      WRITE(82,'(1pg24.16)') deneut00
      WRITE(82,'(100(1pg24.16))') ddeneut(1:nparmx)
      IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,9)
      WRITE(82,'(a)') asec//'E/A-neut'
      WRITE(82,'(1pg24.16)') eneut00
      WRITE(82,'(100(1pg24.16))') deneut0(1:nparmx)
      IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,8)
      WRITE(82,'(a)') asec//'E/A-neut12'
      WRITE(82,'(1pg24.16)') eneut120
      WRITE(82,'(100(1pg24.16))') deneut12(1:nparmx)
      IF(tsecord) 
     &   write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,10)
      IF(nntmax>0) THEN
        open(89,file='errorneut.'//id2)
        rhostep = rntmax/nntmax
        DO i=1,nntmax
         WRITE(89,'(a,1pg12.4)') ' E/A-neut',rhostep*i
         WRITE(89,'(1pg24.16)') eneutsav(i)
         WRITE(89,'(100(1pg24.16))') deneutsav(1:nparmx,i)
        END DO
        CLOSE(89)
      END IF
      IF(nsymax>0) THEN
        open(89,file='errorsym.'//id2)
        rhostep = rsymax/nsymax
        DO i=1,nsymax
         WRITE(89,'(a,f5.3)') ' E/A-sym',rhostep*i
         WRITE(89,'(1pg24.16)') esymsav(i)
         WRITE(89,'(100(1pg24.16))') desymsav(1:nparmx,i)
        END DO
        CLOSE(89)
      END IF
C
C     ACCUMULATE EXTRAPOLATION ERRORS BY FOLDING WITH 'ZARRAY'.
C
      ERRRHO = 0.0
      ERRCMP = 0.0
      ERRMAS = 0.0
      ERRENR = 0.0
      ERRSYM = 0.0
      errdsy = 0.0
      errska = 0.0
      errent = 0.0
      errden = 0.0
      errent12 = 0.0
C
      IF(NPARMX.GT.0) THEN
       DO I=1,NPARMX
         ERRRHO = DRHONM(I)*DRHONM(I) + ERRRHO
         ERRCMP = DCOMPR(I)*DCOMPR(I) + ERRCMP
         ERRMAS = DEFMAS(I)*DEFMAS(I) + ERRMAS
         ERRENR = DEOVRA(I)*DEOVRA(I) + ERRENR
         ERRSYM = DSYMEN(I)*DSYMEN(I) + ERRSYM
         ERRDSY = DDSYME(I)*DDSYME(I) + ERRDSY
         ERRSka = DSkapp(I)*DSkapp(I) + ERRSka
         errent = DENEUT0(I)*DENEUT0(I) + errent
         errden = DDENEUT(I)*DDENEUT(I) + errden
         errent12 = DENEUT12(I)*DENEUT12(I) + errent12
       enddo
       ERRRHO = SQRT(ERRRHO)/parmnatfac
       ERRCMP = SQRT(ERRCMP)/parmnatfac
       ERRMAS = SQRT(ERRMAS)/parmnatfac
       ERRENR = SQRT(ERRENR)/parmnatfac
       ERRSYM = SQRT(ERRSYM)/parmnatfac
       ERRDSY = SQRT(ERRDSY)/parmnatfac
       ERRSka = SQRT(ERRSka)/parmnatfac
       errent = SQRT(errent)/parmnatfac
       errden = SQRT(errden)/parmnatfac
       errent12 = SQRT(errent12)/parmnatfac
      ENDIF
C
C     FINAL PRINT
C
      write(02,'(a/a,a/1x,f9.4,f9.3,f9.1,f9.4,2f9.2,f9.4)')
     &     ' the nuclear matter parameters (natural evaluation):',
     &     '      rho      e/a     compr.  effmas     a4    d/drho_a4',
     &     '    skappa',
     &     rhonm0,eovera,compr,effmas,symenr,dsymme,skappa
      if(nparmx.gt.0) write(02,'(1x,f9.4,f9.2,f9.1,f9.4,2f9.2,f9.4)')
     &     errrho,errenr,errcmp,errmas,errsym,errdsy,errska
      write(22,'(a,f9.4,f9.2,f9.1,f9.4,2f9.2,f9.4)')
     &     commen(1:20),
     &     rhonm0,eovera,compr,effmas,symenr,dsymme,skappa
      write(02,'(a,2(1pg12.4)/a,2(1pg12.4)/a,2(1pg12.4))') 
     &     ' Neutron E/A  :',eneut00,errent,
     &     '         slope:',deneut00,errden,
     &     ' Neutron EA12 :',eneut120,erremt12
c
c
c     estimate curvature contributions
c
      IF(tsecord) THEN
        WRITE(2,'(a)') 'E/A-sym'
        CALL estim_curv(deovra,work2(:,:,1),nparmx,tstsecord)
        WRITE(2,'(a)') 'rho0-sym'
        CALL estim_curv(drhonm,work2(:,:,2),nparmx,tstsecord)
        WRITE(2,'(a)') 'incompr. K'
        CALL estim_curv(dcompr,work2(:,:,3),nparmx,tstsecord)
        WRITE(2,'(a)') 'm*/m'
        CALL estim_curv(defmas,work2(:,:,4),nparmx,tstsecord)
        WRITE(2,'(a)') 'a_sym'
        CALL estim_curv(dsymen,work2(:,:,5),nparmx,tstsecord)
        WRITE(2,'(a)') 'd_rho a_sym'
        CALL estim_curv(ddsyme,work2(:,:,7),nparmx,tstsecord)
        WRITE(2,'(a)') 'kappa'
        CALL estim_curv(dskapp,work2(:,:,6),nparmx,tstsecord)
        WRITE(2,'(a)') 'E/A-neut'
        CALL estim_curv(deneut0,work2(:,:,8),nparmx,tstsecord)
        WRITE(2,'(a)') 'd_rho E/A-neut' 
        CALL estim_curv(ddeneut,work2(:,:,9),nparmx,tstsecord)
        WRITE(2,'(a)') 'E/A-neut12' 
        CALL estim_curv(ddeneut,work2(:,:,10),nparmx,tstsecord)
      END IF
      CALL flush(2)
c
c     print error file for nclear and neutron matter
c
      IF(nntmax>0) DEALLOCATE(eneutsav,deneutsav)
      IF(nsymax>0) DEALLOCATE(esymsav,desymsav)
c
      return
      end
c-----extobs2_all--------------------------------------------------------
c
      subroutine extobs2_all(par,workdv,nparmx,iextobs)
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
      include 'constn.inc'
      dimension par(kparmx),workdv(kparmx,1)
c
c-------------------------------------------------------------------
c
      nextobs = 0            
      nextobs(iextobs) = 1
!      WRITE(2,'(a,100i3)') ' extobs: iextobs=',iextobs,nextobs

      IF(iextobs <= 20) THEN
        call extobs2(par,workdv,nparmx)
      ELSE
        call extobsmf2(par,workdv,nparmx)
      END IF
      CALL flush(2)

      RETURN
      END

c-----extobs--------------------------------------------------------
c
      subroutine extobs2(par,workdv,nparmx)
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
c      idipzr     = dipole polarizability in 90Zr
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
      REAL(8) ::   par(kparmx),pars(kparmx),workdv(kparmx,1)
      real*8  work2(kparmx,kparmx) !,vect(kparmx,kparmx),eig(kparmx)
c
      logical ttest
      parameter (ttest=.true.)
      logical tsymderiv
      parameter (tsymderiv=.true.)
      logical tsecord,tstsecord
      parameter (tsecord=.false.,tstsecord=.false.)
      character*1 asec
      character(3) :: idi,idj,idparm
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
     &   .ne.1) stop "EXTOBS2: only one option can be active"
ctest      ifplot=3
ctest      write(2,'(a)') 'EXTOBS entered'
ctest      write(6,'(a)') 'EXTOBS entered'
c
c     loop over variation of parameters
c
      IF(tsymderiv) THEN
        nparlow = -nparmx
      ELSE
        nparlow = 0
      END IF
      pars = par
      do i=nparlow,nparmx
      IF(tsecord) THEN
        nparm2=abs(i)
        nparlw2=-abs(i)
      ELSE
        nparm2=0
        nparlw2=0
      END IF
      do j=nparlw2,nparm2
        if(j.ne.0 .and. i.eq.0) CYCLE
        if(j.ne.0 .and. j.eq.-i) CYCLE
        if(j.eq.0 .and. i.eq.0) par(1:nparmx) = pars(1:nparmx)
        if(j.eq.0 .and. i.ne.0) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        if(j.ne.0 .and. j.ne.i) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))/2D0
     &             + SIGN(1,j)*parmnatfac*parmnat(1:nparmx,abs(j))/2D0
        if(j.ne.0 .and. j.eq.i) par(1:nparmx) = pars(1:nparmx)
     &             - SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        call attach(par,jp)
        call writ05()
        IF(igmrpb.EQ.1 .AND. ABS(parmnatfac-0.1D0)<1D-7) THEN
          CALL prep_ident(parmnatfac,i,j,idparm,idi,idj)
          call system('cp for005 for005.'//idparm//'.'//idi//'.'//idj)
          CALL  writsky(idparm,idi,idj)
        END IF
        call system('./haforpa.exe < for005 > for006')
        resenact = resenergy()
c
c       collect array of deviations
c
        CALL  collectw(i,j,resenact,resensav,workdv,work2,tsymderiv)
c
        par(1:nparmx) = pars(1:nparmx)
      enddo
      enddo
!      IF(igmrpb.EQ.1) STOP ' input files written'
c
c     accumulate extrapolation errors by folding with 'zarray'.
c
      erresen   = 0.0
c  
      if(nparmx.gt.0) then
       do i=1,nparmx
         erresen   = workdv(i,1)*workdv(i,1)
     &              + erresen
       enddo
       erresen   = sqrt(erresen)/parmnatfac
      endif

c
      IF(tsecord) THEN
        asec='#'
      ELSE
        asec=' '
      END IF
      if(igmrpb.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' GMR in 208PB:  goal,energy,uncertainty=',
     &    13.7,resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' GMR in 208PB:  goal,energy,uncertainty=',
     &    13.7,resensav,erresen
        write(82,'(a)') asec//'GMR 208Pb'
      endif
      if(igdrpb.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' GDR in 208PB:  goal,energy,uncertainty=',
     &    13.6,resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' GDR in 208PB:  goal,energy,uncertainty=',
     &    13.6,resensav,erresen
        write(82,'(a)') asec//'GDR 208Pb'
      endif
      if(igqrpb.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' GQR in 208PB:  goal,energy,uncertainty=',
     &    10.9,resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' GQR in 208PB:  goal,energy,uncertainty=',
     &    10.9,resensav,erresen
        write(82,'(a)') asec//'GQR 208Pb'
      endif
      if(igdro.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' GDR in 16O:  goal,energy,uncertainty=',
     &    23.0,resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' GDR in 16O:  goal,energy,uncertainty=',
     &    23.0,resensav,erresen
        write(82,'(a)') asec//'GDR 16O'
      endif
      if(idipol.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' polariz. in 208Pb: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' polariz. in 208Pb: value,uncertainty=',
     &       resensav,erresen
        write(82,'(a)') asec//'polariz. 208Pb'
      elseif(idipca.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' polariz. in 48Ca: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' polariz. in 48Ca: value,uncertainty=',
     &       resensav,erresen
        write(82,'(a)') asec//'polariz. 48Ca'
      elseif(idipsn.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' polariz. in 132Sn: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' polariz. in 132Sn: value,uncertainty=',
     &       resensav,erresen
        write(82,'(a)') asec//'polariz. 132Sn'
      elseif(idipni.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' polariz. in 68Ni: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' polariz. in 68Ni: value,uncertainty=',
     &       resensav,erresen
        write(82,'(a)') asec//'polariz. 68Ni'
      elseif(idipzr.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' polariz. in 90Zr: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' polariz. in 90Zr: value,uncertainty=',
     &       resensav,erresen
        write(82,'(a)') asec//'polariz. 90Zr'
      elseif(ipygmy.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' pygmy res. 208Pb: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' pygmy res. 208Pb: value,uncertainty=',
     &       resensav,erresen
        write(82,'(a)') asec//'pygmy res. 208Pb'
      elseif(ipygsn.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' pygmy res. 132Sn: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' pygmy res. 132Sn: value,uncertainty=',
     &       resensav,erresen
        write(82,'(a)') asec//'pygmy res. 132Sn'
      elseif(imonop.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' monop. polariz. in 208Pb: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' monop. polariz. in 208Pb: value,uncertainty=',
     &       resensav,erresen
        write(82,'(a)') asec//'monop.pol. 208Pb'
      elseif(imonsn.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' monop. polariz. in 132Sn: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' monop. polariz. in 132Sn: value,uncertainty=',
     &       resensav,erresen
        write(82,'(a)') asec//'monop.pol. 132Sn'
      elseif(imonni.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' monop. polariz. in 68Ni: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' monop. polariz. in 68Ni: value,uncertainty=',
     &       resensav,erresen
        write(82,'(a)') asec//'monop.pol. 68Ni'
      elseif(imonca.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' monop. polariz. in 48Ca: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' monop. polariz. in 48Ca: value,uncertainty=',
     &       resensav,erresen
        write(82,'(a)') asec//'monop.pol. 48Ca'
      elseif(iquadp.eq.1) then
        write(02,'(a,3f8.3)')
     &    ' quadr. polariz. in 208Pb: value,uncertainty=',
     &    resensav,erresen
        if(ttest) write(06,'(a,3f8.3)')
     &    ' quadr. polariz. in 208Pb: value,uncertainty=',
     &       resensav,erresen
        write(82,'(a)') asec//'quadr.pol. 208Pb'
      endif
      CALL flush(2)
c
c     print collection
c
      write(82,'(1pg24.16)') 
     &    resensav
      write(82,'(100(1pg24.16))') 
     &    (workdv(i,1),i=1,nparmx)
      IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx)
      CALL flush(82)
c
c     estimate curvature contributions
c
      IF(tsecord) CALL estim_curv(workdv,work2,nparmx,tstsecord)
c
      return
      end

c-----extobsmf2--------------------------------------------------------
c
      subroutine extobsmf2(par,workdv,nparmx)
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
c      ifrmne48    = weak-charge formfactor in 48Ca
c      ibindexo    = binding energy of very neutron rich 148Sn
c      irdfcoul    = charge diffraction radius in 208Pb
c      isigcoul    = charge surface thickness in 208Pb
c      irmscoul    = charge  rms radius in 208Pb
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
      REAL(8) ::   par(kparmx),pars(kparmx),workdv(kparmx,1)
      real*8  work2(kparmx,kparmx) !,vect(kparmx,kparmx),eig(kparmx)
      parameter (qtestneut=0.475D0)   ! Q [fm^-1] for testing neutron formfactor
      parameter (qtestne48=0.778D0)   ! Q [fm^-1] for testing neutron formfactor
      logical tsymderiv
      parameter (tsymderiv=.true.)
      logical tsecord,tstsecord
      parameter (tsecord=.false.,tstsecord=.false.)
      character*1 asec
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
     &       + ifrmne48-ifrmne48*ifrmne48
     &       + ibindexo-ibindexo*ibindexo
     &       + irdfcoul-irdfcoul*irdfcoul
     &       + isigcoul-isigcoul*isigcoul
     &       + irmscoul-irmscoul*irmscoul
     &       + isurfe-isurfe*isurfe
     &       + isurfas-isurfas*isurfas
     &       + ifrmapv-ifrmapv*ifrmapv
      if(itest.ne.0) stop 'EXTOBS: options not 0 or 1'
      if(iskinpb+ilevelsn+ishe+iqalpha+is2n+is2p+iisotop
     &   +iskinsn74+iskinsn82+iskinsn90+iskinca+iskinni
     &   +iskindf+irmsneut+irdfneut+ifrmneut+ifrmne48+ibindexo
     &   +irdfcoul+isigcoul+irmscoul+isurfe+isurfas+ifrmapv
     &   .ne.1) 
     &   stop "EXTOBSMF: only one option can be active"
c
c     loop over variation of parameters
c
      IF(tsymderiv) THEN
        nparlow = -nparmx
      ELSE
        nparlow = 0
      END IF
      pars = par
      do i=nparlow,nparmx
      IF(tsecord) THEN
        nparm2=abs(i)
        nparlw2=-abs(i)
      ELSE
        nparm2=0
        nparlw2=0
      END IF
      do j=nparlw2,nparm2
        if(j.ne.0 .and. i.eq.0) CYCLE
        if(j.ne.0 .and. j.eq.-i) CYCLE
        if(j.eq.0 .and. i.eq.0) par(1:nparmx) = pars(1:nparmx)
        if(j.eq.0 .and. i.ne.0) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        if(j.ne.0 .and. j.ne.i) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))/2D0
     &             + SIGN(1,j)*parmnatfac*parmnat(1:nparmx,abs(j))/2D0
        if(j.ne.0 .and. j.eq.i) par(1:nparmx) = pars(1:nparmx)
     &             - SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        call attach(par,jp)
c
        if(iskinpb.eq.1 .or. iskindf.eq.1 .or. irmsneut.eq.1
     &     .or. irdfneut.eq.1 .or. ifrmneut.eq.1 
     &     .or. irdfcoul.eq.1 .or. isigcoul.eq.1
     &     .or. irmscoul.eq.1 .or. ifrmapv.eq.1) then
          nprot  = 82
          nneut  = 126
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          if(iskinpb.eq.1) observ = rmsn-rmsp
!     &          rmsn-sqrt(0.6*rdmsn*rdmsn+3.0*surfn*surfn) !tmp rmsn-rmsp
          if(iskindf.eq.1) observ = rdmsn-rdmsp   !tmp surfn
          if(irmsneut.eq.1) observ = rmsn 
          if(irdfneut.eq.1) observ = rdmsn
          if(ifrmneut.eq.1) then
            observ = rmsw - 5.503 !fbint(qtestneut,fow) 
          endif
          if(ifrmapv.eq.1) observ = apv_weak
          if(irdfcoul.eq.1) observ = rdmsc
          if(isigcoul.eq.1) observ = surfc
          if(irmscoul.eq.1) observ = rmsc
        elseif(iskinca.eq.1 .or. ifrmne48.eq.1) then
          nprot  =  20
          nneut  =  28
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          if(iskinca.eq.1)  observ = rmsn -rmsp
          if(ifrmne48.eq.1) observ = apv_weak !fbint(qtestne48,fow)
!          if(ifrmne48.eq.1) observ = rms4c
!          WRITE(2,'(a,i3,1pg18.10)') 'test:',i,observ
!          CALL flush(2)
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
          nprot  = 20  !50
          nneut  = 20  !90
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = rmsn-rmsp
        elseif(ibindexo.eq.1) then
          nprot  = 50
          nneut  = 98
          IF(tsimple) nneut=90
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          IF(tsimple) THEN
            CALL esimple()
          ELSE   
            call hafosu
          END IF
          observ = etotal
        elseif(ilevelsn.eq.1) then  ! misuse for weak skin ! misuse for gap 44Ca
          nprot  =  50
          nneut  =  82
          nnumb  = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          call hafosu
          observ = spenrg(nmaxp+16)-spenrg(nmaxp+14)
!          observ = deltn
!          observ = rmsw-rmsc
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
        elseif(isurfe.eq.1) then
          if(abs(b4p).gt.1.0E-7) T4 = (t4+b4p)*0.666666667
          call ASYMT(T0,T1,T2,T3,t4,X0,X1,X2,X3,power,
     &               (amfc0(1)+amfc0(2))/2D0,
     +               YRHO0,eovera,xKOMPr,a4,YEFFM,YXKF,YELLO,
     +               surfe,surfas,YH,YP,YG,YQT1,YQT2)
          observ = surfe
        elseif(isurfas.eq.1) then
          if(abs(b4p).gt.1.0E-7) T4 = (t4+b4p)*0.666666667
          call ASYMT(T0,T1,T2,T3,t4,X0,X1,X2,X3,power,
     &               (amfc0(1)+amfc0(2))/2D0,
     +               YRHO0,eovera,xKOMPr,a4,YEFFM,YXKF,YELLO,
     +               surfe,surfas,YH,YP,YG,YQT1,YQT2)
          observ = surfas
        endif
cold      ngrid  = 0
c
c       collect array of deviations
c
        CALL  collectw(i,j,observ,observsav,workdv,work2,tsymderiv)
c
        par(1:nparmx) = pars(1:nparmx)
      enddo
      enddo
c
c     accumulate extrapolation errors by folding with 'zarray'.
c
      errobs   = 0.0
c
      if(nparmx.gt.0) then
        do i=1,nparmx
          errobs  = workdv(i,1)*workdv(i,1) + errobs
        enddo
        errobs  = sqrt(errobs)/parmnatfac
      endif
c
c     final print
c
      IF(tsecord) THEN
        asec='#'
      ELSE
        asec=' '
      END IF
        if(iskinpb.eq.1) then
          write(02,'(a,2f8.3)')
     &     'rms neutron skin 208Pb: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'skin 208Pb'
        elseif(iskinca.eq.1) then
          write(02,'(a,2f8.3)')
     &     'rms neutron skin 48Ca: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'skin 48Ca'
        elseif(iskinni.eq.1) then
          write(02,'(a,2f8.3)')
     &     'rms neutron skin 68Ni: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'skin 68Ni'
        elseif(iskinsn74.eq.1) then
          write(02,'(a,2f8.3)')
     &     'rms neutron skin 124Sn: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'skin 124Sn'
        elseif(iskinsn82.eq.1) then
          write(02,'(a,2f8.3)')
     &     'rms neutron skin 132Sn: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'skin 132Sn'
        elseif(iskinsn90.eq.1) then
          write(02,'(a,2f8.3)')
     &     'rms neutron skin 140Sn: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'skin 140Sn'
        elseif(ilevelsn.eq.1) then
          write(02,'(a,2f8.3)')
     &     'neutron level diff. 132Sn: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'level 132Sn'
        elseif(ishe.eq.1) then
          write(02,'(a,f10.3,f8.3)')
     &     'binding energy 120/182: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'EB(120/182)'
        elseif(iqalpha.eq.1) then
          write(02,'(a,2f8.3)')
     &     'Q_alpha 120/182: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'Qalpha(120/182)'
        elseif(is2n.eq.1) then
          write(02,'(a,2f8.3)')
     &     'S_2n 120/182: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'S2n(120/182)'
        elseif(is2p.eq.1) then
          write(02,'(a,2f8.3)')
     &     'S_2p 120/182: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'S2p(120/182)'
        elseif(iisotop.eq.1) then
          write(02,'(a,2f8.3)')
     &     'iso-shift 208Pb: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'r(214)-r(208)'
        elseif(iskindf.eq.1) then
          write(02,'(a,2f8.3)')
     &     'diff.skin 208Pb: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'dskin 208Pb'
        elseif(irmsneut.eq.1) then
          write(02,'(a,2f8.3)')
     &     'neutron rms: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'rmsn 208Pb'
        elseif(irdfneut.eq.1) then
          write(02,'(a,2f8.3)')
     &     'neutron diff.rad.: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'rdmsn 208Pb'
        elseif(ifrmneut.eq.1) then
          write(02,'(a,2f8.3)')
     &     'weak-charge F_W(0.4): value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'formw 208Pb'
        elseif(ifrmapv.eq.1) then
          write(02,'(a,2f8.3)')
     &     'weak-charge A_PV(0.4): value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'A_PV 208Pb'
        elseif(ifrmne48.eq.1) then
          write(02,'(a,2f8.3)')
     &     'weak-charge 48Ca F_W(0.778): value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'formw 48Ca'
        elseif(ibindexo.eq.1) then
          write(02,'(a,f12.3,f8.3)')
     &     'binding energy 148Sn: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'E_B 148Sn '
        elseif(irdfcoul.eq.1) then
          write(02,'(a,f12.3,f8.3)')
     &     ' charge R_diffr 208Pb: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'ch.R_diffr 208Pb '
        elseif(isigcoul.eq.1) then
          write(02,'(a,f12.3,f8.3)')
     &     ' charge surfac.th. 208Pb: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'ch.surf.th 208Pb '
        elseif(irmscoul.eq.1) then
          write(02,'(a,f12.3,f8.3)')
     &     ' charge rms radius 208Pb: value,uncertainty=',
     &     observsav,errobs
          write(82,'(a)') asec//'ch.rms rad. 208Pb '
        elseif(isurfe.eq.1) then
          write(02,'(a,f12.3,f8.3)')
     &     ' semi-classical surface energy=',
     &     observsav,errobs
          write(82,'(a)') asec//'asurf SC         '
        elseif(isurfas.eq.1) then
          write(02,'(a,f12.3,f8.3)')
     &     ' semi-classical surface asymmetry=',
     &     observsav,errobs
          write(82,'(a)') asec//'asurf-symm SC    '
        endif
        write(82,'(1pg24.16)')  observsav
        write(82,'(100(1pg24.16))') workdv(1:nparmx,1)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx)
        CALL flush(2)
        CALL flush(82)
c
c     estimate curvature contributions
c
      IF(tsecord) CALL estim_curv(workdv,work2,nparmx,tstsecord)
c
      return
      end

c-----extobsFq--------------------------------------------------------
c
      subroutine extobsFq2(par,workdv,nparmx)
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
      REAL(8) ::   par(kparmx),pars(kparmx),workdv(kparmx,1)
      parameter (qwmax=2D0,nqwmax=40)
      double precision formfw(0:kparmx,nqwmax),surfsav(0:kparmx)
      double precision rmssav(0:kparmx),rdmssav(0:kparmx)
      double precision halosav(0:kparmx)
      character(3) :: a(1:4)
      integer,save :: icount=1
c
c-----------------------------------------------------------------------
c
!      WRITE(2,'(a)') ' extobsFw: '
      a(1)='010'
      a(2)='025'
      a(3)='050'
      a(4)='100'
c
c     loop over variation of parameters
c
      if(icount>4) stop 'EXTOBSFQ2: icount exhausted'
      open(281,file='error-Fw.'//a(icount))
      icount = 1+ icount
      qwstep = qwmax/nqwmax
      pars = par
      do i=0,nparmx
        if(i.gt.0) par(1:nparmx) = pars(1:nparmx)
     &               + parmnatfac*parmnat(1:nparmx,i)
        call attach(par,jp)
c
        nprot  = 82
        nneut  = 126
        nnumb  = nprot+nneut
        anumb  = nnumb
        nmaxp  = 0
        nmaxn  = 0
        call hafosu
        do iqw=1,nqwmax
          formfw(i,iqw) = fbint(qwstep*iqw,fow)
        enddo
        if(i.EQ.0) then
          rmssav(i) = rmsw
          rdmssav(i) = rdmsw
          surfsav(i) = surfw
          halosav(i) = rmsw-sqrt(0.6*rdmsw*rdmsw+3.0*surfw*surfw)
        else
          rmssav(i) = rmsw-rmssav(0)
          rdmssav(i) = rdmsw-rdmssav(0)
          surfsav(i) = surfw-surfsav(0)
          halosav(i) = rmsw-sqrt(0.6*rdmsw*rdmsw+3.0*surfw*surfw)
     &                 -halosav(0)
        endif
        if(i.gt.0) par(1:nparmx) = pars(1:nparmx)
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
            workdv(i,1) = (formfw(i,iqw)-formfw(0,iqw))
        enddo
c
c       accumulate extrapolation errors by folding with 'zarray'.
c
        observsav = formfw(0,iqw)
        errobs   = 0.0
        if(nparmx.gt.0) then
          do i=1,nparmx
            errobs  = workdv(i,1)*workdv(i,1) + errobs
          enddo
          errobs  = sqrt(errobs)/parmnatfac
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
c
c     evaluate errors for the three basic form-parameters
c
      errobs   = 0.0
      if(nparmx.gt.0) then
        do i=1,nparmx
          errobs  = rmssav(i)**2 + errobs
        enddo
        errobs  = sqrt(errobs)/parmnatfac
      endif
      write(02,'(a,2f8.3)')
     &     'rms_W value,uncertainty=',rmssav(0),errobs
      write(281,'(a)') ' rms_W'
      write(281,'(1pg24.16)') rmssav(0)
      write(281,'(100(1pg24.16))') (rmssav(i),i=1,nparmx)
c
      errobs   = 0.0
      if(nparmx.gt.0) then
        do i=1,nparmx
          errobs  = rdmssav(i)**2 + errobs
        enddo
        errobs  = sqrt(errobs)/parmnatfac
      endif
      write(02,'(a,2f8.3)')
     &     'rdms_W value,uncertainty=',rdmssav(0),errobs
      write(281,'(a)') ' rdms_W'
      write(281,'(1pg24.16)') rdmssav(0)
      write(281,'(100(1pg24.16))') (rdmssav(i),i=1,nparmx)
c
      errobs   = 0.0
      if(nparmx.gt.0) then
        do i=1,nparmx
          errobs  = surfsav(i)**2 + errobs
        enddo
        errobs  = sqrt(errobs)/parmnatfac
      endif
      write(02,'(a,2f8.3)')
     &     'surf_W value,uncertainty=',surfsav(0),errobs
      write(281,'(a)') ' surf_W'
      write(281,'(1pg24.16)') surfsav(0)
      write(281,'(100(1pg24.16))') (surfsav(i),i=1,nparmx)
c
      errobs   = 0.0
      if(nparmx.gt.0) then
        do i=1,nparmx
          errobs  = halosav(i)**2 + errobs
        enddo
        errobs  = sqrt(errobs)/parmnatfac
      endif
      write(02,'(a,2f8.3)')
     &     'halo_W value,uncertainty=',halosav(0),errobs
      write(281,'(a)') ' halo_W'
      write(281,'(1pg24.16)') halosav(0)
      write(281,'(100(1pg24.16))') (halosav(i),i=1,nparmx)
c
c     all done
c
      close(281)
c
      return
      end
c-----extobsFq--------------------------------------------------------
c
      subroutine extobsLanPar(par,workdv,nparmx)
c
c     Extrapolate to Landau parameters including extrapolation error using
c     matrix 'zarray' which is transferred via common /fitpar/ .
c     Eight different Landau parameters are printed F0s,F0v,F1s,F1v
c     from the rho-dependent normalized class and correspondingly
c     f0s,f0v,f1s,f1c from the class with fixed normalization.
c     The label 's' in  the naming means isoscalar parameters and
c     the label 'v' means isovector parameters.
c
c     The densities at which the errors are computed are given
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
      REAL(8) ::   par(kparmx),pars(kparmx),workdv(kparmx,1)
      parameter (rhomax=0.2D0,nmaxrho=20)
      double precision flanpar(0:kparmx,nmaxrho)
      character(3) :: a(1:4)
      integer,save :: icount=1
      INTEGER :: icase, numcases=10
      CHARACTER(3) :: names(10)=(/'f0s','f0v','f1s','f1v',
     &                            'P0s','P0v','F1s','F1v','asy','inK'/)
!      CHARACTER(3) :: names(10)=(/'f0s','f0v','f1s','f1v',
!     &                            'F0s','F0v','F1s','F1v','asy','inK'/)
!      CHARACTER(3),DIMENSION(8) :: names='F0sF0vF1sF1vf0sf0vf1sf1v'
c
c-----------------------------------------------------------------------
c
!      WRITE(2,'(a)') ' extobsFw: '
      a(1)='010'
      a(2)='025'
      a(3)='050'
      a(4)='100'
c
c     loop over variation of parameters
c
!      if(icount>4) stop 'EXTOBSLANPAR: icount exhausted'
      if(icount>1) stop 'EXTOBSLANPAR: icount exhausted'
      rhostep = rhomax/nmaxrho
c
      casesloop: DO icase=1,numcases
      IF(icase==7 .OR. icase==8) CYCLE
      open(281,file='error-'//names(icase)(1:3)//'.'//a(icount))
      WRITE(2,*) 'icase=',icase,' opened: ',
     &     'error-'//names(icase)(1:3)//'.'//a(icount)
      pars = par
      do i=0,nparmx
        if(i.gt.0) par(1:nparmx) = pars(1:nparmx)
     &               + parmnatfac*parmnat(1:nparmx,i)
        call attach(par,jp)
        call initfc()
        DO irho=1,nmaxrho
          rhoact = irho*rhostep
          call lanpar(rhoact,flan0,flan0p,glan0,glan0p,f0,f0p,g0,g0p,
     &                f1,f1p,c0fac)
          SELECT CASE(icase)
            CASE(1)
              flanpar(i,irho) = f0
            CASE(2)
              flanpar(i,irho) = f0p
            CASE(3)
              flanpar(i,irho) = f1
            CASE(4)
              flanpar(i,irho) = f1p
            CASE(5)
              flanpar(i,irho) = (f0-f1/3D0)/(1+f1/3D0)  !f0*c0fac/150D0
            CASE(6)
              flanpar(i,irho) = (f0p-f1/3D0)/(1+f1/3D0) !f0p*c0fac/150D0
            CASE(7)
              flanpar(i,irho) = f1*c0fac/150D0
            CASE(8)
              flanpar(i,irho) = f1p*c0fac/150D0
            CASE(9)
              call enernm(rhoact,half-0.001,zero,energ2,efmasp,efmasn)
              call enernm(rhoact,half,zero,energ0,efmasp,efmasn)
              call enernm(rhoact,half+0.001,zero,energ1,efmasp,efmasn)
              asymact = (energ1-2.0*energ0+energ2)/(0.001*0.001)/8.0
              flanpar(i,irho) = asymact
            CASE(10)
              call enernm(rhoact+0.01,half,zero,eplus,efmasp,efmasn)
              call enernm(rhoact-0.01,half,zero,eminus,efmasp,efmasn)
              call enernm(rhoact     ,half,zero,energ0,efmasp,efmasn)
              compract = 9.0*rhoact*rhoact*(eminus-energ0-energ0+eplus)
     *                          /(0.01*0.01)
              flanpar(i,irho) = compract
          END SELECT
        END DO
        if(i.gt.0) par(1:nparmx) = pars(1:nparmx)
      enddo
cold      ngrid  = 0
c
c     loop over all momenta and compute error relations
c
      DO irho=1,nmaxrho
c
c       collect array of deviations
c
        do i=1,nparmx
            workdv(i,1) = (flanpar(i,irho)-flanpar(0,irho))
        enddo
c
c       accumulate extrapolation errors by folding with 'zarray'.
c
        observsav = flanpar(0,irho)
        errobs   = 0D0
        if(nparmx.gt.0) then
          do i=1,nparmx
            errobs  = workdv(i,1)*workdv(i,1) + errobs
          enddo
          errobs  = sqrt(errobs)/parmnatfac
        endif
c
c     final print
c
        write(02,'(a,f5.2,a,2f8.3)')
     &     names(icase)//'(rho=',irho*rhostep,') value,uncertainty=',
     &     observsav,errobs
        write(281,'(a,f4.2,a)') names(icase)//'(',irho*rhostep,')'
        write(281,'(1pg24.16)')  observsav
        write(281,'(100(1pg24.16))') 
     &      (workdv(i,1),i=1,nparmx)
c
      END DO
      close(281)
      END DO casesloop
      icount = 1+ icount
c
      return
      end

c-----extbarr--------------------------------------------------------
c
      subroutine extbarr2(par,nparmx)
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
      dimension   par(kparmx),pars(kparmx)
c
      character*3 ident,id2
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
      pars = par
      do i=0,nparmx
        if(i.gt.0) par(1:nparmx) = pars(1:nparmx) 
     &             + parmnatfac*parmnat(1:nparmx,i)
        call attach(par,jp)
c
c       prepare identifier
c
        open(31,status='scratch')
        write(31,'(i3)') i+100
        rewind(31)
        read(31,'(a3)') ident
        rewind(31)
        write(31,'(i4)') nint(parmnatfac*100+1000)
        rewind(31)
        read(31,'(1x,a3)') id2
        close(31)
        ident(1:1) = 's' 
c
c       prepare input for NAMELIST 'force' in 'SkyAx'
c
        open(31,file='force2D.'//id2//'.'//ident)
        WRITE(31,'(a)') "&force name='none,'"
        WRITE(31,'(a,f5.2,3(a,i3)/3(a,L1),3(a,f6.2)/a,2(1pg18.10))')
     &   'chg=',acoul,
     &   ',ex=',ifex,
     &   ',zpe=',ifzpe,
     &   ',ipair=',ipair,
     &   ',tcoulfold=',ifrhoc, 
     &   ',tcutoff_uvspace=',tphpair,
     &   ',tregpair=',tregpair,
     &   ',pairstab=',-ecstab,
     &   ',pair_cutoff=',cutfac,
     &   ',hfbwidth=',hfbwidth,
     &   ',h2m=',amfc0
        IF(ioptsk==13) THEN
          WRITE(31,'(a,L1/6(a,1pg20.12)/6(a,1pg20.12)/
     &              6(a,1pg20.12)/5(a,1pg20.12)/a,i2)')
     &   ',tFayans=',.TRUE.,
     &   ',avp=',avp,
     &   ',h1vp=',h1vp,
     &   ',h2vp=',h2vp,
     &   ',avm=',avm,
     &   ',h1vm=',h1vm,
     &   ',h2vm=',h2vm,
     &   ',asp=',asp,
     &   ',hsp=',hsp,
     &   ',hgradp=',hgradp,
     &   ',asm=',asm,
     &   ',hsm=',hsm,
     &   ',hgradm=',hgradm, 
     &   ',fxi=',fxi,
     &   ',h1xi=',h1xi,
     &   ',hgradxi=',hgradxi,
     &   ',delfxi=',delfxi,
     &   ',delh1xi=',delh1xi,
     &   ',delhgradxi=',delhgradxi,
     &   ',c0tau=',c0tau,
     &   ',c1tau=',c1tau,
     &   ',c0nabj=',c0nabj,
     &   ',c1nabj=',c1nabj,
     &   ',hcoul=',hcoul, 
     &   ',isovectormodel=',isovectormodel
       ELSE
         WRITE(31,'(5(a,1pg20.12)/5(a,1pg20.12)/
     &              5(a,1pg20.12))')
     &   't0=',t0,
     &   ',t1=',t1,
     &   ',t2=',t2,
     &   ',t3=',t3,
     &   ',t4=',t4,
     &   ',x0=',x0,
     &   ',x1=',x1,
     &   ',x2=',x2,
     &   ',x3=',x3,
     &   ',b4p=',b4p,
     &   ',power=',power,
     &   ',riftm=',riftm, 
     &   ',v0prot=',v0prot,
     &   ',v0neut=',v0neut,
     &   ',rho0pr=',rho0pr
      END IF
      WRITE(31,'(a)') '&end'
      CLOSE(31)
!         isTsym                                         !NEW ODD                                                                         



c
c       write for005.*
c
        nprot = 120
        nneut = 182
c                                         for the ground state
        open(31,file='skyinpgs.'//id2//'.'//ident)
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
        open(31,file='skyinpbar.'//id2//'.'//ident)
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
        open(31,file='skyinphs.'//id2//'.'//ident)
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
        if(ioptsk.eq.12) then
          WRITE(31,'(a,1pg17.9)') ' Gsigma=',par(jp(1))  
          WRITE(31,'(a,1pg17.9)') ' Gomega=',par(jp(2))  
          WRITE(31,'(a,1pg17.9)') ' Grho=',par(jp(3))    
          WRITE(31,'(a,1pg17.9)') ' isoscalarMu=',par(jp(4))
          WRITE(31,'(a,1pg17.9)') ' Msigma=',par(jp(6))     
          WRITE(31,'(a,1pg17.9)') ' Momega=',par(jp(7))     
          WRITE(31,'(a,1pg17.9)') ' Mrho=',par(jp(8))       
          WRITE(31,'(a,1pg17.9)') ' isovectorMu=',par(jp(9))
          WRITE(31,'(a,1pg17.9)') ' dpol=',par(jp(10))      
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
        if(i.gt.0) par(1:nparmx) = pars(1:nparmx)
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
      pars = par
      do i=0,nparmx
        if(i.gt.0) par(1:nparmx) = pars(1:nparmx)
     &               + parmnatfac*parmnat(1:nparmx,i)
        call attach(par,jp)
        open(31,status='scratch')
        write(31,'(i3)') i+100
        rewind(31)
        read(31,'(a3)') ident
        rewind(31)
        write(31,'(i4)') nint(parmnatfac*100+1000)
        rewind(31)
        read(31,'(1x,a3)') id2
        close(31)
        ident(1:1) = 's' 
        open(31,file='for005.'//id2//'.'//ident)
c
c     the mean field part
c
!        WRITE(2,'(3a,i4)') '#protoc for005:',id2,ident,ioptsk
        write(31,'(a)')
     &   '1,,    IWUNIT',
     &   "' resonances ion 208 Pb'",
     &   '0,0.3,              NGRID,RSTEP',
     &   '-1,1E-10,1E-10,300,15000,3,-1, IWR,DABR/E,ITPR/MAX,IFPL',
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
          iu=31
          WRITE(iu,'(3(1pg20.12),a)') avp,h1vp,h2vp,' avp,h1vp,h2vp'
          WRITE(iu,'(5(1pg20.12))') avm,h1vm,h2vm
          WRITE(iu,'(5(1pg20.12))') sigmavol,rho0fy
          WRITE(iu,'(5(1pg20.12))') asp,hsp,hgradp 
          WRITE(iu,'(5(1pg20.12))') asm,hsm,hgradm 
          WRITE(iu,'(5(1pg20.12))') hcoul,sigmacoul
          WRITE(iu,'(5(1pg20.12))') fxi,h1xi,hgradxi,gammaxi,rho0xi
          WRITE(iu,'(4(1pg20.12),i2)') delfxi,ecstab,delh1xi,delhgradxi,
     &            isovectormodel
          WRITE(iu,'(5(1pg20.12))') c0tau,c1tau 
          WRITE(iu,'(5(1pg20.12))') c0nabj,c1nabj,c0tens,c1tens
      ELSE
        write(31,'(5(1pg19.11,a)/5(g19.11,a))')
     &  t0,',',t1,',',t2,',',t3,',',t4,',      t0,t1,t2,t3,t4',
     &  x0,',',x1,',',x2,',',x3,',',b4p,',       x0,x1,x2,x3,b4p'
        if(ioptsk.eq.5 .or. ioptsk.eq.6 .or. ioptsk.eq.11) then
          write(31,'(2(1pg16.8,a)/2(1pg16.8,a))')
     &      c1,',',c4,',                      c1 c4',
     &      c1p,',',c4p,',                    c1p,c4p'
        endif
        if(ioptsk.eq.8 .or. ioptsk.eq.11) then
          write(31,'(1(1pg16.8,a)/2(1pg16.8,a))')
     &      b3tilde,',                      b3tilde',
     &      b12,',',b12p,',                     b12,b12p'
        endif
        if(ioptsk.eq.9 .or. ioptsk.eq.11) then
          write(31,'(2(1pg16.8,a)/1pg16.8,a)')
     &      b7,',',b7p,',               b7,b7p',
     &      pow7,',               pow7'
        endif
        if(ioptsk.eq.10 .or. ioptsk.eq.11) then
          write(31,'(8(1pg19.11),a)')
     &      riftm,rho0pr,powerp,ecstab,alpinv,betinv,
     &      alpin2,betin2,
     &    ',      riftm,rho0pr,powerp,ecstab,alpinv'
        endif
        if(ioptsk.le.4 .or. ioptsk.eq.7) then
          write(31,'(5(f10.3,a)/4(f10.6,a))')
     &    t5,',',t6,',',t7,',',t8,',',t9,',    t5,t6,t7,t8,t9',
     &    x5,',',x6,',',x7,',',x8,',           x5,x6,x7,x8'
        endif
        write(31,'((f10.3,a),2(f10.6,a))')
     &    tyuk,',',xyuk,',',ayuk,',           tyuk,xyuk,ayuk'
      END IF
        write(31,'(f15.11,a,i2,2(1pg20.12),a,2(1pg19.11),a)')
     &     power,',1.0,',ipair,pairfp,pairfn,
     &  ',0,,',cutfac,hfbwidth,',, POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,'
c        write(31,'(f15.11,a,a)')
c     &     power,',0.0,0,0.2,0.2',
c     &  ',0,,,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,...'
        if(i.gt.0) par(1:nparmx) = pars(1:nparmx)
        close(31)
      enddo
c
      return
      end

      subroutine collectw(i,j,observ,observsav,workdv,work2,tsymderiv)
c
c       collect array of deviations
c
      include 'params.inc'
      include 'comfit.inc'
c
      logical tsymderiv
      integer i,j
      real*8  observ,observsav,workdv(kparmx,1),work2(kparmx,kparmx)
c
c
!        IF(j==0) WRITE(2,*) ' iparam,obs=',i,observ
        if(i.eq.0 .and. j.eq.0) then
          observsav = observ
        elseif(j.eq.0) then
          IF(i<0) THEN
            workdv(abs(i),1) = observ
          ELSE IF(tsymderiv .AND. i>0) THEN
            workdv(i,1) = (observ-workdv(i,1))/2D0
          ELSE
            workdv(i,1) = (observ-observsav)
          END IF
        elseif(j.ne.0 .and. abs(j).ne.abs(i)) then
          IF(j<0 .AND. i<0) THEN
            work2(abs(j),abs(i)) = observ
          ELSE IF(j>0 .AND. i<0) THEN
            work2(j,abs(i)) = work2(j,abs(i))-observ
          ELSE IF(j<0 .AND. i>0) THEN
            work2(abs(j),i) = work2(abs(j),i)-observ
          ELSE IF(i>0 .AND. j>0) THEN
            work2(j,i) = (work2(j,i)+observ)
            work2(i,j) = work2(j,i)
          END IF
        elseif(j.ne.0 .and. j.eq.i) then
          IF(i<0) THEN
            work2(abs(i),abs(i)) = observ
!            WRITE(2,*) ' i<0,obs=',i,observ
          ELSE IF(tsymderiv .AND. i>0) THEN
!            WRITE(2,*) ' i,work2=',i,work2(i,i)
            work2(i,i) = (work2(i,i)+observ-2D0*observsav)
!            WRITE(2,*) ' i,obs,work=',i,observ,work2(i,i),observsav
          ELSE
            work2(i,i) = (workdv(i,1)+observ-observsav)
          END IF
        endif
c
!        IF(j==0) WRITE(2,*)  ' iparam,slope=',i,workdv(i,1)
      return
      end

      subroutine estim_curv(workdv,work2,nparmx,tstsecord)
c
c     estimate curvature contributions
c
      include 'params.inc'
      include 'comfit.inc'
c
      logical tstsecord
      integer nparmx
      real*8 workdv(kparmx,1)
      real*8 work2(kparmx,kparmx),vect(kparmx,kparmx),eig(kparmx)
c
c
!        WRITE(2,'(a)') 'work2:'
        derivfac=1D0/parmnatfac**2
        do i=1,nparmx
        do j=1,nparmx
          IF(i.ne.j) THEN
!            work2(i,j) = (work2(i,j)-workdv(i,1)-workdv(j,1))*derivfac
            work2(i,j) = work2(i,j)*derivfac
          ELSE
!            work2(i,i) = (workdv(i,1)-work2(i,i))*derivfac
            work2(i,i) = work2(i,i)*derivfac
          ENDIF
        enddo
!        WRITE(2,'(200(1pg13.5))') (work2(i,j),j=1,nparmx)
        enddo
        workdv(1:nparmx,1) = workdv(1:nparmx,1)/parmnatfac
        IF(tstsecord) THEN
        WRITE(2,'(a)') 'diagonals:'
          do i=1,nparmx
            WRITE(2,'(i5,3(1pg13.5))') 
     &        i,workdv(i,1)*2D0,work2(i,i),work2(i,i)/(workdv(i,1)*2D0)
          enddo
        END IF
        vmax1 = SQRT(SUM(workdv(1:nparmx,1)**2))
        CALL diag_real_symm(work2,eig,vect,nparmx,kparmx,kparmx*kparmx)
        IF(tstsecord) WRITE(2,'(a,38(1pg13.5))') 'eig:',eig(1:nparmx)
        vmax2 = 0D0
        do i=1,nparmx
          vmax2 = max(vmax2,abs(eig(i)))        
        enddo
!        prederr=0D0
        trace1 = 0D0
        trace2 = 0D0
        do i=1,nparmx
!          prederr = prederr +  workdv(i,1)**2+workdv(i,1)*work2(i,i)
          trace1 = trace1 + work2(i,i)
          trace2 = trace2 + SUM(work2(i,1:nparmx)*work2(1:nparmx,i))
        enddo
        stderr = SUM(workdv(1:nparmx,1)**2)
        prederr = SQRT(stderr + trace1**2 + trace2)
        stderr = SQRT(stderr)
        WRITE(2,'(a,4(1pg13.5))')
     &         'curvature:',vmax2/(2*vmax1),vmax1,vmax2,prederr,stderr
c
      return
      end
c-----prep_ident-----------------------------------------------------
c
      subroutine prep_ident(parmnatfac,i,j,idparm,idi,idj)
c
c     prepare identifier
c
      INTEGER, INTENT(IN) :: i,j
      REAL(8), INTENT(IN) :: parmnatfac
      CHARACTER(3), INTENT(OUT) :: idparm,idi,idj
c
c--------------------------------------------------------------------
c
      open(31,status='scratch')
      write(31,'(i3)') i+100
      rewind(31)
      read(31,'(a3)') idi
      rewind(31)
      write(31,'(i3)') j+100
      rewind(31)
      read(31,'(a3)') idj
      rewind(31)
      write(31,'(i4)') nint(parmnatfac*100+1000)
      rewind(31)
      read(31,'(1x,a3)') idparm
      close(31)
      idi(1:1) = 's' 
      idj(1:1) = 't' 
c
      return
      end


c-----writsky-----------------------------------------------------
c
      subroutine writsky(idparm,idi,idj)
c
c     prepare identifier
c
      include 'params.inc'
      include 'constn.inc'
      include 'skypar.inc'
c
      CHARACTER(3), INTENT(in) :: idparm,idi,idj
c
c--------------------------------------------------------------------
c
c   
c
        open(31,file='skyinpgs.'//idparm//'.'//idi//'.'//idj)
        write(31,'(a/a/2i5,a,2(/4(1pg17.9),a),(/5(1pg17.9),a)
     &             /i4,0pf12.8,i4,f6.1,i4,a)')
     &   " 100000          ITMXT   max. total iterations",
     &   "'"//idi//":SHE '",
     &   108,158,"  nprot,nneut",
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
c
      return
      end

c-----extobsiso--------------------------------------------------------
c
      subroutine extobsiso(par,workdv,nparmx,npr1,nne1,npr2,nne2)
c
c     Prepare output for computing correlations with isotopic shifts.
c       npr1,npr2  proton numbers for the two nuclei
c       nne1,nne2  neutron numbers for the two nuclei
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
      dimension   par(kparmx),pars(kparmx),workdv(kparmx,1)
      real*8  work2(kparmx,kparmx) !,vect(kparmx,kparmx),eig(kparmx)
      logical tsymderiv
      parameter (tsymderiv=.true.)
      logical tsecord,tstsecord
      parameter (tsecord=.false.,tstsecord=.false.)
      character*1 asec
      CHARACTER(26) :: nametyp
c
c
c-----------------------------------------------------------------------
c
c      WRITE(2,'(a)') ' extobsiso: '
c
      IF(npr1==npr2 .AND. nne1==nne2) RETURN
c
c     loop over variation of parameters
c
      IF(tsymderiv) THEN
        nparlow = -nparmx
      ELSE
        nparlow = 0
      END IF
      pars = par
      do i=nparlow,nparmx
      IF(tsecord) THEN
        nparm2=abs(i)
        nparlw2=-abs(i)
      ELSE
        nparm2=0
        nparlw2=0
      END IF
      do j=nparlw2,nparm2
        if(j.ne.0 .and. i.eq.0) CYCLE
        if(j.ne.0 .and. j.eq.-i) CYCLE
        if(j.eq.0 .and. i.eq.0) par(1:nparmx) = pars(1:nparmx)
        if(j.eq.0 .and. i.ne.0) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        if(j.ne.0 .and. j.ne.i) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))/2D0
     &             + SIGN(1,j)*parmnatfac*parmnat(1:nparmx,abs(j))/2D0
        if(j.ne.0 .and. j.eq.i) par(1:nparmx) = pars(1:nparmx)
     &             - SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        call attach(par,jp)
c
        nprot  = npr1
        nneut  = nne1
        nnumb  = nprot+nneut
        anumb  = nnumb
        nmaxp  = 0
        nmaxn  = 0
        call hafosu
!        observ = etotal !rmsc**2
        observ = rmsc !**2
!        IF(nisotyp>1) observ=rmsc                  ! mod
        nprot  = npr2
        nneut  = nne2
        nnumb  = nprot+nneut
        anumb  = nnumb
        nmaxp  = 0
        nmaxn  = 0
        call hafosu
        IF(nisotyp==6) THEN
           observ=rmsc-observ
        ELSE
           observ = rmsc-observ !rmsc**2-observ
!           observ = etotal-observ
        END IF
c
c       collect array of deviations
c
        CALL  collectw(i,j,observ,observsav,workdv,work2,tsymderiv)
c
        par(1:nparmx) = pars(1:nparmx)
      enddo
      enddo

c
c     accumulate extrapolation errors by folding with 'zarray'.
c
      errobs   = 0.0
c
      if(nparmx.gt.0) then
        do i=1,nparmx
          errobs  = workdv(i,1)*workdv(i,1) + errobs
        enddo
        errobs  = sqrt(errobs)/parmnatfac
      endif
c
c     final print
c
      IF(tsecord) THEN
        asec='#'
      ELSE
        asec=' '
      END IF
      WRITE(nametyp,'(a,2(i2,a,i3,a))') 
     &   'isoshift (',npr1,'/',nne1,')-(',npr2,'/',nne2,')'
      write(02,'(a,2f9.4)')
     &     nametyp//': value,uncertainty=',observsav,errobs
        write(82,'(a)') asec//nametyp
        write(82,'(1pg24.16)')  observsav
        write(82,'(100(1pg24.16))') workdv(1:nparmx,1)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx)
      CALL flush(2)
      CALL flush(82)
c
c     estimate curvature contributions
c
      IF(tsecord) CALL estim_curv(workdv,work2,nparmx,tstsecord)
c
      return
      end

c-----extobsiso--------------------------------------------------------
c
      subroutine extobsiso_old(par,workdv,nparmx,nisotyp)
c
c     Prepare output for computing correlations with isotopic shifts.
c     'nisotyp' selects the isotopic shift to be computed, see
c     tabulated values in 'nprottyp' and 'nneuttyp'.
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
      dimension   par(kparmx),pars(kparmx),workdv(kparmx,1)
      real*8  work2(kparmx,kparmx) !,vect(kparmx,kparmx),eig(kparmx)
      logical tsymderiv
      parameter (tsymderiv=.true.)
      logical tsecord,tstsecord
      parameter (tsecord=.false.,tstsecord=.false.)
      character*1 asec
      REAL(8) :: nprottyp(2,5),nneuttyp(2,5)
      DATA nprottyp/20,20,50,50,82,82,62,50,88,82/
      DATA nneuttyp/24,28,74,82,132,126,82,82,126,126/
!      DATA nprottyp/20,20,20,20,20,16,20,20,20,20/
!      DATA nneuttyp/28,20,24,28,16,20,16,20,18,20/
!      DATA nprottyp/20,20,20,20,20,20,26,26,82,82/
!      DATA nneuttyp/28,20,24,28,32,28,32,28,132,126/
!      DATA nprottyp/20,20,28,22,28,24,28,26,16,20/
!      DATA nneuttyp/28,20,22,28,24,28,26,28,20,16/
      CHARACTER(18) :: nametyp(5)
!      DATA nametyp/'isoshift 48Ca-40Ca','mirror 50Ni-50Ti  ',
!     &             'mirror 52Ni-52Cr  ','mirror 54Ni-54Fe  ',
!     &             'mirror 36S-36Ca   '/
!      DATA nametyp/'isoshift 48Ca-40Ca','isoshift 44Ca-48Ca',
!     &             'isoshift 52Ca-48Ca','isoshift 58Fe-54Fe',
!     &             'isoshift 214-208Pb'/
!      DATA nametyp/'isoshift 48Ca-40Ca','isoshift 44Ca-48Ca',
!     &             'isoshift  36Ca-36S','isoshift 36Ca-40Ca',
!     &             'isoshift 38Ca-40Ca'/
      DATA nametyp/'isoshift 44Ca-48Ca','isoshift 124-132Sn',
     &             'isoshift 214-208Pb','isoshift N144-132 ',
     &             'isoshift 214-208N '/
c
c
c-----------------------------------------------------------------------
c
c      WRITE(2,'(a)') ' extobsiso: '
c
c     loop over variation of parameters
c
      IF(tsymderiv) THEN
        nparlow = -nparmx
      ELSE
        nparlow = 0
      END IF
      pars = par
      do i=nparlow,nparmx
      IF(tsecord) THEN
        nparm2=abs(i)
        nparlw2=-abs(i)
      ELSE
        nparm2=0
        nparlw2=0
      END IF
      do j=nparlw2,nparm2
        if(j.ne.0 .and. i.eq.0) CYCLE
        if(j.ne.0 .and. j.eq.-i) CYCLE
        if(j.eq.0 .and. i.eq.0) par(1:nparmx) = pars(1:nparmx)
        if(j.eq.0 .and. i.ne.0) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        if(j.ne.0 .and. j.ne.i) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))/2D0
     &             + SIGN(1,j)*parmnatfac*parmnat(1:nparmx,abs(j))/2D0
        if(j.ne.0 .and. j.eq.i) par(1:nparmx) = pars(1:nparmx)
     &             - SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        call attach(par,jp)
c
        nprot  = nprottyp(1,nisotyp)
        nneut  = nneuttyp(1,nisotyp)
        nnumb  = nprot+nneut
        anumb  = nnumb
        nmaxp  = 0
        nmaxn  = 0
        call hafosu
        observ = rmsc**2
!        IF(nisotyp>1) observ=rmsc                  ! mod
        nprot  = nprottyp(2,nisotyp)
        nneut  = nneuttyp(2,nisotyp)
        nnumb  = nprot+nneut
        anumb  = nnumb
        nmaxp  = 0
        nmaxn  = 0
        call hafosu
        IF(nisotyp==6) THEN
           observ=rmsc-observ
        ELSE
           observ = rmsc**2-observ
        END IF
c
c       collect array of deviations
c
        CALL  collectw(i,j,observ,observsav,workdv,work2,tsymderiv)
c
        par(1:nparmx) = pars(1:nparmx)
      enddo
      enddo

c
c     accumulate extrapolation errors by folding with 'zarray'.
c
      errobs   = 0.0
c
      if(nparmx.gt.0) then
        do i=1,nparmx
          errobs  = workdv(i,1)*workdv(i,1) + errobs
        enddo
        errobs  = sqrt(errobs)/parmnatfac
      endif
c
c     final print
c
      IF(tsecord) THEN
        asec='#'
      ELSE
        asec=' '
      END IF
      write(02,'(a,2f8.3)')
     &     nametyp(nisotyp)//': value,uncertainty=',
     &     observsav,errobs
        write(82,'(a)') asec//nametyp(nisotyp)
        write(82,'(1pg24.16)')  observsav
        write(82,'(100(1pg24.16))') workdv(1:nparmx,1)
        IF(tsecord) 
     &     write(82,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx)
c
c     estimate curvature contributions
c
      IF(tsecord) CALL estim_curv(workdv,work2,nparmx,tstsecord)
c
      return
      end
c-----extobs3gap--------------------------------------------------------
c
      subroutine extobs3gap(par,nparmx,nprmin,nnemin,nprmax,nnemax,
     &                      ifrms)
c
c     Prepare output for computing correlations with 3-point gaps
c     along isotopic or isotonic chains. The chain is selected by
c        nprmin  = minimum proton number
c        nnemin  = minimum neutron number
c        nprmax  = maximum proton number
c        nnemax  = maximum neutron number
c        ifrms   = 'ifrms=1' switches to odd-even staggering of radii
c                  'ifrms=0' switches to odd-even staggering of energies
c                  'ifrms=3' switches to even-even staggering of radii
c                  'ifrms=2' switches to even-even staggering of energies
c      Typical chains are
c        nprmin,nnemin,nprmax,nnemax  =  20,22,20,32, 
c                                        50,56,50,80,
c                                        82,120,82,124,
c                                        14,20,22,20,
c                                        34,50,44,50,
c                                        52,82,66,82.
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
      REAL(8) ::   par(kparmx)
      REAL(8) ::   pars(kparmx)
      REAL(8), ALLOCATABLE :: workdv(:,:),work2(:,:,:),esav(:),errobs(:)
      REAL(8), ALLOCATABLE :: observsav(:)
      INTEGER, ALLOCATABLE :: nprotchain(:),nneutchain(:)
!      real*8  work2(kparmx,kparmx) !,vect(kparmx,kparmx),eig(kparmx)
      LOGICAL, PARAMETER :: tsymderiv=.true.
      LOGICAL, PARAMETER :: tsecord=.false.,tstsecord=.false.
      INTEGER, PARAMETER :: nmagic(6)=(/20,28,50,82,92,126/)
      INTEGER :: nstr=2        ! stride
      LOGICAL :: topen83
      CHARACTER(1) :: asec
      CHARACTER(3) :: id2
      CHARACTER(7) :: nametyp(2,0:3)
c
c
c-----------------------------------------------------------------------
c
c
      IF(ifrms<0 .OR. ifrms>3) STOP 'EXTOBS3GAP: invalid IFRMS'
      nametyp(1,0) = 'p-3gap1'
      nametyp(2,0) = 'n-3gap1'
      nametyp(1,1) = 'p-3rms1'
      nametyp(2,1) = 'n-3rms1'
      nametyp(1,2) = 'p-3gap2'
      nametyp(2,2) = 'n-3gap2'
      nametyp(1,3) = 'p-3rms2'
      nametyp(2,3) = 'n-3rms2'
      nstr=ifrms/2+1
!      WRITE(2,'(a,i5)') '3GAP: nstr=',nstr
c
c     prepare output file
c
      INQUIRE(83,opened=topen83)
      IF(.NOT.topen83) THEN     
        open(39,status='scratch')
        write(39,'(i4)') nint(parmnatfac*100+1000)
        rewind(39)
        read(39,'(1x,a3)') id2
        close(39)
        open(83,file='error-gaps.'//id2)
        WRITE(2,*) 'EXTOBS3GAP: file 83 opened'
      END IF
c
c     prepare bookkeeping fields
c
      lenprot =  nprmax-nprmin
      lenneut =  nnemax-nnemin
      IF(lenprot*lenneut.NE.0) 
     &  STOP 'EXTOBS3GAP: not an isotopic or isotonic  chain'
      lenchain = MAX(lenprot,lenneut)
!      WRITE(2,'(a,2i4)') 'EXTOBS3GAP: chains=',lenprot,lenneut
      ALLOCATE(workdv(nparmx,0:lenchain))
      ALLOCATE(work2(nparmx,nparmx,0:lenchain))
      ALLOCATE(esav(-nstr:lenchain+nstr),errobs(0:lenchain))
      ALLOCATE(observsav(0:lenchain))
      ALLOCATE(nprotchain(-nstr:lenchain+nstr))
      ALLOCATE(nneutchain(-nstr:lenchain+nstr))
      IF(lenprot+lenneut==0) THEN
        isochain=0
        DO icheck=1,6
          IF(nmagic(icheck)==nprmin) THEN
            isochain=2
            DO i=-nstr,lenchain+nstr
              nprotchain(i) = nprmin
              nneutchain(i) = i+nnemin
            END DO
            EXIT
          ELSE IF(nmagic(icheck)==nnemin) THEN
            isochain=1
            DO i=-nstr,lenchain+nstr
              nprotchain(i) = i+nprmin
              nneutchain(i) = nnemin
            END DO
          END IF
        END DO
        IF(isochain==0) STOP '3GAP: undefined direction'
      ELSE IF(lenprot>0) THEN
        isochain=1
        DO i=-nstr,lenchain+nstr
          nprotchain(i) = i+nprmin
          nneutchain(i) = nnemin
        END DO
      ELSE
        isochain=2
        DO i=-nstr,lenchain+nstr
          nprotchain(i) = nprmin
          nneutchain(i) = i+nnemin
        END DO
      END IF
!      WRITE(2,'(a,i3)') 'isochain=',isochain
c
c     loop over variation of parameters
c
      IF(tsymderiv) THEN
        nparlow = -nparmx
      ELSE
        nparlow = 0
      END IF
      pars = par
      do i=nparlow,nparmx
      IF(tsecord) THEN
        nparm2=abs(i)
        nparlw2=-abs(i)
      ELSE
        nparm2=0
        nparlw2=0
      END IF
      do j=nparlw2,nparm2
        if(j.ne.0 .and. i.eq.0) CYCLE
        if(j.ne.0 .and. j.eq.-i) CYCLE
        if(j.eq.0 .and. i.eq.0) par(1:nparmx) = pars(1:nparmx)
        if(j.eq.0 .and. i.ne.0) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        if(j.ne.0 .and. j.ne.i) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))/2D0
     &             + SIGN(1,j)*parmnatfac*parmnat(1:nparmx,abs(j))/2D0
        if(j.ne.0 .and. j.eq.i) par(1:nparmx) = pars(1:nparmx)
     &             - SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        call attach(par,jp)
c
        DO io=-nstr,lenchain+nstr,nstr
          nprot = nprotchain(io)
          nneut = nneutchain(io) 
          nnumb = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          CALL blocking
          IF(MOD(ifrms,2)==1) THEN
            esav(io) = rmsc
          ELSE IF(MOD(ifrms,2)==0) THEN
            esav(io) = etotal   
          END IF
!          WRITE(2,'(a,4i4,1pg13.5)')'HaFo :',i,io,nprot,nneut,esav(io)
!          WRITE(6,'(a,4i4,1pg13.5)')'HaFo :',i,io,nprot,nneut,esav(io)
        END DO
c
c       collect array of deviations
c
        DO io=0,lenchain,nstr
!          observ = 0.5D0*(esav(io-1)-2D0*esav(io)+esav(io+1))
          observ = 0.5D0*(esav(io-nstr)-2D0*esav(io)+esav(io+nstr))
!          WRITE(2,'(a,i4,1pg13.5)')'observ:',io,observ
          CALL  collectw(i,j,observ,observsav(io),workdv(:,io),
     &                   work2(:,:,io),tsymderiv)
        END DO
c
        par(1:nparmx) = pars(1:nparmx)
      enddo
      enddo

c
c     accumulate extrapolation errors by folding with 'zarray'.
c
      errobs   = 0D0
      if(nparmx.gt.0) then
        do io=0,lenchain,nstr
          do i=1,nparmx
            errobs(io)  = workdv(i,io)*workdv(i,io) + errobs(io)
          enddo
          errobs(io)  = sqrt(errobs(io))/parmnatfac
        enddo
      endif
c
c     final print
c
      IF(tsecord) THEN
        asec='#'
      ELSE
        asec=' '
      END IF
      DO io=0,lenchain,nstr
!      WRITE(2,*) 'EXTOBS3GAP: before writing. io=',io
      write(02,'(a,i3,a,i3,a,f10.5,f10.5)')
     &     nametyp(isochain,ifrms)//'(',nprotchain(io),',',
     &     nneutchain(io),'): value,uncertainty=',
     &     observsav(io),errobs(io)
        write(83,'(a,i3,a,i3,a)') asec//nametyp(isochain,ifrms)//'(',
     &                  nprotchain(io),',',nneutchain(io),')'
        write(83,'(1pg24.16)')  observsav(io)
        write(83,'(100(1pg24.16))') workdv(1:nparmx,io)
        IF(tsecord) 
     &     write(83,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,io)
        CALL flush(2)
        CALL flush(83)
      END DO
      DEALLOCATE(workdv)
      DEALLOCATE(work2)
      DEALLOCATE(esav,errobs)
      DEALLOCATE(observsav)
      DEALLOCATE(nprotchain,nneutchain)
c
c     estimate curvature contributions
c
!      IF(tsecord) CALL estim_curv(workdv,work2,nparmx,tstsecord)
c
      return
      end
c-----extobs5gap--------------------------------------------------------
c
      subroutine extobs5gap(par,nparmx,nprmin,nnemin,nprmax,nnemax,
     &                      ifrms)
c
c     Prepare output for computing correlations with 3-point gaps
c     along isotopic or isotonic chains. The chain is selected by
c        nprmin  = minimum proton number
c        nnemin  = minimum neutron number
c        nprmax  = maximum proton number
c        nnemax  = maximum neutron number
c        ifrms   = 'ifrms=1' switches to odd-even staggering of radii
c                  'ifrms=0' switches to odd-even staggering of energies
c      Typical chains are
c        nprmin,nnemin,nprmax,nnemax  =  20,22,20,32, 
c                                        50,56,50,80,
c                                        82,120,82,124,
c                                        14,20,22,20,
c                                        34,50,44,50,
c                                        52,82,66,82.
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
      REAL(8) ::   par(kparmx)
      REAL(8) ::   pars(kparmx)
      REAL(8), ALLOCATABLE :: workdv(:,:),work2(:,:,:),esav(:),errobs(:)
      REAL(8), ALLOCATABLE :: observsav(:)
      INTEGER, ALLOCATABLE :: nprotchain(:),nneutchain(:)
!      real*8  work2(kparmx,kparmx) !,vect(kparmx,kparmx),eig(kparmx)
      LOGICAL, PARAMETER :: tsymderiv=.true.
      LOGICAL, PARAMETER :: tsecord=.false.,tstsecord=.false.
      INTEGER, PARAMETER :: nmagic(6)=(/20,28,50,82,92,126/)
      INTEGER :: nstr=1        ! stride
      LOGICAL :: topen83
      CHARACTER(1) :: asec
      CHARACTER(3) :: id2
      CHARACTER(6) :: nametyp(2,0:1)
c
c
c-----------------------------------------------------------------------
c
c
      IF(ifrms<0 .OR. ifrms>1) STOP 'EXTOBS5GAP: invalid IFRMS'
      nametyp(1,0) = 'p-5gap'
      nametyp(2,0) = 'n-5gap'
      nametyp(1,1) = 'p-5rms'
      nametyp(2,1) = 'n-5rms'
c
c     prepare output file
c
      INQUIRE(83,opened=topen83)
      IF(.NOT.topen83) THEN     
        open(39,status='scratch')
        write(39,'(i4)') nint(parmnatfac*100+1000)
        rewind(39)
        read(39,'(1x,a3)') id2
        close(39)
        open(83,file='error-gaps.'//id2)
        WRITE(2,*) 'EXTOBS5GAP: file 83 opened'
      END IF
c
c     prepare bookkeeping fields
c
      lenprot =  nprmax-nprmin
      lenneut =  nnemax-nnemin
      IF(lenprot*lenneut.NE.0) 
     &  STOP 'EXTOBS5GAP: not an isotopic or isotonic  chain'
      lenchain = MAX(lenprot,lenneut)
!      IF(lenchain<5) STOP 'EXTOBS5GAP: chain too short for 5-gap'
!      WRITE(2,*) 'EXTOBS5GAP: lenchain=',lenchain
      ALLOCATE(workdv(nparmx,0:lenchain))
      ALLOCATE(work2(nparmx,nparmx,0:lenchain))
      ALLOCATE(esav(-2*nstr:lenchain+2*nstr),errobs(0:lenchain))
      ALLOCATE(observsav(0:lenchain))
      ALLOCATE(nprotchain(-2*nstr:lenchain+2*nstr))
      ALLOCATE(nneutchain(-2*nstr:lenchain+2*nstr))
      IF(lenprot+lenneut==0) THEN
        isochain=0
        DO icheck=1,6
          IF(nmagic(icheck)==nprmin) THEN
            isochain=2
            DO i=-2*nstr,lenchain+2*nstr
              nprotchain(i) = nprmin
              nneutchain(i) = nstr*i+nnemin
            END DO
          ELSE IF(nmagic(icheck)==nnemin) THEN
            isochain=1
            DO i=-2*nstr,lenchain+2*nstr
              nprotchain(i) = nstr*i+nprmin
              nneutchain(i) = nnemin
            END DO
          END IF
        END DO
        IF(isochain==0) STOP '5GAP: undefined direction'
      ELSE IF(lenprot>0) THEN
        isochain=1
        DO i=-2*nstr,lenchain+2*nstr
          nprotchain(i) = i+nprmin
          nneutchain(i) = nnemin
        END DO
      ELSE
        isochain=2
        DO i=-2*nstr,lenchain+2*nstr
          nprotchain(i) = nprmin
          nneutchain(i) = i+nnemin
        END DO
      END IF
c
c     loop over variation of parameters
c
      IF(tsymderiv) THEN
        nparlow = -nparmx
      ELSE
        nparlow = 0
      END IF
      pars = par
      do i=nparlow,nparmx
      IF(tsecord) THEN
        nparm2=abs(i)
        nparlw2=-abs(i)
      ELSE
        nparm2=0
        nparlw2=0
      END IF
      do j=nparlw2,nparm2
        if(j.ne.0 .and. i.eq.0) CYCLE
        if(j.ne.0 .and. j.eq.-i) CYCLE
        if(j.eq.0 .and. i.eq.0) par(1:nparmx) = pars(1:nparmx)
        if(j.eq.0 .and. i.ne.0) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        if(j.ne.0 .and. j.ne.i) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))/2D0
     &             + SIGN(1,j)*parmnatfac*parmnat(1:nparmx,abs(j))/2D0
        if(j.ne.0 .and. j.eq.i) par(1:nparmx) = pars(1:nparmx)
     &             - SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        call attach(par,jp)
c
        DO io=-2*nstr,lenchain+2*nstr,nstr
          nprot = nprotchain(io)
          nneut = nneutchain(io) 
          nnumb = nprot+nneut
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
!          WRITE(2,'(a,4i4)') 'HaFo for:',i,io,nprot,nneut
          CALL blocking
          IF(MOD(ifrms,2)==1) THEN
            esav(io) = rmsc
          ELSE IF(MOD(ifrms,2)==0) THEN
            esav(io) = etotal   
          END IF
        END DO
c
c       collect array of deviations
c
        DO io=0,lenchain,nstr
!          observ = 0.5D0*(esav(io-1)-2D0*esav(io)+esav(io+1))
          observ = 0.125D0*(esav(io-2*nstr)-4D0*esav(io-nstr)
     &              +6D0*esav(io)-4D0*esav(io+nstr)+esav(io+2*nstr))
          CALL  collectw(i,j,observ,observsav(io),workdv(:,io),
     &                   work2(:,:,io),tsymderiv)
        END DO
c
        par(1:nparmx) = pars(1:nparmx)
      enddo
      enddo

c
c     accumulate extrapolation errors by folding with 'zarray'.
c
      errobs   = 0D0
      if(nparmx.gt.0) then
        do io=0,lenchain,nstr
          do i=1,nparmx
            errobs(io)  = workdv(i,io)*workdv(i,io) + errobs(io)
          enddo
          errobs(io)  = sqrt(errobs(io))/parmnatfac
        enddo
      endif
c
c     final print
c
      IF(tsecord) THEN
        asec='#'
      ELSE
        asec=' '
      END IF
      DO io=0,lenchain,nstr
!      WRITE(2,*) 'EXTOBS5GAP: before writing. io=',io
      write(02,'(a,i3,a,i3,a,f10.5,f10.5)')
     &     nametyp(isochain,ifrms)//'(',nprotchain(io),',',
     &     nneutchain(io),'): value,uncertainty=',
     &     observsav(io),errobs(io)
        write(83,'(a,i3,a,i3,a)') asec//nametyp(isochain,ifrms)//'(',
     &                  nprotchain(io),',',nneutchain(io),')'
        write(83,'(1pg24.16)')  observsav(io)
        write(83,'(100(1pg24.16))') workdv(1:nparmx,io)
        IF(tsecord) 
     &     write(83,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,io)
        CALL flush(83)
        CALL flush(2)
      END DO
      DEALLOCATE(workdv)
      DEALLOCATE(work2)
      DEALLOCATE(esav,errobs)
      DEALLOCATE(observsav)
      DEALLOCATE(nprotchain,nneutchain)
c
c     estimate curvature contributions
c
!      IF(tsecord) CALL estim_curv(workdv,work2,nparmx,tstsecord)
c
      return
      end

c-----extobssgap--------------------------------------------------------
c
      subroutine extobssgap(par,nparmx,nprmin,nnemin,nprmax,nnemax)
c
c     Prepare output for computing correlations with spectral u*v-gaps
c     along isotopic or isotonic chains. The chain is selected by
c        nprmin  = minimum proton number
c        nnemin  = minimum neutron number
c        nprmax  = maximum proton number
c        nnemax  = maximum neutron number
c      Typical chains are
c        nprmin,nnemin,nprmax,nnemax  =  20,22,20,32, 
c                                        50,56,50,80,
c                                        82,120,82,124,
c                                        14,20,22,20,
c                                        34,50,44,50,
c                                        52,82,66,82.
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
      REAL(8) ::   par(kparmx)
      REAL(8) ::   pars(kparmx)
      REAL(8), ALLOCATABLE :: workdv(:,:),work2(:,:,:),errobs(:)
      REAL(8), ALLOCATABLE :: observsav(:)
      INTEGER, ALLOCATABLE :: nprotchain(:),nneutchain(:)
      INTEGER, PARAMETER :: nmagic(6)=(/20,28,50,82,92,126/)
!      real*8  work2(kparmx,kparmx) !,vect(kparmx,kparmx),eig(kparmx)
      LOGICAL, PARAMETER :: tsymderiv=.true.
      LOGICAL, PARAMETER :: tsecord=.false.,tstsecord=.false.
      CHARACTER(1) :: asec
      CHARACTER(6), PARAMETER :: nametyp(2)= (/'p-sgap','n-sgap'/)
c
c
c-----------------------------------------------------------------------
c
c     prepare bookkeeping fields
c
      lenprot =  (nprmax-nprmin)/2
      lenneut =  (nnemax-nnemin)/2
      IF(lenprot*lenneut.NE.0) 
     &  STOP 'EXTOBSSGAP: not an isotopic or isotonic  chain'
      lenchain = MAX(lenprot,lenneut)
      ALLOCATE(workdv(nparmx,0:lenchain))
      ALLOCATE(work2(nparmx,nparmx,0:lenchain))
      ALLOCATE(errobs(0:lenchain))
      ALLOCATE(observsav(0:lenchain))
      ALLOCATE(nprotchain(0:lenchain),nneutchain(0:lenchain))
      IF(lenprot+lenneut==0) THEN
        isochain=0
        DO icheck=1,6
          IF(nmagic(icheck)==nprmin) THEN
            isochain=2
            DO i=0,lenchain
              nprotchain(i) = nprmin
              nneutchain(i) = 2*i+nnemin
            END DO
          ELSE IF(nmagic(icheck)==nnemin) THEN
            isochain=1
            DO i=0,lenchain
              nprotchain(i) = 2*i+nprmin
              nneutchain(i) = nnemin
            END DO
          END IF
        END DO
        IF(isochain==0) STOP '5GAP: undefined direction'
      ELSE IF(lenprot>0) THEN
        isochain=1
        DO i=0,lenchain
          nprotchain(i) = 2*i+nprmin
          nneutchain(i) = nnemin
        END DO
      ELSE
        isochain=2
        DO i=0,lenchain
          nprotchain(i) = nprmin
          nneutchain(i) = 2*i+nnemin
        END DO
      END IF
c
c     loop over variation of parameters
c
      IF(tsymderiv) THEN
        nparlow = -nparmx
      ELSE
        nparlow = 0
      END IF
      pars = par
      do i=nparlow,nparmx
      IF(tsecord) THEN
        nparm2=abs(i)
        nparlw2=-abs(i)
      ELSE
        nparm2=0
        nparlw2=0
      END IF
      do j=nparlw2,nparm2
        if(j.ne.0 .and. i.eq.0) CYCLE
        if(j.ne.0 .and. j.eq.-i) CYCLE
        if(j.eq.0 .and. i.eq.0) par(1:nparmx) = pars(1:nparmx)
        if(j.eq.0 .and. i.ne.0) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        if(j.ne.0 .and. j.ne.i) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))/2D0
     &             + SIGN(1,j)*parmnatfac*parmnat(1:nparmx,abs(j))/2D0
        if(j.ne.0 .and. j.eq.i) par(1:nparmx) = pars(1:nparmx)
     &             - SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        call attach(par,jp)
c
c       collect array of deviations
c
        DO io=0,lenchain
          nprot = nprotchain(io)
          nneut = nneutchain(io) 
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          CALL blocking
          IF(isochain==1) THEN
            observ = deltp
          ELSE
            observ = deltn
          END IF
          CALL  collectw(i,j,observ,observsav(io),workdv(:,io),
     &                   work2(:,:,io),tsymderiv)
        END DO
c
        par(1:nparmx) = pars(1:nparmx)
      enddo
      enddo

c
c     accumulate extrapolation errors by folding with 'zarray'.
c
      errobs   = 0D0
      if(nparmx.gt.0) then
        do io=0,lenchain
          do i=1,nparmx
            errobs(io)  = workdv(i,io)*workdv(i,io) + errobs(io)
          enddo
          errobs(io)  = sqrt(errobs(io))/parmnatfac
        enddo
      endif
c
c     final print
c
      IF(tsecord) THEN
        asec='#'
      ELSE
        asec=' '
      END IF
      DO io=0,lenchain
      write(02,'(a,i3,a,i3,a,2f8.3)')
     &     nametyp(isochain)//'(',nprotchain(io),',',nneutchain(io),
     &     '): value,uncertainty=',observsav(io),errobs(io)
        write(83,'(a,i3,a,i3,a)') asec//nametyp(isochain)//'(',
     &                  nprotchain(io),',',nneutchain(io),')'
        write(83,'(1pg24.16)')  observsav(io)
        write(83,'(100(1pg24.16))') workdv(1:nparmx,io)
        IF(tsecord) 
     &     write(83,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,io)
        CALL flush(83)
        CALL flush(2)
      END DO
      DEALLOCATE(workdv)
      DEALLOCATE(work2)
      DEALLOCATE(errobs)
      DEALLOCATE(observsav)
      DEALLOCATE(nprotchain,nneutchain)
c
c     estimate curvature contributions
c
!      IF(tsecord) CALL estim_curv(workdv,work2,nparmx,tstsecord)
c
      return
      end

c-----extobssgap--------------------------------------------------------
c
      subroutine extobs_centrdep(par,nparmx,nprmin,nnemin,nprmax,nnemax,
     &                           nrobs)

c
c     Prepare output for computing correlations with central depression
c     and related observables along isotopic or isotonic chains. 
c     The chain is selected by
c        nprmin  = minimum proton number
c        nnemin  = minimum neutron number
c        nprmax  = maximum proton number
c        nnemax  = maximum neutron number
c     Typical chains are
c        nprmin,nnemin,nprmax,nnemax  =  20,22,20,32, 
c                                        50,56,50,80,
c                                        82,120,82,124,
c                                        14,20,22,20,
c                                        34,50,44,50,
c                                        52,82,66,82.
c     The observable chosen is
c        nrobs  = 1 --> proton central depression 
c        nrobs  = 2 --> neutron central depression 
c        nrobs  = 3 --> formerly Coulomb energy, now charge rms radius
c
      include 'params.inc'
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'grid.inc'
      include 'comfit.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'extobs.inc'
      include 'constn.inc'
c
c     'workdv' is field for derivatives with respect to parameter-variation.
c
      REAL(8) ::   par(kparmx)
      REAL(8) ::   pars(kparmx)
      REAL(8), ALLOCATABLE :: workdv(:,:),work2(:,:,:),errobs(:)
      REAL(8), ALLOCATABLE :: observsav(:)
      INTEGER, ALLOCATABLE :: nprotchain(:),nneutchain(:)
!      real*8  work2(kparmx,kparmx) !,vect(kparmx,kparmx),eig(kparmx)
      LOGICAL, PARAMETER :: tsymderiv=.true.
      LOGICAL, PARAMETER :: tsecord=.false.,tstsecord=.false.
      CHARACTER(1) :: asec
      CHARACTER(9), PARAMETER :: 
     &    nametyp(3)= (/'rho_C(0) ','rho_n(0) ','c-centrdp'/)
!     &    nametyp(3)= (/'rho_p(0) ','rho_n(0) ','E_coulomb'/)
!     CHARACTER(9), PARAMETER :: 
!    &    nametyp(3)= (/'p-centrep','n-centrep','E_coulomb'/)
      LOGICAL :: topen84
      CHARACTER(3) :: id2,id2old='   '
c
c
c-----------------------------------------------------------------------
c
!      IF(nrobs==3) RETURN
c
c     prepare output file
c
      open(39,status='scratch')
      write(39,'(i4)') nint(parmnatfac*100+1000)
      rewind(39)
      read(39,'(1x,a3)') id2
      close(39)
      IF(id2.NE.id2old) CLOSE(84)
      id2old=id2
      INQUIRE(84,opened=topen84)
      IF(.NOT.topen84) THEN     
        open(84,file='error-rhocentr.'//id2)
!        open(84,file='error-centrdep.'//id2)
        WRITE(2,*) 'EXTOBS_CENTRDEP: file 84 opened'
      END IF
c
c     prepare bookkeeping fields
c
      lenprot =  (nprmax-nprmin)/2
      lenneut =  (nnemax-nnemin)/2
      IF(lenprot*lenneut.NE.0) 
     &  STOP 'EXTOBS_CENTRDEP: not an isotopic or isotonic  chain'
      lenchain = MAX(lenprot,lenneut)
      ALLOCATE(workdv(nparmx,0:lenchain))
      ALLOCATE(work2(nparmx,nparmx,0:lenchain))
      ALLOCATE(errobs(0:lenchain))
      ALLOCATE(observsav(0:lenchain))
      ALLOCATE(nprotchain(0:lenchain),nneutchain(0:lenchain))
      IF(lenprot>0) THEN
        isochain=1
        DO i=0,lenchain
          nprotchain(i) = 2*i+nprmin
          nneutchain(i) = nnemin
        END DO
      ELSE
        isochain=2
        DO i=0,lenchain
          nprotchain(i) = nprmin
          nneutchain(i) = 2*i+nnemin
        END DO
      END IF
c
c     loop over variation of parameters
c
      IF(tsymderiv) THEN
        nparlow = -nparmx
      ELSE
        nparlow = 0
      END IF
      pars = par
      do i=nparlow,nparmx
      IF(tsecord) THEN
        nparm2=abs(i)
        nparlw2=-abs(i)
      ELSE
        nparm2=0
        nparlw2=0
      END IF
      do j=nparlw2,nparm2
        if(j.ne.0 .and. i.eq.0) CYCLE
        if(j.ne.0 .and. j.eq.-i) CYCLE
        if(j.eq.0 .and. i.eq.0) par(1:nparmx) = pars(1:nparmx)
        if(j.eq.0 .and. i.ne.0) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        if(j.ne.0 .and. j.ne.i) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))/2D0
     &             + SIGN(1,j)*parmnatfac*parmnat(1:nparmx,abs(j))/2D0
        if(j.ne.0 .and. j.eq.i) par(1:nparmx) = pars(1:nparmx)
     &             - SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        call attach(par,jp)
c
c       collect array of deviations
c
        DO io=0,lenchain
          nprot = nprotchain(io)
          nneut = nneutchain(io) 
          anumb  = nnumb
          nmaxp  = 0
          nmaxn  = 0
          ngrid  = 81 !0 !81
          CALL blocking
          IF(nrobs==1) THEN
            observ = rho(1) !15.3*(rdmsp/rsecp-1.0) ! rho(1) !
          ELSE IF(nrobs==2) THEN
            observ = rho(ngridp) ! 15.3*(rdmsn/rsecn-1.0) ! rho(ngridp)  !
          ELSE IF(nrobs==3) THEN
            observ = rmsc ! 15.3*(rdmsc/rsecc-1.0) !ecoul
          ELSE
            STOP 'EXTOBS_CENTRDEP: invalid observable'
          END IF
          CALL  collectw(i,j,observ,observsav(io),workdv(:,io),
     &                   work2(:,:,io),tsymderiv)
        END DO
c
        par(1:nparmx) = pars(1:nparmx)
      enddo
      enddo

c
c     accumulate extrapolation errors by folding with 'zarray'.
c
      errobs   = 0D0
      if(nparmx.gt.0) then
        do io=0,lenchain
          do i=1,nparmx
            errobs(io)  = workdv(i,io)*workdv(i,io) + errobs(io)
          enddo
          errobs(io)  = sqrt(errobs(io))/parmnatfac
        enddo
      endif
c
c     final print
c
      IF(tsecord) THEN
        asec='#'
      ELSE
        asec=' '
      END IF
      DO io=0,lenchain
      write(02,'(a,i3,a,i3,a,2f8.3,i4)')
     &     nametyp(nrobs)//'(',nprotchain(io),',',nneutchain(io),
     &     '): value,uncertainty=',observsav(io),errobs(io)
        write(84,'(a,i3,a,i3,a)') asec//nametyp(nrobs)//'(',
     &                  nprotchain(io),',',nneutchain(io),')'
        write(84,'(1pg24.16)')  observsav(io)
        write(84,'(100(1pg24.16))') workdv(1:nparmx,io)
        IF(tsecord) 
     &     write(84,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx,io)
      END DO
      DEALLOCATE(workdv)
      DEALLOCATE(work2)
      DEALLOCATE(errobs)
      DEALLOCATE(observsav)
      DEALLOCATE(nprotchain,nneutchain)
c
c     estimate curvature contributions
c
!      IF(tsecord) CALL estim_curv(workdv,work2,nparmx,tstsecord)
c
      return
      end



      include 'asymt_sub.f'      
