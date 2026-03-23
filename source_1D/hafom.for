c-----hafomn------------------------------------------------------
c
        program hafomn
c
      include 'params.inc'
c
c     program hafo: direct entry to the skyrme-hartree-fock package.
c     reads parameters and calls the package via 'hafosu'.
c
c************** sample input: file 05  ***************************
c
c 1,                                iwunit
c 'Skyrme 3'                        comment
c 00,0.3,                           ngrid,rstep
c -1,1.0E-6,1.0E-6,18,240,3,-1,     iwr,dabr,dabe,itpr,itmax,ifplot,idisc
c 0.30,50.0,0.3,-1.0,               x0damp,e0damp,addnew,wall
c 1,0,0,1,0,0,0,0,,                 if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold
c 1,                                ifform
c 0.64,0.0,0.0,0.0,0.0,
c 0.64,0.0,0.0,0.0,0.0,
c -1128.75,395,-95,14000.0,120.0,   t0,t1,t2,t3,t4
c 0.45,0.,0.,1.,,                   x0,x1,x2,x3,b4p
c 0.0,0.0,0.0,0.0,,                 t5,t6,t7,t8
c 0.0,0.0,0.0,0.0,                  x5,x6,x7,x8
c 0.0,0.0,0.0,                      tyuk,xyuk,ayuk
c 1.0,1.0,3,0,0,1,,                 power,acoul,ipair,pairfp,pairfn,ilnrsp,temprt
c 40,50,0,0,2,                      nprot,nneut,nmaxp,nmaxn,iparsh
c 0,                                ifgsc
c 1,,                            iwunit
c 'SkI1'                         comment
c 0,0.3,                         ngrid,rstep
c -1,1.E-7,1.E-7,18,280,3,-1,    iwr,dabr,dabe,itpr,itmax,ifplot,idisc
c 0.3,50.,0.3,-1.0,,             x0damp,e0damp,addnew,wall
c 1,0,1,1,1,0,0,0,,              if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold
c 1,,                            ifform
c 0.64,0.0,0.0,0.0,0.0,
c 0.64,0.0,0.0,0.0,0.0,
c -1913.620,439.808,2697.590,10592.300,124.260,      t0,t1,t2,t3,t4
c -0.954500,-5.782400,-1.287400,-1.561400,62.13,     x0,x1,x2,x3,b4p
c 0.0,0.0,0.0,0.0,,              t5,t6,t7,t8
c 0.0,0.0,0.0,0.0,               x5,x6,x7,x8
c 0.0,0.0,0.0,                   tyuk,xyuk,ayuk
c 0.2500000,1.0,2,0,0,0,,        power,acoul,ipair,pairfp,pairfn,ilnrsp,temprt
c 8,8,0,0,,,,,,,,,               nprot,nneut,nmaxp,nmaxn,iparsh
c 0,                             ifgsc
c 1,,                            iwunit
c 'SkI4'                         comment
c 0,0.3,                         ngrid,rstep
c -1,1.E-6,1.E-6,18,280,3,-1,    iwr,dabr,dabe,itpr,itmax,ifplot,idisc
c 0.3,50.,0.3,-1.0,,             x0damp,e0damp,addnew,wall
c 1,0,1,1,1,0,0,0,,              if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold
c 1,,                            ifform
c 0.64,0.0,0.0,0.0,0.0,
c 0.64,0.0,0.0,0.0,0.0,
c -1855.83,473.830,1006.85,9703.63,366.194,  t0,t1,t2,t3,t4
c 0.40512,-2.8891,-1.3252,1.1452,-180.351,   x0,x1,x2,x3,b4p
c 0.0,0.0,0.0,0.0,,              t5,t6,t7,t8
c 0.0,0.0,0.0,0.0,               x5,x6,x7,x8
c 0.0,0.0,0.0,                   tyuk,xyuk,ayuk
c 0.2500000,1.0,0,0,0,0,,        power,acoul,ipair,pairfp,pairfn,ilnrsp,temprt
c 82,126,0,0,0,                  nprot,nneut,nmaxp,nmaxn,iparsh
c 0,                             ifgsc
c 1,,                                  iwunit
c 'Skymre M* '                         comment
c 0,0.3,,,                             ngrid,rstep
c -1,1.E-6,1.E-6,18,280,3,-1,          iwr,dabr,dabe,itpr,itmax,ifplot,idisc
c 0.3,50.,0.3,-1.0,,                   x0damp,e0damp,addnew,wall
c 1,0,1,1,1,0,0,0,,                    if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold
c 1,,                                  ifform
c 0.64,0.0,0.0,0.0,0.0,
c 0.64,0.0,0.0,0.0,0.0,
c -2645.0,410.0,-135.0,15595.0,130.0,, t0,t1,t2,t3,t4
c  0.09,0.000000,0.000000,0.00,,,      x0,x1,x2,x3,b4p
c 0.0,0.0,0.0,0.0,,                    t5,t6,t7,t8
c 0.0,0.0,0.0,0.0,,                    x5,x6,x7,x8
c 0.0,0.0,0.0,                         tyuk,xyuk,ayuk
c 0.16666667,1.0,2,0,0,0,,             power,acoul,ipair,pairfp,pairfn,ilnrsp,temprt
c 50,72,0,0,2,                         nprot,nneut,nmaxp,nmaxn,iparsh
c 0,                                   ifgsc
c -1,,,,        iwunit<0 --> stop
c
c  other "old" skyrme forces are:
c  skyrme 1:
c -1057.3,235.9,-100.,14463.5,120.,       t0,t1,t2,t3,t4
c .56,0.,0.,1.,                           x0,x1,x2,x3
c  skyrme 2:
c -1169.9,585.6,-27.1,9331.1,105.,        t0,t1,t2,t3,t4
c .34,0.,0.,1.,                           x0,x1,x2,x3
c  skyrme 3:
c -1128.75,395.,-95.,14000.,120.,         t0,t1,t2,t3,t4
c .45,0.,0.,1.,                           x0,x1,x2,x3
c  skyrme 4:
c -1205.6,765.,35.,5000.,150.,            t0,t1,t2,t3,t4
c .05,0.,0.,1.,                           x0,x1,x2,x3
c  skyrme 5:
c -1248.29,970.56,107.22,0.,150.,         t0,t1,t2,t3,t4
c -.17,0.,0.,1.,                          x0,x1,x2,x3
c  skyrme 6:
c -1101.81,271.67,-138.33,17000.,115.,    t0,t1,t2,t3,t4
c .583,0.,0.,1.,                          x0,x1,x2,x3
c  in all cases is 'power'=1.0.
c
c******************************************************************
c
      logical*4  top11
c
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'rmsnuc.inc'
      include 'constn.inc'
!
      CHARACTER(20) :: ain
c
c------------------------------------------------------------------------------
c
c       open(unit=05,status='old',file='for005')
c
c     read switch 'iwunit' and open output unit accordingly
c
       iwuold=-1
       do 1000 irep=1,999
c
c      initialize values
c
       ngrid=0
       nmaxp=0
       nmaxn=0
       b4p   = zero
       t5 = zero
       t6 = zero
       t7 = zero
       t8 = zero
       t9 = zero
       ioptsk = 0
!
!      starts reading sequence
!
       read(05,'(a)',err=1001,end=1001)   ain
       WRITE(6,*) ain
!       rewind(05)
       BACKSPACE(05)
!       STOP
       IF(INDEX(ain,'&MAIN')>0 .OR. INDEX(ain,'&main')>0) THEN
         CALL nameinp()
!         WRITE(11,*) 'CUTFAC=',cutfac
       ELSE

       WRITE(6,*) 'in old branch'
       read(05,'(a)',err=1001,end=1001)   ain
       WRITE(6,*) ain
       read(ain,*,end=1001)   iwunit
       WRITE(6,*) iwunit
       if(iwunit.lt.0) stop ! " terminate at users request"
       if(iwunit.ne.iwuold) then
        if(iwuold.ne.-1) close(unit=06)
        if(iwunit.eq.0) then
         open(unit=06,file='for006',status='unknown')
c       else
c         open(unit=06,file='con')
        endif
       endif
       iwuold=iwunit
!       WRITE(6,'(a)') 'HAFOM:  next nucleus'
c
c      read all data
c
       read(05,*,err=992,end=992) title
c       write(6,*) title
       read(05,*,err=1001,end=1001)   ngrid,rstep
       read(05,*)  iwr,dabr,dabe,itpr,itmax,ifplot,idisc
c       write(6,*)  iwr,dabr,dabe,itpr,itmax,ifplot,idisc
       if(ifplot.gt.0) then
        inquire(11,opened=top11)
        if(.not.top11) then
         open(unit=11,status='unknown',file='for011')
        endif
c        write(11,'(/40a)')  (' *',i=1,40)
c        write(11,*) title
       endif
       read(05,*)  x0damp,e0damp,addnew,wall
c       write(6,*)  x0damp,e0damp,addnew,wall
       read(05,*)  ifex,iftm,ifenrg,ifrso,ifzpe,ifdfms,ifrhoc,
     &             ifenav,imfold,ioptsk
c       write(6,*)  ifex,iftm,ifenrg,ifrso,ifzpe,ifdfms,ifrhoc,
c     &             ifenav,imfold,ioptsk
!       if(ifzpe.eq.1) ifenrg=1
       read(05,*)  ifform
c       if(ifform.lt.0) stop ' this option IFFORM<0 not valid'
c       write(6,*)  ifform
       read(05,*)  apro0,apro11,apro12,bpro11,bpro12
       read(05,*)  aneu0,aneu11,aneu12,bneu11,bneu12
c       write(6,*)  apro0,apr011,apro12,bpro11,bpro12
c       write(6,*)  aneu0,aneu11,aneu12,bneu11,bneu12
       IF(ioptsk==13) THEN            ! Fayans functional
         tpairFy=.TRUE.
         READ(5,*) avp,h1vp,h2vp     
         WRITE(6,*) avp,h1vp,h2vp     
         READ(5,*) avm,h1vm,h2vm
         WRITE(6,*) avm,h1vm,h2vm
         READ(5,*) sigmavol,rho0fy
         WRITE(6,*) sigmavol,rho0fy
         READ(5,*) asp,hsp,hgradp 
         WRITE(6,*) asp,hsp,hgradp 
         READ(5,*) asm,hsm,hgradm 
         WRITE(6,*) asm,hsm,hgradm 
         READ(5,*) hcoul,sigmacoul
         WRITE(6,*) hcoul,sigmacoul
         READ(5,*) fxi,h1xi,hgradxi,gammaxi,rho0xi
         WRITE(6,*) fxi,h1xi,hgradxi,gammaxi,rho0xi
         READ(5,*) delfxi,ecstab,delh1xi,delhgradxi !,isovectormodel
         WRITE(6,*) delfxi,ecstab,delh1xi,delhgradxi,isovectormodel
!         delfxi=0D0
!
         READ(5,*) c0tau,c1tau 
         WRITE(6,*) c0tau,c1tau 
         READ(5,*) c0nabj,c1nabj,c0tens,c1tens
         WRITE(6,*) c0nabj,c1nabj,c0tens,c1tens
       ELSE IF(ioptsk==12) THEN
         read(05,*)     b0,b1,b2,b3,t4
         write(6,*)     'b0,b1,b2,b3=',b0,b1,b2,b3
         read(05,*)     b0p,b1p,b2p,b3p,b4p
         write(6,*)     'b0p,b1p,b2p,b3p=',b0p,b1p,b2p,b3p
       ELSE
         read(05,*)     t0,t1,t2,t3,t4
c         write(6,*)     t0,t1,t2,t3,t4
         read(05,*)     x0,x1,x2,x3,b4p
c         write(6,*)     x0,x1,x2,x3,b4p
       END IF
       if(ioptsk.eq.5 .or. ioptsk.eq.6 .or. ioptsk.eq.11) then
         read(05,*) c1,c4
         read(05,*) c1p,c4p
         write(*,*) ' c1,c1p,c4,c4p=',c1,c1p,c4,c4p
       endif
       if(ioptsk.eq.8 .or. ioptsk.eq.11) then
         read(05,*)  b3tilde,tte,tto
         read(05,*)  b12,b12p
         write(*,*) ' b3tilde,b12,b12p=', b3tilde,b12,b12p
       endif
       if(ioptsk.eq.9 .or. ioptsk.eq.11) then
         read(05,*)  b7,b7p
         read(05,*)  pow7
         write(*,*) ' b7,b7p,pow7=',b7,b7p,pow7
       endif
       if(ioptsk.eq.10) then
         read(05,*)  riftm,rho0pr,powerp,ecstab
         write(*,*) 'riftm,rho0pr,powerp=',riftm,rho0pr,powerp
       endif
       if(ioptsk.eq.11) then
         read(05,*)  riftm,rho0pr,powerp,ecstab,alpinv,betinv,
     &               alpin2,betin2
         write(*,'(a,8(1pg15.7))') 
     &     'riftm,rho0pr,powerp,ecstab,alpinv=',
     &      riftm,rho0pr,powerp,ecstab,alpinv,betinv,
     &       alpinv,betinv
       endif
       if(ioptsk.le.4 .or. ioptsk.eq.7) then
         read(05,*)     t5,t6,t7,t8,t9
c        write(6,*)     t5,t6,t7,t8,t9
         read(05,*)     x5,x6,x7,x8
c        write(6,*)     x5,x6,x7,x8
       endif
       if(t7.ne.zero) stop ' old parameter entry T7 disabled '
       IF(ioptsk.NE.13) THEN
         read(05,*)     tyuk,xyuk,ayuk
         write(*,*)     'tyuk,xyuk,ayuk=',tyuk,xyuk,ayuk
         ayukp  = ayuk                            ! restore global folding
         ayukm  = ayuk
         ayukls = 0.0D0  ! ayuk
       END IF
       read(05,*,err=991,end=991)
     &              power,acoul,ipair,pairfp,pairfn,ilnrsp,temprt,cutfac
     &              ,hfbwidth
       WRITE(6,*) acoul,ipair,cutfac,hfbwidth
       read(05,*)     nprot,nneut,nmaxp,nmaxn,xparsh
!       WRITE(6,'(a,4i4,1pg13.5)')'HAFOM:',nprot,nneut,nmaxp,nmaxn,xparsh
       read(05,*)     ifgsc
       iparsh = NINT(xparsh)

       END IF

!      prepare a couple of internal Fayans parameters
         IF(isovectormodel==0) THEN
           fxip=fxi
           h1xip=h1xi
           hgradxip=hgradxi
           fxin=fxi
           h1xin=h1xi
           hgradxin=hgradxi
         ELSE IF(isovectormodel==2) THEN
           fxip=fxi+delfxi
           h1xip=h1xi
           hgradxip=hgradxi
           fxin=fxi-delfxi
           h1xin=h1xi
           hgradxin=hgradxi
         ELSE IF(isovectormodel==1) THEN
           xipfac=(fxi+delfxi)/fxi
           fxip=fxi*xipfac
           h1xip=h1xi*xipfac
           hgradxip=hgradxi*xipfac
           delh1xi = h1xip-h1xi
           delhgradxi = hgradxip-hgradxi
           xinfac=(fxi-delfxi)/fxi
           fxin=fxi*xinfac
           h1xin=h1xi*xinfac
           hgradxin=hgradxi*xinfac
         ELSE IF(isovectormodel==3) THEN
           fxip=fxi+delfxi
           h1xip=h1xi+delh1xi
           hgradxip=hgradxi+delhgradxi
           fxin=fxi-delfxi
           h1xin=h1xi-delh1xi
           hgradxin=hgradxi-delhgradxi
         END IF

c
c      run the hartree-fock package
c
!       ecstab=0.3
!       call hafosu
       call blocking
c       write(21,*) '# call ORTHPOL: rstep,ngrid,langmo,maxorder=',
c     &   rstep,ngrid,1,8
c       call  orthopol(rho,rho(ngridp),rstep,ngridp-1,1,8)
c
 1000  continue
c
 1001  stop " error in reading "
  991  stop ' error in reading pairing '
  992  stop ' error in reading title '
       end
c      include 'hafbody.for'
c      include 'sumdum.for'
      include 'nameinp.f'
