c-----map1738----------------------------------------------------
c
c     Maps input with 17 parameters (ioptsk=10) 
c     to 38 parameters (ioptsk=11).
c
      implicit double precision (a-h,o-z)
      parameter (kparmx1=17,kparmx2=38)
      parameter (kzerolines=14)
      real*8 par1(kparmx1),par2(kparmx2)
      real*8 dpar1(kparmx1),dpar2(kparmx2)
      integer jp1(kparmx1),jp2(kparmx2)
      character*128 inline
      character*6 zuord(kparmx2)
c
c-----------------------------------------------------------------
c
      call naming(zuord,11)
!      write(*,*) zuord
      par2 = 0D0
      dpar2 = 0D0
      open(unit=1,file='for001.17')
      do j=1,kzerolines-1
        read(1,'(a)') inline
        lenline = len_trim(inline)
        write(6,'(a)') inline(1:lenline)
      enddo
      read(1,*) NPARIN,NUMAX,IVAR,FREEP,ENDCON,IFTINP
      write(6,'(3(i2,a),f4.1,a,f7.4,a)') 
     &   NPARIN,',',NUMAX,',',IVAR,',',FREEP,',',ENDCON,
     &   ',3,,11,, NPARIN,NUMAX,IVAR,FREEP,ENDCON,IFTINP'
      do j=1,kparmx1
        read(1,*) par1(j),dpar1(j),jp1(j)
      enddo
c
c  old
c
!         riftm  = par(jpact(11))
!         rho0pr = par(jpact(12))      ! default 0.16
!         powerp = par(jpact(13))
!       exmass = par(jpact(14))        ! switch to exchange functional
c
c  new
c
!         riftm  = par(jpact(14))
!         rho0pr = par(jpact(18))      ! default 0.16
!         powerp = par(jpact(19))
c  remaining
!       1...10
!       b4p    = par(jpact(15))
!       pairfp = par(jpact(16))
!       pairfn = par(jpact(17))
      do j=1,kparmx1
        IF(jp1(j)==11) THEN
          jp2(j) = 14
        ELSE IF(jp1(j)==12) THEN
          jp2(j) = 18
        ELSE IF(jp1(j)==13) THEN
          jp2(j) = 19
        ELSE IF(jp1(j)==14 .and. J.ne.KPARMX1) THEN
         stop "EXMASS must be the last entry in input list"
        ELSE
          jp2(j) = jp1(j)
        END IF
        par2(j) = par1(j)
        dpar2(j) = dpar1(j)
      enddo
      do j=kparmx1-1,kparmx1+1
        par2(j) = 0D0
        dpar2(j) = 0D0
        jp2(j) = 12+j-kparmx1
      enddo
      do j=kparmx1+2,kparmx2
        par2(j) = 0D0
        dpar2(j) = 0D0
        jp2(j) = j
      enddo
c
      do j=1,kparmx2
        write(6,'(1pg20.11,a,1pg13.5,a,i2,a)') 
     &   par2(j),',',dpar2(j),',',jp2(j),','//zuord(jp2(j))
      enddo
c      
      do j=1,2
        read(1,'(a)') inline
        lenline = len_trim(inline)
        write(6,'(a)') inline(1:lenline)
      enddo
c
      stop
      end

c-----naming----------------------------------------------------------
c
      subroutine naming(zuord,ioptsk)
c
c     genererates array 'zuord' which carries the names of the
c     associated Skyrme parameters
c
c      zuord    = output array of names
c      iftinp   = input variant in 'for001'
c      ioptsk   = option for Skyrme force
c
!      include 'params.inc'
!      include 'constn.inc'
c      parameter (zero=0.0)
c
!      include 'comfit.inc'
      parameter (kparmx=38,iftinp=3)
c
      character*6 zuord(kparmx)
      character*6 namepool(17,0:2)
      character*6 optpool(11:38,0:11)
      data ((namepool(i,j),i=1,15),j=0,2)
     &    /'    b0','    b1','    b2','    b3','    t4',
     &     '   b0p','   b1p','   b2p','   b3p',' power',
     &     '  byuk',' byukp','  ayuk','exmass','   b4p',
     &     '    t0','    t1','    t2','    t3','    t4',
     *     '    x0','    x1','    x2','    x3',' power',
     *     '  tyuk','  xyuk','  ayuk','exmass','   b4p',
     &     '    c0','    b1','    b2','    c3','    t4',
     &     '   c0p','   b1p','   b2p','   c3p',' power',
     &     '  byuk',' byukp','  ayuk','exmass','   b4p'/
      data ((namepool(i,j),i=16,17),j=0,2)
     &    /'pairfp','pairfn','pairfp','pairfn','pairfp','pairfn'/
      data ((optpool(i,j),i=11,14),j=5,6)
     &    /'    c1','   c1p','    c4','   c4p',
     &     '    c1','   c1p','    c4','   c4p'/
      data ((optpool(i,j),i=11,14),j=8,8)
     &    /'b3tild','   b12','  b12p',' undef'/
      data ((optpool(i,j),i=11,14),j=9,9)
     &    /'    b7','   b7p','  pow7',' undef'/
      data ((optpool(i,j),i=11,14),j=10,10)
     &    /' riftm','rho0pr','powerp',' undef'/
      data ((optpool(i,j),i=14,14),j=11,11)
     &    /' riftm'/
      data ((optpool(i,j),i=18,38),j=11,11)
     &    /'rho0pr','powerp','b3tild','   b12','  b12p',
     &     '    b7','   b7p','  pow7',
     &     '    c1','   c1p','    c4','   c4p','ecstab',
     &     'alpinv','betinv','alpin2','betin2','riftm1',
     &     '   bl2','  bl2p','str_sq'/
c
c----------------------------------------------------------------------
c
      if(kparmx.lt.17)
     &  stop ' Fit code requires KPARMX >= 17'
      if(ioptsk.eq.11 .and. kparmx.lt.30)
     &  stop ' IOPTSK=11 requires KPARMX >= 30'
c
      do ipar=1,17
        zuord(ipar) = namepool(ipar,min(iftinp,2))
      enddo
      if(ioptsk.eq.5 .or. ioptsk.eq.6 .or. ioptsk.eq.8) then
        do ipar=11,14
          zuord(ipar) = optpool(ipar,ioptsk)
         enddo
      elseif(ioptsk.eq.9 .or. ioptsk.eq.10) then
        do ipar=11,14
          zuord(ipar) = optpool(ipar,ioptsk)
        enddo
      elseif(ioptsk.eq.11) then
        zuord(14) = optpool(14,ioptsk)
        do ipar=18,kparmx
          zuord(ipar) = optpool(ipar,ioptsk)
        enddo
      endif
c
      return
      end
