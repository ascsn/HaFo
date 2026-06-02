c-----blocking------------------------------------------
c
      SUBROUTINE blocking()
c
c     Routine to compute odd nuclei within spherical 
c     blocking approximation. It checks particle number.
c     If both are even, a standard HF calculation is run.
c     If one odd nucleon number is identified, all possible
c     blocked s.p.states near the Fermi surface are 
c     computed and the one with lowest total energy  is
c     taken as ground state.
c
      include 'params.inc'
      include 'constn.inc'
      include 'wavfun.inc'
      include 'result.inc'
      include 'skypar.inc'
      REAL(8), PARAMETER :: blockband0=25D0 
      REAL(8),ALLOCATABLE :: eblock(:),spsav(:),occsav(:),
     &           ebl2(:,:),rbl2(:,:),
     &           gapbl2(:,:,:),gwsav(:)
      LOGICAL,PARAMETER :: testprintb=.TRUE.,tprispe=.TRUE.
      LOGICAL :: tfirst=.TRUE.
!
!     check odd nucleon sort
!
      IF(MOD(nprot,2)==1 .AND. MOD(nneut,2)==1) THEN
        isoodd = 3
      ELSE IF(MOD(nprot,2)==1) THEN
        isoodd = 1
      ELSE IF(MOD(nneut,2)==1) THEN
        isoodd = 2
      ELSE
        isoodd =0
      END IF
!      WRITE(2,*) ' isoodd:',isoodd,tfirst
!      CALL flush(2) 
      IF(isoodd.NE.0 .AND.  !tfirst) THEN
     &   (tfirst .OR. .NOT.ALLOCATED(eblock))) THEN
        ALLOCATE(eblock(0:ksl2),spsav(ksl2),occsav(ksl2),
     &           ebl2(0:ksl2,0:ksl2),rbl2(0:ksl2,0:ksl2),
     &           gapbl2(2,0:ksl2,0:ksl2),gwsav(ksl2))
        tfirst=.FALSE.
        eblock = 0D0
        ebl2=0D0
!        WRITE(2,*) ' EBLOCK etc allocated:',ksl2,eblock(0)
!        CALL flush(2) 
      END IF
      blockband = blockband0*(nprot+nneut*1D0)**(-0.333333D0)
      nblockp=0
      nblockn=0
      IF(testprintb) 
     &  WRITE(6,'(a,i2)') 'BLOCKING: isoodd=',isoodd
c
c     first a run without blocking and subsequent selection of 
c     blocked states
c
!      IF(isoodd.NE.0) WRITE(20,'(a)') 'initial run for BLOCKING'
      CALL hafosu
      IF(isoodd==0) THEN
        OPEN(unit=962,file='result.000',form='unformatted')
!        OPEN(unit=962,file='result.000',form='formatted')
        WRITE(962) nprot,nneut,nmaxp,nmaxn,etotal,rdmsc,surfc,
     &              rmsc,rmsn,deltp,deltn,rhochr,spenrg
        WRITE(11,'(a,2i8,f10.3,f7.4,2(f7.3))')
     &      'final even even: ',
     &       nprot,nneut,etotal,rmsc,epp,epn !deltp,deltn
        WRITE(12,'(a,2i8,4(1pg13.5))') 'summary even-even:',
     &    nprot,nneut,etotal,rmsc,epp,epn
        WRITE(6,'(a,200(1pg13.5))') 'spenrg:',spenrg(1:nmaxto)
        CLOSE(962)
        RETURN
      END IF
!      WRITE(2,*) 'before:',ALLOCATED(eblock)
!      WRITE(2,*) 'before:',eblock(0),etotal
!      CALL flush(2)
      eblock(0) = etotal
      spsav=spenrg
      occsav=weight
      efepsav=efep
      efensav=efen
      IF(isoodd==0) THEN
        RETURN
      ELSE IF(isoodd==1) THEN
        efact =efepsav
        DO nblockp=1,nmaxp
          IF(ABS(spsav(nblockp)-efepsav)<blockband) THEN
            CALL hafosu
            eblock(nblockp) = etotal
            rbl2(nblockp,1) = rmsc
!            gapbl2(:,nblockp,1) = (/deltp,deltn/)       
!            gapbl2(:,nblockp,1) = (/-epp/deltp,-epn/deltn/)            
            gapbl2(:,nblockp,1) = (/epp,epn/)            
            occsav(nblockp) = v2occ(nblockp)
            gwsav(nblockp) = weight(nblockp)
            WRITE(6,'(a,2i3,a,4(1pg13.5))') 'proton BLOCKING: ',nblockp,
     &        jsp(nblockp),'/2',etotal,rmsc,deltap,deltan
            WRITE(11,'(a,2i3,a,4(1pg13.5))')'proton BLOCKING: ',nblockp,
     &        jsp(nblockp),'/2',etotal,rmsc,deltap,deltan
          END IF
        END DO
        WRITE(11,'(a,2i8,1pg13.5)')
     &    'Summary of blocking calculations: nprot,nneut,E_unblock=',
     &       nprot,nneut,eblock(0)
        WRITE(12,'(a,2i8,1pg13.5)')
     &    'Summary of blocking calculations: nprot,nneut,E_unblock=',
     &       nprot,nneut,eblock(0)
        eblockmin = +1D4
        nblockmin = 0
        DO nblockp=1,nmaxp
          IF(eblock(nblockp).NE.0D0) THEN
            WRITE(6,'(2i4,4(1pg13.4))') nblockp,jsp(nblockp),
     &      spsav(nblockp),occsav(nblockp),eblock(nblockp),efact
            WRITE(11,'(2i4,a,4(1pg14.6))') nblockp,jsp(nblockp),'/2',
     &     eblock(nblockp),rbl2(nblockp,1),v2occ(nblockp),gwsav(nblockp)
            WRITE(12,'(2i4,a,4(1pg14.5))') nblockp,jsp(nblockp),'/2',
     &      eblock(nblockp),rbl2(nblockp,1),gapbl2(:,nblockp,1)
!     &     eblock(nblockp),rbl2(nblockp,1),v2occ(nblockp),gwsav(nblockp)
            IF(eblock(nblockp)<eblockmin) THEN
              eblockmin = eblock(nblockp)
              nblockmin = nblockp
            END IF
          END IF
        END DO
        nblockp = nblockmin
        CALL hafosu
        WRITE(11,'(a,2i8,2i3,a,f9.3,f7.4,2(f7.3))')
     &      'final proton BLOCKING: ',
     &       nprot,nneut,nblockp,jsp(nblockp),'/2   .    . ',etotal,rmsc
     &       ,epp,epn  !deltp,deltn
        WRITE(12,'(a,2i8,2i3,a,f9.3,f7.4,2(f7.3))')
     &      'final proton BLOCKING: ',
     &       nprot,nneut,nblockp,jsp(nblockp),'/2   .    . ',etotal,rmsc
     &       ,epp,epn  !deltp,deltn
      ELSE IF(isoodd==2) THEN
        efact =efensav
        DO nblock=1,nmaxn
          nblockn = nblock + nmaxp
          IF(testprintb) 
     &      WRITE(6,'(a,2i3,2(1pg13.5))') 'before neutron BLOCKING: ',
     &      nblock,nblockn,spenrg(nblockn),efen
          IF(testprintb) 
     &      WRITE(11,'(a,2i3,2(1pg13.5))') 'before neutron BLOCKING: ',
     &      nblock,nblockn,spenrg(nblockn),efen
          IF(ABS(spsav(nblockn)-efensav)<blockband) THEN
            CALL hafosu
            eblock(nblockn) = etotal            
            rbl2(nblockn,1) = rmsc
!            gapbl2(:,nblockn,1) = 2D0*(/sumuvp,sumuvn/)            
            gapbl2(:,nblockn,1) = (/epp,epn/)            
            occsav(nblockn) = v2occ(nblockn)
            gwsav(nblockn) = weight(nblockn)
            IF(testprintb) 
     &        WRITE(6,'(a,i3,1pg13.5)') 'neutron BLOCKING: ',
     &        nblockn,etotal
            IF(testprintb) 
     &         WRITE(11,'(a,i3,3(1pg13.5))') 'neutron BLOCKING: ',
     &         nblockn,etotal,v2occ(nblockn),gwsav(nblockn)
          END IF
        END DO
        WRITE(11,'(a,2i8,1pg13.5)')
     &    'Summary of blocking calculations: nprot,nneut,E_unblock=',
     &       nprot,nneut,eblock(0)
        WRITE(12,'(a,2i8,1pg13.5)')
     &    'Summary of blocking calculations: nprot,nneut,E_unblock=',
     &       nprot,nneut,eblock(0)
        eblockmin = +1D4
        nblockmin = 0
        DO nblock=1,nmaxn
          nblockn = nblock + nmaxp
          IF(eblock(nblockn).NE.0D0) THEN
            WRITE(11,'(2i4,a,6(1pg14.6))') nblockn,jsp(nblock),'/2',
     &      spsav(nblockn),occsav(nblockn),eblock(nblockn),
     &      rbl2(nblockn,1),v2occ(nblockn),gwsav(nblockn)
            WRITE(12,'(2i4,a,6(1pg14.5))') nblockn,jsp(nblock),'/2',
     &      spsav(nblockn),occsav(nblockn),eblock(nblockn),
     &      rbl2(nblockn,1),gapbl2(:,nblockn,1)
!     &      rbl2(nblockn,1),v2occ(nblockn),gwsav(nblockn)
            IF(eblock(nblockn)<eblockmin) THEN
              eblockmin = eblock(nblockn)
              nblockmin = nblockn
            END IF
          END IF
        END DO
        nblockn = nblockmin
!        WRITE(20,'(a)') 'final neutron BLOCKING: '
        CALL hafosu
!        WRITE(11,'(a,4i3,a,4(1pg13.5))')
        WRITE(11,'(a,2i8,2i3,a,f9.3,f8.4,2(f7.3))')
     &      'final neutron BLOCKING: ',
     &       nprot,nneut,nblockn,jsp(nblockn-nmaxp),'/2 . . ',etotal,
     &       rmsc,epp,epn !,deltp,deltn
        WRITE(12,'(a,2i8,2i3,a,f9.3,f8.4,2(f7.3))')
     &      'final neutron BLOCKING: ',
     &       nprot,nneut,nblockn,jsp(nblockn-nmaxp),'/2 . . ',etotal,
     &       rmsc,gapbl2(:,nblockn,1) !,deltp,deltn
      ELSE IF(isoodd==3) THEN
!        efact =efepsav
        DO nblockp=1,nmaxp
          IF(ABS(spsav(nblockp)-efepsav)<blockband) THEN
            DO nblock=1,nmaxn
              nblockn = nblock + nmaxp
              IF(ABS(spsav(nblockn)-efensav)<blockband) THEN
                IF(testprintb) 
     &          WRITE(6,'(a,2i3,2(1pg13.5))')'before double BLOCKING: ',
     &          nblockp,nblock,spenrg(nblockp),spenrg(nblockn)
                IF(testprintb) 
     &          WRITE(11,'(a,2i3,2(1pg13.5))')'before doub. BLOCKING: ',
     &          nblockp,nblock,spenrg(nblockp),spenrg(nblockn)
                CALL hafosu
                ebl2(nblockp,nblockn) = etotal            
                rbl2(nblockp,nblockn) = rmsc            
                gapbl2(:,nblockp,nblockn) = (/deltp,deltn/)              
                WRITE(6,'(a,2(2i3,a),4(1pg13.5))') 'double BLOCKING: ',
     &           nblockp,jsp(nblockp),'/2 ',nblock,
     &           jsp(nblock),'/2 ',etotal,rmsc,deltp,deltn
                WRITE(11,'(a,2(2i3,a),4(1pg14.6))') 'double BLOCKING: ',
     &           nblockp,jsp(nblockp),'/2 ',nblock,
     &           jsp(nblock),'/2 ',etotal,rmsc,deltp,deltn
              END IF
            END DO
          END IF
        END DO
        WRITE(6,'(a)') 'summary double blocking:'
        WRITE(11,'(a)') 'summary double blocking:'
        emin=1D3
        DO nblockp=1,nmaxp; DO nblock=1,nmaxn
          nblockn = nblock + nmaxp
          IF(ebl2(nblockp,nblockn).NE.0D0) THEN
            WRITE(6,'(a,2i3,1pg14.6)') 'nblockp,nblockn,energy:',
     &        nblockp,nblockn,ebl2(nblockp,nblockn)
            WRITE(11,'(2(2i3,a),4(1pg14.6))') 
     &           nblockp,jsp(nblockp),'/2 ',nblockn-nmaxp,
     &           jsp(nblockn-nmaxp),'/2 ',ebl2(nblockp,nblockn),
     &           rbl2(nblockp,nblockn),gapbl2(:,nblockp,nblockn)
            IF(ebl2(nblockp,nblockn)<emin) THEN
              emin=ebl2(nblockp,nblockn)
              nblpmin=nblockp
              nblnmin=nblockn
            END IF
          END IF
        END DO; END DO        
        nblockp=nblpmin
        nblockn=nblnmin
        WRITE(11,'(2i8,a,2i3,a,f9.3,3f7.4)')
        WRITE(11,'(/1x/a/1x/a,2i8,2(2i3,a),f9.3,3f7.4)')
     &     '******************', 'final double BLOCKING: ',
     &    nprot,nneut,nblockp,jsp(nblockp),'/2 ',nblockn-nmaxp,
     &           jsp(nblockn-nmaxp),'/2 ',emin,rbl2(nblockp,nblockn),
     &           gapbl2(:,nblockp,nblockn)
        WRITE(12,'(/1x/a/1x/a,2i8,2(2i3,a),f9.3,f7.4,2(f7.3))')
     &     '******************', 'final double BLOCKING: ',
     &    nprot,nneut,nblockp,jsp(nblockp),'/2 ',nblockn-nmaxp,
     &           jsp(nblockn-nmaxp),'/2 ',emin,rbl2(nblockp,nblockn),
     &           gapbl2(:,nblockp,nblockn)
        CALL hafosu
      END IF
      IF(tprispe) THEN
        WRITE(13,'(2i8,87(1pg11.3))')
     &   nprot,nneut,epp,(spenrg(i),i=1,nmaxp)
        CALL flush(13)
        WRITE(14,'(2i8,86(1pg11.3))')
     &   nprot,nneut,epn,(spenrg(i),i=nmaxp+1,nmaxto)
        CALL flush(14)
      END IF
        OPEN(unit=962,file='result.000',form='unformatted')
!        OPEN(unit=962,file='result.000',form='formatted')
        WRITE(962) nprot,nneut,nmaxp,nmaxn,etotal,rdmsc,surfc,
     &              rmsc,rmsn,deltp,deltn,rhochr,spenrg
        CLOSE(962)
      RETURN
      END
