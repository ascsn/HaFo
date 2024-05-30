c to recouple Fayans parameters
      include 'params.inc'
      include 'comfit.inc'
      include 'skypar.inc'
      CHARACTER(128) :: aline
      character(6) :: zuord(kparmx)
      REAL(8) :: a(kparmx),b(kparmx)

      CALL  naming_fy(zuord)
      DO i=1,17
        READ(5,'(a)') aline
        lena = len_trim(aline)
        WRITE(6,'(a)') aline(1:lena)
      END DO
      read(5,*) nparin,numax,ivar,freep,endcon,iftinp,ib2fix,ioptsk
      WRITE(6,'(3(i2,a1),f3.1,a,f7.5,a,3(i2,a))') 
     & nparin,',',numax,',',ivar,',',freep,',',endcon,',',2,',',
     & ib2fix,',',ioptsk,
     & ', nparin,numax,ivar,freep,endcon,iftinp,ib2fix,ioptsk'
      DO j=1,38
        READ(5,*) a(j),dparmo(j),iparva(j)
      END DO
      b=a
      do j=1,nparin
      do l=1,kparmx
      if(iparva(j).eq.l) jp(l)=j
      enddo
      enddo
      CALL recouple_fy(a,jp)
      DO j=1,38
        WRITE(6,'(g15.8,a,e12.5,a,i2,a,a,a,g15.8)')  
     &        a(j),',',dparmo(j),',',iparva(j),',',zuord(iparva(j)),
     &        ',',b(j)
      END DO
      WRITE(6,*) '-11,'      

      STOP
      END
 
c-----recouple_fy--------------------------------------------------------------------
c
      subroutine recouple_fy(par,jpact)
c
c     the explicit Fayans model parameters 
c     are attached according to the values on the field 'par'
c     and the pointers on the field 'jp' .
c
      include 'params.inc'
c
      include 'skypar.inc'
      include 'comfit.inc'
c
      dimension par(kparmx),jpact(kparmx)
c
      data zero/0D0/
c
c--------------------------------------------------------------------
c
    
      avp    = par(jpact(1))
      avm    = par(jpact(6))
      asp    = par(jpact(3))
      asm    = par(jpact(8))
      hsp = par(jpact(23))
      IF(hsp<-1D3) hsp=par(jpact(11))
      hsm = par(jpact(24))
      IF(hsm<-1D3) hsm=par(jpact(12))
        h1vp   = par(jpact(4))*avp
        par(jpact(4)) = h1vp
        h2vp   = par(jpact(11))*avp
        par(jpact(11)) = h2vp
        h1vm   = par(jpact(9))*avm
        par(jpact(9)) = h1vm
        h2vm   = par(jpact(12))*avm
        par(jpact(12)) = h2vm
        IF(ABS(asp)>1D-10) THEN
          hgradp = par(jpact(13))*asp
          par(jpact(13)) = hgradp
          IF(hsp>-1D3) THEN
            hsp = hsp*asp
            par(jpact(23)) = hsp
          END IF
        END IF
        IF(ABS(asm)>1D-10) THEN
          hgradm = par(jpact(14))*asm
          par(jpact(14)) = hgradm
          IF(hsm>-1D3) THEN
            hsm = hsm*asm
            par(jpact(24)) = hsm
          END IF
        END IF
c
      RETURN
      END
c-----naming_fy----------------------------------------------------------
c
      subroutine naming_fy(zuord)
c
c     genererates array 'zuord' which carries the names of the
c     associated Fayans parameters
c
c      zuord    = output array of names
c
      include 'params.inc'
      include 'constn.inc'
c
      include 'comfit.inc'
c
      character*6 zuord(kparmx)
c
      zuord(:) = '      '
c
       zuord(1)='avp   '
       zuord(4)='h1vp  '
       zuord(11)='h2vp  '
       zuord(6)='avm   '
       zuord(9)='h1vm  '
       zuord(12)='h2vm  '
       zuord(10)='sigmav'
       zuord(22)='rho0fy'
       zuord(3)='asp   '
       zuord(23)='hsp   '
       zuord(13)='hgradp'
       zuord(8)='asm   '
       zuord(24)='hsm   '
       zuord(14)='hgradm'
       zuord(5)='c0nabj'
       zuord(15)='c1nabj'
       zuord(2)='c0tau '
       zuord(7)='c1tau '
       zuord(36)='c0tens'
       zuord(37)='c1tens'
       zuord(21)='hcoul '
       zuord(16)='fxi   '
       zuord(17)='h1xi  '
       zuord(19)='hgradx'
       zuord(18)='rho0xi'
       zuord(20)='gammax'
       zuord(30)='ecstab'
c
      return
      end
