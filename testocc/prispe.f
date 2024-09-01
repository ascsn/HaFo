c-----prispe------------ part of the hartree-fock package ----------------------
c
      subroutine prispe
c
      include 'params.inc'
c
c     printing header and single particle energies on file 'for027'
c
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'result.inc'
      include 'constn.inc'
c
      dimension esort(kmaxn),isort(kmaxn)
      data iop27/0/
c
c------------------------------------------------------------------------------
c
      if(iop27.eq.0) then
        open(unit=27,file='for027')
        iop27 = 1
      endif
c
c     some global infomation as a header
c
      ifile = 27
      ainfac=one/anumb**(one/3.0)
      write(ifile,'(14x,a,4(/a,4(1pg13.5)))')
     &   '     rms          rdms       surface        halo',
     &   '   mass radii:',rms,rdms,surf,
     &    rms-sqrt(0.6*rdms*rdms+3.0*surf*surf),
     &   ' charge radii:',rmsc,rdmsc,surfc,
     &    rmsc-sqrt(0.6*rdmsc*rdmsc+3.0*surfc*surfc),
     &   '              ',rmsc*ainfac,rdmsc*ainfac
       write(ifile,'(2(a,4(1pg13.5)/),2(a,2(1pg13.5)/))')
     &   ' proton radii:',rmsp,rdmsp,surfp,
     &    rmsp-sqrt(0.6*rdmsp*rdmsp+3.0*surfp*surfp),
     &   ' neutr. radii:',rmsn,rdmsn,surfn,
     &    rmsn-sqrt(0.6*rdmsn*rdmsn+3.0*surfn*surfn),
     &   '              ',rmsn*ainfac,rdmsn*ainfac,
     &   ' diffe. radii:',rmsn-rmsp,rdmsn-rdmsp
      write(ifile,'(a,1(1pg16.9),2(a,1pg12.5)/a,1pg12.5/
     &              2(a,1pg12.4),
     &              2(/a,1pg12.5,a,2(1pg12.5)))')
     &       '  total energy=',etotal,
     &               '   zpe=',zpe,'  epair=',epair,
     &       '           e/a=',etotal/anumb,
     &       '  without ZPE: E=',etotal+zpe,'  E/A=',(etotal+zpe)/anumb,
     &       '    total ekin=',ekin,'     ekin_p/n=',ekinpr,ekinne,
     &       '        ekin/A=',ekin/anumb,
     &       ' rel. ekin_p/n=',ekinpr/aprot,ekinne/aneut
      write(ifile,'(a,3(1pg16.9))')
     &        '  shell energies=',eshl(1),eshl(2),eshl(1)+eshl(2)
      sqa = sqrt(anumb)
      write(ifile,'(2(a,2(1pg12.4)))')
     &  ' pairing:  delta_p/n',deltp*sqa,deltn*sqa,
     &           '  gapfc_p/n',gpp*anumb,gpn*anumb,
     &  '          efermi_p/n',efep,efen,
     &           '  epair_p/n',epp,epn
c
c     print single particle energies and occupation weights for protons
c
c     sort energies
c
      do n=1,nmaxp
        esort(n) = spenrg(n)
        isort(n) = n
      enddo
      do n=1,nmaxp
        emin = esort(n)
        nmin = n
        do n2=n,nmaxp
          if(esort(n2).lt.emin) then
            emin        = esort(n2)
            esort(n2)   = esort(imin)
            esort(nmin) = emin
            imin        = isort(n2)
            isort(n2)   = isort(nmin)
            isort(nmin) = imin
            nmin        = n2
          endif
        enddo
      enddo
c
c     print
c
      write(ifile,'(/a/2a)')
     & '  proton spectrum:',
     & '  2j    l    energy       weight   rms     multipl.   '
      do n2=1,nmaxp
        n  = isort(n2)
        js=jsp(n)
        ll=js/2
        if(2*ll.eq.js) js=js-1
        npl = (n-1)*ngrid
        sum0 = 0.0
        sum2 = 0.0
        do ir=1,ngrid
          r  = (ir-1)*rstep
          r2 = r*r
          vol= wf(npl+ir)*wf(npl+ir)*r2
          sum0 = vol + sum0
          sum2 = vol*r2 + sum2
        enddo
        sum2 = sqrt(abs(sum2/max(sum0,xsmall)))
        mmult = min((js+1)/2,multmx)                   ! half the multiplicity
        if(spenrg(n).ne.zero .or. sum2.ne.zero) then
          do imult=1,mmult
            write(ifile,'(1x,i3,i5,3g13.5,i5)') js,ll,
     &         spenrg(n),weight(n),sum2,js+1
          enddo
        endif
      enddo
      backspace 27
      write(27,'(18x,i3,i5,3g13.5,i5)') -js,ll,
     &       spenrg(n),weight(n),sum2,js+1
c
c
c     and now printout for the neutrons
c
c     sort energies
c
      do n=1,nmaxn
        esort(n) = spenrg(n+nmaxp)
        isort(n) = n
      enddo
      do n=1,nmaxn
        emin = esort(n)
        nmin = n
        do n2=n,nmaxn
          if(esort(n2).lt.emin) then
            emin        = esort(n2)
            esort(n2)   = esort(imin)
            esort(nmin) = emin
            imin        = isort(n2)
            isort(n2)   = isort(nmin)
            isort(nmin) = imin
            nmin        = n2
          endif
        enddo
      enddo
c
c     print
c
      write(ifile,'(/a/2a)')
     &  '  neutron spectrum:',
     &  '  2j    l    energy       weight   rms     multipl.   '
      do nn=1,nmaxn
        n=isort(nn)+nmaxp
        js=jsp(nn)
        ll=js/2
        if(2*ll.eq.js) js=js-1
        npl = (n-1)*ngrid
        sum0 = 0.0
        sum2 = 0.0
        do ir=1,ngrid
          r  = (ir-1)*rstep
          r2 = r*r
          vol= wf(npl+ir)*wf(npl+ir)*r2
          sum0 = vol + sum0
          sum2 = vol*r2 + sum2
        enddo
        sum2 = sqrt(abs(sum2/max(sum0,xsmall)))
        if(2*ll.eq.js) js=js-1
        mmult = min((js+1)/2,multmx)                   ! half the multiplicity
        if(spenrg(n).ne.zero .or. sum2.ne.zero) then
          do imult=1,mmult
            write(ifile,'(1x,i3,i5,3g13.5,i5)') js,ll,
     &       spenrg(n),weight(n),sum2,js+1
          enddo
        endif
      enddo
      backspace 27
      write(27,'(18x,i3,i5,3g13.5,i5)') -js,ll,
     &       spenrg(n),weight(n),sum2,js+1
c
      return
c
      end
