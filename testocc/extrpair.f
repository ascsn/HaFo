c-----extrpair----------------------------------------------------------
c
      program extrpair
c
c     to extract effective pairing  energy from file 'for027' of SHF output
c
      parameter (KSTATE=500)
c
      character*80 aline
      character*10 title
      dimension wgt(KSTATE),energ(KSTATE),mult(KSTATE),jtwice(KSTATE),
     &          langmo(KSTATE)    
c
c-----------------------------------------------------------------------
c
      write(6,'(2a)')
     & 'force    Z    etotal   epair_f  e_corr  epair_c',
     & ' e_corrected'
      do i0=1,100
      do i=1,19999
        read(5,'(a)',err=99,end=99) aline
        if(aline(1:4).eq.' * *') then
          read(5,'(a)',err=99,end=99) title
          if(title(1:4).eq.'free') then
            do iskip=1,60
              read(5,'(1x)')
            enddo
          else
c            write(6,'(a)') title
            read(5,'(1x)')
            read(5,'(17x,i5)') iz
c            write(6,'(a,i5)') ' nr. protons=',iz
c            write(6,'(a)') '  2*j   l  l_eff   s.p.energy '
          endif
        elseif(aline(1:7).eq.' diffe.') then
          read(5,'(1x)')
          read(5,'(15x,f12.5,31x,e12.4)') etotal,epair
        elseif(aline(3:17).eq.'proton spectrum') then
          read(5,'(1x)')
          do nst=1,999
            read(5,*,err=19,end=19) j2,l,energy,wgtact,rms,multi
            do j=2,multi/2
              read(5,*,err=99,end=99) j2p,lp,energp,wgtp,rmsp,multip
            enddo
            if(nst.gt.KSTATE) stop ' too many states in input '
            energ(nst)     = energy
            wgt(nst)       = wgtact
            mult(nst)      = multi
            jtwice(nst)    = j2
            langmo(KSTATE) = l
            nstmax     = nst
c            write(6,'(1x,2i4,f8.3,1pg13.5)') 
c     &         j2,l,l+0.4*(j2/2.0-l),energy
          enddo
 19       continue
          epcorr = 0.0
          partnm = 0.0
          pardyn = 0.0
          goal   = iz*1.0
c          write(6,'(a)') 
c     &    '  nstate  energy   weight   shellwgt   multiplicity '
          do nst=1,nstmax
            paradd = mult(nst)
            pardyn = mult(nst)*wgt(nst) + pardyn
            if(partnm.gt.goal) then
              shlwgt = 0.0
            elseif(partnm+paradd.lt.goal) then
              shlwgt = 1.0
            else
              shlwgt = (goal-partnm)/(paradd)
            endif
            partnm = paradd + partnm
            epcadd = energ(nst)*mult(nst)*(wgt(nst)-shlwgt)
            epcorr = epcadd + epcorr
c            write(6,'(1x,i5,f8.2,2f8.5,i5,2f7.2,f8.3)')
c     &       nst,energ(nst),wgt(nst),shlwgt,mult(nst),
c     &       partnm,pardyn,epcadd
          enddo
c          write(6,'(a/1x,5(1pg14.6))')
c     &    '    etotal    epair_f   e_corr    epair_c    e_corrected ',
c     &    etotal,epair,epair+epcorr,epcorr,etotal-epair-epcorr        
          write(6,'(a,i4,f10.3,3f8.3,f10.3)')
     &    title(1:7),iz,etotal,epair,epair+epcorr,epcorr,
     &    etotal-epair-epcorr
        endif
      enddo
      enddo
c
 99   stop
      end



