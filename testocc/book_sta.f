c-----book_sta----------------------------------------------------------
c
      program book_sta
c
c     generates bookkeeping of oscillator states for the spherical
c     H.F. package
c
      parameter (KMAXSH=25)
      parameter (KMAXJS=2*KMAXSH-3)
      parameter (KMAXN=(KMAXSH-1)/2)
      parameter (KDATA=KMAXN*KMAXJS)
      parameter (KDATFL=2*(KMAXSH+1)-10)
c
      dimension nsp(KMAXN,KMAXJS)
      dimension nmaxfl(2,0:KMAXSH)
      dimension nact(KMAXJS)
      dimension jsp(9999)
c
      data nact/KMAXJS*0/
      data nsp/KDATA*999/
      data nmaxfl/0,0, 1,2, 2,6, 3,8, 6,20,KDATFL*0/
      data jsp/9999*0/
c
c-----------------------------------------------------------------------
c
      nsp(1,1) = 1
      nsp(1,2) = 3
      nsp(1,3) = 2
      nsp(1,4) = 5
      nsp(1,5) = 4
      nsp(2,1) = 6
      nact(1)  = 2
      nact(2)  = 1
      nact(3)  = 1
      nact(4)  = 1
      nact(5)  = 1
      jsp(1)   = 1
      jsp(2)   = 3
      jsp(3)   = 2
      jsp(4)   = 5
      jsp(5)   = 4
      jsp(6)   = 1
      nstate = 6                                ! counter for nr. of state
      nnumb  = 20                               ! counter for part.number
      do nprinc=3,KMAXSH-2                      ! loop over principle q.n.
        do l=nprinc,0,-2                        ! orbital angular momentum 
          do jph=l+1,max(l,1),-1                ! total ang.mom. + 1/2
            js     = l+jph
            numsh  = jph+jph                    ! nr. states in shell
            nstate = 1 + nstate
            nnumb  = numsh + nnumb
            nact(js) = 1 + nact(js)
            nsp(nact(js),js) = nstate
            jsp(nstate) = js
            if(l.eq.nprinc .and. jph.eq.l+1) then
              nmaxfl(1,nprinc+2) = nstate
              nmaxfl(2,nprinc+2) = nnumb
              if(nprinc.eq.KMAXSH-2) goto 99
            endif
          enddo
        enddo
      enddo
c
 99   continue
      write(6,'(a)') ' NMAXFL:'
      write(6,101) nmaxfl
 101  format(20(5x,'&',5(1x,i4,',',i6,',')/))
      write(6,'(a)') ' NSP:'
      do js=1,KMAXJS
        write(6,102) (nsp(n,js),n=1,KMAXN)
 102    format(5x,'&',1x,20(i4,','))
      enddo
      write(6,'(a)') ' ZWJ:'
      write(6,103) (float(1+(js-1)/2)/6.283185307,js=1,KMAXJS)
 103  format(20(5x,'&',6(f10.7,',')/))
      write(6,'(a,i5)') ' JSP: up to state nr.',nstate
      write(6,104) (jsp(n),n=1,nstate)
 104  format(20(5x,'&',15(i3,',')/))
c
      stop
      end




