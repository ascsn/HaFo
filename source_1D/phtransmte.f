c-----dosph -------------------------------------------------------------------
c
      subroutine dosph()
c
      include 'params.inc'
c
c     determines the density of ph states in the valence space
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'result.inc'
!      include 'skypar.inc'
      include 'grid.inc'
!      include 'linres.inc'
      include 'constn.inc'
c
      REAL(8),PARAMETER :: xminrd=1.0D-9
      REAL(8),PARAMETER :: ephlow=0D0,ephup=50D0
      INTEGER,PARAMETER :: kmaxphe=500
      REAL(8),PARAMETER :: deleph=(ephup-ephlow)/kmaxphe
      REAL(8),PARAMETER :: sigphwid=0.4D0
      REAL(8),PARAMETER :: normphwid=1D0/(SQRT(3.141926D0)*sigphwid)
      REAL(8) :: dosphe(0:kmaxphe)

c
c--------------------------------------------------------------------
c
c
c     initialises the standard loop structure in the
c     linear response package. occupies the e.f.n. 10, 20 and 30.
c
c     isospin loop
c
      lrpa = 1
      alrpa=lrpa
      mlrpa=2*lrpa+1
      dosphe=0D0
      do 10 iso=1,2
      if(iso.eq.1) then
       nmax   = nmaxp
       nstack = 0
       isoadd = 0
       faciso = aneut/anumb
      else
       nmax   = nmaxn
       nstack = nmaxp*ngrid
       isoadd = ngrid
       faciso = -aprot/anumb
      endif
c
c     loop over the a-states
c
      do 20 na=1,nmax
      jsa   = jsp(na)
      la    = jsa/2
      xla   = la
      xja   = (jsa-1)/2+0.5
      wgt   = weight(na+nstack/ngrid)
!      if(wgt.lt.0.5d0) goto 20
      nadd  = nstack + (na-1)*ngrid
      subfac = spenrg(na+nstack/ngrid)
c
c     loop over the j,l to which they couple
c
      do 30 jdif=1,mlrpa
      xj     = xja+jdif-lrpa-1
      if(xj.lt.0.0d0) goto 30
      ll     = (xj-0.5)
      if(mod(lrpa+la+ll,2).ne.0) then
       ll     = (xj+0.5)
      endif
      xl     = ll
      if(ll.gt.lrpa+la .or. ll.lt.abs(la-lrpa)) go to 30
c
c     determine combined quantum number to 'xj', 'xl'.
c
      if(xj.gt.xl) then
       js     = xj+xj
      else
       js     = xj+xj+one
      endif
c
c     now all quantum numbers are provided:
c      'iso' and according stack pointers  'nstack','isoadd'
c            and according number of states 'nmax'.
c      'jsa' and according q.n. 'xja','xla','la'
c            and occupation weight 'wgt'
c            and single particle energy 'subfac'
c            and wave function stack pointer 'nadd'.
c      'jdif' and according q.n. 'xj','xl','ll' and 'js'
c

      redmat = redume(xj,alrpa,xja)
!      WRITE(519,'(a,4i4,4(1pg13.5))') 
!     &   'redmat:',jsa,js,jdif,na,redmat,xj,xja
      if(redmat**2.le.xminrd) goto 30
c
c     loop over the b-states
c
      do 40 nb=1,nmax
      jsb   = jsp(nb)
      if(jsb.ne.js) goto 40
      lb    = jsb/2
      xlb   = lb
      xjb   = (jsb-1)/2+0.5D0
      wgtb  = weight(nb+nstack/ngrid)
!      if(wgtb.ge.wgt) goto 40
      naddb = nstack + (nb-1)*ngrid
      subfcb = spenrg(nb+nstack/ngrid)
c
c     compute dipole matrix elements inclusive weights
c
      wgtph = ABS(wgtb-wgt)
      edif = abs(subfac-subfcb)
!     WRITE(519,'(a,2i4,2(1pg13.5))') 'edif:',jsb,nb,edif,wgtph
      if(edif.lt.ephup .and. edif.gt.ephlow
     &   .and. subfac.ne.zero .and. subfcb.ne.zero) then
!        acc   = zero
!        do ir=2,ngrid
!          acc   = wf(naddb+ir)*wf(nadd+ir)*((ir-1)*rstep) + acc
!        enddo
!        acc2 = (redmat*acc*faciso)**2*wgtph
        acc2 = (redmat*faciso)**2*wgtph
!        WRITE(519,'(4i4,4(1pg13.5))') 
!     &      na,jsa,nb,jsb,edif,wgtph,acc,acc2
        do ie=0,kmaxphe
          enph = ie*deleph
          dosphe(ie) = dosphe(ie) + 
     &           acc2*normphwid*EXP(-((enph-edif)/sigphwid)**2)          
        enddo
      endif
c
c
   40 continue

c
c     close the loops over the quantum number (complement of 'initdo').
c
   30 continue
   20 continue
   10 continue
c end include block

c
c     print results

      OPEN(unit=519,file='dosph.res')
      WRITE(519,'(a)') '# DoS for isovector dipole 1ph states',
     &     '#  A   Z    e_ph     DoS '
      nnumb = nprot+nneut
      do ie=0,kmaxphe
        enph = ie*deleph
        if(dosphe(ie)<1D-99) dosphe(ie)=0D0
          WRITE(519,'(2i4,f6.2,1pg13.5)') nnumb,nprot,enph,dosphe(ie)
        enddo
      WRITE(519,'(1x)')
      CLOSE(519)
c
      
      return
      end


c-----dosph -------------------------------------------------------------------
c
      subroutine dos()
c
      include 'params.inc'
c
c     determines the density of states 
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'result.inc'
!      include 'skypar.inc'
      include 'grid.inc'
!      include 'linres.inc'
      include 'constn.inc'
c
      REAL(8),PARAMETER :: elow=-50D0,eup=20D0
      INTEGER,PARAMETER :: kmaxe=700
      REAL(8),PARAMETER :: dele=(eup-elow)/kmaxe
      REAL(8),PARAMETER :: sigwid=0.4D0
      REAL(8),PARAMETER :: normwid=1D0/(SQRT(3.141926D0)*sigwid)
      REAL(8) :: dose(0:kmaxe)

c
c--------------------------------------------------------------------
c
c     isospin loop
c
      do iso=1,2
        if(iso.eq.1) then
          nmax   = nmaxp
          OPEN(519,file='dosprot.res')
        else
          nmax   = nmaxn
          OPEN(519,file='dosneut.res')
        endif
        ngw    = (iso-1)*nmaxp
        js     = jsp(n)
        phasew = zwj(js)
        dose = 0D0
        do ie=0,kmaxe
          ensp = ie*dele+elow
            do nshl=1,nmax
              nst = nshl+ngw
              IF(spenrg(nst)>0D0) THEN
                sigact = sigwid*2D0
              ELSE
                sigact = sigwid
              END IF
              dose(ie) = dose(ie) + phasew*
     &                   normwid*EXP(-((ensp-spenrg(nst))/sigact)**2)          
            enddo
        enddo
        nnumb = nprot+nneut
        do ie=0,kmaxe
          ensp = ie*dele+elow
          if(dose(ie)<1D-99) dose(ie)=0D0
          WRITE(519,'(2i4,2(1pg13.5))') nnumb,nprot,ensp,dose(ie)
        enddo
        WRITE(519,'(1x)')
        CLOSE(519)
      enddo
 
      return
      end
