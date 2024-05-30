c-----introbs----------------------------------------------------------
c
      subroutine introbs(par,nparmx)
c
      include 'params.inc'
c
c     Computes the parameter variations for all fit variables
c     and writes it to file for later use in correlation analysis.
c
      parameter (kobs=12)
      character*27 heatyp(12)
      character*3 id2
c
      include 'wavfun.inc'
c      include 'grid.inc'
      include 'skypar.inc'
      include 'numpar.inc'
      include 'result.inc'
      include 'comfit.inc'
c      include 'constn.inc'
c
      integer npoipr(5),npoine(6),maxdat(kobs)
      double precision par(kparmx),pars(kparmx)
      double precision varobs(kobs,kdat,0:nparmx)
c
      data zero,x10lw/0.0,1.e-10/
cold      data zero,x10up,x10lw/0.0,1.e+10,1.e-10/
      data heatyp/'  binding-energy',
     *            '  diffr.-radius',
     *            '  surface thickness',
     *            '  rms radius',
     *            '  l*s splitting proton 16O',
     *            '  l*s splitting neutron 16O',
     *            '  (rho(1.8)+rho(0.))/2.',
     *            '  neutron rms radius',
     *            '  gap-proton',
     *            '  gap-neutron',
     *            '  rho(7.2) or 2. l*s neutr.',
     *            '  m*/m symm. nuclear matter'/
c
c      pointer fields to the last filled shells in 208pb.
c
      data npoipr/12,13,14,15,16/
      data npoine/17,18,19,20,21,22/
c
c---------------------------------------------------------------------
c
      WRITE(2,*) 'INTROBS entered. parmnatfac=',parmnatfac
      pars = par
      do n=0,nparmx
        if(n.gt.0) par(1:nparmx) = pars(1:nparmx) 
     &            + parmnatfac*parmnat(1:nparmx,n)
        do i=1,npts                ! loop over nuclei
        call attach(par,jp)
c
c     mean-field calculation: first set the the parameters for
c        the actual nucleus which are then transferred to 'rmf'
c        via common's.
c
        nra    = iay(i)
        nnumb  = nra
        anumb  = nnumb
        nprot  = iprot(i)
        nneut  = nra-nprot
        nmaxp  = 0
        nmaxn  = 0
        ngrid  = 0
        if(ipair.lt.5 .or. kparmx.le.15) then
          pairfp = 0.0
          pairfn = 0.0
        endif
        call hafosu
        if(nprot.eq.108 .and. nneut.eq.156 .or.
     &     nprot.eq.94 .and. nneut.eq.152 .or.
     &     nprot.eq.66 .and. nneut.eq.98) then    ! special deformed
          etotal = energy_deform(nprot,nneut)
        endif
c
c      now extract the results which are transferred back to
c      this routine via common's and collect it on the appropriate
c      arrays
c
          varobs(lifit(1),i,n)=etotal
          if(ifitr>0) then
            varobs(lifit(2),i,n)=(rdmsc+0.011-0.323*surfc)/0.975
          else   
            varobs(lifit(2),i,n)=rdmsc
          endif
          varobs(lifit(3),i,n)=surfc
          varobs(lifit(4),i,n)=rmsc
c
          varobs(lifit(7),i,n)=(rhochr(7)+rhochr(1))/2.*100.
          varobs(lifit(8),i,n)=rmsn   !  -rmsp
          varobs(lifit(9),i,n)=deltp
          varobs(lifit(10),i,n)=deltn
          varobs(lifit(11),i,n)=rhochr(25)*100.
          if(nra.eq.16.and.nprot.eq.8) then
c           beitrag zu chi**2 von schalenenergie-differenzen
            soprot=spenrg(3)-spenrg(2)
            soneut=spenrg(nmaxp+3)-spenrg(nmaxp+2)
            varobs(lifit(5),i,n)=soprot
            varobs(lifit(6),i,n)=soneut
          elseif(nra.eq.40.and.nprot.eq.20) then
            varobs(lifit(5),i,n)=spenrg(5)-spenrg(4)
            varobs(lifit(6),i,n)=spenrg(nmaxp+5)-spenrg(nmaxp+4)
          elseif(nra.eq.48.and.nprot.eq.20) then
            varobs(lifit(5),i,n)=spenrg(5)-spenrg(4)
          elseif(nra.eq.90.and.nprot.eq.40) then
            varobs(lifit(5),i,n)=spenrg(12)-spenrg(11)
            varobs(lifit(6),i,n)=spenrg(nmaxp+12)-spenrg(nmaxp+11)
          elseif(nra.eq.132.and.nprot.eq.50) then
            varobs(lifit(5),i,n)=spenrg(10)-spenrg(9)                 ! 2p prot
            varobs(lifit(6),i,n)=spenrg(nmaxp+14)-spenrg(nmaxp+13)    ! 2d neut
          elseif(nra.eq.208.and.nprot.eq.82) then
            varobs(lifit(5),i,n)=spenrg(14)-spenrg(13)                ! 2d prot
            varobs(lifit(11),i,n)=spenrg(nmaxp+19)-spenrg(nmaxp+18)   ! 1fneut  
            varobs(lifit(6),i,n)=spenrg(nmaxp+21)-spenrg(nmaxp+20)    ! 3p neut
          elseif(nra.eq.207 .and. nprot.eq.82) then
c          beitrag zu chi**2 von level-width in 208pb
            xavera = zero
            eavera = zero
            e2aver = zero
            do m=1,5
              npl    = npoipr(m)
              xmult  = 2.0*float((jsp(npl)+1)/2)
              xavera = xmult + xavera
              eavera = xmult*spenrg(npl) + eavera
              e2aver = xmult*spenrg(npl)**2 +e2aver
            enddo
            widpr  = sqrt(e2aver/xavera-(eavera/xavera)**2)
            varobs(lifit(5),i,n) = widpr
            xavera = zero
            eavera = zero
            e2aver = zero
            do m=1,6
              npl    = npoine(m)+16
              xmult  = 2.0*float((jsp(npl)+1)/2)
              xavera = xmult + xavera
              eavera = xmult*spenrg(npl) + eavera
              e2aver = xmult*spenrg(npl)**2 +e2aver
            enddo
            widne  = sqrt(e2aver/xavera-(eavera/xavera)**2)
            varobs(lifit(6),i,n) = widne
          endif
c
          do l=1,kobs
            if(abs(yexp(lifit(l),i)).lt.x10lw) varobs(lifit(l),i,n)=0D0
          enddo
c
        enddo
      enddo
c
c     compute maximum nr. of entries per data block and
c     then print in contigous blocks
c
      maxdat = 0
      do k=1,ndats
        do i=1,npts
          if(lifit(k)<12 .and. varobs(lifit(k),i,0).NE.0D0)
     &      maxdat(k)=1+maxdat(k)
        enddo
      enddo
      open(82,status='scratch')
      write(82,'(i4)') nint(parmnatfac*100+1000)
      rewind(82)
      read(82,'(1x,a3)') id2 
      close(82)
      open(82,file='error-introbs.'//id2)
      write(82,'(a,i5,f8.2)') 'natural parameters: nparmx=',
     &     nparmx,parmnatfac
      do n=1,nparmx
        write(82,'(100(1pg24.16))') parmnat(1:nparmx,n)
      enddo
      do k=1,ndats
        if(maxdat(k)==0) CYCLE
        write(82,'(a,2i4,2x,a)')
     *       'data block: nr,maxdata,name=',k,maxdat(k),heatyp(ifit1(k))
        do i=1,npts
          if(lifit(k)<12 .and. varobs(lifit(k),i,0).NE.0D0) then
             write(82,'(1pg24.16)') varobs(lifit(k),i,0)
             write(82,'(100(1pg24.16))') 
     &         (varobs(lifit(k),i,n)-varobs(lifit(k),i,0),n=1,nparmx)
          endif
        enddo
      enddo
      close(82)
c
      return
      end
