c=====================================================================
c
c        program shell-gigan
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c      shell-gigan calculates the strutinsky smoothed  energy.      c
c            the single particle energies belonging to H and H_0    c
c            are read in from file swbinp                           c
c                                                                   c      
c                                                                   c
c                                                                   c
c      written by w. nazarewicz, 1985, 1994                         c
c  modified by A. T. Kruppa and T. Vertse 1999 in order to include  c
c  continuum effect  and generalized smoothing procedure            c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit double precision (a-h,o-z)
       parameter(maxlev=1000,ndmain=500,nrmax=80)
       character*72 swbinp,shellout,shellres,shellrdens
       character*20 title,titdum,titsav
       dimension cg(2000),tab(20),dtab(20),sig(ndmain),ilf(20),
     *  esm(ndmain),alam(ndmain)
       dimension enipr(25000),enine(25000),eni0pr(25000),eni0ne(25000)
     * ,ps(25000),ps0(25000),rdensp(0:nrmax,25000),
     & rdensn(0:nrmax,25000),rde(0:nrmax,25000),
     & rdensp0(0:nrmax,25000),rdensn0(0:nrmax,25000),rx(0:nrmax),
     &rhop(0:nrmax),
     & rhon(0:nrmax),rhop0(0:nrmax),rhon0(0:nrmax)
       dimension imlp(25000),imlp0(25000),imln(25000),imln0(25000)
      common
     *       /strut/ her(40),dher(40),cm(40),gammaq,
     *               ip,ilev,inucl
       common/bcs/ eni(25000)/bcs0/eni0(25000)/old/iold
     :       /ome/ aa
       common/fura/rdensp,rdensn,rdensp0,rdensn0,rx
       data nlevm/25000/
       data nf0,nf1,nf2/10,11,12/
c
c
c         i n p u t      i n p u t      i n p u t
c
c         swbinp    -   file which stores sp energies
c
c         shellres  -   file containing results of the shellmesh code
c
c         shellout  -   output file
c
c
c         iz        -   proton number
c
c         in        -   neutron number
c
c         iprot     -   0  (neutron) 1 (proton)
c
c
c         ecut      -   cut-off energy above lambda (in hbar*omega)
c
c         iendd     -   0 (continue), 1 (end)
c
c
c         ipart     -   1  calculate smoothed particle number
c         idens     -   1  calculate smoothed densities
c
       namelist /strutek/ swbinp,shellout,shellres,shellrdens,iz,in,
     *                    gam0,gam1,dgam,ip0,ip1,ipart,idens,alg,felg,
     *                    ecut,enei,enef,dene,ilev,iold,iendd

       enei=-80.d0
       enef=30.d0
       dene=0.1d0
       ecut=60.2
       ilev=0
       nstep=70
       iold=0
       idens=0
       iendd=0
       title = '                    '
       open(unit=5,file='levdenew.inp',status='old')  
c
       ifirst=1
 5     continue
       read(5,strutek)
       write(*,strutek)
c
       if (iendd.eq.1) stop 'allright'
c
       if (ifirst.eq.1) then
         open(nf2,file=shellres,form='formatted',status='unknown')
         open(16,file=shellout,form='formatted',status='unknown')
	 if(idens.eq.1) open(9,file=shellrdens,form='formatted',
     & status='unknown')
     
         write(nf2,'(a,i4,a,i4)')
     &  ' shell corrections for RMF:'
         ifirst=0
       endif
c
      do i=72,1,-1
        lentit = i
        if(swbinp(i:i).ne.' ') goto 998
      enddo
 998   continue
c
       title = swbinp(1:20)
       open(nf1,file=title,form='formatted',status='unknown')
       call swbrea(enipr,imlp,ngorp,nf1,title)  ! read protons
c       write(*,997)(imlp(i),enipr(i),i=1,ngorp)
997    format(1x,i5,2x,f10.3)
       call swbrea(enine,imln,ngorn,nf1,title)  ! read neutrons
c       write(nf2,'(a)')title
       call swbrea(eni0pr,imlp0,ngor0p,nf1,titdum)  ! read proton gas 
       call swbrea(eni0ne,imln0,ngor0n,nf1,titdum)  ! read neutron gas 
       if(idens.eq.1)then

       call swbdens(rde,rx,nrmax,nstep,nsph,9,imlp)
       write(*,*)'nsph=',nsph
       call expand(rde,nrmax,nstep,nsph,imlp,rdensp,nrdp)
       write(*,*)'nrdp=',nrdp
c here rde still hold the spherical proton density
       iiiz=0
       do i=1,nrdp
         iorb=imlp(i)
         iorbn=imlp(i+1)
         iiiz=iiiz+iorb
         if(iz.ge.iiiz.and.iz.lt.iiiz+iorbn)then
           write(*,*)' i,iiiz',i,iiiz
           ifermi=i
           ival=iz-iiiz
         endif
       enddo
       write(*,*)' ifermi,ival=',ifermi,ival
c  calculation of the non-smoothed density for protons
c
       do i=0,nstep
       rhop(i)=0.d0
	do n=1,ifermi
          rhop(i)=rhop(i)+rde(i,n)*imlp(n)	
	enddo 
	if(ival.ne.0)then
	  ifermi1=ifermi+1
          rhop(i)=rhop(i)+rde(i,ifermi1)*ival
        endif
       enddo       



       call swbdens(rde,rx,nrmax,nstep,nsph,9,imln)

       write(*,*)'nsph=',nsph
       call expand(rde,nrmax,nstep,nsph,imln,rdensn,nrdn)
       write(*,*)'nrdn=',nrdn
       call swbdens(rde,rx,nrmax,nstep,nsph,9,imlp0)
       write(*,*)'nsph=',nsph
       call expand(rde,nrmax,nstep,nsph,imlp0,rdensp0,nrdp0)
       write(*,*)'nrdp0=',nrdp0
       call swbdens(rde,rx,nrmax,nstep,nsph,9,imln0)
       write(*,*)'nsph=',nsph
       call expand(rde,nrmax,nstep,nsph,imln0,rdensn0,nrdn0)
       write(*,*)'nrdn0=',nrdn0
       close(unit=9)       
       endif
       open(unit=9, file='density.plt',status='unknown')
       write(9,199)(rx(i),rhop(i),i=0,nstep)
199    format(1x,f10.3,2x,e12.4)       
       close(9)
       if(idens.eq.1)then
       open(unit=9, file='smoothdensity.plt',status='unknown')
       endif
c
c       write(*,*)' iz,in',iz,in
       write(nf2,*)'iz in',iz,in
       aa=iz+in
         write(nf2,1230)
c 1230   format('  N_Z       ',8x,
c     1 '  ip   gam    gamma    Eshell   Esmooth    xlamb   Esp',/
c     2,20x,'      [hom0]  [MeV]    [MeV]     [MeV]     [MeV]   [MeV]')   
 1230   format(' t_z  ',
     1 '  ip   gam    gamma    Eshell   Esmooth    xlamb   Esp',/,6x
     2,'      [hom0]  [MeV]    [MeV]     [MeV]     [MeV]   [MeV]')   
c
       h=41.d0/aa**0.33333333333d0
c
       if(ilev.gt.0)open(unit=2,file='levden.plt',status='unknown')
100    format(1x,8a4)
c  loop for iprot is commented!!
c       do iprot=1,0,-1
c
c             strutinski smooth energy
c
       if(iprot.eq.0) then
c         gam   = gamneu
         ngor  = ngorn
c         ip    = ippro
         do i=1,ngor
           eni(i)  = enine(i)
           eni0(i) = eni0ne(i)
         enddo
       else
c         gam   = gampro
         ngor  = ngorp
c         ip    = ipneu
         do i=1,ngor
           eni(i)  = enipr(i)
           eni0(i) = eni0pr(i)
         enddo
       endif
c
       ngam=(gam1-gam0)/dgam+1
       npp=1
       np=(ip1-ip0)/2+1
       ng0=1
       sala=0.d0
c  loop for p values
       do ip=ip0,ip1,2
       ip2=ip/2
c
c     gamm=gam
c      
c   starting point of the best fit in gamma
       idec=-1
       sig0=1.d+10
c  loop for gamma (gamm is in units h here)
       do ng=ng0,ngam
       gamm=gam0+(ng-1)*dgam      
c
c        write(6,'(a,i3)') 'call SHCORR: iprot=',iprot
        call shcorr(iz,in,iprot,esmoo,ngor,esngp,xlamb,ip,ecut,gamm)
c
       if(ilev.gt.0) then
         nen=(enef-enei)/dene+1
         if(nen.gt.2000)then
           write(*,130)nen
 130       format('    nen=',i6,'  is > 2000')
           stop
         endif
       endif
       gammq=gamm*h
c*****************
c  define energy range in which the best linear fit is required [alg,felg]
c
        xhalf=eni(1)*0.5d0
c        felg=xhalf+0.5D0*h
c        felg=xhalf+fel*h
c        alg=felg-sav*h
c        felg=-20.d0
c	alg=-50.d0
        iel=(alg-enei)/dene+1
        ief=(felg-enei)/dene+1
c  fill array cg with smoothed level density in the energy range [alg,felg]
c
c       
       i=0
c       do ie=1,nen
       do ie=iel,ief
       i=i+1
       ene=enei+(ie-1)*dene 
       zz=hg(ene)
       if (iold.eq.0) then
          z0=hg0(ene)
          cg(i)=zz-z0
c          write(2,105)ene,zz-z0
       else
         cg(i)=zz
c          write(2,105)ene,zz
       endif
105    format(1x,f12.4,3e12.4)
c
c      end of the loop for the energy
c
       enddo
c       write(2,*)
c**********************       
        gm=gammq
c********
          oldslope=slope
c  finding the best straight line fit to the smoothed level density and the
c  khi square of the line: sigma
        call lfif(2000,cg,i,sigma,slope,dene)
c        write(*,*)' sigma,slope=',sigma,slope
        sig(ng)=sigma
        if(sigma.lt.sig0)then
          sig0=sigma
          if(ng.gt.ng0)idec=1
        else
        if(idec.gt.0)then
c  minimum of sigma(gamma) has been found. Restore the position of the minimum
        gamm=gam0+(ng-2)*dgam
c  calculates smoothed energy at the minimum
        call shcorr(iz,in,iprot,esmoo,ngor,esngp,xlamb,ip,ecut,gamm)
        gammq=gamm*h
        ala=xlamb
        sala=sala+ala
        gm=gammq
c  shell correction at the minimum
        eshell=esngp-esmoo
        esm(ng)=eshell
        alam(ng)=xlamb
c********
c        write(nf2,115)gm,eshell,xlamb,oldslope
c115     format(1x,f8.4,f12.4,f10.4,f10.4)
c        write(6,114)ip,gm,eshell,esmoo,ala,oldslope
c114     format(1x,i3,f8.4,2f10.3,2f10.3,f10.3,f9.3)
c*********
        if(iprot.eq.0)then
c	write(nf2,'(a)')title
c	write(*,'(a)')title
c        write(nf2,124)title,ip,gamm,gm,eshell,esmoo,xlamb,esngp
        write(nf2,124)ip,gamm,gm,eshell,esmoo,xlamb,esngp
        else
c	write(nf2,'(a)')title
	write(*,'(a)')title
c        write(nf2,114)title,ip,gamm,gm,eshell,esmoo,xlamb,esngp
        write(nf2,114)ip,gamm,gm,eshell,esmoo,xlamb,esngp
        endif
115     format(1x,f8.4,f12.4,f10.4,f10.4)
c124     format(1x,a,'neutron ',i3,f7.2,f8.4,2f10.3,2f10.3,
c     1 f10.3,f9.3,f10.3)
c114     format(1x,a,'proton  ',i3,f7.2,f8.4,2f10.3,2f10.3,
c     1 f10.3,f9.3,f10.3)
124     format(1x,'neu ',i5,f8.3,f8.4,2f10.3,2f10.3,
     1 f10.3,f9.3,f10.3)
114     format(1x,'pro ',i5,f8.3,f8.4,2f10.3,2f10.3,
     1 f10.3,f9.3,f10.3)

c********
      if(ilev.gt.0)then
c  calculation of the smoothed level density at the minimum

c MB print smoothed level density
c
        do 9876 ener=-80.0,20.0,0.2 
          hggg=hg(ener)
          hggg0=hg0(ener)
          write(16,'(4f12.6)') ener, hggg, hggg0, hggg-hggg0 
          write(2,'(4f12.6)') ener,  hggg-hggg0 
9876    continue
        write(2,*)
        if(ipart.eq.1)then
c  calculate smoothed particle numbers for each orbits	
	  call partsm
     &    (xlamb,eni,eni0,her,dher,ip,cm,ngor,inucl,gm,ps,ps0)
          sumnpart=0.d0
	  do i=1,ngor
            sumnpart=sumnpart+ps(i)-ps0(i)
          enddo
          write(*,*)' sumnpart=',sumnpart
        endif
	
      endif
       
       goto 123
       
       endif
       
       endif
       aa=iz+in
       part=in
       if (iprot.eq.1) part=iz
c       
c       pairing contribution is left out for the moment
c       call hhepcba(xlamb,in,iz,iprot,robar,gpair)
c
       eshell=esngp-esmoo
       esm(ng)=eshell
       alam(ng)=xlamb
       x=eshell
       write(3,'(5f11.5)')gammq,x
       write(16,113) ip,gamm,esmoo,esngp,eshell,xlamb
113    format(1x,i4,f10.5,2f12.5,f10.5,f10.5)       
c
       enddo
c  end of the loop for gamma       
c
       if(ilev.gt.0) then
        write(*,*)'  no minimum found in khi2'
        goto 55
       endif
c  storing minimum values of the smooth energy in the given p
123    tab(npp)=esmoo
       npp=npp+1
       ng0=ng-1
112    format(1x,i3,f8.4,5f10.3,f9.3)
       enddo
c     end loop for p
      if (ilev.gt.0)close(2)
      npp=npp-1
c  calculation of the average values of the smoothed  and Fermi energies 
c   and their variations for the different p values
      call avrg(tab,dtab,npp,sala,esngp,esmooav,var)
 55   continue
c loop for iprot is commented!!
c      enddo
c  end of the loop for iprot
c  go to the starting point      
      goto 5

  2   format(1x,'z=',i3,'  n=',i3,4x,a)
  3   format(1x,'z0=',i3,'  n0=',i3,/,
     :       1x,'beta2=',f7.4,2x,
     :       '  beta3=',f7.4,2x,'beta4=',f7.4,2x,'beta5=',f7.4,2x,
     :       'beta6=',f7.4,2x,/,1x,'nmax=',i2,' p=',i2,' gamma=',f6.2,
     :       '  ecut=',f7.2,/,1x,'esmooth=',f10.3,'  esp=',f10.3,
     :       '  eshell=',f10.4,'  lambda=',f10.4,/,1x,
     :       'robar/a=',f10.5,'  gpair=',f10.5,/,1x,75(1h*),/)

      end
c
c
      subroutine partsm
     &(e,eni,eni0,her,dher,ipher,cm,ilev,inucl,gammaq,ps,ps0)
c  calculates smoothed particle number for each orbits (i) for the
c  HF : ps(i) and for the gas: ps0(i) separately
      implicit double precision (a-h,o-z)
c
      common /old/iold
      dimension
     :          eni(*),eni0(*),her(*),dher(*),cm(*),ps(*),ps0(*)
      data
     :     pi /3.1415926536d0/
c
c         xnn calculates f and df. called from xn.
c
c         eni       -  single particle spectrum
c         her,dher  -  hermite polynomials and their derivatives
c         ipher     -  the order of smoothing in strutinski method
c         ilev      -  the number of proton or  neutron states summed
c         cm        -  auxiliary coefficients (calculated in desc)
c         ps        -  smoothed particle number (HF)
c         ps0       -  smoothed particle number (free)
      sqpi=sqrt(pi)

      do n=1,ilev
      un=(e-eni(n))/gammaq
      un2=un*un
      expun=0.d0
c
      if (un2.gt.30.0d0) go to 10
c
      expun=exp(-un2)
c
      call dhep(un,ipher,her,dher)
c
      sump=0.0d0
c
      do m=2,ipher,2
      sump=sump+cm(m)* her(m)
      enddo
c
   10 continue
c
      skl =+0.5d0*(1.d0+errf(un))-expun*sump/sqpi
      ps(n) =2.d0*skl
c
      ps0(n)=0.d0
      enddo
c
      if (iold.eq.1) return

c***************************free       (gas)

      do n=1,ilev
c
      un=(e-eni0(n))/gammaq
      un2=un*un
      expun=0.d0
c
      if (un2.gt.30.0d0) go to 100
c
      expun=exp(-un2)
c
      call dhep(un,ipher,her,dher)
c
      sump=0.0d0
c
      do m=2,ipher,2
      sump=sump+cm(m)* her(m)
      enddo
c
  100 continue
c
      skl =+0.5d0*(1.d0+errf(un))-expun*sump/sqpi
      ps0(n)=2.d0*skl
c
      enddo
c
c
      return
      end
c

*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      subroutine swbdens(rdens,rx,nrmax,nstep,ilelev,nr,iml)
*
*-----------------------------------------------------------------------
*
*     this routine reads energies produced by swbgigan (mamut version).
*     this is all in double precision. one set of data is read.
*        input  :  nr     - unit number of swbgigan result
*                           (make sure ener is properly dimensioned
*                            in the calling program)
*        output :  rdens   - radial densities 
*                  ilelev - number of levels read. ilelev=0 means that
*                           end of file has been encountered
*                  + characeristics of the spectrum, see header
*
*-----------------------------------------------------------------------
*--------------------------------------------- t.r.werner, february 1994
*-----------------------------------------------------------------------
*
      implicit double precision (a-h,o-z)
c     integer    dom,domm
      dimension  rdens(0:nrmax,*),rx(0:nrmax),iml(*)
      character*6 aline
*
      ilelev=0
      ind=1
c       write(*,*)' nr=',nr
9      do iread=1,9999
        read(nr,102,err=2,end=94) aline,i1,i2,e
	if(aline(2:6).eq.'state')then
	n=i1
	j=i2
        write(*,*)n,j,e
	endif
2        if(aline(1:4).eq.'  r ') goto 999
        enddo


102   format(a6,14x,i2,3x,i2,18x,f7.2)
      write(*,*)' r is not found'
      goto 98

 999  continue
       write(*,9123)aline(1:6)
123    format(f6.2,e13.5)
9123	format((a))

        if(n.eq.1.and.ind.gt.1)then
	 goto 19
	endif

       read(nr,123)(rx(i),rdens(i,ind),i=0,nstep)
c       write(*,123)(rx(i),rdens(i,ind),i=0,nstep)
       iml(ind)=j+1
c       do i=0,nstep
c       rdens(i,ind)=rdens(i,ind)*(j+1.d0)
c       enddo
       ind=ind+1
       goto 9	
 19    continue
       ilelev=ind-1
       write(*,*)'ind ',ilelev
*
c
*
   1   return
  94   ilelev=ind-1
       write(*,*)' end of file reached'
       return
 100   format(i4,i5,f9.3)
 101   format(1x,i4,i5,f9.3)
*
 98    continue
       stop ' EoF encountered before densies read in'
       end
       
      subroutine expand(rdens,nrmax,nstep,ilelev,iml,rdensnew,iveg)
*
*-----------------------------------------------------------------------
*
*     this routine  expands the degeneracy of the spherical symmetric
c      case to the double degeneracy of the deformed case
c     arrays are filled until 1:ilelev
c     iml(i)=2j+1 rdens(ir,i) is the radial density of the i-th
c     spherically symmetric orbit at r=ir*0.3
c     rdensnew(ir,i) stores the multiplicated densities from
c     i=1,..,iveg
c
      implicit double precision (a-h,o-z)
      dimension  rdens(0:nrmax,*),iml(*),rdensnew(0:nrmax,*)
      ind=1

      do isph=1,ilelev
       im2=iml(isph)/2
       im1=im2-1
c  multiply the densities im2 times
       do j=0,im1
         do k=0,nstep
          rdensnew(k,ind+j)=rdens(k,isph)
         enddo
       enddo
c         write(*,100)rdens(nstep,isph)
c         write(*,101)(rdensnew(nstep,ind+j),j=0,im1)
         ind=ind+im2
      enddo
      iveg=ind
c      write(*,*)'iveg ',iveg
      
      return
 100   format(1x,e12.4)
 101   format(1x,4(e12.4))
*
 98    continue
       stop ' EoF encountered before densies read in'
      end
      
      subroutine swbrea(ener,iml,ilelev,nr,title)
*
*-----------------------------------------------------------------------
*
*     this routine reads energies produced by swbgigan (mamut version).
*     this is all in double precision. one set of data is read.
*        input  :  nr     - unit number of swbgigan result
*                           (make sure ener is properly dimensioned
*                            in the calling program)
*        output :  ener   - ordered vector of energies
*                  ilelev - number of levels read. ilelev=0 means that
*                           end of file has been encountered
*                  + characeristics of the spectrum, see header
*                  title  - marker in the inptu file to transfer a comment
*
*-----------------------------------------------------------------------
*--------------------------------------------- t.r.werner, february 1994
*-----------------------------------------------------------------------
*
      implicit double precision (a-h,o-z)
c     integer    dom,domm
      dimension  ener(*),iml(*)
      character*4 aline
      character*14 title
*
      ilelev=0
*
      do iread=1,9999
        read(nr,'(a)',err=98,end=98) aline
c       write(*,9123)aline(1:4)
9123	format((a))
c       write(*,*)' eddig eljott a ciklusban'

        if(aline(1:4).eq.'  2j') goto 999
        if(aline(1:4).eq.' * *') then
          read(nr,'(a)',err=98,end=98) title
	  write(*,'(a)') title
        endif
      enddo
 999  continue
*
      ind=0
cold 9    read(nr,101)j2,l,energ
c 9    read(nr,123)j2,l,energ,xn,imult
 9    read(nr,*)j2,l,energ,xn,rms,imult
123   format(2i5,f9.4,e17.5,i5)
c      write(*,*)j2,l,energ,xn,rms,imult
      ind=ind+1
      ener(ind)=energ
      iml(ind)=imult
      if (j2.lt.0) goto 19
      goto 9
 19   continue
      ilelev=ind
      write(*,*)'ind ',ind
*
      call ord(ener,iml,ilelev)
c
*
    1 return
 100   format(i4,i5,f9.3)
 101   format(1x,i4,i5,f9.3)
*
 98    continue
       stop ' EoF encountered before spectral data'
      end
*
*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      subroutine ord(ra,iml,n)
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*-----------                                                 -----------
*----------- the routine rearranges vector ra of dimension n -----------
*-----------   so that it is in ascending order on return.   -----------
*-----------                                                 -----------
*-----------------------------------------------------------------------
*-----------                                                 -----------
*-----------       method used: heapsort of williams         -----------
*-----------                                                 -----------
*-----------------------------------------------------------------------
*--------------------------------------------- t.r.werner, february 1991
*--- from "numerical  recipes" -----------------------------------------
*-----------------------------------------------------------------------
*
      implicit double precision (a-h,o-z)
      dimension ra(*),iml(*)
*
      if (n.lt.1) stop 'heap00'
      if (n.eq.1) return
*
      l=n/2+1
      ir=n
*
    1 continue
*
         if (l.gt.1) then
            l=l-1
            rra=ra(l)
            irra=iml(l)
         else
            rra=ra(ir)
            ra(ir)=ra(1)
            irra=iml(ir)
            iml(ir)=iml(1)
            ir=ir-1
            if (ir.eq.1) then
               ra(1)=rra
               iml(1)=irra
               return
            end if
         end if
         i=l
         j=l+l
*
    2    continue
*
         if (j.le.ir) then
            if (j.lt.ir.and.ra(j).lt.ra(j+1)) j=j+1
            if (rra.lt.ra(j)) then
               ra(i)=ra(j)
               iml(i)=iml(j)
               i=j
               j=j+j
            else
               j=ir+1
            end if
            go to 2
         end if
         ra(i)=rra
         iml(i)=irra
*
      go to 1
*
      end
c
c--------------------------------------------------------------------
c
      subroutine shcorr(iz,in,iprot,esmoo,infit,esp,xlmbar,ip,ahom,gam)
c
      implicit double precision (a-h,o-z)
c
      common
     *       /strut/ her(40),dher(40),cm(40),gammaq,
     *               ipher,ilev,inucl
       common/bcs/ eni(25000)/bcs0/eni0(25000)/old/iold
     *       /ome/   aa
c
c          shcorr calculates the shell correction to the liquid drop
c          model formula
c
c          h is the oscillator freqeuncy
c
c      write(6,'(a,i3)') ' enter SHCORR: iprot=',iprot
      h=41.d0/aa**(1.d0/3.d0)
      gammaq=gam*h
      ipher=ip
c
      inucl=iz
      if (iprot.eq.0) inucl=in
c
      inupol=inucl/2
c
c       all levels up to (e_fermi + ahom*h) are considered
c
      if (iold.eq.1) then
      porr=eni(inupol)+ahom*h
c
      do i=inupol,infit
      ilev=i
      if (eni(i).gt.porr) go to 41
      enddo
      else
c
       ilev=infit
       endif
   41 continue
c
      call desc(iz,in,xlmbar,esmoo,esp)
c
c      write(6,'(a,i3)') ' leave SHCORR: iprot=',iprot
      return
      end
c
c-------------------------------------------------------------------
c
      subroutine desc(iz,in,xlmbar,esmoo,esp)
c
      implicit double precision (a-h,o-z)
c
c                 ****  called from shcorr  ****
c
      common
     *       /strut/ her(40),dher(40),cm(40),gammaq,
     *               ip,ilev,inucl
       common/bcs/ eni(25000)/bcs0/eni0(25000)/old/iold
      external
     :         xn
      data
     :     pi,iter,eps /3.1415926536d0,200,1.0d-7/

c
c         desc calculates strutinsky shell correction
c
c         gammaq   -  smoothing constant (in mev, bentzwale)
c         ip       -  smoothing order
c         ilev     -  the number of states in the summations
c         xlmbar   -  resulting fermi energy
c
      sqpi=sqrt(pi)
      sil=1.
c
      do k=2,ip,2
      sil=sil*k/2.
      cm(k)=(-1)**(k/2)/2.0**k/sil
      enddo
c
      ipher=ip
      inupol=inucl/2
      first=0.5*(eni(inupol)+eni(inupol+1))
c
      call rtni(xlmbar,f,df,xn,first,eps,iter,ier)
c
      if (ier.ne.0) stop 'desc/rtni'
c
      sum =0.0
      esp =0.0
      inu =inupol
c
      do n=1,ilev
c
      un=(xlmbar-eni(n))/gammaq
c
      call dhep(un,ipher,her,dher)
c
      expun=0.
      un2=un*un
c
      if (un2.gt.30.0d0) go to 10
c
      expun=exp(-un2)
      sump=0.0
      socc=0.0
c
      do m=2,ipher,2
      sump=sump+cm(m)*(0.5*gammaq*her(m+1)+eni(n)*her(m)
     *         +m*gammaq*her(m-1))
      socc=socc+cm(m)*her(m)
      enddo
c
   10 continue
c
      errofn=errf(un)
      skl =   +0.5*eni(n)*(1.+errofn)-0.5*gammaq/sqpi*expun
     *        -expun*sump/sqpi
      sum=sum+skl+skl
c
      if (n.le.inupol) esp=esp+eni(n)*2.d0
      enddo
c
      esmoo=sum
c
c          correction for the odd particle number
c
      if (mod(inucl,2).ne.0) esp=esp+eni(inupol+1)
      if (iold.eq.1) return
c
      sum =0.0d0
c      esp =0.0
c      inu =inupol
c
      do n=1,ilev
c
      un=(xlmbar-eni0(n))/gammaq
c
      call dhep(un,ipher,her,dher)
c
      expun=0.
      un2=un*un
c
      if (un2.gt.30.0d0) go to 100
c
      expun=exp(-un2)
      sump=0.0d0
      socc=0.0d0
c
      do m=2,ipher,2
      sump=sump+cm(m)*(0.5*gammaq*her(m+1)+eni0(n)*her(m)
     *         +m*gammaq*her(m-1))
      socc=socc+cm(m)*her(m)
      enddo
c
100   continue
c
      errofn=errf(un)
      skl =   +0.5d0*eni0(n)*(1.d0+errofn)-0.5d0*gammaq/sqpi
     & *expun-expun*sump/sqpi
      sum=sum+skl+skl
c
c      if (n.le.inupol) esp=esp+eni(n)*2.d0
      enddo
c
      esmoo=esmoo-sum
      return
      end
c
c-------------------------------------------------------------------
c
      subroutine xn(ue,f,df)
c
      implicit double precision (a-h,o-z)
c
      common
     *       /strut/ her(40),dher(40),cm(40),gammaq,
     *               ipher,ilev,inucl
       common/bcs/ eni(25000)/bcs0/eni0(25000)/old/iold
c
c         called from rtni (for calculating lambda bar)
c
      call xnn(ue,eni,eni0,her,dher,ipher,cm,ilev,inucl,gammaq,f,df)
c
      return
      end
c
c--------------------------------------------------------------------
c
      subroutine xnn
     &(e,eni,eni0,her,dher,ipher,cm,ilev,inucl,gammaq,f,df)
c
      implicit double precision (a-h,o-z)
c
      common /old/iold
      dimension
     :          eni(*),eni0(*),her(*),dher(*),cm(*)
c     :          eni(1),eni0(1),her(1),dher(1),cm(1)
      data
     :     pi /3.1415926536d0/
c
c         xnn calculates f and df. called from xn.
c
c         eni       -  single particle spectrum
c         her,dher  -  hermite polynomials and their derivatives
c         ipher     -  the order of smoothing in strutinski method
c         ilev      -  the number of proton or  neutron states summed
c         cm        -  auxiliary coefficients (calculated in desc)
c
      sqpi=sqrt(pi)
      sum=0.0d0
      dsum=0.0d0
c
      do n=1,ilev
c
      un=(e-eni(n))/gammaq
      un2=un*un
      expun=0.d0
c
      if (un2.gt.30.0d0) go to 10
c
      expun=exp(-un2)
c
      call dhep(un,ipher,her,dher)
c
      sumd=0.0d0
      sump=0.0d0
c
      do m=2,ipher,2
      sump=sump+cm(m)* her(m)
      sumd=sumd+cm(m)*dher(m)
      enddo
c
   10 continue
c
      skl =+0.5d0*(1.d0+errf(un))-expun*sump/sqpi
      dskl=expun*(1.d0+2.d0*un*sump-sumd)
      sum =sum + skl +skl
      dsum=dsum+dskl+dskl
c
      enddo
c
      f=inucl-sum
      df=-dsum/gammaq/sqpi

      if (iold.eq.1) return

c***************************free       

      sum=0.0d0
      dsum=0.0d0
c
      do n=1,ilev
c
      un=(e-eni0(n))/gammaq
      un2=un*un
      expun=0.d0
c
      if (un2.gt.30.0d0) go to 100
c
      expun=exp(-un2)
c
      call dhep(un,ipher,her,dher)
c
      sumd=0.0d0
      sump=0.0d0
c
      do m=2,ipher,2
      sump=sump+cm(m)* her(m)
      sumd=sumd+cm(m)*dher(m)
      enddo
c
  100 continue
c
      skl =+0.5d0*(1.d0+errf(un))-expun*sump/sqpi
      dskl=expun*(1.d0+2.d0*un*sump-sumd)
      sum =sum + skl +skl
      dsum=dsum+dskl+dskl
c
      enddo
c
      f=f+sum
      df=df+dsum/gammaq/sqpi
c
      return
      end
c
c--------------------------------------------------------------------
c
      subroutine dhep(x,n,her,dher)
      implicit double precision (a-h,o-z)
      dimension
     :          her(*),dher(*)
c     :          her(1),dher(1)
*
*         dhep calculates hermite polynomials at point x up to n-th
*         degree. it also calculates  derivatives  dher(i) of these
*         polynomials.
*                            dher(i)=2*i*her(i-1)
*
      her (1)=1.d0
      dher(1)=0.d0
*
      if (n) 1,1,2
*
   1  return
*
   2  her(2)=x+x
      dher(2)=2.d0
*
      if (n-1) 1,1,3
*
   3  continue
*
      do 4 i=2,n
*
      f=x*her(i)-(i-1)*her(i-1)
      her(i+1)=f+f
      df=i*her(i)
      dher(i+1)=df+df
*
   4  continue
*
      end
c
c-------------------------------------------------------------------
c
      function errf(y)
c
      implicit double precision (a-h,o-z)
c
      data
     *     a1,a2,a3,a4,a5,a6
     *    /.0705230784d0,.0422820123d0,.0092705272d0,
     *     .0001520143d0,.0002765672d0,.0000430638d0/
c
      x=abs(y)
      zr=(x*(a1+x*(a2+x*(a3+x*(a4+x*(a5+x*a6)))))+1.)
c
      if (abs(zr).lt.100.d0) go to 1
c
      zrf=1.0d0
c
      go to 2
c
   1  zrf=1.0d0-1.0d0/zr**16
c
   2  continue
c
      errf=sign(zrf,y)
c
      return
      end

c=======================================================================
      subroutine hhepcba(xlambr,in,iz,iprot,robar,gpair)
      implicit double precision (a-h,o-z)
c
c         hepcba calculates the smoothed value of the .hepc.
c
c         xlambr   -  strutinskys lambda bar
c         gpair    -  pairing force strength being calculated
c
      aa=iz+in
c
      if (iprot.eq.0) delavr=8.6d0/sqrt(aa)
      if (iprot.eq.1) delavr=9.9d0/sqrt(aa)
c
      robar=hg(xlambr)
c
      if (iprot.eq.1) then
         ntot=iz
         nsign=1
      else
         ntot=in
         nsign=-1
      end if
c
      iodd=0
      if (mod(ntot,2).eq.1) iodd=1
c
      if (iodd.eq.1) ntot=ntot-1
c
      nnd=1
      nng=ntot
c
      xntot=ntot
c
      y1=(-xntot/2.d0+nnd-1)/robar
      y2=(-xntot/2.d0+nng  )/robar
c
      gpair=2.d0/robar/(log(sqrt(y2**2+delavr**2)+y2)-
     *                log(sqrt(y1**2+delavr**2)+y1))
c
      return
      end
c=======================================================================
      function hg(e)
      implicit double precision (a-h,o-z)
      common
     *       /strut/ her(40),dher(40),cm(40),gammaq,
     *               ipher,ilev,inucl
     *       /bcs/   eni(25000)/bcs0/eni0(25000)
      data
     :     pi /3.1415926536d0/
c
c         hg calculates strutinsky smoothing function
c
      sum =0.0d0
c
      do 1 n=1,ilev
c
      un=(e-eni(n))/gammaq
      expun=0.d0

      un2=un*un
c
      if (un2.gt.30.0d0) go to 1
c
      expun=exp(-un2)
      sump=1.0d0
      call dhep(un,ipher,her,dher)
c
      do  i=2,ipher,2
      sump=sump+cm(i)*her(i+1)
      enddo
c
      sum =sum +expun*sump
 1    continue
c
      hg  =sum /gammaq/sqrt(pi)
      hg=hg+hg
c
      return
      end
c
c=======================================================================
c
      function hg0(e)
      implicit double precision (a-h,o-z)
      common
     *       /strut/ her(40),dher(40),cm(40),gammaq,
     *               ipher,ilev,inucl
     *       /bcs/   eni(25000)/bcs0/eni0(25000)
      data
     :     pi /3.1415926536d0/
c
c         hg calculates strutinsky smoothing function
c
      sum =0.0d0
c
      do 1 n=1,ilev
c
      un=(e-eni0(n))/gammaq
      expun=0.d0

      un2=un*un
c
      if (un2.gt.30.0d0) go to 1
c
      expun=exp(-un2)
      sump=1.0d0
      call dhep(un,ipher,her,dher)
c
      do  i=2,ipher,2
      sump=sump+cm(i)*her(i+1)
      enddo
c
      sum =sum +expun*sump
 1    continue
c
      hg0  =sum /gammaq/sqrt(pi)
      hg0=hg0+hg0
c
      return
      end


       subroutine lfif(ndm,y,n,sig,slope,dene)
       implicit real*8 (a-h,o-z)
       dimension y(ndm),xk(2000),c(0:5)
      do i=1,n
      xk(i)=i*dene
      enddo
      nagyn=1
      call normal(xk,y,nagyn,n,c)
      slope=c(1)
cd      write(6,102)(i,c(i),i=0,1)
102   format('  i',10x,'c(i)',/(i5,f15.4))
      sig=delt(xk,y,nagyn,n,c)
cd      write(6,103) sig
103   format('   delta=',e12.5)
       return
       end
      subroutine avrg(tab,dtab,npp,sala,esngp,esmooav,var)
      implicit real*8 (a-h,o-z)
      dimension dtab(npp),tab(npp)
c      complex*16 eshell,esngp
      sala=sala/npp
      esmooav=0.d0
      do i=1,npp
      esmooav=esmooav+tab(i)
      dt=dtab(i)
      enddo
      esmooav=esmooav/npp
      var=0.d0
      do i=1,npp
      dv=(tab(i)-esmooav)**2
      var=var+dv
      enddo
      var=sqrt(var/npp)
      write(12,125) esmooav,var
      eshell=esngp-esmooav
      write(12,126)eshell,var,sala
126   format(' average shell correction=',f10.3,' +- ',f8.3,/
     1' average Fermi level',f10.3)
125   format(' average value of esmooth=',f10.3,' var=',f7.3)
      return
      end
