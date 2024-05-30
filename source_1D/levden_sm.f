c====================================================================
c====================================================================
c        program shell-gigan
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c      shell-gigan calculates the strutinsky smoothed  energy.      c
c            the single particle energies belonging to H and H_0    c
c            are read in from file swbinp                           c
c                                                                   c      
c                                                                   c
c      written by w. nazarewicz, 1985, 1994                         c
c  modified by A. T. Kruppa and T. Vertse 1998 in order to include  c
c  continuum effect                                                 c                                                            c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit double precision (a-h,o-z)
       parameter(maxlev=1000,ndmain=500)
       character*72 swbinp,shellout,shellres,shellsmooth,shellrho
       character*20 title,titdum,titsav
       dimension cg(2000),tab(20),dtab(20),sig(ndmain),ilf(20),
     *           esm(ndmain),alam(ndmain)
       dimension enipr(25000),enine(25000),eni0pr(25000),eni0ne(25000)
       dimension j2pr(25000),j2ne(25000),j2dum(25000)
       dimension ps(25000),ps0(25000)
    
       common/strut/ her(40),dher(40),cm(40),gammaq,
     *               ip,ilev,inucl
       common/bcs/ eni(25000),j2act(25000)
     &       /bcs0/eni0(25000)
     &       /old/iold
     :       /ome/ aa

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
c         gampro    -   smoothing constant for protons (in hbar*omega)
c         gamneu    -   smoothing constant for neutrons (in hbar*omega)
c
c         ippro     -   smoothing order for protons
c         ipneu     -   smoothing order for neutrons
c
c         ecut      -   cut-off energy above lambda (in hbar*omega)
c
c         iendd     -   0 (continue), 1 (end)
c         ipart     -   1  calculate smoothed particle number
c         idens     -   1  calculate smoothed densities
c
c
       namelist /strutek/ swbinp,shellout,shellres,shellrho,iz,in,
     *                    gampro,gamneu,ippro,ipneu,
     *                    ecut,enei,enef,dene,ilev,iold,iendd,ipart
      common /radii/rmspsm,rmspqm,rmsnsm,rmsnqm,rmssm,rmsqm,
     &              rdmspsm,rdmspqm,rdmsnsm,rdmsnqm,rdmssm,rdmsqm
c
       enei=-50.d0
       enef=30.d0
       dene=0.1d0
       ecut=60.2
       ilev=0
       iold=0
       ippro=12
       ipneu=12
       gampro=2.1
       gamneu=2.1
       iendd=0
       ipart = 0
       shellout = 'shell.out'
       shellres = 'shell.res'
       shellsmooth = 'shell.smooth'
       shellrho = 'shell.rho'
       title = '            '
       open(unit=5,file='levden.inp',status='old')  
c
       ifirst=1
 5     continue
       read(5,strutek)
       write(*,strutek)
c
       if (iendd.eq.1) stop 'allright'
c
       if (ifirst.eq.1) then
         open(nf1,file=swbinp,status='old')
         open(nf2,file=shellres,form='formatted',status='unknown')
         open(6,file=shellout,form='formatted',status='unknown')
	 if(ipart.eq.1) 
     &      open(9,file=shellsmooth,form='formatted',status='unknown')
         open(6,file=shellout,form='formatted',status='unknown')
         write(nf2,'(a,i4,a,i4)')
     &  ' shell corrections for Z=',iz,' and N=',in
         write(nf2,123)
 123   format(' force      ',8x,
     1 '  ip   gam    gamma    Eshell   Esmooth    xlamb   Esp',/
     2,20x,'      [hom0]  [MeV]    [MeV]     [MeV]     [MeV]   [MeV]')   
       ifirst=0
       else
         rewind(nf1)      
       endif
c
c       rewind nf1
c
       titsav = title
       call swbrea(enipr,j2pr,ngorp,nf1,tit le)  ! read protons
       call swbrea(enine,j2ne,ngorn,nf1,title)  ! read neutrons
       call swbrea(eni0pr,j2dum,ngor0p,nf1,titdum)  ! read proton gas 
       call swbrea(eni0ne,j2dum,ngor0n,nf1,titdum)  ! read neutron gas 
       rewind(nf1)
c
       write(6,*)'iz in',iz,in
       aa=iz+in
c
        h=41.d0/aa**0.33333333333d0
c
       if(ilev.gt.0)open(unit=2,file='levden.plt',status='unknown')
100    format(1x,8a4)
c
       do iprot=1,0,-1
c
c             strutinski smooth energy
c
       if(iprot.eq.0) then
         gam   = gamneu
         ngor  = ngorn
         ip    = ippro
         do i=1,ngor
           eni(i)  = enine(i)
           eni0(i) = eni0ne(i)
           j2act(i) = j2ne(i)
         enddo
       else
         gam   = gampro
         ngor  = ngorp
         ip    = ipneu
         do i=1,ngor
           eni(i)  = enipr(i)
           eni0(i) = eni0pr(i)
           j2act(i) = j2pr(i)
         enddo
       endif
       gamm=gam      
       write(6,'(a,i3)') ' call SHCORR: iprot=',iprot
       call shcorr(iz,in,iprot,esmoo,ngor,esngp,xlamb,ip,ecut,gamm)
c
       gammq=gamm*h
       if(ilev.gt.0) then
         nen=(enef-enei)/dene+1
         if(nen.gt.2000)then
           write(*,130)nen
 130       format('    nen=',i6,'  is > 2000')
           stop
         endif
         if(ipart.eq.1)then
c          calculate smoothed particle numbers for each orbits	
	   call partsm
     &       (xlamb,eni,eni0,her,dher,ip,cm,ngor,inucl,gammq,ps,ps0)
           sumnpart=0.d0
           sumtest = 0.0D0
	   do i=1,ngor
             write(*,*) ' iprot,igor,energy,ps=',
     &       iprot,i,eni(i),ps(i),ps0(i)
             sumnpart=sumnpart+ps(i)-ps0(i)
             sumtest = sumtest+ps(i)
           enddo
           write(*,*)' sumnpart,sumtest=',sumnpart,sumtest
           call smodens(ps,ps0,ngor,iprot,iz,in,shellrho)
         endif
       endif
       i=0
       do ie=1,nen
       i=i+1
       ene=enei+(ie-1)*dene 
       zz=hg(ene)
       if (iold.eq.0) then
          z0=hg0(ene)
          cg(i)=zz-z0
          write(2,105)ene,zz-z0
       else
         cg(i)=zz
          write(2,105)ene,zz
       endif
105    format(1x,f12.4,3e12.4)
c
c      end of the loop for the energy
c
       enddo
       write(2,*)
        gm=gammq
        eshell=esngp-esmoo
        if(iprot.eq.0)then
        write(nf2,124)title,ip,gam,gm,eshell,esmoo,xlamb,esngp
        else
        write(nf2,114)title,ip,gam,gm,eshell,esmoo,xlamb,esngp
        endif
115     format(1x,f8.4,f12.4,f10.4,f10.4)
124     format(1x,a,'neutron ',i3,f7.2,f8.4,2f10.3,2f10.3,
     1 f10.3,f9.3,f10.3)
114     format(1x,a,'proton  ',i3,f7.2,f8.4,2f10.3,2f10.3,
     1 f10.3,f9.3,f10.3)
       if(iprot.EQ.0) write(nf2,'(a,3(/a,4(1pg13.5)))')
     &  '   radii:   rms_smooth    rdms_smooth     rms_qm      rdms_qm',
     &  '  proton:',rmspsm,rdmspsm,rmspqm,rdmspqm,
     &  ' neutron:',rmsnsm,rdmsnsm,rmsnqm,rdmsnqm,
     &  '   total:',rmssm,rdmssm,rmsqm,rdmsqm
       aa=iz+in
       part=in
       if (iprot.eq.1) part=iz
c
c       call hhepcba(xlamb,in,iz,iprot,robar,gpair)
c
c       x=eshell
c       write(3,'(5f11.5)')gammq,x
c
c      write(nf2,'(5i7)')        iz,in,iz0,in0,iprot
c      write(nf2,'(5f10.5)')     beta2,beta3,beta4,beta5,beta6
c      write(nf2,'(2i7,2f11.5)') nmax,ip,gamm,ecut
c      write(nf2,'(4f11.5)')     esmoo,esngp,eshell,xlamb
       write(6,113) ip,gamm,esmoo,esngp,eshell,xlamb
113    format(1x,i4,f10.5,2f12.5,f10.5,f10.5)       
112    format(1x,i3,f8.4,5f10.3,f9.3)
      if (ilev.gt.0)close(2)
  2   format(1x,'z=',i3,'  n=',i3,4x,a)
  3   format(1x,'z0=',i3,'  n0=',i3,/,
     :       1x,'beta2=',f7.4,2x,
     :       '  beta3=',f7.4,2x,'beta4=',f7.4,2x,'beta5=',f7.4,2x,
     :       'beta6=',f7.4,2x,/,1x,'nmax=',i2,' p=',i2,' gamma=',f6.2,
     :       '  ecut=',f7.2,/,1x,'esmooth=',f10.3,'  esp=',f10.3,
     :       '  eshell=',f10.4,'  lambda=',f10.4,/,1x,
     :       'robar/a=',f10.5,'  gpair=',f10.5,/,1x,75(1h*),/)
      enddo
c
      goto 5
c
      end
c
*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
cold      subroutine swbrea(ener,ilelev,nlevm,ndmain,nr,title)
      subroutine swbrea(ener,j2num,ilelev,nr,title)
*
*-----------------------------------------------------------------------
*
*     this routine reads energies produced by swbgigan (mamut version).
*     this is all in double precision. one set of data is read.
*        input  :  nr     - unit number of swbgigan result
*                           (make sure ener is properly dimensioned
*                            in the calling program)
*        output :  ener   - ordered vector of energies
*                  j2num  - 2*j of the state
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
      dimension  ener(*),j2num(*)
      character*20 aline
      character*12 title
*
      ilelev=0
*
      do iread=1,9999
        read(nr,'(a)',err=98,end=98) aline
        if(aline(1:4).eq.'  2j') goto 999
        if(aline(1:4).eq.' * *') then
          read(nr,'(a)',err=98,end=98) title
        endif
      enddo
 999  continue
*
      ind=0
cold 9    read(nr,101)j2,l,energ
 9    read(nr,*)j2,l,energ
      ind=ind+1
      ener(ind)=energ
      j2num(ind)=abs(j2)+1000*l
c      write(6,*) 'j2,l,energ=',j2,l,energ,j2num(ind)
      if (j2.lt.0) goto 19
      goto 9
 19   continue
      ilelev=ind
      write(*,*)'ind ',ind
*
      call ord(ener,ilelev)
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
      subroutine ord(ra,n)
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
      dimension ra(1)
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
         else
            rra=ra(ir)
            ra(ir)=ra(1)
            ir=ir-1
            if (ir.eq.1) then
               ra(1)=rra
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
               i=j
               j=j+j
            else
               j=ir+1
            end if
            go to 2
         end if
         ra(i)=rra
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
       common/bcs/ eni(25000),j2act(25000)
     &       /bcs0/eni0(25000)/old/iold
     *       /ome/   aa
c
c          shcorr calculates the shell correction to the liquid drop
c          model formula
c
c          h is the oscillator freqeuncy
c
      write(6,'(a,i3)') ' enter SHCORR: iprot=',iprot
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
      write(6,'(a,2i5)') ' leave SHCORR: iprot,ilev=',iprot,ilev
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
       common/bcs/ eni(25000),j2act(25000)
     &       /bcs0/eni0(25000)/old/iold
      external
     :         xn
      data
     :     pi,iter,eps /3.1415926536d0,2000,1.0d-7/

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
       common/bcs/ eni(25000),j2act(25000)
     &       /bcs0/eni0(25000)/old/iold
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
     :          eni(1),eni0(1),her(1),dher(1),cm(1)
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
c---------------------------------------------------------------------
c
      subroutine rtni(x,f,derf,fct,xst,eps,iend,ier)
      implicit double precision (a-h,o-z)
      external fct
c
c        x     - root of equation f(x)=0
c
c        f     function value at  x
c        derf  value of derivative at  x
c
c        fct   external subroutine fct(x,f,derf).  it computes
c              function value f  and derivative  derf  at  x
c
c        xst   initial guess of  x
c        eps   the upper bound of the error of  0
c        iend  max. number of iterations
c
c        ier   resultant error parameter
c              ier=0 - no errors
c                  1 - no convergence after iend iterations
c                  2 - at any iteration derf was equal to 0.
c
c
c        prepare iteration
c
      ier=0
      x=xst
      tdl=x
c
      call fct(tdl,f,derf)
c
      tdlf=100.*eps
c
c         start iteration loop
c
      do 6 i=1,iend
c
      if (f) 1,7,1
c
c        equation is not satisfied by x
c
    1 if (derf) 2,8,2
c
c        iteration is possible
c
    2 dx=f/derf
      x=x-dx
      tdl=x
c
      call fct(tdl,f,derf)
c
c        test on satisfactory accuracy
c
      tdl=eps
      a=abs(x)
c
      if (a-1.) 4,4,3
c
    3 tdl=tdl*a
c
    4 if (abs(dx)-tdl) 5,5,6
    5 if (abs(f)-tdlf) 7,7,6
c
    6 continue
c
c         end of iteration loop
c
c         no convergence after iend iteration steps, error return.
c
      ier=1
c
    7 return
c
c        error return in case of zero divisor
c
    8 ier=2
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
       common/bcs/ eni(25000),j2act(25000)
     *       /bcs0/eni0(25000)
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
      function hg0(e)
      implicit double precision (a-h,o-z)
      common
     *       /strut/ her(40),dher(40),cm(40),gammaq,
     *               ipher,ilev,inucl
       common/bcs/ eni(25000),j2act(25000)
     *       /bcs0/eni0(25000)
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

c------------------------------------------------------------------------
      subroutine partsm
     &    (e,eni,eni0,her,dher,ipher,cm,ilev,inucl,gammaq,ps,ps0)
c
c  calculates smoothed particle number for each orbits (i) for the
c  HF : ps(i) and for the gas: ps0(i) separately
c
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



c-----smodens-------------------------------------------------------------
c
      subroutine smodens(ps,ps0,ngor,iprot,nprot,nneut,shellrho)
c
      implicit double precision (a-h,o-z)
c
c     printing smoothed occupation number,
c     accumulating smoothed density,
c     and computing radii from that.
c
      parameter (KGRID2=400)    ! double size of coordinate space grid
      dimension ps(25000),ps0(25000)
      dimension rho(KGRID2),rho0(KGRID2),rhoq(KGRID2),sj(KGRID2)
      dimension formf(KGRID2),rhochr(KGRID2)
      save rho,rho0,rhoq,sj
c
      character*72 shellrho
c    
      common/bcs/ eni(25000),j2act(25000)
      common /radii/rmspsm,rmspqm,rmsnsm,rmsnqm,rmssm,rmsqm,
     &              rdmspsm,rdmspqm,rdmsnsm,rdmsnqm,rdmssm,rdmsqm
c
c-----------------------------------------------------------------------
c
c     preparatory: open and position file, reset accumulators
c
      open(28,file='for028')
      call spooliso(iprot,ngrid,rstep)
      if(2*ngrid.GT.KGRID2) 
     &    stop ' SMODENS: too large NGRID. enhance KGRID2.'
      iso = 2-iprot
      isoadd = ngrid*iso-ngrid
      do ir=isoadd+1,isoadd+ngrid
        rho(ir)  = 0.0D0
        rho0(ir) = 0.0D0
        rhoq(ir) = 0.0D0
        sj(ir)   = 0.0D0
      enddo
      if(iprot.EQ.1) then
        write(9,'(a)') ' smoothed occupations for protons:'
      else
        write(9,'(a)') ' smoothed occupations for neutrons:'
      endif
      multi = 1
      iniv  = 0
      ensav = -1000.0D0
      j2sav = -1
      do i=1,ngor
ctest        write(*,*) ' igor,energy=',i,eni(i)
        if(ensav.NE.eni(i) .or. j2sav.NE.j2act(i)) then
          if(i.GT.1) then
            iniv = 1+iniv
            occup  = ps(i-1)/2.0D0
            occup0 = ps0(i-1)/2.0D0
            write(9,'(2i5,3(1pg13.5))') 
     &         iniv,multi,eni(i-1),occup,occup0
            call accdens(rho,rho0,rhoq,ngrid,iso,multi,occup,occup0)
          endif
          multi = 1
          ensav = eni(i)
          j2sav = j2act(i)
        elseif(i.EQ.ngor) then
          multi = 1+ multi
          iniv = 1+iniv
          occup  = ps(i)/2.0D0
          occup0 = ps0(i)/2.0D0
          write(9,'(2i5,3(1pg13.5))') 
     &       iniv,multi,eni(i),occup,occup0
          call accdens(rho,rho0,rhoq,ngrid,iso,multi,occup,occup0)
        else
          multi = 1+ multi
        endif
      enddo
      close(28)
c
      if(iprot.NE.0) return
c
c     final print
c
      open(28,file=shellrho)
c
      write(28,'(a,2i5)')
     &  '# smoothed density for nprot,nneut=',nprot,nneut
      call radius(rho,rho(ngrid+1),sj,sj(ngrid+1),
     &            nprot,nneut,rstep,ngrid,20.75D0,20.75D0,0.0D0,0,
     &            -1,0,0,formf,rhochr,
     &                  rmsc,rmspsm,rmsnsm,rmssm,
     &                  rdmsc,rdmspsm,rdmsnsm,rdmssm,
     &                  surfc,surfp,surfn,surf,
     &                  smomc,smomp,smomn,smom)
      write(28,'(a/a,3(/a,3(1pg13.5)))')
     &  '# smoothed density:',
     &  '#             rms      rdms      surf   ',
     &  '# proton: ',rmspsm,rdmspsm,surfp,
     &  '# neutron:',rmsnsm,rdmsnsm,surfn,
     &  '# total:  ',rmssm,rdmssm,surf
      call radius(rhoq,rhoq(ngrid+1),sj,sj(ngrid+1),
     &            nprot,nneut,rstep,ngrid,20.75D0,20.75D0,0.0D0,0,
     &            -1,0,0,formf,rhochr,
     &                  rmsc,rmspqm,rmsnqm,rmsqm,
     &                  rdmsc,rdmspqm,rdmsnqm,rdmsqm,
     &                  surfc,surfp,surfn,surf,
     &                  smomc,smomp,smomn,smom)
      write(28,'(a/a,3(/a,3(1pg13.5)))')
     &  '# quantum density:',
     &  '#             rms      rdms      surf   ',
     &  '# proton: ',rmspqm,rdmspqm,surfp,
     &  '# neutron:',rmsnqm,rdmsnqm,surfn,
     &  '# total:  ',rmsqm,rdmsqm,surf
c
      write(28,'(a)') '#  ir  r  rhosm_p  rhoqm_p  rhosm_n  rhoqm_n '
      do ir=1,ngrid
        write(28,'(1x,i5,f8.2,6(1pg13.5))')
     &   ir,(ir-1)*rstep,rho(ir),rho0(ir),rhoq(ir),
     &    rho(ir+ngrid),rho0(ir+ngrid),rhoq(ir+ngrid)
      enddo
      close(28)
c
      return
      end

c-----accdens-------------------------------------------------------------
c
      subroutine accdens(rhosm,rhosm0,rhoqm,ngrid,iso,multi,
     &                   occup,occup0)
c
      implicit double precision (a-h,o-z)
c
c     Accumulate densities:
c      rhosm   = smoothed density, weighted with smoot occupation 'occup'
c      rhoqm   = q.m. density, weighted with smoot occupation 'occup0'
c      ngrid   = wanted grid size 
c      iso     = type of nucleon
c      multi   = multiplicity of state
c      occup   = smoothed occupation
c      occup0  = q.m. occupation
c     The s.p. density is read from 'for028'. Quantum numbers and
c     grid size are checked. The file 'for028' needs to be sorted
c     in the same order as 'for027' where the s.p. energies came from.
c
      logical ttest
      parameter (ttest=.true.)
      parameter (KGRID2=400)      ! size of coordinate space grid
      dimension rhosm(KGRID2),rhosm0(KGRID2),rhoqm(KGRID2)
c
      character*80 inline
c
c-----------------------------------------------------------------------
c
      do iread=1,9999
        read(28,'(a)') inline
        if(inline(1:7).EQ.' state:') then
          read(28,*) isoin,n,xj,l,spenerg,occupin,ngridin
          goto 99
        endif
      enddo
      stop ' ACCDENS: file count exhausted'
c
 99   continue
      if(ttest) then
        write(*,*) ' isoin,n,xj,l,spenerg,occupin,ngridin=',
     &             isoin,n,xj,l,spenerg,occupin,ngridin
        write(*,*) ' iso,ngrid,multi,occup,occup0=',
     &             iso,ngrid,multi,occup,occup0
      endif
      if(isoin.NE.iso) stop ' ACCDENS: isospin does not match '
      if(ngridin.NE.ngrid) stop ' ACCDENS: grid does not match '
c      if(occupin.NE.occup0) stop ' ACCDENS: occup does not match '
      multiin = xj+0.500000001D0    !    2.0D0*xj+1.000000000001D0
      if(multiin.NE.multi) 
     &   stop ' ACCDENS: multiplicity does not match '
c
      isoadd = ngrid*iso-ngrid
      fac    = occup*multi*2
      fac0   = (occup-occup0)*multi*2
      facq   = occupin*multi*2
ctest      sum = 0.0
      read(28,'(1x)')
      do ir=isoadd+1,isoadd+ngrid
        read(28,*) rin,rhosp
        rhosm(ir)  = rhosm(ir) + fac*rhosp        
        rhosm0(ir) = rhosm0(ir) + fac0*rhosp        
        rhoqm(ir)  = rhoqm(ir) + facq*rhosp        
ctest        sum = rin*rin*rhosp + sum
      enddo
      if(ttest) then
        write(6,'(a,2i5,f5.1,i5,3g13.5)')
     &  ' density accumulated for iso,n,j,l,spenerg,occup=',
     &  isoin,n,xj,l,spenerg,occupin,occup
ctest        write(6,'(a,g12.4)') ' partial p.num.=',sum*1.2*3.1415926
c        write(*,'(5g12.4)')
c     &   (rhosm(ir),ir=isoadd+1,isoadd+ngrid)
      endif 
c
      return
      end

c-----spooliso-------------------------------------------------------
c
      subroutine spooliso(iprot,ngrid,rstep)
c
      implicit double precision (a-h,o-z)
c
      character*80 inline
c
c---------------------------------------------------------------------
c
      do iread=1,99999
        read(28,'(a)') inline
        if(inline(1:7).EQ.' state:') then
          read(28,*) isoin,n,xj,l,spenerg,occup,ngrid
          if(isoin.EQ.2-iprot) goto 99
          do i=1,ngrid+1
            read(28,'(1x)')
          enddo
        endif
      enddo
      stop ' SPOOLISO: file count exhausted '
c
 99   continue
      read(28,'(1x)')
      read(28,*) r0,dummy
      read(28,*) r1,dummy
      rstep = r1-r0
      do iback=1,5
        backspace(28)
      enddo
c
      return
      end


c-----radius-----------------------------------------------------------
c
      subroutine radius(rhop,rhon,sjsp,sjsn,
     &                  nprot,nneut,rstep,ngrid,amfpr,amfne,zpe,ifzpe,
     &                  iwr,ifrso,ifform,
     &                  formf,rhochr,
     &                  rmsc,rmsp,rmsn,rms,
     &                  rdmsc,rdmsp,rdmsn,rdms,
     &                  surfc,surfp,surfn,surf,
     &                  smomc,smomp,smomn,smom)
c
      include 'params.inc'
c      implicit double precision (a-h,o-z)
c
c     * * * * *    rms-radius, dms-radius and surface-thickness    * * * * *
c     * * * *    for proton-, neutron-, total- and charge-density    * * * *
c     * * * * * * *    given on an equidistant radial grid     * * * * * * *
c
c     the charge-density includes folding with intrinsic electric and
c     magentic formfactors of proton and neutron.
c     the contribution of the anomalous magnetic moment is included
c     only if ifrso=1.
c
c     Input variables are:
c      rhop    = array of radial proton density (rhop(1) --> at r=0)
c      rhon    = array of radial neutron density
c      sjsp    = Nabla of spin-orbit current for protons (used for ifrso=1)
c      sjsn    = Nabla of spin-orbit current for neutrons (used for ifrso=1)
c      nprot   = number of protons
c      nneut   = number of neutrons
c      rstep   = spacing of the (equidistant) radial grid
c      ngrid   = number of radial grid points
c      amfpr   = hbar^2/2m_proton
c      amfne   = hbar^2/2m_neutron
c      zpe     = c.m. corection energy (>0), in case of <P^2>/2mA recipe
c      ifzpe   = 1 -->  <P^2>/2mA recipe, requires input of 'zpe'
c                0 -->  m-m/A recipe
c      iwr     = print level
c                -2 suppresses any printing
c                -1 gives reasonable amount
c                > -1 only for testing
c      ifrso   = 1 --> chrarge density contains contribtion from
c                      spin-orbit current
c      ifform  = 1 --> nucleon formfactors are folded in
c                0 --> mere proton density is analyzed
c
c     Output variables are:
c      formf   = charge formfactor
c      rhochr  = charge density inluding all nucleon formfactors
c      rmsc    = charge r.m.s. radius
c      rmsp    = proton r.m.s. radius
c      rmsn    = neutron r.m.s. radius
c      rms     = total r.m.s. radius
c      rdmsc   = charge diffraction radius (sharp- or box-radius)
c      rdmsp   = proton diffraction radius
c      rdmsn   = neutron diffraction radius
c      rdms    = total diffraction radius
c      surfc   = charge surface thickness
c      surfp   = proton surface thickness
c      surfn   = neutron surface thickness
c      surf    = total surface thickness
c
c     the following fields 'foc' etc are workspace for formfactors:
c      foc0   = charge formfactor, including folding only with electric
c               nucleon distribution
c      foc    = complete charge formfactor, including also the folding with
c               the anomalous magnetic moment
c      fop    = formfactor of the proton distribution
c      fon    = formfactor of the neutorn distribution
c      fo     = formfactor of the total nucleon distribution
c
c      parameter (kfrm=50)        ! maximum dimension of formfactor fields
c
      dimension foc(kfrm),fop(kfrm),fon(kfrm),fo(kfrm),foc0(kfrm)
c
      dimension rhop(ngrid),rhon(ngrid)         ! input information
      dimension sjsp(ngrid),sjsn(ngrid)         ! input information
      dimension rhochr(ngrid),formf(kfrm)       ! results
c
      common /radpar/delq,rleng,iqmax           ! to communicate to subr.
c
c      data  zero,one,pi/0.0,1.0,3.1415926/
      data  zero,half,one,pi/0.0,0.5,1.0,3.1415926/
c
c-------------------------------------------------------------------------------
c
c     preparation of auxiliary parameters
c
      anumb  = one*(nprot+nneut)
c      kqmax  = min(kfrm,ngrid*3/4)       ! ???
c      iqmax  = min(kfrm,ngrid*3/4)
      iqmax  = ngrid-1
      rleng  = (ngrid-1)*rstep
      delq   = 3.1415926/rleng
      frpirs = 4*pi*rstep
c
c     compute 'cmwid' as width for for c.m.correction on the formfactor
c
      if(iwr.ge.-1) write(6,'(a,g13.5)') ' ZPE=',zpe
      if(ifzpe.eq.0) then
       cmwid  = 1.58/(6.0*anumb**0.66666667)
      else
       cmwid  = 0.375*0.5*(amfpr+amfne)/(anumb*zpe)
      endif
c
c      formfactor by fourier-bessel transformation
c
c                                                 q=0.0 separately
      sump   = zero
      sumn   = zero
      sum2p = zero
      sum2n = zero
c      sumsp  = zero
      r      = zero
      do ir=2,ngrid
        r     = rstep+r
        sump   = r*r*rhop(ir)+sump
        sumn   = r*r*rhon(ir)+sumn
        sum2p   = r**4*rhop(ir)+sum2p
        sum2n   = r**4*rhon(ir)+sum2n
      enddo
      write(6,'(a,2(1pg13.5))')
     &  ' direct: rmsp,rmsn=',sqrt(sum2p/sump),sqrt(sum2n/sumn)
      foc0(1) = frpirs*sump
      fop(1)  = foc0(1)
      fon(1)  = frpirs*sumn
      fo(1)   = fop(1)+fon(1)
      foc(1)  = foc0(1)
c                                                      all other q
      do iq=2,iqmax
        q1     = iq*delq-delq
        call fourf(q1,
     &             rhop,rhon,sjsp,sjsn,
     &             rstep,ngrid,ifrso,ifform,cmwid,
     &             foc0(iq),fop(iq),fon(iq),fo(iq),foc(iq))
      enddo
c
c     finding first diffraction radii 'rdms..' at first zero
c     of the formfactor,
c     and second diffraction radii 'rsec..' at second zero.
c
c     the endings indicate the type of radius:
c      ....c  = radius of the complete charge distribution
c      ....c0 = radius of the charge distribution without spin-orbit effects
c      ....p  = radius of the proton distribution
c      ....n  = radius of the neutron distribution
c      ....   = radius of the total mass distribution
c
      rdmsc  = 4.493/qzero(foc,1)
      rsecc  = 7.7253/qzero(foc,int(4.493/(delq*rdmsc))+2)
      if(iwr.ge.-1) then
       rdmsc0 = 4.493/qzero(foc0,1)
       rsecc0 = 7.7253/qzero(foc0,int(4.493/(delq*rdmsc0))+2)
       rdms   = 4.493/qzero(fo,1)
       rsec   = 7.7253/qzero(fo,int(4.493/(delq*rdms))+2)
       if(nprot.gt.0) then
        rdmsp  = 4.493/qzero(fop,1)
        rsecp  = 7.7253/qzero(fop,int(4.493/(delq*rdmsp))+2)
       else
        rdmsp  = zero
        rsecp  = zero
       endif
       if(nneut.gt.0) then
        rdmsn  = 4.493/qzero(fon,1)
        rsecn  = 7.7253/qzero(fon,int(4.493/(delq*rdmsn))+2)
       else
        rdmsn  = zero
        rsecn  = zero
       endif
      endif
c      rdmsca=rdmsc/anumb**.3333333
c
c      evaluating surface width from height of first maximum.
c      the 'q1'="momentum_of_first_maximum" is estimated from
c      'rdmsc' as in a square-well density-distribution.
c
      q1     = 5.600/rdmsc
      if(nprot.gt.0) then
       qr     = q1*rdmsc
       surfc  = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &              *3.*nprot/(qr*qr*fbint(q1,foc))))))
       surfc0 = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &            *3.*nprot/(qr*qr*fbint(q1,foc0))))))
      else
       surfc  = zero
       surfc0 = zero
      endif
      if(iwr.ge.-1) then
       qr     = q1*rdms
       surf   = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &           *3.*(nprot+nneut)/(qr*qr*fbint(q1,fo))))))
       if(nprot.gt.0) then
        qr     = q1*rdmsp
        surfp  = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &                *3.*nprot/(qr*qr*fbint(q1,fop))))))
       else
        surfp  = zero
       endif
       if(nneut.gt.0) then
        qr     = q1*rdmsn
        surfn  = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &             *3.*nneut/(qr*qr*fbint(q1,fon))))))
       else
        surfn  = zero
       endif
      endif
c
c     rms-radii from curvature of formfactor at q=0
c     ---> extremely dangerous procedure due to insufficient interpolations
c     ----> values will be overwritten later by direct integration
c
      dqt    = delq/10.0     !  delq/4.0
      dqtd2  = 6.0/(dqt*dqt)
      rmsc   = sqrt(abs(dqtd2*(1.25*foc(1)
     &                  -1.33333333*fbint(dqt,foc)
     &                  +fbint((dqt+dqt),foc)/12.0)/foc(1)))
      if(iwr.ge.-1) then
       rmsc0  = sqrt(abs(dqtd2*(1.25*foc0(1)
     &                   -1.33333333*fbint(dqt,foc0)
     &                +fbint((dqt+dqt),foc0)/12.0)/foc0(1)))
       if(nprot.gt.0) then
        rmsp   = sqrt(abs(dqtd2*(1.25*fop(1)
     &                   -1.33333333*fbint(dqt,fop)
     &                +fbint((dqt+dqt),fop)/12.0)/fop(1)))
       else
        rmsp   = zero
       endif
       if(nneut.gt.0) then
        rmsn   = sqrt(abs(dqtd2*(1.25*fon(1)
     &                   -1.33333333*fbint(dqt,fon)
     &                +fbint((dqt+dqt),fon)/12.0)/fon(1)))
       else
        rmsn   = zero
       endif
        rms    = sqrt(abs(dqtd2*(1.25*fo(1)
     &                   -1.33333333*fbint(dqt,fo)
     &                +fbint((dqt+dqt),fo)/12.0)/fo(1)))
      endif
c
c     save complete charge formfactor on 'formf'.
c
      do iq=1,iqmax
        formf(iq)=foc(iq)
      enddo
c
c     testprints of spin-orbit current
c
      if(ifrso.eq.1 .and. iwr.ge.0) then
       write(06,'(a,20(/1x,10g12.4))') ' foc0:',
     &                 (foc0(iq),iq=1,iqmax)
       write(06,'(a,20(/1x,10g12.4))') ' foc:',
     &                 (foc(iq),iq=1,iqmax)
      endif
c
c     charge density by fourier back-transformation
c
      dpi3   = half/(pi*pi)
c                                                         r=0.0 separately
      sum    = zero
      do iq=2,iqmax
        qa     = delq*iq-delq
        sum    = qa*qa*formf(iq)+sum
      enddo
      rhochr(1) = abs(sum*delq*dpi3)
c                                                         now all other r
      r      = zero
      do ir=2,ngrid
        r      = rstep+r
c
        dc     = cos(r*delq)
        ds     = sin(r*delq)
        c      = one
        s      = zero
        sum    = zero
        qa     = zero
        do iq=2,iqmax
          qa     = delq+qa
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
          sum    = qa*s*foc(iq)+sum
        enddo
        rhochr(ir) = abs(sum*delq*dpi3/r)
c
      enddo
c
c     surface width from moments
c
      if(iwr.ge.-1) 
     &  write(6,'(a,4g13.5)')
     &  ' from F(q): rms,rmsp,rmsn,rmsc=', rms,rmsp,rmsn,rmsc
c
c        rmom0p = zero
c        rmom0n = zero
      rmom2p = zero
      rmom4p = zero
      rmom2n = zero
      rmom4n = zero
      rmom2c = zero
      rmom4c = zero
      r      = zero
      do ir=2,ngrid
        r     = rstep+r
        r2    = r*r
        r4    = r2*r2
        r6    = r4*r2
c        rmom0p = r2*rhop(ir) + rmom0p
        rmom2p = r4*rhop(ir) + rmom2p
        rmom4p = r6*rhop(ir) + rmom4p
c        rmom0n = r2*rhon(ir) + rmom0n
        rmom2n = r4*rhon(ir) + rmom2n
        rmom4n = r6*rhon(ir) + rmom4n
        rmom2c = r4*rhochr(ir) + rmom2c
        rmom4c = r6*rhochr(ir) + rmom4c
      enddo
      rmom2  = frpirs*(rmom2p+rmom2n)/fo(1)
      rmom4  = frpirs*(rmom4p+rmom4n)/fo(1)
      rmom2p = frpirs*rmom2p/fop(1)
      rmom4p = frpirs*rmom4p/fop(1)
      rmom2n = frpirs*rmom2n/fon(1)
      rmom4n = frpirs*rmom4n/fon(1)
      rmom2c = frpirs*rmom2c/foc(1)
      rmom4c = frpirs*rmom4c/foc(1)
      rms    = sqrt(rmom2)
      rmsp   = sqrt(rmom2p)
      rmsn   = sqrt(rmom2n)
      rmsc   = sqrt(rmom2c)
c
      if(iwr.ge.-1) then
        write(6,'(a,6g13.5)')
     &   ' rmom2p,rmom4p,rmom2n,rmom4n,rmom2c,rmom4c=',
     &     rmom2p,rmom4p,rmom2n,rmom4n,rmom2c,rmom4c
        smomp  = (rmom4p- rmom2p*rmom2p*25.0/21.0)**0.25
        smomn  = (rmom4n- rmom2n*rmom2n*25.0/21.0)**0.25
        smomc  = (rmom4c- rmom2c*rmom2c*25.0/21.0)**0.25
        smom   = (rmom4 - rmom2*rmom2*25.0/21.0)**0.25
c
c       print results
c
ctest        write(6,'(a,4g13.5)')
ctest     &   ' fo(1)=',fop(1),fon(1),fo(1),foc(1),
ctest     &   ' dir. =',frpirs*rmom0p,frpirs*rmom0n,frpirs*(rmom0p+rmom0p)
ctest          rmspt = (fop(1)-fop(2))/fop(1)*6.0/delq**2
ctest          rmsnt = (fon(1)-fon(2))/fon(1)*6.0/delq**2
ctest        write(6,'(a,2g13.5)') ' rmpst,rmsnt=',sqrt(rmspt),sqrt(rmsnt)
ctest          rmsnt = sqrt(abs((6.0/delq**2)*
ctest     &                     (1.25*fon(1)-1.33333333*fon(2)
ctest     &                                 +fon(3)/12.0      )/fon(1)))
ctest          rmspt = sqrt(abs((6.0/delq**2)*
ctest     &                     (1.25*fop(1)-1.33333333*fop(2)
ctest     &                                 +fop(3)/12.0      )/fop(1)))
ctest        write(6,'(a,2g13.5)') ' rmpst,rmsnt=',rmspt,rmsnt
ctest        do it=0,8
ctest          dqt    = delq/2**it
ctest          dqtd2  = 6.0/(dqt*dqt)
ctest          rmspt3 = sqrt(abs(dqtd2*(1.25*fop(1)
ctest     &                   -1.33333333*fbint(dqt,fop)
ctest     &                +fbint((dqt+dqt),fop)/12.0)/fop(1)))
ctest          rmspt2  = sqrt(dqtd2*(fop(1)-fbint(dqt,fop))/fop(1))
ctest          write(6,'(a,i3,g12.4,2g13.5)')
ctest     &      ' it,delq,rmspt=',it,dqt,rmspt2,rmspt3
ctest        enddo
c
        write(06,'(16x,a,5(/1x,a,5g13.4))')
     &   ' diff.radius    surface     sec.radius   rms-radius surf.mom',
     &   ' charge density:',rdmsc,surfc,rsecc,rmsc,smomc,
     &   '   mass density:',rdms,surf,rsec,rms,smom,
     &   ' proton density:',rdmsp,surfp,rsecp,rmsp,smomp,
     &   'neutron density:',rdmsn,surfn,rsecn,rmsn,smomn,
     &   ' charge no s.o.:',rdmsc0,surfc0,rsecc0,rmsc0
       write(6,'(a/4(a,g13.5))')
     &   ' rms-radii by direct integration:',
     &   ' charge=',sqrt(rmom2c),'  mass=',sqrt(rmom2),
     &   '  proton=',sqrt(rmom2p),'  neutron=',sqrt(rmom2n)
       rmsn   = sqrt(rmom2n)
       rmsp   = sqrt(rmom2p)
c
        write(06,'(1h0,a,50(/1x,10g13.5))')
     &     'testprint rhochr:', (rhochr(i),i=1,ngrid)
c
      endif
c
      return
      end
c-----fourf ------------ part of the hartree-fock package ----------------------
c
      subroutine fourf(q1,
     &                 rhop,rhon,sjsp,sjsn,
     &                 rstep,ngrid,ifrso,ifform,cmwid,
     &                 foc0,fop,fon,fo,foc)
c
      implicit double precision (a-h,o-z)
c
c     Input variables are:
c      q1    : radial momentum at which the fourier-transformation
c              is to be evaluated
c      rhop    = array of radial proton density (rhop(1) --> at r=0)
c      rhon    = array of radial neutron density
c      sjsp    = Nabla of spin-orbit current for protons (used for ifrso=1)
c      sjsn    = Nabla of spin-orbit current for neutrons (used for ifrso=1)
c      rstep   = spacing of the (equidistant) radial grid
c      ngrid   = number of radial grid points
c      ifrso   = 1 --> chrarge density contains contribtion from
c                      spin-orbit current
c      ifform  = 1 --> nucleon formfactors are folded in
c                0 --> mere proton density is analyzed
c      cmwid : parameter for center of mass correction
c     output paramaters:
c      foc0 :   charge formfactor
c      fop  :   proton "
c      fon  :   neutron    "
c      fo   :   total      "
c      foc  :   charge     "   including spin-orbit current
c
      dimension rhop(ngrid),rhon(ngrid)         ! input information
      dimension sjsp(ngrid),sjsn(ngrid)         ! input information
c
c      data  zero,one,pi/0.0,1.0,3.1415926/
      data  zero,half,one,pi/0.0,0.5,1.0,3.1415926/
c
c      the parameters in data are for the sachs-formfactors:
c       sa..,sm.. are for the isoscalar electric formfactor
c       va..,vm.. are for the isovector electric formfactor
c       gp..,gm.. are for the magnetic formfactors of proton and neutron
c
c       isoscalar and isovector formfactors are recoupled to
c       the electric formfactors of proton and neutron.
c
c     the sachs-formfactors are then recoupled to the nucleon
c     vector formfactor (form1..) and the nucleon tensor
c     formfactor (form2..).
c
c
c     parameters for scalar and vector electric formfactor
c     (walther, priv.comm., sept. 86)
c
      data    sa1,    sa2,    sa3,   sa4,   sm1,  sm2,  sm3,  sm4
     *    /2.2907,-0.6777,-0.7923,0.1793, 15.75,26.68,41.01,134.2/
      data    va1,    va2,    va3,   va4,   vm1,  vm2,  vm3,  vm4
     *    /0.3681, 1.2263,-0.6316,0.0372,  5.00,15.02,44.08,154.2/
c
c     magnetic formfactor , simon et al n.p. a333 (1980) 381
c
      data    gp1,    gp2,    gp3,   gp4,   gm1,  gm2,  gm3,  gm4
     *    / 0.694,  0.719,-0.418 ,0.005 ,  8.50,15.02,44.08,355.4/
c
c     old version of the
c     parameters for the magnetic formfactor (arenhoevel 1980)
c
c      data gp1,gp2,gp3,gp4,gm1,gm2,gm3,gm4/.794,.594,-.393,.005,
c     *                      8.73,15.02,44.08,355.40/
c
c     the anomanlous magnetic moments and the darwin factor 'd4m'
c
      data apro0,apro11,apro12,bpro11
     *    / 0.74, 0.0  ,  0.0 ,  0.0  /
      data aneu0,aneu11,aneu12,bneu11
     *    / 0.74, 0.0  ,  0.0 ,  0.0  /
      data xmgne,xmgpr,d4m/.04197,-.06131,.010989/
c
c----------------------------------------------------------------
c
      frpirs = 4*pi*rstep
c
c     fourier-bessel transformation for proton- and neutron-density
c     given on field 'rho' on common /wavfun/.
c
      dc     = cos(q1*rstep)
      ds     = sin(q1*rstep)
      r      = zero
      sump   = zero
      sumn   = zero
      c      = one
      s      = zero
       do 10 i=2,ngrid
       r      = rstep+r
       x      = c*dc-s*ds
       s      = s*dc+c*ds
       c      = x
       sump   = r*s*rhop(i)+sump
       sumn   = r*s*rhon(i)+sumn
   10  continue
      sump   = frpirs*sump/q1
      sumn   = frpirs*sumn/q1
c
c     fourier-bessel transformation for spin-orbit density of
c     protons and neutrons given on field 'sjs'.
c
      if(ifrso.eq.1) then
       sumsp  = zero
       sumsn  = zero
       r      = zero
       c      = one
       s      = zero
        do 20 i=2,ngrid
        r      = rstep+r
        x      = c*dc-s*ds
        s      = s*dc+c*ds
        c      = x
        sumsp  = r*s*sjsp(i)+sumsp
        sumsn  = r*s*sjsn(i)+sumsn
   20   continue
       sumsp=frpirs*sumsp/q1
       sumsn=frpirs*sumsn/q1
      else
       sumsp  = zero
       sumsn  = zero
      endif
c
c     three cases for the nucleon folding:
c      ifform = 0  --> no folding is done
c      ifform = 1  --> folding with free nucleon formfactor parametrised
c                      by sum of lorentzians, see data on sa.., va..
c                      and gp..
c      ifform =-1  --> density dependent folding for charge density
c
      if(ifform.eq.0) then
       schsmp = 0.0
c       schcmn = 0.0
       sumc   = sump
      elseif(ifform.gt.0) then
c
c      isoscalar and isovector electrical formactors on 'schses'
c      and 'schsev', recoupled then to proton and neutron electric
c      formfactors  'schsep' and 'schsen'. proton and neutron
c      magnetic formfactor on 'schsmp' and 'schsmn'.
c
       qa     = q1*q1
       schses = sa1/(1.+qa/sm1)+sa2/(1.+qa/sm2)+sa3/(1.+qa/sm3)
     &         +sa4/(1.+qa/sm4)
       schsev = va1/(1.+qa/vm1)+va2/(1.+qa/vm2)+va3/(1.+qa/vm3)
     &         +va4/(1.+qa/vm4)
       schsep = half*(schses+schsev)
       schsen = half*(schses-schsev)
       schsm  = ( gp1/(one+qa/gm1)+gp2/(one+qa/gm2)+gp3/(one+qa/gm3)
     &           +gp4/(one+qa/gm4) )
       schsmp = schsm*xmgpr - schsep*d4m
       schsmn = schsm*xmgne - schsen*d4m
c
       sumc   = sump*schsep+sumn*schsen
c
      else
c
c      charge fromfactor with density-dependent intrinsic nucleon radii,
c      density dependence is given by the parameters 'apro..','bpro..',
c      'aneu..' and 'bneu..' transferred via common /rmsnuc/.
c
c       iradd = ngrid
c       if(ifform.eq.-2) iradd=0
       fq2    = q1*q1
       dc     = cos(q1*rstep)
       ds     = sin(q1*rstep)
       r      = zero
       sumc   = zero
       c      = one
       s      = zero
       do 30 i=2,ngrid
       r      = rstep+r
       x      = c*dc-s*ds
       s      = s*dc+c*ds
       c      = x
c                                   proton-radius as fkt. of density:
       rhotot = rhop(i)+rhon(i)
       rmspi2 = ( apro0+apro11*rhotot
     &           +apro12*rhotot*rhotot
     &           +bpro11*(rhop(i)-rhon(i)) )/6.0
c                                   neutron-radius as fkt. of density:
       rmsni2 = aneu0+aneu11*rhotot
     &         +aneu12*rhotot*rhotot
     &         +bneu11*(rhop(i)-rhon(i))
c                                   accumlate folding:
       fintp  = exp(-rmspi2*fq2)
       fintn  = 12./(12.+(0.64+rmsni2)*fq2)-12./(12.+(0.64-rmsni2)*fq2)
       sumc   = r*s*(rhop(i)*fintp + rhon(i)*fintn) + sumc
   30  continue
       sumc   = frpirs*sumc/q1
c
      endif
c
c      factor for darwin correction 'darfac'
c      and the factor for the c.m. correction 'exphf'.
c
       qa     = q1*q1
       if(ifrso.eq.1) then
        darfac = one/sqrt(one+qa*d4m)
       else
        darfac = one
       endif
       exphf  = exp(cmwid*qa)
c
c     finally assembling the formfactors
c
      foc0   = exphf*sumc
      fop    = exphf*sump
      fon    = exphf*sumn
      fo     = exphf*(sump+sumn)
      foc    = exphf*(sumc*darfac + (schsmp*sumsp+schsmn*sumsn))
c
      return
      end
c-----fbint ------------ part of the hartree-fock package ----------------------
c
      function fbint(qin,formfa)
c
      implicit double precision (a-h,o-z)
c
c     evaluates and returns formfactor at the momentum 'qeff=q*rleng' by
c     fourier-bessel interpolation
c     from formfactors at fourier-grid points given on array 'formfa'
c     in the span '1..iqmax'.
c
      parameter (kfrm=199)        ! maximum dimension of formfactor fields
      dimension formfa(kfrm)
c
      common /radpar/delq,rleng,iqmax           ! to communicate to subr.
c
c     'precis' defines the neighbourhood of a pole-term for which the
c     routine switches to an explicit laurent expansion.
c
      data precis/1.d-9/          ! /1.e-2/
      data  zero,one,pi/0.0d0,1.0d0,3.14159265358979d0/
c
c-----------------------------------------------------------------------
c
      if(iqmax.gt.kfrm) stop ' too large grid in FBINT '
c
      qeff   = qin*rleng
      xfx    = (qeff/pi)
      xfx2   = xfx*xfx
      ixfx   = (xfx+precis)
c
c     the fourier-bessel interpolation runs into problems with round-off
c     errors if 'qeff' comes too close to a multiple of 'pi'. we switch
c     to an explicit treatment of the pole-term in that case (i.e. first
c     case in the if-closure).
c
      if(abs(xfx-ixfx) .lt. precis) then
c                                             pole term separately
       ipole  = ixfx+1
       signum = -one
       sum    =  zero
       do i=2,iqmax
         signum = -signum
         if(i.ne.ipole)
     *    sum    = formfa(i)*signum/(one-xfx2/float((i-1)*(i-1))) + sum
       enddo
       signum = (mod(ipole,2)*2-1)
       qdif   = (xfx-ixfx)
       sum2   = signum*formfa(ipole)*(ixfx*ixfx)/((ixfx+ixfx)+qdif)
       qdif   = qdif*pi
       sinq   = (one-qdif*qdif*0.166666667)*signum
       fbint  = 2.0*(sum*qdif+sum2*pi)*sinq/qeff
      else
c                                             normal case
       signum = -one
       sum    = zero
       do i=2,iqmax
         signum = -signum
         sum    = formfa(i)*signum/(one-xfx2/float((i-1)*(i-1))) + sum
       enddo
       fbint  = (sum+sum)*sin(qeff)/qeff
      endif
c
      return
      end
c-----qzero ------------ part of the hartree-fock package ----------------------
c
      function qzero(foin,iqstrt)
c
      implicit double precision (a-h,o-z)
c
c     evaluates and returns the momentum at which the formfactor
c     given on array 'foin' has a zero.
c     scan for zero starts at the field position 'iqstrt'.
c     'iqmax' is the actual length of the array 'foin'.
c     'delq' is the stepsize in the momentum.
c
      parameter (kfrm=50)        ! maximum dimension of formfactor fields
      dimension foin(kfrm)
c
      common /radpar/delq,rleng,iqmax           ! to communicate to subr.
c
c     some numerical parameters for termination of the search
c
      data itq0mx,endq0/30,1.0D-9/
      data  zero,one/0.0D0,1.0D0/
c
c---------------------------------------------------------------------
c
c     scan changing sign in 'foin'.
c      'iq1' is position before the zero and 'iq1' after the zero.
c
      iqmaxm = iqmax-1
      do 10 iscan=iqstrt,iqmaxm
      if(foin(iscan)*foin(iscan+1).le.zero) goto 19
   10 continue
   19 continue
      iq1   = iscan
      iq2   = iq1+1
      if(foin(iq1).le.zero) iq2   = iq1-1
      q1    = iq1*delq-delq
      q2    = iq2*delq-delq
      f1    = foin(iq1)
      f2    = foin(iq2)
c
c     come closer to 'foin'=0.0 by repeated linear, interpolation.
c
      q3     = q1-f1*(q2-q1)/(f2-f1)
      do 20 itq0=1,itq0mx
      f3     = fbint(q3,foin)
      if(f1*f3.gt.zero) then
       q1     = q3
       f1     = f3
      else
       q2     = q3
       f2     = f3
      endif
      q3old  = q3
      q3     = q1-f1*(q2-q1)/(f2-f1)
      if(abs(q3-q3old).lt.endq0) goto 29
   20 continue
   29 continue
c
c     improve position of the zero by quadratic interpolation:
c
c                         the coefficients of the interpolating binomial
c
      f1     = fbint(q1,foin)
      f2     = fbint(q2,foin)
      f3     = fbint(q3,foin)
      ca     = one/((q1-q2)*(q2-q3)*(q3-q1))
      c11     = ca*(q3*(q3-q2)*q2*f1+q1*(q1-q3)*q3*f2+q2*(q2-q1)*q1*f3)
      c2     = ca*((q2*q2-q3*q3)*f1+(q3*q3-q1*q1)*f2+(q1*q1-q2*q2)*f3)
      c3     = ((q3-q2)*f1+(q1-q3)*f2+(q2-q1)*f3)*ca
c
c                         the two zeroes of the interpolating binomial
c
      qa1    = ( sqrt(abs(c2*c2-4.*c11*c3))-c2)/(c3+c3)
      qa2    = (-sqrt(abs(c2*c2-4.*c11*c3))-c2)/(c3+c3)
c
c                         select the zero which is closer to q3
c
      q      = qa1
      if(abs(qa1-q3).gt.abs(qa2-q3)) q=qa2
c
c     return zero thus found.
c
      qzero  = q
c
      return
      end
