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
       character*72 swbinp,shellout,shellres
       character*12 title,titdum,titsav
       dimension cg(2000)
       dimension enipr(25000),enine(25000),eni0pr(25000),eni0ne(25000)

       common/bcs/ eni(25000)/bcs0/eni0(25000)/old/iold
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
c
c
       namelist /strutek/ swbinp,shellout,shellres,iz,in,
     *                    gampro,gamneu,ippro,ipneu,
     *                    ecut,enei,enef,dene,ilev,iold,iendd
c
       enei=-50.d0
       enef=30.d0
       dene=0.1d0
       ecut=60.2
       ilev=0
       iold=0
       ippro=10
       ipneu=10
       gampro=1.54
       gamneu=1.66
       iendd=0
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
       call swbrea(enipr,ngorp,nf1,title)  ! read protons
       call swbrea(enine,ngorn,nf1,title)  ! read neutrons
       call swbrea(eni0pr,ngor0p,nf1,titdum)  ! read proton gas 
       call swbrea(eni0ne,ngor0n,nf1,titdum)  ! read neutron gas 
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
         enddo
       else
         gam   = gampro
         ngor  = ngorp
         ip    = ipneu
         do i=1,ngor
           eni(i)  = enipr(i)
           eni0(i) = eni0pr(i)
         enddo
       endif
       gamm=gam      
       write(6,'(a,i3)') ' call SHCORR: iprot=',iprot
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
      subroutine swbrea(ener,ilelev,nr,title)
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
      dimension  ener(1)
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
c      write(6,101)j2,l,energ
      ind=ind+1
      ener(ind)=energ
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
       common/bcs/ eni(25000)/bcs0/eni0(25000)/old/iold
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
      write(6,'(a,i3)') ' leave SHCORR: iprot=',iprot
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

