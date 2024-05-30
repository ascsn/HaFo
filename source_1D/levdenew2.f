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
c
c      written by w. nazarewicz, 1985, 1994                         c
c  modified by A. T. Kruppa and T. Vertse 1999 in order to include  c
c  continuum effect  and generalized smoothing procedure            c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit double precision (a-h,o-z)
       parameter(maxlev=1000,ndmain=500)
       character*72 swbinp,shellout,shellres,shellsmooth
       character*20 title,titdum,titsav
       dimension cg(2000),tab(20),dtab(20),sig(ndmain),ilf(20),
     *           esm(ndmain),alam(ndmain)
       dimension enipr(25000),enine(25000),eni0pr(25000),eni0ne(25000)
       dimension ps(25000),ps0(25000)
    
       common/strut/ her(40),dher(40),cm(40),gammaq,
     *               ip,ilev,inucl
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
c
c         ecut      -   cut-off energy above lambda (in hbar*omega)
c
c         iendd     -   0 (continue), 1 (end)
c         ipart     -   1  calculate smoothed particle number
c         idens     -   1  calculate smoothed densities
c
c
       namelist /strutek/ swbinp,shellout,shellres,shellsmooth,iz,in,
     *                    gam0,gam1,dgam,ip0,ip1,alg,felg,
     *                    ecut,enei,enef,dene,ilev,iold,iendd,ipart

       enei=-80.d0
       enef=30.d0
       dene=0.1d0
       ecut=60.2
c       fel=0.75d0
c       sav=1.5d0
       fel=2.d0
       sav=4.d0
       ilev=0
       iold=0
       iendd=0
       ipart = 0
       shellout = 'shell.out'
       shellres = 'shell.res'
       shellsmooth = 'shell.smooth'
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
	 if(ipart.eq.1) 
     &      open(9,file=shellsmooth,form='formatted',status='unknown')
         write(nf2,'(a,i4,a,i4)')
     &  ' shell corrections for effective mean field:'
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
       call swbrea(enipr,ngorp,nf1,title)  ! read protons
c       write(*,*)'after reading protons'
       call swbrea(enine,ngorn,nf1,title)  ! read neutrons
c       write(*,*)'after reading neutrons'
       write(nf2,'(a)')title
       call swbrea(eni0pr,ngor0p,nf1,titdum)  ! read proton gas 
       call swbrea(eni0ne,ngor0n,nf1,titdum)  ! read neutron gas 

c
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
c
       do iprot=1,0,-1
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
c        write(*,*)' alg,felg=',alg,felg
c	write(*,*)' iel,ief',iel,ief
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
9876  continue
        write(2,*)
       endif
        if(ipart.eq.1)then
c  calculate smoothed particle numbers for each orbits	
	  call partsm
     &    (xlamb,eni,eni0,her,dher,ip,cm,ngor,inucl,gm,ps,ps0)
          sumnpart=0.d0
          sumtest = 0.0D0
	  do i=1,ngor
            sumnpart=sumnpart+ps(i)-ps0(i)
            sumtest = sumtest+ps(i)
          enddo
          write(*,*)' sumnpart,sumtest=',sumnpart,sumtest
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
      if(iprot.EQ.1) then
        write(9,'(a)') ' smoothed occupations for protons:'
      else
        write(9,'(a)') ' smoothed occupations for neutrons:'
      endif
      multi = 0
      do i=1,ngor
        if(ensav.NE.eni(i)) then
          multi=1+multi
          iniv = 1+iniv
          write(9,'(2i5,2(1pg13.5))')  iniv,multi,eni(i),ps(i)/2.0D0
          multi = 0
          ensav = eni(i)
        else
          multi = 1+ multi
        endif
      enddo

      enddo
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
      function delt(xk,yk,n,nagyn,c)
      implicit real*8 (a-h,o-z)
c
c    subrutine for calculating delta value
c
      dimension xk(nagyn),yk(nagyn),c(0:n)
      delt=0.d0
      do i=1,nagyn
      delt=delt+(yk(i)-horn(n,c,xk(i)))**2
      enddo
      return
      end

      function horn(n,a,x)
c
c function for the Horner schem
c  calculates value of the polinomial of degree n
c  horn(x)=a(0)+a(1)*x+...+a(n)*x**n
c
      implicit real*8 (a-h,x-z)
      dimension a(0:n)
      b=a(n)
      do i=n,1,-1
      b=a(i-1)+x*b
      enddo
      horn=b
      return
      end
      subroutine normal(xk,yk,n,nagyn,c)
      parameter (nmax=49)
      implicit real*8 (a-h,o-z)
c
c    subrutine for solution of the normal equations
c
      dimension xk(nagyn),yk(nagyn),c(0:n)
      dimension a(0:nmax,0:nmax),b(0:nmax)
      do j=0,n
        b(j)=0.d0
        do i=j,n
          s=0.d0
          ij=i+j
          do k=1,nagyn
            if(ij.eq.0)then
              s=s+1.d0
              b(j)=b(j)+yk(k)
            else 
              s=s+xk(k)**ij
              if(i.eq.j) b(j)=b(j)+yk(k)*xk(k)**j
            endif
          enddo
          a(i,j)=s
          if(i.ne.j) a(j,i)=s
        enddo
      enddo
      np1=n+1
      nmaxp=nmax+1
c
c     solution of the system of lin. equations with np1 unknowns
c     a*x=b
c
      call gaussj(a,np1,nmaxp,b,1,1)
c     solution x is in array b at exit 
      do i=0,n
      c(i)=b(i)
      enddo
      return
      end
      SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)
      implicit real*8 (a-h,o-z)
c
c     solution of the system of lin. equations with n unknowns
c     a*x=b
c     using Gauss-Jordan elimination with full pivoting.
c     input quantities: matrix a of the coefficients and matrix b
c      which containes the (multiple) right hand side of the system
c     at the exit it will be overwritten by the solution i.e.
c     solution x is in array b at exit 
c
      PARAMETER (NMAX=50)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                PAUSE 'Singular matrix'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        IF (A(ICOL,ICOL).EQ.0.) PAUSE 'Singular matrix.'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
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
      character*4 aline
      character*14 title
*
      ilelev=0
*
      do iread=1,9999
c      write(*,*)' before'
        read(nr,'(a)',err=98,end=98) aline
c      write(*,*)' after'
c       write(*,9123)aline(1:4)
9123	format((a))
        if(aline(1:4).eq.'  2j') goto 999
        if(aline(1:4).eq.' * *') then
          read(nr,'(a)',err=98,end=98) title
c	  write(*,'(a)') title
        endif
      enddo
 999  continue
*
c        write(*,*)' before'
      ind=0
cold 9    read(nr,101)j2,l,energ
c 9    read(nr,123)j2,l,energ,xn,imult
 9    read(nr,*)j2,l,energ,xn,rms,imult
c        write(*,*)' after'
123   format(2i5,f9.4,e17.5,i5)
c      write(*,*)j2,l,energ,xn,rms,imult
      ind=ind+1
      ener(ind)=energ
      if (j2.lt.0) goto 19
      goto 9
 19   continue
      ilelev=ind
c      write(*,*)'ind ',ind
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

