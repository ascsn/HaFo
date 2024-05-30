      SUBROUTINE LN_iterate(e,delta,ph,gw,pairgw,partnm,nmax,ndim,iter,
     &                     iab,eferm,parnum,iwfoff)
c
c     input via list:
c      e       = array(1..nmax) of single particle energies.
c      delta   = array(1..nmax) of gap parameters.
c      gw      = array(1..nmax) of occupation probabilities,
c                input is used as initial guess for the iteration.
c      ph      = array(1..nmax) of the multiplicities of the single
c                particle levels.
c                This multiplicities INCLUDES here the spin degeneracy,
c                e.g. a 1s1/2 state has ph = 2 .
c      pairgw  = space for dynamical cut-off factor
c      partnm  = desired particle number
c      nmax    = maximum number of single particle states
c                (needs to be smaller than the parameter 'ndim')
c      ndim    = dimensioning of the arrays, particularly of 'vpairf',
c                in the calling routine
c      iter    = number of pairing iterations,
c                abs(iter) is the number of iterations
c      iab     = write level
c      iwfoff  = offset for levels and wavefunctions 
c
c      output via list:
c      eferm   = Fermi energy
c      parnum  = actual particle number after iteration
c
      INCLUDE 'params.inc'

      REAL(8),INTENT(IN) :: e(ndim),ph(ndim),delta(ndim)
      REAL(8),INTENT(IN OUT) :: pairgw(ndim),gw(ndim)

      REAL(8),INTENT(IN OUT) :: parnum,eferm
      INTEGER,INTENT(IN) :: iwfoff

      REAL(8) :: vln(ndim),suv2rs              ! workspaces for LN
C      data slpmin/2.0/          ! minimum slope in dN/dE
      REAL(8),PARAMETER :: xl2stp=1D0,dxllim=0.5d0,dx2lim=0.1D0
      REAL(8),PARAMETER :: detmin=1D-4,dampin=0.5D0,delefe=0.1D0
      INTEGER :: iftest=1,ifirst=0
      REAL(8),SAVE :: xlam2sav(2)=(/-1D99,-1D99/)
C     
c
c       the LN case where p.n. and p.n. fluctuations are iterated
c     
      IF(iwfoff>1) THEN
        IF(xlam2sav(2)<1D-98) THEN
          xlam2=0D0
        ELSE
          xlam2=xlam2sav(2)
        END IF
      ELSE
        IF(xlam2sav(1)<1D-98) THEN
          xlam2=0D0
        ELSE
          xlam2=xlam2sav(1)
        END IF
      END IF
      iftest  = iab   !  + 1
      iterh   = iter/2
      iterp   = iter/10
      do it=1,iter
        if(iftest.gt.0) write(6,'(a,i5,a)') ' iteration nr.',it,' :'
c     
c         initial call to pairing equations
c
        npart  = partnm+0.000001
        if(npart.eq.50 .and. ifirst.eq.1) then
          ifirst = 0
          suv2rs = MAX(SUM(gw(1:nmax)-gw(1:nmax)**2),1D-9)
          WRITE(6,'(a)') 'multiplicity:'
          WRITE(6,'(6(1pg12.4))') ph(1:nmax)
          WRITE(6,'(a)') 's.p. energies:'
          WRITE(6,'(6(1pg12.4))') e(1:nmax)
          WRITE(6,'(a)') 'gaps:'
          WRITE(6,'(6(1pg12.4))') delta(1:nmax)
          write(6,'(a)')  
     &      '  lambda     lambda_2            f1           f2    '
          do il2=0,20
            tlam2 = -0.5 + il2*0.05
            do il1=0,20
              tlam  = -25.0 + il1*1D0
                call lm2def(e,delta,ph,pairgw,1,nmax,iwfoff,npart,
     &                      tlam,tlam2,gw,gw,iftest,suv2rs, f1,f2)
              write(6,'(1x,2f8.2,3(1pg12.4))')
     &          tlam,tlam2,f1,f2,f1**2+5*f2**2
            enddo
          enddo
          stop
        endif
        suv2rs = MAX(SUM(gw(1:nmax)-gw(1:nmax)**2),1D-9)
        call lm2def(e,delta,ph,pairgw,1,nmax,iwfoff,npart,
     &        eferm,xlam2,gw,gw,iftest,suv2rs, f1,f2)
         if(iftest.gt.0)
     &        write(6,'(a,2g12.3,a,2g12.3)')
     &        ' xlam,xlam2=',eferm,xlam2,'   f1,f2=',f1,f2
c
c         test termination
c
         if(it.gt.1 .and. abs(f1).lt.eps 
     &        .and. abs(f2).lt.eps) goto 1099
c     
c     derivatives by finite differences
c
         call lm2def(e,delta,ph,pairgw,1,nmax,iwfoff,npart,
     &        eferm+delefe,xlam2,gw,vln,iftest,suv2rs, f1xl,f2xl)
         call lm2def(e,delta,ph,pairgw,1,nmax,iwfoff,npart,
     &        eferm,xlam2+xl2stp,gw,vln,iftest,suv2rs, f1xl2,f2xl2)
      if(iftest.ge.0) write(6,'(a,i4,6(1pg11.3))')
     &     ' LN: iter,xlam,xlam2,f1,f2=',it,eferm,xlam2,f1,f2
         dlf1   = (f1xl-f1)/delefe
         dlf2   = (f2xl-f2)/delefe
         dl2f1  = (f1xl2-f1)/xl2stp
         dl2f2  = (f2xl2-f2)/xl2stp
         if(iftest.ge.1) 
     &        write(6,'(a,5(1pg12.4))')
     &        ' dlf1,dlf2,dl2f1,dl2f2=',dlf1,dlf2,dl2f1,dl2f2,
     &        (dlf1*dl2f2-dlf2*dl2f1)
c
c         propagation to improved solution by inverse matrix
c
         if(it.lt.iterh) then
           det    = (dlf1*dl2f2-dlf2*dl2f1) !   - sign ?
           det    = sign(max(abs(det),detmin),det) ! lower limit for 'det'
           dxlam  = (f2*dl2f1-f1*dl2f2)/det
           dxlam2 = (f1*dlf2-f2*dlf1)/det
           if(iftest.ge.0) 
     &        write(6,'(a,4g12.4)')
     &        ' dir: dxlam,dxlam2=',dxlam,dxlam2
c
c       propagation by minimzation of least squares .
c       this as such is equivalent to direct gradient,
c       but here we can employ the "damping" by systematically
c       enhancing the diagonal elements (see Bevington step).
c
         else
           dchi1  = f1*dlf1 + 5.0*f2*dlf2
           dchi2  = f1*dl2f1 + 5.0*f2*dl2f2
           dchi11 = dlf1*dlf1 + 5.0*dlf2*dlf2
           dchi12 = dl2f1*dlf1 + 5.0*dl2f2*dlf2
           dchi22 = dl2f1*dl2f1 + 5.0*dl2f2*dl2f2
           inhanc = (it-iterh)/iterp+1
           enhanc = 2.0**inhanc
           dchi11 = dchi11*enhanc
           dchi22 = dchi22*enhanc
           det    = dchi11*dchi22-dchi12*dchi12
           det    = sign(max(abs(det),detmin),det) ! lower limit for 'det'
           dinv11 = dchi22/det
           dinv22 = dchi11/det
           dinv12 = -dchi12/det
           dxlam  = -(dinv11*dchi1+dinv12*dchi2)
           dxlam2 = -(dinv12*dchi1+dinv22*dchi2)
           if(iftest.ge.0) 
     &        write(6,'(a,4g12.4)')
     &        ' chi: dxlam,dxlam2=',dxlam,dxlam2
         endif
c
!         if(iftest.ge.1) 
!     &        write(6,'(a,4g12.4)')
!     &        ' dxlam,dxlam2=',dxlam,dxlam2
         quenf1 = dxllim/max(abs(dampin*dxlam ),dxllim)
         quenf2 = dx2lim/max(abs(dampin*dxlam2),dx2lim)
         quenf  = min(quenf1,quenf2)
         dxlam  = quenf*dampin*dxlam
         dxlam2 = quenf*dampin*dxlam2
cold         dxlam  = sign(min(abs(dampin*dxlam ),dxllim),dxlam )
cold         dxlam2 = sign(min(abs(dampin*dxlam2),dx2lim),dxlam2)
         eferm  = eferm + dxlam
         xlam2  = xlam2 + dxlam2
         if(iftest.ge.0)  write(6,'(a,6(1pg12.4))')
     &       ' dxlam etc=',dxlam,dxlam2,eferm,dxlam2,quenf1,quenf2
      enddo
c     
      if(iftest.ge.0) write(6,'(a)')
     &   ' no convergence of pairing parameters in PAIRDN, LN-branch'
c
 1099 continue
      parnum = f1+npart
      if(iftest.ge.0) write(6,'(a,i4,6(1pg11.3))')
     &     ' PAIRDN-LN: iter,xlam,xlam2,f1,f2=',it,eferm,xlam2,f1,f2
      IF(iwfoff>1) THEN
        xlam2sav(2)=xlam2
      ELSE
        xlam2sav(1)=xlam2
      END IF

      RETURN
      END


c-----lm2def-------------------------------------------------------------------
c
      subroutine lm2def(spen,delta,phm,pairgw,nspd,nspg,iwfoff,npart,
     &                  xlam,xlam2,vv2in,vv2out,iftest,suv2rs, f1,f2)
c
c---------------------------------------------------------------------
c
      include 'params.inc'
      parameter (kmaxln=KSL2)
c
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'grid.inc'
c
c---------------------------------------------------------------------
c
c     routine to compute the second order multiplies 'lambda_2`
c     for the LN scheme in connection with the delta-force pairing.
c
c     input:
c       spen   = array of s.p. energies
c       delta  = array of pairing gaps
c       phm    = array of multiplicities, e.g. a j=7/2 state has phm = 4 .
c       pairgw = cutoff weight complementing gap 'delta'
c       nspd   = lowest level in pairing
c       nspg   = highest level in pairing
c       iwfoff = offset for addressing wavefunctions,
c                note that wf's require appropriate
c                absolute addresses because information
c                on wavefunctions is entered via common /wavefr/ .
c       npart  = particle number
c       vpair  = strength of pairing force
c       xlam   = Fermi energy
c       xlam2  = initial value for fluctuation constraint
c       vv2in  = array of initial internal occupation probabilities
c       iftest = switch for testprints
c       ilipng = switch between BCS and Lipkin-Nogami,
c                 1 --> f3=0
c     the following two parameters are entered via common
c     whereby proton <--> neutron is switched by 'iwfoff'
c       erspln = Skyrme response energy for 'lambda_2'
c       suv2ra =   u^2*v^2 from response routine, for normalization
c       vpair  = pairing strength
c     output:
c       f1     = sum(2v**2)-n         =  particle number condition
c       f2     = xlam2-vpair/4*(sum(u**3*v)sum(v**3*u)-
c                               sum(u**4*v**4))/((sum(u**2*v**2))**2 -
c                               sum(u**4*v**4))
c              = LN condition
c       vv2out = new occupation weights
c
c---------------------------------------------------------------------
c
c
      REAL(8),INTENT(IN) :: spen(kmaxln),vv2in(kmaxln),phm(kmaxln)
      REAL(8),INTENT(OUT) :: vv2out(kmaxln)
      REAL(8),INTENT(IN) :: delta(kmaxln),suv2rs
      REAL(8),INTENT(IN OUT) :: pairgw(kmaxln)      
      INTEGER,INTENT(IN) :: iftest
c
c       iterv2 = max.nr. iterations in 'vv2'
c       epsv2  = precision on 'vv2'
c       epsdv2 = limit on change of energy in determining 'vv2'
c       dampin = general reduction of the step, top stabilize
c
      INTEGER,PARAMETER :: iterv2=50
      REAL(8),PARAMETER :: epsv2=1D-7,epsdv2=1D-7,dampin=0.5D0
      LOGICAL :: tfirst=.TRUE.,texis
      REAL(8),SAVE :: gpairfac=1D0
c
c---------------------------------------------------------------------
c
c     preparation
c
      IF(tfirst) THEN
        INQUIRE(FILE='LNparams',EXIST=texis)
        IF(texis) THEN
          OPEN(112,file='LNparams')
          READ(112,*) gpairfac
        END IF
        tfirst=.FALSE.
      END IF
      if(nspg.gt.kmaxln) stop ' last level > "kmaxln"'
      fnn    = npart
      if(iwfoff.eq.0) then
        erspln = rspln(1)
        suv2ra = suv2rs
        vpair  = pairfp
        gpair  = gpp
      else
        erspln = rspln(2)
        suv2ra = suv2rs
        vpair  = pairfn
        gpair  = gpn
      endif
      if(iftest.ge.3) 
     &  write(6,'(a,i5,3g12.4)')
     &   ' iwfoff,vpair,erspln,suv2ra=',iwfoff,vpair,erspln,suv2ra
c
c     reset accumulators for modified pairing currents
c
      aux1(1:ngrid)=0D0
      aux2(1:ngrid)=0D0
c
c     reset accumulators for combinations of 'u' and 'v'.
c
      sv2   = 0D0                                  ! sum v^2
      suv   = 0D0                                  ! sum u*v
      suv4  = 0D0                                  ! sum u^4 v^4
      suv2  = 0D0                                  ! sum u^2 v^2
      su3v  = 0D0                                  ! sum u*v u^2
      suv3  = 0D0                                  ! sum u*v v^2
c
      alpha  = 4D0*xlam2
c
      do k=nspd,nspg
c
c     iteration loop to compute occupation 'vv2'
c
        d2     = (delta(k)*pairgw(k))**2
        alphd2 = alpha*0.5D0*d2
        e     = spen(k)-xlam
        v2    = vv2in(k)                                ! start from previous value
        x     = e+alpha*v2
        if(iftest.ge.3)
     &   write(6,'(a,i5,a,3(1pg12.4))')
     &    ' state nr.',k,' :  v2in,e,delta=',v2,e,delta(k)*pairgw(k)
c
c                            iterate to s.c. occup.numbers and lambda_2
c       
!        WRITE(*,*) 'alpha etc:',alpha,alphd2
        do it=1,iterv2
!          WRITE(6,*) ' it etc:',it,x,d2,derf
          eq    = sqrt(d2+x*x)
          v2    = 0.5D0-0.5D0*x/eq
          f     = e+alpha*v2-x
          derf  = -alphd2/(eq*eq*eq)-1D0
          dx    = dampin*f/derf         ! optionally switch on if unstable
!          dx    = MAX(MIN(dx,1D0),-1D0)
CTEST          dx    = f/derf
          if(abs(f).lt.epsv2 .and. abs(dx).lt.epsdv2) go to 7
          x     = x-dx
        enddo
c
ctst        write(6,'(a)')
ctst     &    ' maximal number of iterations exceeded in "LM2DEF"'
        if(iftest.ge.3) write(7,'(a)')
     &    ' maximal number of iterations exceeded in "LM2DEF"'
c
    7   continue
        v2   = max(sqrt(1D-20),v2)
        vv2out(k)=v2
c
c       accumulate sums with new 'v2'
c
        xmult = phm(k)*2D0
        if(iftest.ge.3 .AND. it==1)
     &   write(7,'(a,4i4,4(1pg12.4))')
     &    ' nr,state,it,iwf,v2,ph,',npart,k,it,iwfoff,v2,xmult
        u2    = 1D0-v2
!        u2    = max(sqrt(1D-20),u2)
!        WRITE(*,*) 'uv etc:',k,u2,v2,u2*v2
        uv    = sqrt(MAX(u2*v2,1D-20))
        sv2   = sv2  + xmult*v2
        suv   = suv  + xmult*uv
        suv4  = suv4 + xmult*u2*u2*v2*v2
        suv2  = suv2 + xmult*u2*v2
        su3v  = su3v + xmult*uv*u2
        suv3  = suv3 + xmult*uv*v2
c
c       accumulate modified pair currents
c
!        istate = k+iwfoff
        nadd = (k+iwfoff-1)*ngrid
!        WRITE(*,*) 'k,iwfoff,ngrid,nadd=',k,iwfoff,ngrid,nadd
        aux1(1:ngrid) = aux1(1:ngrid) + xmult*
     &               (uv*u2)*(wf(nadd+1:nadd+ngrid)*rinvs(1:ngrid))**2
        aux2(1:ngrid) = aux2(1:ngrid) + xmult*
     &               (uv*v2)*(wf(nadd+1:nadd+ngrid)*rinvs(1:ngrid))**2
!        call sdenti(psi(0,0,1,istate),wrk1,uv*u2,istate)
!        call sdenti(psi(0,0,1,istate),wrk2,uv*v2,istate)
c
      enddo
c
c     first result: particle number 
c
      f1     = sv2-fnn          
!      f1     = 2D0*sv2-fnn          
c
c     second result: the LN equation
c     
c     but the numerator from Skyrme response is a bit special.
c     it was computed with an older "v2". The corresponding
c     "sum V2*U2" is stored in "suv2ra". In order to avoid
c     all to frequent recomputation of that expensive piece,
c     we rescale it here just by the new "suv2".
c
      denom  = 32.0*(suv2*suv2-suv4)                   !  denominator
!     IF(tpairFy) THEN
      IF(ioptsk==13) THEN
        pint = 2D0*fy_pairint(aux1,aux2,iwfoff)
!        WRITE(*,*) 'pint=',pint
      ELSE IF(ipair==5 .OR. ipair==9) THEN
        pint = SUM(aux1(1:ngrid)*aux2(1:ngrid)*rsquar(1:ngrid))*frpirs
        pint = -vpair*2.0*pint
      ELSE IF(ipair==7 .OR. ipair==10) THEN
        pint = SUM(aux1(1:ngrid)*aux2(1:ngrid)*rsquar(1:ngrid)*
     &          (1D0-(rho(1:ngrid)+rho(ngridp:ngrid2))/rho0pr))*frpirs
        pint = -vpair*2.0*pint
      ELSE
        STOP 'LM2DEF: wrong option for IPAIR'
      END IF
!       write(6,'(a,6g12.4)') ' suv2 etc=',suv2,suv4,suv2ra,erspln,vpair,pint
!     xlam2 is augmented by still flexible re-normalization 'gpairfac'
!     xlam2 from Jacek's recipe with effective pairing strength
!      enume1 = -gpair/(gpairfac*4D0)*(su3v*suv3-suv4) !  1. numerator - average
!     xlam2 from explicit pairing functional (Skyrme only)
      enume1 = pint/(256D0*gpairfac)       !  1. numerator - detail
      enume2 = suv2*suv2*erspln/(suv2ra*suv2ra)        !  2. numerator
!      enume2=0D0     ! override self-consistent fead back
      f2     = xlam2+(enume1+enume2)/denom             ! p.n. fluctuation eq.
c
      if(iftest.ge.3) 
     &  write(6,'(a,6g12.4)')
     &   ' suv2,suv4,suv2ra,erspln,vpair=',
     &   suv2,suv4,suv2ra,erspln,vpair
c
      return
      end
