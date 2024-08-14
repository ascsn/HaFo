c-----pairol------------ part of the hartree-fock package ----------------------
c
      subroutine pairol(e,gw,ph,nz,nmax,gp,eferm,epair,delta,ipair,
     &                  partnm,iter,eps,iab)
c
c     this is a pairing pairing routine. but it is presently used
c     as a routine to determine the thermal occupation of the
c     single particle states. this is achieved in the "constant gap"
c     option, using 'iter<0' as a switch.
c
c     the include file 'params.inc' communicates the dimensioning of
c     the fields 'e' and 'gw' with the parameter 'ksl2'. it allows
c     furthermore to make the whole subroutine double precision by
c     transfering a line 'implicit double precision (a-h,o-z)'.
c
      include 'params.inc'
c
c     * * * * * * * * * *     pairing iteration      * * * * * * * * * * * *
c
c     determines the fermi-energy 'eferm' and the gap 'delta'
c     such that the "gap-equation" and the "particle number condition"
c     are fulfilled simultanously. this is achieved by iteration
c     with a gradient scheme (newton's tangential formula)
c     where the gradient is evaluated by finite differences.
c
c     there are two schemes possible:
c     1. the "constant force" approach, where the pairing force is given
c                                       in 'gp'.
c     2. the "constant gap" approach,   where the gap is kept fixed and
c                                       is given in 'delta'.
c     there is only 'eferm' to be determined in the "constant gap"
c     approach. the gap-equation is used at the end to determine an
c     according force 'gp'.
c
c     input via list:
c      e       = array(1..nmax) of single particle energies.
c      gw      = array(1..nmax) of occupation probabilities,
c                input is used as initila guess for the iteration.
c      ph      = array(1..nmax) of the multiplicities of the single
c                particle levels, to be provided by the calling routine.
c      nz      = desired particle number
c      nmax    = maximum number of single particle states
c                (needs to be smaller than the parameter 'ksl2'
c                 and smaller than 46, for dimensioning reasons).
c      gp      = size of pairing-force, is input in case of
c                "constant force" approach.
c      eferm   = fermi-energy, input is used as initial guess;
c                if eferm=0.0 is input then a suitable initial guess
c                is computed from input 'gw'.
c      epair   = pairing energy
c      delta   = pairing gap, is input in case of "constant gap" approach.
c      ipair   = switch for type of pairing
c                 0 --> constant Delta
c                 1 --> constant Force
c                 4 --> Fermi distribuion using T=Delta
c      iter    = number of pairing iterations
c      eps     = end condition for pairing iteration;
c                iteration terminates if error in particle number
c                is smaller than 'eps' and if 'eferm' and/or 'delta'
c                are smaller than 'eps'.
c      iab     = print level; significant values are -2, -1, 0 and 1;
c                larger value prints more.
c
c     output via list:
c      gw      = array(1..nmax) of the occupation probabilities
c                finally achieved.
c      gp      = size of pairing-force, is output in case of
c                "constant gap" approach.
c      delta   = pairing gap, is output in case of
c                "constant force" approach.
c      eferm   = fermi-energy.
c      sum     = particle number, from adding 'gw' with degeneracy 'ph'.
c
c     some auxiliary fields are:
c      dparn   = array(1..3) of variation in the particle number as
c                consequence of 'deferm' and 'ddelta'.
c      dgap    = array(1..3) of variation in the gap equation as
c                consequence of 'deferm' and 'ddelta'.
c      gr      = double array of variation of particle number and of
c                gap equation, used for construction of the tangantial
c                step
c
      dimension e(ksl2),gw(ksl2),ph(ksl2)
      dimension dparn(3),dgap(3),gr(2,2)
c
cold      data ph/1.,2.,1.,3.,2.,1.,4.,3.,2.,1.,5.,4.,3.,2.,1.,6.,5.,
cold     &        4.,3.,2.,1.,7.,6.,5.,4.,3.,2.,1.,8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        9.,8.,7.,6.,5.,4.,3.,2.,1.,10./
      data gr,dparn,dgap/10*0./
c
c     further numerical variables are:
c      deferm  = step in the fermi-energy, to evaluate
c                numerically the derivatives.
c      ddelta  = step in the gap, to evaluate numerically the
c                derivatives; used only in the "constant force" approach.
c      dampin  = damping factor for the tangential step
c      efmax,dlmax,detmin,defmin,ddlmin  = limitations for the step
c
c     note: these parameters are adjusted for an application with
c           the nuclear shell model. it is assumed that the single
c           particle energies are given in units of [mev] and display
c           a reasonable nuclear spectrum.
c
      data dampin,deferm,ddelta/0.9,1.0,0.1/
      data efmax,dlmax,detmin,defmin,ddlmin/3.,.6,.1e-3,.1e-3,.1e-3/
c
      data xmaxlw,x10lw/1.0e-30,1.0e-10/
      data one,zero/1.0,0.0/
c
c-------------------------------------------------------------------------------
c
c      set switch variables
c
      iwarng = 0
c
      if(iab.ge.0)
     & write(6,'(/a,i4,4(a,f9.4)/a,10(/1x,10f8.3))')
     &   ' *** pair entered with: iter=',iter,'  delta=',delta,
     &   '  eferm=',eferm,'  gp=',gp,'  deferm=',deferm,
     &   '                        single particle energies:',
     &   (e(i),i=1,nmax)
c
c      prepare fermi-energy 'eferm', if not yet initialised.
c
      if(eferm.eq.zero) then
       isum   = 0
       do 120 it=1,nmax
       itsum  = it
       isum   = int(ph(it)+.000001) +isum
       if(isum+isum-nz.ge.0) goto 110
  120  continue
       stop ' not enough states too reach particle number'
  110  continue
       eferm  = e(itsum)+deferm
       if(iab.ge.1) write(06,'(5x,a,i5,a,g11.4)')
     &   'initialistion: isum=',isum,'  eferm=',eferm
c
      endif
c
      if(ipair.eq.1) then
c
c-------------------------------------------------------------------------------
c
c     case of constant pairing-force, given by 'gp'.
c      gap 'delta', fermi-energy 'eferm', and occupations 'gw'
c      are computed.
c
c-------------------------------------------------------------------------------
c
       if(ipair.ne.0 .and. ipair.ne.1 .and. ipair.ne.4)
     &   stop ' invalid option "ipair" in subroutine "pairol"'
       delti  = delta
       if(abs(delti).lt.x10lw) delti = one
       do 200 it=1,abs(iter)
c
c       gap-equation & particle-number at 'eferm' & 'delti' (for k=1)
c       and at 'eferm' varied (k=2) and at 'delta' varied (k=3).
c
        do 210 ka=1,3
c                                          actual fermi-energy and gap
        k      = 4-ka
        if(k.eq.1) then
         elam   = eferm
         delact = delti
        elseif(k.eq.2) then
         elam   = eferm+deferm
         delact = delti
        else
         elam   = eferm
         delact = delti+ddelta
        endif
c                                         compute actual quasiparticle energies
c                                         and occupation-weights,
c                                         accumulate gap-equation and
c                                         particle number
        parnum = zero
        gapeq  = zero
         do 220 i=1,nmax
         equasi = sqrt((e(i)-elam)*(e(i)-elam)+delact*delact)
         gw(i)  = max(0.5-0.5*(e(i)-elam)/equasi,zero)
         gapeq  = ph(i)/equasi   + gapeq
         parnum = gw(i)*ph(i)    + parnum
  220    continue
c                                         store actual part.number and
c                                         gap-equation
        dgap(k)  = 0.5*gp*gapeq-one
        dparn(k) = parnum+parnum-nz
  210   continue
c
c       construction of the gradient-matrix and its inverse
c        ( determinant 'det' is limited in value to avoid singularities )
c
        gr(1,1) = ( dgap(2)- dgap(1))/deferm
        gr(1,2) = (dparn(2)-dparn(1))/deferm
        gr(2,1) = ( dgap(3)- dgap(1))/ddelta
        gr(2,2) = (dparn(3)-dparn(1))/ddelta
        det     = gr(1,1)*gr(2,2)-gr(1,2)*gr(2,1)
        det     = sign(max(abs(det),detmin),det)
        if(abs(det).lt.xmaxlw)  det    = sign(xmaxlw,det)
c
        if(iab.ge.1) then
         write(6,'(a,i4,3(a,g11.3)/a,3g11.3,a,3g11.3/a,4g11.3)')
     &    ' iteration at it=',it,' :  delti=',delti,'  eferm=',eferm,
     &    '  det=',det,
     &    ' dgap=',dgap,' ,  dparn=',dparn,
     &    '   gr=',gr
         write(6,'(a,10f7.4,10(/6x,10f7.4))')
     &    '   gw=',gw
        endif
c
        grsav   =  gr(1,1)
        gr(1,1) =  gr(2,2)/det
        gr(1,2) = -gr(1,2)/det
        gr(2,1) = -gr(2,1)/det
        gr(2,2) = grsav/det
c
c       tangential step  (including some step limiters)
c
        efold  = eferm
        delold = delti
        efstep = dampin*(gr(1,1)*dgap(1)+gr(2,1)*dparn(1))
        dlstep = dampin*(gr(1,2)*dgap(1)+gr(2,2)*dparn(1))
        eferm  = eferm-sign(min(abs(efstep),efmax),efstep)
        delti  = abs(delti-sign(min(abs(dlstep),dlmax),dlstep) )
        deferm = efold-eferm
        ddelta = delold-delti
        if(abs(deferm).lt.defmin) deferm = defmin
        if(abs(ddelta).lt.ddlmin) ddelta = ddlmin
c
        if(iab.eq.0)
     &   write(6,'(a,i4,2(a,g11.3))')
     &    ' pairing iteration it=',it,'  :  delti=',delti,
     &    '  eferm=',eferm
c
        if(abs(dgap(1)).lt.eps .and. abs(dparn(1)).lt.eps)   goto 299
  200   continue
        iwarng = 1
c
  299   continue
        delta  = delti
c
      else
c
c-------------------------------------------------------------------------------
c
c      case of constant paring-gap, given by 'delta'.   
c      or case of Temperature, given by Delta
c
c-------------------------------------------------------------------------------
c
       do 310 it=1,abs(iter)
c
c      particle number at 'eferm' (k=1) and 'eferm' varied (k=2).
c
        do 320 ka=1,2
c                                         set actual fermi-energy 'elam'
        k      = 3-ka
        if(k.eq.2) then
         elam   = eferm+deferm
        else
         elam   = eferm
        endif
c                                         compute actual particle number
        parnum = zero
	if(ipair.eq.0) then                          ! pairing
         do 330 i=1,nmax
         equasi = sqrt( (e(i)-elam)*(e(i)-elam) + delta*delta )
         gw(i)  = 0.5-0.5*(e(i)-elam)/equasi
         parnum = gw(i)*ph(i)    + parnum
  330    continue
        elseif(ipair.eq.4) then                      ! thermal
         do 331 i=1,nmax
         equasi = sqrt( (e(i)-elam)*(e(i)-elam) + delta*delta )
         gw(i)  = 1.0/(1.0+exp((e(i)-elam)/delta))
         parnum = gw(i)*ph(i)    + parnum
  331    continue
	endif
c                                         store particle number
        dparn(k) = parnum+parnum-nz
  320   continue
c
c      gradient step to next 'eferm'   (including step limiters)
c
        gr(1,1) = (dparn(2)-dparn(1))/deferm
        det    = sign(max(abs(gr(1,1)),detmin),gr(1,1))
        efold  = eferm
        efstep = dampin*dparn(1)/det
        eferm  = eferm-sign(min(abs(efstep),efmax),efstep)
        deferm = efold-eferm
        deferm = sign(max(abs(deferm),defmin),deferm)
        if(iab.eq.0)
     &   write(6,'(a,i4,2(a,g11.3))')
     &    ' pairing iteration it=',it,'  :  eferm=',eferm,
     &    '  dparn(1)=',dparn(1)
        if(iab.ge.1)
     &   write(6,'(a,i4,2(a,g11.3)/a,2g11.3,a,g11.3/
     &             a,10f7.4,10(/6x,10f7.4))')
     &    ' iteration at it=',it,' :   eferm=',eferm,
     &    '  det=',det,
     &    ' dparn=',dparn(1),dparn(2),' ;   gr=',gr(1,1),
     &    '   gw=',gw
        if(abs(dparn(1)).lt.eps)   goto 350
  310   continue
        iwarng = 1
  350  continue
c
c      determine force 'gp' according to given 'delta' and the
c      equilibrium reached.
c
       gapeq  = zero
        do 360 i=1,nmax
  360   gapeq = ph(i)/sqrt((e(i)-eferm)*(e(i)-eferm)+delta*delta)
     *           + gapeq
       gp     = 2.0/gapeq
c
      endif
c
c-------------------------------------------------------------------------------
c
c     finally the gap-equation and particle-number are checked and printed
c
c-------------------------------------------------------------------------------
c
      if(ipair.eq.0 .or. ipair.eq.1) then
       gphalf = gp*0.5
       del2   = max(xmaxlw,delta*delta)
       gapeq  = zero
       partnm = zero
        do 410 i=1,nmax
        partnm = gw(i)*ph(i) + partnm
  410   gapeq  = gphalf*ph(i)/sqrt(del2+(e(i)-eferm)**2) + gapeq
c
        if(iab.ge.-1)
     &   write(6,'(a,i4,2(a,g13.5))')
     &   ' pair finished with ',it,' iterations: gap-eq.=',gapeq,
     &   '  part.nr=',partnm
c                                  the pairing energy
        epair    = -delta*delta/(gp+xmaxlw)
      elseif(ipair.eq.4) then
        entrop = zero
        partnm = zero
        do 411 i=1,nmax
        partnm = gw(i)*ph(i) + partnm
        occ    = max(gw(i),x10lw)
        occd   = max(1.0-occ,x10lw)
        entrop = 2.0*ph(i)*(occ*log(occ)+occd*log(occd)) + entrop
  411   continue 
c
        if(iab.ge.-1)
     &   write(6,'(a,i4,2(a,g13.5))')
     &   ' pair finished with ',it,' iterations: entropy=',entrop,
     &   '  part.nr=',partnm
c                                 the T*S
	epair   = Delta*entrop
      endif
      IF(iab.GE.0) WRITE(6,'(a,50(/5(1pg13.5)))') ' weights:',gw(1:nmax)
c
c     in case of incomplete convergence print a warning
c
      if(iab.ge.-1 .and. iwarng.eq.1)
     & write(6,'(2(a,g11.3),2(a,i3))')
     &  ' ---> no convergence in pair! err(gap-eq)=',dgap(1),
     &  '  err(part.nr)=',dparn(1),'  nz=',nz,'  nmzx=',nmax
c
      return
      end

c========================================================================
c========================================================================
c    
c     Here follows the Pairing and Lipkin-Nogami package of
c      -->   G.A.Leander and W.Nazarewicz, march 1985   <--
c     revised by P.-G. Reinhard, june 1994
c
c     BCS is standard, but includes correction 'gpair*v^2' to
c     the s.p. energies.
c
c     The method of Lipkin and Nogami and coworkers is explained,
c     e.g. in Nucl. Phys. A201 (1973) 357 .
c
c========================================================================
c-----bcslip-------------------------------------------------------------------
c
      subroutine bcslip(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,
     &                  vv2,vln,epair,ecorr,entot,erspln,suv2rs,
     &                  ilipng,iftest)
c
ctest      implicit double precision (a-h,o-z)
c
ctest      parameter (kmaxln=46)
      include 'params.inc'                       ! activate that for
      parameter (kmaxln=ksl2)                    ! HaFo package
c
c------------------------------------------------------------------------------
c
c     Program to solve the Lipkin-Nogami, or BCS, equations
c       gap equation:      0 = f1 = delta-gpair*sum(u*v)
c       particle number:   0 = f2 = sum(2v**2)-n
c       p.n. fluctuation:  0 = f3 = xlam2-gpair/4*(sum(u**3*v)sum(v**3*u)-
c                                   sum(u**4*v**4))/((sum(u**2*v**2))**2 -
c                                   sum(u**4*v**4))
c     the f3 is omitted for BCS.
c
c     The 3 equations are solved by Newton's method.
c     The matrix of derivatives is computed by finite differences.
c     The internal numerical parameters are explained below.
c
c------------------------------------------------------------------------------
c
c     The input parameters are:
c
c       spen   = array of s.p. energies
c       ph     = array of multiplicities, e.g. a j=7/2 state has ph = 4 .
c       nspd   = lowest level in pairing
c       nspg   = highest level in pairing
c       npart  = particle number
c       gpair  = strength of pairing force
c       delta  = initial value for pairing gap
c       xlam   = initial value for Fermi energy
c       xlam2  = initial value for fluctuation constraint
c       vv2    = array of internal occupation probabilities (initial)
c       erspln = Skyrme response energy for 'lambda_2', needed for
c                Skyrme H.F. codes, to be provided in a separte routine.
c                In all other cases set to 0.0 .
c       suv2rs = u^2*v^2 from Skyrme response routine, for normalization.
c                In all other cases set to 1.0 .
c       ilipng = switch between BCS and Lipkin-Nogami
c                 1 --> Lipkin-Nogami
c                 0 --> BCS
c       iftest = switch for testprints
c
c     The output parameters are:
c       delta  = result for pairing gap
c       xlam   = result for Fermi energy
c       xlam2  = result for fluctuation constraint
c       vv2    = array of internal occupation probabilities
c       vlb    = array of occupation probabilities for external observables
c       epair  = naive pairing energy , -delta^2/gpair
c       ecorr  = consistent pairing energy as correction for the hf energy:
c                =-delta**2/gpair
c                 -gpair*sum{v(i)**4-v(i)**2}
c                 -4*xlam2*sum_(u(i)*v(i))**2
c                = -(4*xlam2-gpair)*sum_(u(i)*v(i))**2
c       entot  = total s.q.p. energy (= s.p. sum - g * npart/2 for normal state)
c
c------------------------------------------------------------------------------
c
c
      dimension spen(kmaxln),vv2(kmaxln),vln(kmaxln),ph(kmaxln)
c
ctest      dimension derv(3,3)
c
c------------------------------------------------------------------------------
c
c     some useful constants, allow easy switch single <--> double precision
c
      data zero,one/0.0,1.0/
c
c------------------------------------------------------------------------------
c
c     iteration parameters:
c      itbcsl  = max.nr. iterations in 3 pairing parameters
c      iterv2  = max.nr. iterations in 'vv2'
c      epsf1   = precision on pairing equation 'f1'
c      epsf2   = precision on pairing equation 'f2'
c      epsf3   = precision on pairing equation 'f3'
c      epsv2   = precision on 'vv2'
c      epsdv2  = limit on change of energy in determining 'vv2'
c      detmin  = limiting minimum value for the determinant in linear eq.
c      ddllim  = limit on stepsize in 'delta'
c      dxllim  = limit on stepsize in 'xlam'
c      dx2lim  = limit on stepsize in 'xlam2'
c      dampin  = general reduction of the step, top stabilize
c      delstp  = step in 'delta' for first derivative in finite difference
c      xlstp   = step in 'xlam' for first derivative in finite difference
c      xl2stp  = step in 'xlam2' for first derivative in finite difference
c
      data itbcsl,iterv2,epsf1,epsf2,epsf3,epsv2,epsdv2
     &    /    50,    50,     5*1.0D-8                 /
      data detmin,ddllim,dxllim,dx2lim,dampin
     &    /1.0D-1,  0.2,    3.0,   1.0, 0.8/
      data delstp,xlstp,xl2stp
     &    /  0.01,  0.01,  0.01/
c
c---------------------------------------------------------------------
c
c     preparation and tests
c
      if(nspg.gt.kmaxln) stop ' last level > "KMAXLN" in BCSLIP'
      if(ilipng.ne.1) then
       xlam2 = zero
       dxlam2 = zero
      endif
c
c    iteration of the three pairing parameters 'delta', 'xlam', and 'xlam2' .
c
      do 1000 itbcs=1,itbcsl
      if(iftest.gt.0)
     &   write(6,'(a,i5,a)')  ' iteration nr.',itbcs,' :'
c
c     initial call to pairing equations
c
      call eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,
     &           vv2,vv2,erspln,suv2rs,iftest,ilipng,
     &           iterv2, epsv2,epsdv2,dampin,
     &           f1,f2,f3)
      if(iftest.gt.0) then
        write(6,'(a,3g12.3/a,3g12.3)')
     &   ' delta,xlam,xlam2=',delta,xlam,xlam2,
     &   '         f1,f2,f3=',f1,f2,f3
      endif
c
c     test termination
c
      if(itbcs.gt.1 .and. abs(f1).lt.epsf1 .and. abs(f2).lt.epsf2
     &   .and. abs(f3).lt.epsf3) goto 1099
c
c     derivatives by finite differences
c
      call eqdef(spen,ph,nspd,nspg,npart,gpair,delta*(one+delstp),
     &            xlam,xlam2,vv2,vln,erspln,suv2rs,iftest,ilipng,  ! vln = aux
     &            iterv2, epsv2,epsdv2,dampin,
     &            f1del,f2del,f3del)
      call eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam+xlstp,
     &            xlam2,vv2,vln,erspln,suv2rs,iftest,ilipng,       ! vln = aux
     &            iterv2, epsv2,epsdv2,dampin,
     &            f1xl,f2xl,f3xl)
      if(ilipng.eq.1)
     & call eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2+xl2stp,
     &            vv2,vln,erspln,suv2rs,iftest,ilipng,             ! vln = aux
     &            iterv2, epsv2,epsdv2,dampin,
     &            f1xl2,f2xl2,f3xl2)
      ddelac = delta*delstp
      ddf1   = (f1del-f1)/ddelac
      ddf2   = (f2del-f2)/ddelac
      dlf1   = (f1xl-f1)/xlstp
      dlf2   = (f2xl-f2)/xlstp
      if(ilipng.eq.1) then
        dl2f1  = (f1xl2-f1)/xl2stp
        dl2f2  = (f2xl2-f2)/xl2stp
        ddf3   = (f3del-f3)/ddelac
        dlf3   = (f3xl-f3)/xlstp
        dl2f3  = (f3xl2-f3)/xl2stp
c
c     solve the linear equations for the case of 3 steps
c
c   define matrix in orders
c      derv(1,1) = ddf1
c      derv(2,1) = ddf2
c      derv(3,1) = ddf3
c      derv(1,2) = dlf1
c      derv(2,2) = dlf2
c      derv(3,2) = dlf3
c      derv(1,3) = dl2f1
c      derv(2,3) = dl2f2
c      derv(3,3) = dl2f3
c
c    matrix of inverse is computed explicitely and stored on bij
c
        det    =(ddf1*dlf2*dl2f3-ddf1*dl2f2*dlf3-ddf2*dlf1*dl2f3
     #          +ddf2*dl2f1*dlf3+ddf3*dlf1*dl2f2-ddf3*dl2f1*dlf2)
        det    = sign(max(abs(det),detmin),det)        ! lower limit for 'det'
        detinv = one/det
        b11 = (dlf2*dl2f3-dl2f2*dlf3)*detinv
        b12 = -(dlf1*dl2f3-dl2f1*dlf3)*detinv
        b13 = (dlf1*dl2f2-dl2f1*dlf2)*detinv
        b21 = -(ddf2*dl2f3-dl2f2*ddf3)*detinv
        b22 = (ddf1*dl2f3-dl2f1*ddf3)*detinv
        b23 = -(ddf1*dl2f2-dl2f1*ddf2)*detinv
        b31 = (ddf2*dlf3-dlf2*ddf3)*detinv
        b32 = -(ddf1*dlf3-dlf1*ddf3)*detinv
        b33 = (ddf1*dlf2-dlf1*ddf2)*detinv
ctest      if(iftest.ge.1) write(6,'(a,3g12.4,2(/19x,3g12.4))')
ctest     &   ' derivative matrix:',derv
ctest      call  MATIN2 (derv,3,3,NERROR,DET)
ctest      if(iftest.ge.1)
ctest     &  write(6,'(a,9g12.4)')
ctest     &   'd:',derv,
ctest     &   'b:',b11,b21,b31,b12,b22,b32,b13,b23,b33
ctest      ddelta = -derv(1,1)*f1-derv(1,2)*f2-derv(1,3)*f3
ctest      dxlam  = -derv(2,1)*f1-derv(2,2)*f2-derv(2,3)*f3
ctest      dxlam2 = -derv(3,1)*f1-derv(3,2)*f2-derv(3,3)*f3
        ddelta = -b11*f1-b12*f2-b13*f3
        dxlam  = -b21*f1-b22*f2-b23*f3
        dxlam2 = -b31*f1-b32*f2-b33*f3
        if(iftest.gt.0)
     &   write(6,'(a,4g12.3)')
     &    ' det,ddelta,dxlam,dxlam2=',det,ddelta,dxlam,dxlam2
      else
c
c     case of BCS, step in two quantities
c
        det    = -(dlf1*ddf2-dlf2*ddf1)
        det    = sign(max(abs(det),detmin),det)      ! lower limit for 'det'
        ddelta = (f2*dlf1-f1*dlf2)/det
        dxlam  = (f1*ddf2-f2*ddf1)/det
      endif
c
c     limit stepsize if necessary
c
      ddelta = dampin*sign(min(abs(ddelta/delta),ddllim),ddelta)
      dxlam  = dampin*sign(min(abs(dxlam ),dxllim),dxlam )
      if(ilipng.eq.1)
     &   dxlam2 = dampin*sign(min(abs(dxlam2),dx2lim),dxlam2)
c
c     step forward
c
      delta  = delta*exp(ddelta)
      xlam   = xlam  + dxlam
      if(ilipng.eq.1)
     &    xlam2  = xlam2 + dxlam2
c
 1000 continue
c
      if(iftest.ge.0) write(6,'(a)')
     & ' no convergence of pairing parameters in BCSLIP'
c
 1099 continue
      if(iftest.ge.0) write(6,'(a,i4,6(1pg11.3))')
     & ' BCSLIP: iter,delta,xlam,xlam2,f1,f2,f3=',
     &   itbcs,delta,xlam,xlam2,f1,f2,f3
c
c     compute some observables after convergence
c
      suv    = zero
      sv2    = zero
      sv4    = zero
      suv2   = zero
      sspe   = zero
      do 200 k=nspd,nspg
      v2     = vv2(k)
      uv     = sqrt(v2-v2*v2)
      vln(k) = uv                                     ! auxiliary storage
      xmult  = ph(k)
      suv    = suv    + xmult*uv
      sv2    = sv2    + xmult*v2
      sv4    = sv4    + xmult*v2*v2
      suv2   = suv2   + xmult*(v2 - v2*v2)
      sspe   = sspe   + xmult*spen(k)*v2
  200 continue
cold      ecorr  = -delta*delta/gpair-gpair*sv4-4.0*xlam2*suv2+gpair*sv2
      ecorr  = -delta*delta/gpair-(4.0*xlam2-gpair)*suv2
      entot  =  sspe - gpair * npart/2
      epair1 = -gpair*sv2
      efluc  = -4.*xlam2*suv2
      epair  = -gpair*suv
cold      etot   = 2.0*sspe+epair1+efluc+epair
cold      epairr = epair
      epair1 = epair1+efluc+epair
c
c           corrected coefficients, vln(i)=<a^+_i a_i>_ln
c
      if(ilipng.eq.1) then
        sum2=0.
        sum4=0.
        sumc=0.
        sumc2=0.
c
        do 120 k=nspd,nspg
        v2    = vv2(k)
        xmult = ph(k)
cold        uv    = sqrt(v2-v2*v2)
        uv    = vln(k)
        s     = uv+uv
        c     = v2+v2-one
        sum2  = sum2  + xmult*s**2
        sum4  = sum4  + xmult*s**4
        sumc  = sumc  + xmult*s*s*c
        sumc2 = sumc2 + xmult*(s*c)**2
  120   continue
c
        od=sum2*0.25/(sumc2*sum2-sumc*sumc+0.5*sum2*(sum2*sum2-sum4))
c
        do 130 k=nspd,nspg
        v2    = vv2(k)
        uv    = vln(k)
cold        uv    = sqrt(v2-v2*v2)
        s     = uv+uv
        c     = v2+v2-one
        vln(k)= v2+od*s*s*(c*sum2-sumc)
  130   continue
      else
        do 135 k=nspd,nspg
        vln(k)= vv2(k)
  135   continue
      endif
c
      return
      end
c-----eqdef-------------------------------------------------------------------
c
      subroutine eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,
     &                 vv2in,vv2out,erspln,suv2rs,iftest,ilipng,
     &                 iterv2, epsv2,epsdv2,dampin,
     &                 f1,f2,f3)
c
c---------------------------------------------------------------------
c
ctest      implicit double precision (a-h,o-z)
c
ctest      parameter (kmaxln=46)
      include 'params.inc'
      parameter (kmaxln=ksl2)
c
c---------------------------------------------------------------------
c
c     routine to compute the three pairing equations, as explained
c     in the header of 'BCSLIP'.
c
c     input:
c       spen   = array of s.p. energies
c       ph     = array of multiplicities, e.g. a j=7/2 state has ph = 4 .
c       nspd   = lowest level in pairing
c       nspg   = highest level in pairing
c       npart  = particle number
c       gpair  = strength of pairing force
c       delta  = pairing gap
c       xlam   = initial value for Fermi energy
c       xlam2  = initial value for fluctuation constraint
c       vv2in  = array of initial internal occupation probabilities
c       erspln - Skyrme response energy for 'lambda_2'
c       suv2rs -   u^2*v^2 from response routine, for normalization
c       iftest = switch for testprints
c       ilipng = switch between BCS and Lipkin-Nogami,
c                 1 --> f3=0
c       iterv2 = max.nr. iterations in 'vv2'
c       epsv2  = precision on 'vv2'
c       epsdv2 = limit on change of energy in determining 'vv2'
c       dampin = general reduction of the step, top stabilize
c     output:
c       f1     = delta-gpair*sum(u*v) =  gap equation
c       f2     = sum(2v**2)-n         =  particle number condition
c       f3     = xlam2-gpair/4*(sum(u**3*v)sum(v**3*u)-
c                               sum(u**4*v**4))/((sum(u**2*v**2))**2 -
c                               sum(u**4*v**4))
c              = p.n. fluctuation
c       vv2out = new occupation weights
c
c---------------------------------------------------------------------
c
c
      dimension spen(kmaxln),vv2in(kmaxln),vv2out(kmaxln),ph(kmaxln)
c
      data zero,one/0.0,1.0/
      data itbisc/25/                   ! max. number of bisection
c      data zero,half,one/0.0,0.5,1.0/
c
c---------------------------------------------------------------------
c
c     preparation
c
      fnn    = npart
      if(nspg.gt.kmaxln) stop ' last level > "kmaxln"'
      d2     = delta*delta
      alpha  = 4.*xlam2-gpair
      alphd2 = alpha*half*d2
c
c     reset accumulators for combinations of 'u' and 'v'.
c
      sv2   = zero                                  ! sum v^2
      suv   = zero                                  ! sum u*v
      if(ilipng.eq.1) then
       suv4  = zero                                  ! sum u^4 v^4
       suv2  = zero                                  ! sum u^2 v^2
       su3v  = zero                                  ! sum u*v u^2
       suv3  = zero                                  ! sum u*v v^2
      endif
c
      do 500 k=nspd,nspg                            ! each state separately
c
c     iteration loop to compute occupation 'vv2'
c
cold      e     = spen(k)-xlam
cold      v2    = vv2in(k)                                ! start from previous value
cold      x     = e+alpha*v2
coldc
cold      do 6 it=1,iterv2
cold      eq    = sqrt(d2+x*x)
cold      v2    = half-half*x/eq
cold      f     = e+alpha*v2-x
cold      derf  = -alphd2/(eq*eq*eq)-one
cold      dx    = dampin*f/derf         ! optionally switch on if unstable
coldc      dx    = f/derf
cold      if(abs(f).lt.epsv2 .and. abs(dx).lt.epsdv2) go to 7
cold      x     = x-dx
cold    6 continue
c                                      one tangential as startup
      e     = spen(k)-xlam
      x     = e
      f     = e+alpha*(0.5-0.5*x/sqrt(d2+x*x))-x
      dfdx  = -alpha*0.5*d2/(d2+x*x)**1.5-1.0
      dx    = f/dfdx
      x     = x-dx
c                                      save first point, with sign
      if(f.gt.zero) then
        fpos  = f
        xpos  = x+dx
        fneg  = zero
      elseif(f.lt.zero) then
        fneg  = f
        xneg  = x+dx
        fpos  = zero
      endif
c                                      secant step, embracing f=0
      do 6 it=1,iterv2
      f     = e+alpha*(0.5-0.5*x/sqrt(d2+x*x))-x
      if(abs(f).lt.epsv2 .and. abs(dx).lt.epsdv2) go to 7
      if(f.gt.zero) then
        if(fneg.ne.zero) then
          dfdx  = (fneg-f)/(xneg-x)
        else
          dfdx  = (fpos-f)/(xpos-x)
        endif
        fpos  = f
        xpos  = x
      else
        if(fpos.ne.zero) then
          dfdx  = (fpos-f)/(xpos-x)
        else
          dfdx  = (fneg-f)/(xneg-x)
        endif
        fneg  = f
        xneg  = x
      endif
      dx    = f/dfdx
      x     = x-dx
c      write(6,'(a,i4,4f14.4)')
c     & ' state #k: it,e,x,f,dx=',it,e,x,f,dx
    6 continue
c
c      the secant step did not succeed. use bisection
c
      do it=1,itbisc
        xnew  = 0.5*(xpos+xneg)
        f     = e+alpha*(0.5-0.5*x/sqrt(d2+x*x))-x
c       if(iab.ge.0)
c     &   write(6,'(a,i4,2(a,g11.3))')
c     &    ' bisection it=',it,'  :  xnew=',xnew,'  f=',f
        if(abs(f).lt.epsv2) then
           go to 7
        elseif(f.lt.zero) then
          xneg  = xnew
          fneg  = f
        elseif(f.gt.zero) then
          xpos  = xnew
          fpos  = f
        endif
      enddo
c
      write(6,'(a,4(1pg12.4))')
     & ' EQDEF: max. iter. exceeded, brackets(f,x)=',fneg,xneg,fpos,xpos
c
    7 continue
c
      v2    = 0.5-0.5*x/sqrt(d2+x*x)
      if(iftest.gt.1)
     & write(6,'(a,i5,a,i5,g12.3)')
     &  ' state nr.',k,' :  iterv,v2',it,v2
      vv2out(k)=v2
c
c     accumulate sums with new 'v2'
c
      xmult = ph(k)
      u2    = one-v2
      uv    = sqrt(u2*v2)
      sv2   = sv2  + xmult*v2
      suv   = suv  + xmult*uv
      if(ilipng.eq.1) then
       suv4  = suv4 + xmult*u2*u2*v2*v2
       suv2  = suv2 + xmult*u2*v2
       su3v  = su3v + xmult*uv*u2
       suv3  = suv3 + xmult*uv*v2
      endif
c
  500 continue
c
c     final assembling of the three pairing equations
c
      f1     = delta-gpair*suv                    ! gap equation
      f2     = sv2+sv2-fnn                        ! particle number
c     
c     the numerator from Skyrme response is a bit special.
c     it was computed with an older "v2". The corresponding
c     "sum V2*U2" is stored in "suv2rs". In order to avoid
c     all to frequent recomputation of that expensive piece,
c     we rescale it here just by the new "suv2".
c
      if(ilipng.eq.1) then
       denom  = (suv2*suv2-suv4)                   !  denominator/32
       enume1 = -gpair/4.0*(su3v*suv3-suv4)        !  1. numerator/32
       enume2 = suv2*suv2*erspln/(suv2rs*suv2rs*32.0)   !  2. numerator/32
       f3     = xlam2+(enume1-enume2)/denom        ! p.n. fluctuations
cold   f3=xlam2-gpair/4.*(su3v*suv3-suv4)/(suv2*suv2-suv4)
       endif
c
      return
      end
c================================================================
      function gavrg(ibcs,iz,in,iprot,rhosc,n1,n2,xmeff,a0)
c
c         gavrg calculates the strength of the monopole pairing
c         interaction using the average gap method. the average
c         gap model is based on moller and nix, nucl. phys.a536
c         (1992) 20. the average pair density is computed using
c         the semi-classical approximation  for   the  harmonic
c         oscillator model modified  by  the  isospin-dependent
c         factor derived based on the wigner-kirkwood expansion
c         for the woods-saxon potential (+ s.o. + coulomb).
c
c         input:    iz:    proton  number
c                   in:    neutron number
c                   iprot: 0 (neutrons) or 1 (protons)
c                   n1:    number of first h.f. level used
c                          in the bcs (ln) calculations.
c                          (preferably n1=1)
c                   n2:    number of last h.f. level used
c                          in the bcs (ln) calculations.
c                          (preferably n2=npart, i.e., iz or in)
c                   xmeff: effective mass (=1 in the woods
c                          saxon model)
c                   ibcs:  1 (bsc) or .ne.1 (ln)
c                   a0:    surface parameter in a Saxon Woods shape,
c                          is 0.58 * surface width from formfactor
c
c        output:    gavrg: pairing strength
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c      written by w. nazarewicz, 1994                               c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      implicit real*8 (a-h,o-z)
      include 'params.inc'
c
      dimension alphap(2),bethap(2)

      data alphap /5.13,-2.35/,
     :     bethap /-1.28e-2,-2.38e-4/
      data zero/0.0/
cold      data edr/0.33333333/
c
c----------------------------------------------------------------------
c
      a=iz+in
      ntot=iz
      if (iprot.eq.0) ntot=in
c
      if(rhosc.eq.zero) then
c
c          pair density at the fermi level
c
ctest      write(6,'(a/4i4,1f8.3,2i4,2f8.3)')
ctest     &   ' ibcs,iz,in,iprot,rhosc,n1,n2,xmeff,a0',
ctest     &     ibcs,iz,in,iprot,rhosc,n1,n2,xmeff,a0
      alpha=alphap(iprot+1)
      betha=bethap(iprot+1)
      fac=(2.*ntot/a)**(1.d0/3.)
      part3=(3.*ntot)**(2.0/3.)
      hom=41/a**(1.0/3)*fac
      tz=(in-iz)/a
      tz2=tz*tz
cold      rhosc=part3/hom/(1.-0.25/part3)/2.*(1.+alpha*tz2+betha*a*tz2)
c4      rhosc=part3/hom/(1.-4.0/part3)/2.*(1.+alpha*tz2+betha*a*tz2)
      rhosc=part3/hom/(1.-1.0/part3)/2.*(1.+alpha*tz2+betha*a*tz2)
c
c                here: correction for the effective mass
c
       r=1.14*a**(1.d0/3)
       eta=a0/r
c
       xm=log(xmeff)
       pieta=(3.1415926*eta)**2
       xmm=xmeff*(1.-eta*xm)**3*(1.+2*eta*xm*pieta/(1.+pieta))
      rhosc=rhosc*xmm**1.5
c
      endif
      robar=rhosc
c

      if (ibcs.eq.1) then
c
c              bcs variant
c
         if (iprot.eq.0) delavr=13.3/sqrt(a)
         if (iprot.eq.1) delavr=13.9/sqrt(a)

      else
c
c              ln variant
c
         if (iprot.eq.0) delavr=8.6/sqrt(a)
         if (iprot.eq.1) delavr=9.9/sqrt(a)
caltern         if (iprot.eq.0) delavr=3.7*a**(-edr)
caltern         if (iprot.eq.1) delavr=4.3*a**(-edr)
c
      endif
c
      np=ntot/2
c
      y1=(-np+n1-1)/robar
      y2=(-np+n2  )/robar
c
      gavrg=2./robar/(log(sqrt(y2**2+delavr**2)+y2)-
     *                log(sqrt(y1**2+delavr**2)+y1))
c
      return
      end

c-----lnresp------------------------------------------------------
c
      subroutine lnresp
c
      include 'params.inc'
c
c     the contribution from the Skyrme residual interaction to
c     the "lambda_2" in Lipkin-Nogami Pairing. 
c
c     computes the response-densities to the particle number operator
c     using pairing weights from 'weight'. 
c     uses the densties 'rho','tau','sjs' as intermediate storage
c     (note that these fields are destroyed).
c     the response Hamiltonian is constructed as in the rpa-package,
c     for the details see 'Annalen der Physik, 1 (1992) 598 & 632'.
c
c     information is handled via common:
c       'weight'       <-- /shells/
c       'rho,tau,...'  <-- /wavfun/
c
c     results are returned via common /shells/
c      rspln    = U^2*V^2 response energy to be used in Lipk.Nogami
c      suv2r    = sum U^2*V^2, to be used in Lipkin-Nogami
c
c     the Lipkin-Nogami does only proton-proton pairing and
c     neutron-neutron pairing. thus only the Z^2 and N^2 response
c     is required. We use nonetheless the typed formulae from the
c     rpa routine setting explicitely the conjugate particle type
c     to zero. this is computational overhead but safer to program.
c     The procedure is extensivelycommented for the first term
c     'xpotp' and 'xpotn'. The other terms run analogously.
c
c
      include 'wavfun.inc'
      include 'potinv.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
      dimension wgtres(ksl2)             ! to carry the response weights
c
c     data xminlw/1.0e-20/
c
c--------------------------------------------------------------------
c
c     compute the response weights (2uv)**2 and according densities
c
      suv2p  = zero                            ! first the protons
      do 10 n=1,nmaxp
	v2 = weight(n)
	wgtres(n) = 4.0*(v2-v2*v2)    / 2.0    !   factor 1/2 to compensate
c                                              !   factor 2 in densit
        js=jsp(n)
        ll=js/2
        if(2*ll.eq.js) js=js-1
        mmult = (js+1)/2                      ! half the multiplicity
	suv2p  = suv2p  +  mmult*(v2-v2*v2)
   10 continue
      suv2n  = zero                            ! now the neutrons
      do 20 n=1,nmaxn
	nc = n+nmaxp
	v2 = weight(nc)
	wgtres(nc) = 4.0*(v2-v2*v2)    / 2.0   !   factor 1/2 to compensate
c                                              !   factor 2 in densit
        js=jsp(n)
        ll=js/2
        if(2*ll.eq.js) js=js-1
        mmult = (js+1)/2                      ! half the multiplicity
	suv2n  = suv2n  +mmult*(v2-v2*v2)
   20 continue
      suv2r(1) = suv2p
      suv2r(2) = suv2n
c
      do 30 i=1,ngrid+ngrid                   ! save old densities
	 aux1(i) = rho(i)
	 aux2(i) = tau(i)
	 aux3(i) = sjs(i)
   30 continue
      call densit(wgtres,v2occ)
c
c     terms depending on 'rho'
c
      t0x    = t0+0.5*t0*x0                                   ! b0
      t0y    = 0.5*t0+t0*x0                                   !b0p
      t3x    = (t3+t3*0.5*x3)*(power+one)*(power+2.0)/12.0
      t3y1   = (0.5*t3+t3*x3)/6.0
      t3y2   = 0.5*t3y1*power*(power-one)
      xpotp  = zero
      xpotn  = zero
      do 110 ir=1,ngrid
c                                    first the proton contribution
      irq    = ir                        ! position of proton densities
      irco   = ir+ngrid                  ! position of neutron densities
      rq     = aux1(irq)                 ! density for like sort
      rc     = aux1(irco)                ! density for other sort, yet needed
      rs     = rq+rc
      rpowm2 = rs**(power-2.0)
      prq    = rho(irq)                  ! response dens. for like sort
      prc    = zero                      ! response dens. for other sort --> 0
      prs    = prq+prc
      xpotp  = xpotp +   rsquar(ir)*prq*
     &          ( t0x*prs-t0y*prq
     &           +t3x*rpowm2*(rs*rs)*prs
     &           -t3y1*rpowm2*rs*(power*(rq*prs+rq*prq+rc*prc)
     &                            +rs*prq)  
     &           -t3y2*rpowm2*prs*(rq*rq+rc*rc)  )
c                                    now the neutron contribution
      irq    = ir+ngrid                  ! position of neutron densities
      irco   = ir                        ! position of proton densities
      rq     = aux1(irq)                 ! density for like sort
      rc     = aux1(irco)                ! density for other sort, yet needed
      rs     = rq+rc
      rpowm2 = rs**(power-2.0)
      prq    = rho(irq)                  ! response dens. for like sort
      prc    = zero                      ! response dens. for other sort --> 0
      prs    = prq+prc
      xpotn  = xpotn +   rsquar(ir)*prq*
     &          ( t0x*prs-t0y*prq
     &           +t3x*rpowm2*(rs*rs)*prs
     &           -t3y1*rpowm2*rs*(power*(rq*prs+rq*prq+rc*prc)
     &                            +rs*prq)  
     &           -t3y2*rpowm2*prs*(rq*rq+rc*rc)  )
  110 continue
      xpotp  = frpirs*xpotp
      xpotn  = frpirs*xpotn
c
c     terms depending on 'tau'
c
      t12co  = 0.25*(t1+t1*0.5*x1+t2+t2*0.5*x2)                     ! b1
      t12q   = t12co-0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)               ! b1-b1p
      xkinp  = zero
      xkinn  = zero
      do 120 ir=1,ngrid
      irq    = ir                                  ! protons
      pt     = tau(irq)
      ptc    = zero    
      xkinp  = xkinp +  rsquar(ir)*(
     &              rho(irq)*(t12q*pt+t12co*ptc)  )
c    &              +rho(irc)*(t12q*ptc+t12co*pt) )    ! other dens. --> zero
      irq    = ir+ngrid                            ! neutrons
      pt     = tau(irq)
      ptc    = zero    
      xkinn  = xkinn +  rsquar(ir)*(
     &              rho(irq)*(t12q*pt+t12co*ptc)  )
c    &              +rho(irc)*(t12q*ptc+t12co*pt) )    ! other dens. --> zero
  120 continue
      xkinp  = 2.0*frpirs*xkinp      ! this term occurs twice:
      xkinn  = 2.0*frpirs*xkinn      ! due to doubel variation  rho <--> tau
c
c     coulomb contribution
c     (the fields 'sqb' and 'sqm' are used as work-space inside 'coulmb'.
c      the field 'binv' is used as workspace for carrying the coulomb
c      potential.)
c
      if(acoul.ne.zero) then
       call coulmb(rho(1),binv,acoul)
       xcol   = zero
       do 240 ir=1,ngrid
  240  xcol   = xcol + rsquar(ir)*rho(ir)*binv(ir)
       xcol   = frpirs*xcol
c
c      coulomb exchange term in local density approximation
c
       xcex   = zero
       if(ifex.eq.1) then
        excfac = -1.41803*abs(acoul)/3.0
c       excfac = zero
        twdthr = -0.66666666667
        do 250 ir=1,ngrid
        rp     = aux1(ir)
        if(rp.gt.xsmall)
     &    xcex = excfac*rp**twdthr*rho(ir)*rho(ir) * rsquar(ir)
     &          + xcex
  250   continue
	xcex   = frpirs*xcex
       endif
      endif
c
c     terms depending on laplace->'rho'
c
      if(imfold.eq.1) then
	stop ' option "imfold.eq.1" not valid with Lipkin-Nogami'
      endif
      t12co  = -0.125*(3.0*(t1+t1*0.5*x1)-t2-t2*0.5*x2)             ! -b2
      t12q   = t12co+0.125*(3.0*(0.5*t1+t1*x1)+0.5*t2+t2*x2)        ! -b2+b2p
      xpd2p  = zero
      xpd2n  = zero
      do 130 ir=2,ngrid
      irq    = ir                             !  the proton contribution
      dpr    = d2r(irq)
      dprc   = zero
      xpd2p  = xpd2p +  rsquar(ir)*(
     &                   rho(irq)*(t12q*dpr+t12co*dprc )   )
c    &                   +rho(irc)*(t12q*dprc+t12co*dpr)     ) ! other --> 0
      irq    = ngrid+ir                       !  the neutron contribution
      dpr    = d2r(irq)
      dprc   = zero
      xpd2n  = xpd2n +  rsquar(ir)*(
     &                   rho(irq)*(t12q*dpr+t12co*dprc )   )
c    &                   +rho(irc)*(t12q*dprc+t12co*dpr)     ) ! other --> 0
  130 continue
      xpd2p  = frpirs*xpd2p
      xpd2n  = frpirs*xpd2n
c
c     terms depending on laplace->'tau'           
c   
      if(t9.ne.zero) then
	stop ' option "t9.ne.0" not valid with Lipkin-Nogami'
      endif
c
      if(t4.ne.zero) then  
c
c     terms depending on 'sjs'
c
      t4co   = -b4
      t4q    = -b4-b4p
      xspop  = zero
      xspon  = zero
      do 40 ir=2,ngrid
      irq    = ir                               ! the proton contribution
      ps     = sjs(irq)
      psc    = zero
      xspop  = xspop +  rsquar(ir)*(
     &                          rho(irq)*(t4q*ps+t4co*psc)  )
c    &                         +rho(irc)*(t4q*psc+t4co*ps)  ) ! other --> 0
      irq    = ngrid+ir                         ! the neutron contribution
      ps     = sjs(irq)
      psc    = zero
      xspon  = xspon +  rsquar(ir)*(
     &                          rho(irq)*(t4q*ps+t4co*psc)  )
c    &                         +rho(irc)*(t4q*psc+t4co*ps)  ) ! other --> 0
   40 continue
      xspop  = 2.0*frpirs*xspop                 ! appers twice:
      xspon  = 2.0*frpirs*xspon                 ! due to rho <--> sjs
      endif
c
c     now the final accumulation
c     and optional testprint
c
      rspln(1) = xpotp+xkinp+xpd2p+xcol+xcex+xspop
      rspln(2) = xpotn+xkinn+xpd2n+xspon
c      if(iwr.ge.-1) then
      if(iwr.gt.-1) then
        write(6,'(a,6g11.3/a,g11.3,a,g11.3)')
     &  ' Lip.Nog. response  protons: pot,kin,dp2,col,cex,spo=',
     &   xpotp,xkinp,xpd2p,xcol,xcex,xspop,
     &  '  total=',rspln(1),' suv2=',suv2r(1)
        write(6,'(a,4g11.3/a,g11.3,a,g11.3)')
     &  ' Lip.Nog. response neutrons: pot,kin,dp2,spo=',
     &   xpotp,xkinp,xpd2p,xspop,
     &  '  total=',rspln(2),' suv2=',suv2r(2) 
c       write(6,'(a,30(/1x,10f8.3))')
c    &   ' weights:',(weight(n),n=1,nmaxto)
      endif
c
      return
      end

c-----semrho------------ part of the hartree-fock package ----------------------
c
      function semrho(eref,potmf,bmasmf)
c
      include 'params.inc'
c
c     semiclassical level density
c     for the potential 'potmf'
c     with effective mass 'bmasmf'
c     taken at the reference energy 'eref'.
c     grid parameters enter via 'grid.inc'
c
c
      include 'skypar.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
      dimension potmf(kg1),bmasmf(kg1)
c
c
c-------------------------------------------------------------------------------
c
      sum    = zero
      do 10 ir=2,ngrid
      bb     = one/bmasmf(ir)
      eover  = eref-potmf(ir)
      if(eover.gt.zero)
     &   sum    = rsquar(ir)*bb*sqrt(bb*eover) + sum
   10 continue
      semrho = frpirs*sum/(4.0*pi*pi)
c
      return
      end
clflpar      include 'fullpair.for'
cold      include 'hfb.for'
c-----hfbpot------------------------------------------------------
c
      subroutine hfbpot
c
      include 'params.inc'
c
c     computes the pair density
c
c       sum_n  u_n*v_n*|phi_n(r)|^2  ,
c
c     subsequently the pair potential on 'pairpt' from /wavfun/,
c     and finally the dynamic gaps on 'delatf' from /shells/
c     as expectation value of 'pairpt'.
c
c     'ipair' is used here to indicate the switch to surface pairing
c             5,6 --> volume pairing
c             7,8 --> surface pairing
c
c     Note: computation of 'pptmax' for HFB iteration with 
c           Tajima damping yet erroneous
c
c
      include 'wavfun.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'constn.inc'
c
      dimension epa(2)
c
      data chipro,chineu/10000.02,10000.02/
      LOGICAL, PARAMETER :: tpaircur=.FALSE.  ! .TRUE. to revive old version
      REAL(8), PARAMETER :: pairoffset=1D-100  !1D-20  !xsmall
      LOGICAL, PARAMETER :: folding=.TRUE.
      REAL(8), DIMENSION(ngrid) :: fpsi
c
c--------------------------------------------------------------------
c
c
c     zero field before accumulation
c
      IF(tpaircur) THEN
c       write(6,'(a)') ' HFB: weights for protons:'
c       write(6,'(5(1pg12.4))') (weight(i),i=1,nmaxp)
c       write(6,'(a)') ' HFB; weights for neutrons:'
c       write(6,'(5(1pg12.4))') (weight(i+nmaxp),i=1,nmaxn)
        do i=1,kgr2
          paircr(i) = zero
        enddo
c
c       accumulate:
c
c       iso-loop, initialise "shift constants" 'im,ia,ie' for the
c       densities, 'ngw' for the weights, and 'nadd' for the
c       wave functions.
c
        do 20 iso=1,2
        ie     = ngrid*iso
        im     = ie-ngrid+1
        ia     = im+1
        if(iso.eq.1) then
         nmax   = nmaxp
        else
         nmax   = nmaxn
        endif
        ngw    = (iso-1)*nmaxp
        nadd   = ngw*ngrid-ie
c
c       loop over states for given isospin.
c        'nshl' counts the states in harmonic oscillator ordering.
c        'js' is the according composed ang.mom.,
c        'll' is the orbital ang.mom.,
c        'wl' is the spin-orbit splitting,
c        'om' is the actual phase-space factor
c
        do 30 nshl=1,nmax
        js     = jsp(nshl)
        ll     = js/2
c        xlp    = float(ll*ll+ll+1)
c        wl     = wls(js)
        if(temprt.eq.zero) then
          v2     = weight(nshl+ngw)
!          uv     = sqrt(max(v2-v2*v2,1D-16))
          uv     = phpair(nshl+ngw)*sqrt(max(v2-v2*v2,pairoffset))
!          uv     = phpair(nshl+ngw)*sqrt(max(v2-v2*v2,xsmall**2))
        else
          ttan   = ttthrm(nshl+ngw)
          v2     = (weight(nshl+ngw)-0.5*(one-ttan))/ttan
          uv     = ttan*phpair(nshl+ngw)*sqrt(max(v2-v2*v2,xsmall**2))
        endif
        om     = zwj(js)*uv*0.5              ! factor 1/2 because "state>0"
c        write(*,*) ' iso,nshl,js,uv,om=',iso,nshl,js,uv,om
        nadd   = nadd+ngrid
        if(om.lt.xsmall) goto 30
c                                    the case r=0.0 needs special treatment
        if(ll.lt.1) then
         x      = wfs(im+nadd)
         IF(.NOT.thfbfold) paircr(im) = (x*x)*om + paircr(im)
        endif
c                                    now all r>0.0
        r      = zero
        IF (thfbfold) THEN
!            WRITE(6,'(a,i4)') 'HFBPOT js=',js
            call Gauss_folding_m(wf(nadd+1:nadd+ngrid), 
     &       fpsi(nadd+1:nadd+ngrid), rstep, ngrid, js, hfbwidth, 
     &            MAXVAL(jsmax(1:2)))
            DO ir=ia,ie
              r      = r+rstep
              wfa    = wffold(nadd+ir)/r
              paircr(ir) = om*(wfa*wfa) + paircr(ir)
            END DO
        ELSE
            DO ir=ia,ie
              r      = r+rstep
              wfa    = wf(nadd+ir)/r
              paircr(ir) = om*(wfa*wfa) + paircr(ir)
            END DO
        END IF
   30   continue

!    extrapolate to r=0.0
        IF(thfbfold) CALL extrap(paircr(im),1D0)

   20   continue
      END IF
c
c      write(*,*) 'proton paircr:',(paircr(ir),ir=1,ngrid)
c
c     now multiply with strength to obtain local pair-potential
c
      IF(tpairFy) THEN
        CALL fy_hfbpot()
!        CALL fy_pairpot()
      ELSE
        if(ipair.le.6) then
          v0prot = -pairfp
          v0neut = -pairfn
          IF(tregpair .AND. iteration_number>50) THEN
            do i=1,ngrid
              geff = geff_regpair(efep,-cutfac,pot(i),bmass(i),v0prot) 
              pairpt(i)  = paircr(i) * geff
              regpairadd(i) = 1D0 !(geff/v0prot)**2
              in = i+ngrid
              geff =geff_regpair(efen,-cutfac,pot(in),bmass(in),v0neut) 
              pairpt(in)  = paircr(in) * geff
              regpairadd(in) = 1D0 !(geff/v0neut)**2
            enddo
          ELSE
            do i=1,ngrid
              pairpt(i)  = v0prot * paircr(i)
              in = i+ngrid
              pairpt(in) = v0neut * paircr(in)
            enddo
          END IF
        else
ctest          rho0pr = 1000.0                                 ! for tests
          v0prot = -pairfp
          v0neut = -pairfn
!          WRITE(6,'(a,3(1pg13.5))') 'V0pair=',v0prot,v0neut,rho0pr
          do i=1,ngrid
            in = i+ngrid
            rhofac = (rho(i)+rho(in))/rho0pr
            rhofac = rhofac**pairgm
            rhfacp = one-rhofac                               !new
     &           -2.0*paircr(i)**2/chipro**2                  !new
            rhfacn = one-rhofac                               !new
     &           -2.0*paircr(in)**2/chineu**2                 !new
            IF(tregpair .AND. iteration_number>50) THEN
              geff = geff_regpair(efep,efaddp,pot(i),bmass(i),
     &               v0prot*rhfacp) 
              pairpt(i)  = paircr(i) * geff
              regpairadd(i) = (geff/(v0prot*rhfacp))**2
              in = i+ngrid
              geff = geff_regpair(efen,efaddn,pot(in),bmass(in),
     &               v0neut*rhfacn) 
              pairpt(in)  = paircr(in) * geff
              regpairadd(in) = (geff/(v0neut*rhfacn))**2
!              STOP "regularized pairing not yet for DDDI pairing"
            ELSE
              pairpt(i)  = v0prot * paircr(i)  * rhfacp      !new
              pairpt(in) = v0neut * paircr(in) * rhfacn      !new 
            END IF
cold            rhofac = one-(rho(i)+rho(in))/rho0pr
cold            pairpt(i)  = v0prot * paircr(i)  * rhofac
cold            pairpt(in) = v0neut * paircr(in) * rhofac
          enddo
        endif
      END IF
c
c     finally compute the actual gaps as s.p. expectation values with
c     the pair potential
c
      nadd   = -ngrid
      do istate=1,nmaxto
        if(istate.le.nmaxp) then
          isoadd = 0
          iso    = 1
        else
          isoadd = ngrid
          iso    = 2
        endif
        nadd   = ngrid+nadd
        acc    = zero
        pptacc = zero
        IF(thfbfold) THEN
          do ir=1,ngrid
            wfa    = wffold(nadd+ir)
            pptacc = max(abs(pairpt(isoadd+ir)),pptacc)
            acc    = (wfa*wfa)*pairpt(isoadd+ir) + acc
          enddo
          pptmax(iso) = pptacc            
        ELSE
          do ir=1,ngrid
            wfa    = wf(nadd+ir)
            pptacc = max(abs(pairpt(isoadd+ir)),pptacc)
            acc    = (wfa*wfa)*pairpt(isoadd+ir) + acc
          enddo
          pptmax(iso) = pptacc            
        END IF
!        deltaf(istate) = (1D0-addnew)*deltaf(istate)-addnew*rstep*acc
        deltaf(istate) = MAX(-rstep*acc,1D-6)
      enddo
c
      if(iwr.ge.-1 .or. ecstab.NE.0D0) then
!      if(iwr.gt.-1 .or. ecstab.gt.0D0) then
!      if(iwr.gt.-3 .or. ecstab.ge.0D0 .and. .NOT.tpairFY) then
        epp = 0.0D0
        epn = 0.0D0
        do ir=1,ngrid
          vol = (ir-1)*rstep
          vol = vol*vol
          epp = paircr(ir)*pairpt(ir)*vol + epp
          ic = ir+ngrid
          epn = paircr(ic)*pairpt(ic)*vol + epn
        enddo
        epp = frpirs*epp
        epn = frpirs*epn
        epair = epp+epn
        epp0=epp
        epn0=epn
!        WRITE(6,'(a,5(1pg13.5))') 'HFBepair=',epair,epp,epn
!        CALL pair_stab(spenrg,deltaf,weight,ph,phpair,v2occ,epp3,
!     &                     0D0,nmaxp,kmaxn,1)
!        CALL pair_stab(spenrg(nminn),deltaf(nminn),weight(nminn),
!     &                 ph(nminn),phpair(nminn),v2occ(nminn),epn3,
!     &                     0D0,nmaxn,kmaxn,1)
!        WRITE(6,'(a,5(1pg13.5))') 'HFBepai2=',epp3+epn3,epp3,epn3
!        WRITE(6,'(100(5(1pg13.5)))') (pairpt(i),i=ngrid+1,2*ngrid)
!        WRITE(6,'(100(5(1pg13.5)))') (paircr(i),i=ngrid+1,2*ngrid)

      endif
c
c       augment by correction for stabilized pairing
c
c      write(*,*) ' ECSTAB=',ecstab
      if(ecstab.NE.0D0) then
ctest        write(*,*) ' EPP.EPN=',epp,epn
!        WRITE(6,'(a,200(1pg13.5))') 'deltaf:',deltaf(1:nmaxp)
        do istate=1,nmaxto
          if(istate.le.nmaxp) then
            epact = epp
          else
            epact = epn
          endif
          deltaf(istate) = deltaf(istate)
     &                    *(1D0+(ecstab/epact)**2)
          deltaf(istate) = min(deltaf(istate),5D0)
        enddo
!        WRITE(6,'(a,200(1pg13.5))') 'deltaf:',deltaf(1:nmaxp)
        epp2 = epp+ecstab**2/max(-epp,1D-8)
        epn2 = epn+ecstab**2/max(-epn,1D-8)
!        write(89,*) ' EPP.EPN=',epp,epp2,epn,epn2
      endif
      DO istate=1,nmaxto 
!         deltaf(istate)=MAX(deltaf(istate),1D-16)
         deltaf(istate)=MAX(deltaf(istate),1D-8)
      END DO
c
      if(iwr.gt.-1) then
        write(6,'(a,3(1pg12.4))')
     &     ' epair p/n/t=',epp,epn,epp+epn
        write(6,'(a)') ' dynamics gaps for protons, directly:'
        write(6,'(5(1pg12.4))') (deltaf(i),i=1,nmaxp)
        write(6,'(a)') ' dynamics gaps for neutrons, directly:'
        write(6,'(5(1pg12.4))') (deltaf(i+nmaxp),i=1,nmaxn)
        write(6,'(a)') ' pair current for protons:'
        write(6,'(5(1pg12.4))') (paircr(i),i=1,ngrid)
        write(6,'(a)') ' pair current for neutrons:'
        write(6,'(5(1pg12.4))') (paircr(i),i=ngrid+1,2*ngrid)
        write(6,'(a)') ' pair potential for protons:'
        write(6,'(5(1pg12.4))') (pairpt(i),i=1,ngrid)
        write(6,'(a)') ' pair potential for neutrons:'
        write(6,'(5(1pg12.4))') (pairpt(i),i=ngrid+1,2*ngrid)
        write(6,'(a)') ' contrib. to epair for neutrons:'
        write(6,'(5(1pg12.4))') (deltaf(i)*
     &     sqrt(weight(i)-weight(i)**2)*phpair(i),i=nmaxp+1,nmaxp+nmaxn)
        acc=0D0
        DO i=nmaxp+1,nmaxp+nmaxn
          acc=acc+deltaf(i)*sqrt(weight(i)-weight(i)**2)*phpair(i)*ph(i)
        END DO
        WRITE(6,*) 'summed to epair:',acc
      endif
c
      return
      end
c-----pairdn------------ part of the hartree-fock package ----------------------
c                        changed to fit into the axial H.F.code
c
      subroutine pairdn(e,delta,gw,ph,pairgw,wgthrm,v2occ,cutfac,
     &                  partnm,eferm,epair,avdelt,avdlt2,eqmin,
     &                  avg,efadd,rtmp,nmax,ndim,isosh,nblck,
     &                  iterin,eps,iab,tlipnog)
c
      include 'params.inc'
       COMMON /explicit/v2explicit,nexplicit,aprotnum,aneutnum
!      include 'skypar.inc'
c
c
c     the include file 'params.inc' allows
c     to make the whole subroutine double precision by
c     transfering a line 'implicit double precision (a-h,o-z)'.
c     Furthermore it carries the common block /CONSTN/ from which
c     the pseudo-constants ONE=1.0 and ZERO=0.0 are taken here.
c
c
c     * * * * * * * * * *   pairing with full interaction   * * * * * * * * * *
c
c     determines the fermi-energy 'eferm'
c     such that the "particle number condition"
c     is fulfilled. this is achieved by iteration.
c
c     input via list:
c      e       = array(1..nmax) of single particle energies.
c      delta   = array(1..nmax) of gap parameters.
c      gw      = array(1..nmax) of occupation probabilities,
c                input is used as initial guess for the iteration.
c      ph      = array(1..nmax) of the multiplicities of the single
c                particle levels.
c                This multiplicities excludes (!) the spin degeneracy,
c                e.g. a 1s1/2 state has ph = 1 .
c      pairgw  = space for dynamical cut-off factor
c      wgthrm  = thermal*pairing occupation factor tanh(...)*u*v
c      cutfac  = cutoff factor for soft pairing space
c                >0: factor for definition in terms of particle number
c                <0: abs.value defines energy band above Fermim energy
c      partnm  = desired particle number
c      eferm   = fermi-energy, input is used as initial guess;
c                if eferm=0.0 is input then a suitable initial guess
c                is computed from input 'gw'.
c      epair   = pairing energy
c      avdelt  = average gap  u*v weighted
c      avdltw  = average gap  v**2 weighted
c      eqmin   = minimum of quasiparticle energies (a measure for gap)
c      avg     = averag pairing force
c      efadd   = pairing band above 'eferm', in case of dynamical pair space
c                see routine  CUTOFF
c      rtmp    = temperature (in cooperation with BCS)
c      nmax    = maximum number of single particle states
c                (needs to be smaller than the parameter 'ndim')
c      ndim    = dimensioning of the arrays, particularly of 'vpairf',
c                in the calling routine
c      nblck   = number of blocked state (=0 if none is blocked)
c      isosh   = offset in spatial array according to isospin particle sort
c      iter    = number of pairing iterations,
c                abs(iter) is the number of iterations
c      eps     = end condition for pairing iteration;
c                iteration terminates if error in particle number
c                is smaller than 'eps' and if 'eferm' and/or 'delta'
c                are smaller than 'eps'.
c      iab     = print level; significant values are -2, -1, 0 and 1;
c                larger value prints more.
c
c     output via list:
c      gw      = array(1..nmax) of the occupation probabilities
c                finally achieved.
c      eferm   = fermi-energy.
c
      REAL(8) :: e(ndim),gw(ndim),ph(ndim),delta(ndim),pairgw(ndim)
      REAL(8) :: wgthrm(ndim),v2occ(ndim)
      LOGICAL,INTENT(IN) :: tlipnog

      REAL(8) :: gwsav(kmaxn),pwsav(kmaxn),v2sav(kmaxn),esav(ksl2)
c
      common /ctest/edirec(ksl2),efkincut,ifkincut
c
c     several numerical parameters, here for nuclei
c
c      data dltini/1.0/          ! initial gaps
c      data stepd /0.4/          ! stepsize for iteration of gap
      data stepl  /0.1/  ! 0.2   ! stepsize for iteration of Fermi energy
      data delefe /2.0/  ! 3.0   ! difference of Fermi energy for derivative
      data ifcut  /1/            ! switch for flexible cutoff
      data slpmin /1.0/  ! 2.0   ! minimum slope in dN/dE
c
c      data one,zero,xsmall/1.0,0.0,1.0E-20/
      data one,zero,half,xsmall/1D0,0D0,0.5D0,1D-100/  !1.0D-20/
c      data cutfac/2.00/         ! active pairing space , default 1.65
c
      data icount/0/
c      data ifedir/3/           ! switches sort of energy in cutoff
      data ifedir/0/           ! switches sort of energy in cutoff
!      data ifkincut/1/        ! switches to kinetic-energy cutoff
      data addpgw/1D0/        ! mixing of pairing weights
      data addmix/1D0/        ! mixing of occupation numbers
c
      logical,parameter ::  TPOLD=.FALSE.
      logical,parameter ::  tstcut=.TRUE.
      REAL(8), PARAMETER :: efermkin=35D0    ! assumed kinetic Fermi energy
      LOGICAL :: tswitchexplicit

c
c-------------------------------------------------------------------------------
c
      tswitchexplicit = nexplicit>0 .AND. ABS(partnm-aneutnum)<1D-6
      IF(tswitchexplicit) nexplicit=nexplicit-isosh+1
!      WRITE(6,'(a,i2,2(1pg13.5),L2)') 
!     &    'EXPLICIT:',nexplicit,partnm,aneutnum,tswitchexplicit
      IF(nblck>0 .AND. nexplicit>0) 
     &   STOP "NBLCK and NEXPLICIT incompatible"
      ifkincut=0      ! switches to kinetic-energy cutoff if 1
      esav=edirec
!
      IF(TPOLD .AND. .NOT.tphpair) 
     &   STOP 'incompatible TPOLD and TPHPAIR'
      IF(iterin<0) addpgw=0D0
      iter=ABS(iterin)
      if(iab.gt.-1) then
        write(6,'(/a,i4,(a,f9.4),a,i4/a,10(/1x,10f8.3))')
     &   ' *** pairing with dynamic gap entered with: iter=',iter,
     &   '  eferm=',eferm,' blocking state=',nblck,
     &   '  single particle energies:',
     &   (e(i),i=1,nmax)
        write(6,'(a,10(/1x,10f8.4))') '  delta:',(delta(i),i=1,nmax)
        write(6,'(a,10(/1x,10f8.4))') '    ph:',(ph(i),i=1,nmax)
        write(6,'(a,10(/1x,10f8.4))') 'pairgw:',(pairgw(i),i=1,nmax)
        write(6,'(a,10(/1x,10f8.4))') 
     &         ' spkin:',(edirec(isosh+i),i=0,nmax-1)
      endif
      IF(nblck>0) THEN
        partnm = partnm - 1D0                 ! blocking
        IF(eferm.NE.0D0) THEN
          IF(ph(nblck)==1D0) THEN
            gw(nblck) = 0D0
          ELSE
            v2b = v2occ(nblck)     ! old:  gw(nblck)-0.5D0/ph(nblck)
            gw(nblck) = v2b        ! old: *ph(nblck)/(ph(nblck)-1D0)
          END IF
        END IF
        ph(nblck) = ph(nblck) - 1D0           ! blocking
      END IF
      IF(tswitchexplicit) THEN
        partnm = partnm - v2explicit*2D0                 ! one level explicit
        WRITE(6,'(a,2(1pg13.5))') 
     &    'PART.GOAL:',partnm,partnm+v2explicit*2D0
        IF(eferm.NE.0D0) THEN
          IF(ph(nexplicit)==1D0) THEN
            gw(nexplicit) = 0D0
          ELSE
            v2b = v2occ(nexplicit)
            gw(nexplicit) = v2b   
          END IF
        END IF
        ph(nexplicit) = ph(nexplicit) - 1D0
      END IF
      icount = 1 + icount
      if(eferm.ne.zero) then
        gwsav(1:nmax)=gw(1:nmax)
        v2sav(1:nmax)=v2occ(1:nmax)
        ifsav = 1
      else
        icount=0
        ifsav = 0
      endif
c
c                               preset pairing weights
c
ctest      if(eferm.ne.zero .and. ifcut.eq.1 .and. icount.gt.120) then
      pwsav=1D0
      if(eferm.ne.zero .and. ifcut.eq.1 .AND. cutfac>0D0) then
       actnum = partnm + cutfac*partnm**0.6666667 
       delcut = 0.1
       efadd  = 20.0/partnm**0.33333333
       do n=1,nmax
         pwsav(n) = pairgw(n) 
       enddo
       if(ifedir.eq.0) then
         call cutoff(e,pairgw,ph,actnum,nmax,eferm,efadd,delcut,
     &                tstnum,iter,eps,iab)
!         WRITE(6,'(a,4(1pg13.5))') 'CUTOFF: actnum,efadd,tstnum=',
!     &        actnum,efadd,tstnum,cutfac,partnm
       elseif(ifedir.eq.1) then
         call cutoff(edirec(isosh),pairgw,ph,actnum,nmax,eferm,efadd,
     &                delcut,tstnum,iter,eps,iab)
       elseif(ifedir.eq.2) then
         call sprad(edirec(isosh),rmsact,gw,ph,isosh)
         write(6,'(a,2f8.4)')
     &      ' rmsact,guess=',rmsact,0.9*partnm**0.333333
         rmsadd = 0.8*(actnum**0.33333-partnm**0.33333)
         call cutoff(edirec(isosh),pairgw,ph,actnum,nmax,rmsact,rmsadd,
     &                delcut,tstnum,iter,eps,iab)
       elseif(ifedir.eq.3) then
         call cutpot(e,edirec,isosh)
         call cutoff(edirec,pairgw,ph,actnum,nmax,eferm,efadd,
     &                delcut,tstnum,iter,eps,iab)
       elseif(ifedir.eq.4) then
         do n=1,nmax
           pairgw(n) = 1.0
         enddo
       elseif(ifedir.eq.5) then
         call cutang(e,edirec,isosh)
         call cutoff(edirec,pairgw,ph,actnum,nmax,eferm,efadd,
     &                delcut,tstnum,iter,eps,iab)
       else
         stop 'IFEDIR > 5 not yet provided'
       endif
       if(icount.gt.300 .and. addpgw.lt.one) then
         do n=1,nmax
           pairgw(n) = (1D0-addpgw)*pwsav(n)+addpgw*pairgw(n)
         enddo
       endif
      else if(eferm.ne.0D0 .and. ifcut.eq.1 .AND. cutfac<0D0) then
        IF(ifkincut==1) THEN
!          efcut=-cutfac
!          efkincut=efermkin-cutfac
!          edirec(1:nmax)=esav(isosh:isosh+nmax-1)
          edirec(1:nmax)=e(1:nmax)   ! restore s.p. energy
          efkincut=eferm-cutfac
        ELSE
          edirec(1:nmax)=e(1:nmax)   ! restore s.p. energy
          efcut=eferm-cutfac
        END IF
        acc = 00D0
        IF(ifkincut==1) THEN
          pairgw(1:nmax) =
     &     1D0/(1D0+EXP((edirec(1:nmax)-efkincut)/(-0.1D0*cutfac)))
        ELSE
          pairgw(1:nmax) =
     &     1D0/(1D0+EXP((e(1:nmax)-eferm+cutfac)/(-0.1D0*cutfac)))
          IF(tlowercut) pairgw(1:nmax) = pairgw(1:nmax)/
     &      (1D0+EXP(-(e(1:nmax)-eferm-cutfac)/(-0.1D0*cutfac)))
        END IF
        acc = SUM(pairgw(1:nmax)*ph(1:nmax))
!        WRITE(6,'(a,2i4,1pg13.5)') 'pairing space=',iterin,nmax,acc
!        WRITE(6,'(a,1pg13.5)') 'PAIRGW updated before: efcut=',efcut
!        WRITE(6,*) '2. s.p.energies:'
!        WRITE(6,'(5(1pg13.5))') e(1:nmax)
!        WRITE(6,*) 'pairing space weights:'
!        WRITE(6,'(5(1pg13.5))') pairgw(1:nmax)
      else
       do n=1,nmax
        pairgw(n) = pwsav(n) ! 1.0D0
       enddo
      endif
      IF(tstcut) THEN
        finnum=0D0
        do n=1,nmax
          finnum = finnum + (pairgw(n)-pairgw(n)**2)*ph(n)
        enddo
!       	WRITE(6,*) 'nr. states in cutoff:',actnum,finnum
      END IF

c
c      prepare fermi-energy 'eferm', if not yet initialised.
c
      if(eferm.eq.zero) then
       icount = 0
       partsm = zero
       do it=1,nmax
        itsum  = it
        partsm = 2.0*ph(it) + partsm
        if(partsm-partnm.ge.zero) goto 110
       enddo
       stop ' not enough states too reach particle number'
  110  continue
       eferm  = 0.5*(e(itsum)+e(itsum+1))
       if(iab.ge.1) write(06,'(1x,a,2g11.4,a,g11.4)')
     &   'initialization: partsm,partnm=',partsm,partnm,' eferm=',eferm
       do n=1,nmax
         if(n.le.itsum) then
           gw(n) = one
         else
           gw(n) = zero
         endif
       enddo
c
      endif
c
      IF(iab>11) THEN
        WRITE(6,'(a)') 'PARTNM0:'
        DO ie=5,10
          elam= (ie-15)*1D0
          WRITE(6,'(2(1pg13.5))') elam,
     &       parnmt0(E,DELTA,PH,v2occ,PAIRGW,elam,NDIM,NMAX,cutfac)
!     &        parnmt0(E,DELTA,PH,v2occ,PAIRGW,elam,NDIM,NMAX,cutfac)
        END DO
      END IF
c
c     the coupled gradient iteration,
c     misusing 'gw' as workspace
c
C#########################
C     CR 03-4-1998
C
!      IF (RTMP.EQ.ZERO) THEN
      IF (RTMP.LT.ZERO) THEN   ! old iteration scheme
c
         if(.not.TPOLD) then
            efermc = 0D0
!            call mixed(E,DELTA,PH,v2occ,PAIRGW,EFERM,
!     &                 NDIM,NMAX,partnm,parnum,50,cutfac,efermc)
           call bisection(E,DELTA,PH,v2occ,PAIRGW,EFERM,
     &                     NDIM,NMAX,partnm,parnum,50,cutfac,efermc)
!     &                     NDIM,NMAX,partnm,parnum,10,cutfac,efermc)
            efermsav=eferm
            error  = parnum-partnm
            IF(iab>-1) WRITE(6,'(a,5(1pg13.5))') 
     &            '1. BISECT: eferm,efermc,error=',eferm,efermc,error
!            IF(ABS(error)<eps) GOTO 299
            GOTO 299
         endif
C
         DO  it=1,iter
c
c     
            elam   = eferm
            elam2  = eferm+delefe ! varied to compute derivative
            parnum = zero
            parnm2 = zero
C  
            if(TPOLD) then   
              DO na=1,nmax
                 edif   = e(na)-elam
                 delta2 = delta(na)*pairgw(na)
ctest                 delta2 = delta(na)
                 delta2 = delta2*delta2
                 v2occ(na) = half-half*edif/sqrt(edif*edif+delta2)
                 parnum = 2.0*v2occ(na)*ph(na)    + parnum
                 edif2  = e(na)-elam2
                 parnm2 = 2.0*ph(na)
     &              *(half-half*edif2/sqrt(edif2*edif2+delta2))
     &              + parnm2
              parnumsav = parnum
              parnm2sav = parnm2
              ENDDO
            else
              parnm2 = parnmt0(E, DELTA, PH, V2OCC, PAIRGW,
     &                         elam2, NDIM, NMAX, cutfac)
              parnum = parnmt0(E, DELTA, PH, V2OCC, PAIRGW,
     &                         elam, NDIM, NMAX, cutfac)
            endif
c            write(6,'(a,i5,4(1pg13.5))')
c     &       ' iter,parnmums:',iter,parnum,parnumsav,
c     &         parnm2-parnum,parnm2sav-parnumsav
C
            slope  = (parnm2-parnum)/delefe
            slope  = sign(max(abs(slope),slpmin),slope)
            error  = parnum-partnm
            eferm  = eferm - stepl*error/slope
c
            if(iab.ge.0)
     &           write(6,'(a,i5,a,5g12.4)')
     &           ' iter=',it,' eferm,parnum,parnm2,slope,error_partnm=',
     &           eferm,parnum,parnm2,slope,error
            if(iab.ge.2)
     &           write(6,'(a,20(/1x,5g12.4))')
     &           ' occupation v**2:',(v2occ(na),na=1,nmax)
            if(abs(error).lt.eps) goto 299
c
         ENDDO
         if(.not.TPOLD) then
            eferm = MIN(efermsav,efermc)
            eferm2 = MAX(efermsav,efermc)
            efermc=eferm2
            IF(iab.GE.-1) WRITE(6,'(a,3(1pg13.5))')
     &            'BISECTION: before=',eferm,efermc
            call bisection(E,DELTA,PH,V2OCC,PAIRGW,EFERM,
     &                     NDIM,NMAX,partnm,parnum,50,cutfac,efermc)
            error  = parnum-partnm
            IF(iab.GE.-1 .AND. ABS(error)>eps) WRITE(6,'(a,3(1pg13.5))')
     &            'BISECTION: efermi,err=',eferm,efermc,error
         endif
         if(iab.ge.-1) then
            IF(TPOLD) THEN
              WRITE(6,'(2(a,g11.3),1(a,i3))')
     &          ' ---> no convergence in PAIRDN!  err(part.nr)=',error,
     &          '  nz=',partnm,'  nmzx=',nmax
            ELSE
              parnumact = SUM(v2occ(1:nmax)*ph(1:nmax))*2D0
              IF(ABS(parnumact-partnm)>eps) WRITE(6,'(a,5(1pg13.5))')
     &          'BISECTION mismatch:',partnm,parnumact-partnm,
     &           eferm,efermc,efermc-eferm
            END IF
            if(abs(error).gt.eps) then
               write(6,'(2(a,g11.3),1(a,i3))')
     &          ' ---> no convergence in BISECTION!  err(part.nr)=',
     &          error,'  nz=',partnm,'  nmzx=',nmax
!               write(6,'(a,20(/1x,5g12.4))')
!     &             ' s.p. energy v**2:',(e(na),na=1,nmax)
!               write(6,'(a,20(/1x,5g12.4))')
!     &             ' occupation v**2:',(v2occ(na),na=1,nmax)
            endif
         endif
c         return
c     
 299     continue

      ELSE IF (RTMP.EQ.ZERO) THEN   ! new case T=0
        IF(tlipnog) THEN
!          IF(ioptsk==13) STOP 'LN not yet for Fayans pairing'
          CALL LN_iterate(e,delta,ph,v2occ,pairgw,partnm,nmax,ndim,
     &                    ABS(iter),iab,eferm,parnum,isosh)
         WRITE(*,*) 'after LN: ',partnm,parnum,eferm
        ELSE
          EFERM = RBRENT0(E, DELTA, PH, V2OCC, PAIRGW, PARTNM, 
     &                    NDIM, NMAX, IAB,cutfac)
          parnum = parnmt0(E, DELTA, PH, V2OCC, PAIRGW,
     &                     EFERM, NDIM, NMAX,cutfac)
!         WRITE(*,*) 'test particle number: ',partnm,parnum
C
          IT = 0               ! iteration indicator obsolete
        END IF
      ELSE
C
         EFERM = RBRENT(E, DELTA, PH, V2OCC, PAIRGW, RTMP, PARTNM, 
     &                  NDIM, NMAX, IAB)
         parnum = parnmt(E, DELTA, PH, V2OCC, PAIRGW, RTMP, 
     &                   EFERM, NDIM, NMAX)
C
         IT = 0
C
      ENDIF
C########################
c
c     mix old and new weights
c
!        write(6,'(a,200(1pg13.5))') 'V2SAV:',v2sav(1:nmax)
!        write(6,'(a,200(1pg13.5))') 'V2OCC:',v2occ(1:nmax)
      if(ifsav.eq.1 .and. icount.ge.100 .and. addmix.lt.one) then
        do n=1,nmax
c          addmix = (gwsav(n)-gwsav(n)*gwsav(n))**0.5
          v2occ(n) = addmix*v2occ(n) + (1.0-addmix)*v2sav(n)
        enddo
      endif
c
c     update pairing space
c
      if(tphpair .AND. ifcut.eq.1 .AND. cutfac<0D0) then
        IF(ifkincut==1) THEN
!          efkincut=efermkin-cutfac
          efadd=-cutfac
!          edirec(1:nmax)=esav(isosh:isosh+nmax-1)
          edirec(1:nmax)=e(1:nmax)   ! restore s.p. energy
          efkincut=eferm-cutfac
        ELSE
          edirec(1:nmax)=e(1:nmax)   ! restore s.p. energy
          efcut=eferm-cutfac
          efadd=-cutfac
        END IF
!        do n=1,nmax
        IF(ifkincut==1) THEN
          pairgw(1:nmax) = 
     &      1D0/(1D0+EXP((edirec(1:nmax)-efkincut)/(-0.1D0*cutfac)))
        ELSE
          pairgw(1:nmax) = 
     &      1D0/(1D0+EXP((e(1:nmax)-eferm+cutfac)/(0.100D0*efadd)))
!     &      1D0/(1D0+EXP((edirec(1:nmax)-efcut)/(0.100D0*efadd)))
          IF(tlowercut) pairgw(1:nmax) = pairgw(1:nmax)/
     &      (1D0+EXP(-(e(1:nmax)-eferm-cutfac)/(0.100D0*efadd)))
!     &      (1D0+EXP(-(edirec(1:nmax)-eferm-cutfac)/(0.100D0*efadd)))
         END IF
         IF(iab>-1) THEN
           WRITE(6,*) '2. s.p.energies:'
           WRITE(6,'(5(1pg13.5))') edirec(1:nmax)
           WRITE(6,*) 'pairing space weights:'
           WRITE(6,'(5(1pg13.5))') pairgw(1:nmax)
!        enddo
           WRITE(6,'(a,1pg13.5)') 'PAIRGW updated after: efcut=',efcut
         END IF
      endif
c
c     finally compute average gap and pairing energy
c
      sumuv  = 0D0
      sumduv = 0D0
      sumv2  = 0D0
      sumdv2 = 0D0
      eqmin  = 10000D0
      sumdev = 0D0
C
      IF (RTMP.EQ.ZERO) THEN
         DO  na=1,nmax
            IF(tphpair) THEN
              deltaa = delta(na)*pairgw(na)   !  old cutoff
            ELSE
              deltaa = delta(na)  ! new cutoff
            END IF
!            deltab = delta(na)
!            edif   = e(na)-eferm
            edif   = e(na)-eferm
            equasi = sqrt(edif*edif+deltaa*deltaa)
            eqmin  = min(eqmin,equasi)
            v2     = (half-half*edif/equasi)
            sumdev = (v2-v2occ(na))**2*ph(na) + sumdev
            uv     = SQRT(MAX(v2-v2*v2,xsmall)) 
            IF(.NOT.tphpair) THEN
              v2     = v2*pairgw(na)   ! new cutoff
              uv     = uv*pairgw(na)   ! new cutoff
            END IF
            vol    = ph(na)*uv
            sumuv  = vol        + sumuv
            sumduv = vol*deltaa + sumduv
!            sumduv = vol*deltab + sumduv
!            sumepa = ph(na)*half*(deltaa*deltaa)/equasi + sumepa
            vol    = ph(na)*v2
            sumv2  = vol        + sumv2
            sumdv2 = vol*deltaa + sumdv2
         ENDDO
      ELSE
         DO  na=1,nmax
            deltaa = delta(na)*pairgw(na)
!            edif   = e(na)-eferm
            edif   = e(na)-eferm
            equasi = sqrt(edif*edif+deltaa*deltaa)
            v2     = half-half*edif/equasi
            TEXP   = 1D0 / (EXP(EQUASI/RTMP)+1D0)
            TTAN   = 1D0 - 2D0 * TEXP
            wgthrm(na) = TTAN
            vol    = ph(na)*sqrt(max(v2-v2*v2,xsmall))*TTAN
            sumuv  = vol        + sumuv
            sumduv = vol*deltaa + sumduv
!            sumepa = ph(na)*half*(deltaa*deltaa*TTAN*TTAN)/equasi 
!     &           + sumepa
         ENDDO
      ENDIF
C
      sumdev = SQRT(MAX(sumdev/nmax,1D-40))
      sumuv  = max(sumuv,xsmall)
      avdelt = sumduv/sumuv
      avdlt2 = sumdv2/max(sumv2,xsmall)
      epair  = -sumduv
      avg    = epair/sumuv**2
      if(iab.ge.-1)
     &   write(6,'(a,i4,a,f8.2,a,2f8.3,4f6.2,1pg13.5)')
     &  ' PAIRDN: it.=',it,' #part.=',parnum,
     &  ' eferm,epair,gaps=',
     &    eferm,epair,avdelt,avdlt2,eqmin,sumuv,sumdev
      IF(sumdev>eps) THEN
        WRITE(6,'(a,5(1pg13.5))') 'PAIRDN: V2 deviates:',
     &   sumdev,sumv2*2D0,eferm,efermc,efermc-eferm     
!        WRITE(6,'(a,i5)') 'V2occ:',isosh
!        WRITE(6,'(6(1pg13.5))') V2occ(1:nmax)
      ENDIF
      if(iab.gt.0) write(6,'(5x,a,3g12.4)')
     &   ' av.Delta,av.G,pair.E=',avdelt,avg,epair
!
!     recoover from blocking
!
      gw(1:nmax)=v2occ(1:nmax)
      IF(nblck>0) THEN
        partnm = partnm + 1D0                     ! add blocked particle
        ph(nblck) = ph(nblck) + 1D0
!        gw(nblck) = gw(nblck)*(ph(nblck)-1D0)/ph(nblck)
!     &              +0.5D0/ph(nblck)
        gw(nblck) = (v2occ(nblck)*(ph(nblck)-1D0)+0.5D0)/ph(nblck)
        IF(iab>0) WRITE(6,'(a,1pg13.5)') ' actual part.nr.=',
     &    SUM(ph(1:nmax)*gw(1:nmax))*2D0
      END IF
c
      IF(tswitchexplicit) THEN
        IF(isosh>1) WRITE(6,'(a,3(1pg13.5))') 
     &      'PART.NR.:',partnm,parnum,parnum + v2explicit*2D0
        partnm = partnm + v2explicit*2D0                     ! add blocked particle
        ph(nexplicit) = ph(nexplicit) + 1D0
!        gw(nexplicit) = gw(nexplicit)*(ph(nexplicit)-1D0)/ph(nexplicit)
!     &              +0.5D0/ph(nexplicit)
        gw(nexplicit) = (v2occ(nexplicit)*(ph(nexplicit)-1D0)
     &                   +v2explicit)/ph(nexplicit)
!        IF(tswitchexplicit) WRITE(6,'(a,3(1pg13.5))')
!     &    'GW(NEXPLICIT)=',gw(nexplicit),v2occ(nexplicit),ph(nexplicit)
!        IF(iab>0) WRITE(6,'(a,1pg13.5)') ' actual part.nr.=',
        v2occ=gw        
        sumduv=0D0
        DO  na=1,nmax
           IF(tphpair) THEN
             deltaa = delta(na)*pairgw(na)   !  old cutoff
           ELSE
             deltaa = delta(na)  ! new cutoff
           END IF
!           edif   = e(na)-eferm
!           equasi = sqrt(edif*edif+deltaa*deltaa)
!           eqmin  = min(eqmin,equasi)
           v2 = v2occ(na)
!           v2     = (half-half*edif/equasi)
           uv     = SQRT(MAX(v2-v2*v2,xsmall)) 
           IF(.NOT.tphpair) THEN
             v2     = v2*pairgw(na)   ! new cutoff
             uv     = uv*pairgw(na)   ! new cutoff
           END IF
            vol    = ph(na)*uv
            sumuv  = vol        + sumuv
            sumduv = vol*deltaa + sumduv
            vol    = ph(na)*v2
            sumv2  = vol        + sumv2
            sumdv2 = vol*deltaa + sumdv2
         ENDDO
         epair=-sumduv
        WRITE(6,'(a,2(1pg13.5))') ' actual part.nr.,epair=',
     &    SUM(ph(1:nmax)*gw(1:nmax))*2D0,epair
      END IF
      edirec=esav
!      if(iab.gt.-1) then
!        write(6,'(/a,i4,(a,f9.4),a,i4/a,10(/1x,10f8.3))')
!     &   ' *** pairing with dynamic gap exited with: iter=',iter,
!     &    ' cutfac=',cutfac
!        write(6,'(a,10(/1x,10f8.4))') 'pairgw:',(pairgw(i),i=1,nmax)
!        write(6,'(a,10(/1x,10f8.4))') 'pairgw:',
!     &    (1D0/(1D0+EXP((edirec(isosh+i)+cutfac)/(-0.05D0*cutfac)))
!     &     ,i=0,nmax-1)
!        write(6,'(a,10(/1x,10f8.4))') 
!     &         ' spkin:',(edirec(isosh+i),i=0,nmax-1)
!      endif
      IF(tswitchexplicit) nexplicit=nexplicit+isosh-1
      return
      end
c-----lagpar------------ part of the hartree-fock package ----------------------
c
      subroutine lagpar(iter)
c
      include 'params.inc'
c
c     determines the matrix of Lagrange parameters for HFB iteration
c     and stores result on PAIRMT on common /shells/.
c     the action of the HFB Hamiltonian and the wavefunction is
c     stored on the workspace SAVLR from common /linres/
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'linres.inc'
      include 'constn.inc'
c
      INTEGER,INTENT(IN) :: iter

      REAL(8), ALLOCATABLE :: symatrix(:,:)
      REAL(8), ALLOCATABLE :: workt1(:),workt2(:)
      LOGICAL,PARAMETER :: tblocknew=.FALSE.  ! corrected blocking factor
c
      dimension savpri(KNMX)
!      REAL(8) :: symmcond(2)
c
      common /ctest/edirec(ksl2),efkincut,ifkincut
!      REAL(8), PARAMETER :: pairoffset=xsmall
c
c      data xminv/1.0E-05/
c
c-------------------------------------------------------------------------------
c
c     reset matrix of Lagrange parameters
c
      do na=1,nmaxp
        do nb=1,nmaxp
          pairmp(na,nb) = zero
        enddo
      enddo
      do na=1,nmaxn
        do nb=1,nmaxn
          pairmn(na,nb) = zero
        enddo
      enddo
c     iso-loop
c      'isoadd' is a pointer for the potential and mass fields;
c      'isnpad' is a pointer for the wave function pointer 'npla';
c      on 'egx' is accumulated lateron the shell-energy.
c
      ALLOCATE(workt1(1:ngrid),workt2(1:ngrid))
      do 1000 iso=1,2
      if(iso.eq.1) then
        isoadd = 0
        isnpad = 0
        eferm  = efep
        nmaxa=nmaxp
        nnumb=nprot
        nblock=nblockp
      else
        isoadd = ngrid
        isnpad = kjmx
        eferm  = efen
        nmaxa=nmaxn
        nnumb=nneut
        nblock = nblockn !-nmaxp
      endif
      ALLOCATE(symatrix(nmaxa,nmaxa))
      symatrix=0D0
c
c     square root of inverse free mass on 'sqb';
c
      call derv15(bmass(isoadd+1),d1b(isoadd+1),one,ngrid,rstep)
c
c  FIRST BIG LOOP: compute HFB-Hamiltonian on wavefunction and
c                  matrix of Lagrange parameters
c
c     loop over composed ang.mom. 'js'
c
      do 40 js=1,jsmax(iso)
      jsnpl = js+isnpad
      nmax  = npla(1,jsnpl)
      if(nmax.gt.0) then
c
c      loop over radial quantum number 'n'.
c       accumulating matrix elements with HFB Hamiltonian (mix of
c       mean field and gap potential)
c       'npl' is wave function pointer,
c       'nplshl' is number of state in harmonic oscillator ordering.
c
       do 70 n=1,nmax
       npl    = npla(n+1,jsnpl)
       nplm   = npl-1
       nplshl = npl/ngrid+1
c
c      action of the mean-field hamiltonian
c
       call hpsi(wf(npl),aux3,wfs(npl),aux2,0,js,iso,zero)
       IF(ifkincut==1) THEN
!         call tpsi(wf(npl),aux2,wfs(npl),aux2,0,js,iso,zero)
         aux2=aux3
         fackincut = (1D0-phpair(nplshl))*deltaf(nplshl)/(-0.1D0*cutfac)
       END IF
!       savlr(npl:nplm+ngrid,2)=wf(npl:nplm+ngrid)
!       savlr(npl:nplm+ngrid,3)=aux3(1:ngrid)
ctest       if(iwr.ge.0) then
       if(iwr.ge.-2) then
         testsp = zero
         do ir=2,ngrid
           testsp   = aux3(ir)*wf(nplm+ir)    + testsp
ctest           testsp   = ((ir-1)*rstep*wf(nplm+ir))**2    + testsp
         enddo
         delcm(nplshl) = testsp*rstep
!         edirec(nplshl) = testsp*rstep
       endif
c
c      complement to actual HFB Hamiltonian
c
       v2  = max(weight(nplshl),xsmall)
       uv     = sqrt(max(v2-v2*v2,1D-200)) !xsmall))   
       IF(nplshl==nblock) THEN
         v2b = v2occ(nplshl)
!         v2b = v2-0.5D0/ph(nblock)
!         v2b = v2b*ph(nblock)/(ph(nblock)-1D0)
         uv  = sqrt(max(v2b-v2b*v2b,1D-200))
       END IF
!       IF(v2b==0D0 .AND. uv>xsmall) THEN
!         WRITE(11,'(a,i4,5(1pg13.5))') ' strange blocking:',
!     &     nplshl,v2b,weight(nplshl),uv
!       END IF
       IF(nplshl==nblock .AND. tblocknew) THEN
         uv = uv*(ph(nblock)-1D0)/ph(nblock)
       END IF
       IF(tphpair) THEN
         uv     = uv * phpair(nplshl) ! old cutoff
       END IF
!       v2 = MAX(weight(nplshl),xsmall)
!       v2 = weight(nplshl)
       accel  = v2**(-atajima)
       v2a    = v2*accel
       uva    = uv*accel
       IF(thfbfold) THEN
         IF(ifkincut==1) STOP "kinetic cutoff incompatible with folding"
!         CALL gapoperator(wffold(nplm+1:nplm+ngrid),aux2(1:ngrid),
!     &                    rstep,ngrid,js,pairpt(isoadd+1:isoadd+ngrid))
!         aux1(1:ngrid)=pairpt(isoadd+1:isoadd+ngrid)*
!     &                 wffold(nplm+1:nplm+ngrid)
         workt1(1:ngrid)=pairpt(isoadd+1:isoadd+ngrid)*
     &                 wffold(nplm+1:nplm+ngrid)
!         WRITE(6,'(a,4i4)') 'before Call Gaussian:',js,jsmax(iso),nplshl
         CALL Gauss_folding_m(workt1(1:ngrid),workt2(1:ngrid),
     &                       rstep,ngrid,js,hfbwidth,MAXVAL(jsmax(1:2)))
!         CALL Gauss_folding_m(aux1(1:ngrid),aux2(1:ngrid),rstep,ngrid,js,
!     &                        hfbwidth,jsmax(iso))
!         wfdiff=rstep*SUM((aux2(1:ngrid)-aux1(1:ngrid))**2)
!         IF(wfdiff>1D-5 .AND. uva>1D-3) THEN
!           WRITE(6,'(a,2i4,2(1pg13.5))') 'WFDIFF:',js,nplshl,uva,wfdiff
!         END IF
         do ir=1,ngrid
           aux3(ir) = v2a*aux3(ir)+uva*workt2(ir)  !aux2(ir)
           savlr(nplm+ir,1) = aux3(ir)                    ! save action of HFB
         enddo
       ELSE
         do ir=1,ngrid
cold         aux3(ir) = v2*aux3(ir)+uv*pairpt(ir+isoadd)*wf(nplm+ir)
!         IF(nplshl.NE.nblock) 
           aux3(ir) = v2a*aux3(ir)+uva*pairpt(ir+isoadd)*wf(nplm+ir)
           IF(ifkincut==1) aux3(ir) = aux3(ir) +uva*fackincut*aux2(ir)
           savlr(nplm+ir,1) = aux3(ir)                    ! save action of HFB
         enddo
       END IF
       if(iwr.ge.0) then
         testsp = zero
         do ir=1,ngrid
           testsp   = aux3(ir)*wf(nplm+ir)    + testsp
         enddo
         drkin(nplshl) = rstep*testsp
       endif
c
c      accumulate to matrix of Lagrange parameters, thereby symmetrising
c
       do 71 nb=1,nmax
       nplb   = npla(nb+1,jsnpl)
       nplmb  = nplb-1
       nplshb = nplb/ngrid+1
       acc    = zero
       do ir=2,ngrid
         acc    = wf(nplmb+ir)*aux3(ir) + acc
       enddo
cold       acc    = 0.5*rstep*acc
       acc    = rstep*acc
       facmix = 1.0/max(weight(nplshb)**(atajima),xsmall)
     &         +1.0/max(weight(nplshl)**(atajima),xsmall)
       if(iso.eq.1) then
         pairmp(nplshb,nplshl) = acc/facmix + pairmp(nplshb,nplshl)
         pairmp(nplshl,nplshb) = acc/facmix + pairmp(nplshl,nplshb)
         symatrix(nplshb,nplshl) = acc
       else
         nxa = nplshl-nmaxp
         nxb = nplshb-nmaxp
         pairmn(nxb,nxa) = acc/facmix + pairmn(nxb,nxa)
         pairmn(nxa,nxb) = acc/facmix + pairmn(nxa,nxb)
         symatrix(nxb,nxa) = acc
       endif
   71  continue
c
   70  continue
      endif
   40 continue
c
c  END OF FIRST BIG LOOP
c
c
c
c  SECONG BIG LOOP: subract Lagrange parameters from stored wavefunctions
c
c
c     loop over composed ang.mom. 'js'
c
      do 140 js=1,jsmax(iso)
      jsnpl = js+isnpad
      nmax  = npla(1,jsnpl)
      if(nmax.gt.0) then
c
c      double loop over radial quantum number 'n', 'nb'.
c
       if(iwr.gt.0)  write(6,'(a,2i5)')
     &   ' Lagrange matrix for: iso,js=',iso,js
       do 170 n=1,nmax
       npl    = npla(n+1,jsnpl)
       nplm   = npl-1
       nplshl = npl/ngrid+1
c
       do 171 nb=1,nmax
       nplb   = npla(nb+1,jsnpl)
       nplmb  = nplb-1
       nplshb = nplb/ngrid+1
       if(iso.eq.1) then
         fac    = pairmp(nplshb,nplshl)
       else
         fac    = pairmn(nplshb-nmaxp,nplshl-nmaxp)
       endif
       if(iwr.gt.0) savpri(nb) = fac
       do ir=2,ngrid
         savlr(nplm+ir,1) = savlr(nplm+ir,1) - fac*wf(nplmb+ir)
       enddo
  171  continue
       if(iwr.gt.0) write(6,'(1x,5(1pg12.4))') (savpri(nb),nb=1,nmax)
cold       vpact = 1.0/max(sqrt(weight(nplshl)),xminv)
cold       do ir=2,ngrid
cold         savlr(nplm+ir,1) = vpact*savlr(nplm+ir,1)
cold       enddo
c
  170  continue
      endif
  140 continue
c
c  END OF SECOND BIG LOOP
c
c    test prints
c
      if(iwr.ge.0) then
        if(iso.eq.1) then
         npr1 = 1
         npr2 = nmaxp
        else
         npr1 = nmaxp+1
         npr2 = nmaxn+nmaxp
        endif
        write(6,'(a,10(/1x,10f8.3))')
     &        ' spe direct:',(delcm(i),i=npr1,npr2)
        write(6,'(a,10(/1x,10f8.3))') ' spe HFB:',(drkin(i),i=npr1,npr2)
      endif
c
c     compute symmetry condition
c
      symmcond(iso)=0D0      
      DO na=1,nmaxa
      DO nb=1,nmaxa
        symmcond(iso) = symmcond(iso) + 
     &    (symatrix(na,nb)-symatrix(nb,na))**2
      END DO
      END DO
      symmcond(iso) = SQRT(symmcond(iso)/nnumb)
!      IF(iso==2) THEN
!        OPEN(97,file='for097')
!        DO na=1,nmaxa
!          WRITE(97,'(i3,50(1pg13.5))') na,symatrix(na,:)-symatrix(:,na)
!        END DO
!        CLOSE(97)
!      END IF
      DEALLOCATE(symatrix)
 1000 continue
      DEALLOCATE(workt1,workt2)
c
!      WRITE(6,'(a,i5,2(1pg13.5))') 'symmcond:',iter,symmcond
!      WRITE(96,'(a,i5,2(1pg13.5))') 'HFB symmcond:',iter,symmcond
      return
      end
c-----cutoff------------ part of the hartree-fock package ----------------------
c
      subroutine cutoff(e,gw,ph,anz,nmax,eferm,efadd,delta,
     &                  partnm,iter,eps,iab)
c
c     this routine determines the effective "Fermi-energy" and
c     cut-off distribution for the pairing space.
c
c     the include file 'params.inc' communicates the dimensioning of
c     the fields 'e' and 'gw' with the parameter 'ksl2'. it allows
c     furthermore to make the whole subroutine double precision by
c     transfering a line 'implicit double precision (a-h,o-z)'.
c
      include 'params.inc'
c
c     input via list:
c      e       = array(1..nmax) of single particle energies.
c      gw      = array(1..nmax) of cut-off weights,
c                input is used as initial guess for the iteration.
c      ph      = array(1..nmax) of the multiplicities of the single
c                particle levels, to be provided by the calling routine.
c      anz     = desired effective particle number
c      nmax    = maximum number of single particle states
c                (needs to be smaller than the parameter 'ksl2'
c                 and smaller than 46, for dimensioning reasons).
c      eferm   = true fermi-energy as given by the pairing routine.
c      efadd   = additional "fermi-energy" for the pairing-space weights
c      delta   = relative width of the distribution
c      iter    = number of pairing iterations
c      eps     = end condition for iteration;
c                iteration terminates if error in effective particle number
c                is smaller than 'eps' and if 'eferm' and/or 'delta'
c                are smaller than 'eps'.
c      iab     = print level; significant values are -2, -1, 0 and 1;
c                larger value prints more.
c
c     output via list:
c      gw      = array(1..nmax) of the occupation probabilities
c                finally achieved.
c      efadd   = additional "fermi-energy" for the pairing-space weights
c      partnm  = particle number, from adding 'gw' with degeneracy 'ph'.
c
c     some auxiliary fields are:
c      dparn   = array(1..3) of variation in the particle number as
c                consequence of 'deferm' and 'ddelta'.
c      dgap    = array(1..3) of variation in the gap equation as
c                consequence of 'deferm' and 'ddelta'.
c      gr      = double array of variation of particle number and of
c                gap equation, used for construction of the tangantial
c                step
c
      dimension e(ksl2),gw(ksl2),ph(ksl2)
      dimension dparn(2),gr(1,1)
c
cold      data ph/1.,2.,1.,3.,2.,1.,4.,3.,2.,1.,5.,4.,3.,2.,1.,6.,5.,
cold     &        4.,3.,2.,1.,7.,6.,5.,4.,3.,2.,1.,8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        9.,8.,7.,6.,5.,4.,3.,2.,1.,10./
      data gr,dparn/3*0./
c
c     further numerical variables are:
c      deferm  = step in the fermi-energy, to evaluate
c                numerically the derivatives.
c      ddelta  = step in the gap, to evaluate numerically the
c                derivatives; used only in the "constant force" approach.
c      dampin  = damping factor for the tangential step
c      efmax,dlmax,detmin,defmin,ddlmin  = limitations for the step
c
c     note: these parameters are adjusted for an application with
c           the nuclear shell model. it is assumed that the single
c           particle energies are given in units of [mev] and display
c           a reasonable nuclear spectrum.
c
      data dampin,deferm/0.9,1.0/
cold      data efmax,detmin,defmin/3.,.1e-3,.1e-3/
      data efmax,detmin/3.,.1e-3/
c
c      data xmaxlw,x10lw/1.0e-30,1.0e-10/
      data one,zero,fifty/1.0,0.0,50.0/
      data iftang/0/
c
c-------------------------------------------------------------------------------
c
c      set switch variables
c
      iwarng = 0
c
      if(iab.ge.0)
     & write(6,'(/a,2i4,a,4f9.4)')
     &   ' *** cutoff entered with: iter,nmax=',iter,nmax,
     &   '  delta,efadd,eferm,partnm=',delta,efadd,eferm,anz
      if(iab.ge.1)
     & write(6,'(a,10(/1x,10f8.3))')
     &   ' single particle energies:',(e(i),i=1,nmax)
c
c      prepare fermi-energy 'eferm', if not yet initialised.
c
      if(efadd.eq.zero) then
       isum   = 0
       do it=1,nmax
        itsum  = it
        isum   = int(ph(it)+.000001) +isum
        if((one+one)*isum-anz.ge.0) goto 110
       enddo
       stop ' not enough states too reach particle number'
  110  continue
       efadd  = e(itsum)+deferm-eferm
       if(iab.ge.1) write(06,'(5x,a,i5,a,g11.4)')
     &   'initialistion: isum=',isum,'  efadd=',efadd
c
      endif
c
c
c-------------------------------------------------------------------------------
c
c     iterate Fermi energy and the distribution
c
c-------------------------------------------------------------------------------
c
c     tangent step
c
      if(iftang.eq.1) then
        do it=1,iter
c
c       particle number at 'eferm' (k=1) and 'eferm' varied (k=2).
c
         do ka=1,2
c                                        set actual fermi-energy 'elam'
          k      = 3-ka
          deltac  = efadd*delta
          if(k.eq.2) then
           elam   = eferm+efadd+deferm*deltac
          else
           elam   = eferm+efadd
          endif
c                                         compute actual particle number
          parnum = zero
          do i=1,nmax
           arg    = (e(i)-elam)/deltac
           if(arg.gt.fifty) then
             gw(i)   = zero
           elseif(arg.lt.-fifty) then
            gw(i)   = one
           else
             gw(i)  = 1.0/(1.0+exp(arg))
           endif
           parnum = gw(i)*ph(i)    + parnum
          enddo
c                                         store particle number
          dparn(k) = parnum+parnum-anz
         enddo
c
c      gradient step to next 'eferm'   (including step limiters)
c
         gr(1,1) = (dparn(2)-dparn(1))/(deferm*deltac)
         det    = sign(max(abs(gr(1,1)),detmin),gr(1,1))
c         efold  = efadd
         efstep = dampin*dparn(1)/det
         efadd  = efadd-sign(min(abs(efstep),efmax),efstep)
c         deferm = efold-efadd
c         deferm = sign(max(abs(deferm),defmin),deferm)
         if(iab.eq.0)
     &    write(6,'(a,i4,2(a,g11.3))')
     &     ' pairing-weight iteration it=',it,'  :  efadd=',efadd,
     &     '  dparn(1)=',dparn(1)
         if(iab.ge.1)
     &    write(6,'(a,i4,2(a,g11.3)/a,2g11.3,a,g11.3/
     &              a,10f7.4,10(/6x,10f7.4))')
     &     ' iteration at it=',it,' :   efadd=',efadd,
     &     '  det=',det,
     &     ' dparn=',dparn(1),dparn(2),' ;   gr=',gr(1,1),
     &     '   gw=',(gw(n),n=1,nmax)
         if(abs(dparn(1)).lt.eps)   goto 350
        enddo
        iwarng = 1
  350   continue
      else
c
c     secant step
c
c       eferm1 = eferm                          !   +(eferm-e(1))/200.0
       eferm1 = 0.0D0
       eferm2 = e(1)
       do i=2,nmax
         eferm1 = min(e(i),eferm1)
         eferm2 = max(e(i),eferm2)
       enddo
       do it=1,iter
         if(it.eq.1) then
           eferm3 = eferm1
         elseif(it.eq.2) then
           eferm3 = eferm2
         else
           eferm3 = eferm1 - dparn1*(eferm2-eferm1)/(dparn2-dparn1)
         endif
         efadd   = eferm3-eferm
         deltac  = max(efadd*delta,delta)
         parnum  = 0.0
         do i=1,nmax
          arg    = (e(i)-eferm3)/deltac
          if(arg.gt.fifty) then
            gw(i)   = zero
          elseif(arg.lt.-fifty) then
            gw(i)   = one
          else
            gw(i)  = 1.0/(1.0+exp(arg))
          endif
          parnum = gw(i)*ph(i)    + parnum
         enddo
         dparnm = parnum+parnum-anz
         if(it.eq.1) then
           dparn1  = dparnm
           if(dparn1.gt.zero) then
             write(6,'(a,3(1pg13.5))')
     &         ' eferm1,dparnm,e(1)=', eferm1,dparnm,e(1)
             write(11,'(a,3(1pg13.5))')
     &         ' eferm1,dparnm,e(1)=', eferm1,dparnm,e(1)
             stop ' lower state too high in CUTOFF'
           endif
         elseif(it.eq.2) then
           dparn2  = dparnm
           if(dparn2.lt.zero) stop ' upper state too low in CUTOFF'
         else
           dparn3  = dparnm
         endif
         if(iab.ge.0)
     &    write(6,'(a,i4,2(a,g11.3))')
     &    ' secant iteration it=',it,'  :  eferm=',eferm3,
     &    '  dparn3=',dparn3
         if(it.gt.2) then
           if(abs(dparn3).lt.eps) then
             goto 399
           elseif(dparn3.gt.zero) then
             eferm2 = eferm3
             dparn2 = dparn3
           else
             eferm1 = eferm3
             dparn1 = dparn3
           endif
         endif
       enddo
  399  continue
      endif
c
c-------------------------------------------------------------------------------
c
c     finally the particle-number is checked and printed
c
c-------------------------------------------------------------------------------
c
      partnm = zero
      do i=1,nmax
       partnm = gw(i)*ph(i) + partnm
      enddo
c
      if(iab.ge.-1)
     &  write(6,'(a,i4,2(1pg13.5))')
     &   ' CUTOFF finished with iter,part.nr.,efadd=',it,2*partnm,efadd
      if(iab.ge.0)
     & write(6,'(a,10(/1x,10f8.3))')
     &   ' pairing-weights:',(gw(i),i=1,nmax)
c
c     in case of incomplete convergence print a warning
c
      if(iab.ge.-1 .and. iwarng.eq.1)
     & write(6,'(2(a,g11.3),2(a,i3))')
     &  ' ---> no convergence in cutoff! err(part.nr)=',dparn(1)
c
      return
      end
c-----sprad ------------ part of the hartree-fock package ---------------------
c
      subroutine sprad(rmssta,rmsav,gw,phin,isosh)
c
      include 'params.inc'
c
c     computes rms radii of separate s.p. states for type 'iso',
c     result is returned on 'rmssta'.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
      dimension rmssta(ksl2),gw(ksl2),phin(ksl2)
      dimension rmsmin(0:kjmx)
c
      data ifcort/0/
c
c------------------------------------------------------------------------------
c
       do l=0,kjmx
         rmsmin(l) = 10000.0
       enddo
       if(isosh.eq.1) then
         nmax = nmaxp
         nadd = 0
       else
         nmax = nmaxn
         nadd = nmaxp
       endif
       lmax = 0
       rmsav = 0.0
       avnum = 0.0
       do nn=1,nmax
         n = nn + nadd
         js=jsp(nn)
         ll=js/2
         lmax = max(ll,lmax)
         if(2*ll.eq.js) js=js-1
         npl = (n-1)*ngrid
         sum0 = 0.0
         sum2 = 0.0
         do ir=1,ngrid
           r  = (ir-1)*rstep
           r2 = r*r
c           vol= wf(npl+ir)*wf(npl+ir)*r2
           vol= wf(npl+ir)*wf(npl+ir)
           sum0 = vol + sum0
           sum2 = vol*r2 + sum2
         enddo
c         write(6,'(a,i5,2g12.4)') ' nn,sum2,sum0=',nn,sum2,sum0
         rmssta(nn) = sqrt(abs(sum2/max(sum0,xsmall)))
         vol   = gw(nn)*phin(nn)
         rmsav = vol*rmssta(nn)**2 + rmsav
         avnum = vol + avnum
         rmsmin(ll) = min(rmssta(nn),rmsmin(ll))
       enddo
       rmsav = sqrt(rmsav/avnum)        
       write(6,'(a,i5,a,20(/1x,10f8.4))') 
     &   ' maximum L=',lmax,' RMSMIN:',(rmsmin(ll),ll=0,lmax)
c
c      correct for trend with angular momentum 'l'
c
       if(ifcort.eq.1) then
         avnor = zero
         avrms = zero
         do nn=1,nmax
           n = nn + nadd
           js=jsp(nn)
           ll=js/2
           rmssta(nn) = rmssta(nn) - rmsmin(ll) + rmsmin(0)
           vol   = gw(nn)*phin(nn)
           avnor = vol + avnor
           avrms = vol*(rmssta(nn)*rmssta(nn)) + avrms
        enddo
        rmsav = sqrt(avrms/avnor)
        write(6,'(a,20(/1x,10f8.4))') ' RMS:',(rmssta(nn),nn=1,nmax)
      endif
c
      return
      end
c-----cutpot------------ part of the hartree-fock package ---------------------
c
      subroutine cutpot(ein,eout,iso)
c
      include 'params.inc'
c
c     computes expectation values of the "penalty potential" for
c     continuum states.
c     ein   = array of s.p. energies
c     eout  = resulting s.p. energies plus penalty energies 
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'result.inc'
      include 'constn.inc'
c
      dimension ein(kg1),eout(kg1)
c
      data iftest/0/
      data ifirst/1/
      data rhocut/0.0016d0/
      data potfac/-3.0d0/
c
c------------------------------------------------------------------------------
c
c  testprint log(rho) as additional potential
c
      if(ifirst.eq.0) then
        write(6,'(/a)') '# r log(rho_p)   log(rho_n)  log(rho):'
        do ir=1,ngrid
          irn = ir+ngrid
          r = (ir-1)*rstep
          rhoc2 = rhocut/2.0
          write(6,'(1x,f8.1,3(1pg13.5))')
     &     r,log(rho(ir)/rhoc2),log(rho(irn)/rhoc2),
     &     log((rho(ir)+rho(irn))/rhocut)
        enddo
        ifirst=1
      endif
c
c     loop over states and computation of averages
c
      if(isosh.eq.1) then
        nmax   = nmaxp
        npladd = 0
      else
        nmax   = nmaxn
        npladd = nmaxp
      endif
      do n=1,nmax
        js=jsp(n)
        ll=js/2
        if(2*ll.eq.js) js=js-1
        npl = (n+npladd-1)*ngrid
        sum  = 0.0
        sumn = 0.0
        do ir=1,ngrid
          r  = (ir-1)*rstep
          sumn = (wf(npl+ir)*wf(npl+ir))  + sumn
          sum  = (wf(npl+ir)*wf(npl+ir))*
     &           max(zero,potfac*log((rho(ir)+rho(ir+ngrid))/rhocut))
     &          + sum
        enddo
        eout(n) = ein(n) + sum/max(sumn,xsmall)
c        write(6,'(a,i5,1pg13.5)') ' test norm:',n,sumn
      enddo
c
c     testprint
c
      if(iftest.eq.1) then
        write(6,'(/a,i1,a/a)') 
     &   ' corrected s.p. energies for iso=',iso,':',
     &   '  state    s.p.    &penalty    diff. '
        do n=1,nmax
          write(6,'(1x,i5,3(1pg13.5))')
     &     n,ein(n),eout(n),eout(n)-ein(n)
        enddo
      endif
c
      return
      end
c-----cutang------------ part of the hartree-fock package ---------------------
c
      subroutine cutang(ein,eout,iso)
c
      include 'params.inc'
c
c     computes s.p. energies which are shifted by the
c     minimum of the effective l-dependent potential relative to 1s.
c     ein   = array of s.p. energies
c     eout  = resulting s.p. energies plus shift
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'result.inc'
      include 'constn.inc'
c
      dimension ein(kg1),eout(kg1)
c
      data iftest/0/
c
c------------------------------------------------------------------------------
c
c     loop over states and computation of averages
c
      if(isosh.eq.1) then
        nmax   = nmaxp
        npladd = 0
        isoadd = 0
      else
        nmax   = nmaxn
        npladd = nmaxp
        isoadd = ngrid
      endif
      call derv15(bmass(isoadd+1),d1b(1),one,ngrid,rstep)
      d1b(1)=(16.0*bmass(2+isoadd)-15.0*bmass(1+isoadd)
     *        -bmass(3+isoadd))/(6.0*rstep**rstep)
      do i=2,ngrid
        d1b(i) = d1b(i)*rinvs(i)
      enddo
c      write(6,'(a,i3)') ' isoadd=',isoadd
c      write(6,'(1x,4g12.4)')
c     &  (pot(i),potls(i),bmass(i),d1b(i),i=ngrid+1,2*ngrid)
      do n=1,nmax
        js=jsp(n)
        ll=js/2
        npl = (n+npladd-1)*ngrid
        potlmn = 1.0D30
        xlp    = ll*ll+ll
        xl0    = wls(js)
c        xl0    = 0.0
        do ir=2,ngrid
         j      = ir+isoadd
         potact = pot(j)+xl0*potls(j)+bmass(j)*xlp*zentrf(ir)+d1b(ir)
c         write(6,'(a,i3,g12.4)') ' ir,potact=',ir,potact
         if(potlmn.ge.potact) then
           potlmn = potact
         else
           goto 99
         endif
        enddo
   99   continue
        if(n.eq.1) then
          absmin = potlmn
          shift  = 0.0
        else
          eshift  = potlmn - absmin
        endif
        eout(n) = ein(n) + eshift
c        write(6,'(a,i5,1pg13.5)')  ' n,shift=',n,eshift
      enddo

c     testprint
c
      if(iftest.eq.1) then
        write(6,'(/a,i1,a/a)') 
     &   ' corrected s.p. energies for iso=',iso,':',
     &   '  state    s.p.    &penalty    diff. '
        do n=1,nmax
          write(6,'(1x,i5,3(1pg13.5))')
     &     n,ein(n),eout(n),eout(n)-ein(n)
        enddo
      endif
c
      return
      end



C#########################
C     CR 03-4-1998
C
C     rbrent calculates efermi in pairdn for finite temperature using
C     zerofkt
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Van Wijngaarden-Dekker-Brent Method for Root Finding
C
C W.H. Press et. al., Numerical Recipes in C, Sec. Ed, Cambridge 1995
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      FUNCTION RBRENT(E, DELTA, PH, GW, PAIRGW, RTMP, PARTNM, 
     &     NDIM, NMAX, IAB)
C
      include 'params.inc'
C
C----------------------------------------------------------------------------
C 
C      rtmp     = real temperature in MeV
C      e        = array(1..nmax) of single particle energies.
C      delta    = array(1..nmax) of gap parameters.
C      ph       = array(1..nmax) of the multiplicities of the single
C                 particle levels.
C                 This multiplicities excludes (!) the spin degeneracy,
C                 e.g. a 1s1/2 state has ph = 1 .
C      pairgw   = space for dynamical cut-off factor
C      partnm   = desired particle number
C
C      iab      = print level; significant values are -2, -1, 0 and 1;
C                 larger value prints more.
C
C      rxa,rxb,rxc,rlna,rlnb,rda,rdb are constants for finding roots
C                  of the zerofkt 
C
C
C      DATA FOR Wijngaarden-Dekker-Brent Method for Root
C
      DIMENSION E(NDIM),PH(NDIM),GW(NDIM),DELTA(NDIM),PAIRGW(NDIM)
C
      DATA     RA0,   RB0,   REPS,IITMAX,   RTOL
     &     /-100.0, 100.0,1.0e-14,   100,1.0e-16/
C
C-----------------------------------------------------------------------------
C
C     SETTING UP
C
      RA  = RA0
      RB  = RB0
      RFA = partnm 
     &     - parnmt(E, DELTA, PH, GW, PAIRGW, RTMP, RA, NDIM, NMAX)
      RFB = partnm 
     &     - parnmt(E, DELTA, PH, GW, PAIRGW, RTMP, RB, NDIM, NMAX)
      RFC = RFB
      RC  = RB
      RD  = RB-RA
      RE  = RD
C     
C     TESTING
C
      IF ((RFA .GT. 0.0).AND.(RFB .GT. 0.0)) THEN
         IF (RFA.GT.RFB) THEN
            RB = 1.0
         ELSE
            RB = 0.0
         ENDIF
      ELSE IF ((RFA .LT. 0.0).AND.(RFB .LT. 0.0)) THEN
         IF (RFA.GT.RFB) THEN
            RB = 0.0
         ELSE
            RB = 1.0
         ENDIF
      ELSE
C
C     ITERATION
C
         DO II=1,IITMAX
C     RENANME RA,RB,RC AND AJUST BOUNDING INTERVAL RD
            IF (((RFB .GT. 0.0).AND.(RFC .GT. 0.0))
     &           .OR.((RFB .LT. 0.0).AND.(RFC .LT. 0.0))) THEN
               RC  = RA
               RFC = RFA
               RD  = RB-RA
               RE  = RD
            ENDIF
            IF (ABS(RFC).LT.ABS(RFB)) THEN
               RA  = RB
               RB  = RC
               RC  = RA
               RFA = RFB
               RFB = RFC
               RFC = RFA
            ENDIF
C     CONVERGENCE CHECK
            RTOL1 = 2.0*REPS*ABS(RB) + 0.5*RTOL
            RXM   = 0.5*(RC-RB)
            IF ((ABS(RXM).LE.RTOL1).OR.((ABS(RFB) .EQ. 0.0))) GOTO 2999
C            IF ((ABS(RXM).LE.RTOL1)
C     &           .OR.((ABS(RXM).LE.(10*RTOL1))
C     &           .AND.(ABS(RFB) .LE. RTOL2)))
C     &           GOTO 2999
            IF ((ABS(RE).GE.RTOL1).OR.(ABS(RFA).GT.ABS(RFB))) THEN
C     ATTEMPT INVERSE QUADRATIC INTERPOLATION
               RS=RFB / RFA
               IF (RA .EQ. RC) THEN
                  RP = 2.0*RXM*RS
                  RQ = 1.0-RS
               ELSE
                  RQ = RFA / RFC
                  RR = RFB / RFC
                  RP = RS * (2.0*RXM*RQ*(RQ-RR)-(RB-RA)*(RR-1.0))
                  RQ = (RQ-1.0)*(RR-1.0) *(RS-1.0)
               ENDIF
C     CHECK WHETHER IN BOUNDS
               IF (RP .GT. 0.0) THEN
                  RQ = -RQ
               ENDIF
               RP    = ABS(RP)
               RMIN1 = 3.0*RXM*RQ-ABS(RTOL1*RQ)
               RMIN2 = ABS(RQ*RE)
               IF (2.0*RP .LT. MIN(RMIN1,RMIN2)) THEN
C     ACCEPT INTERPOLATION
                  RE = RD
                  RD = RP / RQ
C
                  IF(IAB .GT. 0) THEN
                     WRITE(6,'(A)') ' QUADRATIC INTERPOLATION '
                  ENDIF
C
               ELSE
C     INTERPOLATION FAILED -> BISECTION
                  RD = RXM
                  RE = RD
C
                  IF(IAB .GT. 0) THEN
                     WRITE(6,'(A)') ' FAILED -> BISECTION '
                  ENDIF
C
               ENDIF     
C     BOUNDS DECREASING TOO SLOWLY -> BISECTION
            ELSE
               RD = RXM
               RE = RD
C     
               IF(IAB .GT. 0) THEN
                  WRITE(6,'(A)') ' TOO SLOWLY -> BISECTION '
               ENDIF
C
            ENDIF
C     MOVE LAST BEST GUESS TO RA
            RA  = RB
            RFA = RFB
C     EVALUATE NEW TRIAL ROOT
            IF (ABS(RD) .GT. RTOL1) THEN
               RB = RB + RD
            ELSE
               RB = RB + SIGN(RTOL1,RXM)
            ENDIF
            RFBO = RFB
            RFB  = partnm 
     &            - parnmt( E, DELTA, PH, GW, PAIRGW, RTMP, RB, 
     &                      NDIM,NMAX)
C
            IF(IAB .GT. 0) THEN
               WRITE (6,'(I3,A,F28.18,A,F28.18,A,F28.18)')
     &              II,' :  rb = ',RB,'  rfb = ',RFB,'  drfb = ',
     &              (RFBO-RFB) 
            ENDIF
C
         ENDDO
         STOP ':( LOST'
 2999 ENDIF
C
      IF(IAB .GT. 0) THEN
         WRITE (6,'(A,F28.18)') 'EFERMI=NST= ', RB
      ENDIF
C
      RBRENT = RB
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CR 03-4-1998
C
C     'parnmt' computes the particle number for given 'efermi'
C              in case of thermal BCS at temperature 'RTMP'
C
c
      FUNCTION parnmt(E, DELTA, PH, GW, PAIRGW, RTMP, EFERMI, 
     &                NDIM, NMAX)
C
      include 'params.inc'
C
C----------------------------------------------------------------------------
C 
C      rtmp     = real temperature in MeV
C      e        = array(1..nmax) of single particle energies.
C      delta    = array(1..nmax) of gap parameters.
C      ph       = array(1..nmax) of the multiplicities of the single
C                 particle levels.
C                 This multiplicities excludes (!) the spin degeneracy,
C                 e.g. a 1s1/2 state has ph = 1 .
C      pairgw   = space for dynamical cut-off factor
C      partnm   = desired particle number
C      efermi   = fermi-energy
C
C      iab      = print level; significant values are -2, -1, 0 and 1;
C                 larger value prints more.
C
      DIMENSION E(NDIM),PH(NDIM),GW(NDIM),DELTA(NDIM),PAIRGW(NDIM)
C
C-----------------------------------------------------------------------------
C
      PARTN = 0.0
C
      DO NA = 1, NMAX
C
         DELTA2 = DELTA(NA)*PAIRGW(NA)
         DELTA2 = DELTA2*DELTA2
         EDIF   = E(NA)-EFERMI
         EQUASI = SQRT(EDIF*EDIF+DELTA2)
         TEXP   = 1.0 / (EXP(EQUASI/RTMP) + 1.0)
         GW(NA) = TEXP + (0.5 - TEXP) * (1.0 - EDIF/EQUASI)
C
         PARTN  = PH(NA)*GW(NA) + PARTN 
C
      ENDDO
C
      parnmt  = PARTN+PARTN
C
      RETURN
      END
C#########################
C     CR 03-4-1998
C
C     rbrent calculates efermi in pairdn for finite temperature using
C     zerofkt
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Van Wijngaarden-Dekker-Brent Method for Root Finding
C
C W.H. Press et. al., Numerical Recipes in C, Sec. Ed, Cambridge 1995
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      FUNCTION RBRENT0(E, DELTA, PH, GW, PAIRGW, PARTNM, 
     &     NDIM, NMAX, IAB,cutfac)
C
      include 'params.inc'
C
C----------------------------------------------------------------------------
C 
C      e        = array(1..nmax) of single particle energies.
C      delta    = array(1..nmax) of gap parameters.
C      ph       = array(1..nmax) of the multiplicities of the single
C                 particle levels.
C                 This multiplicities excludes (!) the spin degeneracy,
C                 e.g. a 1s1/2 state has ph = 1 .
C      pairgw   = space for dynamical cut-off factor
C      partnm   = desired particle number
C
C      iab      = print level; significant values are -2, -1, 0 and 1;
C                 larger value prints more.
C
C      rxa,rxb,rxc,rlna,rlnb,rda,rdb are constants for finding roots
C                  of the zerofkt 
C
C
C      DATA FOR Wijngaarden-Dekker-Brent Method for Root
C
      DIMENSION E(NDIM),PH(NDIM),GW(NDIM),DELTA(NDIM),PAIRGW(NDIM)
C
      DATA     RA0,   RB0,   REPS,IITMAX,   RTOL
     &     /-100.0, 100.0,1.0e-14,   100,1.0e-16/
C
C-----------------------------------------------------------------------------
C
C     SETTING UP
C
      RA  = RA0
      RB  = RB0
      RFA = partnm 
     &     - parnmt0(E, DELTA, PH, GW, PAIRGW, RA, NDIM, NMAX,cutfac)
      RFB = partnm 
     &     - parnmt0(E, DELTA, PH, GW, PAIRGW, RB, NDIM, NMAX,cutfac)
      RFC = RFB
      RC  = RB
      RD  = RB-RA
      RE  = RD
C     
C     TESTING
C
      IF ((RFA .GT. 0.0).AND.(RFB .GT. 0.0)) THEN
         IF (RFA.GT.RFB) THEN
            RB = 1.0
         ELSE
            RB = 0.0
         ENDIF
      ELSE IF ((RFA .LT. 0.0).AND.(RFB .LT. 0.0)) THEN
         IF (RFA.GT.RFB) THEN
            RB = 0.0
         ELSE
            RB = 1.0
         ENDIF
      ELSE
C
C     ITERATION
C
         DO II=1,IITMAX
C     RENANME RA,RB,RC AND AJUST BOUNDING INTERVAL RD
            IF (((RFB .GT. 0.0).AND.(RFC .GT. 0.0))
     &           .OR.((RFB .LT. 0.0).AND.(RFC .LT. 0.0))) THEN
               RC  = RA
               RFC = RFA
               RD  = RB-RA
               RE  = RD
            ENDIF
            IF (ABS(RFC).LT.ABS(RFB)) THEN
               RA  = RB
               RB  = RC
               RC  = RA
               RFA = RFB
               RFB = RFC
               RFC = RFA
            ENDIF
C     CONVERGENCE CHECK
            RTOL1 = 2.0*REPS*ABS(RB) + 0.5*RTOL
            RXM   = 0.5*(RC-RB)
            IF ((ABS(RXM).LE.RTOL1).OR.((ABS(RFB) .EQ. 0.0))) GOTO 2999
C            IF ((ABS(RXM).LE.RTOL1)
C     &           .OR.((ABS(RXM).LE.(10*RTOL1))
C     &           .AND.(ABS(RFB) .LE. RTOL2)))
C     &           GOTO 2999
            IF ((ABS(RE).GE.RTOL1).OR.(ABS(RFA).GT.ABS(RFB))) THEN
C     ATTEMPT INVERSE QUADRATIC INTERPOLATION
               RS=RFB / RFA
               IF (RA .EQ. RC) THEN
                  RP = 2.0*RXM*RS
                  RQ = 1.0-RS
               ELSE
                  RQ = RFA / RFC
                  RR = RFB / RFC
                  RP = RS * (2.0*RXM*RQ*(RQ-RR)-(RB-RA)*(RR-1.0))
                  RQ = (RQ-1.0)*(RR-1.0) *(RS-1.0)
               ENDIF
C     CHECK WHETHER IN BOUNDS
               IF (RP .GT. 0.0) THEN
                  RQ = -RQ
               ENDIF
               RP    = ABS(RP)
               RMIN1 = 3.0*RXM*RQ-ABS(RTOL1*RQ)
               RMIN2 = ABS(RQ*RE)
               IF (2.0*RP .LT. MIN(RMIN1,RMIN2)) THEN
C     ACCEPT INTERPOLATION
                  RE = RD
                  RD = RP / RQ
C
                  IF(IAB .GT. 0) THEN
                     WRITE(6,'(A)') ' QUADRATIC INTERPOLATION '
                  ENDIF
C
               ELSE
C     INTERPOLATION FAILED -> BISECTION
                  RD = RXM
                  RE = RD
C
                  IF(IAB .GT. 0) THEN
                     WRITE(6,'(A)') ' FAILED -> BISECTION '
                  ENDIF
C
               ENDIF     
C     BOUNDS DECREASING TOO SLOWLY -> BISECTION
            ELSE
               RD = RXM
               RE = RD
C     
               IF(IAB .GT. 0) THEN
                  WRITE(6,'(A)') ' TOO SLOWLY -> BISECTION '
               ENDIF
C
            ENDIF
C     MOVE LAST BEST GUESS TO RA
            RA  = RB
            RFA = RFB
C     EVALUATE NEW TRIAL ROOT
            IF (ABS(RD) .GT. RTOL1) THEN
               RB = RB + RD
            ELSE
               RB = RB + SIGN(RTOL1,RXM)
            ENDIF
            RFBO = RFB
            RFB  = partnm 
     &            - parnmt0( E, DELTA, PH, GW, PAIRGW, RB, 
     &                      NDIM,NMAX,cutfac)
C
            IF(IAB .GT. 0) THEN
               WRITE (6,'(I3,A,F28.18,A,F28.18,A,F28.18)')
     &              II,' :  rb = ',RB,'  rfb = ',RFB,'  drfb = ',
     &              (RFBO-RFB) 
            ENDIF
C
         ENDDO
         STOP ':( LOST'
 2999 ENDIF
C
      IF(IAB .GT. 0) THEN
         WRITE (6,'(A,F28.18)') 'EFERMI=NST= ', RB
      ENDIF
C
      RBRENT0 = RB
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

**************************************************************************
C     PGR 12.01.09
C
C     'parnmt0' computes the particle number for given 'efermi'
C               for BCS in case temperature zero
C
c
      FUNCTION parnmt0(E,DELTA,PH,GW,PAIRGW,EFERMI,NDIM,NMAX,cutfac)
C
      include 'params.inc'
C
C----------------------------------------------------------------------------
C 
C      e        = array(1..nmax) of single particle energies.
C      delta    = array(1..nmax) of gap parameters.
C      ph       = array(1..nmax) of the multiplicities of the single
C                 particle levels.
C                 This multiplicities excludes (!) the spin degeneracy,
C                 e.g. a 1s1/2 state has ph = 1 .
C      pairgw   = space for dynamical cut-off factor   ! old cutoff
C      efermi   = fermi-energy
c      cutfac   = cutoff band
C
C
      common /ctest/edirec(ksl2),efkincut,ifkincut
      DIMENSION E(NDIM),PH(NDIM),GW(NDIM),DELTA(NDIM),PAIRGW(NDIM)
C
C-----------------------------------------------------------------------------
C
!      IF(tphpair.AND.cutfac<0D0) pairgw(1:nmax) =
!      WRITE(*,*) 'eferm:',nmax,efermi
      IF(cutfac<0D0) THEN
        IF(ifkincut==1) THEN
          pairgw(1:nmax) =
     &     1D0/(1D0+EXP((edirec(1:nmax)-efkincut)/(-0.1D0*cutfac)))
        ELSE
          pairgw(1:nmax) =
     &     1D0/(1D0+EXP((e(1:nmax)-efermi+cutfac)/(-0.1D0*cutfac)))
          IF(tlowercut) pairgw(1:nmax) = pairgw(1:nmax)/
     &      (1D0+EXP(-(e(1:nmax)-efermi-cutfac)/(-0.1D0*cutfac)))
        END IF
      END IF

      PARTN = 0.0
C
      DO NA = 1, NMAX
C
         IF(tphpair) THEN
           DELTA2 = DELTA(NA)*PAIRGW(NA)  !old cutoff
         ELSE
           DELTA2 = DELTA(NA)  ! new cutoff
         END IF
         DELTA2 = DELTA2*DELTA2
         EDIF   = E(NA)-EFERMI
         EQUASI = SQRT(MAX(EDIF*EDIF+DELTA2,1D-200))
         GW(NA) = (0.5D0 - 0.5D0*EDIF/EQUASI) 
C
         IF(tphpair) THEN
           PARTN  = PH(NA)*GW(NA) + PARTN 
         ELSE
           PARTN  = PH(NA)*GW(NA)*pairgw(na) + PARTN   ! new cutoff
         END IF
C
      ENDDO
C
      parnmt0  = PARTN+PARTN
C
      RETURN
      END

**************************************************************************
C     PGR 12.01.09
C
C     'bisection' determines the Fermi energy brute force 
C                 by the bisection method
C
c
      subroutine bisection(E,DELTA,PH,GW,PAIRGW,EFERMI,NDIM,NMAX,
     &                     partnm,parnmr,ibisect,cutfac,efermic)
C
      include 'params.inc'
C
C----------------------------------------------------------------------------
C 
C      e        = array(1..nmax) of single particle energies.
C      delta    = array(1..nmax) of gap parameters.
C      ph       = array(1..nmax) of the multiplicities of the single
C                 particle levels.
C                 This multiplicities excludes (!) the spin degeneracy,
C                 e.g. a 1s1/2 state has ph = 1 .
C      pairgw   = space for dynamical cut-off factor
C      partnm   = desired particle number
C      efermi   = fermi-energy (result)
C      parnmr   = achieved poarticle number (result)
c      ibisect  = number of bisection steps
c      cutfac   = cutoff band (used in case of cutfac<0)
c      efermic  = the other side of the final interval 
c                 input: 0D0 starts from scratch
c                        .NEQ.0D0 starts from band
C
C
      DIMENSION E(NDIM),PH(NDIM),GW(NDIM),DELTA(NDIM),PAIRGW(NDIM)
C
      INTEGER, INTENT(IN) ::  IBISECT
      REAL(8),PARAMETER :: precnum=1D-12  ! precision in particle number
      REAL(8),PARAMETER :: precnumfac=1D0+precnum  ! precision factor
      REAL(8),PARAMETER :: precferm=3D-14  ! precision in Fermi energy
!      parameter (IBISECT=90)
C
C----------------------------------/-------------------------------------------
C
ctest               write(6,'(a,20(/1x,5g12.4))')
ctest     &             ' s.p. energy v**2:',(e(na),na=1,nmax)
ctest               write(6,'(a,20(/1x,5g12.4))')
ctest     &             ' occupation v**2:',(gw(na),na=1,nmax)
ctest      write(6,*) ' partnm=',partnm
C
C     determine upper and lower end of energy interval
C
       IF(efermic.NE.0D0) THEN
         efmax = MAX(efermi,efermic)
         efmin = MIN(efermi,efermic)
         eferm1 = efmin-1D-3
         eferm2 = efmax+1D-3
       ELSE
         eferm1 = MINVAL(e(1:NMAX))-1D-5
         eferm2 = MAXVAL(e(1:NMAX))+1D-5
       END IF
!       WRITE(*,*) 'BIS.:',efermi,efermic,eferm1,eferm2
!       WRITE(*,'(a,200(1pg13.5))') 'spe:',E(1:NMAX)
C
C      iterate
C
       parnm1 = parnmt0(E,DELTA,PH,GW,PAIRGW,eferm1,NDIM,NMAX,cutfac)
!       write(6,*) ' eferm1,parnm1=',eferm1,parnm1,eferm2
       if(parnm1.gt.partnm*precnumfac) THEN
         write(6,*) ' eferm1,eferm2,parnm1,parnm2=',
     &                eferm1,eferm2,parnm1,parnm2
         stop " BISECTION: lower limit too high"
       endif
       parnm2 = parnmt0(E,DELTA,PH,GW,PAIRGW,eferm2,NDIM,NMAX,cutfac)
       if(parnm2.lt.partnm/precnumfac) then
         write(6,*) ' eferm1,eferm2,parnm1,parnm2=',
     &                eferm1,eferm2,parnm1,parnm2
         stop " BISECTION: upper limit too low"
       endif
       do iter=1,IBISECT
         eferm3 = 0.5D0*(eferm1+eferm2)
         parnm3 = parnmt0(E,DELTA,PH,GW,PAIRGW,eferm3,NDIM,NMAX,cutfac)
!         write(6,*)' BIS: iter,eferm3,error=',iter,eferm3,parnm3-partnm
         if(parnm3.gt.partnm) then
           eferm2 = eferm3
           parnm2 = parnm3
         else
           eferm1 = eferm3
           parnm1 = parnm3
         endif
         IF(ABS(parnm1-partnm)<precnum*partnm) THEN
           efermi=eferm1
           efermic=eferm2
           parnmr =parnmt0(E,DELTA,PH,GW,PAIRGW,efermi,NDIM,NMAX,cutfac)
           RETURN
         ELSE IF(ABS(parnm2-partnm)<precnum*partnm) THEN
           efermi=eferm2
           efermic=eferm1
           parnmr=parnm2
           RETURN
         END IF
         IF(ABS(eferm1-eferm2)<precferm*eferm1) CYCLE
       enddo
C
C      final choice of optimum 
C
       if(abs(parnm1-partnm).lt.abs(parnm2-partnm)) then
         efermi = eferm1
         parnmr = parnm1
         efermic = eferm2
       else
         efermi = eferm2
         parnmr = parnm2
         efermic = eferm1
       endif
C
      RETURN
      END

**************************************************************************
C     PGR 03.06.19
C
C     'mixed' determines the Fermi energy and occupation numbers
C             by a mix of bisection and tangential step
C
C     !!  testing stage, 10% gain, probably not worth the risk
C
c
      subroutine mixed(E,DELTA,PH,GW,PAIRGW,EFERMI,NDIM,NMAX,
     &                     partnm,parnmr,ibisect,cutfac,efermic)
C
      include 'params.inc'
C
C----------------------------------------------------------------------------
C 
C      e        = array(1..nmax) of single particle energies.
C      delta    = array(1..nmax) of gap parameters.
C      ph       = array(1..nmax) of the multiplicities of the single
C                 particle levels.
C                 This multiplicities excludes (!) the spin degeneracy,
C                 e.g. a 1s1/2 state has ph = 1 .
C      gw       = occupation numbers v**2 (result)
C      pairgw   = space for dynamical cut-off factor
C      partnm   = desired particle number
C      efermi   = fermi-energy (result)
C      parnmr   = achieved poarticle number (result)
c      isecant  = number of bisection steps
c      cutfac   = cutoff band (used in case of cutfac<0)
c      efermic  = the other side of the final interval 
c                 input: 0D0 starts from scratch
c                        .NEQ.0D0 starts from band
C
C
      DIMENSION E(NDIM),PH(NDIM),GW(NDIM),DELTA(NDIM),PAIRGW(NDIM)
C
      INTEGER, INTENT(IN) ::  IBISECT
      REAL(8),PARAMETER :: precnum=1D-12  ! precision in particle number
      REAL(8),PARAMETER :: precnumfac=1D0+precnum  ! precision factor
      REAL(8),PARAMETER :: precferm=5D-14  ! precision in Fermi energy
      REAL(8),PARAMETER :: deferm=0.5D0   ! rough slope to cope with steep steps
!      parameter (IBISECT=90)
C
C----------------------------------/-------------------------------------------
C
ctest               write(6,'(a,20(/1x,5g12.4))')
ctest     &             ' s.p. energy v**2:',(e(na),na=1,nmax)
ctest               write(6,'(a,20(/1x,5g12.4))')
ctest     &             ' occupation v**2:',(gw(na),na=1,nmax)
ctest      write(6,*) ' partnm=',partnm
C
C     determine upper and lower end of energy interval
C
       IF(efermic.NE.0D0) THEN
         efmax = MAX(efermi,efermic)
         efmin = MIN(efermi,efermic)
         eferm1 = efmin-1D-3
         eferm2 = efmax+1D-3
       ELSE
         eferm1 = MINVAL(e(1:NMAX))-1D-5
         eferm2 = MAXVAL(e(1:NMAX))+1D-5
       END IF
!       WRITE(*,*) 'BIS.:',efermi,efermic,eferm1,eferm2
!       WRITE(*,'(a,200(1pg13.5))') 'spe:',E(1:NMAX)
C
C      iterate
C
       parnm1 = parnmt0(E,DELTA,PH,GW,PAIRGW,eferm1,NDIM,NMAX,cutfac)
     &         -partnm
!       write(6,*) ' eferm1,parnm1=',eferm1,parnm1,eferm2
       if(parnm1+partnm.gt.partnm*precnumfac) THEN
         write(6,*) ' eferm1,eferm2,parnm1,parnm2=',
     &                eferm1,eferm2,parnm1,parnm2
         stop " SECANT: lower limit too high"
       endif
       parnm2 = parnmt0(E,DELTA,PH,GW,PAIRGW,eferm2,NDIM,NMAX,cutfac)
     &         -partnm
       if(parnm2+partnm.lt.partnm/precnumfac) then
         write(6,*) ' eferm1,eferm2,parnm1,parnm2=',
     &                eferm1,eferm2,parnm1,parnm2
         stop " SECANT: upper limit too low"
       endif
       do iter=1,IBISECT

!       WRITE(6,'(a,i4,6(1pg13.5))')
!     &  'MIXED: iter,eferms,partnms=',
!     &    iter,eferm1,eferm2,parnm1,parnm2

!        check if sufficient precison has been achieved
         IF(ABS(parnm1)<precnum*partnm) THEN
           efermi=eferm1
           efermic=eferm2
           parnmr =parnmt0(E,DELTA,PH,GW,PAIRGW,efermi,NDIM,NMAX,cutfac)
           RETURN
         ELSE IF(ABS(parnm2)<precnum*partnm) THEN
           efermi=eferm2
           efermic=eferm1
           parnmr=parnm2
           RETURN
         END IF

!        check bracketing conditions
         IF(parnm1*parnm2>0) STOP 'SECANT: solution not bracketed'
         IF(parnm1>parnm2) STOP 'SECANT: wrong order of brackets'

         IF(parnm2<-parnm1) THEN
!          tangential step from above - reference point eferm2
           eferm3 = eferm2-deferm
           parnm3 = parnmt0(E,DELTA,PH,GW,PAIRGW,eferm3,NDIM,NMAX,
     &                      cutfac)  - partnm
           slope  = -(parnm3-parnm2)/deferm
           eferm3 = eferm2 - parnm2/slope
!           WRITE(6,'(a,5(1pg13.5))')'slope2:',eferm2,slope,parnm2,eferm3
!     &              ,(eferm1+eferm2)/2D0
           IF(eferm3.GE.eferm2) STOP 'SECANT2: slope in wrong direction'
           eferm3=MAX(eferm3,(eferm1+eferm2)/2D0)
           parnm3 = parnmt0(E,DELTA,PH,GW,PAIRGW,eferm3,NDIM,NMAX,
     &                      cutfac)  - partnm
!           WRITE(6,'(a,2(1pg13.5))') '       ',eferm3,parnm3
         ELSE
!          tangential step from below - reference point eferm1
           eferm3 = eferm1+deferm
           parnm3 = parnmt0(E,DELTA,PH,GW,PAIRGW,eferm3,NDIM,NMAX,
     &                      cutfac)  - partnm
           slope  = (parnm3-parnm1)/deferm
           eferm3 = eferm1 - parnm1/slope
!           WRITE(6,'(a,5(1pg13.5))')'slope1:',eferm1,slope,parnm1,eferm3
!     &           ,(eferm1+eferm2)/2D0
           IF(eferm3.LE.eferm1) STOP 'SECANT1: slope in wrong direction'
           eferm3=MIN(eferm3,(eferm1+eferm2)/2D0)
           parnm3 = parnmt0(E,DELTA,PH,GW,PAIRGW,eferm3,NDIM,NMAX,
     &                      cutfac)  - partnm
!           WRITE(6,'(a,2(1pg13.5))') '       ',eferm3,parnm3
         END IF

!        define new brackets    
         if(parnm3>0) then
           eferm2 = eferm3
           parnm2 = parnm3
         else
           eferm1 = eferm3
           parnm1 = parnm3
         endif

 !       Fermi energy converged (particle number perhaps not fully exact)
         IF(ABS(eferm1-eferm2)<precferm*eferm1) THEN
           IF(parnm2<-parnm1) THEN
             eferm  = eferm2
             parnmr = parnm2+partnm
             efermc = eferm1
           ELSE
             eferm  = eferm1
             parnmr = parnm1+partnm
             efermc = eferm2
           END IF
           WRITE(6,'(a,3(1pg13.5))')
     &          'SECANT Fermi energy converged: eferm,partnum,error=',
     &           eferm,parnmr,parnmr-partnm
          IF(eferm.NE.eferm3)
     &     parnm3 = parnmt0(E,DELTA,PH,GW,PAIRGW,eferm3,NDIM,NMAX,
     &                      cutfac) 
           RETURN
         END IF
       enddo
C
C      final choice of optimum 
C
       IF(parnm2<-parnm1) THEN
         eferm  = eferm2
         parnmr = parnm2+partnm
         efermc = eferm1
       ELSE
         eferm  = eferm1
         parnmr = parnm1+partnm
         efermc = eferm2
       END IF
       IF(eferm.NE.eferm3)
     &   parnm3 = parnmt0(E,DELTA,PH,GW,PAIRGW,eferm3,NDIM,NMAX,
     &                    cutfac) 
C
      RETURN
      END

c-----pair_stab------------ part of the hartree-fock package -----------------
c  
c     modify gap according to stabilized pairing
c
      subroutine pair_stab(e,delta,gw,ph,pairgw,v2occ,epair,
     &                     ecstab,nmax,ndim,iab)
c
c     input via list:
c      e       = array(1..nmax) of single particle energies.
c      delta   = array(1..nmax) of gap parameters.
c      gw      = array(1..nmax) of occupation probabilities,
c      ph      = array(1..nmax) of the multiplicities of the single
c                particle levels.
c                This multiplicities excludes (!) the spin degeneracy,
c                e.g. a 1s1/2 state has ph = 1 .
c      pairgw  = space for dynamical cut-off factor
c      ecstab  = stabilizing energy
c      nmax    = maximum number of single particle states
c                (needs to be smaller than the parameter 'ndim')
c      ndim    = dimensioning of the arrays, particularly of 'vpairf',
c                in the calling routine
c      iab     = print level; significant values are -2, -1, 0 and 1;
c                larger value prints more.
c     output via list:
c      delta   = array(1..nmax) of gap parameters.
c      epair   = pairing energy
c
      include 'params.inc'
      REAL(8) :: e(ndim),gw(ndim),ph(ndim),delta(ndim),pairgw(ndim)
      REAL(8) :: v2occ(ndim)
      LOGICAL,PARAMETER :: tprint=.FALSE.
c
!      IF(ecstab==0D0) RETURN
c
      IF(tprint) THEN
        WRITE(6,'(a,i4,200(1pg13.5))') 'PAIR_STAB entry: ',nmax
        WRITE(6,'(a,200(1pg13.5))') ' delta: ',delta(1:nmax)
        WRITE(6,'(a,200(1pg13.5))') ' v2occ: ',v2occ(1:nmax)
        WRITE(6,'(a,200(1pg13.5))') 'pairgw: ',pairgw(1:nmax)
        WRITE(6,'(a,200(1pg13.5))') '    ph: ',ph(1:nmax)
      END IF
c
c     compute pairing energy
c
      sumduv = 0D0
      DO  na=1,nmax
         IF(tphpair) THEN
           deltaa = delta(na)*pairgw(na)   !  old cutoff
         ELSE
           deltaa = delta(na)  ! new cutoff
         END IF
!         deltab = delta(na)
!         edif   = e(na)-eferm
!         equasi = sqrt(edif*edif+deltaa*deltaa)
!         v2     = (0.5D0-0.5D0*edif/equasi)
         v2     = v2occ(na)
         uv     = SQRT(MAX(v2-v2*v2,xsmall)) 
         IF(.NOT.tphpair) THEN
           v2     = v2*pairgw(na)   ! new cutoff
           uv     = uv*pairgw(na)   ! new cutoff
         END IF
         sumduv = ph(na)*uv*deltaa + sumduv
      ENDDO
C
      epair  = -sumduv
      if(iab.gt.1) write(6,'(a,3g12.4)') 'PAIR_STAB epair:',epair

      IF(ecstab==0D0) RETURN
c
c     modify gaps
c
      DO na=1,nmax
        delta(na) = delta(na)*(1D0+(ecstab/epair)**2)
        delta(na) = min(delta(na),5D0)
      ENDDO

c
      RETURN
      END


**************************************************************************



      include 'blocking.for'
      include 'LN.f'
