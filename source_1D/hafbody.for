c-----block ------------ part of the hartree-fock package ---------------------c
      blockdata  block
c
      include 'params.inc'
!      integer counter 
      parameter (kmaxp2=kmaxp*kmaxp)
      parameter (kmaxn2=kmaxn*kmaxn)
c
      include 'grid.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'wavfun.inc'
      include 'rmsnuc.inc'
      include 'constn.inc'
c      include 'linres.inc'
      REAL(8) :: pairenhance=2D0,v2limh=1D-6,delmix=0.8D0
      INTEGER :: iterhfb=400,iternat=0,modnat=999999
      INTEGER :: iterbcs=200
!      LOGICAL :: tregpair=.FALSE.  ! .TRUE. for regularized pairing
!      LOGICAL :: tphpair = .TRUE.  ! .TRUE. for uv-cutoff
!      LOGICAL :: tnamelist=.FALSE.
c
c-------------------------------------------------------------------------------
c
c     data in common /grid/
c
      data ngrid,rstep
     &    /    0,  0.3/
c
c     data in common /numpar/
c
      data itmax,itpr,  dabr,  dabe,addnew,x0damp,e0damp, epspr,wall
     &    /  400,  50,0.0001,0.0001,   0.5,   0.1, 110.0,1.0e-5,-1.0/
!     &    /  400,   5,0.0001,0.0001,   0.5,   0.2, 110.0,1.0e-5,-1.0/
      data iwr,ifplot,ifrhoc,ifanlr,ifgrho
     &    / -1,     0,     0,     0,     0/
      data ifex,iftm,ifenrg,ifrso,ifgsc,ifzpe,ifdfms,ifenav,ifsave
     &    /   1,   0,     1,    1,    0,    0,     0,     0,     1/
      data itrpa,epsrpa
     &    /    0,   0.5/
c
c     data in common /skypar/
c
      data  t5, t6, t7, t8, t9, x5, x6, x7, x8, b4p,amfc0
     &    /0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, 0D0,0D0,0D0/
      data b5p, b6,b6p, b7,b7p
     &    /0D0,0D0,0D0,0D0,0D0/
      data  ipair,ilnrsp,iparsh,pairfp, pairfn,acoul,beta20,temprt
     &    /     1,     1,     0,   0D0,    0D0,  1.0,   0D0,   0D0/
      data   c1,c1p, c4,c4p,tte,tto
     &     /0D0,0D0,0D0,0D0,0D0,0D0/
      data cutfac,atajima,ecstab,alpinv,betinv,alpin2,betin2,riftm1
     &    / 1.65 ,    0.5,   0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,   0D0/
!     &    / 1.65 ,    0.0,   0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,   0D0/
      data b3tilde,b12,b12p,b7,b7p,pow7
     &    /    0D0,0D0,0D0,0D0,0D0, 0D0/
      data rho0pr,pairgm,rhoswm,exmass,powerp,hfbwidth
     &    /  0.16,   1.0,   0.0,   0.0, 0.0D0, 0.0D0/
c     &    /  0.16,0.1666,   0.0,   0.0/
      data str_squark
     &    /  -0.24/
      data imfold,ioptsk
     &    /     0,     0/
c
c     data in common /wavfun/
c
      data bmass/kgr2*0.0/
      data pot/kgr2*0.0/
      data potls/kgr2*0.0/
      data ttthrm/ksl2*1.0/
      data gapgoal,xmugap,xlamgap,c0gap/0D0,1D0,0D0,0.5D0/ ! add to 'for005' later
c
c     data in common /shells/
c
      data nprot,nneut,nmaxp,nmaxn
     &    /    8,    8,    0,    0/
cold      data jsp/                 1,                         3,2,
cold     &                      5,4,1,                     7,6,3,2,
cold     &                  9,8,5,4,1,               11,10,7,6,3,2,
cold     &            13,12,9,8,5,4,1,         15,14,11,10,7,6,3,2,
cold     &      17,16,13,12,9,8,5,4,1,   19,18,15,14,11,10,7,6,3,2,
cold     & 21,20,17,16,13,12,9,8,5,4,1,
cold     &                         23,22,19,18,15,14,11,10,7,6,3,2,
cold     & 25,24,21,20,17,16,13,12,9,8,5,4,1,
cold     & 27/
cold      data wls/0.,-2.,1.,-3.,2.,-4.,3.,-5.,4.,-6.,5.,-7.,6.,-8.,7.,
cold     *         -9.,8.,-10.,9., -11., 10.,
cold     &         -12.,11.,-13.,12.,-14.,13.,
cold     &         -15.,14., -16.,15., -17.,16.,  -18.,17., -19.,18.,
cold     &         -20.,19., -21.,20., -22.,21.,  -23.,22., -24.,23./
cold      data zwj/.15915494,.15915494,.31830989,.31830989,.47746483,
cold     &         .47746483,.63661977,.63661977,.79577472,.79577472,
cold     &         .95492966,.95492966,1.1140846,1.1140846,1.2732395,
cold     &         1.2732395,1.4232945,1.4232945,1.5915494,1.5915494,
cold     &         1.7507044,1.7507044,1.9098593,1.9098593,2.0690143,
cold     &         2.0690143,2.2281692/
cold      data zwj/
cold     & 0.1591549, 0.1591549, 0.3183099, 0.3183099, 0.4774648, 0.4774648,
cold     & 0.6366197, 0.6366197, 0.7957747, 0.7957747, 0.9549296, 0.9549296,
cold     & 1.1140846, 1.1140846, 1.2732395, 1.2732395, 1.4323945, 1.4323945,
cold     & 1.5915494, 1.5915494, 1.7507043, 1.7507043, 1.9098593, 1.9098593,
cold     & 2.0690143, 2.0690143, 2.2281692, 2.2281692, 2.3873241, 2.3873241,
cold     & 2.5464790, 2.5464790, 2.7056339, 2.7056339, 2.8647890, 2.8647890,
cold     & 3.0239439, 3.0239439, 3.1830988, 3.1830988, 3.3422537, 3.3422537,
cold     & 3.5014086, 3.5014086, 3.6605637, 3.6605637, 3.8197186/
c      data pairmp/kmaxp2*0.0/
c      data pairmn/kmaxn2*0.0/
c
c      data ck/10*1.e-35/
c
c     data in common /rmsnuc/
c
      data ifform/1/
      data apro0,apro11,apro12,bpro11,bpro12
     *    / 0.74, 0.0  ,  0.0 ,  0.0 ,  0.0 /
      data aneu0,aneu11,aneu12,bneu11,bneu12
     *    / 0.74, 0.0  ,  0.0 ,  0.0 ,  0.0 /
c
c     data for common /constn/
c
cold      data pi,         sq2,        sqfrpi,     d4pi
cold     *    /3.141592654,1.414213562,3.544907702,.0795774155/
      data pi,         sqfrpi,     d4pi
     *    /3.14159265358979D0,3.54490770181103D0,7.957747154594767D-02/
      data       fourpi,  xsmall
     *    /12.5663706143592D0, 1.0D-40/
c
c     defaults for Fayans functional
c
      DATA sigmavol,rho0fy/0.3333333333D0,0.08D0/
      DATA hsp,hgradp,asm,hsm,hgradm/0D0,0D0,0D0,0D0,0D0/
      DATA hcoul,sigmacoul/0D0,0.33333333333D0/
      DATA h1xi,hgradxi,gammaxi,rho0xi/0D0,0D0,0.66666666667D0,0.16D0/
      DATA c0tau,c1tau/0D0,0D0/
      DATA c0nabj,c1nabj,c0tens,c1tens/0D0,0D0,0D0,00D0/
!      DATA tpairFy/.true./
      DATA tregpair/.FALSE./
      DATA tphpair/.FALSE./
      DATA isovectormodel/3/
      DATA tnamelist/.FALSE./
c
c     data for common /linres/
c
c      data ifsepi,resdmp
c     &    /     0,   0.0/
c
       end
c-----hafosu------------ part of the hartree-fock package ----------------------
c
      subroutine hafosu
c
      include 'params.inc'
c
c     * * * * * * *    skyrme-hartree-fock with pairing    * * * * * * *
c
c     solution of the spherical hartree-fock equations with schematic
c     pairing for the occupation probabilities in open shell nuclei.
c
c     input parameters via common /grid/:
c      rstep  = radial grid size (in fm).
c               recommended and preoccupied is 0.3fm.
c      ngrid  = number of radial grid points. if ngrid=0 is given
c               then the program determines a value which is
c               optimal in most cases.
c
c     input parameters via common /shells/:
c      nprot  = number of protons
c      nneut  = number of neutrons
c      nmaxp  = number of proton shells
c      nmaxn  = number of neutron shells
c               (note that the shells are counted in the ordering
c                of the harmonic oscillator shell model.)
c
c     input parameters via common /numpar/:
c      dabr   = "radius-condition" for terminating the iteration;
c               difference in r.m.s. radii should come below that
c               value. a reasonable value is preoccupied in blockdata.
c      dabe   = "energy-condition" for terminating the iteration;
c               difference in energy/particle should come below that
c               value. a reasonable value is preoccupied in blockdata.
c      addnew = weight with which the new h.f. potentials are mixed
c               to the old ones. this mixing is to ensure stable
c               iteration. a reasonable value is addnew=0.5 .
c      x0damp = switch for the step
c               and size of gradient step kinetic energy damping.
c               switch: 'x0damp'=0 --> damping by direct h-inversion
c                       'x0damp'>0 --> kinetic energy damping
c      e0damp = energy parameter for kinetic energy damping
c      ifex   = for 'ifex=1' an approximate coulomb-exchange force
c               is taken into account.
c      iftm   = for 'iftm=1' the contribution from the kinetic terms
c               (skyrme paramters t1 & t2) to the spin-orbit potential
c               is taken into account.
c               note: the option 'ifex' and 'iftm' are explained in
c                     phys.rev. c33 (1986) 335 .
c      itmax  = maximum number of iterations.
c      itpr   = maximum number of subiterations in the pairing routine.
c      iwr    = print level: -2 --> nothing is printed, except warnings
c                            -1 --> standard output.
c                             0 --> more info. during iteration
c                             1 --> wave functions are printed (!)
c      ifdfms = switches the nucleon masses:
c                0 --> equal masses 20.7525 are used
c                1 --> proton and neutron mass are different
c                2 --> equal masses 20.73553 are used
c                3 --> equal masses 20.73621 are used (Pearson, Farine)
c               the detailed values are given in subroutine 'initfc'.
c      ifenrg = skyrme potentials and energies are computed from
c               scratch (i.e. without mixing 'addnew') after the
c               iteration.
c      ifrso  = for 'ifrso=1' the contribution of the spin-orbit
c               current to the charge density is included (see
c               friar & negele, adv.nucl.phys.).
c      ifgsc  = switches handling of theories 'on top' of hartree-fock:
c               0 --> no further calculation;
c               1 --> calls rpa&gsc program 'sumrul';
c               2 --> calls continuum rpa program 'extrsp';
c               3 --> calls polarization correction program 'specor':
c               4 --> calls printout of separable matrix elements 'resmte'
c               7 --> explicit time evolution of linear perturbation
c      ifzpe  = switches type of c.m.correction:
c               -4 --> nothing for energy but correction for rho
c               -3 --> energy as <p**2>/(2*a*m) and nothing for rho
c               -2 --> no correction at all
c               -1 --> factor (1-1/a) on the mass is used for correction
c                      c.m. correction from actual ZPE
c                0 --> factor (1-1/a) on the mass is used for correction
c                      c.m. correction uses average guess
c                1 --> correction is computed explicitely and a
c                      posteriori as <p**2>/(2*a*m).
c                2 --> correction is computed explicitely and before
c                      variation as <p**2>/(2*a*m).
c      idisc  = for 'idisc=+1' wave functions and other results are
c               stored on disk and read in from disk in subsequent
c               calls; usage is recommended e.g. in a fitting routine
c               where 'hafosu' is called many times.
c               to disable put 'idisc=-1'.
c      ifrhoc = switch for charge- or proton-density in coulomb energy:
c                0 --> proton-density is used for coulomb potential
c                1 --> charge-density is used for coulomb potential
c      ifenav = a switch for energy-averages in width
c                0 --> no energy averages are calculated
c                1 --> energy averages are calculated
c      ifplot = switch for compressed information on file 11,
c                     e.g. for plotting
c               0  -->  nothing is written on file 11
c               1  -->  global results
c               2  -->  densities in addition
c               3  -->  and furthermore the single particle energies
c     ioptsk = switch for variants of the Skyrme functional
c               0 --> standard
c               1 --> b5 as input and density dependent mass (surface?)
c                     (only important for 'hffit.for')
c               2 --> rational approximant for isoscalar density dep. term
c               3 --> rational approximant for isoscalar and isovector term
c                     (isovector requires "exmass" = parameter(14) as input)
c               4 --> MF like kinetic&ls term (use 'b1kin' and 'b5kin')
c               5 --> rho-dependent kinetic and l*s, direct version
c               6 --> rho-dependent kinetic and l*s, inverse version
c               7 --> BHF functional for rho-dependent terms, needs B1,B1P=0
c               8 --> alternative isovector b3 and dens.dep. gradient
c               9 --> extra density dependent term with own 'power'
c              10 --> flexible tensor terms (parameter 'riftm')
c              11 --> accumulates the options 5, 8, 9, and 10
c              12 --> QMC model
c              13 --> Fayans functional
c
c     input parameters via common /skypar/:
c      t0-t4,x0-x4 = skyrme force parameters for the standard model
c                    (t4=0 signals RMF lke kinetic&ls functional)
c      t5-t8,x5-x8 = parameters for the extended model
c      tyuk,xyuk,ayuk = parameters for the yukawa-folding force
c                       the parameters are scaled to correspond
c                       to t0 and x0 in the zero-range limit.
c                       the ayuk is the range of the yukawa potential
c                       in units of [fm].
c                 
c                       note: the zero-range equivalent strength of
c                             the yukawa potential is subtracted from
c                             the zer0-range two-body force. thus the
c                             parameter t0 still describes the overall
c                             strength.
c                       note: ayuk serves as switch to folded SHF
c                             ayuk.ne.0.0 and tyuk.eq.0.0
c                             signal the folded SHF
c       ipair          = type of pairing:
C                        -1 = filling approximation
C                         0 = constant Delta approximation
C                         1 = constant Force approximation (BCS)
C                         2 = constant Force approximation with
c                             s.p. correction (BCS)
C                         3= Lipkin-Nogami
C                         4 = thermal HF, temperature handled by "DELTP"
C                         5 = BCS with full pairing matrix from delta force
C                         6 = HFB with full pairing matrix from delta force
C                         7 = BCS with surface pairing force
C                                V_0 delta(1-2) (1-rho(r)/rho0pr0)
C                         8 = HFB with surface pairing force
c                         9 = Lipkin-Nogami BCS with VDI 
c                        10 = Lipkin-Nogami BCS with DDDI 
c       iparsh         = number of extra shells above automatic setting
c      pairfp,pairfn   = for proton, or neutron, pairing;
c                       if 'pairfp/n'=0 then a default strenght
c                                 is taken,
c                              ipair=0 ->  delta=11.2mev/A**.5
c                              ipair=1,2,3 -> from Moeller/Nix in subr. GAVRG
c                              ipair=4 ->  forbidden
c                              ipair=5,6 ->  standard V_0 = 260.0/300.0
c                              ipair=7,8 ->  standard V_0 = 1000.0, rho0pr=0.16
c                       else 'pairfp/n' is the input for the strength, i.e.
c                              ipair=0 ->  Deltap/n
c                              ipair=1,2,3 -> pairing Force G*A
c                              ipair=4 ->  temperature T
c                              ipair=5,6 ->  V_0 in delta-pairing-force
c                              ipair=8,7 ->  V_0 in surface-pairing-force
c      acoul  = factor on the coulomb force. standard is 'acoul=1.0'.
c               acoul=0.0 disables coulomb force.
c
c     input parameters via common /rmsnuc/:
c      ifform = inclusion of folding with intrinsic nucleon formfactor.
c                1  ==> free nucleon formfactor according to mainz data;
c                0  ==> no intrinsic formfactor;
c              the following option disabled:
c                -1 ==> density-dependent nucleon formfactor parmetrised
c                       by the parameters
c                       'apro0,apro11,apro12,bpro11,bpro12,
c                        aneu0,aneu11,aneu12,bneu11,bneu12'
c                       in common /rmsnuc/.
c
c
c     output parameters via common /result/:
c      rms..  = r.m.s.radii; ending means:
c            .. = nothing --> total radius
c            .. = c       --> charge radius
c            .. = c0      --> charge radius before gr.st. correlations
c            .. = p       --> proton radius
c            .. = n       --> neutron radius
c      rdms.. = diffraction radii; ending as above.
c      surf.. = surface thickness; ending as above.
c      etotal = total binding energy
c      etot0  = total binding energy before correlations
c      eshl(iso) = proton/neutron-shell part of the energies
c      enrear = rearrangment energy
c      epair  = pairing contribution to the energy
c      zpe    = c.m.correction to the energy (for 'ifzpe=1')
c      spenrg = single particle energies
c      rhochr = charge density, including all nucleon foldings
c      formf  = charge formfactor, including all foldings
c
c     resulting fields via common /wavfun/
c      wf     = single particle wave functions
c               (the modus for storage of wave functions is explained in
c                subroutine 'initgr')
c      bmass  = inverse effective mass
c      pot    = potential
c      potls  = spin-orbit potential
c      rho    = density
c               (all densities and potentials are stored in a linear
c                array. protons occupy the first interval [1..'ngrid']
c                and neutrons the second interval ['ngrid'+1...2*'ngrid']).
c
c     dimensions of bookkeeping fields (e.g. 'jsp') are
c     are specified such that nuclei up to 208-pb can be calculated.
c
c     character*1 tdum
      integer counter 
      logical*4   top11,top12,top13,top14
      logical     tnopot       ! switch to case without Skyrme potential
      LOGICAL, PARAMETER :: tcentrdep=.FALSE.
      LOGICAL, PARAMETER :: thfbnatprint=.TRUE.,thfbnatural=.TRUE.
      LOGICAL :: tprint96
!      LOGICAL :: thfbfac=.FALSE.
      REAL(8) :: pairaddfacp=0D0,pairaddfacn=0D0
      REAL(8) :: hfbfacact=1D0
      REAL(8) :: symcondsav(2)
cold      dimension ph(277)
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'rmsnuc.inc'
      include 'linres.inc'
      include 'constn.inc'
c
cold      data ph/1.,
cold     &        2.,1.,
cold     &        3.,2.,1.,
cold     &        4.,3.,2.,1.,
cold     &        5.,4.,3.,2.,1.,
cold     &        6.,5.,4.,3.,2.,1.,
cold     &        7.,6.,5.,4.,3.,2.,1.,
cold     &        8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        9.,8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        10.,9.,8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        14.,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        15.,14.,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        16.,15.,14.,13.,12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        17.,16.,15.,14.,13.,12.,11.,10.,9.,8.,7.,6.,
cold     &            5.,4.,3.,2.,1.,
cold     &        18.,17.,16.,15.,14.,13.,12.,11.,10.,9.,8.,7.,6.,
cold     &            5.,4.,3.,2.,1.,
cold     &        19.,18.,17.,16.,15.,14.,13.,12.,11.,10.,9.,8.,7.,6.,
cold     &            5.,4.,3.,2.,1.,
cold     &        20.,19.,18.,17.,16.,15.,14.,13.,12.,11.,10.,9.,
cold     &            8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        21.,20.,19.,18.,17.,16.,15.,14.,13.,12.,11.,10.,9.,
cold     &            8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        22.,21.,20.,19.,18.,17.,16.,15.,14.,13.,12.,11.,10.,9.,
cold     &            8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        23.,22.,21.,20.,19.,18.,17.,16.,15.,14.,13.,12.,11.,10.,
cold     &            9.,8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        24./
c
c      data ifop11/0/
      data ifpoic/1/                 ! switch to point Coulomb
c
c-----------------------------------------------------------------------
c
c      write(2,*)
c     &  ' tyuk,xyuk,ayuk=',tyuk,xyuk,ayuk
c     optionally read previous solutions for nucleus with
c     given 'nprot' and 'nneut' from disk.
c
!      ecstab=0.001
!     set auxiliary field to default (=inactive)
!      WRITE(6,*) 'iternat=',iternat
      regpairadd=1D0
!     optionally read parameters for pairing-energy constraint
!     --> attention: not finally tested for Fayans pairing
      INQUIRE(FILE='pairaddfac',EXIST=tprint96)
      IF(tprint96) THEN
        OPEN(111,file='pairaddfac')
          READ(111,*) pairaddfacp,pairaddfacn
        CLOSE(111)
      END IF
! optionally read parameters for HFB iteration
      INQUIRE(FILE='hfbparams',EXIST=tprint96)
      tprint96=thfbnatprint .AND. thfbnatural .AND. tprint96
      IF(tprint96 .AND. .NOT.tnamelist) THEN
        !  preliminary test switch, implement folding in HFB sub.
!        IF(thfbfold.AND.iternat>0) 
!     &      STOP 'folding and HFB subiterations yet incompatible'
        OPEN(111,file='hfbparams')
        READ(111,*) iterhfb,iternat,modnat,pairenhance,iterbcs,thfbfac
        WRITE(96,'(1x/1x)')
        WRITE(96,'(a,3i5,f7.2,i5,l3)')
     &    '# iterhfb,iternat,modnat,pairenhance,iterbcs=',
     &    iterhfb,iternat,modnat,pairenhance,iterbcs,thfbfac
        CLOSE(111)
      END IF 
!      idisc=-1
      if(idisc.ge.0) then
        call buffio('i')
        deltps = deltp              !  save it for later use
        deltns = deltn
      endif
      if(TMONOPC) then
        constr0 = wall
        wall = -1D0
      else
        constr0 = 0D0
      endif
!     optionally vary pairing strenght
      IF(pairaddfacn.NE.0D0) THEN
        pairfn0=pairfn
        pairfn=pairfn*(1D0+pairaddfacn)
        fxin0=fxin
        h1xin0=h1xin
        hgradxin0=hgradxin
        fxin=fxin*(1D0+pairaddfacn)
        h1xin=h1xin*(1D0+pairaddfacn)
        hgradxin=hgradxin*(1D0+pairaddfacn)
      END IF 
      IF(pairaddfacp.NE.0D0) THEN
        pairfp0=pairfp
        pairfp=pairfp*(1D0+pairaddfacp)
        fxip0=fxip
        h1xip0=h1xip
        hgradxip0=hgradxip
        fxip=fxip*(1D0+pairaddfacp)
        h1xip=h1xip*(1D0+pairaddfacp)
        hgradxip=hgradxip*(1D0+pairaddfacp)
      END IF 
      IF(hfbwidth.NE.0D0) THEN
        thfbfold = .TRUE.
      ELSE
        thfbfold = .FALSE.
      END IF
c
c     check consistency of pairing options
c
      if(ipair.gt.10) stop ' option "ipair">10 not supplied '
      if((ipair.eq.6.or.ipair.eq.8) .and. x0damp.eq.zero)
     &   stop ' HFB requires iteration with finite X0DAMP'
      if(ipair.lt.-1) stop ' option "ipair"<-1 not supplied '
c      if(kmaxn.gt.277)
c     & stop ' "KMAXN" to large for provided multiplicities "PH"'
c      if(kmaxp.gt.277)
c     & stop ' "KMAXP" to large for provided multiplicities "PH"'
      if((khfbn.eq.1 .or. khfbp.eq.1) .and.
     &   (ipair.eq.6 .or. ipair.eq.8))
     &  stop ' to run HFB compile with KIFHFB=1 in PARAMS.INC'
      if(ioptsk.eq.3 .and. exmass.eq.zero)
     &  stop ' IOPTSK=3 requires nonzero EXMASS'
!      IF(ioptsk.EQ.13 .AND. ipair.NE.7)
!      IF(ioptsk.EQ.13 .AND. ipair.LT.7)
!     &  stop ' IOPTSK=13 requires IPAIR=7'
!      WRITE(6,*) 'iternat=',iternat
      IF(iternat.NE.0 .AND. (ipair.eq.8 .or. ipair.eq.6) .AND.
     &   ecstab<0D0)  STOP 'HFB sub-iterations require ECSTAB=0.0'
!     &   ecstab.NE.0D0)  STOP 'HFB sub-iterations require ECSTAB=0.0'
c
c     check wheter free nucleon gas is required --> 'tnopot'
c     this case should only run after a full Skyrme iteration
c     has been done before and delivered a density on 'rho'
c
      if(abs(t0)+abs(t1)+abs(t2)+abs(t3).le.xsmall
     &   .AND. ioptsk.NE.12  .AND. ioptsk.NE.13) then
        tnopot = .true.
      else
        tnopot = .false.
      endif
c
c     initialise default grid size and check bounds
c
      anumb  = nprot+nneut
      if(ngrid.eq.0) then
!        ngrid=int((1.2D0*anumb**.3333333D0+15D0)/rstep)+1  ! includes exotic n.
        ngrid=int((1.2*anumb**.3333333+6.4)/rstep)+1  ! old estimate
!        ngrid = 1.4*ngrid+1
        if(ifgsc.eq.1) ngrid = 1.15*ngrid
c        if(ifgsc.eq.1) ngrid = 2*ngrid
      endif
c
      if(ngrid.gt.kg1) stop 'too large "ngrid"'
      if(nmaxp.gt.kmaxp) stop 'too large "nmaxp"'
      if(nmaxn.gt.kmaxn) stop 'too large "nmaxn"'
c
c     initialise force, grid, wave functions and occupation weights
c
!      WRITE(6,'(6(1pg13.5))') b0,b1,b2,b3,b0p,b1p,b2p,b3p,power
      call initgr
      call initfc
      if(idisc.lt.0 .or. ifread.ne.1) then
       call initwf
ctest       write(*,*) ' INITWF done'
      else
!       call initgw
       deltp = deltps            ! restore read value
       deltn = deltns
      endif
      IF(thfbfold) CALL test_Gauss_folding(ngrid,rstep,hfbwidth)
!      WRITE(6,*) 'after initialization: tpair,ipair=',tpair,tpairFy,ipair
!      WRITE(6,'(a,2i3)') 'IFPAIR:',ifpair
c
c     test nuclear matter
c
!      WRITE(6,'(6(1pg13.5))') b0,b1,b2,b3,b0p,b1p,b2p,b3p,power
!      WRITE(6,'(a)') 'test: nuclear matter EoS'
!      DO i=1,60
!        call enernm(0.005D0*i,half,zero,e0,efmasp,efmasn)
!        WRITE(6,'(f6.3,1pg14.6)') 0.005D0*i,e0
!      END DO
      call nucmat(rhonm,eovera,compr,effmas,f0,f0p,g0,g0p,
     &            asymme,dasymm,d2asym,skappa,0,-2)
      enucmat = eovera
      if(iwr.gt.-2) then
        OPEN(39,file='nucmat.dat')
        write(39,'(a/6(1pg20.12))')
     &    'NUCMAT: rhonm,eovera,compr,effmas,asymme,dasymm,kappa=',
     &     rhonm,eovera,compr,effmas,asymme,dasymm,skappa
c        CALL asymm_nucmat(rhonm,39)
        CLOSE(39)
      endif
c
c     printing the header
c
      if(iwr.gt.-2) then
        WRITE(6,'(2a//a,i4,a,i4)')
     &    '* * * * * * * * * * *      ',
     &    'hartree-fock with skyrme-force    * * * * * * * * * * * ',
     &    ' nucleus:  protons=',nprot,'   neutrons=',nneut
        WRITE(6,'(a,2i4,a,i5,a,2f10.5,a,f10.5,a,i5)')
     &   ' pairing: nmaxp/n=',nmaxp,nmaxn,'  ipair=',ipair,
     &   '  pairfp/n=',pairfp,pairfn,'  temprt=',temprt,
     &   ' ilnrsp=',ilnrsp
        IF(ioptsk==13) THEN            ! Fayans functional
          WRITE(6,'(8(a,1pg15.7))') 'avp=',avp,'h1v=',h1vp,'h2vp=',h2vp     
          WRITE(6,'(8(a,1pg15.7))') 'avm=',avm,'h1v=',h1vm,'h2vm=',h2vm
          WRITE(6,'(8(a,1pg15.7))')
     &           'sigmavol=',sigmavol,'rho0fy=',rho0fy
          WRITE(6,'(8(a,1pg15.7))')
     &           'asp=',asp,'hsp=',hsp,'hgradp=',hgradp
          WRITE(6,'(8(a,1pg15.7))')
     &           'asm=',asm,'hsm=',hsm,'hgradm=',hgradm
          WRITE(6,'(8(a,1pg15.7))')
     &           'hcoul=',hcoul,'sigmacoul=',sigmacoul
          WRITE(6,'(8(a,1pg15.7))') 'fxi=',fxi,'h1xi=',h1xi,
     &           'hgradxi=',hgradxi,'gammaxi=',gammaxi,'rho0xi=',rho0xi
          WRITE(6,'(2(a,1pg15.7),a,i5)') 'delfxi=',delfxi,
     &             ' ecstab=',ecstab,' isovectormodel=',isovectormodel
          WRITE(6,'(8(a,1pg15.7))') 'c0tau=',c0tau,'c1tau=',c1tau
          WRITE(6,'(8(a,1pg15.7))') 'c0nabj=',c0nabj,'c1nabj=',c1nabj,
     &           'c0tens=',c0tens,'c1tens=',c1tens
        ELSE 
          WRITE(6,'(a,9f10.2/a,9f10.4,a,f10.4/2(a,f10.3),a,f10.5)') 
     &      ' t0-t8=',t0,t1,t2,t3,t4,t5,t6,t7,t8,
     &      ' x0-x8=',x0,x1,x2,x3,x4,x5,x6,x7,x8,' power=',power,
     &      ' tyuk=',tyuk,'  xyuk=',xyuk,'  ayuk=',ayuk
          IF(hgradxi.NE.0D0) WRITE(6,'(a,4(1pg13.5))')
     &     ' fxi,h1xi,gammaxi,hgradxi=',fxi,h1xi,gammaxi,hgradxi
        END IF
        WRITE(6,'(a,i3,a,f8.5,a,i3,a,i4,a,2(1pg9.2),a,i2,a,i3)') 
     &       ' numeric parameters:  ngrid=',ngrid,'  rstep=',rstep,
     &          '  itpr=',itpr,'  itmax=',itmax,'  dabr,dabe=',
     &          dabr,dabe,' iwr=',iwr,'  ifplot=',ifplot
        WRITE(6,'(a,0pf4.1,10(a,i2))') 
     &       '    coulomb-factor=',acoul,'  ifex=',ifex,' iftm=',iftm,
     &       '  ifrhoc=',ifrhoc,'  ifenrg=',ifenrg,'  ifrso=',ifrso,
     &       '  ifgsc=',ifgsc,'  ifzpe=',ifzpe,
     &       '  ifdfms=',ifdfms,'  imfold=',imfold,'  ioptsk=',ioptsk
        WRITE(6,'(a,0pf4.1,a,f9.2,a,f9.2)') 
     &     '   wave function iteration with: x0damp=',x0damp,
     &     '  e0damp=',e0damp,'artificial wall=',wall
!       write(06,1) nprot,nneut,
!     &             nmaxp,nmaxn,ipair,pairfp,pairfn,temprt,ilnrsp,
!     &             t0,t1,t2,t3,
!     &             t4,t5,t6,t7,t8,x0,x1,x2,x3,x4,x5,x6,x7,x8,power,
!     &             tyuk,xyuk,ayuk,
!     &             ngrid,rstep,itpr,itmax,dabr,dabe,iwr,ifplot,
!     &             acoul,ifex,iftm,ifrhoc,ifenrg,ifrso,ifgsc,ifzpe,
!     &             ifdfms,imfold,ioptsk,
!     &             x0damp,e0damp,wall
!    1 format(1h1,20('* '),5x,'hartree-fock with skyrme-force',5x,
!     &       20(' *')///1h0,36x,'protons=',i4,32x,'neutrons=',i4/
!     &       1h0,4x,'pairing: nmaxp/n=',2i4,'  ipair=',i5,
!     &              '  pairfp/n=',2f10.5,'  temprt=',f10.5,
!     &              ' ilnrsp=',i5,//
!     &       5x,' t0-t8=',9f10.2/
!     &       1x,' x0-x8=',9f10.4,' power=',f10.4/
!     &       1x,' tyuk=',f10.3,'  xyuk=',f10.3,'  ayuk=',f10.5//
!     &       5x,'numeric parameter:','  ngrid=',i3,'  rstep=',f8.5,
!     &          '  itpr=',i3,'  itmax=',i3,'  dabr=',
!     &          1pg9.2,' dabe=',1pg9.2,' iwr=',i2,'  ifplot=',i3/
!     &       5x,'coulomb-factor=',0pf4.1,'  ifex=',i2,' iftm=',i2,
!     &          '  ifrhoc=',i2,'  ifenrg=',i2,'  ifrso=',i2,
!     &          '  ifgsc=',i2,'  ifzpe=',i2,
!     &          '  ifdfms=',i2,'  imfold=',i2,'  ioptsk=',i2/
!     &       5x,'wave function iteration with: x0damp=',0pf9.4,
!     &          '  e0damp=',f9.2/
!     &       5x,'artificial wall with height=',f9.2/)
       if(ifgsc.eq.0) write(06,'(5x,a)')
     &   'the linear response part is not required'
       if(ioptsk.eq.5 .or. ioptsk.eq.6 .or. ioptsk.eq.11) 
     &   write(06,'(4(a,1pg15.7))')
     &    '  c1=',c1,'  c1p=',c1p,'  c4=',c4,'  c4p=',c4p
       if(ioptsk.eq.8 .or. ioptsk.eq.11) 
     &   write(06,'(4(a,f12.4))')
     &    '  b3tilde=',b3tilde,
     &    '  b12=',b12,'  b12p=',b12p
       if(ioptsk.eq.10 .or. ioptsk.eq.11) 
     &   write(06,'(4(a,f12.4))') 
     &    '  riftm=',riftm,'  rho0pr=',rho0pr,
     &    ' powerp=',powerp,'  ecstab=',ecstab
       if(ifform.lt.0)
     & write(06,'(a,i3,a,2(/5x,5(a,f10.6)))')
     &   '    ifform=',ifform,
     &   '. thus density-dependent nucleon folding is used with:',
     &   '   apro0=', apro0,'  apro11=',apro11,'  apro12=',apro12,
     &   '  bpro11=',bpro11,'  bpro12=',bpro12,
     &   '   aneu0=', aneu0,'  aneu11=',aneu11,'  aneu12=',aneu12,
     &   '  bneu11=',bneu11,'  bneu12=',bneu12
      endif
      if(ifplot.gt.0) then
       inquire(11,opened=top11)
       if(.not.top11) then
        open(unit=11,status='unknown',file='for011')
        open(unit=23,status='unknown',file='for023')
       endif
!       IF(ipair==6 .OR. ipair==8) WRITE(11,'(a,3i5,f6.2)')
!     &       '# iterhfb,iternat,modnat,pairenhance=',
!     &       iterhfb,iternat,modnat,pairenhance
       inquire(12,opened=top12)
       if(.not.top12) then
         counter=0
         open(unit=12,status='unknown',file='for012')   ! testprint iteration
c         write(12,'(a,a)')
c     &   '  Z    N         E_B        asurf      R   E_pair   ',
c     &   '  ZPE   r / E/A  R/A**1/3'
c         write(12,'(a,a)')
c     &   '  Z    N     A^{-1/3}   I    esurf2   esurf3  '
         write(12,'(a,a)')
     &   '  iter E    r  dE   dr   delta_p   delta_n'
       endif 
       inquire(13,opened=top13)
       if(.not.top13) then
         open(unit=13,status='unknown',file='for013')   ! testprint iteration
         write(13,'(a)')
     &   ' iter    spe_protons   '
       endif
       inquire(14,opened=top14)
       if(.not.top14) then
         open(unit=14,status='unknown',file='for014')   ! testprint iteration
         write(14,'(a)')
     &    '  iter    spe_neutrons   '
       endif
       write(11,'(/40a)')  (' *',i=1,40)
       write(11,*) title
c      write(12,'(/40a//4(a,i5)/a,3f10.3,f10.1,f10.3/a,4f10.5,a,f8.5)')
c    &   (' *',i=1,40),
c    &   ' nucleus:  nprot=', nprot,'  nneut=',nneut,
c    &            '  nmaxp=', nmaxp,'  nmaxn=',nmaxn,
c    &   ' force:    t0-t4=',t0,t1,t2,t3,t4,
c    &   '           x0-x3=',x0,x1,x2,x3,'  power=',power
c      write(12,'(a,f6.3,a,f6.2,a,f6.2)')
c    &   ' step:    x0damp=',x0damp,' e0damp=',e0damp,' addnew=',addnew
c      write(12,'(a)') '  iter  energy       ewid       detot'
       write(11,'(/4(a,i5))')
     &   ' nucleus:  nprot=', nprot,'  nneut=',nneut,
     &            '  nmaxp=', nmaxp,'  nmaxn=',nmaxn
       IF(ioptsk==13) THEN            ! Fayans functional
         WRITE(11,'(8(a,1pg15.7))') 'avp=',avp,'h1v=',h1vp,'h2vp=',h2vp     
         WRITE(11,'(8(a,1pg15.7))') 'avm=',avm,'h1v=',h1vm,'h2vm=',h2vm
         WRITE(11,'(8(a,1pg15.7))')
     &           'sigmavol=',sigmavol,'rho0fy=',rho0fy
         WRITE(11,'(8(a,1pg15.7))')
     &           'asp=',asp,'hsp=',hsp,'hgradp=',hgradp
         WRITE(11,'(8(a,1pg15.7))')
     &           'asm=',asm,'hsm=',hsm,'hgradm=',hgradm
         WRITE(11,'(8(a,1pg15.7))')
     &           'hcoul=',hcoul,'sigmacoul=',sigmacoul
         WRITE(11,'(8(a,1pg15.7))') 'fxi=',fxi,'h1xi=',h1xi,
     &           'hgradxi=',hgradxi,'gammaxi=',gammaxi,'rho0xi=',rho0xi
         WRITE(11,'(3(a,1pg15.7))') 'delfxi=',delfxi,
     &             ' delh1xi=',delh1xi,' delhgradxi=',delhgradxi
         WRITE(11,'(3(a,1pg15.7))') '  fxip=',fxip,
     &             '   h1xip=',h1xip,'   hgradxip=',hgradxip
         WRITE(11,'(3(a,1pg15.7))') '  fxin=',fxin,
     &             '   h1xin=',h1xin,'   hgradxin=',hgradxin
         WRITE(11,'(1(a,1pg15.7),a,i5)') 
     &             ' ecstab=',ecstab,' isovectormodel=',isovectormodel
         WRITE(11,'(8(a,1pg15.7))') 'c0tau=',c0tau,'c1tau=',c1tau
         WRITE(11,'(8(a,1pg15.7))') 'c0nabj=',c0nabj,'c1nabj=',c1nabj,
     &           'c0tens=',c0tens,'c1tens=',c1tens
         WRITE(11,'(10x,a,2f10.5)')  ' h2m=',amfc0
       ELSE
         write(11,'(a,3f10.3,f10.1,f10.3/a,5f10.5/a,f8.5,a,2f10.5)')
     &   ' force:    t0-t4=    ',t0,t1,t2,t3,t4,
     &   '           x0-x3,b4p=',x0,x1,x2,x3,b4p,
     &   '           power=    ',power,'  amfc0=',amfc0
         if(tyuk.ne.zero) write(11,'(12x,a,f10.3,a,f10.5,a,f8.5)')
     &               'tyuk=',tyuk,'   xyuk=',xyuk,'   ayuk=',ayuk
         if(abs(t5)+abs(t6)+abs(t7)+abs(t8).gt.zero)
     &     write(11,'(10x,a,4f10.3/10x,a,4f10.5)')
     &             ' t5-t8=',t5,t6,t7,t8,
     &             ' x5-x8=',x5,x6,x7,x8
         if(ioptsk.eq.5 .or. ioptsk.eq.6 .or. ioptsk.eq.11) 
     &     write(11,'(4(a,f8.5))')
     &    '             c1=',c1,'  c1p=',c1p,
     &    '  c4=',c4,'  c4p=',c4p
         if(ioptsk.eq.8 .or. ioptsk.eq.11) 
     &     write(11,'(4(a,f12.4))')
     &    '        b3tilde=',b3tilde,
     &    '  b12=',b12,'  b12p=',b12p
         if(ioptsk.eq.11)  write(11,'(4(a,f12.4))')
     &    '            tte=',tte,'  tto=',tto
         if(ioptsk.eq.10 .or. ioptsk.eq.11) 
     &     write(11,'(a,f6.2,2(a,f12.4),a,f7.3)') 
     &     '           riftm=',riftm,'   rho0pr=',rho0pr,
     &     ' powerp=',powerp,' ecstab=',ecstab
       END IF
       write(11,'(10x,a,i3,a,2f8.3,a,f8.3,a,i3,/
     &            10x,a,l1,3(a,f7.3)/5(a,i5)/10x,5(a,i5)/9x,5(a,i5))')
     &            ' ipair=',ipair,' pairfp/n=',pairfp,pairfn,
     &            ' temprt=',temprt,
     &            ' ilnrsp=',ilnrsp,' oldcut=',tphpair,
     &            ' cutfac=',cutfac,' hfbwid=',hfbwidth,' acoul=',acoul,
     &   ' cases:     ifex=',ifex,'  iftm=',iftm,'  ifzpe=',ifzpe,
     &            ' ifdfms=',ifdfms,' ifrhoc=',ifrhoc,
     &            ' ifrso=',ifrso,' ifgsc=',ifgsc,' ifenrg=',ifenrg,
     &            ' ifplot=',ifplot,' imfold=',imfold,' ioptsk=',ioptsk
       write(11,'(a,f8.3,a,i5/a,f6.3,a,f6.2,2(a,1pg9.2)/10x,2(a,i5))')
     &   ' grid:     rstep=',rstep,' ngrid=',ngrid,
     &   ' step:     x0damp=',x0damp,' e0damp=',e0damp,' dabr=',dabr,
     &              ' dabe=',dabe,
     &              ' itmax=',itmax,' itpr=',itpr
       IF(pairaddfacp.NE.0D0) WRITE(11,'(a,1pg13.5)')
     &   ' proton pairing constraint:',pairaddfacp
       IF(pairaddfacn.NE.0D0) WRITE(11,'(a,1pg13.5)')
     &   ' neutron pairing constraint:',pairaddfacn
       IF(ipair==6 .OR. ipair==8) write(11,'(a,3i5,f8.2,l3)')
     &   ' HFB: iterhfb,iternat,modnat,pairenhance,thfbfac=',
     &   iterhfb,iternat,modnat,pairenhance,thfbfac
       if(ifform.lt.0)
     &  write(11,'(a,i3,a,2(/5x,5(a,f10.6)))')
     &   '    ifform=',ifform,
     &   '. thus density-dependent nucleon folding is used with:',
     &   '   apro0=', apro0,'  apro11=',apro11,'  apro12=',apro12,
     &   '  bpro11=',bpro11,'  bpro12=',bpro12,
     &   '   aneu0=', aneu0,'  aneu11=',aneu11,'  aneu12=',aneu12,
     &   '  bneu11=',bneu11,'  bneu12=',bneu12
c       write(11,'(//)')
      endif
      if(iwr.gt.-2) then
       t0eff = 1.333333333333D0*(t0a-0.5*t0b)
       x0eff = 1.333333333333D0*(t0b-0.5*t0a)
       write(06,'(a,2g13.6)') ' effective t0,x0*t0=',t0eff,x0eff
      endif
c
c     optionaly testprint of initial wave functions
c
      if(iwr.ge.1) call priwf(06)
c
c     test densities for 'itmax=0'
c
      if(itmax.eq.0)  then
       call densit(weight,v2occ)
       write(06,91) '  rho:',(rho(i),i=1,ngrid)
       write(06,92) (rho(i),i=ngridp,ngrid2)
       write(06,91) '  tau:',(tau(i),i=1,ngrid)
       write(06,92) (tau(i),i=ngridp,ngrid2)
       write(06,91) '  d1r:',(d1r(i),i=1,ngrid)
       write(06,92) (d1r(i),i=ngridp,ngrid2)
       write(06,91) '  d2r:',(d2r(i),i=1,ngrid)
       write(06,92) (d2r(i),i=ngridp,ngrid2)
       write(06,91) '   sj:',( sj(i),i=1,ngrid)
       write(06,92) ( sj(i),i=ngridp,ngrid2)
       write(06,91) '  sjs:',(sjs(i),i=1,ngrid)
       write(06,92) (sjs(i),i=ngridp,ngrid2)
   91  format(1h0,a,20(/1x,10g13.5))
   92  format(20(/1x,10g13.5))
       return
      endif
c
c     initialise counter for stopping iteration
c
      ebind  = zero
      rmst   = zero
      ierror = 0
c
c     initialize second pairing constraint for Lipkin-Nogami
c
      xlamp2 = gpp/4.0
      xlamn2 = gpn/4.0
c
c     here starts the iteration
c
ctest      if(ipair.eq.3 .and. ifread.ne.1) then
ctest      if((ipair.eq.3 .or. ipair.ge.5) .and. ifread.ne.1) then
      iparac = ipair
      if((ipair.eq.3 .or. ipair.ge.6) .and. ifread.ne.1) then
!         iparac = 0       ! start Delta-appr. to save time & to stabilize
         if(ipair.eq.6) iparac = 5
         if(ipair.eq.8) iparac = 7
         gppsav = gpp
         gpnsav = gpn
!      else
!         iparac = ipair
      endif
c
      thfb = ipair==6 .OR. ipair==8
!      WRITE(11,*) 'before iteration: tpair,ipair=',
!     &    tpair,tpairFy,ipair,thfb
      IF(thfb) THEN
        fxisav=fxi
        pairfpsav=pairfp
        pairfnsav=pairfn
!        WRITE(11,'(a,3i5,f7.2)')'# iterhfb,iternat,modnat,pairenhance=',
!     &   iterhfb,iternat,modnat,pairenhance
      END IF
!      WRITE(6,*) ' thfb:',
!     &         thfb,iparac,fxisav,pairfpsav,pairfnsav,gppsav,gpnsav
      phpair = 1D0
      efep   = zero
      efen   = zero
      efaddp = zero
      efaddn = zero
c
      do 1000 iter=1,itmax
      itersav=iter
      iteration_number=iter
      IF(ipair>4 .AND. iter>iterbcs) THEN
        iparac=ipair
      ELSE
        iparac=ipair-1
      END IF
      IF(thfb) THEN
        IF(iter<iterhfb) THEN
          pairfacact=1D0+pairenhance*MAX(0,iterhfb-iter)/iterhfb
          fxi = fxisav*pairfacact
          pairfp=pairfpsav*pairfacact
          pairfn=pairfnsav*pairfacact
          gpp = gppsav*pairfacact
          gpn = gpnsav*pairfacact
          IF(thfbfac) hfbfacact=  MIN((1D0*iter)/iterhfb,1D0) ! 1D0-MAX(0,iterhfb-iter)/iterhfb
        ELSE
          fxi = fxisav
          pairfp=pairfpsav
          pairfn=pairfnsav
          gpp = gppsav
          gpn = gpnsav
          hfbfacact=1D0
        END IF
      END IF
!      WRITE(6,'(a,i5,3(1pg13.5))') ' iparac,fxi...=',
!     &         iparac,fxi,pairfp,pairfn
c
c     iterative steps
c
c                                           compute densities
      IF(.not.tnopot) THEN
        call densit(weight,v2occ)   
      END IF
!      IF(MOD(iter,100).EQ.0) WRITE(*,*) 'PAIRCR-1:',paircr(1:2*ngrid)
!      WRITE(6,*) 'before SKYPOT: iter,tpair,ipair=',
!     &    iter,tpair,tpairFy,ipair
c                                           mean-field hamiltonian
      if(.not.tnopot .or. ifpoic.eq.0) then
        call hfbpot   !???
        call skypot(iter)
      else
c                                           point Coulomb only
        call pure_C()
      endif
c                                           stepping wave function
      if(x0damp.eq.zero) then
       call welfu
      else
       if(iparac.eq.6.or.iparac.eq.8) then
         call hfbpot
         call lagpar(iter)
         IF(thfbnatural .AND. thfbnatprint .AND.
     &      iter>iterhfb .AND. MOD(iter,modnat)==0) THEN
           CALL HFB_natural(1,efep,1D0*nprot,-1)
           CALL HFB_natural(2,efen,1D0*nneut,-1)
           symcondsav = symcondi
         END IF
       endif
!       WRITE(6,'(a,5(1pg13.5))') 'iter,pots:',iter,pot(2),potls(2),
!     &      bmass(2),d1b(2)
       call welfdr(iter)
      endif
!        call hfbpot
c                                           pairing iteration
      IF(MOD(iter,100)==100000) THEN
        efep   = -20D0 !zero
        efen   = -20D0 !zero
        efaddp = zero
        efaddn = zero
      END IF
      epp = zero
      epn = zero
      if(iparac.eq.3 .and. ilnrsp.eq.1) then
	call lnresp                         ! prepare response terms for L.N.
      else
        rspln(1) = zero
        rspln(2) = zero
        suv2r(1) = one
        suv2r(2) = one
      endif
      if(ifpair(1).ne.0) then
!       WRITE(6,'(a)') 'proton pairing'
       if(ipair.eq.4 .or. ipair.le.1 .or. iparac.eq.0) then
         call pairol(spenrg(1),weight(1),ph(1),nprot,nmaxp,gpp,
     &               efep,epp,deltp,iparac,sp,itpr,epspr,iwr)
         if(ipair.eq.4) epp  = zero
       elseif(ipair.ge.2 .and. ipair.le.3) then
         if(iter.gt.4) then
          rhoscp = semrho(efep,pot(1),bmass(1))
          gpp  = gavrg(3-ipair,nprot,nneut,1,rhoscp,1,numstp,
     &                 effmas,0.58*0.960d0)              ! estimate surfp=0.96
         endif
         call bcslip(spenrg(1),ph(1),1,nmaxp,nprot,gpp,deltp,efep,
     &               xlamp2,weight(1),wgtobs(1),eppold,epp,entotp,
     &               rspln(1),suv2r(1),ipair-2,iwr)
       elseif(ipair.ge.5) then
clflpar         v0prot = pairfp
clflpar         call parmat(spenrg(1),pairmp,v0prot,efermp,1,kmaxn,iwr)
clflpar         call flpair(spenrg,deltaf(1),weight(1),ph,pairmp,one*nprot,
clflpar     &               efermp,epp,nmaxp,kmaxn,abs(itpr),epspr,iwr)
        if(iparac.eq.5.or.iparac.eq.7) then
          tpair = (iparac.ge.7)
          call hfbpot
!      IF(MOD(iter,100).EQ.0) WRITE(*,*) 'PAIRCR-2:',paircr(1:2*ngrid)
          if(iter.le.2) then   !   kick off from pairing breakdown
            do n=1,nmaxp+nmaxn
              deltaf(n) = max(deltaf(n),1D-4)
            enddo
          endif
        endif
c        efermp = zero
        itpract=ABS(itpr)
        IF((iparac.eq.6.or.iparac.eq.8).AND.iter>20000) itpract=-itpract
        call pairdn(spenrg,deltaf(1),weight(1),ph(1),phpair(1),
     &              ttthrm(1),v2occ(1),
     &              cutfac,one*nprot,efep,epp,deltp,deltp2,eqminp,gpp,
     &              efaddp,temprt,nmaxp,kmaxn,1,nblockp,
     &              itpract,epspr,iwr,tlipnog)
       endif
      endif
      if(ifpair(2).ne.0) then
!       WRITE(6,'(a)') 'neutron pairing'
       if(ipair.eq.4 .or. ipair.le.1 .or. iparac.eq.0) then
         call pairol(spenrg(nminn),weight(nminn),ph(nminn),nneut,
     &               nmaxn,gpn,efen,epn,deltn,iparac,sn,itpr,epspr,iwr)
         if(ipair.eq.4) epn  = zero
       elseif(ipair.ge.2 .and. ipair.le.3) then
         if(iter.gt.4) then
          rhoscn = semrho(efen,pot(ngridp),bmass(ngridp))
          gpn  = gavrg(3-ipair,nprot,nneut,0,rhoscn,1,numstn,
     &                 effmas,0.58*0.960d0)              ! estimate surfp=0.96
         endif
         call bcslip(spenrg(nminn),ph(nminn),1,nmaxn,nneut,gpn,deltn,
     &             efen,xlamn2,
     &             weight(nminn),wgtobs(nminn),epnold,epn,entotn,
     &             rspln(2),suv2r(2),ipair-2,iwr)
       elseif(ipair.ge.5) then
clflpar         v0neut = pairfn
clflpar         call parmat(spenrg(nminn),pairmn,v0neut,efermn,2,kmaxn,iwr)
clflpar         call flpair(spenrg(nminn),deltaf(nminn),weight(nminn),ph,
clflpar     &               pairmn,one*nneut,
clflpar     &               efermn,epn,nmaxn,kmaxn,abs(itpr),epspr,iwr)
c         efermn = zero
         call pairdn(spenrg(nminn),deltaf(nminn),weight(nminn),
     &               ph(nminn),phpair(nminn),ttthrm(nminn),v2occ(nminn),
     &               cutfac,one*nneut,efen,epn,deltn,deltn2,eqminn,gpn,
     &               efaddn,temprt,nmaxn,kmaxn,nminn,nblockn-nmaxp,
     &               itpract,epspr,iwr,tlipnog)
         IF(MOD(iter,10000)==100001) THEN
           DO isub=1,5
           call pairdn(spenrg(nminn),deltaf(nminn),weight(nminn),
     &               ph(nminn),phpair(nminn),ttthrm(nminn),v2occ(nminn),
     &               cutfac,one*nneut,efen,epn,deltn,deltn2,eqminn,gpn,
     &               efaddn,temprt,nmaxn,kmaxn,nminn,nblockn-nmaxp,
     &               itpract,epspr,iwr,tlipnog)
           END DO
         END IF
       endif
      endif
      IF(.NOT.tpairFy) THEN
        epair  = epp+epn
        IF(tepair .AND. ipair>4 .AND. ecstab.NE.0D0) epair=epp0+epn0
      END IF
!      WRITE(6,*) 'after PAIRDN: tpair=',tpair,tpairFy
      if(iwr.ge.0 .and. iparac.ge.2)
     & write(6,'(a,4g12.4)')
     &  ' semicl. rhoscp/n, efep/n=',rhoscp,rhoscn,efep,efen
c
      eshell = eshl(1)+eshl(2)
cold      etotal = eshell+epair+enrear
      etotal = ekin+epotsk+epair
ctest      write(6,'(a,3g12.4)') ' epairs=', epair,epp,epn
      ebiold = ebind
      ebind  = etotal/anumb
      detot  = ebiold-ebind
c
c     optionally HFB iteration in configuration  space
c
!      IF((ipair==6 .OR. ipair==8) .AND. iter>iterhfb .AND. 
      IF((iparac==6 .OR. iparac==8) .AND. iter>iterhfb .AND. 
     &    MOD(iter,modnat)==0 .AND. thfbnatural) THEN
        CALL HFB_natural(1,efep,1D0*nprot,iternat)
        CALL HFB_natural(2,efen,1D0*nneut,iternat)
      END IF
c
c     check termination
c
      if (iter.gt.1) then
c
       rmsold = rmst
       rmst   = zero
       r      = zero
       do i=2,ngrid
         r      = rstep+r
         rvol   = (r*r)*(r*r)
         rmst   = rvol*(rho(i)+rho(i+ngrid)) + rmst
       enddo
       rmst   = sqrt(fourpi*rmst*rstep/anumb)
       drmst  = abs(rmst-rmsold)
c
       IF(iparac==6 .OR. iparac==8 .AND. thfbnatural 
     &    .AND. thfbnatprint .AND. tprint96) THEN
          WRITE(96,'(2i5,10(1pg13.5))')
     &    iter,iparac,detot,drmst,symmcond,symcondi,symcondf,symcondsav
         FLUSH(96)
       END IF
c                                               set finally lip.nog. scheme
       if(iparac.ne.ipair .and.
     &    ((drmst.lt.200.0*dabr .and. abs(detot).lt.200.0*dabe)
     &     .or. iter.gt.iterbcs)  ) then
         iparac = ipair
         gpp = gppsav
         gpn = gpnsav
       elseif(drmst.lt.dabr.and.abs(detot).lt.dabe.and.iter.gt.40) then
         IF(.NOT.thfb .OR. (thfb .AND. iter>iterhfb)) goto 1009
       endif
      endif
c
c     optionally print intermediate results
c
      if(iwr.ge.0) write(06,'(//a,i4,a/)')
     &              '   results of the ',iter,' th. iteration:'
      if(iwr.ge.-1) then
        write(06,'(1x,i4,2(a,g11.4,a,g12.5),a,g11.4)')
     &              iter,' th. iteration: e-diff.=',detot,
     &              ' energy=',etotal,'  r-diff.=',drmst,
     &              ' radius=',rmst,'  enrear=',enrear
      endif
!      write(12,'(1x,i5,5(1pg12.4))')
!     &  iter,etotal,ewid,detot,epp,epn
      if (iwr.ge.0)
     & write(06,90) eshell,eshl(1),eshl(2),epair,epp,epn,
     &              efep,efen,etotal,ebind
      if (iwr.ge.1) call drpot(06)
      if(iwr>0) write(12,'(1x,i5,6(1pg15.7))')
     &   iter,etotal,rmst,abs(detot),abs(drmst),deltp,deltn
      if(iwr.ge.0 .and. mod(iter,10).eq.0) then
c      if(iwr.ge.0 .and. mod(iter,10).eq.100) then
       write(13,'(1x,i5,137(1pg11.3))')
     &   iter,(spenrg(i),i=1,nmaxp)
c     &   iter,(wfnorm(savlr(i*ngrid+1,1)),i=1,nmaxp)
       write(14,'(1x,i5,146(1pg13.5))')
     &   iter,(spenrg(i),i=nmaxp+1,nmaxto)
c     &   iter,(wfnorm(savlr(i*ngrid+1,1)),i=nmaxp+1,nmaxto)
      endif
!
!      IF(iter==634) CALL test_hpsi()

!        IF(tprint96) WRITE(23,'(3i4,87(1pg11.3))')
!     &   iter,nprot,nneut,(spenrg(i),i=1,nmaxp)
!        IF(tprint96) WRITE(24,'(3i4,87(1pg11.3))')
!     &   iter,nprot,nneut,(spenrg(i),i=nmaxp+1,nmaxp+nmaxn)
c
 1000 continue
c
c     no convergence if the routine comes to that point
c
      ierror = 2
      write(06,'(a,2i4,a,i5,a)') 'Z,N=',nprot,nneut,
     &  ': no convergence within ',itmax,
     &                     ' iterations !'
      WRITE(6,'(a,5(1pg13.5))') 'etotal,drmst,detot=',
     & etotal,drmst,detot
c
c     the end of the iteration
c
 1009 continue
!      WRITE(6,'(a,3(1pg13.5))') ' I. epairs:',epair,epp,epn
      IF(tprint96) WRITE(23,'(1x/1x)')
      IF(tprint96) WRITE(24,'(1x/1x)')
ctest      write(6,*) ' HAFOSU: end of iteration'
      IF(ipair==6 .OR. ipair==8) THEN
        WRITE(11,'(11x,a,3i5,2(1pg13.5))') 'nr. iterations,symm.cond.=',
     &         nprot,nneut,itersav,symcondi
      ELSE
        WRITE(11,'(11x,a,i5)') 'nr. iterations=',itersav
      END IF
      if(iwr.ge.-2 .and. ierror.ne.2)
     & write(06,'(a,2i4,a,i4,a)') ' z,n=',nprot,nneut,
     &   ', normal termination after',iter,' iterations'
      if(iwr.ge.1) call priwf(06)
!      IF(thfbnatural .AND. KMAXN<300) CALL hfb_matrix2()
!      CALL HFB_natural(2,efermhfb,1D0*nneut)
c
c     optionally saving on file 19
c
      if(idisc.ge.0 .and. ifsave.eq.1) call buffio('o')
c
c     compute c.m. correction
c
      if(ifzpe.gt.0 .or. ifzpe.le.-3) then
       call zpecor
       call rvariance
       if(iwr.ge.-1) then
         write(6,'(a,100(/i5,3g12.4))')
     &    ' cm corrections for protons:',
     &     (jsp(n),spenrg(n),specm(n),delcm(n),n=1,nmaxp)
         write(6,'(a,100(/i5,3g12.4))')
     &    ' cm corrections for neutrons:',
     &     (jsp(n),spenrg(n),specm(n),delcm(n),n=nminn,nmaxto)
       endif
       if(iwr.ge.-1) write(6,'(a,g13.5)') ' ZPE=',zpe
       if(ifzpe.gt.-4) etotal = etotal-zpe
       ebind  = etotal/anumb
      else
       zpe    = zero
      endif
c
c     compute kinetic energy
c
      ekinpr  = zero
      ekinne  = zero
      r       = zero
      do i=2,ngrid
        r       = rstep+r
        r2      = r*r
        ekinpr  = r2*amfc(1)*tau(i)       + ekinpr
        ekinne  = r2*amfc(2)*tau(ngrid+i) + ekinne
      enddo
      ekinpr  = frpirs*ekinpr
      ekinne  = frpirs*ekinne
      ekin    = ekinpr+ekinne
      if(iwr.ge.0) write(6,'(a,2g12.4)')
     &  ' kinetic energy, diagonal zpe=',ekin,ekin/anumb
c
c    test penalty potential
c
c      call cutpot(spenrg(1),hpswid(1),1)
c      call cutpot(spenrg(nmaxp+1),hpswid(nmaxp+1),2)
c
c    test folding
c
      call nucmat(rhonm,eovera,compr,effmas,f0,f0p,g0,g0p,
     &            asymme,dasymm,d2asym,skappa,0,-2)
ctest      write(6,*) ' HAFOSU: after NUCMAT'
c      if(iwr.ge.0) then
c        call filter0(rho,aux2)     
c        write(38,'(a)') '#  test filtering',
c     &      '#   ir    r     in   out '
c        do ir=1,ngrid
c          write(38,'(1x,i5,f8.2,2(1pg13.5))') 
c     &      ir,rfact(ir),rho(ir),aux2(ir)
c        enddo
c      endif
c
c     compute surface energy
c
      if(iwr.ge.-1 .and. .not.tnopot) call enrlda
!      if(iwr.ge.-2 .and. .not.tnopot) call enrlda
c
c     final print of fields
c
!      WRITE(6,'(a,3(1pg13.5))') ' II. epairs:',epair,epp,epn
      if(ipair.ge.5 .and. iwr.ge.-1) then
!      if(ipair.ge.5 .and. iwr.ge.-2) then
clflpar        v0prot = 260.0
clflpar        call parmat(spenrg(1),pairmp,v0prot,efermp,1,kmaxn,iwr+1)
clflpar        v0neut = 300.0
clflpar        call parmat(spenrg(nminn),pairmn,v0neut,efermn,2,kmaxn,iwr+1)
        call hfbpot
       if(ifpair(1).ne.0)
     &  call pairdn(spenrg,deltaf(1),weight(1),ph(1),
     &              phpair(1),ttthrm(1),v2occ(1),
     &              cutfac,one*nprot,efep,epp,deltp,deltp2,eqminp,gpp,
     &              efaddp,temprt,nmaxp,kmaxn,1,nblockp,
     &              abs(itpr),epspr,0,tlipnog)
       if(ifpair(2).ne.0)
     &  call pairdn(spenrg(nminn),deltaf(nminn),weight(nminn),ph(nminn),
     &              phpair(nminn),ttthrm(nminn),v2occ(nminn),
     &              cutfac,one*nneut,efen,epn,deltn,deltn2,eqminn,gpn,
     &              efaddn,temprt,nmaxn,kmaxn,nminn,nblockn-nmaxp,
     &              abs(itpr),epspr,0,tlipnog)
!        write(6,'(a,10(/1x,5g12.4))')
!     &   ' dynamic gaps for protons:',(deltaf(n),n=1,nmaxp)
!        write(6,'(a,10(/1x,5g12.4))')
!     &   ' dynamic gaps for neutrons:',(deltaf(n),n=nminn,nmaxto)
!        write(6,'(a,10(/1x,5g12.4))')
!     &   ' pairing space for protons:',(phpair(n),n=1,nmaxp)
!        write(6,'(a,10(/1x,5g12.4))')
!     &   ' pairing space for neutrons:',(phpair(n),n=nminn,nmaxto)
!        write(6,'(a,2g12.4)') ' test: Fermi energies=',efep,efen
      endif
!      WRITE(6,'(a,3(1pg13.5))') ' III. epairs:',epair,epp,epn
!      WRITE(6,'(a,4(1pg13.5))') 
!     &          'test1 pair:',epair,epp+epn,epp,epn
!      WRITE(6,'(a,40(/5(1pg13.5)))') 'add potential:',ppadd(1:ngrid)
c
c     test Fayans pairing
c
!      if(ipair.ge.5 .and. (ioptsk==13 .and. .NOT.tpairFy)) then
      if(ipair.ge.5 .and. (ioptsk==13)) then
        WRITE(*,*) 'V0 pairing potential:'
        WRITE(6,'(5(1pg13.5))') ppadd(1:ngrid)
        WRITE(6,'(5(1pg13.5))') ppadd(ngrid+1:2*ngrid)
        tpairFy= .TRUE.
        WRITE(6,*) '!! now with Fayans pairing:'
        call skypot(0)
        call hfbpot
       if(ifpair(1).ne.0)
     &  call pairdn(spenrg,deltaf(1),weight(1),ph(1),
     &              phpair(1),ttthrm(1),v2occ(1),
     &              cutfac,one*nprot,efep,epp2,deltp,deltp2,eqminp,gpp,
     &              efaddp,temprt,nmaxp,kmaxn,1,nblockp,
     &              abs(itpr),epspr,0,tlipnog)
       if(ifpair(2).ne.0)
     &  call pairdn(spenrg(nminn),deltaf(nminn),weight(nminn),ph(nminn),
     &              phpair(nminn),ttthrm(nminn),v2occ(nminn),
     &              cutfac,one*nneut,efen,epn2,deltn,deltn2,eqminn,gpn,
     &              efaddn,temprt,nmaxn,kmaxn,nminn,nblockn-nmaxp,
     &              abs(itpr),epspr,0,tlipnog)
!        write(6,'(a,10(/1x,5g12.4))')
!     &   ' dynamic gaps for protons:',(deltaf(n),n=1,nmaxp)
!        write(6,'(a,10(/1x,5g12.4))')
!     &   ' dynamic gaps for neutrons:',(deltaf(n),n=nminn,nmaxto)
!        write(6,'(a,10(/1x,5g12.4))')
!     &   ' pairing space for protons:',(phpair(n),n=1,nmaxp)
!        write(6,'(a,10(/1x,5g12.4))')
!     &   ' pairing space for neutrons:',(phpair(n),n=nminn,nmaxto)
        write(6,'(a,2g12.4)') ' test: Fermi energies=',efep,efen
        WRITE(6,'(a,4(1pg13.5))') 
     &          'test2 pair:',epair,epp2+epn2,epp2,epn2
        WRITE(6,'(a,40(/5(1pg13.5)))')'add potential FY:',ppadd(1:ngrid)
        WRITE(6,'(5(1pg13.5))') ppadd(ngrid+1:2*ngrid)
        epn=epn2
        epp=epp2
        tpairFy=.FALSE.
        call skypot(0)
      endif
!      WRITE(6,'(a,3(1pg13.5))') ' IV. epairs:',epair,epp,epn

      if(iwr.ge.+1) call priwf(06)
      if(iwr.ge.-1) then
       call drpot(06)
       WRITE(6,'(a,4(1pg13.5))') 'E_kin,E_pot,E_pair,E_cor=',
     &       ekin,epotsk,epair,ecored
       write(06,90) eshell,eshl(1),eshl(2),epair,epp,epn,
     &              efep,efen,etotal,ebind,zpe
90     format(//5x,'shell-energies',3g17.9/5x,'pairing-energies',
     &        3g13.5,' fermi-energies',2g13.5//5x,'the total energy',
     &        g13.5,5x,'per particle',g13.5,'   zpe=',g12.4)
c
      if(iwr.ge.-1 .and. ipair.ge.2)
     & write(6,'(a,4g12.4)')
     &  ' semicl. rhoscp/n, efep/n=',rhoscp,rhoscn,efep,efen
      CALL FLUSH(6)
c
      endif
c
c     separate energy contributions
c
      call sepenr()
c
c     formfactors and radii
c
      if(.not.tnopot) call densit(weight,v2occ)      ! from raw density
c      call radius
      CALL FLUSH(6)
      call radius(rho,rho(1+ngrid),sjs,sjs(1+ngrid),
     &                  nprot,nneut,rstep,ngrid,
     &                  amfc(1),amfc(2),zpe,ekin,ifzpe,
     &                  iwr,ifrso,ifform,
     &                  formf,rhochr,  foc,fop,fon,fo,foc0,fow,
     &                  rmsc,rmsp,rmsn,rms,
     &                  rdmsc,rdmsp,rdmsn,rdms,
     &                  surfc,surfp,surfn,surf,
     &                  rsecc,rsecp,rsecn,rsec,
     &                  smomc,smomp,smomn,smom, rhominmax, str_squark)
!      write(6,'(a,2i5,6(1pg13.5))') ' test rms: Z,N,radii=',nprot,nneut,
!     &    rmsp,rmsn,rmsn-rmsp
      IF(tcentrdep) THEN
        OPEN(47,file='centrdep',position='append')
        rhominmax=0D0
        CALL evaluate_dip(rhochr,foc,rhomin,rhomax,ngrid,rstep,delq)
        rhominmax= (rhomax-rhomin)/rhomax
        rhoaver = foc(1)/(4D0*pi/3D0*rdmsc**3)
        rhominaver=(rhoaver-rhomin)/rhoaver
        rhoaver = fop(1)/(4D0*pi/3D0*rdmsp**3)
        rho0averp =(rhoaver-rho(1))/rhoaver
        rhoaver = fon(1)/(4D0*pi/3D0*rdmsn**3)
        rho0avern =(rhoaver-rho(1+ngrid))/rhoaver
        WRITE(47,'(2i4,14(1pg24.14))') nprot,nneut,
     &    15.3*(rdmsp/rsecp-1.0),15.3*(rdmsn/rsecn-1.0),
     &    15.3*(rdmsc/rsecc-1.0),ecoul,
     &    rhominmax,rhominaver,rho0averp,rho0avern,rho(1),rho(1+ngrid),
     &    rho(1)+rho(1+ngrid),rho(1+ngrid)-rho(1),
     &    rho(1+ngrid)/nneut-rho(1)/nprot,rmsc
        CLOSE(47)
      END IF
      IF(TFORMW) THEN
        IF(ifzpe>0) THEN
          cmwid  = 0.375*0.5*(amfc(1)+amfc(2))/((nprot+nneut)*zpe)
        ELSE
          cmwid=0.0
        END IF
!        WRITE(*,*) ' before RADII'
        OPEN(131,file='APVdensities.dat')
          WRITE(131,'(a)') title(1:30)
          WRITE(131,'(1pg20.12,a)') cmwid,'   cmwid'
          WRITE(131,'(2a)') '#  r   rho_prot   rho_neut  rhols_prot',
     &        'rhols_neut'
          WRITE(131,'((5(1pg20.12)))') ((ir-1)*rstep, 
     &       rho(ir),rho(ir+ngrid),sjs(ir),sjs(ir+ngrid),ir=1,ngrid)
        CLOSE(131)
        CALL FLUSH(6)
        call radii(rho,rho(1+ngrid),sjs,sjs(1+ngrid),
     &             rstep,ngrid,ifrso,0,cmwid,aux1,
     &             rmsw,rmsc2,rmsp2,rmsn2,rms2,
     &             rdmsw,rdmsc2,rdmsp2,rdmsn2,rdms2,
     &             surfw,surfc2,surfp2,surfn2,surf2,fow,str_squark)
        IF(nprot==82 .AND. nneut==126) THEN
          apv_weak = apv_funct(0.3978D0)*1D9
        ELSE IF(nprot==20 .AND. nneut==28) THEN
!          apv_weak = apv_funct(0.873)*1D9
          apv_weak = apv_funct(0.8733D0)*1D9
        ELSE
          apv_weak = 0D0
        END IF
        WRITE(*,*) 'A_PV=',apv_weak
      sump0 = 0D0
      sump2 = 0D0
      sumn0 = 0D0
      sumn2 = 0D0
      do ir=1,ngrid
        r = (ir-1)*rstep
        vol = (r*r)*rho(ir)
        sump0 = vol + sump0
        sump2 = vol*(r*r) + sump2
        vol = (r*r)*rho(ir+ngrid)
        sumn0 = vol + sumn0
        sumn2 = vol*(r*r) + sumn2
      end do
      write(6,'(a,2i5,6(1pg13.5))') ' test rms: Z,N,radii=',nprot,nneut,
     &    sqrt(sump2/sump0),sqrt(sumn2/sumn0),
     &    sqrt(sumn2/sumn0)-sqrt(sump2/sump0)
      END IF
ctest      write(*,*) ' HAFOSU: after RADIUS'
        CALL barrad(nprot,nneut)
      if(ifplot.gt.0) then
c        epair = epp+epn
!        IF(ifenrg.NE.1) CALL plot(11)
        CALL plot(11)
!        CALL test_hpsi()
        if(ipair.gt.4 .AND. ifplot.gt.0)
     &   write(11,'(a,2(1pg12.4))')
     &   ' pairing band EFADDP/N=',efaddp,efaddn
        IF(ifplot.GE.4) call prispe
c         counter=counter+1
         a13   = anumb**(1D0/3D0)
         asurf = (etotal/anumb-enucmat)*a13
        IF(iwr>0) THEN
c           write(12,'(2x,2i6,2(1pg20.12),4(1pg13.5),3(1pg13.5))')
c     &      nprot,nneut,etotal,asurf,
c     &      rdms,epair,zpe,rms,rdms/a13,rms/a13/sqrt(0.6D0)
          write(13,'(1x,i5,37(1pg11.4))')
     &      counter,(spenrg(i),i=1,nmaxp)
          write(14,'(1x,i5,46(1pg11.4))')
     &      counter,(spenrg(i),i=nmaxp+1,nmaxto)
        END IF
!        call momden
c         call fbtpot()
c         call enrwig()
      endif
      if(ifplot.ge.4) call prispd
ctest      write(*,*) ' HAFOSU: after PRISPD'
      if(ipair.eq.3) then
        if(.not.tnopot) call densit(wgtobs,v2occ)           ! from Lipkin-Nogami
c        call radius
         call radius(rho,rho(1+ngrid),sjs,sjs(1+ngrid),
     &                  nprot,nneut,rstep,ngrid,
     &                  amfc(1),amfc(2),zpe,ekin,ifzpe,
     &                  iwr,ifrso,ifform,
     &                  formf,rhochr,  foc,fop,fon,fo,foc0,fow,
     &                  rmsc,rmsp,rmsn,rms,
     &                  rdmsc,rdmsp,rdmsn,rdms,
     &                  surfc,surfp,surfn,surf,
     &                  rsecc,rsecp,rsecn,rsec,
     &                  smomc,smomp,smomn,smom, rhominmax, str_squark)
         if(ifplot.gt.0) call plot(11)
      endif
c
!      CALL dosph()
!      CALL dos()
c      CALL structure_factor
c      call testpnp
!       CALL pnp_overlap(v2occ,ph,phpair,spenrg,nmaxp,nmaxp,1,
!     &                       nprot,nneut,rho,paircr)
!       WRITE(73,'(a,1pg13.5)') 'epair:',
!     &     SUM(paircr(1:ngrid)**2*rfact(1:ngrid)**2)*frpirs

!       CALL pnshift_overlap(v2occ(nminn),ph(nminn),phpair(nminn),
!     &           spenrg(nminn),deltaf(nminn),nmaxn,nmaxn,2,nprot,nneut)
!       CALL pnp_overlap(v2occ(nminn),ph(nminn),phpair(nminn),
!     &           spenrg(nminn),nmaxn,nmaxn,2,nprot,nneut,
!     &           rho(ngrid+1),paircr(ngrid+1))
!       WRITE(73,'(a,1pg13.5)') 'epair:',
!     & SUM(paircr(ngrid+1:ngrid+ngrid)**2*rfact(1:ngrid)**2)*frpirs
!       DO j=0,30
!         CALL dos2ph(0D0,16D0,j)
!       END DO
c
c     optionally ground-state-correlations
c
ctest      write(*,*) ' HAFOSU: before GSC'
      ecored = zero
      etot0  = etotal
      if(ifgsc.eq.1) then
       call sumrul
      elseif(ifgsc.eq.2) then
       call conrpa
      elseif(ifgsc.eq.3) then
       call specor
      elseif(ifgsc.eq.4) then
       call resmte
      elseif(ifgsc.eq.5) then
       call cormom
      elseif(ifgsc.eq.6) then
       call spreso
      elseif(ifgsc.eq.7) then
       call lintdhf
      endif
c
c     test HFB
c
!      CALL HFB_iterat()
!      CALL HFB_natural(2,efermhfb,1D0*nneut)
!      efephfb=efep
!      partnump=aprot
!      CALL hfb_efermi(1,efephfb,partnump)
!      efenhfb=efen
!      partnumn=aneut
!      CALL hfb_efermi(2,efenhfb,partnumn)
c
c
c     for ifenrg.eq.1 compute energies ab initio
c
      if(ifenrg.eq.1) then
       tpairFy = ipair.ge.5 .and. (ioptsk==13)
       IF(pairaddfacp.NE.0D0) THEN
         pairfp=pairfp0
         fxip=fxip0
         h1xip=h1xip0
         hgradxip=hgradxip0
       END IF
       IF(pairaddfacn.NE.0D0) THEN
         pairfn=pairfn0
         fxin=fxin0
         h1xin=h1xin0
         hgradxin=hgradxin0
       END IF
       if(.not.tnopot) call densit(weight,v2occ)
       WRITE(11,'(a,3(1pg13.5))') 'before HFBPOT:',epp,epn,epair
!       CALL hfbpot
       WRITE(11,'(a,3(1pg13.5))') 'after HFBPOT:',epp,epn,epair
       call skypot(1)
       if(x0damp.eq.zero) then
        call welfu
       else
!        x0damp=1D-20
        call welfdr(itmax+1)
       endif
       eshell = eshl(1)+eshl(2)
cold       if(ifpair(1).eq.1 .and. ipair.lt.2) then
cold        epp   = -deltp*deltp/(gpp+xsmall)
cold       else
cold        epp   = zero
cold       endif
cold       if(ifpair(2).eq.1 .and. ipair.lt.2) then
cold        epn    = -deltn*deltn/(gpn+xsmall)
cold       else
cold        epn    = zero
cold       endif
       IF(IOPTSK.NE.13) THEN
         epair  = epp+epn
         IF(tepair .AND. ipair>4 .AND. ecstab.NE.0D0) epair=epp0+epn0
       END IF
c       epair  = epp+epn
c       if(ifzpe.eq.1.and.ifgsc.ne.1)
c    *                call width(1,1,deltp,deltn,-2)
c       zpe=0.e0
c       if(ifzpe.eq.1) zpe=11.67195*(nprot+nneut)/(wids0(1,2)+
c    *                      widp0(1,2)+widp0(1,2))
       detot  = etotal
cold       etotal = eshell+epair+enrear-zpe+ecored
       etotal = ekin+epotsk+epair-zpe+ecored
       etotal2 = ekin+epotsk+epair/2-zpe+ecored
       detot  = detot-etotal
       ebind  = etotal/anumb
       if(iwr.ge.-1) then
        WRITE(6,'(a,4(1pg13.5))') 'E_kin,E_pot,E_pair,E_cor=',
     &       ekin,epotsk,epair,ecored
        write(06,90) eshell,eshl(1),eshl(2),epair,epp,epn,
     &               efep,efen,etotal,ebind,zpe
        write(06,'(a,g12.5)') ' energy-difference (ifenrg):',detot
       endif
       if(ifplot.gt.0) call plot(11)
       IF(ifplot>3) THEN
         open(15,file='for015',position='append')
         write(15,'(2i5,f8.4,3(1pg20.12),1pg13.5)') 
     &     nprot,nneut,constr0,etotal,etotal2,anumb*rms**2,rms
         close(15)
       END IF
      endif
!   &  15.3*(rdmsn/rsecn-1.0)*rdmsn*fon(1)*3D0/(4D0*pi*rdmsn**3),
c
c     reset some values
c
ctest      write(*,*) ' HAFOSU: before RETURN'
      ngrid  = 0
!      WRITE(6,'(a,3(1pg13.5))') ' V. epairs:',epair,epp,epn
c
      return
      end
c-----initfc------------ part of the hartree-fock package ----------------------
c
      subroutine initfc
c
      include 'params.inc'
c
c     * * * * * * *      initialise force parameters      * * * * * * * *
c
c     initialises auxiliary force parameters from the force parameters
c     given via common /skypar/.
c     the given parameters are the skyrme parameters
c     't0-t8,tyuk', 'x0-x8,xyuk', 'ayuk' and 'power';
c     furthermore the pairing parameters 'pairfp' and 'pairfn',
c     and the coulomb strength 'acoul'
c     all the other parameters in common /skypar/ are derived quantities.
c
      include 'skypar.inc'
      include 'wavfun.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
      data iprib1/0/
c
cold      data x10up/1.0e10/
c
c-------------------------------------------------------------------------------
c
c     determine effective masses according to options 'ifdfms'
c     and 'ifzpe'. in fact, 'amfc' is half the inverse mass in
c     units mev/fm**2.
c
      if(ioptsk==12) then
       amfc0(1) = 197.3285879D0/4.753458673288399D0/2D0
       amfc0(2) = amfc0(1)
      elseif(ifdfms.eq.0) then
       amfc0(1) = 20.7525
       amfc0(2) = 20.7525
!       amfc0(1) = 20.73
!       amfc0(2) = 20.73
      elseif(ifdfms.eq.1) then
       amfc0(1) = 20.750
       amfc0(2) = 20.722
      elseif(ifdfms.eq.2) then
       amfc0(1) = 20.73553 ! 20.7339829978 !
       amfc0(2) = 20.73553 ! 20.7339829978 !
cskp       amfc0(1) = 20.73
cskp       amfc0(2) = 20.73
      elseif(ifdfms.eq.3) then
c       Mail from K.Rutz 15.4.1998 gives NIST-values as most precise:
c         hbar**2/(2m_n) = 20.721 259(6) MeV*fm**2 (NIST)
c         hbar**2/(2m_p) = 20.749 821(6) MeV*fm**2 (NIST)
c
       amfc0(1) = 20.749811 ! 20.74982   ! 20.7505   !20.747
       amfc0(2) = 20.721249 ! 20.72126   ! 20.7219   !20.719
      elseif(ifdfms.eq.4) then
c
c      the Pearson & Farine recipe
c
       amfc0(1) =  20.7505   
       amfc0(2) =  20.7219  
c       amfc0(1) = 20.73621
c       amfc0(2) = 20.73621
      elseif(SUM(amfc0)==0D0) then
        stop ' IFDFMS must be out of [0..3]'
      endif
      if(ifzpe.lt.1 .and. ifzpe.gt.-2) then
        if(ifdfms.ge.3) then
          shellN  = (1.5*anumb)**0.33333333
          effmfc  = 2.0/(shellN+0.33333333/shellN)      ! N.P. A422 (1984) 157
          amfc(1) = amfc0(1)-amfc0(1)*effmfc/anumb
          amfc(2) = amfc0(2)-amfc0(2)*effmfc/anumb
        else
          amfc(1) = amfc0(1)-amfc0(1)/anumb
          amfc(2) = amfc0(2)-amfc0(2)/anumb
        endif
      else
          amfc(1) = amfc0(1)
          amfc(2) = amfc0(2)
      endif
c
c
c     combinations of skyrme-parameters which occur frequently:
c     the simple assigments are (hopefully) self-explaining.
c     note: in case of yukawa folding, i.e. tyuk.ne.0, the
c           zero-range equivalent strength of the yukawa
c           potential is subtracted from the zer0-range two-
c           body force. thus the t0 still describes the overall
c           strength.
c
cold      aysm3  = one/(fourpi*ayuk*ayuk*ayuk)
c
      IF(ioptsk==13) THEN
        CALL  ini_fyfunc()
      ELSE IF(.NOT.ioptsk==12) THEN     ! override
        if((ioptsk.eq.10 .or. ioptsk.eq.11)
     &    .and. powerp.gt.1D-20) then
          tpowerp = .true.
        else
          tpowerp= .false.
        endif
        if(ioptsk==11 .and. riftm1==0D0) riftm1 = riftm
        IF(ABS(tte)+ABS(tto)>1D-10) THEN
          riftm1=1D0
          riftm=1D0
          iftm=1
        END IF
        if(ayuk.eq.zero) then
          tyuk = zero
          xyuk = zero
        endif
        t0b    = t0*(0.5+x0)-tyuk*(0.5+xyuk)
        t3b    = t3*(0.5+x3)/6.0
        t0a    = t0*(1.0+0.5*x0)-tyuk*(1.0+0.5*xyuk)
        if(ioptsk.eq.2 .or. ioptsk.eq.3) then
          t3a    = 0.0
          b3dir  = t3*(1.0+0.5*x3)/12.0
c          write(2,'(a,2g12.4)') ' INIT: b0dir,gamdir=',b0dir,gamdir
        else
          t3a    = t3*(1.0+0.5*x3)/6.0
          b3dir  = 0.0
        endif
c        write(*,*) ' t0a,...=',t0a,t0b,t3a,t3b
        if(exmass.ne.zero) then
          coefx  = t3b
          t3b    = 0.0                   ! to disable old isovector part
        endif
c        t7a    = t7*(1.0+0.5*x7)/24.0
c        t7b    = t7*(0.5+x7)/24.0
cold        t3v    = 0.25*t3
        t3as   = power*0.5*t3a
        t3bs   = power*0.5*t3b
        if(tpowerp) then
          write(6,'(a,1pg13.5)') ' isovector power=',powerp
          b3scal = t3/8.0D0 * 0.5D0      !  t3*3.0D0/16.0D0 *0.5D0
          b3vect = (t3+t3*x3*2.0D0)/24.0D0 *0.5D0  !  (t3+t3*x3*2.0D0)/16.0D0 *0.5D0  
          b3eff  = b3scal*(2.0D0+power)
          b3peff = powerp*b3vect
          b3peff2 = 2.0D0*b3vect
        endif
        t8zw   = t8/12.0
        t1p2   = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)     !    =  b1
        t1m2   = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)     !    =  b1p
        if(ioptsk.ne.4) then
          b4     = 0.5*t4
        else
          b4     = -0.5*b4p
          t4     = -b4p
          t1p2   = 0.5*t1m2
        endif
        if(b4p.eq.zero) b4p = b4
        if(ioptsk.eq.4 .and. iprib1.le.10) then
          write(6,'(a,1pg13.5)') ' B1=',t1p2-0.5*t1m2
  
          iprib1 = 1 + iprib1
        endif
c      if(tyuk.eq.zero) then
         t13p2  = 0.125*(1.5*t1+t1*x1*3.0+0.5*t2+t2*x2)
         t13m2  = 0.125*(3.0*t1+t1*x1*1.5-t2-t2*x2*0.5)
c        else
c         t13p2  = 0.0
c         t13m2  = 0.0
c        endif
         t1m20  = (t1-t2)/8D0 + 5D0*(tte-tto)/8D0   ! term on spin-orbit density
!      c0tens = (t1-2D0*t1*x1-t2-2D0*t2*x2)/16D0+ 5D0*(tte+3D0*tto)/16D0
!      c1tens = t1m20/2
!      t1p2x  = -c0tens+c1tens
        t1p2x  = (x1*t1+x2*t2)/8.0 - 5D0*tto/4D0
        t5p6   = (t5+t5*x5*0.5+t6+t6*x6*0.5)/24.0
        t5m6   = (0.5*t5+t5*x5-0.5*t6-t6*x6)/24.0
        t53p6  = (1.5*t5+t5*x5*3.0+0.5*t6+t6*x6)/96.0
        t53m6  = (3.0*t5+t5*x5*1.5-t6-t6*x6*0.5)/96.0
        t5m60  = (t5-t6)/48.0
        t5p6x  = (x5*t5+x6*t6)/48.0
        t5m60h = 0.5*t5m60
        t5p6xh = 0.5*t5p6x
        t53m6d = 2.0*t53m6
        if(ayuk.ne.zero) then
         ayms2  = 1.0/(ayuk*ayuk)
        else
         ayms2  = zero
        endif
        tya    = tyuk*(1.0+0.5*xyuk)
        tyb    = tyuk*(0.5+xyuk)
        IF(hgradxi.NE.0D0 .AND. hgradxi<HUGE(1D0)) CALL  ini_fyfunc()
      END IF                            ! end override
c
c     setting pairing parameters:
c      the input quantities are 'ipair', 'pairfp' and 'pairfn'.
c      their meaning is explained in the header of 'hafbody'.
c
c     assignement of the pairing parameters for the different
c     cases of 'ipair'
c
      IF(ipair==9) THEN
        tlipnog=.TRUE.
        ipair=5
      ELSE IF(ipair==10) THEN
        tlipnog=.TRUE.
        ipair=7
      ELSE
        tlipnog=.FALSE.
      END IF
      IF(tlipnog.AND.ecstab.NE.0D0) 
     &   STOP 'LN excludes stabilized pairing'
      if(nprot.ne.0 .and. nneut.ne.0) then        ! check real nucleus
        if(ipair.eq.3) then
          ibcs  = 0                           ! Lip.Nog.for GAVRG
        else
          ibcs  = 1                           ! BCS.for GAVRG
        endif
        if(ipair.eq.0) then                   !  constant Delta
          if(pairfp.eq.zero) then             ! protons
            deltp = 11.2/sqrt(anumb)
cpow3            deltp = 5.4*anumb**(-one/3)
          elseif(pairfp.gt.zero) then
            deltp = pairfp
          else
            stop ' "pairfp" < 0 is nonsense.'
          endif
          if(pairfn.eq.zero) then             ! neutrons
            deltn = 11.2/sqrt(anumb)
cpow3            deltn = 5.4*anumb**(-one/3)
          elseif(pairfn.gt.zero) then
            deltn = pairfn
          else
            stop ' "pairfn" < 0 is nonsense.'
          endif
c
        elseif(ipair.ge.1 .and. ipair.le.3) then   ! constant force G
          call nucmat(rhonm,eovera,compr,effmas,f0,f0p,g0,g0p,
     &                asymme,dasymm,d2asym,skappa,0,-2)
c         if(iwr.ge.0)
          if(iwr.ge.-1)
     &       write(6,'(a,g12.4)') ' in NUCMAT: m*/m=',effmas
          if(pairfp.eq.zero)  then                   ! protons
            numsta  = 0                            ! compute number of states
            do 2011 n=1,nmaxp
             js=jsp(n)
             ll=js/2
             if(2*ll.eq.js) js=js-1
             mmult = (js+1)/2
             numsta = mmult + numsta
 2011        continue
             numstp = numsta
            rhoscp = zero
            gpp  = gavrg(ibcs,nprot,nneut,1,rhoscp,1,numstp,
     &                   effmas,0.58*0.960d0)              ! estimate surfp=0.96
            if(ipair.eq.1) gpp  = 17.0/anumb               !! override
            deltp = 11.2/sqrt(anumb)                       ! initial value
            if(iwr.ge.-1)
     &        write(6,'(a,3g13.4)')
     &     ' GAVRG: numstp,gpp*A,rhoscp=',numstp,gpp*anumb,rhoscp
          elseif(pairfp.gt.zero) then
            gpp  = pairfp/anumb
          else
            stop ' "pairfp" < 0 is nonsense.'
          endif
          if(pairfn.eq.zero)  then                   ! neutrons
            numsta  = 0                            ! compute number of states
            do n=1,nmaxn
              js=jsp(n)
              ll=js/2
              if(2*ll.eq.js) js=js-1
              mmult = (js+1)/2
              numsta = mmult + numsta
            enddo
            numstn = numsta
            rhoscn = zero
            gpn  = gavrg(ibcs,nprot,nneut,0,rhoscn,1,numstn,
     &                   effmas,0.58*0.960d0)              ! estimate surfp=0.96
            if(ipair.eq.1)  gpn  = 23.0/anumb             !! override
            deltn = 11.2/sqrt(anumb)                       ! initial value
            if(iwr.ge.-1)
     &        write(6,'(a,3g13.4)')
     &     ' GAVRG: numstn,gpn*A,rhoscn=',numstn,gpn*anumb,rhoscn
          elseif(pairfn.gt.zero) then
            gpn  = pairfn/anumb
          else
            stop ' "pairfn" < 0 is nonsense.'
          endif
          if(iwr.ge.-1) write(6,'(5x,a,2g12.4)')
     &     ' Pairing force now GP,GN=',gpp,gpn
          if(ifplot.gt.0) write(11,'(5x,a,2g12.4)')
     &     ' Pairing force now GP,GN=',gpp,gpn
        elseif(ipair.eq.4) then                ! temperature
          if(temprt.le.zero)
     &        stop ' negative or zero temperature !'
          if(pairfp.ne.pairfn)
     &      write(6,'(a)')
     &       ' !! warning: T_prot .ne. T_neut ??'
          Deltp  = temprt             !       pairfp
          Deltn  = temprt             !       pairfn
        elseif(ipair.ge.5) then               !  full pairing force
          deltp = 11.2/sqrt(anumb)           !  fixed gap for first iteration
          deltn = 11.2/sqrt(anumb)
cpow3         deltp = 5.4*anumb**(-one/3)
cpow3         deltn = 5.4*anumb**(-one/3)
          if(ipair.ge.7) then
            if(pairfp.eq.zero) then             ! protons
              pairfp = 700.0                    !  1000.0
            endif
            if(pairfn.eq.zero) then             ! neutrons
              pairfn = 700.0                    !  1000.0
            endif
c            rho0pr = 0.16
          else
            if(pairfp.eq.zero) then             ! protons
              pairfp = 280.0                      !  260.0
            endif
            if(pairfn.eq.zero) then             ! neutrons
              pairfn = 310.0                      !  300.0
            endif
          endif
          do n=1,ksl2
            deltaf(n) = deltp
          enddo
        endif
        epspr  = dabe
      endif
c
c     switch force-option: normal parametrisation --> isky=8
c                          extend. parametrisation --> isky=10
c                          folded yukawa param. --> isky=9
c
      if(abs(t5)+abs(t6)+abs(t7)+abs(t8).gt.zero) then
       if(tyuk.eq.zero) then
        isky   = 10
       else
        write(06,'(a,a)') ' incompatible parameters: ',
     *            ' tyuk and t5-t8 not simultanously treated.'
        stop ' tyuk and t5-t8 not simultanously treated.'
       endif
      elseif(tyuk.ne.zero .and. ayuk.ne.zero) then
       isky   = 9
      else
       isky   = 8
      endif
      if(isky.gt.8 .and. ipair.eq.3)
     & stop ' extended Skyrme not compatible with Lipkin-Nogami'
c
c     option for the force translated to logical variables
c     for faster switching inside the do-loop.
c      'texcha' = exchange contribution to coulomb force
c      'tyukaw' = finite range two body force with yukawa folding
c      'texten' = extended parametrisation ( using t5-t8)
c
      texcha= (ifex.ne.0)
      tyukaw= (isky.eq.9)
      texten= (isky.gt.9)
      tssmas= (t9.ne.zero)
      tpair = (iparac.ge.7 .OR. ioptsk==13) 
      tpairFy = ioptsk==13 .or. (hgradxi.NE.0D0.AND.hgradxi<HUGE(1D0))
      IF(.NOT.tpairFy .AND. ioptsk==13) THEN
        pairfp = -fxi*pfacfy   
        pairfn = -fxi*pfacfy
        rho0pr = -rho0xi*fxi/h1xi
        WRITE(6,'(a,4(1pg13.5))') 'test Pair:',pfacfy,pairfp,rho0pr
      END IF
      texmas= (exmass.ne.zero)
      tb3til= (b3tilde.ne.zero)
      tb12  = (b12.ne.zero) .or. (b12p.ne.zero)
      tb11  = (b11.ne.zero) .or. (b11p.ne.zero)
      talpinv= alpinv.ne.0D0
      talpin2= alpin2.ne.0D0
      tdirpa= (ioptsk.eq.2) .or. (ioptsk.eq.3)
      tdens2= ioptsk.eq.9 .or. 
     &        (ioptsk.eq.11 .and. abs(pow7).gt.1D-20 
     &         .and. .not.talpin2  .and. betin2==0D0)
      tdenslow= ioptsk.eq.11 .and. abs(pow7).gt.1D-20 
     &         .and. .not.talpin2 .and. betin2.ne.0D0
      trmfki= ioptsk.eq.4
      tkindir= (ioptsk.eq.5 .or. ioptsk.eq.11)
     &        .and.  (abs(c1)+abs(c1p)+abs(c4)+abs(c4p)).gt.1D-20
      tkininv= ioptsk.eq.6
c      write(2,'(a,4(1pg13.5)/10x,4(1pg13.5))') 
c     &  ' TKINDIR: tkindir,c1,c1p,c4,c4p=',
c     &    tkindir,c1,c1p,c4,c4p,
c     &    c1*t1p2,c1p*t1m2,c4*b4,c4p*b4p
c      write(2,*)
c     &  ' tyuk,xyuk,ayuk,tyukaw=',tyuk,xyuk,ayuk,tyukaw
c
      return
      end
c-----initgr------------ part of the hartree-fock package --------------------
c
      subroutine initgr
c
      include 'params.inc'
c
c     * * * * * * * *      initialise grid parameters      * * * * * * * * *
c
c     initialisation of grid properties on common /grid/
c     and of shell ordering on common /shells/.
c
c     the shells are ordered as in the harmonic oscillator shell model.
c     within a shell, the ordering goes from highest total angular
c     momentum to lowest. thus the ordering starts as follows:
c      1s1/2, 1p3/2, 1p1/2, 1d5/2, 1d3/2, 2s1/2, 1f7/2, .....
c
c     the pointer field 'npla' for the linear storage of the wave functions
c     is also provided here:
c       the wave functions are stacked one above another in the field 'wf'
c       in portions of 'ngrid'.
c       they are ordered according to the oscillator shell model ordering,
c       precomputed in 'nsp'.
c
c
c     fields for the bookkeeping of the shell structure:
c       'nmaxfl' : describes the magic shells. number of magic shell is
c                  the second index. nmaxfl(1,nr.magic) gives the highes
c                  oscillator state in the shell and nmaxfl(2,nr.magic)
c                  gives the according maximum number of nucleons filled
c                  up the shell.
c       'nsp'    = nsp(n,js) where 'n' is the radial quantum number
c                  and 'js' the combined ang.mom. quantum number. the
c                  value of 'nsp' gives the number of the state (n,js)
c                  in harmonic oscillator ordering (see data statment).
c       'jsp'    = 'jsp(nho)' where 'nho' is the number of a state in
c                  harmonic oscillator ordering. the value 'jsp' gives
c                  the according combined ang.mom. 'js'.
c     the combined ang.mom. quantum number 'js' means:
c       js = 2*l   for  j=l-1/2
c       js = 2*l+1 for  j=l+1/2
c     where j=total ang.mom. and l=orbital ang.mom.
c
c     the magic numbers according to this h.o. ordering are
c
c                 (n+1)(n+2)(n+3)
c       N   (n) = ---------------*2  +  (2*n+4)
c        mag             6
c
c     where n=0,1,2,... is the principle quantum number of the
c     last fully closed shell. The extra (2n+4) comes from the
c     largest 'j' of the nect shell.
c
c
c     parameters for Nilsson model initialization
c
      parameter (hbom0=41.5, xkap0=0.05, xmu0=0.15)   ! parameters for Nillson
      parameter (edr=0.33333333333333D0)
c
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'constn.inc'
c
      dimension nsp(12,47),nmaxfl(2,0:25),jspsav(277)
c
      logical tdyninit                  ! switch to energy initialization
      data tdyninit/.FALSE./
c
cold      data nmaxfl/0,0, 1,2, 2,6, 3,8, 6,20, 7,28, 11,50, 16,82, 22,126,
cold     &            29,184, 37,258, 46,350, 56,462, 67,596, 79,754,
cold     &            92,938/
      data nmaxfl/
     &    0,     0,    1,     2,    2,     6,    3,     8,    6,    20,
     &    7,    28,   11,    50,   16,    82,   22,   126,   29,   184,
     &   37,   258,   46,   350,   56,   462,   67,   596,   79,   754,
     &   92,   938,  106,  1150,  121,  1392,  137,  1666,  154,  1974,
     &  172,  2318,  191,  2700,  211,  3122,  232,  3586,  254,  4094,
     &  277,  4648/
c
c     nsp  = max. # of radial wavefunctions for given JSP,
c            where JSP is the combined q.n. of J and L.
c                                           JSP   L    #states
cold      data nsp/1, 6,15,28,45,66,91,    !     1   0       2
cold     *         3,10,21,36,55,78,  99,  !     2   1       2
cold     *         2, 9,20,35,54,77,  99,  !     3   1       4
cold     *         5,14,27,44,65,90,  99,  !     4   2       4
cold     *         4,13,26,43,64,89,  99,  !     5   2       6
cold     *         8,19,34,53,76,   2*99,  !     6   3       6
cold     *         7,18,33,52,75,   2*99,  !     7   3       8
cold     *        12,25,42,63,88,   2*99,  !     8   4       8
cold     *        11,24,41,62,87,   2*99,  !     9   4      10
cold     *        17,32,51,74,      3*99,  !    10   5      10
cold     *        16,31,50,73,      3*99,  !    11   5      12
cold     *        23,40,61,86,      3*99,  !    12   6      12
cold     *        22,39,60,85,      3*99,  !    13   6      14
cold     *        30,49,72,         4*99,  !    14   7      14
cold     *        29,48,71,         4*99,  !    15   7      16
cold     *        38,59,84,         4*99,  !    16   8      16
cold     *        37,58,83,         4*99,  !    17   8      18
cold     *        47,70,            5*99,  !    18   9      18
cold     *        46,69,            5*99,  !    19   9      20
cold     *        57,82,            5*99,  !    20  10      20
cold     *        56,81,            5*99,  !    21  10      22
cold     *        68,               6*99,  !    22  11      22
cold     *        67,               6*99,  !    23  11      24
cold     *        80,               6*99,  !    24  12      24
cold     *        79,               6*99,  !    25  12      26
cold     *                          7*99,  !    26  13      26
cold     *        92,               6*99/  !    27  14      28
      data nsp/
     &    1,   6,  15,  28,  45,  66,  91, 120, 153, 190, 231, 276,
     &    3,  10,  21,  36,  55,  78, 105, 136, 171, 210, 253, 999,
     &    2,   9,  20,  35,  54,  77, 104, 135, 170, 209, 252, 999,
     &    5,  14,  27,  44,  65,  90, 119, 152, 189, 230, 275, 999,
     &    4,  13,  26,  43,  64,  89, 118, 151, 188, 229, 274, 999,
     &    8,  19,  34,  53,  76, 103, 134, 169, 208, 251, 999, 999,
     &    7,  18,  33,  52,  75, 102, 133, 168, 207, 250, 999, 999,
     &   12,  25,  42,  63,  88, 117, 150, 187, 228, 273, 999, 999,
     &   11,  24,  41,  62,  87, 116, 149, 186, 227, 272, 999, 999,
     &   17,  32,  51,  74, 101, 132, 167, 206, 249, 999, 999, 999,
     &   16,  31,  50,  73, 100, 131, 166, 205, 248, 999, 999, 999,
     &   23,  40,  61,  86, 115, 148, 185, 226, 271, 999, 999, 999,
     &   22,  39,  60,  85, 114, 147, 184, 225, 270, 999, 999, 999,
     &   30,  49,  72,  99, 130, 165, 204, 247, 999, 999, 999, 999,
     &   29,  48,  71,  98, 129, 164, 203, 246, 999, 999, 999, 999,
     &   38,  59,  84, 113, 146, 183, 224, 269, 999, 999, 999, 999,
     &   37,  58,  83, 112, 145, 182, 223, 268, 999, 999, 999, 999,
     &   47,  70,  97, 128, 163, 202, 245, 999, 999, 999, 999, 999,
     &   46,  69,  96, 127, 162, 201, 244, 999, 999, 999, 999, 999,
     &   57,  82, 111, 144, 181, 222, 267, 999, 999, 999, 999, 999,
     &   56,  81, 110, 143, 180, 221, 266, 999, 999, 999, 999, 999,
     &   68,  95, 126, 161, 200, 243, 999, 999, 999, 999, 999, 999,
     &   67,  94, 125, 160, 199, 242, 999, 999, 999, 999, 999, 999,
     &   80, 109, 142, 179, 220, 265, 999, 999, 999, 999, 999, 999,
     &   79, 108, 141, 178, 219, 264, 999, 999, 999, 999, 999, 999,
     &   93, 124, 159, 198, 241, 999, 999, 999, 999, 999, 999, 999,
     &   92, 123, 158, 197, 240, 999, 999, 999, 999, 999, 999, 999,
     &  107, 140, 177, 218, 263, 999, 999, 999, 999, 999, 999, 999,
     &  106, 139, 176, 217, 262, 999, 999, 999, 999, 999, 999, 999,
     &  122, 157, 196, 239, 999, 999, 999, 999, 999, 999, 999, 999,
     &  121, 156, 195, 238, 999, 999, 999, 999, 999, 999, 999, 999,
     &  138, 175, 216, 261, 999, 999, 999, 999, 999, 999, 999, 999,
     &  137, 174, 215, 260, 999, 999, 999, 999, 999, 999, 999, 999,
     &  155, 194, 237, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     &  154, 193, 236, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     &  173, 214, 259, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     &  172, 213, 258, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     &  192, 235, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     &  191, 234, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     &  212, 257, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     &  211, 256, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     &  233, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     &  232, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     &  255, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     &  254, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     &  999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999,
     &  277, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999/
c
c JSP: up to state nr.  277
      data jspsav/
     &  1,  3,  2,  5,  4,  1,  7,  6,  3,  2,  9,  8,  5,  4,  1,
     & 11, 10,  7,  6,  3,  2, 13, 12,  9,  8,  5,  4,  1, 15, 14,
     & 11, 10,  7,  6,  3,  2, 17, 16, 13, 12,  9,  8,  5,  4,  1,
     & 19, 18, 15, 14, 11, 10,  7,  6,  3,  2, 21, 20, 17, 16, 13,
     & 12,  9,  8,  5,  4,  1, 23, 22, 19, 18, 15, 14, 11, 10,  7,
     &  6,  3,  2, 25, 24, 21, 20, 17, 16, 13, 12,  9,  8,  5,  4,
     &  1, 27, 26, 23, 22, 19, 18, 15, 14, 11, 10,  7,  6,  3,  2,
     & 29, 28, 25, 24, 21, 20, 17, 16, 13, 12,  9,  8,  5,  4,  1,
     & 31, 30, 27, 26, 23, 22, 19, 18, 15, 14, 11, 10,  7,  6,  3,
     &  2, 33, 32, 29, 28, 25, 24, 21, 20, 17, 16, 13, 12,  9,  8,
     &  5,  4,  1, 35, 34, 31, 30, 27, 26, 23, 22, 19, 18, 15, 14,
     & 11, 10,  7,  6,  3,  2, 37, 36, 33, 32, 29, 28, 25, 24, 21,
     & 20, 17, 16, 13, 12,  9,  8,  5,  4,  1, 39, 38, 35, 34, 31,
     & 30, 27, 26, 23, 22, 19, 18, 15, 14, 11, 10,  7,  6,  3,  2,
     & 41, 40, 37, 36, 33, 32, 29, 28, 25, 24, 21, 20, 17, 16, 13,
     & 12,  9,  8,  5,  4,  1, 43, 42, 39, 38, 35, 34, 31, 30, 27,
     & 26, 23, 22, 19, 18, 15, 14, 11, 10,  7,  6,  3,  2, 45, 44,
     & 41, 40, 37, 36, 33, 32, 29, 28, 25, 24, 21, 20, 17, 16, 13,
     & 12,  9,  8,  5,  4,  1, 47/
c
c------------------------------------------------------------------------------
c
c     quick and dirty swicth to closed shells
c
c      IF(nprot==20) nmaxp=6
c      IF(nneut==20) nmaxn=6
c      IF(nprot==28) nmaxp=7
c      IF(nneut==28) nmaxn=7
c
c     prepare book-keeping fields 'wls' and 'zwj'
c
      do js=1,KJMX
        ll = js/2
        j2 = (js+1)/2
        j2 = j2+j2                       ! mutiplicity of 'js'
        zwj(js) = j2*d4pi
        IF(js.GT.ll+ll) THEN
          wls(js) = ll        
        ELSE
          wls(js) = -ll-1
        ENDIF
      enddo
      do n=1,MIN(MAX(KMAXN,KMAXP),277)
        jsp(n) = jspsav(n)
      enddo
      IF(KJMX.gt.47 .OR. KNMX.gt.12) THEN
        write(6,'(a)') ' using energy initialization '
        write(7,'(a)') ' using energy initialization '
        tdyninit = .TRUE.
c        IF(nprot.NE.nneut) 
c     &    STOP ' energy initialization only for NNEUT=NPROT'
      ENDIF
c
c     automatic setting the necessary number of s.p.states according
c     to given particle number 'nprot' or 'nneut'. it is done in case
c     of 'nmaxp'=0 or 'nmaxn'=0   (using the bookkeeping field 'nmaxfl').
c
      IF(.NOT.tdyninit) THEN
!        if(iparsh.lt.0 .or. iparsh.gt.3)
        if(iparsh.lt.0 .or. iparsh.gt.6)
     &    stop ' iparsh out of range 0..3 '
        if(nmaxp.eq.0) then
         do i=0,11
           isav   = i
           nmaxp  = nmaxfl(1,i)
           if(nprot.le.nmaxfl(2,i)) goto 19
         enddo
         stop ' table of "nmaxp" exhausted in subroutine "initgr"'
   19    continue
         nmaxp  = nmaxfl(1,isav+iparsh)
ccut         nmaxp  = nmaxfl(1,isav+iparsh+1)
ctst         if(iparsh.eq.2 .and. 2*nprot.le.nmaxfl(2,isav+1))
         if(ipair.lt.4 .and .
     &      iparsh.eq.2 .and. 1.8*nprot.le.1.0*nmaxfl(2,isav+1))  ! 1.78
     &    nmaxp  = nmaxfl(1,isav+1)
         if(iwr.ge.-1) write(6,'(a,i5)')
     &    ' automatic setting: NMAXP=',nmaxp
        endif
        if(nmaxp.gt.kmaxp)
     &    stop ' in initgr: more proton states than "kmaxp" allows '
        if(nmaxn.eq.0) then
         do i=0,11
           isav   = i
           nmaxn  = nmaxfl(1,i)
           if(nneut.le.nmaxfl(2,i)) goto 29
         enddo
         stop ' table of "nmaxn" exhausted in subroutine "initgr"'
   29    continue
         nmaxn  = nmaxfl(1,isav+iparsh)
ccut         nmaxn  = nmaxfl(1,isav+iparsh+1)
ctst         if(iparsh.eq.2 .and. 2*nneut.le.nmaxfl(2,isav+1))
         if(ipair.lt.4 .and.
     &      iparsh.eq.2 .and. 1.8*nneut.le.1.0*nmaxfl(2,isav+1))   ! 1.78
     &    nmaxn  = nmaxfl(1,isav+1)
ctst         if(iparsh.eq.2 .and. nneut.eq.132) nmaxn=37            ! special case
         if(iwr.ge.-1) write(6,'(a,i5)')
     &    ' automatic setting: NMAXN=',nmaxn
        endif
        if(nmaxn.gt.kmaxn)
     &    stop ' in initgr: more proton states than "kmaxp" allows '
      ENDIF
c
c     derived numeric parameters for shells and grid
c
      ngridp = ngrid+1
      ngrid2 = ngrid+ngrid
      nmaxto = nmaxp+nmaxn
      anumb  = nprot+nneut
      aprot  = nprot
      aprot  = max(aprot,xsmall)
      aneut  = nneut
      nminn  = nmaxp+1
      frpirs = fourpi*rstep
      rleng  = (ngrid-1)*rstep
      delq   = 3.1415926/rleng
      iqmax  = min(kfrm,ngrid)       ! ???  ngrid*3/4
c      iqmax  = ngrid-1
c
c     zentrifugal force and other functions of 'r'
c
      zentrf(1) = zero
      rinvs(1)  = zero
      rsquar(1) = zero
      rfact(1)  = zero
      do i=2,ngrid
        x=(i-1)*rstep
        rfact(i)  = x
        zentrf(i) = one/(x*x)
        rinvs(i)  = one/x
        rsquar(i) = x*x
      enddo
c
c     setting up the pointer field 'npla' for the wave functions.
c     npla(n+1,js) describes the first position in the wave function
c     field 'wf' (and 'wfs') for a state with radial quantum
c     number 'n' and combined ang.mom. 'js';
c     npla(1,js) carries the information on the maximal number
c     of radial 'n' for given 'js'.
c
!      IF(ifread==1) RETURN
c                                             zero pointer field
      do jx=1,kjmx+kjmx
      do n=1,knmx+1
        npla(n,jx) = 0
      enddo
      enddo
      IF(.NOT.tdyninit) THEN
c                                             compute pointer field
        do iso=1,2
        do js=1,kjmx
          jx     = (iso-1)*kjmx+js
          isoadd = ngrid*nmaxp*(iso-1)
          if(iso.eq.1) then
           nmax   = nmaxp
          else
           nmax   = nmaxn
          endif
          do n=1,knmx
            nact   = n-1
            nrn    = nsp(n,js)
            if(nmax-nrn.lt.0) goto 139
            npla(n+1,jx) = ngrid*(nrn-1)+1+isoadd
          enddo
  139     continue
          npla(1,jx) = nact
        enddo
        enddo
c                                               compute maximum 'js' and 'ph'
        do iso=1,2
          if(iso.eq.1) then
           nmax   = nmaxp
           isoadd = 0
          else
           nmax   = nmaxn
           isoadd = nmaxp
          endif
          jsm    = 0
          do i=1,nmax
            if(jsm-jsp(i).lt.0) jsm  = jsp(i)
            j2 = jsp(i)
            j2 = (j2+1)/2
            ph(i+isoadd) = j2*1D0
          enddo
          jsmax(iso) = jsm
        enddo
      ELSE
c
c       Here follows the initialization according to the Nilsson model
c       preferably used for huge nuclei.
c       States belonging to active space are estimated with the 
c       s.p.energy of the Nilsson model.
c       The cutoff energy 'eocc' is estimated from the Fermi energy in
c       an oscillator model plus a margin of 'iparsh' oscillator quanta.
c       States are sorted in order of 'js' the combined ang.mom. label.
c       For each 'js', as much states are taken as fit below the Fermi energy.
c
        write(6,'(a,3i5)')
     &   ' old params: nmaxp,nmaxn,nmaxto=',nmaxp,nmaxn,nmaxto
        hbom = hbom0*anumb**(-edr)
        xkm  = xkap0*xmu0
        nshell = 0
        eocc = max((1.4D0*aprot**edr + xparsh)*hbom,
     &             (1.4D0*aneut**edr + xparsh)*hbom)
        do iso=1,2
c          IF(iso.EQ.1) THEN
c            eocc = (1.4D0*aprot**edr + xparsh)*hbom
c          ELSE
c            eocc = (1.4D0*aneut**edr + xparsh)*hbom
c          ENDIF
          write(6,'(a,2(1pg13.5))') ' hbom,eocc=',hbom,eocc
          partnum = zero
          jslast=0
          do js=1,kjmx
            IF(js-jslast>5) EXIT
            jx     = (iso-1)*kjmx+js
            isoadd = ngrid*nmaxp*(iso-1)
            if(iso.eq.1) then
             nmax   = nmaxp
             nnumb  = nprot
            else
             nmax   = nmaxn
             nnumb  = nneut
            endif
            ll = js/2
            do n=1,knmx+1
              energy = hbom*( (2*(n-1)+ll) - xkap0*wls(js)
     &                       - xkm*(ll+ll*ll)/(1D0+xkm*ll)    )
!     &                       - xkm*(ll+ll*ll)                )
                write(6,'(a,3i5,1pg13.5)')
     &           ' iso,js,n,energy=',iso,js,n,energy
              IF(energy.LE.eocc) THEN
                IF(n.GT.knmx) THEN
                  write(6,'(a,3i5)')
     &             ' KNMX too small: iso,js,n=',iso,js,n
                  STOP ' KNMX too small '
                ENDIF
                jslast=js
                jsmax(iso) = js
                nshell = 1 + nshell
                IF(nshell.GT.KSL2) 
     &            STOP ' too many states required. enhance KSL2'
                npla(1,jx) = n
                npla(n+1,jx) = (nshell-1)*ngrid+1
                jsp(nshell) = js
                spenrg(nshell) = energy
                ph(nshell) = (js+1)/2  !D0
                partnum = ph(nshell)*2D0 + partnum
                write(6,'(a,2i5,1pg12.4)') 
     &           ' update: jsmax,nshell,partnum=',
     &                     jsmax(iso),nshell,partnum
              ENDIF
            enddo
          enddo
          IF(iso.EQ.1) THEN
            nmaxp  = nshell
          ELSE
            nmaxto = nshell
            nmaxn  = nmaxto-nmaxp
          ENDIF
          IF(partnum.LT.nnumb*one) THEN
            write(6,'(a,3i5,1pg13.5)')
     &        ' iso,jsmax,nshell,partnum=',
     &          iso,jsmax(iso),nshell,partnum
            STOP ' phase space not large enough to match part.number'
          ELSEIF(partnum.GT.nnumb*one) THEN
            ifpair(iso) = 1
          ELSE
            ifpair(iso) = 0
          ENDIF
          do js=1,jsmax(iso)
            jx     = (iso-1)*kjmx+js
            write(6,'(a,2i4,a)')  ' npla for iso,js=',iso,js,' :'
            write(6,'(20i5)')   (npla(n,jx),n=1,npla(1,jx)+1)
          enddo
        enddo
        nminn  = nmaxp+1
        write(6,'(a,3i5)')
     &   ' new params: nmaxp,nmaxn,nmaxto=',nmaxp,nmaxn,nmaxto
        write(6,'(a,100(/1x,10f6.1))')
     &   ' ph:',(ph(n),n=1,nmaxto)
        write(6,'(a,100(/1x,20i5))')
     &   ' jsp:',(jsp(n),n=1,nmaxto)
c  
c       determine a first occupation according to oscillator energies.
c  
        efep   =  zero
        deltini = 0.1D0
        call pairol(spenrg(1),weight(1),ph(1),nprot,nmaxp,gpp,
     &              efep,epp,deltini,0,sp,itpr,epspr,iwr)
        efen   =  zero
        deltini = 0.1D0
        call pairol(spenrg(nminn),weight(nminn),ph(nminn),nneut,
     &            nmaxn,gpn,efen,epn,deltini,iparac,sn,itpr,epspr,iwr)
        wgtobs=weight
        v2occ=weight
c
      ENDIF
c
      return
      end
c-----initwf------------ part of the hartree-fock package --------------------
c
      subroutine initwf
c
      include 'params.inc'
c
c     * * * * * * * *      initialise wave functions      * * * * * * * * *
c
c     initialises the wave function field 'wf' with harmonic oscillator
c     wave functions.
c     and initialises the pairing weights 'weight'.
c
c     information on the nuclear shell structure for initilaising the
c     pairing occupation weights is carried in the following fields:
c       nzahl(n)   = the maximum number of particles filled up to a given
c                    oscillator shell 'n'.
c       nmag(n)    = the boundaries of a magic shell which contains a
c                    given oscillator state 'n'.
c                    the last two digits are the upper bound of the shell
c                    and the first two digits the lower bound.
c
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'constn.inc'
c
      dimension nzahl(0:154),nmag(154)
c
c
c     nzahl(N) = # particles up to shell nr. N
c
      data nzahl/0,
     &           2,
     &           6,8,
     &           14,18,20,
     &           28,
     &           34,38,40,50,
     &           58,64,68,70,82,
     &           92,100,106,110,112,126,
     &           138,148,156,162,166,168,184,
     &           198,210,220,228,234,238,240,258,
     &           274,288,300,310,318,324,328,330,350,
     &           368,384,398,410,420,428,434,438,440,462,
     &           482,500,516,530,542,552,560,566,570,572,596,
     &           618,638,656,672,686,698,708,716,722,726,728,754,
     &           778,800,820,838,854,868,880,890,898,904,908,910,938,
     &           964,988,1010,1030,1048,1064,1078,1090,1100,1108,
     &               1114,1118,1120,1150,
     &           1178,1204,1228,1250,1270,1288,1304,1318,1330,1340,
     &               1348,1354,1358,1360,1392,
     &           1422,1450,1476,1500,1522,1542,1560,1576,1590,1602,
     &               1612,1620,1626,1630,1632,1666,
     &           1698,1728,1756,1782,1806,1828,1848,1866,1882,1896,
     &               1908,1918,1926,1932,1936,1938,1974/

c
c     nmag   = boundaries of magic shells
c              last three digits = upper end of interval
c              first three digits = lower end of interval
c
      data nmag/1001,2002,3003,3*4006,7007,4*8011,5*12016,6*17022,
     &          7*23029, 8*30037, 9*38046, 10*47056, 11*57067, 12*68079, 
     &          13*80092, 14*93106, 15*107121, 16*122137, 17*138154/
c
      data x9999/0.9999/
c
c     'omeg' is the standard shell spacing in mev,
c     'hmc'  is nucleon mass divided by hbar**2
c
      data omeg,hmc/41.,.024109348/
c
c-------------------------------------------------------------------------------
c
c     harm.osc.-w.f. as initial w.f.
c
c     iso-loop
c
      d2w    = rstep*sqrt(.5*hmc*omeg/anumb**.333333333)
      do 1010 iso=1,2
      isnpad = (iso-1)*kjmx
c
c     loop over combined angular momentum q.n. 'js'
c
      jsm    = jsmax(iso)
      do 10 js=1,jsm
      nmax   = npla(1,js+isnpad)
      if(nmax.ge.1) then
c                                          set qu.numbers
       ll     = js/2
       llp    = ll+1
       sig    = float(2*mod(ll,2)-1)
       npa1   = npla(2,js+isnpad)
c                                          w.f., unorthogonalised
       nold   = npa1
       wf(npa1) = zero
       na     = npa1+1
       ne     = npa1+ngrid-1
       x      = zero
       do i=na,ne  
         x      = x+d2w
         wf(i)  = exp(-x*x)*(x**llp)
       enddo
       if(nmax.gt.0) then
        do n=2,nmax
          na     = npla(n+1,js+isnpad)
          nold   = na-nold
          ne     = na+ngrid-1
          x      = -d2w
          do i=na,ne
            x      = x+d2w
            wf(i)  = wf(i-nold)*x*x
          enddo
          nold   = na
        enddo
       endif
c
c      schmidt-orthonormalisation
c
       call orthon(npla(2,js+isnpad),nmax)
c
c      derivatives of the wave functions on 'wfs'.
c
       do n=1,nmax
         ia     = npla(n+1,js+isnpad)
         call derv15(wf(ia),wfs(ia),sig,ngrid,rstep)
         IF(thfbfold) CALL Gauss_folding_m(wf(ia),wffold(ia),rstep,
     &                         ngrid,js,hfbwidth,MAXVAL(jsmax(1:2)))
       enddo
c
      endif
   10 continue
c
c     initialise the s.p. energies
c
      do n=1,nmaxp+nmaxn
        spenrg(n) = zero
      enddo
c
      IF(KJMX.LE.47 .AND. KNMX.LE.12) THEN
c
c       occupy 'weight' as pairing initialisation.
c       in case of filled shells pairing is switched off using the
c       field 'ifpair'.
c       in case of 'ipair'=-1 the field 'weight' is read
c       from input and pairing is then switched off using 'ifpair'.
c
        if(iso.eq.1) then
         nmax   = nmaxp
         nwgt   = 0
         nz     = nprot
        else
         nmax   = nmaxn
         nwgt   = nmaxp
         nz     = nneut
        endif
        if(nz.gt.1974) stop 
     &     ' too large particle number, enhance NMAG in INITWF'
        if(nmax.gt.0) then
         do  i=1,nmax
           wgtobs(i+nwgt) = zero
           weight(i+nwgt) = zero
         enddo
         ie     = nmax
         do n=1,ie
           if(nzahl(n)-nz.ge.0) goto 329
           wgtobs(n+nwgt) = one
           weight(n+nwgt) = one
         enddo
         WRITE(6,'(a,100i5)') ' iso,nprot,nneut,nmaxp,nmaxn,nzahl:',
     &               iso,nprot,nneut,nmaxp,nmaxn,nzahl(1:ie)
         stop ' too large proton or neutron number'
        endif
c  
  329   continue
        ia     = nmag(n)
        ie     = mod(ia,1000)
        ia     = ia/1000
        nold   = nzahl(ia-1)
        x      = float(nz-nold)/float(nzahl(ie)-nold)
        do n=ia,ie
          weight(n+nwgt) = x
        enddo
        if(weight(nmax+nwgt).lt.x9999 .and. ipair.ge.0) then
         ifpair(iso) = 1
        else
         ifpair(iso) = 0
        endif
        if(ipair.eq.-1) then
         weight(1+nwgt:nmax+nwgt)=0D0
         read(5,*) (weight(i+nwgt),i=1,nmax)
         write(06,'(a,i3,10(/5x,13f9.5))')
     &         '0    weights read in from input data for iso=',i,
     &                  (weight(i+nwgt),i=1,nmax)
        endif
      ENDIF
      v2occ=weight
c
c     initial guess for fermi energies
c
      efep   =  zero
      efen   =  zero
c
 1010 continue
c
      return
      end
c-----initgw------------ part of the hartree-fock package --------------------
c
      subroutine initgw
c
      include 'params.inc'
c
c     initialises the weights by reading from file 5
c     if 'gpp' or 'gpn' is zero.
c     else 'weight' remains unchanged.
c     this routine is used only in connection with restart from
c     scratch file 19.
c
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c-------------------------------------------------------------------------------
c
      do 1010 iso=1,2
c
      if(iso.eq.1) then
       nmax   = nmaxp
       nwgt   = 0
cold       nz     = nprot
      else
       nmax   = nmaxn
       nwgt   = nmaxp
cold       nz     = nneut
      endif
c
      if(ipair.eq.-1) then
       read(5,*) (weight(i+nwgt),i=1,nmax)
       write(06,'(a,i3,10(/5x,13f9.5))')
     &        '0    weights read in from input data for iso=',i,
     &                (weight(i+nwgt),i=1,nmax)
      endif
c
 1010 continue
c
      return
      end


c-----nucmat------- part of the Hartree-Fock package ------------------
c
      subroutine nucmat(rhonm,eovera,compr,effmas,f0,f0p,g0,g0p,
     &                  asymme,dasymm,d2asym,skappa,ifland,ifprnt)
c
c     evaluates the stationary state for nuclear matter from
c     minimising the energy as function of the density.
c
c     the force parameters are communicated via common/skypar/.
c
c     then a few nuclear matter parameters at the ninimum are
c     computed:
c       rhonm   = density
c       eovera  = energy per particle
c       compr   = compressibility
c       effmas  = effective mass
c       asymme  = asymmetry energy
c       dasymm  = d/drho of asymmetry energy
c       eosneut = neutron EoS at rho=0.01
c       deosneut= slope of neutron EoS at rho=0.01
c
c     computation stops after 'rhonm' is found for "ifland.eq.-1".
c     the Landau parameters are computed for "ifland.eq.1"
c       glan0   = landau paramer    g0
c       glan0p  = landau paramer    g0'
c       flan0   = landau paramer    f0
c       flan0p  = landau paramer    f0'
c
c     original version from joerg friedrich, rewritten 20.6.87
c
      include 'params.inc'
      include 'skypar.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'wavfun.inc'
      include 'constn.inc'
c
c      data pi/3.141592654/
c      data zero/0.0/
cold      data rhoini,drhoit,epsil,itrmax
cold     *    / 0.16 , 0.01 ,1.e-5,    10/
      data rhoini,drhoit,itrmax,drhocmp
     *    / 0.16 , 0.00001 ,    200, 1D-5/
      data deltpc/0.01/            ! variation of proton content
      data betab /5.742467999/     ! (3/5)*(3*pi^2)^*(2/3)  factor for tau_n/p
!      teovera(x) = 0.5D0*BC0*x+0.5D0*BC3*x**(1D0+power)
!     &            +b1*betat*x**1.666666667D0
!     &            +amfact*betat*x**0.666666667D0
      logical ttestparams
      data ttestparams/.false./
      real(8),parameter :: delrho=1D-4
c
c-----------------------------------------------------------------------
c
c     iteration to find the minimum energy
c
ctest      write(*,*) ' NUCMAT entered'
      eold   = 1.e10
      rhonm0 = rhoini
      do iter=1,itrmax
        rhonmp = rhonm0+drhoit
        call enernm(rhonmp,half,zero,eplus,efmasp,efmasn)
        call enernm(rhonm0,half,zero,energy,efmasp,efmasn)
        rhonmm = rhonm0-drhoit
        call enernm(rhonmm,half,zero,eminus,efmasp,efmasn)
        rhonm0 = -0.5*drhoit*(eplus-eminus)/(eplus-energy-energy+eminus)
     *           +rhonm0
        if(ifprnt.ge.0)
     &   write(06,'(a,i3,a,2g11.4,a,3g16.9)') ' it=',iter,' rho=',
     &      rhonmm+drhoit,rhonm0,' e+0-=',energy,eplus-eminus,
     &      (eplus-energy-energy+eminus)
        if(abs(eold-energy).lt.dabe) goto 29
        eold   = energy
      enddo
      write(06,'(a,i5,a)') ' no convergency after ',itrmax,
     *            ' steps in subroutine nucmat !'
c
   29 continue
      rhonm  = rhonm0
      fermp2 = (3.*pi*pi/2.*rhonm0)**(0.66666667)
ctest      write(*,*) ' NUCMAT before return'
      if(ifland.eq.-1) return
c
c     now evaluate the nuclear matter parameters at the minimum
c     (e/a nach ergaenzt gemaess krewald, diss., p. 65, gl (25) )
c
      call enernm(rhonm,half,zero,eovera,efmasp,efmasn)
      call enernm(rhonm+drhocmp,half,zero,eplus,efmasp,efmasn)
      call enernm(rhonm-drhocmp,half,zero,eminus,efmasp,efmasn)
      call enernm(rhonm     ,half,zero,energ0,efmasp,efmasn)
      compr = 9.0*rhonm0*rhonm0*(eminus-energy-energy+eplus)
     *                          /drhocmp**2
      effmas = efmasp
      hqqd2m = 0.5*(amfc0(1)+amfc0(2))
!
! test 'alpha'=power of density dependence
!
      IF(ttestparams) THEN
!        betat = betab*0.5D0**0.6666667D0
!        BC0=t0a-t0b/2D0
!        BC3=t3a-t3b/2D0
!        amfact = (amfc0(1)+amfc0(2))/2D0
!        termmas = b1*betat*rhonm0**1.666666667D0
!        termkin = amfact*betab*(rhonm0/2D0)**0.666666667D0
!        testleft = 9D0*(1D0+power)*eovera+compr
!        testright = (10D0-6D0-6D0*power)*termmas
!     &             +(3D0*(1D0+power)-2D0)*termkin
!        testeff = amfact/(amfact+b1*rhonm0)
!        t2eovera = 0.5D0*b0*rhonm0+0.5D0*BC3*rhonm0**(1D0+power)
!     &           +termmas+termkin
!        tzero = 0.5D0*b0*rhonm0
!     &         +(1D0+power)*0.5D0*BC3*rhonm0**(1D0+power)
!     &         +1.666666667*termmas+0.666666667*termkin
!        tcompr = 4.5D0*power*(1D0+power)*BC3*rhonm0**(1D0+power)
!     &           +10D0*termmas-2D0*termkin
!        counter = 10*b1*betat*rhonm0**1.666666667D0
!     &           -2D0*amfact*betat*rhonm0**0.666666667D0-compr
!        denom  = 9*eovera+6*b1*betat*rhonm0**1.666667D0
!     &          -3D0*amfact*betat*rhonm0**0.666666667D0
!        alphatest = counter/denom-1D0
!        write(6,'(a,4(1pg13.5))') '# test alpha:',
!     &      alphatest,compr,eovera,rhonm0
!        write(6,'(a,4(1pg13.5))') '# test NMP:  ',
!     &      t2eovera,tzero,tcompr
!        write(6,'(a,4(1pg13.5))') '# test pot/k:',
!     &      t2eovera-termkin,termkin
!        write(6,'(a,4(1pg13.5))') '# test b0,b3:',
!     &      b0,bC0,b3,b3/BC3
!        write(6,'(a,4(1pg13.5))') '# sub tests: ',
!     &      testleft,testright,testeff
        WRITE(6,'(a,6(1pg13.5))') '# c0 etc:',b0,b3,c0,c3,t0,t0*x0
!        CALL  params_from_NMP(alphatest,BD0,BD3,BD1,
!     &                        rhonm0,eovera,amfact,effmas,compr)
        write(6,'(a,4(1pg13.5))') '# params:   ',
     &      alphatest,BD0,BD3,BD1
      END IF
!      f0dir = effmas/hqqd2m/fermp2/6D0*compr-1D0
!      WRITE(6,'(a,1pg13.5)') ' F0 direct=',f0dir
!      WRITE(2,'(a,1pg13.5)') ' F0 direct=',f0dir
c
      call enernm(rhonm0,half-delrho,zero,energ2,efmasp,efmasn)
      call enernm(rhonm0,half,zero,energ0,efmasp,efmasn)
      call enernm(rhonm0,half+delrho,zero,energ1,efmasp,efmasn)
      asymme = (energ1-2.0*energ0+energ2)/(delrho*delrho)/8.0
!   kinetic symmetry energy
!      energ2 = ekinnm(rhonm0,half-delrho)
!      energ0 = ekinnm(rhonm0,half)
!      energ1 = ekinnm(rhonm0,half+delrho)
!      asymkn = (energ1-2.0*energ0+energ2)/(delrho*delrho)/8.0
!      WRITE(2,'(a,2(1pg13.5))') 'kine.symm.energy=',rhonm0,asymkn
!
      if(iwr.ge.01)
     & write(6,'(a,3g12.4/a,g12.4)')
     &  ' test ENERNM: eovera,energ0,energ1=',eovera,energ0,energ1,
     &  '              ekin=',hqqd2m*betab*rhonm0**(2./3.)
c       sum rule enhancement 'kappa' from slightly asymmetric matter
      amfcpp = amfc0(1)/efmasp
      amfcpn = amfc0(2)/efmasn
      amfcav = 0.5*(amfcpp+amfcpn)
      amfc0a = 0.5*(amfc0(1)+amfc0(2))
      damfcp = (amfcpn-amfcpp)/delrho                  !  /deltpc
      skappa = (amfcav-amfc0a+damfcp/4.0)/amfc0a
c
c     density dependence of asymmetry energy
c
      rhonmp = rhonm0 + delrho*10D0
      call enernm(rhonmp,half-delrho,zero,energ2,efmasp,efmasn)
      call enernm(rhonmp,half,zero,energ0,efmasp,efmasn)
      call enernm(rhonmp,half+delrho,zero,energ1,efmasp,efmasn)
      asymmp = (energ1-2.0*energ0+energ2)/(delrho*delrho)/8.0
      rhonmm = rhonm0 - delrho*10D0
      call enernm(rhonmm,half-delrho,zero,energ2,efmasp,efmasn)
      call enernm(rhonmm,half,zero,energ0,efmasp,efmasn)
      call enernm(rhonmm,half+delrho,zero,energ1,efmasp,efmasn)
      asymmm = (energ1-2.0*energ0+energ2)/(delrho*delrho)/8.0
      dasymm = (asymmp-asymmm)/(delrho+delrho)/10D0
! re-scale to comply with L_symm
      IF(TLSYMM) dasymm = 3D0*rhonm0*dasymm  !/asymmm
      d2asym = 9.0*rhonm0*rhonm0*(asymmp-2.0*asymme+asymmm)/
     &         delrho**2/100.0D0
c
c
      if(ifland.lt.1) return
c
      call lanpar(rhonm0,flan0,flan0p,glan0,glan0p,f0,f0p,g0,g0p,
     &            f1,f1p,c0fac)
c
c     the type 1 Landau parameters are still done explicitely here
c
      a1     = t1*(1.0+0.5*x1)
      a1p    = t1*(0.5+x1)
      a2     = t2*(1.0+0.5*x2)
      a2p    = t2*(0.5+x2)
      bs1     = 0.25*(a1+a2)
      b1p    = 0.25*(a1p-a2p)
      b1     = bs1-0.5*b1p
c      f1     = -3.0*rhonm*b1/(hqqd2m/effmas)
c      f1p    = 1.5*rhonm*b1p/(hqqd2m/effmas)
c      skappo = (f1p-f1)/(3.0+f1)
c
c     now print the results
c
c      if(ifprnt.gt.-2) then
      if(ifprnt.gt.-1) then
        f1b=3D0*(effmas-1D0)
        f1pb=3D0*effmas*(skappa+1D0)-3D0
        write(02,'(a,4(/1x,4(a,g12.5)))')
     &     ' the landau parameters (volume):',
     &     '      f0=', flan0,'     f0s=',flan0p,
     &     '      g0=', glan0,'     g0s=',glan0p,
     &     '      f0=',    f0,'     f0s=',   f0p,
     &     '      g0=',    g0,'     g0s=',   g0p,
     &     '      f1=',    f1,'     f1s=',   f1p,
     &     '     f1b=',   f1b,'    f1sb=',   f1pb,
     &     '  rat.f1=',   f1b/f1,' rat.f1s=',   f1pb/f1p
      endif
c
      return
      end
c-----nucmat------- part of the Hartree-Fock package ------------------
c
      subroutine anucmat(rhonm,procon,eovera)
c
c     evaluates the stationary state for asymmetric nuclear matter from
c     minimising the energy as function of the density.
c
c     the force parameters are communicated via common/skypar/.
c
c     then a few nuclear matter parameters at the ninimum are
c     computed:
c       rhonm   = density
c       procon =  proton content (0.5 for symm.matter)
c       eovera  = energy per particle
c
      include 'params.inc'
      include 'skypar.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'wavfun.inc'
      include 'constn.inc'
c
c      data pi/3.141592654/
c      data zero/0.0/
cold      data rhoini,drhoit,epsil,itrmax
cold     *    / 0.16 , 0.01 ,1.e-5,    10/
      data rhoini,drhoit,itrmax,drhocmp
     *    / 0.16 , 0.00001 ,    200, 1D-5/
      data deltpc/0.01/            ! variation of proton content
      data betab /5.742467999/     ! (3/5)*(3*pi^2)^*(2/3)  factor for tau_n/p
!      teovera(x) = 0.5D0*BC0*x+0.5D0*BC3*x**(1D0+power)
!     &            +b1*betat*x**1.666666667D0
!     &            +amfact*betat*x**0.666666667D0
      logical ttestparams
      data ttestparams/.false./
      real(8),parameter :: delrho=1D-4
c
c-----------------------------------------------------------------------
c
c     iteration to find the minimum energy
c
ctest      write(*,*) ' NUCMAT entered'
      eold   = 1.e10
      rhonm0 = rhoini
      IF(ipair==4) THEN
        tmpact=temprt
      ELSE
        tmpact=0D0
      END IF
      do iter=1,itrmax
        rhonmp = rhonm0+drhoit
        call enernm(rhonmp,procon,tmpact,eplus,efmasp,efmasn)
        call enernm(rhonm0,procon,tmpact,energy,efmasp,efmasn)
        rhonmm = rhonm0-drhoit
        call enernm(rhonmm,procon,tmpact,eminus,efmasp,efmasn)
        rhonm0 = -0.5*drhoit*(eplus-eminus)/(eplus-energy-energy+eminus)
     *           +rhonm0
        if(ifprnt.ge.0)
     &   write(06,'(a,i3,a,2g11.4,a,3g11.4)') ' it=',iter,' rho=',
     &      rhonmm+drhoit,rhonm0,' e+0-=',eplus,energy,eminus
        if(abs(eold-energy).lt.dabe) goto 29
        eold   = energy
      enddo
      write(06,'(a,i5,a)') ' no convergency after ',itrmax,
     *            ' steps in subroutine anucmat !'
c
   29 continue
      rhonm  = rhonm0
      eovera = energy
c
      return
      end
c-----lanpar------- part of the Hartree-Fock package ------------------
c
      subroutine lanpar(rhonm0,flan0,flan0p,glan0,glan0p,
     &                  f0,f0p,g0,g0p,f1,f1p,c0fac)
c
c     evaluates the Landau parameters for symmetric nuclear matter
c     at given density 'rhonm0'.
c
c     !!!  This routine works only in connection with old standard SHF
c
c     the force parameters are communicated via common/skypar/.
c
c     the Landau parameters are 
c       glan0   = landau paramer    g0
c       glan0p  = landau paramer    g0'
c       flan0   = landau paramer    f0
c       flan0p  = landau paramer    f0'
c     their dimensionless analogues are 'f0','f0p','g0', and 'g0p'
c       c0fac   = level density factor to make the parameters
c                 dimensionless
c
      include 'params.inc'
      include 'skypar.inc'
c?      include 'numpar.inc'
c      include 'grid.inc'
c?      include 'wavfun.inc'
c     include 'constn.inc'
c
      data pi/3.141592654D0/
      data zero,one/0.0,1.0/
c
c-----------------------------------------------------------------------
c
!      write(2,*) 'LANPAR entered'
c
c     preparation of recoupled parameters
c     t12p = cm(krewald), t45p = cm1(krewald), diss. p.64, gl.(22,23)
c
      powp1  = power+1.0
      betab  = 3./5.*(3.*pi*pi/2.)**(2./3.)
      hqqd2m = 0.5*(amfc0(1)+amfc0(2))
c?      if(ifzpe.ne.1) then
c?        hqqd2m = hqqd2m/(1.0-1.0/anumb)
c?      endif
c     hqqd2m = 197.3**2/2./936.
      t12p   = 3.*t1+5.*t2*(1.+4./5.*x2)
      t45p   = (3.*t5+5.*t6*(1.+4./5.*x6))/6.0
      fermp2 = (3.*pi*pi/2.*rhonm0)**(0.66666667)
      if(rhoswm.gt.zero) then
        tanfac = tanh(rhonm0/rhoswm)
      else
        tanfac = 1.0
      endif
      effmiv = one+rhonm0/hqqd2m
     &               *((t1p2-0.5*t1m2)*tanfac
     &                 +(t45p+(b6+0.5*b6p)*0.6*fermp2)*rhonm0/16.0)
cold      effmiv = one+rhonm0/hqqd2m
cold     &               *(t12p/16.0*tanfac
cold     &                 +(t45p+(b6+0.5*b6p)*0.6*fermp2)*rhonm0/16.0)
      effmas = one/effmiv
c
      hvonti = 0.5*t0*x0
     &        +1./12.*t3*x3*rhonm0**power
     &        +0.25*t1*x1*fermp2
     &        +0.25*(t5*x5/6.0)*rhonm0*fermp2
      flan0  = 3./4.*t0
     &        +t3/16.0*powp1*(power+2.0)*rhonm0**power
     &        +1./8.*fermp2*(t12p+9.*rhonm0*t5/6.0)
      glan0p =-t0/4.-rhonm0**power*t3/24.
     &        -fermp2/8.*(t1-t2+rhonm0*t5/6.0)
      flan0p = glan0p+t2*x2*fermp2/4.-hvonti
      glan0  = glan0p+t2*x2*fermp2/4.+hvonti
c
c     use formulae from RMP for the g-parameters 
c
ctest      betag = (3.*pi*pi/2.)**(2./3.)
      rhonmpow = rhonm0**power
ctest      rhonm23 = rhonm0**0.6666666666667D0
ctest      c0s = -t0*(0.5D0-x0)/4D0 -t3*(0.5D0-x3)/24D0*rhonmpow
ctest      c0T = riftm*(-t1*(0.5D0-x1)/8D0 +t2*(0.5D0+x2)/8D0)
ctest      c1s = -t0/8D0 -t3/48D0*rhonmpow
ctest      c1T = riftm*(-t1/16D0 +t2/16D0)
ctest      glan0n = 2D0*(c0s+c0T*betag*rhonm23)
ctest      glan0pn = 2D0*(c1s+c1T*betag*rhonm23)
      glan0d = -t0*(1-2*x0)/4D0-t3*(1-2*x3)*rhonmpow/24D0
      glan0d2=(-t1*(1-2*x1)+t2*(1+2*x2))*fermp2/8D0
      glan0pd=-t0/4D0-t3*rhonmpow/24D0
      glan0pd2=(-t1+t2)*fermp2/8D0
c
c     the f1 parameters
c
      c0tau  = 3D0*t1/16D0+0.25D0*t2*(1.25D0+x2)
      flan1  = -2D0*c0tau*fermp2
      c1tau  = -0.125D0*t1*(0.5D0+x1)+0.125D0*t2*(0.5D0+x2)
      flan1p = -2D0*c1tau*fermp2
      write(6,'(4(1pg13.5))') c0tau,flan1,c1tau,flan1p
c
      wfac   = sqrt(fermp2)*effmas/(pi*pi*hqqd2m)
      c0fac  = 1D0/wfac
      f0     = wfac*flan0
      f0p    = wfac*flan0p
      f1     = wfac*flan1
      f1p    = wfac*flan1p
      g0     = wfac*(glan0d+riftm*glan0d2)      !glan0
      g0p    = wfac*(glan0pd+riftm*glan0pd2)    !glan0p
      write(6,'(2a)')
     & '  g0_notensor    g0_tensor      g0     ',
     & '  g0p_notensor   g0p_tensor     g0p     '
      write(6,'(7f10.4)')
     &    wfac*glan0d,wfac*glan0d2,g0,
     &    wfac*glan0pd,wfac*glan0pd2,g0p
!      WRITE(2,*) 1D0/wfac
c
      return
      end

c-----enrlda--------------------------------------------------------
c
      subroutine enrlda
c
      include 'params.inc'
c
c     energy in LDA:
c     uses densities given in common /wavfun/
c     and computes the LDA-energy by integrating with 'enernm'
c     the surface energy is deduce by comparison with exact energy
c
      include 'numpar.inc'
      include 'grid.inc'
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'constn.inc'
c
c
c-----------------------------------------------------------------------
c
c     the surface energy by comparison with LDA
c
      elda    = 0.0
      do ir=2,ngrid
        rhop   = rho(ir)
        rhon   = rho(ir+ngrid)
        rhot   = rhop+rhon
        procon = rhop/rhot
        call enernm(rhot,procon,zero,energm,efmasp,efmasn)
        ract   = (ir-1)*rstep
        elda    = (ract*ract)*rhot*energm+ elda
      enddo
      elda     = elda*frpirs
      anumb3   = anumb**(one/3.0)
      esurf1   = (etotal+zpe-elda)/anumb * anumb3
c
c     the surface energy by comparison with infinite matter
c
!      if(iwr.ge.-1) then
      if(iwr.ge.-2) then
        call nucmat(rhonm,eovera,compr,effmas,f0,f0p,g0,g0p,
     &              asymme,dasymm,d2asym,skappa,0,-2)
        procon = nprot/anumb
        call anucmat(rhonma,procon,eoveraa)
        call enernm(rhonm,half,pairfp,energm,efmasp,efmasn)
c        esurf3   = ((etotal+zpe)/anumb-eovera)*anumb3
        esurf2   = (etotal/anumb-eovera)*anumb3
        esurf3   = (etotal/anumb-eoveraa)*anumb3
        IF(ifplot>0) write(11,'(a,4(1pg12.4))')
     &      ' E/A(2),T,rhonm=',eovera,energm,pairfp,rhonm
c
c     print
c
        write(6,'(a,f7.4,a,f8.3,a,3f7.2)')
     &  ' N^-1/3=',one/anumb3,
     &  ' LDA-E/A=',elda/anumb,
     &  '  surface energies=',esurf1,esurf2
        IF(ifplot>0) write(11,'(a,f7.4,a,f8.3,a,3f7.2)')
     &  ' N^-1/3=',one/anumb3,
     &  ' LDA-E/A=',elda/anumb,
     &  '  surface energies=',esurf1,esurf2,esurf3
c        WRITE(12,'(2i7,4(1pg13.5))')
c     &    nprot,nneut,one/anumb3,1-2D0*procon,esurf2,esurf3
c
        call xkappa(skapp2)
      endif
c
      return
      end
c-----xkappa--------------------------------------------------------
c
      subroutine xkappa(skapp2)
c
      include 'params.inc'
c
c     Sum rule enhancement 'kappa' via effective masses
c       rhonm  = density at which 'kappa' is computed
c
c     Uses the following way to evaluate 'kappa':
c
c       hbar^2         hbar^2   hbar^2   1  d  hbar^2   hbar^2
c       ------ kappa = ------ - ------ + - -- (------ - ------)
c         2m            2m^*      2m     4 dy  2m_n^*   2m_p^*
c
c     where 'y' is a isospin-deformation defined as
c
c       rho_p -->  rho_nm (0.5 + y)   .
c
      include 'skypar.inc'
      include 'constn.inc'
c
      data deltpc/0.01/            ! variation of proton content
c
c-----------------------------------------------------------------------
c
c     equilibrium point
c
      call nucmat(rhonm,eovera,compr,effmas,f0,f0p,g0,g0p,
     &            asymme,dasymm,d2asym,skappa,0,-2)
c
c     + variation of proton content
c
      procon = 0.5+deltpc
      rhop   = procon*rhonm
      rhon   = rhonm-rhop
      amfcpp = amfc0(1)
     &          +t1p2*rhonm-t1m2*rhop
     &          +rhonm*(t5p6*rhonm-t5m6*rhop)
      amfcpn = amfc0(2)
     &          +t1p2*rhonm-t1m2*rhon
     &          +rhonm*(t5p6*rhonm-t5m6*rhon)
c
c     0 variation = average
c
      amfcav = 0.5*(amfcpp+amfcpn)
c
c     - variation of proton content
c
c      procon = 0.5-deltpc
c      rhop   = procon*rhonm
c      rhon   = rhonm-rhop
c      amfcmp = amfc0(1)
c     &          +t1p2*rhonm-t1m2*rhop
c     &          +rhonm*(t5p6*rhonm-t5m6*rhop)
c      amfcmn = amfc0(2)
c     &          +t1p2*rhonm-t1m2*rhon
c     &          +rhonm*(t5p6*rhonm-t5m6*rhon)
c
c     derivate and kappa
c
      amfc0a = 0.5*(amfc0(1)+amfc0(2))
      damfcp = (amfcpn-amfcpp)/(deltpc)
c      damfcn = (amfcpn-amfcmn)/(2.0*deltpc*rhonm)
      skapp2 = (amfcav-amfc0a+damfcp/4.0)/amfc0a
c
c     protocol
c
      IF(ifplot>0) write(11,'(a,2(1pg12.4))')
     &  ' XKAPPA: m*/m,skappa=',amfc0a/amfcav,skapp2
      IF(ifplot>0) write(11,'(a,2(1pg12.4))')
     &  ' NUCMAT: m*/m,skappa=',effmas,skappa
      return
      end
c-----hotnm --------------------------------------------------------
c
      subroutine hotnm(rhonmp,taunmp,temper)
c
      include 'params.inc'
c
c     solves Fermi equilibrium for hot nuclear matter:
c      rhonmp = input: total proton density
c      taunmp = kinetic energy density
c      temper = temperature
c
      include 'skypar.inc'
      include 'constn.inc'
c
      data f1d3  /0.333333333/
      data itmax/50/
      data delfp/0.01/
      data epsfp/1.0E-5/
      data ipgrid/200/
      data delp/0.01/
c
c-----------------------------------------------------------------------
c
      fermp  = (3.*pi*pi*rhonmp)**f1d3            ! first guess Fermi momentum
      amfca0 = 0.5*(amfc0(1)+amfc0(2))              ! average hbar^2/2m
      do iter=1,itmax
        sum0   = 0.0
        sum1   = 0.0
        fermp2 = fermp*fermp
        fermpp = (fermp+delfp)*(fermp+delfp)
        do ip=1,ipgrid
          pact2  = ip*delp
          pact2  = pact2*pact2
c          write(6,'(a,i5,4(1pg12.4))')
c     &    ' ip,pact,amfca0,frmp2,sum0=',ip,pact,amfca0,frmp2,sum0
          sum0   = pact2/(1.0+exp(amfca0*(pact2-fermp2)/temper)) + sum0
          sum1   = pact2/(1.0+exp(amfca0*(pact2-fermpp)/temper)) + sum1
        enddo
        sum0   = delp*sum0/(pi*pi)
        sum1   = delp*sum1/(pi*pi)
        drhodf = (sum1-sum0)/delfp
        fpstep = -(sum0-rhonmp)/sign(max(abs(drhodf),xsmall),drhodf)
        fpstep = sign(min(abs(fpstep),10.0*delfp),fpstep)
        fermp  = fermp+fpstep
c        write(6,'(a,i5,3(1pg12.4))')
c     &   ' HOTNM: iter,fermp,fpstep,rhonmp=',iter,fermp,fpstep,sum0
        if(abs(fpstep).lt.epsfp) goto 99
      enddo
c
c     Fermi momentum is determined. Now compute kinetic energy.
   99 continue
      sumk   = 0.0
      sum0   = 0.0
      fermp2 = fermp*fermp
      do ip=1,ipgrid
        pact2  = ip*delp
        pact2  = pact2*pact2
        vol    = pact2/(1.0+exp(amfca0*(pact2-fermp2)/temper))
        sum0   = vol       + sum0
        sumk   = vol*pact2 + sumk
      enddo
      sum0   = delp*sum0/(pi*pi)
      sumk   = delp*sumk/(pi*pi)
      taunmp = sumk
      averen = sumk/sum0*amfca0
      fermp  = (3.*pi*pi*rhonmp)**f1d3            ! first guess Fermi momentum
      avere0  = 0.6*fermp*fermp*amfca0
      write(6,'(a,5(1pg12.4))')
     &   ' HOTNM: rhonmp(1..2),taunmp,averen,avere0=',
     &    rhonmp,sum0,sumk,averen,avere0
      IF(ifplot>0) write(11,'(a,5(1pg12.4))')
     &   ' HOTNM: rhonmp(1..2),taunmp,averen,avere0=',
     &    rhonmp,sum0,sumk,averen,avere0
c
      return
      end

c-----enernm--------------------------------------------------------
c
      subroutine enernm(rhonm,procon,temper,energm,efmasp,efmasn)
c
      include 'params.inc'
c
c     E/A for nuclear matter
c      rhonm  = input: total density
c      procon = input: proton content (0.5 for symm.matter)
c      temper = temperature
c      energm = output: E/A in MeV
c      effmas = output: effective mass m*/m
c     The force paxrameters are entered via common /skypar/
c
      include 'skypar.inc'
      include 'wavfun.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
      data betab /5.742467999/     ! (3/5)*(3*pi^2)^*(2/3)  factor for tau_n/p
      data f2d3  /0.6666666666667D0/
      LOGICAL :: texchasav  ! to switch off Coulomb exchange in matter
c
c-----------------------------------------------------------------------
c
ctest      write(*,*) 'ENERNM entered'
      rhop   = procon*rhonm
      rhon   = rhonm-rhop
      if(temper.eq.zero) then
        taup   = betab*rhop**f2d3*rhop
        taun   = betab*rhon**f2d3*rhon
      elseif(procon.eq.half) then
        call hotnm(rhop,taup,temper)
        rhon   = rhop
        taun   = taup
      else
        stop ' non zero T only for symmetric nuclear matter '
      endif
!      IF(procon==half) WRITE(911,'(7(1pg16.7))')  rhop,taup,
!     &   rhop**f2d3,(amfc0(1)*taup+amfc0(2)*taun)/rhonm,
!     &   (amfc0(1)+amfc0(2))/2D0
      zpefac = one
      IF(ioptsk==13) THEN
!        texchasav = texcha
!        texcha = .FALSE.
        call fyfunc(rhop,rhon,zero,zero,zero,zero,zero,
     &              taup,taun,zero,zero,zero,zero,zero,zero,
     &              rhop,rhon,
     &              potp,potn,bmasp,bmasn,zero,zero,potlsp,potlsn,
     &              energd)
!        texcha = texchasav
!        IF(procon==half) THEN
!        de = 1D-6
!          call fyfunc(rhop+de,rhon-de,zero,zero,zero,zero,zero,
!     &              taup,taun,zero,zero,zero,zero,zero,zero,
!     &              rhop,rhon,
!     &              potp,potn,bmasp,bmasn,zero,zero,potlsp,potlsn,
!     &              energdpp)
!          call fyfunc(rhop-de,rhon+de,zero,zero,zero,zero,zero,
!     &              taup,taun,zero,zero,zero,zero,zero,zero,
!     &              rhop,rhon,
!     &              potp,potn,bmasp,bmasn,zero,zero,potlsp,potlsn,
!     &              energdpm)
!          WRITE(911,*) ' energs:',rhonm,energd/rhonm,
!     &            rhonm*(energdpp-2D0*energd+energdpm)/(8D0*de**2)
!        END IF
      ELSE 
        call skfunc(rhop,rhon,zero,zero,zero,zero,zero,
     &              taup,taun,zero,zero,zero,zero,zero,zero,
     &              rhop,rhon,
     &              potp,potn,bmasp,bmasn,bmaslp,bmasln,potlsp,potlsn,
     &              energd)
      END IF
ctest      write(*,*) ' ENERNM after skfunk'
      enrg   = zpefac*(amfc0(1)*taup+amfc0(2)*taun) + energd
c        if(iwr.gt.-2)
c     &   write(6,'(a,g12.4)') ' ekin=',(amfc0(1)*taup+amfc0(2)*taun)/rhonm
      energm = enrg/rhonm
c      write(2,'(a,3g12.4)') 
c     &   ' rhonm,procon,energm=',rhonm,procon,energm
      efmasp = zpefac*amfc0(1)/(zpefac*amfc0(1)+bmasp)
      efmasn = zpefac*amfc0(2)/(zpefac*amfc0(1)+bmasn)
      IF(iwr>-1) 
     & WRITE(6,'(a,7(1pg13.5))') 'ENERNM:',rhop,rhon,taup,taun,
     &   zpefac*(amfc0(1)*taup+amfc0(2)*taun)/rhonm,energd/rhonm,energm
!      WRITE(6,'(7x,6(1pg13.5))') 'ENERNM:',rhop,rhon,taup,taun,amfc0
      
c
      return
      end

c-----ekinnm--------------------------------------------------------
c
      function ekinnm(rhonm,procon)
c
      include 'params.inc'
c
c     E_kin/A for nuclear matter
c      rhonm  = input: total density
c      procon = input: proton content (0.5 for symm.matter)
c      ekinnm = kinetic energy per particle
c     The force paxrameters are entered via common /skypar/
c
      include 'skypar.inc'
      include 'wavfun.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
      data betab /5.742467999/     ! (3/5)*(3*pi^2)^*(2/3)  factor for tau_n/p
      data f2d3  /0.666666667/
c
c-----------------------------------------------------------------------
c
      rhop   = procon*rhonm
      rhon   = rhonm-rhop
      taup   = betab*rhop**f2d3*rhop
      taun   = betab*rhon**f2d3*rhon
      ekinnm = (amfc0(1)*taup+amfc0(2)*taun)/rhonm
c
      return
      end

c-----densit------------ part of the hartree-fock package ----------------------
c
      subroutine densit(occwgt,v2wgt)
c
      include 'params.inc'
c
c     * * * * * * * * *      densities and currents      * * * * * * * * * *
c
c     accumulates the density 'rho',
c                 the kinetic density 'tau' and
c                 the spin-orbit density 'sj'
c     from the radial single particle wave functions given on 'wf'
c     and their derivatives given on 'wfs'.
c
c     'occwgt'  is the actual occupation number of the s.p. states
c               entered by the calling list

c
c     each density occurs as a proton- and a neutron-density.
c     the proton density ('iso'=1) is stored first in the
c     interval [1,...,'ngrid']; directly after that the neutron
c     density is stored in the interval [ngrid+1,...,2*ngrid].
c
c     for this task the routine needs some information on the shell
c     structure; this is provided by the following fields:
c      jsp   = total and orbital angular momentum for the state 'nshl'.
c              (for a more detailed explanation see 'initgr')
c      zwj   = "nr. of states in shell with combined ang.mom. 'js'",
c              i.e. (2*j+1)/4*pi where j is total angular momentum.
c      wls   = strength of spin-orbit splitting for ang.mom. 'js'
c     (note that 'js' and also 'jsp' is a combined quantum number
c      for total and orbital ang.mom.)
c
      include 'skypar.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'constn.inc'
      dimension occwgt(ksl2),v2wgt(ksl2)
c
c      dimension zwj(21)
c
c-------------------------------------------------------------------------------
c
c     zero fields before accumulation
c
      do i=1,kgr2
        rho(i) = zero
        tau(i) = zero
        taul(i) = zero
        sj(i)  = zero
        paircr(i) = zero
      enddo
c
c     iso-loop, initialise "shift constants" 'im,ia,ie' for the
c     densities, 'ngw' for the weights, and 'nadd' for the
c     wave functions.
c
      do 20 iso=1,2
      ie     = ngrid*iso
      im     = ie-ngrid+1
      ia     = im+1
      if(iso.eq.1) then
       nmax   = nmaxp
       nblock = nblockp
      else
       nmax   = nmaxn
       nblock = nblockn-nmaxp
      endif
      ngw    = (iso-1)*nmaxp
      nadd   = ngw*ngrid-ie
c
c     loop over states for given isospin.
c      'nshl' counts the states in harmonic oscillator ordering.
c      'js' is the according composed ang.mom.,
c      'll' is the orbital ang.mom.,
c      'wl' is the spin-orbit splitting,
c      'om' is the actual phase-space factor
c
      do 30 nshl=1,nmax
      js     = jsp(nshl)
      ll     = js/2
      xlp    = float(ll*ll+ll+1)
      wl     = wls(js)
      IF(tphpair) THEN
        om     = zwj(js)*occwgt(nshl+ngw) 
      ELSE
        om     = zwj(js)*occwgt(nshl+ngw)*phpair(nshl+ngw) 
      END IF
      v2     = v2wgt(nshl+ngw)
      uv = phpair(nshl+ngw)*sqrt(max(v2-v2*v2,xsmall))
      nadd   = nadd+ngrid
      IF(nshl==nblock) THEN
!        uv = sqrt(max(v2b-v2b*v2b,xsmall))
        ompair = zwj(js)*uv*0.5D0*(ph(nblock+ngw)-1D0)/ph(nblock+ngw)
      ELSE
!        uv     = sqrt(max(v2-v2*v2,xsmall))
        ompair = zwj(js)*uv*0.5D0              ! factor 1/2 because "state>0"
      END IF
      if(om.lt.xsmall) goto 30
c                                  the case r=0.0 needs special treatment
      if(ll.lt.1) then
       x      = wfs(im+nadd)
       rho(im) = (x*x)*om + rho(im)
       IF(.NOT.thfbfold) paircr(im) = (x*x)*ompair + paircr(im)
      elseif(ll.eq.1) then
       x      = (16.0*wf(ia+nadd)-wf(ia+nadd+1))/(6.0*rstep*rstep)
       tau(im) =.75*om*(x*x) + tau(im)
      endif
c                                  now all r>0.0
      r      = zero
      do i=ia,ie
        r      = r+rstep
        ri     = one/r
        phas   = (ri*ri)*om
c
        wfa    = wf(nadd+i)
        wfa2   = wfa*wfa
c       if(wfa2.gt.xsmall) then
         rho(i) = (phas*wfa2) + rho(i)
         taul(i) = (phas*wfa2)*xlp + taul(i)
         sj(i)  = (phas*wfa2)*wl*ri + sj(i)
c
         wfas   = wfs(nadd+i)
         tau(i) = (wfas*wfas - ((wfa*wfas)+(wfa*wfas)-xlp*wfa2*ri)*ri)
     &           *phas + tau(i)
         IF(thfbfold) wfa2=wffold(nadd+i)**2
         paircr(i) = ompair*(ri*ri)*wfa2 + paircr(i)
c     endif
      enddo
c
   30 continue
c
c     some derivatives: d/dr 'rho' --> 'd1r'
c                       d/dr**2 'rho' + 2/r d/dr 'rho' --> 'd2r'
c                       d/dr 'sj' + 2/r 'sj' --> 'sjs'
c     'im' is the initial position for given 'iso'.
c
      zwh    = 2.0/rstep
      call derv15(rho(im),d1r(im), one,ngrid,rstep)
      call derv25(rho(im),d2r(im), one,ngrid,rstep)
      d2r(im) = d2r(im) + (8.0*d1r(im+1)-d1r(im+2))/(3.0*rstep)
      do i=im+1,ie
        d2r(i) = d2r(i) + zwh*d1r(i)/(i-im)
      enddo
      if(t9.ne.zero) then
       call derv15(tau(im),d1t(im), one,ngrid,rstep)
       call derv25(tau(im),d2t(im), one,ngrid,rstep)
       d2t(im) = d2t(im) + (8.0*d1t(im+1)-d1t(im+2))/(3.0*rstep)
       do i=im+1,ie
         d2t(i) = d2t(i) + zwh*d1t(i)/(i-im)
       enddo
      endif
c
      call derv15(sj(im),sjs(im),-one,ngrid,rstep)
      sjs(im) = 3.0*sjs(im)
      do i=im+1,ie
        sjs(i) = sjs(i) + zwh*sj(i)/(i-im)
      enddo

!    extrapolate to r=0.0
      IF(thfbfold) CALL extrap(paircr(im),1D0)
c
c     end iso-loop
c
   20 continue
c
      return
      end
c-----prispd------------ part of the hartree-fock package ----------------------
c
      subroutine prispd
c
      include 'params.inc'
c
c     computing and printing the single particle density on unit 'for028'.
c     the densities are printed in columns to allow easy input to standard
c     plotting software.
c
c     the printed information is transferred via common's.
c
      include 'wavfun.inc'
      include 'result.inc'
      include 'grid.inc'
c
      logical top28
      character*5 nucleo(2)
      data nucleo/' prot','neutr'/
c
c-------------------------------------------------------------------------------
c
      inquire(28,opened=top28)
      if(.not.top28)
     &  open(28,file='spdens.res')
c     &  open(28,file='for028')
      write(28,'(/a/a,2i4/a,1pg13.5)')
     & '#**********************************************************',
     & '# Single nucleon densities and w.f.s for the nucleus Z,N=',  
     &   nprot,nneut,'# having E_B=',etotal
c
c     loop over isospin
c
      ntotal=0
      phasej = zwj(1)*0.5
      do iso=1,2
        isadd  = ngrid*nmaxp*(iso-1)
        if(iso.eq.1) then
         nend   = nmaxp
        isoadd = 0
        else
         nend   = nmaxn
         isoadd = nmaxp
        endif
c
c       loop over states in oscillator ordering
c
        do n2=1,nend
          ntotal=1+ntotal
          n      = nsort(n2+isoadd)-isoadd
          js     = jsp(n)
          ll     = js/2
          j2     = ((js+1)/2)*2-1
          ie     = n*ngrid+isadd
          ia     = ie-ngrid+1
          occfac = weight(isoadd+n)*zwj(js)/phasej
          if(ll.eq.0) then
            rhosp  = wfs(ia)**2*phasej
          else
            rhosp  = 0.0
          endif
c
        write(28,'(/1x,a,a)') 
c          write(28,'(3a,i3,2(a,i2),a,f8.3,a,f5.3/a)')
c     &     ' state: ',nucleo(iso),'ons, n=',n,',j=',j2,'/2,l=',ll,
c     &     ',  s.p.en.=',spenrg(isoadd+n),
c     &     ',  occup.=',weight(isoadd+n),
c     &     '  r      rho_state  r*wavefunction '
          write(28,'(3a/a,3i5,f5.1,i5,2g13.5,i5/a)')
     &     '# state: ',nucleo(iso),
     &     'ons, info in order: nr,iso,n,j,l,spenerg,occup,ngrid=','#',
     &     ntotal,
     &     iso,n,j2*0.5d0,ll,spenrg(isoadd+n),weight(isoadd+n),ngrid,
     &     '#  r      rho_state   r*wavefunction: '
c     &     '  r   r**2   rho_state   r**2*rho_state   rho_stat*occup: '
c          write(28,'(1x,f5.2,f7.2,3(1pg13.5))') 
c     &      0.0,0.0,rhosp,0.0,rhosp*occfac
          write(28,'(1x,f5.2,3(1pg13.5))')    0.0,rhosp,0.0
          do i=ia+1,ie
            r     = (i-ia)*rstep
            rhosp = (wf(i)/r)**2*phasej
c            write(28,'(1x,f5.2,f7.2,3(1pg13.5))') 
c     &        r,(r*r),rhosp,(r*r)*rhosp,rhosp*occfac
            write(28,'(1x,f5.2,1pg13.5,1pg20.12)')  r,rhosp,wf(i)
          enddo
c
        enddo
c
      enddo
c
      return
      end
c-----skypot------------ part of the hartree-fock package ----------------------
c
      subroutine skypot(iter)
c
      include 'params.inc'
c
c     * * * * * * *     skyrme-hartree-fock potentials    * * * * * * * *
c
c     composes single particle potential on 'pot',
c              spin-orbit potential on 'potls', and
c              inverse mass on 'bmass'
c     using the various densities accumulated in subroutine 'densit'.
c     proton potentials are stored first, in the interval [1,...,'ngrid'];
c     neutron potentials follow in the interval ['ngrid'+1,...,2*'ngrid'].
c
c     the emerging potentials are mixed with the previous potentials
c     with a mixing factor 'addnew'; this should stabilise the iteration.
c
      logical tpospr,tposne
cold      logical tcount
c
      include 'skypar.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'result.inc'
      include 'potinv.inc'
      include 'constn.inc'
c
      dimension worki(kg1)
      REAL(8) :: workold(kgr2),worknew(kgr2)
      logical ttest
      logical testen
c
c     some constants:
c      spiba    = coulomb contribution to spin-orbit potential
c      iflinw   = switch for artificial wall on the local potential
c                 ( 1 --> linear extrapolation wall,
c                   0 --> square-well .             )
c
      data  spiba/0.0308/
      data iflinw/1/
      data itste/0/
      data ttest/.false./
      data testen/.false./

      logical TPAIRON
      parameter (TPAIRON=.true.)
c
c    reduced Fermi-momentum for exchange functional as internal statement fct.
c      (29.6088132=3*pi**2)
c
!      eta(rhoact) = (29.6088132*rhoact)**0.333333333d0*2.0/exmass
c
c------------------------------------------------------------------------------
c
c     setting addition weights:
c      for 'iter'=1 the potentials are computed from scratch,
c      in all other cases they are mixed between new and old
c      with the weight 'addact' for the new contribution.
c
      if(iter.eq.1) then
!      if(iter.eq.1 .OR. iter==900) then
       addact=one
      else
       addact=addnew
      endif
      addold=one-addact
c
c     folding the densities with yukawa of length ayuk using the
c     inversion 'invert'.
c
c                                            save charge density before folding
c      if(ifdfms.ne.1 .and. ifdfms.ne.3) then
      if(ifrhoc.eq.0) then
        do ir=1,ngrid
          rhochr(ir) = rho(ir)
        enddo
      else
       call gaufol(rho,rhochr,one*0.8)
      endif
c
      if(tyukaw) then
       do i=1,ngrid
         d1t(i)       = rho(i)                       ! use 'd1t' as workspace
         d1t(i+ngrid) = rho(i+ngrid)
       enddo
       call lapfil(d1t,worki,ayms2,-1)
       call lapfil(d1t(ngridp),worki,ayms2,-1)
       if(addact.eq.one .and. iwr.ge.-1) then
        write(06,'(a,20(/1x,10g13.5))')
     &        'testprint yukawa-potential:',(d1t(i),i=1,ngrid)
        write(06,'(20(/1x,10g13.5))')
     &        (d1t(i),i=ngridp,ngrid2)
       endif
c
c      folding of all densities for thourough finite range potentials
      elseif(ayuk.ne.zero) then
        call nucmat(rhosnm,eovera,compr,effmas,f0,f0p,g0,g0p,
     &              asymme,dasymm,d2asym,skappa,-1,-2)
ctoinv       dampay = ayuk*ayuk
        write(6,'(a,3g12.4)')
     &   ' SKYPOT: delq,rhosnm,ayuk=',delq,rhosnm,ayuk
       if(ttest) then
         do ir=1,ngrid
           aux1(ir) = rho(ir)
           aux2(ir) = rho(ir)
         enddo
         call lapfil(aux1,worki,ayms2,-1)
         call t0inv(aux2,aux2,worki,dampay,0)
         ttest = .false.
         write(06,'(a/a,3g13.5/a,500(/1x,i3,3g13.5))')
     &        'testprint folding:',
     &        '  ayuk,ayms2,dampay=',ayuk,ayms2,dampay,
     &        '  i    rho    lapfil    t0inv  ',
     &        (i,rho(i),aux1(i),aux2(i),i=1,ngrid)
       endif
       do iso=1,2
         im = ngrid*(iso-1)+1
         ie = im+ngrid-1
c         call lapfil(rho(im),worki,ayms2,-1)
c         call lapfil(tau(im),worki,ayms2,-1)
c         call lapfil(d2r(im),worki,ayms2,-1)
         if(ayukp.ne.zero) then
ctoinv          call t0inv(rho(im),rho(im),worki,ayukp*ayukp,0)
ctoinv          call t0inv(d1r(im),d1r(im),worki,ayukp*ayukp,1)
ctoinv          call t0inv(d2r(im),d2r(im),worki,ayukp*ayukp,0)
          if(iter.eq.1 .and. im.eq.1) then
            write(6,'(a,30(/1x,5g12.4))') 
     &           ' rho before:',(rho(ir),ir=im,ie)
          endif
          call filter0(rho(im),rho(im))
          if(iter.eq.1 .and. im.eq.1) then
            write(6,'(a,30(/1x,5g12.4))') 
     &           ' rho after:',(rho(ir),ir=im,ie)
          endif
          call filter0(d1r(im),d1r(im))
          call filter0(d2r(im),d2r(im))
         endif
         if(ayukm.ne.zero)
     &    call filter0(tau(im),tau(im))
ctoinv     &    call t0inv(tau(im),tau(im),worki,ayukm*ayukm,0)
         if(tssmas) call lapfil(d2t(im),worki,ayms2,-1)
c         call lapfil(sjs(im),worki,ayms2,-1)  
         if(ayukls.ne.zero) then
           if(iter.eq.1 .and. im.eq.1) then
             write(6,'(a,30(/1x,5g12.4))') 
     &           ' sjs before:',(sjs(ir),ir=im,ie)
           endif
           call filter0(sjs(im),sjs(im))
           if(iter.eq.1 .and. im.eq.1) then
             write(6,'(a,30(/1x,5g12.4))') 
     &           ' sjs after:',(sjs(ir),ir=im,ie)
           endif
           call filter0(sj(im),sj(im))
          write(6,'(A,i2)') ' l*s branch. IM=',im
ctoinv           call t0inv(sjs(im),sjs(im),worki,ayukls*ayukls,0)
ctoinv           call t0inv(sj(im),sj(im),worki,ayukls*ayukls,1)
         endif
         if(imfold.eq.1) then
           write(6,'(A,i2)') ' IMFOLD branch'
           call lapfil(tau(im),worki,ayms2,-1)
           call lapfil(d2r(im),worki,ayms2,-1)
           if(tssmas) then
            call lapfil(d2t(im),worki,ayms2,-1)
            call lapfil(d2t(im),worki,ayms2,-1)
           endif
           call lapfil(sjs(im),worki,ayms2,-1)
           call t0inv(sj(im),sj(im),worki,dampay,1)
           call t0inv(d1r(im),d1r(im),worki,dampay,1)
         endif
       enddo
       do ir=1,ngrid
         rhochr(ir) = rho(ir)         ! want possibly folded charge density
       enddo
      endif
c
c      compose the potential 'pot', l*s-potential 'potls' and
c      and inverse mass 'bmass' from densities.
c      at the same time the rearrangment energy is accumulated
c      on 'enrear'.
c
c                                  save old potentials
      do ir=1,ngrid+ngrid
        aux1(ir) = pot(ir)
        aux2(ir) = bmass(ir)
        aux3(ir) = potls(ir)
      enddo
c                                  prepare some auxiliary variables
      if(tpair .and. TPAIRON) parfac = pairgm/rho0pr         ! factor on pairing potential
c                                  iso-loop,
c                                  'im' is shifter to conjugate isospin
      enrear = zero
      epotsk = zero               ! accumulator for Skyrme potential E
      x      = zero
!      WRITE(6,*) 'testprint rho,pot:'
!      WRITE(6,*) 'rho:',rho(1:2*ngrid)
      do i=2,ngrid
        x      = x+rstep
        j      = i+ngrid
        IF(ioptsk==13) THEN
          call fyfunc(rho(i),rho(j),rhochr(i),d1r(i),d1r(j),
     &                d2r(i),d2r(j),
     &                tau(i),tau(j),taul(i),taul(j),sj(i),sj(j),
     &                sjs(i),sjs(j),d1t(i),d1t(j),
     &                pot(i),pot(j),bmass(i),bmass(j),
     &                bmassl(i),bmassl(j),potls(i),potls(j),energd)
!          WRITE(6,'(i3,f5.1,4(1pg13.5))')i,x,rho(i),rho(j),pot(i),pot(j)
        ELSE
          call skfunc(rho(i),rho(j),rhochr(i),d1r(i),d1r(j),
     &                d2r(i),d2r(j),
     &                tau(i),tau(j),taul(i),taul(j),sj(i),sj(j),
     &                sjs(i),sjs(j),d1t(i),d1t(j),
     &                pot(i),pot(j),bmass(i),bmass(j),
     &                bmassl(i),bmassl(j),potls(i),potls(j),energd)
        END IF
        if(itste.eq.1) 
     &     write(6,'(a,i4,3g12.4)') 'i,E=',i,energd,rho(i),rho(j)
        epotsk = (x*x)*energd  + epotsk
!        IF(i==3) WRITE(6,'(a,i3,5(1pg13.5))') 
!     &    'testpot:',i,energd,rho(i),rho(j),pot(i),pot(j)      
c
c       contribution from density-dependent pairing
c
!        tpair=.TRUE.
!        WRITE(6,*) 'SKYPOT: tpair=',tpair
        if(tpair .and. TPAIRON) then
!          WRITE(6,*) 'tpair switch'
          IF(.not.tpairFy) THEN
            if(pairgm.eq.1D0) then
              adparp = (pairfp*paircr(i)*paircr(i)
     &                  +pairfn*paircr(j)*paircr(j))*parfac
            else
              adparp = (pairfp*paircr(i)*paircr(i)
     &                +pairfn*paircr(j)*paircr(j))
     &              *parfac*((rho(i)+rho(j))/rho0pr)**(pairgm-one)
            endif
            IF(tregpair) THEN
              pot(i) = adparp*regpairadd(i) + pot(i)
              pot(j) = adparp*regpairadd(j) + pot(j)
            ELSE
              pot(i) = adparp + pot(i)
              pot(j) = adparp + pot(j)
            END IF
!            ppadd(i) = adparp
!            ppadd(j) = paircr(i)**2+paircr(j)**2
!            enrear = (x*x)*(rs)*adparp  + enrear    ! ???
!            enrear = (x*x)*(rs)*adparp  + enrear    ! ???
            enrear = (x*x)*(rho(i)*regpairadd(i)+
     &                      rho(j)*regpairadd(j))*adparp  + enrear   
          END IF
        endif
c
c       add external Woods-Saxon potential in case of neutron drop
c
        if(weight(1).eq.zero .and. nneut.le.8) then
          potext = -20.0/(one+exp((x-3.0)/0.65))   ! neutron stabilization
          pot(j) = potext + pot(j)
c          enrear = -(x*x)*rn*potext  + enrear
c???          epotsk = (x*x)*rn*potext  + epotsk
        endif
c
c       1/r factor on spin-orbit potential
c
        potls(i) = potls(i)/x
        potls(j) = potls(j)/x
      enddo
      enrear = -frpirs*half*enrear
      epotsk = frpirs*epotsk
      itste = 0
!
!     Fayans functional: add gradient contribution to potential 
!                        and pairing contribution to potential
!
      IF(tpairFy) CAll fy_pairpot()
      IF(ioptsk==13) THEN
!        CALL fyfunc_grad_linres(worknew,1D0,1D0,1D0,1D0,1D0,1D0)
!        workold = pot
        CALL fyfunc_grad()
!        workold = pot-workold
!        OPEN(111,file='testgrad')
!        WRITE(111,'(200(i4,3(1pg13.5)/))')
!     &   (i,worknew(i),workold(i),worknew(i)-workold(i),i=1,ngrid+ngrid)
!        CLOSE(111)
!        pot = pot + worknew
      END IF
c
      if(.not.tyukaw .and. ayuk.ne.zero) then
c       call lapfil(pot(1),worki,ayms2,-1)
c       call lapfil(pot(ngrid+1),worki,ayms2,-1)
       if(ayukp.ne.zero) then
         write(6,'(A,i2)') ' potential branch'
           if(iter.eq.1) then
             write(6,'(a,30(/1x,5g12.4))') 
     &           ' pot before:',(pot(ir),ir=1,ngrid)
           endif
         call filter0(pot(1),pot(1))
           if(iter.eq.1) then
             write(6,'(a,30(/1x,5g12.4))') 
     &           ' pot after:',(pot(ir),ir=1,ngrid)
           endif
         call filter0(pot(ngrid+1),pot(ngrid+1))
ctoinv         call t0inv(pot(1),pot(1),worki,ayukp*ayukp,0)
ctoinv         call t0inv(pot(ngrid+1),pot(ngrid+1),worki,ayukp*ayukp,0)
       endif
       if(ayukls.ne.zero) then
         write(6,'(A,i2)') ' l*s branch'
         call filter0(potls(1),potls(1))
         call filter0(potls(ngrid+1),potls(ngrid+1))
ctoinv         call t0inv(potls(1),potls(1),worki,ayukls*ayukls,1)
ctoinv         call t0inv(potls(ngrid+1),potls(ngrid+1),worki,ayukls*ayukls,1)
       endif
c       call lapfil(bmass(1),worki,ayms2,-1)
c       call lapfil(bmass(ngrid+1),worki,ayms2,-1)
       if(ayukm.ne.zero) then
         write(6,'(A,i2)') ' mass branch'
           if(iter.eq.1) then
             write(6,'(a,30(/1x,5g12.4))') 
     &           ' bmass before:',(bmass(ir),ir=1,ngrid)
           endif
         call filter0(bmass(1),bmass(1))
           if(iter.eq.1) then
             write(6,'(a,30(/1x,5g12.4))') 
     &           ' bmass after:',(bmass(ir),ir=1,ngrid)
           endif
         call filter0(bmass(ngrid+1),bmass(ngrid+1))
ctoinv         call t0inv(bmass(1),bmass(1),worki,ayukm*ayukm,0)
ctoinv         call t0inv(bmass(ngrid+1),bmass(ngrid+1),worki,ayukm*ayukm,0)
       endif
       if(imfold.eq.1) then
        call lapfil(bmass(1),worki,ayms2,-1)
        call lapfil(bmass(ngrid+1),worki,ayms2,-1)
       endif
      else
       pot(1)=(4.0*pot(2)-pot(3))*0.33333333
       potls(1)=(4.0*potls(2)-potls(3))*0.33333333
       bmass(1)=(4.0*bmass(2)-bmass(3))*0.33333333
       pot(1+ngrid)=(4.0*pot(2+ngrid)-pot(3+ngrid))*0.33333333
       potls(1+ngrid)=(4.0*potls(2+ngrid)-potls(3+ngrid))*0.33333333
       bmass(1+ngrid)=(4.0*bmass(2+ngrid)-bmass(3+ngrid))*0.33333333
      endif
c                                         mix old and new potentials
      do ir=1,ngrid
        pot(ir)   = addold*aux1(ir) + addact*pot(ir)
        bmass(ir) = addold*aux2(ir) + addact*(bmass(ir)+amfc(1))
        potls(ir) = addold*aux3(ir) + addact*potls(ir)
        ic        = ir+ngrid
        pot(ic)   = addold*aux1(ic) + addact*pot(ic)
        bmass(ic) = addold*aux2(ic) + addact*(bmass(ic)+amfc(2))
        potls(ic) = addold*aux3(ic) + addact*potls(ic)
      enddo
c
c      coulomb-energy
c       the fields 'sqb' and 'sqm' are used as work-space inside 'coulmb'.
c       the field 'binv' is used as workspace for carrying the coulomb
c       potential
c
      if(acoul.ne.zero) then
        call coulmb(rhochr(1),binv,acoul)
        IF(TADSIC) THEN
          binv(1:ngrid) = binv(1:ngrid)*(1D0-1D0/nprot)
        END IF
       ecoul = zero
       do i=1,ngrid
         pot(i) = binv(i)*addact + pot(i)
         ecoul  = rsquar(i)*rhochr(i)*binv(i) + ecoul
       enddo
       ecoul  = half*frpirs*ecoul
       epotsk = ecoul + epotsk
!       WRITE(*,*) 'ecoul=',ecoul
c
c      el.-magn. spin-orbit potential for protons
c
       if(ifex.gt.1) then
        sqm(1) = (8.0*sqb(2)-sqb(3))/(6.0*rstep)
        do i=2,ngrid
          sqm(i) = sqb(i)/((i-1)*rstep)
        enddo
        vor    = half*addact*spiba/(rstep*rstep)
        potls(1) = 4.0*vor*(sqm(2)-sqm(1)) + potls(1)
        addenr = zero
        do i=2,ngrid-1
          potls(i) = vor*(sqm(i+1)-sqm(i-1))/(i-1) + potls(i)
          addenr = rsquar(i)*sj(i)*potls(i) + addenr             ! ????
        enddo
        potls(ngrid) = 2.0*vor*(sqm(ngrid)-sqm(ngrid-1))/(ngrid-1)
     &                + potls(ngrid)
        addenr = frpirs*addenr
        epotsk = addenr+epotsk
       endif
c
c       for ifrhoc.eq.1 a correction to the coulomb-energy
c       is calculated coming form the finite extension of
c       proton.
c       parameters : 4*pi*e**2 = 18.09557
c                    12/radprot**2 = 18.75
c       the field 'yuk' is used as workspace for the laplacian
c       inversion. the work-fields 'potl', 'sqb' and 'd2b'
c       are preoccupied as to describe the folding equation
c         ( -laplacian + 18.75 ) 'yuk' = 'rho' .
c
cold       if(ifrhoc.eq.1) then
cold        umin   = zero
cold        do 270 i=1,ngrid
cold        aux3(i)  = (i-1)*rstep*(-18.09557*rho(i))
cold        sqb(i)  = one
cold        d2b(i)  = zero
cold        potl(i) = 18.75
cold  270   continue
cold        call invert(aux3(1),aux3(1),-1)
coldc
cold        aux3(1) = (8.0*aux3(2)-aux3(3))/(6.0*rstep)
cold        pot(1) = aux3(1)*addact + pot(1)
cold        do 280 i=2,ngrid
cold        aux3(i) = aux3(i)/((i-1)*rstep)
cold        pot(i) = aux3(i)*addact + pot(i)
cold  280   continue
c
cold        vor    = one/(18.09557*rstep*rstep)
cold        eccorr = zero
cold        do 290 i=2,ngrid-1
cold        d2ba   = (i*aux3(i+1)-2.0*(i-1)*aux3(i)+(i-2)*aux3(i-1))
cold     &          *vor/(i-1)
cold 290    eccorr = ((i-1)*(i-1))*aux3(i)*(rho(i)+d2ba) + eccorr
cold        rst2 = rstep*rstep          !   ???? this is uncertain !
cold        eccorr = rst2*frpirs*eccorr
cold        enrear = eccorr+enrear
coldc
cold        if(iwr.ge.-1 .and. iter.eq.1) then
cold         write(06,'(a,g12.4)')' correction to coulomb-correction=',
cold     &                        eccorr
cold         write(06,'(1h0,a,20(/1x,10g12.4))')
cold     &        'testprint coul.cor-potential:', (aux3(i),i=1,ngrid)
cold        endif
cold       endif
c
      endif
c
c     build an artificial wall to stabilise states with
c     positive energy.
c
c
c     test occurence of positive energy states
c
      if(wall.eq.zero) then
       spmxpr = -999.0
       do nshl=1,nmaxp
         if(spenrg(nshl).gt.spmxpr) spmxpr = spenrg(nshl)
       enddo
       tpospr = (spmxpr.gt.zero)
       spmxne = -999.0
       do nshl=nmaxp+1,nmaxto
         if(spenrg(nshl).gt.spmxne) spmxne = spenrg(nshl)
       enddo
       tposne = (spmxne.gt.zero)
      elseif(wall.gt.zero) then
       tpospr = .true.
       tposne = .true.
       spmxne = wall-3.0
       spmxpr = wall-3.0
      else
       tpospr = .false.
       tposne = .false.
      endif
c
c     determine position of maximum of the proton potential
c
      if(tpospr .or. tposne) then
cold       potmax = zero
cold       tcount = .false.
       do ir=2,ngrid
         if(pot(ir).gt.zero) then
           irpmax = ir
           goto 339
         endif
       enddo
  339  continue
       if(irpmax.le.2) irpmax=ngrid/2
       if(iwr.ge.-1)
     &     write(6,'(a,i5,2(a,g12.4))')
     &     ' irpmax=',irpmax,' spmxpr=',spmxpr,' spmxne=',spmxne
      endif
c
c     add a wall to 'pot'
c
      spmxpr = spmxpr+3.0
      if(tpospr) then
       if(iflinw.eq.1) then
       if(iter.eq.1 .and. iwr.ge.-1)
c       if( iwr.ge.-1)
     &     write(6,'(a)') ' proton potential with linear wall'
        do ir=irpmax+1,ngrid
          potadd = min((ir-irpmax)*10.0*rstep,spmxpr)
          pot(ir) = potadd
        enddo
       else
        if(iter.eq.1 .and. iwr.ge.-1)
     &     write(6,'(a)') ' proton potential with square wall'
        do ir=irpmax+1,ngrid
          pot(ir) = spmxpr
        enddo
       endif
      endif
c
      spmxne = spmxne+3.0
      if(tposne) then
       if(iflinw.eq.1) then
        if(iter.eq.1 .and. iwr.ge.-1)
     &     write(6,'(a)') ' neutron potential with linear wall'
        do ir=irpmax+1,ngrid
          potadd = min((ir-irpmax-1)*10.0*rstep,spmxne)
          pot(ir+ngrid) = potadd
        enddo
       else
        if(iter.eq.1 .and. iwr.ge.-1)
     &     write(6,'(a)') ' neutron potential with square wall'
        do ir=irpmax+1,ngrid
          pot(ir+ngrid) = spmxne
        enddo
       endif
      endif
c
c     optionally test separate energies
c
      if(testen) then
        sumt0 = 0.0
        sumt1 = 0.0
        sumls = 0.0
        sumkin = 0.0
        do i=2,ngrid
          j      = i+ngrid
          r      = (i-1)*rstep
          vol    = r*r
          rp     = rho(i)
          rm     = rho(j)
          rs     = rp+rm
          tp     = tau(i)
          tm     = tau(j)
          ts     = tp+tm
          ssp    = sjs(i)
          ssm    = sjs(j)
          sss    = ssp+ssm
          sp     = sj(i)
          sm     = sj(j)
          ss     = sp+sm
          d1p    = d1r(i)
          d1m    = d1r(j)
          d1s    = d1p+d1m
          d2p    = d2r(i)
          d2m    = d2r(j)
          d2s    = d2p+d2m
          sumt0  = vol*(t0a*rs*rs-t0b*(rp*rp+rm*rm))     + sumt0  
          sumt1  = vol*(2.0*(t1p2*ts*rs-t1m2*(tp*rp+tm*rm))
     &                  -t13m2*d2s*rs+t13p2*(d2p*rp+d2m*rm))
     &            + sumt1
          sumls  = vol*( -b4*sss*rs- b4p*(ssp*rp+ssm*rm)) + sumls
          sumkin = vol*(amfc(1)*tp+amfc(2)*tm)    + sumkin
        enddo
        sumt0  = frpirs*sumt0*0.5
        sumt1  = frpirs*sumt1*0.5
        sumls  = frpirs*sumls
        sumkin = frpirs*sumkin
!        write(6,'(a,4g13.5)')
!     &    ' epot_0,epot_1,epot_ls,ekin=',sumt0,sumt1,sumls,sumkin
      endif
c
      if(constr0.NE.0D0) then
        do ir=1,ngrid
          addcnstr = addact * constr0*((ir-1)*rstep)**2
          pot(ir)   = pot(ir) + addcnstr
          pot(ir+ngrid)   = pot(ir+ngrid) + addcnstr
        enddo
      endif
c
      return
      end
c-----coulmb------------ part of the hartree-fock package ----------------------
c
      subroutine coulmb(crho,coulpt,chpart)
c
      include 'params.inc'
c
c     * * * * * * * * * *      coulomb potential         * * * * * * * * *
c
c     evaluates the coulomb potential for given charge density.
c     the parameters in the calling list mean:
c      crho    = density of particles from which the coulomb
c                potential is to be computed ( in [fm**(-3)] ).
c      coulpt  = emerging coulomb potnetial ( in [mev]).
c      chpart  = fractional charge of the particles which constitute
c                'crho'. e.g. for protons set 'chpart'=1 .
c     grid properties are transferred via common /grid/.
c
c      the fields 'sqb' and 'sqm' of the common /potinv/ are used as
c      work-space for solving the laplacian equation.
c
c     the resulting coulomb potential is returned on 'coulpt'.
c
      include 'skypar.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'potinv.inc'
      include 'constn.inc'
c
      dimension crho(kg1),coulpt(kg1)
c
c     some constants needed below
c
      data     vpeq3,        ezw,        fsx
     &    /6.0317561,0.083333333,0.833333333/
c
c-------------------------------------------------------------------------------
c
c     compute total charge from density 'crho' and fractio 'chpart',
c     and prepare r-weighted charge on workspace 'sqb'.
c
      if(chpart.eq.zero) then
       do ir=1,ngrid
         coulpt(ir) = zero
       enddo
      else
       vor    = 3.0*vpeq3*rstep*chpart
       charge = zero
       do i=1,ngrid
         charge = ((i-1)*(i-1))*crho(i) + charge
         sqb(i) = vor*crho(i)*(i-1)
       enddo
       charge = fourpi*rstep**3*charge
       if(iwr.ge.0)
     &  write(6,'(a,f11.5)') ' total charge from rhochr=',charge
c
c      part of green's function from inside to outside on 'sqm'
c
       intm   = ngrid-1
       sqm(1) = zero
       do i=2,intm
         sqm(i) = ezw*(sqb(i+1)+sqb(i-1))+fsx*sqb(i)
     &         +sqm(i-1)*(i-2)/(i-one)
       enddo
c
c      loop backward and collect r-weighted coulomb potential on 'sqb'
c
       rstep2 = rstep*rstep
       sqb(ngrid) = ( ezw*sqb(intm)+fsx*sqb(ngrid)
     &               +(sqm(intm)*(ngrid-2))/intm   )
     &                         *(rstep2*intm)/ngrid
     &             +  (1.440000*charge*ngrid)/ngridp
       do i=intm,2,-1
         sqb(i)=(sqm(i)*rstep2+sqb(i+1))*(i-one)/i
       enddo
       sqb(1) = zero
c
c      remove r-weight to obtain the final result
c
       coulpt(1) = (8.0*sqb(2)-sqb(3))/(6.0*rstep)
       do i=2,ngrid
         coulpt(i) = sqb(i)/((i-1)*rstep)
       enddo
      endif
c
      return
      end
c-----welfu ------------ part of the hartree-fock package ----------------------
c
      subroutine welfu
c
      include 'params.inc'
c
c     * * * *   iteration step for the single particle wave functions  * * * *
c
c     the iteration consists in applying the inverse mean-field
c     hamiltonian:
c        |new-'wf'> = (h - 'umin')**(-1) |old-'wf'>
c     where h is the mean-field hamiltonian and 'umin' the
c     minimum of the actual potential (including zentrifugal force).
c
c     the routine loops through all single particle states with isospin 'iso',
c     composed ang.mom. 'js' and radial quantum number 'n'.
c     for given 'iso' the auxiliary fields 'sqb' and 'd2b' are
c     occupied. for additionally given 'js' the effective potential
c     is computed on 'potl'.
c
c     these pieces constitute the single particle hartree-fock hamiltonian
c     from which the single particle energies 'spenrg' are computed.
c     the single particle energies and the single particle kinetic energies
c     are added up to the shell-part of the total energy on 'eshl'.
c
c     the effective hamiltonian goes via common /potinv/ into the
c     routine 'invert' where one iteration step is performed by acting
c     with the inverse hamiltonian on the single particle wave function.
c     the iteration is completed then by orthogonalising the new set
c     of wave functions.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'numpar.inc'
      include 'result.inc'
      include 'potinv.inc'
      include 'constn.inc'
      common /ctest/edirec(ksl2),efkincut,ifkincut    ! to communicate 'spkin' (obsolete)
c
      data ifolds/0/
c
c-------------------------------------------------------------------------------
c
c     iso-loop
c      'isoadd' is a pointer for the potential and mass fields;
c      'isnpad' is a pointer for the wave function pointer 'npla';
c      on 'egx' is accumulated lateron the shell-energy.
c
      do 1000 iso=1,2
      isoadd = (iso-1)*ngrid
      isnpad = (iso-1)*kjmx
      egx    = zero
c
c     square root of inverse mass 'bmass' on 'sqb';
c     d/dr**2 'sqb' / 'sqb' on 'd2b';
c     d/dr('bmass')/r on 'd1b'.
c
      do i=1,ngrid
        sqb(i) = sqrt(bmass(i+isoadd))
      enddo
c
      call derv25(sqb(1),d2b(1),one,ngrid,rstep)
      do i=1,ngrid
        d2b(i) = d2b(i)/sqb(i)
      enddo
c
      call derv15(bmass(isoadd+1),d1b(1),one,ngrid,rstep)
      d1b(1)=(16.0*bmass(2+isoadd)-15.0*bmass(1+isoadd)
     *        -bmass(3+isoadd))/(6.0*rstep**rstep)
      do i=2,ngrid
        d1b(i) = d1b(i)*rinvs(i)
      enddo
c
c     loop over composed ang.mom. 'js'
c      'll' is the orbital ang.mom.;
c      'sig' is the parity of the wave function 'wf' (it is -parity
c      of the state because 'wf' is multiplied by r);.
c      'phasjs' is the degeneracy of the shell.
c      'jsnpl' is pointer to to nucleon shell in pointer field 'npla'.
c
      do 40 js=1,jsmax(iso)
      jsnpl = js+isnpad
      nmax  = npla(1,jsnpl)
      if(nmax.gt.0) then
       ll     = js/2
       sig    = float(mod(ll,2)*2-1)
       xlp    = float(ll*ll+ll)
       xl0    = wls(js)
       phasjs = float(js-ll)
c
c      the effective potential for given 'js' and 'iso' on 'potl';
c      its minimum on 'umin';
c      the zentrifugal force on 'aux3'.
c
       aux3(1) = zero
       potl(1) = pot(1+isoadd)+xl0*potls(1+isoadd)+d1b(1)
       umin   = 1000.0
       do i=2,ngrid
         n      = i+isoadd
         y      = xlp*zentrf(i)
         aux3(i) = y
         y       = pot(n)+xl0*potls(n)+bmass(n)*y+d1b(i)
         umin   = min(y,umin)
         potl(i) = y
       enddo
c
c      loop over radial quantum number 'n'.
c       accumulating energies and performing one w.f. iteration.
c       'npl' is wave function pointer,
c       'nplshl' is the number of the state in harmonic oscillator ordering.
c
       do 70 n=1,nmax
       npl    = npla(n+1,jsnpl)
       nplm   = npl-1
       nplshl = npl/ngrid+1
       if(ifolds.eq.1) then
c                                         single particle energy
        spe    = (wfs(npl)*wfs(npl))*bmass(1+isoadd)*0.5
        do i=2,ngrid
          spe    = bmass(i+isoadd)*(wfs(nplm+i)*wfs(nplm+i))
     &            +potl(i)*(wf(nplm+i)*wf(nplm+i))          + spe
        enddo
        spe    = rstep*spe
        spenrg(nplshl) = spe
        IF(tphpair) THEN
          egx    = spe*(weight(nplshl)*phasjs) + egx
        ELSE
          egx    = spe*(weight(nplshl)*phpair(nplshl)*phasjs) + egx
        END IF
c                                         s.p. kinetic energy
        spekin = (wfs(npl)*wfs(npl))*0.5
        do i=2,ngrid
          spekin = (wfs(i+nplm)*wfs(i+nplm))
     &            +aux3(i)*(wf(i+nplm)*wf(i+nplm)) + spekin
        enddo
        spekin = rstep*spekin*amfc(iso)
        spkin(nplshl) = spekin
!        IF(ifkincut==1) edirec(nplshl) = spekin  ! communication to 'pairdn'
        IF(tphpair) THEN
!          egx    = spekin*amfc(iso)*(weight(nplshl)*phasjs) + egx
          egx    = spekin*(weight(nplshl)*phasjs) + egx
        ELSE
          egx    = spekin*(weight(nplshl)*phpair(nplshl)*phasjs) + egx
        END IF
       endif
c                                         w.f. iteration by h-inversion
       if(n.eq.1) then
        irep   = -1
       else
        irep   = +1
       endif
       call invert(aux3,wf(npl),irep)
       if(ifolds.eq.1) then
        do i=2,ngrid
          wf(nplm+i) = aux3(i)
        enddo
       else
        spe2   = zero
        do i=2,ngrid
          spe2   = aux3(i)*wf(nplm+i) + spe2
          wf(nplm+i) = aux3(i)
        enddo
        spe2   = umin+one/(rstep*spe2)
        spenrg(nplshl) = spe2
        IF(tphpair) THEN
          egx    = spe2*(weight(nplshl)*phasjs) + egx
        ELSE
          egx    = spe2*(weight(nplshl)*phpair(nplshl)*phasjs) + egx
        END IF
       endif
c
   70  continue
c
c      schmidt-orthonormalisation
c
       call orthon(npla(2,jsnpl),nmax)
c
c      derivative of 'wf' on 'wfs'
c
       do n=1,nmax
         npl    = npla(n+1,jsnpl)
         call derv15(wf(npl),wfs(npl),sig,ngrid,rstep)
         IF(thfbfold) CALL Gauss_folding_m(wf(npl),wffold(npl),rstep,
     &                           ngrid,js,hfbwidth,MAXVAL(jsmax(1:2)))
       enddo
c
      endif
   40 continue
c
c     store total shell-energy on 'eshl'
c
      if(ifolds.ne.1) then
       ekino   = zero
       r       = zero
       do i=2,ngrid
         r       = rstep+r
         ekino   = r*r*amfc(iso)*tau(isoadd+i) + ekino
       enddo
       ekino   = half*frpirs*ekino
       egx     = egx + ekino
      endif
c
      eshl(iso)=egx
c
 1000 continue
c
      return
      end
c-----welfdr------------ part of the hartree-fock package ----------------------
c
      subroutine welfdr(iter)
c
      include 'params.inc'
c
c     * * *   wave function iteration with direct action of 'hpsi'   * * *
c
c     the routine uses the gradient step with kinetic-energy damping.
c     it performs one iteration step as:
c
c                                'x0damp'
c      |wf_new> = |wf_old> - ----------------- [h-'spenrg(wf)'] |wf_old>
c                             t_kin + 'e0damp'
c
c     where t_kin = p**2/2m is the operator of the kinetic energy,
c     h is the mean-field hamiltonian, 'x0damp' and 'e0damp' are step
c     parameters, and 'spenrg(wf)' is the single particle energy of
c     the wave function, 'spenrg' = <wf|h|wf> .
c
c     the single particle energies and the single particle kinetic energies
c     are added up to the shell-part of the total energy on 'eshl'.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'potinv.inc'
      include 'linres.inc'
      include 'constn.inc'
      common /ctest/edirec(ksl2),efkincut,ifkincut   ! to communicate 'spkin'
c
!      data v2limh/1.0E-6/                 ! lower limit in HFB step
ctest      data v2limh/1.0/                 ! lower limit in HFB step
!      data delmix/0.8/                   ! factor on gap potential in damping
ctest      data delmix/0.5/                   ! factor on gap potential in damping
c      data v2damp/0.0/
c
c-------------------------------------------------------------------------------
c
c     iso-loop
c      'isoadd' is a pointer for the potential and mass fields;
c      'isnpad' is a pointer for the wave function pointer 'npla';
c      on 'egx' is accumulated lateron the shell-energy.
c
      ekin   = zero
      do 1000 iso=1,2
      isoadd = (iso-1)*ngrid
      isnpad = (iso-1)*kjmx
      egx    = zero
c      write(6,'(a,g12.4)') ' pptmax=',pptmax(iso)
c
c     square root of inverse free mass on 'sqb';
c     d/dr**2 'sqb' / 'sqb' on 'd2b';
c     d/dr 'bmass' /r on 'd1b'.
c
      call derv15(bmass(isoadd+1),d1b(isoadd+1),one,ngrid,rstep)
      sqfree = sqrt(amfc(iso))
      do i=1,ngrid
        d2b(i) = zero
        sqb(i) = sqfree
      enddo
c
c     loop over composed ang.mom. 'js'
c      'll' is the orbital ang.mom.;
c      'sig' is the parity of the wave function 'wf' (it is -parity
c      of the state because 'wf' is multiplied by r);.
c      'phasjs' is the degeneracy of the shell.
c
      do 40 js=1,jsmax(iso)
      jsnpl = js+isnpad
      nmax  = npla(1,jsnpl)
      if(nmax.gt.0) then
       ll     = js/2
       sig    = float(mod(ll,2)*2-1)
       xlp    = float(ll*ll+ll)
cold       xl0    = wls(js)
       phasjs = float(js-ll)
c
c      prepare effective potential for kinetic-energy inversion
c
       xlfac  = sqfree**2*xlp
       do i=1,ngrid
         potl(i) = e0damp + xlfac*zentrf(i)
       enddo
       umin   = zero
c
c
c      loop over radial quantum number 'n'.
c       accumulating energies and performing one w.f. iteration.
c       'npl' is wave function pointer,
c       'nplshl' is number of state in harmonic oscillator ordering.
c
       do 70 n=1,nmax
       npl    = npla(n+1,jsnpl)
       nplm   = npl-1
       nplshl = npl/ngrid+1
c
c      s.p. kinetic energy
c
!       call hpsi(wf(npl),aux3,wfs(npl),aux2,0,js,iso,zero)
       call tpsi(wf(npl),aux3,aux2,js,iso,zero)
       xnorm  = zero
       spekin = zero
       spekin2 = 0D0
       do i=2,ngrid
         ri     = rinvs(i)
         wfa    = wf(nplm+i)
         wfas   = wfs(nplm+i)
         wfa2   = wfa*wfa
         xnorm  = wfa*wfa*rho(isoadd+i) + xnorm
         spekin = wf(i+nplm)*aux3(i) + spekin
!         spekin2 = (wfas*wfas - ((wfa*wfas)+(wfa*wfas)-xlp*wfa2*ri)*ri)
         spekin2 = (wfas*wfas + xlp*wfa2*ri*ri)*rho(isoadd+i)
     &            + spekin2
       enddo
       xnorm  = rstep*xnorm
       spekin = rstep*spekin
       spkin(nplshl)=spekin !/xnorm
       ifkincut=0   ! test
       IF(ifkincut==1) THEN
         edirec(nplshl) = rstep*spekin2*amfc(iso)/xnorm  ! communication to 'pairdn'
       ELSE
         edirec(nplshl) = 0D0
       END IF
       IF(tphpair) THEN
         ekin   = spekin*(weight(nplshl)*phasjs) + ekin
         egx    = spekin*(weight(nplshl)*phasjs) + egx
       ELSE
         ekin   = spekin*(weight(nplshl)*phpair(nplshl)*phasjs) + ekin
         egx    = spekin*(weight(nplshl)*phpair(nplshl)*phasjs) + egx
       END IF
c
c      action of the mean-field hamiltonian
c
       call hpsi(wf(npl),aux3,wfs(npl),aux2,0,js,iso,zero)
c
c      single particle energy and energy-fluctuation
c
       xnorm  = zero
       ewid   = zero
       spe    = zero
       do i=2,ngrid
         xnorm  = wf(nplm+i)*wf(nplm+i) + xnorm
         ewid   = aux3(i)*aux3(i) + ewid
         spe    = wf(nplm+i)*aux3(i) + spe
       enddo
       xnorm  = rstep*xnorm
       ewid   = rstep*ewid/xnorm
       spe    = rstep*spe/xnorm
       spenrg(nplshl) = spe
       hpswid(nplshl) = SQRT(MAX(ewid-spe*spe,1D-180))
       IF(tphpair) THEN
         egx    = spe*(weight(nplshl)*phasjs) + egx
       ELSE
         egx    = spe*(weight(nplshl)*phpair(nplshl)*phasjs) + egx
       END IF
       if(nplshl.eq.1 .and. iwr.gt.-1) then
        ewid   = sign(sqrt(abs(ewid - spe*spe)),ewid-spe*spe)
        write(6,'(a,g12.5,a,g12.5)') '     testprint ewid=',ewid,
     &     '  ,  xnorm=',xnorm
       endif
c
c      optionally variance of c.m. operator
c
ctest       if(iflast.eq.1) 
ctest     &   call cmact(wf(npl),aux3,wfs(npl),aux2,js,iso,
ctest     &              weight(nplshl))
       if(ifzpe.eq.2) then
         call cmact(wf(npl),aux3,wfs(npl),aux2,js,iso,
     &              weight(nplshl),specm(nplshl),delcm(nplshl))
         xnorm  = zero
         spe    = zero
         do i=2,ngrid
           xnorm  = wf(nplm+i)*wf(nplm+i) + xnorm
           spe    = wf(nplm+i)*aux3(i) + spe
         enddo
         spe    = spe/xnorm
       endif
c
c      w.f. iteration by t-inversion of [h-'spenrg']|wf>
c
       if(iparac.le.5 .or. iparac.eq.7) then
         do ir=1,ngrid
           aux3(ir) = aux3(ir) - spe*wf(nplm+ir)
         enddo
         x0act  = x0damp
       else
!         v2act  = max(v2occ(nplshl),v2limh)
         v2act  = max(weight(nplshl),v2limh)
         e0act  = e0damp+sqrt((one-v2act)/v2act)*pptmax(iso)*delmix
cold         x0act  = x0damp/(v2act)**(0.5-v2damp)
         x0act = x0damp/v2act**(1.0-ATAJIMA)
         xlfac  = sqfree**2*xlp
         speact=0D0
         do ir=1,ngrid
           aux3(ir) = savlr(nplm+ir,1)          ! retrieved from HFB code
           potl(i) = e0act + xlfac*zentrf(i)
           speact = speact + wf(nplm+ir)*aux3(ir)
         enddo
         speact=speact*rstep
         aux3(1:ngrid)=aux3(1:ngrid)-speact*wf(npl:nplm+ngrid)
       endif
       irep   = 2*min(2,n)-3
       call invert(aux3,aux3,irep)
       do i=1,ngrid
         wf(nplm+i) = wf(nplm+i) - x0act*aux3(i)
       enddo
c
   70  continue
c
c      schmidt-orthonormalisation
c
       call orthon(npla(2,jsnpl),nmax)
c
c      derivative of 'wf' on 'wfs'
c
       do n=1,nmax
         npl    = npla(n+1,jsnpl)
         call derv17(wf(npl),wfs(npl),sig,ngrid,rstep)
         IF(thfbfold) CALL Gauss_folding_m(wf(npl),wffold(npl),rstep,
     &                         ngrid,js,hfbwidth,MAXVAL(jsmax(1:2)))
       enddo
c
      endif
   40 continue
c
c     store total shell-energy on 'eshl'
c
      eshl(iso)=egx
c
 1000 continue
      ekin = ekin+ekin
c
      return
      end
c-----orthon------------ part of the hartree-fock package ----------------------
c
      subroutine orthon(npa,nmax)
c
c     schmidt-orthonormalisation of the wave functions belonging to
c     one angular momentum and given by the pointer field 'npa'.
c     'nmax' is the number of wave function to be orthonormalised.
c
      include 'params.inc'
c
      dimension npa(knmx)
c
      include 'wavfun.inc'
      include 'grid.inc'
c
c-------------------------------------------------------------------------------
c
c     loop over the states to be orthonormalised
c
       do 110 n=1,nmax
       npl    = npa(n)
       nplm   = npl-1
       nple   = npl+ngrid-1
       if(n.gt.1) then
c
c      loop over the previous states, already orthonormalised
c
        do 120 n2=1,n-1
        npl2   = npa(n2)
        npl2m  = npl2-1
c
c       overlap <wf_n|wf_n2>
c
        wfovlp = wf(npl)*wf(npl2)*0.5
        do i=2,ngrid
          wfovlp = wf(i+nplm)*wf(i+npl2m)+wfovlp
        enddo
        wfovlp = wfovlp*rstep
c
c       project piece of |wf_n2> out of |wf_n>
c
        xnorm  = -wfovlp
        ndif   = npl2-npl
        do i=npl,nple
          wf(i)  = xnorm*wf(i+ndif) + wf(i)
        enddo
  120   continue
       endif
c
c      normalise |wf_n>
c
       wfnor  = wf(npl)*wf(npl)*0.5
       do i=2,ngrid
         wfnor  = wf(i+nplm)*wf(i+nplm)+wfnor
       enddo
       wfnor  = wfnor*rstep
c
       xnorm  = 1./sqrt(wfnor)
       do i=npl,nple
         wf(i)  = xnorm*wf(i)
       enddo
c
c
  110  continue
c
       return
       end
c-----zpecor------------ part of the hartree-fock package ----------------------
c
      subroutine zpecor
c
      include 'params.inc'
c
c     * * * * * *    center-of-mass correction for the energy    * * * * * *
c
c     all necessary information is handled via common's.
c     input is:
c      ifzpe    is the switch, entered via 'numpar'
c      wf       are the wave functions, eneterd via 'wavfun'
c     and further information on shell ordering, grid etc, comes
c     via the according common's (see there for further explanations).
c
c     output is:
c      zpe      the zero-point energy, returned via 'result'
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'result.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c-----------------------------------------------------------------------
c
c     check whether 'zpe' is wanted
c
      if(ifzpe.lt.1 .and. ifzpe.gt.-3) then
       zpe    = zero
      else
c
c-------------------------------------------------------------------------------
c
c      full evaluation of momentum width <p**2> and subsequent classical
c      zpe correction.
c
c-------------------------------------------------------------------------------
c
c      isospin loop. the same isospin applies for the 'a'-state and
c      the 'b'-state, because the p-operator does not mix isospin.
c
       diag   = 0.0
       offdig = 0.0
       do 10 iso=1,2
       isnpad = (iso-1)*kjmx
c
c      'a'-loop over the angular momenta 'jsa'.
c        the derived quantities are:
c         'xja'   = total angular momentum
c         'xla'   = orbital angular momentum
c         'la'    =    "       "       "      as integer
c
       do 20 jsa=1,jsmax(iso)
       nmaxa  = npla(1,jsa+isnpad)
       if(nmaxa.le.0) goto 20
       xja    = float((jsa-1)/2)+0.5
       la     = jsa/2
       xla    = one*la
       xla2   = xla*xla+xla
c
c      the diagonal contribution from p**2.
c
       do 30 nshla=1,nmaxa
       npta   = npla(nshla+1,jsa+isnpad)
       nptma  = npta-1
       sum    = 0.0
       do i=2,ngrid
         wfsa = wfs(i+nptma)-wf(i+nptma)/((i-1)*rstep)
c        sum    = ( (wfs(i+nptma)*wfs(i+nptma))
         sum    = ( wfsa*wfsa
     *             +(wf(i+nptma)*wf(i+nptma))*xla2*zentrf(i) )
     *           + sum
       enddo
       sum    = rstep*sum
       v2a    = weight(npta/ngrid+1)
       if(temprt.gt.zero) then
        ttana  = ttthrm(npta/ngrid+1)
        v2a    = (v2a-0.5*(one-ttana))/ttana
        v2a    = (one+ttan)*0.5*v2a
c         write(6,'(a,3(1pg12.4))')
c     &    ' v2a,newfac,thfac=',weight(npta/ngrid+1),v2a,
c     &    (one+ttthrm(npta/ngrid+1))*0.5         
       endif
       IF(tphpair) THEN
         diag   = v2a*(xja+xja+one)*sum + diag
       ELSE
         diag   = v2a*(xja+xja+one)*sum*phpair(npta/ngrid+1) + diag
       END IF
   30  continue
c
c      the loop over the 'b'-angular momenta 'jsb'.
c        the derived quantities are:
c         'xjb'   = total angular momentum
c         'xlb'   = orbital angular momentum
c         'lb'    =    "       "       "      as integer
c
       do 40 jsb=1,jsa
       nmaxb  = npla(1,jsb+isnpad)
       if(nmaxb.le.0) goto 40
       xjb    = float((jsb-1)/2)+0.5
       lb     = jsb/2
       xlb    = one*lb
c
c      check triangular relation for angular momentum selection rules
c
       if(mod(la+lb,2).ne.1) goto 40
       if(iabs(la-lb).gt.1) goto 40
       if(abs(xja-xjb).gt.one) goto 40
c
c      compute factors for reduced matrix element.
c
       if(iabs(la-lb).eq.1) then
        if(xja-xla.eq.xjb-xlb) then
         sixfcu = sqrt((xla+xjb+2.5)*(xla+xjb-0.5)
     *                 /((xla+xla+one)*(xlb+xlb+one)))
        else
         sixfcu = sqrt((xjb-xla+1.5)*(xla-xjb+1.5)
     *                 /((xla+xla+one)*(xlb+xlb+one)))
     *            *2.0*(xjb-xlb)
        endif
        if(xla.gt.xlb) then
         facdru = sqrt(xla)
         facctu = -xla
        else
         facdru = -sqrt(xlb)
         facctu = xlb
        endif
       else
        sixfcu = zero
        facdru = zero
        facctu = zero
       endif
c
c      the radial matrix elements for given jsa and jsb, sums running
c      over all nshla (to jsa) and nshlb (to jsb).
c
       do 50 nshla=1,nmaxa
       npta   = npla(nshla+1,jsa+isnpad)
       nptma  = npta-1
       v2a    = weight(npta/ngrid+1)
       if(temprt.gt.zero) then
        ttana  = ttthrm(npta/ngrid+1)
        v2a    = (v2a-0.5*(one-ttana))/ttana
c        v2a    = (one+ttan)*0.5*v2a
       endif
!                    sloppy in case of blocking
       uva    = sqrt(max(v2a-v2a*v2a,xsmall))
c
       do 50 nshlb=1,nmaxb
       nptb   = npla(nshlb+1,jsb+isnpad)
       nptmb  = nptb-1
       v2b    = weight(nptb/ngrid+1)
       if(temprt.gt.zero) then
        ttanb  = ttthrm(nptb/ngrid+1)
        v2b    = (v2b-0.5*(one-ttanb))/ttanb
c        v2a    = (one+ttan)*0.5*v2a
       endif
       uvb    = sqrt(max(v2b-v2b*v2b,xsmall))
c
       radmtu = zero
       do i=2,ngrid
         ra     = (i-1)*rstep
         radmtu = wf(i+nptmb)
     *                   *(-wfs(i+nptma)+facctu*wf(i+nptma)/ra)
     *            + radmtu
       enddo
       radmtu = rstep*radmtu
       if(temprt.gt.zero) then
         wtha  = (one-ttthrm(npta/ngrid+1))*0.5
         wthb  = (one-ttthrm(nptb/ngrid+1))*0.5
         wthca = one-wtha
         wthcb = one-wthb
         wfac = wtha*wthb*(-0.5*(v2a+v2b)+(v2a*v2b+uva*uvb))
     &          +wthb*wthca*v2a
     &          +wthca*wthcb*(v2a*v2b+uva*uvb)
       else
         IF(tphpair) THEN
           wfac = (v2a*v2b+uva*uvb)
         ELSE
           pha    = phpair(npta/ngrid+1)
           phb    = phpair(nptb/ngrid+1)
           wfac = (v2a*v2b+uva*uvb)*pha*phb
         END IF
       endif
       offdig = (radmtu*facdru*sixfcu)**2*wfac
     *         + offdig
c
   50  continue
   40  continue
   20  continue
   10  continue
c
c      finally assembling the zpe-correction  (in MeV).
c
       pwidth = (diag-offdig-offdig)
       zpe    = pwidth*half*(amfc(1)+amfc(2))/anumb
c
c      optionally comparision with simple estimates for 'zpe'
c
       if(iwr.ge.-1) then
        write(06,'(5(a,1pg12.5))') ' zpe=',zpe,' pwidth=',pwidth,
     &                ' diag=',diag*half*(amfc(1)+amfc(2))/anumb,
     &                ' offdig=',offdig*half*(amfc(1)+amfc(2))/anumb,
     &                ' ekin/A=',ekin/anumb
c
c      average zpe adjusted to skyrme-hartree-fock results
c       (reproduces trends in particle number, misses isotopic effects)
c
        zpe1  = 17.3*anumb**(-0.2)
c
c      estimate from the harmonic oscillator shell model
c       (only a rough approach)
c
        zpe2  = .75*41.*anumb**(-0.3333333)
        write(06,'(2(a,g12.5))') ' schematic:  zpe(power 1/3)=',zpe2,
     *            ' zpe(power 0.2)=',zpe1
       endif
c
      endif
c
      return
      end
c     include 'radiusnw.for'
c-----radius-----------------------------------------------------------
c
      subroutine radius(rhop,rhon,sjsp,sjsn,
     &                  nprot,nneut,rstep,ngrid,amfpr,amfne,zpe,ekin,
     &                  ifzpe,iwr,ifrso,ifform,
     &                  formf,rhochr,  foc,fop,fon,fo,foc0,fow,
     &                  rmsc,rmsp,rmsn,rms,
     &                  rdmsc,rdmsp,rdmsn,rdms,
     &                  surfc,surfp,surfn,surf,
     &                  rsecc,rsecp,rsecn,rsec,
     &                  smomc,smomp,smomn,smom, rhominmax,str_squark)
c
      include 'params.inc'
!      include 'skypar.inc'
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
c      ifzpe   =-1 --> factor (1-1/a) on the mass is used for correction
c                      c.m. correction from actual ZPE
c                0 --> factor (1-1/a) on the mass is used for correction
c                      c.m. correction uses average guess
c                1 --> correction is computed explicitely and a
c                      posteriori as <p**2>/(2*a*m).
c                2 --> correction is computed explicitely and before
c                      variation as <p**2>/(2*a*m).
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
      dimension foc(kfrm),fop(kfrm),fon(kfrm),fo(kfrm),foc0(kfrm),
     &          fow(kfrm)
c
      dimension rhop(ngrid),rhon(ngrid)         ! input information
      dimension sjsp(ngrid),sjsn(ngrid)         ! input information
      dimension rhochr(ngrid),formf(kfrm)       ! results
!      REAL(8) :: aux(ngrid)
c
      common /radpar/delq,rleng,iqmax           ! to communicate to subr.
c
c      data  zero,one,pi/0.0,1.0,3.1415926/
      data  zero,half,one,pi/0D0,0.5D0,1D0,3.141592653589793D0/
      LOGICAL, PARAMETER :: ttestrad=.FALSE.
c
c-------------------------------------------------------------------------------
c
c     preparation of auxiliary parameters
c
      anumb  = one*(nprot+nneut)
c      kqmax  = min(kfrm,ngrid*3/4)       ! ???
      iqmax  = min(kfrm,ngrid) !*3/4)
c      iqmax  = ngrid-1
      rleng  = (ngrid-1)*rstep
      delq   = pi/rleng
      frpirs = 4*pi*rstep
c
c     compute 'cmwid' as width for for c.m.correction on the formfactor
c
      if(ifzpe.eq.0) then
        cmwid  = 1.58/(6.0*anumb**0.66666667)
        if(iwr.ge.-1) 
     &     write(6,'(a,2(1pg13.5))') ' cmwid=',cmwid
      elseif(ifzpe.eq.-1) then
        zpedig = ekin/anumb                    ! effective diagonal ZPE
        cmwid  = 0.375*0.5*(amfpr+amfne)/(anumb*zpedig)
        if(iwr.ge.-1) 
     &     write(6,'(a,2(1pg13.5))') ' ZPE,cmwid=',zpedig,cmwid
      elseif(ifzpe.le.-2 .and. ifzpe.ne.-4) then
        cmwid = 0.0D0
      else
        cmwid  = 0.375*0.5*(amfpr+amfne)/(anumb*zpe)
        if(iwr.ge.-1) 
     &     write(6,'(a,2(1pg13.5))') ' ZPE,cmwid=',zpe,cmwid
      endif
      cmwidsav=cmwid
!      write(6,'(a,4(1pg13.5))') 'ZPE etc=',zpe,cmwid,rleng,amfpr+amfne
c
c      formfactor by fourier-bessel transformation
c
c                                                 q=0.0 separately
      call fourf(0D0,
     &             rhop,rhon,sjsp,sjsn,
     &             rstep,ngrid,ifrso,ifform,cmwid,
     &             foc0(1),fop(1),fon(1),fo(1),foc(1))
      sump   = zero
      sumn   = zero
c      sumsp  = zero
      r      = zero
      do ir=2,ngrid
        r     = rstep+r
        sump   = r*r*rhop(ir)+sump
        sumn   = r*r*rhon(ir)+sumn
      enddo
      foc0(1) = frpirs*sump
      fop(1)  = foc0(1)
      fon(1)  = frpirs*sumn
      fo(1)   = fop(1)+fon(1)
      foc(1)  = foc0(1)
      sump   = zero
      sumn   = zero
c      sumsp  = zero
      r      = zero
      do ir=2,ngrid
        r     = rstep+r
        sump   = r*r*sjsp(ir)+sump
        sumn   = r*r*sjsn(ir)+sumn
      enddo
      WRITE(6,'(a,2(1pg13.5))') 'average ls-density:',
     &      frpirs*sump,frpirs*sumn
c                                                      all other q
      do iq=2,iqmax
        q1     = iq*delq-delq
        call fourf(q1,
     &             rhop,rhon,sjsp,sjsn,
     &             rstep,ngrid,ifrso,ifform,cmwid,
     &             foc0(iq),fop(iq),fon(iq),fo(iq),foc(iq))
      enddo
c
c     compute weak-charge form factor
c
!      WRITE(2,'(a,1pg13.5)') 'RADIUS: str_squark=',str_squark
!      CALL formf_weak(fop,fon,str_squark,delq,iqmax,fow)
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
!      rsecc  = 7.7253/qzero(foc,int(4.493/(delq*rdmsc))+2)
!      istart2 = nint(5.6D0/(delq*rdmsc))
      istart2 = int(4.493D0/(delq*rdmsc))+1
      rsecc  = 7.7253/qzero(foc,istart2+1)
!      WRITE(6,'(a,i5,3(1pg13.5))') 'NEXT INT:',
!     &   istart2,delq,qzero(foc,1),rsecc
!      rsecc  = 7.7253/qzero(foc,nint(5.6D0/(delq*rdmsc)))
      if(iwr.ge.-1) then
!      if(iwr.ge.-2) then
       rdmsc0 = 4.493/qzero(foc0,1)
       rsecc0 = 7.7253/qzero(foc0,int(4.493/(delq*rdmsc0))+2)
!       rsecc0  = 7.7253/qzero(foc0,nint(5.6/(delq*rdmsc0)))
       rdms   = 4.493/qzero(fo,1)
       rsec   = 7.7253/qzero(fo,int(4.493/(delq*rdms))+2)
!       rsec  = 7.7253/qzero(fo,nint(5.6/(delq*rdms)))
       if(nprot.gt.0) then
        rdmsp  = 4.493/qzero(fop,1)
        rsecp  = 7.7253/qzero(fop,int(4.493/(delq*rdmsp))+2)
!        rsecp  = 7.7253/qzero(fop,nint(5.6/(delq*rdmsp)))
       else
        rdmsp  = zero
        rsecp  = zero
       endif
       if(nneut.gt.0) then
        rdmsn  = 4.493/qzero(fon,1)
        rsecn  = 7.7253/qzero(fon,int(4.493/(delq*rdmsn))+2)
!        rsecn  = 7.7253/qzero(fon,nint(5.6/(delq*rdmsn)))
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
      if(nprot.gt.0) then
       q1     = 5.600/rdmsc
       qr     = q1*rdmsc
       surfc  = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &              *3.*nprot/(qr*qr*fbint(q1,foc))))))
       surfc0 = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &            *3.*nprot/(qr*qr*fbint(q1,foc0))))))
!       surfe  = acosh(abs((sin(qr)/qr-cos(qr))*3D0*nprot/
!     &                      (qr*qr*fbint(q1,foc))))/q1
       helmpower=2D0
!       CALL find_surfe(fbint(q1,foc)/nprot,q1,rdmsc,0.1D0,3D0,surfe,
!     &                 helmpower)
      else
       surfc  = zero
       surfc0 = zero
      endif
      if(iwr.ge.-1) then
!      if(iwr.ge.-2) then
       q1     = 5.600/rdms
       qr     = q1*rdms
       surf   = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &           *3.*(nprot+nneut)/(qr*qr*fbint(q1,fo))))))
       if(nprot.gt.0) then
        q1     = 5.600/rdmsp
        qr     = q1*rdmsp
        surfp  = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &                *3.*nprot/(qr*qr*fbint(q1,fop))))))
       else
        surfp  = zero
       endif
       if(nneut.gt.0) then
        q1     = 5.600/rdmsn
        qr     = q1*rdmsn
        surfn  = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &             *3.*nneut/(qr*qr*fbint(q1,fon))))))
       else
        surfn  = zero
       endif
      endif
c
c     rms-radii from curvature of formfactor at q=0
c     ---> dangerous procedure, step 'dqt' should not be too small
c     ---> values will be overwritten later by direct integration
c
      dqt    = delq/2D0 !/4.0 !10.0     !  delq/4.0
      dqtd2  = 6D0/(dqt*dqt)
      rmsc   = sqrt(abs(dqtd2*(1.25D0*foc(1)
     &                  -1.33333333D0*fbint(dqt,foc)
     &                  +fbint((dqt+dqt),foc)/12D0)/foc(1)))
!        WRITE(78,*) 0D0,foc(1)
!        DO i=1,1000
!        call fourf(i*dqt,
!     &             rhop,rhon,sjsp,sjsn,
!     &             rstep,ngrid,ifrso,ifform,cmwid,
!     &             xfoc0,xfop,xfon,xfo,xfoc)
!        WRITE(78,*) i*dqt,xfoc,abs((dqtd2/i**2)*(foc(1)-xfoc)/foc(1))
!        END DO
!        call fourf(dqt,
!     &             rhop,rhon,sjsp,sjsn,
!     &             rstep,ngrid,ifrso,ifform,cmwid,
!     &             xfoc0,xfop,xfon,xfo,xfoc)
!        call fourf(1D-10,
!     &             rhop,rhon,sjsp,sjsn,
!     &             rstep,ngrid,ifrso,ifform,cmwid,
!     &             xfoc0,xfop,xfon,xfo,xfoc2)
!      WRITE(77,*) 'testprint F(0):',foc(1),xfoc,xfoc2
!!      rmsc   = sqrt(abs(dqtd2*(foc(1)-fbint(dqt,foc))/foc(1)))
!      rmsc   = sqrt(abs(dqtd2*(foc(1)-xfoc)/foc(1)))
      if(iwr.ge.-1) then
!      if(iwr.ge.-2) then
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
      rmscsav = rmsc
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
        sum    = qa*qa*foc(iq)+sum
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
!          s = SIN(r*qa)
!          c = COS(qa*r)
          sum    = qa*s*foc(iq)+sum
        enddo
        rhochr(ir) = (sum*delq*dpi3/r)
!        rhochr(ir) = abs(sum*delq*dpi3/r)
c
      enddo
c
c     test FOURB
c      
!      CALL fourb(foc,aux,delq,rstep,ngrid,iqmax)
!      acc=0D0
!      DO ir=2,ngrid
!        acc=(rhochr(ir)-aux(ir))**2+acc
!      END DO
!      WRITE(11,*) '# average deviation FOURB:',sqrt(acc/ngrid)
c
c     surface width from moments
c
      if(iwr.ge.-1) 
     &  write(6,'(a,4g13.5)')
     &  ' from F(q): rms,rmsp,rmsn,rmsc=', rms,rmsp,rmsn,rmsc
c
      rmomt  = zero
      rmom0c = zero
      rmom0p = zero
      rmom0n = zero
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
        rmomt  = rhochr(ir)*r2*6D0*(1D0-SIN(dqt*r)/(dqt*r))+rmomt
        rmom0p = r2*rhop(ir) + rmom0p
        rmom2p = r4*rhop(ir) + rmom2p
        rmom4p = r6*rhop(ir) + rmom4p
        rmom0n = r2*rhon(ir) + rmom0n
        rmom2n = r4*rhon(ir) + rmom2n
        rmom4n = r6*rhon(ir) + rmom4n
        rmom0c = r2*rhochr(ir) + rmom0c
        rmom2c = r4*rhochr(ir) + rmom2c
        rmom4c = r6*rhochr(ir) + rmom4c
      enddo
      rmomt  = rmomt/rmom0c
      rmom2  = frpirs*(rmom2p+rmom2n)/fo(1)
      rmom4  = frpirs*(rmom4p+rmom4n)/fo(1)
      rmom2p = rmom2p/rmom0p  !frpirs*rmom2p/fop(1)
      rmom4p = frpirs*rmom4p/fop(1)
      rmom2n = rmom2n/rmom0n  !frpirs*rmom2n/fon(1)
      rmom4n = frpirs*rmom4n/fon(1)
      rmom0c = frpirs*rmom0c
      rmom2c = frpirs*rmom2c/rmom0c !foc(1)
      rmom4c = frpirs*rmom4c/foc(1)
      rms    = sqrt(rmom2)
      rmsp   = sqrt(rmom2p)
      rmsn   = sqrt(rmom2n)
      rmsc   = sqrt(rmom2c)
      WRITE(77,*) foc(1),rmom0c,rmsc,rmom4c**0.25D0 !rmomt/dqt**2
c
      if(iwr.ge.-1) then
        write(6,'(a,6g13.5)')
     &   ' rmom2p,rmom4p,rmom2n,rmom4n,rmom2,rmom4=',
     &     rmom2p,rmom4p,rmom2n,rmom4n,rmom2,rmom4
        smomp  = (rmom4p- rmom2p*rmom2p*25.0/21.0)**0.25
        smomn  = (rmom4n- rmom2n*rmom2n*25.0/21.0)**0.25
        smomc  = (rmom4c- rmom2c*rmom2c*25.0/21.0)**0.25
        smom   = (rmom4 - rmom2*rmom2*25.0/21.0)**0.25
c
c       print results
c
        haloe  = rmsc-SQRT(0.6D0*rdmsc**2+3D0*surfe**2)
        haloe4 = rmsc-(3D0*rdmsc**4/7D0+6D0*rdmsc**2*surfe**2
     &                     +21.94D0*surfe**4)**0.25D0
        halo4  = rmsc-(3D0*rdmsc**4/7D0+6D0*rdmsc**2*surfc**2
     &                     +15D0*surfc**4)**0.25D0
        apar = 15D0
        afac = 13.5D0-0.7D0*apar
        rms3 = rmsc**2/afac
        drms4 = (25D0*rmsc**4-21D0*rmom4c)/(30D0*afac)
        sigma4 = SQRT(rms3-SQRT(rms3**2+drms4))
!        apar = 21.96D0
        apar = 5D0*(2D0*helmpower+3D0)/helmpower
        afac = 13.5D0-0.7D0*apar
        rms3 = rmsc**2/afac
        drms4 = (25D0*rmsc**4-21D0*rmom4c)/(30D0*afac)
        sigmae4 = SQRT(rms3-SQRT(rms3**2+drms4))
        WRITE(311,'(2i4,8(1pg13.5))') nprot,nneut,rmsc,rdmsc,surfc,
     &        surfe,haloe,rmsc-SQRT(0.6D0*rdmsc**2+3D0*surfc**2),
     &        sigma4,sigmae4
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
        write(06,'(16x,a,5(/1x,a,5g13.5))')
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
!       rmsn   = sqrt(rmom2n)
!       rmsp   = sqrt(rmom2p)
c
        write(06,'(1h0,a,50(/1x,10g13.5))')
     &     'testprint rhochr:', (rhochr(i),i=1,ngrid)
c
c      nucleon folding for rms radius directly
c
       rmomt = 0D0
       rlsmomp0 = 0D0      
       rlsmomn0 = 0D0      
       rlsmomp2 = 0D0      
       rlsmomn2 = 0D0      
       r      = zero
       DO ir=2,ngrid
         r     = rstep+r
         r2    = r*r
         r4    = r2*r2
         rmomt  = sjsp(ir)*r2*SIN(dqt*r)/(dqt*r)+rmomt
         rlsmomp0 = rlsmomp0 + sjsp(ir)*r2
         rlsmomn0 = rlsmomn0 + sjsn(ir)*r2
         rlsmomp2 = rlsmomp2 + sjsp(ir)*r4
         rlsmomn2 = rlsmomn2 + sjsn(ir)*r4
       END DO
       rmomt = -frpirs*(rmomt-rlsmomp0)*6D0
       rlsmomp0 = frpirs*rlsmomp0
       rlsmomn0 = frpirs*rlsmomn0
       rlsmomp2 = frpirs*rlsmomp2
       rlsmomn2 = frpirs*rlsmomn2
       WRITE(77,'(a,6(1pg13.5))') 'testprint ls-moments:',
     &    rlsmomp0,rlsmomn0,rlsmomp2,rlsmomn2,rmomt/dqt**2,dqt
       WRITE(6,'(a,6(1pg13.5))') 'testprint ls-moments:',
     &    rlsmomp0,rlsmomn0,rlsmomp2,rlsmomn2,rmomt/dqt**2,dqt
       r2intrp = 0.7297D0
       r2intrn = (-0.1174D0*nneut)/nprot
       darcor = 3D0*0.010989D0
       xmgne = -0.04197D0      ! magn.moment * 2 * darcor
       xmgpr = 0.06131D0
       r2lsp = -(xmgpr+0.010989D0)*rlsmomp2/nprot !(0.010989D0-xmgpr)*rlsmomp2/nprot
       r2lsn = -xmgne*rlsmomn2/nprot
       r2cmcor = -6D0*cmwidsav
       WRITE(6,'(a,a)') 'rmsc**2     rmsp**2     rintr2p    rtint2n   ',
     &     '   darfac    r2lsp       r2lsn    r2cmcor     rmscold**2'
       WRITE(6,'(a,9(1pg13.5))') 'rms_C^2:',
     &    rmsp*rmsp+r2intrp+r2intrn+darcor+r2lsp+r2lsn+r2cmcor,
     &    rmsp*rmsp,r2intrp,r2intrn,darcor,r2lsp,r2lsn,r2cmcor,rmsc**2
       WRITE(6,'(a,7(1pg13.5))') ' rms_C1:',
     &    SQRT(rmsp*rmsp+r2intrp+r2intrn+darcor+r2lsp+r2lsn+r2cmcor),
     &    rmsc,rmscsav,rmsc0,rmsp,
     &    SQRT(rmsp*rmsp+r2intrp+r2intrn+r2cmcor)
       r2lsp = -(xmgpr-0.010989D0)*rlsmomp2/nprot !(0.010989D0-xmgpr)*rlsmomp2/nprot
       WRITE(6,'(a,a)') 'rmsc**2     rmsp**2     rintr2p    rtint2n   ',
     &     '   darfac    r2lsp       r2lsn    r2cmcor     rmscold**2'
       WRITE(6,'(a,9(1pg13.5))') 'rms_C^2:',
     &    rmsp*rmsp+r2intrp+r2intrn+darcor+r2lsp+r2lsn+r2cmcor,
     &    rmsp*rmsp,r2intrp,r2intrn,darcor,r2lsp,r2lsn,r2cmcor,rmsc**2
       WRITE(6,'(a,7(1pg13.5))') ' rms_C2:',
     &    SQRT(rmsp*rmsp+r2intrp+r2intrn+darcor+r2lsp+r2lsn+r2cmcor),
     &    rmsc,rmscsav,rmsc0,rmsp,
     &    SQRT(rmsp*rmsp+r2intrp+r2intrn+r2cmcor)
       WRITE(6,'(a,2i4,7(1pg13.5))') ' rms_C3:',nprot,nneut,
     &    rmsc,
     &    SQRT(rmsp*rmsp+r2intrp+r2intrn+darcor+r2lsp+r2lsn+r2cmcor),
     &    rmsc0,rmsp
       r2lsp = -(xmgpr-2D0*0.010989D0)*rlsmomp2/nprot !(0.010989D0-xmgpr)*rlsmomp2/nprot
       WRITE(6,'(a,7(1pg13.5))') ' rms_Cr:',
     &    SQRT(rmsp*rmsp+r2intrp+r2intrn+darcor+r2lsp+r2lsn+r2cmcor),
     &    rmsc,rmscsav,rmsc0,rmsp,
     &    SQRT(rmsp*rmsp+r2intrp+r2intrn+r2cmcor)
c
c
      endif
!
!      IF(ttestrad) THEN
!!       compte sigma*l expectation values
!        sl2p = 0D0
!        DO nshl=1,nmaxp
!          js     = jsp(nshl)
!          sl2p = sl2p+weight(nshl)*zwj(js)*wls(js)
!        END DO
!        sl2n = 0D0
!        DO nshl=1,nmaxn
!          js     = jsp(nshl)
!          sl2n = sl2n+weight(nshl+nmaxp)*zwj(js)*wls(js)
!        END DO
!        WRITE(6,'(a,4(1pg13.5))') '<sigma*l>=',sl2p,sl2n,
!     &      rlsmomp2,rlsmomn2
!      END IF
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
c      ifrso   =  1 --> charge density with contribution from ls current - old
c                -1 --> charge density with contribution from ls current - new
c                 0 --> charge density without contribution from ls current 
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
      LOGICAL,PARAMETER :: ttestterm=.FALSE.      ! test of l*s contribution
      REAL(8),SAVE :: sump0,sumn0,sumsp0,sumsn0,
     &      schsmp0,schsmn0,darfac0,sumsachsp0,sumsachsn0

c
c----------------------------------------------------------------
c
!      WRITE(*,*) 'FOURF: cmwid=',cmwid
      frpirs = 4*pi*rstep
      IF(ttestterm) cmwid=0D0
      IF(abs(ifrso)>1) STOP 'IFRSO too large'
c
c     fourier-bessel transformation for proton- and neutron-density
c     given on field 'rho' on common /wavfun/.
c
      IF(q1==0D0) THEN
        r      = zero
        sump   = zero
        sumn   = zero
          do i=2,ngrid
            r      = rstep+r
            sump   = (r*r)*rhop(i)+sump
            sumn   = (r*r)*rhon(i)+sumn
          enddo
        sump   = frpirs*sump
        sumn   = frpirs*sumn
      ELSE
        dc     = cos(q1*rstep)
        ds     = sin(q1*rstep)
        r      = zero
        sump   = zero
        sumn   = zero
        c      = one
        s      = zero
          do i=2,ngrid
            r      = rstep+r
            x      = c*dc-s*ds
            s      = s*dc+c*ds
            c      = x
            sump   = (r*s)*rhop(i)+sump
            sumn   = (r*s)*rhon(i)+sumn
          enddo
        sump   = frpirs*sump/q1
        sumn   = frpirs*sumn/q1
      END IF
c
c     fourier-bessel transformation for spin-orbit density of
c     protons and neutrons given on field 'sjs'.
c
      if(ifrso.NE.0) then
        IF(q1==0D0) THEN
          sumsp  = zero
          sumsn  = zero
          r      = zero
          do i=2,ngrid
            r      = rstep+r
            sumsp  = (r*r)*sjsp(i)+sumsp
            sumsn  = (r*r)*sjsn(i)+sumsn
          enddo
          sumsp=frpirs*sumsp
          sumsn=frpirs*sumsn
        ELSE
          sumsp  = zero
          sumsn  = zero
          r      = zero
          c      = one
          s      = zero
          do i=2,ngrid
            r      = rstep+r
            x      = c*dc-s*ds
            s      = s*dc+c*ds
            c      = x
            sumsp  = (r*s)*sjsp(i)+sumsp
            sumsn  = (r*s)*sjsn(i)+sumsn
          enddo
          sumsp=frpirs*sumsp/q1
          sumsn=frpirs*sumsn/q1
        END IF
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
       IF(ifrso==1) THEN
         schsmp = schsm*xmgpr - schsep*d4m     ! yields mu+1/2 = old
         schsmn = schsm*xmgne - schsen*d4m
       ELSE 
         schsmp = schsm*xmgpr + schsep*d4m     ! yields mu-1/2 = correct
         schsmn = schsm*xmgne + schsen*d4m
       END IF
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
       if(ifrso.NE.0) then
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
      IF(ttestterm) THEN
!        sumc   = sump*(schsep-xmgpr*0.5D0*d4m*schsm)
!     &          + sumn*(schsen-xmgne*0.5D0*d4m*schsm)
!        schsmp = schsm*xmgpr + 2D0*schsep*d4m
!        schsmn = schsm*xmgne + 2D0*schsen*d4m
        sumc   =  sump*schsep + sumn*schsen
        schsmp = schsm*xmgpr - schsep*d4m
        schsmn = schsm*xmgne !- schsen*d4m
!        foc    = exphf*(sumc*darfac + (schsmp*sumsp+schsmn*sumsn))
        foc    = (sump + schsmp*sumsp)
        IF(q1==0D0) THEN
          sump0=sump
          sumn0=sumn
          sumsp0=sumsp
          sumsn0=sumsn
          schsmp0=schsmp
          schsmn0=schsmn
          darfac0=darfac
          sumsachsp0=schsmp*sumsp
          sumsachsn0=schsmn*sumsn
          WRITE(199,'(a,a)') 'q1,sump,sumn,sumsp,sumsn,',
     &      'schsmp,schsmn,darfac,schsmp*sumsp,schsmn*sumsn'
        ELSE
        END IF
        WRITE(199,'(10(1pg13.5))') q1,sump,sumn,sumsp,sumsn,
     &      schsmp,schsmn,darfac,schsmp*sumsp,schsmn*sumsn
        WRITE(199,'(13x,10(1pg13.5))') 
     &     -6D0*(sump-sump0)/qa/sump0,
     &     -6D0*(sumn-sumn0)/qa/sump0,
     &     -6D0*(sumsp-sumsp0)/qa/sump0,
     &     -6D0*(sumsn-sumsn0)/qa/sump0,
     &     -6D0*(schsmp-schsmp0)/qa/sump0,
     &     -6D0*(schsmn-schsmn0)/qa/sump0,
     &     -6D0*(darfac-darfac0)/qa/sump0,
     &     -6D0*(schsmp*sumsp-sumsachsp0)/qa/sump0,
     &     -6D0*(schsmn*sumsn-sumsachsn0)/qa/sump0
      ELSE
!       schsmp = schsm*xmgpr - schsep*d4m
!       schsmn = schsm*xmgne - schsen*d4m
        foc    = exphf*(sumc*darfac + (schsmp*sumsp+schsmn*sumsn))
      END IF
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
      parameter (kfrm=1999)        ! maximum dimension of formfactor fields
!      dimension formfa(kfrm)
      dimension formfa(*)
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
     *    sum    = formfa(i)*signum/(one-xfx2/DBLE((i-1)*(i-1))) + sum
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
         sum    = formfa(i)*signum/(one-xfx2/DBLE((i-1)*(i-1))) + sum
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
      parameter (kfrm=999)        ! maximum dimension of formfactor fields
!      dimension foin(kfrm)
      dimension foin(*)
c
      common /radpar/delq,rleng,iqmax           ! to communicate to subr.
c
c     some numerical parameters for termination of the search
c
      data itq0mx,endq0/20,1D-12/
      data itqbis,endqbis/50,1D-16/
      data  zero,one/0.0,1.0/
      LOGICAL, PARAMETER :: tbisect=.TRUE.
c
c---------------------------------------------------------------------
c
c     scan changing sign in 'foin'.
c      'iq1' is position before the zero and 'iq1' after the zero.
c
      iqmaxm = iqmax-1
      do 10 iscan=iqstrt,iqmaxm
      iq1 = iscan
      if(foin(iscan)*foin(iscan+1).le.zero) goto 19
   10 continue
   19 continue
!      iq1   = iscan
      iq2   = iq1+1
      if(.NOT.tbisect .AND. foin(iq1).le.zero) iq2   = iq1-1
      q1    = iq1*delq-delq
      q2    = iq2*delq-delq
      f1    = foin(iq1)
      f2    = foin(iq2)
c
c     branch to bisection
c
      IF(tbisect) THEN
        DO itq0=1,itqbis
          q3 = (q1+q2)/2D0
          f3     = fbint(q3,foin)
          IF(ABS(f3)<endqbis) THEN
            qzero = q3
            RETURN
          ELSE IF(f1*f3>0D0) THEN
            f1 = f3
            q1 = q3
          ELSE
            f2 = f3
            q2 = q3
          END IF
        END DO
        qzero = q3
      ELSE
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
      ENDIF
c
      return
      end


c-----fourb ------------ part of the hartree-fock package ----------------------
c
      subroutine fourb(formfin,rhoout,delq,rstep,ngrid,iqmax)
c
      IMPLICIT REAL(8) (a-h,o-z)
c
c     Input variables are:
c      fromfin = given formfactor
c      rhoout  = back transformed density
c      delq    = step size in q space
c      rstep   = spacing of the (equidistant) radial grid
c      ngrid   = number of radial grid points
c      iqmax   = maximal number of grid points in q space
c
      REAL(8),INTENT(IN) :: formfin(ngrid)
      REAL(8),INTENT(OUT) :: rhoout(ngrid)
      REAL(8),INTENT(IN) :: delq,rstep
      INTEGER,INTENT(IN) :: ngrid,iqmax

      REAL(8),PARAMETER :: pi=3.141592654D0
c
c
c----------------------------------------------------------------
c
c
c     charge density by fourier back-transformation
c
      dpi3   = 0.5D0/(pi*pi)
c                                                         r=0.0 separately
      ACC    = 0D0
      do iq=2,iqmax
        qa     = delq*iq-delq
        ACC    = qa*qa*formfin(iq)+ACC
      enddo
      rhoout(1) = abs(ACC*delq*dpi3)
c                                                         now all other r
      r      = 0D0
      do ir=2,ngrid
        r      = rstep+r
c
        dc     = cos(r*delq)
        ds     = sin(r*delq)
        c      = 1D0
        s      = 0D0
        ACC    = 0D0
        qa     = 0D0
        do iq=2,iqmax
          qa     = delq+qa
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
!          s = SIN(r*qa)
!          c = COS(qa*r)
          ACC    = qa*s*formfin(iq)+ACC
        enddo
        rhoout(ir) = (ACC*delq*dpi3/r)
c
      enddo

      RETURN
      END
c-----rms ------------ part of the hartree-fock package ----------------------
c
      REAL(8) FUNCTION rmsact(rhoin,rstep,ngrid)
c
      IMPLICIT REAL(8) (a-h,o-z)
c
c     Input variables are:
c      rhoin  =  density
c      rstep  = spacing of the (equidistant) radial grid
c      ngrid   = number of radial grid points
c
      REAL(8),INTENT(IN) :: rhoin(ngrid)
      REAL(8),INTENT(IN) :: rstep
      INTEGER,INTENT(IN) :: ngrid
c
c
c----------------------------------------------------------------
c
c
c     charge density by fourier back-transformation
c
      acc0 = 0D0
      acc2 = 0D0
      r    = 0D0
      do ir=2,ngrid
        r    = rstep+r
        r2   = r*r 
c
        acc0 = r2*rhoin(ir) + acc0
        acc2 = r2*r2*rhoin(ir) + acc2
      enddo
      rmsact = SQRT(acc2/acc0)

      RETURN
      END

c-----overlp------------ part of the hartree-fock package ----------------------
c
      function wfnorm(wfa)
c
c     norm of the wave function 'wf'
c
      include 'params.inc'
c
      dimension wfa(kg1)
c
      include 'grid.inc'
c
c-------------------------------------------------------------------------------
c
c
      wfovlp = wfa(1)*wfa(1)*0.5
      do i=2,ngrid
        wfovlp = wfa(i)*wfa(i)+wfovlp
      enddo
      wfnorm = wfovlp*rstep
c
      return
      end
c-----drpot ------------ part of the hartree-fock package ----------------------
c
      subroutine drpot(ifile)
c
      include 'params.inc'
c
c     prints densities and potentials on unit 'ifile'.
c
c     the printed fields and length of the arrays are trsnferred
c     via common's, see the include statements.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'result.inc'
c
c-------------------------------------------------------------------------------
c
c     potentials, only on standard write-unit '6'
c
      if(ifile.eq.6) then
       write(ifile,38) 'inverse mass',(bmass(i),i=1,ngrid)
       write(ifile,39)   (bmass(i),i=ngridp,ngrid2)
       write(ifile,38) '   potential',(   pot(i),i=1,ngrid)
       write(ifile,39)      (pot(i),i=ngridp,ngrid2)
       write(ifile,38) '    l*s-pot.',(    potls(i),i=1,ngrid)
       write(ifile,39)      ( potls(i),i=ngridp,ngrid2)
      endif
c
c     the proton and neutron density
c
      write(ifile,38) '     density',(rho(i),i=1,ngrid)
      write(ifile,39)   (rho(i),i=ngridp,ngrid2)
c
c     further densities, only on standard write-unit '6'
c
      if(itmax.le.1 .and. ifile.eq.6) then
       write(ifile,38) ' kin.-densi.',(tau(i),i=1,ngrid)
       write(ifile,39)   (tau(i),i=ngridp,ngrid2)
       write(ifile,38) ' spin-densi.',( sj(i),i=1,ngrid)
       write(ifile,39)   ( sj(i),i=ngridp,ngrid2)
       write(ifile,38) 'spins-densi.',(sjs(i),i=1,ngrid)
       write(ifile,39)   (sjs(i),i=ngridp,ngrid2)
       write(ifile,38) ' d/dr-densi.',(d1r(i),i=1,ngrid)
       write(ifile,39)   (d1r(i),i=ngridp,ngrid2)
       write(ifile,38) ' d/dr**2-de.',(d2r(i),i=1,ngrid)
       write(ifile,39)   (d2r(i),i=ngridp,ngrid2)
      endif
c
c     single particle properties
c
      write(ifile,40)
   40 format(1h0,4x,'single particle energies')
      write(ifile,45) (spenrg(i),i=1,nmaxp)
      write(ifile,39) (spenrg(i),i=nminn,nmaxto)
      write(ifile,41)
   41 format(1h0,4x,'single particle occupation weights')
      write(ifile,45) (weight(i),i=1,nmaxp)
      write(ifile,39) (weight(i),i=nminn,nmaxto)
c
      return
c
   38 format(1h0,/5x,a//5x,'protons',/(1x,10g13.5))
   39 format(/5x,'neutrons',/(1x,10g13.5))
   45 format (/5x,'protons',/(1x,10g13.5))
c
      end
c-----priwf ------------ part of the hartree-fock package ----------------------
c
      subroutine priwf(ifile)
c
      include 'params.inc'
c
c     printing the single particle wave functions on unit 'ifile'.
c
c     the printed information is transferred via common's.
c
      include 'wavfun.inc'
      include 'result.inc'
      include 'grid.inc'
c
      character*5 nucleo(2)
      data nucleo/' prot','neutr'/
c
c-------------------------------------------------------------------------------
c
c     loop over isospin
c
      write(ifile,'(//2a/a,f5.2,a)')
     &   ' the single particle wavefunctions in the form',
     &   ' "sqrt(4*pi)*r*wf(r)"',' and in steps of ',rstep,' [fm]:'
      do 40 iso=1,2
      isadd  = ngrid*nmaxp*(iso-1)
      write(ifile,'(/1x,a,a)') nucleo(iso),'ons:'
      if(iso.eq.1) then
       nend   = nmaxp
       isoadd = 0
      else
       nend   = nmaxn
       isoadd = nmaxp
      endif
c
c     loop over states in oscillator ordering
c
      do n=1,nend
        js     = jsp(n)
        ll     = js/2
        j2     = ((js+1)/2)*2-1
        ie     = n*ngrid+isadd
        ia     = ie-ngrid+1
c
        write(ifile,'(3(a,i2),a,f7.2,a,f6.3,60(/1x,5g13.4))')
     &  ' state: n=',n,',j=',j2,'/2,l=',ll,
     &  ',  s.p.en.=',spenrg(isoadd+n),
     &  ',  occup.=',weight(isoadd+n),(wf(i),i=ia,ie)
      enddo
c
   40 continue
c
      return
      end
c-----plot  ------------ part of the hartree-fock package ----------------------
c
      subroutine plot(ifile)
c
      include 'params.inc'
c
c     printing densities and single particle energies in columns
c     such that the printout can easily be read by many plotting
c     systems.
c     output is written on unit 'ifile'.
c     printed information is transferred via common's.
c
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'result.inc'
      include 'constn.inc'
      common /ctest/edirec(ksl2),efkincut,ifkincut   ! to communicate 'spkin'
c
c      data one/1.0/
      data ifexpr/0/              ! switch for computing and printing exch.
c
      INTEGER, ALLOCATABLE :: isort(:)
      LOGICAL,PARAMETER :: TSCANPH=.FALSE.
      LOGICAL,PARAMETER :: TFOLDENS=.FALSE.
!      CHARACTER(7),PARAMETER :: gform='1pg13.5',gform6='1pg14.6'
      CHARACTER(8),PARAMETER :: gform='1pg20.12',gform6='1pg20.12'
c
c------------------------------------------------------------------------------
c
c
c     optional c.m. unfolding
c
      IF(tfoldens) THEN
        CALL fourb(fop,aux1,delq,rstep,ngrid,iqmax)
        CALL fourb(fon,aux2,delq,rstep,ngrid,iqmax)
        aux3(1:ngrid)=aux1(1:ngrid)+aux2(1:ngrid)
        WRITE(20,'(a,2i4)') 'nprot,nneut=',nprot,nneut
        WRITE(20,'(a,4(1pg13.5))') 'rms_c,rms_p,rms_n,rms_tot=',
     &        rmsc,rmsact(aux1,rstep,ngrid),rmsact(aux2,rstep,ngrid)
     &        ,rmsact(aux3,rstep,ngrid)
        WRITE(20,'(a/2a)') 'densities with m. correction:',
     &     '  i      r          rho-coul        rho-prot   ',
     &     '    rho-neut       rho-mass'
        DO i=1,ngrid
          r=(i-1)*rstep
          WRITE(20,'(i5,f7.3,2x,6('//gform//'))')  i,r,rhochr(i),
     &             aux1(i),aux2(i),aux3(i)
        END DO
        WRITE(20,'(1x)')
!      ELSE
!        aux1(1:ngrid)=rho(1:ngrid)
!        aux2(1:ngrid)=rho(1+ngrid:ngrid+ngrid)
!        aux3(1:ngrid)=aux1(1:ngrid)+aux2(1:ngrid)
      END IF
c     short print on 23
c
c      write(23,'(2i4,1pg13.5)')   nprot,nneut,etotal
c
c     Coulomb exchange energy
c  
      if(ifexpr.eq.1) then
        IF(TADSIC) THEN
          excfac = -1.418033*acoul*(1D0-(1D0/nprot)**0.33333333333D0)
          write(ifile,'(a)') ' ADSIC invoked'
        ELSE
          excfac = -1.418033*acoul
        END IF
c        excfac = -1.418011569*acoul  ! e^2 =  1.4399784  !  1.418033 for 1.44
        sum   = 0.0
        ract  = 0.0
        do ir=2,ngrid
          ract   = rstep + ract
          sum    = (ract*ract)*rho(ir)**1.33333333 + sum
        enddo
        write(ifile,'(1x,a,1pg13.5)') 
     &   ' Coulomb exchange energy=',excfac*0.75D0*frpirs*sum
        write(ifile,'(1x,a,1pg13.5)') 
     &   '   Coulomb direct energy=',ecoul
      endif
c
c     some global infomation as a header
c
      IF(tfoldens) WRITE(ifile,'(a)') 'radii all with c.m. correction'
      ainfac=one/anumb**(one/3.0)
      armfac=ainfac/sqrt(0.6D0)
      write(ifile,'(14x,a,a,4(/a,5('//gform//')))')
     &   '     rms          rdms       surface        halo  ',
     &   '     central dep.',
     &   '   mass radii:',rms,rdms,surf,
     &    rms-sqrt(0.6*rdms*rdms+3.0*surf*surf),
     &    15.3*(rdms/rsec-1.0),
     &   ' charge radii:',rmsc,rdmsc,surfc,
     &    rmsc-sqrt(0.6*rdmsc*rdmsc+3.0*surfc*surfc),
     &    15.3*(rdmsc/rsecc-1.0),
     &   '              ',rmsc*armfac,rdmsc*ainfac
       write(ifile,'(2(a,5('//gform//')/),2(a,2('//gform//')))')
     &   ' proton radii:',rmsp,rdmsp,surfp,
     &    rmsp-sqrt(0.6*rdmsp*rdmsp+3.0*surfp*surfp),
!     &    15.3*(rdmsp/rsecp-1.0),
     &   rho(1),
     &   ' neutr. radii:',rmsn,rdmsn,surfn,
     &    rmsn-sqrt(0.6*rdmsn*rdmsn+3.0*surfn*surfn),
!     &    15.3*(rdmsn/rsecn-1.0),
     &   rho(1+ngrid),
     &   '              ',rmsn*armfac,rdmsn*ainfac,
     &   ' diffe. radii:',rmsn-rmsp,rdmsn-rdmsp
       write(ifile,'(a,5('//gform//')/)') 
     &   '  c.m. rms**2:',
     &   -0.375D0*3D0*(amfc(1)+amfc(2))/((nprot+nneut)*zpe)
!     &-6D0*cmwid
c      write(23,'(2i4,1x,f9.2,f8.2,f7.2,3f6.2,4f7.3)')
c     &  nprot,nneut,etotal,epair,efaddn,deltn,deltn2,eqminn,
c     &  rmsn,rdmsn,surfn,
c     &  rmsn-sqrt(0.6*rdmsn*rdmsn+3.0*surfn*surfn)
c      write(23,'(2i4,2x,f10.2,3(2x,4(1pg13.5)))')
c     &  nprot,nneut,etotal,
c     &  rmsc*ainfac,rdmsc*ainfac,surfc,
c     &  rmsc-sqrt(0.6*rdmsc*rdmsc+3.0*surfc*surfc),
c     &  rmsn*ainfac,rdmsn*ainfac,surfn,
c     &  rmsn-sqrt(0.6*rdmsn*rdmsn+3.0*surfn*surfn),
c     &  rmsn-rmsp,rdmsn-rdmsp
!      write(ifile,'(a,2i4,1pg20.12,1pg13.5,0pf18.6,4f8.3/)')
!     &     ' Z,N,r**2,E,eferm,gap=',
!     &   nprot,nneut,rmsc*rmsc,
!     &   rmsc*rmsc-0.6*1.44*anumb**((one+one)/3.0),
!     &   etotal,efep,efen,deltp,deltn
      write(ifile,'(a,2i4,f10.6,1pg14.6,6(0pf8.3)/)')
     &     ' Z,N,r_C,E,epair,gap,sumuv=',
     &   nprot,nneut,rmsc,
     &   etotal,epp,epn,deltp,deltn,-epp/deltp,-epn/deltn
      OPEN(26,file='for026',position='APPEND')
      WRITE(26,'(2i4,2(1pg13.5))') nprot,nneut,deltp,deltn
      CLOSE(26)
!     &   rmsc*rmsc-0.6*1.44*anumb**((one+one)/3.0),
!     &   etotal,efep,efen
!      write(ifile,'(2i5,5f9.4,a)')
!     &   nprot,nneut,rms,rmsn-rmsp,rmsp,rmsn,rmsc,' RMS-RADII'
      IF(isoodd==2) THEN
         nblock=nblockn-nmaxp
        WRITE(2,'(a,5i5)') 'test: nprot,nneut,nmaxp,nblockn,nblock',
     &    nprot,nneut,nmaxp,nblockn,nblock
        write(ifile,'(2i4,4(1pg20.12),i3,a)')
     &   nprot,nneut,etotal,rmsc,rdmsc,surfc,jsp(nblock),'/2 BASIC'
      ELSE
        write(ifile,'(2i4,4(1pg20.12),a)')
     &   nprot,nneut,etotal,rmsc,rdmsc,surfc,' BASIC'
      END IF
c      if(nneut.eq.126 .and. nprot.eq.82)
c     &  write(ifile,'(a,1pg12.4)') 
c     &   ' shell oscillations in charge density:',
c     &   rhochr(1)-rhochr(7)
      write(ifile,'(a,1(1pg16.9),2(a,1pg12.5)/a,1pg12.5/
     &              2(a,1pg12.4),
     &              2(/a,1pg12.5,a,2(1pg12.5)))')
     &       '  total energy=',etotal,
     &               '   zpe=',zpe,'  epair=',epair,
     &       '           e/a=',etotal/anumb,
     &       '  without ZPE: E=',etotal+zpe,'  E/A=',(etotal+zpe)/anumb,
     &       '    total ekin=',ekin,'     ekin_p/n=',ekinpr,ekinne,
     &       '        ekin/A=',ekin/anumb,
     &       ' rel. ekin_p/n=',ekinpr/aprot,ekinne/aneut
      write(ifile,'(a,3(1pg16.9))')
     &        '  shell energies=',eshl(1),eshl(2),eshl(1)+eshl(2)
      sqa = sqrt(anumb)
      write(ifile,'(2(a,2(1pg12.4)))')
     &  ' pairing:  delta_p/n',deltp*sqa,deltn*sqa,
     &           '  gapfc_p/n',gpp*anumb,gpn*anumb,
     &  '          efermi_p/n',efep,efen,
     &           '  epair_p/n',epp,epn,
     &  '           delta_u*v',deltp,deltn,
     &           '  delta_v*v',deltp2,deltn2,
     &  '          equasi_min',eqminp,eqminn
c
c     print densities in plot format
c
      if(ifplot.ge.2) then
      IF(tfoldens) WRITE(ifile,'(a)') 
     &   'densities all with c.m. correction'
       write(ifile,'(a/2a)')  ' plot the densities:',
     &     '  i      r          rho-coul        rho-prot   ',
     &     '      rho-neut       rho-mass'
       if(ipair.lt.5) then
         do i=1,ngrid
           r=(i-1)*rstep
           write(ifile,'(i5,f7.3,2x,4('//gform//'))')  i,r,rhochr(i),
     &             rho(i),rho(i+ngrid),(rho(i)+rho(i+ngrid))
         enddo
       else
         acc=0D0
         do i=1,ngrid
           r=(i-1)*rstep
           write(ifile,'(i5,f7.3,2x,8('//gform//'))')  i,r,rhochr(i),
     &             rho(i),rho(i+ngrid),(rho(i)+rho(i+ngrid)),
     &             sjs(i),sjs(i+ngrid),
     &             (d1r(i)+d1r(i+ngrid))**2,
     &             paircr(i)**2+paircr(i+ngrid)**2
           acc=acc+(r*r)*(d1r(i)+d1r(i+ngrid))**2
         enddo
         WRITE(11,'(a,1pg13.5)') ' integrated (grad rho)**2:',
     &      4D0*pi*rstep*acc
       endif
       write(ifile,'(a/a)')  ' plot the formfactor:',
     &     '  i     q       f-coul   f-prot    f-neut   f-total'//
     &     '    f-weak'
       delqh = 0.5*delq
       write(ifile,'(i5,f7.3,5('//gform6//'))')  
     &     1,zero,aprot,aprot,aneut,aprot+aneut,fow(1)
       do iq=2,iqmax
         q1     = iq*delqh-delqh
         write(ifile,'(i5,f7.3,5('//gform6//'))') iq,q1,fbint(q1,formf),
     &     fbint(q1,fop),fbint(q1,fon),fbint(q1,fo),fbint(q1,fow)
       enddo
       write(ifile,'(a/2a)')  ' plot the potentials:',
     &     '  i     r        pot-prot    ',
     &     '    m*/m_p     pot-neut      m*/m_n   pair.pot_p/n'
       potmxp = zero
       potmxn = zero
       DO i=1,ngrid
         r=(i-1)*rstep
         write(ifile,'(i5,f7.3,2x,6g14.5)')
     &      i,r,pot(i),amfc(1)/bmass(i),pot(i+ngrid),
     &      amfc(2)/bmass(i+ngrid),pairpt(i),pairpt(i+ngrid)
         potmxp = max(pot(i),potmxp)
         potmxn = max(pot(i+ngrid),potmxn)
       END DO
       write(ifile,'(a/2a)')  ' plot the l*s potentials:',
     &     '  i     r        wls-prot    ',
     &     '   wls-neut      div*J-prot   div*J-neut'
       DO i=1,ngrid
         r=(i-1)*rstep
         write(ifile,'(i5,f7.3,2x,6g14.5)')
     &      i,r,potls(i),potls(i+ngrid),sjs(i),sjs(i+ngrid)
       END DO
      endif
c
c     print single particle energies and occupation weights
c
      if(ifplot.ge.3) then
c       multmx = 999                            ! set large to get multiplicity
       multmx = 1                              ! set 1 to have shorter print
       write(ifile,'(/a/2a)')
     &  '  proton spectrum:',
     &  '  2j    l    energy       weight    e-variance  rms    ',
     &  '  multipl.   dyn.gap  wgtherm s.p.kinetic'
       ALLOCATE(isort(nmaxp))
       DO n=1,nmaxp
         isort(n) = n
       END DO
       DO n=1,nmaxp
         emin=1D3
         imin=n
         DO n2=n,nmaxp
           nact=isort(n2)
           IF(spenrg(nact)<emin) THEN
             emin=spenrg(nact)
             isn=isort(n)
             isort(n)=isort(n2)
             isort(n2)=isn
           END IF
         END DO
       END DO
       sumuvp = zero
       sumcon = zero
       sencon = zero
       do 20 n2=1,nmaxp
       n=isort(n2)
       js=jsp(n)
       ll=js/2
       xlp    = float(ll*ll+ll+1)
       if(2*ll.eq.js) js=js-1
       IF(tphpair) THEN
         equasi = sqrt((spenrg(n)-efep)**2+phpair(n)*deltaf(n)**2)
       ELSE
         equasi = sqrt((spenrg(n)-efep)**2+deltaf(n)**2)
       END IF
       npl = (n-1)*ngrid
       sum0 = 0.0
       sum2 = 0.0
       sumk = 0D0
       sumk0 = 0D0
       do ir=2,ngrid
         r  = (ir-1)*rstep
         r2 = r*r
         vol= wf(npl+ir)*wf(npl+ir)  ! *r2
         sum0 = vol + sum0
         sum2 = vol*r2 + sum2
!         sumk = wfs(npl+ir)*wfs(npl+ir)*(rho(ir)+rho(ir+ngrid)) + sumk
!         sumk0 = wf(npl+ir)*wf(npl+ir)*(rho(ir)+rho(ir+ngrid)) + sumk0
         wfas=wfs(npl+ir)
         wfa=wf(npl+ir)
         ri=1D0/r
         sumk0 = sumk0 + wfa*wfa*(rho(ir)+rho(ir+ngrid))**3
         sumk = (wfas*wfas - ((wfa*wfas)+(wfa*wfas)-xlp*wfa*wfa*ri)*ri)
     &          *(rho(ir)+rho(ir+ngrid))**3 + sumk
       enddo
       sum2 = sqrt(abs(sum2/max(sum0,xsmall)))
       sumk = sumk/sumk0
       sumuvp = SQRT(weight(n)-weight(n)**2)*(js+1)/2 + sumuvp
       rmssp(n) = sum2
       if(spenrg(n).gt.zero) then
         vol    = weight(n)*(js+1)
         sumcon = vol           + sumcon
         sencon = vol*spenrg(n) + sencon
       endif
       mmult = min((js+1)/2,multmx)                   ! half the multiplicity
       do imult=1,mmult
         spact  = spenrg(n)
         wgact  = weight(n)
         if(spact.ne.zero) then
           if(ipair.lt.5) then
             write(ifile,'(1x,i3,i5,4g13.5,i5)') js,ll,
     &         spenrg(n),weight(n),hpswid(n),sum2,js+1
           else
             if(temprt.eq.zero) then
c               write(ifile,'(1x,i3,i5,4g13.5,i5,2g12.4)') js,ll,
c     &           spenrg(n),weight(n),hpswid(n),sum2,js+1,
c     &           deltaf(n),phpair(n)
               write(ifile,'(1x,i3,i5,4g13.5,i5,4g12.4)') js,ll,
     &           spenrg(n),weight(n),hpswid(n),sum2,js+1,
     &           deltaf(n),phpair(n),spkin(n),edirec(n) ! sumk
!     &           deltaf(n),equasi
             else
               write(ifile,'(1x,i3,i5,4g13.5,i5,2g12.4)') js,ll,
     &           spenrg(n),weight(n),hpswid(n),sum2,js+1,
     &           deltaf(n),ttthrm(n)
             endif
           endif
         endif
       enddo
   20  continue
       if(multmx.gt.1 .and. ifile.eq.11) then
         backspace 11
         write(11,'(18x,i3,i5,3g13.5,i5)') -js,ll,
     &       spact,wgact,sum2,js+1
       endif
!       if(sumcon.gt.zero) then
         write(ifile,'(/2(a,1pg12.3)/a,1pg12.3)')
     &   ' nr. protons in continuum=',sumcon,
     &   ' average E in continuum=',sencon/MAX(sumcon,1D-99),
     &   ' maximum potential:',potmxp
!       endif
!        WRITE(17,'(2i4,86(1pg11.3))')
!     &   nprot,nneut,(rmssp(i),i=1,nmaxp)
!        CALL flush(17)
        WRITE(16,'(2i4,86(1pg11.3))')
     &   nprot,nneut,(SQRT(v2occ(i)-v2occ(i)**2),i=1,nmaxp)
        CALL flush(16)
       DEALLOCATE(isort)
c
       ALLOCATE(isort(nmaxn))
       DO n=1,nmaxn
         isort(n) = n+nmaxp
       END DO
       DO n=1,nmaxn
         nn=n+nmaxp
         emin=1D3
         imin=n
         DO n2=n,nmaxn
           nact=isort(n2)
           IF(spenrg(nact)<emin) THEN
             emin=spenrg(nact)
             isn=isort(n)
             isort(n)=isort(n2)
             isort(n2)=isn
           END IF
         END DO
       END DO
       write(ifile,'(/a/2a)')
     &  '  neutron spectrum:',
     &  '  2j    l    energy       weight    e-variance  rms    ',
     &  '  multipl.   dyn.gap  wgtherm   s.p.kinetic'
       sumcon = zero
       sencon = zero
       sumuvn=zero
       do 30 nn=1,nmaxn
       n=isort(nn)  !       n=nn+nmaxp
       js=jsp(n-nmaxp)
       ll=js/2
       xlp    = float(ll*ll+ll+1)
       if(2*ll.eq.js) js=js-1
       npl = (n-1)*ngrid
       sum0 = 0D0
       sum2 = 0D0
       sum4 = 0D0
       sumk = 0D0
       sumk0 = 0D0
       do ir=2,ngrid
         r  = (ir-1)*rstep
         r2 = r*r
         vol= wf(npl+ir)*wf(npl+ir)  !*r2
         sum0 = vol + sum0
         sum2 = vol*r2 + sum2
         sum4 = vol*r2*r2 + sum4
!         sumk = wfs(npl+ir)*wfs(npl+ir)*(rho(ir)+rho(ir+ngrid)) + sumk
!         sumk0 = wf(npl+ir)*wf(npl+ir)*(rho(ir)+rho(ir+ngrid)) + sumk0
         wfas=wfs(npl+ir)
         wfa=wf(npl+ir)
         ri=1D0/r
         sumk0 = sumk0 + wfa*wfa*(rho(ir)+rho(ir+ngrid))**3
         sumk = (wfas*wfas - ((wfa*wfas)+(wfa*wfas)-xlp*wfa*wfa*ri)*ri)
     &          *(rho(ir)+rho(ir+ngrid))**3 + sumk
!         sumk0 = sumk0 + wfa*wfa
!         sumk = (wfas*wfas - ((wfa*wfas)+(wfa*wfas)-xlp*wfa*wfa*ri)*ri)
!     &          + sumk
       enddo
       sum2 = sqrt(abs(sum2/max(sum0,xsmall)))
       sum4 = (abs(sum4/max(sum0,xsmall)))
       sum4 = (sum4-sum2**4)**0.25D0
       sumk = sumk/sumk0*amfc(2)
       IF(nn>10 .AND. nn<17) 
     &   sumuvn = SQRT(weight(n)-weight(n)**2)*(js+1)/2 + sumuvn
       if(spenrg(n).gt.zero) then
         vol    = weight(n)*(js+1)
         sumcon = vol           + sumcon
         sencon = vol*spenrg(n) + sencon
       endif
       if(2*ll.eq.js) js=js-1
       mmult = min((js+1)/2,multmx)                   ! half the multiplicity
       do imult=1,mmult
         spact  = spenrg(n)
         wgact  = weight(n)
         if(spact.ne.zero) then
           if(ipair.lt.5) then
             write(ifile,'(1x,i3,i5,4g13.5,i5)') js,ll,
     &         spenrg(n),weight(n),hpswid(n),sum2,js+1
           else
             if(temprt.eq.zero) then
c               write(ifile,'(1x,i3,i5,4g13.5,i5,2g12.4)') js,ll,
c     &           spenrg(n),weight(n),hpswid(n),sum2,js+1,
c     &           deltaf(n),phpair(n)
               write(ifile,'(1x,i3,i5,4g13.5,i5,8g12.4)') js,ll,
     &           spenrg(n),weight(n),hpswid(n),sum2,js+1,
     &           deltaf(n),phpair(n),spkin(n),edirec(n),sum4
!     &           1D0/(1D0+EXP((edirec(n)+cutfac)/(-0.05D0*cutfac))) !sum4 !sumk
             else
               write(ifile,'(1x,i3,i5,3g13.5,f6.3,i3,4g12.4)') 
     &           js,ll,
     &           spenrg(n),weight(n),hpswid(n),sum2,js+1,
     &           deltaf(n),ttthrm(n)
             endif
           endif
         endif
       enddo
   30  continue
       if(multmx.gt.1 .and. ifile.eq.11) then
	 backspace 11
         write(11,'(18x,i3,i5,3g13.5,i5)') -js,ll,
     &       spact,wgact,sum2,js+1
       endif
!       WRITE(11,*) 'sumuv=',sumuv
!       if(sumcon.gt.zero) then
         write(ifile,'(/2(a,1pg12.3)/a,1pg12.3)')
     &   ' nr. neutrons in continuum=',sumcon,
     &   ' average E in continuum=',sencon/MAX(sumcon,1D-99),
     &   ' maximum potential:',potmxn
!       endif
      endif
       DEALLOCATE(isort)
c
c     quick hack to print s.p. energies in columns (for Pb isotopes)
c
c      write(32,'(2i4,11f7.2)')
c     &   nprot,nneut,(spenrg(n),n=11,19)
c      write(33,'(2i4,11f7.2)')
c     &   nprot,nneut,(spenrg(n+nmaxp),n=12,22)
c
c     evaluate the average effective mass in the nucleus
c
      if(ifplot.ge.4) then
c                                            'rho' integral over solution
        bmaspr = zero
        bmasne = zero
        do i=1,ngrid
	  vol = rsquar(i)
	  bmaspr = vol*bmass(i)*rho(i) + bmaspr
	  bmasne = vol*bmass(i+ngrid)*rho(i+ngrid) + bmasne
        enddo
        bmaspr = frpirs*bmaspr/(aprot*amfc(1))
        bmasne = frpirs*bmasne/(aneut*amfc(2))
        bmasav = (aprot*bmaspr+aneut*bmasne)/anumb
c                                            'tau' integral over solution
        ataupr = zero
        ataune = zero
        tmaspr = zero
        tmasne = zero
        do i=1,ngrid
	  vol = rsquar(i)
	  ataupr = vol*tau(i) + ataupr
	  ataune = vol*tau(i+ngrid) + ataune
	  tmaspr = vol*bmass(i)*tau(i) + tmaspr
	  tmasne = vol*bmass(i+ngrid)*tau(i+ngrid) + tmasne
        enddo
        tmaspr = tmaspr/(amfc(1))
        tmasne = tmasne/(amfc(2))
        tmasav = (tmaspr+tmasne)/(ataupr+ataune)
        tmaspr = tmaspr/ataupr
        tmasne = tmasne/ataune
c                                            Saxon-Woods estimate
        call nucmat(rhonm,eovera,compr,effmas,f0,f0p,g0,g0p,
     &              asymme,dasymm,d2asym,skappa,0,-2)
	etapr  = 0.58*surfp/rdmsp
	eta2   = (etapr*pi)**2
	bswpr  = one+
     &           (one/effmas-one)*(one-etapr*(3+eta2)/(1+eta2))
	etane  = 0.58*surfn/rdmsn
	eta2   = (etane*pi)**2
	bswne  = one+
     &           (one/effmas-one)*(one-etane*(3+eta2)/(1+eta2))
        bswav = (aprot*bswpr+aneut*bswne)/anumb
c                                            print all
	write(11,'(/a/10x,3(8x,a)/4a/10f8.3)')
     &    ' effective masses:',
     &    'estimate','rho-average','tau-average',
     &    '  m*/m  ',('proton neutron total ',i=1,3),
     &    effmas,
     &     one/bswpr,one/bswne,one/bswav,
     &     one/bmaspr,one/bmasne,one/bmasav,
     &     one/tmaspr,one/tmasne,one/tmasav
c
c       print robust(?) combination of densities and/or kinetic densities
c
        write(11,'(/a/a)')
     &   ' protons:',
     &   ' ir   r   rho*tau-d1r**2/4 ...'
        do ir=1,ngrid
          write(11,'(1x,i5,f6.2,6(1pg12.4))')
     &    ir,(ir-1)*rstep,
     &    rho(ir)*tau(ir)-0.25*d1r(ir)**2,
     &    rho(ir)*tau(ir),d1r(ir)**2,
     &    rho(ir),tau(ir),d1r(ir)
        enddo
        write(11,'(/a/a)')
     &   ' neutrons:',
     &   ' ir   r   d1r**2-rho*d2r  rho*tau-d1r**2/4 ...'
        do ir=ngrid+1,ngrid+ngrid
          write(11,'(1x,i5,f6.2,6(1pg12.4))')
     &    ir,(ir-ngrid-1)*rstep,
     &    rho(ir)*tau(ir)-0.25*d1r(ir)**2,
     &    rho(ir)*tau(ir),d1r(ir)**2,
     &    rho(ir),tau(ir),d1r(ir)
        enddo
      endif
c
      if(TSCANPH .AND. ifplot.ge.3 .AND. nprot.EQ.82
     &    .AND. nneut.EQ.126) then
        OPEN(36,file='scanph.res')
        call scanph(0D0,30D0,12,-1)
        call scanph(0D0,30D0,14,-1)
        CLOSE(36)
      endif      
c
      return
c
      end
c-----prispe------------ part of the hartree-fock package ----------------------
c
      subroutine prispe
c
      include 'params.inc'
      parameter (ksort=kmaxn+1)
c
c     printing header and single particle energies on file 'for027'
c
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'result.inc'
      include 'constn.inc'
c
      dimension esort(ksort),isort(ksort)
      data iop27/0/
      logical tmult
      parameter (tmult=.true.)
      data mmult/1/               ! default
c
c------------------------------------------------------------------------------
c
      if(iop27.eq.0) then
        open(unit=27,file='for027')
        iop27 = 1
      endif
c
c     some global infomation as a header
c
      ifile = 27
      write(27,'(35(a2))') (' *',i=1,35)
      write(27,'(a)') title
       write(27,'(/4(a,i5)/a,3f10.3,f10.1,f10.3/a,4f10.5,a,f8.5)')
     &   ' nucleus:  nprot=', nprot,'  nneut=',nneut,
     &            '  nmaxp=', nmaxp,'  nmaxn=',nmaxn,
     &   ' force:    t0-t4=',t0,t1,t2,t3,t4,
     &   '           x0-x3=',x0,x1,x2,x3,'  power=',power
c       if(tyuk.ne.zero) write(11,'(12x,a,f10.3,a,f10.5,a,f8.5)')
c     &               'tyuk=',tyuk,'   xyuk=',xyuk,'   ayuk=',ayuk
       if(abs(t5)+abs(t6)+abs(t7)+abs(t8).gt.zero)
     &  write(27,'(10x,a,4f10.3/10x,a,4f10.5)')
     &             ' t5-t8=',t5,t6,t7,t8,
     &             ' x5-x8=',x5,x6,x7,x8
       write(27,
     & '(10x,a,i3,a,2f8.3,a,f8.3,a,i3,/10x,a,f6.3/5(a,i5)/10x,5(a,i5))')
     &             ' ipair=',ipair,' pairfp/n=',pairfp,pairfn,
     &             ' temprt=',temprt,
     &             ' ilnrsp=',ilnrsp,' acoul=',acoul,
     &   ' cases:     ifex=',ifex,'  iftm=',iftm,'  ifzpe=',ifzpe,
     &            ' ifdfms=',ifdfms,' ifrhoc=',ifrhoc,
     &             ' ifrso=',ifrso,' ifgsc=',ifgsc,' ifenrg=',ifenrg,
     &            ' ifplot=',ifplot,' imfold=',imfold
       write(27,'(a,f8.3,a,i5/a,f6.3,a,f6.2,2(a,1pg9.2)/10x,2(a,i5))')
     &   ' grid:     rstep=',rstep,' ngrid=',ngrid,
     &   ' step:    x0damp=',x0damp,' e0damp=',e0damp,' dabr=',dabr,
     &              ' dabe=',dabe,
     &              ' itmax=',itmax,' itpr=',itpr
      ainfac=one/anumb**(one/3.0)
      write(ifile,'(14x,a,4(/a,4(1pg13.5)))')
     &   '     rms          rdms       surface        halo',
     &   '   mass radii:',rms,rdms,surf,
     &    rms-sqrt(0.6*rdms*rdms+3.0*surf*surf),
     &   ' charge radii:',rmsc,rdmsc,surfc,
     &    rmsc-sqrt(0.6*rdmsc*rdmsc+3.0*surfc*surfc),
     &   '              ',rmsc*ainfac,rdmsc*ainfac
       write(ifile,'(2(a,4(1pg13.5)/),2(a,2(1pg13.5)/))')
     &   ' proton radii:',rmsp,rdmsp,surfp,
     &    rmsp-sqrt(0.6*rdmsp*rdmsp+3.0*surfp*surfp),
     &   ' neutr. radii:',rmsn,rdmsn,surfn,
     &    rmsn-sqrt(0.6*rdmsn*rdmsn+3.0*surfn*surfn),
     &   '              ',rmsn*ainfac,rdmsn*ainfac,
     &   ' diffe. radii:',rmsn-rmsp,rdmsn-rdmsp
      write(ifile,'(a,1(1pg16.9),2(a,1pg12.5)/a,1pg12.5/
     &              2(a,1pg12.4),
     &              2(/a,1pg12.5,a,2(1pg12.5)))')
     &       '  total energy=',etotal,
     &               '   zpe=',zpe,'  epair=',epair,
     &       '           e/a=',etotal/anumb,
     &       '  without ZPE: E=',etotal+zpe,'  E/A=',(etotal+zpe)/anumb,
     &       '    total ekin=',ekin,'     ekin_p/n=',ekinpr,ekinne,
     &       '        ekin/A=',ekin/anumb,
     &       ' rel. ekin_p/n=',ekinpr/aprot,ekinne/aneut
      write(ifile,'(a,3(1pg16.9))')
     &        '  shell energies=',eshl(1),eshl(2),eshl(1)+eshl(2)
      sqa = sqrt(anumb)
      write(ifile,'(2(a,2(1pg12.4)))')
     &  ' pairing:  delta_p/n',deltp*sqa,deltn*sqa,
     &           '  gapfc_p/n',gpp*anumb,gpn*anumb,
     &  '          efermi_p/n',efep,efen,
     &           '  epair_p/n',epp,epn
c
c     print single particle energies and occupation weights for protons
c
c     sort energies
c
      do n=1,nmaxp
        esort(n) = spenrg(n)
        if(esort(n).eq.zero) esort(n) = 99999.0
        isort(n) = n
      enddo
      do n=1,nmaxp
        emin = esort(n)
        nmin = n
        do n2=n,nmaxp
          if(esort(n2).lt.emin) then
            emin     = esort(n2)
            esort(n2)= esort(n)
            esort(n) = emin
            imin     = isort(n2)
            isort(n2)= isort(n)
            isort(n) = imin
c            nmin        = n2
          endif
        enddo
      enddo
      do n=1,nmaxp
        nsort(n) = isort(n)                 ! save for 'prispd'
      enddo
c      write(27,'(1x,10i5)') (isort(n),n=1,nmaxp)
c      write(27,'(1x,5g12.4)') (esort(n),n=1,nmaxp)
c
c     print
c
      write(ifile,'(/a/2a)')
     & '  proton spectrum:',
     & '  2j    l    energy       weight   rms     multipl.   '
      do n2=1,nmaxp
        n  = isort(n2)
        js=jsp(n)
        ll=js/2
        if(2*ll.eq.js) js=js-1
        npl = (n-1)*ngrid
        sum0 = 0.0
        sum2 = 0.0
        do ir=1,ngrid
          r  = (ir-1)*rstep
          r2 = r*r
          vol= wf(npl+ir)*wf(npl+ir)*r2
          sum0 = vol + sum0
          sum2 = vol*r2 + sum2
        enddo
        sum2 = sqrt(abs(sum2/max(sum0,xsmall)))
        if(tmult) mmult = (js+1)/2                   ! half the multiplicity
        spact  = spenrg(n)
        wgact  = weight(n)
        if(spact.ne.zero .or. sum2.ne.zero) then
          do imult=1,mmult
            write(ifile,'(1x,i3,i5,3g15.7,i5)') js,ll,
     &         spact,wgact,sum2,js+1
          enddo
        endif
      enddo
      backspace 27
      read(27,*) js,ll,spact,wgact,sum2
      backspace 27
      write(27,'(1x,i4,i5,3g15.7,i5)') -js,ll,
     &         spact,wgact,sum2,js+1
c
c
c     and now printout for the neutrons
c
c     sort energies
c
      do n=1,nmaxn
        esort(n) = spenrg(n+nmaxp)
        if(esort(n).eq.zero) esort(n) = 99999.0
        isort(n) = n
      enddo
      do n=1,nmaxn
        emin = esort(n)
        nmin = n
        do n2=n,nmaxn
          if(esort(n2).lt.emin) then
            emin     = esort(n2)
            esort(n2)= esort(n)
            esort(n) = emin
            imin     = isort(n2)
            isort(n2)= isort(n)
            isort(n) = imin
c            nmin        = n2
          endif
        enddo
      enddo
      do n=1,nmaxn
        nsort(n+nmaxp) = isort(n)+nmaxp             ! save for 'prispd'
      enddo
c
c     print
c
      write(ifile,'(/a/2a)')
     &  '  neutron spectrum:',
     &  '  2j    l    energy       weight   rms     multipl.   '
      do nn=1,nmaxn
        n=isort(nn)+nmaxp
        js=jsp(n-nmaxp)
        ll=js/2
        if(2*ll.eq.js) js=js-1
        npl = (n-1)*ngrid
        sum0 = 0.0
        sum2 = 0.0
        do ir=1,ngrid
          r  = (ir-1)*rstep
          r2 = r*r
          vol= wf(npl+ir)*wf(npl+ir)*r2
          sum0 = vol + sum0
          sum2 = vol*r2 + sum2
        enddo
        sum2 = sqrt(abs(sum2/max(sum0,xsmall)))
c??        if(2*ll.eq.js) js=js-1
        if(tmult) mmult = (js+1)/2                   ! half the multiplicity
        spact  = spenrg(n)
        wgact  = weight(n)
        if(spact.ne.zero .or. sum2.ne.zero) then
          do imult=1,mmult
            write(ifile,'(1x,i3,i5,3g15.7,i5)') js,ll,
     &         spact,wgact,sum2,js+1
          enddo
        endif
      enddo
      backspace 27
      read(27,*) js,ll,spact,wgact,sum2
      backspace 27
      write(27,'(1x,i4,i5,3g15.7,i5)') 
     &         -js,ll,spact,wgact,sum2,js+1
c
      return
c
      end
c-----barrad------------ part of the hartree-fock package ----------------------
c
      subroutine barrad(nprot,nneut)
c
      include 'params.inc'
c
c     computes the barret-radius 'rbarr' which is defined by
c      3/rbarr^3*int(r^xk*exp(-alpha*r)),r=0..rbarr) = <r^k*exp(-alpha*r)>.
c
c     input via list is:
c      nprot  = proton number
c     input via common /result/ is:
c      rhochr = array(1..ngrid) carrying charge density
c               (the point r=0 is the first grid point)
c      rdmsc  = diffraction radius of the charge distribution
c               (but an estimate 'rdmsc' = 1.2*a**1/3 would do as well)
c     input via common /grid/ is:
c      ngrid  = number of grid points on which 'rhochr' is given.
c      rstep  = step in r-space of the grid for 'rhochr'
c      frpirs = 4*pi*rstep
c
      include 'grid.inc'
      include 'result.inc'
      include 'constn.inc'
c
c     table of parameters for the barrett-moment
c
      dimension ncharg(0:9),xkapbr(0:9),alphbr(0:9)
      data ncharg/
     &   20,  40,   42,   50,   52,   56,   62,   64,   80,   82/
      data xkapbr/
     &  2.0914, 2.150,2.155,2.187,2.195,2.211,2.338,2.304,2.320,2.333/
      data alphbr/
     &  0.0595, 0.093,0.096,0.108,0.110,0.116,0.139,0.135,0.150,0.153/
c
c     accuracy of the root-finding and max.number of iterations
c
      data accura/0.001/,maxit/50/
c
c     internal statement function for
c     3/rbarr^3*int(r^xk*exp(-alpha*r)),r=0..rbarr)
c
c      barint(rbarr,xk,alpha) =   3.0*rbarr**xk*          ! old maple result
c     &  (.9999999920/(xk+3.0)
c     &  -.9999998944*(rbarr*alpha)**(1)/(xk+4.0)
c     &  +.4999993609*(rbarr*alpha)**(2)/(xk+5.0)
c     &  -.1666643685*(rbarr*alpha)**(3)/(xk+6.0)
c     & + .04166122691*(rbarr*alpha)**(4)/(xk+7.0)
c     &  -.008324461018*(rbarr*alpha)**(5)/(xk+8.0)
c     & + .001378760880 *(rbarr*alpha)**(6)/(xk+9.0)
c     & - .0001903759047*(rbarr*alpha)**(7)/(xk+10.0)
c     &  + .00002050202051*(rbarr*alpha)**(8)/(xk+11.0)
c     & - .000001301715588*(rbarr*alpha)**(9)/(xk+12.0) )
      barint(rbarr,xk,alpha) =   3.0*rbarr**xk*
     &  (1.0/(xk+3.0)
     &  -(rbarr*alpha)**(1)/(xk+4.0)
     &  +(rbarr*alpha)**(2)/(2.0*(xk+5.0))
     &  -(rbarr*alpha)**(3)/(6.0*(xk+6.0))
     &  +(rbarr*alpha)**(4)/(24.0*(xk+7.0))
     &  -(rbarr*alpha)**(5)/(120.0*(xk+8.0))
     &  +(rbarr*alpha)**(6)/(720.0*(xk+9.0))
     &  -(rbarr*alpha)**(7)/(5040.0*(xk+10.0))
     &  +(rbarr*alpha)**(8)/(40320.0*(xk+11.0))
     &  -(rbarr*alpha)**(9)/(362880.0*(xk+12.0)) )
c
c-------------------------------------------------------------------------------
c
c     search table of barrett-parameters
c
      do itry=0,9
        ibar = itry
        if(ncharg(itry).eq.nprot) goto 19
      enddo
      write(11,'(a)') ' no parameters for the barrett-moment.'
      return
c
   19 continue
      xkact  = xkapbr(ibar)
      alphac = alphbr(ibar)
c
c     compute barrett-moment
c
      sum   = zero
      sum2  = zero
      sum4  = zero
      sum6  = zero
      do ir=1,ngrid
        r     = (ir-1)*rstep
        sum   = r**(xkact+2.0)*exp(-alphac*r)*rhochr(ir) + sum
        sum2  = r**(4.0)*rhochr(ir) + sum2
        sum4  = r**(6.0)*rhochr(ir) + sum4
        sum6  = r**(8.0)*rhochr(ir) + sum6
      enddo
      barmom = sum*frpirs/nprot
      rms2   = sqrt(sum2*frpirs/nprot)
      rms4   = (sum4*frpirs/nprot)**(0.25)
      rms6   = (sum6*frpirs/nprot)**(one/6.0)
c      write(6,'(a,g12.4)') ' barrett-moment=',barmom
c
c     search the barrett-radius with the secant method
c
      x1    = 0.5*rdmsc
      x2    = 2.0*rdmsc
c
      fl    = barint(x1,xkact,alphac)-barmom
      fh    = barint(x2,xkact,alphac)-barmom
      if(fl*fh.gt.zero) then
       write(11,'(a)')
     &   ' error in "barrad". first guess too narrow.'
       return
      endif
      if(fl.lt.zero) then
        xl     = x1
        xh     = x2
      else
        xl     = x2
        xh     = x1
        swap   = fl
        fl     = fh
        fh     = swap
      endif
c
      dx     = xh-xl
      do j=1,maxit
        rtflsp = xl+dx*fl/(fl-fh)
        f      = barint(rtflsp,xkact,alphac)-barmom
        if(f.lt.zero) then
          del    = xl-rtflsp
          xl     = rtflsp
          fl     = f
        else
          del    = xh-rtflsp
          xh     = rtflsp
          fh     = f
        endif
        dx     = xh-xl
        if(abs(del).lt.accura.or.f.eq.zero) goto 39
      enddo
      write(11,'(a,g12.4)')
     &   ' warning: barret-radius not converged. delrad=',del
c
   39 continue
      write(6,'(a,g14.6,a,i4,3g12.4)')
     &  ' barrett-radius=',rtflsp,
     &  ' for nprot,kappa,alpha,bar.moment=',nprot,xkact,alphac,barmom
      write(11,'(a,a/a,2i4,7f13.10)')
     & ' nprot nneut  r_bar  <r^2>^1/2',
     & ' <r^4>^1/4 <r^6>^1/6   r_diff    surface',
     & 'Bmom.:',nprot,nneut,rtflsp,rms2,rms4,rms6,rdmsc,surfc
!,rdmsc,surfc
      rms4c=rms4
c
      return
      end
c-----apv_func------------ part of the hartree-fock package ----------------------
c
      REAL(8) FUNCTION apv_funct(q)
c
      include 'params.inc'
c
c     Computes weak-charge asymmetry A_PV(q) from data prepared
c     on the file 'dcs.out'.
c     Requires Z=82 and N=126 and a previious call to 'radii'.
c     'q' is the momentum at which A_PV is computed.
c
      include 'grid.inc'
      include 'wavfun.inc'

      OPEN(131,file='apv.inp')
      WRITE(131,'(2i4)') nneut,nprot
      IF(nprot==82) THEN
        WRITE(131,'(f8.5)') 0.953D0       ! electron energy [GeV]
      ELSE IF(nprot==20 .AND. nneut==28) THEN
        WRITE(131,'(f8.5)') 2.1825D0       ! electron energy [GeV]
      END IF
      WRITE(131,'(2i6)') ngrid,ngrid
      CLOSE(131)
      CALL system(
     &  '~/hafo/tex/weak-charge/DWBA/main < apv.inp > /dev/null ')      

 
      OPEN(131,file='dcs.out')
      READ(131,'(1x/1x/1x)')
      qold=-1D9
      apv_act=0D0
      DO i=1,9999
        READ(131,*,ERR=19,END=19) dum,qact,dum,dum,apv_act
        IF(qact>q) EXIT
        qold=qact
        apv_old=apv_act
      END DO

      CLOSE(131)
      IF(qold<0D0) STOP 'error in APV_FUNCT: input not read'
!     linear interpolation
      apv = apv_old + (q-qold)*(apv_act-apv_old)/(qact-qold)
      WRITE(*,*) 'APV_FUNCT:',q,qold,qact,apv_old,apv_act,apv
      apv_funct = apv

      RETURN
 19   STOP "error in APV_FUNCT: reading DCS.OUT failed"
      END
c-----invert------------ part of the hartree-fock package ----------------------
c
      subroutine invert(r,rhs,irep)
c
      include 'params.inc'
c
      dimension r(kg1),rhs(kg1)
c
      include 'potinv.inc'
      include 'grid.inc'
      include 'constn.inc'
c
       data fsx,ezw/.83333333333333,.8333333333333e-1/
c
c
c      determines 'r' as solution of the inhomogenous linear diff. eq.:
c       " (-d/dx (sqb(x)**2) d/dx + (potl(x) - umin ) ) r = rhs  "
c      for given 'rhs' where 'sqb**2' is x-dependent inverse mass,
c      'potl' a potential and 'umin' a shift constant.
c
c      the arrays in the calling list mean:
c       r      = the solution of the equation
c       rhs    = the right hand side (inhomogeneity) of the equation.
c                r and rhs may be the same arrays in the calling routine
c                rhs is then overwritten.
c      the parameters in the calling list mean:
c       irep   = switch: irep=0 means new potl,sqb such that full
c                        inversion has to be performed from scratch.
c                        irep=1 means that previous potl,sqb still apply
c                        than parts of the first pass can be taken from
c                        the previous calls, whose results have been
c                        saved on the fields ainv and binv.
c
c      information on the hamiltonian is transferred via common /potinv/:
c       sqb    = inverse square-root of x-dependent effective mass
c       potl   = effective potential before mass transformation,
c                to be prepared in subroutine 'welfu'.
c       umin   = constant energy shift
c
c      the other arrays in /potinv/ are used as workspace for this subroutine
c      they need to be kept undestroyed if 'invert' is called subsequent
c      with 'irep'>0 .
c
c      information on the grid is transferred via common /grid/.
c
c      note: the program is written for solving the schroedinger equatio
c            in radially symmetric problems. thus it implies the rhs and
c            poteff to be such that r(x) vanishes at the origin.
c
c      the solution proceeds as follows:
c
c      by introducing scaled wave functions
c       " reff = r/sqb "
c      and
c       " rhseff = rhs/sqb ".
c      the equations is transformed to a laplacian equation
c       " (-d/dx**2 + poteff(x) ) reff = rhseff "
c      with constant mass and an effective potential
c       " poteff = (potl(x)-umin)/(sqb(x)*sqb(x)) + d2b(x) "
c      where
c       " d2b = (d/dx**2 sqb)/sqb " .
c      after solving the laplacian equation, the solution is transformed
c      back by
c       " r   = reff*sqb" .
c
c      the d/dx**2 is treated effectively in 5-point precision.
c      in addition to a standard 3-point formula for d/dx**2
c      a complementing 3-point averaging of ueff*reff and rhseff
c      is used which both together make up the 5-point precision.
c      (see R. Zurmuehl, "Praktische Mathematik", ch. VII, sect. 3,
c       Springer Verlag, Berlin 1963)
c
c
c-----------------------------------------------------------------------
c
c      square-root of inverse eff. mass on 'sqb',
c      second derivative of 'sqb' on 'd2b'.
c      finally assembling  the effective potential on 'poteff'.
c       this part can be overridden if the routine is called for
c       the same 'potl' and 'sqb' again and if the workspaces
c       'sqm' and 'poteff' are not destroyed in the meantime;
c       this situation should be indicated by 'irep>0' in the call.
c
      if(irep.le.0) then
       sqm(1) = one/sqb(1)
        do i=2,ngrid
          sqm(i)    = one/sqb(i)
          poteff(i) = d2b(i)+(potl(i)-umin)*(sqm(i)*sqm(i))
        enddo
      endif
c
c      transform to effective r.h.s., store result on r.
c
       do i=1,ngrid
         r(i)   = rhs(i)*sqm(i)
       enddo
c
c     prepare parameters for kinetic energy and loop-ends.
c
        offdig = one/(rstep*rstep)
        offdg2 = offdig+offdig
        nmaxm  = ngrid-1
c
c      check repetition switch. if new potl,sqm are indicated be irep=0
c      then begin first pass eliminating the lower off-diagonal element.
c      new diagonal and upper off-diagonal elemnts are saved on a and b.
c      the second pass, forward elimination of the r.h.s., is done withi
c      the same loop; result is stored on field t.
c
       if(irep.le.0) then
c
        binv(1) = zero
        ainv(1) = one
        x      = zero
        y      = one
        rim1   = r(1)
c
        do i=2,nmaxm
c                                         forward elim. for r.h.s.
          x      = ezw*(r(i+1)+rim1)+fsx*r(i)- x*(binv(i-1)/y)
          rim1   = r(i)
          r(i)   = x
c                                         forward elim. for matrix
          ui     = poteff(i)
          z      = -offdig+ezw*ui
          binv(i) = z
          y      = offdg2+fsx*ui-z*(binv(i-1)/y)
          ainv(i) = y
        enddo
c                                         last elim. step explicitely
c                                         because r(i+1) is not availabl
        r(ngrid) = ezw*rim1+fsx*r(ngrid)-x*(binv(nmaxm)/y)
        ui     = poteff(ngrid)
        z      = -offdig+ezw*ui
        binv(ngrid) = z
        y      = offdg2+fsx*ui-z*(binv(nmaxm)/y)
        ainv(ngrid) = y
c
c      if unchanged potl,sqm are indicated by irep=1 then perform
c      in the second pass the forward elimination for the r.h.s.
c      using the diagonal and off-diagonal elements of the previous
c      call.
c
       else
        x      = zero
        rim1   = r(1)
        do i=2,nmaxm
          x    = ezw*(r(i+1)+rim1)+fsx*r(i)- x*binv(i-1)/ainv(i-1)
          rim1 = r(i)
          r(i) = x
        enddo
        r(ngrid)  = (ezw*rim1+fsx*r(ngrid)-x*binv(nmaxm)/ainv(nmaxm))
       endif
c
c      in the third pass resolve the remaining upper triangular equation
c      backward.
c
       x      = r(ngrid)/ainv(ngrid)
       r(ngrid) = x*sqm(ngrid)
       do i=nmaxm,2,(-1)
         x      = (r(i)-x*binv(i+1))/ainv(i)
         r(i)   = x*sqm(i)
       enddo
       r(1)   = zero
c
       return
       end
c-----hpsi  ------------ part of the hartree-fock package ----------------------
c
      subroutine hpsi(wfin,hwfout,d1wfin,d2wfin,ifd1wf,jsa,iso,subfac)
c
      include 'params.inc'
c
c     action of the mean-field hamiltonian
c
c                     ^
c         |hwfout> = (h - subfac) |wfin>
c
c     where:
c       wfin    = input wave function field ,
c       hwfout  = output wave function field ,
c       d1wfin  = first derivative of 'wfin'.
c       d2wfin  = second derivative of 'wfin'.
c       ifd1wf  = decides on the computation of 'd1wfin':
c                 ifd1wf>0 --> the 'd1wfin' is interpreted as
c                              workspace and the first derivative
c                              is computed anew;
c                 ifd1wf=0 --> the 'd1wfin' contains already the
c                              first derivative if 'wfin'.
c       jsa     = combined angular momentum quantum number
c                 of wfin and hwfout,
c       iso     = isospin of the wave functions,
c       subfac  = constant energy subtracted from hamiltonian.
c
c     the fields for the mean-field hamiltonian are transferred
c     via common /wavfun/:
c       pot     = potential
c       potls   = spin-orbit potential
c       bmass   = inverse effective mass
c       d1b     = derivative of inverse effective mass
c
c     grid properties are transferred via common /grid/
c
      dimension wfin(kg1),hwfout(kg1),d1wfin(kg1),d2wfin(kg1)
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'constn.inc'
c
c-------------------------------------------------------------------------------
c
c     the quantum numbers of the wave function
c      la      = orbital angular momentum of state  wfin>
c      parity  = parity for the radial wave function (note that it
c                contains an extra factor 'r' such that the 1s1/2
c                state has parity=-1 etc)
c      isoadd  = stack-pointer for potential fields pot, potls
c                and isoadd
c      xl0     = spin-orbit splitting for state  wfin>
c
      la     = jsa/2
      xlp    = (la*la+la)
      parity = (mod(la,2)*2-1)
      isoadd = (iso-1)*ngrid
      xl0    = wls(jsa)
c
c     the derivatives of the wave function
c
      if(ifd1wf.gt.0)
     &   call derv17(wfin(1),d1wfin(1),parity,ngrid,rstep)
      call derv25(wfin(1),d2wfin(1),parity,ngrid,rstep)
c
c     the action of the hamiltonian
c
      hwfout(1) = zero
      do ir=2,ngrid
        irs    = ir+isoadd
        potact = pot(irs)+xl0*potls(irs)
     &          +(bmass(irs)*zentrf(ir)+bmassl(irs))*xlp
     &          +d1b(irs)*rinvs(ir)
        hwfout(ir) = (potact-subfac)*wfin(ir)
     &              -d1b(irs)*d1wfin(ir)-bmass(irs)*d2wfin(ir)
      enddo
c
      return
      end
c-----tpsi  ------------ part of the hartree-fock package ---------------------
c
      subroutine tpsi(wfin,hwfout,d2wfin,jsa,iso,subfac)
c
      include 'params.inc'
c
c     action of the kinetic energy
c
c     where:
c       wfin    = input wave function field ,
c       hwfout  = output wave function field ,
c       jsa     = combined angular momentum q.n. of wfin and hwfout,
c       iso     = isospin of the wave functions,
c       subfac  = constant energy subtracted from kinetic energy operator
c
c     grid properties are transferred via common /grid/
c
      dimension wfin(kg1),hwfout(kg1),d2wfin(kg1)
c
      include 'grid.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c-------------------------------------------------------------------------------
c
c     the quantum numbers of the wave function
c      la      = orbital angular momentum of state  wfin>
c      parity  = parity for the radial wave function (note that it
c                contains an extra factor 'r' such that the 1s1/2
c                state has parity=-1 etc)
c      isoadd  = stack-pointer for potential fields pot, potls
c                and isoadd
c
      la     = jsa/2
      xlp    = (la*la+la)
      parity = (mod(la,2)*2-1)
cold      isoadd = (iso-1)*ngrid
c
c     the derivatives of the wave function
c
      call derv25(wfin,d2wfin,parity,ngrid,rstep)
c
c     the action of the kinetic energy
c
      amfcac = amfc(iso)
      hwfout(1) = amfcac*d2wfin(1)
      do 10 ir=2,ngrid
c      r      = (ir-1)*rstep
      potact = amfcac*xlp*zentrf(ir)
      hwfout(ir) = (potact-subfac)*wfin(ir) - amfcac*d2wfin(ir)
   10 continue
c
      return
      end
c-----k2psi  ------------ part of the hartree-fock package ---------------------
c
      subroutine k2psi(wfin,hwfout,d2wfin,jsa)
c
      include 'params.inc'
c
c     action of the operator -nabla**2 
c
c     where:
c       wfin    = input wave function field ,
c       hwfout  = output wave function field ,
c       jsa     = combined angular momentum q.n. of wfin and hwfout,
c
c     grid properties are transferred via common /grid/
c
      dimension wfin(kg1),hwfout(kg1),d2wfin(kg1)
c
      include 'grid.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c-------------------------------------------------------------------------------
c
c     the quantum numbers of the wave function
c      la      = orbital angular momentum of state  wfin>
c      parity  = parity for the radial wave function (note that it
c                contains an extra factor 'r' such that the 1s1/2
c                state has parity=-1 etc)

c
      la     = jsa/2
      xlp    = (la*la+la)
      parity = (mod(la,2)*2-1)
cold      isoadd = (iso-1)*ngrid
c
c     the derivatives of the wave function
c
      call derv25(wfin,d2wfin,parity,ngrid,rstep)
c
c     action of -nabla**2
c
      hwfout(1) = d2wfin(1)
      DO ir=2,ngrid
        hwfout(ir) = xlp*zentrf(ir)*wfin(ir) - d2wfin(ir)
      END DO
c
      return
      end
c-----derv13------------ part of the hartree-fock package ----------------------
c
      subroutine derv13(wfin,dwfout,parity,ngrid,delr)
c
      include 'params.inc'
c
c     first derivative of 'wfin' on 'dwfout',
c     using three point formula (which effectively is a two point formula).
c
c     input via calling list:
c      wfin    = array for which derivative is to be computed
c      parity  = the radial parity of the field 'wfin'
c                ( > 0 if the field starts with even power of r, < 0 else)
c      ngrid   = length of the arrays
c      delr    = distance between grid points
c     output via calling list:
c      dwfout  = array containing the derivative
c
      dimension wfin(ngrid),dwfout(ngrid)
c
c---------------------------------------------------------------------------
c
c     compute frequently occuring parameters
c
      delind = 2.0/delr
c
c     derivative at the origin, use 'parity' to continue the wave function
c
      if(parity.gt.0.) then
       dwfout(1) = 0.0
      else
       dwfout(1) = (delind+delind)*wfin(2)
      endif
c
c     the straightforward "mid-region"
c
      do i=2,ngrid-1
        dwfout(i) = delind*(wfin(i+1)-wfin(i-1))
      enddo
c
c     right-handed two point formula at the boundary
c
      i=ngrid
      dwfout(i) = delind*(wfin(i)-wfin(i-1))
c
      return
      end
c-----derv23------------ part of the hartree-fock package ----------------------
c
      subroutine derv23(wfin,dwfout,parity,ngrid,delr)
c
      include 'params.inc'
c
c     second derivative of 'wfin' on 'dwfout',
c     using three point formula.
c
c     input via calling list:
c      wfin    = array for which 2.derivative is to be computed
c      parity  = the radial parity of the field 'wfin'
c                ( > 0 if the field starts with even power of r, < 0 else)
c      ngrid   = length of the arrays
c      delr    = distance between grid points
c     output via calling list:
c      dwfout  = array containing the 2. derivative
c
      dimension wfin(ngrid),dwfout(ngrid)
c
      data one,zero/1.0,0.0/
c
c----------------------------------------------------------------------
c
c     compute frequently occuring parameters
c
      delin2 = one/(delr*delr)
      delid2 = delin2+delin2
c
c     derivative at the origin, use 'parity' to continue the wave function
c
      if(parity.gt.0.) then
       dwfout(1) = delid2*(wfin(2)-wfin(1))
      else
       dwfout(1) = zero
      endif
c
c     the straightforward "mid-region"
c
      do i=2,ngrid-1
        dwfout(i) = -delid2*wfin(i)+delin2*(wfin(i-1)+wfin(i+1))
      enddo
c
c     extrapolate linearly to grid boundary
c
      dwfout(ngrid) = 2.0*dwfout(ngrid-1)-dwfout(ngrid-2)
c
      return
      end
c-----derv15------------ part of the hartree-fock package ----------------------
c
      subroutine derv15(wfin,dwfout,parity,ngrid,delr)
c
      include 'params.inc'
c
c     first derivative of 'wfin' on 'dwfout',
c     using five point formula.
c
c     this version uses loop unrolling in a way which inhibits
c     vectorization. a vectorizable variant is stored on the
c     lines which are initiated by 'cvec'
c
c     input via calling list:
c      wfin    = array for which 2.derivative is to be computed
c      parity  = the radial parity of the field 'wfin'
c                ( > 0 if the field starts with even power of r, < 0 else)
c      ngrid   = length of the arrays
c      delr    = distance between grid points
c     output via calling list:
c      dwfout  = array containing the 2. derivative
c
      dimension wfin(ngrid),dwfout(ngrid)
c
      data              vdr,               esx,               zdr
     &    /1.33333333333333,0.1666666666666667,0.6666666666666667/
      data                 ezw
     &    /0.08333333333333333/
cold      data                 ezw,               fsx,            fzwt
cold     &    /0.08333333333333333,0.8333333333333333,2.08333333333333/
      data one,zero/1.0,0.0/
c
c------------------------------------------------------------------------------
c
c     compute frequently occuring parameters
c
      delinv = one/delr
      zdrh   = zdr*delinv
      ezwh   = ezw*delinv
c
c     fill the first set of local wavefunctions
c
      wf1    = wfin(1)
      wf2    = wfin(2)
      wf3    = wfin(3)
      wf4    = wfin(4)
c
c     derivatives near the origin, use 'parity' to continue the wave function
c
      if(parity.gt.0.) then
       dwfout(1) = zero
      else
       dwfout(1) = delinv*(vdr*wf2-esx*wf3)
      endif
      dwfout(2) = zdrh*(wf3-wf1)
     *            -ezwh*(wf4-parity*wf2)
cvec      if(parity.gt.0.) then
cvec       dwfout(1) = zero
cvec      else
cvec       dwfout(1) = delinv*(vdr*wfin(2)-esx*wfin(3))
cvec      endif
cvec      dwfout(2) = zdrh*(wfin(3)-wfin(1))
cvec     *            -ezwh*(wfin(4)-parity*wfin(2))
c
c     the straightforward "mid-region"
c
      ngm2 = ngrid-2
      do i=3,ngm2,5
        wf5 = wfin(i+2)
        dwfout(i) = ezwh*(wf1-wf5) +zdrh*(wf4-wf2)
        ia = i+1
        if(ia.gt.ngm2) goto 90
        wf1 = wfin(ia+2)
        dwfout(ia) = ezwh*(wf2-wf1) +zdrh*(wf5-wf3)
        ia = ia+1
        if(ia.gt.ngm2) goto 90
        wf2 = wfin(ia+2)
        dwfout(ia) = ezwh*(wf3-wf2) +zdrh*(wf1-wf4)
        ia = ia+1
        if(ia.gt.ngm2) goto 90
        wf3 = wfin(ia+2)
        dwfout(ia) = ezwh*(wf4-wf3) +zdrh*(wf2-wf5)
        ia = ia+1
        if(ia.gt.ngm2) goto 90
        wf4 = wfin(ia+2)
        dwfout(ia) = ezwh*(wf5-wf4) +zdrh*(wf3-wf1)
      enddo
cvec      do 10 i=3,ngrid-2
cvec      dwfout(i) = ezwh*(wfin(i-2)-wfin(i+2))
cvec     *           +zdrh*(wfin(i+1)-wfin(i-1))
cvec   10 continue
c
c     near the boundaries lower precision and extrapolate
c
   90 continue
      i=ngrid-1
      dwfout(i) = 0.5*delinv*(wfin(i+1)-wfin(i-1))
      dwfout(ngrid) = 2.0*dwfout(ngrid-1)-dwfout(ngrid-2)
c
      return
      end
c-----derv25------------ part of the hartree-fock package ----------------------
c
      subroutine derv25(wfin,dwfout,parity,ngrid,delr)
c
      include 'params.inc'
c
c
c     second derivative of 'wfin' on 'dwfout',
c     using five point formula.
c
c     this version uses loop unrolling in a way which inhibits
c     vectorization. a vectorizable variant is stored on the
c     lines which are initiated by 'cvec'
c
c     input via calling list:
c      wfin    = array for which 2.derivative is to be computed
c      parity  = the radial parity of the field 'wfin'
c                ( > 0 if the field starts with even power of r, < 0 else)
c      ngrid   = length of the arrays
c      delr    = distance between grid points
c     output via calling list:
c      dwfout  = array containing the 2. derivative
c
      dimension wfin(ngrid),dwfout(ngrid)
c
      data              vdr,                  ezw
     &    /1.33333333333333,0.0833333333333333333/
      data one,zero/1.0,0.0/
c
c----------------------------------------------------------------------
c
c     compute frequently occuring parameters
c
      delin2 = one/(delr*delr)
      fhhh   = 2.5*delin2
      vdrhh  = vdr*delin2
      ezwhh  = ezw*delin2
c
c     fill the first set of local wavefunctions
c
      wf1    = wfin(1)
      wf2    = wfin(2)
      wf3    = wfin(3)
      wf4    = wfin(4)
c
c     derivatives near the origin, use 'parity' to continue the wave function
c
      if(parity.gt.0.) then
       dwfout(1) = 2.0*(vdrhh*wf2-ezwhh*wf3)-fhhh*wf1
      else
       dwfout(1) = zero
      endif
      dwfout(2) = -fhhh*wf2+vdrhh*(wf3+wf1)
     *            -ezwhh*(wf4+parity*wf2)
cvec      if(parity.gt.0.) then
cvec       dwfout(1) = 2.0*(vdrhh*wfin(2)-ezwhh*wfin(3))-fhhh*wfin(1)
cvec      else
cvec       dwfout(1) = zero
cvec      endif
cvec      dwfout(2) = -fhhh*wfin(2)+vdrhh*(wfin(3)+wfin(1))
cvec     *            -ezwhh*(wfin(4)+parity*wfin(2))
c
c     the straightforward "mid-region"
c
      ngm2 = ngrid-2
      do 10 i=3,ngm2,5
      wf5 = wfin(i+2)
      dwfout(i) = -fhhh*wf3+vdrhh*(wf2+wf4)-ezwhh*(wf1+wf5)
      ia = i+1
      if(ia.gt.ngm2) goto 90
      wf1 = wfin(ia+2)
      dwfout(ia) = -fhhh*wf4+vdrhh*(wf3+wf5)-ezwhh*(wf2+wf1)
      ia = ia+1
      if(ia.gt.ngm2) goto 90
      wf2 = wfin(ia+2)
      dwfout(ia) = -fhhh*wf5+vdrhh*(wf4+wf1)-ezwhh*(wf3+wf2)
      ia = ia+1
      if(ia.gt.ngm2) goto 90
      wf3 = wfin(ia+2)
      dwfout(ia) = -fhhh*wf1+vdrhh*(wf5+wf2)-ezwhh*(wf4+wf3)
      ia = ia+1
      if(ia.gt.ngm2) goto 90
      wf4 = wfin(ia+2)
      dwfout(ia) = -fhhh*wf2+vdrhh*(wf1+wf3)-ezwhh*(wf5+wf4)
   10 continue
cvec      do 10 i=3,ngrid-2
cvec      dwfout(i) = -fhhh*wfin(i)+vdrhh*(wfin(i-1)+wfin(i+1))
cvec     *                      -ezwhh*(wfin(i-2)+wfin(i+2))
cvec   10 continue
c
c     use only 3-point formula one grid-point before boundary
c
   90 continue
cold      dwfout(ngrid-1)= (wfin(ngrid)-2.0*wfin(ngrid-1)+wfin(ngrid-2))
cold     *                *delin2
coldc
coldc     extrapolate linearly to grid boundary
coldc
cold      dwfout(ngrid) = 2.0*dwfout(ngrid-1)-dwfout(ngrid-2)
c
c     boundary condition of vanishing wave function
c
      dwfout(ngrid-1) = -fhhh*wfin(ngrid-1)
     &                  +vdrhh*(wfin(ngm2)+wfin(ngrid))
     *                  -ezwhh*(wfin(ngrid-3)    )
      dwfout(ngrid) = -fhhh*wfin(ngrid)
     &                  +vdrhh*(wfin(ngrid-1)    )
     *                  -ezwhh*(wfin(ngm2)    )
c
      return
      end
c-----derv17------------ part of the hartree-fock package ----------------------
c
      subroutine derv17(wfin,dwfout,parity,ngrid,delr)
c
      include 'params.inc'
c
c
c     first derivative of 'wfin' on 'dwfout',
c     using seven point formula.
c
c     this version uses loop unrolling in a way which inhibits
c     vectorization. a vectorizable variant is stored on the
c     lines which are initiated by 'cvec'
c
c     input via calling list:
c      wfin    = array for which 2.derivative is to be computed
c      parity  = the radial parity of the field 'wfin'
c                ( > 0 if the field starts with even power of r, < 0 else)
c      ngrid   = length of the arrays
c      delr    = distance between grid points
c     output via calling list:
c      dwfout  = array containing the 2. derivative
c
      dimension wfin(ngrid),dwfout(ngrid)
c
      data                 ezw,                zdr,                esxt
     &    /0.08333333333333333,0.66666666666666667,0.016666666666666667/
      data one,zero/1.0,0.0/
c
c----------------------------------------------------------------------
c
c     compute frequently occuring parameters
c
      delinv = one/delr
      zdrh   = zdr*delinv
      ezwh   = ezw*delinv
      esxth  = esxt*delinv
      xnsxh  = 9.0*esxth
      f5sxh  = 45.0*esxth
c
c     fill the first set of local wavefunctions
c
      wf1    = wfin(1)
      wf2    = wfin(2)
      wf3    = wfin(3)
      wf4    = wfin(4)
      wf5    = wfin(5)
      wf6    = wfin(6)
c
c     derivatives near the origin, use 'parity' to continue the wave function
c
      if(parity.gt.0.) then
       dwfout(1)= zero
      else
       dwfout(1)= 2.0*(esxth*wf4-xnsxh*wf3
     *                 +f5sxh*wf2)
      endif
      dwfout(2)= esxth*(wf5-parity*wf3)
     *          -xnsxh*(wf4-parity*wf2)
     *          +f5sxh*(wf3-wf1)
      dwfout(3)= esxth*(wf6-parity*wf2)
     *          -xnsxh*(wf5-wf1)
     *          +f5sxh*(wf4-wf2)
cvec      if(parity.gt.0.) then
cvec       dwfout(1)=0.
cvec      else
cvec       dwfout(1)= 2.0*(esxth*wfin(4)-xnsxh*wfin(3)
cvec     *                 +f5sxh*wfin(2))
cvec      endif
cvec      dwfout(2)= esxth*(wfin(5)-parity*wfin(3))
cvec     *          -xnsxh*(wfin(4)-parity*wfin(2))
cvec     *          +f5sxh*(wfin(3)-wfin(1))
cvec      dwfout(3)= esxth*(wfin(6)-parity*wfin(2))
cvec     *          -xnsxh*(wfin(5)-wfin(1))
cvec     *          +f5sxh*(wfin(4)-wfin(2))
c
c     the straightforward "mid-region"
c
      ngm3 = ngrid-3
      do 10 i=4,ngm3,7
      wf7 = wfin(i+3)
      dwfout(i) = esxth*(wf7-wf1) -xnsxh*(wf6-wf2) +f5sxh*(wf5-wf3)
      ia = i+1
      if(ia.gt.ngm3) goto 90
      wf1 = wfin(ia+3)
      dwfout(ia)= esxth*(wf1-wf2) -xnsxh*(wf7-wf3) +f5sxh*(wf6-wf4)
      ia = ia+1
      if(ia.gt.ngm3) goto 90
      wf2 = wfin(ia+3)
      dwfout(ia)= esxth*(wf2-wf3) -xnsxh*(wf1-wf4) +f5sxh*(wf7-wf5)
      ia = ia+1
      if(ia.gt.ngm3) goto 90
      wf3 = wfin(ia+3)
      dwfout(ia)= esxth*(wf3-wf4) -xnsxh*(wf2-wf5) +f5sxh*(wf1-wf6)
      ia = ia+1
      if(ia.gt.ngm3) goto 90
      wf4 = wfin(ia+3)
      dwfout(ia)= esxth*(wf4-wf5) -xnsxh*(wf3-wf6) +f5sxh*(wf2-wf7)
      ia = ia+1
      if(ia.gt.ngm3) goto 90
      wf5 = wfin(ia+3)
      dwfout(ia)= esxth*(wf5-wf6) -xnsxh*(wf4-wf7) +f5sxh*(wf3-wf1)
      ia = ia+1
      if(ia.gt.ngm3) goto 90
      wf6 = wfin(ia+3)
      dwfout(ia)= esxth*(wf6-wf7) -xnsxh*(wf5-wf1) +f5sxh*(wf4-wf2)
   10 continue
cvec      do 10 i=4,ngrid-3
cvec      dwfout(i)= esxth*(wfin(i+3)-wfin(i-3))
cvec     *          -xnsxh*(wfin(i+2)-wfin(i-2))
cvec     *          +f5sxh*(wfin(i+1)-wfin(i-1))
cvec   10 continue
c
c     near the boundaries lower precision and extrapolate
c
   90 continue
      i=ngrid-2
      dwfout(i) = ezwh*(wfin(i-2)-wfin(i+2))
     *           +zdrh*(wfin(i+1)-wfin(i-1))
      i=ngrid-1
      dwfout(i) = 0.5*delinv*(wfin(i+1)-wfin(i-1))
      dwfout(ngrid) = 2.0*dwfout(ngrid-1)-dwfout(ngrid-2)
c
      return
      end
c-----lapfil------------------------------------------------------------
c
       subroutine lapfil(rfilt,awork,xmas2,irep)
c
      include 'params.inc'
c
       dimension rfilt(kg1),awork(kg1)
c
       include 'grid.inc'
c
       data fsx,ezw/.83333333,.83333333e-1/
ctest       data fsx,ezw/1.0,0.0/
c
c      the "laplace-folding-filter" for the skyrme-hartree-fock
c      package.
c      solves the inhomogenous linear diff. eq.:
c       " xmas2^-1 * (-d/dr**2 - 2.0*d/dr + xmas2) rfield = rfield  "
c
c      the d/dx**2 is treated effectively in 5-point precision.
c      in addition to a standard 3-point formula for d/dx**2
c      a complementing 3-point averaging of xmas2*reff and rhseff
c      is used which both together make up the 5-point precision.
c
c
c      the arrays in the call mean:
c       rfilt  = field to be filtered at input and result at output.
c       awork  = work space to keep the off-diagonal element for inversion.
c      the parameters in the call mean:
c       xmas2  = parameter in the laplacian equation, means mass**2 in
c                case of meson folding
c       rstep  = stepsize for the grid
c       ngrid   = maximum number of grid points
c       irep   = switch: irep=0 means new potl,sqb such that full
c                        inversion has to be performed from scratch.
c                        irep=1 means that previous potl,sqb still apply
c                        than parts of the first pass can be taken from
c                        the previous calls, whose results have been
c                        saved on the fields a and b.
c
c-----------------------------------------------------------------------
c
c      transform to cartesian volume element. 'rfilt' is overwritten!
c      and save the normalization.
c
       xnorm = 0.0
       do ir=1,ngrid
         xnorm = rsquar(ir)*rfilt(ir) + xnorm
         rfilt(ir) = rfact(ir)*rfilt(ir)
       enddo
       xnorm = frpirs*xnorm
c
c      check repetition switch. if new xmas2 is indicated be irep=0
c      then begin first pass eliminating the lower off-diagonal element.
c      new diagonal and upper off-diagonal elemnts are saved on a and b.
c      the second pass, forward elimination of the r.h.s., is done withi
c      the same loop; result is stored on field t.
c
c    prepare some parameters for kinetic energy and loop-ends.
c
c       the original on- and off-diagonal elements
c
        rinv2  = 1./(rstep*rstep)
        offdig = -rinv2+ezw*xmas2
        ondiag = rinv2+rinv2+fsx*xmas2
c
c       rescaling by 1/offdig acts on diagonal element and on r.h.s.
c
        offdgi = 1.0/offdig
        ondgrs = ondiag*offdgi
        ezwrs  = ezw*offdgi
        fsxrs  = fsx*offdgi
        nmaxm  = ngrid-1
c
       if(irep.le.0) then
c                                 if(irep.le.0) then
c
        awork(1) = 1.0
        x      = 0.
        y      = 1.0
c
c                                         i=2 explicitely
        x      = ezwrs*(rfilt(3)+rfilt(1))+fsxrs*rfilt(2)
        rfsav  = rfilt(2)
        rfilt(2)   = x
        y      = 1.0/ondgrs
        awork(2)   = y
c                                         now the most other i
        do 13 i=3,nmaxm
c                                         forward elim. for r.h.s.
        x      = ezwrs*(rfilt(i+1)+rfsav)+fsxrs*rfilt(i)- x*y
        rfsav  = rfilt(i)
        rfilt(i)   = x
c                                         forward elim. for matrix
        y      = 1.0/(ondgrs-y)
        awork(i) = y
   13   continue
c                                         last elim. step explicitely
c                                         because rfilt(i+1) is not availabl
        rfilt(ngrid) = ezwrs*rfsav+fsxrs*rfilt(ngrid)-x*y
cold        y      = 1.0/(ondiag-y)
        y      = 1.0/(ondgrs-y)
        awork(ngrid) = y
c
c      if unchanged xmas2 is indicated by irep=1 then perform
c      in the second pass the forward elimination for the r.h.s.
c      using the diagonal and off-diagonal elements of the previous
c      call.
c
       else
        x      = 0.
        rfsav  = rfilt(1)
        do 14 i=2,nmaxm
        x      = ezwrs*(rfilt(i+1)+rfsav)+fsxrs*rfilt(i)- x*awork(i-1)
        rfsav  = rfilt(i)
        rfilt(i)   = x
   14   continue
        rfilt(ngrid)  = (ezwrs*rfsav+fsxrs*rfilt(ngrid)-x*awork(nmaxm))
       endif
c
c      in the third pass resolve the remaining upper triangular equation
c      backward.
c
       x      = rfilt(ngrid)*awork(ngrid)
       rfilt(ngrid) = x
       do 15 i=nmaxm,2,(-1)
       x      = (rfilt(i)-x)*awork(i)
       rfilt(i)   = x
   15  continue
c
c      final division by r and normalization
c
       sum = 0.0
       do ir=2,ngrid
         rfilt(ir) = xmas2*rinvs(ir)*rfilt(ir)
         sum = rsquar(ir)*rfilt(ir) + sum
       enddo
       xnorm = xnorm/(frpirs*sum)
       do ir=2,ngrid
         rfilt(ir) = xnorm*rfilt(ir)
       enddo
       rfilt(1)=(4.0*rfilt(2)-rfilt(3))*0.33333333
c7       rfilt(1) = 1.5*rfilt(2)-0.6*rfilt(3)+0.1*rfilt(4)
c
       return
       end
c-----t0inv ------------------------------------------------------
c
      subroutine t0inv(wfin,wfout,workt,dampfc,ll)
c
      include 'params.inc'
c
c     action of the inverse kinetic energy for the single wavefunction
c     'wfin'. result is on 'wfout'. both wavefunctions are without
c     r-factor. orbital angular momentum is 'll'.
c     'workt' is a workspace.
c
c     some auxiliary fields are needed: 'npa' for intermediate storage
c     of the wave function pointers 'npl'; 'wls' = spin-orbit splitting
c     for given 'js'; 'pkin' = the zentrifugal force for given orbital
c     ang.mom.; 'd1b' = d/dr 'bmass' / r  .
c
      dimension wfin(kg1),wfout(kg1),workt(kg1)
c
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'grid.inc'
cold      include 'linres.inc'
      include 'numpar.inc'
      include 'potinv.inc'
      include 'constn.inc'
c
c-------------------------------------------------------------------------------
c
c     the effective potential, mass and 'umin' for 'invert'.
c
      sqd    = sqrt(dampfc)
      xlpfc  = (ll*ll+ll)*dampfc
      r      = -rstep
      do 60 ir=1,ngrid
      r      = rstep+r
      potl(ir) = one+xlpfc*zentrf(ir)
      workt(ir) = r*wfin(ir)
      sqb(ir) = sqd
      d2b(ir) = zero
   60 continue
      umin   = zero
c       write(6,'(a,20(/1x,5g12.4))')
c     &  ' wfin:',(wfin(i),i=1,ngrid)
c       write(6,'(a,20(/1x,5g12.4))')
c     &  ' work1:',(work1(i),i=1,ngrid)
c
      call invert(workt,workt,-1)
c
c     test inversion
c
c      ifinvt = 1
c      if(ifinvt.eq.1) then
c       write(6,'(3(a,i3),a,g12.4)')
c     &   ' test inversion at iso=',iso,' jsa=',jsa,'  js=',js
c       call tpsi(work1(1),work3(1),aux1(1),2*ll,1,-amfc(1)/dampfc)
c       afac   = dampfc/amfc(1)
c       sumin  = zero
c       sumout = zero
c       sumio  = zero
c       r      = zero
c       do 180 ir=2,ngrid
c       r      = rstep+r
c       work3(ir) = afac*work3(ir)
c       sumin  = wfin(ir)*wfin(ir)*r*r + sumin
c       sumio  = wfin(ir)*r*work3(ir)  + sumio
c       sumout = work3(ir)*work3(ir)   + sumout
c  180  continue
c       angle  = one - sumio*sumio/max(sumin*sumout,xsmall)
c       write(6,'(a,3g12.4,a,g12.4)')
c     &   '  in,out,in-out=',sumin,sumout,sumio,
c     &   '  angle**2=',angle
c      endif
c      ifwfts = 0
c      if(ifwfts.eq.1) then
c       write(6,'(a,20(/1x,5g12.4))')
c     &  ' wfin:',(wfin(i)*(i-1)*rstep,i=1,ngrid)
c       write(6,'(a,20(/1x,5g12.4))')
c     &  ' invert:',(work1(i),i=1,ngrid)
c       write(6,'(a,20(/1x,5g12.4))')
c     &  ' hpsi*invert:',(work3(i),i=1,ngrid)
c     endif
c
c     remove r-factor
c
      do ir=2,ngrid
        wfout(ir) = rinvs(ir)*workt(ir)
      enddo
      parl = (1-mod(ll,2)*2)
      if(parl.gt.zero) then
c7       wfout(1) = 1.5*wfout(2)-0.6*wfout(3)+0.1*wfout(4)
       wfout(1)=(4.0*wfout(2)-wfout(3))*0.33333333
      else
       wfout(1) = zero
      endif
c       write(6,'(a,20(/1x,5g12.4))')
c     &  ' wfout:',(wfout(i),i=1,ngrid)
c       write(6,'(a,20(/1x,5g12.4))')
c     &  ' work1:',(work1(i),i=1,ngrid)
c
c
      return
      end
c-----filter0------------------------------------------------------------------
c
      subroutine filter0(fldin,fldout)
c
      include 'params.inc'
c
c     folding operation using Fourier-Bessel transformation for L=0
c     input on 'fldin'. output on 'fldout'. ('fldin' and 'fldout'
c     can be the same fields).
c     the workspaces 'aux1' and 'aux2' are needed.
c     The folding function does enhance weights of high q.
c
      dimension fldin(kg1),fldout(kg1)
      dimension auxfou(kgr2)
c
      include 'wavfun.inc'
      include 'linres.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c      equivalence (aux1(1),auxfou(1))
c
      logical tposdef
c
c     anti-folding function, in Fourier space enhancing high 'q'.
c     parameters via common /skyapr/:
c       'ayuk'   = amplitude of enahniong high 'q'
c       'rhosnm' = nuclear matter density of given force.
c
      antifold(q) = 1.0D0+ayuk*(TANH((q-fermp)/0.027D0)+1.0D0)
c
c-------------------------------------------------------------------
c
c
c     formfactor of 'fldin',
c     augmented by folding weight,
c     and stored on workspace 'auxfou'.
c
c      write(38,'(a,2g12.4)')
c     &  '  RSTEP,DELQ=',rstep,delq
c      write(38,'(30(/1x,5g12.4))') fldin
      fermp = (1.5D0*pi*pi*rhosnm)**(0.333333333)
      sum = 0.0D0
      r   = 0.0D0
      tposdef=.TRUE.
      do ir=2,ngrid
        r      = rstep+r
        sum    = (r*r)*fldin(ir)+sum
        if(tposdef) tposdef = fldin(ir).ge.0.0D0
      enddo
      auxfou(1) = rstep*sum 
      do iq=2,ngrid
        q1     = iq*delq-delq
        sum    = zero
        dc     = cos(q1*rstep)
        ds     = sin(q1*rstep)
        r      = 0.0D0
        c      = one
        s      = zero
        do ir=2,ngrid
          r      = rstep+r
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
          sum    = r*s*fldin(ir)+sum
        enddo
        auxfou(iq) = rstep*sum/q1   * antifold(q1)
c        write(38,'(i5,f8.3,3(1pg13.5))') 
c     &     iq,q1,abs(rstep*sum/q1),abs(auxfou(iq)),antifold(q1)
      enddo
c
c     fourier-bessel backtransformation
c
      dpifac   = 2.0D0/pi
      sum    = zero
      qa     = zero
      do iq=2,ngrid
        qa     = delq+qa
        sum    = (qa*qa)*auxfou(iq)+sum
      enddo
      fldout(1) = sum*delq*dpifac
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
        do iq=2,ngrid
          qa     = delq+qa
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
          sum    = qa*s*auxfou(iq)+sum
        enddo
        if(tposdef) then
          fldout(ir) = abs(sum*delq*dpifac/r)
        else
          fldout(ir) = (sum*delq*dpifac/r)
        endif
c
      enddo
c
      sum = 0.0D0
      r   = 0.0D0
      do ir=2,ngrid
        r      = rstep+r
        sum    = (r*r)*fldout(ir)+sum
      enddo
      write(38,'(a,2g12.4)')
     &   '# sum IN,OUT=',auxfou(1),rstep*sum
      if(tposdef) then      
        renorm = auxfou(1)/(rstep*sum)
c      write(38,'(a,g12.4)') ' renomr=',renorm
        do ir=2,ngrid
          fldout(ir) = fldout(ir)*renorm
c        write(38,'(i3,1pg12.4)') ir,fldout(ir)
        enddo
      endif
c
      return
      end

c-----buffio------------ part of the hartree-fock package ----------------------
c
      subroutine buffio(direc)
c
      include 'params.inc'
c
c     buffering the single particle wave functions on scratch file 19.
c
c     'direc' determines the direction of the action:
c      direc = 'i' --> wave functions for the actual nucleus are read
c                      from file 19. before that, the file-table 'nfile'
c                      is checked to compute the pointer.
c      direc = 'o' --> wave functions of the actual nucleaus are written
c                      on file 19. the file table 'nfile' is checked to
c                      compute the pointer. the file table 'nfile' is
c                      extended, if the nucleus is not yet contained.
c
c     the wave functions and other information about the actual nucleus
c     are transferred via common's.
c
c     the following parameter statement determines the partitioning
c     of the i/o:
c      kwrdln  = number of words in one records read written;
c                this number needs to be adjusted to the actual
c                computer for most efficient buffering.
c      kspfil  = length of a record in machine units (usually bytes);
c                the factor 'kvarln' in 'params.inc' needs to be set
c                anew for each computer.
c      kfileb  = factor of proton-number in the combined key-values
c                'ntest' or 'nfile(1,i)'.
c
!      parameter (kwrdln=9192,kfileb=200,
!      parameter (kwrdln=131072,kfileb=256,
      parameter (kwrdln=65536,kfileb=256,
     &           kspfil=kvarln*kwrdln)
c
      character*1 direc
      integer     ntest,nfile(2,200)
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'result.inc'
c
c     teh following variables are "initialised" in data to make sure that
c     they are saved after leaving the routine.
c
      data ifopen,markfl,maxfl,ntest,nfile/0,1,1,0,400*0/
      LOGICAL,PARAMETER :: ttest=.TRUE.
c
c-------------------------------------------------------------------------------
c
c     check disc-switch 'idisc'
c
      if(idisc.lt.0) return
c
c     check length of first communication field
c
      if(KNMXP*KJMX2+4*KSL2+56+15.gt.kwrdln)
     &  stop ' insufficient record length in BUFFIO'
c
c     open file 19 if necessary, 'ifopen' serves as indicator
c
      if(ifopen.eq.0) then
       open(unit=19,access='direct',form='unformatted',recl=kspfil,
     &      status='scratch')
       write(06,'(a)') ' unit 19 opened'
       ifopen = 1
      endif
c
      if(direc.eq.'i') then
c
c-------------------------------------------------------------------------------
c
c     read the wave function from file 19, or extend file-table
c     if necessary.
c
c-------------------------------------------------------------------------------
c
c      first, check file table for occurence of 'nprot' and 'nneut'
c
       ifread = 0
       if(markfl.le.0) return
       ntest  = kfileb*nprot+nneut
       if(iwr.ge.-1)
     &  write(06,'(a,12i5)')
     &  '  ntest,maxfl,nfile=',ntest,maxfl,((nfile(i,j),j=1,5),i=1,2)
       do 281 mark=1,0500
       mar2  = mark
       if(ntest.eq.nfile(1,mark) .or. nfile(1,mark).eq.0) goto 282
  281  continue
       stop 'file table exhausted'
c
  282  continue
       if(nfile(1,mar2).eq.0) then
c
c       nucleus was not on file-table 'nfile'; extend the table
c       by the new nucleus and return 'ifread=0' to indicate to
c       the calling routine that nothing has been read in.
c
        markfl = maxfl
        nfile(2,mar2) = maxfl
        nfile(1,mar2) = ntest
        ifread = 0
        IF(TTEST) WRITE(6,'(a)') 
     &     'BUFFIO: file table extended, nothing read'
       else
c
c       nucleus was found on file-table 'nfile'; extract pointer
c       to first record from 'nfile(2,...)' and store it on 'markfl'.
c       'markal' is the dynamic pointer which scans through the
c       sequence of records.
c
        markfl = nfile(2,mar2)
        markal = markfl
        read(19,rec=markal) rstep,ngrid,
     &                      nprot,nneut,nmaxp,nmaxn,
     &                      weight,v2occ,deltaf,spenrg,
     &                      efen,efep,deltp,deltn,xlamp2,xlamn2,
     &                      jsmax,ifpair,npla
c
c       compute partitions of the wave functions
c
        nmaxto = nmaxp+nmaxn
        nmaxln = nmaxto*ngrid
        nmaxbl = nmaxln/kwrdln
        nkspm  = kwrdln-1
c
c       read the wave functions
c
        do i1=1,nmaxln-kwrdln,kwrdln
          markal = 1+markal
          read(19,rec=markal) (wf(i),i=i1,i1+nkspm)
        enddo
        i1    =nmaxbl*kwrdln+1
        markal=1+markal
        read(19,rec=markal) (wf(i),i=i1,nmaxln)
c
c       reconstruct derivative of the wave functions
c
        do iso=1,2
          isnpad=(iso-1)*kjmx
          do js=1,jsmax(iso)
            nmax=npla(1,js+isnpad)
            if(nmax.ge.1) then
              ll=js/2
              sig=float(2*mod(ll,2)-1)
              do n=1,nmax
                ia=npla(n+1,js+isnpad)
                call derv15(wf(ia),wfs(ia),sig,ngrid,rstep)
                IF(thfbfold) CALL Gauss_folding_m(wf(ia),wffold(ia),
     &                       rstep,ngrid,js,hfbwidth,MAXVAL(jsmax(1:2)))
              enddo
            endif
          enddo
        enddo
c
c       test nucleus from file
c
        if(kfileb*nprot+nneut.ne.ntest) then
         npt=ntest/kfileb
         nnt=mod(ntest,kfileb)
         write(06,9299) nprot,nneut,npt,nnt
 9299    format(1h0,'-----> file ftn19 did not match: nprot,nneut=',
     *        2i4,'   but from file nprot,nneut=',2i4)
         stop "file ftn19 did not match"
        endif
c
c       give notice of succesful reading
c
        ifread = 1
        IF(TTEST) write(06,'(a,2i4)')
     &     'BUFFIO: succesful read: nprot,nneut=',nprot,nneut
c
       endif
c
      elseif(direc.eq.'o') then
c
c-------------------------------------------------------------------------------
c
c     write result on scratch file 19.
c       writing goes exactly the same way as reading.
c       the file-table is not checked anew. the last value for 'markfl',
c       the pointer to the first record, is taken as it is, assuming
c       that the same nucleus is now written as was read in just before.
c
c-------------------------------------------------------------------------------
c
       if(markfl.gt.0 .and. idisc.eq.1) then
        markal = markfl
        write(19,rec=markal) rstep,ngrid,
     &                      nprot,nneut,nmaxp,nmaxn,
     &                      weight,v2occ,deltaf,spenrg,
     &                      efen,efep,deltp,deltn,xlamp2,xlamn2,
     &                      jsmax,ifpair,npla
        nmaxln = nmaxto*ngrid
        nmaxbl = nmaxln/kwrdln
        nkspm  = kwrdln-1
        do i1=1,nmaxln-kwrdln,kwrdln
          markal = 1+markal
          write(19,rec=markal) (wf(i),i=i1,i1+nkspm)
        enddo
        i1     = nmaxbl*kwrdln+1
        markal = 1+markal
        write(19,rec=markal) (wf(i),i=i1,nmaxln)
        markfl = 2+nmaxbl+markfl
        maxfl  = max(markfl,maxfl)
        IF(TTEST) write(06,'(a,2i4)')
     &     'BUFFIO: succesful write: nprot,nneut=',nprot,nneut
       endif
      else
       stop ' inappropriate option in buffio'
      endif
c
      return
      end

      include 'matin2old.f'

c-----gaufol------------ part of the hartree-fock package ----------------------
c
      subroutine gaufol(rhoin,rhoout,width)
c
      include 'params.inc'
c
c     Gaussian folding using Fourier-Bessel transformation:
c      rhoin   = input density
c      rhoout  = folded output density
c      width   = width of folding (r.m.s. of Gaussian)
c     The grid parameters are transferred via commons.
c
      dimension rhoin(kg1),rhoout(kg1)
c
      include 'grid.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'constn.inc'
c
      data itest/0/
c
c
c----------------------------------------------------------------
c
c     preparation of auxiliary parameters
c
      anumb  = one*(nprot+nneut)
c      iqmax  = ngrid-1
c      rleng  = (ngrid-1)*rstep
c      delq   = pi/rleng
c      frpirs = 4*pi*rstep
      xlambd = 1.5/(width*width)
      qlambd = 0.25/xlambd
c
c     fourier-bessel transformation forward
c
      sum    = zero
      r      = zero
      do ir=2,ngrid
        r     = rstep+r
        sum   = r*r*rhoin(ir)+sum
      enddo
      aux1(1) = frpirs*sum
c
      do iq=2,iqmax
        q1     = iq*delq-delq
        dc     = cos(q1*rstep)
        ds     = sin(q1*rstep)
        r      = zero
        sum   = zero
        c      = one
        s      = zero
        do i=2,ngrid
          r      = rstep+r
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
          sum    = r*s*rhoin(i)+sum
        enddo
        sum        = frpirs*sum/q1
        aux1(iq)   = exp(-qlambd*(q1*q1))*sum
      enddo
c
c     fourier-bessel transformation backward
c
      dpi3   = half/(pi*pi)
c                                                         r=0.0 separately
      sum    = zero
      do iq=2,iqmax
        qa     = delq*iq-delq
        sum    = qa*qa*aux1(iq)+sum
      enddo
      rhoout(1) = abs(sum*delq*dpi3)
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
          sum    = qa*s*aux1(iq)+sum
        enddo
        rhoout(ir) = abs(sum*delq*dpi3/r)
c
      enddo
c
c     optionally test
c
      if(itest.eq.1) then
        sum1   = zero
        sum2   = zero
        sumr1  = zero
        sumr2  = zero
        r      = zero
        do ir=2,ngrid
          r     = rstep+r
          r2    = r*r
          sum1  = r2*rhoin(ir)+sum1
          sum2  = r2*rhoout(ir)+sum2
          r2    = r2*r2
          sumr1 = r2*rhoin(ir)+sumr1
          sumr2 = r2*rhoout(ir)+sumr2
        enddo
        sum1 = frpirs*sum1
        sum2 = frpirs*sum2
        sumr1 = frpirs*sumr1/sum1
        sumr2 = frpirs*sumr2/sum2
        write(6,'(a,2g12.4,a,g12.4)')
     &  ' rms in-out=',sqrt(sumr1),sqrt(sumr2),' diff. r^2=',sumr2-sumr1
      endif
c
      return
      end
      
cmomdenc-----momden------------ part of the hartree-fock package ----------------------
cmomdenc
cmomden      subroutine momden
cmomdenc
cmomden      include 'params.inc'
cmomdenc
cmomdenc     * * * * * * * * *  density in momentum space   * * * * * * * * * *
cmomdenc
cmomdenc     accumulates the momentum density and prints it.
cmomdenc
cmomdenc     for this task the routine needs some information on the shell
cmomdenc     structure; this is provided by the following fields:
cmomdenc      jsp   = total and orbital angular momentum for the state 'nshl'.
cmomdenc              (for a more detailed explanation see 'initgr')
cmomdenc      zwj   = "nr. of states in shell with combined ang.mom. 'js'",
cmomdenc              i.e. (2*j+1)/4*pi where j is total angular momentum.
cmomdenc      wls   = strength of spin-orbit splitting for ang.mom. 'js'
cmomdenc     (note that 'js' and also 'jsp' is a combined quantum number
cmomdenc      for total and orbital ang.mom.)
cmomdenc
cmomden      include 'skypar.inc'
cmomden      include 'grid.inc'
cmomden      include 'numpar.inc'
cmomden      include 'wavfun.inc'
cmomden      include 'constn.inc'
cmomdenc
cmomdenc     provide space for accumulation
cmomdenc
cmomden      dimension auxmom(3*kg1)
cmomden      equivalence (aux1,auxmom)
cmomdenc
cmomdenc      dimension zwj(19)
cmomdenc      data zwj/.15915494,.15915494,.31830989,.31830989,.47746483,
cmomdenc     &         .47746483,.63661977,.63661977,.79577472,.79577472,
cmomdenc     &         .95492966,.95492966,1.1140846,1.1140846,1.2732395,
cmomdenc     &         1.2732395,1.4232945,1.4232945,1.5915494/
cmomdenc
cmomdenc-------------------------------------------------------------------------------
cmomdenc
cmomden      delqmd = 0.1
cmomden      nkqmax = 50
cmomden      if(nkqmax.gt.3*kg1) stop ' too large nlqmax in momden'
cmomdenc
cmomdenc     iso-loop, initialise "shift constants" 'im,ia,ie' for the
cmomdenc     densities, 'ngw' for the weights, and 'nadd' for the
cmomdenc     wave functions.
cmomdenc
cmomden      do 20 iso=1,2
cmomden      ie     = ngrid*iso
cmomden      im     = ie-ngrid+1
cmomden      ia     = im+1
cmomden      if(iso.eq.1) then
cmomden       nmax   = nmaxp
cmomden      else
cmomden       nmax   = nmaxn
cmomden      endif
cmomden      ngw    = (iso-1)*nmaxp
cmomden      nadd   = ngw*ngrid-ngrid
cmomdenc
cmomdenc     zero fields before accumulation
cmomdenc
cmomden      do 10 i=1,nkqmax
cmomden      auxmom(i) = zero
cmomden   10 continue
cmomdenc
cmomdenc     loop over states for given isospin.
cmomdenc      'nshl' counts the states in harmonic oscillator ordering.
cmomdenc      'js' is the according composed ang.mom.,
cmomdenc      'll' is the orbital ang.mom.,
cmomdenc      'wl' is the spin-orbit splitting,
cmomdenc      'om' is the actual phase-space factor
cmomdenc
cmomden      do 30 nshl=1,nmax
cmomden      js     = jsp(nshl)
cmomden      ll     = js/2
cmomden      om     = zwj(js)*weight(nshl+ngw)/17.4532
cmomden      nadd   = nadd+ngrid
cmomden      if(om.lt.xsmall) goto 30
cmomden      write(6,'(a,2i4,a,g12.5)')
cmomden     & ' momentum density for state iso,na=',iso,nshl
cmomdenc      sumtst = zero
cmomdencc      iqmax  = kqmax
cmomdenc                                                 q=0.0 separately
cmomden      if(ll.eq.0) then
cmomden       sum    = zero
cmomden       r      = zero
cmomden       do 500 ir=2,ngrid
cmomden       r     = rstep+r
cmomden       sum   = r*wf(nadd+ir)+sum
cmomden  500  continue
cmomdenc       sumtst  = sum*sum + sumtst
cmomden       write(6,'(1x,f8.3,1pg12.4)') 0.0,sum*sum
cmomden       auxmom(1) = om*sum*sum + auxmom(1)
cmomden      endif
cmomdenc                                                      all other q
cmomden      do 530 iq=2,nkqmax
cmomden      q1     = iq*delqmd-delqmd
cmomdenc      dc     = cos(q1*rstep)
cmomdenc      ds     = sin(q1*rstep)
cmomden      r      = zero
cmomden      sum    = zero
cmomdenc      c      = one
cmomdenc      s      = zero
cmomden       do 510 i=2,ngrid
cmomden       r      = rstep+r
cmomdenc       x      = c*dc-s*ds
cmomdenc       s      = s*dc+c*ds
cmomdenc       c      = x
cmomden       sum    = r*bessel(r*q1,ll)*wf(nadd+i)+sum
cmomden  510  continue
cmomdenc      sumtst   = (sum*q1)**2 + sumtst
cmomden      write(6,'(1x,f8.3,1pg12.4)') q1,sum*sum
cmomden      auxmom(iq) = om*sum*sum + auxmom(iq)
cmomden  530 continue
cmomdenc      write(6,'(a,2i4,a,g12.5)')
cmomdenc     & ' state iso,na=',iso,nshl,': norm=',sumtst
cmomdenc
cmomden   30 continue
cmomdenc
cmomdenc     now print
cmomdenc
cmomden      write(6,'(/a,i2,a/a,66(/f9.3,g12.4))')
cmomden     & ' momentum density for iso=',iso,':',
cmomden     & '     q     density',
cmomden     & (delqmd*(iq-1),auxmom(iq),iq=1,nkqmax)
cmomdenc
cmomdenc     test total particle munber
cmomdenc
cmomden      summom = zero
cmomden      do 590 iq=1,nkqmax
cmomden      q1     = iq*delqmd-delqmd
cmomden      summom = auxmom(iq)*q1*q1 + summom
cmomden  590 continue
cmomden      write(6,'(a,g12.5)')
cmomden     & ' test total density: ',fourpi*delqmd*summom
cmomdenc
cmomden   20 continue
cmomdenc
cmomdenc     test bessel
cmomdenc
cmomdenc      do 310 ll=0,2
cmomdenc      do 310 ir=1,ngrid,2
cmomdenc      r     = rstep*(ir-1)
cmomdenc      sum   = zero
cmomdenc      do 320 iq=1,2*nkqmax
cmomdenc      q1     = iq*delqmd-delqmd
cmomdenc      sum   = (bessel(r*q1,ll)*q1)**2 + sum
cmomdenc  320 continue
cmomdenc      write(6,'(a,i3,2g12.5)')
cmomdenc     & ' test bessel: ll,r,bessel=',ll,r,fourpi*delqmd*sum*r*r
cmomdenc  310 continue
cmomdenc
cmomden      return
cmomden      end
c-----sixjot--------------------------------------------------------------------
c
c     this is a package for a.m.coupling coefficients from ikp mainz.
c
c     note: routines are changed to work with single precision calls
c           (but internal double precision) , 8.8.88 p.-g.r.
c
c received: by derrze0 via bitnet with njf for mpt218@derrze0; 2 aug 88 11:20:16
c date:     tue,  2 aug 88 11:16 n
c original_from:  vkpmza::bermuth
c from:     postmaster@dmznat51
c subject:  wigner und co.
c to:       mpt218@derrze0
c
c lieber prof. reinhard !
c achtung: die routinen sind sehr empfindlich auf real*8,
c d.h. bei eingabe von 0.5d0 funktionieren sie richtig,
c mit 0.5 kommt aber schrott raus.
c zur zeit schreibe ich bessere routinen. diese sind aber noch nicht fertig.
c                 viele gruesse
c                           klaus bermuth
c-----xinjot--------------------------------------------------------------------
c
      function xinjot (a,b,c,d,e,f,g,h,p)
c
c entnommen aus ibmlib1.for
c
c                       / a b c \
c     the 9-j symbol    | d e f |
c                       \ g h p /
c
      include 'params.inc'
c
      z=max(abs(a-p),abs(h-d),abs(b-f))
      zmax=min(a+p,h+d,b+f)+0.01
      sum=0.
    1 k=((z+z)+0.1)
      vz=1-2*mod(k,2)
      sum=sum+(z+z+1.)*sixjot(a,d,g,h,p,z)
     1 *sixjot(b,e,h,d,z,f)*sixjot(c,f,p,z,a,b)*vz
      z=z+1.
      if (zmax-z) 2,1,1
    2 xinjot=sum
      return
      end
c-----sixjot--------------------------------------------------------------------
c
      function  sixjot ( a,b,e,d,c,f)
c
c entnommen aus ibmlib1.for
c ueberarbeitet fuer external single precision, p.-g.r.
c
      include 'params.inc'
c
      double precision s,u,x,z
c
c     'sqr' stores the square roots of integer numbers
c
      dimension sqr(101),s1(47),s2(54)
c
c     'fak'  stores the log of facuulty of 'i'
c
      double precision fak(101),f1(47),f2(54)
c
      equivalence (fak(1),f1(1)),(fak(48),f2(1)),(sqr(1),s1(1)),
     &            (sqr(48),s2(1))
c
      data f1/0.d0,0.d0,
     &6.93147180559945d-01,1.79175946922806d 00,3.17805383034795d 00,
     &4.78749174278205d 00,6.57925121201010d 00,8.52516136106541d 00,
     &1.06046029027453d 01,1.28018274800815d 01,1.51044125730755d 01,
     &1.75023078458739d 01,1.99872144956619d 01,2.25521638531234d 01,
     &2.51912211827387d 01,2.78992713838409d 01,3.06718601060807d 01,
     &3.35050734501369d 01,3.63954452080331d 01,3.93398841871995d 01,
     &4.23356164607535d 01,4.53801388984769d 01,4.84711813518352d 01,
     &5.16066755677644d 01,5.47847293981123d 01,5.80036052229805d 01,
     &6.12617017610020d 01,6.45575386270063d 01,6.78897431371815d 01,
     &7.12570389671680d 01,7.46582363488302d 01,7.80922235533153d 01,
     &8.15579594561150d 01,8.50544670175815d 01,8.85808275421977d 01,
     &9.21361756036871d 01,9.57196945421432d 01,9.93306124547874d 01,
     &1.02968198614514d 02,1.06631760260643d 02,1.10320639714757d 02,
     &1.14034211781462d 02,1.17771881399745d 02,1.21533081515439d 02,
     &1.25317271149357d 02,1.29123933639127d 02,1.32952575035616d 02/
      data f2/
     &1.36802722637326d 02,1.40673923648234d 02,1.44565743946345d 02,
     &1.48477766951773d 02,1.52409592584497d 02,1.56360836303079d 02,
     &1.60331128216631d 02,1.64320112263195d 02,1.68327445448428d 02,
     &1.72352797139163d 02,1.76395848406997d 02,1.80456291417544d 02,
     &1.84533828861449d 02,1.88628173423672d 02,1.92739047287845d 02,
     &1.96866181672890d 02,2.01009316399282d 02,2.05168199482641d 02,
     &2.09342586752537d 02,2.13532241494563d 02,2.17736934113954d 02,
     &2.21956441819130d 02,2.26190548323728d 02,2.30439043565777d 02,
     &2.34701723442818d 02,2.38978389561834d 02,2.43268849002983d 02,
     &2.47572914096187d 02,2.51890402209723d 02,2.56221135550010d 02,
     &2.60564940971863d 02,2.64921649798553d 02,2.69291097651020d 02,
     &2.73673124285694d 02,2.78067573440366d 02,2.82474292687630d 02,
     &2.86893133295427d 02,2.91323950094270d 02,2.95766601350761d 02,
     &3.00220948647014d 02,3.04686856765669d 02,3.09164193580147d 02,
     &3.13652829949879d 02,3.18152639620209d 02,3.22663499126726d 02,
     &3.27185287703775d 02,3.31717887196928d 02,3.36261181979198d 02,
     &3.40815058870799d 02,3.45379407062267d 02,3.49954118040770d 02,
     &3.54539085519441d 02,3.59134205369575d 02,3.63739375555563d 02/
      data s1/1.00000000000000e 00,1.41421356237310e 00,
     &1.73205080756888e 00,2.00000000000000e 00,2.23606797749979e 00,
     &2.44948974278318e 00,2.64575131106459e 00,2.82842712474619e 00,
     &3.00000000000000e 00,3.16227766016838e 00,3.31662479035540e 00,
     &3.46410161513775e 00,3.60555127546399e 00,3.74165738677394e 00,
     &3.87298334620742e 00,4.00000000000000e 00,4.12310562561766e 00,
     &4.24264068711929e 00,4.35889894354067e 00,4.47213595499958e 00,
     &4.58257569495584e 00,4.69041575982343e 00,4.79583152331272e 00,
     &4.89897948556636e 00,5.00000000000000e 00,5.09901951359278e 00,
     &5.19615242270663e 00,5.29150262212918e 00,5.38516480713450e 00,
     &5.47722557505166e 00,5.56776436283002e 00,5.65685424949238e 00,
     &5.74456264653803e 00,5.83095189484530e 00,5.91607978309962e 00,
     &6.00000000000000e 00,6.08276253029822e 00,6.16441400296898e 00,
     &6.24499799839840e 00,6.32455532033676e 00,6.40312423743285e 00,
     &6.48074069840786e 00,6.55743852430200e 00,6.63324958071080e 00,
     &6.70820393249937e 00,6.78232998312527e 00,6.85565460040104e 00/
      data s2/
     &6.92820323027551e 00,7.00000000000000e 00,7.07106781186548e 00,
     &7.14142842854285e 00,7.21110255092798e 00,7.28010988928052e 00,
     &7.34846922834953e 00,7.41619848709566e 00,7.48331477354788e 00,
     &7.54983443527075e 00,7.61577310586391e 00,7.68114574786861e 00,
     &7.74596669241483e 00,7.81024967590665e 00,7.87400787401181e 00,
     &7.93725393319377e 00,8.00000000000000e 00,8.06225774829855e 00,
     &8.12403840463596e 00,8.18535277187245e 00,8.24621125123532e 00,
     &8.30662386291807e 00,8.36660026534076e 00,8.42614977317636e 00,
     &8.48528137423857e 00,8.54400374531753e 00,8.60232526704263e 00,
     &8.66025403784439e 00,8.71779788708135e 00,8.77496438739212e 00,
     &8.83176086632785e 00,8.88819441731559e 00,8.94427190999916e 00,
     &9.00000000000000e 00,9.05538513813742e 00,9.11043357914430e 00,
     &9.16515138991168e 00,9.21954445729289e 00,9.27361849549570e 00,
     &9.32737905308882e 00,9.38083151964686e 00,9.43398113205660e 00,
     &9.48683298050514e 00,9.53939201416946e 00,9.59166304662544e 00,
     &9.64365076099295e 00,9.69535971483266e 00,9.74679434480896e 00,
     &9.79795897113271e 00,9.84885780179610e 00,9.89949493661167e 00,
     &9.94987437106620e 00,1.00000000000000e 01,1.00498756211209e 01/
c
      data zero,one,two/0.0D0,1.0D0,2.0D0/
c
c     (-1.0)**xl as internal statement function
c
      phase(xl)=one-2.0D0*mod(abs(xl),two)
c
c     check triangular relations and integer matching of related a.m.
c
c           a - b - e
      if( e.gt.a+b .or. e.lt.abs(a-b) .or. mod(a+b+e,one).ne.zero) then
       sixjot = zero
       return
      endif
c           a - c - f
      if( f.gt.a+c .or. f.lt.abs(a-c) .or. mod(a+c+f,one).ne.zero) then
       sixjot = zero
       return
      endif
c           d - b - f
      if( f.gt.d+b .or. f.lt.abs(d-b) .or. mod(d+b+f,one).ne.zero) then
       sixjot = zero
       return
      endif
c           d - c - e
      if( e.gt.d+c .or. e.lt.abs(d-c) .or. mod(e+d+c,one).ne.zero) then
       sixjot = zero
       return
      endif
c
c     select simple cases for faster computation
c
      if (min(a,b,c,d,e,f)) 30,4,15
    4 if (min(d,c,f)) 30,5,10
    5 if (min(c,f)) 30,6,9
    6 if (f) 30,7,8
    7 j1=a+a+one
      j2=b+b+one
      sixjot = phase(a+b+e) / (sqr(j1)*sqr(j2))
      return
    8 j1=a+a+one
      j2=e+e+one
      sixjot = phase(a+b+e) / (sqr(j1)*sqr(j2))
      return
    9 j1=b+b+one
      j2=e+e+one
      sixjot = phase(a+b+e) / (sqr(j1)*sqr(j2))
      return
c
   10 if (min(a,b)) 30,12,11
   11 j1=a+a+one
      j2=c+c+one
      sixjot = phase(a+c+f) / (sqr(j1)*sqr(j2))
      return
   12 if (a) 30,13,14
   13 j1=c+c+one
      j2=e+e+one
      sixjot = phase(d+c+e) / (sqr(j1)*sqr(j2))
      return
   14 j1=d+d+one
      j2=e+e+one
      sixjot = phase(d+c+e) / (sqr(j1)*sqr(j2))
      return
c
c     all remaining cases
c
   15 continue
c
c     triangular forefactors
c
      z=0.0d0
      do 22 k = 1,4
      go to (16,17,18,19), k
   16 x1=a
      x2=b
      x3=e
      go to 20
   17 x1=c
      x2=d
      x3=e
      go to 20
   18 x1=a
      x2=c
      x3=f
      go to 20
   19 x1=b
      x2=d
      x3=f
   20 x=x1+x2+x3
      i1=x1+x2-x3
      i2=x1-x2+x3
      i3=x2+x3-x1
      if (max(i1,i2,i3)-100) 21,21,30
   21 i4=x
   22 z=z+fak(i1+1)+fak(i2+1)+fak(i3+1)-fak(i4+2)
      z=0.5*z
c
c     prepare intermedieate a.m. for the main sum of the 6-j symbol
c
      i1=a+b+e
      i2=c+d+e
      i3=a+c+f
      i4=b+d+f
      i5=a+b+c+d
      i6=a+d+e+f
      i7=b+c+e+f
      if (max(i1,i2,i3,i4,i5,i6,i7)-100) 23,23,30
c
c     main sum of the 6-j symbol
c
   23 k1= max(i1,i2,i3,i4)
      k2= min(i5,i6,i7)
      if (k1-k2) 24,24,30
   24 s=0.0d0
      u=1.0d0
      if(mod(k1,2)) 30,26,25
   25 u=-1.0d0
   26 do 29 i = k1,k2
      l1=i-i1
      l2=i-i2
      l3=i-i3
      l4=i-i4
      l5=i5-i
      l6=i6-i
      l7=i7-i
   27 x=fak(i+2)-fak(l1+1)-fak(l2+1)-fak(l3+1)-fak(l4+1)-fak(l5+1)
     &   -fak(l6+1)-fak(l7+1)+z
      if (x-85.0d0) 28,28,30
   28 s=s+u* dexp(x)
   29 u=-u
c
c     all done
c
      sixjot = s
      return
c
   30 sixjot = zero
      return
c
      end
c-----wigner--------------------------------------------------------------------
c
      function wigner(x1,x2,x3,y1,y2,y3)
c
c berechnet 3j-symbole
c ueberarbeitet aus ibmlib1.for von k.b.
c nochmals ueberarbeitet fuer external single precision, p.-g.r.
c
      include 'params.inc'
c
      double precision fctor,factsm,phas,vcc,done
c
c     'fact' stores the log of faculty of i
c
      double precision fact,facdat,facda2
      dimension fact(0:100),facdat(0:52),facda2(1:48)
      equivalence (fact(0),facdat(0)),(fact(53),facda2(1))
      data facdat/0.d0,0.d0,
     &6.93147180559945d-01,1.79175946922806d 00,3.17805383034795d 00,
     &4.78749174278205d 00,6.57925121201010d 00,8.52516136106541d 00,
     &1.06046029027453d 01,1.28018274800815d 01,1.51044125730755d 01,
     &1.75023078458739d 01,1.99872144956619d 01,2.25521638531234d 01,
     &2.51912211827387d 01,2.78992713838409d 01,3.06718601060807d 01,
     &3.35050734501369d 01,3.63954452080331d 01,3.93398841871995d 01,
     &4.23356164607535d 01,4.53801388984769d 01,4.84711813518352d 01,
     &5.16066755677644d 01,5.47847293981123d 01,5.80036052229805d 01,
     &6.12617017610020d 01,6.45575386270063d 01,6.78897431371815d 01,
     &7.12570389671680d 01,7.46582363488302d 01,7.80922235533153d 01,
     &8.15579594561150d 01,8.50544670175815d 01,8.85808275421977d 01,
     &9.21361756036871d 01,9.57196945421432d 01,9.93306124547874d 01,
     &1.02968198614514d 02,1.06631760260643d 02,1.10320639714757d 02,
     &1.14034211781462d 02,1.17771881399745d 02,1.21533081515439d 02,
     &1.25317271149357d 02,1.29123933639127d 02,1.32952575035616d 02,
     &1.36802722637326d 02,1.40673923648234d 02,1.44565743946345d 02,
     &1.48477766951773d 02,1.52409592584497d 02,1.56360836303079d 02/
      data facda2/1.60331128216631d 02,
     &1.64320112263195d 02,1.68327445448428d 02,
     &1.72352797139163d 02,1.76395848406997d 02,1.80456291417544d 02,
     &1.84533828861449d 02,1.88628173423672d 02,1.92739047287845d 02,
     &1.96866181672890d 02,2.01009316399282d 02,2.05168199482641d 02,
     &2.09342586752537d 02,2.13532241494563d 02,2.17736934113954d 02,
     &2.21956441819130d 02,2.26190548323728d 02,2.30439043565777d 02,
     &2.34701723442818d 02,2.38978389561834d 02,2.43268849002983d 02,
     &2.47572914096187d 02,2.51890402209723d 02,2.56221135550010d 02,
     &2.60564940971863d 02,2.64921649798553d 02,2.69291097651020d 02,
     &2.73673124285694d 02,2.78067573440366d 02,2.82474292687630d 02,
     &2.86893133295427d 02,2.91323950094270d 02,2.95766601350761d 02,
     &3.00220948647014d 02,3.04686856765669d 02,3.09164193580147d 02,
     &3.13652829949879d 02,3.18152639620209d 02,3.22663499126726d 02,
     &3.27185287703775d 02,3.31717887196928d 02,3.36261181979198d 02,
     &3.40815058870799d 02,3.45379407062267d 02,3.49954118040770d 02,
     &3.54539085519441d 02,3.59134205369575d 02,3.63739375555563d 02/
c
      data zero/0.0/
      data done/1.0d0/
c
c     phase of integer a.m. 'l' as internal statement function
c
      iphas(l)=1-2*mod(abs(l),2)
c
c     check consistency of m-quantum numbers
c
      msu=2*(y1+y2+y3)
      if(msu.ne.0) then
       wigner = zero
       return
      endif
c
      vcc=0.0d0
      jche=x1-x2-y3
      j1=x1*2
      j2=x2*2
      j3=x3*2
      m1=y1*2
      m2=y2*2
c
c     check consistency of j and according m quantum numbers
c
      if ((mod(j1+m1,2).ne.0) .or.
     &    (mod(j2+m2,2).ne.0)    ) then
          write(06,*) 'wigner-fehler',x1,x2,x3,y1,y2,y3
          write(10,*) 'wigner-fehler',x1,x2,x3,y1,y2,y3
          wigner = zero
          return
      endif
      if(j1.lt.iabs(m1)) then
       wigner = zero
       return
      endif
      if(j2.lt.iabs(m2)) then
       wigner = zero
       return
      endif
      if(j3.lt.iabs(m1+m2)) then
       wigner = zero
       return
      endif
c
c     arrange otpimal ordering
c
      if(j1.lt.j2) go to 20
      if(j3.lt.j2) go to 30
      icntr=0
      go to 40
   20 if(j3.lt.j1) go to 30
      icntr=-1
      it=j1
      j1=j2
      j2=it
      it=m1
      m1=m2
      m2=it
      go to 40
   30 icntr=1
      it=j2
      j2=j3
      j3=it
      m2=-m1-m2
   40 continue
c
c     composed angular momenta and check of triangular relations
c
      jz1=(j1+j2-j3)/2
      if(jz1.lt.0) then
       wigner = zero
       return
      endif
c
      jz2=(j1+j3-j2)/2
      if(jz2.lt.0) then
       wigner = zero
       return
      endif
c
      jz3=(j2+j3-j1)/2
      if(jz3.lt.0) then
       wigner = zero
       return
      endif
c
      jt1=(j1-j3+m2)/2
      jt2=(j2-j3-m1)/2
      numin=max(jt1,jt2,0)
      jt3=(j1-m1)/2
      jt4=(j2+m2)/2
      numax=min(jt3,jt4,jz1)
      jt5=(j2-m2)/2
      if(numax.lt.numin) then
       wigner = zero
       return
      endif
      j4=j1/2
      j5=j3/2
      phas=iphas(numin)
      do 100 nu=numin,numax
      factsm = fact(j4)+fact(j5)-fact(jt3-nu)-fact(nu-jt2)
     &        -fact(jt4-nu)-fact(nu-jt1)-fact(jz1-nu)-fact(nu)
      vcc = vcc+phas*dexp(factsm)
      phas=-phas
  100 continue
      fctor = fact((j1+m1)/2)+fact((j3+m1+m2)/2)+fact((j3-m1-m2)/2)
     &       +fact(jt3)+fact(jt4)+fact(jt5)+fact(jz1)+fact(jz2)
     &       +fact(jz3)
     &       -2.d0*fact(j4)-2.d0*fact(j5)-fact((j1+j2+j3)/2+1)
      fctor = (j3+1)*dexp(fctor)
      vcc=dsqrt(fctor)*vcc
      if(icntr)120,150,110
  110 vcc=vcc*dsqrt((j2+done)/(j3+done))*iphas(jt3)
      go to 150
  120 vcc=vcc*iphas(jz1)
  150 wigner=iphas(jche)*vcc/dsqrt(x3+x3+1.0d0)
      return
      end
c-----extrap------------------------------------------------------
c
      subroutine extrap(field,parity)
c
      include 'params.inc'
c
c     extrapolates the 'field' given at 'i'=2...ngrid  to its
c     value at 'i'=1 according to the 'parity'
c
      dimension field(kg1)
c
      data zero,four,dthree/0D0,4.0,0.333333333333/
      data trhalf,trfift,oneten/1.5,0.6,0.1/
c
      if(parity.gt.zero) then
c       five point quality
c      field(1) = (four*field(2)-field(3))*dthree
c       seven point quality
       field(1) = trhalf*field(2)-trfift*field(3)+oneten*field(4)
      else
       field(1) = zero
      endif
c
      return
      end
c
      include 'pairing.for'
      include 'exmass.f'
      include 'skfunc.f'
      include 'fyfunc.f'
c      include 'rmfekin.f'
      include 'pure_C.f'
      include 'cmact.f'
      include 'momden.f'
      include 'scanph.f'
      include 'asymm_nucmat.f'
      include 'evaluate_dip.f'
      include 'pnp_overlap.f'
      include 'sepenerg.f'
!      include 'test_hpsi.f'
!      include 'phtransmte.f'
