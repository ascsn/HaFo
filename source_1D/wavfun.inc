C INCLUDE BLOCK: WAVFUN
C
C     THE BLOCK CARRIES INFORMATION ON THE S.P. WAVE FUNCTIONS,
C     ON THE POTENTIALS AND ON THE VARIOUS DENSITIES.
C     THE ARRAYS MEAN:
C      WF      = WAVE FUNCTIONS, ORDERED IN SLICES OF LENGTH 'NGRID'
C                IN OSCILLATOR SHELL MODEL ORDERING. THE FIELD
C                POINTER 'NPLA' IS CARRIED IN COMMON /SHELLS/
C                ANF COMPUTED IN ROUTINE 'INITGR'.
C      WFS     = DERIVATIVE OF THE WAVE FUNCTIONS
C      BMASS   = INVERSE EFFECTIVE MASS
C      POT     = POTENTIAL
C      POTLS   = SPIN-ORBIT POTENTIAL
C      D1B     = WORK SPACE FOR DERIVATIVE OF EFFECTIVE MASS
C      RHO     = DENSITY
C      TAU     = KINETIC ENERGY DENSITY
C      SJ      = SPIN-ORBIT DENSITY
C      SJS     = DERIVATIVE OF SPIN-ORBIT DENSITY
C      D2R     = 2. DERIVATIVE OF DENSITY
C      D1R     = 1. DERIVATIVE OF DENSITY
C      AUX1-3  = WORKSPACE
C     ALL POTENTIALS AND DENSITIES STORE INFORMATION FOR PROTONS
C     IN THE PLACES [1..'NGRID'] AND FOR NEUTRONS IN THE PLACES
C     ['NGRID'+1...2*'NGRID'].
C
      COMMON /WAVFUN/ WF(KX11),WFS(KX11),
     *                BMASS(KGR2),POT(KGR2),POTLS(KGR2),D1B(KGR2),
     *                paircr(kgr2),       ! pair current for HFB
     &                pairpt(kgr2),       ! pair potential for HFB
     &                bmassl(kgr2),       ! inverse "mass" for L**2 term
     *                RHO(KGR2),
     *                TAU(KGR2),SJ(KGR2),SJS(KGR2),D2R(KGR2),D1R(KGR2),
     *                d2t(kgr2),d1t(kgr2),taul(kgr2),                 !new!
     *                AUX1(KGR2),AUX2(KGR2),AUX3(KGR2),
     &                ppadd(kgr2),regpairadd(kgr2)
C
C     THE VARIABLES AND ARRAYS IN COMMON /SHELL/ MEAN:
C      APROT    = PROTON NUMBER (AS REAL NUMBER)
C      ANEUT    = NEUTRON NUMBER (AS REAL NUMBER)
C      ANUMB    = TOTAL NUCLEON NUMBER (AS REAL NUMBER)
C      WEIGHT   = OCCUPATION WEIGHT FOR THE S.P. STATES
c      wgtobs   = weight for computing observables, differs from
c                 "WEIGHT" in case of Lipkin-Nogami pairing
c      v2occ    = the v**2 occupation coefficients from pairing
c                 (they differ from WEIGHT for the blocked state)
c      uvthrm   = thernmal*pairing occupation factor 'TTAN*u*V'
c      deltaf   = array of pairing gaps for the case of full pairing matrix
C      EFEP/N   = FERMI ENERGY FOR PROTONS OR NEUTRONS
C      WLS      = SPIN-ORBIT SPLITTING FOR THE S.P. STATES
C      zwj      = multiplicity factor/2*pi, for accumulating densities
c      pairmp/n = matrix workspace for full pairing and HFB
c      nsort    = pointer field for energy ordering (see 'prispe')
C      NPROT    = PROTON NUMBER (AS INTEGER NUMBER) ( IS INPUT PARAMETER )
C      NNEUT    = NEUTRON NUMBER (AS INTEGER NUMBER) ( IS INPUT PARAMETER )
C      NNUMB    = TOTAL NUCLEON NUMBER (AS INTEGER NUMBER)
C      NMAXP    = MAXIMUM NUMBER OF PROTON STATES ( IS INPUT PARAMETER )
C      NMAXN    = MAXIMUM NUMBER OF NEUTRON STATES ( IS INPUT PARAMETER )
C      NMINN    = POINTER TO FIRST NEUTRON STATE IN 'WEIGHT' ETC
C      NMAXTO   = TOTAL NUMBER OF STATES
c      nblock   = blocked state
c      isoodd   = isosping of blocked nucleon state
C      JSP      = COMBINED 'JS'-QUANTUM NUMBER OFR OSCILLATOR STATE NR. 'N'
C      JSMAX    = MAXIMUM 'JS' FOR PROTONS OR NEUTRONS
C      NPLA     = POINTER FIELD FOR WAVE FUNCTIONS (SEE ROUTINE 'INITGR')
c      xlamgap  = Lagrange parameter for linear gap constraint
c      xmugap   = Lagrange parameter for quadratic gap constraint
c      gapgoal  = wanted value of average gap (set =0 to switch off)
c      c0gap    = stepsize for update of Lagrangian parm. gap constraint
C
      COMMON /SHELLS/ APROT,ANEUT,ANUMB, EFEP,EFEN,xlamp2,xlamn2,
     &                efaddp,efaddn,
     &                xlamgap,xmugap,gapgoal,c0gap,
     &                WEIGHT(KSL2),wgtobs(ksl2),deltaf(ksl2),
     &                v2occ(KSL2),WLS(KJMX),zwj(KJMX),
     &                ph(KSL2),phpair(ksl2),ttthrm(ksl2),
     &                nsort(ksl2),
     *                JSP(KSL2),JSMAX(2),NPLA(KNMXP,KJMX2),
     *                NPROT,NNEUT,NNUMB,NMAXP,NMAXN,NMINN,NMAXTO,
     &                nblockp,nblockn,isoodd
cold     &                WLS(47),zwj(47),

       common /hfb/
     &                pairmn(khfbn,khfbn),pairmp(khfbp,khfbp)
C END INCLUDE BLOCK





