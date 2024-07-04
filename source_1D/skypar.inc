C INCLUDE BLOCK: SKYPAR
C
C      THE VARIABLES IN COMMON /SKYPAR/ MEAN:
C       T0-T4 & X0-X3  = PARAMETERS FOR THE SKYRME FORCE, STANDARD PART
C       T5-T8 & X5-X8  = PARAMETERS FOR THE EXTENDED VERSION
C                        (DENSITY DEPENDENT KINETIC TERMS)
c       tte,tto        = parameters of the tensor terms
C       TYUK  & XYUK   = PARAMETERS FOR A POSSIBLE FINITE RANGE TWO BODY FORCE
C       AYUK           = RANGE OF YUKAWA FOLDING IN FINITE RANGE FORCE
C       POWER          = POWER OF DENSITY DEPENDENCE
C       PAIRFP         = PARAMETER FOR PROTON PAIRING FORCE
C       PAIRFN         = PARAMETER FOR NEUTRON PAIRING FORCE
C       ACOUL          = FRACTIONAL CHARGE FOR PROTONS (USUALLY 1.0)
C       BETA20         = EMPIRICAL WIDTH OF QUADRUPOLE VIBRATIONS
C       T0A.....AYMS2  = FREQUENTLY OCCURING COMBINATIONS OF THE
C                        GIVEN FORCE PARAMETERS
c       b...           = new parameters of extended Skyrme forces
c       b11,b11p       = for "(nabla rho)**2 * tau" term
C       AMFC           = HBAR**2/(2*NUCLEON_MASS) FOR PROTON OR NEUTRON
c       rspln    = U^2*V^2 response energy for protons/neutrons,
c                  used for 'lambda_2' in Lipkin-Nogami
c       suv2r    = sum U^2*V^2 for protons/neutrons, to be used in Lipkin-Nogami
c       ilnrsp         = switch for residual interaction in Lipkin-Nogami,
c                        is to be set in BLOCK DATA
C       GPP/GPN        = STRENGTH OF PAIRING FORCE
c       pairfp/pairfn  = pairing strength in delta-force pairing
c       rho0pr         = reference density for surface pairing
C       DELTP/DELTN    = PAIRING GAP IN MEV   u*v weighted
c       deltp2/deltn2  = pairing gap  v**2 weighted
c       eqminp,eqminn  = min. quasiparticle energy
c       ecstab         = cut off parameter in stabilized pairing
c       alpinv         = parameter for Pade apprixmant density dep. 
c       betinv         = parameter for Pade apprixmant density dep. 
c       atajima        = exponent in Tajima's accelerated HFB step
c                        value is set in 'block data', typically 0.5
C       TEMPRT         = temperature,
C                        with 'ipair.eq.4'  --> pure thermal calculations,
C                        with 'ipair.ge.5'  --> BCS & temperature
c       rhoscp/n       = level density from GAVRG
c       rhoswm         = switch density in rho-dependent effective mass
c                        .eq.0 --> normal Skyrme model
c                        .ne.0 --> switch point in tanh-parametrization
c       exmass         = fictitous meson mass in explicit exchange functional
c       coefx          = strength of explicit exchange functional
c       ipair          = type of pairing:
C                        -1 = filling approximation
C                         0 = constant Delta approximation
C                         1 = constant Force approximation (BCS)
C                         2 = constant Force approximation with
c                             s.p. correction (BCS)
C                         3= Lipkin-Nogami
C                         4 = thermal HF, temperature handled by "DELTP"
C                         5 = BCS pairing with delat force (possibly T)
C       IFPAIR         = SIGNALS WHETHER THE PROTONS OR NEUTRONS NEED PAIRING
c       iparac         = actual version of IPAIR (changes during iteration!)
c       iparsh   = nr. of extra shells above minimum needed,
c                  0 -> magic is really magic, also useful for const. Gap
c                          and BCS (if strength is properly adjusted)
c                  1 -> recommended for Lipkin-Nogami, useful for BCS
c                  2 -> necessary in connection with delta-function pairing
c       xparsh   = real input for 'iparsh', used in case of energy occupation
c       cutfac   = factor determining dynamic cut off for pairing space
C       ISKY           = TYPE OF SKYRME FORCE (DERIVED FROM THE PARAMETERS)
C       IMFOLD   = SWITCH FOR DOUBLE-FOLDING OF KINETIC CURRENTS
c       ioptsk   = switch options for energy functional
c                  0 --> standard 
c                  1 --> IFB5IN = switch option in ATTACH yielding b5 term
c                  2 --> b0&b3-terms replaced by Pade approximant
c                  3 --> rational approximant for isoscalar and isovector term
c                        (isovector requires "exmass" = parameter(14) as input)
c                  4 --> RMF like kinetic&ls term (use 'b1kin' and 'b5kin')
c                  5 --> rho-dependent kinetic and l*s, direct version
c                  6 --> rho-dependent kinetic and l*s, inverse version
c                  7 --> BHF functional for rho-dependent terms, needs B1,B1P=0
c                  8 --> alternative isovector b3 and dens.dep. gradient
c                  9 --> extra density dependent term with own 'power'
c                 10 --> more flexible tensor terms ('riftm')
C
       REAL(8) :: T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,TYUK,tte,tto,
     *            X0,X1,X2,X3,X4,X5,X6,X7,X8,XYUK, AYUK, POWER,
     *            ayukp,ayukm,ayukls,
     *            ACOUL,  BETA20,
     *            T0A,T0B,T3A,T3B,T7A,T7B, T1P2,T1M2,T5P6,T5M6,
     *            T13P2,T13M2,T53P6,T53M6, T1M20,T1P2X,T5M60,T5P6X,
     *            B4,B4P,T8ZW,T3AS,T3BS,T53M6D,T5M60H,T5P6XH,    ! T3V,
     *            b3scal,b3eff,b3vect,b3peff,b3peff2,powerp,
     &            b5p,b6,b6p,
     *            TYA,TYB,AYMS2,
     &            rhoswm,exmass,coefx,b3dir,   b1kin,b5kin,
     &            b0,b1,b2,b3,b0p,b1p,b2p,b3p,
     &            c1,c1p,c4,c4p,bl2,bl2p,
     &            b3tilde,b11,b11p,b12,b12p, b7,b7p,pow7, 
     &            riftm,riftm1,
     *            AMFC(2),amfc0(2),
     &            rspln(2),suv2r(2),
     *            PAIRFP,PAIRFN,rho0pr,pptmax(2), ! pairing strength & cut dens.
     &            xparsh,cutfac,atajima,pairgm,
     &            hfbwidth,           ! folding width    
     &            ecstab,             ! cut-off energy for stabil. pairing
     &            alpinv,betinv,      ! parameters for pade-appr. densit.dep.
     &            alpin2,betin2,      ! parameters for pade-appr. densit.dep.
     *            GPP,GPN,DELTP,DELTN,deltp2,deltn2,eqminp,eqminn,
     &            temprt,                         ! temperature
     &            constr0,                        ! monpole force
     &            rhoscp,rhoscn,                  ! semiclassical level dens.
     &            rhosnm,              ! density symm.nucl.matter
     &            str_squark          ! strength s-quark --> weak charge
      LOGICAL thfb

       COMMON /SKYPAR/ T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,TYUK,tte,tto,
     *                 X0,X1,X2,X3,X4,X5,X6,X7,X8,XYUK, AYUK, POWER,
     *                 ayukp,ayukm,ayukls,
     *                 ACOUL,  BETA20,
     *                 T0A,T0B,T3A,T3B,T7A,T7B, T1P2,T1M2,T5P6,T5M6,
     *                 T13P2,T13M2,T53P6,T53M6, T1M20,T1P2X,T5M60,T5P6X,
     *                 B4,B4P,T8ZW,T3AS,T3BS,T53M6D,T5M60H,T5P6XH,    ! T3V,
     *                 b3scal,b3eff,b3vect,b3peff,b3peff2,powerp,
     &                 b5p,b6,b6p,
     *                 TYA,TYB,AYMS2,
     &                 rhoswm,exmass,coefx,b3dir,   b1kin,b5kin,
     &                 b0,b1,b2,b3,b0p,b1p,b2p,b3p,
     &                 c1,c1p,c4,c4p,bl2,bl2p,
     &                 b3tilde,b11,b11p,b12,b12p, b7,b7p,pow7, 
     &                 riftm,riftm1,
     *                 AMFC,amfc0,
     &                 rspln,suv2r,
     *                 PAIRFP,PAIRFN,rho0pr,pptmax, ! pairing strength & cut dens.
     &                 xparsh,cutfac,atajima,pairgm,
     &                 hfbwidth,           ! folding width    
     &                 ecstab,             ! cut-off energy for stabil. pairing
     &                 alpinv,betinv,      ! parameters for pade-appr. densit.dep.
     &                 alpin2,betin2,      ! parameters for pade-appr. densit.dep.
     *                 GPP,GPN,DELTP,DELTN,deltp2,deltn2,eqminp,eqminn,
     &                 temprt,                         ! temperature
     &                 constr0,                        ! monpole force
     &                 rhoscp,rhoscn,                  ! semiclassical level dens.
     &                 rhosnm,              ! density symm.nucl.matter
     &                 str_squark,          ! strength s-quark --> weak charge
     &                 numstp,numstn,
     &                 IPAIR,IFPAIR(2),ISKY,ilnrsp,IMFOLD,iparsh,iparac,
     &                 ioptsk,ipairsav,iteration_number,
     &                 thfb
c
c     switches for the various variants of the force
c      'texcha' = exchange contribution to coulomb force
c      'tyukaw' = finite range two body force with yukawa folding
c      'texten' = extended parametrisation ( using t5-t8)
c      'tpowerp'= different power in isovector part
c      'tkindir'= rho-dependent kinetic and l*s (itopsk=5 and 11)
c      'tkininv'= inv. rho-dependent kinetic and l*s (itopsk=6)
c      'talpinv'= rho dependence rho/(1+alpinv*rho)
c      'tpairFy'= pairing switch for Fayans functional
c      'tlipnog'= activates Lipkin-Nogami
c      'thfbfold'= folding of the apiring density in HFB
c
      logical tyukaw,texten,texcha,tssmas,tpair,texmas,tdirpa,trmfki,
     &        tbhfen,tb3til,tb11,tb12,tdens2,tpowerp,tkindir,tkininv,
     &        talpinv,talpin2,tdenslow,tpairFy,tlipnog,thfbfold
      common /tswitc/tyukaw,texten,texcha,tssmas,tpair,tpowerp,
     &               texmas,tdirpa,trmfki,tbhfen,tb3til,tb11,tb12,
     &               tdens2,tkindir,tkininv,talpinv,talpin2,tdenslow,
     &               tlipnog,thfbfold
c
c     Parameters for the Fayans functional
c
      REAL(8) ::avp,h1vp,h2vp,avm,h1vm,h2vm,sigmavol,rho0fy,
     &          asp,hsp,hgradp,asm,hsm,hgradm,hcoul,sigmacoul,
     &          fxi,h1xi,hgradxi,gammaxi,rho0xi,
     &          fxip,h1xip,hgradxip,fxin,h1xin,hgradxin,
     &          c0tau,c1tau,c0nabj,c1nabj,c0tens,c1tens,
     &          r0fy,r0fy2,rho0inv,efermfy,efacfy,ufacfy,pfacfy,
     &          clsp,clsm,ctensp,ctensm,ufacfypr,pgradfac
!      INTEGER,PARAMETER :: isovectormodel=1
      INTEGER,PARAMETER :: isovectormodel=2
!      INTEGER,PARAMETER :: isovectormodel=3
      common /faypar/avp,h1vp,h2vp,avm,h1vm,h2vm,sigmavol,rho0fy,
     &               asp,hsp,hgradp,asm,hsm,hgradm,hcoul,sigmacoul,
     &               fxi,h1xi,hgradxi,gammaxi,rho0xi,
     &               delfxi,delh1xi,delhgradxi,
     &               fxip,h1xip,hgradxip,fxin,h1xin,hgradxin,
     &               c0tau,c1tau,c0nabj,c1nabj,c0tens,c1tens,
     &               r0fy,r0fy2,rho0inv,efermfy,efacfy,ufacfy,pfacfy,
     &               clsp,clsm,ctensp,ctensm,ufacfypr,pgradfac,
     &               tpairFy
c
C END INCLUDE BLOCK





