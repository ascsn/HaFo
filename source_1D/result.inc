C INCLUDE BLOCK: RESULT
C
C      THE VARIABLES AND ARRAYS IN COMMON /RESULT/ MEAN:
C       RMS..   = ROOT-MEAN-SQUARE RADIUS OF THE ..-DENSITY;
C                 THE .. ARISE AS
C                   "BLANK"  = TOTAL NUCLEON DENSITY
C                   C        = CHARGE DENSITY WITH ALL FOLDINGS
C                   C0       = AS 'C' BUT WITHOUT MAGNETIC FORMFACTOR
C                              FOR THE NUCLEONS
C                   P        = PROTON DENSITY
C                   N        = NEUTRON DENSITY
C       RDMS..  = DIFFRACTION RADIUS DETERMINED FROM THE FIRST ZERO
C                 IN THE FORMFACTOR, FOR THE ..-DENSITY (AS ABOVE)
C       SURF..  = SURFACE WITDTH FOR THE ..-DENSITY (AS ABOVE)
C       FORMF   = ARRAY FOR THE CHARGE FORMFACTOR (I.E. IN MOMENTUM SPACE)
C       RHOCHR  = ARRAY FOR THE CHARGE DENSITY (IN COORDINATE SPACE)
C       ETOTAL  = TOTAL ENERGY
C       ETOT0   = TOTAL ENERGY BEFORE GROUND STATE CORRELATIONS
C       ESHL    = SHELL ENERGY FOR PROTONS OR NEUTRONS
C       ENREAR  = REARRANGEMENT ENERGY (FROM DENS.DEP. INTERACTION)
C       EPAIR   = PAIRING ENERGY
c       epp     = pairing energy for protons
c       epn     = pairing energy for neutrons
C       ZPE     = ZERO-POINT ENERGY
C       EKIN    = TOTAL KINETIC ENERGY
C       ECORTO  = CORRELATION ENERGY, STRAIGHTFORWARD WAY
C       ECORED  = CORRELATION ENERGY WITH DOUBLE COUNTING CORRECTION
C       EWID    = TOTAL WIDTH OF THE ENERGY
C       HPSIWID = ARRAY FOR THE ENERGY-FLUCTUATION OF THE S.P. STATES
C       SPENRG  = ARRAY FOR THE SINGLE PARTICLE ENERGIES
C
       COMMON /RESULT/ RMS,RMSC,RMSP,RMSN,RDMS,RDMSC,RDMSP,RDMSN,
     *                 SURF,SURFC,SURFP,SURFN,RDMSC0,SURFC0,SCORTO,
     &                 rmsw,rdmsw,surfw, rms4c,cmwid, apv_weak,
     &                 smom,smomc,smomp,smomn,
     &                 rsec,rsecn,rsecp,rsecc,rsec0,rhominmax,
     *                 FORMF(KFRM),RHOCHR(KG1),
     &                 foc(kfrm),fop(kfrm),fon(kfrm),fo(kfrm),
     &                 foc0(kfrm),fow(kfrm),epp0,epn0,sumuvp,sumuvn,
     *                 ETOTAL,ETOT0,ESHL(2),ENREAR,EPAIR,epp,epn,ZPE,
     *                 EKIN,epotsk,ekinpr,ekinne,pwidth,ecoul,enucmat,
     *                 ECORTO,ECORED,EWID,
     &                 drmst,detot,symcondi(2),symcondf(2),symmcond(2),
     *               HPSWID(KSL2),SPENRG(KSL2),spkin(KSL2),specm(KSL2),
     &                 delcm(KSL2),rmssp(ksl2)
C END INCLUDE BLOCK
