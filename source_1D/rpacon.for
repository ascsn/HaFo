C-----CORMOM-----------------------------------------------------
C
      SUBROUTINE CORMOM
C
      INCLUDE 'PARAMS.INC'
      PARAMETER(KRPMAT=KRPMAX*(KRPMAX+1)/2)
      PARAMETER (KMOMAX=20)
C
C     ESTIMATES VARIANCE OF MULTIPOLE FROM INTERPOLATING ODD MOMENTS.
C     THREE INTERPOLATING SCHEMES ARE IMPLEMENTED:
C      1. LAGUERRE INTERPOLATION
C      2. CUBIC SPLINE INTERPOLATION
C      3. ADJUSTMENT OF ESPECIALLY ADAPTED FUNTIONAL
C              MOM_K = SUM_N B(EL)_N ENERG_N**K
C         WITH THE FREE PARAMETERS B(EL)_N AND ENERG_N TO THE
C         GIVEN MOMENTS BY A LEAST SQUARES FIT. THE DESIRED
C         VARIANCE IS THEN MOM_0.
C     SPLINE AND LAGUERRE INTERPOLATION FAIL IN CRITICAL SPECTRA
C     WHICH CONTAIN VERY DISTINCT AND REMOTE PEAKS. ONLY METHOD 3
C     GIVES RELIABLE ESTIMATES IN ALL CASES WITH ERRORS WITHIN 5%.
C
C     THE PROGRAM IS CALLED FROM 'HAFOSU' FOR 'IFGSC'=5.
C
C     THE FIRST TWO INPUT LINES READ FROM THIS PROGRAM ARE
C     THE SAME AS IN THE PROGRAM 'SUMRUL'. FOR A DETAILED EXPLANATION
C     SEE THE COMMENTS THERE.
C     THE THIRD INPUT LINE IS SPECIAL TO THIS ROUTINES. IT CONTAINS:
C      ISORPA  = ISOSPIN OF THE MULTIPOLE OPERATOR WHOSE VARIANCE
C                IS TO BE COMPUTED;
C                 1 --> ISOSCALAR
C                 2 --> ISOVECTOR
C                 3 --> PROTON
C      IWGMIN  = MINIMUM ENERGY WEIGHT OF INITIAL 'QLROUT'; THE LOWEST
C                MOMENT IS THEN K_MIN = 2*IWGMIN-1.
C      IWGMAX  = MAXIMUM ENERGY WEIGHT OF INITIAL 'QLROUT'; THE HIGHEST
C                MOMENT IS THEN K_MAX = 2*IWGMAX-1.
C     RECOMMENDED INPUT: IWGMIN=-2 & IWGMAX=4 ; LOWER VALUES GIVE TOO
C     LITTLE INFORMATION IN CRITICAL CASES; HIGHER VALUES ARE COSTLY AND
C     RUN INTO NUMERICAL DIFFICULTIES (ROUND-OFF, MAXIMUM ENERGY).
C
      INCLUDE 'SKYPAR.INC'
      INCLUDE 'WAVFUN.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'RESULT.INC'
      INCLUDE 'NUMPAR.INC'
      INCLUDE 'CONSTN.INC'
      INCLUDE 'LINRES.INC'
c
      common /mompar/nbel,nmomin,nmomax,nrmom,iwrte
C
      DIMENSION SAVMOM(KMOMAX)
      DIMENSION NPARSV(KRPMAX),NHOLSV(KRPMAX),ISOPHS(KRPMAX)
C
      CHARACTER*1 ICHR(62)
      CHARACTER*16 TYPMOD
      CHARACTER*12 JSELEC
      CHARACTER*50 TOPER
      CHARACTER*6 QOP,POP,QBOP,QROP,PROP
      DATA ICHR/'1','2','3','4','5','6','7','8','9','0',
     &          'A','B','C','D','E','F','G','H','I','J',
     &          'K','L','M','N','O','P','Q','R','S','T',
     &          'U','V','W','X','Y','Z',
     &          'a','b','c','d','e','f','g','h','i','j',
     &          'k','l','m','n','o','p','q','r','s','t',
     &          'u','v','w','x','y','z'/
      DATA     QOP,      POP,    QBOP,    QROP,    PROP
     &    /'QIN   ','P     ','[H,P ]','QRES  ','PRES  '/
      DATA IFPRIM,CUTMOD,RESFAC,CMFAC
     &    /     0,1.0E-5,   1.0,  0.0/
C
      DATA IFOP14,IFOP17/0,0/
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C
C      READ PARAMETERS FOR 'CORMOM' SIMILAR AS IN 'SUMRUL'
C
       READ(5,*)     LRPA,RESFAC,CMFAC,CUTMOD
       READ(5,*)     ITRPA,EPSRPA,DABRPA,IFANLR,IFIRPA,IFGRHO
C
C     OPEN FILE 14 FOR CONCENTRATED OUTPUT OF MATRIX ELEMENTS
C
      IF(IFOP14.EQ.0) THEN
       OPEN(UNIT=14,STATUS='UNKNOWN',FORM='FORMATTED',FILE='FOR014')
       IFOP14 = 1
      ENDIF
C
C     SUPPLY DERIVATIVE OF 'BMASS' ON 'D1B' FOR USE IN 'HPSI'
C     AND DERIVATIVE OF 'RHO' ON 'D1R'.
C
      CALL DERV17(BMASS(1),D1B(1),ONE,NGRID,RSTEP)
      CALL DERV17(BMASS(NGRID+1),D1B(NGRID+1),ONE,NGRID,RSTEP)
      CALL DERV17(RHO(1),D1R(1), ONE,NGRID,RSTEP)
      CALL DERV17(RHO(ngrid+1),D1R(ngrid+1), ONE,NGRID,RSTEP)
C
C     WRITE A HEADER
C
      WRITE(14,'(A/1X,2A/1X,2(A,I4),5X,2(A,I4))')
     &  '1',
     &  '  TITLE:',TITLE,
     &  'NUCLEUS:  NPROT=',NPROT,'  NNEUT=',NNEUT,
     &  ' SHELLS:  NMAXP=',NMAXP,'  NAMXN=',NMAXN
      WRITE(14,'(1X,A,5G13.6/7X,A,4G13.6/7X,A,G13.6)')
     &  'FORCE: T0-T4=',T0,T1,T2,T3,T4,
     &        ' X0-X3=',X0,X1,X2,X3,
     &        ' POWER=',POWER
      Write(14,'(7x,3(a,g13.6)/)')
     &        ' TYUK=',tyuk,' XYUK=',xyuk,' AYUK=',ayuk
      WRITE(14,'(9X,A,2(/A,3G12.5,5X,G12.5))')
     &  '   R_RMS      R_DIFFR     SURFACE          ENERGY',
     &  ' MASS:   ', RMS, RDMS, SURF,ETOTAL,
     &  ' CHARGE: ',RMSC,RDMSC,SURFC
C
C     WRITE OPTIONS AND PARAMETERS FOR RPA PART
C
      WRITE(6,'(/1X,A,I3,3(A,g10.3))')
     &   ' RPA WITH: LRPA=',LRPA,'  Res.Int.Fac=',RESFAC,
     &   '  CmInv.Fac=',CMFAC,'  CutMode=',CUTMOD
      WRITE(6,'(21X,A,I5,A,F7.3,A,G12.4)')
     &    '  ITRPA=',ITRPA,' EPSRPA=',EPSRPA,' DAPRPA=',DABRPA
      WRITE(14,'(/1X,A,I3,3(A,g10.3))')
     &   ' RPA WITH: LRPA=',LRPA,'  Res.Int.Fac=',RESFAC,
     &   '  CmInv.Fac=',CMFAC,'  CutMode=',CUTMOD
c
c      set parameters inherited from 'extrsp'
c
c       IFCMPL = 0
c       ENEXT  = ZERO
c       ENWID  = ZERO
C
C     Set the angular momentum parameters according to given 'LRPA'
C
      ALRPA=LRPA
      MLRPA=2*LRPA+1
      CALL SETRED
      if(lrpa.ne.0) then
       FACMOM = (FOURPI/MLRPA)**2
      else
       FACMOM = (FOURPI/5)**2
      endif
C
C     READ INPUT SPECIAL TO THIS ROUTINE 'CORMOM'
C     AND DEDUCE BOUNDS OF MOMENTS FROM IT:
C
      READ(5,*) ISORPA,IWGMIN,IWGMAX,dellam
      IF(MOD(IWGMIN,2).NE.0) STOP ' IWGMIN NOT EVEN'
      IF(MOD(IWGMAX,2).NE.0) STOP ' IWGMAX NOT EVEN'
      IF(IWGMIN.GT.-2) STOP ' IWGMIN TOO LARGE'
      IF(IWGMAX.LT.0) STOP ' IWGMAX TOO SMALL'
      MINMOM = 2*IWGMIN-1
      nmomin = minmom
      MAXMOM = 2*IWGMAX-1
      nmomax = maxmom
      NRMOM  = (MAXMOM-MINMOM)/2+1
      IF(NRMOM.GT.KMOMAX) STOP ' TOO MANY MOMENTS'
C
      WRITE(6,'(3(A,i5))')
     & ' RUNNING CORMOM WITH: ISORPA=',ISORPA,
     & ' IWGMIN=',IWGMIN,' IWGMAX=',IWGMAX
      WRITE(14,'(3(A,i5))')
     & ' RUNNING CORMOM WITH: ISORPA=',ISORPA,
     & ' IWGMIN=',IWGMIN,' IWGMAX=',IWGMAX
C
C     program is run at least twice:
C      with RESFAC=1.0 to compute correlated variances, and
C      with RESFAC=0.0 to compute uncorrelated variances; sometimes also
C      with RESFAC=DELLAM to compute linear variances (this option
C      is run for the correlation energy, but still invalid).
C
      if(dellam.eq.zero) then
       ireslm = 0
      else
       ireslm = -1
      endif
      DO 1010 IRESF=1,ireslm,-1
      if(iresf.ge.0) then
       RESACT = IRESF*RESFAC
      else
       RESACT = DELLAM*RESFAC
      endif
C
C     INITIALISE BUFFERING FOR RESPONSE WAVE FUNCTIONS
C
      CALL BUFFQL(QLROUT,'C','      ')
C
C     INITIALIZE COUNTERS FOR 'QINIT'
C
      IFLAGR = 0
      IMODE  = 0
C
C     INITIALIZE MODE-SELECTION PARAMETERS FRO 'QINIT'
C
      TYPMOD = 'MOM=0'
      JSELEC = 'NOJJ'
      ISPLIT = 0
      IPRSHL = 0
      IH0WGT = 0
      IFPRIM = 1
C
C     LOOP OVER MODES FROM H**IWGMIN TO H**IWGMAX:
C
      DO 30 IHWGT=IWGMIN,IWGMAX,2
C                                            FIRST OPERATOR IS COMPUTED WITH
C                                            INVERSE MOMENTUM P=H**-1 Q
C                                            TO ACCESS THE LOWEST MOMENT.
      IF(IHWGT.EQ.IWGMIN) THEN
       IPWGT  = -1
      ELSE
       IPWGT  = 1
      ENDIF
C                                            SETUP THE OPERATOR
      CALL QINIT(TYPMOD,JSELEC,ISPLIT,IHWGT,IH0WGT,IPWGT,
     &           IPRSHL,IMODE,IFPRIM,IFLAGR,XKQRD,RESACT,CMFAC,
     &           NPARSV,NHOLSV,ISOPHS,TOPER)
C                                            STORE PRELIMINARY Q-OPERATOR
      QOP(4:4 ) = ICHR(IMODE)
      CALL BUFFQL(QLROUT,'O',QOP)
C                                            COMPUTE POINTER NAMES
      QOP(4:4 ) = ICHR(IMODE)
      POP(2:2 ) = ICHR(IMODE)
      QBOP(5:5) = ICHR(IMODE)
C                                            SET LEVEL OF ANALYSIS IN 'MOMENT'
      IFANLY = IFANLR
C                                            COMPUTE & STORE P, [H,P]
      CALL MOMENT(OVQQ,OVQHQ,OVQP,OVPP,OVPQB,OVQH0Q,QQINV,WIDRPA,
     &            QRLPR,QRLNE,RESACT,CMFAC,
     &            ESHIFQ,ESHIFT,ESHIFI,
     &            QOP,POP,QBOP,IPWGT,IFANLY)
      write(6,'(a,g12.5)') ' <Q*Q>=',ovqq
      write(14,'(a,g12.5)') ' <Q*Q>=',ovqq
      CALL BUFFQL(PLROUT,'I','[H,Q] ')
      QOP(4:4 ) = ICHR(imode)
      CALL BUFFQL(QLROUT,'I',QOP)
C                                            SAVE DATA FOR CORR. ENERGY
      if(ihwgt.eq.0) then
       BCOLL  = OVQHQ
       CCOLL  = OVPQB/ovqp**2
      endif
C                                            MOMENTS FROM OVERLAPS OF
C                                            PRESENT AND PREVIOUS OPERATORS
      if(imode.eq.1) then
       CALL BUFFQL(savlr,'I',POP)
       savmom(1) = QPOVER(QLROUT,savlr)*FACMOM
       savmom(2) = QPOVER(QLROUT,PLROUT)*FACMOM
      else
       imom   = 2*imode
       savmom(imom) = QPOVER(QLROUT,PLROUT)*FACMOM
       QOP(4:4 ) = ICHR(imode-1)
       CALL BUFFQL(QLROUT,'I',QOP)
       savmom(imom-1) = QPOVER(QLROUT,PLROUT)*FACMOM
      endif
C
   30 CONTINUE
      MAXMOD = IMODE
C
C     COMPUTE SIMPLE ESTIMATE FROM M_1 AND M_-1 TO DEFINE TYPICAL
C     DIMENSIONS AND TO MEASURE THE MOMENTS IN THESE NATURAL UNITS.
C
      imom1 = (1-nmomin)/2+1
      XMOM1  = savmom(IMOM1)
      XMOM1M = savmom(IMOM1-1)
      XMOM0P = SQRT(XMOM1*XMOM1M)
      AVENER = SQRT(XMOM1/XMOM1M)
      WRITE(6,'(/a,g12.5/2(a,g12.5))')
     &  ' RUNNING MOMENTS FOR RESFAC=',RESACT,
     &  ' M1*M-1 ESTIMATE: MOM0=',XMOM0P,' Av.E=',AVENER
      WRITE(14,'(/a,g12.5/2(a,g12.5))')
     &  ' RUNNING MOMENTS FOR RESFAC=',RESACT,
     &  ' M1*M-1 ESTIMATE: MOM0=',XMOM0P,' Av.E=',AVENER
C
      DO 110 IMOM=1,NRMOM
      XMOM   = savmom(IMOM)
      MOM    = MINMOM+2*(IMOM-1)
      SAVMOM(IMOM) = LOG(XMOM/(XMOM0P*AVENER**MOM))
      WRITE(6,'(a,i5,a,g12.5,a,i5,a,g12.5)')
     &  ' MOMENTS AT IMOM=',IMOM,': XMOM=',XMOM,' MOM=',MOM,
     &  ' SAVMOM=',SAVMOM(IMOM)
      WRITE(14,'(a,i5,a,g12.5,a,i5,a,g12.5)')
     &  ' MOMENTS AT IMOM=',IMOM,': XMOM=',XMOM,' MOM=',MOM,
     &  ' SAVMOM=',SAVMOM(IMOM)
  110 CONTINUE
C
C     SET COUNTERS FOR INTERPOLATION
C
      XSPL     = IMOM1-0.5
      IMOM   = NRMOM-1
      ISTART = 1
      npoint = imom-istart+1
      nbel   = npoint/2
c
c     interpolation by least squares fit to B(El)-form.
c       using AUX1 as workspace here, further workspaces are
c       employed in subroutine MOMFIT.
c
      open(unit=2,form='formatted',status='unknown',file='con')
      write(2,'(a,i5)')  ' nbel=',nbel
      call momfit(AUX1,AUX1(nbel+1),savmom,smom0f,0,2,chisqr)
      write(14,'(a,2i5,2g12.5)')
     &  ' MOMFIT with: npoint,nbel,smom0f,chisqr=',
     &  npoint,nbel,smom0f,chisqr
      close(unit=2)
      if(chisqr.gt.ONE*0.01) write(14,'(20x,a,g12.5)')
     &     ' !! CHISQR=',chisqr
      XMOM0F = EXP(SMOM0F)*XMOM0P
      IF(IRESF.EQ.0) THEN
       VAR0   = XMOM0F
       BCOLL0 = BCOLL
       CCOLL0 = CCOLL
      ELSEIF(IRESF.EQ.1) THEN
       VARRPA = XMOM0F
       BCOLLR = BCOLL
       CCOLLR = CCOLL
      ELSE
       VARDL  = XMOM0F
       BCOLLD = BCOLL
       CCOLLD = CCOLL
      ENDIF
c
c     interpolation by Laguerre polynomials
c
      XMOMLG = XLAGIN(SAVMOM,ISTART,IMOM,XSPL)
      XMOM0  = EXP(XMOMLG)*XMOM0P
c
c     print results
c
      do 550 iwu=6,14,8
  550 write(iwu,'(/a,i2,a/17x,a,2(/a,2g12.5)/)')
     &  ' Interpolation for ',npoint,' given points in SAVMOM:',
     &  ' LAGUERRE     SPECIAL',
     &  ' LOG(REN.MOMENT)=',XMOMLG,SMOM0F,
     &  '          MOMENT=',XMOM0,XMOM0F
C
 1010 CONTINUE
C
C     COMPUTE CORRECTION ON SURFACE THICKNESS, BASED ON VARINCES
C     FROM FIT TO SPECIAL FUNCTIONAL.
C
      CALL SCLCOR(SMSCLP,SMSCLN,SMSCLT)
      CORRPA = VARRPA-VAR0
      write(14,'(a,3g12.5)')
     &  ' Scaling: SMSCLP,SMSCLN,SMSCLT=',SMSCLP,SMSCLN,SMSCLT,
     &  ' moments: VARRPA,  VAR0,CORRPA=',VARRPA,  VAR0,CORRPA
      DELSGP = ZERO
      DELSGS = ZERO
      IF(ISORPA.eq.1) THEN
       SCALT  = ANUMB*SMSCLT/FOURPI
       DELSGS = CORRPA/(SCALT*SCALT)
c      elseif(ISORPA.eq.2) then
c       write(14,'(a)')
c     &  ' coupling of isovector mode to DELSIG not implemented.'
      elseif(NPROT.GT.0) then
       SCALP  = APROT*SMSCLP/FOURPI
       DELSGP = CORRPA/(SCALP*SCALP)
      endif
      WRITE(14,'(2(A,G12.4))')
     &    ' OLD DELTA-SIGMA**2:  PROTON=',DELSGP,
     &    '  TOTAL=',DELSGS
      WRITE(6,'(2(A,G12.4))')
     &    ' OLD DELTA-SIGMA**2:  PROTON=',DELSGP,
     &    '  TOTAL=',DELSGS
      if(DELSGS.ne.ZERO) DELSGS = (CORRPA/VARRPA)*DELSGS
      if(DELSGP.ne.ZERO) DELSGP = (CORRPA/VARRPA)*DELSGP
      WRITE(14,'(2(A,G12.4))')
     &    ' NEW DELTA-SIGMA**2:  PROTON=',DELSGP,
     &    '  TOTAL=',DELSGS
      WRITE(6,'(2(A,G12.4))')
     &    ' NEW DELTA-SIGMA**2:  PROTON=',DELSGP,
     &    '  TOTAL=',DELSGS
C
C     COMPUTE CORRECTIONS TO ENERGY, assuming approximate redundance
C
      ifetst = 0
      if(ifetst.ne.1) return
c
      qsqrpa = varrpa/facmom
      qsq0   = var0/facmom
      ECORR  = 0.25*(CCOLLR*qsqrpa+0.25*BCOLLR/qsqrpa)
      ECOR0  = 0.25*(CCOLL0*qsq0+0.25*BCOLL0/qsq0)
      IF(DELLAM.NE.ZERO) THEN
       qsqdl  = vardl/facmom
       ECORD  = 0.25*(CCOLLD*qsqdl+0.25*BCOLLD/qsqdl)
      ELSE
       qsqdl  = zero
       ECORD  = zero
      ENDIF
      write(14,'(/a/a,3(/2x,F4.2,1x,g12.5,3g11.4))')
     & ' Energy Correlations:',
     & ' RESFAC     ECOR         B         C      <Q**2>',
     & RESFAC,ECORR,BCOLLR,CCOLLR,qsqrpa,
     & ZERO,ECOR0,BCOLL0,CCOLL0,qsq0,
     & DELLAM*RESFAC,ECORD,BCOLLD,CCOLLD,qsqdl
      ECOR   = -( ECORR-ECOR0)*MLRPA
      IF(DELLAM.NE.ZERO) THEN
       ECORED = -( ECORR-ECOR0 - (ECORD-ECOR0)/DELLAM**2 )*MLRPA
      ELSE
       ECORED = zero
      ENDIF
      write(14,'(a,g12.5,a,g12.5)')
     & ' Corr.Energy=',ECOR,' Reduced Corr.Energy=',ECORED
C
      RETURN
      END
C-----XLAGIN-----------------------------------------------------
C
      FUNCTION XLAGIN(FIELD,ISTART,IEND,X)
C
      INCLUDE 'PARAMS.INC'
C
C     LAGUERRE INTERPOLATION
C
      DIMENSION FIELD(IEND)
C
C------------------------------------------------------------------------------
C
      SUM   = 0.0
      DO 10 IFUNCT=ISTART,IEND
      FACT  = 1.0
      DO 20 IFACT=ISTART,IEND
      IF(IFACT.NE.IFUNCT) FACT  = FACT*(X-IFACT)/(IFUNCT-IFACT)
   20 CONTINUE
      SUM   = FIELD(IFUNCT)*FACT+SUM
   10 CONTINUE
      XLAGIN = SUM
C
      RETURN
      END
c-----momfit-------------------------------------------------------------------
c
      subroutine momfit(bel,energ,xmomen,xmom0,iweigh,IWU,CHISQR)
c
      INCLUDE 'PARAMS.INC'
      parameter (KPARMX=12, KDATMX=12, KTESTC=1)
c
c     the routine adjusts 'nbel' representative B(El) values 'bel(n)'
c     and according energies 'energ(n)' such that the given moments
c     'xmomen' are reproduced by
c
c        Moment_k = sum_n B(El)_n energ_n**k  ,   k = odd .
c
c
c     As one special result, the zeroth moment 'xmom0' is returned explcitely.
c
c     The list parameters are:
c      bel       = array(1..nbel) for fitted B(El)-values.
c      energ     = array(1..nbel) for fitted energies.
c      xmomen    = array(1..nrmom) for given odd moments from 'nmomin' til
c                  'nmomax' in steps of 2.
c      nbel      = number of representatives for B(El)-values and
c                  according energies.
c      xmom0     = zeroth moment extrapolated from the given odd moments
c                  via the representation in terms of B(El) and energies.
c     transfer via common:
c      nmomin    = smallest (odd) moment on 'xmomen'.
c      nrmom     = number of (odd) moments on 'xmomen' from 'nmomin' on.
c
      INCLUDE 'LINRES.INC'
      INCLUDE 'CONSTN.INC'
c
      dimension bel(nbel),energ(nbel)
      dimension xmomen(nrmom)
      DIMENSION Y(KDATMX),SIGMAY(KDATMX),YFIT(KDATMX),
     &          A(KPARMX),ARES(KPARMX),SIGMAA(KPARMX),SIGMA0(KPARMX),
     &          DPARMO(KPARMX)
      DIMENSION ALPHA(KPARMX,KPARMX),zarray(KPARMX,KPARMX),
     &          B(KPARMX),beta(KPARMX)
      DIMENSION PARMOM(KPARMX),CHIMOM(KPARMX),PARSAV(KPARMX,3),
     &          YFITMO(KDATMX),FDERIJ(KDATMX,KPARMX)
c
      dimension IRAND(30)
c
      common /mompar/nbel,nmomin,nmomax,nrmom,iwrte
c
c     access workspaces by equivalence
c
      EQUIVALENCE (Y,HLINRS),(SIGMAY,WORK1),(YFIT,WORK2),(A,WORK3),
     &            (ARES,QRAD),(SIGMAA,QCUR(1,1)),(SIGMA0,QCUR(1,2)),
     &            (DPARMO,DQCUR),(PARMOM,RQCUR),(CHIMOM,PRHO),
     &            (PARSAV,PSJ),(YFITMO,DPRHO),(B,D2PRHO),
     &            (BETA,PTAU),
     &            (ALPHA,QLROUT),(ZARRAY,PLROUT),(FDERIJ,SAVLR)
c
      data itmax,epsrel,epsabs,dampit/220,1.e-6,1.e-4,0.8/
c      data dparmo/KPARMX*0.01/
      data ISTCHM/5/
      DATA IRAND/1,1,1,0,1,1,0,1,0,1,1,0,1,0,0,1,1,0,1,1,1,0,0,
     *           1,1,0,1,0,1,0/
c
c------------------------------------------------------------------------------
c
c     check input
c
      if(KPARMX+ISTCHM.gt.30) stop ' field IRAND in MOMFIT too small'
      write(IWU,'(a,3i5/a,10(/1x,5g12.5))')
     &  ' momfit entered with: nbel,nmomin,nrmom=',nbel,nmomin,nrmom,
     &  ' xmomen=',xmomen
      if(mod(nmomin,2).eq.0) stop ' NMOMIN not odd in routine MOMFIT'
      nmomax = nmomin+2*(nrmom-1)
c      if(2*nbel.ne.nrmom) stop ' NRMOM should be NRMOM/2'
      nmomax = nmomin+2*(nrmom-1)
      if(nmomin.gt.-1) stop ' NMOMIN should be smaller than 0'
      if(nmomax.lt.1) stop ' NMOMAX should be larger than 0'
c
c     first guess from +-1 moment
c
      imom1  = (-nmomin)/2+2
      write(IWU,'(a,2i5)')
     &  ' imom1,nmomax=',imom1,nmomax
c      bel1   = sqrt(xmomen(imom1)*xmomen(imom1-1))
c      energ1 = sqrt(xmomen(imom1)/xmomen(imom1-1))
c      write(6,'(2(a,g12.5))')
c     &  ' first estimate: average energy=',energ1,
c     &  ' zeroth moment=',bel1
c
c     initialize 'bel' and 'energ' with first guess.
c
      xnbel  = nbel
      xnbel  = 1.0/sqrt(xnbel)
      estep  = 1.8/(nbel)
      do 10 ibel=1,nbel
      bel(ibel)    = xnbel
      energ(ibel)  = ((ibel-1-nbel/2)*estep+1.0)
      a(ibel)      = log(bel(ibel))
      a(ibel+nbel) = log(energ(ibel))
   10 continue
c
c     rescale 'xmomen' in natural units defined from first guess
c
      do 20 imom=1,nrmom
      kmomen = nmomin+2*(imom-1)
c      xmomen(imom) = xmomen(imom)/(bel1*energ1**kmomen)
      Y(imom)    = xmomen(imom)
      if(iweigh.eq.1) then
       dparmo(imom) = 0.01*sqrt(ONE*abs(kmomen))
      else
       dparmo(imom) = 0.01
      endif
   20 continue
c
c
c  iterate to full solution by least squares
c
c
      NWDR   = itmax
      ISTOCH = 0
      IFCURV = 20
      IFLAST = 0
      IFALT2 = 0
      NPARMX = 2*nbel
      NDATMX = nrmom
      X      = 0.0
      do 90 imom=1,nrmom
   90 SIGMAY(imom) = 1.0
      chisqr = 1.0e30
      IF(NWDR.LT.2) THEN
       IST1 =1
      ELSE
       IST1 =0
      ENDIF
      do 100 iter=1,itmax
      if(ist1.gt.0) write(6,'(/a,i4)')  ' Iteration Nr.',iter
      CHIA=CHISQR
      WRITE(IWU,'(A,I3)')
     &  ' Next Iteration starts now. Iteration Nr.=',ITER
      call FIT(A,dparmo,Y,YFIT,SIGMAY,SIGMAA,SIGMA0,
     &         NPARMX,NDATMX,IFALT2,NWDR,IWARN,IFLAST,IFCURV,
     &         chisqr,IWU,IST1,
     &         ALPHA,zarray,B,beta,PARMOM,CHIMOM,PARSAV,
     &         YFITMO,FDERIJ)
      chitst = chisqr
C
C      IWARN=99 INDICATES THAT THE GRADIENT STEP RUN INTO
C      BAD CONDITIONS. THE PARAMETERS ARE CHANGED STOCHASTICALLY
C      AND THE FIT IS CONTINUED.
C
       IF(IWARN.EQ.99) THEN
        DO 359 J=1,NPARMX
  359   A(J) = A(J) + (2*IRAND(ISTOCH+J)-1)*DPARMO(J)*2.0
        ISTOCH = 1+ISTOCH
        IFALT2 = 0
        CHITST=.1E30
        WRITE(IWU,'(/A/)') ' IWARN=99. STOCHASTIC CHANGE OF PARAMS.'
       ENDIF
       IF(ISTOCH.GT.ISTCHM) THEN
        WRITE(IWU,'(/A)') ' NUMBER OF STOCHASTIC ATTEMPTS EXHAUSTED\'
        GOTO 199
       ENDIF
      chimax = .1e30
      if((chia.ge.chitst .and. chia/chitst.lt.1.+epsrel
     &    .and. chitst.lt.chimax)
     &   .or. chitst.lt.epsabs) goto 199
c
  100 continue
c
  199 continue
      if(iter.lt.itmax) write(IWU,'(a,i5,a,g12.5)')
     &  ' Convergence after ',iter,' iterations. CHISQR=',chisqr
c
c     unscale 'xmomen' back to original units
c
c      do 30 imom=1,nrmom
c      kmomen = nmomin+2*(imom-1)
c      xmomen(imom) = xmomen(imom)*(bel1*energ1**kmomen)
c   30 continue
c      do 40 ibel=1,nbel
c      bel(ibel) = bel1*bel(ibel)
c      energ(ibel) = energ1*energ(ibel)
c   40 continue
c
c     compute zeroth moment
c
      sum    = 0.0
      do 50 ibel=1,nbel
      if(a(ibel).lt.-50.0*ONE) then
       expmom = 0.0
      else
       expmom = exp(a(ibel))
      endif
   50 sum    = expmom+sum
      xmom0  = log(sum)
      write(IWU,'(a,g12.5)')
     &  ' final zeroth moment=',xmom0
c
      return
      end
C-----FIT   ---------------------------------------------------------
C
      SUBROUTINE FIT(A,dparmo,Y,YFIT,SIGMAY,SIGMAA,SIGMA0,
     &               NPARMX,NDATMX,IFALT2,NWDR,IWARN,IFLAST,IFCURV,
     &               chisqr,IWU,IST1,
     &               ALPHA,zarray,B,beta,PARMOM,CHIMOM,PARSAV,
     &               YFITMO,FDERIJ)
C
C     CURFIT FROM BEVINGTON , COMPLEMENTED BY A CURVED STEP.
C     THIS VERSION IS ADAPTED TO THE RMF-PACKAGE.
C
C      THE PARAMETERS IN THE CALLING SEQUENCE ARE:
C       A       = INITIAL PARAMETERS AS INPUT, NEW PARAMETERS AFTER
C                 ONE ITERATION STEP GIVEN BACK AS OUTPUT.
C       DPARMO  = VARIATION OF PARAMETERS FOR DERIVATIVES
C       Y       = FIELD CONTAINING THE EXPERIMENTAL DATA,.
C       YFIT    = FIELD CONTAINING THE THEORETICAL RESULTS IN THE
C                  SAME MANNER AS THE EXPERIMENTAL DATA ON Y.
C       SIGMAY  = FIELD CONTAINING THE CHOOSEN UNCERTAINTY ON THE
C                 EXPERIMENTAL DATA.
C       SIGMAA  = FIELD CONTAINING THE CORRELATED ERRORS ON THE
C                  MODEL PARAMETERS, IT IS A 1-DIMENSIONAL FIELD
C                  SIGMAA(PARAMETER).
C       SIGMA0  = FIELD CONTAINING THE UNCORRELATED ERRORS.
C       NPARMX  = NUMBER OF PARAMETERS (ON 'A') TO BE FITTED
C       NDATMX  = NUMBER OF DATA POINTS IN THE FIT
C       IFALT2  = SWITCH FOR INTERNAL INITIALISATION.
C                 IS EQUAL TO 0 IF THE
C                 ROUTINE IS CALLED THE FIRST TIME AT THE
C                 START OF AN ITERATION. IS THEN RETURNES WITH
C                 IFALT2=1 WHICH INDICATES THAT EVERYTHING IS
C                 PROPERLY INITIALISED NOW.
C       NWDR    = MAXIMUM NUMBER OF ITERATIONS; USED HERE ONLY TO
C                 CUT-OFF CALCULATION OF ERROR MATRIX IF NWDR.LE.0  .
C       IWARN   = RETURNS THE NUMBER OF ATTEMPTS TO FIND A LOWER
C                 CHI**2 . IWARN=1 MEANS A DIRECT SUCCES. IWARN=99
C                 MEANS THAT NO LOWER CHI**2 HAS BEEN FOUND.
C       IFLAST  = 1 TO INDICATE THE FINAL CALL WHICH PRODUCES A
C                   DETAILED OUTPUT, 0 IN ALL OTHER CASES.
C       IFCURV  = NUMBER OF NORMAL BEVINGTON STEPS WHICH ARE TO BE
C                 RUN BEFORE THE CURVED STEP IS SWITCHED ON. IF NO
C                 CURVED STEP IS DESIRED MAKE IFCURVE LARGER THAN
C                 NWDR.
C       CHISQR  = CHI**2 AFTER ONE STEP OF 'FIT'
C       IWU     = WRITE UNIT FOR CONTROL OUTPUT
C       IST1    = SWITCHES PRINT LEVEL IN THE ROUTINE FIT.
C
C      THE FOLLOWING LIST PARAMETERS ARE WORKSPACES WHOSE SIZE IS TO
C      BE READ OFF FROM THE SUBSEQUENT DIMENSION STATEMENT:
C       ALPHA, ZARRAY, B, BETA, PARMOM, CHIMOM, PARSAV, FDERIJ,
C       YFITM0   = AUXILIARY FIELD FOR INTERMEDIATE SAVING AND
C                  TRANSFER OF THEORETICAL RESULTS TO FIELD YFIT.
C
C      THE FOLLOWING PARAMETERS ARE OPTIONAL SWITCHES FOR THE ROUTINES
C      CALLED FROM 'FIT'. THEY MAY BE TRANSFERRED VIA AN APPROPRIATE
C      COMMON BLOCK:
C       IFSAVE  = SWITCHES THE SAVING OF HARTREE-FOCK INFORMATION
C                  IN THE HAFO-PACKAGE. IFSAVE=1 MEANS SAVING ON
C                  FILE 19. IT IS SWITCHED ON ONLY FOR THE "CENTRAL"
C                  PARAMETER SETS AND SWITCHED OF DURING THE
C                  VARIATION OF THE PARAMETERS.
C       IWR     = SWITCHES PRINT LEVEL IN THE HAFO-PACKAGE.
C       IWRTE   = SWITCHES PRINT LEVEL IN THE ROUTINE CHIFCT.
C
C     DIMENSION STATEMENT VALID FOR NPARMX UP TO KPARMX
C
C
      INCLUDE 'PARAMS.INC'
      parameter (KPARMX=12, KDATMX=12, KTESTC=1)
C
      CHARACTER*4 PARBRA
      CHARACTER*2 PARKET
C     SAVE PARSAV,CHISAV,ITERA
      DIMENSION Y(KDATMX),SIGMAY(KDATMX),YFIT(KDATMX),
     &          A(KPARMX),ARES(KPARMX),SIGMAA(KPARMX),SIGMA0(KPARMX),
     &          DPARMO(KPARMX)
c      DIMENSION ALPHA(KPARMX,KPARMX),zarray(KPARMX,KPARMX),
c     &          DERIV(KPARMX),B(KPARMX),beta(KPARMX)
      DIMENSION ALPHA(KPARMX,KPARMX),zarray(KPARMX,KPARMX),
     &          B(KPARMX),beta(KPARMX)
      DIMENSION PARMOM(KPARMX),CHIMOM(KPARMX),PARSAV(KPARMX,3),
     &          CHISAV(3),
     &          YFITMO(KDATMX),FDERIJ(KDATMX,KPARMX)
      DIMENSION CHIS00(KTESTC),chisu0(KTESTC),SUMCHI(KTESTC)
C
c      INCLUDE 'NUMPAR.INC'
c      INCLUDE 'COMFIT.INC'
c
      common /mompar/nbel,nmomin,nmomax,nrmom,iwrte
C
      DATA NDUMMY,IREDUC,ITRMAX,IFREPT
     *    /   0  ,   1  ,   3  ,   0  /
C
C     THESE PARAMETERS ARE OPTIONS FOR TESTING ETC :
C       IREDUC   = 0 --> DAMPING OF STEP BY INCREASING 'FLAMDA' AS
C                        PROPSED BY BEVINGTON
C                  1 --> DAMPING OF STEP BY REDUCING IT BY FACTOR 2
C       ITRMAX   = NUMBER OF ATTEMPTS WITH REDUCED STEP IN CURVED STEP
C       IFREPT   = 0 --> CHI**2 AND SITUATION FROM PREVIOUS FINAL
C                        CALCULATION IS USED
C                  1 --> NEW CALCULATION AT THE BEGINNING
C                  2 --> EQUILIBRISED STEP
C
      DATA PARBRA,PARKET/'PAR(',')='/
      DATA XLOWLM/1.E-30/
      data one/1.0/
C
C--------------------------------------------------------------------
C
      IF(IFALT2.EQ.0) THEN
       ITERA = 1
      ELSE
       ITERA = 1+ITERA
      ENDIF
       WRITE(IWU,'(A,I5,A,I5)') '0 FIT ENTERED WITH IFALT2=',
     *                   IFALT2,'. NOW ITERA=',ITERA
       WRITE(IWU,'(3(A,I3))') ' RUNNING WITH THE OPTIONS: IREDUC=',
     *            IREDUC,'     IFREPT=',IFREPT
      NFREE  = NDATMX - NPARMX
      RFREE  = NFREE+1.e-10
      CHIIN  = CHISQR
      IF(NFREE.lt.0 .AND. IFLAST.NE.1) THEN
       WRITE(IWU,'(2A)') ' NOT ENOUGH DATA FOR FITTING THE PARAMETERS.',
     *                  ' ====> STOP \'
       STOP
      ENDIF
      if(ist1.ne.0)  write(6,'(a,9g12.5)')
     &  ' FIT for  YFIT=',(Y(imom),imom=1,NDATMX)
C
C   BERECHNUNG DER THEORETISCHEN WERTE MIT ANFANGSPARAMETERN,
C   D.H. WENN IFALT2 = 0 , ODER
C   HA-FO-RECHNUNGEN FUER EINZELNE KERNE AUSGEBEN BEI LETZTER
C   RECHNUNG, D.H. WENN IFLAST = 0   .
C
      IF(IFLAST.EQ.1 .OR. IFALT2.EQ.0 .OR. IFREPT.EQ.1) THEN
       IFALT2=1
       IF(IFLAST.EQ.1) IWR=-1
       IWRTE=ist1+1
       IFSAVE=1
       WRITE(IWU,'(A)') ' INITIAL OR FINAL CHI**2'
       CALL CHIFCT(A,Y,SIGMAY,YFIT,NPARMX,CHISU,SUMCHI)
       IF(IERROR.NE.0) THEN
        WRITE(IWU,'(A,I5,A)') ' THE ERROR ',IERROR,
     *     ' OCCURED AT THE INITIAL OR FINAL EVALUATION! STOP.'
        STOP
       ENDIF
       IWRTE=ist1
       IWR=-2
       IF(IFLAST.NE.1) THEN
        CHIIN=CHISU
        CHIINI=CHIIN
        WRITE(IWU,'(//A,E12.4,A//)') '  CHI**2 INPUT = ',CHIIN,
     *                             '  ***********'
        if(KTESTC.eq.KDATMX) then
         DO 10 KL=1,NDATMX
         CHISU0(KL)=SUMCHI(KL)
   10    CHIS00(KL)=CHISU0(KL)
         WRITE(IWU,'(15(1X,F7.2))') (CHISU0(KL),KL=1,NDATMX)
        endif
       ENDIF
      ENDIF
      IF(NWDR.LE.0) RETURN
C
C   BERECHNUNG MESSWERTE MIT UM DPARMO GEAENDERETEN PARAMETERN
C   UND DAMIT ABLEITUNG NACH DEN PARAMETERN :
C
C       VORBESETZUNGEN
C
      DO 110 I=1,KPARMX
      B(I)=A(I)
      ARES(I)=A(I)
  110 PARMOM(I)=A(I)
      IF(IST1.GT.1) WRITE(IWU,'(A)')
     *       '  VARIATION DER PARAMETER'
      CHIRES=CHIIN
      IRESER=0
C
C        LOOP UEBER DIE PARAMETER, JEDER WIRD DER REIHE NACH
C        EINMAL GEAENDERT.
C
      DO 100 JPAR=1,NPARMX
      PARMOM(JPAR)=A(JPAR)+DPARMO(JPAR)
      IFSAVE=0
      if(ist1.ge.1)
     &  WRITE(IWU,'(A,I4)') ' CHI**2 FOR VARIED PARAMETER NR.',JPAR
      CALL CHIFCT(PARMOM,Y,SIGMAY,YFITMO,NPARMX,CHISU,SUMCHI)
       IF(IERROR.NE.0) THEN
        WRITE(IWU,'(A,I5,A)') ' THE ERROR ',IERROR,
     *     ' OCCURED AT THE PARAMETER VARIATION! STOP.'
        STOP
       ENDIF
C
C     RESERVE PARAMETER ABSPEICHERN , WENN SCHON BEI DER
C     PARAMETERVARIATION EIN KLEINERES CHI**2 HERAUSGEKOMMEN IST.
C
C
      CHIMOM(JPAR)=CHISU
      IF((CHIMOM(JPAR)/CHIRES).LE.0.999) THEN
       DO 120 L=1,NPARMX
  120  ARES(L)=PARMOM(L)
       CHIRES=CHIMOM(JPAR)
       IRESER=1
      ENDIF
C
C     ABLEITUNG DER MESSWERTE NACH DEN PARAMETERN. DAS RESULTAT WIRD
C     AUF DEM FELD 'FDERIJ' GESPEICHERT.
C
      if(ist1.eq.1) then
       write(6,'(a,i2,a/a,9g12.5)')
     &   ' test FDERIJ for JPAR=',JPAR,' :',
     &   '   YFIT=',(YFIT(imom),imom=1,NDATMX)
       write(6,'(a,9g12.5)')
     &   ' YFITMO=',(YFITMO(imom),imom=1,NDATMX)
      endif
       DO 130 IDAT=1,NDATMX
       FDERIJ(IDAT,JPAR)=(YFIT(IDAT)-YFITMO(IDAT))/
     *                    (A(JPAR)-PARMOM(JPAR))
  130  CONTINUE
C
C      ZWISCHENAUSDRUCKE
C
      IF(IST1.GT.1) THEN
       WRITE(IWU,'(I3,8E14.6,E13.3/3X,7E14.6,E12.4,I3)')
     *      JPAR,(PARMOM(KL),KL=1,8),DPARMO(JPAR),
     *      (PARMOM(KL),KL=9,KPARMX),CHIMOM(JPAR),ICHI
       if(KTESTC.eq.KDATMX) then
        WRITE(IWU,'(15(1X,F7.2))') (SUMCHI(KL),KL=1,NDATMX)
        DO 140 KL=1,NDATMX
  140   SUMCHI(KL)=(SUMCHI(KL)-CHISU0(KL))
        WRITE(IWU,'(15(1X,F7.2))') (SUMCHI(KL),KL=1,NDATMX)
        DO 150 KL=1,NDATMX
  150   SUMCHI(KL)=SUMCHI(KL)/DPARMO(JPAR)
        WRITE(IWU,'(15(1X,F7.2))') (SUMCHI(KL),KL=1,NDATMX)
       endif
      ENDIF
C
      PARMOM(JPAR)=A(JPAR)
  100 CONTINUE
C
      IF(IRESER.NE.0) WRITE(IWU,'(A,2(/1X,10G13.6))')
     *             '   RESERVE PARAMETER:',(ARES(J),J=1,NPARMX)
C
C     SECOND DERIVATIVE OF CHI**2 ON THE FIELD 'ALPHA' ,
C     THE GRADIENT OF CHI**2 ON THE FIELD 'BETA' .
C
      DO 210 J=1,NPARMX
      BETA(J)=0.
      DO 210 K=1,J
      ALPHA(K,J)=0.
  210 ALPHA(J,K)=0.
C
      DO 200 J=1,NPARMX
C
      DO 220 IDAT=1,NDATMX
      BETA(J)=BETA(J)+(Y(IDAT)-YFIT(IDAT))*
     +           FDERIJ(IDAT,J)/SIGMAY(IDAT)**2
  220 CONTINUE
C
      DO 230 K=1,J
      DO 240 IDAT=1,NDATMX
      ALPHA(J,K)=ALPHA(J,K)
     &           +FDERIJ(IDAT,J)*FDERIJ(IDAT,K)/SIGMAY(IDAT)**2
  240 CONTINUE
      ALPHA(K,J)=ALPHA(J,K)
  230 CONTINUE
  200 CONTINUE
C
C  NICHT-KORRELIERTE FEHLER
C
      DO 250 J=1,NPARMX
      IF(ALPHA(J,J).GT.0.) THEN
       SIGMA0(J)=1./SQRT(ALPHA(J,J))
      ELSE
       SIGMA0(J)=0.
      ENDIF
  250 CONTINUE
C
C    CHI**2 MIT ANFANGSWERTEN:
C
      CHISQ1=CHIIN
C
C        INVERT MODIFIED CURVATURE MATRIX TO FIND NEW PARAMETERS
C
      DO 300 IWARN2=1,6
      IF(IWARN2.EQ.1 .OR. IREDUC.EQ.0) THEN
       FLAMDA = .0001*10.**IWARN2
       IF(IWARN2.GT.1)
     *   WRITE(IWU,'(//A,I1,A/A,E10.2)')
     *             '  LAMDA ERHOEHT UM 10**',IWARN,
     *     '   OHNE VERBESSERUNG IN CHI**2','  LAMBDA JETZT ',FLAMDA
       DO 310 J=1,NPARMX
       DO 320 K=1,NPARMX
  320  ZARRAY(J,K)=ALPHA(J,K)/SQRT(MAX(XLOWLM,ALPHA(J,J)*ALPHA(K,K)))
  310  ZARRAY(J,J)=1.+FLAMDA
       CALL MATIN2(ZARRAY,NPARMX,KPARMX,NERR,DET)
       if(ist1.ge.1) write(6,'(a,i5,g12.5)')
     &   ' invert CHI**2 matrix with NERROR,DET=',NERR,DET
       IF(IFLAST.EQ.1) GOTO 399
       DO 330 J=1,NPARMX
       B(J)=0.0
       DO 330 K=1,NPARMX
  330  B(J)=B(J)+BETA(K)*ZARRAY(J,K)
     *         /SQRT(MAX(XLOWLM,ALPHA(J,J)*ALPHA(K,K)))
       bmax   = 0.0
       do 331 J=1,NPARMX
  331  bmax  = max(abs(b(J))/dparmo(J),bmax)
       bfac  = min(one,10.0/max(xlowlm,bmax))
       do 332 J=1,NPARMX
  332  B(J) = A(J)+bfac*B(J)
       WRITE(IWU,'(15F8.1)') (B(K),K=1,NPARMX)
      ELSE
       DO 315 J=1,NPARMX
  315  B(J)=A(J)+(B(J)-A(J))*.5
       WRITE(IWU,'(A)') ' STEP SHORTENED BY FACTOR 2. TRY AGAIN.'
      ENDIF
C
C     BERECHNE ABSCHAETZUNG DES NEUEN CHI**2 AUS SCHRITT UND
C     KRUEMMUNGSMATRIX.
C
      SUM1=0.
      SUM2=0.
      DO 360 J=1,NPARMX
      SUM1=(B(J)-A(J))*BETA(J)+SUM1
      DO 360 K=1,NPARMX
  360 SUM2=(B(J)-A(J))*ALPHA(J,K)*(B(K)-A(K))+SUM2
      CHITST=CHISQR-SUM1-SUM1+SUM2
      WRITE(IWU,'(4(A,G12.5))') ' ESTIMATE NEW CHI**2:',CHITST,
     *       '      FOR CHISQR=',CHISQR,'   LINEAR CONTR.=',SUM1,
     *       '   QUADR. CONTR.=',SUM2
C
C    BERECHNE CHI**2 MIT DEN NEUEN PARAMETERN (FELD B)
C
      IFSAVE=1
      IF(IST1.GE.1) IWRTE=ist1+1
      WRITE(IWU,'(A)') ' EVALUATE NEW CHI**2'
      CALL CHIFCT(B,Y,SIGMAY,YFITMO,NPARMX,CHISU,SUMCHI)
      IF(IERROR.NE.0) THEN
        WRITE(IWU,'(A,I5,A)') ' THE ERROR ',IERROR,
     *     ' OCCURED AT THE NEW CHI**2 EVALUATION! STOP.'
        CHISQR=1./XLOWLM
      ELSE
       CHISQR=CHISU
       if(KTESTC.eq.KDATMX) then
        DO 340 KL=1,NDATMX
  340   CHIS00(KL)=SUMCHI(KL)
       endif
      ENDIF
       IWRTE=0
       WRITE(IWU,'(/A,G12.4,A,G12.4,A,I3/)')
     *          '  CHI**2:   ALT: ',CHISQ1,'     NEU: ',CHISQR,
     *         '     IWARN2 = ',IWARN2
       IF(CHISQR.LE.CHISQ1.AND.CHISQR.LE.CHIRES) GO TO 399
       IF(IRESER.NE.0) GOTO 309
  300 CONTINUE
C
C   WAR CHI*2 BEI VARIATION EINES PARAMETERS KLEINER GEWORDEN?
C   IN DIESEM FALL IST IRESER=1 GESETZT
C
  309 CONTINUE
      IF(IRESER.NE.0) THEN
       DO 350 J=1,NPARMX
       B(J)=ARES(J)
       A(J)=ARES(J)
  350  CONTINUE
       WRITE(IWU,'(/A/A/A/A//)')
     *          '  CHI**2 BEI VARIATION EINES PARAMETERS KLEINER',
     * '  DAHER ENTSPRECHENDE PARAMETER GEWAEHLT FUER FORTSETZUNG',
     * '  DA DIE BERECHNETE VARIATION NICHT ZU VERBESSERUNG IN CHI**2',
     * '  FUEHRTE'
       IRESER=0
       IWARN  = 1
       CHISQR=CHIRES
      ELSE
       IWARN = 99
      ENDIF
C
C        EVALUATE PARAMETERS AND UNCERTAINTIES
C
  399 CONTINUE
C
C   HAT SICH CHI**2 VERBESSERT?
C
      VERGLC=CHIIN-CHISQR
      IF(VERGLC.LT.0.0) THEN
       CHISQR=CHIIN
       if(KTESTC.eq.KDATMX) then
        DO 410 KL=1,NDATMX
  410   CHIS00(KL)=CHISU0(KL)
       endif
       DO 420 J=1,NPARMX
  420  B(J)=A(J)
       ITERA = ITERA - 1
       IWARN = 99
      ELSE
       IWARN = 1
      ENDIF
C
C    AUSDRUCK DER RESULTATE DES BEVINGTON STEPS
C
      IF(IFLAST.NE.1) THEN
       WRITE(IWU,'(/A)') '  NEUE PARAMETERWERTE:'
       IF(IWARN.GE.99) WRITE(IWU,'(A/)')
     *        '  = ALTE PARAMETER (KEINE VERBESSERUNG IN CHI**2'
       WRITE(IWU,'(5(2X,A,I2,A,G14.7))')
     *        (PARBRA,J,PARKET,B(J),J=1,NPARMX)
       CHIWRI=CHISQR/RFREE
       WRITE(IWU,'(2(A,E12.4))') ' CHI**2 PRO FREIHEITSGRAD = ',CHIWRI,
     *                ',      CHI**2 GESAMT            = ',CHISQR
      ENDIF
C
      IF(IFLAST.NE.1 .AND. IWARN.LT.99 .AND. NWDR.NE.1) THEN
       DO 440 J=1,NPARMX
       A(J)=B(J)
  440  CONTINUE
       DO 450 IDAT=1,NDATMX
  450  YFIT(IDAT)=YFITMO(IDAT)
      ENDIF
C
      if(KTESTC.eq.KDATMX) then
       DO 460 KL=1,NDATMX
  460  CHISU0(KL)=CHIS00(KL)
      endif
C
C     ERROR RETURN IF NO IMRPOVEMENT HAS BEEN ACHIEVED
C
      IF(IWARN.GE.99) THEN
       WRITE(IWU,'(A)') ' NO IMPROVMENT IN CHI**2. RETURN EARLY.'
       RETURN
      ENDIF
C
C  FEHLERMATRIX AUF 'ZARRAY'
C
      DO 470 J=1,NPARMX
      DO 470 K=1,NPARMX
  470 ZARRAY(J,K)=ALPHA(J,K)/SQRT(MAX(XLOWLM,ALPHA(J,J)*ALPHA(K,K)))
      CALL MATIN2(ZARRAY,NPARMX,KPARMX,NERR,DET)
      DO 480 J=1,NPARMX
  480 SIGMAA(J)
     &   =SQRT(max(xlowlm**1.5,ZARRAY(J,J))/MAX(XLOWLM,ALPHA(J,J)))
      DO 490 J=1,NPARMX
      DO 490 K=1,NPARMX
  490 ZARRAY(J,K)=ZARRAY(J,K)/SQRT(MAX(XLOWLM,ALPHA(J,J)*ALPHA(K,K)))
C
C     PRINT CHI**2 MATRIX
C
      IF(IST1.EQ.1) THEN
       WRITE(IWU,'(A)')   ' CHI**2 MATRIX:'
       do 390 K=1,NPARMX
  390  WRITE(IWU,'(10F8.3)')
     &  (ALPHA(J,K)/SQRT(MAX(XLOWLM,ALPHA(J,J)*ALPHA(K,K))),J=1,NPARMX)
      ENDIF
C
C     SAVE PARAMETERS AND CORRESPONDING CHI**2 FOR CURVED STEP.
C
      IF(ITERA.LE.MAX0(3,IFCURV)) THEN
       DO 510 I=1,KPARMX
       PARSAV(I,1) = PARSAV(I,2)
       PARSAV(I,2) = PARSAV(I,3)
  510  PARSAV(I,3) = A(I)
       CHISAV(1) = CHISAV(2)
       CHISAV(2) = CHISAV(3)
       CHISAV(3) = CHISQR
      ELSE
       DO 520 I=1,KPARMX
  520  PARSAV(I,3) = B(I)
       CHISAV(3)=CHISQR
      ENDIF
      CHIVRG = CHISAV(3)
C
C
C        NOW THE CURVED STEP
C
      IF(ITERA.GE.MAX0(3,IFCURV) .AND. IFLAST.NE.1) THEN
       WRITE(IWU,'(A,3G14.7)') ' CHISAV:',(CHISAV(J),J=1,3)
       WRITE(IWU,'(A)') ' PARSAV:'
       DO 710 J=1,10
  710  WRITE(IWU,'(8X,3G14.7)') (PARSAV(J,I),I=1,3)
       WRITE(IWU,'(A)') ' RUNNING CURVED STEP '
C
CVA1    X3=-1.
CVA1    X1=+1.
C
CVA2        X3=0.
CVA2        DO 720 I=1,KPARMX
CVA2  720   IF(SIGMAA(I).GT.1.E-20)
CVA2     *  X3=X3-(PARSAV(I,3)-PARSAV(I,2))**2/SIGMAA(I)
CVA2        X1=0.
CVA2        DO 730 I=1,KPARMX
CVA2  730   IF(SIGMAA(I).GT.1.E-20)
CVA2     *  X1=X1+(PARSAV(I,1)-PARSAV(I,2))**2/SIGMAA(I)
C
CVA3        X3=0.
CVA3        DO 740 I=1,KPARMX
CVA3  740   IF(SIGMA0(I).GT.1.E-20)
CVA3     *  X3=X3-(PARSAV(I,3)-PARSAV(I,2))**2/SIGMA0(I)
CVA3        X1=0.
CVA3        DO 750 I=1,KPARMX
CVA3  750   IF(SIGMA0(I).GT.1.E-20)
CVA3     *  X1=X1+(PARSAV(I,1)-PARSAV(I,2))**2/SIGMA0(I)
C
C     CURVE IN PARAMETER SPACE IS FITTED BY A PARABOLA, USING
C     THE CURVATURE MATRIX OF CHI**2 AS A METRIC.
C
        X3=0.
        DO 760 I=1,NPARMX
        DO 760 J=1,NPARMX
  760   X3=X3+(PARSAV(I,3)-PARSAV(I,2))*ALPHA(I,J)
     *               *(PARSAV(J,3)-PARSAV(J,2))
        X1=0.
        DO 770 I=1,NPARMX
        DO 770 J=1,NPARMX
  770   X1=X1+(PARSAV(I,1)-PARSAV(I,2))*ALPHA(I,J)
     *               *(PARSAV(J,1)-PARSAV(J,2))
        X1=SQRT(X1)
        X3=-SQRT(X3)
C
C     EQUILIBRISED CURVED STEP:
C      A THIRD POINT ON THE TANGENTIAL CURVE IS GENERATED WITH
C      SAME DISTANCE TO POINT 3 AS THE DISTANCE POINT3 TO POINT 2.
C      THE MINIMUM CHI**2 IS THEN EXTRAPOLATED BETWEEN THIS THREE
C      POINTS.
C
       X0=X3+X3
       DO 780 I=1,KPARMX
  780  PARMOM(I)= PARSAV(I,1)*(X0-X3)*X0/((X1-X3)*X1)
     *           +PARSAV(I,2)*(X0-X3)*(X0-X1)/(X3*X1)
     *           +PARSAV(I,3)*X0*(X0-X1)/(X3*(X3-X1))
       IFSAVE = 0
       WRITE(IWU,'(A)') ' CHI**2 FOR CURVED STEP'
       CALL CHIFCT(PARMOM,Y,SIGMAY,YFITMO,NPARMX,CHISUS,SUMCHI)
       WRITE(IWU,'(5(2X,A,I2,A,G14.7))')
     *       (PARBRA,J,PARKET,PARMOM(J),J=1,NPARMX)
       WRITE(IWU,'(4(A,G12.5))') ' EQUILIBRISED STEP: X1=',X1,
     *         '   X3=',X3,'   X0=',X0,'   CHISUS=',CHISUS
       IF(IERROR.NE.0) THEN
        WRITE(IWU,'(A,I5,A)') ' THE ERROR ',IERROR,
     *     ' OCCURED AT EQUILIBRISED CURVED STEP! RETURN.'
        RETURN
       ENDIF
       IF(CHISUS.LT.CHIVRG) CHIVRG = CHISUS
       CHISAV(1) =CHISAV(2)
       CHISAV(2) = CHISAV(3)
       CHISAV(3) = CHISUS
       DO 790 J=1,NPARMX
       PARSAV(J,1) = PARSAV(J,2)
       PARSAV(J,2) = PARSAV(J,3)
  790  PARSAV(J,3) = PARMOM(J)
       X1     =  1.
       X3     = -1.
C
       D1CHI  = -X3*CHISAV(1)/((X1-X3)*X1)
     *          -(X1+X3)*CHISAV(2)/(X3*X1)
     *          -X1*CHISAV(3)/((X3-X1)*X3)
       D2CHI  = 2.*( CHISAV(3)/((X1-X3)*X1)
     *              +CHISAV(2)/(X3*X1)
     *              +CHISAV(1)/((X3-X1)*X3) )
       IF(D2CHI.GT.0.) THEN
        X0     = -D1CHI/D2CHI
       ELSE
         X0    = -2.
       ENDIF
C      IF(IST1.EQ.1) THEN
        WRITE(IWU,'(5(A,G12.5))') ' CURVED STEP: X1 = ',X1,
     *          '   X3 = ',X3,'   D1CHI = ',D1CHI,'   D2CHI = ',
     *          D2CHI,'   X0 = ',X0
C      ENDIF
       XMIN  = X3
       IF(CHISAV(3).GT.CHISAV(2)) XMIN  = 0.0
       DO 600 ITRY=1,ITRMAX
       DO 610 I=1,KPARMX
  610  PARMOM(I)= PARSAV(I,1)*(X0-X3)*X0/((X1-X3)*X1)
     *           +PARSAV(I,2)*(X0-X3)*(X0-X1)/(X3*X1)
     *           +PARSAV(I,3)*X0*(X0-X1)/(X3*(X3-X1))
       IFSAVE = 1
       WRITE(IWU,'(A)') ' TRIAL CHI**2 AFTER CURVED STEP'
       CALL CHIFCT(PARMOM,Y,SIGMAY,YFITMO,NPARMX,CHISU,SUMCHI)
       IF(IERROR.NE.0) THEN
        WRITE(IWU,'(A,I5,A)') ' THE ERROR ',IERROR,
     *     ' OCCURED AT EXTENDED CURVED STEP! PROCEED.'
       ENDIF
       IF(CHISU.LT.CHIVRG) GOTO 699
       WRITE(IWU,'(5(2X,A,I2,A,G14.7))')
     *        (PARBRA,J,PARKET,PARMOM(J),J=1,NPARMX)
       WRITE(IWU,'(A,I2,A,G12.5)') ' CURVED STEP: ITRY=',ITRY,
     *                           '   CHI**2=',CHISU
       X0=XMIN+(X0-XMIN)*.5
  600  CONTINUE
C
C     CURVED FIT FAILED TO PRODUCE LOWER CHI**2 , RESET SAVED
C     PARAMETERS  FOR NEXT BEVINGTON STEP.
C
       WRITE(IWU,'(A)') '0 NO SUCCESS IN CURVED STEP'
       IF(CHISUS.EQ.CHIVRG) THEN
        WRITE(IWU,'(A)') ' BUT SUCCESS IN TRIAL STEP. KEEP THIS.'
        DO 620 I=1,KPARMX
        PARSAV(I,1)= PARSAV(I,2)
        PARSAV(I,2)= PARSAV(I,3)
  620   A(I)       = PARSAV(I,3)
        CHISAV(1) =  CHISAV(2)
        CHISAV(2) =  CHISAV(3)
        CHISQR    =  CHISAV(3)
       ENDIF
       RETURN
C
C     CURVED FIT SUCCEEDED
C
  699  CONTINUE
       WRITE(IWU,'(A,I3)') '0 SUCCESS IN CURVED STEP AT ITRY=',ITRY
       DO 630 I=1,KPARMX
       PARSAV(I,2)=PARMOM(I)
  630  A(I)       =PARMOM(I)
       DO 650 IDAT=1,NDATMX
  650  YFIT(IDAT)=YFITMO(IDAT)
       CHISAV(2) = CHISU
       CHISQR    = CHISU
       WRITE(IWU,'(//A/)') '  NEUE PARAMETERWERTE'
       WRITE(IWU,'(5(2X,A,I2,A,G14.7))')
     *        (PARBRA,J,PARKET,A(J),J=1,NPARMX)
       CHIWRI=CHISU/RFREE
       WRITE(IWU,'(A,E12.4)') ' CHI**2 PRO FREIHEITSGRAD = ',CHIWRI
       WRITE(IWU,'(A,E12.4)') ' CHI**2 GESAMT            = ',CHISU
      ENDIF
C
      RETURN
      END
C-----CHIFCT----------------------------------------------------------
C
      SUBROUTINE CHIFCT(PAR,FEXP,DFEXP,FBER,NPARMX,CHISU,SUMCHI)
C
      INCLUDE 'PARAMS.INC'
      parameter (KPARMX=12, KDATMX=12, KTESTC=1)
C
      DIMENSION FEXP(KDATMX),DFEXP(KDATMX),FBER(KDATMX),
     *          PAR(KPARMX)
      DIMENSION SUMCHI(KTESTC)
c
      INCLUDE 'CONSTN.INC'
c
      common /mompar/nbel,nmomin,nmomax,nrmom,iwrte
      data xlowmn/1.0e-30/
      data fifty/50.0/
      data xupmin/1.0e37/
c
c------------------------------------------------------------------------------
c
      if(NPARMX.ne.2*nbel) stop ' wrong NPARMX in CHIFCT'
      do 110 imom=1,nrmom
      mom    = nmomin+2*(imom-1)
      sum    = 0.0
      do 120 iq=1,nbel
      if(par(iq).lt.-fifty) then
       expmom = 0.0
      elseif(par(iq).gt.fifty) then
       expmom = xupmin
      else
       expmom = exp(par(iq))
      endif
      if(mom*par(iq+nbel).lt.-fifty) then
       expenk = 0.0
      elseif(mom*par(iq+nbel).gt.fifty) then
       expenk = xupmin
      else
       expenk = exp(par(iq+nbel)*mom)
      endif
      if(expmom.ge.xupmin .or. expenk.ge.xupmin) then
       sum = xupmin
       goto 129
      else
       sum    = expmom*expenk + sum
      endif
  120 continue
  129 continue
      sum    = min(sum,xupmin)
      sum    = max(sum,xlowmn)
c      if(sum.gt.xupmin .or. sum.le.zero) then
c       fber(imom) = xupmin
c      else
       fber(imom) = log(sum)
c      endif
  110 continue
      if(IWRTE.ge.2) then
        write(6,'(a/a,6g12.5,10(/5x,6g12.5))')
     &  ' in CHIFCT:',' par=',(par(iq),iq=1,2*nbel)
        write(6,'(a,6g12.5,10(/5x,6g12.5))')
     &  ' FBER:',(fber(imom),imom=1,nrmom)
      endif
c
      sum    = 0.0
      do 10 imom=1,nrmom
      addchi = ((fexp(imom)-fber(imom))/dfexp(imom))**2
      if(KTESTC.eq.KDATMX) sumchi(imom) = addchi
   10 sum    = addchi+sum
      chisu  = sum
c
      return
      end
C-----SPRESO-----------------------------------------------------
C
      SUBROUTINE SPRESO
C
      INCLUDE 'PARAMS.INC'

C
C     SPRESO COMPUTES THE MULTIPOLE STRENGTH FUNCTION OF THE
C     MEAN-FIELD HAMILTONIAN IN THE CONTINUUM.
C
      INCLUDE 'SKYPAR.INC'
      INCLUDE 'WAVFUN.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'RESULT.INC'
      INCLUDE 'NUMPAR.INC'
      INCLUDE 'CONSTN.INC'
      INCLUDE 'LINRES.INC'
C
      DIMENSION NPARSV(KRPMAX),NHOLSV(KRPMAX),ISOPHS(KRPMAX)
C
      CHARACTER*16 TYPMOD
      CHARACTER*12 JSELEC
      CHARACTER*50 TOPER
      character*2  ncol
C
      DATA IFOP14,IFOP17/0,0/
      data ncol/'N:'/
cc      data ncol/'  '/
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C
C      READ PARAMETERS FOR 'SUMRUL'
C
       READ(05,*)   LRPA
       READ(05,*)   ESTART,ESTEP,EEND,ENWID
       NUMENR = (EEND-ESTART)/ESTEP
C
C      INITIALISE SOME FIELDS WITH ZEROS
C
       WRITE(6,'(A)') ' SPRESO ENTERED '
C
C     OPEN FILE 14 FOR CONCENTRATED OUTPUT OF MATRIX ELEMENTS
C
      IF(IFOP14.EQ.0) THEN
       OPEN(UNIT=14,STATUS='UNKNOWN',FORM='FORMATTED',FILE='FOR014')
       IFOP14 = 1
      ENDIF
C
C     WRITE A HEADER
C
      if(ncol.eq.'N:') then
       write(14,'(4(a/))')
     &   'NEXT:',
     &   'S.P.Reson.',
     &   'X:E',
     &   'Y:Im(Resp)'
      else
       write(14,'(a)') '1'
      endif
      WRITE(14,'(a,2A/a,2(A,I4),5X,2(A,I4))')
     &  ncol,'  TITLE:',TITLE,
     &  ncol,'NUCLEUS:  NPROT=',NPROT,'  NNEUT=',NNEUT,
     &  ' SHELLS:  NMAXP=',NMAXP,'  NAMXN=',NMAXN
      WRITE(14,'(a,A,5G13.6/a,5X,A,4G13.6/a,5X,A,G13.6)')
     &  ncol,'FORCE: T0-T4=',T0,T1,T2,T3,T4,
     &  ncol,      ' X0-X3=',X0,X1,X2,X3,
     &  ncol,      ' POWER=',POWER
      Write(14,'(a,5x,3(a,g13.6)/)')
     &  ncol,      ' TYUK=',tyuk,' XYUK=',xyuk,' AYUK=',ayuk
      WRITE(14,'(a,2(A,F10.3)/a,16X,2(A,F10.3))')
     &  ncol,'SOME PROPERTIES: RMS-RADIUS=',RMS,
     &       '    DIFFR.RADIUS=',RDMS,
     &  ncol,'   SURFACE=',SURF,'          ENERGY=',ETOTAL
c
c     provide derivative of "rho" and "bmass"
c
      CALL DERV17(BMASS(1),D1B(1),ONE,NGRID,RSTEP)
      CALL DERV17(BMASS(NGRID+1),D1B(NGRID+1),ONE,NGRID,RSTEP)
      CALL DERV17(RHO(1),D1R(1), ONE,NGRID,RSTEP)
      CALL DERV17(RHO(ngrid+1),D1R(ngrid+1), ONE,NGRID,RSTEP)
C
C     WRITE OPTIONS AND PARAMETERS FOR RPA PART
C
      WRITE(6,'(/1X,A,I3,3(A,g10.3))')
     &   ' S.P. Resonances with: LRPA=',LRPA
      WRITE(14,'(/a,A,I3,3(A,g10.3))')
     &   ncol,' S.P. Resonances with: LRPA=',LRPA
      Write(14,'(a,4(a,f8.3))')
     &   ncol, '  ESTART=',ESTART,'  ESTEP=',ESTEP,'  EEND=',EEND,
     &    '  ENWID=',ENWID
      Write(6,'(6(a,f8.3))')
     &    '  ESTART=',ESTART,'  ESTEP=',ESTEP,'  EEND=',EEND,
     &    '  ENWID=',ENWID
C
       IFCMPL = 1
       ALRPA  = LRPA
       MLRPA  = 2*LRPA+1
C
       CALL BUFFQL(QLROUT,'C','      ')
C
C      read type of external operator in RPA-standard form,
C      and set-up he operator in QINIT
C
      READ(5,*) TYPMOD,ISORPA,JSELEC,ISPLIT,IHWGT,IH0WGT,IPWGT,IPRSHL
      write(6,'(a,i5,a,5i5)')
     &   TYPMOD,ISORPA,JSELEC,ISPLIT,IHWGT,IH0WGT,IPWGT,IPRSHL
      CALL QINIT(TYPMOD,JSELEC,ISPLIT,IHWGT,IH0WGT,IPWGT,
     &           IPRSHL,IMODE,IFPRIM,IFLAGR,XKQRD,RESFAC,CMFAC,
     &           NPARSV,NHOLSV,ISOPHS,TOPER)
      call buffql(qlrout,'O',' QINIT')
      ovqq = qpover(qlrout,qlrout)
C
      write(14,'(/2a)')
     &   ncol,' ENERGY  Re(Resp)    Im(Resp)'
      if(ncol.eq.'N:')
     &   write(14,'(a)')  'H:  X      N       Y1,L'
      DO 888 IENR=0,NUMENR
      ENEXT  = ESTART+ESTEP*IENR
      call copywf(qlrout,plrout)
      call DAMPH1(plrout,savlr,enext,enwid,1)
      ovqp = qpover(savlr,qlrout)*enwid/ovqq
      ovqpr = qpover(plrout,qlrout)*enwid/ovqq
      write(6,'(3(a,g12.4))')  ' Energy=',enext,
     &    ' Re(R)=',ovqpr,' Im(R)=',ovqp
      write(14,'(1x,f7.3,2g12.4)') enext,ovqpr,ovqp
  888 CONTINUE
C
      if(ncol.eq.'N:')
     &  write(14,'(/2(a,f6.2),a,a)')
     &  'S: DIST(4.0) THICKN(0.7) XMIN(',ESTART,') XMAX(',EEND,')',
     &  ' YMIN(0.0)'
C
      RETURN
      END
C-----DAMPH1------------------------------------------------------
C
      SUBROUTINE DAMPH1(QWF,QIWF,EREAL,EIMAG,NRINV)
C
      INCLUDE 'PARAMS.INC'
C
C     COMPUTES DAMPING BY  (H0-E)**-1
C     WITH COMPLEX ENERGY 'E' = 'EREAL'-i*'EIMAG'
C
      DIMENSION QWF(KX11,KLRPA)
      DIMENSION QIWF(KX11,KLRPA)
C
      INCLUDE 'WAVFUN.INC'
      INCLUDE 'POTINV.INC'
      INCLUDE 'SKYPAR.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'RESULT.INC'
      INCLUDE 'LINRES.INC'
      INCLUDE 'NUMPAR.INC'
      INCLUDE 'CONSTN.INC'
C
C-------------------------------------------------------------------------------
C
      iftsin = 0
C
C     PROVIDE NECESSARY QUANTUM NUMBERS AND POINTERS WITH
C     DO LOOPS 10, 20 AND 30
C
      INCLUDE 'INITDO.INC'
C
C     ACTION OF LINEAR RESPONSE HAMILTONIAN, result on 'hlinrs'
C
C     PREPARE POTENATIAL- AND MASS-FIELDS FOR INVERSION
C
      DO 110 I=1,NGRID
  110 SQB(I) = SQRT(BMASS(I+ISOADD))
C
      CALL DERV25(SQB(1),D2B(1),ONE,NGRID,RSTEP)
      DO 120 I=1,NGRID
  120 D2B(I) = D2B(I)/SQB(I)
C
C      THE EFFECTIVE POTENTIAL FOR GIVEN 'JSA' AND 'ISO' ON 'POTL';
C
C      'PHASJS' IS THE DEGENERACY OF THE SHELL.
C      'JSNPL' IS POINTER TO TO NUCLEON SHELL IN POINTER FIELD 'NPLA'.
C
      XLP    = FLOAT(LL*LL+LL)
      XL0    = WLS(JS)
C
      D1B0=(16.0*BMASS(2+ISOADD)-15.0*BMASS(1+ISOADD)
     *        -BMASS(3+ISOADD))/(6.0*RSTEP**RSTEP)
      POTL(1) = POT(1+ISOADD)+XL0*POTLS(1+ISOADD)+D1B0
      DO 160 I=2,NGRID
      N      = I+ISOADD
      POTL(I) = POT(N)+XL0*POTLS(N)+XLP*BMASS(N)*ZENTRF(I)
     &          +D1B(I+isoadd)*rinvs(i)
  160 CONTINUE
C
C     INVERT WITH MEAN-FIELD HAMILTONIAN MINUS EXTERNAL ENERGY
C     AND PROPAGATE TO IMPROVED GUESS FOR 'PLROUT'
C     'QWF' IS USED AS WORKSPACE IN THE INVERSION STEPS
C
c      CALL PROJEC(QWF(NADD+1,JDIF),JS,NSTACK,NMAX)
      DO 140 IR=1,NGRID
      WORK1(IR)=QWF(NADD+IR,JDIF)*((IR-1)*RSTEP)
      WORK3(IR)=WORK1(IR)
 140  CONTINUE
      CALL INVCOM(WORK1,WORK2,WORK1,SUBFAC+EREAL,EIMAG,LL)
       DO 151 IR=2,NGRID
       QWF(IR+NADD,JDIF) = WORK1(IR)*RINVS(IR)
  151  QIWF(IR+NADD,JDIF) = WORK2(IR)*RINVS(IR)
       if(iftsin.eq.1) then
         write(6,'(a,i5,a,4i5)')
     &    ' test INVCOM for state nr.',NA,
     &    ' JS,LL,JSA,LA=',js,ll,jsa,la
        CALL HPSI(WORK1(1),AUX1(1),AUX2(1),AUX3(1),1,JS,ISO,
     &            SUBFAC+EREAL)
        anorm  = zero
        onorm  = zero
        oaovlp = zero
        DO 710 IR=1,NGRID
        anorm  = aux1(ir)*aux1(ir) + anorm
        onorm  = work3(ir)*work3(ir) + onorm
        oaovlp = aux1(ir)*work3(ir) + oaovlp
  710   CONTINUE
        cosang = oaovlp*oaovlp/max(anorm*onorm,xsmall)
        write(6,'(a,3g11.4,A,G11.4)')
     &     ' norms: a,o,oa=',anorm,onorm,oaovlp,' cosang=',cosang
        write(6,'(a,99(/1x,f8.3,3g12.4))')
     &   '   R       before       after        error',
     &   ((ir-1)*rstep,work3(ir),aux1(ir),
     &     (work3(ir)-aux1(ir))/work3(ir),ir=2,ngrid)
        CALL HPSI(WORK2(1),AUX1(1),AUX2(1),AUX3(1),1,JS,ISO,
     &            SUBFAC+EREAL)
        anorm  = zero
        DO 720 IR=1,NGRID
c       AUX1(IR) = AUX1(IR)-EIMAG*WORK1(IR)
        anorm  = (AUX1(IR)-EIMAG*WORK1(IR))**2 + anorm
  720   CONTINUE
        write(6,'(a,g12.4)')
     &     ' norm imaginary part:',anorm
       endif
c
      PARLL  = (1-MOD(LL,2)*2)
      CALL EXTRAP(QWF(NADD+1,JDIF),PARLL)
      CALL PROJEC(QWF(NADD+1,JDIF),JS,NSTACK,NMAX)
      CALL EXTRAP(QIWF(NADD+1,JDIF),PARLL)
      CALL PROJEC(QIWF(NADD+1,JDIF),JS,NSTACK,NMAX)
C
      INCLUDE 'ENDDO.INC'
C
      RETURN
      END
C-----INVCOM------------ PART OF THE QHQ PACKAGE -------------------------------
C
      SUBROUTINE INVCOM(RREAL,RIMAG,RHS,EINVR,EINVI,Langmo)
C
      INCLUDE 'PARAMS.INC'
C
      DIMENSION RREAL(KG1),RIMAG(KG1),RHS(KG1)
C
      INCLUDE 'WAVFUN.INC'
      INCLUDE 'POTINV.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'CONSTN.INC'
C
       DATA FSX,EZW/.83333333333333,.8333333333333E-1/
cc       DATA FSX,EZW/1.0,0.0/
C
C
C      DETERMINES 'R' AS SOLUTION OF THE INHOMOGENOUS LINEAR DIFF. EQ.:
C       " (-D/DX (SQB(X)**2) D/DX + (POTL(X) - EINV ) ) R = RHS  "
C      FOR GIVEN 'RHS' WHERE 'SQB**2' IS X-DEPENDENT INVERSE MASS,
C      'POTL' A POTENTIAL AND 'EINV' A SHIFT CONSTANT.
C      THE ENERGY-SHIFT 'EINV' IS COMPLEX WITH REAL PART 'EINVR' AND
C      IMAGINARY PART 'EINVI'.
C      ACCORDINGLY THE RESPONSE FIELD 'R' IS COMPLEX WITH REAL PART 'RREAL'
C      AND IMAGINARY PART 'RIMAG'.
C
C      THE ARRAYS IN THE CALLING LIST MEAN:
C       RREAL  = REAL PART OF THE THE SOLUTION OF THE EQUATION
C       RIMAG  = IMAGINARY PART OF THE THE SOLUTION OF THE EQUATION
C       RHS    = THE RIGHT HAND SIDE (INHOMOGENEITY) OF THE EQUATION.
C                R AND RHS MAY BE THE SAME ARRAYS IN THE CALLING ROUTINE
C                RHS IS THEN OVERWRITTEN.
C      THE PARAMETERS IN THE CALLING LIST MEAN:
C       EINVR  = REAL PART OF THE ENERGY-SHIFT
C       EINVI  = IMAGINARY PART OF THE ENERGY-SHIFT
C
C      INFORMATION ON THE HAMILTONIAN IS TRANSFERRED VIA COMMON /POTINV/:
C       SQB    = INVERSE SQUARE-ROOT OF X-DEPENDENT EFFECTIVE MASS
C       POTL   = EFFECTIVE POTENTIAL BEFORE MASS TRANSFORMATION,
C                TO BE PREPARED IN SUBROUTINE 'WELFU'.
C
C      INFORMATION ON THE GRID IS TRANSFERRED VIA COMMON /GRID/.
C
C      NOTE: THE PROGRAM IS WRITTEN FOR SOLVING THE SCHROEDINGER EQUATIO
C            IN RADIALLY SYMMETRIC PROBLEMS. THUS IT IMPLIES THE RHS AND
C            POTEFF TO BE SUCH THAT R(X) VANISHES AT THE ORIGIN.
C
C      THE SOLUTION PROCEEDS AS FOLLOWS:
C
C      BY INTRODUCING SCALED WAVE FUNCTIONS
C       " REFF = R/SQB "
C      AND
C       " RHSEFF = RHS/SQB ".
C      THE EQUATIONS IS TRANSFORMED TO A LAPLACIAN EQUATION
C       " (-D/DX**2 + POTEFF(X) ) REFF = RHSEFF "
C      WITH CONSTANT MASS AND AN EFFECTIVE POTENTIAL
C       " POTEFF = (POTL(X)-UMIN)/(SQB(X)*SQB(X)) + D2B(X) "
C      WHERE
C       " D2B = (D/DX**2 SQB)/SQB " .
C      AFTER SOLVING THE LAPLACIAN EQUATION, THE SOLUTION IS TRANSFORMED
C      BACK BY
C       " R   = REFF*SQB" .
C
C      THE D/DX**2 IS TREATED EFFECTIVELY IN 5-POINT PRECISION.
C      IN ADDITION TO A STANDARD 3-POINT FORMULA FOR D/DX**2
C      A COMPLEMENTING 3-POINT AVERAGING OF UEFF*REFF AND RHSEFF
C      IS USED WHICH BOTH TOGETHER MAKE UP THE 5-POINT PRECISION.
C      (SEE R.ZURMUEHL, "PRAKTISCHE MATHEMATIK", CH. VII, SECT. 3,
C       SPRINGER VERLAG, BERLIN 1963)
C
C
C      THE WORKSPACES FOR THE HAMILTONIAN, 'POTEFF', 'AINV' AND 'BINV' ARE
C      NOW BECOMING COMPLEX. WE AVOID EXPLICIT COMPLEX ARITHMETIC
C      AND ASSIGN THE WORKSPACES 'AUX..' AS FOLLOWS:
C         AUX1   = IMAGINARY PART OF 'AINV'
C         AUX2   = IMAGINARY PART OF 'BINV'
C         AUX3   = IMAGINARY PART OF 'POTEFF'
C
C-----------------------------------------------------------------------
C
C      SQUARE-ROOT OF INVERSE EFF. MASS ON 'SQB',
C      SECOND DERIVATIVE OF 'SQB' ON 'D2B'.
C      FINALLY ASSEMBLING  THE EFFECTIVE POTENTIAL ON 'POTEFF'.
C
       SQM(1) = ONE/SQB(1)
       DO 80 I=2,NGRID
       SQM(I)    = ONE/SQB(I)
       x         = SQM(I)*SQM(I)
       POTEFF(I) = D2B(I)+(POTL(I)-EINVR)*x
       AUX3(I)   = (-EINVI)*x
   80  CONTINUE
C
C      TRANSFORM TO EFFECTIVE R.H.S., STORE RESULT ON R.
C
       DO 17 I=1,NGRID
   17  RREAL(I)   = RHS(I)*SQM(I)
C
C     PREPARE PARAMETERS FOR KINETIC ENERGY AND LOOP-ENDS.
C
       OFFDIG = ONE/(RSTEP*RSTEP)
       OFFDG2 = OFFDIG+OFFDIG
       NMAXM  = NGRID-1
c
c      the centrifugal force at R -> 0 requires a special treatment
c      of the first few places. This is handled with the variable XLL.
c
       if(langmo.eq.1 .and. fsx.lt.one) then
cc        xll    = (langmo*langmo+langmo)*offdig/144.0
        xll    = (langmo*langmo+langmo)*offdig/12.0
       else
        xll    = zero
       endif
C
C      BEGIN FIRST PASS ELIMINATING THE LOWER OFF-DIAGONAL ELEMENT.
C      NEW DIAGONAL AND UPPER OFF-DIAGONAL ELEMNTS ARE SAVED ON A AND B.
C      THE SECOND PASS, FORWARD ELIMINATION OF THE R.H.S., IS DONE WITHI
C      THE SAME LOOP; RESULT IS STORED ON FIELD T.
C
C
       BINV(1) = ZERO
       AUX1(1) = ZERO
       AINV(1) = ONE
       AUX2(1) = ZERO
       XR      = ZERO
       XI      = ZERO
       YR      = ONE
       YI      = ZERO
       RIM1    = RREAL(1)
C
ctest
ctest       write(6,'(a)')'  r        ar          ai          br          bi'
ctest
       DO 13 I=2,NMAXM
C                                         FORWARD ELIM. FOR R.H.S.
       BDYI    = YR*YR+YI*YI
       BDYR    = (BINV(I-1)*YR+AUX2(I-1)*YI)/BDYI
       BDYI    = (AUX2(I-1)*YR-BINV(I-1)*YI)/BDYI
       X       = EZW*(RREAL(I+1)+RIM1)+FSX*RREAL(I)
     &          - (XR*BDYR-XI*BDYI)
       XI      = -XR*BDYI-XI*BDYR
       XR      = X
       RIM1    = RREAL(I)
       RREAL(I)= XR
       RIMAG(I)= XI
C                                         FORWARD ELIM. FOR MATRIX
       UR      = POTEFF(I)
       UI      = AUX3(I)
       ZR      = -OFFDIG+EZW*UR
       ZI      = EZW*UI
c       if(i.eq.3)  then
c        zr = -offdig
c        zi = zero
c       endif
cc       if(i.eq.3)  zr = zr-xll
       BINV(I) = ZR
       AUX2(I) = ZI
       YR      = OFFDG2+FSX*UR-(ZR*BDYR-ZI*BDYI)
       YI      = FSX*UI-(ZR*BDYI+ZI*BDYR)
cc       if(i.eq.3)  zr = zr+xll
cc       if(i.eq.2) yr = yr+xll*16.0
       if(i.eq.2) yr = yr+xll
c       if(i.eq.2) then
c        yr = yr+(ezw+ezw)*ur
c        yi = yi+(ezw+ezw)*ui
c       endif
       AINV(I) = YR
       AUX1(I) = YI
ctest
ctest       write(6,'(1x,f4.2,4g12.4)') (I-1)*rstep,yr,yi,zr,zi
ctest
   13  CONTINUE
C                                         LAST ELIM. STEP EXPLICITELY
C                                         BECAUSE R(I+1) IS NOT AVAILABL
       BDYI    = YR*YR+YI*YI
       BDYR    = (BINV(NMAXM)*YR+AUX2(NMAXM)*YI)/BDYI
       BDYI    = (AUX2(NMAXM)*YR-BINV(NMAXM)*YI)/BDYI
       RREAL(NGRID) = EZW*RIM1+FSX*RREAL(NGRID)
     &               - (XR*BDYR-XI*BDYI)
       RIMAG(NGRID) = -XR*BDYI-XI*BDYR
       UR      = POTEFF(NGRID)
       UI      = AUX3(NGRID)
       ZR      = -OFFDIG+EZW*UR
       ZI      = EZW*UI
       BINV(NGRID) = ZR
       AUX2(NGRID) = ZI
       YR      = OFFDG2+FSX*UR-(ZR*BDYR-ZI*BDYI)
       YI      = FSX*UI-(ZR*BDYI+ZI*BDYR)
       AINV(NGRID) = YR
       AUX1(NGRID) = YI
c
c      determinant for testing
c
       detr   = one
       deti   = zero
       do 9010 ir=1,ngrid
       x      = detr*ainv(ir)-deti*aux1(ir)
       deti   = detr*aux1(ir)+deti*ainv(ir)
       detr   = x
 9010  continue
ctest       write(6,'(a,2g12.4))') ' Determ.=',detr,deti
C
C      IN THE THIRD PASS RESOLVE THE REMAINING UPPER TRIANGULAR EQUATION
C      BACKWARD.
C
       DENOM  = AINV(NGRID)*AINV(NGRID)+AUX1(NGRID)*AUX1(NGRID)
       XR     = (RREAL(NGRID)*AINV(NGRID)+RIMAG(NGRID)*AUX1(NGRID))
     &         /DENOM
       XI     = (RIMAG(NGRID)*AINV(NGRID)-RREAL(NGRID)*AUX1(NGRID))
     &         /DENOM
       RREAL(NGRID) = XR*SQM(NGRID)
       RIMAG(NGRID) = XI*SQM(NGRID)
       DO 15 I=NMAXM,2,(-1)
       YR     = RREAL(I) -(XR*BINV(I+1)-XI*AUX2(I+1))
       YI     = RIMAG(I) -(XI*BINV(I+1)+XR*AUX2(I+1))
       DENOM  = AINV(I)*AINV(I)+AUX1(I)*AUX1(I)
       XR     = (YR*AINV(I)+YI*AUX1(I))/DENOM
       XI     = (YI*AINV(I)-YR*AUX1(I))/DENOM
       RREAL(I)   = XR*SQM(I)
       RIMAG(I)   = XI*SQM(I)
   15  CONTINUE
       RREAL(1)   = ZERO
       RIMAG(1)   = ZERO
C
       RETURN
       END
C-----CONRPA-----------------------------------------------------
C
      SUBROUTINE CONRPA
C
      INCLUDE 'PARAMS.INC'

C
C     CONRPA COMPUTES MULTIPOLE STRENGTH FUNCTION IN TE CONTINUUM
C     BY LINEAR RESPONSE:
C
      INCLUDE 'SKYPAR.INC'
      INCLUDE 'WAVFUN.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'RESULT.INC'
      INCLUDE 'NUMPAR.INC'
      INCLUDE 'CONSTN.INC'
      INCLUDE 'LINRES.INC'
C
      DIMENSION NPARSV(KRPMAX),NHOLSV(KRPMAX),ISOPHS(KRPMAX)
C
      CHARACTER*16 TYPMOD
      CHARACTER*12 JSELEC
      CHARACTER*50 TOPER
C
      DATA IFOP14,IFOP17/0,0/
C
C  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
C
C      CHECK BOUNDS
C
       IF(KRPMAX.GT.62) STOP ' KRPMAX TOO LARGE'
C
C      READ PARAMETERS FOR 'SUMRUL'
C
       READ(05,*)   LRPA,RESFAC,CMFAC
       READ(05,*)   ITRPA,EPSRPA,DABRPA
       READ(05,*)   ESTART,ESTEP,EEND,ENWID,ERPRDM,ERPIDM
       NUMENR = (EEND-ESTART)/ESTEP
       CALL SETRED
C
C      INITIALISE SOME FIELDS WITH ZEROS
C
       WRITE(6,'(A)') ' CONRPA ENTERED '
C
C     OPEN FILE 14 FOR CONCENTRATED OUTPUT OF MATRIX ELEMENTS
C
      IF(IFGRHO.NE.0)
     &  OPEN(UNIT=17,STATUS='SCRATCH',FORM='UNFORMATTED',
     &      ACCESS='DIRECT',RECL=5*KVARLN*NGRID2)
      IF(IFOP14.EQ.0) THEN
       OPEN(UNIT=14,STATUS='UNKNOWN',FORM='FORMATTED',FILE='FOR014')
       IFOP14 = 1
      ENDIF
C
C     WRITE A HEADER
C
      WRITE(14,'(A/1X,2A/1X,2(A,I4),5X,2(A,I4))')
     &  '1',
     &  '  TITLE:',TITLE,
     &  'NUCLEUS:  NPROT=',NPROT,'  NNEUT=',NNEUT,
     &  ' SHELLS:  NMAXP=',NMAXP,'  NAMXN=',NMAXN
      WRITE(14,'(1X,A,5G13.6/7X,A,4G13.6/7X,A,G13.6)')
     &  'FORCE: T0-T4=',T0,T1,T2,T3,T4,
     &        ' X0-X3=',X0,X1,X2,X3,
     &        ' POWER=',POWER
      Write(14,'(7x,3(a,g13.6)/)')
     &        ' TYUK=',tyuk,' XYUK=',xyuk,' AYUK=',ayuk
      WRITE(14,'(1X,2(A,F10.3)/18X,2(A,F10.3))')
     &  'SOME PROPERTIES: RMS-RADIUS=',RMS,'    DIFFR.RADIUS=',RDMS,
     &  '   SURFACE=',SURF,'          ENERGY=',ETOTAL
c
c     cancel mass and potential for testing
c
c      do 9010 ir=1,ngrid
c      pot(ir) = zero
c      bmass(ir) = amfc(1)
c      ira = ngrid+ir
c      pot(ira) = zero
c      bmass(ira) = amfc(2)
c 9010 continue
C
      CALL DERV17(BMASS(1),D1B(1),ONE,NGRID,RSTEP)
      CALL DERV17(BMASS(NGRID+1),D1B(NGRID+1),ONE,NGRID,RSTEP)
c
c     provide derivative of "rho"
c
      CALL DERV17(RHO(1),D1R(1), ONE,NGRID,RSTEP)
      CALL DERV17(RHO(ngrid+1),D1R(ngrid+1), ONE,NGRID,RSTEP)
C
C     WRITE OPTIONS AND PARAMETERS FOR RPA PART
C
      WRITE(6,'(/1X,A,I3,3(A,g10.3))')
     &   ' Cont.-RPA WITH: LRPA=',LRPA,'  Res.Int.Fac=',RESFAC,
     &   '  CmInv.Fac=',CMFAC
      WRITE(14,'(/1X,A,I3,3(A,g10.3))')
     &   ' Cont.-RPA WITH: LRPA=',LRPA,'  Res.Int.Fac=',RESFAC,
     &   '  CmInv.Fac=',CMFAC
      WRITE(6,'(21X,A,I5,A,F7.3,A,G12.4)')
     &    '  ITRPA=',ITRPA,' EPSRPA=',EPSRPA,' DAPRPA=',DABRPA
      Write(6,'(6(a,f8.3))')
     &    '  ESTART=',ESTART,'  ESTEP=',ESTEP,'  EEND=',EEND,
     &    '  ENWID=',ENWID,'  ERPRDM=',ERPRDM,'  ERPIDM=',ERPIDM
      Write(14,'(6(a,f8.3))')
     &    '  ESTART=',ESTART,'  ESTEP=',ESTEP,'  EEND=',EEND,
     &    '  ENWID=',ENWID,'  ERPRDM=',ERPRDM,'  ERPIDM=',ERPIDM
C
       IFCMPL = 1
       ALRPA  = LRPA
       MLRPA  = 2*LRPA+1
C
       WRITE(6,'(/1X,A,I3,A,I3)')
     &   'LINEAR RESPONSE FOR RPA MODE WITH L=',LRPA
       WRITE(14,'(/1X,A,I3,A,I3)')
     &   'LINEAR RESPONSE FOR RPA MODE WITH L=',LRPA
C
       CALL BUFFQL(QLROUT,'C','      ')
C
C      read type of external operator in RPA-standard form,
C      and set-up he operator in QINIT
C
      READ(5,*) TYPMOD,ISORPA,JSELEC,ISPLIT,IHWGT,IH0WGT,IPWGT,IPRSHL
      write(6,'(a,i5,a,5i5)')
     &   TYPMOD,ISORPA,JSELEC,ISPLIT,IHWGT,IH0WGT,IPWGT,IPRSHL
      CALL QINIT(TYPMOD,JSELEC,ISPLIT,IHWGT,IH0WGT,IPWGT,
     &           IPRSHL,IMODE,IFPRIM,IFLAGR,XKQRD,RESFAC,CMFAC,
     &           NPARSV,NHOLSV,ISOPHS,TOPER)
      call buffql(qlrout,'O',' QINIT')
C
       DO 888 IENR=0,NUMENR
       ENEXT  = ESTART+ESTEP*IENR
       CALL EXTRSP(ENEXT,ENWID,ERPRDM,ERPIDM,RESFAC,CMFAC,OVQP)
       WRITE(6,'(3(A,G12.4))')
     &   ' ENERGY=',ENEXT,',  WIDTH=',ENWID,
     &   ' YIELDS <Q H**-1 Q>=',OVQP
       WRITE(14,'(3(A,G12.4))')
     &   ' ENERGY=',ENEXT,',  WIDTH=',ENWID,
     &   ' YIELDS <Q H**-1 Q>=',OVQP
  888  CONTINUE
C
      RETURN
      END
C-----EXTRSP------------------------------------------------------
C
      SUBROUTINE EXTRSP(ENEXT,ENWID,ERPRDM,ERPIDM,RESFAC,CMFAC,OVQP)
C
      INCLUDE 'PARAMS.INC'
C
C     COMPUTES RESPONSE TO EXTERNAL FIELD 'Q' GIVEN ON 'QLROUT'
C     WITH COMPLEX ENERGY 'ENEXT'-i*'ENWID', THUS SOLVING ITERATIVELY
C
C            [H-'ENEXT'+i*'ENWID',P] = Q
C
C     THE RESPONSE 'P' IS HANDLED ON 'PLROUT' (REAL PART) AND
C     'SAVLR' (IMAGINARY PART).
C
C     THE ITERATION PROCEEDS AS
C
C       P := (H_hf-'ENEXT'+i*'ENWID')**(-1) [ Q - H_P |PHI> ]
C
C     THE 'ERPDMP' IS AN ADDITIONAL ENERGY WHICH MAY SERVE TO STABILIZE
C     THE ITERATIONS (COMPARABLE TO 'E0DAMP' IN THE H.F. CASE). FASTEST
C     CONVERGENCE (IF STABLE) IS ACHIEVED FOR 'ERPDMP'=0.0 .
C
      INCLUDE 'WAVFUN.INC'
      INCLUDE 'POTINV.INC'
      INCLUDE 'SKYPAR.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'RESULT.INC'
      INCLUDE 'LINRES.INC'
      INCLUDE 'NUMPAR.INC'
      INCLUDE 'CONSTN.INC'
C
C-------------------------------------------------------------------------------
C
C     CHECK OPTION AND SET PARAMETERS
C
      IF(ITRPA.LE.0) RETURN
C
      call buffql(qlrout,'I',' QINIT')
      IFPRJ  = 1
      ADDACT = ONE-EPSRPA
      QWID   = QNORM(QLROUT)
      write(6,'(a,g12.4)') ' INVMOM: QWID=',qwid
c
c     initialise PLROUT , REAL PART ON LABEL 'PR' AND
C     IMAGINARY PART ON LABEL 'PI'.
c
      DO 110 JDIF=1,MLRPA
      DO 110 I=1,NMAXTO*NGRID
  110 SAVLR(I,JDIF) = ZERO
      CALL BUFFQL(SAVLR,'O','RES2QR')
      CALL BUFFQL(SAVLR,'O','RES2QI')
      if(itrpa.gt.2) call DAMPH2(qlrout,savlr,enext,enwid,2)
      CALL BUFFQL(QLROUT,'O','RES2PI')
      CALL BUFFQL(SAVLR,'O','RES2PR')
c
c     test damph2 and invcmp
c
       iftsin = 1
       if(iftsin.eq.1 .and. itrpa.le.1) then
        call copywf(qlrout,plrout)
        if(itrpa.eq.1) then
         call DAMPH2(plrout,savlr,ERPRDM,ERPIDM,1)
         stop
        elseif(itrpa.eq.2) then
         call DAMPH2(plrout,savlr,enext,enwid,1)
         ovqp = qpover(savlr,qlrout)
         ovqpr = qpover(plrout,qlrout)
         write(6,'(3(a,g12.4))')  ' Energy=',enext,
     &    ' Re(R)=',ovqpr,' Im(R)=',ovqp
         write(14,'(1x,f7.3,2g12.4)') enext,ovqpr,ovqp
         CALL BUFFQL(SAVLR,'I','RES2PR')
         return
        endif
       endif
C
C     COMPUTE RIGHT-HAND-SIDE
C
      CALL BUFFQL(SAVLR,'O','RHSPR ')
      CALL COMMHQ(OVQP,OVQP0,OVPP,OVP0P0,IFPRJ,RESFAC,0)
      CALL MULTWF(PLROUT,-ONE)
      CALL BUFFQL(PLROUT,'O','RHSPI ')
      IF(ENWID.EQ.ZERO) THEN
       CALL BUFFQL(SAVLR,'O','RHSQI ')
       CALL MULTWF(QLROUT,-ENEXT)
       CALL BUFFQL(QLROUT,'O','RHSQR ')
      ELSE
       CALL MULTWF(QLROUT,-ENWID)
       CALL BUFFQL(QLROUT,'O','RHSQI ')
       CALL MULTWF(QLROUT,-ENEXT/ENWID)
       CALL BUFFQL(QLROUT,'O','RHSQR ')
      ENDIF
C
C     ITERATE RESPONSE
C
      DO 1000 ITER=1,ITRPA
c
c     set switch for testing. First iteration is only for tests.
c
      if(iter.eq.1) then
ctst   iftest = 1   , --> change test-program before switching on 'iftest'
       iftest = 0
      else
       iftest = 0
      endif
C
C     1. ACTION OF HAMILTONIAN, PRODUCING Q-LIKE REAL PART
C
      CALL BUFFQL(PLROUT,'I','RES2PI')
      CALL COMMHP(OVPQ,OVPQ0,OVQQ,OVQ0Q0,RESFAC,CMFAC,0)
      CALL BUFFQL(SAVLR,'I','RES2QI')
      CALL BUFFQL(PLROUT,'I','RES2QR')
      CALL ADDWF3(PLROUT,SAVLR,QLROUT,-ENEXT, ENWID,-ONE)
      CALL BUFFQL(PLROUT,'O','RESPQR')
C
C     1. ACTION OF HAMILTONIAN, PRODUCING Q-LIKE IMAGINARY PART
C
      CALL BUFFQL(PLROUT,'I','RES2PR')
      CALL COMMHP(OVPQ,OVPQ0,OVQQ,OVQ0Q0,RESFAC,CMFAC,0)
      CALL BUFFQL(SAVLR,'I','RES2QI')
      CALL BUFFQL(PLROUT,'I','RES2QR')
      CALL ADDWF3(PLROUT,SAVLR,QLROUT,-ENWID,-ENEXT, ONE)
      CALL BUFFQL(PLROUT,'O','RESPQI')
C
C     1. ACTION OF HAMILTONIAN, PRODUCING P-LIKE REAL PART
C
      CALL BUFFQL(QLROUT,'I','RES2QI')
      CALL COMMHQ(OVQP,OVQP0,OVPP,OVP0P0,IFPRJ,RESFAC,0)
      CALL BUFFQL(SAVLR,'I','RES2PR')
      CALL BUFFQL(QLROUT,'I','RES2PI')
      CALL ADDWF3(PLROUT,SAVLR,QLROUT,ONE,-ENEXT, ENWID)
      CALL BUFFQL(PLROUT,'O','RESPPR')
C
C     1. ACTION OF HAMILTONIAN, PRODUCING P-LIKE IMAGINARY PART
C
      CALL BUFFQL(QLROUT,'I','RES2QR')
      CALL COMMHQ(OVQP,OVQP0,OVPP,OVP0P0,IFPRJ,RESFAC,0)
      CALL BUFFQL(SAVLR,'I','RES2PR')
      CALL BUFFQL(QLROUT,'I','RES2PI')
      CALL ADDWF3(PLROUT,SAVLR,QLROUT,-ONE,-ENWID,-ENEXT)
      CALL BUFFQL(PLROUT,'O','RESPPI')
C
C     2. ACTION OF HAMILTONIAN, PRODUCING Q-LIKE REAL PART
C
      CALL BUFFQL(PLROUT,'I','RESPPI')
      CALL COMMHP(OVPQ,OVPQ0,OVQQ,OVQ0Q0,RESFAC,CMFAC,0)
      CALL BUFFQL(SAVLR,'I','RESPQI')
      CALL BUFFQL(PLROUT,'I','RESPQR')
      CALL ADDWF3(PLROUT,SAVLR,QLROUT,-ENEXT,-ENWID,-ONE)
C
C     COMPUTE MISMATCH, DAMP IT, AND STEP TO IMPROVED Q-LIKE REAL PART
C
      CALL BUFFQL(SAVLR,'I','RHSQR ')
      CALL ADDWF2(PLROUT,SAVLR,ONE,-ONE)
      ERRBQR = QNORM(PLROUT)
      CALL DAMPH2(PLROUT,SAVLR,ERPRDM,ERPIDM,2)
      ERRDQR = QNORM(PLROUT)
      ERRIQR = QNORM(SAVLR)
      CALL BUFFQL(SAVLR,'I','RES2QR')
      CALL ADDWF2(SAVLR,PLROUT,ONE,-EPSRPA)
      CALL BUFFQL(SAVLR,'O','RES2QR')
C
C     2. ACTION OF HAMILTONIAN, PRODUCING Q-LIKE IMAGINARY PART
C
      CALL BUFFQL(PLROUT,'I','RESPPR')
      CALL COMMHP(OVPQ,OVPQ0,OVQQ,OVQ0Q0,RESFAC,CMFAC,0)
      CALL BUFFQL(SAVLR,'I','RESPQI')
      CALL BUFFQL(PLROUT,'I','RESPQR')
      CALL ADDWF3(PLROUT,SAVLR,QLROUT, ENWID,-ENEXT, ONE)
C
C     COMPUTE MISMATCH, DAMP IT, AND STEP TO IMPROVED Q-LIKE IMAGINARY PART
C
      CALL BUFFQL(SAVLR,'I','RHSQI ')
      CALL ADDWF2(PLROUT,SAVLR,ONE,-ONE)
      ERRBQI = QNORM(PLROUT)
      CALL DAMPH2(PLROUT,SAVLR,ERPRDM,ERPIDM,2)
      ERRDQI = QNORM(PLROUT)
      ERRIQI = QNORM(SAVLR)
      CALL BUFFQL(SAVLR,'I','RES2QI')
      CALL ADDWF2(SAVLR,PLROUT,ONE,-EPSRPA)
      CALL BUFFQL(SAVLR,'O','RES2QI')
C
C     2. ACTION OF HAMILTONIAN, PRODUCING P-LIKE REAL PART
C
      CALL BUFFQL(QLROUT,'I','RESPQI')
      CALL COMMHQ(OVQP,OVQP0,OVPP,OVP0P0,IFPRJ,RESFAC,0)
      CALL BUFFQL(SAVLR,'I','RESPPR')
      CALL BUFFQL(QLROUT,'I','RESPPI')
      CALL ADDWF3(PLROUT,SAVLR,QLROUT, ONE,-ENEXT,-ENWID)
C
C     COMPUTE MISMATCH, DAMP IT, AND STEP TO IMPROVED Q-LIKE REAL PART
C
      CALL BUFFQL(SAVLR,'I','RHSPR ')
      CALL ADDWF2(PLROUT,SAVLR,ONE,-ONE)
      ERRBPR = QNORM(PLROUT)
      CALL DAMPH2(PLROUT,SAVLR,ERPRDM,ERPIDM,2)
      ERRDPR = QNORM(PLROUT)
      ERRIPR = QNORM(SAVLR)
      CALL BUFFQL(SAVLR,'I','RES2PR')
      CALL ADDWF2(SAVLR,PLROUT,ONE,-EPSRPA)
      CALL BUFFQL(SAVLR,'O','RES2PR')
C
C     2. ACTION OF HAMILTONIAN, PRODUCING P-LIKE IMAGINARY PART
C
      CALL BUFFQL(QLROUT,'I','RESPQR')
      CALL COMMHQ(OVQP,OVQP0,OVPP,OVP0P0,IFPRJ,RESFAC,0)
      CALL BUFFQL(SAVLR,'I','RESPPR')
      CALL BUFFQL(QLROUT,'I','RESPPI')
      CALL ADDWF3(PLROUT,SAVLR,QLROUT,-ONE, ENWID,-ENEXT)
C
C     COMPUTE MISMATCH, DAMP IT, AND STEP TO IMPROVED Q-LIKE IMAGINARY PART
C
      CALL BUFFQL(SAVLR,'I','RHSPI ')
      CALL ADDWF2(PLROUT,SAVLR,ONE,-ONE)
      ERRBPI = QNORM(PLROUT)
      CALL DAMPH2(PLROUT,SAVLR,ERPRDM,ERPIDM,2)
      ERRDPI = QNORM(PLROUT)
      ERRIPI = QNORM(SAVLR)
      CALL BUFFQL(SAVLR,'I','RES2PI')
      CALL ADDWF2(SAVLR,PLROUT,ONE,-EPSRPA)
      CALL BUFFQL(SAVLR,'O','RES2PI')
C
C     CHECK CONVERGENCE AND PRINT
C
      WRITE(6,'(A,I4,(A,4G10.2),2(/24X,A,4G10.2))')
     & ' AT ITER=',ITER,
     & '  MISMATCH: BEFORE=',ERRBQR,ERRBQI,ERRBPR,ERRBPI,
     &            ' DAMPED=',ERRDQR,ERRDQI,ERRDPR,ERRDPI,
     &            'D.IMAG.=',ERRIQR,ERRIQI,ERRIPR,ERRIPI
      WRITE(14,'(A,I4,A,4G10.2)')
     & ' AT ITER=',ITER,
     & '  MISMATCH: BEFORE=',ERRBQR,ERRBQI,ERRBPR,ERRBPI
      ALLERR = SQRT(ERRBQR**2+ERRBQI**2+ERRBPR**2+ERRBPI**2)
      IF(ALLERR.LT.DABRPA) GOTO 1009
 1000 CONTINUE
      WRITE(6,'(A,I4,A)')
     & ' WARNING: NO CONVERGENCE WITHIN ITER=',ITER,' STEPS!.'
 1009 CONTINUE
C
C     COMPUTE NORM OF RESPONSE
C
      CALL BUFFQL(QLROUT,'I','RESPQR')
      RESPNO = QNORM(QLROUT)
      CALL BUFFQL(QLROUT,'I','RESPQI')
      RESPNO = QNORM(QLROUT)  + RESPNO
      CALL BUFFQL(QLROUT,'I','RESPPR')
      RESPNO = QNORM(QLROUT)  + RESPNO
      CALL BUFFQL(QLROUT,'I','RESPPI')
      RESPNO = QNORM(QLROUT)  + RESPNO
      RESPNO = SQRT(RESPNO)
      WRITE(6,'(A,G13.4)') '  NORM OF RESPONSE=',RESPNO
      OVQP   = RESPNO
C
      RETURN
      END
ckinC-----DAMPH2------------------------------------------------------
ckinC
ckin      SUBROUTINE DAMPH2(QWF,QIWF,EREAL,EIMAG,NRINV)
ckinC
ckin      INCLUDE 'PARAMS.INC'
ckinC
ckinC     COMPUTES DAMPING BY (H0-E*)**-1 (H0-E)**-1
ckinC     WITH COMPLEX ENERGY 'E' = 'EREAL'-i*'EIMAG'
ckinC
ckin      DIMENSION QWF(KX11,KLRPA)
ckin      DIMENSION QIWF(KX11,KLRPA)
ckin      dimension anorsv(21)
ckinC
ckin      INCLUDE 'WAVFUN.INC'
ckin      INCLUDE 'POTINV.INC'
ckin      INCLUDE 'SKYPAR.INC'
ckin      INCLUDE 'GRID.INC'
ckin      INCLUDE 'RESULT.INC'
ckin      INCLUDE 'LINRES.INC'
ckin      INCLUDE 'NUMPAR.INC'
ckin      INCLUDE 'CONSTN.INC'
ckinC
ckinC-------------------------------------------------------------------------------
ckinC
ckinC     CHECK OPTION AND SET PARAMETERS
ckinC
ckin      iftsin = 0
ckinC
ckinC     PREPARE POTENATIAL- AND MASS-FIELDS FOR kinetic energy INVERSION
ckinC
ckin      averms = SQRT(0.5*(amfc(1)+amfc(2)))
ckin      DO 110 I=1,NGRID
ckin      d2b(i) = zero
ckin      potl(i) = zero
ckin  110 SQB(I) = averms
ckinC
ckinC     PROVIDE NECESSARY QUANTUM NUMBERS AND POINTERS WITH
ckinC     DO LOOPS 10, 20 AND 30
ckinC
ckin      INCLUDE 'INITDO.INC'
ckinC
ckinC     ACTION OF LINEAR RESPONSE HAMILTONIAN, result on 'hlinrs'
ckinC
ckinC     INVERT WITH MEAN-FIELD HAMILTONIAN MINUS EXTERNAL ENERGY
ckinC     AND PROPAGATE TO IMPROVED GUESS FOR 'PLROUT'
ckinC     'QWF' IS USED AS WORKSPACE IN THE INVERSION STEPS
ckinC
ckin      CALL PROJEC(QWF(NADD+1,JDIF),JS,NSTACK,NMAX)
ckin      DO 140 IR=1,NGRID
ckin      WORK1(IR)=QWF(NADD+IR,JDIF)*((IR-1)*RSTEP)
ckin      WORK3(IR)=WORK1(IR)
ckin 140  CONTINUE
ckin      umin = ereal
ckin      CALL invert(WORK1,WORK1,-1)
ckin      IF(NRINV.EQ.1) THEN
ckin       DO 151 IR=2,NGRID
ckin       QWF(IR+NADD,JDIF) = WORK1(IR)*RINVS(IR)
ckin  151  QIWF(IR+NADD,JDIF) = zero
ckin      ELSE
ckin       CALL invert(WORK1,WORK1,1)
ckin       DO 150 IR=2,NGRID
ckin       QIWF(IR+NADD,JDIF) = zero
ckin  150  QWF(IR+NADD,JDIF) = WORK1(IR)*RINVS(IR)
ckin      ENDIF
ckin      PARLL  = (1-MOD(LL,2)*2)
ckin      CALL EXTRAP(QWF(NADD+1,JDIF),PARLL)
ckin      CALL PROJEC(QWF(NADD+1,JDIF),JS,NSTACK,NMAX)
ckin      CALL EXTRAP(QIWF(NADD+1,JDIF),PARLL)
ckin      CALL PROJEC(QIWF(NADD+1,JDIF),JS,NSTACK,NMAX)
ckinC
ckin      INCLUDE 'ENDDO.INC'
ckinC
ckin      RETURN
ckin      END
C-----DAMPH2------------------------------------------------------
C
      SUBROUTINE DAMPH2(QWF,QIWF,EREAL,EIMAG,NRINV)
C
      INCLUDE 'PARAMS.INC'
C
C     COMPUTES DAMPING BY (H0-E*)**-1 (H0-E)**-1
C     WITH COMPLEX ENERGY 'E' = 'EREAL'-i*'EIMAG'
C
      DIMENSION QWF(KX11,KLRPA)
      DIMENSION QIWF(KX11,KLRPA)
      dimension anorsv(21)
C
      INCLUDE 'WAVFUN.INC'
      INCLUDE 'POTINV.INC'
      INCLUDE 'SKYPAR.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'RESULT.INC'
      INCLUDE 'LINRES.INC'
      INCLUDE 'NUMPAR.INC'
      INCLUDE 'CONSTN.INC'
C
C-------------------------------------------------------------------------------
C
C     CHECK OPTION AND SET PARAMETERS
C
      if(itrpa.le.1) then
       iftsin = 1
      else
       iftsin = 0
      endif
C
C     PROVIDE NECESSARY QUANTUM NUMBERS AND POINTERS WITH
C     DO LOOPS 10, 20 AND 30
C
      INCLUDE 'INITDO.INC'
C
C     ACTION OF LINEAR RESPONSE HAMILTONIAN, result on 'hlinrs'
C
C     PREPARE POTENATIAL- AND MASS-FIELDS FOR INVERSION
C
      DO 110 I=1,NGRID
  110 SQB(I) = SQRT(BMASS(I+ISOADD))
C
      CALL DERV25(SQB(1),D2B(1),ONE,NGRID,RSTEP)
      DO 120 I=1,NGRID
  120 D2B(I) = D2B(I)/SQB(I)
C
c      CALL DERV15(BMASS(ISOADD+1),D1B(1),ONE,NGRID,RSTEP)
c      D1B(1)=(16.0*BMASS(2+ISOADD)-15.0*BMASS(1+ISOADD)
c     *        -BMASS(3+ISOADD))/(6.0*RSTEP**RSTEP)
c      DO 130 I=2,NGRID
c  130 D1B(I) = D1B(I)*RINVS(I)
C
C      THE EFFECTIVE POTENTIAL FOR GIVEN 'JSA' AND 'ISO' ON 'POTL';
C
C      'PHASJS' IS THE DEGENERACY OF THE SHELL.
C      'JSNPL' IS POINTER TO TO NUCLEON SHELL IN POINTER FIELD 'NPLA'.
C
      XLP    = FLOAT(LL*LL+LL)
      XL0    = WLS(JS)
C
      D1B0=(16.0*BMASS(2+ISOADD)-15.0*BMASS(1+ISOADD)
     *        -BMASS(3+ISOADD))/(6.0*RSTEP**RSTEP)
      POTL(1) = POT(1+ISOADD)+XL0*POTLS(1+ISOADD)+D1B0
      DO 160 I=2,NGRID
      N      = I+ISOADD
      POTL(I) = POT(N)+XL0*POTLS(N)+XLP*BMASS(N)*ZENTRF(I)
     &          +D1B(I+isoadd)*rinvs(i)
  160 CONTINUE
C
C     INVERT WITH MEAN-FIELD HAMILTONIAN MINUS EXTERNAL ENERGY
C     AND PROPAGATE TO IMPROVED GUESS FOR 'PLROUT'
C     'QWF' IS USED AS WORKSPACE IN THE INVERSION STEPS
C
c     BNORM  = ZERO
      CALL PROJEC(QWF(NADD+1,JDIF),JS,NSTACK,NMAX)
      DO 140 IR=1,NGRID
      WORK1(IR)=QWF(NADD+IR,JDIF)*((IR-1)*RSTEP)
      WORK3(IR)=WORK1(IR)
      WORK2(IR)=WORK1(IR)
c     BNORM  = WORK1(IR)*WORK1(IR)+BNORM
 140  CONTINUE
ct1        CALL HPSI(WORK2(1),work1(1),AUX2(1),AUX3(1),1,JS,ISO,
ct1     &            SUBFAC+EREAL)
      CALL INVCMP(WORK1,WORK2,WORK1,SUBFAC+EREAL,EIMAG,LL)
c      CALL PROJEC(work1(1),JS,NSTACK,NMAX)
c      CALL PROJEC(work2(1),JS,NSTACK,NMAX)
      IF(NRINV.EQ.1) THEN
       DO 151 IR=2,NGRID
c      ANORM  = (WORK1(IR)+WORK2(IR))**2+ANORM
       QWF(IR+NADD,JDIF) = WORK1(IR)*RINVS(IR)
  151  QIWF(IR+NADD,JDIF) = WORK2(IR)*RINVS(IR)
       if(iftsin.eq.1) then
         write(6,'(a,i5,a,4i5)')
     &    ' test INVCMP for state nr.',NA,
     &    ' JS,LL,JSA,LA=',js,ll,jsa,la
        CALL HPSI(WORK1(1),AUX1(1),AUX2(1),AUX3(1),1,JS,ISO,
     &            SUBFAC+EREAL)
        anorm  = zero
        onorm  = zero
        oaovlp = zero
        DO 710 IR=1,NGRID
ct1        aux1(ir) = work1(ir)
c        AUX1(IR) = AUX1(IR)+EIMAG*WORK2(IR)
        anorm  = aux1(ir)*aux1(ir) + anorm
        onorm  = work3(ir)*work3(ir) + onorm
        oaovlp = aux1(ir)*work3(ir) + oaovlp
  710   CONTINUE
        cosang = oaovlp*oaovlp/max(anorm*onorm,xsmall)
        write(6,'(a,3g11.4,A,G11.4)')
     &     ' norms: a,o,oa=',anorm,onorm,oaovlp,' cosang=',cosang
        write(6,'(a,99(/1x,f8.3,3g12.4))')
     &   '   R       before       after        error',
     &   ((ir-1)*rstep,work3(ir),aux1(ir),
     &     (work3(ir)-aux1(ir))/work3(ir),ir=2,ngrid)
        CALL HPSI(WORK2(1),AUX1(1),AUX2(1),AUX3(1),1,JS,ISO,
     &            SUBFAC+EREAL)
c        CALL HPSI(WORK2(1),AUX1(1),AUX2(1),AUX3(1),1,JS,ISO,
c     &            SUBFAC+EREAL)
        anorm  = zero
        DO 720 IR=1,NGRID
c       AUX1(IR) = AUX1(IR)-EIMAG*WORK1(IR)
        anorm  = (AUX1(IR)-EIMAG*WORK1(IR))**2 + anorm
  720   CONTINUE
        write(6,'(a,g12.4)')
     &     ' norm imaginary part:',anorm
c        write(6,'(a,99(/1x,f8.3,3g12.4))')
c     &   '   R        work1       work2        aux1',
c     &   ((ir-1)*rstep,work1(ir),work2(ir),aux1(ir),ir=1,ngrid)
       endif
      ELSE
       CALL INVCMP(WORK1,QIWF(NADD+1,JDIF),WORK1,
     &             SUBFAC+EREAL,-EIMAG,LL)
       CALL INVCMP(WORK3,WORK2,WORK2,SUBFAC+EREAL,-EIMAG,LL)
c      ANORM  = ZERO
       DO 150 IR=2,NGRID
c      ANORM  = (WORK1(IR)+WORK2(IR))**2+ANORM
       QIWF(IR+NADD,JDIF) =
     &                     (QIWF(IR+NADD,JDIF)+WORK3(IR))*RINVS(IR)
       WORK3(IR)=QWF(NADD+IR,JDIF)*((IR-1)*RSTEP)
  150  QWF(IR+NADD,JDIF) = (WORK1(IR)-WORK2(IR))*RINVS(IR)
c      WRITE(6,'(A,3G12.4)')
c    &  ' IN DAMPH2: BEFORE,AFTER,E**-2=',BNORM,ANORM,
c    &  BNORM/MAX(ANORM,XSMALL)
       if(iftsin.eq.1) then
        DO 240 IR=1,NGRID
        WORK1(IR)=QWF(NADD+IR,JDIF)*((IR-1)*RSTEP)
 240    CONTINUE
        CALL HPSI3(WORK1(1),work2(1),AUX2(1),AUX3(1),1,JS,ISO,
     &            SUBFAC+EREAL)
        CALL HPSI3(work2(1),aux1(1),AUX2(1),AUX3(1),1,JS,ISO,
     &            SUBFAC+EREAL)
c        CALL HPSI(WORK1(1),work2(1),AUX2(1),AUX3(1),1,JS,ISO,
c     &            SUBFAC+EREAL)
c        CALL HPSI(WORK2(1),AUX1(1),AUX2(1),AUX3(1),1,JS,ISO,
c     &            SUBFAC+EREAL)
        ainorm = zero
        anorm  = zero
        onorm  = zero
        oaovlp = zero
        DO 210 IR=2,NGRID
        AUX1(IR) = aux1(IR)+(EIMAG*EIMAG)*WORK3(IR)
        anorm  = aux1(ir)*aux1(ir) + anorm
        onorm  = work3(ir)*work3(ir) + onorm
        oaovlp = aux1(ir)*work3(ir) + oaovlp
        ainorm =(QIWF(NADD+IR,JDIF)*((IR-1)*RSTEP))**2 + ainorm
  210   CONTINUE
        cosang = oaovlp*oaovlp/max(anorm*onorm,xsmall)
        if(anorm.gt.zero) then
         write(6,'(a,i5,a,4i5)')
     &    ' test INVCMP*INVCMP for state nr.',NA,
     &    ' JS,LL,JSA,LA=',js,ll,jsa,la
         write(6,'(a,3g11.4,A,G11.4)')
     &      ' norms: a,o,oa=',anorm,onorm,oaovlp,' cosang=',cosang
         write(6,'(a)')
     &    '   R     before       after       raw1         raw2'
         write(6,'(a,g11.4)') ' norm of imag. part=',ainorm
         do 9020 ir=1,ngrid
         R      = (ir-1)*rstep
 9020    write(6,'(1x,f5.2,4g12.4)')
     &       R,work3(ir),aux1(ir),work1(ir),work2(ir)
        endif
       endif
      ENDIF
      PARLL  = (1-MOD(LL,2)*2)
      CALL EXTRAP(QWF(NADD+1,JDIF),PARLL)
      CALL PROJEC(QWF(NADD+1,JDIF),JS,NSTACK,NMAX)
      CALL EXTRAP(QIWF(NADD+1,JDIF),PARLL)
      CALL PROJEC(QIWF(NADD+1,JDIF),JS,NSTACK,NMAX)
C
      INCLUDE 'ENDDO.INC'
C
      RETURN
      END
coldC-----EXTOLD------------------------------------------------------
coldC
cold      SUBROUTINE EXTOLD(ENEXT,ENWID,OVQP)
coldC
cold      INCLUDE 'PARAMS.INC'
coldC
coldC     COMPUTES RESPONSE TO EXTERNAL FIELD 'Q' GIVEN ON 'QLROUT'
coldC     WITH COMPLEX ENERGY 'ENEXT'-i*'ENWID', THUS SOLVING ITERATIVELY
coldC
coldC            [H-'ENEXT'+i*'ENWID',P] = Q
coldC
coldC     THE RESPONSE 'P' IS HANDLED ON 'PLROUT' (REAL PART) AND
coldC     'SAVLR' (IMAGINARY PART).
coldC
coldC     THE ITERATION PROCEEDS AS
coldC
coldC       P := (H_hf-'ENEXT'+i*'ENWID')**(-1) [ Q - H_P |PHI> ]
coldC
cold      INCLUDE 'WAVFUN.INC'
cold      INCLUDE 'POTINV.INC'
cold      INCLUDE 'SKYPAR.INC'
cold      INCLUDE 'GRID.INC'
cold      INCLUDE 'RESULT.INC'
cold      INCLUDE 'LINRES.INC'
cold      INCLUDE 'NUMPAR.INC'
cold      INCLUDE 'CONSTN.INC'
coldC
coldC-------------------------------------------------------------------------------
coldC
coldC     CHECK OPTION AND SET PARAMETERS
coldC
cold      IF(ITRPA.LE.0) RETURN
coldC
cold      ADDACT = ONE-EPSRPA
cold      QWID   = QNORM(QLROUT)
cold      write(6,'(a,g12.4)') ' INVMOM: QWID=',qwid
coldc
coldc     initialise PLROUT
coldc
cold      CALL COPYWF(QLROUT,PLROUT)
coldC
coldC     ITERATE
coldC
cold      OVQPOL = ZERO
cold      DO 1000 ITER=1,ITRPA
coldC
coldC     PREPARE THE RESPONSE CURRENT 'QCUR' AND ITS DERIVATIVE 'DQCUR'
coldC
cold      CALL PDENSI(one)
coldC
coldC     PREPARE THE KINETIC RESPONSE POTENTIAL FOR GIVEN CURRENT
coldC
cold      CALL PPOT(one,zero)
coldC
coldC     PROVIDE NECESSARY QUANTUM NUMBERS AND POINTERS WITH
coldC     DO LOOPS 10, 20 AND 30
coldC
cold      OVQP   = ZERO
cold      OVPPR  = ZERO
cold      OVPPI  = ZERO
cold      OVPP   = ZERO
coldC
cold      INCLUDE 'INITDO.INC'
coldC
coldC     ACTION OF LINEAR RESPONSE HAMILTONIAN, result on 'hlinrs'
coldC
cold      CALL HRESPP(ISO,LA,XJA,NADD,XJ,XL,JDIF,WF(NADD+1),WFS(NADD+1))
cold      DO 40 IR=2,NGRID
cold   40 HLINRS(IR) = ((IR-1)*RSTEP)*(QLROUT(NADD+IR,JDIF)-HLINRS(IR))
cold      HLINRS(1)  = ZERO
coldC
coldC     PREPARE POTENATIAL- AND MASS-FIELDS FOR INVERSION
coldC
cold      DO 110 I=1,NGRID
cold  110 SQB(I) = SQRT(BMASS(I+ISOADD))
coldC
cold      CALL DERV25(SQB(1),D2B(1),ONE,NGRID,RSTEP)
cold      DO 120 I=1,NGRID
cold  120 D2B(I) = D2B(I)/SQB(I)
coldC
cold      CALL DERV15(BMASS(ISOADD+1),D1B(1),ONE,NGRID,RSTEP)
cold      D1B(1)=(16.0*BMASS(2+ISOADD)-15.0*BMASS(1+ISOADD)
cold     *        -BMASS(3+ISOADD))/(6.0*RSTEP**RSTEP)
cold      DO 130 I=2,NGRID
cold  130 D1B(I) = D1B(I)*RINVS(I)
coldC
coldC      THE EFFECTIVE POTENTIAL FOR GIVEN 'JSA' AND 'ISO' ON 'POTL';
coldC
coldC      'PHASJS' IS THE DEGENERACY OF THE SHELL.
coldC      'JSNPL' IS POINTER TO TO NUCLEON SHELL IN POINTER FIELD 'NPLA'.
coldC
cold      XLP    = FLOAT(Ll*Ll+Ll)
cold      XL0    = WLS(JS)
coldC
cold      POTL(1) = POT(1+ISOADD)+XL0*POTLS(1+ISOADD)+D1B(1)
cold      DO 160 I=2,NGRID
cold      N      = I+ISOADD
cold      POTL(I) = POT(N)+XL0*POTLS(N)+XLP*BMASS(N)*ZENTRF(I)+D1B(I)
cold  160 CONTINUE
coldC
coldC     INVERT WITH MEAN-FIELD HAMILTONIAN MINUS EXTERNAL ENERGY
coldC     AND PROPAGATE TO IMPROVED GUESS FOR 'PLROUT'
coldC
cold      CALL INVCMP(WORK1,WORK2,HLINRS,SUBFAC+ENEXT,-ENWID)
cold      DO 150 IR=2,NGRID
cold      PLROUT(IR+NADD,JDIF) = ADDACT*PLROUT(IR+NADD,JDIF)
cold     &                      +EPSRPA*WORK1(IR)*RINVS(IR)
cold  150 SAVLR(IR+NADD,JDIF)  = WORK2(IR)*RINVS(IR)
cold      PARLL  = (1-MOD(LL,2)*2)
cold      CALL EXTRAP(PLROUT(NADD+1,JDIF),PARLL)
cold      CALL EXTRAP(SAVLR(NADD+1,JDIF),PARLL)
cold      CALL PROJEC(PLROUT(NADD+1,JDIF),JS,NSTACK,NMAX)
cold      CALL PROJEC(SAVLR(NADD+1,JDIF),JS,NSTACK,NMAX)
coldC
cold      OVQP   = OVERLP(QLROUT(NADD+1,JDIF),PLROUT(NADD+1,JDIF),
cold     &                RSTEP,NGRID)    +  OVQP
cold      PPR    = OVERLP(PLROUT(NADD+1,JDIF),PLROUT(NADD+1,JDIF),
cold     &                RSTEP,NGRID)
cold      PPI    = OVERLP(SAVLR(NADD+1,JDIF),SAVLR(NADD+1,JDIF),
cold     &                RSTEP,NGRID)
cold      OVPPR  = PPR + OVPPR
cold      OVPPI  = PPI + OVPPI
cold      OVPP   = PPR+PPI + OVPP
coldC
cold      INCLUDE 'ENDDO.INC'
coldC                                      PRINTS
cold      EM1    = QWID/OVQP
cold      WRITE(6,'(1X,A,I3,A,4G11.3,A,G12.4)')
cold     &   'INV-IT=',ITER,' OVQP,OVPP,OVPPR,OVPPI=',OVQP,OVPP,OVPPR,OVPPI,
cold     &   ' E(-1)=',EM1
cold      IF(ITER.GE.2 .AND. ABS(OVQP-OVQPOL).LT.DABRPA) GOTO 1009
cold      OVQPOL = OVQP
coldC
cold 1000 CONTINUE
cold      WRITE(6,'(A,A,I5,A)') ' INVERSE MOMENTUM NOT FULLY CONVERGED',
cold     &   ' WITHIN',I4,' ITERATIONS!'
cold 1009 CONTINUE
coldC
cold      WRITE(6,'(2(A,G12.4))')
cold     &  ' AFTER INVMRQ: OVQP=',OVQP,'  OVPP=',OVPP
coldC
cold      RETURN
cold      END
ckinC-----INVKCM------------ PART OF THE QHQ PACKAGE -------------------------------
ckinC
ckin      SUBROUTINE INVKCM(RREAL,RIMAG,RHS,BMASS0,EINVR,EINVI)
ckinC
ckin      INCLUDE 'PARAMS.INC'
ckinC
ckin      DIMENSION RREAL(KG1),RIMAG(KG1),RHS(KG1)
ckinC
ckin      INCLUDE 'WAVFUN.INC'
ckin      INCLUDE 'POTINV.INC'
ckin      INCLUDE 'GRID.INC'
ckin      INCLUDE 'CONSTN.INC'
ckinC
ckin       DATA FSX,EZW/.83333333333333,.8333333333333E-1/
ckinC
ckinC
ckinC      DETERMINES 'R' AS SOLUTION OF THE INHOMOGENOUS LINEAR DIFF. EQ.:
ckinC       " (-D/DX BMASS0 D/DX + ( - EINV ) ) R = RHS  "
ckinC      FOR GIVEN 'RHS' WHERE 'BMASS' IS X-DEPENDENT INVERSE MASS,
ckinC      AND 'EINV' A SHIFT CONSTANT.
ckinC      THE ENERGY-SHIFT 'EINV' IS COMPLEX WITH REAL PART 'EINVR' AND
ckinC      IMAGINARY PART 'EINVI'.
ckinC      ACCORDINGLY THE RESPONSE FIELD 'R' IS COMPLEX WITH REAL PART 'RREAL'
ckinC      AND IMAGINARY PART 'RIMAG'.
ckinC
ckinC      THE ARRAYS IN THE CALLING LIST MEAN:
ckinC       RREAL  = REAL PART OF THE THE SOLUTION OF THE EQUATION
ckinC       RIMAG  = IMAGINARY PART OF THE THE SOLUTION OF THE EQUATION
ckinC       RHS    = THE RIGHT HAND SIDE (INHOMOGENEITY) OF THE EQUATION.
ckinC                R AND RHS MAY BE THE SAME ARRAYS IN THE CALLING ROUTINE
ckinC                RHS IS THEN OVERWRITTEN.
ckinC      THE PARAMETERS IN THE CALLING LIST MEAN:
ckinC       EINVR  = REAL PART OF THE ENERGY-SHIFT
ckinC       EINVI  = IMAGINARY PART OF THE ENERGY-SHIFT
ckinC
ckinC      INFORMATION ON THE HAMILTONIAN IS TRANSFERRED VIA COMMON /POTINV/:
ckinC       SQB    = INVERSE SQUARE-ROOT OF X-DEPENDENT EFFECTIVE MASS
ckinC       POTL   = EFFECTIVE POTENTIAL BEFORE MASS TRANSFORMATION,
ckinC                TO BE PREPARED IN SUBROUTINE 'WELFU'.
ckinC
ckinC      INFORMATION ON THE GRID IS TRANSFERRED VIA COMMON /GRID/.
ckinC
ckinC      NOTE: THE PROGRAM IS WRITTEN FOR SOLVING THE SCHROEDINGER EQUATIO
ckinC            IN RADIALLY SYMMETRIC PROBLEMS. THUS IT IMPLIES THE RHS AND
ckinC            POTEFF TO BE SUCH THAT R(X) VANISHES AT THE ORIGIN.
ckinC
ckinC      THE SOLUTION PROCEEDS AS FOLLOWS:
ckinC
ckinC      BY INTRODUCING SCALED WAVE FUNCTIONS
ckinC       " REFF = R/SQB "
ckinC      AND
ckinC       " RHSEFF = RHS/SQB ".
ckinC      THE EQUATIONS IS TRANSFORMED TO A LAPLACIAN EQUATION
ckinC       " (-D/DX**2 + POTEFF(X) ) REFF = RHSEFF "
ckinC      WITH CONSTANT MASS AND AN EFFECTIVE POTENTIAL
ckinC       " POTEFF = (POTL(X)-UMIN)/(SQB(X)*SQB(X)) + D2B(X) "
ckinC      WHERE
ckinC       " D2B = (D/DX**2 SQB)/SQB " .
ckinC      AFTER SOLVING THE LAPLACIAN EQUATION, THE SOLUTION IS TRANSFORMED
ckinC      BACK BY
ckinC       " R   = REFF*SQB" .
ckinC
ckinC      THE D/DX**2 IS TREATED EFFECTIVELY IN 5-POINT PRECISION.
ckinC      IN ADDITION TO A STANDARD 3-POINT FORMULA FOR D/DX**2
ckinC      A COMPLEMENTING 3-POINT AVERAGING OF UEFF*REFF AND RHSEFF
ckinC      IS USED WHICH BOTH TOGETHER MAKE UP THE 5-POINT PRECISION.
ckinC      (SEE R.ZURMUEHL, "PRAKTISCHE MATHEMATIK", CH. VII, SECT. 3,
ckinC       SPRINGER VERLAG, BERLIN 1963)
ckinC
ckinC
ckinC      THE WORKSPACES FOR THE HAMILTONIAN, 'POTEFF', 'AINV' AND 'BINV' ARE
ckinC      NOW BECOMING COMPLEX. WE AVOID EXPLICIT COMPLEX ARITHMETIC
ckinC      AND ASSIGN THE WORKSPACES 'AUX..' AS FOLLOWS:
ckinC         AUX1   = IMAGINARY PART OF 'AINV'
ckinC         AUX2   = IMAGINARY PART OF 'BINV'
ckinC         AUX3   = IMAGINARY PART OF 'POTEFF'
ckinC
ckinC-----------------------------------------------------------------------
ckinC
ckinC      SQUARE-ROOT OF INVERSE EFF. MASS ON 'SQB',
ckinC      SECOND DERIVATIVE OF 'SQB' ON 'D2B'.
ckinC      FINALLY ASSEMBLING  THE EFFECTIVE POTENTIAL ON 'POTEFF'.
ckinC
ckin       SQB0   = SQRT(BMASS0)
ckin       SQM0   = ONE/SQB0
ckin       XMASS0 = SQM0*SQM0
ckin       DO 80 I=2,NGRID
ckin       POTEFF(I) = D2B(I)+(POTL(I)-EINVR)*XMASS0
ckin       AUX3(I)   = (-EINVI)*XMASS0
ckin   80  CONTINUE
ckinC
ckinC      TRANSFORM TO EFFECTIVE R.H.S., STORE RESULT ON R.
ckinC
ckin       DO 17 I=1,NGRID
ckin   17  RREAL(I)   = RHS(I)*SQM0
ckinC
ckinC     PREPARE PARAMETERS FOR KINETIC ENERGY AND LOOP-ENDS.
ckinC
ckin       OFFDIG = ONE/(RSTEP*RSTEP)
ckin       OFFDG2 = OFFDIG+OFFDIG
ckin       NMAXM  = NGRID-1
ckinC
ckinC      BEGIN FIRST PASS ELIMINATING THE LOWER OFF-DIAGONAL ELEMENT.
ckinC      NEW DIAGONAL AND UPPER OFF-DIAGONAL ELEMNTS ARE SAVED ON A AND B.
ckinC      THE SECOND PASS, FORWARD ELIMINATION OF THE R.H.S., IS DONE WITHI
ckinC      THE SAME LOOP; RESULT IS STORED ON FIELD T.
ckinC
ckinC
ckin       BINV(1) = ZERO
ckin       AUX1(1) = ZERO
ckin       AINV(1) = ONE
ckin       AUX2(1) = ZERO
ckin       XR      = ZERO
ckin       XI      = ZERO
ckin       YR      = ONE
ckin       YI      = ZERO
ckin       RIM1    = RREAL(1)
ckinC
ckin       DO 13 I=2,NMAXM
ckinC                                         FORWARD ELIM. FOR R.H.S.
ckin       BDYI    = YR*YR+YI*YI
ckin       BDYR    = (BINV(I-1)*YR+AUX2(I-1)*YI)/BDYI
ckin       BDYI    = (AUX2(I-1)*YR-BINV(I-1)*YI)/BDYI
ckin       XR      = EZW*(RREAL(I+1)+RIM1)+FSX*RREAL(I)
ckin     &          - (XR*BDYR-XI*BDYI)
ckin       XI      = -XR*BDYI-XI*BDYR
ckin       RIM1    = RREAL(I)
ckin       RREAL(I)= XR
ckin       RIMAG(I)= XI
ckinC                                         FORWARD ELIM. FOR MATRIX
ckin       UR      = POTEFF(I)
ckin       UI      = AUX3(I)
ckin       ZR      = -OFFDIG+EZW*UR
ckin       ZI      = EZW*UI
ckin       BINV(I) = ZR
ckin       AUX2(I) = ZI
ckin       YR      = OFFDG2+FSX*UR-(ZR*BDYR-ZI*BDYI)
ckin       YI      = FSX*UI-(ZR*BDYI+ZI*BDYR)
ckin       AINV(I) = YR
ckin       AUX1(I) = YI
ckin   13  CONTINUE
ckinC                                         LAST ELIM. STEP EXPLICITELY
ckinC                                         BECAUSE R(I+1) IS NOT AVAILABL
ckin       BDYI    = YR*YR+YI*YI
ckin       BDYR    = (BINV(NMAXM)*YR+AUX2(NMAXM)*YI)/BDYI
ckin       BDYI    = (AUX2(NMAXM)*YR-BINV(NMAXM)*YI)/BDYI
ckin       RREAL(NGRID) = EZW*RIM1+FSX*RREAL(NGRID)
ckin     &               - (XR*BDYR-XI*BDYI)
ckin       RIMAG(NGRID) = -XR*BDYI-XI*BDYR
ckin       UR      = POTEFF(NGRID)
ckin       UI      = AUX3(NGRID)
ckin       ZR      = -OFFDIG+EZW*UR
ckin       ZI      = EZW*UI
ckin       BINV(NGRID) = ZR
ckin       AUX2(NGRID) = ZI
ckin       YR      = OFFDG2+FSX*UR-(ZR*BDYR-ZI*BDYI)
ckin       YI      = FSX*UI-(ZR*BDYI+ZI*BDYR)
ckin       AINV(NGRID) = YR
ckin       AUX1(NGRID) = YI
ckinC
ckinC      IN THE THIRD PASS RESOLVE THE REMAINING UPPER TRIANGULAR EQUATION
ckinC      BACKWARD.
ckinC
ckin       DENOM  = AINV(NGRID)*AINV(NGRID)+AUX1(NGRID)*AUX1(NGRID)
ckin       XR     = (RREAL(NGRID)*AINV(NGRID)+RIMAG(NGRID)*AUX1(NGRID))
ckin     &         /DENOM
ckin       XI     = (RIMAG(NGRID)*AINV(NGRID)-RREAL(NGRID)*AUX1(NGRID))
ckin     &         /DENOM
ckin       RREAL(NGRID) = XR*SQM0
ckin       RIMAG(NGRID) = XI*SQM0
ckin       DO 15 I=NMAXM,2,(-1)
ckin       YR     = RREAL(I) -(XR*BINV(I+1)-XI*AUX2(I+1))
ckin       YI     = RIMAG(I) -(XI*BINV(I+1)+XR*AUX2(I+1))
ckin       DENOM  = AINV(I)*AINV(I)+AUX1(I)*AUX1(I)
ckin       XR     = (YR*AINV(I)+YI*AUX1(I))/DENOM
ckin       XI     = (YI*AINV(I)-YR*AUX1(I))/DENOM
ckin       RREAL(I)   = XR*SQM0
ckin       RIMAG(I)   = XI*SQM0
ckin   15  CONTINUE
ckin       RREAL(1)   = ZERO
ckin       RIMAG(1)   = ZERO
ckinC
ckin       RETURN
ckin       END
C-----INVCMP------------ PART OF THE QHQ PACKAGE -------------------------------
C
      SUBROUTINE INVCMP(RREAL,RIMAG,RHS,EINVR,EINVI,Langmo)
C
      INCLUDE 'PARAMS.INC'
C
      DIMENSION RREAL(KG1),RIMAG(KG1),RHS(KG1)
C
      INCLUDE 'WAVFUN.INC'
      INCLUDE 'POTINV.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'CONSTN.INC'
C
       DATA FSX,EZW/.83333333333333,.8333333333333E-1/
cc       DATA FSX,EZW/1.0,0.0/
C
C
C      DETERMINES 'R' AS SOLUTION OF THE INHOMOGENOUS LINEAR DIFF. EQ.:
C       " (-D/DX (SQB(X)**2) D/DX + (POTL(X) - EINV ) ) R = RHS  "
C      FOR GIVEN 'RHS' WHERE 'SQB**2' IS X-DEPENDENT INVERSE MASS,
C      'POTL' A POTENTIAL AND 'EINV' A SHIFT CONSTANT.
C      THE ENERGY-SHIFT 'EINV' IS COMPLEX WITH REAL PART 'EINVR' AND
C      IMAGINARY PART 'EINVI'.
C      ACCORDINGLY THE RESPONSE FIELD 'R' IS COMPLEX WITH REAL PART 'RREAL'
C      AND IMAGINARY PART 'RIMAG'.
C
C      THE ARRAYS IN THE CALLING LIST MEAN:
C       RREAL  = REAL PART OF THE THE SOLUTION OF THE EQUATION
C       RIMAG  = IMAGINARY PART OF THE THE SOLUTION OF THE EQUATION
C       RHS    = THE RIGHT HAND SIDE (INHOMOGENEITY) OF THE EQUATION.
C                R AND RHS MAY BE THE SAME ARRAYS IN THE CALLING ROUTINE
C                RHS IS THEN OVERWRITTEN.
C      THE PARAMETERS IN THE CALLING LIST MEAN:
C       EINVR  = REAL PART OF THE ENERGY-SHIFT
C       EINVI  = IMAGINARY PART OF THE ENERGY-SHIFT
C
C      INFORMATION ON THE HAMILTONIAN IS TRANSFERRED VIA COMMON /POTINV/:
C       SQB    = INVERSE SQUARE-ROOT OF X-DEPENDENT EFFECTIVE MASS
C       POTL   = EFFECTIVE POTENTIAL BEFORE MASS TRANSFORMATION,
C                TO BE PREPARED IN SUBROUTINE 'WELFU'.
C
C      INFORMATION ON THE GRID IS TRANSFERRED VIA COMMON /GRID/.
C
C      NOTE: THE PROGRAM IS WRITTEN FOR SOLVING THE SCHROEDINGER EQUATIO
C            IN RADIALLY SYMMETRIC PROBLEMS. THUS IT IMPLIES THE RHS AND
C            POTEFF TO BE SUCH THAT R(X) VANISHES AT THE ORIGIN.
C
C      THE SOLUTION PROCEEDS AS FOLLOWS:
C
C      BY INTRODUCING SCALED WAVE FUNCTIONS
C       " REFF = R/SQB "
C      AND
C       " RHSEFF = RHS/SQB ".
C      THE EQUATIONS IS TRANSFORMED TO A LAPLACIAN EQUATION
C       " (-D/DX**2 + POTEFF(X) ) REFF = RHSEFF "
C      WITH CONSTANT MASS AND AN EFFECTIVE POTENTIAL
C       " POTEFF = (POTL(X)-UMIN)/(SQB(X)*SQB(X)) + D2B(X) "
C      WHERE
C       " D2B = (D/DX**2 SQB)/SQB " .
C      AFTER SOLVING THE LAPLACIAN EQUATION, THE SOLUTION IS TRANSFORMED
C      BACK BY
C       " R   = REFF*SQB" .
C
C      THE D/DX**2 IS TREATED EFFECTIVELY IN 5-POINT PRECISION.
C      IN ADDITION TO A STANDARD 3-POINT FORMULA FOR D/DX**2
C      A COMPLEMENTING 3-POINT AVERAGING OF UEFF*REFF AND RHSEFF
C      IS USED WHICH BOTH TOGETHER MAKE UP THE 5-POINT PRECISION.
C      (SEE R.ZURMUEHL, "PRAKTISCHE MATHEMATIK", CH. VII, SECT. 3,
C       SPRINGER VERLAG, BERLIN 1963)
C
C
C      THE WORKSPACES FOR THE HAMILTONIAN, 'POTEFF', 'AINV' AND 'BINV' ARE
C      NOW BECOMING COMPLEX. WE AVOID EXPLICIT COMPLEX ARITHMETIC
C      AND ASSIGN THE WORKSPACES 'AUX..' AS FOLLOWS:
C         AUX1   = IMAGINARY PART OF 'AINV'
C         AUX2   = IMAGINARY PART OF 'BINV'
C         AUX3   = IMAGINARY PART OF 'POTEFF'
C
C-----------------------------------------------------------------------
C
C      SQUARE-ROOT OF INVERSE EFF. MASS ON 'SQB',
C      SECOND DERIVATIVE OF 'SQB' ON 'D2B'.
C      FINALLY ASSEMBLING  THE EFFECTIVE POTENTIAL ON 'POTEFF'.
C
       SQM(1) = ONE/SQB(1)
       DO 80 I=2,NGRID
       SQM(I)    = ONE/SQB(I)
       x         = SQM(I)*SQM(I)
       POTEFF(I) = D2B(I)+(POTL(I)-EINVR)*x
       AUX3(I)   = (-EINVI)*x
   80  CONTINUE
C
C      TRANSFORM TO EFFECTIVE R.H.S., STORE RESULT ON R.
C
       DO 17 I=1,NGRID
   17  RREAL(I)   = RHS(I)*SQM(I)
C
C     PREPARE PARAMETERS FOR KINETIC ENERGY AND LOOP-ENDS.
C
       OFFDIG = ONE/(RSTEP*RSTEP)
       OFFDG2 = OFFDIG+OFFDIG
       NMAXM  = NGRID-1
       if(langmo.eq.1 .and. fsx.lt.one) then
cc        xll    = (langmo*langmo+langmo)*offdig/144.0
        xll    = (langmo*langmo+langmo)*offdig/12.0
       else
        xll    = zero
       endif
C
C      BEGIN FIRST PASS ELIMINATING THE LOWER OFF-DIAGONAL ELEMENT.
C      NEW DIAGONAL AND UPPER OFF-DIAGONAL ELEMNTS ARE SAVED ON A AND B.
C      THE SECOND PASS, FORWARD ELIMINATION OF THE R.H.S., IS DONE WITHI
C      THE SAME LOOP; RESULT IS STORED ON FIELD T.
C
C
       BINV(1) = ZERO
       AUX1(1) = ZERO
       AINV(1) = ONE
       AUX2(1) = ZERO
       XR      = ZERO
       XI      = ZERO
       YR      = ONE
       YI      = ZERO
       RIM1    = RREAL(1)
C
ctest
ctest       write(6,'(a)')'  r        ar          ai          br          bi'
ctest
       DO 13 I=2,NMAXM
C                                         FORWARD ELIM. FOR R.H.S.
       BDYI    = YR*YR+YI*YI
       BDYR    = (BINV(I-1)*YR+AUX2(I-1)*YI)/BDYI
       BDYI    = (AUX2(I-1)*YR-BINV(I-1)*YI)/BDYI
       X       = EZW*(RREAL(I+1)+RIM1)+FSX*RREAL(I)
     &          - (XR*BDYR-XI*BDYI)
       XI      = -XR*BDYI-XI*BDYR
       XR      = X
       RIM1    = RREAL(I)
       RREAL(I)= XR
       RIMAG(I)= XI
C                                         FORWARD ELIM. FOR MATRIX
       UR      = POTEFF(I)
       UI      = AUX3(I)
       ZR      = -OFFDIG+EZW*UR
       ZI      = EZW*UI
c       if(i.eq.3)  then
c        zr = -offdig
c        zi = zero
c       endif
cc       if(i.eq.3)  zr = zr-xll
       BINV(I) = ZR
       AUX2(I) = ZI
       YR      = OFFDG2+FSX*UR-(ZR*BDYR-ZI*BDYI)
       YI      = FSX*UI-(ZR*BDYI+ZI*BDYR)
cc       if(i.eq.3)  zr = zr+xll
cc       if(i.eq.2) yr = yr+xll*16.0
       if(i.eq.2) yr = yr+xll
c       if(i.eq.2) then
c        yr = yr+(ezw+ezw)*ur
c        yi = yi+(ezw+ezw)*ui
c       endif
       AINV(I) = YR
       AUX1(I) = YI
ctest
ctest       write(6,'(1x,f4.2,4g12.4)') (I-1)*rstep,yr,yi,zr,zi
ctest
   13  CONTINUE
C                                         LAST ELIM. STEP EXPLICITELY
C                                         BECAUSE R(I+1) IS NOT AVAILABL
       BDYI    = YR*YR+YI*YI
       BDYR    = (BINV(NMAXM)*YR+AUX2(NMAXM)*YI)/BDYI
       BDYI    = (AUX2(NMAXM)*YR-BINV(NMAXM)*YI)/BDYI
       RREAL(NGRID) = EZW*RIM1+FSX*RREAL(NGRID)
     &               - (XR*BDYR-XI*BDYI)
       RIMAG(NGRID) = -XR*BDYI-XI*BDYR
       UR      = POTEFF(NGRID)
       UI      = AUX3(NGRID)
       ZR      = -OFFDIG+EZW*UR
       ZI      = EZW*UI
       BINV(NGRID) = ZR
       AUX2(NGRID) = ZI
       YR      = OFFDG2+FSX*UR-(ZR*BDYR-ZI*BDYI)
       YI      = FSX*UI-(ZR*BDYI+ZI*BDYR)
       AINV(NGRID) = YR
       AUX1(NGRID) = YI
c
c      determinant for testing
c
       detr   = one
       deti   = zero
       do 9010 ir=1,ngrid
       x      = detr*ainv(ir)-deti*aux1(ir)
       deti   = detr*aux1(ir)+deti*ainv(ir)
       detr   = x
 9010  continue
ctest       write(6,'(a,2g12.4))') ' Determ.=',detr,deti
C
C      IN THE THIRD PASS RESOLVE THE REMAINING UPPER TRIANGULAR EQUATION
C      BACKWARD.
C
       DENOM  = AINV(NGRID)*AINV(NGRID)+AUX1(NGRID)*AUX1(NGRID)
       XR     = (RREAL(NGRID)*AINV(NGRID)+RIMAG(NGRID)*AUX1(NGRID))
     &         /DENOM
       XI     = (RIMAG(NGRID)*AINV(NGRID)-RREAL(NGRID)*AUX1(NGRID))
     &         /DENOM
       RREAL(NGRID) = XR*SQM(NGRID)
       RIMAG(NGRID) = XI*SQM(NGRID)
       DO 15 I=NMAXM,2,(-1)
       YR     = RREAL(I) -(XR*BINV(I+1)-XI*AUX2(I+1))
       YI     = RIMAG(I) -(XI*BINV(I+1)+XR*AUX2(I+1))
       DENOM  = AINV(I)*AINV(I)+AUX1(I)*AUX1(I)
       XR     = (YR*AINV(I)+YI*AUX1(I))/DENOM
       XI     = (YI*AINV(I)-YR*AUX1(I))/DENOM
       RREAL(I)   = XR*SQM(I)
       RIMAG(I)   = XI*SQM(I)
   15  CONTINUE
       RREAL(1)   = ZERO
       RIMAG(1)   = ZERO
C
       RETURN
       END
cpivC-----INVCMP------------ PART OF THE QHQ PACKAGE -------------------------------
cpivC
cpiv      SUBROUTINE INVCMP(RREAL,RIMAG,RHS,EINVR,EINVI)
cpivC
cpivC     This is the version of INVCMP which does pivoting.
cpivC     It is, of course, a bit more involved in coding, and
cpivC     sometimes in running.
cpivC     Switch to that variant if the simple version of INVCMP
cpivC     runs in obvious trouble.
cpivC
cpiv      INCLUDE 'PARAMS.INC'
cpiv      parameter (koffsv=20)
cpivC
cpiv      logical tchang(kg1)
cpiv      dimension cr(koffsv),ci(koffsv)
cpivc
cpiv      DIMENSION RREAL(KG1),RIMAG(KG1),RHS(KG1)
cpivC
cpiv      INCLUDE 'WAVFUN.INC'
cpiv      INCLUDE 'POTINV.INC'
cpiv      INCLUDE 'GRID.INC'
cpiv      INCLUDE 'CONSTN.INC'
cpivC
cpivcc       DATA FSX,EZW/.83333333333333,.8333333333333E-1/
cpiv       DATA FSX,EZW/1.0,0.0/
cpivC
cpivC
cpivC      DETERMINES 'R' AS SOLUTION OF THE INHOMOGENOUS LINEAR DIFF. EQ.:
cpivC       " (-D/DX (SQB(X)**2) D/DX + (POTL(X) - EINV ) ) R = RHS  "
cpivC      FOR GIVEN 'RHS' WHERE 'SQB**2' IS X-DEPENDENT INVERSE MASS,
cpivC      'POTL' A POTENTIAL AND 'EINV' A SHIFT CONSTANT.
cpivC      THE ENERGY-SHIFT 'EINV' IS COMPLEX WITH REAL PART 'EINVR' AND
cpivC      IMAGINARY PART 'EINVI'.
cpivC      ACCORDINGLY THE RESPONSE FIELD 'R' IS COMPLEX WITH REAL PART 'RREAL'
cpivC      AND IMAGINARY PART 'RIMAG'.
cpivC
cpivC      THE ARRAYS IN THE CALLING LIST MEAN:
cpivC       RREAL  = REAL PART OF THE THE SOLUTION OF THE EQUATION
cpivC       RIMAG  = IMAGINARY PART OF THE THE SOLUTION OF THE EQUATION
cpivC       RHS    = THE RIGHT HAND SIDE (INHOMOGENEITY) OF THE EQUATION.
cpivC                R AND RHS MAY BE THE SAME ARRAYS IN THE CALLING ROUTINE
cpivC                RHS IS THEN OVERWRITTEN.
cpivC      THE PARAMETERS IN THE CALLING LIST MEAN:
cpivC       EINVR  = REAL PART OF THE ENERGY-SHIFT
cpivC       EINVI  = IMAGINARY PART OF THE ENERGY-SHIFT
cpivC
cpivC      INFORMATION ON THE HAMILTONIAN IS TRANSFERRED VIA COMMON /POTINV/:
cpivC       SQB    = INVERSE SQUARE-ROOT OF X-DEPENDENT EFFECTIVE MASS
cpivC       POTL   = EFFECTIVE POTENTIAL BEFORE MASS TRANSFORMATION,
cpivC                TO BE PREPARED IN SUBROUTINE 'WELFU'.
cpivC
cpivC      INFORMATION ON THE GRID IS TRANSFERRED VIA COMMON /GRID/.
cpivC
cpivC      NOTE: THE PROGRAM IS WRITTEN FOR SOLVING THE SCHROEDINGER EQUATIO
cpivC            IN RADIALLY SYMMETRIC PROBLEMS. THUS IT IMPLIES THE RHS AND
cpivC            POTEFF TO BE SUCH THAT R(X) VANISHES AT THE ORIGIN.
cpivC
cpivC      THE SOLUTION PROCEEDS AS FOLLOWS:
cpivC
cpivC      BY INTRODUCING SCALED WAVE FUNCTIONS
cpivC       " REFF = R/SQB "
cpivC      AND
cpivC       " RHSEFF = RHS/SQB ".
cpivC      THE EQUATIONS IS TRANSFORMED TO A LAPLACIAN EQUATION
cpivC       " (-D/DX**2 + POTEFF(X) ) REFF = RHSEFF "
cpivC      WITH CONSTANT MASS AND AN EFFECTIVE POTENTIAL
cpivC       " POTEFF = (POTL(X)-UMIN)/(SQB(X)*SQB(X)) + D2B(X) "
cpivC      WHERE
cpivC       " D2B = (D/DX**2 SQB)/SQB " .
cpivC      AFTER SOLVING THE LAPLACIAN EQUATION, THE SOLUTION IS TRANSFORMED
cpivC      BACK BY
cpivC       " R   = REFF*SQB" .
cpivC
cpivC      THE D/DX**2 IS TREATED EFFECTIVELY IN 5-POINT PRECISION.
cpivC      IN ADDITION TO A STANDARD 3-POINT FORMULA FOR D/DX**2
cpivC      A COMPLEMENTING 3-POINT AVERAGING OF UEFF*REFF AND RHSEFF
cpivC      IS USED WHICH BOTH TOGETHER MAKE UP THE 5-POINT PRECISION.
cpivC      (SEE R.ZURMUEHL, "PRAKTISCHE MATHEMATIK", CH. VII, SECT. 3,
cpivC       SPRINGER VERLAG, BERLIN 1963)
cpivC
cpivC
cpivC      THE WORKSPACES FOR THE HAMILTONIAN, 'POTEFF', 'AINV' AND 'BINV' ARE
cpivC      NOW BECOMING COMPLEX. WE AVOID EXPLICIT COMPLEX ARITHMETIC
cpivC      AND ASSIGN THE WORKSPACES 'AUX..' AS FOLLOWS:
cpivC         AUX1   = IMAGINARY PART OF 'AINV'
cpivC         AUX2   = IMAGINARY PART OF 'BINV'
cpivC         AUX3   = IMAGINARY PART OF 'POTEFF'
cpivC
cpivC-----------------------------------------------------------------------
cpivC
cpivC      SQUARE-ROOT OF INVERSE EFF. MASS ON 'SQB',
cpivC      SECOND DERIVATIVE OF 'SQB' ON 'D2B'.
cpivC      FINALLY ASSEMBLING  THE EFFECTIVE POTENTIAL ON 'POTEFF'.
cpivC
cpiv       SQM(1) = ONE/SQB(1)
cpiv       DO 80 I=2,NGRID
cpiv       SQM(I)    = ONE/SQB(I)
cpiv       POTEFF(I) = D2B(I)+(POTL(I)-EINVR)*(SQM(I)*SQM(I))
cpiv       AUX3(I)   = (-EINVI)*(SQM(I)*SQM(I))
cpiv   80  CONTINUE
cpivC
cpivC      TRANSFORM TO EFFECTIVE R.H.S., STORE RESULT ON R.
cpivC
cpiv       DO 17 I=1,NGRID
cpiv   17  RREAL(I)   = RHS(I)*SQM(I)
cpivC
cpivC     PREPARE PARAMETERS FOR KINETIC ENERGY AND LOOP-ENDS.
cpivC
cpiv       OFFDIG = ONE/(RSTEP*RSTEP)
cpiv       OFFDG2 = OFFDIG+OFFDIG
cpiv       NMAXM  = NGRID-1
cpivc
cpivc      setup of tridiagonal elements of the l.h.s. matrix and of the
cpivc      mixing for the r.h.s.
cpivc
cpiv       rim1   = rreal(1)
cpiv       do 50 i=2,nmaxm
cpiv       rim1   = ezw*(rreal(i+1)+rim1)+fsx*rreal(i)
cpiv       rreal(i) = rim1
cpiv       rimag(i) = zero
cpiv   50  continue
cpiv       rreal(ngrid) = ezw*rim1+fsx*rreal(ngrid)
cpiv       rimag(ngrid) = zero
cpiv       BINV(1) = ZERO
cpiv       AUX1(1) = ZERO
cpiv       AINV(1) = ONE
cpiv       AUX2(1) = ZERO
cpiv       do 60 i=2,ngrid
cpiv       ur       = poteff(i)
cpiv       ainv(i)  = offdg2+fsx*ur
cpiv       binv(i)  = -offdig+ezw*ur
cpiv       ui       = aux3(i)
cpiv       aux1(i)  = fsx*ui
cpiv       aux2(i)  = ezw*ui
cpiv   60  continue
cpivC
cpivC      BEGIN FIRST PASS ELIMINATING THE LOWER OFF-DIAGONAL ELEMENT.
cpivC      NEW DIAGONAL AND UPPER OFF-DIAGONAL ELEMNTS ARE SAVED ON A AND B.
cpivC      THE SECOND PASS, FORWARD ELIMINATION OF THE R.H.S., IS DONE WITHI
cpivC      THE SAME LOOP; RESULT IS STORED ON FIELD T.
cpivC
cpivC
cpiv       XR      = ZERO
cpiv       XI      = ZERO
cpiv       YR      = ONE
cpiv       YI      = ZERO
cpiv       zr      = zero
cpiv       zi      = zero
cpiv       numoff  = 1
cpivC
cpivctest
cpivctest       write(6,'(a)')'  r        ar          ai          br          bi'
cpivctest
cpiv       DO 13 I=2,ngrid
cpivctest
cpivctest       write(6,'(1x,f4.2,6g12.4)') (I-2)*rstep,yr,yi,zr,zi,binv(i-1),
cpivctest     &                             aux2(i-1)
cpivctest
cpivc                                         search pivot element
cpivc       if(abs(zr).gt.abs(yr).and.I.lt.ngrid) then
cpiv       if(zr*zr+zi*zi.gt.yr*yr+yi*yi.and.I.lt.ngrid) then
cpivc                                         exchange row I and I-1
cpiv        ainv(i-1) = zr
cpiv        aux1(i-1) = zi
cpiv        x       = ainv(i)
cpiv        ainv(i) = binv(i)
cpiv        binv(i) = x
cpiv        x       = aux1(i)
cpiv        aux1(i) = aux2(i)
cpiv        aux2(i) = x
cpiv        x       = rreal(i)
cpiv        rreal(i) = rreal(i-1)
cpiv        rreal(i-1) = x
cpiv        rimag(i) = rimag(i-1)
cpiv        rimag(i-1) = zero
cpivC                                         PREPARE B/A
cpiv        BDYI    = zr*zr+zi*zi
cpiv        BDYR    = (zr*YR+zi*YI)/BDYI
cpiv        BDYI    = (-zi*YR+zr*YI)/BDYI
cpivc                                         set counters and extra element
cpiv        if(numoff.gt.koffsv) stop ' too much pivoting in INVCMP'
cpiv        cr(numoff) = binv(i+1)
cpiv        ci(numoff) = aux2(i+1)
cpiv        tchang(i)  = .true.
cpiv        numoff     = 1+numoff
cpiv       else
cpivC                                         PREPARE B/A
cpiv        BDYI    = YR*YR+YI*YI
cpiv        BDYR    = (zr*YR+zi*YI)/BDYI
cpiv        BDYI    = (zi*YR-zr*YI)/BDYI
cpiv        tchang(i)  = .false.
cpiv       endif
cpivC                                         FORWARD ELIM. FOR MATRIX
cpiv       ZR      = BINV(I)
cpiv       ZI      = AUX2(I)
cpiv       YR      = ainv(i) - (ZR*BDYR-ZI*BDYI)
cpiv       YI      = aux1(i) - (ZR*BDYI+ZI*BDYR)
cpiv       AINV(I) = YR
cpiv       AUX1(I) = YI
cpiv       if(tchang(i)) then
cpivc                                         modified off diogonal piece
cpiv        zr        = binv(i+1)
cpiv        zi        = aux2(i+1)
cpiv        binv(i+1) = -(zr*bdyr-zi*bdyi)
cpiv        aux2(i+1) = -(zr*bdyi+zi*bdyr)
cpivC                                         FORWARD ELIM. FOR R.H.S.
cpiv        XR      = xr - (rreal(i-1)*BDYR)
cpiv        XI      = xi - rreal(i-1)*BDYI
cpiv        RREAL(I)= XR
cpiv        RIMAG(I)= XI
cpiv       else
cpivC                                         FORWARD ELIM. FOR R.H.S.
cpiv        X       = RREAL(I) - (XR*BDYR-XI*BDYI)
cpiv        XI      = rimag(i) - (XR*BDYI+XI*BDYR)
cpiv        xr      = x
cpiv        RREAL(I)= XR
cpiv        RIMAG(I)= XI
cpiv       endif
cpivc                                         prepare B(I-1) for next step
cpiv       zr = -offdig+ezw*poteff(i)
cpiv       zi = ezw*aux3(i)
cpiv   13  CONTINUE
cpivc
cpivc      determinant for testing
cpivc
cpiv       detr   = one
cpiv       deti   = zero
cpiv       do 9010 ir=1,ngrid
cpiv       x      = detr*ainv(ir)-deti*aux1(ir)
cpiv       deti   = detr*aux1(ir)+deti*ainv(ir)
cpiv       detr   = x
cpiv 9010  continue
cpivctest       write(6,'(a,2g12.4))') ' Determ.=',detr,deti
cpivC
cpivC      IN THE THIRD PASS RESOLVE THE REMAINING UPPER TRIANGULAR EQUATION
cpivC      BACKWARD.
cpivC
cpiv       DENOM  = AINV(NGRID)*AINV(NGRID)+AUX1(NGRID)*AUX1(NGRID)
cpiv       XR     = (RREAL(NGRID)*AINV(NGRID)+RIMAG(NGRID)*AUX1(NGRID))
cpiv     &         /DENOM
cpiv       XI     = (RIMAG(NGRID)*AINV(NGRID)-RREAL(NGRID)*AUX1(NGRID))
cpiv     &         /DENOM
cpiv       RREAL(NGRID) = XR*SQM(NGRID)
cpiv       RIMAG(NGRID) = XI*SQM(NGRID)
cpiv       DO 15 I=NMAXM,2,(-1)
cpiv       YR     = RREAL(I) -(XR*BINV(I+1)-XI*AUX2(I+1))
cpiv       YI     = RIMAG(I) -(XI*BINV(I+1)+XR*AUX2(I+1))
cpiv       if(tchang(i+1)) then
cpiv        xr = rreal(i+2)/sqm(i+2)
cpiv        xi = rimag(i+2)/sqm(i+2)
cpiv        numoff = numoff-1
cpiv        if(numoff.lt.1) stop ' second off-diagonals exhausted'
cpiv        yr = yr - (xr*cr(numoff)-xi*ci(numoff))
cpiv        yi = yi - (xr*ci(numoff)+xi*cr(numoff))
cpiv       endif
cpiv       DENOM  = AINV(I)*AINV(I)+AUX1(I)*AUX1(I)
cpiv       XR     = (YR*AINV(I)+YI*AUX1(I))/DENOM
cpiv       XI     = (YI*AINV(I)-YR*AUX1(I))/DENOM
cpiv       RREAL(I)   = XR*SQM(I)
cpiv       RIMAG(I)   = XI*SQM(I)
cpiv   15  CONTINUE
cpiv       RREAL(1)   = ZERO
cpiv       RIMAG(1)   = ZERO
cpivC
cpiv       RETURN
cpiv       END
C-----ADDWF2------------------------------------------------------
C
      SUBROUTINE ADDWF2(QWF1,QWF2,FAC1,FAC2)
C
      INCLUDE 'PARAMS.INC'
C
      INCLUDE 'WAVFUN.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'RESULT.INC'
      INCLUDE 'LINRES.INC'
      INCLUDE 'CONSTN.INC'
C
      DIMENSION QWF1(KX11,KLRPA),QWF2(KX11,KLRPA)
C
C------------------------------------------------------------------------------
C
      INCLUDE 'INITDO.INC'
C
      DO 45 IR=1,NGRID
   45 QWF1(NADD+IR,JDIF) = FAC2*QWF2(NADD+IR,JDIF)
     &                    +FAC1*QWF1(NADD+IR,JDIF)
C
      INCLUDE 'ENDDO.INC'
C
      RETURN
      END
C-----ADDWF3------------------------------------------------------
C
      SUBROUTINE ADDWF3(QWF1,QWF2,QWF3,FAC1,FAC2,FAC3)
C
      INCLUDE 'PARAMS.INC'
C
      INCLUDE 'WAVFUN.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'RESULT.INC'
      INCLUDE 'LINRES.INC'
      INCLUDE 'CONSTN.INC'
C
      DIMENSION QWF1(KX11,KLRPA),QWF2(KX11,KLRPA),QWF3(KX11,KLRPA)
C
C------------------------------------------------------------------------------
C
      INCLUDE 'INITDO.INC'
C
      DO 45 IR=1,NGRID
   45 QWF1(NADD+IR,JDIF) = FAC2*QWF2(NADD+IR,JDIF)
     &                    +FAC3*QWF3(NADD+IR,JDIF)
     &                    +FAC1*QWF1(NADD+IR,JDIF)
C
      INCLUDE 'ENDDO.INC'
C
      RETURN
      END
C-----HPSI3 ------------ PART OF THE HARTREE-FOCK PACKAGE ----------------------
C
      SUBROUTINE HPSI3(WFIN,HWFOUT,D1WFIN,D2WFIN,IFD1WF,JSA,ISO,SUBFAC)
C
      INCLUDE 'PARAMS.INC'
C
C     ACTION OF THE MEAN-FIELD HAMILTONIAN
C
C                     ^
C         |HWFOUT> = (H - SUBFAC) |WFIN>
C
C     WHERE:
C       WFIN    = INPUT WAVE FUNCTION FIELD ,
C       HWFOUT  = OUTPUT WAVE FUNCTION FIELD ,
C       D1WFIN  = FIRST DERIVATIVE IF 'WFIN'.
C       D2WFIN  = SECOND DERIVATIVE IF 'WFIN'.
C       IFD1WF  = DECIDES ON THE COMPUTATION OF 'D1WFIN':
C                 IFD1WF>0 --> THE 'D1WFIN' IS INTERPRETED AS
C                              WORKSPACE AND THE FIRST DERIVATIVE
C                              IS COMPUTED ANEW;
C                 IFD1WF=0 --> THE 'D1WFIN' CONTAINS ALREADY THE
C                              FIRST DERIVATIVE IF 'WFIN'.
C       JSA     = COMBINED ANGULAR MOMENTUM QUANTUM NUMBER
C                 OF WFIN AND HWFOUT,
C       ISO     = ISOSPIN OF THE WAVE FUNCTIONS,
C       SUBFAC  = CONSTANT ENERGY SUBTRACTED FROM HAMILTONIAN.
C
C     THE FIELDS FOR THE MEAN-FIELD HAMILTONIAN ARE TRANSFERRED
C     VIA COMMON /WAVFUN/:
C       POT     = POTENTIAL
C       POTLS   = SPIN-ORBIT POTENTIAL
C       BMASS   = INVERSE EFFECTIVE MASS
C       D1B     = DERIVATIVE OF INVERSE EFFECTIVE MASS
C
C     GRID PROPERTIES ARE TRANSFERRED VIA COMMON /GRID/
C
      DIMENSION WFIN(KG1),HWFOUT(KG1),D1WFIN(KG1),D2WFIN(KG1)
C
      INCLUDE 'WAVFUN.INC'
      INCLUDE 'GRID.INC'
      INCLUDE 'CONSTN.INC'
C
C-------------------------------------------------------------------------------
C
C     THE QUANTUM NUMBERS OF THE WAVE FUNCTION
C      LA      = ORBITAL ANGULAR MOMENTUM OF STATE  WFIN>
C      PARITY  = PARITY FOR THE RADIAL WAVE FUNCTION (NOTE THAT IT
C                CONTAINS AN EXTRA FACTOR 'R' SUCH THAT THE 1S1/2
C                STATE HAS PARITY=-1 ETC)
C      ISOADD  = STACK-POINTER FOR POTENTIAL FIELDS POT, POTLS
C                AND ISOADD
C      XL0     = SPIN-ORBIT SPLITTING FOR STATE  WFIN>
C
      LA     = JSA/2
      XLP    = (LA*LA+LA)
      PARITY = (MOD(LA,2)*2-1)
      ISOADD = (ISO-1)*NGRID
      XL0    = WLS(JSA)
C
C     THE DERIVATIVES OF THE WAVE FUNCTION
C
      IF(IFD1WF.GT.0)
     &   CALL DERV15(WFIN(1),D1WFIN(1),PARITY,NGRID,RSTEP)
      CALL DERV23(WFIN(1),D2WFIN(1),PARITY,NGRID,RSTEP)
C
C     THE ACTION OF THE HAMILTONIAN
C
      HWFOUT(1) = ZERO
      DO 10 IR=2,NGRID
      IRS    = IR+ISOADD
      POTACT = POT(IRS)+XL0*POTLS(IRS)
     &        +BMASS(IRS)*XLP*ZENTRF(IR)+D1B(IRS)*RINVS(IR)
      HWFOUT(IR) = (POTACT-SUBFAC)*WFIN(IR)
     &            -D1B(IRS)*D1WFIN(IR)-BMASS(IRS)*D2WFIN(IR)
   10 CONTINUE
C
      RETURN
      END
