C-----RECOUP ----------------------------------------------------------
C
       PROGRAM RECOUP
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   !
C
C   RECOUPLES SKYRME PARAMETERS FROM t0,.... to B0,....
C   AND VICE VERSA.
C   INPUT AND OUTPUT IN STANDARD FORM FOR FIT-PROGRAM.
C
C     IFTINP = SWITCH FOR INPUT OF SKYRME PARAMETERS
C               1 --> STANDARD FORM T0,T1,...
C               0 --> RECOUPLED FORM B0,B0P,...
C
C***********   SAMPLE INPUT FILE, UNIT 01  **************************
C
C 011,'  WALECKA-MODEL   ' ,ISTEUR,COMMEN
C 03,        'NWDR'
C 2,1,3,0,0,0,0,0,0,0,0,0,0,0,0,0,        'AUSWAHL DER MESSGROESSEN'
C 2.0,1.0,5.0,0.0,0.0,0.,0.,0.,0.,0.,0.,0.,0.,0., 'FEHLER AN DATEN'
C 00.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,  ' % ODER ABSOLUT
C 00.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., ' KORRREKTUR '
C 010,220,        'IAUNT,IAOB,
C 0,0,0,1,1,1,0,  'NMINV,IF:EX,TM,ENRG,RSO,ZPE,TOND'
C 1.0,,,,         'ACOUL,ADDNEW,DABE,DABR '
C 0,0,1,,         'IWRTE,IST1,IFCURV,ITMAX'
C 15,4,0,1.,,     'NPARIN,NUMAX,IVAR,FREEP,ENDCON'
C -1128.75,1.000,1, T0
C +395.   ,1.000,2, T1
C -095.   ,0.200,3, T2
C +14000. ,0.100,4, T3
C +120.   ,0.010,5, T4
C +0.45   ,0.001,6, X0
C +0.0    ,1.000,7, X1
C +0.0    ,1.000,8, X2
C +1.0    ,1.000,9, X3
C +1.0    ,1.000,10, POWER
C +0.0    ,1.000,11, TYUK
C +0.0    ,1.000,12, XYUK
C +0.5    ,1.000,13, AYUK
C +000.000,1.000,14,  N.N.
C +000.000,1.000,15,  N.N.
C -11,'  STOP THE FIT HERE '
C
      DIMENSION A(15),PAR(15),DPARMO(15),
     &          JP(15),LIFIT(12),IFIT1(12),IPARVA(15)
      DIMENSION SIGMAA(15),SIGMA0(15),DPARM0(15)
C
      CHARACTER*6 ZUORD(15,0:2)
      CHARACTER*80 LINE
C
      DATA ZUORD/'    B0','    B1','    B2','    B3','    T4',
     &           '   B0P','   B1P','   B2P','   B3P',' POWER',
     &           '  BYUK',' BYUKP','  AYUK',' N.N. ',' N.N. ',
     &           '    T0','    T1','    T2','    T3','    T4',
     *           '    X0','    X1','    X2','    X3',' POWER',
     *           '  TYUK','  XYUK','  AYUK',' N.N. ',' N.N. ',
     &           '    C0','    B1','    B2','    C3','    T4',
     &           '   C0P','   B1P','   B2P','   C3P',' POWER',
     &           '  BYUK',' BYUKP','  AYUK',' N.N. ',' N.N. '/
      DATA A/15*0.01/
      DATA RHONM/0.159/
C
C--------------------------------------------------------------------
C
       OPEN(UNIT=1,FORM='FORMATTED',STATUS='OLD',FILE='RECOUP.INP')
       OPEN(UNIT=2,FORM='FORMATTED',STATUS='UNKNOWN',FILE='RECOUP.OUT')
c       OPEN(UNIT=6,FORM='FORMATTED',STATUS='UNKNOWN',FILE='FOR006')
C
C    DEFINE ACTION FROM USER INPUT
C
      WRITE(6,'(A)')
     & ' GIVE "IFTINP" OF OUTPUT (0,1,2):'
      READ(5,*) IFTINN
C    INITIALISE AUXILIARY FIELDS FOR FIT
C    (LIFIT=12 DEACTIVATES THE FIT-POINTER)
C
      DO 9001 I=1,12
 9001 LIFIT(I)=12
C
C      OUTER LOOP. FOR EACH 'IRUN' A NEW COMPUTATION IS STARTED.
C      LOOP STOPS IF 'ISTEUR<0' IS READ IN.
C
      DO 5 IRUN=1,1000
C
C     PRESET SOME CONTROL PARAMETERS.
C
      IFTINP = 1
      ZERO   = 0.0
C
C**********************************************************************
C
C   FIRST BLOCK: READING CONTROL PARAMETERS, FIT OPTIONS,
C                EXPERIMENTAL DATA, AND FORCE PARAMETERS.
C
C**********************************************************************
C
      DO 110 IR=1,10
      READ(1,'(A)',END=9999) LINE
      WRITE(2,'(A)') LINE
      write(6,'(a)') line
  110 CONTINUE
      READ(1,*) NPARIN,NUMAX,IVAR,FREEP,ENDCON,IFTINP
      WRITE(2,'(1X,3(I2,A),F4.1,A,F6.4,A,I1,A)')
     &       NPARIN,',',NUMAX,',',IVAR,',',FREEP,',',ENDCON,',',
     &       IFTINN,',    NPARIN,NUMAX,IVAR,FREEP,ENDCON,IFTINP'
C
C   EINGABE DER ANFANGSPARAMETER
C
      DO 310 J=1,NPARIN
      READ(1,*) A(J),DPARMO(J),IPARVA(J)
      WRITE(6,'(2G13.6,I6)') A(J),DPARMO(J),IPARVA(J)
      DPARM0(J)=DPARMO(J)
      PAR(J)=A(J)
      IF(J.GT.NUMAX) DPARMO(J)=0.
      DO 310 L=1,15
      IF(IPARVA(J).EQ.L) JP(L)=J
  310 CONTINUE
      IF(PAR(JP(1)).EQ.ZERO .AND. PAR(JP(6)).EQ.ZERO)
     &  WRITE(6,'(A,A)') ' --> VERSION: ZERO-RANGE',
     &     ' TWO-BODY FORCE SWITCHED OFF.'
C
C    STORE SOME INITIAL PARAMETER FOR LATER PROTOCOLL ON FILE 13
C
C**********************************************************************
C
C   SECOND BLOCK: PERFORMING THE FIT BY REPEATED CALLS TO THE
C                FITTING SUBROUTINE 'XD1'.
C
C**********************************************************************
C
C     UNFOLD PARAMETERS
C
       IF(IFTINP.EQ.1) THEN
C                    ASSIGN T0,... Values
        T0     = PAR(JP(1))
        T1     = PAR(JP(2))
        T2     = PAR(JP(3))
        T3     = PAR(JP(4))
        T4     = PAR(JP(5))
        X0     = PAR(JP(6))
        X1     = PAR(JP(7))
        X2     = PAR(JP(8))
        X3     = PAR(JP(9))
        POWER  = PAR(JP(10))
        TYUK   = PAR(JP(11))
        XYUK   = PAR(JP(12))
        AYUK   = PAR(JP(13))
C                   RECOUPLE TO B0,...
        BS0     = T0*(1.0+0.5*X0)
        B0P    = T0*(0.5+X0)
        BS3     = T3*(1.0+0.5*X3)
        B3P    = T3*(0.5+X3)
        A1     = T1*(1.0+0.5*X1)
        A1P    = T1*(0.5+X1)
        A2     = T2*(1.0+0.5*X2)
        A2P    = T2*(0.5+X2)
        BS1     = 0.25*(A1+A2)
        B1P    = 0.25*(A1P-A2P)
        BS2     = 0.125*(3.0*A1-A2)
        B2P    = 0.125*(3.0*A1P+A2P)
        BSYUK   = TYUK*(1.0+0.5*XYUK)
        BYUKP  = TYUK*(0.5+XYUK)
        B0     = BS0-0.5*B0P
        B1     = BS1-0.5*B1P
        B2     = BS2-0.5*B2P
        B3     = BS3-0.5*B3P
        POWFAC = RHONM**POWER/6.0
        C0     = B0+POWFAC*B3
        C3     = B0-POWFAC*B3
        C0P    = B0P+POWFAC*B3P
        C3P    = B0P-POWFAC*B3P
       ELSE
        IF(IFTINP.EQ.0) THEN
         B0     = PAR(JP(1))
         B3     = PAR(JP(4))
         B0P    = PAR(JP(6))
         B3P    = PAR(JP(9))
         POWER  = PAR(JP(10))
         POWFAC = RHONM**POWER/6.0
         C0     = B0+POWFAC*B3
         C3     = B0-POWFAC*B3
         C0P    = B0P+POWFAC*B3P
         C3P    = B0P-POWFAC*B3P
        ELSE
         C0     = PAR(JP(1))
         C3     = PAR(JP(4))
         C0P    = PAR(JP(6))
         C3P    = PAR(JP(9))
         POWER  = PAR(JP(10))
         POWFAC = 3.0/RHONM**POWER
         B0     = 0.5*(C0+C3)
         B3     = POWFAC*(C0-C3)
         B0P    = 0.5*(C0P+C3P)
         B3P    = POWFAC*(C0P-C3P)
        ENDIF
        B1     = PAR(JP(2))
        B2     = PAR(JP(3))
        T4     = PAR(JP(5))
        B1P    = PAR(JP(7))
        B2P    = PAR(JP(8))
        BYUK   = PAR(JP(11))
        BYUKP  = PAR(JP(12))
        AYUK   = PAR(JP(13))
C               Recouple Parameters to T0 etc
        BS0    = B0+0.5*B0P
        BS1    = B1+0.5*B1P
        BS2    = B2+0.5*B2P
        BS3    = B3+0.5*B3P
        fourtr = 4.0/3.0
        T0     = fourtr*(BS0-0.5*B0P)
        X0     = B0P/T0-0.5
        T3     = fourtr*(BS3-0.5*B3P)
        X3     = B3P/T3-0.5
        A1     = BS1+BS2+BS2
        A1P    = B1P+B2P+B2P
        A2     = 3.0*BS1-2.0*BS2
        A2P    = -3.0*B1P+2.0*B2P
        T1     = fourtr*(A1-0.5*A1P)
        X1     = A1P/T1-0.5
        T2     = fourtr*(A2-0.5*A2P)
        X2     = A2P/T2-0.5
        TYUK   = fourtr*(BSYUK-0.5*BYUKP)
        IF(TYUK.NE.ZERO) THEN
         XYUK   = BYUKP/TYUK-0.5
        ELSE
         XYUK   = ZERO
        ENDIF
       ENDIF
C
C      WRITE IN NEW FORM
C
       IF(IFTINN.EQ.1) THEN
C                Reassign to JP
        PAR(JP(1))  = T0
        PAR(JP(2))  = T1
        PAR(JP(3))  = T2
        PAR(JP(4))  = T3
        PAR(JP(5))  = T4
        PAR(JP(6))  = X0
        PAR(JP(7))  = X1
        PAR(JP(8))  = X2
        PAR(JP(9))  = X3
        PAR(JP(10)) = POWER
        PAR(JP(11)) = TYUK
        PAR(JP(12)) = XYUK
        PAR(JP(13)) = AYUK
        DPARMO(JP(1))  = 0.2
        DPARMO(JP(2))  = 0.2
        DPARMO(JP(3))  = 0.2
        DPARMO(JP(4))  = 2.0
        DPARMO(JP(5))  = 0.2
        DPARMO(JP(6))  = 0.0002
        DPARMO(JP(7))  = 0.0002
        DPARMO(JP(8))  = 0.0002
        DPARMO(JP(9))  = 0.0002
        DPARMO(JP(11))  = 0.2
        DPARMO(JP(12))  = 0.0002
       ELSE
C                   Reassign to JP
        IF(IFTINN.EQ.0) THEN
         PAR(JP(1))  = B0
         PAR(JP(4))  = B3
         PAR(JP(6))  = B0P
         PAR(JP(9))  = B3P
        ELSE
         PAR(JP(1))  = C0
         PAR(JP(4))  = C3
         PAR(JP(6))  = C0P
         PAR(JP(9))  = C3P
        ENDIF
        PAR(JP(2))  = B1
        PAR(JP(3))  = B2
        PAR(JP(5))  = T4
        PAR(JP(7))  = B1P
        PAR(JP(8))  = B2P
        PAR(JP(10)) = POWER
        PAR(JP(11)) = BYUK
        PAR(JP(12)) = BYUKP
        PAR(JP(13)) = AYUK
        DPARMO(JP(1)) = 0.2
        DPARMO(JP(2)) = 0.2
        DPARMO(JP(3)) = 0.2
        DPARMO(JP(4)) = 2.0
        DPARMO(JP(6)) = 2.0
        DPARMO(JP(7)) = 2.0
        DPARMO(JP(8)) = 2.0
        DPARMO(JP(9)) = 20.0
        DPARMO(JP(11)) = 0.2
        DPARMO(JP(12)) = 2.0
       ENDIF
C
C      OUTPUT OF RESULTING PARAMETERS IN FTN01-INPUT FORMAT
C
       DO 1019 J=1,NPARIN
       WRITE(02,'(G13.6,A,G12.5,A,I2,A,A)')  PAR(J),',',DPARMO(J),',',
     &                 IPARVA(J),',',ZUORD(IPARVA(J),IFTINN)
 1019  CONTINUE
C
    5 CONTINUE
C
 9999 STOP
      END
