      character*80 a
      real*8 h2min
      common /mass/h2min
      open(11,form='formatted',file='asymt.011')
      do 100 iread=1,99
      read(5,'(a)',err=99,end=99) a
      write(6,'(a/a)')
     &  ' running semi-infinite nuclear matter for the force:',
     &  a
      write(11,'(a/a)')
     &  ' running semi-infinite nuclear matter for the force:',
     &  a
      read(5,*,err=99,end=99) t0,t1,t2,t3,t4
      read(5,*,err=99,end=99) x0,x1,x2,x3,b4p
      read(5,*,err=99,end=99) power,h2min
c
c     reshape extended l*s to effective T4
c
c      T4_eff = (t4+b4p)*2/3
c
      if(abs(b4p).gt.1.0E-7) T4 = (t4+b4p)*0.666666667
C
      call ASYMT(T0,T1,T2,T3,t4,X0,X1,X2,X3,power,
     +YRHO0,eovera,xKOMPr,a4,YEFFM,YXKF,YELLO,
     +surfe,surfas,YH,YP,YG,YQT1,YQT2)
C    OUTPUT PARAMETERS ARE:
C     ES0    : SURFACE ENERGY COEFFICIENT
C     XJO    : ASYMMETRY ENERGY J
C     Q      : SURFACE ASYMMETRY
C     H,P,G  : AUXILIARY EXPANSION COEFFICIENTS DEFINING Q, SEE
C              PHYSICS REPORTS "SEMICLASSICAL ..."
C     QT1,QT2 : ALTERNATIVE Q-COEFFICIENTS;
C               THE MOST RELEVANT OF THE COEFFICIENTS Q, QT1 AND QT2
C               IS THE COEFFICIENT QT1 !
C               THE DEVIATION BETWEEN Q, QT1 AND QT2 INDICATES SOME
C               INSUFFICIENCY OF THE PARAMETRISATION OF THE
C               DENSITY PROFILE ( PROBLEMS YET TO BE SETTLED).
C     ELLO   : EXPANSION COEFFICIENT L, STANDING IN FRONT OF THE
C              TERM WITH DELTA AND EPSILON IN THE LDM EXPANSION.
C
      write(6,'(6(/a,1pg12.4)/a,3(1pg12.4))')
     &  ' density=',yrho0,
     &  '     E/A=',eovera,
     &  ' incomp.=',xkompr,
     &  ' asymm.E=',a4,
     &  '    m*/m=',yeffm,
     &  '  surf.E=',surfe,
     &  ' surf.As=',surfas,YQT1,YQT2
      write(11,'(6(/a,1pg12.4)/a,3(1pg12.4))')
     &  ' density=',yrho0,
     &  '     E/A=',eovera,
     &  ' incomp.=',xkompr,
     &  ' asymm.E=',a4,
     &  '    m*/m=',yeffm,
     &  '  surf.E=',surfe,
     &  ' surf.As=',surfas,YQT1,YQT2
  100 continue
   99 stop
      end
C-----ASYMT ---------------------------------------------------------
C
C ****  ASYMT  ****     VERSION  10.5.84  ****
C
C     THIS IS ONLY THE MASTER ROUTINE, HANDLING INPUT DATA;
C     IT TRANSFERS CONTROL QUICKLY TO THE SUBROUTINE 'ASYMSB'
C     WHICH DOES ALL OF THE FURTHER WORK.
C     THE INPUT PARAMETERS ARE ALL NAMED ENDING WITH AN 'L';
C     THEY ARE TRANSFERRED VIA LIST AND PUT ONTO THE SAME NAMES,
C     WITHOUT 'L', IN THE VARIOUS COMMON BLOCKS.
C     INPUT AND OUTPUT PARAMETERS ARE EXPLAINED IN DETAIL IN
C     SUBROUTINE 'ASYMSB' .
C
C     IMPLICIT REAL*8(A-H),REAL*8(O-Z)
C
      SUBROUTINE ASYMT(YT0L,YT1L,YT2L,YT3L,YW0L,YX0L,
     +YX1L,YX2L,YX3L,YALFAL,
     +YRHO0,YELAM,YBKOMP,YXJO,YEFFM,YXKF,YELLO,
     +YES0,YQ,YH,YP,YG,YQT1,YQT2)
      IMPLICIT DOUBLE PRECISION (A-H,O-X,Z)
      T0L=YT0L
      T1L=YT1L
      T2L=YT2L
      T3L=YT3L
      W0L=YW0L
      X0L=YX0L
      X1L=YX1L
      X2L=YX2L
      X3L=YX3L
      ALFAL=YALFAL
C
C
C
C
C--------------------------------------------------------------------
C
C     OPEN(UNIT=05,ACCESS='SEQUENTIEL',FORM='FORMATTED',
C    *     FILE='ASYMDAT')
C     OPEN(UNIT=06,ACCESS='SEQUENTIEL',FORM='FORMATTED',
C    *     FILE='RESULT')
C 100 CONTINUE
CW    WRITE(02,8)
C
C     ISKY : NUMBER FOR SKYRME-FORCE (SEE SUBROUTINE SKYRM (LINE 5510))
C            0      FOR STOP
C            -1     FOR NEW SKYRME-FORCES
C     N    : SWITCH PARAMETER                              USUAL -2
C     Z1   : LEFT INTEGRATION BOUNDARY                     USUAL -20
C     Z2   : RIGHT INTEGRATION BOUNDARY                    USUAL +20
C     DZ   : INTEGRATION STEP WIDTH                        USUAL 0.5
C     APS  : TOTAL ACCURACY (MEV)                          USUAL 1.E-4 /
C     TSK  : ONLY NEEDED FOR TEMP<>0                       USUAL 0.0
C     DDD  : DELTA                                         USUAL 0.01
C     DDT  : TAU                                           USUAL 0.01
C     TEMP : TEMPERATURE                                   USUAL 0.0
C
      ISKYL=-1
      NL=-2
      Z1L=-20.
      Z2L=+20.
      DZL=0.5
      APSL=1.E-4
      TSKL=0.0
      DDDL=0.01
      DDTL=0.01
      TEMPL=0.0
      AJL=0.
      BJL=0.
      APL=0.5
      ANL=0.5
      GAPL=1.0
      GANL=1.0
      DELTAL=1.0
C     READ(05,*) ISKYL,NL,Z1L,Z2L,DZL,APSL,TSKL,DDDL,DDTL,TEMPL
C     IF(ISKYL.EQ.0) STOP
CW    WRITE(02,9)
C     IF(ISKYL.LT.0)
C    *   READ(05,*) T0L,T1L,T2L,T3L,X0L,X1L,X2L,X3L,W0L,ALFAL,AJL,BJL
CW    WRITE(02,7)
C
C    INITIALISATION FOR ITERATION IS READ IN:
C     AP    : ALFA PROTON      STARTING VALUE              USUAL 0.5
C     AN    : ALFA NEUTRON     STARTING VALUE              USUAL 0.5
C     GAP   : GAMMA PROTON                              FIXED TO 1.0
C     GAN   : GAMMA NEUTRON                             FIXED TO 1.0
C     DELTA : DELTA            STARTING VALUE           FIXED TO 1.0
C
C     READ(05,*) APL,ANL,GAPL,GANL,DELTAL
      CALL ASYMSB(ISKYL,NL,Z1L,Z2L,DZL,APSL,TSKL,DDDL,DDTL,TEMPL,
     *            T0L,T1L,T2L,T3L,X0L,X1L,X2L,X3L,W0L,ALFAL,AJL,BJL,
     *            APL,ANL,GAPL,GANL,DELTAL,
     +            RHO0,ELAM,BKOMP,EFFM,XKF,
     *            ES0,XJO,Q,H,P,G,QT1,QT2,ELLO)
CW    WRITE(02,'(6(A,G12.5)/11X,5(A,G12.5))')
CW   *    ' TESTPRINT:   ES0=',   ES0,'   XJO=',   XJO,'     Q=',Q,
CW   *               '     H=',     H,'     P=',     P,'     G=',G,
CW   *               '   QT1=',   QT1,'   QT2=',   QT2,'  ELLO=',ELLO
      YRHO0=RHO0
      YELAM=ELAM
      YBKOMP=BKOMP
      YXJO=XJO
      YEFFM=EFFM
      YXKF=XKF
      YELLO=ELLO
      YES0=ES0
      YQ=Q
      YH=H
      YP=P
      YG=G
      YQT1=QT1
      YQT2=QT2
C     GOTO 100
C
    7 FORMAT(' READ: AP,AN,GAP,GAN,DELTA')
    8 FORMAT(' READ: ISKY,N,Z1,Z2,DZ,APS,TSK,DDD,DDT,TEMP')
    9 FORMAT(' READ: T0,T1,T2,T3,X0,X1,X2,X3,W0,ALFA,AJ,BJ')
C
      RETURN
      END
C-----ASYMSB---------------------------------------------------------
C
C     IMPLICIT REAL*8(A-H),REAL*8(O-Z)
C
      SUBROUTINE ASYMSB(ISKYL,NL,Z1L,Z2L,DZL,APSL,TSKL,DDDL,DDTL,TEMPL,
     *                  T0L,T1L,T2L,T3L,X0L,X1L,X2L,X3L,W0L,ALFAL,
     *                  AJL,BJL,  APL,ANL,GAPL,GANL,DELTAL,
     +           ARHO0,AELAM,ABKOMP,AEFFM,AXKF,
     *                  ES0,XJO,Q,H,P,G,QT1,QT2,ELLO)
C
C     ALL LIST PARAMETERS END WITH AN 'L'. THEY ARE ASSIGNED TO
C     THE COMMON-PARAMETERS WIDHT THE SAME NAME, BUT WITHOUT 'L',
C     AT THE BEGINNING OF THIS ROUTINE.
C     MANY PARAMETERS ARE NUMERICAL PARAMETERS AND NEED NOT MUCH
C     CHANGE. TYPICAL VALUES ARE:
C     ISKY : NUMBER FOR SKYRME-FORCE (SEE SUBROUTINE SKYRM (LINE 5510))
C            0      FOR STOP
C            -1     FOR NEW SKYRME-FORCES
C     N    : SWITCH PARAMETER                              USUAL -2
C     Z1   : LEFT INTEGRATION BOUNDARY                     USUAL -20
C     Z2   : RIGHT INTEGRATION BOUNDARY                    USUAL +20
C     DZ   : INTEGRATION STEP WIDTH                        USUAL 0.5
C     APS  : TOTAL ACCURACY (MEV)                          USUAL 1.E-4 /
C     TSK  : ONLY NEEDED FOR TEMP<>0                       USUAL 0.0
C     DDD  : DELTA                                         USUAL 0.01
C     DDT  : TAU                                           USUAL 0.01
C     TEMP : TEMPERATURE                                   USUAL 0.0
C    INITIALISATION FOR ITERATION IS READ IN:
C     AP    : ALFA PROTON      STARTING VALUE              USUAL 0.5
C     AN    : ALFA NEUTRON     STARTING VALUE              USUAL 0.5
C     GAP   : GAMMA PROTON                              FIXED TO 1.0
C     GAN   : GAMMA NEUTRON                             FIXED TO 1.0
C     DELTA : DELTA            STARTING VALUE           FIXED TO 1.0
C
C    OUTPUT PARAMETERS ARE:
C     ES0    : SURFACE ENERGY COEFFICIENT
C     XJO    : ASYMMETRY ENERGY J
C     Q      : SURFACE ASYMMETRY
C     H,P,G  : AUXILIARY EXPANSION COEFFICIENTS DEFINING Q, SEE
C              PHYSICS REPORTS "SEMICLASSICAL ..."
C     QT1,QT2 : ALTERNATIVE Q-COEFFICIENTS;
C               THE MOST RELEVANT OF THE COEFFICIENTS Q, QT1 AND QT2
C               IS THE COEFFICIENT QT1 !
C               THE DEVIATION BETWEEN Q, QT1 AND QT2 INDICATES SOME
C               INSUFFICIENCY OF THE PARAMETRISATION OF THE
C               DENSITY PROFILE ( PROBLEMS YET TO BE SETTLED).
C     ELLO   : EXPANSION COEFFICIENT L, STANDING IN FRONT OF THE
C              TERM WITH DELTA AND EPSILON IN THE LDM EXPANSION.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      EXTERNAL ES
C
      COMMON /CONST1/ TSK,YTA,AVI,RHO0,ELAM,ELL,RAA0,BKOMP,XJ,R0,
     +EFFM,XKF
      COMMON /SKY/ T0,T1,T2,T3,X0,X1,X2,X3,ALFA,W0,ISKY,TEMP,AJ,BJ
      COMMON/RANGE/Z1,Z2,DZ,IPOINT,N
C
      DIMENSION PA(20),DD(20),P1(5),P2(5),P3(5),P4(5)
      DATA PA/20*0.0/
C
      ISKY   = ISKYL
      N      = NL
      Z1     = Z1L
      Z2     = Z2L
      DZ     = DZL
      APS    = APSL
      TSK    = TSKL
      DDD    = DDDL
      DDT    = DDTL
      TEMP   = TEMPL
      T0     = T0L
      T1     = T1L
      T2     = T2L
      T3     = T3L
      X0     = X0L
      X1     = X1L
      X2     = X2L
      X3     = X3L
      W0     = W0L
      ALFA   = ALFAL
      AJ     = AJL
      BJ     = BJL
      PA(1)  = APL
      PA(2)  = ANL
      PA(3)  = GAPL
      PA(4)  = GANL
      PA(5)  = DELTAL
c      write(6,'(a)') ' in ASYMTB before INITIA'
      CALL INITIA
c      write(6,'(a)') ' in ASYMTB after INITIA'
      WRITE(02,11) T0,T1,T2,T3,X0,X1,X2,X3,W0,ALFA,AJ,BJ
      JUMP=0
      IF(N.EQ.-5) JUMP=1
      IF(N.EQ.-2) JUMP=2
      N=IABS(N)
      IPOINT=((Z2-Z1)/DZ+1.00001)
CW    WRITE(02,12) ISKY,N,Z1,Z2,DZ,APS,TSK,TEMP
CW    WRITE(02,13) (PA(I),I=1,5)
C
      IF(JUMP.EQ.1) GO TO 30
      DO 1 I=1,20
    1 DD(I)=0.0005
      DD(3)=-.05
      DD(4)=-.05
      AL=1.0
      BE=0.5
      GA=2.0
      IF(JUMP.EQ.2) GO TO 40
      IF(N.LE.4) CALL LOKLAM(PA(5))
      CALL SIMPLE(N,PA,DD,AL,BE,GA,APS,ESU,JZAHL,ES)
c      write(6,'(a)') ' in ASYMTB after 1. SIMPLE'
CW    WRITE(02,14) JZAHL,ESU,TSK,(PA(I),I=1,5)
      RETURN
   30 CONTINUE
      CALL LOKLAM(PA(5))
      CALL SEMIC(PA(1),PA(2),PA(3),PA(4),PA(5),ETOT2,ETOT4)
CW    WRITE(02,17) ETOT2,ETOT4
      RETURN
   40 CONTINUE
      TSK=0.
      PA(5)=0.
      CALL LOKLAM(PA(5))
      CALL SIMPLE(N,PA,DD,AL,BE,GA,APS,ES0,JZAHL,ES)
c      write(6,'(a)') ' in ASYMTB after 2. SIMPLE'
CW    WRITE(02,14) JZAHL,ES0,TSK,(PA(I),I=1,5)
      TSK=DDT
      PA(5)=DDD
      CALL LOKLAM (PA(5))
      CALL SIMPLE(N,PA,DD,AL,BE,GA,APS,ES3,JZAHL,ES)
c      write(6,'(a)') ' in ASYMTB after 3. SIMPLE'
CW    WRITE(02,14) JZAHL,ES3,TSK,(PA(I),I=1,5)
      TSK=0.
      PA(5)=DDD
      CALL LOKLAM(PA(5))
      CALL SIMPLE(N,PA,DD,AL,BE,GA,APS,ES2,JZAHL,ES)
c      write(6,'(a)') ' in ASYMTB after 4. SIMPLE'
CW    WRITE(02,14) JZAHL,ES2,TSK,(PA(I),I=1,5)
      TSK=DDT
      PA(5)=0.
      ALMI=0.5*(PA(1)+PA(2))
      PA(1)=ALMI
      PA(2)=ALMI
      CALL LOKLAM(PA(5))
      CALL SIMPLE(N,PA,DD,AL,BE,GA,APS,ES1,JZAHL,ES)
c      write(6,'(a)') ' in ASYMTB after 5. SIMPLE'
CW    WRITE(02,14) JZAHL,ES1,TSK,(PA(I),I=1,5)
      TAU=DDT/R0
      H=(ES1-ES0)/TAU**2
      G=(ES0-ES2)/DDD**2
      P=(ES3-ES2-ES1+ES0)/(2.*TAU*DDD)
      Q=H+P**2/G
      QT1=H/(1.-2.*P/(3.*XJ))
      QT2=1.5*XJ*P/G
CW    WRITE(02,18) ES0,XJ,Q,H,P,G,QT1,QT2
      XJO    = XJ
      ELLO   = ELL
C
      ARHO0=RHO0
      AELAM=ELAM
      ABKOMP=BKOMP
      AEFFM=EFFM
      AXKF=XKF
      RETURN
C  10 FORMAT(80B)
   11 FORMAT(1H0,'SKYRME PARAMETERS: T0,T1,T2,T3=',4F10.3/
     1 20X,'X0,X1,X2,X3=',4F10.3/20X,'W0,ALFA,AJ,BJ=',4F10.3/)
   12 FORMAT(1H0,'ISKY,N,Z1,Z2,DZ,EPS,TSK,TEMP=',2I3,6F10.5)
   13 FORMAT(1H0,'STARTING VALUES = ',10F10.5/)
   14 FORMAT(' ITER ESU :',I4,F10.6//' TSK=',F10.5,
     1 ' PARAMETERS:',10F10.5,/)
   15 FORMAT(1X,4F10.5,5X,4F10.5)
   16 FORMAT(1H0)
   17 FORMAT(1H0,'ETOT2,ETOT4:',2F10.4)
   18 FORMAT(1H0,' DROPLET PARAMETERS:  AS=',F6.2,' J=',F6.2,
     1' Q=',F6.2,' H,P,G=',3F6.2,'TEST: Q1,Q2=',2F6.2,/)
      END
C-----INITIA---------------------------------------------------------
C
      SUBROUTINE INITIA
C     IMPLICIT REAL*8(A-H),REAL*8(O-Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /SKY/ T0,T1,T2,T3,X0,X1,X2,X3,ALFA,W0,ISKY,TEMP,AJ,BJ
      COMMON /CONST1/ TSK,YTA,AVI,RHO0,ELAM,ELL,RAA0,BKOMP,XJ,R0,
     +EFFM,XKF
      common /mass/h2min
      COMMON /CONST2/ RKAPP,H2M,C53,C83,TOND,CT,PI4,T,ETA,ETAP,ETAN
      COMMON /RANGE/ Z1,Z2,DZ,IPOINT,N
      COMMON IR,IW
      PI=3.1415926536
      PIPI=PI*PI
      PI4=4.0*PI
      PY=(PIPI*1.5)**(2.0/3.0)
      CALL SKYRM
      IF(ISKY.NE.4) AJ=T2-T1
      IF(ISKY.NE.4) BJ=X1*T1+X2*T2
cold      H2M=20.73
      H2M=h2min
      H22M=1.0/H2M
      EINS=1.
      ZERO=0.
      C13=1.0/3.0
      C53=5.0/3.0
      C23=2.0/3.0
      C83=8.0/3.0
      CM13=-1.0/3.0
      C1A=1.0/ALFA
      C35=1.0/C53
      C85=8.0/5.0
      RKAPP=H2M*C35*(3.0*PIPI)**C23
      TOND=0.0
      IF(ISKY.EQ.16) TOND=1.0
      CT=67.766
      HJ1=H2M*PY*C35
      HJ2=T0/C83
      HJ3=T3/16.0
      HJ4=ALFA+2.0
      HJ5=1.+ALFA
      AETA=H2M/36.0
      BETA=H22M*(3.0*T1+5.0*T2+4.0*T2*X2)/16.0
      ZETA=H2M/12.0*BETA**2
      BBB=(9.0*T1-5.0*T2)/64.0
C     IF(ISKY.EQ.7) BETA=0.0
      SIGMA=BBB-11.0*AETA*BETA
      XI=H2M*3.0/5.0*PY
      BS1=3.0*T0/8.0
      BS2=T3*HJ5/16.
      BS3=BS2*HJ4
      CSO=(0.75*W0)**2*0.5*H22M
      CTE=SQRT(PY)*TEMP**2/(H2M*18.)
      ETA=0.
      DDR=0.
      T=TEMP
      IF(T.LE.0.001) GO TO 21
      FACST=5./(3.*T)
      FACJT=T*(T/H2M)**1.5/PIPI
      ETA=46./T
C
C*****DETERMINATION OF RHO AND LAMBDA
C
   21 SIG=1.0
      IF(DZ.LT.0.0) SIG=-1.0
      IT=0
      R1=0.16
    1 IT=IT+1
      EFF=1.+BETA*R1
      IF(T.LE.0.001) GO TO 22
      CALL FINDET(EINS,R1,R1,BETA,ZERO,T,ETA,FJ1,FJ3,NITT,DDR,DDRM,YX,
     1YXXX)
      TK=FACJT*FJ3/EFF**1.5
      X=BS1*R1**2+BS2*R1**HJ4+TK*(C53-1./EFF)
      XP=2.*BS1*R1+BS3*R1**HJ5+2.5*BETA*TK*(1.-EFF)/EFF**2
     1  +DDR*T*R1*(2.5-1.5/EFF)
      GO TO 23
   22 TK=XI*R1**C53*EFF
      X=C23*TK+XI*BETA*R1**C83+BS1*R1**2+BS2*R1**HJ4
      XP=10.*(TK/(9.*R1)+C13*XI*BETA*R1**C53)+2.*BS1*R1+BS3*R1**HJ5
   23 R2=R1-X/XP
CW    WRITE(02,33) IT,R2,X,XP,NITT,ETA,DDR
      IF(ABS(R1-R2).LT.1.0E-6) GO TO 2
      IF(IT.GT.50) STOP
   33 FORMAT(' IT,R2,X,XP=',I3,F9.6,2E12.4,
     1' NITT,ETA,DDR=',I3,2E15.6)
      R1=R2
      GO TO 1
    2 CONTINUE
      RHO=R2
      RHO23=RHO**C23
      DDLA=0.
      IF(T.GT.0.001) DDLA=-C53*TK/RHO+T*ETA
      ELL=10.0/9.0*H2M*PY*C35*RHO23-0.75*T0*(0.5+X0)*RHO+
     +25.0/3.0*C35*PY*(1.0/6.0*T2*(1.0+1.25*X2)-
     -1.0/8.0*X1*T1)*RHO**C53-(1.0+ALFA)/16.0*T3*
     *(1.0+2.0*X3)*RHO**(1.0+ALFA)
      BKK=-2.*TK/RHO+12.*XI*BETA*RHO**C53
      IF(T.GT.0.001) BKK=
     1  TK*(-30.+18./EFF)/RHO+9.*T*RHO*DDR*(2.5-1.5/EFF)
      BKOMP=BKK+9.*ALFA*BS2*RHO**HJ5
      ELAM=TK/RHO+BS1*RHO+HJ3*RHO**HJ5+DDLA
      F0=1.0+BETA
    3 CONTINUE
      XJ=H2M*PY/3.0*RHO23-0.25*T0*(0.5+X0)*RHO+PY*RHO**C53*
     *(1.0/6.0*T2*(1.0+1.25*X2)-X1*T1/8.0)-T3/48.0*(1.0+2.0*X3)*
     *RHO**(1.0+ALFA)-DDLA/9.
      IF(ISKY.EQ.16) XJ=RHO**C23*(CT+5.0/9.0*PY*H2M*C35)-DDLA/9.
      IF(ISKY.EQ.16) ELL=2.0*XJ
      RHO0=RHO
      EFFM=1.0/EFF
      XKF=(1.5*PIPI*RHO)**(1.0/3.0)
      R0=(3.0/PI4/RHO)**C13
      YTA=PI4*R0*R0
CW    WRITE(02,4) RHO,ELAM,DDLA,BKOMP,XJ,EFFM,XKF,ELL
    4 FORMAT(1H0,' RHOINF,LAM,DLAT,K,J,M*,KF,L=',9F10.5/)
      RETURN
      END
C-----SIMPLE---------------------------------------------------------
C
      SUBROUTINE SIMPLE(N,X,D,ALPHA,BETA,GAMMA,EPS,FNEU,JZAHL,FKT)
C     IMPLICIT REAL*8(A-H),REAL*8(O-Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(20),XX(20,21),YY(20),XNULL(20),XR(20),XE(20),XC(20)
      DIMENSION F(20),D(N )
      COMMON IR,IW
      EXTERNAL FKT
      YY(4)=X(4)
      YY(5)=X(5)
      XR(4)=X(4)
      XR(5)=X(5)
      XE(4)=X(4)
      XE(5)=X(5)
      XC(4)=X(4)
      XC(5)=X(5)
      YY(3)=X(3)
      XR(3)=X(3)
      XE(3)=X(3)
      XC(3)=X(3)
      NP1=N+1
      ALPH1=1+ALPHA
      BETA1=1-BETA
      GAMM1=1-GAMMA
C
C     ERSTELLEN DER KOEFFIZIENTENMATRIX
C
      JZAHL=0
      DO 1 K=1,NP1
      DO 1 I=1,N
      XX(I,K)=X(I)
      IF(K.EQ.I) XX(I,K)=XX(I,K)+D(I)
    1 CONTINUE
C
C     BERECHNUNG DER DAZUGEHOERIGEN FUNKTIONSWERTEN
C
      DO 3 I=1,NP1
      DO 4 K=1,N
    4 YY(K)=XX(K,I)
    3 F(I)=FKT(N,YY)
C
C     MINIMUM
   55 LL=1
      FL=F(1)
      DO 5 I=2,NP1
      IF((F(I)-FL).GE.0.) GO TO 5
      FL=F(I)
      LL=I
    5 CONTINUE
C
C     MAXIMUM
   66 LH=1
      FH=F(1)
      DO 6 I=2,NP1
      IF((F(I)-FH).LE.0.) GO TO 6
      FH=F(I)
      LH=I
    6 CONTINUE
C
C     ZWEITES MAXIMUM
   77 II=2
      LH2=1
      FH2=F(1)
      IF(LH-1) 74,75,74
   75 LH2=2
      II=3
      FH2=F(2)
   74 DO 7 I=II,NP1
      IF(I.EQ.LH) GO TO 7
      IF((F(I)-FH2).LE.0.) GO TO 7
      FH2=F(I)
      LH2=I
    7 CONTINUE
      JZAHL=JZAHL+1
      FQUER=0.
      DO 8 K=1,NP1
    8 FQUER=FQUER+F(K)
      FQUER=FQUER/NP1
      S=0.
      DO 9 I=1,NP1
    9 S=S+(F(I)-FQUER)**2
      S=SQRT(S/N)
      IF(S.LE.EPS)   GO TO 216
      GO TO 213
216   DO 11 I=1,N
11    X(I)=XX(I,LL)
      FNEU=F(LL)
       RETURN
C
C
  213 DO 46 I=1,N
      XNULL(I)=0.
      DO 47 K=1,NP1
      IF(K.EQ.LH) GO TO 47
      XNULL(I)=XNULL(I)+XX(I,K)
   47 CONTINUE
   46 XNULL(I)=XNULL(I)/N
      DO 13 K=1,N
   13 XR(K)=ALPH1*XNULL(K)-ALPHA*XX(K,LH)
      FR=FKT(N,XR)
      IF(FR.GE.F(LL)) GO TO 26
      DO 16 K=1,N
   16 XE(K)=GAMMA*XR(K)+GAMM1*XNULL(K)
      FE=FKT(N,XE)
      IF(FE.GE.F(LL)) GO TO 23
      DO 19 I=1,N
   19 XX(I,LH)=XE(I)
      F(LH)=FE
   20 LL=LH
      LH=LH2
      GO TO 77
   23 DO 24 I=1,N
   24 XX(I,LH)=XR(I)
      F(LH)=FR
      GO TO 20
   26 IF(FR.GT.F(LH2)) GO TO 30
      DO 27 I=1,N
   27 XX(I,LH)=XR(I)
      F(LH)=FR
      LH=LH2
      GO TO 77
   30 IF(FR.GE.F(LH)) GO TO 32
      DO 31 I=1,N
   31 XX(I,LH)=XR(I)
      F(LH)=FR
   32 DO 33 I=1,N
   33 XC(I)=BETA*XX(I,LH)+BETA1*XNULL(I)
      FC=FKT(N,XC)
      IF(FC.GE.F(LH)) GO TO 36
      DO 35 I=1,N
   35 XX(I,LH)=XC(I)
      F(LH)=FC
      GO TO 55
   36 DO 42 I=1,N
      DO 42 K=1,NP1
      IF(K.EQ.LL) GO TO  42
       XX(I,K)=(XX(I,K)+XX(I,LL))*0.5
   42 CONTINUE
      DO 37 I=1,NP1
      DO 38 K=1,N
   38 YY(K)=XX(K,I)
   37 F(I)=FKT(N,YY)
      GO TO  66
      END
C-----SEMIC ---------------------------------------------------------
C
      SUBROUTINE SEMIC(AP,AN,GAP,GAN,DELTA,ETOT2,ETOT4)
C     IMPLICIT REAL*8(A-H),REAL*8(O-Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONST1/ TSK,YTA,AVI,RHO0,ELAM,ELL,RAA0,BKOMP,XJ,R0,
     +EFFM,XKF
      COMMON /CONST2/ RKAPP,H2M,C53,C83,TOND,CT,PI4,T,ETA,ETAP,ETAN
      COMMON /SKY/ T0,T1,T2,T3,X0,X1,X2,X3,ALFA,W0,ISKY,TEMP,AJ,BJ
      COMMON /RANGE/ Z1,Z2,DZ,IPOINT,N
      COMMON IR,IW
      DIMENSION F(30)
      CORC=0.0
      PI=3.141592654
      H2M=20.73
      DO 2 I=1,30
    2 F(I)=0.0
      X11=1.+X1/2.
      X12=0.5+X1
      X21=1.+X2/2.
      X22=0.5+X2
      X31=1.+X3/2.
      X32=0.5+X3
      T5=(T1*(1.-X1)+3.*T2*(1.+X2))/8./H2M
      T6=(T1*X11+T2*X21)/4./H2M
      T52=T5*T5
      T62=T6*T6
      T56=2.*T5*T6
      W02=W0*W0
      T01=(1.+0.5*X0)*T0*0.5
      T02=-(0.5+X0)*T0*0.5
      T7=(T2*X21-3.*T1*X11)/16.
      T8=(T2*X22+3.*T1*X12)/16.
      T31=T3/12.
      FKAP=(3.*PI*PI)**(2./3.)
      FKAP2=SQRT(FKAP)/12./H2M
      CTEM=TEMP**2*FKAP2
      TOND=0.0
      IF(ISKY.EQ.16) TOND=1.0
      CT=67.766
      C23=2.0/3.0
      FAKT=H2M*FKAP*0.6
      TWO=2.
      FACJT=0.5*T*(T/H2M)**1.5/PI**2
      Z=Z1-DZ
      DO 1 IZ=1,IPOINT
      Z=Z+DZ
      CALL DERI(Z,AP,AN,DELTA,GAP,GAN,DP,DN,DRP,DRN,DLP,DLN)
      D=DP+DN
      FN=1.+T5*DN+T6*DP
      FP=1.+T5*DP+T6*DN
      FNI=1./FN
      FPI=1./FP
      DNI=1./DN
      DPI=1./DP
      FN2=FN**2
      FP2=FP**2
      DP2=DP**2
      DN2=DN**2
      DNP=DN*DP
      D2=D**2
      DALFA=D**ALFA
      DL=DLP+DLN
      GN=DRN**2
      GP=DRP**2
      HN=GN/DN2
      HP=GP/DP2
      GNP=DRN*DRP
      GFN=T52*GN+T62*GP+GNP*T56
      GFP=T52*GP+T62*GN+GNP*T56
      GHN=T5*GN+T6*GNP
      GHP=T5*GP+T6*GNP
      FLN=T5*DLN+T6*DLP
      FLP=T5*DLP+T6*DLN
      SN=W02*(GN+0.25*GP+GNP)
      SP=W02*(GP+0.25*GN+GNP)
      DSN=W0*(GN+0.5*GNP)
      DSP=W0*(GP+0.5*GNP)
      FSN=W0*(T5*GN+0.5*T6*GP+(T5/2.+T6)*GNP)
      FSP=W0*(T5*GP+0.5*T6*GN+(T5/2.+T6)*GNP)
      GSN=W0*(DLN+0.5*DLP)
      GSP=W0*(DLP+0.5*DLN)
      DNP13=DN**0.33333333
      DPP13=DP**0.33333333
C****  CENTRAL POTENTIAL SKYRME ENERGY
      F1=T01*D2+T02*(DN2+DP2)+T31*DALFA*(X31*D2-X32*(DN2+DP2))
     1   +T7*D*DL+T8*(DN*DLN+DP*DLP)+TOND*CT*(DN-DP)**2/D**(1./3.)
C**   KINETIC ENERGY 0TH ORDER
      DN53=DN**(5./3.)
      DP53=DP**(5./3.)
      TKP=FAKT*(FN*DN53+FP*DP53-CORC*(DN53+DP53))
      IF(T.LE.0.001) GOTO 22
      CALL FINDET(TWO,DP,DN,T5,T6,T,ETAP,F1P,F3P,NITP,WP,WWP,YNP,ZEP)
      CALL FINDET(TWO,DN,DP,T5,T6,T,ETAN,F1N,F3N,NITN,WN,WWN,YNN,ZEN)
      GETP=ZEP-YNP
      GETN=ZEN-YNN
      DETP=(9.*ZEP-7./12.)/4.
      DETN=(9.*ZEN-7./12.)/4.
      TKP=-2.*FACJT*(F3N*FNI**1.5+F3P*FPI**1.5)/3.+T*(ETAP*DP+ETAN*DN)
  22  F(2)=F(2)+TKP+F1-AVI*D
C**   KINETIC ENERGY 2ND ORDER
      F(3)=F(3)+(FN*GN*DNI+FP*GP*DPI)/36.0-(GHN+GHP)/3.0-(GFN*DN*FNI
     1      +GFP*DP*FPI)/12.0
      IF(T.LE.0.001) GO TO 23
      F(3)=FN*GN*DNI*ZEN+FP*GP*DPI*ZEP+3.*(ZEN*GHN+ZEP*GHP)
     1 +5.*(DLN*FN+DLP*FP)/12. + (GFN*DN*FNI*DETN+GFP*DP*FPI*DETP)
C**  KINETIC ENERGY 4TH ORDER
  23  DN23=DNP13*DNI
      DP23=DPP13*DPI
      DP53=DP23*DPI
      DN53=DN23*DNI
      DN83=DN53*DNI
      DP83=DP53*DPI
      DN113=DN83*DNI
      DP113=DP83*DPI
      F(4)=F(4)+(DNP13*FNI)*(-FLN**2+GFN*(2.0*FLN-GFN*FNI)*FNI)
     1+(DPP13*FPI)*(-FLP**2+GFP*(2.0*FLP-GFP*FPI)*FPI)
      F(5)=F(5)+(DN23*FNI)*(DLN*(FN*FLN*2.-GFN)+FLN*GHN)
     1+(DP23*FPI)*(DLP*(FP*FLP*2.-GFP)+FLP*GHP)
      F(6)=F(6)+DN53*FNI*(7.*GHN*(GHN+2.*FN*DLN)+GN*(4.*GFN-7.*FN*FLN)
     1+8.*(FN*DLN)**2)
     1  +DP53*FPI*(7.*GHP*(GHP+2.*FP*DLP)+GP*(4.*GFP-7.*FP*FLP)
     1+8.*(FP*DLP)**2)
      F(7)=F(7)+
     +DN83*GN*(27.*FN*DLN+22.*GHN)+DP83*GP*(27.*FP*DLP+22.*GHP)
      F(8)=F(8)+DN113*FN*GN*GN+DP113*FP*GP*GP
      IF(ISKY.EQ.16) GO TO 20
C **   SPIN-ORBIT CONTRIBUTION TO 2ND ORDER **
      F(9)=F(9)+DN*SN*FNI+DP*SP*FPI
C **   SPIN-ORBIT CONTRIBUTIONS TO 4TH ORDER  **
      F(10)=F(10)+
     =DNP13/FN2*(GSN*(4.*GSN*FN-6.*FSN)+SN*FLN+2.*FSN*FSN*FNI)
     1  +DPP13/FP2*(GSP*(4.*GSP*FP-6.*FSP)+SP*FLP+2.*FSP*FSP*FPI)
      F(11)=F(11)+DN23/FN2*(DSN*2.*(GSN*FN-0.*FSN)+SN*(DLN*FN-0.*GHN))
     1  +DP23/FP2*(DSP*2.*(GSP*FP-0.*FSP)+SP*(DLP*FP-0.*GHP))
      F(12)=F(12)+DN53*FNI*SN*GN+DP53*FPI*SP*GP
      F(13)=F(13)+DNP13*FNI**3*SN*SN+DPP13*FPI**3*SP*SP
      GO TO 30
   20 CONTINUE
C****  S.-O. FOR THE TONDEUR FUNCTIONAL
      F(9)=F(9)+D*(GN+GP+2.0*GNP)
      F(10)=F(10)+DL*DL*(DNP13+DPP13)
      F(11)=F(11)+DL*(DN23*(GN+GNP)+DP23*(GP+GNP))
      F(12)=F(12)+(GN+GP+2.0*GNP)*(DLN*DN23+DLP*DP23-C23*(GN
     **DN53+GP*DP53))
      F(13)=F(13)+(GN+GP+2.0*GNP)**2*(DNP13+DPP13)
   30 CONTINUE
      F(14)=F(14)+(AJ+BJ)*(GN*DN2/FN2+GP*DP2/FP2)
     1 +2.*BJ*GNP*DNP/(FN*FP)
    1 CONTINUE
      DO 13 I=1,14
   13 F(I)=F(I)*DZ
      TFK=F(2)
      WZK=H2M*F(3)
      T4K1=F(4)/240.
      T4K2=F(5)/360.
      T4K3=F(6)/2160.
      T4K4=-F(7)/6480.
      T4K5=F(8)/810.
      T4K=(T4K1+T4K2+T4K3+T4K4+T4K5)*H2M/FKAP
      IF(ISKY.EQ.16) GO TO 40
      SO2=-(0.5/H2M)*F(9)
      SO41=F(10)/16.
      SO42=F(11)/48.
      SO43=-F(12)/72.
      SO44=0.5*F(13)
      SO4=(0.5/H2M)/FKAP*(SO41+SO42+SO43+SO44*(0.5/H2M)**2)
      GO TO 41
   40 CONTINUE
      SO2=-W02*0.5/H2M*F(9)
      SO4=W02/(FKAP*H2M)*0.5*(0.25*F(10)+F(11)/24.+F(12)/48.
     1 +F(13)*W02/(8.*H2M**2))
   41 CONTINUE
      SJ2=-F(14)*(0.75*W0/H2M)**2/16.
      IF(ISKY.EQ.4) SO2=SO2+SJ2
      ETOT2=(TFK+WZK+SO2)*YTA
      ETOT4=ETOT2+(T4K+SO4)*YTA
      SJ2=SJ2*YTA
CW    WRITE(02,3) DELTA,AP,AN,GAP,GAN,ETOT4,SJ2
    3 FORMAT(1X,'DELTA,AP,AN,GAP,GAN,ETOT4=',5F8.4,F11.6,
     1 ' SJ2=',F9.4)
      RETURN
      END
C-----SKYRM ---------------------------------------------------------
C
      SUBROUTINE SKYRM
C     IMPLICIT REAL*8(A-H),REAL*8(O-Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /SKY/ T0,T1,T2,T3,X0,X1,X2,X3,ALFA,W0,ISKY,TEMP,AJ,BJ
      IF(ISKY.LT.0) RETURN
      X1=0.
      X2=0.
      X3=1.
      ALFA=1.
      AJ=0.
      BJ=0.
      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), ISKY
    1 CONTINUE
      RETURN
    3 T0=-1128.75
      T1=395.0
      T2=-95.0
      T3=14000.
      X0=0.45
      W0=120.
      RETURN
    2 T0=-2645.0
C**   SKM*
      T1=410.
      T2=-135.
      T3=15595.
      X0=0.09
      X3=0.
      W0=130.
      ALFA=1./6.
      RETURN
    4 T0=-1788.9
      T1=301.5
      T2=502.5
      T3=12764.
      X0=0.353
      W0=130.
      X1=-2.5
      X2=-1.7
      X3=0.475
      ALFA=1./3.
      BJ=-575.
      RETURN
    5 T0=-2160.
C**   RATP
      T1=513.
      T2=121.
      T3=11600
      X0=0.418
      X1=-0.36
      X2=-2.29
      X3=0.586
      ALFA=0.2
      W0=120.
      RETURN
    6 CONTINUE
      RETURN
    7 CONTINUE
C**   K* FORCE OF P.G.REINHARD AND FRIEDRICH (1984)
      T0=-1998.23
      T1=375.992
      T2=-103.327
      T3=11931.6
      X0=-0.7418
      X3=-1.4985
      W0=126.515
      ALFA=0.25
      RETURN
    8 T0=-2645.
      T1=391.8
      T2=-124.1
      T3=15595.
      X0=0.09
      X3=0.
      W0=130.
      ALFA=1./6.
      RETURN
    9 T0=-2645.
C**   SKM
      T1=385.
      T2=-120.
      T3=15595.
      X0=0.09
      X3=0.
      W0=130.
      ALFA=1.0/6.0
      RETURN
 10   T0=-2563.0714
      T1=437.11642
      T2=-409.24024
      T3=14725.436
      X0=-0.22212679
      X1=-0.11986273
      X2=-0.86447481
      X3=-0.36613607
      ALFA=1./6.
      W0=120.0
      RETURN
C  FORCES OF TONDEUR  M*=M
   11 CONTINUE
      T0=-1794.0
      T1=298.0
      T2=-298.0
      T3=12812.0
      X0=0.154
      X1=-0.5
      X2=-0.5
      X3=0.089
      W0=110.0
      ALFA=1.0/3.0
      RETURN
   12 CONTINUE
      T0=-1791.6
      T1=300.0
      T2=-300.0
      T3=12792.0
      X0=0.154
      X1=-0.5
      X2=-0.5
      X3=0.089
      W0=120.0
      ALFA=1.0/3.0
      RETURN
   13 CONTINUE
      T0=-1791.8
      T1=298.5
      T2=-99.5
      T3=12794.0
      X0=0.138
      X1=-1.0
      X2=+1.0
      X3=0.075
      W0=126.0
      ALFA=1.0/3.0
      RETURN
   14 CONTINUE
      T0=-1808.8
      T1=303.4
      T2=-303.4
      T3=12980.0
      X0=-0.177
      X1=-0.5
      X2=-0.5
      X3=-0.5
      W0=113
      ALFA=1.0/3.0
      RETURN
   15 CONTINUE
      T0=-2917.1
      T1=328.2
      T2=-328.2
      T3=18584.0
      X0=-0.295
      X1=-0.5
      X2=-0.5
      X3=-0.5
      W0=114
      ALFA=1.0/6.0
      RETURN
   16 CONTINUE
C     FUNCTIONAL DUE TO TONDEUR NP A315(1978) 353
      X0=-0.5
      X1=-0.5
      X2=-0.5
      X3=-0.5
      T0=8.0/3.0*(-678.9)
      T1=16.0/3.0*58.25
      T2=-T1
      T3=16.0*812.07
      ALFA=1.0/3.0
      D=82.0
      W0=D
      RETURN
   17 CONTINUE
C     FORCE OF FRIEDRICH + REINHARD (MAINZ, MARCH 1982)
      T0=-1218.8
      T1=106.1
      T2=-183.
      T3=14100.
      X0=0.03
      X3=-0.43
      W0=110.4
      ALFA=2./3.
      RETURN
      END
C-----DERI  ---------------------------------------------------------
C
      SUBROUTINE DERI(Z,AP,AN,DELTA,GAP,GAN,DP,DN,DRP,DRN,DLP,DLN)
C     IMPLICIT REAL*8(A-H),REAL*8(O-Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONST1/ TSK,YTA,AVI,RHO0,ELAM,ELL,RAA0,BKOMP,XJ,R0,
     +EFFM,XKF
      TT=0.5*TSK
      TP=TT/AP
      TN=-TT/AN
      RHO0P=0.5*RAA0*(1.0-DELTA)
      RHO0N=0.5*RAA0*(1.0+DELTA)
      XIP=-GAP/AP
      XIN=-GAN/AN
      EP=EXP(Z/AP)*EXP(TP)
      EN=EXP(Z/AN)*EXP(TN)
      ETAP=EP/(1.0+EP)
      ETAN=EN/(1.0+EN)
      ETAPK=ETAP**2
      ETANK=ETAN**2
      ETAP1=(ETAP-ETAPK)/AP
      ETAN1=(ETAN-ETANK)/AN
      DP=RHO0P/(1.0+EP)**GAP
      DN=RHO0N/(1.0+EN)**GAN
      DRP=XIP*DP*ETAP
      DRN=XIN*DN*ETAN
      DLP=XIP*DP*(XIP*ETAPK+ETAP1)
      DLN=XIN*DN*(XIN*ETANK+ETAN1)
      RETURN
      END
C-----ES    ---------------------------------------------------------
C
      FUNCTION ES(N,Y)
C     IMPLICIT REAL*8(A-H),REAL*8(O-Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(20)
      COMMON /CONST1/ TSK,YTA,AVI,RHO0,ELAM,ELL,RAA0,BKOMP,XJ,R0,
     +EFFM,XKF
      IF(N.EQ.5) CALL LOKLAM(Y(5))
      CALL SEMIC(Y(1),Y(2),Y(3),Y(4),Y(5),ETOT2,ETOT4)
      ES=ETOT4
      RETURN
      END
C-----LOKLAM---------------------------------------------------------
C
      SUBROUTINE LOKLAM(DELTA)
C     IMPLICIT REAL*8(A-H),REAL*8(O-Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CONST1/ TSK,YTA,AVI,RHO0,ELAM,ELL,RAA0,BKOMP,XJ,R00,
     +EFFM,XKF
      COMMON /CONST2/ RKAPP,H2M,C53,C83,TOND,CT,PI4,T,ETA,ETAP,ETAN
      COMMON /SKY/ T0,T1,T2,T3,X0,X1,X2,X3,ALFA,W0,ISKY,TEMP,AJ,BJ
      COMMON IR,IW
      EPS=ELL/BKOMP*DELTA**2
      C13=1./3.
      TWO=2.
      PI=3.14159265
      FPI=(3.*PI**2)**C13
      CTE=TEMP**2*FPI/(12.*H2M)
      BETAQ=(T1*(1.-X1)+3.*T2*(1.+X2))/(8.*H2M)
      BETAMQ=(T1*(1.+0.5*X1)+T2*(1.+0.5*X2))/(4.*H2M)
      CM13=-1.0/3.0
      ETAP=0.
      ETAN=0.
      IF(T.LE.0.001) GOTO 21
      FACJT=0.5*T*(T/H2M)**1.5/PI**2
      ETAP=38./T
      ETAN=50./T
   21 RAA0=RHO0*(1.0-3.0*EPS)
      RHON=0.5*RAA0*(1.0+DELTA)
      RHOP=0.5*RAA0*(1.0-DELTA)
      FP=1.+(BETAQ*RHOP+BETAMQ*RHON)
      FN=1.+(BETAQ*RHON+BETAMQ*RHOP)
      RR53=RHON**C53+RHOP**C53
      RR2=RHON**2+RHOP**2
      TKP=RKAPP*(RR53+(0.25*(T1*(1.+0.5*X1)+T2*(1.+0.5*X2))*RR53*RAA0
     1   +0.25*(T2*(X2+0.5)-T1*(X1+0.5))*(RHON**C83+RHOP**C83))/H2M)
      IF(T.LE.0.001) GO TO 22
      CALL FINDET(TWO,RHOP,RHON,BETAQ,BETAMQ,T,ETAP,F1P,F3P,NITP,ZX,ZXX,
     1YNX,ZEX)
      CALL FINDET(TWO,RHON,RHOP,BETAQ,BETAMQ,T,ETAN,F1N,F3N,NITN,ZX,ZXX,
     1YNX,ZEX)
      TKP=-2.*FACJT*(F3P/FP**1.5+F3N/FN**1.5)/3.+T*(RHON*ETAN+RHOP*ETAP)
   22 AVI=TKP+0.5*T0*(RAA0**2*(1.0+0.5*X0)-RR2*(X0+0.5))
     1   +T3/12.0*RAA0**ALFA*(RAA0**2*(1.0+0.5*X3)-RR2*(X3+0.5))
     2   +TOND*CT*(RHON-RHOP)**2*RAA0**CM13
      AVI=AVI/RAA0
      R0=(3.0/PI4/RAA0)**0.3333333
CW    WRITE(02,100) RAA0,RHOP,RHON,ETAP,ETAN,AVI
  100 FORMAT(' LOKLAM:  RHO0,RHOP,RHON=',3F8.5,' ETAP,ETAN,AVI=',
     1 2E12.4,F10.5,/)
      YTA=PI4*R0**2
      RETURN
      END
C-----FINDET---------------------------------------------------------
C
      SUBROUTINE FINDET(FAC,RHO,RHOM,BET,GAM,T,ETA,F1,F3,NIT,DEDR,DEDRM,
     1YN,ZE)
C     IMPLICIT REAL*8(A-H),REAL*8(O-Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DATA H2M,PI  / 20.73, 3.141592653 /
      IPRIT=0
      RTEST=1.E-07
      FQ=1.+RHO*BET+RHOM*GAM
      F=RHO*FAC*PI*PI*(SQRT(H2M*FQ/T))**3
      NIT=0
    1 NIT=NIT+1
      P1=0.5
      F1=FERDIR(ETA,P1)
      DIFR=(F-F1)/F
      IF(ABS(DIFR).LE.RTEST) GO TO 2
      P2=-0.5
      FM1=FERDIR(ETA,P2)
      ETA=ETA+(F-F1)*2./FM1
      IF(IPRIT.LE.0) GO TO 1
CW    WRITE(02,10) NIT,ETA,F,F1
   10 FORMAT(' NIT,ETA,F,F1=',I3,3E14.6)
      GO TO 1
    2 P3=1.5
      F3=FERDIR(ETA,P3)
      FFMI=2.*F/FM1
      DEDR=FFMI*(1./RHO+1.5*BET/FQ)
      DEDRM=FFMI*1.5*GAM/FQ
      DD=0.01
      ETAP=ETA+DD
      ETAM=ETA-DD
      FM1=FERDIR(ETA,P2)
      FM1P=FERDIR(ETAP,P2)
      FM1M=FERDIR(ETAM,P2)
      FM3=(FM1M-FM1P)/DD
      FM5=4.*(FM1P+FM1M-2.*FM1)/(3.*DD**2)
      X=F1*FM3/FM1**2
      Y=F1**2*FM5/FM1**3
      YN=(X*(1.+2.*X)-3.*Y)/8.
      ZE=-X/12.
      RETURN
      END
C-----FERDIR---------------------------------------------------------
C
C/  ADD    NAME=FERDIR
      FUNCTION FERDIR(X,P)
C     IMPLICIT REAL*8(A-H),REAL*8(O-Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION P1(5),P2(5),P3(5),P4(5),P5(5),P6(5),P7(5),P8(5),P9(5)
      DIMENSION Q1(5),Q2(5),Q3(5),Q4(5),Q5(5),Q6(5),Q7(5),Q8(5),Q9(5)
      DATA P1
     1/-1.25331 41288 20E+0, -1.72366 35577 01E+0, -6.55904 57292 58E-1,
     2 -6.34228 31976 82E-2, -1.48838 31061 16E-5/
      DATA Q1
     1/+1.00000 00000 00E+0, +2.19178 09259 80E+0, +1.60581 29554 06E+0,
     2 +4.44366 95274 81E-1, +3.62423 22881 12E-2/
      DATA P2
     1/-3.13328 53055 70E-1, -4.16187 38522 93E-1, -1.50220 84005 88E-1,
     2 -1.33957 93751 73E-2, -1.51335 07001 38E-5/
      DATA Q2
     1/+1.00000 00000 00E+0, +1.87260 86759 02E+0, +1.14520 44465 78E+0,
     2 +2.57022 55875 73E-1, +1.63990 25435 68E-2/
      DATA P3
     1/-2.34996 39854 06E-1, -2.92737 36375 47E-1, -9.88309 75887 38E-2,
     2 -8.25138 63795 51E-3, -1.87438 41532 23E-5/
      DATA Q3
     1/+1.00000 00000 00E+0, +1.60859 71091 46E+0, +8.27528 95308 80E-1,
     2 +1.52232 23828 50E-1, +7.69512 04750 64E-3/
      DATA P4
     1/+1.07381 27694E+0, +5.60033 03660E+0, +3.68822 11270E+0,
     2 +1.17433 92816E+0, +2.36419 35527E-1/
      DATA Q4
     1/+1.00000 00000E+0, +4.60318 40667E+0, +4.30759 10674E-1,
     2 +4.21511 32145E-1, +1.18326 01601E-2/
      DATA P5
     1/+6.78176 62666 0E-1, +6.33124 01791 0E-1, +2.94479 65177 2E-1,
     2 +8.01320 71141 9E-2, +1.33918 21294 0E-2/
      DATA Q5
     1/+1.00000 00000 0E+0, +1.43740 40039 7E-1, +7.08662 14845 0E-2,
     2 +2.34579 49473 5E-3, -1.29449 92883 5E-5/
      DATA P6
     1/+1.15302 13402E+0, +1.05915 58972E+0, +4.68988 03095E-1,
     2 +1.18829 08784E-1, +1.94387 55787E-2/
      DATA Q6
     1/+1.00000 00000E+0, +3.73489 53841E-2, +2.32484 58137E-2,
     2 -1.37667 70874E-3, +4.64663 92781E-5/
      DATA P7
     1/-8.22255 9330E-1, -3.62036 9345E+1, -3.01538 5410E+3,
     2 -7.04987 1579E+4, -5.69814 5924E+4/
      DATA Q7
     1/+1.00000 0000E+0, +3.93568 9841E+1, +3.56875 6266E+3,
     2 +4.18189 3625E+4, +3.38513 8907E+5/
      DATA P8
     1/+8.22449 97626E-1, +2.00463 03393E+1, +1.82680 93446E+3,
     2 +1.22265 30374E+4, +1.40407 50092E+5/
      DATA Q8
     1/+1.00000 00000E+0, +2.34862 07659E+1, +2.20134 83743E+3,
     1 +1.14426 73596E+4, +1.65847 15900E+5/
      DATA P9
     1/+2.46740 02368 4E+0, +2.19167 58236 8E+2, +1.23829 37907 5E+4,
     2 +2.20667 72496 8E+5, +8.49442 92003 4E+5/
      DATA Q9
     1/+1.00000 00000 0E+0, +8.91125 14061 9E+1, +5.04575 66966 7E+3,
     2 +9.09075 94630 4E+4, +3.89960 91564 1E+5/
      IF(P .EQ. 1.5) GO TO 30
      IF(P .EQ. 0.5) GO TO 20
      IF(P .EQ.-0.5) GO TO 10
      FERDIR=0.
CW    WRITE(02,100) P
  100 FORMAT(1X,23HFERDIR ... ILLEGAL P =  ,F10.6)
      RETURN
   10 IF(X .GT. 4.0) GO TO 11
      IF(X .GT. 1.0) GO TO 12
      Y=EXP(X)
      FERDIR=Y*(1.772453850905516+Y*
     1       (P1(1)+Y*(P1(2)+Y*(P1(3)+Y*(P1(4)+Y*P1(5)))))/
     2       (Q1(1)+Y*(Q1(2)+Y*(Q1(3)+Y*(Q1(4)+Y*Q1(5))))))
      RETURN
   12 FERDIR=(P4(1)+X*(P4(2)+X*(P4(3)+X*(P4(4)+X*P4(5)))))/
     1       (Q4(1)+X*(Q4(2)+X*(Q4(3)+X*(Q4(4)+X*Q4(5)))))
      RETURN
   11 Y=1.0/X**2
      FERDIR=SQRT(X)*(2.0+Y*
     1       (P7(1)+Y*(P7(2)+Y*(P7(3)+Y*(P7(4)+Y*P7(5)))))/
     2       (Q7(1)+Y*(Q7(2)+Y*(Q7(3)+Y*(Q7(4)+Y*Q7(5))))))
      RETURN
   20 IF(X .GT. 4.0) GO TO 21
      IF(X .GT. 1.0) GO TO 22
      Y=EXP(X)
      FERDIR=Y*(0.886226925452758+Y*
     1       (P2(1)+Y*(P2(2)+Y*(P2(3)+Y*(P2(4)+Y*P2(5)))))/
     2       (Q2(1)+Y*(Q2(2)+Y*(Q2(3)+Y*(Q2(4)+Y*Q2(5))))))
      RETURN
   22 FERDIR=(P5(1)+X*(P5(2)+X*(P5(3)+X*(P5(4)+X*P5(5)))))/
     1       (Q5(1)+X*(Q5(2)+X*(Q5(3)+X*(Q5(4)+X*Q5(5)))))
      RETURN
   21 Y=1.0/X**2
      FERDIR=X*SQRT(X)*(0.666666666666667+Y*
     1       (P8(1)+Y*(P8(2)+Y*(P8(3)+Y*(P8(4)+Y*P8(5)))))/
     2       (Q8(1)+Y*(Q8(2)+Y*(Q8(3)+Y*(Q8(4)+Y*Q8(5))))))
      RETURN
   30 IF(X .GT. 4.0) GO TO 31
      IF(X .GT. 1.0) GO TO 32
      Y=EXP(X)
      FERDIR=Y*(1.329340388179137+Y*
     1       (P3(1)+Y*(P3(2)+Y*(P3(3)+Y*(P3(4)+Y*P3(5)))))/
     2       (Q3(1)+Y*(Q3(2)+Y*(Q3(3)+Y*(Q3(4)+Y*Q3(5))))))
      RETURN
   32 FERDIR=(P6(1)+X*(P6(2)+X*(P6(3)+X*(P6(4)+X*P6(5)))))/
     1       (Q6(1)+X*(Q6(2)+X*(Q6(3)+X*(Q6(4)+X*Q6(5)))))
      RETURN
   31 XS=X**2
      Y=1.0/XS
      FERDIR=XS*SQRT(X)*(0.4+Y*
     1       (P9(1)+Y*(P9(2)+Y*(P9(3)+Y*(P9(4)+Y*P9(5)))))/
     2       (Q9(1)+Y*(Q9(2)+Y*(Q9(3)+Y*(Q9(4)+Y*Q9(5))))))
      RETURN
      END
 
