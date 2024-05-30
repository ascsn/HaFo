c-----testrd-----------------------------------------------------------
c
      program testrd
c
c     tests occurence of reduced matrix elements
c
      INCLUDE 'PARAMS.INC'
c
      data LRPA/2/
      data one/1.0/
c
c----------------------------------------------------------------------
c
      write(6,'(a,i5/3a)')
     & ' Test Reduced Matrix Elements for LRPA=',LRPA,
     & ' XJ  XL XJA XLA    REDMAT  SPMDER  SPMRAD ',
     & ' SPPDER  SPPRAD  CRMDER  CRMRAD  CRPDER ',
     & ' CRPRAD  REDMSG'
      MLRPA = 2*LRPA+1
      ALRPA = LRPA
c
      do 10 jsa = 1,KJMX
      LA    = JSA/2
      XLA   = LA
      XJA   = (JSA-1)/2+0.5
C
C     LOOP OVER THE J,L TO WHICH THEY COUPLE
C
      DO 30 JDIF=1,MLRPA
      XJ     = XJA+JDIF-LRPA-1
c      IF(XJ.LT.0.0) GOTO 30
      LL     = (XJ-0.5)
      IF(MOD(LRPA+LA+LL,2).NE.0) THEN
       LL     = (XJ+0.5)
      ENDIF
      XL     = LL
c      IF(LL.GT.LRPA+LA .OR. LL.LT.ABS(LA-LRPA)) GO TO 30
C
C     DETERMINE COMBINED QUANTUM NUMBER TO 'XJ', 'XL'.
C
      IF(XJ.GT.XL) THEN
       JS     = XJ+XJ
      ELSE
       JS     = XJ+XJ+ONE
      ENDIF
c
c     compute and print matrix elements
c
      REDMAT = REDUME(XJ,ALRPA,XJA)
      CALL REDSPO(XJ,XL,ALRPA,ALrpa+one,XJA,XLA,SPPDER,SPPRAD)
      CALL REDSPO(XJ,XL,ALRPA,ALrpa-one,XJA,XLA,SPMDER,SPMRAD)
      CALL REDCUR(XJ,XL,ALRPA,ALrpa+one,XJA,XLA,CRPDER,CRPRAD)
      CALL REDCUR(XJ,XL,ALRPA,ALrpa-one,XJA,XLA,CRMDER,CRMRAD)
      REDMSG = REDSIG(XJ,XL,ALRPA,ALRPA,XJA,XLA)
      write(6,'(4f4.1,2x,10f8.3)') xj,xl,xja,xla,
     &   redmat,spmder,spmrad,sppder,spprad,crmder,crmrad,
     &   crpder,crprad,redmsg
c
   30 continue
   10 continue
c
      write(6,'(/a)')
     &   ' Reverse relations:'
c
      do 110 jsa = 1,KJMX
      LA    = JSA/2
      XLA   = LA
      XJA   = (JSA-1)/2+0.5
C
C     LOOP OVER THE J,L TO WHICH THEY COUPLE
C
      DO 130 JDIF=1,MLRPA
      XJ     = XJA+JDIF-LRPA-1
c      IF(XJ.LT.0.0) GOTO 30
      LL     = (XJ-0.5)
      IF(MOD(LRPA+LA+LL,2).NE.0) THEN
       LL     = (XJ+0.5)
      ENDIF
      XL     = LL
c      IF(LL.GT.LRPA+LA .OR. LL.LT.ABS(LA-LRPA)) GO TO 30
C
C     DETERMINE COMBINED QUANTUM NUMBER TO 'XJ', 'XL'.
C
      IF(XJ.GT.XL) THEN
       JS     = XJ+XJ
      ELSE
       JS     = XJ+XJ+ONE
      ENDIF
c
c     compute and print matrix elements
c
      REDMAT = REDUME(XJa,ALRPA,XJ)
      CALL REDSPO(XJa,XLa,ALRPA,ALrpa+one,XJ,XL,SPPDER,SPPRAD)
      CALL REDSPO(XJa,XLa,ALRPA,ALrpa-one,XJ,XL,SPMDER,SPMRAD)
      CALL REDCUR(XJa,XLa,ALRPA,ALrpa+one,XJ,XL,CRPDER,CRPRAD)
      CALL REDCUR(XJa,XLa,ALRPA,ALrpa-one,XJ,XL,CRMDER,CRMRAD)
      REDMSG = REDSIG(XJa,XLa,ALRPA,ALRPA,XJ,XL)
      write(6,'(4f4.1,2x,10f8.3)') xj,xl,xja,xla,
     &   redmat,spmder,spmrad,sppder,spprad,crmder,crmrad,
     &   crpder,crprad,redmsg
c
  130 continue
  110 continue
c
      stop
      end
