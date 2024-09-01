c-----curfit---------------------------------------------------------
c
      program curfit
      IMPLICIT REAL*8(A-H,O-Z)
c
c     Fit a curve to given data.
c     The data come from the first and second column in standard input.
c     The input is copied to output and the fitted parameters are appended.
c     The model to be fitted is contained in the function 'model'.
c     A sample input is provided in 'curfit_sample.inp'.
c
c     The routine 'model' has to be provided separately and combined
c     with 'curfit'  at compile time. The callinc structure of 'model'
c     is 'call model(parameters,nparmx,yfit,ndatmx)' with
c        parameters   =  the actual model parameters
c        nparmx       =  maximum number of parameters
c        yfit         =  the 'y'-values resulting from the fit
c        ndatmx       =  maxmimum number of data points
c     The 'x'-values for the data are communicated via 
c        common /times/tim(KDATMX),delom
c     where 'tim' contains the sequency of 'x'-values.
c     (The 'delom' is obsolete.) This common block needs to be included
c     in the subroutine 'model'.
c
c     The maximum number of data point is KDATMX to set at compile time.
c     The maximum number of model parameters is KPARMX to set at compile 
c     time.
C
C
C
      parameter (KPARMX=20, KDATMX=20)
      parameter (KTESTC=1)
C
      CHARACTER*80 aline
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
      common /datinp/rhodat(KDATMX)
      dimension IRAND(30)
      DATA IRAND/1,1,1,0,1,1,0,1,0,1,1,0,1,0,0,1,1,0,1,1,1,0,0,
     *           1,1,0,1,0,1,0/
      data r0,r1/0.0,1.0/
C
c
c    initial settings
c
      IWU = 7
      itmax  = 30000
      IFCURV = 500
      epsrel = 1D-15
      epsabs = r1
      open(unit=7,form='formatted',status='unknown',file='for007')
      open(unit=8,form='formatted',status='unknown',file='curfit.cmp')
c
c     read parameters and experimental data from FOR005
c
      read(5,'(a)',err=999,end=999) aline
      write(IWU,'(a)') aline
      write(8,'(a)') aline
      do 1111 iread1=1,9999
      do 1111 iread2=1,9999
      read(5,'(a)',err=999,end=999) aline
      write(IWU,'(a)') aline
      if(aline(1:1).eq.'.') then
       read(5,'(a)') aline
       do 10 iread=1,KDATMX
       read(5,*,err=19,end=19) rhodat(iread),Y(iread),sigmay(iread)
       write(6,'(3f12.4)')
       NDATMX = iread
   10  continue
       write(IWU,'(a)') ' !! KDATMX exhausted.'
   19  continue
c
c      setting initial values
c
       do iread=1,KDATMX
         a(iread) = 0.0
       enddo
       do iread=1,KDATMX
         read(5,*,err=29,end=29) a(iread),dparmo(iread)
         nparmx = iread
       enddo
       stop ' file count exhausted'
 29    continue
c       nparmx = 3
c       a(1) = 10.0
c       dparmo(1) = 0.01
c       a(2) = -5.0
c       dparmo(2) = 0.01
c       a(3) = -1.0
c       dparmo(3) = 0.01 
c       a(4) = -1
c       dparmo(4) = 0.02
       write(IWU,'(a,i5/a,i5,a,2f10.4)')
     &   ' Nr. of data points=',ndatmx,
     &   ' Nr. of parameters=',nparmx,', parameters=',a(1),a(2)
c
      NWDR   = itmax    ! = maximum number of iterations !
      ISTOCH = 0
      ISTCHM = 3
      IFLAST = 0
      IFALT2 = 0
      chisqr = 1.0e30
      IF(NWDR.LT.2) THEN
       IST1 =1
      ELSE
       IST1 =0
      ENDIF
      do 100 iter=1,itmax
      if(ist1.gt.0) write(IWU,'(/a,i4)')  ' Iteration Nr.',iter
      CHIA=CHISQR
      WRITE(IWU,'(A,I3)')
     &  ' Next Iteration starts now. Iteration Nr.=',ITER
c      write(6,'(a,i5,3f6.1)')
c     &  'iter,parameters=',iter,a(1),a(2),a(3)
      call FIT(A,dparmo,Y,YFIT,SIGMAY,SIGMAA,SIGMA0,
     &         NPARMX,NDATMX,IFALT2,NWDR,IWARN,IFLAST,IFCURV,
     &         chisqr,IWU,IST1,
     &         ALPHA,zarray,B,beta,PARMOM,CHIMOM,PARSAV,
     &         YFITMO,FDERIJ,sumchi)
       write(6,'(a,i4,a,1pg12.4)')
     &   ' Iteration Nr.',iter,'  Chi^2=',chisqr
       write(7,'(a,i4,a,1pg12.4)')
     &   ' Iteration Nr.',iter,'  Chi^2=',chisqr
       chitst = chisqr
C
C      IWARN=99 INDICATES THAT THE GRADIENT STEP RUN INTO
C      BAD CONDITIONS. THE PARAMETERS ARE CHANGED STOCHASTICALLY
C      AND THE FIT IS CONTINUED.
C
       IF(IWARN.EQ.99) THEN
        ISTOCH = 1+ISTOCH
        IF(ISTOCH.GT.ISTCHM) THEN
        WRITE(IWU,'(/A)') ' NUMBER OF STOCHASTIC ATTEMPTS EXHAUSTED'
         GOTO 199
        ENDIF
        DO 359 J=1,NPARMX
  359   A(J) = A(J) + (2*IRAND(ISTOCH+J)-1)*DPARMO(J)*2.0
        IFALT2 = 0
        CHITST=.1E30
        WRITE(IWU,'(/A/)') ' IWARN=99. STOCHASTIC CHANGE OF PARAMS.'
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
C
C     This ends the iteration section
C
      write(IWU,'(/2a/a,i4,2(a,1pg12.4),a,1pg12.4)')
     &  ' Fit for case=',Titel,
     &  ' NPARMX=',nparmx,'  SIGMAY=',sigmay(1),'   CHI**2=',chisqr,
     &  '  partial CHI**2=',sumchi(1)
      write(IWU,'(a/a)')  ' Results:',
     &  '  i-1   parameter   corr.err.   unc.err.'
      do i=1,nparmx
      write(IWU,'(i5,1pg14.7,2(1pg11.3))') i-1,A(i),Sigmaa(i),Sigma0(i)
      enddo
      write(iwu,'(a/a)')  ' Fitted data:',
     &  '     rho      Yfit      Ydata   '
      do i=1,ndatmx
        write(iwu,'(1x,3g12.4)') rhodat(i),yfit(i),y(i)
      enddo
      write(8,'(a,i5,a,g12.5)')
     &  ' Convergence after ',iter,' iterations. CHISQR=',chisqr
      write(8,'(a/a)')  ' Results:',
     &  '  i-1   parameter   corr.err.   unc.err.'
      do i=1,nparmx
       write(8,'(i5,1pg14.7,2(1pg11.3))') i-1,A(i),Sigmaa(i),Sigma0(i)
      enddo
c
      endif
c
 1111 continue
c
  999 continue
      stop ' normal end'
      end
c-----chifct-----------------------------------------------------------
c
      subroutine CHIFCT(A,Y,SIGMAY,YFIT,NPARMX,NDATMX,CHISU,SUMCHI)
      IMPLICIT REAL*8(A-H,O-Z)
c
c     computation if chi**2
c
c     Input:
c      A      = array of model parameters (1...nparmx)
c      NPARMX = number of model parameters
c      Y      = array of experimental data
c      SIGMAY = array of adopted errors on data
c      NDATMX = number of data points
c      
c     Output:
c      YFIT   = fitted data
c      CHISU  = total chi**2
c 
c      INCLUDE 'PARAMS.INC'
      parameter (KPARMX=20, KDATMX=20)
      parameter (KTESTC=1)
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
      common /datinp/rhodat(KDATMX)
c
c     sum chi
c
c      write(6,'(a,4f6.1)') 
c     &  'chifct: ',(A(i),i=1,nparmx)
      call model(A,nparmx,yfit,ndatmx)
      sum   = 0.0
      do 10 i=1,ndatmx
      sum   = ((y(i)-yfit(i))/sigmay(i))**2  + sum
   10 continue
      CHISU = sum
      sumchi(1) = sum
c
      return
      end
c-----chifct-----------------------------------------------------------
c
      subroutine model(A,nparmx,yfit,ndatmx)
      IMPLICIT REAL*8(A-H,O-Z)
c
      parameter (KPARMX=20, KDATMX=20)
      DIMENSION YFIT(KDATMX),A(KPARMX)
      common /datinp/rhodat(KDATMX)
c
c-------------------------------------------------------------------
c
      do irho=1,ndatmx
        rhoact = rhodat(irho)
        sum = A(1)*rhoact
        do ipar=2,nparmx
          sum = sum + A(ipar)*rhoact**(ipar)
        enddo
        yfit(irho) = sum
      enddo

      return
      end

c--------------------------------------------------------------------
c
c     input from Baldo:
c#  Baldo data in Symmetric Nuclear Matter
c.
c.
c   0.04    -8.40315  2
c   0.08   -11.7278   1
c   0.12   -14.5467   0.05
c   0.16   -15.4002   0.02
c   0.20   -14.2367   0.05
c   0.24   -11.1535   0.1
c   0.28   -6.10037   0.4
c   0.32   0.956474   1
c   0.36    10.0410   5
c   0.40    21.1712   10
c   0.44    34.3608   15
c   0.48    49.6205   25
c   0.52    66.9592   50
c...
c   -100  1
c   200  2
c   0    4
c...




C-----FIT   ---------------------------------------------------------
C
      SUBROUTINE FIT(A,dparmo,Y,YFIT,SIGMAY,SIGMAA,SIGMA0,
     &               NPARMX,NDATMX,IFALT2,NWDR,IWARN,IFLAST,IFCURV,
     &               chisqr,IWU,IST1,
     &               ALPHA,zarray,B,beta,PARMOM,CHIMOM,PARSAV,
     &               YFITMO,FDERIJ,sumchi)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     CURFIT FROM BEVINGTON , COMPLEMENTED BY A CURVED STEP.
C
C     The routine is one step in minimizing a least-squares functional
C                             theory            exp        2
C         2                / f      (ndat,par)-f   (ndat) \
C      chi (par) = sum     | ---------------------------- |
C                     ndat \        Delta_f (ndat)        /
C
C     where 'par' is an array of model paramaters,
C     'ndat' runs over a given set of data points, and
C     'Delta_f' is the adopted error on data point number 'ndat'.
C
C     The chi^2 is to be computed in a soubroutine
C       CHIFCT(PARMOM,Y,SIGMAY,YFITMO,NPARMX,CHISU,SUMCHI)
C     where PARMOM is the input string for the model parameters 'par',
C     Y for the experimental data 'f^exp', and SIGMAY for the error on the
C     data 'Delta_f'. The number of model parameters is input as 'NPARMX'.
C     The routine should give as output: the model result 'f^theory' for
C     the given paramaters on the array YFITMO, the final 'chi^2' on
C     CHISU (!), and optionally the chi^2 from the single data points
C     on SUMCHI.
C
C
C
C
C      THE PARAMETERS IN THE CALLING SEQUENCE ARE:
C       A       = 'par'. INITIAL PARAMETERS AS INPUT, NEW PARAMETERS AFTER
C                 ONE ITERATION STEP GIVEN BACK AS OUTPUT.
C       DPARMO  = VARIATION OF PARAMETERS FOR DERIVATIVES
C       Y       = 'f^exp'. FIELD CONTAINING THE EXPERIMENTAL DATA,.
C       YFIT    = 'f^theory'. FIELD CONTAINING THE THEORETICAL RESULTS IN THE
C                  SAME MANNER AS THE EXPERIMENTAL DATA ON Y.
C       SIGMAY  = 'Delta_f'. FIELD CONTAINING THE CHOOSEN UNCERTAINTY ON THE
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
C
C
C     The full iteration is to be done in the calling routine with a
C     typical sequence like this:
C       -->  initialise the model dependent input
C              A, DPARMO, Y, SIGMAY, NPARMX, and NDATMX.
CC
CC     Initialise numerical parameters for the routine CHIFIT
CC
C      NWDR   = itmax    ! = maximum number of iterations !
C      ISTOCH = 0
C      IFCURV = 20
C      IFLAST = 0
C      IFALT2 = 0
C      chisqr = 1.0e30
C      IF(NWDR.LT.2) THEN
C       IST1 =1
C      ELSE
C       IST1 =0
C      ENDIF
C      do 100 iter=1,itmax
C      if(ist1.gt.0) write(6,'(/a,i4)')  ' Iteration Nr.',iter
C      CHIA=CHISQR
C      WRITE(IWU,'(A,I3)')
C     &  ' Next Iteration starts now. Iteration Nr.=',ITER
C      call FIT(A,dparmo,Y,YFIT,SIGMAY,SIGMAA,SIGMA0,
C     &         NPARMX,NDATMX,IFALT2,NWDR,IWARN,IFLAST,IFCURV,
C     &         chisqr,IWU,IST1,
C     &         ALPHA,zarray,B,beta,PARMOM,CHIMOM,PARSAV,
C     &         YFITMO,FDERIJ)
C     chitst = chisqr
CC
CC      IWARN=99 INDICATES THAT THE GRADIENT STEP RUN INTO
CC      BAD CONDITIONS. THE PARAMETERS ARE CHANGED STOCHASTICALLY
CC      AND THE FIT IS CONTINUED.
CC
C       IF(IWARN.EQ.99) THEN
C        DO 359 J=1,NPARMX
C  359   A(J) = A(J) + (2*IRAND(ISTOCH+J)-1)*DPARMO(J)*2.0
C        ISTOCH = 1+ISTOCH
C        IFALT2 = 0
C        CHITST=.1E30
C        WRITE(IWU,'(/A/)') ' IWARN=99. STOCHASTIC CHANGE OF PARAMS.'
C       ENDIF
C       IF(ISTOCH.GT.ISTCHM) THEN
C        WRITE(IWU,'(/A)') ' NUMBER OF STOCHASTIC ATTEMPTS EXHAUSTED\'
C        GOTO 199
C       ENDIF
C      chimax = .1e30
C      if((chia.ge.chitst .and. chia/chitst.lt.1.+epsrel
C     &    .and. chitst.lt.chimax)
C     &   .or. chitst.lt.epsabs) goto 199
Cc
C  100 continue
Cc
C  199 continue
C      if(iter.lt.itmax) write(IWU,'(a,i5,a,g12.5)')
C     &  ' Convergence after ',iter,' iterations. CHISQR=',chisqr
CC
CC     This ends the iteration section
CC
C
C      The stochastic change of parameters requires a field IRAND.
C      This can be given e.g. as
C
C
      parameter (KPARMX=20, KDATMX=10)
      parameter (KTESTC=1)
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
      data r0,r1/0.0,1.0/
      data ierror/0/
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
     *                  ' ====> STOP '
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
       CALL CHIFCT(A,Y,SIGMAY,YFIT,NPARMX,NDATMX,CHISU,SUMCHI)
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
      CALL CHIFCT(PARMOM,Y,SIGMAY,YFITMO,NPARMX,NDATMX,CHISU,SUMCHI)
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
       write(6,'(a,i2,a,a,(/10g12.5))')
     &   ' test FDERIJ for JPAR=',JPAR,' :'
c     &   ,'   YFIT=',(YFIT(imom),imom=1,NDATMX)
c       write(6,'(a,(/10g12.5))')
c     &   ' YFITMO=',(YFITMO(imom),imom=1,NDATMX)
       write(6,'(a,(/10g12.5))')  ' a,a_var=',a(jpar),parmom(jpar)
      endif
       DO 130 IDAT=1,NDATMX
       FDERIJ(IDAT,JPAR)=(YFIT(IDAT)-YFITMO(IDAT))/
     *                    (A(JPAR)-PARMOM(JPAR))
  130  CONTINUE
C
C      ZWISCHENAUSDRUCKE
C
      IF(IST1.GT.1) THEN
c       WRITE(IWU,'(I3,8E14.6,E13.3/3X,7E14.6,E12.4,I3)')
c     *      JPAR,(PARMOM(KL),KL=1,8),DPARMO(JPAR),
c     *      (PARMOM(KL),KL=9,KPARMX),CHIMOM(JPAR),ICHI
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
       B(J)=r0
       DO 330 K=1,NPARMX
  330  B(J)=B(J)+BETA(K)*ZARRAY(J,K)
     *         /SQRT(MAX(XLOWLM,ALPHA(J,J)*ALPHA(K,K)))
       bmax   = r0
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
      CALL CHIFCT(B,Y,SIGMAY,YFITMO,NPARMX,NDATMX,CHISU,SUMCHI)
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
      IF(VERGLC.LT.r0) THEN
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
       DO 710 J=1,KPARMX
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
       CALL CHIFCT(PARMOM,Y,SIGMAY,YFITMO,NPARMX,NDATMX,CHISUS,SUMCHI)
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
       CALL CHIFCT(PARMOM,Y,SIGMAY,YFITMO,NPARMX,NDATMX,CHISU,SUMCHI)
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
C-----MATIN2----------------------------------------------------------
C
      SUBROUTINE  MATIN2 (A,NMAX,NDIM,NERROR,DETERM)
      IMPLICIT REAL*8(A-H,O-Z)
C
c      INCLUDE 'PARAMS.INC'
      parameter (KPARMX=20, KDATMX=10)
      PARAMETER (KINDEX=KPARMX*KPARMX)
C
C     INVERTS THE MATRIX 'C'.
C
C     THE PARAMETERS ON THE I/O LIST ARE:
C
C       A       = AT INPUT:  MATRIX TO BE INVERTED
C                 AT OUTPUT: INVERSE MATRIX
C       NMAX    = NUMBER OF ACTUAL ROWS AND COLUMNS
C       NDIM    = DIMENSION OF MATRIX IN THE CALLING ROUTINE
C       NERROR  = ERROR FLAG
C       DETERM  = DETERMINANT OF INPUT MATRIX 'A'
C
C       THE ROUTINE IS TAKEN FROM 'CERNLIB' AND MODIFIED FOR THIS
C       PACKAGE.
C
      DIMENSION A(KINDEX),INDEX(KINDEX)
      data r0,r1/0.0,1.0/
C
C-----------------------------------------------------------------------------
C
      IF(NMAX.GT.KINDEX) STOP ' TOO LARGE DIMENSION IN MATIN2'
C
      DETER=1.00
      N=NMAX
C     IEMAT=N+N2
      IEMAT=N+0
      IDIM=NMAX
      NMIN1=N-1
C
C   UMSETZEN DES NMAX*NMAX-ANTEILS DER MATRIX A IN GEPACKTE
C   EIN-DIMENSIONALE FORM.
C
      DO 400 I=2,IDIM
      DO 400 J=1,IDIM
      IPACK  = (I-1)*NMAX + J
      IMATR  = (I-1)*NDIM + J
      A(IPACK) = A(IMATR)
  400 CONTINUE
C        THE ROUTINE DOES ITS OWN EVALUATION FOR DOUBLE SUBSCRIPTING OF
C        ARRAY A.
      IPIVC=1-IDIM
C        MAIN LOOP TO INVERT THE MATRIX
      DO 11 MAIN=1,N
      PIVOT=0.00
      IPIVC=IPIVC+IDIM
C        SEARCH FOR NEXT PIVOT IN COLUMN MAIN.
      IPIVC1=IPIVC+MAIN-1
      IPIVC2=IPIVC +NMIN1
      DO 2 I1=IPIVC1,IPIVC2
      IF(ABS(A(I1))-ABS(PIVOT)) 2,2,1
    1 PIVOT=A(I1)
      LPIV=I1
    2 CONTINUE
C        IS PIVOT DIFFERENT FROM ZERO
      IF(PIVOT) 3,15,3
C        GET THE PIVOT-LINE INDICATOR AND SWAP LINES IF NECESSARY
    3 ICOL=LPIV-IPIVC+1
      INDEX(MAIN)=ICOL
      IF(ICOL-MAIN) 6,6,4
C        COMPLEMENT THE DETERMINANT
    4 DETER=-DETER
C        POINTER TO LINE PIVOT FOUND
      ICOL=ICOL-IDIM
C        POINTER TO EXACT PIVOT LINE
      I3=MAIN-IDIM
      DO 5 I=1,IEMAT
      ICOL=ICOL+IDIM
      I3=I3+IDIM
      SWAP=A(I3)
      A(I3)=A(ICOL)
    5 A(ICOL)=SWAP
C        COMPUTE DETERMINANT
    6 DETER=DETER*PIVOT
      PIVOT=1./PIVOT
C        TRANSFORM PIVOT COLUMN
      I3=IPIVC+NMIN1
      DO 7 I=IPIVC,I3
    7 A(I)=-A(I)*PIVOT
      A(IPIVC1)=PIVOT
C        PIVOT ELEMENT TRANSFORMED
C
C        NOW CONVERT REST OF THE MATRIX
      I1=MAIN-IDIM
C        POINTER TO PIVOT LINE ELEMENTS
      ICOL=1-IDIM
C        GENERAL COLUMN POINTER
      DO 10 I=1,IEMAT
      ICOL=ICOL+IDIM
      I1=I1+IDIM
C        POINTERS MOVED
      IF(I-MAIN) 8,10,8
C        PIVOT COLUMN EXCLUDED
    8 JCOL=ICOL+NMIN1
      SWAP=A(I1)
      I3=IPIVC-1
      DO 9 I2=ICOL,JCOL
      I3=I3+1
    9 A(I2)=A(I2)+SWAP*A(I3)
      A(I1)=SWAP*PIVOT
   10 CONTINUE
   11 CONTINUE
C        NOW REARRANGE THE MATRIX TO GET RIGHT INVERS
      DO 14 I1=1,N
      MAIN=N+1-I1
      LPIV=INDEX(MAIN)
      IF(LPIV-MAIN) 12,14,12
   12 ICOL=(LPIV-1)*IDIM+1
      JCOL=ICOL+NMIN1
      IPIVC=(MAIN-1)*IDIM+1-ICOL
      DO 13 I2=ICOL,JCOL
      I3=I2+IPIVC
      SWAP=A(I2)
      A(I2)=A(I3)
   13 A(I3)=SWAP
   14 CONTINUE
      DETERM=DETER
      NERROR=0
C   RUECKSETZEN DER GEPACKTEN EINDIMENSIONALEN FORM IN DEN
C   NMAX*NMAX-ANTEIL DER MATRIX C
C
      DO 500 I=IDIM,2,-1
      DO 500 J=IDIM,1,-1
      IPACK  = (I-1)*NMAX + J
      IMATR  = (I-1)*NDIM + J
  500 A(IMATR) = A(IPACK)
C
      RETURN
   15 NERROR=MAIN
      DETERM=DETER
      RETURN
      END
c
c      include 'lev8.f'
c      include 'sp.f'
