c-----nucvar----------------------------------------------------------
c
      subroutine nucvar(par,dparmo,jp,nparmx)
c
c     computes extrapolated errors on nuclear matter properties
c
c      par    = array of model parameters
c      dparmo = array of variations of model parameters
c      nparmx = number of varied parameters
c               0 --> no evaluation of extrapolated errors
c
      include 'params.inc'
c
      DIMENSION DEOVRA(15),DRHONM(15),DCOMPR(15),DEFMAS(15),
     *          DSYMEN(15),dskapp(15)
      dimension par(15),dparmo(15),jp(15)
c
      include 'comfit.inc'
c
c----------------------------------------------------------------------
c
      DO 10 I=0,NPARMX
      IF(I.GT.0) PAR(I) = PAR(I) + DPARMO(I)
      CALL ATTACH(PAR,JP)
      call nucmat(rhonm0,eovra0,compr0,efmas0,glan0,symen0,d2asy0,
     &                skapp0,1,-2)
C
      IF(I.EQ.0) THEN
       EOVERA = EOVRA0
       RHONM  = RHONM0
       COMPR  = COMPR0
       EFFMAS = EFMAS0
       SYMENR = SYMEN0
       skappa = skapp0
      ELSE
       DEOVRA(I) = (EOVRA0-EOVERA)/DPARMO(I)
       DRHONM(I) = (RHONM0-RHONM )/DPARMO(I)
       DCOMPR(I) = (COMPR0-COMPR )/DPARMO(I)
       DEFMAS(I) = (EFMAS0-EFFMAS)/DPARMO(I)
       DSYMEN(I) = (SYMEN0-SYMENR)/DPARMO(I)
       dskapp(i) = (skapp0-skappa)/dparmo(i)
      ENDIF
C
      IF(I.GT.0) PAR(I) = PAR(I) - DPARMO(I)
   10 CONTINUE
C
C     ACCUMULATE EXTRAPOLATION ERRORS BY FOLDING WITH 'ZARRAY'.
C
      ERRRHO = 0.0
      ERRCMP = 0.0
      ERRMAS = 0.0
      ERRENR = 0.0
      ERRSYM = 0.0
      errska = 0.0
C
      IF(NPARMX.GT.0) THEN
       DO 40 I=1,NPARMX
       DO 40 J=1,NPARMX
       ERRRHO = DRHONM(I)*DRHONM(J)*ZARRAY(I,J) + ERRRHO
       ERRCMP = DCOMPR(I)*DCOMPR(J)*ZARRAY(I,J) + ERRCMP
       ERRMAS = DEFMAS(I)*DEFMAS(J)*ZARRAY(I,J) + ERRMAS
       ERRENR = DEOVRA(I)*DEOVRA(J)*ZARRAY(I,J) + ERRENR
       ERRSYM = DSYMEN(I)*DSYMEN(J)*ZARRAY(I,J) + ERRSYM
       ERRSka = DSkapp(I)*DSkapp(J)*ZARRAY(I,J) + ERRSka
   40  CONTINUE
       ERRRHO = SQRT(ERRRHO)
       ERRCMP = SQRT(ERRCMP)
       ERRMAS = SQRT(ERRMAS)
       ERRENR = SQRT(ERRENR)
       ERRSYM = SQRT(ERRSYM)
       ERRSka = SQRT(ERRSka)
      ENDIF
C
C     FINAL PRINT
C
      write(02,'(a/a/1x,f8.4,f8.2,f8.1,f8.4,f8.2,f8.4)')
     &     ' the nuclear matter parameters:',
     &     '     rho     e/a     compr.     effmas      a4    skappa',
     &     rhonm0,eovera,compr,effmas,asymme,skappa
      if(nparmx.gt.0) write(02,'(1x,f8.4,f8.2,f8.1,f8.4,f8.2,f8.4)')
     &     errrho,errenr,errcmp,errmas,errsym,errska
c
c     final call to print the Landau parameters
c
      call nucmat(rhonm0,eovra0,compr0,efmas0,glan0,symen0,d2asy0,
     &                skapp0,1,-1)
c
      return
      end
