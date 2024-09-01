C INCLUDE BLOCK: CONSTN
C     SOME CONSTANT WHICH OCCUR FREQUENTLY:
C      PI     = 3.141592654 = AREA OF A UNIT CIRCLE
C      SQ2    = SQUARE ROOT OF 2.0
C      FOURPI = 4*PI = 12.56637061
C      D4PI   = 1.0/(4*PI)
C      SQFRPI = SQUARE ROOT OF 4*PI
C      XSMALL = VERY SMALL NUMBER ( 1.0E-20 )
C      ZERO   = 0.0
C      ONE    = 1.0
C      HALF   = 0.5
C
       parameter (one=1.0D0)
       parameter (zero=0.0D0)
       parameter (half=0.5D0)
       parameter (sq2=1.414213562D0)
c       parameter (xsmall=1D-20)
       COMMON /CONSTN/ PI,FOURPI,D4PI,SQFRPI,XSMALL
C
C END INCLUDE BLOCK
