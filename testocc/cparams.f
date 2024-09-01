c-----cparams--------------------------------------------------------------------
c
      subroutine cparams()
c
c     compute the 'C' parameters according to the appendic of the RMP
c     take input  t0,t1,... from common /skypar/
c     result is printed on 'for002'
c
      include 'params.inc'
c
      include 'skypar.inc'
      include 'comfit.inc'
c
c--------------------------------------------------------------------
c
      crho0 = 3.0/8.0*t0+3.0/48.0*t3
