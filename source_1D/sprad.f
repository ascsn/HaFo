c-----sprad ------------ part of the hartree-fock package ---------------------
c
      subroutine sprad(rmssta,iso)
c
      include 'params.inc'
c
c     computes rms radii of separate s.p. states for type 'iso',
c     result is returned on 'rmssta'.
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
       if(iso.eq.1) then
         nmax = nmaxp
         nadd = 0
       else
         nmax = nmaxn
         nadd = nmaxp
       endif
       do nn=1,nmax
         n = nn + nadd
         js=jsp(nn)
         ll=js/2
         if(2*ll.eq.js) js=js-1
         npl = (n-1)*ngrid
         sum0 = 0.0
         sum2 = 0.0
         do ir=1,ngrid
           r  = (ir-1)*rstep
           r2 = r*r
           vol= wf(npl+ir)*wf(npl+ir)*r2
           sum0 = vol + sum0
           sum2 = vol*r2 + sum2
         enddo
         rmssta(nn) = sqrt(abs(sum2/max(sum0,xsmall)))
      enddo
c
      return
      end


