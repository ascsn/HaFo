c-----lapfil------------------------------------------------------------
c
       subroutine lapfil(rfilt,awork,xmas2,irep)
c
      include 'params.inc'
c
       dimension rfilt(kg1),awork(kg1)
c
       include 'grid.inc'
c
       data fsx,ezw/.83333333,.83333333e-1/
ctest       data fsx,ezw/1.0,0.0/
c
c      the "laplace-folding-filter" for the skyrme-hartree-fock
c      package.
c      solves the inhomogenous linear diff. eq.:
c       " xmas2^-1 * (-d/dr**2 - 2.0*d/dr + xmas2) rfield = rfield  "
c
c      the d/dx**2 is treated effectively in 5-point precision.
c      in addition to a standard 3-point formula for d/dx**2
c      a complementing 3-point averaging of xmas2*reff and rhseff
c      is used which both together make up the 5-point precision.
c
c
c      the arrays in the call mean:
c       rfilt  = field to be filtered at input and result at output.
c       awork  = work space to keep the off-diagonal element for inversion.
c      the parameters in the call mean:
c       xmas2  = parameter in the laplacian equation, means mass**2 in
c                case of meson folding
c       rstep  = stepsize for the grid
c       ngrid   = maximum number of grid points
c       irep   = switch: irep=0 means new potl,sqb such that full
c                        inversion has to be performed from scratch.
c                        irep=1 means that previous potl,sqb still apply
c                        than parts of the first pass can be taken from
c                        the previous calls, whose results have been
c                        saved on the fields a and b.
c
c-----------------------------------------------------------------------
c
c      transform to cartesian volume element. 'rfilt' is overwritten!
c      and save the normalization.
c
       xnorm = 0.0
       do ir=1,ngrid
         xnorm = rsquar(ir)*rfilt(ir) + xnorm
         rfilt(ir) = rfact(ir)*rfilt(ir)
       enddo
       xnorm = frpirs*xnorm
c
c      check repetition switch. if new xmas2 is indicated be irep=0
c      then begin first pass eliminating the lower off-diagonal element.
c      new diagonal and upper off-diagonal elemnts are saved on a and b.
c      the second pass, forward elimination of the r.h.s., is done withi
c      the same loop; result is stored on field t.
c
c    prepare some parameters for kinetic energy and loop-ends.
c
c       the original on- and off-diagonal elements
c
        rinv2  = 1./(rstep*rstep)
        offdig = -rinv2+ezw*xmas2
        ondiag = rinv2+rinv2+fsx*xmas2
c
c       rescaling by 1/offdig acts on diagonal element and on r.h.s.
c
        offdgi = 1.0/offdig
        ondgrs = ondiag*offdgi
        ezwrs  = ezw*offdgi
        fsxrs  = fsx*offdgi
        nmaxm  = ngrid-1
c
       if(irep.le.0) then
c                                 if(irep.le.0) then
c
        awork(1) = 1.0
        x      = 0.
        y      = 1.0
c
c                                         i=2 explicitely
        x      = ezwrs*(rfilt(3)+rfilt(1))+fsxrs*rfilt(2)
        rfsav  = rfilt(2)
        rfilt(2)   = x
        y      = 1.0/ondgrs
        awork(2)   = y
c                                         now the most other i
        do 13 i=3,nmaxm
c                                         forward elim. for r.h.s.
        x      = ezwrs*(rfilt(i+1)+rfsav)+fsxrs*rfilt(i)- x*y
        rfsav  = rfilt(i)
        rfilt(i)   = x
c                                         forward elim. for matrix
        y      = 1.0/(ondgrs-y)
        awork(i) = y
   13   continue
c                                         last elim. step explicitely
c                                         because rfilt(i+1) is not availabl
        rfilt(ngrid) = ezwrs*rfsav+fsxrs*rfilt(ngrid)-x*y
cold        y      = 1.0/(ondiag-y)
        y      = 1.0/(ondgrs-y)
        awork(ngrid) = y
c
c      if unchanged xmas2 is indicated by irep=1 then perform
c      in the second pass the forward elimination for the r.h.s.
c      using the diagonal and off-diagonal elements of the previous
c      call.
c
       else
        x      = 0.
        rfsav  = rfilt(1)
        do 14 i=2,nmaxm
        x      = ezwrs*(rfilt(i+1)+rfsav)+fsxrs*rfilt(i)- x*awork(i-1)
        rfsav  = rfilt(i)
        rfilt(i)   = x
   14   continue
        rfilt(ngrid)  = (ezwrs*rfsav+fsxrs*rfilt(ngrid)-x*awork(nmaxm))
       endif
c
c      in the third pass resolve the remaining upper triangular equation
c      backward.
c
       x      = rfilt(ngrid)*awork(ngrid)
       rfilt(ngrid) = x
       do 15 i=nmaxm,2,(-1)
       x      = (rfilt(i)-x)*awork(i)
       rfilt(i)   = x
   15  continue
c
c      final division by r and normalization
c
       sum = 0.0
       do ir=2,ngrid
         rfilt(ir) = xmas2*rinvs(ir)*rfilt(ir)
         sum = rsquar(ir)*rfilt(ir) + sum
       enddo
       xnorm = xnorm/(frpirs*sum)
       do ir=2,ngrid
         rfilt(ir) = xnorm*rfilt(ir)
       enddo
       rfilt(1)=(4.0*rfilt(2)-rfilt(3))*0.33333333
c7       rfilt(1) = 1.5*rfilt(2)-0.6*rfilt(3)+0.1*rfilt(4)
c
       return
       end
