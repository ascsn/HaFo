c-----testci--------------------------------------------------------------------
c
      program testci
c
c     test compex inversion
c
      include 'params.inc'
c
      complex*16 r(kg1),cumin,diff,test
      dimension  rhs(kg1)
c
      include 'potinv.inc'
      include 'grid.inc'
      include 'wavfun.inc'
      include 'constn.inc'
c
c-------------------------------------------------------------------------------
c
      one   = 1.0
      zero  = 0.0
      rstep = 0.3
      ngrid = 35
      h2m   = 20.75
      do i=1,ngrid
        rad    = rstep*(i-1)
        sqb(i)  = sqrt(h2m)
        potl(i) = 0.0
        rhs(i)  = rad*exp(-rad*rad/4.0)
        d2b(i)  = 0.0
      enddo
      cumin   = cmplx(5*one,one/2)
c
      call cinver(r,rhs,cumin)
c
      do i=1,ngrid
        wf(i)  = real(r(i))
        wf(i+ngrid)  = aimag(r(i))
      enddo
      call derv25(wf(1),wfs(1),-one,ngrid,rstep)
      call derv25(wf(1+ngrid),wfs(1+ngrid),-one,ngrid,rstep)
      do i=1,ngrid
        test = (potl(i)-cumin)*r(i)-h2m*cmplx(wfs(i),wfs(i+ngrid))
        diff = test-cmplx(rhs(i),zero)
        write(6,'(i5,f8.2,3(2g12.4,2x))')
     &   i,rstep*(i-1),test,diff,r(i)
      enddo
c
      stop
      end
c-----cinver------------ part of the hartree-fock package ----------------------
c
      subroutine cinver(r,rhs,cumin)
c
      include 'params.inc'
c
      complex*16 r(kg1)
      complex*16 cainv(kg1),cbinv(kg1),cpotef(kg1)
      complex*16 cumin,x,y,rim1,ui,z
      dimension  rhs(kg1)
c
      include 'potinv.inc'
      include 'grid.inc'
      include 'constn.inc'
c
       data fsx,ezw/.83333333333333,.8333333333333e-1/
c
c
c      determines 'r' as solution of the inhomogenous linear diff. eq.:
c       " (-d/dx (sqb(x)**2) d/dx + (potl(x) - cumin ) ) r = rhs  "
c      for given 'rhs' where 'sqb**2' is x-dependent inverse mass,
c      'potl' a potential and 'cumin' a shift constant.
c
c      the arrays in the calling list mean:
c       r      = the solution of the equation
c       rhs    = the right hand side (inhomogeneity) of the equation.
c                r and rhs may be the same arrays in the calling routine
c                rhs is then overwritten.
c      the parameter in the calling list:
c       cumin   = constant energy shift
c
c      information on the hamiltonian is transferred via common /potinv/:
c       sqb    = inverse square-root of x-dependent effective mass
c       potl   = effective potential before mass transformation,
c                to be prepared in subroutine 'welfu'.
c
c      information on the grid is transferred via common /grid/.
c
c      note: the program is written for solving the schroedinger equatio
c            in radially symmetric problems. thus it implies the rhs and
c            cpotef to be such that r(x) vanishes at the origin.
c
c      the solution proceeds as follows:
c
c      by introducing scaled wave functions
c       " reff = r/sqb "
c      and
c       " rhseff = rhs/sqb ".
c      the equations is transformed to a laplacian equation
c       " (-d/dx**2 + cpotef(x) ) reff = rhseff "
c      with constant mass and an effective potential
c       " cpotef = (potl(x)-cumin)/(sqb(x)*sqb(x)) + d2b(x) "
c      where
c       " d2b = (d/dx**2 sqb)/sqb " .
c      after solving the laplacian equation, the solution is transformed
c      back by
c       " r   = reff*sqb" .
c
c      the d/dx**2 is treated effectively in 5-point precision.
c      in addition to a standard 3-point formula for d/dx**2
c      a complementing 3-point averaging of ueff*reff and rhseff
c      is used which both together make up the 5-point precision.
c      (see r.zurmuehl, "praktische mathematik", ch. vii, sect. 3,
c       springer verlag, berlin 1963)
c
c
c-----------------------------------------------------------------------
c
c      square-root of inverse eff. mass on 'sqb',
c      second derivative of 'sqb' on 'd2b'.
c      finally assembling  the effective potential on 'cpotef'.
c       this part can be overridden if the routine is called for
c       the same 'potl' and 'sqb' again and if the workspaces
c       'sqm' and 'cpotef' are not destroyed in the meantime;
c       this situation should be indicated by 'irep>0' in the call.
c
c      if(irep.le.0) then
       sqm(1) = one/sqb(1)
        do 80 i=2,ngrid
        sqm(i)    = one/sqb(i)
        cpotef(i) = d2b(i)+(potl(i)-cumin)*(sqm(i)*sqm(i))
   80   continue
c      endif
c
c      transform to effective r.h.s., store result on r.
c
       do 17 i=1,ngrid
   17  r(i)   = rhs(i)*sqm(i)
c
c     prepare parameters for kinetic energy and loop-ends.
c
        offdig = one/(rstep*rstep)
        offdg2 = offdig+offdig
        nmaxm  = ngrid-1
c
c      check repetition switch. if new potl,sqm are indicated be irep=0
c      then begin first pass eliminating the lower off-diagonal element.
c      new diagonal and upper off-diagonal elemnts are saved on a and b.
c      the second pass, forward elimination of the r.h.s., is done withi
c      the same loop; result is stored on field t.
c
c      if(irep.le.0) then
c
        cbinv(1) = zero
        cainv(1) = one
        x      = zero
        y      = one
        rim1   = r(1)
c
        do 13 i=2,nmaxm
c                                         forward elim. for r.h.s.
        x      = ezw*(r(i+1)+rim1)+fsx*r(i)- x*(cbinv(i-1)/y)
        rim1   = r(i)
        r(i)   = x
c                                         forward elim. for matrix
        ui     = cpotef(i)
        z      = -offdig+ezw*ui
        cbinv(i) = z
        y      = offdg2+fsx*ui-z*(cbinv(i-1)/y)
        cainv(i) = y
   13   continue
c                                         last elim. step explicitely
c                                         because r(i+1) is not availabl
        r(ngrid) = ezw*rim1+fsx*r(ngrid)-x*(cbinv(nmaxm)/y)
        ui     = cpotef(ngrid)
        z      = -offdig+ezw*ui
        cbinv(ngrid) = z
        y      = offdg2+fsx*ui-z*(cbinv(nmaxm)/y)
        cainv(ngrid) = y
c
c      if unchanged potl,sqm are indicated by irep=1 then perform
c      in the second pass the forward elimination for the r.h.s.
c      using the diagonal and off-diagonal elements of the previous
c      call.
c
c       else
c        x      = zero
c        rim1   = r(1)
c        do 14 i=2,nmaxm
c        x      = ezw*(r(i+1)+rim1)+fsx*r(i)- x*cbinv(i-1)/cainv(i-1)
c        rim1   = r(i)
c        r(i)   = x
c   14   continue
c        r(ngrid)  = (ezw*rim1+fsx*r(ngrid)-x*cbinv(nmaxm)/cainv(nmaxm))
c       endif
c
c      in the third pass resolve the remaining upper triangular equation
c      backward.
c
       x      = r(ngrid)/cainv(ngrid)
       r(ngrid) = x*sqm(ngrid)
       do 15 i=nmaxm,2,(-1)
       x      = (r(i)-x*cbinv(i+1))/cainv(i)
       r(i)   = x*sqm(i)
   15  continue
       r(1)   = zero
c
       return
       end
c-----derv25------------ part of the hartree-fock package ----------------------
c
      subroutine derv25(wfin,dwfout,parity,ngrid,delr)
c
      include 'params.inc'
c
c
c     second derivative of 'wfin' on 'dwfout',
c     using five point formula.
c
c     this version uses loop unrolling in a way which inhibits
c     vectorization. a vectorizable variant is stored on the
c     lines which are initiated by 'cvec'
c
c     input via calling list:
c      wfin    = array for which 2.derivative is to be computed
c      parity  = the radial parity of the field 'wfin'
c                ( > 0 if the field starts with even power of r, < 0 else)
c      ngrid   = length of the arrays
c      delr    = distance between grid points
c     output via calling list:
c      dwfout  = array containing the 2. derivative
c
      dimension wfin(ngrid),dwfout(ngrid)
c
      data              vdr,                  ezw
     &    /1.33333333333333,0.0833333333333333333/
      data one,zero/1.0,0.0/
c
c----------------------------------------------------------------------
c
c     compute frequently occuring parameters
c
      delin2 = one/(delr*delr)
      fhhh   = 2.5*delin2
      vdrhh  = vdr*delin2
      ezwhh  = ezw*delin2
c
c     fill the first set of local wavefunctions
c
      wf1    = wfin(1)
      wf2    = wfin(2)
      wf3    = wfin(3)
      wf4    = wfin(4)
c
c     derivatives near the origin, use 'parity' to continue the wave function
c
      if(parity.gt.0.) then
       dwfout(1) = 2.0*(vdrhh*wf2-ezwhh*wf3)-fhhh*wf1
      else
       dwfout(1) = zero
      endif
      dwfout(2) = -fhhh*wf2+vdrhh*(wf3+wf1)
     *            -ezwhh*(wf4+parity*wf2)
cvec      if(parity.gt.0.) then
cvec       dwfout(1) = 2.0*(vdrhh*wfin(2)-ezwhh*wfin(3))-fhhh*wfin(1)
cvec      else
cvec       dwfout(1) = zero
cvec      endif
cvec      dwfout(2) = -fhhh*wfin(2)+vdrhh*(wfin(3)+wfin(1))
cvec     *            -ezwhh*(wfin(4)+parity*wfin(2))
c
c     the straightforward "mid-region"
c
      ngm2 = ngrid-2
      do 10 i=3,ngm2,5
      wf5 = wfin(i+2)
      dwfout(i) = -fhhh*wf3+vdrhh*(wf2+wf4)-ezwhh*(wf1+wf5)
      ia = i+1
      if(ia.gt.ngm2) goto 90
      wf1 = wfin(ia+2)
      dwfout(ia) = -fhhh*wf4+vdrhh*(wf3+wf5)-ezwhh*(wf2+wf1)
      ia = ia+1
      if(ia.gt.ngm2) goto 90
      wf2 = wfin(ia+2)
      dwfout(ia) = -fhhh*wf5+vdrhh*(wf4+wf1)-ezwhh*(wf3+wf2)
      ia = ia+1
      if(ia.gt.ngm2) goto 90
      wf3 = wfin(ia+2)
      dwfout(ia) = -fhhh*wf1+vdrhh*(wf5+wf2)-ezwhh*(wf4+wf3)
      ia = ia+1
      if(ia.gt.ngm2) goto 90
      wf4 = wfin(ia+2)
      dwfout(ia) = -fhhh*wf2+vdrhh*(wf1+wf3)-ezwhh*(wf5+wf4)
   10 continue
cvec      do 10 i=3,ngrid-2
cvec      dwfout(i) = -fhhh*wfin(i)+vdrhh*(wfin(i-1)+wfin(i+1))
cvec     *                      -ezwhh*(wfin(i-2)+wfin(i+2))
cvec   10 continue
c
c     use only 3-point formula one grid-point before boundary
c
   90 continue
cold      dwfout(ngrid-1)= (wfin(ngrid)-2.0*wfin(ngrid-1)+wfin(ngrid-2))
cold     *                *delin2
coldc
coldc     extrapolate linearly to grid boundary
coldc
cold      dwfout(ngrid) = 2.0*dwfout(ngrid-1)-dwfout(ngrid-2)
c
c     boundary condition of vanishing wave function
c
      dwfout(ngrid-1) = -fhhh*wfin(ngrid-1)
     &                  +vdrhh*(wfin(ngm2)+wfin(ngrid))
     *                  -ezwhh*(wfin(ngrid-3)    )
      dwfout(ngrid) = -fhhh*wfin(ngrid)
     &                  +vdrhh*(wfin(ngrid-1)    )
     *                  -ezwhh*(wfin(ngm2)    )
c
      return
      end

