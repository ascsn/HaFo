c-----k2psi  ------------ part of the hartree-fock package ---------------------
c
      subroutine k2psi(wfin,hwfout,d2wfin,jsa)
c
!      include 'params.inc'
      IMPLICIT REAL(8) (a-h,o-z)
c
c     action of the operator -nabla**2 
c
c     where:
c       wfin    = input wave function field ,
c       hwfout  = output wave function field ,
c       jsa     = combined angular momentum q.n. of wfin and hwfout,
c
c     grid properties are transferred via common /grid/
c
      dimension wfin(*),hwfout(*),d2wfin(*)
c
      include 'grid.inc'
!      include 'skypar.inc'
!      include 'constn.inc'
c
c-------------------------------------------------------------------------------
c
c     the quantum numbers of the wave function
c      la      = orbital angular momentum of state  wfin>
c      parity  = parity for the radial wave function (note that it
c                contains an extra factor 'r' such that the 1s1/2
c                state has parity=-1 etc)

c
      la     = jsa/2
      xlp    = (la*la+la)
      parity = (mod(la,2)*2-1)
cold      isoadd = (iso-1)*ngrid
c
c     the derivatives of the wave function
c
      call derv25(wfin,d2wfin,parity,ngrid,rstep)
c
c     action of -nabla**2
c
      hwfout(1) = d2wfin(1)
      DO ir=2,ngrid
        hwfout(ir) = xlp*wfin(ir)/((ir-1)*rstep)**2 - d2wfin(ir)
      END DO
c
      return
      end


c-----derv25------------ part of the hartree-fock package ----------------------
c
      subroutine derv25(wfin,dwfout,parity,ngrid,delr)
c
!      include 'params.inc'
      IMPLICIT REAL(8) (a-h,o-z)
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
