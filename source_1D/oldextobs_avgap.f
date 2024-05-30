c-----extobsiso--------------------------------------------------------
c
      subroutine extobs_avgap(par,workdv,nparmx,npr,nne,isodelt)
c
c     Prepare output for computing correlations with average gaps
c       npr  proton number
c       nne  neutron number
c       isodelt isospin for gap
c
      include 'params.inc'
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'comfit.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'extobs.inc'
      include 'constn.inc'
c
c     'workdv' is field for derivatives with respect to parameter-variation.
c
      dimension   par(kparmx),pars(kparmx),workdv(kparmx,1)
      real*8  work2(kparmx,kparmx) !,vect(kparmx,kparmx),eig(kparmx)
      logical tsymderiv
      parameter (tsymderiv=.true.)
      logical tsecord,tstsecord
      parameter (tsecord=.false.,tstsecord=.false.)
      logical t182
      character*1 asec
      character*3 id2
c
c
c-----------------------------------------------------------------------
c
      INQUIRE(83,opened=t182)
      IF(.NOT.t182) THEN     
        open(39,status='scratch')
        write(39,'(i4)') nint(parmnatfac*100+1000)
        rewind(39)
        read(39,'(1x,a3)') id2
        close(39)
        OPEN(unit=182,status='unknown',file='error-avgap.'//id2)
      END IF
c
c     loop over variation of parameters
c
      IF(tsymderiv) THEN
        nparlow = -nparmx
      ELSE
        nparlow = 0
      END IF
      pars = par
      do i=nparlow,nparmx
      IF(tsecord) THEN
        nparm2=abs(i)
        nparlw2=-abs(i)
      ELSE
        nparm2=0
        nparlw2=0
      END IF
      do j=nparlw2,nparm2
        if(j.ne.0 .and. i.eq.0) CYCLE
        if(j.ne.0 .and. j.eq.-i) CYCLE
        if(j.eq.0 .and. i.eq.0) par(1:nparmx) = pars(1:nparmx)
        if(j.eq.0 .and. i.ne.0) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        if(j.ne.0 .and. j.ne.i) par(1:nparmx) = pars(1:nparmx)
     &             + SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))/2D0
     &             + SIGN(1,j)*parmnatfac*parmnat(1:nparmx,abs(j))/2D0
        if(j.ne.0 .and. j.eq.i) par(1:nparmx) = pars(1:nparmx)
     &             - SIGN(1,i)*parmnatfac*parmnat(1:nparmx,abs(i))
        call attach(par,jp)
c
        nprot  = npr
        nneut  = nne
        nnumb  = nprot+nneut
        anumb  = nnumb
        nmaxp  = 0
        nmaxn  = 0
        call hafosu
        IF(isodelt==1) THEN
          observ = deltp
        ELSE
          observ = deltn
        END IF
c
c       collect array of deviations
c
        CALL  collectw(i,j,observ,observsav,workdv,work2,tsymderiv)
c
        par(1:nparmx) = pars(1:nparmx)
      enddo
      enddo

c
c     accumulate extrapolation errors by folding with 'zarray'.
c
      errobs   = 0.0
c
      if(nparmx.gt.0) then
        do i=1,nparmx
          errobs  = workdv(i,1)*workdv(i,1) + errobs
        enddo
        errobs  = sqrt(errobs)/parmnatfac
      endif
c
c     final print
c
      IF(tsecord) THEN
        asec='#'
      ELSE
        asec=' '
      END IF
      IF(isodelt==1) THEN
        write(02,'(a,2f9.4)')
     &  'proton gap: value,uncertainty=',observsav,errobs
        write(182,'(a)') asec//' proton gap'
      ELSE
        write(02,'(a,2f9.4)')
     &  'neutron gap: value,uncertainty=',observsav,errobs
        write(182,'(a)') asec//' neutron gap'
      END IF
        write(182,'(1pg24.16)')  observsav
        write(182,'(100(1pg24.16))') workdv(1:nparmx,1)
        IF(tsecord) 
     &     write(182,'(1000(1pg24.16))') work2(1:nparmx,1:nparmx)
      CALL flush(2)
      CALL flush(182)
c
c     estimate curvature contributions
c
      IF(tsecord) CALL estim_curv(workdv,work2,nparmx,tstsecord)
c
      return
      end
