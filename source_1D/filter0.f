c-----filter0-------------------------------------------------------------------
c
      subroutine filter0(fldin,fldout)
c
      include 'params.inc'
c
c     folding operation using Fourier-Bessel transformation for L=0
c     input on 'fldin'. output on 'fldout'. ('fldin' and 'fldout'
c     can be the same fields).
c     the workspaces 'aux1' and 'aux2' are needed.
c     The folding function does enhance weights of high q.
c
      dimension fldin(kg1),fldout(kg1)
      dimension auxfou(kgr2)
c
      include 'wavfun.inc'
      include 'linres.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
      equivalence (aux1(1),auxfou(1))
c
      antifold(q) = 1.0D0+0.3D0*TANH((q-2.7D0)/0.27D0)
c
c-------------------------------------------------------------------
c
c
c     formfactor of 'fldin',
c     augmented by folding weight,
c     and stored on workspace 'auxfou'.
c
      do iq=1,ngrid
        q1     = iq*delq-delq
        sum    = zero
        dc     = cos(q1*rstep)
        ds     = sin(q1*rstep)
        r      = 0.0D0
        c      = one
        s      = zero
        do ir=2,ngrid
          r      = rstep+r
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
          sum    = r*s*fldin(i)+sum
        enddo
        auxfou(iq) = rstep*sum  ! * antifold(q1)
      enddo
c
c     fourier-bessel backtransformation
c
      r      = zero
      do ir=2,ngrid
        r      = rstep+r
c
        dc     = cos(r*delq)
        ds     = sin(r*delq)
        c      = one
        s      = zero
        sum    = zero
        qa     = zero
        do iq=2,iqmax
          qa     = delq+qa
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
          sum    = qa*s*auxfou(iq)+sum
        enddo
        fldout(ir) = (sum*delq*dpi3/r)
c
      enddo
c
c
      return
      end
