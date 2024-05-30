c-----gaufol------------ part of the hartree-fock package ----------------------
c
      subroutine gaufol(rhoin,rhoout,width)
c
      include 'params.inc'
c
c     Gaussian folding using Fourier-Bessel transformation:
c      rhoin   = input density
c      rhoout  = folded output density
c      width   = width of folding (r.m.s. of Gaussian)
c     The grid parameters are transferred via commons.
c
      dimension rhoin(kg1),rhoout(kg1)
c
      include 'grid.inc'
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'constn.inc'
c
c
c----------------------------------------------------------------
c
c     preparation of auxiliary parameters
c
      anumb  = one*(nprot+nneut)
      iqmax  = ngrid-1
      rleng  = (ngrid-1)*rstep
      delq   = pi/rleng
      frpirs = 4*pi*rstep
      xlambd = 0.5/(width*width)
      qlambd = 0.25/xlambd
c
c     fourier-bessel transformation forward
c
      sum    = zero
      r      = zero
      do ir=2,ngrid
        r     = rstep+r
        sum   = r*r*rhoin(ir)+sum
      enddo
      aux1(1) = frpirs*sum
c
      do iq=2,iqmax
        q1     = iq*delq-delq
        dc     = cos(q1*rstep)
        ds     = sin(q1*rstep)
        r      = zero
        sum   = zero
        c      = one
        s      = zero
        do i=2,ngrid
          r      = rstep+r
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
          sum    = r*s*rhop(i)+sum
        enddo
        sum        = frpirs*sum/q1
        aux1(iq)   = exp(-qlambd*(q1*q1))*sum
      enddo
c
c     fourier-bessel transformation backward
c
      dpi3   = half/(pi*pi)
c                                                         r=0.0 separately
      sum    = zero
      do iq=2,iqmax
        qa     = delq*iq-delq
        sum    = qa*qa*aux1(iq)+sum
      enddo
      rhoout(1) = abs(sum*delq*dpi3)
c                                                         now all other r
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
          sum    = qa*s*aux1(iq)+sum
        enddo
        rhoout(ir) = abs(sum*delq*dpi3/r)
c
      enddo
c
      return
      end
