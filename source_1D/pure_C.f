c-----pure_C------------ part of the hartree-fock package ---------------------
c
      subroutine pure_C()
c
      include 'params.inc'
c
c     puts point-Coulomb on proton potential and sets all other
c     potentials, spin-orbit potential, effective masses to zero.
c     to be used for the 'Fermi gas' calculations.
c
      include 'grid.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'wavfun.inc'
      include 'result.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------------
c
      chgfac = nprot*1.44*acoul
      do ir=1,ngrid
        pot(ir)   = chgfac/(rstep*max(ir-1,1))
        potls(ir) = 0.0
        bmass(ir) = amfc(1)
        iq = ir+ngrid
        pot(iq)   = 0.0
        potls(iq)  = 0.0
        bmass(iq) = amfc(2)
      enddo
      enrear = 0.0
c
      return
      end
