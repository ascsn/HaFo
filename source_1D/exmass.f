c-----exfunc------------ part of the hartree-fock package --------------------
c
      double precision function exfunc(rhoact)
c
c     the explicit exchange functional. Input is
c      'rhoact'  via list
c      'exmass'  via  common /skypar/
c      'coefx'   via  common /skypar/
c
      include 'params.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------
c
c
c     reduced Fermi-momentum 'eta' employs 29.6088132=3*pi**2
c
      eta    = (29.6088132*rhoact)**0.333333333d0*2.0/exmass
      eta2   = eta*eta
      if(ioptsk.ne.3) then
        exfunc = coefx/24.0 * (
     &           1.5d0*eta2*eta2-eta2
     &           -4.0*eta2*eta*atan(eta)
     &           +(3.0d0*eta2+1.0d0)*log(1.0d0+eta2)      )
      else
        exfunc = coefx*eta**6/(one+0.11111111111*eta2)
      endif
c
      return
      end
c-----expot ------------ part of the hartree-fock package --------------------
c
      double precision function expot(rhoact)
c
c     the explicit mean field potental from the exchange functional. Input is
c      'rhoact'  via list
c      'exmass'  via  common /skypar/
c      'coefx'   via  common /skypar/
c
      include 'params.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------
c
c
c     reduced Fermi-momentum 'eta' employs 29.6088132=3*pi**2
c
      eta    = (29.6088132*rhoact)**0.333333333d0*2.0/exmass
      eta2   = eta*eta
      if(ioptsk.ne.3) then
        expot  = coefx*eta/(3.0d0*rhoact) * (
     &           0.5d0*eta2*eta
     &           -eta2*atan(eta)
     &           +0.5d0*eta*log(1.0d0+eta2)      )     * 0.5
      else
        expot  = coefx*eta/(3.0d0*rhoact)
     &          *eta**5*(6.0+0.444444444444*eta2)
     &          /(one+0.11111111111*eta2)**2
      endif
c
      return
      end
c-----exenre------------ part of the hartree-fock package --------------------
c
      double precision function exenre(rhoact)
c
c     the rearrangement energy from the exchange functional. Input is
c      'rhoact'  via list
c      'exmass'  via  common /skypar/
c      'coefx'   via  common /skypar/
c
      include 'params.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------
c
c
c     reduced Fermi-momentum 'eta' employs 29.6088132=3*pi**2
c
      eta    = (29.6088132*rhoact)**0.333333333d0*2.0/exmass
      eta2   = eta*eta
      atane  = eta2*atan(eta)
      xloge  = log(1.0d0+eta2)
      if(ioptsk.ne.3) then
        expot  = coefx*eta/(3.0d0*rhoact) * (
     &           0.5d0*eta2*eta-atane+0.5d0*eta*xloge      )    *   0.5
        exfunc = coefx/24.0 * (
     &           1.5d0*eta2*eta2-eta2-4.0*eta*atane
     &           +(3.0d0*eta2+1.0d0)*xloge      )
        exenre = exfunc - 0.5d0*rhoact*expot
      else
        exenre = coefx/27.0 * eta**8 /(one+0.11111111111*eta2)**2
      endif
c
      return
      end
c-----dirfun------------ part of the hartree-fock package --------------------
c
      double precision function dirfun(rhoact)
c
c     the direct-part functional in Pade approximation. Input is
c      'rhoact'  via list
c      'b3dir'  via  common /skypar/
c      'power'   via  common /skypar/
c
      include 'params.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------
c
c
cold      dirfun = 0.5*b0dir*rhoact*rhoact/(one+gamdir*rhoact)
      dirfun = b3dir*rhoact**3/(one+power*rhoact)
c
      return
      end
c-----dirpot------------ part of the hartree-fock package --------------------
c
      double precision function dirpot(rhoact)
c
c     the m.f. potential for direct-part functional in Pade approximation. 
c     Input is
c      'rhoact'  via list
c      'b3dir'  via  common /skypar/
c      'power'   via  common /skypar/
c
      include 'params.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------
c
c
cold      dirpot = 0.5*b0dir*rhoact*(2.0+gamdir*rhoact)
cold     &         /(one+gamdir*rhoact)**2
      dirpot = b3dir*rhoact*rhoact*(3.0+2.0*power*rhoact)
     &         /(one+power*rhoact)**2
c      write(2,'(a,2g12.4)') ' DIRPOT: b0dir,gamdir=',b0dir,gamdir
c
      return
      end
c-----direnr------------ part of the hartree-fock package --------------------
c
      double precision function direnr(rhoact)
c
c     the rearr. energy for direct-part functional in Pade approximation. 
c     Input is
c      'rhoact'  via list
c      'b3dir'  via  common /skypar/
c      'power'   via  common /skypar/
c
      include 'params.inc'
      include 'skypar.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------
c
c
      direnr = -0.5*b3dir*rhoact**3/(one+power*rhoact)**2
c
      return
      end









