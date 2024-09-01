c-----hoteos------- part of the Hartree-Fock package ------------------
c
      subroutine hoteos(temper)
c
c     evaluates E/A and pessure for symmetric nuclear matter.

c     the force parameters are communicated via common/skypar/.
c
c     list arameters are:
c       temper  = input of temperature
c
      include 'params.inc'
      include 'skypar.inc'
      include 'numpar.inc'
      include 'grid.inc'
      include 'wavfun.inc'
      include 'result.inc'
      include 'constn.inc'
c
c      data pi/3.141592654/
c      data zero/0.0/
      data rhoini,rhomax,drhoit
     *    / 0.0, 0.3, 0.01/
      data deltpc/1.0D-3/            ! variation of proton content
c
c-----------------------------------------------------------------------
c
      itrmax = (rhomax-rhoini)/drhoit
      drhoit = (rhomax-rhoini)/max(itrmax,1)
      write(12,'(a,f6.2/a)') 
     &  '#  EoS at temperature =',temper,
     &  '#  rho  E/A   pressure '
      do irho=1,itrmax
        rhonm  = rhoini + irho*drhoit
c
        call enernm(rhonm,half,temper,eovera,efmasp,efmasn)
        fovera = eovera - 2*temper*entropy/rhonm   ! -efermh
        call enernm(rhonm+deltpc,half,temper,eplus,efmasp,efmasn)
        fplus  = eplus - 2*temper*entropy/rhonm   ! -efermh
        call enernm(rhonm-deltpc,half,temper,eminus,efmasp,efmasn)
        fminus = eminus - 2*temper*entropy/rhonm   ! -efermh
        press = rhonm*rhonm*(eplus-eminus)/(deltpc+deltpc)
        pressf = rhonm*rhonm*(fplus-fminus)/(deltpc+deltpc)
        write(12,'(1x,6(1pg13.5))')   
     &    rhonm,eovera,press,fovera,pressf,efermh
c
      enddo
      write(12,'(a)') '     ','     '
      return
      end

