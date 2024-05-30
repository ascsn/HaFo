c-----lanpar------- part of the Hartree-Fock package ------------------
c
      subroutine lanpar(rhonm0,flan0,flan0p,glan0,glan0p,
     &                  f0,f0p,g0,g0p)
c
c     evaluates the Landau parameters for symmetric nulcear matter
c     at given density 'rhonm0'.
c
c     the force parameters are communicated via common/skypar/.
c
c     the Landau parameters are 
c       glan0   = landau paramer    g0
c       glan0p  = landau paramer    g0'
c       flan0   = landau paramer    f0
c       flan0p  = landau paramer    f0'
c     their dimensionless analogues are 'f0','f0p','g0', and 'g0p'
c
      include 'params.inc'
      include 'skypar.inc'
c      include 'numpar.inc'
c      include 'grid.inc'
c      include 'wavfun.inc'
c     include 'constn.inc'
c
      data pi/3.141592654/
      data zero/0.0/
cold      data rhoini,drhoit,epsil,itrmax
cold     *    / 0.16 , 0.01 ,1.e-5,    10/
      data rhoini,drhoit,itrmax
     *    / 0.16 , 0.01 ,    20/
c
c-----------------------------------------------------------------------
c
c     preparation of recoupled parameters
c     t12p = cm(krewald), t45p = cm1(krewald), diss. p.64, gl.(22,23)
c
      powp1  = power+1.0
      betab  = 3./5.*(3.*pi*pi/2.)**(2./3.)
      hqqd2m = 0.5*(amfc(1)+amfc(2))
      if(ifzpe.ne.1) then
        hqqd2m = hqqd2m/(1.0-1.0/anumb)
      endif
c     hqqd2m = 197.3**2/2./936.
      t12p   = 3.*t1+5.*t2*(1.+4./5.*x2)
      t45p   = (3.*t5+5.*t6*(1.+4./5.*x6))/6.0
      fermp2 = (3.*pi*pi/2.*rhonm0)**(0.66666667)
      if(rhoswm.gt.zero) then
        tanfac = tanh(rhonm/rhoswm)
      else
        tanfac = 1.0
      endif
      effmas = 1./(1.+rhonm0/16./hqqd2m*(t12p*tanfac+t45p*rhonm0
     &               +(b6+0.5*b6p)*0.6*fermp2*rhonm0))
c
      hvonti = 0.5*t0*x0
     &        +1./12.*t3*x3*rhonm0**power
     &        +0.25*t1*x1*fermp2
     &        +0.25*(t5*x5/6.0)*rhonm0*fermp2
      flan0  = 3./4.*t0
     &        +t3/16.0*powp1*(power+2.0)*rhonm0**power
     &        +1./8.*fermp2*(t12p+9.*rhonm0*t5/6.0)
      glan0p =-t0/4.-rhonm0**power*t3/24.
     &        -fermp2/8.*(t1-t2+rhonm0*t5/6.0)
      flan0p = glan0p+t2*x2*fermp2/4.-hvonti
      glan0  = glan0p+t2*x2*fermp2/4.+hvonti
c
      wfac   = sqrt(fermp2)*effmas/(pi*pi*hqqd2m)
      f0     = wfac*flan0
      f0p    = wfac*flan0p
      g0     = wfac*glan0
      g0p    = wfac*glan0p
c
      return
      end

