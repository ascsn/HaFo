c-----hffit ----------------------------------------------------------
c
       program hffit
c
c    !! modification: array 'chi' in calling list of 'xd1' and 'ffverg'.
c
c   Routine for fitting Skyrme-Hartree-Fock model parameters to
c   empirical data on ground state properties and differences
c   thereof.
c
c   Started: 28.12.1982
c
c
c    the input parameters read in from file 01 are:
c     isteur = switches reading of parameters.
c              isteur < 0  --> stop
c              isteur > 0  --> read and run
c              isteur >10  -->  "    "   " , include 'korrektur' in
c                              the input data
c              recommendation: use only isteur=0 or isteur=11
c              together with the input card sequence below
c     commen = chance to put a comment which is printed as a header
c     nwdr   = maximum number of iteration steps in the fit
c              nwdr=0 just evaluates the chi**2 in one path
c              typical values for a fit are between 8 and 20
c              (then look at the results and restart if needed)
c nwdr,ityp,nparot,nwdro,chitmp
c     ityp   = method for chi**2 minimisation
c               0 --> standard Bevington
c               1 --> simplex
c               2 --> Powell
c               3 --> interlaced Bevington 
c                     (dividing into outer and inner parameters,
c                      inner parameters are fully optimised for each
c                      given outer set,
c                      outer parameters are varied to provide gradient
c                      and full curvature matrix,
c                      the next set is determined by estimated minimum
c                      from gradient and curvature)
c                4 --> mix of Monte-Carle and Bevington
c                      (inner paremeters are again optmised with Bevington,
c                       outer parameters are driven by Monte Carlo)
c     nparot  = number of outer parameters  (in case of ityp=3 or 4)
c                (number of inner parameters is again given by 'nparin')
c                outer parameters follow immediately after inner parameters
c                in input string.
c     nwdro   = number of out iteration steps (in case of ityp=3 or 4)
c     chitmp  = "temperature" for Monte-Carlo decisions (case ityp=4),
c               dimensionless, determines weight in exponential
c                   exp(-delta_chi**2/chitmp) .
c
c     ifit1 = 'choice of observables', is an integer field ifit1(12).
c              allows to select the observables to be fitted by
c              inserting their key values
c                 1 <--> binding-energy
c                            bind.-energy/particle for nuc.matter
c                 2 <--> diffraction radius
c                          density in case of nuclear matter
c                 3 <--> surface thickness
c                          incompressibility in case of nuclear matter
c                 4 <--> energy/particle for the equation of state
c                        of symmetric nuclear matter
c                        or the asymmetry energy at the ground state
c                        of nuclear matter
c                 5 <--> proton  spin-orbit splitting in 16o
c                        Landau parameter g0 in case of nuclear matter
c                 6 <--> neutron spin-orbit splitting in 16o
c                        Landau parameter g0' in case of nuclear matter
c                 7 <--> (rho(0.fm)+rho(1.8fm))/2.
cold                 8 <--> (rho(0.fm)-rho(1.8fm))*100. 'shell fluct.'
c                 8 <--> rms neutron skin
c                        or
c                        d/drho a_4 in case of nuclear matter
c                 9 <--> rho(3.6fm)
c                        or average proton gap in non magic nuclei
c                10 <--> rho(5.4fm)
c                        or average neutron gap in non magic nuclei
c                        or strength of strange quark in case of n.m.
c                11 <--> rho(7.2fm)     or
c                        sumrule 'kappa' in case of nuclear matter
c                12 <--> the effective mass in case of nuclear matter
c              in any desired order. collection of observables stops
c              as soon as the first zero is detected.
c              e.g.  ifit1=2,1,3,0,0,0,0,...  means that radius,
c              energy and surfaces are included in the fit, for a
c              selection of nuclei read in from file 10.
c              note that some key definitions have been overlayed
c              in order to allow for simultanous fit of some
c              nuclear matter properties. nuclear matter is denoted
c              by reading nra=999 from file 10; the nrz=0 then means
c              the ground state of nuclear matter having the
c              'observables' e/a, density, compressibility and
c              effective mass, whereas nrz.ne.0 specifies the fermi
c              momentum by fer.mom.=nrz/100 and the 'observable' is
c              then the according e/a.
c     fehler = 'error in data', is the adopted error for the observables
c              specified in the line above (ifit1), to be placed in the
c              order as the observables in ifit1. this error is diffe-
c              reht from the experimental error; it merely adjusts the
c              relative weight of the observables in the fit in the
c              expression for the chi**2
c              in case of nuclear matter fit 'fehler' is multiplied
c              with a correction factor on 'dyy' which allows to
c              vary the error weight with each fermi-momentum.
c     fehlpr = ' % or absolut', determines whether the absolut error
c              or the relative error (in %) is considered in ch**2.
c                0.0 --> 'fehler' is to be relative error in chi**2:
c                 chi**2 = sum<((exp.obs.-calc.obs.)/erreff)**2!
c                          with erreff= (error/100.)*exp.obs
c                1.0 --> 'fehler' is to be absolute error in chi**2:
c                 chi**2 = sum<((exp.obs.-calc.obs.)/error)**2!
c     faksub = ' correction ', allows to change the values of the
c              observables by the amount given. this is for testing
c              purposes and is better left unused in actual runs.
c              just put a 0.0 in each of the 12 entries (corresponding
c              to the 12 possible observables in ifit1).
c     iaunt  = lower limit for for the mass number window for the
c              nuclei to be included in the fit ( out of the set
c              for which data is read in from file 10)
c     iaob   = upper limit for for the mass number window
c     maxidf = max. nr. of isotopic differences
c              maxdif >0 then read maxdif sequences of
c                 nuc1a,nuc1z,nuc2a,nuc2z,iobstp,errdif(idif)
c                 nuc1a   = mass of first nucleus
c                 nuc1z   = charge of first nucleus
c                 nuc2a   = mass of second nucleus
c                 nuc2z   = charge of second nucleus
c                 iobstp  = type of observable in numbering as explained
c                           above
c                 errdif  = error on this isotopic difference
c     if3gap = read nuclei for three-point differences from 'for009'
c     nminv  = key-value (out of 1,...,12 from the above list for
c              ifits1). all nuclei out of the selection provided in
c              file 10 for which key-value is missing experimental
c              data are excluded from the fit. this option is dis-
c              abled if nminv=0 (recommended value).
c     acoul  = factor which is multiplying the coulomb force.
c              typical values are
c                 1.0 --> normal coulomb force switched on
c                 0.0 --> coulomb force switched off
c     baryms(1) = proton  mass in mev
c     baryms(2) = neutron mass in mev
c     ifenrg = cleaning up energy in hartree caclculation.
c                 1 --> one path more in hartree calculation after
c                       iteration has finished (recommended for fit)
c                 0 --> just takes the last current value of energy
c                       after iteration has finished
c     ifrso  = switches evaluation of tensor current for crho.
c                 1 --> tensor current contributes to crho
c                 0 --> tensor current left out of crho
c              yet to be implemented ( only 0 is working )
c     ifform = switches evaluation of charge-formfactor
c                 0 --> does the minimum needed to determine
c                       radius and surface
c                 1 --> does full evaluation (recommended)
c     ifzpe  = switches zero-point  correction to the energy
c                -1 --> correction estimated by fit to skyrme h.f.
c                 0 --> correction estimated by harm.osc. shell model
c                 1 --> subtracts the spurious kinetic energy of the
c                       center-of-mass (recommended)
c     ifprpr = number of pairing shells above the minimum needed.
c              abs(ifprpr) switches:
c                 0 --> pairing just up to the valence shell
c                 1 --> one more shell above the valence shell
c                 2 --> two more shells above the valence shell
c              sign(ifprpr) switches:
c               <=0 --> constant gap pairing
c               >0  --> constant force pairing
c     iwrte  = print level in the chi**2-subroutine 'ffverg'
c     ist1   = print level in the fitting-subroutine 'xd1'
c                 0 --> minimum necessary information printed
c                 1 --> some more intermediate quantities are
c                       printed in addition ( for checking the
c                       fitting path )
c     iextan = switch for final analysis of other nuclei
c     nparin = number of starting model parameters to be read in below.
c              must be 'kparmx' in the present code.
c     numax  = number of parameters to be varied in the fit.
c              the first numax parameters out of the list read in below
c              will be varied in the fit. numax must be greater than 0.
c     freep  = determines length of the trial steps in the chi**2
c              surfaces to get the gradient and curvature. use values
c              between 1. and 4. ( 1. is in principle more precise,
c              4. is more stable in presence of round-off errors).
c     iftinp = switch for input of skyrme parameters
c               1 --> standard form t0,t1,...
c               0 --> recoupled form b0,b0p,...
c     ib2fix = switch for "automatic" b2                             !new
c               0 --> b2 as given in for001                          !new
c               1 --> b2 = b1/2 and b2p = 1p/2 fixed                !new
c     ioptsk = switch for variants of the Skyrme functional
c               0 --> standard
c               1 --> b5 as input and density dependent mass (surface?)
c               2 --> rational approximant for isoscalar density dep. term
c               3 --> rational approximant for isoscalar and isovector term
c                     (isovector requires "exmass" = parameter(14) as input)
c               4 --> RMF like kinetic&ls term (use 'b1kin' from input 'b1'
c                     and 'b5kin' as input from parameter(14))
c               5 --> rho-dependent kinetic and l*s, direct version
c               6 --> rho-dependent kinetic and l*s, inverse version
c               7 --> BHF functional for rho-dependent terms, needs B1,B1P=0
c               8 --> alternative isovector b3 and dens.dep. gradient
c               9 --> extra density dependent term with own 'power'
c              10 --> flexible tensor couplings
c                     use this option also if Fayans pairing is wanted
c              11 --> accumulates the options 5, 8, 9, and 10
c              12 --> QMC model, requires IFTINP=0
c
c     now come 'kparmx' lines to be read which contain the initial values
c     of the force parameters, the length of the first trial step,
c     and the key-number, in the order
c        parameter, trial-step, key
c     ( changing the initial values may change the final fit. one
c       should experiment a little bit to make sure that one has
c       found a minimum ).
c     place the parameters which are going to be varied on top of
c     the list ( the first numax parameters are varied, the others
c     are kept fixed ).
c
c
c***********   sample input file, unit 01  **************************
c
c 011,'  walecka-model   ' ,isteur,commen
c 03,        'nwdr'
c 2,1,3,0,0,0,0,0,0,0,0,0,0,0,0,0,        'auswahl der messgroessen'
c 2.0,1.0,5.0,0.0,0.0,0.,0.,0.,0.,0.,0.,0.,0.,0., 'fehler an daten'
c 00.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,  ' % oder absolut
c 00.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., ' korrrektur '
c 010,220,0,        'iaunt,iaob,maxdif,'
c 0,0,0,1,1,1,0,  'nminv,if:ex,tm,enrg,rso,zpe,tond'
c 1.0,,,,         'acoul,addnew,dabe,dabr '
c 0,0,1,,         'iwrte,ist1,ifcurv,itmax'
c 15,4,0,1.,,     'nparin,numax,ivar,freep,endcon'
c -1128.75,1.000,1, t0
c +395.   ,1.000,2, t1
c -095.   ,0.200,3, t2
c +14000. ,0.100,4, t3
c +120.   ,0.010,5, t4
c +0.45   ,0.001,6, x0
c +0.0    ,1.000,7, x1
c +0.0    ,1.000,8, x2
c +1.0    ,1.000,9, x3
c +1.0    ,1.000,10, power
c +0.0    ,1.000,11, tyuk
c +0.0    ,1.000,12, xyuk
c +0.5    ,1.000,13, ayuk
c +000.000,1.000,14,  n.n.
c +000.000,1.000,15,  n.n.
c -11,'  stop the fit here '
c
c***********************************************************************
c
c     experimental data are read in from file 10, using two lines
c     per nucleus. data are read in as
c     mass, charge, binding-energy [mev], charge radius [fm],
c     surface thickness [fm], unused, proton s.o. splitting,
c     neutron s.o. splitting, average crho at center, shell fluctuation
c     in crho at center, rho(3.6fm), rho(5.4fm), rho(7.2fm), unused;
c     the second line contains possible shift values ( switched by the
c     option faksub).
c     reading of data stops if mass<0 is encountered.
c     if experimental data is missing zero value should be typed in.
c
c***********   sample input file, unit 10  **************************
c
c 016 08 -0127.6 2.777 0.839 0.00 6.3 6.1 7.50 6.79 0.96 0.00 0.00 0.00
c         0000.0 0.000 0.100 0.00 0.0 0.0 0.00 0.00 0.00 0.00 0.00 0.00
c 040 20 -0342.1 3.845 0.978 0.00 0.0 0.0 8.84 7.87 4.12 0.32 0.00 0.00
c         0000.0 0.000 0.131 0.00 0.0 0.0 0.00 0.00 0.00 0.00 0.00 0.00
c 048 20 -0416.0 3.964 0.881 0.00 0.0 0.0 7.88 7.33 4.30 0.32 0.00 0.00
c         0000.0 0.000 0.125 0.00 0.0 0.0 0.00 0.00 0.00 0.00 0.00 0.00
c 208 82 -1636.4 6.806 0.900 0.00 0.0 0.0 6.17 5.94 6.23 5.79 1.70 0.00
c         0000.0 0.000 0.084 0.00 0.0 0.0 0.00 0.00 0.00 0.00 0.00 0.00
c -58 28 -0506.5 4.356 0.911 0.00 0.0 0.0 8.10 7.74 5.85 0.78 0.00 0.00
c         0000.0 0.000 0.209 0.00 0.0 0.0 0.00 0.00 0.00 0.00 0.00 0.00
c 090 40 -0783.9 5.040 0.957 0.00 0.0 0.0 0.00 0.00 0.00 0.00 0.00 0.00
c         0000.0 0.000 0.084 0.00 0.0 0.0 0.00 0.00 0.00 0.00 0.00 0.00
c 116 50 -0988.7 5.537 0.947 0.00 0.0 0.0 6.60 6.97 6.82 3.36 0.23 0.00
c         0000.0 0.000 0.128 0.00 0.0 0.0 0.00 0.00 0.00 0.00 0.00 0.00
c 124 50 -1050.0 5.640 0.908 0.00 0.0 0.0 6.11 6.55 6.57 3.49 0.25 0.00
c         0000.0 0.000 0.107 0.00 0.0 0.0 0.00 0.00 0.00 0.00 0.00 0.00
c-118 50 -1005.0 5.571 0.931 0.00 0.0 0.0 0.00 0.00 0.00 0.00 0.00 0.00
c         0000.0 0.000 0.084 0.00 0.0 0.0 0.00 0.00 0.00 0.00 0.00 0.00
c 999 120 0.0 0.0 0.0 -14.2829 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
c         0.0 0.0 0.0   0.0    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
c 999 130 0.0 0.0 0.0 -15.3053 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
c         0.0 0.0 0.0   0.0    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
c 999 140 0.0 0.0 0.0 -15.5711 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
c         0.0 0.0 0.0   0.0    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
c 999 150 0.0 0.0 0.0 -14.4938 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
c         0.0 0.0 0.0   0.0    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
c 999 160 0.0 0.0 0.0 -10.8203 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
c         0.0 0.0 0.0   0.0    0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
c 1. spalte: a
c 2. spalte: z
c spalte 01: bindungs-energie  mev
c spalte 02: diffraktions-radius   fm
c spalte 03: hautdicke  fm
c spalte 04: z.zt. leer
c spalte 05: l-s-aufspaltung proton in 16-o
c spalte 06: l-s-aufspaltung neutron in 16-o
c spalte 07: rho(r=0) * 100.   fm**(-3)   gemessene ladungsverteilung
c spalte 08: rho(r=1.8)*100.       "            "
c spalte 09: rho(r=3.6)*100.       "            "
c spalte 10: rho(r=5.4)*100.       "            "
c spalte 11: rho(r=7.2)*100.       "            "
c     the "untere zeile" supplies the relative
c     error weight for each fermi-momentum (it multiplies the
c     error given by 'fehler').
c
c********************************************************************
c
c
c
      include 'params.inc'
      include 'comfit.inc'
c
      dimension rmserr(12),rmserra(12)
      dimension a(kparmx),par(kparmx)
      dimension yy(12),fehler(12),fehlpr(12),faksub(12),dyy(12),
     *          sigmaa(kparmx),sigma0(kparmx),dparm0(kparmx),irand(30)
      REAL(8),PARAMETER :: parmfseries(4)=(/0.1D0,0.25D0,0.5D0,1.0D0/)
c
      character*6 zuord(kparmx)
      character*6 dummy
      character*27 heatyp(12),headact,headact2
      logical*4  top11
      logical tasym,trpaobs
      parameter (tasym=.false.)
      parameter (trpaobs=.true.)
      REAL(8) :: datapointnumber(12)      ! counter for effective nrr. of data points
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'numpar.inc'
      include 'result.inc'
      include 'constn.inc'
c
cold      data pi,sq2,h2m,hbarc,d4pi,hmc,omeg,sqfrpi,zero,one,x30lw
cold     *    /3.1415926,1.4142,20.75,197.3,0.07957742,0.2411,41.,3.54491,
cold     *     0.0,1.0,1.e-30/
      data x30lw/1.0D-30/
c
      data heatyp/'  binding-energy           ',
     *            '  diffr.-radius            ',
     *            '  surface thickness        ',
     *            '  rms radius               ',
     *            '  l*s splitting proton 16O ',
     *            '  l*s splitting neutron 16O',
     *            '  (rho(1.8)+rho(0.))/2.    ',
     *            '  neutron rms radius       ',
     *            '  gap-proton               ',
     *            '  gap-neutr.               ',
     *            '  rho(7.2) or 2. l*s neutr.',
     *            '  m*/m symm. nuclear matter'/
c
      data istchm,numax/3,14/
      data iextan,if3gap/0,0/
      data irand/1,1,1,0,1,1,0,1,0,1,1,0,1,0,0,1,1,0,1,1,1,0,0,
     *           1,1,0,1,0,1,0/
      data a/kparmx*0.01/
      data nlanmx/1/
c      data yexp/360*0.0/
c
c--------------------------------------------------------------------
c
       tsimple = .FALSE.
       parmnatfac = parmfseries(1)
c
       open(unit=1,form='formatted',status='old',file='for001')
       open(unit=2,form='formatted',status='unknown',file='for002')
       open(unit=3,form='formatted',status='unknown',file='for003')
       open(unit=7,form='formatted',status='unknown',file='for007')
       open(unit=10,form='formatted',status='old',file='for010')
       open(unit=21,form='formatted',file='for021')
       open(unit=22,form='formatted',file='for022')
c
c    initialise auxiliary fields for fit
c    (lifit=12 deactivates the fit-pointer)
c
      lifit=12
      yexp=0D0
c
c      outer loop. for each 'irun' a new computation is started.
c      loop stops if 'isteur<0' is read in.
c
      do 5 irun=1,1000
c
c     preset some control parameters.
c
      ngrid_explicit  = 0
      iftinp = 1
      ioptsk = 0
      ib2fix = 0
      ifrhoc = 0
      ifgsc  = 0
!      idisc=1
      idisc=-1
      endcon = 0.001
      ifalt2=0
      kama=12
      datapointnumber = 0D0
      iwr=  -2
      chisqr=.1e30
      chimin = chisqr
      istoch=1
      widpr  = 0D0
      widne  = 0D0
      wall   = -one
c
c**********************************************************************
c
c   first block: reading control parameters, fit options,
c                experimental data, and force parameters.
c
c**********************************************************************
c
      read(1,*) isteur,commen
      write(6,'(i6,a60)') isteur,commen
      if(isteur.lt.0) stop ' normal exit'
      write(02,'(43a)') (' * ',i=1,43)
      write(02,'(2x,a,2x,a)') 
     &  'fit of Skyrme-Hartree-Fock model to experimental data  --',
     &  'comment:'//commen
      read(1,*) nwdr,ityp,nparot,nwdro,chitmp
      write(6,'(a,i6)') 
     &   ' nwdr,ityp,nparot,nwdro=',nwdr,ityp,nparot,nwdro
c
c   zahl der datensaetze
c
      read(1,*) (ifit1(k),k=1,kama)
      read(1,*) (fehler(k),k=1,kama)
      read(1,*) (fehlpr(k),k=1,kama)
      read(1,*) (faksub(k),k=1,kama)
      write(6,'(12i6)') (ifit1(k),k=1,kama)
      write(6,'(12f6.3)') (fehler(k),k=1,kama)
      write(6,'(12f6.3)') (fehlpr(k),k=1,kama)
c
c   logik fuer anpassung von fermi-entf. radius:
c
      ifitr=0
      do k=1,kama
        if(ifit1(k).lt.0) ifitr=1
        ifit1(k)=iabs(ifit1(k))
      enddo
      if(ifitr.gt.0) write(2,'(//a//)')
     *          '  fit an fermi-entfaltete radien'
c
c   welcher datensatz enthaelt welche physikalische groesse?
c
      do l=1,kama
        if(ifit1(l).gt.0) lifit(ifit1(l))=l
      enddo
c
c    accumulate number of data. result  on 'ndats'.
c
      ndats=0
      do k=1,kama
        if(ifit1(k).ne.0) ndats=1+ndats
      enddo
      write(2,9) ndats,(ifit1(k),k=1,kama),(fehler(k),k=1,kama),
     +           (fehlpr(k),k=1,kama),(faksub(k),k=1,kama)
    9 format(1h0,'  fit an ',i2,' datensaetze'/
     +'  daten aus spalten ',12i6/
     +'                    ',12f6.2,'  fehler'/
     +'                    ',12f6.2,'  0: prozent,  1: absolut'/
     +'                    ',12f6.2,' * var. add. zum exp. wert')
      rewind 10
c
      iwrte=0
      if(nwdr.eq.0) then
!      if(nwdr.ge.0) then
       iwr=-1
       inquire(11,opened=top11)
       if(.not.top11) then
        open(unit=11,status='unknown',file='for011')
c        write(11,'(/40a)')  (' *',i=1,40)
       endif
       title = commen(1:50)
       ifplot=3
      endif
c
      do idat=1,kama
       dyy(idat)=0.
       do  j=1,kdat
        yexp(idat,j)=0D0
        sigmay(idat,j)=one
        yfit(idat,j)=0D0
       enddo
      enddo
      do j=1,kparmx
       sigmaa(j)=0.
       sigma0(j)=0.
       a(j)=.01
      enddo
c
      read(1,*) iaunt,iaob,maxdif,if3gap
c
c     limits 'iaunt' and 'iaob' now defined --> read data
c
      write(2,'(/a/)')
     *         '  **********  input  -  daten  ***********  '
c
c     loop 30: einlesen der exp. kerndaten , diesmal von file 10.
c
      read(10,'(a)') dummy
      if(dummy(1:1).ne.'#') backspace(10)
      nminv = 0                 ! overrides any later setting
      j=0
      do 30 i=1,kdat
      read(10,*) nra,zrprot,(yy(k),k=1,kama)
      nrprot = NINT(zrprot)
      if(nra.lt.0) go to 40
c
      iold8 = 0
      if(iold8.eq.1) then
        yyy7=yy(7)
        yyy8=yy(8)
        yy7=(yy(8)+yy(7))/2.
        yy8=(yy(8)-yy(7))
        yy(7)=yy7
        yy(8)=yy8
      endif
      read(10,*) (dyy(k),k=1,kama)
      write(6,'(2i5,1x,3g11.3,9f7.3)') nra,nrprot,(yy(k),k=1,kama)
      write(06,'(9x,f7.3,2f6.3,9f6.3)')  (dyy(k),k=1,kama)
      DO kact=1,kama
        IF(ABS(yy(kact))>1D-10) datapointnumber(kact)=
     &      datapointnumber(kact) + 1D0/dyy(kact)
      END DO
c
c   nra=a, nrprot=z,iint=zahl der radialen gitterpunkte,
c   nmaxp=zahl der proton-zustaende, nmaxn=zahl der neutronzust.,
c   yy(1)=bindungsenergie
c   yy(2)=dms-radius
c   yy(3)=sigma
c   yy(4)=rms-radius
c   yy(5)=schalenenergie-differenz protonen
c   yy(6)=        "          "     neutronen
c   yy(7)= rho(0.0)
c   yy(8)= rho(1.8)
c       yy(7) and yy(8) are recoupled to sum and difference.
c
c
c       skip if this observable is missing from the list
c
      if(nminv.gt.0) then
       if(abs(yy(nminv)).lt.x30lw) go to 30
      endif
c
c     check if we want to exclude any nucleus
c
      if(nra.lt.iaunt.or.nra.gt.iaob) go to 30
c
      j=j+1
      iay(j)=nra
      iprot(j)=nrprot
      zprot(j)=zrprot
cold      rnra=nra
cold      rnrpro=nrprot
c
      do 280 l=1,kama
      if(ifit1(l).eq.0) then
       yexp(l,j)=0D0
      elseif(faksub(ifit1(l)).ne.0D0 .and. yy(ifit1(l)).ne.0D0) then
!       yexp(l,j)=yy(ifit1(l))+dyy(ifit1(l))*faksub(ifit1(l))
       yexp(l,j)=yy(ifit1(l))
     &        +faksub(ifit1(l))*fehlpr(l)
     &        +faksub(ifit1(l))*yy(ifit1(l))*(1D0-fehlpr(l))/1D2
      else
       yexp(l,j)=yy(ifit1(l))
      endif
c
c   fermi-entfaltung der radien: p.b.9, seite 79
c
      if(ifitr.eq.1.and.l.eq.lifit(2))
     * yexp(l,j)=(yexp(l,j)+0.011-0.323*yexp(lifit(3),j))/0.975
      if(l.ne.lifit(8) .or. iold8.ne.1) then
c
c     prepare error weights
c
       if(ifit1(l).gt.0) then
         if(dyy(ifit1(l)).ne.zero) then
           fehlac=fehler(l)*dyy(ifit1(l))
         else
           fehlac=fehler(l)
         endif
       else
         fehlac=fehler(l)
       endif
       sigmay(l,j)=abs(fehlpr(l)*fehlac-(fehlpr(l)-1.)/100.*
     +                  fehlac*abs(yexp(l,j)))
      else
c       mittelwert r=0/r=1.8
       fehl1=fehlpr(lifit(7))*fehler(lifit(7))-(fehlpr(lifit(7))-1.)
     + /100.*fehler(lifit(7))*yyy7
       fehl2=fehlpr(lifit(8))*fehler(lifit(8))-(fehlpr(lifit(8))-1.)
     + /100.*fehler(lifit(8))*yyy8
       sigmay(l,j)=abs(fehl1)+abs(fehl2)
      endif
      if(abs(yexp(l,j)).lt.0.1e-10) sigmay(l,j)=0.1e-10
  280 continue
c
c   print experimental input data
c
      do ipri=1,(ndats-1)/8+1
        kpr1 = ipri*8-7
        kpr2  = min0(ndats,ipri*8)
        write(2,'(i3,i4,i3,8(f9.3,6x))')
     *             i,nra,nrprot,(yy(ifit1(k)),k=kpr1,kpr2)
        write(2,'(10x,8(f9.3,f6.3))') 
     &        (yexp(k,j),sigmay(k,j),k=kpr1,kpr2)
      enddo
c
   30 continue
      write(2,'(a)')
     & ' ! ATTENTION: more data (>50) in for010 than run here !'
   40 continue
      npts=j
      WRITE(2,'(a,12f7.2)') 'Nr. data points:',datapointnumber
c
c     optionally read isotopic differences or 3gaps and compute pointers
c
      IF(maxdif>0) THEN                     ! read isotopic differnces
        IF(maxdif>kdifmax) STOP 'nr. of isotopic differences too large'
        DO idif=1,maxdif
          READ(1,*) nuc1a,nuc1z,nuc2a,nuc2z,nucdif(3,idif),errdif(idif)
          nucdif(1,idif) = iposnuc(nuc1a,nuc1z)
          nucdif(2,idif) = iposnuc(nuc2a,nuc2z)
          IF(nucdif(1,idif)==0 .OR. nucdif(2,idif)==0) THEN
            WRITE(2,'(a,2i4,3x,2i4)')
     &      'one of 2 nuclei not in list, A,Z=',nuc1a,nuc1z,nuc2a,nuc2z
            STOP
          END IF
        END DO
      END IF
      IF(if3gap>0) THEN                     ! read 3gaps
        OPEN(9,file='for009')
        max3gap = 0
        DO i=1,k3gapmax
          READ(9,*,ERR=991,END=19) 
     &       nuc1a,nuc1z,nuc2a,nuc2z,nuc3a,nuc3z,err3gap(i),itypobs
          IF(itypobs<1 .AND. itypobs>12) 
     &         STOP 'ITYPOBS out of range. Check "for009"'
          n3gap(1,i) = iposnuc(nuc1a,nuc1z)
          n3gap(2,i) = iposnuc(nuc2a,nuc2z)
          n3gap(3,i) = iposnuc(nuc3a,nuc3z)
          n3gap(5,i) = itypobs
!          WRITE(6,*)'3gap',i,nuc1a,nuc1z,nuc2a,nuc2z,nuc3a,nuc3z,
!     &    itypobs,
!     &    iposnuc(nuc1a,nuc1z),iposnuc(nuc2a,nuc2z),iposnuc(nuc3a,nuc3z)          
          IF(n3gap(1,i)==0 .OR. n3gap(2,i)==0 .OR. n3gap(3,i)==0) THEN
            WRITE(2,'(a,2i4,3x,2i4,3x,2i4)')
     &      'one of 3 nuclei not in list, A,Z=',
     &      nuc1a,nuc1z,nuc2a,nuc2z,nuc3a,nuc3z
            STOP
          END IF
          IF(nuc1z==nuc2z .AND. nuc2z==nuc3z) THEN
            n3gap(4,i) = 9      ! proton gaps -- only for print
          ELSE IF(nuc1a-nuc1z==nuc2a-nuc2z .AND. 
     &            nuc2a-nuc2z==nuc3a-nuc3z) THEN
            n3gap(4,i) = 10     ! neutron gaps -- only for print
          ELSE
            STOP '3gaps do not match protons or neutrons'
          END IF          
          max3gap = i
        END DO
 19     CONTINUE
        CLOSE(9)
      END IF
c
c     continue reading fit options
c
      read(1,*) nminv,ifex,iftm,ifenrg,ifrso,ifzpe,ifrhoc,ifdfms,imfold
      read(1,*) acoul,addnew,dabe,dabr,ipair,ilnrsp,iparsh,cutfac,
     &          hfbwidth
      ipairsav=ipair
      IF(cutfac>0D0) itpr=10
      itpr=200
      IF(ipair>8) itpr=500
cold      if(ipair.eq.2) iparol = 0
      if(ipair.ge.5) x0damp = 0.2
      write(6,'(12i6)') iaunt,iaob
      write(6,'(8i10)')
     *          nminv,ifex,iftm,ifenrg,ifrso,ifzpe,ifrhoc,ifdfms
      write(6,'(3g12.5)') acoul
      write(2,223) acoul,ifex,iftm,ifenrg,ifrso,
     *             ifzpe,ifrhoc,ifdfms,ipair,ilnrsp,iparsh,imfold,
     *             cutfac,hfbwidth,isovectormodel
  223 format(/5x,
     +     'coulomb-kraft multipliziert mit', f6.3/
     +5x,'ifex  = ',i2,'  ( 0: ohne     , 1: mit coulomb exch.)'/
     +5x,'iftm  = ',i2,'  ( 0: ohne     , 1: mit t1-2 cont. ls)'/
     +5x,'if energie = ',i1,' : "1" rechne e-bindung ab initio,',
     +'  "0" gehe von bisherigem aus'/
     +5x,'ifrso = ',i2,'  ( 0: ohne     , 1: mit ls cont. rho )'/
     +5x,'ifzpe = ',i2,'  ( 0: m*(1-1/a) , 1: m but zpe corr. )'/
     +5x,'ifrhoc= ',i2,'  ( 0: rho_p in Coul., 1: rho_C in Coul. )'/
     +5x,'ifdfms= ',i2,'  ( 0: equal mass, 1: p-n mass diff.  )'/
     +5x,' ipair= ',i2,'  ( 0: const.gap; 1&2: BCS; 3: Lip.Nog. )'/
     +5x,'ilnrsp= ',i2,'  ( 1: response in Lipkin Nogami )'/
     +5x,'iparsh= ',i2,'  ( nr. of extra shells for BCS )'/
     +5x,'imfold= ',i2,'  ( 0: sing. fold. 1: doub.fold.  )'/
     +5x,'cutfac= ',f7.2, ' (cutoff factor in pairing)'/
     +5x,'hfbwidth= ',f6.2, ' (folding width in pairing)'/
     +5x,'isovectormodel= ',i2/)
      write(2,'(a,3(a,g12.3))') '     convergence parameters:',
     &          ' addnew=',addnew,' dabe=',dabe,' dabr=',dabr
c
c     for nuclear matter fit give compressed output on file 13
c
      if(iaunt.gt.900 .and. iaob.gt.900) then
ctest  do 9301 iread=1,9999
c9301  read(13,'(a)',err=9309,end=9309)  dummy
c9309  continue
       open(unit=13,form='formatted',status='unknown',file='for013')
       write(13,'(a/a/2a/1x)') '@',
     *   ' protocoll of nuclear matter fit :',
     *   '    comment=',commen
      endif
c
      if(nminv.gt.0) write(2,'(/a,i1,a)')
     *         '  nur die daten, fuer die auch ein wert in der',
     *         nminv,'. spalte vorhanden ist'
c
c   einlesen anfangs-parameter etc
c
      read(1,*) RSYMAX,NSYMAX,RNTMAX,NNTMAX                ! EoS analysis
      read(1,*) iwrte,ist1,ifcurv,itmax,iextan,ngrid_explicit
      read(1,*) nparin,numax,ivar,freep,endcon,iftinp,ib2fix,ioptsk
      nparfit=numax
!      if(nparin.ne.kparmx)
!     & stop ' nr. of parameters does not match KPARMX in the code'
      if(nparin.gt.kparmx)
     & stop ' nr. of parameters does larger than KPARMX'
      if(iftinp.eq.1) ib2fix=0
      if(ioptsk.eq.4 .and. iftinp.eq.1)
     &   stop ' IFTINP=1 not possible with IOPTSK=4'
      write(6,'(12i6)') iwrte,ist1
      write(6,'(3i6,f6.3,g12.5,3i3)')
     &       nparin,numax,ivar,freep,endcon,iftinp,ib2fix,ioptsk
      write(02,'(5x,a,g12.5)')
     *         'end condition for hartree iteration: endcon=',endcon
      if(ib2fix.eq.1)
     & write(2,'(5x,a)') 'b2=b1/2 and b2p=b1p/2 fixed.'
      call  naming(zuord,ioptsk)
!      IF(ioptsk==13 .AND. iftinp==6) x0damp=0.1     ! preliminary
c
!      IF(ioptsk==13) THEN
!        e0damp=160D0
!      END IF
c
c   eingabe der anfangsparameter
c
      a=0D0
      write(2,'(/a/)')
     *         '  nr. anfangswert  variation    nummer'
      do j=1,nparin
        read(1,*) a(j),dparmo(j),iparva(j)
        IF(nwdr==1) dparmo(j)=dparmo(j)*freep                        ! new
        write(6,'(2g13.6,i6)') a(j),dparmo(j),iparva(j)
        dparm0(j)=dparmo(j)
        par(j)=a(j)
        write(2,'(i5,2e14.7,i5)') j,a(j),dparmo(j),iparva(j)
        do l=1,kparmx
          if(iparva(j).eq.l) jp(l)=j
        enddo
      enddo
      if(par(jp(1)).eq.zero .and. par(jp(6)).eq.zero)
     &  write(2,'(a,a)') ' --> version: zero-range',
     &     ' two-body force switched off.'
c
c**********************************************************************
c
c   second block: performing the fit by repeated calls to the
c                fitting subroutine 'xd1'.
c
c**********************************************************************
c
c     determine number of datas and degrees of freedom for the fit
c
      free=0
      do idat=1,ndats
        free1=0.
        do l=1,npts
          if(abs(yexp(idat,l)).ge.1.e-30) free1=one+free1
        enddo
        write(2,'(a,i1,a,f4.0,a)')
     *       ' im ',idat,'. satz ',free1,' daten von null verschieden'
        free=free+free1
      enddo
      nptsg=free
      free=free-numax
      write(2,342) nptsg,numax,free,freep
  342 format(/i4,' messwerte verschieden von null'/
     +i4,' parameter anzupassen'/f5.0,' freiheitsgrade'/
     +'  parameter variation = sigma0 * ',f8.3//)
      flamda=0.001
!
!       call attach(a,jp)
!       WRITE(6,'(a)') 'TEST PAR:'
!       WRITE(6,'(6(1pg13.5))') avp,h1vp,h2vp,avm,h1vm,h2vm
!       WRITE(6,'(6(1pg13.5))') asp,hsp,hgradp,asm,hsm,hgradm
c
c     prepare parameters for nuclear matter analysis in case of nwdr<0
c
      if(nwdr.lt.0) then
       call attach(a,jp)
       call initfc
c
c   einmalige rechnung ohne sigmaa (nwdr=0) mit sigmaa (nwdr=1)
c
      elseif(nwdr.eq.0) then
       numax = 1                            ! to override check
       chia=chisqr
       iwrte=1
       itera=1
       call xd1(a,sigmaa,sigma0,numax,nwdr,iwarn,1,ifcurv,zuord)
       ifplot = 3!0
c
c   iteration with gradient in chi**2-landscape
c
      else
       ifplot = 0
       if(ityp.eq.1) then
         call simplex(a,numax,nwdr)
       elseif(ityp.eq.2) then
         iwr=-2
         ifplot=0
         call powell(a,numax,nwdr)
       elseif(ityp.eq.3) then
         call splitb(zuord,a,sigmaa,sigma0,numax,nparot,
     *               nwdr,nwdro,ifcurv,jp)
       elseif(ityp.eq.4) then
         call mcsear(zuord,a,sigmaa,sigma0,numax,nparot,
     *         nwdr,nwdro,ifcurv,chitmp)
       else
         call beving(zuord,a,sigmaa,sigma0,numax,nwdr,ifcurv)
       endif
c
c      immediate print of result, in order of reading.
c
       write(2,'(a//3a)')
     *         ' printout of result in the order of reading:',
     *         ' read nr.  par.nr.   type      par.value   ',
     *         '   corr.error     uncorr.error  '
       do j=1,nparin
         itypa = iparva(j)
         write(2,'(1x,2i9,1x,a6,1x,3g16.7)')
     *    j,itypa,zuord(itypa),a(j),sigmaa(j),sigma0(j)
       enddo
c
      endif
c
c**********************************************************************
c
c   third block: printing and nuclear matter analysis.
c
c**********************************************************************
c
      iwr = -1                     ! enable printing
!      ifplot=3
      do i=1,nparin
        par(i)=a(i)
      enddo
      call attach(par,jp)
c
c     optional ensemble of "reasonable" parametrizations
      IF(iextan==4) CALL sample_reasonable(a,sigmaa,nparin,numax)
c
c
c      output of resulting parameters in ftn01-input format
c
       if(nwdr.ne.0)
     &   call xd1(a,sigmaa,sigma0,numax,0,iwarn,1,ifcurv,zuord)
       write(02,'(a)') ' resulting param.s in hffit input format:'
       do j=1,nparin
         write(02,'(g13.6,a,g12.5,a,i2,a,a)')  
     &     a(j),',',dparmo(j),',',
     &     iparva(j),',',zuord(iparva(j))
       enddo
c
c      output in ftn05 input format
c
      write(02,'(a)') ' resulting param.s in hafo input format:'
      open(unit=15,file='for015')
      do iu=2,15,13
       write(iu,'(a)')
     &   '1,,    IWUNIT',
     &   "'"//commen(1:60)//"'",
     &   '0,0.3,              NGRID,RSTEP',
     &   '-1,1E-8,1E-9,278,5780,3,-1, IWR,DABR/E,ITPR/MAX,IFPL',
     &   '0.2,90.,0.2,-1,,               X0DAMP,E0DAMP,ADDNEW,WALL'
       write(iu,'(9(i1,a),i2,a,5x,a)') 
     &    ifex,',',iftm,',',ifenrg,',',ifrso,',',ifzpe,',',
     &    ifdfms,',',ifrhoc,',',ifenav,',',imfold,',',ioptsk,
     &    ', if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk'
       write(iu,'(a)')
     &   '1,,,                       ifform',
     &   '0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12',
     &   '0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12'
c       write(iu,'(5(f10.3,a)/4(f10.6,a),f11.5,a)')
       write(iu,'(5(1pg18.10,a)/4(1pg18.10,a),0pf12.6,a)')
     &  t0,',',t1,',',t2,',',t3,',',t4,',      t0,t1,t2,t3,t4',
     &  x0,',',x1,',',x2,',',x3,',',b4p,',       x0,x1,x2,x3,b4p'
       if(ioptsk.eq.5 .or. ioptsk.eq.6 .or. ioptsk.eq.11) then
         write(iu,'(2(f10.6,a)/2(f10.6,a))')
     &      c1,',',c4,',                      c1 c4',
     &      c1p,',',c4p,',                    c1p,c4p'
       endif
       if(ioptsk.eq.8 .or. ioptsk.eq.11) then
         write(iu,'(1(f10.3,a)/2(f10.3,a))')
     &      b3tilde,',                      b3tilde',
     &      b12,',',b12p,',                     b12,b12p'
       endif
       if(ioptsk.eq.9 .or. ioptsk.eq.11) then
         write(iu,'(2(f10.3,a)/f10.6,a)')
     &      b7,',',b7p,',               b7,b7p',
     &      pow7,',               pow7'
       endif
       if(ioptsk.eq.10) then
         write(iu,'(3f10.5,a)')
     &      riftm,rho0pr,powerp,',      riftm,rho0pr,powerp'
       endif
       if(ioptsk.eq.11) then
         write(iu,'(8f10.6,a)')
     &      riftm,rho0pr,powerp,ecstab,alpinv,betinv,
     &      alpin2,betin2,
     &     ',  riftm,rho0pr,powerp,ecstab,alpinv'
       endif
       if(csurf.ne.zero.and.ioptsk.ne.1) then
         write(iu,'(5(f10.3,a)/4(f10.6,a))')
     &    t5,',',t6,',',t7,',',t8,',',t9,',    t5,t6,t7,t8,t9',
     &    x5,',',x6,',',x7,',',x8,',           x5,x6,x7,x8'
       endif
       write(iu,'((f10.3,a),2(f10.6,a))')
     &    tyuk,',',xyuk,',',ayuk,',           tyuk,xyuk,ayuk'
       write(iu,'(f14.10,a,i2,a,f12.5,a,f12.5,a)')
     &     power,',1.0,',ipair,',',pairfp,',',pairfn,
     &  ',0,,,,,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,...'
      enddo
      close(15)
c
c     table of chi**2 finally achieved
c
       write(02,'(/a)')
     &  'table of chi**2 finally achieved:'
        write(2,400)
  400   format(3x,' a  z',6x,'f(exp)',7x,'f(fit)',6x,'difr(%)',
     :  1x,'stdev(%) ',2x,'chi',4x,'sum chi**2 ')
       WRITE(16,'(a)') '#table of errors:'
       WRITE(16,'(a)') '#    a  z    value  deviat   extrap.error'
c
        totchi=0
        do idat=1,ndats
         write(2,'(1x,a)') heatyp(ifit1(idat))
         IF(idat>1) WRITE(16,'(1x/1x)')
         write(16,'(a)') '# '//heatyp(ifit1(idat))
c
         chi1=0D0
         sumerr = 0D0    ! for % error
         averr  = 0D0    ! for % error
         numerr = 0      ! for % error
         sumerra = 0D0    ! for abs. error
         averra  = 0D0    ! for abs. error
         xnumerra = 0D0   ! for abs. error
         varerror = 0D0   ! mean square extrapolation error
         iprotold = iprot(1)
         do j=1,npts
          chi0=0.
          difr=0.
          sigpct=0.
          if(abs(yexp(idat,j)).ge.x30lw) then
c           chi=(yexp(idat,j)-yfit(idat,j))/sigmay(idat,j)
           chi0=chi(idat,j)
           chi1=chi1+chi0**2
           totchi=totchi+chi0**2
           if(yexp(idat,j).ne.zero .and. iay(j).lt.900) then
            difr=(yexp(idat,j)-yfit(idat,j))*100./yexp(idat,j)
            averr  = difr      + averr
            sumerr = difr*difr + sumerr
            numerr = 1 + numerr
            sigpct=abs(sigmay(idat,j)*100./yexp(idat,j))
            difr=(yexp(idat,j)-yfit(idat,j))
            vol = 1.0/sigmay(idat,j)**2
            averra  = vol*difr      + averra
            sumerra = vol*difr*difr + sumerra
            xnumerra = vol + xnumerra
            varerror = vol*varextrap(idat,j) + varerror
c            WRITE(2,'(a,2i5,3f8.2)') 
c     &         ' IDAT,J=',idat,j,difr**2*vol,xnumerra/vol,1D0/SQRT(vol)
           endif
          endif
          IF(abs(yexp(idat,j)).gt.x30lw) THEN
c            WRITE(2,'(2i4,3x,2(1x,g12.5),5(1x,f7.2),e12.5,f9.4)')
            WRITE(2,'(2i4,3x,2(1x,g14.6),5(1x,f7.2),e12.5,f9.4)')
     &            iay(j),iprot(j),yexp(idat,j),yfit(idat,j),difr,sigpct,
     &            chi0,totchi,chi1
            IF(iprot(j).NE.iprotold) WRITE(16,'(1x)')
            WRITE(16,'(2i4,3(1pg13.5))') iay(j),iprot(j),
     &            yfit(idat,j),ABS(difr),SQRT(varextrap(idat,j))
            iprotold=iprot(j)
          END IF
c
         enddo
         if(numerr.gt.0)  then
            rmserr(idat) = sqrt(sumerr/numerr)
            write(2,'(21x,a,f7.2,a,f7.2)') 
     &        ' rms error in %=',rmserr(idat),
     &        ' average error in %=',averr/numerr
            rmserra(idat) = sqrt(sumerra/xnumerra)
            write(2,'(21x,a,1pg12.4,a,1pg12.4)') 
     &        ' r.m.s. error (abs)=',rmserra(idat),
     &        ' average error (abs)=',averra/xnumerra
            varerror = sqrt(varerror/xnumerra)
            write(2,'(21x,a,1pg12.4,a,1pg12.4)') 
     &        ' rms extrap. error =',varerror,
     &        ' ratio to rms error =',varerror/rmserra(idat)
         endif
        enddo
!        write(21,'(a,12(1pg12.3))')
!     &     commen(1:12),(rmserr(idat),idat=1,ndats)
!
!       protocol chi**2 for isotopic differences
!
        if(maxdif.gt.0) then
          write(2,'(a/a)') '    isotopic differences',
     &    ' observable  A_1   Z_1   A_2   Z_2     chi    sumchi  subchi'
          idat = 13
          chi1 = 0D0
          do idif=1,maxdif
            inuc1 = nucdif(1,idif)
            inuc2 = nucdif(2,idif)
c?          iobs = lifit(nucdif(3,idif))
            chi0=chi(idat,idif)
            chi1=chi1+chi0**2
            totchi=totchi+chi0**2
            headact=heatyp(nucdif(3,idif))
            write(2,'(a,2i4,2x,2i4,2x,3(1pg12.4))')
     &      headact(1:12),
     &      iay(inuc1),iprot(inuc1),iay(inuc2),iprot(inuc2),
     &      chi0,totchi,chi1
            CALL FLUSH(2)
          enddo
       endif
!
!       protocol chi**2 for 3gaps
!
        IF(max3gap>0) THEN
          WRITE(2,'(a/a)') ' 3-point gaps:',
     &    ' observable  A_1,Z_1  A_2,Z_2 A3,Z3   chi    sumchi  subchi'
          idat = 14
          chi1 = 0D0
          erraver = 0D0
          errrms = 0D0
          DO i3=1,max3gap
            inuc1 = n3gap(1,i3)
            inuc2 = n3gap(2,i3)
            inuc3 = n3gap(3,i3)
            chi0=chi(idat,i3)
            erract = chi0*err3gap(i3)
            erraver = erraver+erract
            errrms = errrms + erract**2
            chi1=chi1+chi0**2
            totchi=totchi+chi0**2
            headact = heatyp(n3gap(4,i3))
            headact2 = heatyp(n3gap(5,i3))
            write(2,'(a,3(2i4,2x),3(1pg12.4))')
     &      headact(1:12)//headact2(1:12),iay(inuc1),iprot(inuc1),
     &      iay(inuc2),iprot(inuc2),iay(inuc3),iprot(inuc3),
     &      chi0,totchi,chi1
            CALL FLUSH(2)
          END DO
       END IF
       erraver = erraver/max3gap
       errrms = errrms/max3gap
       WRITE(2,'(a,1pg13.5,a,1pg13.5)')
     &   ' r.m.s. error=',SQRT(errrms-erraver**2),
     &   ' average error=',erraver
c
c       write(02,'(a,2(a,g12.5))') ' the level width in pb208',
c     *          ': for protons =',widpr,
c     *          ',  for neutrons=',widne
c       write(02,'(a,2(a,g12.5))') ' the shell fluctuation in',
c     *          ' pb208 (times 100) =',fluc
c
c    in case of nuclear matter fit give protocoll on file 13
c
      if(iaunt.gt.900 .and. iaob.gt.900) then
       write(13,'(2a)')
     *   ' nr. fer.mom. input e/a   adop.err.    result e/a  ',
     *   ' res. err.    chi**2 '
       totchi=zero
       do  j=1,npts
         difr=yexp(1,j)-yfit(1,j)
         chi0=difr/sigmay(1,j)
         chi0=chi0*chi0
         totchi=chi0+totchi
         write(13,'(1x,i3,3x,f5.2,2(2x,f8.2),4x,3(2x,f8.2))')
     &      j,iprot(j)/100.,yexp(1,j),sigmay(1,j),yfit(1,j),difr,chi0
       enddo
       write(13,'(58x,f8.2)') totchi
       istchm = istoch-1
       write(13,'(1x,i5,a,i5,a,i3,a)')
     *   itera,' nr. of iterations, out of',nwdr,'. thereby',
     *   istchm,' stochastic restarts.'
       write(13,'(a)') ' resulting param.s in hffit input format:'
       do 1039 j=1,nparin
       write(13,'(g13.6,a,g12.5,a,i2,a,a)')  a(j),',',dparmo(j),',',
     &                 iparva(j),',',zuord(iparva(j))
 1039  continue
cold       ipri13=1
cold      else
cold       ipri13=0
      endif
c
c     test print surface energies
c
      IF(.NOT.tsimple) THEN
      call attach(par,jp)
      if(abs(b4p).gt.1.0E-7) T4 = (t4+b4p)*0.666666667
      IF(tasym) THEN
        call ASYMT(T0,T1,T2,T3,t4,X0,X1,X2,X3,power,
     &           (amfc0(1)+amfc0(2))/2D0,
     +           YRHO0,eovera,xKOMPr,a4,YEFFM,YXKF,YELLO,
     +           surfe,surfas,YH,YP,YG,YQT1,YQT2)
        WRITE(2,'(a,4(1pg13.5))') 'testprint surface energies:',
     &     surfe,surfas,YQT1,YQT2
      END IF
!c
!c         early analysis
!c
!          parmnatfac = parmfseries(1)
!          call nucvarnat(par,numax)
!          stop 'early'

c

      WRITE(2,'(a,1pg18.10)') 'final chi**2=',totchi

c
c     nuclear matter analysis
c
      iwr =-1
      if(nwdr.eq.0) then
       call nucvar(par,0,zuord)
      else
       call nucvar(par,numax,zuord)
      endif
      if(nsymax.gt.0) then
        rhostp = rsymax/nsymax
        write(2,'(a/2a)')
     &   ' EoS in symmetric nuclear matter:',
     &   '  rho   k_F      E/A      variance    Ekin/A    Epot/A',
     &   '    a_sym '
        do irho=1,nsymax
c         do irho = 1,20
c          fkf = irho*0.1
          rhoact = irho*rhostp
          fkf = (rhoact/0.0675475)**0.333333333D0
c          rhoact=0.0675475*fkf*fkf*fkf
           call nucvar2(par,numax,rhoact,half,energm,errenm)
c          call enernm(rhoact,half,zero,energm,efmasp,efmasn)
          ekinoa = ekinnm(rhoact,half)
!          write(2,'(2f8.4,7(1pg13.5))') 
          write(2,'(2f8.4,7(1pg16.8))') 
     &      rhoact,fkf,energm,errenm,ekinoa,energm-ekinoa,
     &      asymrho(rhoact)
        enddo
      endif
      if(nntmax.gt.0) then
        rhostp = rntmax/nntmax
        write(2,'(2a)')
     &   ' EoS in neutron matter:','  rho     E/A  variance '
        write(28,'(2a)')
     &   ' EoS in neutron matter:','  rho     E/A'
        do irho=1,nntmax
          rhoact = irho*rhostp
          call nucvar2(par,numax,rhoact,zero,energm,errenm)
c          call enernm(rhoact,zero,zero,energm,efmasp,efmasn)
          ekinoa = ekinnm(rhoact,zero)
          write(2,'(6(1pg12.4))') 
     &         rhoact,energm,errenm,ekinoa,energm-ekinoa
          write(28,'(4(1pg12.4))') rhoact,energm
        enddo
      endif
      if(nlanmx.gt.0) then
        nlanmx = max(nsymax,1)
        if(nlanmx.gt.1) then
          delrho = rsymax/nlanmx
          write(2,'(a/2a)')
     &    ' Landau parameter:',
     &    '    rho         f0        f0p        g0         g0p  ',
     &    '    c0fac    f1     f1p   F0    F0p'
          do nl=0,nlanmx
            call attach(par,jp) 
            rhoact = nl*delrho
!            WRITE(2,*) 'LANPAR called'
!            c0fac2=c0fac/150D0
            call lanpar(rhoact,flan0,flan0p,glan0,glan0p,f0,f0p,g0,g0p,
     &                  f1,f1p,c0fac)
            write(2,'(1x,f8.4,9(1pg12.4))') rhoact,f0,f0p,g0,g0p,c0fac,
     &      f1,f1p,(f0-f1/3D0)/(1D0-f1/3D0),(f0p-f1/3D0)/(1D0-f1/3D0)
          enddo
!        else
        endif
          call nucmat(rhonm,eovera,compr,effmas,f0,f0p,g0,g0p,
     &            asymme,dasymm,d2asymm,sk,0,-1)
          call lanpar(rhonm,flan0,flan0p,glan0,glan0p,f0,f0p,g0,g0p,
     &                f1,f1p,c0fac)
!          c0fac2=c0fac/150D0
          write(2,'(a/2a)')
     &    ' Landau parameter at equilibrium:',
     &    '    rho         f0        f0p        g0         g0p    ',
     &    '    c0fac    f1    f1p'
          write(2,'(1x,f8.4,9(1pg12.4))') rhonm,f0,f0p,g0,g0p,c0fac,
     &      f1,f1p,(f0-f1/3D0)/(1D0-f1/3D0),(f0p-f1/3D0)/(1D0-f1/3D0)
!          write(2,'(a,3(1pg13.5))') ' kappa=',(1D0+f1p/3D0)/effmas,
!     &       (f1p-f1)/(1D0+f1/3D0)/3D0
          write(2,'(1x)')
!        endif
      endif
      END IF
!          CALL  extos3gap(par,1,34,50,34,50,0)
!          CALL  extobs3gap(par,1,34,50,34,50,1)
!          CALL  extobs3gap(par,1,34,50,34,50,2)
!          CALL  extobs3gap(par,1,20,24,20,24,0)
!          CALL  extobs3gap(par,1,20,24,20,24,1)
!          CALL  extobs3gap(par,1,20,24,20,24,2)
c
c     asymmetric matter
c
c      CALL asym_equi()
c
c     other nuclei analysis
c
c      iwr =-2
c        call extobsmf(par,yfit,numax*min(nwdr,1),1,0,0,0,0,0,0,0,0,0,
c     &                0,0,0,0,0,0,0)
      if(iextan.eq.1) call extnuc(par,yfit,numax*min(nwdr,1))
      IF(tsimple) THEN
        DO ifac=1,4
          parmnatfac = parmfseries(ifac)
          WRITE(2,'(/a,f7.3)') ' parameter variation now:',
     &      parmnatfac
          call extobs2_all(par,yfit,numax*min(nwdr,1),37)
       END DO
      ELSE
      if(iextan.eq.2) then
!        IF(iextan<0) then
        call extbarr(par,numax*min(nwdr,1))
        CALL attach(par,jp)
        IF(trpaobs) THEN
!          call extobs_all(par,yfit,numax*min(nwdr,1),1)
!          call extobs_all(par,yfit,numax*min(nwdr,1),2)
!          call extobs_all(par,yfit,numax*min(nwdr,1),3)
!          call extobs_all(par,yfit,numax*min(nwdr,1),4)
!          call extobs_all(par,yfit,numax*min(nwdr,1),5)
!          call extobs_all(par,yfit,numax*min(nwdr,1),6)
!          call extobs_all(par,yfit,numax*min(nwdr,1),11)
        END IF
        CALL flush(2)
        call extobs_all(par,yfit,numax*min(nwdr,1),21)
        call extobs_all(par,yfit,numax*min(nwdr,1),31)
        call extobs_all(par,yfit,numax*min(nwdr,1),33)
        call extobs_all(par,yfit,numax*min(nwdr,1),34)
        call extobs_all(par,yfit,numax*min(nwdr,1),35)
        call extobs_all(par,yfit,numax*min(nwdr,1),36)
        call extobs_all(par,yfit,numax*min(nwdr,1),38)
        call extobs_all(par,yfit,numax*min(nwdr,1),39)
        call extobs_all(par,yfit,numax*min(nwdr,1),43)
        call extobs_all(par,yfit,numax*min(nwdr,1),27)
        call extobs_all(par,yfit,numax*min(nwdr,1),22)
        call extobs_all(par,yfit,numax*min(nwdr,1),23)
        call extobs_all(par,yfit,numax*min(nwdr,1),24)
        call extobs_all(par,yfit,numax*min(nwdr,1),25)
        call extobs_all(par,yfit,numax*min(nwdr,1),26)
        call extobs_all(par,yfit,numax*min(nwdr,1),28)
!        call extobs_all(par,yfit,numax*min(nwdr,1),29)
!        call extobs_all(par,yfit,numax*min(nwdr,1),30)
!        call extobs_all(par,yfit,numax*min(nwdr,1),6)
!        call extobs_all(par,yfit,numax*min(nwdr,1),7)
!        call extobs_all(par,yfit,numax*min(nwdr,1),8)
!        call extobs_all(par,yfit,numax*min(nwdr,1),9)
!        call extobs_all(par,yfit,numax*min(nwdr,1),32)
!        call extobs_all(par,yfit,numax*min(nwdr,1),10)
        call extobs_all(par,yfit,numax*min(nwdr,1),37)
!        call extobs_all(par,yfit,numax*min(nwdr,1),40)
!        call extobs_all(par,yfit,numax*min(nwdr,1),12)
!        call extobs_all(par,yfit,numax*min(nwdr,1),13)
!        call extobs_all(par,yfit,numax*min(nwdr,1),14)
!        call extobs_all(par,yfit,numax*min(nwdr,1),15)
!        call extobsFq(par,yfit,numax*min(nwdr,1))
!        endif
!        DO ifac=-5,5
!          parmnatfac = ifac*0.2D0
!          WRITE(2,'(/a,f7.3)') ' parameter variation now:',
!     &      parmnatfac
!          call nucvarnat(par,numax)
!          call extobs2_all(par,yfit,numax*min(nwdr,1),1)
!          call extobs2_all(par,yfit,numax*min(nwdr,1),2)
!          call extobs2_all(par,yfit,numax*min(nwdr,1),3)
!          call extobs2_all(par,yfit,numax*min(nwdr,1),23)
!        END DO
!        IF(numax>1) THEN
        CALL flush(2)
        IF(numax>0) THEN
        DO ifac=1,4,3  !4,1,-3 !1,4
          parmnatfac = parmfseries(ifac)
          WRITE(2,'(/a,f7.3)') ' parameter variation now:',
     &      parmnatfac
!          WRITE(2,'(5(1pg13.5))') parmnat(1:numax,1:numax)
!          call introbs(par,numax*min(nwdr,1))
!          DO nn=16,32,2     ! 44
!            CALL  extobsiso(par,yfit,numax*min(nwdr,1),20,28,20,nn)
!          END DO
!          CALL  extobs3gap(par,numax*min(nwdr,1),20,24,20,24,1)
!          CALL  extobs3gap(par,numax*min(nwdr,1),50,82,50,82,3)
!          CALL  extobs3gap(par,numax*min(nwdr,1),82,126,82,126,3)
          call nucvarnat(par,numax)
          call extbarr2(par,numax*min(nwdr,1))
          call extobs2_all(par,yfit,numax*min(nwdr,1),41)
          call extobs2_all(par,yfit,numax*min(nwdr,1),42)
          IF(trpaobs) THEN
            call extobs2_all(par,yfit,numax*min(nwdr,1),1)
            call extobs2_all(par,yfit,numax*min(nwdr,1),2)
            call extobs2_all(par,yfit,numax*min(nwdr,1),3)
            call extobs2_all(par,yfit,numax*min(nwdr,1),4)
            call extobs2_all(par,yfit,numax*min(nwdr,1),5)
            call extobs2_all(par,yfit,numax*min(nwdr,1),6)  !  prel.
            call extobs2_all(par,yfit,numax*min(nwdr,1),7)
          END IF
          CALL flush(2)
!          call extobs2_all(par,yfit,numax*min(nwdr,1),11)
          call extobs2_all(par,yfit,numax*min(nwdr,1),21)
          call extobs2_all(par,yfit,numax*min(nwdr,1),31)
          call extobs2_all(par,yfit,numax*min(nwdr,1),30)
          call extobs2_all(par,yfit,numax*min(nwdr,1),22)
          call extobs2_all(par,yfit,numax*min(nwdr,1),33)
          call extobs2_all(par,yfit,numax*min(nwdr,1),34)
          call extobs2_all(par,yfit,numax*min(nwdr,1),38)
          call extobs2_all(par,yfit,numax*min(nwdr,1),39)
          call extobs2_all(par,yfit,numax*min(nwdr,1),43)
          call extobs2_all(par,yfit,numax*min(nwdr,1),37)
          call extobs2_all(par,yfit,numax*min(nwdr,1),23)
          call extobs2_all(par,yfit,numax*min(nwdr,1),24)

!tmp          call extobs2_all(par,yfit,numax*min(nwdr,1),4)
!tmp          call extobs2_all(par,yfit,numax*min(nwdr,1),35)
          IF(tformw) THEN
            call extobs2_all(par,yfit,numax*min(nwdr,1),36)
            call extobs2_all(par,yfit,numax*min(nwdr,1),34)
            call extobs2_all(par,yfit,numax*min(nwdr,1),44)
            call extobs2_all(par,yfit,numax*min(nwdr,1),40)
          END IF
!          call extobs2_all(par,yfit,numax*min(nwdr,1),31)
          call extobs2_all(par,yfit,numax*min(nwdr,1),32)
!tmp          call extobs2_all(par,yfit,numax*min(nwdr,1),25)
!tmp          call extobs2_all(par,yfit,numax*min(nwdr,1),26)
!          call extobs2_all(par,yfit,numax*min(nwdr,1),28)
!          call extobs2_all(par,yfit,numax*min(nwdr,1),29)
!          call extobs2_all(par,yfit,numax*min(nwdr,1),30)
!tmp          call extobs2_all(par,yfit,numax*min(nwdr,1),6)
!tmp          call extobs2_all(par,yfit,numax*min(nwdr,1),7)
!tmp          call extobs2_all(par,yfit,numax*min(nwdr,1),8)
!tmp          call extobs2_all(par,yfit,numax*min(nwdr,1),16)
!          call extobs2_all(par,yfit,numax*min(nwdr,1),9)
!          call extobs2_all(par,yfit,numax*min(nwdr,1),32)
!          call extobs2_all(par,yfit,numax*min(nwdr,1),10)
!          call extobs2_all(par,yfit,numax*min(nwdr,1),12)
!          call extobs2_all(par,yfit,numax*min(nwdr,1),13)
!          call extobs2_all(par,yfit,numax*min(nwdr,1),14)
          call extobs2_all(par,yfit,numax*min(nwdr,1),15)
!tmp          call extobs2_all(par,yfit,numax*min(nwdr,1),27)
!tmp          call extobs2_all(par,yfit,numax*min(nwdr,1),11)
!tmp          call extobsFq2(par,yfit,numax*min(nwdr,1))
!          call extobsLanPar(par,yfit,numax*min(nwdr,1))
!          DO ir=30,40,2
!            CALL  extobsiso(par,yfit,numax*min(nwdr,1),28,ir-2,28,ir)
!          END DO
!          DO ir=32,40,2
!            CALL  extobsiso(par,yfit,numax*min(nwdr,1),28,30,28,ir)
!          END DO
!          mirror series
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),28,26,26,28)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),28,24,24,28)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),28,22,22,28)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),28,20,20,28)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),20,16,16,20)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),21,20,21,19)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),21,19,20,19)
          CALL  extobsiso(par,yfit,numax*min(nwdr,1),20,20,20,16)
          CALL  extobsiso(par,yfit,numax*min(nwdr,1),20,20,20,24)
          CALL  extobsiso(par,yfit,numax*min(nwdr,1),20,24,20,28)
          CALL  extobsiso(par,yfit,numax*min(nwdr,1),20,20,20,28)
          CALL  extobsiso(par,yfit,numax*min(nwdr,1),50,74,50,82)
          CALL  extobsiso(par,yfit,numax*min(nwdr,1),50,84,50,82)
          CALL  extobsiso(par,yfit,numax*min(nwdr,1),82,120,82,126)
          CALL  extobsiso(par,yfit,numax*min(nwdr,1),82,132,82,126)
          CALL flush(2)
          IF(1==1) THEN 
          CALL  extobs3gap(par,numax*min(nwdr,1),20,24,20,24,2)
          CALL  extobs3gap(par,numax*min(nwdr,1),50,68,50,72,2)
!          CALL  extobs3gap(par,numax*min(nwdr,1),14,16,18,16,2)
          CALL  extobs3gap(par,numax*min(nwdr,1),38,50,38,50,2)
          CALL  extobs3gap(par,numax*min(nwdr,1),44,50,44,50,2)
          CALL  extobs3gap(par,numax*min(nwdr,1),52,82,52,82,2)
          CALL  extobs3gap(par,numax*min(nwdr,1),58,82,58,82,2)
          CALL  extobs3gap(par,numax*min(nwdr,1),88,126,88,126,2)
          CALL  extobs3gap(par,numax*min(nwdr,1),20,24,20,24,3)
          CALL  extobs3gap(par,numax*min(nwdr,1),50,68,50,72,3)
          CALL  extobs3gap(par,numax*min(nwdr,1),82,120,82,120,3)
          CALL  extobs3gap(par,numax*min(nwdr,1),38,50,38,50,3)
          CALL  extobs3gap(par,numax*min(nwdr,1),44,50,44,50,3)
          CALL  extobs3gap(par,numax*min(nwdr,1),52,82,52,82,3)
          CALL  extobs3gap(par,numax*min(nwdr,1),58,82,58,82,3)
          CALL  extobs3gap(par,numax*min(nwdr,1),88,126,88,126,3)
          CALL flush(2)
          CALL  extobs3gap(par,numax*min(nwdr,1),20,44,20,44,0)
          CALL  extobs3gap(par,numax*min(nwdr,1),50,73,50,73,0)
          CALL  extobs3gap(par,numax*min(nwdr,1),82,101,82,101,0)
          CALL  extobs3gap(par,numax*min(nwdr,1),40,50,40,50,0)
          CALL  extobs3gap(par,numax*min(nwdr,1),56,82,56,82,0)
          CALL  extobs3gap(par,numax*min(nwdr,1),88,126,88,126,0)
          CALL  extobs3gap(par,numax*min(nwdr,1),20,44,20,44,1)
          CALL  extobs3gap(par,numax*min(nwdr,1),50,73,50,73,1)
          CALL  extobs3gap(par,numax*min(nwdr,1),82,101,82,101,1)
          CALL  extobs3gap(par,numax*min(nwdr,1),40,50,40,50,1)
          CALL  extobs3gap(par,numax*min(nwdr,1),56,82,56,82,1)
          CALL  extobs3gap(par,numax*min(nwdr,1),88,126,88,126,1)
          CALL flush(2)
!          DO ir=0,1
!            CALL  extobs3gap(par,numax*min(nwdr,1),50,73,50,73,ir)
!            CALL  extobs3gap(par,numax*min(nwdr,1),82,121,82,121,ir)
!*            CALL  extobs3gap(par,numax*min(nwdr,1),57,82,57,82,ir)
!            CALL  extobs3gap(par,numax*min(nwdr,1),82,126,82,126,3)
!            CALL flush(2)
!          END DO
!          DO ir=0,3
!            CALL  extobs3gap(par,numax*min(nwdr,1),20,24,20,24,ir)
!            CALL  extobs3gap(par,numax*min(nwdr,1),50,72,50,72,ir)
!            CALL  extobs3gap(par,numax*min(nwdr,1),82,122,82,122,ir)
!            CALL  extobs3gap(par,numax*min(nwdr,1),92,132,92,132,ir)
!            CALL  extobs3gap(par,numax*min(nwdr,1),14,20,14,20,ir)
!            CALL  extobs3gap(par,numax*min(nwdr,1),34,50,34,50,ir)
!            CALL  extobs3gap(par,numax*min(nwdr,1),56,82,56,82,ir)
!            CALL  extobs3gap(par,numax*min(nwdr,1),88,126,88,126,ir)
!            CALL flush(2)
!          END DO
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),20,20,20,16)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),20,28,20,20)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),20,28,20,24)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),20,32,20,28)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),28,40,28,30)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),30,42,30,32)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),30,50,30,32)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),50,82,50,74)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),50,82,50,80)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),50,84,50,82)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),82,126,82,120)
!          CALL  extobsiso(par,yfit,numax*min(nwdr,1),82,132,82,126)
          END IF
!          DO nn=16,32,2     ! 44
!            CALL  extobsiso(par,yfit,numax*min(nwdr,1),20,20,20,nn)
!          END DO
!          DO nn=16,32,2     ! 44
!            CALL  extobsiso(par,yfit,numax*min(nwdr,1),20,28,20,nn)
!          END DO
!          IF(2 < 1) THEN
!          DO nn=62,88,2     ! 74
!            CALL  extobsiso(par,yfit,numax*min(nwdr,1),50,82,50,nn)
!          END DO
!          DO nn=118,132,2   ! 132
!            CALL  extobsiso(par,yfit,numax*min(nwdr,1),82,126,82,nn)
!          END DO
!          DO np=14,18,2     ! 16
!            CALL  extobsiso(par,yfit,numax*min(nwdr,1),20,20,np,20)
!          END DO
!          DO np=36,42,2     ! ?
!            CALL  extobsiso(par,yfit,numax*min(nwdr,1),40,50,np,50)
!          END DO
!          DO np=52,70,2     ! 62
!            CALL  extobsiso(par,yfit,numax*min(nwdr,1),50,82,np,82)
!          END DO
!          DO np=80,88,2     ! 88
!            CALL  extobsiso(par,yfit,numax*min(nwdr,1),82,126,np,126)
!          END DO
!          END IF
          close(82)
!          STOP '3gap done'
!          WRITE(2,*) ' Isotopic shifts finished. Start 3gap:'
!
          CALL  extobs5gap(par,numax*min(nwdr,1),20,24,20,24,0)
          CALL  extobs5gap(par,numax*min(nwdr,1),50,72,50,72,0)
          CALL  extobs5gap(par,numax*min(nwdr,1),82,122,82,122,0)
!          CALL  extobs5gap(par,numax*min(nwdr,1),92,132,92,132,0)
          CALL  extobs5gap(par,numax*min(nwdr,1),34,50,34,50,0)
          CALL  extobs5gap(par,numax*min(nwdr,1),56,82,56,82,0)
          CALL  extobs5gap(par,numax*min(nwdr,1),88,126,88,126,0)
!
          CALL  extobssgap(par,numax*min(nwdr,1),20,24,20,24)
          CALL  extobssgap(par,numax*min(nwdr,1),50,72,50,72)
          CALL  extobssgap(par,numax*min(nwdr,1),82,122,82,122)
!          CALL  extobssgap(par,numax*min(nwdr,1),92,132,92,132)
!          CALL  extobssgap(par,numax*min(nwdr,1),14,20,18,20)
          CALL  extobssgap(par,numax*min(nwdr,1),34,50,36,50)
          CALL  extobssgap(par,numax*min(nwdr,1),56,82,58,82)
          CALL  extobssgap(par,numax*min(nwdr,1),88,126,90,126)
!          STOP '5gap done'
!          END IF
!          CALL  extobs3gap(par,numax*min(nwdr,1),50,66,50,76,1)
!          CALL  extobs3gap(par,numax*min(nwdr,1),82,122,82,122,1)
!          CALL  extobs3gap(par,numax*min(nwdr,1),92,130,92,132,1)
!          CALL  extobs3gap(par,numax*min(nwdr,1),20,24,20,24,1)
!          CALL  extobs3gap(par,numax*min(nwdr,1),28,36,28,36,1)
!          CALL  extobs3gap(par,numax*min(nwdr,1),50,66,50,66,1)
!          CALL  extobs3gap(par,numax*min(nwdr,1),50,74,50,74,1)
!          CALL  extobs3gap(par,numax*min(nwdr,1),82,122,82,122,1)
!          CALL  extobs3gap(par,numax*min(nwdr,1),92,130,92,132,1)
!          CALL  extobs3gap(par,numax*min(nwdr,1),20,22,20,32,0)
!          CALL  extobs3gap(par,numax*min(nwdr,1),28,28,28,36,0)
!          CALL  extobs3gap(par,numax*min(nwdr,1),50,56,50,80,0)
!          CALL  extobs3gap(par,numax*min(nwdr,1),82,120,82,124,0)
!          CALL  extobssgap(par,numax*min(nwdr,1),20,22,20,32)
!          CALL  extobssgap(par,numax*min(nwdr,1),50,56,50,80)
!          CALL  extobssgap(par,numax*min(nwdr,1),82,120,82,124)
!          CALL  extobssgap(par,numax*min(nwdr,1),14,20,22,20)
!          CALL  extobssgap(par,numax*min(nwdr,1),34,50,44,50)
!          CALL  extobssgap(par,numax*min(nwdr,1),56,82,60,82)
!          CALL  extobs3gap(par,numax*min(nwdr,1),20,22,20,32,1)
!          CALL  extobs3gap(par,numax*min(nwdr,1),28,28,28,36,1)
!          CALL  extobs3gap(par,numax*min(nwdr,1),50,56,50,80,1)
!          CALL  extobs3gap(par,numax*min(nwdr,1),82,120,82,124,1)
          close(83)
!           DO nn=22,26,2
!             CALL extobs_avgap(par,workdv,nparmx,20,nn,2)
!           END DO
!           DO nn=62,80,2
!             CALL extobs_avgap(par,workdv,nparmx,50,nn,2)
!           END DO
!           DO nn=122,130,2
!             IF(nn.NE.126) 
!     &         CALL extobs_avgap(par,workdv,nparmx,50,nn,2)
!           END DO
!           DO np=14,18,2
!             CALL extobs_avgap(par,workdv,nparmx,np,20,1)
!           END DO
!           DO np=34,44,2
!             CALL extobs_avgap(par,workdv,nparmx,np,50,1)
!           END DO
!           DO np=52,66,2
!             CALL extobs_avgap(par,workdv,nparmx,np,82,1)
!           END DO
          nparact=numax*min(nwdr,1)
!          CALL extobs_centrdep(par,nparact,20,16,20,34,3)    ! rms radii
!          CALL extobs_centrdep(par,nparact,14,14,14,14,1)
!          CALL extobs_centrdep(par,nparact,20,20,20,20,1)
!          CALL extobs_centrdep(par,nparact,20,28,20,28,1)
!          CALL extobs_centrdep(par,nparact,28,30,28,30,1)
!          CALL extobs_centrdep(par,nparact,28,32,28,32,1)
!          CALL extobs_centrdep(par,nparact,28,34,28,34,1)
!          CALL extobs_centrdep(par,nparact,28,36,28,36,1)
!          CALL extobs_centrdep(par,nparact,40,50,40,50,1)
!          CALL extobs_centrdep(par,nparact,40,52,40,52,1)
!          CALL extobs_centrdep(par,nparact,40,54,40,54,1)
!          CALL extobs_centrdep(par,nparact,50,66,50,66,1)
!          CALL extobs_centrdep(par,nparact,50,74,50,74,1)
!          CALL extobs_centrdep(par,nparact,62,86,62,86,1)
!          CALL extobs_centrdep(par,nparact,82,122,82,122,1)
!          CALL extobs_centrdep(par,nparact,82,124,82,126,1)
!          CALL extobs_centrdep(par,nparact,82,126,82,126,1)
!          CALL extobs_centrdep(par,nparact,14,14,14,14,3)
!          CALL extobs_centrdep(par,nparact,20,20,20,20,3)
!          CALL extobs_centrdep(par,nparact,20,28,20,28,3)
!          CALL extobs_centrdep(par,nparact,28,30,28,30,3)
!          CALL extobs_centrdep(par,nparact,28,32,28,32,3)
!          CALL extobs_centrdep(par,nparact,28,34,28,34,3)
!          CALL extobs_centrdep(par,nparact,28,36,28,36,3)
!          CALL extobs_centrdep(par,nparact,40,50,40,50,3)
!          CALL extobs_centrdep(par,nparact,40,52,40,52,3)
!          CALL extobs_centrdep(par,nparact,40,54,40,54,3)
!          CALL extobs_centrdep(par,nparact,50,66,50,66,3)
!          CALL extobs_centrdep(par,nparact,50,74,50,74,3)
!          CALL extobs_centrdep(par,nparact,62,86,62,86,3)
!          CALL extobs_centrdep(par,nparact,82,122,82,122,3)
!          CALL extobs_centrdep(par,nparact,82,124,82,126,3)
!          CALL extobs_centrdep(par,nparact,82,126,82,126,3)
!          CALL extobs_centrdep(par,nparact,50,82,66,82,1)
!          CALL extobs_centrdep(par,nparact,76,126,94,126,1)
!          CALL extobs_centrdep(par,nparact,94,184,122,184,1)
!          CALL extobs_centrdep(par,nparact,20,28,20,28,2)
!          CALL extobs_centrdep(par,nparact,50,82,66,82,2)
!          CALL extobs_centrdep(par,nparact,76,126,94,126,2)
!          CALL extobs_centrdep(par,nparact,94,184,122,184,2)
!          CALL extobs_centrdep(par,nparact,50,82,66,82,3)
!          CALL extobs_centrdep(par,nparact,76,126,94,126,3)
!          CALL extobs_centrdep(par,nparact,94,184,122,184,3)
        END DO
        END IF
      endif
      if(iextan.eq.3) call extobschain(par,yfit,numax*min(nwdr,1))
      END IF
c
    5 continue
      STOP ' normal exit'
 991  STOP ' error in reading FOR009'
c
c
      end

c-----naming----------------------------------------------------------
c
      subroutine naming(zuord,ioptsk)
c
c     genererates array 'zuord' which carries the names of the
c     associated Skyrme parameters
c
c      zuord    = output array of names
c      iftinp   = input variant in 'for001'
c      ioptsk   = option for Skyrme force
c
      include 'params.inc'
      include 'constn.inc'
c      parameter (zero=0D0)
c
      include 'comfit.inc'
c
      character*6 zuord(kparmx)
      character*6 namepool(17,0:4)
      character*6 optpool(11:38,0:11)
      data ((namepool(i,j),i=1,17),j=0,4)
     &    /'    b0','    b1','    b2','    b3','    t4',
     &     '   b0p','   b1p','   b2p','   b3p',' power',
     &     '  byuk',' byukp','  ayuk','exmass','   b4p',
     &     'pairfp','pairfn',
     &     '    t0','    t1','    t2','    t3','    t4',
     *     '    x0','    x1','    x2','    x3',' power',
     *     '  tyuk','  xyuk','  ayuk','exmass','   b4p',
     &     'pairfp','pairfn',
     &     '    c0','    b1','    b2','    c3','    t4',
     &     '   c0p','   b1p','   b2p','   c3p',' power',
     &     '  byuk',' byukp','  ayuk','exmass','   b4p',
     &     'pairfp','pairfn',
     &     '    c0','    b1','    b2','    c3','    t4',
     &     '   c0p','   b1p','   b2p','   c3p',' power',
     &     '  byuk',' byukp','  ayuk','exmass','   b4p',
     &     'pairfp','pairfn',
     &     'EoverA','effmas','c0delr','rho_nm','    t4',
     &     ' asymm',' kappa','c1delr','dasymm',' compr',
     &     '  byuk',' byukp','  ayuk','exmass','   b4p',
     &     'pairfp','pairfn'/
      data ((optpool(i,j),i=11,14),j=5,6)
     &    /'    c1','   c1p','    c4','   c4p',
     &     '    c1','   c1p','    c4','   c4p'/
      data ((optpool(i,j),i=11,14),j=8,8)
     &    /'b3tild','   b12','  b12p',' undef'/
      data ((optpool(i,j),i=11,14),j=9,9)
     &    /'    b7','   b7p','  pow7',' undef'/
      data ((optpool(i,j),i=11,14),j=10,10)
     &    /' riftm','rho0pr','powerp',' undef'/
      data ((optpool(i,j),i=14,14),j=11,11)
     &    /' riftm'/
      data ((optpool(i,j),i=18,38),j=11,11)
     &    /'rho0pr','powerp','b3tild','   b12','  b12p',
     &     '    b7','   b7p','  pow7',
     &     '    c1','   c1p','    c4','   c4p','ecstab',
     &     'alpinv','betinv','alpin2','betin2','riftm1',
     &     '   tte','   tto','str_sq'/
c
c----------------------------------------------------------------------
c
      zuord = '      '
      if(kparmx.lt.17)
     &  stop ' Fit code requires KPARMX >= 17'
      if(ioptsk.eq.11 .and. kparmx.lt.30)
     &  stop ' IOPTSK=11 requires KPARMX >= 30'
      IF(ioptsk==13) THEN
        CALL naming_fy(zuord)
        RETURN
      END IF
c
      izuord = iftinp
      IF(iftinp==5) izuord = 4
      do ipar=1,17
        zuord(ipar) = namepool(ipar,izuord)
      enddo
      if(ioptsk.eq.5 .or. ioptsk.eq.6 .or. ioptsk.eq.8) then
        do ipar=11,14
          zuord(ipar) = optpool(ipar,ioptsk)
         enddo
      elseif(ioptsk.eq.9 .or. ioptsk.eq.10) then
        do ipar=11,14
          zuord(ipar) = optpool(ipar,ioptsk)
        enddo
      elseif(ioptsk.eq.11) then
        zuord(14) = optpool(14,ioptsk)
        do ipar=18,kparmx
          zuord(ipar) = optpool(ipar,ioptsk)
        enddo
      elseif(ioptsk.eq.12) then
        zuord(1)  = 'Gsigma'
        zuord(2)  = 'Gomega'
        zuord(3)  = 'Grho  '
        zuord(4)  = 'iscalm'
        zuord(6)  = 'Msigma'
        zuord(7)  = 'Momega'
        zuord(8)  = 'Mrho  '
        zuord(9)  = 'ivectm'
        zuord(10) = 'dpol  '
        zuord(18) = 'rho0pr'
      endif
c
      return
      end

c-----attach--------------------------------------------------------------------
c
      subroutine attach(par,jpact)
c
      USE NMP_map
c
c     the explicit model parameters ,
c       t0,t1,...
c     are attached according to the values on the field 'par'
c     and the pointers on the field 'jp' .
c
      include 'params.inc'
c
      include 'skypar.inc'
      include 'numpar.inc'
      include 'comfit.inc'
c
      dimension par(kparmx),jpact(kparmx)
c
c     for mass-recoupling according to 'recoupmass.tex':
      data betatau2/7.23505631D0/  ! (3*pi**2/2)**(2/3)*(6/5)
c
      data zero/0D0/
c
c--------------------------------------------------------------------
c
ctest      write(6,*) ' ATTACH with IFTINP:',iftinp
      if(ioptsk.eq.11 .and. kparmx.lt.30)
     &  stop ' IOPTSK=11 requires KPARMX >= 30'
      IF(ioptsk==13) THEN
        pairfp=0D0
        pairfn=0D0
        CALL attach_fy(par,jpact)
        RETURN
      END IF
      IF(iftinp>4) STOP "ATTACH: IFTINP too large for Skyrme branch"
c
c      prepare choice of rhonm for transformation C0 C3
c
c      if(iftinp.eq.2 .or. iftinp.eq.3) then
      if(iftinp.eq.2) then
        rhonm = 0.159
      else
        rhonm = 0.160
      endif
c
c     set defaults
c
      powerp = 0.0D0
      rhonm7 = rhonm
      riftm1 = 0D0
      str_squark = -0.24D0
!      pairgm = 1D0
c
        t5     = zero
        x5     = zero
        t6     = zero
        x6     = zero
         b5     = 0D0
         b6     = 0D0
         b6p    = 0D0
         b5p    = 0D0
         tyuk   = 0D0
         ayuk   = 0D0
         byukp  = 0D0
         csurf  = 0D0
c
c     first deal with the NMP entry
c
      IF(iftinp==4) THEN
c
       eoain   = par(jpact(1))
       rhoin   = par(jpact(4))
       compin = par(jpact(10))
       asymin = par(jpact(6))
       dasyin = par(jpact(9))
       efmain = par(jpact(2))
       skapin = par(jpact(7))
       c0delrhoin = par(jpact(3))
       c1delrhoin = par(jpact(8))
       IF(ifdfms.eq.0) THEN
         amfc0(1) = 20.7525
         amfc0(2) = 20.7525
       ELSE IF(ifdfms.eq.2) THEN
         amfc0(1) = 20.73553 ! 20.7339829978 !
         amfc0(2) = 20.73553 ! 20.7339829978 !
       ELSE IF(ifdfms.eq.3) THEN
         amfc0(1) = 20.749811 ! 20.74982   ! 20.7505   !20.747
         amfc0(2) = 20.721249 ! 20.72126   ! 20.7219   !20.719
       ELSE
         STOP "IFDFMS must be 1,2 or 3 for NMP entry"
       END IF
       CALL NMP_Skyrme(t0,t1,t2,t3,x0,x1,x2,x3,power,
     &                 rhoin,eoain,compin,efmain,
     &                 asymin,dasyin,skapin,
     &                 amfc0(1),amfc0(2),c0delrhoin,c1delrhoin)
c
       b2=c0delrhoin          ! only b2 & b2p together correct
       b2p=c1delrhoin
       t4 = par(jpact(5))
       b4p = par(jpact(15))
       pairfp = par(jpact(16))
       pairfn = par(jpact(17))
       IF(ioptsk==0 .OR. ioptsk==10) THEN
         riftm  = par(jpact(11))
         rho0pr = par(jpact(12))
         powerp = par(jpact(13))
         IF(kparmx>29) ecstab = par(jpact(30))
         IF(kparmx>30) hgradxi = par(jpact(31))
       ELSE IF(ioptsk==11) THEN
         tyuk   = par(jpact(11))
         xyuk   = par(jpact(12))
         ayuk   = par(jpact(13))
         riftm  = par(jpact(14))
         rho0pr = par(jpact(18))      ! default 0.16
         powerp = par(jpact(19))
         b3tilde = par(jpact(20))
         b12    = par(jpact(21))
         b12p   = par(jpact(22))
         pow7   = par(jpact(25))
         b7     = par(jpact(23))/rhonm**pow7
         b7p    = par(jpact(24))/rhonm**pow7
         c1     = par(jpact(26))
         c1p    = par(jpact(27))
         c4     = par(jpact(28))
         c4p    = par(jpact(29))
         ecstab = par(jpact(30))
         alpinv = par(jpact(31))
         betinv = par(jpact(32))
         if(kparmx.ge.34) then
           alpin2 = par(jpact(33))
           betin2 = par(jpact(34))
         endif
         if(kparmx.ge.37) then
           riftm1 = par(jpact(35))
           if(riftm1==0D0) riftm1 = riftm
!           WRITE(2,*) ' RIFTM1=',riftm1
           tte    = par(jpact(36))       ! even tensor coupling
           tto    = par(jpact(37))       ! odd tensor coupling
         endif
         IF(kparmx.ge.38) THEN
!           str_squark = par(jpact(38))
           IF(par(jpact(38))>1D-2) pairgm = par(jpact(38))
         END IF
       ELSE
         STOP " NMP entry needs  IOPSK=1,10, or 11"
       END IF
       RETURN
      END IF
c
c     now all the rest
c    
c
      IF(ioptsk==12) THEN
        IF(iftinp.NE.0) STOP 'ATTACH: QMC requires IFTINP=0'
        CALL attach_QMC(par,jpact)
        WRITE(6,'(a,3(1pg13.5))') 'PAIRA:',pairfp,pairfn,rho0pr
        pairfp = par(jpact(16))
        pairfn = par(jpact(17))
        rho0pr = par(jpact(18))       ! default 0.16
        RETURN
      END IF
c
      if(iftinp.eq.1) then
       t0     = par(jpact(1))
       t1     = par(jpact(2))
       t2     = par(jpact(3))
       t3     = par(jpact(4))
       t4     = par(jpact(5))
       x0     = par(jpact(6))
       x1     = par(jpact(7))
       x2     = par(jpact(8))
       x3     = par(jpact(9))
       power  = par(jpact(10))
!       rho0pr= par(jpact(12))
!       xyuk=0D0
       if(ioptsk.eq.1) then
         b6     = par(jpact(11))
         b6p    = par(jpact(12))
         b5p    = par(jpact(13))
       elseif(ioptsk.eq.5 .or. ioptsk.eq.6) then
         c1     = par(jpact(11))
         c1p    = par(jpact(12))
         c4     = par(jpact(13))
         c4p    = par(jpact(14))
       elseif(ioptsk.eq.8) then
         b3tilde = par(jpact(11))
         b12     = par(jpact(12))
         b12p    = par(jpact(13))
       elseif(ioptsk.eq.9) then
         pow7   = par(jpact(13))
         b7     = par(jpact(11))/rhonm**pow7
         b7p    = par(jpact(12))/rhonm**pow7
       elseif(ioptsk.eq.10) then
         riftm  = par(jpact(11))
         rho0pr = par(jpact(12))
         powerp = par(jpact(13))
         IF(kparmx>29) ecstab = par(jpact(30))
         IF(kparmx>30) hgradxi = par(jpact(31))
       elseif(ioptsk.eq.11) then
         tyuk   = par(jpact(11))
         xyuk   = par(jpact(12))
         ayuk   = par(jpact(13))
         riftm  = par(jpact(14))
         rho0pr = par(jpact(18))      ! default 0.16
         powerp = par(jpact(19))
         b3tilde = par(jpact(20))
         b12    = par(jpact(21))
         b12p   = par(jpact(22))
         pow7   = par(jpact(25))
         b7     = par(jpact(23))/rhonm**pow7
         b7p    = par(jpact(24))/rhonm**pow7
         c1     = par(jpact(26))
         c1p    = par(jpact(27))
         c4     = par(jpact(28))
         c4p    = par(jpact(29))
         ecstab = par(jpact(30))
         alpinv = par(jpact(31))
         betinv = par(jpact(32))
         if(kparmx.ge.34) then
           alpin2 = par(jpact(33))
           betin2 = par(jpact(34))
         endif
         if(kparmx.ge.37) then
           riftm1 = par(jpact(35))
           if(riftm1==0D0) riftm1 = riftm
!           WRITE(2,*) ' RIFTM1=',riftm1
           tte    = par(jpact(36))       ! even tensor coupling
           tto    = par(jpact(37))       ! odd tensor coupling
         endif
         IF(kparmx.ge.38) THEN
!           str_squark = par(jpact(38))
           IF(par(jpact(38))>1D-2) pairgm = par(jpact(38))
         END IF
       else
         tyuk   = par(jpact(11))
         xyuk   = par(jpact(12))
         ayuk   = par(jpact(13))
       endif
       exmass = par(jpact(14))        ! switch to exchange functional
       rhoswm = 0D0                   ! tanh switching for nonl. mass disabled
       csurf  = 0D0
cold       csurf  = par(jpact(14))
       b4p    = par(jpact(15))
c
c       if(b4p.eq.zero) b4p = 1.0e-20
c        t9     = par(jpact(15))
       if(csurf.ne.zero)
     &  stop  ' option csurf.ne.0 invalid with iftinp.eq.1'
      else
       b4p    = par(jpact(15))
       b1     = par(jpact(2))
       b2     = par(jpact(3))
       if(iftinp.eq.0) then
        b0     = par(jpact(1))
        b3     = par(jpact(4))
        b0p    = par(jpact(6))
        b3p    = par(jpact(9))
        power  = par(jpact(10))
       elseif(iftinp.eq.2) then
        c0     = par(jpact(1))
        c3     = par(jpact(4))
        c0p    = par(jpact(6))
        c3p    = par(jpact(9))
        power  = par(jpact(10))
        powfac = 3.0/rhonm**power
        b0     = 0.5*(c0+c3)
        b3     = powfac*(c0-c3)
        b0p    = 0.5*(c0p+c3p)
        b3p    = powfac*(c0p-c3p)
       elseif(iftinp.eq.3) then
        d0     = par(jpact(1))
        d3     = par(jpact(4))
        power  = par(jpact(10))
        powfac = rhonm**0.66666666666666667D0
        c0     = d0 - b1*betatau2*powfac
        c3     = d3 - b1*betatau2*powfac
     &                *(1.0D0-20.0D0/(9.0D0*power*(power+1.0D0)))
c        write(*,*) ' c0,c3=',c0,c3
        c0p    = par(jpact(6))
        c3p    = par(jpact(9))
        powfac = 3.0/rhonm**power
        b0     = 0.5*(c0+c3)
        b3     = powfac*(c0-c3)
        b0p    = 0.5*(c0p+c3p)
        b3p    = powfac*(c0p-c3p)
!        write(2,'(a,4(1pg16.8))') ' b0...=',b0,b3,b0p,b3p
       else
         stop " this IFTINP not yet provided"
       endif
       if(iftinp.eq.3) then
         b4sum  = par(jpact(5))
         if(b4p.ne.ZERO) then
           t4     = b4sum - par(jpact(15))
         else
           t4     = b4sum*2.0D0/3.0D0
         endif
       else
         t4     = par(jpact(5))
       endif
       b1p    = par(jpact(7))
       b2p    = par(jpact(8))
       if(ioptsk.eq.1) then
         b6     = par(jpact(11))
         b6p    = par(jpact(12))
         b5p    = par(jpact(13))
       elseif(ioptsk.eq.5 .or. ioptsk.eq.6) then
         c1     = par(jpact(11))
         c1p    = par(jpact(12))
         c4     = par(jpact(13))
         c4p    = par(jpact(14))
       elseif(ioptsk.eq.8) then
         b3tilde = par(jpact(11))
         b12     = par(jpact(12))
         b12p    = par(jpact(13))
       elseif(ioptsk.eq.9) then
         pow7   = par(jpact(13))
         b7     = par(jpact(11))/rhonm**pow7
         b7p    = par(jpact(12))/rhonm**pow7
       elseif(ioptsk.eq.10) then
         riftm  = par(jpact(11))
         rho0pr = par(jpact(12))      ! default 0.16
         powerp = par(jpact(13))
         IF(kparmx>29) ecstab = par(jpact(30))
         IF(kparmx>30) hgradxi = par(jpact(31))
       elseif(ioptsk.eq.11) then
         byuk   = par(jpact(11))
         byukp   = par(jpact(12))
         ayuk   = par(jpact(13))
         riftm  = par(jpact(14))
         rho0pr = par(jpact(18))      ! default 0.16
         powerp = par(jpact(19))
         b3tilde = par(jpact(20))
         b12    = par(jpact(21))
         b12p   = par(jpact(22))
         pow7   = par(jpact(25))
         b7     = par(jpact(23))/rhonm**pow7
         b7p    = par(jpact(24))/rhonm**pow7
         c1     = par(jpact(26))
         c1p    = par(jpact(27))
         c4     = par(jpact(28))
         c4p    = par(jpact(29))
         ecstab = par(jpact(30))
         alpinv = par(jpact(31))
         betinv = par(jpact(32))
         if(kparmx.ge.34) then
           IF(jpact(33)>0) alpin2 = par(jpact(33))
           IF(jpact(34)>0) betin2 = par(jpact(34))
         endif
         if(kparmx.ge.37) then
           IF(jpact(35)>0) riftm1 = par(jpact(35))
           if(riftm1==0D0) riftm1 = riftm
!           WRITE(2,*) ' RIFTM1=',riftm1
           IF(jpact(36)>0) tte    = par(jpact(36))       ! even tensor coupling
           IF(jpact(37)>0) tto    = par(jpact(37))       ! odd tensor coupling
         endif
         IF(kparmx.ge.38) THEN
!           IF(jpact(38)>0) str_squark = par(jpact(38))
           IF(par(jpact(38))>1D-2) pairgm = par(jpact(38))
         END IF
       else
         byuk   = par(jpact(11))
         byukp  = par(jpact(12))
         ayuk   = par(jpact(13))
       endif
cnew       csurf  = 0D0
       rhoswm = 0D0                   ! tanh switching for nonl. mass disabled
       if(ioptsk.eq.1) then
         csurf  = par(jpact(14))      ! 'csurf' neede for this extended nonl.
         exmass = 0D0
       elseif(ioptsk.eq.2 .or. ioptsk.eq.3) then
         exmass = par(jpact(14))      ! switch to exchange functional
         csurf  = 0D0
       elseif(ioptsk.eq.4) then
         b1kin  = b1
         b1     = zero
         b5kin  = par(jpact(14))
       endif
       b4p    = par(jpact(15))
c
c       if(b4p.eq.zero) b4p = 1.0e-20
c       t9     = par(jpact(15))
cs       if(csurf.ne.zero) then
cs        c1     = b1
cs        b1     = csurf-sqrt(csurf*(csurf-2.0*c1))
cs        b5     = (c1-b1)/rhonm
cs        b5p    = zero
cs       endif
c               evaluate optionally b2 and b2p
       if(ib2fix.eq.1) then
        b2     = b1*0.5
        b2p    = b1p*0.5
       endif
c              recouple parameters to t0 etc
       bs0    = b0+0.5*b0p
       bs1    = b1+0.5*b1p
       bs2    = b2+0.5*b2p
       bs3    = b3+0.5*b3p
       fourtr = 4.0/3.0
       t0     = fourtr*(bs0-0.5*b0p)
       if(t0.ne.zero) then
         x0     = b0p/t0-0.5
       else
         x0     = zero
       endif
       t3     = fourtr*(bs3-0.5*b3p)
       x3     = b3p/t3-0.5
       a1     = bs1+bs2+bs2
       a1p    = b1p+b2p+b2p
       a2     = 3.0*bs1-2.0*bs2
       a2p    = -3.0*b1p+2.0*b2p
       t1     = fourtr*(a1-0.5*a1p)
       x1     = a1p/t1-0.5
       t2     = fourtr*(a2-0.5*a2p)
       x2     = a2p/t2-0.5
       if(ioptsk.ne.1) then
         bsyuk  = byuk + 0.5*byukp
         tyuk   = fourtr*(bsyuk-0.5*byukp)
         if(tyuk.ne.zero) then
          xyuk   = byukp/tyuk-0.5
         else
          xyuk   = zero
         endif
       endif
       if(csurf.ne.zero.and.ioptsk.ne.1) then
        b5p    = zero
        b5     = csurf                           ! misuse entry 'csurf'
        bs5    = b5+0.5*b5p
        t5     = 6.0*fourtr*(bs5-0.5*b5p)
        x5     = b5p/t5-0.5
        t6     = 6.0*fourtr*(3.0*bs5+1.5*b5p)
        x6     = -3.0*b5p/t6-0.5
       else
        t5     = zero
        x5     = zero
        t6     = zero
        x6     = zero
       endif
      endif
c
      if(t0.eq.zero .and. x0.eq.zero) then
       t0=tyuk
       x0=xyuk
      endif
c
c     set actual folding model
c
c      write(2,*)
c     &  ' byuk,byukp,ayuk=',byuk,byukp,ayuk
c      write(2,*)
c     &  ' tyuk,xyuk,ayuk=',tyuk,xyuk,ayuk
      ayukp  = 0D0               ! ayuk     
      ayukm  = ayuk
      ayukls = ayuk
c
      if(kparmx.ge.17) then
       pairfp = par(jpact(16))
       pairfn = par(jpact(17))
       IF(pairfn<-1D3) pairfn = pairfp
      endif
c
c     map to Fayans pairing parameters
c
      IF(hgradxi.NE.0D0) THEN
        fxi = pairfp
        rho0xi = 0.16D0
        rho0fy = 0.08D0
        h1xi = pairfn
        gammaxi = rho0pr
!        WRITE(2,*) ' Fayans: ',fxi,h1xi,gammaxi
      END IF
c
!      WRITE(2,*) 'b:',b0,b1,b2,b3,t4
!      WRITE(2,*) 'bp:',b0p,b1p,b2p,b3p,b4p
!      WRITE(2,*) 'power:',power

      IF(ABS(tte)+ABS(tto)>0D0 .AND. ABS(riftm1)>0D0) THEN
        STOP ' active TTE&TTO not compatible with active RIFTM1'
      END IF
c
      return
      end

c-----beving ---------------------------------------------------------
c
      subroutine beving(zuord,a,sigmaa,sigma0,numax,nwdr,ifcurv)
c
c     driver routine for bevington step complemented by a
c     stochastic change of parameters if step is frozen.
c
c     bevington step is iteration with gradient in chi**2-landscape
c     with a simple estimate for the curvature.
c     the step is complmented by a curved step.
c
c
c      the parameters in the calling sequence are:
c       nwdr    = maximum number of iterations
c       iwarn   = returns the number of attempts to find a lower
c                 chi**2 . iwarn=1 means a direct succes. iwarn=99
c                 means that no lower chi**2 has been found.
c       iflast  = 1 to indicate the final call which produces a
c                   detailed output, 0 in all other cases.
c       ifcurv  = number of normal bevington steps which are to be
c                 run before the curved step is switched on. if no
c                 curved step is desired make ifcurve larger than
c                 nwdr.
c       zuord   = names of parameters for printing of protocol 'for003'
c       a       = initial parameters as input, new parameters after
c                 one iteration step given back as output.
c       sigmaa   = field containing the correlated errors on the
c                  model parameters, it is a 1-dimensional field
c                  sigmaa(parameter).
c       sigma0   = field containing the uncorrelated errors.
c
c
c     dimension statement valid for nterms up to 'kparmx'
c
      include 'params.inc'
      include 'comfit.inc'
      parameter (kparsv=3*kparmx)
c
      character*6 zuord(kparmx)
      dimension a(kparmx),sigmaa(kparmx),sigma0(kparmx)
      dimension amin(kparmx)
c
cold      include 'numerc.inc'
cold      include 'shells.inc'
cold      include 'forces.inc'
      include 'result.inc'
c
      data istchm,iseed/3,19067/
c
c--------------------------------------------------------------------
c
      istoch=1
      ifalt2=0
      chisqr = 1.0D30
      chisqm = chisqr                    !  local minimum in this call
      amin(1:nparin) = a(1:nparin)
c
      do itera=1,nwdr
        chia=chisqr
        if(itera.eq.1.and.nwdr.lt.2) then
         iwrte=1
        else
         iwrte=-2
        endif
c
        write(06,'(/a,i3,a,g12.5)')
     &   ' iteration nr.',itera,' , chi**2=',chia
!        iflast=0
!        IF(itera==nwdr) iflast=1
        call xd1(a,sigmaa,sigma0,numax,nwdr,iwarn,0,
     &           ifcurv,zuord)
c
        if(chisqr.lt.chisqm) then
          write(6,'(a,2g12.4)') 
     &     ' update local minimum chi old/new=',chisqm,chisqr
          write(2,'(a,2g12.4)') 
     &     ' update local minimum chi old/new=',chisqm,chisqr
          chisqm = chisqr
          amin(1:nparin) = a(1:nparin)
          write(6,'(a,2g12.4)') 
     &     ' update FOR003. chi old/new=',chimin,chisqr
          chimin = chisqr
          rewind(3)
          do j=1,nparin
           write(3,'(1pg18.10,a,1pg12.4,a,i2,a,a)')  
     &        a(j),',',dparmo(j),',',iparva(j),',',zuord(iparva(j))
          enddo
          write(3,'(a/a,i4,1pg13.5)')
     &     '-11, "to stop"',' iteration,chi**2=',itera,chimin
          call flush(3)
        endif
c
c      iwarn=99 indicates that the gradient step run into
c      bad conditions. the parameters are changed stochastically
c      and the fit is continued.
c
        if(iwarn.eq.99 .and. itera.lt.nwdr) then
!        if(iwarn.eq.99) then
         if(istoch.gt.istchm) then
          write(02,'(/a)') ' number of stochastic attempts exhausted'
          goto 3601
         endif
         do j=1,numax
cold           a(j) = a(j) + (2*irand(istoch+j)-1)*dparmo(j)*freep
           a(j) = a(j) + (2.0D0*ran1(iseed)-1.0D0)*dparmo(j)
         enddo
         istoch = 1+istoch
         ifalt2 = 0
         chisqr=.1e30
         write(02,'(/a/)') ' iwarn=99. stochastic change of params.'
         write(02,'(a,6(/1x,4g14.6))')
     &     ' new parameters:',(a(j),j=1,numax)
        endif
c
c
        write(2,'(a,i4,10x,a,i4,10x,a,e10.3,10x,a,i3)') 
     &            '   nwdr = ',nwdr,'itera = ',itera,
     &            'flamda = ',flamda,'    iwarn = ',iwarn
        call flush(2)
c
        do i=1,numax
          dparmo(i)=sigma0(i)*freep
creset          if(dparmo(i).gt.dparm0(i)) dparmo(i)=dparm0(i)
        enddo
c         if((chia/chisqr).lt.1.001 .and. iwarn.lt.99) goto 3601
        endcp  = 1.0 + endcon
        if((chia/chisqr).lt.endcp .and. iwarn.lt.99) goto 3601 
        write(2,'(/2x,2f10.4,3i6/)') chia,chisqr,itera,nwdr
        call flush(2)
      enddo
c
 3601 continue
c
      if(chisqm.lt.chisqr) then
        write(2,'(a)')  
     &     ' reset to minimum chi*2: chi**2,chimin=',chisqr,chisqm
        write(6,'(a)')  
     &     ' reset to minimum chi*2: chi**2,chimin=',chisqr,chisqm
        chisqr = chisqm
        do i=1,nparin
          a(i) = amin(i)
        enddo
      endif
      call flush(2)
c
      return
      end




c-----xd1   ---------------------------------------------------------
c
      subroutine xd1(a,sigmaa,sigma0,nterms,nwdr,iwarn,iflast,ifcurv,
     &               zuord)
c
c     curfit from bevington , complemented by a curved step.
c     this version is adapted to the rmf-package.
c
c      the parameters in the calling sequence are:
c       nwdr    = maximum number of iterations; used here only to
c                 cut-off calculation of error matrix if nwdr.le.0  .
c       iwarn   = returns the number of attempts to find a lower
c                 chi**2 . iwarn=1 means a direct succes. iwarn=99
c                 means that no lower chi**2 has been found.
c       iflast  = 1 to indicate the final call which produces a
c                   detailed output, 0 in all other cases.
c       ifcurv  = number of normal bevington steps which are to be
c                 run before the curved step is switched on. if no
c                 curved step is desired make ifcurve larger than
c                 nwdr.
c       a       = initial parameters as input, new parameters after
c                 one iteration step given back as output.
c       yfitm0   = auxiliary field for intermediate saving and
c                  transfer of theoretical results to field yfit.
c       sigmaa   = field containing the correlated errors on the
c                  model parameters, it is a 1-dimensional field
c                  sigmaa(parameter).
c       sigma0   = field containing the uncorrelated errors.
c
c      some important parameters are transferred via common:
c       ifalt2   = in common /parart/ .  is equal to 0 if the
c                  routine is called the first time at the
c                  start of an iteration. is then returnes with
c                  ifalt2=1 which indicates that everything is
c                  properly initialised now.
c       ifsave   = switches the saving of hartree-fock information
c                  in the hafo-package. ifsave=1 means saving on
c                  file 19. it is switched on only for the "central"
c                  parameter sets and switched of during the
c                  variation of the parameters.
c       iwrte    = switches print level in the routine ffverg.
c       iwr      = switches print level in the hafo-package.
c       yfit     = field containing the theoretical results in the
c                  same manner as the experimental data on yexp.
c       yexp    = field containing the experimental data, it is
c                 2-dimensional field yexp(observable,nucleus).
c       sigmay  = field containing the choosen uncertainty on the
c                 experimental data.
c
c     dimension statement valid for nterms up to 'kparmx'
c
c
      include 'params.inc'
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'comfit.inc'
      parameter (kparsv=3*kparmx)
c
      logical tdiagz,tinvert
      parameter (tdiagz=.true.,tinvert=.TRUE.)
c
      character*6 zuord(kparmx)
      character*4 parbra
      character*2 parket
      character*16 errpri(2)
      dimension parsav(kparmx,3),chisav(3)
      dimension a(kparmx),ares(kparmx),sigmaa(kparmx),sigma0(kparmx)
      dimension alpha(kparmx,kparmx),b(kparmx),beta(kparmx)
      dimension vectz(kparmx,kparmx),rootz(kparmx)
      dimension alphdr(kparmx),betadr(kparmx)
      dimension parmom(kparmx),yfitmo(12,kdat),chis00(12)
c
c
      data ireduc,itrmax,ifrept
     *    /   0  ,   3  ,   1  /
      data errpri/'"sqb negative"  ','"too largeiiter"'/
c
c     these parameters are options for testing etc :
c       ireduc   = 0 --> damping of step by increasing 'flamda' as
c                        propsed by bevington
c                  1 --> damping of step by reducing it by factor 2
c       itrmax   = number of attempts with reduced step in curved step
c       ifrept   = 0 --> chi**2 and situation from previous final
c                        calculation is used
c                  1 --> new calculation at the beginning
c                  2 --> equilibrised step
c
      data parbra,parket/'par(',')='/
      data xlowlm,one/1.e-30,1.0/
      data parsav/kparsv*0D0/,chisav/3*0D0/,itera/0/
      data ifdrc/0/                         ! switch to direct 2. derivative
      LOGICAL, PARAMETER :: tonestep=.TRUE.  ! evaluation of chi**2 landscape
      LOGICAL, PARAMETER :: tpreigen=.TRUE.  ! print details eigenvetors
      LOGICAL, PARAMETER :: tsensitivity=.TRUE.
c
c--------------------------------------------------------------------
c
      if(ifalt2.eq.0) then
       flamda= 0.001
       itera = 1
      else
       itera = 1+itera
      endif
       write(02,'(a,i5,a,i5)') '0 xd1 entered with ifalt2=',
     *                   ifalt2,'. now itera=',itera
       write(02,'(3(a,i3))') ' running with the options: ireduc=',
     *            ireduc,'     ifrept=',ifrept
      call flush(2)
      nfree  = nptsg - nterms
      rfree  = nfree
cold      numax  = nterms
      chiin  = chisqr
      if(nfree.le.0 .and. iflast.ne.1) then
       write(2,'(2a)')
     * ' not enough data for fitting the parameters',
     *                  ' ====> stop '
       stop
      endif
c
c   berechnung der theoretischen werte mit anfangsparametern,
c   d.h. wenn ifalt2 = 0 , oder
c   ha-fo-rechnungen fuer einzelne kerne ausgeben bei letzter
c   rechnung, d.h. wenn iflast = 0   .
c
      if(iflast.eq.1 .or. ifalt2.eq.0 .or. ifrept.eq.1) then
       ifalt2=1
       if(iflast.eq.1) then
         iwr=-1
         ifplot=3
       else
!         iwr=-1 !-2
         iwr=-2
         ifplot=0
       endif
       iwrte=1
       ifsave=1
       write(06,'(a)') ' initial or final chi**2'
       chisu = ffverg(a)
       if(ierror.ne.0) then
        write(02,'(3a)') ' the error ',errpri(ierror),
     *     ' occured at the initial or final evaluation! stop.'
        if(nwdr.gt.0) stop
       endif
       iwrte=0
       iwr=-2
       if(iflast.ne.1) then
        iwr=-2
!        iwr=-1
        chiin=chisu
cold        chiini=chiin
        write(2,'(//a,e12.4,a//)') '  chi**2 input = ',chiin,
     *                             '  ***********'
        write(6,'(//a,e12.4,a//)') '  chi**2 input = ',chiin
       call flush(2)
       chisu0(1:ndats)=sumchi(1:ndats)
       chis00(1:ndats)=chisu0(1:ndats)
        write(2,'(15(1x,f7.2))') chisu0(1:ndats)
       endif
      endif
      if(nwdr.le.0) return
c
c   berechnung messwerte mit um dparmo geaendereten parametern
c   und damit ableitung nach den parametern :
c
c       vorbesetzungen
c
      b=a
      ares=a
      parmom=a
      if(ist1.gt.1) write(2,'(a)') '  variation der parameter'
      chires=chiin
      ireser=0
c
c        loop ueber die parameter, jeder wird der reihe nach
c        einmal geaendert.
c
      IF(tdiffobs) THEN
!        ALLOCATE(fjacobian(kderij,0:kparmx))
        ndatsd=ndats
        IF(maxdif>0) ndatsd=ndatsd+1 
        IF(if3gap>0) ndatsd=ndatsd+1 
        ndatsd=14
!        WRITE(2,'(a,i3)') 'FJACOBIAN-0: ndatsd=',ndatsd
        ilauf=0
        do idat=1,ndatsd
        do ipar=1,npts
          ilauf=ilauf+1
          fjacobian(ilauf,0)=chi(idat,ipar)                ! store base point
!          WRITE(2,'(3i3,1pg13.5)') idat,ipar,ilauf,fjacobian(ilauf,0)
        enddo
        enddo
      END IF
      yfitmo(1:ndats,1:npts) = yfit(1:ndats,1:npts)        ! save central values
      chi0sav(1:ndatsd,1:npts) = chi(1:ndatsd,1:npts)        ! save central values
      do 100 jpar=1,nterms
      parmom(jpar)=a(jpar)+dparmo(jpar)
      ifsave=0
      write(06,'(a,i4)') ' chi**2 for varied parameter nr.',jpar
      chisu = ffverg(parmom)
      write(06,'(a,g12.4)') ' chi**2=',chisu
      write(02,'(a,i4,a,g12.4,3a,1pg14.6)') ' parameter nr.',jpar,
     &     ' chi**2=',chisu,' param= ',zuord(iparva(jpar)),
     &     ' value=',parmom(jpar)
      call flush(2)
       if(ierror.ne.0) then
        write(02,'(3a)') ' the error ',errpri(ierror),
     *     ' occured at the parameter variation! stop.'
        write(02,'(a,2i5)') ' nucleus was: nprot,nneut=',nprot,nneut
        stop
       endif
c
c     reserve parameter abspeichern , wenn schon bei der
c     parametervariation ein kleineres chi**2 herausgekommen ist.
c
c
      chimom(jpar)=chisu
      if((chimom(jpar)/chires).le.0.999) then
        ares(1:nterms)=parmom(1:nterms)
        chires=chimom(jpar)
        ireser=jpar
      endif
c
c     opposite point in case of direct 2. derivatives
c
      if(ifdrc.eq.1) then
        parmom(jpar)=a(jpar)-dparmo(jpar)
        write(06,'(a,i4)') ' chi**2 for -varied parameter nr.',jpar
        chismn = ffverg(parmom)
        write(06,'(a,g12.4)') ' chi**2=',chismn
        write(02,'(a,i4,a,g12.4)')
     &        ' parameter nr.',-jpar,' chi**2=',chismn
        betadr(jpar) = chismn
        if((chismn/chires).le.0.999) then
          ares(1:nterms)=parmom(1:nterms)
          chires=chismn
          ireser=jpar
        endif
        alphdr(jpar) = (chimom(jpar)-2.0*chiin+chismn)/dparmo(jpar)**2
        betadr(jpar) = (chimom(jpar)-chismn)/(2.0*dparmo(jpar))
      endif
c
c     ableitung der messwerte nach den parametern. das resultat wird
c     auf dem feld 'fderij' gespeichert.
c
      ilauf=0
      icount=0
      do idat=1,ndats
      do ipar=1,npts
        ilauf=ilauf+1
        fderij(ilauf,jpar)=(yfitmo(idat,ipar)-yfit(idat,ipar))/
     *                     (a(jpar)-parmom(jpar))
        ftest=fderij(ilauf,jpar)/sigmay(idat,ipar)
!        IF(ABS(ftest)>1D-20) THEN
        IF(ABS(ftest)>1D20) THEN
          icount=1+icount
          WRITE(2,'(a,4i5,1pg13.5)') 'FDERIJ:',
     &     idat,ipar,icount,jpar,ftest
        END IF
      enddo
      enddo
      IF(tdiffobs) THEN
!        WRITE(2,'(a,i3)') 'FJACOBIAN-jpar: jpar=',jpar
        ilauf=0
        icount=0
        do idat=1,ndatsd
        do ipar=1,npts
          ilauf=ilauf+1
          fjacobian(ilauf,jpar)=
     &     -(chi(idat,ipar)-chi0sav(idat,ipar))/(a(jpar)-parmom(jpar))
!          IF(ABS(fjacobian(ilauf,jpar))>1D-20) THEN
          IF(ABS(fjacobian(ilauf,jpar))>1D20) THEN
            icount=icount+1
            WRITE(2,'(a,4i5,1pg13.5)') 'JACOBIAN:',
     &        idat,ipar,icount,jpar,fjacobian(ilauf,jpar)
          END IF
        enddo
        enddo
      END IF
c
c      zwischenausdrucke   ??: meaning not clear
c
      if(ist1.gt.1) then
       write(2,'(i3,8e14.6,e13.3/3x,9e14.6,e12.4)')
     *      jpar,(parmom(kl),kl=1,8),dparmo(jpar),
     *      (parmom(kl),kl=9,kparmx),chimom(jpar)
       write(2,'(15(1x,f7.2))') sumchi(1:ndats)
       sumchi(1:ndats)=(sumchi(1:ndats)-chisu0(1:ndats))
       write(2,'(15(1x,f7.2))') sumchi(1:ndats)
       sumchi(1:ndats)=sumchi(1:ndats)/dparmo(jpar)
       write(2,'(15(1x,f7.2))') sumchi(1:ndats)
      endif
c
      parmom(jpar)=a(jpar)
  100 continue
c
c
      if(ist1.gt.0) then
        write(2,'(a)') ' FDERIJ:'
        do i=1,ilauf
         write(2,'(1x,5g12.4)') fderij(i,1:nterms)
        enddo
      endif
      yfit(1:ndats,1:npts) = yfitmo(1:ndats,1:npts)   ! restore central values
c
      if(ireser.ne.0) write(2,'(a,2(/1x,10g13.6))')
     *             '   reserve parameter:',(ares(j),j=1,nterms)
c
c     second derivative of chi**2 on the field 'alpha' ,
c     the gradient of chi**2 on the field 'beta' .
c
      do j=1,nterms
        beta(j)=0D0
        do k=1,j
          alpha(j,k)=0D0
        enddo
      enddo
c
      DO j=1,nterms
c
        i=0
        do idat=1,ndats
        do ilauf=1,npts
          i=i+1
          beta(j)=beta(j)+(yexp(idat,ilauf)-yfitmo(idat,ilauf))*
     +            fderij(i,j)/sigmay(idat,ilauf)**2
        enddo
        enddo
c
        do k=1,j
          IF(tdiffobs) THEN
            i=0
            do idat=1,ndatsd
            do ilauf=1,npts
              i=i+1
!              WRITE(2,'(a,3i5)') 'ALPHA:',idat,ilauf,i
              alpha(j,k)=alpha(j,k)+fjacobian(i,j)*fjacobian(i,k)
            enddo
            enddo
          ELSE
            i=0
            do idat=1,ndats
            do ilauf=1,npts
              i=i+1
              alpha(j,k)=alpha(j,k)+fderij(i,j)*fderij(i,k)
     +                   /sigmay(idat,ilauf)**2
            enddo
            enddo
          END IF
          alpha(k,j)=alpha(j,k)
        enddo
      END DO
      curvature(1:nterms,1:nterms) =  alpha(1:nterms,1:nterms) ! save curvature
c
c     optional correction of diagonal elements
c
      if(ifdrc.eq.1) then
        do j=1,nterms
        write(2,'(a,i5,2x,2(1pg12.4),3x,2(1pg12.4))')
     &  ' 1. deriavtive: exact,approx,ratio for j=',
     &   j,betadr(j),beta(j),betadr(j)/beta(j)
        write(2,'(a,i5,2x,2(1pg12.4),3x,2(1pg12.4))')
     &  ' 2. deriavtive: exact,approx,ratio for j=',
     &   j,alphdr(j),alpha(j,j),alphdr(j)/alpha(j,j)
        alpha(j,j) = alphdr(j)
        beta(j)    = betadr(j)
        enddo
      endif
c
c  nicht-korrelierte fehler
c
      do j=1,nterms
        if(alpha(j,j).gt.0D0) then
         sigma0(j)=1./sqrt(alpha(j,j))
        else
         sigma0(j)=0D0
        endif
      enddo
c
c    chi**2 mit anfangswerten:
c
      chisq1=chiin
c
c        invert modified curvature matrix to find new parameters
c
      flamda = 0.1*flamda
      iwarnmax = 6
      if(nwdr.eq.0) iwarnmax=1
      IF(tonestep .AND. nwdr==1) iwarnmax=1 
      do 300 iwarn2=1,iwarnmax
      iwarns = iwarn2
      if(iwarn2.eq.1 .or. ireduc.eq.0) then
cold       flamda = .0001*10.**iwarn2
       flamda = 10.0*flamda
       if(iwarn2.gt.1)
     *   write(2,'(//a,i1,a/a,e10.2)')
     *             '  lamda erhoeht um 10**',iwarn,
     *     '   ohne verbesserung in chi**2','  lambda jetzt ',flamda
       do j=1,nterms
         do k=1,nterms
          zarray(j,k)=alpha(j,k)/sqrt(max(xlowlm,alpha(j,j)*alpha(k,k)))
         enddo
         zarray(j,j)=1.+flamda
       enddo
       call matin2(zarray,nterms,kparmx,nerr,det)
       if(iflast.eq.1) goto 399
       do j=1,nterms
         b(j)=a(j)
         do  k=1,nterms
            b(j)=b(j)+beta(k)*zarray(j,k)
     *         /sqrt(max(xlowlm,alpha(j,j)*alpha(k,k)))
       enddo
       enddo
       write(2,'(a/17f8.1)') 'B:',b(1:kparmx)
      else
        b(1:nterms)=a(1:nterms)+(b(1:nterms)-a(1:nterms))*.5
        write(02,'(a)') ' step shortened by factor 2. try again.'
      endif
      call flush(2)
c
c     berechne abschaetzung des neuen chi**2 aus schritt und
c     kruemmungsmatrix.
c
      sum1=0D0
      sum2=0D0
      do j=1,nterms
        sum1=(b(j)-a(j))*beta(j)+sum1
        do k=1,nterms
          sum2=(b(j)-a(j))*alpha(j,k)*(b(k)-a(k))+sum2
        enddo
      enddo
      chitst=chisqr-sum1-sum1+sum2
      write(02,'(4(a,g12.5))') ' estimate new chi**2:',chitst,
     *       '      for chisqr=',chisqr,'   linear contr.=',sum1,
     *       '   quadr. contr.=',sum2
!      WRITE(2,'(a,8(/5(1pg13.5)))') 'old parameters:',a(1:kparmx)
!      WRITE(2,'(a,8(/5(1pg13.5)))') 'new parameters:',b(1:kparmx)
c
c    berechne chi**2 mit den neuen parametern (feld b)
c
      ifsave=1
      if(ist1.ge.1) iwrte=1
      write(06,'(a)') ' evaluate new chi**2'
      chisu = ffverg(b)
      write(06,'(a,g12.4)') ' chi**2=',chisu
      if(ierror.ne.0) then
        write(02,'(3a)') ' the error ',errpri(ierror),
     *     ' occured at the new chi**2 evaluation! stop.'
        write(02,'(a,2i5)') ' nucleus was: nprot,nneut=',nprot,nneut
        chisqr=1./xlowlm
      else
       chisqr=chisu
       do kl=1,ndats
          chis00(kl)=sumchi(kl)
       enddo
      endif
       iwrte=0
       write(2,'(/a,g12.4,a,g12.4,a,i3/)')
     *          '  chi**2:   alt: ',chisq1,'     neu: ',chisqr,
     *         '     iwarn2 = ',iwarn2
       call flush(2)
       if(chisqr.le.chisq1.and.chisqr.le.chires) go to 399
       if(ireser.ne.0) goto 309
  300 continue
      call flush(2)
c
c   war chi*2 bei variation eines parameters kleiner geworden?
c   in diesem fall ist ireser=1 gesetzt
c
  309 continue
      if(ireser.ne.0) then
       write(2,'(/a,a/)')
     *  '  chi**2 advantagous at parameter variation ',
     *  '  --> quadratic step along this parameter .'
       do j=1,nterms
         parmom(j) = a(j)
       enddo
       parmom(ireser)=a(ireser)-dparmo(ireser)
       ifsave=0
       write(06,'(a,i4)') ' chi**2 for varied parameter nr.',-ireser
       chisu = ffverg(parmom)
       write(6,'(a,1g12.4)') ' chi**2=',chisu
       write(2,'(a,1g12.4)') ' chi**2=',chisu
       x1  = (chires-chisu)/(chiin+chiin-chires-chisu)
       x1  = sign(min(abs(x1),10*one),x1)
       write(2,'(a,4g12.4)')
     &   ' chires,chiin,chisu,x1=',chires,chiin,chisu,x1
       do 352 itry=1,4
       parmom(ireser) = a(ireser)+x1*dparmo(ireser)
       ifsave=1
       write(06,'(a)') ' try chi**2 in quadratic step.'
       chisu = ffverg(parmom)
       write(6,'(a,1g12.4)') ' chi**2=',chisu
       write(2,'(a,2g12.4)') ' x1,chi**2=',x1,chisu
       if(chisu.lt.chires) goto 353
       write(2,'(a)') ' no improvement in quadratic step, shorten '
       x1  = 0.7*x1
  352  continue
       write(2,'(a)') ' no improvement all quadratic steps '
       do j=1,nterms
         b(j)=ares(j)
         a(j)=ares(j)
       enddo
       ireser=0
       iwarn  = 1
       chisqr=chires
       goto 399
c
  353  continue
       write(2,'(a)') ' succesful quadratic step '
       do j=1,nterms
         b(j)=parmom(j)
         a(j)=parmom(j)
       enddo
       chisqr = chisu
      else
       iwarn = 99
      endif
c
c        evaluate parameters and uncertainties
c
  399 continue
      call flush(2)
c
c   hat sich chi**2 verbessert?
c
      verglc=chiin-chisqr
      if(verglc.lt.0.000) then
       chisqr=chiin
       do kl=1,ndats
         chis00(kl)=chisu0(kl)
       enddo
       do j=1,nterms
         b(j)=a(j)
       enddo
       itera = itera - 1
      else
       iwarn = 1
       if(iwarns.eq.1) flamda = max(0.001*one,flamda/10.0)
      endif
c
c    ausdruck der resultate des bevington steps
c
      if(iflast.ne.1) then
       write(2,'(/a)') '  neue parameterwerte:'
       if(iwarn.ge.99) write(2,'(a/)')
     *        '  = alte parameter (keine verbesserung in chi**2'
       write(02,'(5(2x,a,i2,a,g14.7))')
     *        (parbra,j,parket,b(j),j=1,nterms)
       chiwri=chisqr/rfree
       write(2,'(2(a,e12.4))') ' chi**2 pro freiheitsgrad = ',chiwri,
     *                ',      chi**2 gesamt            = ',chisqr
      endif
c
      if(iflast.ne.1 .and. iwarn.lt.99 .and. nwdr.ne.1) then
       do j=1,nterms
         a(j)=b(j)
       enddo
       do idat=1,ndats
       do ilauf=1,npts
         yfitmo(idat,ilauf)=yfit(idat,ilauf)
       enddo
       enddo
      endif
c
      do kl=1,ndats
        chisu0(kl)=chis00(kl)
      enddo
      call flush(2)
c
c     error return if no improvement has been achieved
c
      if(iwarn.ge.99) then
       write(02,'(a)') ' no improvment in chi**2. return early.'
       IF(.NOT.tonestep .AND. nwdr>1) return
      endif
c
c  fehlermatrix auf 'zarray'
c
      if(tdiagz) then
        IF(tpreigen) THEN
          WRITE(198) a(1:nterms)
          WRITE(198) (alpha(i,i),i=1,nterms)
          WRITE(199) a(1:nterms)
          WRITE(2,'(a)') 'Hessian matrix (unconditioned):'
          do n=1,nterms   
            WRITE(2,'(i3,30(1pg20.12))') 
     &                n,alpha(1:nterms,n)
          end do
          WRITE(199) alpha(1:nterms,1:nterms)
        END IF
        zarray(1:nterms,1:nterms)=alpha(1:nterms,1:nterms)
        call diag_real_symm(zarray,rootz,vectz,
     &                      nterms,kparmx,(nterms*nterms+nterms)/2)
        write(2,'(a)') 
     &  ' roots of Hessian matrix (unconditioned):'
        write(2,'(5(1pg20.12))') (rootz(i),i=1,nterms)
        IF(tpreigen) THEN
          write(199) rootz(1:nterms)
          WRITE(2,'(a)') 'Eigenvectors of Hessian matrix:'
          do n=1,nterms   
            WRITE(2,'(i3,30(1pg20.12))') 
     &                n,vectz(1:nterms,n)
          end do
          WRITE(199) vectz(1:nterms,1:nterms)
        END IF
      endif
      do j=1,nterms
      do k=1,nterms
       zarray(j,k)=alpha(j,k)/sqrt(alpha(j,j)*alpha(k,k))
      enddo
      enddo
        IF(tpreigen) THEN
          WRITE(2,'(a)') 'Hessian matrix (conditioned):'
          do n=1,nterms   
            WRITE(2,'(i3,30(1pg15.7))') 
     &                n,zarray(1:nterms,n)
          end do
          WRITE(198) zarray(1:nterms,1:nterms)
          WRITE(2,'(a)') 'uncorrelated variances for conditioning:'
          write(2,'(5(1pg13.5))') (1D0/SQRT(alpha(i,i)),i=1,nterms)
        END IF
      if(tinvert) then
        open(34,form='unformatted',file='alpha.dat')
        write(34) nterms
        write(34) ((alpha(j,k),j=1,nterms),k=1,nterms)
        close(34)
      endif      
      if(tdiagz) then
!        IF(tpreigen) THEN
!          WRITE(2,'(a)') 'Covariance matrix (conditioned):'
!          do n=1,nterms   
!            WRITE(2,'(i3,30(1pg15.7))') 
!     &                n,zarray(1:nterms,n)
!          end do
!        END IF
        call diag_real_symm(zarray,rootz,vectz,
     &                      nterms,kparmx,(nterms*nterms+nterms)/2)
        write(2,'(a)') 
     & ' roots of Hessian (condited):'
        write(2,'(5(1pg13.5))') (rootz(i),i=1,nterms)
        WRITE(198) rootz(1:nterms)
!
        IF(tsensitivity) THEN
          IF(tdiffobs) THEN
            CALL sensitivity(alpha,vectz,rootz,nterms,ndatsd)
          ELSE
            CALL sensitivity(alpha,vectz,rootz,nterms,ndats)
          END IF
        END IF

c             the matrix of natural parameters
        IF(tpreigen) WRITE(2,'(a)') 
     &   'eigenvectors of conditioned Hessian:'
        IF(tpreigen)  WRITE(198) vectz(1:nterms,1:nterms)
        do n=1,nterms     ! loop over eigenstates with length 1
          IF(tpreigen) WRITE(2,'(i3,1pg13.5,3x,30(1pg13.5))') 
     &                n,rootz(n),vectz(1:nterms,n)
          do i=1,nterms
            parmnat(i,n) = vectz(i,n)/sqrt(alpha(i,i)*rootz(n))
!     &                    *parmnatfac
          end do
        end do
        IF(tpreigen) THEN
          WRITE(2,'(a)') 'natural parameters:'
          do n=1,nterms   
            WRITE(2,'(i3,30(1pg15.7))') 
     &                n,parmnat(1:nterms,n)
          end do
        END IF
      endif
      call matin2(zarray,nterms,kparmx,nerr,det)
      IF(tpreigen) THEN
        WRITE(2,'(a,1pg13.5)') 
     &    'covariance matrix (conditioned): det=',det
        do n=1,nterms   
          WRITE(2,'(i3,30(1pg15.7))') 
     &                n,zarray(1:nterms,n)
        end do
        WRITE(2,'(a,1pg13.5)') 
     &    'covariance matrix (conditioned, normalized): det=',det
        do n=1,nterms   
          WRITE(2,'(i3,30(1pg15.7))') 
     &    n,(zarray(i,n)/SQRT(zarray(i,i)*zarray(n,n)),i=1,nterms)
        end do
      END IF
       if(tdiagz) then
         call diag_real_symm(zarray,rootz,vectz,
     &                       nterms,kparmx,(nterms*nterms+nterms)/2)
         write(2,'(a)') 
     &     ' roots of covariance matrix (conditioned):'
         write(2,'(5(1pg13.5))') (rootz(i),i=1,nterms)
         write(2,'(a)') 
     &     ' inverse roots of covariance matrix (conditioned):'
         write(2,'(5(1pg13.5))') (1D0/rootz(i),i=1,nterms)
         IF(tpreigen) THEN
           WRITE(2,'(a)') 
     &      'eigenvectors covariance matrix (conditioned):'
           do n=1,nterms  
             WRITE(2,'(i3,1pg13.5,3x,30(0pf9.5))') 
     &                n,rootz(n),vectz(1:nterms,n)
           enddo
         END IF
       endif
      do j=1,nterms
        sigmaa(j)=sqrt(zarray(j,j)/alpha(j,j))
      enddo
      do j=1,nterms
      do k=1,nterms
        zarray(j,k)=zarray(j,k)/sqrt(alpha(j,j)*alpha(k,k))
      enddo
      enddo
      if(tdiagz) then
        IF(tpreigen) THEN
          WRITE(2,'(a)') 'covariance matrix (unconditioned):'
          do n=1,nterms   
            WRITE(2,'(i3,30(1pg15.7))') 
     &                n,zarray(1:nterms,n)
         end do
!         CALL correl_data(nterms) 
        END IF
         call diag_real_symm(zarray,rootz,vectz,
     &                       nterms,kparmx,(nterms*nterms+nterms)/2)
         call diag_real_symm(zarray,rootz,vectz,
     &                       nterms,kparmx,(nterms*nterms+nterms)/2)
         write(2,'(a)')
     &     ' roots of covariance matrix (unconditioned):'
         write(2,'(5(1pg13.5))') (rootz(i),i=1,nterms)
         write(2,'(a)')
     &     ' inverse roots of coariance matrix (unconditioned):'
         write(2,'(5(1pg13.5))') (1D0/rootz(i),i=1,nterms)
      endif
!
!     test of 'parmnat'
      IF(tpreigen) THEN
        DO i=1,nterms; DO j=1,nterms
          vectz(i,j)=SUM(parmnat(i,:)*parmnat(j,:))
        END DO; END DO
        WRITE(2,'(a)') 'restored Hessian matrix (unconditioned):'
        do n=1,nterms   
          WRITE(2,'(i3,30(1pg15.7))') 
     &                n,vectz(1:nterms,n)
        end do
      END IF
c
      if(tinvert) then
        write(2,'(a)') ' test inversion by multiplication:'
        errmax = 0D0
        do j=1,nterms
          do k=1,nterms
            sum1 = 0D0
            do l=1,nterms
              sum1 = sum1 + zarray(j,l)*alpha(l,k)
            enddo
            rootz(k) = sum1
            if(k.eq.j) sum1 = sum1-1D0
            errmax = max(errmax,abs(sum1))
          enddo
          write(2,'(5(1pg14.6))') (rootz(k),k=1,nterms)
        enddo
        write(2,'(a,1pg13.5)') ' maximum deviation=',errmax
      endif
c
c     the "extrapolation errors" for fit observables
c
      ilauf=0
      do idat=1,ndats 
      do ipar=1,npts
        ilauf=ilauf+1
        varextrap(idat,ipar) = 0D0
        DO jpar1=1,nterms
        DO jpar2=1,nterms
          varextrap(idat,ipar) = varextrap(idat,ipar) +
     &       fderij(ilauf,jpar1)*zarray(jpar1,jpar2)*
     &       fderij(ilauf,jpar2)
        END DO
        END DO
      enddo
      enddo

c
c     print chi**2 matrix
c
      if(ist1.eq.1) then
       write(02,'(a)')   ' chi**2 matrix:'
       do k=1,nterms
        write(02,'(10f8.3)')
     &  (alpha(j,k)/sqrt(max(xlowlm,alpha(j,j)*alpha(k,k))),j=1,nterms)
       enddo
      endif
      call flush(2)
c
c     save parameters and corresponding chi**2 for curved step.
c
      if(itera.le.3) then
       do i=1,kparmx
         parsav(i,itera)=a(i)
       enddo
       write(2,'(a,10(1pg12.4))')
     & ' parsav=',(parsav(i,itera),i=1,kparmx)
       chisav(itera)=chisqr
      else
       do i=1,kparmx
        parsav(i,1) = parsav(i,2)
        parsav(i,2) = parsav(i,3)
        parsav(i,3)=b(i)
       enddo
       chisav(1) = chisav(2)
       chisav(2) = chisav(3)
       chisav(3) = chisqr
      endif
      chivrg = chisav(3)
c
c
c        now the curved step
c
      if(itera.ge.max0(3,ifcurv) .and. iflast.ne.1
     &   .and. chisav(1).ne.chisav(2)) then
       write(02,'(a,3g14.7)') ' chisav:',(chisav(j),j=1,3)
       write(02,'(a)') ' parsav:'
       do j=1,kparmx
         write(02,'(8x,3g14.7)') (parsav(j,i),i=1,3)
       enddo
       write(02,'(a)') ' running curved step '
c
cva1    x3=-1.
cva1    x1=+1.
c
cva2        x3=0.
cva2        do 720 i=1,15
cva2  720   if(sigmaa(i).gt.1.e-20)
cva2     *  x3=x3-(parsav(i,3)-parsav(i,2))**2/sigmaa(i)
cva2        x1=0.
cva2        do 730 i=1,15
cva2  730   if(sigmaa(i).gt.1.e-20)
cva2     *  x1=x1+(parsav(i,1)-parsav(i,2))**2/sigmaa(i)
c
cva3        x3=0.
cva3        do 740 i=1,15
cva3  740   if(sigma0(i).gt.1.e-20)
cva3     *  x3=x3-(parsav(i,3)-parsav(i,2))**2/sigma0(i)
cva3        x1=0.
cva3        do 750 i=1,15
cva3  750   if(sigma0(i).gt.1.e-20)
cva3     *  x1=x1+(parsav(i,1)-parsav(i,2))**2/sigma0(i)
c
c     curve in parameter space is fitted by a parabola, using
c     the curvature matrix of chi**2 as a metric.
c
        x3=0.
        do i=1,nterms
        do j=1,nterms
          x3=x3+(parsav(i,3)-parsav(i,2))*alpha(i,j)
     *               *(parsav(j,3)-parsav(j,2))
        enddo
        enddo
        x1=0.
        do i=1,nterms
        do j=1,nterms
          x1=x1+(parsav(i,1)-parsav(i,2))*alpha(i,j)
     *               *(parsav(j,1)-parsav(j,2))
        enddo
        enddo
        x1=sqrt(x1)
        x3=-sqrt(x3)
c
c     equilibrised curved step:
c      a third point on the tangential curve is generated with
c      same distance to point 3 as the distance point3 to point 2.
c      the minimum chi**2 is then extrapolated between this three
c      points.
c
       x0=x3+x3
       do i=1,kparmx
        parmom(i)= parsav(i,1)*(x0-x3)*x0/((x1-x3)*x1)
     *            +parsav(i,2)*(x0-x3)*(x0-x1)/(x3*x1)
     *            +parsav(i,3)*x0*(x0-x1)/(x3*(x3-x1))
       enddo
       ifsave = 0
       write(06,'(a)') ' chi**2 for curved step'
       chisus = ffverg(parmom)
       write(02,'(5(2x,a,i2,a,g14.7))')
     *       (parbra,j,parket,parmom(j),j=1,nterms)
       write(02,'(4(a,g12.5))') ' equilibrised step: x1=',x1,
     *         '   x3=',x3,'   x0=',x0,'   chisus=',chisus
      call flush(2)
       if(ierror.ne.0) then
        write(02,'(3a)') ' the error ',errpri(ierror),
     *     ' occured at equilibrised curved step! return.'
        write(02,'(a,2i5)') ' nucleus was: nprot,nneut=',nprot,nneut
        return
       endif
       if(chisus.lt.chivrg) then
         chivrg = chisus
         chisav(1) =chisav(2)
         chisav(2) = chisav(3)
         chisav(3) = chisus
         do j=1,kparmx
           parsav(j,1) = parsav(j,2)
           parsav(j,2) = parsav(j,3)
           parsav(j,3) = parmom(j)
         enddo
       endif
       x1     =  1.
       x3     = -1.
c
       d1chi  = -x3*chisav(1)/((x1-x3)*x1)
     *          -(x1+x3)*chisav(2)/(x3*x1)
     *          -x1*chisav(3)/((x3-x1)*x3)
       d2chi  = 2.*( chisav(3)/((x1-x3)*x1)
     *              +chisav(2)/(x3*x1)
     *              +chisav(1)/((x3-x1)*x3) )
       if(d2chi.gt.0.) then
cold        x0     = -d1chi/d2chi
        xlimit = 3.0
        x0     = -min(d1chi/d2chi,xlimit)
       else
         x0    = -2.
       endif
c      if(ist1.eq.1) then
        write(2,'(5(a,g12.5))') ' curved step: x1 = ',x1,
     *          '   x3 = ',x3,'   d1chi = ',d1chi,'   d2chi = ',
     *          d2chi,'   x0 = ',x0
c      endif
       xmin  = x3
       if(chisav(3).gt.chisav(2)) xmin  = 0D0
       do 600 itry=1,itrmax
       do i=1,kparmx
         parmom(i)= parsav(i,1)*(x0-x3)*x0/((x1-x3)*x1)
     *             +parsav(i,2)*(x0-x3)*(x0-x1)/(x3*x1)
     *             +parsav(i,3)*x0*(x0-x1)/(x3*(x3-x1))
       enddo
       ifsave = 1
       write(06,'(a)') ' trial chi**2 after curved step'
       chisu = ffverg(parmom)
       if(ierror.ne.0) then
        write(02,'(3a)') ' the error ',errpri(ierror),
     *     ' occured at extended curved step! proceed.'
        write(02,'(a,2i5)') ' nucleus was: nprot,nneut=',nprot,nneut
       endif
       if(chisu.lt.chivrg) goto 699
       write(02,'(5(2x,a,i2,a,g14.7))')
     *        (parbra,j,parket,parmom(j),j=1,nterms)
       write(2,'(a,i2,a,g12.5)') ' curved step: itry=',itry,
     *                           '   chi**2=',chisu
       x0=xmin+(x0-xmin)*.5
  600  continue
c
c     curved fit failed to produce lower chi**2 , reset saved
c     parameters  for next bevington step.
c
       write(2,'(a)') '0 no success in curved step'
       if(chisus.eq.chivrg) then
        write(2,'(a)') ' but success in trial step. keep this.'
        do i=1,kparmx
c        parsav(i,1)= parsav(i,2)
c        parsav(i,2)= parsav(i,3)
         a(i)       = parsav(i,3)
        enddo
c        chisav(1) =  chisav(2)
c        chisav(2) =  chisav(3)
        chisqr    =  chisav(3)
       endif
       return
c
c     curved fit succeeded
c
  699  continue
       write(2,'(a,i3)') '0 success in curved step at itry=',itry
       do i=1,kparmx
        parsav(i,3)=parmom(i)
        a(i)       =parmom(i)
       enddo
       chisav(3) = chisu
       chisqr    = chisu
c       do 650 idat=1,ndats
c       do 650 ilauf=1,npts
c       yfit(idat,ilauf)=yfitmo(idat,ilauf)
c  650  continue
       write(2,'(//a/)') '  neue parameterwerte'
       write(02,'(5(2x,a,i2,a,g14.7))')
     *        (parbra,j,parket,a(j),j=1,nterms)
       chiwri=chisu/rfree
       write(2,'(a,e12.4)') ' chi**2 pro freiheitsgrad = ',chiwri
       write(2,'(a,e12.4)') ' chi**2 gesamt            = ',chisu
      endif
c
      return
      end
c-----ffverg----------------------------------------------------------
c
      function ffverg(par)
c
#if(runomp)
      USE omp_lib
#endif
c
      include 'params.inc'
c
      logical tnucmat              ! flag for print n.m. properties in 31
      parameter (tnucmat=.true.)
c
      character*27 heatyp(12)
c
      include 'wavfun.inc'
      include 'grid.inc'
      include 'skypar.inc'
      include 'numpar.inc'
      include 'result.inc'
      include 'comfit.inc'
c      include 'constn.inc'
c
      dimension    par(kparmx),npoipr(5),npoine(6)
      REAL(8),ALLOCATABLE :: spetmp(:),rhotmp(:)
      REAL(8),ALLOCATABLE :: workpol(:,:)
c
      data zero,x10lw/0D0,1.e-10/
cold      data zero,x10up,x10lw/0D0,1.e+10,1.e-10/
      data heatyp/'  binding-energy',
     *            '  diffr.-radius',
     *            '  surface thickness',
     *            '  rms radius',
     *            '  l*s splitting proton 16O',
     *            '  l*s splitting neutron 16O',
     *            '  (rho(1.8)+rho(0.))/2.',
     *            '  neutron rms radius',
     *            '  gap-proton',
     *            '  gap-neutron',
     *            '  rho(7.2) or 2. l*s neutr.',
     *            '  m*/m symm. nuclear matter'/
c
ciso      data isav82/0/   ! storage of pointer for isotope shifts in z=82
ciso      data isav38/0/   ! storage of pointer for isotope shifts in z=38
ciso      data isav20/0/   ! storage of pointer for isotope shifts in z=20
c
c
c      pointer fields to the last filled shells in 208pb.
c
      data npoipr/12,13,14,15,16/
      data npoine/17,18,19,20,21,22/
      LOGICAL, PARAMETER :: tswitchca=.FALSE.
c      logical tsimple                      ! switch to simple test model
c      data tsimple/.true./
      LOGICAL,PARAMETER :: ttest=.FALSE.
      LOGICAL,PARAMETER :: tswitchpolariz=.TRUE.
#if(runomp)
      LOGICAL,PARAMETER :: tsubcall=.TRUE.
#else
      LOGICAL,PARAMETER :: tsubcall=.FALSE.
#endif
      CHARACTER(1) :: nodename(0:7)=(/'0','1','2','3','4','5','6','7'/)
      LOGICAL :: tfirst=.TRUE.
c
c---------------------------------------------------------------------
c
#if(runomp)
      IF(tfirst) THEN
        IF(.NOT.tsubcall) STOP 'OPENMP requires TSUBCALL=.T.'
        numthr = 8 !OMP_GET_MAX_THREADS()
        call OMP_SET_NUM_THREADS(numthr)
        maxthr  = numthr-1
        DO myn=0,maxthr
          CALL system('mkdir '//nodename(myn))
        END DO
        IF(ttest) WRITE(2,*) 'number of threads=',numthr
        tfirst=.FALSE.
      END IF
#endif
c
      IF(nwdr>0) ifplot=0
c
c      write(6,*) ' FFVERG: IWR=',iwr
cold      do 90 i=1,nparin
cold      par1(i)=par(i)
cold   90 continue
c
      do i=1,12
      do j=1,kdat
        yfit(i,j)=zero
      enddo
      enddo
c
      if(ist1.ge.1) write(2,105) t0,t1,t2,t3,t4,x0,x1,x2,x3,power,
     +                         tyuk,xyuk,ayuk
  105 format(/5f8.2,5f8.5,f8.2,2f8.5)
c
      do k=1,ndats
       sumchi(k)=zero
      enddo
c
      if(tnucmat) 
     &  write(31,'(a)') '     ','      ','#'//commen(1:20)
#if(runomp)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(nra,nnumb,anumb,nprot,zrprot,
!$OMP&   nneut,nmaxp,nmaxn,ngrid,etotal,rdmsc,surfc,rmsc,rmsn,deltp,
!$OMP&   deltn,rhotmp,spetmp,yy7,yy8,soprot,soneut,myn,handshake,iu,
!$OMP&   nprotc,nneutc,nmaxpc,nmaxnc)
#endif
      ALLOCATE(spetmp(ksl2),rhotmp(kg1))
!    WRITE(7,*) ' PSISTATE allocated, thread=',OMP_GET_THREAD_NUM()
#if(runomp)
!$OMP DO SCHEDULE(DYNAMIC)
#endif
      do 1000 i=1,npts
#if(runomp)
      myn=OMP_GET_THREAD_NUM()
#else
      myn=0
#endif
      IF(ttest) WRITE(2,*) ' MYN=',myn,i
      call attach(par,jp)
!      WRITE(6,'(6(1pg13.5))') b0,b1,b2,b3,b0p,b1p,b2p,b3p,power
c
c     mean-field calculation: first set the the parameters for
c        the actual nucleus which are then transferred to 'rmf'
c        via common's.
c
      nra    = iay(i)
      nnumb  = nra
cold      rnra   = nra
      anumb  = nnumb
      nprot  = iprot(i)
      zrprot = zprot(i)
c                            switch to finite nuclei
      if(nra.lt.999) then
       nneut  = nra-nprot
       nmaxp  = 0
       nmaxn  = 0
       ngrid  = 0
       IF(ngrid_explicit>0) THEN
         ngrid=ngrid_explicit
       END IF
!       WRITE(02,'(a,2i5)') 'FF: ngrid_explicit=',ngrid_explicit,ngrid
       IF(tswitchca .AND. nprot==20 .AND. nneut==20) THEN
         nmaxp = 6
         nmaxn = 6
       ELSE IF(tswitchca .AND. nprot==20 .AND. nneut==28) THEN
         nmaxp = 6
         nmaxn = 7
       ELSE IF(tswitchca .AND. nprot==20) THEN
         nmaxp = 6
       END IF
       if(ipair.lt.5 .or. kparmx.le.15) then
         pairfp = 0D0
         pairfn = 0D0
       endif
c
       IF(tsimple) THEN
c        WRITE(6,'(a,2(1pg13.5))') ' before ESIMPLE:',par(jp(1)),t0
         call esimple
       ELSE
!        WRITE(6,'(6(1pg13.5))') b0,b1,b2,b3,b0p,b1p,b2p,b3p,power
!        call hafosu
         IF(tsubcall) THEN         
          iu=930+myn
          open(unit=iu,file=nodename(myn)//'/for005')
          CALL handshake05(iu,nprot,nneut,ngrid,rstep)
          close(unit=iu)
          IF(ttest) WRITE(2,'(a,2i4)') 
     &          'file for005: node,nucleus='//nodename(myn),
     &             nprot,nneut
          CALL SYSTEM('cd '//nodename(myn)//
     &               '; ../haforpa.exe < for005 > /dev/null;
     &                 cd .. ')
          IF(ttest) WRITE(2,'(a)') nodename(myn)//': haforpa done'
          OPEN(unit=iu,file=nodename(myn)//'/result.000',
     &         form='unformatted')
          READ(iu) nprotc,nneutc,nmaxpc,nmaxnc,etotal,rdmsc,surfc,
     &              rmsc,rmsn,deltp,deltn,rhotmp,spetmp
!          WRITE(6,'(a,200(1pg13.5))') 'spetmp:',spetmp(1:100)
          CALL flush(6)
          CLOSE(iu)
          IF(ttest) WRITE(2,'(a,4i4)') nodename(myn)//':result read',
     &             nprot,nneut,nprotc,nneutc
          IF(nprotc.NE.nprot .OR. nneutc.NE.nneut) 
     &      STOP 'reading RESULT.000: nucleus not matching'
         ELSE
          ipair=ipairsav
          call blocking
          spetmp=spenrg
          rhotmp=rhochr
          nmaxpc=nmaxp
         END IF
       ENDIF
       IF(ierror.EQ.2) THEN
         write(02,'(a,2i5)') ' too large nr. iterations: nprot,nneut=',
     &    nprot,nneut
         WRITE(2,'(a,5(1pg13.5))') 'etotal,drmst,detot=',
     &            etotal,drmst,detot
!         do j=1,nparfit
!           write(2,'(g15.8,a,e12.5,a,i2,a)')  
!     &        par(j),',',dparmo(j),',',iparva(j),','
!         enddo
         write(2,'(50(1g15.8))') (par(j),j=1,nparfit)
         CALL flush(2)
       END IF
       IF(.NOT.tsimple) THEN
       if(nprot.eq.108 .and. nneut.eq.156 .or.
     &    nprot.eq.94 .and. nneut.eq.152 .or.
     &    nprot.eq.66 .and. nneut.eq.98) then    ! special deformed
         etotal = energy_deform(nprot,nneut)
       endif
       ENDIF
c
c      now extract the results which are transferred back to
c      this routine via common's
c
 1001  continue
c
       yfit(lifit(1),i)=etotal
!       WRITE(2,*) 'ENERGY:',i,lifit(1),nra,etotal
c
       yfit(lifit(2),i)=rdmsc
       if(ifitr.gt.0) yfit(lifit(2),i)=(rdmsc+0.011-0.323*surfc)/0.975
c
       yfit(lifit(3),i)=surfc
c
       yfit(lifit(4),i)=rmsc
c
       yy7=(rhotmp(7)+rhotmp(1))/2.*100.
       yy8=(rhotmp(7)-rhotmp(1))*100.
       yfit(lifit(7),i)=yy7
cold       yfit(lifit(8),i)=yy8
       yfit(lifit(8),i)= apv_weak !rmsn   !  -rmsp
c
cold       yfit(lifit(9),i)=rhotmp(13)*100.
       yfit(lifit(9),i)=deltp
c   kern ist nicht 16-o:
cold       yfit(lifit(10),i)=rhotmp(19)*100.
       yfit(lifit(10),i)=deltn
c
       yfit(lifit(11),i)=rhotmp(25)*100.
       if(nra.eq.16.and.nprot.eq.8) then
c       beitrag zu chi**2 von schalenenergie-differenzen
         soprot=spetmp(3)-spetmp(2)
         soneut=spetmp(nmaxpc+3)-spetmp(nmaxpc+2)
         yfit(lifit(5),i)=soprot
         yfit(lifit(6),i)=soneut
       elseif(nra.eq.40.and.nprot.eq.20) then
         yfit(lifit(5),i)=spetmp(5)-spetmp(4)
         yfit(lifit(6),i)=spetmp(nmaxpc+5)-spetmp(nmaxpc+4)
       elseif(nra.eq.48.and.nprot.eq.20) then
         yfit(lifit(5),i)=spetmp(5)-spetmp(4)
       elseif(nra.eq.90.and.nprot.eq.40) then
         yfit(lifit(5),i)=spetmp(12)-spetmp(11)
         yfit(lifit(6),i)=spetmp(nmaxpc+12)-spetmp(nmaxpc+11)
       elseif(nra.eq.132.and.nprot.eq.50) then
         yfit(lifit(5),i)=spetmp(10)-spetmp(9)                 ! 2p prot
         yfit(lifit(6),i)=spetmp(nmaxpc+14)-spetmp(nmaxpc+13)    ! 2d neut
       elseif(nra.eq.208.and.nprot.eq.82) then
         yfit(lifit(5),i)=spetmp(14)-spetmp(13)                ! 2d prot
         yfit(lifit(11),i)=spetmp(nmaxpc+19)-spetmp(nmaxpc+18)   ! 1fneut  
         yfit(lifit(6),i)=spetmp(nmaxpc+21)-spetmp(nmaxpc+20)    ! 3p neut
       elseif(nra.eq.207 .and. nprot.eq.82) then
c       elseif(nra.eq.208 .and. nprot.eq.82) then
c        beitrag zu chi**2 von level-width in 208pb
        xavera = zero
        eavera = zero
        e2aver = zero
        do n=1,5
          npl    = npoipr(n)
          xmult  = 2.0*float((jsp(npl)+1)/2)
          xavera = xmult + xavera
          eavera = xmult*spetmp(npl) + eavera
          e2aver = xmult*spetmp(npl)**2 +e2aver
        enddo
        widpr  = sqrt(e2aver/xavera-(eavera/xavera)**2)
        yfit(lifit(5),i) = widpr
        xavera = zero
        eavera = zero
        e2aver = zero
        do n=1,6
          npl    = npoine(n)+16
          xmult  = 2.0*float((jsp(npl)+1)/2)
          xavera = xmult + xavera
          eavera = xmult*spetmp(npl) + eavera
          e2aver = xmult*spetmp(npl)**2 +e2aver
        enddo
        widne  = sqrt(e2aver/xavera-(eavera/xavera)**2)
        yfit(lifit(6),i) = widne
        fluc   = yy8
       endif
c
c     for nra.ge.999 compute nuclear matter.
c      first option is nprot.ge.999 
c       --> enter total density via data column 2 and
c                 proton fraction via data column 3
c                 resulting E/A via data column 4
c                 note: columns 2 & 3 are (mis)used as dummy data
      elseif(nprot.ge.999) then
       rhot=yexp(lifit(2),i)
       procon=yexp(lifit(3),i)
       call attach(par,jp)
       call initfc
       call enernm(rhot,procon,zero,energm,efmasp,efmasn)
       yfit(lifit(4),i)=energm
       yfit(lifit(3),i)=yexp(lifit(3),i)  
       yfit(lifit(2),i)=yexp(lifit(2),i)
c       write(6,'(a,5(1pg12.4))')
c     &  'rhot,procon,energm,efmasp,efmasn=',
c     &   rhot,procon,energm,efmasp,efmasn
c      second option is nprot.gt.0 for which the equation of state
c      is computed and compared with data.
c
      elseif(nprot.gt.0) then
       fkf=nprot/100.
       rhonm=0.033773727*fkf**3   !0.0675475*fkf*fkf*fkf
       call attach(par,jp)
       call initfc
       call enernm(rhonm,0.5d0,zero,energm,efmasp,efmasn)
       yfit(lifit(4),i)=energm
c
c      third option is nprot=0 for which the ground state of nuclear
c      matter is computed and compared with pseudodata.
c
      elseif(ABS(zrprot)<1D-11) then
       call attach(par,jp)
       call initfc
       call nucmat(rhonm,eovera,compr,effmas,f0,f0p,g0,g0p,
     &             asymme,dasymm,d2asym,skappa,0,-2)
       if(ioptsk.eq.4) call xkappa(skappa)
       yfit(lifit(1),i)=eovera
       yfit(lifit(2),i)=rhonm
       yfit(lifit(3),i)=compr
       yfit(lifit(4),i)=asymme
! misuse 'g=' to communicate polarizability
       IF(tswitchpolariz) THEN
         ALLOCATE(workpol(12,kdat))
         IF(ABS(yexp(lifit(5),i))> 0D0) THEN
           call extobs_all(par,workpol,1,5)   ! polarizability in 208Pb
!         call extobs_all(par,workpol,1,7)   ! presently low E_2 in 48Ca
           WRITE(2,*) 'POLARIZ 208Pb:',i,lifit(5),nra,yfit(lifit(5),i),
     &                 yfit(lifit(1),1)
           yfit(lifit(5),i)=resenergy()
         ELSE
           yfit(lifit(5),i)=0D0
         END IF
         IF(ABS(yexp(lifit(6),i))> 0D0) THEN
           call extobs_all(par,workpol,1,7)   ! polarizability in 208Pb
           yfit(lifit(6),i)=resenergy()
           WRITE(2,*) 'POLARIZ  48Pb:',i,lifit(5),nra,yfit(lifit(5),i),
     &                yfit(lifit(1),1)
         ELSE
           yfit(lifit(6),i)=0D0
         END IF
         DEALLOCATE(workpol)
       ELSE         
         yfit(lifit(5),i)=g0
         yfit(lifit(6),i)=g0p
       END IF
       yfit(lifit(8),i)=dasymm
c       write(02,'(a,2i5,2(1pg13.5))')
c     &  ' dasymm:',i,lifit(8),dasymm,yfit(lifit(8),i)
       yfit(lifit(10),i)=str_squark               ! strength of strange quark
       yfit(lifit(11),i)=skappa
       yfit(lifit(12),i)=effmas
c
c      fourth option is nprot.lt.0 for which the equation of state
c      is computed and compared with data.
c
      elseif(zrprot.lt.0D0) then
       fkf=-zrprot/100D0
       rhonm=0.0337737*fkf*fkf*fkf
       call initfc
       call enernm(rhonm,zero,zero,energm,efmasp,efmasn)
!       WRITE(2,'(a,3(1pg13.5))') 'EoS: ',fkf,rhonm,energm
       yfit(lifit(4),i)=energm
      endif
c
      do 50 l=1,12
      if(abs(yexp(lifit(l),i)).lt.x10lw) yfit(lifit(l),i)=zero
   50 continue
c
      if(tnucmat) 
     &  write(31,'(2i4,f5.0,f5.1,f6.0,f5.2,f6.2,2(1pg12.4))')
     &  iprot(i),iay(i)-iprot(i),
     &  compr,asymme,dasymm,skappa,effmas,
     &  rmsn-rmsp,rdmsn-rdmsp
c     &  yfit(lifit(3),i),yfit(lifit(4),i),
c     &  yfit(lifit(8),i),yfit(lifit(11),i),yfit(lifit(12),i),
c
 1000 continue
#if(runomp)
!$OMP END DO
#endif
      DEALLOCATE(spetmp,rhotmp)
#if(runomp)
!$OMP END PARALLEL
#endif
c
c     finally sum up the chi**2. optionally write table of chi**2.
c
      if(iwrte.eq.1.and.ist1.gt.1) write(2,588)
  588 format(///,3x,'  a  z   n       f(exp)',8x,'d-f(exp)',8x,
     +'f(ber)    (yfit-yexp)/df      s u m m e  c h i  ** 2',
     +'     f(ber)0   f(ber)0 -',
     +/,65x,'chi          dieser satz     gesamt',18x,'f(ber)',/)
c
      chisu=0.
!      call attach(par,jp)
!      chisu = (str_squark/1.64D0)**2               ! offset from strange-quark
!      write(2,'(a,2(1pg13.5))') 
!     &    ' str_squark,chi_offset=',str_squark,chisu
      do 600 k=1,ndats
      if(iwrte.eq.1.and.ist1.gt.1) write(2,'(a,i2/a)')
     *       '   daten satz ',k,heatyp(ifit1(k))
c
      do 600 i=1,npts
      chi(k,i)=0.
      if(sigmay(k,i).gt.x10lw) then
       chi(k,i)=(yfit(k,i)-yexp(k,i))/sigmay(k,i)
c       write(02,'(a,2i5,3(1pg13.5))')
c     &  ' chi:',k,i,yfit(k,i),yexp(k,i),sigmay(k,i),chi(k,i)
cisoc
cisoc     extra check for two special isotope shifts
cisoc
cisoc       write(2,'(a,3i5,a,2i5,g12.4)')
cisoc     &  ' k,i,lifit=',k,i,lifit(k),
cisoc     &  ' : a,z,sigmay=',iay(i),iprot(i),sigmay(k,i)
ciso      if(k.eq.lifit(4).and.iprot(i).eq.82) then
ciso         if(iay(i).eq.208) then
cisoc            write(2,'(a)') ' pb208 seen'
ciso            isav82 = i
ciso         elseif(iay(i).eq.214 .and. isav82.gt.0) then
cisoc            write(2,'(a,g11.3)') ' pb214 seen. chi=',chi(k,i)
cisoc            write(2,'(a,g11.3)') ' old yfit=',yfit(k,i)
ciso            shift  = yfit(k,i)**2-yfit(k,isav82)**2
ciso            shifte = yexp(k,i)**2-yexp(k,isav82)**2
ciso            chi(k,i) = (shift-shifte)/(0.01*sigmay(k,i))
ciso            if(iwr.lt.-1) then
ciso             yfit(k,i) = sigmay(k,i)*chi(k,i)+yexp(k,i)   ! pseudo 'yfit' for xd1
cisoc             write(2,'(a,2(1pg12.4))')
cisoc     &        ' isotope shifts pb214-pb208, exact, fit=',shifte,shift
ciso            else
ciso             write(2,'(a,2(1pg12.4))')
ciso     &        ' isotope shifts pb214-pb208, exact, fit=',shifte,shift
ciso            endif
cisoc            write(2,'(a,g11.3)') ' new chi=',chi(k,i)
cisoc            write(2,'(a,g11.3)') ' new yfit=',yfit(k,i)
ciso            isav82 = 0
ciso         endif
ciso      endif
ciso      if(k.eq.lifit(4).and.iprot(i).eq.38) then
ciso         if(iay(i).eq.88) then
cisoc            write(2,'(a)') ' sr88 seen'
ciso            isav38 = i
ciso         elseif(iay(i).eq.84 .and. isav38.gt.0) then
cisoc            write(2,'(a,g11.3)') ' sr80 seen. chi=',chi(k,i)
cisoc            write(2,'(a,g11.3)') ' old yfit=',yfit(k,i)
ciso            shift  = yfit(k,i)**2-yfit(k,isav38)**2
ciso            shifte = yexp(k,i)**2-yexp(k,isav38)**2
ciso            chi(k,i) = (shift-shifte)/(0.01*sigmay(k,i))
ciso            if(iwr.lt.-1) then
ciso             yfit(k,i) = sigmay(k,i)*chi(k,i)+yexp(k,i)   ! pseudo 'yfit' for xd1
cisoc             write(2,'(a,2(1pg12.4))')
cisoc     &        ' isotope shifts sr84-sr88, exact, fit=',shifte,shift
ciso            else
ciso             write(2,'(a,2(1pg12.4))')
ciso     &        ' isotope shifts sr84-sr88, exact, fit=',shifte,shift
ciso            endif
cisoc            write(2,'(a,g11.3)') ' new chi=',chi(k,i)
cisoc            write(2,'(a,g11.3)') ' new yfit=',yfit(k,i)
ciso            isav38 = 0
ciso         endif
ciso      endif
ciso      if(k.eq.lifit(4).and.iprot(i).eq.20) then
ciso         if(iay(i).eq.40) then
cisoc            write(2,'(a)') ' ca40 seen'
ciso            isav20 = i
ciso         elseif(iay(i).eq.48 .and. isav20.gt.0) then
cisoc            write(2,'(a,g11.3)') ' ca48 seen. chi=',chi(k,i)
cisoc            write(2,'(a,g11.3)') ' old yfit=',yfit(k,i)
ciso            shift  = yfit(k,i)**2-yfit(k,isav20)**2
ciso            shifte = yexp(k,i)**2-yexp(k,isav20)**2
ciso            chi(k,i) = (shift-shifte)/(0.01*sigmay(k,i))
ciso            if(iwr.lt.-1) then
ciso             yfit(k,i) = sigmay(k,i)*chi(k,i)+yexp(k,i)   ! pseudo 'yfit' for xd1
cisoc             write(2,'(a,2(1pg12.4))')
cisoc     &        ' isotope shifts ca48-ca40, exact, fit=',shifte,shift
ciso            else
ciso             write(2,'(a,2(1pg12.4))')
ciso     &        ' isotope shifts ca48-ca40, exact, fit=',shifte,shift
ciso            endif
cisoc            write(2,'(a,g11.3)') ' new chi=',chi(k,i)
cisoc            write(2,'(a,g11.3)') ' new yfit=',yfit(k,i)
ciso            isav20 = 0
ciso         endif
ciso      endif
cisoc
cisoc     end of extra check
cisoc
      endif
      sumchi(k)=sumchi(k)+chi(k,i)**2
      if(k.le.ndats) chisu=chisu+chi(k,i)**2
      nneut=iay(i)-iprot(i)
      if(iwrte.eq.1.and.ist1.gt.1) then
       relf=(yexp(k,i)-yfit(k,i))/yexp(k,i)
       write(02,'(i3,i4,i3,i4,6e15.6,1x,e15.6,f10.6)')
     *              i,iay(i),iprot(i),nneut,yexp(k,i),sigmay(k,i),
     *              yfit(k,i),chi(k,i),sumchi(k),chisu,relf
      endif
  600 continue
c
c     evalute contribution from isotopic differences
c
      if(maxdif.gt.0) then
      do 630 idif=1,maxdif
      k     = 13                         ! storage number for isotopoc shifts
      inuc1 = nucdif(1,idif)
      inuc2 = nucdif(2,idif)
      iobs  = lifit(nucdif(3,idif))
!      if(nucdif(3,idif).ne.4) then
      if(nucdif(3,idif).ne.4 .and. nucdif(3,idif).ne.1) then
c        relf1 =(yexp(iobs,inuc1)-yfit(iobs,inuc1))/sigmay(iobs,inuc1)
c        relf2 =(yexp(iobs,inuc2)-yfit(iobs,inuc2))/sigmay(iobs,inuc2)
c        chi(k,idif) = (relf1-relf2)/errdif(idif)
        errf1 =(yexp(iobs,inuc1)-yfit(iobs,inuc1))
        errf2 =(yexp(iobs,inuc2)-yfit(iobs,inuc2))
        chi(k,idif) = (errf1-errf2)/errdif(idif)
      else
        IF(nucdif(3,idif)==4) THEN
          shift  = yfit(iobs,inuc1)**2-yfit(iobs,inuc2)**2
          shifte = yexp(iobs,inuc1)**2-yexp(iobs,inuc2)**2
        ELSE
          shift  = yfit(iobs,inuc1)-yfit(iobs,inuc2)
          shifte = yexp(iobs,inuc1)-yexp(iobs,inuc2)
        END IF
        chi(k,idif) = (shift-shifte)/(errdif(idif)*sigmay(iobs,inuc1))
        if(iwr.ge.-1)
!        if(iwr.ge.0)
     &    write(2,'(4(a,i4),a,5(1pg12.4))')
     &        ' isotope shifts (',iay(inuc1),',',iprot(inuc1),
     &        ')-(',iay(inuc2),',',iprot(inuc2),
     &        '): exact, fit=',shifte,shift,(shift-shifte),
     &        errdif(idif)*sigmay(iobs,inuc1),chi(k,idif)
        endif
      chisu = chi(k,idif)**2 + chisu
c      write(6,'(a,4i5,g12.4)')
c     &  ' idif,inuc1,inuc2,iobs,chi(k,idif)=',
c     &    idif,inuc1,inuc2,iobs,chi(k,idif)
  630 continue
      endif
c  610 continue
c
c     evaluate three-point gaps
c
      IF(max3gap>0) THEN
        DO i3=1,max3gap
          inuc1 = n3gap(1,i3)
          inuc2 = n3gap(2,i3)
          inuc3 = n3gap(3,i3)
          iobs  = n3gap(5,i3) 
          gap3theo = 0.5D0*(yfit(iobs,inuc1)-2D0*yfit(iobs,inuc2)
     &                      +yfit(iobs,inuc3))
          gap3exp = 0.5D0*(yexp(iobs,inuc1)-2D0*yexp(iobs,inuc2)
     &                      +yexp(iobs,inuc3))
          chi(14,i3) = (gap3theo-gap3exp)/err3gap(i3)
          chisu = chi(14,i3)**2 + chisu
          IF(iwr.ge.-1) THEN
!          if(iwr.ge.0)
            WRITE(2,'(2(a,i4),a,2(1pg12.4))')
     &        ' 3gap(',iay(inuc2),',',iprot(inuc2),
     &        '): exact, fit=',gap3exp,gap3theo
!            WRITE(2,'(a,5i4)') 'details:',i3,n3gap(1:4,i3)
          END IF
        END DO
      END IF
c
      ffverg = chisu
c
      return
      end
c-----extnuc--------------------------------------------------------
c
      subroutine extnuc(par,workdv,nparmx)
c
c     extrapolate to nuclei including extrapolation error using
c     matrix 'zarray' which is transferred via common /fitpar/ .
c
      include 'params.inc'
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'comfit.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'constn.inc'
c
c     'workdv' is field for derivatives with respect to parameter-variation.
c
      character*1 procen
      dimension   par(kparmx),workdv(kparmx,6)
      dimension   nuclei(2,3)
c
      data procen/'%'/
      data nuclei/ 28,50, 50,50, 50,82/
c      logical tsimple
c      data tsimple/.true./
cold      data x10up/1.e+10/
c
c-----------------------------------------------------------------------
c
      write(2,'(/a/3a)')
     & ' other nuclei including extrapolation errors:',
     & '   z   n    ',
     & '   e         rdms        surfc        rmsc',
     & '      rdms_halo    surf_halo'
c
c     loop over the test nuclei
c
      serre=0D0
      serrrd=0D0
      serrs=0D0
      serrrms=0D0
      do 1010 inuc=1,npts+3
      if(inuc.le.npts) then
       if(iay(inuc).gt.220) goto 1010
      endif
c
c     loop over variation of parameters
c
      do 10 i=0,nparmx
      if(i.gt.0) par(i) = par(i) + dparmo(i)
      call attach(par,jp)
c
      if(inuc.le.npts) then
       nprot  = iprot(inuc)
       nneut  = iay(inuc)-nprot
      else
       nprot  = nuclei(1,inuc-npts)
       nneut  = nuclei(2,inuc-npts)
      endif
      nnumb  = nprot+nneut
      anumb  = nnumb
cold      gapp   = -x10up
cold      gapn   = -x10up
      nmaxp  = 0
      nmaxn  = 0
cold      ngrid  = 0
c
      IF(tsimple) THEN
        call esimple
      ELSE
!        call hafosu
        ipair=ipairsav
        call blocking
      ENDIF
      rhalo  = rdmsn-rdmsp
      shalo  = surfn-surfp
c
      if(i.eq.0) then
       write(2,'(2(1x,i3),2x,6(g12.4))')
     &    nprot,nneut,etotal,rdmsc,surfc,rmsc,rhalo,shalo
       esav   = etotal
       rdsav  = rdmsc
       ssav   = surfc
       rmssav = rmsc
       rhalos = rhalo
       shalos = shalo
      else
       workdv(i,1) = (etotal-esav)/dparmo(i)
       workdv(i,2) = (rdmsc-rdsav)/dparmo(i)
       workdv(i,3) = (surfc-ssav)/dparmo(i)
       workdv(i,4) = (rmsc-rmssav)/dparmo(i)
       workdv(i,5) = (rhalo-rhalos)/dparmo(i)
       workdv(i,6) = (shalo-shalos)/dparmo(i)
      endif
c
      if(i.gt.0) par(i) = par(i) - dparmo(i)
   10 continue
c
c     accumulate extrapolation errors by folding with 'zarray'.
c
      erre   = 0D0
      errrd  = 0D0
      errs   = 0D0
      errrms = 0D0
      errrha = 0D0
      errsha = 0D0
c
      if(nparmx.gt.0) then
       do 40 i=1,nparmx
       do 40 j=1,nparmx
       erre   = workdv(i,1)*workdv(j,1)*zarray(i,j) + erre
       errrd  = workdv(i,2)*workdv(j,2)*zarray(i,j) + errrd
       errs   = workdv(i,3)*workdv(j,3)*zarray(i,j) + errs
       errrms = workdv(i,4)*workdv(j,4)*zarray(i,j) + errrms
       errrha = workdv(i,5)*workdv(j,5)*zarray(i,j) + errrha
       errsha = workdv(i,6)*workdv(j,6)*zarray(i,j) + errsha
   40  continue
       erre   = sqrt(erre)
       errrd  = sqrt(errrd)
       errs   = sqrt(errs)
       errrms = sqrt(errrms)
       errrha = sqrt(errrha)
       errsha = sqrt(errsha)
      endif
c
c     final print
c
      if(nparmx.gt.0) then
       write(2,'(10x,6(g12.4))')
     &    erre,errrd,errs,errrms,errrha,errsha
       write(2,'(10x,6(f9.3,a1,2x))')
     &    100.0*erre/esav,procen,100.0*errrd/rdsav,procen,
     &    100.0*errs/ssav,procen,100.0*errrms/rmssav,procen
      IF(inuc.LE.npts) THEN
        serre=serre+erre**2
        serrrd=serrrd+errrd**2
        serrs=serrs+errs**2
        serrrms=serrrms+errrms**2
      END IF
      endif
c
 1010 continue
       write(2,'(a,6(g12.4))')  'average:  ',
     &    SQRT(serre/npts),SQRT(serrrd/npts),SQRT(serrs/npts),
     &    SQRT(serrrms/npts)
c
      return
      end
c-----nucvar----------------------------------------------------------
c
      subroutine nucvar(par,nparmx,zuord)
c
c     computes extrapolated errors on nuclear matter properties
c
c      par    = array of model parameters
c      dparmo = array of variations of model parameters
c      nparmx = number of varied parameters
c               0 --> no evaluation of extrapolated errors
c
      include 'params.inc'
      include 'constn.inc'
c      parameter (zero=0D0)
c
      include 'comfit.inc'
      include 'skypar.inc'
c
      character*6 zuord(kparmx)
      DIMENSION DEOVRA(kparmx),DRHONM(kparmx),DCOMPR(kparmx),
     *          DSYMEN(kparmx),dskapp(kparmx),DDSYME(kparmx)
      REAL(8) :: DD2SYME(kparmx)
      dimension DEFMAS(kparmx),DDENEUT(kparmx),DENEUT0(kparmx)
      dimension par(kparmx),dtto(kparmx),dtte(kparmx)
      dimension db4(kparmx),db4p(kparmx),dvpp(kparmx),dvpn(kparmx)
      dimension drho0pr(kparmx)
c
c----------------------------------------------------------------------
c
      DO 10 I=0,NPARMX
      IF(I.GT.0) PAR(I) = PAR(I) + DPARMO(I)
      CALL ATTACH(PAR,JP)
      call initfc
      call nucmat(rhonm0,eovra0,compr0,efmas0,f0,f0p,g0,g0p,
     &            symen0,dasym0,d2asy0,skapp0,0,-2)
      if(ioptsk.eq.4)  call xkappa(skapp0)
c
c     neutron matter at rho=0.01 fm**-3 and slope thereof
c
      call enernm(0.10d0,zero,zero,eneut0,efmasp,efmasn)
      call enernm(0.1d0+0.0005,zero,zero,eneutp,efmasp,efmasn)
      call enernm(0.1d0-0.0005,zero,zero,eneutm,efmasp,efmasn)
      deneut = (eneutp-eneutm)/(0.001)
!      deneut = compr0-6D0*dasym0             ! override for test
c      call enernm(0.16d0,zero,zero,eneut1,efmasp,efmasn)
c      write(2,'(/a,2g13.5)')
c     &    ' E/A in neutron matter at rho=0.10,slop:',eneut0,deneut
c
c     now compute variances
C
      IF(I.EQ.0) THEN
       EOVERA = EOVRA0
       RHONM  = RHONM0
       COMPR  = COMPR0
       EFFMAS = EFMAS0
       SYMENR = SYMEN0
       skappa = skapp0
       dsymme = dasym0
       d2symme = d2asy0
       eneut00 = eneut0
       deneut00 = deneut
       b40=t4/2D0
       b4p0=b4p
       vpp0=pairfp
       vpn0=pairfn
       rho0pr0=rho0pr
       tto0=tto
       tte0=tte
      ELSE
       DEOVRA(I) = (EOVRA0-EOVERA)/DPARMO(I)
       DRHONM(I) = (RHONM0-RHONM )/DPARMO(I)
       DCOMPR(I) = (COMPR0-COMPR )/DPARMO(I)
       DEFMAS(I) = (EFMAS0-EFFMAS)/DPARMO(I)
       DSYMEN(I) = (SYMEN0-SYMENR)/DPARMO(I)
       DDSYME(I) = (dasym0-dsymme)/DPARMO(I)
       DD2SYME(I) = (d2asy0-d2symme)/DPARMO(I)
       dskapp(i) = (skapp0-skappa)/dparmo(i)
       DENEUT0(i)= (eneut00-eneut0)/dparmo(i)
       DDENEUT(i)= (deneut00-deneut)/dparmo(i)
       db4(i) = (t4/2D0-b40)/dparmo(i)
       db4p(i) = (b4p-b4p0)/dparmo(i)
       dvpp(i) = (pairfp-vpp0)/dparmo(i)
       dvpn(i) = (pairfn-vpn0)/dparmo(i)
       drho0pr(i) = (rho0pr-rho0pr0)/dparmo(i)
       dtto(i) = (tto-tto0)/dparmo(i)
       dtte(i) = (tte-tte0)/dparmo(i)
      ENDIF
C
      IF(I.GT.0) PAR(I) = PAR(I) - DPARMO(I)
   10 CONTINUE
c
c     write protokol of derivatives (CONTINUED IN EXTOBS)
c
      open(81,file='erroranalysis')
      write(81,'(a,i5)') 'nparmx=',nparmx
      WRITE(81,'(a/100a)') 'model parameters:',
     &    (zuord(iparva(j))//' ',j=1,nparmx)
      WRITE(81,'(a/100(1pg16.8))') 'values:',par(1:nparmx)
      WRITE(81,'(a/100(1pg16.8))') 'variation:',dparmo(1:nparmx)
      WRITE(81,'(a)') 'covariance matrix "zarray":'
      do i=1,nparmx
        write(81,'(100(1pg24.16))') (zarray(j,i),j=1,nparmx)
      enddo
!      write(81,'(a)') 
!     &  ' E/A rhonm  K  m*/m  asym  dasym kappa  E/Aneut slope'
!        write(81,'(100(1pg24.16))') 
!     &    EOVRA0,RHONM0,COMPR0,EFMAS0,SYMEN0,dasym0,skapp0,
!     &    eneut00,deneut00
!      do i=1,nparmx
!        write(81,'(100(1pg24.16))') 
!     &    DEOVRA(I),DRHONM(I),DCOMPR(I),DEFMAS(I),DSYMEN(I),
!     &    DDSYME(I),Dskapp(i),DENEUT0(I),DDENEUT(I)
!      enddo
      WRITE(81,'(a)') ' E/A-sym'
      WRITE(81,'(1pg24.16)') eovra0
      WRITE(81,'(100(1pg24.16))') deovra(1:nparmx)
      WRITE(81,'(a)') ' rho0-sym'
      WRITE(81,'(1pg24.16)') rhonm0
      WRITE(81,'(100(1pg24.16))') drhonm(1:nparmx)
      WRITE(81,'(a)') ' incompr. K'
      WRITE(81,'(1pg24.16)') compr0
      WRITE(81,'(100(1pg24.16))') dcompr(1:nparmx)
      WRITE(81,'(a)') ' m*/m'
      WRITE(81,'(1pg24.16)') efmas0
      WRITE(81,'(100(1pg24.16))') defmas(1:nparmx)
      WRITE(81,'(a)') ' a_sym'
      WRITE(81,'(1pg24.16)') symen0
      WRITE(81,'(100(1pg24.16))') dsymen(1:nparmx)
      WRITE(81,'(a)') ' d_rho a_sym'
      WRITE(81,'(1pg24.16)') dasym0
      WRITE(81,'(100(1pg24.16))') ddsyme(1:nparmx)
      WRITE(81,'(a)') ' d2_rho a_sym'
      WRITE(81,'(1pg24.16)') d2asym0
      WRITE(81,'(100(1pg24.16))') dd2syme(1:nparmx)
      WRITE(81,'(a)') ' kappa'
      WRITE(81,'(1pg24.16)') skapp0
      WRITE(81,'(100(1pg24.16))') dskapp(1:nparmx)
      WRITE(81,'(a)') ' d_rho E/A-neut'
      WRITE(81,'(1pg24.16)') deneut00
      WRITE(81,'(100(1pg24.16))') ddeneut(1:nparmx)
      WRITE(81,'(a)') ' E/A-neut'
      WRITE(81,'(1pg24.16)') eneut00
      WRITE(81,'(100(1pg24.16))') deneut0(1:nparmx)
      WRITE(81,'(a)') ' b4'
      WRITE(81,'(1pg24.16)') b40
      WRITE(81,'(100(1pg24.16))') db4(1:nparmx)
      WRITE(81,'(a)') ' b4p'
      WRITE(81,'(1pg24.16)') b4p0
      WRITE(81,'(100(1pg24.16))') db4p(1:nparmx)
      WRITE(81,'(a)') ' vpairp'
      WRITE(81,'(1pg24.16)') vpp0
      WRITE(81,'(100(1pg24.16))') dvpp(1:nparmx)
      WRITE(81,'(a)') ' vpairn'
      WRITE(81,'(1pg24.16)') vpn0
      WRITE(81,'(100(1pg24.16))') dvpn(1:nparmx)
      WRITE(81,'(a)') ' rho0pr'
      WRITE(81,'(1pg24.16)') rho0pr0
      WRITE(81,'(100(1pg24.16))') drho0pr(1:nparmx)
      IF(ABS(tto0)>1D-5) THEN
        WRITE(81,'(a)') ' tto'
        WRITE(81,'(1pg24.16)') tto0
        WRITE(81,'(100(1pg24.16))') dtto(1:nparmx)
        WRITE(81,'(a)') ' tte'
        WRITE(81,'(1pg24.16)') tte0
        WRITE(81,'(100(1pg24.16))') dtte(1:nparmx)
      END IF
C
C     ACCUMULATE EXTRAPOLATION ERRORS BY FOLDING WITH 'ZARRAY'.
C
      ERRRHO = 0D0
      ERRCMP = 0D0
      ERRMAS = 0D0
      ERRENR = 0D0
      ERRSYM = 0D0
      errdsy = 0D0
      errd2sy = 0D0
      errska = 0D0
      errent = 0D0
      errden = 0D0
C
      IF(NPARMX.GT.0) THEN
       DO I=1,NPARMX
       DO J=1,NPARMX
         ERRRHO = DRHONM(I)*DRHONM(J)*ZARRAY(I,J) + ERRRHO
         ERRCMP = DCOMPR(I)*DCOMPR(J)*ZARRAY(I,J) + ERRCMP
         ERRMAS = DEFMAS(I)*DEFMAS(J)*ZARRAY(I,J) + ERRMAS
         ERRENR = DEOVRA(I)*DEOVRA(J)*ZARRAY(I,J) + ERRENR
         ERRSYM = DSYMEN(I)*DSYMEN(J)*ZARRAY(I,J) + ERRSYM
         ERRDSY = DDSYME(I)*DDSYME(J)*ZARRAY(I,J) + ERRDSY
         ERRD2SY = DD2SYME(I)*DD2SYME(J)*ZARRAY(I,J) + ERRD2SY
         ERRSka = DSkapp(I)*DSkapp(J)*ZARRAY(I,J) + ERRSka
         errent = DENEUT0(I)*DENEUT0(I)*ZARRAY(I,J) + errent
         errden = DDENEUT(I)*DDENEUT(I)*ZARRAY(I,J) + errden
       enddo
       enddo
       ERRRHO = SQRT(ERRRHO)
       ERRCMP = SQRT(ERRCMP)
       ERRMAS = SQRT(ERRMAS)
       ERRENR = SQRT(ERRENR)
       ERRSYM = SQRT(ERRSYM)
       ERRDSY = SQRT(ERRDSY)
       ERRD2SY = SQRT(ERRD2SY)
       ERRSka = SQRT(ERRSka)
       errent = SQRT(errent)
       errden = SQRT(errden)
      ENDIF
C
C     FINAL PRINT
C
      write(02,'(a/a,a/1x,f9.4,f9.3,f9.1,f9.4,2f9.2,f9.4,f9.2)')
     &     ' the nuclear matter parameters:',
     &     '      rho      e/a     compr.  effmas     a4    d/drho_a4',
     &     '    skappa     K_sym',
     &     rhonm0,eovera,compr,effmas,symenr,dsymme,skappa,d2symme
      if(nparmx.gt.0) write(02,'(1x,f9.4,f9.2,f9.1,f9.4,2f9.2,f9.4,
     &                           f9.2)')
     &     errrho,errenr,errcmp,errmas,errsym,errdsy,errska,errd2sy
!      write(22,'(a,f9.4,f9.2,f9.1,f9.4,2f9.2,f9.4)')
!     &     commen(1:20),
!     &     rhonm0,eovera,compr,effmas,symenr,dsymme,skappa
      write(02,'(a,2(1pg12.4)/a,2(1pg12.4))') 
     &     ' Neutron E/A(0.1)  :',eneut00,errent,
     &     '         slope(0.1):',deneut00,errden
c
c     separate energy contributions
c
      if(ioptsk.eq.5 .or. ioptsk.eq.6) then
       call sepenu(rhonm0,half,zero,eovera,enr0sc,enr0vc,enr1sc,
     &              enr1vc,enr3sc,enr3vc,ekinA)
       call sepenu(0.16d0,zero,zero,eneut1,enm0sc,enm0vc,enm1sc,
     &              enm1vc,enm3sc,enm3vc,ekinmA)
       write(2,'(/a/2a)')
     & " the separate contribution to the energy in infinite matter:",
     & "   typ       b_0    b'_0     b_1    b'_1     b_3    b'_3  ",
     & "   kin    total"
       write(2,'(1x,a,1x,8f8.2)')
     & ' symm. ',enr0sc,enr0vc,enr1sc,enr1vc,enr3sc,enr3vc,ekinA,eovera,
     & 'neutron',enm0sc,enm0vc,enm1sc,enm1vc,enm3sc,enm3vc,ekinmA,eneut1
c      else
c       write(2,'(a)') 
c     &  ' separate energy contributions not yet available'
      endif
c
c     print EoS
c

c
c     final call to print the Landau parameters
c
c      call nucmat(rhonm0,eovra0,compr0,efmas0,f0,f0p,g0,g0p,
c     &            symen0,dasymm,d2asym,skapp0,1,-1)
c
      return
      end

c-----nucvar2----------------------------------------------------------
c
      subroutine nucvar2(par,nparmx,rhonm,procon,energnm,errenm)
c
c     computes extrapolated errors on nuclear matter energy
c
c      par    = array of model parameters
c      nparmx = number of varied parameters
c               0 --> no evaluation of extrapolated errors
c      rhonm  = average n.m. density at which E/A is evaluated
c      procon = proton/neutron relation
c      energnm = resulting E/A
c      errenm = error on E/A
c
      include 'params.inc'
      include 'constn.inc'
c      parameter (zero=0D0)
c
      include 'comfit.inc'
      include 'skypar.inc'
c
      DIMENSION denergnm(kparmx)
      dimension par(kparmx)
c
c----------------------------------------------------------------------
c
      DO I=0,NPARMX
        IF(I.GT.0) PAR(I) = PAR(I) + DPARMO(I)
        CALL ATTACH(PAR,JP)
        call initfc
        call enernm(rhonm,procon,zero,energnm0,efmasp,efmasn)
c
c     now compute variances
C
        IF(I.EQ.0) THEN
         energnm = energnm0
        ELSE
         denergnm(I) = (energnm0-energnm)/DPARMO(I)
        ENDIF
C
        IF(I.GT.0) PAR(I) = PAR(I) - DPARMO(I)
      ENDDO
C
C     ACCUMULATE EXTRAPOLATION ERRORS BY FOLDING WITH 'ZARRAY'.
C
      ERRenm = 0D0
      IF(NPARMX.GT.0) THEN
       DO I=1,NPARMX
       DO J=1,NPARMX
        errenm = Denergnm(I)*Denergnm(J)*ZARRAY(I,J) + errenm
       enddo
       enddo
       ERRenm = SQRT(ERRenm)
      ENDIF
c
      return
      end
c-----asymrho------- part of the Hartree-Fock package ------------------
c
      double precision function asymrho(rhonm)
c
c     Computes the symmtry energy of symmetric nuclear matter
c     for given density 'rhonm'.
c 
      include 'params.inc'
      include 'skypar.inc'
      include 'numpar.inc'
c      include 'grid.inc'
c      include 'wavfun.inc'
      include 'constn.inc'
c
c------------------------------------------------------------------------
c
      call enernm(rhonm,half-0.001,zero,energ2,efmasp,efmasn)
      call enernm(rhonm,half,zero,energ0,efmasp,efmasn)
      call enernm(rhonm,half+0.001,zero,energ1,efmasp,efmasn)
      asymrho = (energ1-2D0*energ0+energ2)/(0.001*0.001)/8D0
c
      return
      end



c-----sepenu------------------------------------------------------------
c
      subroutines epenu(rhonm,procon,temper,energm,enr0sc,enr0vc,enr1sc,
     &                  enr1vc,enr3sc,enr3vc,ekinA)     
c
      include 'params.inc'
c
c     soubroutine evaluating the separate energy contributions
c     in the case of symetric nuclear matter:
c      rhonm  = input: total density
c      procon = input: proton content (0.5 for symm.matter)
c      temper = temperature
c      energm = output: E/A in MeV
c      enr0sc = b_0 part
c      enr0vc = b'_0 part
c      enr1sc = b_1 part
c      enr1vc = b'_1 part
c      enr3sc = b_3 part
c      enr3vc = b'_3 part
c      ekinA  = kinetic energy
c     The force parameters are entered via common /skypar/
c
      include 'skypar.inc'
      include 'wavfun.inc'
      include 'numpar.inc'
      include 'constn.inc'
c
      data betab /5.742467999/     ! (3/5)*(3*pi^2)^*(2/3)  factor for tau_n/p
      data f2d3  /0.666666667/
c
c-----------------------------------------------------------------------
c
      rhop   = procon*rhonm
      if(temper.eq.zero) then
       taup   = betab*rhop**f2d3*rhop
      else
       call hotnm(rhop,taup,temper)
      endif
c      if(ifzpe.ne.1) then
c        zpefac = one/(one-one/anumb)
c      else
        zpefac = one
c      endif
      rhopow = rhonm**power *half
      if(rhoswm.gt.zero) then
        tanfac = tanh(rhonm/rhoswm)
      else
        tanfac = 1.0
      endif
      if(procon.eq.half) then
        taunm  = taup+taup
        enr0sc = half*t0a*rhonm*rhonm        
        enr0vc = -half*half*t0b*rhonm*rhonm
        enr1sc = t1p2*rhonm*taunm*tanfac
        enr1vc = -half*t1m2*rhonm*taunm*tanfac
        enr3sc = rhopow*t3a*rhonm*rhonm
        enr3vc = -rhopow*half*t3b*rhonm*rhonm
        ekinA  = zpefac*(amfc(1)+amfc(2))*taup
        if(exmass.ne.zero) enr3vc = 2.0d0*exfunc(rhop)
        if(ioptsk.eq.2 .or. ioptsk.eq.3)    
     &        enr3sc = dirfun(rhonm)   
      else
        rhon   = rhonm-rhop
        taun   = betab*rhon**f2d3*rhon
        taunm  = taup+taun
        enr0sc = half*t0a*rhonm*rhonm
        enr0vc = -half*t0b*(rhop*rhop+rhon*rhon)
        enr1sc = t1p2*rhonm*taunm
        enr1vc = -t1m2*(rhop*taup+rhon*taun)
        enr3sc = rhopow*t3a*rhonm*rhonm
        enr3vc = -rhopow*t3b*(rhop*rhop+rhon*rhon)
        ekinA  = zpefac*(amfc(1)*taup+amfc(2)*taun)
        if(exmass.ne.zero) enr3vc = exfunc(rhop)+exfunc(rhon)
        if(ioptsk.eq.2 .or. ioptsk.eq.3)    
     &        enr3sc = dirfun(rhonm) 
      endif
c
      enr0sc = enr0sc/rhonm                
      enr0vc = enr0vc/rhonm     
      enr1sc = enr1sc/rhonm               
      enr1vc = enr1vc/rhonm
      enr3sc = enr3sc/rhonm
      enr3vc = enr3vc/rhonm
      ekinA  = ekinA/rhonm
      energm = ekinA+enr0sc+enr0vc+enr1sc+enr1vc+enr3sc+enr3vc
c
      return
      end
c-----ran1  ------------ part of hafo package -------------------------------
c
      function ran1(idum)
c
c     random number generator with reshuffling two linear congruential
c     generators, see 'press et al, numerical recipes', sect. 7.1.
c
      include 'params.inc'
c
      parameter (m1=259200,ia1=7141,ic1=54773,rm1=1.0/m1)
      parameter (m2=134456,ia2=8121,ic2=28411,rm2=1.0/m2)
      parameter (m3=243000,ia3=4561,ic3=51349)
c      integer*4 m,ia,ic,idum,iy,ir
      dimension r(97)
      data ix1,ix2,ix3/3*0/,r/97*0D0/
      data iff/0/
c                                        initialization part
      if(idum.lt.0 .or. iff.eq.0) then
        iff   = 1
        ix1   = mod(ic1-idum,m1)
        ix1   = mod(ia1*ix1+ic1,m1)
        ix2   = mod(ix1,m2)
        ix1   = mod(ia1*ix1+ic1,m1)
        ix3   = mod(ix1,m3)
        do 11 j=1,97
        ix1   = mod(ia1*ix1+ic1,m1)
        ix2   = mod(ia2*ix2+ic2,m2)
        r(j)  = (ix1+ix2*rm2)*rm1
   11   continue
        idum  = 1
      endif
c                                        compute random number
      ix1   = mod(ia1*ix1+ic1,m1)
      ix2   = mod(ia2*ix2+ic2,m2)
      ix3   = mod(ia3*ix3+ic3,m3)
      j     = 1 + (97*ix3)/m3
      if(j.gt.97 .or. j.lt.1) stop ' error in ran1'
      ran1  = r(j)
      r(j)  = (ix1+ix2*rm2)*rm1
c
      return
      end

c-----esimple----------------------------------------------------------
c
      subroutine esimple()
c
c     simple model for energy: E = (E/A) * A
c     the volume energy E/A is transferred as parameter t0
c     the surface energy is transferred as parameter t1
c     the symmetry energy is transferred as x_0
c
      include 'params.inc'
c
      include 'wavfun.inc'
      include 'skypar.inc'
      include 'result.inc'

      data ac/0.691/               !  coefficient of Coulomb energy 

      etotal = t0*anumb+t1*anumb**0.66666666667D0
     &        +t2*anumb**0.333333333D0
     &        +x0*(nprot-nneut)**2/anumb
      etotal = etotal + ac*nprot**2/anumb**0.333333333D0
c      WRITE(6,'(a,3(1pg13.5))') ' ESIMPLE:',t0,anumb,etotal

      return
      end

!
      INTEGER FUNCTION iposnuc(nuca,nucz)
!
!     find the position of a nucleus A=nuc1, Z=nuc2 in the table
!     of fit nuclei
!
      include 'params.inc'
      include 'comfit.inc'
!
      iposnuc = 0
      DO inuc=1,npts  
!        WRITE(2,'(a,5i4))') 'iposnuc:',
!     &      inuc,iay(inuc),iprot(inuc),nuca,nucz
        IF(iay(inuc).eq.nuca .and. iprot(inuc).eq.nucz) THEN
          iposnuc = inuc
          EXIT
        END IF
      END DO
      RETURN
      END


c-----handshake--------------------------------------------------------
c
      subroutine handshake05(iu,nprota,nneuta,ngrida,rstepa)
c
c     Prepare input files for subsequent call of 1D HF program.
c     Requires that all detailed force and model parameters are
c     already defined in commons.
c
      include 'params.inc'
c
      include 'numpar.inc'
      include 'wavfun.inc'
      include 'comfit.inc'
      include 'skypar.inc'
      include 'result.inc'
      include 'constn.inc'
c
!      CHARACTER(1),INTENT(IN) :: ident
c
c-----------------------------------------------------------------------
c
!      WRITE(2,'(a)') 'HANDSHAKE: '//ident//'/for005'
c
c     the mean field part
c
        write(iu,'(a/a/i4,f8.3,a)')
     &   '1,,    IWUNIT',
     &   "' 1D SHF input'",
     &   ngrida,rstepa,'     NGRID,RSTEP'
      WRITE(iu,'(i4,2(1pg13.5),2i7,i4,a)')
     &   iwr,dabr,dabe,itpr,itmax,ifplot,
     &   ' -1    iwr,dabr,dabe,itpr,itmax,ifplot'
      WRITE(iu,'(4f9.2,a)')
     &   x0damp,e0damp,addnew,wall,'  x0damp,e0damp,addnew,wall'
        write(iu,'(9(i2,a),i2,a,5x,a)') 
     &    ifex,',',iftm,',',ifenrg,',',ifrso,',',ifzpe,',',
     &    ifdfms,',',ifrhoc,',',ifenav,',',imfold,',',ioptsk,
     &    ', if:ex,tm,enrg,rso,zpe,dfms,rhoc,enav,mfold,optsk'
        write(iu,'(a)')
     &   '1,,,                       ifform',
     &   '0.64,0,0,0,0,   apro0,apro11,apro12,bpro11,bpro12',
     &   '0.64,0,0,0,0,   aneu0,aneu11,aneu12,bneu11,bneu12'
        IF(ioptsk==13) THEN            ! Fayans functional
          WRITE(iu,'(3(1pg20.12),a)') 
     &         avp,h1vp,h2vp,' avp,h1vp,h2vp'
          WRITE(iu,'(5(1pg20.12))') avm,h1vm,h2vm
          WRITE(iu,'(5(1pg20.12))') sigmavol,rho0fy
          WRITE(iu,'(5(1pg20.12))') asp,hsp,hgradp 
          WRITE(iu,'(5(1pg20.12))') asm,hsm,hgradm 
          WRITE(iu,'(5(1pg20.12))') hcoul,sigmacoul
          WRITE(iu,'(5(1pg20.12))') fxi,h1xi,hgradxi,gammaxi,rho0xi
          WRITE(iu,'(4(1pg20.12),i2)') delfxi,ecstab,delh1xi,
     &         delhgradxi,isovectormodel
          WRITE(iu,'(5(1pg20.12))') c0tau,c1tau 
          WRITE(iu,'(5(1pg20.12))') c0nabj,c1nabj,c0tens,c1tens
        ELSE 
          write(iu,'(5(1pg20.12,a)/5(g20.12,a))')
     &    t0,',',t1,',',t2,',',t3,',',t4,',      t0,t1,t2,t3,t4',
     &    x0,',',x1,',',x2,',',x3,',',b4p,',       x0,x1,x2,x3,b4p'
        END IF
        if(ioptsk.eq.5 .or. ioptsk.eq.6 .or. ioptsk.eq.11) then
          write(iu,'(2(1pg20.12,a)/2(1pg20.12,a))')
     &      c1,',',c4,',                      c1 c4',
     &      c1p,',',c4p,',                    c1p,c4p'
        endif
        if(ioptsk.eq.8 .or. ioptsk.eq.11) then
          write(iu,'(3(1pg20.12),a/2(1pg20.12,a))')
     &      b3tilde,tte,tto,',   b3tilde,tte,tto',
     &      b12,',',b12p,',                     b12,b12p'
        endif
        if(ioptsk.eq.9 .or. ioptsk.eq.11) then
          write(iu,'(2(1pg20.12,a)/1pg20.12,a)')
     &      b7,',',b7p,',               b7,b7p',
     &      pow7,',               pow7'
        endif
        if(ioptsk.eq.10 .or. ioptsk.eq.11) then
          write(iu,'(8(1pg20.12),a)')
     &      riftm,rho0pr,powerp,ecstab,alpinv,betinv,
     &      alpin2,betin2,
     &    ',      riftm,rho0pr,powerp,ecstab,alpinv'
        endif
        if(ioptsk.le.4 .or. ioptsk.eq.7) then
          write(iu,'(5(1pg20.12,a)/4(1pg20.12,a))')
     &    t5,',',t6,',',t7,',',t8,',',t9,',    t5,t6,t7,t8,t9',
     &    x5,',',x6,',',x7,',',x8,',           x5,x6,x7,x8'
        endif
        IF(ioptsk.NE.13) 
     &     write(iu,'((1pg20.12,a),2(1pg20.12,a))')
     &      tyuk,',',xyuk,',',ayuk,',           tyuk,xyuk,ayuk'
        write(iu,
     &        '(g21.12,a,i2,a,1pg20.12,a,1pg20.12,a,1pg20.12,f8.2,a)')
     &     power,',1.0,',ipair,',',pairfp,',',pairfn,
     &  ',0,,',cutfac,hfbwidth,
     &  ',,   POWER,ACOUL,IPAIR,PAIRFP/N,ILNRESP,...'
        WRITE(iu,'(4i5,f8.2,a)') 
     &    nprota,nneuta,nmaxp,nmaxn,iparsh*1D0,
     &    '  nprot,nneut,nmaxp,nmaxn,iparsh'
        WRITE(iu,'(a)') '0, ifgsc'
        WRITE(iu,'(a)') '-1, to stop'

c
      return
      end


      include 'sensitivity.f'
      include 'simplex.for'
      include 'powell.for'
      include 'driver.f'
      include 'extobs.f'
      include 'extobs2.f'
      include 'introbs.f'
      include 'givens.f'
      include 'sample_reasonable.f'
      include 'correl-data.f'
