c-----nameinp ----------------------------------------------------------
c
      subroutine nameinp()
c
c     Reads namelist-format input for hafomn (1D Hartree-Fock code).
c     Called from hafomn when the first line of unit 05 contains '&MAIN'.
c     Reads namelists MAIN, GRIDS, STATIC, FORCE, THERMAL from unit 05.
c
c     Variable mapping (namelist keyword -> common block variable):
c       &MAIN:    title      -> TITLE   (NUMPAR)
c                 iwunit     -> local (opens output unit 06)
c                 itmax,iwr,ifplot,iterbcs,iterhfb,pairenhance,
c                 modnat,iternat,v2limh,delmix  -> NUMPAR (direct)
c       &GRIDS:   ngrid,rstep          -> GRID (direct)
c       &STATIC:  nprot,nneut,nmaxp,nmaxn -> SHELLS (direct)
c                 deocc(2)   -> iparsh (SKYPAR): max(NINT(deocc))
c                 e0inv      -> e0damp (NUMPAR)
c                 x0inv      -> x0damp (NUMPAR)
c                 addrho     -> addnew (NUMPAR)
c                 errenr     -> dabe   (NUMPAR)
c                 errrad     -> dabr   (NUMPAR)
c                 npiter     -> itpr   (NUMPAR)
c                 trestart   -> idisc  (NUMPAR)
c       &FORCE:   t0..t4,x0..x3,b4p,power,riftm,ipair,rho0pr -> SKYPAR
c                 h2m(2)     -> amfc(2) (SKYPAR)
c                 chg        -> acoul   (SKYPAR)
c                 ex         -> ifex    (NUMPAR)
c                 zpe        -> ifzpe   (NUMPAR)
c                 v0prot     -> pairfp  (SKYPAR)
c                 v0neut     -> pairfn  (SKYPAR)
c                 pair_cutoff -> ecstab (SKYPAR)
c                 pairstab   -> ecstab  = -pairstab (SKYPAR, if pair_cutoff=0)
c                 hfbwidth   -> hfbwidth (SKYPAR)
c                 tregpair   -> tregpair (/switches/)
c                 tcutoff_uvspace -> tphpair (/switches/)
c                 tcoulfold  -> ifrhoc  (NUMPAR)
c                 tFayans    -> tpairFy (/faypar/)
c                 isovectormodel -> isovectormodel (/faypar/)
c                 Fayans force parameters -> /faypar/ (direct)
c       &THERMAL: temprt -> temprt (SKYPAR)
c
      include 'params.inc'
      include 'grid.inc'
      include 'numpar.inc'
      include 'skypar.inc'
      include 'wavfun.inc'
c
c     SAVE variable to track output unit state across calls (loop over nuclei)
      integer, save :: iwuold_nml
      data iwuold_nml / -1 /
c
c     Local variables for namelist keywords that differ from common block names
c
      integer :: iwunit_nml
      real(8) :: e0inv, x0inv, addrho, errenr, errrad
      integer :: npiter_nml
      real(8) :: deocc(2)
      logical :: trestart_nml
      character*20 :: name
      real(8) :: h2m(2)
      real(8) :: chg
      integer :: ex_nml, zpe_nml
      real(8) :: v0prot, v0neut
      real(8) :: pair_cutoff, pairstab_nml
      real(8) :: hfbwidth_nml
      logical :: tcoulfold_nml
      logical :: tFayans_nml
      logical :: tcutoff_uvspace
c
c     Namelists
c
      namelist /main/ title, iwunit_nml,
     &     itmax, iwr, ifplot,
     &     iterbcs, iterhfb, pairenhance,
     &     modnat, iternat, v2limh, delmix
c
      namelist /grids/ ngrid, rstep
c
      namelist /static/ nprot, nneut, nmaxp, nmaxn, deocc,
     &     e0inv, x0inv, addrho,
     &     errenr, errrad, npiter_nml, trestart_nml
c
      namelist /force/ name,
     &     t0, t1, t2, t3, t4,
     &     x0, x1, x2, x3, b4p, power,
     &     h2m, chg, ex_nml, riftm, zpe_nml,
     &     ipair, v0prot, v0neut, rho0pr,
     &     pairstab_nml, pair_cutoff, hfbwidth_nml,
     &     tregpair, tcutoff_uvspace, tcoulfold_nml, tFayans_nml,
     &     isovectormodel,
     &     avp, h1vp, h2vp, avm, h1vm, h2vm,
     &     asp, hsp, hgradp, asm, hsm, hgradm,
     &     hcoul,
     &     fxi, h1xi, hgradxi,
     &     delfxi, delh1xi, delhgradxi,
     &     c0tau, c1tau, c0nabj, c1nabj, c0tens, c1tens,
     &     gammaxi, rho0xi, sigmavol, rho0fy
c
      namelist /thermal/ temprt
c
c --- set defaults before reading ---
c
c     &MAIN defaults
      iwunit_nml  = 1
      title       = ' '
c     remaining MAIN variables keep blockdata defaults (itmax,iwr,ifplot etc.)
c
c     &GRIDS defaults (from blockdata: ngrid=0, rstep=0.3)
c     nothing to override; blockdata values are already set
c
c     &STATIC defaults
      deocc(1)    = 0D0
      deocc(2)    = 0D0
      e0inv       = e0damp
      x0inv       = x0damp
      addrho      = addnew
      errenr      = dabe
      errrad      = dabr
      npiter_nml  = itpr
      trestart_nml = .false.
      nmaxp       = 0
      nmaxn       = 0
c
c     &FORCE defaults (use current common block values as defaults)
      name        = 'none'
      h2m(1)      = amfc(1)
      h2m(2)      = amfc(2)
      chg         = acoul
      ex_nml      = ifex
      zpe_nml     = ifzpe
      v0prot      = pairfp
      v0neut      = pairfn
      pairstab_nml = 0D0
      pair_cutoff = ecstab
      hfbwidth_nml = hfbwidth
      tcoulfold_nml = (ifrhoc .ne. 0)
      tFayans_nml = tpairFy
      tcutoff_uvspace = tphpair
c
c --- read all namelist groups ---
c
      read(05, nml=main)
c
c     handle termination request
      if(iwunit_nml .lt. 0) stop ' nameinp: stop requested by user'
c
c     handle output unit 06 (same logic as old-format branch)
      if(iwunit_nml .ne. iwuold_nml) then
        if(iwuold_nml .ne. -1) close(unit=06)
        if(iwunit_nml .eq. 0) then
          open(unit=06, file='for006', status='unknown')
        endif
      endif
      iwuold_nml = iwunit_nml
c
      read(05, nml=grids)
      read(05, nml=static)
      read(05, nml=force)
      read(05, nml=thermal)
c
c --- map local variables to common block variables ---
c
c     NUMPAR
      e0damp = e0inv
      x0damp = x0inv
      addnew = addrho
      dabe   = errenr
      dabr   = errrad
      itpr   = npiter_nml
      if(trestart_nml) then
        idisc = 0
      else
        idisc = -1
      endif
      ifex   = ex_nml
      ifzpe  = zpe_nml
      if(tcoulfold_nml) then
        ifrhoc = 1
      else
        ifrhoc = 0
      endif
c
c     SHELLS: deocc(1/2) = extra proton/neutron shells above minimum
      iparsh = max(nint(deocc(1)), nint(deocc(2)))
c
c     SKYPAR
      amfc(1) = h2m(1)
      amfc(2) = h2m(2)
      acoul   = chg
      pairfp  = v0prot
      pairfn  = v0neut
c     pair_cutoff takes precedence; pairstab is -ecstab convention from output
      if(pair_cutoff .ne. 0D0) then
        ecstab = pair_cutoff
      elseif(pairstab_nml .ne. 0D0) then
        ecstab = -pairstab_nml
      endif
      hfbwidth = hfbwidth_nml
c
c     /switches/
      tphpair = tcutoff_uvspace
c     tregpair is direct (already in NAMELIST from /switches/ common)
c
c     /faypar/
      tpairFy = tFayans_nml
c     isovectormodel and Fayans force parameters are set directly via NAMELIST
c
c     mark that namelist input was used
      tnamelist = .true.
c
      return
      end
