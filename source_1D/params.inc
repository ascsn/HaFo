C INCLUDE BLOCK: PARAMS
C
C       Contains all settings for dimensions etc.
C       Check lines with ! for desired settings.
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                   !
        PARAMETER (KVARLN=8)                                                  !
C
C       KVARLN = LENGTH OF A VARIABLE,
C                SET TO 4 FOR SINGLE PRECISION IBM,  PC, ...
C                SET TO 8 FOR DOUBLE PRECISION IBM,  PC, ...
C                SET TO 1 FOR SINGLE PRECISION VAX
C                SET TO 2 FOR DOUBLE PRECISION VAX
C                SET TO 1 FOR SINGLE PRECISION CDC
C                SET TO 2 FOR DOUBLE PRECISION CDC
C
C       KG1   = MAX.NUMBER OF GRID POINTS
C       KGR2  = 2 * KG1
C       KNMX  = MAX. NUMBER OF RADIAL NODES 
C       KJMX  = MAX. NUMBER OF JS
C              KJMX and KNMX are also switches to energy initialization:
C              KJMX>47 or KNMX>12 signal that switch
C       KX11  = SIZE OF WF-FIELD = KG1 * KSL2
C       KSLM  = MAX. OCCURING OSCILLATOR SHELL
C       KNJX  = (KNMX + 1) * 2 KJMX
C       KSL2  = MAX. NUMBER OF STATES (PROT. + NEUT.)
C       KFILE = SHIFT IN FILE STORAGE NUMBER
C       KFRM  = LENGTH OF FORM-FACTOR FIELD
C       KQMAX = MAXIMUM OF FORM-FACTOR STEPS
C       KMAXP = MAX.NUMBER OF PROTON SHELLS
C       KMAXN = MAX.NUMBER OF NEUTRON SHELLS
C       KJRPMX= MAX.NUMBER OF 'JS' IN STORAGE OF REDUCED M.E.
C       KRESWF= Size of Arrays for the separable residual interaction:
C               set to 1 if separable interaction is disenabled,
C               set KX11 if you want to use separable interaction.
C       KTIMWF= Size of Arrays for the separable residual interaction:
C               set to 1 if separable interaction is disenabled,
C               set KX11 if you want to use separable interaction.
C       KICHR = size of array ICHR which contains the keys for BUFFIO.
C       KSCDIP= switch to r**3 for isoscalar dipole modes
C       KCMCORR = approximate c.m. subtraction in r**3 mode and toroidal str.
C       TMONOPC = switch to monopole constraint (using 'wall' variable)
c       TLSYMM  = re-scale slope of asymmetry energy --> 3*rho*dasymm
c       TEPAIR  = direct pairing energy for stabilized pairing
c       TPHPAIR = switch to old handling of pairing space (weight only on
c                   Delta, not on v**2)
c       TLOWERCUT = activates also pairing cutoff below Fermi surface
C
      logical TMONOPC,TADSIC,TLSYMM,TFORMW,TEPAIR,TLOWERCUT
      LOGICAL :: tregpair
      LOGICAL :: tphpair 
csmall       PARAMETER (KG1=61,KNMX=7,KJMX=27,KMAXN=56,KMAXP=46,                 !
cmed       PARAMETER (KG1=199,KNMX=11,KJMX=47,KMAXN=56,KMAXP=56,                 !
      PARAMETER (KG1=301,KNMX=11,KJMX=47,KMAXN=277,KMAXP=277,              !
cold       PARAMETER (KG1=99,KNMX=11,KJMX=47,KMAXN=56,KMAXP=56,                 !
!      PARAMETER (KG1=2289,KNMX=11,KJMX=47,KMAXN=277,KMAXP=277,              !
chuge       PARAMETER (KG1=230,KNMX=30,KJMX=140,KMAXN=3000,KMAXP=3000,              !
!      PARAMETER (KG1=199,KNMX=11,KJMX=47,KMAXN=137,KMAXP=121,                 !
     *            KFRM=KG1,KFILE=400,KQMAX=KG1,KLRPA=25,KRPMAX=599,          !
     *            KICHR=KRPMAX,
     *            KIFHFB=1,
     *            KSCDIP=0,KCMCORR=1,
     *            KHFBN= (KMAXN-1)*KIFHFB+1,
     *            KHFBP= (KMAXP-1)*KIFHFB+1,
     *            KNMXP=KNMX+1,KJMX2=2*KJMX,
     *            KNJX=KNMXP*KJMX2,KSL2=KMAXP+KMAXN,KSLM=KMAXN,
     *            KGR2=2*KG1,KGR2P=KGR2+1,KX11=KSL2*KG1,
     *            KJRPMX=KJMX+KLRPA/2,KSTART=KNMXP*KJMX+1,
     *            KRESWF=1,
     *            KTIMWF=KX11,TMONOPC=.false.,TADSIC=.false.,
     &            TLSYMM=.true.,TFORMW=.FALSE.,TEPAIR=.TRUE.,
     &            TLOWERCUT=.FALSE.)
c
c      parameters for switches in the RPA part
c
      logical oddoff,evenoff
      parameter (oddoff=.false.,evenoff=.false.)
!
      COMMON /switches/tregpair,tphpair
C END INCLUDE BLOCK


