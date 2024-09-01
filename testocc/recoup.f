c-----recoup ----------------------------------------------------------
c
       program recoup
c
        implicit double precision (a-h,o-z)                                   !
c
c   recouples skyrme parameters from t0,.... to b0,....
c   and vice versa.
c   input and output in standard form for fit-program.
c
c     iftinp = switch for input of skyrme parameters
c               1 --> standard form t0,t1,...
c               0 --> recoupled form b0,b0p,...
c
c***********   sample input file, unit 01  **************************
c
c 011,'  walecka-model   ' ,isteur,commen
c 03,        'nwdr'
c 2,1,3,0,0,0,0,0,0,0,0,0,0,0,0,0,        'auswahl der messgroessen'
c 2.0,1.0,5.0,0.0,0.0,0.,0.,0.,0.,0.,0.,0.,0.,0., 'fehler an daten'
c 00.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,  ' % oder absolut
c 00.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., ' korrrektur '
c 010,220,        'iaunt,iaob,
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
      dimension a(17),par(17),dparmo(17),
     &          jp(17),lifit(12),ifit1(12),iparva(17)
      dimension sigmaa(17),sigma0(17),dparm0(17)
c
      character*6 zuord(15,0:2)
      character*80 line
c
      data zuord/'    b0','    b1','    b2','    b3','    t4',
     &           '   b0p','   b1p','   b2p','   b3p',' power',
     &           '  byuk',' byukp','  ayuk',' n.n. ',' n.n. ',
     &           '    t0','    t1','    t2','    t3','    t4',
     *           '    x0','    x1','    x2','    x3',' power',
     *           '  tyuk','  xyuk','  ayuk',' n.n. ',' n.n. ',
     &           '    c0','    b1','    b2','    c3','    t4',
     &           '   c0p','   b1p','   b2p','   c3p',' power',
     &           '  byuk',' byukp','  ayuk',' n.n. ',' n.n. '/
      data a/17*0.01/
      data rhonm/0.16/
c
c--------------------------------------------------------------------
c
       open(unit=1,form='formatted',status='old',file='recoup.inp')
       open(unit=2,form='formatted',status='unknown',file='recoup.out')
c       open(unit=6,form='formatted',status='unknown',file='for006')
c
c    define action from user input
c
      write(6,'(a)')
     & ' give "iftinp" of output (0,1,2):'
      read(5,*) iftinn
c
      if(iftinp.eq.2) then
        rhonm = 0.159
      else
        rhonm = 0.160
      endif
c    initialise auxiliary fields for fit
c    (lifit=12 deactivates the fit-pointer)
c
      do 9001 i=1,12
 9001 lifit(i)=12
c
c      outer loop. for each 'irun' a new computation is started.
c      loop stops if 'isteur<0' is read in.
c
      do 5 irun=1,1000
c
c     preset some control parameters.
c
      iftinp = 1
      zero   = 0.0
c
c**********************************************************************
c
c   first block: reading control parameters, fit options,
c                experimental data, and force parameters.
c
c**********************************************************************
c
      do 110 ir=1,11
      read(1,'(a)',end=9999) line
      write(2,'(a)') line
  110 continue
      read(1,*) nparin,numax,ivar,freep,endcon,iftinp
      write(2,'(1x,3(i2,a),f4.1,a,f6.4,a,i1,a)')
     &       nparin,',',numax,',',ivar,',',freep,',',endcon,',',
     &       iftinn,',    nparin,numax,ivar,freep,endcon,iftinp'
c
c   eingabe der anfangsparameter
c
      do 310 j=1,nparin
      read(1,*) a(j),dparmo(j),iparva(j)
      write(6,'(2g13.6,i6)') a(j),dparmo(j),iparva(j)
      dparm0(j)=dparmo(j)
      par(j)=a(j)
      if(j.gt.numax) dparmo(j)=0.
      do 310 l=1,15
      if(iparva(j).eq.l) jp(l)=j
  310 continue
      if(par(jp(1)).eq.zero .and. par(jp(6)).eq.zero)
     &  write(6,'(a,a)') ' --> version: zero-range',
     &     ' two-body force switched off.'
c
c    store some initial parameter for later protocoll on file 13
c
c**********************************************************************
c
c   second block: performing the fit by repeated calls to the
c                fitting subroutine 'xd1'.
c
c**********************************************************************
c
c     unfold parameters
c
       if(iftinp.eq.1) then
c                    assign t0,... values
        t0     = par(jp(1))
        t1     = par(jp(2))
        t2     = par(jp(3))
        t3     = par(jp(4))
        t4     = par(jp(5))
        x0     = par(jp(6))
        x1     = par(jp(7))
        x2     = par(jp(8))
        x3     = par(jp(9))
        power  = par(jp(10))
        tyuk   = par(jp(11))
        xyuk   = par(jp(12))
        ayuk   = par(jp(13))
c                   recouple to b0,...
        bs0     = t0*(1.0+0.5*x0)
        b0p    = t0*(0.5+x0)
        bs3     = t3*(1.0+0.5*x3)
        b3p    = t3*(0.5+x3)
        a1     = t1*(1.0+0.5*x1)
        a1p    = t1*(0.5+x1)
        a2     = t2*(1.0+0.5*x2)
        a2p    = t2*(0.5+x2)
        bs1     = 0.25*(a1+a2)
        b1p    = 0.25*(a1p-a2p)
        bs2     = 0.125*(3.0*a1-a2)
        b2p    = 0.125*(3.0*a1p+a2p)
        bsyuk   = tyuk*(1.0+0.5*xyuk)
        byukp  = tyuk*(0.5+xyuk)
        b0     = bs0-0.5*b0p
        b1     = bs1-0.5*b1p
        b2     = bs2-0.5*b2p
        b3     = bs3-0.5*b3p
        powfac = rhonm**power/6.0
        c0     = b0+powfac*b3
        c3     = b0-powfac*b3
        c0p    = b0p+powfac*b3p
        c3p    = b0p-powfac*b3p
       else
        if(iftinp.eq.0) then
         b0     = par(jp(1))
         b3     = par(jp(4))
         b0p    = par(jp(6))
         b3p    = par(jp(9))
         power  = par(jp(10))
         powfac = rhonm**power/6.0
         c0     = b0+powfac*b3
         c3     = b0-powfac*b3
         c0p    = b0p+powfac*b3p
         c3p    = b0p-powfac*b3p
        else
         c0     = par(jp(1))
         c3     = par(jp(4))
         c0p    = par(jp(6))
         c3p    = par(jp(9))
         power  = par(jp(10))
         powfac = 3.0/rhonm**power
         b0     = 0.5*(c0+c3)
         b3     = powfac*(c0-c3)
         b0p    = 0.5*(c0p+c3p)
         b3p    = powfac*(c0p-c3p)
        endif
        b1     = par(jp(2))
        b2     = par(jp(3))
        t4     = par(jp(5))
        b1p    = par(jp(7))
        b2p    = par(jp(8))
        byuk   = par(jp(11))
        byukp  = par(jp(12))
        ayuk   = par(jp(13))
c               recouple parameters to t0 etc
        bs0    = b0+0.5*b0p
        bs1    = b1+0.5*b1p
        bs2    = b2+0.5*b2p
        bs3    = b3+0.5*b3p
        fourtr = 4.0/3.0
        t0     = fourtr*(bs0-0.5*b0p)
        x0     = b0p/t0-0.5
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
        tyuk   = fourtr*(bsyuk-0.5*byukp)
        if(tyuk.ne.zero) then
         xyuk   = byukp/tyuk-0.5
        else
         xyuk   = zero
        endif
       endif
c
c      write in new form
c
       if(iftinn.eq.1) then
c                reassign to jp
        par(jp(1))  = t0
        par(jp(2))  = t1
        par(jp(3))  = t2
        par(jp(4))  = t3
        par(jp(5))  = t4
        par(jp(6))  = x0
        par(jp(7))  = x1
        par(jp(8))  = x2
        par(jp(9))  = x3
        par(jp(10)) = power
        par(jp(11)) = tyuk
        par(jp(12)) = xyuk
        par(jp(13)) = ayuk
        dparmo(jp(1))  = 0.2
        dparmo(jp(2))  = 0.2
        dparmo(jp(3))  = 0.2
        dparmo(jp(4))  = 2.0
        dparmo(jp(5))  = 0.2
        dparmo(jp(6))  = 0.0002
        dparmo(jp(7))  = 0.0002
        dparmo(jp(8))  = 0.0002
        dparmo(jp(9))  = 0.0002
        dparmo(jp(11))  = 0.2
        dparmo(jp(12))  = 0.0002
       else
c                   reassign to jp
        if(iftinn.eq.0) then
         par(jp(1))  = b0
         par(jp(4))  = b3
         par(jp(6))  = b0p
         par(jp(9))  = b3p
        else
         par(jp(1))  = c0
         par(jp(4))  = c3
         par(jp(6))  = c0p
         par(jp(9))  = c3p
        endif
        par(jp(2))  = b1
        par(jp(3))  = b2
        par(jp(5))  = t4
        par(jp(7))  = b1p
        par(jp(8))  = b2p
        par(jp(10)) = power
        par(jp(11)) = byuk
        par(jp(12)) = byukp
        par(jp(13)) = ayuk
        dparmo(jp(1)) = 0.2
        dparmo(jp(2)) = 0.2
        dparmo(jp(3)) = 0.2
        dparmo(jp(4)) = 2.0
        dparmo(jp(6)) = 2.0
        dparmo(jp(7)) = 2.0
        dparmo(jp(8)) = 2.0
        dparmo(jp(9)) = 20.0
        dparmo(jp(11)) = 0.2
        dparmo(jp(12)) = 2.0
       endif
c
c      output of resulting parameters in ftn01-input format
c
       do 1019 j=1,nparin
       write(02,'(g13.6,a,g12.5,a,i2,a,a)')  par(j),',',dparmo(j),',',
     &                 iparva(j),',',zuord(iparva(j),iftinn)
 1019  continue
c
    5 continue
c
 9999 stop
      end
