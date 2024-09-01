c-----input -----------------------------------------------------------
c
      program input
c
c     reformats the simple input 'input' for SHF or RMF into the
c     more detailed input files 'for005.SHF' and 'for005.RMF'.
c     numerical options are all chosen by default values.
c
      character*80 a,blank,title
      data blank/'                                                  '/
c
c----------------------------------------------------------------------
c
c     open files
c
      open(19,form='formatted',status='scratch')
      open(11,form='formatted',file='for005.RMF')
      open(12,form='formatted',file='for005.SHF')
c
c     loop for reading sequences
c
      do 1000 read=1,999
c
c     read terminal selector and nucleus
c
      read(5,*,err=999,end=999) iwunit
      if(iwunit.lt.0) goto 1999
      read(5,*,err=999,end=999) nprot,nneut
c
c     read force
c
      a    = blank
      read(5,*,err=999,end=999) a
      lena  = length(a)
      if(lena.gt.8) then
        write(6,'(a/2a)')
     &  ' the name of the force was to long!',
     &  ' you entered:',a(1:lena)
        stop ' too long name for force'
      endif
c                                                     read heading line in force
      title = a
      open(13,form='formatted',status='old',
     & file='/usr/people/reinhard/HFbase/spher1D/'//a(1:lena)//'.005')
      read(13,'(a)',err=998,end=998) a
      write(6,'(3a)')
     & ' compare force labels:',title(1:lena),a
c                                                     analyze RMF <--> SHF
      read(13,'(a)',err=998,end=998) a
      backspace 13
      do 10 itext=1,80
      if(a(itext:itext).eq.'t') then
        imodel = 1
        goto 19
      elseif(a(itext:itext).eq.'g') then
        imodel = 2
        goto 19
      endif
   10 continue
      stop ' neither SHF nor RMF indicated in force file'
   19 continue
c
c     continue with SHF
c
      if(imodel.eq.1) then
c                                                     read force parameters and opt.
        read(13,*,err=998,end=998) t0,t1,t2,t3
        read(13,*,err=998,end=998) x0,x1,x2,x3
        read(13,*,err=998,end=998) b4,b4p,alpha,amfc
        if(amfc.eq.20.7525) then
          ifdfms = 0
        elseif(amfc.eq.20.750) then
          ifdfms = 1                            ! Tondeur option M_n \neq M_p
        elseif(amfc.eq.20.73553) then
          ifdfms = 2
        elseif(amfc.eq.20.7505) then
          ifdfms = 3                            ! Pearson option, --> special cm
        else
          stop ' this value for inverse mass not yet provided'
        endif
        read(13,*,err=998,end=998) ifex,iftm,ifzpe
        if(ifzpe.eq.2) then
          write(6,'(2a)')
     &     ' IFZPE=2, i.e. variational <P^2>/2m, not yet implemented',
     &     ' --> reset to IFZPE=1'
          ifzpe  = 1
        endif
        if(ifzpe.ne.0 .and. ifzpe.ne.1)
     &   stop ' invalid IFZPE in SHF (only 0 or 1 possible)'
        read(13,*,err=998,end=998) ipair,pairdp,pairdn
        close(13)
c                                                     read pairing
        read(5,*,err=999,end=999) ipair,pairfp,pairfn
        if(ipair.eq.5) then
          if(pairfp.eq.0.0) pairfp = pairdp
          if(pairfn.eq.0.0) pairfn = pairdn
          iparsh = 2
        elseif(ipair.ge.2) then
          iparsh = 2
        else
          iparsh = 0
        endif
        if(ipair.eq.3) then
          ilnrsp = 1
        else
          ilnrsp = 0
        endif
c                                                     print options
        read(5,*,err=999,end=999) iwr,ifplot
c                                                     prepare 'for005.SHF'
        write(12,'(i2,a/3a/a/i2,a,i1,a/a/i2,3(a,i1),a,3(/a))')
     &  iwunit,',                   IWUNIT',
     &  ' ''',title(1:lena),'''',
     &  ' 0,0.3,                  NGRID,RSTEP',
     &  iwr,',1.E-7,1.E-7,278,780,',ifplot,
     &  ',-1,   IWR,DABR,DABE,ITPR,ITMAX,IRAD,IFPLOT,IDISC ',
     &  ' 0.2,60.,0.2,-1.0,,              XODAMP,E0DAMP,ADDNEW,WALL',
     &  ifex,',',iftm,',1,1,',ifzpe,',',ifdfms,
     &   ',0,0,,      IF:EX,TM,ENRG,RSO,ZPE,DFMS,RHOC,ENAV',
     &  ' 1,,                      IFFORM',
     &  ' 0.64,0.0,0.0,0.0,0.0,    APRO0,APR011,APRO12,BPRO11,BPRO12',
     &  ' 0.64,0.0,0.0,0.0,0.0,    ANEU0,ANEU11,ANEU12,BNEU11,BNEU12'
        write(12,'(5(f10.3,a)/5(f10.5,a),3(/a)/
     &             f10.7,a,i1,a,f8.3,a,f8.3,a,i1,a)')
     &   t0,',',t1,',',t2,',',t3,',',2*b4,'  ,t0,t1,t2,t3,t4',
     &   x0,',',x1,',',x2,',',x3,',',b4p,',  x0,x1,x2,x3,b4p',
     &   ' 0.000,0.000,0.000,0.000,0.000,    t5,t6,t7,t8,t9',
     &   ' 0.000,0.000,0.000,0.000000,       x5,x6,x7,x8',
     &   ' 0.000,0.000,0.000,                tyuk,xyuk,ayuk',
     &   alpha,',1.0,',ipair,',',pairfp,',',pairfn,',',ilnrsp,
     &    ', power,acoul,ipair,pairfp,pairfn,ilnrsp'
        write(12,'(i3,a,i3,a,i1,a/a)')
     &   nprot,',',nneut,',0,0,',iparsh,
     &   ',   nprot,nneut,nmaxp,nmaxn,iparsh',
     &   ' 0,     ifgsc'
c
c     continue with rmf
c
      else
        read(13,*,err=998,end=998)
     &        smescp,vmescp,rmescp,dmescp,smsb2c,smsb3c
        read(13,*,err=998,end=998) smesms,vmesms,rmesms,dmesms
        if(dmesms.eq.0.0) dmesms = 980.0
        if(rmesms.eq.0.0) rmesms = 763.0
        read(13,*,err=998,end=998) barms
        smsfn2 = 0.0
        if(title(1:1).eq.'P') read(13,*,err=998,end=998) smsfn2
        read(13,*,err=998,end=998) ifzpe
        if(ifzpe.eq.3) ifzpe = 0
        if(ifzpe.eq.4) ifzpe = -1
        read(13,*,err=998,end=998) ipair,pairdp,pairdn
        close(13)
c                                                     read pairing
        read(5,*,err=999,end=999) ipair,pairfp,pairfn
        if(ipair.eq.5) then
          if(pairfp.eq.0.0) pairfp = pairdp
          if(pairfn.eq.0.0) pairfn = pairdn
          iparsh = 2
        elseif(ipair.ge.2) then
          stop ' IPAIR=2,3 not valid in RMF'
        else
          iparsh = 0
        endif
c                                                     print options
        read(5,*,err=999,end=999) iwr,ifplot
c                                                     prepare 'for005.RMF'
        write(11,'(3a)') '''',title(1:lena),''''
        write(11,'(a/4i5/a/a,7(/a)/2i5,2f9.3,i5,a/a/7i5)')
     &  '  NPROT,NNEUT,NSHLPR,NSHLNE',nprot,nneut,0,0,
     &  '  RSTEP,NGRID','  0.3,   0, ',
     &  '  ENDCON,ADDNEW,ITERMX,REDSTP,MODSTP,NUMSTP',
     &  '   1.E-6,    0.2,  299,   0.7,    1,    1,',
     &  '                         XMSSTP,                  EPSSTP,',
     &  '  938.9,1400.,3800.,600.0, 60.0,   0.2,1.6,0.3,0.004,0.1,',
     &  '                          ENINV,                   ISTEP,',
     &  '        50.,100.,100., 4.0,40.0,               1,2,2,3,1,',
     &  '  IPAIR,ITPAIR,gappr,gapne,IWRIT,IFRSO,IDISC,IWRMAX',
     &     ipair,250,pairfp,pairfn,iwr,',  1,  -1,  10,',
     &  '  IFPAIR(2),IFPLOT,IFENRG,IFZPE,Iparsh,IFEINV',
     &     0,0,ifplot,1,ifzpe,iparsh,1
        write(11,'(a/2f8.2,a/a/4f10.3,5(/a)/8f10.5/a/f4.1,f10.5,a)')
     &  '  BARYMS(2)  ,BRMASS(2)  ,BARMSE(2)  ,',
     &     barms,barms,',  0.0,  0.0,  0.0,  0.0,',
     &  '  SMESMS ,VMESMS ,RMESMS,DMESMS,',smesms,vmesms,rmesms,dmesms,
     &  '  SMASS2,VMASS2,RMASS2,DMASS2,',
     &  '     0.0,   0.0,   0.0,   0.0,',
     &  '  SCPL2  ,VCPL2  ,VTNSFC,RCPL2  ,RTNSFC,DCPL2,SMESB2,SMESB3,',
     &  '     0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,',
     &  '  SMESCP,VMESCP,vtnscp,RMESCP,RTNSCP,DMESCP,SMSB2C,SMSB3C,',
     &     SMESCP,VMESCP,0.0,RMESCP,0.0,DMESCP,SMSB2C,SMSB3C,
     &  '  COULOM,SMSFN2,SPOTFC,SMES0,CUTOFF,',
     &     1.0,smsfn2,', 0.000,  0.0,   0.0,'
        write(11,'(6(a/),a/102a)')
     &  '  INPEFF,IFSCND,IFBRMS,ifnlgs,IFBAG,IFVACP',
     &  '       1,    -1,     1,    0,    0,     0,',
     &  '  SCALPT,VECTPT,OSCAMP,OSCLEN,RADIPT,SURFPT,RADINI',
     &  '    -0.0,   0.0,  0.00,  1.35,   1.3,   0.5,   1.0,',
     &  '  itreps,itrone,epsflm,epsfac,itersw,keepsp',
     &  '       1,    30,  1.03, 1.414,     0,     1,',
     &  '  speini',(',',i=1,102)
      endif
 1000 continue
 1999 continue
      write(11,'(a/a/4i5)')
     &  '''stop''',
     &  '  NPROT,NNEUT,NSHLPR,NSHLNE',-nprot,nneut,0,0
      write(12,'(a)') ' -1,  to stop'
      stop ' reformatting of input succesfully completed'
  999 stop ' error or premature EoF in input file'
  998 stop ' error or premature EoF in force file'
      end
c-----length-----------------------------------------------------------
c
      function length(text)
c
      character*80 text
c
      do 10 l=80,1,-1
      ileng = l
      if(text(l:l).ne.' ') goto 19
   10 continue
      ileng = 0
   19 continue
      length = ileng
c
      return
      end

