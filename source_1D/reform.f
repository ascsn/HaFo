c
c     reads input from discusson file 'testsky.tex'
c     and reformats that to 'for001' of the hffit package
c
      character*80 a
c
      do 10 iread=1,9999
      read(5,'(a)',err=999,end=999) a
      if(a(1:1).eq.'"') then
        write(6,'(2a)')  '011,',a(1:20)
        read(5,*) t0,t1,t2,t3
        read(5,*) x0,x1,x2,x3
        read(5,*) b4,b4p,alpha,xmass
        read(5,*) iex,ils,icm
c          write for001
        write(6,'(a/4(2a/),a/a,i1,a,i1,a,i1,a/2a,5(/a))')
     &   '0,        NWDR',
     &   '1,  2,  3,  5,  6,  4,  11,   9,  10,0,0,0,0,0,0,0, ',
     &   '       AUSWAHL DER MESSGROESSEN',
     &   '0.2,0.5,1.5,0.1,0.1,0.5,0.01,0.05,0.05,0.,0.,0.,0.,0., ',
     &   'FEHLER AN DATEN',
     &   '0.0,0.0,0.0,1.0,1.0,0.0,1.00, 1.0, 1.0,0.,0.,0.,0.,0.,0,0,0,',
     &   '   % ODER ABSOLUT',
     &   '00.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,',
     &   ' KORRREKTUR ',
     &   '016,1000,0,      IAUNT,IAOB,',
     &   '0,',iex,',',ils,',1,1,',icm,
     &   ',0,0,0,  NMINV,IF:EX,TM,ENRG,RSO,ZPE,GSC,TOND,MFOLD',
     &   '1.0,0.2,.1E-6,.1E-6,1,0,0,   ',
     &   'ACOUL,ADDNEW,DABE,DABR,ipair,ilnrsp,iparsh',
     &   '0.0074,0,0.0099,0,            RSYMAX,NSYMAX,RNTMAX,NNTMAX',
     &   '0,0,0,820,0,                IWRTE,IST1,IFCURV,ITMAX,IEXTAN',
     &   '17, 1,0,2.0,0.001,1,,  NPARIN,NUMAX,IVAR,FREEP,ENDCON,IFTINP',
     &   '0, 5.0        ,16, pairfp',
     &   '0, 5.0        ,17, pairfn'
      write(6,'(1pg14.6,a)')
     &    t0,', 0.20000    , 1,    T0',
     &    t1,', 0.20000    , 2,    T1',
     &    t2,', 0.20000    , 3,    T2',
     &    t3,',  2.0000    , 4,    T3',
     &    x0,', 0.20000E-03, 6,    X0',
     &    x1,', 0.20000E-03, 7,    X1',
     &    x2,', 0.20000E-03, 8,    X2',
     &    x3,', 0.20000E-03, 9,    X3',
     &    2.0*b4,', 0.20000    , 5,    T4',
     &    b4p,', 0.00000E+00,15, B4P',
     &    alpha,', 0.00000E+00,10, POWER'
      write(6,'(a)')
     &  '0.000000E+00, 0.20000    ,11,  TYUK',
     &  '0.000000E+00, 0.20000E-03,12,  XYUK',
     &  '0.000000E+00, 0.00000E+00,13,  AYUK',
     &  '0.000000E+00, 0.00000E+00,14, N.N.'
      endif
   10 continue
c
  999 stop
      end
