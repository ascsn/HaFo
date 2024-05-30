USE  NMP_map
IMPLICIT REAL(8) (a-h,o-z)
REAL(8), PARAMETER :: PI=3.14159265358979D0
REAL(8), PARAMETER :: betap=0.6D0*(1.5D0*pi**2)**(2D0/3D0)

! nomass nokappa
  rhonm= 0.1595 
  eovera= -15.921
  compr=     234.5   
  asymme=   31.00 
  dasymm=   47.80
  effmas=  1.0
  skappa=  0.0
  bmassp= 20.749811
  bmassn= 20.721249

    t0= -1814.012
    t1=   271.920
    t2=   -14.079
    t3= 12964.276
    x0=0.21909140   
    x1=-1.2111678   
    x2= 13.235426   
    x3=0.20775555   
    t4=   102.348
   b4p=   52.30659
 power= 0.33


! use SkM* for test
  rhonm= 0.160212158804
  eovera= -15.7522909704
  compr= 216.407689816
  asymme= 30.0356268908
  dasymm= 45.8014713202
  effmas= 0.788772752789
  skappa= 0.530759457139
  bmassp= 20.7525
  bmassn= 20.7525

    t0= -2645.00    
    t1=  410.000    
    t2= -135.000    
    t3=  15595.0    
    x0= 0.900000E-01
    x1=  0.00000    
    x2=  0.00000    
    x3=  0.00000    
    t4=  130.000    
   b4p=  65.0000    
 power= 0.166667    



  CALL mapSk_TtoC(c0rho,c0alpha,c1rho,c1alpha,c0tau,c1tau,&
                  c0delrho,c1delrho,&
                  t0,t1,t2,t3,x0,x1,x2,x3)

  WRITE(6,'(a,5(1pg13.5))') &
   'c0rho,c0alpha,c1rho,c1alpha=',c0rho,c0alpha,c1rho,c1alpha
  WRITE(6,'(a,5(1pg13.5))') &
   'c0tau,c0delrho,c1tau,c1delrho=',c0tau,c0delrho,c1tau,c1delrho
  WRITE(6,'(a,5(1pg13.5))') &
   'test mass:',bmassp/(bmassp+c0tau*rhonm),&
       1D0/effmas-1D0-c1tau*rhonm/bmassp

  rho23=rhonm**0.66666666667D0
  eta=-1D0
  asymmtest=bmassp*betap*rho23*5D0/9D0+c0tau*betap*rho23*rhonm*5D0/9D0 &
            -eta*c1tau*betap*rho23*rhonm*5D0/3D0 &
            +c1rho*rhonm+c1alpha*rhonm**(1D0+power)
  dasymtest=bmassp*betap*rho23*10D0/9D0+c0tau*betap*rho23*rhonm*25D0/9D0 &
            -eta*c1tau*betap*rho23*rhonm*25D0/3D0 &
            +3D0*c1rho*rhonm+3D0*(1D0+power)*c1alpha*rhonm**(1D0+power)

  WRITE(6,'(a,5(1pg13.5))') &
    ' asymme,dasymm=',asymmtest,dasymtest

  CALL mapSk_CtoT(c0rho,c0alpha,c1rho,c1alpha,c0tau,c1tau,&
                  c0delrho,c1delrho,&
                  t0ot,t1ot,t2ot,t3ot,x0ot,x1ot,x2ot,x3ot)

  WRITE(6,'(a,2(1pg13.5))') 't0=   ',t0,t0ot
  WRITE(6,'(a,2(1pg13.5))') 't1=   ',t1,t1ot
  WRITE(6,'(a,2(1pg13.5))') 't2=   ',t2,t2ot
  WRITE(6,'(a,2(1pg13.5))') 't3=   ',t3,t3ot
  WRITE(6,'(a,2(1pg13.5))') 'x0=   ',x0,x0ot
  WRITE(6,'(a,2(1pg13.5))') 'x1=   ',x1,x1ot
  WRITE(6,'(a,2(1pg13.5))') 'x2=   ',x2,x2ot
  WRITE(6,'(a,2(1pg13.5))') 'x3=   ',x3,x3ot



  CALL NMP_Skyrme(t0ot,t1ot,t2ot,t3ot,x0ot,x1ot,x2ot,x3ot,alphaot,&
                      rhonm,eovera,compr,effmas,asymme,dasymm,skappa,&
                      bmassp,bmassn,c0delrho,c1delrho)


  WRITE(6,'(a,2(1pg13.5))') 't0=   ',t0,t0ot
  WRITE(6,'(a,2(1pg13.5))') 't1=   ',t1,t1ot
  WRITE(6,'(a,2(1pg13.5))') 't2=   ',t2,t2ot
  WRITE(6,'(a,2(1pg13.5))') 't3=   ',t3,t3ot
  WRITE(6,'(a,2(1pg13.5))') 'x0=   ',x0,x0ot
  WRITE(6,'(a,2(1pg13.5))') 'x1=   ',x1,x1ot
  WRITE(6,'(a,2(1pg13.5))') 'x2=   ',x2,x2ot
  WRITE(6,'(a,2(1pg13.5))') 'x3=   ',x3,x3ot
  WRITE(6,'(a,2(1pg13.5))') 'power=',power,alphaot

STOP
END
