USE  NMP_map
IMPLICIT REAL(8) (a-h,o-z)




! test m*/m=1
  rhonm= 0.162272505679
  eovera=  -15.8622266582
  compr= 217.174773932
  asymme= 30.3917468969
  dasymm= 58.9132153373
  effmas= 1.0
  skappa= 0.0
  bmassp= 20.749811
  bmassn= 20.721249

c0tau= 0.0
c1tau= 0.0
  avp= -9.49577
  h1vp= -5.95543/avp
  h2vp= -1.37930/avp 
  avm=12.4287
  h1vm= -17.177052/avm
  h2vm= 240.97291/avm
  sigmavol= 0.33333333
  rho0fy= 0.800000E-01



! with m*/m & kappa
rhonm=  0.164031091053
eovera= -15.8396151174     
compr=  204.705629623    
effmas=  0.863480536454 
asymme=  29.9582023189     
dasymm=  45.5365840450
skappa=  0.396220201041
  bmassp= 20.749811
  bmassn= 20.721249


c0tau= 20.0
c1tau=-30.0
avp= -10.436376    
h1vp= 0.54178836    
h2vp= 0.42032398    
avm=  8.5390869    
h1vm= -2.0862222    
h2vm=  12.523734    
sigmavol= 0.33333333    
rho0fy= 0.80000000E-01




  CALL NMP_Fayans(avpout,h1vpout,h2vpout,avmout,h1vmout,h2vm,sigmavol,&
                      c0tout,c1tout,rho0fy,&
                      rhonm,eovera,compr,effmas,asymme,dasymm,skappa,&
                      bmassp,bmassn)
       WRITE(2,'(a,7(1pg13.5))') 'rho,eoa,comp,efma,asym,dasy,skap=',&
              rhonm,eovera,compr,effmas,asymme,dasymm,skappa
        WRITE(2,'(a,7(1pg13.5))') 'avp,h1vp,h2vp,avm,h1vm,h2vm=',&
         avpout,h1vpout,h2vpout,avmout,h1vmout,h2vm
 

  WRITE(6,'(a,2(1pg13.5))') 'avp=   ',avp,avpout
  WRITE(6,'(a,2(1pg13.5))') 'h1vp=  ',h1vp,h1vpout
  WRITE(6,'(a,2(1pg13.5))') 'h2vp=  ',h2vp,h2vpout
  WRITE(6,'(a,2(1pg13.5))') 'avm=   ',avm,avmout
  WRITE(6,'(a,2(1pg13.5))') 'h1vm=  ',h1vm,h1vmout
  WRITE(6,'(a,2(1pg13.5))') 'c0tau=  ',c0tau,c0tout
  WRITE(6,'(a,2(1pg13.5))') 'c1tau=  ',c1tau,c1tout


STOP
END
