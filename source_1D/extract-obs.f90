! Program to extract final chi**2, nuclear matter properties
! and other observables from 'for002' of the fitting routine.
!
! The program is written in free format and uses Fortran 90.
!
program extract_obs
!
  implicit none
!
!
  integer LENLINE
  parameter (LENLINE=132)
!
  character(LENLINE) :: line
  integer iread
  integer il
  real d1,d2,d3,d4,d5,d6,d7,d8,d9
  real chi2,rho,ea,compr,effmas,asym,dasym,kappa,gmr,gdr,gqr,gdr16
  real skin,levelsn,qalpha,eb120182,ebsn148
  real vchi2,vrho,vea,vcompr,veffmas,vasym,vdasym,vkappa,vgmr,vgdr,vgqr,vgdr16
  real vskin,vlevelsn,vqalpha,veb120182,vebsn148

!
!----------------------------------------------------------------
!
  do iread=1,99999
    read(5,'(a)',err=1,end=2) line
    if(index(line,'subchi').gt.0) then
!      read(5,'(1x)',err=1,end=2)
      read(5,'(a)',err=1,end=2) line
      read(line(14:lenline),*,err=3,end=3) d1,d2,d3,d4,d5,chi2
    elseif(index(line,'effmas').gt.0) then
      read(5,*,err=1,end=2) rho,ea,compr,effmas,asym,dasym,kappa
      read(5,*,err=1,end=2) vrho,vea,vcompr,veffmas,vasym,vdasym,vkappa
    elseif(index(line,'GMR').gt.0) then
      il=index(line,'=')+1
      read(line(il:lenline),*,err=3,end=3) d1,gmr,vgmr
    elseif(index(line,'GDR in 208PB').gt.0) then
      il=index(line,'=')+1
      read(line(il:lenline),*,err=3,end=3) d1,gdr,vgdr
    elseif(index(line,'GQR').gt.0) then
      il=index(line,'=')+1
      read(line(il:lenline),*,err=3,end=3) d1,gqr,vgqr
    elseif(index(line,'n skin 2').gt.0) then
      il=index(line,'=')+1
      read(line(il:lenline),*,err=3,end=3) skin,vskin
    elseif(index(line,'level diff').gt.0) then
      il=index(line,'=')+1
      read(line(il:lenline),*,err=3,end=3) levelsn,vlevelsn
    elseif(index(line,'Q_alpha').gt.0) then
      il=index(line,'=')+1
      read(line(il:lenline),*,err=3,end=3) qalpha,vqalpha
    elseif(index(line,'gy 120').gt.0) then
      il=index(line,'=')+1
      read(line(il:lenline),*,err=3,end=3) eb120182,veb120182
    elseif(index(line,'gy 148').gt.0) then
      il=index(line,'=')+1
      read(line(il:lenline),*,err=3,end=3) ebsn148,vebsn148
    endif
  enddo
  stop ' ATTENTION: line count exhausted'
!
 1    stop ' error while reding file'
!
 2 continue
   write(6,'(2a)') '  chi2    rho    E/A    K     m*/m  asym  dasym   kappa',&
                   '  GMR   GDR   GQR   skin lev.Sn Qalpha E_302  E_148'
   write(6,'(f7.2,f7.3,2f7.2,f6.3,2f7.2,f7.3,3f6.2,3f6.2,2f8.1)') &
    chi2,rho,ea,compr,effmas,asym,dasym,kappa,gmr,gdr,gqr, &
    skin,levelsn,qalpha,eb120182,ebsn148
   write(6,'(f7.2,f7.3,2f7.2,f6.3,2f7.2,f7.3,3f6.2,3f6.2,2f8.1)') &
    vchi2,vrho,vea,vcompr,veffmas,vasym,vdasym,vkappa,vgmr,vgdr,vgqr, &
    vskin,vlevelsn,vqalpha,veb120182,vebsn148
   stop ' normal termination'
!
 3    stop ' error while analyzing line'
!
end
