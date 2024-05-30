PROGRAM extract_obs

! To extract basic observables from 'for002'.
! The programm requires a special 'for010' input
! with only the four nuclei 48Ca, 68Ni, 1332Sn, and 208Pb.

IMPLICIT NONE
INTEGER, PARAMETER :: numnuc=4
INTEGER, PARAMETER :: numrho=20

INTEGER :: na(numnuc),nz(numnuc)
REAL :: ebind(numnuc),rd(numnuc),sigma(numnuc),rms(numnuc),rmsn(numnuc)

REAL :: fermom(numrho),rho(numrho),eossym(numrho),eosn(numrho)

REAL :: rho0,ea,compr,effmas,asym,dasymm,trkappa,eosn01
REAL :: gmr,gdr,gqr
REAL :: pol208,pol132,pol68,pol48
REAL :: skin208,skin132,skin124,skin140,skin68,skin48


CHARACTER*80 :: inline
CHARACTER*8 :: force

INTEGER :: iread,i
REAL :: dum


DO iread=1,99999
  READ(5,'(a)',ERR=99,END=19) inline
  IF(inline(22:28) == 'comment') THEN
    READ(inline,'(29x,a8)',ERR=98,END=98) force
  END IF
  IF(inline(4:17) == 'binding-energy') THEN
    DO i=1,numnuc
      READ(5,*,ERR=98,END=98) na(i),nz(i),dum,ebind(i)
    END DO
  END IF
  IF(inline(4:16) == 'diffr.-radius') THEN
    DO i=1,numnuc
      READ(5,*,ERR=98,END=98) dum,dum,dum,rd(i)
    END DO
  END IF
  IF(inline(4:20) == 'surface thickness') THEN
    DO i=1,numnuc
      READ(5,*,ERR=98,END=98) dum,dum,dum,sigma(i)
    END DO
  END IF
  IF(inline(4:13) == 'rms radius') THEN
    DO i=1,numnuc
      READ(5,*,ERR=98,END=98) dum,dum,dum,rms(i)
    END DO
  END IF
  IF(inline(4:21) == 'neutron rms radius') THEN
    DO i=1,numnuc
      READ(5,*,ERR=98,END=98) dum,dum,dum,rmsn(i)
    END DO
  END IF
  IF(inline(2:12) == 'the nuclear') THEN
    READ(5,'(1x)')
    READ(5,*,ERR=98,END=98) rho0,ea,compr,effmas,asym,dasymm,trkappa
    READ(5,'(15x,f10.3)') eosn01
  END IF
  IF(inline(2:12) == 'EoS in symm') THEN
    READ(5,'(1x)')
    DO i=1,numrho
      READ(5,*,ERR=97,END=97) rho(i),fermom(i),eossym(i)
    END DO
  END IF
  IF(inline(2:12) == 'EoS in neut') THEN
    DO i=1,numrho
      READ(5,*,ERR=97,END=97) dum,eosn(i)
    END DO
  END IF
  IF(inline(2:4) == 'GMR') THEN
    READ(inline,'(48x,f8.3)',ERR=96,END=96) gmr
    READ(5,'(48x,f8.3)',ERR=96,END=96) gdr
    READ(5,'(48x,f8.3)',ERR=96,END=96) gqr
    READ(5,'(1x)')
    READ(5,'(38x,f8.3)',ERR=96,END=96) pol208
    READ(5,'(42x,f8.3)',ERR=96,END=96) skin208
  END IF
  IF(inline(1:22) == 'rms neutron skin 124Sn') THEN
    READ(inline,'(42x,f8.3)',ERR=96,END=96) skin124
    READ(5,'(42x,f8.3)',ERR=96,END=96) skin132
    READ(5,'(42x,f8.3)',ERR=96,END=96) skin140
    READ(5,'(41x,f8.3)',ERR=96,END=96) skin48
    READ(5,'(38x,f8.3)',ERR=96,END=96) pol132
    READ(5,'(37x,f8.3)',ERR=96,END=96) pol48
    READ(5,'(37x,f8.3)',ERR=96,END=96) pol68
    READ(5,'(1x)')
    READ(5,'(41x,f8.3)',ERR=96,END=96) skin68
  END IF
END DO



STOP 'file count exhausted'

19 CONTINUE

WRITE(6,'(a)') ' Summary of results for the Skyrme parameterization:'//force
WRITE(6,'(a)') '   A   Z     E_B  r_diffr sigma r_rms  r_rms,n '
DO i=1,numnuc
  WRITE(6,'(2i4,f10.3,4(f6.3))') &
  na(i),nz(i),ebind(i),rd(i),sigma(i),rms(i),rmsn(i)
END DO
WRITE(6,'(a)')  '  rho_eq   E/A   compr.   m*/m   a_sym  '//&
   'd/drho(a_sym)  TRK-kappa    EoS_neut '
WRITE(6,'(f8.4,f8.3,f7.1,f7.4,2f7.2,f7.4,f7.3)') &
  rho0,ea,compr,effmas,asym,dasymm,trkappa,eosn01
WRITE(6,'(a)')  '       rho        k_f      EoS_symm       EoS_neut '
DO i=1,numrho
  WRITE(6,'(4f12.4)') rho(i),fermom(i),eossym(i),eosn(i)
END DO
WRITE(6,'(a)')  '   GMR      GDR      GQR   in 208Pb'
WRITE(6,'(5f8.3)') gmr,gdr,gqr
WRITE(6,'(a)')  '  polarizabilities in 48Ca, 68Ni, 132Sn, 208 Pb [fm**2/MeV]'
WRITE(6,'(5f8.3)') pol48*8.38,pol68*8.38,pol132*8.38,pol208*8.38
WRITE(6,'(a)') 'rms skins: 48Ca, 68Ni, 124Sn, 132Sn, 140Sn, 208Pb [fm]'
WRITE(6,'(6f8.3)') skin48,skin68,skin124,skin132,skin140,skin208

OPEN(2,file='summary.res')
WRITE(2,'(a8,4f10.3,16f6.3,f8.4,f8.3,f7.1,f7.4,2f7.2,f7.4,f7.3,13f8.3)') &
  force, &
  ebind,rd,sigma,rms,rmsn,rho0,ea,compr,effmas,asym,dasymm,trkappa,eosn01,&
  gmr,gdr,gqr,pol48*8.38,pol68*8.38,pol132*8.38,pol208*8.38, &
  skin48,skin68,skin124,skin132,skin140,skin208
CLOSE(2)

STOP 'regular end'

96 STOP 'error or EoF in reading GR etc'
97 STOP 'error or EoF in reading EoS'
98 STOP 'error or EoF in reading nucleus'
99 STOP 'error in reading text'

END PROGRAM extract_obs


