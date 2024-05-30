PROGRAM lsmooth_torsumr
!
! Reads E and B(E1), computes cumulated strength with erf smoothing
!
IMPLICIT NONE
!
REAL(8), PARAMETER :: sqpi=1.772454D0


INTEGER,PARAMETER :: kstates=199
LOGICAL,PARAMETER :: ttest=.false.
REAL(8),PARAMETER :: width=0.5D0,span=40D0,delen=0.5D0

! bel(1) = T=1 dipole strength
! bel(2) = T=0 compressional strength
! bel(3) = T=1 compressional strength
REAL(8) :: omega(kstates),bel(3,kstates) 
INTEGER :: nstates,nst,idum,nbin,ibin
REAL(8) :: strcum(3,3),dum(11),add(3),energy,omegact


OPEN(1,file='tor_strength.res')
READ(1,'(1x)')
READ(1,'(1x)')
READ(1,'(1x)')

DO nst=1,kstates
  READ(1,*,ERR=99,END=19) omega(nst),bel(2,nst),bel(3,nst),dum,bel(1,nst)
  nstates = nst
  IF(omega(nst)<1D0) bel(:,nst)=0D0    ! remove c.m.
END DO

19 CONTINUE

IF(ttest) THEN
  WRITE(6,*)  '# nstates:',nstates
  WRITE(6,*) '# omega:',omega(1:nstates)
  WRITE(6,*) '# BE:',bel(:,1:nstates)
END IF

! 
! cumulated strengths
!
nbin = NINT(span/delen)

DO ibin=1,nbin
  energy = ibin*delen
  strcum = 0D0
  DO nst=1,nstates
    omegact = omega(nst)
    add = bel(:,nst)*0.5D0*(erf((energy-omegact)/width)+1D0)
    strcum(:,1) = strcum(:,1) + add/omegact
    strcum(:,2) = strcum(:,2) + add
    strcum(:,3) = strcum(:,3) + add*omegact
  END DO
  WRITE(6,'(f8.2,9(1pg13.5))') energy,strcum
END DO

strcum = 0D0
DO nst=1,nstates
  energy = omega(nst)
  strcum(:,1) = strcum(:,1) + bel(:,nst)/energy
  strcum(:,2) = strcum(:,2) + bel(:,nst)
  strcum(:,3) = strcum(:,3) + bel(:,nst)*energy
END DO
WRITE(6,'(a,9(1pg13.5))') '# total:',strcum

!strinteg = 0D0
!DO nst=1,nstates
!  strinteg = strinteg + bel(nst)
!  WRITE(6,'(f10.3,2(1pg13.5))')   omega(nst),strinteg,strinteg/strtotal
!END DO


STOP 'regular end'
99 STOP ' error in input file'

END PROGRAM lsmooth_torsumr
