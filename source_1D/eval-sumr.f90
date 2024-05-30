PROGRAM eval_sumr
!
! Reads dipole and coompressional strengths from 'tor_strength.res'
! and computes cumulated strengths (Gaussain folded with 'width') 
! as function of energy.
! Does the same for transition formfactor at two selected 'q'.
! Prints all in columns.
! All strengths are cumulated with E weight and E**(-1) weight.
!
IMPLICIT NONE
!
REAL(8), PARAMETER :: sqpi=1.772454D0


INTEGER,PARAMETER :: kstates=399
LOGICAL,PARAMETER :: ttest=.false.
REAL(8),PARAMETER :: width=1D0,span=40D0,delen=1D0

! bel(1) = T=1 dipole strength
! bel(2) = T=0 compressional strength
! bel(3) = T=1 compressional strength
REAL(8) :: omega(kstates),bel(3,kstates) 
INTEGER :: nstates,nst,idum,nbin,ibin,ngrid,n,nstin
REAL(8) :: strcum(3,2),dum(11),add(3),energy,omegact
REAL(8) :: formf(2,4,kstates),qact,addf(4)

REAL(8),ALLOCATABLE :: strcollect(:,:,:)

INTEGER,PARAMETER :: iq1=10,iq2=15  ! selects q-values -- needs adjustment
REAL(8),ALLOCATABLE :: frmfcollect(:,:,:)
CHARACTER(80) :: line

nbin = NINT(span/delen)
ALLOCATE(strcollect(3,2,0:nbin))

!
! The data from 'tor_strength.res'
!

OPEN(1,file='tor_strength.res')
READ(1,'(1x)')
READ(1,'(1x)')
READ(1,'(1x)')

DO nst=1,kstates
  READ(1,*,ERR=99,END=19) omega(nst),bel(2,nst),bel(3,nst),dum,bel(1,nst)
  nstates = nst
  IF(omega(nst)<1D0) bel(:,nst)=0D0    ! remove c.m.
  IF(ttest) WRITE(6,*) ' tor_strength.res read. nst=',nst
END DO

19 CONTINUE

IF(ttest) THEN
  WRITE(6,*)  '# nstates:',nstates
  WRITE(6,*) '# omega:',omega(1:nstates)
  WRITE(6,*) '# BE:',bel(:,1:nstates)
END IF
CLOSE(1)

! 
! cumulated strengths
!

DO ibin=1,nbin
  energy = ibin*delen
  strcum = 0D0
  DO nst=1,nstates
    omegact = omega(nst)
    add = bel(:,nst)*0.5D0*(erf((energy-omegact)/width)+1D0)
    strcum(:,1) = strcum(:,1) + add/omegact
    strcum(:,2) = strcum(:,2) + add*omegact
  END DO
!  WRITE(6,'(f8.2,9(1pg13.5))') energy,strcum
  strcollect(:,:,ibin) = strcum
END DO

strcum = 0D0
DO nst=1,nstates
  energy = omega(nst)
  strcum(:,1) = strcum(:,1) + bel(:,nst)/energy
  strcum(:,2) = strcum(:,2) + bel(:,nst)*energy
END DO
!WRITE(6,'(a,9(1pg13.5))') '# total:',strcum
strcollect(:,:,0) = strcum

!
! The data from 'transform.res'
!
ALLOCATE(frmfcollect(4,2,0:nbin))
OPEN(1,file='transform.res')

READ(1,'(1x)') 
READ(1,'(8x,i5)') ngrid
IF(ttest) WRITE(*,*) ' input from transform.res with ngrid=',ngrid
READ(1,'(1x)') 


DO nst=1,kstates

  READ(1,'(10x,i5,9x,e13.5,5x,4e13.5)',ERR=98,END=97) &
    nstin,omega(nst),bel(:,nst)
  IF(ttest) WRITE(*,*) ' nst,nstin,omega=',nst,nstin,omega(nst)
  IF(nst /= nstin) STOP ' states do not match'

  DO n=1,iq1-1
   READ(1,'(1x)') 
  END DO
  READ(1,*,ERR=98,END=97) qact,formf(1,:,nst)
  IF(ttest) WRITE(*,*) 'q_1=',qact
  DO n=iq1+1,iq2-1
   READ(1,'(1x)') 
  END DO
  READ(1,*,ERR=98,END=97) qact,formf(2,:,nst)
  IF(ttest) WRITE(*,*) 'q_2=',qact
  DO n=iq2+1,ngrid
   READ(1,'(1x)') 
  END DO
  IF(omega(nst)<1D0) formf(:,:,nst)=0D0          ! supress c.m.
!  READ(1,'(1x)') 
!  READ(1,'(1x)') 

  nstates = nst


  WRITE(7,'(i5,f10.3,4(1pg13.5))')  nst,omega(nst),formf(1,:,nst)
  WRITE(7,'(i5,f10.3,4(1pg13.5))')  nst,omega(nst),formf(2,:,nst)

END DO
STOP ' too many states'
97 CONTINUE
CLOSE(1)

frmfcollect = 0D0
DO ibin=1,nbin
  energy = ibin*delen
  DO nst=1,nstates
    omegact = omega(nst)
    addf(1) = formf(1,3,nst)**2*0.5D0*(erf((energy-omegact)/width)+1D0)
    addf(2) = formf(1,4,nst)**2*0.5D0*(erf((energy-omegact)/width)+1D0)
    addf(3) = formf(2,3,nst)**2*0.5D0*(erf((energy-omegact)/width)+1D0)
    addf(4) = formf(2,4,nst)**2*0.5D0*(erf((energy-omegact)/width)+1D0)
    frmfcollect(:,1,ibin) = frmfcollect(:,1,ibin) + addf/omegact
    frmfcollect(:,2,ibin) = frmfcollect(:,2,ibin) + addf*omegact
  END DO
END DO

!
! now write on file
!
OPEN(1,file='sumr-all.res')
WRITE(1,'(a,i3)') '#nr.bins:',nbin
WRITE(1,'(a)') '# strengths cumulated with inverse energy weight'
WRITE(1,'(a)') '# ibin energy   dipole   compr0    compr1   F_0(0.65)  &
   F_1(0.65)   F_0(1.05)  F_1(1.05)'
DO ibin=1,nbin
  energy = ibin*delen
  WRITE(1,'(i5,f8.2,7(1pg13.5))') &
    ibin,energy,strcollect(:,1,ibin),frmfcollect(:,1,ibin)
END DO
WRITE(1,'(1x)')
WRITE(1,'(1x)')
WRITE(1,'(a)') '# strengths cumulated with direct energy weight'
WRITE(1,'(a)') '# ibin energy   dipole   compr0    compr1   F_0(0.65)  &
   F_1(0.65)   F_0(1.05)  F_1(1.05)'
DO ibin=1,nbin
  energy = ibin*delen
  WRITE(1,'(i5,f8.2,7(1pg13.5))') &
    ibin,energy,strcollect(:,2,ibin),frmfcollect(:,2,ibin)
END DO


CLOSE(1)

! postludium

STOP 'regular end'
99 STOP ' error in input file tor_strength.res'
98 CONTINUE
BACKSPACE(1)
READ(1,'(a)') line
WRITE(*,*) line
STOP ' error in input file transformf.res'

END PROGRAM eval_sumr
