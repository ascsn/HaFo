PROGRAM integstrength
!
! Reads E and B(E1), computes integrated strength
!
IMPLICIT NONE
!
REAL(8), PARAMETER :: sqpi=1.772454D0


INTEGER,PARAMETER :: kstates=199
LOGICAL,PARAMETER :: ttest=.false.

REAL(8) :: omega(kstates),bel(kstates)
INTEGER :: nstates,nst,idum
REAL(8) :: strtotal,strinteg,dum(4)


OPEN(1,file='lsmooth.inp')
READ(1,'(1x)')
READ(1,'(1x)')
READ(1,'(1x)')

DO nst=1,kstates
  READ(1,*,ERR=99,END=19) idum,omega(nst),dum,bel(nst)
  nstates = nst
END DO

19 CONTINUE

IF(ttest) THEN
  WRITE(6,*)  '# nstates:',nstates
  WRITE(6,*) '# omega:',omega(1:nstates)
  WRITE(6,*) '# BE:',bel(1:nstates)
END IF

strtotal = 0D0
DO nst=1,nstates
  strtotal = strtotal + bel(nst)
END DO

strinteg = 0D0
DO nst=1,nstates
  strinteg = strinteg + bel(nst)
  WRITE(6,'(f10.3,2(1pg13.5))')   omega(nst),strinteg,strinteg/strtotal
END DO


STOP 'regular end'
99 STOP ' error in input file'

END
