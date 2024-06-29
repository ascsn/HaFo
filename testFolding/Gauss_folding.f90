SUBROUTINE Gauss_folding_m(psi_in,fpsi_out,rstep,ngrid,js,foldwid,jsmax)

! Subroutine to apply Gaussian folding to the s.p. wavefunction
! 'psi_in' having the combined angular momentum quantum number 'js'.
! The resulting wavefunction is returned on 'fpsi_out'.
! The combined angular momentum 'js' is resolved in to spatial
! angular momentum 'll' and twice total angular momentum 'j2'
! as exemplified below.
! The grid parameters are entered as stepsize 'rstep' and 
! number of grid points 'ngrid'.
! The maximum number of 'js' values as well as the folding width
! 'foldwid' is given for the case of intializaton.
!
! The folding kernel is generated at the first call of this routine
! and stored internally on teh array 'foldwid'. To limit the summations
! for folding an array 'nmin' for lower and another one 'nmax' for
! upper bound is also computed and stored during initialization  of
! the folding kernel.
!
! The routine is written in Fortran 90 and *.f90 foramt. This is not
! compatible with the *.inc files for communcation via common blocks.
! Thus all communcation goes through list parameters.

IMPLICIT NONE

REAL(8), INTENT(IN) :: psi_in(1:ngrid)    ! s.p. wavefunction to be folded
REAL(8), INTENT(OUT) :: fpsi_out(1:ngrid) ! result of folding
REAL(8), INTENT(IN) :: rstep              ! stepsize of radial grid
INTEGER, INTENT(IN) :: ngrid              ! number of grid points
INTEGER, INTENT(IN) :: js                 ! combinded ang.mom. quantum number
REAL(8), INTENT(IN) :: foldwid            ! spatial width of the folding 
INTEGER, INTENT(IN) :: jsmax              ! maximum 'js' supported

REAL(8), PARAMETER :: folding_limit=1D-6  ! threshold above which folding acts
REAL(8), ALLOCATABLE, SAVE :: fold1D(:,:,:) ! array containing the folding matrix
INTEGER, ALLOCATABLE, SAVE :: nmin(:)     ! minimal grid point above threshold
INTEGER, ALLOCATABLE, SAVE :: nmax(:)     ! maximal grid point above threshold
LOGICAL, SAVE :: tfirst=.TRUE.            ! triggers initializaton at first call

INTEGER :: ll      ! spatial angular momentum of 'psi_in'
INTEGER :: j2      ! twice total angular momentum of 'psi_in'
INTEGER :: i

REAL(8), DIMENSION(ngrid) :: hwfout, d2wfin

IF(tfirst) THEN
  IF(jsmax<0) STOP 'folding must be intialized already in that call'
  CALL Init_Gauss_folding
  tfirst = .FALSE.
END IF

ll=js/2
IF(2*ll == js) THEN 
  j2=js-1
ELSE
  j2=js
END IF


!  ...
!      code here the folding
!      optionally exploit limited range in r-space
!
!      Note that the 'psi_in' contains the radial wavefunction R multplied
!      by r, i.e. 'psi_in' <-->  r*R(r). This has to be undone before
!      folding and restored after.
!  ...

fpsi_out = 0d0
do i=1,ngrid
  fpsi_out(i) = SUM(fold1D(i,:,js)*psi_in(:))
end do

CONTAINS

! Initializes the matrix for the folding operation on the radial grid.
! Folding width 'foldwidth' is given as list parameter of the calling
! routine.
 
SUBROUTINE Init_Gauss_folding


integer :: ndimlin, i, la, j, k, lang, l
integer :: xlp, parity
real(8) :: jn
!real(8), dimension(ngrid) :: r
real(8), dimension(SIZE(psi_in), SIZE(psi_in)) :: hwfArray
real(8), ALLOCATABLE :: exp_uv(:), d2psi_uv(:,:), eigVals(:), eigVects(:,:)
real(8), ALLOCATABLE :: psi_u(:)
logical,parameter :: ttestfold=.FALSE.

! initializes folding matrix and selection of active elements
ALLOCATE(fold1D(1:ngrid,1:ngrid,jsmax),nmin(1:ngrid),nmax(1:ngrid))
ALLOCATE(exp_uv(1:ngrid), d2psi_uv(1:ngrid, 1:ngrid))
ALLOCATE(eigVals(1:ngrid), eigVects(1:ngrid, 1:ngrid))
ALLOCATE(psi_u(1:ngrid))



!  ...
!      Define here the folding matrix for each 'js=1,jsmax'.
!      Evaluate the bounds 'nmin' and 'nmax' according to 'folding_limit'.
!
!      The idea is to compute the folding matrix via Fourier-Bessel basis.
!      The basis is generated on the grid using by diagonalizing the
!      matrix of the -nabla**2 operator using 'diag_real_symm' in 'givens.f'.
!      The matrix of -nabla**2 can be composed using the routine 'k2psi' in
!      'hafbody.for' which computes the action of -nabla**2.
!     
!  ...


! Define bessel-fourier grid for the same matrix
ndimlin = ngrid*(ngrid+1)/2
do l = 1, jsmax             ! loop over total angular momenta
  lang=l/2                  ! orbital angular momentum

! compose matrix of -partial_r^2 on the grid points 
  d2psi_uv(:,1)=0D0
  do j = 2, ngrid     
    psi_u=0D0
    psi_u(j) = 1D0

!    Calculate the matrix of -nabla**2
    call k2psi(psi_u, hwfout, d2wfin, l)
    d2psi_uv(:,j) = hwfout(:)
  end do
  IF(ttestfold) THEN
    WRITE(6,'(a,i3)') 'The matrix of -partial_r^2 for js=',js
    do i=1,ngrid
      write(6,'(200(1pg13.5))') d2psi_uv(:,i)
      WRITE(6,'(/1x)')
    end do
  END IF

! diagonalize to have "numerical Bessel functions" on 'eigVects'
  call diag_real_symm(d2psi_uv, eigVals, eigVects, ngrid, ngrid, ndimlin)
  exp_uv = exp(-eigVals*(0.5D0*foldwid)**2)
  IF(ttestfold) THEN
    WRITE(6,'(/a)') 'The eigenvalues of -partial_r^2:'
    write(6,'(200(1pg13.5))') eigVals
    WRITE(6,'(/a,i3)') 'The eigenvectors'
    do i=1,ngrid
      write(6,'(200(1pg13.5))') eigVects(i,:)
!      WRITE(6,'(/1x)')
    end do
  END IF

! finally compose folding matrix 
  do i=1,ngrid; do j = 1, ngrid
!   sum over "Bessel" eigenstates  
    fold1D(i,j,l)=0D0
    do k = 2, ngrid
      fold1D(i,j,l) = fold1D(i,j,l) + eigVects(i,k)*exp_uv(k)*eigVects (j, k)
    end do
  end do; end do
  IF(ttestfold) THEN
    WRITE(6,'(/a,i3)') 'The folding matrix for js=',l
    do i=1,ngrid
      write(6,'(200(1pg13.5))') fold1D(:,i,l)
      write(6,'(1x)')
    end do
  END IF

end do                     ! end loop over total angular momenta

END SUBROUTINE Init_Gauss_folding

! General Bessel basis function
real(8) function bessel_basis(order, r, theta)
    implicit none
    real(8), intent(in) :: r, theta
    integer, intent(in) :: order
    real(8) :: jn

    ! Compute the nth-order Bessel function of the first kind
    call bessel_calc(order, r, jn)
    bessel_basis = jn * cos(theta)
end function bessel_basis

! Subroutine to compute the nth-order Bessel function of the first kind
subroutine bessel_calc(n, x, result)
    implicit none
    integer :: n
    real(8), intent(in) :: x
    real(8), intent(out) :: result

    ! Use a built-in function or an approximation
    result = bessel_jn(n, x)
end subroutine bessel_calc


END SUBROUTINE Gauss_folding_m


SUBROUTINE gapoperator(fpsi_in,gappsi_out,rstep,ngrid,js,pairpt)

! Subroutine to compute the action of the gap operator "Delta".
! Input is the already folded s.p. wavefunction 'fpsi_in'.
! Output is 'gappsi_out'. 
! The 'js' is the combined angular momentum quantum numbero  of
! the actual state. 
! The grid parameters are enterd as stepsize 'rstep' and 
! number of grid points 'ngrid'.
! The local part of the gap potential is enetered via'pairpt'.
!
! The routine is written in Fortran 90 and *.f90 foramt. This is not
! compatible with the *.inc files for communcation via common blocks.
! Thus all communcation goes through list parameters.

IMPLICIT NONE

REAL(8), INTENT(IN) :: fpsi_in(1:ngrid)    ! folded s.p. wavefunction 
REAL(8), INTENT(OUT) :: gappsi_out(1:ngrid) ! result of the action
REAL(8), INTENT(IN) :: rstep              ! stepsize of radial grid
INTEGER, INTENT(IN) :: ngrid              ! number of grid points
INTEGER, INTENT(IN) :: js                 ! combinded ang.mom. quantum number
REAL(8), INTENT(IN) :: pairpt(1:ngrid)    ! local gap potential

REAL(8), ALLOCATABLE :: psiwork(:)

ALLOCATE(psiwork(1:ngrid))

! action of local gap potential
psiwork = pairpt*fpsi_in

! action of the second folding (last two entries dummy)

CALL  Gauss_folding_m(psiwork,gappsi_out,rstep,ngrid,js,1D0,-1)

DEALLOCATE(psiwork)


END SUBROUTINE gapoperator




REAL(8) FUNCTION Delta_state_folded(fpsi_in,rstep,ngrid,js,pairpt)

! Function computes the gap 'Delta_alpha' of a state.
! Input is the already folded s.p. wavefunction 'fpsi_in'.
! The 'js' is the combined angular momentum quantum numbero  of
! the actual state. 
! The grid parameters are enterd as stepsize 'rstep' and 
! number of grid points 'ngrid'.
! The local part of the gap potential is enetered via'pairpt'.
!
! The routine is written in Fortran 90 and *.f90 foramt. This is not
! compatible with the *.inc files for communcation via common blocks.
! Thus all communcation goes through list parameters.
!
! ==> Actually, this routine is trivial. Its call my be replaced by
!     inserting this one instruction into the calling routine.

IMPLICIT NONE

REAL(8), INTENT(IN) :: fpsi_in(1:ngrid)    ! folded s.p. wavefunction 
REAL(8), INTENT(IN) :: rstep              ! stepsize of radial grid
INTEGER, INTENT(IN) :: ngrid              ! number of grid points
INTEGER, INTENT(IN) :: js                 ! combinded ang.mom. quantum number
REAL(8), INTENT(IN) :: pairpt(1:ngrid)    ! local gap potential

INTEGER :: i       ! running index for loop over grid 
REAL(8), ALLOCATABLE :: psiwork(:)

Delta_state_folded = SUM(pairpt*fpsi_in**2)*rstep


END FUNCTION Delta_state_folded
