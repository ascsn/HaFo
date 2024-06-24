PROGRAM test_folding
    ! Test program for Gaussian folding subroutine.
    use spherical_harmonics
    IMPLICIT NONE
    ! external hfbpot

    INTEGER, PARAMETER :: n_grid = 5
    INTEGER :: i
    REAL(8) :: js
    ! REAL(8), parameter :: pi = 3.141592653589793d0
    ! REAL(8) :: phi
    REAL(8), DIMENSION(n_grid) :: psi_in, r, deltar, deltar2, fpsi_out
    COMPLEX(8) :: Yll
    theta = pi/6.0
    phi = 1.d0 ! phi contribution isn't included
    l = 1
    js = 2

    ! Define the wavefunction to be folded
    r = 0d0
    do i = 1, n_grid
        r(i) = i*0.2d0
    end do
    deltar(1) = 0d0
    deltar2(1) = 0d0


    ! Calculate Yll for l = 1
    Yll = spherical_harmonic(l, l, theta, phi)
    print*, "YLL", Yll

    do i = 2, n_grid
        deltar(i) = r(i) - r(i-1)
        deltar2(i) = deltar(i)**2
        psi_in(i) = r(i)**l*REAL(Yll)*EXP(-r(i)**2/deltar2(i))
    end do
    print*, "Psi_IN", psi_in
    call Gauss_folding_m(psi_in, fpsi_out, 0.2, n_grid, js, 1d0, 4)
    print*, fpsi_out

    ! Test the norm, widths of the new function wrt old function

END PROGRAM test_folding