PROGRAM test_folding
    ! Test program for Gaussian folding subroutine.
    use spherical_harmonics
    IMPLICIT NONE
    ! external hfbpot

    INTEGER, PARAMETER :: n_grid = 15
    INTEGER :: i, js
    ! REAL(8), parameter :: pi = 3.141592653589793d0
    REAL(8) :: t1, t2, rstep
    REAL(8), DIMENSION(n_grid) :: psi_in, r, deltar, deltar2, fpsi_out, phi_exact
    REAL(8), DIMENSION(n_grid) :: redPsi_in
    COMPLEX(8) :: Yll
    theta = pi/6.0
    phi = 1.d0 ! phi contribution isn't included
    l = 1
    js = 2
    rstep = 0.3d0

    ! Define the wavefunction to be folded
    r = 0d0
    do i = 1, n_grid
        r(i) = i*rstep
    end do
    deltar(1) = rstep
    deltar2(1) = deltar(1)**2


    ! Calculate Yll for l = 1
    Yll = spherical_harmonic(l, l, theta, phi)
    print*, "YLL", Yll

    do i = 2, n_grid
        deltar(i) = r(i) - r(i-1)
        deltar2(i) = deltar(i)**2
        psi_in(i) = r(i)**l*REAL(Yll)*EXP(-r(i)**2/deltar2(i))
    end do
    
    ! SUBROUTINE Gauss_folding_m(psi_in,fpsi_out,rstep,ngrid,js,foldwid,jsmax)
    call Gauss_folding_m(psi_in, fpsi_out, rstep, n_grid, js, 3d0, 4)

    ! Test the norm, widths of the new function wrt old function
    do i = 1, n_grid
        t1 = (-deltar2(i)/(2*deltar2(i)))**l*r(i)**(l+1)
        t2  = exp(- r(i)**2/(2*deltar2(i)))
        phi_exact(i) = (-1)**l*(2*deltar2(i))**(3/2)/deltar(i)**3*t1*t2
    end do

    print*, "phi_exact", phi_exact
    print*, "fpsi", fpsi_out
    print*, phi_exact - fpsi_out

    redPsi_in = psi_in*r/Yll

    ! Write data to file
    open(unit=10, file="testFoldingData.dat", status="unknown")
    do i = 1, n_grid
        write(10, '(5es15.6)') r(i), redPsi_in(i), fpsi_out(i),&
        phi_exact(i)
    end do
    close(10)

END PROGRAM test_folding