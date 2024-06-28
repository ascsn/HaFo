PROGRAM test_folding
    ! Test program for Gaussian folding subroutine.
!    use spherical_harmonics
    IMPLICIT NONE

!   common 'grid' to communicate 'ngrid'
    include 'grid.inc'

!   set the model parameters here by PARAMETER statements
!   (may be changed to be read in later)
    INTEGER, PARAMETER :: n_grid = 31      ! number of radial grid points
    INTEGER, PARAMETER :: js=1             ! angular momentum
    INTEGER, PARAMETER :: jsmax=4          ! maximum angular momentum
    REAL(8), PARAMETER :: rstep_in=0.3D0   ! grid spacing
    REAL(8), PARAMETER :: foldwidth=2.0D0    ! Gaussian folding width
    REAL(8), PARAMETER :: wpwidth=3D0      ! with of test wavepacket

    INTEGER :: i, ll
    ! REAL(8), parameter :: pi = 3.141592653589793d0
    REAL(8) :: t1, t2
    REAL(8), DIMENSION(n_grid) :: psi_in, r, deltar, deltar2, fpsi_out, phi_exact
    REAL(8), DIMENSION(n_grid) :: redPsi_in
    COMPLEX(8) :: Yll
!    theta = pi/6.0
!    phi = 1.d0 ! phi contribution isn't included
    ngrid = n_grid            ! to be available in 'k2psi'
    rstep = rstep_in          ! to be available in 'k2psi'
    ll = js/2  ! orbital angular momentum to combined a.m. quantum number 'js'

    ! Define the wavefunction to be folded
    do i = 1, n_grid
        r(i) = (i-1)*rstep
    end do
!    deltar(1) = rstep
!    deltar2(1) = deltar(1)**2


!    ! Calculate Yll for l = 1
!    Yll = spherical_harmonic(l, l, theta, phi)
!    print*, "YLL", Yll

    do i = 1, n_grid
!        deltar(i) = r(i) - r(i-1)
!        deltar2(i) = deltar(i)**2
        psi_in(i) = r(i)**(ll+1)*EXP(-(r(i)/wpwidth)**2)
    end do
    
    ! SUBROUTINE Gauss_folding_m(psi_in,fpsi_out,rstep,ngrid,js,foldwid,jsmax)
    call Gauss_folding_m(psi_in, fpsi_out, rstep, n_grid, js, foldwidth, jsmax)

    ! Test the norm, widths of the new function wrt old function
    do i = 1, n_grid
        t1 = 1D0/(wpwidth**2+foldwidth**2) !(-deltar2(i)/(2*deltar2(i)))**l*r(i)**(l+1)
        t2  = wpwidth**2*t1
        phi_exact(i) = ((-1)**ll)*t2**(ll+1.5D0)*r(i)**(ll+1)*EXP(-r(i)**2*t1)
!        phi_exact(i) = ((-1)**ll)*t2**(ll)*r(i)**(ll+1)*EXP(-r(i)**2*t1)
    end do

    print*, "phi_exact", phi_exact
    print*, "fpsi", fpsi_out
    print*, phi_exact - fpsi_out

!    redPsi_in = psi_in*r/Yll

    ! Write data to file
    open(unit=10, file="testFoldingData.dat", status="unknown")
    WRITE(10,'(a,2(1pg13.5))') &
     '# test Gaussian folding for w.p.width,fold.width=',wpwidth,foldwidth
    WRITE(10,'(a)') '# r   init w.p.   folded w.p.    exact folding'
    do i = 1, n_grid
        write(10, '(6(1pg13.5))') r(i), psi_in(i), fpsi_out(i),&
        phi_exact(i),fpsi_out(i)/phi_exact(i)
    end do
    close(10)

END PROGRAM test_folding
