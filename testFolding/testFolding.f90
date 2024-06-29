PROGRAM test_folding
    ! Test program for Gaussian folding subroutine.
!    use spherical_harmonics
    IMPLICIT NONE

!   common 'grid' to communicate 'ngrid'
    include 'grid.inc'

!   set the model parameters here by PARAMETER statements
!   (may be changed to be read in later)
    INTEGER, PARAMETER :: n_grid = 31      ! number of radial grid points
    ! INTEGER, PARAMETER :: js=2             ! angular momentum
    INTEGER, PARAMETER :: jsmax=15          ! maximum angular momentum
    INTEGER :: js
    REAL(8), PARAMETER :: rstep_in=0.3D0   ! grid spacing
    ! REAL(8), PARAMETER :: foldwidth=2.0D0    ! Gaussian folding width
    ! REAL(8), PARAMETER :: wpwidth=3D0      ! with of test wavepacket
    REAL(8) :: foldwidth, wpwidth

    INTEGER :: i, ll, k
    INTEGER :: fw, wp
    ! REAL(8), parameter :: pi = 3.141592653589793d0
    REAL(8) :: t1, t2           
    REAL(8), DIMENSION(n_grid) :: psi_in, r, deltar, deltar2, fpsi_out, phi_exact
    REAL(8), DIMENSION(n_grid) :: redPsi_in
    COMPLEX(8) :: Yll
    CHARACTER(len=100) :: filename
    CHARACTER(*), parameter :: filepath="./results_data/"
!    theta = pi/6.0
!    phi = 1.d0 ! phi contribution isn't included
    ngrid = n_grid            ! to be available in 'k2psi'
    rstep = rstep_in          ! to be available in 'k2psi'

    ! Define a loop over js, foldwidth, wpwidth
    jsLoop: do js = 1, jsmax
    fwLoop: do fw = 1, 5
        foldwidth = 0.3d0 + dble(fw - 1)*0.3d0
    wpLoop: do wp = 1, 8
        wpwidth = 1d0 + dble(wp-1)*0.5d0
        ll = js/2  ! orbital angular momentum to combined a.m. quantum number 'js'

        ! Define the wavefunction to be folded
        do i = 1, n_grid
            r(i) = (i-1)*rstep
        end do

        do i = 1, n_grid
            psi_in(i) = r(i)**(ll+1)*EXP(-(r(i)/wpwidth)**2)
        end do
        
        ! SUBROUTINE Gauss_folding_m(psi_in,fpsi_out,rstep,ngrid,js,foldwid,jsmax)
        call Gauss_folding_m(psi_in, fpsi_out, rstep, n_grid, js, foldwidth, jsmax)

        ! Test the norm, widths of the new function wrt old function
        do i = 1, n_grid
            t1 = 1D0/(wpwidth**2+foldwidth**2) 
            t2  = wpwidth**2*t1
            phi_exact(i) = ((-1)**(2*ll))*t2**(ll+1.5D0)*r(i)**(ll+1)*EXP(-r(i)**2*t1)
        end do

        ! print*, "phi_exact", phi_exact
        ! print*, "fpsi", fpsi_out
        ! print*, phi_exact - fpsi_out


        ! Write data to file
        k = js + 10
        WRITE(filename, '("test_js_",I2.2,"_fw_",f3.1,"_wp_",f3.1,".dat")') js, foldwidth, wpwidth  
        ! open(unit=10, file="testFoldingData.dat", status="unknown")
        OPEN(unit=k, file=filepath//filename, status='unknown')
        WRITE(k,'(a,2(1pg13.5))') &
        '# test Gaussian folding for w.p.width,fold.width=',wpwidth,foldwidth
        WRITE(k,'(a)') '# r   init w.p.   folded w.p.    exact folding'
        do i = 1, n_grid
            write(k, '(6(1pg13.5))') r(i), psi_in(i), fpsi_out(i),&
            phi_exact(i),fpsi_out(i)/phi_exact(i)
        end do
        close(k)
    end do wpLoop
    end do fwLoop
    end do jsLoop

END PROGRAM test_folding
