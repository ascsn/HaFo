PROGRAM test_folding
    ! Test program for Gaussian folding subroutine.
!    use spherical_harmonics
    IMPLICIT NONE

!   common 'grid' to communicate 'ngrid'
    include 'grid.inc'

!   set the model parameters here by PARAMETER statements
!   (may be changed to be read in later)
    INTEGER, PARAMETER :: n_grid = 50      ! number of radial grid points
    ! INTEGER, PARAMETER :: js=2             ! angular momentum
    INTEGER, PARAMETER :: jsmax=15          ! maximum angular momentum
    INTEGER :: js
    REAL(8), PARAMETER :: rstep_in=0.3D0   ! grid spacing
    ! REAL(8), PARAMETER :: foldwidth=3.0D0    ! Gaussian folding width
    ! REAL(8), PARAMETER :: wpwidth=3D0      ! with of test wavepacket
    REAL(8) :: foldwidth, wpwidth
    REAL(8) :: psi_norm,psi_width,fpsi_norm,fpsi_width,exact_norm,exact_width

    INTEGER :: i, ll, k
    INTEGER :: fw, wp
    REAL(8), parameter :: pi = 3.141592653589793d0
    REAL(8) :: t1, t2
    REAL(8), DIMENSION(n_grid) :: psi_in, r, deltar, deltar2, fpsi_out, phi_exact
    REAL(8), DIMENSION(n_grid) :: redPsi_in
    COMPLEX(8) :: Yll
    CHARACTER(len=100) :: filename
    CHARACTER(*), parameter :: filepath="./test_foldwidth_15fm/"
!    CHARACTER(*), parameter :: filepath="./results_data2/"
    ngrid = n_grid            ! to be available in 'k2psi'
    rstep = rstep_in          ! to be available in 'k2psi'

    ! Define a loop over js, foldwidth, wpwidth
    jsLoop: do js = 1, jsmax
    fwLoop: do fw = 1, 20
        foldwidth = 0.3d0 + dble(fw - 1)*0.3d0
    wpLoop: do wp = 1, 8
        wpwidth = 1.0d0 + dble(wp-1)*0.5d0
        ll = js/2  ! orbital angular momentum to combined a.m. quantum number 'js'

        ! Define the wavefunction to be folded
        do i = 1, n_grid
            r(i) = (i-1)*rstep
        end do

        do i = 1, n_grid
            psi_in(i) = r(i)**(ll+1)*EXP(-(r(i)/wpwidth)**2)
        end do
        psi_in=psi_in/(pi**1.5D0*wpwidth**3)*4D0*pi
        
        ! SUBROUTINE Gauss_folding_m(psi_in,fpsi_out,rstep,ngrid,js,foldwid,jsmax)
        call Gauss_folding_m(psi_in, fpsi_out, rstep, n_grid, js, foldwidth, jsmax)

        ! Test the norm, widths of the new function wrt old function
        do i = 1, n_grid
            t1 = 1D0/(wpwidth**2+foldwidth**2) 
            t2  = wpwidth**2*t1
            phi_exact(i) = ((-1)**(2*ll))*t2**(ll+1.5D0)*r(i)**(ll+1)*EXP(-r(i)**2*t1)
        end do
        phi_exact=phi_exact/(pi**1.5D0*wpwidth**3)*4D0*pi

        ! print*, "phi_exact", phi_exact
        ! print*, "fpsi", fpsi_out
        ! print*, phi_exact - fpsi_out

        ! Calculate norms



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
        psi_norm = rstep*SUM(r*psi_in)
        fpsi_norm = rstep*SUM(r*fpsi_out)
        exact_norm = rstep*SUM(r*phi_exact)
        psi_width = SQRT(rstep*SUM(r**3*psi_in)/psi_norm)
        fpsi_width = SQRT(rstep*SUM(r**3*fpsi_out)/fpsi_norm)
        exact_width = SQRT(rstep*SUM(r**3*phi_exact)/exact_norm)
!        write(k,'(a,5(1pg13.5))') '# norms initial w.p.,folded,exact,ratio=', &
!          rstep*SUM(r*psi_in),rstep*SUM(r*fpsi_out),rstep*SUM(r*phi_exact), &
!          SUM(r*fpsi_out)/SUM(r*psi_in)
        write(k,'(a,5(1pg13.5))') '# norms initial,folded,exact,ratio=', &
          psi_norm,fpsi_norm,exact_norm,fpsi_norm/MAX(psi_norm,1D-90)
        write(k,'(a,5(1pg13.5))') '#widths initial,folded,exact,ratio=', &
          psi_width,fpsi_width,exact_width,fpsi_width/MAX(exact_width,1D-90)
        IF(js==1) write(k,'(a,5(1pg13.5))') '#widths analytical                =', &
          SQRT(1.5D0)*wpwidth,SQRT(1.5D0*(wpwidth**2+foldwidth**2))
        close(k)
    end do wpLoop
    end do fwLoop
    end do jsLoop

END PROGRAM test_folding
