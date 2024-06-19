module spherical_harmonics
    implicit none
    integer :: l, m
    real(8) :: theta, phi
    complex(8) :: Ylm
    complex(8), parameter :: pi = 3.141592653589793d0

contains

    ! Function to calculate the spherical harmonics
    complex(8) function spherical_harmonic(l, m, theta, phi)
        integer, intent(in) :: l, m
        real(8), intent(in) :: theta, phi
        real(8) :: norm_factor, P_lm
        complex(8) :: e_imphi

        ! Calculate the normalization factor
        norm_factor = sqrt((2.0d0*l+1.0d0)/(4.0d0*pi) * factorial(l-m)/factorial(l+m))

        ! Calculate the associated Legendre polynomial
        P_lm = legendre_polynomial(l, m, cos(theta))

        ! Calculate the exponential factor
        e_imphi = 1.d0 !cmplx(cos(m*phi), sin(m*phi))

        ! Combine all parts to get the spherical harmonic
        spherical_harmonic = norm_factor * P_lm * e_imphi
    end function spherical_harmonic

    ! Function to calculate the associated Legendre polynomials
    real(8) function legendre_polynomial(l, m, x)
        integer, intent(in) :: l, m
        real(8), intent(in) :: x
        integer :: k
        real(8) :: pmm, pmmp1, pll

        if (m < 0) then
            print *, "Error: m must be non-negative."
            stop
        end if

        pmm = 1.0d0
        if (m > 0) then
            pmmp1 = sqrt(1.0d0 - x*x)
            pmm = pmmp1
            do k = 2, m
                pmmp1 = sqrt(1.0d0 - x*x)
                pmm = pmm * pmmp1
            end do
        end if

        if (l == m) then
            legendre_polynomial = pmm
            return
        end if

        pmmp1 = x * (2*m+1) * pmm
        if (l == m+1) then
            legendre_polynomial = pmmp1
            return
        end if

        do k = m+2, l
            pll = ((2*k-1)*x*pmmp1 - (k+m-1)*pmm) / (k-m)
            pmm = pmmp1
            pmmp1 = pll
        end do

        legendre_polynomial = pll
    end function legendre_polynomial

    ! Function to calculate the factorial of an integer
    real(8) function factorial(n)
        integer, intent(in) :: n
        integer :: i
        factorial = 1.0d0
        do i = 2, n
            factorial = factorial * i
        end do
    end function factorial

end module spherical_harmonics
