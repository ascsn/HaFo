SUBROUTINE params_from_NMP(power_from_K,B0_deduced,B3_deduced,&
                           B1_deduced,&
                           rhotot,eovera,hbar2m,effmas,compr)
IMPLICIT NONE
!
! computes isoscalar Skyrme parameters from NMP.
! Input parameters (NMP) are:
!  rhotot = equilibrium density of symmetric matter,
!  eovera = binding energy per nucleon at equilibrium,
!  hbar2m = hbar**2/(2*nucleon_mass),
!  effmas = effective mass m*/m
!  compr  = incompressibility at equilibrium.
! Output parameters are:
!  power_from_K = power of density dependence,
!  B0_deduced   = B0 parameter,
!  B1_deduced   = B1 parameter,
!  B3_deduced   = B3 parameter.
!
!  Note: The B3 parameter deduced here is factor 6 smaller than
!        the equivalent used in the 1D code.
!
  REAL(8), INTENT(IN) :: rhotot,eovera,hbar2m,effmas,compr
  REAL(8), INTENT(OUT) :: power_from_K,B0_deduced,B3_deduced,B1_deduced

! (3/5)*(3*pi^2)^*(2/3)  factor for tau_n/p
  REAL(8), PARAMETER :: f2d3= 2D0/3D0, f5d3=5D0/3D0 
  REAL(8), PARAMETER :: betab=5.742467999D0 
  REAL(8), PARAMETER :: betat=betab*0.5D0**f2d3
  REAL(8) :: counter,denom,b1,power
  REAL(8) :: hbmterm      ! kinetic energy term
  REAL(8) :: b1term       ! kinetic interaction term
  REAL(8) :: rhopowp1     ! rho**(1+power)/2

  b1       = hbar2m*(1D0/effmas-1D0)/rhotot
  B1_deduced = b1

  hbmterm  = hbar2m*betat*rhotot**f2d3
  b1term   = b1*betat*rhotot**f5d3

  counter = 10D0*b1term-2D0*hbmterm-compr
  denom   = 9D0*eovera+6D0*b1term-3D0*hbmterm
  power   = counter/denom-1D0
  power_from_K = power


  rhopowp1 = rhotot**(1D0+power)/2D0
  B3_deduced = (hbmterm/3D0-f2d3*b1term-eovera)/(power*rhopowp1)
  B0_deduced = 2D0*(eovera-hbmterm-b1term-B3_deduced*rhopowp1)/rhotot

  RETURN
END SUBROUTINE params_from_NMP
