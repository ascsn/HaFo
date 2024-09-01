REAL(8) FUNCTION geff_regpair(efermi_a,ecut_a,pot_a,bmas_a,g_pair)
!
!  Computes the regularization kernel g_eff for
!  regularized pairing.
!  
IMPLICIT NONE
REAL(8),INTENT(IN) :: efermi_a     ! Fermi energy
REAL(8),INTENT(IN) :: ecut_a       ! Cutoff energy
REAL(8),INTENT(IN) :: pot_a        ! local potential U(r)
REAL(8),INTENT(IN) :: bmas_a       ! local inverse effective mass B*(r)=hbar**2/(2M*)
REAL(8),INTENT(IN) :: g_pair       ! local pairing strength from pairing model

  geff_regpair = g_pair/(1-g_pair*Keff_regpair())
!  geff_regpair = g_pair/(1-g_pair*Keff_regpair2())

  RETURN

CONTAINS

REAL(8) FUNCTION Keff_regpair
!
!  Computes the regularization kernel K_eff for
!  regularized pairing.
!  

REAL(8),PARAMETER :: fourpi=12.56637061D0   ! 4*pi
REAL(8),PARAMETER :: fourpi2=39.4784176D0   ! 4*pi**2
REAL(8) :: kf2                     ! square of Fermi momentum
REAL(8) :: kf                      ! absolute value of Fermi momentum
REAL(8) :: kcut                    ! cutoff momentum, deduced from ecut_a
REAL(8) :: ecut_eff                ! cutoff energy relatice to cont. threshold

kf2= (efermi_a-pot_a)/bmas_a
kf = SQRT(ABS(kf2))
ecut_eff = ecut_a  !+efermi_a
kcut = (ecut_eff+efermi_a-pot_a)/bmas_a
IF(kcut<0D0) THEN
  WRITE(6,'(a,3(1pg13.5))') ' e_cut,e_fermi,pot=',ecut_a,efermi_a,pot_a
  STOP "'kcut' is imaginary"
END IF
kcut=SQRT(kcut)

IF(kf2.GE.0D0) THEN
  Keff_regpair = (kcut/(fourpi2*bmas_a)) * &
                 (1D0 - (kf/(kcut+kcut))*LOG((kcut+kf)/(kcut-kf))  )
ELSE
  Keff_regpair = (kcut/(fourpi2*bmas_a)) * &
                 (1D0 + (kf/(kcut))*ATAN((kf/(kcut)))  )
END IF

RETURN

END FUNCTION Keff_regpair


REAL(8) FUNCTION Keff_regpair2
!
!  Computes the regularization kernel K_eff for
!  regularized pairing.
!  

REAL(8),PARAMETER :: fourpi=12.56637061D0   ! 4*pi
REAL(8),PARAMETER :: fourpi2=39.4784176D0   ! 4*pi**2
REAL(8) :: kf2                     ! square of Fermi momentum
REAL(8) :: kf                      ! absolute value of Fermi momentum
REAL(8) :: kcut                    ! cutoff momentum, deduced from ecut_a
REAL(8) :: lcut                    ! lower cutoff momentum, deduced from ecut_a
REAL(8) :: ecut_eff                ! cutoff energy relatice to cont. threshold

REAL(8) :: a,b                     ! auxiliary

ecut_eff = ecut_a  !+efermi_a

IF ( pot_a + ecut_eff - efermi_a < 0D0 ) THEN
   !
   kf = SQRT( ( efermi_a - pot_a ) / bmas_a )
   kcut = SQRT( ( efermi_a + ecut_eff - pot_a ) / bmas_a )
   lcut = SQRT( ( efermi_a - ecut_eff - pot_a ) / bmas_a )
   a = kcut / fourpi2 / bmas_a &
      * ( 1 - kf / 2 / kcut * log( ( kcut + kf ) / ( kcut - kf ) ) )
   b = lcut / fourpi2 / bmas_a &
      * ( 1 - kf / 2 / lcut * log( ( lcut + kf ) / ( kf - lcut ) ) )
!         rega(p,it) = a
!         regb(p,it) = b
   !........................ geff(p,it) = 1 / ( 1 / g - a - b )
ELSE
   IF ( pot_a - efermi_a < 0D0 ) THEN
      !
      kf = SQRT( ( efermi_a - pot_a ) / bmas_a )
      kcut = SQRT( ( efermi_a + ecut_eff - pot_a ) / bmas_a )
      a = kcut / fourpi2 / bmas_a * &
         ( 1 - kf / 2 / kcut * log( ( kcut + kf ) / ( kcut - kf ) ) )
      b = 0D0
!            rega(p,it) = a
!            regb(p,it) = 0D0
      !..................... geff(p,it) = 1 / ( 1 / g - a )
   ELSE
      IF ( pot_a - ecut_eff - efermi_a < 0D0 ) THEN
         !
         kf = SQRT( - ( efermi_a - pot_a ) / bmas_a )
         kcut = SQRT( ( efermi_a + ecut_eff - pot_a ) / bmas_a )
         a = kcut / fourpi2 / bmas_a &
            * ( 1 + kf / kcut * atan( kf / kcut ) )
         b = 0D0
!               rega(p,it) = a
!               regb(p,it) = 0D0
         !.................. geff(p,it) = 1 / ( 1 / g - a )
      ELSE
         !
         a = 0D0
         b = 0D0
!               rega(p,it) = 0D0
!               regb(p,it) = 0D0
         !.................. geff(p,it) = g
      END IF
   END IF
END IF
Keff_regpair2 = a+b

RETURN

END FUNCTION Keff_regpair2

END FUNCTION geff_regpair
