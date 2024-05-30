!     common block to supply the number of observable for which
!  'extobs' is to be run actually.
      parameter (numextobs=44)
      common /extobsc/nextobs(numextobs)

!     variables in 'extobs':
      equivalence (nextobs(1),igmrpb)      ! GMR in 208Pb
      equivalence (nextobs(2),igdrpb)      ! GDR in 208Pb
      equivalence (nextobs(3),igqrpb)      ! GQR in 208Pb
      equivalence (nextobs(4),igdro)       ! GDR in 16O
      equivalence (nextobs(5),idipol)      ! dipole polarizability in 208Pb
      equivalence (nextobs(6),idipsn)      ! dipole polarizability in 132Sn
      equivalence (nextobs(7),idipca)      ! dipole polarizability in 48Ca
      equivalence (nextobs(8),idipni)      ! dipole polarizability in 68Ni
      equivalence (nextobs(9),ipygmy)      ! pygmy in 208Pb 
      equivalence (nextobs(10),ipygsn)     ! pygmy in 132Sn 
      equivalence (nextobs(11),imonop)     ! monopole polarizability in 208Pb
      equivalence (nextobs(12),imonsn)     ! monopole polarizability in 132Sn
      equivalence (nextobs(13),imonni)     ! monopole polarizability in 68Ni
      equivalence (nextobs(14),imonca)     ! monopole polarizability in 488Ca
      equivalence (nextobs(15),iquadp)     ! quadrupole polarizability in 208Pb
      equivalence (nextobs(16),idipzr)     ! dipole polarizability in 90Zr

!     variables in 'extobsmf':
      equivalence (nextobs(21),iskinpb)    ! neutron rsm skin in 208Pb
      equivalence (nextobs(22),ilevelsn)   ! neutron level sequence in 132Sn
      equivalence (nextobs(23),ishe)       ! E_B  for 120/182
      equivalence (nextobs(24),iqalpha )   ! Q_alpha for 120/182
      equivalence (nextobs(25),is2n)       ! S_2n for 120/182
      equivalence (nextobs(26),is2p)       ! S_2p for 120/182
      equivalence (nextobs(27),iisotop)    ! isotope shift 214Pb-208Pb
      equivalence (nextobs(28),iskinsn74)  ! neutron rms skin in 124Sn
      equivalence (nextobs(29),iskinsn82)  ! neutron rms skin in 132Sn
      equivalence (nextobs(30),iskinsn90)  ! neutron rms skin in 140Sn
      equivalence (nextobs(31),iskinca)    ! neutron rms skin in 48Ca
      equivalence (nextobs(32),iskinni)    ! neutron rms skin in 68Ni
      equivalence (nextobs(33),iskindf)    ! neutron diffraction skin in 208Pb
      equivalence (nextobs(34),irmsneut)   ! neutron rms radius in 208Pb
      equivalence (nextobs(35),irdfneut)   ! neutron diffr. radius in 208Pb
      equivalence (nextobs(36),ifrmneut)   ! weak formfactor in 208Pb
      equivalence (nextobs(37),ibindexo)   ! binding energy of 148Sn
      equivalence (nextobs(38),irdfcoul)   ! R_diffr_charge 208Pb
      equivalence (nextobs(39),isigcoul)   ! sigma_charge 208Pb
      equivalence (nextobs(40),ifrmne48)   ! weak charge FF for 48Ca
      equivalence (nextobs(41),isurfe)     ! surface energy from SC
      equivalence (nextobs(42),isurfas)    ! surface asymmetry from SC
      equivalence (nextobs(43),irmscoul)   ! rms_charge 208Pb
      equivalence (nextobs(44),ifrmapv)    ! weak A_PV in 208Pb


