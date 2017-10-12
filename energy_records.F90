
! Module that holds the energy records for PME and GB to avoid cyclic
! dependencies.

module energy_records_mod

  implicit none

  type gb_pot_ene_rec
    sequence
    double precision    :: total
    double precision    :: vdw_tot
    double precision    :: elec_tot
    double precision    :: gb
    double precision    :: surf
    double precision    :: bond
    double precision    :: angle
    double precision    :: dihedral
    double precision    :: vdw_14
    double precision    :: elec_14
    double precision    :: restraint
    double precision    :: angle_ub !CHARMM Urey-Bradley Energy
    double precision    :: imp      !CHARMM Improper Energy
    double precision    :: cmap     !CHARMM CMAP
    double precision    :: dvdl
    double precision    :: amd_boost ! AMD boost energy
    double precision    :: emap ! EMAP restraint energy
  end type gb_pot_ene_rec

  integer, parameter    :: gb_pot_ene_rec_size = 17

  type(gb_pot_ene_rec), parameter      :: null_gb_pot_ene_rec = &
    gb_pot_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                   0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)

  ! Potential energies, with breakdown, from pme.  This is intended to be the
  ! external interface to potential energy data produced by this module, in
  ! particular by subroutine pme_force().

  type pme_pot_ene_rec
    sequence
    double precision    :: total
! Added by Lijiang Yang
    double precision    :: total_slu
    double precision    :: total_sol
    double precision    :: total_inter
!------------------------------------------------------
    double precision    :: vdw_tot      ! total of dir, recip
    double precision    :: vdw_dir
! Added by Lijiang Yang
    double precision    :: vdw_dir_slu
    double precision    :: vdw_dir_sol
    double precision    :: vdw_dir_inter
! ---------------------------------------------
    double precision    :: vdw_recip
    double precision    :: elec_tot     ! total of dir, recip, nb_adjust, self
    double precision    :: elec_dir
! Added by Lijiang Yang
    double precision    :: elec_dir_slu
    double precision    :: elec_dir_sol
    double precision    :: elec_dir_inter
! ---------------------------------------------
    double precision    :: elec_recip
    double precision    :: elec_nb_adjust
    double precision    :: elec_self
    double precision    :: hbond
! Added by Lijiang Yang
    double precision    :: hbond_slu
    double precision    :: hbond_sol
    double precision    :: hbond_inter
! ---------------------------------------------
    double precision    :: bond
! Added by Lijiang Yang
    double precision    :: bond_slu
    double precision    :: bond_sol
    double precision    :: bond_inter
! ---------------------------------------------
    double precision    :: angle
! Added by Lijiang Yang
    double precision    :: angle_slu
    double precision    :: angle_sol
    double precision    :: angle_inter
! ---------------------------------------------
    double precision    :: dihedral
! Added by Lijiang Yang
    double precision    :: dihedral_slu
    double precision    :: dihedral_sol
    double precision    :: dihedral_inter
! ---------------------------------------------
    double precision    :: vdw_14
    double precision    :: elec_14
! Added by Lijiang Yang
    double precision    :: vdw_14_slu
    double precision    :: elec_14_slu
    double precision    :: vdw_14_sol
    double precision    :: elec_14_sol
    double precision    :: vdw_14_inter
    double precision    :: elec_14_inter
! ---------------------------------------------
    double precision    :: restraint
    double precision    :: angle_ub !CHARMM Urey-Bradley Energy
    double precision    :: imp      !CHARMM Improper Energy
    double precision    :: cmap     !CHARMM CMAP
    double precision    :: amd_boost ! AMD boost energy
    double precision    :: emap ! EMAP restraint energy
  end type pme_pot_ene_rec

! Modified by Lijiang Yang
  integer, parameter    :: pme_pot_ene_rec_size = 48 

  type(pme_pot_ene_rec), parameter      :: null_pme_pot_ene_rec = &
    pme_pot_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                    0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                    0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                    0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                    0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                    0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
! ---------------------------------------------
                         

  ! Virials, with breakdown, from pme.  This is intended to be the
  ! structure for storing virials; it currently is not exported but could be.

  type pme_virial_rec
    sequence
    double precision    :: molecular(3, 3)
    double precision    :: atomic(3, 3)
    double precision    :: elec_direct(3, 3)
    double precision    :: elec_nb_adjust(3, 3)
    double precision    :: elec_recip(3, 3)
    double precision    :: elec_recip_vdw_corr(3, 3)
    double precision    :: elec_recip_self(3, 3)
    double precision    :: elec_14(3, 3)
    double precision    :: ep_frame(3, 3)       ! for extra points...
    double precision    :: eedvir               ! used in Darden's error est.
  end type pme_virial_rec

  integer, parameter    :: pme_virial_rec_size = 82

  type(pme_virial_rec), parameter      :: null_pme_virial_rec = &
    pme_virial_rec(9*0.d0, 9*0.d0, 9*0.d0, &
                   9*0.d0, 9*0.d0, 9*0.d0, &
                   9*0.d0, 9*0.d0, 9*0.d0, &
                   0.d0)

end module energy_records_mod
