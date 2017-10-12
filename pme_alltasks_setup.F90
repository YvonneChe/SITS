#include "copyright.i"

!*******************************************************************************
!
! Module: pme_alltasks_setup_mod
!
! Description:  Setup of data structures for uniprocessor code as well as
!               mpi master and slave processors.
!              
!*******************************************************************************

module pme_alltasks_setup_mod

use file_io_dat_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  pme_alltasks_setup
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine pme_alltasks_setup(num_ints, num_reals)

  use angles_mod
  use angles_ub_mod
  use bonds_mod
  use cit_mod
  use constraints_mod
  use dihedrals_mod
  use dihedrals_imp_mod
  use cmap_mod
  use dynamics_mod
  use dynamics_dat_mod
  use extra_pnts_nb14_mod
  use ene_frc_splines_mod
  use loadbal_mod
  use pme_direct_mod
  use pme_force_mod
  use pme_recip_dat_mod
  use pme_blk_recip_mod
  use pme_slab_recip_mod
  use mol_list_mod
  use prfs_mod
  use img_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use nb_exclusions_mod
  use nb_pairlist_mod
  use parallel_dat_mod
  use parallel_mod
  use pbc_mod
  use prfs_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use random_mod
  use charmm_mod, only : charmm_active
  use mdin_debugf_dat_mod
#ifdef CUDA
  use nmr_calls_mod
#endif

  implicit none

! Formal arguments:

! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed

#ifdef MPI
#else
  integer               :: use_atm_map(natom)

  use_atm_map(:) = 0
#endif


#ifdef MPI
! Send and receive common blocks from the master node.  In the non-masters,
! memory is allocated as needed.

  call bcast_mdin_ewald_dat
  call bcast_img_dat(natom)                              
  call bcast_nb_pairlist_dat(natom, vdw_cutoff + skinnb)
  call bcast_nb_exclusions_dat(natom, next)
  call bcast_pme_force_dat(ntypes)
  call bcast_pme_direct_dat
  call bcast_ene_frc_splines_dat

  ! The mask lists made here are used in nonbonded pairlist processing
  ! to keep nonbonded calcs from being done on excluded atoms.  The lists are
  ! static (atom-based).

  call make_atm_excl_mask_list(natom, atm_numex, gbl_natex)

  if (ntb .ne. 0) then
    call bcast_pbc
    if (.not. master) then
      call set_cit_tbl_dims(pbc_box, vdw_cutoff + skinnb, cut_factor)
    end if
  end if

  ! Need new molecule list structures to support no_intermolecular_bonds option.

  if ((ntp .gt. 0 .and. imin .eq. 0) .or. &
      (emil_do_calc .gt. 0)          .or. &
      (master .and. iwrap .ne. 0)) then
    call setup_molecule_lists(num_ints, num_reals, natom, nspm, atm_nsp, &
                              nbona, gbl_bond(bonda_idx), &
                              no_intermolecular_bonds)
  end if

  ! Old molecule lists no longer valid...
  if  ( emil_do_calc .eq. 0 ) then !...but still useful for emil calc.
     num_ints = num_ints - size(atm_nsp)
     deallocate(atm_nsp)
     nspm = 0
  end if

  call dynamics_dat_init(natom, ntp, imin, atm_mass, atm_crd, &
                         num_ints, num_reals)

  call setup_prfs(num_ints, num_reals, natom, nres, gbl_res_atms, &
                  nbonh, gbl_bond, gbl_frame_cnt, ep_frames)

  ! The frc-ene task data structures must be set up before we can do reciprocal
  ! space data structures set up.  This is a facility that allows for assigning
  ! no frc-ene workload to some processes; it is not used in this implementation
  ! and will likely be superceded by a facility allowing a partial frc-ene
  ! workload on some tasks to enable them to arrive at communications
  ! rendevous first.

  call setup_frc_ene_task_dat(num_ints, num_reals)

  ! We must know about fft slab or block allocations before we do atom division
  ! in parallel_setup():

  call pme_recip_dat_setup(num_ints, num_reals)

  if (block_fft .eq. 0) then
    call pme_slab_recip_setup(num_ints, num_reals)
  else
    call pme_blk_recip_setup(frc_ene_task_cnt, gbl_frc_ene_task_lst, &
                             num_ints, num_reals)
  end if

  ! Set up the pseudo-residue fragment data now; it will also be needed when
  ! atom division is done in parallel_setup().

  if (ntp .gt. 0 .and. imin .eq. 0) then
    call setup_mol_prf_data(num_ints, num_reals, natom, gbl_prf_cnt, &
                            gbl_prf_listdata, gbl_prf_lists)

    call setup_fragmented_molecules(natom, num_ints, num_reals)
  end if

  ! Divide atoms up among the processors.  The atom division is redone
  ! periodically under cit, and is either residue or molecule-based, with
  ! locality.  In other words, under cit a contiguous block of atoms owned by
  ! each process is a thing of the past.  This is done along with various
  ! allocations specific to the pme parallel implementation in parallel_setup():

  call parallel_setup(num_ints, num_reals)

  ! Master needs space to receive data for load balancing.  The 5 values
  ! for each task are the "direct force time" due to image nonbonded calcs,
  ! the "reciprocal force time" due to pme nonbonded reciprocal force calcs,
  ! the "atom owner time" due to computations required of atom owners, and the
  ! "atom user time" due to managing used atom and image data structures

  if (master) then
    allocate(gbl_loadbal_node_dat(4, 0:numtasks - 1), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
    num_ints = num_ints + size(gbl_loadbal_node_dat)
    gbl_loadbal_node_dat(:,:) = 0
  end if

  call bcast_mdin_debugf_dat()

#else /* begin non-MPI code */

  ! The mask lists made here are used in nonbonded pairlist processing
  ! to keep nonbonded calcs from being done on excluded atoms.  The lists are
  ! static (atom-based).

  call make_atm_excl_mask_list(natom, atm_numex, gbl_natex)

  ! Need new molecule list structures to support no_intermolecular_bonds option.
  if ( (ntp .gt. 0 .and. imin .eq. 0) .or. &
                       (iwrap .ne. 0) .or. &
                (emil_do_calc .gt. 0) ) then
    call setup_molecule_lists(num_ints, num_reals, natom, nspm, atm_nsp, &
                              nbona, gbl_bond(bonda_idx), &
                              no_intermolecular_bonds)
  end if

  ! Old molecule lists no longer valid...
  if  ( emil_do_calc .eq. 0 ) then !...but still useful for emil calc.
     num_ints = num_ints - size(atm_nsp)
     deallocate(atm_nsp)
     nspm = 0
  end if

  ! Initialize dynamics data.

  call dynamics_dat_init(natom, ntp, imin, atm_mass, atm_crd, &
                         num_ints, num_reals)

  ! Uniprocessor values:

  my_mol_cnt = gbl_mol_cnt

  ! We set up reciprocal force data structures here in parallel with
  ! where it has to be done for mpi code:

  call pme_recip_dat_setup(num_ints, num_reals)
  call pme_slab_recip_setup(num_ints, num_reals)

  call bonds_setup(num_ints, num_reals, use_atm_map)
  call angles_setup(num_ints, num_reals, use_atm_map)
  call dihedrals_setup(num_ints, num_reals, use_atm_map)
#ifdef CUDA 
  if (nmropt .ne. 0) then
    call cuda_nmr_setup()         
  end if
#endif
  call nb14_setup(num_ints, num_reals, use_atm_map)

  if (charmm_active) then
    call angles_ub_setup(num_ints, num_reals, use_atm_map)
    call dihedrals_imp_setup(num_ints, num_reals, use_atm_map)
    call cmap_setup(num_ints, num_reals, use_atm_map)
  endif

  ! gbl_bond is still needed for shake setup and resetup, but gbl_angle and
  ! gbl_dihed can be deallocated .

  !Added by Lijiang Yang
  num_ints = num_ints - size(gbl_bondslu) * bond_rec_ints
  num_ints = num_ints - size(gbl_bondsol) * bond_rec_ints
  num_ints = num_ints - size(gbl_bondinter) * bond_rec_ints
  num_ints = num_ints - size(gbl_angle) * angle_rec_ints
  num_ints = num_ints - size(gbl_angleslu) * angle_rec_ints
  num_ints = num_ints - size(gbl_anglesol) * angle_rec_ints
  num_ints = num_ints - size(gbl_angleinter) * angle_rec_ints
  num_ints = num_ints - size(gbl_dihed) * dihed_rec_ints
  num_ints = num_ints - size(gbl_dihedslu) * dihed_rec_ints
  num_ints = num_ints - size(gbl_dihedsol) * dihed_rec_ints
  num_ints = num_ints - size(gbl_dihedinter) * dihed_rec_ints
  !---------------------------------------------------------
!  num_ints = num_ints - size(gbl_dihed) * dihed_rec_ints
  if (charmm_active) then
    num_ints = num_ints - size(gbl_angle_ub) * angle_ub_rec_ints
    num_ints = num_ints - size(gbl_dihed_imp) * dihed_imp_rec_ints
    num_ints = num_ints - size(gbl_cmap) * cmap_rec_ints

    deallocate(gbl_angle_ub, gbl_dihed_imp, gbl_cmap)
  endif

  deallocate(gbl_bondslu, gbl_bondsol, gbl_bondinter)
  deallocate(gbl_angle, gbl_angleslu, gbl_anglesol, gbl_angleinter)
  deallocate(gbl_dihed, gbl_dihedslu, gbl_dihedsol, gbl_dihedinter)
!  deallocate(gbl_angle, gbl_dihed)

  call make_nb_adjust_pairlst(natom, use_atm_map, &
                              gbl_nb_adjust_pairlst, &
                              atm_nb_maskdata, atm_nb_mask)

#endif /* not MPI */

! Initialize random number generator at same point in all processors. Then, if
! random initial velocities are needed, generate them in all processors. 
! In general, we must be careful to generate the same sequence of random
! numbers in all processors.

  call amrset(ig)

  if (ntx .eq. 1 .or. ntx .eq. 2) then
    call all_atom_setvel(natom, atm_vel, atm_mass_inv, tempi)
#ifdef CUDA
    call gpu_upload_vel(atm_vel)
#endif
  end if

  if (ibelly .gt. 0) then
    call all_atom_belly(natom, atm_igroup, atm_vel)
  end if
  
#ifdef CUDA
  call gpu_pme_alltasks_setup(nfft1, nfft2, nfft3, gbl_prefac1, gbl_prefac2, gbl_prefac3, ew_coeff, ips)
#endif

  return

end subroutine pme_alltasks_setup

end module pme_alltasks_setup_mod
