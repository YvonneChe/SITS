#include "copyright.i"

!*******************************************************************************
!
! Module: barostats_mod
!
! Description: Monte-Carlo barostat
!              
!*******************************************************************************

module barostats_mod

  use random_mod, only : random_state, amrset_gen, amrand_gen
#ifdef CUDA
  use pbc_mod
#endif

  implicit none

  private

  double precision, save :: dvmax = 0.02d0

  integer, save :: total_mcbar_attempts = 0
  integer, save :: total_mcbar_successes = 0

  integer, save      :: mcbar_attempts = 0
  integer, save      :: mcbar_successes = 0
  integer, parameter :: dvmax_interval = 10

  type(random_state), save :: mcbar_gen

! Variable explanations
!
! Public variables
!   dvmax           : Size of trial dV move
!   total_mcbar_attempts  : # of trial moves we've attempted with MC barostat
!   total_mcbar_successes : # of trial moves we've accepted with MC barostat
!
!   mcbar_attempts  : # of trial moves; used to adjust size of volume move
!   mcbar_successes : # of successes; used to adjust size of volume move
!   dvmax_interval  : # of exchange attempts to do before deciding whether or
!                     not to change the size of the volume move

  public mcbar_trial, mcbar_setup, mcbar_summary

contains

!*******************************************************************************
!
! Subroutine:  mcbar_setup
!
! Description: Sets up the random number generator for the MC barostat
!              
!*******************************************************************************

#ifdef MPI
! Make sure all threads on the same replica have the exact same random # stream
subroutine mcbar_setup(ig)

  use parallel_dat_mod, only : master, err_code_mpi, pmemd_comm

  implicit none
  include 'mpif.h'
  integer, intent(in) :: ig
  integer :: local_ig

  local_ig = ig + 10

  call mpi_bcast(local_ig, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)

  call amrset_gen(mcbar_gen, local_ig)

  return

end subroutine mcbar_setup
#else
! Non-MPI case below
subroutine mcbar_setup(ig)

  implicit none

  integer, intent(in) :: ig

  call amrset_gen(mcbar_gen, ig + 10)

  return

end subroutine mcbar_setup
#endif

!*******************************************************************************
!
! Subroutine:  mcbar_trial
!
! Description: Perform the trial move for the MC barostat
!              
!*******************************************************************************

#ifdef MPI
subroutine mcbar_trial(atm_cnt, crd, frc, mass, my_atm_lst, new_list, &
                       orig_pot_ene, verbose, pme_err_est, vel, numextra)
#else
subroutine mcbar_trial(atm_cnt, crd, frc, mass, new_list, orig_pot_ene, &
                       verbose, pme_err_est)
#endif

  use cit_mod, only : set_cit_tbl_dims
  use file_io_dat_mod, only : mdout
  use img_mod, only : gbl_img_atm_map, gbl_atm_img_map
  use mdin_ctrl_dat_mod, only : vdw_cutoff, ntp, temp0, pres0, nmropt
  use nmr_calls_mod, only : nmrdcp, skip_print
#ifdef MPI
  use extra_pnts_nb14_mod, only : ep_frames, gbl_frame_cnt, zero_extra_pnts_vec
  use loadbal_mod, only : atm_redist_needed
  use parallel_mod, only : mpi_allgathervec
#endif
  use parallel_dat_mod, only : master
  use pbc_mod, only : uc_volume
  use pme_force_mod, only : pme_pot_ene_rec, pme_force
  use mol_list_mod, only : gbl_mol_cnt

  implicit none

! Passed parameters
  
  integer, intent(in)     :: atm_cnt
  integer, intent(in)     :: verbose
  logical, intent(in out) :: new_list
#ifdef MPI
  integer, intent(in)     :: numextra
  integer, intent(in)     :: my_atm_lst(*)

  double precision, intent(in out)      :: vel(3, atm_cnt)
#endif
  double precision, intent(in out)      :: crd(3, atm_cnt)
  double precision, intent(in out)      :: frc(3, atm_cnt)
  double precision, intent(in)          :: mass(atm_cnt)
  double precision, intent(in out)      :: pme_err_est

  type(pme_pot_ene_rec), intent(in out) :: orig_pot_ene

! Local parameters

  double precision, dimension(3) :: dv
  double precision, dimension(3) :: rmu
  double precision, dimension(3) :: virial
  double precision, dimension(3) :: ekcmt
  double precision, dimension(3, atm_cnt) :: tmp_frc

  logical :: scaled_new_list

  double precision :: randval
  double precision :: pv_work
  double precision :: orig_vol
  double precision :: local_pme_err_est
  double precision :: expfac
  double precision :: nbeta

  type(pme_pot_ene_rec) :: new_pot_ene_rec

! Constants
  double precision, parameter :: ONE_THIRD = 1.d0 / 3.d0
  double precision, parameter :: ONE_HALF = 0.5d0

  ! This is another mcbar attempt
  total_mcbar_attempts = total_mcbar_attempts + 1
  mcbar_attempts = mcbar_attempts + 1

  ! Back up the original volume

  orig_vol = uc_volume

  ! Get the dV move we plan on doing

  if (ntp .eq. 1) then
    ! Isotropic
    call amrand_gen(mcbar_gen, randval)
    dv(1) = (randval - ONE_HALF) * dvmax
    dv(2) = dv(1)
    dv(3) = dv(1)
  else if (ntp .eq. 2) then
    ! Non-isotropic
    call amrand_gen(mcbar_gen, randval)
    dv(1) = (randval - ONE_HALF) * dvmax
    call amrand_gen(mcbar_gen, randval)
    dv(2) = (randval - ONE_HALF) * dvmax
    call amrand_gen(mcbar_gen, randval)
    dv(3) = (randval - ONE_HALF) * dvmax
  else if (ntp .eq. 3) then
    ! Semi-isotropic
    call amrand_gen(mcbar_gen, randval)
    dv(1) = (randval - ONE_HALF) * dvmax
    call amrand_gen(mcbar_gen, randval)
    dv(2) = (randval - ONE_HALF) * dvmax
    dv(3) = dv(2)
  end if

  rmu(1) = (1.d0 + dv(1)) ** ONE_THIRD
  rmu(2) = (1.d0 + dv(2)) ** ONE_THIRD
  rmu(3) = (1.d0 + dv(3)) ** ONE_THIRD

#ifdef MPI
  call scale_system_volume(rmu, verbose, atm_cnt, crd, mass, scaled_new_list, &
                           my_atm_lst)
#else
  call scale_system_volume(rmu, verbose, atm_cnt, crd, mass, scaled_new_list)
#endif

  new_list = new_list .or. scaled_new_list

#ifdef MPI
  ! redistribute coordinates if we're rebuilding the list. (always do this for
  ! now since something appears to be wrong in parallel...)
  ! TODO: Figure out how to avoid global redistributions every time we attempt a
  !       volume change
! if (new_list) then
    call mpi_allgathervec(atm_cnt, crd)
!   if (atm_redist_needed) then
      ! Also need velocities and forces, since we haven't propagated dynamics
      ! yet and the next force call will automagically redistribute atoms...
      call mpi_allgathervec(atm_cnt, vel)
      call mpi_allgathervec(atm_cnt, frc)
      if (numextra .gt. 0) &
        call zero_extra_pnts_vec(vel, ep_frames, gbl_frame_cnt)
!   end if
! end if
#endif

#ifdef CUDA
  call gpu_download_frc(frc)
#endif
  ! Now we need to calculate our next energy, but do not print any nmropt values
  if (nmropt .gt. 0) skip_print = .true.
#ifdef MPI
!Modified by Lijiang Yang
  call pme_force(atm_cnt, crd, tmp_frc, tmp_frc, tmp_frc, tmp_frc, gbl_img_atm_map, gbl_atm_img_map, &
!  call pme_force(atm_cnt, crd, tmp_frc, gbl_img_atm_map, gbl_atm_img_map, &
!                 my_atm_lst, new_list, .true., .false., &
                 my_atm_lst, new_list, .true., .true., .false., &
!------------------------------------------------------------------------
                 new_pot_ene_rec, virial, ekcmt, local_pme_err_est)
#else
!Modified by Lijiang Yang
  call pme_force(atm_cnt, crd, tmp_frc, tmp_frc, tmp_frc, tmp_frc, gbl_img_atm_map, gbl_atm_img_map, &
!  call pme_force(atm_cnt, crd, tmp_frc, gbl_img_atm_map, gbl_atm_img_map, &
!                 new_list, .true., .false., &
                 new_list, .true., .true., .false., &
!------------------------------------------------------------------------
                 new_pot_ene_rec, virial, ekcmt, local_pme_err_est)
#endif /* MPI */

  ! Decrement the NMR opt counter since this force call should not 'count'

  if (nmropt .gt. 0) call nmrdcp

  ! p*dV (6.02204d-2 / 4184.0d0 converts from bar*A^3/particle to kcal/mol)
  pv_work = pres0 * (uc_volume - orig_vol) * 6.02204d-2 / 4184.0d0

  nbeta = -4.184d0 / (8.31441d-3 * temp0)

  expfac = new_pot_ene_rec%total - orig_pot_ene%total + pv_work + &
           gbl_mol_cnt * log(rmu(1) * rmu(2) * rmu(3)) / nbeta

  call amrand_gen(mcbar_gen, randval)

  ! Monte carlo decision

  ! DEBUG
  if (master) &
    write(mdout, '(a)', advance='NO') '| Attempting MC barostat change: '

  if (randval .lt. exp(nbeta * expfac)) then
    ! Accept -- update force array
    frc(1,:) = tmp_frc(1,:)
    frc(2,:) = tmp_frc(2,:)
    frc(3,:) = tmp_frc(3,:)
    pme_err_est = local_pme_err_est
    total_mcbar_successes = total_mcbar_successes + 1
    mcbar_successes = mcbar_successes + 1
    orig_pot_ene = new_pot_ene_rec
    ! DEBUG
    if(master) write(mdout, '(a)') 'Succeeded'
  else
    ! Reject -- rescale everything and put back original forces
    rmu(1) = 1.d0 / rmu(1)
    rmu(2) = 1.d0 / rmu(2)
    rmu(3) = 1.d0 / rmu(3)
#ifdef MPI
    call scale_system_volume(rmu, verbose, atm_cnt, crd, mass, scaled_new_list, &
                             my_atm_lst)
#else
    call scale_system_volume(rmu, verbose, atm_cnt, crd, mass, scaled_new_list)
#endif
    new_list = new_list .or. scaled_new_list
#ifdef CUDA
    ! For GPUs, reload the old forces back to the GPU
    call gpu_upload_frc(frc)
#endif
    ! DEBUG
    if(master) write(mdout, '(a)') 'Failed'
  end if

  if (mod(mcbar_attempts, dvmax_interval) .eq. 0) then

    ! If our success fraction is too large or too small, adjust dvmax

    if (mcbar_successes .ge. 0.75d0 * mcbar_attempts) then

      dvmax = dvmax * 1.1d0
      mcbar_attempts = 0
      mcbar_successes = 0
      if (master) &
        write(mdout, '(a)') '| MC Barostat: Increasing size of volume moves'

    else if (mcbar_successes .le. 0.25d0 * mcbar_attempts) then

      dvmax = dvmax / 1.1d0
      mcbar_attempts = 0
      mcbar_successes = 0
      if (master) &
        write(mdout, '(a)') '| MC Barostat: Decreasing size of volume moves'

    end if

  end if
  
end subroutine mcbar_trial

!*******************************************************************************
!
! Subroutine:  scale_system_volume
!
! Description: Scales the system volume and checks if a new pairlist is needed
!              
!*******************************************************************************

#ifdef MPI
subroutine scale_system_volume(rmu, verbose, atm_cnt, crd, mass, new_list, &
                               my_atm_lst)
#else
subroutine scale_system_volume(rmu, verbose, atm_cnt, crd, mass, new_list)
#endif
  
  use cit_mod, only : set_cit_tbl_dims
  use constraints_mod, only : natc, atm_xc
  use dynamics_dat_mod, only : gbl_mol_mass_inv, gbl_my_mol_lst, gbl_mol_com
  use mdin_ctrl_dat_mod, only : vdw_cutoff, ntp, ntr
  use mdin_ewald_dat_mod, only : skinnb
#ifdef MPI
  use nb_pairlist_mod, only : gbl_atm_saved_crd, check_my_atom_movement
#else
  use nb_pairlist_mod, only : gbl_atm_saved_crd, check_all_atom_movement
#endif
  use parallel_dat_mod ! for MPI functions and MPI-related variables
  use pbc_mod, only : pressure_scale_pbc_data, pressure_scale_crds, &
                      pressure_scale_restraint_crds, pbc_box, cut_factor
  implicit none

! Formal arguments
  
  integer, intent(in)  :: atm_cnt
  integer, intent(in)  :: verbose
#ifdef MPI
  integer, intent(in)  :: my_atm_lst(*)
#endif
  logical, intent(out) :: new_list

  double precision, intent(in out) :: crd(3, atm_cnt)
  double precision, intent(in)     :: mass(atm_cnt)
  double precision, intent(in)     :: rmu(3)

  integer, dimension(2) :: new_list_buf

  new_list = .false.
  new_list_buf(:) = 0

! Local variables
  
  call pressure_scale_pbc_data(rmu, vdw_cutoff + skinnb, verbose)
#ifdef CUDA
  call gpu_pressure_scale(ucell, recip, uc_volume)
#else
  call set_cit_tbl_dims(pbc_box, vdw_cutoff + skinnb, cut_factor)
#endif

#ifdef MPI
  call pressure_scale_crds(crd, mass, gbl_mol_mass_inv, gbl_my_mol_lst, &
                           gbl_mol_com)
#else
  call pressure_scale_crds(crd, mass, gbl_mol_mass_inv, gbl_mol_com)
#endif

  ! If we have restraints, we need to scale restraint coordinates
  if (ntr .gt. 0 .and. natc .gt. 0) then
#ifdef MPI
    call pressure_scale_restraint_crds(atm_xc, mass, gbl_mol_mass_inv, &
                                       gbl_my_mol_lst)
#else
    call pressure_scale_restraint_crds(atm_xc, mass, gbl_mol_mass_inv)
#endif
  end if

  ! Now we have to see if we need to rebuild our pairlist

#ifdef CUDA
  call gpu_skin_test()
#else

#ifdef MPI
  call check_my_atom_movement(crd, gbl_atm_saved_crd, my_atm_lst, &
                              skinnb, ntp, new_list)
  if (new_list) &
    new_list_buf(1) = 1
  
  ! Make sure everyone knows if we need to rebuild a list according to any
  ! worker thread
  call mpi_allreduce(new_list_buf(1), new_list_buf(2), 1, mpi_integer, &
                     mpi_sum, pmemd_comm, err_code_mpi)
  
  new_list = new_list_buf(2) .gt. 0
#else
  call check_all_atom_movement(atm_cnt, crd, gbl_atm_saved_crd, skinnb, ntp, &
                               new_list)
#endif /* MPI */

#endif /* CUDA */

end subroutine scale_system_volume

!*******************************************************************************
!
! Subroutine: mcbar_summary
!
! Description: Print out a summary of the MC barostat statistics
!
!*******************************************************************************

subroutine mcbar_summary()

  use parallel_dat_mod, only : master

  implicit none

  double precision :: success_percent

  if (.not. master) return

  success_percent = &
     dble(total_mcbar_successes) / dble(total_mcbar_attempts) * 100.d0

  write(6, '(("| MC Barostat: ",i10," volume changes attempted."))') &
     total_mcbar_attempts
  
  write(6, '(("| MC Barostat: ",i10," changes successful (",1f6.2,"%)"))') &
     total_mcbar_successes, success_percent
  write(6, '(t2,78("-"),/)')

  return

end subroutine mcbar_summary

end module barostats_mod
