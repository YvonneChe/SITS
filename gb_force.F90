#include "copyright.i"

!*******************************************************************************
!
! Module: gb_force_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module gb_force_mod

  use energy_records_mod, only : gb_pot_ene_rec, gb_pot_ene_rec_size, &
                                 null_gb_pot_ene_rec
  use gbl_datatypes_mod

  implicit none

  ! Potential energies, with breakdown, from GB.  This is intended to be the
  ! external interface to potential energy data produced by this module, in
  ! particular by subroutine gb_force().

contains

!*******************************************************************************
!
! Subroutine:  gb_force
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine gb_force(atm_cnt, crd, frc, pot_ene, irespa, need_pot_ene)

  use angles_mod
  use angles_ub_mod
  use bonds_mod
  use charmm_mod, only : charmm_active
  use constantph_dat_mod, only : on_cpstep, proposed_qterm, cnstph_frc
  use constraints_mod
  use cmap_mod
  use dihedrals_mod
  use dihedrals_imp_mod
  use emap_mod,only:emapforce
  use extra_pnts_nb14_mod
  use gb_ene_mod
  use gb_parallel_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_debugf_dat_mod, only : do_charmm_dump_gold
  use nmr_calls_mod
  use parallel_dat_mod
  use parallel_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use runfiles_mod
  use timers_mod
  use charmm_mod, only : charmm_active
  use mdin_debugf_dat_mod, only : do_charmm_dump_gold
  use amd_mod
  use scaledMD_mod
  use gbsa_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt
  double precision              :: crd(3, atm_cnt)
  double precision              :: frc(3, atm_cnt)
  type(gb_pot_ene_rec)          :: pot_ene
  type(gb_pot_ene_rec)          :: pot_ene_prop ! CpH proposed state energy rec.
  integer, intent(in)           :: irespa
  logical, intent(in)           :: need_pot_ene

! Local variables:

  double precision              :: enmr(3)
  double precision              :: enemap
  integer                       :: i, j
  double precision              :: temp_tot_dih
  double precision              :: totdih
  integer                       :: buf_size

  ! Constant pH energy terms for proposed state

  double precision              :: cph_egb
  double precision              :: cph_eel
  double precision              :: cph_eel14
  double precision              :: cph_vdw
  double precision              :: cph_vdw14
  double precision              :: cph_surf

  ! End constant pH energy terms for proposed state


#ifdef CUDA
  ! Update any weight modifications and zero-out the potential energy record
  if (nmropt .ne. 0) then
    call nmr_weight(atm_cnt, crd, 6)
    pot_ene = null_gb_pot_ene_rec !Zero's entire structure
    enmr(:) = 0
  end if
  if (need_pot_ene) then !DO NOT CHANGE THIS OR ADD ADDITIONAL CLAUSES, THIS IS SET ONCE IN RUNMD.F90.
#ifdef MPI    
    if((iamd.eq.2).or.(iamd.eq.3))then
      totdih = 0.0
      if(num_amd_lag .eq. 0)then
        temp_tot_dih = 0.0
        !Note this double calculates the dihedral energy since we do it
        !here for AMD and then repeat it later - consider merging these
        !at some point for efficiency.
        call gpu_calculate_gb_amd_dihedral_energy(totdih)
        buf_size = 1
        call mpi_allreduce(totdih, temp_tot_dih, &
                       buf_size, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)
        totdih = temp_tot_dih
      end if
      call gpu_calculate_amd_dihedral_weight(totdih)
    end if
#else    
    if((iamd.eq.2).or.(iamd.eq.3))then
      call gpu_calculate_gb_amd_dihedral_energy_weight()
    endif
#endif    
! BEGIN GBSA implementation to work with the GPU code.
! It runs on the CPU for now, needs to write a kernel for this.
! Regular CPU call is located in gb_ene
! pmemd.cuda only takes nrespa=1 for now, so GBSA contribution must ALWAYS
! be calculated (regardless of need_pot_ene)
    if (gbsa .eq. 1) then
      frc=0.d0
      gbsafrc=0.d0
      call gpu_download_crd(crd)
      call gbsa_ene(crd, gbsafrc, pot_ene%surf ,atm_cnt, jj, r2x, belly_atm_cnt) 
    end if
! END GBSA implementation to work with the GPU code.
    if (on_cpstep) then
      call gpu_refresh_charges(cit_nb14, gbl_one_scee, proposed_qterm)
      call gpu_gb_ene(pot_ene_prop, enmr)
      call gpu_refresh_charges(cit_nb14, gbl_one_scee, atm_qterm)
      call gpu_gb_ene(pot_ene, enmr)
      pot_ene%dvdl = pot_ene_prop%total - pot_ene%total
    else
      call gpu_gb_ene(pot_ene, enmr)
    end if
! BEGIN GBSA implementation to work with the GPU code.
! force elements sent to the GPU to be added to the other force contributions
    if (gbsa .eq. 1) then
#ifdef MPI    
      buf_size = atm_cnt *3
        call mpi_allreduce(gbsafrc, frc, &
                       buf_size, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)
      call gpu_gbsa_frc_add(frc)
#else
      call gpu_gbsa_frc_add(gbsafrc)
#endif    
      pot_ene%total =  pot_ene%total +  pot_ene%surf
    end if
! END GBSA implementation to work with the GPU code.
#ifdef MPI    
    call gb_distribute_enes(pot_ene)
#endif    
  else
! BEGIN GBSA implementation to work with the GPU code.
! It runs on the CPU for now, needs to write a kernel for this.
! Regular CPU call is located in gb_ene
! pmemd.cuda only takes nrespa=1 for now, so GBSA contribution must ALWAYS
! be calculated (regardless of need_pot_ene)
    if (gbsa .eq. 1) then
      frc=0.d0
      gbsafrc=0.d0
      call gpu_download_crd(crd)
      call gbsa_ene(crd, gbsafrc, pot_ene%surf, atm_cnt, jj, r2x, belly_atm_cnt) 
    end if
! END GBSA implementation to work with the GPU code.
    call gpu_gb_forces()
! BEGIN GBSA implementation to work with the GPU code.
! force elements sent to the GPU to be added to the other force contributions
    if (gbsa .eq. 1) then
#ifdef MPI    
      buf_size = atm_cnt *3
        call mpi_allreduce(gbsafrc, frc, &
                       buf_size, mpi_double_precision, &
                       mpi_sum, pmemd_comm, err_code_mpi)
     call gpu_gbsa_frc_add(frc)
#else
     call gpu_gbsa_frc_add(gbsafrc)
#endif    
    end if
! END GBSA implementation to work with the GPU code.
  end if
  if (nmropt .ne. 0) then
    call nmr_calc(crd, frc, enmr, 6)
  end if
  if(iamd.gt.0)then
!AMD calculate weight and scale forces
    call gpu_calculate_and_apply_amd_weights(pot_ene%total, pot_ene%dihedral, &
                                             pot_ene%amd_boost, num_amd_lag)
     ! Update total energy
     pot_ene%total = pot_ene%total + pot_ene%amd_boost
  endif

  if(scaledMD.gt.0)then
!scaledMD scale forces
    call gpu_scaledmd_scale_frc(pot_ene%total)
    pot_ene%total = pot_ene%total * scaledMD_lambda
  endif

  call update_time(nonbond_time)
  return
#else /* CUDA */

! Zero energies that are stack or call parameters:

  pot_ene = null_gb_pot_ene_rec
  pot_ene_prop = null_gb_pot_ene_rec

! Zero internal energies:

  enmr(:) = 0.d0

  if (ti_mode .ne. 0) then
    call ti_zero_arrays
    ti_nb_frc(:,:,:) = 0.d0
  end if

! Do weight changes, if requested.

  if (nmropt .ne. 0) call nmr_weight(atm_cnt, crd, 6)

! If no force calcs are to be done, clear the frc array and bag out now.

  if (ntf .eq. 8) then
    frc(:,:) = 0.d0
    return
  end if

  call zero_time()
  call zero_gb_time()

! Calculate the non-bonded contributions:

  frc(:,:) = 0.d0

  if (ti_mode .eq. 0) then
    call gb_ene(crd, frc, atm_gb_radii, atm_gb_fs, atm_qterm, atm_iac, &
                typ_ico, atm_numex, gbl_natex, atm_cnt, belly_atm_cnt, &
                pot_ene%gb, pot_ene%elec_tot, pot_ene%vdw_tot,         &
                pot_ene%surf, irespa)
  else
    ! NOTE: By calling gb_ene twice, this gives us at most 50% efficiency. 
    ! However, if we try to split the function (as with the
    ! other energy calculations), we'll have problems if the rborn of all of
    ! the atoms in the TI regions are not the same. This causes the reff
    ! of common atoms to be different in V0/V1, so we lose most of the speed
    ! advantage in the general case. 
    ti_mask_piece = 2
    call gb_ene(crd, frc, atm_gb_radii, atm_gb_fs, atm_qterm, atm_iac, &
                typ_ico, atm_numex, gbl_natex, atm_cnt, belly_atm_cnt, &
                ti_pot_gb(1, 1), ti_pot_gb(1, 2), ti_pot_gb(1, 3),     &
                ti_pot_gb(1, 4), irespa)
    ti_nb_frc(1,:,:) = frc(:,:) * ti_weights(1)
    frc(:,:) = 0.d0
    ti_mask_piece = 1
    call gb_ene(crd, frc, atm_gb_radii, atm_gb_fs, atm_qterm, atm_iac, &
                typ_ico, atm_numex, gbl_natex, atm_cnt, belly_atm_cnt, &
                ti_pot_gb(2, 1), ti_pot_gb(2, 2), ti_pot_gb(2, 3),     &
                ti_pot_gb(2, 4), irespa)
    frc(:,:) = ti_nb_frc(1,:,:) + frc(:,:) * ti_weights(2)

    call ti_update_ene_all(ti_pot_gb(:,1), si_hbond_ene, pot_ene%gb)
    call ti_update_ene_all(ti_pot_gb(:,2), si_elect_ene, pot_ene%elec_tot)
    call ti_update_ene_all(ti_pot_gb(:,3), si_vdw_ene,   pot_ene%vdw_tot)
    call ti_update_ene_all(ti_pot_gb(:,4), si_surf_ene,  pot_ene%surf)
  end if

! Calculate the 1-4 vdw and electrostatics contributions:

  if (charmm_active) then
!    call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
!                         gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
!                         pot_ene%elec_14, pot_ene%vdw_14)
  else
!    call get_nb14_energy(atm_qterm, crd, frc, atm_iac, typ_ico, &
!                         gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt, &
!                         pot_ene%elec_14, pot_ene%vdw_14)
  endif 

! For constant pH, get the EEL/VDW energies for the proposed state. This allows
! for easy extension to VDW modifications

  if (on_cpstep) then

    cnstph_frc(:,:) = 0.d0
    cph_egb         = 0.d0
    cph_eel         = 0.d0
    cph_eel14       = 0.d0
    cph_vdw         = 0.d0
    cph_vdw14       = 0.d0
    cph_surf        = 0.d0

    ! The last .true. argument forces gb_ene to skip calculating the effective
    ! radii, thereby saving time (since they will be unchanged from the previous
    ! gb_ene call)

    call gb_ene(crd, cnstph_frc, atm_gb_radii, atm_gb_fs, proposed_qterm, &
                atm_iac, typ_ico, atm_numex, gbl_natex, atm_cnt,          &
                belly_atm_cnt, cph_egb, cph_eel, cph_vdw, cph_surf, irespa, &
                .true.)

    if (charmm_active) then
!      call get_nb14_energy(proposed_qterm, crd, cnstph_frc, atm_iac, typ_ico, &
!                           gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt,      &
!                           cph_eel14, cph_vdw14)
    else
!      call get_nb14_energy(proposed_qterm, crd, cnstph_frc, atm_iac, typ_ico, &
!                           gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt,          &
!                           cph_eel14, cph_vdw14)
    end if

    ! Now compute dvdl = proposed state - current state

    pot_ene%dvdl = cph_egb + cph_eel + cph_vdw + cph_surf + &
                   cph_eel14 + cph_vdw14 - &
                   pot_ene%gb - pot_ene%elec_tot - pot_ene%vdw_tot - &
                   pot_ene%surf - pot_ene%elec_14 - pot_ene%vdw_14

  end if ! on_cpstep

  call update_time(nonbond_time)
  
  ! Calculate the other contributions:

  call gb_bonded_force(crd, frc, pot_ene)

  ! --- calculate the EMAP constraint energy ---

   if(iemap>0) then   ! (EMAP restraints)
       call emapforce(natom,enemap,atm_mass,crd,frc )
       pot_ene%emap = enemap
       pot_ene%restraint = pot_ene%restraint + enemap
   end if

  ! Sum up total potential energy for this task:

  pot_ene%total = pot_ene%vdw_tot + &
                  pot_ene%elec_tot + &
                  pot_ene%gb + &
                  pot_ene%surf + &
                  pot_ene%bond + &
                  pot_ene%angle + &
                  pot_ene%dihedral + &
                  pot_ene%vdw_14 + &
                  pot_ene%elec_14 + &
                  pot_ene%restraint + &
                  pot_ene%angle_ub + &
                  pot_ene%imp + &
                  pot_ene%cmap 


#ifdef MPI
! Distribute energies to all processes.

  call gb_distribute_enes(pot_ene)

  if (ti_mode .ne. 0) then
    call ti_dist_enes_virs_netfrcs(.true.,.false.)
  end if

#endif /* MPI */


 if(iamd.gt.1)then
 ! Calculate the boosting weight for dihedrals amd
   call calculate_amd_dih_weights(atm_cnt,pot_ene%dihedral,frc,crd)
 endif

#ifdef MPI
! Distribute forces to atom owners; this is a reduce sum operation:

  call gb_frcs_distrib(atm_cnt, frc, dbl_mpi_recv_buf)

  ! Clear the unused forces
  if (ti_mode .ne. 0) then
    do i = 1, atm_cnt
      if (gbl_atm_owner_map(i) .ne. mytaskid) frc(:,i) = 0.d0
    end do
  end if
#endif /* MPI */

  if (charmm_active .and. do_charmm_dump_gold == 1) then
#ifdef MPI
    call gb_mpi_gathervec(atm_cnt, frc)
    if (master) then
      call gb_charmm_dump_gold(atm_cnt,frc,pot_ene)
      write(mdout, '(a)') 'charmm_gold() completed. Exiting'
    endif
#else
    call gb_charmm_dump_gold(atm_cnt,frc,pot_ene)
    write(mdout, '(a)') 'charmm_gold() completed. Exiting'
#endif
    call mexit(6, 0)
  endif

! Calculate the NMR restraint energy contributions, if requested.
 
  if (nmropt .ne. 0) then
    call nmr_calc(crd, frc, enmr, 6)
    pot_ene%restraint = pot_ene%restraint + enmr(1) + enmr(2) + enmr(3)
    pot_ene%total = pot_ene%total + enmr(1) + enmr(2) + enmr(3)
  end if

  if (ti_mode .ne. 0) then
    call ti_calc_dvdl
  end if

!AMD DUAL BOOST CALC START
   if(iamd.gt.0)then
!calculate totboost and apply weight to frc. frc=frc*fwgt
     call calculate_amd_total_weights_gb(atm_cnt,pot_ene%total,pot_ene%dihedral,&
                                         pot_ene%amd_boost,frc,crd)
     ! Update total energy
     pot_ene%total = pot_ene%total + pot_ene%amd_boost
   end if

!scaledMD 
   if(scaledMD.gt.0)then
!apply scale to frc. frc=frc*scaledMD_lambda
     call scaledMD_scale_frc_gb(atm_cnt,pot_ene%total,frc,crd)
   end if

! If belly is on then set the belly atom forces to zero:

  if (ibelly .gt. 0) call bellyf(atm_cnt, atm_igroup, frc)
#endif /* CUDA */
  return

end subroutine gb_force

!*******************************************************************************
!
! Subroutine:  gb_bonded_force
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine gb_bonded_force(crd, frc, pot_ene)

  use angles_mod
  use angles_ub_mod
  use bonds_mod
  use constraints_mod
  use dihedrals_mod
  use dihedrals_imp_mod
  use cmap_mod
  use dynamics_dat_mod
  use mdin_ctrl_dat_mod
  use nmr_calls_mod
  use timers_mod

  implicit none

! Formal arguments:

  double precision              :: crd(3, *)
  double precision              :: frc(3, *)
  type(gb_pot_ene_rec)          :: pot_ene

! Local variables:

  ! These energy variables are temporaries, for summing. DON'T use otherwise!

  double precision              :: bond_ene
  double precision              :: angle_ene
  double precision              :: dihedral_ene
  double precision              :: ub_ene
  double precision              :: dihedral_imp_ene
  double precision              :: cmap_ene

  double precision              :: molvir(3, 3)         ! dummy argument
  double precision              :: e14_vir(3, 3)        ! dummy argument

! Bond energy contribution:

! The ebdev/eadev stuff currently only is output under nmr_calls for non-mpi
! code, so we basically drop it here under mpi.

#ifdef MPI
#else
  ebdev = 0.d0
#endif
  if (ntf .le. 1) then
!comment out temporarily by Lijiang Yang
!    if (cit_nbonh .gt. 0) then
!      call get_bond_energy(cit_nbonh, cit_h_bond, crd, frc, bond_ene)
!      pot_ene%bond = bond_ene
!    end if
  end if

  if (ntf .le. 2) then
!    if (cit_nbona .gt. 0) then
!comment out temporarily by Lijiang Yang
!      call get_bond_energy(cit_nbona, cit_a_bond, crd, frc, bond_ene)
!      pot_ene%bond = pot_ene%bond + bond_ene
!    end if
  end if
#ifdef MPI
#else
!    if (cit_nbonh + cit_nbona .gt. 0) &
!      ebdev = sqrt(ebdev / (cit_nbonh + cit_nbona))
#endif

  call update_time(bond_time)

! UB Angle energy contribution
  if (cit_nthet_ub .gt. 0) then
    call get_angle_ub_energy(cit_nthet_ub, cit_angle_ub, crd, frc, ub_ene)
    pot_ene%angle_ub = ub_ene
  endif

! Angle energy contribution:

#ifdef MPI
#else
  eadev = 0.d0
#endif
  if (ntf .le. 3) then
!comment out temporarily by Lijiang Yang
!    if (cit_ntheth .gt. 0) then
!      call get_angle_energy(cit_ntheth, cit_angle, crd, frc, angle_ene)
!      pot_ene%angle = angle_ene
!    end if
  end if

  if (ntf .le. 4) then
!comment out temporarily by Lijiang Yang
!    if (cit_ntheta .gt. 0) then
!      call get_angle_energy(cit_ntheta, cit_angle(cit_ntheth+1), &
!                            crd, frc, angle_ene)
!      pot_ene%angle = pot_ene%angle + angle_ene
!    end if
  end if
#ifdef MPI
#else
!  if (cit_ntheth + cit_ntheta .gt. 0) &
!    eadev = 57.296 * sqrt(eadev / (cit_ntheth + cit_ntheta))
#endif

  call update_time(angle_time)

! Improper dihedral energy contributions:

!comment out temporarily by Lijiang Yang
!  if (cit_nimphi .gt. 0) then
!    call get_dihed_imp_energy(cit_nimphi, cit_dihed_imp, crd, frc, dihedral_imp_ene)
!    pot_ene%imp = dihedral_imp_ene
!  endif

!CMAP
!  if (cit_cmap_term_count .gt. 0) then
!    call get_cmap_energy(cit_cmap_term_count, cit_cmap, crd, frc, cmap_ene)
!    pot_ene%cmap = cmap_ene
!  endif

  

! Dihedral energy contribution:

 if(iamd.gt.1)then
!comment out temporarily by Lijiang Yang
!  if (ntf .le. 5) then
!    if (cit_nphih .gt. 0) then
!      call get_dihed_energy_amd(cit_nphih, cit_dihed, crd, dihedral_ene)
!      pot_ene%dihedral = dihedral_ene
!    end if
!  end if

  if (ntf .le. 6) then
!comment out temporarily by Lijiang Yang
!    if (cit_nphia .gt. 0) then
!      call get_dihed_energy_amd(cit_nphia, cit_dihed(cit_nphih + 1), crd, &
!                            dihedral_ene)
!      pot_ene%dihedral = pot_ene%dihedral + dihedral_ene
!    end if
  end if
 else
  if (ntf .le. 5) then
!comment out temporarily by Lijiang Yang
!    if (cit_nphih .gt. 0) then
!      call get_dihed_energy(cit_nphih, cit_dihed, crd, frc, dihedral_ene)
!      pot_ene%dihedral = dihedral_ene
!    end if
  end if

  if (ntf .le. 6) then
!comment out temporarily by Lijiang Yang
!    if (cit_nphia .gt. 0) then
!      call get_dihed_energy(cit_nphia, cit_dihed(cit_nphih + 1), crd, frc, &
!                            dihedral_ene)
!      pot_ene%dihedral = pot_ene%dihedral + dihedral_ene
!    end if
  end if
 endif

  call update_time(dihedral_time)

! Calculate the position constraint energy:

  if (natc .gt. 0) then
     call get_crd_constraint_energy(natc, pot_ene%restraint, atm_jrc, &
                                    crd, frc, atm_xc, atm_weight)
  endif
  return

end subroutine gb_bonded_force

!*******************************************************************************
!
! Subroutine:  gb_cph_ene
!
! Description: Gets CpHMD dvdl term
!              
!*******************************************************************************

subroutine gb_cph_ene(atm_cnt, crd, qterm, dvdl, skipradii)

  use charmm_mod, only          : charmm_active
  use constantph_dat_mod, only  : cnstph_frc
  use constraints_mod, only     : belly_atm_cnt
  use extra_pnts_nb14_mod, only : cit_nb14, cit_nb14_cnt, get_nb14_energy
  use gb_ene_mod, only          : gb_ene
  use parallel_dat_mod
  use prmtop_dat_mod, only      : atm_gb_radii, atm_gb_fs, atm_iac, typ_ico, &
                                  atm_numex, gbl_natex, gbl_cn114, gbl_cn214, &
                                  gbl_cn1, gbl_cn2

  implicit none

  ! Passed arguments

  integer, intent(in)             :: atm_cnt
  double precision, intent(in)    :: crd(3, atm_cnt)
  double precision, intent(in)    :: qterm(atm_cnt)
  double precision, intent(out)   :: dvdl
  logical, intent(in)             :: skipradii

  ! Local arguments

  double precision :: egb = 0.d0
  double precision :: eel = 0.d0
  double precision :: dvdl_buf = 0.d0
  double precision :: surf = 0.d0
  double precision :: vdw = 0.d0
  double precision :: eel14 = 0.d0
  double precision :: vdw14 = 0.d0

  ! Compute the energy of the current state

  call gb_ene(crd, cnstph_frc, atm_gb_radii, atm_gb_fs, qterm,   &
              atm_iac, typ_ico, atm_numex, gbl_natex, atm_cnt,   &
              belly_atm_cnt, egb, eel, vdw, surf, 0, skipradii)

  if (charmm_active) then
!    call get_nb14_energy(qterm, crd, cnstph_frc, atm_iac, typ_ico, &
!                         gbl_cn114, gbl_cn214, cit_nb14, cit_nb14_cnt, &
!                         eel14, vdw14)
  else
!    call get_nb14_energy(qterm, crd, cnstph_frc, atm_iac, typ_ico, &
!                         gbl_cn1, gbl_cn2, cit_nb14, cit_nb14_cnt,     &
!                         eel14, vdw14)
  end if

  ! Store the constant pH-relevant portion of the current state energy

  dvdl = egb + eel + vdw + surf + eel14 + vdw14

#ifdef MPI
  ! We need to reduce our total energy in parallel
  call mpi_allreduce(dvdl, dvdl_buf, 1, mpi_double_precision, mpi_sum, &
                     pmemd_comm, err_code_mpi)

  dvdl = dvdl_buf
#endif

! if (master) write (666, '(5(a,f15.8,x))') &
!     'EGB = ', egb, 'EEL = ', eel, 'VDW = ', vdw, 'EEL 1-4 = ', eel14, &
!     'VDW 1-4 = ', vdw14

  return

end subroutine gb_cph_ene

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  gb_distribute_enes
!
! Description: <TBS>
!
!*******************************************************************************

subroutine gb_distribute_enes(pot_ene)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  type(gb_pot_ene_rec) :: pot_ene

! Local variables:

  type(gb_pot_ene_rec), save    :: dat_in, dat_out

  dat_in = pot_ene

  call mpi_allreduce(dat_in%total, dat_out%total, &
                     gb_pot_ene_rec_size, mpi_double_precision, &
                     mpi_sum, pmemd_comm, err_code_mpi)

  pot_ene = dat_out

  return

end subroutine gb_distribute_enes
#endif

end module gb_force_mod
