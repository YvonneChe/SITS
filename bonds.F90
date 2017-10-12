#include "copyright.i"

!*******************************************************************************
!
! Module:  bonds_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module bonds_mod

  use gbl_datatypes_mod

  implicit none

! The following are derived from prmtop bond info:

! Added by Lijiang Yang
  integer, save                         :: cit_nbonhslu, cit_nbonaslu
  integer, save                         :: cit_nbonhsol, cit_nbonasol
  integer, save                         :: cit_nbonhinter, cit_nbonainter
!-----------------------------------------------------------------

! Added by Lijiang Yang
  type(bond_rec), allocatable, save     :: cit_h_bondslu(:)
  type(bond_rec), allocatable, save     :: cit_a_bondslu(:)
  type(bond_rec), allocatable, save     :: cit_h_bondsol(:)
  type(bond_rec), allocatable, save     :: cit_a_bondsol(:)
  type(bond_rec), allocatable, save     :: cit_h_bondinter(:)
  type(bond_rec), allocatable, save     :: cit_a_bondinter(:)
!-----------------------------------------------------------------

contains

!*******************************************************************************
!
! Subroutine:  bonds_setup
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine bonds_setup(num_ints, num_reals, use_atm_map)

  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals
  integer                       :: use_atm_map(natom)

! Local variables:

  integer               :: alloc_failed
! Added by Lijiang Yang
  type(bond_rec)        :: bonds_copyslu(nbonhslu + nbonaslu)
  type(bond_rec)        :: bonds_copysol(nbonhsol + nbonasol)
  type(bond_rec)        :: bonds_copyinter(nbonhinter + nbonainter)
!-----------------------------------------------------------------

  ! This routine can handle reallocation, and thus can be called multiple
  ! times.

! Modified by Lijiang Yang
  call find_my_bonds(nbonhslu, gbl_bondslu, cit_nbonhslu, bonds_copyslu, use_atm_map)
  call find_my_bonds(nbonhsol, gbl_bondsol, cit_nbonhsol, bonds_copysol, use_atm_map)
  call find_my_bonds(nbonhinter, gbl_bondinter, cit_nbonhinter, bonds_copyinter, use_atm_map)

  if (cit_nbonhslu .gt. 0) then
    if (allocated(cit_h_bondslu)) then
      if (size(cit_h_bondslu) .lt. cit_nbonhslu) then
        num_ints = num_ints - size(cit_h_bondslu) * bond_rec_ints
        deallocate(cit_h_bondslu)
       allocate(cit_h_bondslu(cit_nbonhslu), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_h_bondslu) * bond_rec_ints
      end if
    else
      allocate(cit_h_bondslu(cit_nbonhslu), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_h_bondslu) * bond_rec_ints
    end if
    cit_h_bondslu(1:cit_nbonhslu) = bonds_copyslu(1:cit_nbonhslu)
  end if
  if (cit_nbonhsol .gt. 0) then
    if (allocated(cit_h_bondsol)) then
      if (size(cit_h_bondsol) .lt. cit_nbonhsol) then
       num_ints = num_ints - size(cit_h_bondsol) * bond_rec_ints
        deallocate(cit_h_bondsol)
        allocate(cit_h_bondsol(cit_nbonhsol), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_h_bondsol) * bond_rec_ints
      end if
    else
      allocate(cit_h_bondsol(cit_nbonhsol), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_h_bondsol) * bond_rec_ints
    end if
    cit_h_bondsol(1:cit_nbonhsol) = bonds_copysol(1:cit_nbonhsol)
  end if
  if (cit_nbonhinter .gt. 0) then
    if (allocated(cit_h_bondinter)) then
      if (size(cit_h_bondinter) .lt. cit_nbonhinter) then
        num_ints = num_ints - size(cit_h_bondinter) * bond_rec_ints
        deallocate(cit_h_bondinter)
        allocate(cit_h_bondinter(cit_nbonhinter), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_h_bondinter) * bond_rec_ints
      end if
    else
      allocate(cit_h_bondinter(cit_nbonhinter), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_h_bondinter) * bond_rec_ints
    end if
    cit_h_bondinter(1:cit_nbonhinter) = bonds_copyinter(1:cit_nbonhinter)
  end if

! Added by Lijiang Yang
  call find_my_bonds(nbonaslu, gbl_bondslu(bonda_idx_slu), cit_nbonaslu, bonds_copyslu, use_atm_map)
  call find_my_bonds(nbonasol, gbl_bondsol(bonda_idx_sol), cit_nbonasol, bonds_copysol, use_atm_map)
  call find_my_bonds(nbonainter, gbl_bondinter(bonda_idx_inter), cit_nbonainter, bonds_copyinter, use_atm_map)
!-----------------------------------------------------------------

  if (cit_nbonaslu .gt. 0) then
    if (allocated(cit_a_bondslu)) then
      if (size(cit_a_bondslu) .lt. cit_nbonaslu) then
        num_ints = num_ints - size(cit_a_bondslu) * bond_rec_ints
        deallocate(cit_a_bondslu)
        allocate(cit_a_bondslu(cit_nbonaslu), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_a_bondslu) * bond_rec_ints
      end if
    else
      allocate(cit_a_bondslu(cit_nbonaslu), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_a_bondslu) * bond_rec_ints
    end if
    cit_a_bondslu(1:cit_nbonaslu) = bonds_copyslu(1:cit_nbonaslu)
  end if
  if (cit_nbonasol .gt. 0) then
    if (allocated(cit_a_bondsol)) then
      if (size(cit_a_bondsol) .lt. cit_nbonasol) then
        num_ints = num_ints - size(cit_a_bondsol) * bond_rec_ints
        deallocate(cit_a_bondsol)
        allocate(cit_a_bondsol(cit_nbonasol), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_a_bondsol) * bond_rec_ints
      end if
    else
      allocate(cit_a_bondsol(cit_nbonasol), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_a_bondsol) * bond_rec_ints
    end if
    cit_a_bondsol(1:cit_nbonasol) = bonds_copysol(1:cit_nbonasol)
  end if
  if (cit_nbonainter .gt. 0) then
    if (allocated(cit_a_bondinter)) then
      if (size(cit_a_bondinter) .lt. cit_nbonainter) then
        num_ints = num_ints - size(cit_a_bondinter) * bond_rec_ints
        deallocate(cit_a_bondinter)
        allocate(cit_a_bondinter(cit_nbonainter), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_a_bondinter) * bond_rec_ints
      end if
    else
      allocate(cit_a_bondinter(cit_nbonainter), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_a_bondinter) * bond_rec_ints
    end if
    cit_a_bondinter(1:cit_nbonainter) = bonds_copyinter(1:cit_nbonainter)
  end if
!-----------------------------------------------------------------

#ifdef CUDA
  call gpu_bonds_setup(cit_nbona, cit_a_bond, cit_nbonh, cit_h_bond, gbl_req, gbl_rk);    
#endif

  return

end subroutine bonds_setup

!*******************************************************************************
!
! Subroutine:  find_my_bonds
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine find_my_bonds(bond_cnt, bonds, my_bond_cnt, my_bonds, &
                             use_atm_map)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: bond_cnt
  type(bond_rec)        :: bonds(bond_cnt)
  integer               :: my_bond_cnt
  type(bond_rec)        :: my_bonds(*)
  integer               :: use_atm_map(*)

! Local variables:

  integer               :: atm_i, atm_j, bonds_idx

! Find all bonds for which this process owns either atom:

  my_bond_cnt = 0

  do bonds_idx = 1, bond_cnt

    atm_i = bonds(bonds_idx)%atm_i
    atm_j = bonds(bonds_idx)%atm_j

#if defined(MPI) && !defined(CUDA)
    if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
      my_bond_cnt = my_bond_cnt + 1
      my_bonds(my_bond_cnt) = bonds(bonds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
      my_bond_cnt = my_bond_cnt + 1
      my_bonds(my_bond_cnt) = bonds(bonds_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    end if
#else
    my_bond_cnt = my_bond_cnt + 1
    my_bonds(my_bond_cnt) = bonds(bonds_idx)
    use_atm_map(atm_i) = 1
    use_atm_map(atm_j) = 1
#endif

  end do

  return

end subroutine find_my_bonds


!*******************************************************************************
!
! Subroutine:  get_bond_energy
!
! Description:
!              
! Routine to get bond energy and forces for the potential of cb*(b-b0)**2.
!
!*******************************************************************************

subroutine get_bond_energy(bond_cnt, bond, x, frc, subfrc, bond_energy)

  use nmr_calls_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: bond_cnt
  type(bond_rec)        :: bond(*)
  double precision      :: x(3, *)
  double precision      :: frc(3, *)
  double precision      :: subfrc(3, *)
  double precision      :: bond_energy

! Local variables:

  double precision      :: da
  double precision      :: df
  double precision      :: dfw
  integer               :: i, j, ic, jn
  double precision      :: lcl_bond_energy
  double precision      :: xa, ya, za
  double precision      :: rij
  double precision      :: xij, yij, zij

  lcl_bond_energy = 0.0d0

! Grand loop for the bond stuff:

  do jn = 1, bond_cnt

    i = bond(jn)%atm_i
    j = bond(jn)%atm_j

! Calculation of the bond vector:

    xij = x(1, i) - x(1, j)
    yij = x(2, i) - x(2, j)
    zij = x(3, i) - x(3, j)

    rij = sqrt(xij * xij + yij * yij + zij * zij)

! Calculation of the energy and deriv:

    ic = bond(jn)%parm_idx
    da = rij - gbl_req(ic)

#ifdef MPI
! If ebdev is ever supported under mpi, you will need to treat it like the
! energy...
#else
    ebdev = ebdev + da * da ! For rms deviation from ideal bonds:
#endif
    
    df = gbl_rk(ic) * da

    if (ti_mode .ne. 0) then
      ti_region = 0
      if (ti_lst(1,i) .eq. 1 .or. ti_lst(1,j) .eq. 1) then
        ti_region = 1
      else if (ti_lst(2,i) .eq. 1 .or. ti_lst(2,j) .eq. 1) then
        ti_region = 2
      end if
      if (ti_region .gt. 0) then
        if (ti_mode .ne. 1) then
          ti_lscale = (ti_sc_lst(i) .lt. 2 .and. &
                       ti_sc_lst(j) .lt. 2) 
        end if
#ifdef MPI
        if (gbl_atm_owner_map(i) .eq. mytaskid) then
#endif    
          if (ti_mode .eq. 1 .or. ti_lscale) then !linear scaling      
            call ti_update_ene(df * da, si_bond_ene, ti_region)
            !no need to scale da as we always scale df
          else
            ti_ene(ti_region,si_bond_ene) = &
              ti_ene(ti_region,si_bond_ene) + df * da
            da = 0.d0
          end if
#ifdef MPI
        end if
#endif
        if (ti_mode .eq. 1 .or. ti_lscale) df = df * ti_weights(ti_region)
      end if
    end if

    dfw = (df + df) / rij

! Calculation of the force:

    xa = dfw * xij
    ya = dfw * yij
    za = dfw * zij

#ifdef MPI
    ! We use atm_i to determine who sums up the energy...
    if (gbl_atm_owner_map(i) .eq. mytaskid) then
      lcl_bond_energy = lcl_bond_energy + df * da
!      frc(1, i) = frc(1, i) - xa
!      frc(2, i) = frc(2, i) - ya
!      frc(3, i) = frc(3, i) - za
      subfrc(1, i) = subfrc(1, i) - xa
      subfrc(2, i) = subfrc(2, i) - ya
      subfrc(3, i) = subfrc(3, i) - za
    end if

    if (gbl_atm_owner_map(j) .eq. mytaskid) then
!      frc(1, j) = frc(1, j) + xa
!      frc(2, j) = frc(2, j) + ya
!      frc(3, j) = frc(3, j) + za
      subfrc(1, j) = subfrc(1, j) + xa
      subfrc(2, j) = subfrc(2, j) + ya
      subfrc(3, j) = subfrc(3, j) + za
    end if
#else
      lcl_bond_energy = lcl_bond_energy + df * da
!      frc(1, i) = frc(1, i) - xa
!      frc(2, i) = frc(2, i) - ya
!      frc(3, i) = frc(3, i) - za
!      frc(1, j) = frc(1, j) + xa
!      frc(2, j) = frc(2, j) + ya
!      frc(3, j) = frc(3, j) + za
      subfrc(1, i) = subfrc(1, i) - xa
      subfrc(2, i) = subfrc(2, i) - ya
      subfrc(3, i) = subfrc(3, i) - za
      subfrc(1, j) = subfrc(1, j) + xa
      subfrc(2, j) = subfrc(2, j) + ya
      subfrc(3, j) = subfrc(3, j) + za
#endif

  end do

  bond_energy = lcl_bond_energy

  return

end subroutine get_bond_energy

end module bonds_mod
