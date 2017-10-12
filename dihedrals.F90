#include "copyright.i"

!*******************************************************************************
!
! Module:  dihedrals_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module dihedrals_mod

  use gbl_datatypes_mod

  implicit none

! The following are derived from prmtop dihedral info:

  integer, save                         :: cit_nphihslu, cit_nphiaslu
  integer, save                         :: cit_nphihsol, cit_nphiasol
  integer, save                         :: cit_nphihinter, cit_nphiainter

!AMD
  double precision, save                :: fwgtd

  type(dihed_rec), allocatable, save    :: cit_dihedslu(:)
  type(dihed_rec), allocatable, save    :: cit_dihedsol(:)
  type(dihed_rec), allocatable, save    :: cit_dihedinter(:)

contains

!*******************************************************************************
!
! Subroutine:  dihedrals_setup
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine dihedrals_setup(num_ints, num_reals, use_atm_map)

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
  type(dihed_rec)       :: diheds_copyslu(nphihslu + nphiaslu)
  type(dihed_rec)       :: diheds_copysol(nphihsol + nphiasol)
  type(dihed_rec)       :: diheds_copyinter(nphihinter + nphiainter)
  integer               :: atm_i, atm_j, atm_k, atm_l, diheds_idx
  integer               :: my_dihed_cnt_slu, my_dihed_cnt_sol, my_dihed_cnt_inter

! This routine can handle reallocation, and thus can be called multiple
! times.

! Find all diheds for which this process owns either atom:

  my_dihed_cnt_slu = 0

  if ( nphihslu .gt. 0 ) then
    do diheds_idx = 1, nphihslu

      atm_i = gbl_dihedslu(diheds_idx)%atm_i
      atm_j = gbl_dihedslu(diheds_idx)%atm_j
      atm_k = iabs(gbl_dihedslu(diheds_idx)%atm_k)
      atm_l = iabs(gbl_dihedslu(diheds_idx)%atm_l)

#if defined(MPI) && !defined(CUDA)
      if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
        my_dihed_cnt_slu = my_dihed_cnt_slu + 1
        diheds_copyslu(my_dihed_cnt_slu) = gbl_dihedslu(diheds_idx)
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_k = gbl_dihedslu(diheds_idx)%atm_k
        !diheds_copyslu(my_dihed_cnt_slu)%atm_l = gbl_dihedslu(diheds_idx)%atm_l
        !diheds_copyslu(my_dihed_cnt_slu)%parm_idx = gbl_dihedslu(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
        my_dihed_cnt_slu = my_dihed_cnt_slu + 1
        diheds_copyslu(my_dihed_cnt_slu) = gbl_dihedslu(diheds_idx)
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_k = gbl_dihedslu(diheds_idx)%atm_k
        !diheds_copyslu(my_dihed_cnt_slu)%atm_l = gbl_dihedslu(diheds_idx)%atm_l
        !diheds_copyslu(my_dihed_cnt_slu)%parm_idx = gbl_dihedslu(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
        my_dihed_cnt_slu = my_dihed_cnt_slu + 1
        diheds_copyslu(my_dihed_cnt_slu) = gbl_dihedslu(diheds_idx)
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_k = gbl_dihedslu(diheds_idx)%atm_k
        !diheds_copyslu(my_dihed_cnt_slu)%atm_l = gbl_dihedslu(diheds_idx)%atm_l
        !diheds_copyslu(my_dihed_cnt_slu)%parm_idx = gbl_dihedslu(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_l) .eq. mytaskid) then
        my_dihed_cnt_slu = my_dihed_cnt_slu + 1
        diheds_copyslu(my_dihed_cnt_slu) = gbl_dihedslu(diheds_idx)
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_k = gbl_dihedslu(diheds_idx)%atm_k
        !diheds_copyslu(my_dihed_cnt_slu)%atm_l = gbl_dihedslu(diheds_idx)%atm_l
        !diheds_copyslu(my_dihed_cnt_slu)%parm_idx = gbl_dihedslu(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      end if
#else
      my_dihed_cnt_slu = my_dihed_cnt_slu + 1
      diheds_copyslu(my_dihed_cnt_slu) = gbl_dihedslu(diheds_idx)
      !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
      !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
      !diheds_copyslu(my_dihed_cnt_slu)%atm_k = gbl_dihedslu(diheds_idx)%atm_k
      !diheds_copyslu(my_dihed_cnt_slu)%atm_l = gbl_dihedslu(diheds_idx)%atm_l
      !diheds_copyslu(my_dihed_cnt_slu)%parm_idx = gbl_dihedslu(diheds_idx)%parm_idx
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
#endif

    end do

    cit_nphihslu = my_dihed_cnt_slu
  endif

  if ( nphiaslu .gt. 0 ) then
    do diheds_idx = diheda_idx_slu, diheda_idx_slu + nphiaslu - 1

      atm_i = gbl_dihedslu(diheds_idx)%atm_i
      atm_j = gbl_dihedslu(diheds_idx)%atm_j
      atm_k = iabs(gbl_dihedslu(diheds_idx)%atm_k)
      atm_l = iabs(gbl_dihedslu(diheds_idx)%atm_l)

#if defined(MPI) && !defined(CUDA)  
      if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
        my_dihed_cnt_slu = my_dihed_cnt_slu + 1
        diheds_copyslu(my_dihed_cnt_slu) = gbl_dihedslu(diheds_idx)
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_k = gbl_dihedslu(diheds_idx)%atm_k
        !diheds_copyslu(my_dihed_cnt_slu)%atm_l = gbl_dihedslu(diheds_idx)%atm_l
        !diheds_copyslu(my_dihed_cnt_slu)%parm_idx = gbl_dihedslu(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
        my_dihed_cnt_slu = my_dihed_cnt_slu + 1
        diheds_copyslu(my_dihed_cnt_slu) = gbl_dihedslu(diheds_idx)
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_k = gbl_dihedslu(diheds_idx)%atm_k
        !diheds_copyslu(my_dihed_cnt_slu)%atm_l = gbl_dihedslu(diheds_idx)%atm_l
        !diheds_copyslu(my_dihed_cnt_slu)%parm_idx = gbl_dihedslu(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
        my_dihed_cnt_slu = my_dihed_cnt_slu + 1
        diheds_copyslu(my_dihed_cnt_slu) = gbl_dihedslu(diheds_idx)
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_k = gbl_dihedslu(diheds_idx)%atm_k
        !diheds_copyslu(my_dihed_cnt_slu)%atm_l = gbl_dihedslu(diheds_idx)%atm_l
        !diheds_copyslu(my_dihed_cnt_slu)%parm_idx = gbl_dihedslu(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_l) .eq. mytaskid) then
        my_dihed_cnt_slu = my_dihed_cnt_slu + 1
        diheds_copyslu(my_dihed_cnt_slu) = gbl_dihedslu(diheds_idx)
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
        !diheds_copyslu(my_dihed_cnt_slu)%atm_k = gbl_dihedslu(diheds_idx)%atm_k
        !diheds_copyslu(my_dihed_cnt_slu)%atm_l = gbl_dihedslu(diheds_idx)%atm_l
        !diheds_copyslu(my_dihed_cnt_slu)%parm_idx = gbl_dihedslu(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      end if
#else
      my_dihed_cnt_slu = my_dihed_cnt_slu + 1
      diheds_copyslu(my_dihed_cnt_slu) = gbl_dihedslu(diheds_idx)
      !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
      !diheds_copyslu(my_dihed_cnt_slu)%atm_j = gbl_dihedslu(diheds_idx)%atm_j
      !diheds_copyslu(my_dihed_cnt_slu)%atm_k = gbl_dihedslu(diheds_idx)%atm_k
      !diheds_copyslu(my_dihed_cnt_slu)%atm_l = gbl_dihedslu(diheds_idx)%atm_l
      !diheds_copyslu(my_dihed_cnt_slu)%parm_idx = gbl_dihedslu(diheds_idx)%parm_idx
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
#endif

    end do

    cit_nphiaslu = my_dihed_cnt_slu - cit_nphihslu
  endif

  if (my_dihed_cnt_slu .gt. 0) then
    if (allocated(cit_dihedslu)) then
      if (size(cit_dihedslu) .lt. my_dihed_cnt_slu) then
        num_ints = num_ints - size(cit_dihedslu) * dihed_rec_ints
        deallocate(cit_dihedslu)
        allocate(cit_dihedslu(my_dihed_cnt_slu), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_dihedslu) * dihed_rec_ints
      end if
    else
      allocate(cit_dihedslu(my_dihed_cnt_slu), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_dihedslu) * dihed_rec_ints
    end if
    cit_dihedslu(1:my_dihed_cnt_slu) = diheds_copyslu(1:my_dihed_cnt_slu)
    !do diheds_idx = 1, my_dihed_cnt_slu 
    !   cit_dihedslu(diheds_idx)%atm_i = diheds_copyslu(diheds_idx)%atm_i
    !   cit_dihedslu(diheds_idx)%atm_j = diheds_copyslu(diheds_idx)%atm_j
    !   cit_dihedslu(diheds_idx)%atm_k = diheds_copyslu(diheds_idx)%atm_k
    !   cit_dihedslu(diheds_idx)%atm_l = diheds_copyslu(diheds_idx)%atm_l
    !   cit_dihedslu(diheds_idx)%parm_idx = diheds_copyslu(diheds_idx)%parm_idx
    !end do
  end if

  my_dihed_cnt_sol = 0

  if ( nphihsol .gt. 0 ) then
    do diheds_idx = 1, nphihsol

      atm_i = gbl_dihedsol(diheds_idx)%atm_i
      atm_j = gbl_dihedsol(diheds_idx)%atm_j
      atm_k = iabs(gbl_dihedsol(diheds_idx)%atm_k)
      atm_l = iabs(gbl_dihedsol(diheds_idx)%atm_l)

#if defined(MPI) && !defined(CUDA)
      if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
        my_dihed_cnt_sol = my_dihed_cnt_sol + 1
        !Modified by Lijiang Yang: it is not right if we write the assignment
        !statement like the following. This is very wierd!
        !diheds_copysol(my_dihed_cnt_sol) = gbl_dihedsol(diheds_idx) 
        !-----------------------------------------------------------------------
        diheds_copysol(my_dihed_cnt_sol)%atm_i = gbl_dihedsol(diheds_idx)%atm_i
        diheds_copysol(my_dihed_cnt_sol)%atm_j = gbl_dihedsol(diheds_idx)%atm_j
        diheds_copysol(my_dihed_cnt_sol)%atm_k = gbl_dihedsol(diheds_idx)%atm_k
        diheds_copysol(my_dihed_cnt_sol)%atm_l = gbl_dihedsol(diheds_idx)%atm_l
        diheds_copysol(my_dihed_cnt_sol)%parm_idx = gbl_dihedsol(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
        my_dihed_cnt_sol = my_dihed_cnt_sol + 1
        !diheds_copysol(my_dihed_cnt_sol) = gbl_dihedsol(diheds_idx)
        diheds_copysol(my_dihed_cnt_sol)%atm_i = gbl_dihedsol(diheds_idx)%atm_i
        diheds_copysol(my_dihed_cnt_sol)%atm_j = gbl_dihedsol(diheds_idx)%atm_j
        diheds_copysol(my_dihed_cnt_sol)%atm_k = gbl_dihedsol(diheds_idx)%atm_k
        diheds_copysol(my_dihed_cnt_sol)%atm_l = gbl_dihedsol(diheds_idx)%atm_l
        diheds_copysol(my_dihed_cnt_sol)%parm_idx = gbl_dihedsol(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
        my_dihed_cnt_sol = my_dihed_cnt_sol + 1
        !diheds_copysol(my_dihed_cnt_sol) = gbl_dihedsol(diheds_idx)
        diheds_copysol(my_dihed_cnt_sol)%atm_i = gbl_dihedsol(diheds_idx)%atm_i
        diheds_copysol(my_dihed_cnt_sol)%atm_j = gbl_dihedsol(diheds_idx)%atm_j
        diheds_copysol(my_dihed_cnt_sol)%atm_k = gbl_dihedsol(diheds_idx)%atm_k
        diheds_copysol(my_dihed_cnt_sol)%atm_l = gbl_dihedsol(diheds_idx)%atm_l
        diheds_copysol(my_dihed_cnt_sol)%parm_idx = gbl_dihedsol(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_l) .eq. mytaskid) then
        my_dihed_cnt_sol = my_dihed_cnt_sol + 1
        !diheds_copyslu(my_dihed_cnt_sol) = gbl_dihedsol(diheds_idx)
        diheds_copysol(my_dihed_cnt_sol)%atm_i = gbl_dihedsol(diheds_idx)%atm_i
        diheds_copysol(my_dihed_cnt_sol)%atm_j = gbl_dihedsol(diheds_idx)%atm_j
        diheds_copysol(my_dihed_cnt_sol)%atm_k = gbl_dihedsol(diheds_idx)%atm_k
        diheds_copysol(my_dihed_cnt_sol)%atm_l = gbl_dihedsol(diheds_idx)%atm_l
        diheds_copysol(my_dihed_cnt_sol)%parm_idx = gbl_dihedsol(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      end if
#else
      my_dihed_cnt_sol = my_dihed_cnt_sol + 1
      !diheds_copysol(my_dihed_cnt_sol) = gbl_dihedsol(diheds_idx)
      diheds_copysol(my_dihed_cnt_sol)%atm_i = gbl_dihedsol(diheds_idx)%atm_i
      diheds_copysol(my_dihed_cnt_sol)%atm_j = gbl_dihedsol(diheds_idx)%atm_j
      diheds_copysol(my_dihed_cnt_sol)%atm_k = gbl_dihedsol(diheds_idx)%atm_k
      diheds_copysol(my_dihed_cnt_sol)%atm_l = gbl_dihedsol(diheds_idx)%atm_l
      diheds_copysol(my_dihed_cnt_sol)%parm_idx = gbl_dihedsol(diheds_idx)%parm_idx
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
#endif

    end do

    cit_nphihsol = my_dihed_cnt_sol
  endif

  if ( nphiasol .gt. 0 ) then
    do diheds_idx = diheda_idx_sol, diheda_idx_sol + nphiasol - 1

      atm_i = gbl_dihedsol(diheds_idx)%atm_i
      atm_j = gbl_dihedsol(diheds_idx)%atm_j
      atm_k = iabs(gbl_dihedsol(diheds_idx)%atm_k)
      atm_l = iabs(gbl_dihedsol(diheds_idx)%atm_l)

#if defined(MPI) && !defined(CUDA)  
      if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
        my_dihed_cnt_sol = my_dihed_cnt_sol + 1
        !diheds_copysol(my_dihed_cnt_sol) = gbl_dihedsol(diheds_idx)
        diheds_copysol(my_dihed_cnt_sol)%atm_i = gbl_dihedsol(diheds_idx)%atm_i
        diheds_copysol(my_dihed_cnt_sol)%atm_j = gbl_dihedsol(diheds_idx)%atm_j
        diheds_copysol(my_dihed_cnt_sol)%atm_k = gbl_dihedsol(diheds_idx)%atm_k
        diheds_copysol(my_dihed_cnt_sol)%atm_l = gbl_dihedsol(diheds_idx)%atm_l
        diheds_copysol(my_dihed_cnt_sol)%parm_idx = gbl_dihedsol(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
        my_dihed_cnt_sol = my_dihed_cnt_sol + 1
        !diheds_copysol(my_dihed_cnt_sol) = gbl_dihedsol(diheds_idx)
        diheds_copysol(my_dihed_cnt_sol)%atm_i = gbl_dihedsol(diheds_idx)%atm_i
        diheds_copysol(my_dihed_cnt_sol)%atm_j = gbl_dihedsol(diheds_idx)%atm_j
        diheds_copysol(my_dihed_cnt_sol)%atm_k = gbl_dihedsol(diheds_idx)%atm_k
        diheds_copysol(my_dihed_cnt_sol)%atm_l = gbl_dihedsol(diheds_idx)%atm_l
        diheds_copysol(my_dihed_cnt_sol)%parm_idx = gbl_dihedsol(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
        my_dihed_cnt_sol = my_dihed_cnt_sol + 1
        !diheds_copysol(my_dihed_cnt_sol) = gbl_dihedsol(diheds_idx)
        diheds_copysol(my_dihed_cnt_sol)%atm_i = gbl_dihedsol(diheds_idx)%atm_i
        diheds_copysol(my_dihed_cnt_sol)%atm_j = gbl_dihedsol(diheds_idx)%atm_j
        diheds_copysol(my_dihed_cnt_sol)%atm_k = gbl_dihedsol(diheds_idx)%atm_k
        diheds_copysol(my_dihed_cnt_sol)%atm_l = gbl_dihedsol(diheds_idx)%atm_l
        diheds_copysol(my_dihed_cnt_sol)%parm_idx = gbl_dihedsol(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_l) .eq. mytaskid) then
        my_dihed_cnt_sol = my_dihed_cnt_sol + 1
        !diheds_copysol(my_dihed_cnt_sol) = gbl_dihedsol(diheds_idx)
        diheds_copysol(my_dihed_cnt_sol)%atm_i = gbl_dihedsol(diheds_idx)%atm_i
        diheds_copysol(my_dihed_cnt_sol)%atm_j = gbl_dihedsol(diheds_idx)%atm_j
        diheds_copysol(my_dihed_cnt_sol)%atm_k = gbl_dihedsol(diheds_idx)%atm_k
        diheds_copysol(my_dihed_cnt_sol)%atm_l = gbl_dihedsol(diheds_idx)%atm_l
        diheds_copysol(my_dihed_cnt_sol)%parm_idx = gbl_dihedsol(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      end if
#else
      my_dihed_cnt_sol = my_dihed_cnt_sol + 1
      !diheds_copysol(my_dihed_cnt_sol) = gbl_dihedsol(diheds_idx)
      diheds_copysol(my_dihed_cnt_sol)%atm_i = gbl_dihedsol(diheds_idx)%atm_i
      diheds_copysol(my_dihed_cnt_sol)%atm_j = gbl_dihedsol(diheds_idx)%atm_j
      diheds_copysol(my_dihed_cnt_sol)%atm_k = gbl_dihedsol(diheds_idx)%atm_k
      diheds_copysol(my_dihed_cnt_sol)%atm_l = gbl_dihedsol(diheds_idx)%atm_l
      diheds_copysol(my_dihed_cnt_sol)%parm_idx = gbl_dihedsol(diheds_idx)%parm_idx
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
#endif

    end do

    cit_nphiasol = my_dihed_cnt_sol - cit_nphihsol
  endif

  if (my_dihed_cnt_sol .gt. 0) then
    if (allocated(cit_dihedsol)) then
      if (size(cit_dihedsol) .lt. my_dihed_cnt_sol) then
        num_ints = num_ints - size(cit_dihedsol) * dihed_rec_ints
        deallocate(cit_dihedsol)
        allocate(cit_dihedsol(my_dihed_cnt_sol), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_dihedsol) * dihed_rec_ints
      end if
    else
      allocate(cit_dihedsol(my_dihed_cnt_sol), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_dihedsol) * dihed_rec_ints
    end if
    !cit_dihedsol(1:my_dihed_cnt_sol) = diheds_copysol(1:my_dihed_cnt_sol)
    do diheds_idx = 1, my_dihed_cnt_sol 
       cit_dihedsol(diheds_idx)%atm_i = diheds_copysol(diheds_idx)%atm_i
       cit_dihedsol(diheds_idx)%atm_j = diheds_copysol(diheds_idx)%atm_j
       cit_dihedsol(diheds_idx)%atm_k = diheds_copysol(diheds_idx)%atm_k
       cit_dihedsol(diheds_idx)%atm_l = diheds_copysol(diheds_idx)%atm_l
       cit_dihedsol(diheds_idx)%parm_idx = diheds_copysol(diheds_idx)%parm_idx
    end do
  end if

  my_dihed_cnt_inter = 0

  if ( nphihinter .gt. 0 ) then
    do diheds_idx = 1, nphihinter

      atm_i = gbl_dihedinter(diheds_idx)%atm_i
      atm_j = gbl_dihedinter(diheds_idx)%atm_j
      atm_k = iabs(gbl_dihedinter(diheds_idx)%atm_k)
      atm_l = iabs(gbl_dihedinter(diheds_idx)%atm_l)

#if defined(MPI) && !defined(CUDA)
      if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
        my_dihed_cnt_inter = my_dihed_cnt_inter + 1
        diheds_copyinter(my_dihed_cnt_inter) = gbl_dihedinter(diheds_idx)
        !diheds_copyinter(my_dihed_cnt_inter)%atm_i = gbl_dihedinter(diheds_idx)%atm_i
        !diheds_copyinter(my_dihed_cnt_inter)%atm_j = gbl_dihedinter(diheds_idx)%atm_j
        !diheds_copyinter(my_dihed_cnt_inter)%atm_k = gbl_dihedinter(diheds_idx)%atm_k
        !diheds_copyinter(my_dihed_cnt_inter)%atm_l = gbl_dihedinter(diheds_idx)%atm_l
        !diheds_copyinter(my_dihed_cnt_inter)%parm_idx = gbl_dihedinter(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
        my_dihed_cnt_inter = my_dihed_cnt_inter + 1
        diheds_copyinter(my_dihed_cnt_inter) = gbl_dihedinter(diheds_idx)
        !diheds_copyinter(my_dihed_cnt_inter)%atm_i = gbl_dihedinter(diheds_idx)%atm_i
        !diheds_copyinter(my_dihed_cnt_inter)%atm_j = gbl_dihedinter(diheds_idx)%atm_j
        !diheds_copyinter(my_dihed_cnt_inter)%atm_k = gbl_dihedinter(diheds_idx)%atm_k
        !diheds_copyinter(my_dihed_cnt_inter)%atm_l = gbl_dihedinter(diheds_idx)%atm_l
        !diheds_copyinter(my_dihed_cnt_inter)%parm_idx = gbl_dihedinter(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
        my_dihed_cnt_inter = my_dihed_cnt_inter + 1
        diheds_copyinter(my_dihed_cnt_inter) = gbl_dihedinter(diheds_idx)
        !diheds_copyinter(my_dihed_cnt_inter)%atm_i = gbl_dihedinter(diheds_idx)%atm_i
        !diheds_copyinter(my_dihed_cnt_inter)%atm_j = gbl_dihedinter(diheds_idx)%atm_j
        !diheds_copyinter(my_dihed_cnt_inter)%atm_k = gbl_dihedinter(diheds_idx)%atm_k
        !diheds_copyinter(my_dihed_cnt_inter)%atm_l = gbl_dihedinter(diheds_idx)%atm_l
        !diheds_copyinter(my_dihed_cnt_inter)%parm_idx = gbl_dihedinter(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_l) .eq. mytaskid) then
        my_dihed_cnt_inter = my_dihed_cnt_inter + 1
        diheds_copyslu(my_dihed_cnt_inter) = gbl_dihedinter(diheds_idx)
        !diheds_copyinter(my_dihed_cnt_inter)%atm_i = gbl_dihedinter(diheds_idx)%atm_i
        !diheds_copyinter(my_dihed_cnt_inter)%atm_j = gbl_dihedinter(diheds_idx)%atm_j
        !diheds_copyinter(my_dihed_cnt_inter)%atm_k = gbl_dihedinter(diheds_idx)%atm_k
        !diheds_copyinter(my_dihed_cnt_inter)%atm_l = gbl_dihedinter(diheds_idx)%atm_l
        !diheds_copyinter(my_dihed_cnt_inter)%parm_idx = gbl_dihedinter(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      end if
#else
      my_dihed_cnt_inter = my_dihed_cnt_inter + 1
      diheds_copyinter(my_dihed_cnt_inter) = gbl_dihedinter(diheds_idx)
      !diheds_copyinter(my_dihed_cnt_inter)%atm_i = gbl_dihedinter(diheds_idx)%atm_i
      !diheds_copyinter(my_dihed_cnt_inter)%atm_j = gbl_dihedinter(diheds_idx)%atm_j
      !diheds_copyinter(my_dihed_cnt_inter)%atm_k = gbl_dihedinter(diheds_idx)%atm_k
      !diheds_copyinter(my_dihed_cnt_inter)%atm_l = gbl_dihedinter(diheds_idx)%atm_l
      !diheds_copyinter(my_dihed_cnt_inter)%parm_idx = gbl_dihedinter(diheds_idx)%parm_idx
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
#endif

    end do

    cit_nphihinter = my_dihed_cnt_inter
  endif

  if ( nphiainter .gt. 0 ) then
    do diheds_idx = diheda_idx_inter, diheda_idx_inter + nphiainter - 1

      atm_i = gbl_dihedinter(diheds_idx)%atm_i
      atm_j = gbl_dihedinter(diheds_idx)%atm_j
      atm_k = iabs(gbl_dihedinter(diheds_idx)%atm_k)
      atm_l = iabs(gbl_dihedinter(diheds_idx)%atm_l)

#if defined(MPI) && !defined(CUDA)  
      if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
        my_dihed_cnt_inter = my_dihed_cnt_inter + 1
        diheds_copyinter(my_dihed_cnt_inter) = gbl_dihedinter(diheds_idx)
        !diheds_copyinter(my_dihed_cnt_inter)%atm_i = gbl_dihedinter(diheds_idx)%atm_i
        !diheds_copyinter(my_dihed_cnt_inter)%atm_j = gbl_dihedinter(diheds_idx)%atm_j
        !diheds_copyinter(my_dihed_cnt_inter)%atm_k = gbl_dihedinter(diheds_idx)%atm_k
        !diheds_copyinter(my_dihed_cnt_inter)%atm_l = gbl_dihedinter(diheds_idx)%atm_l
        !diheds_copyinter(my_dihed_cnt_inter)%parm_idx = gbl_dihedinter(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
        my_dihed_cnt_inter = my_dihed_cnt_inter + 1
        diheds_copyinter(my_dihed_cnt_inter) = gbl_dihedinter(diheds_idx)
        !diheds_copyinter(my_dihed_cnt_inter)%atm_i = gbl_dihedinter(diheds_idx)%atm_i
        !diheds_copyinter(my_dihed_cnt_inter)%atm_j = gbl_dihedinter(diheds_idx)%atm_j
        !diheds_copyinter(my_dihed_cnt_inter)%atm_k = gbl_dihedinter(diheds_idx)%atm_k
        !diheds_copyinter(my_dihed_cnt_inter)%atm_l = gbl_dihedinter(diheds_idx)%atm_l
        !diheds_copyinter(my_dihed_cnt_inter)%parm_idx = gbl_dihedinter(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
        my_dihed_cnt_inter = my_dihed_cnt_inter + 1
        diheds_copyinter(my_dihed_cnt_inter) = gbl_dihedinter(diheds_idx)
        !diheds_copyinter(my_dihed_cnt_inter)%atm_i = gbl_dihedinter(diheds_idx)%atm_i
        !diheds_copyinter(my_dihed_cnt_inter)%atm_j = gbl_dihedinter(diheds_idx)%atm_j
        !diheds_copyinter(my_dihed_cnt_inter)%atm_k = gbl_dihedinter(diheds_idx)%atm_k
        !diheds_copyinter(my_dihed_cnt_inter)%atm_l = gbl_dihedinter(diheds_idx)%atm_l
        !diheds_copyinter(my_dihed_cnt_inter)%parm_idx = gbl_dihedinter(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      else if (gbl_atm_owner_map(atm_l) .eq. mytaskid) then
        my_dihed_cnt_inter = my_dihed_cnt_inter + 1
        diheds_copyinter(my_dihed_cnt_inter) = gbl_dihedinter(diheds_idx)
        !diheds_copyinter(my_dihed_cnt_inter)%atm_i = gbl_dihedinter(diheds_idx)%atm_i
        !diheds_copyinter(my_dihed_cnt_inter)%atm_j = gbl_dihedinter(diheds_idx)%atm_j
        !diheds_copyinter(my_dihed_cnt_inter)%atm_k = gbl_dihedinter(diheds_idx)%atm_k
        !diheds_copyinter(my_dihed_cnt_inter)%atm_l = gbl_dihedinter(diheds_idx)%atm_l
        !diheds_copyinter(my_dihed_cnt_inter)%parm_idx = gbl_dihedinter(diheds_idx)%parm_idx
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
        use_atm_map(atm_l) = 1
      end if
#else
      my_dihed_cnt_inter = my_dihed_cnt_inter + 1
      diheds_copyinter(my_dihed_cnt_inter) = gbl_dihedinter(diheds_idx)
      !diheds_copyinter(my_dihed_cnt_inter)%atm_i = gbl_dihedinter(diheds_idx)%atm_i
      !diheds_copyinter(my_dihed_cnt_inter)%atm_j = gbl_dihedinter(diheds_idx)%atm_j
      !diheds_copyinter(my_dihed_cnt_inter)%atm_k = gbl_dihedinter(diheds_idx)%atm_k
      !diheds_copyinter(my_dihed_cnt_inter)%atm_l = gbl_dihedinter(diheds_idx)%atm_l
      !diheds_copyinter(my_dihed_cnt_inter)%parm_idx = gbl_dihedinter(diheds_idx)%parm_idx
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
      use_atm_map(atm_l) = 1
#endif

    end do

    cit_nphiainter = my_dihed_cnt_inter - cit_nphihinter
  endif

  if (my_dihed_cnt_inter .gt. 0) then
    if (allocated(cit_dihedinter)) then
      if (size(cit_dihedinter) .lt. my_dihed_cnt_inter) then
        num_ints = num_ints - size(cit_dihedinter) * dihed_rec_ints
        deallocate(cit_dihedinter)
        allocate(cit_dihedinter(my_dihed_cnt_inter), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_dihedinter) * dihed_rec_ints
      end if
    else
      allocate(cit_dihedinter(my_dihed_cnt_inter), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_dihedinter) * dihed_rec_ints
    end if
    cit_dihedinter(1:my_dihed_cnt_inter) = diheds_copyinter(1:my_dihed_cnt_inter)
    !do diheds_idx = 1, my_dihed_cnt_inter 
    !   cit_dihedinter(diheds_idx)%atm_i = diheds_copyinter(diheds_idx)%atm_i
    !   cit_dihedinter(diheds_idx)%atm_j = diheds_copyinter(diheds_idx)%atm_j
    !   cit_dihedinter(diheds_idx)%atm_k = diheds_copyinter(diheds_idx)%atm_k
    !   cit_dihedinter(diheds_idx)%atm_l = diheds_copyinter(diheds_idx)%atm_l
    !   cit_dihedinter(diheds_idx)%parm_idx = diheds_copyinter(diheds_idx)%parm_idx
    !end do
  end if
#ifdef CUDA
   call gpu_dihedrals_setup(my_dihed_cnt, cit_nphih, cit_dihed, gbl_ipn, gbl_pn, gbl_pk, gbl_gamc, gbl_gams); 
#endif

!BEGIN DBG
! write(0,*)'task,cit_nphih,cit_nphia', mytaskid, &
!            cit_nphih, cit_nphia,
! write(0,*)'task,    nphih,    nphia', mytaskid, &
!                nphih,     nphia
!END DBG

  return

end subroutine dihedrals_setup

!*******************************************************************************
!
! Subroutine:  get_dihed_energy
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine get_dihed_energy(dihed_cnt, dihed, x, frc, subfrc, ep)

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer                       :: dihed_cnt
  type(dihed_rec)               :: dihed(*)
  double precision              :: x(3, *)
  double precision              :: frc(3, *)
  double precision              :: subfrc(3, *)
  double precision              :: ep

! Local variables:

  double precision      :: ap
  double precision      :: cosnp
  double precision      :: cphi
  double precision      :: ct, ct0
  double precision      :: dc1, dc2, dc3, dc4, dc5, dc6
  double precision      :: df
  double precision      :: dr1, dr2, dr3, dr4, dr5, dr6
  double precision      :: drx, dry, drz
  double precision      :: dums
  double precision      :: dx, dy, dz
  double precision      :: epl
  double precision      :: epw
  double precision      :: f1, f2
  double precision      :: fxi, fyi, fzi
  double precision      :: fxj, fyj, fzj
  double precision      :: fxk, fyk, fzk
  double precision      :: fxl, fyl, fzl
  double precision      :: g
  double precision      :: gmul(10)
  double precision      :: gx, gy, gz
  integer               :: ic, ic0
  integer               :: inc
  integer               :: i, j, k, kt, l, lt
  integer               :: jn
  double precision      :: one
  double precision      :: s
  double precision      :: sinnp
  double precision      :: sphi
  double precision      :: tenm3
  double precision      :: tm06, tm24
  double precision      :: xa, ya, za
  double precision      :: xij, yij, zij
  double precision      :: xkj, ykj, zkj
  double precision      :: xkl, ykl, zkl
  double precision      :: z1, z2
  double precision      :: z11, z12, z22
  double precision      :: zero

  data gmul /0.d0, 2.d0, 0.d0, 4.d0, 0.d0, 6.d0, 0.d0, 8.d0, 0.d0, 10.d0/

  data tm24, tm06, tenm3/1.d-18, 1.d-06, 1.d-03/

  data zero, one /0.d0, 1.d0/

! Arrays gbl_gamc = gbl_pk * cos(phase) and gbl_gams = gbl_pk * sin(phase)

  epl = zero

! Grand loop for the dihedral stuff:

  do jn = 1, dihed_cnt

    i = dihed(jn)%atm_i
    j = dihed(jn)%atm_j
    kt = dihed(jn)%atm_k
    lt = dihed(jn)%atm_l
    k = iabs(kt)
    l = iabs(lt)

! Calculation of ij, kj, kl vectors:

    xij = x(1, i) - x(1, j)
    yij = x(2, i) - x(2, j)
    zij = x(3, i) - x(3, j)
    xkj = x(1, k) - x(1, j)
    ykj = x(2, k) - x(2, j)
    zkj = x(3, k) - x(3, j)
    xkl = x(1, k) - x(1, l)
    ykl = x(2, k) - x(2, l)
    zkl = x(3, k) - x(3, l)                                   

! Get the normal vector:

    dx = yij * zkj - zij * ykj
    dy = zij * xkj - xij * zkj
    dz = xij * ykj - yij * xkj
    gx = zkj * ykl - ykj * zkl
    gy = xkj * zkl - zkj * xkl
    gz = ykj * xkl - xkj * ykl

    fxi = sqrt(dx * dx + dy * dy + dz * dz + tm24)
    fyi = sqrt(gx * gx + gy * gy + gz * gz + tm24)
    ct = dx * gx + dy * gy + dz * gz

! Branch if linear dihedral:

    if (tenm3 .le. fxi) then
      z1 = one / fxi
    else
      z1 = zero
    endif

    if (tenm3 .le. fyi) then
      z2 = one / fyi
    else
      z2 = zero
    endif

    z12 = z1 * z2

    if (z12 .ne. zero) then
      fzi = one
    else
      fzi = zero
    end if

    s = xkj * (dz * gy - dy * gz) + ykj * (dx * gz - dz * gx) + &
        zkj * (dy * gx - dx * gy)

    ap = PI - sign(acos(max(-one, min(one, ct * z12))), s)

    cphi = cos(ap)
    sphi = sin(ap)

! Calculate the energy and the derivatives with respect to cosphi:

    ic = dihed(jn)%parm_idx
    inc = gbl_ipn(ic)
    ct0 = gbl_pn(ic) * ap
    cosnp = cos(ct0)
    sinnp = sin(ct0)
    epw = (gbl_pk(ic) + cosnp * gbl_gamc(ic) + sinnp * gbl_gams(ic)) * fzi
    
    dums = sphi + sign(tm24, sphi)

    if (tm06 .gt. abs(dums)) then
      df = fzi * gbl_gamc(ic) * (gbl_pn(ic) - gmul(inc) + gmul(inc) * cphi)
    else
      df = fzi * gbl_pn(ic) * (gbl_gamc(ic) * sinnp - gbl_gams(ic) * cosnp)/dums
    end if

    if (ti_mode .ne. 0) then
      ti_region = 0
      if (ti_lst(1,i) .eq. 1 .or. ti_lst(1,j) .eq. 1 &
          .or. ti_lst(1,k) .eq. 1 .or. ti_lst(1,l) .eq. 1) then
        ti_region = 1
      else if (ti_lst(2,i) .eq. 1 .or. ti_lst(2,j) .eq. 1 &
               .or. ti_lst(2,k) .eq. 1 .or. ti_lst(2,l) .eq. 1) then
        ti_region = 2
      end if

      if (ti_region .gt. 0) then
        if (ti_mode .ne. 1) then
          ti_lscale = (ti_sc_lst(i) .lt. 2 .and. &
                       ti_sc_lst(j) .lt. 2 .and. &
                       ti_sc_lst(k) .lt. 2 .and. &
                       ti_sc_lst(l) .lt. 2) 
        end if
#ifdef MPI
        if (gbl_atm_owner_map(i) .eq. mytaskid) then
#endif     
          if (ti_mode .eq. 1 .or. ti_lscale) then !linear scaling       
            call ti_update_ene(epw, si_dihedral_ene, ti_region)
            epw = epw * ti_weights(ti_region)
          else
            ti_ene(ti_region,si_dihedral_ene) = &
               ti_ene(ti_region,si_dihedral_ene) + epw
            epw = 0.d0    
          end if
#ifdef MPI
        end if
#endif
        if (ti_mode .eq. 1 .or. ti_lscale) df = df * ti_weights(ti_region)
      end if
    end if

! Now do torsional first derivatives:

! Now, set up array dc = 1st der. of cosphi w/respect to cartesian differences:

    z11 = z1 * z1
    z12 = z1 * z2
    z22 = z2 * z2
    dc1 = -gx * z12 - cphi * dx * z11
    dc2 = -gy * z12 - cphi * dy * z11
    dc3 = -gz * z12 - cphi * dz * z11
    dc4 =  dx * z12 + cphi * gx * z22
    dc5 =  dy * z12 + cphi * gy * z22
    dc6 =  dz * z12 + cphi * gz * z22

! Update the first derivative array:

    dr1 = df * ( dc3 * ykj - dc2 * zkj)
    dr2 = df * ( dc1 * zkj - dc3 * xkj)
    dr3 = df * ( dc2 * xkj - dc1 * ykj)
    dr4 = df * ( dc6 * ykj - dc5 * zkj)
    dr5 = df * ( dc4 * zkj - dc6 * xkj)
    dr6 = df * ( dc5 * xkj - dc4 * ykj)
    drx = df * (-dc2 * zij + dc3 * yij + dc5 * zkl - dc6 * ykl)
    dry = df * ( dc1 * zij - dc3 * xij - dc4 * zkl + dc6 * xkl)
    drz = df * (-dc1 * yij + dc2 * xij + dc4 * ykl - dc5 * xkl)
    fxi = -dr1
    fyi = -dr2
    fzi = -dr3
    fxj = -drx + dr1
    fyj = -dry + dr2
    fzj = -drz + dr3
    fxk = +drx + dr4
    fyk = +dry + dr5
    fzk = +drz + dr6
    fxl = -dr4
    fyl = -dr5
    fzl = -dr6

#ifdef MPI
! We use atm_i to determine who sums up the energy...

    if (gbl_atm_owner_map(i) .eq. mytaskid) then
      epl = epl  + epw
  !    frc(1, i) = frc(1, i) + fxi
  !    frc(2, i) = frc(2, i) + fyi
  !    frc(3, i) = frc(3, i) + fzi
      subfrc(1, i) = subfrc(1, i) + fxi
      subfrc(2, i) = subfrc(2, i) + fyi
      subfrc(3, i) = subfrc(3, i) + fzi
    end if

    if (gbl_atm_owner_map(j) .eq. mytaskid) then
  !    frc(1, j) = frc(1, j) + fxj 
  !    frc(2, j) = frc(2, j) + fyj
  !    frc(3, j) = frc(3, j) + fzj
      subfrc(1, j) = subfrc(1, j) + fxj 
      subfrc(2, j) = subfrc(2, j) + fyj
      subfrc(3, j) = subfrc(3, j) + fzj
    end if

    if (gbl_atm_owner_map(k) .eq. mytaskid) then
  !    frc(1, k) = frc(1, k) + fxk
  !    frc(2, k) = frc(2, k) + fyk
  !    frc(3, k) = frc(3, k) + fzk
      subfrc(1, k) = subfrc(1, k) + fxk
      subfrc(2, k) = subfrc(2, k) + fyk
      subfrc(3, k) = subfrc(3, k) + fzk
    end if

    if (gbl_atm_owner_map(l) .eq. mytaskid) then
  !    frc(1, l) = frc(1, l) + fxl
  !    frc(2, l) = frc(2, l) + fyl
  !    frc(3, l) = frc(3, l) + fzl
      subfrc(1, l) = subfrc(1, l) + fxl
      subfrc(2, l) = subfrc(2, l) + fyl
      subfrc(3, l) = subfrc(3, l) + fzl
    end if
#else
    epl = epl + epw

 !   frc(1, i) = frc(1, i) + fxi
 !   frc(2, i) = frc(2, i) + fyi
 !   frc(3, i) = frc(3, i) + fzi
 !   frc(1, j) = frc(1, j) + fxj 
 !   frc(2, j) = frc(2, j) + fyj
 !   frc(3, j) = frc(3, j) + fzj
 !   frc(1, k) = frc(1, k) + fxk
 !   frc(2, k) = frc(2, k) + fyk
 !   frc(3, k) = frc(3, k) + fzk
 !   frc(1, l) = frc(1, l) + fxl
 !   frc(2, l) = frc(2, l) + fyl
 !   frc(3, l) = frc(3, l) + fzl
    subfrc(1, i) = subfrc(1, i) + fxi
    subfrc(2, i) = subfrc(2, i) + fyi
    subfrc(3, i) = subfrc(3, i) + fzi
    subfrc(1, j) = subfrc(1, j) + fxj 
    subfrc(2, j) = subfrc(2, j) + fyj
    subfrc(3, j) = subfrc(3, j) + fzj
    subfrc(1, k) = subfrc(1, k) + fxk
    subfrc(2, k) = subfrc(2, k) + fyk
    subfrc(3, k) = subfrc(3, k) + fzk
    subfrc(1, l) = subfrc(1, l) + fxl
    subfrc(2, l) = subfrc(2, l) + fyl
    subfrc(3, l) = subfrc(3, l) + fzl
#endif
    
  end do

  ep = epl

  return

end subroutine get_dihed_energy

!*******************************************************************************
!
! Subroutine:  get_dihed_energy_amd
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine get_dihed_energy_amd(dihed_cnt, dihed, x, ep)

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer                       :: dihed_cnt
  type(dihed_rec)               :: dihed(*)
  double precision              :: x(3, *)
  double precision              :: ep

! Local variables:

  double precision      :: ap
  double precision      :: cosnp
  double precision      :: cphi
  double precision      :: ct, ct0
  double precision      :: dx, dy, dz
  double precision      :: epl
  double precision      :: epw
  double precision      :: fxi, fyi, fzi
  double precision      :: gmul(10)
  double precision      :: gx, gy, gz
  integer               :: ic
  integer               :: inc
  integer               :: i, j, k, kt, l, lt
  integer               :: jn
  double precision      :: one
  double precision      :: s
  double precision      :: sinnp
  double precision      :: sphi
  double precision      :: tenm3
  double precision      :: tm06, tm24
  double precision      :: xij, yij, zij
  double precision      :: xkj, ykj, zkj
  double precision      :: xkl, ykl, zkl
  double precision      :: z1, z2
  double precision      :: z12
  double precision      :: zero

  data gmul /0.d0, 2.d0, 0.d0, 4.d0, 0.d0, 6.d0, 0.d0, 8.d0, 0.d0, 10.d0/

  data tm24, tm06, tenm3/1.d-18, 1.d-06, 1.d-03/

  data zero, one /0.d0, 1.d0/

! Arrays gbl_gamc = gbl_pk * cos(phase) and gbl_gams = gbl_pk * sin(phase)

  epl = zero

! Grand loop for the dihedral stuff:

  do jn = 1, dihed_cnt

    i = dihed(jn)%atm_i
    j = dihed(jn)%atm_j
    kt = dihed(jn)%atm_k
    lt = dihed(jn)%atm_l
    k = iabs(kt)
    l = iabs(lt)

! Calculation of ij, kj, kl vectors:

    xij = x(1, i) - x(1, j)
    yij = x(2, i) - x(2, j)
    zij = x(3, i) - x(3, j)
    xkj = x(1, k) - x(1, j)
    ykj = x(2, k) - x(2, j)
    zkj = x(3, k) - x(3, j)
    xkl = x(1, k) - x(1, l)
    ykl = x(2, k) - x(2, l)
    zkl = x(3, k) - x(3, l)                                   

! Get the normal vector:

    dx = yij * zkj - zij * ykj
    dy = zij * xkj - xij * zkj
    dz = xij * ykj - yij * xkj
    gx = zkj * ykl - ykj * zkl
    gy = xkj * zkl - zkj * xkl
    gz = ykj * xkl - xkj * ykl

    fxi = sqrt(dx * dx + dy * dy + dz * dz + tm24)
    fyi = sqrt(gx * gx + gy * gy + gz * gz + tm24)
    ct = dx * gx + dy * gy + dz * gz

! Branch if linear dihedral:

    if (tenm3 .le. fxi) then
      z1 = one / fxi
    else
      z1 = zero
    endif

    if (tenm3 .le. fyi) then
      z2 = one / fyi
    else
      z2 = zero
    endif

    z12 = z1 * z2

    if (z12 .ne. zero) then
      fzi = one
    else
      fzi = zero
    end if

    s = xkj * (dz * gy - dy * gz) + ykj * (dx * gz - dz * gx) + &
        zkj * (dy * gx - dx * gy)

    ap = PI - sign(acos(max(-one, min(one, ct * z12))), s)

    cphi = cos(ap)
    sphi = sin(ap)

! Calculate the energy and the derivatives with respect to cosphi:

    ic = dihed(jn)%parm_idx
    inc = gbl_ipn(ic)
    ct0 = gbl_pn(ic) * ap
    cosnp = cos(ct0)
    sinnp = sin(ct0)
    epw = (gbl_pk(ic) + cosnp * gbl_gamc(ic) + sinnp * gbl_gams(ic)) * fzi

#ifdef MPI
! We use atm_i to determine who sums up the energy...

    if (gbl_atm_owner_map(i) .eq. mytaskid) then
      epl = epl  + epw
    end if

#else
    epl = epl + epw

#endif
    
  end do

  ep = epl

  return

end subroutine get_dihed_energy_amd

!*******************************************************************************
!
! Subroutine:  get_dihed_forces_amd
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine get_dihed_forces_amd(dihed_cnt, dihed, x, frc)

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer                       :: dihed_cnt
  type(dihed_rec)               :: dihed(*)
  double precision              :: x(3, *)
  double precision              :: frc(3, *)

! Local variables:

  double precision      :: ap
  double precision      :: cosnp
  double precision      :: cphi
  double precision      :: ct, ct0
  double precision      :: dc1, dc2, dc3, dc4, dc5, dc6
  double precision      :: df
  double precision      :: dr1, dr2, dr3, dr4, dr5, dr6
  double precision      :: drx, dry, drz
  double precision      :: dums
  double precision      :: dx, dy, dz
  double precision      :: fxi, fyi, fzi
  double precision      :: fxj, fyj, fzj
  double precision      :: fxk, fyk, fzk
  double precision      :: fxl, fyl, fzl
  double precision      :: gmul(10)
  double precision      :: gx, gy, gz
  integer               :: ic
  integer               :: inc
  integer               :: i, j, k, kt, l, lt
  integer               :: jn
  double precision      :: one
  double precision      :: s
  double precision      :: sinnp
  double precision      :: sphi
  double precision      :: tenm3
  double precision      :: tm06, tm24
  double precision      :: xij, yij, zij
  double precision      :: xkj, ykj, zkj
  double precision      :: xkl, ykl, zkl
  double precision      :: z1, z2
  double precision      :: z11, z12, z22
  double precision      :: zero

  data gmul /0.d0, 2.d0, 0.d0, 4.d0, 0.d0, 6.d0, 0.d0, 8.d0, 0.d0, 10.d0/

  data tm24, tm06, tenm3/1.d-18, 1.d-06, 1.d-03/

  data zero, one /0.d0, 1.d0/

! Arrays gbl_gamc = gbl_pk * cos(phase) and gbl_gams = gbl_pk * sin(phase)

! Calculate dihedral forces using the boosting weight for amd

! Grand loop for the dihedral stuff:


  do jn = 1, dihed_cnt

    i = dihed(jn)%atm_i
    j = dihed(jn)%atm_j
    kt = dihed(jn)%atm_k
    lt = dihed(jn)%atm_l
    k = iabs(kt)
    l = iabs(lt)

! Calculation of ij, kj, kl vectors:

    xij = x(1, i) - x(1, j)
    yij = x(2, i) - x(2, j)
    zij = x(3, i) - x(3, j)
    xkj = x(1, k) - x(1, j)
    ykj = x(2, k) - x(2, j)
    zkj = x(3, k) - x(3, j)
    xkl = x(1, k) - x(1, l)
    ykl = x(2, k) - x(2, l)
    zkl = x(3, k) - x(3, l)                                   

! Get the normal vector:

    dx = yij * zkj - zij * ykj
    dy = zij * xkj - xij * zkj
    dz = xij * ykj - yij * xkj
    gx = zkj * ykl - ykj * zkl
    gy = xkj * zkl - zkj * xkl
    gz = ykj * xkl - xkj * ykl

    fxi = sqrt(dx * dx + dy * dy + dz * dz + tm24)
    fyi = sqrt(gx * gx + gy * gy + gz * gz + tm24)
    ct = dx * gx + dy * gy + dz * gz

! Branch if linear dihedral:

    if (tenm3 .le. fxi) then
      z1 = one / fxi
    else
      z1 = zero
    endif

    if (tenm3 .le. fyi) then
      z2 = one / fyi
    else
      z2 = zero
    endif

    z12 = z1 * z2

    if (z12 .ne. zero) then
      fzi = one
    else
      fzi = zero
    end if

    s = xkj * (dz * gy - dy * gz) + ykj * (dx * gz - dz * gx) + &
        zkj * (dy * gx - dx * gy)

    ap = PI - sign(acos(max(-one, min(one, ct * z12))), s)

    cphi = cos(ap)
    sphi = sin(ap)

! Calculate the energy and the derivatives with respect to cosphi:

    ic = dihed(jn)%parm_idx
    inc = gbl_ipn(ic)
    ct0 = gbl_pn(ic) * ap
    cosnp = cos(ct0)
    sinnp = sin(ct0)

    dums = sphi + sign(tm24, sphi)

    if (tm06 .gt. abs(dums)) then
      df = fzi * gbl_gamc(ic) * (gbl_pn(ic) - gmul(inc) + gmul(inc) * cphi)
    else
      df = fzi * gbl_pn(ic) * (gbl_gamc(ic) * sinnp - gbl_gams(ic) * cosnp)/dums
    end if

! Now do torsional first derivatives:

! Now, set up array dc = 1st der. of cosphi w/respect to cartesian differences:

    z11 = z1 * z1
    z12 = z1 * z2
    z22 = z2 * z2
    dc1 = -gx * z12 - cphi * dx * z11
    dc2 = -gy * z12 - cphi * dy * z11
    dc3 = -gz * z12 - cphi * dz * z11
    dc4 =  dx * z12 + cphi * gx * z22
    dc5 =  dy * z12 + cphi * gy * z22
    dc6 =  dz * z12 + cphi * gz * z22

! Update the first derivative array:

    df = df * fwgtd
    dr1 = df * ( dc3 * ykj - dc2 * zkj)
    dr2 = df * ( dc1 * zkj - dc3 * xkj)
    dr3 = df * ( dc2 * xkj - dc1 * ykj)
    dr4 = df * ( dc6 * ykj - dc5 * zkj)
    dr5 = df * ( dc4 * zkj - dc6 * xkj)
    dr6 = df * ( dc5 * xkj - dc4 * ykj)
    drx = df * (-dc2 * zij + dc3 * yij + dc5 * zkl - dc6 * ykl)
    dry = df * ( dc1 * zij - dc3 * xij - dc4 * zkl + dc6 * xkl)
    drz = df * (-dc1 * yij + dc2 * xij + dc4 * ykl - dc5 * xkl)
    fxi = -dr1
    fyi = -dr2
    fzi = -dr3
    fxj = -drx + dr1
    fyj = -dry + dr2
    fzj = -drz + dr3
    fxk = +drx + dr4
    fyk = +dry + dr5
    fzk = +drz + dr6
    fxl = -dr4
    fyl = -dr5
    fzl = -dr6

#ifdef MPI
! We use atm_i to determine who sums up the energy...

    if (gbl_atm_owner_map(i) .eq. mytaskid) then
      frc(1, i) = frc(1, i) + fxi
      frc(2, i) = frc(2, i) + fyi
      frc(3, i) = frc(3, i) + fzi
    end if

    if (gbl_atm_owner_map(j) .eq. mytaskid) then
      frc(1, j) = frc(1, j) + fxj 
      frc(2, j) = frc(2, j) + fyj
      frc(3, j) = frc(3, j) + fzj
    end if

    if (gbl_atm_owner_map(k) .eq. mytaskid) then
      frc(1, k) = frc(1, k) + fxk
      frc(2, k) = frc(2, k) + fyk
      frc(3, k) = frc(3, k) + fzk
    end if

    if (gbl_atm_owner_map(l) .eq. mytaskid) then
      frc(1, l) = frc(1, l) + fxl
      frc(2, l) = frc(2, l) + fyl
      frc(3, l) = frc(3, l) + fzl
    end if
#else

    frc(1, i) = frc(1, i) + fxi
    frc(2, i) = frc(2, i) + fyi
    frc(3, i) = frc(3, i) + fzi
    frc(1, j) = frc(1, j) + fxj 
    frc(2, j) = frc(2, j) + fyj
    frc(3, j) = frc(3, j) + fzj
    frc(1, k) = frc(1, k) + fxk
    frc(2, k) = frc(2, k) + fyk
    frc(3, k) = frc(3, k) + fzk
    frc(1, l) = frc(1, l) + fxl
    frc(2, l) = frc(2, l) + fyl
    frc(3, l) = frc(3, l) + fzl
#endif
    
  end do

  return

end subroutine get_dihed_forces_amd

end module dihedrals_mod
