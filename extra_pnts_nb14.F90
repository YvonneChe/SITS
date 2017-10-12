#include "copyright.i"

!*******************************************************************************!
! Module:  extra_pnts_nb14_mod
!
! Description: <TBS>
!
!*******************************************************************************

module extra_pnts_nb14_mod

  use gbl_datatypes_mod

  implicit none

! Frames are centered on an atom parent_atm that "owns" extra points
! ep_cnt are number of extra points attached to parent_atm
! extra_pnt are pointers to (at most 2) extra points attached to parent_atm
! type is 1 if atom has at least 2 other bonds, pref to heavy atoms
! type is 2 if atom has only one other bond e.g. O6 in guanine
! frame_atm1 frame_atm2 and frame_atm3 are atom nums defining frame
! loc_frame are local coords

! Original code by Tom Darden based on Jim Caldwell's ideas.
!
! All atom names beginning with EP are considered extra points.
!
! If you want bond angle and dihedral forces as usual for lone pairs
! (i.e. as if they are amber atoms) then set frameon = 0.
!
! If frameon .eq. 1, (DEFAULT) the bonds, angles and dihedral interactions
! involving the lone pairs / extra points are removed except for constraints
! added during parm. The lone pairs are kept in ideal geometry relative to
! local atoms, and resulting torques are transferred to these atoms.
!
! If chngmask .eq. 1 (DEFAULT), new 1-1, 1-2, 1-3 and 1-4 interactions are
! calculated. An extra point belonging to an atom has a 1-1 interaction with it,
! and participates in any 1-2, 1-3 or 1-4 interaction that atom has.
! For example, suppose (excusing the geometry)
! C1,C2,C3,C4 form a dihedral and each has 1 extra point attached as below
!
!           C1------C2------C3---------C4
!           |        |       |         |
!           |        |       |         |
!          Ep1      Ep2     Ep3       Ep4
!
! The 1-4 interactions include  C1&C4, Ep1&C4, C1&Ep4, and Ep1&Ep4
!
! To see a printout of all 1-1, 1-2, 1-3 and 1-4 interactions set verbose = 1.
! These interactions are masked out of nonbonds. Thus the amber mask list is
! rebuilt from these 1-1, 1-2, 1-3 and 1-4 pairs. Pairs that aren't in the
! union of these are not masked.
!
! A separate list of 1-4 nonbonds is then compiled. This list does not agree
! in general with the above 1-4, since a 1-4 could also be a 1-3 if its
! in a ring. The rules in ephi() are used to see who is included:
!
! Here is that code:
!
!             DO 700 JN = 1,MAXLEN
!               I3 = IP(JN + IST)
!               K3T = KP(JN + IST)
!               L3T = LP(JN + IST)
!               IC0 = ICP(JN + IST)
!               IDUMI = ISIGN(1,K3T)
!               IDUML = ISIGN(1,L3T)
!               KDIV = (2 + IDUMI + IDUML) / 4
!               L3 = IABS(L3T)
!               FMULN = FLOAT(KDIV) * FMN(IC0)
!   C
!               II = (I3 + 3) / 3
!               JJ = (L3 + 3) / 3
!               IA1 = IAC(II)
!               IA2 = IAC(JJ)
!               IBIG = MAX0(IA1,IA2)
!               ISML = MIN0(IA1,IA2)
!               IC = IBIG * (IBIG - 1) / 2 + ISML
!   C
!   C             ----- CALCULATE THE 14-EEL ENERGY -----
!   C
!               R2 = FMULN / CT(JN)
!               R1 = SQRT(R2)
!       ...........
!
! So a pair is included in the 1-4 list if kdiv .gt. 0 and FMN(ic0) .gt. 0.
! This is decided at startup. This decision logic is applied to the parent
! atoms, and if they are included, so are the extra points attached:
!
! That is, in the above situation, if C1 and C4 pass the test then
! C1&C4, Ep1&C4, C1&Ep4, and Ep1&Ep4 are included. The dihedrals involving the
! extra points are not tested since the decision is based solely on parent
! atoms.
!
! The list of 1-4 nonbonds is also spit out if verbose = 1.
!
! To scale 1-4 charge-dipole and dipole-dipole interactions the same as
! 1-4 charge-charge (i.e. divided by scee) set scaldip = 1 (DEFAULT).
! If scaldip .eq. 0 the 1-4 charge-dipole and dipole-dipole interactions
! are treated the same as other dipolar interactions (i.e. divided by 1)
!
!-----------------------------------------------------------------------------
! Using the bond list (the one not involving hydrogens), find the number
! of neighbors (heavy atom, hydrogens and extra points) attached to each
! atom. For example if atom i is heavy, numnghbr(1,i) is the number of heavy
! atoms attached to atom i, while numnghbr(2,i) is the number of hydrogens
! attached, and numnghbr(3,i) is the number of extra points attached to i.
! The identities of neighbors are packed back to back in nghbrlst. The attyp
! array holds the atom types, usded to distinguish extra points or lone pairs
! from regular atoms.
!-----------------------------------------------------------------------------

  integer, parameter    :: extra_pnts_nb14_int_cnt = 5

  integer                     gbl_nb14_cnt, gbl_frame_cnt
! Added by Lijiang Yang 
  integer                     gbl_nb14_cnt_slu, gbl_nb14_cnt_sol, gbl_nb14_cnt_inter
!-----------------------------------------------------------------------------------

! Modified by Lijiang Yang 
  common / extra_pnts_nb14_int / gbl_nb14_cnt, gbl_frame_cnt, gbl_nb14_cnt_slu, gbl_nb14_cnt_sol, gbl_nb14_cnt_inter
!-----------------------------------------------------------------------------------

  ! Stuff below here used for 1-4 nonbonded interactions:

  integer, save                                 :: cit_nb14_cnt = 0
! Added by Lijiang Yang 
  integer, save                                 :: cit_nb14_cnt_slu = 0
  integer, save                                 :: cit_nb14_cnt_sol = 0
  integer, save                                 :: cit_nb14_cnt_inter = 0
!-----------------------------------------------------------------------------------

  integer,              allocatable, save       :: gbl_nb14(:,:)
! Added by Lijiang Yang 
  integer,              allocatable, save       :: gbl_nb14slu(:,:)
  integer,              allocatable, save       :: gbl_nb14sol(:,:)
  integer,              allocatable, save       :: gbl_nb14inter(:,:)
!-----------------------------------------------------------------------------------

  ! This is allocated / reallocated as needed in nb14_setup():

  integer,              allocatable, save       :: cit_nb14(:,:)
! Added by Lijiang Yang 
  integer,              allocatable, save       :: cit_nb14slu(:,:)
  integer,              allocatable, save       :: cit_nb14sol(:,:)
  integer,              allocatable, save       :: cit_nb14inter(:,:)
!-----------------------------------------------------------------------------------

  ! Nonbonded 1-4 interactions possibly increase by:
  integer, parameter                            :: nb14_mult_fac = 9

  ! Extra Points Frame datatype:

  type ep_frame_rec
    sequence
    integer             :: extra_pnt(2)
    integer             :: ep_cnt
    integer             :: type
    integer             :: parent_atm
    integer             :: frame_atm1
    integer             :: frame_atm2
    integer             :: frame_atm3
  end type ep_frame_rec

  integer, parameter    :: ep_frame_rec_ints = 8 ! don't use for allocation!

  type(ep_frame_rec), parameter :: null_ep_frame_rec = &
                                   ep_frame_rec(2*0,0,0,0,0,0,0)

  ! Stuff below here used only for extra points:

  type(ep_frame_rec),   allocatable, save       :: ep_frames(:)
  double precision,     allocatable, save       :: ep_lcl_crd(:,:,:)
#ifdef MPI
  integer, save                                 :: my_ep_frame_cnt
  integer, allocatable, save                    :: gbl_my_ep_frame_lst(:)
#endif /* MPI */

  character(len=4), parameter                   :: ep_symbl = 'EP  '
  integer, private                              :: maxa

  private       get_nghbrs, &
                define_frames, &
                fill_bonded, &
                redo_masked, &
                build_nb14, &
                copy_nb14, &
                trim_bonds, &
                trim_theta, &
                trim_phi, &
                do_bond_pairs, &
                do_angle_pairs, &
                do_dihed_pairs, &
                add_one_list_iblo, &
                add_one_list_inb, &
                do_14pairs

contains

!*******************************************************************************
!
! Subroutine:  nb14_setup
!
! Description:  Handle workload subdivision for nb14 list.  This is also called
!               once in uniprocessor code, but this is really not necessary -
!               just the convention borrowed from bond-angle-dihedral code. One
!               advantage, though, is that gbl_nb14 can then be deallocated in
!               uniprocessor code, and it is probably bigger than it needs to
!               be.
!*******************************************************************************

subroutine nb14_setup(num_ints, num_reals, use_atm_map)

  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
#ifdef CUDA
  use charmm_mod, only : charmm_active
#endif

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals
  integer                       :: use_atm_map(natom)

! Local variables:

  integer               :: alloc_failed
!  integer               :: nb14_copy(3, gbl_nb14_cnt)
  integer               :: nb14_copyslu(3, gbl_nb14_cnt_slu)
  integer               :: nb14_copysol(3, gbl_nb14_cnt_sol)
  integer               :: nb14_copyinter(3, gbl_nb14_cnt_inter)
  integer               :: atm_i, atm_j, nb14_idx

! This routine can handle reallocation, and thus can be called multiple
! times.

! Find all diheds for which this process owns either atom:

  cit_nb14_cnt_slu = 0

  do nb14_idx = 1, gbl_nb14_cnt_slu

    atm_i = gbl_nb14slu(1, nb14_idx)
    atm_j = gbl_nb14slu(2, nb14_idx)
!    nb14_parm_idx = gbl_nb14(3, nb14_idx)

#if defined(MPI) && !defined(CUDA)
    if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
      cit_nb14_cnt_slu = cit_nb14_cnt_slu + 1
      nb14_copyslu(:, cit_nb14_cnt_slu) = gbl_nb14slu(:, nb14_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
      cit_nb14_cnt_slu = cit_nb14_cnt_slu + 1
      nb14_copyslu(:, cit_nb14_cnt_slu) = gbl_nb14slu(:, nb14_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    end if
#else
    cit_nb14_cnt_slu = cit_nb14_cnt_slu + 1
    nb14_copyslu(:, cit_nb14_cnt_slu) = gbl_nb14slu(:, nb14_idx)
    use_atm_map(atm_i) = 1
    use_atm_map(atm_j) = 1
#endif

  end do

  cit_nb14_cnt_sol = 0

  do nb14_idx = 1, gbl_nb14_cnt_sol

    atm_i = gbl_nb14sol(1, nb14_idx)
    atm_j = gbl_nb14sol(2, nb14_idx)
!    nb14_parm_idx = gbl_nb14(3, nb14_idx)

#if defined(MPI) && !defined(CUDA)
    if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
      cit_nb14_cnt_sol = cit_nb14_cnt_sol + 1
      nb14_copysol(:, cit_nb14_cnt_sol) = gbl_nb14sol(:, nb14_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
      cit_nb14_cnt_sol = cit_nb14_cnt_sol + 1
      nb14_copysol(:, cit_nb14_cnt_sol) = gbl_nb14sol(:, nb14_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    end if
#else
    cit_nb14_cnt_sol = cit_nb14_cnt_sol + 1
    nb14_copysol(:, cit_nb14_cnt_sol) = gbl_nb14sol(:, nb14_idx)
    use_atm_map(atm_i) = 1
    use_atm_map(atm_j) = 1
#endif

  end do

  cit_nb14_cnt_inter = 0

  do nb14_idx = 1, gbl_nb14_cnt_inter

    atm_i = gbl_nb14inter(1, nb14_idx)
    atm_j = gbl_nb14inter(2, nb14_idx)
!    nb14_parm_idx = gbl_nb14(3, nb14_idx)

#if defined(MPI) && !defined(CUDA)
    if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
      cit_nb14_cnt_inter = cit_nb14_cnt_inter + 1
      nb14_copyinter(:, cit_nb14_cnt_inter) = gbl_nb14inter(:, nb14_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
      cit_nb14_cnt_inter = cit_nb14_cnt_inter + 1
      nb14_copyinter(:, cit_nb14_cnt_inter) = gbl_nb14inter(:, nb14_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
    end if
#else
    cit_nb14_cnt_inter = cit_nb14_cnt_inter + 1
    nb14_copyinter(:, cit_nb14_cnt_inter) = gbl_nb14inter(:, nb14_idx)
    use_atm_map(atm_i) = 1
    use_atm_map(atm_j) = 1
#endif

  end do

  if (cit_nb14_cnt_slu .gt. 0) then
    if (allocated(cit_nb14slu)) then
      if (size(cit_nb14slu) .lt. cit_nb14_cnt_slu * 3) then
        num_ints = num_ints - size(cit_nb14slu) * 3
        deallocate(cit_nb14slu)
        allocate(cit_nb14slu(3, cit_nb14_cnt_slu), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_nb14slu) * 3
      end if
    else
      allocate(cit_nb14slu(3, cit_nb14_cnt_slu), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_nb14slu) * 3
    end if
    cit_nb14slu(:, 1:cit_nb14_cnt_slu) = nb14_copyslu(:, 1:cit_nb14_cnt_slu)
  end if

  if (cit_nb14_cnt_sol .gt. 0) then
    if (allocated(cit_nb14sol)) then
      if (size(cit_nb14sol) .lt. cit_nb14_cnt_sol * 3) then
        num_ints = num_ints - size(cit_nb14sol) * 3
        deallocate(cit_nb14sol)
        allocate(cit_nb14sol(3, cit_nb14_cnt_sol), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_nb14sol) * 3
      end if
    else
      allocate(cit_nb14sol(3, cit_nb14_cnt_sol), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_nb14sol) * 3
    end if
    cit_nb14sol(:, 1:cit_nb14_cnt_sol) = nb14_copysol(:, 1:cit_nb14_cnt_sol)
  end if

  if (cit_nb14_cnt_inter .gt. 0) then
    if (allocated(cit_nb14inter)) then
      if (size(cit_nb14inter) .lt. cit_nb14_cnt_inter * 3) then
        num_ints = num_ints - size(cit_nb14inter) * 3
        deallocate(cit_nb14inter)
        allocate(cit_nb14inter(3, cit_nb14_cnt_inter), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_nb14inter) * 3
      end if
    else
      allocate(cit_nb14inter(3, cit_nb14_cnt_inter), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_nb14inter) * 3
    end if
    cit_nb14inter(:, 1:cit_nb14_cnt_inter) = nb14_copyinter(:, 1:cit_nb14_cnt_inter)
  end if

#ifdef CUDA
  if (charmm_active) then
    call gpu_nb14_setup(cit_nb14_cnt, cit_nb14, gbl_one_scee, gbl_one_scnb, ntypes, atm_iac, typ_ico, gbl_cn1, gbl_cn2, &
      gbl_cn114, gbl_cn214)
  else
    call gpu_nb14_setup(cit_nb14_cnt, cit_nb14, gbl_one_scee, gbl_one_scnb, ntypes, atm_iac, typ_ico, gbl_cn1, gbl_cn2, &
      gbl_cn1, gbl_cn2)
  end if
#endif 
  
  

!BEGIN DBG
! write(0,*)'task, cit_nb14_cnt =', mytaskid, cit_nb14_cnt
!END DBG

  return

end subroutine nb14_setup

!*******************************************************************************!
! Subroutine:  alloc_nb14_mem_only
!
! Description: <TBS>
!
!*******************************************************************************

subroutine alloc_nb14_mem_only(max_nb14, num_ints, num_reals)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: max_nb14

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed

! Modified by Lijiang Yang
  allocate(gbl_nb14(3, max_nb14), &
           gbl_nb14slu(3, max_nb14), &
           gbl_nb14sol(3, max_nb14), &
           gbl_nb14inter(3, max_nb14), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(gbl_nb14) * 4

  gbl_nb14(:,:) = 0
  gbl_nb14slu(:,:) = 0
  gbl_nb14sol(:,:) = 0
  gbl_nb14inter(:,:) = 0
!--------------------------------------------- 
  return

end subroutine alloc_nb14_mem_only 

!*******************************************************************************!
! Subroutine:  alloc_extra_pnts_nb14_mem
!
! Description: <TBS>
!
!*******************************************************************************

subroutine alloc_extra_pnts_nb14_mem(max_frames, max_nb14, num_ints, num_reals)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: max_frames
  integer, intent(in)           :: max_nb14

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer               :: alloc_failed

  allocate(ep_frames(max_frames), &
           ep_lcl_crd(3, 2, max_frames), &
           gbl_nb14(3, max_nb14), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(ep_frames) * ep_frame_rec_ints + &
                        size(gbl_nb14)

  num_reals = num_reals + size(ep_lcl_crd)

  ep_frames(:) = null_ep_frame_rec
  gbl_nb14(:,:) = 0
  ep_lcl_crd(:,:,:) = 0.d0
  
  return

end subroutine alloc_extra_pnts_nb14_mem 

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_extra_pnts_nb14_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_extra_pnts_nb14_dat

  use parallel_dat_mod
  use pmemd_lib_mod
  use prmtop_dat_mod

  implicit none

! Local variables:

  integer       :: num_ints, num_reals  ! returned values discarded
  integer       :: bytes_per_unit
  integer       :: alloc_failed

  call mpi_bcast(gbl_nb14_cnt, extra_pnts_nb14_int_cnt, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  if (.not. master) then
    num_ints = 0
    num_reals = 0
    if (gbl_frame_cnt .eq. 0) then
      call alloc_nb14_mem_only(gbl_nb14_cnt, num_ints, num_reals)
    else
      call alloc_extra_pnts_nb14_mem(numextra, gbl_nb14_cnt, &
                                     num_ints, num_reals)
    end if
  end if

  ! All callers need the nb14 stuff...

! Modified by Lijiang Yang
  call mpi_bcast(gbl_nb14, gbl_nb14_cnt * 3, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)
  call mpi_bcast(gbl_nb14slu, gbl_nb14_cnt_slu * 3, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)
  call mpi_bcast(gbl_nb14sol, gbl_nb14_cnt_sol * 3, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)
  call mpi_bcast(gbl_nb14inter, gbl_nb14_cnt_inter * 3, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)


  if (gbl_frame_cnt .ne. 0) then

    call get_bytesize(ep_frames(1), ep_frames(2), bytes_per_unit)

    call mpi_bcast(ep_frames, gbl_frame_cnt * bytes_per_unit, mpi_byte, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(ep_lcl_crd, 3 * 2 * gbl_frame_cnt, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    
    ! Here we allocate the process-local frame list in all tasks, including
    ! the master:

    allocate(gbl_my_ep_frame_lst(gbl_frame_cnt), stat = alloc_failed)
    if (alloc_failed .ne. 0) call setup_alloc_error
    num_ints = num_ints + size(gbl_my_ep_frame_lst)
    gbl_my_ep_frame_lst(:) = 0
    my_ep_frame_cnt = 0 ! until first atom division...

  end if

  return

end subroutine bcast_extra_pnts_nb14_dat
#endif

!*******************************************************************************!
! Subroutine:  init_nb14_only
!
! Description: <TBS>
!
!*******************************************************************************

subroutine init_nb14_only(num_ints, num_reals)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables

  integer                       :: dihed_idx

  call alloc_nb14_mem_only(2 * (nphih + nphia), num_ints, num_reals)
   
  gbl_nb14_cnt = 0
  gbl_frame_cnt = 0

  do dihed_idx = 1, nphih + nphia

    if (gbl_dihed(dihed_idx)%atm_k .gt. 0 .and. &
        gbl_dihed(dihed_idx)%atm_l .gt. 0 .and. &
        gbl_fmn(gbl_dihed(dihed_idx)%parm_idx) .ne. 0.d0) then
        gbl_nb14_cnt = gbl_nb14_cnt + 1
        gbl_nb14(1, gbl_nb14_cnt) = gbl_dihed(dihed_idx)%atm_i
        gbl_nb14(2, gbl_nb14_cnt) = gbl_dihed(dihed_idx)%atm_l
        gbl_nb14(3, gbl_nb14_cnt) = gbl_dihed(dihed_idx)%parm_idx
    end if
  end do  !  n = 1, nphih + nphia
! Modified by Lijiang Yang
  gbl_nb14_cnt_slu = 0
  do dihed_idx = 1, nphihslu + nphiaslu

    if (gbl_dihedslu(dihed_idx)%atm_k .gt. 0 .and. &
        gbl_dihedslu(dihed_idx)%atm_l .gt. 0 .and. &
        gbl_fmn(gbl_dihedslu(dihed_idx)%parm_idx) .ne. 0.d0) then
        gbl_nb14_cnt_slu = gbl_nb14_cnt_slu + 1
        gbl_nb14slu(1, gbl_nb14_cnt_slu) = gbl_dihedslu(dihed_idx)%atm_i
        gbl_nb14slu(2, gbl_nb14_cnt_slu) = gbl_dihedslu(dihed_idx)%atm_l
        gbl_nb14slu(3, gbl_nb14_cnt_slu) = gbl_dihedslu(dihed_idx)%parm_idx
    end if
  end do  !  n = 1, nphihslu + nphiaslu

  gbl_nb14_cnt_sol = 0
  do dihed_idx = 1, nphihsol + nphiasol

    if (gbl_dihedsol(dihed_idx)%atm_k .gt. 0 .and. &
        gbl_dihedsol(dihed_idx)%atm_l .gt. 0 .and. &
        gbl_fmn(gbl_dihedsol(dihed_idx)%parm_idx) .ne. 0.d0) then
        gbl_nb14_cnt_sol = gbl_nb14_cnt_sol + 1
        gbl_nb14sol(1, gbl_nb14_cnt_sol) = gbl_dihedsol(dihed_idx)%atm_i
        gbl_nb14sol(2, gbl_nb14_cnt_sol) = gbl_dihedsol(dihed_idx)%atm_l
        gbl_nb14sol(3, gbl_nb14_cnt_sol) = gbl_dihedsol(dihed_idx)%parm_idx
    end if
  end do  !  n = 1, nphihsol + nphiasol

  gbl_nb14_cnt_inter = 0
  do dihed_idx = 1, nphihinter + nphiainter

    if (gbl_dihedinter(dihed_idx)%atm_k .gt. 0 .and. &
        gbl_dihedinter(dihed_idx)%atm_l .gt. 0 .and. &
        gbl_fmn(gbl_dihedinter(dihed_idx)%parm_idx) .ne. 0.d0) then
        gbl_nb14_cnt_inter = gbl_nb14_cnt_inter + 1
        gbl_nb14inter(1, gbl_nb14_cnt_inter) = gbl_dihedinter(dihed_idx)%atm_i
        gbl_nb14inter(2, gbl_nb14_cnt_inter) = gbl_dihedinter(dihed_idx)%atm_l
        gbl_nb14inter(3, gbl_nb14_cnt_inter) = gbl_dihedinter(dihed_idx)%parm_idx
    end if
  end do  !  n = 1, nphihinter + nphiainter

  return

end subroutine init_nb14_only 

!*******************************************************************************!
! Subroutine:  init_extra_pnts_nb14
!
! Description: <TBS>
!
!*******************************************************************************

subroutine init_extra_pnts_nb14(num_ints, num_reals)

  use mdin_ewald_dat_mod
  use prmtop_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer, allocatable  :: epbtyp(:,:)
  integer, allocatable  :: nghbrs(:,:)
  integer, allocatable  :: hnghbrs(:,:)
  integer, allocatable  :: enghbrs(:,:)
  integer, allocatable  :: numnghbr(:,:)
  integer, allocatable  :: epowner(:)
  integer, allocatable  :: offset(:)
  integer, allocatable  :: test(:)
  integer, allocatable  :: i11(:,:)
  integer, allocatable  :: i12(:,:)
  integer, allocatable  :: i13(:,:)
  integer, allocatable  :: i14(:,:)
  integer, allocatable  :: nb_14_list(:,:)
  integer, allocatable  :: s3(:)
   
  integer               :: n
  integer               :: numextra_test
  integer               :: num11, num12, num13, num14
  integer               :: max11, max12, max13, max14
  integer               :: alloc_failed
  integer               :: i
   
  gbl_nb14_cnt = 0
  gbl_frame_cnt = 0

  call alloc_extra_pnts_nb14_mem(numextra, nb14_mult_fac * (nphih + nphia), &
                                 num_ints, num_reals)
   
  numextra_test = 0
  do n = 1, natom
    if (atm_isymbl(n) .eq. ep_symbl) numextra_test = numextra_test + 1
  end do

  if (numextra_test .ne. numextra) then
    write(6, *) 'Error in numextra_test'
    call mexit(6, 1)
  end if

  if (numextra .eq. 0) frameon = 0

  max11 = natom + numextra
  max12 = 3 * (nbonh + nbona)
  max13 = 3 * (ntheth + ntheta)
  max14 = nb14_mult_fac * (nphih + nphia)
  maxa = max(max11, max12, max13, max14)

  allocate(epbtyp(5, natom), &
           nghbrs(5, natom), &
           hnghbrs(5, natom), &
           enghbrs(5, natom), &
           numnghbr(3, natom), &
           epowner(natom), &
           offset(natom), &
           test(natom), &
           i11(2, (natom + numextra)), &
           i12(2, 3 * (nbonh + nbona)), &
           i13(2, 3 * (ntheth + ntheta)), &
           i14(2, nb14_mult_fac * (nphih + nphia)), &
           nb_14_list(3,nb14_mult_fac*(nphih + nphia)), &
           s3(maxa), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  nghbrs(:,:) = 0
  hnghbrs(:,:) = 0
  enghbrs(:,:) = 0
  numnghbr(:,:) = 0
  epowner(:) = 0
  epbtyp(:,:) = 0

  i11(:,:) = 0
  i12(:,:) = 0
  i13(:,:) = 0
  i14(:,:) = 0
  nb_14_list(:,:) = 0

  call get_nghbrs(nghbrs, hnghbrs, enghbrs, numnghbr, epowner, epbtyp)

  call define_frames(natom, atm_isymbl, nghbrs, hnghbrs, enghbrs, &
                     numnghbr, ep_frames, ep_lcl_crd, &
                     gbl_frame_cnt, epbtyp, gbl_req, verbose)
   
  call fill_bonded(max11, max12, max13, max14, num11, num12, num13, num14, &
                   i11, i12, i13, i14, enghbrs, numnghbr, epowner, verbose)
   
  if (chngmask .eq. 1 ) then
    call redo_masked(natom, atm_numex, gbl_natex, next, &
                     num11, num12, num13, num14, &
                     i11, i12, i13, i14, offset, test)
  end if
   
  ! Use nb_14_list and then copy to permanent:
   
  call build_nb14(nb_14_list, gbl_nb14_cnt, max14, epowner, numnghbr, enghbrs, &
                  natom, chngmask, verbose)

  call copy_nb14(nb_14_list, gbl_nb14, gbl_nb14_cnt)
   
  ! if frameon = 1 use frames and forces to define ep's;
  ! else use their masses and amber force params
   
  if (frameon .eq. 1) then
      
    ! Zero out mass and massinv for extra points:
      
    call fix_masses(natom, atm_mass, epowner)
      
    ! Now remove bonds etc involving extra points:
      
    call trim_bonds(nbonh, gbl_bond, epowner)
    call trim_bonds(nbona, gbl_bond(bonda_idx), epowner)

    ! Make the bond arrays sequential for shake and force routines:

    do i = 1, nbona
      gbl_bond(nbonh + i) = gbl_bond(bonda_idx + i - 1)
    end do

    bonda_idx = nbonh + 1

    call trim_theta(ntheth, gbl_angle, epowner)
    call trim_theta(ntheta, gbl_angle(anglea_idx), epowner)

    ! Make the angle arrays sequential:

    do i = 1, ntheta
      gbl_angle(ntheth + i) = gbl_angle(anglea_idx + i - 1)
    end do

    anglea_idx = ntheth + 1

    call trim_phi(nphih, gbl_dihed, epowner)
    call trim_phi(nphia, gbl_dihed(diheda_idx), epowner)

    do i = 1, nphia
      gbl_dihed(nphih + i) = gbl_dihed(diheda_idx + i - 1)
    end do

    diheda_idx = nphih + 1

  end if
   
  deallocate(epbtyp, nghbrs, hnghbrs, enghbrs, numnghbr, epowner, offset, &
             test, i11, i12, i13, i14, nb_14_list, s3)

  return

end subroutine init_extra_pnts_nb14 

!*******************************************************************************!
! Internal Subroutine:  get_nghbrs
!
! Description: find atoms involved in each center with an EP
!
!*******************************************************************************

subroutine get_nghbrs(nghbrs, hnghbrs, enghbrs, numnghbr, epowner, epbtyp)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: nghbrs(5, *)
  integer               :: hnghbrs(5, *)
  integer               :: enghbrs(5, *)
  integer               :: numnghbr(3, *)
  integer               :: epowner(*)
  integer               :: epbtyp(5, *)

! Local variables:

  integer               :: bond_idx, atm_i, atm_j

  do bond_idx = bonda_idx, bonda_idx + nbona - 1

    atm_i = gbl_bond(bond_idx)%atm_i
    atm_j = gbl_bond(bond_idx)%atm_j

    if (atm_isymbl(atm_i) .eq. ep_symbl) then
      numnghbr(3, atm_j) = numnghbr(3, atm_j) + 1
      enghbrs(numnghbr(3, atm_j), atm_j) = atm_i
      epowner(atm_i) = atm_j
      epbtyp(numnghbr(3, atm_j), atm_j) = gbl_bond(bond_idx)%parm_idx
    else if (atm_isymbl(atm_j) .eq. ep_symbl) then
      numnghbr(3, atm_i) = numnghbr(3, atm_i) + 1
      enghbrs(numnghbr(3, atm_i), atm_i) = atm_j
      epowner(atm_j) = atm_i
      epbtyp(numnghbr(3, atm_i), atm_i) = gbl_bond(bond_idx)%parm_idx
    else
      numnghbr(1, atm_i) = numnghbr(1, atm_i) + 1
      numnghbr(1, atm_j) = numnghbr(1, atm_j) + 1
      nghbrs(numnghbr(1, atm_i), atm_i) = atm_j
      nghbrs(numnghbr(1, atm_j), atm_j) = atm_i
    end if

  end do

  do bond_idx = 1, nbonh

    atm_i = gbl_bond(bond_idx)%atm_i
    atm_j = gbl_bond(bond_idx)%atm_j

    numnghbr(2, atm_i) = numnghbr(2, atm_i) + 1
    numnghbr(2, atm_j) = numnghbr(2, atm_j) + 1
    hnghbrs(numnghbr(2, atm_i), atm_i) = atm_j
    hnghbrs(numnghbr(2, atm_j), atm_j) = atm_i

  end do

  return

end subroutine get_nghbrs 

!*******************************************************************************!
! Internal Subroutine:  define_frames
!
! Description:  Fix the positions of EP in the local frame / coord depending
!               on the kind of frame and atom types.
!
!*******************************************************************************

subroutine define_frames(natom, isymbl, nghbrs, hnghbrs, enghbrs, numnghbr, &
                         frames, ep_lcl_crd, frame_cnt, epbtyp, req, verbose)

  use gbl_constants_mod, only : DEG_TO_RAD
  use pmemd_lib_mod

  implicit none
  
! Formal arguments:

  integer               :: natom
  character(len=4)      :: isymbl(*)
  integer               :: nghbrs(5, *)
  integer               :: hnghbrs(5, *)
  integer               :: enghbrs(5, *)
  integer               :: numnghbr(3, *)
  type(ep_frame_rec)    :: frames(*)
  double precision      :: ep_lcl_crd(3, 2, *)
  integer               :: frame_cnt
  integer               :: epbtyp(5, *)
  double precision      :: req(*)
  integer               :: verbose

! Local variables:

  integer               :: k, m, n, l
  double precision      :: tetcos, tetsin, angle, scos, ssin
  character(len=4)      :: sulf, sulfh

  sulf = 'S   '
  sulfh = 'SH  '
   
  ! Get half-angle for tetrahedral:
   
  angle = 54.735d0
  angle = angle * DEG_TO_RAD
  tetcos = cos(angle)
  tetsin = sin(angle)
   
  ! Get cos, sin for 60:
   
  scos = 0.5d0
  ssin = sqrt(1.d0 - scos * scos)
   
  frame_cnt = 0

  do n = 1, natom
      
    if (numnghbr(3, n) .gt. 0) then
         
      if (numnghbr(1, n) + numnghbr(2, n) .gt. 2) then
        write(6, *) 'EXTRA_PTS: too many nghbrs!!'
        call mexit(6, 1)
      end if

      frame_cnt = frame_cnt + 1
      frames(frame_cnt)%parent_atm = n
      frames(frame_cnt)%ep_cnt = numnghbr(3, n)

      frames(frame_cnt)%extra_pnt(:) = 0

      do k = 1, frames(frame_cnt)%ep_cnt
        frames(frame_cnt)%extra_pnt(k) = enghbrs(k, n)
      end do
         
      if (numnghbr(1, n) .eq. 0  .and. &
          numnghbr(2, n) .eq. 2 .and. &
          numnghbr(3, n) .eq. 1) then
            
        ! TIP4P water:
        ! (temporarily assign type to 3, will set to 1 in next section)
            
        frames(frame_cnt)%type = 3
        frames(frame_cnt)%frame_atm1 = hnghbrs(1, n)
        frames(frame_cnt)%frame_atm2 = n
        frames(frame_cnt)%frame_atm3 = hnghbrs(2, n)
            
      else if (numnghbr(1, n) .eq. 0  .and. &
               numnghbr(2, n) .eq. 2 .and. &
               numnghbr(3, n) .eq. 2) then
            
        ! TIP5P water:
            
        frames(frame_cnt)%type = 1
        frames(frame_cnt)%frame_atm1 = hnghbrs(1, n)
        frames(frame_cnt)%frame_atm2 = n
        frames(frame_cnt)%frame_atm3 = hnghbrs(2, n)
            
      else if (numnghbr(1, n) .gt. 1) then
            
        ! "ordinary" type of frame defined by two other heavy atoms:
            
        frames(frame_cnt)%type = 1
        frames(frame_cnt)%frame_atm1 = nghbrs(1, n)
        frames(frame_cnt)%frame_atm2 = n
        frames(frame_cnt)%frame_atm3 = nghbrs(2, n)
            
      else if (numnghbr(1, n) .eq. 1 .and. numnghbr(2, n) .eq. 1) then
            
        ! frame defined by one heavy atom and one hydrogen:
            
        frames(frame_cnt)%type = 1
        frames(frame_cnt)%frame_atm1 = nghbrs(1, n)
        frames(frame_cnt)%frame_atm2 = n
        frames(frame_cnt)%frame_atm3 = hnghbrs(1, n)
            
      else if (numnghbr(1, n) .eq. 1 .and. numnghbr(2, n) .eq. 0) then
            
        ! Assume this is CARBONYL oxygen:
        ! (Need to use midpoints of other two bonds of the carbon
        ! for frame_atm1 and frame_atm3 in orient force; thus (mis)use
        ! frame_atm1 frame_atm2 frame_atm3 and parent_atm to store 4 atoms for
        ! this special case.)
            
        frames(frame_cnt)%type = 2
        m = nghbrs(1, n)
        frames(frame_cnt)%frame_atm2 = m
        if (numnghbr(1, m) .ne. 3 .or.numnghbr(2, m) .gt. 0) then
          write(6, *) 'EXTRA_PTS: frame type 2 Should not be here'
          write(6, *) n, m, numnghbr(1, m), numnghbr(2, m)
          call mexit(6, 1)
        end if
            
        ! numnghbr(1, m) = 3. One is n (the oxygen) and other 2 are
        ! other bonding partners of carbon:

        frames(frame_cnt)%frame_atm1 = 0
        frames(frame_cnt)%frame_atm3 = 0

        k = 1
        do while (k .lt. 4 .and. frames(frame_cnt)%frame_atm1 .eq. 0)

          if (nghbrs(k, m) .ne. n) then
            frames(frame_cnt)%frame_atm1 = nghbrs(k, m)
          end if
          k = k + 1

        end do

        k = 1
        do while (k .lt. 4 .and. frames(frame_cnt)%frame_atm3 .eq. 0)

          if (nghbrs(k, m) .ne. n .and. &
              nghbrs(k, m) .ne. frames(frame_cnt)%frame_atm1) then
            frames(frame_cnt)%frame_atm3 = nghbrs(k, m)
          end if
          k = k + 1

        end do

        if (frames(frame_cnt)%frame_atm1.eq. 0 .or. &
            frames(frame_cnt)%frame_atm3 .eq. 0) then
          write(6, *) 'EXTRA_PTS: cannot find first or third frame point '
          write(6, *) 'define: ', n, numnghbr(1, n), numnghbr(2, n), &
            numnghbr(3, n), frames(frame_cnt)%frame_atm1, &
            frames(frame_cnt)%frame_atm3
          call mexit(6, 1)
        endif

      else

        write(6, *) 'EXTRA_PTS: unexpected numnghbr array: '
        write(6, *) 'define: ', n, numnghbr(1, n), numnghbr(2, n), &
          numnghbr(3, n)
        call mexit(6, 1)

      end if  !  ( numnghbr(1, n) .eq. 0  .and. numnghbr(2, n) .eq. 2

    end if  !  ( numnghbr(1, n) .eq. 0  .and. numnghbr(2, n) .eq. 2

  end do  ! ( numnghbr(3, n) .gt. 0 )
   
  ! Get the local coords:
   
  ep_lcl_crd(1:3, 1:2, 1:frame_cnt) = 0.d0

  do n = 1, frame_cnt

    l = frames(n)%parent_atm

    if (frames(n)%type .eq. 1 .or. frames(n)%type .eq. 3) then
         
      ! Z axis is along symmetry axis of second atom opposite
      ! bisector of frame_atm1, frame_atm3;
      ! X axis is along the diff vector frame_atm3 minus frame_atm1;
      ! Y axis is cross product.
         
      if (frames(n)%ep_cnt .eq. 1) then
            
        ! Extra point is along the z-direction: positive for ordinary
        ! lone pair, negative for TIP4P water extra point:
            
        ep_lcl_crd(3, 1, n) = req(epbtyp(1, l))
        if (frames(n)%type .eq. 3) then
          ep_lcl_crd(3, 1, n) = -req(epbtyp(1, l))
          frames(n)%type = 1
        end if
            
      else if (frames(n)%ep_cnt .eq. 2) then
            
        ! Extra points are in the z, y plane, tetrahedrally
        ! (unless frame_atm2 atom is sulfur, in which case they
        ! are opposite along y):
            
        m = frames(n)%frame_atm2
        if (isymbl(m) .eq. sulf .or. isymbl(m) .eq. sulfh) then
          ep_lcl_crd(2, 1, n) = req(epbtyp(1, l))
          ep_lcl_crd(2, 2, n) = -req(epbtyp(2, l))
        else
          ep_lcl_crd(3, 1, n) = tetcos * req(epbtyp(1, l))
          ep_lcl_crd(2, 1, n) = tetsin * req(epbtyp(1, l))
          ep_lcl_crd(3, 2, n) = tetcos * req(epbtyp(2, l))
          ep_lcl_crd(2, 2, n) = -tetsin * req(epbtyp(2, l))
        end if
            
      else

        write(6, *) 'EXTRA_PTS: unexpected ep_cnt value: ', frames(n)%ep_cnt
        call mexit(6, 1)

      end if  ! ( frames(n)%ep_cnt .eq. 1 )
         
    else if (frames(n)%type .eq. 2) then
         
      ! Z axis is along bond from frame_atm2 to parent_atm;
      ! X axis in plane of parent_atm and midpoints of frame_atm1, frame_atm2
      ! and frame_atm2, frame_atm3.
         
      if (frames(n)%ep_cnt .eq. 1) then
        ep_lcl_crd(3, 1, n) = req(epbtyp(1, l))
      else if (frames(n)%ep_cnt .eq. 2) then
        ep_lcl_crd(3, 1, n) = scos * req(epbtyp(1, l))
        ep_lcl_crd(1, 1, n) = ssin * req(epbtyp(1, l))
        ep_lcl_crd(3, 2, n) = scos * req(epbtyp(2, l))
        ep_lcl_crd(1, 2, n) = -ssin * req(epbtyp(2, l))
      else
        write(6, *) 'EXTRA_PTS: unexpected ep_cnt value: ', frames(n)%ep_cnt
        call mexit(6, 1)
      end if

    else
      write(6, *) 'EXTRA_PTS: unexpected frame type value: ', frames(n)%type
      call mexit(6, 1)
    end if  ! ( frames(n)%type .eq. 1 .or. frames(n)%type .eq. 3 )

  end do  !  n = 1, frame_cnt
   
   if (verbose .gt. 3) then
     write(6, *) 'frames:'
     do n = 1, frame_cnt
       write(6, 666) n, frames(n)%parent_atm, isymbl(frames(n)%parent_atm), &
         frames(n)%ep_cnt, &
         frames(n)%extra_pnt(1), frames(n)%extra_pnt(2), frames(n)%type, &
         frames(n)%frame_atm1, frames(n)%frame_atm2, frames(n)%frame_atm3
       write(6, 667) ep_lcl_crd(1,1,n), ep_lcl_crd(2,1,n), ep_lcl_crd(3,1,n), &
         ep_lcl_crd(1,2,n), ep_lcl_crd(2,2,n), ep_lcl_crd(3,2,n)
     end do
   end if
   
  return

666 format(1x, 2i7, 1x, a4, 7i7)
667 format(1x, 6(1x, f10.4))

end subroutine define_frames 

!*******************************************************************************!
! Internal Subroutine:  fill_bonded
!
! Description: <TBS>
!
!*******************************************************************************

subroutine fill_bonded(max11, max12, max13, max14, &
                       num11, num12, num13, num14, &
                       list11, list12, list13, list14, &
                       enghbrs, numnghbr, epowner, verbose)

  use pmemd_lib_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer       :: max11
  integer       :: max12
  integer       :: max13
  integer       :: max14
  integer       :: num11
  integer       :: num12
  integer       :: num13
  integer       :: num14
  integer       :: list11(2, max11)
  integer       :: list12(2, max12)
  integer       :: list13(2, max13)
  integer       :: list14(2, max14)
  integer       :: enghbrs(5, *)
  integer       :: numnghbr(3, *)
  integer       :: epowner(*)
  integer       :: verbose

! Local variables:

  integer       :: n, ifail

  num11 = 0
  num12 = 0
  num13 = 0
  num14 = 0

  do n = 1, natom

    if (numnghbr(3, n) .eq. 1) then
      num11 = num11 + 1
      if (num11 .gt. max11) goto 100
      list11(1, num11) = n
      list11(2, num11) = enghbrs(1, n)
    else if (numnghbr(3, n) .eq. 2) then
      if (num11  + 3 .gt. max11) goto 100
      num11 = num11 + 1
      if (num11 .gt. max11) goto 100
      list11(1, num11) = n
      list11(2, num11) = enghbrs(1, n)
      num11 = num11 + 1
      if (num11 .gt. max11) goto 100
      list11(1, num11) = n
      list11(2, num11) = enghbrs(2, n)
      num11 = num11 + 1
      if (num11 .gt. max11) goto 100
      list11(1, num11) = enghbrs(1, n)
      list11(2, num11) = enghbrs(2, n)
    else if (numnghbr(3, n) .eq. 3) then
      write(6, *) 'fill_bonded: should not be here!'
      call mexit(6, 1)
    end if

  end do

  ! Bonds:

  call do_bond_pairs(list12, num12, max12, epowner, numnghbr, enghbrs, ifail)

  if (ifail .eq. 1) then
    write(6, *) 'fill_bonded: max12 exceeded!!'
    call mexit(6, 1)
  end if

  call sort_pairs(list12, num12, natom)

  ! Angles:

  call do_angle_pairs(list13, num13, max13, epowner, numnghbr, enghbrs, ifail)

  if (ifail .eq. 1) then
    write(6, *) 'fill_bonded: max13 exceeded!!'
    call mexit(6, 1)
  end if

  call sort_pairs(list13, num13, natom)

  ! Dihedrals:

  call do_dihed_pairs(list14, num14, max14, epowner, numnghbr, enghbrs, ifail)

  if (ifail .eq. 1) then
    write(6, *) 'fill_bonded: max14 exceeded!!'
    call mexit(6, 1)
  end if

  call sort_pairs(list14, num14, natom)
   
  if (verbose .gt. 0) &
    write(6, '(a, 4i6)') '| EXTRA PNTS fill_bonded: num11-14 = ', &
      num11, num12, num13, num14

  if (verbose .gt. 3) then
    write(6, *) '$$$$$$$$$$$$$$$$$ 1-1 pairs $$$$$$$$$$$$$$$$$$$$$$$$'
    do n = 1, num11
      write(6, 666) n, list11(1, n), list11(2, n)
    end do
    write(6, *) '$$$$$$$$$$$$$$$$$ 1-2 pairs $$$$$$$$$$$$$$$$$$$$$$$$'
    do n = 1, num12
      write(6, 666) n, list12(1, n), list12(2, n)
    end do
    write(6, *) '$$$$$$$$$$$$$$$$$ 1-3 pairs $$$$$$$$$$$$$$$$$$$$$$$$'
    do n = 1, num13
      write(6, 666) n, list13(1, n), list13(2, n)
    end do
    write(6, *) '$$$$$$$$$$$$$$$$$ 1-4 pairs $$$$$$$$$$$$$$$$$$$$$$$$'
    do n = 1, num14
      write(6, 666) n, list14(1, n), list14(2, n)
    end do
    write(6, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
  end if

  return

100 write(6, *) 'fill_bonded: max11 exceeded!!'
  call mexit(6, 1)

666 format(1x, i5, ':', 3x, 'i, j = ', i5, 1x, i5)

end subroutine fill_bonded 

!*******************************************************************************!
! Internal Subroutine:  redo_masked
!
! Description: <TBS>
!
!*******************************************************************************

subroutine redo_masked(natom, iblo, inb, nnb, &
                       num11, num12, num13, num14, &
                       list11, list12, list13, list14, offset, test)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer       :: natom
  integer       :: iblo(*)
  integer       :: inb(*)
  integer       :: nnb
  integer       :: num11
  integer       :: num12
  integer       :: num13
  integer       :: num14
  integer       :: list11(2, *)
  integer       :: list12(2, *)
  integer       :: list13(2, *)
  integer       :: list14(2, *)
  integer       :: offset(*)
  integer       :: test(*)
   
! Local variables:

  integer       :: j, n, m, ntot

  ! Build the mask list from list11-14. Make sure no duplication:
   
  do n = 1, natom
    iblo(n) = 0
    offset(n) = 0
    test(n) = 0
  end do
   
  ! Pass 1:  Fill iblo.
   
  call add_one_list_iblo(iblo, list11, num11)
  call add_one_list_iblo(iblo, list12, num12)
  call add_one_list_iblo(iblo, list13, num13)
  call add_one_list_iblo(iblo, list14, num14)
   
  ! Check totals while finding offsets, resetting iblo:
   
  ntot = 0
  do n = 1, natom
      offset(n) = ntot
      ntot = ntot + iblo(n)
      iblo(n) = 0
  end do

  if (ntot .gt. 2 * nnb) then
      write(6, *) 'EXTRA POINTS: nnb too small! '
      write(6, *) 'nnb, ntot = ', nnb, ntot
      call mexit(6, 1)
  end if
   
  ! Pass 2 fill inb, redo iblo:
   
  call add_one_list_inb(iblo, inb, offset, list11, num11)
  call add_one_list_inb(iblo, inb, offset, list12, num12)
  call add_one_list_inb(iblo, inb, offset, list13, num13)
  call add_one_list_inb(iblo, inb, offset, list14, num14)
   
  ! Pass 3 filter inb, remove duplicate entries:
   
  do n = 1, natom - 1
    do m = 1, iblo(n)
      j = inb(offset(n) + m)
      if (test(j) .ne. n ) then
        test(j) = n
      else
        inb(offset(n) + m) = 0
      end if
    end do
  end do

  return

end subroutine redo_masked 

!*******************************************************************************!
! Internal Subroutine:  build_nb14
!
! Description: <TBS>
!
!*******************************************************************************

subroutine build_nb14(nb14, nb14_cnt, maxnb14, epowner, numnghbr, enghbrs, &
                      atm_cnt, chngmask, verbose)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer               :: nb14(3, *)
  integer               :: nb14_cnt
  integer               :: maxnb14
  integer               :: epowner(*)
  integer               :: numnghbr(3, *)
  integer               :: enghbrs(5, *)
  integer               :: atm_cnt
  integer               :: chngmask
  integer               :: verbose

! Local variables:

  integer               :: ifail, n

  nb14_cnt = 0

  call do_14pairs(nb14, nb14_cnt, maxnb14, epowner, numnghbr, enghbrs, &
                  ifail, chngmask)

  if (ifail .eq. 1) then
    write(6, *) 'exceeded maxnb14 in build14: check extra_pnts_nb14.fpp'
    call mexit(6, 1)
  end if

  call sort_pairs_14_nb(nb14, nb14_cnt, atm_cnt)

  if (verbose .gt. 0) write(6, '(a, i6)') &
    '| EXTRA_PTS, build_nb14: num of 14 terms = ', nb14_cnt

  if (verbose .gt. 3) then
    write(6, *) '$$$$$$$$$$$$$$$$$$$$$$$  1-4 nb list $$$$$$$$$$'
    do n = 1, nb14_cnt
      write(6, 666) n, nb14(1, n), nb14(2, n)
    end do
    write(6, *) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
  end if

  return

666 format(1x, i5, ':', 3x, 'i, j = ', i5, 1x, i5)

end subroutine build_nb14 

!*******************************************************************************!
! Internal Subroutine:  copy_nb14
!
! Description: <TBS>
!
!*******************************************************************************

subroutine copy_nb14(from14, to14, nb14_cnt)

  implicit none

! Formal arguments:

  integer       :: from14(3, *)
  integer       :: to14(3, *)
  integer       :: nb14_cnt

! Local variables:

  integer       :: n

  do n = 1, nb14_cnt
    to14(1, n) = from14(1, n)
    to14(2, n) = from14(2, n)
    to14(3, n) = from14(3, n)
  end do

  return

end subroutine copy_nb14 

!*******************************************************************************!
! Internal Subroutine:  trim_bonds
!
! Description: <TBS>
!
!*******************************************************************************

subroutine trim_bonds(bond_cnt, bonds, epowner)

  implicit none

! Formal arguments:

  integer, intent(in out)               :: bond_cnt
  type(bond_rec), intent(in out)        :: bonds(*)
  integer, intent(in)                   :: epowner(*)

! Local variables:

  integer       :: atm_i, atm_j
  integer       :: n, m
   
  write(6, '(a, 2i6)') &
    '|      EXTRA_PTS, trim_bonds: num bonds BEFORE trim =', bond_cnt, 0
   
  m = 0

  do n = 1, bond_cnt

    atm_i = bonds(n)%atm_i
    atm_j = bonds(n)%atm_j
      
    ! Only keep if neither is extra:
      
    if (epowner(atm_i) .eq. 0 .and. epowner(atm_j) .eq. 0) then
      m = m + 1
      bonds(m) = bonds(n)
    end if

  end do

  bond_cnt = m
   
  write(6, '(a, 2i6)') &
    '|      EXTRA_PTS, trim_bonds: num bonds AFTER  trim =', bond_cnt, 0

  return

end subroutine trim_bonds 

!*******************************************************************************!
! Internal Subroutine:  trim_theta
!
! Description: <TBS>
!
!*******************************************************************************

subroutine trim_theta(angle_cnt, angles, epowner)

  implicit none

! Formal arguments:

  integer, intent(in out)               :: angle_cnt
  type(angle_rec), intent(in out)       :: angles(*)
  integer, intent(in)                   :: epowner(*)

! Local variables:

  integer       :: atm_i, atm_k
  integer       :: n, m

  m = 0

  write(6, '(a, 2i6)') &
    '|      EXTRA_PTS, trim_theta: num angle BEFORE trim =', angle_cnt, 0
   
  do n = 1, angle_cnt

    atm_i = angles(n)%atm_i
    atm_k = angles(n)%atm_k
      
    ! Only keep if neither is extra:
      
    if (epowner(atm_i) .eq. 0 .and. epowner(atm_k) .eq. 0) then
      m = m + 1
      angles(m) = angles(n)
    end if

  end do

  angle_cnt = m
   
  write(6, '(a, 2i6)') &
    '|      EXTRA_PTS, trim_theta: num angle AFTER  trim =', angle_cnt, 0

  return

end subroutine trim_theta 

!*******************************************************************************!
! Internal Subroutine:  trim_phi
!
! Description: <TBS>
!
!*******************************************************************************

subroutine trim_phi(dihed_cnt, dihedrals, epowner)

  implicit none

! Formal arguments:

  integer, intent(in out)               :: dihed_cnt
  type(dihed_rec), intent(in out)       :: dihedrals(*)
  integer, intent(in)                   :: epowner(*)

! Local variables:

  integer       :: atm_i, atm_l
  integer       :: n, m

  m = 0

  write(6, '(a, 2i6)') &
    '|      EXTRA_PTS, trim_phi:  num diheds BEFORE trim =', dihed_cnt, 0
   
  do n = 1, dihed_cnt

    atm_i = dihedrals(n)%atm_i
    atm_l = iabs(dihedrals(n)%atm_l)
      
    ! Only keep if neither is extra:
      
    if (epowner(atm_i) .eq. 0 .and. epowner(atm_l) .eq. 0) then
      m = m + 1
      dihedrals(m) = dihedrals(n)
    end if

  end do

  dihed_cnt = m
   
  write(6, '(a, 2i6)') &
    '|      EXTRA_PTS, trim_phi:  num diheds AFTER  trim =', dihed_cnt, 0

  return

end subroutine trim_phi 

!*******************************************************************************!
! Internal Subroutine:  do_bond_pairs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine do_bond_pairs(list, num, maxp, epowner, numnghbr, enghbrs, ifail)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer       :: list(2, *)
  integer       :: num
  integer       :: maxp
  integer       :: epowner(*)
  integer       :: numnghbr(3, *)
  integer       :: enghbrs(5, *)
  integer       :: ifail

! Local variables:

  integer       :: atm_i, atm_j
  integer       :: i, j, k, l, n

  ifail = 0

  do n = 1, nbonh + nbona
      
    atm_i = gbl_bond(n)%atm_i
    atm_j = gbl_bond(n)%atm_j
      
    ! Check neither is extra; Count this bond and also extra attached:
      
    if (epowner(atm_i) .eq. 0 .and. epowner(atm_j) .eq. 0) then
      k = numnghbr(3, atm_i)
      l = numnghbr(3, atm_j)
      if (num + 1 + k + l + k * l .gt. maxp) then
        ifail = 1
        return
      end if
      num = num + 1
      list(1, num) = atm_i
      list(2, num) = atm_j
      do i = 1, k
        num = num + 1
        list(1, num) = enghbrs(i, atm_i)
        list(2, num) = atm_j
      end do
      do j = 1, l
        num = num + 1
        list(1, num) = atm_i
        list(2, num) = enghbrs(j, atm_j)
      end do
      do i = 1, k
        do j = 1, l
          num = num + 1
          list(1, num) = enghbrs(i, atm_i)
          list(2, num) = enghbrs(j, atm_j)
        end do
      end do
    end if
  end do  !  n = 1, nbonh + nbona

  return

end subroutine do_bond_pairs 

!*******************************************************************************!
! Internal Subroutine:  do_angle_pairs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine do_angle_pairs(list, num, maxp, epowner, numnghbr, enghbrs, ifail)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer       :: list(2, *)
  integer       :: num
  integer       :: maxp
  integer       :: epowner(*)
  integer       :: numnghbr(3, *)
  integer       :: enghbrs(5, *)
  integer       :: ifail

! Local variables:

  integer       :: atm_i, atm_k
  integer       :: i, j, k, l, n

  ifail = 0

  do n = 1, ntheth + ntheta
      
    atm_i = gbl_angle(n)%atm_i
    atm_k = gbl_angle(n)%atm_k
      
    ! Check neither is extra; Count this bond and also extra attached:
      
    if (epowner(atm_i) .eq. 0 .and. epowner(atm_k) .eq. 0) then
      k = numnghbr(3, atm_i)
      l = numnghbr(3, atm_k)
      if (num + 1 + k + l + k * l .gt. maxp) then
        ifail = 1
        return
      end if
      num = num + 1
      list(1, num) = atm_i
      list(2, num) = atm_k
      do i = 1, k
        num = num + 1
        list(1, num) = enghbrs(i, atm_i)
        list(2, num) = atm_k
      end do
      do j = 1, l
        num = num + 1
        list(1, num) = atm_i
        list(2, num) = enghbrs(j, atm_k)
      end do
      do i = 1, k
        do j = 1, l
          num = num + 1
          list(1, num) = enghbrs(i, atm_i)
          list(2, num) = enghbrs(j, atm_k)
        end do
      end do
    end if
  end do  !  n = 1, ntheth + ntheta

  return

end subroutine do_angle_pairs 

!*******************************************************************************!
! Internal Subroutine:  do_dihed_pairs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine do_dihed_pairs(list, num, maxp, epowner, numnghbr, enghbrs, ifail)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer       :: list(2, *)
  integer       :: num
  integer       :: maxp
  integer       :: epowner(*)
  integer       :: numnghbr(3, *)
  integer       :: enghbrs(5, *)
  integer       :: ifail

! Local variables:

  integer       :: atm_i, atm_l
  integer       :: i, j, k, l, n

  ifail = 0

  do n = 1, nphih + nphia
      
    ! Sometimes second index negative (improper dihedrals):
      
    atm_i = gbl_dihed(n)%atm_i
    atm_l = iabs(gbl_dihed(n)%atm_l)
      
    ! Check neither is extra; Count this bond and also extra attached:
      
    if (epowner(atm_i) .eq. 0 .and. epowner(atm_l) .eq. 0) then
      k = numnghbr(3, atm_i)
      l = numnghbr(3, atm_l)
      if (num + 1 + k + l + k * l .gt. maxp) then
        ifail = 1
        return
      end if
      num = num + 1
      list(1, num) = atm_i
      list(2, num) = atm_l
      do i = 1, k
        num = num + 1
        list(1, num) = enghbrs(i, atm_i)
        list(2, num) = atm_l
      end do
      do j = 1, l
        num = num + 1
        list(1, num) = atm_i
        list(2, num) = enghbrs(j, atm_l)
      end do
      do i = 1, k
        do j = 1, l
          num = num + 1
          list(1, num) = enghbrs(i, atm_i)
          list(2, num) = enghbrs(j, atm_l)
        end do
      end do
    end if
  end do  !  n = 1, nphih + nphia

  return

end subroutine do_dihed_pairs 

!*******************************************************************************!
! Internal Subroutine:  sort_pairs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine sort_pairs(list, num, atm_cnt)

  implicit none

! Formal arguments:

  integer       :: list(2, *)
  integer       :: num
  integer       :: atm_cnt

! Local variables:

  integer       :: i, j, k, m, n, ntot
  integer       :: scr1(atm_cnt), scr2(atm_cnt), scr3(maxa)

  do n = 1, num
    i = list(1, n)
    j = list(2, n)
    if (i .lt. j) then
      list(1, n) = i
      list(2, n) = j
    else
      list(1, n) = j
      list(2, n) = i
    end if
  end do
   
  ! Now get rid of duplicates:
   
  ! First pass:
   
  do n = 1, atm_cnt
    scr1(n) = 0
  end do

  do n = 1, num
    i = list(1, n)
    scr1(i) = scr1(i) + 1
  end do

  scr2(1) = 0
  do n = 2, atm_cnt
    scr2(n) = scr2(n - 1) + scr1(n - 1)
    scr1(n - 1) = 0
  end do
  scr1(atm_cnt) = 0
   
  ! Second pass:
   
  do n = 1, num
    i = list(1, n)
    j = list(2, n)
    scr1(i) = scr1(i) + 1
    scr3(scr1(i) + scr2(i)) = j
  end do

  scr2(1:atm_cnt) = 0
   
  ! Now trim them:
   
  ntot = 0
  k = 0
  do n = 1, atm_cnt
    do m = 1, scr1(n)
      j = scr3(ntot + m)
      if (scr2(j) .ne. n) then
        k = k + 1
        list(1, k) = n
        list(2, k) = j
        scr2(j) = n
      end if
    end do
    ntot = ntot + scr1(n)
  end do
  num = k
   
  return
   
end subroutine sort_pairs 

!*******************************************************************************!
! Internal Subroutine:  sort_pairs_14_nb
!
! Description: <TBS>
!
!*******************************************************************************

subroutine sort_pairs_14_nb(nb_14_list, num, atm_cnt)

  implicit none

! Formal arguments:

  integer       :: nb_14_list(3, *)
  integer       :: num
  integer       :: atm_cnt

! Local variables:

  integer       :: i, j, k, m, n, ntot, parm_idx
  integer       :: scr1(atm_cnt), scr2(atm_cnt), scr3(maxa), scr4(maxa)

 !Indexes of nb_14_list = 1,j,parm_index
 !First pass, make sure the first index contains the lowest number. Ignore
 !index 3 here since this won't change anything.

  do n = 1, num
    i = nb_14_list(1, n)
    j = nb_14_list(2, n)
    if (i .lt. j) then
      nb_14_list(1, n) = i
      nb_14_list(2, n) = j
    else
      nb_14_list(1, n) = j
      nb_14_list(2, n) = i
    end if
  end do
   
  ! Now get rid of duplicates:
   
  ! First pass:
  scr1(1:atm_cnt) = 0 
  do n = 1, num
    i = nb_14_list(1, n)
    scr1(i) = scr1(i) + 1
  end do

  scr2(1) = 0
  do n = 2, atm_cnt
    scr2(n) = scr2(n - 1) + scr1(n - 1)
    scr1(n - 1) = 0
  end do
  scr1(atm_cnt) = 0
   
  ! Second pass:
   
  do n = 1, num
    i = nb_14_list(1, n)
    j = nb_14_list(2, n)
    parm_idx = nb_14_list(3, n)
    scr1(i) = scr1(i) + 1
    scr3(scr1(i) + scr2(i)) = j
    scr4(scr1(i) + scr2(i)) = parm_idx
  end do

  scr2(1:atm_cnt) = 0
   
  ! Now trim them:
   
  ntot = 0
  k = 0
  do n = 1, atm_cnt
    do m = 1, scr1(n)
      j = scr3(ntot + m)
      parm_idx = scr4(ntot + m)
      if (scr2(j) .ne. n) then
        k = k + 1
        nb_14_list(1, k) = n
        nb_14_list(2, k) = j
        nb_14_list(3, k) = parm_idx
        scr2(j) = n
      end if
    end do
    ntot = ntot + scr1(n)
  end do
  num = k
   
  return

end subroutine sort_pairs_14_nb

!*******************************************************************************!
! Internal Subroutine:  add_one_list_iblo
!
! Description: <TBS>
!
!*******************************************************************************

subroutine add_one_list_iblo(iblo, list, num)

  implicit none

! Formal arguments:

  integer       :: iblo(*)
  integer       :: list(2, *)
  integer       :: num

! Local variables:

  integer       :: n, i

  do n = 1, num
    i = list(1, n)
    iblo(i) = iblo(i) + 1
  end do

  return

end subroutine add_one_list_iblo 

!*******************************************************************************!
! Internal Subroutine:  add_one_list_inb
!
! Description: <TBS>
!
!*******************************************************************************

subroutine add_one_list_inb(iblo, inb, offset, list, num)

  implicit none

! Formal arguments:

  integer       :: iblo(*)
  integer       :: inb(*)
  integer       :: offset(*)
  integer       :: list(2, *)
  integer       :: num

! Local variables:

  integer       :: n, i, j, m

  do n = 1, num
    i = list(1, n)
    j = list(2, n)
    m = offset(i)
    iblo(i) = iblo(i) + 1
    inb(m + iblo(i)) = j
  end do

  return

end subroutine add_one_list_inb 

!*******************************************************************************!
! Internal Subroutine:  do_14pairs
!
! Description: <TBS>
!
!*******************************************************************************

subroutine do_14pairs(nb_14_list, num, maxnb14, epowner, numnghbr, enghbrs, &
                      ifail, chngmask)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer               :: nb_14_list(3, *)
  integer               :: num
  integer               :: maxnb14
  integer               :: epowner(*)
  integer               :: numnghbr(3, *)
  integer               :: enghbrs(5, *)
  integer               :: ifail
  integer               :: chngmask

! Local variables:

  integer               :: atm_i, atm_l, parm_idx
  integer               :: i, j, n, ni, nl

  ifail = 0

  do n = 1, nphih + nphia

    if (gbl_dihed(n)%atm_k .gt. 0 .and. &
        gbl_dihed(n)%atm_l .gt. 0 .and. &
        gbl_fmn(gbl_dihed(n)%parm_idx) .gt. 0.d0) then

      atm_i = gbl_dihed(n)%atm_i
      atm_l = gbl_dihed(n)%atm_l
      parm_idx = gbl_dihed(n)%parm_idx
      if (chngmask .eq. 0) then
        if (num + 1 .gt. maxnb14) then
          ifail = 1
          return
        end if
        num = num + 1
        nb_14_list(1, num) = atm_i
        nb_14_list(2, num) = atm_l
        nb_14_list(3, num) = parm_idx  !used for lookup into scnb and scee arrays.
      else
        ! Check neither is extra. Count this bond and also extra attached:
        if (epowner(atm_i) .eq. 0 .and. epowner(atm_l) .eq. 0) then
          ni = numnghbr(3, atm_i)
          nl = numnghbr(3, atm_l)
          if (num + 1 + ni + nl + ni * nl .gt. maxnb14) then
            ifail = 1
            return
          end if
          num = num + 1
          nb_14_list(1, num) = atm_i
          nb_14_list(2, num) = atm_l
          nb_14_list(3, num) = parm_idx  !used for lookup into scnb and scee arrays.
          do i = 1, ni
            num = num + 1
            nb_14_list(1, num) = enghbrs(i, atm_i)
            nb_14_list(2, num) = atm_l
            nb_14_list(3, num) = parm_idx  !used for lookup into scnb and scee arrays.
          end do
          do j = 1, nl
            num = num + 1
            nb_14_list(1, num) = atm_i
            nb_14_list(2, num) = enghbrs(j, atm_l)
            nb_14_list(3, num) = parm_idx  !used for lookup into scnb and scee arrays.
          end do
          do i = 1, ni
            do j = 1, nl
              num = num + 1
              nb_14_list(1, num) = enghbrs(i, atm_i)
              nb_14_list(2, num) = enghbrs(j, atm_l)
              nb_14_list(3, num) = parm_idx  !used for lookup into scnb and scee arrays.
            end do
          end do
        end if
      end if  ! ( chngmask .eq. 0 )
    end if
  end do  !  n = 1, nphih + nphia

  return

end subroutine do_14pairs

!*******************************************************************************!
! Subroutine:  all_local_to_global
!
! Description:  Put EP in position in world coord system based on the
!               position of the frame and the local coordinates, for all
!               coordinates; done only in master.
!
! Frames from Stone and Alderton Mol Phys. 56, 5, 1047 (1985) plus weird
! modification for carbonyl to symmetrize.
!*******************************************************************************

subroutine all_local_to_global(crd, frames, lcl_crd, frame_cnt)

  implicit none

! Formal arguments:

  double precision      :: crd(3, *)
  type(ep_frame_rec)    :: frames(*)
  double precision      :: lcl_crd(3, 2, *)
  integer               :: frame_cnt
   
! Local variables:

  double precision      :: uvec(3)
  double precision      :: vvec(3)
  double precision      :: ave(3)
  double precision      :: diff(3)
  double precision      :: usiz
  double precision      :: vsiz
  double precision      :: asiz
  double precision      :: dsiz
  double precision      :: f(3, 3)
  double precision      :: a(3), b(3), c(3)
  integer               :: j, k, m, n
   
  if (frame_cnt .eq. 0) return

  do n = 1, frame_cnt

    if (frames(n)%type .eq. 1) then
      do m = 1, 3
        a(m) = crd(m, frames(n)%frame_atm1)
        b(m) = crd(m, frames(n)%frame_atm2)
        c(m) = crd(m, frames(n)%frame_atm3)
      end do
    else if (frames(n)%type .eq. 2) then
      do m = 1, 3
        a(m) = 0.5d0 * (crd(m, frames(n)%frame_atm1) + &
               crd(m, frames(n)%frame_atm2))
        b(m) = crd(m, frames(n)%parent_atm)
        c(m) = 0.5d0 * (crd(m, frames(n)%frame_atm3) + &
               crd(m, frames(n)%frame_atm2))
      end do
    end if
      
    ! Z-axis along symmmetry axis of b midway between
    ! unit vector to a and unit vector to c; points opposite:
      
    usiz = 0.d0
    vsiz = 0.d0
    do m = 1, 3
      uvec(m) = a(m) - b(m)
      usiz = usiz + uvec(m) * uvec(m)
      vvec(m) = c(m) - b(m)
      vsiz = vsiz + vvec(m) * vvec(m)
    end do
    usiz = sqrt(usiz)
    vsiz = sqrt(vsiz)
    asiz = 0.d0
    dsiz = 0.d0
    do m = 1, 3
      uvec(m) = uvec(m) / usiz
      vvec(m) = vvec(m) / vsiz
      ave(m) = (uvec(m) + vvec(m)) / 2.d0
      asiz = asiz + ave(m) * ave(m)
      diff(m) = (vvec(m) - uvec(m)) / 2.d0
      dsiz = dsiz + diff(m) * diff(m)
    end do
    asiz = sqrt(asiz)
    dsiz = sqrt(dsiz)
    do m = 1, 3
      f(m, 3) = -ave(m) / asiz
      f(m, 1) = diff(m) / dsiz
    end do
    f(1, 2) = f(2, 3) * f(3, 1) - f(3, 3) * f(2, 1)
    f(2, 2) = f(3, 3) * f(1, 1) - f(1, 3) * f(3, 1)
    f(3, 2) = f(1, 3) * f(2, 1) - f(2, 3) * f(1, 1)
    do k = 1, frames(n)%ep_cnt
      j = frames(n)%extra_pnt(k)
      do m = 1, 3
        crd(m, j) = crd(m, frames(n)%parent_atm) + &
                    lcl_crd(1, k, n) * f(m, 1) + &
                    lcl_crd(2, k, n) * f(m, 2) + &
                    lcl_crd(3, k, n) * f(m, 3)
      end do
    end do

  end do  !  n = 1, frame_cnt

  return

end subroutine all_local_to_global 

!*******************************************************************************!
! Subroutine:  local_to_global
!
! Description:  Put EP in position in world coord system based on the
!               position of the frame and the local coordinates.
!
! Frames from Stone and Alderton Mol Phys. 56, 5, 1047 (1985) plus weird
! modification for carbonyl to symmetrize.
!*******************************************************************************

subroutine local_to_global(crd, frames, lcl_crd, frame_cnt)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision      :: crd(3, *)
  type(ep_frame_rec)    :: frames(*)
  double precision      :: lcl_crd(3, 2, *)
  integer               :: frame_cnt
   
! Local variables:

  double precision      :: uvec(3)
  double precision      :: vvec(3)
  double precision      :: ave(3)
  double precision      :: diff(3)
  double precision      :: usiz
  double precision      :: vsiz
  double precision      :: asiz
  double precision      :: dsiz
  double precision      :: f(3, 3)
  double precision      :: a(3), b(3), c(3)
  integer               :: i, j, k, m
  integer               :: frame_id
#ifdef MPI
  integer               :: my_lst_idx
#endif
   
  if (frame_cnt .eq. 0) return

#ifdef MPI
  do my_lst_idx = 1, my_ep_frame_cnt
    frame_id = gbl_my_ep_frame_lst(my_lst_idx)
#else
  do frame_id = 1, frame_cnt
#endif

    if (frames(frame_id)%type .eq. 1) then
      do m = 1, 3
        a(m) = crd(m, frames(frame_id)%frame_atm1)
        b(m) = crd(m, frames(frame_id)%frame_atm2)
        c(m) = crd(m, frames(frame_id)%frame_atm3)
      end do
    else if (frames(frame_id)%type .eq. 2) then
      do m = 1, 3
        a(m) = 0.5d0 * (crd(m, frames(frame_id)%frame_atm1) + &
               crd(m, frames(frame_id)%frame_atm2))
        b(m) = crd(m, frames(frame_id)%parent_atm)
        c(m) = 0.5d0 * (crd(m, frames(frame_id)%frame_atm3) + &
               crd(m, frames(frame_id)%frame_atm2))
      end do
    end if
      
    ! Z-axis along symmmetry axis of b midway between
    ! unit vector to a and unit vector to c; points opposite:
      
    usiz = 0.d0
    vsiz = 0.d0
    do m = 1, 3
      uvec(m) = a(m) - b(m)
      usiz = usiz + uvec(m) * uvec(m)
      vvec(m) = c(m) - b(m)
      vsiz = vsiz + vvec(m) * vvec(m)
    end do
    usiz = sqrt(usiz)
    vsiz = sqrt(vsiz)
    asiz = 0.d0
    dsiz = 0.d0
    do m = 1, 3
      uvec(m) = uvec(m) / usiz
      vvec(m) = vvec(m) / vsiz
      ave(m) = (uvec(m) + vvec(m)) / 2.d0
      asiz = asiz + ave(m) * ave(m)
      diff(m) = (vvec(m) - uvec(m)) / 2.d0
      dsiz = dsiz + diff(m) * diff(m)
    end do
    asiz = sqrt(asiz)
    dsiz = sqrt(dsiz)
    do m = 1, 3
      f(m, 3) = -ave(m) / asiz
      f(m, 1) = diff(m) / dsiz
    end do
    f(1, 2) = f(2, 3) * f(3, 1) - f(3, 3) * f(2, 1)
    f(2, 2) = f(3, 3) * f(1, 1) - f(1, 3) * f(3, 1)
    f(3, 2) = f(1, 3) * f(2, 1) - f(2, 3) * f(1, 1)
    do k = 1, frames(frame_id)%ep_cnt
      j = frames(frame_id)%extra_pnt(k)
      do m = 1, 3
        crd(m, j) = crd(m, frames(frame_id)%parent_atm) + &
                    lcl_crd(1, k, frame_id) * f(m, 1) + &
                    lcl_crd(2, k, frame_id) * f(m, 2) + &
                    lcl_crd(3, k, frame_id) * f(m, 3)
      end do
    end do

  end do  !  frame_id = 1, frame_cnt

  return

end subroutine local_to_global 

!*******************************************************************************!
! Subroutine:   orient_frc
!
! Description:  Transfer forces from EP to main atoms in frame.
!
! Frames from Stone and Alderton Mol Phys. 56, 5, 1047 (1985)
!*******************************************************************************

subroutine orient_frc(crd, frc, framevir, frames, frame_cnt)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision      :: crd(3, *)
  double precision      :: frc(3, *)
  double precision      :: framevir(3, 3)
  type(ep_frame_rec)    :: frames(*)
  integer               :: frame_cnt
   
! Local variables:

  double precision      :: u(3), v(3), w(3)
  double precision      :: up(3), vp(3), diff(3)
  double precision      :: usiz, vsiz, wsiz, upsiz, vpsiz
  double precision      :: dotdu, dotdv, dphidu, dphidv, dphidw
  double precision      :: c, s, uvdis, vudis, du(3), dv(3)
  double precision      :: force(3), torque(3), rel(3)
  integer               :: i, j, k, l, j1, j2, j3, j4, m
  integer               :: frame_id
#ifdef MPI
  integer               :: my_lst_idx
#endif
  double precision      :: ap(3), bp(3), cp(3)
   
  ! Motions of the frame can be described in terms of rotations about the
  ! unit vectors u and v from frame_atm2 to frame_atm1, frame_atm3 respecively
  ! and rotation about the unit cross product w.
   
  if (frame_cnt .eq. 0) return
   
  framevir = 0.d0

#ifdef MPI
  do my_lst_idx = 1, my_ep_frame_cnt
    frame_id = gbl_my_ep_frame_lst(my_lst_idx)
#else
  do frame_id = 1, frame_cnt
#endif
    force  = 0.d0
    torque = 0.d0
    i = frames(frame_id)%parent_atm
    do k = 1, frames(frame_id)%ep_cnt
      j = frames(frame_id)%extra_pnt(k)
      force(1:3) = force(1:3) + frc(1:3, j)
      rel(1:3) = crd(1:3, j) - crd(1:3, i)
         
      ! Get transferred force component of virial:
         
      do m = 1, 3
        framevir(:, m) = framevir(:, m) + frc(:, j) * rel(m)
      end do
         
      ! Torque is rel x frc:
         
      torque(1) = torque(1) + rel(2) * frc(3, j) - rel(3) * frc(2, j)
      torque(2) = torque(2) + rel(3) * frc(1, j) - rel(1) * frc(3, j)
      torque(3) = torque(3) + rel(1) * frc(2, j) - rel(2) * frc(1, j)
      frc(:, j) = 0.d0
    end do
      
    if (frames(frame_id)%type .eq. 1) then
      do m = 1, 3
        ap(m) = crd(m, frames(frame_id)%frame_atm1)
        bp(m) = crd(m, frames(frame_id)%frame_atm2)
        cp(m) = crd(m, frames(frame_id)%frame_atm3)
      end do
    else if (frames(frame_id)%type .eq. 2) then
      do m = 1, 3
        ap(m) = 0.5d0 * (crd(m, frames(frame_id)%frame_atm1) + &
                crd(m, frames(frame_id)%frame_atm2))
        bp(m) = crd(m, frames(frame_id)%parent_atm)
        cp(m) = 0.5d0 * (crd(m, frames(frame_id)%frame_atm3) + &
                crd(m, frames(frame_id)%frame_atm2))
      end do
    end if
    usiz = 0.d0
    vsiz = 0.d0
    do m = 1, 3
      u(m) = ap(m) - bp(m)
      usiz = usiz + u(m) * u(m)
      v(m) = cp(m) - bp(m)
      vsiz = vsiz + v(m) * v(m)
    end do
    usiz = sqrt(usiz)
    vsiz = sqrt(vsiz)
    w(1) = u(2) * v(3) - u(3) * v(2)
    w(2) = u(3) * v(1) - u(1) * v(3)
    w(3) = u(1) * v(2) - u(2) * v(1)
    wsiz = sqrt(w(1) * w(1) + w(2) * w(2) + w(3) * w(3))
    dotdu = 0.d0
    dotdv = 0.d0
    do m = 1, 3
      u(m) = u(m) / usiz
      v(m) = v(m) / vsiz
      w(m) = w(m) / wsiz
      diff(m) = v(m) - u(m)
      dotdu = dotdu + u(m) * diff(m)
      dotdv = dotdv + v(m) * diff(m)
    end do
      
    ! Get perps to u, v to get direction of motion of u or v
    ! due to rotation about the cross product vector w:
      
    upsiz = 0.d0
    vpsiz = 0.d0
    do m = 1, 3
      up(m) = diff(m) - dotdu * u(m)
      vp(m) = diff(m) - dotdv * v(m)
      upsiz = upsiz + up(m) * up(m)
      vpsiz = vpsiz + vp(m) * vp(m)
    end do
    upsiz = sqrt(upsiz)
    vpsiz = sqrt(vpsiz)
    do m = 1, 3
      up(m) = up(m) / upsiz
      vp(m) = vp(m) / vpsiz
    end do
      
    ! Negative of dot product of torque with unit vectors
    ! along u, v and w.  Give result of infinitesmal rotation
    ! along these vectors, i.e. dphi / dtheta = dot product.
      
    dphidu = -(torque(1) * u(1) + torque(2) * u(2) + torque(3) * u(3))
    dphidv = -(torque(1) * v(1) + torque(2) * v(2) + torque(3) * v(3))
    dphidw = -(torque(1) * w(1) + torque(2) * w(2) + torque(3) * w(3))
      
    ! Get projected distances between vectors:
      
    c = u(1) * v(1) + u(2) * v(2) + u(3) * v(3)
    s = sqrt(1.d0 - c * c)
    uvdis = usiz * s
    vudis = vsiz * s
      
    !---------------------------------------------------------------------
    ! Frame formed by bisector of u, v, its perp, and w.
    ! movement of u by dz out of plane -> rotation about v of -dz / uvdis
    ! since positive rotation about v move u in negative dir. wrt w
    ! dphi / dz = dphi / dtheta dtheta / dz = -dotvt / uvdis
    ! movement of v by dz out of plane -> rotation about u of dz / vudis
    ! movement of u by dy along up -> rotation about w of 1/2 dy / usiz
    ! since bisector only rotates 1/2 as much as u or v in isolation
    ! movement of v by dy along vperp -> rotation about w of 1/2 dy / vsiz
    ! movement of u by dx along u doesn't change frame
    ! movement of v by dx along v doesn't change frame
    ! So... du_du = 0, du_dw = -dotvt / uvdis, du_dup = dotwt / (2.d0 * usiz)
    ! So... dv_dv = 0, dv_dw = dotut / vudis, du_dup = dotwt / (2.d0 * usiz)
    !---------------------------------------------------------------------
      
    if (frames(frame_id)%type .eq. 1) then
      j1 = frames(frame_id)%frame_atm1
      j2 = frames(frame_id)%frame_atm2
      j3 = frames(frame_id)%frame_atm3
      do m = 1, 3
        du(m) = -w(m) * dphidv / uvdis + up(m) * dphidw / (2.d0 * usiz)
        dv(m) = w(m) * dphidu / vudis + vp(m) * dphidw / (2.d0 * vsiz)
        frc(m, j1) = frc(m, j1) - du(m)
        frc(m, j3) = frc(m, j3) - dv(m)
        frc(m, j2) = frc(m, j2) + dv(m) + du(m) + force(m)
      end do
         
      ! Get torque contribution to virial:
         
      do m = 1, 3
        framevir(:, m) = framevir(:, m) + du(:) * (ap(m) - bp(m)) + &
                           dv(:) * (cp(m) - bp(m))
      end do
         
    else if (frames(frame_id)%type .eq. 2) then
      
      ! Need to transfer forces from midpoints to atoms:
         
      j1 = frames(frame_id)%frame_atm1
      j2 = frames(frame_id)%frame_atm2
      j3 = frames(frame_id)%frame_atm3
      j4 = frames(frame_id)%parent_atm
      do m = 1, 3
        du(m) = -w(m) * dphidv / uvdis + up(m) * dphidw / (2.d0 * usiz)
        dv(m) = w(m) * dphidu / vudis + vp(m) * dphidw / (2.d0 * vsiz)
        frc(m, j1) = frc(m, j1) - 0.5d0 * du(m)
        frc(m, j3) = frc(m, j3) - 0.5d0 * dv(m)
        frc(m, j2) = frc(m, j2) - 0.5d0 * (du(m) + dv(m))
        frc(m, j4) = frc(m, j4) + dv(m) + du(m) + force(m)
      end do
         
      ! Get torque contribution to virial:
         
      do m = 1, 3
        framevir(:, m) = framevir(:, m) + du(:) * (ap(m) - bp(m)) + &
                         dv(:) * (cp(m) - bp(m))
      end do
         
    end if  ! ( frames(frame_id)%type .eq. 1 )
      
    !---------------------------------------------------------------------
    ! OTHER TYPE FRAME; NOT SEEN YET
    ! Frame formed by  u, its perp, and w
    ! movement of v in plane doesn't change frame
    ! movement of u by dz out of plane -> rotation about v of -dz / uvdis
    ! since positive rotation about v move u in negative dir. wrt w
    ! dphi / dz = dphi / dtheta dtheta / dz = -dotvt / uvdis
    ! movement of v by dz out of plane -> rotation about u of dz / vudis
    ! movement of u by dy along up -> rotation about w of dy / usiz
    ! since frame rotates as much as u in isolation
    !---------------------------------------------------------------------
    !  do m = 1, 3
    !    du(m) = -w(m) * dphidv / uvdis + up(m) * dphidw / usiz
    !    dv(m) = w(m) * dphidu / vudis
    !    frc(m, j1) = frc(m, j1) - du(m) + force(m)
    !    frc(m, j3) = frc(m, j3) - dv(m)
    !    frc(m, j2) = frc(m, j2) + dv(m) + du(m)
    !  enddo
    !
    !  get torque contribution to virial:
    !
    !  do m = 1, 3
    !    do l = 1, 3
    !      framevir(l, m) = framevir(l, m) + &
    !                       du(l) * (crd(m, j1) - crd(m, j2)) + &
    !                       dv(l) * (crd(m, j3) - crd(m, j2))
    !    enddo
    !  enddo
      
  end do  !  frame_id = 1, frame_cnt

  return

end subroutine orient_frc 

!*******************************************************************************!
! Subroutine:   zero_extra_pnts_vec
!
! Description:  Set extra points vector entries to 0.0.
!
!*******************************************************************************

subroutine zero_extra_pnts_vec(vec, frames, frame_cnt)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision      :: vec(3, *)
  type(ep_frame_rec)    :: frames(*)
  integer               :: frame_cnt
   
! Local variables:

  integer               :: i
  integer               :: frame_id
#ifdef MPI
  integer               :: my_lst_idx
#endif

  if (frame_cnt .eq. 0) return
   
#ifdef MPI
  do my_lst_idx = 1, my_ep_frame_cnt
    frame_id = gbl_my_ep_frame_lst(my_lst_idx)
#else
  do frame_id = 1, frame_cnt
#endif
    do i = 1, frames(frame_id)%ep_cnt
      vec(:, frames(frame_id)%extra_pnt(i)) = 0.d0
    end do
  end do
  return

end subroutine zero_extra_pnts_vec 

!*******************************************************************************!
! Internal Subroutine:  get_nb14_energy
!
! Description: <TBS>
!
!*******************************************************************************

subroutine get_nb14_energy(charge, crd, frc, subfrc, iac, ico, cn1, cn2, &
                           nb14, nb14_cnt, ee14, enb14, e14vir)

  use prmtop_dat_mod, only : ntypes, gbl_one_scee, gbl_one_scnb
  use parallel_dat_mod
  use ti_mod

  implicit none

! Formal arguments:

  double precision, intent(in)                  :: charge(*)
  double precision, intent(in)                  :: crd(3, *)
  double precision, intent(in out)              :: frc(3, *)
! Added by Lijiang Yang
  double precision, intent(in out)              :: subfrc(3, *)
!----------------------------------------------------------------
  integer, intent(in)                           :: iac(*)
  integer, intent(in)                           :: ico(*)
  double precision, intent(in)                  :: cn1(*)
  double precision, intent(in)                  :: cn2(*)
  integer, intent(in)                           :: nb14(3, *)
  integer, intent(in)                           :: nb14_cnt
  double precision, intent(out)                 :: ee14
  double precision, intent(out)                 :: enb14
  double precision, optional, intent(out)       :: e14vir(3, 3)
   
! Local variables:

  integer               :: n, i, j, ic, parm_idx, ia1, ia2, ibig, isml
  integer               :: ntypes_lcl
  logical               :: do_virial
  double precision      :: dx, dy, dz, r2, r2inv, rinv
  double precision      :: scnb0, g, f6, f12, r6, df

  ee14 = 0.d0
  enb14 = 0.d0

  do_virial = present(e14vir)

  if (do_virial) e14vir(:, :) = 0.d0

  if (nb14_cnt .eq. 0) return

  ntypes_lcl = ntypes

  do n = 1, nb14_cnt
    i = nb14(1, n)
    j = nb14(2, n)
    parm_idx = nb14(3, n)
    scnb0 = gbl_one_scnb(parm_idx)
    dx = crd(1, j) - crd(1, i)
    dy = crd(2, j) - crd(2, i)
    dz = crd(3, j) - crd(3, i)
    r2 = dx * dx + dy * dy + dz * dz
    rinv = sqrt(1.d0 / r2)
    r2inv = rinv * rinv
    r6 = r2inv * r2inv * r2inv
    g = charge(i) * charge(j) * rinv * gbl_one_scee(parm_idx)
    !  always use the 6-12 parameters, even if 10-12 are available:
    ia1 = iac(i)
    ia2 = iac(j)
    ibig = max0(ia1,ia2)
    isml = min0(ia1,ia2)
    ic = ibig*(ibig-1)/2+isml
    f6 = cn2(ic) * r6
    f12 = cn1(ic) * (r6 * r6)
    df = (g + scnb0 * (12.d0 * f12 - 6.d0 * f6)) * r2inv

    if (ti_mode .ne. 0) then
      ti_region = 0
      if (ti_lst(1,i) .eq. 1 .or. ti_lst(1,j) .eq. 1) then
        ti_region = 1
      else if(ti_lst(2,i) .eq. 1 .or. ti_lst(2,j) .eq. 1) then
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
          if(ti_mode .eq. 1 .or. ti_lscale) then !linear scaling
            call ti_update_ene(g, si_elect_14_ene, ti_region)
            call ti_update_ene((f12 - f6) * scnb0, si_vdw_14_ene, ti_region)
            g = g * ti_weights(ti_region)    
            f12 = f12 * ti_weights(ti_region)
            f6 = f6 * ti_weights(ti_region)       
          else
            ti_ene(ti_region,si_vdw_14_ene) = &
               ti_ene(ti_region,si_vdw_14_ene) + (f12 - f6) * scnb0
            ti_ene(ti_region,si_elect_14_ene) = &
               ti_ene(ti_region,si_elect_14_ene) + g
            f6 = 0.d0
            f12 = 0.d0
            g = 0.d0
          end if 
#ifdef MPI
        end if
#endif
        if (ti_mode .eq. 1 .or. ti_lscale) df = df * ti_weights(ti_region) 
      end if
    end if

#ifdef MPI
! We use i to determine who sums up the energy...

    if (gbl_atm_owner_map(i) .eq. mytaskid) then

      ee14 = ee14 + g
      enb14 = enb14 + (f12 - f6) * scnb0

      frc(1, i) = frc(1, i) - df * dx
      frc(2, i) = frc(2, i) - df * dy
      frc(3, i) = frc(3, i) - df * dz
      subfrc(1, i) = - df * dx
      subfrc(2, i) = - df * dy
      subfrc(3, i) = - df * dz

      if (do_virial) then
        e14vir(1, 1) = e14vir(1, 1) - df * dx * dx
        e14vir(1, 2) = e14vir(1, 2) - df * dx * dy
        e14vir(1, 3) = e14vir(1, 3) - df * dx * dz
        e14vir(2, 2) = e14vir(2, 2) - df * dy * dy
        e14vir(2, 3) = e14vir(2, 3) - df * dy * dz
        e14vir(3, 3) = e14vir(3, 3) - df * dz * dz
      end if

      if (ti_mode .ne. 0) then
        call ti_update_nb_frc(-df * dx, -df * dy ,-df * dz, i, ti_region)
      end if

    end if

    if (gbl_atm_owner_map(j) .eq. mytaskid) then

      frc(1, j) = frc(1, j) + df * dx
      frc(2, j) = frc(2, j) + df * dy
      frc(3, j) = frc(3, j) + df * dz
      subfrc(1, j) = df * dx
      subfrc(2, j) = df * dy
      subfrc(3, j) = df * dz

      if (ti_mode .ne. 0) then
        call ti_update_nb_frc(df * dx, df * dy, df * dz, j, ti_region)
      end if

    end if
     
#else

    ee14 = ee14 + g
    enb14 = enb14 + (f12 - f6) * scnb0
    frc(1, i) = frc(1, i) - df * dx
    frc(2, i) = frc(2, i) - df * dy
    frc(3, i) = frc(3, i) - df * dz
    frc(1, j) = frc(1, j) + df * dx
    frc(2, j) = frc(2, j) + df * dy
    frc(3, j) = frc(3, j) + df * dz
    subfrc(1, i) = - df * dx
    subfrc(2, i) = - df * dy
    subfrc(3, i) = - df * dz
    subfrc(1, j) = df * dx
    subfrc(2, j) = df * dy
    subfrc(3, j) = df * dz

    if (do_virial) then
      e14vir(1, 1) = e14vir(1, 1) - df * dx * dx
      e14vir(1, 2) = e14vir(1, 2) - df * dx * dy
      e14vir(1, 3) = e14vir(1, 3) - df * dx * dz
      e14vir(2, 2) = e14vir(2, 2) - df * dy * dy
      e14vir(2, 3) = e14vir(2, 3) - df * dy * dz
      e14vir(3, 3) = e14vir(3, 3) - df * dz * dz
    end if

    if (ti_mode .ne. 0) then
      call ti_update_nb_frc(-df * dx, -df * dy, -df * dz, i, ti_region)
      call ti_update_nb_frc(df * dx, df * dy, df * dz, j, ti_region)
    end if

#endif /* MPI */
  end do

  if (do_virial) then
    e14vir(2, 1) = e14vir(1, 2)
    e14vir(3, 1) = e14vir(1, 3)
    e14vir(3, 2) = e14vir(2, 3)
  end if

  return

end subroutine get_nb14_energy 

!*******************************************************************************!
! Internal Subroutine:  get_nb14_energy_sub
!
! Description: <TBS>
!
!*******************************************************************************

subroutine get_nb14_energy_sub(charge, crd, subfrc, iac, ico, cn1, cn2, &
                           nb14, nb14_cnt, ee14, enb14, e14vir, first)

  use prmtop_dat_mod, only : ntypes, gbl_one_scee, gbl_one_scnb
  use parallel_dat_mod
  use ti_mod

  implicit none

! Formal arguments:

  double precision, intent(in)                  :: charge(*)
  double precision, intent(in)                  :: crd(3, *)
! Added by Lijiang Yang
  double precision, intent(in out)              :: subfrc(3, *)
!----------------------------------------------------------------
  integer, intent(in)                           :: iac(*)
  integer, intent(in)                           :: ico(*)
  double precision, intent(in)                  :: cn1(*)
  double precision, intent(in)                  :: cn2(*)
  integer, intent(in)                           :: nb14(3, *)
  integer, intent(in)                           :: nb14_cnt
  double precision, intent(out)                 :: ee14
  double precision, intent(out)                 :: enb14
  double precision, optional, intent(out)       :: e14vir(3, 3)
! Added by Lijiang Yang
  logical, optional                             :: first
!----------------------------------------------------------------
   
! Local variables:

  integer               :: n, i, j, ic, parm_idx, ia1, ia2, ibig, isml
  integer               :: ntypes_lcl
  logical               :: do_virial, firstcall, lastcall
  double precision      :: dx, dy, dz, r2, r2inv, rinv
  double precision      :: scnb0, g, f6, f12, r6, df

  ee14 = 0.d0
  enb14 = 0.d0

  do_virial = present(e14vir)
  firstcall = present(first)

  if (do_virial .and. firstcall) e14vir(:, :) = 0.d0

  if (nb14_cnt .eq. 0) return

  ntypes_lcl = ntypes

  do n = 1, nb14_cnt
    i = nb14(1, n)
    j = nb14(2, n)
    parm_idx = nb14(3, n)
    scnb0 = gbl_one_scnb(parm_idx)
    dx = crd(1, j) - crd(1, i)
    dy = crd(2, j) - crd(2, i)
    dz = crd(3, j) - crd(3, i)
    r2 = dx * dx + dy * dy + dz * dz
    rinv = sqrt(1.d0 / r2)
    r2inv = rinv * rinv
    r6 = r2inv * r2inv * r2inv
    g = charge(i) * charge(j) * rinv * gbl_one_scee(parm_idx)
    !  always use the 6-12 parameters, even if 10-12 are available:
    ia1 = iac(i)
    ia2 = iac(j)
    ibig = max0(ia1,ia2)
    isml = min0(ia1,ia2)
    ic = ibig*(ibig-1)/2+isml
    f6 = cn2(ic) * r6
    f12 = cn1(ic) * (r6 * r6)
    df = (g + scnb0 * (12.d0 * f12 - 6.d0 * f6)) * r2inv

    if (ti_mode .ne. 0) then
      ti_region = 0
      if (ti_lst(1,i) .eq. 1 .or. ti_lst(1,j) .eq. 1) then
        ti_region = 1
      else if(ti_lst(2,i) .eq. 1 .or. ti_lst(2,j) .eq. 1) then
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
          if(ti_mode .eq. 1 .or. ti_lscale) then !linear scaling
            call ti_update_ene(g, si_elect_14_ene, ti_region)
            call ti_update_ene((f12 - f6) * scnb0, si_vdw_14_ene, ti_region)
            g = g * ti_weights(ti_region)    
            f12 = f12 * ti_weights(ti_region)
            f6 = f6 * ti_weights(ti_region)       
          else
            ti_ene(ti_region,si_vdw_14_ene) = &
               ti_ene(ti_region,si_vdw_14_ene) + (f12 - f6) * scnb0
            ti_ene(ti_region,si_elect_14_ene) = &
               ti_ene(ti_region,si_elect_14_ene) + g
            f6 = 0.d0
            f12 = 0.d0
            g = 0.d0
          end if 
#ifdef MPI
        end if
#endif
        if (ti_mode .eq. 1 .or. ti_lscale) df = df * ti_weights(ti_region) 
      end if
    end if

#ifdef MPI
! We use i to determine who sums up the energy...

    if (gbl_atm_owner_map(i) .eq. mytaskid) then

      ee14 = ee14 + g
      enb14 = enb14 + (f12 - f6) * scnb0

      subfrc(1, i) = subfrc(1, i) - df * dx
      subfrc(2, i) = subfrc(2, i) - df * dy
      subfrc(3, i) = subfrc(3, i) - df * dz

      if (do_virial) then
        e14vir(1, 1) = e14vir(1, 1) - df * dx * dx
        e14vir(1, 2) = e14vir(1, 2) - df * dx * dy
        e14vir(1, 3) = e14vir(1, 3) - df * dx * dz
        e14vir(2, 2) = e14vir(2, 2) - df * dy * dy
        e14vir(2, 3) = e14vir(2, 3) - df * dy * dz
        e14vir(3, 3) = e14vir(3, 3) - df * dz * dz
      end if

      if (ti_mode .ne. 0) then
        call ti_update_nb_frc(-df * dx, -df * dy ,-df * dz, i, ti_region)
      end if

    end if

    if (gbl_atm_owner_map(j) .eq. mytaskid) then

      subfrc(1, j) = subfrc(1, j) + df * dx
      subfrc(2, j) = subfrc(2, j) + df * dy
      subfrc(3, j) = subfrc(3, j) + df * dz

      if (ti_mode .ne. 0) then
        call ti_update_nb_frc(df * dx, df * dy, df * dz, j, ti_region)
      end if

    end if
     
#else

    ee14 = ee14 + g
    enb14 = enb14 + (f12 - f6) * scnb0
    subfrc(1, i) = subfrc(1, i) - df * dx
    subfrc(2, i) = subfrc(2, i) - df * dy
    subfrc(3, i) = subfrc(3, i) - df * dz
    subfrc(1, j) = subfrc(1, j) + df * dx
    subfrc(2, j) = subfrc(2, j) + df * dy
    subfrc(3, j) = subfrc(3, j) + df * dz

    if (do_virial) then
      e14vir(1, 1) = e14vir(1, 1) - df * dx * dx
      e14vir(1, 2) = e14vir(1, 2) - df * dx * dy
      e14vir(1, 3) = e14vir(1, 3) - df * dx * dz
      e14vir(2, 2) = e14vir(2, 2) - df * dy * dy
      e14vir(2, 3) = e14vir(2, 3) - df * dy * dz
      e14vir(3, 3) = e14vir(3, 3) - df * dz * dz
    end if

    if (ti_mode .ne. 0) then
      call ti_update_nb_frc(-df * dx, -df * dy, -df * dz, i, ti_region)
      call ti_update_nb_frc(df * dx, df * dy, df * dz, j, ti_region)
    end if

#endif /* MPI */
  end do

!  if (do_virial) then
!    e14vir(2, 1) = e14vir(1, 2)
!    e14vir(3, 1) = e14vir(1, 3)
!    e14vir(3, 2) = e14vir(2, 3)
!  end if

  return

end subroutine get_nb14_energy_sub 

!*******************************************************************************!
! Internal Subroutine:  get_nb14_energy_inter
!
! Description: <TBS>
!
!*******************************************************************************

subroutine get_nb14_energy_intervir(charge, crd, subfrc, iac, ico, cn1, cn2, &
                           nb14, nb14_cnt, ee14, enb14, e14vir)

  use prmtop_dat_mod, only : ntypes, gbl_one_scee, gbl_one_scnb
  use parallel_dat_mod
  use ti_mod

  implicit none

! Formal arguments:

  double precision, intent(in)                  :: charge(*)
  double precision, intent(in)                  :: crd(3, *)
! Added by Lijiang Yang
  double precision, intent(in out)              :: subfrc(3, *)
!----------------------------------------------------------------
  integer, intent(in)                           :: iac(*)
  integer, intent(in)                           :: ico(*)
  double precision, intent(in)                  :: cn1(*)
  double precision, intent(in)                  :: cn2(*)
  integer, intent(in)                           :: nb14(3, *)
  integer, intent(in)                           :: nb14_cnt
  double precision, intent(out)                 :: ee14
  double precision, intent(out)                 :: enb14
  double precision, optional, intent(out)       :: e14vir(3, 3)
   
! Local variables:

  integer               :: n, i, j, ic, parm_idx, ia1, ia2, ibig, isml
  integer               :: ntypes_lcl
  double precision      :: dx, dy, dz, r2, r2inv, rinv
  double precision      :: scnb0, g, f6, f12, r6, df

  ee14 = 0.d0
  enb14 = 0.d0


  if (nb14_cnt .eq. 0) return

  ntypes_lcl = ntypes

  do n = 1, nb14_cnt
    i = nb14(1, n)
    j = nb14(2, n)
    parm_idx = nb14(3, n)
    scnb0 = gbl_one_scnb(parm_idx)
    dx = crd(1, j) - crd(1, i)
    dy = crd(2, j) - crd(2, i)
    dz = crd(3, j) - crd(3, i)
    r2 = dx * dx + dy * dy + dz * dz
    rinv = sqrt(1.d0 / r2)
    r2inv = rinv * rinv
    r6 = r2inv * r2inv * r2inv
    g = charge(i) * charge(j) * rinv * gbl_one_scee(parm_idx)
    !  always use the 6-12 parameters, even if 10-12 are available:
    ia1 = iac(i)
    ia2 = iac(j)
    ibig = max0(ia1,ia2)
    isml = min0(ia1,ia2)
    ic = ibig*(ibig-1)/2+isml
    f6 = cn2(ic) * r6
    f12 = cn1(ic) * (r6 * r6)
    df = (g + scnb0 * (12.d0 * f12 - 6.d0 * f6)) * r2inv

    if (ti_mode .ne. 0) then
      ti_region = 0
      if (ti_lst(1,i) .eq. 1 .or. ti_lst(1,j) .eq. 1) then
        ti_region = 1
      else if(ti_lst(2,i) .eq. 1 .or. ti_lst(2,j) .eq. 1) then
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
          if(ti_mode .eq. 1 .or. ti_lscale) then !linear scaling
            call ti_update_ene(g, si_elect_14_ene, ti_region)
            call ti_update_ene((f12 - f6) * scnb0, si_vdw_14_ene, ti_region)
            g = g * ti_weights(ti_region)    
            f12 = f12 * ti_weights(ti_region)
            f6 = f6 * ti_weights(ti_region)       
          else
            ti_ene(ti_region,si_vdw_14_ene) = &
               ti_ene(ti_region,si_vdw_14_ene) + (f12 - f6) * scnb0
            ti_ene(ti_region,si_elect_14_ene) = &
               ti_ene(ti_region,si_elect_14_ene) + g
            f6 = 0.d0
            f12 = 0.d0
            g = 0.d0
          end if 
#ifdef MPI
        end if
#endif
        if (ti_mode .eq. 1 .or. ti_lscale) df = df * ti_weights(ti_region) 
      end if
    end if

#ifdef MPI
! We use i to determine who sums up the energy...

    if (gbl_atm_owner_map(i) .eq. mytaskid) then

      ee14 = ee14 + g
      enb14 = enb14 + (f12 - f6) * scnb0

      subfrc(1, i) = subfrc(1, i) - df * dx
      subfrc(2, i) = subfrc(2, i) - df * dy
      subfrc(3, i) = subfrc(3, i) - df * dz

      e14vir(1, 1) = e14vir(1, 1) - df * dx * dx
      e14vir(1, 2) = e14vir(1, 2) - df * dx * dy
      e14vir(1, 3) = e14vir(1, 3) - df * dx * dz
      e14vir(2, 2) = e14vir(2, 2) - df * dy * dy
      e14vir(2, 3) = e14vir(2, 3) - df * dy * dz
      e14vir(3, 3) = e14vir(3, 3) - df * dz * dz

      if (ti_mode .ne. 0) then
        call ti_update_nb_frc(-df * dx, -df * dy ,-df * dz, i, ti_region)
      end if

    end if

    if (gbl_atm_owner_map(j) .eq. mytaskid) then

      subfrc(1, j) = subfrc(1, j) + df * dx
      subfrc(2, j) = subfrc(2, j) + df * dy
      subfrc(3, j) = subfrc(3, j) + df * dz

      if (ti_mode .ne. 0) then
        call ti_update_nb_frc(df * dx, df * dy, df * dz, j, ti_region)
      end if

    end if
     
#else

    ee14 = ee14 + g
    enb14 = enb14 + (f12 - f6) * scnb0
    subfrc(1, i) = subfrc(1, i) - df * dx
    subfrc(2, i) = subfrc(2, i) - df * dy
    subfrc(3, i) = subfrc(3, i) - df * dz
    subfrc(1, j) = subfrc(1, j) + df * dx
    subfrc(2, j) = subfrc(2, j) + df * dy
    subfrc(3, j) = subfrc(3, j) + df * dz

    e14vir(1, 1) = e14vir(1, 1) - df * dx * dx
    e14vir(1, 2) = e14vir(1, 2) - df * dx * dy
    e14vir(1, 3) = e14vir(1, 3) - df * dx * dz
    e14vir(2, 2) = e14vir(2, 2) - df * dy * dy
    e14vir(2, 3) = e14vir(2, 3) - df * dy * dz
    e14vir(3, 3) = e14vir(3, 3) - df * dz * dz

    if (ti_mode .ne. 0) then
      call ti_update_nb_frc(-df * dx, -df * dy, -df * dz, i, ti_region)
      call ti_update_nb_frc(df * dx, df * dy, df * dz, j, ti_region)
    end if

#endif /* MPI */
  end do

    e14vir(2, 1) = e14vir(1, 2)
    e14vir(3, 1) = e14vir(1, 3)
    e14vir(3, 2) = e14vir(2, 3)

  return

end subroutine get_nb14_energy_intervir 

!*******************************************************************************!
! Subroutine:  fix_masses
!
! Description: We just fix up the atom masses; everything else is derived later
!              from them...
!
!*******************************************************************************

subroutine fix_masses(atm_cnt, amass, epowner)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: amass(*)
  integer               :: epowner(*)

! Local variables:

  integer               :: n
   
  ! Zero out mass for extra points;
   
  do n = 1, atm_cnt
    if (epowner(n) .ne. 0) then
      amass(n) = 0.d0
    end if
  end do
   
  return

end subroutine fix_masses 

end module extra_pnts_nb14_mod
