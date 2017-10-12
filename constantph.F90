!*******************************************************************************
!
! Module: constantph_mod
!
! Description: This module houses the constant pH functionality initially
!              implemented in sander in implicit solvent by John Mongan. 
!              Implemented here by Jason Swails
!
!*******************************************************************************

module constantph_mod

use constantph_dat_mod, only : on_cpstep, proposed_qterm
use file_io_mod,     only : amopen
use file_io_dat_mod, only : cpin, cpin_name, cpout, cpout_name, &
                            cprestrt, cprestrt_name, mdout, owrite
use random_mod,      only : random_state, amrand_gen, amrset_gen

implicit none

private ! everything here is private

! Limits on the sizes of the data structures
#define MAX_TITR_RES 50
#define MAX_TITR_STATES 200
#define MAX_ATOM_CHRG 1000
#define MAX_H_COUNT 4

! Constant pH variables

type :: cnstph_info
  sequence
  integer :: num_states
  integer :: first_atom
  integer :: num_atoms
  integer :: first_state
  integer :: first_charge
end type cnstph_info
integer, parameter :: SIZE_CNSTPH_INFO = 5

type(cnstph_info), parameter :: NULL_CPH_INFO = cnstph_info(0,0,0,0,0)
type(random_state), save     :: cnstph_randgen

double precision, save, allocatable :: gbl_chrgdat(:)
double precision, save, allocatable :: gbl_statene(:)

character (len=40), save, allocatable :: gbl_resname(:)

integer, save, allocatable :: gbl_resstate(:)
integer, save, allocatable :: gbl_protcnt(:)
integer, save, allocatable :: mobile_atoms(:)
integer, save              :: trescnt

integer, save, allocatable :: hindex(:)
integer, save, allocatable :: tpair(:)
integer, save, allocatable :: iselres(:)
integer, save, allocatable :: iselstat(:)

integer, save :: cphfirst_sol
integer, save, public :: cph_igb

double precision, save :: cph_intdiel ! not implemented

logical, save              :: first_cpout ! first pass through cpout file
logical, save              :: first_cprestrt ! first pass through cprestrt file
logical, save              :: cph_success ! to keep track of whether or not a
                                          ! protonation state changed

type(cnstph_info), save, allocatable :: gbl_stateinf(:)

#ifdef MPI
public cnstph_setup, cnstph_update_pairs, cnstph_begin_step, cnstph_end_step, &
       cnstph_write_restart, cnstph_write_cpout, cnstph_read, cnstph_bcast, &
       total_protonation, explicit_cnstph_begin_step, cnstph_explicitmd
#else
public cnstph_setup, cnstph_update_pairs, cnstph_begin_step, cnstph_end_step, &
       cnstph_write_restart, cnstph_write_cpout, cnstph_read, &
       total_protonation, explicit_cnstph_begin_step, cnstph_explicitmd
#endif

! tpair and hindex arrays: 
!
!  These arrays contain index information for each titr. res., and have similar
!  structures. The head of each array is a pointer section, from indices 0 to 
!  trescnt.  The residue number (0 to trescnt-1) is an index that shows the 
!  location of the FIRST element belonging to that residue (be it neighbors or
!  titrating protons). You can get the *number* of elements belonging to residue
!  "x" by subtracting array(x+1) - array(x). The trescnt-th element of each
!  array is a tail index so that trick also works for the last titr. res.
!  Let me show you an example using hindex: Consider 2 titr. res.  The 0th res
!  has titrating proton numbers 12, 15, 16, and 17. The 1st res has titrating
!  proton numbers 22 and 23. The hindex array would then be arranged as follows:
!  hindex = 3, 7, 9, 12, 15, 16, 17, 22, 23

contains

!*******************************************************************************
!
! Subroutine: cnstph_read
!
! Description: Reads the constant pH input file and initializes the data
!
!*******************************************************************************

subroutine cnstph_read(num_ints, num_reals)

  use gbl_constants_mod, only : AMBER_ELECTROSTATIC
  use mdin_ctrl_dat_mod
  use pmemd_lib_mod,     only : mexit, get_atomic_number
  use prmtop_dat_mod, only    : natom, atm_gb_fs, atm_atomicnumber, &
            loaded_atm_atomicnumber, atm_igraph, atm_mass, atm_gb_radii

  implicit none

! Passed variables
  
  integer, intent(in out) :: num_ints
  integer, intent(in out) :: num_reals

! Local variables

  ! Namelist variables -- same as gbl_ counterparts above

  double precision  :: chrgdat(0:MAX_ATOM_CHRG-1)
  double precision  :: statene(0:MAX_TITR_STATES-1)
  type(cnstph_info) :: stateinf(0:MAX_TITR_RES-1)

  integer       :: protcnt(0:MAX_TITR_STATES-1)
  integer       :: resstate(0:MAX_TITR_RES-1)
  character(40) :: resname(0:MAX_TITR_RES)

  integer       :: itres        ! residue iterator
  integer       :: itatm        ! atom iterator
  integer       :: last_res     ! last residue (test for overflow)
  integer       :: last_charge  ! last charge  (test for overflow)
  integer       :: last_state   ! last state   (test for overflow)
  integer       :: atomicnumber ! holder for atomic numbers
  integer       :: i            ! counter

  namelist /cnstph/ stateinf, resstate, protcnt, chrgdat, statene, &
                    trescnt, resname, cphfirst_sol, cph_igb, cph_intdiel
  
  ! Allocate the global constant pH data arrays

  call cnstph_allocate(num_ints, num_reals)

  ! Zero out our temporary arrays

  call cnstph_zero(stateinf, protcnt, resname, resstate, chrgdat, statene)

  ! Open the cpin file and read the cnstph namelist

  write(mdout, '(a,a)') '|reading charge increments from file: ', cpin_name

  call amopen(cpin, cpin_name, 'O', 'F', 'R')

  ! Initialize the internal dielectric in case it's not specified

  cph_intdiel = 1.0d0

  read (cpin, nml=cnstph)

  ! Scale the charges to amber internal units

  do itatm = 0, MAX_ATOM_CHRG - 1
    chrgdat(itatm) = chrgdat(itatm) * AMBER_ELECTROSTATIC
  end do

  ! Check to see if any of our fields are overflowing...
  
  if (trescnt .gt. MAX_TITR_RES) then
    write(mdout, '(a)') 'Too many titrating residues. Alter MAX_TITR_RES in &
                        &constantph.fpp and recompile'
    call mexit(mdout, 1)
  end if

  do itres = 0, trescnt - 1
    
    if (stateinf(itres)%first_state + stateinf(itres)%num_states .ge. &
        MAX_TITR_STATES) then
      write(mdout, '(a)') 'Too many titrating states; alter MAX_TITR_STATES in &
                          &constantph.fpp and recompile'
      call mexit(mdout, 1)
    end if

    if (stateinf(itres)%first_charge + stateinf(itres)%num_atoms .ge. &
        MAX_ATOM_CHRG) then
      write(mdout, '(a)') 'Too much charge data; alter MAX_ATOM_CHRG in &
                          &constantph.fpp and recompile'
      call mexit(mdout, 1)
    end if

  end do

  ! Copy the data into the main arrays
  
  gbl_stateinf(:) = stateinf(:)
  gbl_resstate(:) = resstate(:)
  gbl_protcnt(:)  = protcnt(:)
  gbl_chrgdat(:)  = chrgdat(:)
  gbl_statene(:)  = statene(:)
  gbl_resname(:)  = resname(:)

  ! See if we need to set up GB parameters

  if (icnstph .eq. 2) then
    gb_cutoff = 1000.d0  ! GB cutoff is infinite here
    if (cph_igb .eq. 1) then
      gb_alpha = 1.d0
      gb_beta = 0.d0
      gb_gamma = 0.d0
    else if (cph_igb .eq. 2) then
      ! Use our best guesses for Onufriev/Case GB  (GB^OBC I):
      gb_alpha = 0.8d0
      gb_beta = 0.d0
      gb_gamma = 2.909125d0
    else if (cph_igb .eq. 5) then
      ! Use our second best guesses for Onufriev/Case GB (GB^OBC II):
      gb_alpha = 1.d0
      gb_beta = 0.8d0
      gb_gamma = 4.85d0
    else if (cph_igb .eq. 7) then
      ! Use parameters for Mongan et al. CFA GBNECK:
      gb_alpha = 1.09511284d0
      gb_beta = 1.90792938d0
      gb_gamma = 2.50798245d0
      gb_neckscale = 0.361825d0
    else if (cph_igb .eq. 8) then
      gb_alpha_h   = 0.788440d0
      gb_beta_h    = 0.798699d0
      gb_gamma_h   = 0.437334d0
      gb_alpha_c   = 0.733756d0
      gb_beta_c    = 0.506378d0
      gb_gamma_c   = 0.205844d0
      gb_alpha_n   = 0.503364d0
      gb_beta_n    = 0.316828d0
      gb_gamma_n   = 0.192915d0
      gb_alpha_os  = 0.867814d0
      gb_beta_os   = 0.876635d0
      gb_gamma_os  = 0.387882d0
      gb_alpha_p   = 1.0d0
      gb_beta_p    = 0.8d0
      gb_gamma_p   = 4.85d0  
      gb_neckscale = 0.826836d0
      screen_h = 1.425952d0
      screen_c = 1.058554d0
      screen_n = 0.733599d0
      screen_o = 1.061039d0
      screen_s = -0.703469d0
      screen_p = 0.5d0
      offset  = 0.195141d0
    end if

    ! Set gb_kappa as long as saltcon is .ge. 0.d0 (and catch bug below if not).

    if (saltcon .ge. 0.d0) then

      ! Get Debye-Huckel kappa (A**-1) from salt concentration (M), assuming:
      !   T = 298.15, epsext=78.5,
      gb_kappa = sqrt(0.10806d0 * saltcon)

      ! Scale kappa by 0.73 to account(?) for lack of ion exclusions:

      gb_kappa = 0.73d0 * gb_kappa

    end if

    ! We need to set up the data structures for igb == 7 or igb == 8 since this
    ! is done in subroutine init_prmtop_dat by default, and that was already
    ! called so we would know the value of 'natom'

    if (cph_igb .eq. 7) then

      write(mdout,'(a)') &
        ' Replacing prmtop screening parameters with GBn (igb=7) values'

      do i = 1, natom

        if (loaded_atm_atomicnumber) then
          atomicnumber = atm_atomicnumber(i)
        else
          call get_atomic_number(atm_igraph(i), atm_mass(i), atomicnumber)
        end if

        if (atomicnumber .eq. 6) then
          atm_gb_fs(i) = 4.84353823306d-1
        else if (atomicnumber .eq. 1) then
          atm_gb_fs(i) = 1.09085413633d0
        else if (atomicnumber .eq. 7) then
          atm_gb_fs(i) = 7.00147318409d-1
        else if (atomicnumber .eq. 8) then
          atm_gb_fs(i) = 1.06557401132d0
        else if (atomicnumber .eq. 16) then
          atm_gb_fs(i) = 6.02256336067d-1
        else
          atm_gb_fs(i) = 5.d-1 ! not optimized
        end if

      end do
    
    else if (cph_igb .eq. 8) then
      
      write(mdout, '(a)') &
        ' Replacing prmtop screening parameters with GBn2 (igb=8) values'

      do i = 1, natom

        if(loaded_atm_atomicnumber) then
          atomicnumber = atm_atomicnumber(i)
        else
          call get_atomic_number(atm_igraph(i), atm_mass(i), atomicnumber)
        end if

        ! The screen_ variables are found in mdin_ctrl_dat_mod

        if (atomicnumber .eq. 6) then
          atm_gb_fs(i) = screen_c
        else if (atomicnumber .eq. 1) then
          atm_gb_fs(i) = screen_h
        else if (atomicnumber .eq. 7) then
          atm_gb_fs(i) = screen_n
        else if (atomicnumber .eq. 8) then
          atm_gb_fs(i) = screen_o
        else if (atomicnumber .eq. 16) then
          atm_gb_fs(i) = screen_s
        else if (atomicnumber .eq. 15) then
          atm_gb_fs(i) = screen_p
        else
          atm_gb_fs(i) = 5.d-1 ! not optimized
        end if

      end do

    end if ! cph_igb .eq. 7

    if (cph_igb .eq. 7 .or. cph_igb .eq. 8) then

      ! In this case, we need to re-compute atm_gb_fs(i), since we changed its
      ! value above

      gb_fs_max = 0.d0

      do i = 1, natom
        atm_gb_fs(i) = atm_gb_fs(i) * (atm_gb_radii(i) - offset)
        gb_fs_max = max(gb_fs_max, atm_gb_fs(i))
      end do

    end if

    mobile_atoms(1:cphfirst_sol-1) = 0
    mobile_atoms(cphfirst_sol:natom) = 1
  end if ! (icnstph .eq. 2)

  return

end subroutine cnstph_read

!*******************************************************************************
!
! Subroutine: cnstph_allocate
!
! Description: Allocates constant pH data structures
!
!*******************************************************************************

subroutine cnstph_allocate(num_ints, num_reals)

  use constantph_dat_mod, only: allocate_cnstph_dat
  use mdin_ctrl_dat_mod, only : icnstph
  use pmemd_lib_mod,  only    : mexit
  use prmtop_dat_mod, only    : natom

  implicit none

! Passed variables

  integer, intent(in out) :: num_ints
  integer, intent(in out) :: num_reals

! Local variables

  integer :: alloc_failed

  allocate(gbl_chrgdat(0:MAX_ATOM_CHRG-1),    &
           gbl_statene(0:MAX_TITR_STATES-1),  &
           gbl_stateinf(0:MAX_TITR_RES-1),    &
           gbl_protcnt(0:MAX_TITR_STATES-1),  &
           gbl_resname(0:MAX_TITR_RES),       &
           gbl_resstate(0:MAX_TITR_RES-1),    &
           hindex(0:MAX_TITR_RES * (MAX_H_COUNT+1) - 1), &
           tpair(0:MAX_TITR_RES * 4),         &
           iselres(0:MAX_TITR_RES),           &
           iselstat(2),                       &
           stat=alloc_failed )

  if (alloc_failed .ne. 0) then
    write(6,'(a)') 'Error in constant pH allocation!'
    call mexit(6, 1)
  end if

  call allocate_cnstph_dat(natom, num_reals, alloc_failed)
  if (alloc_failed .ne. 0) then
    write(6,'(a)') 'Error in constant pH allocation!'
    call mexit(6, 1)
  end if

  ! Allocate our mobile_atoms array

  if (icnstph .eq. 2) then
    allocate(mobile_atoms(natom), stat=alloc_failed)
    if (alloc_failed .ne. 0) then
      write(6,'(a)') 'Error in constant pH allocation!'
      call mexit(6, 1)
    end if
    num_ints = num_ints + size(mobile_atoms)
  end if

  num_ints = num_ints + size(gbl_stateinf)   &
                      + size(gbl_protcnt)    &
                      + size(gbl_resname)    &
                      + size(gbl_resstate)   &
                      + size(hindex)         &
                      + size(tpair)          &
                      + size(iselres)        &
                      + size(iselstat)

  num_reals = num_reals + size(gbl_chrgdat)  &
                        + size(gbl_statene)

end subroutine cnstph_allocate

!*******************************************************************************
!
! Subroutine: cnstph_zero
!
! Description: Zeroes out all of the main data arrays/types used for constant pH
!
!*******************************************************************************

subroutine cnstph_zero(stateinf, protcnt, resname, resstate, chrgdat, statene)

  use mdin_ctrl_dat_mod, only : icnstph

  implicit none

! Passed Variables

  type(cnstph_info), intent(out)  :: stateinf(0:MAX_TITR_RES-1)
  integer, intent(out)            :: protcnt(0:MAX_TITR_STATES-1)
  integer, intent(out)            :: resstate(0:MAX_TITR_RES-1)
  character (len=40), intent(out) :: resname(0:MAX_TITR_RES)
  double precision, intent(out)   :: chrgdat(0:MAX_ATOM_CHRG-1)
  double precision, intent(out)   :: statene(0:MAX_TITR_STATES-1)

  ! cnstph_info type
  stateinf(:) = NULL_CPH_INFO

  ! Integers
  protcnt(:)  = 0
  resstate(:) = 0
  hindex(:)   = 0
  tpair(:)    = 0
  iselres(:)  = 0
  iselstat(:) = 0

  ! Chars
  resname(:)  = ' '

  ! Reals
  chrgdat(:)    = 0.d0
  statene(:)    = 0.d0

  return

end subroutine cnstph_zero

#ifdef MPI
!*******************************************************************************
!
! Subroutine: cnstph_bcast
!
! Description: Broadcasts all constant pH data to the whole MPI universe
!
!*******************************************************************************

subroutine cnstph_bcast(num_ints, num_reals)

  use mdin_ctrl_dat_mod, only : icnstph
  use prmtop_dat_mod, only : atm_qterm, natom
  use parallel_dat_mod  ! Provides MPI routines/constants, pmemd_comm and master

  implicit none

! Passed variables

  integer, intent(in out) :: num_ints
  integer, intent(in out) :: num_reals

  ! The slave nodes haven't been allocated yet...

  if (.not. master) call cnstph_allocate(num_ints, num_reals)

  call mpi_bcast(gbl_stateinf, MAX_TITR_RES * SIZE_CNSTPH_INFO, mpi_integer, &
                 0, pmemd_comm, err_code_mpi)

  call mpi_bcast(gbl_resstate, MAX_TITR_RES, mpi_integer, 0, pmemd_comm, &
                 err_code_mpi)

  call mpi_bcast(gbl_protcnt, MAX_TITR_STATES, mpi_integer, 0, pmemd_comm, &
                 err_code_mpi)

  call mpi_bcast(trescnt, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)

  call mpi_bcast(gbl_resname, 40 * (MAX_TITR_RES + 1), mpi_character, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(gbl_chrgdat, MAX_ATOM_CHRG, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(gbl_statene, MAX_TITR_STATES, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)
   
  call mpi_bcast(cph_igb, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)
  call mpi_bcast(cphfirst_sol, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)
  call mpi_bcast(cph_intdiel, 1, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)
  if (icnstph .eq. 2) &
    call mpi_bcast(mobile_atoms, natom, mpi_integer, 0, &
                   pmemd_comm, err_code_mpi)
  
end subroutine cnstph_bcast
#endif /* MPI */

!*******************************************************************************
!
! Subroutine: cnstph_setup
!
! Description: 
!
! Sets up constant pH calculation, launching the random number generator,
! setting up the hindex array so we know where all of the hydrogens are,
! and setting up the tpair array for the first time. This should be called by
! all threads.
!
!*******************************************************************************

subroutine cnstph_setup(crd)

  use mdin_ctrl_dat_mod, only : ig
  use parallel_dat_mod
  use prmtop_dat_mod,    only : natom, atm_qterm, gbl_one_scee
  use extra_pnts_nb14_mod

  implicit none

! Passed variables
  
  double precision, intent(in) :: crd(3, natom)

! Local variables

  integer  :: itres   ! residue counter
  integer  :: itstate ! state counter
  integer  :: itatm   ! atom counter
  integer  :: idx_ptr ! index pointer for building hindex pairlist
  integer  :: i       ! loop counter
  integer  :: randseed! Random number seed
  integer  :: fnd_atms(MAX_H_COUNT) ! which protons have been found so far

  logical  :: is_fnd  ! whether current proton is found yet or not

  ! Open up the cpout file if you're the master

  if (master) then
    ! For some reason, owrite is set to 'U' given the -O flag, but this just
    ! appends to the cpout file instead of overwriting it. Passing 'R' to amopen
    ! overwrites the file (this is the behavior in sander). Since this has been
    ! around since Bob first wrote the code, I don't want to remove it. Hence
    ! the kludge here.
    if (owrite .eq. 'U') then
      call amopen(cpout, cpout_name, 'R', 'F', 'W')
    else
      call amopen(cpout, cpout_name, owrite, 'F', 'W')
    end if
  end if

  ! Set initial charges based on initial states

  do itres = 0, trescnt - 1
    call cnstph_update_charge(atm_qterm, itres, gbl_resstate(itres))
  end do
#ifdef CUDA
  call gpu_refresh_charges(cit_nb14, gbl_one_scee, atm_qterm)
#endif

  ! Now copy this set of charges to the proposed charge set

  proposed_qterm(:) = atm_qterm(:)

  ! Set up the random number generator

  randseed = ig

#ifdef MPI
  call mpi_bcast(randseed, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)
#endif

  call amrset_gen(cnstph_randgen, randseed)

  ! Set up our first time through

  first_cpout = .true.
  first_cprestrt = .true.

  ! We don't do an exchange attempt on our first step

  on_cpstep = .false.

  ! The pointer section starts after the header section

  idx_ptr = trescnt + 1

  ! This is not a terribly efficient way of doing this, since it embeds a lot 
  ! of loops inside loops, but it's only done once at the beginning of a 
  ! simulation, and all of the loops are (reasonably) small, so this is not
  ! worth optimizing.

  do itres = 0, trescnt - 1

    hindex(itres) = idx_ptr

    fnd_atms(:) = 0

    do itatm = 0, gbl_stateinf(itres)%num_atoms - 1

      do itstate = 0, gbl_stateinf(itres)%num_states - 1

        ! Look for protons with a charge of 0 in at least 1 state
        if ( gbl_chrgdat(gbl_stateinf(itres)%first_charge + itstate * &
                         gbl_stateinf(itres)%num_atoms + itatm) .eq. 0.d0 ) then
          is_fnd = .false.
          do i = 1, MAX_H_COUNT
            if (fnd_atms(i) .eq. (gbl_stateinf(itres)%first_atom + itatm)) &
              is_fnd = .true.
          end do
          
          if (.not. is_fnd) then
            hindex(idx_ptr) = gbl_stateinf(itres)%first_atom + itatm
            idx_ptr = idx_ptr + 1

            do i = 1, MAX_H_COUNT
              if (fnd_atms(i) .eq. 0) then
                fnd_atms(i) = gbl_stateinf(itres)%first_atom + itatm
                exit
              end if
            end do

          end if

        end if ! gbl_chrgdat(gbl_stateinf...

      end do ! itstate

    end do ! itatm

  end do ! itres

  hindex(trescnt) = idx_ptr

  ! Initial pairlist build

  call cnstph_update_pairs(crd)

end subroutine cnstph_setup

!*******************************************************************************
!
! Subroutine: cnstph_update_pairs
!
! Description: Updates the neighboring titrating residue pairlist
!
!*******************************************************************************

subroutine cnstph_update_pairs(crd)

  use pmemd_lib_mod,  only : mexit
  use prmtop_dat_mod, only : natom

  implicit none

! Passed arguments

  double precision, intent(in) :: crd(3, natom) ! atomic coordinates

! Local variables

  integer :: idx_ptr     ! index into tpair array
  integer :: resi, resj  ! residue counters
  integer :: hi, hj      ! proton counters
  integer :: atmi, atmj  ! atom counters
  
  double precision :: x_ij, y_ij, z_ij  ! distances

  double precision, parameter ::  CUTOFF = 4.0d0 ! square of cutoff distance

  idx_ptr = trescnt + 1

  do resi = 0, trescnt - 1

    tpair(resi) = idx_ptr

    do resj = 0, trescnt - 1

      if (resi .eq. resj) cycle

      hloop: do hi = 0, hindex(resi + 1) - hindex(resi) - 1
        
        atmi = hindex(hi + hindex(resi))

        do hj = 0, hindex(resj + 1) - hindex(resj) - 1

          atmj = hindex(hj + hindex(resj))

          x_ij = crd(1, atmi) - crd(1, atmj)
          y_ij = crd(2, atmi) - crd(2, atmj)
          z_ij = crd(3, atmi) - crd(3, atmj)

          if (CUTOFF .gt. x_ij * x_ij + y_ij * y_ij + z_ij * z_ij) then

            if (idx_ptr < MAX_TITR_STATES * 4) then
              tpair(idx_ptr) = resj
              idx_ptr = idx_ptr + 1
              exit hloop
            else
              write(mdout, '(a)') 'Constant pH pair list overflow. Increase &
                                  &MAX_TITR_RES in constantph.fpp; recompile.'
              call mexit(mdout, 1)
            end if

          end if

        end do ! hj = 0

      end do hloop

    end do ! resj = 0

  end do ! resi = 0

  tpair(trescnt) = idx_ptr

end subroutine cnstph_update_pairs

!*******************************************************************************
!
! Subroutine: cnstph_begin_step
!
! Description: 
!
! Begins the constant pH step -- selects the titratable residue(s) that will be
! titrated as well as the state that will be attempted. There is a 25% chance 
! that we will attempt a multi-site move among neighboring residues. This is 
! done so we can at times observe a proton transfer when that would otherwise be
! a very rare event. The code only supports a 2-site move, but iselres/iselstat
! have been generalized to an arbitrary number of sites (as long as their
! dimensions are appropriately increased).
!
! We will first determine if we do a multisite move, then randomly select a
! residue to titrate. If we tried a multisite move, see if our selected site has
! any neighbors, then select from among those to titrate, also. Update the
! iselstat/iselres arrays to list which residues we're changing to which states.
! If we have no neighbors, don't do multisite.
!              
!*******************************************************************************

subroutine cnstph_begin_step

  use prmtop_dat_mod, only : natom

  implicit none

  integer          :: neighbor_cnt, i
  double precision :: random_value

  ! Do we do multisite?

  call amrand_gen(cnstph_randgen, random_value)
  
  if (random_value .lt. 2.5d-1) then
    iselres(0) = 2
  else
    iselres(0) = 1
  end if

  ! Select a random residue

  call amrand_gen(cnstph_randgen, random_value)
  iselres(1) = int( (random_value * 0.99999999d0) * trescnt)

  ! Select a random, unique state

  call amrand_gen(cnstph_randgen, random_value)
  iselstat(1) = int( (random_value * 0.99999999d0) * &
                     (gbl_stateinf(iselres(1))%num_states - 1) )
  
  if (iselstat(1) .ge. gbl_resstate(iselres(1))) &
    iselstat(1) = iselstat(1) + 1

  ! Update our proposed state charge array

  call cnstph_update_charge(proposed_qterm, iselres(1), iselstat(1))

  ! Determine our second titrating residue if we're doing more than 1

  if (iselres(0) .eq. 2) then
    neighbor_cnt = tpair(iselres(1) + 1) - tpair(iselres(1))

    if (neighbor_cnt .gt. 0) then

      ! New residue

      call amrand_gen(cnstph_randgen, random_value)
      iselres(2) = tpair( tpair(iselres(1)) + &
                          int(random_value * 0.99999999d0 * neighbor_cnt) )
      
      ! New state

      call amrand_gen(cnstph_randgen, random_value)
      iselstat(2) = int( (random_value * 0.99999999d0) * &
                         (gbl_stateinf(iselres(2))%num_states - 1) )
      if (iselstat(2) .ge. gbl_resstate(iselres(2))) &
        iselstat(2) = iselstat(2) + 1

      ! Update charge for this second residue

      call cnstph_update_charge(proposed_qterm, iselres(2), iselstat(2))

    else
      iselres(0) = 1
    end if

  end if

end subroutine cnstph_begin_step

!*******************************************************************************
!
! Subroutine: cnstph_end_step
!
! Description: Ends the constant pH step and accepts/rejects the MC move
!
!*******************************************************************************

subroutine cnstph_end_step(dvdl)

  use gbl_constants_mod, only : KB, LN_TO_LOG
  use mdin_ctrl_dat_mod, only : solvph, temp0
  use prmtop_dat_mod
  use charmm_mod
  use extra_pnts_nb14_mod

  implicit none

! Passed variables

  double precision, intent(in) :: dvdl

! Local variables

  double precision :: random_value
  double precision :: deltae
  double precision :: metrop

  integer          :: statebase
  integer          :: i

  deltae = 0.d0

  do i = 1, iselres(0)
    
    statebase = gbl_stateinf(iselres(i))%first_state

    ! deltae = E (proposed state) - E (current state)

    deltae = deltae + gbl_statene(iselstat(i) + statebase) - &
                      gbl_statene(gbl_resstate(iselres(i)) + statebase)

    ! Adjust for pH (delta protons * pH * LN_TO_LOG * KB * T)

    deltae = deltae - (gbl_protcnt(iselstat(i) + statebase) - &
                       gbl_protcnt(gbl_resstate(iselres(i)) + statebase)) * &
                       solvph * LN_TO_LOG * KB * temp0

  end do

  call amrand_gen(cnstph_randgen, random_value)

  metrop = exp( (deltae - dvdl) / (KB * temp0) )

! write(0, '(a,f16.5)') 'deltae  = ', deltae
! write(0, '(a,f16.5)') 'dvdl    = ', dvdl
! write(0, '(a,f16.5)') 'randval = ', random_value
! write(0, '(a,f16.5)') 'metrop  = ', metrop
! write(0, '(72("="))')

  if (random_value .le. metrop) then

    do i = 1, iselres(0)
      gbl_resstate(iselres(i)) = iselstat(i)
      call cnstph_update_charge(atm_qterm, iselres(i), iselstat(i))
    end do

    cph_success = .true.

#ifdef CUDA
    call gpu_refresh_charges(cit_nb14, gbl_one_scee, proposed_qterm)
#endif
    
  else
    
    ! Revert charges if we rejected the change

    do i = 1, iselres(0)
      call cnstph_update_charge(proposed_qterm, iselres(i), &
                                gbl_resstate(iselres(i)) )
    end do

  end if

end subroutine cnstph_end_step

!*******************************************************************************
!
! Subroutine: cnstph_update_charge
!
! Description: Update the given charge array with proposed states
!
!*******************************************************************************

subroutine cnstph_update_charge(chg_arry, selres, selstat)

  use prmtop_dat_mod, only : natom

  implicit none

! Passed variables

  double precision, intent(in out) :: chg_arry(natom) ! charge array to change
  integer, intent(in)              :: selres          ! selected residue
  integer, intent(in)              :: selstat         ! selected state

! Local variables

  integer :: itatm ! atom iterator

  do itatm = 0, gbl_stateinf(selres)%num_atoms - 1

    chg_arry(gbl_stateinf(selres)%first_atom + itatm) = &
      gbl_chrgdat(gbl_stateinf(selres)%first_charge + &
      selstat * gbl_stateinf(selres)%num_atoms + itatm)

  end do

end subroutine cnstph_update_charge

!*******************************************************************************
!
! Subroutine: cnstph_write_cpout
!
! Description: Writes constant pH titration information to cpout file
!
!*******************************************************************************

subroutine cnstph_write_cpout(nstep, nstlim, time, remd_method)

  use mdin_ctrl_dat_mod, only : ntwx, solvph, ntcnstph
  use parallel_dat_mod, only  : master

  implicit none

! Passed variables
  
  integer, intent(in) :: nstep
  integer, intent(in) :: nstlim
  integer, intent(in) :: remd_method

  double precision, intent(in) :: time

! Local variables

  logical :: full
  integer :: i

  ! Only the master writes cpout info

  if (.not. master) return

  if (ntwx .gt. 0) then
    full = (mod(nstep, ntwx) .eq. 0 .or. nstep .eq. nstlim .or. first_cpout)
  else
    full = (nstep .eq. nstlim .or. first_cpout)
  end if

  if (full) then

    write(cpout, '(a,f8.5)') 'Solvent pH: ', solvph
    write(cpout, '(a,i8)')   'Monte Carlo step size: ', ntcnstph
    write(cpout, '(a,i8)')   'Time step: ', nstep
    write(cpout, '(a,f10.3)') 'Time: ', time

    if (remd_method .eq. 4) then
      do i = 0, trescnt - 1
        write(cpout, '(a,i4,a,i2,a,f7.3)') &
            'Residue ', i, ' State: ', gbl_resstate(i), ' pH: ', solvph
      end do
    else
      do i = 0, trescnt - 1
        write(cpout, '(a,i4,a,i2)') 'Residue ', i, ' State: ', gbl_resstate(i)
      end do
    end if

  else if (remd_method .eq. 4) then

    do i = 1, iselres(0)
      write(cpout, '(a,i4,a,i2,a,f7.3)') 'Residue ', iselres(i), ' State: ', &
            gbl_resstate(iselres(i)), ' pH: ', solvph
    end do

  else

    do i = 1, iselres(0)
      write(cpout, '(a,i4,a,i2)') 'Residue ', iselres(i), ' State: ', &
            gbl_resstate(iselres(i))
    end do

  end if ! full

  write(cpout, '()')

  ! No longer first time around
  
  first_cpout = .false.

end subroutine cnstph_write_cpout

!*******************************************************************************
!
! Subroutine: cnstph_write_restart
!
! Description: writes the cprestrt file
!
!*******************************************************************************

subroutine cnstph_write_restart

  use gbl_constants_mod, only : ONE_AMBER_ELECTROSTATIC
  use parallel_dat_mod, only  : master
  use pmemd_lib_mod,     only : mexit

  implicit none

  ! Namelist variables -- same as gbl_ counterparts above

  double precision  :: chrgdat(0:MAX_ATOM_CHRG-1)
  double precision  :: statene(0:MAX_TITR_RES-1)
  type(cnstph_info) :: stateinf(0:MAX_TITR_RES-1)

  integer       :: protcnt(0:MAX_TITR_STATES-1)
  integer       :: resstate(0:MAX_TITR_RES-1)
  character(40) :: resname(0:MAX_TITR_RES)

  character(7)  :: stat

  integer       :: i

  namelist /cnstph/ stateinf, resstate, protcnt, chrgdat, statene, &
                    trescnt, resname, cphfirst_sol, cph_igb, cph_intdiel

  if (.not. master) return

  chrgdat(:) = gbl_chrgdat(:) * ONE_AMBER_ELECTROSTATIC
  statene(:) = gbl_statene(:)
  stateinf(:) = gbl_stateinf(:)
  protcnt(:) = gbl_protcnt(:)
  resstate(:) = gbl_resstate(:)
  resname(:) = gbl_resname(:)

  if (first_cprestrt) then
    if (owrite .eq. 'N') then
      stat = 'NEW'
    else if (owrite .eq. 'O') then
      stat = 'OLD'
    else if (owrite .eq. 'R') then
      stat = 'REPLACE'
    else if (owrite .eq. 'U') then
      stat = 'UNKNOWN'
    end if
    open(unit=cprestrt, file=cprestrt_name, status=stat, form='FORMATTED', &
         delim='APOSTROPHE', err=666)
    first_cprestrt = .false.
  else
    open(unit=cprestrt, file=cprestrt_name, status='OLD', form='FORMATTED', &
         delim='APOSTROPHE')
  end if

  write(cprestrt, nml=cnstph)

  close(cprestrt)

  return

666 write(mdout, '(a)') 'Error opening ', trim(cprestrt_name)
    call mexit(mdout, 1)

end subroutine cnstph_write_restart

!*******************************************************************************
!
! Function: total_protonation
!
! Description: Finds the total number of 'active' protons
!
!*******************************************************************************

integer function total_protonation()

  implicit none

  integer :: i

  total_protonation = 0

  do i = 0, trescnt - 1
    total_protonation = total_protonation + &
             gbl_protcnt(gbl_stateinf(i)%first_state + gbl_resstate(i))
  end do

  return

end function total_protonation

!*******************************************************************************
!
! Subroutine: explicit_cnstph_begin_step
!
! Description: The beginning of a constant pH step in explicit solvent
!
!*******************************************************************************

subroutine explicit_cnstph_begin_step(dcharge, done_res, num_done)

  implicit none

! Passed variables

  double precision, intent(in out) :: dcharge(1:*)
  integer, intent (in out)         :: done_res(1:MAX_TITR_RES)
  integer, intent(in)              :: num_done

! Local variables

  integer :: i
  integer :: j
  integer :: neighbor_cnt
  integer :: holder_res(1:MAX_TITR_RES)

  double precision :: randval

! This subroutine will randomly select what to titrate, as long as it has not
! yet been titrated. That is the only difference from cnstph_begin_step

  call amrand_gen(cnstph_randgen, randval)
  if(randval .lt. 0.25d0) then
    iselres(0) = 2
  else
    iselres(0) = 1
  end if

  ! Select a random residue
  call amrand_gen(cnstph_randgen, randval)
  iselres(1) = int((randval * 0.9999999d0) * (trescnt - num_done))

  if (num_done .eq. 0) then
    done_res(1) = iselres(1)
  else

    ! Now move the selected residue off of those that have been selected already
    ! THIS ARRAY SHOULD ALREADY BE SORTED (see below)

    do i = 1, num_done
      if (iselres(1) .ge. done_res(i)) iselres(1) = iselres(1) + 1
    end do

    ! If the residue we chose was larger than our last number, add it to the end
    if (iselres(1) .gt. done_res(num_done)) &
      done_res(num_done + 1) = iselres(1)

  end if ! (num_done .eq. 0)

  do i = 1, num_done
    if (iselres(1) .lt. done_res(i)) then
      holder_res(i:num_done) = done_res(i:num_done)
      done_res(i) = iselres(1)
      done_res(i+1:num_done+1) = holder_res(i:num_done)
      exit
    end if
  end do

  ! Select a random, unique state

  call amrand_gen(cnstph_randgen, randval)
  iselstat(1) = int( (randval * 0.99999999d0) * &
                     (gbl_stateinf(iselres(1))%num_states - 1) )
  
  if (iselstat(1) .ge. gbl_resstate(iselres(1))) &
    iselstat(1) = iselstat(1) + 1

  ! Update our proposed state charge array

  call cnstph_update_charge(proposed_qterm, iselres(1), iselstat(1))

  ! Determine our second titrating residue if we're doing more than 1

  if (iselres(0) .eq. 2) then
    neighbor_cnt = tpair(iselres(1) + 1) - tpair(iselres(1))

    if (neighbor_cnt .gt. 0) then

      ! New residue

      call amrand_gen(cnstph_randgen, randval)
      iselres(2) = tpair( tpair(iselres(1)) + &
                          int(randval * 0.99999999d0 * neighbor_cnt) )
      
      ! New state

      call amrand_gen(cnstph_randgen, randval)
      iselstat(2) = int( (randval * 0.99999999d0) * &
                         (gbl_stateinf(iselres(2))%num_states - 1) )
      if (iselstat(2) .ge. gbl_resstate(iselres(2))) &
        iselstat(2) = iselstat(2) + 1

      ! Update charge for this second residue

      call cnstph_update_charge(proposed_qterm, iselres(2), iselstat(2))

    else
      iselres(0) = 1
    end if

  end if

  return

end subroutine explicit_cnstph_begin_step

!*******************************************************************************
!
! Subroutine: cnstph_explicitmd
!
! Description: Do MC cycles in explicit solvent
!
!*******************************************************************************

#ifdef MPI
subroutine cnstph_explicitmd(atm_cnt, crd, mass, frc, vel, last_vel, &
                             nstep, nstlim, time, my_atm_lst)
#else
subroutine cnstph_explicitmd(atm_cnt, crd, mass, frc, vel, last_vel, &
                             nstep, nstlim, time)
#endif

  use gb_force_mod, only       : gb_cph_ene
  use mdin_ctrl_dat_mod, only  : ntp, ntb, ips, igb, nmropt, nscm, ntrelax, iamd
  use nmr_calls_mod, only      : nmrdcp
#ifdef MPI
  use parallel_mod, only       : mpi_allgathervec
#endif
  use parallel_dat_mod, only   : master
  use prmtop_dat_mod, only     : natom, atm_qterm, atm_gb_radii
  use remd_mod, only           : remd_method
  use relaxmd_mod, only        : relaxmd

  implicit none

  ! Passed variables

#ifdef MPI
  integer, intent(in)     :: my_atm_lst(*)
#endif
  integer, intent(in)     :: nstep
  integer, intent(in)     :: nstlim
  integer, intent(in out) :: atm_cnt
  double precision, intent(in)     :: time
  double precision, intent(in out) :: crd(3, atm_cnt)
  double precision, intent(in)     :: mass(atm_cnt)
  double precision, intent(in out) :: frc(3, atm_cnt)
  double precision, intent(in out) :: vel(3, atm_cnt)
  double precision, intent(in out) :: last_vel(3, atm_cnt)

  ! Local variables

  integer :: holder_ntb, holder_natom, holder_ntp, holder_ips, j, i, k, &
             natom3, selres_holder(1:MAX_TITR_RES), holder_nscm, holder_iamd

#ifdef CUDA
  integer :: fat
#endif

  double precision, dimension(3, atm_cnt) :: vtemp

  double precision :: dvdl_current
  double precision :: dvdl_proposed
  double precision :: dvdl

  logical          :: any_success

  ! This subroutine is the main driver for running constant pH MD in explicit
  ! solvent. The first thing it does is call being_step to set up the MC.
  ! Then it removes the periodic boundary conditions (PBC), selects a GB model,
  ! then calls force with nstep = 0 which will force oncpstep to be .true., so
  ! we'll get dvdl back. This is then passed to cnstphendstep to evaluate the
  ! transition. If it fails, then we restore the PBC, turn off the GB model,
  ! restore the CUToff, and return to the calling routine. If it's successful,
  ! we still return the above variables, but then we also turn on belly and fix
  ! the protein while we call runmd to relax the solvent for ntrelax steps.
  ! Some important things we do are:
  !  o  turn off trajectory/energy printing
  !  o  zero-out the v and vold arrays for the solute
  !  o  turn on belly and assign belly arrays

! Collect all of the coordinates

#ifdef MPI
  ! Collect all coordinates
  call mpi_allgathervec(atm_cnt, crd)
#endif

  natom3          = atm_cnt * 3
  holder_natom    = atm_cnt
  holder_ntb      = ntb
  holder_ntp      = ntp
  holder_ips      = ips

  natom   = cphfirst_sol - 1
  atm_cnt = cphfirst_sol - 1
  ntb     = 0
  igb     = cph_igb

  dvdl_current = 0.d0
  dvdl_proposed = 0.d0

  ! Zero-out the holder for the selected residues

  selres_holder(1:trescnt) = 0

  ! Initiate trial moves for protonation states

#ifdef CUDA
  call gpu_download_crd(crd)
#endif
  
  any_success = .false.

  ! Get the initial energy
  call gb_cph_ene(atm_cnt, crd, atm_qterm, dvdl_current, .false.)

  do i = 1, trescnt

    cph_success = .false.

    call explicit_cnstph_begin_step(proposed_qterm, selres_holder, i-1)

    ! Call gb force to get energies

    call gb_cph_ene(atm_cnt, crd, proposed_qterm, dvdl_proposed, .true.)

    dvdl = dvdl_proposed - dvdl_current
!   write(0,'(a,2i5,f20.10)') 'Selected residue, state, dvdl is ', &
!            selres_holder(i), iselstat(1), dvdl

    ! End the step

    call cnstph_end_step(dvdl)

    ! If we succeeded, update our 'current' state energy
    if (cph_success) &
      dvdl_current = dvdl_proposed

    any_success = any_success .or. cph_success

    ! Decrement the nmropt counter

    if (nmropt .ne. 0) call nmrdcp

  end do

! if (master) write(666, '(80("="))')

  ! Restore PBC variables

  natom   = holder_natom
  atm_cnt = holder_natom
  ntb     = holder_ntb
  igb     = 0
  ntp     = holder_ntp
  ips     = holder_ips

  ! Use iselres to stipulate that every residue is printed out

  iselres(0) = trescnt
  do i = 1, trescnt
    iselres(i) = i - 1
  end do

  if (master) call cnstph_write_cpout(nstep, nstlim, time, remd_method)

  ! If we did not succeed, bail out

  if (.not. any_success) return

  holder_nscm     = nscm
  holder_iamd     = iamd
  nscm    = 0
  iamd    = 0

  ! If we succeeded, adjust the variables for relaxation dynamics. Don't print
  ! during relaxation steps. Also zero the velocity array, and don't remove COM
  ! motion

#ifdef CUDA
  fat = cphfirst_sol - 1
  call gpu_download_vel(vel)
  call gpu_set_first_update_atom(fat)
#endif

#ifdef MPI
  call mpi_allgathervec(atm_cnt, vel)
#endif

  do i = 1, cphfirst_sol - 1
    vtemp(:, i) = vel(:, i)
    vel(:, i) = 0.d0
  end do

#ifdef CUDA
  call gpu_upload_vel(vel)
#endif
  ! Call routine to do relaxation dynamics

#ifdef MPI
  call relaxmd(atm_cnt, crd, mass, frc, vel, last_vel, my_atm_lst, &
               mobile_atoms, ntrelax)
  call mpi_allgathervec(atm_cnt, crd)
  call mpi_allgathervec(atm_cnt, vel)
#else
  call relaxmd(atm_cnt, crd, mass, frc, vel, last_vel, &
               mobile_atoms, ntrelax)
#endif

  ! Restore the original velocities

#ifdef CUDA
  call gpu_download_vel(vel)
#endif

  do i = 1, cphfirst_sol - 1
    vel(:,i) = vtemp(:,i)
  end do

#ifdef CUDA
  fat = 0
  call gpu_upload_vel(vel)
  call gpu_set_first_update_atom(fat)
#endif

  nscm = holder_nscm
  iamd = holder_iamd

  return 

end subroutine cnstph_explicitmd

end module constantph_mod
