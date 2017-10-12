#include "copyright.i"

!*******************************************************************************
!
! Module: master_setup_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module master_setup_mod

use file_io_dat_mod

  implicit none

! Hide internal routines:

  private       open_output_files

contains

!*******************************************************************************
!
! Subroutine:  master_setup
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine master_setup(num_ints, num_reals, new_stack_limit, terminal_flag)

  use axis_optimize_mod
  use charmm_mod, only : charmm_active
  use cmap_mod, only : generate_cmap_derivatives
  use cit_mod
  use constantph_mod, only : cnstph_read, cph_igb
  use constraints_mod
  use dynamics_mod
  use dynamics_dat_mod
  use emap_mod, only: emap_options
  use extra_pnts_nb14_mod
  use ene_frc_splines_mod
  use pme_direct_mod
  use pme_force_mod
  use file_io_mod
  use findmask_mod, only : atommask
  use gbl_constants_mod
  use get_cmdline_mod
  use img_mod
  use inpcrd_dat_mod
  use mdin_ctrl_dat_mod
  use mdin_debugf_dat_mod, only : init_mdin_debugf_dat, validate_mdin_debugf_dat
  use mdin_ewald_dat_mod
  !Added by Lijiang Yang for Integrated Tempering Sampling
  use mdin_its_dat_mod
  use mdin_subunit_dat_mod
  !-------------------------------------------------------
  use nb_exclusions_mod
  use nb_pairlist_mod
  use nmr_calls_mod
  use pbc_mod
  use pmemd_lib_mod
  use prmtop_dat_mod
  use parallel_dat_mod
#ifdef MPI
  use parallel_mod
#endif /* MPI */
  use remd_mod
  use ti_mod
  use shake_mod
#ifdef EMIL
  use mdin_emil_dat_mod 
#endif

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out) :: num_ints, num_reals
  integer, intent(in)     :: new_stack_limit
  logical, intent(out)    :: terminal_flag

! Local variables:

  double precision      :: box_alpha, box_beta, box_gamma
  double precision      :: box(3)
  integer               :: i, j
  integer               :: ifind
  integer               :: inerr
  integer               :: inpcrd_natom
  character(80)         :: inpcrd_title
  character(80)         :: mdin_title
  integer               :: ord1, ord2, ord3
  integer               :: itmp(3)
  double precision      :: rtmp(3)
  character(8)          :: date
  character(10)         :: time
  character(512)        :: char_tmp_512

  inerr = 0

  terminal_flag = .false.
  call get_cmdline(terminal_flag)   ! Get the file names.
  
! terminal_flag is set to .true. inside get_cmdline if --version or --help is
! given. If this happened in the groupfile, we quit in error already. We just
! return here if terminal_flag is .true. and let the calling routine exit (to
! avoid MPI hangups)
  if (terminal_flag) return

! Read the control data and open different files. 

  call amopen(mdin, mdin_name, 'O', 'F', 'R')
  call amopen(mdout, mdout_name, owrite, 'F', 'W')

  write(mdout, 1000)
  write(mdout, '(a, /)') '| PMEMD implementation of SANDER, Release 14'

  call date_and_time(DATE=date, TIME=time)

  write(mdout,'(12(a),/)') '| Run on ', date(5:6), '/', date(7:8), '/',  &
        date(1:4), ' at ', time(1:2), ':', time(3:4), ':', time(5:6)
! Write the path of the current executable and working directory
  call get_command_argument(0, char_tmp_512)
  write(mdout,'(a,a)') '|   Executable path: ', trim(char_tmp_512)
  call getcwd(char_tmp_512)
  write(mdout,'(a,a)') '| Working directory: ', trim(char_tmp_512)
! Write the hostname if we can get it from environment variable
! Note: get_environment_variable is part of the F2003 standard but seems
!       to be supported by GNU, Intel, IBM and Portland (2010+) compilers
  call get_environment_variable("HOSTNAME", char_tmp_512, inerr)
  if (inerr .eq. 0) then
    write(mdout,'(a,a,/)')  '|          Hostname: Unknown'
  else
    write(mdout,'(a,a,/)')  '|          Hostname: ', trim(char_tmp_512)
  end if

! Print warning if the stack could not be unlimited
  if (master .and. new_stack_limit .gt. 0) then
    write(mdout, '(a,a,i10,a)') warn_hdr, &
      'Stack usage limited by a hard resource limit of ', new_stack_limit, ' bytes!'
    write(mdout, '(a,a)') extra_line_hdr, &
      'If segment violations occur, get your sysadmin to increase the limit.'

    ! Flush the output in case we segfault - note calls to flush() are not 
    ! reliable and can be machine dependent.
    close(mdout)
    open(unit=mdout, file=mdout_name, status='OLD', position='APPEND')

  end if

  if (owrite .eq. 'U') write(mdout, '(2x,a,/)') '[-O]verwriting output'

! Echo the file assignments to the user:

  write(mdout, 1010) 'MDIN',   mdin_name(1:70),   'MDOUT',  mdout_name(1:70),  &
                     'INPCRD', inpcrd_name(1:70), 'PARM',   prmtop_name(1:70), &
                     'RESTRT', restrt_name(1:70), 'REFC',   refc_name(1:70),   &
                     'MDVEL',  mdvel_name(1:70),  'MDEN',   mden_name(1:70),   &
#ifdef MPI
                     'MDCRD',  mdcrd_name(1:70),  'MDINFO', mdinfo_name(1:70), &
                     'LOGFILE',  logfile_name(1:70), 'MDFRC', mdfrc_name(1:70)
#else
                     'MDCRD',  mdcrd_name(1:70),  'MDINFO', mdinfo_name(1:70), &
                     'MDFRC',  mdfrc_name(1:70)
#endif

  call echoin(5, 6)     ! Echo the input file to the user.

  write(mdout, '(/)')

! Read data characterizing the md-run:

! Read the title line in the mdin file:

  read(mdin, '(a80)') mdin_title   ! BUGBUG - No longer serves a purpose...

! Read the cntrl namelist in the mdin file:

  call init_mdin_ctrl_dat(remd_method)
  call validate_mdin_ctrl_dat(remd_method)

#ifdef EMIL
  ! Read the EMIL namelist (if wanted and one is present)
  if( emil_do_calc .gt. 0 ) then
     call init_emil_dat()
  end if
#endif

#ifdef CUDA
  !RCW: Call to set device is mostly redundant now due to'
  !     removal of -gpu command line argument. But leave for now.
  call gpu_set_device(-1)
  call gpu_init()

  !Write info about CUDA GPU(s)
#ifdef MPI
  call gpu_write_cuda_info(mdout, mytaskid, numtasks, using_pme_potential,iamd)

#else
  call gpu_write_cuda_info(mdout, 1, 1, using_pme_potential,iamd)

#endif

  ! Check for GeForce GTX TITAN or 780 and NVIDIA driver version
  call gpu_check_titan()

#endif


  ! Read the input coordinates or restart file (inpcrd):

  call init_inpcrd_dat(num_ints, num_reals, inpcrd_natom, &
                       box_alpha, box_beta, box_gamma, &
                       box, t, inpcrd_title)



  if (using_pme_potential) then

    ! Note that the box parameters may be overridden by ewald input:

    call init_mdin_ewald_dat(box_alpha, box_beta, box_gamma, box)

    ! After the above call, box() and nfft1..3 may have been flipped to internal
    ! optimized coordinates.  Beware!  Flipping is only done for orthogonal unit
    ! cells, so there is no need to worry with the angles.

    call validate_mdin_ewald_dat(box_alpha, box_beta, box_gamma, box)

  else

    ! Turn off axis flipping, which is only relevant under pme.

    call setup_axis_opt(1.d0, 1.d0, 1.d0)
  end if

! Added by Lijiang Yang: Read in options for Integrated Tempering Sampling 
  if (doits .eq. 1) then
    call init_mdin_its_dat
    call validate_mdin_its_dat
  endif
! Added by Lijiang Yang: Read in options for subunit energy
  if (isub .eq. 1) then
    call init_mdin_subunit_dat(inpcrd_natom)
    call validate_mdin_subunit_dat
  endif

! Print conditional compilation flag information:

  call printdefines()

  if (ntb .ne. 0) then

    ! Initialize stuff associated with periodic boundary conditions:

    call init_pbc(box(1), box(2), box(3), box_alpha, box_beta, box_gamma, &
                  vdw_cutoff + skinnb)
 
    ! Nonisotropic scaling / nonorthorhombic unit cell check:

    if (ntp .gt. 1) then
      if (abs(box_alpha - 90.d0) .gt. 1.d-5 .or. &
          abs(box_beta - 90.d0) .gt. 1.d-5 .or. &
          abs(box_gamma - 90.d0) .gt. 1.d-5) then
        write(mdout, '(a,a,a)') error_hdr, &
                                'Nonisotropic scaling on nonorthorhombic ', &
                                'unit cells is not permitted.'
        write(mdout,'(a,a,a)') extra_line_hdr, &
                               'Please use ntp=1 if unit cell angles are ', &
                               'not 90 degrees.'
        call mexit(6,1)
      end if
    end if

  end if

  if (using_pme_potential) then

    ! Set up coordinate index table dimensions; coordinate flipping is done for
    ! the output values if necessary (pbc_box is already flipped).

    call set_cit_tbl_dims(pbc_box, vdw_cutoff + skinnb, cut_factor)

  end if

! We only support a subset of the debugf namelist stuff, so if there is a debugf
! namelist, issue a warning:

  rewind(mdin)                      ! Insurance against maintenance mods.

  call nmlsrc('debugf', mdin, ifind)

  if (ifind .ne. 0) then
    write(mdout, '(a,/)') warn_hdr, &
      'debugf namelist found in mdin. PMEMD only supports a very small subset &
      &of the debugf options. Unsupported options will be ignored.'

  end if

! Read the parameters and topology file.  This routine reads all parameters
! used by amber pme ff's, or amber generalized Born ff's.

  call init_prmtop_dat(num_ints, num_reals, inpcrd_natom)

! Read the cpin file if we're doing constant pH MD. We need natom, so this must
! appear after we initialize the prmtop data

  if (icnstph .gt. 0) call cnstph_read(num_ints, num_reals)

  if (charmm_active) call generate_cmap_derivatives

!DEBUG
! do i=1,cmap_type_count
!   write(mdout, *)"gbl_cmap_res(i)", gbl_cmap_res(i)
!   write(mdout, *)"gbl_cmap_grid(i)", gbl_cmap_grid(i,1:24,1:24)
!   write(mdout, *)"gbl_cmap_dPhi(i)", gbl_cmap_dPhi(i,1:24,1:24)
!   write(mdout, *)"gbl_cmap_dPsi(i)", gbl_cmap_dPsi(i,1:24,1:24)
!   write(mdout, *)"gbl_cmap_dPhi_dPsi(i)", gbl_cmap_dPhi_dPsi(i,1:24,1:24)
! enddo


! If the user has requested NMR restraints, do a cursory read of the
! restraints file(s) now to determine the amount of memory necessary
! for these restraints, and allocate the memory in the master.

  if (nmropt .ne. 0) call init_nmr_dat(num_ints, num_reals)

! Process EMAP restraints is requested

  if (iemap > 0) call emap_options(mdin)

! Now do axis flip optimization. The way setup works on this is that if it
! was not selected then the axes_flip routine actually does nothing.  At
! present, the decision to do axis flipping occurs in init_mdin_ewald_dat()
! because in this routine we first have access to the final unit cell lengths
! and angles and we need to do flipping if we are going to because nfft1,2,3
! must be set.

  do i = 1, natom
    call axes_flip(atm_crd(1,i), atm_crd(2,i), atm_crd(3,i))
  end do

  ! Check input crds for NaNs.
  ! NaNs have no value therefore can be identified by checking the variable
  ! against itself.  
  do i = 1, natom
     do j = 1, 3
       if (atm_crd(j, i) /= atm_crd(j, i)) then
          write(mdout,'(a,a,a)') error_hdr, "NaN(s) found in input coordinates." 
          write(mdout,'(11x,a)') "This likely means that something went wrong in the previous simulation."
          call mexit(6,1) 
       end if
    enddo       
  end do 

  if (ntx .ne. 1 .and. ntx .ne. 2) then

    do i = 1, natom
      call axes_flip(atm_vel(1,i), atm_vel(2,i), atm_vel(3,i))
    end do

  ! Check input velocities for NaNs.
  ! NaNs have no value therefore can be identified by checking the variable
  ! against itself.  
    do i = 1, natom
       do j = 1, 3
         if (atm_vel(j, i) /= atm_vel(j, i)) then
            write(mdout,*) error_hdr, "NaN(s) found in input velocities. "
            write(mdout,'(11x,a)') "This likely means that something went wrong in the previous simulation."
            call mexit(6,1) 
         end if
       enddo       
    end do 

  end if

  if (using_pme_potential) then

    ord1 = axis_flipback_ords(1)
    ord2 = axis_flipback_ords(2)
    ord3 = axis_flipback_ords(3)
    itmp(1) = cit_tbl_x_dim
    itmp(2) = cit_tbl_y_dim
    itmp(3) = cit_tbl_z_dim

    write(mdout, '(a, 3i5)')'| Coordinate Index Table dimensions: ', &
                            itmp(ord1), itmp(ord2), itmp(ord3)

    rtmp(1) = pbc_box(1)/cit_tbl_x_dim
    rtmp(2) = pbc_box(2)/cit_tbl_y_dim
    rtmp(3) = pbc_box(3)/cit_tbl_z_dim

    write(mdout,'(a, 3f10.4, /)')'| Direct force subcell size = ', &
                                 rtmp(ord1), rtmp(ord2), rtmp(ord3)

  end if

! Init constraints (and belly) data:

  call init_constraints_dat(natom, ibelly, ntr, num_ints, num_reals)

  if (using_pme_potential) then

    ! Set up image dynamic memory:

    call alloc_img_mem(natom, num_ints, num_reals)

    ! Set up pairlist memory:

    call alloc_nb_pairlist_mem(natom, vdw_cutoff + skinnb, &
                               num_ints, num_reals)

    call alloc_nb_exclusions_mem(natom, next, num_ints, num_reals)

    ! Set up ewald variables and memory:

    call alloc_pme_force_mem(ntypes, num_ints, num_reals)

    call init_pme_direct_dat(num_ints, num_reals)

    call init_ene_frc_splines_dat(num_ints, num_reals)

  end if

! Code added to detect the existence of any 10-12 terms that must be
! examined.  If none are found, it speeds up the nonbonded calculations.

#ifdef HAS_10_12
#else
  do i = 1, nphb
    if (gbl_asol(i) .ne. 0.d0 .or. gbl_bsol(i) .ne. 0.d0) then
      write(mdout, '(a,a,a)') error_hdr, &
                              'Found a non-zero 10-12 coefficient, but ',&
                              ' source was not compiled with -DHAS_10_12.'
      write(mdout,'(a,a,a)') extra_line_hdr, &
                             'If you are using a pre-1994 force field, you',&
                             ' will need to re-compile with this flag.'
      call mexit(6,1)
    end if
  end do
#endif

! (ifbox comes from prmtop & is for indicating presence & type of box).

  if (ifbox .eq. 1) write(mdout, '(5x,''BOX TYPE: RECTILINEAR'',/)')
  if (ifbox .eq. 2) write(mdout, '(5x,''BOX TYPE: TRUNCATED OCTAHEDRON'',/)')
  if (ifbox .eq. 3) write(mdout, '(5x,''BOX TYPE: GENERAL'',/)')

! Print data characterizing the md-run. 

! If the axis flipping optimization is in effect, we want to restore the
! box lengths to original values for printout...

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  ! Print control data header:

  write(mdout, 1040)

  ! Strangely enough the prmtop title occurs next:

  write(mdout, '(a80)') prmtop_ititl

  ! Then the &ctrl data:

  call print_mdin_ctrl_dat(remd_method, cph_igb)

  ! Then the &ewald data:

  if (using_pme_potential) &
    call print_mdin_ewald_dat(box_alpha, box_beta, box_gamma, box, es_cutoff)

  ! Then the &subunit data: Added by Lijiang Yang: Read in options for Subunit Energy
  if (isub .eq. 1) then
    call print_mdin_subunit_dat
  endif
  ! Then the &its data: Added by Lijiang Yang: Read in options for Integrated Tempering Sampling 
  if (doits .eq. 1) then
    call print_mdin_its_dat
  endif

! Check if ifbox variable from prmtop file matches actual angles.  This must
! occur immediately after print_mdin_ewald_dat for consistency with sander11
! mdout output.

  if (ifbox .eq. 1) then
    if (abs(box_alpha - 90.0d0 ) .gt. 1.d-5 .or. &
        abs(box_beta - 90.0d0) .gt. 1.d-5 .or. &
        abs(box_gamma - 90.0d0) .gt. 1.d-5) then
      ifbox = 3
      write(mdout,'(a)') '     Setting ifbox to 3 for non-orthogonal unit cell'
    end if
  end if

  if (ifbox .eq. 2) then
    if (abs(box_alpha - 109.4712190d0) .gt. 1.d-5 .or. &
        abs(box_beta - 109.4712190d0) .gt. 1.d-5 .or. &
        abs(box_gamma - 109.4712190d0) .gt. 1.d-5) then
      write(mdout,'(/2x,a)') &
        'Error: ifbox=2 in prmtop but angles are not correct'
      inerr = 1
    end if
  end if

! Consistency checking:

  if (ntb .ne. 0 .and. ntp .ne. 0 .and. ifbox .eq. 0) then
    write(mdout, '(a,a)') error_hdr, &
      'the combination ntb != 0, ntp != 0, ifbox == 0 is not supported!'
    inerr = 1
  end if

  if (using_pme_potential) then

    if (vdw_cutoff .ge. box(1) * 0.5d0 .or. &
        vdw_cutoff .ge. box(2) * 0.5d0 .or. &
        vdw_cutoff .ge. box(3) * 0.5d0) then
      write(mdout, '(a,a)') error_hdr, &
                            'max cut must be < half smallest box dimension!'
      write(mdout, '(a,a)') extra_line_hdr, 'max cut=', vdw_cutoff
      write(mdout, '(a,a)') extra_line_hdr, 'box(1)=', box(ord1)
      write(mdout, '(a,a)') extra_line_hdr, 'box(2)=', box(ord2)
      write(mdout, '(a,a)') extra_line_hdr, 'box(3)=', box(ord3)
      inerr = 1
    end if

  end if

! Warnings:

  if (using_pme_potential .and. ibelly .gt. 0) then
     write(mdout, '(a,/,a,/)') 'Warning: Although EWALD will work with belly', &
           '(for equilibration), it is not strictly correct!'
  end if

  if (inerr .eq. 1) then
    write(mdout, '(a,/)') ' Input errors occurred. Terminating execution.' 
    call mexit(6, 1)
  end if

! Load the constrained (or belly) atoms. If their respective masks are
! not set, they are read as groups:

  natc = 0
  belly_atm_cnt = 0

  if (ntr .gt. 0) then
    write(mdout, '(/4x,a,/)') 'LOADING THE CONSTRAINED ATOMS AS GROUPS'
    call read_restraints(natom, ntrx, atm_xc)

    ! Axis flip optimization...

    do i = 1, natom
      call axes_flip(atm_xc(1,i), atm_xc(2,i), atm_xc(3,i))
    end do

    ! Load constrained atoms as mask if it's provided, or GROUP if it's not
    if (len_trim(restraintmask) .le. 0) then
      call rgroup(natom, natc, nres, gbl_res_atms, gbl_labres, &
                  atm_igraph, atm_isymbl, atm_itree, atm_jrc, atm_weight, &
                  .true., .false., 5)
    else
      call atommask(natom, nres, 0, atm_igraph, atm_isymbl, gbl_res_atms, &
                    gbl_labres, atm_crd, restraintmask, atm_jrc)

      ! Gather constrained atoms together as is done in file_io_mod's rgroup
      natc = 0
      do i = 1, natom
        if (atm_jrc(i) .le. 0) cycle
        natc = natc + 1
        atm_jrc(natc) = i
        atm_weight(natc) = restraint_wt !!!BUGBUGJWK - Don't need to reset the whole array each loop ... 
      end do

      write(6,'(a,a,a,i5,a)') '     Mask ', &
            restraintmask(1:len_trim(restraintmask)), &
            ' matches ', natc, ' atoms'
    end if
  end if

  if (ibelly .gt. 0) then
    write(mdout, '(/4x,a,/)') 'LOADING THE BELLY ATOMS AS GROUPS'

    ! load belly atoms as mask if it's provided, or GROUP if it's not
    if (len_trim(bellymask) .le. 0) then
      call rgroup(natom, belly_atm_cnt, nres, gbl_res_atms, gbl_labres, &
           atm_igraph, atm_isymbl, atm_itree, atm_igroup, atm_weight, &
           .false., .true., 5)
    else
      call atommask(natom, nres, 0, atm_igraph, atm_isymbl, gbl_res_atms, &
                    gbl_labres, atm_crd, bellymask, atm_igroup)
      
      belly_atm_cnt = 0
      do i = 1, natom
        if (atm_igroup(i) .gt. 0) &
          belly_atm_cnt = belly_atm_cnt + 1
      end do

      write(6,'(a,a,a,i5,a)') '     Mask ', bellymask(1:len_trim(bellymask)), &
            ' matches ', belly_atm_cnt, ' atoms'
    end if
  end if

  ti_mode = 0 ! Always initialize
  if(icfe .ne. 0) then
    call ti_alloc_mem(natom, ntypes, num_ints, num_reals)
    call ti_init_mask(natom, nres, atm_igraph, atm_isymbl, &
      gbl_res_atms, gbl_labres, atm_crd, timask1, timask2, &
      scmask1, scmask2)   
    call ti_crgmask(natom, nres, atm_igraph, atm_isymbl, gbl_res_atms, &
      gbl_labres, atm_crd, crgmask, atm_qterm)
    call ti_init_dat(natom, atm_crd, icfe, ifsc, atm_qterm)
  else
    ! Allow for crgmask even w/o TI 
    call ti_crgmask(natom, nres, atm_igraph, atm_isymbl, gbl_res_atms, &
      gbl_labres, atm_crd, crgmask, atm_qterm)
  end if

  call shake_noshakemask(natom, nres, atm_igraph, atm_isymbl, gbl_res_atms, &
                         gbl_labres, atm_crd, noshakemask, num_ints, num_reals)
! All the bond, angle, and dihedral parameters may be changed here as the
! bond, angle, and dihedral arrays are repacked! Note in particular that
! diheda_idx may also be changed.  We also count atoms in the "belly" here,
! which is probably redundant (also done in rgroup() above).

  if (ibelly .gt. 0) then

    call remove_nonbelly_bnd_ang_dihed
    call count_belly_atoms(natom, belly_atm_cnt, atm_igroup)

    if (.not. using_pme_potential) then

      ! The only allowable belly here has just the first belly_atm_cnt atoms
      ! in the moving part.  Confirm this.

      do i = belly_atm_cnt + 1, natom
        if (atm_igroup(i) .ne. 0) then
          write(mdout, *)'When ibelly != 0 and igb != 0, the moving part must'
          write(mdout, *)'  be at the start of the molecule, which seems to'
          write(mdout, *)'  not be the case!'
          call mexit(6, 1)
        end if
      end do

    end if

  end if

  ! Make the bond arrays sequential for shake and force routines:

  do i = 1, nbona
    gbl_bond(nbonh + i) = gbl_bond(bonda_idx + i - 1)
  end do

  bonda_idx = nbonh + 1

  ! Make the angle arrays sequential:

  do i = 1, ntheta
    gbl_angle(ntheth + i) = gbl_angle(anglea_idx + i - 1)
  end do

  anglea_idx = ntheth + 1

  ! Make the dihedrals sequential:

  do i = 1, nphia
    gbl_dihed(nphih + i) = gbl_dihed(diheda_idx + i - 1)
  end do

  diheda_idx = nphih + 1

  if (using_pme_potential) then
    if (numextra .eq. 0) then
      call init_nb14_only(num_ints, num_reals)
    else
      call init_extra_pnts_nb14(num_ints, num_reals)
      ! Make sure input coordinates have correctly placed extra points:
      if (frameon .ne. 0 .and. gbl_frame_cnt .gt. 0) &
        call all_local_to_global(atm_crd, ep_frames, ep_lcl_crd, gbl_frame_cnt)
    end if
  else if (using_gb_potential) then
    call init_nb14_only(num_ints, num_reals)
  end if

! Dump inpcrd output here, to be consistent with sander output:

  write(mdout, 1050)
  write(mdout, '(a80)') inpcrd_title
  write(mdout, '(t2,a,f10.3,a,/)') &
    'begin time read from input coords =', t, ' ps'

! MAINTENANCE WARNING!!! If axis flip optimization is ever supported for
! dipole code, the atm_inddip and atm_dipvel array values need to be flipped.

! If we are reading NMR restraints/weight changes, read them, and then determine
! how many of the torsional parameters are improper:

  if (nmropt .ne. 0) then
    call nmr_read(atm_crd, mdin, mdout)
    call set_num_improp_dihed(nphih, gbl_dihed, &
                              nphia, gbl_dihed(diheda_idx), nptra)
  endif

! Open the data dumping files and position it depending on the type of run:

  call open_output_files

! atm_itree is no longer needed, so deallocate:

  num_ints = num_ints - size(atm_itree)
  deallocate(atm_itree)

! We need atm_isymbl for igb = 8 so we can assign the gb_alpha/beta/gamma params

  if (.not. using_gb_potential .and. icnstph .ne. 2) then
    deallocate(atm_isymbl)
    num_ints = num_ints - size(atm_isymbl)
  end if

!Limited debugf namelist support
  call init_mdin_debugf_dat()
  call validate_mdin_debugf_dat()

  return

! Standard format statements:

! |= screen out in dacdif

 1000 format(/10x, 55('-'), /10x, &
             'Amber 14 SANDER                              2014', &
             /10x, 55('-')/)

#ifdef MPI
 1010 format('File Assignments:', /, 12('|', a7, ': ', a, /))
#else
 1010 format('File Assignments:', /, 11('|', a7, ': ', a, /))
#endif

 1040 format(80('-')/,'   2.  CONTROL  DATA  FOR  THE  RUN',/80('-')/)

 1050 format(/80('-')/,'   3.  ATOMIC COORDINATES AND VELOCITIES',/80('-')/)

end subroutine master_setup

!*******************************************************************************
!
! Subroutine:   open_output_files
!
! Description:  Routine to open the dumping and restart files.
!*******************************************************************************

subroutine open_output_files

  use bintraj_mod
  use file_io_mod
  use mdin_ctrl_dat_mod
  use prmtop_dat_mod

  implicit none

  character(10), parameter      :: file_version = '9.00'
  integer                       :: box_flag
  character(100)                :: bin4_title

! ioutfm .ne. 0 selects binary output, theoretically for all files below. In
! reality though, we never open mden for binary output.

  if (ioutfm .le. 0) then               ! Formatted dumping:

    if (ntwx .gt. 0) then
      call amopen(mdcrd, mdcrd_name, owrite, 'F', 'W')
      write(mdcrd, 1000) prmtop_ititl
    end if

    if (ntwv .gt. 0) then
      call amopen(mdvel, mdvel_name, owrite, 'F', 'W')
      write(mdvel, 1000) prmtop_ititl
    end if

    if (ntwf .gt. 0) then
      call amopen(mdfrc, mdfrc_name, owrite, 'F', 'W')
      write(mdfrc, 1000) prmtop_ititl
    end if

  else if (ioutfm .eq. 1) then

    call open_binary_files

  else if (ioutfm .eq. 2) then  ! The new "bin4" efficiency format...

    if (ntwx .gt. 0) then

      bin4_title = trim(mdcrd_name) // '.bin4'
      call amopen(mdcrd, bin4_title, owrite, 'U', 'W') 
      write(mdcrd) file_version
      write(mdcrd) prmtop_ititl

      if (ntb .gt. 0) then
        box_flag = 1
      else
        box_flag = 0
      end if

      if (ntwprt .ne. 0) then
        write(mdcrd) ntwprt, box_flag
      else
        write(mdcrd) natom, box_flag
      end if

    end if 

    if (ntwv .gt. 0) then

      bin4_title = trim(mdvel_name) // '.bin4'
      call amopen(mdvel, bin4_title, owrite, 'U', 'W')   
      write(mdvel) file_version
      write(mdvel) prmtop_ititl

      box_flag = 0

      if (ntwprt .ne. 0) then
        write(mdvel) ntwprt, box_flag
      else
        write(mdvel) natom, box_flag
      end if

    end if

  end if

! Open the energies file:

  if (ntwe .gt. 0) then
    call amopen(mden, mden_name, owrite, 'F', 'W')
  end if

! Open the restart file

  if (ntxo .le. 0) then
    call amopen(restrt, restrt_name, owrite, 'U', 'W')
  else if (ntxo .eq. 1) then
    call amopen(restrt, restrt_name, owrite, 'F', 'W')
  end if

! If we are doing MD, then the restrt file gets opened and then
! closed to force flushing the file write buffer.  For minimizations,
! however, we keep the file open and only write a restart file at the
! end, since no intermediate restarts are written.  Therefore, keep
! restrt open for minimizations

  if (imin .eq. 0) close(restrt)

! Open the AMD file:

  if (iamd .gt. 0) then
    call amopen(amdlog, amdlog_name, owrite, 'F', 'W') 
  end if

! Open the scaledMD file:

  if (scaledMD .gt. 0) then
    call amopen(scaledMDlog, scaledMDlog_name, owrite, 'F', 'W') 
  end if

! Open the Subunit file: Added by Lijiang Yang
  if (isub .gt. 0) then
    call amopen(sluene, sluene_name, owrite, 'F', 'W')
    call amopen(solene, solene_name, owrite, 'F', 'W')
    call amopen(interene, interene_name, owrite, 'F', 'W')
    write(sluene, 2000)
    write(sluene, '(a, /)') '| Solute Potential Energy by Lijiang Yang, PKU'
    write(solene, 2000)
    write(solene, '(a, /)') '| Solvent Potential Energy by Lijiang Yang, PKU'
    write(interene, 2000)
    write(interene, '(a, /)') '| Interaction Potential Energy by Lijiang Yang, PKU'
  end if

#ifdef MPI
! Open the mpi logfile:

  call amopen(logfile, logfile_name, owrite, 'F', 'W') 
#endif

  return

1000 format(a80)
2000 format(/10x, 55('-'), /10x, &
             'Amber 14 SANDER                              2014', &
             /10x, 55('-')/)

end subroutine open_output_files

!*******************************************************************************
!
! Subroutine:   printdefines
!
! Description:  Routine to print info about conditional compilation defines.
!               We just print defines with significant functional, performance,
!               or configurational significance.
!*******************************************************************************

subroutine printdefines()

  use file_io_mod
  use gbl_constants_mod

  implicit none

  write(mdout,'(a)') '| Conditional Compilation Defines Used:'

#ifdef HAS_10_12
  write(mdout, '(a)') '| HAS_10_12'
#endif

#ifdef MPI
  write(mdout, '(a)') '| MPI'
#endif

#ifdef FFTW_FFT
  write(mdout, '(a)') '| FFTW_FFT'
#endif

#ifdef PUBFFT
  write(mdout, '(a)') '| PUBFFT'
#endif

#ifdef BINTRAJ
  write(mdout, '(a)') '| BINTRAJ'
#endif

#ifdef MKL
  write(mdout, '(a)') '| MKL'
#endif

#ifdef CUDA
  write(mdout, '(a)') '| CUDA'
#endif

#ifdef EMIL
  write(mdout, '(a)') '| EMIL'
#endif

  write(mdout, *)

  return

end subroutine printdefines 

end module master_setup_mod
