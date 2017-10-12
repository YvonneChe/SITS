#include "copyright.i"
!*******************************************************************************
!
! Module: get_cmdline_mod
!
! Description: Module controlling getting command line arguments from either
!              a line in a group file or from the command-line itself
!              
!*******************************************************************************

module get_cmdline_mod

  implicit none

  private       add_suffix

contains

#ifdef USE_PXFGETARG
!*******************************************************************************
!
! Subroutine:  pmemd_getarg
!
! Description: wrapper for intrinsic getarg to make it work with nonstandard
!              systems
!
!*******************************************************************************

subroutine pmemd_getarg(i, c)

  implicit none

! Formal arguments:

  character(*)  :: c
  integer       :: i

! Local variables:

  integer       :: ilen
  integer       :: ierror

! Note that the pxfgetarg subroutine is broken in ifc 7.1; use getarg instead
! for any compiler that supports it!

  call pxfgetarg(i, c, ilen, ierror)

  if (ierror .ne. 0) then
    write(mdout, '(a)') 'ERROR: Command line parsing failed in getarg!'
    call mexit(6, 1)
  end if

  return

end subroutine pmemd_getarg
#else
! The preferred command line subroutine is getarg.
#define pmemd_getarg    getarg
#endif /* USE_PXFGETARG */

!*******************************************************************************
!
! Subroutine:  get_cmdline
!
! Description: Gets command-line arguments either from CL or from groupfile
!              
! OUTPUT: File names and command-line arguments
!
! Author: George Seibel. Modifications by JMS to handle strings for multipmemd
!
!*******************************************************************************

subroutine get_cmdline(terminal_flag)

  use file_io_dat_mod
  use gbl_constants_mod, only : error_hdr, VERSION
  use pmemd_lib_mod
  use parallel_dat_mod
#ifdef MPI
  use remd_mod,         only  : remd_method
#endif

  implicit none

! Passed variables:

  logical, intent(in out)  :: terminal_flag

! Local variables:

  character(max_fn_len) :: arg  ! temp buffer for each of the whitespace-
                        ! delimited command line words.
  integer       :: iarg ! arg pointer, final number of arguments.
  integer       :: indx

#ifdef MPI
  logical       :: in_groupfile ! if we're reading this from the groupfile
                                ! or _really_ on the command-line
#endif

! Flags that determine if/when the -suffix gets used:

  logical       :: mdout_specified
  logical       :: mdinfo_specified
  logical       :: mdcrd_specified
  logical       :: mdvel_specified
  logical       :: mdfrc_specified
  logical       :: mden_specified
  logical       :: restrt_specified
  logical       :: logfile_specified
  logical       :: cpout_specified
  logical       :: cprestrt_specified
  logical       :: outfile_suffix_specified
  logical       :: amdlog_specified
  logical       :: scaledMDlog_specified
  logical       :: always_apply_suffix
! Added by Lijiang Yang for Integrated Tempering Sampling
  logical       :: fbinput_specified
  logical       :: normlinput_specified
  logical       :: fbrestrt_specified
  logical       :: normlrestrt_specified
  logical       :: fbtrj_specified
  logical       :: normltrj_specified
! ------------------------------------------------------
! Added by Lijiang Yang for Subunit Energy
  logical       :: sluene_specified
  logical       :: solene_specified
  logical       :: interene_specified
! ------------------------------------------------------

  integer       :: suffix_len

  mdout_specified = .false.
  mdinfo_specified = .false.
  mdcrd_specified = .false.
  mdvel_specified = .false.
  mdfrc_specified = .false.
  mden_specified = .false.
  restrt_specified = .false.
  logfile_specified = .false.
  cpout_specified = .false.
  cprestrt_specified = .false.
  outfile_suffix_specified = .false.
  amdlog_specified = .false.
  scaledMDlog_specified = .false.
  always_apply_suffix = .false.
! Added by Lijiang Yang for Integrated Tempering Sampling
  fbinput_specified = .false.
  normlinput_specified = .false.
  fbrestrt_specified = .false.
  normlrestrt_specified = .false.
  fbtrj_specified = .false.
  normltrj_specified = .false.
! ------------------------------------------------------
! Added by Lijiang Yang for Subunit Energy
  sluene_specified = .false.
  solene_specified = .false.
  interene_specified = .false.
! ------------------------------------------------------

! Default file names:

  mdin_name   = 'mdin'
  mdout_name  = 'mdout'
  inpcrd_name = 'inpcrd'
  prmtop_name = 'prmtop'
  restrt_name = 'restrt'
  refc_name   = 'refc'
  mdvel_name  = 'mdvel'
  mdfrc_name  = 'mdfrc'
  mden_name   = 'mden'
  mdcrd_name  = 'mdcrd'
  mdinfo_name = 'mdinfo'
  cpin_name   = 'cpin'
  cpout_name  = 'cpout'
  cprestrt_name = 'cprestrt'
  logfile_name = 'logfile'
! ------------------------------------------------------
! Added by Lijiang Yang for Subunit Energy
  sluene_name = 'sluene'
  solene_name = 'solene'
  interene_name = 'interene'
  fbinput_name = 'fbinput'
  normlinput_name = 'normlinput'
  fbrestrt_name = 'fbrestrt'
  normlrestrt_name = 'normlrestrt'
  fbtrj_name = 'fbtrj'
  normltrj_name = 'normltrj'
! ------------------------------------------------------
#ifdef MPI
  proc_map_name = ' '

! To make sure we have unique output file names in all independent replicas,
! activate the outfile_suffix_specified, and build the suffix from the replica
! number (i.e. .001, .002, etc.). Activate the outfile_suffix_specified here,
! and if -suffix is provided, take that as a key to apply it to all files, even
! the ones we DID specify in the groupfile.

  if ( numgroups .gt. 1) then
    outfile_suffix_specified = .true.
    call assign_group_suffix(outfile_suffix)
  end if

  in_groupfile = len_trim(groupline_buffer) .ne. 0
#endif /* MPI */


! Default status of output: new

  owrite = 'N'

! Get com line arguments:

  iarg = 0
  indx = iargc_wrap()

  if (indx .eq. 0) goto 20

  10 continue

  iarg = iarg + 1

  call getarg_wrap(iarg, arg)

  ! Look for a --version or -V flag and dump out the version of pmemd
  if (arg .eq. '-V' .or. arg .eq. '--version') then
#ifdef MPI
    if (in_groupfile) then
      write(0, '(2a)') trim(arg), &
                     ' Bad flag in groupfile. Must be on command-line!'
      call mexit(mdout, 0)
    end if
#endif
    call getarg_wrap(0, arg)
    call basename(arg)
    write(6, '(3a)') trim(arg), ': ', trim(VERSION)
    terminal_flag = .true.
    return ! no need to do anything more here.
  else if (arg .eq. '-O') then
    owrite = 'U'
  else if (arg .eq. '-i') then
    iarg = iarg + 1
    call getarg_wrap(iarg, mdin_name)
  else if (arg .eq. '-o') then
    iarg = iarg + 1
    call getarg_wrap(iarg, mdout_name)
    mdout_specified = .true.
  else if (arg .eq. '-p') then
    iarg = iarg + 1
    call getarg_wrap(iarg, prmtop_name)
  else if (arg .eq. '-c') then
    iarg = iarg + 1
    call getarg_wrap(iarg, inpcrd_name)
  else if (arg .eq. '-r') then
    iarg = iarg + 1
    call getarg_wrap(iarg, restrt_name)
    restrt_specified = .true.
  else if (arg .eq. '-ref' .or. arg .eq.'-z') then
    iarg = iarg + 1
    call getarg_wrap(iarg, refc_name)
  else if (arg .eq. '-e') then
    iarg = iarg + 1
    call getarg_wrap(iarg, mden_name)
    mden_specified = .true.
  else if (arg .eq. '-v') then
    iarg = iarg + 1
    call getarg_wrap(iarg, mdvel_name)
    mdvel_specified = .true.
  else if (arg .eq. '-frc') then
    iarg = iarg + 1
    call getarg_wrap(iarg, mdfrc_name)
    mdfrc_specified = .true.
  else if (arg .eq. '-x' .or. arg .eq.'-t') then
    iarg = iarg + 1
    call getarg_wrap(iarg, mdcrd_name)
    mdcrd_specified = .true.
  else if (arg .eq. '-inf') then
    iarg = iarg + 1
    call getarg_wrap(iarg, mdinfo_name)
    mdinfo_specified = .true.
  else if (arg .eq. '-l') then
    iarg = iarg + 1
    call getarg_wrap(iarg, logfile_name)
    logfile_specified = .true.
! Added by Lijiang Yang for Integrated Tempering Sampling
  else if (arg .eq. '-cfb') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, fbinput_name)
    fbinput_specified = .true.
  else if (arg .eq. '-cnorml') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, normlinput_name)
    normlinput_specified = .true.
  else if (arg .eq. '-rfb') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, fbrestrt_name)
    fbrestrt_specified = .true.
  else if (arg .eq. '-rnorml') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, normlrestrt_name)
    normlrestrt_specified = .true.
  else if (arg .eq. '-xfb') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, fbtrj_name)
    fbtrj_specified = .true.
  else if (arg .eq. '-xnorml') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, normltrj_name)
    normltrj_specified = .true.
! ------------------------------------------------------
! Added by Lijiang Yang for Subunit Energy
  else if (arg .eq. '-sluene') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, sluene_name)
    sluene_specified = .true.
  else if (arg .eq. '-solene') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, solene_name)
    solene_specified = .true.
  else if (arg .eq. '-interene') then
    iarg = iarg + 1
    call pmemd_getarg(iarg, interene_name)
    interene_specified = .true.
! ------------------------------------------------------
  else if (arg .eq. '-suffix') then
    iarg = iarg + 1
    call getarg_wrap(iarg, outfile_suffix)
    if (numgroups .gt. 1) always_apply_suffix = .true.
    outfile_suffix_specified = .true.
  else if (arg .eq. '-cpin') then
    iarg = iarg + 1
    call getarg_wrap(iarg, cpin_name)
  else if (arg .eq. '-cpout') then
    iarg = iarg + 1
    call getarg_wrap(iarg, cpout_name)
    cpout_specified = .true.
  else if (arg .eq. '-cprestrt') then
    iarg = iarg + 1
    call getarg_wrap(iarg, cprestrt_name)
    cprestrt_specified = .true.
  else if (arg .eq. '-amd') then
    iarg = iarg + 1
    call getarg_wrap(iarg, amdlog_name)
    amdlog_specified = .true.
  else if (arg .eq. '-scaledMD') then
    iarg = iarg + 1
    call getarg_wrap(iarg, scaledMDlog_name)
    scaledMDlog_specified = .true.
  else if (arg .eq. '-help' .or. arg .eq. '--help' .or. arg .eq. '-H') then
#ifdef MPI
    if (in_groupfile) then
      write(0, '(2a)') trim(arg), &
                     ' Bad flag in groupfile. Must be on command-line!'
      call mexit(mdout, 0)
    end if
#endif
    write(mdout, 9000)
    terminal_flag = .true.
    return

#ifdef MPI
  else if (arg .eq. '-ng') then
    iarg = iarg + 1
    call getarg_wrap(iarg, arg)
    read(arg, '(1i5)') numgroups
  else if (arg .eq. '-rem') then
    ! Make sure we didn't call 
    iarg = iarg + 1
    call getarg_wrap(iarg, arg)
    read(arg, '(1i5)', err=666) remd_method
    if (remd_method .ne. -1 .and. len_trim(remd_dimension_name) .ne. 0) goto 667
  else if (arg .eq. '-remd-file') then
    iarg = iarg + 1
    call getarg_wrap(iarg, remd_dimension_name)
    if (remd_method .ne. 0) goto 667
    ! Make remd_method -1 here, since that's multi-D REMD indicator
    remd_method = -1
  else if (arg .eq. '-remlog') then
    iarg = iarg + 1
    call getarg_wrap(iarg, remlog_name)
  else if (arg .eq. '-remtype') then
!   Not supported yet!
!   iarg = iarg + 1
!   call getarg_wrap(iarg, remlog_name)
    write(0, *) '-remtype not supported in pmemd.MPI. Use sander.MPI instead!'
    call mexit(mdout, 1)
  else if (arg .eq. '-groupfile') then
    iarg = iarg + 1
    call getarg_wrap(iarg, groupfile_name)
  else if (arg .eq. '-ng-nonsequential') then
    ng_nonsequential = .true.
  else if (arg .eq. '-gpes') then
    iarg = iarg + 1
    call getarg_wrap(iarg, proc_map_name)
  else if (arg .eq. '-p4pg') then
    iarg = iarg + 1
  else if (arg .eq. '-p4wd') then
    iarg = iarg + 1
  else if (arg .eq. '-np') then
    iarg = iarg + 1
  else if (arg .eq. '-mpedbg') then
    continue
  else if (arg .eq. '-dbx') then
    continue
  else if (arg .eq. '-gdb') then
    continue
#endif
  else
    if (arg .eq. ' ') goto 20
    write(mdout, '(/,2x,a,a)') 'unknown flag: ', arg
    write(mdout, 9000)
    call mexit(6, 1)
  end if

  if (iarg .lt. indx) goto 10

  20 continue

  if (outfile_suffix_specified) then

    suffix_len = len_trim(outfile_suffix)

    if (.not. mdout_specified .or. always_apply_suffix) &
      call add_suffix(mdout_name, outfile_suffix, suffix_len)

    if (.not. mdinfo_specified .or. always_apply_suffix) &
      call add_suffix(mdinfo_name, outfile_suffix, suffix_len)

    if (.not. mdcrd_specified .or. always_apply_suffix) &
      call add_suffix(mdcrd_name, outfile_suffix, suffix_len)

    if (.not. mdvel_specified .or. always_apply_suffix) &
      call add_suffix(mdvel_name, outfile_suffix, suffix_len)

    if (.not. mdfrc_specified .or. always_apply_suffix) &
      call add_suffix(mdfrc_name, outfile_suffix, suffix_len)

    if (.not. mden_specified .or. always_apply_suffix) &
      call add_suffix(mden_name, outfile_suffix, suffix_len)

    if (.not. restrt_specified .or. always_apply_suffix) &
      call add_suffix(restrt_name, outfile_suffix, suffix_len)

    if (.not. logfile_specified .or. always_apply_suffix) &
      call add_suffix(logfile_name, outfile_suffix, suffix_len)

    if (.not. cpout_specified .or. always_apply_suffix) &
      call add_suffix(cpout_name, outfile_suffix, suffix_len)

    if (.not. cprestrt_specified .or. always_apply_suffix) &
      call add_suffix(cprestrt_name, outfile_suffix, suffix_len)

    if (.not. amdlog_specified .or. always_apply_suffix) &
      call add_suffix(amdlog_name, outfile_suffix, suffix_len)

    if (.not. scaledMDlog_specified .or. always_apply_suffix) &
      call add_suffix(scaledMDlog_name, outfile_suffix, suffix_len)

  end if

  return
#ifdef MPI
9000 format(/, 2x, &
      'usage: pmemd  [-O] -i mdin -o mdout -p prmtop -c inpcrd -r restrt', &
      /14x, '[-ref refc -x mdcrd -v mdvel -frc mdfrc -e mden &
             &-inf mdinfo -l logfile]', &
      /14x, '[-ng numgroups -groupfile groupfile -rem remd_method]', &
      /14x, '[-amd amdlog_name -scaledMD scaledMDlog_name &
             &-suffix output_files_suffix]',/)

666 write(mdout, '(2a)') error_hdr, 'bad integer read on command-line (-rem)'
    call mexit(mdout, 1)

667 write(mdout, '(2a)') error_hdr, &
      'do not specify -rem # when you supply a -remdim <dimension file>'
    call mexit(mdout, 1)
#else
9000 format(/, 2x, &
      'usage: pmemd  [-O] -i mdin -o mdout -p prmtop -c inpcrd -r restrt', &
      /14x, '[-ref refc -x mdcrd -v mdvel -frc mdfrc -e mden &
             &-inf mdinfo -l logfile]', &
      /14x, '[-amd amdlog_name -scaledMD scaledMDlog_name &
             &-suffix output_files_suffix]',/)
#endif /* MPI */

end subroutine get_cmdline

!*******************************************************************************
!
! Subroutine:  add_suffix
!
! Description: adds a suffix to each file
!              
!*******************************************************************************

subroutine add_suffix(file_name, suffix, suffix_len)

  use file_io_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  character(max_fn_len), intent(in out) :: file_name
  character(max_fn_len), intent(in)     :: suffix
  integer, intent(in)                   :: suffix_len

  if (suffix_len + len_trim(file_name) .lt. max_fn_len) then
    file_name = trim(file_name) // '.' // trim(suffix)
  else
    write(mdout, '(a)') 'ERROR: Filename with suffix too long!'
    call mexit(6, 1)
  end if

  return

end subroutine add_suffix

!*******************************************************************************
!
! Function: iargc_wrap
!
! Description: wrapper for iargc() intrinsic so we can read arguments from a 
!              groupfile. Counts number of arguments in groupfile line or on
!              command-line
!
!*******************************************************************************

integer function iargc_wrap()

#ifdef MPI
  use file_io_dat_mod,  only : groupline_buffer
  use pmemd_lib_mod,    only : get_num_tokens
#endif

  implicit none

#ifdef USE_PXFGETARG
  integer  :: ipxfargc
#else
  integer  :: iargc  ! intrinsic iargc if we're on command-line
#endif

#ifdef MPI
  ! If numgroups > 1 and groupfile has been found, parse groupline

  if (len_trim(groupline_buffer) .ne. 0) then
    call get_num_tokens(groupline_buffer, iargc_wrap)
  else
# ifdef USE_PXFGETARG
    iargc_wrap = ipxfargc()
# else
    iargc_wrap = iargc()
# endif
  end if

#else

# ifdef USE_PXFGETARG
  iargc_wrap = ipxfargc()
# else
  iargc_wrap = iargc()
# endif /* USE_PXFGETARG */

#endif /* MPI */

end function iargc_wrap

!*******************************************************************************
!
! Subroutine: getarg_wrap
!
! Description: Wrapper subroutine for intrinsic getarg to support parsing
!              arguments from a groupfile
!
!*******************************************************************************

subroutine getarg_wrap(iarg, arg)

#ifdef MPI
  use file_io_dat_mod, only : groupline_buffer
  use pmemd_lib_mod, only : get_token
#endif

  implicit none

! Passed variables

  integer, intent(in)       :: iarg  ! argument number

  character(*), intent(out) :: arg   ! argument string

#ifdef MPI
  if (len_trim(groupline_buffer) .gt. 0) then
    call get_token(groupline_buffer, iarg, arg)
  else
    call pmemd_getarg(iarg, arg)
  end if
#else
  call pmemd_getarg(iarg, arg)
#endif /* MPI */

end subroutine getarg_wrap

end module get_cmdline_mod
