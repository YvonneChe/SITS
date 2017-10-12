!*******************************************************************************
!
! Module: file_io_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module file_io_dat_mod

  implicit none

! File-related data.  No need to broadcast this stuff, since only the master
! does i/o:

  integer, parameter            :: max_fn_len = 256
#ifdef MPI
  integer, parameter            :: grpfl_line_len = 4096
#endif

  character(max_fn_len), save   :: mdin_name
  character(max_fn_len), save   :: mdout_name
  character(max_fn_len), save   :: mdinfo_name
  character(max_fn_len), save   :: prmtop_name
  character(max_fn_len), save   :: inpcrd_name
  character(max_fn_len), save   :: refc_name
  character(max_fn_len), save   :: mdcrd_name
  character(max_fn_len), save   :: mdvel_name
  character(max_fn_len), save   :: mdfrc_name
  character(max_fn_len), save   :: mden_name
  character(max_fn_len), save   :: restrt_name
  character(max_fn_len), save   :: logfile_name
! Added by Lijiang Yang for Integrated Tempering Sampling
  character(max_fn_len), save   :: fbinput_name
  character(max_fn_len), save   :: normlinput_name
  character(max_fn_len), save   :: fbrestrt_name
  character(max_fn_len), save   :: normlrestrt_name
  character(max_fn_len), save   :: fbtrj_name
  character(max_fn_len), save   :: normltrj_name
! ------------------------------------------------------
! Added by Lijiang Yang for Subunit Energy
  character(max_fn_len), save   :: sluene_name
  character(max_fn_len), save   :: solene_name
  character(max_fn_len), save   :: interene_name
! ------------------------------------------------------
  character(max_fn_len), save   :: infile_suffix
  character(max_fn_len), save   :: outfile_suffix
  character(max_fn_len), save   :: cpin_name
  character(max_fn_len), save   :: cpout_name
  character(max_fn_len), save   :: cprestrt_name
!AMD file for reweighting values
  character(max_fn_len), save   :: amdlog_name = 'amd.log'
!scaledMD file for reweighting values
  character(max_fn_len), save   :: scaledMDlog_name = 'scaledMD.log'
#ifdef MPI
! Set file names for multipmemd and REMD. proc_map refers
! to the file describing how many processors to assign to
! each communicator
  character(max_fn_len), save   :: groupfile_name = ' '
  character(max_fn_len), save   :: remlog_name = 'rem.log'
  character(max_fn_len), save   :: remtype_name = 'rem.type' ! not supported yet
  character(max_fn_len), save   :: proc_map_name
  character(max_fn_len), save   :: remd_dimension_name = ' '
  character(grpfl_line_len), save :: groupline_buffer = ' '
  character(grpfl_line_len), save :: proc_map_buffer = ' '
#endif
  character, save       :: owrite

! Logical unit numbers for pmemd files.

  integer, parameter    :: mdin        =  5
  integer, parameter    :: mdout       =  6
  integer, parameter    :: mdinfo      =  7
  integer, parameter    :: prmtop      =  8
  integer, parameter    :: inpcrd      =  9
  integer, parameter    :: refc        = 10
  ! who was 11?
  integer, parameter    :: mdcrd       = 12
  integer, parameter    :: mdvel       = 13
  integer, parameter    :: mdfrc       = 14
  integer, parameter    :: mden        = 15
  integer, parameter    :: restrt      = 16
  integer, parameter    :: restrt2     = 17
  integer, parameter    :: logfile     = 18
#ifdef MPI
  integer, parameter    :: groupfile     = 19
  integer, parameter    :: remlog        = 20
  integer, parameter    :: proc_map      = 21
  integer, parameter    :: remd_file     = 22
#endif
  integer, parameter    :: cpin          = 23
  integer, parameter    :: cpout         = 24
  integer, parameter    :: cprestrt      = 25

#ifdef MPI
! Maybe not quite the right place to put this, but it can't be in
! multipmemd_mod, since we'll get a cyclic dependency.

  logical, save         :: ng_nonsequential = .false.
#endif
  integer, parameter    :: amdlog        = 23
  integer, parameter    :: scaledMDlog        = 24

! Added by Lijiang Yang for Integrated Tempering Sampling
  integer, parameter    :: fbinput        = 25
  integer, parameter    :: normlinput     = 26
  integer, parameter    :: fbrestrt       = 27
  integer, parameter    :: normlrestrt    = 28
  integer, parameter    :: fbtrj          = 29
  integer, parameter    :: normltrj       = 30
! ------------------------------------------------------
! Added by Lijiang Yang for Subunit Energy
  integer, parameter    :: sluene         = 31
  integer, parameter    :: solene         = 32
  integer, parameter    :: interene       = 33
! ------------------------------------------------------

end module file_io_dat_mod
