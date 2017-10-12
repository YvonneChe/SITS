!*******************************************************************************
!
! Module: mdin_subunit_dat_mod
!
! Description: Added by Lijiang Yang for Subunit Energy
!              
!*******************************************************************************

module mdin_subunit_dat_mod

use file_io_dat_mod

  implicit none

!  bgsluatm = 1 !The first atom index of solute
!  nsluatm = 0  !Number of atoms of solute
!  bgsolatm = 1 !The first atom index of solvent

  integer, parameter    :: mdin_subunit_int_cnt = 5

  integer                       bgsluatm, nsluatm, bgsolatm, nsolatm,&
                                solene_avg_sampling

  common / mdin_subunit_int /   bgsluatm, nsluatm, bgsolatm, nsolatm,&
                                solene_avg_sampling

  save  :: /mdin_subunit_int /

  private       :: subunit

  namelist /subunit/        bgsluatm, nsluatm, bgsolatm, nsolatm, solene_avg_sampling

! Private subroutines:


contains

!*******************************************************************************
!
! Subroutine:  init_mdin_subunit_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine init_mdin_subunit_dat(natom)

  use file_io_mod
  use gbl_constants_mod
  use mdin_ctrl_dat_mod

  implicit none

  integer               natom
! Local variables:

  integer               :: ifind
  integer               :: inerr

  inerr = 0

  bgsluatm = 1
  nsluatm  = 0
  bgsolatm = 1
  nsolatm  = natom - nsluatm
  solene_avg_sampling = -1      ! -1 means "set to default value"

  rewind(mdin)                  ! Insurance against maintenance mods.

  call nmlsrc('subunit', mdin, ifind)

  if (ifind .ne. 0) read(mdin, nml = subunit)   ! Namelist found. Read it:
  nsolatm  = natom - nsluatm
  if (solene_avg_sampling .eq. -1) then
    if (ntpr .gt. 0) then
      solene_avg_sampling = ntpr
    else
      solene_avg_sampling = 1
    end if
  end if

  if (solene_avg_sampling .le. 0) then
    write(mdout, '(a,a)') error_hdr, 'solene_avg_sampling must be .gt. 0!'
  end if

  return

end subroutine init_mdin_subunit_dat

!*******************************************************************************
!
! Subroutine:  validate_mdin_subunit_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine validate_mdin_subunit_dat

!  use mdin_ctrl_dat_mod

  implicit none

! Formal arguments:

! Local variables

  integer       :: inerr

! Check on bogus data and unsupported options:

  inerr = 0


! Field any errors and bag out.

  if (inerr .eq. 1) then
    write(mdout, '(/,a)') ' Input errors occurred. Terminating execution.'
    call mexit(6, 1)
  else
    write(mdout, '(a)') ' '
  end if

  return

end subroutine validate_mdin_subunit_dat

!*******************************************************************************
!
! Subroutine:  print_mdin_subunit_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine print_mdin_subunit_dat

!  use prmtop_dat_mod

  implicit none

! Formal arguments:

! Local variables:

  write(mdout,'(/a)') 'Subunit Partition:'
  write(mdout,'(5x,4(a,i8))') 'bgsluatm =', bgsluatm, &
          ', nsluatm  =', nsluatm, ', bgsolatm =', bgsolatm, ', nsolatm  =', nsolatm
  write(mdout,'(5x,(a,i8))') 'solene_avg_sampling =', solene_avg_sampling

  return

end subroutine print_mdin_subunit_dat

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_mdin_subunit_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_mdin_subunit_dat

  use parallel_dat_mod

  implicit none

  call mpi_bcast(bgsluatm, mdin_subunit_int_cnt, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  return

end subroutine bcast_mdin_subunit_dat
#endif

end module mdin_subunit_dat_mod
