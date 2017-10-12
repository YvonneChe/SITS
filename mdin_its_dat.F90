!*******************************************************************************
!
! Module: mdin_its_dat_mod
!
! Description: Added by Lijiang Yang for Integrated Tempering Sampling 
!              
!*******************************************************************************

module mdin_its_dat_mod

use file_io_dat_mod

  implicit none

! ntemp = 100 number of temperature used in ITS.
! templ = 283 the lowest temperature used in ITS.
! temph = 600 the highest temperature used in ITS. 
! peshift = 0 the shiftting energy to make the potential energy dealing in the ITS be larger than zero
! start_cycle = 0 how many cycles ITS has been run
! update_step = 100 how many steps to update fb once
! fb_bias = 0.0 factor to enhance the sampling in the certein energy range (negative to the high energy, positive to the low energy)
! fb_rest = 0 whether it is the first time running the ITS
! fb_const = 0 whether using the constant fb
! rb_fac1 = 0.5 rb updating strategy: (rbfb_i + rbfb_i+1)*0.5 
! rb_fac2 = 0.0 rb updating strategy: mcycle*0.0 
! fb_init = -0.005 fb initialization seed: fb0=exp(fb_init*ntemp)

  integer, parameter    :: mdin_its_int_cnt = 5

  integer                       ntemp, start_cycle, update_step, fb_rest, fb_const

  common / mdin_its_int /       ntemp, start_cycle, update_step, fb_rest, fb_const

  save  :: /mdin_its_int /

  integer, parameter    :: mdin_its_dbl_cnt = 7

  double precision              templ, temph, peshift, fb_bias, rb_fac1, rb_fac2, fb_init

  common / mdin_its_dbl /       templ, temph, peshift, fb_bias, rb_fac1, rb_fac2, fb_init

  save  :: /mdin_its_dbl /

  ! 1 cal = 4.184 joules  beta for 300K: 4.184*1000/(6.022*10E23*300*1.380653*10E-23)=1.677
!  double precision, parameter   :: tempcoff = 503.229837949d0

  private       :: its

  namelist /its/        ntemp, start_cycle, update_step, fb_rest, fb_const, &
                        templ, temph, peshift, fb_bias, rb_fac1, rb_fac2, fb_init 

! Private subroutines:


contains

!*******************************************************************************
!
! Subroutine:  init_mdin_its_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine init_mdin_its_dat

  use file_io_mod
  use mdin_ctrl_dat_mod

  implicit none

! Local variables:

  integer               :: ifind
  integer               :: inerr

  inerr = 0

  ntemp = 100
  start_cycle = 0
  update_step = 100
  fb_rest = 0
  fb_const = 0
  templ = 283.0d0
  temph = 600.0d0
  peshift = 0.0d0
  fb_bias = 0.0d0
  rb_fac1 = 0.5d0
  rb_fac2 = 0.0d0
  fb_init = -0.005d0


  rewind(mdin)                  ! Insurance against maintenance mods.

  call nmlsrc('its', mdin, ifind)

  if (ifind .ne. 0) read(mdin, nml = its)   ! Namelist found. Read it:

  return

end subroutine init_mdin_its_dat

!*******************************************************************************
!
! Subroutine:  validate_mdin_its_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine validate_mdin_its_dat

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

end subroutine validate_mdin_its_dat

!*******************************************************************************
!
! Subroutine:  print_mdin_its_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine print_mdin_its_dat

!  use prmtop_dat_mod

  implicit none

! Formal arguments:

! Local variables:

  write(mdout,'(/a)') 'ITS parameters:'
!  write(mdout,'(5x,4(a,i8))') 'ntemp =', ntemp, &
!          ', start_cycle =', start_cycle, ', update_step  =', update_step, &
!          ', fb_rest =', fb_rest
  write(mdout,'(5x,4(a,i8))') 'ntemp =', ntemp, &
          ', update_step  =', update_step, ', fb_rest =', fb_rest, &
          ', fb_const =', fb_const
  write(mdout, '(5x,3(a,f9.3),a,f10.1)') 'templ =', templ, ', temph =', temph, &
           ', fb_bias =', fb_bias, ', peshift =', peshift
  write(mdout, '(5x,(a,f9.3))') 'fb_init =', fb_init

  return

end subroutine print_mdin_its_dat

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_mdin_its_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_mdin_its_dat

  use parallel_dat_mod

  implicit none

  call mpi_bcast(ntemp, mdin_its_int_cnt, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)
  call mpi_bcast(templ, mdin_its_dbl_cnt, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)
  return

end subroutine bcast_mdin_its_dat
#endif

end module mdin_its_dat_mod
