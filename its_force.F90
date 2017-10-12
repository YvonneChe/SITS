#include "copyright.i"

!*******************************************************************************
!
! Module: its_force_mod
!
! Description: Added by Lijiang Yang for Integrated Tempering Sampling 
!              
!*******************************************************************************

module its_force_mod


  implicit none

    integer, save          :: mcycle
    double precision, save :: tempcoff
    double precision, save :: beta0
    double precision, save :: dtemp
    double precision, save :: gfe
    double precision, save :: gff
    double precision, save :: gftt
    double precision, save :: betaht
    double precision, save :: betalt
    double precision, allocatable, save :: mybeta(:)
    double precision, allocatable, save :: fb(:)
    double precision, allocatable, save :: gf(:)
    double precision, allocatable, save :: gft(:)
    double precision, allocatable, save :: ratio(:)
    double precision, allocatable, save :: pratio(:)
    double precision, allocatable, save :: norml(:)
    double precision, allocatable, save :: normlold(:)
    double precision, allocatable, save :: rb(:)
    double precision, allocatable, save :: rbfb(:)


contains

!*******************************************************************************
!
! Subroutine:  its_init
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine its_init(num_reals)

  use mdin_ctrl_dat_mod, only: temp0
  use mdin_its_dat_mod
  use file_io_dat_mod
  use file_io_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:
  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_reals


! Local variables:

  integer                       :: i

! begin of code
  
  call  alloc_its_force_mem(num_reals)

  tempcoff = 503.229837949  ! 1 cal = 4.184 joules  tempcoff=4184.0d0/(6.022*1.380653)
  beta0 = tempcoff/temp0    ! for temp0=300K beta0=4.184*1000/(6.022*10E23*300*1.380653*10E-23)=1.677
  betaht=tempcoff/temph
  betalt=tempcoff/templ
  dtemp=betalt-betaht
  mcycle = start_cycle

  do i=1,ntemp
     mybeta(i)=betalt-(i-1)*dtemp/(ntemp-1)
  enddo
  
  if ( fb_rest .lt. 1 ) then
     do i=1,ntemp-1
        fb(i) = exp(fb_init*ntemp)
        norml(i) = 0.0d0 
     enddo
     fb(ntemp) = exp(fb_init*ntemp)
  else
     call amopen(fbinput, fbinput_name, 'O', 'F', 'R')
     do i = 1, ntemp
        read(fbinput, *) fb(i)
     enddo
     close(fbinput)
     call amopen(normlinput, normlinput_name, 'O', 'F', 'R')
     do i = 1, ntemp-1
        read(normlinput, *) norml(i)
     enddo
     close(normlinput)
  endif
     
  return

1001 continue

  write(mdout, '(a,a)') 'FATAL: Could not read fb from ', &
                        fbinput_name
  call mexit(6, 1)

1002 continue

  write(mdout, '(a,a)') 'FATAL: Could not read normal from ', &
                        normlinput_name
  call mexit(6, 1)

end subroutine its_init

!*******************************************************************************
!
! Subroutine:  alloc_its_force_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_its_force_mem(num_reals)

  use mdin_its_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_reals


! Local variables:

  integer               :: i, alloc_failed

! begin of code

  ! allocate array
  allocate(      fb( ntemp     ), &
              norml( ntemp - 1 ), &
           normlold( ntemp - 1 ), &
             mybeta( ntemp     ), &
                 gf( ntemp     ), &
              ratio( ntemp     ), &
             pratio( ntemp     ), &
                 rb( ntemp     ), &
               rbfb( ntemp     ), &
                gft( ntemp     ), &
             stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals   + size(mybeta)   + size(fb)     + size(gf) +       &
              size(gft)   + size(ratio)    + size(pratio) +                  &
              size(norml) + size(normlold) + size(rb)     +                  &
              size(rbfb)

  return

end subroutine  alloc_its_force_mem

!*******************************************************************************
!
! Subroutine:  its_updatefb
!
! Description: update the neighbering sampling ratio every update_step step 
!              
!*******************************************************************************

subroutine its_updatefb

  use mdin_its_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:


! Local variables:

  integer                       :: i,alloc_failed
  double precision,allocatable  :: value1(:)

! begin of code
  mcycle = mcycle + 1
  allocate(  value1( ntemp ), &
             stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  do i = 1, ntemp-1
     rb(i) = ( rbfb(i) + rbfb(i+1) )*0.5
     ratio(i) = fb(i) - fb(i+1)
     normlold(i) = norml(i)
     value1(i)=rbfb(i+1)-rbfb(i)+fb_bias+rb(i)
  enddo

  if ( fb_rest .eq. 0 .and. mcycle .eq. 1 ) then
     do i = 1, ntemp-1
        norml(i) = rb(i)
        normlold(i) = -10.0d0**10
     enddo
  else 
     do i = 1, ntemp-1
        if ( norml(i) .gt. rb(i) ) then
           norml(i) = norml(i) + log( 1.0d0 + exp( rb(i) -  norml(i) ) )
        else
           norml(i) = rb(i) + log( 1.0d0 + exp( -rb(i) + norml(i) ) )
        endif
     enddo
  endif

  pratio(1) = 0.0d0
  fb(1) = -pratio(1)
  do i = 1, ntemp-1
      if ( normlold(i) .gt. value1(i) ) then
         ratio(i) = ratio(i) + normlold(i) - norml(i) + log( 1.0d0 + exp( -normlold(i) + value1(i) ) )
      else
         ratio(i) = value1(i) + ratio(i) - norml(i) + log( 1.0d0 + exp( normlold(i) - value1(i) ) )
      endif
     pratio( i + 1 ) = pratio(i) + ratio(i)
     fb( i + 1 ) = -pratio( i + 1 )
  enddo        

  do i = 1, ntemp
     rbfb(i) = 0.0d0
  enddo
  
  deallocate (value1)

  return

end subroutine its_updatefb


!*******************************************************************************
!
! Subroutine:  its_force
!
! Description: 
!              
!*******************************************************************************

subroutine its_force(f_its, pot_sluene, pot_interene, nstep) 

  use mdin_its_dat_mod

  implicit none

! Formal arguments:

  integer                       :: nstep
  double precision              :: f_its
  double precision              :: pot_sluene
  double precision              :: pot_interene


! Local variables:

  integer                       :: i

! begin of code
  gfe = pot_sluene + 0.5*pot_interene - peshift  
  do i = 1, ntemp
     gft(i) = -mybeta(i)*gfe + fb(i)
     gf(i) = -mybeta(i)*gfe + fb(i) + log( mybeta(i) )
  enddo
 
  gff = gf(1)
  gftt = gft(1)
  do i = 2, ntemp
     if ( gff .gt. gf(i) ) then
        gff = gff + log( 1.0d0 + exp( gf(i)-gff ) )
     else
        gff = gf(i) + log( 1.0d0 + exp( gff-gf(i) ) )
     endif
 
     if ( gftt .gt. gft(i) ) then
        gftt = gftt + log( 1.0d0 + exp( gft(i)-gftt ) )
     else
        gftt = gft(i) + log( 1.0d0 + exp( gftt-gft(i) ) )
     endif
  enddo

  if ( fb_const .ne. 1) then
     if ( nstep .eq. 1 .or. mod(nstep,update_step) .eq. 0 ) then
        do i = 1, ntemp
           rbfb(i) = gft(i) !from Shao
        enddo
     else
        do i = 1, ntemp
           if ( rbfb(i) .gt. ( gft(i) - gff ) ) then
               rbfb(i) = rbfb(i) + log( 1.0d0 + exp( gft(i) - gff - rbfb(i) ) ) 
           else
               rbfb(i) = gft(i)- gff + log( 1.0d0 + exp( rbfb(i) - gft(i) + gff ) )
           endif
        enddo
     endif
  endif

  f_its = exp(gff-gftt)/beta0

  return

end subroutine its_force

!*******************************************************************************
!
! Subroutine:  its_prntrst
!
! Description: 
!              
!*******************************************************************************
subroutine its_prntrst

  use mdin_its_dat_mod
  use file_io_dat_mod
  use file_io_mod

  implicit none

! Formal arguments:


! Local variables:

  integer                       :: i

! begin of code
  call amopen(fbrestrt, fbrestrt_name, owrite, 'F', 'W')
  call amopen(normlrestrt, normlrestrt_name, owrite, 'F', 'W')
  do i = 1, ntemp-1
     write(fbrestrt, *) fb(i)
     write(normlrestrt, *) norml(i)
  enddo
  write(fbrestrt, *) fb(ntemp)
  close(fbrestrt)
  close(normlrestrt)

  return

end subroutine its_prntrst

!*******************************************************************************
!
! Subroutine:  its_prnttrj
!
! Description: 
!              
!*******************************************************************************
subroutine its_prnttrj

  use mdin_its_dat_mod
  use file_io_dat_mod
  use file_io_mod

  implicit none

! Formal arguments:


! Local variables:

  integer                       :: i

! begin of code
  call amopen(fbtrj, fbtrj_name, owrite, 'F', 'A')
  call amopen(normltrj, normltrj_name, owrite, 'F', 'A')
  do i = 1, ntemp-1
     write(fbtrj,*) mcycle,i,fb(i)
     write(normltrj,*) mcycle,i,norml(i)
  enddo
  write(fbtrj,*) mcycle,ntemp,fb(ntemp)
  close(fbtrj)
  close(normltrj)

  return

end subroutine its_prnttrj

end module its_force_mod
