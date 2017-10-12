#include "copyright.i"

!*******************************************************************************
!
! Module: remd_exchg_mod
!
! Description: 
!
! This module holds all of the subroutines that conduct exchanges.
! 
! This module needed to be introduced external to the REMD module because many
! of the more advanced exchange subroutines introduced dependencies for remd_mod
! that resulted in cyclic dependencies that prevented compiling. It's probably
! good to have them separate, anyway.
!
! IMPORTANT: If you add new exchange types (in the form of subroutines) here,
! you MUST make sure that, if you must call gb/pme_force, you do so with the
! temporary force array (frc_temp), and in the event of a successful exchange,
! you *must* replace frc with frc_temp so our forces are right for our first
! step after this exchange. Note if the exchange attempt fails, DO NOT perform
! this update. Also, any stand-alone exchange type should write data to the
! remlog if and *only* if it is not called from multid_exchange. Return before
! any print-outs are done if remd_method .eq. -1
!
!*******************************************************************************

module remd_exchg_mod

  use parallel_dat_mod
  use random_mod, only : random_state, amrset_gen, amrand_gen
  use remd_mod

  implicit none

#ifdef MPI
! Constants

  double precision, parameter :: ONEKB = 503.01d0 ! 1/Kb in internal units
  double precision, parameter :: TINY = 1.d-6

! To get pH-REMD to agree with sander

  logical, save :: jumpright = .true.

! Random number generator

  type(random_state), save, private :: remd_rand_gen

contains

!*******************************************************************************
!
! Subroutine: multid_exchange
!
! Description:
!
! This subroutine is responsible for coordinating the setting up of the REMD
! comms to facilitate multi-dimensional replica exchange attempts. You'll need
! to edit the select case() sections to add new exchange types
!
!*******************************************************************************

subroutine multid_exchange(atm_cnt, crd, vel, frc, remd_ptot, &
                           actual_temperature, my_atm_lst, gbl_img_atm_map, &
                           gbl_atm_img_map, mdloop)

  use file_io_dat_mod, only : remd_dimension_name, remd_file, max_fn_len
  use pmemd_lib_mod, only   : strip

  implicit none

! Passed variables
  
  integer, intent(in) :: atm_cnt
  integer, intent(in) :: mdloop
  integer             :: my_atm_lst(*)
  integer             :: gbl_img_atm_map(*)
  integer             :: gbl_atm_img_map(*)

  double precision, intent(in) :: remd_ptot
  double precision, intent(in) :: actual_temperature
  double precision             :: crd(3,atm_cnt)
  double precision             :: frc(3,atm_cnt)
  double precision             :: vel(3,atm_cnt)

! Local variables

  double precision :: l_fe
  double precision :: r_fe

  integer   :: my_dim
  integer   :: i, j, k
  integer   :: rep
  integer   :: gid ! group id

  logical   :: group_master

  character :: torf_char

  character(len=max_fn_len) :: filename
  character(len=max_fn_len) :: buf

! Variable descriptions
!
!  my_dim       : which dimension we choose to exchange between. Alternate 
!                 b/w different dimensions
!  i            : counter/iterator
!  buf          : holder string for the remlog filename extension
!  group_master : master of the group_master_comm
!  filename     : the full file name of the remlog we need to write to
!  torf_char    : character for if we succeeded (T) or not (F) in H-REMD output

  ! Determine my_dim, starting with the lowest dimension. We want my_dim to
  ! cycle between 1 - remd_dimension, starting at 1

  my_dim = mod(mdloop - 1, remd_dimension) + 1

  ! Nobody is group master yet, this will be set later

  group_master = .false.
  
  ! Set up REMD communicators, and free them if they're still formed

  if (master) then

    if (remd_comm .ne. mpi_comm_null) &
      call mpi_comm_free(remd_comm, err_code_mpi)
    remd_comm = mpi_comm_null

    call mpi_comm_split(pmemd_master_comm, group_num(my_dim), &
                        master_rank, remd_comm, err_code_mpi)
    call mpi_comm_size(remd_comm, remd_size, err_code_mpi)
    call mpi_comm_rank(remd_comm, remd_rank, err_code_mpi)

    remd_master = remd_rank .eq. 0

    if (group_master_comm .ne. mpi_comm_null) &
      call mpi_comm_free(group_master_comm, err_code_mpi)
    group_master_comm = mpi_comm_null

    if (remd_master) then
      call mpi_comm_split(pmemd_master_comm, 0, master_rank, &
                          group_master_comm, err_code_mpi)
      call mpi_comm_rank(group_master_comm, group_master_rank, err_code_mpi)
      call mpi_comm_size(group_master_comm, num_rem_grp, err_code_mpi)
      group_master = group_master_rank .eq. 0
    else
      ! If we are not remd_master then call split with MPI_UNDEFINEDs.
      call mpi_comm_split(pmemd_master_comm, MPI_UNDEFINED, MPI_UNDEFINED, &
                          group_master_comm, err_code_mpi)
    end if

  end if ! master

  select case (remd_types(my_dim))

    case(1)
      call temperature_exchange(atm_cnt, vel, remd_ptot, my_dim, remd_size, &
                                actual_temperature, .true., mdloop)
    case(3)
      call hamiltonian_exchange(atm_cnt, my_atm_lst, crd, vel, frc, remd_ptot, &
                                gbl_img_atm_map, gbl_atm_img_map, my_dim, &
                                remd_size, .true., mdloop)

  end select

  ! Now collect all of the multid_print_data's for all of the various 
  ! remd_masters to the group_master

  if (remd_master) then
    call mpi_allgather(multid_print_data, SIZE_REMLOG_DATA*numgroups, &
                       mpi_double_precision, multid_print_data_buf, &
                       SIZE_REMLOG_DATA*numgroups, mpi_double_precision, &
                       group_master_comm, err_code_mpi)

   ! Update the free energy data structures if this dimension is H-REMD

    if (remd_types(my_dim) .eq. 3) then

      do i = 1, num_rem_grp

        gid = int(multid_print_data_buf(1,i)%group_num)

        do j = 1, int(multid_print_data_buf(1,i)%num_rep)

          rep = int(multid_print_data_buf(j,i)%repnum)
          num_left_exchg(my_dim,gid,rep) = num_left_exchg(my_dim,gid,rep) &
                 + int(multid_print_data_buf(j,i)%left_exchg)
          num_right_exchg(my_dim,gid,rep) = num_right_exchg(my_dim,gid,rep) &
                 + int(multid_print_data_buf(j,i)%right_exchg)
          total_left_fe(my_dim,gid,rep) = total_left_fe(my_dim,gid,rep) &
                 + multid_print_data_buf(j,i)%left_fe
          total_right_fe(my_dim,gid,rep) = total_right_fe(my_dim,gid,rep) &
                 + multid_print_data_buf(j,i)%right_fe

        end do ! j = 1, num_rep

      end do ! i = 1, num_rem_grp

    end if ! remd_types(my_dim) .eq. 3

  end if

  ! Now that we have all of the data on group_master, have *that* process
  ! open up the appropriate remlog and write out the data to that file. We
  ! call the intrinsic "open" so we can append to this old file (it should
  ! have already been opened via amopen in the setup routine). This is the
  ! approach taken for the hard-flushing of mdout/mdinfo as well (so this has
  ! the side-effect of flushing this rem.log every time we reach this point

  if (group_master) then

    write(buf, '(i5)') my_dim
    call strip(buf)
    filename = trim(remlog_name) // '.' // trim(buf)

    open(unit=remlog, file=filename, status='OLD', position='APPEND')
    
    ! Modify your particular REMD's exchange method's printout here:
  
    select case(remd_types(my_dim))
  
      case(1) ! TEMPERATURE
        do i = 1, num_rem_grp

          gid = int(multid_print_data_buf(1,i)%group_num)

          write(remlog, '(2(a,i8))') '# exchange ', mdloop, ' REMD group ', gid

          do j = 1, int(multid_print_data_buf(1,i)%num_rep)
            write(remlog, '(i2,6f10.2,i8)') &
               j, &
               multid_print_data_buf(j,i)%scaling,       &
               multid_print_data_buf(j,i)%real_temp,     &
               multid_print_data_buf(j,i)%pot_ene_tot,   &
               multid_print_data_buf(j,i)%temp0,         &
               multid_print_data_buf(j,i)%new_temp0,     &
               multid_print_data_buf(j,i)%success_ratio, &
               int(multid_print_data_buf(j,i)%struct_num)
          end do ! j = 1, %num_rep

        end do ! i = 1, num_rem_grp

      case(3) ! HAMILTONIAN
        ! We need to update our running tally of number of left/right exchanges
        ! as well as our total left and right free energies
        do i = 1, num_rem_grp

          gid = int(multid_print_data_buf(1,i)%group_num)

          write(remlog, '(2(a,i8))') '# exchange ', mdloop, ' REMD group ', gid

          do j = 1, int(multid_print_data_buf(1,i)%num_rep)
            
            rep = int(multid_print_data_buf(j,i)%repnum)
            
            if (num_right_exchg(my_dim, i, rep) .eq. 0) then
              r_fe = 0.d0
            else
              r_fe = multid_print_data_buf(j,i)%temp0 / ONEKB * &
                     log(total_right_fe(my_dim, gid, rep) / &
                         num_right_exchg(my_dim, gid, rep))
            end if

            if (num_left_exchg(my_dim, i, rep) .eq. 0) then
              l_fe = 0.d0
            else
              l_fe = multid_print_data_buf(j,i)%temp0 / ONEKB * &
                     log(total_left_fe(my_dim, gid, rep) / &
                         num_left_exchg(my_dim, gid, rep))
            end if

            if (multid_print_data_buf(j,i)%success .eq. 1.d0) then
              torf_char = 'T'
            else
              torf_char = 'F'
            end if

            write(remlog, '(2i6,5f10.2,4x,a,2x,f10.2)')      &
               int(multid_print_data_buf(j,i)%repnum),       &
               int(multid_print_data_buf(j,i)%neighbor_rep), &
               multid_print_data_buf(j,i)%temp0,             &
               multid_print_data_buf(j,i)%pot_ene_tot,       &
               multid_print_data_buf(j,i)%nei_pot_ene,       &
               l_fe, &
               r_fe, &
               torf_char, &
               multid_print_data_buf(j,i)%success_ratio

          end do ! j = 1, %num_rep

        end do ! i = 1, num_rem_grp

    end select

    close(remlog)

  end if ! group_master

  return

end subroutine multid_exchange

!*******************************************************************************
!
! Subroutine: temperature_exchange
!
! Description: Performs the temperature replica exchange attempt. We don't need
!              any extra energies here.
!
!*******************************************************************************

subroutine temperature_exchange(atm_cnt, vel, my_pot_ene_tot, t_dim, &
                                num_replicas, actual_temperature,    &
                                print_exch_data, mdloop)

  use mdin_ctrl_dat_mod, only : temp0
  use pmemd_lib_mod,     only : mexit

  implicit none
 
! Passed variables

  integer, intent(in)             :: atm_cnt
  integer, intent(in)             :: num_replicas
  integer, intent(in)             :: mdloop

  double precision, intent(in)    :: my_pot_ene_tot
  double precision, intent(inout) :: vel(3, atm_cnt)
  double precision, intent(in)    :: actual_temperature
  
  integer, intent(in)             :: t_dim ! which dimension T-exchanges are
  logical, intent(in)             :: print_exch_data

! Local variables

  type :: exchange_data ! facilitates gather for remlog
    sequence
    double precision :: scaling
    double precision :: real_temp
    double precision :: pot_ene_tot
    double precision :: temp0
    double precision :: new_temp0
    double precision :: struct_num
  end type exchange_data

  integer, parameter :: SIZE_EXCHANGE_DATA = 6 ! for mpi_gather

  type (exchange_data) :: exch_buffer(num_replicas) ! gather buffer

  type (exchange_data) :: my_exch_data ! stores data for this replica
  type (exchange_data) :: exch_data_tbl(num_replicas) ! exch_data for all reps

  double precision  :: delta
  double precision  :: metrop
  double precision  :: random_value
  double precision  :: success_ratio

  integer           :: neighbor_rank
  integer           :: i
  integer           :: success_array(num_replicas)
  integer           :: success_buf(num_replicas)

  ! For exchanging replica positions in ladders upon successful exchanges
  integer           :: group_num_buffer(remd_dimension)
  integer           :: replica_indexes_buffer(remd_dimension)

  logical           :: success
  logical           :: i_do_exchg

  integer, dimension(mpi_status_size) :: stat_array ! needed for mpi_recv

! Explanation of local variables:
!
! delta: The calculated DELTA value used as part of the detailed balance calc.
!
! metrop: The actual value compared to a random_value to decide success
!
! success: Did this exchange attempt succeed?
!
! success_array: buffer for keeping track of which temperatures exchanged
!                successfully. Used in a mpi_reduce at the end
!
! success_buf: buffer for mpi_reduce call on success_array
!
! neighbor_rank: the remd_rank of the replica with whom we need to be doing ALL
!                of our communication. Trying to do everything with the
!                'partners' array is a bit too confusing.

! Set the variables that we know at the beginning

  my_exch_data%real_temp   = actual_temperature
  my_exch_data%temp0       = temp0
  my_exch_data%new_temp0   = temp0 ! changed later if exch. succeeds
  my_exch_data%struct_num  = -1    ! not implemented yet
  my_exch_data%pot_ene_tot = my_pot_ene_tot
  my_exch_data%scaling     = -1.d0
  success_array(:)         = 0

  if (master) then

    call set_partners(t_dim, num_replicas)

    ! Call amrand_gen here to generate a random number. Sander does this
    ! to pick out structures from a reservoir which pmemd doesn't implement
    ! (yet). To keep the random #s synched, though, and to ensure proper
    ! regression testing, this must be done to keep the random #s synched
    ! b/w sander and pmemd

    call amrand_gen(remd_rand_gen, random_value)

#ifdef VERBOSE_REMD
    write(mdout, '(26("="),a,26("="))') 'REMD EXCHANGE CALCULATION'
    write(mdout, '(a,i10,a,i1)') 'Exch= ', mdloop, ' RREMD= ', 0
#endif

    ! We alternate exchanges: temp up followed by temp down
    if (even_exchanges(t_dim)) then
      i_do_exchg = even_replica(t_dim)
    else
      i_do_exchg = .not. even_replica(t_dim)
    end if

    even_exchanges(t_dim) = .not. even_exchanges(t_dim)

    ! Find out what rank our neighbor is (-1 for indexing from 0), and
    ! If we are controlling the exchange, then our partner will be our 
    ! 2nd partner (above). Otherwise it will be our 1st partner (below)
    
    if (i_do_exchg) then
      neighbor_rank = partners(2) - 1
    else
      neighbor_rank = partners(1) - 1
    end if

    ! Collect potential energies/temperatures

    call mpi_allgather(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                exch_data_tbl, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                remd_comm, err_code_mpi)

#ifdef VERBOSE_REMD
    write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') 'Replica          Temp= ', &
      my_exch_data%temp0, ' Indx= ', replica_indexes(t_dim), ' Rep#= ', &
      remd_rank+1, ' EPot= ', my_pot_ene_tot
#endif

    ! Calculate the exchange probability. Note that if delta is negative, then
    ! metrop will be > 1, and it will always succeed (so we don't need an extra
    ! conditional)
    
    if (i_do_exchg) then
#ifdef VERBOSE_REMD
      write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') 'Partner          Temp= ', &
        exch_data_tbl(neighbor_rank+1)%temp0, ' Indx= ', neighbor_rank+1, &
        ' Rep#= ', neighbor_rank + 1, ' EPot= ', &
        exch_data_tbl(neighbor_rank+1)%pot_ene_tot
#endif

      delta = (my_exch_data%pot_ene_tot - &
               exch_data_tbl(neighbor_rank+1)%pot_ene_tot) * &
              (my_exch_data%temp0 - &
               exch_data_tbl(neighbor_rank+1)%temp0) * ONEKB /  &
              (my_exch_data%temp0 * exch_data_tbl(neighbor_rank+1)%temp0)

      metrop = exp(-delta)

      call amrand_gen(remd_rand_gen, random_value)

      success = random_value .lt. metrop

      ! Let my neighbor know about our success

      call mpi_send(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, remd_comm, err_code_mpi)

      if (success) then
        success_array(replica_indexes(t_dim)) = 1
        my_exch_data%new_temp0 = exch_data_tbl(neighbor_rank+1)%temp0
        my_exch_data%scaling = sqrt(my_exch_data%new_temp0 / my_exch_data%temp0)
      end if

#ifdef VERBOSE_REMD
      write(mdout,'(a8,E16.6,a8,E16.6,a12,f10.2)') &
               "Metrop= ",metrop," delta= ",delta," o_scaling= ", &
               1 / my_exch_data%scaling
#endif

    else

#ifdef VERBOSE_REMD
      ! Even if we don't control exchange, still print REMD info to mdout

      write(mdout, '(a,f6.2,2(a7,i2),a7,f10.2)') 'Partner          Temp= ', &
        exch_data_tbl(neighbor_rank+1)%temp0, ' Indx= ', neighbor_rank+1, &
        ' Rep#= ', neighbor_rank + 1, ' EPot= ', &
        exch_data_tbl(neighbor_rank+1)%pot_ene_tot
      write(mdout, '(a)') 'Not controlling exchange.'
#endif

      call amrand_gen(remd_rand_gen, random_value) ! to stay synchronized

      ! Get the message from the exchanging replica about our success
      call mpi_recv(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, remd_comm, stat_array, err_code_mpi)

      ! We scale velocities only if we succeed
      if (success) then
        my_exch_data%new_temp0 = exch_data_tbl(neighbor_rank+1)%temp0
        my_exch_data%scaling = sqrt(my_exch_data%new_temp0 / my_exch_data%temp0)
      end if

    end if

#ifdef VERBOSE_REMD
    write(mdout,'(a8,E16.6,a12,f10.2,a10,L1)') &
      'Rand=   ', random_value, ' MyScaling= ', my_exch_data%scaling, &
      ' Success= ', success
    
    write(mdout,'(24("="),a,24("="))') "END REMD EXCHANGE CALCULATION"
#endif

  end if ! master

  ! We'll broadcast all of the my_exch_data here once. From that, we can tell
  ! if the exchange succeeded based on whether or not the temperatures have
  ! changed. Then make sure we update temp0 so the simulation reflects that.

  call mpi_bcast(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                 0, pmemd_comm, err_code_mpi)

  temp0 = my_exch_data%new_temp0

  if (.not. master) &
    success = abs(my_exch_data%temp0 - my_exch_data%new_temp0) .gt. TINY

#ifdef VERBOSE_REMD
  if (master) &
    write(mdout, '(2(a,f6.2))') &
      'REMD: checking to see if bath T has changed: ', &
      my_exch_data%temp0, '->', my_exch_data%new_temp0
#endif

  if (success) call rescale_velocities(atm_cnt, vel, my_exch_data%scaling)
    
  ! We want to keep track of our exchange successes, so gather them
  ! from all of the masters here, and increment them in exchange_successes
  ! We do this whether or not we print so we can keep track of exchange
  ! success ratios. For the index and exchange data, only gather those to
  ! the remd_master if we actually plan to print. This should be much
  ! more efficient for high exchange attempts.

  ! Then recollect the temperatures as some may have changed, and reset our
  ! partners.

  if (master) then

    if (print_exch_data) then
      call mpi_gather(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                      exch_buffer, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                      0, remd_comm, err_code_mpi)
    end if

    call mpi_reduce(success_array, success_buf, num_replicas, mpi_integer, &
                    mpi_sum, 0, remd_comm, err_code_mpi)

    if (remd_master) &
      exchange_successes(t_dim,:) = exchange_successes(t_dim,:) + success_buf(:)

    ! Increment exchange counter, and swap our replica ranks and numbers with 
    ! our neighbor if our attempt succeeded, since we effectively swapped places
    ! in this array. That's because we swap temperatures and not structures due
    ! to performance considerations.

    if (success) then
      call mpi_sendrecv(group_num, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, &
                        group_num_buffer, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, remd_comm, &
                        stat_array, err_code_mpi)
      call mpi_sendrecv(replica_indexes, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, &
                        replica_indexes_buffer, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, remd_comm, &
                        stat_array, err_code_mpi)
      group_num = group_num_buffer
      replica_indexes = replica_indexes_buffer
    end if
  end if
  
  ! If we're doing NMROPT, then we need to make sure we re-read in the TEMP0,
  ! since this subroutine overwrites TEMP0 to the original since you can modify
  ! temperature with NMROPT settings.

  remd_modwt = .true.

  ! Write the data to the remlog

  ! remd_method == -1 means we're doing multi-D REMD.
  if (remd_method .ne. -1) then
    if (print_exch_data .and. remd_master) then
      write(remlog, '(a,i8)') '# exchange ', mdloop
      do i = 1, num_replicas
        success_ratio = dble(exchange_successes(t_dim, index_list(i))) / &
                        dble(mdloop) * dble(remd_dimension) * 2
        write(remlog, '(i2,6f10.2,i8)') &
              i, &
              exch_buffer(i)%scaling, &
              exch_buffer(i)%real_temp, &
              exch_buffer(i)%pot_ene_tot, &
              exch_buffer(i)%temp0, &
              exch_buffer(i)%new_temp0, &
              success_ratio, &
              int(exch_buffer(i)%struct_num)
      end do
    end if
  else if (remd_master) then
    do i = 1, num_replicas
      success_ratio = dble(exchange_successes(t_dim, index_list(i))) / &
                      dble(mdloop) * dble(remd_dimension) * 2
      multid_print_data(i)%scaling       = exch_buffer(i)%scaling
      multid_print_data(i)%real_temp     = exch_buffer(i)%real_temp
      multid_print_data(i)%pot_ene_tot   = exch_buffer(i)%pot_ene_tot
      multid_print_data(i)%temp0         = exch_buffer(i)%temp0
      multid_print_data(i)%new_temp0     = exch_buffer(i)%new_temp0
      multid_print_data(i)%group_num     = group_num(t_dim)
      multid_print_data(i)%success_ratio = success_ratio
      multid_print_data(i)%num_rep       = num_replicas
    end do
  end if ! remd_method .ne. -1

  return

end subroutine temperature_exchange

!*******************************************************************************
!
! Subroutine: hamiltonian_exchange
!
! Description: Performs the Hamiltonian replica exchange attempt. Additional
!              energies need to be computed here.
!
!*******************************************************************************

subroutine hamiltonian_exchange(atm_cnt, my_atm_lst, crd, vel, frc, &
                                my_pot_ene_tot, gbl_img_atm_map, &
                                gbl_atm_img_map, h_dim, &
                                num_replicas, print_exch_data, mdloop)

  use gb_force_mod
  use mdin_ctrl_dat_mod
  use nmr_calls_mod, only : nmrdcp, skip_print
  use pme_force_mod

  implicit none

! Passed variables

  integer          :: atm_cnt
  integer          :: num_replicas
  integer          :: mdloop
  integer          :: my_atm_lst(*)
  integer          :: gbl_img_atm_map(*)
  integer          :: gbl_atm_img_map(*)
  logical          :: print_exch_data
  double precision :: crd(3,atm_cnt)
  double precision :: frc(3,atm_cnt)
  double precision :: vel(3,atm_cnt)
  integer          :: h_dim
  double precision :: my_pot_ene_tot

! Local variables

  type :: ene_temp ! Useful type to limit communication
    sequence
    double precision :: neighbor_rep ! neighbor replica # (for remlog printing)
    double precision :: energy_1     ! energy with MY coordinates
    double precision :: energy_2     ! energy with THEIR coordinates
    double precision :: temperature  ! my temperature
    double precision :: repnum       ! replica number
    double precision :: left_fe      ! My left free energy for this step
    double precision :: right_fe     ! My right free energy for this step
    double precision :: right_exchg  ! 1 if we exchanged to the right, 0 left
    double precision :: left_exchg   ! 1 if we exchanged to the left, 0 right
  end type ene_temp
  integer, parameter :: SIZE_ENE_TEMP = 9

  double precision :: my_pot_ene_tot_2   ! MY pot ene with THEIR coordinates
  double precision :: delta 
  double precision :: metrop
  double precision :: random_value
  double precision :: success_ratio
  double precision :: virial(3)          ! For pme_force_call
  double precision :: ekcmt(3)           ! For pme_force_call
  double precision :: pme_err_est        ! For pme_force_call
  double precision :: vel_scale_coff     ! Velocity scaling factor
  double precision :: l_fe, r_fe         ! Holders for FE calcs

  logical          :: i_do_exchg
  logical          :: rank_success(num_replicas) ! log of all replica's success
  logical          :: success
  integer          :: success_array(num_replicas) ! to keep track of successes
  integer          :: success_buf(num_replicas) ! reduce buffer for successes
  character        :: success_char ! to print T or F if we succeeded in remlog
  integer          :: neighbor_rank
  integer          :: istat(mpi_status_size)
  integer          :: ier, i
  integer          :: rep ! replica counter

  type(gb_pot_ene_rec)  :: my_new_gb_ene  ! gb energy record for THEIR coords
  type(pme_pot_ene_rec) :: my_new_pme_ene ! pme energy record for THEIR coords
  type(ene_temp)   :: my_ene_temp
  type(ene_temp)   :: neighbor_ene_temp
  type(ene_temp)   :: ene_temp_buffer(num_replicas) ! to collect/print REMD data


  ! Get coordinates and forces from GPU
#ifdef CUDA
  call gpu_download_crd(crd)
  call gpu_download_frc(frc)
#endif


  ! Set up some stuff for pme_force call:

  virial(:) = 0.d0
  ekcmt(:) = 0.d0

  ! First determine if we're controlling the exchange or not

  if (master) then

    ! Figure out our partners

    call set_partners(h_dim, num_replicas)

    ! Initialize some data

    success_array(:) = 0
    my_ene_temp%energy_1 = my_pot_ene_tot
    my_ene_temp%temperature = temp0
    my_ene_temp%repnum = replica_indexes(h_dim)
    my_ene_temp%left_fe = 0.d0
    my_ene_temp%right_fe = 0.d0
    my_ene_temp%right_exchg = 0.d0
    my_ene_temp%left_exchg = 0.d0

#ifdef VERBOSE_REMD
    write(mdout, '(24("="),a,24("="))') ' H-REMD EXCHANGE CALCULATION '
    write(mdout, '(a,i10,a,i1)') 'Exch= ', mdloop, ' RREMD= ', 0
#endif
    
    ! Alternate exchanges

    if (even_exchanges(h_dim)) then
      i_do_exchg = even_replica(h_dim)
    else
      i_do_exchg = .not. even_replica(h_dim)
    end if

    even_exchanges(h_dim) = .not. even_exchanges(h_dim)
    
    ! Set partner rank 
    if (i_do_exchg) then
      my_ene_temp%neighbor_rep = partners(2)
      neighbor_rank = partners(2) - 1
    else
      my_ene_temp%neighbor_rep = partners(1)
      neighbor_rank = partners(1) - 1
    end if

    ! Exchange coordinates with your neighbor. crd_temp must be a backup since
    ! the reciprocal space sum needs the global array to be updated with the new
    ! coordinates (and crd is a reference to the global array)
    crd_temp(1,:) = crd(1,:)
    crd_temp(2,:) = crd(2,:)
    crd_temp(3,:) = crd(3,:)
    call mpi_sendrecv_replace(crd, atm_cnt * 3, mpi_double_precision, &
               neighbor_rank, remd_tag, &
               neighbor_rank, remd_tag, &
               remd_comm, istat, err_code_mpi)

  end if ! master

  ! Now broadcast the temporary coordinates to the rest of pmemd_comm. I don't
  ! know if this is necessary or not (i.e. is it done in the load balancing?)

  call mpi_bcast(crd, atm_cnt * 3, mpi_double_precision, &
                 0, pmemd_comm, err_code_mpi)
  call mpi_bcast(crd_temp, atm_cnt * 3, mpi_double_precision, &
                 0, pmemd_comm, err_code_mpi)

#ifdef CUDA
  ! Upload the temporary coordinates
  call gpu_upload_crd(crd)
#endif

  ! We do not want to print out any DUMPAVE values for this force call

  if (nmropt .ne. 0) skip_print = .true.

  ! Now it's time to get the energy of the other structure:
  if (using_gb_potential) then
    call gb_force(atm_cnt, crd, frc_temp, my_new_gb_ene, irespa, .true.)

    my_ene_temp%energy_2 = my_new_gb_ene%total

  else if (using_pme_potential) then
#ifdef CUDA
    ! CUDA pairlist building is not controlled by new_list, it is handled
    ! internally. Since we are changing coordinates, we must build a new
    ! pairlist
    call gpu_force_new_neighborlist()
#endif
    ! First .true. is new_list, next is need_pot_enes, .false. is need_virials
! Modified by Lijiang Yang
    call pme_force(atm_cnt, crd, frc_temp, frc_temp, frc_temp, frc_temp, gbl_img_atm_map, &
!    call pme_force(atm_cnt, crd, frc_temp, gbl_img_atm_map, &
!                   gbl_atm_img_map, my_atm_lst, .true., .true., .false., &
                   gbl_atm_img_map, my_atm_lst, .true., .true., .true., .false., &
!-----------------------------------------------------------
                   my_new_pme_ene, virial, ekcmt, pme_err_est)
    
    my_ene_temp%energy_2 = my_new_pme_ene%total
  end if

  ! Since we made an 'extra' force call just now, decrement the nmropt counter

  if (nmropt .ne. 0) call nmrdcp

  ! Now trade that potential energy with your partner
  if (master) then

    call mpi_sendrecv(my_ene_temp, SIZE_ENE_TEMP, mpi_double_precision, &
                      neighbor_rank, remd_tag, &
                      neighbor_ene_temp, SIZE_ENE_TEMP, mpi_double_precision, &
                      neighbor_rank, remd_tag, &
                      remd_comm, istat, err_code_mpi)

  ! Now it's time to calculate the exchange criteria. This detailed balance was
  ! derived under the assumption that *if* we have different T's, then they are
  ! NOT exchanged. Swapping temperatures requires a slightly different equation
  ! and then requires you to update temp0 at the end with an mpi_sendrecv call.

    if (i_do_exchg) then

      delta = -ONEKB / temp0 * (my_ene_temp%energy_2 - my_ene_temp%energy_1) - &
              ONEKB / neighbor_ene_temp%temperature * &
              (neighbor_ene_temp%energy_2 - neighbor_ene_temp%energy_1)

      metrop = exp(delta)

      call amrand_gen(remd_rand_gen, random_value)

      success = random_value .lt. metrop

      ! Let my neighbor know about our success

      call mpi_send(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, remd_comm, err_code_mpi)
      if (success) then
#ifdef CUDA
        call gpu_download_vel(vel)
#endif
        success_array(replica_indexes(h_dim)) = 1
        ! swapping velocities after exchange succeed
        ! and scaling velocities   (DSD 09/12)
        call mpi_sendrecv_replace(vel, atm_cnt * 3, mpi_double_precision, &
                   neighbor_rank, remd_tag, &
                   neighbor_rank, remd_tag, &
                   remd_comm, istat, err_code_mpi)
        vel_scale_coff = sqrt(temp0 / neighbor_ene_temp%temperature)
      end if

#ifdef VERBOSE_REMD
      write(mdout, '(a,f16.6)') 'My Eptot_1:       ', my_ene_temp%energy_1
      write(mdout, '(a,f16.6)') 'My Eptot_2:       ', my_ene_temp%energy_2
      write(mdout, '(a,f16.6)') 'Neighbor Eptot_1: ', neighbor_ene_temp%energy_1
      write(mdout, '(a,f16.6)') 'Neighbor Eptot_2: ', neighbor_ene_temp%energy_2
      write(mdout, '(3(a,f16.6))') 'Delta= ', delta, ' Metrop= ', metrop, &
                                   ' Random #= ', random_value
      if (success) then
        write(mdout, '(a)') 'Exchange Succeeded!'
      else
        write(mdout, '(a)') 'Exchange Failed!'
      end if
#endif

    else
      
      call amrand_gen(remd_rand_gen, random_value) ! keep in sync

      call mpi_recv(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, remd_comm, istat, err_code_mpi)
      if (success) then
#ifdef CUDA
        call gpu_download_vel(vel)
#endif
        call mpi_sendrecv_replace(vel, atm_cnt * 3, mpi_double_precision, &
                      neighbor_rank, remd_tag, &
                      neighbor_rank, remd_tag, &
                      remd_comm, istat, err_code_mpi)
        vel_scale_coff = sqrt(temp0 / neighbor_ene_temp%temperature)
      end if

#ifdef VERBOSE_REMD
      write(mdout, '(a,f16.6)') 'My Eptot_1:       ', my_ene_temp%energy_1
      write(mdout, '(a,f16.6)') 'My Eptot_2:       ', my_ene_temp%energy_2
      write(mdout, '(a,f16.6)') 'Neighbor Eptot_1: ', neighbor_ene_temp%energy_1
      write(mdout, '(a,f16.6)') 'Neighbor Eptot_2: ', neighbor_ene_temp%energy_2
      write(mdout, '(3(a))') '| Not controlling exchange.'
      if (success) then
        write(mdout, '(a)') 'Exchange Succeeded!'
      else
        write(mdout, '(a)') 'Exchange Failed!'
      end if
#endif

    end if ! i_do_exchg

  end if ! master
  
  ! Now broadcast our success for all to hear! This is the ONLY part that is not
  ! for masters only

  call mpi_bcast(success, 1, mpi_logical, 0, pmemd_comm, err_code_mpi)

#ifdef CUDA
  ! If the exchange succeeded, then the coordinates and forces on the gpu
  ! are already correct and we don't need to do anything. Otherwise, we
  ! need to send back the originals.
  if (.not. success) then
    ! Upload the orginal coordinates again and force a new neighborlist if we're
    ! running PME
    if (using_pme_potential) &
      call gpu_force_new_neighborlist()
    call gpu_upload_crd(crd_temp)
    call gpu_upload_frc(frc)
  else
    call mpi_bcast(vel, atm_cnt * 3, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(vel_scale_coff, 1, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call gpu_upload_vel(vel)
    call rescale_velocities(atm_cnt, vel, vel_scale_coff)
  end if
#else
  ! Update our forces if we succeeded. Note that the exchange probability
  ! equation used currently satisfies detailed balance only for NOT exchanging
  ! temperatures. We also update our force array, since we won't have to run
  ! our force call on the first pass through our next MD loop in runmd
  if (success) then
    frc(1, :) = frc_temp(1, :)
    frc(2, :) = frc_temp(2, :)
    frc(3, :) = frc_temp(3, :)
    call mpi_bcast(vel, atm_cnt * 3, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call mpi_bcast(vel_scale_coff, 1, mpi_double_precision, 0, &
                   pmemd_comm, err_code_mpi)
    call rescale_velocities(atm_cnt, vel, vel_scale_coff)
  else
    ! Restore the coordinates
    crd(1, :) = crd_temp(1, :)
    crd(2, :) = crd_temp(2, :)
    crd(3, :) = crd_temp(3, :)
  end if
#endif /* CUDA */

  ! Calculate the free energies for jumping right and jumping left. This follows
  ! free energy perturbation (FEP). The equation for FEP is
  !
  !                          /        -(Eb - Ea)     \
  ! D G_(a->b) = -Kb * T ln | < exp (-------------) > |
  !                          \          Kb * T       /  A
  !
  ! Where the exponential is averaged. Thus, we have to keep a running total of
  ! the exponential Eb - Ea term (total_right_fe/total_left_fe), which can be
  ! combined into a total D G for this step via the above equation. These are
  ! stored in the ene_temp type to minimize communications. They can be 
  ! calculated regardless of whether we succeeded or not!
  ! For the layout of the rem.log file: left = up, right = down.
  ! Then calculate the 'averaged' free energy right after, keeping in mind that
  ! right/left alternate. Therefore, we introduce a counter for the number of
  ! left-moves and number of right-moves to make sure we divide by the right
  ! number of samples.

  if (master) then
    if (i_do_exchg) then

      my_ene_temp%right_fe = &
                   exp((my_ene_temp%energy_1 - neighbor_ene_temp%energy_2) * &
                   ONEKB / my_ene_temp%temperature)
      my_ene_temp%right_exchg = 1.d0
      my_ene_temp%left_fe = 0.d0
      
    else

      my_ene_temp%left_fe = &
                   exp((my_ene_temp%energy_1 - neighbor_ene_temp%energy_2) * &
                   ONEKB / my_ene_temp%temperature)
      my_ene_temp%left_exchg = 1.d0
      my_ene_temp%right_fe = 0.d0

    end if

    ! Collect data and write it out to the rem.log

    call mpi_reduce(success_array, success_buf, num_replicas, mpi_integer, &
                    mpi_sum, 0, remd_comm, err_code_mpi)
    
    call mpi_gather(success, 1, mpi_logical, rank_success, 1, mpi_logical, &
                    0, remd_comm, err_code_mpi)

    call mpi_gather(my_ene_temp, SIZE_ENE_TEMP, mpi_double_precision, &
                    ene_temp_buffer, SIZE_ENE_TEMP, mpi_double_precision, &
                    0, remd_comm, err_code_mpi)

    if (remd_master) then

      exchange_successes(h_dim,:) = exchange_successes(h_dim,:) + success_buf(:)

      if (remd_method .ne. -1 .and. print_exch_data) then
  
        do i = 1, num_replicas
          num_left_exchg(h_dim, group_num(h_dim), i)  = &
                num_left_exchg(h_dim, group_num(h_dim), i) &
              + int(ene_temp_buffer(i)%left_exchg)
  
          num_right_exchg(h_dim, group_num(h_dim), i) = &
                num_right_exchg(h_dim, group_num(h_dim), i) &
              + int(ene_temp_buffer(i)%right_exchg)
        end do

        write(remlog, '(a,i8)') '# exchange ', mdloop

        do rep = 1, num_replicas
  
          ! Did we succeed or not? We need to print a character
          if (rank_success(rep)) then
            success_char = 'T'
          else
            success_char = 'F'
          end if
          
          i = ene_temp_buffer(rep)%repnum

          total_left_fe(h_dim, group_num(h_dim), i) = &
             total_left_fe(h_dim, group_num(h_dim), i) + &
             ene_temp_buffer(rep)%left_fe

          total_right_fe(h_dim, group_num(h_dim), i) = &
             total_right_fe(h_dim, group_num(h_dim), i) + &
             ene_temp_buffer(rep)%right_fe
  
          ! Calculate the success ratio:
  
          success_ratio = dble(exchange_successes(h_dim, &
                          int(ene_temp_buffer(rep)%repnum))) / &
                          dble(mdloop) * dble(remd_dimension) * 2
  
          ! Print info to the remlog
          
          if (num_left_exchg(h_dim, group_num(h_dim), rep) > 0) then
            l_fe = ene_temp_buffer(rep)%temperature / ONEKB * &
                   log(total_left_fe(h_dim, group_num(h_dim), rep) / &
                   num_left_exchg(h_dim, group_num(h_dim), rep))
          else
            l_fe = 0.d0
          end if
          
          if (num_right_exchg(h_dim, group_num(h_dim), rep) > 0) then
            r_fe = ene_temp_buffer(rep)%temperature / ONEKB * &
                   log(total_right_fe(h_dim, group_num(h_dim), rep) / &
                   num_right_exchg(h_dim, group_num(h_dim), rep))
          else
            r_fe = 0.d0
          end if

          write(remlog, '(2i6,5f10.2,4x,a,2x,f10.2)')   &
                index_list(rep), &
                int(ene_temp_buffer(rep)%neighbor_rep), &
                ene_temp_buffer(rep)%temperature, &
                ene_temp_buffer(rep)%energy_1,    &
                ene_temp_buffer(rep)%energy_2,    &
                l_fe,         &
                r_fe,         &
                success_char, &
                success_ratio
  
        end do
      else
        do rep = 1, num_replicas

          success_ratio = dble(exchange_successes(h_dim, &
                          int(ene_temp_buffer(rep)%repnum))) /  &
                          dble(mdloop) * dble(remd_dimension) * 2
          if (rank_success(rep)) then
            multid_print_data(rep)%success = 1.d0
          else
            multid_print_data(rep)%success = 0.d0
          end if

          multid_print_data(rep)%neighbor_rep= ene_temp_buffer(rep)%neighbor_rep
          multid_print_data(rep)%temp0       = ene_temp_buffer(rep)%temperature
          multid_print_data(rep)%pot_ene_tot = ene_temp_buffer(rep)%energy_1
          multid_print_data(rep)%nei_pot_ene = ene_temp_buffer(rep)%energy_2
          multid_print_data(rep)%left_fe     = ene_temp_buffer(rep)%left_fe
          multid_print_data(rep)%right_fe    = ene_temp_buffer(rep)%right_fe
          multid_print_data(rep)%success_ratio= success_ratio
          multid_print_data(rep)%num_rep     = num_replicas
          multid_print_data(rep)%group_num   = group_num(h_dim)
          multid_print_data(rep)%repnum      = ene_temp_buffer(rep)%repnum
          multid_print_data(rep)%left_exchg  = ene_temp_buffer(rep)%left_exchg
          multid_print_data(rep)%right_exchg = ene_temp_buffer(rep)%right_exchg
        end do
      end if

    end if ! remd_master

#ifdef VERBOSE_REMD
    write(mdout, '(26("="),a,26("="))') 'END H-REMD CALCULATION'
#endif

  end if ! master

  remd_modwt = .true.

  return

end subroutine hamiltonian_exchange

!*******************************************************************************
!
! Subroutine: ph_remd_exchange
!
! Description: Performs the pH replica exchange attempt.
!
!*******************************************************************************

subroutine ph_remd_exchange(ph_dim, num_replicas, mdloop)

  use constantph_mod, only    : total_protonation
  use gbl_constants_mod, only : KB, LN_TO_LOG
  use mdin_ctrl_dat_mod, only : solvph
  use parallel_dat_mod

  implicit none

! Passed variables

  integer, intent(in) :: ph_dim
  integer, intent(in) :: num_replicas
  integer, intent(in) :: mdloop

! Local variables
  
  double precision :: delta                  ! DELTA value for MC transition
  double precision :: randval                ! Random #
  double precision :: all_ph(num_replicas)   ! pH table 
  double precision :: new_ph(num_replicas)   ! table of pHs after exchange
  double precision :: exchfrac(num_replicas) ! fraction of successes
  double precision :: o_ph                   ! pH of neighbor replica
  double precision :: success_ratio          ! success rate

  integer :: i    ! counter
  integer :: prot ! # of protons in this replica
  integer :: suc_arry(num_replicas)   ! array with success values
  integer :: suc_buf(num_replicas)    ! reduce buffer for the above
  integer :: prot_table(num_replicas) ! table of `prot' for each replica
  integer :: my_index                 ! my position in pH table
  integer :: o_index                  ! index of replica we are exchanging with
  integer :: o_repnum                 ! replica number we're exchanging with
  integer :: o_prot                   ! # of protons in other replica

  logical :: success                  ! did we succeed?
  logical :: i_do_exchg               ! Do I do the exchange?

  integer, dimension(mpi_status_size) :: stat_array ! needed for mpi_recv

  ! Initialize
  suc_arry(:) = 0
  suc_buf(:) = 0
  success = .false.
  delta = 0.d0
  randval = 0.d0

  if (master) then

    call set_partners(ph_dim, num_replicas)

#ifdef VERBOSE_REMD
    write(6,'(a)') '=============== REMD ==============='
#endif

    prot = total_protonation()

    call mpi_allgather(prot, 1, mpi_integer, prot_table, &
                       1, mpi_integer, remd_comm, err_code_mpi)
    
    call mpi_allgather(solvph, 1, mpi_double_precision, &
                       all_ph, 1, mpi_double_precision, &
                       remd_comm, err_code_mpi)

    my_index = replica_indexes(ph_dim)

    if (jumpright) then
      
      if (mod(my_index, 2) .eq. 0) then
        o_index = my_index + 1
        if (o_index > num_replicas) o_index = 1
        o_repnum = partners(2)
      else
        o_index = my_index - 1
        if (o_index < 1) o_index = num_replicas
        o_repnum = partners(1)
      end if

    else ! not jumpright

      if (mod(my_index, 2) .eq. 0) then
        o_index = my_index - 1
        if (o_index < 1) o_index = num_replicas
        o_repnum = partners(1)
      else
        o_index = my_index + 1
        if (o_index > num_replicas) o_index = 1
        o_repnum = partners(2)
      end if

    end if ! (jumpright)

    o_ph = all_ph(o_repnum)
    o_prot = prot_table(o_repnum)

#ifdef VERBOSE_REMD
    write(6,'(a,i3,a,f5.2,a,i3)') 'Found partner information: &
             &protonation = ', o_prot, ' pH = ', o_ph, ' repnum = ', o_repnum
#endif
    
    if (mod(my_index, 2) .eq. 0) then

      call amrand_gen(remd_rand_gen, randval)
      delta = LN_TO_LOG * (prot - o_prot) * (o_ph - solvph)
      success = randval < exp(-delta)

      ! Tell our neighbor if we succeeded

      call mpi_send(success, 1, mpi_logical, o_repnum-1, &
                    22, remd_comm, err_code_mpi)

      ! Increment our success

      if (success) then
        if (jumpright) then
          suc_arry(my_index) = 1
        else
          suc_arry(o_index) = 1
        end if
      end if

#ifdef VERBOSE_REMD
       write(6, '(2(a,i3))') 'Proton count: ', prot, ' --> ', o_prot
       write(6, '(2(a,f8.3))') 'pH transition:', solvph, ' --> ', o_ph
       if (success) then
         write(6, '(a,E16.8)') 'Success! delta = ', delta
       else
         write(6, '(a,E16.8)') 'Failure. delta = ', delta
       end if
#endif

    else ! I do not exchange
      
      call mpi_recv(success, 1, mpi_logical, o_repnum-1, &
                    22, remd_comm, stat_array, err_code_mpi)
#ifdef VERBOSE_REMD
         write(6, '(2(a,i3))') 'Proton count: ', prot, ' --> ', o_prot
         write(6, '(2(a,f8.3))') 'pH transition:', solvph, ' --> ', o_ph
         if (success) then
            write(6, '(a)') 'Success!'
         else
            write(6, '(a)') 'Failure.'
         end if
#endif
    end if ! (mod(my_index, 2) .eq. 0)

    ! Update our solvph if necessary

    if (success) solvph = o_ph

      ! Swap the replica_indexes and group_num if our exchange was successful
    if (success .and. master) then
      call mpi_sendrecv_replace(replica_indexes, remd_dimension, &
                mpi_integer, o_repnum-1, 22, o_repnum-1, 22, &
                remd_comm, stat_array, err_code_mpi)
      call mpi_sendrecv_replace(group_num, remd_dimension, mpi_integer, &
                o_repnum-1, 23, o_repnum-1, 23, remd_comm, stat_array, &
                err_code_mpi)
    end if
    
    ! Reduce our success array to master for printing

    call mpi_reduce(suc_arry, suc_buf, num_replicas, mpi_integer, &
                    mpi_sum, 0, remd_comm, err_code_mpi)
    
    ! Generate our new pH table

    call mpi_gather(solvph, 1, mpi_double_precision, new_ph, 1, &
                    mpi_double_precision, 0, remd_comm, err_code_mpi)

  end if ! (master)

  ! Tell the rest of our replica about our (maybe) new pH

  call mpi_bcast(solvph, 1, mpi_double_precision, 0, pmemd_comm, err_code_mpi)

  if (remd_master) &
    exchange_successes(ph_dim,:) = exchange_successes(ph_dim,:) + suc_buf(:)

  ! Write to the remlog

  if (remd_master) then

    write(remlog, '(a,i8)') '# exchange ', mdloop
    do i = 1, num_replicas
      success_ratio = dble(exchange_successes(ph_dim, index_list(i))) / &
                      dble(mdloop) * dble(remd_dimension) * 2
      write(remlog, '(i6,x,i7,x,2f7.3,x,f8.4)') &
            i, &
            prot_table(i), &
            all_ph(i), &
            new_ph(i), &
            success_ratio
    end do

  end if

  jumpright = .not. jumpright

end subroutine ph_remd_exchange

!*******************************************************************************
!
! Subroutine: setup_remd_randgen
!
! Description: 
!
! Initializes the REMD random number stream. This *should* go in remd.fpp
! (subroutine remd_setup), but the Intel compilers fail when it's taken from a
! different module. I think this is caused because the size of random_state is
! not a multiple of its largest member (causes a warning), so ifort confuses the
! types. This seems to me to be a compiler bug, but oh well...
!
!*******************************************************************************

subroutine setup_remd_randgen

  use mdin_ctrl_dat_mod, only : ig

  implicit none

  integer rand_seed
  integer my_ig
  integer orepnum ! other replica number

  integer, dimension(mpi_status_size) :: stat_array ! needed for mpi_recv

  ! Only master nodes do exchanges, so seed the generator to protect against
  ! calling amrand on non-masters, but don't go through the logic below since
  ! the necessary communicator is not set up

  if (.not. master) then
    call amrset_gen(remd_rand_gen, ig + rand_seed)
    return
  end if

  ! sander and pmemd use different strategies when it comes to exchanging. While
  ! sander has only even replicas exchange 'up' and 'down', pmemd has every
  ! replica exchange 'up', and alternates which replica exchanges. Therefore,
  ! the random number streams between adjacent replicas here that are used for 
  ! exchange attempt evaluatinos should be identical so that the random numbers
  ! used in REMD evaluations are synchronized between pmemd and sander. This way
  ! all results will be identical b/w the two programs and they will be easier
  ! to validate

  my_ig = ig

  if (mod(repnum, 2) .eq. 1) then
    ! This is actually an 'even' replica when indexing from 0
    if (repnum .eq. 1) then
      orepnum = numgroups
    else
      orepnum = repnum - 1
    end if
    call mpi_send(my_ig, 1, mpi_integer, orepnum-1, 10, &
                  pmemd_master_comm, err_code_mpi)
  else
    if (repnum .eq. numgroups) then
      orepnum = 1
    else
      orepnum = repnum + 1
    end if
    call mpi_recv(my_ig, 1, mpi_integer, orepnum-1, 10, pmemd_master_comm, &
                  stat_array, err_code_mpi)
  end if
  rand_seed = mod(repnum-1, numgroups) / 2 * 2

  call amrset_gen(remd_rand_gen, my_ig + rand_seed)

  return

end subroutine setup_remd_randgen

!*******************************************************************************
!
! Subroutine: rxsgld_exchange
!
! Description: Performs the SGLD replica exchange attempt. It allow replicas with
!     different guiding effects and/or different temperatures to exchange.  It is
!     an extended version of temperature_exchange.
!
!*******************************************************************************
subroutine rxsgld_exchange(atm_cnt, my_atm_lst, amass,vel, my_pot_ene_tot, t_dim, &
                                num_replicas, actual_temperature,    &
                                print_exch_data, mdloop)

  use mdin_ctrl_dat_mod, only : temp0,sgft,tempsg
  use pmemd_lib_mod,     only : mexit
  use sgld_mod, only :  tsgset, temprxlf, epotlf, avgtlf, &
                    avgeflf, avgefhf, avgcflf, avgcfhf, sgld_exchg, &
                    rxsgld_scale

  implicit none
 
! Passed variables

  integer, intent(in)             :: atm_cnt
  integer                         :: my_atm_lst(*)
  integer, intent(in)             :: num_replicas
  integer, intent(in)             :: mdloop

  double precision, intent(in)    :: my_pot_ene_tot
  double precision, intent(inout) :: amass(atm_cnt)
  double precision, intent(inout) :: vel(3, atm_cnt)
  double precision, intent(in)    :: actual_temperature
  
  integer, intent(in)             :: t_dim ! which dimension T-exchanges are
  logical, intent(in)             :: print_exch_data

! Local variables

  type :: exchange_data ! facilitates gather for remlog
    sequence
    double precision :: temp0
    double precision :: sgft
    double precision :: tempsg
    double precision :: new_temp0
    double precision :: new_sgft
    double precision :: new_tempsg
    double precision :: real_temp
    double precision :: pot_ene_tot
    double precision :: epotlf
    double precision :: avgtlf
    double precision :: avgeflf
    double precision :: avgefhf
    double precision :: avgcflf
    double precision :: avgcfhf
    double precision :: scaling
    double precision :: scalsg
    double precision :: exchange
  end type exchange_data

  integer, parameter :: SIZE_EXCHANGE_DATA = 17 ! for mpi_gather

  type (exchange_data) :: exch_buffer(num_replicas) ! gather buffer

  type (exchange_data) :: my_exch_data ! stores data for this replica
  type (exchange_data) :: exch_data_tbl(num_replicas) ! exch_data for all reps

   double precision  :: myscaling,o_scaling,myscalsg, o_scalsg
   double precision  :: o_sglf, sglf, o_sghf, sghf

   double precision  :: delta
  double precision  :: metrop
  double precision  :: random_value
  double precision  :: success_ratio

  integer           :: neighbor_rank
  integer           :: i,index_base
  integer           :: success_array(num_replicas)
  integer           :: success_buf(num_replicas)

  ! For exchanging replica positions in ladders upon successful exchanges
  integer           :: group_num_buffer(remd_dimension)
  integer           :: replica_indexes_buffer(remd_dimension)

  logical           :: success
  logical           :: i_do_exchg

  integer, dimension(mpi_status_size) :: stat_array ! needed for mpi_recv

! Explanation of local variables:
!
! delta: The calculated DELTA value used as part of the detailed balance calc.
!
! metrop: The actual value compared to a random_value to decide success
!
! success: Did this exchange attempt succeed?
!
! success_array: buffer for keeping track of which temperatures exchanged
!                successfully. Used in a mpi_reduce at the end
!
! success_buf: buffer for mpi_reduce call on success_array
!
! neighbor_rank: the remd_rank of the replica with whom we need to be doing ALL
!                of our communication. Trying to do everything with the
!                'partners' array is a bit too confusing.

! Set the variables that we know at the beginning

  my_exch_data%temp0       = temp0
  my_exch_data%sgft   = sgft 
  my_exch_data%tempsg  = tempsg
  my_exch_data%new_temp0       = temp0
  my_exch_data%new_sgft   = sgft 
  my_exch_data%new_tempsg  = tempsg
  my_exch_data%real_temp   = actual_temperature
  my_exch_data%pot_ene_tot = my_pot_ene_tot
  my_exch_data%epotlf     = epotlf
  my_exch_data%avgtlf     = avgtlf
  my_exch_data%avgeflf     = avgeflf
  my_exch_data%avgefhf     = avgefhf
  my_exch_data%avgcflf     = avgcflf
  my_exch_data%avgcfhf     = avgcfhf
  my_exch_data%scaling     = -1.0d0
  my_exch_data%scalsg     = -1.0d0
  my_exch_data%exchange   = -1.0d0
  success_array(:)         = 0

  if (master) then
    my_exch_data%exchange   = replica_indexes(t_dim)

    call set_partners(t_dim, num_replicas)

    ! Call amrand_gen here to generate a random number. Sander does this
    ! to pick out structures from a reservoir which pmemd doesn't implement
    ! (yet). To keep the random #s synched, though, and to ensure proper
    ! regression testing, this must be done to keep the random #s synched
    ! b/w sander and pmemd

    call amrand_gen(remd_rand_gen, random_value)

#ifdef VERBOSE_REMD
    write(mdout, '(26("="),a,26("="))') 'RXSGLD EXCHANGE CALCULATION'
    write(mdout, '(a,i10,a,i1)') 'Exch= ', mdloop, ' RREMD= ', 0
    write(mdout, *) 'exdata= ', my_exch_data
#endif

    ! We alternate exchanges: temp up followed by temp down
    if (even_exchanges(t_dim)) then
      i_do_exchg = even_replica(t_dim)
    else
      i_do_exchg = .not. even_replica(t_dim)
    end if

    even_exchanges(t_dim) = .not. even_exchanges(t_dim)

    ! Find out what rank our neighbor is (-1 for indexing from 0), and
    ! If we are controlling the exchange, then our partner will be our 
    ! 2nd partner (above). Otherwise it will be our 1st partner (below)
    
    if (i_do_exchg) then
      neighbor_rank = partners(2) - 1
    else
      neighbor_rank = partners(1) - 1
    end if

    ! Collect potential energies/temperatures

    call mpi_allgather(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                exch_data_tbl, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                remd_comm, err_code_mpi)

    ! Using the base replica to calculate temprxlf
    index_base=-1
    do i=1,num_replicas
      if(temp0==exch_data_tbl(i)%temp0.and.(exch_data_tbl(i)%sgft==0.0d0).and.  &
     (exch_data_tbl(i)%tempsg==0.0d0.or.exch_data_tbl(i)%tempsg==temp0))index_base=i
    enddo
    if(index_base>0)then
       temprxlf=exch_data_tbl(index_base)%avgtlf
    else
       temprxlf=0.0d0
    endif
   
#ifdef VERBOSE_REMD
    write(mdout, '(a,f8.2,f8.4,f8.2,f8.2,2(a7,i2),a7,f10.2)') &
    'Temp,sgft,tempsg,temprxlf= ', &
      temp0,sgft,tempsg,temprxlf, ' Indx= ', replica_indexes(t_dim), ' Rep#= ', &
      remd_rank+1, ' EPot= ', my_pot_ene_tot
    write(mdout, *) 'exdatn= ', exch_data_tbl(neighbor_rank+1)
#endif

    ! Calculate the exchange probability. Note that if delta is negative, then
    ! metrop will be > 1, and it will always succeed (so we don't need an extra
    ! conditional)
    
    if (i_do_exchg) then
#ifdef VERBOSE_REMD
      write(mdout, '(a,f6.2,f6.4,f6.2,2(a7,i2),a7,f10.2)') 'Partner Temp,sgft,tempsg= ', &
        exch_data_tbl(neighbor_rank+1)%temp0,exch_data_tbl(neighbor_rank+1)%sgft,&
        exch_data_tbl(neighbor_rank+1)%tempsg, ' Indx= ', neighbor_rank+1, &
        ' Rep#= ', neighbor_rank + 1, ' EPot= ', &
        exch_data_tbl(neighbor_rank+1)%pot_ene_tot
#endif

     ! delta = (my_exch_data%pot_ene_tot - &
     !          exch_data_tbl(neighbor_rank+1)%pot_ene_tot) * &
     !         (my_exch_data%temp0 - &
     !          exch_data_tbl(neighbor_rank+1)%temp0) * ONEKB /  &
     !         (my_exch_data%temp0 * exch_data_tbl(neighbor_rank+1)%temp0)
            ! Replica exchange self-guided Langevin dynamics
            !  Works also for temperature-based replica exchange
           sghf = avgefhf * avgcfhf * 503.01d0 / my_exch_data%temp0
           sglf = avgeflf * avgcflf * 503.01d0 / my_exch_data%temp0 - sghf
           o_sghf = exch_data_tbl(neighbor_rank+1)%avgefhf * &
           exch_data_tbl(neighbor_rank+1)%avgcfhf * &
           503.01d0/ exch_data_tbl(neighbor_rank+1)%temp0
           o_sglf = exch_data_tbl(neighbor_rank+1)%avgeflf *  &
           exch_data_tbl(neighbor_rank+1)%avgcflf * &
                    503.01d0 / exch_data_tbl(neighbor_rank+1)%temp0 - o_sghf
            
           delta =(sglf-o_sglf)*(exch_data_tbl(neighbor_rank+1)%epotlf - epotlf)  &
                 + (sghf - o_sghf) * (exch_data_tbl(neighbor_rank+1)%pot_ene_tot - my_pot_ene_tot)

      metrop = exp(-delta)

      call amrand_gen(remd_rand_gen, random_value)

      success = random_value .lt. metrop

      ! Let my neighbor know about our success

      call mpi_send(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, remd_comm, err_code_mpi)

      if (success) then
        success_array(replica_indexes(t_dim)) = 1
        myscaling = sqrt(exch_data_tbl(neighbor_rank+1)%temp0 / my_exch_data%temp0)
        myscalsg = sqrt(exch_data_tbl(neighbor_rank+1)%avgtlf / my_exch_data%avgtlf)
        my_exch_data%new_temp0 = exch_data_tbl(neighbor_rank+1)%temp0
        my_exch_data%new_sgft = exch_data_tbl(neighbor_rank+1)%sgft
        my_exch_data%new_tempsg = exch_data_tbl(neighbor_rank+1)%tempsg
        my_exch_data%scaling = myscaling
        my_exch_data%scalsg = myscalsg
        o_scaling = 1.0d0 / myscaling
        o_scalsg = 1.0d0 / myscalsg
      end if

#ifdef VERBOSE_REMD
      write(mdout,'(a8,E16.6,a8,E16.6,a12,f10.2,f10.2)') &
               "Metrop= ",metrop," delta= ",delta," o_scaling,o_scalsg= ", &
               o_scaling,o_scalsg
#endif

    else

#ifdef VERBOSE_REMD
      ! Even if we don't control exchange, still print REMD info to mdout

      write(mdout, '(a,f6.2,f6.4,f6.2,2(a7,i2),a7,f10.2)') 'Partner Temp,sgft,tempsg= ', &
        exch_data_tbl(neighbor_rank+1)%temp0,exch_data_tbl(neighbor_rank+1)%sgft, &
        exch_data_tbl(neighbor_rank+1)%tempsg, ' Indx= ', neighbor_rank+1, &
        ' Rep#= ', neighbor_rank + 1, ' EPot= ', &
        exch_data_tbl(neighbor_rank+1)%pot_ene_tot
      write(mdout, '(a)') 'Not controlling exchange.'
#endif

      call amrand_gen(remd_rand_gen, random_value) ! to stay synchronized

      ! Get the message from the exchanging replica about our success
      call mpi_recv(success, 1, mpi_logical, neighbor_rank, &
                    remd_tag, remd_comm, stat_array, err_code_mpi)

      ! We scale velocities only if we succeed
      if (success) then
        myscaling = sqrt(exch_data_tbl(neighbor_rank+1)%temp0 / my_exch_data%temp0)
        myscalsg = sqrt(exch_data_tbl(neighbor_rank+1)%avgtlf / my_exch_data%avgtlf)
        my_exch_data%new_temp0 = exch_data_tbl(neighbor_rank+1)%temp0
        my_exch_data%new_sgft = exch_data_tbl(neighbor_rank+1)%sgft
        my_exch_data%new_tempsg = exch_data_tbl(neighbor_rank+1)%tempsg
        my_exch_data%scaling = myscaling
        my_exch_data%scalsg = myscalsg
        o_scaling = 1.0d0 / myscaling
        o_scalsg = 1.0d0 / myscalsg
      end if

    end if

#ifdef VERBOSE_REMD
    write(mdout,'(a8,E16.6,a12,f10.2,f10.2,a10,L1)') &
      'Rand=   ', random_value, ' MyScaling,myscalsg= ', myscaling, &
      myscalsg,' Success= ', success
    
    write(mdout,'(24("="),a,24("="))') "END RXSGLD EXCHANGE CALCULATION"
#endif

  end if ! master

  ! We'll broadcast all of the my_exch_data here once. From that, we can tell
  ! if the exchange succeeded based on whether or not the temperatures have
  ! changed. Then make sure we update temp0 so the simulation reflects that.

  call mpi_bcast(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                 0, pmemd_comm, err_code_mpi)

  if (master.and.success)then
     ! exchange all SGLD common data block
     call sgld_exchg(neighbor_rank)
  endif
  if (.not. master ) &
    success = my_exch_data%scaling.gt.TINY
  
#ifdef VERBOSE_REMD
  if (master) &
    write(mdout,*)'rank, sucess: ', &
       my_exch_data%exchange,replica_indexes(t_dim),neighbor_rank+1
   if (master) &
     write(mdout, '(2(a,i4,f6.2,f6.4,f6.2))') &
      'RXSGLD: checking to see if TEMP0,SGFT,TEMPSG has changed: ', &
      replica_indexes(t_dim),temp0,sgft,tempsg, '->', &
     int(my_exch_data%exchange),my_exch_data%temp0,my_exch_data%sgft,my_exch_data%tempsg
#endif

  if (success)then
   ! ---=== RXSGLD property SCALING ===---
     temp0 = my_exch_data%new_temp0
     sgft = my_exch_data%new_sgft
     tempsg = my_exch_data%new_tempsg
     call rxsgld_scale(atm_cnt,my_atm_lst, my_exch_data%scaling,my_exch_data%scalsg, amass, vel)
  endif
  ! We want to keep track of our exchange successes, so gather them
  ! from all of the masters here, and increment them in exchange_successes
  ! We do this whether or not we print so we can keep track of exchange
  ! success ratios. For the index and exchange data, only gather those to
  ! the remd_master if we actually plan to print. This should be much
  ! more efficient for high exchange attempts.

  ! Then recollect the temperatures as some may have changed, and reset our
  ! partners.
  if (master) then

    if (print_exch_data) then
      call mpi_gather(my_exch_data, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                      exch_buffer, SIZE_EXCHANGE_DATA, mpi_double_precision, &
                      0, remd_comm, err_code_mpi)
    end if

    call mpi_reduce(success_array, success_buf, num_replicas, mpi_integer, &
                    mpi_sum, 0, remd_comm, err_code_mpi)

    if (remd_master) &
      exchange_successes(t_dim,:) = exchange_successes(t_dim,:) + success_buf(:)

    ! Increment exchange counter, and swap our replica ranks and numbers with 
    ! our neighbor if our attempt succeeded, since we effectively swapped places
    ! in this array. That's because we swap temperatures and not structures due
    ! to performance considerations.

    if (success) then
      call mpi_sendrecv(group_num, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, &
                        group_num_buffer, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, remd_comm, &
                        stat_array, err_code_mpi)
      call mpi_sendrecv(replica_indexes, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, &
                        replica_indexes_buffer, remd_dimension, mpi_integer, &
                        neighbor_rank, remd_tag, remd_comm, &
                        stat_array, err_code_mpi)
      group_num = group_num_buffer
      replica_indexes = replica_indexes_buffer
    end if
  end if
 
  ! If we're doing NMROPT, then we need to make sure we re-read in the TEMP0,
  ! since this subroutine overwrites TEMP0 to the original since you can modify
  ! temperature with NMROPT settings.

  remd_modwt = .true.

  ! Write the data to the remlog

  ! remd_method == -1 means we're doing multi-D REMD.
  if (remd_method .ne. -1) then
    if (print_exch_data .and. remd_master) then
      write(remlog, '(a,i8)') '# exchange ', mdloop
      do i = 1, num_replicas
        success_ratio = dble(exchange_successes(t_dim, index_list(i))) / &
                        dble(mdloop) * dble(remd_dimension) * 2
        write(remlog, '(i4,i4,2f8.4,2f8.2,e14.6,f8.4)') &
              i, int(exch_buffer(i)%exchange),&
              exch_buffer(i)%scaling, &
              exch_buffer(i)%scalsg, &
              exch_buffer(i)%real_temp, &
              exch_buffer(i)%avgtlf, &
              exch_buffer(i)%pot_ene_tot, &
              success_ratio
      end do
    end if
  else if (remd_master) then
    do i = 1, num_replicas
      success_ratio = dble(exchange_successes(t_dim, index_list(i))) / &
                      dble(mdloop) * dble(remd_dimension) * 2
      multid_print_data(i)%scaling       = exch_buffer(i)%scaling
      multid_print_data(i)%real_temp     = exch_buffer(i)%real_temp
      multid_print_data(i)%pot_ene_tot   = exch_buffer(i)%pot_ene_tot
      multid_print_data(i)%temp0         = exch_buffer(i)%avgtlf
      multid_print_data(i)%new_temp0     = exch_buffer(i)%scalsg
      multid_print_data(i)%group_num     = group_num(t_dim)
      multid_print_data(i)%success_ratio = success_ratio
      multid_print_data(i)%num_rep       = num_replicas
    end do
  end if ! remd_method .ne. -1
#ifdef VERBOSE_REMD
  if (master) &
    write(mdout, '(2(a,i8))') &
      'RXSGLD: checking  ', &
      replica_indexes(t_dim), '->', &
     int(my_exch_data%exchange)
#endif

  return

end subroutine rxsgld_exchange

#endif /* MPI */

end module remd_exchg_mod
