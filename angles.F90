#include "copyright.i"

!*******************************************************************************
!
! Module:  angles_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module angles_mod

  use gbl_datatypes_mod

  implicit none

! The following are derived from prmtop angle info:

! Added by Lijiang Yang
  integer, save                         :: cit_nthethslu, cit_nthetaslu
  integer, save                         :: cit_nthethsol, cit_nthetasol
  integer, save                         :: cit_nthethinter, cit_nthetainter
!-----------------------------------------------------------------

! Added by Lijiang Yang
  type(angle_rec), allocatable, save    :: cit_angleslu(:)
  type(angle_rec), allocatable, save    :: cit_anglesol(:)
  type(angle_rec), allocatable, save    :: cit_angleinter(:)
!-----------------------------------------------------------------

contains

!*******************************************************************************
!
! Subroutine:  angles_setup
!
! Description:  <TBS>
!
!*******************************************************************************

subroutine angles_setup(num_ints, num_reals, use_atm_map)

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
! Added by Lijiang Yang
  type(angle_rec)       :: angles_copyslu(nthethslu + nthetaslu)
  type(angle_rec)       :: angles_copysol(nthethsol + nthetasol)
  type(angle_rec)       :: angles_copyinter(nthethinter + nthetainter)
!-----------------------------------------------------------------
  integer               :: atm_i, atm_j, atm_k, angles_idx
  integer               :: my_angle_cnt_slu, my_angle_cnt_sol, my_angle_cnt_inter

! This routine can handle reallocation, and thus can be called multiple
! times.

! Find all angles for which this process owns either atom:

  my_angle_cnt_slu = 0

  if ( nthethslu .gt. 0 ) then
    do angles_idx = 1, nthethslu

      atm_i = gbl_angleslu(angles_idx)%atm_i
      atm_j = gbl_angleslu(angles_idx)%atm_j
      atm_k = gbl_angleslu(angles_idx)%atm_k

#if defined(MPI) && !defined(CUDA)
      if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
        my_angle_cnt_slu = my_angle_cnt_slu + 1
        angles_copyslu(my_angle_cnt_slu) = gbl_angleslu(angles_idx)
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
      else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
        my_angle_cnt_slu = my_angle_cnt_slu + 1
        angles_copyslu(my_angle_cnt_slu) = gbl_angleslu(angles_idx)
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
      else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
        my_angle_cnt_slu = my_angle_cnt_slu + 1
        angles_copyslu(my_angle_cnt_slu) = gbl_angleslu(angles_idx)
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
      end if
#else
      my_angle_cnt_slu = my_angle_cnt_slu + 1
      angles_copyslu(my_angle_cnt_slu) = gbl_angleslu(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
#endif

    end do
  endif

  cit_nthethslu = my_angle_cnt_slu

  my_angle_cnt_sol = 0

  if (nthethsol .gt. 0) then
     do angles_idx = 1, nthethsol

        atm_i = gbl_anglesol(angles_idx)%atm_i
        atm_j = gbl_anglesol(angles_idx)%atm_j
        atm_k = gbl_anglesol(angles_idx)%atm_k

#if defined(MPI) && !defined(CUDA)
        if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
          my_angle_cnt_sol = my_angle_cnt_sol + 1
          angles_copysol(my_angle_cnt_sol) = gbl_anglesol(angles_idx)
          use_atm_map(atm_i) = 1
          use_atm_map(atm_j) = 1
          use_atm_map(atm_k) = 1
        else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
          my_angle_cnt_sol = my_angle_cnt_sol + 1
          angles_copysol(my_angle_cnt_sol) = gbl_anglesol(angles_idx)
          use_atm_map(atm_i) = 1
          use_atm_map(atm_j) = 1
          use_atm_map(atm_k) = 1
        else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
          my_angle_cnt_sol = my_angle_cnt_sol + 1
          angles_copysol(my_angle_cnt_sol) = gbl_anglesol(angles_idx)
          use_atm_map(atm_i) = 1
          use_atm_map(atm_j) = 1
          use_atm_map(atm_k) = 1
        end if
#else
        my_angle_cnt_sol = my_angle_cnt_sol + 1
        angles_copysol(my_angle_cnt_sol) = gbl_anglesol(angles_idx)
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
#endif

      end do
  endif

  cit_nthethsol = my_angle_cnt_sol

  my_angle_cnt_inter = 0

  if (nthethinter .gt. 0) then
     do angles_idx = 1, nthethinter

        atm_i = gbl_angleinter(angles_idx)%atm_i
        atm_j = gbl_angleinter(angles_idx)%atm_j
        atm_k = gbl_angleinter(angles_idx)%atm_k

#if defined(MPI) && !defined(CUDA)
        if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
          my_angle_cnt_inter = my_angle_cnt_inter + 1
          angles_copyinter(my_angle_cnt_inter) = gbl_angleinter(angles_idx)
          use_atm_map(atm_i) = 1
          use_atm_map(atm_j) = 1
          use_atm_map(atm_k) = 1
        else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
          my_angle_cnt_inter = my_angle_cnt_inter + 1
          angles_copyinter(my_angle_cnt_inter) = gbl_angleinter(angles_idx)
          use_atm_map(atm_i) = 1
          use_atm_map(atm_j) = 1
          use_atm_map(atm_k) = 1
        else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
          my_angle_cnt_inter = my_angle_cnt_inter + 1
          angles_copyinter(my_angle_cnt_inter) = gbl_angleinter(angles_idx)
          use_atm_map(atm_i) = 1
          use_atm_map(atm_j) = 1
          use_atm_map(atm_k) = 1
        end if
#else
        my_angle_cnt_inter = my_angle_cnt_inter + 1
        angles_copyinter(my_angle_cnt_inter) = gbl_angleinter(angles_idx)
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
#endif

      end do
  endif

  cit_nthethinter = my_angle_cnt_inter

  if ( nthetaslu .gt. 0 ) then
    do angles_idx = anglea_idx_slu, (anglea_idx_slu + nthetaslu - 1 ) 
       atm_i = gbl_angleslu(angles_idx)%atm_i
       atm_j = gbl_angleslu(angles_idx)%atm_j
       atm_k = gbl_angleslu(angles_idx)%atm_k

#if defined(MPI) && !defined(CUDA)
       if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
         my_angle_cnt_slu = my_angle_cnt_slu + 1
         angles_copyslu(my_angle_cnt_slu) = gbl_angleslu(angles_idx)
         use_atm_map(atm_i) = 1
         use_atm_map(atm_j) = 1
         use_atm_map(atm_k) = 1
       else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
         my_angle_cnt_slu = my_angle_cnt_slu + 1
         angles_copyslu(my_angle_cnt_slu) = gbl_angleslu(angles_idx)
         use_atm_map(atm_i) = 1
         use_atm_map(atm_j) = 1
         use_atm_map(atm_k) = 1
       else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
         my_angle_cnt_slu = my_angle_cnt_slu + 1
         angles_copyslu(my_angle_cnt_slu) = gbl_angleslu(angles_idx)
         use_atm_map(atm_i) = 1
         use_atm_map(atm_j) = 1
         use_atm_map(atm_k) = 1
       end if
#else
       my_angle_cnt_slu = my_angle_cnt_slu + 1
       angles_copyslu(my_angle_cnt_slu) = gbl_angleslu(angles_idx)
       use_atm_map(atm_i) = 1
       use_atm_map(atm_j) = 1
       use_atm_map(atm_k) = 1
#endif

     end do
  endif

  cit_nthetaslu = my_angle_cnt_slu - cit_nthethslu

  if ( nthetasol .gt. 0 ) then
    do angles_idx = anglea_idx_sol, (anglea_idx_sol + nthetasol - 1)

      atm_i = gbl_anglesol(angles_idx)%atm_i
      atm_j = gbl_anglesol(angles_idx)%atm_j
      atm_k = gbl_anglesol(angles_idx)%atm_k

#if defined(MPI) && !defined(CUDA)
      if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
        my_angle_cnt_sol = my_angle_cnt_sol + 1
        angles_copysol(my_angle_cnt_sol) = gbl_anglesol(angles_idx)
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
      else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
        my_angle_cnt_sol = my_angle_cnt_sol + 1
        angles_copysol(my_angle_cnt_sol) = gbl_anglesol(angles_idx)
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
      else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
        my_angle_cnt_sol = my_angle_cnt_sol + 1
        angles_copysol(my_angle_cnt_sol) = gbl_anglesol(angles_idx)
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
      end if
#else
      my_angle_cnt_sol = my_angle_cnt_sol + 1
      angles_copysol(my_angle_cnt_sol) = gbl_anglesol(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
#endif

    end do
  endif

  cit_nthetasol = my_angle_cnt_sol - cit_nthethsol

  if ( nthetainter .gt. 0 ) then
    do angles_idx = anglea_idx_inter, (anglea_idx_inter + nthetainter - 1)

      atm_i = gbl_angleinter(angles_idx)%atm_i
      atm_j = gbl_angleinter(angles_idx)%atm_j
      atm_k = gbl_angleinter(angles_idx)%atm_k

#if defined(MPI) && !defined(CUDA)
      if (gbl_atm_owner_map(atm_i) .eq. mytaskid) then
        my_angle_cnt_inter = my_angle_cnt_inter + 1
        angles_copyinter(my_angle_cnt_inter) = gbl_angleinter(angles_idx)
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
      else if (gbl_atm_owner_map(atm_j) .eq. mytaskid) then
        my_angle_cnt_inter = my_angle_cnt_inter + 1
        angles_copyinter(my_angle_cnt_inter) = gbl_angleinter(angles_idx)
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
      else if (gbl_atm_owner_map(atm_k) .eq. mytaskid) then
        my_angle_cnt_inter = my_angle_cnt_inter + 1
        angles_copyinter(my_angle_cnt_inter) = gbl_angleinter(angles_idx)
        use_atm_map(atm_i) = 1
        use_atm_map(atm_j) = 1
        use_atm_map(atm_k) = 1
      end if
#else
      my_angle_cnt_inter = my_angle_cnt_inter + 1
      angles_copyinter(my_angle_cnt_inter) = gbl_angleinter(angles_idx)
      use_atm_map(atm_i) = 1
      use_atm_map(atm_j) = 1
      use_atm_map(atm_k) = 1
#endif

    end do
  endif

  cit_nthetainter = my_angle_cnt_inter - cit_nthethinter

  if (my_angle_cnt_slu .gt. 0) then
    if (allocated(cit_angleslu)) then
      if (size(cit_angleslu) .lt. my_angle_cnt_slu) then
        num_ints = num_ints - size(cit_angleslu) * angle_rec_ints
        deallocate(cit_angleslu)
        allocate(cit_angleslu(my_angle_cnt_slu), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_angleslu) * angle_rec_ints
      end if
    else
      allocate(cit_angleslu(my_angle_cnt_slu), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_angleslu) * angle_rec_ints
    end if
    cit_angleslu(1:my_angle_cnt_slu) = angles_copyslu(1:my_angle_cnt_slu)
  end if

  if (my_angle_cnt_sol .gt. 0) then
    if (allocated(cit_anglesol)) then
      if (size(cit_anglesol) .lt. my_angle_cnt_sol) then
        num_ints = num_ints - size(cit_anglesol) * angle_rec_ints
        deallocate(cit_anglesol)
        allocate(cit_anglesol(my_angle_cnt_sol), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_anglesol) * angle_rec_ints
      end if
    else
      allocate(cit_anglesol(my_angle_cnt_sol), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_anglesol) * angle_rec_ints
    end if
    cit_anglesol(1:my_angle_cnt_sol) = angles_copysol(1:my_angle_cnt_sol)
  end if

  if (my_angle_cnt_inter .gt. 0) then
    if (allocated(cit_angleinter)) then
      if (size(cit_angleinter) .lt. my_angle_cnt_inter) then
        num_ints = num_ints - size(cit_angleinter) * angle_rec_ints
        deallocate(cit_angleinter)
        allocate(cit_angleinter(my_angle_cnt_inter), stat = alloc_failed)
        if (alloc_failed .ne. 0) call setup_alloc_error
        num_ints = num_ints + size(cit_angleinter) * angle_rec_ints
      end if
    else
      allocate(cit_angleinter(my_angle_cnt_inter), stat = alloc_failed)
      if (alloc_failed .ne. 0) call setup_alloc_error
      num_ints = num_ints + size(cit_angleinter) * angle_rec_ints
    end if
    cit_angleinter(1:my_angle_cnt_inter) = angles_copyinter(1:my_angle_cnt_inter)
  end if

#ifdef CUDA
    call gpu_angles_setup(my_angle_cnt, cit_ntheth, cit_angle, gbl_teq, gbl_tk)
#endif

!BEGIN DBG
! write(0,*)'task,cit_ntheth,cit_ntheta,', mytaskid, &
!            cit_ntheth, cit_ntheta
! write(0,*)'task,    ntheth,    ntheta', mytaskid, &
!                ntheth,     ntheta
!END DBG
  return

end subroutine angles_setup

!*******************************************************************************
!
! Subroutine:  get_angle_energy
!
! Description:  Routine to get the bond energies and forces for potentials of
!               the type ct*(t-t0)**2.
!
!*******************************************************************************

subroutine get_angle_energy(angle_cnt, angle, x, frc, subfrc, angle_energy)

  use nmr_calls_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: angle_cnt
  type(angle_rec)       :: angle(*)
  double precision      :: x(3, *)
  double precision      :: frc(3, *)
  double precision      :: subfrc(3, *)
  double precision      :: angle_energy

! Local variables:

  double precision, parameter   :: pt999 = 0.9990d0

  double precision      :: ant
  double precision      :: cst
  double precision      :: dfw
  double precision      :: eaw
  double precision      :: rij, rik, rkj
  double precision      :: xij, yij, zij
  double precision      :: xkj, ykj, zkj
  double precision      :: cii, cik, ckk
  double precision      :: da
  double precision      :: df
  double precision      :: dt1, dt2, dt3, dt4, dt5, dt6
  double precision      :: sth
  integer               :: i, j, k, ic, jn

  angle_energy = 0.0d0

! Grand loop for the angle stuff:

  do jn = 1, angle_cnt

    i = angle(jn)%atm_i
    j = angle(jn)%atm_j
    k = angle(jn)%atm_k
    ic = angle(jn)%parm_idx

! Calculation of the angle:

    xij = x(1, i) - x(1, j)
    xkj = x(1, k) - x(1, j)

    yij = x(2, i) - x(2, j)
    ykj = x(2, k) - x(2, j)

    zij = x(3, i) - x(3, j)
    zkj = x(3, k) - x(3, j)

    rij = xij * xij + yij * yij + zij * zij
    rkj = xkj * xkj + ykj * ykj + zkj * zkj

    rik = sqrt(rij * rkj)

    cst = min(pt999, max(-pt999, (xij * xkj + yij * ykj + zij * zkj) / rik))

    ant = acos(cst)

! Calculation of the energy and deriv:

    da = ant - gbl_teq(ic)

#ifdef MPI
! If eadev is ever supported under mpi, you will need to split the contribution.
#else
    eadev = eadev + da * da         ! For rms deviation from ideal angles.
#endif

    df = gbl_tk(ic) * da
    eaw = df * da

    if (ti_mode .ne. 0) then
      ti_region = 0
      if (ti_lst(1,i) .eq. 1 .or. ti_lst(1,j) .eq. 1 &
          .or. ti_lst(1,k) .eq. 1) then
        ti_region = 1
      else if (ti_lst(2,i) .eq. 1 .or. ti_lst(2,j) .eq. 1 &
               .or. ti_lst(2,k) .eq. 1) then
        ti_region = 2
      end if
      if (ti_region .gt. 0) then 
        if (ti_mode .ne. 1) then
          ti_lscale = (ti_sc_lst(i) .lt. 2 .and. &
                       ti_sc_lst(j) .lt. 2 .and. &
                       ti_sc_lst(k) .lt. 2) 
        end if
#ifdef MPI
        if (gbl_atm_owner_map(i) .eq. mytaskid) then
#endif        
          if (ti_mode .eq. 1 .or. ti_lscale) then
            call ti_update_ene(eaw, si_angle_ene, ti_region)
            eaw = eaw * ti_weights(ti_region)
          else
            ti_ene(ti_region,si_angle_ene) = &
              ti_ene(ti_region,si_angle_ene) + eaw
            eaw = 0.d0
          end if
#ifdef MPI
        end if
#endif
        if (ti_mode .eq. 1 .or. ti_lscale) df = df * ti_weights(ti_region)
      end if
    end if

    dfw = -(df + df) / sin(ant)

! Calculation of the force: 

    cik = dfw / rik
    sth = dfw * cst
    cii = sth / rij
    ckk = sth / rkj

    dt1 = cik * xkj - cii * xij
    dt2 = cik * ykj - cii * yij
    dt3 = cik * zkj - cii * zij
    dt4 = cik * xij - ckk * xkj
    dt5 = cik * yij - ckk * ykj
    dt6 = cik * zij - ckk * zkj
    
#ifdef MPI
    ! We use atm_i to determine who sums up the energy...
    if (gbl_atm_owner_map(i) .eq. mytaskid) then

      angle_energy = angle_energy + eaw

!      frc(1, i) = frc(1, i) - dt1
!      frc(2, i) = frc(2, i) - dt2
!      frc(3, i) = frc(3, i) - dt3
      subfrc(1, i) = subfrc(1, i) - dt1
      subfrc(2, i) = subfrc(2, i) - dt2
      subfrc(3, i) = subfrc(3, i) - dt3

    end if

    if (gbl_atm_owner_map(j) .eq. mytaskid) then

!      frc(1, j) = frc(1, j) + dt1 + dt4
!      frc(2, j) = frc(2, j) + dt2 + dt5
!      frc(3, j) = frc(3, j) + dt3 + dt6
      subfrc(1, j) = subfrc(1, j) + dt1 + dt4
      subfrc(2, j) = subfrc(2, j) + dt2 + dt5
      subfrc(3, j) = subfrc(3, j) + dt3 + dt6

    end if

    if (gbl_atm_owner_map(k) .eq. mytaskid) then

!      frc(1, k) = frc(1, k) - dt4
!      frc(2, k) = frc(2, k) - dt5
!      frc(3, k) = frc(3, k) - dt6
      subfrc(1, k) = subfrc(1, k) - dt4
      subfrc(2, k) = subfrc(2, k) - dt5
      subfrc(3, k) = subfrc(3, k) - dt6

    end if

#else
!    frc(1, i) = frc(1, i) - dt1
!    frc(2, i) = frc(2, i) - dt2
!    frc(3, i) = frc(3, i) - dt3
!    frc(1, j) = frc(1, j) + dt1 + dt4
!    frc(2, j) = frc(2, j) + dt2 + dt5
!    frc(3, j) = frc(3, j) + dt3 + dt6
!    frc(1, k) = frc(1, k) - dt4
!    frc(2, k) = frc(2, k) - dt5
!    frc(3, k) = frc(3, k) - dt6
    subfrc(1, i) = subfrc(1, i) - dt1
    subfrc(2, i) = subfrc(2, i) - dt2
    subfrc(3, i) = subfrc(3, i) - dt3
    subfrc(1, j) = subfrc(1, j) + dt1 + dt4
    subfrc(2, j) = subfrc(2, j) + dt2 + dt5
    subfrc(3, j) = subfrc(3, j) + dt3 + dt6
    subfrc(1, k) = subfrc(1, k) - dt4
    subfrc(2, k) = subfrc(2, k) - dt5
    subfrc(3, k) = subfrc(3, k) - dt6

    angle_energy = angle_energy + eaw
#endif
        
  end do

  return

end subroutine get_angle_energy

end module angles_mod
