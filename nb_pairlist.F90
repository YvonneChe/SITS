#include "copyright.i"

!*******************************************************************************
!
! Module:  nb_pairlist_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module nb_pairlist_mod

  implicit none

! Data that is not broadcast:

  double precision, save, allocatable   :: gbl_atm_saved_crd(:,:)
  double precision, save, allocatable   :: gbl_saved_imgcrd(:,:)
! Added by Lijiang Yang for subunit energy and force
  integer, allocatable, save            :: gbl_sluipairs(:)
  integer, allocatable, save            :: gbl_solipairs(:)
  integer, allocatable, save            :: gbl_interipairs(:)
! -------------------------------------------------------------
  double precision, save                :: gbl_saved_box(3)
  integer, save                         :: ipairs_maxsize
  ! Pairlist definitions for TI
  integer, parameter    :: nb_list_null = 0 
  integer, parameter    :: nb_list_common = 1
  integer, parameter    :: nb_list_linear_V0 = 2
  integer, parameter    :: nb_list_linear_V1 = 3  
  integer, parameter    :: nb_list_sc_common_V0 = 4
  integer, parameter    :: nb_list_sc_common_V1 = 5
  integer, parameter    :: nb_list_sc_sc_V0 = 6
  integer, parameter    :: nb_list_sc_sc_V1 = 7
  integer, parameter    :: nb_list_cnt = 7
contains

!*******************************************************************************
!
! Subroutine:  alloc_nb_pairlist_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_nb_pairlist_mem(atm_cnt, cutlist, num_ints, num_reals)

  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  double precision, intent(in)  :: cutlist

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals

! Local variables:

  integer                       :: alloc_failed
  integer, parameter            :: ipairs_size_min = 1000
!  double precision, parameter   :: ipairs_size_coef = 0.225d0 
  double precision, parameter   :: ipairs_size_coef = 0.5d0 

  allocate(gbl_atm_saved_crd(3, atm_cnt), &
           gbl_saved_imgcrd(3, atm_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(gbl_atm_saved_crd) + &
                          size(gbl_saved_imgcrd)

! Allocate pairs list:

! The following ipairs_maxsize calc is a crude heuristic assuming that the
! number of nonbonded pairs roughly scales with the cutoff volume and
! the number of atoms in the system.  If MPI is running, the total is
! divided up among the processes.  The 2 or 3 is for the counters and flags at
! the head of each atom sublist.

  ipairs_maxsize = int((ipairs_size_coef * cutlist**3 + 3.d0) * dble(atm_cnt))

#ifdef MPI
  if (numtasks .gt. 0) ipairs_maxsize = ipairs_maxsize / numtasks
#endif /* MPI */

  if (ipairs_maxsize .lt. ipairs_size_min) ipairs_maxsize = ipairs_size_min

! Modified by Lijiang Yang for subunit energy and force -----
  allocate(gbl_sluipairs(ipairs_maxsize), &
           gbl_solipairs(ipairs_maxsize), &
           gbl_interipairs(ipairs_maxsize), &
           stat = alloc_failed)
!  ----------------------------------------------------------
!  Original form --------------------------------------------
!  allocate(gbl_ipairs(ipairs_maxsize), stat = alloc_failed)
!  ----------------------------------------------------------

! Added by Lijiang Yang (I think it should be presented)
  num_ints = num_ints + size(gbl_sluipairs) + &
                        size(gbl_solipairs) + &
                        size(gbl_interipairs)

  if (alloc_failed .ne. 0) call setup_alloc_error

  return

end subroutine alloc_nb_pairlist_mem
  
#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_nb_pairlist_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_nb_pairlist_dat(atm_cnt, cutlist)

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  double precision, intent(in)  :: cutlist

! Local variables:

  integer               :: num_ints, num_reals  ! returned values discarded
  
  if (.not. master) then
    num_ints = 0
    num_reals = 0
    call alloc_nb_pairlist_mem(atm_cnt, cutlist, num_ints, num_reals)
  end if

  ! The allocated data is not initialized from the master node.

  return

end subroutine bcast_nb_pairlist_dat
#endif

!*******************************************************************************
!
! Subroutine:   map_pairlist_imgs
!
! Description:  TBS
!
!*******************************************************************************

#ifdef MPI
subroutine map_pairlist_imgs(atm_cnt, flat_cit, fraction, charge, iac, &
                             img_crd, img_qterm, img_iac, mapped_img_cnt, &
                             mapped_img_lst, img_atm_map)
#else
subroutine map_pairlist_imgs(atm_cnt, fraction, charge, iac, &
                             img_crd, img_qterm, img_iac, img_atm_map)
#endif

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
#ifdef MPI
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)
#endif /* MPI */
  double precision      :: fraction(3, atm_cnt)
  double precision      :: charge(atm_cnt)
  integer               :: iac(atm_cnt)
#ifdef MPI
  integer               :: mapped_img_cnt
  integer               :: mapped_img_lst(*)
#endif /* MPI */
  double precision      :: img_crd(3, atm_cnt)
  double precision      :: img_qterm(atm_cnt)
  integer               :: img_iac(atm_cnt)
  integer               :: img_atm_map(atm_cnt)

! Local variables:

! "bkt" refers to a bucket or cell in the flat cit table.
! "bkts" refers to the bucket mapping tables.

  double precision      :: f1, f2, f3
  double precision      :: x_box, y_box, z_box
  double precision      :: ucell_stk(3, 3)
#ifdef MPI
  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
#endif /* MPI */
  integer               :: atm_idx, img_idx
  integer               :: is_orthog_stk
#ifdef MPI
  integer               :: lst_idx
  integer               :: i_bkt_x_idx, i_bkt_y_idx, i_bkt_z_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: img_i_lo, atm_i_lo
  integer               :: i_bkt_lo, i_bkt_hi
  integer               :: i_bkt, cur_bkt, yz_bkt, z_bkt
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
#endif /* MPI */

  x_box = pbc_box(1)
  y_box = pbc_box(2)
  z_box = pbc_box(3)

  ucell_stk(:,:) = ucell(:,:)

  is_orthog_stk = is_orthog

#ifdef MPI

  do lst_idx = 1, mapped_img_cnt
    img_idx = mapped_img_lst(lst_idx)
    img_iac(img_idx) = 0
  end do
  mapped_img_cnt = 0

  if (my_img_lo .gt. my_img_hi) return

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

! Set up the bucket mapping arrays.  Redoing this is only necessary if the
! box dimensions change, but it is cheap and there is a locality benefit
! to having it on the stack.

  call setup_cit_tbl_bkts(x_bkts, y_bkts, z_bkts)

! Get the low and high flat cit bucket indexes for this task:

  call get_flat_cit_idx(my_img_lo, img_atm_map, fraction, i_bkt_lo)
  call get_flat_cit_idx(my_img_hi, img_atm_map, fraction, i_bkt_hi)

! Map all the images you will reference in building the pairlist.
! Note that "mapping" is different than claiming the images as "used".
! When mapped, this just means all the image's data structures are valid.

  do i_bkt = i_bkt_lo, i_bkt_hi

    img_i_lo = flat_cit(i_bkt)%img_lo
    if (img_i_lo .eq. 0) cycle

    ! Get the x,y,z 3d cit indexes to be able to get the bucket ranges
    ! to search:

    atm_i_lo = img_atm_map(img_i_lo)
    i_bkt_x_idx = int(fraction(1, atm_i_lo) * scale_fac_x)
    i_bkt_y_idx = int(fraction(2, atm_i_lo) * scale_fac_y)
    i_bkt_z_idx = int(fraction(3, atm_i_lo) * scale_fac_z)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = i_bkt_x_idx + cit_tbl_x_dim - cit_bkt_delta
    x_bkts_hi = i_bkt_x_idx + cit_tbl_x_dim + cit_bkt_delta
    y_bkts_lo = i_bkt_y_idx + cit_tbl_y_dim - cit_bkt_delta
    y_bkts_hi = i_bkt_y_idx + cit_tbl_y_dim + cit_bkt_delta
    z_bkts_lo = i_bkt_z_idx + cit_tbl_z_dim - cit_bkt_delta
    z_bkts_hi = i_bkt_z_idx + cit_tbl_z_dim
                                        
    z_loop: &
    do z_bkts_idx = z_bkts_lo, z_bkts_hi
      z_bkt = z_bkts(z_bkts_idx)
      y_loop: &
      do y_bkts_idx = y_bkts_lo, y_bkts_hi
        yz_bkt = z_bkt + y_bkts(y_bkts_idx)
        x_loop: &
        do x_bkts_idx = x_bkts_lo, x_bkts_hi
          cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
          img_i_lo = flat_cit(cur_bkt)%img_lo
          if (img_i_lo .eq. 0) cycle
          if (img_iac(img_i_lo) .ne. 0) cycle
          if (is_orthog_stk .ne. 0) then
            do img_idx = img_i_lo, flat_cit(cur_bkt)%img_hi
              atm_idx = img_atm_map(img_idx)
              img_crd(1, img_idx) = fraction(1, atm_idx) * x_box
              img_crd(2, img_idx) = fraction(2, atm_idx) * y_box
              img_crd(3, img_idx) = fraction(3, atm_idx) * z_box
              img_qterm(img_idx) = charge(atm_idx)
              img_iac(img_idx) = iac(atm_idx)
              mapped_img_cnt = mapped_img_cnt + 1
              mapped_img_lst(mapped_img_cnt) = img_idx
            end do
          else
            do img_idx = img_i_lo, flat_cit(cur_bkt)%img_hi
              atm_idx = img_atm_map(img_idx)
              f1 = fraction(1, atm_idx)
              f2 = fraction(2, atm_idx)
              f3 = fraction(3, atm_idx)
              ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
              ! we can simplify the expression in this critical inner loop
              img_crd(1, img_idx) = f1 * ucell_stk(1, 1) + &
                                    f2 * ucell_stk(1, 2) + &
                                    f3 * ucell_stk(1, 3)
              img_crd(2, img_idx) = f2 * ucell_stk(2, 2) + &
                                    f3 * ucell_stk(2, 3)
              img_crd(3, img_idx) = f3 * ucell_stk(3, 3)
              img_qterm(img_idx) = charge(atm_idx)
              img_iac(img_idx) = iac(atm_idx)
              mapped_img_cnt = mapped_img_cnt + 1
              mapped_img_lst(mapped_img_cnt) = img_idx
            end do
          end if
          if (cur_bkt .eq. i_bkt) exit z_loop ! Done with cit buckets
                                             ! image i claims pairs from.
        end do x_loop
      end do y_loop
    end do z_loop

  end do

#else

! For the uniprocessor, just map it all...

  if (is_orthog_stk .ne. 0) then
    do img_idx = 1, atm_cnt
      atm_idx = img_atm_map(img_idx)
      img_crd(1, img_idx) = fraction(1, atm_idx) * x_box
      img_crd(2, img_idx) = fraction(2, atm_idx) * y_box
      img_crd(3, img_idx) = fraction(3, atm_idx) * z_box
      img_qterm(img_idx) = charge(atm_idx)
      img_iac(img_idx) = iac(atm_idx)
    end do
  else
    do img_idx = 1, atm_cnt
      atm_idx = img_atm_map(img_idx)
      f1 = fraction(1, atm_idx)
      f2 = fraction(2, atm_idx)
      f3 = fraction(3, atm_idx)
      ! ucell(2,1), ucell(3,1), and ucell(3,2) are always 0.d0, so
      ! we can simplify the expression in this critical inner loop
      img_crd(1, img_idx) = f1 * ucell_stk(1, 1) + &
                            f2 * ucell_stk(1, 2) + &
                            f3 * ucell_stk(1, 3)
      img_crd(2, img_idx) = f2 * ucell_stk(2, 2) + &
                            f3 * ucell_stk(2, 3)
      img_crd(3, img_idx) = f3 * ucell_stk(3, 3)
      img_qterm(img_idx) = charge(atm_idx)
      img_iac(img_idx) = iac(atm_idx)
    end do
  end if

#endif /* MPI */

  return

end subroutine map_pairlist_imgs
#ifdef MIC_offload
!*******************************************************************************
!
! Subroutine:   get_nb_list_offlaod
!
! Description:  Main routine for pairlist building with MIC offloading.  
!
!*******************************************************************************

subroutine get_nb_list_offload(atm_cnt, flat_cit, img_crd, img_atm_map, &
                       used_img_map, used_img_cnt, used_img_lst,ifail_dev, &
                       ico, fraction, tranvec, atm_maskdata, atm_mask, &
                       atm_img_map, excl_img_flags, &
                       img_iac, igroup, ntypes, ibelly, &
                       es_cutoff, vdw_cutoff, skinnb, verbose, ifail)

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod
  use ti_mod
  use omp_lib
! using pre-allocated mic data
  use offload_allocation_mod 

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)

  double precision      :: img_crd(3, atm_cnt)
  integer               :: img_atm_map(atm_cnt)
#ifdef MPI
  integer(byte)         :: used_img_map(atm_cnt)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
#endif
  integer               :: ico(*)
  double precision      :: fraction(3, atm_cnt)
  double precision      :: tranvec(1:3, 0:17)
  !type(listdata_rec)    :: atm_maskdata(*)
  type(listdata_rec)    :: atm_maskdata(atm_cnt)
  integer               :: atm_mask(*)
  integer               :: atm_img_map(atm_cnt)
  integer               :: excl_img_flags(atm_cnt)
  integer               :: img_iac(atm_cnt)
  integer               :: igroup(atm_cnt)
  integer               :: ntypes
  integer               :: ibelly
  double precision      :: es_cutoff
  double precision      :: vdw_cutoff
  double precision      :: skinnb
  integer               :: verbose
  integer               :: ifail

  integer               :: ifail_dev ! first defined in pme_direct


! Local variables:

!DIR$ OPTIONS /offload_attribute_target=mic

  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
  double precision      :: cutlist_sq
  double precision      :: x_i, y_i, z_i
  integer               :: atm_i, img_i
  integer               :: i
  integer               :: atm_mask_idx

  integer               :: i_bkt_x_idx, i_bkt_y_idx, i_bkt_z_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: i_bkt, saved_i_bkt_img_hi
  integer               :: img_i_lo, img_i_hi
  integer               :: i_bkt_lo, i_bkt_hi
  logical               :: common_tran
  logical               :: dont_skip_belly_pairs
  integer               :: iaci, ntypes_stk
  integer               :: ee_eval_cnt
  integer               :: full_eval_cnt
  integer               :: is_orthog_stk

! Variables use in the included files.  We would put this stuff in 
! contained subroutines, but it really whacks performance for itanium/ifort.

  double precision      :: dx, dy, dz
  integer               :: cur_bkt, yz_bkt, z_bkt
  integer               :: cur_tran, yz_tran, z_tran
  integer               :: img_j
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  double precision      :: r_sq
  integer               :: enc_img
  double precision      :: i_tranvec(1:3, 0:17)

! Variables used for double cutoffs:
  
  double precision      :: es_cut_sq
  logical               :: cutoffs_equal

  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: x_trans(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: y_trans(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
  integer               :: z_trans(0 : cit_tbl_z_dim * 2 - 1)

  integer               :: total_pairs  ! DBG

  integer               :: img_j_ee_eval(atm_cnt)
  integer               :: img_j_full_eval(atm_cnt)
 /*The following are for OpenMP in MIC offload region*/
#ifdef MPI
  integer 		:: j 

  integer		:: tmp_flat_cit_img_hi, hi_range ! rpalcing flat_cit()%img_hi use in the loop, produced race cond
  integer		:: my_thread_id
 ! integer		:: thread_ipairs(0:ipairs_maxsize) 
 ! integer		:: thread_num_packed(0:1)

 ! integer(byte)         :: used_img_map_thread(atm_cnt)
 ! integer		:: ipairs_maxsize_thread
  integer		:: ifail_thread(0:1)
!offload version
  integer		:: ipairs_maxsize_thread_dev

!offload related
  integer		::start_img_crd_dev,end_img_crd_dev
  integer		::start_img_crd_host,end_img_crd_host
  integer		::i_bkt_lo_dev,i_bkt_hi_dev
  integer		::i_bkt_lo_host,i_bkt_hi_host
  !double precision	::mic_fraction                          
  integer		:: cit_bkt_delta_dev
  integer		:: sig1 !for asynchronous offload
  integer		:: thread_offset
  double precision	:: rsqr
  
  integer		:: cit_tbl_x_dim_dev, cit_tbl_y_dim_dev, cit_tbl_z_dim_dev
!dir$ end options
  mytaskid_dev = mytaskid
!I have to use the following since they are originally static
	cit_tbl_x_dim_dev = cit_tbl_x_dim
	cit_tbl_y_dim_dev = cit_tbl_y_dim
	cit_tbl_z_dim_dev = cit_tbl_z_dim
  cit_bkt_delta_dev = cit_bkt_delta !cit_bkt_delta_dev to be inside offload
  

if(mytaskid .eq. numtasks/2 .or. mytaskid .eq. numtasks/2 - 1) then 
  start_img_crd_host = my_img_lo+ int(dble(my_img_hi - my_img_lo)*mic_fraction) 
  end_img_crd_host = my_img_hi ! end crd for host
  start_img_crd_dev = my_img_lo ! start crd for mic
  end_img_crd_dev = start_img_crd_host -1
else
  start_img_crd_host = my_img_lo 
  end_img_crd_host = my_img_hi 
end if

	ifail_thread(0) =0 ! one thread on host

        excl_img_flags = 0 ! arrays sets to zero
#endif
/*Done with the Openmp variables*/
! "bkt" refers to a bucket or cell in the flat cit table.
! "bkts" refers to the bucket mapping tables.
! "trans" refers to the translation mapping tables.

#ifdef MPI
  if (my_img_lo .gt. my_img_hi) return
#endif /* MPI */
  num_packed = 0        ! For "regular" pairlist
  num_packed_dev = 0        ! For "regular" pairlist in offload
  total_pairs = 0       ! DBG

  dont_skip_belly_pairs = ibelly .eq. 0 .and. ti_mode .eq. 0
  ntypes_stk = ntypes
  is_orthog_stk = is_orthog

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

  es_cut_sq = (es_cutoff + skinnb) * (es_cutoff + skinnb)
  cutoffs_equal = (es_cutoff .eq. vdw_cutoff)
  cutlist_sq = (vdw_cutoff + skinnb) * (vdw_cutoff + skinnb)

  ! Set up the bucket and translation index mapping arrays.  We really only need
  ! to do this each time if the cit table dimensions are changing, but there
  ! may be advantages to having only local variables anyway.

  call setup_cit_tbl_bkts(x_bkts, y_bkts, z_bkts, x_trans, y_trans, z_trans)

  ! Get the low and high flat cit bucket indexes for this task:

#ifdef MPI

/*the following 4 calls are offload related to get the boundaries for i_bkt_lo and i_bkt_hi for offload*/
if(mytaskid .eq. numtasks/2 .or. mytaskid .eq. numtasks/2 - 1) then !changed
  call get_flat_cit_idx(start_img_crd_dev, img_atm_map, fraction, i_bkt_lo_dev)
  call get_flat_cit_idx(end_img_crd_dev, img_atm_map, fraction, i_bkt_hi_dev)
end if
  call get_flat_cit_idx(start_img_crd_host, img_atm_map, fraction, i_bkt_lo_host)
  call get_flat_cit_idx(end_img_crd_host, img_atm_map, fraction, i_bkt_hi_host)
#else
  i_bkt_lo = 0
  i_bkt_hi = cit_tbl_x_dim * cit_tbl_y_dim * cit_tbl_z_dim - 1
#endif /* MPI */

  ! Loop over the atoms you own, finding their nonbonded partners...
  ! img_i is the atom for which you are finding the partner atoms, img_j

!introduced mic offload with openmp threading below, now it only works with MPI
!1st part on device

#if 1
!offload version
if(mytaskid .eq. numtasks/2 .or. mytaskid .eq. numtasks/2 - 1) then
!dir$ offload begin mandatory signal(sig1) target (mic:0)  &

 in(flat_cit(0:cit_tbl_x_dim*cit_tbl_y_dim*cit_tbl_z_dim-1) : alloc_if(.false.) free_if(.false.) into(flat_cit_dev(:)))&
 in(fraction(:,1:atm_cnt) : alloc_if(.false.) free_if(.false.) into(fraction_dev(:,:)))&
 in(x_bkts(0 : cit_tbl_x_dim * 3 - 1) : alloc_if(.false.) free_if(.false.) into(x_bkts_dev(:)))&
 in(x_trans(0 : cit_tbl_x_dim * 3 - 1) : alloc_if(.false.) free_if(.false.) into(x_trans_dev(:)))& 
 in(y_trans(0 : cit_tbl_y_dim * 3 - 1) : alloc_if(.false.) free_if(.false.) into(y_trans_dev(:)))&
 in(y_bkts(0 : cit_tbl_y_dim * 3 - 1) : alloc_if(.false.) free_if(.false.) into(y_bkts_dev(:)))&
 in(z_bkts(0 : cit_tbl_z_dim * 2 - 1) : alloc_if(.false.) free_if(.false.) into(z_bkts_dev(:)))&
 in(z_trans(0 : cit_tbl_z_dim * 2 - 1) : alloc_if(.false.) free_if(.false.) into(z_trans_dev(:)))&
 out(used_img_lst_thread_dev(1:atm_cnt) : alloc_if(.false.) free_if(.false.))&
 nocopy(thread_num_packed_dev(0:num_threads_offload) : alloc_if(.false.) free_if(.false.))&
 nocopy(ipairs_dev(0:ipairs_maxsize_dev) : alloc_if(.false.) free_if(.false.))& !it could be ipairs_maxsize/fraction
 nocopy(thread_ipairs_dev(0:ipairs_maxsize_dev) : alloc_if(.false.) free_if(.false.))& !it could be ipairs_maxsize/fraction
 in(img_iac(1:atm_cnt) : alloc_if(.false.) free_if(.false.) into(gbl_img_iac_dev(:)))& !can we pre in at offload_alloc
 nocopy(excl_img_flags_dev(1:atm_cnt) : alloc_if(.false.) free_if(.false.))&
 in(atm_img_map(1:atm_cnt) : alloc_if(.false.) free_if(.false.) into(atm_img_map_dev(:)))& !can we pre in at offload_alloc
 in(img_atm_map(1:atm_cnt) : alloc_if(.false.) free_if(.false.) into(img_atm_map_dev(:)))& !can we pre in at offload_alloc
 nocopy(atm_maskdata_dev(1:atm_cnt) : alloc_if(.false.) free_if(.false.))& ! already "in" in the offload allocation, amt_maskdata is originally allocated and built in nb_exclusion.F90
 nocopy(atm_mask_dev(1:local_next*local_next_mult_fac) : alloc_if(.false.) free_if(.false.)) & ! local_* is defined in the offload_alloc, and offload_alloc got it from prmtop_dat.F90 and nb_exclusion.F90
 in(img_crd(:,1:atm_cnt) : alloc_if(.false.) free_if(.false.) into(img_crd_dev(:,1:atm_cnt)))&
 in(tranvec(:,0:17) : alloc_if(.false.) free_if(.false.) into(tranvec_dev(:,0:17)))&
 nocopy(typ_ico_dev(1:ntypes*ntypes) : alloc_if(.false.) free_if(.false.))& 
 nocopy(ifail_thread_dev(0:num_threads_offload) : alloc_if(.false.) free_if(.false.))&
 nocopy(img_j_full_eval_dev(1:max_nb_per_atm) : alloc_if(.false.) free_if(.false.))&
 nocopy(img_j_ee_eval_dev(1:max_nb_per_atm) : alloc_if(.false.) free_if(.false.))&
 nocopy(i_tranvec_dev(:,0:17) : alloc_if(.false.) free_if(.false.))&
 nocopy(start_indexes_dev(1:(1+end_img_crd_dev-start_img_crd_dev)): alloc_if(.false.) free_if(.false.)) &
 nocopy(i,ipairs_maxsize_thread_dev,my_thread_id, full_eval_cnt, ee_eval_cnt,i_bkt, img_i_lo, img_i_hi, saved_i_bkt_img_hi, atm_i, i_bkt_x_idx, common_tran,i_bkt_y_idx, i_bkt_z_idx,x_bkts_lo,x_bkts_hi,y_bkts_lo,y_bkts_hi,z_bkts_lo,z_bkts_hi,img_i,iaci,atm_mask_idx,x_i, y_i, z_i, z_bkts_idx, z_bkt,  y_bkts_idx, yz_bkt,x_bkts_idx, cur_bkt, img_j,dx, dy, dz,z_tran,yz_tran,cur_tran, enc_img, r_sq, tmp_flat_cit_img_hi, hi_range,total_pairs,j)&
 in(atm_cnt,num_threads_offload, ntypes_stk, start_img_crd_dev,end_img_crd_dev,i_bkt_lo_dev, i_bkt_hi_dev,ipairs_maxsize_dev,scale_fac_x, scale_fac_y,scale_fac_z,cit_tbl_x_dim_dev,cit_tbl_y_dim_dev,cit_tbl_z_dim_dev,cit_bkt_delta_dev,cutoffs_equal,cutlist_sq,dont_skip_belly_pairs,es_cut_sq,thread_offset,rsqr,mytaskid_dev)&
 out(ifail_dev)! , num_packed_dev) ! &
!$omp parallel
  my_thread_id = omp_get_thread_num()
  thread_num_packed_dev(my_thread_id)=0 !pre_allocation
  ifail_thread_dev(my_thread_id) =0 !pre_allocation
  !$omp do
! !dir$ vector nontemporal 	
  do i = 1, atm_cnt
         excl_img_flags_dev(i) = 0 
         used_img_lst_thread_dev(i) = 0 !experimental
end do
!$omp end parallel

!$omp parallel do default(shared) private(i_bkt, my_thread_id, &
!$omp img_i_lo, img_i_hi, saved_i_bkt_img_hi, atm_i, &  
!$omp i_bkt_x_idx, common_tran, i_bkt_y_idx, i_bkt_z_idx, x_bkts_lo,x_bkts_hi, &
!$omp  y_bkts_lo,y_bkts_hi, z_bkts_lo,z_bkts_hi, &
!$omp img_i,full_eval_cnt, ee_eval_cnt, iaci, atm_mask_idx, i,j, &
!$omp x_i, y_i, z_i, z_bkts_idx, z_bkt,  y_bkts_idx, yz_bkt, x_bkts_idx, cur_bkt, img_j, &
!$omp dx, dy, dz, z_tran, yz_tran , &  
!$omp cur_tran, enc_img, r_sq, tmp_flat_cit_img_hi, hi_range, total_pairs, &
!$omp ipairs_maxsize_thread_dev,thread_offset,rsqr, & 
!$omp img_j_ee_eval_dev,i_tranvec_dev , img_j_full_eval_dev) & 
!$omp firstprivate(excl_img_flags_dev) !array 
  do i_bkt = i_bkt_lo_dev,i_bkt_hi_dev  
	/* OpenMP related*/
    	my_thread_id = omp_get_thread_num()
	ipairs_maxsize_thread_dev = ipairs_maxsize_dev * (1) / num_threads_offload ! max possible neigh. limit each thread
  	thread_offset = ipairs_maxsize_dev*my_thread_id/num_threads_offload
    img_i_lo = flat_cit_dev(i_bkt)%img_lo
    img_i_hi = flat_cit_dev(i_bkt)%img_hi

    if (img_i_lo .eq. 0) cycle
    !saved_i_bkt_img_hi = img_i_hi ! dont need since we are using tmp_flat_cit_img_hi
    if (img_i_lo .lt. start_img_crd_dev) img_i_lo = start_img_crd_dev
    if (img_i_hi .gt. end_img_crd_dev) img_i_hi = end_img_crd_dev

    ! Get the x,y,z 3d cit indexes to be able to get the bucket ranges
    ! to search:

    atm_i = img_atm_map_dev(img_i_lo)
    i_bkt_x_idx = int(fraction_dev(1, atm_i) * scale_fac_x)
    i_bkt_y_idx = int(fraction_dev(2, atm_i) * scale_fac_y)
    i_bkt_z_idx = int(fraction_dev(3, atm_i) * scale_fac_z)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = i_bkt_x_idx + cit_tbl_x_dim_dev - cit_bkt_delta_dev
    x_bkts_hi = i_bkt_x_idx + cit_tbl_x_dim_dev + cit_bkt_delta_dev
    y_bkts_lo = i_bkt_y_idx + cit_tbl_y_dim_dev - cit_bkt_delta_dev
    y_bkts_hi = i_bkt_y_idx + cit_tbl_y_dim_dev + cit_bkt_delta_dev
    z_bkts_lo = i_bkt_z_idx + cit_tbl_z_dim_dev - cit_bkt_delta_dev
    z_bkts_hi = i_bkt_z_idx + cit_tbl_z_dim_dev
                                        
    common_tran = x_trans_dev(x_bkts_lo) .eq. 1 .and. &
                  x_trans_dev(x_bkts_hi) .eq. 1 .and. &
                  y_trans_dev(y_bkts_lo) .eq. 3 .and. &
                  y_trans_dev(y_bkts_hi) .eq. 3 .and. &
                  z_trans_dev(z_bkts_lo) .eq. 9 .and. &
                  z_trans_dev(z_bkts_hi) .eq. 9

    do img_i = img_i_lo, img_i_hi
	tmp_flat_cit_img_hi = img_i - 1	
      
      ee_eval_cnt = 0
      full_eval_cnt = 0
      iaci = ntypes_stk * (gbl_img_iac_dev(img_i) - 1)

      atm_i = img_atm_map_dev(img_i)

      ! Mark excluded j images:

      atm_mask_idx = atm_maskdata_dev(atm_i)%offset
      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata_dev(atm_i)%cnt
        excl_img_flags_dev(atm_img_map_dev(atm_mask_dev(i))) = 1
      end do

      ! These are imaged coordinates:

      x_i = img_crd_dev(1, img_i)
      y_i = img_crd_dev(2, img_i)
      z_i = img_crd_dev(3, img_i)


      if (cutoffs_equal) then

        if (common_tran) then

          z_loop1: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts_dev(z_bkts_idx)
            y_loop1: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts_dev(y_bkts_idx)
              x_loop1: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts_dev(x_bkts_idx)
		/*The following conditional is added for OpenMP to work*/
		if (cur_bkt .eq. i_bkt) then
			hi_range = tmp_flat_cit_img_hi
		else
			hi_range = flat_cit_dev(cur_bkt)%img_hi
		end if
		!dir$ noprefetch excl_img_flags_dev,typ_ico_dev,gbl_img_iac_dev
		!dir$ ivdep
                do img_j = flat_cit_dev(cur_bkt)%img_lo, hi_range
		!do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd_dev(1, img_j) - x_i
                  dy = img_crd_dev(2, img_j) - y_i
                  dz = img_crd_dev(3, img_j) - z_i
		  rsqr = dx * dx + dy * dy + dz * dz
                  if (rsqr .lt. cutlist_sq) then
                    if (excl_img_flags_dev(img_j) .eq. 0) then
                      if (typ_ico_dev(iaci + gbl_img_iac_dev(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval_dev(ee_eval_cnt) = img_j
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval_dev(full_eval_cnt) = img_j
                      end if
			used_img_lst_thread_dev(img_j) = 1
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop1 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop1
            end do y_loop1
          end do z_loop1

        else    ! not common_tran

          !dir$ ivdep
          do i = 0, 17
            i_tranvec_dev(1, i) = tranvec_dev(1, i) - x_i 
            i_tranvec_dev(2, i) = tranvec_dev(2, i) - y_i 
            i_tranvec_dev(3, i) = tranvec_dev(3, i) - z_i 
          end do

          z_loop2: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts_dev(z_bkts_idx)
            z_tran = z_trans_dev(z_bkts_idx)
            y_loop2: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts_dev(y_bkts_idx)
              yz_tran = z_tran + y_trans_dev(y_bkts_idx)
              x_loop2: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts_dev(x_bkts_idx)
                cur_tran = x_trans_dev(x_bkts_idx) + yz_tran
		/*The following conditional is added for OpenMP to work*/
		if (cur_bkt .eq. i_bkt) then
			hi_range = tmp_flat_cit_img_hi
		else
			hi_range = flat_cit_dev(cur_bkt)%img_hi
		end if
		!dir$ noprefetch excl_img_flags_dev,typ_ico_dev,gbl_img_iac_dev
		!dir$ ivdep
                do img_j = flat_cit_dev(cur_bkt)%img_lo, hi_range
	!	do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd_dev(1, img_j) + i_tranvec_dev(1, cur_tran)
                  dy = img_crd_dev(2, img_j) + i_tranvec_dev(2, cur_tran)
                  dz = img_crd_dev(3, img_j) + i_tranvec_dev(3, cur_tran)
		  rsqr = dx * dx + dy * dy + dz * dz
                  if (rsqr .lt. cutlist_sq) then
                    if (excl_img_flags_dev(img_j) .eq. 0) then
                      enc_img = ior(img_j, ishft(cur_tran, 27))
                      if (typ_ico_dev(iaci + gbl_img_iac_dev(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval_dev(ee_eval_cnt) = enc_img
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval_dev(full_eval_cnt) = enc_img
                      end if
			used_img_lst_thread_dev(img_j) = 1
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop2 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop2
            end do y_loop2
          end do z_loop2

        end if
      else  ! cutoffs not equal

        if (common_tran) then

          z_loop3: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts_dev(z_bkts_idx)
            y_loop3: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts_dev(y_bkts_idx)
              x_loop3: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts_dev(x_bkts_idx)
		/*The following conditional is added for OpenMP to work*/
		if (cur_bkt .eq. i_bkt) then
			hi_range = tmp_flat_cit_img_hi
		else
			hi_range = flat_cit_dev(cur_bkt)%img_hi
                end if
                !dir$ noprefetch excl_img_flags_dev,typ_ico_dev,gbl_img_iac_dev
                !dir$ ivdep
                do img_j = flat_cit_dev(cur_bkt)%img_lo, hi_range
                !do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd_dev(1, img_j) - x_i
                  dy = img_crd_dev(2, img_j) - y_i
                  dz = img_crd_dev(3, img_j) - z_i
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags_dev(img_j) .eq. 0) then
                      if (typ_ico_dev(iaci + gbl_img_iac_dev(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval_dev(ee_eval_cnt) = img_j
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval_dev(full_eval_cnt) = img_j
                      end if
			used_img_lst_thread_dev(img_j) = 1
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop3 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop3
            end do y_loop3
          end do z_loop3

        else  ! not common_tran

          do i = 0, 17
            i_tranvec_dev(1, i) = tranvec_dev(1, i) - x_i 
            i_tranvec_dev(2, i) = tranvec_dev(2, i) - y_i 
            i_tranvec_dev(3, i) = tranvec_dev(3, i) - z_i 
          end do

          z_loop4: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts_dev(z_bkts_idx)
            z_tran = z_trans_dev(z_bkts_idx)
            y_loop4: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts_dev(y_bkts_idx)
              yz_tran = z_tran + y_trans_dev(y_bkts_idx)
              x_loop4: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts_dev(x_bkts_idx)
                cur_tran = x_trans_dev(x_bkts_idx) + yz_tran
		if (cur_bkt .eq. i_bkt) then
			hi_range = tmp_flat_cit_img_hi
		else
			hi_range = flat_cit_dev(cur_bkt)%img_hi
		end if
		!dir$ noprefetch excl_img_flags_dev,typ_ico_dev,gbl_img_iac_dev
		!dir$ ivdep
                do img_j = flat_cit_dev(cur_bkt)%img_lo, hi_range
		!do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd_dev(1, img_j) + i_tranvec_dev(1, cur_tran)
                  dy = img_crd_dev(2, img_j) + i_tranvec_dev(2, cur_tran)
                  dz = img_crd_dev(3, img_j) + i_tranvec_dev(3, cur_tran)
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags_dev(img_j) .eq. 0) then
                      if (typ_ico_dev(iaci + gbl_img_iac_dev(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval_dev(ee_eval_cnt) = &
                            ior(img_j,ishft(cur_tran, 27))
			used_img_lst_thread_dev(img_j) = 1
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval_dev(full_eval_cnt) = &
                          ior(img_j,ishft(cur_tran,27))
			used_img_lst_thread_dev(img_j) = 1
                      end if
                    end if
                  end if
                end do

                if (cur_bkt .eq. i_bkt) exit z_loop4 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop4
            end do y_loop4
          end do z_loop4
        end if
!cutoffs not equal is not supported

      end if
      !total_pairs = total_pairs + full_eval_cnt + ee_eval_cnt       ! DBG


      if (dont_skip_belly_pairs) then
#ifdef MPI
!dir$ forceinline
	 call pack_nb_list_offload(ee_eval_cnt, img_j_ee_eval_dev, &
                          full_eval_cnt, img_j_full_eval_dev, &
                          thread_ipairs_dev(thread_offset), &
			  thread_num_packed_dev(my_thread_id), common_tran,ipairs_maxsize_thread_dev,ifail_thread_dev, my_thread_id,num_threads_offload)
#else
	call pack_nb_list(ee_eval_cnt, img_j_ee_eval_dev, &
                          full_eval_cnt, img_j_full_eval_dev, &
                          gbl_sluipairs, num_packed_dev)       
#endif /* MPI */


!currently the pack_nb_list_skip_ti_pairs and pack_nb_list_skip_belly_pairs  subroutines are not supported for offload
#if 0
      else if (ti_mode .ne. 0) then
        call pack_nb_list_skip_ti_pairs(ee_eval_cnt, img_j_ee_eval_dev, &
                                        full_eval_cnt, img_j_full_eval_dev, &
                                        gbl_sluipairs, ti_lst, img_atm_map_dev, &
                                        num_packed_dev)
      else
        call pack_nb_list_skip_belly_pairs_offload(ee_eval_cnt, img_j_ee_eval_dev, &
                                           full_eval_cnt, img_j_full_eval_dev, &
                                           gbl_sluipairs, igroup, img_atm_map_dev, &
#ifdef MPI
					num_packed, ipairs_maxsize_thread_dev,ifail_thread_dev, my_thread_id,num_threads_offload)
#else
					num_packed_dev)
#endif


#endif
      end if

      ! Clear excluded j images flags:

      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata_dev(atm_i)%cnt
        excl_img_flags_dev(atm_img_map_dev(atm_mask_dev(i))) = 0
      end do
! we need to get to get rid of the following
#ifndef MPI ! the MPI version for handling the case ifail=1 (re)
      if (ifail .eq. 1) then
        flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi
        return
      end if
#endif
    end do ! end of images in i_bkt


#ifndef MPI
    flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi ! dont need this since I used  tmp_flat_cit_img_hi
#endif

  end do ! end of all cit buckets with i images owned by this task

	ifail_dev=0 ! its in value was 0 in pme_direct
!	!$omp simd Reduction(+:ifail_dev) 
	do i=0,num_threads_offload-1
		ifail_dev=ifail_dev + ifail_thread_dev(i)
	end do

!packing the thread parts of thread_ipairs_dev into ipairs_dev
#if 1
if(ifail_dev .eq. 0) then
  num_packed_dev = 0
  Do i=0, num_threads_offload-1 
  	thread_offset = ipairs_maxsize_dev*i/num_threads_offload
        ipairs_dev((num_packed_dev+1):(num_packed_dev+1)+thread_num_packed_dev(i)) = thread_ipairs_dev(thread_offset : thread_offset + thread_num_packed_dev(i)-1)
	num_packed_dev =num_packed_dev + thread_num_packed_dev(i)
  end do
end if
#endif
#if 0 /* experimental version */
  Do i=0, num_threads_offload-1 
  	thread_offset = ipairs_maxsize_dev*i/num_threads_offload
	!$omp parallel do schedule(dynamic,20)
	Do j= 1, thread_num_packed_dev(i)
		ipairs_dev(num_packed_dev + j) = thread_ipairs_dev(thread_offset+j-1)
	end do	
	num_packed_dev =num_packed_dev + thread_num_packed_dev(i)
  end do
#endif
/* done with the nb list building on device*/
#endif
!dir$ end offload

end if ! mytaskid
/*2nd part on host running parallely with MIC*/
  !do i_bkt = i_bkt_lo, i_bkt_hi 
  do i_bkt = i_bkt_lo_host, i_bkt_hi_host
    img_i_lo = flat_cit(i_bkt)%img_lo
    img_i_hi = flat_cit(i_bkt)%img_hi

    if (img_i_lo .eq. 0) cycle
    saved_i_bkt_img_hi = img_i_hi ! dont need for omp since we are using tmp_flat_cit_img_hi
    !if (img_i_lo .lt. my_img_lo) img_i_lo = my_img_lo
    !if (img_i_hi .gt. my_img_hi) img_i_hi = my_img_hi
    if (img_i_lo .lt. start_img_crd_host ) img_i_lo = start_img_crd_host
    if (img_i_hi .gt. end_img_crd_host ) img_i_hi = end_img_crd_host

    ! Get the x,y,z 3d cit indexes to be able to get the bucket ranges
    ! to search:

    atm_i = img_atm_map(img_i_lo)
    i_bkt_x_idx = int(fraction(1, atm_i) * scale_fac_x)
    i_bkt_y_idx = int(fraction(2, atm_i) * scale_fac_y)
    i_bkt_z_idx = int(fraction(3, atm_i) * scale_fac_z)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = i_bkt_x_idx + cit_tbl_x_dim - cit_bkt_delta
    x_bkts_hi = i_bkt_x_idx + cit_tbl_x_dim + cit_bkt_delta
    y_bkts_lo = i_bkt_y_idx + cit_tbl_y_dim - cit_bkt_delta
    y_bkts_hi = i_bkt_y_idx + cit_tbl_y_dim + cit_bkt_delta
    z_bkts_lo = i_bkt_z_idx + cit_tbl_z_dim - cit_bkt_delta
    z_bkts_hi = i_bkt_z_idx + cit_tbl_z_dim
                                        
    common_tran = x_trans(x_bkts_lo) .eq. 1 .and. &
                  x_trans(x_bkts_hi) .eq. 1 .and. &
                  y_trans(y_bkts_lo) .eq. 3 .and. &
                  y_trans(y_bkts_hi) .eq. 3 .and. &
                  z_trans(z_bkts_lo) .eq. 9 .and. &
                  z_trans(z_bkts_hi) .eq. 9

    do img_i = img_i_lo, img_i_hi
	flat_cit(i_bkt)%img_hi = img_i - 1
      
      ee_eval_cnt = 0
      full_eval_cnt = 0
      iaci = ntypes_stk * (img_iac(img_i) - 1)

      atm_i = img_atm_map(img_i)

      ! Mark excluded j images:

      atm_mask_idx = atm_maskdata(atm_i)%offset
      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 1
      end do

      ! These are imaged coordinates:

      x_i = img_crd(1, img_i)
      y_i = img_crd(2, img_i)
      z_i = img_crd(3, img_i)


      if (cutoffs_equal) then

        if (common_tran) then

          z_loop5: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop5: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop5: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
		do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = img_j
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
                      if (used_img_map(img_j) .eq. 0) then
			used_img_map(img_j) = 1
			used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop5 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop5
            end do y_loop5
          end do z_loop5

        else    ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop6: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop6: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop6: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
		do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      enc_img = ior(img_j, ishft(cur_tran, 27))
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = enc_img
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = enc_img
                      end if
                      if (used_img_map(img_j) .eq. 0) then
			used_img_map(img_j) = 1
			used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop6 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop6
            end do y_loop6
          end do z_loop6

        end if

      else  ! cutoffs not equal

        if (common_tran) then

          z_loop7: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop7: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop7: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
		do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = img_j
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop7 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop7
            end do y_loop7
          end do z_loop7

        else  ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop8: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop8: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop8: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
		do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = &
                            ior(img_j,ishft(cur_tran, 27))
                          if (used_img_map(img_j) .eq. 0) then
                            used_img_map(img_j) = 1
                            used_img_cnt = used_img_cnt + 1
                            used_img_lst(used_img_cnt) = img_j
                          end if
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = &
                          ior(img_j,ishft(cur_tran,27))
                        if (used_img_map(img_j) .eq. 0) then
                          used_img_map(img_j) = 1
                          used_img_cnt = used_img_cnt + 1
                          used_img_lst(used_img_cnt) = img_j
                        end if
                      end if
                    end if
                  end if
                end do

                if (cur_bkt .eq. i_bkt) exit z_loop8 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop8
            end do y_loop8
          end do z_loop8
        end if
      end if

      total_pairs = total_pairs + full_eval_cnt + ee_eval_cnt       ! DBG


      if (dont_skip_belly_pairs) then
        call pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                          full_eval_cnt, img_j_full_eval, &
                          gbl_sluipairs, num_packed)
      else if (ti_mode .ne. 0) then
        call pack_nb_list_skip_ti_pairs(ee_eval_cnt, img_j_ee_eval, &
                                        full_eval_cnt, img_j_full_eval, &
                                        gbl_sluipairs, img_atm_map, &
                                        num_packed)
      else
        call pack_nb_list_skip_belly_pairs(ee_eval_cnt, img_j_ee_eval, &
                                           full_eval_cnt, img_j_full_eval, &
                                           gbl_sluipairs, igroup, img_atm_map, &
                                           num_packed)
      end if

      ! Clear excluded j images flags:

      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 0
      end do
#ifndef MPI ! the MPI version for handling the case ifail=1 (re)
      if (ifail .eq. 1) then
        flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi
        return
      end if
#endif
    end do ! end of images in i_bkt


    flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi 


  end do ! end of all cit buckets with i images owned by this task

/*done with nb list building on host*/


#ifdef MPI 

if(mytaskid .eq. numtasks/2 .or. mytaskid .eq. numtasks/2 - 1) then
 !dir$ offload begin   mandatory target(mic:0) wait(sig1) !changed
 !dir$ end offload
end if ! mytaskid

      if (ifail .gt. 0 .or. ifail_dev .gt. 0) then
!        print *, "From nb_pairlist: nb list needs to grow ", "Host ifail", ifail,"Device ifail_dev",ifail_dev,"for mytaskid", mytaskid
        return ! to allocate larger space for gbl_ipairs
      end if
      ! need to have some code if ifail_dev .gt. 0, increase the array size ipairs_dev
#endif

#ifdef MPI ! now reducing the used_img_lst from all the threads, removing the duplicate img_j

if(mytaskid .eq. numtasks/2 .or. mytaskid .eq. numtasks/2 - 1) then 
  !reduce the dev version,for actual offload we dont need to reduce into gbl_ipairs
  !updata used_img_img_map and used_img_lst for mic nb list
	do j=1, atm_cnt
		if(used_img_lst_thread_dev(j) .ne. 0) then
			if(used_img_map(j) .eq. 0) then
				used_img_map(j)=1		
				used_img_cnt = used_img_cnt +1 
				used_img_lst(used_img_cnt)=j
			end if
		end if
	end do
end if

#endif

! write(0, *) 'DBG: total_pairs, listtot = ', &
!             total_pairs, num_packed
! write(0, *) 'DBG: avg packing = ', &
!             dble(total_pairs)/dble(num_packed)

  if (verbose .gt. 0) then
    if (master) write(mdout, *) 'listtot = ', num_packed
  end if

  return

contains

#include "pack_nb_list_dflt.i"
#ifdef MIC_offload
#include "pack_nb_list_dflt_offload.i"
#endif
end subroutine get_nb_list_offload
#endif /*MIC_offload*/
!*******************************************************************************
!
! Subroutine:   get_nb_list
!
! Description:  Main routine for pairlist building.  
!
!*******************************************************************************

subroutine get_nb_list(atm_cnt, flat_cit, img_crd, img_atm_map, &
#ifdef MPI
                       used_img_map, used_img_cnt, used_img_lst, &
#endif
                       ico, fraction, tranvec, atm_maskdata, atm_mask, &
                       atm_img_map, excl_img_flags, &
                       img_iac, igroup, ntypes, ibelly, &
                       es_cutoff, vdw_cutoff, skinnb, verbose, ifail)

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)

  double precision      :: img_crd(3, atm_cnt)
  integer               :: img_atm_map(atm_cnt)
#ifdef MPI
  integer(byte)         :: used_img_map(atm_cnt)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
#endif
  integer               :: ico(*)
  double precision      :: fraction(3, atm_cnt)
  double precision      :: tranvec(1:3, 0:17)
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
  integer               :: atm_img_map(atm_cnt)
  integer               :: excl_img_flags(atm_cnt)
  integer               :: img_iac(atm_cnt)
  integer               :: igroup(atm_cnt)
  integer               :: ntypes
  integer               :: ibelly
  double precision      :: es_cutoff
  double precision      :: vdw_cutoff
  double precision      :: skinnb
  integer               :: verbose
  integer               :: ifail

! Local variables:

  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
  double precision      :: cutlist_sq
  double precision      :: x_i, y_i, z_i
  integer               :: atm_i, img_i
  integer               :: i
  integer               :: atm_mask_idx
  integer               :: num_packed
  integer               :: i_bkt_x_idx, i_bkt_y_idx, i_bkt_z_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: i_bkt, saved_i_bkt_img_hi
  integer               :: img_i_lo, img_i_hi
  integer               :: i_bkt_lo, i_bkt_hi
  logical               :: common_tran
  logical               :: dont_skip_belly_pairs
  integer               :: iaci, ntypes_stk
  integer               :: ee_eval_cnt
  integer               :: full_eval_cnt
  integer               :: is_orthog_stk

! Variables use in the included files.  We would put this stuff in 
! contained subroutines, but it really whacks performance for itanium/ifort.

  double precision      :: dx, dy, dz
  integer               :: cur_bkt, yz_bkt, z_bkt
  integer               :: cur_tran, yz_tran, z_tran
  integer               :: img_j
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  double precision      :: r_sq
  integer               :: enc_img
  double precision      :: i_tranvec(1:3, 0:17)

! Variables used for double cutoffs:

  double precision      :: es_cut_sq
  logical               :: cutoffs_equal

  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: x_trans(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: y_trans(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
  integer               :: z_trans(0 : cit_tbl_z_dim * 2 - 1)

  integer               :: total_pairs  ! DBG

  integer               :: img_j_ee_eval(atm_cnt)
  integer               :: img_j_full_eval(atm_cnt)
  integer               :: atm_j

! "bkt" refers to a bucket or cell in the flat cit table.
! "bkts" refers to the bucket mapping tables.
! "trans" refers to the translation mapping tables.

#ifdef MPI
  if (my_img_lo .gt. my_img_hi) return
#endif /* MPI */
  num_packed = 0        ! For "regular" pairlist
  total_pairs = 0       ! DBG

  dont_skip_belly_pairs = ibelly .eq. 0 .and. ti_mode .eq. 0
  ntypes_stk = ntypes
  is_orthog_stk = is_orthog

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

  es_cut_sq = (es_cutoff + skinnb) * (es_cutoff + skinnb)
  cutoffs_equal = (es_cutoff .eq. vdw_cutoff)
  cutlist_sq = (vdw_cutoff + skinnb) * (vdw_cutoff + skinnb)

  ! Set up the bucket and translation index mapping arrays.  We really only need
  ! to do this each time if the cit table dimensions are changing, but there
  ! may be advantages to having only local variables anyway.

  call setup_cit_tbl_bkts(x_bkts, y_bkts, z_bkts, x_trans, y_trans, z_trans)

  ! Get the low and high flat cit bucket indexes for this task:

#ifdef MPI
  call get_flat_cit_idx(my_img_lo, img_atm_map, fraction, i_bkt_lo)
  call get_flat_cit_idx(my_img_hi, img_atm_map, fraction, i_bkt_hi)
#else
  i_bkt_lo = 0
  i_bkt_hi = cit_tbl_x_dim * cit_tbl_y_dim * cit_tbl_z_dim - 1
#endif /* MPI */

  ! Loop over the atoms you own, finding their nonbonded partners...
  ! img_i is the atom for which you are finding the partner atoms, img_j
  do i_bkt = i_bkt_lo, i_bkt_hi

    img_i_lo = flat_cit(i_bkt)%img_lo
    img_i_hi = flat_cit(i_bkt)%img_hi

    if (img_i_lo .eq. 0) cycle
    saved_i_bkt_img_hi = img_i_hi
#ifdef MPI
    if (img_i_lo .lt. my_img_lo) img_i_lo = my_img_lo
    if (img_i_hi .gt. my_img_hi) img_i_hi = my_img_hi
#endif /* MPI */

    ! Get the x,y,z 3d cit indexes to be able to get the bucket ranges
    ! to search:

    atm_i = img_atm_map(img_i_lo)
    i_bkt_x_idx = int(fraction(1, atm_i) * scale_fac_x)
    i_bkt_y_idx = int(fraction(2, atm_i) * scale_fac_y)
    i_bkt_z_idx = int(fraction(3, atm_i) * scale_fac_z)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = i_bkt_x_idx + cit_tbl_x_dim - cit_bkt_delta
    x_bkts_hi = i_bkt_x_idx + cit_tbl_x_dim + cit_bkt_delta
    y_bkts_lo = i_bkt_y_idx + cit_tbl_y_dim - cit_bkt_delta
    y_bkts_hi = i_bkt_y_idx + cit_tbl_y_dim + cit_bkt_delta
    z_bkts_lo = i_bkt_z_idx + cit_tbl_z_dim - cit_bkt_delta
    z_bkts_hi = i_bkt_z_idx + cit_tbl_z_dim
                                        
    common_tran = x_trans(x_bkts_lo) .eq. 1 .and. &
                  x_trans(x_bkts_hi) .eq. 1 .and. &
                  y_trans(y_bkts_lo) .eq. 3 .and. &
                  y_trans(y_bkts_hi) .eq. 3 .and. &
                  z_trans(z_bkts_lo) .eq. 9 .and. &
                  z_trans(z_bkts_hi) .eq. 9

    do img_i = img_i_lo, img_i_hi

      flat_cit(i_bkt)%img_hi = img_i - 1

      ee_eval_cnt = 0
      full_eval_cnt = 0
      iaci = ntypes_stk * (img_iac(img_i) - 1)

      atm_i = img_atm_map(img_i)

      ! Mark excluded j images:

      atm_mask_idx = atm_maskdata(atm_i)%offset
      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 1
      end do

      ! These are imaged coordinates:

      x_i = img_crd(1, img_i)
      y_i = img_crd(2, img_i)
      z_i = img_crd(3, img_i)

      if (cutoffs_equal) then

        if (common_tran) then

          z_loop1: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop1: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop1: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = img_j
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop1 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop1
            end do y_loop1
          end do z_loop1

        else    ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop2: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop2: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop2: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                 atm_j=img_atm_map(img_j)
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      enc_img = ior(img_j, ishft(cur_tran, 27))
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = enc_img
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = enc_img
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop2 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop2
            end do y_loop2
          end do z_loop2

        end if

      else  ! cutoffs not equal

        if (common_tran) then

          z_loop3: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop3: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop3: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = img_j
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop3 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop3
            end do y_loop3
          end do z_loop3

        else  ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop4: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop4: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop4: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = &
                            ior(img_j,ishft(cur_tran, 27))
#ifdef MPI
                          if (used_img_map(img_j) .eq. 0) then
                            used_img_map(img_j) = 1
                            used_img_cnt = used_img_cnt + 1
                            used_img_lst(used_img_cnt) = img_j
                          end if
#endif /* MPI */
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = &
                          ior(img_j,ishft(cur_tran,27))
#ifdef MPI
                        if (used_img_map(img_j) .eq. 0) then
                          used_img_map(img_j) = 1
                          used_img_cnt = used_img_cnt + 1
                          used_img_lst(used_img_cnt) = img_j
                        end if
#endif /* MPI */
                      end if
                    end if
                  end if
                end do

                if (cur_bkt .eq. i_bkt) exit z_loop4 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop4
            end do y_loop4
          end do z_loop4
        end if
      end if

      total_pairs = total_pairs + full_eval_cnt + ee_eval_cnt       ! DBG

      if (dont_skip_belly_pairs) then
        call pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                          full_eval_cnt, img_j_full_eval, &
                          gbl_sluipairs, num_packed)
      else if (ti_mode .ne. 0) then
        call pack_nb_list_skip_ti_pairs(ee_eval_cnt, img_j_ee_eval, &
                                        full_eval_cnt, img_j_full_eval, &
                                        gbl_sluipairs, img_atm_map, &
                                        num_packed)
      else
        call pack_nb_list_skip_belly_pairs(ee_eval_cnt, img_j_ee_eval, &
                                           full_eval_cnt, img_j_full_eval, &
                                           gbl_sluipairs, igroup, img_atm_map, &
                                           num_packed)
      end if

      ! Clear excluded j images flags:

      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 0
      end do

      if (ifail .eq. 1) then
        flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi
        return
      end if

    end do ! end of images in i_bkt

    flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi

  end do ! end of all cit buckets with i images owned by this task

! write(0, *) 'DBG: total_pairs, listtot = ', &
!             total_pairs, num_packed
! write(0, *) 'DBG: avg packing = ', &
!             dble(total_pairs)/dble(num_packed)

  if (verbose .gt. 0) then
    if (master) write(mdout, *) 'listtot = ', num_packed
  end if

  return

contains

#include "pack_nb_list_dflt.i"

end subroutine get_nb_list

!*******************************************************************************
!
! Subroutine:   get_nb_slulist
!
! Description:  Main routine for solute pairlist building. Added by Lijiang Yang 
!
!*******************************************************************************

subroutine get_nb_slulist(atm_cnt, bgsluatm, nsluatm, flat_cit, img_crd, img_atm_map, &
#ifdef MPI
                       used_img_map, used_img_cnt, used_img_lst, &
#endif
                       ico, fraction, tranvec, atm_maskdata, atm_mask, &
                       atm_img_map, excl_img_flags, &
                       img_iac, igroup, ntypes, ibelly, &
                       es_cutoff, vdw_cutoff, skinnb, verbose, ifail)

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
! Added by Lijiang Yang for solute pairlist
  integer               :: bgsluatm,nsluatm
! ------------------------------------------
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)

  double precision      :: img_crd(3, atm_cnt)
  integer               :: img_atm_map(atm_cnt)
#ifdef MPI
  integer(byte)         :: used_img_map(atm_cnt)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
#endif
  integer               :: ico(*)
  double precision      :: fraction(3, atm_cnt)
  double precision      :: tranvec(1:3, 0:17)
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
  integer               :: atm_img_map(atm_cnt)
  integer               :: excl_img_flags(atm_cnt)
  integer               :: img_iac(atm_cnt)
  integer               :: igroup(atm_cnt)
  integer               :: ntypes
  integer               :: ibelly
  double precision      :: es_cutoff
  double precision      :: vdw_cutoff
  double precision      :: skinnb
  integer               :: verbose
  integer               :: ifail

! Local variables:

  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
  double precision      :: cutlist_sq
  double precision      :: x_i, y_i, z_i
  integer               :: atm_i, img_i
  integer               :: atm_j !Added by Lijiang Yang
  integer               :: i
  integer               :: atm_mask_idx
  integer               :: num_packed
  integer               :: i_bkt_x_idx, i_bkt_y_idx, i_bkt_z_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: i_bkt, saved_i_bkt_img_hi
  integer               :: img_i_lo, img_i_hi
  integer               :: i_bkt_lo, i_bkt_hi
  logical               :: common_tran
  logical               :: dont_skip_belly_pairs
  integer               :: iaci, ntypes_stk
  integer               :: ee_eval_cnt
  integer               :: full_eval_cnt
  integer               :: is_orthog_stk

! Variables use in the included files.  We would put this stuff in 
! contained subroutines, but it really whacks performance for itanium/ifort.

  double precision      :: dx, dy, dz
  integer               :: cur_bkt, yz_bkt, z_bkt
  integer               :: cur_tran, yz_tran, z_tran
  integer               :: img_j
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  double precision      :: r_sq
  integer               :: enc_img
  double precision      :: i_tranvec(1:3, 0:17)

! Variables used for double cutoffs:

  double precision      :: es_cut_sq
  logical               :: cutoffs_equal

  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: x_trans(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: y_trans(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
  integer               :: z_trans(0 : cit_tbl_z_dim * 2 - 1)

  integer               :: total_pairs  ! DBG

  integer               :: img_j_ee_eval(atm_cnt)
  integer               :: img_j_full_eval(atm_cnt)
! "bkt" refers to a bucket or cell in the flat cit table.
! "bkts" refers to the bucket mapping tables.
! "trans" refers to the translation mapping tables.

#ifdef MPI
  if (my_img_lo .gt. my_img_hi) return
#endif /* MPI */
  num_packed = 0        ! For "regular" pairlist
  total_pairs = 0       ! DBG

  dont_skip_belly_pairs = ibelly .eq. 0 .and. ti_mode .eq. 0
  ntypes_stk = ntypes
  is_orthog_stk = is_orthog

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

  es_cut_sq = (es_cutoff + skinnb) * (es_cutoff + skinnb)
  cutoffs_equal = (es_cutoff .eq. vdw_cutoff)
  cutlist_sq = (vdw_cutoff + skinnb) * (vdw_cutoff + skinnb)

  ! Set up the bucket and translation index mapping arrays.  We really only need
  ! to do this each time if the cit table dimensions are changing, but there
  ! may be advantages to having only local variables anyway.

  call setup_cit_tbl_bkts(x_bkts, y_bkts, z_bkts, x_trans, y_trans, z_trans)

  ! Get the low and high flat cit bucket indexes for this task:

#ifdef MPI
  call get_flat_cit_idx(my_img_lo, img_atm_map, fraction, i_bkt_lo)
  call get_flat_cit_idx(my_img_hi, img_atm_map, fraction, i_bkt_hi)
#else
  i_bkt_lo = 0
  i_bkt_hi = cit_tbl_x_dim * cit_tbl_y_dim * cit_tbl_z_dim - 1
#endif /* MPI */

  ! Loop over the atoms you own, finding their nonbonded partners...
  ! img_i is the atom for which you are finding the partner atoms, img_j

  do i_bkt = i_bkt_lo, i_bkt_hi

    img_i_lo = flat_cit(i_bkt)%img_lo
    img_i_hi = flat_cit(i_bkt)%img_hi

    if (img_i_lo .eq. 0) cycle
    saved_i_bkt_img_hi = img_i_hi
#ifdef MPI
    if (img_i_lo .lt. my_img_lo) img_i_lo = my_img_lo
    if (img_i_hi .gt. my_img_hi) img_i_hi = my_img_hi
#endif /* MPI */

    ! Get the x,y,z 3d cit indexes to be able to get the bucket ranges
    ! to search:

    atm_i = img_atm_map(img_i_lo)
    i_bkt_x_idx = int(fraction(1, atm_i) * scale_fac_x)
    i_bkt_y_idx = int(fraction(2, atm_i) * scale_fac_y)
    i_bkt_z_idx = int(fraction(3, atm_i) * scale_fac_z)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = i_bkt_x_idx + cit_tbl_x_dim - cit_bkt_delta
    x_bkts_hi = i_bkt_x_idx + cit_tbl_x_dim + cit_bkt_delta
    y_bkts_lo = i_bkt_y_idx + cit_tbl_y_dim - cit_bkt_delta
    y_bkts_hi = i_bkt_y_idx + cit_tbl_y_dim + cit_bkt_delta
    z_bkts_lo = i_bkt_z_idx + cit_tbl_z_dim - cit_bkt_delta
    z_bkts_hi = i_bkt_z_idx + cit_tbl_z_dim
                                        
    common_tran = x_trans(x_bkts_lo) .eq. 1 .and. &
                  x_trans(x_bkts_hi) .eq. 1 .and. &
                  y_trans(y_bkts_lo) .eq. 3 .and. &
                  y_trans(y_bkts_hi) .eq. 3 .and. &
                  z_trans(z_bkts_lo) .eq. 9 .and. &
                  z_trans(z_bkts_hi) .eq. 9

    do img_i = img_i_lo, img_i_hi

      flat_cit(i_bkt)%img_hi = img_i - 1

      ee_eval_cnt = 0
      full_eval_cnt = 0
      iaci = ntypes_stk * (img_iac(img_i) - 1)

      atm_i = img_atm_map(img_i)
   
      if(atm_i .gt. nsluatm) cycle

      ! Mark excluded j images:

      atm_mask_idx = atm_maskdata(atm_i)%offset
      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 1
      end do

      ! These are imaged coordinates:

      x_i = img_crd(1, img_i)
      y_i = img_crd(2, img_i)
      z_i = img_crd(3, img_i)

      if (cutoffs_equal) then

        if (common_tran) then

          z_loop1: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop1: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop1: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                 atm_j=img_atm_map(img_j)
                 if(atm_j .ge. bgsluatm .and. atm_j .le. nsluatm) then !Added by Lijiang Yang 
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = img_j
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                 endif !atm_j
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop1 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop1
            end do y_loop1
          end do z_loop1

        else    ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop2: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop2: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop2: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                 atm_j=img_atm_map(img_j)
                 if(atm_j .ge. bgsluatm .and. atm_j .le. nsluatm) then !Added by Lijiang Yang 
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      enc_img = ior(img_j, ishft(cur_tran, 27))
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = enc_img
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = enc_img
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                 end if !atm_j
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop2 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop2
            end do y_loop2
          end do z_loop2

        end if

      else  ! cutoffs not equal

        if (common_tran) then

          z_loop3: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop3: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop3: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                 atm_j=img_atm_map(img_j)
                 if(atm_j .ge. bgsluatm .and. atm_j .le. nsluatm) then !Added by Lijiang Yang 
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = img_j
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                 endif
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop3 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop3
            end do y_loop3
          end do z_loop3

        else  ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop4: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop4: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop4: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                 atm_j=img_atm_map(img_j)
                 if(atm_j .ge. bgsluatm .and. atm_j .le. nsluatm) then !Added by Lijiang Yang 
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = &
                            ior(img_j,ishft(cur_tran, 27))
#ifdef MPI
                          if (used_img_map(img_j) .eq. 0) then
                            used_img_map(img_j) = 1
                            used_img_cnt = used_img_cnt + 1
                            used_img_lst(used_img_cnt) = img_j
                          end if
#endif /* MPI */
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = &
                          ior(img_j,ishft(cur_tran,27))
#ifdef MPI
                        if (used_img_map(img_j) .eq. 0) then
                          used_img_map(img_j) = 1
                          used_img_cnt = used_img_cnt + 1
                          used_img_lst(used_img_cnt) = img_j
                        end if
#endif /* MPI */
                      end if
                    end if
                  end if
                 end if !atm_j
                end do

                if (cur_bkt .eq. i_bkt) exit z_loop4 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop4
            end do y_loop4
          end do z_loop4
        end if
      end if

      total_pairs = total_pairs + full_eval_cnt + ee_eval_cnt       ! DBG

      if (dont_skip_belly_pairs) then
        call pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                          full_eval_cnt, img_j_full_eval, &
                          gbl_sluipairs, num_packed)
      else if (ti_mode .ne. 0) then
        call pack_nb_list_skip_ti_pairs(ee_eval_cnt, img_j_ee_eval, &
                                        full_eval_cnt, img_j_full_eval, &
                                        gbl_sluipairs, img_atm_map, &
                                        num_packed)
      else
        call pack_nb_list_skip_belly_pairs(ee_eval_cnt, img_j_ee_eval, &
                                           full_eval_cnt, img_j_full_eval, &
                                           gbl_sluipairs, igroup, img_atm_map, &
                                           num_packed)
      end if

      ! Clear excluded j images flags:

      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 0
      end do

      if (ifail .eq. 1) then
        flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi
        return
      end if

    end do ! end of images in i_bkt

    flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi

  end do ! end of all cit buckets with i images owned by this task

! write(0, *) 'DBG: total_pairs, listtot = ', &
!             total_pairs, num_packed
! write(0, *) 'DBG: avg packing = ', &
!             dble(total_pairs)/dble(num_packed)

  if (verbose .gt. 0) then
    if (master) write(mdout, *) 'listslu = ', num_packed
  end if

  return

contains

#include "pack_nb_list_dflt.i"

end subroutine get_nb_slulist

!*******************************************************************************
!
! Subroutine:   get_nb_sollist
!
! Description:  Main routine for solvent pairlist building. Added by Lijiang Yang 
!
!*******************************************************************************

subroutine get_nb_sollist(atm_cnt, bgsolatm, flat_cit, img_crd, img_atm_map, &
#ifdef MPI
                       used_img_map, used_img_cnt, used_img_lst, &
#endif
                       ico, fraction, tranvec, atm_maskdata, atm_mask, &
                       atm_img_map, excl_img_flags, &
                       img_iac, igroup, ntypes, ibelly, &
                       es_cutoff, vdw_cutoff, skinnb, verbose, ifail)

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
! Added by Lijiang Yang for solvent pairlist
  integer               :: bgsolatm
! ------------------------------------------
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)

  double precision      :: img_crd(3, atm_cnt)
  integer               :: img_atm_map(atm_cnt)
#ifdef MPI
  integer(byte)         :: used_img_map(atm_cnt)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
#endif
  integer               :: ico(*)
  double precision      :: fraction(3, atm_cnt)
  double precision      :: tranvec(1:3, 0:17)
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
  integer               :: atm_img_map(atm_cnt)
  integer               :: excl_img_flags(atm_cnt)
  integer               :: img_iac(atm_cnt)
  integer               :: igroup(atm_cnt)
  integer               :: ntypes
  integer               :: ibelly
  double precision      :: es_cutoff
  double precision      :: vdw_cutoff
  double precision      :: skinnb
  integer               :: verbose
  integer               :: ifail

! Local variables:

  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
  double precision      :: cutlist_sq
  double precision      :: x_i, y_i, z_i
  integer               :: atm_i, img_i
  integer               :: atm_j !Added by Lijiang Yang
  integer               :: i
  integer               :: atm_mask_idx
  integer               :: num_packed
  integer               :: i_bkt_x_idx, i_bkt_y_idx, i_bkt_z_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: i_bkt, saved_i_bkt_img_hi
  integer               :: img_i_lo, img_i_hi
  integer               :: i_bkt_lo, i_bkt_hi
  logical               :: common_tran
  logical               :: dont_skip_belly_pairs
  integer               :: iaci, ntypes_stk
  integer               :: ee_eval_cnt
  integer               :: full_eval_cnt
  integer               :: is_orthog_stk

! Variables use in the included files.  We would put this stuff in 
! contained subroutines, but it really whacks performance for itanium/ifort.

  double precision      :: dx, dy, dz
  integer               :: cur_bkt, yz_bkt, z_bkt
  integer               :: cur_tran, yz_tran, z_tran
  integer               :: img_j
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  double precision      :: r_sq
  integer               :: enc_img
  double precision      :: i_tranvec(1:3, 0:17)

! Variables used for double cutoffs:

  double precision      :: es_cut_sq
  logical               :: cutoffs_equal

  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: x_trans(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: y_trans(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
  integer               :: z_trans(0 : cit_tbl_z_dim * 2 - 1)

  integer               :: total_pairs  ! DBG

  integer               :: img_j_ee_eval(atm_cnt)
  integer               :: img_j_full_eval(atm_cnt)

! "bkt" refers to a bucket or cell in the flat cit table.
! "bkts" refers to the bucket mapping tables.
! "trans" refers to the translation mapping tables.

#ifdef MPI
  if (my_img_lo .gt. my_img_hi) return
#endif /* MPI */
  num_packed = 0        ! For "regular" pairlist
  total_pairs = 0       ! DBG

  dont_skip_belly_pairs = ibelly .eq. 0 .and. ti_mode .eq. 0
  ntypes_stk = ntypes
  is_orthog_stk = is_orthog

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

  es_cut_sq = (es_cutoff + skinnb) * (es_cutoff + skinnb)
  cutoffs_equal = (es_cutoff .eq. vdw_cutoff)
  cutlist_sq = (vdw_cutoff + skinnb) * (vdw_cutoff + skinnb)

  ! Set up the bucket and translation index mapping arrays.  We really only need
  ! to do this each time if the cit table dimensions are changing, but there
  ! may be advantages to having only local variables anyway.

  call setup_cit_tbl_bkts(x_bkts, y_bkts, z_bkts, x_trans, y_trans, z_trans)

  ! Get the low and high flat cit bucket indexes for this task:

#ifdef MPI
  call get_flat_cit_idx(my_img_lo, img_atm_map, fraction, i_bkt_lo)
  call get_flat_cit_idx(my_img_hi, img_atm_map, fraction, i_bkt_hi)
#else
  i_bkt_lo = 0
  i_bkt_hi = cit_tbl_x_dim * cit_tbl_y_dim * cit_tbl_z_dim - 1
#endif /* MPI */

  ! Loop over the atoms you own, finding their nonbonded partners...
  ! img_i is the atom for which you are finding the partner atoms, img_j

  do i_bkt = i_bkt_lo, i_bkt_hi

    img_i_lo = flat_cit(i_bkt)%img_lo
    img_i_hi = flat_cit(i_bkt)%img_hi

    if (img_i_lo .eq. 0) cycle
    saved_i_bkt_img_hi = img_i_hi
#ifdef MPI
    if (img_i_lo .lt. my_img_lo) img_i_lo = my_img_lo
    if (img_i_hi .gt. my_img_hi) img_i_hi = my_img_hi
#endif /* MPI */

    ! Get the x,y,z 3d cit indexes to be able to get the bucket ranges
    ! to search:

    atm_i = img_atm_map(img_i_lo)
    i_bkt_x_idx = int(fraction(1, atm_i) * scale_fac_x)
    i_bkt_y_idx = int(fraction(2, atm_i) * scale_fac_y)
    i_bkt_z_idx = int(fraction(3, atm_i) * scale_fac_z)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = i_bkt_x_idx + cit_tbl_x_dim - cit_bkt_delta
    x_bkts_hi = i_bkt_x_idx + cit_tbl_x_dim + cit_bkt_delta
    y_bkts_lo = i_bkt_y_idx + cit_tbl_y_dim - cit_bkt_delta
    y_bkts_hi = i_bkt_y_idx + cit_tbl_y_dim + cit_bkt_delta
    z_bkts_lo = i_bkt_z_idx + cit_tbl_z_dim - cit_bkt_delta
    z_bkts_hi = i_bkt_z_idx + cit_tbl_z_dim
                                        
    common_tran = x_trans(x_bkts_lo) .eq. 1 .and. &
                  x_trans(x_bkts_hi) .eq. 1 .and. &
                  y_trans(y_bkts_lo) .eq. 3 .and. &
                  y_trans(y_bkts_hi) .eq. 3 .and. &
                  z_trans(z_bkts_lo) .eq. 9 .and. &
                  z_trans(z_bkts_hi) .eq. 9

    do img_i = img_i_lo, img_i_hi

      flat_cit(i_bkt)%img_hi = img_i - 1

      ee_eval_cnt = 0
      full_eval_cnt = 0
      iaci = ntypes_stk * (img_iac(img_i) - 1)

      atm_i = img_atm_map(img_i)

      if(atm_i .lt. bgsolatm) cycle

      ! Mark excluded j images:

      atm_mask_idx = atm_maskdata(atm_i)%offset
      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 1
      end do

      ! These are imaged coordinates:

      x_i = img_crd(1, img_i)
      y_i = img_crd(2, img_i)
      z_i = img_crd(3, img_i)

      if (cutoffs_equal) then

        if (common_tran) then

          z_loop1: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop1: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop1: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                 atm_j=img_atm_map(img_j)
                 if(atm_j .ge. bgsolatm .and. atm_j .le. atm_cnt) then !Added by Lijiang Yang 
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = img_j
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                 endif !atm_j
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop1 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop1
            end do y_loop1
          end do z_loop1

        else    ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop2: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop2: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop2: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                 atm_j=img_atm_map(img_j)
                 if(atm_j .ge. bgsolatm .and. atm_j .le. atm_cnt) then !Added by Lijiang Yang 
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      enc_img = ior(img_j, ishft(cur_tran, 27))
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = enc_img
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = enc_img
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                 end if !atm_j
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop2 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop2
            end do y_loop2
          end do z_loop2

        end if

      else  ! cutoffs not equal

        if (common_tran) then

          z_loop3: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop3: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop3: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                 atm_j=img_atm_map(img_j)
                 if(atm_j .ge. bgsolatm .and. atm_j .le. atm_cnt) then !Added by Lijiang Yang 
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = img_j
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                 end if !atm_j
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop3 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop3
            end do y_loop3
          end do z_loop3

        else  ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop4: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop4: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop4: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                 atm_j=img_atm_map(img_j)
                 if(atm_j .ge. bgsolatm .and. atm_j .le. atm_cnt) then !Added by Lijiang Yang 
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = &
                            ior(img_j,ishft(cur_tran, 27))
#ifdef MPI
                          if (used_img_map(img_j) .eq. 0) then
                            used_img_map(img_j) = 1
                            used_img_cnt = used_img_cnt + 1
                            used_img_lst(used_img_cnt) = img_j
                          end if
#endif /* MPI */
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = &
                          ior(img_j,ishft(cur_tran,27))
#ifdef MPI
                        if (used_img_map(img_j) .eq. 0) then
                          used_img_map(img_j) = 1
                          used_img_cnt = used_img_cnt + 1
                          used_img_lst(used_img_cnt) = img_j
                        end if
#endif /* MPI */
                      end if
                    end if
                  end if
                 end if !atm_j
                end do

                if (cur_bkt .eq. i_bkt) exit z_loop4 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop4
            end do y_loop4
          end do z_loop4
        end if
      end if

      total_pairs = total_pairs + full_eval_cnt + ee_eval_cnt       ! DBG

      if (dont_skip_belly_pairs) then
        call pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                          full_eval_cnt, img_j_full_eval, &
                          gbl_solipairs, num_packed)
      else if (ti_mode .ne. 0) then
        call pack_nb_list_skip_ti_pairs(ee_eval_cnt, img_j_ee_eval, &
                                        full_eval_cnt, img_j_full_eval, &
                                        gbl_solipairs, img_atm_map, &
                                        num_packed)
      else
        call pack_nb_list_skip_belly_pairs(ee_eval_cnt, img_j_ee_eval, &
                                           full_eval_cnt, img_j_full_eval, &
                                           gbl_solipairs, igroup, img_atm_map, &
                                           num_packed)
      end if

      ! Clear excluded j images flags:

      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 0
      end do

      if (ifail .eq. 1) then
        flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi
        return
      end if

    end do ! end of images in i_bkt

    flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi

  end do ! end of all cit buckets with i images owned by this task

! write(0, *) 'DBG: total_pairs, listtot = ', &
!             total_pairs, num_packed
! write(0, *) 'DBG: avg packing = ', &
!             dble(total_pairs)/dble(num_packed)

  if (verbose .gt. 0) then
    if (master) write(mdout, *) 'listsol = ', num_packed
  end if

  return

contains

#include "pack_nb_list_dflt.i"

end subroutine get_nb_sollist

!*******************************************************************************
!
! Subroutine:   get_nb_interlist
!
! Description:  Main routine for interaction pairlist building. Added by Lijiang Yang 
!
!*******************************************************************************

subroutine get_nb_interlist(atm_cnt, bgsluatm, nsluatm, bgsolatm, flat_cit, img_crd, img_atm_map, &
#ifdef MPI
                       used_img_map, used_img_cnt, used_img_lst, &
#endif
                       ico, fraction, tranvec, atm_maskdata, atm_mask, &
                       atm_img_map, excl_img_flags, &
                       img_iac, igroup, ntypes, ibelly, &
                       es_cutoff, vdw_cutoff, skinnb, verbose, ifail)

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod
  use ti_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
! Added by Lijiang Yang for solute pairlist
  integer               :: bgsluatm,nsluatm,bgsolatm
! ------------------------------------------
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)

  double precision      :: img_crd(3, atm_cnt)
  integer               :: img_atm_map(atm_cnt)
#ifdef MPI
  integer(byte)         :: used_img_map(atm_cnt)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
#endif
  integer               :: ico(*)
  double precision      :: fraction(3, atm_cnt)
  double precision      :: tranvec(1:3, 0:17)
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
  integer               :: atm_img_map(atm_cnt)
  integer               :: excl_img_flags(atm_cnt)
  integer               :: img_iac(atm_cnt)
  integer               :: igroup(atm_cnt)
  integer               :: ntypes
  integer               :: ibelly
  double precision      :: es_cutoff
  double precision      :: vdw_cutoff
  double precision      :: skinnb
  integer               :: verbose
  integer               :: ifail

! Local variables:

  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
  double precision      :: cutlist_sq
  double precision      :: x_i, y_i, z_i
  integer               :: atm_i, img_i
  integer               :: atm_j !Added by Lijiang Yang
  integer               :: i
  integer               :: atm_mask_idx
  integer               :: num_packed
  integer               :: i_bkt_x_idx, i_bkt_y_idx, i_bkt_z_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: i_bkt, saved_i_bkt_img_hi
  integer               :: img_i_lo, img_i_hi
  integer               :: i_bkt_lo, i_bkt_hi
  logical               :: common_tran
  logical               :: dont_skip_belly_pairs
  integer               :: iaci, ntypes_stk
  integer               :: ee_eval_cnt
  integer               :: full_eval_cnt
  integer               :: is_orthog_stk

! Variables use in the included files.  We would put this stuff in 
! contained subroutines, but it really whacks performance for itanium/ifort.

  double precision      :: dx, dy, dz
  integer               :: cur_bkt, yz_bkt, z_bkt
  integer               :: cur_tran, yz_tran, z_tran
  integer               :: img_j
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  double precision      :: r_sq
  integer               :: enc_img
  double precision      :: i_tranvec(1:3, 0:17)

! Variables used for double cutoffs:

  double precision      :: es_cut_sq
  logical               :: cutoffs_equal

  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: x_trans(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: y_trans(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
  integer               :: z_trans(0 : cit_tbl_z_dim * 2 - 1)

  integer               :: total_pairs  ! DBG

  integer               :: img_j_ee_eval(atm_cnt)
  integer               :: img_j_full_eval(atm_cnt)

! "bkt" refers to a bucket or cell in the flat cit table.
! "bkts" refers to the bucket mapping tables.
! "trans" refers to the translation mapping tables.

#ifdef MPI
  if (my_img_lo .gt. my_img_hi) return
#endif /* MPI */
  num_packed = 0        ! For "regular" pairlist
  total_pairs = 0       ! DBG

  dont_skip_belly_pairs = ibelly .eq. 0 .and. ti_mode .eq. 0
  ntypes_stk = ntypes
  is_orthog_stk = is_orthog

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

  es_cut_sq = (es_cutoff + skinnb) * (es_cutoff + skinnb)
  cutoffs_equal = (es_cutoff .eq. vdw_cutoff)
  cutlist_sq = (vdw_cutoff + skinnb) * (vdw_cutoff + skinnb)

  ! Set up the bucket and translation index mapping arrays.  We really only need
  ! to do this each time if the cit table dimensions are changing, but there
  ! may be advantages to having only local variables anyway.

  call setup_cit_tbl_bkts(x_bkts, y_bkts, z_bkts, x_trans, y_trans, z_trans)

  ! Get the low and high flat cit bucket indexes for this task:

#ifdef MPI
  call get_flat_cit_idx(my_img_lo, img_atm_map, fraction, i_bkt_lo)
  call get_flat_cit_idx(my_img_hi, img_atm_map, fraction, i_bkt_hi)
#else
  i_bkt_lo = 0
  i_bkt_hi = cit_tbl_x_dim * cit_tbl_y_dim * cit_tbl_z_dim - 1
#endif /* MPI */

  ! Loop over the atoms you own, finding their nonbonded partners...
  ! img_i is the atom for which you are finding the partner atoms, img_j

  do i_bkt = i_bkt_lo, i_bkt_hi

    img_i_lo = flat_cit(i_bkt)%img_lo
    img_i_hi = flat_cit(i_bkt)%img_hi

    if (img_i_lo .eq. 0) cycle
    saved_i_bkt_img_hi = img_i_hi
#ifdef MPI
    if (img_i_lo .lt. my_img_lo) img_i_lo = my_img_lo
    if (img_i_hi .gt. my_img_hi) img_i_hi = my_img_hi
#endif /* MPI */

    ! Get the x,y,z 3d cit indexes to be able to get the bucket ranges
    ! to search:

    atm_i = img_atm_map(img_i_lo)
    i_bkt_x_idx = int(fraction(1, atm_i) * scale_fac_x)
    i_bkt_y_idx = int(fraction(2, atm_i) * scale_fac_y)
    i_bkt_z_idx = int(fraction(3, atm_i) * scale_fac_z)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = i_bkt_x_idx + cit_tbl_x_dim - cit_bkt_delta
    x_bkts_hi = i_bkt_x_idx + cit_tbl_x_dim + cit_bkt_delta
    y_bkts_lo = i_bkt_y_idx + cit_tbl_y_dim - cit_bkt_delta
    y_bkts_hi = i_bkt_y_idx + cit_tbl_y_dim + cit_bkt_delta
    z_bkts_lo = i_bkt_z_idx + cit_tbl_z_dim - cit_bkt_delta
    z_bkts_hi = i_bkt_z_idx + cit_tbl_z_dim
                                        
    common_tran = x_trans(x_bkts_lo) .eq. 1 .and. &
                  x_trans(x_bkts_hi) .eq. 1 .and. &
                  y_trans(y_bkts_lo) .eq. 3 .and. &
                  y_trans(y_bkts_hi) .eq. 3 .and. &
                  z_trans(z_bkts_lo) .eq. 9 .and. &
                  z_trans(z_bkts_hi) .eq. 9

    do img_i = img_i_lo, img_i_hi

      flat_cit(i_bkt)%img_hi = img_i - 1

      ee_eval_cnt = 0
      full_eval_cnt = 0
      iaci = ntypes_stk * (img_iac(img_i) - 1)

      atm_i = img_atm_map(img_i)
!      if(atm_i .gt. nsluatm) cycle

      ! Mark excluded j images:

      atm_mask_idx = atm_maskdata(atm_i)%offset
      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 1
      end do

      ! These are imaged coordinates:

      x_i = img_crd(1, img_i)
      y_i = img_crd(2, img_i)
      z_i = img_crd(3, img_i)

      if (cutoffs_equal) then

        if (common_tran) then

          z_loop1: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop1: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop1: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                 atm_j=img_atm_map(img_j)
                 if( (atm_i .ge. bgsluatm .and. atm_i .le. nsluatm .and. atm_j &
                    .ge. bgsolatm .and. atm_j .le. atm_cnt) .or. &
                     (atm_i .ge. bgsolatm .and. atm_i .le. atm_cnt .and. atm_j &
                      .ge. bgsluatm .and. atm_j .le. nsluatm) ) then !Added by Lijiang Yang 
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = img_j
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                 endif !atm_i,atm_j
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop1 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop1
            end do y_loop1
          end do z_loop1

        else    ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop2: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop2: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop2: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                 atm_j=img_atm_map(img_j)
                 if( (atm_i .ge. bgsluatm .and. atm_i .le. nsluatm .and. atm_j &
                    .ge. bgsolatm .and. atm_j .le. atm_cnt) .or. &
                     (atm_i .ge. bgsolatm .and. atm_i .le. atm_cnt .and. atm_j &
                      .ge. bgsluatm .and. atm_j .le. nsluatm) ) then !Added by Lijiang Yang 
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      enc_img = ior(img_j, ishft(cur_tran, 27))
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = enc_img
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = enc_img
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                 end if !atm_i,atm_j
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop2 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop2
            end do y_loop2
          end do z_loop2

        end if

      else  ! cutoffs not equal

        if (common_tran) then

          z_loop3: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop3: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop3: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                 atm_j=img_atm_map(img_j)
                 if( (atm_i .ge. bgsluatm .and. atm_i .le. nsluatm .and. atm_j &
                    .ge. bgsolatm .and. atm_j .le. atm_cnt) .or. &
                     (atm_i .ge. bgsolatm .and. atm_i .le. atm_cnt .and. atm_j &
                      .ge. bgsluatm .and. atm_j .le. nsluatm) ) then !Added by Lijiang Yang 
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = img_j
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                 endif
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop3 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop3
            end do y_loop3
          end do z_loop3

        else  ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop4: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop4: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop4: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                 atm_j=img_atm_map(img_j)
                 if( (atm_i .ge. bgsluatm .and. atm_i .le. nsluatm .and. atm_j &
                    .ge. bgsolatm .and. atm_j .le. atm_cnt) .or. &
                     (atm_i .ge. bgsolatm .and. atm_i .le. atm_cnt .and. atm_j &
                      .ge. bgsluatm .and. atm_j .le. nsluatm) ) then !Added by Lijiang Yang 
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = &
                            ior(img_j,ishft(cur_tran, 27))
#ifdef MPI
                          if (used_img_map(img_j) .eq. 0) then
                            used_img_map(img_j) = 1
                            used_img_cnt = used_img_cnt + 1
                            used_img_lst(used_img_cnt) = img_j
                          end if
#endif /* MPI */
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = &
                          ior(img_j,ishft(cur_tran,27))
#ifdef MPI
                        if (used_img_map(img_j) .eq. 0) then
                          used_img_map(img_j) = 1
                          used_img_cnt = used_img_cnt + 1
                          used_img_lst(used_img_cnt) = img_j
                        end if
#endif /* MPI */
                      end if
                    end if
                  end if
                 end if !atm_i,atm_j
                end do

                if (cur_bkt .eq. i_bkt) exit z_loop4 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop4
            end do y_loop4
          end do z_loop4
        end if
      end if

      total_pairs = total_pairs + full_eval_cnt + ee_eval_cnt       ! DBG

      if (dont_skip_belly_pairs) then
        call pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                          full_eval_cnt, img_j_full_eval, &
                          gbl_interipairs, num_packed)
      else if (ti_mode .ne. 0) then
        call pack_nb_list_skip_ti_pairs(ee_eval_cnt, img_j_ee_eval, &
                                        full_eval_cnt, img_j_full_eval, &
                                        gbl_interipairs, img_atm_map, &
                                        num_packed)
      else
        call pack_nb_list_skip_belly_pairs(ee_eval_cnt, img_j_ee_eval, &
                                           full_eval_cnt, img_j_full_eval, &
                                           gbl_interipairs, igroup, img_atm_map, &
                                           num_packed)
      end if

      ! Clear excluded j images flags:

      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 0
      end do

      if (ifail .eq. 1) then
        flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi
        return
      end if

    end do ! end of images in i_bkt

    flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi

  end do ! end of all cit buckets with i images owned by this task

! write(0, *) 'DBG: total_pairs, listtot = ', &
!             total_pairs, num_packed
! write(0, *) 'DBG: avg packing = ', &
!             dble(total_pairs)/dble(num_packed)

  if (verbose .gt. 0) then
    if (master) write(mdout, *) 'listinter = ', num_packed
  end if

  return

contains

#include "pack_nb_list_dflt.i"

end subroutine get_nb_interlist

!*******************************************************************************
!
! Subroutine:   get_nb_ips_list
!
! Description:  Main routine for pairlist building.  
!
!*******************************************************************************

subroutine get_nb_ips_list(atm_cnt, flat_cit, img_crd, img_atm_map, &
#ifdef MPI
                       used_img_map, used_img_cnt, used_img_lst, &
                       used_img_ips_map, used_img_ips_cnt, used_img_ips_lst, &
#endif
                       ico, fraction, tranvec, atm_maskdata, atm_mask, &
                       atm_img_map, excl_img_flags, &
                       img_iac, igroup, ntypes, ibelly, &
                       es_cutoff, vdw_cutoff, skinnb, verbose, ifail)

  use cit_mod
  use gbl_datatypes_mod
  use img_mod
  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  type(cit_tbl_rec)     :: flat_cit(0 : cit_tbl_x_dim * &
                                        cit_tbl_y_dim * &
                                        cit_tbl_z_dim - 1)

  double precision      :: img_crd(3, atm_cnt)
  integer               :: img_atm_map(atm_cnt)
#ifdef MPI
  integer(byte)         :: used_img_map(atm_cnt)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
  integer(byte)         :: used_img_ips_map(atm_cnt)
  integer               :: used_img_ips_cnt
  integer               :: used_img_ips_lst(*)
#endif
  integer               :: ico(*)
  double precision      :: fraction(3, atm_cnt)
  double precision      :: tranvec(1:3, 0:17)
  type(listdata_rec)    :: atm_maskdata(*)
  integer               :: atm_mask(*)
  integer               :: atm_img_map(atm_cnt)
  integer               :: excl_img_flags(atm_cnt)
  integer               :: img_iac(atm_cnt)
  integer               :: igroup(atm_cnt)
  integer               :: ntypes
  integer               :: ibelly
  double precision      :: es_cutoff
  double precision      :: vdw_cutoff
  double precision      :: skinnb
  integer               :: verbose
  integer               :: ifail

! Local variables:

  double precision      :: scale_fac_x, scale_fac_y, scale_fac_z
  double precision      :: cutlist_sq
  double precision      :: x_i, y_i, z_i
  integer               :: atm_i, img_i
  integer               :: i
  integer               :: atm_mask_idx
  integer               :: num_packed
  integer               :: i_bkt_x_idx, i_bkt_y_idx, i_bkt_z_idx
  integer               :: x_bkts_lo, x_bkts_hi
  integer               :: y_bkts_lo, y_bkts_hi
  integer               :: z_bkts_lo, z_bkts_hi
  integer               :: i_bkt, saved_i_bkt_img_hi
  integer               :: img_i_lo, img_i_hi
  integer               :: i_bkt_lo, i_bkt_hi
  logical               :: common_tran
  logical               :: dont_skip_belly_pairs
  integer               :: iaci, ntypes_stk
  integer               :: ee_eval_cnt
  integer               :: full_eval_cnt
  integer               :: is_orthog_stk

! Variables use in the included files.  We would put this stuff in 
! contained subroutines, but it really whacks performance for itanium/ifort.

  double precision      :: dx, dy, dz
  integer               :: cur_bkt, yz_bkt, z_bkt
  integer               :: cur_tran, yz_tran, z_tran
  integer               :: img_j
  integer               :: x_bkts_idx, y_bkts_idx, z_bkts_idx
  double precision      :: r_sq
  integer               :: enc_img
  double precision      :: i_tranvec(1:3, 0:17)

! Variables used for double cutoffs:

  double precision      :: es_cut_sq
  logical               :: cutoffs_equal

  integer               :: x_bkts(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: x_trans(0 : cit_tbl_x_dim * 3 - 1)
  integer               :: y_bkts(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: y_trans(0 : cit_tbl_y_dim * 3 - 1)
  integer               :: z_bkts(0 : cit_tbl_z_dim * 2 - 1)
  integer               :: z_trans(0 : cit_tbl_z_dim * 2 - 1)

  integer               :: total_pairs  ! DBG

  integer               :: img_j_ee_eval(atm_cnt)
  integer               :: img_j_full_eval(atm_cnt)

! "bkt" refers to a bucket or cell in the flat cit table.
! "bkts" refers to the bucket mapping tables.
! "trans" refers to the translation mapping tables.

#ifdef MPI
  if (my_img_lo .gt. my_img_hi) return
#endif /* MPI */
  num_packed = 0        ! For "regular" pairlist
  total_pairs = 0       ! DBG

  dont_skip_belly_pairs = ibelly .eq. 0
  ntypes_stk = ntypes
  is_orthog_stk = is_orthog

  scale_fac_x = dble(cit_tbl_x_dim)
  scale_fac_y = dble(cit_tbl_y_dim)
  scale_fac_z = dble(cit_tbl_z_dim)

  es_cut_sq = (es_cutoff + skinnb) * (es_cutoff + skinnb)
  cutoffs_equal = (es_cutoff .eq. vdw_cutoff)
  cutlist_sq = (vdw_cutoff + skinnb) * (vdw_cutoff + skinnb)

  ! Set up the bucket and translation index mapping arrays.  We really only need
  ! to do this each time if the cit table dimensions are changing, but there
  ! may be advantages to having only local variables anyway.

  call setup_cit_tbl_bkts(x_bkts, y_bkts, z_bkts, x_trans, y_trans, z_trans)

  ! Get the low and high flat cit bucket indexes for this task:

#ifdef MPI
  call get_flat_cit_idx(my_img_lo, img_atm_map, fraction, i_bkt_lo)
  call get_flat_cit_idx(my_img_hi, img_atm_map, fraction, i_bkt_hi)
#else
  i_bkt_lo = 0
  i_bkt_hi = cit_tbl_x_dim * cit_tbl_y_dim * cit_tbl_z_dim - 1
#endif /* MPI */

  ! Loop over the atoms you own, finding their nonbonded partners...
  ! img_i is the atom for which you are finding the partner atoms, img_j

  do i_bkt = i_bkt_lo, i_bkt_hi

    img_i_lo = flat_cit(i_bkt)%img_lo
    img_i_hi = flat_cit(i_bkt)%img_hi

    if (img_i_lo .eq. 0) cycle
    saved_i_bkt_img_hi = img_i_hi
#ifdef MPI
    if (img_i_lo .lt. my_img_lo) img_i_lo = my_img_lo
    if (img_i_hi .gt. my_img_hi) img_i_hi = my_img_hi
#endif /* MPI */

    ! Get the x,y,z 3d cit indexes to be able to get the bucket ranges
    ! to search:

    atm_i = img_atm_map(img_i_lo)
    i_bkt_x_idx = int(fraction(1, atm_i) * scale_fac_x)
    i_bkt_y_idx = int(fraction(2, atm_i) * scale_fac_y)
    i_bkt_z_idx = int(fraction(3, atm_i) * scale_fac_z)

    ! Determine the bucket ranges that need to be searched:

    x_bkts_lo = i_bkt_x_idx + cit_tbl_x_dim - cit_bkt_delta
    x_bkts_hi = i_bkt_x_idx + cit_tbl_x_dim + cit_bkt_delta
    y_bkts_lo = i_bkt_y_idx + cit_tbl_y_dim - cit_bkt_delta
    y_bkts_hi = i_bkt_y_idx + cit_tbl_y_dim + cit_bkt_delta
    z_bkts_lo = i_bkt_z_idx + cit_tbl_z_dim - cit_bkt_delta
    z_bkts_hi = i_bkt_z_idx + cit_tbl_z_dim
                                        
    common_tran = x_trans(x_bkts_lo) .eq. 1 .and. &
                  x_trans(x_bkts_hi) .eq. 1 .and. &
                  y_trans(y_bkts_lo) .eq. 3 .and. &
                  y_trans(y_bkts_hi) .eq. 3 .and. &
                  z_trans(z_bkts_lo) .eq. 9 .and. &
                  z_trans(z_bkts_hi) .eq. 9

    do img_i = img_i_lo, img_i_hi

      flat_cit(i_bkt)%img_hi = img_i - 1

      ee_eval_cnt = 0
      full_eval_cnt = 0
      iaci = ntypes_stk * (img_iac(img_i) - 1)

      atm_i = img_atm_map(img_i)

      ! Mark excluded j images:

      atm_mask_idx = atm_maskdata(atm_i)%offset
      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 1
      end do

      ! These are imaged coordinates:

      x_i = img_crd(1, img_i)
      y_i = img_crd(2, img_i)
      z_i = img_crd(3, img_i)

      if (cutoffs_equal) then

        if (common_tran) then

          z_loop1: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop1: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop1: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = img_j
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
                    else
                      if (used_img_ips_map(img_j) .eq. 0) then
                        used_img_ips_map(img_j) = 1
                        used_img_ips_cnt = used_img_ips_cnt + 1
                        used_img_ips_lst(used_img_ips_cnt) = img_j                      
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop1 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop1
            end do y_loop1
          end do z_loop1

        else    ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop2: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop2: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop2: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  if (dx * dx + dy * dy + dz * dz .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      enc_img = ior(img_j, ishft(cur_tran, 27))
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        ee_eval_cnt = ee_eval_cnt + 1
                        img_j_ee_eval(ee_eval_cnt) = enc_img
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = enc_img
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
                    else 
                      if (used_img_ips_map(img_j) .eq. 0) then
                        used_img_ips_map(img_j) = 1
                        used_img_ips_cnt = used_img_ips_cnt + 1
                        used_img_ips_lst(used_img_ips_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop2 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop2
            end do y_loop2
          end do z_loop2

        end if

      else  ! cutoffs not equal

        if (common_tran) then

          z_loop3: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            y_loop3: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              x_loop3: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) - x_i
                  dy = img_crd(2, img_j) - y_i
                  dz = img_crd(3, img_j) - z_i
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = img_j
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = img_j
                      end if
#ifdef MPI
                      if (used_img_map(img_j) .eq. 0) then
                        used_img_map(img_j) = 1
                        used_img_cnt = used_img_cnt + 1
                        used_img_lst(used_img_cnt) = img_j
                      end if
                    else
                      if (used_img_ips_map(img_j) .eq. 0) then
                        used_img_ips_map(img_j) = 1
                        used_img_ips_cnt = used_img_ips_cnt + 1
                        used_img_ips_lst(used_img_ips_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do
                if (cur_bkt .eq. i_bkt) exit z_loop3 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop3
            end do y_loop3
          end do z_loop3

        else  ! not common_tran

          do i = 0, 17
            i_tranvec(1, i) = tranvec(1, i) - x_i 
            i_tranvec(2, i) = tranvec(2, i) - y_i 
            i_tranvec(3, i) = tranvec(3, i) - z_i 
          end do

          z_loop4: &
          do z_bkts_idx = z_bkts_lo, z_bkts_hi
            z_bkt = z_bkts(z_bkts_idx)
            z_tran = z_trans(z_bkts_idx)
            y_loop4: &
            do y_bkts_idx = y_bkts_lo, y_bkts_hi
              yz_bkt = z_bkt + y_bkts(y_bkts_idx)
              yz_tran = z_tran + y_trans(y_bkts_idx)
              x_loop4: &
              do x_bkts_idx = x_bkts_lo, x_bkts_hi
                cur_bkt = yz_bkt + x_bkts(x_bkts_idx)
                cur_tran = x_trans(x_bkts_idx) + yz_tran
                do img_j = flat_cit(cur_bkt)%img_lo, flat_cit(cur_bkt)%img_hi
                  dx = img_crd(1, img_j) + i_tranvec(1, cur_tran)
                  dy = img_crd(2, img_j) + i_tranvec(2, cur_tran)
                  dz = img_crd(3, img_j) + i_tranvec(3, cur_tran)
                  r_sq = dx * dx + dy * dy + dz * dz
                  if (r_sq .lt. cutlist_sq) then
                    if (excl_img_flags(img_j) .eq. 0) then
                      if (ico(iaci + img_iac(img_j)) .eq. 0) then
                        if (r_sq .lt. es_cut_sq) then
                          ee_eval_cnt = ee_eval_cnt + 1
                          img_j_ee_eval(ee_eval_cnt) = &
                            ior(img_j,ishft(cur_tran, 27))
#ifdef MPI
                          if (used_img_map(img_j) .eq. 0) then
                            used_img_map(img_j) = 1
                            used_img_cnt = used_img_cnt + 1
                            used_img_lst(used_img_cnt) = img_j
                          end if
#endif /* MPI */
                        end if
                      else
                        full_eval_cnt = full_eval_cnt + 1
                        img_j_full_eval(full_eval_cnt) = &
                          ior(img_j,ishft(cur_tran,27))
#ifdef MPI
                        if (used_img_map(img_j) .eq. 0) then
                          used_img_map(img_j) = 1
                          used_img_cnt = used_img_cnt + 1
                          used_img_lst(used_img_cnt) = img_j
                        end if
#endif /* MPI */
                      end if
#ifdef MPI
                    else
                      if (used_img_ips_map(img_j) .eq. 0) then
                        used_img_ips_map(img_j) = 1
                        used_img_ips_cnt = used_img_ips_cnt + 1
                        used_img_ips_lst(used_img_ips_cnt) = img_j
                      end if
#endif /* MPI */
                    end if
                  end if
                end do

                if (cur_bkt .eq. i_bkt) exit z_loop4 ! Done with cit buckets
                                                     ! img i claims pairs from.
              end do x_loop4
            end do y_loop4
          end do z_loop4
        end if
      end if

      total_pairs = total_pairs + full_eval_cnt + ee_eval_cnt       ! DBG

      if (dont_skip_belly_pairs) then
        call pack_nb_list(ee_eval_cnt, img_j_ee_eval, &
                          full_eval_cnt, img_j_full_eval, &
                          gbl_sluipairs, num_packed)
      else
        call pack_nb_list_skip_belly_pairs(ee_eval_cnt, img_j_ee_eval, &
                                           full_eval_cnt, img_j_full_eval, &
                                           gbl_sluipairs, igroup, img_atm_map, &
                                           num_packed)
      end if

      ! Clear excluded j images flags:

      do i = atm_mask_idx + 1, atm_mask_idx + atm_maskdata(atm_i)%cnt
        excl_img_flags(atm_img_map(atm_mask(i))) = 0
      end do

      if (ifail .eq. 1) then
        flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi
        return
      end if

    end do ! end of images in i_bkt

    flat_cit(i_bkt)%img_hi = saved_i_bkt_img_hi

  end do ! end of all cit buckets with i images owned by this task

! write(0, *) 'DBG: total_pairs, listtot = ', &
!             total_pairs, num_packed
! write(0, *) 'DBG: avg packing = ', &
!             dble(total_pairs)/dble(num_packed)

  if (verbose .gt. 0) then
    if (master) write(mdout, *) 'listtot = ', num_packed
  end if

  return

contains

#include "pack_nb_list_dflt.i"

end subroutine get_nb_ips_list

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  check_my_atom_movement
!
! Description:  This routine checks if any atom on this task's atom list has
!               moved more than half the skin distance.  This is intended for
!               use only with parallel MD.
!
!*******************************************************************************

subroutine check_my_atom_movement(crd, saved_crd, my_atm_lst, skinnb, ntp, &
                                  new_list)

  use parallel_dat_mod
  use pbc_mod

  implicit none

! Formal arguments:

  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  integer               :: my_atm_lst(*)
  double precision      :: skinnb
  integer               :: ntp
  logical               :: new_list

! Local variables:

  integer               :: atm_lst_idx
  double precision      :: dx, dy, dz
  double precision      :: box_dx, box_dy, box_dz
  double precision      :: distance_limit
  integer               :: n

! If anybody moves more than half the nbskin added to cutoff in generating the
! pairlist, update the pairlist.

  distance_limit = 0.25d0 * skinnb * skinnb

  if (ntp .eq. 0) then  ! Constant volume

    do atm_lst_idx = 1, my_atm_cnt
      n = my_atm_lst(atm_lst_idx)

      dx = crd(1, n) - saved_crd(1, n)
      dy = crd(2, n) - saved_crd(2, n)
      dz = crd(3, n) - saved_crd(3, n)

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  else  ! Constant pressure scaling.

    box_dx = pbc_box(1) / gbl_saved_box(1)
    box_dy = pbc_box(2) / gbl_saved_box(2)
    box_dz = pbc_box(3) / gbl_saved_box(3)

    do atm_lst_idx = 1, my_atm_cnt
      n = my_atm_lst(atm_lst_idx)

      dx = crd(1, n) - saved_crd(1, n) * box_dx
      dy = crd(2, n) - saved_crd(2, n) * box_dy
      dz = crd(3, n) - saved_crd(3, n) * box_dz

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  end if

  new_list = .false.

  return

end subroutine check_my_atom_movement
#endif /* MPI */

!*******************************************************************************
!
! Subroutine:  check_all_atom_movement
!
! Description:  This routine checks if any atom has moved more than half the
!               skin distance.  This is intended for all minimizations and
!               for single processor MD.
!*******************************************************************************

subroutine check_all_atom_movement(atm_cnt, crd, saved_crd, skinnb, ntp, &
                                   new_list)

  use pbc_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  double precision      :: skinnb
  integer               :: ntp
  logical               :: new_list

! Local variables:

  double precision      :: dx, dy, dz
  double precision      :: box_dx, box_dy, box_dz
  double precision      :: distance_limit
  integer               :: n

! If anybody moves more than half the nbskin added to cutoff in generating the
! pairlist, update the pairlist.

  distance_limit = 0.25d0 * skinnb * skinnb

  if (ntp .eq. 0) then  ! Constant volume

    do n = 1, atm_cnt

      dx = crd(1, n) - saved_crd(1, n)
      dy = crd(2, n) - saved_crd(2, n)
      dz = crd(3, n) - saved_crd(3, n)

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  else  ! Constant pressure scaling.

    box_dx = pbc_box(1) / gbl_saved_box(1)
    box_dy = pbc_box(2) / gbl_saved_box(2)
    box_dz = pbc_box(3) / gbl_saved_box(3)

    do n = 1, atm_cnt

      dx = crd(1, n) - saved_crd(1, n) * box_dx
      dy = crd(2, n) - saved_crd(2, n) * box_dy
      dz = crd(3, n) - saved_crd(3, n) * box_dz

      if (dx * dx + dy * dy + dz * dz .ge. distance_limit) then
        new_list = .true.
        return
      end if

    end do

  end if

  new_list = .false.
  return

end subroutine check_all_atom_movement

!*******************************************************************************
!
! Subroutine:  save_all_atom_crds
!
! Description:  This is needed for the skin test for buffered pairlists.
!*******************************************************************************

subroutine save_all_atom_crds(atm_cnt, crd, saved_crd)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: saved_crd(3, atm_cnt)

! Local variables:

  integer               :: i

  do i = 1, atm_cnt
    saved_crd(1, i) = crd(1, i)
    saved_crd(2, i) = crd(2, i)
    saved_crd(3, i) = crd(3, i)
  end do

  return

end subroutine save_all_atom_crds

!*******************************************************************************
!
! Subroutine:  save_imgcrds
!
! Description:  <TBS> 
!
!*******************************************************************************

#ifdef MPI
subroutine save_imgcrds(img_crd, saved_imgcrd, used_img_cnt, used_img_lst)
#else
subroutine save_imgcrds(img_cnt, img_crd, saved_imgcrd)
#endif

  use img_mod

  implicit none

! Formal arguments:

#ifdef MPI
  double precision      :: img_crd(3, *)
  double precision      :: saved_imgcrd(3, *)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
#else
  integer               :: img_cnt
  double precision      :: img_crd(3, *)
  double precision      :: saved_imgcrd(3, *)
#endif /* MPI */

! Local variables:

#ifdef MPI
  integer               :: lst_idx
#endif /* MPI */
  integer               :: img_idx

#ifdef MPI
  do lst_idx = 1, used_img_cnt
    img_idx = used_img_lst(lst_idx)
    saved_imgcrd(:, img_idx) = img_crd(:, img_idx)
  end do
#else
  do img_idx = 1, img_cnt
    saved_imgcrd(:,img_idx) = img_crd(:,img_idx)
  end do
#endif

  return

end subroutine save_imgcrds

!*******************************************************************************
!
! Subroutine:  adjust_imgcrds
!
! Description:  <TBS>
!              
!*******************************************************************************

#ifdef MPI
subroutine adjust_imgcrds(img_cnt, img_crd, img_atm_map, &
                          used_img_cnt, used_img_lst, &
                          saved_imgcrd, crd, saved_crd, ntp)
#else
subroutine adjust_imgcrds(img_cnt, img_crd, img_atm_map, &
                          saved_imgcrd, crd, saved_crd, ntp)
#endif /* MPI */

  use gbl_datatypes_mod
  use img_mod
  use pbc_mod

  implicit none

! Formal arguments:

#ifdef MPI
  integer               :: img_cnt
  double precision      :: img_crd(3, *)
  integer               :: img_atm_map(*)
  integer               :: used_img_cnt
  integer               :: used_img_lst(*)
  double precision      :: saved_imgcrd(3, *)
  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  integer               :: ntp
#else
  integer               :: img_cnt
  double precision      :: img_crd(3, *)
  integer               :: img_atm_map(*)
  double precision      :: saved_imgcrd(3, *)
  double precision      :: crd(3, *)
  double precision      :: saved_crd(3, *)
  integer               :: ntp
#endif /* MPI */

! Local variables:

  integer               :: atm_idx
  integer               :: img_idx
  double precision      :: box_del(3)
#ifdef MPI
  integer               :: lst_idx
#endif /* MPI */

#ifdef MPI
  if (ntp .eq. 0) then  ! Constant volume

      do lst_idx = 1, used_img_cnt
        img_idx = used_img_lst(lst_idx)
        atm_idx = img_atm_map(img_idx)
        img_crd(:,img_idx) = crd(:,atm_idx) + saved_imgcrd(:,img_idx) - &
                             saved_crd(:,atm_idx)
      end do

  else  ! Constant pressure scaling.

    box_del(:) = pbc_box(:) / gbl_saved_box(:)

      do lst_idx = 1, used_img_cnt
        img_idx = used_img_lst(lst_idx)
        atm_idx = img_atm_map(img_idx)
        img_crd(:,img_idx) = crd(:,atm_idx) + (saved_imgcrd(:,img_idx) - &
                             saved_crd(:,atm_idx)) * box_del(:)
      end do

  end if

#else

  if (ntp .eq. 0) then  ! Constant volume

    do img_idx = 1, img_cnt
        atm_idx = img_atm_map(img_idx)
        img_crd(:,img_idx) = crd(:,atm_idx) + saved_imgcrd(:,img_idx) - &
                             saved_crd(:,atm_idx)
    end do

  else  ! Constant pressure scaling.

    box_del(:) = pbc_box(:) / gbl_saved_box(:)

    do img_idx = 1, img_cnt
      atm_idx = img_atm_map(img_idx)
      img_crd(:,img_idx) = crd(:,atm_idx) + (saved_imgcrd(:,img_idx) - &
                           saved_crd(:,atm_idx)) * box_del(:)
    end do

  end if

#endif /* MPI */

  return

end subroutine adjust_imgcrds

end module nb_pairlist_mod
