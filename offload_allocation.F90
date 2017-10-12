#ifdef MIC_offload
/*eed_cub, efs_tbll, and fs_tbl are not trasnferred inside, since now only erfc is used instead of lookup table*/

module offload_allocation_mod
use cit_mod
use gbl_datatypes_mod
use omp_lib
use file_io_dat_mod

!dec$ options /offload_attribute_target=mic
!for direct force
  double precision,allocatable      ::  img_frc_thread_dev(:,:)
  double precision,dimension(:,:),allocatable 	:: img_frc_dev 
  double precision,dimension(:,:),allocatable 	:: img_crd_dev 
  integer,dimension(:),allocatable 	:: ipairs_dev 
  integer,dimension(:),allocatable 	:: start_indexes_dev 
  double precision,dimension(:,:),allocatable 	:: x_tran_dev
  double precision,dimension(:,:),allocatable 	:: tranvec_dev
  double precision,dimension(:),allocatable 	:: img_qterm_dev
  double precision,dimension(:),allocatable 	:: eed_cub_dev
  double precision,dimension(:),allocatable 	:: gbl_cn1_dev
  double precision,dimension(:),allocatable 	:: gbl_cn2_dev
  double precision,dimension(:),allocatable 	:: efs_tbl_dev
  double precision,dimension(:),allocatable 	:: fs_tbl_dev
  integer,dimension(:),allocatable 	:: typ_ico_dev
  integer,dimension(:),allocatable 	:: gbl_img_iac_dev
  integer                :: num_threads_offload
!for nb list
  type(cit_tbl_rec),dimension(:),allocatable  :: flat_cit_dev
  double precision,dimension(:,:),allocatable  :: fraction_dev
  integer,dimension(:),allocatable  :: x_bkts_dev
  integer,dimension(:),allocatable  :: x_trans_dev
  integer,dimension(:),allocatable  :: y_bkts_dev
  integer,dimension(:),allocatable  :: y_trans_dev
  integer,dimension(:),allocatable  :: z_bkts_dev
  integer,dimension(:),allocatable  :: z_trans_dev
!  integer,dimension(:),allocatable  :: used_img_cnt_thread_dev
  integer(byte),dimension(:),allocatable  :: used_img_lst_thread_dev
  integer,dimension(:),allocatable  :: thread_num_packed_dev
  integer,dimension(:),allocatable  :: thread_ipairs_dev
  integer,dimension(:),allocatable  :: ifail_thread_dev
  integer,dimension(:),allocatable  :: excl_img_flags_dev
  integer,dimension(:),allocatable  :: atm_img_map_dev
  integer,dimension(:),allocatable  :: img_atm_map_dev
!  integer(byte), dimension(:),allocatable  :: used_img_map_thread_dev
  type(listdata_rec),dimension(:),allocatable  :: atm_maskdata_dev
  integer,dimension(:),allocatable  :: atm_mask_dev
  integer,dimension(:),allocatable  :: img_j_full_eval_dev
  integer,dimension(:),allocatable  :: img_j_ee_eval_dev
  double precision,dimension(:,:),allocatable  :: i_tranvec_dev
  integer                                          :: local_next,local_next_mult_fac
  integer                                    :: num_packed, num_packed_dev
  integer				:: max_nb_per_atm
  double precision                      :: mic_fraction
  !Reuse: img_iac use gbl_img_iac_dev
  !tranvec,typ_ico_dev for ico,img_crd_dev 

  integer  :: ipairs_maxsize_dev
  integer  :: used_img_cnt_dev
  integer (kind=kmp_affinity_mask_kind) ::  mask
  integer :: tnum       
  integer :: proc       
  integer :: threads_per_core       
  integer :: mytaskid_dev , temp, offload_numtasks      

!dec$ end options

contains

subroutine offload_allocate

  use img_mod
  use cit_mod
  use timers_mod
  use mdin_ctrl_dat_mod
  use mdin_ewald_dat_mod
  use parallel_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use ene_frc_splines_mod
!  use pmemd_lib_mod
!  use pme_alltasks_setup_mod
  use nb_exclusions_mod
  use prmtop_dat_mod
  use omp_lib
  integer :: mxeedtab
#if 1 
  integer, parameter            :: ipairs_size_min = 1000
  double precision, parameter   :: ipairs_size_coef = 0.225d0 
  
  mic_fraction = 0.80
  local_next = next !it is used here and in the nb_pairlist
  local_next_mult_fac = next_mult_fac !it is used here and in the nb_pairlist
  max_nb_per_atm = ( (vdw_cutoff+skinnb)**3 ) * 2 
  ipairs_maxsize_dev = int(ipairs_size_coef * (max_nb_per_atm + 3.d0) * dble(natom))
#ifdef MPI
  if (numtasks .gt. 0) ipairs_maxsize_dev = ipairs_maxsize_dev/4 
#endif /* MPI */

  if (ipairs_maxsize_dev .lt. ipairs_size_min) ipairs_maxsize_dev = ipairs_size_min


#endif
offload_numtasks = 2 /*Now 2 MPI in the middle (such as MPI 11 and 12 of 24 MPI) are offloading in a single node*/
mytaskid_dev = numtasks/2 - mytaskid /*mytaskid_dev is the transformation of the offloading mpi from 0 to offload_numtasks-1  to calculate the openmp affinity in the below*/

!dir$  offload begin mandatory target(mic:0) &
in(mytaskid_dev, offload_numtasks) &
nocopy(tnum, threads_per_core, proc, temp)&
out(num_threads_offload)
        !$omp parallel
                num_threads_offload = omp_get_num_threads()
        !$omp end parallel
        if(num_threads_offload .gt. 120) then
               !Write statement needs to be fed to master in order to write to mdout
                ! write(mdout, '(a)') "OpenMP thread numbers are not set or too high. Please set it <= 120 with proper affinity and MIC partitions, and re-rerun. Now Forcing # of threads to 28 per MPI on offload_allocate routine"
                call OMP_SET_NUM_THREADS(30) 
        end if
                 !$omp parallel private(mask, tnum, thread_per_core, proc,temp)
                        num_threads_offload = omp_get_num_threads() /*total threads per MPI*/
                        tnum = omp_get_thread_num()  /*this thread number*/
                        call kmp_create_affinity_mask(mask) /*initializing mask*/
                        /*scatter affinity*/    
                        proc = (240/offload_numtasks) * mytaskid_dev
                        proc = proc + (tnum - (tnum/30)*30)*4 + (tnum/30) + 1
                        temp = kmp_set_affinity_mask_proc(proc, mask) /*mapping=affinitize of omp threads to HW threads*/
                        ! Write statement needs to be fed to master in order to
                        ! write to mdout
                        !if (kmp_set_affinity(mask) .ne.  0) then
                        !        write(mdout, '(a)') "Could not set affinity mask in offload allocation"
                        !end if
                !$omp end parallel
!dir$ end offload

mxeedtab = int(ew_coeff * eedtbdns * es_cutoff * 1.5d0) 

allocate(img_frc_thread_dev(1:3,1:natom*num_threads_offload))
allocate(img_frc_dev(1:3,1:natom))
allocate(img_crd_dev(1:3,1:natom))
allocate(ipairs_dev(0:ipairs_maxsize_dev)) 
allocate(start_indexes_dev(1:natom))!natom replace by natom/#MPI 
allocate(x_tran_dev(1:3,0:17))
allocate(tranvec_dev(1:3,0:17))
allocate(img_qterm_dev(1:natom)) 
!allocate(eed_cub_dev(1:mxeedtab)) 
allocate(gbl_cn1_dev(1:nttyp)) 
allocate(gbl_cn2_dev(1:nttyp)) 
!allocate(efs_tbl_dev(1:8*efs_tbl_siz)) 
!allocate(fs_tbl_dev(1:4*efs_tbl_siz)) 
allocate(typ_ico_dev(1:ntypes*ntypes)) 
allocate(gbl_img_iac_dev(1:natom)) 
!nb list related
allocate(flat_cit_dev(0 :cit_tbl_x_dim*cit_tbl_y_dim*cit_tbl_z_dim - 1))
allocate(fraction_dev(1:3,1:natom)) 
allocate(x_bkts_dev(0: cit_tbl_x_dim*3 -1)) 
allocate(x_trans_dev(0: cit_tbl_x_dim*3 -1)) 
allocate(y_bkts_dev(0: cit_tbl_y_dim*3 -1)) 
allocate(y_trans_dev(0: cit_tbl_y_dim*3 -1)) 
allocate(z_bkts_dev(0: cit_tbl_z_dim*2 -1)) 
allocate(z_trans_dev(0: cit_tbl_z_dim*2 -1)) 
!allocate(used_img_cnt_thread_dev(0:num_threads_offload)) 
allocate(used_img_lst_thread_dev(1:natom)) 
allocate(thread_num_packed_dev(0:num_threads_offload)) 
allocate(thread_ipairs_dev(0:ipairs_maxsize_dev)) 
allocate(ifail_thread_dev(0:  num_threads_offload)) 
allocate(excl_img_flags_dev(1: natom)) 
allocate(atm_img_map_dev(1: natom)) 
allocate(img_atm_map_dev(1: natom)) 
!allocate(used_img_map_thread_dev(1: natom)) 
allocate(atm_maskdata_dev(1: natom)) 
!allocate(atm_mask_dev(1: 291192)) 
allocate(atm_mask_dev(1: next*next_mult_fac)) 
!allocate(igroup_dev(1: natom)) 
allocate(img_j_full_eval_dev(1:max_nb_per_atm)) 
allocate(img_j_ee_eval_dev(1:max_nb_per_atm)) 
allocate(i_tranvec_dev(1:3,0:17)) 
!allocate(ico_dev(1:ntypes*ntypes)) 

!dir$  offload begin   mandatory target(mic:0) &
!direct force related
nocopy(img_frc_thread_dev(:,1:natom*num_threads_offload): alloc_if(.true.) free_if(.false.) ) & 
nocopy(img_frc_dev(:,1:natom) : alloc_if(.true.) free_if(.false.) ) & 
nocopy(img_crd_dev(:,1:natom) : alloc_if(.true.) free_if(.false.) ) & 
!nocopy(ipairs_dev(1:natom/numtasks*400) : alloc_if(.true.) free_if(.false.) ) & 
nocopy(ipairs_dev(0:ipairs_maxsize_dev) : alloc_if(.true.) free_if(.false.) ) & 
nocopy(start_indexes_dev(1:natom) : alloc_if(.true.) free_if(.false.)) &  
nocopy(x_tran_dev(:,0:17) : alloc_if(.true.) free_if(.false.)) &  
nocopy(tranvec_dev(:,0:17) : alloc_if(.true.) free_if(.false.)) &  
nocopy(img_qterm_dev(1:natom) : alloc_if(.true.) free_if(.false.))  &
!nocopy(eed_cub_dev(1:mxeedtab) : alloc_if(.true.) free_if(.false.))  &
in(gbl_cn1(1:nttyp) : alloc_if(.true.) free_if(.false.) into(gbl_cn1_dev))  &
in(gbl_cn2(1:nttyp) : alloc_if(.true.) free_if(.false.) into(gbl_cn2_dev))  &
!in(efs_tbl(1:8*efs_tbl_siz) : alloc_if(.true.) free_if(.false.) into(efs_tbl_dev))  &
!in(fs_tbl(1:4*efs_tbl_siz) : alloc_if(.true.) free_if(.false.) into(fs_tbl_dev))  &
in(typ_ico(1:ntypes*ntypes) : alloc_if(.true.) free_if(.false.) into(typ_ico_dev))  &
nocopy(gbl_img_iac_dev(1:natom) : alloc_if(.true.) free_if(.false.))& 
!nblist related
nocopy(flat_cit_dev(0 :cit_tbl_x_dim*cit_tbl_y_dim*cit_tbl_z_dim - 1) : alloc_if(.true.) free_if(.false.))& 
nocopy(fraction_dev(:,1:natom) : alloc_if(.true.) free_if(.false.))& 
nocopy(x_bkts_dev(0: cit_tbl_x_dim*3 -1) : alloc_if(.true.) free_if(.false.)) &
nocopy(x_trans_dev(0: cit_tbl_x_dim*3 -1) : alloc_if(.true.) free_if(.false.)) &
nocopy(y_bkts_dev(0: cit_tbl_y_dim*3 -1) : alloc_if(.true.) free_if(.false.)) &
nocopy(y_trans_dev(0: cit_tbl_y_dim*3 -1) : alloc_if(.true.) free_if(.false.)) & 
nocopy(z_bkts_dev(0: cit_tbl_z_dim*2 -1) : alloc_if(.true.) free_if(.false.)) &
nocopy(z_trans_dev(0: cit_tbl_z_dim*2 -1) : alloc_if(.true.) free_if(.false.)) & 
!nocopy(used_img_cnt_thread_dev(0:num_threads_offload) : alloc_if(.true.) free_if(.false.)) &
nocopy(used_img_lst_thread_dev(1:natom) : alloc_if(.true.) free_if(.false.))& 
nocopy(thread_num_packed_dev(0:num_threads_offload) : alloc_if(.true.) free_if(.false.)) &
nocopy(thread_ipairs_dev(0:ipairs_maxsize_dev) : alloc_if(.true.) free_if(.false.))& 
nocopy(ifail_thread_dev(0:  num_threads_offload) : alloc_if(.true.) free_if(.false.))& 
nocopy(excl_img_flags_dev(1: natom) : alloc_if(.true.) free_if(.false.)) &
nocopy(atm_img_map_dev(1: natom) : alloc_if(.true.) free_if(.false.)) &
nocopy(img_atm_map_dev(1: natom) : alloc_if(.true.) free_if(.false.))& 
!nocopy(used_img_map_thread_dev(1: natom) : alloc_if(.true.) free_if(.false.)) &
!nocopy(atm_maskdata_dev(1: natom) : alloc_if(.true.) free_if(.false.)) &
in(atm_nb_maskdata(1: natom) : alloc_if(.true.) free_if(.false.) into(atm_maskdata_dev(:))) &
!nocopy(atm_mask_dev(1: 291192) : alloc_if(.true.) free_if(.false.))& 
in(atm_nb_mask(1: next*next_mult_fac) : alloc_if(.true.) free_if(.false.) into(atm_mask_dev(:)))& 
!nocopy(igroup_dev(1: natom) : alloc_if(.true.) free_if(.false.)) &
nocopy(img_j_full_eval_dev(1:max_nb_per_atm) : alloc_if(.true.) free_if(.false.))& 
nocopy(img_j_ee_eval_dev(1:max_nb_per_atm) : alloc_if(.true.) free_if(.false.)) &
nocopy(i_tranvec_dev(:,0:17) : alloc_if(.true.) free_if(.false.))
!nocopy(ico_dev(1:ntypes*ntypes) : alloc_if(.true.) free_if(.false.))
!dir$ end  offload
end subroutine offload_allocate


subroutine offload_deallocate


!dir$  offload begin   mandatory target(mic:0) &
/*direct force related*/
nocopy(img_frc_thread_dev(:,:): alloc_if(.false.) free_if(.true.) ) & 
nocopy(img_frc_dev(:,:) : alloc_if(.false.) free_if(.true.) ) & 
nocopy(img_crd_dev(:,:) : alloc_if(.false.) free_if(.true.) ) & 
!nocopy(ipairs_dev(:) : alloc_if(.false.) free_if(.true.) ) & 
nocopy(ipairs_dev(0:) : alloc_if(.false.) free_if(.true.) ) & 
nocopy(start_indexes_dev(:) : alloc_if(.false.) free_if(.true.)) &  
nocopy(x_tran_dev(:,:) : alloc_if(.false.) free_if(.true.)) &  
nocopy(tranvec_dev(:,:) : alloc_if(.false.) free_if(.true.)) &  
nocopy(img_qterm_dev(:) : alloc_if(.false.) free_if(.true.))  &
!nocopy(eed_cub_dev(:) : alloc_if(.false.) free_if(.true.))  &
nocopy(gbl_cn1_dev(:) : alloc_if(.false.) free_if(.true.))  &
nocopy(gbl_cn2_dev(:) : alloc_if(.false.) free_if(.true.))  &
!nocopy(efs_tbl_dev(:) : alloc_if(.false.) free_if(.true.))  &
!nocopy(fs_tbl_dev(:) : alloc_if(.false.) free_if(.true.))  &
nocopy(typ_ico_dev(:) : alloc_if(.false.) free_if(.true.))  &
nocopy(gbl_img_iac_dev(:) : alloc_if(.false.) free_if(.true.)) & 
/*nblist related*/
nocopy(flat_cit_dev(:) : alloc_if(.false.) free_if(.true.))& 
nocopy(fraction_dev(:,:) : alloc_if(.false.) free_if(.true.))& 
nocopy(x_bkts_dev(0:) : alloc_if(.false.) free_if(.true.)) &
nocopy(x_trans_dev(:) : alloc_if(.false.) free_if(.true.)) &
nocopy(y_bkts_dev(:) : alloc_if(.false.) free_if(.true.)) &
nocopy(y_trans_dev(:) : alloc_if(.false.) free_if(.true.)) & 
nocopy(z_bkts_dev(:) : alloc_if(.false.) free_if(.true.)) &
nocopy(z_trans_dev(:) : alloc_if(.false.) free_if(.true.)) & 
!nocopy(used_img_cnt_thread_dev(:) : alloc_if(.false.) free_if(.true.)) &
nocopy(used_img_lst_thread_dev(:) : alloc_if(.false.) free_if(.true.))& 
nocopy(thread_num_packed_dev(:) : alloc_if(.false.) free_if(.true.)) &
nocopy(thread_ipairs_dev(:) : alloc_if(.false.) free_if(.true.))& 
nocopy(ifail_thread_dev(:) : alloc_if(.false.) free_if(.true.))& 
nocopy(excl_img_flags_dev(:) : alloc_if(.false.) free_if(.true.)) &
nocopy(atm_img_map_dev(:) : alloc_if(.false.) free_if(.true.)) &
nocopy(img_atm_map_dev(:) : alloc_if(.false.) free_if(.true.))& 
!nocopy(used_img_map_thread_dev(:) : alloc_if(.false.) free_if(.true.)) &
!nocopy(atm_maskdata_dev(:) : alloc_if(.false.) free_if(.true.)) &
nocopy(atm_maskdata_dev(:) : alloc_if(.false.) free_if(.true.)) &
nocopy(atm_mask_dev(:) : alloc_if(.false.) free_if(.true.))& 
nocopy(img_j_full_eval_dev(:) : alloc_if(.false.) free_if(.true.))& 
nocopy(img_j_ee_eval_dev(:) : alloc_if(.false.) free_if(.true.)) &
nocopy(i_tranvec_dev(:,:) : alloc_if(.false.) free_if(.true.))
!dir$ end  offload

deallocate(img_frc_thread_dev)
deallocate(img_frc_dev)
deallocate(img_crd_dev)
deallocate(ipairs_dev) 
deallocate(start_indexes_dev) 
deallocate(x_tran_dev)
deallocate(tranvec_dev)
deallocate(img_qterm_dev) 
!deallocate(eed_cub_dev) 
deallocate(gbl_cn1_dev) 
deallocate(gbl_cn2_dev) 
!deallocate(efs_tbl_dev) 
!deallocate(fs_tbl_dev) 
deallocate(typ_ico_dev) 
deallocate(gbl_img_iac_dev) 
deallocate(flat_cit_dev)
deallocate(fraction_dev) 
deallocate(x_bkts_dev) 
deallocate(x_trans_dev) 
deallocate(y_bkts_dev) 
deallocate(y_trans_dev) 
deallocate(z_bkts_dev) 
deallocate(z_trans_dev) 
!deallocate(used_img_cnt_thread_dev) 
deallocate(used_img_lst_thread_dev) 
deallocate(thread_num_packed_dev) 
deallocate(thread_ipairs_dev) 
deallocate(ifail_thread_dev) 
deallocate(excl_img_flags_dev) 
deallocate(atm_img_map_dev) 
deallocate(img_atm_map_dev) 
!deallocate(used_img_map_thread_dev) 
deallocate(atm_maskdata_dev) 
!deallocate(atm_mask_dev) 
deallocate(atm_mask_dev) 
!deallocate(igroup_dev) 
deallocate(img_j_full_eval_dev) 
deallocate(img_j_ee_eval_dev) 
deallocate(i_tranvec_dev) 

end subroutine offload_deallocate

end module offload_allocation_mod

#endif /*MIC_offload*/
