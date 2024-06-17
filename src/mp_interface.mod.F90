#include "cpmd_global.h"

! ==================================================================
! == INTERFACE FOR MESSAGE PASSING CALLS:                         ==
! ==================================================================

MODULE mp_interface
  USE, INTRINSIC :: iso_c_binding,     ONLY: c_loc,&
                                             c_f_pointer,&
                                             c_ptr
  USE error_handling,                  ONLY: stopgm
  USE gpu,                             ONLY: comm_buffers_on_host
  USE kinds,                           ONLY: int_1,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE machine,                         ONLY: m_flush,&
                                             m_walltime
  USE para_global,                     ONLY: para_use_mpi_in_place,&
                                             para_buff,&
                                             il_para_buff,&
                                             buff_size_in_bytes
  USE parac,                           ONLY: parai,&
                                             paral
  USE pstat
  USE reshaper,                        ONLY: reshape_inplace
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif
  USE system,                          ONLY: cnti
  USE zeroing_utils,                   ONLY: zeroing

#ifdef __PARALLEL
    USE mpi_f08
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mp_start
  PUBLIC :: mp_end
  PUBLIC :: mp_is_comm_null
  PUBLIC :: mp_environ
  PUBLIC :: mp_comm_free
  PUBLIC :: mp_comm_dup
  PUBLIC :: mp_comm_compare
  PUBLIC :: mp_comm_info
  PUBLIC :: mp_task_query
  PUBLIC :: mp_cart
  PUBLIC :: mp_mpi_error_assert
  PUBLIC :: mp_bcast
  PUBLIC :: mp_all2all
  PUBLIC :: mp_sync
  PUBLIC :: mp_sum
  PUBLIC :: mp_prod
  PUBLIC :: mp_max
  PUBLIC :: mp_min
  PUBLIC :: mp_dims_create
  PUBLIC :: mp_split
  PUBLIC :: mp_split_type
  PUBLIC :: mp_group
  PUBLIC :: mp_send
  PUBLIC :: mp_recv
  PUBLIC :: mp_sendrecv
  PUBLIC :: mp_get_processor_name
  PUBLIC :: mp_get_node_env
  PUBLIC :: mp_get_version
  PUBLIC :: mp_get_library_version
  PUBLIC :: mp_win_alloc_shared_mem
  PUBLIC :: mp_win_sync
  PUBLIC :: mp_win_lock_all_shared
  PUBLIC :: mp_win_unlock_all_shared
  !
  ! interfaces
  !
  INTERFACE mp_bcast
     MODULE PROCEDURE mp_bcast_char_r0
     MODULE PROCEDURE mp_bcast_char_r1
     MODULE PROCEDURE mp_bcast_int4_r0
     MODULE PROCEDURE mp_bcast_int4_r1
     MODULE PROCEDURE mp_bcast_int4_r2
     MODULE PROCEDURE mp_bcast_int8_r0
     MODULE PROCEDURE mp_bcast_int8_r1
     MODULE PROCEDURE mp_bcast_real8_r0
     MODULE PROCEDURE mp_bcast_real8_r1
     MODULE PROCEDURE mp_bcast_real8_r2
     MODULE PROCEDURE mp_bcast_real8_r3
     MODULE PROCEDURE mp_bcast_real8_r4
     MODULE PROCEDURE mp_bcast_real8_r5
     MODULE PROCEDURE mp_bcast_complex8_r0
     MODULE PROCEDURE mp_bcast_complex8_r1
     MODULE PROCEDURE mp_bcast_complex8_r2
     MODULE PROCEDURE mp_bcast_complex8_r3
     MODULE PROCEDURE mp_bcast_complex8_r4
     MODULE PROCEDURE mp_bcast_logical_r0
     MODULE PROCEDURE mp_bcast_logical_r1
     MODULE PROCEDURE mp_bcast_logical_r2
  END INTERFACE mp_bcast


  INTERFACE mp_all2all
     MODULE PROCEDURE mp_all2all_int4_r1
     MODULE PROCEDURE mp_all2all_int8_r1
     MODULE PROCEDURE mp_all2all_real8_r1
     MODULE PROCEDURE mp_all2all_complex4_r1
     MODULE PROCEDURE mp_all2all_complex8_r1
  END INTERFACE mp_all2all


  INTERFACE mp_send
     MODULE PROCEDURE mp_send_char_r0
     MODULE PROCEDURE mp_send_char_r1
     MODULE PROCEDURE mp_send_int4_r0
     MODULE PROCEDURE mp_send_int4_r1
     MODULE PROCEDURE mp_send_int4_r2
     MODULE PROCEDURE mp_send_int8_r0
     MODULE PROCEDURE mp_send_int8_r1
     MODULE PROCEDURE mp_send_real8_r0
     MODULE PROCEDURE mp_send_real8_r1
     MODULE PROCEDURE mp_send_real8_r2
     MODULE PROCEDURE mp_send_real8_r3
     MODULE PROCEDURE mp_send_real8_r4
     MODULE PROCEDURE mp_send_real8_r5
     MODULE PROCEDURE mp_send_complex4_r0
     MODULE PROCEDURE mp_send_complex4_r1
     MODULE PROCEDURE mp_send_complex4_r2
     MODULE PROCEDURE mp_send_complex4_r3
     MODULE PROCEDURE mp_send_complex4_r4
     MODULE PROCEDURE mp_send_complex8_r0
     MODULE PROCEDURE mp_send_complex8_r1
     MODULE PROCEDURE mp_send_complex8_r2
     MODULE PROCEDURE mp_send_complex8_r3
     MODULE PROCEDURE mp_send_complex8_r4
     MODULE PROCEDURE mp_send_logical_r0
     MODULE PROCEDURE mp_send_logical_r1
     MODULE PROCEDURE mp_send_logical_r2
  END INTERFACE mp_send


  INTERFACE mp_recv
     MODULE PROCEDURE mp_recv_char_r0
     MODULE PROCEDURE mp_recv_char_r1
     MODULE PROCEDURE mp_recv_int4_r0
     MODULE PROCEDURE mp_recv_int4_r1
     MODULE PROCEDURE mp_recv_int4_r2
     MODULE PROCEDURE mp_recv_int8_r0
     MODULE PROCEDURE mp_recv_int8_r1
     MODULE PROCEDURE mp_recv_real8_r0
     MODULE PROCEDURE mp_recv_real8_r1
     MODULE PROCEDURE mp_recv_real8_r2
     MODULE PROCEDURE mp_recv_real8_r3
     MODULE PROCEDURE mp_recv_real8_r4
     MODULE PROCEDURE mp_recv_real8_r5
     MODULE PROCEDURE mp_recv_complex4_r0
     MODULE PROCEDURE mp_recv_complex4_r1
     MODULE PROCEDURE mp_recv_complex4_r2
     MODULE PROCEDURE mp_recv_complex4_r3
     MODULE PROCEDURE mp_recv_complex4_r4
     MODULE PROCEDURE mp_recv_complex8_r0
     MODULE PROCEDURE mp_recv_complex8_r1
     MODULE PROCEDURE mp_recv_complex8_r2
     MODULE PROCEDURE mp_recv_complex8_r3
     MODULE PROCEDURE mp_recv_complex8_r4
     MODULE PROCEDURE mp_recv_logical_r0
     MODULE PROCEDURE mp_recv_logical_r1
     MODULE PROCEDURE mp_recv_logical_r2
  END INTERFACE mp_recv


  INTERFACE mp_sendrecv
     MODULE PROCEDURE mp_sendrecv_int4_r1
     MODULE PROCEDURE mp_sendrecv_int4_r2
     MODULE PROCEDURE mp_sendrecv_int8_r1
     MODULE PROCEDURE mp_sendrecv_real8_r1
     MODULE PROCEDURE mp_sendrecv_real8_r2
     MODULE PROCEDURE mp_sendrecv_real8_r3
     MODULE PROCEDURE mp_sendrecv_real8_r4
     MODULE PROCEDURE mp_sendrecv_complex8_r1
     MODULE PROCEDURE mp_sendrecv_complex8_r2
     MODULE PROCEDURE mp_sendrecv_complex8_r3
     MODULE PROCEDURE mp_sendrecv_complex8_r4
  END INTERFACE mp_sendrecv


  INTERFACE mp_sum
     MODULE PROCEDURE mp_sum_int4_r0
     MODULE PROCEDURE mp_sum_int4_r1
     MODULE PROCEDURE mp_sum_int4_r2
     MODULE PROCEDURE mp_sum_int4_r3
     MODULE PROCEDURE mp_sum_int4_r4
     MODULE PROCEDURE mp_sum_real8_r0
     MODULE PROCEDURE mp_sum_real8_r1
     MODULE PROCEDURE mp_sum_real8_r2
     MODULE PROCEDURE mp_sum_real8_r3
     MODULE PROCEDURE mp_sum_real8_r4
     MODULE PROCEDURE mp_sum_real4_r0
     MODULE PROCEDURE mp_sum_real4_r1
     MODULE PROCEDURE mp_sum_real4_r2
     MODULE PROCEDURE mp_sum_real4_r3
     MODULE PROCEDURE mp_sum_real4_r4
     MODULE PROCEDURE mp_sum_complex8_r0
     MODULE PROCEDURE mp_sum_complex8_r1
     MODULE PROCEDURE mp_sum_complex8_r2
     MODULE PROCEDURE mp_sum_complex8_r3
     MODULE PROCEDURE mp_sum_complex8_r4
     MODULE PROCEDURE mp_sum_root_int4_r0
     MODULE PROCEDURE mp_sum_root_int4_r1
     MODULE PROCEDURE mp_sum_root_int4_r2
     MODULE PROCEDURE mp_sum_root_real8_r0
     MODULE PROCEDURE mp_sum_root_real8_r1
     MODULE PROCEDURE mp_sum_root_real8_r2
     MODULE PROCEDURE mp_sum_root_complex8_r0
     MODULE PROCEDURE mp_sum_root_complex8_r1
     MODULE PROCEDURE mp_sum_root_complex8_r2
     MODULE PROCEDURE mp_sum_root_in_place_int4_r1
     MODULE PROCEDURE mp_sum_root_in_place_int4_r2
     MODULE PROCEDURE mp_sum_root_in_place_int8_r1
     MODULE PROCEDURE mp_sum_root_in_place_real8_r1
     MODULE PROCEDURE mp_sum_root_in_place_real8_r2
     MODULE PROCEDURE mp_sum_root_in_place_real8_r3
     MODULE PROCEDURE mp_sum_root_in_place_real8_r4
     MODULE PROCEDURE mp_sum_root_in_place_real8_r5
     MODULE PROCEDURE mp_sum_root_in_place_complex8_r1
     MODULE PROCEDURE mp_sum_root_in_place_complex8_r2
     MODULE PROCEDURE mp_sum_root_in_place_complex8_r3
     MODULE PROCEDURE mp_sum_root_in_place_complex8_r4
     MODULE PROCEDURE mp_sum_in_place_int1_r1
     MODULE PROCEDURE mp_sum_in_place_int4_r0
     MODULE PROCEDURE mp_sum_in_place_int4_r1
     MODULE PROCEDURE mp_sum_in_place_int4_r2
     MODULE PROCEDURE mp_sum_in_place_int8_r0
     MODULE PROCEDURE mp_sum_in_place_int8_r1
     MODULE PROCEDURE mp_sum_in_place_real8_r0
     MODULE PROCEDURE mp_sum_in_place_real8_r1
     MODULE PROCEDURE mp_sum_in_place_real8_r2
     MODULE PROCEDURE mp_sum_in_place_real8_r3
     MODULE PROCEDURE mp_sum_in_place_real8_r4
     MODULE PROCEDURE mp_sum_in_place_real8_r5
     MODULE PROCEDURE mp_sum_in_place_complex8_r0
     MODULE PROCEDURE mp_sum_in_place_complex8_r1
     MODULE PROCEDURE mp_sum_in_place_complex8_r2
     MODULE PROCEDURE mp_sum_in_place_complex8_r3
     MODULE PROCEDURE mp_sum_in_place_complex8_r4
  END INTERFACE mp_sum


  INTERFACE mp_max
     MODULE PROCEDURE mp_max_int4_r0
     MODULE PROCEDURE mp_max_int8_r0
     MODULE PROCEDURE mp_max_int4_r1
     MODULE PROCEDURE mp_max_real8_r0
     MODULE PROCEDURE mp_max_real8_r1
     MODULE PROCEDURE mp_max_complex8_r0
     MODULE PROCEDURE mp_max_complex8_r1
  END INTERFACE mp_max

  INTERFACE mp_min
     MODULE PROCEDURE mp_min_int4_r0
     MODULE PROCEDURE mp_min_int4_r1
     MODULE PROCEDURE mp_min_real8_r0
     MODULE PROCEDURE mp_min_real8_r1
     MODULE PROCEDURE mp_min_complex8_r0
     MODULE PROCEDURE mp_min_complex8_r1     
  END INTERFACE mp_min
  
  INTERFACE mp_prod
     MODULE PROCEDURE mp_prod_int4_r0
     MODULE PROCEDURE mp_prod_int4_r1
     MODULE PROCEDURE mp_prod_int4_r2
     MODULE PROCEDURE mp_prod_int4_r3
     MODULE PROCEDURE mp_prod_int4_r4
     MODULE PROCEDURE mp_prod_real8_r0
     MODULE PROCEDURE mp_prod_real8_r1
     MODULE PROCEDURE mp_prod_real8_r2
     MODULE PROCEDURE mp_prod_real8_r3
     MODULE PROCEDURE mp_prod_real8_r4
     MODULE PROCEDURE mp_prod_complex8_r0
     MODULE PROCEDURE mp_prod_complex8_r1
     MODULE PROCEDURE mp_prod_complex8_r2
     MODULE PROCEDURE mp_prod_complex8_r3
     MODULE PROCEDURE mp_prod_complex8_r4
  END INTERFACE mp_prod


  ! alias for MPI_COMM_WORLD to be used outside this module
#ifdef __PARALLEL
  type(MPI_COMM), SAVE, PUBLIC :: mp_comm_world
#else
  INTEGER, SAVE, PUBLIC :: mp_comm_world
#endif
  INTEGER, SAVE, PUBLIC :: mp_max_processor_name = HUGE(0)
#ifdef __PARALLEL
  type(MPI_COMM), SAVE, PUBLIC :: mp_comm_null
#else
  INTEGER, SAVE, PUBLIC :: mp_comm_null
#endif

  INTEGER, SAVE, PRIVATE :: mp_int1_in_bytes, mp_int2_in_bytes, mp_int4_in_bytes, mp_int8_in_bytes
  INTEGER, SAVE, PRIVATE :: mp_real_in_bytes, mp_double_in_bytes
  INTEGER, SAVE, PRIVATE :: mp_complex_in_bytes, mp_double_complex_in_bytes
  INTEGER, SAVE, PRIVATE :: mp_char_in_bytes
  INTEGER, SAVE, PRIVATE :: mp_logical_in_bytes
  !TK, FAU Erlangen Nuernberg, 03.19 shared memory windows
#ifdef __PARALLEL
  type(MPI_WIN), SAVE, PRIVATE :: mpi_window(2)
#else
  INTEGER, SAVE, PRIVATE :: mpi_window(2)
#endif


CONTAINS


  SUBROUTINE mp_start
    ! ==--------------------------------------------------------------==
    ! == Initialisation of message passing calls                      ==
    ! ==--------------------------------------------------------------==
    ! Variables
    INTEGER :: i
    INTEGER, PARAMETER :: output_unit       =  6
    CHARACTER(*),PARAMETER::procedureN='mp_start'
#ifdef __PARALLEL
    INTEGER :: ierr, provided
    LOGICAL :: externalInit
    ! ==--------------------------------------------------------------==
    !   check if MPI is already initialized (needed for library QM/MM interface mode)
    CALL MPI_INITIALIZED(externalInit,ierr)
    IF(.NOT.externalInit) THEN
       ! if not, call  mpi_init
       !$ provided = 0
       !$ IF ( .TRUE. ) THEN
       !$omp master
#ifdef _HASNT_MULTITHREAD_MPI_SUPPORT
       !$ CALL mpi_init_thread ( MPI_THREAD_FUNNELED, provided, ierr )
       !$ IF (paral%io_parent) THEN
       !$ WRITE(output_unit,*) 'WARNING: GPU not suported without MPI_THREAD_MULTIPLE'
       !$ ENDIF
#else
       !$ CALL mpi_init_thread ( MPI_THREAD_MULTIPLE, provided, ierr )
#endif 
       !$ IF(provided.LT.MPI_THREAD_FUNNELED) THEN
       !$     CALL stopgm(procedureN,'Thread support not provided',&
       !$       __LINE__,__FILE__)
       !$ ENDIF
       !$omp end master
       !$ ELSE
       CALL mpi_init(ierr)
       !$ ENDIF
       IF (ierr.NE.0) CALL stopgm('MPI_INIT','IERR.NE.0',& 
            __LINE__,__FILE__)
       ! an external interface has to set mp_comm_world
    ENDIF
    CALL mpi_errhandler_set ( MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr )
    IF (ierr/=0) CALL stopgm('MPI_INIT','mpi_errhandler_set',&
         __LINE__,__FILE__)
#else
    mp_comm_world = 0
#endif
    CALL zeroing(cmlen)!,cmpar)
    DO i=1,cmpar
       cmcal(i)=1.e-20_real_8
       cmtim(i)=1.e-20_real_8
    ENDDO
    CALL mp_set_atomic_type_sizes()
    CALL mp_set_global_consts()
       mp_comm_world = MPI_COMM_WORLD
  END SUBROUTINE mp_start


  SUBROUTINE mp_end
    ! ==--------------------------------------------------------------==
    ! == End of message passing calls. Display some information       ==
    ! ==--------------------------------------------------------------==
    ! Variables
#ifdef __PARALLEL
    INTEGER :: ierr
#endif
    CHARACTER(*),PARAMETER::procedureN='mp_end'
    CHARACTER (len=50), PARAMETER :: fformat1='(1X,"=",A,T27,F10.0," BYTES",8X,F12.0,2X,"=")' 
    CHARACTER (len=60), PARAMETER :: fformat2='(1X,"=",A,T27,F10.3,"  MB/S",6X,F10.3," SEC",2X,"=")' 
    INTEGER :: i
    REAL(real_8) :: xln
    ! ==--------------------------------------------------------------==
    IF (paral%io_parent.AND.parai%nproc.GT.1) THEN
       DO i=1,cmpar
          IF (cmcal(i).LT.1) THEN
             cmlen(i)=0
             cmcal(i)=1
          ENDIF
          IF (cmtim(i).LT.1.e-8_real_8) THEN
             cmlen(i)=0
             cmtim(i)=1
          ENDIF
       ENDDO
       xln=LOG(REAL(parai%nproc,kind=real_8))/LOG(2._real_8)
       !   ==------------------------------------------------------------==
       WRITE(6,'(/,1X,64("="))')
       WRITE(6,'(1X,"=",A,2X,A,2X,A,2X,"=")') ' COMMUNICATION TASK',&
            'AVERAGE MESSAGE LENGTH','NUMBER OF CALLS'
       WRITE(6,fformat1) ' SEND/RECEIVE',&
            cmlen(ipar_send)/cmcal(ipar_send),cmcal(ipar_send)
       WRITE(6,fformat1) ' BROADCAST',&
            cmlen(ipar_cast)/cmcal(ipar_cast),cmcal(ipar_cast)
       WRITE(6,fformat1) ' GLOBAL SUMMATION',&
            cmlen(ipar_gsum)/cmcal(ipar_gsum),cmcal(ipar_gsum)
       WRITE(6,fformat1) ' GLOBAL MULTIPLICATION',&
            cmlen(ipar_gmul)/cmcal(ipar_gmul),cmcal(ipar_gmul)
       WRITE(6,fformat1) ' ALL TO ALL COMM ',&
            cmlen(ipar_aall)/cmcal(ipar_aall),cmcal(ipar_aall)
       WRITE(6,fformat1) ' ALLGATHERV ',&
            cmlen(ipar_agav)/cmcal(ipar_agav),cmcal(ipar_agav)
       !   ==------------------------------------------------------------==
       WRITE(6,'(1X,"=",29X,A11,10X,A10,2X,"=")')&
            'PERFORMANCE','TOTAL TIME'
       WRITE(6,fformat2) ' SEND/RECEIVE',&
            cmlen(ipar_send)/cmtim(ipar_send)*1.e-3_real_8,&
            cmtim(ipar_send)*1.e-3_real_8
       WRITE(6,fformat2) ' BROADCAST',&
            cmlen(ipar_cast)/cmtim(ipar_cast)*1.e-3_real_8,&
            cmtim(ipar_cast)*1.e-3_real_8
       WRITE(6,fformat2) ' GLOBAL SUMMATION',&
            xln*cmlen(ipar_gsum)/cmtim(ipar_gsum)*1.e-3_real_8,&
            cmtim(ipar_gsum)*1.e-3_real_8
       WRITE(6,fformat2) ' GLOBAL MULTIPLICATION',&
            xln*cmlen(ipar_gmul)/cmtim(ipar_gmul)*1.e-3_real_8,&
            cmtim(ipar_gmul)*1.e-3_real_8
       WRITE(6,fformat2) ' ALL TO ALL COMM ',&
            cmlen(ipar_aall)/cmtim(ipar_aall)*1.e-3_real_8,&
            cmtim(ipar_aall)*1.e-3_real_8
       WRITE(6,fformat2) ' ALLGATHERV ',&
            cmlen(ipar_agav)/cmtim(ipar_agav)*1.e-3_real_8,&
            cmtim(ipar_agav)*1.e-3_real_8
       WRITE(6,'(1X,"=",A,17X,A,6X,F10.3,A,2X,"=")')&
            ' SYNCHRONISATION ','      ',cmtim(ipar_sync)*1.e-3_real_8,' SEC'
       WRITE(6,'(1X,64("="))')
       CALL m_flush(6)
#if defined(__ES)
       CLOSE(6)
#endif
    ENDIF
#ifdef __PARALLEL
    !IPHIGENIE/CPMD interface finalizes in IPHIGENIE
    call mp_sync(parai%cp_grp)
    IF (cnti%iftype.NE.3) THEN
       call mp_win_dealloc_shared_mem(parai%node_nproc,parai%node_me,parai%node_grp)
       CALL mpi_comm_free(parai%cp_inter_grp,ierr)
       CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
       CALL mpi_comm_free(parai%node_grp,ierr)
       CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
       CALL mpi_comm_free(parai%cp_grp,ierr)
       CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
       CALL mpi_finalize(ierr)
       CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
       IF (ierr.NE.0) CALL stopgm('MPI_FINALIZE','error in mpi_finalize',& 
            __LINE__,__FILE__)
    ENDIF
#endif
    ! ==--------------------------------------------------------------==
  END SUBROUTINE mp_end


  SUBROUTINE mp_mpi_error_assert(ierr,procN,line,file)
    INTEGER, INTENT(in) :: ierr,line
    CHARACTER(*), INTENT(in) :: file,procN

#ifdef __PARALLEL
    INTEGER :: err_msg_len, ierr_err
    CHARACTER(MPI_MAX_ERROR_STRING) :: err_msg

    IF(ierr/=MPI_SUCCESS) THEN
       CALL mpi_error_string(ierr,err_msg,err_msg_len,ierr_err)
       CALL stopgm(procN,err_msg,line,file)
    ENDIF
#endif

  END SUBROUTINE mp_mpi_error_assert


  FUNCTION mp_is_comm_null(comm)
    !     arguments
#ifdef __PARALLEL
    type(MPI_COMM) ::  comm
#else
    INTEGER :: comm
#endif
    LOGICAL :: mp_is_comm_null
#ifdef __PARALLEL
    !     ==--------------------------------------------------------------==
    mp_is_comm_null=mpi_comm_null==comm
#else
    mp_is_comm_null=.FALSE.
#endif
    !     ==--------------------------------------------------------------==
  END FUNCTION mp_is_comm_null


  SUBROUTINE mp_environ(gid,numtask,taskid)
    ! ==--------------------------------------------------------------==
    ! Arguments
#ifdef __PARALLEL
    INTEGER :: numtask,taskid
    type(MPI_COMM) :: gid
#else
    INTEGER :: gid,numtask,taskid
#endif
    CHARACTER(*),PARAMETER::procedureN='mp_environ'
#ifdef __PARALLEL
    ! Variables
    INTEGER :: ierr
    ! ==--------------------------------------------------------------==
    CALL mpi_comm_rank(gid,taskid,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    CALL mpi_comm_size(gid,numtask,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
    ! ==--------------------------------------------------------------==
    numtask = 1
    taskid  = 0
#endif
    ! ==--------------------------------------------------------------==
  END SUBROUTINE mp_environ


  SUBROUTINE mp_comm_free(comm)
    ! ==--------------------------------------------------------------==
    ! Arguments
#ifdef __PARALLEL
    type(MPI_COMM) :: comm
#else
    INTEGER :: comm
#endif
    ! Variables
    CHARACTER(*),PARAMETER::procedureN='mp_comm_free'
#ifdef __PARALLEL
    ! Variables
    INTEGER :: ierr
    ! ==--------------------------------------------------------------==
    CALL mpi_comm_free(comm,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
    ! ==--------------------------------------------------------------==
  END SUBROUTINE mp_comm_free


  SUBROUTINE mp_comm_dup(comm,new_comm)
    ! ==--------------------------------------------------------------==
    ! Arguments
#ifdef __PARALLEL
    type(MPI_COMM) :: comm,new_comm
#else
    INTEGER :: comm,new_comm
#endif
    ! Variables
    CHARACTER(*),PARAMETER::procedureN='mp_comm_dup'
#ifdef __PARALLEL
    ! Variables
    INTEGER :: ierr
    ! ==--------------------------------------------------------------==
    CALL mpi_comm_dup(comm,new_comm,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
    new_comm=comm
#endif
    ! ==--------------------------------------------------------------==
  END SUBROUTINE mp_comm_dup


  SUBROUTINE mp_comm_compare(comm1,comm2)
    ! ==--------------------------------------------------------------==
    ! Arguments
#ifdef __PARALLEL
    type(MPI_COMM) :: comm1,comm2
#else
    INTEGER :: comm1,comm2
#endif
    ! Variables
    CHARACTER(*),PARAMETER::procedureN='mp_comm_compare'
#ifdef __PARALLEL
    ! Variables
    INTEGER :: ierr,RESULT
    ! ==--------------------------------------------------------------==
    CALL mpi_comm_compare(comm1,comm2,RESULT,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    SELECT CASE(RESULT)
    CASE(MPI_IDENT); WRITE(6,*) 'identical'
    CASE(MPI_CONGRUENT); WRITE(6,*) 'congruent'
    CASE(MPI_SIMILAR); WRITE(6,*) 'similar'
    CASE(MPI_UNEQUAL); WRITE(6,*) 'unequal'
    CASE DEFAULT; WRITE(6,*) 'unknown'
    END SELECT
#endif
    ! ==--------------------------------------------------------------==
  END SUBROUTINE mp_comm_compare


  SUBROUTINE mp_comm_info(comm)
    ! ==--------------------------------------------------------------==
    ! Arguments
#ifdef __PARALLEL
    type(MPI_COMM) :: comm
#else
    INTEGER :: comm
#endif
    ! Variables
    INTEGER :: numtask,taskid
    CHARACTER(*),PARAMETER::procedureN='mp_comm_info'
#ifdef __PARALLEL
    ! Variables
    INTEGER :: ierr
    ! ==--------------------------------------------------------------==
    CALL mpi_comm_rank(comm,taskid,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    CALL mpi_comm_size(comm,numtask,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
    ! ==--------------------------------------------------------------==
    numtask = 1
    taskid  = 0
#endif
    WRITE(6,*) ' TASKID=',taskid,' NUMTASK=',numtasK
    ! ==--------------------------------------------------------------==
  END SUBROUTINE mp_comm_info


  SUBROUTINE mp_get_processor_name ( name )
    ! ==--------------------------------------------------------------==
    ! Arguments
    CHARACTER(*), INTENT(OUT) :: name
    CHARACTER(*),PARAMETER::procedureN='mp_get_processor_name'
#ifdef __PARALLEL
    ! Variables
    CHARACTER( MPI_MAX_PROCESSOR_NAME) :: name_tmp
    INTEGER :: ierr, resultlen
    ! ==--------------------------------------------------------------==
    CALL mpi_get_processor_name ( name_tmp , resultlen, ierr )
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    name = name_tmp(1:resultlen)
#else
    ! ==--------------------------------------------------------------==
    name='serial'
#endif
    ! ==--------------------------------------------------------------==
  END SUBROUTINE mp_get_processor_name


  SUBROUTINE mp_task_query(cp_grp)
    ! ==--------------------------------------------------------------==
    ! Arguments
#ifdef __PARALLEL
    type(MPI_COMM) :: cp_grp
#else
    INTEGER :: cp_grp
#endif
    CHARACTER(*),PARAMETER::procedureN='mp_task_query'
#ifdef __PARALLEL
    ! Variables
    INTEGER :: ierr
    ! ==--------------------------------------------------------------==
    CALL mpi_comm_dup(mp_comm_world,cp_grp,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
    ! ==--------------------------------------------------------------==
    cp_grp=0
#endif
    ! ==--------------------------------------------------------------==
  END SUBROUTINE mp_task_query


  SUBROUTINE mp_cart(comm,nogrp,npgrp,meogrp,mepgrp)
    ! ==--------------------------------------------------------------==
    ! == Forms cartesian groups of processors                         ==
    ! ==--------------------------------------------------------------==
#ifdef __PARALLEL
    INTEGER                                  :: nogrp, npgrp
    type(MPI_COMM)                           :: comm, meogrp, mepgrp
#else
    INTEGER                                  :: comm, nogrp, npgrp, &
                                                meogrp, mepgrp
#endif
    CHARACTER(*), PARAMETER                  :: procedureN = 'mp_cart'

#ifdef __PARALLEL
    type(MPI_GROUP)                          :: world, ogr, pgr
#else
    INTEGER                                  :: world, ogr, pgr
#endif

#ifdef __PARALLEL
    ! Variables
#ifdef __PARALLEL
    INTEGER :: i,io,np(2),ierr,ioo,nnodes
    type(MPI_COMM) :: cart_grp
#else
    INTEGER :: cart_grp,i,io,np(2),ierr,ioo,nnodes
#endif
    LOGICAL :: period(2),order,s_col(2),s_row(2)
    INTEGER :: err_msg_len, ierr_err
    ! ==--------------------------------------------------------------==
    order=.TRUE.
    period(1)=.FALSE.
    period(2)=.FALSE.
    s_col(2)=.TRUE.
    s_col(1)=.FALSE.
    s_row(1)=.TRUE.
    s_row(2)=.FALSE.
    np(1)=npgrp
    np(2)=nogrp

    ! if needed: get the 2d grid
    CALL mpi_comm_size ( comm, nnodes, ierr )
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    IF ( ANY( np == 0 ) ) THEN
       CALL mpi_dims_create ( nnodes, 2, np, ierr )
       CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
       npgrp=np(1)
       nogrp=np(2)
    ENDIF

    CALL mpi_cart_create(comm,2,np,period,order,&
         CART_GRP,IERR)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    CALL mpi_cart_sub(cart_grp,s_row,mepgrp,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    CALL mpi_cart_sub(cart_grp,s_col,meogrp,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    CALL mpi_comm_group(comm,world,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    CALL mpi_comm_group(mepgrp,pgr,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    CALL mpi_comm_group(meogrp,ogr,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
    ! ==--------------------------------------------------------------==
    meogrp=0
    mepgrp=0
    IF (nogrp==0) nogrp=1
    IF (npgrp==0) npgrp=1
#endif
  END SUBROUTINE mp_cart

  SUBROUTINE mp_set_global_consts()
    IMPLICIT NONE
#ifdef __PARALLEL
    mp_max_processor_name = MPI_MAX_PROCESSOR_NAME
    mp_comm_null = MPI_COMM_NULL
#else
    mp_max_processor_name = 256
    mp_comm_null = -1
#endif
  END SUBROUTINE mp_set_global_consts

  SUBROUTINE mp_set_atomic_type_sizes()
    ! ==--------------------------------------------------------------==
    IMPLICIT NONE
#ifdef __PARALLEL
    ! Variables
    INTEGER :: ierr
    ! ==--------------------------------------------------------------==
    CALL mpi_type_size(mpi_integer1,         mp_int1_in_bytes,           ierr)
    CALL mpi_type_size(mpi_integer2,         mp_int2_in_bytes,           ierr)
    CALL mpi_type_size(mpi_integer,          mp_int4_in_bytes,           ierr)
    CALL mpi_type_size(mpi_integer8,         mp_int8_in_bytes,           ierr)
    CALL mpi_type_size(mpi_real,             mp_real_in_bytes,           ierr)
    CALL mpi_type_size(mpi_double_precision, mp_double_in_bytes,         ierr)
    CALL mpi_type_size(mpi_complex,          mp_complex_in_bytes,        ierr)
    CALL mpi_type_size(mpi_double_complex,   mp_double_complex_in_bytes, ierr)
    CALL mpi_type_size(mpi_character,        mp_char_in_bytes,           ierr)
    CALL mpi_type_size(mpi_logical,          mp_logical_in_bytes,        ierr)
#else
    ! so far we dont care those sizes for serial compilation
    mp_int1_in_bytes=0
    mp_int2_in_bytes=0
    mp_int4_in_bytes=0
    mp_int8_in_bytes=0
    mp_real_in_bytes=0
    mp_double_in_bytes=0
    mp_complex_in_bytes=0
    mp_double_complex_in_bytes=0
    mp_char_in_bytes=0
    mp_logical_in_bytes=0
#endif
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mp_set_atomic_type_sizes


  SUBROUTINE mp_sync(gid)
    ! ==--------------------------------------------------------------==
    ! == SYNCHRONISATION OF ALL PROCESSORS                            ==
    ! ==--------------------------------------------------------------==
#ifdef __PARALLEL
    type(MPI_COMM)                           :: gid
#else
    INTEGER                                  :: gid
#endif

    CHARACTER(*), PARAMETER                  :: procedureN = 'mp_sync'

#ifdef __PARALLEL
    ! Variables
    INTEGER :: ierr
    REAL(real_8) :: tim1,tim2
    ! ==--------------------------------------------------------------==
    tim1=m_walltime()
    CALL mpi_barrier(gid,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
    ! Variables
    REAL(real_8) :: tim1,tim2
    ! ==--------------------------------------------------------------==
    tim1=m_walltime()
#endif
    tim2=m_walltime()
    cmtim(ipar_sync)=cmtim(ipar_sync)+tim2-tim1
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mp_sync



  ! ==================================================================
  SUBROUTINE mp_sum_in_place_int1_r1(DATA,n,comm)
    ! ==--------------------------------------------------------------==
    ! == Wrapper to mpi_allreduce (summation)                         ==
    ! ==--------------------------------------------------------------==
    INTEGER(int_1), INTENT(INOUT), TARGET    :: DATA(*)
#ifdef __PARALLEL
    INTEGER, INTENT(IN)                      :: n
    type(MPI_COMM), INTENT(IN)               :: comm
#else
    INTEGER, INTENT(IN)                      :: n, comm
#endif

    CHARACTER(*), PARAMETER :: procedureN = 'mp_sum_in_place_int1_r1'

    ! Variables
    INTEGER(int_1), POINTER __CONTIGUOUS :: para_buff_user(:), data_ptr(:)
    TYPE(c_ptr) :: loc_x

#ifdef __PARALLEL
    INTEGER :: ierr,stat,i
    REAL(real_8) :: tim1,tim2
    TYPE(MPI_Op), SAVE :: op_sum_i1
    LOGICAL, SAVE      :: init=.FALSE.
    ! ==--------------------------------------------------------------==
    ! INT*1 not always supported, so we define an op
    IF(.NOT.init)THEN
       op_sum_i1 = mpi_op_null
       CALL mpi_op_create(sum_i1,.TRUE.,op_sum_i1,ierr)
       init=.TRUE.
    END IF
    cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
    cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_int1_in_bytes
    tim1=m_walltime()
    IF(para_use_mpi_in_place)THEN
       CALL mpi_allreduce(MPI_IN_PLACE,data,n,mpi_integer1,&
            op_sum_i1,comm,ierr)
       CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    ELSE
#if !defined _USE_SCRATCHLIBRARY
       IF(.NOT.ALLOCATED(para_buff))THEN
          ALLOCATE(para_buff(il_para_buff(1)),stat=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate para_buff', &
               __LINE__,__FILE__)
       END IF
#endif
       IF(n*mp_int1_in_bytes.GT.buff_size_in_bytes*il_para_buff(1))THEN
          il_para_buff(1)=CEILING(REAL(n*mp_int1_in_bytes,real_8)/&
               buff_size_in_bytes)
#if !defined _USE_SCRATCHLIBRARY
          DEALLOCATE(para_buff,stat=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
               __LINE__,__FILE__)
          ALLOCATE(para_buff(il_para_buff(1)),stat=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate para_buff', &
               __LINE__,__FILE__)
#endif
       END IF
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate para_buff', &
            __LINE__,__FILE__)
#endif
       CALL mpi_allreduce(data,para_buff,n,mpi_integer1,&
            op_sum_i1,comm,ierr)
       CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
            __LINE__,__FILE__)
#endif
       loc_x = C_LOC(para_buff)
       CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
       loc_x = C_LOC(data)
       CALL C_F_POINTER(loc_x, data_ptr, (/n/))
       !$omp parallel do simd
       do i=1,n
          data_ptr(i)=para_buff_user(i)
       end do
    END IF
    tim2=m_walltime()
    cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
  END SUBROUTINE mp_sum_in_place_int1_r1

  SUBROUTINE sum_i1(invec,inoutvec,len,TYPE)
    USE, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer
    INTEGER                                  :: len
    type(C_PTR), VALUE                       :: inoutvec, invec
    type(MPI_Datatype)                       :: TYPE
    INTEGER(int_1), POINTER                  :: invec_f(:), inoutvec_f(:)

    INTEGER                                  :: i
    
    CALL c_f_pointer(invec,invec_f,(/ len /))
    CALL c_f_pointer(inoutvec,inoutvec_f,(/ len /))
    CALL sumi1(invec_f,inoutvec_f,len)
  END SUBROUTINE sum_i1
  
  SUBROUTINE sumi1(in,inout,len)
    INTEGER, INTENT(IN) :: len
    INTEGER(int_1), INTENT(IN) :: in(*)
    INTEGER(int_1), INTENT(INOUT) :: inout(*)
    INTEGER :: i

    DO i=1,len
       inout(i)=in(i)+inout(i)
    END DO

  END SUBROUTINE sumi1

  SUBROUTINE mp_dims_create(nodes,ndims,dims)
    ! ==--------------------------------------------------------------==
    ! == Wrapper to mpi_dims_create                                   ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nodes, ndims, dims(ndims)

    CHARACTER(*), PARAMETER                  :: procedureN = 'mp_dims_create'

! Variables

#ifdef __PARALLEL
    INTEGER :: ierr
    INTEGER :: err_msg_len, ierr_err
    ! ==--------------------------------------------------------------==
    CALL mpi_dims_create(nodes,ndims,dims,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
    dims = 1
#endif
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mp_dims_create


  SUBROUTINE mp_split(oldcomm,color,key,gid,newrank)
    ! ==--------------------------------------------------------------==
    ! == Splits processors into groups                                ==
    ! ==--------------------------------------------------------------==
#ifdef __PARALLEL
    INTEGER                                  :: color, key, newrank
    type(MPI_COMM)                           :: oldcomm, gid
#else
    INTEGER                                  :: oldcomm, color, key, gid, &
                                                newrank
#endif
    CHARACTER(*), PARAMETER                  :: procedureN = 'mp_split'

#ifdef __PARALLEL
    ! Variables
    INTEGER :: ierr
    CHARACTER (len=60) :: strout
    ! ==--------------------------------------------------------------==
    CALL mpi_comm_split(oldcomm,color,key,gid,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)

    CALL mpi_comm_rank(gid,newrank,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
    ! ==--------------------------------------------------------------==
    gid=0
    newrank=0
#endif
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mp_split

  SUBROUTINE mp_split_type(key,comm,newcomm)
#ifdef __PARALLEL
    INTEGER,INTENT(IN)                       :: key
    type(MPI_COMM),INTENT(IN)                :: comm
    type(MPI_COMM),INTENT(OUT)               :: newcomm
#else
    INTEGER,INTENT(IN)                       :: key,comm
    INTEGER,INTENT(OUT)                      :: newcomm
#endif
    INTEGER                                  :: ierr
    CHARACTER(*), PARAMETER                  :: procedureN = 'mp_comm_split_type'
    ! Author: Tobias Kloeffel, FAU Erlangen Nuernberg, March 2019

#ifdef __PARALLEL
    CALL mpi_comm_split_type(comm,MPI_COMM_TYPE_SHARED,key,MPI_INFO_NULL,newcomm,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
    newcomm=comm
#endif
  END SUBROUTINE mp_split_type


  SUBROUTINE mp_group(gsize,glist,gid,pcomm)
    ! ==--------------------------------------------------------------==
    ! == Forms groups of processors                                   ==
    ! ==--------------------------------------------------------------==
#ifdef __PARALLEL
    INTEGER                                  :: gsize, glist(*)
    type(MPI_COMM)                           :: gid, pcomm
#else
    INTEGER                                  :: gsize, glist(*), gid, pcomm
#endif
    CHARACTER(*), PARAMETER                  :: procedureN = 'mp_group'

#ifdef __PARALLEL
    ! Variables
#ifdef __PARALLEL
    type(MPI_GROUP) :: world,newgroup
#else
    INTEGER :: world,newgroup
#endif
    INTEGER :: ierr
    ! ==--------------------------------------------------------------==
    CALL mpi_comm_group(pcomm,world,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    CALL mpi_group_incl(world,gsize,glist,newgroup,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    CALL mpi_comm_create(pcomm,newgroup,gid,ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
    ! ==--------------------------------------------------------------==
    gid=0
#endif
    ! ==--------------------------------------------------------------==
  END SUBROUTINE mp_group


  SUBROUTINE mp_get_node_env ( comm, node_numtasks, node_taskid )
#ifdef __PARALLEL
    type(MPI_COMM), INTENT(IN) :: comm
#else
    INTEGER, INTENT(IN) :: comm
#endif
    INTEGER, INTENT(OUT) :: node_numtasks, node_taskid
    CHARACTER(*), PARAMETER                  :: procedureN = 'mp_get_node_env'
#ifdef __PARALLEL
    INTEGER :: numtasks , ierr, namelength, taskid, itask, stat
    CHARACTER(MPI_MAX_PROCESSOR_NAME), ALLOCATABLE, DIMENSION(:) :: names
    CHARACTER(MPI_MAX_PROCESSOR_NAME) :: name
#else

#endif

#ifdef __PARALLEL

    CALL mpi_comm_size( comm , numtasks, ierr )
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    CALL mpi_comm_rank( comm, taskid, ierr )
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)

    ALLOCATE( names(0:numtasks-1), stat=stat )
    IF (stat.NE.0) CALL stopgm(procedureN,'Allocation problem',&
         __LINE__,__FILE__)
    names(:)=''
    name=''

    CALL mpi_get_processor_name( name, namelength, ierr )
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)

    CALL mpi_allgather( name, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
         names, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, comm, ierr )
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)

    ! count how many tasks own the same processor name
    node_numtasks = 0
    DO itask = 0, numtasks-1
       IF( itask == taskid ) node_taskid = node_numtasks
       IF( name == names( itask ) ) node_numtasks = node_numtasks + 1
    ENDDO

    DEALLOCATE( names, stat=stat )
    IF (stat.NE.0) CALL stopgm(procedureN,'Deallocation problem',&
         __LINE__,__FILE__)

#else
    node_numtasks = 1
    node_taskid = 0

#endif

  END SUBROUTINE mp_get_node_env


  SUBROUTINE mp_get_version(version, subversion)
    INTEGER, INTENT(OUT) :: version, subversion
    CHARACTER(*), PARAMETER                  :: procedureN = 'mp_get_version'

#ifdef __PARALLEL
    INTEGER :: ierr
#else

#endif

#ifdef __PARALLEL
    CALL mpi_get_version(version, subversion, ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
    version = 0
    subversion = 0
#endif
  END SUBROUTINE mp_get_version


  SUBROUTINE mp_get_library_version( version )
    CHARACTER(*), INTENT(OUT) :: version
    CHARACTER(*), PARAMETER                  :: procedureN = 'mp_get_library_version'

#if defined(__PARALLEL) && ! defined(_HASNT_MPI_30)
    INTEGER :: ierr, resulten
    CHARACTER(MPI_MAX_LIBRARY_VERSION_STRING) :: v
#endif

#if defined(__PARALLEL) && ! defined(_HASNT_MPI_30)

    CALL mpi_get_library_version ( v, resulten, ierr )
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)

    version = v(1:MIN(resulten,LEN(version)))

#else

    version = 'not available'

#endif

  END SUBROUTINE mp_get_library_version

  SUBROUTINE mp_win_lock_all_shared(comm)
    ! ==--------------------------------------------------------------==
    ! == Wrapper for mpi_win_fence                                    ==
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, FAU Erlangen Nuernberg, March 2019
#ifdef __PARALLEL
    type(MPI_COMM),INTENT(IN) :: comm
#else
    INTEGER,INTENT(IN) :: comm
#endif
#ifdef __PARALLEL
    INTEGER:: index,ierr
    CHARACTER(*),PARAMETER::procedureN='mp_win_sync'

    IF(comm.EQ.parai%node_grp)THEN
       index=1
    ELSEIF(comm.EQ.parai%cp_inter_node_grp)THEN
       index=2
    ELSE
       CALL stopgm(procedureN,'Unsupported mpi communicator',&
            __LINE__,__FILE__)
    END IF
    CALL mpi_win_lock_all(MPI_MODE_NOCHECK,mpi_window(index),ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  END SUBROUTINE mp_win_lock_all_shared

  SUBROUTINE mp_win_unlock_all_shared(comm)
    ! ==--------------------------------------------------------------==
    ! == Wrapper for mpi_win_fence                                    ==
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, FAU Erlangen Nuernberg, March 2019
#ifdef __PARALLEL
    type(MPI_COMM),INTENT(IN) :: comm
#else
    INTEGER,INTENT(IN) :: comm
#endif
#ifdef __PARALLEL
    INTEGER:: index,ierr
    CHARACTER(*),PARAMETER::procedureN='mp_win_sync'

    IF(comm.EQ.parai%node_grp)THEN
       index=1
    ELSEIF(comm.EQ.parai%cp_inter_node_grp)THEN
       index=2
    ELSE
       CALL stopgm(procedureN,'Unsupported mpi communicator',&
            __LINE__,__FILE__)
    END IF
    CALL mpi_win_unlock_all(mpi_window(index),ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  END SUBROUTINE mp_win_unlock_all_shared

  SUBROUTINE mp_win_sync(comm)
    ! ==--------------------------------------------------------------==
    ! == Wrapper for mpi_win_fence                                    ==
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, FAU Erlangen Nuernberg, March 2019
#ifdef __PARALLEL
    type(MPI_COMM),INTENT(IN) :: comm
#else
    INTEGER,INTENT(IN) :: comm
#endif
#ifdef __PARALLEL
    INTEGER:: index,ierr
    CHARACTER(*),PARAMETER::procedureN='mp_win_sync'

    IF(comm.EQ.parai%node_grp)THEN
       index=1
    ELSEIF(comm.EQ.parai%cp_inter_node_grp)THEN
       index=2
    ELSE
       CALL stopgm(procedureN,'Unsupported mpi communicator',&
            __LINE__,__FILE__)
    END IF
    CALL mpi_win_sync(mpi_window(index),ierr)
    CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    CALL mp_sync(comm)
#endif
  END SUBROUTINE mp_win_sync

  
  SUBROUTINE mp_win_alloc_shared_mem(type,lda,n,baseptr,nproc,mypos,comm)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
    ! ==--------------------------------------------------------------==
    ! == Return baseptr to shared memory window                       ==
    ! == Currently two shared memory windows are maintained, one for  ==
    ! == comm.eq.node_grp (index.eq.1) and one for                    ==
    ! == comm.eq.cp_inter_node_grp (index.eq.2)                       ==
    ! == Deallocation for type .eq. A or a, reallocation possible     ==
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, FAU Erlangen Nuernberg, March 2019
#ifdef __PARALLEL
    INTEGER,INTENT(IN) :: lda,n,nproc,mypos
    type(MPI_COMM),INTENT(IN) :: comm
#else
    INTEGER,INTENT(IN) :: lda,n,nproc,mypos,comm
#endif
    CHARACTER(1),INTENT(IN) :: type
    TYPE(C_PTR),INTENT(OUT) :: baseptr(0:nproc-1)
#ifdef __PARALLEL
    INTEGER(int_8), SAVE :: allocated_size(2) = 0
    INTEGER(int_8) :: requested_size
    INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
    type(MPI_Info) :: info
    INTEGER :: proc,ierr,index,displ
    LOGICAL, SAVE :: alloc(2)=.FALSE.
    CHARACTER(*),PARAMETER::procedureN='mp_win_alloc_shared_mem'

    IF(type.EQ.'C'.OR.type.EQ.'c')THEN
       requested_size=lda*n*mp_double_complex_in_bytes
    ELSEIF(type.EQ.'R'.OR.type.EQ.'r')THEN
       requested_size=lda*n*mp_double_in_bytes
    ELSEIF(type.EQ.'D'.OR.type.EQ.'d')THEN
       DO index=1,size(alloc)
          IF(alloc(index))THEN
             CALL MPI_WIN_FREE(mpi_window(index), IERR)
             CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
             alloc(index)=.FALSE.
             allocated_size(index)=0
          END IF
       END DO
       RETURN
    END IF

    IF(comm.EQ.parai%node_grp)THEN
       index=1
    ELSEIF(comm.EQ.parai%cp_inter_node_grp)THEN
       index=2
    ELSE
       CALL stopgm(procedureN,'Unsupported mpi communicator',&
            __LINE__,__FILE__)
    END IF

    CALL mp_max(requested_size,comm)

    !check if allocated window is big enough...
    IF (allocated_size(index).LT.requested_size) THEN
       IF(alloc(index)) THEN
          CALL MPI_WIN_FREE(mpi_window(index), IERR)
          CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
          alloc(index)=.FALSE.
          allocated_size(index)=0
       END IF
       CALL MPI_INFO_CREATE(info,ierr)
       CALL MPI_INFO_SET(INFO, 'same_disp_unit', 'true',ierr)
       CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
       CALL MPI_INFO_SET(INFO, 'alloc_shared_noncontig', 'true',ierr)
       CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
       CALL mpi_info_set(info, 'same_size', 'true',ierr)
       CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
       windowsize=requested_size
       allocated_size(index)=requested_size
       displ=8
       CALL MPI_WIN_allocate_shared(windowsize , displ, info,  &
            comm, baseptr(mypos), mpi_window(index), ierr)
       CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
       CALL mpi_info_free(info,ierr)
       CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
       alloc(index)=.TRUE.
    END IF
    !set baseptrs for all procs
    DO proc=0,nproc-1
       CALL mpi_win_shared_query(mpi_window(index), proc, windowsize, &
            displ, baseptr(proc), ierr)
       CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
    END DO
#endif
  END SUBROUTINE mp_win_alloc_shared_mem

    SUBROUTINE mp_win_dealloc_shared_mem(nproc,mypos,comm)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
    ! ==--------------------------------------------------------------==
    ! == Return baseptr to shared memory window                       ==
    ! == Currently two shared memory windows are maintained, one for  ==
    ! == comm.eq.node_grp (index.eq.1) and one for                    ==
    ! == comm.eq.cp_inter_node_grp (index.eq.2)                       ==
    ! == Deallocation for type .eq. A or a, reallocation possible     ==
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, FAU Erlangen Nuernberg, March 2019
#ifdef __PARALLEL
    INTEGER,INTENT(IN) :: nproc,mypos
    type(MPI_COMM),INTENT(IN) :: comm
#else
    INTEGER,INTENT(IN) :: nproc,mypos,comm
#endif
    TYPE(C_PTR) :: baseptr(0:nproc-1)
#ifdef __PARALLEL
    CHARACTER(*),PARAMETER::procedureN='mp_win_dealloc_shared_mem'

    CALL mp_win_alloc_shared_mem('d',1,1,baseptr,nproc,mypos,comm)
#endif
  END SUBROUTINE mp_win_dealloc_shared_mem

  !
  ! include file for the interfaces
  !
#include "bcast.inc"
#include "bcast_r0.inc"
#include "all2all.inc"

#include "send_and_recv.inc"

#include "sendrecv.inc"

#include "sum.inc"
#include "sum_r0.inc"
#include "sum_root.inc"
#include "sum_root_in_place.inc"
#include "sum_in_place.inc"
#include "sum_scalar_in_place.inc"

#include "prod.inc"
#include "prod_r0.inc"

#include "max_scalar.inc"
#include "max.inc"
#include "min_scalar.inc"
#include "min.inc"


END MODULE mp_interface
