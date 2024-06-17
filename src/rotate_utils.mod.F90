#include "cpmd_global.h"

MODULE rotate_utils
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes,&
                                             cp_grp_redist_array_f
  USE distribution_utils,              ONLY: dist_entity
  USE error_handling,                  ONLY: stopgm
  USE gpu
  USE kinds,                           ONLY: real_8,&
                                             int_8,&
                                             int_4
  USE mp_interface,                    ONLY: mp_win_alloc_shared_mem,&
                                             mp_win_sync,&
                                             mp_sync,&
                                             mp_win_lock_all_shared,&
                                             mp_win_unlock_all_shared
  USE nvtx_utils
  !$ USE omp_lib,                      ONLY: omp_set_nested,&
  !$                                         omp_get_thread_num,&
  !$                                         omp_set_num_threads
  USE parac,                           ONLY: parai
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE, INTRINSIC :: ISO_C_BINDING,     ONLY: C_PTR,&
                                             C_F_POINTER
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif
#ifdef _INTEL_MKL
  use mkl_service
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rotate
  !!public :: rotate_c
  !!public :: rotate_da
  PUBLIC :: rottr
  PUBLIC :: rottr_fnl
  PUBLIC :: rotate_fnl
  PUBLIC :: rottr_c0_fnl
  PUBLIC :: rotate_c0_fnl

CONTAINS
  ! ==================================================================
  SUBROUTINE rotate_fnl(ldf,fnl_p,fnlgam_p,nstate,gam)
    ! ==--------------------------------------------------------------==
    ! ==         FNLGAM <= A*FNL*GAM                                  ==
    ! == USPP needs fnl*gam two times. This operation is rather cheap ==
    ! == and is only parallelizable on the node level. Hence we create==
    ! == a mpi shared memory window and perform the dgemm directly in ==
    ! == this window.                                                 ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate, ldf
    REAL(real_8),INTENT(IN)                  :: gam(nstate,*), fnl_p(ldf,nstate)
    REAL(real_8),INTENT(OUT), TARGET         :: fnlgam_p(ldf,nstate)
    INTEGER                                  :: isub, loc_work, start_work,&
                                                i, arrayshape(2), proc, nmin,&
                                                nchunk,nmin_2,nchunk_2, lda,&
                                                work_proc(2,0:parai%node_nproc-1),&
                                                len_proc(0:parai%node_nproc-1),&
                                                ierr
    TYPE(C_PTR)                              :: baseptr(0:parai%node_nproc-1)
    CHARACTER(*), PARAMETER                  :: procedureN = 'rotate_fnl'
    REAL(real_8), POINTER __CONTIGUOUS       :: loc(:,:)
    TYPE CONTAINER
       REAL(real_8),  POINTER __CONTIGUOUS   :: temp(:,:)
    END TYPE CONTAINER
    TYPE(CONTAINER),ALLOCATABLE :: iproc(:)
! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)

    IF(parai%node_nproc.GT.1)THEN
       !parallelize dsymms at node level
       CALL dist_entity(ldf,parai%node_nproc,work_proc)
       DO proc=0,parai%node_nproc-1
          len_proc(proc)=work_proc(2,proc)-work_proc(1,proc)+1
       END DO
       loc_work=len_proc(parai%node_me)
       start_work=work_proc(1,parai%node_me)
       lda=MAXVAL(len_proc)
       ALLOCATE(iproc(0:parai%node_nproc-1), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate iproc',&
            __LINE__,__FILE__)

       CALL mp_win_alloc_shared_mem('r',lda,nstate,baseptr,parai%node_nproc,&
            parai%node_me,parai%node_grp)
       CALL mp_win_lock_all_shared(parai%node_grp)
       arrayshape(1)=lda
       arrayshape(2)=nstate
       DO proc=0,parai%node_nproc-1
          CALL C_F_POINTER(baseptr(proc), iproc(proc)%temp,arrayshape)
       END DO
       loc=>iproc(parai%node_me)%temp(:,:)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target enter data map(alloc:loc)
#endif
    ELSE
       loc_work=ldf
       lda=ldf
       start_work=1
       loc=>fnlgam_p
    END IF

    !set up spin settings
    IF (cntl%tlsd) THEN
       nmin=1
       nchunk=spin_mod%nsup
       nmin_2=spin_mod%nsup+1
       nchunk_2=spin_mod%nsdown
    ELSE
       nmin=1
       nchunk=nstate
       nmin_2=0
       nchunk_2=0
    ENDIF

    CALL cpmd_dsymm('R','U',loc_work,nchunk,1.0_real_8,&
         gam(nmin,nmin),nstate,&
         fnl_p(start_work,nmin),ldf,&
         0.0_real_8,loc(1,nmin),lda)
    IF(cntl%tlsd)THEN
       CALL cpmd_dsymm('R','U',loc_work,nchunk_2,1.0_real_8,&
            gam(nmin_2,nmin_2),nstate,&
            fnl_p(start_work,nmin_2),ldf,&
            0.0_real_8,loc(1,nmin_2),lda)
    END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target update from(loc)
#endif
    IF(parai%node_nproc.GT.1)THEN
       !sync shared memory window
       CALL mp_win_sync(parai%node_grp)
       !get data from all procs
       !$omp parallel do private(proc,i)
       DO i=1,nstate
          DO proc=0,parai%node_nproc-1
             CALL copy(fnlgam_p(work_proc(1,proc):work_proc(2,proc),i),&
                  iproc(proc)%temp(:,i),len_proc(proc))
          END DO
       END DO
       DEALLOCATE(iproc, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate iproc',&
            __LINE__,__FILE__)
       CALL mp_win_unlock_all_shared(parai%node_grp)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target update to(fnlgam_p)
       !$omp target exit data map(delete:loc)
#endif
    END IF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rotate_fnl
  ! ==================================================================
  SUBROUTINE rottr_fnl(ldf,fnl_p,nstate,gam)
    ! ==--------------------------------------------------------------==
    ! ==         FNL <= A*FNL*GAM                                     ==
    ! ==   SPECIAL CASE FOR GAM UPPER TRIAGONAL                       ==
    ! ==   Usage: rgsvan                                              ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: ldf, nstate
    REAL(real_8),INTENT(INOUT)               :: fnl_p(ldf,*)
    REAL(real_8),INTENT(IN)                  :: gam(nstate,nstate)

    INTEGER                                  :: isub, loc_work, start_work,&
                                                i, arrayshape(2), proc, nmin,&
                                                nchunk,nmin_2,nchunk_2, lda,&
                                                work_proc(2,0:parai%node_nproc-1),&
                                                len_proc(0:parai%node_nproc-1),&
                                                ierr
    CHARACTER(*), PARAMETER                  :: procedureN = 'rottr_fnl'
    TYPE(C_PTR)                              :: baseptr(0:parai%node_nproc-1)
    TYPE CONTAINER
       REAL(real_8),  POINTER __CONTIGUOUS   :: temp(:,:)
    END TYPE CONTAINER
    TYPE(CONTAINER),ALLOCATABLE :: iproc(:)

! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    IF(parai%node_nproc.GT.1)THEN
       !parallelize dtrmms at node level
       CALL dist_entity(ldf,parai%node_nproc,work_proc)
       DO proc=0,parai%node_nproc-1
          len_proc(proc)=work_proc(2,proc)-work_proc(1,proc)+1
       END DO
       loc_work=len_proc(parai%node_me)
       start_work=work_proc(1,parai%node_me)
       lda=MAXVAL(len_proc)
       ALLOCATE(iproc(0:parai%node_nproc-1), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate iproc',&
            __LINE__,__FILE__)

       CALL mp_win_alloc_shared_mem('r',lda,nstate,baseptr,parai%node_nproc,&
            parai%node_me,parai%node_grp)
       CALL mp_win_lock_all_shared(parai%node_grp)
       arrayshape(1)=lda
       arrayshape(2)=nstate
       DO proc=0,parai%node_nproc-1
          CALL C_F_POINTER(baseptr(proc), iproc(proc)%temp,arrayshape)
       END DO
    ELSE
       loc_work=ldf
       lda=ldf
       start_work=1
    END IF

    !set up spin settings
    IF (cntl%tlsd) THEN
       nmin=1
       nchunk=spin_mod%nsup
       nmin_2=spin_mod%nsup+1
       nchunk_2=spin_mod%nsdown
    ELSE
       nmin=1
       nchunk=nstate
       nmin_2=0
       nchunk_2=0
    ENDIF
    
    IF(loc_work.GT.0)THEN
       CALL cpmd_dtrmm('R','U','N','N',loc_work,nchunk,1.0_real_8,gam(nmin,nmin),nstate, &
            fnl_p(start_work,nmin),ldf)
       IF(cntl%tlsd)THEN
          CALL cpmd_dtrmm('R','U','N','N',loc_work,nchunk_2,1.0_real_8,gam(nmin_2,nmin_2),nstate,&
               fnl_p(start_work,nmin_2),ldf)
       END IF
    END IF
    
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target update from(fnl_p(:,:nstate)) if(.NOT.update_result_to_host.OR.(parai%node_nproc.EQ.1))
#endif
    
    IF(parai%node_nproc.GT.1)THEN
       !copy local chunk into shared mem
       proc=parai%node_me
       !$omp parallel do private (i)
       DO i=1,nstate
          CALL copy(iproc(proc)%temp(:,i),fnl_p(work_proc(1,proc),i),len_proc(proc))
       END DO
       !sync shared memory window
       CALL mp_win_sync(parai%node_grp)
       !get data from other procs
       !$omp parallel do private (i,proc)
       DO i=1,nstate
          DO proc=0,parai%node_nproc-1
             IF(proc.EQ.parai%node_me)CYCLE
             CALL copy(fnl_p(work_proc(1,proc):work_proc(2,proc),i),&
                  iproc(proc)%temp(:,i),len_proc(proc))
          END DO
       END DO
       CALL mp_win_unlock_all_shared(parai%node_grp)
       DEALLOCATE(iproc, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate iproc',&
            __LINE__,__FILE__)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target update to(fnl_p(:,:nstate))
#endif
    END IF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rottr_fnl
  ! ==================================================================
  PURE SUBROUTINE copy (out,in,len)
    !should be inlined, get rid of pointer attribute
    INTEGER,INTENT(IN)                       :: len
    REAL(real_8),INTENT(IN)                  :: in(len)
    REAL(real_8), INTENT(OUT)                :: out(len)

    out=in
  END SUBROUTINE  copy
  ! ==================================================================
  SUBROUTINE rotate_c0_fnl(ldc,c0,c2,ldf,fnl_p,fnlgam_p,nstate,gam,redist)
    ! ==--------------------------------------------------------------==
    ! ==         FNLGAM <= A*FNL*GAM                                  ==
    ! ==         C2 <= 1.0*C2-1.0*C0*GAM                              ==
    ! == USPP also needs fnl*gam. This operation in rather cheap,     ==
    ! == especially when it is parallelized across all procs. So we   ==
    ! == do the allgather step while doing c2=c2-c0*gam               ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate, ldf, ldc
    REAL(real_8),INTENT(IN)                  :: gam(nstate,*),fnl_p(ldf,nstate)
    REAL(real_8),INTENT(OUT)                 :: fnlgam_p(ldf,nstate)
    COMPLEX(real_8),INTENT(INOUT)            :: c2(ldc,*)
    COMPLEX(real_8),INTENT(IN)               :: c0(ldc,*)
    LOGICAL,INTENT(IN)                       :: redist
    CHARACTER(*), PARAMETER                  :: procedureN = 'rotate_c0_fnl'

    INTEGER                                  :: isub, loc_work, start_work, ierr,&
                                                proc, nmin, nchunk, ibeg_c0, &
                                                nmin_2, nchunk_2, ngw_local, &
                                                work_proc(2,0:parai%nproc-1),&
                                                methread, nthreads, nested_threads, nbmax
    INTEGER(int_8)                           :: il_temp(3)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: temp(:,:,:)
#else
    REAL(real_8), ALLOCATABLE                :: temp(:,:,:)
#endif
    CALL tiset(procedureN,isub)
    !set up spin settings
    IF (cntl%tlsd) THEN
       nmin=1
       nchunk=spin_mod%nsup
       nmin_2=spin_mod%nsup+1
       nchunk_2=spin_mod%nsdown
    ELSE
       nmin=1
       nchunk=nstate
       nmin_2=0
       nchunk_2=0
    ENDIF
    !end spin settings
    !get local fnl chunk
    CALL dist_entity(ldf,parai%nproc,work_proc,nbmax=nbmax,nblocal=loc_work,&
         iloc=parai%me)
    il_temp(1)=INT(nbmax,kind=int_8)
    start_work=work_proc(1,parai%me)
    !end local fnl chunk
    !allocate fnlgam buffer
    proc=parai%me+1
    il_temp(2)=nstate
    il_temp(3)=parai%nproc
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_temp,temp,procedureN//'_temp',ierr)
#else
    ALLOCATE(temp(il_temp(1),il_temp(2),il_temp(3)), stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate temp',&
         __LINE__,__FILE__)
    CALL cpmd_dsymm('R','U',loc_work,nchunk,1.0_real_8,&
         gam(nmin,nmin),nstate,fnl_p(start_work,nmin),ldf,&
         0.0_real_8,temp(1,nmin,proc),INT(il_temp(1),kind=int_4))
    IF(cntl%tlsd)THEN
       CALL cpmd_dsymm('R','U',loc_work,nchunk_2,1.0_real_8,&
            gam(nmin_2,nmin_2),nstate,fnl_p(start_work,nmin_2),ldf,&
            0.0_real_8,temp(1,nmin_2,proc),INT(il_temp(1),kind=int_4))
    END IF
    IF(cntl%overlapp_comm_comp)THEN
       nthreads=MIN(2,parai%ncpus)
       nested_threads=(MAX(parai%ncpus-1,1))
#if !defined(_INTEL_MKL)
       CALL stopgm(procedureN, 'Overlapping communication and computation: Behavior of BLAS &
            routine inside parallel region not checked',&
            __LINE__,__FILE__)
#endif
    ELSE
       nthreads=1
       nested_threads=parai%ncpus
    END IF
    methread=0
    !$omp parallel if(nthreads.eq.2) num_threads(nthreads) &
    !$omp private (methread,ngw_local,ibeg_c0) proc_bind(close)
    !$ methread = omp_get_thread_num()
    !$ IF(methread.EQ.1)THEN
    !$    CALL omp_set_max_active_levels(2)
    !$    CALL omp_set_num_threads(nested_threads)
#if defined(_INTEL_MKL) && !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$    CALL mkl_set_dynamic(0)
    !$    ierr = mkl_set_num_threads_local(nested_threads)
#endif
    !$ END IF
    
    IF(methread.EQ.1.OR.nthreads.EQ.1)THEN
       CALL cp_grp_get_sizes(ngw_l=ngw_local,first_g=ibeg_c0)
       CALL cpmd_dsymm('R','U',ngw_local*2,nchunk,-1.0_real_8,&
            gam(nmin,nmin),nstate,c0(ibeg_c0,nmin),ldc*2,&
            1.0_real_8,c2(ibeg_c0,nmin),ldc*2)
       IF(cntl%tlsd)THEN
          CALL cpmd_dsymm('R','U',ngw_local*2,nchunk_2,-1.0_real_8,&
               gam(nmin_2,nmin_2),nstate,c0(ibeg_c0,nmin_2),ldc*2,&
               1.0_real_8,c2(ibeg_c0,nmin_2),ldc*2)
       END IF
    END IF
    IF(methread.EQ.0.OR.nthreads.EQ.1)THEN
       CALL my_concat_inplace(temp,il_temp(1)*nstate,parai%allgrp)
    END IF

    !$ IF (methread.EQ.1) THEN
    !$    CALL omp_set_max_active_levels(1)
    !$    CALL omp_set_num_threads(parai%ncpus)
#if defined(_INTEL_MKL) && !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$    CALL mkl_set_dynamic(1)
    !$    ierr = mkl_set_num_threads_local(0)
#endif
    !$ END IF

    !$omp end parallel
    
    CALL copy_back(fnlgam_p,temp,work_proc,ldf,nstate,parai%nproc)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_temp,temp,procedureN//'_temp',ierr)
#else
    DEALLOCATE(temp, stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate temp',&
         __LINE__,__FILE__)
    IF(redist)CALL cp_grp_redist_array_f(c2,ldc,nstate)
    CALL tihalt(procedureN,isub)
  END SUBROUTINE rotate_c0_fnl
  ! ==================================================================
  SUBROUTINE rottr_c0_fnl(ldc,c0,ldf,fnl,gam,nstate,redist)
    ! ==--------------------------------------------------------------==
    ! ==         FNL <= FNL*GAM                                       ==
    ! ==          C0 <= C0*GAM                                        ==
    ! == USPP also needs fnl*gam. This operation in rather cheap,     ==
    ! == especially when it is parallelized across all procs. So we   ==
    ! == do the allgather step while doing c0=c0*gam                  ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate, ldf, ldc
    REAL(real_8),INTENT(INOUT)               :: fnl(ldf,nstate)
    REAL(real_8),INTENT(IN)                  :: gam(nstate,*)
    COMPLEX(real_8),INTENT(INOUT)            :: c0(ldc,nstate)
    LOGICAL,INTENT(IN)                       :: redist

    CHARACTER(*), PARAMETER                  :: procedureN = 'rottr_c0_fnl'

    INTEGER                                  :: isub, loc_work, start_work, ierr,&
                                                proc, nmin, nchunk, ibeg_c0, &
                                                nmin_2, nchunk_2, ngw_local, &
                                                work_proc(2,0:parai%nproc-1),&
                                                methread, nthreads, nested_threads, nbmax
    INTEGER(int_8)                           :: il_temp(3)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: temp(:,:,:)
#else
    REAL(real_8), ALLOCATABLE                :: temp(:,:,:)
#endif
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    !set up spin settings
    IF (cntl%tlsd) THEN
       nmin=1
       nchunk=spin_mod%nsup
       nmin_2=spin_mod%nsup+1
       nchunk_2=spin_mod%nsdown
    ELSE
       nmin=1
       nchunk=nstate
       nmin_2=0
       nchunk_2=0
    ENDIF
    !end spin settings
    !get local fnl chunk
    CALL dist_entity(ldf,parai%nproc,work_proc,nbmax=nbmax,nblocal=loc_work,&
         iloc=parai%me)
    il_temp(1)=INT(nbmax,KIND=int_4)
    start_work=work_proc(1,parai%me)
    !end local fnl chunk
    !allocate fnlgam buffer
    proc=parai%me+1
    il_temp(2)=nstate
    il_temp(3)=parai%nproc
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_temp,temp,procedureN//'_temp',ierr)
#else
    ALLOCATE(temp(il_temp(1),il_temp(2),il_temp(3)), stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate temp',&
         __LINE__,__FILE__)
    IF(loc_work.GT.0)THEN
       CALL cpmd_dtrmm('R','U','N','N',loc_work,nchunk,1.0_real_8,gam(nmin,nmin),nstate, &
            fnl(start_work,nmin),ldf)
       IF(cntl%tlsd)THEN
          CALL cpmd_dtrmm('R','U','N','N',loc_work,nchunk_2,1.0_real_8,gam(nmin_2,nmin_2),nstate,&
               fnl(start_work,nmin_2),ldf)
       END IF
    END IF
    CALL copy_in(temp(:,:,parai%me+1),fnl,work_proc(:,parai%me),INT(il_temp(1),KIND=int_4),nstate)
    IF(cntl%overlapp_comm_comp)THEN
       nthreads=MIN(2,parai%ncpus)
       nested_threads=(MAX(parai%ncpus-1,1))
#if !defined(_INTEL_MKL)
       CALL stopgm(procedureN, 'Overlapping communication and computation: Behavior of BLAS &
            routine inside parallel region not checked',&
            __LINE__,__FILE__)
#endif
    ELSE
       nthreads=1
       nested_threads=parai%ncpus
    END IF
    methread=0
    !$omp parallel if(nthreads.eq.2) num_threads(nthreads) &
    !$omp private (methread,ngw_local,ibeg_c0) proc_bind(close)
    !$ methread=omp_get_thread_num()
    !$ IF(methread.EQ.1)THEN
    !$    CALL omp_set_max_active_levels(2)
    !$    CALL omp_set_num_threads(nested_threads)
#if defined(_INTEL_MKL) && !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$    CALL mkl_set_dynamic(0)
    !$    ierr = mkl_set_num_threads_local(nested_threads)
#endif
    !$ END IF

    IF(methread.EQ.1.OR.nthreads.EQ.1)THEN
       CALL cp_grp_get_sizes(ngw_l=ngw_local,first_g=ibeg_c0)
       IF(ngw_local.GT.0)THEN
          CALL cpmd_dtrmm('R','U','N','N',ngw_local*2,nchunk,1.0_real_8,gam(nmin,nmin),nstate,&
               c0(ibeg_c0,nmin),ldc*2)
          IF(cntl%tlsd)THEN
             CALL cpmd_dtrmm('R','U','N','N',ngw_local*2,nchunk_2,1.0_real_8,gam(nmin_2,nmin_2),&
                  nstate,c0(ibeg_c0,nmin_2),ldc*2)
          END IF
       END IF
    END IF
    IF(methread.EQ.0.OR.nthreads.EQ.1)THEN
       CALL my_concat_inplace(temp,il_temp(1)*nstate,parai%allgrp)
    END IF
    !$ IF (methread.EQ.1) THEN
    !$    CALL omp_set_max_active_levels(1)
    !$    CALL omp_set_num_threads(parai%ncpus)
#if defined(_INTEL_MKL) && !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$    CALL mkl_set_dynamic(1)
    !$    ierr = mkl_set_num_threads_local(0)
#endif
    !$ END IF
    !$omp end parallel

    CALL copy_back(fnl,temp,work_proc,ldf,nstate,parai%nproc)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_temp,temp,procedureN//'_temp',ierr)
#else
    DEALLOCATE(temp, stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate temp',&
         __LINE__,__FILE__)
    IF(redist)CALL cp_grp_redist_array_f(c0,ldc,nstate)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rottr_c0_fnl
  ! ==================================================================
  SUBROUTINE copy_back(fnl_p,temp,work_proc,ldf,nstate,n)
    INTEGER, INTENT(IN)                      :: n,nstate,ldf,work_proc(2,n)
    REAL(real_8), INTENT(IN) __CONTIGUOUS    :: temp(:,:,:)
    REAL(real_8), INTENT(OUT)                :: fnl_p(ldf,nstate)
    INTEGER                                  :: proc,i

    !$omp parallel do private(proc,i)
    DO i=1,nstate
       DO proc=1,n
          CALL copy(fnl_p(work_proc(1,proc):work_proc(2,proc),i),&
               temp(:,i,proc),work_proc(2,proc)-work_proc(1,proc)+1)
       END DO
    END DO
  END SUBROUTINE copy_back
  ! ==================================================================
  SUBROUTINE copy_in(temp,fnl_p,work_proc,ldt,nstate)
    INTEGER, INTENT(IN)                      :: nstate, ldt ,work_proc(2)
    REAL(real_8), INTENT(IN) __CONTIGUOUS    :: fnl_p(:,:)
    REAL(real_8), INTENT(OUT)                :: temp(ldt,nstate)
    INTEGER                                  :: i
    !$omp parallel do private(i)
    DO i=1,nstate
       CALL copy(temp(:,i),fnl_p(work_proc(1):work_proc(2),i),work_proc(2)-work_proc(1)+1)
    END DO
  END SUBROUTINE  copy_in
  ! ==================================================================
  SUBROUTINE rottr(a,c1,gam,transa,nstate,n,tlsd,na,nb,redist,use_cp)
    ! ==--------------------------------------------------------------==
    ! ==         C1 <= A*C1*GAM                                       ==
    ! ==   SPECIAL CASE FOR GAM UPPER TRIAGONAL                       ==
    ! ==--------------------------------------------------------------==
    REAL(real_8),INTENT(IN)                  :: a
    COMPLEX(real_8),INTENT(INOUT)            :: c1(:,:)
    REAL(real_8),INTENT(IN)                  :: gam(:,:)
    CHARACTER(len=*),INTENT(IN)              :: transa
    INTEGER,INTENT(IN)                       :: nstate, n, na, nb
    LOGICAL,INTENT(IN)                       :: tlsd
    LOGICAL,INTENT(IN),OPTIONAL              :: redist,use_cp

    INTEGER                                  :: isub, naa,ibeg_c0,ngw_local
    LOGICAL                                  :: rdst,cp_active
    CHARACTER(*), PARAMETER                  :: procedureN = 'rottr'

!(nstate,nstate)
! Variables
! ==--------------------------------------------------------------==
    IF (n.EQ.0 .OR. nstate.EQ.0) RETURN
    CALL tiset(procedureN,isub)
    !can we use cp_grp tricks?
    IF(PRESENT(use_cp))THEN
       cp_active=use_cp
    ELSE
       cp_active=.FALSE.
    END IF
    IF(cp_active)THEN
       CALL cp_grp_get_sizes(ngw_l=ngw_local,first_g=ibeg_c0)
       IF(PRESENT(redist))THEN
          rdst=redist
       ELSE
          rdst=.TRUE.
       END IF
    ELSE
       ngw_local=n
       ibeg_c0=1
       rdst=.FALSE.
    END IF
    IF (tlsd) THEN
       naa=na+1
       CALL cpmd_dtrmm('R','U',transa,'N',2*ngw_local,na,a,gam,nstate,&
            c1(ibeg_c0,1),2*n)
       CALL cpmd_dtrmm('R','U',transa,'N',2*ngw_local,nb,a,gam(naa,naa),&
            nstate,c1(ibeg_c0,naa),2*n)
    ELSE
       CALL cpmd_dtrmm('R','U',transa,'N',2*ngw_local,nstate,a,gam,&
            nstate,c1(ibeg_c0,1),2*n)
    ENDIF
    IF(rdst) CALL cp_grp_redist_array_f(c1,n,nstate)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rottr
  ! ==================================================================
  SUBROUTINE rotate(a,c1,b,c2,gam,nstate,n,tlsd,na,nb,symm,redist,use_cp)
    ! ==--------------------------------------------------------------==
    ! ==         C2 <= B*C2 + A*C1*GAM                                ==
    ! ==   ALSO FOR LSD                                               ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate,n,na,nb
    REAL(real_8),INTENT(IN)                  :: a,b
    COMPLEX(real_8),INTENT(IN)               :: c1(:,:)
    COMPLEX(real_8),INTENT(INOUT)            :: c2(:,:)
    REAL(real_8),INTENT(IN)                  :: gam(nstate,*)
    LOGICAL,INTENT(IN)                       :: tlsd
    LOGICAL,INTENT(IN),OPTIONAL              :: redist,use_cp,symm

    CHARACTER(*), PARAMETER                  :: procedureN = 'rotate'

    INTEGER                                  :: isub, naa, ibeg_c0, ngw_local, &
                                                first_g
    LOGICAL                                  :: rdst,cp_active,use_dsymm
! Variables
! ==--------------------------------------------------------------==
    IF (n.EQ.0 .OR. nstate.EQ.0) RETURN

    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )

    !can we use cp_grp tricks?
    IF(PRESENT(use_cp))THEN
       cp_active=use_cp
    ELSE
       cp_active=.FALSE.
    END IF
    IF(PRESENT(symm))THEN
       use_dsymm=symm
    ELSE
       use_dsymm=.FALSE.
    END IF
    IF(cp_active)THEN
       CALL cp_grp_get_sizes(ngw_l=ngw_local,first_g=ibeg_c0)
       IF(PRESENT(redist))THEN
          rdst=redist
       ELSE
          rdst=.TRUE.
       END IF
    ELSE
       ngw_local=n/2
       ibeg_c0=1
       rdst=.FALSE.
    END IF
    IF(use_dsymm)THEN
       IF (n>0) THEN
          IF (tlsd) THEN
             naa=na+1
             CALL cpmd_dsymm('R','U',2*ngw_local,na,a,gam(1,1),nstate,&
                  c1(ibeg_c0,1),n,b,c2(ibeg_c0,1),n)
             CALL cpmd_dsymm('R','U',2*ngw_local,nb,a,gam(naa,naa),nstate,&
                  c1(ibeg_c0,naa),n,b,c2(ibeg_c0,naa),n)
          ELSE
             CALL cpmd_dsymm('R','U',2*ngw_local,nstate,a,gam,nstate,&
                  c1(ibeg_c0,1),n,b,c2(ibeg_c0,1),n)
          ENDIF
       ENDIF
    ELSE
       IF (n>0) THEN
          IF (tlsd) THEN
             naa=na+1
             CALL cpmd_dgemm('N','N',2*ngw_local,na,na,a,c1(ibeg_c0,1),n,&
                  gam(1,1),nstate,b,c2(ibeg_c0,1),n)
             CALL cpmd_dgemm('N','N',2*ngw_local,nb,nb,a,c1(ibeg_c0,naa),n,&
                  gam(naa,naa),nstate,b,c2(ibeg_c0,naa),n)
          ELSE
             CALL cpmd_dgemm('N','N',2*ngw_local,nstate,nstate,a,c1(ibeg_c0,1),&
                  n,gam,nstate,b,c2(ibeg_c0,1),n)
          ENDIF
       ENDIF
    END IF
    IF(rdst)CALL cp_grp_redist_array_f(c2,n/2,nstate)

    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE rotate
  ! ==================================================================
END MODULE rotate_utils

! ==================================================================
SUBROUTINE rotate_da(a,c1,b,c2,gam,ld,n,nstate,ndd1,ndd2,nddx,mepos,pgroup,nproc,grp,tlsd,na,nb)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_bcast
  USE cuda_types,                      ONLY: cuda_memory_t
  USE cuda_utils,                      ONLY: cuda_alloc_bytes, cuda_dealloc, &
       cuda_memcpy_host_to_device, cuda_memcpy_device_to_host, cuda_d_points_to, cuda_mem_zero_bytes
  USE cublas_types,                    ONLY: cublas_handle_t
  USE cublas_utils,                    ONLY: cublas_create, cublas_destroy, cublas_dgemm, cublas_dger
  USE sizeof_kinds,                    ONLY: sizeof_real_8
  USE cp_cuwfn_types, ONLY: cp_cuwfn_device_get_ptrs, cp_cuwfn, cp_cuwfn_get
  USE nvtx_utils
#ifdef __PARALLEL
  USE mpi_f08
#endif
  IMPLICIT NONE
  REAL(real_8)                               :: a, b
  INTEGER                                    :: ld, n, nstate
  REAL(real_8)                               :: c1(ld,nstate), c2(ld,nstate), &
                                                gam(nstate,nstate)
#ifdef __PARALLEL
  INTEGER                                    :: ndd1(0:*), ndd2(0:*), nddx, &
                                                mepos, pgroup(0:*), nproc
  type(MPI_COMM)                             :: grp
#else
  INTEGER                                    :: ndd1(0:*), ndd2(0:*), nddx, &
                                                mepos, pgroup(0:*), nproc, grp
#endif
  LOGICAL                                    :: tlsd
  INTEGER                                    :: na, nb

  CHARACTER(*), PARAMETER                    :: procedureN = 'rotate_da'
  LOGICAL, PARAMETER                         :: use_gpu = .FALSE.

  INTEGER                                    :: i, i1, ierr, ii, ip, isub, j, &
                                                n1, nmax, nmin
  INTEGER(int_8)                             :: n_bytes
  REAL(real_8), ALLOCATABLE                  :: aux(:)
  REAL(real_8), ALLOCATABLE, DIMENSION(:, :) :: c2_tmp
  TYPE(cublas_handle_t)                      :: blas_handle
  TYPE(cuda_memory_t)                        :: c1_d, c2_d, gam_d

! ==--------------------------------------------------------------==

  CALL tiset(procedureN,isub)
  __NVTX_TIMER_START ( procedureN )

  IF (tlsd) THEN
     IF (ndd1(mepos).LT.na) THEN
        nmax=MIN(na,ndd2(mepos))
        DO i=ndd1(mepos),nmax
           ii=i-ndd1(mepos)+1
           DO j=na+1,nstate
              gam(j,ii)=0._real_8
           ENDDO
        ENDDO
     ENDIF
     IF (ndd2(mepos).GT.na) THEN
        nmin=MAX(na+1,ndd1(mepos))
        DO i=nmin,ndd2(mepos)
           ii=i-ndd1(mepos)+1
           DO j=1,na
              gam(j,ii)=0._real_8
           ENDDO
        ENDDO
     ENDIF
  ENDIF
  CALL dscal(nstate*n,b,c2,1)
  ALLOCATE(aux(nstate*nddx),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  CALL dcopy(nstate*nddx,gam,1,aux,1)

  IF( use_gpu ) THEN
     WRITE(*,*) 'rotate_utils.mod.F90: use_gpu',use_gpu
     CALL cublas_create ( blas_handle, 0 )

     n_bytes = SIZE( c1, KIND=int_8 ) * INT( sizeof_real_8, int_8 )
     CALL cuda_alloc_bytes ( c1_d, n_bytes, 0 )
     CALL cuda_memcpy_host_to_device ( c1, c1_d )

     n_bytes = SIZE( c1, KIND=int_8 ) * INT( sizeof_real_8, int_8 )
     CALL cuda_alloc_bytes ( c2_d, n_bytes, 0 )
     CALL cuda_mem_zero_bytes ( c2_d, n_bytes )

     n_bytes = SIZE( gam, KIND=int_8 ) * INT( sizeof_real_8, int_8 )
     CALL cuda_alloc_bytes ( gam_d, n_bytes, 0 )
  ENDIF

  DO ip=0,nproc-1
     i1=ndd1(ip)
     n1=ndd2(ip)-ndd1(ip)+1
     IF (n1.GT.0) THEN
        CALL dcopy(nstate*nddx,aux,1,gam,1)
        CALL mp_bcast(gam,nstate*nddx,pgroup(ip+1),grp)
        IF (ld>0) THEN
           IF( use_gpu ) THEN
              CALL cuda_memcpy_host_to_device ( gam, gam_d )
              CALL cublas_dgemm ( blas_handle, 'N', 'N', n, n1, nstate, &
                   & a, c1_d, ld, &
                   & gam_d, nstate, &
                   & 1.0_real_8, cuda_d_points_to(c2_d,(i1-1)*ld+1), ld )
           ELSE
              CALL cpmd_dgemm('N','N',n,n1,nstate,a,c1(1,1),ld,&
                   gam(1,1),nstate,1._real_8,c2(1,i1),ld)
           ENDIF
        ENDIF
     ENDIF
  ENDDO

  IF( use_gpu ) THEN
     ALLOCATE(c2_tmp(SIZE(c2,1),SIZE(c2,2)),STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
          __LINE__,__FILE__)
     CALL cuda_memcpy_device_to_host ( c2_d, c2_tmp )
     CALL daxpy ( SIZE(c2), 1.0_real_8, c2_tmp, 1, c2, 1 )
     DEALLOCATE(c2_tmp,STAT=ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
          __LINE__,__FILE__)
     CALL cuda_dealloc ( gam_d )
     CALL cuda_dealloc ( c2_d )
     CALL cuda_dealloc ( c1_d )
     CALL cublas_destroy ( blas_handle )
  ENDIF

  CALL dcopy(nstate*nddx,aux,1,gam,1)
  DEALLOCATE(aux,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)

  __NVTX_TIMER_STOP
  CALL tihalt(procedureN,isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE rotate_da
! ==================================================================
SUBROUTINE rotate_c(a,c1,b,c2,gam,nstate)
  ! ==--------------------------------------------------------------==
  ! ==         C2 <= B*C2 + A*C1*GAM                                ==
  ! ==   ALSO FOR LSD                                               ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE system , ONLY:cntl,nkpt
  USE parac, ONLY : paral,parai
  USE spin , ONLY:spin_mod
  IMPLICIT NONE
  COMPLEX(real_8)                            :: a, c1(nkpt%ngwk,*), b, &
                                                c2(nkpt%ngwk,*)
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: gam(nstate,*)

  INTEGER                                    :: isub, naa

  IF (nkpt%ngwk.EQ.0 .OR. nstate.EQ.0) RETURN
  CALL tiset('  ROTATE_C',isub)
  IF (cntl%tlsd) THEN
     naa=spin_mod%nsup+1
     IF (nkpt%ngwk>0) THEN
        CALL cpmd_zgemm('N','N',nkpt%ngwk,spin_mod%nsup,spin_mod%nsup,a,c1(1,1),nkpt%ngwk,gam(1,1),&
             nstate,b,c2(1,1),nkpt%ngwk)
        CALL cpmd_zgemm('N','N',nkpt%ngwk,spin_mod%nsdown,spin_mod%nsdown,a,c1(1,naa),nkpt%ngwk,&
             gam(naa,naa),nstate,b,c2(1,naa),nkpt%ngwk)
     ENDIF
  ELSE
     IF (nkpt%ngwk>0) THEN
        CALL cpmd_zgemm('N','N',nkpt%ngwk,nstate,nstate,a,c1(1,1),nkpt%ngwk,&
             gam(1,1),nstate,b,c2(1,1),nkpt%ngwk)
     ENDIF
  ENDIF
  CALL tihalt('  ROTATE_C',isub)
  ! ==--------------------------------------------------------------==
END SUBROUTINE rotate_c
