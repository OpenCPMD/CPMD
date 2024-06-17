#include "cpmd_global.h"

MODULE crotwf_utils
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes
  USE error_handling,                  ONLY: stopgm
  USE gpu
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE part_1d,                         ONLY: part_1d_get_blk_bounds
  USE pslo,                            ONLY: pslo_com
  USE nort,                            ONLY: nort_ovlap,&
                                             nort_com
  USE rotate_utils,                    ONLY: rotate
  USE spin,                            ONLY: spin_mod
  USE summat_utils,                    ONLY: summat
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE utils,                           ONLY: dsyevx_driver,&
                                             dsyevd_driver,&
                                             elpa_driver
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: crotwf
CONTAINS
  ! ==================================================================
  SUBROUTINE crotwf(c0,cm,c2,sc0,nstate,gam,use_cp_grps)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate
    COMPLEX(real_8),INTENT(INOUT)            :: c2(ncpw%ngw,nstate), &
                                                cm(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    COMPLEX(real_8),INTENT(OUT)              :: sc0(ncpw%ngw,nstate)
    REAL(real_8),INTENT(OUT)                 :: gam(nstate,nstate)
    LOGICAL,INTENT(IN)                       :: use_cp_grps

    CHARACTER(*), PARAMETER                  :: procedureN = 'crotwf'

    INTEGER                                  :: i, ierr, j, isub, &
                                                ibeg, iend, ig, isub1
    INTEGER(int_8)                           :: il_eigval(1), il_temp(2)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: eigval(:),temp(:,:),temp1(:,:)
#else
    REAL(real_8), ALLOCATABLE                :: eigval(:),temp(:,:),temp1(:,:)
#endif
    LOGICAL                                  :: nopara
    CALL tiset(procedureN,isub)

    IF(cntl%tlsd)THEN
       il_eigval(1)=max(spin_mod%nsup,spin_mod%nsdown)
    ELSE
       il_eigval(1)=nstate
    END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    update_first_to_gpu  =.FALSE.
    update_second_to_gpu =.FALSE.
    update_third_to_gpu  =.FALSE.
    update_result_to_host=.FALSE.
    comm_buffers_on_host =.FALSE.
#endif
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_eigval,eigval,procedureN//'_eigval',ierr)
#else
    ALLOCATE(eigval(il_eigval(1)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    nopara=(nort_com%scond.LT.1.e-9_real_8.OR.parai%cp_nproc.LT.17.OR.nstate.LT.1000)

    CALL tiset(procedureN//'solve',isub1)
    IF (.NOT.cntl%tlsd) THEN
       CALL solve_eigenvector(nstate,eigval,nort_ovlap,gam,nopara)
    ELSE
       !for the moment we copy each spins out in a temporary buffer...
       !spin up
       il_temp(1)=spin_mod%nsup
       il_temp(2)=spin_mod%nsup
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_temp,temp,procedureN//'_temp',ierr)
#else
       ALLOCATE(temp(il_temp(1),il_temp(2)),STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_temp,temp1,procedureN//'_temp1',ierr)
#else
       ALLOCATE(temp1(il_temp(1),il_temp(2)),STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       !$omp parallel do private(i,j)
       DO i=1,spin_mod%nsup
          DO j=1,i
             temp(j,i)=nort_ovlap(j,i)
          END DO
       END DO
       CALL solve_eigenvector(spin_mod%nsup,eigval,temp,temp1,nopara)
       !$omp parallel do private(i,j)
       DO i=1,spin_mod%nsup
          DO j=1,spin_mod%nsup
             nort_ovlap(j,i)=temp(j,i)
          END DO
       END DO
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_temp,temp1,procedureN//'_temp1',ierr)
#else
       DEALLOCATE(temp1,STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_temp,temp,procedureN//'_temp',ierr)
#else
       DEALLOCATE(temp,STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       !spin down
       il_temp(1)=spin_mod%nsdown
       il_temp(2)=spin_mod%nsdown
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_temp,temp,procedureN//'_temp',ierr)
#else
       ALLOCATE(temp(il_temp(1),il_temp(2)),STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_temp,temp1,procedureN//'_temp1',ierr)
#else
       ALLOCATE(temp1(il_temp(1),il_temp(2)),STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       !$omp parallel do private(i,j)
       DO i=spin_mod%nsup+1,nstate
          DO j=spin_mod%nsup+1,i
             temp(j-spin_mod%nsup,i-spin_mod%nsup)=nort_ovlap(j,i)
          END DO
       END DO
       CALL solve_eigenvector(spin_mod%nsdown,eigval,temp,temp1,nopara)
       !$omp parallel do private(i,j)
       DO i=spin_mod%nsup+1,nstate
          DO j=spin_mod%nsup+1,nstate
             nort_ovlap(j,i)=temp(j-spin_mod%nsup,i-spin_mod%nsup)
          END DO
       END DO
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_temp,temp1,procedureN//'_temp1',ierr)
#else
       DEALLOCATE(temp1,STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_temp,temp,procedureN//'_temp',ierr)
#else
       DEALLOCATE(temp,STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    CALL tihalt(procedureN//'solve',isub1)

    IF(use_cp_grps)THEN
       CALL cp_grp_get_sizes(first_g=ibeg,last_g=iend)
    ELSE
       ibeg=1
       iend=ncpw%ngw
    END IF
    CALL rotate(1.0_real_8,c0,0.0_real_8,sc0,nort_ovlap,nstate,2*ncpw%ngw,cntl%tlsd,&
         spin_mod%nsup,spin_mod%nsdown,redist=.NOT.use_cp_grps,use_cp=use_cp_grps)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute parallel do collapse(2)
#else
    !$omp parallel do private(ig,i)
#endif
    DO i=1,nstate
       DO ig=ibeg,iend
          c0(ig,i)=sc0(ig,i)
       END DO
    END DO
    CALL rotate(1.0_real_8,cm,0.0_real_8,sc0,nort_ovlap,nstate,2*ncpw%ngw,cntl%tlsd,&
         spin_mod%nsup,spin_mod%nsdown,redist=.NOT.use_cp_grps,use_cp=use_cp_grps)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute parallel do collapse(2)
#else
    !$omp parallel do private(ig,i)
#endif
    DO i=1,nstate
       DO ig=ibeg,iend
          cm(ig,i)=sc0(ig,i)
       END DO
    END DO
    CALL rotate(1.0_real_8,c2,0.0_real_8,sc0,nort_ovlap,nstate,2*ncpw%ngw,cntl%tlsd,&
         spin_mod%nsup,spin_mod%nsdown,redist=.NOT.use_cp_grps,use_cp=use_cp_grps)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute parallel do collapse(2)
#else
    !$omp parallel do private(ig,i)
#endif
    DO i=1,nstate
       DO ig=ibeg,iend
          c2(ig,i)=sc0(ig,i)
       END DO
    END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    update_first_to_gpu  =.TRUE.
    update_second_to_gpu =.TRUE.
    update_third_to_gpu  =.TRUE.
    update_result_to_host=.TRUE.
    comm_buffers_on_host =.TRUE.
#endif
    ! ==--------------------------------------------------------------==
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_eigval,eigval,procedureN//'_eigval',ierr)
#else
    DEALLOCATE(eigval,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==

    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE crotwf
  ! ==================================================================
  SUBROUTINE solve_eigenvector(nstate,eigval,ovlap,gam,nopara)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate
    REAL(real_8),INTENT(OUT) __CONTIGUOUS    :: eigval(:)
    REAL(real_8),INTENT(OUT)                 :: gam(nstate,nstate)
    REAL(real_8),INTENT(INOUT)               :: ovlap(nstate,nstate)
    LOGICAL,INTENT(IN)                       :: nopara

    INTEGER                                  :: i,iopt,first,last,&
                                                chunks(2,0:parai%cp_nproc-1),&
                                                recvcnt(0:parai%cp_nproc-1),&
                                                displ(0:parai%cp_nproc-1),isub

    iopt=21
    IF(cntl%use_elpa)THEN
       CALL elpa_driver(ovlap,eigval,nstate)
       CALL mp_bcast(ovlap,nstate**2,parai%io_source,parai%cp_grp)
    ELSE
       IF(nopara) THEN
          !parallelization using multiple dsyevx/r does not work =>
          !fall back to dsyevd on root
          IF(paral%io_parent) CALL dsyevd_driver(iopt,ovlap,eigval,nstate)
          CALL tiset('bcast',isub)
          CALL mp_bcast(ovlap,nstate**2,parai%io_source,parai%cp_grp)
          CALL tihalt('bcast',isub)
       ELSE
          !we can distribute the eigenvalue problem by using multiple
          !dsyevx/r instances
          !Not very efficient parallelization though
          !For numerical stability we should consider to serialize and broadcast the
          !tridiagonilization step
          recvcnt=-1
          displ=0
          DO i = 0,parai%cp_nproc-1
             CALL part_1d_get_blk_bounds(nstate,i,parai%cp_nproc,chunks(1,i),chunks(2,i))
             recvcnt(i)=(chunks(2,i)-chunks(1,i)+1)*nstate
             IF (i.GT.0) displ(i)=displ(i-1)+recvcnt(i-1)
          END DO
          first=chunks(1,parai%cp_me)
          last=chunks(2,parai%cp_me)
          CALL dsyevx_driver(iopt,ovlap,gam,eigval,nstate,first,last,-1.0_real_8)
          CALL my_concatv(gam,ovlap,(last-first+1)*nstate,recvcnt,displ,parai%cp_grp)
       END IF
    END IF
  END SUBROUTINE solve_eigenvector
  ! ==================================================================

END MODULE crotwf_utils
