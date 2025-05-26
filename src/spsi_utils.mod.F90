#include "cpmd_global.h"

MODULE spsi_utils
  USE cppt,                            ONLY: twnl_nghtol
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms,&
                                             cp_grp_get_sizes
  USE beta_utils,                      ONLY: build_beta
  USE cvan,                            ONLY: qq
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8,&
                                             int_8,&
                                             int_4
  USE nlps,                            ONLY: nghtol,&
                                             nlps_com,&
                                             imagp
  !$ USE omp_lib,                      ONLY: omp_set_nested,&
  !$                                         omp_get_thread_num,&
  !$                                         omp_set_num_threads
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: eigr,&
                                             fnl,&
                                             il_fnl_packed
  USE system,                          ONLY: cntl,&
                                             ncpw,&
                                             maxsys
  USE timer,                           ONLY: tihalt,&
                                             tiset
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif
#ifdef _INTEL_MKL
  use mkl_service
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: spsi

CONTAINS

  ! ==================================================================
  SUBROUTINE spsi(nstate,sc0,fnl_p,redist)
    ! == Rewritten: Tobias Kloeffel, FAU Erlangen-Nuernberg May 2019  ==
    ! == Introduce overlapping communication computation algorithm    ==
    ! == enable via USE_OVERLAPPING_COMM_COMP in the &CPMD section    ==
    ! == relies on fnl_packed calculated by rottr_c0_fnl and/or       ==
    ! == rnlsm1, rotate_c0_fnl                                        ==
    ! == cp_grp distribution of fnl is taken care of                  ==
    ! == full performance only with saved arrays or scratch_library   ==
    ! == technically speaking this is just a copy of nlforce, with    ==
    ! == the exception of a different build_dai routine               ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate
    COMPLEX(real_8),INTENT(INOUT)            :: sc0(ncpw%ngw,nstate)
    REAL(real_8),INTENT(IN)                  :: fnl_p(il_fnl_packed(1),nstate)
    LOGICAL,INTENT(IN)                       :: redist

    character(*), PARAMETER                  :: proceduren = 'spsi'
    INTEGER                                  :: i, offset_fnl, offset_dai, isa0, ibeg, ngw_local,&
                                                is, ia_sum, fnl_start, ia_fnl, isub, ierr, &
                                                ld_grp(0:parai%cp_nogrp-1), grp,&
                                                 nthreads, nested_threads, methread
    INTEGER(int_8)                           :: il_dai(3), il_eiscr(2)
    INTEGER, ALLOCATABLE                     :: na_grp(:,:,:), na_fnl(:,:), na(:,:)

#ifdef _USE_SCRATCHLIBRARY
    COMPLEX(real_8),POINTER __CONTIGUOUS     :: eiscr(:,:)
    REAL(real_8),POINTER __CONTIGUOUS        :: dai(:,:,:)
#else
    REAL(real_8),ALLOCATABLE                 :: dai(:,:,:)
    COMPLEX(real_8),ALLOCATABLE              :: eiscr(:,:)
#endif

    CALL tiset(procedureN,isub)
    IF (imagp.EQ.2) call stopgm(procedureN,'k-point not implemented',&
         __LINE__,__FILE__)
    IF (cntl%tfdist) call stopgm(procedureN,'fnl dist. not implemented',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ALLOCATE(na(2,ions1%nsp),na_fnl(2,ions1%nsp),na_grp(2,ions1%nsp,0:parai%cp_nogrp-1)&
         , stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'allocation problem',&
         __LINE__,__FILE__)

    CALL cp_grp_split_atoms(na_grp)
    na_fnl(:,:)=na_grp(:,:,parai%cp_inter_me)
    !if cp_groups are used we exchange the dia arrays between cp_groups
    !this way we can do the exchange phase during the local dgemm calculation
    !if we need the full sc0 we apply all dia arrays to the full set of
    !ngws, else we just update the cp_grp local ngws
    IF(redist)THEN
       ngw_local=ncpw%ngw
       ibeg=1
    ELSE
       CALL cp_grp_get_sizes(ngw_l=ngw_local,first_g=ibeg)
    END IF
    IF (pslo_com%tivan) THEN
       ! ==--------------------------------------------------------==
       ! ==  vanderbilt pp                                         ==
       ! ==--------------------------------------------------------==
       ld_grp=0
       DO grp=0,parai%cp_nogrp-1
          DO is=1,ions1%nsp
             IF (pslo_com%tvan(is) ) THEN
                ld_grp(grp)=ld_grp(grp)+(na_grp(2,is,grp)-na_grp(1,is,grp)+1)&
                     *nlps_com%ngh(is)
             ELSE
                !filter out non uspp atoms
                na_grp(1,is,grp)=0
                na_grp(2,is,grp)=-1
             END IF
          END DO
       END DO
       na(:,:)=na_grp(:,:,parai%cp_inter_me)
       il_dai(1)=MAXVAL(ld_grp)
       il_dai(2)=nstate
       il_dai(3)=parai%cp_nogrp
       il_eiscr(1)=ngw_local
       il_eiscr(2)=MAXVAL(ld_grp)
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_dai,dai,procedureN//'_dai',ierr)
#else
       ALLOCATE(dai(il_dai(1),il_dai(2),il_dai(3)), stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate dai',&
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_eiscr,eiscr,procedureN//'_eiscr',ierr)
#else
       ALLOCATE(eiscr(il_eiscr(1),il_eiscr(2)), stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate eiscr',&
            __LINE__,__FILE__)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute  &
#else
       !$omp parallel do &
#endif
       !$omp& private (i,offset_fnl,offset_dai,isa0,is,ia_fnl,ia_sum,fnl_start)
       DO i=1,nstate
          !offset for packed fnl/fnlgam
          offset_fnl=1
          !offset for dai
          isa0=0
          offset_dai=1
          !fill local part of dai
          DO is=1,ions1%nsp
             ia_fnl=na_fnl(2,is)-na_fnl(1,is)+1
             ia_sum=na(2,is)-na(1,is)+1
             fnl_start=na(1,is)-na_fnl(1,is)
             IF(ia_sum.GT.0)THEN
                      !DIR$ forceinline
                CALL build_dai(dai(offset_dai:,i,&
                     parai%cp_inter_me+1),&
                     fnl_p(offset_fnl,i),&
                     qq(:,:,is),nlps_com%ngh(is),ia_sum,ia_fnl,fnl_start,&
                     maxsys%nhxs)
                offset_dai=offset_dai+nlps_com%ngh(is)*ia_sum
             END IF
             offset_fnl=offset_fnl+nlps_com%ngh(is)*ia_fnl
             isa0=isa0+ions0%na(is)
          END DO
       END DO
!#if defined(_HAS_OMP_TARGET_OFFLOAD)
!       !$omp target update from(dai)
!#endif
       IF(cntl%overlapp_comm_comp)THEN
          nthreads=MIN(2,parai%ncpus)
          nested_threads=(MAX(parai%ncpus-1,1))
#if !defined(_INTEL_MKL)
          CALL stopgm(procedureN, 'Overlapping communication and computation: Behavior of &
               BLAS routine inside parallel region not checked',&
               __LINE__,__FILE__)
#endif
       ELSE
          nthreads=1
          nested_threads=parai%ncpus
       END IF
       methread=0
       IF(parai%cp_nogrp.EQ.1)THEN
          nthreads=1
          nested_threads=parai%ncpus
       END IF
       methread=0
       !$omp parallel if(nthreads.EQ.2) num_threads(nthreads) &
       !$omp private(methread,grp) proc_bind(close)
       !$ methread = omp_get_thread_num()
       IF(methread.EQ.0.AND.parai%cp_nogrp.GT.1)THEN
          !get data from other cp_grp other threads build local beta and perform dgemms
          CALL my_concat_inplace(dai,il_dai(1)*nstate,parai%cp_inter_grp)
       END IF
       IF(methread.EQ.1.OR.nthreads.EQ.1)THEN
          !$ IF(methread.EQ.1)THEN
          !$    CALL omp_set_max_active_levels(2)
          !$    CALL omp_set_num_threads(nested_threads)
#if defined(_INTEL_MKL) && !defined(_HAS_OMP_TARGET_OFFLOAD)
          !$    CALL mkl_set_dynamic(0)
          !$    ierr = mkl_set_num_threads_local(nested_threads)
#endif
          !$ END IF

          grp=parai%cp_inter_me
          CALL build_beta(na_grp(:,:,grp),eigr,eiscr,ncpw%ngw,ibeg,&
               ngw_local,twnl_nghtol=twnl_nghtol(:,:,:,1))
          IF(ld_grp(grp).GT.0)THEN
             CALL CPMD_DGEMM('N','N',2*ngw_local,nstate,ld_grp(grp)&
                  ,1._real_8,eiscr(1,1),2*ngw_local&
                  ,dai(1,1,grp+1),INT(il_dai(1),kind=int_4),1.0_real_8,sc0(ibeg,1),2*ncpw%ngw)
          END IF
          
          !$ IF (methread.EQ.1) THEN
          !$    CALL omp_set_max_active_levels(1)
          !$    CALL omp_set_num_threads(parai%ncpus)
#if defined(_INTEL_MKL) && !defined(_HAS_OMP_TARGET_OFFLOAD)
          !$    CALL mkl_set_dynamic(1)
          !$    ierr = mkl_set_num_threads_local(0)
#endif
          !$ END IF
       END IF
       !$omp end parallel
       !now we apply the dai arrays of all other groups
       IF(parai%cp_nogrp.GT.1)THEN
          DO grp=0,parai%cp_nogrp-1
             IF(grp.EQ.parai%cp_inter_me)CYCLE
             CALL build_beta(na_grp(:,:,grp),eigr,eiscr,ncpw%ngw,ibeg,&
                  ngw_local,twnl_nghtol=twnl_nghtol(:,:,:,1))
             IF(ld_grp(grp).GT.0)THEN
                CALL CPMD_DGEMM('N','N',2*ngw_local,nstate,ld_grp(grp)&
                     ,1._real_8,eiscr(1,1),2*ngw_local&
                     ,dai(1,1,grp+1),INT(il_dai(1),kind=int_4),1.0_real_8,sc0(ibeg,1),2*ncpw%ngw)
             END IF
          END DO
       END IF
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_eiscr,eiscr,procedureN//'_eiscr',ierr)
#else
    DEALLOCATE(eiscr, stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate eiscr',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_dai,dai,procedureN//'_dai',ierr)
#else
    DEALLOCATE(dai, stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate dai',&
         __LINE__,__FILE__)
    DEALLOCATE(na_grp, na, na_fnl, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'deallocation problem)',&
         __LINE__,__FILE__)

    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE spsi
  ! ==================================================================
  !DIR$ ATTRIBUTES FORCEINLINE::build_dai
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  SUBROUTINE build_dai(mat,fnl_p,qq_,ngh,ia_sum,ia_fnl,fnl_start,maxngh)
#else
  PURE SUBROUTINE build_dai(mat,fnl_p,qq_,ngh,ia_sum,ia_fnl,fnl_start,maxngh)
#endif
    INTEGER,INTENT(IN)                       :: ngh,ia_sum,ia_fnl,fnl_start,&
                                                maxngh
    REAL(real_8),INTENT(IN)                  :: fnl_p(ia_fnl,ngh,*),&
                                                qq_(maxngh,*)
    REAL(real_8),INTENT(OUT)                 :: mat(ia_sum,ngh,*)
    INTEGER                                  :: iv,ia,jv,isa
    real(real_8) :: fac1
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp parallel do private(iv,ia,jv,isa,fac1)
#endif
    DO iv=1,ngh
       DO ia=1,ia_sum
          mat(ia,iv,1)=0.0_real_8
       END DO
       DO jv=1,ngh
          isa=fnl_start
          fac1=qq_(jv,iv)
          IF (ABS(fac1).GT.1.e-5_real_8) THEN
             DO ia=1,ia_sum
                isa=isa+1
                mat(ia,iv,1)=mat(ia,iv,1)&
                     +qq_(jv,iv)*fnl_p(isa,jv,1)
             END DO
          END IF
       END DO
    END DO
  END SUBROUTINE build_dai
  ! ==================================================================

END MODULE spsi_utils
