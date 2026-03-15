#include "cpmd_global.h"

MODULE csmat_utils
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms
  USE cvan,                            ONLY: qq
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE distribution_utils,              ONLY: dist_atoms
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE gpu
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: nlps_com
  USE nort,                            ONLY: nort_com,&
                                             nort_ovlap
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: fnl_packed
  USE spin,                            ONLY: spin_mod
  USE summat_utils,                    ONLY: summat
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: csmat

CONTAINS

  ! ==================================================================
  ! == FOR TKPNT=.TRUE. FNL IS COMPLEX -- CSMAT_C WILL BE WRITTEN   ==
  ! ==================================================================
  SUBROUTINE csmat(a,c0,nstate,ikind,full,store_nonort,only_parent)
    ! ==--------------------------------------------------------------==
    ! ==         COMPUTES THE OVERLAP MATRIX A = < C0 |S| C0 >        ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: ikind, nstate
    COMPLEX(real_8),INTENT(IN)               :: c0(ncpw%ngw,nstate)
    REAL(real_8),INTENT(OUT)                 :: a(nstate,*)
    LOGICAL,INTENT(IN)                       :: full,only_parent,store_nonort
    INTEGER                                  :: i, j, ia, is, ispin, nspin,&
                                                ia_sum, tot_work, ns(2), nmin(2), &
                                                off_i, off_mat, off_fnl, ia_fnl, &
                                                start_fnl, end_fnl, fnl_start, start_mat, &
                                                end_mat, isub, ierr, n
    INTEGER(int_8)                           :: il_fnlat(2), il_fnlatj(2)
    REAL(real_8)                             :: fractions(parai%nproc),selem, temp
    CHARACTER(*), PARAMETER                  :: procedureN = 'csmat'
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8),POINTER __CONTIGUOUS        :: fnlat(:,:), fnlatj(:,:)
#else
    REAL(real_8),ALLOCATABLE                 :: fnlat(:,:), fnlatj(:,:)
#endif
    INTEGER,ALLOCATABLE,SAVE                 :: na_buff(:,:,:)
    INTEGER,ALLOCATABLE                      :: na(:,:), na_fnl(:,:), na_grp(:,:,:)
    CALL tiset(procedureN,isub)
    IF (cntl%tfdist) CALL stopgm('CSMAT','TFDIST NOT IMPLEMENTED',&
         __LINE__,__FILE__)

    !we only need the upper part here
    CALL ovlap(nstate,a,c0,c0,redist=.FALSE.,full=.FALSE.)
    IF(store_nonort) CALL store_ovlap(a,nstate)

    IF (pslo_com%tivan) THEN
       ALLOCATE(na(2,ions1%nsp),na_fnl(2,ions1%nsp),na_grp(2,ions1%nsp,0:parai%cp_nogrp-1)&
            , stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'allocation problem',&
            __LINE__,__FILE__)
       !get cp_grp atom distribution
       CALL cp_grp_split_atoms(na_grp)
       na_fnl(:,:)=na_grp(:,:,parai%cp_inter_me)
       IF (.NOT.ALLOCATED(na_buff))THEN
          !distribute atoms between procs
          ALLOCATE(na_buff(2,ions1%nsp,0:parai%nproc-1), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate na_buff',&
               __LINE__,__FILE__)
          fractions=1.0_real_8/REAL(parai%nproc,KIND=real_8)
          CALL dist_atoms(parai%nproc,fractions,na_fnl,na_buff,only_uspp=.TRUE.)
       END IF
       na(:,:)=na_buff(:,:,parai%mepos)

       ! VANDERBILT PP
       ns(1)=nstate
       ns(2)=0
       nmin(1)=1
       nmin(2)=0
       nspin=1
       IF (cntl%tlsd) THEN
          ns(1)=spin_mod%nsup
          ns(2)=spin_mod%nsdown
          nmin(1)=1
          nmin(2)=spin_mod%nsup+1
          nspin=2
       ENDIF
       !end preparation

       !calculate local part of total work
       tot_work=0
       DO is=1,ions1%nsp
          tot_work=tot_work+(na(2,is)-na(1,is)+1)*nlps_com%ngh(is)
       END DO

       il_fnlat(1)=tot_work
       il_fnlat(2)=MAXVAL(ns)
       il_fnlatj(1)=tot_work
       il_fnlatj(2)=MAXVAL(ns)
       IF(tot_work.GT.0)THEN
#ifdef _USE_SCRATCHLIBRARY
          CALL request_scratch(il_fnlat,fnlat,procedureN//'_fnlat',ierr)
#else
          ALLOCATE(fnlat(il_fnlat(1),il_fnlat(2)), stat=ierr)
#endif
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlat',&
               __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
          CALL request_scratch(il_fnlatj,fnlatj,procedureN//'_fnlatj',ierr)
#else
          ALLOCATE(fnlatj(il_fnlatj(1),il_fnlatj(2)), stat=ierr)
#endif
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlatj',&
               __LINE__,__FILE__)

          DO ispin=1,nspin
             off_i=0
             IF(ispin.EQ.2) off_i=spin_mod%nsup
             n=ns(ispin)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp target teams distribute &
#else
             !$omp parallel do &
#endif
             !$omp& private(i,off_mat,off_fnl,is,ia_fnl,ia_sum,fnl_start)
             DO i=1,n
                off_mat=1
                off_fnl=1
                DO is=1,ions1%nsp
                   ia_fnl=na_fnl(2,is)-na_fnl(1,is)+1
                   ia_sum=na(2,is)-na(1,is)+1
                   !nc alreday filtered out in dist_atoms
                   IF(ia_sum.GT.0)THEN
                      !starting index fnl_packed
                      fnl_start=na(1,is)-na_fnl(1,is)
                      !DIR$ forceinline
                      CALL prepare_matrix(fnl_packed(off_fnl:,i+off_i),&
                           fnlat(off_mat:,i),&
                           fnlatj(off_mat:,i),qq(:,:,is),nlps_com%ngh(is),&
                           ia_sum,ia_fnl,fnl_start,maxsys%nhxs)
                      off_mat=off_mat+ia_sum*nlps_com%ngh(is)
                   END IF
                   off_fnl=off_fnl+ia_fnl*nlps_com%ngh(is)
                END DO
             END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp target update from(fnlat) if(update_first_to_gpu)
             !$omp target update from(fnlatj) if(update_second_to_gpu)
#endif
#ifdef _HAS_DGEMMT
             CALL cpmd_dgemmt('U','T','N',n,tot_work,1.0_real_8,&
                  fnlat(1,1),tot_work,fnlatj(1,1),tot_work,1.0_real_8,&
                  a(nmin(ispin),nmin(ispin)),nstate)
#else
             CALL cpmd_dgemm('T','N',n,n,tot_work,1.0_real_8,&
                  fnlat(1,1),tot_work,fnlatj(1,1),tot_work,1.0_real_8,&
                  a(nmin(ispin),nmin(ispin)),nstate)
#endif
          END DO
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_fnlatj,fnlatj,procedureN//'_fnlatj',ierr)
#else
          DEALLOCATE(fnlatj, stat=ierr)
#endif
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnlatj',&
               __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_fnlat,fnlat,procedureN//'_fnlat',ierr)
#else
          DEALLOCATE(fnlat, stat=ierr)
#endif
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnlat',&
               __LINE__,__FILE__)
       END IF
       
       DEALLOCATE(na,na_fnl,na_grp, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'deallocation problem',&
            __LINE__,__FILE__)

    END IF
    CALL summat(a,nstate,symmetrization=full,lsd=.TRUE.,gid=parai%cp_grp,&
         parent=only_parent)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE csmat
  ! ==================================================================
  SUBROUTINE store_ovlap(a,nstate)
    ! ==--------------------------------------------------------------==
    ! ==         STORE OVLAP MATRIX FOR LATER USE IN CROTWF           ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate
    REAL(real_8),INTENT(IN)                  :: a(nstate,*)
    INTEGER                                  :: i,j, ierr
    REAL(real_8)                             :: temp, selem
    CHARACTER(*), PARAMETER                  :: procedureN = 'store_nort'
    ! ==--------------------------------------------------------------==
    IF (.NOT. ALLOCATED(nort_ovlap)) THEN
       ALLOCATE(nort_ovlap(nstate,nstate), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate nort_ovlap',&
            __LINE__,__FILE__)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target enter data map(alloc:nort_ovlap)
#endif
    ELSE
       IF(SIZE(nort_ovlap).NE.nstate*nstate)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target exit data map(delete:nort_ovlap)
#endif
          DEALLOCATE(nort_ovlap, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate nort_ovlap',&
               __LINE__,__FILE__)
          ALLOCATE(nort_ovlap(nstate,nstate), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate nort_ovlap',&
               __LINE__,__FILE__)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target enter data map(alloc:nort_ovlap)
#endif
       END IF
    END IF
    IF(cntl%tlsd)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute parallel do simd private(i,j)
#else
       !$omp parallel private(i,j)
       !$omp do
#endif
       DO i=1,spin_mod%nsup
          DO j=1,i
             nort_ovlap(j,i)=a(j,i)
          END DO
       END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute parallel do simd private(i,j)
#else
       !$omp end do nowait
       !$omp  do
#endif
       DO i=spin_mod%nsup+1,nstate
          DO j=spin_mod%nsup+1,i
             nort_ovlap(j,i)=a(j,i)
          END DO
       END DO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp end parallel
#endif
    ELSE
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute parallel do simd private(i,j)
#else
       !$omp parallel do private(i,j)
#endif
       DO i=1,nstate
          DO j=1,i
             nort_ovlap(j,i)=a(j,i)
          END DO
       END DO
    END IF
    !lapack only needs upper or lower part
    CALL summat(nort_ovlap,nstate,symmetrization=cntl%use_elpa,lsd=.TRUE.,gid=parai%cp_grp,&
         parent=.FALSE.)
    !always check threshold.
    temp=0.0_real_8
    IF(cntl%tlsd)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute parallel do simd private(i,j,selem)&
       !$omp& reduction(max:temp)
#else
       !$omp parallel private(i,j,selem)reduction(max:temp)
       !$omp do
#endif
       DO i=1,spin_mod%nsup
          DO j=1,i-1
             selem=ABS(nort_ovlap(j,i))
             IF (temp.LT.selem) temp=selem
          ENDDO
       ENDDO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute parallel do simd private(i,j,selem)&
       !$omp& reduction(max:temp)
#else
       !$omp end do nowait
       !$omp do
#endif
       DO i=spin_mod%nsup+1,nstate
          DO j=spin_mod%nsup+1,i-1
             selem=ABS(nort_ovlap(j,i))
             IF (temp.LT.selem) temp=selem
          ENDDO
       ENDDO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp end parallel
#endif
    ELSE
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute parallel do simd private(i,j,selem)&
       !$omp& reduction(max:temp)
#else
       !$omp parallel do private(i,j,selem)reduction(max:temp)
#endif
       DO i=1,nstate
          DO j=1,i-1
             selem=ABS(nort_ovlap(j,i))
             IF (temp.LT.selem) temp=selem
          ENDDO
       ENDDO
    END IF
    nort_com%scond=temp

  END SUBROUTINE store_ovlap
  ! ==================================================================
  !DIR$ ATTRIBUTES FORCEINLINE::prepare_matrix
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  SUBROUTINE prepare_matrix(fnl_p,fnli,fnlj,qq_,ngh,ia_sum,ia_fnl,fnl_start,maxngh)
#else
  PURE SUBROUTINE prepare_matrix(fnl_p,fnli,fnlj,qq_,ngh,ia_sum,ia_fnl,fnl_start,maxngh)
#endif
    INTEGER,INTENT(IN)                       :: ngh,ia_sum,ia_fnl,fnl_start,maxngh
    REAL(real_8),INTENT(IN)                  :: fnl_p(ia_fnl,ngh,*),qq_(maxngh,*)
    REAL(real_8),INTENT(OUT)                 :: fnli(ia_sum,ngh,*),fnlj(ia_sum,ngh,*)
    INTEGER                                  :: iv,ia,jv
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp parallel do private(iv,ia)
#endif
    DO iv=1,ngh
       DO ia=1,ia_sum
          fnli(ia,iv,1)=fnl_p(ia+fnl_start,iv,1)
          fnlj(ia,iv,1)=0.0_real_8
       END DO
    END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp parallel do private(iv,jv,ia)
#endif
    DO iv=1,ngh
       DO jv=1,ngh
          IF (ABS(qq_(jv,iv)).GT.1.e-5_real_8) THEN
             DO ia=1,ia_sum
                fnlj(ia,iv,1)=fnlj(ia,iv,1)&
                     +qq_(jv,iv)*fnli(ia,jv,1)
             END DO
          END IF
       END DO
    END DO
  END SUBROUTINE prepare_matrix
  ! ==================================================================
END MODULE csmat_utils
