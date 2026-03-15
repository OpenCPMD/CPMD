#include "cpmd_global.h"
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
  include "mkl_lapack_omp_offload_lp64.f90"
#endif

  SUBROUTINE cpmd_dtrtri( UPLO, DIAG, N, A, LDA, INFO )
    USE kinds,                          ONLY: real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
    USE onemkl_lapack_omp_offload_lp64
    USE gpu
#endif
    IMPLICIT NONE

    CHARACTER(1)       DIAG, UPLO
    INTEGER            INFO, LDA, N
    REAL(real_8)       A( LDA, * )
    integer, save :: info_save=huge(0)
    if(info_save.eq.huge(0))then
       !$omp target enter data map(alloc:info_save)
    end if
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target enter data map(alloc:A(1:LDA,1:N))
    !$omp target update to(A(1:LDA,1:N)) if(update_first_to_gpu)
    !$omp target data use_device_addr(a,info_save)
    !$omp dispatch
#endif
    CALL dtrtri( UPLO, DIAG, N, A, LDA, INFO_save )
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end target data
    !$omp target update from(info_save)  
    !$omp target update from(A(1:LDA,1:N)) if(update_result_to_host)
    !$omp target exit data map(release:A(1:LDA,1:N))
#endif
    info=info_save
  END SUBROUTINE cpmd_dtrtri
  
  SUBROUTINE cpmd_dpotrf( UPLO, N, A, LDA, INFO )
    USE kinds,                          ONLY: real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
    USE onemkl_lapack_omp_offload_lp64
    USE gpu
#endif
    IMPLICIT NONE

    CHARACTER(1)       UPLO
    INTEGER            INFO, LDA, N
    REAL(real_8)       A( LDA, * )
    integer, save :: info_save=huge(0)
    if(info_save.eq.huge(0))then
       !$omp target enter data map(alloc:info_save)
    end if
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target enter data map(alloc:A(1:LDA,1:N))
    !$omp target update to(A(1:LDA,1:N)) if(update_first_to_gpu)
    !$omp target data use_device_addr(a,info_save)
    !$omp dispatch
#endif
    CALL dpotrf( UPLO, N, A, LDA, INFO_SAVE )
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end target data
    !$omp target update from(info_save)
    info=info_save
    !$omp target update from(A(1:LDA,1:N)) if(update_result_to_host)
    !$omp target exit data map(release:A(1:LDA,1:N))
#endif
  END SUBROUTINE cpmd_dpotrf

  SUBROUTINE cpmd_dsyevd(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO)
    USE kinds,                          ONLY: real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
    USE onemkl_lapack_omp_offload_lp64
    USE gpu
#endif
    IMPLICIT NONE
    CHARACTER(1)          jobz, uplo
    INTEGER,intent(out)::               info
    INTEGER,intent(in)::               lda, liwork, lwork, n
    INTEGER  ::             iwork( * )
    DOUBLE PRECISION::      a( lda, * ), w( * ), work( * )
    LOGICAL               work_space_query
    integer, save :: info_save=huge(0)

    if(info_save.eq.huge(0))then
       !$omp target enter data map(alloc:info_save)
    end if
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    work_space_query=.FALSE.
    IF(lwork.EQ.-1) work_space_query=.TRUE.
    !$omp target enter data map(alloc:A(1:LDA,1:N),W(1:N))
    IF(work_space_query)THEN
       !$omp target enter data map(alloc:work(1),iwork(1))
    ELSE
       !$omp target enter data map(alloc:work(1:lwork),iwork(1:liwork))
    END IF
    !$omp target update to(A(1:LDA,1:N)) if(update_first_to_gpu)
    !$omp target data use_device_addr(a,w,work,iwork,info_save)
    !$omp dispatch
#endif   
    CALL dsyevd(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO_save)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end target data
    !$omp target update from(info_save)
    IF(work_space_query)THEN
       !$omp target update from(work(1),iwork(1))
       !$omp target exit data map(release:work(1),iwork(1))
    ELSE
       !$omp target exit data map(release:work(1:lwork),iwork(1:liwork))
    END IF
    !$omp target update from(A(1:LDA,1:N),W(1:N)) if(update_result_to_host)
    !$omp target exit data map(release:A(1:LDA,1:N),W(1:N))
#endif
    info=info_save

  END SUBROUTINE cpmd_dsyevd

  
  SUBROUTINE cpmd_dsyevx(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, &
       WORK, LWORK, IWORK, IFAIL, INFO)
    USE kinds,                          ONLY: real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    USE onemkl_lapack_omp_offload_lp64
    USE gpu
#endif
    IMPLICIT NONE
    CHARACTER(1)          jobz, uplo, range
    INTEGER,intent(out)::               info, m
    INTEGER,intent(in)::               lda, lwork, n, ldz, il, iu
    INTEGER  ::             iwork( * ), ifail( * )
    DOUBLE PRECISION::      a( lda, * ), w( * ), work( * ), z( ldz, *)
    DOUBLE PRECISION,intent(in)::      vl, vu, abstol
    LOGICAL               work_space_query
    integer, save :: info_save=huge(0), m_save=huge(0)

    if(info_save.eq.huge(0))then
       !$omp target enter data map(alloc:info_save,m_save)
    end if
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    work_space_query=.FALSE.
    IF(lwork.EQ.-1) work_space_query=.TRUE.
    !$omp target enter data map(alloc:A(1:LDA,1:N),W(1:N),Z(1:LDZ,1:N),IFAIL(1:N))
    IF(work_space_query)THEN
       !$omp target enter data map(alloc:work(1),iwork(1))
    ELSE
       !$omp target enter data map(alloc:work(1:lwork),iwork(1:5*N))
    END IF
    !$omp target update to(A(1:LDA,1:N)) if(update_first_to_gpu)

    !$omp target data use_device_addr(a,w,work,iwork,ifail,z,m_save,info_save)
    !$omp dispatch
#endif
    CALL dsyevx(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M_SAVE, W, Z, LDZ, &
       WORK, LWORK, IWORK, IFAIL, INFO_SAVE)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end target data
    !$omp target update from(info_save,m_save)
    IF(work_space_query)THEN
       !$omp target update from(work(1),iwork(1))
       !$omp target exit data map(release:work(1),iwork(1))
    ELSE
       !$omp target exit data map(release:work(1:lwork),iwork(1:5*N))
    END IF
    !$omp target update from(W(1:N),Z(1:LDZ,1:N),IFAIL(1:N)) if(update_result_to_host)
    !$omp target exit data map(release:A(1:LDA,1:N),W(1:N),Z(1:LDZ,1:N),IFAIL(1:N))
#endif
    info=info_save
    m=m_save

  END SUBROUTINE cpmd_dsyevx
