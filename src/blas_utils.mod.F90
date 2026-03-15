#include "cpmd_global.h"
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
  include "mkl_blas_omp_offload_lp64.f90"
#endif


SUBROUTINE cpmd_dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  USE kinds,                          ONLY: real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
  USE onemkl_blas_omp_offload_lp64
  USE gpu
#endif
  IMPLICIT NONE
  CHARACTER(1)                :: TRANSA, TRANSB
  INTEGER                     :: M, N, K, LDA, LDB, LDC
  REAL(real_8), INTENT(IN)    :: ALPHA, A(*), B(*), BETA
  REAL(real_8), INTENT(INOUT) :: C(*)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  INTEGER                     :: ka, kb, la, lb
  INTEGER                     :: len_a, len_b, len_c
  IF( TRANSA .EQ. 'N' .OR. TRANSA .eq. 'n' )THEN
     ka = k
     la = m
  ELSE
     ka = m
     la = k
  END IF

  IF( TRANSB .EQ. 'N' .OR. TRANSB .eq. 'n' )THEN
     kb = n
     lb = k
  ELSE
     kb = k
     lb = n
  END IF

  len_a = LDA * (ka - 1) + la
  len_b = LDB * (kb - 1) + lb
  len_c = LDC * (n  - 1) + m
  !$omp target enter data map(alloc:A(1:len_a),B(1:len_b),C(1:len_c))
  !$omp target update to(A(1:len_a)) if(update_first_to_gpu)
  !$omp target update to(B(1:len_b)) if(update_second_to_gpu)
  !$omp target update to(C(1:len_c)) if(update_third_to_gpu)
  !$omp dispatch
#endif
  CALL dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target update from (C(1:len_c)) if(update_result_to_host)
  !$omp target exit data map(release:A(1:len_a),B(1:len_b),C(1:len_c))
#endif
END SUBROUTINE cpmd_dgemm

SUBROUTINE cpmd_dgemmt(UPLO, TRANSA, TRANSB, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  USE kinds,                          ONLY: real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
  USE onemkl_blas_omp_offload_lp64
  USE gpu
#endif
  IMPLICIT NONE
  CHARACTER(1)                :: UPLO, TRANSA, TRANSB
  INTEGER                     :: N, K, LDA, LDB, LDC
  REAL(real_8), INTENT(IN)    :: ALPHA, A(*), B(*), BETA
  REAL(real_8), INTENT(INOUT) :: C(*)
#if defined(_HAS_DGEMMT)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  INTEGER                     :: ka, kb, la, lb
  INTEGER                     :: len_a, len_b, len_c
  IF( TRANSA .EQ. 'N' .OR. TRANSA .eq. 'n' )THEN
     ka = k
     la = n
  ELSE
     ka = n
     la = k
  END IF

  IF( TRANSB .EQ. 'N' .OR. TRANSB .eq. 'n' )THEN
     kb = n
     lb = k
  ELSE
     kb = k
     lb = n
  END IF
  len_a = LDA * (ka - 1) + la
  len_b = LDB * (kb - 1) + lb
  len_c = LDC * (n  - 1) + n
  !$omp target enter data map(alloc:A(1:len_a),B(1:len_b),C(1:len_c))
  !$omp target update to(A(1:len_a)) if(update_first_to_gpu)
  !$omp target update to(B(1:len_b)) if(update_second_to_gpu)
  !$omp target update to(C(1:len_c)) if(update_third_to_gpu)
  !$omp dispatch
#endif
  CALL dgemmt(UPLO, TRANSA, TRANSB, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target update from (C(1:len_c)) if(update_result_to_host)
  !$omp target exit data map(release:A(1:len_a),B(1:len_b),C(1:len_c))
#endif
#endif
END SUBROUTINE cpmd_dgemmt

SUBROUTINE cpmd_dsymm(SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  USE kinds,                          ONLY: real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
  USE onemkl_blas_omp_offload_lp64
  USE gpu
#endif
  IMPLICIT NONE
  CHARACTER(1)                :: SIDE, UPLO
  INTEGER                     :: M, N, LDA, LDB, LDC
  REAL(real_8), INTENT(IN)    :: ALPHA, A(*), B(*), BETA
  REAL(real_8), INTENT(INOUT) :: C(*)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  INTEGER                     :: ka, la
  INTEGER                     :: len_a, len_b, len_c
  IF( SIDE .EQ. 'L' .OR. SIDE .eq. 'l' )THEN
     ka = m
     la = m
  ELSE
     ka = n
     la = n
  END IF
  len_a = LDA * (ka - 1) + la
  len_b = LDB * (n  - 1) + m
  len_c = LDC * (n  - 1) + m 

  !$omp target enter data map(alloc:A(1:len_a),B(1:len_b),C(1:len_c))
  !$omp target update to(A(1:len_a)) if(update_first_to_gpu)
  !$omp target update to(B(1:len_b)) if(update_second_to_gpu)
  !$omp target update to(C(1:len_c)) if(update_third_to_gpu)
  !$omp dispatch
#endif
  CALL dsymm(SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target update from (C(1:len_c)) if(update_result_to_host)
  !$omp target exit data map(release:A(1:len_a),B(1:len_b),C(1:len_c))
#endif
END SUBROUTINE cpmd_dsymm

SUBROUTINE cpmd_dsyrk(UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC)
  USE kinds,                          ONLY: real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
  USE onemkl_blas_omp_offload_lp64
  USE gpu
#endif
  IMPLICIT NONE
  CHARACTER(1)                :: UPLO, TRANS
  INTEGER                     :: N, K, LDA, LDC
  REAL(real_8), INTENT(IN)    :: ALPHA, A(*), BETA
  REAL(real_8), INTENT(INOUT) :: C(*)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  INTEGER                     :: ka, la
  INTEGER                     :: len_a, len_c

  IF( TRANS .EQ. 'N' .OR. TRANS .eq. 'n' )THEN
     ka = k
     la = n
  ELSE
     ka = n
     la = k
  END IF

  len_a = LDA * (ka - 1) + la
  len_c = LDC * (n  - 1) + n
  !$omp target enter data map(alloc:A(1:len_a),C(1:len_c))
  !$omp target update to(A(1:len_a)) if(update_first_to_gpu)
  !$omp target update to(C(1:len_c)) if(update_third_to_gpu)
  !$omp dispatch
#endif
  CALL dsyrk(UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target update from (C(1:len_c)) if(update_result_to_host)
  !$omp target exit data map(release:A(1:len_a),C(1:len_c))
#endif
END SUBROUTINE cpmd_dsyrk

SUBROUTINE cpmd_dtrmm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
  USE kinds,                          ONLY: real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
  USE onemkl_blas_omp_offload_lp64
  USE gpu
#endif
  IMPLICIT NONE
  CHARACTER(1)               :: SIDE, UPLO, TRANSA, DIAG
  INTEGER                    :: M, N, LDA, LDB
  REAL(real_8),INTENT(IN)    :: ALPHA, A(*)
  REAL(real_8),INTENT(INOUT) :: B(*)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  INTEGER                    :: k
  INTEGER                    :: len_a, len_b
  IF( SIDE .EQ. 'L' .OR. SIDE .eq. 'l' )THEN
     k = m
  ELSE
     k = n
  END IF

  len_a = LDA * (k - 1) + k
  len_b = LDB * (n - 1) + m
  !$omp target enter data map(alloc:A(1:len_a),B(1:len_b))
  !$omp target update to(A(1:len_a)) if(update_first_to_gpu)
  !$omp target update to(B(1:len_b)) if(update_second_to_gpu)
  !$omp dispatch
#endif
  CALL dtrmm(SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target update from (B(1:len_b)) if(update_result_to_host)
  !$omp target exit data map(release:A(1:len_a),B(1:len_b))
#endif
END SUBROUTINE cpmd_dtrmm

SUBROUTINE cpmd_dger(M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
  USE kinds,                          ONLY: real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  USE onemkl_blas_omp_offload_lp64
  USE gpu
#endif
  INTEGER                    :: M, N, INCX, INCY, LDA
  REAL(real_8),INTENT(IN)    :: ALPHA, X(*), Y(*)
  REAL(real_8),INTENT(INOUT) :: A(*)
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
  INTEGER                    :: len_x, len_y, len_a

  len_x = 1 + ( m - 1 ) * abs( INCX )
  len_y = 1 + ( n - 1 ) * abs( INCY )
  len_a = LDA * (n - 1) + m
  !$omp target enter data map(alloc:A(1:len_a),X(1:len_x),Y(1:len_y))
  !$omp target update to(X(1:len_x)) if(update_first_to_gpu)
  !$omp target update to(Y(1:len_y)) if(update_second_to_gpu)
  !$omp target update to(A(1:len_a)) if(update_third_to_gpu)
#endif
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp dispatch
#endif
  CALL dger(M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target update from (A(1:len_a)) if(update_result_to_host)
  !$omp target exit data map(release:X(1:len_x),Y(1:len_y),A(1:len_a))
#endif
END SUBROUTINE cpmd_dger

SUBROUTINE cpmd_dcopy(N,DX,INCX,DY,INCY)
  USE kinds,                          ONLY: real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
  USE onemkl_blas_omp_offload_lp64
  USE gpu
#endif
  IMPLICIT NONE
  INTEGER                      :: N, INCX, INCY
  REAL(real_8), INTENT(IN)     :: DX(*)
  REAL(REAL_8), INTENT(OUT)    :: DY(*)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target enter data map(alloc:DX(1:INCX*N),DY(1:INCY*N))
  !$omp target update to(DX(1:INCX*N)) if(update_first_to_gpu)
  !$omp dispatch
#endif
  call dcopy(n,dx,incx,dy,incy)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target update from(DY(1:INCY*N)) if(update_result_to_host)
  !$omp target exit data map(release:DX(1:INCX*N),DY(1:INCY*N))
#endif
END SUBROUTINE cpmd_dcopy
SUBROUTINE cpmd_zgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
  USE kinds,                          ONLY: real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
  USE onemkl_blas_omp_offload_lp64
  USE gpu
#endif
  IMPLICIT NONE
  CHARACTER(1)                  :: TRANSA, TRANSB
  INTEGER                       :: M, N, K, LDA, LDB, LDC
  COMPLEX(real_8),INTENT(IN)    :: ALPHA, A(*), B(*), BETA
  COMPLEX(real_8),INTENT(INOUT) :: C(*)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  INTEGER                       :: ka, kb, la, lb
  INTEGER                       :: len_a, len_b, len_c
  IF( TRANSA .EQ. 'N' .OR. TRANSA .eq. 'n' )THEN
     ka = k
     la = m
  ELSE
     ka = m
     la = k
  END IF

  IF( TRANSB .EQ. 'N' .OR. TRANSB .eq. 'n' )THEN
     kb = n
     lb = k
  ELSE
     kb = k
     lb = n
  END IF
  len_a = LDA * (ka - 1) + la
  len_b = LDB * (kb - 1) + lb
  len_c = LDC * (n  - 1) + m

  !$omp target enter data map(alloc:A(1:len_a),B(1:len_b),C(1:len_c))
  !$omp target update to(A(1:len_a)) if(update_first_to_gpu)
  !$omp target update to(B(1:len_b)) if(update_second_to_gpu)
  !$omp target update to(C(1:len_c)) if(update_third_to_gpu)
  !$omp dispatch
#endif
  CALL zgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target update from (C(1:len_c)) if(update_result_to_host)
  !$omp target exit data map(release:A(1:len_a),B(1:len_b),C(1:len_c))
#endif
END SUBROUTINE cpmd_zgemm

SUBROUTINE cpmd_zgemv(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
  USE kinds,                          ONLY: real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
  USE onemkl_blas_omp_offload_lp64
  USE gpu
#endif
  IMPLICIT NONE
  CHARACTER(1)                  :: TRANS
  INTEGER                       :: M, N, LDA, INCX, INCY
  COMPLEX(real_8),INTENT(IN)    :: ALPHA, X(*), BETA, A(*)
  COMPLEX(real_8),INTENT(INOUT) :: Y(*)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  INTEGER                       :: len_x, len_y, len_a

  IF( TRANS .EQ. 'N' .OR. TRANS .eq. 'n' )THEN
     len_x = 1 + ( n - 1 )*abs( INCX )
     len_y = 1 + ( m - 1 )*abs( INCY )
  ELSE
     len_x = 1 + ( m - 1 )*abs( INCX )
     len_y = 1 + ( n - 1 )*abs( INCY ) 
  END IF

  len_a = LDA * ( m - 1 ) * n + m
  !$omp target enter data map(alloc:A(1:len_a),X(1:len_x),Y(1:len_y))
  !$omp target update to(A(1:len_a)) if(update_first_to_gpu)
  !$omp target update to(X(1:len_x)) if(update_second_to_gpu)
  !$omp target update to(Y(1:len_y)) if(update_third_to_gpu)
  !$omp dispatch
#endif  
  CALL zgemv(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target update from (Y(1:len_y)) if(update_result_to_host)
  !$omp target exit data map(release:X(1:len_x),A(1:len_a),Y(1:len_y))
#endif
END SUBROUTINE cpmd_zgemv

!DIR$ ATTRIBUTES FORCEINLINE::cpmd_ddot_2
FUNCTION cpmd_ddot_2(N,x,y) RESULT(ddot)
  USE kinds, ONLY: real_8
  USE omp_lib
  !, ONLY: omp_is_inital_device
  INTEGER, INTENT(IN) :: n
  REAL(real_8),INTENT(IN)  :: x(*), y(*)
  REAL(real_8) :: ddot
!$omp declare target  
  ddot=0._real_8
  DO i=1,n
     ddot=ddot+x(i)*y(i)
  END DO
END FUNCTION cpmd_ddot_2

!DIR$ ATTRIBUTES FORCEINLINE::cpmd_ddot_1
FUNCTION cpmd_ddot_1(N,x) RESULT(ddot)
  USE kinds, ONLY: real_8
  USE omp_lib
  !, ONLY: omp_is_inital_device
  INTEGER, INTENT(IN) :: n
  REAL(real_8),INTENT(IN)  :: x(*)
  REAL(real_8) :: ddot
  !$omp declare target
  ddot=0._real_8
  DO i=1,n
     ddot=ddot+x(i)**2
  END DO
END FUNCTION cpmd_ddot_1
