#include "cpmd_global.h"

MODULE utils
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE gpu
  USE kinds,                           ONLY: int_8,&
                                             real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE reshaper,                        ONLY: reshape_inplace
  USE timer,                           ONLY: tiset,&
                                             tihalt
  USE system,                          ONLY: parap,&
                                             cntl,&
                                             maxsys
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: unitmx
  PUBLIC :: simpsn
  PUBLIC :: zclean
  PUBLIC :: icopy
  PUBLIC :: symma
  PUBLIC :: symmat_pack
  PUBLIC :: symmat_unpack
  PUBLIC :: numcpus
  PUBLIC :: invmat
  PUBLIC :: inversemat
  PUBLIC :: elpa_driver
  PUBLIC :: dsyevd_driver
  PUBLIC :: dsyevx_driver
  PUBLIC :: dspevy
  PUBLIC :: zgthr
  PUBLIC :: zgthr_no_omp
  PUBLIC :: zsctr_no_omp
  PUBLIC :: zclean_k
  PUBLIC :: setkwf
  PUBLIC :: fskip
  PUBLIC :: determ
  PUBLIC :: rtrans
  PUBLIC :: rmatmov
  PUBLIC :: izamax
  PUBLIC :: idamin
  PUBLIC :: dgive
  PUBLIC :: zgive
  PUBLIC :: nxxfun
  PUBLIC :: sbes0
  PUBLIC :: fc4
  PUBLIC :: print_debug_ions
#if defined(__SR11KIBM)
  PUBLIC :: dzamax
#endif
#if defined(__HAS_EXTERNAL_IZAMAX)
  INTEGER, EXTERNAL :: izamax
#endif
#if defined(__HAS_EXTERNAL_IDAMIN)
  INTEGER, EXTERNAL :: idamin
#endif

CONTAINS

  ! ==================================================================
  SUBROUTINE unitmx(a,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: a(n,n)

    INTEGER                                  :: i

    CALL zeroing(a)!,n*n)
    !$omp parallel do private(I)
    DO i=1,n
       a(i,i)=1._real_8
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE unitmx
  ! ==================================================================
  SUBROUTINE simpsn(n,inte,sum)
    ! ==--------------------------------------------------------------==
    ! == COMPUTES ONE-DIMENSIONAL INTEGRALS BY THE SIMPSON METHOD     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: inte(n), sum

    REAL(real_8), PARAMETER :: c1 = 109._real_8/48._real_8, &
      c2 = -5._real_8/48._real_8, c3 = 63._real_8/48._real_8, &
      c4 = 49._real_8/48._real_8

    INTEGER                                  :: i

! ==--------------------------------------------------------------==

    inte(1)   = inte(1)*c1
    inte(2)   = inte(2)*c2
    inte(3)   = inte(3)*c3
    inte(4)   = inte(4)*c4
    inte(n-1) = inte(n-1)*c1
    inte(n-2) = inte(n-2)*c2
    inte(n-3) = inte(n-3)*c3
    inte(n-4) = inte(n-4)*c4
    sum=0._real_8
    !$omp parallel do private(I) reduction(+:SUM)
    DO i=1,n-1
       sum=sum+inte(i)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE simpsn
  ! ==================================================================
  SUBROUTINE zclean(a,n,ngw)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: n, ngw
    COMPLEX(real_8), TARGET, INTENT(INOUT)   :: a(ngw,n)

    INTEGER                                  :: i
    REAL(real_8), POINTER __CONTIGUOUS       :: pa(:,:,:)

    CALL reshape_inplace(a, (/2, ngw, n/), pa)

    IF (ngw.GT.0) THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target update to(a) if(update_first_to_gpu)
       !$omp target teams distribute parallel do simd &
#else
       !$omp parallel do simd &
       
#endif
       !$omp& private(I)
       DO i=1,n
          pa(2,1,i)=0._real_8
       ENDDO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target update from(a) if(update_result_to_host)
#endif
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE zclean

  ! ==================================================================
  ! LSA#endif
  SUBROUTINE symma(a,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: a(n,n)

    INTEGER                                  :: i, j

    !$omp parallel do private(I,J) schedule(static,1)
    DO i=1,n
       DO j=i+1,n
          a(i,j)=0.5_real_8*(a(i,j)+a(j,i))
          a(j,i)=a(i,j)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE symma
  ! ==================================================================
  SUBROUTINE symmat_pack(a,packeda,n,n1,n2)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: n,n1,n2
    REAL(real_8),INTENT(IN)                  :: a(n,n)
    REAL(real_8),INTENT(OUT) __CONTIGUOUS    :: packeda(:)

    INTEGER                                  :: offset, i, j, k
    !for lsd=>n=nstate,n1=spind_mod%nsup,n2=spin_mod%nsdown
    !else=>n=nstate,n1=nstate,n2=0
    ! ==--------------------------------------------------------------==
    offset=n1*(n1+1)/2
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target enter data map(alloc:a,packeda)
    !$omp target update to(a) if(update_first_to_gpu)
    !$omp target teams &
#else
    !$omp parallel &
#endif
    !$omp& private(i,j,k)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp distribute 
#else
    !$omp do schedule(static)
#endif
    DO i=0,n1-1
       k=i*(i-1)/2+i+1
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do
#endif
       DO j=k,k+i
          packeda(j)=a(j-k+1,i+1)
       END DO
    END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp distribute 
#else
    !$omp end do nowait
    !$omp do schedule(static)
#endif
    DO i=0,n2-1
       k=i*(i-1)/2+i+1+offset
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do
#endif
       DO j=k,k+i
          packeda(j)=a(n1+j-k+1+i+1,n1+i+1)
       END DO
    END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end target teams
    !$omp target update from(packeda) if(update_result_to_host)
    !$omp target exit data map(release:a,packeda)
#else
    !$omp end parallel
#endif
    RETURN
  END SUBROUTINE symmat_pack
  ! ==================================================================
  SUBROUTINE symmat_unpack(a,packeda,n,n1,n2,full)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: n,n1,n2
    REAL(real_8),INTENT(OUT)                 :: a(n,n)
    REAL(real_8),INTENT(IN)  __CONTIGUOUS    :: packeda(:)
    LOGICAL,INTENT(IN)                       :: full

    INTEGER                                  :: offset, i, j, k
    !for lsd=>n=nstate,n1=spind_mod%nsup,n2=spin_mod%nsdown
    !else=>n=nstate,n1=nstate,n2=0
    ! ==--------------------------------------------------------------==
    offset=n1*(n1+1)/2
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target enter data map(alloc:a,packeda)
    !$omp target update to(packeda) if(update_first_to_gpu)
    !$omp target teams &
#else
    !$omp& parallel &
#endif
    !$omp& private(i,j,k)

#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp distribute 
#else
    !$omp do schedule(static)
#endif
    DO i=0,n1-1
       k=i*(i-1)/2+i+1
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do
#endif
       DO j=k,k+i
          a(j-k+1,i+1)=packeda(j)
       END DO      
    END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp distribute
#else
    !$omp end do nowait
    !$omp do schedule(static)
#endif
    DO i=0,n2-1
       k=i*(i-1)/2+i+1+offset
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do
#endif
       DO j=k,k+i
          a(n1+j-k+1+i+1,n1+i+1)=packeda(j)
       END DO
    END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end target teams
#else
    !$omp end parallel
#endif
    IF(full)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams&
#else
       !$omp parallel &
#endif
       !$omp& private(i,j)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp distribute parallel do
#else
       !$omp do schedule(static)
#endif
       DO i=1,n1
          DO j=i,n1
             a(j,i)=a(i,j)
          END DO
       END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp distribute parallel do
#else
       !$omp end do nowait
       !$omp do schedule(static)
#endif
       DO i=n1+1,n
          DO j=i,n
             a(j,i)=a(i,j)
          END DO
       END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp end target teams
#else
       !$omp end parallel
#endif
    END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target update from(a) if(update_result_to_host)
    !$omp target exit data map(release:a,packeda)
#endif
  END SUBROUTINE symmat_unpack
  ! ==================================================================
  SUBROUTINE numcpus(ncpus)
    ! ==--------------------------------------------------------------==
    ! == Number of CPU for a distributed machine                      ==
    ! == Used with OpenMP version.                                    ==
    ! == First call Give the number of CPU and store locally in MCPUS ==
    ! == Other calls Use saved value in MCPUS variable                ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ncpus

    INTEGER                                  :: idummy
    INTEGER, SAVE                            :: mcpus = 0

    !$  integer :: omp_get_num_threads

    ! ==--------------------------------------------------------------==
    IF (mcpus.NE.0) THEN
       ncpus=mcpus
    ELSE
       ! OpenMP version
       ! We forcedly disable nested paralellism.
       ncpus=1
       !$omp parallel shared(NCPUS)
       !$      NCPUS=OMP_GET_NUM_THREADS()
#ifndef _HASNT_OMP_SET_NESTED
       !$      CALL OMP_SET_NESTED(.FALSE.)
#endif
       !$      CALL omp_set_max_active_levels( 1 )
       !$omp end parallel
       mcpus=ncpus
#if defined(_HAS_FFT_FFTW3)
       !$      CALL DFFTW_INIT_THREADS(IDUMMY)
       !$      CALL DFFTW_PLAN_WITH_NTHREADS(MCPUS)
#endif
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE numcpus
  ! ==================================================================
  SUBROUTINE invmat(n,a,b,info)
    ! ==--------------------------------------------------------------==
    ! == Inverse Matrix A(N,N)                                        ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:  N Dimension                                          ==
    ! ==         A(N,N) Matrix                                        ==
    ! == OUTPUT: A(N,N) inverse matrix of A                           ==
    ! ==         B is a work array                                    ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: a(n,n), b(n,n)
    INTEGER                                  :: info

    INTEGER                                  :: lb

    IF (n.EQ.1) THEN
       a(1,1)=1._real_8/a(1,1)
    ELSE
       lb=(n-1)*n
       ! Compute an LU factorization
       CALL dgetrf(n,n,a,n,b(1,1),info)
       IF (info.EQ.0) THEN
          ! Compute the inverse of a matrix using the LU factorization
          CALL dgetri(n,a,n,b(1,1),b(1,2),lb,info)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE invmat
  ! ==================================================================
  SUBROUTINE inversemat(n,a,lda,b,info)
    ! ==--------------------------------------------------------------==
    ! == Inverse Matrix A(LDA,N)                                      ==
    ! ==--------------------------------------------------------------==
    ! == INPUT:  N   Dimension                                        ==
    ! ==         LDA Leading dimension of A                           ==
    ! ==         A(LDA,N) Matrix                                      ==
    ! == OUTPUT: A(LDA,N) inverse matrix of A                         ==
    ! ==         B is a work array                                    ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n, lda
    REAL(real_8)                             :: a(lda,n), b(n,n)
    INTEGER                                  :: info

    INTEGER                                  :: lb

    IF (n.EQ.1) THEN
       a(1,1)=1._real_8/a(1,1)
    ELSE
       lb=(n-1)*n
       ! Compute an LU factorization
       CALL dgetrf(n,n,a,lda,b(1,1),info)
       IF (info.EQ.0) THEN
          ! Compute the inverse of a matrix using the LU factorization
          CALL dgetri(n,a,lda,b(1,1),b(1,2),lb,info)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE inversemat
  ! ==----------------------------------------------------------------==
  SUBROUTINE elpa_driver(a,w,n)
    ! ==--------------------------------------------------------------==
    ! == PARALLEL ELPA EIGENVALUE SOLVER                              ==
    ! == ASUMES GLOBALLY REPLICATED SYMMETRIC MATRIX A, UPPER AND     ==
    ! == LOWER PART OF A MUST BE SET                                  ==
    ! == RETURNS GLOBALLY REPLICATED EIGENVECTORS OF A IN A           ==
    ! == CAN HANDLE TO SIZES (MAINLY FOR SPINUP AND DOWN)             ==
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, Erlangen
    ! Date August 2019
#ifdef _HAS_LIBELPA
    USE elpa
    USE elpa_utils,                           ONLY: elpa_ob_create,&
                                                    elpa_get_eigvect,&
                                                    elpa_set_mat,&
                                                    elpa_autotune
    IMPLICIT NONE
#endif
    INTEGER,INTENT(IN)                       :: n
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: w(:), a(:,:)
#ifdef _HAS_LIBELPA
    INTEGER                                  :: i, j, id, success, isub
    INTEGER(int_8)                           :: il_aux(2), il_c(2)
    LOGICAL                                  :: create
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: aux(:,:), c(:,:)
#else
    REAL(real_8), ALLOCATABLE                :: aux(:,:), c(:,:)
#endif
    CLASS(ELPA_T),SAVE,POINTER               :: el1, el2
    INTEGER, SAVE, ALLOCATABLE               :: idx1(:,:,:), idx2(:,:,:)
    INTEGER, SAVE                            :: na_row(2)=0, na_col(2)=0, n_loc(2)=0
#endif
    CHARACTER(*),PARAMETER                   :: procedureN='elpa_driver'
#ifdef _HAS_LIBELPA
    LOGICAL, SAVE                            :: autotune1_done=.FALSE., autotune2_done=.FALSE.
    create=.FALSE.
    IF(n_loc(1).EQ.0)THEN
       id=1
       create=.TRUE.
    ELSE IF(n.EQ.n_loc(1))THEN
       id=1
    ELSE IF(n_loc(2).EQ.0) THEN
       id=2
       create=.TRUE.
    ELSE IF(n.EQ.n_loc(2))THEN
       id=2
    ELSE
       CALL stopgm(procedureN,'problem with ELPA driver', &
            __LINE__,__FILE__)
    END IF
    n_loc(id)=n
    IF(create)THEN
       IF(id.EQ.1)THEN
          CALL elpa_ob_create(n_loc(id),el1,idx1,na_row(id),na_col(id))
       ELSEIF(id.EQ.2)THEN
          CALL elpa_ob_create(n_loc(id),el2,idx2,na_row(id),na_col(id))
       ELSE
          CALL stopgm(procedureN,'problem with ELPA driver', &
               __LINE__,__FILE__)
       END IF
    END IF
    IF(na_row(id).GT.0)THEN
       il_c(1)=na_row(id)
       il_c(2)=na_col(id)
       il_aux=il_c
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_c,c,procedureN//'_c',ierr)
#else
       ALLOCATE(c(il_c(1),il_c(2)),STAT = ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_aux,aux,procedureN//'_aux',ierr)
#else
       ALLOCATE(aux(il_aux(1),il_aux(2)),STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       IF(id.EQ.1)THEN
          CALL elpa_set_mat(a,c,idx1,na_col(id),na_row(id))
          IF(cntl%use_elpa_autotune.AND..NOT.autotune1_done) THEN
             call tiset(procedureN//'_autotune',isub)
             call elpa_autotune(c,w,aux,el1)
             call tihalt(procedureN//'_autotune',isub)
             autotune1_done=.TRUE.
          END IF
          CALL el1%eigenvectors(c, w, aux, success)
          CALL elpa_get_eigvect(aux,a,idx1,na_col(id),na_row(id))
       ELSEIF(id.EQ.2)THEN
          CALL elpa_set_mat(a,c,idx2,na_col(id),na_row(id))
          IF(cntl%use_elpa_autotune.AND..NOT.autotune2_done) THEN
             call tiset(procedureN//'_autotune',isub)
             call elpa_autotune(c,w,aux,el1)
             call tihalt(procedureN//'_autotune',isub)
             autotune2_done=.TRUE.
          END IF
          CALL el2%eigenvectors(c, w, aux, success)
          CALL elpa_get_eigvect(aux,a,idx2,na_col(id),na_row(id))
       ELSE
          CALL stopgm(procedureN,'problem with ELPA driver', &
               __LINE__,__FILE__)
       END IF
       IF(success.NE.0) CALL stopgm(procedureN,'ELPA: diagonalization failed', &
            __LINE__,__FILE__)

#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_aux,aux,procedureN//'_aux',ierr)
#else
       DEALLOCATE(aux,STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_c,c,procedureN//'_c',ierr)
#else
       DEALLOCATE(c1,STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
#else
    CALL stopgm(procedureN,'ELPA not installed', &
         __LINE__,__FILE__)
#endif
  END SUBROUTINE elpa_driver
  ! ==================================================================
  SUBROUTINE dsyevd_driver(iopt,a,w,n)
    ! ==--------------------------------------------------------------==
    ! == DIAGONALIZATION ROUTINE: FOLLOW THE ESSL CONVENTION          ==
    ! == DIVIDE AND CONQUER UNPACKED MEMORY INTENSE VARIANT OF DSPEVY ==
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, Erlangen
    ! Date March 2019

    INTEGER,INTENT(IN)                       :: iopt, n
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: w(:), a(:,:)
    !local
    INTEGER                                  :: dummy_int(1)
    INTEGER(int_8)                           :: il_work(1), il_iwork(1)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: work(:)
    INTEGER, POINTER __CONTIGUOUS            :: iwork(:)
#else
    REAL(real_8), ALLOCATABLE                :: work(:)
    INTEGER, ALLOCATABLE                     :: iwork(:)
#endif

    REAL(real_8)                             :: dummy_real(1)
    CHARACTER(1)                             :: jobz, uplo
    CHARACTER(*),PARAMETER                   :: procedureN='dsyevd_driver'
    INTEGER                                  :: info, ierr
    integer, save                            :: n_save=0, il_work_save=0, il_iwork_save=0, iopt_save=0
    IF(iopt.EQ.0)THEN
       jobz='N'
       uplo='L'
    ELSEIF(iopt.EQ.1)THEN
       jobz='V'
       uplo='L'
    ELSEIF(iopt.EQ.20)THEN
       jobz='N'
       uplo='U'
    ELSEIF(iopt.EQ.21)THEN
       jobz='V'
       uplo='U'
    END IF

    il_work(1)=-1
    il_iwork=-1
    IF(iopt_save.eq.iopt) then
       if(n_save.eq.n)then
          il_work(1)=il_work_save
          il_iwork(1)=il_iwork_save
       end if
    end IF
    if(il_work(1).eq.-1)then
       !workspace query
       CALL cpmd_dsyevd(jobz,uplo,n,a,n,w,dummy_real,INT(il_work(1)),dummy_int,INT(il_iwork(1))&
            ,info)
       il_work(1)=INT(dummy_real(1))
       il_iwork=dummy_int(1)
       il_work_save=il_work(1)
       il_iwork_save=il_iwork(1)
       n_save=n
       iopt_save=iopt
    end if
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_work,work,procedureN//'_work',ierr)
#else
    ALLOCATE(work(il_work(1)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_iwork,iwork,procedureN//'_iwork',ierr)
#else
    ALLOCATE(iwork(il_iwork(1)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    !actual calculation
    CALL cpmd_dsyevd(jobz,uplo,n,a,n,w,work,INT(il_work(1)),iwork,INT(il_iwork(1)),info)
    IF (info.NE.0) CALL stopgm(procedureN,'FAILED TO DIAGONALIZE',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_iwork,iwork,procedureN//'_iwork',ierr)
#else
    DEALLOCATE(work,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_work,work,procedureN//'_work',ierr)
#else
    DEALLOCATE(work,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
  END SUBROUTINE dsyevd_driver
  ! ==================================================================
  SUBROUTINE dsyevx_driver(iopt,a,z,w,n,il,iu,abstol)
    ! ==--------------------------------------------------------------==
    ! == DIAGONALIZATION ROUTINE: FOLLOW THE ESSL CONVENTION          ==
    ! == UNPACKED MEMORY INTENSE VARIANT OF DSPEVY, PARALLELIZATION   ==
    ! == POSSIBLE SELECTING DIFFERENT EIGEVECTOR INDIZES ON EACH PROC ==
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, Erlangen
    ! Date March 2019

    INTEGER,INTENT(IN)                       :: iopt, n, il, iu
    REAL(real_8),INTENT(IN)                  :: abstol
    REAL(real_8),INTENT(INOUT)               :: w(:), a(:,:), z(:,:)
    !local
    INTEGER                                  :: m
    INTEGER(int_8)                           :: il_work(1), il_iwork(1), il_ifail(1)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: work(:)
#else
    REAL(real_8), ALLOCATABLE                :: work(:)
#endif
#ifdef _USE_SCRATCHLIBRARY
    INTEGER, POINTER , CONTIGUOUS            :: iwork(:), ifail(:)
#else
    INTEGER, ALLOCATABLE                     :: iwork(:), ifail(:)
#endif
    REAL(real_8)                             :: dummy_real(1)
    CHARACTER(1)                             :: jobz, uplo, range
    CHARACTER(*),PARAMETER                   :: procedureN='dsyevr_driver'
    INTEGER                                  :: info, ierr
    INTEGER, SAVE                            :: n_save=0, il_work_save=0, iopt_save=0
    
    IF(iu-il+1.LE.0)RETURN
    IF(iopt.EQ.0)THEN
       jobz='N'
       uplo='L'
       range='I'
       IF(il.EQ.1.AND.iu.EQ.n)range='A'
    ELSEIF(iopt.EQ.1)THEN
       jobz='V'
       uplo='L'
       range='I'
       IF(il.EQ.1.AND.iu.EQ.n)range='A'
    ELSEIF(iopt.EQ.20)THEN
       jobz='N'
       uplo='U'
       range='I'
       IF(il.EQ.1.AND.iu.EQ.n)range='A'
    ELSEIF(iopt.EQ.21)THEN
       jobz='V'
       uplo='U'
       range='I'
       IF(il.EQ.1.AND.iu.EQ.n)range='A'
    END IF
    il_ifail=n
    il_iwork=5*n
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_ifail,ifail,procedureN//'_ifail',ierr)
#else
    ALLOCATE(ifail(il_ifail(1)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_iwork,iwork,procedureN//'_iwork',ierr)
#else
    ALLOCATE(iwork(il_iwork(1)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    !workspace query
    il_work(1)=-1_int_8
    IF(iopt_save.eq.iopt) then
       if(n_save.eq.n)then
          il_work(1)=il_work_save
       end if
    end IF
    IF(il_work(1).EQ.-1_int_8)THEN
       !workspace query
       CALL cpmd_dsyevx(jobz,range,uplo,n,a,n,1.0_real_8,1.0_real_8,il,iu,0.0_real_8,m,w,&
            z,n,dummy_real,INT(il_work(1)),iwork,ifail,info)
       il_work(1)=INT(dummy_real(1))
       il_work_save=il_work(1)
       n_save=n
       iopt_save=iopt
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_work,work,procedureN//'_work',ierr)
#else
    ALLOCATE(work(il_work(1)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    !actual calculation
    CALL cpmd_dsyevx(jobz,range,uplo,n,a,n,1.0_real_8,1.0_real_8,il,iu,-1_real_8,m,w,&
         z,n,work,il_work(1),iwork,ifail,info)
    IF (info.NE.0) CALL stopgm(procedureN,'FAILED TO DIAGONALIZE',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_work,work,procedureN//'_work',ierr)
#else
    DEALLOCATE(work,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_iwork,iwork,procedureN//'_iwork',ierr)
#else
    DEALLOCATE(iwork,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_ifail,ifail,procedureN//'_ifail',ierr)
#else
    DEALLOCATE(ifail,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
  END SUBROUTINE dsyevx_driver
    ! ==--------------------------------------------------------------==
  SUBROUTINE dspevy(iopt,ap,w,z,ldz,n,aux,naux)
    ! ==--------------------------------------------------------------==
    ! == DIAGONALIZATION ROUTINE: FOLLOW THE ESSL CONVENTION          ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iopt, ldz, n
    REAL(real_8)                             :: z(ldz,n), w(n), ap(n*(n+1)/2)
    INTEGER                                  :: naux
    REAL(real_8)                             :: aux(naux)

    INTEGER                                  :: info

#if defined ( ESSL )
    CALL dspev(iopt,ap,w,z,ldz,n,aux,naux)
#else
    info=0
    IF (iopt.EQ.0) THEN
       CALL dspev('N','L',n,ap,w,z,ldz,aux,info)
    ELSEIF (iopt.EQ.1) THEN
       CALL dspev('V','L',n,ap,w,z,ldz,aux,info)
    ELSEIF (iopt.EQ.20) THEN
       CALL dspev('N','U',n,ap,w,z,ldz,aux,info)
    ELSEIF (iopt.EQ.21) THEN
       CALL dspev('V','U',n,ap,w,z,ldz,aux,info)
    ENDIF
    IF (info.NE.0) CALL stopgm('DSPEVY','FAILED TO DIAGONALIZE',&
         __LINE__,__FILE__)
#endif
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dspevy
  ! ==================================================================
#if ! (defined(ESSL))
  ! ==================================================================
  SUBROUTINE zgthr(n,a,b,ind)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: a(*), b(n)
    INTEGER                                  :: ind(n)

    INTEGER                                  :: i
    if(update_first_to_gpu.or.update_second_to_gpu.or.update_result_to_host) write(*,*) 'zgthr'
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target update to(a(1:maxval(ind))) IF(update_first_to_gpu)
    !$omp target update to(b) IF(update_second_to_gpu)
    !$omp target teams distribute parallel do &
#else
    !$omp parallel do &
#endif
    !$omp&private(I)
    DO i=1,n
       b(i)=a(ind(i))
    ENDDO
    !$omp target update from(b) IF(update_result_to_host)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE zgthr
  ! ==================================================================
  SUBROUTINE zgthr_no_omp(n,a,b,ind)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: a(*), b(n)
    INTEGER                                  :: ind(n)

    INTEGER                                  :: i

    DO i=1,n
       b(i)=a(ind(i))
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE zgthr_no_omp
  ! ==================================================================
  SUBROUTINE zsctr_no_omp(n,a,ind,b)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: a(n)
    INTEGER                                  :: ind(n)
    COMPLEX(real_8)                          :: b(*)

    INTEGER                                  :: i

    DO i=1,n
       b(ind(i))=a(i)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE zsctr_no_omp
!! cmb - The following function was defined to overcome a library
!! cmb - problem in IBM-like machines (Hitachi SR11000 inlcuded)
!! cmb - we keep it here just for emergency
#if defined(__SR11KIBM)
  ! ==================================================================
  FUNCTION dzamax(n,a,inc)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: a(*)
    INTEGER                                  :: inc
    REAL(real_8)                             :: dzamax

    INTEGER                                  :: i, ii
    REAL(real_8)                             :: y, ymax

    ymax=-1.0_real_8
    !$omp parallel do private(I,II,Y) reduction(MAX:YMAX)
    DO i=1,n
       ii=(i-1)*inc+1
       y=ABS(a(ii))
       IF (y.GT.ymax) ymax=y
    ENDDO
    dzamax=ymax
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION dzamax
#endif
  ! ==================================================================
#if   ! (defined(__HAS_EXTERNAL_IZAMAX))
  FUNCTION izamax(n,a,inc)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: a(*)
    INTEGER                                  :: inc, izamax

    INTEGER                                  :: i, ii
    REAL(real_8)                             :: xmax, y

    izamax=0
    xmax=-1.0_real_8
#if defined(__SR8000) ||defined(__SR11000) || defined(_vpp_)
    DO i=1,n
       ii=(i-1)*inc+1
       y=a(ii)*CONJG(a(ii))
       IF (y.GT.xmax) THEN
          xmax=y
          izamax=i
       ENDIF
    ENDDO
    xmax=SQRT(xmax)
#else
    DO i=1,n
       ii=(i-1)*inc+1
       y=ABS(a(ii))
       IF (y.GT.xmax) THEN
          xmax=y
          izamax=i
       ENDIF
    ENDDO
#endif
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION izamax
  ! ==================================================================
#endif
#endif
#if ! (defined(__HAS_EXTERNAL_IDAMIN))
  FUNCTION idamin(n,x,inx)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: x(n)
    INTEGER                                  :: inx, idamin

    INTEGER                                  :: i, ii, ij
    REAL(real_8)                             :: xmin

    xmin=x(1)
    ij=1
    DO i=2,n
       ii=(i-1)*inx+1
       IF (xmin.GE.ABS(x(ii))) ij=i
    ENDDO
    idamin=ij
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION idamin
  ! ==================================================================
#endif
  ! =================================================================
  SUBROUTINE zclean_k(a,n,ngw)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n, ngw
    COMPLEX(real_8)                          :: a(2*ngw,n)

    INTEGER                                  :: i

    IF (ngw.GT.0) THEN
       !$omp parallel do private(I) shared(NGW)
       DO i=1,n
          a(ngw+1,i)=CMPLX(0._real_8,0._real_8,kind=real_8)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE zclean_k
  ! ==================================================================
  SUBROUTINE setkwf(ngw,nstate,a)
    ! ==--------------------------------------------------------------==
    ! == A USES INVERSION SYMMETRY.                                   ==
    ! == WE TRANSFORM A IN ARRAY WITHOUT INVERSION SYMMETRY.          ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ngw, nstate
    COMPLEX(real_8)                          :: a(ngw,2*nstate)

    INTEGER                                  :: i, ii, j

    DO i=nstate,2,-1
       CALL dcopy(2*ngw,a(1,i),1,a(1,2*i-1),1)
    ENDDO
    DO i=1,nstate
       ii=2*i
       CALL dcopy(2*ngw,a(1,ii-1),1,a(1,ii),1)
       DO j=1,ngw
          a(j,ii)=CONJG(a(j,ii))
       ENDDO
    ENDDO
    IF (geq0) CALL zclean_k(a,nstate,ngw)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE setkwf
  ! ==================================================================
  SUBROUTINE fskip(nf,nrec)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nf, nrec

    INTEGER                                  :: i

    IF (nrec.LT.0) THEN
       WRITE(6,*) ' FSKIP! NREC=',NREC
       CALL stopgm('FSKIP','BAD NUMBER OF RECORDS (NREC) TO SKIP',&
            __LINE__,__FILE__)
    ENDIF
    DO i=1,nrec
       READ(nf,END=20,err=30)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
20  CONTINUE
    CALL stopgm('FSKIP','END OF FILE',&
         __LINE__,__FILE__)
30  CONTINUE
    CALL stopgm('FSKIP','READ ERROR',&
         __LINE__,__FILE__)
  END SUBROUTINE fskip
  ! ==================================================================
  FUNCTION dgive(a,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: a(n), dgive

! ==--------------------------------------------------------------==

    dgive=a(n)
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION dgive
  ! ==================================================================
  FUNCTION zgive(a,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: a(n), zgive

! ==--------------------------------------------------------------==

    zgive=a(n)
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION zgive

  SUBROUTINE determ(dmat,ld,n,dd)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ld
    COMPLEX(real_8)                          :: dmat(ld,*)
    INTEGER                                  :: n
    COMPLEX(real_8)                          :: dd

    CHARACTER(*), PARAMETER                  :: procedureN = 'determ'

    INTEGER                                  :: i, ierr, info
    INTEGER, ALLOCATABLE                     :: ipiv(:)

    ALLOCATE(ipiv(n),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    CALL zgetrf(n,n,dmat,ld,ipiv,info)
    IF (info.NE.0) CALL stopgm('DETERM','ILLEGAL RESULTS',&
         __LINE__,__FILE__)
    DEALLOCATE(ipiv,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    dd=CMPLX(1._real_8,0._real_8,kind=real_8)
    DO i=1,n
       dd=dd*dmat(i,i)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE determ
  ! ==================================================================
  FUNCTION nxxfun(nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate, nxxfun

    INTEGER                                  :: i, n1, n2, ns

    n1=0
    n2=0
    DO i=0,parai%nproc-1
       ns=parap%nst12(i,2)-parap%nst12(i,1)+1
       n1=MAX(n1,ns)
       n2=MAX(n2,parap%sparm(3,i))
    ENDDO
    nxxfun=n1*n2
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION nxxfun
  ! ==================================================================
  FUNCTION sbes0(x)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: x, sbes0

    REAL(real_8), PARAMETER                  :: small = 1.e-35_real_8

! ==--------------------------------------------------------------==
! ==  CALCULATE SPHERICAL BESSEL FUNCTION OF THE FIRST KIND       ==
! ==--------------------------------------------------------------==

    IF (ABS(x).LT.small) THEN
       sbes0=1.0_real_8 - x*x/6.0_real_8 + x*x*x*x/120.0_real_8
    ELSE
       sbes0=SIN(x)/x
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION sbes0
  ! ==================================================================
  FUNCTION fc4(x,y)
    REAL(real_8)                             :: x, y, fc4

    REAL(real_8), PARAMETER                  :: small = 1.e-10_real_8

    REAL(real_8)                             :: diff, dxy

! ..CALCULATE (X*SIN(Y) - Y*SIN(X))/(X*Y*(Y^2 - X^2))

    diff = y*y - x*x
    IF (ABS(diff).GT.small) THEN
       fc4=(x*SIN(y)-y*SIN(x))/(x*y*(y*y-x*x))
    ELSE
       dxy = x - y
       fc4 = -0.5_real_8*sb1(y) + 0.25_real_8*sb2(y)*dxy&
            -1.0_real_8/120.0_real_8*(7.0_real_8*sb3(y) - 3.0_real_8*sb1(y))*dxy*dxy&
            +1.0_real_8/336.0_real_8*(3.0_real_8*sb4(y) - 4.0_real_8*sb2(y))*dxy*dxy*dxy&
            +1.0_real_8/30240.0_real_8*(-31.0_real_8*sb5(y)+77.0_real_8*sb3(y)&
            -18.0_real_8*sb1(y))*dxy*dxy*dxy*dxy
    ENDIF
    RETURN
  END FUNCTION fc4
  ! ======================================================================C
  FUNCTION sb1(y)
    REAL(real_8)                             :: y, sb1

    REAL(real_8), PARAMETER                  :: small = 1.e-6_real_8

    IF (ABS(y).GT.small) THEN
       sb1=SIN(y)/y**3 - COS(y)/y**2
    ELSE
       sb1=(1._real_8-y**2/10.0_real_8+y**4/280.0_real_8)/3.0_real_8
    ENDIF
    RETURN
  END FUNCTION sb1
  ! ======================================================================C
  FUNCTION sb2(y)
    REAL(real_8)                             :: y, sb2

    REAL(real_8), PARAMETER                  :: small = 1.e-6_real_8

    IF (ABS(y).GT.small) THEN
       sb2=SIN(y)*(3.0_real_8/y**4-1.0_real_8/y**2) - COS(y)*(3.0_real_8/y**3)
    ELSE
       sb2=(1._real_8-y**2/14.0_real_8+y**4/504.0_real_8)*y/15.0_real_8
    ENDIF
    RETURN
  END FUNCTION sb2
  ! ======================================================================C
  FUNCTION sb3(y)
    REAL(real_8)                             :: y, sb3

    REAL(real_8), PARAMETER                  :: small = 1.e-6_real_8

    IF (ABS(y).GT.small) THEN
       sb3=SIN(y)*(15.0_real_8/y**5-6.0_real_8/y**3)&
            + COS(y)*(-15.0_real_8/y**4+1.0_real_8/y**2)
    ELSE
       sb3=(1._real_8-y**2/18.0_real_8)*y**2/105.0_real_8
    ENDIF
    RETURN
  END FUNCTION sb3
  ! ======================================================================C
  FUNCTION sb4(y)
    REAL(real_8)                             :: y, sb4

    REAL(real_8), PARAMETER                  :: small = 1.e-6_real_8

    IF (ABS(y).GT.small) THEN
       sb4=SIN(y)*(105.0_real_8/y**6-45.0_real_8/y**4+1.0_real_8/y**2)&
            + COS(y)*(-105.0_real_8/y**5+10.0_real_8/y**3)
    ELSE
       sb4=(1._real_8-y**2/22.0_real_8)*y**3/945.0_real_8
    ENDIF
    RETURN
  END FUNCTION sb4
  ! ======================================================================C
  FUNCTION sb5(y)
    REAL(real_8)                             :: y, sb5

    REAL(real_8), PARAMETER                  :: small = 1.e-6_real_8

    IF (ABS(y).GT.small) THEN
       sb5=SIN(y)*(945.0_real_8/y**7-420.0_real_8/y**5+15.0_real_8/y**3)&
            + COS(y)*(-945.0_real_8/y**6+105.0_real_8/y**4-1.0_real_8/y**2)
    ELSE
       sb5=y**4/10395.0_real_8
    ENDIF
    RETURN
  END FUNCTION sb5
  ! ==================================================================
  SUBROUTINE rtrans(a,n)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n
    REAL(real_8)                             :: a(n,n)

    INTEGER                                  :: i, j
    REAL(real_8)                             :: aij

! ==--------------------------------------------------------------==

    !$omp parallel do private(I,J,AIJ) shared(N,A) schedule(static,1)
    DO i=1,n
       DO j=i+1,n
          aij=a(i,j)
          a(i,j)=a(j,i)
          a(j,i)=aij
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rtrans
  ! ==================================================================
  SUBROUTINE rmatmov(n,m,a,lda,b,ldb)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n, m, lda
    REAL(real_8)                             :: a(lda,*)
    INTEGER                                  :: ldb
    REAL(real_8)                             :: b(ldb,*)

    INTEGER                                  :: i, j

    !$omp parallel do private(I,J) shared(M,N,A,B)
    DO i=1,m
       DO j=1,n
          b(j,i)=a(j,i)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rmatmov
  ! ==================================================================
  SUBROUTINE icopy(n,a,ia,b,ib)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: n, ia, a(ia*n), ib, b(ib*n)

    INTEGER                                  :: i, iax, ibx

! Variables
! ==--------------------------------------------------------------==

    !$omp parallel do private(I,IAX,IBX)
    DO i=1,n
       iax=(i-1)*ia + 1
       ibx=(i-1)*ib + 1
       b(ibx) = a(iax)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE icopy

  ! ==================================================================

  ! ==================================================================
  SUBROUTINE print_debug_ions(name,array)
    ! ==--------------------------------------------------------------==
    CHARACTER(*)                             :: name
    REAL(real_8)                             :: array(:,:,:)

    INTEGER                                  :: is, ia, ierr
    REAL(real_8), ALLOCATABLE                :: debug(:,:,:)
    CHARACTER(*), PARAMETER                  :: procedureN='print_debug_ions'

    ALLOCATE(debug(3,maxsys%nax,maxsys%nsx), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate dbg_forces',& 
         __LINE__,__FILE__)
    !$omp target enter data map(alloc:debug,array) if(.not.comm_buffers_on_host)
    !$omp target update to(array) if(.not.comm_buffers_on_host)
    CALL mp_sum(array,debug,3*maxsys%nax*maxsys%nsx,parai%allgrp)
    !$omp target update from(debug) if(.NOT.comm_buffers_on_host)
    IF (paral%io_parent) THEN
       WRITE(6,*) "===================================="
       WRITE(6,*) TRIM(ADJUSTL(name))
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             WRITE(6,*) debug(1:3,ia,is),ia,is
          END DO
       END DO
    END IF
    !$omp target exit data map(delete:debug) if(.not.comm_buffers_on_host)
    !$omp target exit data map(release:array) if(.not.comm_buffers_on_host)
    DEALLOCATE(debug,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE print_debug_ions
  ! ==================================================================

END MODULE utils


! ==================================================================
SUBROUTINE dcopy_s(n,dx,incx,dy,incy)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  IMPLICIT NONE
  INTEGER                                    :: n
  REAL(real_4)                               :: dx(*)
  INTEGER                                    :: incx
  REAL(real_8)                               :: dy(*)
  INTEGER                                    :: incy

  INTEGER                                    :: i, ix, iy, m, mp1

  IF (n.LE.0) RETURN
  IF (incx.EQ.1 .AND. incy.EQ.1) THEN
     ! ==--------------------------------------------------------------==
     m = MOD(n,7)
     IF (m.NE.0) THEN
        DO i = 1,m
           dy(i) = dx(i)
        ENDDO
        IF (n.LT.7) RETURN
     ENDIF
     mp1 = m + 1
     !$omp parallel do private(I) shared(N,MP1,DX,DY)
     DO i = mp1,n,7
        dy(i) = dx(i)
        dy(i+1) = dx(i+1)
        dy(i+2) = dx(i+2)
        dy(i+3) = dx(i+3)
        dy(i+4) = dx(i+4)
        dy(i+5) = dx(i+5)
        dy(i+6) = dx(i+6)
     ENDDO
  ELSE
     ix = 1
     iy = 1
     IF (incx.LT.0) ix = (-n+1)*incx + 1
     IF (incy.LT.0) iy = (-n+1)*incy + 1
     DO i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
     ENDDO
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE dcopy_s
! ==================================================================
SUBROUTINE cgthr_z(n,a,b,ind)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  IMPLICIT NONE
  INTEGER                                    :: n
  COMPLEX(real_8)                            :: a(*)
  COMPLEX(real_4)                            :: b(n)
  INTEGER                                    :: ind(n)

  INTEGER                                    :: i, m

! ==--------------------------------------------------------------==
! Variables
! ==--------------------------------------------------------------==

  m = MOD(n,7)
  IF (m.NE.0) THEN
     DO i = 1,m
        b(i) = a(ind(i))
     ENDDO
     IF (n.LT.7) RETURN
  ENDIF
  !$omp parallel do private(I) shared(N,M,A,B,IND)
  DO i = m+1,n,7
     b(i  ) = a(ind(i  ))
     b(i+1) = a(ind(i+1))
     b(i+2) = a(ind(i+2))
     b(i+3) = a(ind(i+3))
     b(i+4) = a(ind(i+4))
     b(i+5) = a(ind(i+5))
     b(i+6) = a(ind(i+6))
  ENDDO
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE cgthr_z

! ==================================================================
SUBROUTINE scopy_d(n,dx,incx,dy,incy)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  IMPLICIT NONE
  INTEGER                                    :: n
  REAL(real_8)                               :: dx(*)
  INTEGER                                    :: incx
  REAL(real_4)                               :: dy(*)
  INTEGER                                    :: incy

  INTEGER                                    :: i, ix, iy, m, mp1

! ==--------------------------------------------------------------==
! Variables
! ==--------------------------------------------------------------==

  IF (n.LE.0) RETURN
  IF (incx.EQ.1 .AND. incy.EQ.1) THEN
     ! ==--------------------------------------------------------------==
     m = MOD(n,7)
     IF (m.NE.0) THEN
        DO i = 1,m
           dy(i) = dx(i)
        ENDDO
        IF (n.LT.7) RETURN
     ENDIF
     mp1 = m + 1
     !$omp parallel do private(I) shared(N,MP1,DX,DY)
     DO i = mp1,n,7
        dy(i) = dx(i)
        dy(i+1) = dx(i+1)
        dy(i+2) = dx(i+2)
        dy(i+3) = dx(i+3)
        dy(i+4) = dx(i+4)
        dy(i+5) = dx(i+5)
        dy(i+6) = dx(i+6)
     ENDDO
  ELSE
     ix = 1
     iy = 1
     IF (incx.LT.0) ix = (-n+1)*incx + 1
     IF (incy.LT.0) iy = (-n+1)*incy + 1
     DO i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
     ENDDO
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE scopy_d
! ==================================================================
SUBROUTINE zsctr_c(n,a,ind,b)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  IMPLICIT NONE
  INTEGER                                    :: n
  COMPLEX(real_4)                            :: a(n)
  INTEGER                                    :: ind(n)
  COMPLEX(real_8)                            :: b(*)

  INTEGER                                    :: i, m

! ==--------------------------------------------------------------==
! Variables
! ==--------------------------------------------------------------==

  m = MOD(n,7)
  IF (m.NE.0) THEN
     DO i = 1,m
        b(ind(i)) = a(i)
     ENDDO
     IF (n.LT.7) RETURN
  ENDIF
  !$omp parallel do private(I) shared(N,M,A,B,IND)
  DO i = m+1,n,7
     b(ind(i  )) = a(i  )
     b(ind(i+1)) = a(i+1)
     b(ind(i+2)) = a(i+2)
     b(ind(i+3)) = a(i+3)
     b(ind(i+4)) = a(i+4)
     b(ind(i+5)) = a(i+5)
     b(ind(i+6)) = a(i+6)
  ENDDO
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE zsctr_c


#if ! (defined(ESSL))
! ==================================================================
SUBROUTINE zgetmo(a,lda,m,n,b,ldb)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  IMPLICIT NONE
  INTEGER                                    :: lda, m, n
  COMPLEX(real_8)                            :: a(lda,n)
  INTEGER                                    :: ldb
  COMPLEX(real_8)                            :: b(ldb,m)

  INTEGER                                    :: i, j

! ==--------------------------------------------------------------==
! Variables
! ==--------------------------------------------------------------==

  !$omp parallel do private(I,J) shared(M,N,A,B)
  DO i=1,m
     DO j=1,n
        b(j,i)=a(i,j)
     ENDDO
  ENDDO
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE zgetmo
#endif

! ==================================================================
SUBROUTINE matmov(n,m,a,lda,b,ldb)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  IMPLICIT NONE
  INTEGER                                    :: n, m, lda
  COMPLEX(real_8)                            :: a(lda,*)
  INTEGER                                    :: ldb
  COMPLEX(real_8)                            :: b(ldb,*)

! ==--------------------------------------------------------------==
! Variables

#if defined(__VECTOR)
  INTEGER :: i,j
  ! ==--------------------------------------------------------------==
  !$omp parallel do private(I,J) shared(M,N,A,B)
  DO i=1,m
     DO j=1,n
        b(j,i)=a(j,i)
     ENDDO
  ENDDO
#elif defined(__HP)
  INTEGER :: i,j
  ! ==--------------------------------------------------------------==
  DO i=1,m
     DO j=1,n
        b(j,i)=a(j,i)
     ENDDO
  ENDDO
#else
  INTEGER :: i,j,jj
  ! ==--------------------------------------------------------------==
  jj=MOD(n,4)
  !$omp parallel do private(I,J) shared(M,N,A,B,JJ)
  DO i=1,m
     DO j=1,jj
        b(j,i)=a(j,i)
     ENDDO
     DO j=1+jj,n,4
        b(j  ,i)  =  a(j  ,i)
        b(j+1,i)  =  a(j+1,i)
        b(j+2,i)  =  a(j+2,i)
        b(j+3,i)  =  a(j+3,i)
     ENDDO
  ENDDO
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE matmov
