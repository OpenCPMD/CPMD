SUBROUTINE MP_MAX_int4_r1(data,n,comm)
  IMPLICIT NONE
  ! Arguments
#ifdef __PARALLEL
  INTEGER, INTENT(IN) :: N
  TYPE(MPI_COMM), INTENT(IN) :: COMM
#else
  INTEGER, INTENT(IN) :: N,COMM
#endif
  integer(int_4), INTENT(INOUT), TARGET , dimension(*) :: data
  ! Variables
  integer(int_4), POINTER __CONTIGUOUS :: para_buff_user(:), data_ptr(:)
  TYPE(c_ptr) :: loc_x

  CHARACTER(*),PARAMETER :: procedureN='mp_max_int4_r1'
  INTEGER :: ierr, i
  REAL(REAL_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_int4_in_bytes
  tim1=m_walltime()
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#else
     !$omp target update from(data(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
     CALL mpi_allreduce(MPI_IN_PLACE,data,n,MPI_INTEGER,&
          MPI_MAX,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
     !$omp end target data
#else
     !$omp target update to(data(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  ELSE
#if !defined _USE_SCRATCHLIBRARY
     IF(.NOT.ALLOCATED(para_buff))THEN
        ALLOCATE(para_buff(il_para_buff(1)),stat=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate para_buff', &
             __LINE__,__FILE__)
     END IF
#endif
     IF(n*mp_int4_in_bytes.GT.buff_size_in_bytes*il_para_buff(1))THEN
	il_para_buff(1)=CEILING(REAL(n*mp_int4_in_bytes,real_8)/&
             buff_size_in_bytes)
#if !defined _USE_SCRATCHLIBRARY
        DEALLOCATE(para_buff,ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
             __LINE__,__FILE__)
        ALLOCATE(para_buff(il_para_buff(1)),ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate para_buff', &
             __LINE__,__FILE__)
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL request_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate para_buff', &
          __LINE__,__FILE__)
#endif
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#else
     !$omp target update from(data(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
     CALL mpi_allreduce(data,para_buff,n,MPI_INTEGER,&
          MPI_MAX,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
     !$omp end target data
#else
     !$omp target update to(para_buff(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     loc_x = C_LOC(para_buff)
     CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
     loc_x = C_LOC(data)
     CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     IF(.NOT.comm_buffers_on_host)THEN
        !$omp target teams distribute parallel do simd
#else
        !$omp parallel do simd
#endif
        DO i=1,n
           data_ptr(i)=para_buff_user(i)
        END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     ELSE
        !$omp parallel do simd
        do i=1,n
           data_ptr(i)=para_buff_user(i)
        END DO
     END IF
#endif
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
END SUBROUTINE MP_MAX_int4_r1
SUBROUTINE MP_MAX_real8_r1(data,n,comm)
  IMPLICIT NONE
  ! Arguments
#ifdef __PARALLEL
  INTEGER, INTENT(IN) :: N
  TYPE(MPI_COMM), INTENT(IN) :: COMM
#else
  INTEGER, INTENT(IN) :: N,COMM
#endif
  real(real_8), INTENT(INOUT), TARGET , dimension(*) :: data
  ! Variables
  real(real_8), POINTER __CONTIGUOUS :: para_buff_user(:), data_ptr(:)
  TYPE(c_ptr) :: loc_x

  CHARACTER(*),PARAMETER :: procedureN='mp_max_real8_r1'
  INTEGER :: ierr, i
  REAL(REAL_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_in_bytes
  tim1=m_walltime()
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#else
     !$omp target update from(data(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
     CALL mpi_allreduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_PRECISION,&
          MPI_MAX,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
     !$omp end target data
#else
     !$omp target update to(data(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  ELSE
#if !defined _USE_SCRATCHLIBRARY
     IF(.NOT.ALLOCATED(para_buff))THEN
        ALLOCATE(para_buff(il_para_buff(1)),stat=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate para_buff', &
             __LINE__,__FILE__)
     END IF
#endif
     IF(n*mp_double_in_bytes.GT.buff_size_in_bytes*il_para_buff(1))THEN
	il_para_buff(1)=CEILING(REAL(n*mp_double_in_bytes,real_8)/&
             buff_size_in_bytes)
#if !defined _USE_SCRATCHLIBRARY
        DEALLOCATE(para_buff,ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
             __LINE__,__FILE__)
        ALLOCATE(para_buff(il_para_buff(1)),ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate para_buff', &
             __LINE__,__FILE__)
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL request_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate para_buff', &
          __LINE__,__FILE__)
#endif
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#else
     !$omp target update from(data(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
     CALL mpi_allreduce(data,para_buff,n,MPI_DOUBLE_PRECISION,&
          MPI_MAX,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
     !$omp end target data
#else
     !$omp target update to(para_buff(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     loc_x = C_LOC(para_buff)
     CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
     loc_x = C_LOC(data)
     CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     IF(.NOT.comm_buffers_on_host)THEN
        !$omp target teams distribute parallel do simd
#else
        !$omp parallel do simd
#endif
        DO i=1,n
           data_ptr(i)=para_buff_user(i)
        END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     ELSE
        !$omp parallel do simd
        do i=1,n
           data_ptr(i)=para_buff_user(i)
        END DO
     END IF
#endif
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
END SUBROUTINE MP_MAX_real8_r1
SUBROUTINE MP_MAX_complex8_r1(data,n,comm)
  IMPLICIT NONE
  ! Arguments
#ifdef __PARALLEL
  INTEGER, INTENT(IN) :: N
  TYPE(MPI_COMM), INTENT(IN) :: COMM
#else
  INTEGER, INTENT(IN) :: N,COMM
#endif
  complex(real_8), INTENT(INOUT), TARGET , dimension(*) :: data
  ! Variables
  complex(real_8), POINTER __CONTIGUOUS :: para_buff_user(:), data_ptr(:)
  TYPE(c_ptr) :: loc_x

  CHARACTER(*),PARAMETER :: procedureN='mp_max_complex8_r1'
  INTEGER :: ierr, i
  REAL(REAL_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_complex_in_bytes
  tim1=m_walltime()
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#else
     !$omp target update from(data(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
     CALL mpi_allreduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_COMPLEX,&
          MPI_MAX,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
     !$omp end target data
#else
     !$omp target update to(data(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  ELSE
#if !defined _USE_SCRATCHLIBRARY
     IF(.NOT.ALLOCATED(para_buff))THEN
        ALLOCATE(para_buff(il_para_buff(1)),stat=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate para_buff', &
             __LINE__,__FILE__)
     END IF
#endif
     IF(n*mp_double_complex_in_bytes.GT.buff_size_in_bytes*il_para_buff(1))THEN
	il_para_buff(1)=CEILING(REAL(n*mp_double_complex_in_bytes,real_8)/&
             buff_size_in_bytes)
#if !defined _USE_SCRATCHLIBRARY
        DEALLOCATE(para_buff,ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
             __LINE__,__FILE__)
        ALLOCATE(para_buff(il_para_buff(1)),ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate para_buff', &
             __LINE__,__FILE__)
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL request_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate para_buff', &
          __LINE__,__FILE__)
#endif
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#else
     !$omp target update from(data(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
     CALL mpi_allreduce(data,para_buff,n,MPI_DOUBLE_COMPLEX,&
          MPI_MAX,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
     !$omp end target data
#else
     !$omp target update to(para_buff(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     loc_x = C_LOC(para_buff)
     CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
     loc_x = C_LOC(data)
     CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     IF(.NOT.comm_buffers_on_host)THEN
        !$omp target teams distribute parallel do simd
#else
        !$omp parallel do simd
#endif
        DO i=1,n
           data_ptr(i)=para_buff_user(i)
        END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     ELSE
        !$omp parallel do simd
        do i=1,n
           data_ptr(i)=para_buff_user(i)
        END DO
     END IF
#endif
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
END SUBROUTINE MP_MAX_complex8_r1
