subroutine mp_sum_root_in_place_int4_r1(data,n,  root, comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n,  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n,  root, comm
#endif
  integer(int_4), intent(inout), target , dimension(*) :: data
  ! Variables
  integer(int_4), POINTER __CONTIGUOUS :: para_buff_user(:),data_ptr(:)
  TYPE(c_ptr) :: loc_x
  
  character(*),parameter :: procedureN='mp_sum_in_place_int4_r1'
  integer :: ierr,me,i
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_int4_in_bytes
  tim1=m_walltime()      
  CALL mpi_comm_rank(comm,me, ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#endif
     IF (root.EQ.me) THEN
        CALL mpi_reduce(MPI_IN_PLACE,data,n,MPI_INTEGER,&
             MPI_SUM,root,comm,ierr)
     ELSE
        CALL mpi_reduce(data,data,n,MPI_INTEGER,&
             MPI_SUM,root,comm,ierr)
     END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
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
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#endif
     CALL mpi_reduce(data,para_buff,n,MPI_INTEGER,&
          MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     IF(root.EQ.me) THEN
        loc_x = C_LOC(para_buff)
        CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
        loc_x = C_LOC(data)
        CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        if(.NOT.comm_buffers_on_host)then
           !$omp target teams distribute parallel do simd
#else
           !$omp parallel do simd
#endif
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        else
           !$omp parallel do simd
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
        end if
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_root_in_place_int4_r1
subroutine mp_sum_root_in_place_int4_r2(data,n,  root, comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n,  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n,  root, comm
#endif
  integer(int_4), intent(inout), target , dimension(1,*) :: data
  ! Variables
  integer(int_4), POINTER __CONTIGUOUS :: para_buff_user(:),data_ptr(:)
  TYPE(c_ptr) :: loc_x
  
  character(*),parameter :: procedureN='mp_sum_in_place_int4_r2'
  integer :: ierr,me,i
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_int4_in_bytes
  tim1=m_walltime()      
  CALL mpi_comm_rank(comm,me, ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#endif
     IF (root.EQ.me) THEN
        CALL mpi_reduce(MPI_IN_PLACE,data,n,MPI_INTEGER,&
             MPI_SUM,root,comm,ierr)
     ELSE
        CALL mpi_reduce(data,data,n,MPI_INTEGER,&
             MPI_SUM,root,comm,ierr)
     END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
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
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#endif
     CALL mpi_reduce(data,para_buff,n,MPI_INTEGER,&
          MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     IF(root.EQ.me) THEN
        loc_x = C_LOC(para_buff)
        CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
        loc_x = C_LOC(data)
        CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        if(.NOT.comm_buffers_on_host)then
           !$omp target teams distribute parallel do simd
#else
           !$omp parallel do simd
#endif
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        else
           !$omp parallel do simd
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
        end if
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_root_in_place_int4_r2
subroutine mp_sum_root_in_place_int8_r1(data,n, root, comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, root, comm
#endif
  integer(int_8), intent(inout), target , dimension(*) :: data
  ! Variables
  integer(int_8), POINTER __CONTIGUOUS :: para_buff_user(:),data_ptr(:)
  TYPE(c_ptr) :: loc_x
  
  character(*),parameter :: procedureN='mp_sum_in_place_int8_r1'
  integer :: ierr,me,i
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_int8_in_bytes
  tim1=m_walltime()      
  CALL mpi_comm_rank(comm,me, ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#endif
     IF (root.EQ.me) THEN
        CALL mpi_reduce(MPI_IN_PLACE,data,n,MPI_INTEGER8,&
             MPI_SUM,root,comm,ierr)
     ELSE
        CALL mpi_reduce(data,data,n,MPI_INTEGER8,&
             MPI_SUM,root,comm,ierr)
     END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
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
     IF(n*mp_int8_in_bytes.GT.buff_size_in_bytes*il_para_buff(1))THEN
        il_para_buff(1)=CEILING(REAL(n*mp_int8_in_bytes,real_8)/&
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
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#endif
     CALL mpi_reduce(data,para_buff,n,MPI_INTEGER8,&
          MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     IF(root.EQ.me) THEN
        loc_x = C_LOC(para_buff)
        CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
        loc_x = C_LOC(data)
        CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        if(.NOT.comm_buffers_on_host)then
           !$omp target teams distribute parallel do simd
#else
           !$omp parallel do simd
#endif
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        else
           !$omp parallel do simd
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
        end if
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_root_in_place_int8_r1
subroutine mp_sum_root_in_place_real8_r1(data,n,  root, comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n,  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n,  root, comm
#endif
  real(real_8), intent(inout), target , dimension(*) :: data
  ! Variables
  real(real_8), POINTER __CONTIGUOUS :: para_buff_user(:),data_ptr(:)
  TYPE(c_ptr) :: loc_x
  
  character(*),parameter :: procedureN='mp_sum_in_place_real8_r1'
  integer :: ierr,me,i
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_in_bytes
  tim1=m_walltime()      
  CALL mpi_comm_rank(comm,me, ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#endif
     IF (root.EQ.me) THEN
        CALL mpi_reduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_PRECISION,&
             MPI_SUM,root,comm,ierr)
     ELSE
        CALL mpi_reduce(data,data,n,MPI_DOUBLE_PRECISION,&
             MPI_SUM,root,comm,ierr)
     END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
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
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#endif
     CALL mpi_reduce(data,para_buff,n,MPI_DOUBLE_PRECISION,&
          MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     IF(root.EQ.me) THEN
        loc_x = C_LOC(para_buff)
        CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
        loc_x = C_LOC(data)
        CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        if(.NOT.comm_buffers_on_host)then
           !$omp target teams distribute parallel do simd
#else
           !$omp parallel do simd
#endif
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        else
           !$omp parallel do simd
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
        end if
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_root_in_place_real8_r1
subroutine mp_sum_root_in_place_real8_r2(data,n,  root, comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n,  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n,  root, comm
#endif
  real(real_8), intent(inout), target , dimension(1,*) :: data
  ! Variables
  real(real_8), POINTER __CONTIGUOUS :: para_buff_user(:),data_ptr(:)
  TYPE(c_ptr) :: loc_x
  
  character(*),parameter :: procedureN='mp_sum_in_place_real8_r2'
  integer :: ierr,me,i
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_in_bytes
  tim1=m_walltime()      
  CALL mpi_comm_rank(comm,me, ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#endif
     IF (root.EQ.me) THEN
        CALL mpi_reduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_PRECISION,&
             MPI_SUM,root,comm,ierr)
     ELSE
        CALL mpi_reduce(data,data,n,MPI_DOUBLE_PRECISION,&
             MPI_SUM,root,comm,ierr)
     END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
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
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#endif
     CALL mpi_reduce(data,para_buff,n,MPI_DOUBLE_PRECISION,&
          MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     IF(root.EQ.me) THEN
        loc_x = C_LOC(para_buff)
        CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
        loc_x = C_LOC(data)
        CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        if(.NOT.comm_buffers_on_host)then
           !$omp target teams distribute parallel do simd
#else
           !$omp parallel do simd
#endif
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        else
           !$omp parallel do simd
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
        end if
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_root_in_place_real8_r2
subroutine mp_sum_root_in_place_real8_r3(data,n,  root, comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n,  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n,  root, comm
#endif
  real(real_8), intent(inout), target , dimension(1,1,*) :: data
  ! Variables
  real(real_8), POINTER __CONTIGUOUS :: para_buff_user(:),data_ptr(:)
  TYPE(c_ptr) :: loc_x
  
  character(*),parameter :: procedureN='mp_sum_in_place_real8_r3'
  integer :: ierr,me,i
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_in_bytes
  tim1=m_walltime()      
  CALL mpi_comm_rank(comm,me, ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#endif
     IF (root.EQ.me) THEN
        CALL mpi_reduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_PRECISION,&
             MPI_SUM,root,comm,ierr)
     ELSE
        CALL mpi_reduce(data,data,n,MPI_DOUBLE_PRECISION,&
             MPI_SUM,root,comm,ierr)
     END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
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
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#endif
     CALL mpi_reduce(data,para_buff,n,MPI_DOUBLE_PRECISION,&
          MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     IF(root.EQ.me) THEN
        loc_x = C_LOC(para_buff)
        CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
        loc_x = C_LOC(data)
        CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        if(.NOT.comm_buffers_on_host)then
           !$omp target teams distribute parallel do simd
#else
           !$omp parallel do simd
#endif
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        else
           !$omp parallel do simd
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
        end if
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_root_in_place_real8_r3
subroutine mp_sum_root_in_place_real8_r4(data,n,  root, comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n,  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n,  root, comm
#endif
  real(real_8), intent(inout), target , dimension(1,1,1,*) :: data
  ! Variables
  real(real_8), POINTER __CONTIGUOUS :: para_buff_user(:),data_ptr(:)
  TYPE(c_ptr) :: loc_x
  
  character(*),parameter :: procedureN='mp_sum_in_place_real8_r4'
  integer :: ierr,me,i
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_in_bytes
  tim1=m_walltime()      
  CALL mpi_comm_rank(comm,me, ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#endif
     IF (root.EQ.me) THEN
        CALL mpi_reduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_PRECISION,&
             MPI_SUM,root,comm,ierr)
     ELSE
        CALL mpi_reduce(data,data,n,MPI_DOUBLE_PRECISION,&
             MPI_SUM,root,comm,ierr)
     END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
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
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#endif
     CALL mpi_reduce(data,para_buff,n,MPI_DOUBLE_PRECISION,&
          MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     IF(root.EQ.me) THEN
        loc_x = C_LOC(para_buff)
        CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
        loc_x = C_LOC(data)
        CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        if(.NOT.comm_buffers_on_host)then
           !$omp target teams distribute parallel do simd
#else
           !$omp parallel do simd
#endif
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        else
           !$omp parallel do simd
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
        end if
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_root_in_place_real8_r4
subroutine mp_sum_root_in_place_real8_r5(data,n,  root, comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n,  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n,  root, comm
#endif
  real(real_8), intent(inout), target , dimension(1,1,1,1,*) :: data
  ! Variables
  real(real_8), POINTER __CONTIGUOUS :: para_buff_user(:),data_ptr(:)
  TYPE(c_ptr) :: loc_x
  
  character(*),parameter :: procedureN='mp_sum_in_place_real8_r5'
  integer :: ierr,me,i
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_in_bytes
  tim1=m_walltime()      
  CALL mpi_comm_rank(comm,me, ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#endif
     IF (root.EQ.me) THEN
        CALL mpi_reduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_PRECISION,&
             MPI_SUM,root,comm,ierr)
     ELSE
        CALL mpi_reduce(data,data,n,MPI_DOUBLE_PRECISION,&
             MPI_SUM,root,comm,ierr)
     END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
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
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#endif
     CALL mpi_reduce(data,para_buff,n,MPI_DOUBLE_PRECISION,&
          MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     IF(root.EQ.me) THEN
        loc_x = C_LOC(para_buff)
        CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
        loc_x = C_LOC(data)
        CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        if(.NOT.comm_buffers_on_host)then
           !$omp target teams distribute parallel do simd
#else
           !$omp parallel do simd
#endif
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        else
           !$omp parallel do simd
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
        end if
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_root_in_place_real8_r5
subroutine mp_sum_root_in_place_complex8_r1(data,n,  root, comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n,  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n,  root, comm
#endif
  complex(real_8), intent(inout), target , dimension(*) :: data
  ! Variables
  complex(real_8), POINTER __CONTIGUOUS :: para_buff_user(:),data_ptr(:)
  TYPE(c_ptr) :: loc_x
  
  character(*),parameter :: procedureN='mp_sum_in_place_complex8_r1'
  integer :: ierr,me,i
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_complex_in_bytes
  tim1=m_walltime()      
  CALL mpi_comm_rank(comm,me, ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#endif
     IF (root.EQ.me) THEN
        CALL mpi_reduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_COMPLEX,&
             MPI_SUM,root,comm,ierr)
     ELSE
        CALL mpi_reduce(data,data,n,MPI_DOUBLE_COMPLEX,&
             MPI_SUM,root,comm,ierr)
     END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
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
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#endif
     CALL mpi_reduce(data,para_buff,n,MPI_DOUBLE_COMPLEX,&
          MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     IF(root.EQ.me) THEN
        loc_x = C_LOC(para_buff)
        CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
        loc_x = C_LOC(data)
        CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        if(.NOT.comm_buffers_on_host)then
           !$omp target teams distribute parallel do simd
#else
           !$omp parallel do simd
#endif
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        else
           !$omp parallel do simd
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
        end if
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_root_in_place_complex8_r1
subroutine mp_sum_root_in_place_complex8_r2(data,n,  root, comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n,  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n,  root, comm
#endif
  complex(real_8), intent(inout), target , dimension(1,*) :: data
  ! Variables
  complex(real_8), POINTER __CONTIGUOUS :: para_buff_user(:),data_ptr(:)
  TYPE(c_ptr) :: loc_x
  
  character(*),parameter :: procedureN='mp_sum_in_place_complex8_r2'
  integer :: ierr,me,i
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_complex_in_bytes
  tim1=m_walltime()      
  CALL mpi_comm_rank(comm,me, ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#endif
     IF (root.EQ.me) THEN
        CALL mpi_reduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_COMPLEX,&
             MPI_SUM,root,comm,ierr)
     ELSE
        CALL mpi_reduce(data,data,n,MPI_DOUBLE_COMPLEX,&
             MPI_SUM,root,comm,ierr)
     END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
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
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#endif
     CALL mpi_reduce(data,para_buff,n,MPI_DOUBLE_COMPLEX,&
          MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     IF(root.EQ.me) THEN
        loc_x = C_LOC(para_buff)
        CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
        loc_x = C_LOC(data)
        CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        if(.NOT.comm_buffers_on_host)then
           !$omp target teams distribute parallel do simd
#else
           !$omp parallel do simd
#endif
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        else
           !$omp parallel do simd
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
        end if
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_root_in_place_complex8_r2
subroutine mp_sum_root_in_place_complex8_r3(data,n,  root, comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n,  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n,  root, comm
#endif
  complex(real_8), intent(inout), target , dimension(1,1,*) :: data
  ! Variables
  complex(real_8), POINTER __CONTIGUOUS :: para_buff_user(:),data_ptr(:)
  TYPE(c_ptr) :: loc_x
  
  character(*),parameter :: procedureN='mp_sum_in_place_complex8_r3'
  integer :: ierr,me,i
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_complex_in_bytes
  tim1=m_walltime()      
  CALL mpi_comm_rank(comm,me, ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#endif
     IF (root.EQ.me) THEN
        CALL mpi_reduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_COMPLEX,&
             MPI_SUM,root,comm,ierr)
     ELSE
        CALL mpi_reduce(data,data,n,MPI_DOUBLE_COMPLEX,&
             MPI_SUM,root,comm,ierr)
     END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
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
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#endif
     CALL mpi_reduce(data,para_buff,n,MPI_DOUBLE_COMPLEX,&
          MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     IF(root.EQ.me) THEN
        loc_x = C_LOC(para_buff)
        CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
        loc_x = C_LOC(data)
        CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        if(.NOT.comm_buffers_on_host)then
           !$omp target teams distribute parallel do simd
#else
           !$omp parallel do simd
#endif
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        else
           !$omp parallel do simd
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
        end if
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_root_in_place_complex8_r3
subroutine mp_sum_root_in_place_complex8_r4(data,n,  root, comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n,  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n,  root, comm
#endif
  complex(real_8), intent(inout), target , dimension(1,1,1,*) :: data
  ! Variables
  complex(real_8), POINTER __CONTIGUOUS :: para_buff_user(:),data_ptr(:)
  TYPE(c_ptr) :: loc_x
  
  character(*),parameter :: procedureN='mp_sum_in_place_complex8_r4'
  integer :: ierr,me,i
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_complex_in_bytes
  tim1=m_walltime()      
  CALL mpi_comm_rank(comm,me, ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  
  IF(para_use_mpi_in_place)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data) if(.NOT.comm_buffers_on_host)
#endif
     IF (root.EQ.me) THEN
        CALL mpi_reduce(MPI_IN_PLACE,data,n,MPI_DOUBLE_COMPLEX,&
             MPI_SUM,root,comm,ierr)
     ELSE
        CALL mpi_reduce(data,data,n,MPI_DOUBLE_COMPLEX,&
             MPI_SUM,root,comm,ierr)
     END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
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
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp target data use_device_addr(data,para_buff) if(.NOT.comm_buffers_on_host)
#endif
     CALL mpi_reduce(data,para_buff,n,MPI_DOUBLE_COMPLEX,&
          MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
     !$omp end target data
#endif
     CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
     IF(root.EQ.me) THEN
        loc_x = C_LOC(para_buff)
        CALL C_F_POINTER(loc_x, para_buff_user, (/n/))
        loc_x = C_LOC(data)
        CALL C_F_POINTER(loc_x, data_ptr, (/n/))
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        if(.NOT.comm_buffers_on_host)then
           !$omp target teams distribute parallel do simd
#else
           !$omp parallel do simd
#endif
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
#if defined(_HAS_OMP_TARGET_OFFLOAD)
        else
           !$omp parallel do simd
           do i=1,n
              data_ptr(i)=para_buff_user(i)
           end do
        end if
#endif
     END IF
#ifdef _USE_SCRATCHLIBRARY
     CALL free_scratch(il_para_buff,para_buff,procedureN//'_para_buff',ierr)
     IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate para_buff', &
          __LINE__,__FILE__)
#endif
  END IF
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#endif
end subroutine mp_sum_root_in_place_complex8_r4
