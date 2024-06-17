subroutine mp_all2all_int4_r1(data_in, data_out, n, comm )
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n
  type(MPI_COMM), intent(in) ::  comm
#else
  integer, intent(in) :: n, comm
#endif
  integer(int_4), intent(in) , dimension(*) ::  data_in
  integer(int_4), intent(out) , dimension(*) :: data_out
  ! Variables
  integer :: ierr,nnodes
  character(*),parameter :: procedureN='mp_all2all_int4_r1'
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL  
  call mpi_comm_size ( comm, nnodes, ierr )
  cmcal(ipar_aall)=cmcal(ipar_aall)+1.0d0
  cmlen(ipar_aall)=cmlen(ipar_aall)+n*nnodes*mp_int4_in_bytes
  tim1=m_walltime()     
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target data use_device_addr(data_in,data_out) if(.NOT.comm_buffers_on_host)
#endif
  CALL mpi_alltoall(data_in,n,MPI_INTEGER,data_out,n,MPI_INTEGER,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp end target data
#endif
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_aall)=cmtim(ipar_aall)+tim2-tim1
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_int4_r1
subroutine mp_all2all_int8_r1(data_in, data_out, n, comm )
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n
  type(MPI_COMM), intent(in) ::  comm
#else
  integer, intent(in) :: n, comm
#endif
  integer(int_8), intent(in) , dimension(*) ::  data_in
  integer(int_8), intent(out) , dimension(*) :: data_out
  ! Variables
  integer :: ierr,nnodes
  character(*),parameter :: procedureN='mp_all2all_int8_r1'
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL  
  call mpi_comm_size ( comm, nnodes, ierr )
  cmcal(ipar_aall)=cmcal(ipar_aall)+1.0d0
  cmlen(ipar_aall)=cmlen(ipar_aall)+n*nnodes*mp_int8_in_bytes
  tim1=m_walltime()     
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target data use_device_addr(data_in,data_out) if(.NOT.comm_buffers_on_host)
#endif
  CALL mpi_alltoall(data_in,n,MPI_INTEGER8,data_out,n,MPI_INTEGER8,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp end target data
#endif
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_aall)=cmtim(ipar_aall)+tim2-tim1
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_int8_r1
subroutine mp_all2all_real8_r1(data_in, data_out, n, comm )
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n
  type(MPI_COMM), intent(in) ::  comm
#else
  integer, intent(in) :: n, comm
#endif
  real(real_8), intent(in) , dimension(*) ::  data_in
  real(real_8), intent(out) , dimension(*) :: data_out
  ! Variables
  integer :: ierr,nnodes
  character(*),parameter :: procedureN='mp_all2all_real8_r1'
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL  
  call mpi_comm_size ( comm, nnodes, ierr )
  cmcal(ipar_aall)=cmcal(ipar_aall)+1.0d0
  cmlen(ipar_aall)=cmlen(ipar_aall)+n*nnodes*mp_double_in_bytes
  tim1=m_walltime()     
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target data use_device_addr(data_in,data_out) if(.NOT.comm_buffers_on_host)
#endif
  CALL mpi_alltoall(data_in,n,MPI_DOUBLE_PRECISION,data_out,n,MPI_DOUBLE_PRECISION,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp end target data
#endif
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_aall)=cmtim(ipar_aall)+tim2-tim1
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_real8_r1
subroutine mp_all2all_complex4_r1(data_in, data_out, n, comm )
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n
  type(MPI_COMM), intent(in) ::  comm
#else
  integer, intent(in) :: n, comm
#endif
  complex(real_4), intent(in) , dimension(*) ::  data_in
  complex(real_4), intent(out) , dimension(*) :: data_out
  ! Variables
  integer :: ierr,nnodes
  character(*),parameter :: procedureN='mp_all2all_complex4_r1'
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL  
  call mpi_comm_size ( comm, nnodes, ierr )
  cmcal(ipar_aall)=cmcal(ipar_aall)+1.0d0
  cmlen(ipar_aall)=cmlen(ipar_aall)+n*nnodes*mp_complex_in_bytes
  tim1=m_walltime()     
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target data use_device_addr(data_in,data_out) if(.NOT.comm_buffers_on_host)
#endif
  CALL mpi_alltoall(data_in,n,MPI_COMPLEX,data_out,n,MPI_COMPLEX,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp end target data
#endif
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_aall)=cmtim(ipar_aall)+tim2-tim1
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_complex4_r1
subroutine mp_all2all_complex8_r1(data_in, data_out, n, comm )
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n
  type(MPI_COMM), intent(in) ::  comm
#else
  integer, intent(in) :: n, comm
#endif
  complex(real_8), intent(in) , dimension(*) ::  data_in
  complex(real_8), intent(out) , dimension(*) :: data_out
  ! Variables
  integer :: ierr,nnodes
  character(*),parameter :: procedureN='mp_all2all_complex8_r1'
  real(real_8) :: tim1,tim2
#ifdef __PARALLEL  
  call mpi_comm_size ( comm, nnodes, ierr )
  cmcal(ipar_aall)=cmcal(ipar_aall)+1.0d0
  cmlen(ipar_aall)=cmlen(ipar_aall)+n*nnodes*mp_double_complex_in_bytes
  tim1=m_walltime()     
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target data use_device_addr(data_in,data_out) if(.NOT.comm_buffers_on_host)
#endif
  CALL mpi_alltoall(data_in,n,MPI_DOUBLE_COMPLEX,data_out,n,MPI_DOUBLE_COMPLEX,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp end target data
#endif
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_aall)=cmtim(ipar_aall)+tim2-tim1
#else
  data_out(1:n)=data_in(1:n)
#endif
end subroutine mp_all2all_complex8_r1
