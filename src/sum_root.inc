subroutine mp_sum_root_int4_r0(data_in,data_out, root, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  root, comm
#endif
#if defined(_HAS_GPU_AWARE_MPI) || !defined(_HAS_OMP_TARGET_OFFLOAD)
  integer(int_4), intent(in)  ::  data_in
#else
  integer(int_4), intent(inout)  ::  data_in
#endif
  integer(int_4), intent(inout)  :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_int4_r0'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+1*mp_int4_in_bytes
  tim1=m_walltime()      
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp target data use_device_addr(data_in,data_out) if(.NOT.comm_buffers_on_host)
#else
  !$omp target update from(data_in) if(.not.comm_buffers_on_host)
#endif
#endif
  call mpi_reduce(data_in,data_out,1,MPI_INTEGER,&
       MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp end target data
#else
  !$omp target update to(data_out) if(.not.comm_buffers_on_host)
#endif
#endif
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out =data_in 
#endif
end subroutine mp_sum_root_int4_r0
subroutine mp_sum_root_int4_r1(data_in,data_out,n, root, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, root, comm
#endif
#if defined(_HAS_GPU_AWARE_MPI) || !defined(_HAS_OMP_TARGET_OFFLOAD)
  integer(int_4), intent(in) , dimension(*) ::  data_in
#else
  integer(int_4), intent(inout) , dimension(*) ::  data_in
#endif
  integer(int_4), intent(inout) , dimension(*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_int4_r1'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_int4_in_bytes
  tim1=m_walltime()      
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp target data use_device_addr(data_in,data_out) if(.NOT.comm_buffers_on_host)
#else
  !$omp target update from(data_in(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
  call mpi_reduce(data_in,data_out,n,MPI_INTEGER,&
       MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp end target data
#else
  !$omp target update to(data_out(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out (1:n)=data_in (1:n)
#endif
end subroutine mp_sum_root_int4_r1
subroutine mp_sum_root_int4_r2(data_in,data_out,n, root, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, root, comm
#endif
#if defined(_HAS_GPU_AWARE_MPI) || !defined(_HAS_OMP_TARGET_OFFLOAD)
  integer(int_4), intent(in) , dimension(1,*) ::  data_in
#else
  integer(int_4), intent(inout) , dimension(1,*) ::  data_in
#endif
  integer(int_4), intent(inout) , dimension(1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_int4_r2'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_int4_in_bytes
  tim1=m_walltime()      
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp target data use_device_addr(data_in,data_out) if(.NOT.comm_buffers_on_host)
#else
  !$omp target update from(data_in(1,1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
  call mpi_reduce(data_in,data_out,n,MPI_INTEGER,&
       MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp end target data
#else
  !$omp target update to(data_out(1,1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out (1,1:n)=data_in (1,1:n)
#endif
end subroutine mp_sum_root_int4_r2
subroutine mp_sum_root_real8_r0(data_in,data_out, root, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  root, comm
#endif
#if defined(_HAS_GPU_AWARE_MPI) || !defined(_HAS_OMP_TARGET_OFFLOAD)
  real(real_8), intent(in)  ::  data_in
#else
  real(real_8), intent(inout)  ::  data_in
#endif
  real(real_8), intent(inout)  :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_real8_r0'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+1*mp_double_in_bytes
  tim1=m_walltime()      
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp target data use_device_addr(data_in,data_out) if(.NOT.comm_buffers_on_host)
#else
  !$omp target update from(data_in) if(.not.comm_buffers_on_host)
#endif
#endif
  call mpi_reduce(data_in,data_out,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp end target data
#else
  !$omp target update to(data_out) if(.not.comm_buffers_on_host)
#endif
#endif
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out =data_in 
#endif
end subroutine mp_sum_root_real8_r0
subroutine mp_sum_root_real8_r1(data_in,data_out,n, root, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, root, comm
#endif
#if defined(_HAS_GPU_AWARE_MPI) || !defined(_HAS_OMP_TARGET_OFFLOAD)
  real(real_8), intent(in) , dimension(*) ::  data_in
#else
  real(real_8), intent(inout) , dimension(*) ::  data_in
#endif
  real(real_8), intent(inout) , dimension(*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_real8_r1'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_in_bytes
  tim1=m_walltime()      
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp target data use_device_addr(data_in,data_out) if(.NOT.comm_buffers_on_host)
#else
  !$omp target update from(data_in(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
  call mpi_reduce(data_in,data_out,n,MPI_DOUBLE_PRECISION,&
       MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp end target data
#else
  !$omp target update to(data_out(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out (1:n)=data_in (1:n)
#endif
end subroutine mp_sum_root_real8_r1
subroutine mp_sum_root_real8_r2(data_in,data_out,n, root, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, root, comm
#endif
#if defined(_HAS_GPU_AWARE_MPI) || !defined(_HAS_OMP_TARGET_OFFLOAD)
  real(real_8), intent(in) , dimension(1,*) ::  data_in
#else
  real(real_8), intent(inout) , dimension(1,*) ::  data_in
#endif
  real(real_8), intent(inout) , dimension(1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_real8_r2'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_in_bytes
  tim1=m_walltime()      
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp target data use_device_addr(data_in,data_out) if(.NOT.comm_buffers_on_host)
#else
  !$omp target update from(data_in(1,1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
  call mpi_reduce(data_in,data_out,n,MPI_DOUBLE_PRECISION,&
       MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp end target data
#else
  !$omp target update to(data_out(1,1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out (1,1:n)=data_in (1,1:n)
#endif
end subroutine mp_sum_root_real8_r2
subroutine mp_sum_root_complex8_r0(data_in,data_out, root, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  root, comm
#endif
#if defined(_HAS_GPU_AWARE_MPI) || !defined(_HAS_OMP_TARGET_OFFLOAD)
  complex(real_8), intent(in)  ::  data_in
#else
  complex(real_8), intent(inout)  ::  data_in
#endif
  complex(real_8), intent(inout)  :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_complex8_r0'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+1*mp_double_complex_in_bytes
  tim1=m_walltime()      
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp target data use_device_addr(data_in,data_out) if(.NOT.comm_buffers_on_host)
#else
  !$omp target update from(data_in) if(.not.comm_buffers_on_host)
#endif
#endif
  call mpi_reduce(data_in,data_out,1,MPI_DOUBLE_COMPLEX,&
       MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp end target data
#else
  !$omp target update to(data_out) if(.not.comm_buffers_on_host)
#endif
#endif
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out =data_in 
#endif
end subroutine mp_sum_root_complex8_r0
subroutine mp_sum_root_complex8_r1(data_in,data_out,n, root, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, root, comm
#endif
#if defined(_HAS_GPU_AWARE_MPI) || !defined(_HAS_OMP_TARGET_OFFLOAD)
  complex(real_8), intent(in) , dimension(*) ::  data_in
#else
  complex(real_8), intent(inout) , dimension(*) ::  data_in
#endif
  complex(real_8), intent(inout) , dimension(*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_complex8_r1'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_complex_in_bytes
  tim1=m_walltime()      
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp target data use_device_addr(data_in,data_out) if(.NOT.comm_buffers_on_host)
#else
  !$omp target update from(data_in(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
  call mpi_reduce(data_in,data_out,n,MPI_DOUBLE_COMPLEX,&
       MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp end target data
#else
  !$omp target update to(data_out(1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out (1:n)=data_in (1:n)
#endif
end subroutine mp_sum_root_complex8_r1
subroutine mp_sum_root_complex8_r2(data_in,data_out,n, root, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, root, comm
#endif
#if defined(_HAS_GPU_AWARE_MPI) || !defined(_HAS_OMP_TARGET_OFFLOAD)
  complex(real_8), intent(in) , dimension(1,*) ::  data_in
#else
  complex(real_8), intent(inout) , dimension(1,*) ::  data_in
#endif
  complex(real_8), intent(inout) , dimension(1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_root_complex8_r2'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: tim1,tim2
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+n*mp_double_complex_in_bytes
  tim1=m_walltime()      
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp target data use_device_addr(data_in,data_out) if(.NOT.comm_buffers_on_host)
#else
  !$omp target update from(data_in(1,1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
  call mpi_reduce(data_in,data_out,n,MPI_DOUBLE_COMPLEX,&
       MPI_SUM,root,comm,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
#if defined(_HAS_GPU_AWARE_MPI)
  !$omp end target data
#else
  !$omp target update to(data_out(1,1:n)) if(.not.comm_buffers_on_host)
#endif
#endif
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out (1,1:n)=data_in (1,1:n)
#endif
end subroutine mp_sum_root_complex8_r2
