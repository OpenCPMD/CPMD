subroutine mp_min_int4_r0(data,comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: comm
#endif
  integer(int_4), intent(inout) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_min_int4_r0'
#ifdef __PARALLEL
  integer :: ierr
  integer(int_4) :: stack_buff
  call mpi_allreduce(data,stack_buff,1,MPI_INTEGER,&
         MPI_MIN,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  data=stack_buff
#endif
end subroutine mp_min_int4_r0
subroutine mp_min_int8_r0(data,comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: comm
#endif
  integer(int_8), intent(inout) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_min_int8_r0'
#ifdef __PARALLEL
  integer :: ierr
  integer(int_8) :: stack_buff
  call mpi_allreduce(data,stack_buff,1,MPI_INTEGER8,&
         MPI_MIN,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  data=stack_buff
#endif
end subroutine mp_min_int8_r0
subroutine mp_min_real8_r0(data,comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: comm
#endif
  real(real_8), intent(inout) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_min_real8_r0'
#ifdef __PARALLEL
  integer :: ierr
  real(real_8) :: stack_buff
  call mpi_allreduce(data,stack_buff,1,MPI_DOUBLE_PRECISION,&
         MPI_MIN,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  data=stack_buff
#endif
end subroutine mp_min_real8_r0
subroutine mp_min_complex8_r0(data,comm)
  implicit none
  ! Arguments
#ifdef __PARALLEL
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: comm
#endif
  complex(real_8), intent(inout) :: data
  ! Variables
  character(*),parameter :: procedureN='mp_min_complex8_r0'
#ifdef __PARALLEL
  integer :: ierr
  complex(real_8) :: stack_buff
  call mpi_allreduce(data,stack_buff,1,MPI_DOUBLE_COMPLEX,&
         MPI_MIN,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  data=stack_buff
#endif
end subroutine mp_min_complex8_r0
