subroutine mp_prod_int4_r0(data_in,data_out, comm)
  ! Arguments
#ifdef __PARALLEL
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: comm
#endif
  integer(int_4), intent(in)  :: data_in
  integer(int_4), intent(inout)  :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_int4_r0'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,1,MPI_INTEGER,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out = data_in 
#endif
end subroutine mp_prod_int4_r0
subroutine mp_prod_real8_r0(data_in,data_out, comm)
  ! Arguments
#ifdef __PARALLEL
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: comm
#endif
  real(real_8), intent(in)  :: data_in
  real(real_8), intent(inout)  :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_real8_r0'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,1,MPI_DOUBLE_PRECISION,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out = data_in 
#endif
end subroutine mp_prod_real8_r0
subroutine mp_prod_complex8_r0(data_in,data_out, comm)
  ! Arguments
#ifdef __PARALLEL
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: comm
#endif
  complex(real_8), intent(in)  :: data_in
  complex(real_8), intent(inout)  :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_complex8_r0'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,1,MPI_DOUBLE_COMPLEX,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out = data_in 
#endif
end subroutine mp_prod_complex8_r0
