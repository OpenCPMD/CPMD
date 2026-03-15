subroutine mp_prod_int4_r1(data_in,data_out,n, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n 
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, comm
#endif
  integer(int_4), intent(in), dimension(*) :: data_in
  integer(int_4), intent(inout), dimension(*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_int4_r1'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,n,MPI_INTEGER,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1:n)=data_in (1:n)
#endif
end subroutine mp_prod_int4_r1
subroutine mp_prod_int4_r2(data_in,data_out,n, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n 
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, comm
#endif
  integer(int_4), intent(in), dimension(1,*) :: data_in
  integer(int_4), intent(inout), dimension(1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_int4_r2'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,n,MPI_INTEGER,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1,1:n)=data_in (1,1:n)
#endif
end subroutine mp_prod_int4_r2
subroutine mp_prod_int4_r3(data_in,data_out,n, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n 
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, comm
#endif
  integer(int_4), intent(in), dimension(1,1,*) :: data_in
  integer(int_4), intent(inout), dimension(1,1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_int4_r3'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,n,MPI_INTEGER,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1,1,1:n)=data_in (1,1,1:n)
#endif
end subroutine mp_prod_int4_r3
subroutine mp_prod_int4_r4(data_in,data_out,n, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n 
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, comm
#endif
  integer(int_4), intent(in), dimension(1,1,1,*) :: data_in
  integer(int_4), intent(inout), dimension(1,1,1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_int4_r4'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,n,MPI_INTEGER,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1,1,1,1:n)=data_in (1,1,1,1:n)
#endif
end subroutine mp_prod_int4_r4
subroutine mp_prod_real8_r1(data_in,data_out,n, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n 
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, comm
#endif
  real(real_8), intent(in), dimension(*) :: data_in
  real(real_8), intent(inout), dimension(*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_real8_r1'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,n,MPI_DOUBLE_PRECISION,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1:n)=data_in (1:n)
#endif
end subroutine mp_prod_real8_r1
subroutine mp_prod_real8_r2(data_in,data_out,n, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n 
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, comm
#endif
  real(real_8), intent(in), dimension(1,*) :: data_in
  real(real_8), intent(inout), dimension(1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_real8_r2'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,n,MPI_DOUBLE_PRECISION,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1,1:n)=data_in (1,1:n)
#endif
end subroutine mp_prod_real8_r2
subroutine mp_prod_real8_r3(data_in,data_out,n, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n 
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, comm
#endif
  real(real_8), intent(in), dimension(1,1,*) :: data_in
  real(real_8), intent(inout), dimension(1,1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_real8_r3'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,n,MPI_DOUBLE_PRECISION,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1,1,1:n)=data_in (1,1,1:n)
#endif
end subroutine mp_prod_real8_r3
subroutine mp_prod_real8_r4(data_in,data_out,n, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n 
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, comm
#endif
  real(real_8), intent(in), dimension(1,1,1,*) :: data_in
  real(real_8), intent(inout), dimension(1,1,1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_real8_r4'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,n,MPI_DOUBLE_PRECISION,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1,1,1,1:n)=data_in (1,1,1,1:n)
#endif
end subroutine mp_prod_real8_r4
subroutine mp_prod_complex8_r1(data_in,data_out,n, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n 
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, comm
#endif
  complex(real_8), intent(in), dimension(*) :: data_in
  complex(real_8), intent(inout), dimension(*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_complex8_r1'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,n,MPI_DOUBLE_COMPLEX,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1:n)=data_in (1:n)
#endif
end subroutine mp_prod_complex8_r1
subroutine mp_prod_complex8_r2(data_in,data_out,n, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n 
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, comm
#endif
  complex(real_8), intent(in), dimension(1,*) :: data_in
  complex(real_8), intent(inout), dimension(1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_complex8_r2'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,n,MPI_DOUBLE_COMPLEX,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1,1:n)=data_in (1,1:n)
#endif
end subroutine mp_prod_complex8_r2
subroutine mp_prod_complex8_r3(data_in,data_out,n, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n 
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, comm
#endif
  complex(real_8), intent(in), dimension(1,1,*) :: data_in
  complex(real_8), intent(inout), dimension(1,1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_complex8_r3'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,n,MPI_DOUBLE_COMPLEX,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1,1,1:n)=data_in (1,1,1:n)
#endif
end subroutine mp_prod_complex8_r3
subroutine mp_prod_complex8_r4(data_in,data_out,n, comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n 
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, comm
#endif
  complex(real_8), intent(in), dimension(1,1,1,*) :: data_in
  complex(real_8), intent(inout), dimension(1,1,1,*) :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_prod_complex8_r4'
#ifdef __PARALLEL
  integer :: ierr
  call mpi_allreduce(data_in,data_out,n,MPI_DOUBLE_COMPLEX,&
       MPI_PROD,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  data_out (1,1,1,1:n)=data_in (1,1,1,1:n)
#endif
end subroutine mp_prod_complex8_r4
