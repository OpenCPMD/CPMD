subroutine mp_sum_int4_r0(data_in,data_out,comm)
  ! Arguments
#ifdef __PARALLEL
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  comm
#endif
  integer(int_4), intent(in) :: data_in
  integer(int_4), intent(inout)  :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_int4_r0'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  integer :: ierr,count
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+1*mp_int4_in_bytes
  tim1=m_walltime()      
  call mpi_allreduce(data_in,data_out,1,MPI_INTEGER,&
       MPI_SUM,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out=data_in
#endif
end subroutine mp_sum_int4_r0
subroutine mp_sum_real8_r0(data_in,data_out,comm)
  ! Arguments
#ifdef __PARALLEL
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  comm
#endif
  real(real_8), intent(in) :: data_in
  real(real_8), intent(inout)  :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_real8_r0'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  integer :: ierr,count
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+1*mp_double_in_bytes
  tim1=m_walltime()      
  call mpi_allreduce(data_in,data_out,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out=data_in
#endif
end subroutine mp_sum_real8_r0
subroutine mp_sum_real4_r0(data_in,data_out,comm)
  ! Arguments
#ifdef __PARALLEL
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  comm
#endif
  real(real_4), intent(in) :: data_in
  real(real_4), intent(inout)  :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_real4_r0'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  integer :: ierr,count
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+1*mp_real_in_bytes
  tim1=m_walltime()      
  call mpi_allreduce(data_in,data_out,1,MPI_REAL,&
       MPI_SUM,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out=data_in
#endif
end subroutine mp_sum_real4_r0
subroutine mp_sum_complex8_r0(data_in,data_out,comm)
  ! Arguments
#ifdef __PARALLEL
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  comm
#endif
  complex(real_8), intent(in) :: data_in
  complex(real_8), intent(inout)  :: data_out
  ! Variables
  character(*),parameter :: procedureN='mp_sum_complex8_r0'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  integer :: ierr,count
  cmcal(ipar_gsum)=cmcal(ipar_gsum)+1.0d0
  cmlen(ipar_gsum)=cmlen(ipar_gsum)+1*mp_double_complex_in_bytes
  tim1=m_walltime()      
  call mpi_allreduce(data_in,data_out,1,MPI_DOUBLE_COMPLEX,&
       MPI_SUM,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_gsum)=cmtim(ipar_gsum)+tim2-tim1
#else
  data_out=data_in
#endif
end subroutine mp_sum_complex8_r0
