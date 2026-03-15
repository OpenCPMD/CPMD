subroutine mp_sendrecv_int4_r1(data_in,n_in,dest,data_out,n_out,src,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n_in,n_out,src,dest
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n_in,n_out,src,dest,comm
#endif
  integer(int_4), intent(in) , dimension(*) :: data_in
  integer(int_4), intent(out) , dimension(*) :: data_out
  ! Variables
  integer :: ierr,tag_in,tag_out
  character(*),parameter::procedureN='mp_sendrecv_int4_r1'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+(n_in+n_out)*mp_int4_in_bytes
  tim1=m_walltime()
  tag_in=1;tag_out=1
  call mpi_sendrecv(data_in,n_in,MPI_INTEGER,dest,tag_in,&
       data_out,n_out,MPI_INTEGER,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#else
  call mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
end subroutine mp_sendrecv_int4_r1
subroutine mp_sendrecv_int4_r2(data_in,n_in,dest,data_out,n_out,src,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n_in,n_out,src,dest
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n_in,n_out,src,dest,comm
#endif
  integer(int_4), intent(in) , dimension(1,*) :: data_in
  integer(int_4), intent(out) , dimension(1,*) :: data_out
  ! Variables
  integer :: ierr,tag_in,tag_out
  character(*),parameter::procedureN='mp_sendrecv_int4_r2'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+(n_in+n_out)*mp_int4_in_bytes
  tim1=m_walltime()
  tag_in=1;tag_out=1
  call mpi_sendrecv(data_in,n_in,MPI_INTEGER,dest,tag_in,&
       data_out,n_out,MPI_INTEGER,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#else
  call mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
end subroutine mp_sendrecv_int4_r2
subroutine mp_sendrecv_int8_r1(data_in,n_in,dest,data_out,n_out,src,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n_in,n_out,src,dest
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n_in,n_out,src,dest,comm
#endif
  integer(int_8), intent(in) , dimension(*) :: data_in
  integer(int_8), intent(out) , dimension(*) :: data_out
  ! Variables
  integer :: ierr,tag_in,tag_out
  character(*),parameter::procedureN='mp_sendrecv_int8_r1'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+(n_in+n_out)*mp_int8_in_bytes
  tim1=m_walltime()
  tag_in=1;tag_out=1
  call mpi_sendrecv(data_in,n_in,MPI_INTEGER8,dest,tag_in,&
       data_out,n_out,MPI_INTEGER8,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#else
  call mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
end subroutine mp_sendrecv_int8_r1
subroutine mp_sendrecv_real8_r1(data_in,n_in,dest,data_out,n_out,src,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n_in,n_out,src,dest
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n_in,n_out,src,dest,comm
#endif
  real(real_8), intent(in) , dimension(*) :: data_in
  real(real_8), intent(out) , dimension(*) :: data_out
  ! Variables
  integer :: ierr,tag_in,tag_out
  character(*),parameter::procedureN='mp_sendrecv_real8_r1'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+(n_in+n_out)*mp_double_in_bytes
  tim1=m_walltime()
  tag_in=1;tag_out=1
  call mpi_sendrecv(data_in,n_in,MPI_DOUBLE_PRECISION,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_PRECISION,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#else
  call mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
end subroutine mp_sendrecv_real8_r1
subroutine mp_sendrecv_real8_r2(data_in,n_in,dest,data_out,n_out,src,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n_in,n_out,src,dest
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n_in,n_out,src,dest,comm
#endif
  real(real_8), intent(in) , dimension(1,*) :: data_in
  real(real_8), intent(out) , dimension(1,*) :: data_out
  ! Variables
  integer :: ierr,tag_in,tag_out
  character(*),parameter::procedureN='mp_sendrecv_real8_r2'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+(n_in+n_out)*mp_double_in_bytes
  tim1=m_walltime()
  tag_in=1;tag_out=1
  call mpi_sendrecv(data_in,n_in,MPI_DOUBLE_PRECISION,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_PRECISION,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#else
  call mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
end subroutine mp_sendrecv_real8_r2
subroutine mp_sendrecv_real8_r3(data_in,n_in,dest,data_out,n_out,src,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n_in,n_out,src,dest
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n_in,n_out,src,dest,comm
#endif
  real(real_8), intent(in) , dimension(1,1,*) :: data_in
  real(real_8), intent(out) , dimension(1,1,*) :: data_out
  ! Variables
  integer :: ierr,tag_in,tag_out
  character(*),parameter::procedureN='mp_sendrecv_real8_r3'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+(n_in+n_out)*mp_double_in_bytes
  tim1=m_walltime()
  tag_in=1;tag_out=1
  call mpi_sendrecv(data_in,n_in,MPI_DOUBLE_PRECISION,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_PRECISION,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#else
  call mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
end subroutine mp_sendrecv_real8_r3
subroutine mp_sendrecv_real8_r4(data_in,n_in,dest,data_out,n_out,src,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n_in,n_out,src,dest
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n_in,n_out,src,dest,comm
#endif
  real(real_8), intent(in) , dimension(1,1,1,*) :: data_in
  real(real_8), intent(out) , dimension(1,1,1,*) :: data_out
  ! Variables
  integer :: ierr,tag_in,tag_out
  character(*),parameter::procedureN='mp_sendrecv_real8_r4'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+(n_in+n_out)*mp_double_in_bytes
  tim1=m_walltime()
  tag_in=1;tag_out=1
  call mpi_sendrecv(data_in,n_in,MPI_DOUBLE_PRECISION,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_PRECISION,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#else
  call mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
end subroutine mp_sendrecv_real8_r4
subroutine mp_sendrecv_complex8_r1(data_in,n_in,dest,data_out,n_out,src,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n_in,n_out,src,dest
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n_in,n_out,src,dest,comm
#endif
  complex(real_8), intent(in) , dimension(*) :: data_in
  complex(real_8), intent(out) , dimension(*) :: data_out
  ! Variables
  integer :: ierr,tag_in,tag_out
  character(*),parameter::procedureN='mp_sendrecv_complex8_r1'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+(n_in+n_out)*mp_double_complex_in_bytes
  tim1=m_walltime()
  tag_in=1;tag_out=1
  call mpi_sendrecv(data_in,n_in,MPI_DOUBLE_COMPLEX,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_COMPLEX,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#else
  call mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
end subroutine mp_sendrecv_complex8_r1
subroutine mp_sendrecv_complex8_r2(data_in,n_in,dest,data_out,n_out,src,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n_in,n_out,src,dest
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n_in,n_out,src,dest,comm
#endif
  complex(real_8), intent(in) , dimension(1,*) :: data_in
  complex(real_8), intent(out) , dimension(1,*) :: data_out
  ! Variables
  integer :: ierr,tag_in,tag_out
  character(*),parameter::procedureN='mp_sendrecv_complex8_r2'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+(n_in+n_out)*mp_double_complex_in_bytes
  tim1=m_walltime()
  tag_in=1;tag_out=1
  call mpi_sendrecv(data_in,n_in,MPI_DOUBLE_COMPLEX,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_COMPLEX,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#else
  call mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
end subroutine mp_sendrecv_complex8_r2
subroutine mp_sendrecv_complex8_r3(data_in,n_in,dest,data_out,n_out,src,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n_in,n_out,src,dest
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n_in,n_out,src,dest,comm
#endif
  complex(real_8), intent(in) , dimension(1,1,*) :: data_in
  complex(real_8), intent(out) , dimension(1,1,*) :: data_out
  ! Variables
  integer :: ierr,tag_in,tag_out
  character(*),parameter::procedureN='mp_sendrecv_complex8_r3'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+(n_in+n_out)*mp_double_complex_in_bytes
  tim1=m_walltime()
  tag_in=1;tag_out=1
  call mpi_sendrecv(data_in,n_in,MPI_DOUBLE_COMPLEX,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_COMPLEX,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#else
  call mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
end subroutine mp_sendrecv_complex8_r3
subroutine mp_sendrecv_complex8_r4(data_in,n_in,dest,data_out,n_out,src,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n_in,n_out,src,dest
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n_in,n_out,src,dest,comm
#endif
  complex(real_8), intent(in) , dimension(1,1,1,*) :: data_in
  complex(real_8), intent(out) , dimension(1,1,1,*) :: data_out
  ! Variables
  integer :: ierr,tag_in,tag_out
  character(*),parameter::procedureN='mp_sendrecv_complex8_r4'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+(n_in+n_out)*mp_double_complex_in_bytes
  tim1=m_walltime()
  tag_in=1;tag_out=1
  call mpi_sendrecv(data_in,n_in,MPI_DOUBLE_COMPLEX,dest,tag_in,&
       data_out,n_out,MPI_DOUBLE_COMPLEX,src,tag_out,comm,&
       MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#else
  call mp_mpi_error_assert(1,procedureN,__LINE__,__FILE__)
#endif
end subroutine mp_sendrecv_complex8_r4
