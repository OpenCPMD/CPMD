subroutine mp_send_char_r0(data, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) ::  dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  dest,tag,comm
#endif
  character(*), intent(in)  :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_char_r0'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+len(data)*mp_char_in_bytes
  tim1=m_walltime()
  call mpi_send(data,len(data),MPI_CHARACTER,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_char_r0

subroutine mp_recv_char_r0(data, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  src,tag,comm
#endif
  character(*), intent(out)  :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_char_r0'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+len(data)*mp_char_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,len(data),MPI_CHARACTER,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_char_r0
subroutine mp_send_char_r1(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  character(*), intent(in) , dimension(*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_char_r1'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+len(data)*n*mp_char_in_bytes
  tim1=m_walltime()
  call mpi_send(data,len(data)*n,MPI_CHARACTER,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_char_r1

subroutine mp_recv_char_r1(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  character(*), intent(out) , dimension(*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_char_r1'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+len(data)*n*mp_char_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,len(data)*n,MPI_CHARACTER,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_char_r1
subroutine mp_send_int4_r0(data, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) ::  dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  dest,tag,comm
#endif
  integer(int_4), intent(in)  :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_int4_r0'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+1*mp_int4_in_bytes
  tim1=m_walltime()
  call mpi_send(data,1,MPI_INTEGER,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_int4_r0

subroutine mp_recv_int4_r0(data, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  src,tag,comm
#endif
  integer(int_4), intent(out)  :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_int4_r0'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+1*mp_int4_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,1,MPI_INTEGER,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_int4_r0
subroutine mp_send_int4_r1(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  integer(int_4), intent(in) , dimension(*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_int4_r1'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_int4_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_INTEGER,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_int4_r1

subroutine mp_recv_int4_r1(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  integer(int_4), intent(out) , dimension(*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_int4_r1'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_int4_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_INTEGER,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_int4_r1
subroutine mp_send_int4_r2(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  integer(int_4), intent(in) , dimension(1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_int4_r2'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_int4_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_INTEGER,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_int4_r2

subroutine mp_recv_int4_r2(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  integer(int_4), intent(out) , dimension(1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_int4_r2'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_int4_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_INTEGER,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_int4_r2
subroutine mp_send_int8_r0(data, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) ::  dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  dest,tag,comm
#endif
  integer(int_8), intent(in)  :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_int8_r0'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+1*mp_int8_in_bytes
  tim1=m_walltime()
  call mpi_send(data,1,MPI_INTEGER8,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_int8_r0

subroutine mp_recv_int8_r0(data, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  src,tag,comm
#endif
  integer(int_8), intent(out)  :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_int8_r0'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+1*mp_int8_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,1,MPI_INTEGER8,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_int8_r0
subroutine mp_send_int8_r1(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  integer(int_8), intent(in) , dimension(*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_int8_r1'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_int8_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_INTEGER8,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_int8_r1

subroutine mp_recv_int8_r1(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  integer(int_8), intent(out) , dimension(*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_int8_r1'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_int8_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_INTEGER8,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_int8_r1
subroutine mp_send_real8_r0(data, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) ::  dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  dest,tag,comm
#endif
  real(real_8), intent(in)  :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_real8_r0'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+1*mp_double_in_bytes
  tim1=m_walltime()
  call mpi_send(data,1,MPI_DOUBLE_PRECISION,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_real8_r0

subroutine mp_recv_real8_r0(data, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  src,tag,comm
#endif
  real(real_8), intent(out)  :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_real8_r0'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+1*mp_double_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,1,MPI_DOUBLE_PRECISION,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_real8_r0
subroutine mp_send_real8_r1(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  real(real_8), intent(in) , dimension(*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_real8_r1'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_DOUBLE_PRECISION,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_real8_r1

subroutine mp_recv_real8_r1(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  real(real_8), intent(out) , dimension(*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_real8_r1'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_DOUBLE_PRECISION,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_real8_r1
subroutine mp_send_real8_r2(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  real(real_8), intent(in) , dimension(1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_real8_r2'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_DOUBLE_PRECISION,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_real8_r2

subroutine mp_recv_real8_r2(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  real(real_8), intent(out) , dimension(1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_real8_r2'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_DOUBLE_PRECISION,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_real8_r2
subroutine mp_send_real8_r3(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  real(real_8), intent(in) , dimension(1,1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_real8_r3'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_DOUBLE_PRECISION,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_real8_r3

subroutine mp_recv_real8_r3(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  real(real_8), intent(out) , dimension(1,1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_real8_r3'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_DOUBLE_PRECISION,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_real8_r3
subroutine mp_send_real8_r4(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  real(real_8), intent(in) , dimension(1,1,1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_real8_r4'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_DOUBLE_PRECISION,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_real8_r4

subroutine mp_recv_real8_r4(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  real(real_8), intent(out) , dimension(1,1,1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_real8_r4'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_DOUBLE_PRECISION,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_real8_r4
subroutine mp_send_real8_r5(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  real(real_8), intent(in) , dimension(1,1,1,1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_real8_r5'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_DOUBLE_PRECISION,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_real8_r5

subroutine mp_recv_real8_r5(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  real(real_8), intent(out) , dimension(1,1,1,1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_real8_r5'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_DOUBLE_PRECISION,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_real8_r5
subroutine mp_send_complex4_r0(data, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) ::  dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  dest,tag,comm
#endif
  complex(real_4), intent(in)  :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_complex4_r0'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+1*mp_complex_in_bytes
  tim1=m_walltime()
  call mpi_send(data,1,MPI_COMPLEX,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_complex4_r0

subroutine mp_recv_complex4_r0(data, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  src,tag,comm
#endif
  complex(real_4), intent(out)  :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_complex4_r0'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+1*mp_complex_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,1,MPI_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_complex4_r0
subroutine mp_send_complex4_r1(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  complex(real_4), intent(in) , dimension(*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_complex4_r1'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_complex_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_COMPLEX,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_complex4_r1

subroutine mp_recv_complex4_r1(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  complex(real_4), intent(out) , dimension(*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_complex4_r1'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_complex_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_complex4_r1
subroutine mp_send_complex4_r2(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  complex(real_4), intent(in) , dimension(1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_complex4_r2'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_complex_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_COMPLEX,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_complex4_r2

subroutine mp_recv_complex4_r2(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  complex(real_4), intent(out) , dimension(1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_complex4_r2'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_complex_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_complex4_r2
subroutine mp_send_complex4_r3(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  complex(real_4), intent(in) , dimension(1,1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_complex4_r3'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_complex_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_COMPLEX,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_complex4_r3

subroutine mp_recv_complex4_r3(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  complex(real_4), intent(out) , dimension(1,1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_complex4_r3'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_complex_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_complex4_r3
subroutine mp_send_complex4_r4(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  complex(real_4), intent(in) , dimension(1,1,1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_complex4_r4'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_complex_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_COMPLEX,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_complex4_r4

subroutine mp_recv_complex4_r4(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  complex(real_4), intent(out) , dimension(1,1,1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_complex4_r4'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_complex_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_complex4_r4
subroutine mp_send_complex8_r0(data, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) ::  dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  dest,tag,comm
#endif
  complex(real_8), intent(in)  :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_complex8_r0'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+1*mp_double_complex_in_bytes
  tim1=m_walltime()
  call mpi_send(data,1,MPI_DOUBLE_COMPLEX,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_complex8_r0

subroutine mp_recv_complex8_r0(data, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  src,tag,comm
#endif
  complex(real_8), intent(out)  :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_complex8_r0'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+1*mp_double_complex_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,1,MPI_DOUBLE_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_complex8_r0
subroutine mp_send_complex8_r1(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  complex(real_8), intent(in) , dimension(*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_complex8_r1'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_complex_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_DOUBLE_COMPLEX,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_complex8_r1

subroutine mp_recv_complex8_r1(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  complex(real_8), intent(out) , dimension(*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_complex8_r1'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_complex_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_DOUBLE_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_complex8_r1
subroutine mp_send_complex8_r2(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  complex(real_8), intent(in) , dimension(1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_complex8_r2'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_complex_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_DOUBLE_COMPLEX,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_complex8_r2

subroutine mp_recv_complex8_r2(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  complex(real_8), intent(out) , dimension(1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_complex8_r2'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_complex_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_DOUBLE_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_complex8_r2
subroutine mp_send_complex8_r3(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  complex(real_8), intent(in) , dimension(1,1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_complex8_r3'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_complex_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_DOUBLE_COMPLEX,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_complex8_r3

subroutine mp_recv_complex8_r3(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  complex(real_8), intent(out) , dimension(1,1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_complex8_r3'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_complex_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_DOUBLE_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_complex8_r3
subroutine mp_send_complex8_r4(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  complex(real_8), intent(in) , dimension(1,1,1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_complex8_r4'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_complex_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_DOUBLE_COMPLEX,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_complex8_r4

subroutine mp_recv_complex8_r4(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  complex(real_8), intent(out) , dimension(1,1,1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_complex8_r4'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_double_complex_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_DOUBLE_COMPLEX,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_complex8_r4
subroutine mp_send_logical_r0(data, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) ::  dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  dest,tag,comm
#endif
  logical, intent(in)  :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_logical_r0'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+1*mp_logical_in_bytes
  tim1=m_walltime()
  call mpi_send(data,1,MPI_LOGICAL,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_logical_r0

subroutine mp_recv_logical_r0(data, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  src,tag,comm
#endif
  logical, intent(out)  :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_logical_r0'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+1*mp_logical_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,1,MPI_LOGICAL,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_logical_r0
subroutine mp_send_logical_r1(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  logical, intent(in) , dimension(*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_logical_r1'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_logical_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_LOGICAL,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_logical_r1

subroutine mp_recv_logical_r1(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  logical, intent(out) , dimension(*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_logical_r1'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_logical_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_LOGICAL,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_logical_r1
subroutine mp_send_logical_r2(data,n, dest,tag,comm)
  ! Arguments 
#ifdef __PARALLEL
  integer, intent(in) :: n, dest,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, dest,tag,comm
#endif
  logical, intent(in) , dimension(1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_isend_logical_r2'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_logical_in_bytes
  tim1=m_walltime()
  call mpi_send(data,n,MPI_LOGICAL,dest,tag,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_send_logical_r2

subroutine mp_recv_logical_r2(data,n, src,tag,comm)
  use pstat
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) :: n, src,tag
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) :: n, src,tag,comm
#endif
  logical, intent(out) , dimension(1,*) :: data
  ! Variables
  integer :: ierr
  real(real_8) :: tim1,tim2
  character(*),parameter::procedureN='mp_recv_logical_r2'
#ifdef __PARALLEL
  cmcal(ipar_send)=cmcal(ipar_send)+1.0d0
  cmlen(ipar_send)=cmlen(ipar_send)+n*mp_logical_in_bytes
  tim1=m_walltime()
  call mpi_recv(data,n,MPI_LOGICAL,src,tag,comm,&
        MPI_STATUS_IGNORE,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_send)=cmtim(ipar_send)+tim2-tim1
#endif
end subroutine mp_recv_logical_r2
