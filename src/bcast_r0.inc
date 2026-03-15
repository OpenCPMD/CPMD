subroutine mp_bcast_char_r0(data, root,comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  root,comm
#endif
  character(*), intent(inout)  :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_char_r0'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+len(data)*mp_char_in_bytes
  tim1=m_walltime()      
  call mpi_bcast(data,len(data),MPI_CHARACTER,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_char_r0
subroutine mp_bcast_int4_r0(data, root,comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  root,comm
#endif
  integer(int_4), intent(inout)  :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_int4_r0'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+1*mp_int4_in_bytes
  tim1=m_walltime()      
  call mpi_bcast(data,1,MPI_INTEGER,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_int4_r0
subroutine mp_bcast_int8_r0(data, root,comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  root,comm
#endif
  integer(int_8), intent(inout)  :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_int8_r0'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+1*mp_int8_in_bytes
  tim1=m_walltime()      
  call mpi_bcast(data,1,MPI_INTEGER8,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_int8_r0
subroutine mp_bcast_real8_r0(data, root,comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  root,comm
#endif
  real(real_8), intent(inout)  :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_real8_r0'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+1*mp_double_in_bytes
  tim1=m_walltime()      
  call mpi_bcast(data,1,MPI_DOUBLE_PRECISION,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_real8_r0
subroutine mp_bcast_complex8_r0(data, root,comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  root,comm
#endif
  complex(real_8), intent(inout)  :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_complex8_r0'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+1*mp_double_complex_in_bytes
  tim1=m_walltime()      
  call mpi_bcast(data,1,MPI_DOUBLE_COMPLEX,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_complex8_r0
subroutine mp_bcast_logical_r0(data, root,comm)
  ! Arguments
#ifdef __PARALLEL
  integer, intent(in) ::  root
  type(MPI_COMM), intent(in) :: comm
#else
  integer, intent(in) ::  root,comm
#endif
  logical, intent(inout)  :: data
  ! Variables
  integer :: ierr
  character(*),parameter :: procedureN='mp_bcast_logical_r0'
#ifdef __PARALLEL
  real(real_8) :: tim1,tim2
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+1*mp_logical_in_bytes
  tim1=m_walltime()      
  call mpi_bcast(data,1,MPI_LOGICAL,root,comm,ierr)
  call mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
end subroutine mp_bcast_logical_r0
