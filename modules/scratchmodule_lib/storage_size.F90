program stor
  use iso_fortran_env
  implicit none

  complex(real64) :: a
  integer(int32), parameter :: size_of_real=storage_size(cmplx(1.0_real64,1.0_real64,kind=real64))/storage_size('C')

  write(*,*) size_of_real

end program stor
