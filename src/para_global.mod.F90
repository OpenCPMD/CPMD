#include "cpmd_global.h"

MODULE para_global
  USE kinds,                           ONLY: real_8,&
                                             int_8

  IMPLICIT NONE

  PRIVATE

  ! allows MPI to allocate memory for reduction operations 
  LOGICAL, PUBLIC :: para_use_mpi_in_place=.FALSE.

  ! maximum buffer size for reduction operations (require an allocation)
  INTEGER(int_8), PUBLIC     :: il_para_buff(1)=2**16
  INTEGER, PUBLIC            :: buff_size_in_bytes=STORAGE_SIZE(para_buff)/STORAGE_SIZE('A')
#ifdef _USE_SCRATCHLIBRARY
  COMPLEX(real_8), PUBLIC, POINTER __CONTIGUOUS  :: para_buff(:)
#else
  COMPLEX(real_8), PUBLIC, ALLOCATABLE, TARGET   :: para_buff(:)
#endif

END MODULE para_global
