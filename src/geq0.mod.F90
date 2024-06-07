#include "cpmd_global.h"

MODULE geq0mod
  IMPLICIT NONE

  ! ==================================================================
  ! == Used in parallel version otherwise GEQ0 is always .TRUE.     ==
  ! == GEQ0 = .TRUE. if the processors are the 0 component          ==
  ! ==================================================================
  LOGICAL :: geq0
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp declare target(geq0)
#endif
  ! ==================================================================

END MODULE geq0mod
