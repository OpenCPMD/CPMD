MODULE cvan
  USE kinds,                           ONLY: real_8

  IMPLICIT NONE

  ! ==================================================================
  ! == VANDERBILD PSEUDOPOTENTIALS                                  ==
  ! ==--------------------------------------------------------------==
  INTEGER, ALLOCATABLE :: nelev(:)

  REAL(real_8), ALLOCATABLE :: qq(:,:,:)
  REAL(real_8), ALLOCATABLE :: dvan(:,:,:)
  REAL(real_8), ALLOCATABLE :: deeq(:,:,:,:)
  REAL(real_8), ALLOCATABLE :: deeq_fnl_hfx(:,:,:)
  COMPLEX(real_8), ALLOCATABLE, TARGET :: qg(:,:)
  COMPLEX(real_8), ALLOCATABLE, TARGET :: qg_dipole(:,:)
  ! ==================================================================

END MODULE cvan
