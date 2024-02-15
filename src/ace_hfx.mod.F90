MODULE ace_hfx
  USE kinds,                           ONLY: real_8
  !
  IMPLICIT NONE
  SAVE
  ! ==================================================================
  ! == INCLUDE FILE FOR ACE+MTS CALCULATIONS                        ==
  ! ==================================================================
  ! == XI         :: projectors of ACE formulation                  ==
  ! == use_ace    :: Using ACE formulation or not                   ==
  ! == switch_ace :: decide what to do in HFX calculation           ==
  ! == status_ace :: low level or high level force calculation      ==
  ! ==--------------------------------------------------------------==
  LOGICAL :: use_ace, switch_ace, status_ace
  !
  COMPLEX(real_8), ALLOCATABLE  ::      XI(:,:)
! ==================================================================
  LOGICAL      :: HFX_SCDM_STATUS
!----------------------------------------------
  REAL(real_8) :: SCDM_CUTOFF

!-------------------------------------------------------------------
END MODULE ace_hfx
