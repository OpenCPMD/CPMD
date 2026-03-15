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
  LOGICAL      :: HFX_SCDM_STATUS, NEW_SCDM
!----------------------------------------------
  REAL(real_8) :: SCDM_CUTOFF, de_cutoff
  integer      :: n_loop
  REAL(real_8), ALLOCATABLE  ::  rho_scdm(:,:), grad_scdm(:,:)
!-------------------------------------------------------------------
  ! ==================================================================
  ! == INCLUDE FILE FOR LANGEVIN DYNAMICS WITH SCDM APPROACH        ==
  ! ==================================================================
  ! == GAMMA    VALUE OF LANGEVIN FRICTION COEFFICIENT IN A.U.      ==
  ! == T_BATH   TARGET TEMPERATURE IN KELVIN                        ==
  ! ==--------------------------------------------------------------==
  LOGICAL :: LANG_DYN
  !
  REAL(real_8) :: GAMMA,T_BATH,C_0,C_1,C_2,&
       C_RV,SIGMA_R,SIGMA_V
! ==================================================================
END MODULE ace_hfx
