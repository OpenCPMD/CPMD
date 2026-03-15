#include "cpmd_global.h"

MODULE rnlsm_utils
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nlm
  USE rnlsm1_utils,                    ONLY: rnlsm1
  USE rnlsm2_utils,                    ONLY: rnlsm2
  USE rnlsmd_utils,                    ONLY: rnlsmd
  USE system,                          ONLY: cntl
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlsm

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlsm(c0,nstate,ikpt,ikind,tfor,only_dfnl,unpack_dfnl_fnl)
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8),INTENT(IN) __CONTIGUOUS  :: c0(:,:)
    INTEGER,INTENT(IN)                       :: nstate, ikpt, ikind
    LOGICAL,INTENT(IN)                       :: tfor
    LOGICAL,INTENT(IN),OPTIONAL              :: only_dfnl, unpack_dfnl_fnl

    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlsm'

    LOGICAL                                  :: dfnl, unpack
    INTEGER                                  :: isub

! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    IF(PRESENT(only_dfnl))THEN
       dfnl=only_dfnl
    ELSE
       dfnl=.FALSE.
    END IF
    IF(PRESENT(unpack_dfnl_fnl))THEN
       unpack=unpack_dfnl_fnl
    ELSE
       unpack=.TRUE.
    END IF
    IF (cntl%tfdist) THEN
       CALL rnlsmd(c0,nstate,ikind)
    ELSE
       IF(.NOT.dfnl) CALL rnlsm1(c0,nstate,ikind,fnl_unpack=unpack)
    ENDIF
    IF (tfor) CALL rnlsm2(c0,nstate,ikpt,ikind,dfnl_unpack=unpack)
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)

  END SUBROUTINE rnlsm
  ! ==================================================================


END MODULE rnlsm_utils
