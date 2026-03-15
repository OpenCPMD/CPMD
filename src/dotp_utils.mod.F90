#include "cpmd_global.h"

MODULE dotp_utils

  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dotp
  PUBLIC :: dotp_c1_cp
  PUBLIC :: dotp_c2_cp

  INTERFACE dotp
     MODULE PROCEDURE dotp_c
  END INTERFACE dotp

CONTAINS

  ! ==================================================================
  FUNCTION dotp_c(n,a,b) RESULT(dotp)
    INTEGER, INTENT(IN)                      :: n
    COMPLEX(real_8), INTENT(IN) __CONTIGUOUS :: a(:), b(:)
    REAL(real_8)                             :: dotp

    CHARACTER(*), PARAMETER                  :: procedureN = 'dotp_c'

    REAL(real_8), EXTERNAL                   :: ddot

! ==--------------------------------------------------------------==
    dotp = 0._real_8
    IF (n.EQ.0) RETURN

    IF (n>SIZE(a).OR.n>SIZE(b)) CALL stopgm(procedureN,&
         'wrong dimension',&
         __LINE__,__FILE__)

    IF (geq0) THEN
       dotp = REAL(a(1),KIND=real_8) * REAL(b(1),KIND=real_8)
    ELSE
       dotp=2.0_real_8*( REAL(a(1),KIND=real_8)*REAL(b(1),KIND=real_8) + &
            AIMAG(a(1))*AIMAG(b(1)))
    ENDIF
    IF (n.GT.1) THEN
       dotp=dotp+2.0_real_8*ddot(2*n-2,a(2),1,b(2),1)
    ENDIF
  END FUNCTION dotp_c
  ! ==================================================================
  FUNCTION dotp_c1_cp(n,a,geq0_loc) RESULT(dotp_cp)
    INTEGER,INTENT(IN)                         :: n
    COMPLEX(real_8),INTENT(IN)                 :: a(*)
    REAL(real_8)                               :: dotp_cp
    LOGICAL,INTENT(IN)                         :: geq0_loc

    REAL(real_8), EXTERNAL                     :: ddot
    ! ==--------------------------------------------------------------==
    dotp_cp = 0._real_8
    IF (n.EQ.0) RETURN

    IF (geq0_loc) THEN
       dotp_cp = REAL(a(1),KIND=real_8) * REAL(a(1),KIND=real_8)
    ELSE
       dotp_cp=2.0_real_8*( REAL(a(1),KIND=real_8)*REAL(a(1),KIND=real_8) + &
            AIMAG(a(1))*AIMAG(a(1)))
    ENDIF
    IF (n.GT.1) THEN
       dotp_cp=dotp_cp+2.0_real_8*ddot(2*n-2,a(2),1,a(2),1)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION dotp_c1_cp
  ! ==================================================================
  FUNCTION dotp_c2_cp(n,a,b,geq0_loc) RESULT(dotp_cp)
    INTEGER,INTENT(IN)                         :: n
    COMPLEX(real_8),INTENT(IN)                 :: a(*), b(*)
    REAL(real_8)                               :: dotp_cp
    LOGICAL,INTENT(IN)                         :: geq0_loc

    REAL(real_8), EXTERNAL                     :: ddot
    ! ==--------------------------------------------------------------==
    dotp_cp = 0._real_8
    IF (n.EQ.0) RETURN

    IF (geq0_loc) THEN
       dotp_cp = REAL(a(1),KIND=real_8) * REAL(b(1),KIND=real_8)
    ELSE
       dotp_cp=2.0_real_8*( REAL(a(1),KIND=real_8)*REAL(b(1),KIND=real_8) + &
            AIMAG(a(1))*AIMAG(b(1)))
    ENDIF
    IF (n.GT.1) THEN
       dotp_cp=dotp_cp+2.0_real_8*ddot(2*n-2,a(2),1,b(2),1)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END FUNCTION dotp_c2_cp
  ! ==================================================================
END MODULE dotp_utils

! ==================================================================
FUNCTION dotp(n,a,b)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_8
  USE geq0mod , ONLY:geq0
  IMPLICIT NONE
  INTEGER                                    :: n
  REAL(real_8)                               :: a(2*n), b(2*n), dotp

#if defined (_vpp_) || defined (__SR8000) || defined (__ES)
  INTEGER                                    :: ii
  REAL(real_8)                               :: ddo
#else
  REAL(real_8), EXTERNAL                     :: ddot
#endif
  ! ==--------------------------------------------------------------==
  IF (n.EQ.0) THEN
     dotp=0._real_8
     RETURN
  ENDIF
  IF (geq0) THEN
     dotp=a(1)*b(1)
  ELSE
     dotp=2.0_real_8*(a(1)*b(1)+a(2)*b(2))
  ENDIF
  IF (n.GT.1) THEN
!! cmb - We keep the following for platforms with no ddot
#if defined (_vpp_) || defined (__SR8000)  || defined (__ES)
     ddo=0._real_8
     !$omp parallel do private(II) reduction(+:DDO)
     DO ii=3,2*n
        ddo=ddo+a(ii)*b(ii)
     ENDDO
     dotp=dotp+ddo+ddo
#else
     dotp=dotp+2.0_real_8*ddot(2*n-2,a(3),1,b(3),1)
#endif
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
END FUNCTION dotp
! ==================================================================
