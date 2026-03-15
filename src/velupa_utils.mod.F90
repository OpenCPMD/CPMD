#include "cpmd_global.h"

MODULE velupa_utils
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes
  USE harm,                            ONLY: dtan2c,&
                                             dtan2w,&
                                             freq,&
                                             xmu
  USE kinds,                           ONLY: real_8
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             ncpw,&
                                             nkpt
  USE reshaper,                        ONLY: reshape_inplace
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt_elec,&
                                             dtb2me

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: velupa

CONTAINS

  ! ==================================================================
  SUBROUTINE velupa(c0,cm,c2,nstate,nnlst,use_cp_grps)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(IN)               :: c0(nkpt%ngwk,*), &
                                                c2(nkpt%ngwk,*)
    COMPLEX(real_8),INTENT(INOUT)            :: cm(nkpt%ngwk,*)
    INTEGER,INTENT(IN)                       :: nstate, nnlst
    LOGICAL,INTENT(IN),OPTIONAL              :: use_cp_grps

    LOGICAL                                  :: cp_active
    INTEGER                                  :: i, ig, isub, &
                                                ibeg_c0, iend_c0, ngwk_local
    REAL(real_8)                             :: dtx, hfo, hgi
    REAL(real_8),POINTER __CONTIGUOUS        :: c2_r(:,:),cm_r(:,:)
    CHARACTER(*), PARAMETER                  :: procedureN = 'velupa'

    CALL tiset(procedureN,isub)
    IF(PRESENT(use_cp_grps))THEN
       cp_active=use_cp_grps
    ELSE
       cp_active=.FALSE.
    END IF
    IF(cp_active)THEN
       CALL cp_grp_get_sizes(ngwk_l=ngwk_local,firstk_g=ibeg_c0,lastk_g=iend_c0)
    ELSE
       ngwk_local=nkpt%ngwk
       ibeg_c0=1
       iend_c0=nkpt%ngwk
    END IF
    ! WAVEFUNCTION VELOCITY UPDATE
    IF (cntl%tharm) THEN
       !$omp parallel do private(i,ig,hfo,hgi)
       DO i=1,nstate
          DO ig=ibeg_c0,iend_c0
             IF (cntl%tmass) THEN
                hfo=-xmu(ig)*freq(ig)**2
             ELSE
                hfo=-cntr%emass*freq(ig)**2
             END IF
             IF (nnlst.EQ.1) THEN
                hgi=dtan2w(ig)
                cm(ig,i)=cm(ig,i)+hgi*(c2(ig,i)-hfo*c0(ig,i))
             ELSE
                hgi=dtan2c(ig)
                cm(ig,i)=cm(ig,i)+hgi*c2(ig,i)
             END IF
          END DO
       END DO
    ELSE IF (cntl%tmass) THEN
       !$omp parallel do private(i,ig,hgi)
       DO i=1,nstate
          DO ig=ibeg_c0,iend_c0
             hgi=REAL(nnlst,kind=real_8)*0.5_real_8*dt_elec/xmu(ig)
             cm(ig,i)=cm(ig,i)+hgi*c2(ig,i)
          END DO
       END DO
    ELSE
       dtx=REAL(nnlst,kind=real_8)*dtb2me
       CALL reshape_inplace(cm,(/2*ncpw%ngw,nstate/),cm_r)
       CALL reshape_inplace(c2,(/2*ncpw%ngw,nstate/),c2_r)
       CALL update_cm(cm_r,c2_r,nstate,dtx,ibeg_c0*2-1,iend_c0*2)
    END IF

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE velupa
  ! ==================================================================
  SUBROUTINE update_cm(cm_r,c2_r,nstate,dtx,ibeg_c0,iend_c0)
    ! ==--------------------------------------------------------------==
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: cm_r(:,:)
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: c2_r(:,:)
    REAL(real_8),INTENT(IN)                  :: dtx
    INTEGER,INTENT(IN)                       :: nstate, ibeg_c0, iend_c0

    INTEGER                                  :: i,ig

#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute parallel do simd collapse(2)&
#else
    !$omp parallel do &
#endif
    !$omp& private(i,ig)
    DO i=1,nstate
       DO ig=ibeg_c0,iend_c0
          cm_r(ig,i)=cm_r(ig,i)+dtx*c2_r(ig,i)
       END DO
    END DO

    ! ==--------------------------------------------------------------==
  END SUBROUTINE update_cm
  ! ==================================================================

END MODULE velupa_utils
