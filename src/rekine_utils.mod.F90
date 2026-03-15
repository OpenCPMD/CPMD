#include "cpmd_global.h"

MODULE rekine_utils
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes
  USE dotp_utils,                      ONLY: dotp_c1_cp
  USE geq0mod,                         ONLY: geq0
  USE harm,                            ONLY: xmu
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai
  USE reshaper,                        ONLY: reshape_inplace
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
#ifdef __PARALLEL
  USE mpi_f08
#endif
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rekine

CONTAINS

  ! ==================================================================
  SUBROUTINE rekine(cm,nstate,ekinc,use_cp_grps)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate
    COMPLEX(real_8),INTENT(IN)               :: cm(ncpw%ngw,nstate)
    REAL(real_8),INTENT(OUT)                 :: ekinc
    LOGICAL,INTENT(IN),OPTIONAL              :: use_cp_grps

#ifdef __PARALLEL
    INTEGER                                  :: i, ig, isub, ngw_local,&
                                                ibeg_c0, iend_c0
    type(MPI_COMM)                           :: gid
#else
    INTEGER                                  :: i, ig, isub, ngw_local,&
                                                ibeg_c0, iend_c0, gid
#endif
    REAL(real_8),POINTER __CONTIGUOUS        :: cm_r(:,:)
    LOGICAL                                  :: geq0_local, cp_active
    REAL(real_8)                             :: ax, bx, pf
    CHARACTER(*), PARAMETER                  :: procedureN = 'rekine'
! ==--------------------------------------------------------------==
! ==  COMPUTE FICTITIOUS KINETIC ENERGY OF THE ELECTRONS          ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF(PRESENT(use_cp_grps))THEN
       cp_active=use_cp_grps
    ELSE
       cp_active=.FALSE.
    END IF
    IF(cp_active)THEN
       CALL cp_grp_get_sizes(ngw_l=ngw_local,first_g=ibeg_c0,last_g=iend_c0,geq0_l=geq0_local)
       gid=parai%cp_grp
    ELSE
       ngw_local=ncpw%ngw
       ibeg_c0=1
       iend_c0=ncpw%ngw
       gid=parai%allgrp
       geq0_local=geq0
    END IF

    ekinc=0._real_8
    IF (cntl%tmass) THEN
       !$omp parallel do private(i,ig,pf,ax,bx) reduction(+:ekinc)
       DO i=1,nstate
          DO ig=ibeg_c0,iend_c0
             pf=2.0_real_8*xmu(ig)
             ax=REAL(cm(ig,i))
             bx=AIMAG(cm(ig,i))
             ekinc=ekinc+pf*(ax*ax+bx*bx)
          ENDDO
          IF (geq0_local) ekinc=ekinc-xmu(1)*REAL(cm(1,i))*REAL(cm(1,i))
       ENDDO
    ELSE
       CALL reshape_inplace(cm,(/ncpw%ngw*2,nstate/),cm_r)
       CALL calc_ekinc(cm_r,ekinc,ibeg_c0*2-1,iend_c0*2,nstate,geq0_local)
       ekinc=ekinc*cntr%emass
    ENDIF
    CALL mp_sum(ekinc,gid)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rekine
  ! ==================================================================
  SUBROUTINE calc_ekinc(cm_r,ekinc,ibeg_c0,iend_c0,nstate,geq0_local)
    ! ==--------------------------------------------------------------==
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: cm_r(:,:)
    REAL(real_8),INTENT(OUT)                 :: ekinc
    INTEGER,INTENT(IN)                       :: nstate, ibeg_c0, iend_c0
    LOGICAL,INTENT(IN)                       :: geq0_local

    INTEGER                                  :: i,ig
    REAL(real_8)                             :: ekinc_s
    
    ekinc=0.0_real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute &
#else
    !$omp parallel do &
#endif
    !$omp& private(i,ekinc_s,ig) reduction(+:ekinc)
    DO i=1,nstate
       IF(geq0_local)THEN
          ekinc_s=cm_r(ibeg_c0,i)**2*0.5_real_8
       ELSE
          ekinc_s=cm_r(ibeg_c0,i)**2
          ekinc_s=ekinc_s+cm_r(ibeg_c0+1,i)**2
       END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do simd private(ig)  &
#else
       !$omp simd &
#endif
       !$omp& reduction(+:ekinc_s)
       DO ig=ibeg_c0+2,iend_c0
          ekinc_s=ekinc_s+cm_r(ig,i)**2
       END DO
       ekinc=ekinc+ekinc_s*2.0_real_8
    END DO

  END SUBROUTINE calc_ekinc
  ! ================================================================== 
END MODULE rekine_utils
