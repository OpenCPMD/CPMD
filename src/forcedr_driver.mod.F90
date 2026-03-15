#include "cpmd_global.h"

MODULE forcedr_driver
  USE bsym,                            ONLY: bsclcs
  USE cnstfc_utils,                    ONLY: restfc
  USE cotr,                            ONLY: cotr007,&
                                             gsrate,&
                                             resval,&
                                             resval_dest
  USE error_handling,                  ONLY: stopgm
  USE forces_driver,                   ONLY: forces
  USE kinds,                           ONLY: real_8
  USE mm_dim_utils,                    ONLY: mm_dim
  USE mm_dimmod,                       ONLY: mm_go_mm,&
                                             mm_go_qm,&
                                             mm_revert
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: ipcurr,&
                                             np_local,&
                                             np_low
  USE reshaper,                        ONLY: reshape_inplace
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt_ions
  USE utils,                           ONLY: unitmx
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: forcedr

CONTAINS

  ! ==================================================================
  SUBROUTINE forcedr(c0,c2,sc0,rhoe,psi,tau0,fion,eigv,&
       nstate,nkpoint,lproj,tfor,tcsize)
    ! ==--------------------------------------------------------------==
    ! ==             DRIVER FOR FORCE CALCULATION                     ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(IN),TARGET&
         __CONTIGUOUS                        :: c0(:,:)
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: c2(:,:)
    REAL(real_8),INTENT(OUT) __CONTIGUOUS    :: rhoe(:,:)
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: psi(:,:)
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: tau0(:,:,:)
    REAL(real_8),INTENT(OUT) __CONTIGUOUS    :: fion(:,:,:)
    !eigv not used
    REAL(real_8),INTENT(INOUT)               :: eigv(*)
    INTEGER,INTENT(IN)                       :: nstate
    ! sc0 not used
    COMPLEX(real_8),INTENT(INOUT)            :: sc0(ncpw%ngw,nstate)
    INTEGER,INTENT(IN)                       :: nkpoint
    LOGICAL,INTENT(IN)                       :: lproj, tfor, tcsize

    CHARACTER(*), PARAMETER                  :: procedureN = 'forcedr'

    COMPLEX(real_8), POINTER __CONTIGUOUS    :: c0_ptr(:,:,:)
    INTEGER                                  :: i, ierr, ipx, isub, numx
    LOGICAL                                  :: oldstatus, statusdummy

    CALL tiset(procedureN,isub)
    CALL mm_dim(mm_go_mm,oldstatus)
    CALL zeroing(fion)!,3*maxsys%nax*maxsys%nsx)
    CALL mm_dim(mm_go_qm,statusdummy)
    NULLIFY(c0_ptr)
    CALL reshape_inplace(c0, (/SIZE(c0,1),SIZE(c0,2),1/), c0_ptr)
    CALL forces(c0_ptr,c2,tau0,fion,rhoe,psi,&
         nstate,nkpoint,lproj,tfor,tcsize)
    ! Forces from restraints (may be in the MM part as well)
    IF ( tfor ) THEN
       CALL mm_dim(mm_go_mm,statusdummy)
       IF (paral%parent) THEN
          DO i=1,cotr007%mrestr
             resval(i)=resval(i)+gsrate(i)*dt_ions
             IF (resval_dest(i).NE.-999._real_8) THEN
                IF (gsrate(i).GT.0._real_8.AND.resval(i).GT.resval_dest(i))&
                     resval(i)=resval_dest(i) ! increase
                IF (gsrate(i).LT.0._real_8.AND.resval(i).LT.resval_dest(i))&
                     resval(i)=resval_dest(i) ! decrease
             ENDIF
          ENDDO
          CALL restfc(tau0,fion)
       ENDIF
       !CALL mp_bcast(fion,3*maxsys%nax*maxsys%nsx,parai%source,parai%allgrp)
       CALL mp_bcast(fion,3*maxsys%nax*maxsys%nsx,parai%io_source,parai%cp_grp)
    ENDIF
    CALL mm_dim(mm_revert,oldstatus)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE forcedr

END MODULE forcedr_driver
