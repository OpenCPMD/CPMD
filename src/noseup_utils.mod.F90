#include "cpmd_global.h"

MODULE noseup_utils
  USE bsym,                            ONLY: bsclcs
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes
  USE enosmove_utils,                  ONLY: enosmove
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: real_8
  USE mp_interface,                    ONLY: mp_bcast
  USE nose,                            ONLY: dtsuz,&
                                             ncalls,&
                                             nit,&
                                             tnosepc
  USE parac,                           ONLY: parai,&
                                             paral
  USE pimd,                            ONLY: pimd1,&
                                             pma0s
  USE pnosmove_utils,                  ONLY: pnosmove
  USE prcnosmove_utils,                ONLY: prcnosmove,&
                                             prcnosmove_iso
  USE prcp,                            ONLY: prcpl
  USE prpcmove_utils,                  ONLY: prpcmove,&
                                             prpcmove_iso
  USE prpcnosmove_utils,               ONLY: prpcnosmove,&
                                             prpcnosmove_iso
  USE prpnosmove_utils,                ONLY: prpnosmove,&
                                             prpnosmove_iso
  USE rekine_utils,                    ONLY: rekine
  USE reshaper,                        ONLY: reshape_inplace
  USE rmas,                            ONLY: rmass
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE utils,                           ONLY: print_debug_ions
  USE timer,                           ONLY: tihalt,&
                                             tiset
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: noseup

CONTAINS

  ! ==================================================================
  SUBROUTINE noseup(velp,cm,nstate,ip,use_cp_grps)
    ! ==--------------------------------------------------------------==
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: velp(:,:,:)
    COMPLEX(real_8),INTENT(INOUT)            :: cm(ncpw%ngw,*)
    INTEGER,INTENT(IN)                       :: nstate, ip
    LOGICAL,INTENT(IN),OPTIONAL              :: use_cp_grps
    CHARACTER(*), PARAMETER                  :: procedureN = 'noseup'

    LOGICAL                                  :: cp_active
    INTEGER                                  :: i, isub, j, &
                                                ig, ibeg_c0, iend_c0
    REAL(real_8)                             :: ekinc, sctot
    REAL(real_8),POINTER __CONTIGUOUS        :: cm_r(:,:)

    CALL tiset(procedureN,isub)
    IF(PRESENT(use_cp_grps))THEN
       cp_active=use_cp_grps
    ELSE
       cp_active=.FALSE.
    END IF
    IF(cp_active)THEN
       CALL cp_grp_get_sizes(first_g=ibeg_c0,last_g=iend_c0)
    ELSE
       ibeg_c0=1
       iend_c0=ncpw%ngw
    END IF
    IF(cntl%tverbosevel)THEN
       CALL print_debug_ions('DEBUG VELOCITIES '//procedureN, velp)
    END IF

    IF (cntl%tnosee) THEN
       CALL rekine(cm,nstate,ekinc,use_cp_grps=cp_active)
       IF (paral%io_parent) THEN
          sctot=1.0_real_8
          DO i=1,nit
             DO j=1,ncalls
                CALL enosmove(ekinc,dtsuz(j),sctot,ip)
             ENDDO
          ENDDO
       ENDIF
       CALL mp_bcast(sctot,parai%io_source,parai%cp_grp)
       CALL reshape_inplace(cm,(/2*ncpw%ngw,nstate/),cm_r)
       CALL update_cm(cm_r,sctot,nstate,ibeg_c0*2-1,iend_c0*2)
    ENDIF
    ! FOR BS_CPMD: IF WF IS HS, SKIP THE REST
    IF (cntl%bsymm.AND.(bsclcs.EQ.2))THEN
       CALL tihalt(procedureN,isub)
       RETURN
    ENDIF
    !
    IF (cntl%tnosep.AND.paral%io_parent.AND..NOT.cntl%tnosec) THEN
       DO i=1,nit
          DO j=1,ncalls
             IF (cntl%tpath.AND.cntl%tpimd) THEN
                IF (.NOT.((pimd1%tcentro.OR.pimd1%tringp).AND.ip.EQ.1.AND..NOT.tnosepc)) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL pnosmove(velp,dtsuz(j),pma0s(1,ip),ip)
                   ELSE
                      CALL pnosmove(velp,dtsuz(j),rmass%pma,ip)
                   ENDIF
                ENDIF
             ELSE
                CALL pnosmove(velp,dtsuz(j),rmass%pma,ip)
             ENDIF
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! AK: FIXME  
    IF (cntl%tprcp.AND.paral%io_parent.AND.prcpl%tzflex) THEN
       CALL stopgm(procedureN,'NO NOSE THERMOSTAT FOR ZFLEXIBLE CELL',& 
            __LINE__,__FILE__)
    ENDIF
    ! 
    IF (cntl%tprcp.AND.ip.EQ.1.AND.paral%io_parent.AND.prcpl%tisot) THEN
       IF ((cntl%tnosep.AND.tnosepc).AND..NOT.cntl%tnosec) THEN
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prpnosmove_iso(velp,dtsuz(j),pma0s(1,ip),ip)
                   ELSE
                      CALL prpnosmove_iso(velp,dtsuz(j),rmass%pma,ip)
                   ENDIF
                ELSE
                   CALL prpnosmove_iso(velp,dtsuz(j),rmass%pma,ip)
                ENDIF
             ENDDO
          ENDDO
       ELSEIF (.NOT.(cntl%tnosep.AND.tnosepc).AND.cntl%tnosec) THEN
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prcnosmove_iso(velp,dtsuz(j),pma0s(1,ip),ip)
                   ELSE
                      CALL prcnosmove_iso(velp,dtsuz(j),rmass%pma,ip)
                   ENDIF
                ELSE
                   CALL prcnosmove_iso(velp,dtsuz(j),rmass%pma,ip)
                ENDIF
             ENDDO
          ENDDO
       ELSEIF ((cntl%tnosep.AND.tnosepc).AND.cntl%tnosec) THEN
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prpcnosmove_iso(velp,dtsuz(j),pma0s(1,ip),ip)
                   ELSE
                      CALL prpcnosmove_iso(velp,dtsuz(j),rmass%pma,ip)
                   ENDIF
                ELSE
                   CALL prpcnosmove_iso(velp,dtsuz(j),rmass%pma,ip)
                ENDIF
             ENDDO
          ENDDO
       ELSE
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prpcmove_iso(velp,dtsuz(j),pma0s(1,ip))
                   ELSE
                      CALL prpcmove_iso(velp,dtsuz(j),rmass%pma)
                   ENDIF
                ELSE
                   CALL prpcmove_iso(velp,dtsuz(j),rmass%pma)
                ENDIF
             ENDDO
          ENDDO
       ENDIF
       ! ==--------------------------------------------------------------==
    ELSEIF (cntl%tprcp.AND.ip.EQ.1.AND.paral%io_parent.AND..NOT.prcpl%tisot.AND..NOT.cntl%tshock) THEN
       ! ==--------------------------------------------------------------==
       IF ((cntl%tnosep.AND.tnosepc).AND..NOT.cntl%tnosec) THEN
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prpnosmove(velp,dtsuz(j),pma0s(1,ip),ip)
                   ELSE
                      CALL prpnosmove(velp,dtsuz(j),rmass%pma,ip)
                   ENDIF
                ELSE
                   CALL prpnosmove(velp,dtsuz(j),rmass%pma,ip)
                ENDIF
             ENDDO
          ENDDO
       ELSEIF (.NOT.(cntl%tnosep.AND.tnosepc).AND.cntl%tnosec) THEN
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prcnosmove(velp,dtsuz(j),pma0s(1,ip),ip)
                   ELSE
                      CALL prcnosmove(velp,dtsuz(j),rmass%pma,ip)
                   ENDIF
                ELSE
                   CALL prcnosmove(velp,dtsuz(j),rmass%pma,ip)
                ENDIF
             ENDDO
          ENDDO
       ELSEIF ((cntl%tnosep.AND.tnosepc).AND.cntl%tnosec) THEN
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prpcnosmove(velp,dtsuz(j),pma0s(1,ip),ip)
                   ELSE
                      CALL prpcnosmove(velp,dtsuz(j),rmass%pma,ip)
                   ENDIF
                ELSE
                   CALL prpcnosmove(velp,dtsuz(j),rmass%pma,ip)
                ENDIF
             ENDDO
          ENDDO
       ELSE
          DO i=1,nit
             DO j=1,ncalls
                IF (cntl%tpath.AND.cntl%tpimd) THEN
                   IF (pimd1%tpinm.OR.pimd1%tstage) THEN
                      CALL prpcmove(velp,dtsuz(j),pma0s(1,ip))
                   ELSE
                      CALL prpcmove(velp,dtsuz(j),rmass%pma)
                   ENDIF
                ELSE
                   CALL prpcmove(velp,dtsuz(j),rmass%pma)
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    !TK syncronize cp_grps
    IF(paral%parent) CALL mp_bcast(velp,SIZE(velp),parai%io_source,&
         parai%cp_inter_grp)
    IF(cntl%tverbosevel)THEN
       CALL print_debug_ions('DEBUG VELOCITIES '//procedureN, velp)
    END IF
    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE noseup
  ! ==================================================================
  SUBROUTINE update_cm(cm_r,sctot,nstate,ibeg_c0,iend_c0)
    ! ==--------------------------------------------------------------==
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: cm_r(:,:)
    INTEGER,INTENT(IN)                       :: nstate, ibeg_c0, iend_c0
    REAL(real_8),INTENT(IN)                  :: sctot
    INTEGER                                  :: i,ig

    !$omp parallel do private (i,ig)
    DO i=1,nstate
       DO ig=ibeg_c0,iend_c0
          cm_r(ig,i)=cm_r(ig,i)*sctot
       END DO
    END DO

    ! ==--------------------------------------------------------------==
  END SUBROUTINE update_cm
  ! ==================================================================



END MODULE noseup_utils
