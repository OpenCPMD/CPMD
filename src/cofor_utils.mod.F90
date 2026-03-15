#include "cpmd_global.h"

MODULE cofor_utils
  USE cppt,                            ONLY: gk,&
                                             inyh
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes
  USE ions,                            ONLY: ions1,&
                                             ions0
  USE mp_interface,                    ONLY: mp_sum
  USE kinds,                           ONLY: real_8
  USE nlcc,                            ONLY: corel,&
                                             rhoc,&
                                             vnlcc,&
                                             vnlt
  USE parac,                           ONLY: parai
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             ncpw,&
                                             parm,&
                                             maxsys
  USE utils,                           ONLY: print_debug_ions
  USE timer,                           ONLY: tihalt,&
                                             tiset


  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cofor

CONTAINS

  ! ==================================================================
  SUBROUTINE cofor(fion,vpot)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == FORCES ON IONS DUE TO CORE CHARGES                           ==
    ! ==--------------------------------------------------------------==
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    COMPLEX(real_8)                          :: vpot(ncpw%nhg,*)

    COMPLEX(real_8)                          :: ei123, vxc
    INTEGER                                  :: ia, ig, is, isa, k, isa0, &
                                                ig_start, ig_end, isub
    REAL(real_8)                             :: omtp, vcgs, ft1, ft2, ft3
    CHARACTER(*), PARAMETER                  :: procedureN = 'cofor'

    CALL tiset(procedureN,isub)

    IF(parai%cp_nogrp.GT.1.AND.parai%cp_inter_me.GT.0)THEN
       !$omp parallel do private(is,ia)
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             fion(:,ia,is)=0._real_8
          END DO
       END DO
    END IF
    CALL cp_grp_get_sizes(first_nhg=ig_start,last_nhg=ig_end)
    omtp=2._real_8*parm%omega*parm%tpiba
    IF(cntl%bigmem)THEN
       IF(cntl%tlsd)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target teams private(isa0,is,ft1,ft2,ft3,ia,isa)
#else
          !$omp parallel private(isa0,is,ia,ft1,ft2,ft3,isa,ig,vxc,vcgs,k)
#endif
          isa0=0
          DO is=1,ions1%nsp
             IF(corel%tnlcc(is))THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
                !$omp distribute
#else
                !$omp do
#endif
                DO ia=1,ions0%na(is)
                   ft1=0._real_8
                   ft2=0._real_8
                   ft3=0._real_8
                   isa=isa0+ia
#if defined(_HAS_OMP_TARGET_OFFLOAD)
                   !$omp parallel do reduction (+:ft1,ft2,ft3) private(vxc,vcgs)
#endif
                   DO ig=ig_start,ig_end
                      vxc=0.5_real_8*(vpot(ig,1)+vpot(ig,2)-2*vnlt(ig)-&
                           vnlcc(ig,1)-vnlcc(ig,2))
                      vcgs=-AIMAG(CONJG(vxc)*eigrb(ig,isa)*rhoc(ig,is))
                      ft1=ft1+gk(1,ig)*vcgs
                      ft2=ft2+gk(2,ig)*vcgs
                      ft3=ft3+gk(3,ig)*vcgs
                   END DO
                   fion(1,ia,is)=fion(1,ia,is)+ft1*omtp
                   fion(2,ia,is)=fion(2,ia,is)+ft2*omtp
                   fion(3,ia,is)=fion(3,ia,is)+ft3*omtp
                END DO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
                !$omp end do nowait
#endif
             END IF
             isa0=isa0+ions0%na(is)
          END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp end target teams
#else
          !$omp end parallel
#endif
       ELSE
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target teams private(isa0,is,ft1,ft2,ft3,ia,isa)
#else
          !$omp parallel private(isa0,is,ia,ft1,ft2,ft3,isa,ig,vxc,vcgs,k)
#endif
          isa0=0
          DO is=1,ions1%nsp
             IF(corel%tnlcc(is))THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
                !$omp distribute
#else
                !$omp do
#endif
                DO ia=1,ions0%na(is)
                   ft1=0._real_8
                   ft2=0._real_8
                   ft3=0._real_8
                   isa=isa0+ia
#if defined(_HAS_OMP_TARGET_OFFLOAD)
                   !$omp parallel do reduction (+:ft1,ft2,ft3) private(vxc,vcgs)
#endif
                   DO ig=ig_start,ig_end
                      vxc=vpot(ig,1)-vnlt(ig)-vnlcc(ig,1)
                      vcgs=-AIMAG(CONJG(vxc)*eigrb(ig,isa)*rhoc(ig,is))
                      ft1=ft1+gk(1,ig)*vcgs
                      ft2=ft2+gk(2,ig)*vcgs
                      ft3=ft3+gk(3,ig)*vcgs
                   END DO
                   fion(1,ia,is)=fion(1,ia,is)+ft1*omtp
                   fion(2,ia,is)=fion(2,ia,is)+ft2*omtp
                   fion(3,ia,is)=fion(3,ia,is)+ft3*omtp
                END DO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
                !$omp end do nowait
#endif
             END IF
             isa0=isa0+ions0%na(is)
          END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp end target teams
#else
          !$omp end parallel
#endif
       END IF
    ELSE
       !$omp parallel do private(isa,ia,is,ig,ei123,vxc,vcgs)
       DO isa=1,ions1%nat
          ia=iatpt(1,isa)
          is=iatpt(2,isa)
          IF (corel%tnlcc(is)) THEN
             DO ig=ig_start,ig_end
                ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                     ei3(isa,inyh(3,ig))
                vxc=vpot(ig,1)-vnlt(ig)-vnlcc(ig,1)
                IF (cntl%tlsd) vxc=0.5_real_8*(vxc+vpot(ig,2)-vnlt(ig)-vnlcc(ig,2))
                vcgs=-AIMAG(CONJG(vxc)*ei123*rhoc(ig,is))
                fion(1,ia,is)=fion(1,ia,is)+gk(1,ig)*vcgs*omtp
                fion(2,ia,is)=fion(2,ia,is)+gk(2,ig)*vcgs*omtp
                fion(3,ia,is)=fion(3,ia,is)+gk(3,ig)*vcgs*omtp
             ENDDO
          ENDIF
       ENDDO
    END IF

    IF(parai%cp_nogrp.GT.1)THEN
       CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%cp_inter_grp)
    END IF

    IF(cntl%tverbosefor)THEN
       CALL print_debug_ions('DEBUG FORCES '//procedureN, fion)
    END IF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE cofor
  ! ==================================================================

END MODULE cofor_utils
