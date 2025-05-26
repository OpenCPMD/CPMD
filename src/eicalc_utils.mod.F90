#include "cpmd_global.h"

MODULE eicalc_utils
  USE cppt,                            ONLY: gk,&
                                             inyh,&
                                             rhops,&
                                             vps
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nvtx_utils
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             maxsys,&
                                             ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: eicalc
  PUBLIC :: eicalc1

CONTAINS

  ! ==================================================================
  SUBROUTINE eicalc(eivps,eirop,update_to_host)
    ! ==--------------------------------------------------------------==
    ! == EIVPS : phase factor times local pseudopotential  (VPS)      ==
    ! == EIROP : phase factor times Gaussian charge distributions     ==
    ! ==         which replaced ionic point charges (RHOPS)           ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: eivps(:), eirop(:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'eicalc'

    COMPLEX(real_8)                          :: ei123, eivps_s, eirop_s
    INTEGER                                  :: ia, ig, is, isa, isub, isa0,&
                                                ind1, ind2, ind3
    REAL(real_8)                             :: ei, er,  vps_s, rhops_s
    LOGICAL, OPTIONAL, INTENT(IN)            :: update_to_host
    LOGICAL                                  :: update
    
    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )
    IF(PRESENT(update_to_host))THEN
       update=update_to_host
    ELSE
       update=.TRUE.
    END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute num_teams(ncpw%nhg) &
#else
    !$omp parallel do &
#endif
    !$omp& private(IG,ISA,IA,IS,ER,EI,EI123,ISA0,eivps_s,eirop_s,&
    !$omp& ind1,ind2,ind3,vps_s,rhops_s) shared(EIVPS,EIROP)
    DO ig=1,ncpw%nhg
       eivps_s=(0.0_real_8,0.0_real_8)
       eirop_s=(0.0_real_8,0.0_real_8)
       ind1=inyh(1,ig)
       ind2=inyh(2,ig)
       ind3=inyh(3,ig)
       isa0=0
       DO is=1,ions1%nsp
          rhops_s=rhops(is,ig)
          vps_s=vps(is,ig)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp parallel do reduction(+:eivps_s,eirop_s), &
          !$omp& private(ia,isa,ei123,er,ei) num_threads(16)
#endif
          DO ia=1,ions0%na(is)
             isa=isa0+ia
             ei123=ei1(isa,ind1)*ei2(isa,ind2)*&
                  ei3(isa,ind3)
             er=REAL(ei123)
             ei=AIMAG(ei123)
             eivps_s=eivps_s+&
                  CMPLX(er*vps_s,ei*vps_s,kind=real_8)
             eirop_s=eirop_s+&
                  CMPLX(er*rhops_s,ei*rhops_s,kind=real_8)
          END DO
          isa0=isa0+ions0%na(is)
       END DO
       eivps(ig)=eivps_s
       eirop(ig)=eirop_s
    END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target update from(eivps,eirop) IF(update)
#endif

    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE eicalc
  ! ==================================================================
  SUBROUTINE eicalc1(k,is,isa,eivps1,eirop1)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: k, is, isa
    COMPLEX(real_8)                          :: eivps1(ncpw%nhg), &
                                                eirop1(ncpw%nhg)

    COMPLEX(real_8)                          :: ei123, g
    INTEGER                                  :: ig, isub
    REAL(real_8)                             :: ei, er

    CALL tiset('   EICALC1',isub)
    IF (cntl%bigmem) THEN
       !$omp parallel do private(IG,G)
       DO ig=1,ncpw%nhg
          g=CMPLX(0._real_8,-gk(k,ig),kind=real_8)*parm%tpiba*eigrb(ig,isa)
          eivps1(ig)=vps(is,ig)*g
          eirop1(ig)=rhops(is,ig)*g
       ENDDO
    ELSE
       !$omp parallel do private(IG,G,EI123,ER,EI)
       DO ig=1,ncpw%nhg
          g=CMPLX(0._real_8,-gk(k,ig),kind=real_8)*parm%tpiba
          ei123=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
               ei3(isa,inyh(3,ig))
          er=REAL(ei123)
          ei=AIMAG(ei123)
          eivps1(ig)=CMPLX(er*vps(is,ig),ei*vps(is,ig),kind=real_8)*g
          eirop1(ig)=CMPLX(er*rhops(is,ig),ei*rhops(is,ig),kind=real_8)*g
       ENDDO
    ENDIF
    CALL tihalt('   EICALC1',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE eicalc1
  ! ==================================================================

END MODULE eicalc_utils
