#include "cpmd_global.h"

MODULE kin_energy_utils
  USE cppt,                            ONLY: hg
  USE dotp_utils,                      ONLY: dotp
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_c,&
                                             ener_com
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE prcp,                            ONLY: prcp_com
  USE reshaper,                        ONLY: reshape_inplace
  USE special_functions,               ONLY: cp_erf
  USE spin,                            ONLY: clsd,&
                                             lspin2
  USE system,                          ONLY: ncpw,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: kin_energy

CONTAINS

  ! ==================================================================
  SUBROUTINE kin_energy(c0,nstate,rsum)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE KINETIC ENERGY EKIN. IT IS DONE IN RECIPROCAL SPACE     ==
    ! ==  WHERE THE ASSOCIATED OPERATORS ARE DIAGONAL.                ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(IN) __CONTIGUOUS  :: c0(:,:)
    INTEGER,INTENT(IN)                       :: nstate
    REAL(real_8),INTENT(OUT)                 :: rsum

    REAL(real_8), PARAMETER                  :: deltakin = 1.e-10_real_8 
    REAL(real_8),POINTER __CONTIGUOUS        :: c0_r(:,:)
    INTEGER                                  :: ig, isub
    REAL(real_8)                             :: ima, imb, ra, rb, &
                                                sk1, sk2, xkin
    CHARACTER(*), PARAMETER                  :: procedureN = 'kin_energy'
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ! Accumulate the charge and kinetic energy
    rsum=0._real_8
    xkin=0._real_8
    CALL reshape_inplace(c0,(/2*ncpw%ngw,nstate/),c0_r)
    CALL kin_energy_r(c0_r,hg,xkin,rsum,nstate,ncpw%ngw)
    ener_com%ekin=xkin*parm%tpiba2
    ! Other kinetic energies for CAS22 method
    IF (lspin2%tlse .AND. (lspin2%tcas22.OR.lspin2%tpenal)) THEN
       IF (prcp_com%akin.GT.deltakin) CALL stopgm('RHOOFR',&
            'CAS22 and AKIN not implemented',& 
            __LINE__,__FILE__)
       sk1=0._real_8
       DO ig=1,ncpw%ngw
          ra=REAL(c0(ig,clsd%ialpha))
          rb=REAL(c0(ig,clsd%ibeta))
          ima=AIMAG(c0(ig,clsd%ialpha))
          imb=AIMAG(c0(ig,clsd%ibeta))
          sk1=sk1+hg(ig)*(ra*rb+ima*imb)
       ENDDO
       ener_c%ekin_ab = sk1 * parm%tpiba2
    ENDIF
    IF (lspin2%tlse .AND. lspin2%tcas22) THEN
       sk1=0._real_8
       sk2=0._real_8
       DO ig=1,ncpw%ngw
          sk1=sk1+REAL(hg(ig)*CONJG(c0(ig,clsd%ialpha))*c0(ig,clsd%ialpha))
          sk2=sk2+REAL(hg(ig)*CONJG(c0(ig,clsd%ibeta))*c0(ig,clsd%ibeta))
       ENDDO
       ener_c%ekin_a = ener_com%ekin - parm%tpiba2 * ( sk2 - sk1 )
       ener_c%ekin_2 = ener_com%ekin + parm%tpiba2 * ( sk2 - sk1 )
    ENDIF
    ! 
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE kin_energy
  ! ==================================================================
  SUBROUTINE kin_energy_r(c0_r,hg,xkin,rsum,nstate,ngw)
    INTEGER,INTENT(IN)                       :: nstate, ngw
    REAL(real_8),INTENT(IN)                  :: c0_r(2,ngw,nstate),hg(ngw)
    REAL(real_8),INTENT(OUT)                 :: rsum,xkin
    REAL(real_8)                             :: xskin,temp,sk1,arg,g2
    INTEGER                                  :: is1,ig,i
    REAL(real_8), PARAMETER                  :: deltakin = 1.e-10_real_8
    
    rsum=0.0_real_8
    xkin=0.0_real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute &
#else
    !$omp parallel do &
#endif
    !$omp& private(I,SK1,XSKIN,IS1,ARG,G2,IG,TEMP) &
    !$omp  reduction(+:RSUM,XKIN)
    DO i=1,nstate
       IF (crge%f(i,1).NE.0._real_8) THEN! TODO check F(I,1) === F(I)

          sk1=0.0_real_8
          temp=0._real_8
          is1=1
          IF (prcp_com%akin.GT.deltakin) THEN
             xskin=1._real_8/prcp_com%gskin
             IF (geq0) THEN
                is1=2
                arg=-prcp_com%gckin*xskin
                g2=0.5_real_8*prcp_com%gakin*(1._real_8+cp_erf(arg))
                sk1=sk1+g2*(c0_r(1,1,i)**2+c0_r(2,1,i)**2)
                temp=c0_r(1,1,i)**2*0.5_real_8
             ENDIF
#if defined(_HAS_OMP_TARGET_OFFLOAD)             
             !$omp parallel do reduction(+:temp,sk1)
#endif
             DO ig=is1,ngw
                arg=(hg(ig)-prcp_com%gckin)*xskin
                g2=hg(ig)+prcp_com%gakin*(1._real_8+cp_erf(arg))
                sk1=sk1+g2*(c0_r(1,ig,i)**2+c0_r(2,ig,i)**2)
                temp=temp+c0_r(1,ig,i)**2+c0_r(2,ig,i)**2
             ENDDO
          ELSE
             IF(geq0)THEN
                temp=c0_r(1,1,i)**2*0.5_real_8
                is1=2
             END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)             
             !$omp parallel do reduction(+:temp,sk1)
#endif
             DO ig=is1,ngw
                sk1=sk1+hg(ig)*(c0_r(1,ig,i)**2+c0_r(2,ig,i)**2)
                temp=temp+c0_r(1,ig,i)**2+c0_r(2,ig,i)**2
             ENDDO
          ENDIF
          rsum=rsum+temp*crge%f(i,1)*2.0_real_8
          xkin=xkin+crge%f(i,1)*sk1
       ENDIF
    ENDDO
  END SUBROUTINE kin_energy_r
  
END MODULE kin_energy_utils
