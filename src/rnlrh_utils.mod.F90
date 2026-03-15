#include "cpmd_global.h"

MODULE rnlrh_utils
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms
  USE cvan,                            ONLY: dvan
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_c,&
                                             ener_com
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE kpnt,                            ONLY: wk
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlm,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: fnl,&
                                             fnl2,&
                                             fnlgp,&
                                             fnl_packed
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE spin,                            ONLY: clsd,&
                                             lspin2
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             ipept,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlrh
  !public :: rnlcas
  PUBLIC :: rnlrhg

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlrh(enl,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==     THE NON-LOCAL POTENTIAL CONTRIBUTION TO THE TOTAL        ==
    ! ==     ENERGY, I.E. ENL                                         ==
    ! ==     K-POINTS IMPLEMENTED (FNL IS COMPLEX -> IMAGP=2)         ==
    ! ==--------------------------------------------------------------==
    REAL(real_8), INTENT(OUT)                :: enl
    INTEGER, INTENT(IN)                      :: nstate, nkpoint

    INTEGER                                  :: i, ia, ii, ik, is, isa, isa0, &
                                                isub, iv, jv, ki, kj, l, l2, &
                                                li, lj ,ia_sum, na(2,ions1%nsp), &
                                                offset, ierr
    INTEGER, ALLOCATABLE                     :: na_grp(:,:,:)
    REAL(real_8)                             :: sum, weight
    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlrh'

    enl=0._real_8
    ! If no non-local components -> return.
    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    ! == Compute the non-local contribution to the total energy (ENL) ==
    ! ==--------------------------------------------------------------==
    ALLOCATE(na_grp(2,ions1%nsp,0:parai%cp_nogrp-1), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate na_grp',&
         __LINE__,__FILE__)

    CALL cp_grp_split_atoms(na_grp)
    na(:,:)=na_grp(:,:,parai%cp_inter_me)
    IF(pslo_com%tivan)THEN
       ! VANDERBILT PSEUDOPOTENTIAL
       IF (imagp.EQ.2)&
            CALL stopgm(procedureN,'K-POINTS NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       !$omp parallel private(is,i,offset,ia_sum) reduction(+:enl)
       offset=0
       DO is=1,ions1%nsp
          ia_sum=na(2,is)-na(1,is)+1
          IF(ia_sum.EQ.0)CYCLE
          IF (pslo_com%tvan(is)) THEN
             ! VANDERBILT PSEUDOPOTENTIAL
             !$omp do
             DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                enl=enl+calc_enl_uspp(ia_sum,nlps_com%ngh(is),&
                     fnl_packed(offset+1:offset+ia_sum*nlps_com%ngh(is),i),dvan(:,:,is))*&
                     wk(1)*crge%f(i,1)
             END DO
             !$omp end do nowait
          END IF
          offset=offset+ia_sum*nlps_com%ngh(is)
       END DO
       !$omp end parallel
    END IF
    DO ik=1,nkpoint
       isa0=0
       DO is=1,ions1%nsp
          IF(pslo_com%tivan)THEN
          ELSEIF (sgpp1%tsgp(is)) THEN
             ! Stefan Goedecker pp
             DO iv=1,nlps_com%ngh(is)
                l=nghtol(iv,is)+1
                ki=sgpp2%lfval(iv,is)
                li=sgpp2%lpval(iv,is)
                DO jv=iv,nlps_com%ngh(is)
                   l2=nghtol(jv,is)+1
                   lj=sgpp2%lpval(jv,is)
                   IF (l2.EQ.l.AND.li.EQ.lj) THEN
                      kj=sgpp2%lfval(jv,is)
                      sum=0.0_real_8
                      DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                         weight=wk(ik)*crge%f(i,ik)
                         IF (weight.EQ.0._real_8) GOTO 2000
                         ii=i-parap%nst12(parai%mepos,1)+1
                         IF (imagp.EQ.1) THEN
                            IF (cntl%tfdist) THEN
                               DO ia=na(1,is),na(2,is)
                                  isa=isa0+ia
                                  sum=sum+weight*fnl2(1,isa,iv,ii,ik)*&
                                       fnl2(1,isa,jv,ii,ik)
                               ENDDO
                            ELSE
                               DO ia=na(1,is),na(2,is)
                                  isa=isa0+ia
                                  sum=sum+weight*fnl2(1,isa,iv,i,ik)*&
                                       fnl2(1,isa,jv,i,ik)
                               ENDDO
                            ENDIF
                         ELSE
                            IF (cntl%tfdist) THEN
                               DO ia=na(1,is),na(2,is)
                                  isa=isa0+ia
                                  sum=sum+weight*&
                                       (fnl2(1,isa,iv,ii,ik)*fnl2(1,isa,jv,ii,ik)&
                                       +fnl2(2,isa,iv,ii,ik)*fnl2(2,isa,jv,ii,ik))
                               ENDDO
                            ELSE
                               DO ia=na(1,is),na(2,is)
                                  isa=isa0+ia
                                  sum=sum+weight*&
                                       (fnl2(1,isa,iv,i,ik)*fnl2(1,isa,jv,i,ik)&
                                       +fnl2(2,isa,iv,i,ik)*fnl2(2,isa,jv,i,ik))
                               ENDDO
                            ENDIF
                         ENDIF
2000                     CONTINUE
                      ENDDO
                      IF (iv.NE.jv) sum=2._real_8*sum
                      enl=enl+sum*sgpp2%hlsg(ki,kj,l,is)
                   ENDIF
                ENDDO
             ENDDO
          ELSE
             ! BHS AND RELATED 
             DO iv=1,nlps_com%ngh(is)
                sum=0.0_real_8
                DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                   ii=i-parap%nst12(parai%mepos,1)+1
                   IF (imagp.EQ.1) THEN
                      IF (cntl%tfdist) THEN
                         DO ia=na(1,is),na(2,is)
                            isa=isa0+ia
                            sum=sum+crge%f(i,ik)*wk(ik)*&
                                 fnl2(1,isa,iv,ii,ik)*fnl2(1,isa,iv,ii,ik)
                         ENDDO
                      ELSE
                         DO ia=na(1,is),na(2,is)
                            isa=isa0+ia
                            sum=sum+crge%f(i,ik)*wk(ik)*&
                                 fnl2(1,isa,iv,i,ik)*fnl2(1,isa,iv,i,ik)
                         ENDDO
                      ENDIF
                   ELSE
                      IF (cntl%tfdist) THEN
                         DO ia=na(1,is),na(2,is)
                            isa=isa0+ia
                            sum=sum+crge%f(i,ik)*wk(ik)*&
                                 (fnl2(1,isa,iv,ii,ik)*fnl2(1,isa,iv,ii,ik)&
                                 +fnl2(2,isa,iv,ii,ik)*fnl2(2,isa,iv,ii,ik))
                         ENDDO
                      ELSE
                         DO ia=na(1,is),na(2,is)
                            isa=isa0+ia
                            sum=sum+crge%f(i,ik)*wk(ik)*&
                                 (fnl2(1,isa,iv,i,ik)*fnl2(1,isa,iv,i,ik)&
                                 +fnl2(2,isa,iv,i,ik)*fnl2(2,isa,iv,i,ik))
                         ENDDO
                      ENDIF
                   ENDIF
                ENDDO
                enl=enl+wsg(is,iv)*sum
             ENDDO
          ENDIF
          isa0=isa0+ions0%na(is)
       ENDDO
    ENDDO
    IF (lspin2%tlse .AND. lspin2%tcas22) CALL rnlcas
    IF(parai%cp_nogrp.GT.1) CALL mp_sum(enl,parai%cp_inter_grp)
    DEALLOCATE(na_grp, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate na_grp',&
         __LINE__,__FILE__)

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlrh
  ! ==================================================================
  PURE FUNCTION calc_enl_uspp(ia_sum,ngh,fnl_,dvan_) RESULT(enl)
    INTEGER,INTENT(IN)                       :: ia_sum, ngh
    REAL(real_8),INTENT(IN)                  :: fnl_(ia_sum,ngh,*)
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: dvan_(:,:)
    REAL(real_8)                             :: sum
    REAL(real_8)                             :: enl
    INTEGER                                  :: iv,ia,jv
    enl=0._real_8
    DO iv=1,ngh
       sum=0._real_8
       DO ia=1,ia_sum
          sum=sum+fnl_(ia,iv,1)*fnl_(ia,iv,1)
       END DO
       enl=enl+sum*dvan_(iv,iv)
       DO jv=iv+1,ngh
          sum=0._real_8
          DO ia=1,ia_sum
             sum=sum+fnl_(ia,iv,1)*fnl_(ia,jv,1)
          END DO
          enl=enl+sum*dvan_(jv,iv)*2.0_real_8
       END DO
    END DO
  END FUNCTION calc_enl_uspp
  ! ==================================================================
  SUBROUTINE rnlcas
    ! ==--------------------------------------------------------------==
    ! == Calculates NL PP contribution to CAS22 Energies              ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ia, iat, iatl, is, iv, jv, &
                                                ki, kj, l, l2, li, lj
    REAL(real_8)                             :: sm

! ==--------------------------------------------------------------==

    ener_c%enl_ab = 0._real_8
    ener_c%enl_a  = ener_com%enl
    DO iat=ipept(1,parai%mepos),ipept(2,parai%mepos)
       ia=iatpt(1,iat)
       is=iatpt(2,iat)
       IF (pslo_com%tvan(is)) THEN
          ! VANDERBILT PSEUDOPOTENTIAL
       ELSEIF (sgpp1%tsgp(is)) THEN
          ! Stefan Goedecker pp
          DO iv=1,nlps_com%ngh(is)
             l=nghtol(iv,is)+1
             li=sgpp2%lpval(iv,is)
             ki=sgpp2%lfval(iv,is)
             DO jv=iv,nlps_com%ngh(is)
                IF (iv.EQ.jv) THEN
                   sm=1._real_8
                ELSE
                   sm=2._real_8
                ENDIF
                l2=nghtol(jv,is)+1
                lj=sgpp2%lpval(jv,is)
                IF (l2.EQ.l.AND.li.EQ.lj) THEN
                   kj=sgpp2%lfval(jv,is)
                   IF (cntl%tfdist) THEN
                      iatl=iat-ipept(1,parai%mepos)+1
                   ELSE
                      iatl=iat
                   ENDIF
                   ener_c%enl_ab=ener_c%enl_ab+sgpp2%hlsg(ki,kj,l,is)*0.5_real_8*sm*&
                        (fnl(1,iatl,iv,clsd%ialpha,1)*fnl(1,iatl,jv,clsd%ibeta,1)+&
                        fnl(1,iatl,jv,clsd%ialpha,1)*fnl(1,iatl,iv,clsd%ibeta,1))
                   ener_c%enl_a=ener_c%enl_a+sgpp2%hlsg(ki,kj,l,is)*sm*&
                        (fnl(1,iatl,iv,clsd%ialpha,1)*fnl(1,iatl,jv,clsd%ialpha,1)-&
                        fnl(1,iatl,iv,clsd%ibeta,1)*fnl(1,iatl,jv,clsd%ibeta,1))
                ENDIF
             ENDDO
          ENDDO
       ELSE
          ! BHS AND RELATED
          DO iv=1,nlps_com%ngh(is)
             IF (cntl%tfdist) THEN
                iatl=iat-ipept(1,parai%mepos)+1
             ELSE
                iatl=iat
             ENDIF
             ener_c%enl_ab=ener_c%enl_ab+wsg(is,iv)*fnl(1,iatl,iv,clsd%ialpha,1)*&
                  fnl(1,iatl,iv,clsd%ibeta,1)
             ener_c%enl_a=ener_c%enl_a+wsg(is,iv)*(fnl(1,iatl,iv,clsd%ialpha,1)**2&
                  -fnl(1,iatl,iv,clsd%ibeta,1)**2)
          ENDDO
       ENDIF
    ENDDO
    ener_c%enl_2  = ener_com%enl - (ener_c%enl_a - ener_com%enl)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlcas
  ! EHR[
  ! ==================================================================
  SUBROUTINE rnlrhg(enl,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==     THE NON-LOCAL POTENTIAL CONTRIBUTION TO THE TOTAL        ==
    ! ==     ENERGY, I.E. ENL                                         ==
    ! ==     K-POINTS IMPLEMENTED (FNL IS COMPLEX -> IMAGP=2)         ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: enl
    INTEGER                                  :: nstate, nkpoint

    INTEGER                                  :: i, ia, ii, ik, is, isa, isa0, &
                                                isub, iv
    REAL(real_8)                             :: sum

    enl=0._real_8
    ! If no non-local components -> return.
    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('     RNLRH',isub)
    ! ==--------------------------------------------------------------==
    ! == Compute the non-local contribution to the total energy (ENL) ==
    ! ==--------------------------------------------------------------==
    DO ik=1,nkpoint
       isa0=0
       DO is=1,ions1%nsp
          ! BHS AND RELATED
          DO iv=1,nlps_com%ngh(is)
             sum=0.0_real_8
             DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                ii=i-parap%nst12(parai%mepos,1)+1
                IF (imagp.EQ.1) THEN
                   IF (cntl%tfdist) THEN
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         sum=sum+crge%f(i,ik)*wk(ik)*&
                              fnl2(1,isa,iv,ii,ik)*fnl2(1,isa,iv,ii,ik)
                      ENDDO
                   ELSE
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         sum=sum+crge%f(i,ik)*wk(ik)*&
                              fnl2(1,isa,iv,i,ik)*fnl2(1,isa,iv,i,ik)
                      ENDDO
                   ENDIF
                ELSE
                   IF (cntl%tfdist) THEN
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         sum=sum+crge%f(i,ik)*wk(ik)*&
                              (fnl2(1,isa,iv,ii,ik)*fnl2(1,isa,iv,ii,ik)&
                              +fnl2(2,isa,iv,ii,ik)*fnl2(2,isa,iv,ii,ik))
                      ENDDO
                   ELSE
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         sum=sum+crge%f(i,ik)*wk(ik)* 2._real_8 *&
                              ( (fnl(1,isa,iv,i,ik)*fnlgp(1,isa,iv,i,ik)&
                              +fnl(2,isa,iv,i,ik)*fnlgp(2,isa,iv,i,ik))&
                              )
                      ENDDO
                   ENDIF
                ENDIF
             ENDDO
             enl=enl+wsg(is,iv)*sum
          ENDDO
          isa0=isa0+ions0%na(is)
       ENDDO
    ENDDO
    IF (lspin2%tlse .AND. lspin2%tcas22) CALL rnlcas
    CALL tihalt('     RNLRH',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlrhg
  ! ==================================================================
  ! EHR]

END MODULE rnlrh_utils
