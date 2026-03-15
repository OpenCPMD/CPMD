#include "cpmd_global.h"

MODULE phfac_utils
  USE cppt,                            ONLY: inyh
  USE error_handling,                  ONLY: stopgm
  USE gvec,                            ONLY: gvec_com
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE kpnt,                            ONLY: eigkr,&
                                             eikr,&
                                             rk
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  !$ USE omp_lib,                         ONLY: omp_get_thread_num,&
  !$                                            omp_get_team_num
  USE parac,                           ONLY: paral, &
                                             parai
  USE prmem_utils,                     ONLY: prmem
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigr,&
                                             eigrb,&
                                             eigrb_dipole,&
                                             natx
  USE system,                          ONLY: cntl,&
                                             iatpt,&
                                             kpbeg,&
                                             ncpw,&
                                             nkpbl,&
                                             nkpt,&
                                             parm,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: phfac
  PUBLIC :: phfac_dipole
  PUBLIC :: calc_eigkr

CONTAINS

  SUBROUTINE phfac_dipole(aug_ind,aug_cntrl)
    ! ==--------------------------------------------------------------==
    ! ==  CALCULATION OF ARRAY EIGRB_DIPOLE                           ==
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'phfac_dipole'

    INTEGER, INTENT(IN)          :: aug_ind(6)
    LOGICAL, INTENT(IN)          :: aug_cntrl(6)
    INTEGER                      :: isub, k, ierr
    CALL tiset(procedureN,isub)
    IF(cntl%bigmem)THEN
       IF(.not. allocated(eigrb_dipole))THEN
          ALLOCATE(eigrb_dipole(6,ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       END IF
       eigrb_dipole=cmplx(0.0_real_8,0.0_real_8)
       !$omp parallel do private(k)
       DO k=1, 6
          IF (aug_cntrl(k)) THEN
             eigrb_dipole(k,:)=eigrb(aug_ind(k),:)
          ENDIF
       END DO
       CALL mp_sum(eigrb_dipole,6*ions1%nat,parai%allgrp)
    END IF
    CALL tihalt(procedureN,isub)
    RETURN
    
  END SUBROUTINE phfac_dipole
  
  ! ==================================================================
  SUBROUTINE phfac(tau0,force_update)
    ! ==--------------------------------------------------------------==
    REAL(real_8), INTENT(IN) __CONTIGUOUS    :: tau0(:,:,:)
    LOGICAL, OPTIONAL                        :: force_update

    CHARACTER(*), PARAMETER                  :: procedureN = 'phfac'

    COMPLEX(real_8)                          :: ctem1, ctem2, ctem3, ctep1, &
                                                ctep2, ctep3, ei10, ei20, &
                                                ei30, svtmpm, svtmpp, zsum
    INTEGER                                  :: i, ia, ierr, ig, ik, ikk, &
                                                ikpt, is, isa, isub, j, k, &
                                                nh1, nh2, nh3,  methread, nthreads
    INTEGER(int_8)                           :: il_ei2t(2), il_ei3t(2), il_ei1t(2)
    INTEGER, SAVE                            :: ifirst = 0
    REAL(real_8)                             :: ar1, ar2, ar3, sum, sum1, &
                                                sum2, sum3
#ifdef _USE_SCRATCHLIBRARY
    COMPLEX(real_8),POINTER __CONTIGUOUS     :: ei1t(:,:), ei2t(:,:), ei3t(:,:)
#else
    COMPLEX(real_8),ALLOCATABLE              :: ei1t(:,:), ei2t(:,:), ei3t(:,:)
#endif
! ==--------------------------------------------------------------==

    IF (spar%nr1s.LT.3) THEN
       CALL stopgm('PHFAC',' PHFAC: NR1 TOO SMALL ',& 
            __LINE__,__FILE__)
    ENDIF
    IF (spar%nr2s.LT.3) THEN
       CALL stopgm('PHFAC',' PHFAC: NR2 TOO SMALL ',& 
            __LINE__,__FILE__)
    ENDIF
    IF (spar%nr3s.LT.3) THEN
       CALL stopgm('PHFAC',' PHFAC: NR3 TOO SMALL ',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    IF (ifirst.EQ.0) THEN
       !TK either eigrb or ei1-3
       IF (cntl%bigmem) THEN
          ALLOCATE(eigrb(ncpw%nhg,ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target enter data map(alloc:eigrb)
#endif
       ENDIF
       ALLOCATE(ei1(natx,(2*spar%nr1s-1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ei2(natx,(2*spar%nr2s-1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE(ei3(natx,(2*spar%nr3s-1)),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

       ALLOCATE(eigr(ncpw%ngw,ions1%nat,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)! FIXME deallocate missing
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target enter data map(alloc:eigr,ei1,ei2,ei3)
#endif
       IF (tkpts%tkpnt) THEN
          ALLOCATE(eikr(nkpt%nkpts,ions1%nat),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)

          ALLOCATE(eigkr(nkpt%ngwk,ions1%nat,nkpt%nkpnt),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
       ifirst  = 1
       IF (paral%parent) THEN
          CALL prmem('     PHFAC')
       ENDIF
    ENDIF
    !TK temporary arrays
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    nthreads=ions1%nat
#else
    nthreads=parai%ncpus
#endif
    il_ei1t(1)=2*spar%nr1s-1
    il_ei1t(2)=nthreads
    il_ei2t(1)=2*spar%nr2s-1
    il_ei2t(2)=nthreads
    il_ei3t(1)=2*spar%nr3s-1
    il_ei3t(2)=nthreads
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_ei1t,ei1t,procedureN//'_ei1t',ierr)
#else
    ALLOCATE(ei1t(il_ei1t(1),il_ei2t(2)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_ei2t,ei2t,procedureN//'_ei2t',ierr)
#else
    ALLOCATE(ei2t(il_ei2t(1),il_ei2t(2)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_ei3t,ei3t,procedureN//'_ei3t',ierr)
#else
    ALLOCATE(ei3t(il_ei3t(1),il_ei3t(2)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    nh1=spar%nr1s/2
    nh2=spar%nr2s/2
    nh3=spar%nr3s/2
    methread=1
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams num_teams(ions1%nat) &
    !$omp& map (to:tau0,gvec_com,iatpt,parm) &
#else
    !$omp parallel &
#endif
    !$omp& private(isa,ia,is,sum1,sum2,sum3,ar1,ar2,ar3) &
    !$omp private(ctep1,ctep2,ctep3,ctem1,ctem2,ctem3) &
    !$omp private(svtmpp,svtmpm,i,j,k,ei10,ei20,ei30,methread) &
    !$omp shared(nh1,nh2,nh3)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$ methread=omp_get_team_num()+1
    !$omp distribute
#else
    !$ methread=omp_get_thread_num()+1
    !$omp do
#endif
    DO isa=1,ions1%nat
       ia=iatpt(1,isa)
       is=iatpt(2,isa)
       sum1=gvec_com%b1(1)*tau0(1,ia,is)+gvec_com%b1(2)*tau0(2,ia,is)+&
            gvec_com%b1(3)*tau0(3,ia,is)
       sum2=gvec_com%b2(1)*tau0(1,ia,is)+gvec_com%b2(2)*tau0(2,ia,is)+&
            gvec_com%b2(3)*tau0(3,ia,is)
       sum3=gvec_com%b3(1)*tau0(1,ia,is)+gvec_com%b3(2)*tau0(2,ia,is)+&
            gvec_com%b3(3)*tau0(3,ia,is)
       ar1=parm%tpiba*sum1
       ar2=parm%tpiba*sum2
       ar3=parm%tpiba*sum3
       ctep1=CMPLX(COS(ar1),-SIN(ar1),kind=real_8)
       ctep2=CMPLX(COS(ar2),-SIN(ar2),kind=real_8)
       ctep3=CMPLX(COS(ar3),-SIN(ar3),kind=real_8)
       ctem1=CONJG(ctep1)
       ctem2=CONJG(ctep2)
       ctem3=CONJG(ctep3)

       ei10=ctep1**(-nh1)
       ei1t(1,methread)=CMPLX(1.0_real_8,0.0_real_8,kind=real_8)*ei10
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do private(i)
#endif     
       DO i=2,spar%nr1s
          ei1t(i,methread)=ctep1**(i-1)*ei10
          ei1t(spar%nr1s+i-1,methread)=ctem1**(i-1)*ei10
       END DO

       ei20=ctep2**(-nh2)
       ei2t(1,methread)=CMPLX(1.0_real_8,0.0_real_8,kind=real_8)*ei20
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do private(j)
#endif
       DO j=2,spar%nr2s
          ei2t(j,methread)=ctep2**(j-1)*ei20
          ei2t(spar%nr2s+j-1,methread)=ctem2**(j-1)*ei20
       END DO

       ei30=ctep3**(-nh3)
       ei3t(1,methread)=CMPLX(1.0_real_8,0.0_real_8,kind=real_8)*ei30
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do private(k)
#endif
       DO k=2,spar%nr3s
          ei3t(k,methread)=ctep3**(k-1)*ei30
          ei3t(spar%nr3s+k-1,methread)=ctem3**(k-1)*ei30
       END DO

       ! ==--------------------------------------------------------------==
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do private(ig)
#endif
       DO ig=1,2*spar%nr1s-1
          ei1(isa,ig)=ei1t(ig,methread)
       END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do private(ig)
#endif
       DO ig=1,2*spar%nr2s-1
          ei2(isa,ig)=ei2t(ig,methread)
       END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do private(ig)
#endif
       DO ig=1,2*spar%nr3s-1
          ei3(isa,ig)=ei3t(ig,methread)
       END DO

       IF (cntl%bigmem) THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp parallel do private(ig)
#endif
          DO ig=1,ncpw%ngw
             eigr(ig,isa,1)=ei1t(inyh(1,ig),methread)*ei2t(inyh(2,ig),methread)&
                  *ei3t(inyh(3,ig),methread)
             eigrb(ig,isa)=ei1t(inyh(1,ig),methread)*ei2t(inyh(2,ig),methread)&
                  *ei3t(inyh(3,ig),methread)
          END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp parallel do private(ig)
#endif
          DO ig=ncpw%ngw+1,ncpw%nhg
             eigrb(ig,isa)=ei1t(inyh(1,ig),methread)*ei2t(inyh(2,ig),methread)&
                  *ei3t(inyh(3,ig),methread)
          END DO
       ELSE
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp parallel do private(ig)
#endif
          DO ig=1,ncpw%ngw
             eigr(ig,isa,1)=ei1t(inyh(1,ig),methread)*ei2t(inyh(2,ig),methread)&
                  *ei3t(inyh(3,ig),methread)
          END DO
       END IF
    END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end target teams
    IF(PRESENT(force_update))THEN
       !$omp target update from(eigr,eigrb,ei1,ei2,ei3) IF(force_update)
    END IF
#else
    !$omp end do nowait
    !$omp end parallel
#endif
    ! ==--------------------------------------------------------------==
    IF (tkpts%tkpnt) THEN
       DO ikpt=1,nkpt%nblkp
         DO ik=1,nkpbl(ikpt)
             ikk=kpbeg(ikpt)+ik
             !$omp parallel do private(ISA,IA,IS,SUM1,SUM2,SUM3,SUM,ZSUM,IG)
             DO isa=1,ions1%nat
                ia=iatpt(1,isa)
                is=iatpt(2,isa)
                sum1=gvec_com%b1(1)*tau0(1,ia,is)+gvec_com%b1(2)*tau0(2,ia,is)+&
                     gvec_com%b1(3)*tau0(3,ia,is)
                sum2=gvec_com%b2(1)*tau0(1,ia,is)+gvec_com%b2(2)*tau0(2,ia,is)+&
                     gvec_com%b2(3)*tau0(3,ia,is)
                sum3=gvec_com%b3(1)*tau0(1,ia,is)+gvec_com%b3(2)*tau0(2,ia,is)+&
                     gvec_com%b3(3)*tau0(3,ia,is)
                sum=rk(1,ikk)*sum1+rk(2,ikk)*sum2+rk(3,ikk)*sum3
                zsum=CMPLX(COS(sum),SIN(sum),kind=real_8)
                eikr(ikk,isa)=zsum
                DO ig=1,ncpw%ngw
                   eigkr(ig,isa,ik)=eigr(ig,isa,1)*zsum
                   eigkr(ig+ncpw%ngw,isa,ik)=CONJG(eigr(ig,isa,1))*zsum
                END DO
             END DO
          END DO
          IF (tkpts%tkblock) THEN
             CALL wkpt_swap(eigkr,1,ikpt,'EIGKR')
          ENDIF
       END DO
    ELSE
       ! kpt..not to duplicate code.
       ! NB: it is impossible to achieve this using Fortran 90, because
       ! EIGKR has rank 3, and EIGR has rank 2. Even if 
       ! EIGKR has dimensions (:,:,1) it cannot point to (:,:)
       ! To refactor this we will have to either make EIGR rank 3 with last dim = 1,
       ! or allocate memory for EIGKR and copy...
       ! 
       ! semantically this should be the following
       ! EIGKR => EIGR
       ! in principle, EIGKR should be just removed (?), and EIGR should have rank 3
       ! where last dim is NKPNT

       ! TODO: fix this!
       eigkr => eigr
    ENDIF

#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_ei3t,ei3t,procedureN//'_ei3t',ierr)
#else
    DEALLOCATE(ei3t,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_ei2t,ei2t,procedureN//'_ei2t',ierr)
#else
    DEALLOCATE(ei2t,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_ei1t,ei1t,procedureN//'_ei1t',ierr)
#else
    DEALLOCATE(ei1t,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE phfac
  ! ==================================================================
  SUBROUTINE calc_eigkr(ikpt)
    ! ==--------------------------------------------------------------==
    ! == CALCULATES EIGKR (USED WITH THE OPTION TKBLOCK)              ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ikpt

    COMPLEX(real_8)                          :: zsum
    INTEGER                                  :: ia, ig, ik, ikind, is, isa, &
                                                isub

    CALL tiset('CALC_EIGKR',isub)
    DO ikind=1,nkpbl(ikpt)
       ik=kpbeg(ikpt)+ikind
       isa=0
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             isa=isa+1
             zsum=eikr(ik,isa)
             DO ig=1,ncpw%ngw
                eigkr(ig,isa,ikind)=eigr(ig,isa,1)*zsum
                eigkr(ig+ncpw%ngw,isa,ikind)=CONJG(eigr(ig,isa,1))*zsum
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL tihalt('CALC_EIGKR',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE calc_eigkr
  ! ==================================================================

END MODULE phfac_utils
