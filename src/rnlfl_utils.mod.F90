#include "cpmd_global.h"

MODULE rnlfl_utils
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms  
  USE cvan,                            ONLY: qq
  USE cvan,                            ONLY: deeq,&
                                             dvan
  USE distribution_utils,              ONLY: dist_entity
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8,&
       int_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nlps_com
  !$ USE omp_lib,                       ONLY: omp_get_thread_num,&
  !$                                          omp_get_team_num
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: il_fnl_packed,&
                                             il_dfnl_packed
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             maxsys
  USE utils,                           ONLY: print_debug_ions
  USE timer,                           ONLY: tihalt,&
                                             tiset
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlfl

CONTAINS

  ! ==================================================================
  SUBROUTINE rnlfl(fion,nstate,nkpoint,fnl_p,fnlgam_p,dfnl_p)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE FORCE ON THE IONIC DEGREES OF FREEDOM DUE TO THE        ==
    ! ==  ORTHOGONALITY CONSTRAINED IN VANDERBILT PSEUDO-POTENTIALS   ==
    ! ==--------------------------------------------------------------==
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    INTEGER,INTENT(IN)                       :: nstate, nkpoint
    REAL(real_8),INTENT(IN)                  :: dfnl_p(il_dfnl_packed(1),nstate), &
                                                fnl_p(il_fnl_packed(1),nstate), &
                                                fnlgam_p(il_fnl_packed(1),nstate)

    INTEGER                                  :: i, ia, is, isa0, k, ia_sum, iac, &
                                                start_ia, end_ia, &
                                                offset_fnl, offset_dfnl, offset_base, &
                                                methread, ispin, start_isa, &
                                                isub, ierr
    INTEGER(int_8)                           :: il_fiont(4),il_ft(3)
    INTEGER, ALLOCATABLE                     :: na_grp(:,:,:), na(:,:)

#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: fiont(:,:,:,:),ft(:,:,:)
#else
    REAL(real_8), ALLOCATABLE                :: fiont(:,:,:,:),ft(:,:,:)
#endif
    REAL(real_8)                             :: weight,fac
    INTEGER, PARAMETER                       :: nteams=1024
    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlfl'
    
    CALL tiset(procedureN,isub)

    IF (cntl%tfdist) CALL stopgm(procedureN,'TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (imagp.EQ.2)CALL stopgm(procedureN,'K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)

    IF(pslo_com%tivan) THEN
       ALLOCATE(na_grp(2,ions1%nsp,0:parai%cp_nogrp-1), na(2,ions1%nsp), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'allocation problem',&
            __LINE__,__FILE__)

       ! split atoms between cp groups
       CALL cp_grp_split_atoms(na_grp)
       na(:,:)=na_grp(:,:,parai%cp_inter_me)
       IF(parai%cp_nogrp.GT.1.AND.parai%cp_inter_me.GT.0)THEN
          !$omp parallel do private(is,ia)
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                fion(1:3,ia,is)=0._real_8
             END DO
          END DO
       END IF

       il_fiont(1)=3
       il_fiont(2)=maxsys%nax
       il_fiont(3)=ions1%nsp+10 !padding
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       il_fiont(4)=1!nteams
#else
       il_fiont(4)=parai%ncpus
#endif

       il_ft(1)=maxsys%nax
       il_ft(2)=3
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       il_ft(3)=nteams
#else
       il_ft(3)=parai%ncpus
#endif
       
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_fiont,fiont,procedureN//'_fiont',ierr)
#else
       ALLOCATE(fiont(il_fiont(1),il_fiont(2),il_fiont(3),il_fiont(4)), stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fiont',& 
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_ft,ft,procedureN//'_ft',ierr)
#else
       ALLOCATE(fiont(il_ft(1),il_ft(2),il_ft(3)), stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ft',& 
            __LINE__,__FILE__)

       methread=1
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute parallel do collapse(3)
       DO is=1,ions1%nsp
          DO ia=1,maxsys%nax
             DO k=1,3
                fiont(k,ia,is,1)=0._real_8
             END DO
          END DO
       END DO
       !$omp target teams num_teams(nteams) thread_limit(32) &
#else
       !$omp parallel &
#endif
       !$omp& private (methread,is,ia,k,isa0,offset_base,ia_sum,&
       !$omp& i,offset_fnl,offset_dfnl,weight,start_isa,ispin) !reduction(+:fiont)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$ methread=omp_get_team_num()+1
#else
       !$ methread=omp_get_thread_num()+1
       fiont(:,:,1:ions1%nsp,methread)=0.0_real_8
#endif

#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp distribute
#else
       !$omp do
#endif
       do i=1,nstate
          weight=crge%f(i,1)
          ispin=1
          IF (cntl%tlsd.AND.i.GT.spin_mod%nsup) ispin=2

          isa0=0
          offset_base=0

          DO is=1,ions1%nsp
             ia_sum=na(2,is)-na(1,is)+1
             start_isa=isa0+na(1,is)-1
             IF( pslo_com%tvan(is))THEN
                IF (ia_sum.GT.0) THEN            

                   offset_fnl=offset_base+1
                   offset_dfnl=offset_base*3+1
                   !DIR$ forceinline
                   CALL add_force(ft(:,:,methread),dfnl_p(offset_dfnl,i),fnl_p(offset_fnl,i),&
                        fnlgam_p(offset_fnl,i),qq(:,:,is),deeq(:,:,:,ispin),dvan(:,:,is), &
                        ia_sum,weight,nlps_com%ngh(is),start_isa,maxsys%nax)
                   DO ia=1,ia_sum
                      DO k=1,3
#if defined(_HAS_OMP_TARGET_OFFLOAD)
                         !$omp atomic
                         fiont(k,ia,is,1)=fiont(k,ia,is,1)+ft(ia,k,methread)
#else
                         fiont(k,ia,is,methread)=fiont(k,ia,is,methread)+ft(ia,k,methread)
#endif
                      END DO
                   END DO
                END IF
             END IF
             offset_base=offset_base+nlps_com%ngh(is)*ia_sum
             isa0=isa0+ions0%na(is)
          ENDDO
       END do
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp end target  teams
       !$omp target update from(fiont)
       DO methread=1,1
#else
       !$omp end parallel
       DO methread=1,parai%ncpus
#endif
          DO is=1,ions1%nsp
             DO ia=na(1,is),na(2,is)
                iac=ia-na(1,is)+1
                DO k=1,3
                   fion(k,ia,is)=fion(k,ia,is)-fiont(k,iac,is,methread)
                END DO
             END DO
          END DO
       END DO

#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_ft,ft,procedureN//'_ft',ierr)
#else
       DEALLOCATE(fiont, stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate ft',& 
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_fiont,fiont,procedureN//'_fiont',ierr)
#else
       DEALLOCATE(fiont, stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fiont',& 
            __LINE__,__FILE__)
       IF (parai%cp_nogrp.GT.1 ) THEN
          CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%cp_inter_grp)
       END IF

       DEALLOCATE(na_grp, na, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'deallocation problem',&
            __LINE__,__FILE__)

    END IF
    IF(cntl%tverbosefor)THEN
       CALL print_debug_ions('DEBUG FORCES '//procedureN, fion)
    END IF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlfl
  ! ==================================================================
  !DIR$ ATTRIBUTES FORCEINLINE::add_force
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  SUBROUTINE add_force(fion,dfnl,fnl,fnlgam,qq_,deeq_,dvan_,ia_sum,weight,ngh,start_isa,ld_f)
#else
  PURE SUBROUTINE add_force(fion,dfnl,fnl,fnlgam,qq_,deeq_,dvan_,ia_sum,weight,ngh,start_isa,ld_f)
#endif
    INTEGER,INTENT(IN)                       :: ia_sum, ngh, start_isa, ld_f
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: qq_(:,:), deeq_(:,:,:), dvan_(:,:)
    REAL(real_8),INTENT(OUT)                 :: fion(ld_f,3)
    REAL(real_8),INTENT(IN)                  :: dfnl(ia_sum,3,ngh,*), fnl(ia_sum,ngh,*), &
                                                fnlgam(ia_sum,ngh,*), weight
    REAL(real_8)                             :: t1, fac1, fac2
    INTEGER                                  :: iv, jv, ia, k, isa

    fion=0.0_real_8
    fac2=weight*2.0_real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp parallel private(iv,fac1,t1,k,ia,isa,jv)
#endif
    DO iv=1,ngh
       !special case iv.eq.jv=> qq always .gt. 0 
       fac1=qq_(iv,iv)*2.0_real_8
       t1=dvan_(iv,iv)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp do collapse(2)
#endif
       DO k=1,3
          DO ia=1,ia_sum
             isa=start_isa+ia
             fion(ia,k)=fion(ia,k)+&
                  (fnlgam(ia,iv,1)*fac1+&
                  fac2*fnl(ia,iv,1)*&
                  (deeq_(isa,iv,iv)+t1))*&
                  dfnl(ia,k,iv,1)
          END DO
       END DO
       DO jv=iv+1,ngh
          fac1=qq_(jv,iv)*2.0_real_8
          t1=dvan_(jv,iv)
          IF (ABS(fac1).GT.2.e-5_real_8) THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp do collapse(2)
#endif
             DO k=1,3
                DO ia=1,ia_sum
                   isa=start_isa+ia
                   fion(ia,k)=fion(ia,k)+&
                        (fnlgam(ia,iv,1)*fac1+&
                        fac2*fnl(ia,iv,1)*&
                        (deeq_(isa,jv,iv)+t1))*&
                        dfnl(ia,k,jv,1)+&
                        (fnlgam(ia,jv,1)*fac1+&
                        fac2*fnl(ia,jv,1)*&
                        (deeq_(isa,jv,iv)+t1))*&
                        dfnl(ia,k,iv,1)
                END DO
             END DO
          ELSE
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp do collapse(2)
#endif
             DO k=1,3
                DO ia=1,ia_sum
                   isa=start_isa+ia
                   fion(ia,k)=fion(ia,k)+&
                        (fnl(ia,iv,1)*dfnl(ia,k,jv,1)+&
                        fnl(ia,jv,1)*dfnl(ia,k,iv,1))*&
                        fac2*(deeq_(isa,jv,iv)+t1)
                END DO
             END DO
          END IF
       END DO
    END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end parallel
#endif
  END SUBROUTINE add_force
  ! ==================================================================
END MODULE rnlfl_utils
