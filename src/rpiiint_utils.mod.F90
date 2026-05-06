#include "cpmd_global.h"

MODULE rpiiint_utils
  USE cnst,                            ONLY: dsqrtpi_inv
  USE eam,                             ONLY: tieam
  USE eam_pot_utils,                   ONLY: eam_pot
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_sum
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE pbc_utils,                       ONLY: pbc
  !$ USE omp_lib,                       ONLY: omp_get_thread_num
  USE ragg,                            ONLY: raggio
  USE special_functions,               ONLY: cp_erfc
  USE system,                          ONLY: iatpe_cp,&
                                             parm, maxsys, iatpt, cntl
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: print_debug_ions
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rpiiint

CONTAINS
  ! ==================================================================
  SUBROUTINE rpiiint(esr,tau0,fion,iesr,tfor)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! == THE RESIDUAL PART OF THE ION-ION INTERACTION (ESR)           ==
    ! == DUE TO THE OVERLAP OF THE SMEARED IONIC CHARGE DENSITIES     ==
    ! == (ionic point charges replaced by gaussian charge distrib.)   ==
    ! == CORRESPONDING TO DIFFERENT ATOMS.                            ==
    ! == ESR depends only on TAU0 (RAGGIO and VALENCE CHARGES)        ==
    ! ==--------------------------------------------------------------==
    ! Modified: Tobias Kloeffel, Erlangen
    ! Date May 2019
    ! cp_grp parallelization, openmp, vectorization(??)
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: tau0(:,:,:)
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    REAL(real_8),INTENT(OUT)                 :: esr
    INTEGER,INTENT(IN)                       :: iesr
    LOGICAL,INTENT(IN)                       :: tfor

    REAL(real_8), PARAMETER                  :: argmax = 20._real_8

    INTEGER                                  :: iat, inf, ishft, isub, ix, &
                                                iy, iz, j, k, l, lax, m, ierr, methread, &
                                                ia, is, thread, new_atom_positions
    INTEGER(int_8)                           :: il_ftmp(4), il_rxlm(3), il_ht(2), il_erre2(3)
    INTEGER, SAVE                            :: iflag = 0
    LOGICAL                                  :: tzero
    REAL(real_8) :: addesr, addpre, arg,  esrtzero, rckj, repand, rlm, &
      xlm, ylm, zlm, xlm_, ylm_, zlm_ , zv2, fiont(6), &
        thresh, rckj_inv
    REAL(real_8),ALLOCATABLE, SAVE           :: tau0_save(:,:,:)
    REAL(real_8),SAVE                        :: esr_save
    LOGICAL,SAVE                             :: first
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8),POINTER __CONTIGUOUS        :: ftmp(:,:,:,:),rxlm(:,:,:),erre2(:,:,:),&
                                                ht(:,:)
#else
    REAL(real_8),ALLOCATABLE                 :: ftmp(:,:,:,:),rxlm(:,:,:),erre2(:,:,:),&
                                                ht(:,:)
#endif
    INTEGER                                  :: iesr_arr((iesr*2+1)**3,3), num_ind, &
                                                tot_ind,ind
    CHARACTER(*), PARAMETER                  :: procedureN='rpiiint'

    ! ==--------------------------------------------------------------==
    IF (iflag.EQ.0.AND.paral%parent) THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,T52,I2,A,I2,A,I2,A)')&
            ' EWALD| SUM IN REAL SPACE OVER ',&
            2*iesr+1,'*',2*iesr+1,'*',2*iesr+1,' CELLS'
       iflag=1
    ENDIF
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)

    IF(.NOT.ALLOCATED(tau0_save))THEN
       ALLOCATE(tau0_save(3,maxsys%nax,maxsys%nsx),STAT=ierr)
       tau0_save=-HUGE(0._real_8)
       esr_save=-HUGE(0._real_8)
       first=.TRUE.
    END IF
    IF(first.OR.&
         (.NOT.ANY(tau0_save(1:3,1,1).EQ.tau0(1:3,1,1))))THEN
       first=.FALSE.
       new_atom_positions=1
       !$omp parallel do private(ia,is)
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             tau0_save(1:3,ia,is)=tau0(1:3,ia,is)
          END DO
       END DO
    ELSE
       new_atom_positions=0
       !$omp parallel do private(ia,is) reduction(+:new_atom_positions)
       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             IF(.NOT.ANY(tau0_save(1:3,ia,is).EQ.tau0(1:3,ia,is)))THEN
                new_atom_positions=1
                tau0_save(1:3,ia,is)=tau0(1:3,ia,is)
             END IF
          END DO
       END DO
    END IF

    IF(new_atom_positions.GT.0.OR.TFOR)THEN
       ind=0
       DO IX=-IESR,IESR
          DO IY=-IESR,IESR
             DO IZ=-IESR,IESR
                ind=ind+1
                iesr_arr(ind,1)=ix
                iesr_arr(ind,2)=iy
                iesr_arr(ind,3)=iz
             END DO
          END DO
       END DO
       tot_ind=ind
       !bring ix=iy=iz=0 at the end of the array iesr_arr
       iesr_arr((tot_ind+1)/2,:)=iesr_arr(tot_ind,:)
       iesr_arr(tot_ind,:)=0

       il_ftmp(1)=3
       il_ftmp(2)=maxsys%nax
       il_ftmp(3)=ions1%nsp+10 !padding
       il_ftmp(4)=parai%ncpus
       il_ht(1)=tot_ind
       il_ht(2)=3
       il_rxlm(1)=tot_ind
       il_rxlm(2)=30 !padding
       il_rxlm(3)=parai%ncpus
       il_erre2(1)=tot_ind
       il_erre2(2)=10 !padding
       il_erre2(3)=parai%ncpus

#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_ftmp,ftmp,procedureN//'_ftmp',ierr)
#else
       ALLOCATE(ftmp(il_ftmp(1),il_ftmp(2),il_ftmp(3),il_ftmp(4)), stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ftmp',&
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_ht,ht,procedureN//'_ht',ierr)
#else
       ALLOCATE(ht(il_ht(1),il_ht(2)), stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ht',&
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_rxlm,rxlm,procedureN//'_rxlm',ierr)
#else
       ALLOCATE(rxlm(il_rxlm(1),il_rxlm(2),il_rxlm(3)), stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate rxlm',&
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_erre2,erre2,procedureN//'_erre2',ierr)
#else
       ALLOCATE(erre2(il_erre2(1),il_erre2(2),il_erre2(3)), stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate erre2',&
            __LINE__,__FILE__)
       !$omp simd
       DO ind=1,tot_ind
          ht(ind,1)=iesr_arr(ind,1)*metr_com%ht(1,1)&
               +iesr_arr(ind,2)*metr_com%ht(2,1)+iesr_arr(ind,3)*metr_com%ht(3,1)
          ht(ind,2)=iesr_arr(ind,1)*metr_com%ht(1,2)&
               +iesr_arr(ind,2)*metr_com%ht(2,2)+iesr_arr(ind,3)*metr_com%ht(3,2)
          ht(ind,3)=iesr_arr(ind,1)*metr_com%ht(1,3)&
               +iesr_arr(ind,2)*metr_com%ht(2,3)+iesr_arr(ind,3)*metr_com%ht(3,3)
       END DO

       esr=0._real_8
       methread=1
       !$omp parallel private(thread,is,ia,iat,k,j,zv2,rckj,rckj_inv,thresh,lax,l,&
       !$omp inf,m,methread,xlm,ylm,zlm,xlm_,ylm_,zlm_,tzero,fiont,num_ind,ind,&
       !$omp rlm,arg,esrtzero,addesr,addpre,repand) reduction(+:esr)
       !$ methread=omp_get_thread_num()+1

       IF(parai%cp_nogrp.GT.1.AND.parai%cp_inter_me.GT.0.AND.tfor)THEN
          !$omp do
          DO is=1,ions1%nsp
             DO ia=1,ions0%na(is)
                fion(1:3,ia,is)=0._real_8
             END DO
          END DO
          !$omp end do nowait
       END IF

       DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
             ftmp(1:3,ia,is,methread)=0._real_8
          END DO
       END DO

       iat=0
       DO k=1,ions1%nsp
          DO j=k,ions1%nsp
             zv2=ions0%zv(k)*ions0%zv(j)
             IF (ABS(zv2).LT.1.e-10_real_8) GOTO 2000
             rckj=SQRT(raggio(k)*raggio(k)+raggio(j)*raggio(j))
             rckj_inv=1.0_real_8/rckj
             thresh=(argmax*rckj)*(argmax*rckj)
             lax=ions0%na(k)
             DO l=1,lax
                IF (iatpe_cp(iat+l,parai%cp_inter_me).NE.parai%mepos) GOTO 1000
                inf=1
                IF (k.EQ.j)inf=l
                !$omp do
                DO M=INF,ions0%NA(J)
                   IF(L.EQ.M.AND.K.EQ.J) THEN
                      xlm=0.d0
                      ylm=0.d0
                      zlm=0.d0
                      TZERO=.TRUE.
                      ESRTZERO=0.5D0
                   ELSE
                      TZERO=.FALSE.
                      ESRTZERO=1.D0
                      xlm_=tau0(1,l,k)-tau0(1,m,j)
                      ylm_=tau0(2,l,k)-tau0(2,m,j)
                      zlm_=tau0(3,l,k)-tau0(3,m,j)
                      CALL pbc(xlm_,ylm_,zlm_,xlm,ylm,zlm,1,parm%apbc,parm%ibrav)
                   ENDIF
                   IF(TFOR) THEN
                      fiont(1:6)=0.0d0
                   ENDIF
                   num_ind=tot_ind
                   !skip last element if TZERO)
                   IF(TZERO)num_ind=num_ind-1
                   !$omp simd
                   DO ind=1,num_ind
                      rxlm(ind,1,methread)=xlm+ht(ind,1)
                      rxlm(ind,2,methread)=ylm+ht(ind,2)
                      rxlm(ind,3,methread)=zlm+ht(ind,3)
                      erre2(ind,1,methread)=&
                           rxlm(ind,1,methread)*rxlm(ind,1,methread)+&
                           rxlm(ind,2,methread)*rxlm(ind,2,methread)+&
                           rxlm(ind,3,methread)*rxlm(ind,3,methread)
                   END DO
                   DO ind=1,num_ind
                      IF(erre2(ind,1,methread).LE.thresh) THEN !ADDESR,ADDPRE /= 0
                         RLM=SQRT(ERRE2(ind,1,methread))
                         ARG=RLM*RCKJ_inv
                         ADDESR=ZV2*ERFC(ARG)/RLM
                         ESR=ESR+ADDESR*ESRTZERO
                         IF(TFOR) THEN
                            ADDPRE=(2.D0*ZV2*DSQRTPI_inv)*DEXP(-ARG*ARG)*RCKJ_inv
                            REPAND=ESRTZERO*(ADDESR+ADDPRE)/ERRE2(ind,1,methread)
                            fiont(1)=fiont(1)+REPAND*RXLM(ind,1,methread)
                            fiont(2)=fiont(2)+REPAND*RXLM(ind,2,methread)
                            fiont(3)=fiont(3)+REPAND*RXLM(ind,3,methread)
                            fiont(4)=fiont(4)-REPAND*RXLM(ind,1,methread)
                            fiont(5)=fiont(5)-REPAND*RXLM(ind,2,methread)
                            fiont(6)=fiont(6)-REPAND*RXLM(ind,3,methread)
                         ENDIF
                      ENDIF
                   END DO

                   IF(TFOR) THEN
                      Ftmp(1,L,K,methread) =Ftmp(1,L,K,methread)+fiont(1)
                      Ftmp(2,L,K,methread) =Ftmp(2,L,K,methread)+fiont(2)
                      Ftmp(3,L,K,methread) =Ftmp(3,L,K,methread)+fiont(3)
                      Ftmp(1,M,J,methread) =Ftmp(1,M,J,methread)+fiont(4)
                      Ftmp(2,M,J,methread) =Ftmp(2,M,J,methread)+fiont(5)
                      Ftmp(3,M,J,methread) =Ftmp(3,M,J,methread)+fiont(6)
                   ENDIF
                ENDDO
                !$omp end do nowait
1000            CONTINUE
             ENDDO
2000         CONTINUE
          ENDDO
          IAT=IAT+ions0%NA(K)
       ENDDO
       !$omp barrier
       !$omp do
       DO is=1,ions1%nsp
          DO thread=1,parai%ncpus
             DO ia=1,ions0%na(is)
                fion(1:3,ia,is)=fion(1:3,ia,is)+ftmp(1:3,ia,is,thread)
             END DO
          END DO
       END DO
       !$omp end do nowait
       !$omp end parallel
       ! 
       ! Embedded Atom Model
       ! 
       IF (tieam) THEN
          CALL eam_pot(esr,tau0,iesr,fion,tfor)
       ENDIF
       !
       CALL mp_sum(esr,parai%cp_grp)
       IF (parai%cp_nogrp.GT.1.AND.TFOR) THEN
          CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%cp_inter_grp)
       END IF
       IF (.NOT.paral%parent) esr=0._real_8
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_erre2,erre2,procedureN//'_erre2',ierr)
#else
       DEALLOCATE(erre2, stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate erre2',&
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_rxlm,rxlm,procedureN//'_rxlm',ierr)
#else
       DEALLOCATE(rxlm, stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate rxlm',&
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_ht,ht,procedureN//'_ht',ierr)
#else
       DEALLOCATE(ht, stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate ht',&
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_ftmp,ftmp,procedureN//'_ftmp',ierr)
#else
       DEALLOCATE(ftmp, stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate ftmp',&
            __LINE__,__FILE__)
       esr_save=esr
    ELSE
       esr=esr_save
    END IF

    IF(tfor.AND.cntl%tverbosefor)THEN
       CALL print_debug_ions('DEBUG FORCES '//procedureN, fion)
    END IF
    CALL tihalt(procedureN,isub)

    ! ==================================================================    
    RETURN
  END SUBROUTINE rpiiint
  ! ==================================================================

END MODULE rpiiint_utils
