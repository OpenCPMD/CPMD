#include "cpmd_global.h"

MODULE rnlfor_utils
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms
  USE cvan,                            ONLY: deeq,&
                                             deeq_fnl_hfx, &
                                             dvan
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8, &
                                             int_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlm,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: dfnl,&
                                             fnl,&
                                             fnl2, &
                                             fnla,&
                                             !dfnla,&
                                             !dfnl_packed,&
                                             il_dfnl_packed
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             ipept,&
                                             maxsys,&
                                             parap
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: print_debug_ions
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch !, &
                                             !request_saved_scratch, &
                                             !save_scratch
#endif


!$ USE omp_lib,                         ONLY: omp_get_thread_num

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rnlfor
  !public :: rcasfor
  PUBLIC :: rnlfor_hfx
CONTAINS

  ! ==================================================================
  SUBROUTINE rnlfor(fion,f,wk,nstate,nkpoint)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NON-LOCAL POTENTIAL CONTRIBUTION TO THE FORCE ON THE    ==
    ! ==  IONIC DEGREES OF FREEDOM                                    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    INTEGER                                  :: nstate, nkpoint
    REAL(real_8)                             :: f(nstate,nkpoint), wk(nkpoint)

    INTEGER                                  :: i, ia, ii, ik, is, isa, isa0, &
                                                ispin, isub, iv, jv, ki, kj, &
                                                l, l2, li, lj
    REAL(real_8)                             :: tdbl, temp, tt, weight, &
                                                wk1_1, wk1_2, wk1_3, wk2_1, &
                                                wk2_2, wk2_3
    REAL(real_8),ALLOCATABLE                 :: fnlt(:,:,:),dfnlt(:,:,:,:)
    CHARACTER(*), PARAMETER                  :: procedureN = 'rnlfor'

! Variables
! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    DO ik=1,nkpoint
       isa0=0
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is)) THEN
             !TK USPP case moved to rnlfl to access dfnl only once
             !This also allows us to call rnlsm2 just in time
          ELSEIF (sgpp1%tsgp(is)) THEN
             ! Stefan Goedecker pp
             DO iv=1,nlps_com%ngh(is)
                l=nghtol(iv,is)+1
                ki=sgpp2%lfval(iv,is)
                li=sgpp2%lpval(iv,is)
                DO jv=1,nlps_com%ngh(is)
                   l2=nghtol(jv,is)+1
                   lj=sgpp2%lpval(jv,is)
                   IF (l2.EQ.l.AND.li.EQ.lj) THEN
                      kj=sgpp2%lfval(jv,is)
                      DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                         weight=wk(ik)*f(i,ik)
                         IF (ABS(weight).GT.1.e-12_real_8) THEN
                            tt=2.0_real_8*weight*sgpp2%hlsg(ki,kj,l,is)
                            ii=i-parap%nst12(parai%mepos,1)+1
                            IF (imagp.EQ.2) THEN
                               IF (cntl%tfdist) THEN
                                  !$omp parallel do private(IA,ISA)
                                  DO ia=1,ions0%na(is)
                                     isa=isa0+ia
                                     fion(1,ia,is)=fion(1,ia,is)-tt*dfnl(1,isa,&
                                          iv,1,ii,ik)*fnl2(1,isa,jv,  ii,ik)-tt*&
                                          dfnl(2,isa,iv,1,ii,ik)* fnl2(2,isa,jv,&
                                          ii,ik)
                                     fion(2,ia,is)=fion(2,ia,is)-tt*dfnl(1,isa,&
                                          iv,2,ii,ik)*fnl2(1,isa,jv,  ii,ik)-tt*&
                                          dfnl(2,isa,iv,2,ii,ik)* fnl2(2,isa,jv,&
                                          ii,ik)
                                     fion(3,ia,is)=fion(3,ia,is)-tt*dfnl(1,isa,&
                                          iv,3,ii,ik)*fnl2(1,isa,jv,  ii,ik)-tt*&
                                          dfnl(2,isa,iv,3,ii,ik)* fnl2(2,isa,jv,&
                                          ii,ik)
                                  ENDDO
                               ELSE
                                  !$omp parallel do private(IA,ISA)
                                  DO ia=1,ions0%na(is)
                                     isa=isa0+ia
                                     fion(1,ia,is)=fion(1,ia,is)-tt*dfnl(1,isa,&
                                          iv,1,ii,ik)*fnl2(1,isa,jv,   i,ik)-tt*&
                                          dfnl(2,isa,iv,1,ii,ik)* fnl2(2,isa,jv,&
                                          i,ik)
                                     fion(2,ia,is)=fion(2,ia,is)-tt*dfnl(1,isa,&
                                          iv,2,ii,ik)*fnl2(1,isa,jv,   i,ik)-tt*&
                                          dfnl(2,isa,iv,2,ii,ik)* fnl2(2,isa,jv,&
                                          i,ik)
                                     fion(3,ia,is)=fion(3,ia,is)-tt*dfnl(1,isa,&
                                          iv,3,ii,ik)*fnl2(1,isa,jv,   i,ik)-tt*&
                                          dfnl(2,isa,iv,3,ii,ik)* fnl2(2,isa,jv,&
                                          i,ik)
                                  ENDDO
                               ENDIF
                            ELSE
                               IF (cntl%tfdist) THEN
                                  !$omp parallel do private(IA,ISA)
                                  DO ia=1,ions0%na(is)
                                     isa=isa0+ia
                                     fion(1,ia,is)=fion(1,ia,is)-tt*dfnl(1,isa,&
                                          iv,1,ii,ik)*fnl2(1,isa,jv,  ii,ik)
                                     fion(2,ia,is)=fion(2,ia,is)-tt*dfnl(1,isa,&
                                          iv,2,ii,ik)*fnl2(1,isa,jv,  ii,ik)
                                     fion(3,ia,is)=fion(3,ia,is)-tt*dfnl(1,isa,&
                                          iv,3,ii,ik)*fnl2(1,isa,jv,  ii,ik)
                                  ENDDO
                               ELSE
                                  !$omp parallel do private(IA,ISA)
                                  DO ia=1,ions0%na(is)
                                     isa=isa0+ia
                                     fion(1,ia,is)=fion(1,ia,is)-tt*dfnl(1,isa,&
                                          iv,1,ii,ik)*fnl2(1,isa,jv,   i,ik)
                                     fion(2,ia,is)=fion(2,ia,is)-tt*dfnl(1,isa,&
                                          iv,2,ii,ik)*fnl2(1,isa,jv,   i,ik)
                                     fion(3,ia,is)=fion(3,ia,is)-tt*dfnl(1,isa,&
                                          iv,3,ii,ik)*fnl2(1,isa,jv,   i,ik)
                                  ENDDO
                               ENDIF
                            ENDIF
                         ENDIF
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ELSE
             ! Other pp (numeric)
             DO iv=1,nlps_com%ngh(is)
                temp=2._real_8*wsg(is,iv)
                DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                   weight=wk(ik)*f(i,ik)
                   IF (ABS(weight).GT.1.e-12_real_8) THEN
                      ii=i-parap%nst12(parai%mepos,1)+1
                      IF (cntl%tfdist) THEN
                         IF (imagp.EQ.2) THEN
                            !$omp parallel do private(IA,ISA,wk1_1,wk1_2,wk1_3)
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               wk1_1=dfnl(1,isa,iv,1,ii,ik)*fnl2(1,isa,iv, ii,&
                                    ik)+dfnl(2,isa,iv,1,ii,ik)*fnl2(2,isa,iv,&
                                    ii,ik)
                               wk1_2=dfnl(1,isa,iv,2,ii,ik)*fnl2(1,isa,iv, ii,&
                                    ik)+dfnl(2,isa,iv,2,ii,ik)*fnl2(2,isa,iv,&
                                    ii,ik)
                               wk1_3=dfnl(1,isa,iv,3,ii,ik)*fnl2(1,isa,iv, ii,&
                                    ik)+dfnl(2,isa,iv,3,ii,ik)*fnl2(2,isa,iv,&
                                    ii,ik)
                               fion(1,ia,is)=fion(1,ia,is)-temp*weight*wk1_1
                               fion(2,ia,is)=fion(2,ia,is)-temp*weight*wk1_2
                               fion(3,ia,is)=fion(3,ia,is)-temp*weight*wk1_3
                            ENDDO
                         ELSE
                            !$omp parallel do private(IA,ISA,wk1_1,wk1_2,wk1_3)
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               wk1_1=dfnl(1,isa,iv,1,ii,ik)*fnl2(1,isa,iv,&
                                    ii,ik)
                               wk1_2=dfnl(1,isa,iv,2,ii,ik)*fnl2(1,isa,iv,&
                                    ii,ik)
                               wk1_3=dfnl(1,isa,iv,3,ii,ik)*fnl2(1,isa,iv,&
                                    ii,ik)
                               fion(1,ia,is)=fion(1,ia,is)-temp*weight*wk1_1
                               fion(2,ia,is)=fion(2,ia,is)-temp*weight*wk1_2
                               fion(3,ia,is)=fion(3,ia,is)-temp*weight*wk1_3
                            ENDDO
                         ENDIF
                      ELSE
                         IF (imagp.EQ.2) THEN
                            !$omp parallel do private(IA,ISA,wk1_1,wk1_2,wk1_3)
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               wk1_1=dfnl(1,isa,iv,1,ii,ik)*fnl2(1,isa,iv, i,&
                                    ik)+dfnl(2,isa,iv,1,ii,ik)*fnl2(2,isa,iv,&
                                    i,ik)
                               wk1_2=dfnl(1,isa,iv,2,ii,ik)*fnl2(1,isa,iv, i,&
                                    ik)+dfnl(2,isa,iv,2,ii,ik)*fnl2(2,isa,iv,&
                                    i,ik)
                               wk1_3=dfnl(1,isa,iv,3,ii,ik)*fnl2(1,isa,iv, i,&
                                    ik)+dfnl(2,isa,iv,3,ii,ik)*fnl2(2,isa,iv,&
                                    i,ik)
                               fion(1,ia,is)=fion(1,ia,is)-temp*weight*wk1_1
                               fion(2,ia,is)=fion(2,ia,is)-temp*weight*wk1_2
                               fion(3,ia,is)=fion(3,ia,is)-temp*weight*wk1_3
                            ENDDO
                         ELSE
                            !$omp parallel do private(IA,ISA,wk1_1,wk1_2,wk1_3)
                            DO ia=1,ions0%na(is)
                               isa=isa0+ia
                               wk1_1=dfnl(1,isa,iv,1,ii,ik)*fnl2(1,isa,iv,&
                                    i,ik)
                               wk1_2=dfnl(1,isa,iv,2,ii,ik)*fnl2(1,isa,iv,&
                                    i,ik)
                               wk1_3=dfnl(1,isa,iv,3,ii,ik)*fnl2(1,isa,iv,&
                                    i,ik)
                               fion(1,ia,is)=fion(1,ia,is)-temp*weight*wk1_1
                               fion(2,ia,is)=fion(2,ia,is)-temp*weight*wk1_2
                               fion(3,ia,is)=fion(3,ia,is)-temp*weight*wk1_3
                            ENDDO
                         ENDIF
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
          isa0 = isa0 + ions0%na(is)
       ENDDO
    ENDDO
    IF(cntl%tverbosefor)THEN
       CALL print_debug_ions('DEBUG FORCES '//procedureN, fion)
    END IF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlfor
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE rnlfor_hfx(fion,f,wk,nstate,nkpoint,dfnl_packed)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NON-LOCAL POTENTIAL CONTRIBUTION TO THE FORCE ON THE    ==
    ! ==  IONIC DEGREES OF FREEDOM                                    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)
    INTEGER                                  :: nstate, nkpoint
    REAL(real_8)                             :: f(nstate,nkpoint), wk(nkpoint)
    REAL(real_8),INTENT(IN)                  :: dfnl_packed(il_dfnl_packed(1),nstate)
    INTEGER                                  :: i, ia, ii, ik, is, isa, isa0, &
                                                ispin, isub, iv, jv, ki, kj, &
                                                l, l2, li, lj, k, ierr, &
                                                methread, j,&
                                                na_grp(2,ions1%nsp,parai%cp_nogrp),iac,&
                                                ia_sum, start_ia,end_ia,offset1,offset2,&
                                                offset3,offset4,offset_base, tot_work, &
                                                offset, start_isa,ikind
   INTEGER(INT_8)                           :: il_fiont(4)
   REAL(real_8)                             :: tdbl, tt, weight, fac, &
                                                wk1_1, wk1_2, wk1_3, wk2_1, &
                                                wk2_2, wk2_3, wk_t(3), ft(maxsys%nax,3)
   LOGICAL                                  :: need_old
   real(real_8), allocatable :: deeq21(:,:,:,:)
#ifdef _USE_SCRATCHLIBRARY   
   REAL(real_8), POINTER __CONTIGUOUS       :: fiont(:,:,:,:),temp(:,:)
#else
   REAL(real_8),ALLOCATABLE                 :: fiont(:,:,:,:),temp(:,:)
#endif
    CHARACTER(*), PARAMETER                 :: procedureN = 'rnlfor_hfx' 
    ! split atoms between cp groups

    call cp_grp_split_atoms(na_grp)

! Variables
! ==--------------------------------------------------------------==
! If no non-local components -> return.
    IF (nlm.EQ.0) RETURN

    CALL tiset(procedureN,isub)

    need_old=.FALSE.
    IF(imagp.EQ.2)need_old=.TRUE.
    IF(il_dfnl_packed(1).EQ.0)need_old=.TRUE.
    DO is=1,ions1%nsp
       IF(.NOT.pslo_com%tvan(is))need_old=.TRUE.
    END DO
    IF (il_dfnl_packed(1).GT.0.AND.pslo_com%tivan.AND.imagp.EQ.1.AND..NOT.cntl%tfdist) THEN
       !Vanderbilt optimized: rnlsm2 glosum is skipped, instead rnlfor runs over all states
       !and only forces need to be summed up. Glosum of dfnl is always more expensive 
       CALL cp_grp_split_atoms(na_grp)
   
       if (parai%cp_nogrp.gt.1 .and. parai%cp_inter_me .gt. 0) then
          !$omp parallel do private(is,ia)
          do is=1,ions1%nsp
             do ia=1,ions0%na(is)
                fion(1:3,ia,is)=0._real_8
             end do
          end do
       end if
       
       il_fiont(1)=3
       il_fiont(2)=maxsys%nax
       il_fiont(3)=ions1%nsp+10
       il_fiont(4)=parai%ncpus
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_fiont,fiont,procedureN//'_fiont',ierr)
!       CALL request_saved_scratch(il_dfnl_packed,dfnl_packed,'DFNL_packed',ierr)
#else
       ALLOCATE(fiont(il_fiont(1),il_fiont(2),il_fiont(3),il_fiont(4)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate dfnlt',& 
            __LINE__,__FILE__)
#endif
       isa0=0
       ALLOCATE(deeq21(ions1%nat,maxsys%nhxs,1,nstate)&
            ,stat=ierr)
       deeq21=0.0d0
       methread=1
       !$omp parallel private(methread,is,ia,k,ik,isa0,offset_base,start_ia,end_ia,ia_sum,&
       !$omp ft,i,weight,ispin,iv,jv,tdbl,isa,offset1,offset2,iac,start_isa,j,fac)

       !$ methread=omp_get_thread_num()+1
       fiont(:,:,:,methread)=0.0_real_8
       DO ik=1,nkpoint
          isa0=0
          offset_base=0
          DO is=1,ions1%nsp
             IF (pslo_com%tvan(is)) THEN
                !Vanderbilt optimized
                start_ia=na_grp(1,is,parai%cp_inter_me +1)
                end_ia=na_grp(2,is,parai%cp_inter_me +1)
                ia_sum=end_ia-start_ia+1
                start_isa=isa0+na_grp(1,is,parai%cp_inter_me +1)-1
!                !$omp do
!                DO i=1,nstate                   
!                   ispin=1
!                   IF (cntl%tlsd.AND.i.GT.spin_mod%nsup) ispin=2
!                   DO iv=1,nlps_com%ngh(is)
!                      DO jv=1,nlps_com%ngh(is)
!                         isa=start_isa
!                         do j=1,nstate
!                            fac=weight
!                            !$omp simd
!                            DO ia=1,ia_sum
!                               isa=start_isa+ia
!                               deeq21(isa,iv,ispin,i)=deeq21(isa,iv,ispin,i)+&
!                                    fnla(isa,jv,j)*&
!                                    deeq2(isa,jv,iv,ispin,j,i)
!                            END DO
!                         end do
!                      END DO
!                   END DO
!                END DO
!
                ft=0.0_real_8
                !$omp do
                DO i=1,nstate                   
                   weight=wk(ik)*f(i,ik)
                   IF (ABS(weight).GT.1.e-12_real_8) THEN
                      ispin=1
                      IF (cntl%tlsd.AND.i.GT.spin_mod%nsup) ispin=2
                      DO iv=1,nlps_com%ngh(is)
                         DO k=1,3
                            offset1=offset_base+(iv-1)*ia_sum*3+(k-1)*ia_sum
                            isa=start_isa
                            fac=weight
                            !$omp simd
                            DO ia=1,ia_sum
                               isa=start_isa+ia
                               ft(ia,k)=ft(ia,k)+&
                                    dfnl_packed(offset1+ia,i)*&
                                    deeq_fnl_hfx(isa,iv,i)*fac
                            END DO
                         END DO
                      END DO
                   END IF
                END DO
                !$omp end do nowait

                DO ia=1,ia_sum
                   DO k=1,3
                      fiont(k,ia,is,methread)=ft(ia,k)*2.0d0
                   END DO
                END DO
             END IF
             isa0 = isa0 + ions0%na(is)
             offset_base=offset_base+nlps_com%ngh(is)*ia_sum*3
          ENDDO
       ENDDO
       !$omp end parallel
       DO methread=1,parai%ncpus
          DO is=1,ions1%nsp
             DO ia=na_grp(1,is,parai%cp_inter_me+1),na_grp(2,is,parai%cp_inter_me+1)
                iac=ia-na_grp(1,is,parai%cp_inter_me+1)+1
                DO k=1,3
                   fion(k,ia,is)=fion(k,ia,is)-fiont(k,iac,is,methread)
                END DO
             END DO
          END DO
       END DO
       IF (parai%cp_nogrp.GT.1) THEN
          CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%cp_inter_grp)
       END IF
       IF(cntl%tverbosefor)THEN
          CALL print_debug_ions('DEBUG FORCES '//procedureN, fion)
       END IF
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_fiont,fiont,procedureN//'_fiont',ierr)
!       CALL save_scratch(il_dfnl_packed,dfnl_packed,'DFNL_packed',ierr)
#else
       DEALLOCATE(fiont, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fiont',& 
            __LINE__,__FILE__)
#endif
       deallocate(deeq21)
    END IF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rnlfor_hfx
  ! ==================================================================


  SUBROUTINE rcasfor(fion)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE NON-LOCAL POTENTIAL CONTRIBUTION TO THE FORCE ON THE    ==
    ! ==  IONIC DEGREES OF FREEDOM                                    ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'rcasfor'

    INTEGER                                  :: ia, iatl, ierr, ii, is, isa, &
                                                isa0, isub, iv, jv, k, ki, &
                                                kj, l, l2, li, lj
    REAL(real_8)                             :: tt
    REAL(real_8), ALLOCATABLE                :: dfab(:,:,:)

! Variables
! ==--------------------------------------------------------------==
! If no non-local components -> return.

    IF (nlm.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset('   RCASFOR',isub)
    ALLOCATE(dfab(ions1%nat,maxsys%nhxs,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO k=1,3
       CALL zeroing(dfab)!,2*ions1%nat)
       IF (clsd%ialpha.GE.parap%nst12(parai%mepos,1).AND.clsd%ialpha.LE.parap%nst12(parai%mepos,2) ) THEN
          ii=clsd%ialpha-parap%nst12(parai%mepos,1)+1
          CALL dcopy(maxsys%nhxs*ions1%nat,dfnl(1,1,1,k,ii,1),1,dfab(1,1,1),1)
       ENDIF
       IF (clsd%ibeta.GE.parap%nst12(parai%mepos,1).AND.clsd%ibeta.LE.parap%nst12(parai%mepos,2) ) THEN
          ii=clsd%ibeta-parap%nst12(parai%mepos,1)+1
          CALL dcopy(maxsys%nhxs*ions1%nat,dfnl(1,1,1,k,ii,1),1,dfab(1,1,2),1)
       ENDIF
       CALL mp_sum(dfab,2*ions1%nat*maxsys%nhxs,parai%allgrp)
       isa0=0
       DO is=1,ions1%nsp
          IF (pslo_com%tvan(is)) THEN
             ! Vanderbild pp
             CALL stopgm("RCASFOR","VDB NOT IMPLEMENTED",& 
                  __LINE__,__FILE__)
          ELSEIF (sgpp1%tsgp(is)) THEN
             ! Stefan Goedecker pp
             DO iv=1,nlps_com%ngh(is)
                l=nghtol(iv,is)+1
                ki=sgpp2%lfval(iv,is)
                li=sgpp2%lpval(iv,is)
                DO jv=1,nlps_com%ngh(is)
                   l2=nghtol(jv,is)+1
                   lj=sgpp2%lpval(jv,is)
                   IF (l2.EQ.l.AND.li.EQ.lj) THEN
                      kj=sgpp2%lfval(jv,is)
                      tt=sgpp2%hlsg(ki,kj,l,is)
                      DO ia=1,ions0%na(is)
                         isa=isa0+ia
                         IF (cntl%tfdist) THEN
                            iatl=isa-ipept(1,parai%mepos)+1
                            fion(k,ia,is)=fion(k,ia,is)-0.5_real_8*tt*(dfab(isa,iv,&
                                 1)*fnl(1,iatl,jv,clsd%ibeta,1)+dfab(isa,iv,2)*&
                                 fnl(1,iatl,jv,clsd%ialpha,1)+ dfab(isa,jv,1)*&
                                 fnl(1,iatl,iv,clsd%ibeta,1)+dfab(isa,jv,2)*&
                                 fnl(1,iatl, iv, clsd%ialpha,1))
                         ELSE
                            fion(k,ia,is)=fion(k,ia,is)-tt*(dfab(isa,iv,1)*&
                                 fnl(1,isa,jv,clsd%ibeta,1)+dfab(isa,iv,2)*fnl(1,&
                                 isa, jv, clsd%ialpha, 1)+dfab(isa, jv,1)*fnl(1,&
                                 isa,iv, clsd%ibeta,1)+ dfab(isa,jv,2)*fnl(1,isa,&
                                 iv,clsd%ialpha,1))
                         ENDIF
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ELSE
             ! Other pp (numeric)
             DO iv=1,nlps_com%ngh(is)
                DO ia=1,ions0%na(is)
                   isa=isa0+ia
                   IF (isa.GE.ipept(1,parai%mepos).AND.isa.LE.ipept(2,parai%mepos))&
                        THEN
                      IF (cntl%tfdist) THEN
                         iatl=isa-ipept(1,parai%mepos)+1
                         fion(k,ia,is)=fion(k,ia,is)-wsg(is,iv)*(dfab(isa,iv,&
                              1)*fnl(1,iatl,iv,clsd%ibeta,1)+fnl(1,iatl,iv,clsd%ialpha,&
                              1)* dfab(isa,iv,2))
                      ELSE
                         fion(k,ia,is)=fion(k,ia,is)-wsg(is,iv)*(dfab(isa,iv,&
                              1)*fnl(1,isa,iv,clsd%ibeta,1)+fnl(1,isa,iv,clsd%ialpha,&
                              1)* dfab(isa,iv,2))
                      ENDIF
                   ENDIF
                ENDDO
             ENDDO
          ENDIF
          isa0 = isa0 + ions0%na(is)
       ENDDO
    ENDDO
    DEALLOCATE(dfab,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('   RCASFOR',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rcasfor
  ! ==================================================================

END MODULE rnlfor_utils
