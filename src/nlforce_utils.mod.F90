#include "cpmd_global.h"

MODULE nlforce_utils
  USE cppt,                            ONLY: twnl_nghtol, twnl
  USE cp_grp_utils,                    ONLY: cp_grp_redist,&
                                             cp_grp_split_atoms,&
                                             cp_grp_get_sizes
  USE beta_utils,                      ONLY: build_beta
  USE cvan,                            ONLY: deeq,&
                                             dvan,&
                                             qq
  USE distribution_utils,              ONLY: dist_entity
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8,&
                                             int_8,&
                                             int_4
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlps_com,&
                                             wsg
  !$ USE omp_lib,                      ONLY: omp_set_nested,&
  !$                                         omp_get_thread_num,&
  !$                                         omp_set_num_threads
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: eigr,&
                                             fnl,&
                                             il_fnl_packed
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap
  USE tbxc,                            ONLY: toldcode
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif
#ifdef _INTEL_MKL
  use mkl_service
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: nlforce

CONTAINS


  ! ==================================================================
  SUBROUTINE nlforce(c2,f,fnl_p,fnlgam_p,nstate,redist)
    ! == Rewritten: Tobias Kloeffel, FAU Erlangen-Nuernberg May 2019  ==
    ! == Introduce overlapping communication computation algorithm    ==
    ! == enable via USE_OVERLAPPING_COMM_COMP in the &CPMD section    ==
    ! == relies on fnl_p and fnlgam_p calculated by                   ==
    ! == rottr_c0_fnl and/or rnlsm1, rotate_c0_fnl                    ==
    ! == cp_grp distribution of fnl/fnlgam is taken care of           ==
    ! == full performance only with saved arrays or scratch_library   ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate
    COMPLEX(real_8),INTENT(INOUT)            :: c2(ncpw%ngw,nstate)
    REAL(real_8),INTENT(IN)                  :: f(*), &
                                                fnlgam_p(il_fnl_packed(1),nstate), &
                                                fnl_p(il_fnl_packed(1),nstate)
    LOGICAL,INTENT(IN)                       :: redist

    character(*), PARAMETER                  :: proceduren = 'nlforce'
    INTEGER                                  :: i, ispin, offset_fnl, offset_dai, isa0, &
                                                is, ia_sum, start_isa, ia_fnl, isub, ierr, &
                                                ld_grp(0:parai%cp_nogrp-1), grp, nthreads, &
                                                ibeg, ngw_local, methread, nested_threads
    INTEGER, ALLOCATABLE                     :: na_grp(:,:,:), na(:,:), na_fnl(:,:)
    INTEGER(int_8)                           :: il_eiscr(2), il_dai(3)
    REAL(real_8)                             :: ffi

#ifdef _USE_SCRATCHLIBRARY
    COMPLEX(real_8),POINTER __CONTIGUOUS &
                       , ASYNCHRONOUS        :: eiscr(:,:)
    REAL(real_8),POINTER __CONTIGUOUS &
                       , ASYNCHRONOUS        :: dai(:,:,:)
#else
    REAL(real_8), ALLOCATABLE &
                       , ASYNCHRONOUS        :: dai(:,:,:)
    COMPLEX(real_8),ALLOCATABLE &
                       , ASYNCHRONOUS        :: eiscr(:,:)
#endif
    LOGICAL                                  :: mixed_psp

    CALL tiset(procedureN,isub)
    IF (imagp.EQ.2) call stopgm(procedureN,'k-point not implemented',&
         __LINE__,__FILE__)
    IF (cntl%tfdist) call stopgm(procedureN,'fnl dist. not implemented',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! fixme: ak 2008/02/14: this was broken by changes from 2006/05/12.
    ! not sure whether those were the real cause.
    ! undoing those changes does not make it work
    ! again. so we stop here until somebody fixes it
    ! update ak 2008/05/24: it looks as IF it can be worked around using
    ! oldcode. it seems to be faster, too.
!    IF (pslo_com%tivan.and.cntl%tlsd) THEN
!       IF (.NOT.toldcode)&
!            CALL stopgm(procedureN,'vanderbilt with lsd requires USE of'&
!            // ' "oldcode" flag in &dft section.',&
!            __LINE__,__FILE__)
!    ENDIF
    ! ==--------------------------------------------------------------==
    ALLOCATE(na_grp(2,ions1%nsp,0:parai%cp_nogrp-1), na(2,ions1%nsp), na_fnl(2,ions1%nsp),&
         stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'allocation problem)',&
         __LINE__,__FILE__)

    CALL cp_grp_split_atoms(na_grp)
    na_fnl(:,:)=na_grp(:,:,parai%cp_inter_me)
    !if cp_groups are used we exchange the dai arrays between cp_groups
    !this way we can do the exchange phase during the local dgemm calculation
    !if we need the full c2 we apply all dai arrays to the full set of
    !ngws, else we just update the cp_grp local ngws
    IF(redist)THEN
       ngw_local=ncpw%ngw
       ibeg=1
    ELSE
       CALL cp_grp_get_sizes(ngw_l=ngw_local,first_g=ibeg)
    END IF
    IF (pslo_com%tivan) THEN
       ! ==--------------------------------------------------------==
       ! ==  vanderbilt pp                                         ==
       ! ==--------------------------------------------------------==
       ld_grp=0
       DO grp=0,parai%cp_nogrp-1
          DO is=1,ions1%nsp
             IF (pslo_com%tvan(is) ) THEN
                ld_grp(grp)=ld_grp(grp)+(na_grp(2,is,grp)-na_grp(1,is,grp)+1)&
                     *nlps_com%ngh(is)
             ELSE
                !filter out non uspp atoms
                na_grp(1,is,grp)=0
                na_grp(2,is,grp)=-1
             END IF
          END DO
       END DO
       na(:,:)=na_grp(:,:,parai%cp_inter_me)
       il_dai(1)=MAXVAL(ld_grp)
       il_dai(2)=nstate
       il_dai(3)=parai%cp_nogrp
       il_eiscr(1)=ngw_local
       il_eiscr(2)=MAXVAL(ld_grp)
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_dai,dai,procedureN//'_dai',ierr)
#else
       ALLOCATE(dai(il_dai(1),il_dai(2),il_dai(3)), stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate dai',&
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_eiscr,eiscr,procedureN//'_eiscr',ierr)
#else
       ALLOCATE(eiscr(il_eiscr(1),il_eiscr(2)), stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate eiscr',&
            __LINE__,__FILE__)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute map(to:f(1:nstate),cntl%tlsd)&
#else
       !$omp parallel do &
#endif
       !$omp& private (i,ffi,ispin,offset_fnl,offset_dai,isa0,is,ia_fnl,ia_sum,&
       !$omp& start_isa)
       DO i=1,nstate
          !setup spin settings
          ffi=f(i)
          IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
          ispin=1
          IF (cntl%tlsd.AND.i.GT.spin_mod%nsup) ispin=2

          !offset for packed fnl/fnlgam
          offset_fnl=1
          !offset for dai
          isa0=0
          offset_dai=1
          !fill local part of dai
          DO is=1,ions1%nsp
             ia_fnl=na_fnl(2,is)-na_fnl(1,is)+1
             ia_sum=na(2,is)-na(1,is)+1
             start_isa=isa0+na(1,is)-1
             IF(ia_sum.GT.0)THEN
                !DIR$ forceinline
                CALL build_dai(dai(offset_dai:,i,&
                     parai%cp_inter_me+1),&
                     fnl_p(offset_fnl,i),&
                     fnlgam_p(offset_fnl,i),&
                     ia_sum,ffi,start_isa,nlps_com%ngh(is),qq(:,:,is),dvan(:,:,is),&
                     deeq(:,:,:,ispin),maxsys%nhxs,ions1%nat)
                offset_dai=offset_dai+nlps_com%ngh(is)*ia_sum
             END IF
             offset_fnl=offset_fnl+nlps_com%ngh(is)*ia_fnl
             isa0=isa0+ions0%na(is)
          END DO
       END DO
       IF(cntl%overlapp_comm_comp)THEN
          nthreads=MIN(2,parai%ncpus)
          nested_threads=(MAX(parai%ncpus-1,1))
#if !defined(_INTEL_MKL)
          CALL stopgm(procedureN, 'Overlapping communication and computation: Behavior of &
               BLAS routine inside parallel region not checked',&
               __LINE__,__FILE__)
#endif
       ELSE
          nthreads=1
          nested_threads=parai%ncpus
       END IF
       methread=0
       IF(parai%cp_nogrp.EQ.1)THEN
          nthreads=1
          nested_threads=parai%ncpus
       END IF
       methread=0
       !$omp parallel if(nthreads.EQ.2) num_threads(nthreads) &
       !$omp private(methread,grp) proc_bind(close)
       !$ methread = omp_get_thread_num()
       IF(methread.EQ.0.AND.parai%cp_nogrp.GT.1)THEN
          !get data from other cp_grp other threads build local beta and perform dgemms
          CALL my_concat_inplace(dai,INT(il_dai(1),kind=int_4)*nstate,parai%cp_inter_grp)
       END IF
       IF(methread.EQ.1.OR.nthreads.EQ.1)THEN
          !$ methread = omp_get_thread_num()
          !$ IF(methread.EQ.1)THEN
          !$    CALL omp_set_max_active_levels(2)
          !$    CALL omp_set_num_threads(nested_threads)
#if defined(_INTEL_MKL) && !defined(_HAS_OMP_TARGET_OFFLOAD)
          !$    CALL mkl_set_dynamic(0)
          !$    ierr = mkl_set_num_threads_local(nested_threads)
#endif
          !$ END IF
          grp=parai%cp_inter_me
          CALL build_beta(na_grp(:,:,grp),eigr,eiscr,ncpw%ngw,ibeg,&
                  ngw_local,twnl_nghtol=twnl_nghtol(:,:,:,1))
          IF(ld_grp(grp).GT.0)THEN
             CALL cpmd_dgemm('N','N',2*ngw_local,nstate,ld_grp(grp)&
                  ,1._real_8,eiscr(1,1),2*ngw_local&
                  ,dai(1,1,grp+1),INT(il_dai(1),kind=int_4),1.0_real_8,c2(ibeg,1),2*ncpw%ngw)
          END IF

          !$ IF (methread.EQ.1) THEN
          !$    CALL omp_set_max_active_levels(1)
          !$    CALL omp_set_num_threads(parai%ncpus)
#if defined(_INTEL_MKL) && !defined(_HAS_OMP_TARGET_OFFLOAD)
          !$    CALL mkl_set_dynamic(1)
          !$    ierr = mkl_set_num_threads_local(0)
#endif
          !$ END IF
       END IF
       !$omp end parallel
       !now we apply the dai arrays of all other groups
       IF(parai%cp_nogrp.GT.1)THEN
          DO grp=0,parai%cp_nogrp-1
             IF(grp.EQ.parai%cp_inter_me)CYCLE
             CALL build_beta(na_grp(:,:,grp),eigr,eiscr,ncpw%ngw,ibeg,&
                  ngw_local,twnl_nghtol=twnl_nghtol(:,:,:,1))
             IF(ld_grp(grp).GT.0)THEN
                CALL cpmd_dgemm('N','N',2*ngw_local,nstate,ld_grp(grp)&
                     ,1._real_8,eiscr(1,1),2*ngw_local&
                     ,dai(1,1,grp+1),INT(il_dai(1),kind=int_4),1.0_real_8,c2(ibeg,1),2*ncpw%ngw)
             END IF
          END DO
       END IF
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_eiscr,eiscr,procedureN//'_eiscr',ierr)
#else
    DEALLOCATE(eiscr, stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate eiscr',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_dai,dai,procedureN//'_dai',ierr)
#else
    DEALLOCATE(dai, stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate dai',&
         __LINE__,__FILE__)
    DEALLOCATE(na_grp, na, na_fnl, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'deallocation problem)',&
         __LINE__,__FILE__)

    CALL tihalt(procedureN,isub)
    !check for mixed psp and call nlforce old
    mixed_psp=.FALSE.
    DO is=1,ions1%nsp
       IF(pslo_com%tvan(is))CYCLE
       mixed_psp=.TRUE.
    END DO
    IF(mixed_psp) CALL nlforce_old(c2,f,nstate)
    RETURN
  END SUBROUTINE nlforce
  ! ==================================================================
  !DIR$ ATTRIBUTES FORCEINLINE::build_dai
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  SUBROUTINE build_dai(dai,fnl_,fnlgam_,ia_,ffi,start_isa,ngh_is,qq_,dvan_,deeq_,maxnhg,nat)
#else
  PURE SUBROUTINE build_dai(dai,fnl_,fnlgam_,ia_,ffi,start_isa,ngh_is,qq_,dvan_,deeq_,maxnhg,nat)
#endif
    INTEGER,INTENT(IN)                       :: ia_,start_isa,ngh_is,maxnhg,nat
    REAL(real_8),INTENT(IN)                  :: ffi
    REAL(real_8),INTENT(OUT)                 :: dai(ia_,ngh_is,*)
    REAL(real_8),INTENT(IN)                  :: fnl_(ia_,ngh_is,*),&
                                                fnlgam_(ia_,ngh_is,*),&
                                                qq_(maxnhg,*),dvan_(maxnhg,*),deeq_(nat,maxnhg,*)
    REAL(real_8)                             :: t1, fac1
    INTEGER                                  :: iv,jv,ia,isa
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp parallel do private(iv,jv,ia,isa,t1,fac1)
#endif
    DO iv=1,ngh_is
       DO ia=1,ia_
          dai(ia,iv,1)=0.0_real_8
       END DO
       DO jv=1,ngh_is
          fac1=qq_(jv,iv)
          t1=dvan_(jv,iv)
          isa=start_isa
          IF (ABS(fac1).GT.1.e-5_real_8) THEN
             DO ia=1,ia_
                isa=isa+1
                dai(ia,iv,1)=dai(ia,iv,1)-&
                     fnl_(ia,jv,1)*ffi*&
                     (deeq_(isa,jv,iv)+t1)-&
                     fac1*fnlgam_(ia,jv,1)
             END DO
          ELSE
             DO ia=1,ia_
                isa=isa+1
                dai(ia,iv,1)=dai(ia,iv,1)-&
                     fnl_(ia,jv,1)*ffi*&
                     (deeq_(isa,jv,iv)+t1)
             END DO
          END IF
       END DO
    END DO

    RETURN
  END SUBROUTINE build_dai
  ! ==================================================================

  SUBROUTINE nlforce_old(c2,f,nstate)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(INOUT)            :: c2(ncpw%ngw,*)
    REAL(real_8),INTENT(IN)                  :: f(*)
    INTEGER,INTENT(IN)                       :: nstate

    COMPLEX(real_8)                          :: ct, ctm
    INTEGER                                  :: i, ia, ig, is, isa, isa0, &
                                                ispin, isub, iv, j, jv, ki, &
                                                kj, l, l2, li, lj, nmax, nmin, &
                                                ierr, na(2,ions1%nsp), &
                                                na_grp(2,ions1%nsp,0:parai%cp_nogrp-1)
    REAL(real_8)                             :: dd, fac, ffi
    character(*), PARAMETER                  :: proceduren = 'nlforce_old'
    REAL(real_8),ALLOCATABLE                 :: ddia(:,:)
    COMPLEX(real_8),ALLOCATABLE              :: auxc(:,:)

! Variables
! ==--------------------------------------------------------------==
! == Compute the non-local contributions to the electron gradient ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF (imagp.EQ.2) CALL stopgm(procedureN,'K-POINT NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    IF (cntl%tfdist) CALL stopgm(procedureN,'FNL DIST. NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    ALLOCATE(ddia(maxsys%nax,nstate), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ddai',&
         __LINE__,__FILE__)
    ALLOCATE(auxc(ncpw%ngw,nstate), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate auxc',&
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    ! FIXME: AK 2008/02/14: this was broken by changes from 2006/05/12.
    ! not sure whether those were the real cause.
    ! undoing those changes does not make it work
    ! again. so we stop here until somebody fixes it
    ! UPDATE AK 2008/05/24: it looks as if it can be worked around using
    ! OLDCODE. it seems to be faster, too.
    CALL cp_grp_split_atoms(na_grp)
    na(:,:)=na_grp(:,:,parai%cp_inter_me)

!    IF (pslo_com%tivan.AND.cntl%tlsd) THEN
!       IF (.NOT.toldcode)&
!            CALL stopgm(procedureN,'VANDERBILT WITH LSD REQUIRES USE OF'&
!            // ' "OLDCODE" FLAG IN &DFT SECTION.',&
!            __LINE__,__FILE__)
!    ENDIF
    ! ==--------------------------------------------------------------==
    isa0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          ! ==--------------------------------------------------------==
          ! ==  VANDERBILT PP                                         ==
          ! ==--------------------------------------------------------==
          !see new code
       ELSEIF (sgpp1%tsgp(is)) THEN
          ! ==--------------------------------------------------------==
          ! ==  Stefan Goedecker PP                                   ==
          ! ==--------------------------------------------------------==
          DO iv=1,nlps_com%ngh(is)
             CALL zeroing(auxc)!,ngw*nstate)
             CALL zeroing(ddia)!,maxsys%nax*nstate)
             DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                ffi=f(i)
                IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
                DO ia=na(1,is),na(2,is)
                   isa=isa0+ia
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   dd=0.0_real_8
                   DO jv=1,nlps_com%ngh(is)
                      l2=nghtol(jv,is)+1
                      lj=sgpp2%lpval(jv,is)
                      IF (l.EQ.l2.AND.li.EQ.lj) THEN
                         kj=sgpp2%lfval(jv,is)
                         ddia(ia,i)=ddia(ia,i)+&
                              fnl(1,isa,jv,i,1)*sgpp2%hlsg(ki,kj,l,is)
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             CALL mp_sum(ddia,nstate*maxsys%nax,parai%cp_grp)
             IF (ncpw%ngw.GT.0)&
                  CALL cpmd_dgemm("N","N",2*ncpw%ngw,nstate,ions0%na(is),1._real_8,&
                  eigr(1,isa0+1,1),2*ncpw%ngw,ddia(1,1),maxsys%nax,0._real_8,auxc,2*ncpw%ngw)
             ct=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)
             !$omp parallel do private(IG,I,FFI,CTM)
             DO i=1,nstate
                ffi=f(i)
                IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
                ctm=-ffi*ct
                DO ig=1,ncpw%ngw
                   auxc(ig,i)=ctm*auxc(ig,i)
                   c2(ig,i)=c2(ig,i)+twnl(ig,iv,is,1)*auxc(ig,i)
                ENDDO
             ENDDO
          ENDDO
       ELSE
          ! ==--------------------------------------------------------==
          ! ==  BACHELET HAMANN SCHLUTER                              ==
          ! ==--------------------------------------------------------==
          DO iv=1,nlps_com%ngh(is)
             CALL zeroing(ddia)!,maxsys%nax*nstate)
             DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
                DO ia=na(1,is),na(2,is)
                   isa=isa0+ia
                   ddia(ia,i)=fnl(1,isa,iv,i,1)
                END DO
             END DO
             CALL mp_sum(ddia,nstate*maxsys%nax,parai%cp_grp)
             IF (ncpw%ngw.GT.0)&
                  CALL cpmd_dgemm('N','N',2*ncpw%ngw,nstate,ions0%na(is),1._real_8,&
                  eigr(1,isa0+1,1),2*ncpw%ngw,ddia,maxsys%nax,0.0_real_8,auxc,2*ncpw%ngw)
             ct=(0.0_real_8,-1.0_real_8)**nghtol(iv,is)*wsg(is,iv)
             !$omp parallel do private(IG,I,FFI,CTM)
             DO i=1,nstate
                ffi=f(i)
                IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
                ctm=-ffi*ct
                DO ig=1,ncpw%ngw
                   auxc(ig,i)=ctm*auxc(ig,i)
                   c2(ig,i)=c2(ig,i)+twnl(ig,iv,is,1)*auxc(ig,i)
                ENDDO
             ENDDO
          ENDDO
          ! ==--------------------------------------------------------------==
       ENDIF
       isa0=isa0+ions0%na(is)
    ENDDO
    DEALLOCATE(ddia, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate ddai',&
         __LINE__,__FILE__)
    DEALLOCATE(auxc, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate auxc',&
         __LINE__,__FILE__)

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE nlforce_old
  ! ==================================================================
END MODULE nlforce_utils
