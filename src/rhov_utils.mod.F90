#include "cpmd_global.h"

MODULE rhov_utils
  USE cppt,                            ONLY: indz,&
                                             inyh,&
                                             nzh
  USE cp_grp_utils,                    ONLY: cp_grp_redist_array,&
                                             cp_grp_get_sizes,&
                                             cp_grp_split_atoms
  USE cvan,                            ONLY: qg,qg_dipole
  USE distribution_utils,              ONLY: dist_entity
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fftmain_utils,                   ONLY: invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE fft_maxfft,                      ONLY: maxfftn
  USE geq0mod,                         ONLY: geq0
  USE gpu
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nlps_com
  !$ USE omp_lib,                       ONLY: omp_get_thread_num,&
  !$                                          omp_set_max_active_levels,&
  !$                                          omp_get_max_active_levels,&
  !$                                          omp_get_dynamic,&
  !$                                          omp_set_dynamic
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE qvan2_utils,                     ONLY: qvan2
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb,eigrb_dipole,&
                                             fnl, &
                                             fnl_packed,&
                                             il_fnl_packed
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             ncpw, &
                                             cnti,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rhov
CONTAINS

  ! ==================================================================
  SUBROUTINE rhov(nstates,rsumv,psi,hfx,dipole)
    ! ==--------------------------------------------------------------==
    ! ==         THIS ROUTINE CALCULATES THE VANDERBILT DENSITY       ==
    ! ==                                                              ==
    ! == N_V(G) = SUM_I,IJ RHO_I,IJ Q_I,JI(G) E^-IG.R_I               ==
    ! == RHO_I,IJ = SUM_N < BETA_I,I | PSI_N >< PSI_N | BETA_I,J >    ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN) __CONTIGUOUS          :: nstates(:,:)
    REAL(real_8),INTENT(OUT), OPTIONAL       :: rsumv
    COMPLEX(real_8),INTENT(INOUT), TARGET,  &
         OPTIONAL                            :: psi(*)
    LOGICAL, INTENT(IN), OPTIONAL            :: hfx, dipole

    LOGICAL                                  :: do_hfx, do_dipole
    INTEGER                                  :: ig, isub, isub1, ig_start, nhg_loc, is, &
                                                ierr, isub2, nstates_local(2,SIZE(nstates,2))
    INTEGER, ALLOCATABLE                     :: na_grp(:,:,:), na(:,:), nst(:,:)
    INTEGER(int_8)                           :: il_deltar(2)
#ifdef _USE_SCRATCHLIBRARY
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: deltar(:,:)
#else
    COMPLEX(real_8), ALLOCATABLE, TARGET     :: deltar(:,:)
#endif
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: deltar_ptr(:,:)
    CHARACTER(*), PARAMETER                  :: procedureN='rhov'
    IF(cntl%fft_tune_batchsize)THEN
       CALL tiset(procedureN//'tuning',isub2)
    ELSE
       CALL tiset(procedureN,isub)
    END IF

    IF (cntl%tfdist) CALL stopgm(procedureN,'TFDIST NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    IF (imagp.EQ.2)CALL stopgm(procedureN,'K-POINT NOT IMPLEMENTED',&
         __LINE__,__FILE__)

    IF(PRESENT(hfx))THEN
       do_hfx=hfx
    ELSE
       do_hfx=.FALSE.
    END IF
    IF(PRESENT(dipole))THEN
       do_dipole=dipole
    ELSE
       do_dipole=.FALSE.
    END IF

    CALL setfftn(0)

    IF(do_hfx.OR.do_dipole)THEN
       IF (.NOT.PRESENT(PSI)) CALL stopgm(procedureN,'PSI needed',&
            __LINE__,__FILE__)

       IF (.NOT.cntl%bigmem) CALL stopgm(procedureN,'USE BIGMEM FOR HFX',&
            __LINE__,__FILE__)

       !In HFX calculations we do not split atoms between CP groups but distribute
       !state pairs between CP groups
       !Also we just need rho only for a state pair, so no sense to distribute the
       !calculation of states in calc_rho
       !HFX might call with more than one state pair
       ALLOCATE(na(2,ions1%nsp),&
            stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'allocation problem)',&
            __LINE__,__FILE__)

       DO is=1, ions1%nsp
          na(1,is)=1
          na(2,is)=ions0%na(is)
       END DO

       ig_start=1
       nhg_loc=ncpw%nhg

       nstates_local=nstates
       IF(do_dipole)THEN
          nhg_loc=6
       END IF
       deltar_ptr(1:nhg_loc,1:SIZE(nstates,2))=>psi(1:nhg_loc*SIZE(nstates,2))
    ELSE
       
       ALLOCATE(na_grp(2,ions1%nsp,0:parai%cp_nogrp-1), na(2,ions1%nsp), nst(2,0:parai%nproc-1),&
            stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'allocation problem)',&
            __LINE__,__FILE__)

       CALL cp_grp_get_sizes(first_nhg=ig_start,nhg_l=nhg_loc)
       CALL dist_entity(nstates(2,1)-nstates(1,1)+1,parai%nproc,nst)
       !shift everything to starting state
       nstates_local(1,1)=nst(1,parai%me)+nstates(1,1)-1
       nstates_local(2,1)=nst(2,parai%me)+nstates(1,1)-1
       CALL cp_grp_split_atoms(na_grp)
       na(:,:)=na_grp(:,:,parai%cp_inter_me)
       il_deltar(1)=ncpw%nhg
       il_deltar(2)=1
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_deltar,deltar,procedureN//'_deltar',ierr)
#else
       ALLOCATE(deltar(il_deltar(1),il_deltar(2)),STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate deltar', &
            __LINE__,__FILE__)

       deltar_ptr=>deltar
    END IF
   ! ==--------------------------------------------------------------==
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    update_first_to_gpu  =.FALSE.
    update_second_to_gpu =.FALSE.
    update_third_to_gpu  =.FALSE.
    update_result_to_host=.FALSE.
    comm_buffers_on_host =.FALSE.
#endif
    IF (cntl%bigmem) THEN
       CALL prep_bigmem_rhov(nstates_local,ig_start,nhg_loc,na,deltar_ptr,fnl_packed,do_hfx,do_dipole)
    ELSE
       CALL prep_smallmem_rhov(nstates_local,ig_start,nhg_loc,na,deltar_ptr,fnl_packed)
    ENDIF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    update_first_to_gpu  =.TRUE.
    update_second_to_gpu =.TRUE.
    update_third_to_gpu  =.TRUE.
    update_result_to_host=.TRUE.
    comm_buffers_on_host =.TRUE.
#endif
    IF(.NOT.do_hfx.AND..NOT.do_dipole)THEN
       IF (parai%cp_nogrp.GT.1) THEN
          CALL TISET(procedureN//'_grpsb',isub1)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target update from(deltar)
#endif         
          CALL cp_grp_redist_array(deltar(:,1),ncpw%nhg)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target update to(deltar)
#endif         
          CALL TIHALT(procedureN//'_grpsb',isub1)
       END IF
       IF(PRESENT(rsumv))THEN
          IF (geq0) THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp target update from(deltar(1,1))
#endif
             rsumv=REAL(deltar(1,1))
          ELSE
             rsumv=0.0_real_8
          ENDIF
       END IF
    ! ==--------------------------------------------------------------==

       IF (PRESENT(psi)) THEN

#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target teams distribute parallel do
          do ig=1,maxfftn
             psi(ig) = cmplx(0.0_real_8,0.0_real_8)
          end do
          !$omp target teams distribute parallel do
#else
          CALL zeroing(psi(:maxfftn))!,maxfft)
          !$omp parallel do private(IG) shared(PSI)
#endif
          DO ig=1,ncpw%nhg
             psi(nzh(ig))=deltar(ig,1)
             psi(indz(ig))=CONJG(deltar(ig,1))
          ENDDO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target update from(psi(:maxfftn))
#endif
       END IF
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_deltar,deltar,procedureN//'_deltar',ierr)
#else
       DEALLOCATE(deltar,STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate deltar', &
            __LINE__,__FILE__)

       IF(PRESENT(psi)) CALL invfftn(psi(:ncpw%nhg), .FALSE.,parai%allgrp)
    END IF
    IF(do_hfx.OR.do_dipole)THEN
       DEALLOCATE(na, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'deallocation problem)',&
            __LINE__,__FILE__)
    ELSE
       DEALLOCATE(na_grp, na, nst, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'deallocation problem)',&
            __LINE__,__FILE__)
    END IF

    ! ==--------------------------------------------------------------==
    IF(cntl%fft_tune_batchsize)THEN
       CALL tihalt(procedureN//'tuning',isub2)
    ELSE
       CALL tihalt(procedureN,isub)
    END IF
    ! ==--------------------------------------------------------------==
!    if(do_hfx) stop
    RETURN
  END SUBROUTINE rhov
  ! ==================================================================
  SUBROUTINE prep_bigmem_rhov(nst,ig_start,nhg,na,deltar,fnl_p,do_hfx,do_dipole)
    !driver routine for bigmem case
    !handles memory allocations and offset counters
    COMPLEX(real_8),INTENT(OUT)              :: deltar(ncpw%nhg,*)
    INTEGER,INTENT(IN) __CONTIGUOUS          :: nst(:,:)
    INTEGER,INTENT(IN)                       :: ig_start, nhg, na(2,ions1%nsp)
    REAL(real_8),INTENT(IN)                  :: fnl_p(il_fnl_packed(1),*)
    LOGICAL, INTENT(IN)                      :: do_hfx, do_dipole

    INTEGER                                  :: is, isa0, nhh0, offset_fnl0, nhh, num_orb, ipair, &
                                                blocksize, blocks, last_block, ig, ierr, num_pairs, nthreads
    INTEGER(int_8)                           :: il_dia(2), il_ctmp(2), il_fnlt(2)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: dia(:,:), fnlt(:,:)
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: ctmp(:,:)
#else
    REAL(real_8), ALLOCATABLE                :: dia(:,:), fnlt(:,:)
    COMPLEX(real_8), ALLOCATABLE             :: ctmp(:,:)
#endif
    CHARACTER(*), PARAMETER                  :: procedureN='prep_bigmem_rhov'


#if defined(_HAS_OMP_TARGET_OFFLOAD)
    nthreads=1
    blocksize=nhg
#else
    nthreads=parai%ncpus
    blocksize=cnti%blocksize_uspp
#endif
    blocksize=(((blocksize * 8 + 511) / 512) * 512 + 64) / 8
    blocksize=blocksize*nthreads
    IF(blocksize.GT.nhg)blocksize=nhg
    blocks=nhg/blocksize
    last_block=MOD(nhg,blocksize)
    
    IF(do_hfx.OR.do_dipole)THEN
       num_orb=1
    ELSE
       num_orb=nst(2,1)-nst(1,1)+1
    END IF
    num_pairs=SIZE(nst,2)
    
    il_dia(1)=0
    DO is=1,ions1%nsp
       IF(pslo_com%tvan(is))THEN
          il_dia(1)=MAX(il_dia(1),ions0%na(is)*(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2)
       END IF
    END DO
    il_dia(2)=num_pairs
    
    il_ctmp(1)=ceiling(real(blocksize,real_8)/real(nthreads,real_8)) &
         *(maxsys%nhxs*(maxsys%nhxs+1))/2 &
         *num_pairs
    il_ctmp(1)=(((il_ctmp(1) * 8 + 511) / 512) * 512 + 64) / 8
    il_ctmp(2)=nthreads
    IF(do_hfx.OR.do_dipole)THEN
       il_fnlt=1
    ELSE
       il_fnlt(1)=num_orb*maxsys%nhxs*maxsys%nax
       il_fnlt(1)=(((il_fnlt(1) * 8 + 511) / 512) * 512 + 64) / 8
       il_fnlt(2)=nthreads
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_dia,dia,procedureN//'_dia',ierr)
#else
    ALLOCATE(dia(il_dia(1),il_dia(2)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate dia', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_fnlt,fnlt,procedureN//'_fnlt',ierr)
#else
    ALLOCATE(fnlt(il_fnlt(1),il_fnlt(2)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate fnlt', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_ctmp,ctmp,procedureN//'_ctmp',ierr)
#else
    ALLOCATE(ctmp(il_ctmp(1),il_ctmp(2)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'cannot allocate ctmp', &
         __LINE__,__FILE__)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute parallel do simd collapse(2) private(ipair,ig)
#endif
    DO ipair=1,num_pairs
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do private(ig)
#endif
       DO ig=ig_start,ig_start+nhg-1
          deltar(ig,ipair)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
       END DO
    END DO
    isa0=1
    nhh0=1
    offset_fnl0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          nhh=(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2
          CALL evaluate_bigmem_rhov(nst,isa0,ions0%na(is),na(2,is)-na(1,is)+1,nlps_com%ngh(is),&
               nhh,offset_fnl0,na(1,is)-1,num_orb,ig_start,blocksize,last_block,blocks,&
               crge%f,dia,deltar,fnlt,fnl_p,ctmp,nhh0,do_hfx,do_dipole,nthreads)
          nhh0=nhh0+nhh
       END IF
       isa0=isa0+ions0%na(is)
       offset_fnl0=offset_fnl0+(na(2,is)-na(1,is)+1)*nlps_com%ngh(is)
    END DO

#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_ctmp,ctmp,procedureN//'_ctmp',ierr)
#else
    DEALLOCATE(ctmp,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate ctmp', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_fnlt,fnlt,procedureN//'_fnlt',ierr)
#else
    DEALLOCATE(fnlt,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate fnlt', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_dia,dia,procedureN//'_dia',ierr)
#else
    DEALLOCATE(dia,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate dia', &
         __LINE__,__FILE__)

  END SUBROUTINE prep_bigmem_rhov
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE evaluate_bigmem_rhov(nst,isa0,ia_sum,ia_fnl,ngh,nhh,offset_fnl0,offset_ylm,num_orb,&
       ig_start,blocksize,last_block,blocks,f,dia,deltar,fnlt,fnl_p,ctmp,nhh0,do_hfx,do_dipole,nthreads)
    INTEGER,INTENT(IN)                       :: isa0, ia_sum, ia_fnl, ngh, nhh, offset_fnl0,&
                                                offset_ylm, num_orb, ig_start, blocksize, &
                                                last_block, blocks, nhh0, nthreads
    INTEGER,INTENT(IN) __CONTIGUOUS          :: nst(:,:)
    REAL(real_8), INTENT(OUT)                :: dia(ia_sum,nhh,*)
    REAL(real_8), INTENT(OUT) __CONTIGUOUS   :: fnlt(:,:)
    REAL(real_8), INTENT(IN)                 :: fnl_p(il_fnl_packed(1),*), f(*)
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: ctmp(:,:)
    COMPLEX(real_8),INTENT(INOUT)            :: deltar(*)
    LOGICAL, INTENT(IN)                      :: do_hfx, do_dipole
    
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: eigrb_ptr(:,:), qg_ptr(:,:)
    INTEGER                                  :: methread, istart, qgstart, iblock, num_pairs, &
                                                loc_block, my_start, my_end, my_size,ierr

    num_pairs=SIZE(nst,2)
    IF(do_dipole)THEN
       eigrb_ptr => eigrb_dipole
       qg_ptr => qg_dipole
    ELSE
       eigrb_ptr => eigrb
       qg_ptr => qg
    END IF
    methread=1

#if !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp parallel private(methread) num_threads(nthreads)
    !$ methread=omp_get_thread_num()+1
#endif
    IF(do_hfx.OR.do_dipole)THEN
       IF(num_pairs.EQ.1)THEN
          CALL calc_rho12(nst(1,1),nst(2,1),ngh,ia_fnl,offset_fnl0,dia,fnl_p)
       ELSE IF(num_pairs.EQ.2)THEN
          CALL calc_rho13(nst(1,1),nst(2,1),nst(2,2),ngh,ia_fnl,offset_fnl0,dia,fnl_p)
       END IF
    ELSE
       CALL calc_rho(nst(1,1),nst(2,1),ngh,ia_fnl,offset_fnl0,num_orb,offset_ylm,ia_sum,dia,&
            fnlt(:,methread),f,fnl_p,nhh)
    END IF
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end parallel
#endif
    IF(.NOT.do_hfx.AND..NOT.do_dipole) THEN
       CALL mp_sum(dia,nhh*ia_sum,parai%cp_grp)
    END IF
    methread=1
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp parallel private(loc_block,istart,qgstart,methread,iblock,my_start,my_end,my_size)
#endif
    loc_block=ceiling(real(blocksize,real_8)/real(nthreads,real_8))
    istart=ig_start
    qgstart=ig_start-1
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$ methread=omp_get_thread_num()+1
#endif
    
    DO iblock=1,blocks
       my_start=istart+loc_block*(methread-1)
       my_end=istart+loc_block*methread-1
       IF(my_end.GE.istart+blocksize-1)my_end=0
       IF(methread.EQ.nthreads)my_end=blocksize+istart-1
       my_size=my_end-my_start+1

       IF(my_size.GT.0)THEN
          CALL rho_evaluate(my_size,my_start,ia_sum,nhh,dia,ctmp(:,methread),&
               deltar,eigrb_ptr(:,isa0:),qg_ptr(:,nhh0:),num_pairs)
       END IF
       qgstart=qgstart+blocksize
       istart=istart+blocksize
    END DO
    IF(last_block.gt.0)THEN
       loc_block=CEILING(REAL(last_block,real_8)/REAL(nthreads,real_8))
       my_start=istart+loc_block*(methread-1)
       my_end=istart+loc_block*methread-1
       IF(my_end.GE.istart+last_block-1)my_end=0
       IF(methread.EQ.nthreads)my_end=last_block+istart-1
       my_size=my_end-my_start+1
       IF(my_size.GT.0)THEN
          CALL rho_evaluate(my_size,my_start,ia_sum,nhh,dia,ctmp(:,methread),&
               deltar,eigrb_ptr(:,isa0:),qg_ptr(:,nhh0:),num_pairs)
       END IF
    END IF
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end parallel
#endif
  END SUBROUTINE evaluate_bigmem_rhov
  ! ==================================================================
  SUBROUTINE prep_smallmem_rhov(nst,ig_start,nhg,na,deltar,fnl_p)
    COMPLEX(real_8),INTENT(OUT)              :: deltar(ncpw%nhg)

    INTEGER,INTENT(IN)                       :: ig_start, nhg, na(2,ions1%nsp)
    INTEGER,INTENT(IN) __CONTIGUOUS          :: nst(:,:)
    REAL(real_8),INTENT(IN)                  :: fnl_p(il_fnl_packed(1),*)

    INTEGER                                  :: isa0, is, iv, jv, ia, isa, i, ig, nhh, ijv, &
                                                methread, offset_fnl0, offset_dia, &
                                                il_ctmp(1), il_dia(1), il_fnlt(3), ierr
    REAL(real_8)                             :: sum
    REAL(real_8), ALLOCATABLE                :: dia(:), fnlt(:,:,:)
    COMPLEX(real_8), ALLOCATABLE             :: ctmp(:)
    CHARACTER(*), PARAMETER                  :: procedureN='prep_smallmem_rhov'


    il_dia(1)=0
    DO is=1,ions1%nsp
       IF(pslo_com%tvan(is))THEN
          il_dia(1)=MAX(il_dia(1),ions0%na(is)*(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2)
       END IF
    END DO
    il_fnlt(1)=nst(2,1)-nst(1,1)+1
    il_fnlt(2)=maxsys%nhxs
    il_fnlt(3)=parai%ncpus
    il_ctmp(1)=ncpw%nhg

    ALLOCATE(dia(il_dia(1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(fnlt(il_fnlt(1),il_fnlt(2),il_fnlt(3)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(ctmp(il_ctmp(1)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    ALLOCATE(qg(ncpw%nhg,1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    !$omp parallel do private(ig)
    DO ig=ig_start,ig_start+nhg-1
       deltar(ig)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
    END DO

    isa0=0
    offset_fnl0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          nhh=(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2
          methread=1
          !$omp parallel private(methread)
          !$ methread=omp_get_thread_num()+1
          CALL calc_rho(nst(1,1),nst(2,1),nlps_com%ngh(is),na(2,is)-na(1,is)+1,offset_fnl0,&
               nst(2,1)-nst(1,1)+1,na(1,is)-1,ions0%na(is),dia,fnlt(:,:,methread),crge%f,fnl_p,nhh)
          !$omp end parallel
          CALL mp_sum(dia,nhh*ions0%na(is),parai%cp_grp)
          ijv=0
          DO iv=1,nlps_com%ngh(is)
             DO jv=iv,nlps_com%ngh(is)
                ijv=ijv+1
                CALL qvan2(iv,jv,is,qg(:,1))
                DO ia=1,ions0%na(is)
                   isa=isa0+ia
                   !$omp parallel do private(IG)
                   DO ig=ig_start,ig_start+nhg-1
                      ctmp(ig)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                           ei3(isa,inyh(3,ig))
                   ENDDO
                   offset_dia=ions0%na(is)*(ijv-1)+ia
                   !$omp parallel do private(IG)
                   DO ig=ig_start,ig_start+nhg-1
                      deltar(ig)=deltar(ig)+qg(ig,1)*dia(offset_dia)*ctmp(ig)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       offset_fnl0=offset_fnl0+(na(2,is)-na(1,is)+1)*nlps_com%ngh(is)
       isa0=isa0+ions0%na(is)
    ENDDO
    DEALLOCATE(qg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(ctmp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(dia,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    DEALLOCATE(fnlt,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)

  END SUBROUTINE prep_smallmem_rhov
  ! ==================================================================
  SUBROUTINE calc_rho(n1,n2,ngh,ia_fnl,offset_fnl0,norb,offset_ylm,ld_ylm,ylm,temp,f,fnl_p,nhh)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                         :: ngh,ia_fnl,offset_fnl0,norb,offset_ylm,&
                                                  ld_ylm,n1,n2,nhh
    REAL(real_8),INTENT(IN)                    :: fnl_p(il_fnl_packed(1),*),f(*)
    REAL(real_8),INTENT(OUT)                   :: ylm(ld_ylm,*),temp(norb,ngh,*)
    INTEGER                                    :: ia,iv,i,ijv,jv,fnl_ia,ii,ijv_t
    REAL(real_8)                               :: ftmp,tmp
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute parallel do private(ia,ijv)
#else
    !$omp do
#endif
    DO ijv=1,nhh
       DO ia=1,ld_ylm
          ylm(ia,ijv)=0.0_real_8
       END DO
    END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute &
    !$omp& private(ia,iv,i,ijv,jv,fnl_ia,ii,ijv_t,ftmp,tmp) map(to:f(n1:n2)) &
    !$omp& firstprivate(offset_fnl0,norb,n1,n2,ngh,offset_ylm,ia_fnl)
#else
    !$omp do
#endif
    DO ia=1,ia_fnl
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do collapse(2) private (iv,i,ii,fnl_ia) 
#endif
       DO iv=1,ngh
          DO i=n1,n2
             fnl_ia=ia+(iv-1)*ia_fnl+offset_fnl0
             ii=i-n1+1
             temp(ii,iv,ia)=fnl_p(fnl_ia,i)
          END DO
       END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do private (ijv,ijv_t,iv,jv,tmp,i,ii,ftmp)
#endif
       do ijv=1,nhh
          ijv_t=0
          iv_loop:do iv=1,ngh
             do jv=iv,ngh
                ijv_t=ijv_t+1
                if(ijv.eq.ijv_t) exit iv_loop
             end do
          end do iv_loop
          tmp=0._real_8
          DO i=1,norb
             ii=i-1+n1
             tmp=tmp+f(ii)*temp(i,iv,ia)*temp(i,jv,ia)
          END DO
          ftmp=1.0_real_8
          IF(iv.NE.jv) ftmp=2._real_8
          YLM(ia+offset_ylm,ijv)=tmp*ftmp
       END DO
    END DO
  END SUBROUTINE calc_rho
  ! ==================================================================
  SUBROUTINE calc_rho12(i,j,ngh,ia_fnl,offset_fnl0,ylm,fnl_p)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                         :: ngh,ia_fnl,offset_fnl0,i,j
    REAL(real_8),INTENT(IN)                    :: fnl_p(il_fnl_packed(1),*)
    REAL(real_8),INTENT(OUT)                   :: ylm(ia_fnl,*)
    INTEGER                                    :: offset_fnl_iv,offset_fnl_jv,ia,iv,ijv,jv,fnl_ia,&
                                                  ijv_t, nijv
    REAL(real_8)                               :: ftmp


    nijv=ngh*(ngh+1)/2
    !$omp do
    do ijv=1,nijv
       ijv_t=0
       iv_loop:do iv=1,ngh
          do jv=iv,ngh
             ijv_t=ijv_t+1
             if(ijv.eq.ijv_t) exit iv_loop
          end do
       end do iv_loop
       offset_fnl_iv=offset_fnl0+(iv-1)*ia_fnl          
       offset_fnl_jv=offset_fnl0+(jv-1)*ia_fnl
       ftmp=1.0_real_8
       if(iv.EQ.jv) ftmp=0.5_real_8
       !$omp simd
       do ia=1,ia_fnl
          ylm(ia,ijv)=ftmp*(fnl_p(offset_fnl_iv+ia,i)*fnl_p(offset_fnl_jv+ia,j)+&
               fnl_p(offset_fnl_jv+ia,i)*fnl_p(offset_fnl_iv+ia,j))
       end do
    end do
  END SUBROUTINE calc_rho12
  ! ==================================================================
  SUBROUTINE calc_rho13(i,j,k,ngh,ia_fnl,offset_fnl0,ylm,fnl_p)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                         :: ngh,ia_fnl,offset_fnl0,i,j,k
    REAL(real_8),INTENT(IN)                    :: fnl_p(il_fnl_packed(1),*)
    REAL(real_8),INTENT(OUT)                   :: ylm(ia_fnl,ngh*(ngh+1)/2,*)
    INTEGER                                    :: offset_fnl_iv,offset_fnl_jv,ia,iv,ijv,jv,fnl_ia,&
                                                  ijv_t, nijv
    REAL(real_8)                               :: ftmp


    nijv=ngh*(ngh+1)/2
    !$omp do
    do ijv=1,nijv
       ijv_t=0
       iv_loop:do iv=1,ngh
          do jv=iv,ngh
             ijv_t=ijv_t+1
             if(ijv.eq.ijv_t) exit iv_loop
          end do
       end do iv_loop
       offset_fnl_iv=offset_fnl0+(iv-1)*ia_fnl
       offset_fnl_jv=offset_fnl0+(jv-1)*ia_fnl
       ftmp=1.0_real_8
       if(iv.EQ.jv) ftmp=0.5_real_8
       !$omp simd
       do ia=1,ia_fnl
          ylm(ia,ijv,1)=ftmp*(fnl_p(offset_fnl_iv+ia,i)*fnl_p(offset_fnl_jv+ia,j)+&
               fnl_p(offset_fnl_jv+ia,i)*fnl_p(offset_fnl_iv+ia,j))
          ylm(ia,ijv,2)=ftmp*(fnl_p(offset_fnl_iv+ia,i)*fnl_p(offset_fnl_jv+ia,k)+&
               fnl_p(offset_fnl_jv+ia,i)*fnl_p(offset_fnl_iv+ia,k))
       end do
    end do
  END SUBROUTINE calc_rho13
  SUBROUTINE rho_evaluate(block_size,block_start,na_is,nhh,dia,ctmp,deltar,eigrb_local,qg_local,num_rho)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                         :: block_size, block_start, na_is, nhh, &
                                                  num_rho
    REAL(real_8),INTENT(IN)                    :: DIA(na_is,*)
    COMPLEX(real_8),INTENT(IN) , CONTIGUOUS    :: eigrb_local(:,:),qg_local(:,:)
    COMPLEX(real_8),INTENT(OUT)                :: ctmp(block_size,nhh,*)
    COMPLEX(real_8),INTENT(INOUT)              :: deltar(size(eigrb_local,1),*)
    INTEGER                                    :: ijv,nhh0,ig,ig2,irho

    IF(block_size.GT.0)THEN
       CALL cpmd_dgemm('n','n',2*block_size,nhh*num_rho,na_is &
            ,1.0_real_8,EIGRB_local(block_start,1),2*size(EIGRB_local,1) &
            ,DIA,na_is,0.0_real_8,CTMP,2*size(ctmp,1))
       DO irho=1,num_rho
          DO ijv=1,nhh
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp target teams distribute parallel do simd private(ig2)
#else
             !$omp simd
#endif
             DO ig=1, block_size
                ig2=ig -1 + block_start
                deltar(ig2,irho)=deltar(ig2,irho)+qg_local(ig2,ijv)*ctmp(ig,ijv,irho)
             ENDDO

          END DO
       END DO
    END IF
  END SUBROUTINE rho_evaluate
  ! ==================================================================



END MODULE rhov_utils
