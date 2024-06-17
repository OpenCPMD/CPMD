#include "cpmd_global.h"

MODULE newd_utils
  USE cppt,                            ONLY: gk,&
                                             gk_trans,&
                                             inyh
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes,&
                                             cp_grp_split_atoms
  USE cvan,                            ONLY: qg
  USE distribution_utils,              ONLY: dist_entity
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  USE gpu
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: int_8,&
                                             real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: nlps_com
  !$ USE omp_lib,                       ONLY: omp_get_thread_num
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE qvan2_utils,                     ONLY: qvan2
  USE sfac,                            ONLY: ei1,&
                                             ei2,&
                                             ei3,&
                                             eigrb,&
                                             fnl_packed,&
                                             il_fnl_packed
  USE system,                          ONLY: cntl,&
                                             iatpe,&
                                             ipept,&
                                             maxsys,&
                                             ncpw,&
                                             parm,&
                                             cnti
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: print_debug_ions
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: newd

CONTAINS

    ! ==--------------------------------------------------------------==
    ! ==         THIS ROUTINE CALCULATES THE ARRAY DEEQ               ==
    ! ==                                                              ==
    ! ==  DEEQ_I,IJ = OMEGA( V(G=0) Q_I,IJ(G=0) +                     ==
    ! ==       2 SUM_G> RE[ V*(G) Q_I,IJ(G) E^-IG.R_I ] )             ==
    ! ==                                                              ==
    ! ==--------------------------------------------------------------==
  SUBROUTINE newd(deeq,f,vpot,fion,tfor,nstates,hfx)
    REAL(REAL_8),INTENT(INOUT) __CONTIGUOUS  :: deeq(:,:,:)
    COMPLEX(real_8),INTENT(IN)               :: vpot(ncpw%nhg,*)
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    REAL(real_8),INTENT(IN)                  :: f(*)
    LOGICAL,INTENT(IN)                       :: tfor
    INTEGER,INTENT(IN) __CONTIGUOUS          :: nstates(:,:)
    LOGICAL,INTENT(IN), OPTIONAL             :: hfx
    LOGICAL                                  :: do_hfx
    LOGICAL, ALLOCATABLE                     :: sym(:)
    INTEGER                                  :: is, ia, isub, ig_start, nhg_loc, ierr, i, &
                                                nstates_local(SIZE(nstates,1),SIZE(nstates,2))
    INTEGER, ALLOCATABLE                     :: na_grp(:,:,:), na(:,:), nst(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN='newd'

    CALL tiset(procedureN,isub)
    IF (cntl%tfdist) CALL stopgm(procedureN,'TFDIST NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    
    IF(PRESENT(hfx))THEN
       do_hfx=hfx
    ELSE
       do_hfx=.FALSE.
    END IF

    IF(do_hfx)THEN
       IF (.NOT.cntl%bigmem) CALL stopgm(procedureN,'USE BIGMEM FOR HFX',&
            __LINE__,__FILE__)

       !In HFX calculations we do not split atoms between CP groups but distribute
       !state pairs between CP groups
       !Also we just need rho only for a state pair, so no sense to distribute the
       !calculation of states in calc_rho
       !HFX might call with more than one potential and more than one state pair
       ALLOCATE(na(2,ions1%nsp), sym(SIZE(nstates,2)),&
            stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'allocation problem)',&
            __LINE__,__FILE__)

       DO is=1, ions1%nsp
          na(1,is)=1
          na(2,is)=ions0%na(is)
       END DO

       ig_start=1
       nhg_loc=ncpw%nhg

       DO i=1, SIZE(nstates,2)
          IF(nstates(1,i).EQ.nstates(2,i)) THEN
             sym(i)=.TRUE.
          ELSE
             sym(i)=.FALSE.
          END IF
       END DO
       nstates_local=nstates
    ELSE

       ALLOCATE(na_grp(2,ions1%nsp,0:parai%cp_nogrp-1), na(2,ions1%nsp), nst(2,0:parai%nproc-1),&
            sym(1), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'allocation problem)',&
            __LINE__,__FILE__)

       CALL cp_grp_split_atoms(na_grp)
       na(:,:)=na_grp(:,:,parai%cp_inter_me)
       CALL cp_grp_get_sizes(first_nhg=ig_start,nhg_l=nhg_loc)
       !generate 'new' nst12 mapping -> nwa12
       CALL dist_entity(nstates(2,1)-nstates(1,1)+1,parai%nproc,nst)
       !shift to starting state
       nstates_local(1,1)=nst(1,parai%me)+nstates(1,1)-1
       nstates_local(2,1)=nst(2,parai%me)+nstates(1,1)-1

       sym=.FALSE.

       
       IF(tfor)THEN
          IF (parai%cp_nogrp.GT.1 .AND. parai%cp_inter_me .GT. 0) THEN
             !$omp parallel do private(is,ia)
             DO is=1,ions1%nsp
                DO ia=1,ions0%na(is)
                   fion(1:3,ia,is)=0._real_8
                END DO
             END DO
          END IF
       END IF

    END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    update_first_to_gpu  =.FALSE.
    update_second_to_gpu =.FALSE.
    update_third_to_gpu  =.FALSE.
    update_result_to_host=.FALSE.
    comm_buffers_on_host =.FALSE.
#endif
    IF (cntl%bigmem) THEN
       CALL prep_bigmem_newd(nstates_local,SIZE(nstates,2),ig_start,nhg_loc,na,deeq,f,fion,vpot,&
            tfor,do_hfx,sym)
    ELSE
       CALL prep_smallmem_newd(nstates_local,SIZE(nstates,2),ig_start,nhg_loc,na,deeq,f,fnl_packed,&
            fion,vpot,tfor,do_hfx,sym)
    ENDIF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    update_first_to_gpu  =.TRUE.
    update_second_to_gpu =.TRUE.
    update_third_to_gpu  =.TRUE.
    update_result_to_host=.TRUE.
    comm_buffers_on_host =.TRUE.
#endif

    IF(tfor.AND..NOT.do_hfx)THEN
       IF (parai%cp_nogrp.gt.1 ) then
          CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%cp_inter_grp)
       END IF
    END IF
    IF(do_hfx)THEN
       DEALLOCATE(na, sym, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'deallocation problem)',&
            __LINE__,__FILE__)
    ELSE
       DEALLOCATE(na_grp, na, nst, sym, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'deallocation problem)',&
            __LINE__,__FILE__)
    END IF
    IF(tfor.AND.cntl%tverbosefor)THEN
       CALL print_debug_ions('DEBUG FORCES '//procedureN, fion)
    END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target update from(deeq)
#endif
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE newd
  ! ==================================================================
  SUBROUTINE prep_bigmem_newd(nst,num_pot,ig_start,nhg,na,deeq,f,fion,vpot,tfor,do_hfx,sym)
    !driver routine for bigmem case
    !handles memory allocations and offset counters
    INTEGER,INTENT(IN)                       :: num_pot, nst(2,num_pot), ig_start, nhg, na(2,ions1%nsp)
    REAL(REAL_8),INTENT(INOUT) __CONTIGUOUS    :: deeq(:,:,:)
    COMPLEX(real_8),INTENT(IN)               :: vpot(ncpw%nhg,num_pot)
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    REAL(real_8),INTENT(IN)                  :: f(*)
    LOGICAL,INTENT(IN)                       :: tfor,do_hfx
    LOGICAL,INTENT(IN) __CONTIGUOUS          :: sym(:)

#ifdef _USE_SCRATCHLIBRARY
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: qg1(:,:), vtmp(:)
    REAL(real_8), POINTER __CONTIGUOUS       :: fnlt(:,:),ylm(:,:)
#else
    REAL(real_8), ALLOCATABLE                :: fnlt(:,:),ylm(:,:)
    COMPLEX(real_8), ALLOCATABLE             :: qg1(:,:), vtmp(:)
#endif

    INTEGER                                  :: i, is, ia, blocksize, blocks, last_block, &
                                                isa0, nhh0, nhh, num_orb, offset_fnl0, ierr,&
                                                ld_fnlt, nthreads
    INTEGER(int_8)                           :: il_qg1(2), il_ylm(2),&
                                                il_fnlt(2)
    CHARACTER(*), PARAMETER                  :: procedureN='prep_bigmem'

    !blocking parameters
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
    !
    IF(DO_HFX)THEN
       num_orb=1
    ELSE
       num_orb=nst(2,1)-nst(1,1)+1
    END IF

    !Memory allocation
    il_ylm(1)=0
    DO is=1,ions1%nsp
       IF(pslo_com%tvan(is))THEN
          il_ylm(1)=MAX(il_ylm(1),ions0%na(is)*(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2)
       END IF
    END DO
    il_qg1(1)=CEILING(REAL(blocksize,real_8)/REAL(nthreads,real_8)) &
         *((maxsys%nhxs*(maxsys%nhxs+1))/2)*num_pot
    if(tfor)il_qg1(1)=il_qg1(1)*4
    il_qg1(1)=(((il_qg1(1) * 8 + 511) / 512) * 512 + 64) / 8
    il_qg1(2)=nthreads

    il_ylm(2)=num_pot
    !if ionic forces are needed we combine both glosums:
    !ylm(:,:,1) contains rho_nm
    !ylm(:,:,2) contains local part of deeq
    !ylm(:,:,3-5) contains derivative of qg
    !ylm(:,:,2-5) enter dgemm
    il_fnlt(1)=il_ylm(1)*num_pot
    IF (tfor) THEN
       il_ylm(2)=5*num_pot
       il_fnlt(1)=il_fnlt(1)*4
       IF(.NOT.do_hfx)THEN
          il_fnlt(1)=max(il_fnlt(1),maxsys%nax*num_orb*maxsys%nhxs)
       END IF
    END IF
    il_fnlt(1)=(((il_fnlt(1) * 8 + 511) / 512) * 512 + 64) / 8
    il_fnlt(2)=nthreads
    IF(DO_HFX)THEN
       num_orb=1
    ELSE
       num_orb=nst(2,1)-nst(1,1)+1
    END IF

#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_fnlt,fnlt,procedureN//'_fnlt',ierr)
#else
    ALLOCATE(fnlt(il_fnlt(1),il_fnlt(2)), stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlt',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_qg1,qg1,procedureN//'_qg1',ierr)
#else
    ALLOCATE(qg1(il_qg1(1),il_qg1(2)), stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate qg1',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_ylm,ylm,procedureN//'_ylm',ierr)
#else
    ALLOCATE(ylm(il_ylm(1),il_ylm(2)), stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ylm',&
         __LINE__,__FILE__)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target update to(vpot(:,:num_pot),fnl_packed)
#endif
    isa0=0
    nhh0=1
    offset_fnl0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is)) THEN
          IF(do_hfx)THEN
             ld_fnlt=na(2,is)-na(1,is)+1
          ELSE
             ld_fnlt=num_orb
          END IF
          
          nhh=(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2
          CALL evaluate_bigmem_newd(isa0,ions0%na(is),na(2,is)-na(1,is)+1,nlps_com%ngh(is),&
               nhh,offset_fnl0,na(1,is)-1,num_orb,nst,ig_start,blocksize,last_block,blocks,tfor, &
               do_hfx,sym,f,vpot,deeq,qg(1,nhh0),ld_fnlt,fnlt,qg1,ylm,&
               fnl_packed,fion(:,:,is),nthreads)
          nhh0=nhh0+nhh
       END IF
       isa0=isa0+ions0%na(is)
       offset_fnl0=offset_fnl0+(na(2,is)-na(1,is)+1)*nlps_com%ngh(is)
    END DO

#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_ylm,ylm,procedureN//'_ylm',ierr)
#else
    DEALLOCATE(ylm, stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate ylm)',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_qg1,qg1,procedureN//'_qg1',ierr)
#else
    DEALLOCATE(qg1, stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate qg1)',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_fnlt,fnlt,procedureN//'_fnlt',ierr)
#else
    DEALLOCATE(fnlt, stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnlt',&
         __LINE__,__FILE__)


  END SUBROUTINE prep_bigmem_newd
  ! ==================================================================
  SUBROUTINE evaluate_bigmem_newd(isa0,ia_sum,ia_fnl,ngh,nhh,offset_fnl0,offset_ylm,num_orb,nst, &
       ig_start,blocksize,last_block,blocks,tfor,do_hfx,sym,f,vpot,deeq,qg_,ld_fnlt,fnlt,qg1,ylm, &
       fnl_p,fion,nthreads)
    INTEGER,INTENT(IN)                       :: isa0, ia_sum, ia_fnl, ngh, nhh, offset_fnl0,&
                                                 offset_ylm, ig_start, blocksize, &
                                                 last_block, blocks, num_orb, ld_fnlt, nthreads
    INTEGER,INTENT(IN) __CONTIGUOUS          :: nst(:,:)
    LOGICAL,INTENT(IN)                       :: tfor, do_hfx
    LOGICAL,INTENT(IN) __CONTIGUOUS          :: sym(:)
    REAL(real_8),INTENT(IN)                  :: f(*), fnl_p(il_fnl_packed(1),*)
    REAL(real_8),INTENT(OUT)                 :: ylm(ia_sum,nhh,*)
    REAL(real_8), INTENT(INOUT) __CONTIGUOUS :: deeq(:,:,:), fnlt(:,:), fion(:,:)
    COMPLEX(real_8), INTENT(OUT) __CONTIGUOUS:: qg1(:,:)
    COMPLEX(real_8), INTENT(IN)              :: qg_(ncpw%nhg,*), vpot(*)


    INTEGER                                  :: istart, iblock, qgstart, ijv, methread, ia, &
                                                isa, iv, jv, len_qg1, start_ylm, num_pot, &
                                                ipot, my_size, my_start, my_end, loc_block,ierr

    REAL(real_8)                             :: fac, omtpiba, ft1, ft2, ft3, otr

    num_pot=SIZE(nst,2)

    len_qg1=num_pot
    start_ylm=1
    IF(tfor)THEN
       len_qg1=len_qg1*4
       start_ylm=start_ylm+num_pot
    END IF

    omtpiba=parm%omega*parm%tpiba
    methread=1
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp parallel private(loc_block,istart,qgstart,methread,iblock,my_start,my_end,my_size,fac) &
    !$omp& num_threads(nthreads)
#endif
    loc_block=CEILING(REAL(blocksize,real_8)/REAL(nthreads,real_8))
    istart=ig_start
    qgstart=ig_start-1
    methread=1
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$ methread=omp_get_thread_num()+1
#endif
    fac=0.0_real_8

    DO iblock=1,blocks
       my_start=istart+loc_block*(methread-1)
       my_end=istart+loc_block*methread-1
       IF(my_end.GE.istart+blocksize-1)my_end=0
       IF(methread.EQ.nthreads)my_end=blocksize+istart-1
       my_size=my_end-my_start+1

       CALL evaluate_block_newd(my_size,my_start,isa0,ia_sum, &
            nhh,fac,tfor,qg1(:,methread),fnlt(:,methread),vpot,qg_,num_pot,gk_trans)

       IF(iblock.EQ.1)THEN
          fac=1.0_real_8
          !remove double counting of geq0
          IF(geq0.AND.istart.EQ.1.AND.my_start.EQ.istart) THEN
             CALL cpmd_dger(ia_sum,len_qg1*nhh,-1.0_real_8,&
                  eigrb(1,isa0+1),2*ncpw%nhg, &
                  qg1(:,methread),2*my_size,&
                  fnlt(:,methread),ia_sum)
          END IF
       END IF
       istart=istart+blocksize
       qgstart=qgstart+blocksize
    END DO
    IF(last_block.GT.0)THEN
       loc_block=CEILING(REAL(last_block,real_8)/REAL(nthreads,real_8))
       my_start=istart+loc_block*(methread-1)
       my_end=istart+loc_block*methread-1
       IF(my_end.GE.istart+last_block-1)my_end=0
       IF(methread.EQ.nthreads)my_end=last_block+istart-1
       my_size=my_end-my_start+1

       CALL evaluate_block_newd(my_size,my_start,isa0,ia_sum, &
            nhh,fac,tfor,qg1(:,methread),fnlt(:,methread),vpot,qg_,num_pot,gk_trans)
    END IF
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end parallel
#endif
    
    CALL collect_ylm(ia_sum*nhh*len_qg1,ylm(1,1,start_ylm),fnlt,nthreads)
    IF(tfor)THEN
       ijv=0
       methread=1
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel private(methread) num_threads(nthreads)
       !$ methread=omp_get_thread_num()+1
#endif
       IF(do_hfx)THEN
          IF(num_pot.EQ.1)THEN
             CALL calc_rho12(nst(1,1),nst(2,1),ngh,ia_sum,offset_fnl0,&
                  ylm,fnl_p)
          ELSE IF(num_pot.EQ.2)THEN
             CALL calc_rho13(nst(1,1),nst(2,1),nst(2,2),ngh,ia_sum,offset_fnl0,&
                  ylm,fnl_p)
          END IF
       ELSE
          CALL calc_rho(nst(1,1),nst(2,1),ngh,ia_fnl,offset_fnl0,num_orb,offset_ylm,ia_sum,&
               ylm,fnlt(:,methread),f,fnl_p,nhh)
       END IF
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp end parallel
#endif
    END IF

    IF(do_hfx)then
       CALL mp_sum(ylm(:,:,start_ylm),nhh*ia_sum*num_pot,parai%cp_grp)
    ELSE
       CALL mp_sum(ylm,nhh*ia_sum*start_ylm,parai%cp_grp)
    END IF
    IF(do_hfx)THEN
       CALL build_deeq_hfx(ngh, ia_sum, isa0, parm%omega, fnl_p, nst, deeq, &
            ylm(1,1,start_ylm), offset_fnl0,fnlt,nthreads)
    ELSE
       CALL build_deeq_gga(ngh, ia_sum, isa0, parm%omega, deeq, ylm(1,1,start_ylm),nthreads)
    END IF

    IF (tfor) THEN
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel private(ipot,fac,ia,ft1,ft2,ft3,ijv,otr) num_threads(nthreads)
#endif
       DO ipot=1,num_pot
          IF(do_hfx)THEN
             fac=2.0_real_8
             IF(sym(ipot)) fac=1.0_real_8
          ELSE
             fac=1.0_real_8
          END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target teams distribute private(ia,ft1,ft2,ft3,ijv,otr)
#else
          !$omp do schedule(static)
#endif
          DO ia=1,ia_sum
             ft1=0._real_8
             ft2=0._real_8
             ft3=0._real_8
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp parallel do private(ijv,otr) reduction(+:ft1,ft2,ft3)
#endif
             DO ijv=1,nhh
                OTR=YLM(ia,ijv,ipot)*omtpiba
                ft1=ft1+OTR*YLM(ia,ijv,num_pot*2+(ipot-1)*3+1)
                ft2=ft2+OTR*YLM(ia,ijv,num_pot*2+(ipot-1)*3+2)
                ft3=ft3+OTR*YLM(ia,ijv,num_pot*2+(ipot-1)*3+3)
             END DO
             fion(1,ia)=fion(1,ia)+ft1*fac
             fion(2,ia)=fion(2,ia)+ft2*fac
             fion(3,ia)=fion(3,ia)+ft3*fac
          END DO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp end do nowait
#endif
       END DO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp end parallel
#endif

    ENDIF
  END SUBROUTINE evaluate_bigmem_newd

  SUBROUTINE collect_ylm(ld_ylm,ylm,ylm_temp,nthreads)
    INTEGER, INTENT(IN)                  :: ld_ylm,nthreads
    REAL(real_8),INTENT(OUT)             :: ylm(*)
    REAL(real_8),INTENT(IN) __CONTIGUOUS :: ylm_temp(:,:)
    INTEGER                              :: ithread,i
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute parallel do simd
#else
    !$omp parallel private(i,ithread) num_threads(nthreads)
    !$omp do simd schedule(static)
#endif
    DO i=1,ld_ylm
       ylm(i)=0.0_real_8
    END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute private(ithread,i)
#else
    !$omp end do simd nowait
#endif
    DO ithread=1,nthreads
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do private(i)
#else
       !$omp do simd schedule(static)
#endif
       DO i=1,ld_ylm
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp atomic
#endif
          ylm(i)=ylm(i)+ylm_temp(i,ithread)
       END DO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp end do simd nowait
#endif
    END DO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end parallel
#endif

  END SUBROUTINE collect_ylm
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
       ftmp=2.0_real_8
       if(iv.EQ.jv) ftmp=1.0_real_8
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
       ftmp=2.0_real_8
       if(iv.EQ.jv) ftmp=1.0_real_8
       !$omp simd
       do ia=1,ia_fnl
          ylm(ia,ijv,1)=ftmp*(fnl_p(offset_fnl_iv+ia,i)*fnl_p(offset_fnl_jv+ia,j)+&
               fnl_p(offset_fnl_jv+ia,i)*fnl_p(offset_fnl_iv+ia,j))
          ylm(ia,ijv,2)=ftmp*(fnl_p(offset_fnl_iv+ia,i)*fnl_p(offset_fnl_jv+ia,k)+&
               fnl_p(offset_fnl_jv+ia,i)*fnl_p(offset_fnl_iv+ia,k))
       end do
    end do
  END SUBROUTINE calc_rho13

  SUBROUTINE evaluate_block_newd(blocksize,block_start,isa0,ia_sum, &
       nhh,fac,tfor,qg1,ylm,vpot,qg_local,num_pot,gk_trans_)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                         :: blocksize, block_start, isa0, &
                                                  ia_sum, nhh, num_pot
    LOGICAL,INTENT(IN)                         :: tfor
    REAL(real_8),INTENT(IN)                    :: fac
    REAL(real_8),INTENT(IN) __CONTIGUOUS       :: gk_trans_(:,:)
    COMPLEX(real_8),INTENT(IN)                 :: vpot(ncpw%nhg,*), qg_local(ncpw%nhg,*)
    COMPLEX(real_8),INTENT(OUT)                :: qg1(blocksize,nhh,*)
    REAL(real_8)                               :: ylm(ia_sum,nhh,*)
    INTEGER                                    :: ijv,ig,k,ig2,len_qg1,ipot,ind
    IF(num_pot.NE.2)THEN
       IF(tfor)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target teams distribute parallel do simd collapse(3) &
          !$omp& private(ijv,ipot,ig,ind,ig2,k)
#endif
          DO ijv=1,nhh
             DO ipot=1, num_pot
                DO ig=1,blocksize
                   ind=(ipot-1)*3+num_pot
                   ig2=block_start-1+ig
                   qg1(ig,ijv,ipot)=CONJG(qg_local(ig2,ijv))*vpot(ig2,ipot)
                   DO k=1,3
                      qg1(ig,ijv,ind+k)=&
                           CMPLX(gk_trans_(ig2,k)*AIMAG(qg1(ig,ijv,ipot)),&
                           gk_trans_(ig2,k)*(-REAL(qg1(ig,ijv,ipot),kind=real_8)))
                   END DO
                END DO
             END DO
          END DO
       ELSE
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target teams distribute parallel do simd collapse(3) &
          !$omp& private(ijv,ipot,ig,ig2)
#endif
          DO ijv=1,nhh
             DO ipot=1, num_pot
                DO ig=1,blocksize
                   ig2=block_start-1+ig
                   qg1(ig,ijv,ipot)=CONJG(qg_local(ig2,ijv))*vpot(ig2,ipot)
                END DO
             END DO
          END DO
       END IF
    ELSE
       IF(tfor)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target teams distribute parallel do simd collapse(2) &
          !$omp& private(ijv,ipot,ig,ind,ig2,k)
#endif
          DO ijv=1,nhh
             DO ig=1,blocksize
                ig2=block_start-1+ig
                qg1(ig,ijv,1)=CONJG(qg_local(ig2,ijv))*vpot(ig2,1)
                qg1(ig,ijv,2)=CONJG(qg_local(ig2,ijv))*vpot(ig2,2)
                DO k=1,3
                   qg1(ig,ijv,2+k)=&
                        CMPLX(gk_trans_(ig2,k)*AIMAG(qg1(ig,ijv,1)),&
                        gk_trans_(ig2,k)*(-REAL(qg1(ig,ijv,1),kind=real_8)))
                   qg1(ig,ijv,5+k)=&
                        CMPLX(gk_trans_(ig2,k)*AIMAG(qg1(ig,ijv,2)),&
                        gk_trans_(ig2,k)*(-REAL(qg1(ig,ijv,2),kind=real_8)))
                END DO
             END DO
          END DO
       ELSE
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target teams distribute parallel do simd collapse(2) &
          !$omp& private(ijv,ipot,ig,ind,ig2)
#endif
          DO ijv=1,nhh
             DO ig=1,blocksize
                ig2=block_start-1+ig
                qg1(ig,ijv,1)=CONJG(qg_local(ig2,ijv))*vpot(ig2,1)
                qg1(ig,ijv,2)=CONJG(qg_local(ig2,ijv))*vpot(ig2,2)
             END DO
          END DO
       END IF
    END IF
    len_qg1=num_pot
    IF(tfor)len_qg1=4*num_pot
    CALL cpmd_dgemm('T','N',ia_sum,nhh*len_qg1,2*blocksize,2.0_real_8,&
         eigrb(block_start,isa0+1),2*ncpw%nhg,&
         qg1,2*blocksize,&
         fac,ylm,ia_sum)
    RETURN
  END SUBROUTINE evaluate_block_newd
  ! ==================================================================
  SUBROUTINE build_deeq_gga(ngh, ia_sum, isa0, fac,deeq, ylm,nthreads)
    INTEGER, INTENT(IN) :: ngh, ia_sum, isa0, nthreads
    REAL(real_8), INTENT(OUT) __CONTIGUOUS :: deeq(:,:,:)
    REAL(real_8), INTENT(IN)  :: ylm(ia_sum,ngh*(ngh+1)/2)
    REAL(real_8), INTENT(IN)  :: fac
    INTEGER                   :: ijv_c, iv, jv, ia, isa, ijv
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute &
#else
    !$omp parallel do num_threads(nthreads)&
#endif
    !$omp& private(ijv_c, iv, jv, ia, isa, ijv)
    do ijv=1, ngh*(ngh+1)/2
       ijv_c=0
       count: DO iv=1,ngh
          DO jv=iv,ngh
             ijv_c=ijv_c+1
             IF(ijv_c.EQ.ijv) EXIT count
          END DO
       END DO count
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do private(ia,isa)
#endif
       DO ia=1,ia_sum
          isa=isa0+ia
          deeq(isa,iv,jv)=ylm(ia,ijv)*fac !parm%omega
          deeq(isa,jv,iv)=ylm(ia,ijv)*fac !parm%omega
       ENDDO
    ENDDO
  END SUBROUTINE build_deeq_gga
  
  SUBROUTINE build_deeq_hfx(ngh, ia_sum, isa0, fac, fnl_p, nst,deeq,ylm,offset_fnl0,temp,nthreads)
    INTEGER, INTENT(IN) :: ngh, ia_sum, isa0, offset_fnl0, nthreads
    INTEGER, INTENT(IN) __CONTIGUOUS :: nst(:,:)
    REAL(real_8), INTENT(OUT) :: deeq(ions1%nat, maxsys%nhxs,*),temp(ia_sum,ngh,nthreads)
    REAL(real_8), INTENT(IN)  :: ylm(ia_sum,ngh*(ngh+1)/2,*), fnl_p(il_fnl_packed(1), &
                                 il_fnl_packed(2)), fac
    INTEGER :: i,j, ijv_c, iv, jv, ia, isa, ipot, ijv, offset_fnl_iv, offset_fnl_jv, methread, ithread
    real(real_8) :: fac_l

    methread=1
    !$omp parallel private(ipot,i,j,ijv,ijv_c,iv,jv,fac_l,offset_fnl_iv,offset_fnl_jv,ia,isa,methread),&
    !$omp& num_threads(nthreads)
    !$ methread=omp_get_thread_num()+1
    do ipot= 1, SIZE(nst,2)
       i=nst(1,ipot)
       j=nst(2,ipot)
       temp(:,:,methread)=0.0_real_8
       !$omp do 
       do ijv=1, ngh*(ngh+1)/2
          ijv_c=0
          count: DO iv=1,ngh
             DO jv=iv,ngh
                ijv_c=ijv_c+1
                IF(ijv_c.EQ.ijv) EXIT count
             END DO
          END DO count
          fac_l=fac
          IF(iv.EQ.jv)fac_l=fac*0.5_real_8
          offset_fnl_iv=offset_fnl0+(iv-1)*ia_sum
          offset_fnl_jv=offset_fnl0+(jv-1)*ia_sum
          DO ia=1,ia_sum
             temp(ia,jv,methread)=temp(ia,jv,methread)+ylm(ia,ijv,ipot)*fac_l*&
                  fnl_p(offset_fnl_iv+ia,i)
             temp(ia,iv,methread)=temp(ia,iv,methread)+ylm(ia,ijv,ipot)*fac_l*&
                  fnl_p(offset_fnl_jv+ia,i)
          ENDDO
       ENDDO
       !$omp do
       DO ia=1,ia_sum
          isa=isa0+ia
          DO iv=1,ngh
             DO ithread=1,nthreads
                deeq(isa,iv,j)=deeq(isa,iv,j)+temp(ia,iv,ithread)
             END DO
          END DO
       END DO
       IF(I.NE.J)THEN
          temp(:,:,methread)=0.0_real_8
          !$omp do
          DO ijv=1, ngh*(ngh+1)/2
             ijv_c=0
             count2: DO iv=1,ngh
                DO jv=iv,ngh
                   ijv_c=ijv_c+1
                   IF(ijv_c.EQ.ijv) EXIT count2
                END DO
             END DO count2
             fac_l=fac
             IF(iv.EQ.jv)fac_l=fac*0.5_real_8
             offset_fnl_iv=offset_fnl0+(iv-1)*ia_sum
             offset_fnl_jv=offset_fnl0+(jv-1)*ia_sum
             DO ia=1,ia_sum
                temp(ia,jv,methread)=temp(ia,jv,methread)+ylm(ia,ijv,ipot)*fac_l*&
                     fnl_p(offset_fnl_iv+ia,j)
                temp(ia,iv,methread)=temp(ia,iv,methread)+ylm(ia,ijv,ipot)*fac_l*&
                     fnl_p(offset_fnl_jv+ia,j)
             ENDDO
          ENDDO
          !$omp do
          DO ia=1,ia_sum
             isa=isa0+ia
             DO iv=1,ngh
                DO ithread=1,nthreads
                   deeq(isa,iv,i)=deeq(isa,iv,i)+temp(ia,iv,ithread)
                END DO
             END DO
          END DO
       END IF
    END DO
    !$omp end parallel
  END SUBROUTINE build_deeq_hfx
  
  SUBROUTINE prep_smallmem_newd(nst,num_pot,ig_start,nhg,na,deeq,f,fnl_p,fion,vpot,tfor,do_hfx,sym)
    REAL(REAL_8),INTENT(OUT)                 :: deeq(ions1%nat,maxsys%nhxs,*)
    COMPLEX(real_8),INTENT(IN)               :: vpot(ncpw%nhg)
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    INTEGER,INTENT(IN) __CONTIGUOUS          :: nst(:,:)
    INTEGER,INTENT(IN)                       :: num_pot, na(2,ions1%nsp), ig_start, nhg
    REAL(real_8),INTENT(IN)                  :: f(*), fnl_p(il_fnl_packed(1),*)
    LOGICAL,INTENT(IN)                       :: tfor, do_hfx
    LOGICAL,INTENT(IN) __CONTIGUOUS          :: sym(:)

    COMPLEX(real_8), ALLOCATABLE             :: vtmp(:)
    REAL(real_8), ALLOCATABLE                :: ylm(:), fnlt(:,:,:)
    INTEGER                                  :: isa0, offset_fnl0, is, ia_sum, ijv, iv, jv, &
                                                offset_iv, offset_jv, ig, ia, isa, i, k, &
                                                nhh, methread, offset_ylm, il_vtmp(1), &
                                                il_ylm(1), il_fnlt(3), ierr
    CHARACTER(*), PARAMETER                  :: procedureN='smallmem_newd'
    REAL(real_8)                             :: ftmp, otr
    REAL(real_8), EXTERNAL                   :: dotp
    LOGICAL                                  :: geq0_back
    ! ==--------------------------------------------------------------==

    il_ylm(1)=0
    IF(tfor)THEN
       DO is=1,ions1%nsp
          IF(pslo_com%tvan(is))THEN
             il_ylm(1)=MAX(il_ylm(1),ions0%na(is)*(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2)
          END IF
       END DO
    END IF
    il_vtmp(1)=nhg

    IF(tfor)THEN
       IF(.NOT.do_hfx)THEN
          il_fnlt(1)=nst(2,1)-nst(1,1)+1
          il_fnlt(2)=maxsys%nhxs
          il_fnlt(3)=parai%ncpus
       END IF
       ALLOCATE(fnlt(il_fnlt(1),il_fnlt(2),il_fnlt(3)), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlt)',&
            __LINE__,__FILE__)
    END IF

    ALLOCATE(qg(ncpw%nhg,1), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate qg)',&
         __LINE__,__FILE__)
    ALLOCATE(vtmp(il_vtmp(1)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate vtmp)',&
         __LINE__,__FILE__)
    ALLOCATE(ylm(il_ylm(1)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ylm)',&
         __LINE__,__FILE__)

    !dotp function is not cp_grp aware... temporary redefine geq0
    geq0_back=geq0
    geq0=(ig_start.EQ.1.AND.geq0)
    !
    isa0=0
    offset_fnl0=0
    DO is=1,ions1%nsp
       ia_sum=na(2,is)-na(1,is)+1
       IF (pslo_com%tvan(is)) THEN
          IF(tfor)THEN
             nhh=(nlps_com%ngh(is)*(nlps_com%ngh(is)+1))/2
             methread=1
             !$omp parallel private(methread)
             !$ methread=omp_get_thread_num()+1
             CALL calc_rho(nst(1,1),nst(2,1),nlps_com%ngh(is),ia_sum,offset_fnl0,nst(2,1)-nst(1,1)+1,&
                  na(1,is)-1,ions0%na(is),ylm,fnlt(:,:,methread),f,fnl_p,nhh)
             !$omp end parallel
             CALL mp_sum(ylm,nhh*ions0%na(is),parai%cp_grp)
          END IF

          ijv=0
          DO iv=1,nlps_com%ngh(is)
             DO jv=iv,nlps_com%ngh(is)
                ijv=ijv+1
                CALL qvan2(iv,jv,is,qg(:,1))
                !$omp parallel do private(IG) shared(QG) schedule(static)
                DO ig=ig_start,ig_start+nhg-1
                   qg(ig,1)=CONJG(qg(ig,1))*vpot(ig)
                ENDDO
                DO ia=1,ions0%na(is)
                   isa=isa0+ia
                   !$omp parallel do private(IG) shared(VTMP) schedule(static)
                   DO ig=ig_start,ig_start+nhg-1
                      vtmp(ig-ig_start+1)=ei1(isa,inyh(1,ig))*ei2(isa,inyh(2,ig))*&
                           ei3(isa,inyh(3,ig))
                   ENDDO
                   deeq(isa,iv,jv)=parm%omega*dotp(nhg,qg(ig_start,1),vtmp)
                   deeq(isa,jv,iv)=deeq(isa,iv,jv)
                   IF (tfor) THEN
                      offset_ylm=ions0%na(is)*(ijv-1)+ia
                      otr=parm%omega*parm%tpiba*ylm(offset_ylm)
                      DO k=1,3
                         CALL cftemp(nhg,qg(ig_start,1),vtmp,gk(1,ig_start),k,ftmp)
!                         IF (iv.NE.jv) ftmp=2.0_real_8*ftmp
                         fion(k,ia,is)=fion(k,ia,is)+ftmp*otr
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDIF
       offset_fnl0=offset_fnl0+(na(2,is)-na(1,is)+1)*nlps_com%ngh(is)
       isa0=isa0+ions0%na(is)
    ENDDO
    geq0=geq0_back
    DEALLOCATE(qg, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate qg)',&
         __LINE__,__FILE__)
    DEALLOCATE(vtmp, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate vtmp)',&
         __LINE__,__FILE__)
    DEALLOCATE(ylm, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ylm)',&
         __LINE__,__FILE__)
    IF(tfor)THEN
       DEALLOCATE(fnlt, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlt)',&
            __LINE__,__FILE__)
    END IF
    CALL mp_sum(deeq,ions1%nat*maxsys%nhxs*maxsys%nhxs,parai%cp_grp)
  END SUBROUTINE prep_smallmem_newd





  ! ==================================================================
END MODULE newd_utils

! ==================================================================
SUBROUTINE cftemp(nhg,qg,vtmp,gk,k,ftmp)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: nhg
  REAL(real_8)                               :: qg(2,nhg), vtmp(2,nhg), &
                                                gk(3,nhg)
  INTEGER                                    :: k
  REAL(real_8)                               :: ftmp

  INTEGER                                    :: ig

! Variables
! ==--------------------------------------------------------------==
! G=0 TERM IS ZERO !

  ftmp=0.0_real_8
  !$omp parallel do private(IG) reduction(+:FTMP) schedule(static)
  DO ig=1,nhg
     ftmp=ftmp+gk(k,ig)*(vtmp(1,ig)*qg(2,ig)-vtmp(2,ig)*qg(1,ig))
  ENDDO
  ftmp=2.0_real_8*ftmp
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE cftemp
