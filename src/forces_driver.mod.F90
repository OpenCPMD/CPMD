#include "cpmd_global.h"

MODULE forces_driver
  USE autotune_utils,                  ONLY: autotune
  USE benc,                            ONLY: ibench
  USE cp_grp_utils,                    ONLY: cp_grp_get_cp_rank,&
                                             cp_grp_get_sizes,&
                                             cp_grp_redist,&
                                             cp_grp_redist_array_f
  USE csize_utils,                     ONLY: csize
  USE dotp_utils,                      ONLY: dotp
  USE efield_utils,                    ONLY: extfield
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fft_maxfft,                      ONLY: maxfft
  USE fft,                             ONLY: batch_fft,&
                                             fft_tune_max_it
  USE fftprp_utils,                    ONLY: autotune_fftbatchsize
  USE fnonloc_utils,                   ONLY: fnonloc
  USE func,                            ONLY: func1
  USE geq0mod,                         ONLY: geq0
  USE gsize_utils,                     ONLY: gsize
  USE hfx_drivers,                     ONLY: hfx
  USE hfxmod,                          ONLY: hfxc3
  USE hnlmat_utils,                    ONLY: hnlmat
  USE hubbardu,                        ONLY: c2u0,hubbu
  USE hubbardu_utils,                  ONLY: add_hubbardu
  USE ions,                            ONLY: ions1
  USE jrotation_utils,                 ONLY: set_orbdist
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE kpts,                            ONLY: tkpts
  USE kpnt,                            ONLY: wk
  USE mp_interface,                    ONLY: mp_sum
  USE nlforce_utils,                   ONLY: nlforce
  USE nlps,                            ONLY: imagp
  USE nlsl_utils,                      ONLY: nlsl
  USE norm,                            ONLY: cnorm,&
                                             gemax,&
                                             gnmax,&
                                             gnorm
  USE nvtx_utils
  USE opeigr_utils,                    ONLY: give_scr_opeigr
  USE ovlap_utils,                     ONLY: ovlap,&
                                             ovlap2_dist
  USE parac,                           ONLY: parai,&
                                             paral
  USE printp_utils,                    ONLY: printp3
  USE pslo,                            ONLY: pslo_com
  USE puttau_utils,                    ONLY: taucl
  USE reigs_utils,                     ONLY: reigs
  USE reshaper,                        ONLY: reshape_inplace
  USE rgsvan_utils,                    ONLY: rgsvan
  USE rnlfl_utils,                     ONLY: rnlfl
  use rnlfor_utils,                    ONLY: rnlfor_hfx
  USE rnlsm_utils,                     ONLY: rnlsm
  USE ropt,                            ONLY: ropt_mod
  USE rotate_utils,                    ONLY: rotate,&
                                             rottr,&
                                             rotate_fnl,&
                                             rotate_c0_fnl
  USE rscpot_utils,                    ONLY: rscpot
  USE rswfmod,                         ONLY: rsactive
  USE sfac,                            ONLY: fnl_packed,&
                                             il_fnl_packed,&
                                             dfnl_packed
  USE store_types,                     ONLY: cprint
  USE spin,                            ONLY: clsd,&
                                             lspin2,&
                                             lspin3,&
                                             spin_mod
  USE summat_utils,                    ONLY: summat
  USE symtrz_utils,                    ONLY: symvec
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             maxsys,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             paraw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean,&
                                             zclean_k
  USE vpsi_utils,                      ONLY: vpsi,&
                                             vpsi_batchfft
!!use rotate_utils, only : rotate_c
!!use ovlap_utils, only : ovlap_c
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: forces
  PUBLIC :: gscal
  PUBLIC :: gscal_c

CONTAINS

  ! ==================================================================
  SUBROUTINE forces(c0,c2,tau0,fion,rhoe,psi,nstate,nkpoint,lproj,tfor,tcsize)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==                     THE TOTAL ENERGY                         ==
    ! ==                  THE ELECTRONIC FORCES                       ==
    ! ==                  THE FORCES ON THE IONS                      ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8),INTENT(IN), TARGET&
                               __CONTIGUOUS  :: c0(:,:,:)
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: c2(:,:)
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: tau0(:,:,:)
    REAL(real_8),INTENT(OUT) __CONTIGUOUS    :: fion(:,:,:), &
                                                rhoe(:,:)
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: psi(:,:)
    INTEGER,INTENT(IN)                       :: nstate, nkpoint
    LOGICAL,INTENT(IN)                       :: lproj, tfor, tcsize

    CHARACTER(*), PARAMETER                  :: procedureN = 'forces'

    CHARACTER(len=30)                        :: tag
    COMPLEX(real_8)                          :: zee
    COMPLEX(real_8), ALLOCATABLE             :: psiab(:)
#ifdef _USE_SCRATCHLIBRARY
    COMPLEX(real_8), ALLOCATABLE, TARGET     :: auxc(:)
#else
    COMPLEX(real_8), ALLOCATABLE, TARGET     :: auxc(:)
#endif
    COMPLEX(real_8), EXTERNAL                :: zdotc
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: cgam(:), c0_ptr(:,:,:)
    INTEGER :: first_g, i, ierr, ik, ikind, il_auxc, il_fsc,  &
      il_psiab, il_scrdip, ipp, isub, isub2, isub3, j, jj, last_g,  &
      naa, ngw_l, nstx, NSTX_grp, is
    INTEGER(int_8)                           :: il_c0_ort(3), il_smat(2), il_gam(2)
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: NWA12_grp
    LOGICAL                                  :: debug, redist_c2, &
                                                case1, case2, case3, case4
    REAL(real_8)                             :: ee, ehfx, fi, fj, vhfx
    REAL(real_8), ALLOCATABLE                :: a1mat(:,:), a2mat(:,:), &
                                                a3mat(:,:), fsc(:), &
                                                scrdip(:)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: smat(:,:), fnlgam_packed(:,:), gam(:,:)
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: c0_ort(:,:,:)
#else
    REAL(real_8), ALLOCATABLE                :: smat(:,:), fnlgam_packed(:,:)
    REAL(real_8), ALLOCATABLE, TARGET        :: gam(:,:)
    COMPLEX(real_8), ALLOCATABLE, TARGET     :: c0_ort(:,:,:)
#endif
    REAL(real_8), EXTERNAL                   :: dasum
    REAL(real_8), POINTER                    :: aux(:)

    INTEGER, SAVE                            :: autotune_it=1
    CALL tiset(procedureN,isub)
    CALL tiset(procedureN//'_a',isub2)
    !TK noforce goes here:
    !differces to forces:
    !orthogonolize c0 in the beginning
    !rotate c2 into nonorthogonal basis in the end
    !lrpoj is ignored == always .true.
    !for pure ncpp nonort is not different at all (smat.EQ.I)
    IF (tkpts%tkpnt.AND.pslo_com%tivan) &
         CALL stopgm(procedureN,'K-POINTS NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    !call rnlsm1 and rnlsm2 in the kpt loop
    case1=.FALSE.
    !call rnlsm1 for rgsvan and orthogonalization and only rnlsm2 in the kpt loop
    case2=.FALSE.
    !call rnlsm1 in the kpt loop and rnlsm2 just before rnlfl
    case3=.FALSE.
    !call rnlsm1 for rgsvan and orthogonalization and only
    !rnlsm2 in the kpt loop but just before rnlfl
    case4=.FALSE.
    IF(.NOT.pslo_com%tivan)THEN
       case1=.TRUE.
    ELSEIF(pslo_com%tivan)THEN
       IF(cntl%nonort.OR..NOT.cntl%tmdcp)THEN
          IF(pslo_com%mixed_psp)THEN
             case2=.TRUE.
          ELSE
             case4=.TRUE.
          END IF
       ELSE
          IF(pslo_com%mixed_psp)THEN
             case1=.TRUE.
          ELSE
             case3=.TRUE.
          END IF
       END IF
    END IF
    IF(case4.OR.case2)THEN
       IF (lspin2%tlse) CALL stopgm('NOFORCE','NO LSE ALLOWED HERE',&
            __LINE__,__FILE__)
       IF (imagp.EQ.2) CALL stopgm('NOFORCE','K-POINT NOT IMPLEMENTED',&
            __LINE__,__FILE__)
       il_c0_ort(1)=ncpw%ngw
       il_c0_ort(2)=nstate
       il_c0_ort(3)=1
       il_smat(1)=nstate
       il_smat(2)=nstate
       !In case of CP-MD we need to preserve the non-orthogonal orbitals to propagate them properly
       IF(cntl%tmdcp)THEN
#ifdef _USE_SCRATCHLIBRARY
          CALL request_scratch(il_c0_ort,c0_ort,procedureN//'_c0_ort',ierr)
#else
          ALLOCATE(c0_ort(il_c0_ort(1),il_c0_ort(2),il_c0_ort(3)),stat=ierr)
#endif
          IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
               __LINE__,__FILE__)
          CALL dcopy(2*ncpw%ngw*nstate,c0,1,c0_ort,1)
          c0_ptr=> c0_ort
       ELSE
          c0_ptr=>c0
       END IF
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_smat,smat,procedureN//'_smat',ierr)
#else
       ALLOCATE(smat(il_smat(1),il_smat(2)),stat=ierr)
#endif
       IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
            __LINE__,__FILE__)

       IF(parai%cp_nogrp.GT.1) CALL cp_grp_redist_array_f(c0_ptr,ncpw%ngw,nstate)
       IF(pslo_com%tivan) CALL rnlsm(c0_ptr(:,:,1),nstate,1,1,.FALSE.,&
            unpack_dfnl_fnl=.FALSE.)
       CALL rgsvan(c0_ptr(:,:,1),nstate,smat,store_nonort=cntl%tmdcp)
    ELSE
       c0_ptr=>c0
    END IF

    !autotuning of rhoofr/vpsi/rnlsm
    IF(cntl%rnlsm_autotune.OR.cntl%fft_tune_batchsize)THEN
       IF(ibench(3).EQ.1) CALL autotune(c0_ptr(:,:,1),c2,rhoe,psi,nstate)
    END IF


    gnmax=0.0_real_8
    gnorm=0.0_real_8
    ehfx =0.0_real_8
    vhfx =0.0_real_8
    ! ==--------------------------------------------------------------==
    redist_c2=.FALSE.
    IF (func1%mhfx.EQ.1.AND..NOT.hfxc3%use_new_hfx) redist_c2=.TRUE.
    IF (pslo_com%tivan.OR.lspin2%tlse) redist_c2=.TRUE.
    IF (cntl%ttau) redist_c2=.TRUE.
    ! ==--------------------------------------------------------------==
    debug=.FALSE.
    il_gam(1) = imagp*nstate
    il_gam(2) = nstate
    il_auxc = 1
    IF (cntl%tfield) CALL give_scr_opeigr(il_scrdip,tag,nstate)
    IF (cntl%tdmal) THEN
       ALLOCATE(NWA12_grp(0:parai%cp_nproc-1,2),stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
            __LINE__,__FILE__)
       CALL set_orbdist(nstate,cnti%nstblk,parai%cp_nproc,NSTX_grp)
       NWA12_grp(0:parai%cp_nproc-1,1:2)=paraw%nwa12(0:parai%cp_nproc-1,1:2)
       CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
       il_gam=1
    ELSE
       il_auxc=MAX(il_auxc,nstate*nstate)! AUXC space for OVLAP (Laio A.)
    ENDIF
    il_fsc=nstate
    il_psiab=2
    IF (lspin2%tlse) il_psiab=2*maxfft
    ! ==--------------------------------------------------------------==
    ! ==   CALCULATE THE POTENTIAL AND THE FORCE ON THE IONS          ==
    ! ==--------------------------------------------------------------==
    DO ikind=1,nkpoint
       IF (nstate >0 ) THEN
          IF(case1)THEN
             CALL rnlsm(c0_ptr(:,:,ikind),nstate,1,ikind,tfor)
          ELSEIF(case2)THEN
             CALL rnlsm(c0_ptr(:,:,ikind),nstate,1,ikind,tfor,only_dfnl=.TRUE.)
          ELSEIF(case3)THEN
             CALL rnlsm(c0_ptr(:,:,ikind),nstate,1,ikind,.FALSE.,&
                  unpack_dfnl_fnl=ropt_mod%calste)
          END IF
       ENDIF
    ENDDO
    rsactive = cntl%krwfn
    IF(cntl%fft_tune_batchsize.AND.batch_fft.AND.(ibench(3).EQ.0))THEN
       IF(autotune_it.LE.fft_tune_max_it+1)THEN
          autotune_it=autotune_it+1
          CALL autotune_fftbatchsize()
       END IF
    END IF
    CALL rscpot(c0_ptr(:,:,1),tau0,fion,rhoe,psi,tfor,ropt_mod%calste,nstate,nkpoint)
    CALL tihalt(procedureN//'_a',isub2)
    ! ==--------------------------------------------------------------==
    ! ==   CALCULATE THE ELECTRONIC FORCE                             ==
    ! ==--------------------------------------------------------------==
    IF (nstate.EQ.0) THEN
       gemax=0._real_8
       cnorm=0._real_8
       GOTO 2000
    ENDIF
    CALL zeroing(c2)!,nkpt%ngwk*nstate)
    ! EHR[
    ! cntl%tmdeh only needs forces acting on the nuclei
    ! skip calculation of C2 (forces on KS orbitals)
    IF (cntl%tmdeh) GOTO 2000
    ! EHR]
    ! ==--------------------------------------------------------------==
    ! == COMPUTE THE FORCE ON THE ELECTRONIC DEGREES OF FREEDOM DUE   ==
    ! == TO THE LOCAL POTENTIAL (STORED IN RHOE)                      ==
    ! ==--------------------------------------------------------------==
    IF (cntl%tdmal) THEN
       ALLOCATE(a1mat(nstate, nstx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(a2mat(nstate, nstx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       ALLOCATE(a3mat(nstate,nstx),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_gam,gam,procedureN//'_gam',ierr)
#else
    ALLOCATE(gam(il_gam(1), il_gam(2)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF(tkpts%tkpnt) CALL reshape_inplace(gam, (/nstate*nstate/), cgam)

    IF(.NOT.pslo_com%tivan)THEN
       ALLOCATE(fsc(il_fsc),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    END IF
    IF (lspin2%tlse) THEN
       ALLOCATE(auxc(il_auxc),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL reshape_inplace(auxc, (/2*il_auxc/), aux)

       ALLOCATE(psiab(il_psiab),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    END IF

    IF (cntl%tfield) THEN
       ALLOCATE(scrdip(il_scrdip),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    ENDIF

    DO ik=1,nkpoint
       IF(batch_fft.AND..NOT.tkpts%tkpnt)THEN
          CALL vpsi_batchfft(c0_ptr(:,:,ik),c2,crge%f(:,1),rhoe,psi(:,1),nstate,ik,&
               clsd%nlsd,redist_c2)
       ELSE
          CALL vpsi(c0_ptr(:,:,ik),c2,crge%f(:,1),rhoe,psi(:,1),nstate,ik,clsd%nlsd,&
               redist_c2)
       END IF
       ! c2u0 is calculated in uprho or rscpot
       IF (hubbu%debug) THEN
            IF (paral%io_parent) write(6,*) procedureN,"| starting add_hubbardu"
       ENDIF
       IF (cntl%thubb) CALL add_hubbardu(c2,c2u0,nstate)
       CALL hfx(c0_ptr(:,:,ik),c2,crge%f(:,1),psi(:,1),nstate,ehfx,vhfx,redist_c2,fion,tfor)
       
       ! McB    IF (cntl%tfield) CALL EFIELD(TAU0,FION,C0_PTR,C2,SCRDIP,IL_SCRDIP,NSTATE)
       IF (cntl%tfield) THEN
          CALL extfield(tau0,fion,c0_ptr,c2,nstate)
       ENDIF
       ener_com%exc=ener_com%exc+ehfx
       ener_com%etot=ener_com%etot+ehfx
       CALL tiset(procedureN//'_b',isub3)
       __NVTX_TIMER_START ( procedureN//'_b' )

       IF (pslo_com%tivan) THEN
          CALL ovlap(nstate,gam,c2,c0_ptr(:,:,ik),redist=.FALSE.,full=.FALSE.)
          CALL hnlmat(gam,crge%f,nstate)
          CALL summat(gam,nstate,lsd=.TRUE.,gid=parai%cp_grp,symmetrization=&
               ropt_mod%prteig.OR.ropt_mod%calste)
#ifdef _USE_SCRATCHLIBRARY
          CALL request_scratch(il_fnl_packed,fnlgam_packed,procedureN//'_fnlgam_packed',ierr)
#else
          ALLOCATE(fnlgam_packed(il_fnl_packed(1),il_fnl_packed(2)),STAT=ierr)
#endif
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          IF(cntl%distribute_fnl_rot)THEN
             CALL rotate_c0_fnl(ncpw%ngw,c0_ptr(:,:,ik),c2,int(il_fnl_packed(1)),fnl_packed,&
                  fnlgam_packed,nstate,gam,redist=.NOT.cntl%nonort)
          ELSE
             CALL rotate(-1.0_real_8,c0_ptr(:,:,ik),1.0_real_8,c2,gam,&
                  nstate,2*nkpt%ngwk,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown,symm=.TRUE.,&
                  use_cp=.TRUE.,redist=.NOT.cntl%nonort)
             CALL rotate_fnl(int(il_fnl_packed(1)),fnl_packed,fnlgam_packed,nstate,gam)
          END IF
          CALL nlforce(c2,crge%f,fnl_packed,fnlgam_packed,nstate,redist=.NOT.cntl%nonort)
          IF (tfor) THEN
             CALL rnlsm(c0_ptr(:,:,ik),nstate,1,ik,tfor,unpack_dfnl_fnl=.FALSE.,&
                  only_dfnl=.TRUE.)
             CALL rnlfl(fion,nstate,nkpoint,fnl_packed,fnlgam_packed,dfnl_packed)
             IF (func1%mhfx.EQ.1.) call rnlfor_hfx(fion,crge%f,wk,nstate,ik,dfnl_packed)
          END IF
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_fnl_packed,fnlgam_packed,procedureN//'_fnlgam_packed',ierr)
#else
          DEALLOCATE(fnlgam_packed, STAT=ierr)
#endif
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
          IF (ropt_mod%calste) CALL nlsl(gam,nstate)
          IF (ropt_mod%prteig.AND.paral%io_parent) THEN
             CALL dscal(nstate*nstate,-1.0_real_8,gam,1)
             CALL reigs(nstate,gam,crge%f)
             CALL dscal(nstate*nstate,-1.0_real_8,gam,1)
          ENDIF
       ELSEIF (lspin2%tlse) THEN
          CALL fnonloc(c2,crge%f,nstate,ik,clsd%nlsd,.TRUE.)
          CALL ovlap(nstate,aux,c0_ptr(:,:,ik),c2)
          IF (lspin2%tros) THEN
             CALL addup2(nstate,gam,aux)
          ELSE
             CALL addup(nstate,gam,aux)
          ENDIF
          IF (lspin2%tmgoe.OR.lspin2%tros) THEN
             CALL mp_sum(gam,nstate*nstate,parai%allgrp)
          ELSE
             CALL summat(gam,nstate)
          ENDIF
          lspin3%hablse=gam(clsd%ialpha,clsd%ibeta)
          lspin3%haolse=dasum(clsd%ialpha-1,gam(1,clsd%ialpha),1)
          lspin3%hbolse=dasum(clsd%ialpha-1,gam(1,clsd%ibeta),1)
          lspin3%rotab=0.5_real_8*ATAN(2._real_8*lspin3%hablse/&
               (gam(clsd%ialpha,clsd%ialpha)-gam(clsd%ibeta,clsd%ibeta)))
          IF (lproj.OR.cntl%nonort)&
               CALL rotate(-0.5_real_8,c0_ptr(:,:,ik),1._real_8,c2,gam,&
               nstate,2*nkpt%ngwk,.FALSE.,spin_mod%nsup,spin_mod%nsdown)
       ELSE
          CALL fnonloc(c2,crge%f,nstate,ik,clsd%nlsd,redist_c2)
          IF (.NOT.redist_c2) CALL cp_grp_redist(c2,nkpt%ngwk,nstate)
          IF (lproj.OR.cntl%nonort) THEN
             DO i=1,nstate
                IF (crge%f(i,ik).EQ.0._real_8) THEN
                   fsc(i)=0._real_8
                ELSE
                   fsc(i)=1._real_8
                ENDIF
             ENDDO
             IF (cnti%iproj.EQ.0) THEN
             ELSEIF (cnti%iproj.EQ.1) THEN
                IF (tkpts%tkpnt) THEN
                   IF (nkpt%ngwk.GT.0) THEN
                      DO i=1,nstate
                         zee=fsc(i)*zdotc(nkpt%ngwk,c0_ptr(1,i,ik),1,c2(1,1),1)
                         CALL mp_sum(zee,parai%allgrp)
                         CALL zaxpy(nkpt%ngwk,-zee,c0_ptr(1,i,ik),1,c2(1,1),1)
                      ENDDO
                   ENDIF
                ELSE
                   IF (ncpw%ngw.GT.0) THEN
                      DO i=1,nstate
                         ee=fsc(i)*dotp(ncpw%ngw,c0_ptr(:,i,ik),c2(:,1))
                         CALL mp_sum(ee,parai%allgrp)
                         CALL daxpy(2*ncpw%ngw,-ee,c0_ptr(1,i,ik),1,c2(1,1),1)
                      ENDDO
                   ENDIF
                ENDIF
             ELSEIF (cnti%iproj.EQ.2.OR.cntl%nonort) THEN
                IF (tkpts%tkpnt) THEN
                   IF (cntl%tlsd) THEN
                      naa=spin_mod%nsup+1
                      CALL zeroing(gam)!,2*nstate*nstate)
                      CALL ovlap_c(spin_mod%nsup,gam,c2,c0_ptr(:,:,ik))
                      CALL ovlap_c(spin_mod%nsdown,gam(2*naa,naa),&
                           c2(1,naa),c0_ptr(1,naa,ik))
                   ELSE
                      CALL ovlap_c(nstate,gam,c2,c0_ptr(:,:,ik))
                   ENDIF
                   ! AK FIXME: can this work? GAM is complex.
                   CALL mp_sum(gam,2*nstate*nstate,parai%allgrp)
                   CALL gscal_c(nstate,cgam,fsc)
                   CALL rotate_c(CMPLX(-1.0_real_8,0._real_8,kind=real_8),c0_ptr(:,:,ik),&
                        CMPLX( 1.0_real_8,0._real_8,kind=real_8),c2,gam,nstate)
                ELSE
                   IF (cntl%tdmal) THEN

                      CALL ovlap2_dist( nstate, a1mat, c2, c0_ptr(:,:,ik) )

                      ipp = cp_grp_get_cp_rank(parai%me,parai%cp_inter_me)
                      IF (parai%cp_me.NE.ipp) CALL stopgm(procedureN,&
                           ' something wrong here',&
                           __LINE__,__FILE__)
                      DO i=1,nstate
                         fi= fsc(i)
                         DO j=NWA12_grp(parai%cp_me,1),nwa12_grp(parai%cp_me,2)
                            jj=j-NWA12_grp(parai%cp_me,1)+1
                            fj= fsc(j)
                            a1mat(i,jj)=(fi*fj)*a1mat(i,jj)
                            ! >>>>
                            ! SDF the overlap above is not taking care of the
                            ! SDF orthogonality of the spins. So we do it here...
                            ! SDF some speedup could be gained in ovlap2 and rotate_da.
                            IF (cntl%tlsd) THEN
                               IF (j.GT.spin_mod%nsup.AND.i.LE.spin_mod%nsup) &
                                    a1mat(i,jj)=0.0_real_8
                               IF (j.LE.spin_mod%nsup.AND.i.GT.spin_mod%nsup) &
                                    a1mat(i,jj)=0.0_real_8
                            ENDIF
                            ! <<<
                         ENDDO
                      ENDDO

                      CALL cp_grp_get_sizes(ngw_l=ngw_l,first_g=first_g,&
                           last_g=last_g)

                      CALL rotate_da(-1._real_8,c0_ptr(first_g,1,ik),1._real_8,&
                           c2(first_g,1),a1mat,2*ncpw%ngw,2*ngw_l,&
                           nstate,NWA12_grp(0,1),nwa12_grp(0,2),&
                           NSTX_grp,parai%cp_me,parap%pgroup,parai%cp_nproc,parai%cp_grp,&
                           cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)

                      CALL cp_grp_zero_g(c2,nkpt%ngwk,nstate,first_g,last_g)
                      CALL cp_grp_redist(c2,nkpt%ngwk,nstate)

                   ELSE
                      CALL ovlap(nstate,gam,c2,c0_ptr(:,:,ik))
                      CALL mp_sum(gam,nstate*nstate,parai%allgrp)
                      CALL gscal(nstate,gam,fsc)
                      CALL rotate(-1.0_real_8,c0_ptr(:,:,ik),1.0_real_8,c2,gam,&
                           nstate,2*nkpt%ngwk,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
       ENDIF
       IF(case4.or.case2) THEN
          ! ==--------------------------------------------------------------==
          ! ==   ROTATE ELECTRONIC FORCE BACK INTO NONORTHOGONAL BASIS      ==
          ! ==--------------------------------------------------------------==
          CALL rottr(1._real_8,c2,smat,"T",nstate,ncpw%ngw,cntl%tlsd,spin_mod%nsup,&
               spin_mod%nsdown,redist=.NOT.cntl%nonort,use_cp=cntl%nonort)
          IF (geq0) CALL zclean(c2,nstate,ncpw%ngw)
       ELSE
          IF (geq0) THEN
             IF (tkpts%tkpnt) THEN
                CALL zclean_k(c2,nstate,ncpw%ngw)
             ELSE
                CALL zclean(c2,nstate,ncpw%ngw)
             ENDIF
          ENDIF
       END IF

       IF(tcsize) CALL csize(c2,nstate,gemax,cnorm,use_cp_grps=cntl%nonort,special=cntl%nonort)

       __NVTX_TIMER_STOP
       CALL tihalt(procedureN//'_b',isub3)
    ENDDO

    IF (cntl%tdmal) THEN
       DEALLOCATE(a1mat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(a2mat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(a3mat,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_gam,gam,procedureN//'_gam',ierr)
#else
    DEALLOCATE(gam,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'cannot deallocate gam', &
         __LINE__,__FILE__)
    IF(.NOT.pslo_com%tivan)THEN
       DEALLOCATE(fsc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    END IF
    IF(lspin2%tlse)THEN
       DEALLOCATE(auxc,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
       DEALLOCATE(psiab,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    END IF
    IF (cntl%tfield) THEN
       DEALLOCATE(scrdip,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDIF
    IF(case4)THEN
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_smat,smat,procedureN//'_smat',ierr)
#else
       DEALLOCATE(smat,stat=ierr)
#endif
       IF(cntl%tmdcp)THEN
          IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',&
               __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_c0_ort,c0_ort,procedureN//'_c0_ort',ierr)
#else
          DEALLOCATE(c0_ort,stat=ierr)
#endif
          IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',&
               __LINE__,__FILE__)
       END IF
    END IF

2000 CONTINUE
    rsactive = .FALSE.
    IF (tfor) THEN
       CALL mp_sum(fion,3*maxsys%nax*maxsys%nsx,parai%allgrp)
       IF (paral%parent) THEN
          CALL symvec(fion)
          IF (cprint%twritefixforcetrajectory) CALL printp3(tau0,fion)
          CALL taucl(fion)
          CALL gsize(fion,gnmax,gnorm)
       ENDIF
    ENDIF
    !
    IF (ALLOCATED(NWA12_grp)) THEN
       DEALLOCATE(NWA12_grp,stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
    !
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==

  END SUBROUTINE forces
  ! ==================================================================
  SUBROUTINE addup(nstate,a1,a2)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: a1(nstate,*), a2(nstate,*)

    INTEGER                                  :: i, j

! ==--------------------------------------------------------------==

    DO i=1,nstate
       DO j=i,nstate
          a1(i,j)=a2(i,j)+a2(j,i)
          a1(j,i)=a1(i,j)
       ENDDO
    ENDDO
    IF (lspin2%tmgoe) THEN
       a1(clsd%ibeta,clsd%ialpha) = 2._real_8*((1._real_8-lspin3%mgmab)*a2(clsd%ibeta,clsd%ialpha)&
            +lspin3%mgmab *A2(clsd%ialpha,clsd%ibeta))
       a1(clsd%ialpha,clsd%ibeta) = 2._real_8*((1._real_8-lspin3%mgmba)*a2(clsd%ialpha,clsd%ibeta)&
            +lspin3%mgmba *A2(clsd%ibeta,clsd%ialpha))
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE addup
  ! ==================================================================
  SUBROUTINE addup2(nstate,a1,a2)
    ! NEW VERSION: MODIFIED GOEDECKER
    INTEGER                                  :: nstate
    REAL(real_8)                             :: a1(nstate,*), a2(nstate,*)

    INTEGER                                  :: i, j

! ==--------------------------------------------------------------==

    DO i=1,nstate-2
       a1(i,i) = 2.0_real_8 * a2(i,i)
       DO j=i+1,nstate-2
          a1(i,j) = a2(i,j) + a2(j,i)
          a1(j,i) = a2(j,i) + a2(i,j)
       ENDDO
       a1(nstate-1,i) = 2.0_real_8*(lspin3%mgab(1)*a2(nstate-1,i)&
            +lspin3%mgab(2)*A2(I,NSTATE-1))
       a1(nstate,i) = 2.0_real_8*(lspin3%mgab(3)*a2(nstate,i)&
            +lspin3%mgab(4)*A2(I,NSTATE))
       a1(i,nstate-1) = 2.0_real_8*(lspin3%mgab(5)*a2(i,nstate-1)&
            +lspin3%mgab(6)*A2(NSTATE-1,I))
       a1(i,nstate) = 2.0_real_8*(lspin3%mgab(7)*a2(i,nstate)&
            +lspin3%mgab(8)*A2(NSTATE,I))
    ENDDO
    a1(nstate,nstate-1) = 2.0_real_8*(lspin3%mgab(9)*a2(nstate,nstate-1)&
         +lspin3%mgab(10)*A2(NSTATE-1,NSTATE))
    a1(nstate-1,nstate) = 2.0_real_8*(lspin3%mgab(11)*a2(nstate-1,nstate)&
         +lspin3%mgab(12)*A2(NSTATE,NSTATE-1))
    a1(nstate-1,nstate-1) = 2.0_real_8 * a2(nstate-1,nstate-1)
    a1(nstate,nstate) = 2.0_real_8 * a2(nstate,nstate)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE addup2
  ! ==================================================================
  SUBROUTINE gscal(nstate,gam,f)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: gam(nstate,nstate), f(nstate)

    INTEGER                                  :: i, j

! ==--------------------------------------------------------------==

    !$omp parallel do private(I,J)
    DO j=1,nstate
       DO i=1,nstate
          gam(i,j)=f(i)*f(j)*gam(i,j)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gscal
  ! ==================================================================
  SUBROUTINE gscal_c(nstate,gam,f)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: gam(nstate,nstate)
    REAL(real_8)                             :: f(nstate)

    INTEGER                                  :: i, j

! ==--------------------------------------------------------------==

    !$omp parallel do private(I,J)
    DO j=1,nstate
       DO i=1,nstate
          gam(i,j)=f(i)*f(j)*gam(i,j)
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE gscal_c
  ! ==================================================================
  

END MODULE forces_driver
