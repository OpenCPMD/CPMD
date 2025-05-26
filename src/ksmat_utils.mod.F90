#include "cpmd_global.h"

MODULE ksmat_utils
  USE atwf,                            ONLY: atwf_mod,&
                                             atwp,&
                                             catom,&
                                             xxmat,&
                                             zxmat
  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnl_set,&
                                             fnlalloc,&
                                             fnldealloc
  USE fnonloc_utils,                   ONLY: fnonloc
  USE hnlmat_utils,                    ONLY: hnlmat
  USE gpu
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE kpts,                            ONLY: tkpts
  USE nlforce_utils,                   ONLY: nlforce
  USE ovlap_utils,                     ONLY: ovlap
  USE pslo,                            ONLY: pslo_com
  USE parac,                           ONLY: parai
  USE rnlsm_utils,                     ONLY: rnlsm
  USE rotate_utils,                    ONLY: rotate, &
                                             rotate_fnl
  USE sfac,                            ONLY: il_fnl_packed,&
                                             fnl,&
                                             fnl_packed
  USE summat_utils,                    ONLY: summat
  USE spin,                            ONLY: lspin2
  USE system,                          ONLY: cntl,&
                                             nkpt,&
                                             maxsys
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE vpsi_utils,                      ONLY: vpsi
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ksmat

CONTAINS

  ! ==================================================================
  SUBROUTINE ksmat(c2,vpot,psi,nstate,ikind)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE KOHN-SHAM MATRIX IN THE ATOMIC WAVEFUNCTION BASIS       ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: c2(:,:)
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: vpot(:,:)
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: psi(:)
    INTEGER,INTENT(IN)                       :: nstate, ikind

    COMPLEX(real_8)                          :: pab(1)
    INTEGER                                  :: ia, is, ierr, &
                                                ist, isub, natst
    INTEGER(int_8)                           :: il_gam(2), il_fnl_save(3), &
                                                il_fnlgam_packed(2)
    LOGICAL                                  :: tfdist2
    REAL(real_8), DIMENSION(20)              :: foc = 1.0_real_8
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: gam(:,:), fnlp_save(:,:), &
                                                fnl_save(:,:,:),&
                                                fnlgam_packed(:,:)
#else
    REAL(real_8), ALLOCATABLE                :: gam(:,:), fnlp_save(:,:), &
                                               fnl_save(:,:,:), &
                                               fnlgam_packed(:,:)
#endif
    CHARACTER(*), PARAMETER                  :: procedureN = 'ksmat'
! ==--------------------------------------------------------------==

    IF (lspin2%tlse) CALL stopgm('KSMAT','NO LSE ALLOWED HERE',& 
         __LINE__,__FILE__)
    ! AK      CALL zeroing(PAB)!,NNR1)
    pab(1)=0.0_real_8
    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    tfdist2=cntl%tfdist
    cntl%tfdist=.FALSE.
    CALL fnl_set('SAVE')
    CALL fnlalloc(atwp%nattot,.FALSE.,.FALSE.)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    update_first_to_gpu  =.FALSE.
    update_second_to_gpu =.FALSE.
    update_third_to_gpu  =.FALSE.
    update_result_to_host=.FALSE.
    !$omp target update to(catom)
    !$omp target enter data map(alloc:xxmat)
    !$omp target update to(xxmat)
#endif
    CALL rnlsm(catom,atwp%nattot,1,ikind,.FALSE.,unpack_dfnl_fnl=.NOT.pslo_com%tivan)
    IF(pslo_com%tivan)THEN
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_fnl_packed,fnlp_save,procedureN//'_fnlp_save',ierr)
#else
       ALLOCATE(fnlp_save(int(il_fnl_packed(1)),int(il_fnl_packed(2))),STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)

       CALL cpmd_dcopy(PRODUCT(INT(il_fnl_packed)),fnl_packed,1,fnlp_save,1)
    ELSE
       il_fnl_save(1)=ions1%nat
       il_fnl_save(2)=maxsys%nhxs
       il_fnl_save(3)=atwp%nattot
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_fnl_save,fnl_save,procedureN//'_fnl_save',ierr)
#else
       ALLOCATE(fnl_save(il_fnl_save(1),il_fnl_save(2),il_fnl_save(3)),STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)

       CALL dcopy(ions1%nat*maxsys%nhxs*atwp%nattot,fnl,1,fnl_save,1)
    END IF
    ist=1
    DO is=1,ions1%nsp
       natst=atwf_mod%numaor(is)
       il_gam(1)=natst
       il_gam(2)=natst
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_gam,gam,procedureN//'_gam',ierr)
#else
       ALLOCATE(gam(il_gam(1),il_gam(2)),STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)

       DO ia=1,ions0%na(is)
          IF(pslo_com%tivan)THEN
             CALL cpmd_dcopy(il_fnl_packed(1)*natst,fnlp_save(1,ist),1,fnl_packed,1)
          ELSE
             CALL dcopy(ions1%nat*maxsys%nhxs*natst,fnl_save(1,1,ist),1,fnl,1)
          END IF
          CALL zeroing(c2(:,1:natst))!,nkpt%ngwk*natst)
          CALL vpsi(catom(:,ist:ist+natst-1),c2,foc,vpot,psi,natst,ikind,1,.TRUE.)
          IF(pslo_com%tivan)THEN
             CALL ovlap(natst,gam,c2,catom(:,ist:ist+natst-1),redist=.FALSE.,full=.FALSE.)
             CALL hnlmat(gam,foc,natst)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             comm_buffers_on_host =.FALSE.
#endif
             CALL summat(gam,natst,lsd=.TRUE.,gid=parai%cp_grp,symmetrization=.FALSE.)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             comm_buffers_on_host =.TRUE.
#endif
             il_fnlgam_packed(1)=il_fnl_packed(1)
             il_fnlgam_packed(2)=natst
#ifdef _USE_SCRATCHLIBRARY
             CALL request_scratch(il_fnlgam_packed,fnlgam_packed,procedureN//'_fnlgam_packed',ierr)
#else
             ALLOCATE(fnlgam_packed(il_fnlgam_packed(1),il_fnlgam_packed(2)),STAT=ierr)
#endif
             IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
                  __LINE__,__FILE__)
             CALL rotate(-1.0_real_8,catom(:,ist:ist+natst-1),1.0_real_8,c2,gam,&
                  natst,2*nkpt%ngwk,.FALSE.,1,0,symm=.TRUE.,&
                  use_cp=.TRUE.,redist=.TRUE.)
             CALL rotate_fnl(int(il_fnl_packed(1)),fnl_packed,fnlgam_packed,natst,gam)
             CALL nlforce(c2,foc,fnl_packed,fnlgam_packed,natst,redist=.TRUE.)
#ifdef _USE_SCRATCHLIBRARY
             CALL free_scratch(il_fnlgam_packed,fnlgam_packed,procedureN//'_fnlgam_packed',ierr)
#else
             DEALLOCATE(fnlgam_packed, STAT=ierr)
#endif
             IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
                  __LINE__,__FILE__)
          ELSE
             CALL fnonloc(c2,foc,natst,ikind,1,.TRUE.)
          END IF
          IF (tkpts%tkpnt) THEN
             CALL ovlap2_c(nkpt%ngwk,atwp%nattot,natst,zxmat(1,ist),catom,c2)
          ELSE
             CALL ovlap2(nkpt%ngwk,atwp%nattot,natst,xxmat(1,ist),catom,c2,.TRUE.)
          ENDIF
          ist=ist+natst
       ENDDO
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_gam,gam,procedureN//'_gam',ierr)
#else
       DEALLOCATE(gam,STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    ENDDO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target update from(xxmat)
    !$omp target exit data map(delete:xxmat)
    update_first_to_gpu  =.true.
    update_second_to_gpu =.true.
    update_third_to_gpu  =.true.
    update_result_to_host=.true.
#endif

    IF(pslo_com%tivan)THEN
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_fnl_packed,fnlp_save,procedureN//'_fnlp_save',ierr)
#else
       DEALLOCATE(fnlp_save,STAT=ierr)
#endif
    ELSE
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_fnl_save,fnl_save,procedureN//'_fnl_save',ierr)
#else
       DEALLOCATE(fnl_save,STAT=ierr)
#endif
    END IF
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL fnldealloc(.FALSE.,.FALSE.)
    CALL fnl_set('RECV')
    cntl%tfdist=tfdist2
    IF (tkpts%tkpnt) THEN
       CALL dscal(2*atwp%nattot*atwp%nattot,-1._real_8,zxmat(1,1),1)
    ELSE
       CALL dscal(  atwp%nattot*atwp%nattot,-1._real_8,xxmat(1,1),1)
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE ksmat
  ! ==================================================================

END MODULE ksmat_utils
