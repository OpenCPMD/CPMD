#include "cpmd_global.h"

MODULE rgsvan_utils
  USE csmat_utils,                     ONLY: csmat
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms,&
                                             cp_grp_redist_dfnl_fnl
  USE error_handling,                  ONLY: stopgm
  USE gpu
  USE geq0mod,                         ONLY: geq0
  USE fnl_utils,                       ONLY: unpack_fnl
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8,&
                                             int_8,&
                                             int_4
  USE mp_interface,                    ONLY: mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE pslo,                            ONLY: pslo_com
  USE rgs_utils,                       ONLY: uinv
  USE ropt,                            ONLY: ropt_mod
  USE rotate_utils,                    ONLY: rottr,&
                                             rottr_c0_fnl,&
                                             rottr_fnl
  USE sfac,                            ONLY: fnla,&
                                             fnl,&
                                             fnl_packed,&
                                             il_fnl_packed
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             ncpw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zclean,&
                                             symmat_pack,&
                                             symmat_unpack

#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rgsvan

CONTAINS

  ! ==================================================================
  SUBROUTINE rgsvan(c0,nstate,smat,store_nonort)
    ! ==--------------------------------------------------------------==
    ! ==  Gram-Schmidt orthogonalization for Vanderbilt pp            ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(INOUT) __CONTIGUOUS &
                                             :: c0(:,:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: smat(nstate,nstate)
    INTEGER                                  :: ierr, isub
    INTEGER, ALLOCATABLE                     :: na_grp(:,:,:)
    INTEGER(int_8)                           :: il_smatpacked(1)
    CHARACTER(*), PARAMETER                  :: procedureN = 'rgsvan'
    LOGICAL, INTENT(IN)                      :: store_nonort
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: smatpacked(:)
#else
    REAL(real_8), ALLOCATABLE                :: smatpacked(:)
#endif
! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    CALL csmat(smat,c0,nstate,1,full=.FALSE.,store_nonort=store_nonort,only_parent=.TRUE.)
    IF(cntl%tlsd) THEN
       il_smatpacked=spin_mod%nsup*(spin_mod%nsup+1)/2+&
            spin_mod%nsdown*(spin_mod%nsdown+1)/2
    ELSE
       il_smatpacked=nstate*(nstate+1)/2
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_smatpacked,smatpacked,procedureN//'_smatpacked',ierr)
#else
    ALLOCATE(smatpacked(il_smatpacked(1)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF(paral%io_parent)THEN
       IF (cntl%tlsd) THEN
          CALL uinv('U',smat(1,1),nstate,spin_mod%nsup)
          CALL uinv('U',smat(spin_mod%nsup+1,spin_mod%nsup+1),nstate,spin_mod%nsdown)
          CALL symmat_pack(smat,smatpacked,nstate,spin_mod%nsup,spin_mod%nsdown)
       ELSE
          CALL uinv('U',smat,nstate,nstate)
          CALL symmat_pack(smat,smatpacked,nstate,nstate,0)
       ENDIF
    END IF
    CALL mp_bcast(smatpacked,INT(il_smatpacked(1),KIND=int_4),parai%io_source,parai%cp_grp)
    IF(cntl%tlsd)THEN
       CALL symmat_unpack(smat,smatpacked,nstate,spin_mod%nsup,spin_mod%nsdown,.FALSE.)
    ELSE
       CALL symmat_unpack(smat,smatpacked,nstate,nstate,0,.FALSE.)
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_smatpacked,smatpacked,procedureN//'_smatpacked',ierr)
#else
    DEALLOCATE(smatpacked,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    IF(pslo_com%tivan)THEN
       IF(cntl%distribute_fnl_rot)THEN
          CALL rottr_c0_fnl(ncpw%ngw,c0,INT(il_fnl_packed(1),KIND=int_4),fnl_packed,smat,nstate,redist=.TRUE.)
       ELSE
          CALL rottr(1._real_8,c0,smat,"N",nstate,ncpw%ngw,cntl%tlsd,spin_mod%nsup,&
               spin_mod%nsdown,use_cp=.TRUE.,redist=.TRUE.)
          CALL rottr_fnl(INT(il_fnl_packed(1),KIND=int_4),fnl_packed,nstate,smat)
       END IF
       IF(pslo_com%mixed_psp.OR.ropt_mod%calste) THEN
          ALLOCATE(na_grp(2,ions1%nsp,0:parai%cp_nogrp-1), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate na_grp',&
               __LINE__,__FILE__)
          CALL cp_grp_split_atoms(na_grp)
          CALL unpack_fnl(na_grp(:,:,parai%cp_inter_me),fnl_packed,unpacked=fnla)
          IF(parai%cp_nogrp.GT.1)CALL cp_grp_redist_dfnl_fnl(.TRUE.,.FALSE.,nstate,1)
          DEALLOCATE(na_grp, stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate na_grp',&
               __LINE__,__FILE__)
       END IF
    ELSE
       IF (ncpw%ngw.GT.0)&
            CALL rottr(1._real_8,c0,smat,"N",nstate,ncpw%ngw,cntl%tlsd,spin_mod%nsup,&
            spin_mod%nsdown,use_cp=.TRUE.,redist=.TRUE.)
    END IF
    IF (geq0) CALL zclean(c0,nstate,ncpw%ngw)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE rgsvan
  ! ==================================================================

END MODULE rgsvan_utils
