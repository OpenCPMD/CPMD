#include "cpmd_global.h"

MODULE fnlalloc_utils
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: imagp,&
                                             ndfnl,&
                                             nlps_com
  USE parac,                           ONLY: parai
  USE sfac,                            ONLY: dfnl,&
                                             dfnla,&
                                             fnl,&
                                             fnla,&
                                             fnl2,&
                                             ldf1,&
                                             ldf2,&
                                             tfnl2,&
                                             fnl_packed,&
                                             dfnl_packed,&
                                             il_dfnl_packed,&
                                             il_fnl_packed
  USE system,                          ONLY: cntl,&
                                             maxsys,&
                                             natpe,&
                                             nkpt,&
                                             norbpe,&
                                             parap
  USE reshaper,                        ONLY: reshape_inplace
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fnlalloc
  PUBLIC :: fnldealloc
  PUBLIC :: fnl_set

CONTAINS

  ! ==================================================================
  SUBROUTINE fnlalloc(nstate,tfor,tstress)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate,is
    LOGICAL                                  :: tfor, tstress

    CHARACTER(*), PARAMETER                  :: procedureN = 'fnlalloc'

    INTEGER                                  :: ierr, ldfnl, lfnl, &
                                                na_grp(2,ions1%nsp,0:parai%cp_nogrp-1),&
                                                na(2,ions1%nsp)

    ! split states between cp groups
    CALL cp_grp_split_atoms(na_grp)
    na(:,:)=na_grp(:,:,parai%cp_inter_me)
    il_fnl_packed(1)=0
    DO is=1,ions1%nsp
       il_fnl_packed(1)=il_fnl_packed(1)+(na(2,is)-na(1,is)+1)*nlps_com%ngh(is)
    END DO
    il_fnl_packed(1)=il_fnl_packed(1)*imagp
    il_fnl_packed(2)=nstate
    il_dfnl_packed(1)=il_fnl_packed(1)*3
    il_dfnl_packed(2)=nstate
    IF (cntl%tfdist) THEN
       tfnl2=.TRUE.

       lfnl=MAX(1,imagp*nstate*natpe*maxsys%nhxs*nkpt%nkpnt)
       IF (lfnl == 1) THEN
          ALLOCATE(fnl(1,1,1,1,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(fnl(imagp,natpe,maxsys%nhxs,nstate,nkpt%nkpnt),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF

       !CALL zeroing(fnl)!,lfnl)

       lfnl=MAX(1,imagp*norbpe*ions1%nat*maxsys%nhxs*nkpt%nkpnt)
       IF (lfnl == 1) THEN
          ALLOCATE(fnl2(1,1,1,1,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(fnl2(imagp,ions1%nat,maxsys%nhxs,norbpe,nkpt%nkpnt),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF

       !CALL zeroing(fnl2)!,lfnl)
       ldf1=natpe
       ldf2=norbpe

    ELSE

       tfnl2=.FALSE.
       lfnl=MAX(1,imagp*nstate*ions1%nat*maxsys%nhxs*nkpt%nkpnt)
       IF (lfnl == 1) THEN
          ALLOCATE(fnl(1,1,1,1,1),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ELSE
          ALLOCATE(fnl(imagp,ions1%nat,maxsys%nhxs,nstate,nkpt%nkpnt),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF

       !CALL zeroing(fnl)!,lfnl)

       fnl2 => fnl

       ldf1=ions1%nat
       ldf2=nstate

    ENDIF
    ! DFNL
    IF ( cntl%tshop ) THEN
       ndfnl=INT(nstate/parai%nproc)+1
    ELSE
       ndfnl=norbpe
    ENDIF
    ldfnl=imagp*3*ions1%nat*maxsys%nhxs*ndfnl*nkpt%nkpnt
    IF (ldfnl.LE.0) ldfnl=1
    IF (ldfnl == 1) THEN
       ALLOCATE(dfnl(1,1,1,1,1,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       ALLOCATE(dfnl(imagp,ions1%nat,maxsys%nhxs,3,ndfnl,nkpt%nkpnt),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    IF(imagp.EQ.1) &
         CALL reshape_inplace(fnl, (/SIZE(fnl,2),SIZE(fnl,3),SIZE(fnl,4)/),fnla)
    IF(imagp.EQ.1) &
         CALL reshape_inplace(dfnl, (/SIZE(dfnl,2),SIZE(dfnl,3),SIZE(dfnl,4),SIZE(dfnl,5)/),dfnla)
    !packed fnl/dfnl arrays
    ALLOCATE(fnl_packed(il_fnl_packed(1),il_fnl_packed(2)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnl_packed',&
         __LINE__,__FILE__)
    ALLOCATE(dfnl_packed(il_dfnl_packed(1),il_dfnl_packed(2)), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate dfnl_packed',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target enter data map(alloc:fnl,dfnl,fnl_packed,dfnl_packed)
#endif
    RETURN
  END SUBROUTINE fnlalloc
  ! ==================================================================
  SUBROUTINE fnldealloc(tfor,tstress)
    ! ==--------------------------------------------------------------==
    LOGICAL                                  :: tfor, tstress

    CHARACTER(*), PARAMETER                  :: procedureN = 'fnldealloc'

    INTEGER                                  :: ierr

! ==--------------------------------------------------------------==
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target exit data map(delete:fnl,dfnl,fnl_packed,dfnl_packed)
#endif
    if (imagp .eq. 1) nullify(fnla,dfnla)
    DEALLOCATE(fnl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    IF (tfnl2) DEALLOCATE(fnl2,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(dfnl,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(fnl_packed, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnl_packed',&
         __LINE__,__FILE__)
    DEALLOCATE(dfnl_packed, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate dfnl_packed',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fnldealloc
  ! ==================================================================
  SUBROUTINE fnl_set(tag)
    ! ==--------------------------------------------------------------==
    ! C.BEKAS

    CHARACTER(len=*)                         :: tag

    INTEGER                                  :: ldfd, ildum(2)
    INTEGER, SAVE                            :: ldf1b, ldf2b, ndfnlb, &
                                                ildfnl(2), ilfnl(2)
    LOGICAL, SAVE                            :: tfnl2b, tfnl2d
    REAL(real_8), POINTER  __CONTIGUOUS      :: xdum(:,:,:,:,:), &
                                                xdum2(:,:,:,:,:,:),&
                                                xdum3(:,:,:),&
                                                xdum4(:,:,:,:),&
                                                xdumfpack(:,:),&
                                                xdumdfpack(:,:)
    REAL(real_8), POINTER, SAVE __CONTIGUOUS :: xdfnl(:,:,:,:,:,:), &
                                                xfnl(:,:,:,:,:), &
                                                xfnl2(:,:,:,:,:),&
                                                xfnlp(:,:),&
                                                xdfnlp(:,:),&
                                                xdfnla(:,:,:,:),&
                                                xfnla(:,:,:)

    IF (INDEX(tag,'SAVE').NE.0) THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target exit data map(delete:fnl,dfnl,fnl_packed,dfnl_packed)
#endif
       xfnl   => fnl
       xfnl2  => fnl2
       xdfnl  => dfnl
       xfnla  => fnla
       xdfnla => dfnla
       xdfnlp => dfnl_packed
       xfnlp => fnl_packed

       ldf1b=ldf1
       ldf2b=ldf2
       ndfnlb=ndfnl
       tfnl2b=tfnl2
       ildfnl=il_dfnl_packed
       ilfnl=il_fnl_packed
    ELSEIF (INDEX(tag,'RECV').NE.0) THEN
       fnl   => xfnl
       fnl2  => xfnl2
       dfnl  => xdfnl
       fnla  => xfnla
       dfnla => xdfnla
       dfnl_packed => xdfnlp
       fnl_packed => xfnlp

       ldf1=ldf1b
       ldf2=ldf2b
       ndfnl=ndfnlb
       tfnl2=tfnl2b
       il_dfnl_packed=ildfnl
       il_fnl_packed=ilfnl
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target enter data map(alloc:fnl,dfnl,fnl_packed,dfnl_packed)
#endif
    ELSEIF (INDEX(tag,'SWITCH').NE.0) THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target exit data map(delete:fnl,dfnl,fnl_packed,dfnl_packed)
#endif
       xdum => fnl
       fnl => xfnl
       xfnl => xdum

       xdum => fnl2
       fnl2 => xfnl2
       xfnl2 => xdum

       xdum2 => dfnl
       dfnl => xdfnl
       xdfnl => xdum2

       xdum3  => fnla
       fnla   => xfnla
       xfnla  => xdum3

       xdum4  => dfnla
       dfnla  => xdfnla
       xdfnla => xdum4

       xdumdfpack => dfnl_packed
       dfnl_packed => xdfnlp
       xdfnlp => xdumdfpack

       xdumfpack => fnl_packed
       fnl_packed => xfnlp
       xfnlp => xdumfpack

       ldfd=ldf1
       ldf1=ldf1b
       ldf1b=ldfd
       ldfd=ldf2
       ldf2=ldf2b
       ldf2b=ldfd
       ldfd=ndfnl
       ndfnl=ndfnlb
       ndfnlb=ldfd
       tfnl2d=tfnl2
       tfnl2=tfnl2b
       tfnl2b=tfnl2d
       ildum=il_dfnl_packed
       il_dfnl_packed=ildfnl
       ildfnl=ildum
       ildum=il_fnl_packed
       il_fnl_packed=ilfnl
       ilfnl=ildum
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target enter data map(alloc:fnl,dfnl,fnl_packed,dfnl_packed)
#endif
    ELSEIF (INDEX(tag,'MIX').NE.0) THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target exit data map(delete:fnl,dfnl,fnl_packed,dfnl_packed)
#endif
       xdum2 => dfnl
       dfnl => xdfnl
       xdfnl => xdum2

       xdum4 => dfnla
       dfnla => xdfnla
       xdfnla => xdum4

       xdumdfpack => dfnl_packed
       dfnl_packed => xdfnlp
       xdfnlp => xdumdfpack

       ldfd=ndfnl
       ndfnl=ndfnlb
       ndfnlb=ldfd
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target enter data map(alloc:fnl,dfnl,fnl_packed,dfnl_packed)
#endif
    ELSE
       CALL stopgm('FNL_SET','INVALID TAG',&
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE fnl_set
  ! ==================================================================

END MODULE fnlalloc_utils
