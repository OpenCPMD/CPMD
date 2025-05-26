#include "cpmd_global.h"

MODULE cp_grp_utils
  USE error_handling,                  ONLY: stopgm
  USE geq0mod,                         ONLY: geq0
  use gpu
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE fnl_utils,                       ONLY: pack_fnl,&
                                             unpack_fnl,&
                                             pack_dfnl,&
                                             unpack_dfnl
  USE kinds,                           ONLY: int_1,&
                                             int_2,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp,&
                                             nlps_com
  USE parac,                           ONLY: parai
  USE part_1d,                         ONLY: part_1d_get_blk_bounds
  USE system,                          ONLY: ncpw,&
                                             nkpt,&
                                             maxsys,&
                                             parap,&
                                             norbpe
  USE sfac,                            ONLY: dfnl,&
                                             fnl,&
                                             dfnla,&
                                             fnla
  USE timer,                           ONLY: tihalt,&
                                             tiset
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cp_grp_get_sizes
  PUBLIC :: cp_grp_redist
  PUBLIC :: cp_grp_redist_array
  PUBLIC :: cp_grp_redist_array_f
  PUBLIC :: cp_grp_split_atoms
  PUBLIC :: cp_grp_redist_dfnl_fnl

  INTERFACE cp_grp_redist
     MODULE PROCEDURE cp_grp_redist_d
     MODULE PROCEDURE cp_grp_redist_z
  END INTERFACE cp_grp_redist
  !TK redistributing arrays distributed along the second dimension
  ! is much cheaper using allgatherv instead of allreduce
  INTERFACE cp_grp_redist_array
     MODULE PROCEDURE cp_grp_redist_array_r1_d
     MODULE PROCEDURE cp_grp_redist_array_r1_z
     MODULE PROCEDURE cp_grp_redist_array_r2_d
     MODULE PROCEDURE cp_grp_redist_array_r2_z
  END INTERFACE cp_grp_redist_array

  INTEGER,DIMENSION(:,:),ALLOCATABLE,PUBLIC :: cp_grp_get_cp_rank

CONTAINS
  ! ==================================================================
  SUBROUTINE cp_grp_get_sizes(ngwk_l,ngw_l,nhg_l,first_nhg,last_nhg,&
       geq0_l,igeq0_l,firstk_g,lastk_g,first_g,last_g,i_g0_l)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(out), OPTIONAL           :: ngwk_l, ngw_l, nhg_l
    LOGICAL, INTENT(out), OPTIONAL           :: geq0_l
    INTEGER, INTENT(out), OPTIONAL           :: igeq0_l, firstk_g, lastk_g, &
                                                first_g, last_g, i_g0_l,&
                                                first_nhg, last_nhg

    INTEGER, SAVE                            :: first, firstk, i_g0, last, &
                                                lastk, igeq0_l_s,&
                                                first_nh, last_nh

    LOGICAL, SAVE                            :: first_call=.TRUE.
    IF(first_call)THEN
       IF(parai%cp_nogrp.GT.1)THEN
          CALL part_1d_get_blk_bounds(nkpt%ngwk,parai%cp_inter_me,parai%cp_nogrp,firstk,&
               lastk)
       ELSE
          firstk=1
          lastk=nkpt%ngwk
       END IF
       !
       ! index of the G=0 element
       i_g0 = 0
       IF (tkpts%tkpnt) THEN
          IF (GEq0) i_g0 = 1 + ncpw%ngw
       ELSE
          IF (GEq0) i_g0 = 1
       ENDIF
       !
       !
       IF(parai%cp_nogrp.GT.1)THEN
          CALL part_1d_get_blk_bounds(ncpw%ngw,parai%cp_inter_me,parai%cp_nogrp,first,last)
       ELSE
          first=1
          last=ncpw%ngw
       END IF
       !
       !
       IF(parai%cp_nogrp.GT.1)THEN
          CALL part_1d_get_blk_bounds(ncpw%nhg,parai%cp_inter_me,parai%cp_nogrp,first_nh,&
               last_nh)
       ELSE
          first_nh=1
          last_nh=ncpw%nhg
       END IF
       !
       !
       igeq0_l_s=0
       IF (i_g0.GE.firstk.AND.i_g0.LE.lastk) igeq0_l_s = parai%cp_me
       CALL mp_sum(igeq0_l_s,parai%cp_grp)

       first_call=.FALSE.
    END IF

    IF (PRESENT(firstk_g)) firstk_g = firstk
    IF (PRESENT(lastk_g)) lastk_g = lastk
    IF (PRESENT(ngwk_l)) ngwk_l = lastk - firstk + 1
    IF (PRESENT(first_g)) first_g = first
    IF (PRESENT(last_g)) last_g = last
    IF (PRESENT(ngw_l)) ngw_l = last - first + 1
    IF (PRESENT(first_nhg)) first_nhg = first_nh
    IF (PRESENT(last_nhg)) last_nhg = last_nh
    IF (PRESENT(nhg_l)) nhg_l = last_nh - first_nh + 1
    !
    !
    IF (PRESENT(geq0_l)) geq0_l = i_g0.GE.firstk.AND.i_g0.LE.lastk

    IF (PRESENT(igeq0_l)) igeq0_l = igeq0_l_s
    IF (PRESENT(i_g0_l)) THEN
       i_g0_l = -HUGE(0)
       IF (i_g0.GE.firstk.AND.i_g0.LE.lastk) i_g0_l = i_g0 - firstk + 1
    ENDIF

       ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_get_sizes
  ! ==================================================================
  SUBROUTINE cp_grp_redist_z(DATA,ld,n)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: DATA(:,:)
    INTEGER                                  :: ld, n

    IF (parai%cp_nogrp.EQ.1) RETURN
    CALL mp_sum(DATA,ld*n,parai%cp_inter_grp)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_z
  ! ==================================================================
  SUBROUTINE cp_grp_redist_d(DATA,ld,n)
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: DATA(:,:)
    INTEGER                                  :: ld, n

    IF (parai%cp_nogrp.EQ.1) RETURN
    CALL mp_sum(DATA,ld*n,parai%cp_inter_grp)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_d

  ! ==--------------------------------------------------------------==
  SUBROUTINE cp_grp_redist_array_r1_z(data,n)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: n
    COMPLEX(real_8)                          :: data(*)
    INTEGER                                  :: n_proc(3,parai%cp_nogrp),&
                                                revcnt(parai%cp_nogrp),&
                                                displ(parai%cp_nogrp),i
    !  ==--------------------------------------------------------------==
    IF(parai%CP_NOGRP.EQ.1) RETURN
    DO i = 1,parai%cp_nogrp
       CALL part_1d_get_blk_bounds(n,i-1,parai%cp_nogrp,n_proc(1,i),n_proc(2,i))
       n_proc(3,i)=n_proc(2,i)- n_proc(1,i)+1
    END DO
    DO i = 1,parai%cp_nogrp
       revcnt(i)=n_proc(3,i)*2
    END DO
    displ=0
    do i= 2,parai%cp_nogrp
       displ(i)=displ(i-1)+n_proc(3,i-1)*2
    END DO
    CALL my_concatv_inplace(data,revcnt,displ,parai%cp_inter_grp)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_array_r1_z
  ! ==================================================================
  SUBROUTINE cp_grp_redist_array_r1_d(data,n)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: n
    REAL(real_8),INTENT(INOUT)               :: data(*)
    INTEGER                                  :: n_proc(3,parai%cp_nogrp),&
                                                revcnt(parai%cp_nogrp),&
                                                displ(parai%cp_nogrp),i
    ! ==--------------------------------------------------------------==
    IF(parai%CP_NOGRP.EQ.1) RETURN
    DO i = 1,parai%cp_nogrp
       CALL part_1d_get_blk_bounds(n,i-1,parai%cp_nogrp,n_proc(1,i),n_proc(2,i))
       n_proc(3,i)=n_proc(2,i)- n_proc(1,i)+1
    END DO
    DO i = 1,parai%cp_nogrp
       revcnt(i)=n_proc(3,i)
    END DO
    displ=0
    do i= 2,parai%cp_nogrp
       displ(i)=displ(i-1)+n_proc(3,i-1)
    END DO
    CALL my_concatv_inplace(data,revcnt,displ,parai%cp_inter_grp)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_array_r1_d
  ! ==================================================================
  SUBROUTINE cp_grp_redist_array_r2_z(data,ld,n)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: ld,n
    COMPLEX(real_8),INTENT(INOUT)            :: data(ld,*)
    INTEGER                                  :: n_proc(3,parai%cp_nogrp),&
                                                revcnt(parai%cp_nogrp),&
                                                displ(parai%cp_nogrp),i
    !  ==--------------------------------------------------------------==
    IF(parai%CP_NOGRP.EQ.1) RETURN
    DO i = 1,parai%cp_nogrp
       CALL part_1d_get_blk_bounds(n,i-1,parai%cp_nogrp,n_proc(1,i),n_proc(2,i))
       n_proc(3,i)=n_proc(2,i)- n_proc(1,i)+1
    END DO
    DO i = 1,parai%cp_nogrp
       revcnt(i)=n_proc(3,i)*ld*2
    END DO
    displ=0
    DO i= 2,parai%cp_nogrp
       displ(i)=displ(i-1)+n_proc(3,i-1)*ld*2
    END DO
    CALL my_concatv_inplace(data,revcnt,displ,parai%cp_inter_grp)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_array_r2_z
  ! ==================================================================
  SUBROUTINE cp_grp_redist_array_r2_d(data,ld,n)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: ld,n
    REAL(real_8),INTENT(INOUT)               :: data(ld,*)
    INTEGER                                  :: n_proc(3,parai%cp_nogrp),&
                                                revcnt(parai%cp_nogrp),&
                                                displ(parai%cp_nogrp),i
    ! ==--------------------------------------------------------------==
    IF(parai%CP_NOGRP.EQ.1) RETURN
    DO i = 1,parai%cp_nogrp
       CALL part_1d_get_blk_bounds(n,i-1,parai%cp_nogrp,n_proc(1,i),n_proc(2,i))
       n_proc(3,i)=n_proc(2,i)- n_proc(1,i)+1
    END DO
    DO i = 1,parai%cp_nogrp
       revcnt(i)=n_proc(3,i)*ld
    END DO
    displ=0
    DO i= 2,parai%cp_nogrp
       displ(i)=displ(i-1)+n_proc(3,i-1)*ld
    END DO
    CALL my_concatv_inplace(data,revcnt,displ,parai%cp_inter_grp)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_array_r2_d
  ! ==--------------------------------------------------------------==
  SUBROUTINE cp_grp_redist_array_f(data,ld,n)
    ! ==------------------------------------------------------------==
    ! redistributes an array distributed along the first dimention
    ! Author: Tobias Kloeffel, FAU Erlangen-Nuernberg April 2019
    ! ==------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: ld,n
    COMPLEX(real_8),INTENT(INOUT)            :: data(ld,*)
    INTEGER                                  :: ld_group(3,parai%cp_nogrp),&
                                                  revcnt,i,ig,ierr,group,my_group
    INTEGER(int_8)                           :: il_buffer(3)
#ifdef _USE_SCRATCHLIBRARY
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: buffer(:,:,:)
#else
    COMPLEX(real_8), ALLOCATABLE             :: buffer(:,:,:)
#endif
    CHARACTER(*), PARAMETER                  :: procedureN = 'cp_grp_redist_array_f'
    !  ==--------------------------------------------------------------==
    IF(parai%CP_NOGRP.EQ.1) RETURN
    DO i = 1,parai%cp_nogrp
       CALL part_1d_get_blk_bounds(ld,i-1,parai%cp_nogrp,ld_group(1,i),ld_group(2,i))
       ld_group(3,i)=ld_group(2,i)-ld_group(1,i)+1
    END DO
    il_buffer(1)=MAXVAL(ld_group(3,:))
    il_buffer(2)=n
    il_buffer(3)=parai%cp_nogrp
    revcnt=il_buffer(1)*n*2
    my_group=parai%cp_inter_me+1
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_buffer,buffer,procedureN//'_buffer',ierr)
#else
    ALLOCATE(buffer(il_buffer(1),il_buffer(2),il_buffer(3)),stat=ierr)
#endif
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
         __LINE__,__FILE__)
    group=parai%cp_inter_me+1
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute &
#else
    !$omp parallel do &
#endif
    !$omp& private(i,ig) 
    DO i=1,n
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do private(ig)
#endif
       DO ig=ld_group(1,group),ld_group(2,group)
          buffer(ig-ld_group(1,group)+1,i,group)=data(ig,i)
       END DO
    END DO
    CALL my_concat_inplace(buffer,revcnt,parai%cp_inter_grp)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams &
#else
    !$omp parallel &
#endif
    !$omp& private(i,ig,group) 
    DO group=1,parai%cp_nogrp
       IF(group.EQ.my_group) CYCLE
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp distribute
#else
       !$omp do
#endif
       DO i=1,n
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp parallel do private(ig)
#else
          !$omp simd 
#endif
          DO ig=ld_group(1,group),ld_group(2,group)
             data(ig,i)=buffer(ig-ld_group(1,group)+1,i,group)
          END DO
       END DO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp end do nowait
#endif
    END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end target teams
#else
    !$omp end parallel
#endif
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_buffer,buffer,procedureN//'_buffer',ierr)
#else
    DEALLOCATE(buffer,stat=ierr)
#endif
    IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    RETURN
    ! ==--------------------------------------------------------------==
  END SUBROUTINE cp_grp_redist_array_f
  ! ==--------------------------------------------------------------==
  SUBROUTINE cp_grp_split_atoms(na_split)
    INTEGER,INTENT(OUT)                    :: na_split(2,ions1%nsp,parai%cp_nogrp)
    INTEGER                                :: is,ia,iv,num(maxsys%nhxs,parai%cp_nogrp+1),&
                                              extra(parai%cp_nogrp), grp, aim, &
                                              distributed(parai%cp_nogrp), max_distributed, &
                                              counter, i, residual, isub, ierr
    LOGICAL,SAVE                           :: first=.true.
    INTEGER,ALLOCATABLE,SAVE               :: na_work(:,:,:)
    CHARACTER(*), PARAMETER                :: procedureN = 'cp_grp_split_atoms'
    ! ==------------------------------------------------------------==
    ! generates custom na mapping, splits atoms between groups
    ! avoids loadimablance by taking betaprojectors into account
    ! Author: Tobias Kloeffel, FAU Erlangen-Nuernberg April 2019
    ! ==------------------------------------------------------------==
    CALL tiset(procedureN,isub)

    IF (first) THEN
       ALLOCATE(na_work(2,ions1%nsp,parai%cp_nogrp), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate na_work',&
            __LINE__,__FILE__)

       !determine the number of atoms with the same nlps_com%ngh(is)
       num=0
       DO iv=1,maxsys%nhxs
          DO is=1,ions1%nsp
             IF(nlps_com%ngh(is).EQ.iv) num(iv,1)=num(iv,1)+ions0%na(is)
          END DO
       END DO

       !divide the sum of atoms with the same nlps_com%ngh(is) into chunks for each group
       DO iv=1,maxsys%nhxs
          DO grp=1,parai%cp_nogrp
             num(iv,grp+1)=FLOOR(DBLE(num(iv,1))/DBLE(parai%cp_nogrp))
          END DO
       END DO

       !check for missing atoms
       distributed=0
       extra=0
       max_distributed=0
       DO iv=1,maxsys%nhxs
          counter=0
          DO grp=1,parai%cp_nogrp
             counter=counter+num(iv,grp+1)
          END DO
          IF (counter.LT.num(iv,1) ) THEN
             residual=num(iv,1)-counter
             !distribute es evenly as possible
             DO i=1,2
                max_distributed=minval(extra,1)
                DO grp=1,parai%cp_nogrp
                   !This group already got an extra atom
                   IF (extra(grp).GT.max_distributed) CYCLE
                   !Some atoms need to be distributed
                   IF (residual.GT.0) THEN
                      extra(grp)=extra(grp)+1
                      residual=residual-1
                      num(iv,grp+1)=num(iv,grp+1)+1
                   END IF
                END DO
             END DO
          END IF
       END DO

       !sanity check
       DO iv=1,maxsys%nhxs
          counter=0
          DO grp=1,parai%cp_nogrp
             counter=counter+num(iv,grp+1)
          END DO
          IF (counter.NE.num(iv,1) ) CALL stopgm(procedureN, 'something went wront',&
            __LINE__,__FILE__)
       END DO

       na_work(1,:,:)=0
       na_work(2,:,:)=-1

       !build up na_work
       DO iv=1,maxsys%nhxs
          DO grp=1,parai%cp_nogrp
             counter=0
             DO is=1,ions1%nsp
                IF (nlps_com%ngh(is).EQ.iv) THEN
                   aim=num(iv,grp+1)
                   max_distributed=0
                   IF (grp.GT.1) THEN
                      DO i=1,grp -1
                         max_distributed=max_distributed+num(iv,i+1)
                      END DO
                   END IF
                   DO ia=1,ions0%na(is)
                      IF(max_distributed.LE.counter.AND.aim.GE.(counter-max_distributed+1))&
                           THEN
                         IF(na_work(1,is,grp).EQ.0) na_work(1,is,grp)=ia
                         na_work(2,is,grp)=ia
                         counter=counter+1
                      ELSE
                         counter=counter+1
                      END IF
                   END DO
                END IF
             END DO
          END DO
       END DO
       first=.FALSE.
    END IF

    na_split=na_work

    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE cp_grp_split_atoms
  ! ==--------------------------------------------------------------==
  SUBROUTINE cp_grp_redist_dfnl_fnl(redist_fnl,redist_dfnl,nstate,ikind)
    LOGICAL,INTENT(IN)                     :: redist_fnl, redist_dfnl
    INTEGER,INTENT(IN)                     :: nstate, ikind
    CHARACTER(*), PARAMETER                :: procedureN = 'cp_grp_redist_dfnl_fnl'
    INTEGER                                :: worksum(parai%cp_nogrp), is,ierr, ii,&
                                              sendcnt, grp, isub
    INTEGER, ALLOCATABLE                   :: na_grp(:,:,:)
    INTEGER(int_8)                         :: il_temp(3)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS     :: temp(:,:,:)
#else
    REAL(real_8), ALLOCATABLE              :: temp(:,:,:)
#endif
    ! ==------------------------------------------------------------==
    ! redistributes dfnl and/or fnl arrays
    ! Author: Tobias Kloeffel, FAU Erlangen-Nuernberg April 2019
    ! ==------------------------------------------------------------==
    IF (parai%cp_nogrp.EQ.1) RETURN
    CALL tiset(procedureN,isub)
    !get cp_grp atomwise (d)fnl distribution
    ALLOCATE(na_grp(2,ions1%nsp,parai%cp_nogrp), stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate na_grp',&
         __LINE__,__FILE__)
    CALL cp_grp_split_atoms(na_grp)

    worksum=0
    !count data sent to other groups
    DO grp=1,parai%cp_nogrp
       DO is=1,ions1%nsp
          worksum(grp)=worksum(grp)+(na_grp(2,is,grp)-na_grp(1,is,grp)+1)*nlps_com%ngh(is)
       END DO
    END DO
    IF(redist_fnl)THEN
       il_temp(1)=imagp*MAXVAL(worksum)
       il_temp(2)=nstate
       il_temp(3)=parai%cp_nogrp
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_temp,temp,procedureN//'_temp',ierr)
#else
       ALLOCATE(temp(il_temp(1),il_temp(2),il_temp(3)), stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate temp',&
            __LINE__,__FILE__)
       !pack local fnl chunk
       grp=parai%cp_inter_me+1
       IF(tkpts%tkpnt)THEN
          CALL pack_fnl(na_grp(:,:,grp),temp(:,:,grp),unpacked_k=fnl(:,:,:,:,ikind))
       ELSE
          CALL pack_fnl(na_grp(:,:,grp),temp(:,:,grp),unpacked=fnla)
       END IF
       sendcnt=il_temp(1)*il_temp(2)
       !exchange data
       CALL my_concat_inplace(temp,sendcnt,parai%cp_inter_grp)
       !unpack back to local fnl
       DO grp=1,parai%cp_nogrp
          IF (grp.EQ.parai%cp_inter_me+1) CYCLE
          IF(tkpts%tkpnt)THEN
             CALL unpack_fnl(na_grp(:,:,grp),temp(:,:,grp),unpacked_k=fnl(:,:,:,:,ikind))
          ELSE
             CALL unpack_fnl(na_grp(:,:,grp),temp(:,:,grp),unpacked=fnla)
          END IF
       END DO
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_temp,temp,procedureN//'_temp',ierr)
#else
       DEALLOCATE(temp, stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate temp',&
            __LINE__,__FILE__)
    END IF

    IF(redist_dfnl)THEN
       il_temp(1)=imagp*MAXVAL(worksum)*3
       il_temp(2)=norbpe
       il_temp(3)=parai%cp_nogrp
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_temp,temp,procedureN//'_temp',ierr)
#else
       ALLOCATE(temp(il_temp(1),il_temp(2),il_temp(3)), stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate temp',&
            __LINE__,__FILE__)
       !pack local dfnl chunk
       grp=parai%cp_inter_me+1
       IF(tkpts%tkpnt)THEN
          CALL pack_dfnl(na_grp(:,:,grp),temp(:,:,grp),unpacked_k=dfnl(:,:,:,:,:,ikind))
       ELSE
          CALL pack_dfnl(na_grp(:,:,grp),temp(:,:,grp),unpacked=dfnla)
       END IF
       sendcnt=il_temp(1)*il_temp(2)
       !exchange data
       CALL my_concat_inplace(temp,sendcnt,parai%cp_inter_grp)
       !unpack to local dfnl
       DO grp=1,parai%cp_nogrp
          IF (grp.EQ.parai%cp_inter_me+1) CYCLE
          IF(tkpts%tkpnt)THEN
             CALL unpack_dfnl(na_grp(:,:,grp),temp(:,:,grp),unpacked_k=dfnl(:,:,:,:,:,ikind))
          ELSE
             CALL unpack_dfnl(na_grp(:,:,grp),temp(:,:,grp),unpacked=dfnla)
          END IF
       END DO
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_temp,temp,procedureN//'_temp',ierr)
#else
       DEALLOCATE(temp, stat=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate temp',&
            __LINE__,__FILE__)
    END IF

    DEALLOCATE(na_grp, stat=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate na_grp',&
         __LINE__,__FILE__)

    CALL tihalt(procedureN,isub)
    RETURN
  END SUBROUTINE cp_grp_redist_dfnl_fnl
  ! ==--------------------------------------------------------------==
END MODULE cp_grp_utils
! ==================================================================

SUBROUTINE cp_grp_zero_g(c0,n,m,first_g,last_g)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: n, m
  COMPLEX(real_8)                            :: c0(n,m)
  INTEGER                                    :: first_g, last_g

  COMPLEX(real_8), PARAMETER :: zzero = (0.0_real_8,0.0_real_8)

  INTEGER                                    :: i, j

  IF (first_g.EQ.1.AND.last_g.EQ.n) RETURN
  !$omp parallel do &
  !$omp            private(i) &
  !$omp            shared(first_g,last_g,N)
  DO j = 1,m
     DO i = 1,first_g-1
        c0(i,j) = zzero
     ENDDO
     DO i = last_g+1,n
        c0(i,j) = zzero
     ENDDO
  ENDDO
  !$omp end parallel do
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_zero_g
! ==================================================================

SUBROUTINE cp_grp_zero_states(c0,n,m,ibeg,iend)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE parac , ONLY:parai
  IMPLICIT NONE
  INTEGER                                    :: n, m
  COMPLEX(real_8)                            :: c0(n,m)
  INTEGER                                    :: ibeg(0:parai%cp_nproc-1), &
                                                iend(0:parai%cp_nproc-1)

  COMPLEX(real_8), PARAMETER :: zzero = (0.0_real_8,0.0_real_8)

  INTEGER                                    :: nn

  IF (ibeg(parai%cp_me).GT.1) THEN
     nn = n * ( ibeg(parai%cp_me) - 1 )
     CALL zcopy(nn,zzero,0,c0(1,1),1)
  ENDIF
  IF (iend(parai%cp_me).LT.m) THEN
     nn = n * ( m - iend(parai%cp_me) - 1 )
     CALL zcopy(nn,zzero,0,c0(1,iend(parai%cp_me)+1),1)
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_zero_states
! ==================================================================
SUBROUTINE cp_grp_copy_wfn_to_local(c0,ld,C0_l,LD_l,ibeg,M_l,n)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: ld
  COMPLEX(real_8)                            :: c0(ld,*)
  INTEGER                                    :: LD_l
  COMPLEX(real_8)                            :: C0_l(LD_l,*)
  INTEGER                                    :: ibeg, M_l, n

  INTEGER                                    :: i

  IF (LD_l.GT.0.AND.ld.GT.0) THEN
     ! copy the right C0 block into the local array
     DO i=1,n
        CALL dcopy(2*M_l,c0(ibeg,i),1,C0_l(1,i),1)
     ENDDO
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_copy_wfn_to_local
! ==================================================================
SUBROUTINE cp_grp_copy_local_to_wfn(C0_l,LD_l,c0,ld,ibeg,M_l,n)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: LD_l
  COMPLEX(real_8)                            :: C0_l(LD_l,*)
  INTEGER                                    :: ld
  COMPLEX(real_8)                            :: c0(ld,*)
  INTEGER                                    :: ibeg, M_l, n

  INTEGER                                    :: i

  IF (LD_l.GT.0.AND.ld.GT.0) THEN
     ! copy back to the right C0 block
     DO i=1,n
        CALL dcopy(2*M_l,C0_l(1,i),1,c0(ibeg,i),1)
     ENDDO
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_copy_local_to_wfn
! ==================================================================
SUBROUTINE cp_grp_add_local_to_wfn(C0_l,LD_l,c0,ld,ibeg,M_l,n)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  IMPLICIT NONE
  INTEGER                                    :: LD_l
  COMPLEX(real_8)                            :: C0_l(LD_l,*)
  INTEGER                                    :: ld
  COMPLEX(real_8)                            :: c0(ld,*)
  INTEGER                                    :: ibeg, M_l, n

  INTEGER                                    :: i, j

  IF (LD_l.GT.0.AND.ld.GT.0) THEN
     ! copy back to the right C0 block
     !$omp parallel do private(i,j) shared(N,M_l,ibeg)
     DO i=1,n
        DO j=1,M_l
           c0(ibeg+j-1,i) = c0(ibeg+j-1,i) + C0_l(j,i)
        ENDDO
     ENDDO
     !$omp end parallel do
  ENDIF
  ! ==--------------------------------------------------------------==
  RETURN
  ! ==--------------------------------------------------------------==
END SUBROUTINE cp_grp_add_local_to_wfn
