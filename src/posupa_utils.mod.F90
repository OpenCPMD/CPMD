#include "cpmd_global.h"

MODULE posupa_utils
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes
  USE dotp_utils,                      ONLY: dotp_c1_cp,&
                                             dotp_c2_cp
  USE error_handling,                  ONLY: stopgm
  USE gpu
  USE geq0mod,                         ONLY: geq0
  USE harm,                            ONLY: dcgt,&
                                             dsgtw,&
                                             dtan2w,&
                                             wdsgt,&
                                             xmu
  USE jrotation_utils,                 ONLY: set_orbdist
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_sum
  USE nort,                            ONLY: nort_com
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE ropt,                            ONLY: ropt_mod
  USE rortog_utils,                    ONLY: give_scr_rortog,&
                                             rortog
  USE rotate_utils,                    ONLY: rotate
  USE reshaper,                        ONLY: reshape_inplace
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             ncpw,&
                                             parap,&
                                             paraw
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpar,                            ONLY: dt2hbe,&
                                             dt_elec,&
                                             dtb2me
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif
#ifdef __PARALLEL
  USE mpi_f08
#endif
!!use rotate_utils, only : rotate
!!se rotate_utils, only : rotate_da

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: posupa
  PUBLIC :: give_scr_posupa

CONTAINS

  ! ==================================================================
  SUBROUTINE posupa(c0,cm,c2,gamx,nstate,use_cp_grps)
    ! ==--------------------------------------------------------------==
    ! ==  UPDATE OF THE POSITIONS FOR VELOCITY VERLET                 ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(OUT)              :: c2(ncpw%ngw,*)
    REAL(real_8),INTENT(OUT)                 :: gamx(*)
    INTEGER,INTENT(IN)                       :: nstate
    COMPLEX(real_8),INTENT(INOUT)            :: cm(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    LOGICAL,INTENT(IN),OPTIONAL              :: use_cp_grps

#ifdef __PARALLEL
    INTEGER                                  :: i, ig, isub, nstx, ibeg_c0,&
                                                iend_c0, ngw_local, ierr
    type(MPI_COMM)                           :: gid
#else
    INTEGER                                  :: i, ig, isub, nstx, ibeg_c0,&
                                                iend_c0, ngw_local, ierr, &
                                                gid
#endif
    INTEGER(int_8)                           :: il_ai(2), il_ai1(2)
    LOGICAL                                  :: prtev, tnon, cp_active, geq0_local
    REAL(real_8)                             :: odt, pf1, pf2, pf3, &
                                                pf4, xi, xi_dt_elec
    REAL(real_8),POINTER __CONTIGUOUS        :: c2_r(:,:),cm_r(:,:),c0_r(:,:)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8),POINTER __CONTIGUOUS        :: ai(:,:), ai1(:,:)
#else
    REAL(real_8),ALLOCATABLE                 :: ai(:,:), ai1(:,:)
#endif
    CHARACTER(*), PARAMETER                  :: procedureN = 'posupa'
! Variables
! ==--------------------------------------------------------------==
! ==  WAVEFUNCTIONS                                               ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF(PRESENT(use_cp_grps))THEN
       cp_active=use_cp_grps
    ELSE
       cp_active=.FALSE.
    END IF
    IF(cp_active)THEN
       CALL cp_grp_get_sizes(ngw_l=ngw_local,geq0_l=geq0_local,&
         first_g=ibeg_c0,last_g=iend_c0)
       gid=parai%cp_grp
    ELSE
       ngw_local=ncpw%ngw
       geq0_local=geq0
       ibeg_c0=1
       iend_c0=ncpw%ngw
       gid=parai%allgrp
    END IF
    IF (tkpts%tkpnt) CALL stopgm('POSUPA','K POINTS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    IF (cntl%nonort) THEN
       IF(nort_com%scond.GT.nort_com%slimit)THEN
          il_ai(1)=3
          il_ai(2)=nstate
       ELSE
          il_ai(1)=nstate
          il_ai(2)=1
       END IF
       il_ai1=il_ai
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_ai,ai,procedureN//'_ai',ierr)
#else
       ALLOCATE(ai(il_ai(1),il_ai(2)),STAT=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ai',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_ai1,ai1,procedureN//'_ai1',ierr)
#else
       ALLOCATE(ai1(il_ai(1),il_ai1(2)),STAT=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate ai1',&
         __LINE__,__FILE__)
       ! ..LINEAR CONSTRAINTS
       IF (cntl%tharm) THEN
          !$omp parallel do private (i,ig,pf1,pf2,pf3,pf4)
          DO i=1,nstate
             DO ig=ibeg_c0,iend_c0
                pf1=dcgt(ig)
                pf2=dsgtw(ig)
                pf3=wdsgt(ig)
                c2(ig,i)=pf1*c0(ig,i)+pf2*cm(ig,i)
                cm(ig,i)=pf1*cm(ig,i)-pf3*c0(ig,i)
                pf4=dsgtw(ig)*dtan2w(ig)/dtb2me
                c0(ig,i)=pf4*c0(ig,i)
             ENDDO
             IF(nort_com%scond.GT.nort_com%slimit)THEN
                ai(:,i)=lagrange_afterrot(c0(ibeg_c0,i),c2(ibeg_c0,i),ngw_local,geq0_local)
             ELSE
                ai(i,1)=lagrange_norot(c2(ibeg_c0,i),ngw_local,geq0_local)
             END IF
          ENDDO
       ELSE IF(cntl%tmass)THEN
          !$omp parallel do private (i,ig,pf4)
          DO i=1,nstate
             DO ig=ibeg_c0,iend_c0
                pf4=cntr%emass/xmu(ig)
                c2(ig,i)=c0(ig,i)+dt_elec*cm(ig,i)
                c0(ig,i)=pf4*c0(ig,i)
             ENDDO
             IF(nort_com%scond.GT.nort_com%slimit)THEN
                ai(:,i)=lagrange_afterrot(c0(ibeg_c0,i),c2(ibeg_c0,i),ngw_local,geq0_local)
             ELSE
                ai(i,1)=lagrange_norot(c2(ibeg_c0,i),ngw_local,geq0_local)
             END IF
          ENDDO
       ELSE
          CALL reshape_inplace(c0,(/2*ncpw%ngw,nstate/),c0_r)
          CALL reshape_inplace(cm,(/2*ncpw%ngw,nstate/),cm_r)
          CALL reshape_inplace(c2,(/2*ncpw%ngw,nstate/),c2_r)
          CALL calc_ai(c0_r,cm_r,c2_r,ai,nstate,dt_elec,ibeg_c0*2-1,iend_c0*2,geq0_local,&
               nort_com%scond,nort_com%slimit)
       ENDIF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       comm_buffers_on_host=.FALSE.
#endif
       CALL mp_sum(ai,ai1,INT(il_ai(1)*il_ai(2)),gid)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       comm_buffers_on_host=.TRUE.
#endif
       CALL update_cm_c0(c0_r,cm_r,c2_r,ai1,nstate,dt_elec,ibeg_c0*2-1,iend_c0*2,geq0_local,&
            nort_com%scond,nort_com%slimit)
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_ai1,ai1,procedureN//'_ai1',ierr)
#else
       DEALLOCATE(ai1,STAT=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate ai1',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_ai,ai,procedureN//'_ai',ierr)
#else
       DEALLOCATE(ai,STAT=ierr)
#endif
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate ai',&
         __LINE__,__FILE__)
    ELSE
       IF (cntl%tharm) THEN
          !$omp parallel do private (I,IG,PF1,PF2,PF3,PF4)
          DO i=1,nstate
             DO ig=1,ncpw%ngw
                pf1=dcgt(ig)
                pf2=dsgtw(ig)
                pf3=wdsgt(ig)
                pf4=pf2*dtan2w(ig)/dt2hbe
                c2(ig,i)=pf1*c0(ig,i)+pf2*cm(ig,i)
                cm(ig,i)=pf1*cm(ig,i)-pf3*c0(ig,i)
                c0(ig,i)=pf4*c0(ig,i)
             ENDDO
          ENDDO
       ELSE
          !$omp parallel do private (I,IG)
          DO i=1,nstate
             DO ig=1,ncpw%ngw
                c2(ig,i)=c0(ig,i)+dt_elec*cm(ig,i)
             ENDDO
          ENDDO
          IF (cntl%tmass) THEN
             !$omp parallel do private (I,IG,PF4)
             DO i=1,nstate
                DO ig=1,ncpw%ngw
                   pf4=cntr%emass/xmu(ig)
                   c0(ig,i)=pf4*c0(ig,i)
                ENDDO
             ENDDO
          ENDIF
       ENDIF
       ! REORTHOGONALIZE THE WAVEFUNCTIONS
       prtev=ropt_mod%prteig .AND. .NOT.pslo_com%tivan
       tnon=cntl%tmass.OR.cntl%tharm
       CALL rortog(c0,c2,gamx,tnon,nstate,prtev)
       ! ADD CONSTRAINTS FOR WAVEFUNCTION VELOCITIES 
       IF (cntl%tharm) THEN
          !$omp parallel do private (I,IG,PF1,PF2,PF4)
          DO i=1,nstate
             DO ig=1,ncpw%ngw
                pf1=dcgt(ig)
                pf2=dsgtw(ig)
                pf4=cntr%delt_elec*(pf1/pf2)
                c0(ig,i)=pf4*c0(ig,i)
             ENDDO
          ENDDO
       ENDIF
       odt=1.0_real_8/dt_elec
       IF (ncpw%ngw.GT.0) THEN
          IF (cntl%tdmal) THEN
             CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
             CALL rotate_da(odt,c0,1._real_8,cm,gamx,2*ncpw%ngw,2*ncpw%ngw,nstate,&
                  paraw%nwa12(0,1),paraw%nwa12(0,2),nstx,parai%mepos,parap%pgroup,parai%nproc,&
                  parai%allgrp,cntl%tlsd,spin_mod%nsup,spin_mod%nsdown)
          ELSE
             CALL rotate(odt,c0,1.0_real_8,cm,gamx,nstate,2*ncpw%ngw,cntl%tlsd,&
                  spin_mod%nsup,spin_mod%nsdown)
          ENDIF
          CALL dcopy(2*nstate*ncpw%ngw,c2(1,1),1,c0(1,1),1)
       ENDIF
    ENDIF
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE posupa
  ! ==================================================================
  FUNCTION lagrange_afterrot(c0,c2,n,geq0_l) RESULT(ai)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: n
    LOGICAL,INTENT(IN)                       :: geq0_l
    COMPLEX(real_8),INTENT(IN)               :: c0(*),c2(*)
    REAL(real_8)                             :: ai(3)

    ai(1)=dotp_c1_cp(n,c0,geq0_l)
    ai(2)=dotp_c2_cp(n,c2,c0,geq0_l)
    ai(3)=dotp_c1_cp(n,c2,geq0_l)

  END FUNCTION lagrange_afterrot
  ! ==================================================================
  FUNCTION lagrange_norot(c2,n,geq0_l) RESULT(ai)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: n
    LOGICAL,INTENT(IN)                       :: geq0_l
    COMPLEX(real_8),INTENT(IN)               :: c2(*)
    REAL(real_8)                             :: ai

    ai=dotp_c1_cp(n,c2,geq0_l)

  END FUNCTION lagrange_norot
  ! ==================================================================
  SUBROUTINE calc_ai(c0_r,cm_r,c2_r,ai,nstate,dt_elec,ibeg_c0,iend_c0,geq0_local,scond,slimit)
    ! ==--------------------------------------------------------------==
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: c0_r(:,:), cm_r(:,:)
    REAL(real_8),INTENT(OUT) __CONTIGUOUS    :: c2_r(:,:), ai(:,:)
    REAL(real_8),INTENT(IN)                  :: dt_elec, scond, slimit
    INTEGER,INTENT(IN)                       :: nstate, ibeg_c0, iend_c0
    LOGICAL,INTENT(IN)                       :: geq0_local

    INTEGER                                  :: i,ig
    REAL(real_8)                             :: ai_1,ai_2,ai_3

#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute&
#else
    !$omp parallel do &
#endif
    !$omp& private (i,ig,ai_1,ai_2,ai_3)
    DO i=1,nstate
       IF(geq0_local)THEN
          c2_r(ibeg_c0,i)=c0_r(ibeg_c0,i)+dt_elec*cm_r(ibeg_c0,i)
          c2_r(ibeg_c0+1,i)=0.0_real_8
          IF(scond.GT.slimit)THEN
             ai_1=c0_r(ibeg_c0,i)**2*0.5_real_8
             ai_2=c0_r(ibeg_c0,i)*c2_r(ibeg_c0,i)*0.5_real_8
             ai_3=c2_r(ibeg_c0,i)**2*0.5_real_8
          ELSE
             ai_1=c2_r(ibeg_c0,i)**2*0.5_real_8
          END IF
       ELSE
          c2_r(ibeg_c0,i)=c0_r(ibeg_c0,i)+dt_elec*cm_r(ibeg_c0,i)
          c2_r(ibeg_c0+1,i)=c0_r(ibeg_c0+1,i)+dt_elec*cm_r(ibeg_c0+1,i)
          IF(scond.GT.slimit)THEN
             ai_1=c0_r(ibeg_c0,i)**2
             ai_2=c0_r(ibeg_c0,i)*c2_r(ibeg_c0,i)
             ai_3=c2_r(ibeg_c0,i)**2
             ai_1=ai_1+c0_r(ibeg_c0+1,i)**2
             ai_2=ai_2+c0_r(ibeg_c0+1,i)*c2_r(ibeg_c0+1,i)
             ai_3=ai_3+c2_r(ibeg_c0+1,i)**2
          ELSE
             ai_1=c2_r(ibeg_c0,i)**2
             ai_1=ai_1+c2_r(ibeg_c0+1,i)**2
          END IF
       END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp parallel do simd private(ig) &
#else
       !$omp simd &
#endif
       !$omp& reduction(+:ai_1,ai_2,ai_3)
       DO ig=ibeg_c0+2,iend_c0
          c2_r(ig,i)=c0_r(ig,i)+dt_elec*cm_r(ig,i)
          IF(scond.GT.slimit)THEN
             ai_1=ai_1+c0_r(ig,i)**2
             ai_2=ai_2+c0_r(ig,i)*c2_r(ig,i)
             ai_3=ai_3+c2_r(ig,i)**2
          ELSE
             ai_1=ai_1+c2_r(ig,i)**2
          END IF
       END DO
       IF(scond.GT.slimit)THEN
          ai(1,i)=ai_1*2.0_real_8
          ai(2,i)=ai_2*2.0_real_8
          ai(3,i)=ai_3*2.0_real_8
       ELSE
          ai(i,1)=ai_1*2.0_real_8
       END IF
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE calc_ai
  ! ==================================================================
  SUBROUTINE update_cm_c0(c0_r,cm_r,c2_r,ai,nstate,dt_elec,ibeg_c0,iend_c0,geq0_local,scond,slimit)
    ! ==--------------------------------------------------------------==
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: c0_r(:,:), cm_r(:,:), ai(:,:)
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: c2_r(:,:)
    REAL(real_8),INTENT(IN)                  :: dt_elec, scond, slimit
    INTEGER,INTENT(IN)                       :: nstate, ibeg_c0, iend_c0
    LOGICAL,INTENT(IN)                       :: geq0_local

    INTEGER                                  :: i,ig
    REAL(real_8)                             :: xi,xi_dt_elec
    IF(scond.GT.slimit)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute &
#else
       !$omp parallel do &
#endif
       !$omp& private(i,ig,xi,xi_dt_elec)
       DO i=1,nstate
          ai(3,i)=ai(3,i)-1.0_real_8
          xi=(-ai(2,i)+SQRT(ai(2,i)*ai(2,i)-ai(1,i)*ai(3,i)))/(ai(1,i))
          xi_dt_elec=xi/dt_elec
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp parallel do
#endif
          DO ig=ibeg_c0,iend_c0
             cm_r(ig,i)=cm_r(ig,i)+c0_r(ig,i)*xi_dt_elec
             c0_r(ig,i)=c2_r(ig,i)+c0_r(ig,i)*xi
          END DO
       END DO
    ELSE
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute &
#else
       !$omp parallel do &
#endif
       !$omp& private(i,ig,xi,xi_dt_elec)
       DO i=1,nstate
          xi=(-1+SQRT(2._real_8-ai(i,1)))
          xi_dt_elec=xi/dt_elec
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp parallel do
#endif
          DO ig=ibeg_c0,iend_c0
             cm_r(ig,i)=cm_r(ig,i)+c0_r(ig,i)*xi_dt_elec
             c0_r(ig,i)=c2_r(ig,i)+c0_r(ig,i)*xi
          END DO
       END DO
    ENDIF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE update_cm_c0
  ! ==================================================================
  SUBROUTINE give_scr_posupa(lposupa,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lposupa
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: lortog, lrotate, nstx

! ==--------------------------------------------------------------==

    IF (cntl%nonort) THEN
       lposupa=0
    ELSE
       IF (cntl%tdmal) THEN
          CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
          lrotate=nstate*nstx
       ELSE
          lrotate=0
       ENDIF
       CALL give_scr_rortog(lortog,tag,nstate)
       lposupa=MAX(lortog,lrotate)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_posupa
  ! ==================================================================

END MODULE posupa_utils
