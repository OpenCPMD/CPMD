#include "cpmd_global.h"

MODULE csize_utils
  USE cp_grp_utils,                    ONLY: cp_grp_get_sizes
  USE dotp_utils,                      ONLY: dotp_c1_cp
  USE elct,                            ONLY: crge
  USE geq0mod,                         ONLY: geq0
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_max,&
                                             mp_sum
  USE parac,                           ONLY: parai
  USE reshaper,                        ONLY: reshape_inplace
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             nkpt,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
#ifdef __PARALLEL
  USE mpi_f08
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: csize

CONTAINS

  ! ==================================================================
  SUBROUTINE csize(c2,nstate,gemax,cnorm,use_cp_grps,special)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate
    COMPLEX(real_8),INTENT(IN)               :: c2(nkpt%ngwk,nstate)
    REAL(real_8),INTENT(OUT)                 :: gemax, cnorm
    LOGICAL,INTENT(IN),OPTIONAL              :: use_cp_grps, special

    CHARACTER(*), PARAMETER                  :: procedureN = 'csize'

#ifdef __PARALLEL
    INTEGER                                  :: i, iiabs, nocc, isub, &
                                                ngwk_local, ibeg_c0
    type(MPI_COMM)                           :: gid
#else
    INTEGER                                  :: i, iiabs, nocc, isub, &
                                                ngwk_local, ibeg_c0, gid
#endif
    REAL(real_8),POINTER __CONTIGUOUS        :: c2_r(:,:)
    LOGICAL                                  :: geq0_local, cp_active, sp
    REAL(real_8), EXTERNAL                   :: ddot

    INTEGER, EXTERNAL                        :: izamax
    ! ==--------------------------------------------------------------==
    !TK in noforce cnorm and gemax were multiplied by the occupation number
    ! special case to get it done in the same way...
    CALL tiset(procedureN,isub)
    IF(PRESENT(special))THEN
       sp=special
    ELSE
       sp=.FALSE.
    END IF
    IF(PRESENT(use_cp_grps))THEN
       cp_active=use_cp_grps
    ELSE
       cp_active=.FALSE.
    END IF
    IF(cp_active)THEN
       CALL cp_grp_get_sizes(ngwk_l=ngwk_local,firstk_g=ibeg_c0,geq0_l=geq0_local)
       gid=parai%cp_grp
    ELSE
       ngwk_local=nkpt%ngwk
       ibeg_c0=1
       geq0_local=geq0
       gid=parai%allgrp
    END IF
    CALL reshape_inplace(c2,(/2*nkpt%ngwk,nstate/),c2_r)
    gemax=0.0_real_8
    cnorm=0.0_real_8
    nocc=0
    IF(ngwk_local.GT.0) &
         CALL csize_r(c2_r,ibeg_c0,nkpt%ngwk,ngwk_local,nstate,geq0_local,sp,gemax,cnorm,nocc)
    CALL mp_sum(cnorm,gid)
    CALL mp_max(gemax,gid)
    cnorm=SQRT(cnorm/REAL(nocc*spar%ngwks,kind=real_8))
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE csize
  ! ==================================================================
  SUBROUTINE csize_r(c2_r,ibeg,ngw,ngw_local,nstate,geq0,sp,gemax,cnorm,nocc)
    INTEGER,INTENT(IN)                       :: nstate, ibeg, ngw, ngw_local
    REAL(real_8),INTENT(IN)                  :: c2_r(2,ngw,*)
    REAL(real_8),INTENT(INOUT)               :: gemax, cnorm
    INTEGER,INTENT(INOUT)                    :: nocc
    LOGICAL,INTENT(IN)                       :: geq0, sp

    INTEGER                                  :: i,ig
    REAL(real_8)                             :: cnorm_l, gemax_l, gemax_ll
    IF (tkpts%tkpnt) THEN
       !$omp parallel do private(i,ig,cnorm_l,gemax_l,gemax_ll) &
       !$omp& reduction(max:gemax) reduction(+:cnorm,nocc)
       DO i=1,nstate
          IF(crge%f(i,1).LE.1.e-5_real_8) CYCLE
          nocc=nocc+1
          cnorm_l=0._real_8
          gemax_l=0._real_8
          gemax_ll=0._real_8
          DO ig=ibeg,ngw_local
             cnorm_l=cnorm_l+c2_r(1,ig,i)**2+c2_r(2,ig,i)**2
             IF(gemax_ll.LT.ABS(c2_r(1,ig,i))+ABS(c2_r(2,ig,i)))THEN
                gemax_ll=ABS(c2_r(1,ig,i))+ABS(c2_r(2,ig,i))
                gemax_l=ABS(c2_r(1,ig,i))**2+ABS(c2_r(2,ig,i))**2
             END IF
          END DO
          IF(sp)THEN
             cnorm=cnorm+cnorm_l*crge%f(i,1)
             gemax=MAX(ABS(crge%f(i,1)*gemax_l),gemax)
          ELSE
             cnorm=cnorm+cnorm_l
             gemax=MAX(gemax_l,gemax)
          END IF
       END DO
    ELSE
       !$omp parallel do private(i,ig,cnorm_l,gemax_l,gemax_ll) &
       !$omp& reduction(max:gemax) reduction(+:cnorm,nocc)
       DO i=1,nstate
          IF(crge%f(i,1).LE.1.e-5_real_8) CYCLE
          nocc=nocc+1
          IF(geq0)THEN
             cnorm_l=c2_r(1,ibeg,i)**2*0.5_real_8
             gemax_l=ABS(c2_r(1,ibeg,i))**2
             gemax_ll=ABS(c2_r(1,ibeg,i))
          ELSE
             cnorm_l=c2_r(1,ibeg,i)**2+c2_r(2,ibeg,i)**2
             gemax_l=ABS(c2_r(1,ibeg,i))**2+ABS(c2_r(2,ibeg,i))**2
             gemax_ll=ABS(c2_r(1,ibeg,i))+ABS(c2_r(2,ibeg,i))
          END IF
          DO ig=ibeg+1,ngw_local
             cnorm_l=cnorm_l+c2_r(1,ig,i)**2+c2_r(2,ig,i)**2
             IF(gemax_ll.LT.ABS(c2_r(1,ig,i))+ABS(c2_r(2,ig,i)))THEN
                gemax_ll=ABS(c2_r(1,ig,i))+ABS(c2_r(2,ig,i))
                gemax_l=ABS(c2_r(1,ig,i))**2+ABS(c2_r(2,ig,i))**2
             END IF
          END DO
          IF(sp)THEN
             cnorm=cnorm+cnorm_l*crge%f(i,1)*2.0_real_8
             gemax=MAX(ABS(crge%f(i,1)*gemax_l),gemax)
          ELSE
             cnorm=cnorm+cnorm_l*2.0_real_8
             gemax=MAX(gemax_l,gemax)
          END IF
       END DO
    END IF
    gemax=SQRT(gemax)
  END SUBROUTINE csize_r

END MODULE csize_utils
