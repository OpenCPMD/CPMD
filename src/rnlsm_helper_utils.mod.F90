#include "cpmd_global.h"

MODULE rnlsm_helper
  USE beta_utils,                      ONLY: build_beta
  USE distribution_utils,              ONLY: dist_atoms
  USE error_handling,                  ONLY: stopgm
  USE gpu
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nlps_com
  USE system,                          ONLY: cntl,&
                                             cnti,&
                                             maxsys,&
                                             parm
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: tune_rnlsm
  PUBLIC :: proj_beta
  PUBLIC :: set_buffers
CONTAINS
  ! Author: Tobias Kloeffel, FAU Erlangen-Nuernberg April 2019
  ! Contains some helper routines for rnlsm1/2 to avoid duplicated code
  ! ==================================================================
  SUBROUTINE tune_rnlsm(count,first,others,timings,work,nstate)
    ! ==--------------------------------------------------------------==
    ! == routine to calculate dgemm chunk sizes for effective         ==
    ! == overlapping of communication and computation                 ==
    ! == currently used in rnlsm1/2                                   ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: work,nstate
    REAL(real_8),INTENT(IN)                  :: timings(2)
    INTEGER,INTENT(OUT)                      :: count
    REAL(real_8),INTENT(OUT)                 :: first,others
    INTEGER                                  :: tot_work

    IF(work.LE.0) THEN
       count=0
       RETURN
    END IF
    !computation .GT. communication?
    IF(timings(1).GT.timings(2))THEN
       tot_work=work
       !work for overlapping part
       tot_work=CEILING(REAL(tot_work,kind=real_8)*timings(2)/timings(1))
       !take at least 300 in each buffer
       count=CEILING(REAL(tot_work,KIND=real_8)/300_real_8)
       !failsave for small amount of work
       IF(count.EQ.1)THEN
          count=2
       ELSEIF(count.EQ.0)THEN
          count=1
       ELSEIF(count.GT.15)THEN
          count=15
       END IF
       IF(count.GT.1) THEN
          first=1.0_real_8-timings(2)/timings(1)
          others=timings(2)/timings(1)/REAL(count-1,kind=real_8)
       ELSE
          first=1.0_real_8
          others=0.0_real_8
       END IF
    ELSE
       !communication has won.. set arbitrary number of buffers
       tot_work=work
       count=CEILING(REAL(tot_work,KIND=real_8)/300_real_8)
       !prevent more than 15 buffers
       IF(count.GT.15) count=15
       first=1.0_real_8/REAL(count,KIND=real_8)
       others=1.0_real_8/REAL(count,KIND=real_8)
    END IF
  END SUBROUTINE tune_rnlsm
  ! ==================================================================
  SUBROUTINE set_buffers(autotune_it,maxbuff,timings,bc,b1,b2,na_in,na_out,ld,start,&
       tot_work,nstate,k,disable)
    ! ==--------------------------------------------------------------==
    ! == returns custom na                                            ==
    ! == handles autotuning stuff                                     ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: na_in(:,:), maxbuff, nstate, k
    INTEGER,INTENT(OUT),ALLOCATABLE          :: na_out(:,:,:)
    INTEGER,INTENT(INOUT)                    :: autotune_it, bc
    REAL(real_8),INTENT(INOUT)               :: timings(2), b1, b2
    INTEGER,INTENT(OUT)                      :: tot_work, ld(maxbuff), start(maxbuff)
    LOGICAL,INTENT(IN)                       :: disable
    INTEGER                                  :: buff
    CHARACTER(*), PARAMETER                  :: procedureN = 'set_buffers'
    INTEGER, ALLOCATABLE                     :: map_isa(:,:)
    REAL(real_8)                             :: fractions(maxbuff)
    INTEGER                                  :: is,ierr
    
    tot_work=0
    DO is=1,ions1%nsp
       tot_work=tot_work+(na_in(2,is)-na_in(1,is)+1)*nlps_com%ngh(is)
    END DO
    tot_work=tot_work*k
    IF(cntl%overlapp_comm_comp)THEN
       !skip splitting if there is too less work...
       IF(tot_work.GT.100.AND..NOT.disable)THEN
          IF (autotune_it.EQ.0.AND.cnti%rnlsm_autotune_maxit.GT.0) THEN
             timings=0.0_real_8
             bc=1
          END IF
          !check if autotuning is active, count iteration
          IF (autotune_it.LT.cnti%rnlsm_autotune_maxit) autotune_it=autotune_it+1
       END IF
    ELSE
       bc=1
    END IF
    IF(disable) bc=1
    IF(bc.GT.1)THEN
       ALLOCATE(map_isa(ions1%nat+1,bc), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate map_isa',& 
            __LINE__,__FILE__)
       ALLOCATE(na_out(2,ions1%nsp,bc), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate na_out',& 
            __LINE__,__FILE__)
       map_isa=0
       fractions(1)=b1
       fractions(2:bc)=b2
       CALL dist_atoms(bc,fractions,na_in,na_out,map_isa)
       start=0
       ld=0
       start(1)=1
       ld(1)=map_isa(ions1%nat+1,1)*k
       DO buff=2,bc
          ld(buff)=map_isa(ions1%nat+1,buff)*k
          start(buff)=start(buff-1)+ld(buff-1)*nstate
       END DO
       DEALLOCATE(map_isa, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate map_isa',& 
            __LINE__,__FILE__)
    ELSE
       ALLOCATE(na_out(2,ions1%nsp,1), stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate na_out',& 
            __LINE__,__FILE__)
       b1=1.0_real_8
       b2=0.0_real_8
       na_out(:,:,1)=na_in(:,:)
       start=0
       start(1)=1
       ld=0
       ld(1)=tot_work
    END IF
    RETURN
  END SUBROUTINE set_buffers
  ! ==================================================================
  SUBROUTINE proj_beta(na,igeq0,nstate,c0,ld_c0,eigkr,eiscr,ld_eiscr,startg,&
       dai,ld_dai,deriv,tkpnt,geq0,twnl_nghtol,twnl_nghtol_gk)
    ! ==--------------------------------------------------------------==
    ! == Calculates the local projection of beta on C0                ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: na(2,*), nstate, startg, ld_eiscr, ld_c0, &
                                                ld_dai, igeq0
    COMPLEX(real_8),INTENT(IN)               :: c0(ld_c0,*),eigkr(ld_c0,*)
    COMPLEX(real_8),INTENT(OUT)              :: eiscr(ld_eiscr,*)
    REAL(real_8),INTENT(OUT)                 :: dai(ld_dai,*)
    REAL(real_8),INTENT(IN),OPTIONAL         :: twnl_nghtol(ld_c0,maxsys%nhxs,ions1%nsp)
    REAL(real_8),INTENT(IN),OPTIONAL         :: twnl_nghtol_gk(ld_c0,3,maxsys%nhxs,ions1%nsp)
    LOGICAL,INTENT(IN)                       :: deriv, tkpnt, geq0
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8,0._real_8) ,&
                                                zzero = (0._real_8,0._real_8)
    COMPLEX(real_8)                          :: zfac
    REAL(real_8)                             :: tfac
    CHARACTER(*), PARAMETER                  :: procedureN = 'proj_beta'
    INTEGER                                  :: isub1,isub2,isub3,isub4
    
    IF(deriv)THEN
       CALL tiset(procedureN//'build_beta_deriv',isub1)
       CALL build_beta(na,eigkr,eiscr,ld_eiscr,startg,ld_c0,ld_dai,igeq0,geq0,tkpnt,twnl_nghtol_gk=twnl_nghtol_gk)
       CALL tihalt(procedureN//'build_beta_deriv',isub1)
    ELSE
       CALL tiset(procedureN//'build_beta',isub2)
       CALL build_beta(na,eigkr,eiscr,ld_eiscr,startg,ld_c0,ld_dai,igeq0,geq0,tkpnt,twnl_nghtol=twnl_nghtol)
       CALL tihalt(procedureN//'build_beta',isub2)
    END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target update from(eiscr(:,:ld_dai)) if(update_first_to_gpu)
#endif
    IF(deriv)THEN
       CALL tiset(procedureN//'proj_beta_deriv',isub3)
    ELSE
       CALL tiset(procedureN//'proj_beta',isub4)
    END IF
    IF(TKPNT) THEN       
       IF(deriv)THEN 
          zfac=parm%tpiba*zone
       ELSE
          zfac=zone
       END IF
       CALL cpmd_zgemm('C','N',ld_dai,nstate,  ld_eiscr,zfac,eiscr,  ld_eiscr,c0(startg,1),  &
              ld_c0, zzero,    dai,ld_dai)
    ELSE
       IF(deriv)THEN
          tfac=2._real_8*parm%tpiba
       ELSE
          tfac=2._real_8
       END IF
       CALL cpmd_dgemm('T','N',ld_dai,nstate,2*ld_eiscr,tfac,eiscr,2*ld_eiscr,c0(startg,1),  &
            2*ld_c0, 0._real_8,dai,ld_dai)
    ENDIF
    IF(deriv)THEN
       CALL tihalt(procedureN//'proj_beta_deriv',isub3)
    ELSE
       CALL tihalt(procedureN//'proj_beta',isub4)
    END IF

  END SUBROUTINE proj_beta

END MODULE rnlsm_helper
