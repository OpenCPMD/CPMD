#include "cpmd_global.h"

MODULE vofrhob_utils
  USE cofor_utils,                     ONLY: cofor
  USE corec_utils,                     ONLY: corec
  USE cppt,                            ONLY: indz,&
                                             nzh
  USE cvan,                            ONLY: deeq
  USE efld,                            ONLY: extf,&
                                             textfld
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE extpotmod,                       ONLY: extpot
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE gcener_utils,                    ONLY: gcener,&
                                             gclsd
  USE geq0mod,                         ONLY: geq0
  USE graden_utils,                    ONLY: graden
  USE isos,                            ONLY: isos3
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE linres,                          ONLY: lrf4
  USE meta_localizespin_utils,         ONLY: localizespin
  USE newd_utils,                      ONLY: newd
  USE nlcc,                            ONLY: corel,&
                                             corer,&
                                             vnlt
  USE nlccstr_utils,                   ONLY: nlccstr
  USE nlps,                            ONLY: imagp
  USE nvtx_utils
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE reshaper,                        ONLY: reshape_inplace
  USE saop_utils,                      ONLY: gllb,&
                                             lb94m,&
                                             saop
  USE sfac,                            ONLY: fnl
  USE spin,                            ONLY: clsd,&
                                             spin_mod
  USE ssic_utils,                      ONLY: ssic
  USE str2,                            ONLY: dqg,&
                                             drhovg
  USE strs,                            ONLY: alpha,&
                                             beta,&
                                             delta,&
                                             dexc
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             ncpw,&
                                             parm,&
                                             spar
  USE tauf,                            ONLY: tau,&
                                             vtau
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE tpot,                            ONLY: c0y,&
                                             eigref,&
                                             eigval,&
                                             foccp,&
                                             mstate
  USE utils,                           ONLY: zgthr
  USE xcener_utils,                    ONLY: xcener
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vofrhob
  PUBLIC :: give_scr_vofrhob

CONTAINS

  ! ==================================================================
  SUBROUTINE vofrhob(fion,rhoe,v,vtemp,tfor,tstress)
    ! ==--------------------------------------------------------------==
    ! ==                        COMPUTES                              ==
    ! ==  THE ONE-PARTICLE POTENTIAL V IN REAL SPACE                  ==
    ! ==  THE TOTAL ENERGY ETOT                                       ==
    ! ==  THE FORCES FION ACTING ON THE IONS                          ==
    ! ==--------------------------------------------------------------==
    ! == RHOE:  in electronic density in real space                   ==
    ! ==       out total electronic potential in real space           ==
    ! ==--------------------------------------------------------------==
    REAL(real_8)                             :: fion(:,:,:), rhoe(:,:)
    COMPLEX(real_8)                          :: v(:,:), vtemp(:,:)
    LOGICAL                                  :: tfor, tstress

    CHARACTER(*), PARAMETER                  :: procedureN = 'vofrhob'

    COMPLEX(real_8), POINTER                 :: dqg_1d(:)
    INTEGER                                  :: ierr, ig, ir, isub, kk, nnrs, &
                                                nstates(2,1)
    INTEGER(int_8)                           :: il_grad(2), il_rhoval(2), il_rvtmp(2),&
                                                il_vpt1(1), il_vpt2(1), il_vtmp(2)
    LOGICAL                                  :: debug
    REAL(real_8)                             :: sgcc, sgcx, sxc, vgc
#ifdef _USE_SCRATCHLIBRARY
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: vtmp(:,:)
    REAL(real_8), POINTER __CONTIGUOUS       :: grad(:,:), rhoval(:,:), &
                                                rvtmp(:,:), vpt1(:), vpt2(:)
#else
    COMPLEX(real_8), ALLOCATABLE             :: vtmp(:,:)
    REAL(real_8), ALLOCATABLE                :: grad(:,:), rhoval(:,:), &
                                                rvtmp(:,:), vpt1(:), vpt2(:)
#endif
    REAL(real_8), DIMENSION(:, :, :), &
      POINTER                                :: fnla, fnladown, fnlaup

    CALL reshape_inplace(dqg, (/SIZE(dqg)/), dqg_1d)
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    __NVTX_TIMER_START ( procedureN )

    CALL setfftn(0)
    ! SCR partition.
    debug=.FALSE.
    il_vtmp(1) = ncpw%nhg
    il_vtmp(2) = MAX(fpar%nnr1,ncpw%nhg*2)*clsd%nlsx/ncpw%nhg      ! GCENER and COREC
    IF (cntl%tgc) THEN
       il_grad(1)=fpar%nnr1
       il_grad(2)=clsd%nlsd*4   ! GRADEN
    ELSE
       il_grad(1)=fpar%nnr1
       il_grad(2)=0
    ENDIF
    IF (corel%tinlc) THEN
       il_rhoval(1)=fpar%nnr1
       il_rhoval(2)=clsd%nlsd   ! XCENER
    ELSE
       il_rhoval(1)=fpar%nnr1
       il_rhoval(2)=0   ! XCENER
    ENDIF
    IF (cntl%tpotential) THEN
       il_vpt1(1)=fpar%nnr1*clsd%nlsd
       il_vpt2(1)=fpar%nnr1*clsd%nlsd
    ELSE
       il_vpt1(1)=0
       il_vpt2(1)=0
    ENDIF
    il_rvtmp(1)=fpar%nnr1
    il_rvtmp(2)=4
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_rvtmp,rvtmp,procedureN//'_rvtmp',ierr)
#else
    ALLOCATE(rvtmp(il_rvtmp(1),il_rvtmp(2)),STAT=ierr)!clsd%nlsx),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_vtmp,vtmp,procedureN//'_vtmp',ierr)
#else
    ALLOCATE(vtmp(il_vtmp(1), il_vtmp(2)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__) ! ! TODO il_vtmp should take into account that vtmp is complex
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_grad,grad,procedureN//'_grad',ierr)
#else
    ALLOCATE(grad(il_grad(1), il_grad(2)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF(corel%tinlc)THEN
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_rhoval,rhoval,procedureN//'_rhoval',ierr)
#else
       ALLOCATE(rhoval(il_rhoval(1), il_rhoval(2)),STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
    END IF
    ! ==--------------------------------------------------------------==
    ! == ADD CORE CHARGE TO RHOE                                      ==
    ! ==--------------------------------------------------------------==
    IF (corel%tinlc) THEN
       CALL dcopy(fpar%nnr1*clsd%nlsd,rhoe(1,1),1,rhoval(1,1),1)
       CALL corec(rhoe,vtmp,v)
       CALL dcopy(2*ncpw%nhg,vtemp(1,1),1,vnlt(1),1)
    ENDIF
    ! ==-------------------------------------------------------------==
    ! == VTEMP (Potential in G-Space) -FFT-> V(R)                     ==
    ! ==--------------------------------------------------------------==
    CALL zeroing(v)
#ifdef __NEC
    !CDIR NODEP
#endif
#ifdef __SR8000
    !poption parallel
#endif
#if defined(_vpp_) || defined(__PRIMERGY) || defined(__PRIMEHPC)
    !ocl novrec(v)
#endif
    !$omp parallel do private(IG)
    DO ig=1,ncpw%nhg
       v(indz(ig),1) = CONJG(vtemp(ig,1))
       v(nzh(ig),1)  = vtemp(ig,1)
    ENDDO
    IF (geq0) v(nzh(1),1) = vtemp(1,1)
    CALL invfftn(v(:,1),.FALSE.,parai%allgrp)
    ! ==--------------------------------------------------------------==
    ! == ADD EXTERNAL POTENTIAL TO V                                  ==
    ! ==--------------------------------------------------------------==
    IF (textfld) THEN
       ! ..potential from classical interface
       !$omp parallel do private(IR) shared(V)
       DO ir=1,fpar%nnr1
          v(ir,1)=v(ir,1)+CMPLX(extf(ir),0._real_8,kind=real_8)
       ENDDO
    ENDIF
    IF (cntl%texpot) THEN
       ! ..static external potential
       !$omp parallel do private(IR) shared(V)
       DO ir=1,fpar%nnr1
          v(ir,1)=v(ir,1)+CMPLX(extpot(ir),0._real_8,kind=real_8)
       ENDDO
    ENDIF
    IF ((cntl%texpot.OR.textfld).AND.corel%tinlc) THEN
       CALL fwfftn(v(:,1),.FALSE.,parai%allgrp)
       !$omp parallel do private(IG) shared(V)
#ifdef __SR8000
       !poption parallel
#endif
#ifdef _vpp_
       !OCL NOALIAS
#endif
       DO ig=1,ncpw%nhg
          vnlt(ig) = v(nzh(ig),1)
       ENDDO
       CALL invfftn(v(:,1),.FALSE.,parai%allgrp)
    ENDIF
    ! == ADD TO THE POTENTIAL THE CONTRIBUTION FROM SPIN LOCALIZATION ==
    CALL localizespin(rhoe,v)
    ! == ADD TO THE POTENTIAL THE SIC CONTRIBUTION ==
    IF (cntl%tsic) THEN
       IF (isos3%ps_type.EQ.1) THEN
          CALL stopgm('SSIC','HOCKNEY POISSON SOLVER NOT IMPLEMENTED',& 
               __LINE__,__FILE__)
       ELSE
          CALL ssic(ener_com%ehsic,rhoe,v)! cmb_ssic
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == COMPUTE EXCHANGE AND CORRELATION ENERGY (EXC)                ==
    ! ==--------------------------------------------------------------==
    IF (corel%tinlc) THEN
       CALL xcener(sxc,ener_com%vxc,rhoval,rhoe,v)
    ELSE
       CALL xcener(sxc,ener_com%vxc,rhoe,rhoe,v)
    ENDIF
    sgcx = 0.0_real_8
    sgcc = 0.0_real_8
    vgc  = 0.0_real_8
    IF (cntl%tgc) THEN
       ! ==------------------------------------------------------------==
       ! == CALCULATE THE GRADIENT OF THE DENSITY                      ==
       ! ==------------------------------------------------------------==
       IF (cntl%tlsd) THEN
          IF (tstress.OR.cntl%tdiag) THEN
             CALL dcopy(2*fpar%nnr1,v(1,1),1,dqg_1d(1),1)
             CALL dcopy(2*fpar%nnr1,v(1,2),1,dqg_1d(fpar%nnr1+1),1)
          ENDIF
          CALL fwfftn(v(:,1),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,v(:,1),vtemp(:,1),nzh)
          CALL fwfftn(v(:,2),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,v(:,2),vtemp(:,2),nzh)
          !$omp parallel do private(IR) shared(RHOE)
          DO ir=1,fpar%nnr1
             rhoe(ir,1)=rhoe(ir,1)-rhoe(ir,2)
          ENDDO
          CALL graden(rhoe(:,1),v,grad(:,1),vtmp)
          CALL graden(rhoe(:,2),v,grad(:,5),vtmp)
       ELSE
          IF (tstress.OR.cntl%tdiag) CALL dcopy(2*fpar%nnr1,v(1,1),1,dqg,1)
          CALL fwfftn(v(:,1),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,v,vtemp,nzh)
          CALL graden(rhoe,v,grad,vtmp)
       ENDIF
       ! ==--------------------------------------------------------------==
       ! == GRADIENT CORRECTION TO THE EXCHANGE ENERGY (EGCX)            ==
       ! ==--------------------------------------------------------------==
       IF (tstress .AND. cntl%ttau)&
            CALL stopgm(procedureN,'NO STRESS WITH TAU FUNCTIONAL',& 
            __LINE__,__FILE__)
       IF (cntl%tpotential) THEN
          IF (tstress) CALL stopgm(procedureN,'NO STRESS WITH TPOTENTIAL'&
               ,& 
               __LINE__,__FILE__)
          IF (pslo_com%tivan) CALL stopgm(procedureN,'NO VDB PP WITH TPOTENTIAL',& 
               __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
          CALL request_scratch(il_vpt1,vpt1,procedureN//'_vpt1',ierr)
#else
          ALLOCATE(vpt1(il_vpt1(1)),STAT=ierr)
#endif
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
          CALL request_scratch(il_vpt2,vpt2,procedureN//'_vpt2',ierr)
#else
          ALLOCATE(vpt2(il_vpt2(1)),STAT=ierr)
#endif
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
               __LINE__,__FILE__)
          CALL vpotential(rhoe,grad,v,vtemp,rvtmp,vpt1,vpt2)
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_vpt2,vpt2,procedureN//'_vpt2',ierr)
#else
          DEALLOCATE(vpt2,STAT=ierr)
#endif
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_vpt1,vpt1,procedureN//'_vpt1',ierr)
#else
          DEALLOCATE(vpt1,STAT=ierr)
#endif
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
               __LINE__,__FILE__)

       ELSEIF (cntl%tlsd) THEN
          ! IF(TSTRESS) CALL STOPGM('VOFRHOB','NO STRESS WITH LSD AND GC')
          CALL gclsd(sgcx,sgcc,rhoe,v,vtemp,rvtmp,grad,tstress)
          IF (tstress.OR.cntl%tdiag) THEN
             CALL daxpy(2*fpar%nnr1,-1.0_real_8,v(1,1),1,dqg_1d(1),1)
             CALL daxpy(2*fpar%nnr1,-1.0_real_8,v(1,2),1,dqg_1d(fpar%nnr1+1),1)
             CALL dscal(4*fpar%nnr1,-1.0_real_8,dqg,1)
             vgc=0.0_real_8
             IF (corel%tinlc) THEN
#if defined(__VECTOR)
                !$omp parallel do private(IR) reduction(+:VGC)
#else
                !$omp parallel do private(IR) reduction(+:VGC) schedule(static)
#endif
                DO ir=1,fpar%nnr1
                   vgc=vgc+(rhoval(ir,1)-rhoval(ir,2))*REAL(dqg_1d(ir))
                   vgc=vgc+rhoval(ir,2)*REAL(dqg_1d(fpar%nnr1+ir))
                ENDDO
             ELSE
#if defined(__VECTOR)
                !$omp parallel do private(IR) reduction(+:VGC)
#else
                !$omp parallel do private(IR) reduction(+:VGC) schedule(static)
#endif
                DO ir=1,fpar%nnr1
                   vgc=vgc+rhoe(ir,1)*REAL(dqg_1d(ir))
                   vgc=vgc+rhoe(ir,2)*REAL(dqg_1d(fpar%nnr1+ir))
                ENDDO
             ENDIF
             IF (cntl%ttau) THEN
#if defined(__VECTOR)
                !$omp parallel do private(IR) reduction(+:VGC)
#else
                !$omp parallel do private(IR) reduction(+:VGC) schedule(static)
#endif
                DO ir=1,fpar%nnr1
                   vgc=vgc+(tau(ir,1)*vtau(ir,1)+tau(ir,2)*vtau(ir,2))
                ENDDO
             ENDIF
          ENDIF
       ELSE
          CALL gcener(sgcx,sgcc,rhoe,v,vtemp,rvtmp(:,1),grad,tstress)
          IF (tstress.OR.cntl%tdiag) THEN
             CALL daxpy(2*fpar%nnr1,-1.0_real_8,v(1,1),1,dqg,1)
             CALL dscal(2*fpar%nnr1,-1.0_real_8,dqg,1)
             vgc=0.0_real_8
             IF (corel%tinlc) THEN
#if defined(__VECTOR)
                !$omp parallel do private(IR) reduction(+:VGC)
#else
                !$omp parallel do private(IR) reduction(+:VGC) schedule(static)
#endif
#ifdef __SR8000
                !poption parallel, tlocal(IR), psum(VGC)
#endif 
                DO ir=1,fpar%nnr1
                   vgc=vgc+rhoval(ir,1)*REAL(dqg_1d(ir))
                ENDDO
             ELSE
#if defined(__VECTOR)
                !$omp parallel do private(IR) reduction(+:VGC)
#else
                !$omp parallel do private(IR) reduction(+:VGC) schedule(static)
#endif
#ifdef __SR8000
                !poption parallel, tlocal(IR), psum(VGC)
#endif 
                DO ir=1,fpar%nnr1
                   vgc=vgc+rhoe(ir,1)*REAL(dqg_1d(ir))
                ENDDO
             ENDIF
             IF (cntl%ttau) THEN
#if defined(__VECTOR)
                !$omp parallel do private(IR) reduction(+:VGC)
#else
                !$omp parallel do private(IR) reduction(+:VGC) schedule(static)
#endif
#ifdef __SR8000
                !poption parallel, tlocal(IR), psum(VGC)
#endif
                DO ir=1,fpar%nnr1
                   vgc=vgc+tau(ir,1)*vtau(ir,1)
                ENDDO
             ENDIF
          ENDIF
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == XC CONTRIBUTION TO STRESS (VANDERBILT CHARGE)                ==
    ! ==--------------------------------------------------------------==
    IF (pslo_com%tivan.AND.tstress) THEN
       CALL xcstr(drhovg,rhoe,dqg_1d(1),dqg_1d(maxfftn+1))
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == NLCC CONTRIBUTION TO STRESS                                  ==
    ! ==--------------------------------------------------------------==
    IF (corel%tinlc.AND.tstress) THEN
       CALL nlccstr(dqg_1d(:),dqg_1d(clsd%nlsd*maxfftn+1:))
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == RELEASE TEMPORARY ARRAYS FROM SCRATCH SPACE                  ==
    ! ==--------------------------------------------------------------==
    IF(corel%tinlc)THEN
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_rhoval,rhoval,procedureN//'_rhoval',ierr)
#else
       DEALLOCATE(rhoval,STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_grad,grad,procedureN//'_grad',ierr)
#else
    DEALLOCATE(grad,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_vtmp,vtmp,procedureN//'_vtmp',ierr)
#else
    DEALLOCATE(vtmp,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_rvtmp,rvtmp,procedureN//'_rvtmp',ierr)
#else
    DEALLOCATE(rvtmp,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    ! == V CONTAINS THE TOTAL POTENTIAL IN R-SPACE                    ==
    ! == MOVE IT TO RHOE                                              ==
    ! ==--------------------------------------------------------------==
#if defined (__VECTOR)
    !$omp parallel do private(IR)
#else
    !$omp parallel do private(IR) schedule(static)
#endif
    DO ir=1,fpar%nnr1
       rhoe(ir,1)=REAL(v(ir,1))
    ENDDO
    IF (cntl%tlsd) THEN
#if defined (__VECTOR)
       !$omp parallel do private(IR)
#else
       !$omp parallel do private(IR) schedule(static)
#endif
       DO ir=1,fpar%nnr1
          rhoe(ir,2)=REAL(v(ir,2))
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == CORRECT ENERGY FOR CORE CHARGE                               ==
    ! ==--------------------------------------------------------------==
    IF (corel%tinlc) THEN
       sxc=sxc-corer%excco
       sgcx=sgcx-corer%egcxco
       sgcc=sgcc-corer%egccco
    ENDIF
    IF (pslo_com%tivan.OR.(corel%tinlc.AND.tfor)) THEN
       ! ==--------------------------------------------------------------==
       ! == TRANSF. TOTAL POTENTIAL IN G-SPACE                           ==
       ! == PUT THE TOTAL POTENTIAL IN G-SPACE INTO VTEMP                ==
       ! ==--------------------------------------------------------------==
       CALL fwfftn(v(:,1),.FALSE.,parai%allgrp)
       CALL zgthr(ncpw%nhg,v,vtemp,nzh)
       IF (cntl%tlsd) THEN
          CALL fwfftn(v(:,2),.FALSE.,parai%allgrp)
          CALL zgthr(ncpw%nhg,v(:,2),vtemp(:,2),nzh)
       ENDIF
    ENDIF
    IF (pslo_com%tivan) THEN
       IF (imagp.EQ.2)CALL stopgm(procedureN,'K-POINT NOT IMPLEMENTED',& 
            __LINE__,__FILE__)
       ! ==--------------------------------------------------------------==
       ! == CALCULATE DEEQ FOR VANDERBILT PP                             ==
       ! == FORCE ON IONS DUE TO THE "VANDERBILT CHARGE"                 ==
       ! ==--------------------------------------------------------------==
       IF (cntl%tlsd) THEN
          nstates(1,1)=1
          nstates(2,1)=spin_mod%nsup
          CALL newd(deeq(:,:,:,1),crge%f,vtemp(:,1),fion,tfor,nstates,.FALSE.)
          nstates(1,1)=spin_mod%nsup+1
          nstates(2,1)=spin_mod%nsup+spin_mod%nsdown
          CALL newd(deeq(:,:,:,2),crge%f,vtemp(:,2),fion,tfor,nstates,.FALSE.)
       ELSE
          nstates(1,1)=1
          nstates(2,1)=crge%n
          CALL newd(deeq(:,:,:,1),crge%f,vtemp,fion,tfor,nstates,.FALSE.)
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! == FORCE ON IONS DUE TO NLCC                                    ==
    ! ==--------------------------------------------------------------==
    IF (corel%tinlc.AND.tfor) CALL cofor(fion,vtemp)
    ! ==--------------------------------------------------------------==
    ! == SINGLE PIECES OF THE ENERGY:                                 ==
    ! ==--------------------------------------------------------------==
    nnrs = spar%nr1s*spar%nr2s*spar%nr3s
    ener_com%egc=(sgcc+sgcx)*parm%omega/REAL(nnrs,kind=real_8)
    ener_com%exc=sxc*parm%omega/REAL(nnrs,kind=real_8) + ener_com%egc
    ener_com%vxc=(ener_com%vxc + vgc)*parm%omega/REAL(nnrs,kind=real_8)
    IF (tstress) THEN
       DO kk=1,6
          dexc(kk) = dexc(kk) + (ener_com%exc-ener_com%vxc)*delta(alpha(kk),beta(kk))
       ENDDO
    ENDIF

    __NVTX_TIMER_STOP
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE vofrhob
  ! ==================================================================
  SUBROUTINE give_scr_vofrhob(lvofrhob,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lvofrhob
    CHARACTER(len=30)                        :: tag

    INTEGER                                  :: ltgc, ltinlc, lvtmp

    IF (cntl%tgc) THEN
       ltgc=fpar%nnr1*clsd%nlsd*4      ! GRADEN
    ELSE
       ltgc=0
    ENDIF
    IF (corel%tinlc) THEN            ! XCENER
       ltinlc=fpar%nnr1*clsd%nlsd
    ELSE
       ltinlc=0
    ENDIF
    lvtmp=MAX(fpar%nnr1,ncpw%nhg*2)*clsd%nlsx
    IF (cntl%tpotential) THEN
       lvtmp=lvtmp+2*fpar%nnr1*clsd%nlsd
    ENDIF
    lvofrhob = lvtmp+ltgc+ltinlc ! GCENER and COREC
    lvofrhob = lvofrhob + 100  ! For boundary checks in SCRPTR
    tag='MAX(NNR1,NHG*2)*NLSX+.....'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_vofrhob
  ! ==================================================================
  SUBROUTINE vpotential(rhoe,grad,psi,vtemp,v1,v2,v3)
    ! ==--------------------------------------------------------------==

    REAL(real_8)                             :: rhoe(:,:), grad(:,:)
    COMPLEX(real_8)                          :: psi(:,:), vtemp(:,:)
    REAL(real_8)                             :: v1(:,:), v2(:), v3(:)

    INTEGER                                  :: ig, ir
    REAL(real_8)                             :: exc

! ==--------------------------------------------------------------==

    IF (cntl%tlsd) THEN
       !$omp parallel do private(IR)
       DO ir=1,fpar%nnr1
          v2(ir) = grad(ir,2)*grad(ir,6)+grad(ir,3)*grad(ir,7)+&
               grad(ir,4)*grad(ir,8)
       ENDDO
       CALL dcopy(fpar%nnr1,grad(1,5),1,grad(1,3),1)
       CALL dcopy(fpar%nnr1,v2,1,grad(1,2),1)
    ENDIF
    CALL zeroing(v1)!,nnr1*clsd%nlsd)
    IF (lrf4%td_functional.EQ.10) THEN
       CALL saop(c0y,eigval,eigref,foccp,mstate,rhoe,grad,exc,v1,v2,v3,&
            grad(:,4),psi(:,1))
    ELSEIF (lrf4%td_functional.EQ.11) THEN
       CALL lb94m(rhoe,grad,exc,v1)
    ELSEIF (lrf4%td_functional.EQ.12) THEN
       CALL gllb(c0y,eigval,eigref,foccp,mstate,rhoe,grad,exc,&
            v1,v2,psi(:,1))
    ELSE
       CALL stopgm("VPOTENTIAL","ILLEGAL TD_FUNCTIONAL",& 
            __LINE__,__FILE__)
    ENDIF
    ! 
    IF (cntl%tlsd) THEN
       CALL zeroing(psi(:,1))!,maxfft)
       !ocl novrec
       !$omp parallel do private(IG) shared(PSI)
       DO ig=1,ncpw%nhg
          psi(nzh(ig),1) = vtemp(ig,1)
          psi(indz(ig),1) = CONJG(vtemp(ig,1))
       ENDDO
       CALL  invfftn(psi(:,1),.FALSE.,parai%allgrp)
       !$omp parallel do private(IR) shared(PSI)
       DO ir=1,fpar%nnr1
          psi(ir,1)=psi(ir,1)+v1(ir,1)
       ENDDO
       CALL zeroing(psi(:,2))!,maxfft)
       !ocl novrec
       !$omp parallel do private(IG) shared(PSI)
       DO ig=1,ncpw%nhg
          psi(nzh(ig),2) = vtemp(ig,2)
          psi(indz(ig),2) = CONJG(vtemp(ig,2))
       ENDDO
       CALL  invfftn(psi(:,2),.FALSE.,parai%allgrp)
       !$omp parallel do private(IR) shared(PSI)
       DO ir=1,fpar%nnr1
          psi(ir,2)=psi(ir,2)+v1(ir,2)
       ENDDO
    ELSE
       CALL zeroing(psi(:,1))!,maxfft)
       !ocl novrec
       !$omp parallel do private(IG) shared(PSI)
       DO ig=1,ncpw%nhg
          psi(nzh(ig),1) = vtemp(ig,1)
          psi(indz(ig),1) = CONJG(vtemp(ig,1))
       ENDDO
       CALL  invfftn(psi(:,1),.FALSE.,parai%allgrp)
       !$omp parallel do private(IR) shared(PSI)
       DO ir=1,fpar%nnr1
          psi(ir,1)=psi(ir,1)+v1(ir,1)
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE vpotential
  ! ==================================================================

END MODULE vofrhob_utils
