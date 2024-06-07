#include "cpmd_global.h"

MODULE fftprp_utils
  USE benc,                            ONLY: ibench
  USE cp_cuda_types,                   ONLY: cp_cuda_devices_fft,&
                                             cp_cuda_env
  USE cp_cufft_types,                  ONLY: cp_cufft
  USE cp_cufft_utils,                  ONLY: cp_cufft_finalize_devices,&
                                             cp_cufft_init_devices
  USE cp_curho_types,                  ONLY: cp_curho
  USE cp_curho_utils,                  ONLY: cp_curho_finalize,&
                                             cp_curho_init
  USE cppt,                            ONLY: indz,&
                                             indzs,&
                                             inyh,&
                                             nzh,&
                                             nzhs
  USE cuda_utils,                      ONLY: cuda_alloc_host,cuda_host_register, cuda_host_unregister, cuda_dealloc_host
  USE elct,                            ONLY: crge
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: &
       fftpool, kr1m, kr2max, kr2min, kr3max, kr3min, lmsqmax, lnzf, lnzs, &
       maxrpt, mg, ms, msp, mxy, mz, ngrm, nhrm, nr1m, nr3m, xf, yf,&
       batch_fft, a2a_msgsize, fft_batchsize, fft_numbatches, fft_residual, &
       fft_total,  fft_tune_max_it, fft_tune_num_it, fft_time_total, fft_batchsizes, fft_min_numbatches
  USE fft_maxfft,                      ONLY: maxfft
  USE fftnew_utils,                    ONLY: addfftnset,&
                                             setfftn
  USE isos,                            ONLY: isos1
  USE kinds,                           ONLY: real_8
  USE kpts,                            ONLY: tkpts
  USE mp_interface,                    ONLY: mp_max,&
                                             mp_sum,&
                                             mp_bcast
  USE parac,                           ONLY: parai,&
                                             paral
  USE part_1d,                         ONLY: part_1d_get_blk_bounds
  USE prmem_utils,                     ONLY: prmem
  USE sizeof_kinds,                    ONLY: sizeof_complex_8
  USE reshaper,                        ONLY: reshape_inplace
  USE rswfmod,                         ONLY: lwdim,&
                                             maxstates,&
                                             rsactive,&
                                             rswf
  USE system,                          ONLY: cntl,&
                                             cntr,&
                                             cnti,&
                                             fpar,&
                                             group,&
                                             ncpw,&
                                             nkpt,&
                                             parap,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: icopy
  USE zeroing_utils,                   ONLY: zeroing
  USE string_utils,                    ONLY: int2str

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: fft_init
  PUBLIC :: fft_finalize
  PUBLIC :: autotune_fftbatchsize

  LOGICAL, PARAMETER, PRIVATE :: cuda_register_memory = .true.

CONTAINS


  SUBROUTINE fft_init ( )
    CHARACTER(*), PARAMETER                  :: procedureN = 'fft_init'

    INTEGER                                  :: cp_ipool, isub

! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)

    CALL fftprp_default_init ( )

    !vw   need to initialized the new fft before they use
    CALL addfftnset(-1._real_8,-1._real_8,cp_ipool)

    ! fftprp_default needs to be called before the cuda initialization
    IF( cp_cuda_env%use_fft ) THEN
       CALL cp_cufft_init_devices ( cp_cuda_env, cp_cuda_devices_fft, cp_cufft )

       CALL cp_curho_init( cp_cuda_env, cp_cuda_devices_fft, cp_curho )
       !CALL cp_curho_alloc_buffers ( cp_curho, fpar%nnr1 ) !vw for the moment only for RKS. need to pass the extra dim.
    ENDIF

    !vw init fft arrays, this should come after fft_init_cufft_*
    CALL setfftn(0)

    ! ==--------------------------------------------------------------==
    CALL tihalt(procedureN,isub)
  END SUBROUTINE fft_init


  SUBROUTINE fft_finalize ( )
    CHARACTER(*), PARAMETER                  :: procedureN = 'fft_finalize'

    INTEGER                                  :: isub

    CALL tiset(procedureN,isub)

    IF( cp_cuda_env%use_fft ) THEN
       CALL cp_cufft_finalize_devices ( cp_cufft )
       !CALL cp_curho_dealloc_buffers ( cp_curho )
       CALL cp_curho_finalize ( cp_curho )
    ENDIF

    CALL fftprp_default_finalize ( )

    CALL tihalt(procedureN,isub)
  END SUBROUTINE fft_finalize


  ! ==================================================================
  SUBROUTINE fftprp_default_init
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'fftprp_default_init'

    INTEGER :: i, ierr, ig, ij, img, iny1, iny2, iny3, ip, ipp, ixf, j, jj, &
      jmg, ldim, len, mxrp, nclu, ngray, nh1, nh2, nh3, nhray, nl1, nl2, nn2, &
      nr3i, nrx, nstate, ny1, ny2, ny3, first, last, lda, it, count
    INTEGER, ALLOCATABLE                     :: my(:)
    REAL(real_8)                             :: rmem, rstate, xmpenm

! ==--------------------------------------------------------------==
! GATHER ARRAY FOR FFT ALONG X
! SPARSITY FOR FFT ALONG Y

    ALLOCATE(mg(fpar%kr2s,fpar%kr3s),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(mz((2*spar%nr3s)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(my((2*spar%nr2s)),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(mg)!,kr2s*kr3s)
    CALL zeroing(mz)!,2*spar%nr3s)
    CALL zeroing(my)!,2*spar%nr2s)
    nh1=spar%nr1s/2+1
    nh2=spar%nr2s/2+1
    nh3=spar%nr3s/2+1
    DO ig=1,ncpw%ngw
       ny2=inyh(2,ig)
       ny3=inyh(3,ig)
       iny2=-ny2+2*nh2
       iny3=-ny3+2*nh3
       mg(ny2,ny3)=mg(ny2,ny3)+1
       mg(iny2,iny3)=mg(iny2,iny3)+1
       my(ny2)=my(ny2)+1
       my(iny2)=my(iny2)+1
       mz(ny3)=mz(ny3)+1
       mz(iny3)=mz(iny3)+1
    ENDDO
    CALL mp_sum(my,my(spar%nr2s+1:),spar%nr2s,parai%allgrp)
    CALL icopy(spar%nr2s,my(spar%nr2s+1),1,my,1)
    CALL mp_sum(mz,mz(spar%nr3s+1:),spar%nr3s,parai%allgrp)
    CALL icopy(spar%nr3s,mz(spar%nr3s+1),1,mz,1)
    kr2min=1
    DO i=1,fpar%kr2s
       IF (my(i).NE.0) THEN
          kr2min=i
          GOTO 50
       ENDIF
    ENDDO
50  CONTINUE
    kr2max=fpar%kr2s
    DO i=fpar%kr2s,1,-1
       IF (my(i).NE.0) THEN
          kr2max=i
          GOTO 51
       ENDIF
    ENDDO
51  CONTINUE
    kr3min=1
    DO i=1,fpar%kr3
       IF (mz(i).NE.0) THEN
          kr3min=i
          GOTO 52
       ENDIF
    ENDDO
52  CONTINUE
    kr3max=fpar%kr3s
    DO i=fpar%kr3,1,-1
       IF (mz(i).NE.0) THEN
          kr3max=i
          GOTO 53
       ENDIF
    ENDDO
53  CONTINUE
    ALLOCATE(nzhs(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(indzs(ncpw%ngw),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target enter data map(alloc:nzhs,indzs)
#endif
    ! MAXIMUM OF NR1, NGRAYS AND NHRAYS FOR MP_INDEX
    nr1m = 0
    nhrm = 0
    ngrm = 0
    DO i=0,parai%nproc-1
       nr1m = MAX(nr1m,parap%sparm(5,i))
       nhrm = MAX(nhrm,parap%sparm(8,i))
       ngrm = MAX(ngrm,parap%sparm(9,i))
    ENDDO
    kr1m=MAX(nr1m+MOD(nr1m+1,2),fpar%kr1)
    img=0
    DO j=1,fpar%kr3s
       DO i=1,fpar%kr2s
          IF (mg(i,j).NE.0) THEN
             img=img+1
             mg(i,j)=img
          ENDIF
       ENDDO
    ENDDO
    ngray=img
    DO ig=ncpw%ngw+1,ncpw%nhg
       ny2=inyh(2,ig)
       ny3=inyh(3,ig)
       iny2=-ny2+2*nh2
       iny3=-ny3+2*nh3
       jmg=mg(ny2,ny3)
       IF (jmg.EQ.0) mg(ny2,ny3)=-1
       jmg=mg(iny2,iny3)
       IF (jmg.EQ.0) mg(iny2,iny3)=-1
    ENDDO
    DO j=1,fpar%kr3s
       DO i=1,fpar%kr2s
          IF (mg(i,j).LT.0) THEN
             img=img+1
             mg(i,j)=img
          ENDIF
       ENDDO
    ENDDO
    nhray=img
    ! SCATTER ARRAY FOR FFT ALONG X
    ALLOCATE(ms(nhrm,2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(ms)!,2*nhrm)
    DO i=1,fpar%kr2s
       DO j=1,fpar%kr3s
          ij=mg(i,j)
          IF (ij.GT.0) THEN
             ms(ij,1)=i
             ms(ij,2)=j
          ENDIF
       ENDDO
    ENDDO
    ! CONCATENATE GATHER/SCATTER ARRAYS
    ALLOCATE(msp(nhrm,2,parai%nproc),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target enter data map(alloc:msp)
#endif
    len = nhrm * 8
    CALL my_concat(ms,msp,len,parai%allgrp)
    ! TRANSLATE I,J TO A SINGLE G/S INDEX
    DO ip=0,parai%nproc-1
       mxrp=parap%sparm(8,ip)
       DO ixf=1,mxrp
          jj=ixf+ip*2*nhrm
          i=msp(ixf,1,ip+1)
          j=msp(ixf,2,ip+1)
          msp(ixf,1,ip+1)=i+(j-1)*fpar%kr2s
          IF (ixf.LE.parap%sparm(9,ip)) msp(ixf,2,ip+1)=i+(j-kr3min)*fpar%kr2s
       ENDDO
    ENDDO
    ! REDEFINE NZH AND INDZ FOR COMPRESSED STORAGE
    nn2=1
    DO ig=1,ncpw%nhg
       ny1=inyh(1,ig)
       ny2=inyh(2,ig)
       ny3=inyh(3,ig)
       iny1=-ny1+2*nh1
       iny2=-ny2+2*nh2
       iny3=-ny3+2*nh3
       nzh(ig)=ny1 + (mg(ny2,ny3)-1)*fpar%kr1s
       indz(ig)=iny1 + (mg(iny2,iny3)-1)*fpar%kr1s
    ENDDO
    !$omp parallel do private(IG)
    DO ig=1,ncpw%ngw
       nzhs(ig)=nzh(ig)
       indzs(ig)=indz(ig)
    ENDDO
    ! Some dimensions used for groups
    fpar%krx=1
    IF (group%nogrp.GT.1) THEN
       group%mpen=0
       DO ip=1,group%nogrp
          ipp=parap%nlink(group%nolist(ip))
          group%mpen=group%mpen+parap%sparm(9,ipp)
       ENDDO
       xmpenm=group%mpen
       CALL mp_max(xmpenm,parai%allgrp)
       group%mpenm=NINT(xmpenm)
       nl1=parap%nlink(group%nolist(1))
       nl2=parap%nlink(group%nolist(group%nogrp))
       nrx=parap%nrxpl(nl1,2)-parap%nrxpl(nl2,1)+1
       fpar%krx=nrx+MOD(nrx+1,2)
       CALL grpgs
    ELSE
       DEALLOCATE(ms,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ENDIF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target update to (nzhs,indzs,msp)
#endif
    ! ARRAY SIZE TO DO A 3D-FFT
    maxfft = MAX(kr1m*fpar%kr2s*fpar%kr3s,parai%nproc*nr1m*nhrm)
    IF (group%nogrp.GT.1) maxfft = MAX(fpar%krx*fpar%kr2s*fpar%kr3s,maxfft)
    IF (isos1%tclust) THEN
       nr3m = 0
       DO i=0,parai%nproc-1
          nr3i = parap%nrzpl(i,2)-parap%nrzpl(i,1)+1
          nr3m = MAX(nr3m,nr3i)
       ENDDO
       nclu = parai%nproc * nr1m*fpar%kr2s*nr3m
       maxfft = MAX(nclu,maxfft)
    ENDIF
    IF( cp_cuda_env%use_fft ) THEN
#if defined(_HAS_CUDA)
       block         
         use sizeof_kinds, ONLY: sizeof_complex_8
         use machine, only: m_getpagesize
         integer :: ps
         ps = m_getpagesize()
         !vw resize maxfft to page size (needed by cuda?)
         maxfft = CEILING(REAL(maxfft,real_8) / REAL( ps / sizeof_complex_8) ) * ps / sizeof_complex_8
       end block
       IF(cuda_register_memory) THEN
          ALLOCATE(xf(maxfft,cp_cuda_env%fft_n_devices_per_task*cp_cuda_env%fft_n_streams_per_device),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
          ALLOCATE(yf(maxfft,cp_cuda_env%fft_n_devices_per_task*cp_cuda_env%fft_n_streams_per_device),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
          CALL cuda_host_register(xf)
          CALL cuda_host_register(yf)
       ELSE
          CALL cuda_alloc_host(xf,[maxfft,cp_cuda_env%fft_n_devices_per_task*cp_cuda_env%fft_n_streams_per_device])
          CALL cuda_alloc_host(yf,[maxfft,cp_cuda_env%fft_n_devices_per_task*cp_cuda_env%fft_n_streams_per_device])
       ENDIF
#else
       CALL stopgm(procedureN,'shall not get to that point', __LINE__,__FILE__)
#endif
#if !defined(_USE_SCRATCHLIBRARY)
    ELSE
       ALLOCATE(xf(maxfft,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
       ALLOCATE(yf(maxfft,1),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
#endif
    ENDIF
#if defined(_HAS_CUDA) || !defined(_USE_SCRATCHLIBRARY)
    CALL zeroing(xf)!,SIZE(xf))
    CALL zeroing(yf)!,SIZE(yf))
#endif
    !
    DEALLOCATE(mg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(mz,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(my,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    !
    ! FFTNPOOLS
    fftpool = 0
    lmsqmax = nhrm
    lnzf    = ncpw%nhg
    lnzs    = ncpw%ngw
    ! ==--------------------------------------------------------------==
    ! ARRAY TO KEEP THE REAL SPACE WAVEFUNCTIONS
    rsactive = .FALSE.
    IF (cntl%krwfn.AND..NOT.batch_fft) THEN
       IF (tkpts%tkblock) CALL stopgm('FFTPRP',&
            'INCOMPATIBLE OPTIONS TKBLOCK AND KRWFN',&
            __LINE__,__FILE__)
       IF (parai%cp_nogrp.EQ.1) THEN
          lwdim = fpar%nnr1
          nstate = crge%n
          IF (tkpts%tkpnt) nstate = crge%n*nkpt%nkpnt
       ELSE
          !
          ! SuperDirtyFix (SDF): to avoid out of bound in the case we use more than 1 group,
          ! we allocate the memory for all the states.

          lwdim = fpar%nnr1
          ! LWDIM = KRX*KR2S*KR3S

          nstate = crge%n + 1
          ! NSTATE = (N+NOGRP-1)/NOGRP

          IF (tkpts%tkpnt) nstate = nstate*nkpt%nkpnt
       ENDIF
       !TK sometimes we have no real space plane...nnr1=0!
       CALL mp_max(lwdim,parai%allgrp)
       !
       ! SuperDirtyFix (SDF): allocate all the memory when CP_NOGRP.GT.1
       IF (cntr%memsize.LT.0.OR.parai%cp_nogrp.GT.1) THEN
          maxstates = nstate
          rstate=1._real_8
       ELSE
          rmem = 0.125_real_8*cntr%memsize*1.e6_real_8/REAL(lwdim,kind=real_8)
          maxstates = ((INT(rmem)+1)/2)*2
          maxstates = MIN(nstate,maxstates)
          rmem = -ABS(maxstates)
          CALL mp_max(rmem,parai%allgrp)
          maxstates = NINT(-rmem)
          rstate= REAL(maxstates,kind=real_8)/REAL(nstate,kind=real_8)
       ENDIF
       ldim  = (maxstates+1)/2 * lwdim
       IF (tkpts%tkpnt) ldim = ldim*2
       ALLOCATE(rswf(lwdim,ldim/lwdim),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       rmem = 16._real_8*ldim*1.e-6_real_8
       IF (paral%parent) THEN
          IF (rstate.GT.0.99999999_real_8) THEN
             IF (paral%io_parent)&
                  WRITE(6,'(A,T59,A)')&
                  ' FFTPRP| ORBITALS KEPT IN REAL SPACE ','    ALL'
          ELSE
             IF (paral%io_parent)&
                  WRITE(6,'(A,T50,F8.1,A)')&
                  ' FFTPRP| ORBITALS KEPT IN REAL SPACE ',&
                  100*rstate,' PERCENT'
          ENDIF
          IF (paral%io_parent)&
               WRITE(6,'(A,T51,F8.3,A)')&
               ' FFTPRP| WAVEFUNCTION TAKES IN REAL SPACE ',rmem,' MBYTES'
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    ! INITIALIZE BATCH FFT ALGORITHM
    IF(batch_fft)THEN
       nstate = crge%n
       CALL part_1d_get_blk_bounds(nstate,parai%cp_inter_me,parai%cp_nogrp,first,last)
       nstate = last - first +1
       IF(tkpts%tkpnt) nstate = nstate*nkpt%nkpnt
       maxstates = nstate
       fft_total=(maxstates+1)/2
       IF (tkpts%tkpnt) fft_total=maxstates+1
       lda=ngrm*nr1m
       !calculate fft_batchsize so that we get something less than a2a_msgsize to send per proc (sparse fft)
       !if autotuning is requested, we use this batchsize to set fft_max_numbatches
       a2a_msgsize=a2a_msgsize*1024/(parai%nproc*16)
       fft_batchsize=FLOOR(REAL(a2a_msgsize,KIND=real_8)/REAL(lda,KIND=real_8))
       IF(fft_batchsize.LE.0) fft_batchsize=1
       fft_residual=MOD(fft_total,fft_batchsize)
       fft_numbatches=(fft_total-fft_residual)/fft_batchsize

       IF(cntl%overlapp_comm_comp)THEN
          IF(cntl%krwfn)THEN
             fft_min_numbatches=2
          ELSE
             fft_min_numbatches=3
          END IF
       ELSE
          fft_min_numbatches=1
       END IF
       IF(fft_numbatches.LT.fft_min_numbatches)THEN
          fft_numbatches=fft_min_numbatches
          fft_batchsize=FLOOR(REAL(fft_total,KIND=real_8)/REAL(fft_numbatches,KIND=real_8))
          fft_residual=MOD(fft_total,fft_batchsize)
       END IF

       CALL mp_bcast(fft_min_numbatches,parai%source,parai%allgrp)
       CALL mp_bcast(fft_numbatches,parai%source,parai%allgrp)
       CALL mp_bcast(fft_residual,parai%source,parai%allgrp)
       CALL mp_bcast(fft_batchsize,parai%source,parai%allgrp)
       CALL mp_bcast(fft_total,parai%source,parai%allgrp)
       
       fft_tune_max_it=0

       IF(cntl%fft_tune_batchsize)THEN
          fft_tune_num_it=0
          fft_tune_max_it=(fft_batchsize*cnti%fft_tune_it_per_batch)+3
          call mp_max(fft_tune_max_it,parai%cp_grp)
          ALLOCATE(fft_time_total(fft_tune_max_it),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
          ALLOCATE(fft_batchsizes(fft_tune_max_it),STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
          fft_batchsizes=1
          !We start with the upper limit and decrease the batchsize in each iteration
          !the first FFT is always very slow so we skip it for the autotuning process.
          fft_batchsizes(1)=fft_batchsize
          fft_batchsizes(2)=fft_batchsize
          DO i=3,cnti%fft_tune_it_per_batch+2
             fft_batchsizes(i)=fft_batchsize
          END DO

          count=cnti%fft_tune_it_per_batch+2
          DO it=fft_tune_max_it-1,1,-1
             fft_batchsize=fft_batchsize-1
             IF(fft_batchsize.EQ.0)fft_batchsize=1
             fft_residual=MOD(fft_total,fft_batchsize)
             fft_numbatches=(fft_total-fft_residual)/fft_batchsize
             IF(fft_numbatches.LT.fft_min_numbatches)THEN
                fft_numbatches=fft_min_numbatches
                fft_batchsize=FLOOR(REAL(fft_total,KIND=real_8)/REAL(fft_numbatches,KIND=real_8))
                fft_residual=MOD(fft_total,fft_batchsize)
             END IF
             IF(fft_batchsize.NE.fft_batchsizes(count))THEN
                DO i=1,cnti%fft_tune_it_per_batch
                   count=count+1
                   fft_batchsizes(count)=fft_batchsize
                END DO
             end IF
          end DO
          fft_tune_max_it=count
          call mp_max(fft_tune_max_it,parai%cp_grp)
          fft_time_total=HUGE(0.0_real_8)
          fft_batchsize=fft_batchsizes(1)
          fft_residual=MOD(fft_total,fft_batchsize)
          fft_numbatches=(fft_total-fft_residual)/fft_batchsize          
       END IF
#if !defined(_USE_SCRATCHLIBRARY)
       IF(ALLOCATED(xf))THEN
          DEALLOCATE(xf,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       END IF
       IF(ALLOCATED(YF))THEN
          DEALLOCATE(yf,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       END IF
       ALLOCATE(xf(MAX(maxfft,fpar%nnr1*fft_batchsize),2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
       ALLOCATE(yf(MAX(maxfft,fpar%nnr1*fft_batchsize),2),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
#endif
       IF (paral%parent) THEN
          WRITE(6,*)&
               ' FFTPRP| BATCH FFT NUMBER ALLTOALL CALLS ',fft_numbatches
          WRITE(6,*)&
               ' FFTPRP| BATCH FFT NUMBER OF STATES PER CALL ',fft_batchsize
       ENDIF
    ENDIF
    ! ==--------------------------------------------------------------==
    IF (paral%parent) CALL prmem('    FFTPRP')
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fftprp_default_init

  ! ==================================================================
  SUBROUTINE autotune_fftbatchsize()
    ! ==--------------------------------------------------------------==
    ! == tunes the number of states per fft batch                     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: ind(1), ierr
    CHARACTER(*), PARAMETER                  :: procedureN = 'autotune_fftbatchsize'
    REAL(real_8)                             :: last, now
    INTEGER                                  :: i
    ! ==--------------------------------------------------------------==
    IF(fft_tune_max_it.EQ.0)RETURN

    IF(fft_tune_num_it.EQ.fft_tune_max_it)THEN
       cntl%fft_tune_batchsize=.FALSE.
    ELSE
       fft_tune_num_it=fft_tune_num_it+1
       fft_batchsize=fft_batchsizes(fft_tune_num_it)
       fft_residual=MOD(fft_total,fft_batchsize)
       fft_numbatches=(fft_total-fft_residual)/fft_batchsize
       IF(fft_tune_num_it.EQ.1)THEN
          RETURN
       ELSEIF(fft_tune_num_it.EQ.2)THEN
          CALL mp_max(fft_time_total(fft_tune_num_it-1),parai%cp_grp)
          RETURN
       END IF
       CALL mp_max(fft_time_total(fft_tune_num_it-1),parai%cp_grp)
    END IF
       
    IF(.NOT.cntl%fft_tune_batchsize)THEN
       DO i=3,fft_tune_max_it,cnti%fft_tune_it_per_batch
          fft_time_total(i+1:i+cnti%fft_tune_it_per_batch-1)=&
               SUM(fft_time_total(i+1:i+cnti%fft_tune_it_per_batch-1))/&
               (REAL(cnti%fft_tune_it_per_batch-1,KIND=real_8))
       END DO
       ind=MINLOC(fft_time_total)
       fft_batchsize=fft_batchsizes(ind(1))
       if(ibench(4).eq.1)then
          if(paral%io_parent)write(6,*) fft_time_total
          if(paral%io_parent)write(6,*) fft_batchsizes
          if(paral%io_parent)write(6,*) ind,fft_time_total(ind),fft_batchsizes(ind)
       end if
       fft_residual=MOD(fft_total,fft_batchsize)
       fft_numbatches=(fft_total-fft_residual)/fft_batchsize
       IF(fft_numbatches.LT.fft_min_numbatches)THEN
          fft_numbatches=fft_min_numbatches
!          fft_batchsize=FLOOR(REAL(fft_total,KIND=real_8)/REAL(fft_numbatches,KIND=real_8))
          fft_batchsize=fft_total/fft_numbatches
          fft_residual=MOD(fft_total,fft_batchsize)
       END IF
    END IF

    CALL mp_bcast(fft_numbatches,parai%source,parai%allgrp)
    CALL mp_bcast(fft_residual,parai%source,parai%allgrp)
    CALL mp_bcast(fft_batchsize,parai%source,parai%allgrp)

#if !defined(_USE_SCRATCHLIBRARY)
    IF(ALLOCATED(xf))THEN
       DEALLOCATE(xf,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    END IF
    IF(ALLOCATED(YF))THEN
       DEALLOCATE(yf,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    END IF
    ALLOCATE(xf(MAX(maxfft,fpar%nnr1*fft_batchsize),2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
    ALLOCATE(yf(MAX(maxfft,fpar%nnr1*fft_batchsize),2),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', __LINE__,__FILE__)
#endif
    IF (.NOT.cntl%fft_tune_batchsize.AND.paral%parent) THEN
       WRITE(6,*)&
            ' autotune_fftbatchsize| BATCH FFT NUMBER ALLTOALL CALLS ',fft_numbatches
       WRITE(6,*)&
            ' autotune_fftbatchsize| BATCH FFT NUMBER OF STATES PER CALL ',fft_batchsize
       WRITE(6,*)&
            ' autotune_fftbatchsize| A2A MESSAGE SIZE ',fft_batchsize*ngrm*nr1m*parai%nproc*16/1024
   ENDIF

    IF(.NOT.cntl%fft_tune_batchsize)THEN
       IF(ALLOCATED(fft_time_total))THEN
          DEALLOCATE(fft_time_total,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)
          DEALLOCATE(fft_batchsizes,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', __LINE__,__FILE__)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE autotune_fftbatchsize
  ! ==================================================================

  ! ==================================================================
  SUBROUTINE fftprp_default_finalize
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'fftprp_default_finalize'
    INTEGER :: ierr

    DEALLOCATE(nzhs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(indzs,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(msp,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF( cp_cuda_env%use_fft ) THEN
       IF(cuda_register_memory) THEN
          CALL cuda_host_unregister(xf)
          CALL cuda_host_unregister(yf)
          DEALLOCATE(xf,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
          DEALLOCATE(yf,STAT=ierr)
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ELSE
#if defined(_HAS_CUDA)
          CALL cuda_dealloc_host(xf)
          CALL cuda_dealloc_host(yf)
#else
       CALL stopgm(procedureN,'shall not get to that point', __LINE__,__FILE__)
#endif
       ENDIF
#if !defined(_USE_SCRATCHLIBRARY)
    ELSE
       DEALLOCATE(xf,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
       DEALLOCATE(yf,STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
#endif
    ENDIF


  END SUBROUTINE fftprp_default_finalize


  ! ==================================================================
  SUBROUTINE grpgs
    ! ==--------------------------------------------------------------==
    CHARACTER(*), PARAMETER                  :: procedureN = 'grpgs'

    INTEGER                                  :: ierr, ii, io, ip, ipro, mxrpt
    INTEGER, POINTER                         :: ms_1d(:)

! ==--------------------------------------------------------------==

    ipro=0
    maxrpt=0
    DO ip=1,group%npgrp
       mxrpt=0
       DO io=1,group%nogrp
          mxrpt=mxrpt+parap%sparm(9,ipro)
          ipro=ipro+1
       ENDDO
       mxy(ip)=mxrpt
       IF (maxrpt.LT.mxrpt) maxrpt=mxrpt
    ENDDO
    DEALLOCATE(ms,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(ms((maxrpt+1),group%npgrp),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL reshape_inplace(ms, (/(maxrpt+1)*group%npgrp/), ms_1d)

    ipro=0
    DO ip=1,group%npgrp
       mxrpt=1
       DO io=1,group%nogrp
          ii=(ip-1)*maxrpt+mxrpt
          CALL icopy(parap%sparm(9,ipro),msp(:,2,ipro+1),1,ms_1d(ii:),1)
          mxrpt=mxrpt+parap%sparm(9,ipro)
          ipro=ipro+1
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
  END SUBROUTINE grpgs
  ! ==================================================================

END MODULE fftprp_utils
