#if defined(__FFT_HAS_LOW_LEVEL_TIMERS)
#define HAS_LOW_LEVEL_TIMERS .TRUE.
#else
#define HAS_LOW_LEVEL_TIMERS .FALSE.
#endif

#include "cpmd_global.h"
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
  include "fftw/offload/fftw3_omp_offload.f90"
#endif

MODULE mltfft_utils

  USE cublas_types,                    ONLY: cublas_handle_t
  USE cublas_utils,                    ONLY: cublas_zdscal
  USE cuda_types,                      ONLY: cuda_memory_t,&
                                             cuda_stream_t
  USE cufft_interfaces,                ONLY: cufft_forward,&
                                             cufft_inverse
  USE cufft_types,                     ONLY: cufft_plan_t
  USE cufft_utils,                     ONLY: cufft_execz2z
  USE cuuser_utils,                    ONLY: CuUser_setblock2zero,&
                                             CuUser_setblock2zero_scale
  USE error_handling,                  ONLY: stopgm
  USE gfft_utils,                      ONLY: ctrig
  USE gpu
  USE kinds,                           ONLY: int_4,&
                                             int_8,&
                                             real_8
  USE utils,                           ONLY: numcpus
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE system,                          ONLY: cntl
  USE fft, ONLY: fft_batchsize
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_NULL_CHAR,&
       C_INT, C_CHAR,&
       C_PTR,&
       C_LONG_DOUBLE_COMPLEX
  USE timer,                           ONLY: tihalt,&
       tiset

#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: mltfft_default
  PUBLIC :: mltfft_essl
  PUBLIC :: mltfft_hp
  PUBLIC :: mltfft_fftw
  PUBLIC :: mltfft_cuda

CONTAINS

  ! ==================================================================
  SUBROUTINE mltfft_default(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,isign,scale)
    CHARACTER(len=1)                         :: transa, transb
    INTEGER                                  :: ldax, lday
    COMPLEX(real_8)                          :: a(ldax,lday)
    INTEGER                                  :: ldbx, ldby
    COMPLEX(real_8)                          :: b(ldbx,ldby)
    INTEGER                                  :: n, m, isign
    REAL(real_8)                             :: scale

    CHARACTER(*), PARAMETER                  :: procedureN = 'mltfft_default'

    INTEGER :: after(20), before(20), i, ic, icpu, ierr, il, inzee, ir, isig, &
      itr, j, length, lot, mh, ml, mm, ncpus, nfft, now(20)
    LOGICAL                                  :: tscal
    REAL(real_8)                             :: trig(2,1024)

#if defined(__NEC)
    INTEGER, PARAMETER :: ncache=1024*1024
#elif defined(__PRIMEHPC)
    INTEGER, PARAMETER :: ncache=1024*32
#elif defined(__VECTOR)
#ifdef __SR8000
    INTEGER, PARAMETER :: ncache=1024*8   ! cmb - for SR8K
#else
    INTEGER, PARAMETER :: ncache=1024*128
#endif
#elif defined(__SGI)
    INTEGER, PARAMETER :: ncache=1024*4
#elif defined(__IBM)
#ifdef __SR11KIBM
    INTEGER, PARAMETER :: ncache=1024*32  ! cmb - for SR11K
#else
    INTEGER, PARAMETER :: ncache=1024*10
#endif
#elif defined(__OSX) || defined(__OSX_IFC)
    INTEGER, PARAMETER :: ncache=1024*10
#elif defined(__i386)
    INTEGER, PARAMETER :: ncache=1024*10
#elif defined(__x86_64)
#ifdef __PACS
    INTEGER, PARAMETER :: ncache=1024*8   ! cmb - for PACS-CS
#else
    INTEGER, PARAMETER :: ncache=1024*6
#endif
#elif defined(__HP)
    INTEGER, PARAMETER :: ncache=1024*64
#elif defined(__ALTIX)
    INTEGER, PARAMETER :: ncache=1024*16  ! Altix3700/3900
#elif defined(_FFT_NCACHE)
    INTEGER, PARAMETER :: ncache=_fft_ncache  ! for convenient overrides in makefile
#else
    INTEGER, PARAMETER :: ncache=1024*2
#endif
#if defined(__SR8000) || defined(__SR11000)
    REAL(real_8) :: z(2,ncache/4+1,2,maxcpu/8)
#else
    REAL(real_8), ALLOCATABLE, SAVE :: z(:,:,:,:)
#endif
    INTEGER, SAVE :: ifirst=0

    ! ==-------------------------------------------------------------==
    CALL numcpus(ncpus)
    IF (ifirst.EQ.0) THEN
       ifirst=1
       length=4*(ncache/4+1)*ncpus
       ALLOCATE(z(2,ncache/4+1,2,ncpus),STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ENDIF
    isig=-isign
    tscal=(ABS(scale-1._real_8).GT.1.e-12_real_8)
    CALL ctrig(n,trig,after,before,now,isig,ic)
    lot=ncache/(4*n)
    lot=lot-MOD(lot+1,2)
    lot=MAX(1,lot)
    IF (4*n*lot.GT.ncache) CALL stopgm('MLFFT_DEFAULT',' NCACHE ',&
         __LINE__,__FILE__)
    mm=MAX(m/ncpus,1)
#if defined(__SR11000) || defined(__BG)
    !$omp parallel do private(nfft,icpu,ml,mh,itr,inzee,i,j,il,ir) &
    !$omp  shared(trig,now,after,before,tscal,scale, &
    !$omp         ldax,lday,ldbx,ldby,isig,n,a,b,lot,mm,m,ncpus, &
    !$omp         transa,transb,ic)
#else
    !$omp parallel do private(nfft,icpu,ml,mh,itr,inzee,i) &
    !$omp  shared(trig,now,after,before,tscal,scale, &
    !$omp         ldax,lday,ldbx,ldby,isig,n,a,b,lot,mm,m,ncpus, &
    !$omp         transa,transb,ic)
#endif
    DO icpu=1,ncpus
       ml=(icpu-1)*mm+1
       mh=MIN(icpu*mm,m)
       IF (icpu.EQ.ncpus) mh=m
       DO itr=ml,mh,lot
          nfft=MIN(mh-itr+1,lot)
          IF (transa.EQ.'N'.OR.transa.EQ.'n') THEN
             CALL fftpre(nfft,nfft,ldax,lot,n,a(1,itr),z(1,1,1,icpu),&
                  trig,now(1),after(1),before(1),isig)
          ELSE
             CALL fftstp(ldax,nfft,n,lot,n,a(itr,1),z(1,1,1,icpu),&
                  trig,now(1),after(1),before(1),isig)
          ENDIF
          ! IF(TSCAL) CALL DSCAL(2*LOT*N,SCALE,Z(1,1,1,ICPU),1)
          IF (tscal) THEN
             IF (lot.EQ.nfft) THEN
#if defined(__SR8000) || defined(__SR11000) || defined(__BG)
                DO j=1,lot*n
                   z(1,j,1,icpu)=scale*z(1,j,1,icpu)
                   z(2,j,1,icpu)=scale*z(2,j,1,icpu)
                ENDDO
#else
                CALL dscal(2*lot*n,scale,z(1,1,1,icpu),1)
#endif
             ELSE
#if defined(__SR8000) || defined(__SR11000) || defined(__BG)
                DO i=1,n
                   DO j=1,nfft
                      il=lot*(i-1)+j
                      ir=il
                      z(1,il,1,icpu)=scale*z(1,ir,1,icpu)
                      z(2,il,1,icpu)=scale*z(2,ir,1,icpu)
                   ENDDO
                ENDDO
#else
                DO i=1,n
                   CALL dscal(2*nfft,scale,z(1,lot*(i-1)+1,1,icpu),1)
                ENDDO
#endif
             ENDIF
          ENDIF
          IF (ic.EQ.1) THEN
             IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
#if defined(__EREG)
                CALL cmplx_transpose(nfft,n,z(1,1,1,icpu),lot,&
                     b(itr,1),ldbx)
#elif defined(__SR8000) || defined(__SR11000) || defined(__BG)
                DO i=1,nfft
                   DO j=1,n
                      il=ldbx*(i-1)+j
                      ir= lot*(j-1)+i
                      b(il,itr)=CMPLX(z(1,ir,1,icpu),z(2,i,ir,icpu),kind=real_8)
                   ENDDO
                ENDDO
#else
                CALL zgetmo(z(1,1,1,icpu),lot,nfft,n,b(1,itr),ldbx)
#endif
             ELSE
#if defined(__SR8000) || defined(__SR11000) || defined(__BG)
                DO i=1,nfft
                   DO j=1,n
                      il=ldbx*(i-1)+j
                      ir= lot*(i-1)+j
                      b(itr,il)=CMPLX(z(1,ir,1,icpu),z(2,i,ir,icpu),kind=real_8)
                   ENDDO
                ENDDO
#else
                CALL matmov(nfft,n,z(1,1,1,icpu),lot,b(itr,1),ldbx)
#endif
             ENDIF
          ELSE
             inzee=1
             DO i=2,ic-1
                CALL fftstp(lot,nfft,n,lot,n,z(1,1,inzee,icpu),&
                     z(1,1,3-inzee,icpu),trig,now(i),after(i),&
                     before(i),isig)
                inzee=3-inzee
             ENDDO
             IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
                CALL fftrot(lot,nfft,n,nfft,ldbx,z(1,1,inzee,icpu),&
                     b(1,itr),trig,now(ic),after(ic),before(ic),isig)
             ELSE
                CALL fftstp(lot,nfft,n,ldbx,n,z(1,1,inzee,icpu),&
                     b(itr,1),trig,now(ic),after(ic),before(ic),isig)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    ! ==--------------------------------------------------------------==
    IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
       !$omp parallel do private (I,J)
       DO j=m+1,ldby
          DO i=1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
       !$omp parallel do private (I,J)
       DO j=1,m
          DO i=n+1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ELSE
       !$omp parallel do private (I,J)
       DO j=n+1,ldby
          DO i=1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
       !$omp parallel do private (I,J)
       DO j=1,n
          DO i=m+1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mltfft_default
  ! ==================================================================
  SUBROUTINE mltfft_essl(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,&
       isign,scale)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: transa, transb
    INTEGER                                  :: ldax, lday
    COMPLEX(real_8)                          :: a(ldax,lday)
    INTEGER                                  :: ldbx, ldby
    COMPLEX(real_8)                          :: b(ldbx,ldby)
    INTEGER                                  :: n, m, isign
    REAL(real_8)                             :: scale

    INTEGER, PARAMETER                       :: naux1 = 100000, naux2 = 100000

    INTEGER                                  :: i, j
    REAL(real_8)                             :: aux1(naux1), aux2(naux2)

#if defined(__HAS_FFT_ESSL)
    IF (transa.EQ.'N'.OR.transa.EQ.'n') THEN
       IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
          CALL dcft(1,a,1,ldax,b,1,ldbx,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
          CALL dcft(0,a,1,ldax,b,1,ldbx,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
       ELSE
          CALL dcft(1,a,1,ldax,b,ldbx,1,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
          CALL dcft(0,a,1,ldax,b,ldbx,1,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
       ENDIF
    ELSE
       IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
          CALL dcft(1,a,ldax,1,b,1,ldbx,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
          CALL dcft(0,a,ldax,1,b,1,ldbx,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
       ELSE
          CALL dcft(1,a,ldax,1,b,ldbx,1,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
          CALL dcft(0,a,ldax,1,b,ldbx,1,n,m,isign,scale,&
               aux1,naux1,aux2,naux2)
       ENDIF
    ENDIF
#endif
    ! ==--------------------------------------------------------------==
    IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
       DO j=m+1,ldby
          DO i=1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
       DO j=1,m
          DO i=n+1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ELSE
       DO j=n+1,ldby
          DO i=1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
       DO j=1,n
          DO i=m+1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mltfft_essl
  ! ==================================================================
  SUBROUTINE mltfft_hp(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,&
       isign,scale)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: transa, transb
    INTEGER                                  :: ldax, lday
    COMPLEX(real_8)                          :: a(ldax,lday)
    INTEGER                                  :: ldbx, ldby
    COMPLEX(real_8)                          :: b(ldbx,ldby)
    INTEGER                                  :: n, m, isign
    REAL(real_8)                             :: scale

    INTEGER(int_4)                           :: i, ier, iopt, j
    REAL(real_8)                             :: ascale

! ==--------------------------------------------------------------==

#if defined(__HAS_FFT_HP)
    IF (transa.EQ.'N'.OR.transa.EQ.'n') THEN
       IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
          ! CALL DCFT(1,A,1,LDAX,B,1,LDBX,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! CALL DCFT(0,A,1,LDAX,B,1,LDBX,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          CALL zgecpy ('N', n, m, a, ldax, b, ldbx)
          IF (isign .GE. 0) THEN
             iopt = -1
          ELSE
             iopt = 1
          ENDIF
          CALL zffts (b, n, 1, m, ldbx, iopt, ier)
          IF (ier .NE. 0) THEN
             CALL stopgm('MLFFT_HP','ERROR IN CALL TO ZFFTS',&
                  __LINE__,__FILE__)
          ENDIF
          IF (isign .GE. 0) THEN
             ascale = scale*REAL(n)
          ELSE
             ascale = scale
          ENDIF
          IF (ascale .NE. 1.0_real_8) THEN
             DO i = 1,m
                CALL zdscal(n, ascale, b(1,i), 1)
             ENDDO
          ENDIF
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ELSE
          ! CALL DCFT(1,A,1,LDAX,B,LDBX,1,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! CALL DCFT(0,A,1,LDAX,B,LDBX,1,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          CALL zgecpy ('T', n, m, a, ldax, b, ldbx)
          IF (isign .GE. 0) THEN
             iopt = -1
          ELSE
             iopt = 1
          ENDIF
          CALL zffts (b, n, ldbx, m, 1, iopt, ier)
          IF (ier .NE. 0) THEN
             CALL stopgm('MLFFT_HP','ERROR IN CALL TO ZFFTS',&
                  __LINE__,__FILE__)
          ENDIF
          IF (isign .GE. 0) THEN
             ascale = scale*REAL(n)
          ELSE
             ascale = scale
          ENDIF
          IF (ascale .NE. 1.0_real_8) THEN
             DO i = 1,n
                CALL zdscal(m, ascale, b(1,i), 1)
             ENDDO
          ENDIF
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ENDIF
    ELSE
       IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
          ! CALL DCFT(1,A,LDAX,1,B,1,LDBX,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! CALL DCFT(0,A,LDAX,1,B,1,LDBX,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          CALL zgecpy ('T', m, n, a, ldax, b, ldbx)
          IF (isign .GE. 0) THEN
             iopt = -1
          ELSE
             iopt = 1
          ENDIF
          CALL zffts (b, n, 1, m, ldbx, iopt, ier)
          IF (ier .NE. 0) THEN
             CALL stopgm('MLFFT_HP','ERROR IN CALL TO ZFFTS',&
                  __LINE__,__FILE__)
          ENDIF
          IF (isign .GE. 0) THEN
             ascale = scale*REAL(n)
          ELSE
             ascale = scale
          ENDIF
          IF (ascale .NE. 1.0_real_8) THEN
             DO i = 1,m
                CALL zdscal(n, ascale, b(1,i), 1)
             ENDDO
          ENDIF
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ELSE
          ! CALL DCFT(1,A,LDAX,1,B,LDBX,1,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! CALL DCFT(0,A,LDAX,1,B,LDBX,1,N,M,ISIGN,SCALE,
          ! *            AUX1,NAUX1,AUX2,NAUX2)
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          CALL zgecpy ('N', m, n, a, ldax, b, ldbx)
          IF (isign .GE. 0) THEN
             iopt = -1
          ELSE
             iopt = 1
          ENDIF
          CALL zffts (b, n, ldbx, m, 1, iopt, ier)
          IF (ier .NE. 0) THEN
             CALL stopgm('MLFFT_HP','ERROR IN CALL TO ZFFTS',&
                  __LINE__,__FILE__)
          ENDIF
          IF (isign .GE. 0) THEN
             ascale = scale*REAL(n)
          ELSE
             ascale = scale
          ENDIF
          IF (ascale .NE. 1.0_real_8) THEN
             DO i = 1,n
                CALL zdscal(m, ascale, b(1,i), 1)
             ENDDO
          ENDIF
          ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ENDIF
    ENDIF
#endif
    ! ==--------------------------------------------------------------==
    IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
       DO j=m+1,ldby
          DO i=1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
       DO j=1,m
          DO i=n+1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ELSE
       DO j=n+1,ldby
          DO i=1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
       DO j=1,n
          DO i=m+1,ldbx
             b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
          ENDDO
       ENDDO
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE mltfft_hp
  ! ==================================================================
  SUBROUTINE mltfft_fftw(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,&
       isign,scale,use_wisdom)
    ! ==--------------------------------------------------------------==
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
    USE FFTW3_OMP_OFFLOAD
#endif
    CHARACTER(len=1)                         :: transa, transb
    INTEGER                                  :: ldax
    COMPLEX(real_8),TARGET                          :: a(ldax,*)
    INTEGER                                  :: lday, ldbx
    COMPLEX(real_8),TARGET                          :: b(ldbx,*)
    INTEGER                                  :: ldby, n, m, isign
    REAL(real_8)                             :: scale
    LOGICAL                                  :: use_wisdom

    INTEGER, PARAMETER                       :: fftw_b = 1, fftw_es = 64, &
                                                fftw_f = -1 , fftw_me = 0, &
                                                fftw_pa = 32

    INTEGER                                  :: fftw_dir, fftw_flags, i, j, &
                                                isub1, isub2, isub3, isub4
    LOGICAL                                  :: tscal
    INTEGER(int_8)                           :: plan
    CHARACTER(*), PARAMETER                  :: procedureN1='dfftw_execute'
    CHARACTER(*), PARAMETER                  :: procedureN2='fftw3 zero'
    complex(real_8), pointer, contiguous    :: a_ptr(:,:), b_ptr(:,:)
#if defined(_HAS_FFT_FFTW3)
    tscal=(ABS(scale-1._real_8).GT.1.e-12_real_8)
    IF (isign.EQ.1) THEN
       fftw_dir=fftw_f
    ELSE
       fftw_dir=fftw_b
    ENDIF
    IF (use_wisdom) THEN
       fftw_flags=fftw_me
    ELSE
       fftw_flags=fftw_es
    ENDIF

    
    CALL get_fftw_plan(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,&
         isign,plan)
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN1//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN1,isub2)
       END IF
    END IF

    a_ptr=> a(1:ldax,1:lday)
    b_ptr=> b(1:ldbx,1:ldby)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target enter data map(alloc:a(1:ldax,1:lday),b(1:ldbx,1:ldby))
    !$omp target update to(a(1:ldax,1:lday)) if(update_first_to_gpu)
#endif
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)    
    !$omp dispatch
#endif
    CALL dfftw_execute_dft(plan,a_ptr,b_ptr)
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN1//'_tuning',isub1)
          CALL tiset(procedureN2//'_tuning',isub3)
       ELSE
          CALL tihalt(procedureN1,isub2)
          CALL tiset(procedureN2,isub4)
       END IF
    END IF
    IF (tscal) THEN
       IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target teams distribute parallel do simd collapse(2)
#else
          !$omp parallel do private(I,J)
#endif
          DO i = 1,m
             DO j = 1,n
                b(j,i)=scale*b(j,i)
             ENDDO
          ENDDO
       ELSE
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target teams distribute parallel do simd collapse(2)
#else
          !$omp parallel do private(I,J)
#endif
          DO i = 1,n
             DO j = 1,m
                b(j,i)=scale*b(j,i)
             ENDDO
          ENDDO
       END IF
    ENDIF
#endif
    IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
       IF(m+1.le.ldby)then
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target teams distribute parallel do simd collapse(2)
#else
          !$omp parallel do private (I,J)
#endif
          DO j=m+1,ldby
             DO i=1,ldbx
                b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
             ENDDO
          ENDDO
       END IF
       IF(n+1.le.ldbx)then
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target teams distribute parallel do simd collapse(2)
#else
          !$omp parallel do private (I,J)
#endif
          DO j=1,m
             DO i=n+1,ldbx
                b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
             ENDDO
          ENDDO
       END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
!       !$omp end target teams
#endif
    ELSE
       IF(n+1.le.ldby)then
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target teams distribute parallel do simd collapse(2)
#else
          !$omp parallel do private (I,J)
#endif
          DO j=n+1,ldby
             DO i=1,ldbx
                b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
             ENDDO
          ENDDO
       END IF
       IF(m+1.le.ldbx)THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
          !$omp target teams distribute parallel do simd collapse(2)
#else
          !$omp parallel do private (I,J)
#endif
          DO j=1,n
             DO i=m+1,ldbx
                b(i,j)=CMPLX(0.0_real_8,0.0_real_8,kind=real_8)
             ENDDO
          ENDDO
       END IF
#if defined(_HAS_OMP_TARGET_OFFLOAD)
!       !$omp end target teams
#endif       
    ENDIF
#if defined(_HAS_OMP_TARGET_OFFLOAD)    
    !$omp target update from(b(1:ldbx,1:ldby)) if(update_result_to_host)
    !$omp target exit data map(release:a(1:ldax,1:lday),b(1:ldbx,1:ldby))
#endif
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN2//'_tuning',isub3)
       ELSE
          CALL tihalt(procedureN2,isub4)
       END IF
    END IF
  END SUBROUTINE mltfft_fftw
  ! ==================================================================
#if defined(_HAS_FFT_FFTW3)
  SUBROUTINE get_fftw_plan(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,&
       isign,plan)
    CHARACTER(len=1), INTENT( IN )           :: transa, transb
    INTEGER, INTENT( IN )                    :: ldax, lday, ldbx, ldby
    COMPLEX(real_8), INTENT( INOUT )         :: a(ldax,*)
    COMPLEX(real_8), INTENT( INOUT )         :: b(ldbx,*)
    INTEGER, INTENT( IN )                    :: n, m, isign
    INTEGER(int_8), INTENT( OUT)             :: plan

    INTEGER, PARAMETER                       :: fftw_b = 1, fftw_f = -1

    INTEGER, ALLOCATABLE, SAVE               :: params(:,:)
    INTEGER(int_8), ALLOCATABLE, SAVE        :: plans(:)
    INTEGER                                  :: params_req(7), fftw_dir, ierr, &
                                                num_plans, i, isub1, isub2
    INTEGER(int_8)                           :: il_asave(2)
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: asave(:,:)
    LOGICAL                                  :: first
    LOGICAL, SAVE                            :: batch_tuning = .FALSE.
    CHARACTER(*), PARAMETER                  :: procedureN = 'get_fftw_plan'
    INTEGER, SAVE                            :: fft_batchsize_plans=0

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF

    IF( cntl%fft_tune_batchsize .AND. .NOT. batch_tuning )THEN
       batch_tuning = .TRUE.
    END IF

    IF( .NOT. cntl%fft_tune_batchsize .AND. batch_tuning )THEN
       batch_tuning = .FALSE.
       CALL clean_fftw_plans(plans, params)
    END IF
    IF(fft_batchsize.NE.fft_batchsize_plans.AND.batch_tuning)THEN
       FFT_batchsize_plans=fft_batchsize
       CALL clean_fftw_plans(plans, params)
    END IF

    IF (isign.EQ.1) THEN
       fftw_dir=fftw_f
    ELSE
       fftw_dir=fftw_b
    ENDIF

    params_req(1)=n
    params_req(2)=m
    CALL get_stride( transa, transb, ldax, ldbx, params_req(3), params_req(4), params_req(5), &
         params_req(6) )
    params_req(7)=fftw_dir
    first = .FALSE.

    IF( .NOT. ALLOCATED( plans ))THEN
       ALLOCATE( params( 7, 1 ), STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       ALLOCATE( plans( 1 ), STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       params = 0
       first = .TRUE.
    END IF
    num_plans = SIZE( plans, 1 )

    DO i = 1, num_plans
       IF(  params_req(1) .EQ. params(1,i) .AND. params_req(2) .EQ. params(2,i) .AND. &
            params_req(3) .EQ. params(3,i) .AND. params_req(4) .EQ. params(4,i) .AND. &
            params_req(5) .EQ. params(5,i) .AND. params_req(6) .EQ. params(6,i) .AND. &
            params_req(7) .EQ. params(7,i) )EXIT
    END DO
    IF( i .LE. num_plans ) THEN
       plan = plans(i)
    ELSE
       IF( .NOT. first )THEN
          num_plans = num_plans + 1
          CALL grow_plan_index(plans,params)
       END IF
       params(:,num_plans)=params_req
       
       il_asave( 1 ) = INT( ldax, KIND = int_8)
       il_asave( 2 ) = INT( lday, KIND = int_8)
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_asave,asave,procedureN//'_asave',ierr)
#else
       ALLOCATE(asave(il_asave(1),il_asave(2)),STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
            __LINE__,__FILE__)
       CALL cpmd_dcopy(ldax*lday*2, a, 1, asave,  1 )
       CALL create_fftw_plan(plan,params(:,num_plans),a,b,ldax,lday,ldbx,ldby)
       CALL cpmd_dcopy(ldax*lday*2, asave, 1, a,  1 )
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_asave,asave,procedureN//'_asave',ierr)
#else
       DEALLOCATE(asave,STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
            __LINE__,__FILE__)

       plans(num_plans) = plan
    END IF

  END SUBROUTINE get_fftw_plan
  ! ==================================================================  
  SUBROUTINE get_stride( transa, transb, ldax, ldbx , istride, idist, ostride, odist )
    CHARACTER(len=1), INTENT( IN )           :: transa, transb
    INTEGER, INTENT( IN )                    :: ldax, ldbx
    INTEGER, INTENT( OUT )                   :: istride, idist, ostride, odist

    IF (transa.EQ.'N'.OR.transa.EQ.'n') THEN
       istride = 1
       idist = ldax
       IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
          ostride = 1
          odist = ldbx
       ELSE
          ostride = ldbx
          odist = 1
       ENDIF
    ELSE
       istride = ldax
       idist = 1
       IF (transb.EQ.'N'.OR.transb.EQ.'n') THEN
          ostride = 1
          odist = ldbx
       ELSE
          ostride = ldbx
          odist = 1
       ENDIF
    ENDIF

  END SUBROUTINE get_stride
  ! ==================================================================
  SUBROUTINE grow_plan_index(plans,params)
    INTEGER(int_8),ALLOCATABLE,INTENT(INOUT) :: plans(:)
    INTEGER,ALLOCATABLE,INTENT(INOUT)        :: params(:,:)

    INTEGER(int_8),ALLOCATABLE               :: swap_plans(:)
    INTEGER,ALLOCATABLE                      :: swap_params(:,:)
    INTEGER                                  :: num_plans
    INTEGER                                  :: i, ierr
    CHARACTER(*),PARAMETER                   :: procedureN = 'grow_plan_index'

    num_plans = SIZE( plans, 1 )
    ALLOCATE( swap_params( 7, num_plans + 1 ), STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    ALLOCATE( swap_plans( num_plans + 1 ), STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)

    DO i = 1, num_plans
       swap_params(:,i) = params(:,i)
       swap_plans(i) = plans(i)
    END DO

    CALL MOVE_ALLOC( swap_params, params)
    CALL MOVE_ALLOC( swap_plans, plans)
  END SUBROUTINE grow_plan_index
  ! ==================================================================
  SUBROUTINE clean_fftw_plans(plans,params)
    INTEGER(int_8),ALLOCATABLE,INTENT(INOUT) :: plans(:)
    INTEGER,ALLOCATABLE,INTENT(INOUT)        :: params(:,:)
    INTEGER                                  :: i, ierr
    CHARACTER(*),PARAMETER                   :: procedureN = 'clean_plans'

    IF(ALLOCATED( params ))THEN
       DEALLOCATE( params, STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    END IF
    IF(ALLOCATED(plans))THEN
       DO i = 1, SIZE( plans, 1 )
          CALL dfftw_destroy_plan(plans(i))
       END DO
       DEALLOCATE( plans, STAT=ierr)
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    END IF
  END SUBROUTINE clean_fftw_plans
  ! ==================================================================
  SUBROUTINE create_fftw_plan(plan,params,a,b,ldax,lday,ldbx,ldby)
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
    USE FFTW3_OMP_OFFLOAD
#endif
    INTEGER(int_8),INTENT(OUT)      :: plan
    INTEGER,INTENT(IN)              :: params(7),ldax,lday,ldbx,ldby
    COMPLEX(real_8),INTENT(INOUT), TARGET &
                                    :: a(ldax,*),b(ldbx,*)

    INTEGER                         :: n, m, istride, idist, ostride, odist, fftw_dir, rank, size(1), &
                                       howmany, inembed(1), onembed(1)
    INTEGER, PARAMETER              :: fftw_exhaustive = 8, fftw_measure=0, fftw_estimate=64
    COMPLEX(real_8), POINTER &
         __CONTIGUOUS                :: a_ptr(:,:), b_ptr(:,:)
    n        = params( 1 )
    m        = params( 2 )
    istride  = params( 3 )
    idist    = params( 4 )
    ostride  = params( 5 )
    odist    = params( 6 )
    fftw_dir = params( 7 )
    rank = 1
    size = n
    howmany = m   
    inembed = n
    onembed = n

    a_ptr=> a(1:ldax,1:lday)
    b_ptr=> b(1:ldbx,1:ldby)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target enter data map(alloc:a(1:ldax,1:lday),b(1:ldbx,1:ldby))
#endif
#if defined(_HAS_OMP_TARGET_OFFLOAD_INTEL)
    !$omp dispatch
#endif
    CALL dfftw_plan_many_dft(plan,rank,size,howmany,a_ptr,inembed,istride,idist,b_ptr,onembed,ostride,odist,&
         fftw_dir,fftw_estimate)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target exit data map(release:a(1:ldax,1:lday),b(1:ldbx,1:ldby))
#endif

  END SUBROUTINE create_fftw_plan
#endif
  ! ==================================================================
  SUBROUTINE mltfft_cuda(transa,transb,a,ldax,lday,b,ldbx,ldby,n,m,isign,scale, plan, blas_handle, stream)
    !     ==--------------------------------------------------------------==
    CHARACTER(len=1)                         :: transa, transb
    TYPE(cuda_memory_t)                      :: a
    INTEGER                                  :: ldax, lday
    TYPE(cuda_memory_t)                      :: b
    INTEGER                                  :: ldbx, ldby, n, m, isign
    REAL(real_8)                             :: scale
    TYPE(cufft_plan_t), INTENT(IN)           :: plan
    TYPE(cublas_handle_t), INTENT(IN)        :: blas_handle
    TYPE(cuda_stream_t), INTENT(IN)          :: stream

    CHARACTER(*), PARAMETER                  :: procedureN = 'mltfft_cuda'

    INTEGER                                  :: cufft_dir
    LOGICAL                                  :: tscal

!     ==--------------------------------------------------------------==

    IF( plan%device /= blas_handle%device .OR. &
         plan%device /= stream%device ) CALL stopgm(procedureN,'devices are not consistent',&
         __LINE__,__FILE__)


    IF (isign.EQ.1) THEN
       cufft_dir=cufft_forward
    ELSE
       cufft_dir=cufft_inverse
    ENDIF

    CALL cufft_execz2z(plan, a, b, cufft_dir)

    tscal=(ABS(scale-1._real_8).GT.1.e-12_real_8)
    IF( .TRUE. ) THEN
       IF(tscal) CALL cublas_zdscal( blas_handle, ldax*lday, scale, b, 1); !acm: should it be ldbx*ldby ?
       CALL CuUser_setblock2zero ( b, transb, n, m, ldbx, ldby, stream )
    ELSE
       IF(tscal) THEN
          CALL CuUser_setblock2zero_scale ( b, scale, transb, n, m, ldbx, ldby, stream )
       ELSE
          CALL CuUser_setblock2zero ( b, transb, n, m, ldbx, ldby, stream )
       ENDIF
    ENDIF

  END SUBROUTINE mltfft_cuda

END MODULE mltfft_utils
