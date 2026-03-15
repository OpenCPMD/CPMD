#include "cpmd_global.h"
#if defined(__FFT_HASNT_THREADED_COPIES)
#define HASNT_THREADED_COPIES .TRUE.
#else
#define HASNT_THREADED_COPIES .FALSE.
#endif

#if defined(__FFT_HAS_LOW_LEVEL_TIMERS)
#define HAS_LOW_LEVEL_TIMERS .TRUE.
#else
#define HAS_LOW_LEVEL_TIMERS .FALSE.
#endif

#if defined(__FFT_HAS_SPECIAL_COPY)
#define HAS_SPECIAL_COPY __FFT_HAS_SPECIAL_COPY
#else
#define HAS_SPECIAL_COPY 0
#endif

#ifndef __COLLAPSE2
#if defined(__FFT_HAS_OMP_COLLAPSE)
#define __COLLAPSE2 collapse(2)
#else
#define __COLLAPSE2
#endif
#endif

MODULE fftutil_utils
  USE fft_maxfft,                      ONLY: maxfft
  USE kinds,                           ONLY: real_4,&
                                             real_8
  USE mp_interface,                    ONLY: mp_all2all
  USE parac,                           ONLY: parai
  USE reshaper,                        ONLY: type_cast,&
                                             reshape_inplace
  USE system,                          ONLY: cntl,&
                                             fpar,&
                                             parap,&
                                             spar
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: zgthr_no_omp,&
                                             zsctr_no_omp
  USE zeroing_utils,                   ONLY: zeroing

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: phase
  PUBLIC :: putz
  PUBLIC :: getz
  PUBLIC :: unpack_y2x
  PUBLIC :: pack_y2x
  PUBLIC :: fft_comm
  PUBLIC :: pack_x2y
  PUBLIC :: unpack_x2y
  PUBLIC :: phasen
!TK special routines for batched fft
  PUBLIC :: getz_n
  PUBLIC :: putz_n
  PUBLIC :: pack_x2y_n
  PUBLIC :: unpack_x2y_n
  PUBLIC :: pack_y2x_n
  PUBLIC :: unpack_y2x_n
!TK
CONTAINS


  ! ==================================================================
  SUBROUTINE phase(f)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8) :: f(fpar%kr1,fpar%kr2s,fpar%kr3s)

    INTEGER                                  :: ii, ijk1, ijk2, isub, j, k

! ==--------------------------------------------------------------==

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset('     PHASE',isub)
    ijk2=parap%nrxpl(parai%mepos,2)-parap%nrxpl(parai%mepos,1)+1 ! jh-mb
    !$omp parallel do private (J,K,II,IJK1) shared(IJK2) __COLLAPSE2
    DO k=1,spar%nr3s
       DO j=1,spar%nr2s
          ii=k+j+parap%nrxpl(parai%mepos,1)
          ijk1=MOD(ii,2)+1
          DO ii=ijk1,ijk2,2
             f(ii,j,k)=-f(ii,j,k)
          ENDDO
       ENDDO
    ENDDO
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt('     PHASE',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE phase
  ! ==================================================================
  SUBROUTINE putz(a,b,krmin,krmax,kr,m)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: krmin, krmax
    COMPLEX(real_8)                          :: a(krmax-krmin+1,*)
    INTEGER                                  :: kr, m
    COMPLEX(real_8)                          :: b(kr,m)

    CHARACTER(*), PARAMETER                  :: procedureN = 'putz'

    INTEGER                                  :: isub, n

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub)
    n=krmax-krmin+1
    CALL zeroing(b)!,m*kr)
    CALL matmov(n,m,a,n,b(krmin,1),kr)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE putz
  ! ==================================================================
  SUBROUTINE getz(a,b,krmin,krmax,kr,m)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: krmin, krmax
    COMPLEX(real_8)                          :: b(krmax-krmin+1,*)
    INTEGER                                  :: kr
    COMPLEX(real_8)                          :: a(kr,*)
    INTEGER                                  :: m

    CHARACTER(*), PARAMETER                  :: procedureN = 'getz'

    INTEGER                                  :: isub, n

! ==--------------------------------------------------------------==

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub)
    n=krmax-krmin+1
    CALL matmov(n,m,a(krmin,1),kr,b,n)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE getz
  ! ==================================================================
  ! ==================================================================
  ! CODE FOR NEW FFT ROUTINES
  ! ==================================================================


  SUBROUTINE unpack_y2x(xf,yf,m,nrays,lda,jrxpl,sp5,maxfft,mproc,&
       tr4a2a)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: xf(*), yf(*)
    INTEGER                                  :: m, nrays, lda, maxfft, mproc, &
                                                sp5(0:mproc-1), &
                                                jrxpl(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: yf4(:)
    INTEGER                                  :: ip, ipp, isub1, k, nrs, nrx

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='UNPACK_Y2X'
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! ..Pack the data for sending
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(yf, maxfft, yf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL dcopy_s(2*nrx,yf4(ipp),1,xf(nrs),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                xf(nrs+k)=yf4(ipp+k)
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL dcopy(2*nrx,yf(ipp),1,xf(nrs),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                xf(nrs+k)=yf(ipp+k)
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE unpack_y2x
  ! ==================================================================
  SUBROUTINE pack_y2x(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,&
       tr4a2a)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: xf(*), yf(*)
    INTEGER                                  :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), maxfft, mproc, &
                                                sp8(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: xf4(:)
    INTEGER                                  :: i, ii, ip, isub1, jj, k, &
                                                mxrp, nrx

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='PACK_Y2X'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    ! ..Pack the data for sending
    ! IF(HAS_LOW_LEVEL_TIMERS) CALL TISET(procedureN,ISUB1)
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(xf, maxfft, xf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL cgthr_z(mxrp,yf(jj),xf4(ii),msp(1,ip+1))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   xf4(ii+k) = yf(jj+msp(k,ip+1))
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default

          !$omp parallel do shared(SP8,NRX,M,LDA) &
          !$omp             private(MXRP,I,II,JJ)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL zgthr_no_omp(mxrp,yf(jj),xf(ii),msp(1,ip+1))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   xf(ii+k) = yf(jj+msp(k,ip+1))
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pack_y2x
  ! ==================================================================
  SUBROUTINE fft_comm(xf,yf,lda,tr4a2a, comm )
    ! ==--------------------------------------------------------------==
#ifdef __PARALLEL
    USE mpi_f08
#endif
    COMPLEX(real_8), TARGET                  :: xf(*), yf(*)
    INTEGER, INTENT(IN)                      :: lda
    LOGICAL, INTENT(IN)                      :: tr4a2a
#ifdef __PARALLEL
    type(MPI_COMM), INTENT(IN)                      :: comm
#else
    INTEGER, INTENT(IN)                      :: comm
#endif

    CHARACTER(*), PARAMETER                  :: procedureN = 'fft_comm'

    COMPLEX(real_4), POINTER                 :: xf4(:), yf4(:)
    INTEGER                                  :: isub1, isub2

! Variables
! ==--------------------------------------------------------------==
! ..All to all communication

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF
    IF (tr4a2a) THEN
       ! is that needed on P?
       CALL type_cast(xf, maxfft, xf4)
       CALL type_cast(yf, maxfft, yf4)       
       CALL mp_all2all( xf4, yf4, lda, comm )
    ELSE
       CALL mp_all2all( xf, yf, lda, comm )
    ENDIF
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE fft_comm
  ! ==================================================================
  SUBROUTINE pack_x2y(xf,yf,nrays,lda,jrxpl,sp5,&
       maxfft,mproc,tr4a2a)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: xf(*), yf(*)
    INTEGER                                  :: nrays, lda, maxfft, mproc, &
                                                sp5(0:mproc-1), &
                                                jrxpl(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: yf4(:)
    INTEGER                                  :: ip, ipp, isub1, k, nrs, nrx

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='PACK_X2Y'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    ! ..Prepare data for sending
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(yf, maxfft, yf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL scopy_d(2*nrx,xf(nrs),1,yf4(ipp),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                yf4(ipp+k)=xf(nrs+k)
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays + 1
             ipp = ip*lda + 1
             CALL dcopy(2*nrx,xf(nrs),1,yf(ipp),1)
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP5,NRAYS,JRXPL,LDA) &
          !$omp             private(NRX,NRS,IPP,K)
          DO ip=0,mproc-1
             nrx = sp5(ip)*nrays
             nrs = (jrxpl(ip)-1)*nrays
             ipp = ip*lda
             DO k=1,nrx
                yf(ipp+k)=xf(nrs+k)
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pack_x2y
  ! ==================================================================
  SUBROUTINE unpack_x2y(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,tr4a2a)
    ! ==--------------------------------------------------------------==
    ! include 'parac.inc'
    COMPLEX(real_8)                          :: xf(*)
    INTEGER                                  :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), maxfft
    COMPLEX(real_8)                          :: yf(maxfft)
    INTEGER                                  :: mproc, sp8(0:mproc-1)
    LOGICAL                                  :: tr4a2a

    COMPLEX(real_4), POINTER                 :: xf4(:)
    INTEGER                                  :: i, ii, ip, isub1, jj, k, mxrp

    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='UNPACK_X2Y'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset(procedureN,isub1)
    ! ..Unpacking the data
    CALL zeroing(yf)!,maxfft)
    !$    IF(HASNT_THREADED_COPIES) THEN
    !$       max_threads = omp_get_max_threads()
    !$       call omp_set_num_threads(1)
    !$    ENDIF
    IF (tr4a2a) THEN
       CALL type_cast(xf, maxfft, xf4)
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL zsctr_c(mxrp,xf4(ii),msp(1,ip+1),yf(jj))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   yf(jj+msp(k,ip+1)) = xf4(ii+k)
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ELSE
       SELECT CASE(HAS_SPECIAL_COPY)
       CASE default
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp + 1
                jj = (i-1)*m + 1
                CALL zsctr_no_omp(mxrp,xf(ii),msp(1,ip+1),yf(jj))
             ENDDO
          ENDDO
       CASE(1)
          ! that seems to work better on the P
          !$omp parallel do shared(SP8,LR1,M,LDA) &
          !$omp             private(MXRP,I,II,JJ,K)
          DO ip=0,mproc-1
             mxrp = sp8(ip)
             DO i=1,lr1
                ii = ip*lda + (i-1)*mxrp
                jj = (i-1)*m
                DO k=1,mxrp
                   yf(jj+msp(k,ip+1)) = xf(ii+k)
                ENDDO
             ENDDO
          ENDDO
       END SELECT
    ENDIF
    !$    IF(HASNT_THREADED_COPIES) call omp_set_num_threads(max_threads)
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt(procedureN,isub1)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE unpack_x2y
  ! ==================================================================
  SUBROUTINE phasen(f,kr1,kr2s,kr3s,n1u,n1o,nr2s,nr3s)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: kr1, kr2s, kr3s
    COMPLEX(real_8)                          :: f(kr1,kr2s,kr3s)
    INTEGER                                  :: n1u, n1o, nr2s, nr3s

    INTEGER                                  :: i, ii, ijk, isub, j, k
    REAL(real_8), DIMENSION(2)               :: pf = (/1._real_8,-1._real_8/)

    IF (HAS_LOW_LEVEL_TIMERS) CALL tiset('     PHASE',isub)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute parallel do collapse(3) &
#else
    !$omp parallel do default(none) __COLLAPSE2 &
#endif
    !$omp             private(K,J,I,II,IJK) &
    !$omp             shared(F,PF,NR3S,NR2S,N1U,N1O)
    DO k=1,nr3s
       DO j=1,nr2s
          DO i=n1u,n1o
             ii=i-n1u+1
             ijk=MOD(k+j+i+1,2)+1
             f(ii,j,k)=f(ii,j,k)*pf(ijk)
          ENDDO
       ENDDO
    ENDDO
    IF (HAS_LOW_LEVEL_TIMERS) CALL tihalt('     PHASE',isub)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE phasen
  ! ==================================================================

  !TK
  ! Rewritten routines taking care of batches of states
  ! no support for A2A in single precision
  ! Author:
  ! Tobias Kloeffel, CCC,FAU Erlangen-Nuernberg tobias.kloeffel@fau.de
  ! Gerald Mathias, LRZ, Garching Gerald.Mathias@lrz.de
  ! Bernd Meyer, CCC, FAU Erlangen-Nuernberg bernd.meyer@fau.de

  SUBROUTINE putz_n(a,b,krmin,krmax,kr,kr1,kr2s,nperbatch)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: krmin, krmax, kr, kr1, kr2s, nperbatch
    COMPLEX(real_8),INTENT(IN)               :: a(krmax-krmin+1,kr1,nperbatch,*)
    COMPLEX(real_8),INTENT(INOUT)            :: b(kr,kr1,kr2s,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'putz_n'
    REAL(real_8), POINTER __CONTIGUOUS       :: b_r(:,:,:,:)
    REAL(real_8), POINTER __CONTIGUOUS       :: a_r(:,:,:,:)
    INTEGER                                  :: isub1,isub2

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF
    
    CALL reshape_inplace(a,(/(krmax-krmin+1)*2,kr1,nperbatch,kr2s/),a_r)
    CALL reshape_inplace(b,(/kr*2,kr1,kr2s,nperbatch/),b_r)
    CALL putz_n_r(a_r,b_r,krmin,krmax,kr,kr1,kr2s,nperbatch)

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE putz_n

  SUBROUTINE getz_n(a,b,krmin,krmax,kr,kr1,kr2s,nperbatch)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: krmin, krmax, kr, kr1, kr2s, nperbatch
    COMPLEX(real_8),INTENT(OUT)              :: b(krmax-krmin+1,kr1,nperbatch,*)
    COMPLEX(real_8),INTENT(IN)               :: a(kr,kr1,kr2s,*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'getz_n'

    INTEGER                                  :: isub1,isub2,n,is,i,j,k,n1

! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF

    n=krmax-krmin+1
    n1=krmin-1
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute parallel do collapse(4) &
#else
    !$omp parallel proc_bind(close) &
#endif
    !$omp& private (i,is,k,j)
    DO is=1,nperbatch
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp do
#endif
       DO i=1,kr2s
          DO k=1,kr1
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp simd
#endif
             DO j=1,n
                b(j,k,is,i)=a(n1+j,k,i,is)
             END DO
          END DO
       END DO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp end do nowait
#endif
    END DO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end parallel
#endif

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE getz_n

  SUBROUTINE pack_y2x_n(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,&
       tr4a2a,nperbatch,offset_in)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(OUT)              :: xf(*)
    COMPLEX(real_8),INTENT(IN)               :: yf(*)
    INTEGER,INTENT(IN)                       :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), mproc, &
                                                sp8(0:mproc-1), maxfft
    LOGICAL,INTENT(IN)                       :: tr4a2a
    INTEGER,INTENT(IN),OPTIONAL              :: nperbatch,offset_in

    INTEGER                                  :: i, ii, ip, isub1, isub2, &
                                                jj, k, mxrp, is, nstate, offset

    CHARACTER(*),PARAMETER :: procedureN='PACK_Y2X_n'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF

    IF(PRESENT(nperbatch))THEN
       nstate=nperbatch
    ELSE
       nstate=1
    END IF
    IF(PRESENT(offset_in))THEN
       offset=offset_in
    ELSE
       offset=0
    END IF

#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute collapse(3) &
#else
    !$omp parallel proc_bind(close) &
#endif
    !$omp& private(is,ip,mxrp,i,ii,jj,k) 
    DO is=1,nstate
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp do
#endif
       DO ip=0,mproc-1
          DO i=1,lr1
             mxrp = sp8(ip)
             ii = ip*lda*nstate + (i-1)*mxrp + (is-1)*lda
             jj = (i-1)*m + (is-1)*offset
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp parallel do simd
#else
!             !$omp simd
#endif
             DO k=1,mxrp
                xf(ii+k) = yf(jj+msp(k,ip+1))
             ENDDO
          ENDDO
       ENDDO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp end do nowait
#endif
    end do
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end parallel
#endif

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pack_y2x_n

  SUBROUTINE unpack_y2x_n(xf,yf,m,nrays,lda,jrxpl,sp5,maxfft,mproc,tr4a2a,count)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(IN)               :: yf(*)
    COMPLEX(real_8),INTENT(OUT)              :: xf(*)
    LOGICAL,INTENT(IN)                       :: tr4a2a
    INTEGER,INTENT(IN)                       :: m, nrays, lda, mproc, &
                                                sp5(0:mproc-1), &
                                                jrxpl(0:mproc-1), maxfft
    INTEGER,INTENT(IN),OPTIONAL              :: count

    INTEGER                                  :: ip, ipp, isub1, isub2, k, nrs, &
                                                nrx, i,is,nstate,sp5_max

    CHARACTER(*),PARAMETER :: procedureN='UNPACK_Y2X_n'
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! ..Pack the data for sending
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF

    If(PRESENT(count))THEN
       nstate=count
    ELSE
       nstate=1
    END IF


#if defined(_HAS_OMP_TARGET_OFFLOAD)
    sp5_max=maxval(sp5)
    !$omp target teams distribute collapse(3) &
#else
    !$omp parallel do proc_bind(close) &
#endif
    !$omp& private(ip,i,is,nrs,ipp,k)
    DO ip=0,mproc-1
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       DO i=1,sp5_max
#else
       DO i=1,sp5(ip)
#endif
          DO is=1,nstate
             IF(i.GT.sp5(ip))CYCLE
             nrs = (jrxpl(ip)-1)*nrays*nstate + (i-1)*nrays*nstate +(is-1)*nrays
             ipp = ip*lda*nstate + (i-1)*nrays + (is-1)*lda
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp parallel do simd
#else
             !$omp simd
#endif
             DO k=1,nrays
                xf(nrs+k)=yf(ipp+k)
             END DO
          END DO
       END DO
    END DO

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE unpack_y2x_n
  ! ==================================================================

  SUBROUTINE pack_x2y_n(xf,yf,nrays,lda,jrxpl,sp5,maxfft,mproc,tr4a2a,count)
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(IN)               :: xf(*)
    COMPLEX(real_8),INTENT(OUT)              :: yf(*)
    INTEGER,INTENT(IN)                       :: nrays, lda, maxfft, mproc, &
                                                sp5(0:mproc-1), &
                                                jrxpl(0:mproc-1)
    LOGICAL,INTENT(IN)                       :: tr4a2a
    INTEGER,INTENT(IN),OPTIONAL              :: count

    INTEGER                                  :: ip, ipp, isub1, isub2, &
                                                 k, nrs, nrx, i,is,nstate,&
                                                 sp5_max
    CHARACTER(*),PARAMETER :: procedureN='PACK_X2Y_n'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF

    IF(PRESENT(count))THEN
       nstate=count
    ELSE
       nstate=1
    END IF

    sp5_max=maxval(sp5)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute collapse(3) &
#else
    !$omp parallel do proc_bind(close) &
#endif
    !$omp& private(ip,i,is,nrs,ipp,k) 
    DO ip=0,mproc-1
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       DO i=1,sp5_max
#else
       DO i=1,sp5(ip)
#endif
          DO is=1,nstate
             IF(i.GT.sp5(ip))CYCLE
             nrs = (jrxpl(ip)-1)*nrays*nstate + (i-1)*nrays*nstate +(is-1)*nrays
             ipp = ip*lda*nstate + (i-1)*nrays + (is-1)*lda
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp parallel do simd
#else
             !$omp simd
#endif
             DO k=1,nrays
                yf(ipp+k)=xf(nrs+k)
             END DO
          END DO
       END DO
    END DO

    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE pack_x2y_n

  ! ==================================================================
  SUBROUTINE unpack_x2y_n(xf,yf,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,tr4a2a,&
       count,offset_in)
    ! ==--------------------------------------------------------------==
    ! include 'parac.inc'
    COMPLEX(real_8),INTENT(IN)               :: xf(*)
    COMPLEX(real_8),INTENT(INOUT)            :: yf(*)
    INTEGER, INTENT(IN)                      :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), maxfft,&
                                                mproc, sp8(0:mproc-1)
    LOGICAL, INTENT(IN)                      :: tr4a2a
    INTEGER, INTENT(IN), OPTIONAL            :: count,offset_in
    INTEGER                                  :: isub1, isub2, nstate, offset
    REAL(real_8), POINTER __CONTIGUOUS       :: yf_r(:),xf_r(:)
    
    !$    INTEGER   max_threads
    !$    INTEGER, EXTERNAL :: omp_get_max_threads
    CHARACTER(*),PARAMETER :: procedureN='UNPACK_X2Y_n'
    ! ==--------------------------------------------------------------==
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tiset(procedureN//'_tuning',isub1)
       ELSE
          CALL tiset(procedureN,isub2)
       END IF
    END IF
    ! ..Unpacking the data
    IF(PRESENT(count))THEN
       nstate=count
    ELSE
       nstate=1
    ENDIF
    IF(PRESENT(offset_in))THEN
       offset=offset_in
    ELSE
       offset=0
    END IF
    CALL reshape_inplace(yf,(/nstate*offset*2/),yf_r)
    CALL reshape_inplace(xf,(/nstate*offset*2/),xf_r)
    CALL unpack_x2y_n_r(xf_r,yf_r,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,&
       offset,nstate)
    IF (HAS_LOW_LEVEL_TIMERS) THEN
       IF(cntl%fft_tune_batchsize) THEN
          CALL tihalt(procedureN//'_tuning',isub1)
       ELSE
          CALL tihalt(procedureN,isub2)
       END IF
    END IF
    ! ==--------------------------------------------------------------==
  END SUBROUTINE unpack_x2y_n
  ! ==================================================================
  SUBROUTINE unpack_x2y_n_r(xf_r,yf_r,m,lr1,lda,msp,lmsp,sp8,maxfft,mproc,&
       offset,nstate)
    ! ==--------------------------------------------------------------==
    ! include 'parac.inc'
    real(real_8),INTENT(IN)               :: xf_r(*)
    real(real_8),INTENT(OUT)              :: yf_r(*)
    INTEGER, INTENT(IN)                      :: m, lr1, lda, lmsp, &
                                                msp(lmsp,*), maxfft,&
                                                mproc, sp8(0:mproc-1)
    INTEGER, INTENT(IN)                    :: offset,nstate
    INTEGER                                  :: i, ii, ip, isub1, jj, k, mxrp,is,kk

#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute parallel do simd
#else
    !$omp parallel private(is,ip,mxrp,i,ii,jj,k,kk) proc_bind(close)
    !$omp do
#endif
    do is=1,nstate*offset*2
       yf_r(is)=0._real_8
    end do
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute collapse(3) private(is,ip,i,mxrp,ii,jj,k,kk)
#endif
    DO is=1,nstate
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp do schedule (static)
#endif
       DO ip=0,mproc-1
          DO i=1,lr1
             mxrp = sp8(ip)
             ii = ip*lda*nstate + (i-1)*mxrp + (is-1)*lda
             jj = (i-1)*m + (is-1)*offset
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp parallel do simd
#endif
             DO k=1,mxrp
                do kk=-1,0
                   yf_r((jj+msp(k,ip+1))*2+kk) = xf_r((ii+k)*2+kk)
                end do
             END DO
          END DO
       END DO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp end do nowait
#endif
    END DO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end parallel
#endif

  END SUBROUTINE unpack_x2y_n_r
  SUBROUTINE putz_n_r(a_r,b_r,krmin,krmax,kr,kr1,kr2s,nperbatch)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: krmin, krmax, kr, kr1, kr2s, nperbatch
    REAL(real_8),INTENT(IN)               :: a_r((krmax-krmin+1)*2,kr1,nperbatch,*)
    REAL(real_8),INTENT(OUT)              :: b_r(kr*2,kr1,kr2s,*)

    INTEGER                                  :: isub,n,n1,n2,n3,n4,is,i,j,k,krmin_loc,kr_loc

    
    n=krmax-krmin+1
    n1=(krmin-1)*2
    n2=(krmin-1+n)*2
    n3=(krmin+n)*2-1
    n4=(krmin-1)*2
    kr_loc=kr*2
    krmin_loc=krmin*2-1

#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute collapse(3) &
#else
    !$omp parallel proc_bind(close) &
#endif
    !$omp& private (i,is,k,j)
    DO is=1,nperbatch
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp do
#endif
       DO i=1,kr2s
          DO k=1,kr1
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp parallel private(j)
             !$omp do simd
#else
             !$omp simd
#endif
             DO j=1,n1
                b_r(j,k,i,is)=0.0_real_8
             END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp end do simd nowait
             !$omp do simd
#else
             !$omp simd
#endif
             DO j=krmin_loc,n2
                b_r(j,k,i,is)=a_r(j-n4,k,is,i)
             END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp end do simd nowait
             !$omp do simd
#else
             !$omp simd
#endif
             DO j=n3,kr_loc
                b_r(j,k,i,is)=0.0_real_8
             END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp end do simd nowait
             !$omp end parallel
#endif
          END DO
       END DO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp end do nowait
#endif
    END DO
#if !defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp end parallel
#endif
    ! ==--------------------------------------------------------------==
  END SUBROUTINE putz_n_r
END MODULE fftutil_utils
