#include "cpmd_global.h"

MODULE fnl_utils
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nlps_com
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: unpack_fnl
  PUBLIC :: unpack_dfnl
  PUBLIC :: pack_fnl
  PUBLIC :: pack_dfnl
  PUBLIC :: sort_fnl
  PUBLIC :: sort_fnl_k
  PUBLIC :: sort_dfnl
  PUBLIC :: sort_dfnl_k

CONTAINS
  ! ==================================================================
  SUBROUTINE unpack_dfnl(na,packed,unpacked_k,unpacked)
    ! ==--------------------------------------------------------------==
    ! == unpack packed_dfnl, relies on custom na                      ==
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, FAU Erlangen-Nuernberg March 2019
    INTEGER,INTENT(IN)                       :: na(:,:)
    REAL(real_8),INTENT(IN)  __CONTIGUOUS    :: packed(:,:)
    REAL(real_8),INTENT(OUT), OPTIONAL&
                             __CONTIGUOUS    :: unpacked_k(:,:,:,:,:),unpacked(:,:,:,:)
    INTEGER                                  :: i,offset,isa0,is,iv,ia,isa,k,im

    IF(PRESENT(unpacked_k))THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute parallel do &
#else
       !$omp parallel do proc_bind(close) &
#endif
       !$omp& private(i,offset,isa0,is,iv,k,ia,isa,im)
       DO i=1,SIZE(packed,2)
          offset=0
          isa0=0
          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                DO k=1,3
                   DO ia=na(1,is),na(2,is)
                      isa=isa0+ia
                      DO im=1,2
                         offset=offset+1
                         unpacked_k(im,isa,iv,k,i)=packed(offset,i)
                      END DO
                   END DO
                END DO
             END DO
             isa0=isa0+ions0%na(is)
          END DO
       END DO
    ELSEIF(PRESENT(unpacked))THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute parallel do &
#else
       !$omp parallel do proc_bind(close) &
#endif
       !$omp& private(i,offset,k,isa0,is,iv,ia,isa) 
       DO i=1,SIZE(packed,2)
          offset=0
          isa0=0
          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                DO k=1,3
                   Do ia=na(1,is),na(2,is)
                      isa=isa0+ia
                      offset=offset+1
                      unpacked(isa,iv,k,i)=packed(offset,i)
                   END DO
                END DO
             END DO
             isa0=isa0+ions0%na(is)
          END DO
       END DO
    END IF
  END SUBROUTINE unpack_dfnl
  ! ==================================================================
  SUBROUTINE unpack_fnl(na,packed,unpacked_k,unpacked)
    ! ==--------------------------------------------------------------==
    ! == unpack packed_fnl, relies on custom na                       ==
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, FAU Erlangen-Nuernberg March 2019
    INTEGER,INTENT(IN)                       :: na(:,:)
    REAL(real_8),INTENT(IN)  __CONTIGUOUS    :: packed(:,:)
    REAL(real_8),INTENT(OUT), OPTIONAL&
                             __CONTIGUOUS    :: unpacked_k(:,:,:,:),unpacked(:,:,:)
    INTEGER                                  :: i,offset,isa0,is,iv,ia,isa,k,im

    IF(PRESENT(unpacked_k))THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute parallel do &
#else       
       !$omp parallel do proc_bind(close) &
#endif
       !$omp& private(i,offset,isa0,is,iv,ia,isa,im)
       DO i=1,SIZE(packed,2)
          isa0=0
          offset=0
          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                DO ia=na(1,is),na(2,is)
                   isa=isa0+ia
                   DO im=1,2
                      offset=offset+1
                      unpacked_k(im,isa,iv,i)=packed(offset,i)
                   END DO
                END DO
             END DO
             isa0=isa0+ions0%na(is)
          END DO
       END DO
    ELSEIF(PRESENT(unpacked))THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute parallel do &
#else       
       !$omp parallel do proc_bind(close) &
#endif
       !$omp& private(i,offset,isa0,is,iv,ia,isa) 
       DO i=1,SIZE(packed,2)
          isa0=0
          offset=0
          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                DO ia=na(1,is),na(2,is)
                   isa=isa0+ia
                   offset=offset+1
                   unpacked(isa,iv,i)=packed(offset,i)
                END DO
             END DO
             isa0=isa0+ions0%na(is)
          END DO
       END DO
    END IF
  END SUBROUTINE unpack_fnl
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE pack_dfnl(na,packed,unpacked_k,unpacked)
    ! ==--------------------------------------------------------------==
    ! == pack dfnl, relies on custom na                               ==
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, FAU Erlangen-Nuernberg March 2019
    INTEGER,INTENT(IN)                       :: na(:,:)
    REAL(real_8),INTENT(OUT)  __CONTIGUOUS   :: packed(:,:)
    REAL(real_8),INTENT(IN), OPTIONAL&
                             __CONTIGUOUS    :: unpacked_k(:,:,:,:,:),unpacked(:,:,:,:)
    INTEGER                                  :: i,offset,isa0,is,iv,ia,isa,k,im

    IF(PRESENT(unpacked_k))THEN
       !$omp parallel do private(i,offset,k,isa0,is,iv,ia,isa,im) proc_bind(close)
       DO i=1,SIZE(packed,2)
          offset=0
          isa0=0
          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                DO k=1,3
                   DO ia=na(1,is),na(2,is)
                      isa=isa0+ia
                      DO im=1,2
                         offset=offset+1
                         packed(offset,i)=unpacked_k(im,isa,iv,k,i)
                      END DO
                   END DO
                END DO
             END DO
             isa0=isa0+ions0%na(is)
          END DO
       END DO
    ELSEIF(PRESENT(unpacked))THEN
       !$omp parallel do private(i,offset,k,isa0,is,iv,ia,isa) proc_bind(close)
       DO i=1,SIZE(packed,2)
          offset=0
          isa0=0
          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                DO k=1,3
                   Do ia=na(1,is),na(2,is)
                      isa=isa0+ia
                      offset=offset+1
                      packed(offset,i)=unpacked(isa,iv,k,i)
                   END DO
                END DO
             END DO
             isa0=isa0+ions0%na(is)
          END DO
       END DO
    END IF
  END SUBROUTINE pack_dfnl
  ! ==================================================================
  SUBROUTINE pack_fnl(na,packed,unpacked_k,unpacked)
    ! ==--------------------------------------------------------------==
    ! == pack dfnl, relies on custom na                               ==
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, FAU Erlangen-Nuernberg March 2019
    INTEGER,INTENT(IN)                       :: na(:,:)
    REAL(real_8),INTENT(OUT)  __CONTIGUOUS   :: packed(:,:)
    REAL(real_8),INTENT(IN), OPTIONAL&
                             __CONTIGUOUS    :: unpacked_k(:,:,:,:),unpacked(:,:,:)
    INTEGER                                  :: i,offset,isa0,is,iv,ia,isa,k,im

    IF(PRESENT(unpacked_k))THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute parallel do &
#else
       !$omp parallel do proc_bind(close) &
#endif
       !$omp& private(i,offset,isa0,is,iv,ia,isa,im) 
       DO i=1,SIZE(packed,2)
          isa0=0
          offset=0
          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                DO ia=na(1,is),na(2,is)
                   isa=isa0+ia
                   DO im=1,2
                      offset=offset+1
                      packed(offset,i)=unpacked_k(im,isa,iv,i)
                   END DO
                END DO
             END DO
             isa0=isa0+ions0%na(is)
          END DO
       END DO
    ELSEIF(PRESENT(unpacked))THEN
#if defined(_HAS_OMP_TARGET_OFFLOAD)
       !$omp target teams distribute parallel do &
#else
       !$omp parallel do proc_bind(close) &
#endif
       !$omp& private(i,offset,isa0,is,iv,ia,isa)
       DO i=1,SIZE(packed,2)
          isa0=0
          offset=0
          DO is=1,ions1%nsp
             DO iv=1,nlps_com%ngh(is)
                DO ia=na(1,is),na(2,is)
                   isa=isa0+ia
                   offset=offset+1
                   packed(offset,i)=unpacked(isa,iv,i)
                END DO
             END DO
             isa0=isa0+ions0%na(is)
          END DO
       END DO
    END IF
  END SUBROUTINE pack_fnl
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE sort_fnl(nchunks,na,chunks,sorted,start_chunks,ld_chunks)
    ! ==--------------------------------------------------------------==
    ! == sorts fnl chunks according to na output: packed fnl          ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nchunks, na(2,ions1%nsp,nchunks),&
                                                start_chunks(nchunks), ld_chunks(nchunks)
    REAL(real_8),INTENT(IN)                  :: chunks(*)
    REAL(real_8),INTENT(OUT) __CONTIGUOUS    :: sorted(:,:)
    INTEGER                                  :: i, offset(nchunks), is, iv, ia, ichunk,&
                                                offset_sort
    ! ==--------------------------------------------------------------==
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute parallel do &
#else
    !$omp parallel do &
#endif
    !$omp& private(i,offset_sort,offset,is,iv,ichunk,ia)
    DO i=1,SIZE(sorted,2)
       offset_sort=0
       offset=(i-1)*ld_chunks+start_chunks-1
       DO is=1,ions1%nsp
          DO iv=1,nlps_com%ngh(is)
             DO ichunk=1,nchunks
                DO ia=na(1,is,ichunk),na(2,is,ichunk)
                   offset(ichunk)=offset(ichunk)+1
                   offset_sort=offset_sort+1
                   sorted(offset_sort,i)=chunks(offset(ichunk))
                END DO
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE sort_fnl
  ! ==================================================================
  SUBROUTINE sort_fnl_k(nchunks,na,chunks,sorted,start_chunks,ld_chunks)
    ! ==--------------------------------------------------------------==
    ! == sorts fnl chunks according to na output: packed fnl          ==
    ! == kpoint version                                               ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nchunks, na(2,ions1%nsp,nchunks), &
                                                start_chunks(nchunks), ld_chunks(nchunks)
    REAL(real_8),INTENT(IN)                  :: chunks(*)
    REAL(real_8),INTENT(OUT) __CONTIGUOUS    :: sorted(:,:)
    INTEGER                                  :: i, offset(nchunks), is, iv, ia, ichunk,&
                                                offset_sort, im
    ! ==--------------------------------------------------------------==
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute parallel do &
#else
    !$omp parallel do &
#endif
    !$omp& private(i,offset_sort,offset,is,iv,ichunk,ia,im)
    DO i=1,SIZE(sorted,2)
       offset_sort=0
       offset=(i-1)*ld_chunks+start_chunks-1
       DO is=1,ions1%nsp
          DO iv=1,nlps_com%ngh(is)
             DO ichunk=1,nchunks
                DO ia=na(1,is,ichunk),na(2,is,ichunk)
                   DO im=1,2
                      offset(ichunk)=offset(ichunk)+1
                      offset_sort=offset_sort+1
                      sorted(offset_sort,i)=chunks(offset(ichunk))
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE sort_fnl_k
  ! ==================================================================
  SUBROUTINE sort_dfnl(nchunks,na,chunks,sorted,start_chunks,ld_chunks)
    ! ==--------------------------------------------------------------==
    ! == sorts dfnl chunks according to na output: packed fnl         ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nchunks, na(2,ions1%nsp,nchunks),&
                                                start_chunks(nchunks), ld_chunks(nchunks)
    REAL(real_8),INTENT(IN)                  :: chunks(*)
    REAL(real_8),INTENT(OUT) __CONTIGUOUS    :: sorted(:,:)
    INTEGER                                  :: i, offset(nchunks), is, iv, ia, ichunk,&
                                                offset_sort, k
    ! ==--------------------------------------------------------------==
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute parallel do &
#else
    !$omp parallel do &
#endif
    !$omp& private(i,offset_sort,offset,is,iv,ichunk,ia,k)
    DO i=1,SIZE(sorted,2)
       offset_sort=0
       offset=(i-1)*ld_chunks+start_chunks-1
       DO is=1,ions1%nsp
          DO iv=1,nlps_com%ngh(is)
             DO k=1,3
                DO ichunk=1,nchunks
                   DO ia=na(1,is,ichunk),na(2,is,ichunk)
                      offset(ichunk)=offset(ichunk)+1
                      offset_sort=offset_sort+1
                      sorted(offset_sort,i)=chunks(offset(ichunk))
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE sort_dfnl
  ! ==================================================================
  SUBROUTINE sort_dfnl_k(nchunks,na,chunks,sorted,start_chunks,ld_chunks)
    ! ==--------------------------------------------------------------==
    ! == sorts dfnl chunks according to na output: packed fnl         ==
    ! == kpoint version                                               ==
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nchunks, na(2,ions1%nsp,nchunks), &
                                                start_chunks(nchunks), ld_chunks(nchunks)
    REAL(real_8),INTENT(IN)                  :: chunks(*)
    REAL(real_8),INTENT(OUT) __CONTIGUOUS    :: sorted(:,:)
    INTEGER                                  :: i, offset(nchunks), is, iv, ia, ichunk,&
                                                offset_sort, im, k
    ! ==--------------------------------------------------------------==
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp target teams distribute parallel do &
#else
    !$omp parallel do &
#endif
    !$omp& private(i,offset_sort,offset,is,iv,ichunk,ia,im,k)
    DO i=1,SIZE(sorted,2)
       offset_sort=0
       offset=(i-1)*ld_chunks+start_chunks-1
       DO is=1,ions1%nsp
          DO iv=1,nlps_com%ngh(is)
             DO k=1,3
                DO ichunk=1,nchunks
                   DO ia=na(1,is,ichunk),na(2,is,ichunk)
                      DO im=1,2
                         offset(ichunk)=offset(ichunk)+1
                         offset_sort=offset_sort+1
                         sorted(offset_sort,i)=chunks(offset(ichunk))
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE sort_dfnl_k
  ! ==================================================================
END MODULE fnl_utils
