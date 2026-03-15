#include "cpmd_global.h"

MODULE summat_utils
  USE error_handling,                  ONLY: stopgm
  USE kinds,                           ONLY: int_4,&
                                             int_8,&
                                             real_8
  USE mp_interface,                    ONLY: mp_sum
  USE nlps,                            ONLY: imagp
  USE parac,                           ONLY: parai,&
                                             paral
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: symmat_pack,&
                                             symmat_unpack
#ifdef __PARALLEL
  USE mpi_f08
#endif
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: summat
  !!public :: sumhmat

CONTAINS

  ! ==================================================================
  SUBROUTINE summat(a,nstate,symmetrization,lsd,gid,parent)
    ! ==--------------------------------------------------------------==
    INTEGER, INTENT(IN)                      :: nstate
    REAL(real_8),INTENT(INOUT)               :: a(nstate,nstate)
#ifdef __PARALLEL
    type(MPI_COMM),INTENT(IN),OPTIONAL       :: gid
#else
    INTEGER,INTENT(IN),OPTIONAL              :: gid
#endif
    LOGICAL,INTENT(IN),OPTIONAL              :: symmetrization,lsd,parent

    CHARACTER(*), PARAMETER                  :: procedureN = 'summat'

    INTEGER                                  :: ierr, isub
    INTEGER(int_8)                           :: il_aux(1),il_aux1(1)
    LOGICAL                                  :: full,lsd_active,is_parent,&
         only_parent
#ifdef __PARALLEL
    type(MPI_COMM)                           :: mpi_com
    INTEGER                                  :: parent_rank
#else
    INTEGER                                  :: mpi_com,parent_rank
#endif
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8), POINTER __CONTIGUOUS       :: aux(:),aux1(:)
#else
    REAL(real_8), ALLOCATABLE                :: aux(:)
#endif
    ! ==--------------------------------------------------------------==
    ! == GLOBAL SUMMATION OF A SYMMETRIC MATRIX                       ==
    ! ==--------------------------------------------------------------==
    ! Modified: Tobias Kloeffel, Erlangen
    ! Date March 2019
    ! symmetrization can be disabled via optional flag symmetrization
    ! communicator can be changed via optional flag gid
    ! summat_parent is obsolete via optional flag parent
    ! able to pack both spins in a single mpi call via optional flag lsd
    ! without optional flags, returns to the original version com=allgrp
    ! symmetrization = .true.
    ! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    !check input parameters
    IF(PRESENT(gid))THEN
       mpi_com=gid
    ELSE
       mpi_com=parai%allgrp
    END IF
    IF(PRESENT(symmetrization))THEN
       full=symmetrization
    ELSE
       full=.TRUE.
    END IF
    IF(PRESENT(lsd))THEN
       lsd_active=lsd
    ELSE
       lsd_active=.FALSE.
    END IF
    IF(PRESENT(parent))THEN
       only_parent=parent
    ELSE
       only_parent=.FALSE.
    END IF
    !end check

    IF (mpi_com .EQ. parai%cp_grp) THEN
       is_parent=paral%io_parent
       parent_rank=parai%io_source
    ELSEIF (mpi_com .EQ. parai%allgrp) THEN
       is_parent=paral%parent
       parent_rank=parai%source
    ELSE
       CALL stopgm(procedureN,'unsupported communicator', &
         __LINE__,__FILE__)
    END IF
    IF(lsd_active.AND.cntl%tlsd)THEN
       il_aux(1)=spin_mod%nsup*(spin_mod%nsup+1)/2+&
         spin_mod%nsdown*(spin_mod%nsdown+1)/2
    ELSE
       il_aux(1)=nstate*(nstate+1)/2
    END IF
    il_aux1=il_aux
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_aux,aux,procedureN//'_aux',ierr)
#else
    ALLOCATE(aux(il_aux(1)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_aux1,aux1,procedureN//'_aux1',ierr)
#else
    ALLOCATE(aux1(il_aux1(1)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
         __LINE__,__FILE__)
    IF(lsd_active.AND.cntl%tlsd)THEN
       CALL symmat_pack(a,aux,nstate,spin_mod%nsup,spin_mod%nsdown)
    ELSE
       CALL symmat_pack(a,aux,nstate,nstate,0)
    END IF
    IF(only_parent)THEN
       CALL mp_sum(aux,aux1,INT(il_aux(1),KIND=int_4),parent_rank,mpi_com)
    ELSE
       CALL mp_sum(aux,aux1,INT(il_aux(1),KIND=int_4),mpi_com)
    END IF
    IF(.NOT.(only_parent.AND..NOT.is_parent))THEN
       IF(lsd_active.AND.cntl%tlsd)THEN
          CALL symmat_unpack(a,aux1,nstate,spin_mod%nsup,spin_mod%nsdown,full)
       ELSE
          CALL symmat_unpack(a,aux1,nstate,nstate,0,full)
       END IF
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_aux1,aux1,procedureN//'_aux1',ierr)
#else
    DEALLOCATE(aux1,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_aux,aux,procedureN//'_aux',ierr)
#else
    DEALLOCATE(aux,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
         __LINE__,__FILE__)

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE summat

END MODULE summat_utils

SUBROUTINE sumhmat(a,nstate)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_sum
  USE parac, ONLY : paral,parai
  IMPLICIT NONE
  INTEGER                                    :: nstate
  COMPLEX(real_8)                            :: a(nstate,nstate)

  CHARACTER(*), PARAMETER                    :: procedureN = 'sumhmat'

  COMPLEX(real_8), ALLOCATABLE               :: aux(:,:)
  INTEGER                                    :: i, ierr, isub, j, k, n2, ntr, &
                                                ntr2

! Variables
! ==--------------------------------------------------------------==
! == GLOBAL SUMMATION OF A HERMITIAN MATRIX                       ==
! ==--------------------------------------------------------------==

  IF (parai%nproc.LE.1) RETURN
  CALL tiset('   SUMHMAT',isub)
  n2 = nstate*nstate
  ntr = (nstate*(nstate+1))/2
  ntr2 = ntr*2
  ALLOCATE(aux(ntr2,1),STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'allocation problem', &
       __LINE__,__FILE__)
  k=0
  DO i=1,nstate
     DO j=i,nstate
        k=k+1
        aux(k,1)=a(i,j)
     ENDDO
  ENDDO
  CALL mp_sum(aux,a,ntr,parai%allgrp)
  CALL dcopy(ntr2,a(1,1),1,aux,1)
  k=0
  DO i=1,nstate
     DO j=i,nstate
        k=k+1
        a(i,j)=aux(k,1)
        IF (i.NE.j) a(j,i)=CMPLX(REAL(aux(k,1)),-AIMAG(aux(k,1)),kind=real_8)
     ENDDO
  ENDDO
  DEALLOCATE(aux,STAT=ierr)
  IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem', &
       __LINE__,__FILE__)
  CALL tihalt('   SUMHMAT',isub)
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE sumhmat
! ==================================================================
