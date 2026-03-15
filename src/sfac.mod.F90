#include "cpmd_global.h"

MODULE sfac
  USE kinds,                           ONLY: real_8,&
                                             int_8

  IMPLICIT NONE

  ! ==--------------------------------------------------------------==
  ! == FNL(IMAGP,NAT,NHXS,NSTATE,NKPNT)                             ==
  ! == DFNL(IMAGP,NAT,NHXS,3,NDFNL,NKPNT)                           ==
  ! ==--------------------------------------------------------------==


  REAL(real_8), POINTER __CONTIGUOUS :: fnl(:,:,:,:,:)
  REAL(real_8), POINTER __CONTIGUOUS :: fnl2(:,:,:,:,:)
  REAL(real_8), POINTER __CONTIGUOUS :: fnla(:,:,:)
  REAL(real_8), POINTER __CONTIGUOUS :: dfnl(:,:,:,:,:,:)
  REAL(real_8), POINTER __CONTIGUOUS :: dfnla(:,:,:,:)
  REAL(real_8), ALLOCATABLE :: ddfnl(:,:,:,:)
  REAL(real_8), POINTER __CONTIGUOUS :: fnl_packed(:,:)
  REAL(real_8), POINTER __CONTIGUOUS :: fnlgam_packed(:,:)
  REAL(real_8), POINTER __CONTIGUOUS :: dfnl_packed(:,:)

  INTEGER(int_8)                     :: il_dfnl_packed(2) = 0,&
                                        il_fnl_packed(2) = 0

  ! NOTE: not clear what happens with this var... 
  REAL(real_8), ALLOCATABLE :: FNLGP(:,:,:,:,:) ! ! FNLGP(IMAGP,NAT,NHXS,LDF2,NKPOINT)

  ! (:,:,1) because EIGR is used only if NKPNT=1. in K-points EIGKR is used instead 
  COMPLEX(real_8), ALLOCATABLE, TARGET :: eigr(:,:,:)

  COMPLEX(real_8), ALLOCATABLE, TARGET :: eigrb(:,:)
  COMPLEX(real_8), ALLOCATABLE, TARGET :: eigrb_dipole(:,:)

  LOGICAL :: tfnl2
  INTEGER :: ldf1,ldf2
  ! ==================================================================
  INTEGER :: natx
  COMPLEX(real_8), ALLOCATABLE :: ei1(:,:)
  COMPLEX(real_8), ALLOCATABLE :: ei2(:,:)
  COMPLEX(real_8), ALLOCATABLE :: ei3(:,:)


  ! ==================================================================

END MODULE sfac
