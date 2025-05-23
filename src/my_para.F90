#include "cpmd_global.h"

SUBROUTINE mp_bcast_byte(DATA,n,root,comm)
  ! ==--------------------------------------------------------------==
  ! == Wrapper to mpi_bcast                                         ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_mpi_error_assert
  USE pstat
  USE machine, ONLY: m_walltime
#ifdef __PARALLEL
  USE mpi_f08
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER :: DATA(*)
  INTEGER :: n,root
#ifdef __PARALLEL
  type(MPI_COMM) :: comm
#else
  INTEGER :: comm
#endif
  ! Variables
  INTEGER :: ierr
  CHARACTER(*),PARAMETER :: procedureN='mp_bcast_byte'
#ifdef __PARALLEL
  INTEGER,SAVE :: size
  LOGICAL :: first=.true.
  REAL(real_8) :: tim1,tim2
  ! ==--------------------------------------------------------------==
  IF(first)THEN
    CALL mpi_type_size(mpi_byte,size,ierr)
    first=.FALSE.
  END IF
  cmcal(ipar_cast)=cmcal(ipar_cast)+1.0d0
  cmlen(ipar_cast)=cmlen(ipar_cast)+n*size
  tim1=m_walltime()
  CALL mpi_bcast(DATA,n,mpi_byte,root,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_cast)=cmtim(ipar_cast)+tim2-tim1
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE mp_bcast_byte


SUBROUTINE my_stopall(code)
  ! ==--------------------------------------------------------------==
  ! == STOP ALL PROCESSORS                                          ==
  ! ==--------------------------------------------------------------==
  USE mp_interface, ONLY: mp_comm_world
#ifdef __PARALLEL
  USE mpi_f08
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER :: code
#ifdef __PARALLEL
  ! Variables
  INTEGER :: ierr
  LOGICAL :: is_init
  ! ==--------------------------------------------------------------==
  CALL mpi_initialized(is_init,ierr)
  if( is_init ) then
     CALL mpi_abort(mp_comm_world,code,ierr)
  else
     stop 1
  endif
#endif
  STOP 1
  ! ==--------------------------------------------------------------==
END SUBROUTINE my_stopall

! ==================================================================
SUBROUTINE my_concat(outmsg,inmsg,blklen,gid)
  ! ==--------------------------------------------------------------==
  ! == Concat all OUTMSG from GID group_ processor into INMSG        ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_mpi_error_assert
  USE utils, ONLY : icopy
#ifdef __PARALLEL
  USE mpi_f08
  USE pstat
  USE machine, ONLY: m_walltime
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER :: outmsg(*),inmsg(*),blklen
#ifdef __PARALLEL
  type(MPI_COMM) :: gid
#else
  INTEGER :: gid
#endif
  CHARACTER(*),PARAMETER::procedureN='my_concat'
#ifdef __PARALLEL
  ! Variables
  INTEGER,SAVE :: size
  LOGICAL :: first=.true.
  INTEGER :: ierr,nnodes
  REAL(real_8) :: tim1,tim2
  ! ==--------------------------------------------------------------==
  IF(first)THEN
    CALL mpi_type_size(mpi_double_precision,size,ierr)
    first=.FALSE.
  END IF
  CALL mpi_comm_size ( gid, nnodes, ierr )
  cmcal(ipar_agav)=cmcal(ipar_agav)+1.0d0
  cmlen(ipar_agav)=cmlen(ipar_agav)+nnodes*size*blklen
  tim1=m_walltime()
  ! ==--------------------------------------------------------------==
  CALL mpi_allgather(outmsg,blklen,mpi_byte,inmsg,blklen,&
       mpi_byte,gid,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_agav)=cmtim(ipar_agav)+tim2-tim1
#else
  ! Variables
  INTEGER :: len
  ! ==--------------------------------------------------------------==
  len = blklen/(8/2)
  CALL icopy(len,outmsg,1,inmsg,1)
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE my_concat
! ==================================================================
SUBROUTINE my_concat_inplace(inmsg,blklen,gid)
  ! ==--------------------------------------------------------------==
  ! == Concat all OUTMSG from GID group_ processor into INMSG        ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_mpi_error_assert
  USE utils, ONLY : icopy
#ifdef __PARALLEL
  USE mpi_f08
  USE pstat
  USE machine, ONLY: m_walltime
#endif
  IMPLICIT NONE
  ! Arguments
  REAL(real_8) :: inmsg(*)
#ifdef __PARALLEL
  integer :: blklen
  type(MPI_COMM) :: gid
#else
  integer :: blklen,gid
#endif
  CHARACTER(*),PARAMETER::procedureN='my_concat'
#ifdef __PARALLEL
  ! Variables
  INTEGER :: ierr,nnodes
  INTEGER,SAVE :: size
  LOGICAL :: first=.true.
  REAL(real_8) :: tim1,tim2
  ! ==--------------------------------------------------------------==
  IF(first)THEN
    CALL mpi_type_size(mpi_double_precision,size,ierr)
    first=.FALSE.
  END IF
  CALL mpi_comm_size ( gid, nnodes, ierr )
  cmcal(ipar_agav)=cmcal(ipar_agav)+1.0d0
  cmlen(ipar_agav)=cmlen(ipar_agav)+nnodes*size*blklen
  tim1=m_walltime()
  CALL mpi_allgather(mpi_in_place,0,MPI_DATATYPE_NULL,inmsg,blklen,&
       mpi_double_precision,gid,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_agav)=cmtim(ipar_agav)+tim2-tim1
#else
  ! Nothing to do, in=out
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE my_concat_inplace
! ==----------------------------------------------------------------==
SUBROUTINE my_concatv(outmsg,inmsg,blklen,recvcnt,recvdispl,gid)
  ! ==--------------------------------------------------------------==
  ! == Concat all OUTMSG of different lengths from GID group_       ==
  ! == processor into INMSG. NOTE DOUBLE PRECISION ARRAYS USED      ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE gpu,                             ONLY: comm_buffers_on_host
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_mpi_error_assert
#ifdef __PARALLEL
  USE mpi_f08
  USE pstat
  USE machine, ONLY: m_walltime
#endif
  IMPLICIT NONE
  ! Arguments
  REAL(real_8) :: outmsg(*),inmsg(*)
  INTEGER :: blklen,recvcnt(*),recvdispl(*)
#ifdef __PARALLEL
  type(MPI_COMM) :: gid
#else
  INTEGER :: gid
#endif
  CHARACTER(*),PARAMETER::procedureN='my_concatv'
#ifdef __PARALLEL
  ! Variables
  INTEGER :: ierr,nnodes,sum,i
  INTEGER,SAVE :: size
  LOGICAL :: first=.true.
  REAL(real_8) :: tim1,tim2
  ! ==--------------------------------------------------------------==
  IF(first)THEN
    CALL mpi_type_size(mpi_double_precision,size,ierr)
    first=.FALSE.
  END IF
  CALL mpi_comm_size ( gid, nnodes, ierr )
  cmcal(ipar_agav)=cmcal(ipar_agav)+1.0d0
  sum=0
  DO i=1,nnodes
     sum=sum+recvcnt(i)
  END DO
  cmlen(ipar_agav)=cmlen(ipar_agav)+sum*size
  tim1=m_walltime()
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp target data use_device_addr(outmsg(1:sum),inmsg(1:sum)) if(.NOT.comm_buffers_on_host)
#endif
  CALL mpi_allgatherv(outmsg,blklen,mpi_double_precision,inmsg,&
       recvcnt,recvdispl,mpi_double_precision&
       ,gid,ierr)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  !$omp end target data
#endif
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_agav)=cmtim(ipar_agav)+tim2-tim1
#else
  ! ==--------------------------------------------------------------==
  CALL dcopy(blklen,outmsg,1,inmsg(1+recvdispl(1)),1)
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE my_concatv
  ! ==--------------------------------------------------------------==
  ! ==--------------------------------------------------------------==
SUBROUTINE my_concatv_inplace(inmsg,recvcnt,recvdispl,gid)
  ! ==--------------------------------------------------------------==
  ! == Concat all OUTMSG of different lengths from GID group_        ==
  ! == processor into INMSG. NOTE DOUBLE PRECISION ARRAYS USED      ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_mpi_error_assert
#ifdef __PARALLEL
  USE mpi_f08
  USE pstat
  USE machine, ONLY: m_walltime
#endif
  IMPLICIT NONE
  ! Arguments
  REAL(real_8) :: inmsg(*)
  INTEGER :: blklen,recvcnt(*),recvdispl(*)
#ifdef __PARALLEL
  type(MPI_COMM) :: gid
#else
  INTEGER :: gid
#endif
  CHARACTER(*),PARAMETER::procedureN='my_concatv'
#ifdef __PARALLEL
  ! Variables
  INTEGER :: ierr,nnodes,sum,i
  INTEGER,SAVE :: size
  LOGICAL :: first=.true.
  REAL(real_8) :: tim1,tim2
  ! ==--------------------------------------------------------------==
  IF(first)THEN
    CALL mpi_type_size(mpi_double_precision,size,ierr)
    first=.FALSE.
  END IF
  CALL mpi_comm_size ( gid, nnodes, ierr )
  cmcal(ipar_agav)=cmcal(ipar_agav)+1.0d0
  sum=0
  DO i=1,nnodes
     sum=sum+recvcnt(i)
  END DO
  cmlen(ipar_agav)=cmlen(ipar_agav)+sum*size
  tim1=m_walltime()
  CALL mpi_allgatherv(mpi_in_place,0,MPI_DATATYPE_NULL,inmsg,&
       recvcnt,recvdispl,mpi_double_precision&
       ,gid,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  tim2=m_walltime()
  cmtim(ipar_agav)=cmtim(ipar_agav)+tim2-tim1
#else
  !nothing to do, in=out
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE my_concatv_inplace
! ==================================================================
! ==================================================================
SUBROUTINE my_source_concatv(outmsg,inmsg,blklen,recvcnt,&
     recvdispl,root,gid)
  ! ==--------------------------------------------------------------==
  ! == Concat all OUTMSG of different lengths from GID group_        ==
  ! == processor into INMSG. NOTE DOUBLE PRECISION ARRAYS USED      ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_mpi_error_assert
#ifdef __PARALLEL
  USE mpi_f08
#endif
  IMPLICIT NONE
  ! Arguments
  REAL(real_8) :: outmsg(*),inmsg(*)
  INTEGER :: blklen,recvcnt(*),recvdispl(*),root
#ifdef __PARALLEL
  type(MPI_COMM) :: gid
#else
  INTEGER :: gid
#endif
  CHARACTER(*),PARAMETER::procedureN='my_source_concatv'
#ifdef __PARALLEL
  ! Variables
  INTEGER :: ierr
  INTEGER :: err_msg_len, ierr_err
  ! ==--------------------------------------------------------------==
  CALL mpi_gatherv(outmsg,blklen,mpi_double_precision,inmsg,&
       recvcnt,recvdispl,mpi_double_precision&
       ,root,gid,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  ! ==--------------------------------------------------------------==
  CALL dcopy(blklen,outmsg,1,inmsg(1+recvdispl(1)),1)
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE my_source_concatv
! ==================================================================
SUBROUTINE my_trans(outmsg,inmsg,blklen,group_)
  ! ==--------------------------------------------------------------==
  ! == Matrix Transpose: group_ 1 : all processors                   ==
  ! ==                   group_ 2 : orbital split                    ==
  ! ==                   group_ 3 : plane wave split                 ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE machine, ONLY: m_walltime
  USE mp_interface, ONLY: mp_mpi_error_assert
  USE system , ONLY:group
  USE parac, ONLY : paral,parai
  USE pstat , ONLY:cmcal,cmlen,cmtim,ipar_aall
  USE utils, ONLY : icopy
#ifdef __PARALLEL
  USE mpi_f08
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER :: outmsg(*),inmsg(*),blklen,group_
#ifdef __PARALLEL
  type(MPI_COMM) :: gid
#else
  INTEGER :: gid
#endif
  ! Variables
  INTEGER :: nn
  REAL(real_8) :: tim1,tim2
  CHARACTER(*),PARAMETER::procedureN='my_trans'
#ifdef __PARALLEL
  INTEGER :: ierr
#else
  INTEGER :: len
#endif
  ! ==--------------------------------------------------------------==
  cmcal(ipar_aall)=cmcal(ipar_aall)+1.0_real_8
  tim1=m_walltime()
  IF (group_.EQ.1) THEN
     gid=parai%allgrp
     nn=parai%nproc
  ELSEIF (group_.EQ.2) THEN
     gid=group%meogrp
     nn=group%nogrp
  ELSEIF (group_.EQ.3) THEN
     gid=group%mepgrp
     nn=group%npgrp
  ENDIF
#ifdef __PARALLEL
  CALL mpi_alltoall(outmsg,blklen,mpi_byte,inmsg,blklen,&
       mpi_byte,gid,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  len = blklen/(8/2)
  CALL icopy(len,outmsg(1),1,inmsg(1),1)
#endif
  tim2=m_walltime()
  cmtim(ipar_aall)=cmtim(ipar_aall)+tim2-tim1
  cmlen(ipar_aall)=cmlen(ipar_aall)+blklen*nn
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE my_trans
! ==================================================================
SUBROUTINE my_allgather_i(DATA,n,comm)
  ! ==--------------------------------------------------------------==
  ! == Wrapper to mpi_allgather                                     ==
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_mpi_error_assert
#ifdef __PARALLEL
  USE mpi_f08
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER :: DATA(*)
  INTEGER :: n
#ifdef __PARALLEL
  type(MPI_COMM) :: comm
#else
  INTEGER :: comm
#endif
  ! Variables
  CHARACTER(*),PARAMETER :: procedureN='MY_ALLGATHER_I'
#ifdef __PARALLEL
  INTEGER :: ierr
  ! ==--------------------------------------------------------------==
  CALL mpi_allgather(mpi_in_place,n,mpi_integer,DATA,n,&
       mpi_integer,comm,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE my_allgather_i
! ==================================================================
SUBROUTINE my_shift(msend,mrecv,msglen,mep,ip,gid)
  ! ==--------------------------------------------------------------==
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  USE timer, ONLY: tiset, tihalt
  USE mp_interface, ONLY: mp_mpi_error_assert
  USE parac, ONLY : paral,parai
  USE pstat , ONLY:cmcal,cmlen,cmtim,ipar_aall, ipar_send
#ifdef __PARALLEL
  USE mpi_f08
#endif
  IMPLICIT NONE
  ! Arguments
  INTEGER :: msglen,mep,ip
#ifdef __PARALLEL
  type(MPI_COMM) :: gid
#else
  INTEGER :: gid
#endif
  REAL(real_8) :: msend(*),mrecv(*)
  CHARACTER(*),PARAMETER::procedureN='my_shift'
#ifdef __PARALLEL
  ! Variables
  INTEGER :: ipsend,iprecv,&
       whoami,howmany,&
       itype,ierr
  type(MPI_Request) :: irequest
  type(MPI_Status) :: status
  ! ==--------------------------------------------------------------==
  CALL mpi_comm_rank(gid,whoami,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  CALL mpi_comm_size(gid,howmany,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  ipsend=MOD(whoami+ip+howmany,howmany)
  iprecv=MOD(whoami-ip+howmany,howmany)

  cmcal(ipar_send)=cmcal(ipar_send)+1.0_real_8
  cmlen(ipar_send)=cmlen(ipar_send)+msglen
  itype=1
  CALL mpi_isend(msend,msglen,mpi_byte,ipsend,itype,&
       gid,irequest,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  CALL mpi_recv(mrecv,msglen,mpi_byte,iprecv,itype,&
       gid,status,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
  CALL mpi_wait(irequest,status,ierr)
  CALL mp_mpi_error_assert(ierr,procedureN,__LINE__,__FILE__)
#else
  CALL stopgm('MY_SHIFT','NOT INSTALLED',& 
       __LINE__,__FILE__)
#endif
  ! ==--------------------------------------------------------------==
  RETURN
END SUBROUTINE my_shift
