#include "cpmd_global.h"

MODULE hnlmat_utils
  USE cp_grp_utils,                    ONLY: cp_grp_split_atoms
  USE cvan,                            ONLY: deeq,&
                                             dvan
  USE distribution_utils,              ONLY: dist_atoms
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE gpu
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE nlps,                            ONLY: imagp,&
                                             nghtol,&
                                             nlps_com,&
                                             wsg
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE sfac,                            ONLY: fnl,&
                                             fnl_packed
  USE sgpp,                            ONLY: sgpp1,&
                                             sgpp2
  USE spin,                            ONLY: spin_mod
  USE system,                          ONLY: cntl,&
                                             parap,&
                                             maxsys
  USE timer,                           ONLY: tihalt,&
                                             tiset
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: hnlmat

CONTAINS

  ! ==================================================================
  SUBROUTINE hnlmat(hmat,f,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER,INTENT(IN)                       :: nstate
    REAL(real_8),INTENT(IN)                  :: f(nstate)
    REAL(real_8),INTENT(INOUT),TARGET        :: hmat(nstate,nstate)
    INTEGER                                  :: i, j, ia, is, isa0, ispin, nspin, isa_start,&
                                                ia_sum, tot_work, ns(2), nmin(2), off_i, &
                                                off_mat, off_fnl, ia_fnl, start_fnl, &
                                                end_fnl, fnl_start, start_mat, end_mat, &
                                                isub, ierr, n
    INTEGER(int_8)                           :: il_fnlat(2), il_fnlatj(2), il_hmat_loc(2)
    LOGICAL                                  :: need_hmat_loc, non_uspp
    REAL(real_8)                             :: ffi, fac, sum , fractions(parai%nproc)
#ifdef _USE_SCRATCHLIBRARY
    REAL(real_8),POINTER __CONTIGUOUS        :: fnlat(:,:), fnlatj(:,:)
#else
    REAL(real_8),ALLOCATABLE                 :: fnlat(:,:), fnlatj(:,:)
#endif
    REAL(real_8),POINTER __CONTIGUOUS        :: hmat_loc(:,:)
    INTEGER,ALLOCATABLE,SAVE                 :: na_buff(:,:,:)
    INTEGER,ALLOCATABLE                      :: na(:,:), na_fnl(:,:),na_grp(:,:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hnlmat'
    ! ==--------------------------------------------------------------==
    ! == Compute the non-local contribution to the Hamilton matrix    ==
    ! == optimized version for uspp only, calls old hnlmat in case of ==
    ! == other pseudo potentials                                      ==
    ! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    IF (cntl%tfdist) CALL stopgm('HNLMAT','TFDIST NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    IF (imagp.EQ.2)CALL stopgm('HNLMAT','K-POINT NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    IF (pslo_com%tivan) THEN
       ! VANDERBILT PP
       ns(1)=nstate
       ns(2)=0
       nmin(1)=1
       nmin(2)=0
       nspin=1
       fac=-2.0_real_8
       IF (cntl%tlsd) THEN
          ns(1)=spin_mod%nsup
          ns(2)=spin_mod%nsdown
          nmin(1)=1
          nmin(2)=spin_mod%nsup+1
          nspin=2
          fac=-1.0_real_8
       ENDIF
       ALLOCATE(na(2,ions1%nsp),na_fnl(2,ions1%nsp),na_grp(2,ions1%nsp,0:parai%cp_nogrp-1)&
            , stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'allocation problem',&
            __LINE__,__FILE__)
       !get cp_grp atom distribution
       CALL cp_grp_split_atoms(na_grp)
       na_fnl(:,:)=na_grp(:,:,parai%cp_inter_me)
       IF (.NOT. ALLOCATED(na_buff))THEN
          !distribute atoms between procs
          ALLOCATE(na_buff(2,ions1%nsp,0:parai%nproc-1), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate na_buff',&
               __LINE__,__FILE__)
          fractions=1.0_real_8/REAL(parai%nproc,KIND=real_8)
          CALL dist_atoms(parai%nproc,fractions,na_fnl,na_buff,only_uspp=.TRUE.)
       END IF
       na(:,:)=na_buff(:,:,parai%mepos)
       !do we have equal occupation numbers?
       need_hmat_loc=.FALSE.
       ffi=f(1)
       DO i=1,nstate
          IF (ffi.NE.f(i)) THEN
             need_hmat_loc=.TRUE.
             EXIT
          ENDIF
       ENDDO
       IF(need_hmat_loc)THEN
          il_hmat_loc(1)=nstate
          il_hmat_loc(2)=nstate

          ALLOCATE(hmat_loc(il_hmat_loc(1),il_hmat_loc(2)), stat=ierr)
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate hmat_loc',&
               __LINE__,__FILE__)
          fac=-1.0_real_8
       ELSE
          hmat_loc=>hmat
       END IF

       !calculate local part of total work
       tot_work=0
       DO is=1,ions1%nsp
          IF(pslo_com%tvan(is))THEN
             tot_work=tot_work+(na(2,is)-na(1,is)+1)*nlps_com%ngh(is)
          END IF
       END DO

       IF(tot_work.GT.0)THEN
          il_fnlat(1)=tot_work
          il_fnlat(2)=MAXVAL(ns)

          il_fnlatj(1)=tot_work
          il_fnlatj(2)=MAXVAL(ns)

#ifdef _USE_SCRATCHLIBRARY
          CALL request_scratch(il_fnlat,fnlat,procedureN//'_fnlat',ierr)
#else
          ALLOCATE(fnlat(il_fnlat(1),il_fnlat(2)), stat=ierr)
#endif
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlat',&
               __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
          CALL request_scratch(il_fnlatj,fnlatj,procedureN//'_fnlatj',ierr)
#else
          ALLOCATE(fnlatj(il_fnlatj(1),il_fnlatj(2)), stat=ierr)
#endif
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot allocate fnlatj',&
               __LINE__,__FILE__)


          DO ispin=1,nspin
             off_i=0
             IF(ispin.EQ.2) off_i=spin_mod%nsup
             n=ns(ispin)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp target teams distribute &
#else
             !$omp parallel do &
#endif
             !$omp& private(i,isa0,off_mat,off_fnl,is,ia_fnl,ia_sum,&
             !$omp& isa_start,fnl_start)
             DO i=1,n
                isa0=0
                off_mat=1
                off_fnl=1
                DO is=1,ions1%nsp
                   ia_fnl=na_fnl(2,is)-na_fnl(1,is)+1
                   ia_sum=na(2,is)-na(1,is)+1
                   !nc alreday filtered out in dist_atoms
                   IF(ia_sum.GT.0)THEN
                      !starting index isa
                      isa_start=isa0+na(1,is)-1
                      !starting index fnl_packed
                      fnl_start=na(1,is)-na_fnl(1,is)
                      !DIR$ forceinline
                      CALL prepare_matrix(fnl_packed(off_fnl:,i+off_i),&
                           fnlat(off_mat:,i),&
                           fnlatj(off_mat:,i),deeq(:,:,:,ispin),dvan(:,:,is),&
                           nlps_com%ngh(is),&
                           ia_sum,ia_fnl,fnl_start,isa_start,maxsys%nhxs,ions1%nat)
                      off_mat=off_mat+ia_sum*nlps_com%ngh(is)
                   END IF
                   off_fnl=off_fnl+ia_fnl*nlps_com%ngh(is)
                   isa0=isa0+ions0%na(is)
                END DO
             END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             !$omp target update from(fnlat) if(update_first_to_gpu)
             !$omp target update from(fnlatj) if(update_second_to_gpu)
#endif
#ifdef _HAS_DGEMMT
             CALL cpmd_dgemmt('U','T','N',n,tot_work,fac,&
                  fnlat(1,1),tot_work,fnlatj(1,1),tot_work,1.0_real_8,&
                  hmat_loc(nmin(ispin),nmin(ispin)),nstate)
#else
             CALL cpmd_dgemm('T','N',n,n,tot_work,fac,&
                  fnlat(1,1),tot_work,fnlatj(1,1),tot_work,1.0_real_8,&
                  hmat_loc(nmin(ispin),nmin(ispin)),nstate)
#endif
          END DO
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_fnlatj,fnlatj,procedureN//'_fnlatj',ierr)
#else
          DEALLOCATE(fnlatj, stat=ierr)
#endif
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnlatj',&
               __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_fnlat,fnlat,procedureN//'_fnlat',ierr)
#else
          DEALLOCATE(fnlat, stat=ierr)
#endif
          IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate fnlat',&
               __LINE__,__FILE__)
          IF(need_hmat_loc)THEN
             !$omp parallel do private(i,j)
             DO i=1,nstate
                DO j=1,nstate
                   IF (f(i) .GE. 1.e-5_real_8) THEN
                      hmat(j,i)=hmat(j,i)+hmat_loc(j,i)*f(i)
                   ELSE
                      hmat(j,i)=hmat(j,i)+hmat_loc(j,i)
                   END IF
                END DO
             END DO
             DEALLOCATE(hmat_loc, stat=ierr)
             IF (ierr /= 0) CALL stopgm(procedureN, 'Cannot deallocate deeqt',&
                  __LINE__,__FILE__)
          END IF
       END IF
       DEALLOCATE(na,na_fnl,na_grp, stat=ierr)
       IF (ierr /= 0) CALL stopgm(procedureN, 'deallocation problem',&
            __LINE__,__FILE__)
    END IF

    !!END VANDERBILT OPTIMIZED!!
    !lets check if there is something else...

    non_uspp=.FALSE.
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is))CYCLE
       non_uspp=.TRUE.
    END DO

    IF (non_uspp) CALL hnlmat_old(hmat,f,nstate)
    CALL tihalt(procedureN,isub)

  END SUBROUTINE hnlmat
  ! ==================================================================
  !DIR$ ATTRIBUTES FORCEINLINE::prepare_matrix
#if defined(_HAS_OMP_TARGET_OFFLOAD)
  SUBROUTINE prepare_matrix(fnl_p,fnli,fnlj,deeq_,dvan_,ngh,ia_sum,ia_fnl,fnl_start,&
       isa_start,maxngh,nat)
#else
  PURE  SUBROUTINE prepare_matrix(fnl_p,fnli,fnlj,deeq_,dvan_,ngh,ia_sum,ia_fnl,fnl_start,&
       isa_start,maxngh,nat)
#endif
    INTEGER,INTENT(IN)                       :: ngh,ia_sum,ia_fnl,isa_start,fnl_start,&
                                                maxngh,nat
    REAL(real_8),INTENT(IN)                  :: fnl_p(ia_fnl,ngh,*),dvan_(maxngh,*),&
                                                deeq_(nat,maxngh,*)
    REAL(real_8),INTENT(OUT)                 :: fnli(ia_sum,ngh,*),fnlj(ia_sum,ngh,*)
    INTEGER                                  :: iv,ia,jv,isa

#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp parallel do private(iv,ia)
#endif
    DO iv=1,ngh
       DO ia=1,ia_sum
          fnli(ia,iv,1)=fnl_p(ia+fnl_start,iv,1)
          fnlj(ia,iv,1)=0.0_real_8
       END DO
    END DO
#if defined(_HAS_OMP_TARGET_OFFLOAD)
    !$omp parallel do private(iv,jv,isa,ia)
#endif
    DO iv=1,ngh
       DO jv=1,ngh
          isa=isa_start
          DO ia=1,ia_sum
             isa=isa+1
             fnlj(ia,iv,1)=fnlj(ia,iv,1)&
                  +(dvan_(jv,iv)+deeq_(isa,jv,iv))*fnli(ia,jv,1)
          END DO
       END DO
    END DO
  END SUBROUTINE prepare_matrix
  ! ==================================================================
  SUBROUTINE hnlmat_old(hmat,f,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate), hmat(nstate,nstate)

    INTEGER                                  :: chunk, i, ia, ind, ind1, is, &
                                                isa, isa0, ispin, isub, iv, &
                                                j, jmax, jv, ki, kj, l, l2, &
                                                li, lj, mxmypos, mxrank, mypos,&
                                                na_grp(2,ions1%nsp,0:parai%cp_nogrp-1),&
                                                na_fnl(2,ions1%nsp)

    REAL(real_8)                             :: dd, fdd, ffi
    CHARACTER(*), PARAMETER                  :: procedureN = 'hnlmat_old'

! Variables
! ==--------------------------------------------------------------==
! == Compute the non-local contribution to the Hamilton matrix    ==
! ==--------------------------------------------------------------==

    CALL tiset(procedureN,isub)
    CALL cp_grp_split_atoms(na_grp)
    na_fnl(:,:)=na_grp(:,:,parai%cp_inter_me)
    IF (cntl%tfdist) CALL stopgm(procedureN,'TFDIST NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (imagp.EQ.2)CALL stopgm(procedureN,'K-POINT NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    DO i=1,nstate
       ffi=f(i)
       IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
       DO j=1,nstate
          hmat(i,j)=hmat(i,j)/ffi
       ENDDO
    ENDDO

    ! AK   FIXME: 20080131 the HNLMAT() optimization does not work for LSD.
    ! AK   FIXME: we drop back to old algorithm.

    IF (cntl%tlsd) THEN
       ! 1.. NSTATE for serial !
       DO i=parap%nst12(parai%mepos,1),parap%nst12(parai%mepos,2)
          ispin=1
          jmax=nstate
          IF (cntl%tlsd.AND.i.LE.spin_mod%nsup) jmax=spin_mod%nsup
          IF (cntl%tlsd.AND.i.GT.spin_mod%nsup) ispin=2
          DO j=i,jmax
             isa0=0
             DO is=1,ions1%nsp
                IF (sgpp1%tsgp(is)) THEN
                   ! Stefan Goedecker PP
                   DO iv=1,nlps_com%ngh(is)
                      l=nghtol(iv,is)+1
                      ki=sgpp2%lfval(iv,is)
                      li=sgpp2%lpval(iv,is)
                      DO jv=1,nlps_com%ngh(is)
                         l2=nghtol(jv,is)+1
                         lj=sgpp2%lpval(jv,is)
                         IF (l2.EQ.l.AND.li.EQ.lj) THEN
                            kj=sgpp2%lfval(jv,is)
                            DO ia=na_fnl(1,is),na_fnl(2,is)
                               isa=isa0+ia
                               fdd=sgpp2%hlsg(ki,kj,l,is)*fnl(1,isa,iv,i,1)*&
                                    fnl(1,isa,jv,j,1)
                               hmat(i,j)=hmat(i,j)-fdd
                               IF (i.NE.j) hmat(j,i)=hmat(j,i)-fdd
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                ELSE
                   ! BACHELET HAMANN SCHLUTER
                   DO iv=1,nlps_com%ngh(is)
                      dd=wsg(is,iv)
                      DO ia=na_fnl(1,is),na_fnl(2,is)
                         isa=isa0+ia
                         fdd=dd*fnl(1,isa,iv,i,1)*fnl(1,isa,iv,j,1)
                         hmat(i,j)=hmat(i,j)-fdd
                         IF (i.NE.j) hmat(j,i)=hmat(j,i)-fdd
                      ENDDO
                   ENDDO
                   ! ==-------------------------------------------------------==
                ENDIF
                isa0=isa0+ions0%na(is)
             ENDDO
          ENDDO
       ENDDO
    ELSE
       chunk = 0
       mxrank = 0
       IF (parap%nst12(parai%mepos,1).LE.parap%nst12(parai%mepos,2)) THEN
          ind1 = 0
          DO i=0, parai%nproc-1
             IF (parap%nst12(i,1).LE.parap%nst12(i,2)) THEN
                ind1 = ind1 + 1
                IF (i.GT.mxrank) mxrank=i
                IF (parai%mepos.EQ.i) mxmypos = ind1
             ENDIF
          ENDDO

          chunk = (nstate*(nstate+1)/2) / ind1
          IF (parai%mepos.EQ.mxrank) chunk=(nstate*(nstate+1)/2)-(ind1-1)*&
               chunk
          ind = 1
          mypos = (mxmypos-1)*chunk+1
          IF (parai%mepos.EQ.mxrank) THEN
             mypos=(ind1-1)*((nstate*(nstate+1)/2)/ind1)+1
          ENDIF
          DO i=1, nstate
             DO j=i, nstate
                IF (ind.EQ.mypos) GOTO 999
                ind = ind +1
             ENDDO
          ENDDO
999       CONTINUE
       ENDIF

       ! 1.. NSTATE for serial !
       isa0 = 0
       DO ind = 1, chunk

          ispin=1
          jmax=nstate
          IF (cntl%tlsd.AND.i.LE.spin_mod%nsup) jmax=spin_mod%nsup
          IF (cntl%tlsd.AND.i.GT.spin_mod%nsup) ispin=2
          DO is=1,ions1%nsp
             IF (sgpp1%tsgp(is)) THEN
                ! Stefan Goedecker PP
                DO iv=1,nlps_com%ngh(is)
                   l=nghtol(iv,is)+1
                   ki=sgpp2%lfval(iv,is)
                   li=sgpp2%lpval(iv,is)
                   DO jv=1,nlps_com%ngh(is)
                      l2=nghtol(jv,is)+1
                      lj=sgpp2%lpval(jv,is)
                      IF (l2.EQ.l.AND.li.EQ.lj) THEN
                         kj=sgpp2%lfval(jv,is)
                         DO ia=na_fnl(1,is),na_fnl(2,is)
                            isa=isa0+ia
                            fdd=sgpp2%hlsg(ki,kj,l,is)*fnl(1,isa,iv,i,1)*&
                                 fnl(1,isa,jv,j,1)
                            hmat(i,j)=hmat(i,j)-fdd
                            IF (i.NE.j) hmat(j,i)=hmat(j,i)-fdd
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
             ELSE
                ! BACHELET HAMANN SCHLUTER
                DO iv=1,nlps_com%ngh(is)
                   dd=wsg(is,iv)
                   DO ia=na_fnl(1,is),na_fnl(2,is)
                      isa=isa0+ia
                      fdd=dd*fnl(1,isa,iv,i,1)*fnl(1,isa,iv,j,1)
                      hmat(i,j)=hmat(i,j)-fdd
                      IF (i.NE.j) hmat(j,i)=hmat(j,i)-fdd
                   ENDDO
                ENDDO
                ! ==-------------------------------------------------------==
             ENDIF
             isa0=isa0+ions0%na(is)
          ENDDO
          j = j+1
          isa0 = 0
          IF (j.GT.nstate) THEN
             j = i+1
             i = i+1
             isa0 = 0
          ENDIF
       ENDDO
       ! AK: end of .NOT.cntl%tlsd branch. FIXME.
    ENDIF

    DO i=1,nstate
       ffi=f(i)
       IF (ffi.LT.1.e-5_real_8) ffi=1.0_real_8
       DO j=1,nstate
          hmat(i,j)=hmat(i,j)*ffi
       ENDDO
    ENDDO
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hnlmat_old
  ! ==================================================================

END MODULE hnlmat_utils
