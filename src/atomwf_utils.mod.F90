#include "cpmd_global.h"

MODULE atomwf_utils
  USE atrho_utils,                     ONLY: atrho,&
                                             give_scr_atrho
  USE atwf,                            ONLY: atwf_mod,&
                                             atwp,&
                                             catom,&
                                             catord,&
                                             xsmat,&
                                             xxmat,&
                                             zsmat,&
                                             zxmat
  USE copot_utils,                     ONLY: copot,&
                                             give_scr_copot
  USE csmat_utils,                     ONLY: csmat
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE fnlalloc_utils,                  ONLY: fnl_set,&
                                             fnlalloc,&
                                             fnldealloc
  USE func,                            ONLY: func1,&
                                             mfxcc_is_pade,&
                                             mfxcc_is_skipped
  USE gpu
  USE gs_disortho_utils,               ONLY: gs_disortho
  USE gsortho_utils,                   ONLY: gs_ortho,&
                                             gs_ortho_c
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8,&
                                             int_8
  USE kpclean_utils,                   ONLY: c_clean
  USE kpts,                            ONLY: tkpts
  USE ksmat_dist_utils,                ONLY: dist_ksmat
  USE ksmat_utils,                     ONLY: ksmat
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sum
  USE nlcc,                            ONLY: corel
  USE nlps,                            ONLY: imagp
  USE ovlap_utils,                     ONLY: ovlap
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE prmem_utils,                     ONLY: prmem
  USE pslo,                            ONLY: pslo_com
  USE randtowf_utils,                  ONLY: randtowf
  USE reshaper,                        ONLY: reshape_inplace
  USE rgs_utils,                       ONLY: rgs,&
                                             rgs_c
  USE rnlsm_utils,                     ONLY: rnlsm
  USE sfac,                            ONLY: fnl
  USE sphe,                            ONLY: tsphere
  USE spin,                            ONLY: clsd,&
                                             lspin2,&
                                             spin_mod
  USE summat_utils,                    ONLY: summat
  USE system,                          ONLY: cntl,&
                                             kpbeg,&
                                             ncpw,&
                                             nkpbl,&
                                             nkpt
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE utils,                           ONLY: setkwf
  USE vofrho_utils,                    ONLY: give_scr_vofrho,&
                                             vofrho
  USE zeroing_utils,                   ONLY: zeroing
#ifdef _USE_SCRATCHLIBRARY
  USE scratch_interface,               ONLY: request_scratch,&
                                             free_scratch
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: atomwf
  PUBLIC :: give_scr_atomwf
  PUBLIC :: dsygvx
  !public :: zsygvx
  !public :: primao

CONTAINS

  ! ==================================================================
  SUBROUTINE atomwf(c0,nstate,tau0,fion,rhoe,psi)
    ! ==--------------------------------------------------------------==
    ! == ATOMIC WAVEFUNCTIONS AS INITIAL GUESS                        ==
    ! == WORKS WITH K POINTS: WITH KPOINT=GAMMA, SAME EIGENVALUES AS  ==
    ! == NO KPOINT VERSION x EXP(I THETA).                            ==
    ! == MOREOVER, IF EIGENSTATES ARE DEGENERATED, IN THE SUBSPACE    ==
    ! == EIGENVECTORS ARE DIFFERENT.                                  ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: c0(:,:,:)
    INTEGER,INTENT(IN)                       :: nstate
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: tau0(:,:,:)
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: rhoe(:,:)
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: psi(:,:)

    CHARACTER(*), PARAMETER                  :: procedureN = 'atomwf'
    COMPLEX(real_8), PARAMETER               :: zone = (1._real_8, 0._real_8),&
                                                zzero = (0._real_8, 0._real_8)

#ifdef _USE_SCRATCHLIBRARY
    COMPLEX(real_8), POINTER __CONTIGUOUS    :: cx(:,:), zorkat(:), zmatat(:,:)
    REAL(real_8), POINTER __CONTIGUOUS       :: eivat(:), rscr(:), workat(:)
#else
    COMPLEX(real_8), ALLOCATABLE             :: cx(:,:), zorkat(:), zmatat(:,:)
    REAL(real_8), ALLOCATABLE                :: eivat(:), rscr(:), workat(:)
#endif
    INTEGER :: ierr, ikind, ikk, ikpt, iopt, isub, kbeg, kend, kinc, nkpoint, nleft, nn, &
         nsmat, il_xmatat
    INTEGER(int_8) :: il_eivar(1), il_workat(1), il_rscr(1), il_xxmat(2), il_xsmat(2), &
         il_catom(2), il_zxmat(2), il_zsmat(2), il_zorkat(1), il_cx(2), il_zmatat(2), &
         il_eivat(1), i
    LOGICAL                                  :: debug, tdmal_tmp, tlsd2, &
                                                tlse2, ttau2
    REAL(real_8)                             :: dummy
    REAL(real_8), POINTER                    :: xmatat(:,:)

    CALL phfac(tau0,force_update=.TRUE.)
    ! ==--------------------------------------------------------------==
    IF (nstate.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)
    ! ==--------------------------------------------------------------==
    tdmal_tmp = cntl%tdmal
    ! vw>>>
    ! vw using always cntl%tdmal = T:
    ! vw    1) still bugs with cntl%cdft = T
    ! vw    2) not yet implemented with TKPNT = T 
    ! TK    3) not yet implemented with USPP
    ! vw so we hack...
    IF (.NOT.tkpts%tkpnt.AND..NOT.cntl%cdft.AND..NOT.pslo_com%tivan) &
         cntl%tdmal = .TRUE.
    ! vw<<<
    IF (cntl%tdmal.AND.paral%io_parent) THEN
       WRITE(6,*)&
            'NOTE: ATOMIC GUESS USING DISTRIBUTED LINALG WITH LANCZOS'
    ENDIF
    ! ==--------------------------------------------------------------==
    ! XMATAT and ZMATAT are matrices NATTOTxNATTOT
    IF (cntl%tdmal) THEN
       il_xmatat=atwp%nattot
    ELSE
       il_xmatat=imagp*atwp%nattot*atwp%nattot
    ENDIF
    il_zmatat(1)=atwp%nattot
    il_zmatat(2)=atwp%nattot
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_zmatat,zmatat,procedureN//'_zmatat',ierr)
#else
    ALLOCATE(zmatat(il_zmatat(1),il_zmatat(2)), stat=ierr)
#endif
    IF (ierr /= 0) CALL stopgm('atomwf.F90','Error allocating ZMATAT',& 
         __LINE__,__FILE__)
    ! XMATAT is real(8) :: alias to ZMATAT (C*16)
    CALL reshape_inplace(zmatat, (/atwp%nattot, atwp%nattot/), xmatat)

    ! EIVAT is eigenvalues.
    il_eivat(1)=imagp*atwp%nattot
    ! WORKAT is an auxiliary array.
    il_workat(1)=3*atwp%nattot
    ! Assignation
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_eivat,eivat,procedureN//'_eivat',ierr)
#else
    ALLOCATE(eivat(il_eivat(1)), stat=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_workat,workat,procedureN//'_workat',ierr)
#else    
    ALLOCATE(workat(il_workat(1)), stat=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    IF (tkpts%tkpnt) THEN
       il_zorkat=il_workat
    ELSE
       il_zorkat=1
    END IF
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_zorkat,zorkat,procedureN//'_zorkat',ierr)
#else    
    ALLOCATE(zorkat(il_zorkat(1)), stat=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    il_cx(1)=nkpt%ngwk
    il_cx(2)=atwp%numaormax
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_cx,cx,procedureN//'_cx',ierr)
#else
    ALLOCATE(cx(il_cx(1),il_cx(2)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(cx)!,nkpt%ngwk*atwp%numaormax)
    ! ==--------------------------------------------------------------==
    il_catom(1)=nkpt%ngwk
    il_catom(2)=atwp%nattot
#ifdef _USE_SCRATCHLIBRARY
    CALL request_scratch(il_catom,catom,procedureN//'_catom',ierr)
#else
    ALLOCATE(catom(il_catom(1),il_catom(2)),STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    CALL zeroing(catom)
    IF (tkpts%tkpnt) THEN
       il_zxmat(1)=atwp%nattot
       il_zxmat(2)=atwp%nattot
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_zxmat,zxmat,procedureN//'_zxmat',ierr)
#else
       ALLOCATE(zxmat(il_zxmat(1),il_zxmat(2)),STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
       il_zsmat(1)=atwp%nattot
       il_zsmat(2)=atwp%nattot
#ifdef _USE_SCRATCHLIBRARY
       CALL request_scratch(il_zsmat,zsmat,procedureN//'_zsmat',ierr)
#else
       ALLOCATE(zsmat(il_zsmat(1),il_zsmat(2)),STAT=ierr)

#endif
       IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)
    ELSE
       IF (.NOT.cntl%tdmal) THEN
          il_xxmat(1)=atwp%nattot
          il_xxmat(2)=atwp%nattot
#ifdef _USE_SCRATCHLIBRARY
          CALL request_scratch(il_xxmat,xxmat,procedureN//'_xxmat',ierr)
#else
          ALLOCATE(xxmat(il_xxmat(1),il_xxmat(2)),STAT=ierr)
#endif
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          il_xsmat(1)=atwp%nattot
          il_xsmat(2)=atwp%nattot
#ifdef _USE_SCRATCHLIBRARY
          CALL request_scratch(il_xsmat,xsmat,procedureN//'_xsmat',ierr)
#else
          ALLOCATE(xsmat(il_xsmat(1),il_xsmat(2)),STAT=ierr)
#endif
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF
    IF (paral%parent) CALL prmem('    ATOMWF')
    ! ==--------------------------------------------------------------==
    ! Compute the potential and energies for XC for the core charges.
    ! Need to have the correct LSD flag.
    IF (corel%tinlc) CALL copot(rhoe,psi,.FALSE.)

    ! ==--------------------------------------------------------------==
    ! Backup Spin polarization flag
    tlsd2=cntl%tlsd
    cntl%tlsd=.FALSE.
    tlse2=lspin2%tlse
    lspin2%tlse=.FALSE.
    ! Backup Meta functional
    ttau2 = cntl%ttau
    cntl%ttau  = .FALSE.
    IF (ttau2) THEN
       func1%mfxcc = mfxcc_is_pade
       cntl%tgc = .FALSE.
    ENDIF
    debug=.FALSE.
    ! ==--------------------------------------------------------------==
    ! Charge density
    CALL atrho(rhoe,psi(:,1),nstate)
    ! Orthogonalize atomic wavefunctions
    IF (cntl%tdmal) THEN
       IF (tkpts%tkpnt) CALL stopgm('ATOMWF','DISTRIBUTED INITIALIZATION '//&
            'NOT AVAILABLE WITH K-POINTS',& 
            __LINE__,__FILE__)
       CALL gs_disortho(atwp%nattot, catom)
    ENDIF
    ! The other kpoints are the same atomic wf (CATOM).
    IF (tkpts%tkpnt) CALL setkwf(ncpw%ngw,atwp%nattot,catom)
    ! Local potential
    CALL vofrho(tau0,fion,rhoe,psi,.FALSE.,.FALSE.)

    ! SUM ETOT, EKIN, EHT, EPSEU, ENL and EXC.
    CALL mp_sum(ener_com%etot,parai%allgrp)
    CALL mp_sum(ener_com%ekin,parai%allgrp)
    CALL mp_sum(ener_com%epseu,parai%allgrp)
    CALL mp_sum(ener_com%enl,parai%allgrp)
    CALL mp_sum(ener_com%eht,parai%allgrp)
    CALL mp_sum(ener_com%ehep,parai%allgrp)

    IF (cntl%tdmal) THEN
       CALL dist_ksmat(cx,rhoe,psi(:,1),nstate,1,c0,&
            MIN(nstate,atwp%nattot),tlsd2,.FALSE.)

       ! >>>
       ! add random vector to C0 and ortho them when NATTOT < NSTATE
       IF (tlsd2) THEN
          nn=MAX(0,spin_mod%nsup-atwp%nattot)+MAX(0,spin_mod%nsdown-atwp%nattot)
          IF (paral%io_parent.AND.nn.GT.0) WRITE(6,'(A,I5,A)')&
               ' ATOMWF| WARNING! RANDOMIZATION OF ',&
               nn,' STATES'
          nleft=spin_mod%nsup-atwp%nattot
          IF (nleft.GT.0) THEN
             CALL randtowf(c0(:,atwp%nattot+1:,1),&
                  spin_mod%nsup-atwp%nattot,ikind,1)
             nsmat=MAX(nleft*atwp%nattot,nleft*nleft)
             CALL gs_ortho(c0(:,:,1),atwp%nattot,c0(:,atwp%nattot+1:,1),&
                  nleft)
          ENDIF
          nleft=spin_mod%nsdown-atwp%nattot
          IF (nleft.GT.0) THEN
             CALL randtowf(c0(:,spin_mod%nsup+atwp%nattot+1:,1),&
                  nstate-(spin_mod%nsup+atwp%nattot),1,1)
             nsmat=MAX(nleft*atwp%nattot,nleft*nleft)
             CALL gs_ortho(c0(:,spin_mod%nsup+1:,1),atwp%nattot,&
                  c0(:,spin_mod%nsup+atwp%nattot+1:,1),nleft)
          ENDIF
       ELSE
          nleft=nstate-atwp%nattot
          IF (paral%io_parent.AND.nleft.GT.0) WRITE(6,'(A,I5,A)')&
               ' ATOMWF| WARNING! RANDOMIZATION OF ',&
               nleft,' STATES'
          IF (nleft.GT.0) THEN
             CALL randtowf(c0(:,atwp%nattot+1:,1),nleft,1,1)
             nsmat=MAX(nleft*atwp%nattot,nleft*nleft)
             CALL gs_ortho(c0,atwp%nattot,c0(:,atwp%nattot+1:,1),&
                  nleft)
          ENDIF
       ENDIF
       ! <<<
       GOTO 100
    ENDIF
    IF (cntl%tipao)THEN
       IF (tlsd2)THEN
          IF (paral%io_parent)&
               WRITE(6,*)"INITIALISING WAVEFUNCTION WITH",&
               " PRIMITIVE PSEUDO ATOMIC ORBITALS"
          CALL primao(c0)
          il_rscr(1)=2*nkpt%ngwk*MAX(nstate,atwp%nattot)
#ifdef _USE_SCRATCHLIBRARY
          CALL request_scratch(il_rscr,rscr,procedureN//'_rscr',ierr)
#else
          ALLOCATE(rscr(il_rscr(1)),STAT=ierr)
#endif
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)
          CALL rgs(c0(:,:,1),spin_mod%nsup,rscr)
          CALL rgs(c0(:,spin_mod%nsup+1:,1),spin_mod%nsdown,rscr)
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_rscr,rscr,procedureN//'_rscr',ierr)
#else
          DEALLOCATE(rscr,STAT=ierr)
#endif
          IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
               __LINE__,__FILE__)

          GOTO 100
       ELSE
          IF (paral%io_parent)&
               WRITE(6,'(A,I5,A)')&
               ' ATOMWF| WARNING! PRIMITIVE AO METHOD ONLY WORKS ',&
               'WITH LSD'
       ENDIF
    ENDIF

    ! ==--------------------------------------------------------------==
    ! Loop over k points
    CALL inq_swap(kbeg,kend,kinc)
    DO ikpt=kbeg,kend,kinc
       nkpoint=nkpbl(ikpt)
       IF (tkpts%tkblock) CALL rkpt_swap(dummy,1,ikpt,&
            'HGKP HGKM TWNL EIGKR MASKGW')
       DO ikind=1,nkpoint
          ikk=kpbeg(ikpt)+ikind
          ! Kohn-Sham and Overlap matrix
          CALL ksmat(cx,rhoe,psi(:,1),nstate,ikind)
          IF (pslo_com%tivan) THEN
             CALL fnl_set('SAVE')
             CALL fnlalloc(atwp%nattot,.FALSE.,.FALSE.)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             update_first_to_gpu  =.FALSE.
             update_second_to_gpu =.FALSE.
             update_third_to_gpu  =.FALSE.
             update_result_to_host=.FALSE.
             !$omp target enter data map (alloc:catom,xsmat,xxmat)
             !$omp target update to(catom,xsmat,xxmat)
#endif
             CALL rnlsm(catom,atwp%nattot,ikpt,ikind,.FALSE.)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             comm_buffers_on_host =.FALSE.
#endif
             CALL csmat(xsmat,catom,atwp%nattot,ikind,full=.FALSE.,store_nonort=.FALSE.,&
                  only_parent=.TRUE.)
             CALL summat(xxmat,atwp%nattot,symmetrization=.false.,parent=.true.)
#if defined(_HAS_OMP_TARGET_OFFLOAD)
             update_first_to_gpu  =.TRUE.
             update_second_to_gpu =.TRUE.
             update_third_to_gpu  =.TRUE.
             update_result_to_host=.TRUE.
             comm_buffers_on_host =.TRUE.
             !$omp target update from(xxmat,xsmat)
             !$omp target exit data map (release:catom,xsmat,xxmat)
#endif

             CALL fnldealloc(.FALSE.,.FALSE.)
             CALL fnl_set('RECV')
          ELSE
             IF (tkpts%tkpnt) THEN
                CALL ovlap_h(atwp%nattot,zsmat,catom,catom)
                CALL sumhmat(zxmat,atwp%nattot)
                CALL sumhmat(zsmat,atwp%nattot)
             ELSE
                CALL ovlap(atwp%nattot,xsmat,catom,catom)
                CALL summat(xxmat,atwp%nattot)
                CALL summat(xsmat,atwp%nattot)
             ENDIF
          ENDIF
          ! XMATAT,EIVAT,WORKAT and ZORKAT are needed only here
          CALL zeroing(xmatat)!,il_xmatat)
          iopt=1
          IF (paral%io_parent) THEN
             IF (tkpts%tkpnt) THEN
                CALL zsygvx(iopt,zxmat,atwp%nattot,zsmat,atwp%nattot,eivat,&
                     zmatat,atwp%nattot,atwp%nattot,zorkat,3*atwp%nattot,workat)
             ELSE
                CALL dsygvx(iopt,xxmat,atwp%nattot,xsmat,atwp%nattot,eivat,&
                     xmatat,atwp%nattot,atwp%nattot,workat,3*atwp%nattot)
             ENDIF
          ENDIF
          CALL mp_bcast(xmatat,il_xmatat,parai%io_source,parai%cp_grp)! vwe need to send that to all
          IF (tkpts%tkpnt) THEN
             IF (tlsd2) THEN
                nn=MIN(spin_mod%nsup,atwp%nattot)
                IF (nn.GT.0) THEN
                   CALL zgemm('N','N',nkpt%ngwk,nn,atwp%nattot,zone,catom(1,1),&
                        nkpt%ngwk,zmatat,atwp%nattot,zzero,c0(:,:,ikind),nkpt%ngwk)
                ENDIF
                nn=MIN(spin_mod%nsdown,atwp%nattot)
                IF (nn.GT.0) THEN
                   CALL zgemm('N','N',nkpt%ngwk,nn,atwp%nattot,zone,catom(1,1),&
                        nkpt%ngwk,zmatat,atwp%nattot,zzero,c0(:,spin_mod%nsup+1:,ikind),nkpt%ngwk)
                ENDIF
             ELSE
                nn=MIN(nstate,atwp%nattot)
                IF (nn.GT.0) THEN
                   CALL zgemm('N','N',nkpt%ngwk,nn,atwp%nattot,zone,catom(1,1),&
                        nkpt%ngwk,zmatat,atwp%nattot,zzero,c0(:,:,ikind),nkpt%ngwk)
                ENDIF
             ENDIF
             ! Spherical cutoff (Clean and Reorthogonalisation).
             IF (tsphere) THEN
                IF (tlsd2) THEN
                   nn=MIN(spin_mod%nsup,atwp%nattot)
                   IF (nn.GT.0) THEN
                      CALL c_clean(c0(:,:,ikind),nn,ikind)
                      CALL rgs_c  (c0(:,:,ikind),nn,zmatat)
                   ENDIF
                   nn=MIN(spin_mod%nsdown,atwp%nattot)
                   IF (nn.GT.0) THEN
                      CALL c_clean(c0(:,spin_mod%nsup+1:,ikind),nn,ikind)
                      CALL rgs_c  (c0(:,spin_mod%nsup+1:,ikind),nn,zmatat)
                   ENDIF
                ELSE
                   nn=MIN(nstate,atwp%nattot)
                   CALL c_clean(c0(:,:,ikind),nn,ikind)
                   CALL rgs_c  (c0(:,:,ikind),nn,zmatat)
                ENDIF
             ENDIF
          ELSE              ! IF (TKPNT) THEN
             IF (tlsd2) THEN
                nn=MIN(spin_mod%nsup,atwp%nattot)
                IF (nn.GT.0 .AND. nkpt%ngwk.GT.0)&
                     CALL dgemm('N','N',2*nkpt%ngwk,nn,atwp%nattot,1.0_real_8,catom(1,1),&
                     2*nkpt%ngwk,xmatat,atwp%nattot,0.0_real_8,c0(:,:, ikind),2*nkpt%ngwk)
                nn=MIN(spin_mod%nsdown,atwp%nattot)
                IF (nn.GT.0 .AND. nkpt%ngwk.GT.0)&
                     CALL dgemm('N','N',2*nkpt%ngwk,nn,atwp%nattot,1.0_real_8,catom(1,1),&
                     2*nkpt%ngwk,xmatat,atwp%nattot,0.0_real_8,c0(:,spin_mod%nsup+ 1:,ikind),&
                     2*nkpt%ngwk)
             ELSE
                nn=MIN(nstate,atwp%nattot)
                IF (nn.GT.0 .AND. nkpt%ngwk.GT.0)&
                     CALL dgemm('N','N',2*nkpt%ngwk,nn,atwp%nattot,1.0_real_8,catom(1,1),&
                     2*nkpt%ngwk,xmatat,atwp%nattot,0.0_real_8,c0(:,:, ikind),2*nkpt%ngwk)
             ENDIF
          ENDIF
          ! Randomize other states and orthogonalize.
          IF ((.NOT.tlsd2 .AND. nstate.GT.atwp%nattot) .OR.&
               (tlsd2 .AND.(spin_mod%nsup.GT.atwp%nattot .OR. spin_mod%nsdown.GT.atwp%nattot) ) ) THEN
             IF (tlsd2) THEN
                nn=MAX(0,spin_mod%nsup-atwp%nattot)+MAX(0,spin_mod%nsdown-atwp%nattot)
                IF (paral%io_parent.AND.ikk.EQ.1)&
                     WRITE(6,'(A,I5,A)')&
                     ' ATOMWF| WARNING! RANDOMIZATION OF ',&
                     nn,' STATES'
                IF (spin_mod%nsup-atwp%nattot .GT. 0) THEN
                   CALL randtowf(c0(:,atwp%nattot+1:,ikind),spin_mod%nsup-atwp%nattot,ikind,&
                        ikk)
                ENDIF
                IF (spin_mod%nsdown-atwp%nattot .GT. 0) THEN
                   CALL randtowf(c0(:,spin_mod%nsup+atwp%nattot+1:,ikind),nstate-(spin_mod%nsup+&
                        atwp%nattot),ikind,ikk)
                ENDIF
             ELSE
                IF (paral%io_parent.AND.ikk.EQ.1)&
                     WRITE(6,'(A,I5,A)')&
                     ' ATOMWF| WARNING! RANDOMIZATION OF ',&
                     nstate-atwp%nattot,' STATES'
                IF (nstate-atwp%nattot .GT. 0) THEN
                   CALL randtowf(c0(:,atwp%nattot+1:,ikind),nstate-atwp%nattot,ikind,&
                        ikk)
                ENDIF
             ENDIF
             IF (tkpts%tkpnt) THEN
                CALL c_clean(c0(:,:,ikind),atwp%nattot,ikind)
                IF (tlsd2) THEN
                   nleft=MAX(spin_mod%nsup-atwp%nattot,spin_mod%nsdown-atwp%nattot)
                   nsmat=MAX(2*nleft*atwp%nattot,2*nleft*nleft)

                   nleft=spin_mod%nsup-atwp%nattot
                   IF (nleft.GT.0)&
                        CALL gs_ortho_c(c0(:,:,ikind),atwp%nattot,&
                        c0(:, atwp%nattot+1,ikind),nleft,zmatat)
                   nleft=spin_mod%nsdown-atwp%nattot
                   IF (nleft.GT.0)&
                        CALL gs_ortho_c(c0(:,spin_mod%nsup+1:,ikind),atwp%nattot,&
                        c0(:,spin_mod%nsup+atwp%nattot+1:,ikind),nleft,zmatat)

                ELSE
                   nleft=nstate-atwp%nattot
                   nsmat=MAX(2*nleft*atwp%nattot,2*nleft*nleft)

                   IF (nleft.GT.0)&
                        CALL gs_ortho_c(c0(:,:,ikind),atwp%nattot,&
                        c0(:,atwp%nattot+1:,ikind),nleft,zmatat)

                ENDIF
             ELSE
                IF (tlsd2) THEN
                   nleft=MAX(spin_mod%nsup-atwp%nattot,spin_mod%nsdown-atwp%nattot)
                   nsmat=MAX(nleft*atwp%nattot,nleft*nleft)

                   nleft=spin_mod%nsup-atwp%nattot
                   IF (nleft.GT.0)&
                        CALL gs_ortho(c0(:,:,ikind),atwp%nattot,&
                        c0(:,atwp%nattot+1:,ikind),nleft)
                   nleft=spin_mod%nsdown-atwp%nattot
                   IF (nleft.GT.0)&
                        CALL gs_ortho(c0(:,spin_mod%nsup+1:,ikind),atwp%nattot,&
                        c0(:,spin_mod%nsup+atwp%nattot+1:,ikind),nleft)

                ELSE
                   nleft=nstate-atwp%nattot
                   nsmat=MAX(nleft*atwp%nattot,nleft*nleft)

                   IF (nleft.GT.0)&
                        CALL gs_ortho(c0(:,:,ikind),atwp%nattot,&
                        c0(:,atwp%nattot+1:,ikind),nleft)
                ENDIF
             ENDIF         ! IF (TKPNT) THEN

          ENDIF
       ENDDO
       IF (tkpts%tkblock) THEN
          CALL wkpt_swap(c0,nstate,ikpt,'C0')
          ! We calculate only for a k point.
          IF (tkpts%tknoswap) GOTO 100
       ENDIF
    ENDDO
100 CONTINUE

    IF (tkpts%tkpnt) THEN
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_zsmat,zsmat,procedureN//'_zsmat',ierr)
#else
       DEALLOCATE(zsmat,STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
       CALL free_scratch(il_zxmat,zxmat,procedureN//'_zxmat',ierr)
#else
       DEALLOCATE(zxmat,STAT=ierr)
#endif
       IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
            __LINE__,__FILE__)
    ELSE
       IF (.NOT.cntl%tdmal) THEN
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_xsmat,xsmat,procedureN//'_xsmat',ierr)
#else
          DEALLOCATE(xsmat,STAT=ierr)
#endif
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
          CALL free_scratch(il_xxmat,xxmat,procedureN//'_xxmat',ierr)
#else
          DEALLOCATE(xxmat,STAT=ierr)
#endif
          IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
               __LINE__,__FILE__)
       ENDIF
    ENDIF

#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_catom,catom,procedureN//'_catom',ierr)
#else
    DEALLOCATE(catom,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_cx,cx,procedureN//'_cx',ierr)
#else
    DEALLOCATE(cx,STAT=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_zorkat,zorkat,procedureN//'_zorkat',ierr)
#else
    DEALLOCATE(zorkat, stat=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_workat,workat,procedureN//'_workat',ierr)
#else
    DEALLOCATE(workat, stat=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_eivat,eivat,procedureN//'_eivat',ierr)
#else
    DEALLOCATE(eivat, stat=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    
#ifdef _USE_SCRATCHLIBRARY
    CALL free_scratch(il_zmatat,zmatat,procedureN//'_zmatat',ierr)
#else
    DEALLOCATE(zmatat, stat=ierr)
#endif
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    ! Spin polarization flag
    cntl%tlsd=tlsd2
    lspin2%tlse=tlse2
    cntl%ttau = ttau2
    IF (cntl%ttau) THEN
       func1%mfxcc = mfxcc_is_skipped
       cntl%tgc = .TRUE.
    ENDIF
    cntl%tdmal = tdmal_tmp
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE atomwf
  ! ==================================================================
  SUBROUTINE give_scr_atomwf(latomwf,tag,nstate)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: latomwf
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate

    INTEGER                                  :: latrho, lcopot, &
                                                lvofrho, &
                                                ngso, nleft
    LOGICAL                                  :: tlsd2, tlse2

    tlsd2=cntl%tlsd
    cntl%tlsd=.FALSE.
    tlse2=lspin2%tlse
    lspin2%tlse=.FALSE.
    lcopot=0
    ! 
    CALL give_scr_atrho(latrho,tag)
    ! COPOT
    IF (corel%tinlc) THEN
       CALL give_scr_copot(lcopot,tag)
    ENDIF
    ! VOFRHO
    CALL give_scr_vofrho(lvofrho,tag)
    IF (cntl%tdmal) THEN
       latomwf=imagp*atwp%nattot +           & ! OBSOLETE: XMATAT or ZMATAT
            imagp*atwp%nattot+                    & ! EIVAT
            imagp*3*atwp%nattot+10                ! WORKAT
    ELSE
       latomwf=imagp*atwp%nattot*atwp%nattot+         & ! XMATAT or ZMATAT
            imagp*atwp%nattot+                    & ! EIVAT
            imagp*3*atwp%nattot+10                ! WORKAT
    ENDIF
    IF (tkpts%tkpnt) latomwf=latomwf+3*atwp%nattot ! ZORKAT
    IF (nstate.GT.clsd%nlsd*atwp%nattot) THEN
       nleft=nstate-clsd%nlsd*atwp%nattot
       ngso=MAX(nleft**2,nleft*atwp%nattot)+nleft
       latomwf=MAX(latomwf,ngso)
    ENDIF
    latomwf=MAX(latrho,latomwf,lcopot,lvofrho)
    ! 
    cntl%tlsd=tlsd2
    lspin2%tlse=tlse2
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_atomwf
  ! ==================================================================
  SUBROUTINE dsygvx(iopt,a,lda,b,ldb,w,z,ldz,n,aux,naux)
    ! ==--------------------------------------------------------------==
    ! == GENERALIZED DIAGONALIZATION ROUTINE:                         ==
    ! == FOLLOW THE ESSL CONVENTION                                   ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iopt, lda, ldb, ldz, n
    REAL(real_8)                             :: z(ldz,n), w(n), b(ldb,n), &
                                                a(lda,n)
    INTEGER                                  :: naux
    REAL(real_8)                             :: aux(naux)

    INTEGER                                  :: i, info

    IF (iopt.EQ.1) THEN
       info=0
       CALL dsygv(1,'V','U',n,a,lda,b,ldb,w,aux,naux,info)
       DO i=1,n
          CALL dcopy(n,a(1,i),1,z(1,i),1)
       ENDDO
    ELSE
       CALL stopgm('DSYGVX','MISSING LIBRARY OPTION',& 
            __LINE__,__FILE__)
    ENDIF
    IF (info.NE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/," DSYGVX| INFO=",I5)') info
       IF (info.LT.0) THEN
          IF (paral%io_parent)&
               WRITE(6,'(" DSYGVX| THE ",I5,A)')&
               info,'TH ARGUMENT HAS AN ILLEGAL VALUE'
       ELSEIF (info.LE.n) THEN
          IF (paral%io_parent)&
               WRITE(6,&
               '(" DSYGVX| FAILED TO CONVERGE WHEN DIAGONALIZING A")')
       ELSEIF (info.GT.n) THEN
          IF (paral%io_parent)&
               WRITE(6,'(" DSYGVX| PB WITH OVERLAP MATRIX FROM CATOM")')
          IF (paral%io_parent)&
               WRITE(6,'(/,"-> CHECK IF YOUR CELL IS NOT TOO SMALL! !!")')
       ENDIF
       CALL stopgm ('DSYGVX','FAILED TO DIAGONALIZE',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dsygvx
  ! ==================================================================
  SUBROUTINE zsygvx(iopt,a,lda,b,ldb,w,z,ldz,n,aux,naux,rwork)
    ! ==--------------------------------------------------------------==
    ! == GENERALIZED DIAGONALIZATION ROUTINE:                         ==
    ! == FOLLOW DSYGVX CONVENTION                                     ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: iopt, lda, ldb, ldz, n
    COMPLEX(real_8)                          :: z(ldz,n)
    REAL(real_8)                             :: w(n)
    COMPLEX(real_8)                          :: b(ldb,n), a(lda,n)
    INTEGER                                  :: naux
    COMPLEX(real_8)                          :: aux(naux)
    REAL(real_8)                             :: rwork(naux)

    INTEGER                                  :: i, info

    info=0
    CALL zhegv(iopt,'V','U',n,a,lda,b,ldb,w,aux,naux,rwork,info)
    DO i=1,n
#ifdef __NEC
       CALL zcopy(n,a(1,i),1,z(1,i),1)
#else
       CALL dcopy(2*n,a(1,i),1,z(1,i),1)
#endif
    ENDDO
    IF (info.NE.0) THEN
       IF (paral%io_parent)&
            WRITE(6,'(/," ZSYGVX| INFO=",I5)') info
       IF (info.LT.0) THEN
          IF (paral%io_parent)&
               WRITE(6,'(" ZSYGVX| THE ",I5,A)')&
               info,'TH ARGUMENT HAS AN ILLEGAL VALUE'
       ELSEIF (info.LE.n) THEN
          IF (paral%io_parent)&
               WRITE(6,&
               '(" ZSYGVX| FAILED TO CONVERGE WHEN DIAGONALIZING A")')
       ELSEIF (info.GT.n) THEN
          IF (paral%io_parent)&
               WRITE(6,'(" ZSYGVX| PB WITH OVERLAP MATRIX FROM CATOM")')
          IF (paral%io_parent)&
               WRITE(6,'(/,"-> CHECK IF YOUR CELL IS NOT TOO SMALL! !!")')
       ENDIF
       CALL stopgm ('ZSYGVX','FAILED TO DIAGONALIZE',& 
            __LINE__,__FILE__)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE zsygvx
  ! ==================================================================
  ! ==================================================================
  SUBROUTINE primao(c0)
    ! ==--------------------------------------------------------------==
    ! == Initialise wavefunctions with atomic orbitals                ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8)                          :: c0(nkpt%ngwk,crge%n,*)

    INTEGER                                  :: ia, iat, idown, is, ish, iup, &
                                                ls, num, o, ocu

! TODO FIGURE OUT LS

    iup=1
    iat=1
    idown=spin_mod%nsup+1
    DO is=1,ions1%nsp
       DO ia=1,ions0%na(is)
          DO ish=1,atwf_mod%nshell(is)
             ocu=atwf_mod%oc(ish,is)
             ls=catord(ish,iat)
             num=MODULO(ocu,2)
             DO o=1,ocu-num,2
                CALL dcopy(2*nkpt%ngwk,catom(1,ls),1,c0(1,iup,1),1)
                iup=iup+1
                CALL dcopy(2*nkpt%ngwk,catom(1,ls),1,c0(1,idown,1),1)
                idown=idown+1
                ls=ls+1
             ENDDO
             IF (num.EQ.1)THEN
                CALL zcopy(nkpt%ngwk,catom(1,ls),1,c0(1,idown,1),1)
                idown=idown+1
             ENDIF
          ENDDO! Shells
          iat=iat+1
       ENDDO! Atoms
    ENDDO ! Species

    IF (spin_mod%nsup.NE.iup-1)THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,I5,A,I5,A)') "PRIMAO| GOT",iup-1,&
            "INSTEAD OF",spin_mod%nsup," UP ELECTRONS"
       CALL stopgm( 'PRIMAO', 'ILLEGAL RESULTS',& 
            __LINE__,__FILE__)
    ELSEIF (crge%n.NE.idown-1)THEN
       IF (paral%io_parent)&
            WRITE(6,'(A,I5,A,I5,A)') "PRIMAO| GOT",idown-1-spin_mod%nsup,&
            "INSTEAD OF",spin_mod%nsdown," DOWN ELECTRONS"
       CALL stopgm( 'PRIMAO', 'ILLEGAL RESULTS',& 
            __LINE__,__FILE__)
    ENDIF

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE primao
  ! ==================================================================

END MODULE atomwf_utils
