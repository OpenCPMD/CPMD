#include "cpmd_global.h"

#define PRINT_GROUP_INFOS .FALSE.

MODULE hfx_utils
  USE cnst,                            ONLY: uimag
  USE cppt,                            ONLY: scgx
  USE dotp_utils,                      ONLY: dotp
  USE error_handling,                  ONLY: stopgm
  USE fft,                             ONLY: inzf,&
                                             inzs,&
                                             jgw,&
                                             jhg,&
                                             llr1,&
                                             nzff,&
                                             nzfs
  USE fft_maxfft,                      ONLY: maxfftn
  USE fftmain_utils,                   ONLY: fwfftn,&
                                             invfftn
  USE fftnew_utils,                    ONLY: setfftn
  USE fnonloc_utils,                   ONLY: fnonloc_hfx
  USE func,                            ONLY: func1,&
                                             func3
  USE geq0mod,                         ONLY: geq0
  USE hfxmod,                          ONLY: hfxc3,&
                                             hfxc4,&
                                             ipoolhfx,&
                                             wcentx
  USE ions,                            ONLY: ions1
  USE kinds,                           ONLY: int_8,&
                                             real_8
  USE kpts,                            ONLY: tkpts
  USE machine,                         ONLY: m_cputime,&
                                             m_flush
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_max,&
                                             mp_sum
  USE newd_utils,                      ONLY: newd
  USE parac,                           ONLY: parai,&
                                             paral
  USE part_1d,                         ONLY: part_1d_elem_to_proc,&
                                             part_1d_get_elem,&
                                             part_1d_nbr_elems,&
                                             part_1d_symm_holds_pair,&
                                             proc_to_grid2d
  USE pbc_utils,                       ONLY: pbc
  USE pslo,                            ONLY: pslo_com
  USE rhov_utils,                      ONLY: rhov
  USE ropt,                            ONLY: infi,&
                                             infw,&
                                             iteropt
  USE rswfmod,                         ONLY: maxstatesx,&
                                             rswfx
  USE spin,                            ONLY: lspin2,&
                                             spin_mod
  USE state_utils,                     ONLY: &
       add_wfn, copy_im_to_im, copy_im_to_re, copy_re_to_im, copy_re_to_re, &
       set_psi_1_state_g, set_psi_2_states_g, set_psi_im_state_r, &
       set_psi_re_state_r, zero_wfn
  USE system,                          ONLY: cntl,&
                                             group,&
                                             ncpw,&
                                             parm, &
                                             spar !SM
  USE timer,                           ONLY: tihalt,&
                                             tiset
  USE zeroing_utils,                   ONLY: zeroing
  !
  USE ace_hfx  !SAGAR
  USE ovlap_utils,                     ONLY: ovlap
  !
  USE cp_grp_utils,                    ONLY: cp_grp_redist !sagar

  IMPLICIT NONE

  PRIVATE


  PUBLIC :: hfx_old
  !public :: check_occupation
  PUBLIC :: get_wannier_separation
  !public :: hfxab
  !public :: hfxab2
  !public :: hfxaa
  PUBLIC :: hfxpsi_old
  !public :: hfxpb
  PUBLIC :: hfxrpa_old
  !public :: hfxrpav
  PUBLIC :: hfx_scdm !SM
  PUBLIC :: hfx_ace !SM

CONTAINS


!===================================================================!
! SM ACE
!================================================================== 
  SUBROUTINE hfx_ace(c0,c2,f,psic,nstate,ehfx,vhfx)
!==--------------------------------------------------------------==
! Apply ACE projectors to calculate EHFX                            
! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:) 
    REAL(real_8)                             :: f(:)
    COMPLEX(real_8)                          :: psic(:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: ehfx, vhfx
!====================================================================    
    CHARACTER(*), PARAMETER                  :: procedureN = 'hfx_ace'
    COMPLEX(real_8), PARAMETER               :: zone = (1.0_real_8,0.0_real_8)

    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C2_ace, cmexx
    
    INTEGER                                  :: i, ierr, isub
!=======================================================================
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: mexx
!=======================================================================
    REAL(real_8)                             :: ehfx_ace
    INTEGER                                  :: info
! ==------------------------------------------------------------------==

    ehfx = 0.0_real_8
    vhfx = 0.0_real_8

    IF (func1%mhfx /= 1) RETURN

    CALL tiset(procedureN,isub)
!-------------------------------------------------------------------
    ALLOCATE(cmexx(nstate,nstate),mexx(nstate,nstate),&
         C2_ace(ncpw%ngw,nstate),&
         stat=ierr)
    IF (ierr /= 0) CALL stopgm( procedureN, 'Allocation problem' ,&
         __LINE__,__FILE__)

!-------------------------------------------------------------------
!   Apply ACE projectors to calculate EHFX 
! ==================================================================
        ! COMPUTES THE OVERLAP MATRIX A = < C1 | C2 > 
        ! ! <xi|phi> ACE potential
        !CALL AZZERO(mexx,NSTATE*NSTATE)
        CALL zeroing(mexx)
        !
        call OVLAP(NSTATE,mexx,xi,C0)
        !
        !CALL MY_SUM_D(mexx,NSTATE*NSTATE,ALLGRP)
        CALL mp_sum(mexx,nstate*nstate,parai%allgrp)
        !
        ! |vv> = |vphi> + (-One) * |xi> * <xi|phi>
        cmexx = (1.d0,0.d0)*mexx
        !
        !CALL ZAZZERO( C2_ace, NGW * NSTATE )
        CALL zeroing(C2_ace)
        !
        CALL ZGEMM ('N','N',ncpw%ngw,nstate,nstate,(1.d0,0.d0),xi, &
                   ncpw%ngw,cmexx,nstate,(1.d0,0.d0),c2_ace,ncpw%ngw)
        !
        !CALL AZZERO(mexx,NSTATE*NSTATE)
        CALL zeroing(mexx)
        !
        call OVLAP(NSTATE,mexx,C0,c2_ace)
        !
        !CALL MY_SUM_D(mexx,NSTATE*NSTATE,ALLGRP)
        CALL mp_sum(mexx,nstate*nstate,parai%allgrp)
        !
        ehfx_ace=0.0d0
        !
        do i=1,nstate
         ehfx_ace=ehfx_ace+mexx(i,i)
        end do
        ehfx_ace=-ehfx_ace*0.5d0
        !
  !      DO IA=1,NSTATE
  !         VHFX=VHFX+DOTP(NGW,C0(1,IA),C2_ace(1,IA))
  !      ENDDO
!---------------------------------------------------------------
!   sagar hack TODO
!    IF (parai%cp_nogrp > 1) THEN
       !CALL tiset(procedureN//'_b',isub6)
       !CALL mp_sum(ehfx,parai%cp_inter_grp)
       !IF (redist_c2)
!        CALL cp_grp_redist(C2_ace,ncpw%ngw,nstate)
       !CALL tihalt(procedureN//'_b',isub6)
!    ENDIF
!---------------------------------------------------------------
        if(parai%cp_inter_me.eq.0) & !hack fixed for CP_GRP sagar
        CALL add_wfn(jgw,nstate,zone,C2_ace,ncpw%ngw,c2,ncpw%ngw)

!  C     ENERGY
  !      EHFX=EHFX*OMEGA
  !      CALL MY_SUM_D(EHFX,1,ALLGRP)
        EHFX  = ehfx_ace

!  C     POTENTIAL
        DO I=1,NSTATE
           !VHFX=VHFX+DOTP(NGW,C0(1,IA),C2_hfx(1,IA))
           VHFX=VHFX+DOTP(ncpw%ngw,C0(:,I),C2(:,I))
        ENDDO
        !CALL MY_SUM_D(VHFX,1,ALLGRP)
        CALL mp_sum(vhfx,parai%allgrp)
        !
        !IF(me.eq.0)write(6,*)"EHFX=",EHFX, "VHFX=",VHFX
!========================================================================
! ACE stuff ends here
!________________________________________________________________________
    DEALLOCATE(c2_ace,cmexx, &
         mexx,STAT=ierr)
    IF (ierr /= 0) CALL stopgm(procedureN,'Deallocation problem',&
         __LINE__,__FILE__)

! ==================================================================

    CALL tihalt(procedureN,isub)
! ==--------------------------------------------------------------==
  END SUBROUTINE hfx_ace
! ==================================================================

  SUBROUTINE hfx_old(c0,c2,f,psia,nstate,ehfx,vhfx,deeq_fnl_hfx,fion,tfor)
    ! ==================================================================
    ! == HARTREE-FOCK EXCHANGE                                        ==
    ! ==--------------------------------------------------------------==

    COMPLEX(real_8),INTENT(INOUT) __CONTIGUOUS :: c2(:,:),c0(:,:)
    REAL(real_8),INTENT(IN) __CONTIGUOUS       :: f(:)
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS    :: fion(:,:,:)
    REAL(real_8),INTENT(inOUT) __CONTIGUOUS      :: deeq_fnl_hfx(:,:,:)
    COMPLEX(real_8),INTENT(INOUT) __CONTIGUOUS :: psia(:)
    INTEGER,INTENT(IN)                         :: nstate
    REAL(real_8),INTENT(OUT)                   :: ehfx, vhfx
    LOGICAL,INTENT(IN)                         :: tfor

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfx_old'
    COMPLEX(real_8), PARAMETER :: zone = (1.0_real_8,0.0_real_8), &
      zzero = (0.0_real_8,0.0_real_8)

    COMPLEX(real_8), ALLOCATABLE             :: psic(:), vpotg(:)
    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C2_hfx, psib
    INTEGER :: geo_iter, i, ia, ib1, ib2, ibin, ibuf, id_states(2,2), ierr, &
      ig, int_i, int_ij, int_j, ispin, isub, isub2, isub3, isub4, isub5, &
      isub6, iwf, j, jb1, jb2, my_pcol, my_prow, n_int, nbr_states(2), &
      npcols, nprows, nspins, nst, nstates(2), st_offst
    INTEGER(int_8) :: nbr_int_skiped_1, nbr_int_skiped_2, nbr_integrals, &
      nbr_rwfn_precomputed, nbr_rwfn_recomputed
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bin_vals, new_vals, &
                                                psi_in_core
    INTEGER, SAVE                            :: int_count = 0, &
                                                prev_geo_id = -HUGE(0)
    LOGICAL                                  :: init_ints, no_more_states
    REAL(real_8) :: bin_max, bin_range, dab, dt, dt_tot, EHFX_loc, &
      EHFX_loc_1, EHFX_loc_2, GHFX_loc, GHFX_loc_1, GHFX_loc_2, max_dab, &
      old_DWFC, pfl, pfx, pfx1, pfx2, t1, t1_tot, t2, t2_tot
    REAL(real_8), ALLOCATABLE                :: vpotr(:)
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: acc_vals, max_vals
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :), SAVE                  :: int_vals

    ehfx=0._real_8
    vhfx=0._real_8
    IF (func1%mhfx.EQ.0) RETURN
    ! ==--------------------------------------------------------------==
    CALL tiset(procedureN,isub)

    ! decide if reinit intergrals
    ! need to use the *right* geo iter counter (in the cpmd spirit)!
    init_ints=.FALSE.
    geo_iter=0
    IF(.NOT.cntl%wfopt) geo_iter=infi
    IF (prev_geo_id/=geo_iter) THEN
       prev_geo_id=geo_iter
       init_ints=.TRUE.
    ENDIF


    IF (paral%io_parent.AND.hfxc3%twscr.AND..FALSE.) THEN
       WRITE(6,*) '>>>>>>>>>>     WFOPT=',cntl%wfopt
       WRITE(6,*) '>>>>>>>>>>       NFI=',iteropt%nfi
       WRITE(6,*) '>>>>>>>>>>      INFI=',infi!<
       WRITE(6,*) '>>>>>>>>>>     IINFI=',iteropt%iinfi
       WRITE(6,*) '>>>>>>>>>>      INFW=',infw
       WRITE(6,*) '>>>>>>>>>> init_ints=',init_ints
    ENDIF

    IF (paral%io_parent.AND.init_ints.AND.hfxc3%twscr.AND..FALSE.) THEN
       WRITE(6,*) 'TWSCR=',hfxc3%twscr
       WRITE(6,*) 'TWFT=',hfxc3%twft,&
            'DWF_INTEGRAL_THRESH=',hfxc4%dwf_integral_thresh
       WRITE(6,*) 'TWFC=',hfxc3%twfc,' DWFC=',hfxc4%dwfc
       WRITE(6,*) 'TDIAW=',hfxc3%tdiaw
    ENDIF

    t1_tot=m_cputime()

    IF (cntl%cdft.AND.parai%cp_nogrp>1) CALL stopgm(procedureN,&
         'CDFT AND CP_NOGRP>1 NOT YET IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (tkpts%tkpnt) CALL stopgm(procedureN,'K-POINTS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (lspin2%tlse) CALL stopgm(procedureN,'NO LSE POSSIBLE',& 
         __LINE__,__FILE__)
    IF (group%nogrp.GT.1) CALL stopgm(procedureN,&
         'TASK GROUPS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    ! 
    CALL setfftn(ipoolhfx)
    ! 
    ! Prepare the indeces for parallelization over the groups
    ! 
    CALL proc_to_grid2d(parai%cp_nogrp,parai%cp_inter_me,nprows,npcols,&
         my_prow,my_pcol)

    ! 
    ! Need to allocate a buffer to accumulate the x part of C2
    ! 
    ALLOCATE(C2_hfx(ncpw%ngw,nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
         __LINE__,__FILE__)

    CALL zero_wfn(jgw,nstate,C2_hfx,ncpw%ngw)
    ! 
    ! For the moment we synchronize the C0 and C2 if the cp groups are used
    ! That should be done else where
    IF (parai%cp_nogrp.GT.1) THEN
       CALL tiset(procedureN//'_a',isub5)
       CALL mp_bcast(c0,ncpw%ngw*nstate,0,parai%cp_inter_grp)
       CALL mp_bcast(c2,ncpw%ngw*nstate,0,parai%cp_inter_grp)
       CALL tihalt(procedureN//'_a',isub5)
    ENDIF

    ! 
    ! randomize the states to improve load balance
    ! 
    ! TODO

    ! 
    t1=m_cputime()

    ! CALL check_occupation(F,NSTATE,cntl%tlsd)

    ! 
    ! set the spins stuff
    ! 
    nstates = 0
    IF (cntl%tlsd) THEN
       nspins = 2
       nstates(1) = spin_mod%nsup
       nstates(2) = spin_mod%nsdown
       pfl = 0.5_real_8
    ELSE
       nspins = 1
       nstates(1) = nstate
       pfl = 0.25_real_8
    ENDIF
    IF (cntl%thybrid) pfl=pfl*func3%phfx
    ALLOCATE(psic(maxfftn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpotg(2*jhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpotr(2*llr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi_in_core(0:nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
         __LINE__,__FILE__)
    psi_in_core(:) = 0

    IF (hfxc3%twscr) THEN
       IF (.NOT.ALLOCATED(int_vals)) THEN
          ALLOCATE(int_vals( nstate*(nstate+1)/2, 2 ),stat=ierr)
          IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
               __LINE__,__FILE__)
          CALL zeroing(int_vals)!,SIZE(int_vals))
       ENDIF
       ALLOCATE(acc_vals( nstate*(nstate+1)/2, 2 ),&
            new_vals( nstate*(nstate+1)/2 ), stat=ierr )
       IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
            __LINE__,__FILE__)
       CALL zeroing(acc_vals)!,SIZE(acc_vals))
       CALL zeroing(new_vals)!,SIZE(new_vals))

       int_count=int_count+1
    ENDIF
    IF (func1%mhfx.EQ.1) THEN
       ! HARTREE-FOCK
       ALLOCATE(psib(maxfftn,2),stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
            __LINE__,__FILE__)
       CALL zeroing(psib)!,2*maxfftn)

       ! ==--------------------------------------------------------------==
       CALL tiset(procedureN//'_outer',isub2)
       ! 
       ! loop over the spin
       ! 
       id_states(:,:) = 0
       nbr_states(:) = 0
       nbr_int_skiped_1 = 0; nbr_int_skiped_2 = 0;
       nbr_integrals = 0
       nbr_rwfn_precomputed = 0
       nbr_rwfn_recomputed = 0
       DO ispin = 1,nspins

          nst = nstates(ispin)
          st_offst = 0
          IF (ispin.EQ.2) st_offst = nstates(1)



          ! ==--------------------------------------------------------------==
          ! ==--------------------------------------------------------------==
          ! 
          ! if possible precompute some real space wavefunctions
          ! 
          IF (cntl%krwfn) THEN
             CALL tiset(procedureN//'_precomp',isub4)

             IF (.FALSE.) THEN

                ! >>>>>>>>>>>>>> OLD
                iwf = 1
                psi_in_core(:) = 0
                DO ia=st_offst+1,st_offst+nst,2
                   IF (iwf.GT.maxstatesx/2) EXIT
                   CALL zeroing(psia)!,maxfftn)
                   IF (ia.EQ.nstate) THEN
                      CALL set_psi_1_state_g(zone,c0(:,ia),psia)
                      psi_in_core(ia) = iwf
                   ELSE
                      CALL set_psi_2_states_g(c0(:,ia),c0(:,ia+1),psia)
                      psi_in_core(ia)   = iwf
                      psi_in_core(ia+1) = -iwf
                   ENDIF
                   nbr_rwfn_precomputed = nbr_rwfn_precomputed + 1
                   CALL invfftn(psia,.TRUE.,parai%allgrp)
                   CALL dcopy(2*llr1,psia(1),1,rswfx(1,iwf),1)
                   iwf = iwf + 1
                ENDDO         ! IA

             ELSE

                ! >>>>>>>>>>>>>> NEW
                iwf = 1
                psi_in_core(:) = 0
                i_loop: DO i = 1,part_1d_nbr_elems(nst,my_prow,nprows)
                   ia = part_1d_get_elem(i,my_prow,nprows)+st_offst
                   IF (f(ia).LT.1.e-6_real_8) CYCLE
                   IF (psi_in_core(ia).EQ.0) THEN
                      ! here compute ia i f needed
                      ! 
                      CALL zeroing(psia)!,maxfftn)
                      ! psi_in_core(ia) ->  iwf ie psi(ia) stored in real at position iwf
                      ! psi_in_core(ia) -> -iwf ie psi(ia) stored in imag at position iwf
                      ! psi_in_core(ia) ->    0 ie psi(ia) not stored
                      CALL set_psi_1_state_g(zone,c0(:,ia),psia)
                      nbr_rwfn_precomputed = nbr_rwfn_precomputed + 1
                      CALL invfftn(psia,.TRUE.,parai%allgrp)
                      IF (ABS(iwf).LE.maxstatesx/2) THEN
                         psi_in_core(ia) = iwf
                         IF (iwf.LT.0) THEN
                            CALL copy_re_to_im(llr1,psia,rswfx(:,-iwf))
                            iwf = -iwf + 1
                         ELSE
                            CALL copy_re_to_re(llr1,psia,rswfx(:, iwf))
                            iwf = -iwf
                         ENDIF
                      ELSE
                         ! We can exit the loops
                         EXIT i_loop
                      ENDIF
                   ENDIF

                   j = 1
                   DO
                      no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol,npcols)
                      IF (no_more_states) EXIT

                      ! find out the correct ib's
                      ib1 = 0; ib2 = 0;
                      DO
                         ! set id1
                         IF (ib1.EQ.0) THEN
                            ib1 = part_1d_get_elem(j,my_pcol,npcols)+st_offst
                            j = j + 1
                         ENDIF
                         IF (ia.EQ.ib1) ib1 = 0
                         IF (.NOT.part_1d_symm_holds_pair(ia,ib1)) ib1 = 0
                         IF (ib1.NE.0) THEN
                            IF (f(ib1).LT.1.e-6_real_8) ib1 = 0
                         ENDIF
                         ! First: check that the separation distance is smaller that the max distance allowed
                         IF (ib1.NE.0.AND.hfxc3%twscr.AND.hfxc3%twfc) THEN
                            CALL get_wannier_separation(ib1,ia,dab)
                            IF (dab.GT.hfxc4%dwfc) THEN
                               int_i=MIN(ia,ib1)
                               int_j=MAX(ia,ib1)
                               int_ij=int_j*(int_j-1)/2+int_i
                               ! Set the int_vals to something that if the 
                               ! pair comes within the radius any time later
                               ! we properly recompute the integral.
                               int_vals(int_ij,:)=1.0_real_8
                               nbr_int_skiped_1 = nbr_int_skiped_1 + 1
                               ib1=0
                            ENDIF
                         ENDIF
                         ! Second: check that the integral is larger than the integral threshold
                         IF (.NOT.init_ints) THEN
                            IF (hfxc3%twscr.AND.ib1.NE.0) THEN
                               int_i=MIN(ia,ib1)
                               int_j=MAX(ia,ib1)
                               int_ij=int_j*(int_j-1)/2+int_i
                               IF (int_vals(int_ij,1)<hfxc4%dwf_integral_thresh) ib1 = 0
                            ENDIF
                         ENDIF
                         IF (psi_in_core(ib1).NE.0) ib1 = 0
                         no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol,&
                              npcols)
                         IF (no_more_states) EXIT
                         ! set ib2
                         IF (ib2.EQ.0) THEN
                            ib2 = part_1d_get_elem(j,my_pcol,npcols)+st_offst
                            j = j + 1
                         ENDIF
                         IF (ia.EQ.ib2) ib2 = 0
                         IF (.NOT.part_1d_symm_holds_pair(ia,ib2)) ib2 = 0
                         IF (ib2.NE.0) THEN
                            IF (f(ib2).LT.1.e-6_real_8) ib2 = 0
                         ENDIF
                         ! First: check that the separation distance is smaller that the max distance allowed
                         IF (ib2.NE.0.AND.hfxc3%twscr.AND.hfxc3%twfc) THEN
                            CALL get_wannier_separation(ib2,ia,dab)
                            IF (dab.GT.hfxc4%dwfc) THEN
                               int_i=MIN(ia,ib2)
                               int_j=MAX(ia,ib2)
                               int_ij=int_j*(int_j-1)/2+int_i
                               ! Set the int_vals to something that if the 
                               ! pair comes within the radius any time later
                               ! we properly recompute the integral.
                               int_vals(int_ij,:)=1.0_real_8
                               nbr_int_skiped_1 = nbr_int_skiped_1 + 1
                               ib2=0
                            ENDIF
                         ENDIF
                         ! Second: check that the integral is larger than the integral threshold
                         IF (.NOT.init_ints) THEN
                            IF (hfxc3%twscr.AND.ib2.NE.0) THEN
                               int_i=MIN(ia,ib2)
                               int_j=MAX(ia,ib2)
                               int_ij=int_j*(int_j-1)/2+int_i
                               IF (int_vals(int_ij,1)<hfxc4%dwf_integral_thresh) ib2 = 0
                            ENDIF
                         ENDIF
                         IF (psi_in_core(ib2).NE.0) ib2 = 0
                         no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol,&
                              npcols)
                         IF (no_more_states) EXIT
                         ! are we set?
                         IF (ib1.NE.0.AND.ib2.NE.0) EXIT
                      ENDDO
                      IF (ib1.EQ.0.AND.ib2.EQ.0) CYCLE
                      IF (ib1.EQ.0) THEN
                         ib1 = ib2; ib2 = 0
                      ENDIF
                      ! compute ib1 and ib2 if needed
                      ! 
                      CALL zeroing(psia)!,maxfftn)
                      ! psi_in_core(ia) ->  iwf ie psi(ia) stored in real at position iwf
                      ! psi_in_core(ia) -> -iwf ie psi(ia) stored in imag at position iwf
                      ! psi_in_core(ia) ->    0 ie psi(ia) not stored
                      IF (ib2.EQ.0) THEN
                         CALL set_psi_1_state_g(zone,c0(:,ib1),psia)
                      ELSE
                         CALL set_psi_2_states_g(c0(:,ib1),c0(:,ib2),psia)
                      ENDIF
                      nbr_rwfn_precomputed = nbr_rwfn_precomputed + 1
                      CALL invfftn(psia,.TRUE.,parai%allgrp)
                      IF (ib2.EQ.0) THEN
                         psi_in_core(ib1) = iwf
                         IF (iwf.LT.0) THEN
                            CALL copy_re_to_im(llr1,psia,rswfx(:,-iwf))
                            iwf = -iwf + 1
                         ELSE
                            CALL copy_re_to_re(llr1,psia,rswfx(:, iwf))
                            iwf = -iwf
                         ENDIF
                      ELSE
                         IF (iwf.LT.0) THEN
                            psi_in_core(ib1) =  iwf
                            CALL copy_re_to_im(llr1,psia,rswfx(:,-iwf))
                            iwf = -iwf + 1
                            IF (ABS(iwf).GT.maxstatesx/2) EXIT i_loop
                            psi_in_core(ib2) =  iwf
                            CALL copy_im_to_re(llr1,psia,rswfx(:, iwf))
                            iwf = -iwf
                         ELSE
                            psi_in_core(ib1) =  iwf
                            psi_in_core(ib2) = -iwf
                            CALL dcopy(2*llr1,psia(1),1,rswfx(1,iwf),1)
                            iwf = iwf + 1
                         ENDIF
                      ENDIF
                      IF (ABS(iwf).GT.maxstatesx/2) EXIT i_loop
                   ENDDO
                ENDDO i_loop ! i
             ENDIF
             CALL tihalt(procedureN//'_precomp',isub4)
          ENDIF
          ! 
          ! ==--------------------------------------------------------------==
          ! ==--------------------------------------------------------------==


          DO i = 1,part_1d_nbr_elems(nst,my_prow,nprows)
             ia = part_1d_get_elem(i,my_prow,nprows)+st_offst
             IF (f(ia).LT.1.e-6_real_8) CYCLE

             IF (psi_in_core(ia).NE.0) THEN
                iwf = psi_in_core(ia)
                IF (iwf.LT.0) THEN
                   iwf=-iwf
                   CALL set_psi_im_state_r(zone,rswfx(:,iwf),&
                        zzero,psia)
                ELSE
                   CALL set_psi_re_state_r(zone,rswfx(:,iwf),&
                        zzero,psia)
                ENDIF
             ELSE
                CALL zeroing(psia)!,maxfftn)
                CALL set_psi_1_state_g(zone,c0(:,ia),psia)
                CALL invfftn(psia,.TRUE.,parai%allgrp)
                nbr_rwfn_recomputed = nbr_rwfn_recomputed + 1
             ENDIF
             ! Handle the diagonal terms
             ! Only when needed....
             IF (part_1d_elem_to_proc(ia-st_offst,nprows).EQ.my_prow.AND.&
                  part_1d_elem_to_proc(ia-st_offst,npcols).EQ.my_pcol) THEN
                IF (parai%me.EQ.0)nbr_integrals=nbr_integrals+1
                pfx=pfl*f(ia)*f(ia)
                CALL hfxaa(EHFX_loc,GHFX_loc,pfx,psia,vpotg,vpotr,psic,&
                     C2_hfx(1,ia),ia,f,deeq_fnl_hfx,fion,tfor)
                ehfx=ehfx+EHFX_loc

                IF (hfxc3%twscr) THEN
                   int_i=MIN(ia,ia)
                   int_j=MAX(ia,ia)
                   int_ij=int_j*(int_j-1)/2+int_i
                   acc_vals( int_ij, 1 ) = EHFX_loc
                   acc_vals( int_ij, 2 ) = GHFX_loc
                   new_vals( int_ij ) = 1
                ENDIF
                IF (hfxc3%twscr.AND.hfxc3%tdiaw) GOTO 101
             ENDIF
             ! Now for the other terms IA-IB
             ! ==--------------------------------------------------------------==
             CALL tiset(procedureN//'_inner',isub3)
             j = 1
             DO
                no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol,&
                     npcols)
                IF (no_more_states) EXIT

                ! 
                ! find out the correct ib's
                ib1 = 0; ib2 = 0;
                DO
                   ! set id1
                   IF (ib1.EQ.0) THEN
                      ib1 = part_1d_get_elem(j,my_pcol,npcols)+st_offst
                      j = j + 1
                   ENDIF
                   IF (ia.EQ.ib1) ib1 = 0
                   IF (.NOT.part_1d_symm_holds_pair(ia,ib1)) ib1 = 0
                   IF (ib1.NE.0) THEN
                      IF (f(ib1).LT.1.e-6_real_8) ib1 = 0
                   ENDIF
                   ! First: check that the separation distance is smaller that the max distance allowed
                   IF (ib1.NE.0.AND.hfxc3%twscr.AND.hfxc3%twfc) THEN
                      CALL get_wannier_separation(ib1,ia,dab)
                      IF (dab.GT.hfxc4%dwfc) THEN
                         int_i=MIN(ia,ib1)
                         int_j=MAX(ia,ib1)
                         int_ij=int_j*(int_j-1)/2+int_i
                         ! Set the int_vals to something that if the 
                         ! pair comes within the radius any time later
                         ! we properly recompute the integral.
                         int_vals(int_ij,:)=1.0_real_8
                         nbr_int_skiped_1 = nbr_int_skiped_1 + 1
                         ib1=0
                      ENDIF
                   ENDIF
                   ! Second: check that the integral is larger than the integral threshold
                   IF (.NOT.init_ints) THEN
                      IF (hfxc3%twscr.AND.ib1.NE.0) THEN
                         int_i=MIN(ia,ib1)
                         int_j=MAX(ia,ib1)
                         int_ij=int_j*(int_j-1)/2+int_i
                         IF (int_vals(int_ij,1)<hfxc4%dwf_integral_thresh) THEN
                            nbr_int_skiped_1 = nbr_int_skiped_1 + 1
                            ib1 = 0
                         ENDIF
                      ENDIF
                   ENDIF
                   no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol,&
                        npcols)
                   IF (no_more_states) EXIT
                   ! set ib2
                   IF (ib2.EQ.0) THEN
                      ib2 = part_1d_get_elem(j,my_pcol,npcols)+st_offst
                      j = j + 1
                   ENDIF
                   IF (ia.EQ.ib2) ib2 = 0
                   IF (.NOT.part_1d_symm_holds_pair(ia,ib2)) ib2 = 0
                   IF (ib2.NE.0) THEN
                      IF (f(ib2).LT.1.e-6_real_8) ib2 = 0
                   ENDIF
                   ! First: check that the separation distance is smaller that the max distance allowed
                   IF (ib2.NE.0.AND.hfxc3%twscr.AND.hfxc3%twfc) THEN
                      CALL get_wannier_separation(ib2,ia,dab)
                      IF (dab.GT.hfxc4%dwfc) THEN
                         int_i=MIN(ia,ib2)
                         int_j=MAX(ia,ib2)
                         int_ij=int_j*(int_j-1)/2+int_i
                         ! Set the int_vals to something that if the 
                         ! pair comes within the radius any time later
                         ! we properly recompute the integral.
                         int_vals(int_ij,:)=1.0_real_8
                         nbr_int_skiped_1 = nbr_int_skiped_1 + 1
                         ib2=0
                      ENDIF
                   ENDIF
                   ! Second: check that the integral is larger than the integral threshold
                   IF (.NOT.init_ints) THEN
                      IF (hfxc3%twscr.AND.ib2.NE.0) THEN
                         int_i=MIN(ia,ib2)
                         int_j=MAX(ia,ib2)
                         int_ij=int_j*(int_j-1)/2+int_i
                         IF (int_vals(int_ij,1)<hfxc4%dwf_integral_thresh) THEN
                            nbr_int_skiped_1 = nbr_int_skiped_1 + 1
                            ib2 = 0
                         ENDIF
                      ENDIF
                   ENDIF
                   no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol,&
                        npcols)
                   IF (no_more_states) EXIT
                   ! are we set?
                   IF (ib1.NE.0.AND.ib2.NE.0) EXIT
                ENDDO

                ! if(me.eq.0) write(6,*)cp_me,'ia,ib1,ib2',ia,ib1,ib2,no_more_states

                IF (.NOT.(ib1.EQ.0.AND.ib2.EQ.0)) THEN
                   IF ((psi_in_core(ib1).NE.0.AND.psi_in_core(ib2).NE.0)&
                        .AND.(ib1.NE.0.OR.ib2.NE.0)) THEN
                      CALL zeroing(psib(:,1))!,maxfftn)
                      IF (ib1.NE.0) THEN
                         iwf = psi_in_core(ib1)
                         IF (iwf.LT.0) THEN
                            iwf=-iwf
                            CALL set_psi_im_state_r(zone,rswfx(:,iwf),&
                                 zone,psib(:,1))
                         ELSE
                            CALL set_psi_re_state_r(zone,rswfx(:,iwf),&
                                 zone,psib(:,1))
                         ENDIF
                      ENDIF
                      IF (ib2.NE.0) THEN
                         iwf = psi_in_core(ib2)
                         IF (iwf.LT.0) THEN
                            iwf=-iwf
                            CALL set_psi_im_state_r(uimag,&
                                 rswfx(:,iwf),zone,psib(:,1))
                         ELSE
                            CALL set_psi_re_state_r(uimag,&
                                 rswfx(:,iwf),zone,psib(:,1))
                         ENDIF
                      ENDIF
                   ELSE
                      IF (ib1.EQ.0) THEN
                         CALL zeroing(psib(:,1))!,maxfftn)
                         CALL set_psi_1_state_g(uimag,c0(:,ib2),psib(:,1))
                         nbr_rwfn_recomputed = nbr_rwfn_recomputed + 1
                      ELSEIF (ib2.EQ.0) THEN
                         CALL zeroing(psib(:,1))!,maxfftn)
                         CALL set_psi_1_state_g(zone,c0(:,ib1),psib(:,1))
                         nbr_rwfn_recomputed = nbr_rwfn_recomputed + 1
                      ELSE
                         CALL zeroing(psib(:,1))!,maxfftn)
                         CALL set_psi_2_states_g(c0(:,ib1),c0(:,ib2),&
                              psib(:,1))
                         nbr_rwfn_recomputed = nbr_rwfn_recomputed + 2
                      ENDIF
                      ! Transform the wavefunction to real space
                      CALL invfftn(psib(:,1),.TRUE.,parai%allgrp)
                   ENDIF
                ENDIF
                ! 
                ! ==--------------------------------------------------------------==
                id_states(1,1) = ib1
                id_states(2,1) = ib2
                nbr_states(1) = 0
                IF (ib1.NE.0) nbr_states(1) = nbr_states(1) + 1
                IF (ib2.NE.0) nbr_states(1) = nbr_states(1) + 1

                IF (nbr_states(1).EQ.1) THEN
                   IF (id_states(1,1).NE.0) THEN
                      IF (id_states(1,2).EQ.0) THEN
                         id_states(1,2) = id_states(1,1)
                         CALL copy_re_to_re(llr1,psib(:,1),psib(:,2))
                      ELSE
                         id_states(2,2) = id_states(1,1)
                         CALL copy_re_to_im(llr1,psib(:,1),psib(:,2))
                      ENDIF
                      id_states(1,1) = 0
                   ELSE
                      IF (id_states(1,2).EQ.0) THEN
                         id_states(1,2) = id_states(2,1)
                         CALL copy_im_to_re(llr1,psib(:,1),psib(:,2))
                      ELSE
                         id_states(2,2) = id_states(2,1)
                         CALL copy_im_to_im(llr1,psib(:,1),psib(:,2))
                      ENDIF
                      id_states(2,1) = 0
                   ENDIF
                   nbr_states(1) = nbr_states(1) - 1
                   nbr_states(2) = nbr_states(2) + 1
                ENDIF

                ! ==--------------------------------------------------------------==

                DO ibuf = 1,2
                   IF (.NOT.no_more_states) THEN
                      IF (nbr_states(ibuf).NE.2) CYCLE
                   ENDIF
                   jb1 = id_states(1,ibuf)
                   jb2 = id_states(2,ibuf)
                   ! 
                   EHFX_loc = 0.0_real_8
                   IF (jb1.NE.0 .AND. jb2.NE.0) THEN
                      nbr_integrals=nbr_integrals+2
                      pfx1=pfl*f(ia)*f(jb1)
                      pfx2=pfl*f(ia)*f(jb2)
                      CALL hfxab2(EHFX_loc_1,EHFX_loc_2,&
                           GHFX_loc_1,GHFX_loc_2,&
                           pfx1,pfx2,psia,&
                           psib(1,ibuf),&
                           vpotg,vpotr,psic,C2_hfx(1,ia),&
                           C2_hfx(1,jb1),c2_hfx(1,jb2),ia,&
                           jb1,jb2,f,deeq_fnl_hfx,fion,tfor)
                      EHFX_loc = EHFX_loc_1 + EHFX_loc_2
                      IF (hfxc3%twscr) THEN
                         int_i=MIN(ia,jb1)
                         int_j=MAX(ia,jb1)
                         int_ij=int_j*(int_j-1)/2+int_i
                         acc_vals( int_ij, 1 ) = EHFX_loc_1
                         acc_vals( int_ij, 2 ) = GHFX_loc_1
                         new_vals( int_ij ) = 1
                         int_i=MIN(ia,jb2)
                         int_j=MAX(ia,jb2)
                         int_ij=int_j*(int_j-1)/2+int_i
                         acc_vals( int_ij, 1 ) = EHFX_loc_2
                         acc_vals( int_ij, 2 ) = GHFX_loc_2
                         new_vals( int_ij ) = 1
                      ENDIF
                   ELSEIF (jb1.NE.0 .AND. jb2.EQ.0) THEN
                      ! Terms IA*IB1
                      nbr_integrals=nbr_integrals+1
                      pfx=pfl*f(ia)*f(jb1)
                      CALL hfxab(EHFX_loc,GHFX_loc,&
                           pfx,psia,psib(1,ibuf),1,&
                           vpotg,vpotr,psic,C2_hfx(1,ia),&
                           C2_hfx(1,jb1),ia,jb1,f,deeq_fnl_hfx,&
                           fion,tfor)
                      IF (hfxc3%twscr) THEN
                         int_i=MIN(ia,jb1)
                         int_j=MAX(ia,jb1)
                         int_ij=int_j*(int_j-1)/2+int_i
                         acc_vals( int_ij, 1 ) = EHFX_loc
                         acc_vals( int_ij, 2 ) = GHFX_loc
                         new_vals( int_ij ) = 1
                      ENDIF
                   ELSEIF (jb1.EQ.0 .AND. jb2.NE.0) THEN
                      ! Terms IA*IB2
                      nbr_integrals=nbr_integrals+1
                      pfx=pfl*f(ia)*f(jb2)
                      CALL hfxab(EHFX_loc,GHFX_loc,pfx,&
                           psia,psib(1,ibuf),2,&
                           vpotg,vpotr,psic,C2_hfx(1,ia),&
                           C2_hfx(1,jb2),ia,jb2,f,deeq_fnl_hfx,&
                           fion,tfor)
                      IF (hfxc3%twscr) THEN
                         int_i=MIN(ia,ib2)
                         int_j=MAX(ia,jb2)
                         int_ij=int_j*(int_j-1)/2+int_i
                         acc_vals( int_ij, 1 ) = EHFX_loc
                         acc_vals( int_ij, 2 ) = GHFX_loc
                         new_vals( int_ij ) = 1
                      ENDIF
                   ENDIF
                   ehfx=ehfx+EHFX_loc

                   ! reset the buffer
                   id_states(:,ibuf) = 0
                   nbr_states(ibuf) = 0
                ENDDO     ! ibuf
             ENDDO           ! IB
             CALL tihalt(procedureN//'_inner',isub3)
             ! ==--------------------------------------------------------------==
101          CONTINUE
          ENDDO                 ! IA
       ENDDO                   ! ispin
       CALL tihalt(procedureN//'_outer',isub2)
       ! ==--------------------------------------------------------------==
    ELSE
       ! HARTREE
       DO ia=1,nstate
          IF (f(ia).GT.1.e-6_real_8) THEN
             CALL zeroing(psia)!,maxfftn)
             !ocl novrec
             DO ig=1,jgw
                psia(nzfs(ig))=c0(ig,ia)
                psia(inzs(ig))=CONJG(c0(ig,ia))
             ENDDO
             IF (geq0) psia(nzfs(1))=c0(1,ia)
             ! Transform the wavefunction to real space
             CALL invfftn(psia,.TRUE.,parai%allgrp)
             pfx=pfl*f(ia)*f(ia)
             CALL hfxaa(EHFX_loc,GHFX_loc,pfx,psia,vpotg,vpotr,&
                  psic,C2_hfx(1,ia),ia,f,deeq_fnl_hfx,fion,tfor)
             ehfx=ehfx+EHFX_loc
          ENDIF
       ENDDO
    ENDIF
    ! free some memory
    ! 
    DEALLOCATE(psic,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpotg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpotr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF (func1%mhfx.EQ.1) THEN
       DEALLOCATE(psib,stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
            __LINE__,__FILE__)
    ENDIF

    t2=m_cputime()
    dt = t2 - t1

    ! 
    ! redistribute EHFX and C2_hfx over the groups if needed
    ! 
    IF (parai%cp_nogrp.GT.1) THEN
       ! call mp_sync(CP_GRP)
       CALL tiset(procedureN//'_b',isub6)
       CALL mp_sum(ehfx,parai%cp_inter_grp)
       ! call cp_grp_redist_z(C2_hfx,NGW,NSTATE)
       CALL mp_sum(C2_hfx,ncpw%ngw*nstate,parai%cp_inter_grp)
       CALL tihalt(procedureN//'_b',isub6)
    ENDIF
    ! 
    ! add up the hfx contribution to C2
    ! 
    CALL add_wfn(jgw,nstate,zone,C2_hfx,ncpw%ngw,c2,ncpw%ngw)
    IF(pslo_com%tivan)CALL fnonloc_hfx(c2,nstate)
    
    IF (hfxc3%twscr) THEN
       CALL mp_sum(acc_vals,SIZE(acc_vals),parai%cp_grp)
       CALL mp_sum(new_vals,SIZE(new_vals),parai%cp_inter_grp)

       bin_range=0.5_real_8
       bin_max=200.0_real_8
       ALLOCATE( max_vals( INT( bin_max / bin_range ) + 1, 2 ),&
            bin_vals( INT( bin_max / bin_range ) + 1 ), stat=ierr )
       IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',& 
            __LINE__,__FILE__)

       max_vals(:,:) = 0.0_real_8
       bin_vals(:) = 0
       n_int=0
       max_dab=0.0_real_8

       ! The following compiler has a problem with the reduction in the omp section
       ! IBM XL Fortran Advanced Edition for Blue Gene/P, V11.1
       ! Version: 11.01.0000.0011
       !$omp parallel do &
       !$omp             default(none) &
       !$omp             private(i,ibin,int_i,int_j,dab) &
       !$omp             shared(new_vals,acc_vals,bin_range) &
       !$omp             shared(int_vals) &
       !$omp             reduction(+:n_int,bin_vals) &
       !$omp             reduction(max:max_vals,max_dab)
       DO i=1,SIZE(int_vals,1)
          IF (new_vals(i)>0) THEN
             int_vals(i,1) = ABS(acc_vals(i,1))
             int_vals(i,2) = acc_vals(i,2)
             n_int=n_int+1
          ENDIF
          int_j=CEILING((-1.0_real_8+SQRT(1.0_real_8+8.0_real_8*REAL(i,kind=real_8)))/2.0_real_8)
          int_i=i-int_j*(int_j-1)/2
          CALL get_wannier_separation(int_i,int_j,dab)
          max_dab=MAX(max_dab,dab)
          ibin=INT(dab/bin_range)+1
          max_vals(ibin,1)=MAX(max_vals(ibin,1),int_vals(i,1))
          max_vals(ibin,2)=MAX(max_vals(ibin,2),int_vals(i,2))
          bin_vals(ibin)=bin_vals(ibin)+1
       ENDDO
       !$omp end parallel do

       ! find new radius      
       old_DWFC=hfxc4%dwfc
       IF (hfxc3%twfc.AND.hfxc3%twft) THEN
          ibin=INT(hfxc4%dwfc/bin_range)+1
          IF (max_vals(ibin  ,1)>hfxc4%dwf_integral_thresh)&
               hfxc4%dwfc=hfxc4%dwfc+bin_range
          IF (max_vals(ibin-1,1)<hfxc4%dwf_integral_thresh.AND.&
               max_vals(ibin-2,1)<hfxc4%dwf_integral_thresh)&
               hfxc4%dwfc=hfxc4%dwfc-bin_range
       ENDIF
       IF (paral%io_parent.AND..FALSE.) THEN
          WRITE(6,*) 'we have ',n_int,' intergrals (re)set'
          WRITE(6,*) 'new DWFC=',hfxc4%dwfc,' old DWFC=',old_DWFC
          DO i=1,INT(max_dab/bin_range)+1
             WRITE(6,'(2F6.1,2E8.2,I8)') (i-1)*bin_range,&
                  i*bin_range,&
                  max_vals(i,:),bin_vals(i)
          ENDDO
       ENDIF

       IF (paral%io_parent.AND..FALSE.) THEN
          DO i=1,SIZE(int_vals,1)
             int_j=CEILING((-1.0_real_8+SQRT(1.0_real_8+8.0_real_8*REAL(i,kind=real_8)))/2.0_real_8)
             int_i=i-int_j*(int_j-1)/2
             CALL get_wannier_separation(int_i,int_j,dab)
             WRITE(6,'(A,I0,5E12.4)') 'INT',int_count,dab,&
                  wcentx(4,int_i),wcentx(4,int_j),&
                  int_vals(i,:)
          ENDDO
       ENDIF

       DEALLOCATE( max_vals,bin_vals,new_vals,acc_vals,stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
            __LINE__,__FILE__)

    ENDIF

    ! ENERGY
    ehfx=ehfx*parm%omega
    CALL mp_sum(ehfx,parai%allgrp)
    ! POTENTIAL
    DO ia=1,nstate
       vhfx=vhfx+dotp(ncpw%ngw,c0(:,ia),c2(:,ia))
    ENDDO
    CALL mp_sum(vhfx,parai%allgrp)
    DEALLOCATE(psi_in_core,C2_hfx,stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',& 
         __LINE__,__FILE__)

    ! 
    ! some infos
    ! 
    t2_tot=m_cputime()
    dt_tot = t2_tot - t1_tot

    IF (PRINT_GROUP_INFOS) THEN
       IF (parai%me.EQ.0) THEN
          WRITE(6,'(1X,6(A,I0),A,F0.2,A,F0.2)')&
               procedureN//'| group ',parai%cp_inter_me,&
               ' computed ',nbr_integrals,&
               ' integrals, precomputed ',nbr_rwfn_precomputed,&
               ', recomputed ',nbr_rwfn_recomputed,&
               ' and skipped ',&
               nbr_int_skiped_1,' + ',nbr_int_skiped_2,&
               ', t_loc ',dt,&
               ', t_per_1k_ints ',1.0e3_real_8*dt/(REAL(nbr_integrals,kind=real_8)+1.0e-6_real_8) ! to avoid NANs
          CALL m_flush(6)
       ENDIF

       IF (parai%cp_nogrp.GT.1) THEN
          CALL mp_sum(nbr_integrals,parai%cp_inter_grp)
          CALL mp_sum(nbr_int_skiped_1,parai%cp_inter_grp)
          CALL mp_sum(nbr_int_skiped_2,parai%cp_inter_grp)
          CALL mp_sum(nbr_rwfn_precomputed,parai%cp_inter_grp)
          CALL mp_sum(nbr_rwfn_recomputed,parai%cp_inter_grp)
          CALL mp_max(dt_tot,parai%cp_inter_grp)
          IF (parai%cp_me.EQ.0) THEN
             WRITE(6,'(1X,5(A,I0),A,F0.2)')&
                  procedureN//'| all the groups computed ',&
                  nbr_integrals,' integrals, precomputed ',&
                  nbr_rwfn_precomputed,&
                  ', recomputed ',nbr_rwfn_recomputed,&
                  ' and skipped ',&
                  nbr_int_skiped_1,' + ',nbr_int_skiped_2,&
                  ', t_tot ',dt_tot
             CALL m_flush(6)
          ENDIF
       ENDIF
    ENDIF

    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfx_old
  ! ==================================================================
  SUBROUTINE check_occupation(f,nstate,tlsd)
    ! ==================================================================
    ! == Check the occupation number                                  ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: nstate
    REAL(real_8)                             :: f(nstate)
    LOGICAL                                  :: tlsd

    REAL(real_8), PARAMETER                  :: one = 1.0_real_8, &
                                                two = 2.0_real_8 , &
                                                zero = 0.0_real_8

    INTEGER                                  :: i
    REAL(real_8)                             :: ff

! ==--------------------------------------------------------------==

    IF (tlsd) THEN
       !$omp parallel do default(none) private(I,FF) shared(F,NSTATE)
       DO i=1,nstate
          ff=ABS(f(i)-one)
          IF (ff.GT.zero.AND.ABS(f(i)).GT.zero) THEN
             CALL stopgm('HFX','ONLY OCCUPATIONS (0,1) ALLOWED',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
       !$omp end parallel do
    ELSE
       !$omp parallel do default(none) private(I,FF) shared(F,NSTATE)
       DO i=1,nstate
          ff=ABS(f(i)-two)
          IF (ff.GT.zero .AND. ABS(f(i)).GT.zero) THEN
             CALL stopgm('HFX','ONLY OCCUPATIONS (0,2) ALLOWED',& 
                  __LINE__,__FILE__)
          ENDIF
       ENDDO
       !$omp end parallel do
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE check_occupation

  ! ==================================================================
  SUBROUTINE get_wannier_separation(icenter,jcenter,dist)
    ! ==================================================================
    ! == Get the distance between two Wannier centers                 ==
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: icenter, jcenter
    REAL(real_8)                             :: dist

    REAL(real_8)                             :: coordab(3),coordab_(3)

! we should protect here
! IF(ICENTER.GT.NSTATE.OR.ICENTER.LT.1.OR.
! &   JCENTER.GT.NSTATE.OR.JCENTER.LT.1)
! &   CALL STOPGM('GET_WANNIER_SEPARATION',' something wrong here...')

    coordab_(1)=wcentx(1,icenter)-wcentx(1,jcenter)
    coordab_(2)=wcentx(2,icenter)-wcentx(2,jcenter)
    coordab_(3)=wcentx(3,icenter)-wcentx(3,jcenter)
    CALL pbc(coordab_(1),coordab_(2),coordab_(3),coordab(1),coordab(2),coordab(3),1,parm%apbc,parm%ibrav)
    dist=SQRT(coordab(1)**2+coordab(2)**2+coordab(3)**2)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE get_wannier_separation

  ! ==================================================================
  SUBROUTINE hfxab(ehfx,ghfx,pf,psia,psib,iran,vpotg,vpotr,psic,&
       c2a,c2b,ia,ib,f,deeq_fnl_hfx,fion,tfor)
    ! ==================================================================
    REAL(real_8),INTENT(OUT)                 :: ehfx, ghfx
    REAL(real_8),INTENT(IN)                  :: pf
    COMPLEX(real_8),INTENT(IN)               :: psia(llr1), psib(llr1)
    INTEGER,INTENT(IN)                       :: iran, ia, ib
    COMPLEX(real_8),INTENT(OUT)              :: vpotg(jhg)
    REAL(real_8),INTENT(OUT)                 :: vpotr(llr1)
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: psic(:)
    COMPLEX(real_8),INTENT(INOUT)            :: c2a(jgw), c2b(jgw)
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: f(:)
    REAL(real_8),INTENT(inOUT) __CONTIGUOUS    :: deeq_fnl_hfx(:,:,:)
    LOGICAL,INTENT(IN)                       :: tfor
    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxab'

    COMPLEX(real_8)                          :: fm, fp, geq0_vpotg
    INTEGER                                  :: ig, ir, isub, nstates(2,1)

    CALL tiset(procedureN,isub)
    ehfx=0.0_real_8
    ghfx=0.0_real_8
    IF (iran.EQ.1) THEN
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=REAL(psia(ir))*REAL(psib(ir))
       ENDDO
    ELSE
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=REAL(psia(ir))*AIMAG(psib(ir))
       ENDDO
    ENDIF
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)

    IF(pslo_com%tivan)THEN
       nstates(1,1)=ia
       nstates(2,1)=ib
       CALL rhov(nstates,psi=vpotg,hfx=.TRUE.)
    ELSE
       !$omp workshare
       vpotg=CMPLX(0.0_real_8,0.0_real_8)
       !$omp end workshare
    END IF
    IF(geq0)geq0_vpotg=vpotg(1)

    !$omp parallel do private(IG,FP) &
    !$omp  reduction(+:EHFX)
    DO ig=1,jhg
       fp=psic(nzff(ig))+vpotg(ig)
       vpotg(ig)=-pf*scgx(ig)*fp
       ehfx=ehfx+REAL(4._real_8*vpotg(ig)*CONJG(fp))
    ENDDO
    IF(geq0) ehfx=ehfx-REAL(2._real_8*vpotg(1)*CONJG(psic(nzff(1))+geq0_vpotg))
    IF(pslo_com%tivan)CALL newd(deeq_fnl_hfx,f,vpotg,fion,tfor,nstates,.TRUE.)

    CALL zeroing(psic)!,maxfftn)
    !$omp parallel do private(IG)
    DO ig=1,jhg
       psic(nzff(ig))=vpotg(ig)
       psic(inzf(ig))=CONJG(vpotg(ig))
    ENDDO
    IF (geq0) psic(nzff(1))=vpotg(1)
    CALL invfftn(psic,.FALSE.,parai%allgrp)
    !$omp parallel do private(IR)
    DO ir=1,llr1
       vpotr(ir)=REAL(psic(ir))
    ENDDO
    IF (iran.EQ.1) THEN
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=vpotr(ir)*(REAL(psia(ir))+uimag*REAL(psib(ir)))
       ENDDO
    ELSE
       !$omp parallel do private(IR)
       DO ir=1,llr1
          psic(ir)=vpotr(ir)*(REAL(psia(ir))+uimag*AIMAG(psib(ir)))
       ENDDO
    ENDIF
    CALL fwfftn(psic,.TRUE.,parai%allgrp)
    !$omp parallel do private(IG,FP,FM) reduction(+:GHFX)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       ghfx=ghfx+ABS(CMPLX(REAL(fp),AIMAG(fm),kind=real_8))&
            +ABS(CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
       c2b(ig)=c2b(ig)-CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       c2a(ig)=c2a(ig)-CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
    ENDDO
    ghfx=ghfx/2.0_real_8
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxab
  ! ==================================================================
  SUBROUTINE hfxab2(ehfx_1,ehfx_2,ghfx_1,ghfx_2,pf1,pf2,psia,psib,&
       vpotg,vpotr,psic,c2a,c2b1,c2b2,ia,ib1,ib2,f,deeq_fnl_hfx,fion,tfor)
    ! ==================================================================
    REAL(real_8),INTENT(OUT)                 :: ehfx_1, ehfx_2, ghfx_1, &
                                                ghfx_2
    REAL(real_8),INTENT(IN)                  :: pf1, pf2
    COMPLEX(real_8),INTENT(IN)               :: psia(llr1), psib(llr1)
    COMPLEX(real_8),INTENT(OUT)              :: vpotg(jhg,2)
    REAL(real_8),INTENT(OUT)                 :: vpotr(llr1,2)
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: psic(:)
    COMPLEX(real_8),INTENT(INOUT)            :: c2a(jgw), c2b1(jgw), c2b2(jgw)
    INTEGER,INTENT(IN)                       :: ia,ib1,ib2
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: f(:)
    REAL(real_8),INTENT(inOUT) __CONTIGUOUS    :: deeq_fnl_hfx(:,:,:)
    LOGICAL,INTENT(IN)                       :: tfor

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxab2'

    COMPLEX(real_8)                          :: fm, fp, geq0_vpotg1, geq0_vpotg2, vp1, vp2
    INTEGER                                  :: ig, ir, isub, nstates(2,2)

    CALL tiset(procedureN,isub)
    ehfx_1=0.0_real_8
    ehfx_2=0.0_real_8
    ghfx_1=0.0_real_8
    ghfx_2=0.0_real_8
    DO ir=1,llr1
       psic(ir)=REAL(psia(ir))*psib(ir)
    ENDDO
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)
    IF(pslo_com%tivan) THEN
       nstates(1,1)=ia
       nstates(2,1)=ib1
       nstates(1,2)=ia
       nstates(2,2)=ib2
       CALL rhov(nstates,psi=vpotg,hfx=.TRUE.)
    ELSE
       !$omp workshare
       vpotg(:,:)=cmplx(0.0_real_8,0.0_real_8)
       !$omp end workshare
    END IF
    IF(geq0)then
       geq0_vpotg1=vpotg(1,1)
       geq0_vpotg2=vpotg(1,2)
    END IF

    !$omp parallel do private(IG,FP,FM,vp1,vp2) &
    !$omp  reduction(+:EHFX_1,EHFX_2)
    DO ig=1,jhg
       fp=psic(nzff(ig))+psic(inzf(ig))
       fm=psic(nzff(ig))-psic(inzf(ig))
       vp1=vpotg(ig,1)
       vp2=vpotg(ig,2)
       vpotg(ig,1)=-pf1*scgx(ig)*(0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)+vp1)
       vpotg(ig,2)=-pf2*scgx(ig)*(0.5_real_8*CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)+vp2)
       ehfx_1=ehfx_1+REAL(2._real_8*vpotg(ig,1)*(CMPLX(REAL(fp),-AIMAG(fm),kind=real_8)&
            +CONJG(vp1)*2.0_real_8))
       ehfx_2=ehfx_2+REAL(2._real_8*vpotg(ig,2)*(CMPLX(AIMAG(fp),REAL(fm),kind=real_8)&
            +CONJG(vp2)*2.0_real_8))
    ENDDO
    IF (geq0) THEN
       fp=psic(nzff(1))+psic(inzf(1))+geq0_vpotg1
       fm=psic(nzff(1))-psic(inzf(1))+geq0_vpotg2
       ehfx_1=ehfx_1-REAL(vpotg(1,1)*CMPLX(REAL(fp),-AIMAG(fm),kind=real_8))
       ehfx_2=ehfx_2-REAL(vpotg(1,2)*CMPLX(AIMAG(fp),REAL(fm),kind=real_8))
    ENDIF
    IF(pslo_com%tivan)CALL newd(deeq_fnl_hfx,f,vpotg,fion,tfor,nstates,.TRUE.)
    CALL zeroing(psic)!,maxfftn)
    !$omp parallel do private(IG)
    DO ig=1,jhg
       psic(nzff(ig))=vpotg(ig,1)+uimag*vpotg(ig,2)
       psic(inzf(ig))=CONJG(vpotg(ig,1))+uimag*CONJG(vpotg(ig,2))
    ENDDO
    IF (geq0) psic(nzff(1))=vpotg(1,1)+uimag*vpotg(1,2)
    CALL invfftn(psic,.FALSE.,parai%allgrp)
    !$omp parallel do private(IR)
    DO ir=1,llr1
       vpotr(ir,1)=REAL(psic(ir))
       vpotr(ir,2)=AIMAG(psic(ir))
    ENDDO
    !$omp parallel do private(IR)
    DO ir=1,llr1
       psic(ir)=vpotr(ir,1)*(REAL(psia(ir))+uimag*REAL(psib(ir)))
    ENDDO
    CALL fwfftn(psic,.TRUE.,parai%allgrp)
    !$omp parallel do private(IG,FP,FM) reduction(+:GHFX_1)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       ghfx_1=ghfx_1+ABS(CMPLX(REAL(fp),AIMAG(fm),kind=real_8))&
            +ABS(CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
       c2b1(ig)=c2b1(ig)-CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       c2a(ig)=c2a(ig)-CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
    ENDDO
    ghfx_1=ghfx_1/2.0_real_8
    !$omp parallel do private(IR)
    DO ir=1,llr1
       psic(ir)=vpotr(ir,2)*(REAL(psia(ir))+uimag*AIMAG(psib(ir)))
    ENDDO
    CALL fwfftn(psic,.TRUE.,parai%allgrp)
    !$omp parallel do private(IG,FP,FM) reduction(+:GHFX_2)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       ghfx_2=ghfx_2+ABS(CMPLX(REAL(fp),AIMAG(fm),kind=real_8))&
            +ABS(CMPLX(AIMAG(fp),-REAL(fm),kind=real_8))
       c2b2(ig)=c2b2(ig)-CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
       c2a(ig)=c2a(ig)-CMPLX(AIMAG(fp),-REAL(fm),kind=real_8)
    ENDDO
    ghfx_2=ghfx_2/2.0_real_8
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxab2
  ! ==================================================================
  SUBROUTINE hfxaa(ehfx,ghfx,pf,psia,vpotg,vpotr,psic,c2a,ia,f,deeq_fnl_hfx,fion,tfor)
    ! ==================================================================
    REAL(real_8),INTENT(OUT)                 :: ehfx, ghfx
    REAL(real_8),INTENT(IN)                  :: pf
    COMPLEX(real_8),INTENT(IN) __CONTIGUOUS  :: psia(:)
    COMPLEX(real_8),INTENT(OUT)              :: vpotg(jhg)
    REAL(real_8)                             :: vpotr(*)
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    REAL(real_8),INTENT(IN) __CONTIGUOUS     :: f(:)
    REAL(real_8),INTENT(inOUT) __CONTIGUOUS    :: deeq_fnl_hfx(:,:,:)
    INTEGER,INTENT(IN)                       :: ia
    LOGICAL,INTENT(IN)                       :: tfor
    COMPLEX(real_8),INTENT(OUT) __CONTIGUOUS :: psic(:)
    COMPLEX(real_8),INTENT(INOUT)            :: c2a(*)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxaa'

    COMPLEX(real_8)                          :: fm, fp, geq0_vpotg
    INTEGER                                  :: ig, ir, isub, nstates(2,1)

    CALL tiset(procedureN,isub)
    ehfx=0.0_real_8
    ghfx=0.0_real_8
    !$omp parallel do private(IR)
    DO ir=1,llr1
       psic(ir)=REAL(psia(ir))*REAL(psia(ir))
    ENDDO
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)

    IF(pslo_com%tivan)THEN
       nstates(1,1)=ia
       nstates(2,1)=ia
       CALL rhov(nstates,psi=vpotg,hfx=.TRUE.)
    ELSE
       !$omp workshare
       vpotg=CMPLX(0.0_real_8,0.0_real_8)
       !$omp end workshare
    END IF
    IF(geq0) geq0_vpotg=vpotg(1)
    !$omp parallel do private(IG,FP) &
    !$omp  reduction(+:EHFX)
    DO ig=1,jhg
       fp=psic(nzff(ig))+vpotg(ig)
       vpotg(ig)=-pf*scgx(ig)*fp
       ehfx=ehfx+REAL(2._real_8*vpotg(ig)*CONJG(fp))
    ENDDO
    IF(geq0) ehfx=ehfx-REAL(vpotg(1)*CONJG(psic(nzff(1))+geq0_vpotg))
    IF(pslo_com%tivan) CALL newd(deeq_fnl_hfx,f,vpotg,fion,tfor,nstates,.TRUE.)
    CALL zeroing(psic)!,maxfftn)
    !$omp parallel do private(IG)
    DO ig=1,jhg
       psic(nzff(ig))=vpotg(ig)
       psic(inzf(ig))=CONJG(vpotg(ig))
    ENDDO
    IF (geq0) psic(nzff(1))=vpotg(1)
    CALL invfftn(psic,.FALSE.,parai%allgrp)
    !$omp parallel do private(IR)
    DO ir=1,llr1
       vpotr(ir)=REAL(psic(ir))
    ENDDO
    !$omp parallel do private(IR)
    DO ir=1,llr1
       psic(ir)=vpotr(ir)*REAL(psia(ir))
    ENDDO
    CALL fwfftn(psic,.TRUE.,parai%allgrp)
    !$omp parallel do private(IG,FP,FM) reduction(+:GHFX)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       ghfx=ghfx+ABS(CMPLX(REAL(fp),AIMAG(fm),kind=real_8))
       c2a(ig)=c2a(ig)-CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
    ENDDO
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxaa
  ! ==================================================================
  SUBROUTINE hfxpsi_old(c0,cpsi,c2,f,sign,psia,nstate,norb)
    ! ==================================================================
    ! == HARTREE-FOCK EXCHANGE                                        ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: c0(:,:)
    REAL(real_8)                             :: f(:), sign
    COMPLEX(real_8)                          :: psia(:)
    INTEGER                                  :: nstate, norb
    COMPLEX(real_8)                          :: c2(ncpw%ngw,norb), &
                                                cpsi(ncpw%ngw,norb)

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxpsi_old'

    COMPLEX(real_8), ALLOCATABLE             :: psib(:), psic(:), vpotg(:)
    INTEGER                                  :: i, ia, ib, ib1, ib2, ierr, &
                                                ig, isub
    REAL(real_8)                             :: pfl, pfx
    REAL(real_8), ALLOCATABLE                :: vpotr(:)

    IF (func1%mhfx.EQ.0) RETURN
    IF (func1%mhfx.EQ.2) CALL stopgm('HFXPSI_OLD',&
         'HARTREE METHOD NOT POSSIBLE',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tiset('    HFXPSI_OLD',isub)
    IF (tkpts%tkpnt) CALL stopgm('HFXPSI_OLD','K-POINTS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (lspin2%tlse) CALL stopgm('HFXPSI_OLD','NO LSE POSSIBLE',& 
         __LINE__,__FILE__)
    IF (pslo_com%tivan) CALL stopgm('HFXPSI_OLD','NO VDB PP POSSIBLE',& 
         __LINE__,__FILE__)
    IF (group%nogrp.GT.1) CALL stopgm('HFXPSI_OLD',&
         'TASK GROUPS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    DO i=1,nstate
       IF (f(1).NE.f(i).AND. f(i).GT.1.e-6_real_8)&
            CALL stopgm('HFXPSI_OLD','OCCUPATION NUMBERS CHANGE',& 
            __LINE__,__FILE__)
    ENDDO
    ! 
    CALL setfftn(ipoolhfx)
    ! 
    pfl=0.5_real_8
    IF (cntl%thybrid) pfl=pfl*func3%phfx
    pfl=pfl*sign
    ALLOCATE(psic(maxfftn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpotg(jhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpotr(llr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! HARTREE-FOCK
    ALLOCATE(psib(maxfftn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO ia=1,norb
       CALL zeroing(psia)!,maxfftn)
       !ocl novrec
       DO ig=1,jgw
          psia(nzfs(ig))=cpsi(ig,ia)
          psia(inzs(ig))=CONJG(cpsi(ig,ia))
       ENDDO
       IF (geq0) psia(nzfs(1))=cpsi(1,ia)
       ! Transform the wavefunction to real space
       CALL invfftn(psia,.TRUE.,parai%allgrp)
       DO ib=1,nstate,2
          ib1=ib
          ib2=ib+1
          IF (f(ib1).LT.1.e-6_real_8) ib1=0
          IF (ib2.GT.nstate) THEN
             ib2=0
          ELSEIF (f(ib2).LT.1.e-6_real_8) THEN
             ib2=0
          ENDIF
          IF (ib1.NE.0 .OR. ib2.NE.0) THEN
             IF (ib1.EQ.0) THEN
                CALL zeroing(psib)!,maxfftn)
                !ocl novrec
                DO ig=1,jgw
                   psib(nzfs(ig))=uimag*c0(ig,ib2)
                   psib(inzs(ig))=uimag*CONJG(c0(ig,ib2))
                ENDDO
                IF (geq0) psib(nzfs(1))=uimag*c0(1,ib2)
             ELSEIF (ib2.EQ.0) THEN
                CALL zeroing(psib)!,maxfftn)
                !ocl novrec
                DO ig=1,jgw
                   psib(nzfs(ig))=c0(ig,ib1)
                   psib(inzs(ig))=CONJG(c0(ig,ib1))
                ENDDO
                IF (geq0) psib(nzfs(1))=c0(1,ib1)
             ELSE
                CALL zeroing(psib)!,maxfftn)
                !ocl novrec
                DO ig=1,jgw
                   psib(nzfs(ig))=c0(ig,ib1)+uimag*c0(ig,ib2)
                   psib(inzs(ig))=CONJG(c0(ig,ib1))+&
                        uimag*CONJG(c0(ig,ib2))
                ENDDO
                IF (geq0) psib(nzfs(1))=c0(1,ib1)+uimag*c0(1,ib2)
             ENDIF
             ! Transform the wavefunction to real space
             CALL invfftn(psib,.TRUE.,parai%allgrp)
             ! Terms IA*IB1
             IF (ib1.NE.0) THEN
                pfx=pfl*f(ib1)
                CALL hfxpb(pfx,psia,psib,1,vpotg,vpotr,psic,c2(1,ia))
             ENDIF
             ! Terms IA*IB2
             IF (ib2.NE.0) THEN
                pfx=pfl*f(ib2)
                CALL hfxpb(pfx,psia,psib,2,vpotg,vpotr,psic,c2(1,ia))
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    DEALLOCATE(psib,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psic,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpotg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpotr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('    HFXPSI_OLD',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxpsi_old
  ! ==================================================================
  SUBROUTINE hfxpb(pf,psia,psib,iran,vpotg,vpotr,psic,c2a)
    ! ==================================================================
    REAL(real_8)                             :: pf
    COMPLEX(real_8)                          :: psia(llr1), psib(llr1)
    INTEGER                                  :: iran
    COMPLEX(real_8)                          :: vpotg(jhg)
    REAL(real_8)                             :: vpotr(llr1)
    COMPLEX(real_8)                          :: psic(llr1), c2a(jgw)

    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: ig, ir

    IF (iran.EQ.1) THEN
       DO ir=1,llr1
          psic(ir)=REAL(psia(ir))*REAL(psib(ir))
       ENDDO
    ELSE
       DO ir=1,llr1
          psic(ir)=REAL(psia(ir))*AIMAG(psib(ir))
       ENDDO
    ENDIF
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)
    DO ig=1,jhg
       vpotg(ig)=-pf*scgx(ig)*psic(nzff(ig))
    ENDDO
    CALL zeroing(psic)!,maxfftn)
    !ocl novrec
    DO ig=1,jhg
       psic(nzff(ig))=vpotg(ig)
       psic(inzf(ig))=CONJG(vpotg(ig))
    ENDDO
    IF (geq0) psic(nzff(1))=vpotg(1)
    CALL invfftn(psic,.FALSE.,parai%allgrp)
    DO ir=1,llr1
       vpotr(ir)=REAL(psic(ir))
    ENDDO
    IF (iran.EQ.1) THEN
       DO ir=1,llr1
          psic(ir)=vpotr(ir)*REAL(psib(ir))
       ENDDO
    ELSE
       DO ir=1,llr1
          psic(ir)=vpotr(ir)*AIMAG(psib(ir))
       ENDDO
    ENDIF
    CALL fwfftn(psic,.TRUE.,parai%allgrp)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       c2a(ig)=c2a(ig)+0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxpb
  ! ==================================================================
  SUBROUTINE hfxrpa_old(c0,c1,c2,psia,nstate,tcis)
    ! ==================================================================
    ! == HARTREE-FOCK EXCHANGE CONTRIBUTION TO RPA                    ==
    ! ==--------------------------------------------------------------==
    COMPLEX(real_8)                          :: psia(:)
    INTEGER                                  :: nstate
    COMPLEX(real_8)                          :: c2(ncpw%ngw,nstate), &
                                                c1(ncpw%ngw,nstate), &
                                                c0(ncpw%ngw,nstate)
    LOGICAL                                  :: tcis

    CHARACTER(*), PARAMETER                  :: procedureN = 'hfxrpa_old'

    COMPLEX(real_8), ALLOCATABLE             :: psib(:), psic(:), vpotg(:)
    INTEGER                                  :: ia, ib, ib1, ib2, ierr, ig, &
                                                isub, nlower, nupper
    REAL(real_8)                             :: pfl
    REAL(real_8), ALLOCATABLE                :: vpotr(:)

! Variables
! (NHG)
! ==--------------------------------------------------------------==

    IF (func1%mhfx.EQ.0) RETURN
    IF (func1%mhfx.EQ.2) CALL stopgm('HFXRPA','HARTREE METHOD NOT POSSIBLE',& 
         __LINE__,__FILE__)
    ! ==--------------------------------------------------------------==
    CALL tiset('    HFXRPA',isub)
    IF (tkpts%tkpnt) CALL stopgm('HFXRPA','K-POINTS NOT IMPLEMENTED',& 
         __LINE__,__FILE__)
    IF (lspin2%tlse) CALL stopgm('HFXRPA','NO LSE POSSIBLE',& 
         __LINE__,__FILE__)
    IF (pslo_com%tivan) CALL stopgm('HFXRPA','NO VDB PP POSSIBLE',& 
         __LINE__,__FILE__)
    IF (group%nogrp.GT.1) CALL stopgm('HFXRPA','TASK GROUPS NOT IMPLEMENTED'&
         ,& 
         __LINE__,__FILE__)
    ! 
    CALL setfftn(ipoolhfx)
    ! 
    pfl=0.5_real_8
    IF (cntl%tlsd) pfl=1.0_real_8
    pfl=pfl*2.0_real_8
    IF (.NOT.tcis) pfl=2._real_8*pfl
    IF (cntl%thybrid) pfl=pfl*func3%phfx
    ALLOCATE(psic(maxfftn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpotg(jhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpotr(llr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ! HARTREE-FOCK
    ALLOCATE(psib(maxfftn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    DO ia=1,nstate
       CALL zeroing(psia)!,maxfftn)
       !ocl novrec
       DO ig=1,jgw
          psia(nzfs(ig))=c1(ig,ia)+uimag*c0(ig,ia)
          psia(inzs(ig))=CONJG(c1(ig,ia))+&
               uimag*CONJG(c0(ig,ia))
       ENDDO
       IF (geq0) psia(nzfs(1))=c1(1,ia)+uimag*c0(1,ia)
       ! Transform the wavefunctions to real space
       CALL invfftn(psia,.TRUE.,parai%allgrp)
       IF (cntl%tlsd) THEN
          IF (ia.LE.spin_mod%nsup) THEN
             nupper=spin_mod%nsup
             nlower=1
          ELSE
             nupper=nstate
             nlower=spin_mod%nsup+1
          ENDIF
       ELSE
          nupper=nstate
          nlower=1
       ENDIF
       DO ib=nlower,nupper,2
          CALL zeroing(psib)!,maxfftn)
          ib1=ib
          ib2=ib+1
          IF (ib2.GT.nupper) ib2=0
          IF (ib2.EQ.0) THEN
             !ocl novrec
             DO ig=1,jgw
                psib(nzfs(ig))=c0(ig,ib1)
                psib(inzs(ig))=CONJG(c0(ig,ib1))
             ENDDO
             IF (geq0) psib(nzfs(1))=c0(1,ib1)
          ELSE
             !ocl novrec
             DO ig=1,jgw
                psib(nzfs(ig))=c0(ig,ib1)+uimag*c0(ig,ib2)
                psib(inzs(ig))=CONJG(c0(ig,ib1))+&
                     uimag*CONJG(c0(ig,ib2))
             ENDDO
             IF (geq0) psib(nzfs(1))=c0(1,ib1)+uimag*c0(1,ib2)
          ENDIF
          ! Transform the wavefunction to real space
          CALL invfftn(psib,.TRUE.,parai%allgrp)
          ! 
          CALL hfxrpav(pfl,psia,psib,1,vpotg,vpotr,psic,c2(1,ib1))
          ! 
          IF (ib2.NE.0) THEN
             CALL hfxrpav(pfl,psia,psib,2,vpotg,vpotr,psic,c2(1,ib2))
          ENDIF
       ENDDO
    ENDDO
    DEALLOCATE(psib,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(psic,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpotg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpotr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    CALL tihalt('    HFXRPA',isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxrpa_old
  ! ==================================================================
  SUBROUTINE hfxrpav(pf,psia,psib,iran,vpotg,vpotr,psic,c2a)
    ! ==================================================================
    REAL(real_8)                             :: pf
    COMPLEX(real_8)                          :: psia(llr1), psib(llr1)
    INTEGER                                  :: iran
    COMPLEX(real_8)                          :: vpotg(jhg)
    REAL(real_8)                             :: vpotr(llr1)
    COMPLEX(real_8)                          :: psic(llr1), c2a(jgw)

    COMPLEX(real_8)                          :: fm, fp
    INTEGER                                  :: ig, ir

    IF (iran.EQ.1) THEN
       DO ir=1,llr1
          psic(ir)=AIMAG(psia(ir))*REAL(psib(ir))
       ENDDO
    ELSE
       DO ir=1,llr1
          psic(ir)=AIMAG(psia(ir))*AIMAG(psib(ir))
       ENDDO
    ENDIF
    CALL dscal(2*llr1,1._real_8/parm%omega,psic,1)
    CALL fwfftn(psic,.FALSE.,parai%allgrp)
    DO ig=1,jhg
       vpotg(ig)=-pf*scgx(ig)*psic(nzff(ig))
    ENDDO
    CALL zeroing(psic)!,maxfftn)
    !ocl novrec
    DO ig=1,jhg
       psic(nzff(ig))=vpotg(ig)
       psic(inzf(ig))=CONJG(vpotg(ig))
    ENDDO
    IF (geq0) psic(nzff(1))=vpotg(1)
    CALL invfftn(psic,.FALSE.,parai%allgrp)
    DO ir=1,llr1
       vpotr(ir)=REAL(psic(ir))
    ENDDO
    DO ir=1,llr1
       psic(ir)=vpotr(ir)*REAL(psia(ir))
    ENDDO
    CALL fwfftn(psic,.TRUE.,parai%allgrp)
    DO ig=1,jgw
       fp=psic(nzfs(ig))+psic(inzs(ig))
       fm=psic(nzfs(ig))-psic(inzs(ig))
       c2a(ig)=c2a(ig)-0.5_real_8*CMPLX(REAL(fp),AIMAG(fm),kind=real_8)
    ENDDO
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE hfxrpav
  ! ==================================================================
  INTEGER FUNCTION dminloc(n,x,incx)
    IMPLICIT NONE
    INTEGER :: n,incx
    REAL(real_8) :: x(*)
    INTEGER :: i,ix
    REAL(real_8) :: dmin
    dminloc = 0
    IF (n.LT.1 .OR. incx.LE.0) RETURN
    ix = 1
    dmin = ABS(x(ix))
    dminloc = ix
    ix = ix + incx
    DO i = 2,n
       IF (ABS(x(ix)).LT.dmin) THEN
          dminloc = i
          dmin = ABS(x(ix))
       ENDIF
       ix = ix + incx
    ENDDO
  END FUNCTION dminloc

  INTEGER FUNCTION dmaxloc(n,x,incx)
    IMPLICIT NONE
    INTEGER :: n,incx
    REAL(real_8) :: x(*)
    INTEGER :: i,ix
    REAL(real_8) :: dmax
    dmaxloc = 0
    IF (n.LT.1 .OR. incx.LE.0) RETURN
    ix = 1
    dmax = ABS(x(ix))
    dmaxloc = ix
    ix = ix + incx
    DO i = 2,n
       IF (ABS(x(ix)).GT.dmax) THEN
          dmaxloc = i
          dmax = ABS(x(ix))
       ENDIF
       ix = ix + incx
    ENDDO
  END FUNCTION dmaxloc
!==========================================================!
! HERE SCDM STUFF IS INCLUDED                              !
! HFX_SCDM ROUTINE CALCULATES HFX PART USING               !
!  SCMD BASED SCREENING                                    !
!----------------------------------------------------------!

!   ==================================================================
!   ==================================================================
    SUBROUTINE HFX_SCDM(C0,C2,F,PSIA,NSTATE,EHFX,VHFX,redist_c2,deeq_fnl_hfx,fion,tfor)
!   ==================================================================
!   == HARTREE-FOCK EXCHANGE                                        ==
!   ==--------------------------------------------------------------==
!   Arguments
    COMPLEX(real_8)                          :: c0(:,:), c2(:,:)
    REAL(real_8)                             :: f(:)
    COMPLEX(real_8)                          :: psia(:)
    INTEGER                                  :: nstate
    REAL(real_8)                             :: ehfx, vhfx
    LOGICAL                                  :: redist_c2
    REAL(real_8),INTENT(INOUT) __CONTIGUOUS  :: fion(:,:,:)
    REAL(real_8),INTENT(inOUT) __CONTIGUOUS    :: deeq_fnl_hfx(:,:,:)
    LOGICAL,INTENT(IN)                       :: tfor
!============================================================
!     Variables
!------------------------------------------------------------
    CHARACTER(*), PARAMETER                  :: procedureN = 'hfx_scdm'
    COMPLEX(real_8), PARAMETER :: zone = (1.0_real_8,0.0_real_8), &
      zzero = (0.0_real_8,0.0_real_8)

    COMPLEX(real_8), ALLOCATABLE             :: psic(:), vpotg(:)
    COMPLEX(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: C2_hfx, psib
    INTEGER :: geo_iter, i, ia, ib1, ib2, ibin, ibuf, id_states(2,2), ierr, &
      ig, int_i, int_ij, int_j, ispin, isub, isub2, isub3, isub4, isub5, &
      isub6, iwf, j, jb1, jb2, my_pcol, my_prow, n_int, nbr_states(2), &
      npcols, nprows, nspins, nst, nstates(2), st_offst
    INTEGER(int_8) :: nbr_int_skiped_1, nbr_int_skiped_2, nbr_integrals, &
      nbr_rwfn_precomputed, nbr_rwfn_recomputed
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: bin_vals, new_vals, &
                                                psi_in_core
    INTEGER, SAVE                            :: int_count = 0, &
                                                prev_geo_id = -HUGE(0)
    LOGICAL                                  :: init_ints, no_more_states
    REAL(real_8) :: bin_max, bin_range, dab, dt, dt_tot, EHFX_loc, &
      EHFX_loc_1, EHFX_loc_2, GHFX_loc, GHFX_loc_1, GHFX_loc_2, max_dab, &
      old_DWFC, pfl, pfx, pfx1, pfx2, t1, t1_tot, t2, t2_tot
    REAL(real_8), ALLOCATABLE                :: vpotr(:)
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :)                        :: acc_vals, max_vals
    REAL(real_8), ALLOCATABLE, &
      DIMENSION(:, :), SAVE                  :: int_vals
!-----------------------------------------------------------------
    REAL(real_8),DIMENSION(:,:),ALLOCATABLE  :: U
    LOGICAL                        :: NEGLECT
!    REAL(real_8)                   :: DIST
    REAL(real_8), ALLOCATABLE      :: RHOIJ(:),PSI(:,:)
    INTEGER                        :: K, KA, KB, IR, ID
    COMPLEX(real_8), ALLOCATABLE   :: C2_HFXKS(:,:)
!C  ====--------------------------------------------------------------==
!------------------------------------------------------------------------!
!   sagar  ACE stuff
    INTEGER                                    :: info
    COMPLEX(real_8), ALLOCATABLE, &
        DIMENSION(:, :)                        :: cmexx
    REAL(real_8), ALLOCATABLE, &
        DIMENSION(:, :)                        :: mexx
!========================================================================!

    ehfx=0._real_8
    vhfx=0._real_8
    IF (func1%mhfx.EQ.0) RETURN
    !
    CALL TISET(procedureN,ISUB)
!-----------------------------------------------------------------------
    !  write(6,*)"here",parai%me,scdm_cutoff,hfx_scdm_status
    !scdm_dist=.false.  !todo
    !scdm_norm=.true.
    !scdm_max=.false.
    !SCDM_CUTOFF=1.D-8
!    do k=1,scdm_number
!      write(6,*)"here",parai%me,scdm_cut(k),scdm_num(k)
!    end do
!    write(6,*)"here hfx_scdm",parai%me,scdm_dist,scdm_norm, & 
!                              scdm_max,hfx_scdm_l,scdm_number
!    write(6,*)"here hfx_scdm_lang",parai%me,lang_dyn,gamma,t_bath
!=======================================================================
    ALLOCATE(C2_hfxks(ncpw%ngw,nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
         __LINE__,__FILE__)

    CALL zero_wfn(jgw,nstate,C2_hfxks,ncpw%ngw)

    ALLOCATE(PSI(LLR1,nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
         __LINE__,__FILE__)

    ALLOCATE(RHOIJ(LLR1),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
         __LINE__,__FILE__)

    ! ==--------------------------------------------------------------==
    ! sagar
    IF(USE_ACE)THEN
      ALLOCATE(cmexx(nstate,nstate),mexx(nstate,nstate),&
           stat=ierr)
      IF (ierr /= 0) CALL stopgm( procedureN, 'Allocation problem' ,&
           __LINE__,__FILE__)
    ENDIF
    !
    ! ==--------------------------------------------------------------==

!    ALLOCATE(SCDM_CENTER(3,nstate),stat=ierr)
!    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
!         __LINE__,__FILE__)

!C     ==--------------------------------------------------------------==
    ! decide if reinit intergrals
    ! need to use the *right* geo iter counter (in the cpmd spirit)!
    init_ints=.FALSE.
    geo_iter=0
    IF(.NOT.cntl%wfopt) geo_iter=infi
    IF (prev_geo_id/=geo_iter) THEN
       prev_geo_id=geo_iter
       init_ints=.TRUE.
    ENDIF
!----------------------------------------
    t1_tot=m_cputime()

    IF(hfxc3%twscr.OR.cntl%krwfn)CALL STOPGM(procedureN, &
          'TWSCR OR KRWFN NOT IMPLEMENTED WITH SCDM METHOD', &
         __LINE__,__FILE__)
    !IF(parai%cp_nogrp > 1)CALL STOPGM(procedureN, &
    !      'CP_NOGRP>1 NOT IMPLEMENTED WITH ACE+SCDM METHOD', &
    !     __LINE__,__FILE__)
!-------------------------------------------------------------

    IF (cntl%cdft.AND.parai%cp_nogrp>1) CALL stopgm(procedureN,&
         'CDFT AND CP_NOGRP>1 NOT YET IMPLEMENTED',&
         __LINE__,__FILE__)
    IF (tkpts%tkpnt) CALL stopgm(procedureN,'K-POINTS NOT IMPLEMENTED',&
         __LINE__,__FILE__)
    IF (lspin2%tlse) CALL stopgm(procedureN,'NO LSE POSSIBLE',&
         __LINE__,__FILE__)
    IF (group%nogrp.GT.1) CALL stopgm(procedureN,&
         'TASK GROUPS NOT IMPLEMENTED',&
         __LINE__,__FILE__)

!===============================================================

!
    CALL SETFFTN(IPOOLHFX)
!
!     Prepare the indeces for parallelization over the groups
!C
    CALL proc_to_grid2d(parai%cp_nogrp,parai%cp_inter_me,nprows,npcols,&
         my_prow,my_pcol)

!C
!C     Need to allocate a buffer to accumulate the x part of C2
!C

    ALLOCATE(C2_hfx(ncpw%ngw,nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
         __LINE__,__FILE__)

    CALL zero_wfn(jgw,nstate,C2_hfx,ncpw%ngw)

!C
!C     For the moment we synchronize the C0 and C2 if the cp groups are used
!C     That should be done else where
!   TODO is it needed?
    !IF (parai%cp_nogrp.GT.1) THEN
    !   CALL tiset(procedureN//'_a',isub5)
    !   CALL mp_bcast(c0,ncpw%ngw*nstate,0,parai%cp_inter_grp)
    !   CALL mp_bcast(c2,ncpw%ngw*nstate,0,parai%cp_inter_grp)
    !   CALL tihalt(procedureN//'_a',isub5)
    !ENDIF

!c
!c     randomize the states to improve load balance
!c
!TODO

!C
    t1=m_cputime()

!      CALL check_occupation(F,NSTATE,TLSD)

!C
!C     set the spins stuff
!C
    nstates = 0
    IF (cntl%tlsd) THEN
       nspins = 2
       nstates(1) = spin_mod%nsup
       nstates(2) = spin_mod%nsdown
       pfl = 0.5_real_8
    ELSE
       nspins = 1
       nstates(1) = nstate
       pfl = 0.25_real_8
    ENDIF

    IF (cntl%thybrid) pfl=pfl*func3%phfx
    ALLOCATE(psic(maxfftn),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpotg(2*jhg),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(vpotr(2*llr1),STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
         __LINE__,__FILE__)
    ALLOCATE(psi_in_core(0:nstate),stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
         __LINE__,__FILE__)
    psi_in_core(:) = 0


     IF (func1%mhfx.EQ.1) THEN
       ! HARTREE-FOCK
       ALLOCATE(psib(maxfftn,2),stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Allocation problem',&
            __LINE__,__FILE__)
       CALL zeroing(psib)!,2*maxfftn)
!---------------------------------------------------------------------==
       CALL G_TO_R_WF(C0,NSTATE,PSI)        !SAGAR HACK TODO
!C     ==--------------------------------------------------------------==
       CALL TISET(procedureN//'_outer',ISUB2)
!C
!C     loop over the spin
!C
       id_states(:,:) = 0
       nbr_states(:) = 0
       nbr_int_skiped_1 = 0; nbr_int_skiped_2 = 0;
       nbr_integrals = 0
       nbr_rwfn_precomputed = 0
       nbr_rwfn_recomputed = 0
       DO ispin = 1,nspins

           nst = nstates(ispin)
           st_offst = 0
           IF(ispin.EQ.2) st_offst = nstates(1)
!=======================================================================
           !
           ALLOCATE(U(NST,NST),STAT=ierr)
           IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
              __LINE__,__FILE__)
           !
           CALL LOCALIZE_SCDM(PSI(1,1+ST_OFFST),NST,U)  ! TODO
           !
!==========================================================
!           IF(SCDM_DIST)THEN
!             DO ID=1+ST_OFFST,NST+ST_OFFST
!               CALL FIND_CENTER(PSI(1,ID), &
!               SCDM_CENTER(1,ID),SCDM_CENTER(2,ID),SCDM_CENTER(3,ID))
!             END DO
!           END IF
!----------------------------------------------------------
           !
           CALL DSCAL(NST*LLR1,DSQRT(DFLOAT(spar%NR1S*spar%NR2S*spar%NR3S)), &
               PSI(1,1+ST_OFFST),1)
           !
!C     ==--------------------------------------------------------------==

        DO i = 1,part_1d_nbr_elems(nst,my_prow,nprows)
           ia = part_1d_get_elem(i,my_prow,nprows)+st_offst
           IF(F(IA).LT.1.e-6_real_8) CYCLE
!======================================================================
!          SAGAR HACK
!----------------------------------------------------------------------
           DO IR=1,LLR1
             PSIA(IR)=DCMPLX(PSI(IR,IA),0.D0)
           ENDDO
!======================================================================
!C     Handle the diagonal terms
!C     Only when needed....
             IF (part_1d_elem_to_proc(ia-st_offst,nprows).EQ.my_prow.AND.&
                  part_1d_elem_to_proc(ia-st_offst,npcols).EQ.my_pcol) THEN
                IF (parai%me.EQ.0)nbr_integrals=nbr_integrals+1
                pfx=pfl*f(ia)*f(ia)
                CALL hfxaa(EHFX_loc,GHFX_loc,pfx,psia,vpotg,vpotr,psic,&
                     C2_hfx(1,ia),ia,f,deeq_fnl_hfx,fion,tfor)
                ehfx=ehfx+EHFX_loc
             ENDIF
!C     Now for the other terms IA-IB
!C     ==--------------------------------------------------------------==
           CALL TISET(procedureN//'_inner',ISUB3)
           j = 1
           DO
              no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol, &
                   npcols)
              IF(no_more_states) EXIT

!c
!c     find out the correct ib's
              ib1 = 0; ib2 = 0;
              DO
!c     set id1
                 IF(ib1.EQ.0) THEN
                    ib1 = part_1d_get_elem(j,my_pcol,npcols)+st_offst
                    j = j + 1
                 ENDIF
                 IF(ia.EQ.ib1) ib1 = 0
                 IF(.NOT.part_1d_symm_holds_pair(ia,ib1)) ib1 = 0
                 IF(ib1.NE.0) THEN
                    IF(F(ib1).LT.1.e-6_real_8) ib1 = 0
                 ENDIF
!---------------------------------------------------------------------------------------
                 IF(IB1.NE.0)THEN    !SAGAR HACK
                  !
!                  IF(.NOT.SCDM_DIST)THEN
                   DO K=1,LLR1
                     RHOIJ(K)=PSI(K,IA)*PSI(K,IB1)
                   END DO
!                  END IF
                  !
                  !IF(SCDM_MAX)THEN
                  !  CALL SELECT_PAIR_1(RHOIJ,NEGLECT)
                  !ELSE IF(SCDM_NORM)THEN
                    CALL SELECT_PAIR(RHOIJ,NEGLECT)
                  !ELSE IF(SCDM_DIST)THEN
                  !  !
                  !  CALL SCDM_SEPARATION &
                  !  (SCDM_CENTER(1,IA),SCDM_CENTER(1,IB1),DIST)
                  !  NEGLECT=.FALSE.
                  !  IF(DIST.GT.SCDM_CUTOFF)NEGLECT=.TRUE.
                  !  !
                  !END IF
                  !
                  !IF(NEGLECT.and.(parai%me.eq.0))write(6,*)ia,ib1
                  IF(NEGLECT)IB1=0
                  !
                 END IF
!---------------------------------------------------------------------------------------
                 no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol,  &
                     npcols)
                 IF(no_more_states) EXIT
!c     set ib2
                 IF(ib2.EQ.0) THEN
                    ib2 = part_1d_get_elem(j,my_pcol,npcols)+st_offst
                    j = j + 1
                 ENDIF
                 IF(ia.EQ.ib2) ib2 = 0
                 IF(.NOT.part_1d_symm_holds_pair(ia,ib2)) ib2 = 0
                 IF(ib2.NE.0) THEN
                    IF(F(ib2).LT.1.e-6_real_8) ib2 = 0
                 ENDIF
!---------------------------------------------------------------------------------------
                 IF(IB2.NE.0)THEN      !SAGAR HACK
                  !
                  !IF(.NOT.SCDM_DIST)THEN
                   DO K=1,LLR1
                     RHOIJ(K)=PSI(K,IA)*PSI(K,IB2)
                   END DO
                  !END IF
                  !
                  !IF(SCDM_MAX)THEN
                  !  CALL SELECT_PAIR_1(RHOIJ,NEGLECT)
                  !ELSE IF(SCDM_NORM)THEN
                    CALL SELECT_PAIR(RHOIJ,NEGLECT)
                  !ELSE IF(SCDM_DIST)THEN
                  !  !
                  !  CALL SCDM_SEPARATION  &
                  ! (SCDM_CENTER(1,IA),SCDM_CENTER(1,IB2),DIST)
                  !  NEGLECT=.FALSE.
                  !  IF(DIST.GT.SCDM_CUTOFF)NEGLECT=.TRUE.
                  !  !
                  !END IF
                  !
                  !IF(NEGLECT.and.(parai%me.eq.0))write(6,*)ia,ib2
                  IF(NEGLECT)IB2=0
                  !
                 END IF
!---------------------------------------------------------------------------------------
                 no_more_states = j.GT.part_1d_nbr_elems(nst,my_pcol, &
                      npcols)
                 IF(no_more_states) EXIT
!c     are we set?
                 IF(ib1.NE.0.AND.ib2.NE.0) EXIT
              ENDDO

!c      if(me.eq.0) write(*,*)cp_me,'ia,ib1,ib2',ia,ib1,ib2,no_more_states
!------------------------------------------------------------------------
!             SAGAR HACK
!------------------------------------------------------------------------
              IF(.NOT.(IB1.EQ.0.AND.IB2.EQ.0)) THEN
                  !CALL ZAZZERO(PSIB(1,1),MAXFFTN)
                  CALL zeroing(psib)
                  IF(IB1.NE.0) THEN
                     DO IR=1,LLR1
                       PSIB(IR,1)=PSIB(IR,1)+DCMPLX(PSI(IR,IB1),0.D0)
                     ENDDO
                  ENDIF
                  !
                  IF(IB2.NE.0) THEN
                     DO IR=1,LLR1
                       PSIB(IR,1)=PSIB(IR,1)+UIMAG*  &
                                 DCMPLX(PSI(IR,IB2),0.D0)
                     ENDDO
                  ENDIF
              ENDIF
!=======================================================================
!C     ==--------------------------------------------------------------==
              id_states(1,1) = ib1
              id_states(2,1) = ib2
              nbr_states(1) = 0
              if(ib1.ne.0) nbr_states(1) = nbr_states(1) + 1
              if(ib2.ne.0) nbr_states(1) = nbr_states(1) + 1

              if(nbr_states(1).eq.1) then
                 if(id_states(1,1).ne.0) then
                    if(id_states(1,2).eq.0) then
                       id_states(1,2) = id_states(1,1)
                       CALL copy_re_to_re(LLR1,psib(:,1),psib(:,2))
                    else
                       id_states(2,2) = id_states(1,1)
                       CALL copy_re_to_im(LLR1,psib(:,1),psib(:,2))
                    endif
                    id_states(1,1) = 0
                 else
                    if(id_states(1,2).eq.0) then
                       id_states(1,2) = id_states(2,1)
                       CALL copy_im_to_re(LLR1,psib(:,1),psib(:,2))
                    else
                       id_states(2,2) = id_states(2,1)
                       CALL copy_im_to_im(LLR1,psib(:,1),psib(:,2))
                    endif
                    id_states(2,1) = 0
                 endif
                 nbr_states(1) = nbr_states(1) - 1
                 nbr_states(2) = nbr_states(2) + 1
              endif

!C     ==--------------------------------------------------------------==

              DO ibuf = 1,2
                 IF(.NOT.no_more_states) THEN
                    IF(nbr_states(ibuf).NE.2) CYCLE
                 ENDIF
                 jb1 = id_states(1,ibuf)
                 jb2 = id_states(2,ibuf) 
!C
                 EHFX_loc = 0.0_real_8
                 IF(JB1.NE.0 .AND. JB2.NE.0) THEN
                    nbr_integrals=nbr_integrals+2
                    PFX1=PFL*F(IA)*F(JB1)
                    PFX2=PFL*F(IA)*F(JB2)
                    CALL HFXAB2(EHFX_loc_1,EHFX_loc_2, &
                        GHFX_loc_1,GHFX_loc_2, &
                        PFX1,PFX2,PSIA, &
                        PSIB(1,ibuf), &
                        VPOTG,VPOTR,PSIC,C2_hfx(1,IA), &
                        C2_hfx(1,JB1),C2_hfx(1,JB2),&
                        ia,jb1,jb2,f,deeq_fnl_hfx,fion,tfor) 
                    EHFX_loc = EHFX_loc_1 + EHFX_loc_2
                 ELSEIF(JB1.NE.0 .AND. JB2.EQ.0) THEN
!C     Terms IA*IB1
                    nbr_integrals=nbr_integrals+1
                    PFX=PFL*F(IA)*F(JB1)
                    CALL HFXAB(EHFX_loc,GHFX_loc, &
                        PFX,PSIA,PSIB(1,ibuf),1, &
                        VPOTG,VPOTR,PSIC,C2_hfx(1,IA), &
                        C2_hfx(1,JB1),ia,jb1,f,deeq_fnl_hfx,&
                        fion,tfor) 
                 ELSEIF(JB1.EQ.0 .AND. JB2.NE.0) THEN
!C     Terms IA*IB2
                    nbr_integrals=nbr_integrals+1
                    PFX=PFL*F(IA)*F(JB2)
                    CALL HFXAB(EHFX_loc,GHFX_loc,PFX, &
                        PSIA,PSIB(1,ibuf),2, &
                        VPOTG,VPOTR,PSIC,C2_hfx(1,IA), &
                        C2_hfx(1,JB2),ia,jb2,f,deeq_fnl_hfx,&
                        fion,tfor)
                 ENDIF
                 EHFX=EHFX+EHFX_loc

!C     reset the buffer
                 id_states(:,ibuf) = 0
                 nbr_states(ibuf) = 0 
              ENDDO             !ibuf
           ENDDO                !IB
           CALL TIHALT(procedureN//'_inner',ISUB3)
!C     ==--------------------------------------------------------------==
 101       CONTINUE
       ENDDO                   !IA
!=======================================================================
!       SAGAR HACK todo change for ACE
!-----------------------------------------------------------------------
        DO KA=1+ST_OFFST,ST_OFFST+NST
         DO K=1,ncpw%ngw
          DO KB=1+ST_OFFST,ST_OFFST+NST
            !
            C2_HFXKS(K,KA)=C2_HFXKS(K,KA)+C2_HFX(K,KB) &
                             *U(KA-ST_OFFST,KB-ST_OFFST)
            !
          END DO
         END DO
        END DO
        !
        DEALLOCATE(U,STAT=ierr)
        IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
           __LINE__,__FILE__)
!=======================================================================
        ENDDO                     !ispin
        CALL TIHALT(procedureN//'_outer',ISUB2)
!C     ==--------------------------------------------------------------==

     ELSE
!C     HARTREE
         DO IA=1,NSTATE
            IF(F(IA).GT.1.e-6_real_8) THEN
               !CALL ZAZZERO(PSIA,MAXFFTN)
               CALL zeroing(psia)
               DO IG=1,JGW
                  PSIA(NZFS(IG))=C0(IG,IA)
                  PSIA(INZS(IG))=DCONJG(C0(IG,IA))
               ENDDO
               IF(GEQ0) PSIA(NZFS(1))=C0(1,IA)
!C     Transform the wavefunction to real space
               !CALL INVFFTN(PSIA,.TRUE.)
               CALL invfftn(psia,.TRUE.,parai%allgrp)
               PFX=PFL*F(IA)*F(IA)
               CALL HFXAA(EHFX_loc,GHFX_loc,PFX,PSIA,VPOTG,VPOTR, &
                    PSIC,C2_hfx(1,IA),ia,f,deeq_fnl_hfx,fion,tfor)
               EHFX=EHFX+EHFX_loc
            ENDIF
         ENDDO
     ENDIF
!C
!C     free some memory
!C

    DEALLOCATE(psic,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpotg,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)
    DEALLOCATE(vpotr,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    IF (func1%mhfx.EQ.1) THEN
       DEALLOCATE(psib,stat=ierr)
       IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',&
            __LINE__,__FILE__)
    ENDIF

    t2=m_cputime()
    dt = t2 - t1


!c
!c     redistribute EHFX and C2_hfx over the groups if needed
!c
    IF (parai%cp_nogrp.GT.1) THEN
       ! call mp_sync(CP_GRP)
       CALL tiset(procedureN//'_b',isub6)
       CALL mp_sum(ehfx,parai%cp_inter_grp)
       ! call cp_grp_redist_z(C2_hfx,NGW,NSTATE)
       !CALL mp_sum(C2_hfx,ncpw%ngw*nstate,parai%cp_inter_grp)
       !CALL mp_sum(C2_hfxks,ncpw%ngw*nstate,parai%cp_inter_grp) !TODO
       IF (redist_c2)CALL cp_grp_redist(C2_hfxks,ncpw%ngw,nstate)
       CALL tihalt(procedureN//'_b',isub6)
    ENDIF
!c
!c     add up the hfx contribution to C2
!c
      !CALL add_wfn(JGW,NSTATE,ZONE,C2_hfx,NGW,C2,NGW)  ! SAGAR HACK
!     CALL add_wfn(jgw,nstate,zone,C2_hfx,ncpw%ngw,c2,ncpw%ngw)
      !
      !if(parai%cp_inter_me.eq.0) & !hack fixed for CP_GRP sagar
      CALL add_wfn(JGW,NSTATE,ZONE,C2_HFXKS,ncpw%NGW,C2,ncpw%NGW)
!=======================================================================
!------------------
    IF ((parai%cp_nogrp > 1).and.use_ace) THEN !bug fixed for CP_GRP sagar TODO
      IF (.not.redist_c2) CALL cp_grp_redist(C2_hfxks,ncpw%ngw,nstate)
    ENDIF
    !write(6,*)redist_c2,parai%cp_inter_me
    if(use_ace)call build_ace_proj()  !sagar
!---------------------------

!=-------------------------------------------------------------------------------------
!C     ENERGY
    ehfx=ehfx*parm%omega
    CALL mp_sum(ehfx,parai%allgrp)
    ! POTENTIAL
    DO ia=1,nstate
       vhfx=vhfx+dotp(ncpw%ngw,c0(:,ia),c2(:,ia))
    ENDDO
    CALL mp_sum(vhfx,parai%allgrp)
    DEALLOCATE(psi_in_core,C2_hfx,stat=ierr)
    IF (ierr.NE.0) CALL stopgm(procedureN,'Deallocation problem',&
         __LINE__,__FILE__)

!C
!C     some infos
!C 
    t2_tot=m_cputime()
    dt_tot = t2_tot - t1_tot
!-------------------------------------------
    IF (PRINT_GROUP_INFOS) THEN
       IF (parai%me.EQ.0) THEN
          WRITE(6,'(1X,6(A,I0),A,F0.2,A,F0.2)')&
               procedureN//'| group ',parai%cp_inter_me,&
               ' computed ',nbr_integrals,&
               ' integrals, precomputed ',nbr_rwfn_precomputed,&
               ', recomputed ',nbr_rwfn_recomputed,&
               ' and skipped ',&
               nbr_int_skiped_1,' + ',nbr_int_skiped_2,&
               ', t_loc ',dt,&
               ', t_per_1k_ints ',1.0e3_real_8*dt/(REAL(nbr_integrals,kind=real_8)+1.0e-6_real_8) ! to avoid NANs
          CALL m_flush(6)
       ENDIF

       IF (parai%cp_nogrp.GT.1) THEN
          CALL mp_sum(nbr_integrals,parai%cp_inter_grp)
          CALL mp_sum(nbr_int_skiped_1,parai%cp_inter_grp)
          CALL mp_sum(nbr_int_skiped_2,parai%cp_inter_grp)
          CALL mp_sum(nbr_rwfn_precomputed,parai%cp_inter_grp)
          CALL mp_sum(nbr_rwfn_recomputed,parai%cp_inter_grp)
          CALL mp_max(dt_tot,parai%cp_inter_grp)
          IF (parai%cp_me.EQ.0) THEN
             WRITE(6,'(1X,5(A,I0),A,F0.2)')&
                  procedureN//'| all the groups computed ',&
                  nbr_integrals,' integrals, precomputed ',&
                  nbr_rwfn_precomputed,&
                  ', recomputed ',nbr_rwfn_recomputed,&
                  ' and skipped ',&
                  nbr_int_skiped_1,' + ',nbr_int_skiped_2,&
                  ', t_tot ',dt_tot
             CALL m_flush(6)
          ENDIF
       ENDIF
    ENDIF

!----------------------------------------------------
    ! sagar
    IF(USE_ACE)THEN
      DEALLOCATE(cmexx,mexx,&
           stat=ierr)
      IF (ierr /= 0) CALL stopgm( procedureN, 'Deallocation problem' ,&
           __LINE__,__FILE__)
    ENDIF
    !
!----------------------------------------------------
    DEALLOCATE(C2_HFXKS,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    DEALLOCATE(PSI,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

    DEALLOCATE(RHOIJ,STAT=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
         __LINE__,__FILE__)

!    DEALLOCATE(SCDM_CENTER,STAT=ierr)
!    IF(ierr/=0) CALL stopgm(procedureN,'deallocation problem',&
!         __LINE__,__FILE__)

!-----------==========================---------------
    CALL TIHALT(procedureN,ISUB)
!C     ==--------------------------------------------------------------==
    !RETURN sagar
    CONTAINS

! ==================================================================
    SUBROUTINE build_ace_proj()
!-------------------------------------------------------------------
!   Building the ACE projectors for use in the later part
!===================================================================
        ! copy c2_hfx to xi
        ! |xi> = Vx[phi]|phi> 
        !CALL DCOPY(2*ncpw%ngw*NSTATE,C2_hfx(1,1),1,xi(1,1),1) !TODO
        CALL DCOPY(2*ncpw%ngw*NSTATE,C2_hfxKS(1,1),1,xi(1,1),1)
        ! !   mexx = <phi|Vx[phi]|phi>
        ! OVLAP(NSTATE,A,C1,C2)
        ! COMPUTES THE OVERLAP MATRIX A = < C1 | C2 >
        !CALL AZZERO(mexx,NSTATE*NSTATE)
        CALL zeroing(mexx)
        !call OVLAP(NSTATE,mexx,C0,xi)
        call OVLAP(NSTATE,mexx,C0,xi)
        !CALL MY_SUM_D(mexx,NSTATE*NSTATE,ALLGRP)
        CALL mp_sum(mexx,nstate*nstate,parai%allgrp)
        !
  !----------------------------------------------------------
  !       A=mexx & n=nstate
  !--------------------
          INFO = -1
          CALL DPOTRF( 'L', nstate, mexx, nstate, INFO )
          !
          !IF(Info.NE.0) &
          ! CALL STOPGM('ACE','ERROR IN DPOTRF')
    IF (info /= 0) CALL stopgm(procedureN,'ERROR IN DPOTRF: ACE',&
         __LINE__,__FILE__)
          !
          !CALL errinfo('DPOTRF','Cholesky failed in invchol.',INFO)
          INFO = -1
          CALL DTRTRI( 'L', 'N', nstate, mexx, nstate, INFO )
          !IF(Info.NE.0) &
          !CALL STOPGM('ACE','ERROR IN DTRTRI')
    IF (info /= 0) CALL stopgm(procedureN,'ERROR IN DTRTRI: ACE',&
         __LINE__,__FILE__)
          !CALL errinfo('DTRTRI','inversion failed in invchol.',INFO)
!==========================================================
        !
        cmexx = (1.d0,0.d0)*mexx
        !
        ! |xi> = -One * Vx[phi]|phi> * rmexx^T
        CALL ZTRMM('R','L','C','N',ncpw%ngw,nstate,(1.d0,0.d0), &
                   cmexx,nstate,xi,ncpw%ngw)
        !
!=========================================================
    END SUBROUTINE build_ace_proj
! ==================================================================

    END SUBROUTINE HFX_SCDM
!=======================================================================

!==============================================================================
!=========================================================================!
      SUBROUTINE LOCALIZE_SCDM(PSI,NORB,U)
!=========================================================================!
!     Nisanth N. Nair, nnair@iitk.ac.in (June, 2015)                      !
!                                                                         !
!     Localized Orbitals from Selected Column Density Matrix Computation  !  
!=========================================================================!
!      IMPLICIT NONE
!      INCLUDE 'system.h'
!      INCLUDE 'parac.inc'
!      INCLUDE 'fft.inc'
!      INCLUDE 'irat.inc'
!
      INTEGER               :: NORB
      REAL(real_8)          :: PSI(LLR1,NORB), u(norb,norb) !U(:,:)
!--------------------------------------------------------------------------
!     local variables
      REAL(real_8),ALLOCATABLE :: C(:,:)
      REAL(real_8)         :: DUM
      INTEGER              :: NGRIDS, LSCR, IERR, IORB, JORB, IR, I, K, &
                              MSGLEN, ISUB,MAXLLR1
      !
      REAL(real_8),ALLOCATABLE :: PSIT(:,:),QMAT(:),MYSCR(:), & 
                                  PCMAT(:,:),SMAT(:,:),PSIC(:,:)
      !
      INTEGER, ALLOCATABLE     :: PIMAT(:), PIMATg(:), ME_LLR1(:) 
!======================================================================
      CALL TISET('LOCALIZE_SCDM',ISUB)
!
      IF(NORB.GT.LLR1)CALL STOPGM('LOCALIZE_SCDM', &
       'NORB>LLR1 NOT IMPLEMENTED',__LINE__,__FILE__)
!
!=====================================================
      ALLOCATE(C(NORB,NORB),STAT=IERR)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','allocation problem',&
               __LINE__,__FILE__)
      !
      CALL zeroing(U)
      CALL zeroing(C)
!-----------------------------------------------------
!#ifdef PARALLEL
#ifdef __PARALLEL
      ALLOCATE(ME_LLR1(PARAI%NPROC),STAT=IERR)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','allocation problem',&
               __LINE__,__FILE__)
      !
      call zeroing(ME_LLR1)

      ME_LLR1(PARAI%ME+1)=LLR1

      CALL mp_sum(me_llr1,SIZE(me_llr1),parai%allgrp)

      MAXLLR1=ME_LLR1(1)

      DO I=1,PARAI%NPROC
        MAXLLR1=MAX(MAXLLR1,ME_LLR1(I))
      END DO

      DEALLOCATE(ME_LLR1,stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
         __LINE__,__FILE__)

#else
      MAXLLR1=LLR1
#endif
!
!
      ALLOCATE(PSIT(NORB,MAXLLR1),stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Allocation problem',&
         __LINE__,__FILE__)
  
      ALLOCATE(PIMAT(MAXLLR1),stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Allocation problem',&
         __LINE__,__FILE__)
!------------------------------------------------------------
!c     Normalize real space orbitals
      NGRIDS=spar%NR1S*spar%NR2S*spar%NR3S
      DUM=DSQRT(1.D0/DFLOAT(NGRIDS))
      CALL DSCAL(NORB*LLR1,DUM,PSI,1)
!C     Transpose 
      DO IR=1,LLR1
        DO IORB=1,NORB
          PSIT(IORB,IR)=PSI(IR,IORB)
        END DO
      END DO
!c     Adjusting the unequal distribution of grid points over processors
!c     (uniform distribution is required for SCALAPACK)
      DO IR=LLR1+1,MAXLLR1
        DO IORB=1,NORB
          PSIT(IORB,IR)=0.D0
        END DO
      END DO
!      call cpu_time(final_time)
!      if(parent)
!     & print *, "Initialization: Timing =",final_time-start_time
!c
!c!check psi^2
!c      do iorb=1,norb
!c        do jorb=1,norb
!c          dum=0.d0
!c          do ir=1,llr1
!c             dum=dum+psi(ir,iorb)*psi(ir,jorb)
!cc          end do
!c          CALL MY_SUM_D(DUM,1,CP_GRP)
!c          if(io_parent)write(6,'(2i8,e14.6)')iorb,jorb,dum
!c        end do
!c      end do
!c
!C     Rank Revealing QR factorization
!
      call zeroing(PIMAT)
!      call cpu_time(start_time)
!#ifdef PARALLEL
#ifdef __PARALLEL
      CALL PRQR(NORB,MAXLLR1,PSIT,PIMAT)
#else
      CALL SRQR(NORB,LLR1,PSIT,PIMAT)
#endif
!====================================================================================
!      call cpu_time(final_time)
!      if(parent)
!     & print *, "QR Factorization: Timing =",final_time-start_time
!
      DEALLOCATE(PSIT,stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
         __LINE__,__FILE__)

!     Constructing the selcted column density matrix using the piovotting information

      ALLOCATE(PCMAT(LLR1,NORB),stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Allocation problem',&
         __LINE__,__FILE__)


!      call cpu_time(start_time)
!#ifdef PARALLEL
#ifdef __PARALLEL
      ALLOCATE(PIMATg(NORB),stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Allocation problem',&
         __LINE__,__FILE__)

      ALLOCATE(PSIC(NORB,NORB),stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Allocation problem',&
         __LINE__,__FILE__)

      call zeroing(PSIC)

      IF(PARAI%ME.EQ.0)THEN
        DO IORB=1,NORB
          PIMATg(IORB)=PIMAT(IORB)
        END DO
      END IF
      !
      CALL mp_bcast(PIMATg,SIZE(PIMATg),parai%source,parai%allgrp)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      call cpu_time(final_time_1)
!      if(parent)
!     & print *, "Pc Mat Construction-1:  =",final_time_1-start_time
!c
      DO IORB=1,NORB
        IF((PIMATg(IORB).GT.PARAI%ME*MAXLLR1).AND. &
          (PIMATg(IORB).LE.(PARAI%ME+1)*MAXLLR1))THEN
           IR=PIMATg(IORB)-PARAI%ME*MAXLLR1
           DO JORB=1,NORB
             PSIC(IORB,JORB)=PSI(IR,JORB)
           END DO
        END IF
      END DO

      CALL mp_sum(PSIC,NORB*NORB,parai%allgrp)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!------------------------------------------------------------------------
!c      call cpu_time(final_time)
!c      if(parent)
!c     & print *, "Pc Mat Construction-2:  =",final_time-final_time_1
!c      DO IR=1,LLR1
!c        DO IORB=1,NORB
!c          DUM=0.D0
!c          DO JORB=1,NORB
!c            DUM=DUM+PSI(IR,JORB)*PSIC(IORB,JORB)
!c          END DO
!c          PCMAT(IR,IORB)=DUM
!c        END DO
!c      END DO
!c      call cpu_time(final_time_1)
!c      if(parent)
!c     & print *, "Pc Mat Construction-3:  =",final_time_1-final_time
!      call cpu_time(start_time)
      CALL DGEMM('N','T',LLR1,NORB,NORB,1.D0,PSI,LLR1,PSIC,NORB, &
                 0.D0,PCMAT,LLR1)
!      call cpu_time(final_time)
!      if(parent)
!     & print *, "Pc Mat 3, Alternative:  =",final_time-start_time

!      CALL FREEM(IP_PSIC)    !SAGAR HACK
#else
!c     FIXME we need to modify using DGEMM like for the parallel version
      DO IR=1,LLR1
        DO IORB=1,NORB
          DUM=0.D0
          DO JORB=1,NORB
            DUM=DUM+PSI(IR,JORB)*PSI(PIMAT(IORB),JORB)
          END DO
          PCMAT(IR,IORB)=DUM
        END DO
      END DO
#endif
!      call cpu_time(final_time)
!      if(parent)
!     &print *, "Pc Mat Construction Total=",final_time-start_time
!c
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ALLOCATE(SMAT(NORB,NORB),stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Allocation problem',&
         __LINE__,__FILE__)
      call zeroing(SMAT)

!c     Constructing the Smatrix
!      call cpu_time(start_time)
!#ifdef PARALLEL
#ifdef __PARALLEL
!c FIXME Only upper diagonal or lower diagonal is required
      DO IORB=1,NORB
        IF((PIMATg(IORB).GT.PARAI%ME*MAXLLR1).AND. &
         (PIMATg(IORB).LE.(PARAI%ME+1)*MAXLLR1))THEN
          IR=PIMATg(IORB)-PARAI%ME*MAXLLR1
          DO JORB=1,NORB
            SMAT(IORB,JORB)=PCMAT(IR,JORB)
          END DO
        END IF
      END DO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DEALLOCATE(PIMATg,stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
         __LINE__,__FILE__)

      CALL mp_sum(SMAT,NORB*NORB,parai%allgrp)
#else
!c FIXME Only upper diagonal or lower diagonal of SMAT is required
      DO IORB=1,NORB
        DO JORB=1,NORB
          SMAT(JORB,IORB)=PCMAT(PIMAT(JORB),IORB)
        END DO
      END DO
#endif
!      call cpu_time(final_time)
!      if(parent)
!     & print *, "S Mat Construction: Timing =",final_time-start_time
!C     Cholesky factorization of SMAT
!
!      call cpu_time(start_time)
!
      CALL DPOTRF('U',NORB,SMAT,NORB,IERR)
!
!      call cpu_time(final_time)
!      if(parent)
!     &print *, "Cholesky : Timing =",final_time-start_time
!
      IF(IERR.NE.0) &
        CALL STOPGM('LOCALIZE_SCDM','ERROR IN DPOTRF', &
         __LINE__,__FILE__)
!------------------------------------------------------------------------
!     HACK SAGAR   CALCULATION OF OVERLAP MATRIX U
      DO IORB=1,NORB
        DO JORB=IORB,NORB
          C(IORB,JORB)=SMAT(IORB,JORB)
        END DO
      END DO
      !
      CALL DTRTRI('U','N',NORB, C,NORB, IERR)
      !
      CALL DGEMM('T','N',NORB,NORB,NORB,1.D0, &
                 PSIC,NORB,C,NORB,0.D0,U,NORB)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DEALLOCATE(PSIC,stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
         __LINE__,__FILE__)
                            !COMMENTED OUT THE PREVIOUS ONE (MUST)
!========================================================================
! FIXME The following is not required and make DTRSM reads lower diagonal
       DO IORB=1,NORB-1
         DO JORB=I+1,NORB
           SMAT(JORB,IORB)=SMAT(IORB,JORB)
         END DO
       END DO
!
!     Building Orthogonal Basis
!      call cpu_time(start_time)
      CALL DTRSM('R','U','N', 'N',LLR1,NORB,1.D0,SMAT,NORB,PCMAT, &
                LLR1)
!      call cpu_time(final_time)
!      if(parent)
!     &print *, "DTRSM : Timing =",final_time-start_time
!     Copying the localized orbitals to PSI
      CALL DCOPY(LLR1*NORB,PCMAT,1,PSI,1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DEALLOCATE(SMAT,stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
         __LINE__,__FILE__)

      DEALLOCATE(PCMAT,stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
         __LINE__,__FILE__)

      DEALLOCATE(PIMAT,stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
         __LINE__,__FILE__)

      DEALLOCATE(C,stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
         __LINE__,__FILE__)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      CALL TIHALT('LOCALIZE_SCDM',ISUB)
      END SUBROUTINE LOCALIZE_SCDM
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
      SUBROUTINE PRQR(M,N,A,PIV)
!     Using Scalapack to do QR factorization with pivotting
!      IMPLICIT NONE
!      INCLUDE 'parac.inc'
      INTEGER             :: M,N,PIV(:)
      REAL(real_8)        :: A(M,N)
!     local variables
      INTEGER             :: ICONTXT, PROW, PCOL, MYROW, MYCOL, &
       MYACOL,MYAROW,ACOLg,AROWg,MB,NB,myunit, IDES_A(9), IA, JA, &
       LSCR, INFO, II,ISUB
      REAL(real_8),ALLOCATABLE        :: QMAT(:),MYSCR(:)

      INTEGER                         :: NUMROC,I, IERR
      INTEGER, ALLOCATABLE            :: IMAP(:,:)
#if defined(_HAS_SCALAPACK)
!=======================================================================
      CALL TISET('PRQR',ISUB)
!      
!     Getting Processor information from BLACS routines
      !CALL BLACS_PINFO(ME,PROCS)
      AROWg=M
      ACOLg=PARAI%NPROC*N
!     Getting Scalapack contxt 
      CALL BLACS_GET(0,0,ICONTXT)
!     Setting processor grid
      PROW=1
      PCOL=PARAI%NPROC
!     Number of rows and columns per block
      MB=M  
      NB=N
!c     Blacs grids initialization
      IF(PARAI%CP_NOGRP.GT.1)THEN

        ALLOCATE(IMAP(1,PCOL),stat=ierr)
        IF (ierr.NE.0)CALL stopgm('LOCALIZE_SCDM','Allocation problem',&
         __LINE__,__FILE__)


        DO I=1,PCOL
          IMAP(1,I)=(I-1)*PARAI%CP_NOGRP+PARAI%CP_INTER_ME
        END DO
        CALL BLACS_GRIDMAP(ICONTXT,IMAP, &
             PROW, PROW, PCOL)
      ELSE
        CALL BLACS_GRIDINIT(ICONTXT,'R',PROW,PCOL) 
      END IF

      CALL BLACS_GRIDINFO(ICONTXT,PROW,PCOL,MYROW,MYCOL)
!c 
      MYAROW=NUMROC(AROWg,MB,MYROW,0,PROW)
      MYACOL=NUMROC(ACOLg,NB,MYCOL,0,PCOL)
!c
      IF((MYAROW.NE.M).OR.(MYACOL.NE.N))THEN
       WRITE(*,*)"MYAROW=",MYAROW, "MYACOL=", MYACOL
       WRITE(*,*)"M=",M, "N=", N
       CALL STOPGM('LOCALIZE_SCDM','ERROR SETTING UP SCALAPCK', &
         __LINE__,__FILE__)
      END IF
      
!c      myunit = 8777+me
!c      write(myunit,*)"Output for processor ",me," to unit ",myunit
!c      write(myunit,*)"Proc ",me,": myrow, mycol in p-array is ", 
!c     &               MYROW,MYCOL
!c      write(myunit,*)"Size of global array is ",AROWg," x ",ACOLg
!c      write(myunit,*)"Size of block is        ",MB," x ",NB
!c      write(myunit,*)"Size of local array is    ",M," x ",N
!c      write(myunit,*)"Size of local array expected",MYAROW," x ",MYACOL

!C     Setting up the matrix description in ides_a array
      CALL DESCINIT(ides_a,AROWg,ACOLg,MB,NB,0,0,icontxt,M,INFO)
!c
      IA=1
      JA=1

      ALLOCATE(QMAT(JA+MIN(M,N)-1),stat=ierr)
      IF (ierr.NE.0)CALL stopgm('LOCALIZE_SCDM','Allocation problem',&
         __LINE__,__FILE__)

!c
      LSCR=-1

      ALLOCATE(MYSCR(2*MYCOL+MYAROW),stat=ierr)
      IF (ierr.NE.0)CALL stopgm('LOCALIZE_SCDM','Allocation problem',&
         __LINE__,__FILE__)

      CALL PDGEQPF(AROWg,ACOLg,A,IA,JA,IDES_A,PIV,QMAT,MYSCR,LSCR,INFO)
      IF(INFO.NE.0)THEN
        PRINT *, 'INFO =', INFO
        CALL STOPGM('LOCALIZE_SCDM','ERROR PSGEQPF', &
         __LINE__,__FILE__)
      END IF
      LSCR=NINT(MYSCR(1))

      DEALLOCATE(MYSCR,stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
         __LINE__,__FILE__)

      ALLOCATE(MYSCR(LSCR),stat=ierr)
      IF (ierr.NE.0)CALL stopgm('LOCALIZE_SCDM','Allocation problem',&
         __LINE__,__FILE__)

      CALL PDGEQPF(AROWg,ACOLg,A,IA,JA,IDES_A,PIV,QMAT,MYSCR,LSCR,INFO)
!c      do ii=1,N
!c         write(myunit,*)ii, PIV(ii)
!c      end do
      IF(INFO.NE.0)THEN
        PRINT *, 'INFO =', INFO
        CALL STOPGM('LOCALIZE_SCDM','ERROR PSGEQPF',  &
         __LINE__,__FILE__)
      END IF


      DEALLOCATE(QMAT,stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
         __LINE__,__FILE__)

      DEALLOCATE(MYSCR,stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
         __LINE__,__FILE__)

!c
      IF(PARAI%CP_NOGRP.GT.1)THEN
        DEALLOCATE(IMAP,stat=ierr)
        IF (ierr.NE.0)CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
          __LINE__,__FILE__)

      END IF
!C
      CALL BLACS_GRIDEXIT(ICONTXT)  !BUG FIXED SAGAR | NNN
!c
      CALL TIHALT('PRQR',ISUB)
#else
      CALL STOPGM('LOCALIZE_SCDM','requires -D__HAS_SCALAPACK', &
           __LINE__,__FILE__)
#endif
!c----------------------------------------------------------------------!
      END SUBROUTINE PRQR 
!c----------------------------------------------------------------------!
!c----------------------------------------------------------------------!
      SUBROUTINE SRQR(M,N,A,PIV)
!C     Using Lapack to do QR factorization with pivotting
!      IMPLICIT NONE
      INTEGER             :: M,N,PIV(:)
      REAL(real_8)        :: A(M,N)
!C    local variables
      REAL(real_8),ALLOCATABLE :: QMAT(:),MYSCR(:)
!c   
      INTEGER             :: IERR, LSCR
!C=======================================================================     
      ALLOCATE(QMAT(MIN(M,N)),stat=ierr)
      IF (ierr.NE.0)CALL stopgm('LOCALIZE_SCDM','Allocation problem',&
         __LINE__,__FILE__)
 
!C
      LSCR=-1

      ALLOCATE(MYSCR(2*M+(M+1)*M),stat=ierr)
      IF (ierr.NE.0)CALL stopgm('LOCALIZE_SCDM','Allocation problem',&
         __LINE__,__FILE__)

      CALL DGEQP3(M,N,A,M,PIV,QMAT,MYSCR,LSCR,IERR)
      LSCR=NINT(MYSCR(1))
      IF(IERR.NE.0) CALL STOPGM('LOCALIZE_SCDM','ERROR in DGEQP3', &
         __LINE__,__FILE__)
!C

      DEALLOCATE(MYSCR,stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
         __LINE__,__FILE__)


      ALLOCATE(MYSCR(LSCR),stat=ierr)
      IF (ierr.NE.0)CALL stopgm('LOCALIZE_SCDM','Allocation problem',&
         __LINE__,__FILE__)

!C
      CALL DGEQP3(M,N,A,M,PIV,QMAT,MYSCR,LSCR,IERR)
      LSCR=NINT(MYSCR(1))


      DEALLOCATE(QMAT,stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
         __LINE__,__FILE__)


      DEALLOCATE(MYSCR,stat=ierr)
      IF (ierr.NE.0) CALL stopgm('LOCALIZE_SCDM','Deallocation problem',&
         __LINE__,__FILE__)

!C      
      END SUBROUTINE SRQR
!c     ==================================================================
!==============================================================================
!==============================================================================
!C     ==================================================================
      SUBROUTINE SELECT_PAIR(RHOIJ,NEGLECT)
!      IMPLICIT NONE
!      INCLUDE 'system.h'
!      INCLUDE 'fft.inc'
!      INCLUDE 'parac.inc'
!      INCLUDE 'scdm.inc'
!     ARGUMENTS
      REAL(real_8)     :: RHOIJ(LLR1)
      LOGICAL          :: NEGLECT
!------------------------------------------------------------------------------
!     LOCAL VARIABLES
      INTEGER          :: K
      REAL(real_8)     :: DUM
!------------------------------------------------------------------------------
      NEGLECT=.FALSE.
      !
      DUM=0.0D0
      !
      DO K=1,LLR1
        DUM=DUM+DABS(RHOIJ(K))
      END DO
      !
      !CALL MY_SUM_D(DUM,1,ALLGRP)
      CALL mp_sum(dum,parai%allgrp)
      !
      DUM=DUM/(DFLOAT(spar%NR1S*spar%NR2S*spar%NR3S))
      !
      IF(DUM.LT.SCDM_CUTOFF)NEGLECT=.TRUE.
      !
      RETURN
      END SUBROUTINE SELECT_PAIR
!==============================================================================
!==============================================================================
      SUBROUTINE G_TO_R_WF(C0,NSTATE,PSI)
!      IMPLICIT NONE
!      INCLUDE 'func.inc'
!      INCLUDE 'system.h'
!      INCLUDE 'fft.inc'
!      INCLUDE 'parac.inc'
!      include 'geq0.inc'
!      include 'cnst.inc'
!     ARGUMENTS
      COMPLEX(real_8)        :: C0(:,:)
      INTEGER                :: NSTATE
      REAL(real_8)           :: PSI(:,:)
!------------------------------------------------------------------------------
!     LOCAL VARIABLES
      INTEGER                :: IA,IG,IWF,IR, ierr
      COMPLEX(real_8),ALLOCATABLE  :: PSIA(:)
      real(real_8)           :: t_1,t_2,dum
!------------------------------------------------------------------------------
      !
      ALLOCATE(psia(maxfftn),STAT=ierr)
      IF(ierr/=0) CALL stopgm('G_TO_R_WF','allocation problem',&
         __LINE__,__FILE__)

      t_1=m_cputime()

      DO IA=1,NSTATE ,2

        CALL zeroing(psia)

        IF(IA.EQ.NSTATE) THEN
          DO IG=1,JGW
             PSIA(NZFS(IG))=C0(IG,IA)
             PSIA(INZS(IG))=DCONJG(C0(IG,IA))
          ENDDO
          IF(GEQ0) PSIA(NZFS(1))=C0(1,IA)
        ELSE
          DO IG=1,JGW
              PSIA(NZFS(IG))=C0(IG,IA)+UIMAG*C0(IG,IA+1)
              PSIA(INZS(IG))=DCONJG(C0(IG,IA))+ &
                             UIMAG*DCONJG(C0(IG,IA+1))
          ENDDO
          IF(GEQ0) PSIA(NZFS(1))=C0(1,IA)+UIMAG*C0(1,IA+1)
        ENDIF
        CALL INVFFTN(PSIA,.TRUE.,parai%allgrp)

!        call phase(psia)       !SAGAR HACK
!--------------------------------------------
        IF(IA.EQ.NSTATE)THEN
          DO IR=1,LLR1
            PSI(IR,IA)=DREAL(PSIA(IR))
          END DO
        ELSE
          DO IR=1,LLR1
            PSI(IR,IA)=DREAL(PSIA(IR))
            PSI(IR,IA+1)=DIMAG(PSIA(IR))
          ENDDO
        ENDIF
!--------------------------------------------
      ENDDO
      !
      t_2=m_cputime()
!     write(6,*)"time 2 of G to R = ",t_2-t_1     

      DEALLOCATE(psia,STAT=ierr)
      IF(ierr/=0) CALL stopgm('G_TO_R_WF','deallocation problem',&
         __LINE__,__FILE__)

      RETURN
      END SUBROUTINE G_TO_R_WF
!==============================================================================
!==========================================================!
END MODULE hfx_utils
