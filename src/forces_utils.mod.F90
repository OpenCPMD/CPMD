MODULE forces_utils
  USE fft_maxfft,                      ONLY: maxfft
  USE jrotation_utils,                 ONLY: set_orbdist
  USE nlps,                            ONLY: imagp
  USE opeigr_utils,                    ONLY: give_scr_opeigr
  USE parac,                           ONLY: parai
  USE pslo,                            ONLY: pslo_com
  USE ropt,                            ONLY: ropt_mod
  USE rscpot_utils,                    ONLY: give_scr_rscpot
  USE spin,                            ONLY: lspin2
  USE symtrz_utils,                    ONLY: give_scr_symvec
  USE system,                          ONLY: cnti,&
                                             cntl

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: give_scr_forces

CONTAINS
  ! ==================================================================
  SUBROUTINE give_scr_forces(lforces,tag,nstate,lproj,tfor)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lforces
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: lproj, tfor

    INTEGER :: il_amat, il_auxc, il_fsc, il_gam, il_psiab, il_scr, &
      il_scrdip, lrscpot, lsymvec, nstx

    lsymvec=0
    lrscpot=0
    IF (tfor) THEN
       CALL give_scr_symvec(lsymvec,tag)
    ENDIF
    CALL give_scr_rscpot(lrscpot,tag,ropt_mod%calste)
    il_auxc=0
    IF (cntl%tdmal) THEN
       il_gam = 1
    ELSE
       il_gam = imagp*nstate*nstate
    ENDIF
    il_scrdip=0
    IF (cntl%tfield) CALL give_scr_opeigr(il_scrdip,tag,nstate)
    IF (cntl%tdmal) THEN
       CALL set_orbdist(nstate,cnti%nstblk,parai%nproc,nstx)
       il_amat=3*nstate*nstx
    ELSE
       il_amat=0
       il_auxc=MAX(il_auxc,nstate**2)! AUXC space for OVLAP (Laio A.)
    ENDIF
    il_fsc=nstate
    il_psiab=2
    IF (lspin2%tlse) il_psiab=2*maxfft
    il_scr=il_gam+il_auxc+il_fsc+il_psiab+il_scrdip+100
    il_scr=il_scr+il_amat
    lforces = MAX(lrscpot,lsymvec,il_scr)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_forces
  ! ==================================================================

END MODULE forces_utils
