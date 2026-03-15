MODULE k_forces_utils
  USE fft_maxfft,                      ONLY: maxfft
  USE nlps,                            ONLY: imagp
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE pslo,                            ONLY: pslo_com
  USE ropt,                            ONLY: ropt_mod
  USE rscpot_utils,                    ONLY: give_scr_rscpot
  USE spin,                            ONLY: lspin2
  USE symtrz_utils,                    ONLY: give_scr_symvec
  USE system,                          ONLY: cnti

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: give_scr_kforces

CONTAINS 
  ! ==================================================================
  SUBROUTINE give_scr_kforces(lforces,tag,nstate,lproj,tfor)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: lforces
    CHARACTER(len=30)                        :: tag
    INTEGER                                  :: nstate
    LOGICAL                                  :: lproj, tfor

    INTEGER :: il_fsc, il_gam, il_ksener, il_psiab, il_scr, &
      lortho, lrscpot, lsymvec

    IF (tfor) THEN
       CALL give_scr_symvec(lsymvec,tag)
    ELSE
       lsymvec=0
    ENDIF
    CALL give_scr_rscpot(lrscpot,tag,ropt_mod%calste)
    il_gam = imagp*nstate*nstate
    il_fsc=nstate
    il_psiab=2
    IF (lspin2%tlse) il_psiab=2*maxfft
    il_scr=il_gam+il_fsc+il_psiab+100
    il_ksener = 2*nstate*nstate + 3*2*nstate + nstate/2 +&
         MAX(1,2*2*nstate-1)
    CALL give_scr_ortho(lortho,tag,nstate)
    lforces = MAX(lrscpot,lsymvec,il_scr,il_ksener,lortho)
    tag='LRSCPOT,LSYMVEC,..'
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_kforces
  ! ==================================================================
END MODULE k_forces_utils
