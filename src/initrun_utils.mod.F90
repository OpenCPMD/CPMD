MODULE initrun_utils
  USE ainitwf_utils,                   ONLY: give_scr_ainitwf
  USE calc_alm_utils,                  ONLY: give_scr_calc_alm
  USE copot_utils,                     ONLY: give_scr_copot
  USE elct,                            ONLY: crge
  USE fint,                            ONLY: fint1
  USE nlcc,                            ONLY: corel
  USE ortho_utils,                     ONLY: give_scr_ortho
  USE pslo,                            ONLY: pslo_com
  USE rinitwf_utils,                   ONLY: give_scr_rinitwf
  USE store_types,                     ONLY: restart1
  USE system,                          ONLY: cntl
  USE timer,                           ONLY: tihalt,&
                                             tiset

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: give_scr_initrun

CONTAINS
  ! ==================================================================
  SUBROUTINE give_scr_initrun(linitrun,tag)
    ! ==--------------------------------------------------------------==
    INTEGER                                  :: linitrun
    CHARACTER(len=30)                        :: tag

    CHARACTER(*), PARAMETER :: procedureN = 'give_scr_initrun'

    INTEGER                                  :: isub, lainitwf, lcalc_alm, &
                                                lcopot, lortho, &
                                                lrinitwf, &
                                                nstate

    CALL tiset(procedureN,isub)
    nstate=crge%n
    lcopot=0
    lortho=0
    lcalc_alm=0
    lainitwf=0
    linitrun=0
    lrinitwf=0
    CALL give_scr_rinitwf(lrinitwf,tag,nstate)
    IF (restart1%restart) THEN
       IF (corel%tinlc) CALL give_scr_copot(lcopot,tag)
       CALL give_scr_ortho(lortho,tag,nstate)
       IF (.NOT. restart1%rnon) CALL give_scr_ainitwf(lainitwf,tag,nstate)
       lortho=MAX(lortho,lainitwf)
    ENDIF
    IF (cntl%tdiag) THEN
       IF (fint1%ttrot) CALL give_scr_calc_alm(lcalc_alm,tag)
    ENDIF
    linitrun=MAX(lrinitwf,lcopot,lortho,&
         lcalc_alm)
    CALL tihalt(procedureN,isub)
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE give_scr_initrun
  ! ==================================================================


END MODULE initrun_utils
