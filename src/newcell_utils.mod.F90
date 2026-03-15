MODULE newcell_utils
  USE initclust_utils,                 ONLY: gf_periodic
  USE nlccset_utils,                   ONLY: nlccset
  USE rggen_utils,                     ONLY: gvector
  USE rinforce_utils,                  ONLY: putps,&
                                             putwnl,&
                                             testspline

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: newcell

CONTAINS

  ! ==================================================================
  SUBROUTINE newcell
    ! ==--------------------------------------------------------------==
    ! ==           INITIALIZE ALL ARRAYS THAT DEPEND ON HT            ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    CALL gvector
    CALL testspline
    CALL putps
    CALL putwnl
    CALL nlccset
    CALL gf_periodic
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE newcell

END MODULE newcell_utils
