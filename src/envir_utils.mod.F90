MODULE envir_utils
  USE envj,                            ONLY: curdir,&
                                             hname,&
                                             my_pid,&
                                             my_uid,&
                                             real_8,&
                                             tjlimit,&
                                             tmpdir,&
                                             user
  USE machine,                         ONLY: m_getcwd,&
                                             m_getenv,&
                                             m_getlog,&
                                             m_getpid,&
                                             m_getuid,&
                                             m_hostname
  USE system,                          ONLY: cnts

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: envir

CONTAINS

  ! ==================================================================
  SUBROUTINE envir
    ! ==--------------------------------------------------------------==
    ! ==  LOOK FOR THE ENVIRONMENT THIS JOB IS RUNNING                ==
    ! ==--------------------------------------------------------------==
    ! ==--------------------------------------------------------------==
    ! ==  LOOK FOR THE ENVIRONMENT THIS JOB IS RUNNING                ==
    ! ==--------------------------------------------------------------==
    CALL m_getpid(my_pid)
    CALL m_getlog(user)
    CALL m_getuid(my_uid)
    CALL m_hostname(hname)
    CALL m_getcwd(curdir)
    tjlimit=0._real_8
    ! cnts%tmpdir is the host's scratch directory. Blank falls through
    ! to TMPDIR, then the working directory.
    IF (LEN_TRIM(cnts%tmpdir).GT.0) THEN
       tmpdir=TRIM(cnts%tmpdir)
    ELSE
       CALL m_getenv('TMPDIR',tmpdir)
    ENDIF
    IF (tmpdir(1:1).EQ.' ') tmpdir=curdir
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE envir
  ! ==================================================================

END MODULE envir_utils
