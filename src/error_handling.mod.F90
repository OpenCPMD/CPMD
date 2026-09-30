MODULE error_handling
  USE parac,                           ONLY: parai
  USE, INTRINSIC :: iso_c_binding,     ONLY: c_funptr, c_null_funptr, &
                                             c_associated, c_f_procpointer, c_int

  !$ USE omp_lib, ONLY: omp_get_thread_num, omp_get_level

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: stopgm

  ! A program that embeds CPMD may install a C function here. stopgm calls
  ! it with the stop code (999) before it prints the call stack. A nonzero
  ! return makes stopgm close the log and return to its caller. The
  ! wavefunction, forces and module state of that call are undefined, so the
  ! host discards them and sets CPMD up again. The function must return.
  ! Unwinding across the Fortran frames that called stopgm is undefined.
  ! Null by default, so cpmd.x stops as it always has.
  TYPE(c_funptr), BIND(C, NAME='cpmd_stopgm_hook'), PUBLIC :: stopgm_hook = c_null_funptr

  ABSTRACT INTERFACE
    FUNCTION stopgm_hook_iface(code) BIND(C) RESULT(handled)
      IMPORT :: c_int
      INTEGER(c_int), VALUE :: code
      INTEGER(c_int) :: handled
    END FUNCTION stopgm_hook_iface
  END INTERFACE
CONTAINS
  ! ==================================================================
  SUBROUTINE stopgm(a,b,line,file)
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*)                         :: a, b
    INTEGER                                  :: line
    CHARACTER(len=*)                         :: file

    CHARACTER(*), PARAMETER                  :: file_name_base = 'LocalError'
    INTEGER, PARAMETER                       :: file_unit = 666

    CHARACTER(100)                           :: buff, file_name
    EXTERNAL                                 :: tistopgm
    PROCEDURE(stopgm_hook_iface), POINTER    :: hook
    INTEGER                                  :: i_level, i_thread, nc

! ==--------------------------------------------------------------==

    !$omp master
    CALL end_swap
    !$omp end master

    i_thread = 0
    i_level = 0
    !$ i_thread = omp_get_thread_num( )
    !$ i_level = omp_get_level( )
    file_name=file_name_base
    WRITE(buff,'(i0,A,i0,A,i0)') parai%cp_me,'-',i_thread,'-',i_level
    file_name=TRIM(file_name)//'-'//TRIM(ADJUSTL(buff))//'.log'
    OPEN(unit=file_unit,file=file_name,action='write')
    WRITE(file_unit,'(5(A,I0))') ' process id''s: ',&
         parai%cp_me,', ',parai%me,', ',parai%cp_inter_me,', ',i_thread,', ',i_level
    WRITE(file_unit,'(A,A)')&
         ' process stops in file: ',TRIM(ADJUSTL(file))
    WRITE(file_unit,'(A,I0)')&
         '               at line: ',line
    WRITE(file_unit,'(A,A)')&
         '               in procedure: ',TRIM(ADJUSTL(a))
    WRITE(file_unit,'(A,A)') ' error message: ',TRIM(ADJUSTL(b))
    nc=999
    IF (c_associated(stopgm_hook)) THEN
      CALL c_f_procpointer(stopgm_hook, hook)
      IF (hook(INT(nc, c_int)) /= 0) THEN
        CLOSE(file_unit)
        RETURN
      END IF
    END IF
    CALL tistopgm(file_unit)
    CLOSE(file_unit)
    CALL my_stopall(nc)
    ! ==--------------------------------------------------------------==
  END SUBROUTINE stopgm
END MODULE error_handling
