MODULE header_utils
  USE envj,                            ONLY: curdir,&
                                             hname,&
                                             my_pid,&
                                             real_8,&
                                             tjlimit,&
                                             tmpdir,&
                                             user
  USE parac,                           ONLY: paral
  USE readsr_utils,                    ONLY: xstring

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: header

CONTAINS

  ! ==================================================================
  SUBROUTINE header(filename)
    ! ==--------------------------------------------------------------==
    ! ==  Writes a header to the standard output                      ==
    ! ==--------------------------------------------------------------==
    CHARACTER(len=*)                         :: filename

    CHARACTER(len=9)                         :: fformat
    INTEGER                                  :: i1, i2

! ==--------------------------------------------------------------==

    IF (paral%io_parent) THEN
!  _____                            ! 
! / ___ \                           ! 
!| |   | |____   ____ ____          ! 
!| |   | |  _ \ / _  )  _ \         ! 
!| |___| | | | ( (/ /| | | |        ! 
! \_____/| ||_/ \____)_| |_|        ! 
!        |_|                        ! 
!  ______ ______  ______  _____     ! 
! / _____|_____ \|  ___ \(____ \    ! 
!| /      _____) ) | _ | |_   \ \   !
!| |     |  ____/| || || | |   | |  !
!| \_____| |     | || || | |__/ /   ! 
! \______)_|     |_||_||_|_____/    !
!                                   !                             
!  ,ad8888ba,                                                  
! d8"'    `"8b                                                 
!d8'        `8b                                                
!88          88  8b,dPPYba,    ,adPPYba,  8b,dPPYba,           
!88          88  88P'    "8a  a8P_____88  88P'   `"8a          
!Y8,        ,8P  88       d8  8PP"""""""  88       88          
! Y8a.    .a8P   88b,   ,a8"  "8b,   ,aa  88       88          
!  `"Y8888Y"'    88`YbbdP"'    `"Ybbd8"'  88       88          
!                88                                            
!                88                                            
!                                                              
!  ,ad8888ba,   88888888ba   88b           d88  88888888ba,    
! d8"'    `"8b  88      "8b  888b         d888  88      `"8b   
!d8'            88      ,8P  88`8b       d8'88  88        `8b  
!88             88aaaaaa8P'  88 `8b     d8' 88  88         88  
!88             88""""""'    88  `8b   d8'  88  88         88  
!Y8,            88           88   `8b d8'   88  88         8P  
! Y8a.    .a8P  88           88    `888'    88  88      .a8P   
!  `"Y8888Y"'   88           88     `8'     88  88888888Y"'    
!                                                              
!                                                              

       WRITE(6,'(/)')
       WRITE(6,'(12X,A)') '_______                             '
       WRITE(6,'(12X,A)') "__  __ \___________________         "
       WRITE(6,'(12X,A)') "_  / / /__  __ \  _ \_  __ \        "
       WRITE(6,'(12X,A)') '/ /_/ /__  /_/ /  __/  / / /        '
       WRITE(6,'(12X,A)') '\____/ _  .___/\___//_/ /_/         '
       WRITE(6,'(12X,A)') '       /_/                          '
       WRITE(6,'(12X,A)') '______________________  __________  '
       WRITE(6,'(12X,A)') '__  ____/__  __ \__   |/  /__  __ \ '
       WRITE(6,'(12X,A)') '_  /    __  /_/ /_  /|_/ /__  / / / '
       WRITE(6,'(12X,A)') '/ /___  _  ____/_  /  / / _  /_/ /  '
       WRITE(6,'(12X,A)') '\____/  /_/     /_/  /_/  /_____/   '
       WRITE(6,'(/,23X,A,A,/)') '   VERSION ', __GIT_REV
#if defined(__GROMOS)
       WRITE(6,'(/,16X,A,/)') 'COMPILED WITH GROMOS-AMBER QM/MM SUPPORT'
#endif
       WRITE(6,'(14X,A)') '    Home Page: http://www.cpmd.org'
       WRITE(6,'(14X,A)') '            FORKED from  CPMD https://github.com/CPMD-code'
       WRITE(6,'(14X,A)') '              COPYRIGHT'
       WRITE(6,'(14X,A)') '        IBM RESEARCH DIVISION (1990-2023)'
       WRITE(6,'(14X,A)') '  MPI FESTKOERPERFORSCHUNG STUTTGART (1994-2001)'
       WRITE(6,'(/)')
       CALL timetag
       ! ==--------------------------------------------------------------==
       CALL xstring(filename,i1,i2)
       WRITE(fformat,'(A,I2,A)') '(A,T',MAX(21,65-(i2-i1)),',A)'
       WRITE(6,fformat) ' THE INPUT FILE IS: ',ADJUSTL(TRIM(filename))
       CALL xstring(hname,i1,i2)
       WRITE(fformat,'(A,I2,A)') '(A,T',MAX(20,65-(i2-i1)),',A)'
       WRITE(6,fformat) ' THIS JOB RUNS ON: ',ADJUSTL(TRIM(hname))
       WRITE(6,'(A)')   ' THE CURRENT DIRECTORY IS: '
       CALL xstring(curdir,i1,i2)
       WRITE(fformat,'(A,I2,A)') '(T',MAX(2,65-(i2-i1)),',A)'
       WRITE(6,fformat) ADJUSTL(TRIM(curdir))
       WRITE(6,'(A)')   ' THE TEMPORARY DIRECTORY IS: '
       CALL xstring(tmpdir,i1,i2)
       WRITE(fformat,'(A,I2,A)') '(T',MAX(2,65-(i2-i1)),',A)'
       WRITE(6,fformat) ADJUSTL(TRIM(tmpdir))
       WRITE(6,'(A,T50,I16)')  ' THE PROCESS ID IS: ',my_piD
       CALL xstring(user,i1,i2)
       IF (i2.NE.0) THEN
          WRITE(fformat,'(A,I2,A)') '(A,T',MAX(28,65-(i2-i1)),',A)'
          WRITE(6,fformat) ' THE JOB WAS SUBMITTED BY: ',&
               ADJUSTL(TRIM(user))
       ENDIF
       IF (tjlimit.NE.0._real_8) THEN
          WRITE(6,'(A,T48,F10.0,A8)')&
               ' THE JOB TIME LIMIT IS:',tjlimit,' SECONDS'
       ENDIF
       WRITE(6,*)
    ENDIF
    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE header
  ! ==================================================================

END MODULE header_utils
