MODULE benc
  IMPLICIT NONE

  ! ==================================================================
  ! ==                 BENCHMARK OPTIONS                            ==
  ! ==                                                              ==
  ! == IBENCH(1)  : 0  default, no special actions                  ==
  ! ==              1  skip writing the restart files               ==
  ! == IBENCH(2)  : 0  default, only print timings up to 1/1000     ==
  ! ==              1  print every timer triggered                  ==
  ! == IBENCH(3)  : 0  default, autotuning during calculation       ==
  ! ==              1  autotuning during intialization, separate    ==
  ! ==                 timers for autotuning                        ==
  ! == IBENCH(4)  : 0  default, do not print detailed FFT tuning    ==
  ! ==                 timings                                      ==
  ! ==              1  print detailed FFT tuning timings            ==
  ! == IBENCH(5)  : 0  default, do not print scratchlib usage       ==
  ! ==              1  print scrtachlib usage                       ==
  ! ==                                                              ==
  ! ==================================================================
  INTEGER, PARAMETER ::   nbentr=10

  INTEGER :: ibench(nbentr)
  ! ==================================================================

END MODULE benc
