MODULE if_parallel
#ifdef __PARALLEL
  USE mpi_f08
#endif
  !     ==--------------------------------------------------------------==
  !     == IFMEPOS : RANK in the Interface group                        ==
  !     == IFNPROC : Number of processes in the interface group         ==
  !     == IFGRP   : Interface group communicator                       ==
  !     == IFSOURCE: IFMEPOS of  process responsible for the Y-Z plane  ==
  !     == IFSRCPOS: RANK of IFSOURCE in the QMMMGROUP                  ==
  !     == IFPARENT: IFMEPOS==IFSOURCE                                  ==
  !     ==--------------------------------------------------------------==
  TYPE ifparai_t
     INTEGER :: ifmepos 
     INTEGER :: ifnproc 
#ifdef __PARALLEL
     type(MPI_COMM) :: ifgrp 
#else
     INTEGER :: ifgrp 
#endif
     INTEGER :: ifsource 
     INTEGER :: ifsrcpos
  END TYPE ifparai_t
  TYPE (ifparai_t) :: ifparai

  TYPE lifpar_t
     LOGICAL :: ifparent
  END TYPE lifpar_t
  TYPE (lifpar_t) :: lifpar
END MODULE if_parallel
