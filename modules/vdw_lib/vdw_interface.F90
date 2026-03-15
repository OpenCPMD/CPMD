! MODULE vdw_grimme_lib
! 
!   USE vdw_param, only: grimvdwin, vdw_init_data
!   USE vdw_calculator, only: vdw_forces
! 
!   IMPLICIT NONE
!   
!   PRIVATE
!   
!   PUBLIC :: vdw_grimme_read_input
!   PUBLIC :: vdw_grimme_init_data_mpi
!   PUBLIC :: vdw_grimme_calc_energy_forces_stress
! 
!   CONTAINS 

    SUBROUTINE vdw_grimme_read_input(io_unit, num_atom_types, atom_names, lattice_constant, functional, error)
      USE vdw_param, only: grimvdwin, vdw_init_data
      USE vdw_calculator, only: vdw_forces
      INTEGER,INTENT(IN) :: io_unit, num_atom_types
      ! io_unit : Input uint
      ! num_atom_types : number of atomic species
      INTEGER,INTENT(OUT) :: error
      ! error =/ 0 failure
      REAL(KIND=8) :: lattice_constant
      ! lattice constant : lattice parameter, has to be constant
      CHARACTER(LEN=3), INTENT(IN) :: atom_names(num_atom_types)
      CHARACTER(LEN=8), INTENT(IN) :: functional
      ! atom_names : array of element names
      ! functional : DFT functional name, currently only PBE and revPBE are allowed

      CALL grimvdwin(io_unit, num_atom_types, atom_names, lattice_constant, functional, error)

    END SUBROUTINE vdw_grimme_read_input

    SUBROUTINE vdw_grimme_init_data_mpi(master_rank, mpi_rank, mpi_com,error)
      USE vdw_param, only: grimvdwin, vdw_init_data
      USE vdw_calculator, only: vdw_forces
      USE mpi_f08
      INTEGER, INTENT(IN) :: master_rank, mpi_rank
      type(MPI_COMM), intent(in)  :: mpi_com
      ! master_rank : mpi rank of IO proc in supplied mpi communicator
      ! mpi_rank : callers mpi rank in supplied mpi communicator
      ! mpi_com : mpi communicator
      INTEGER, INTENT(OUT) :: error

      call vdw_init_data(master_rank, mpi_rank, mpi_com, error)
    END SUBROUTINE VDW_GRIMME_INIT_DATA_MPI

    SUBROUTINE vdw_grimme_calc_energy_forces_stress(lattice_constant, realsp_lattice_vectors, recipro_lattice_vectors,&
               cel_volume, num_atoms, index_atom_types, atomic_coordinates, vdw_energy, forces, vdw_stress,&
               mpi_rank, num_mpi_procs, mpi_com, error)
     USE vdw_param, only: grimvdwin, vdw_init_data
     USE vdw_calculator, only: vdw_forces
     USE mpi_f08
      REAL(KIND=8), INTENT(IN) :: lattice_constant, realsp_lattice_vectors(3,3), recipro_lattice_vectors(3,3),&
                                  cel_volume, atomic_coordinates(3,num_atoms)
      ! lattice_constant : lattice constant
      ! realsp_lattice_vectors : 3x3 array containing the real space lattice vectors
      ! recipro_lattice_vectors : 3x3 array containing the reciprocal lattice vectors
      ! cel_vol : real space cell volume  
      ! atomic_coordinates : 3xnum_atoms array containing the atomic coordinates in xyz style
      INTEGER, INTENT(IN) :: num_atoms, index_atom_types(num_atoms),mpi_rank,num_mpi_procs
      type(MPI_COMM) :: mpi_com
      INTEGER, INTENT(OUT):: error
      ! num_atoms : number of atoms
      ! index_atom_types : array of length num_atoms, containing a map of atom number to element type number
      REAL(KIND=8), INTENT(OUT) :: vdw_energy, forces(3,num_atoms), vdw_stress(3,3)
      ! vdw_energy : calculated vdw_energy in hartree
      ! vdw_forces :  calculated atomic forces in x y z direction in hartree/ bohr^2
      ! vdw_stress : calculated stress in hartree / bohr ^2
      ! mpi_rank : callers mpi rank
      ! num_mpi_procs : size of mpi communicator
      ! mpi_com : mpi communicator
      
      call vdw_forces(lattice_constant, realsp_lattice_vectors, recipro_lattice_vectors,&
               cel_volume, num_atoms, index_atom_types, atomic_coordinates, vdw_energy, forces, vdw_stress,&
               mpi_rank, num_mpi_procs, mpi_com, error)
    END SUBROUTINE vdw_grimme_calc_energy_forces_stress

!   END MODULE vdw_grimme_lib
