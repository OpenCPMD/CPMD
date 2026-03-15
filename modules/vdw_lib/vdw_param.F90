MODULE vdw_param
  !
  ! Module for Grimme-type dispersion correction:
  ! S. Grimme, S. Ehrlich, L. Goerigk, J. Comp. Chem. 32, 1456, (2011)
  ! S. Grimme, J. Antony, S. Ehrlich, H. Krieg, J. Chem. Phys. 132, 154104 (2010)
  ! V. Barone et al., J. Comp. Chem. 30, 934 (2009)
  ! S. Grimme, J. Comp. Chem. 27, 1787 (2006)
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! Some general parameters:
  !
  INTEGER, PRIVATE :: vdw_input
  INTEGER, PUBLIC :: vdw_dir(3)
  INTEGER, ALLOCATABLE, PUBLIC :: vdw_pair(:,:)
  DOUBLE PRECISION, PUBLIC :: rcut2_vdw, rcut2_cn
  CHARACTER(LEN=9) :: version = 'noGrimme '
  CHARACTER(LEN=8) :: vdw_funct = ' '
  !
  ! version   : allowed values: D2, D3org-BJ, D3org-Z, D3mod-BJ, D3mod-Z
  ! vdw_input : allow input of new C6 coefficients (=1) or not (=0)
  ! vdw_dir   : use periodic boundary conditions (=1) or not (=0)
  ! vdw_pair  : apply vdW between pair of species (=1) or not (=0)
  ! vdw_funct : choose Grimme parameters for this functional
  !
  ! rcut2_vdw : cut-off radius squared for vdW interactions (alat units)
  ! rcut2_cn  : cut-off radius squared for coordination number (alat units)
  !
  !
  ! Grimme D2 parameters:
  ! =====================
  !
  DOUBLE PRECISION, PARAMETER, PUBLIC :: beta = 20.d0
  DOUBLE PRECISION, PUBLIC :: s6_D2
  !
  ! beta  : damping function parameter
  ! s6_D2 : global scaling factor
  !
  DOUBLE PRECISION, ALLOCATABLE, PUBLIC  :: C6_ij(:,:), R_sum(:,:)
  !
  ! C6_ij(ityp,jtyp) : C6 coefficients of each atom type pair: sqrt( C6i*C6j )
  ! R_sum(ityp,jtyp) : sum of VdW radii
  !
  !
  ! Grimme D3 parameters:
  ! =====================
  !
  INTEGER, PARAMETER, PUBLIC :: mxCN0 = 8
  DOUBLE PRECISION, PARAMETER, PUBLIC :: rcut_cn0 = 40.d0
  DOUBLE PRECISION, PARAMETER, PUBLIC :: cn_k1 = 16.d0
  DOUBLE PRECISION, PARAMETER, PUBLIC :: cn_k2 = 4.d0/3.d0
  DOUBLE PRECISION, PARAMETER, PUBLIC :: cn_k3 = 4.d0
  !
  ! cn_k3 : reasonable choices are between 3 and 5; this gives smooth
  !         curves with maxima around the integer values; k3=3 give
  !         for CN=0 a slightly smaller value than computed for the
  !         free atom; this also yields to larger CN for atoms in
  !         larger molecules but with the same chemical environment
  !         which is physically not right; values >5 might lead to
  !         bumps in the potential
  !
  DOUBLE PRECISION, PARAMETER, PUBLIC :: s6 = 1.d0
  DOUBLE PRECISION, PARAMETER, PUBLIC :: alp6 = 14.d0
  DOUBLE PRECISION, PARAMETER, PUBLIC :: alp8 = 16.d0
  DOUBLE PRECISION, PUBLIC :: rs6, s8, rs8
  !
  INTEGER, ALLOCATABLE, PUBLIC :: nCN0(:)
  INTEGER                      :: numTypes
  CHARACTER(LEN=8)             :: functional
  DOUBLE PRECISION, ALLOCATABLE, PUBLIC  :: R0ab(:,:), r2r4(:), Rcov(:)
  !
  ! Grimme D3org:
  DOUBLE PRECISION, ALLOCATABLE, PUBLIC  :: C6ab(:,:,:,:,:)
  ! Grimme D3mod:
  DOUBLE PRECISION, ALLOCATABLE, PUBLIC  :: C6a(:,:,:)
#if defined(__HAS_OMP_OFFLOAD)
!$OMP declare target(version,C6ab,C6a,nCN0)
#endif
  !
  CONTAINS
    !---------------------------------------------------------------------------
    ! Initialize parameters
    !---------------------------------------------------------------------------
    !
    SUBROUTINE grimvdwin(iunit, ntyp, namtyp, alat, funct, error)
      !==--------------------------------------------------------------==
      !==  THIS ROUTINE READS THE GRIMME CORRECTION SUBSECTION         ==
      !==--------------------------------------------------------------==
      !==  THE INPUT HAS TO BE IN THE FOLLOWING FORMAT                 ==
      !==    GRIMME CORRECTION                                         ==
      !==      Options                                                 ==
      !==    END GRIMME CORRECTION                                     ==
      !==--------------------------------------------------------------==
      !==  LIST OF OPTIONS                                             ==
      !==                                                              ==
      !==    VDW VERSION                                               ==
      !==      version                                                 ==
      !==    VDW FUNCTIONAL                                            ==
      !==      vdw_funct                                               ==
      !==    VDW RCUT                                                  ==
      !==      rcut_inp                                                ==
      !==    VDW PERIODICITY                                           ==
      !==      vdw_dir(1:3)                                            ==
      !==    VDW INTERACTION PAIRS                                     ==
      !==      vdw_pair(ntype:ntype)                                   ==
      !==    VDW D3 COEFFICIENTS                                       ==
      !==      ncoef                                                   ==
      !==      znuc C6 CN                                              ==
      !==--------------------------------------------------------------==
      !
!For PWscf only:
!      USE io_global, ONLY: ionode, ionode_id
!#if defined __PARA
!      USE mp, ONLY: mp_bcast
!#endif
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: iunit, ntyp
      ! iunit : standard input
      ! ntyp  : number of species
      DOUBLE PRECISION, INTENT(IN) :: alat
      ! alat  : the cell parameter
      CHARACTER(LEN=3), INTENT(IN) :: namtyp(ntyp)
      CHARACTER(LEN=8), INTENT(IN) :: funct
      ! namtyp     : element name
      ! functional : allowed values: PBE, revPBE
      integer,intent(out) :: error
      !
      ! and here the local variables
      !
      INTEGER, ALLOCATABLE :: ielem(:)
      ! ielem : element number; set according to atom name
      INTEGER, EXTERNAL :: get_atom_number
      INTEGER :: i, ityp, itypin, jtyp
      ! ityp,jtyp : counter of atom type
      DOUBLE PRECISION :: rcut_inp
      ! rcut_inp : cut-off radius for vdW interactions
      CHARACTER(LEN=255) :: inpline
!!BM: to be removed:
!      LOGICAL, PARAMETER :: ionode = .true.
!!BM:
      !
      ! set module parameters from input
      functional=funct
      numTypes=ntyp
      ALLOCATE ( vdw_pair(numTypes,numTypes) )
      !
!      IF (IO_PARENT) THEN
        !
        ! Defaults:
        !
        version = ''
        rcut_inp = 0.d0
        vdw_dir(1:3) = -1
        vdw_pair(1:numTypes,1:numTypes) = -1
        !
        ! read additional vdW information
        !
 100    CONTINUE
        READ(iunit,ERR=200,END=200,FMT='(A80)') inpline
        !
        IF ( (INDEX(inpline,'END').NE.0) .AND. (INDEX(inpline,'GRIM').NE.0) ) GOTO 300
        !
        IF ( (INDEX(inpline,'VDW').NE.0) .AND. (INDEX(inpline,'VERSION').NE.0) ) THEN 
          !
          ! read version of Grimme correction scheme
          !
          READ(iunit,ERR=200,END=200,FMT='(a)') inpline
          version = adjustl(inpline)
          GOTO 100
          !
        ELSE IF ( (INDEX(inpline,'VDW').NE.0) .AND. (INDEX(inpline,'FUNCT').NE.0) ) THEN
          !
          ! read functional for Grimme parameters
          !
          READ(iunit,ERR=200,END=200,FMT='(a)') inpline
          vdw_funct = adjustl(inpline)
          GOTO 100
          !
        ELSE IF ( (INDEX(inpline,'VDW').NE.0) .AND. (INDEX(inpline,'RCUT').NE.0) ) THEN
          !
          ! read cutoff radius for vdW interactions
          !
          READ(iunit,ERR=200,END=200,FMT=*) rcut_inp
          GOTO 100
          !
        ELSE IF ( (INDEX(inpline,'VDW').NE.0) .AND. (INDEX(inpline,'PER').NE.0) ) THEN
          !
          ! apply periodic boundary conditions?
          !
          READ(iunit,ERR=200,END=200,FMT=*) vdw_dir(1:3)
          GOTO 100
          !
        ELSE IF ( (INDEX(inpline,'VDW').NE.0) .AND. (INDEX(inpline,'INT').NE.0) ) THEN
          !
          ! specify pairs of species with vdW contribution
          !
          DO ityp = 1 , numTypes
            READ(iunit,ERR=200,END=200,FMT=*) itypin, vdw_pair(ityp,1:numTypes)
            IF (itypin .ne. ityp) then
              CALL VDW_ERROR('GRIMVDWIN','ERROR INPUT VDW PAIRS',error)
              return
           end IF
           ENDDO
          GOTO 100
          !
        ELSE
          ! Dummy line
          GOTO 100
          !
        ENDIF
        !
        ! alternative D2/D3 parameters are read in the respective subroutines
        !
 200    CONTINUE
        CALL VDW_ERROR('GRIMVDWIN','ERROR IN READING INPUT FILE',error)
        return
        !
 300    CONTINUE
        !
        ! Check for reasonalbe input values
        !
        IF ( (version /= 'D2') .AND. &
             (version /= 'D3org-BJ ') .AND. (version /= 'D3org-Z ') .AND. &
             (version /= 'D3surf-BJ') .AND. (version /= 'D3surf-Z') .AND. &
             (version /= 'D3mod-BJ ') .AND. (version /= 'D3mod-Z ') ) then
           CALL VDW_ERROR('GRIMVDWIN','ERROR INPUT VDW VERSION',error)
           return
        endif
        !
        ! Override functional for Grimme parameters
        !
        IF (vdw_funct /= ' ') functional = vdw_funct
        !
        IF (rcut_inp <= 0.d0) then 
           CALL VDW_ERROR('GRIMVDWIN','ERROR INPUT VDW RCUT',error) 
           return
        end if
        !
        DO i = 1, 3
          IF ( (vdw_dir(i) /= 0) .and. (vdw_dir(i) /= 1) ) then
            CALL VDW_ERROR('GRIMVDWIN','ERROR INPUT VDW PERIODICITY',error)
            return
         end if
        ENDDO
        !
        DO ityp = 1 , numTypes
        DO jtyp = ityp , numTypes
          IF ( (vdw_pair(ityp,jtyp) /= 0) .and. (vdw_pair(ityp,jtyp) /= 1) ) then
            CALL VDW_ERROR('GRIMVDWIN','ERROR INPUT VDW PAIRS',error)
            return
         end IF
          IF (vdw_pair(ityp,jtyp) /= vdw_pair(jtyp,ityp)) then
            CALL VDW_ERROR('GRIMVDWIN','ERROR INPUT VDW PAIRS',error)
            return
         end if
        ENDDO
        ENDDO
        !
        ! Printout
        !
        WRITE(6,'(/1x,"Grimme vdW version: ",a9,4x,"functional: ",a)') version, functional
        WRITE(6,'(/5x,"Periodic images: ", 3i3)') vdw_dir(1:3)
        WRITE(6,'(/5x,"Pairwise interaction of species:")')
        DO ityp = 1, numTypes
          WRITE(6,'(5x,i4,":",20i3)') ityp, vdw_pair(ityp,1:numTypes)
        ENDDO
        !
        vdw_input = 0
        IF (version(3:6) == 'surf') THEN
          vdw_input = 1
          version(3:9) = 'org' // version(7:9) // ' '
        ENDIF
        IF (version(3:5) == 'mod') vdw_input = 1
        !
!      ENDIF
      !
!#if defined __PARA   ! --> PWSCF only
!      !
!      ! broadcast data to all processes
!      !
!      CALL mp_bcast(version,  ionode_id)
!      CALL mp_bcast(vdw_funct,ionode_id)
!      CALL mp_bcast(vdw_input,ionode_id)
!      CALL mp_bcast(vdw_dir,  ionode_id)
!      CALL mp_bcast(vdw_pair, ionode_id)
!#endif
      !
!#ifdef PARALLEL
!      CALL MY_CHAR_BCAST(version,IO_SOURCE,CP_GRP)
!      CALL MY_CHAR_BCAST(vdw_funct,IO_SOURCE,CP_GRP)
!      CALL MY_BCAST(vdw_input,1*4,IO_SOURCE,CP_GRP)
!      CALL MY_BCAST(vdw_dir,3*4,IO_SOURCE,CP_GRP)
!      CALL MY_BCAST(vdw_pair,numTypes*numTypes*4,IO_SOURCE,CP_GRP)
!#endif
      !
      ! Set functional dependent parameters
      !
      IF (version(1:2) == 'D2') THEN
        SELECT CASE (functional)
        CASE ('PBE')
          s6_D2 = 0.75d0
        CASE ('revPBE')
          s6_D2 = 1.25d0
        CASE ('BLYP')
          s6_D2 = 1.20d0
        CASE ('BP86')
          s6_D2 = 1.06d0
        CASE ('TPSS')
          s6_D2 = 1.0d0
        CASE ('PBE0')
          s6_D2 = 0.60d0
        CASE ('B3LYP')
          s6_D2 = 1.05d0
        CASE default
          WRITE(6,'(/1x,"Unknown functional: ",a8)') functional
          STOP
        END SELECT
      ELSE IF (version(7:8) == 'BJ') THEN
        SELECT CASE (functional)
        CASE ('PBE')
          rs6 = 0.4289d0
          s8  = 0.7875d0
          rs8 = 4.4407d0
        CASE ('revPBE')
          rs6 = 0.5238d0
          s8  = 2.3550d0
          rs8 = 3.5016d0
        CASE ('RPBE')
          rs6 = 0.8318d0
          s8  = 0.1820d0
          rs8 = 4.0094d0
        CASE ('PBEsol')
          rs6 = 0.4466d0
          s8  = 2.9491d0
          rs8 = 6.1742d0
        CASE ('BLYP')
          rs6 = 0.4298d0
          s8  = 2.6996d0
          rs8 = 4.2359d0
        CASE ('BP86')
          rs6 = 0.3946d0
          s8  = 3.2822d0
          rs8 = 4.8516d0
        CASE ('OPBE')
          rs6 = 0.5512d0
          s8  = 3.3816d0
          rs8 = 2.9444d0
        CASE ('OLYP')
          rs6 = 0.5299d0
          s8  = 2.6205d0
          rs8 = 2.8065d0
        CASE ('HCTC')
          rs6 = 0.3563d0
          s8  = 1.0821d0
          rs8 = 4.3359d0
        CASE ('TPSS')
          rs6 = 0.4535d0
          s8  = 1.9435d0
          rs8 = 4.4752d0
        CASE ('HF')
          rs6 = 0.3385d0
          s8  = 0.9171d0
          rs8 = 2.8830d0
        CASE ('PBE0')
          rs6 = 0.4145d0
          s8  = 1.2177d0
          rs8 = 4.8593d0
        CASE ('revPBE0')
          rs6 = 0.4679d0
          s8  = 1.7588d0
          rs8 = 3.7619d0
        CASE ('HSE06')
          rs6 = 0.383d0
          s8  = 2.310d0
          rs8 = 5.685d0
        CASE ('B3LYP')
          rs6 = 0.3981d0
          s8  = 1.9889d0
          rs8 = 4.4122d0
        CASE ('CAMB3LYP')
          rs6 = 0.3708d0
          s8  = 2.0674d0
          rs8 = 5.4743d0
        CASE default
          WRITE(6,'(/1x,"Unknown functional: ",a8)') functional
          STOP
        END SELECT
      ELSE IF (version(7:7) == 'Z') THEN
        rs8 = 1.d0
        SELECT CASE (functional)
        CASE ('PBE')
          rs6 = 1.217d0
          s8  = 0.722d0
        CASE ('revPBE')
          rs6 = 0.923d0
          s8  = 1.010d0
        CASE ('RPBE')
          rs6 = 0.872d0
          s8  = 0.514d0
        CASE ('PBEsol')
          rs6 = 1.345d0
          s8  = 0.612d0
        CASE ('BLYP')
          rs6 = 1.094d0
          s8  = 1.682d0
        CASE ('BP86')
          rs6 = 1.139d0
          s8  = 1.683d0
        CASE ('OPBE')
          rs6 = 0.837d0
          s8  = 2.055d0
        CASE ('OLYP')
          rs6 = 0.806d0
          s8  = 1.764d0
        CASE ('HCTC')
          rs6 = 1.d0
          s8  = 1.d0
        CASE ('TPSS')
          rs6 = 1.221d0
          s8  = 1.206d0
        CASE ('HF')
          rs6 = 1.158d0
          s8  = 1.746d0
        CASE ('PBE0')
          rs6 = 1.287d0
          s8  = 0.928d0
        CASE ('revPBE0')
          rs6 = 0.949d0
          s8  = 0.792d0
        CASE ('HSE06')
          rs6 = 1.129d0
          s8  = 0.109d0
        CASE ('B3LYP')
          rs6 = 1.261d0
          s8  = 1.703d0
        CASE ('CAMB3LYP')
          rs6 = 1.378d0
          s8  = 1.217d0
        CASE default
          WRITE(6,'(/1x,"Unknown functional: ",a8)') functional
          STOP
        END SELECT
      ELSE IF (version(1:2) .ne. 'D2') THEN
        WRITE(6,'(/1x,"Unknown vdW version: ",a9)') version
        STOP
      ENDIF
      !
      ! Printout
      !
!      IF (IO_PARENT) THEN
        IF (version(1:2) == 'D2') THEN
          WRITE(6,'(/1x,"s6_D2=",f8.4)') s6_D2
        ELSE
          WRITE(6,'(/1x,"s6= ",f8.4,5x,"s8= ",f8.4)') s6, s8
          WRITE(6,'(1x,"rs6=",f8.4,5x,"rs8=",f8.4)') rs6, rs8
        ENDIF
!      ENDIF
      !
      ! Allocate all parameter arrays (on all nodes)
      !
      ALLOCATE ( Rcov(numTypes) )
      !
      IF (version(1:2) == 'D2') THEN
        !
        ALLOCATE ( C6_ij(numTypes,numTypes), R_sum(numTypes,numTypes) )
        !
      ELSE IF (version(1:2) == 'D3') THEN
        !
        ALLOCATE ( R0ab(numTypes,numTypes), r2r4(numTypes), nCN0(numTypes) )
        !
        IF (version(3:5) == 'org') THEN
          !
          ALLOCATE ( C6ab(numTypes,numTypes,mxCN0,mxCN0,3) )
          !
        ELSE IF (version(3:5) == 'mod') THEN
          !
          ALLOCATE ( C6a(numTypes,mxCN0,2) )
          !
        ELSE
          WRITE(6,'(/1x,"Unknown vdW version: ",a9)') version
          STOP
        ENDIF
        !
      ELSE
        WRITE(6,'(/1x,"Unknown vdW version: ",a9)') version
        STOP
      ENDIF
      !
      ! Initialize all parameter arrays (on IO_PARENT only)
      !
!      IF (IO_PARENT) THEN
        !
        ALLOCATE ( ielem(numTypes) )
        !
        DO ityp = 1, numTypes
          ielem(ityp) = get_atom_number( namtyp(ityp) )
        END DO
        !
        ! Coordination numbers are calculated in all cases
        !
        CALL getparam_Rcov(numTypes, ielem, Rcov)
        !
        IF (version(1:2) == 'D2') THEN
          !
          CALL getparam_D2(iunit, numTypes, namtyp, ielem, C6_ij, R_sum)
          !
        ELSE IF (version(1:2) == 'D3') THEN
          !
          CALL getparam_R0ab(numTypes, ielem, R0ab)
          CALL getparam_r2r4(numTypes, ielem, r2r4)
#if defined(__HAS_OMP_OFFLOAD)
!$OMP target update to(r0ab,r2r4)
#endif
          !
          IF (version(3:5) == 'org') THEN
            !
            CALL getparam_D3org(iunit, vdw_input, numTypes, namtyp, ielem, mxCN0, nCN0, C6ab)
            !
#if defined(__HAS_OMP_OFFLOAD)
!$OMP target update to (c6ab,ncn0)
#endif
          ELSE IF (version(3:5) == 'mod') THEN
            !
            CALL getparam_D3mod(iunit, numTypes, namtyp, ielem, mxCN0, nCN0, C6a)
            !
#if defined(__HAS_OMP_OFFLOAD)
!$OMP target update to (c6a,ncn0)
#endif
          ELSE
            WRITE(6,'(/1x,"Unknown vdW version: ",a9)') version
            STOP
          ENDIF
          !
        ELSE
          WRITE(6,'(/1x,"Unknown vdW version: ",a9)') version
          STOP
        ENDIF
        !
        DEALLOCATE ( ielem )
        !
        ! Finally: initialize cutoff radii (in alat units)
        !
        rcut2_vdw = (rcut_inp/alat)**2
        rcut2_cn = (rcut_cn0/alat)**2
        !
        WRITE(6,'(/1x,"rcut_vdw=",f8.2,/1x,"rcut_cn =",f8.2,/)') rcut_inp, rcut_cn0
        !
!      ENDIF
      !
!#if defined __PARA   ! --> PWSCF only
!      !
!      ! broadcast data to all processes
!      !
!      CALL mp_bcast(Rcov,  ionode_id)
!      IF (version(1:2) == 'D2') THEN
!        CALL mp_bcast(C6_ij, ionode_id)
!        CALL mp_bcast(R_sum, ionode_id)
!      ELSE
!        CALL mp_bcast(R0ab,  ionode_id)
!        CALL mp_bcast(r2r4,  ionode_id)
!        CALL mp_bcast(nCN0,  ionode_id)
!        IF (version(3:5) == 'org') THEN
!          CALL mp_bcast(C6ab,  ionode_id)
!        ELSE
!          CALL mp_bcast(C6a,   ionode_id)
!        ENDIF
!      ENDIF
!      CALL mp_bcast(rcut2_vdw, ionode_id)
!      CALL mp_bcast(rcut2_cn,  ionode_id)
!#endif
      !
!#ifdef PARALLEL
!      CALL MY_BCAST(rcov,numTypes*8,IO_SOURCE,CP_GRP)
!      IF (version(1:2) == 'D2') THEN
!        CALL MY_BCAST(C6_ij,numTypes*numTypes*8,IO_SOURCE,CP_GRP)
!        CALL MY_BCAST(R_sum,numTypes*numTypes*8,IO_SOURCE,CP_GRP)
!      ELSE
!        CALL MY_BCAST(R0ab,numTypes*numTypes*8,IO_SOURCE,CP_GRP)
!        CALL MY_BCAST(r2r4,numTypes*8,IO_SOURCE,CP_GRP)
!        CALL MY_BCAST(nCN0,numTypes*4,IO_SOURCE,CP_GRP)
!        IF (version(3:5) == 'org') THEN
!          CALL MY_BCAST(C6ab,numTypes*numTypes*mxCN0*mxCN0*3*8,IO_SOURCE,CP_GRP)
!        ELSE
!          CALL MY_BCAST(C6a,numTypes*mxCN0*8*2,IO_SOURCE,CP_GRP)
!        ENDIF
!      ENDIF
!      CALL MY_BCAST(rcut2_vdw,8,IO_SOURCE,CP_GRP)
!      CALL MY_BCAST(rcut2_cn,8,IO_SOURCE,CP_GRP)
!#endif
      !
      RETURN
      !
    END SUBROUTINE grimvdwin
    !
    !---------------------------------------------------------------------------
    ! Deallocate all arrays
    !---------------------------------------------------------------------------
    !
    SUBROUTINE vdw_dealloc()
      !---------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      IF (ALLOCATED( vdw_pair ) ) DEALLOCATE ( vdw_pair )
      IF (ALLOCATED( C6_ij ) ) DEALLOCATE ( C6_ij )
      IF (ALLOCATED( R_sum ) ) DEALLOCATE ( R_sum )
      IF (ALLOCATED( nCN0  ) ) DEALLOCATE ( nCN0  )
      IF (ALLOCATED( R0ab  ) ) DEALLOCATE ( R0ab  )
      IF (ALLOCATED( r2r4  ) ) DEALLOCATE ( r2r4  )
      IF (ALLOCATED( Rcov  ) ) DEALLOCATE ( Rcov  )
      IF (ALLOCATED( C6ab  ) ) DEALLOCATE ( C6ab  )
      IF (ALLOCATED( C6a   ) ) DEALLOCATE ( C6a   )
      !
      RETURN
      !
    END SUBROUTINE vdw_dealloc
    !

    SUBROUTINE vdw_init_data(master_rank, mpi_rank, mpi_com, error)
      USE MPI_f08
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: master_rank, mpi_rank
      type(MPI_COMM),intent(IN) :: mpi_com
      INTEGER,INTENT(OUT) :: error
      
      
      CALL mpi_bcast(numTypes,1,mpi_integer,MASTER_RANK,MPI_COM,error)
      CALL mpi_bcast(functional,8,mpi_character,MASTER_RANK,MPI_COM,error)
#if defined(__HAS_OMP_OFFLOAD)
!$OMP target enter data map(alloc:numtypes,functional)
!$OMP target update to(numtypes,functional)
#endif
      IF (MASTER_RANK .NE. MPI_RANK) THEN
         ALLOCATE ( vdw_pair(numTypes,numTypes) )
      END IF
      
      CALL mpi_bcast(version,9,mpi_character,MASTER_RANK,MPI_COM,error)
#if defined(__HAS_OMP_OFFLOAD)
!$OMP target enter data map(always,alloc:version)
!$OMP target update to (version)
#endif
      CALL mpi_bcast(vdw_funct,8,mpi_character,MASTER_RANK,MPI_COM,error)
      CALL mpi_BCAST(vdw_input,1,mpi_integer,MASTER_RANK,MPI_COM,error)
      CALL mpi_BCAST(vdw_dir,3,mpi_integer,MASTER_RANK,MPI_COM,error)
      CALL mpi_BCAST(vdw_pair,numTypes*numTypes,mpi_integer,MASTER_RANK,MPI_COM,error)
#if defined(__HAS_OMP_OFFLOAD)
!$OMP target enter data map(alloc:vdw_pair,vdw_dir,vdw_input,vdw_funct)
!$OMP target update to(vdw_pair,vdw_dir,vdw_input,vdw_funct)
#endif
      IF (ERROR .NE. 0) THEN
         call vdw_error('vdw_init_data',' mpi error during bcast',error)
         return
      END IF
      IF (MASTER_RANK .NE. MPI_RANK) THEN
         IF (version(1:2) == 'D2') THEN
            SELECT CASE (functional)
            CASE ('PBE')
               s6_D2 = 0.75d0
            CASE ('revPBE')
               s6_D2 = 1.25d0
            CASE ('BLYP')
               s6_D2 = 1.20d0
            CASE ('BP86')
              s6_D2 = 1.06d0
            CASE ('TPSS')
              s6_D2 = 1.0d0
            CASE ('PBE0')
              s6_D2 = 0.60d0
            CASE ('B3LYP')
              s6_D2 = 1.05d0
            CASE default
               WRITE(6,'(/1x,"Unknown functional: ",a8)') functional
               STOP
            END SELECT
         ELSE IF (version(7:8) == 'BJ') THEN
            SELECT CASE (functional)
            CASE ('PBE')
               rs6 = 0.4289d0
               s8  = 0.7875d0
               rs8 = 4.4407d0
            CASE ('revPBE')
               rs6 = 0.5238d0
               s8  = 2.3550d0
               rs8 = 3.5016d0
            CASE ('RPBE')
               rs6 = 0.8318d0
               s8  = 0.1820d0
               rs8 = 4.0094d0
            CASE ('PBEsol')
               rs6 = 0.4466d0
               s8  = 2.9491d0
               rs8 = 6.1742d0
            CASE ('BLYP')
               rs6 = 0.4298d0
               s8  = 2.6996d0
               rs8 = 4.2359d0
            CASE ('BP86')
              rs6 = 0.3946d0
              s8  = 3.2822d0
              rs8 = 4.8516d0
            CASE ('OPBE')
              rs6 = 0.5512d0
              s8  = 3.3816d0
              rs8 = 2.9444d0
            CASE ('OLYP')
              rs6 = 0.5299d0
              s8  = 2.6205d0
              rs8 = 2.8065d0
            CASE ('HCTC')
              rs6 = 0.3563d0
              s8  = 1.0821d0
              rs8 = 4.3359d0
            CASE ('TPSS')
              rs6 = 0.4535d0
              s8  = 1.9435d0
              rs8 = 4.4752d0
            CASE ('HF')
               rs6 = 0.3385d0
               s8  = 0.9171d0
               rs8 = 2.8830d0
            CASE ('PBE0')
              rs6 = 0.4145d0
              s8  = 1.2177d0
              rs8 = 4.8593d0
            CASE ('revPBE0')
              rs6 = 0.4679d0
              s8  = 1.7588d0
              rs8 = 3.7619d0
            CASE ('HSE06')
              rs6 = 0.383d0
              s8  = 2.310d0
              rs8 = 5.685d0
            CASE ('B3LYP')
              rs6 = 0.3981d0
              s8  = 1.9889d0
              rs8 = 4.4122d0
            CASE ('CAMB3LYP')
              rs6 = 0.3708d0
              s8  = 2.0674d0
              rs8 = 5.4743d0
            CASE default
               WRITE(6,'(/1x,"Unknown functional: ",a8)') functional
               STOP
            END SELECT
         ELSE IF (version(7:7) == 'Z') THEN
            rs8 = 1.d0
            SELECT CASE (functional)
            CASE ('PBE')
               rs6 = 1.217d0
               s8  = 0.722d0
            CASE ('revPBE')
               rs6 = 0.923d0
               s8  = 1.010d0
            CASE ('RPBE')
               rs6 = 0.872d0
               s8  = 0.514d0
            CASE ('PBEsol')
               rs6 = 1.345d0
               s8  = 0.612d0
            CASE ('BLYP')
               rs6 = 1.094d0
               s8  = 1.682d0
            CASE ('BP86')
              rs6 = 1.139d0
              s8  = 1.683d0
            CASE ('OPBE')
              rs6 = 0.837d0
              s8  = 2.055d0
            CASE ('OLYP')
              rs6 = 0.806d0
              s8  = 1.764d0
            CASE ('HCTC')
              rs6 = 1.d0
              s8  = 1.d0
            CASE ('TPSS')
              rs6 = 1.221d0
              s8  = 1.206d0
            CASE ('HF')
               rs6 = 1.158d0
               s8  = 1.746d0
            CASE ('PBE0')
              rs6 = 1.287d0
              s8  = 0.928d0
            CASE ('revPBE0')
              rs6 = 0.949d0
              s8  = 0.792d0
            CASE ('HSE06')
              rs6 = 1.129d0
              s8  = 0.109d0
            CASE ('B3LYP')
              rs6 = 1.261d0
              s8  = 1.703d0
            CASE ('CAMB3LYP')
              rs6 = 1.378d0
              s8  = 1.217d0
            CASE default
               WRITE(6,'(/1x,"Unknown functional: ",a8)') functional
               STOP
            END SELECT
         ELSE IF (version(1:2) .ne. 'D2') THEN
            WRITE(6,'(/1x,"Unknown vdW version: ",a9)') version
            STOP
         END IF
         ! Allocate all parameter arrays (on all nodes)
         !
         ALLOCATE ( Rcov(numTypes) )
         !
         IF (version(1:2) == 'D2') THEN
            !
            ALLOCATE ( C6_ij(numTypes,numTypes), R_sum(numTypes,numTypes) )
            !
         ELSE IF (version(1:2) == 'D3') THEN
            !
            ALLOCATE ( R0ab(numTypes,numTypes), r2r4(numTypes), nCN0(numTypes) )
            !
            IF (version(3:5) == 'org') THEN
               !
               ALLOCATE ( C6ab(numTypes,numTypes,mxCN0,mxCN0,3) )
               !
            ELSE IF (version(3:5) == 'mod') THEN
               !
               ALLOCATE ( C6a(numTypes,mxCN0,2) )
               !
            ELSE
               WRITE(6,'(/1x,"Unknown vdW version: ",a9)') version
               ERROR=1
               RETURN
            ENDIF
            !
         ELSE
            WRITE(6,'(/1x,"Unknown vdW version: ",a9)') version
            ERROR=1
            RETURN
         ENDIF
      ENDIF
      CALL MPI_BCAST(rcov,numTypes*8,MPI_BYTE,MASTER_RANK,MPI_COM,ERROR)
#if defined(__HAS_OMP_OFFLOAD)
!$OMP target enter data map (alloc:rcov)
!$OMP target update to (rcov)
#endif
      IF (version(1:2) == 'D2') THEN
         CALL MPI_BCAST(C6_ij,numTypes*numTypes*8,MPI_BYTE,MASTER_RANK,MPI_COM,ERROR)
         CALL MPI_BCAST(R_sum,numTypes*numTypes*8,MPI_BYTE,MASTER_RANK,MPI_COM,ERROR)
      ELSE
         CALL MPI_BCAST(R0ab,numTypes*numTypes*8,MPI_BYTE,MASTER_RANK,MPI_COM,ERROR)
         CALL MPI_BCAST(r2r4,numTypes*8,MPI_BYTE,MASTER_RANK,MPI_COM,ERROR)
         CALL MPI_BCAST(nCN0,numTypes*4,MPI_BYTE,MASTER_RANK,MPI_COM,ERROR)
#if defined(__HAS_OMP_OFFLOAD)
!$OMP target enter data map (alloc:ncn0,r2r4,r0ab)
!$OMP target enter data map (always,alloc:ncn0)
!$OMP target update to(ncn0,r2r4,r0ab)
#endif
         IF (version(3:5) == 'org') THEN
            CALL MPI_BCAST(C6ab,numTypes*numTypes*mxCN0*mxCN0*3*8,MPI_BYTE,MASTER_RANK,MPI_COM,ERROR)
#if defined(__HAS_OMP_OFFLOAD)
!$OMP target enter data map(always,alloc:c6ab)
!$OMP target update to(c6ab)
#endif
         ELSE
            CALL MPI_BCAST(C6a,numTypes*mxCN0*8*2,MPI_BYTE,MASTER_RANK,MPI_COM,ERROR)
#if defined(__HAS_OMP_OFFLOAD)
!$OMP target enter data map(always,alloc:c6a)
!$OMP target update to(c6a)
#endif
         ENDIF
      ENDIF
      CALL MPI_BCAST(rcut2_vdw,8,MPI_BYTE,MASTER_RANK,MPI_COM,ERROR)
#if defined(__HAS_OMP_OFFLOAD)
!$OMP target enter data map (alloc:rcut2_vdw)
!$OMP target update to(rcut2_vdw)
#endif
      CALL MPI_BCAST(rcut2_cn,8,MPI_BYTE,MASTER_RANK,MPI_COM,ERROR)
#if defined(__HAS_OMP_OFFLOAD)
!$OMP target enter data map (alloc:rcut2_cn)
!$OMP target update to(rcut2_cn)
#endif
      IF (ERROR .NE. 0) then 
         call vdw_error('vdw_init_data', 'mpi error during bcast',error)
         return
      END IF
    END SUBROUTINE vdw_init_data


    SUBROUTINE VDW_ERROR(a,b,err)
      implicit none
      character(*),intent(in) :: a,b
      integer,intent(out) :: err
      write(6,*) a,b
      err=1
      return
    end SUBROUTINE VDW_ERROR

END MODULE vdw_param
!
!---------------------------------------------------------------------------
SUBROUTINE getparam_D2(iunit, ntyp, namtyp, ielem, C6_ij, R_sum,error)
  !---------------------------------------------------------------------
  !
  ! Grimme D2 scheme.
  ! Tabulated C6 values are already in Rydberg units.
  !
  use vdw_param,     only: vdw_error
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iunit, ntyp, ielem(ntyp)
  DOUBLE PRECISION, INTENT(OUT) :: C6_ij(ntyp,ntyp), R_sum(ntyp,ntyp)
  CHARACTER(LEN=3), INTENT(IN) :: namtyp(ntyp)
  !
  ! and here the local variables
  !
  INTEGER, PARAMETER :: mxelem = 54 ! 86
  !
  INTEGER :: ierr, iel, jel, il, ityp, jtyp, ncoef
  INTEGER :: read_input
  INTEGER, INTENT(OUT) :: error
  DOUBLE PRECISION D2_C6_param(mxelem), D2_Rcov_param(mxelem), znuc, C6in, Rin
  CHARACTER(LEN=255) :: inpline
  !
  ! Grimme D2 vdw C6 and radii for the first 54 elements
  !
  DATA D2_C6_param(1:mxelem) / &
         4.857,   2.775,  55.853,  55.853, 108.584,  60.710, & !  H-C
        42.670,  24.284,  26.018,  21.855, 198.087, 198.087, & !  N-Mg
       374.319, 320.200, 271.980, 193.230, 175.885, 159.927, & ! Al-Ar
       374.666, 374.666, 374.666, 374.666, 374.666, 374.666, & !  K-Cr
       374.666, 374.666, 374.666, 374.666, 374.666, 374.666, & ! Mn-Zn
       589.405, 593.221, 567.896, 438.498, 432.600, 416.642, & ! Ga-Kr
       855.833, 855.833, 855.833, 855.833, 855.833, 855.833, & ! Rb-Mo
       855.833, 855.833, 855.833, 855.833, 855.833, 855.833, & ! Tc-Cd
      1294.678,1342.899,1333.532,1101.101,1092.775,1040.391 /  ! In-Xe
!  Conversion factor: 1 J/mol*nm^6 = 34.690531 Ry*bohr^6
!    10937.057,7874.542,6114.275,4880.264,4880.264,4880.264, & ! Cs-Nd
!     4880.264,4880.264,4880.264,4880.264,4880.264,4880.264, & ! Pm-Dy
!     4880.264,4880.264,4880.264,4880.264,4880.264,3646.391, & ! Ho-Hf
!     2818.259,2818.259,2818.259,2818.259,2818.259,2818.259, & ! Ta-Pt
!     2818.259,1989.988,1986.172,2191.123,2204.236,1917.797, & ! Au-Po
!     1983.292,1964.872 /                                    & ! At-Rn
  !
  DATA D2_Rcov_param(1:mxelem) / &
         1.892,   1.912,   1.559,   2.661,   2.806,   2.744, & !  H-C
         2.640,   2.536,   2.432,   2.349,   2.162,   2.578, & !  N-Mg
         3.097,   3.243,   3.222,   3.180,   3.097,   3.014, & ! Al-Ar
         2.806,   2.785,   2.952,   2.952,   2.952,   2.952, & !  K-Cr
         2.952,   2.952,   2.952,   2.952,   2.952,   2.952, & ! Mn-Zn
         3.118,   3.264,   3.326,   3.347,   3.305,   3.264, & ! Ga-Kr
         3.076,   3.035,   3.097,   3.097,   3.097,   3.097, & ! Rb-Mo
         3.097,   3.097,   3.097,   3.097,   3.097,   3.097, & ! Tc-Cd
         3.160,   3.409,   3.555,   3.575,   3.575,   3.555 /  ! In-Xe
!        3.405,   3.330,   3.251,   3.313,   3.313,   3.313, & ! Cs-Nd
!        3.313,   3.313,   3.313,   3.313,   3.313,   3.313, & ! Pm-Dy
!        3.313,   3.313,   3.313,   3.313,   3.313,   3.378, & ! Ho-Hf
!        3.349,   3.349,   3.349,   3.349,   3.349,   3.349, & ! Ta-Pt
!        3.349,   3.322,   3.752,   3.673,   3.586,   3.789, & ! Au-Po
!        3.762,   3.636 /                                    & ! At-Rn
    !
    ! Read alternative C6 coefficients
    !
    ierr = read_input(iunit,'&VDW')
    IF (ierr .NE. 0) then
        CALL VDW_ERROR('GETPARAM_D2','SECTION &VDW NOT FOUND',error)
        return
     end if
    ierr = read_input(iunit,'GRIM')
    IF (ierr .NE. 0) then
       CALL VDW_ERROR('GETPARAM_D2','SECTION GRIMME CORRECTION NOT FOUND',error)
       return
    end if
    !
    ncoef = 0
 20 CONTINUE
    READ(iunit,ERR=60,END=60,FMT='(A80)') inpline
    IF ( (INDEX(inpline,'END').NE.0) .AND. (INDEX(inpline,'GRIM').NE.0) ) GOTO 80
    !
    IF ( (INDEX(inpline,'VDW').NE.0) .AND. &
         (INDEX(inpline,'D2').NE.0) .AND. (INDEX(inpline,'COEF').NE.0) ) THEN
      !
      READ(iunit,ERR=60,END=60,FMT=*) ncoef
      IF (ncoef <= 0) THEN 
         CALL VDW_ERROR('GETPARAM_D2','ILLEGAL NCOEF',error)
         RETURN
      END IF

      DO il = 1, ncoef
        READ(iunit,ERR=60,END=60,FMT=*) znuc, C6in, Rin
        iel = nint(znuc)
        IF ( (iel < 1) .or. (iel > mxelem) .or. (abs(znuc-dble(iel)) > 1.d-6) ) then
          CALL VDW_ERROR('GETPARAM_D2','ERROR WHILE READING D2 COEF',error)
          return
       end if
       IF ( (C6in <= 0.d0) .or. (Rin <= 0.d0) ) then
          CALL VDW_ERROR('GETPARAM_D2','ERROR WHILE READING D2 COEF',error)
          return
       end if
        D2_C6_param(iel) = C6in
        D2_Rcov_param(iel) = Rin
      ENDDO
    ENDIF
    GOTO 20
    !
 60 CONTINUE
    CALL VDW_ERROR('GETPARAM_D2','ERROR WHILE READING D2 COEF',error)
    return
    !
 80 CONTINUE
    !
  ! Here we store C6_ij parameters of each pair of atom types into
  ! a square matrix C6_ij = sqrt ( C6_i * C6_j )
  !
  DO jtyp = 1, ntyp
    jel = ielem(jtyp)
    IF ( (jel < 1) .OR. (jel > mxelem) ) THEN
      WRITE(6,'(/1x,"No D2 parameters for element:",i4)') jel
      STOP
    ENDIF
    !
    DO ityp = 1, ntyp
      iel = ielem(ityp)
      IF ( (iel < 1) .OR. (iel > mxelem) ) THEN
        WRITE(6,'(/1x,"No D2 parameters for element:",i4)') iel
        STOP
      ENDIF
      !
      C6_ij(ityp,jtyp) = sqrt( D2_C6_param(iel)*D2_C6_param(jel) )
      R_sum(ityp,jtyp) = D2_Rcov_param(iel) + D2_Rcov_param(jel)
      !
    ENDDO
  ENDDO
  !
  WRITE(6,'(/5x,"-------------------------------------", &
            /5x,"Parameters for Dispersion Correction:", &
            /5x,"-------------------------------------", &
            /5x,"  atom      VdW radius       C_6     ",/)')
  DO ityp = 1, ntyp
    iel = ielem(ityp)
    WRITE(6,'(8x,a3,6x,f7.3,8x,f7.3)') &
      namtyp(ityp), D2_Rcov_param(iel), D2_C6_param(iel)
  ENDDO
  !
  IF (ncoef > 0) THEN
     WRITE(6,'(/5x,"Number of modified D2 coefficients:",i3)') ncoef
     RETURN
  END IF
  !
END SUBROUTINE getparam_D2
!
!---------------------------------------------------------------------------
SUBROUTINE getparam_D3org(iunit, vdw_input, ntyp, namtyp, ielem, mxCN0, nCN0, C6ab,error)
  !---------------------------------------------------------------------
  !
  ! Original Grimme D3 scheme.
  ! Tabulated C6 values are in Hartree units. Convert to Rydberg.
  !
  use vdw_param,     only: vdw_error
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iunit, vdw_input, ntyp, ielem(ntyp), mxCN0
  INTEGER, INTENT(OUT) :: nCN0(ntyp)
  DOUBLE PRECISION, INTENT(OUT) :: C6ab(ntyp,ntyp,mxCN0,mxCN0,3)
  CHARACTER(LEN=3), INTENT(IN) :: namtyp(ntyp)
  INTEGER, INTENT(OUT) :: error
  !
  ! and here the local variables
  !
  INTEGER, PARAMETER :: mxelem = 94
  INTEGER, PARAMETER :: mxcoef = 32385
  !
  INTEGER :: id, il, ierr, iel, jel, ityp, jtyp, icn, jcn, ncoef
  INTEGER :: icn_replaced(ntyp), nCN0org(ntyp)
  INTEGER :: read_input
  DOUBLE PRECISION d3org_param(5*mxcoef), C6in, C6av, CN0in, CN0av, znuc
  CHARACTER(LEN=255) :: inpline
  !
  include 'vdw_param-D3org.inc'
  !
  nCN0(1:ntyp) = 0
  C6ab(1:ntyp,1:ntyp,1:mxCN0,1:mxCN0,1:3) = -1.D0
  !
  il = 1
  DO id = 1, mxcoef
    iel = nint(d3org_param(il+1))
    jel = nint(d3org_param(il+2))
    icn = iel/100
    jcn = jel/100
    iel = iel - icn*100
    jel = jel - jcn*100
    icn = icn + 1
    jcn = jcn + 1
    !
    DO ityp = 1, ntyp
      IF ( iel == ielem(ityp) ) THEN
        !
        DO jtyp = 1, ntyp
          IF ( jel == ielem(jtyp) ) THEN
            !
            nCN0(ityp) = max( nCN0(ityp), icn )
            nCN0(jtyp) = max( nCN0(jtyp), jcn )
            !
            C6ab(ityp,jtyp,icn,jcn,1) = 2.d0*d3org_param(il)
            C6ab(ityp,jtyp,icn,jcn,2) = d3org_param(il+3)
            C6ab(ityp,jtyp,icn,jcn,3) = d3org_param(il+4)
            !
            C6ab(jtyp,ityp,jcn,icn,1) = 2.d0*d3org_param(il)
            C6ab(jtyp,ityp,jcn,icn,2) = d3org_param(il+4)
            C6ab(jtyp,ityp,jcn,icn,3) = d3org_param(il+3)
            !
          ENDIF
        ENDDO
        !
      ENDIF
    ENDDO
    !
    il = il + 5
  ENDDO
  !
    ! Read additional C6 coefficients
    !
    ierr = read_input(iunit,'&VDW')
    IF (ierr .NE. 0) then
       CALL VDW_ERROR('GETPARAM_D3','SECTION &VDW NOT FOUND',error)
       return
    end if
    ierr = read_input(iunit,'GRIM')
    IF (ierr .NE. 0) then 
       CALL VDW_ERROR('GETPARAM_D3','SECTION GRIMME CORRECTION NOT FOUND',error)
       return
    end if
    !
    ncoef = 0
    nCN0org(1:ntyp) = nCN0(1:ntyp)
    icn_replaced(1:ntyp) = 0
 20 CONTINUE
    READ(iunit,ERR=60,END=60,FMT='(A80)') inpline
    IF ( (INDEX(inpline,'END').NE.0) .AND. (INDEX(inpline,'GRIM').NE.0) ) GOTO 80
    !
    IF ( (INDEX(inpline,'VDW').NE.0) .AND. &
         (INDEX(inpline,'D3').NE.0) .AND. (INDEX(inpline,'COEF').NE.0) ) THEN
      !
      IF ( vdw_input == 0 ) THEN
        WRITE(6,'(/1x,"Modification of D3 coefficients is not allowed in D3org")')
        CALL VDW_ERROR('GETPARAM_D3','ERROR WHILE READING D3 COEF',error)
        return
      ENDIF
      !
      READ(iunit,ERR=60,END=60,FMT=*) ncoef
      IF (ncoef <= 0) then 
         CALL VDW_ERROR('GETPARAM_D3','ILLEGAL NCOEF',error)
         return
      end if
      DO il = 1, ncoef
        READ(iunit,ERR=60,END=60,FMT=*) znuc, C6in, CN0in
        iel = nint(znuc)
        IF ( (iel < 1) .or. (iel > mxelem) .or. (abs(znuc-dble(iel)) > 1.d-6) ) then
           CALL VDW_ERROR('GETPARAM_D3','ERROR WHILE READING D3 COEF',error)
           return
        end if
        IF ( (C6in <= 0.d0) .or. (CN0in <= 0.d0) ) then
          CALL VDW_ERROR('GETPARAM_D3','ERROR WHILE READING D3 COEF',error)
          return
       end if
        !
        DO ityp = 1, ntyp
          IF ( ielem(ityp) == iel ) THEN
            !
            ! Check, if new parameters replace old D3org values or are appended at end of list
            ! Allow replacement only once for each element
            !
            DO icn = 1, nCN0org(ityp)
              IF ( (C6ab(ityp,ityp,icn,icn,2) >= CN0in) ) THEN
                IF ( icn_replaced(ityp) == icn ) CYCLE
                IF ( icn_replaced(ityp) == 0 ) THEN
                  WRITE(6,'(/1x,"WARNING: Removing D3org parameter: ",a3,i2,2f12.4)') &
                    namtyp(ityp), icn-1, c6ab(ityp,ityp,icn,icn,1), c6ab(ityp,ityp,icn,icn,2)
                  icn_replaced(ityp) = icn
                  GOTO 40
                ELSE
                  WRITE(6,'(/1x,"Only one D3org value per element can be replaced",/)')
                  CALL VDW_ERROR('GETPARAM_D3','ERROR WHILE READING D3 COEF',error)
                ENDIF
              ENDIF
            ENDDO
            icn = nCN0(ityp) + 1
            IF ( icn > mxCN0 ) THEN
              WRITE(6,'(/1x,"Too many CN0 values for element:",i4,/)') iel
              CALL VDW_ERROR('GETPARAM_D3','ERROR WHILE READING D3 COEF',error)
              return
            ENDIF
            nCN0(ityp) = icn
            !
 40         CONTINUE
            !
            C6ab(ityp,ityp,icn,icn,1) = C6in
            C6ab(ityp,ityp,icn,icn,2) = CN0in
            C6ab(ityp,ityp,icn,icn,3) = CN0in
            !
            DO jtyp = 1, ntyp
              DO jcn = 1, nCN0(jtyp)
                C6av = sqrt( C6ab(jtyp,jtyp,jcn,jcn,1)*C6in )
                CN0av = C6ab(jtyp,jtyp,jcn,jcn,2)
                !
                C6ab(ityp,jtyp,icn,jcn,1) = C6av
                C6ab(ityp,jtyp,icn,jcn,2) = CN0in
                C6ab(ityp,jtyp,icn,jcn,3) = CN0av
                !
                C6ab(jtyp,ityp,jcn,icn,1) = C6av
                C6ab(jtyp,ityp,jcn,icn,2) = CN0av
                C6ab(jtyp,ityp,jcn,icn,3) = CN0in
              ENDDO
            ENDDO
            !
          ENDIF
        ENDDO
      ENDDO
      !
    ENDIF
    GOTO 20
    !
 60 CONTINUE
    CALL VDW_ERROR('GETPARAM_D3','ERROR WHILE READING D3 COEF',error)
    return
    !
 80 CONTINUE
    !
  !  Check, that we found C6 coefficients for all pairs of elements
  !
  DO jtyp = 1, ntyp
  DO ityp = 1, ntyp
    DO jcn = 1, nCN0(jtyp)
    DO icn = 1, nCN0(ityp)
      IF (C6ab(ityp,jtyp,icn,jcn,1) <= 0.D0) THEN
        WRITE(6,'(/1x,"C6 missing for:")')
        WRITE(6,'(5x,"ityp=",i2,", icn=",i2," of",i2,5x,"jtyp=",i2,", jcn=",i2," of",i2)') &
          ityp, icn, nCN0(ityp), jtyp, jcn, nCN0(jtyp)
        STOP
      ENDIF
    ENDDO
    ENDDO
! DEBUG:
!    write(19,'(4i5)') ityp, nCN0(ityp), jtyp, nCN0(jtyp)
!    write(19,'(6f12.4)') C6ab(ityp,jtyp,1:nCN0(ityp),1:nCN0(jtyp),1)
!    write(19,*)
  ENDDO
  ENDDO
  !
  WRITE(6,'(/5x,"-------------------------------------", &
            /5x,"Parameters for Dispersion Correction:", &
            /5x,"-------------------------------------", &
            /5x,"  atom           C_6        CN0      ",/)')
  DO ityp = 1, ntyp
    do icn = 1, nCN0(ityp)
      WRITE(6,'(8x,a3,i3,5x,f8.3,4x,f8.4)') namtyp(ityp), icn-1, &
        c6ab(ityp,ityp,icn,icn,1), c6ab(ityp,ityp,icn,icn,2)
    ENDDO
  ENDDO
  !
  IF ( (vdw_input == 1) .AND. (ncoef > 0) ) THEN
     WRITE(6,'(/5x,"Number of additional D3 coefficients:",i3)') ncoef
     RETURN
  END IF
  !
END SUBROUTINE getparam_D3org
!
!---------------------------------------------------------------------------
SUBROUTINE getparam_D3mod(iunit, ntyp, namtyp, ielem, mxCN0, nCN0, C6a,error)
  !---------------------------------------------------------------------
  !
  ! Simplyfied version of Grimme D3 scheme.
  ! Tabulated C6 values are in Hartree units. Convert to Rydberg.
  !
  use vdw_param,     only: vdw_error
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iunit, ntyp, ielem(ntyp), mxCN0
  INTEGER, INTENT(OUT) :: nCN0(ntyp)
  DOUBLE PRECISION, INTENT(OUT) :: C6a(ntyp,mxCN0,2)
  CHARACTER(LEN=3), INTENT(IN) :: namtyp(ntyp)
  INTEGER, INTENT(OUT) :: error
  !
  ! and here the local variables
  !
  INTEGER, PARAMETER :: mxelem = 94
  INTEGER, PARAMETER :: mxcoef = 254
  !
  INTEGER :: icn, id, ierr, il, iel, ityp, ncoef
  INTEGER :: read_input
  DOUBLE PRECISION d3mod_param(3*mxcoef), C6in, CN0in, znuc
  CHARACTER(LEN=255) :: inpline
  !
  include 'vdw_param-D3mod.inc'
  !
  C6a(1:ntyp,1:mxCN0,1:2) = 0.D0
  !
  DO ityp = 1, ntyp
    iel = ielem(ityp)
    icn = 0
    il = 1
    DO id = 1, mxcoef
      IF (nint(d3mod_param(il)) == iel) THEN
        icn = icn + 1
        IF (icn > mxCN0) THEN
          WRITE(6,'(/1x,"Too many CN0 values for element:",i4,/)') iel
          CALL VDW_ERROR('GETPARAM_D3','ERROR WHILE READING D3 COEF',error)
          return
        ENDIF
        C6a(ityp,icn,1) = 2.d0*d3mod_param(il+1)
        C6a(ityp,icn,2) = d3mod_param(il+2)
      ENDIF
      il = il + 3
    ENDDO
    IF (icn == 0) THEN
      WRITE(6,'(/1x,"No D3mod parameters for element:",i4)') iel
      STOP
    ENDIF
    nCN0(ityp) = icn
  ENDDO
  !
    ! Read additional C6 coefficients
    !
    ierr = read_input(iunit,'&VDW')
    IF (ierr .NE. 0) then
       CALL VDW_ERROR('GETPARAM_D3','SECTION &VDW NOT FOUND',error)
       return
    end IF
    ierr = read_input(iunit,'GRIM')
    IF (ierr .NE. 0) then 
       CALL VDW_ERROR('GETPARAM_D3','SECTION GRIMME CORRECTION NOT FOUND',error)
       return
    end IF
    !
    ncoef = 0
 20 CONTINUE
    READ(iunit,ERR=60,END=60,FMT='(A80)') inpline
    IF ( (INDEX(inpline,'END').NE.0) .AND. (INDEX(inpline,'GRIM').NE.0) ) GOTO 80
    !
    IF ( (INDEX(inpline,'VDW').NE.0) .AND. &
         (INDEX(inpline,'D3').NE.0) .AND. (INDEX(inpline,'COEF').NE.0) ) THEN
      !
      READ(iunit,ERR=60,END=60,FMT=*) ncoef
      IF (ncoef <= 0) then
         CALL VDW_ERROR('GETPARAM_D3','ILLEGAL NCOEF',error)
         return
      end if
      DO il = 1, ncoef
        READ(iunit,ERR=60,END=60,FMT=*) znuc, C6in, CN0in
        iel = nint(znuc)
        IF ( (iel < 1) .or. (iel > mxelem) .or. (abs(znuc-dble(iel)) > 1.d-6) ) then
          CALL VDW_ERROR('GETPARAM_D3','ERROR WHILE READING D3 COEF',error)
          return
       end IF
        IF ( (C6in <= 0.d0) .or. (CN0in <= 0.d0) ) then
          CALL VDW_ERROR('GETPARAM_D3','ERROR WHILE READING D3 COEF',error) 
          return
       end if
        !
        DO ityp = 1, ntyp
          IF ( ielem(ityp) == iel ) THEN
            icn = nCN0(ityp) + 1
            IF ( icn > mxCN0 ) THEN
              WRITE(6,'(/1x,"Too many CN0 values for element:",i4,/)') iel
              CALL VDW_ERROR('GETPARAM_D3','ERROR WHILE READING D3 COEF',error)
              return
            ENDIF
            nCN0(ityp) = icn
            C6a(ityp,icn,1) = C6in
            C6a(ityp,icn,2) = CN0in
          ENDIF
        ENDDO
      ENDDO
      !
    ENDIF
    GOTO 20
    !
 60 CONTINUE
    CALL VDW_ERROR('GETPARAM_D3','ERROR WHILE READING D3 COEF',error)
    return
    !
 80 CONTINUE
    !
  WRITE(6,'(/5x,"-------------------------------------", &
            /5x,"Parameters for Dispersion Correction:", &
            /5x,"-------------------------------------", &
            /5x,"  atom           C_6        CN0      ",/)')
  DO ityp = 1, ntyp
    do icn = 1, nCN0(ityp)
      WRITE(6,'(8x,a3,i3,5x,f8.3,4x,f8.4)') &
        namtyp(ityp), icn-1, C6a(ityp,icn,1), C6a(ityp,icn,2)
    ENDDO
  ENDDO
  !
  IF (ncoef > 0) THEN
    WRITE(6,'(/5x,"Number of additional D3 coefficients:",i3)') ncoef
    RETURN
 END IF
  !
END SUBROUTINE getparam_D3mod
!
!---------------------------------------------------------------------------
SUBROUTINE getparam_R0ab(ntyp, ielem, R0ab)
  !---------------------------------------------------------------------
  !
  ! Set cut-off radii. Tabulated parameters are in Angstrom.
  ! In parts due to Intel compiler bug.
  !
  use vdw_param,     only: vdw_error
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ntyp, ielem(ntyp)
  DOUBLE PRECISION, INTENT(OUT) :: R0ab(ntyp,ntyp)
  !
  ! and here the local variables
  !
  INTEGER, PARAMETER :: mxelem = 94
  DOUBLE PRECISION, PARAMETER :: Bohr = 0.52917726d0
  !
  INTEGER :: iel, jel, ityp, jtyp, k0, k
  DOUBLE PRECISION R0ab_param(4465)
  !
  include 'vdw_param-R0ab.inc'
  !
  DO jtyp = 1, ntyp
    jel = ielem(jtyp)
    IF ( (jel .le. 0) .OR. (jel .gt. mxelem) ) THEN
      WRITE(6,'(/1x,"No R0ab parameter for element:",i4)') jel
      STOP
    ENDIF
    k0 = ( jel*(jel-1) )/2
    !
    DO ityp = 1, ntyp
      iel = ielem(ityp)
      IF ( (iel .le. 0) .OR. (iel .gt. mxelem) ) THEN
        WRITE(6,'(/1x,"No R0ab parameter for element:",i4)') iel
        STOP
      ENDIF
      !
      IF (iel .le. jel) THEN
        k = k0 + iel
        R0ab(ityp,jtyp) = R0ab_param(k)/Bohr
        R0ab(jtyp,ityp) = R0ab_param(k)/Bohr
      ENDIF
      !
    ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE getparam_R0ab
!
!---------------------------------------------------------------------------
SUBROUTINE getparam_r2r4(ntyp, ielem, r2r4)
  !---------------------------------------------------------------------
  !
  use vdw_param,     only: vdw_error
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ntyp, ielem(ntyp)
  DOUBLE PRECISION, INTENT(OUT) :: r2r4(ntyp)
  !
  ! and here the local variables
  !
  INTEGER, PARAMETER :: mxelem = 94
  !
  INTEGER :: iel, ityp
  DOUBLE PRECISION r2r4_param(mxelem)
  !
  ! PBE0/def2-QZVP atomic values
  !
  !DATA r2r4_param(1:mxelem) / &
  !      8.0589,  3.4698, 29.0974, 14.8517, 11.8799,  7.8715,  5.5588, &
  !      4.7566,  3.8025,  3.1036, 26.1552, 17.2304, 17.7210, 12.7442, &
  !      9.5361,  8.1652,  6.7463,  5.6004, 29.2012, 22.3934, 19.0598, &
  !     16.8590, 15.4023, 12.5589, 13.4788, 12.2309, 11.2809, 10.5569, &
  !     10.1428,  9.4907, 13.4606, 10.8544,  8.9386,  8.1350,  7.1251, &
  !      6.1971, 30.0162, 24.4103, 20.3537, 17.4780, 13.5528, 11.8451, &
  !     11.0355, 10.1997,  9.5414,  9.0061,  8.6417,  8.9975, 14.0834, &
  !     11.8333, 10.0179,  9.3844,  8.4110,  7.5152, 32.7622, 27.5708, &
  !     23.1671, 21.6003, 20.9615, 20.4562, 20.1010, 19.7475, 19.4828, &
  !     15.6013, 19.2362, 17.4717, 17.8321, 17.4237, 17.1954, 17.1631, &
  !     14.5716, 15.8758, 13.8989, 12.4834, 11.4421, 10.2671,  8.3549, &
  !      7.8496,  7.3278,  7.4820, 13.5124, 11.6554, 10.0959,  9.7340, &
  !      8.8584,  8.0125, 29.8135, 26.3157, 19.1885, 15.8542, 16.1305, &
  !     15.6161, 15.1226, 16.1576 /
  !
  ! Scale r4/r2 values of the atoms by sqrt(Z); sqrt is also globally close
  ! to optimum; together with the factor 1/2 this yield reasonable C8 for
  ! He, Ne and Ar; for larger Z, C8 becomes too large, which effectively
  ! mimics higher R^n terms neglected due to stability reasons
  !
  ! r2r4 = sqrt( 0.5*sqrt(i)*r2r4_param(i) ) with i=elementnumber
  !
  ! The large number of digits is just to keep the results consistent with older
  ! versions. They should not imply any higher accuracy than the old values.
  !
  DATA r2r4_param(1:mxelem) / &
        2.00734898,  1.56637132,  5.01986934,  3.85379032,  3.64446594, &
        3.10492822,  2.71175247,  2.59361680,  2.38825250,  2.21522516, &
        6.58585536,  5.46295967,  5.65216669,  4.88284902,  4.29727576, &
        4.04108902,  3.72932356,  3.44677275,  7.97762753,  7.07623947, &
        6.60844053,  6.28791364,  6.07728703,  5.54643096,  5.80491167, &
        5.58415602,  5.41374528,  5.28497229,  5.22592821,  5.09817141, &
        6.12149689,  5.54083734,  5.06696878,  4.87005108,  4.59089647, &
        4.31176304,  9.55461698,  8.67396077,  7.97210197,  7.43439917, &
        6.58711862,  6.19536215,  6.01517290,  5.81623410,  5.65710424, &
        5.52640661,  5.44263305,  5.58285373,  7.02081898,  6.46815523, &
        5.98089120,  5.81686657,  5.53321815,  5.25477007, 11.02204549, &
       10.15679528,  9.35167836,  9.06926079,  8.97241155,  8.90092807, &
        8.85984840,  8.81736827,  8.79317710,  7.89969626,  8.80588454, &
        8.42439218,  8.54289262,  8.47583370,  8.45090888,  8.47339339, &
        7.83525634,  8.20702843,  7.70559063,  7.32755997,  7.03887381, &
        6.68978720,  6.05450052,  5.88752022,  5.70661499,  5.78450695, &
        7.79780729,  7.26443867,  6.78151984,  6.67883169,  6.39024318, &
        6.09527958, 11.79156076, 11.10997644,  9.51377795,  8.67197068, &
        8.77140725,  8.65402716,  8.53923501,  8.85024712 /
  !
  DO ityp = 1, ntyp
    iel = ielem(ityp)
    IF ( (iel .gt. 0) .AND. (iel .le. mxelem) ) THEN
      r2r4(ityp) = r2r4_param(iel)
    ELSE
      WRITE(6,'(/1x,"No r2r4 parameter for element:",i4)') iel
      STOP
    ENDIF
  ENDDO
  !
  RETURN
  !
END SUBROUTINE getparam_r2r4
!
!---------------------------------------------------------------------------
SUBROUTINE getparam_Rcov(ntyp, ielem, Rcov)
  !---------------------------------------------------------------------
  !
  use vdw_param,                only: vdw_error
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ntyp, ielem(ntyp)
  DOUBLE PRECISION, INTENT(OUT) :: Rcov(ntyp)
  !
  ! and here the local variables
  !
  INTEGER, PARAMETER :: mxelem = 94
  !
  INTEGER :: iel, ityp
  DOUBLE PRECISION Rcov_param(mxelem)
  !
  ! Covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
  ! Values for metals decreased by 10%
  !
  !DATA Rcov_param(1:mxelem) / &
  !      0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67, &
  !      1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54, &
  !      1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09, &
  !      1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39, &
  !      1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26, &
  !      1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57, &
  !      1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53, &
  !      1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32, &
  !      1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58, &
  !      1.52, 1.53, 1.54, 1.55 /
  !
  ! These new data are scaled with k2=4/3 and converted from Angstrom to
  ! atomic units via Bohr=0.52917726d0
  !
  DATA Rcov_param(1:mxelem) / &
        0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865, &
        1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527, &
        3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820, &
        2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730, &
        3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923, &
        2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188, &
        2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349, &
        2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216, &
        3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717, &
        2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967, &
        3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625, &
        4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657, &
        3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833, &
        3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098, &
        3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878, &
        2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790, &
        3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584, &
        3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289, &
        3.82984466, 3.85504098, 3.88023730, 3.90543362 /
  !
  DO ityp = 1, ntyp
    iel = ielem(ityp)
    IF ( (iel .gt. 0) .AND. (iel .le. mxelem) ) THEN
      Rcov(ityp) = Rcov_param(iel)
    ELSE
      WRITE(6,'(/1x,"No Rcov parameter for element:",i4)') iel
      STOP
    ENDIF
  ENDDO
  !
  RETURN
  !
END SUBROUTINE getparam_Rcov
!
!-----------------------------------------------------------------------
FUNCTION get_atom_number(nameat)
  !-----------------------------------------------------------------------
  !
  use vdw_param,                only: vdw_error
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(in) :: nameat
  INTEGER :: get_atom_number
  !
  ! and here the local variables
  !
  INTEGER :: iel
  INTEGER, PARAMETER :: ilower = ichar('A') - ichar('a')
  CHARACTER(len=1) :: nchar
  CHARACTER(len=2) :: atom, elements(103)
  !
  DATA elements(1:103) / &
     'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne','Na','Mg', &
     'Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr', &
     'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &
     'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd', &
     'In','Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
     'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf', &
     'Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po', &
     'At','Rn','Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm', &
     'Bk','Cf','Es','Fm','Md','No','Lr' /
  !
  !  Make first character a capital letter (only letters allowed).
  !
  atom = '  '
  nchar = nameat(1:1)
  IF (lge(nchar,'a') .AND. lle(nchar,'z')) THEN
    nchar = char(ichar(nchar) + ilower)
  ELSE IF (llt(nchar,'A') .OR. lgt(nchar,'Z')) THEN
    WRITE(6,*) 'Invalid element name:', nameat
    STOP 999
  ENDIF
  atom(1:1) = nchar
  !
  IF ( len(nameat) > 1 ) THEN
    !
    !  Make sure that second letter is lower case (all characters allowed).
    !
    nchar = nameat(2:2)
    IF (lge(nchar,'A') .AND. lle(nchar,'Z')) nchar = char(ichar(nchar) - ilower)
    atom(2:2) = nchar
  ENDIF
  !
  DO iel = 1, 103
    IF ( atom == elements(iel) ) GOTO 50
  ENDDO
  atom(2:2) = ' '
  DO iel = 1, 103
    IF ( atom == elements(iel) ) GOTO 50
  ENDDO
  WRITE(6,'(/1x,"Atom ",a2," not found")')  atom
  STOP 999
  !
50  CONTINUE
  get_atom_number = iel
  !
  RETURN
  !
END FUNCTION get_atom_number
FUNCTION read_input(iunit,label)
  INTEGER :: iunit,read_input,start
  character(len=*) :: label
  character(len=255) :: string
  read_input=1
  
  rewind(iunit)
20 continue
  read(iunit,end=30,err=30,fmt='(A255)') string
  start=index(string,label)
  if (start .ne. 0) then 
     read_input=0
     goto 30
  else
     goto 20
  end if

30 continue
  return
end FUNCTION read_input
