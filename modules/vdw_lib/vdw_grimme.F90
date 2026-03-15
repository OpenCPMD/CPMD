MODULE vdw_calculator
IMPLICIT NONE
PRIVATE
PUBLIC :: vdw_forces

CONTAINS
!---------------------------------------------------------------------------
! Compute dispersion contribution to atomic forces and stress tensor
!---------------------------------------------------------------------------
!
SUBROUTINE vdw_forces(alat, avec, bvec, celvol, natcel, idx_ityp, coorat, &
                      evdw, forces, stress, mepos, nproc, allgrp)
  !---------------------------------------------------------------------
  !
  USE vdw_param, ONLY: version, rcut2_cn, cn_k1, Rcov, rcut2_vdw, &
                       vdw_dir, vdw_pair, s6_D2, beta, C6_ij, R_sum, &
                       s6, s8, rs6, rs8, alp6, alp8, R0ab, r2r4
!For PWscf only:
!#if defined __PARA
!  USE mp_global, ONLY: me_image, nproc_image, intra_image_comm
!  USE mp,        ONLY: mp_sum
!#endif
#ifdef PARALLEL
  USE mpi_f08
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: natcel, idx_ityp(natcel)
  ! natcel   : number of atoms
  ! idx_ityp : type of each atom
  INTEGER, INTENT(IN) :: mepos, nproc
  type(MPI_COMM), intent(in) :: allgrp
  ! parallelization stuff
  DOUBLE PRECISION, INTENT(IN) :: alat, avec(3,3), bvec(3,3), celvol, &
                                  coorat(3,natcel)
  ! alat   : the cell parameter (in bohr)
  ! avec   : real space lattice vectors (alat units)
  ! bvec   : reciprocal lattice vectors (alat units)
  ! celvol : volume of unit cell
  ! coorat : atomic positions in Cartesian coordinates (alat units)
  DOUBLE PRECISION, INTENT(OUT) :: evdw, forces(3,natcel), stress(3,3)
  ! evdw   : the dispersion energy
  ! forces : the dispersion forces
  ! stress : the dispersion stress
  !
  ! and here the local variables
  !
  INTEGER :: i, ierr, iat, jat, kat, ityp, jtyp, k, it, mxtvec, nrt, nt1, nt2, nt3
  ! iat,jat   : atom counters
  ! ityp,jtyp : type of atom iat and jat
  ! it        : lattice vector counter
  ! mxtvec    : upper bound for number of T-vectors
  ! nrt       : actual number of vectors with tvec^2 < rcut_vdw^2
  ! nt1,nt2,nt3 : size of box containing sphere with radius 'rcut_vdw'
  INTEGER :: jat_first, jat_last, jat_div, jat_rest, &
             revcnt(0:nproc-1), displ(0:nproc-1), jat_nproc(0:nproc-1,2)
  ! parallelization stuff
  DOUBLE PRECISION, ALLOCATABLE :: tvec(:,:), tnrm2(:), CN(:), dE_dCN_sum(:)
  ! tvec   : Cartesian coordinates of ordered T-vectors (alat units)
  ! tnrm2  : length squared of ordered T-vectors (alat units)
  ! CN     : coordination number of all atoms in unit cell
  DOUBLE PRECISION :: frc_c6(3,natcel), str_c6(3,3), dcor(3), ds(3), &
                      tt(3), ttn, ttn2, dist1, dist2, dist4, dist6, &
                      dist8, c6, c8, dc6i, dc6j, arg, tmp, betaRS, &
                      expval, damp, efac6, efac8, efac, cfac, sfac, &
                      ffac6, ffac8, ffac, scal
  ! dcor   : vector R_ij between pair of atoms within unit cell (alat units)
  ! tt     : vector R_ij between current pair of atoms (alat units)
  ! dist1  : distance R_ij between the current pair of atoms (bohr)
  ! dist2  : distance**2
  ! dist6  : distance**6
  ! damp   : damping function
  !
!  DOUBLE PRECISION :: dE_dC6_sum(natcel,natcel), frc2(3,natcel), str2(3,3)
  !
!#if defined __PARA   ! --> PWSCF only
  !
  ! parallelization: divide atoms across processors of this image
  ! (different images have different atomic positions)
  !
!  jat_rest = MOD( natcel, nproc_image )
!  jat_div = natcel / nproc_image
!  IF (me_image+1 <= jat_rest) THEN
!    jat_first = (jat_div + 1) * me_image + 1
!    jat_last  = (jat_div + 1) * (me_image + 1)
!  ELSE
!    jat_first = (jat_div + 1) * jat_rest + jat_div * (me_image - jat_rest) + 1
!    jat_last  = (jat_div + 1) * jat_rest + jat_div * (me_image - jat_rest + 1)
!  ENDIF
!#else
!  jat_first = 1
!  jat_last = natcel
!#endif
  !
#ifdef PARALLEL
  !
  ! parallelization: divide atoms across nproc:
  !
  jat_rest = MOD( natcel, nproc )
  jat_div = natcel / nproc
  jat_nproc(0:nproc-1,1:2) = 1
  k = 0
  DO i = 0, nproc-1
    IF (i .gt. 0) THEN
      jat_nproc(i,1) = jat_nproc(i-1,2) + 1
    ELSE
      jat_nproc(i,1) = 1
    ENDIF
    jat_nproc(i,2) = jat_nproc(i,1) + jat_div - 1
    IF (k .lt. jat_rest) THEN
      k = k+1
      jat_nproc(i,2) = jat_nproc(i,2) + 1
    ENDIF
  ENDDO
  !
  DO i = 0, nproc-1
    IF (i .gt. 0) THEN
      displ(i) = jat_nproc(i-1,2)
    ELSE
      displ(i) = 0
    ENDIF
    revcnt(i) = jat_nproc(i,2) - jat_nproc(i,1) + 1
  ENDDO
  jat_first = jat_nproc(mepos,1)
  jat_last  = jat_nproc(mepos,2)
#else
  jat_first=1
  jat_last=natcel
#endif
  !
  evdw = 0.d0
  forces(1:3,1:natcel) = 0.d0
  stress(1:3,1:3) = 0.d0
  !
  ! Get coordination numbers
  !
  IF (version(1:2) == 'D3') THEN
    !
    ALLOCATE ( CN(natcel), dE_dCN_sum(natcel) )
    !
    CALL get_CN(alat, avec, bvec, natcel, idx_ityp, coorat, CN, &
                jat_first, jat_last)
#ifdef PARALLEL
    CALL MPI_Allgatherv(mpi_in_place, 0, MPI_DATATYPE_NULL, &
                        CN, revcnt, displ, mpi_double_precision, allgrp, ierr)
#endif
    !
    dE_dCN_sum(1:natcel) = 0.d0
!    dE_dC6_sum(1:natcel,1:natcel) = 0.d0
    !
  ENDIF
  !
  ! Get lattice vectors T (sphere 'rcut_vdw')
  !
  CALL get_mxtvec(bvec, rcut2_vdw, mxtvec, nt1, nt2, nt3, vdw_dir)
!  WRITE(6,'(1x,"Upper bound for number of T-vectors (vdW):",i8)') mxtvec
  !
  ALLOCATE ( tvec(3,mxtvec), tnrm2(mxtvec) )
  !
  CALL setup_tvec(avec, rcut2_vdw, nt1, nt2, nt3, nrt, tvec, tnrm2)
  !
  DO jat = jat_first, jat_last
    jtyp = idx_ityp(jat)
    !
!$OMP parallel do private(iat,ityp,dcor,i,k,ds,betaRs,c6,dc6i,dc6j,c8,it,tt,ttn2,scal, &
!$OMP                     dist1,dist2,dist4,dist6,arg,expval,damp,efac,ffac,sfac,dist8, &
!$OMP                     tmp,efac6,ffac6) reduction(+:evdw,forces,stress,de_dcn_sum)
    DO iat = 1, natcel
      ityp = idx_ityp(iat)
      !
      IF (vdw_pair(ityp,jtyp) == 0) CYCLE
      !
      dcor(1:3) = coorat(1:3,jat) - coorat(1:3,iat)
      !
      ! bring 'dcor' into the unit cell centered on the origin (if periodic
      ! boundary conditions are used). This prevents trouble, if atoms are
      ! displaced far away from origin
      ! (remember that translational invariance allows this!)
      !
      DO k = 1, 3
        ds(k) = dcor(1)*bvec(1,k) + dcor(2)*bvec(2,k) + dcor(3)*bvec(3,k)
      ENDDO
      DO k = 1, 3
        IF (vdw_dir(k) == 1) ds(k) = ds(k) - anint(ds(k))
      ENDDO
      DO k = 1, 3
        dcor(k) = ds(1)*avec(k,1) + ds(2)*avec(k,2) + ds(3)*avec(k,3)
      ENDDO
      !
      IF (version(1:2) == 'D2') THEN
        !
        betaRs = beta / R_sum(ityp,jtyp)
        !
      ELSE IF (version(1:2) == 'D3') THEN
        ! Get C6 and C8 coefficients and the derivative dC6/dCN
        CALL get_dC6_dCN(ityp, jtyp, CN(iat), CN(jat), c6, dc6i, dc6j)
        ! Note: r2r4 (Q) is stored as sqrt
        c8 = 3.0d0*c6*r2r4(ityp)*r2r4(jtyp)
      ENDIF
      !
      ! do we need to sort 'tvec+dcor' ????
      !
      DO it = nrt, 1, -1
!      DO it = 1, nrt
        !
        tt(1:3) = tvec(1:3,it) + dcor(1:3)
        ttn2 = tt(1)*tt(1) + tt(2)*tt(2) + tt(3)*tt(3)
!        IF ( (ttn2 > rcut2_vdw) .OR. (abs(ttn2) < 1.d-8) ) CYCLE
        IF ( (abs(ttn2) < 1.d-8) .OR. &
             ( (rcut2_vdw > 0.d0) .AND. (ttn2 > rcut2_vdw) ) ) CYCLE
        !
        ! store grad_Rj(ttn) in tt(1:3)
        !
        ttn = sqrt(ttn2)
        scal = 1.d0 / ttn
        tt(1:3) = scal * tt(1:3)
        !
        dist1 = alat*ttn
        dist2 = alat*alat*ttn2
        dist4 = dist2*dist2
        dist6 = dist2*dist4
        !
        IF (version(1:2) == 'D2') THEN
          !
          arg = beta*( dist1 / R_sum(ityp,jtyp) - 1.D0 )
          expval = exp(-arg)
          damp = 1.d0 / ( 1.d0 + expval )
          efac = s6_D2 * damp * ( C6_ij(ityp,jtyp) / dist6 )
          ffac = efac * ( betaRs*expval*damp - 6.d0/dist1 )
          sfac = ffac * dist1
          !
          evdw = evdw - efac
          IF (iat .ne. jat) forces(1:3,jat) = forces(1:3,jat) + ffac * tt(1:3)
          DO k = 1, 3
          DO i = 1, k
            stress(i,k) = stress(i,k) + sfac * tt(i) * tt(k)
          ENDDO
          ENDDO
          !
        ELSE IF (version(1:2) == 'D3') THEN
          !
          dist8 = dist2*dist6
          !
          IF (version(7:8) == 'BJ') THEN
            !
            ! DFT-D3 BJ damp
            ! Becke-Johnson parameters: a1=rs6, a2=rs8
            !
            tmp = rs6*sqrt(c8/c6) + rs8
            damp = 1.d0 / ( dist6 + tmp**6 )
            efac6 = s6 * c6 * damp
            ffac6 = -6.d0 * dist4 * dist1 * efac6 * damp
            !
            damp = 1.d0 / ( dist8 + tmp**8 )
            efac8 = s8 * c8 * damp
            ffac8 = -8.d0 * dist6 * dist1 * efac8 * damp
            !
            efac = efac6 + efac8
            ffac = ffac6 + ffac8
            sfac = ffac * dist1
            !
            dE_dCN_sum(jat) = dE_dCN_sum(jat) + efac * (dc6j / c6)
!            dE_dC6_sum(iat,jat) = dE_dC6_sum(iat,jat) + efac / c6
            !
            evdw = evdw - efac
            IF (iat .ne. jat) forces(1:3,jat) = forces(1:3,jat) + ffac * tt(1:3)
            DO k = 1, 3
            DO i = 1, k
              stress(i,k) = stress(i,k) + sfac * tt(i) * tt(k)
            ENDDO
            ENDDO
            !
          ELSE IF (version(7:7) == 'Z') THEN
            !
            ! DFT-D3 zero-damp
            !
            arg = R0ab(ityp,jtyp) / dist1
            tmp = 6.d0 * ( rs6 * arg )**alp6
            damp = 1.d0 / ( 1.d0 + tmp )
            efac6 = s6 * damp * ( c6 / dist6 )
            ffac6 = efac6 * ( alp6 * tmp * damp - 6.d0) / dist1
            !
            tmp = 6.d0 * ( rs8 * arg )**alp8
            damp = 1.d0 / ( 1.d0 + tmp )
            efac8 = s8 * damp * ( c8 / dist8 )
            ffac8 = efac8 * ( alp8 * tmp * damp - 8.d0) / dist1
            !
            efac = efac6 + efac8
            ffac = ffac6 + ffac8
            sfac = ffac * dist1
            !
            dE_dCN_sum(jat) = dE_dCN_sum(jat) + efac * (dc6j / c6)
!            dE_dC6_sum(iat,jat) = dE_dC6_sum(iat,jat) + efac / c6
            !
            evdw = evdw - efac
            IF (iat .ne. jat) forces(1:3,jat) = forces(1:3,jat) + ffac * tt(1:3)
            DO k = 1, 3
            DO i = 1, k
              stress(i,k) = stress(i,k) + sfac * tt(i) * tt(k)
            ENDDO
            ENDDO
            !
          ELSE
            STOP ! Unknown damping function
          ENDIF
          !
        ELSE
          STOP ! Unkown Grimme scheme
        ENDIF
        !
      ENDDO
      !
    ENDDO
    !
  ENDDO
  !
  ! Add contribution from dC6/dR (D3 only)
  !
  IF (version(1:2) == 'D3') THEN
    !
    DEALLOCATE ( tvec, tnrm2 )
    !
    ! Get lattice vectors T (sphere 'rcut_cn')
    !
    CALL get_mxtvec(bvec, rcut2_cn, mxtvec, nt1, nt2, nt3, vdw_dir)
!    WRITE(6,'(1x,"Upper bound for number of T-vectors (CN): ",i8)') mxtvec
    !
    ALLOCATE ( tvec(3,mxtvec), tnrm2(mxtvec) )
    !
    CALL setup_tvec(avec, rcut2_cn, nt1, nt2, nt3, nrt, tvec, tnrm2)
    !
    frc_c6(1:3,1:natcel) = 0.d0
    str_c6(1:3,1:3) = 0.d0
    !
!#if defined __PARA   ! --> PWSCF only
!    CALL mp_sum(dE_dCN_sum, intra_image_comm)
!#endif
#ifdef PARALLEL
    CALL MPI_Allgatherv(mpi_in_place, 0, MPI_DATATYPE_NULL, &
                        dE_dCN_sum, revcnt, displ, mpi_double_precision, allgrp, ierr)
#endif
    !
    DO jat = jat_first, jat_last
      jtyp = idx_ityp(jat)
      !
!$OMP parallel do private(iat,ityp,dcor,i,k,ds,it,tt,ttn2,scal,dist1,dist2,tmp,arg, &
!$OMP                     expval,damp,cfac,ffac,sfac) reduction(+:frc_c6,str_c6)
      DO iat = 1, natcel
        ityp = idx_ityp(iat)
        !
        IF (vdw_pair(ityp,jtyp) == 0) CYCLE
        !
        dcor(1:3) = coorat(1:3,jat) - coorat(1:3,iat)
        !
        ! bring 'dcor' into the unit cell centered on the origin (if periodic
        ! boundary conditions are used). This prevents trouble, if atoms are
        ! displaced far away from origin
        ! (remember that translational invariance allows this!)
        !
        DO k = 1, 3
          ds(k) = dcor(1)*bvec(1,k) + dcor(2)*bvec(2,k) + dcor(3)*bvec(3,k)
        ENDDO
        DO k = 1, 3
          IF (vdw_dir(k) == 1) ds(k) = ds(k) - anint(ds(k))
        ENDDO
        DO k = 1, 3
          dcor(k) = ds(1)*avec(k,1) + ds(2)*avec(k,2) + ds(3)*avec(k,3)
        ENDDO
        !
        ! do we need to sort 'tvec+dcor' ????
        !
        DO it = nrt, 1, -1
!        DO it = 1, nrt
          !
          tt(1:3) = tvec(1:3,it) + dcor(1:3)
          ttn2 = tt(1)*tt(1) + tt(2)*tt(2) + tt(3)*tt(3)
          IF ( (ttn2 > rcut2_cn) .OR. (abs(ttn2) < 1.d-8) ) CYCLE
          !
          ! store grad_Rj(ttn) in tt(1:3)
          !
          ttn = sqrt(ttn2)
          scal = 1.d0 / ttn
          tt(1:3) = scal * tt(1:3)
          !
          dist1 = alat*ttn
          dist2 = alat*alat*ttn2
          !
          tmp = Rcov(ityp) + Rcov(jtyp)
          arg = cn_k1*( tmp / dist1 - 1.d0 )
          expval = exp(-arg)
          damp = 1.d0 / ( 1.d0 + expval )
          cfac = -cn_k1 * ( tmp / dist2 ) * expval * damp * damp
          !
          ffac = ( dE_dCN_sum(iat) + dE_dCN_sum(jat) ) * cfac
          sfac = ffac * dist1
          !
          IF (iat .ne. jat) frc_c6(1:3,jat) = frc_c6(1:3,jat) + ffac * tt(1:3)
          DO k = 1, 3
          DO i = 1, k
            str_c6(i,k) = str_c6(i,k) + sfac * tt(i) * tt(k)
          ENDDO
          ENDDO
          !
        ENDDO
        !
      ENDDO
      !
    ENDDO
    !
    forces(1:3,1:natcel) = forces(1:3,1:natcel) + frc_c6(1:3,1:natcel)
    stress(1:3,1:3) = stress(1:3,1:3) + str_c6(1:3,1:3)
    !
!begin debug
!       call debug_dC6_num(alat, avec, bvec, natcel, idx_ityp, &
!                          coorat, dE_dC6_sum, frc2, str2)
!       call debug_dCN_num(alat, avec, bvec, natcel, idx_ityp, &
!                          coorat, dE_dCN_sum, frc2, str2)
!       write(6,9400)
! 9400  format(//1x,'VdW forces: check dC6/dR contribution',/)
!       do iat = 1, natcel
!         write(6,9420) iat, frc2(1:3,iat), frc_c6(1:3,iat)
!       enddo
! 9420  format(i5,3f12.6,5x,3f12.6)
!       write(6,9500)
! 9500  format(//1x,'VdW stress: check dC6/dR contribution',/)
!       do k = 1, 3
!         write(6,9520) str2(1:3,k), str_c6(1:k,k)
!       enddo
! 9520  format(3f12.6,5x,3f12.6)
!end debug
    !
  ENDIF
  !
  evdw = 0.5d0 * evdw
  !
  scal = 0.5d0 / celvol
  DO k = 1, 3
  DO i = 1, k
    stress(i,k) = scal * stress(i,k)
  ENDDO
  ENDDO
  !
  DO k = 1, 2
  DO i = k+1, 3
    stress(i,k) = stress(k,i)
  ENDDO
  ENDDO
!#if defined __PARA   ! --> PWSCF only
!  CALL mp_sum(evdw,   intra_image_comm)
!  CALL mp_sum(forces, intra_image_comm)
!  CALL mp_sum(stress, intra_image_comm)
!#endif
!#ifdef PARALLEL
!  CALL my_sum_d(evdw, 1, allgrp)
!!in CPMD: done in subroutine 'noforce/forces'
!!  CALL my_sum_d(forces, natcel*3, allgrp)
!!in CPMD: done in subroutine 'totstr'
!!  CALL my_sum_d(stress, 9, allgrp)
!#endif
  !
  DEALLOCATE ( tvec, tnrm2 )
  IF (version(1:2) == 'D3') DEALLOCATE ( CN, dE_dCN_sum )
  !
  RETURN
  !
END SUBROUTINE vdw_forces
!
!-----------------------------------------------------------------------
SUBROUTINE get_CN(alat, avec, bvec, natcel, idx_ityp, coorat, CN, &
                  jat_first, jat_last)
  !---------------------------------------------------------------------
  !
  ! Compute coordination numbers by adding an inverse damping function
  !
  ! CN(j) = \sum_i  1 / [ 1 + exp( -k1 * ( (R_i^0+R_j^0) / R_ij - 1 ) ) ]
  !
  USE vdw_param, ONLY: rcut2_cn, cn_k1, Rcov, vdw_dir
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: natcel, idx_ityp(natcel), jat_first, jat_last
  ! natcel   : number of atoms
  ! idx_ityp : type of each atom
  DOUBLE PRECISION, INTENT(IN) :: alat, avec(3,3), bvec(3,3), coorat(3,natcel)
  ! alat   : the cell parameter (in bohr)
  ! avec   : real space lattice vectors (alat units)
  ! bvec   : reciprocal lattice vectors (alat units)
  ! coorat : atomic positions in Cartesian coordinates (alat units)
  DOUBLE PRECISION, INTENT(OUT) :: CN(natcel)
  ! CN : coordination number of all atoms in unit cell
  !
  ! and here the local variables
  !
  INTEGER :: iat, jat, ityp, jtyp, k, it, mxtvec, nrt, nt1, nt2, nt3
  ! iat,jat   : atom counters
  ! ityp,jtyp : type of atom iat and jat
  ! it        : lattice vector counter
  ! mxtvec    : upper bound for number of T-vectors
  ! nrt       : actual number of vectors with tvec^2 < rcut_cn^2
  ! nt1,nt2,nt3 : size of box containing sphere with radius 'rcut_cn'
  DOUBLE PRECISION, ALLOCATABLE :: tvec(:,:), tnrm2(:)
  ! tvec   : Cartesian coordinates of ordered T-vectors (alat units)
  ! tnrm2  : length squared of ordered T-vectors (alat units)
  DOUBLE PRECISION :: dcor(3), ds(3), tt(3), ttn2, dist, cnsum, arg, damp
  ! dcor   : vector R_ij between pair of atoms within unit cell (alat units)
  ! tt     : vector R_ij between current pair of atoms (alat units)
  ! dist   : distance R_ij between the current pair of atoms (bohr)
  !
  CALL get_mxtvec(bvec, rcut2_cn, mxtvec, nt1, nt2, nt3, vdw_dir)
!  WRITE(6,'(1x,"Upper bound for number of T-vectors (CN): ",i8)') mxtvec
  !
  ALLOCATE ( tvec(3,mxtvec), tnrm2(mxtvec) )
  !
  CALL setup_tvec(avec, rcut2_cn, nt1, nt2, nt3, nrt, tvec, tnrm2)
  !
  DO jat = jat_first, jat_last
    jtyp = idx_ityp(jat)
    !
    cnsum = 0.d0
    !
!$OMP parallel do private(iat,ityp,dcor,k,ds,it,tt,ttn2,dist,arg,damp) reduction(+:cnsum)
    DO iat = 1, natcel
      ityp = idx_ityp(iat)
      !
      dcor(1:3) = coorat(1:3,jat) - coorat(1:3,iat)
      !
      ! bring 'dcor' into the unit cell centered on the origin (if periodic
      ! boundary conditions are used). This prevents trouble, if atoms are
      ! displaced far away from origin
      ! (remember that translational invariance allows this!)
      !
      DO k = 1, 3
        ds(k) = dcor(1)*bvec(1,k) + dcor(2)*bvec(2,k) + dcor(3)*bvec(3,k)
      ENDDO
      DO k = 1, 3
        IF (vdw_dir(k) == 1) ds(k) = ds(k) - anint(ds(k))
      ENDDO
      DO k = 1, 3
        dcor(k) = ds(1)*avec(k,1) + ds(2)*avec(k,2) + ds(3)*avec(k,3)
      ENDDO
      !
      ! do we need to sort 'tvec+dcor' ????
      !
      DO it = nrt, 1, -1
!      DO it = 1, nrt
        !
        tt(1:3) = tvec(1:3,it) + dcor(1:3)
        ttn2 = tt(1)*tt(1) + tt(2)*tt(2) + tt(3)*tt(3)
        IF ( (ttn2 <= rcut2_cn) .and. (abs(ttn2) > 1.d-8) ) THEN
          dist = alat*sqrt(ttn2)
          arg = cn_k1*( ( Rcov(ityp) + Rcov(jtyp) ) / dist - 1.d0 )
          damp = 1.d0/( 1.d0 + exp(-arg) )
          cnsum = cnsum + damp
        ENDIF
        !
      ENDDO
      !
    ENDDO
    !
    CN(jat) = cnsum
    !
  ENDDO
  !
  DEALLOCATE ( tvec, tnrm2 )
  !
  RETURN
  !
END SUBROUTINE get_CN
!
!-----------------------------------------------------------------------
SUBROUTINE get_C6(ityp, jtyp, cni, cnj, c6)
  !---------------------------------------------------------------------
  !
  ! Interpolate C6 coefficient.
  !
  ! Warning: interpolation via exponentials is sensitive to numerics
  !          when cni (cnj) is much larger than CN0
  !
  USE vdw_param, ONLY: version, cn_k3, nCN0, C6ab, C6a
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ityp, jtyp
  ! itype,jtyp : type of both atoms
  DOUBLE PRECISION, INTENT(IN) :: cni, cnj
  ! cni,cnj : coordination of both atoms
  DOUBLE PRECISION, INTENT(OUT) :: c6
  ! c6 : coordination-dependend C6 coefficient
  !
  ! and here the local variables
  !
  INTEGER :: icn, jcn
  DOUBLE PRECISION :: c6i, c6j, c6tmp, dist, fexp, &
                      csave, dsave, wsum, zsum
  !
  IF (version(3:5) == 'mod') THEN
    !
    csave = C6a(ityp,1,1)
    dsave = 1.d+99
    zsum = 0.d0
    wsum = 0.d0
    DO icn = 1, nCN0(ityp)
      c6tmp = C6a(ityp,icn,1)
      dist = (cni - C6a(ityp,icn,2))**2
      IF (dist < dsave) THEN
        dsave = dist
        csave = c6tmp
      ENDIF
      fexp = exp(-cn_k3 * dist)
      zsum = zsum + fexp * c6tmp
      wsum = wsum + fexp
    ENDDO
    IF (wsum > 1.d-99) THEN
      c6i = zsum/wsum
    ELSE
      c6i = csave
    ENDIF
    !
    csave = C6a(jtyp,1,1)
    dsave = 1.d+99
    zsum = 0.d0
    wsum = 0.d0
    DO jcn = 1, nCN0(jtyp)
      c6tmp = C6a(jtyp,jcn,1)
      dist = (cnj - C6a(jtyp,jcn,2))**2
      IF (dist < dsave) THEN
        dsave = dist
        csave = c6tmp
      ENDIF
      fexp = exp(-cn_k3 * dist)
      zsum = zsum + fexp * c6tmp
      wsum = wsum + fexp
    ENDDO
    IF (wsum > 1.d-99) THEN
      c6j = zsum/wsum
    ELSE
      c6j = csave
    ENDIF
    !
    c6 = sqrt( c6i * c6j )
    !
  ELSE
    !
    csave = 0.d0
    dsave = 1.d+99
    zsum = 0.d0
    wsum = 0.d0
    DO jcn = 1, nCN0(jtyp)
    DO icn = 1, nCN0(ityp)
      c6tmp = C6ab(ityp,jtyp,icn,jcn,1)
      dist = (cni - C6ab(ityp,jtyp,icn,jcn,2))**2 + &
             (cnj - C6ab(ityp,jtyp,icn,jcn,3))**2
      IF (dist < dsave) THEN
        dsave = dist
        csave = c6tmp
      ENDIF
      fexp = exp(-cn_k3 * dist)
      zsum = zsum + fexp * c6tmp
      wsum = wsum + fexp
    ENDDO
    ENDDO
    IF (wsum > 1.d-99) THEN
      c6 = zsum/wsum
    ELSE
      c6 = csave
    ENDIF
    !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE get_C6
!
!-----------------------------------------------------------------------
SUBROUTINE get_dC6_dCN(ityp, jtyp, cni, cnj, c6, dc6i, dc6j)
  !---------------------------------------------------------------------
  !
  ! Gradient of C6 coefficient with respect to coordination number
  !
  ! Warning: interpolation via exponentials is sensitive to numerics
  !          when cni (cnj) is much larger than CN0
  !
  USE vdw_param, ONLY: version, cn_k3, nCN0, C6ab, C6a
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ityp, jtyp
  ! itype,jtyp : type of both atoms
  DOUBLE PRECISION, INTENT(IN) :: cni, cnj
  ! cni,cnj : coordination of both atoms
  DOUBLE PRECISION, INTENT(OUT) :: c6, dc6i, dc6j
  ! c6 : coordination-dependend C6 coefficient
  !
  ! and here the local variables
  !
  INTEGER :: icn, jcn
  DOUBLE PRECISION :: c6i, c6j, c6tmp, cni_ref, cnj_ref, dist, fexp, &
                      csave, dsave, wsum, zsum, dwsum_i, dwsum_j,    &
                      dzsum_i, dzsum_j, fac
  !
  IF (version(3:5) == 'mod') THEN
    !
    csave = C6a(ityp,1,1)
    dsave = 1.d+99
    zsum = 0.d0
    wsum = 0.d0
    dzsum_i = 0.d0
    dwsum_i = 0.d0
    DO icn = 1, nCN0(ityp)
      c6tmp = C6a(ityp,icn,1)
      cni_ref = C6a(ityp,icn,2)
      dist = (cni - cni_ref)**2
      IF (dist < dsave) THEN
        dsave = dist
        csave = c6tmp
      ENDIF
      fexp = exp(-cn_k3 * dist)
      zsum = zsum + fexp * c6tmp
      wsum = wsum + fexp
      dzsum_i = dzsum_i - 2.d0*cn_k3*(cni - cni_ref) * fexp * c6tmp
      dwsum_i = dwsum_i - 2.d0*cn_k3*(cni - cni_ref) * fexp
    ENDDO
    IF (wsum > 1.d-99) THEN
      c6i = zsum/wsum
      fac = 1.d0/(wsum*wsum)
      dc6i = fac*(wsum*dzsum_i - zsum*dwsum_i)
    ELSE
      c6i = csave
      dc6i = 0.d0
    ENDIF
    !
    csave = C6a(jtyp,1,1)
    dsave = 1.d+99
    zsum = 0.d0
    wsum = 0.d0
    dzsum_j = 0.d0
    dwsum_j = 0.d0
    DO jcn = 1, nCN0(jtyp)
      c6tmp = C6a(jtyp,jcn,1)
      cnj_ref = C6a(jtyp,jcn,2)
      dist = (cnj - cnj_ref)**2
      IF (dist < dsave) THEN
        dsave = dist
        csave = c6tmp
      ENDIF
      fexp = exp(-cn_k3 * dist)
      zsum = zsum + fexp * c6tmp
      wsum = wsum + fexp
      dzsum_j = dzsum_j - 2.d0*cn_k3*(cnj - cnj_ref) * fexp * c6tmp
      dwsum_j = dwsum_j - 2.d0*cn_k3*(cnj - cnj_ref) * fexp
    ENDDO
    IF (wsum > 1.d-99) THEN
      c6j = zsum/wsum
      fac = 1.d0/(wsum*wsum)
      dc6j = fac*(wsum*dzsum_j - zsum*dwsum_j)
    ELSE
      c6j = csave
      dc6j = 0.d0
    ENDIF
    !
    c6 = sqrt( c6i * c6j )
    dc6i = 0.5d0 * ( dc6i / c6 ) * c6j
    dc6j = 0.5d0 * c6i * ( dc6j / c6 )
    !
  ELSE
    !
    csave = 0.d0
    dsave = 1.d+99
    zsum = 0.d0
    wsum = 0.d0
    dzsum_i = 0.d0
    dwsum_i = 0.d0
    dzsum_j = 0.d0
    dwsum_j = 0.d0
    !
    DO jcn = 1, nCN0(jtyp)
    DO icn = 1, nCN0(ityp)
      c6tmp = C6ab(ityp,jtyp,icn,jcn,1)
      cni_ref = C6ab(ityp,jtyp,icn,jcn,2)
      cnj_ref = C6ab(ityp,jtyp,icn,jcn,3)
      dist = (cni - cni_ref)**2 + (cnj - cnj_ref)**2
      IF (dist < dsave) THEN
        dsave = dist
        csave = c6tmp
      ENDIF
      fexp = exp(-cn_k3 * dist)
      zsum = zsum + fexp*c6tmp
      wsum = wsum + fexp
      dzsum_i = dzsum_i - 2.d0*cn_k3*(cni - cni_ref) * fexp * c6tmp
      dwsum_i = dwsum_i - 2.d0*cn_k3*(cni - cni_ref) * fexp
      dzsum_j = dzsum_j - 2.d0*cn_k3*(cnj - cnj_ref) * fexp * c6tmp
      dwsum_j = dwsum_j - 2.d0*cn_k3*(cnj - cnj_ref) * fexp
    ENDDO
    ENDDO
    IF (wsum > 1.d-99) THEN
      c6 = zsum/wsum
      fac = 1.d0/(wsum*wsum)
      dc6i = fac*(wsum*dzsum_i - zsum*dwsum_i)
      dc6j = fac*(wsum*dzsum_j - zsum*dwsum_j)
    ELSE
      c6 = csave
      dc6i = 0.d0
      dc6j = 0.d0
    ENDIF
    !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE get_dC6_dCN
!
!-----------------------------------------------------------------------
SUBROUTINE get_mxtvec(bvec, rcut2, mxtvec, nt1, nt2, nt3, vdw_dir)
  !---------------------------------------------------------------------
  !
  ! Get a save estimate for the number of lattice vectors T with
  ! length smaller than 'rcut'.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: vdw_dir(3)
  ! vdw_dir : use periodic boundary conditions (=1) or not (=0)
  INTEGER, INTENT(out) :: mxtvec, nt1, nt2, nt3
  ! mxtvec  : upper bound for number of T-vectors
  ! nt1,nt2,nt3 : size of box containing sphere with radius 'rcut'
  DOUBLE PRECISION, INTENT(in) :: bvec(3,3), rcut2
  ! bvec  : reciprocal lattice vectors (alat units)
  ! rcut2 : cut-off radius squared (alat units)
  !
  ! and here the local variables
  !
  DOUBLE PRECISION :: blen2
  !
  ! Determine parallelepiped encompassing sphere with radius 'rcut'.
  !
  nt1 = 0
  IF (vdw_dir(1) == 1) THEN
    blen2 = bvec(1,1)*bvec(1,1) + bvec(2,1)*bvec(2,1) + bvec(3,1)*bvec(3,1)
    nt1 = int(sqrt(blen2*rcut2)) + 2
  ENDIF
  nt2 = 0
  IF (vdw_dir(2) == 1) THEN
    blen2 = bvec(1,2)*bvec(1,2) + bvec(2,2)*bvec(2,2) + bvec(3,2)*bvec(3,2)
    nt2 = int(sqrt(blen2*rcut2)) + 2
  ENDIF
  nt3 = 0
  IF (vdw_dir(3) == 1) THEN
    blen2 = bvec(1,3)*bvec(1,3) + bvec(2,3)*bvec(2,3) + bvec(3,3)*bvec(3,3)
    nt3 = int(sqrt(blen2*rcut2)) + 2
  ENDIF
  mxtvec = (2*nt1+1) * (2*nt2+1) * (2*nt3+1)
  !
  RETURN
  !
END SUBROUTINE get_mxtvec
!
!-----------------------------------------------------------------------
SUBROUTINE setup_tvec(avec, rcut2, nt1, nt2, nt3, nrt, tvec, tnrm2)
  !---------------------------------------------------------------------
  !
  ! Generates neighbours shells (cartesian, in units of lattice parameter)
  ! with length < rcut, and returns them in order of increasing length:
  !   tvec(:) = i*a1(:) + j*a2(:) + k*a3(:) ,  tnrm2 = tvec^2
  ! where a1, a2, a3 are primitive lattice vectors.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nt1, nt2, nt3
  ! nt1,nt2,nt3 : size of box containing sphere with radius 'rcut'
  INTEGER, INTENT(out) :: nrt
  ! nrt   : number of vectors with tvec^2 < rcut^2
  DOUBLE PRECISION, INTENT(in) :: avec(3,3), rcut2
  ! avec  : lattice vectors ( a1=avec(:,1), a2=avec(:,2), a3=avec(:,3) )
  ! rcut2 : cut-off radius squared (alat units)
  DOUBLE PRECISION, INTENT(out) :: tvec(3,*), tnrm2(*)
  !
  ! and here the local variables
  !
  INTEGER, ALLOCATABLE :: isort(:)
  INTEGER :: i, j, k, it, itswp, iswap
  DOUBLE PRECISION :: tt(3), ttn2
  !
  nrt = 1
  tvec(1:3,1) = 0.d0
  tnrm2(1) = 0.d0
  IF ( (nt1+nt2+nt3 == 0) .or. (rcut2 == 0.d0) ) RETURN
  !
  ! Collect all lattice vectors inside parallelepiped.
  !
  nrt = 0
  DO k = -nt3, nt3
    DO j = -nt2, nt2
      DO i = -nt1, nt1
        tt(1:3) = dble(i)*avec(1:3,1) + dble(j)*avec(1:3,2) + dble(k)*avec(1:3,3)
        ttn2 = tt(1)*tt(1) + tt(2)*tt(2) + tt(3)*tt(3)
        IF (ttn2 <= rcut2) THEN
          nrt = nrt + 1
!          IF (nrt > mxtvec) CALL errore('setup_tvec', 'too many T-vectors', nrt)
          tvec(1:3,nrt) = tt(1:3)
          tnrm2(nrt) = ttn2
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !
  ! Reorder the vectors in order of increasing length.
  ! Index array is initialized inside of sorting routine.
  !
  ALLOCATE ( isort(nrt) )
  !
  isort(1) = 0
  IF (nrt > 1) CALL hpsort(nrt, tnrm2, isort)
  DO it = 1, nrt-1
20  itswp = isort(it)
    IF (itswp .ne. it) THEN
      tt(1:3) = tvec(1:3,itswp)
      tvec(1:3,itswp) = tvec(1:3,isort(itswp))
      tvec(1:3,isort(itswp)) = tt(1:3)
      iswap = isort(it)
      isort(it) = isort(itswp)
      isort(itswp) = iswap
      GOTO 20
    ENDIF
  ENDDO
  !
  DEALLOCATE (isort)
  !
  RETURN
  !
END SUBROUTINE setup_tvec
!
!---------------------------------------------------------------------
subroutine hpsort (n, ra, ind)  
  !---------------------------------------------------------------------
  ! sort an array ra(1:n) into ascending order using heapsort algorithm.
  ! n is input, ra is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ra).
  ! in case of equal values in the data array (ra) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  implicit none  
  !-input/output variables
  integer :: n  
  integer :: ind (*)  
  double precision :: ra (*)  
  !-local variables
  integer :: i, ir, j, l, iind  
  double precision :: rra  
  ! initialize index array
  if (ind (1) .eq.0) then  
     do i = 1, n  
        ind (i) = i  
     enddo
  endif
  ! nothing to order
  if (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  
  ir = n  
10 continue  
  ! still in hiring phase
  if (l.gt.1) then  
     l = l - 1  
     rra = ra (l)  
     iind = ind (l)  
     ! in retirement-promotion phase.
  else  
     ! clear a space at the end of the array
     rra = ra (ir)  
     !
     iind = ind (ir)  
     ! retire the top of the heap into it
     ra (ir) = ra (1)  
     !
     ind (ir) = ind (1)  
     ! decrease the size of the corporation
     ir = ir - 1  
     ! done with the last promotion
     if (ir.eq.1) then  
        ! the least competent worker at all !
        ra (1) = rra  
        !
        ind (1) = iind  
        return  
     endif
  endif
  ! wheter in hiring or promotion phase, we
  i = l  
  ! set up to place rra in its proper level
  j = l + l  
  !
  do while (j.le.ir)  
     if (j.lt.ir) then  
        ! compare to better underling
        if (ra (j) .lt.ra (j + 1) ) then  
           j = j + 1  
        elseif (ra (j) .eq.ra (j + 1) ) then  
           if (ind (j) .lt.ind (j + 1) ) j = j + 1  
        endif
     endif
     ! demote rra
     if (rra.lt.ra (j) ) then  
        ra (i) = ra (j)  
        ind (i) = ind (j)  
        i = j  
        j = j + j  
     elseif (rra.eq.ra (j) ) then  
        ! demote rra
        if (iind.lt.ind (j) ) then  
           ra (i) = ra (j)  
           ind (i) = ind (j)  
           i = j  
           j = j + j  
        else  
           ! set j to terminate do-while loop
           j = ir + 1  
        endif
        ! this is the right place for rra
     else  
        ! set j to terminate do-while loop
        j = ir + 1  
     endif
  enddo
  ra (i) = rra  
  ind (i) = iind  
  goto 10  
  !
end subroutine hpsort

END MODULE vdw_calculator
