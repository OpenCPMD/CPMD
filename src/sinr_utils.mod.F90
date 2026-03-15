MODULE sinr_utils
USE kinds,                      ONLY : real_8
USE error_handling,             ONLY : stopgm
USE nose,                       ONLY : dtsuz, &
                                       ncalls, &
                                       nit,nedof
USE fileopen_utils,             ONLY : fileopen
                                       
USE fileopenmod,                ONLY : fo_app
USE rmas,                       ONLY : rmass
USE store_types,                ONLY : restart1
USE prng_utils,                 ONLY : repprngu,&
                                       repprngu_vec
USE tpar,                       ONLY : dt_ions
USE ions,                       ONLY : ions0,&
                                       ions1
USE cnst,                       ONLY : pi, &
                                       kboltz, &
                                       factem

USE parac,                      ONLY: parai,&
                                      paral
USE system,                     ONLY : maxsys,iatpt,&
                                       cnti, &
                                       cntr


IMPLICIT NONE


REAL(real_8) :: T_target, sinr_mass_1,sinr_mass_2,lbylp1,gamma_sinr,t_ke,sigma_sinr,&
                esinr,emgt,em2gt,sqe2gt, lambda_sinr, sinr_ke, tau_sinr
REAL(real_8), ALLOCATABLE :: v1_sinr(:,:,:,:), v2_sinr(:,:,:,:)

INTEGER :: l_sinr
LOGICAL :: ferror,ke_eq_LkbT

!!...................................................................................................!!
CONTAINS

  SUBROUTINE invelp_sinr(velp,v1_sinr,v2_sinr)
    REAL(real_8) :: velp(:,:,:),v1_sinr(:,:,:,:),v2_sinr(:,:,:,:)
    INTEGER :: i,ia,is,j,ierr,k
    REAL(real_8) :: width, fact, alfa1,rnr(1),vsq

! INITIALIZE THERMOSTAT VARIABLES RANDOMLY
! SCALE VELOCITIES SO THAT ISOKINETIC CONSTRAINT IS STRICTLY FOLLOWED
 
    width=DSQRT(t_ke/sinr_mass_2)
    DO is=1,ions1%nsp
      DO ia=1,ions0%na(is)
        DO k=1,3
          
          DO j=1,l_sinr
            CALL repprngu_vec(1,rnr)
            alfa1=2.0_real_8*pi*rnr(1)
              v1_sinr(j,k,ia,is)=DSQRT(DLOG(repprngu())*(-2.0_real_8))*DCOS(alfa1)
              v2_sinr(j,k,ia,is)=DSQRT(DLOG(repprngu())*(-2.0_real_8))*DSIN(alfa1)
              v2_sinr(j,k,ia,is)=v2_sinr(j,k,ia,is)*width
          ENDDO
           vsq=velp(k,ia,is)**2.d0
           DO j=1,l_sinr
              vsq=vsq+v1_sinr(j,k,ia,is)**2.d0
           ENDDO

           vsq=DSQRT(vsq)

           velp(k,ia,is)=velp(k,ia,is)/vsq

           DO j=1,l_sinr
              v1_sinr(j,k,ia,is)=v1_sinr(j,k,ia,is)/vsq
           ENDDO

           velp(k,ia,is)=velp(k,ia,is)/(DSQRT(rmass%pma(is)/(REAL(l_sinr,kind=real_8)*t_ke)))
 
           DO j=1,l_sinr
              v1_sinr(j,k,ia,is)=v1_sinr(j,k,ia,is)/(DSQRT(sinr_mass_1/(REAL((l_sinr+1),kind=real_8)*t_ke)))
           ENDDO

        ENDDO
      ENDDO
    ENDDO


  END SUBROUTINE invelp_sinr

!!...............................................................................................!!

  SUBROUTINE sinr_init
    INTEGER :: i,ia,is,j,k
    CHARACTER(*), PARAMETER :: procedureN="sinralloc"
    INTEGER :: ierr

!--------------------- SET THERMOSTAT PARAMETERS-----------------------------!

     T_target=cntr%tempsinr     ! take temperature in kelvin from input

 !   gamma_sinr=0.00025_real_8     
     gamma_sinr=cntr%gammasinr   ! take gamma in atomic unit from input

 !   l_sinr=4 !! change no of thermostat to check
     l_sinr=cnti%lsinr   ! take L from input

     tau_sinr=cntr%tausinr    !take tau in atomic unit from input

!---------------------------------------------------------------------------!

    ALLOCATE(v1_sinr(l_sinr,3,maxsys%nax,maxsys%nsx),v2_sinr(l_sinr,3,maxsys%nax,maxsys%nsx), &
        stat=ierr)
    IF(ierr/=0) CALL stopgm(procedureN,'allocation problem',&
            __LINE__,__FILE__)

     
    lbylp1=REAL(l_sinr,kind=real_8)/REAL((l_sinr+1),kind=real_8)
    t_ke=T_target/factem   ! set target kinetic energy in hartree
    sinr_mass_1=t_ke*tau_sinr**2
    sinr_mass_2=sinr_mass_1
    
    lambda_sinr=REAL(l_sinr,kind=real_8)*t_ke
    sigma_sinr = DSQRT((2.d0 * gamma_sinr * t_ke)/sinr_mass_2)
    emgt=DEXP(-gamma_sinr*dt_ions)
    em2gt=DEXP(-2._real_8*gamma_sinr*dt_ions)
    sqe2gt=DSQRT((1._real_8-em2gt)/(2._real_8*gamma_sinr))

  END SUBROUTINE sinr_init

!!........................................................................................!!

  SUBROUTINE sinrup(velp,v1_sinr,v2_sinr)
    REAL(real_8), INTENT(INOUT) :: velp(:,:,:), v1_sinr(:,:,:,:), &
                                   v2_sinr(:,:,:,:)
  
    INTEGER :: it,isy,i,ia,is,k,j,N
    REAL(real_8) :: step
    REAL(real_8) :: G, H, sum_H, aa

!-----------------------i L_N dt Operator-------------------------------!
        
    DO it=1,nit
      DO isy=1,ncalls
        step=dtsuz(isy)!/REAL(nit,kind=real_8)!*n_inner_step
  
        DO is=1,ions1%nsp
          DO ia=1,ions0%na(is)
          DO k=1,3

!---------------Outer translation operator LN,2-----------------------------!

            DO j=1,l_sinr  ! thermostat chain
               G = (sinr_mass_1*v1_sinr(j,k,ia,is)**2.d0 - t_ke)/sinr_mass_2
               v2_sinr(j,k,ia,is)=v2_sinr(j,k,ia,is)+G * (step/4.d0)
            ENDDO

!---------------Inner translation operator LN,1-----------------------------!

            sum_H = rmass%pma(is)*velp(k,ia,is)**2
            DO j=1,l_sinr
              aa=DEXP(-v2_sinr(j,k,ia,is)*step*0.5_real_8)
              sum_H=sum_H+lbylp1*sinr_mass_1*v1_sinr(j,k,ia,is)**2.d0*aa*aa
            ENDDO
            H = DSQRT(lambda_sinr/sum_H)
            velp(k,ia,is)=velp(k,ia,is)*H
            DO j=1,l_sinr
              aa=DEXP(-v2_sinr(j,k,ia,is)*step*0.5_real_8)
              v1_sinr(j,k,ia,is)=v1_sinr(j,k,ia,is)*H*aa
            ENDDO

!---------------Outer translation operator LN,2-----------------------------!

            DO j=1,l_sinr
              G = (sinr_mass_1*v1_sinr(j,k,ia,is)**2.d0 - t_ke)/sinr_mass_2
              v2_sinr(j,k,ia,is)=v2_sinr(j,k,ia,is)+G * (step/4.d0)
            ENDDO
            
         ENDDO  ! end for k
       ENDDO    ! end for atom
      ENDDO
    ENDDO
    ENDDO
  END SUBROUTINE sinrup

!!................................................................................!!

  SUBROUTINE velupi_sinr(velp,v1_sinr,fion)
    REAL(real_8),INTENT(INOUT) :: velp(:,:,:), v1_sinr(:,:,:,:)
    REAL(real_8),INTENT(IN) :: fion(:,:,:)
    
    REAL(real_8) :: a,b,s,sdot,arg,root_b
    INTEGER :: i,ia,is,k,j,iat
 
!-----------------------i L_v dt Operator-------------------------------!

    DO iat=1,ions1%nat
       ia=iatpt(1,iat)
       is=iatpt(2,iat)
       DO k=1,3

!--------UPDATE SYSTEM VELOCITIES-------------!

        a=(fion(k,ia,is)*velp(k,ia,is))/lambda_sinr
        b=fion(k,ia,is)**2.d0/(rmass%pma(is)*lambda_sinr)

        root_b=DSQRT(b)
        arg=dt_ions*root_b/2.d0

        IF (arg<(1.D-5))THEN
           s=(1.d0/root_b) * sinh_limit(arg)+(a/b) * (cosh_limit(arg)-1.d0)
           sdot=cosh_limit(arg)+(a/root_b) * sinh_limit(arg)
        ELSE
           s=(1.d0/root_b) * DSINH(arg) + (a/b) * (DCOSH(arg)-1.d0)
           sdot=DCOSH(arg) + (a/root_b) * DSINH(arg)
        ENDIF

        velp(k,ia,is) = (velp(k,ia,is)+(fion(k,ia,is)/rmass%pma(is))*s)/sdot

!-----------UPDATE THERMOSTATS----------------!

        DO j=1,l_sinr
          v1_sinr(j,k,ia,is)=v1_sinr(j,k,ia,is)/sdot
        ENDDO
  
     ENDDO
   ENDDO
 
  END SUBROUTINE velupi_sinr

!!...................................................................!!
    FUNCTION sinh_limit(x)result(sinh_l)

    REAL(8) :: x
    REAL(8) :: sinh_l

    sinh_l=x+(x**3.D0)/6.D0

    END FUNCTION sinh_limit
!!...................................................................!!
    FUNCTION cosh_limit(x)result(cosh_l)

    REAL(8) :: x
    REAL(8) :: cosh_l

    cosh_l=1.d0+(x**2.D0)/2.D0+(x**4.d0)/24.d0
    
    END FUNCTION cosh_limit

!!..............................................................................!!

  SUBROUTINE posupi_sinr(tau0,taup,velp,v2_sinr)
    REAL(real_8) :: tau0(:,:,:),velp(:,:,:),taup(:,:,:),v2_sinr(:,:,:,:)
    REAL(real_8) :: RND

    INTEGER :: i,j,k,ia,is,iat

!-------------------------i L_u dt Operator-----------------------------!

    DO iat=1,ions1%nat
       ia=iatpt(1,iat)
       is=iatpt(2,iat)
       DO k=1,3

!            CALL random_force(RND)
            taup(k,ia,is) = tau0(k,ia,is) + velp(k,ia,is) * dt_ions
            DO j=1,l_sinr
                CALL random_force(RND)
                v2_sinr(j,k,ia,is)=v2_sinr(j,k,ia,is)*emgt+sigma_sinr*RND*sqe2gt
            ENDDO
        ENDDO
    ENDDO
  END SUBROUTINE posupi_sinr

!!.....................................................................................!!

  SUBROUTINE random_force(RND)
    REAL(real_8) :: RND
    REAL(real_8) :: ran(1),alfa1

! Create Gaussian Random Force
 
    CALL repprngu_vec(1,ran)
    alfa1=2.0_real_8*pi*ran(1)

    RND = DSQRT(-(2.D0*DLOG(repprngu())))*DCOS(alfa1)

  END SUBROUTINE random_force

!!.....................................................................................!!

  SUBROUTINE sinr_cons(esinr)
    REAL(real_8),INTENT(OUT) :: esinr
    INTEGER :: i,ia,is,k,j


    esinr=0._real_8
    DO is=1,ions1%nsp
      DO ia=1,ions0%na(is)
      DO k=1,3
        DO j=1,l_sinr
          esinr=esinr+sinr_mass_1*v1_sinr(j,k,ia,is)**2
        ENDDO
      ENDDO
    ENDDO
    ENDDO

    esinr=esinr*lbylp1!/REAL(nedof,KIND=real_8)


  END SUBROUTINE sinr_cons

!!....................................................................................!!

  SUBROUTINE print_sinr(nfi,velp,v1_sinr,v2_sinr,sinr_ke)

    INTEGER :: nfi
    REAL(real_8) :: v1_sinr(:,:,:,:), v2_sinr(:,:,:,:),sinr_ke, velp(:,:,:)
    INTEGER :: i,ia,is,k,j

    
    OPEN(unit=500,file='RESTART_SINR',status='unknown')

    CALL check_sinr(velp,v1_sinr,l_sinr,sinr_ke)

    DO is=1,ions1%nsp
      DO ia=1,ions0%na(is)
        DO k=1,3
          DO j=1,l_sinr
            WRITE(500,'(I10,2X,3E16.6)')&
                 NFI,sinr_ke,v1_sinr(j,k,ia,is),v2_sinr(j,k,ia,is)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

   !CALL fileclose(500)


   END SUBROUTINE print_sinr     

!!...................................................................................!!

   SUBROUTINE read_restart_sinr(v1_sinr,v2_sinr,sinr_ke)

    INTEGER :: nfi
    REAL(real_8) :: v1_sinr(:,:,:,:), v2_sinr(:,:,:,:),sinr_ke
    INTEGER :: i,ia,is,k,j
    
    DO is=1,ions1%nsp
      DO ia=1,ions0%na(is)
        DO k=1,3
          DO j=1,l_sinr

            READ(601,*)&
                NFI,sinr_ke,v1_sinr(j,k,ia,is),v2_sinr(j,k,ia,is)

          ENDDO
        ENDDO
      ENDDO
    ENDDO
    END SUBROUTINE read_restart_sinr

!!-----------------------------------------------------------------------------------------!!
   
    SUBROUTINE check_sinr(velp,v1_sinr,l_sinr,sinr_ke)
     
        REAL(real_8), INTENT(IN) :: velp(:,:,:),v1_sinr(:,:,:,:)
        INTEGER, INTENT(IN) :: l_sinr
        INTEGER :: i,ia,is,k,j
        REAL(real_8) :: sinr_ke, ke

        sinr_ke=0._real_8

        DO is=1,ions1%nsp
           DO ia=1,ions0%na(is)
              DO k=1,3
                 ke=rmass%pma(is)*velp(k,ia,is)**2._real_8
                 DO j=1,l_sinr
                    ke=ke+lbylp1*sinr_mass_1*v1_sinr(j,k,ia,is)**2
                 ENDDO
                 
                 sinr_ke=sinr_ke+ke
              ENDDO
           ENDDO
        ENDDO
        sinr_ke=sinr_ke/(3._real_8*REAL(ions1%nat,kind=real_8))

    END SUBROUTINE check_sinr    

END MODULE sinr_utils
