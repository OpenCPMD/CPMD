#include "cpmd_global.h"

MODULE distribution_utils
  USE error_handling,                  ONLY: stopgm
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE kinds,                           ONLY: real_8
  USE nlps,                            ONLY: nlps_com
  USE pslo,                            ONLY: pslo_com
  USE system,                          ONLY: iatpt

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dist_entity
  PUBLIC :: dist_entity2
  PUBLIC :: dist_atoms
CONTAINS

    SUBROUTINE dist_entity(entity,num,dist_array,nblock,nbmax,nblocal,iloc,fw)
    ! ==--------------------------------------------------------------==
    ! == Distributes entity into num equally parts
    INTEGER,INTENT(IN)                       :: entity,num
    INTEGER,INTENT(OUT) __CONTIGUOUS         :: dist_array(:,:)
    INTEGER,INTENT(IN),OPTIONAL              :: nblock, iloc
    LOGICAL,INTENT(IN),OPTIONAL              :: fw 
    INTEGER,INTENT(OUT),OPTIONAL             :: nbmax, nblocal
    INTEGER                                  :: ip, n
    REAL(real_8)                             :: xaim, xentity
    CHARACTER(*), PARAMETER                  :: procedureN = 'distribution'
    LOGICAL                                  :: fw_local
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, FAU Erlangen-Nuernberg March 2019

    IF (PRESENT(nblock)) THEN
       IF (nblock*num .GE. entity) THEN
          xentity=REAL(nblock,kind=real_8)
       ELSE
          xentity=REAL(entity,kind=real_8)
       END IF
    ELSE
       xentity=REAL(entity,kind=real_8)
    END IF
    dist_array(1,:)=0
    dist_array(2,:)=-1

    fw_local=.false.

    IF (PRESENT(FW) ) fw_local=FW

    IF (fw_local) THEN
       DO ip=1,num
          xaim = ip * xentity/num
          dist_array(2,ip)=NINT(xaim)
          IF (ip.EQ.1) THEN
             dist_array(1,ip)=1
          !failsave
          ELSE
             dist_array(1,ip)=dist_array(2,ip-1)+1
          END IF
          IF (NINT(xaim).GT.entity) dist_array(2,ip)=entity
          IF (ip.EQ.num) dist_array(2,ip)=entity
          IF (dist_array(1,ip).GT.entity) dist_array(1,ip)=entity+1
       ENDDO
    ELSE
       DO ip=num,1,-1
          xaim = xentity-(num-ip+1)*xentity/num+1
          dist_array(1,ip)=NINT(xaim)
          IF (ip.EQ.num) THEN
             dist_array(2,ip)=entity
          !failsave
          ELSE
             dist_array(2,ip)=dist_array(1,ip+1)-1
          END IF
          IF (ip.EQ.1) THEN
             dist_array(1,ip)=1
          ENDIF
       ENDDO
    END IF

    IF (PRESENT(nbmax)) THEN
       nbmax=0
       DO ip=0,num-1
          nbmax=MAX(nbmax,dist_array(2,ip+1)-dist_array(1,ip+1)+1)
       END DO
    END IF

    IF (PRESENT(nblocal)) THEN
       IF (PRESENT(iloc)) THEN
          IF (iloc.GE.num) CALL stopgm(procedureN,'Programming error',&
               __LINE__,__FILE__)
          nblocal=dist_array(2,iloc+1)-dist_array(1,iloc+1)+1
       ELSE
          CALL stopgm(procedureN,'Programming error',&
               __LINE__,__FILE__)
       END IF
    END IF

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dist_entity
  ! ==================================================================
  SUBROUTINE dist_entity2(entity,num,dist_array,nblock,nbmax,nblocal,iloc,fw)
    ! ==--------------------------------------------------------------==
    ! == Distributes entity into num equally parts
    INTEGER,INTENT(IN)                       :: entity,num
    INTEGER,INTENT(OUT) __CONTIGUOUS         :: dist_array(:,:)
    INTEGER,INTENT(IN),OPTIONAL              :: nblock, iloc
    LOGICAL,INTENT(IN),OPTIONAL              :: fw 
    INTEGER,INTENT(OUT),OPTIONAL             :: nbmax, nblocal
    INTEGER                                  :: ip, n
    REAL(real_8)                             :: xaim, xentity
    CHARACTER(*), PARAMETER                  :: procedureN = 'distribution'
    LOGICAL                                  :: fw_local
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, FAU Erlangen-Nuernberg March 2019

    IF (PRESENT(nblock)) THEN
       IF (nblock*num .GE. entity) THEN
          xentity=REAL(nblock,kind=real_8)
       ELSE
          xentity=REAL(entity,kind=real_8)
       END IF
    ELSE
       xentity=REAL(entity,kind=real_8)
    END IF

    dist_array(:,1)=0
    dist_array(:,2)=-1

    fw_local=.FALSE.

    IF (PRESENT(FW)) fw_local=FW

    IF (fw_local) THEN
       DO ip=1,num
          xaim = ip * xentity/num
          dist_array(ip,2)=NINT(xaim)
          IF (ip.EQ.1) THEN
             dist_array(ip,1)=1
          !failsave
          ELSE
             dist_array(ip,1)=dist_array(ip-1,2)+1
          END IF
          IF (NINT(xaim).GT.entity) dist_array(ip,2)=entity
          IF (ip.EQ.num) dist_array(ip,2)=entity
          IF (dist_array(ip,1).GT.entity) dist_array(ip,1)=entity+1
       ENDDO
    ELSE
       DO ip=num,1,-1
          xaim = xentity-(num-ip+1)*xentity/num+1
          dist_array(ip,1)=NINT(xaim)
          IF (ip.EQ.num) THEN
             dist_array(ip,2)=entity
          !failsave
          ELSE
             dist_array(ip,2)=dist_array(ip+1,1)-1
          END IF
          IF (ip.EQ.1) THEN
             dist_array(ip,1)=1
          ENDIF
       ENDDO
    END IF

    IF (PRESENT(nbmax)) THEN
       nbmax=0
       DO ip=0,num-1
          nbmax=MAX(nbmax,dist_array(ip+1,2)-dist_array(ip+1,1)+1)
       END DO
    END IF

    IF (PRESENT(nblocal)) THEN
       IF (PRESENT(iloc)) THEN
          IF (iloc.GE.num) CALL stopgm(procedureN,'Programming error',&
               __LINE__,__FILE__)
          nblocal=dist_array(iloc+1,2)-dist_array(iloc+1,1)+1
       ELSE
          CALL stopgm(procedureN,'Programming error',&
               __LINE__,__FILE__)
       END IF
    END IF

    ! ==--------------------------------------------------------------==
    RETURN
  END SUBROUTINE dist_entity2
  ! ==================================================================
  SUBROUTINE  dist_atoms(buffer_count,fractions,na,na_split,map_work,only_uspp)
    ! ==--------------------------------------------------------------==
    ! == Splits atoms into buffer_count chunks according to fractions ==
    ! == Usefull for distribution of work related to beta projectors  ==
    ! == Requieres a custom na array, containing the first and last IA==
    ! == for each species                                             ==
    ! ==--------------------------------------------------------------==
    ! Author: Tobias Kloeffel, FAU Erlangen-Nuernberg March 2019

    INTEGER, INTENT(IN)                      :: buffer_count
    REAL(real_8),INTENT(IN)                  :: fractions(buffer_count)
    INTEGER, INTENT(IN)                      :: na(2,ions1%nsp)
    INTEGER, INTENT(OUT)                     :: na_split(2,ions1%nsp,buffer_count)
    INTEGER, INTENT(OUT), OPTIONAL           :: map_work(ions1%nat+1,buffer_count)
    LOGICAL, INTENT(IN), OPTIONAL            :: only_uspp
    INTEGER                                  :: tot_work, &
                                                div(buffer_count),&
                                                is, iv, ia, isa, isa0, buff, aim, count, &
                                                map(ions1%nat+1,buffer_count)
    LOGICAL                                  :: uspp

    IF(PRESENT(only_uspp))THEN
       uspp=only_uspp
    ELSE
       uspp=.FALSE.
    END IF
    tot_work=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is).OR..NOT.uspp) THEN
          tot_work=tot_work+(na(2,is)-na(1,is)+1)*nlps_com%ngh(is)
       END IF
    END DO

    div=0
    DO buff=1,buffer_count
       div(buff)=CEILING(REAL(tot_work,real_8)*fractions(buff))
    END DO

    map=0
    count=0
    buff=1
    aim=div(buff)
    isa0=0
    DO is=1,ions1%nsp
       IF (pslo_com%tvan(is).OR..NOT.uspp) THEN
          DO ia=na(1,is),na(2,is)
             isa=isa0+ia
             !fail-save last buffer takes any remaining atom
             IF (buff.EQ.buffer_count) THEN
                map(isa,buff)=1
                count=count+nlps_com%ngh(is)
                !accumulate atoms in this buffer
             ELSE IF(count+nlps_com%ngh(is).LE.aim) THEN
                map(isa,buff)=1
                count=count+nlps_com%ngh(is)
                !the current atom would not fit into this buffer, add it to the next buffer
             ELSE IF( count+nlps_com%ngh(is).GT.aim) THEN
                buff=buff+1
                aim=aim+div(buff)
                map(isa,buff)=1
                count=count+nlps_com%ngh(is)
             ENDIF
          ENDDO
       END IF
       !isa0/isa is always the 'real' isa
       isa0=isa0+ions0%na(is)
    ENDDO

    na_split(1,:,:)=0
    na_split(2,:,:)=-1

    DO buff=1,buffer_count
       isa0=0
       DO is=1,ions1%nsp
          DO ia=na(1,is),na(2,is)
             isa=isa0+ia
             !check if this atom is included in this buffer
             IF (map(isa,buff).EQ.1) then
                !Check if we found the first atom of this species
                IF (na_split(1,is,buff).EQ.0) na_split(1,is,buff)=ia
                !Set the last atom of this species
                na_split(2,is,buff)=ia
             ENDIF
          ENDDO
          isa0=isa0+ions0%na(is)
       ENDDO
    ENDDO

    IF (PRESENT(map_work)) THEN
       DO buff=1,buffer_count
          count=0
          DO isa=1,ions1%nat
             IF (map(isa,buff).EQ.1) THEN
                is=iatpt(2,isa)
                count=count+nlps_com%ngh(is)
             END IF
          END DO
          map(ions1%nat+1,buff)=count
       END DO
       map_work=map
    END IF
  END SUBROUTINE dist_atoms

END MODULE distribution_utils
