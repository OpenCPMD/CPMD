#include "scratch_defines.h"

MODULE segment_managment_utils
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  USE data_managment_utils
  IMPLICIT NONE
  PRIVATE

  TYPE short_index
     INTEGER( INT64 ) :: index_i( 6 )
     LOGICAL          :: index_l( 2 )
     CHARACTER( 256 ) :: index_c
  END type short_index


  TYPE, PUBLIC :: segment
     TYPE(scr), PRIVATE                      :: data
     ! all segments are padded to begin at an alignment boundary and have two guard blocks before and after the payload. Since the segments always end at a aligment boundary, the start guard block is exactly one alignment block, whereas the end guard block is at least one alignment block and max one byte less then two alignment blocks in length, this way, every segment is a multiple of aligment blocks in length
     INTEGER( INT64 ), ALLOCATABLE, PRIVATE  :: index_i( : , : ) !< integer index array holding the start, end, length of the actual segments and the padded segments, size(index_i,1)=number of segments, size(index_i,2)=6 => requested(start,end,len),padded(start,end,len), for free segments
     LOGICAL, ALLOCATABLE, PRIVATE           :: index_l( : , : ) !< logical index array for each segment, size(index_l,1)=number of segments, size(index_l,2)=2 => free, moveable
     CHARACTER(256),ALLOCATABLE,PRIVATE     :: index_c( : ) !< character array holding user supplied string for each segment, size(tag,1)= number of segments
  END TYPE segment

  !Segment_managment is the interface and handles all requestes from pool_managment and creates
  !requestes to data_managment
  PUBLIC :: init_segment
  PUBLIC :: finalize_segment
  PUBLIC :: is_free
  PUBLIC :: request_saved_segment
  PUBLIC :: request_free_segment
  PUBLIC :: free_segment
  PUBLIC :: save_segment
  
CONTAINS
  
  SUBROUTINE init_segment( this, requested_len, ierr )
    !just allocate this segment
    TYPE( segment ), INTENT( OUT )  :: this
    INTEGER( INT64 ), INTENT( IN )  :: requested_len
    INTEGER( INT32 ), INTENT( OUT ) :: ierr

    INTEGER( INT32 )                :: id

#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(6I4)' ) selected_int_kind( 15 ), INT64
#endif
    
    id = 1
    CALL allocate_index( this, id, ierr)
    IF( ierr /= 0 )EXIT_ON_ERROR

    !update the index
    this%index_i( id, 1 ) = pad_in_char + ONE_INT64
    this%index_i( id, 2 ) = requested_len + pad_in_char
    this%index_i( id, 3 ) = requested_len
    this%index_i( id, 4 ) = ONE_INT64
    this%index_i( id, 5 ) = requested_len + pad_in_char + alignment_in_char + alignment_in_char &
         - MOD( requested_len, alignment_in_char )
    this%index_i( id, 6 ) = this%index_i( id, 5 ) - this%index_i ( id, 4 ) + ONE_INT64
    
    this%index_l( id, 1 ) = .TRUE.
    this%index_l( id, 2 ) = .TRUE.

    this%index_c( id ) = ''
    
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A)' ) "Initialized new node with segment:"
    WRITE( OUTPUT_UNIT, '(6I17)' ) this%index_i
#endif
!    ALLOCATE( this%data)
    CALL allocate_array( this%data, requested_len, ierr )

  END SUBROUTINE init_segment
  
  SUBROUTINE finalize_segment( this, len, ierr )
    !just allocate this segment
    TYPE( segment ), INTENT( INOUT ) :: this
    INTEGER( INT64 ), INTENT( OUT ) :: len
    INTEGER( INT32 ), INTENT( OUT ) :: ierr

#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A)') "Warning, deallocating segment"
    WRITE( OUTPUT_UNIT, '(A)') this%index_c
#endif
    
    CALL deallocate_array( this%data, len, ierr )
    IF( ierr /= 0 )EXIT_ON_ERROR

    CALL deallocate_index( this, ierr)
    IF( ierr /= 0 )EXIT_ON_ERROR
    
  END SUBROUTINE finalize_segment

  SUBROUTINE is_free( this, free )
    !check only if this segment is free
    TYPE( segment ), INTENT( INOUT ) :: this
    LOGICAL, INTENT( OUT )           :: free

    INTEGER( INT32 )                 :: id
    INTEGER( INT32 )                 :: num_segments
    
    num_segments = SIZE( this%index_l, 1 )
    DO id = 1, num_segments
       IF( .NOT. this%index_l( id, 1 ) ) EXIT
    END DO

    IF( id <= num_segments )THEN
       free = .FALSE.
    ELSE
       free = .TRUE.
    END IF
    
  END SUBROUTINE is_free
  
  SUBROUTINE request_saved_segment( user_ptr, this, requested_len, tag, ierr )
    ! return the id of a free segment with length greater or equal len
    CHARACTER, POINTER, CONTIGUOUS, INTENT( OUT ) :: user_ptr( : )
    TYPE( segment ), INTENT( INOUT )              :: this
    INTEGER( INT64 ), INTENT( IN )                :: requested_len
    CHARACTER( 256 ), INTENT( IN )                :: tag
    INTEGER( INT32 ), INTENT( OUT )               :: ierr
    
    INTEGER( INT32 )                              :: id
    INTEGER( INT64 )                              :: pos( 3 )
    
    CALL search_tagged_segment( this, .TRUE., requested_len, tag, id, ierr )

    IF( ierr < 0 )RETURN
    IF( ierr /= 0 )EXIT_ON_ERROR

    !get payload position
    pos( 1 ) = this%index_i( id, 1 )
    pos( 2 ) = this%index_i( id, 2 )
    pos( 3 ) = this%index_i( id, 3 )

    !mark this segment as locked
    this%index_l( id, 2 ) = .FALSE.

    CALL associate_ptr( user_ptr, this%data, pos, ierr )
    
  END SUBROUTINE request_saved_segment

  SUBROUTINE request_free_segment( user_ptr, this, requested_len, tag, ierr )
    ! return the id of a free segment with length greater or equal len
    CHARACTER, POINTER, CONTIGUOUS, INTENT( OUT ) :: user_ptr( : )
    TYPE( segment ), INTENT( INOUT )              :: this
    INTEGER( INT64 ), INTENT( IN )                :: requested_len
    CHARACTER( 256 ), INTENT( IN )                :: tag
    INTEGER( INT32 ), INTENT( OUT )               :: ierr
    
    INTEGER( INT32 )                              :: id
    INTEGER( INT64 )                              :: pos( 3 )
    
    CALL search_free_segment( this, requested_len, id, ierr )

    IF( ierr <  0 )RETURN
    IF( ierr /= 0 )EXIT_ON_ERROR

    !mark this segment as locked and not free
    this%index_l( id, 1 ) = .FALSE.
    this%index_l( id, 2 ) = .FALSE.
    !set the tag
    this%index_c( id ) = tag

    pos( 1 ) = this%index_i( id, 1 )
    pos( 2 ) = this%index_i( id, 2 )
    pos( 3 ) = this%index_i( id, 3 )
    
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,2I17)' ) "requested len, segment len:", requested_len, this%index_i( id, 3 )
    WRITE( OUTPUT_UNIT, '(3I17)' ) pos
#endif
    !associate the user ptr to the data ptr
    CALL associate_ptr( user_ptr, this%data, pos, ierr )
    IF( ierr /= 0)EXIT_ON_ERROR

  END SUBROUTINE request_free_segment

  SUBROUTINE free_segment( user_ptr, this, requested_len, tag, ierr )
    ! return the id of a free segment with length greater or equal len
    CHARACTER, POINTER, CONTIGUOUS, INTENT( IN ) :: user_ptr( : )
    TYPE( segment ), INTENT( INOUT )             :: this
    INTEGER( INT64 ), INTENT( IN )               :: requested_len
    CHARACTER( 256 ), INTENT( IN )               :: tag
    INTEGER( INT32 ), INTENT( OUT )              :: ierr
    
    INTEGER( INT32 )                             :: id
    INTEGER( INT64 )                             :: pos( 3 )
    
    CALL search_tagged_segment( this, .FALSE., requested_len, tag, id, ierr )
    IF( ierr < 0 )RETURN
    IF( ierr /= 0 )EXIT_ON_ERROR

    pos( 1 ) = this%index_i( id, 1 )
    pos( 2 ) = this%index_i( id, 2 )
    pos( 3 ) = this%index_i( id, 3 )
    
    CALL check_associate_ptr( user_ptr, this%data, pos, ierr )
    IF( ierr /= 0 )EXIT_ON_ERROR

    this%index_l( id, 1 ) = .TRUE.
    this%index_l( id, 2 ) = .TRUE.
    this%index_c( id ) = ''
    
    CALL cleanup_segments( this, id )
    
  END SUBROUTINE free_segment
  
  SUBROUTINE save_segment( user_ptr, this, requested_len, tag, ierr )
    ! return the id of a free segment with length greater or equal len
    CHARACTER, POINTER, CONTIGUOUS, INTENT( IN ) :: user_ptr( : )
    TYPE( segment ), INTENT( INOUT )             :: this
    INTEGER( INT64 ), INTENT( IN )               :: requested_len
    CHARACTER( 256 ), INTENT( IN )               :: tag
    INTEGER( INT32 ), INTENT( OUT )              :: ierr

    INTEGER( INT32 )                             :: id
    INTEGER( INT64 )                             :: pos( 3 )

    CALL search_tagged_segment( this, .FALSE., requested_len, tag, id, ierr )
    IF( ierr < 0 )RETURN
    IF( ierr /= 0 )EXIT_ON_ERROR
    
    pos( 1 ) = this%index_i( id, 1 )
    pos( 2 ) = this%index_i( id, 2 )
    pos( 3 ) = this%index_i( id, 3 )

    CALL check_associate_ptr( user_ptr, this%data, pos, ierr )
    IF( ierr /= 0 )EXIT_ON_ERROR
    
    !mark segment as moveable
    this%index_l( id, 2 ) = .TRUE.
    
  END SUBROUTINE save_segment

  SUBROUTINE search_free_segment( this, requested_len, id, ierr )
    ! return the id of a free segment with length greater or equal len
    TYPE( segment ), INTENT( INOUT ) :: this
    INTEGER( INT64 ), INTENT( IN )   :: requested_len
    INTEGER( INT32 ), INTENT( OUT )  :: id
    INTEGER( INT32 ), INTENT( OUT )  :: ierr

    ierr = 0
    
    DO id = 1, SIZE( this%index_i, 1 )
       !get a free segment
       IF( this%index_l( id, 1 ) ) THEN
          !check if payload len is lower or equal len
          IF( this%index_i( id, 3 ) >= requested_len ) THEN
             EXIT
          END IF
       END IF
    END DO

    IF( id > SIZE( this%index_i, 1 ) ) THEN
#if defined(_DEBUG)
       WRITE(OUTPUT_UNIT, '(A)' ) "No free segment found that might be big enough"
#endif
       id = -1
       ierr = -1
       RETURN
    END IF

    IF( this%index_i( id, 3 ) > requested_len ) THEN

#if defined(_DEBUG)
       WRITE( OUTPUT_UNIT, '(A,2I17,I4)') "found larger segment than needed, found len, &
            requested_len, in segment id:", this%index_i( id, 3 ), requested_len, id
#endif
       CALL shrink_segment( this, id, requested_len, ierr )

       IF( ierr /= 0 )EXIT_ON_ERROR
    END IF

  END SUBROUTINE search_free_segment
  
  SUBROUTINE search_tagged_segment( this, moveable, requested_len, tag, id, ierr )
    ! return the id of the segment specified by len, tag, and moveable
    TYPE( segment ), INTENT( INOUT ) :: this
    LOGICAL, INTENT( IN )            :: moveable
    INTEGER( INT64 ), INTENT( IN )   :: requested_len
    CHARACTER( 256 ), INTENT( IN )   :: tag
    INTEGER( INT32 ), INTENT( OUT )  :: id
    INTEGER( INT32 ), INTENT( OUT )  :: ierr
#if defined(_DEBUG)
    INTEGER( INT32 )                 :: i
#endif

    DO id = 1, SIZE( this%index_i, 1 )
       !get a non free segment
       IF( .NOT. this%index_l( id, 1 ) ) THEN
          !check if request len is matching
          IF( this%index_i( id, 3 ) == requested_len ) THEN
             !check if this segment is moveable or locked
             IF( this%index_l( id, 2) .EQV. moveable ) THEN
                !check if the tag is matching
                IF( this%index_c( id ) == tag ) THEN
                   ierr = 0
                   EXIT
                END IF
             END IF
          END IF
       END IF
    END DO

    IF( id > SIZE( this%index_i, 1 ) ) THEN
       id = -1
       ierr = -1

#if defined(_DEBUG)
       WRITE( OUTPUT_UNIT, '(A,L2,I17,X,A)' ) "No occupied segment found", &
            moveable, requested_len, ADJUSTL(tag)
       DO i = 1, SIZE( this%index_i, 1 )
          WRITE( OUTPUT_UNIT, '(6I17)') this%index_i( i, : )
          WRITE( OUTPUT_UNIT, '(2L2)') this%index_l( i, : )
          WRITE( OUTPUT_UNIT, '(A)') this%index_c( i )
       END DO
#endif

       RETURN
    END IF
    
  END SUBROUTINE search_tagged_segment

  SUBROUTINE shrink_segment( this, new_id, requested_len, ierr )
    !shrink a segment to requested_len
    !first check if can shrink in place
    !if not possible:
    !a)reuse an existing empty segment
    !b)grow the index by one
    !note that on freeing a segment the freed segment will be combined with any other free segment and
    !thus we are sure that there is no free segment close to any other free segment
    TYPE( segment ), INTENT( INOUT ) :: this
    INTEGER( INT32 ), INTENT( IN )   :: new_id
    INTEGER( INT64 ), INTENT( IN )   :: requested_len
    INTEGER( INT32 ), INTENT( OUT )  :: ierr

    INTEGER( INT64 )                 :: old_len
    INTEGER( INT32 )                 :: id
    INTEGER( INT32 )                 :: old_size
    
    old_len = this%index_i( new_id, 3 )
    ierr = 0
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A)' ) "request to add a segment, old segment data"
    old_size = SIZE( this%index_i, 1 )
    WRITE( OUTPUT_UNIT, '(I4)' ) old_size
    DO id = 1, old_size
       WRITE( OUTPUT_UNIT, '(6I17)' ) this%index_i( id, : )
       WRITE( OUTPUT_UNIT, '(2L2)' ) this%index_l( id, : )
       WRITE( OUTPUT_UNIT, '(A)' ) this%index_c( id )
    END DO
#endif

    IF( ( old_len - requested_len ) <= alignment_in_char )THEN
       ! old_len is always a multiple of alignment_in_char and there are always 2 alignment blocks

#if defined(_DEBUG)
       WRITE( OUTPUT_UNIT, '(A)' ) "do not add a new segment, just grow this segment"
#endif

       this%index_i( new_id, 2 ) = this%index_i( new_id, 2 ) - old_len + requested_len
       this%index_i( new_id, 3 ) = requested_len

#if defined(_DEBUG)       
       WRITE( OUTPUT_UNIT, '(A)' ) "request to add a segment, new segment data"
       WRITE( OUTPUT_UNIT, '(6I17)' ) this%index_i( new_id, : )
#endif

       RETURN
    END IF

    old_size = SIZE( this%index_i, 1 )

    IF( this%index_i( old_size, 6 ) == ZERO_INT64 )THEN

#if defined(_DEBUG)
       WRITE( OUTPUT_UNIT, '(A)' ) "reusing empty index"
#endif
       
       !shift all segments by one to the end
       DO id = old_size, new_id + 1 , -1
          this%index_i( id, : ) = this%index_i( id - 1, : )
          this%index_l( id, : ) = this%index_l( id - 1, : )
          this%index_c( id ) = this%index_c( id - 1 )
       END DO
#if defined(_DEBUG)
       WRITE( OUTPUT_UNIT, '(A)' ) "after shift"
       DO id = 1, SIZE( this%index_i, 1 )
          WRITE( OUTPUT_UNIT, '(6I17)' ) this%index_i( id, : )
          WRITE( OUTPUT_UNIT, '(2L2)' ) this%index_l( id, : )
          WRITE( OUTPUT_UNIT, '(A)' ) this%index_c( id )
       END DO
#endif

    ELSE
       CALL grow_index( this, ierr )
       IF( ierr /= 0 )EXIT_ON_ERROR
       !shift all segments by one to the end
       DO id = SIZE( this%index_i, 1 ), new_id + 2 , -1
          this%index_i( id, : ) = this%index_i( id - 1, : )
          this%index_l( id, : ) = this%index_l( id - 1, : )
          this%index_c( id ) = this%index_c( id - 1 )
       END DO
#if defined(_DEBUG)
       WRITE( OUTPUT_UNIT, '(A)' ) "after shift"
       DO id = 1, SIZE( this%index_i, 1 )
          WRITE( OUTPUT_UNIT, '(6I17)' ) this%index_i( id, : )
          WRITE( OUTPUT_UNIT, '(2L2)' ) this%index_l( id, : )
          WRITE( OUTPUT_UNIT, '(A)' ) this%index_c( id )
       END DO
#endif
       
    END IF

    !split segment
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A)' ) "request to add a segment, new segment data 1"
    DO id = 1, SIZE( this%index_i, 1 )
       WRITE( OUTPUT_UNIT, '(6I17)' ) this%index_i( id, : )
       WRITE( OUTPUT_UNIT, '(2L2)' ) this%index_l( id, : )
       WRITE( OUTPUT_UNIT, '(A)' ) this%index_c( id )
    END DO
#endif

    
    !backup the old segment end
    this%index_i( new_id + 1 , 5 ) = this%index_i( new_id, 5 )

    !set the new segment bounds
    this%index_i( new_id, 2 ) = this%index_i( new_id, 1 ) + requested_len - ONE_INT64
    this%index_i( new_id, 3 ) = requested_len
    this%index_i( new_id, 5 ) = this%index_i( new_id, 2 ) + alignment_in_char + alignment_in_char &
         - MOD( requested_len, alignment_in_char )
    this%index_i( new_id, 6 ) = requested_len + pad_in_char + alignment_in_char + alignment_in_char &
         - MOD( requested_len, alignment_in_char )

    !set the remainder segment
    this%index_i( new_id + 1 , 4 ) = this%index_i( new_id, 5 ) + ONE_INT64
    this%index_i( new_id + 1 , 6 ) = this%index_i( new_id + 1, 5 ) - this%index_i( new_id + 1, 4 ) &
         + ONE_INT64
    
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A)' ) "request to add a segment, new segment data 2"
    DO id = 1, SIZE( this%index_i, 1 )
       WRITE( OUTPUT_UNIT, '(6I17)' ) this%index_i( id, : )
       WRITE( OUTPUT_UNIT, '(2L2)' ) this%index_l( id, : )
       WRITE( OUTPUT_UNIT, '(A)' ) this%index_c( id )
    END DO
#endif

    CALL maximize_segment( this, new_id + 1 )
    this%index_l( new_id + 1, 1 ) = .TRUE.
    this%index_l( new_id + 1, 2 ) = .TRUE.
    this%index_c( new_id + 1 ) = ''

#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A)' ) "request to add a segment, new segment data"
    DO id = 1, SIZE( this%index_i, 1 )
       WRITE( OUTPUT_UNIT, '(6I17)' ) this%index_i( id, : )
       WRITE( OUTPUT_UNIT, '(2L2)' ) this%index_l( id, : )
       WRITE( OUTPUT_UNIT, '(A)' ) this%index_c( id )
    END DO
#endif

  END SUBROUTINE shrink_segment
  
  SUBROUTINE grow_index( this, ierr )
    !grow index array by one
    TYPE( segment ), INTENT( INOUT ) :: this
    INTEGER( INT32 ), INTENT( OUT )  :: ierr

    TYPE( segment )                  :: this_backup
    INTEGER( INT32 )                 :: new_size

    new_size = SIZE( this%index_i, 1 ) + 1 !< expand index by 1
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,X2I4)' ) "growing segment to", new_size
#endif
    !backup old index arrays
    CALL backup_index( this, this_backup )

    CALL allocate_index( this, new_size, ierr)
    IF( ierr /= 0 ) EXIT_ON_ERROR

    CALL restore_index( this, this_backup, new_size, ierr )
    IF( ierr /= 0 ) EXIT_ON_ERROR

#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,I4,6L2)') " newsize, index arrays allocated? backup allocated?", &
         SIZE(this%index_i,1 ), ALLOCATED(this%index_i), ALLOCATED(this%index_l), &
         ALLOCATED(this%index_c), ALLOCATED(this_backup%index_i), ALLOCATED(this_backup%index_l), &
         ALLOCATED(this_backup%index_c) 
#endif

  END SUBROUTINE grow_index
  
  SUBROUTINE cleanup_segments( this, start_id )
    TYPE( segment ), INTENT( INOUT ) :: this
    INTEGER( INT32 ), INTENT( IN)    :: start_id

    INTEGER( INT32 )                 :: temp_id
    INTEGER( INT32 )                 :: id_left
    INTEGER( INT32 )                 :: id_right
    INTEGER( INT32 )                 :: num_segments
#if defined(_DEBUG)
    INTEGER( INT32 )                 :: i
#endif

    CALL maximize_segment( this, start_id )

    num_segments = SIZE( this%index_i, 1)
    IF( num_segments == 1 ) RETURN

    id_left = start_id - 1
    id_right = start_id + 1
    temp_id = start_id
    
    !check left and right id for merging
    IF( id_left > 0 )THEN
       IF( this%index_i( id_left, 6 ) > ZERO_INT64 )THEN
          IF( this%index_l( id_left, 1 ) )THEN
             CALL combine_segments( this, id_left, temp_id )
             id_right = id_right - 1
             temp_id = temp_id - 1
          END IF
       END IF
    END IF
    IF( id_right <= num_segments )THEN
       IF( this%index_i( id_right, 6 ) > ZERO_INT64 )THEN
          IF( this%index_l( id_right, 1 ) )THEN
             CALL combine_segments( this, temp_id, id_right )
          END IF
       END IF
    END IF

    !at the moment it is not supported to merge moveable segments
    !this would potentialle cause a lot of overhead and is only reasonable for
    !tight memory constraints

#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A)' ) "cleaning"
    DO i = 1, SIZE( this%index_i, 1 )
       WRITE( OUTPUT_UNIT, '(6I17)') this%index_i( i, : )
       WRITE( OUTPUT_UNIT, '(2L2)') this%index_l( i, : )
       WRITE( OUTPUT_UNIT, '(A)') this%index_c( i )
    END DO
    WRITE( OUTPUT_UNIT, '(A)' ) "end cleaning"
#endif
    
  END SUBROUTINE cleanup_segments

  SUBROUTINE combine_segments( this, left_id, right_id )
    !combines segments id1 and id2, empty segment id2 will be moved to the end of the index array
    TYPE( segment ), INTENT( INOUT ) :: this
    INTEGER( INT32 ), INTENT( IN )   :: left_id
    INTEGER( INT32 ), INTENT( IN )   :: right_id

    INTEGER( INT32 )                 :: id
    INTEGER( INT32 )                 :: segment_size
    
#if defined(_DEBUG)
    IF( right_id - 1 /= left_id )THEN
       WRITE( OUTPUT_UNIT ,'(A1X1I4X1I4)') "Cannot combine nonconsecutive segment ids", &
            left_id, right_id
       EXIT_ON_ERROR
    END IF
#endif
    !segment start remains
    !end position becomes end of right id
    this%index_i( left_id, 5 ) = this%index_i( right_id, 5 )
    this%index_i( left_id, 6 ) = this%index_i( left_id, 5 ) -  this%index_i( left_id, 4 ) + ONE_INT64

    CALL maximize_segment( this, left_id )
    
    segment_size = SIZE( this%index_i, 1 )
    IF( right_id < segment_size )THEN
       DO id = right_id, segment_size - 1 
          this%index_i( id, : ) = this%index_i( id + 1, : )
          this%index_l( id, : ) = this%index_l( id + 1, : )
          this%index_c( id ) = this%index_c( id + 1 )
       END DO
    END IF
    !the freed segment is now placed as the last segment
    this%index_i( segment_size, 1 ) = ZERO_INT64
    this%index_i( segment_size, 2 ) = ZERO_INT64
    this%index_i( segment_size, 3 ) = ZERO_INT64
    this%index_i( segment_size, 4 ) = ZERO_INT64
    this%index_i( segment_size, 5 ) = ZERO_INT64
    this%index_i( segment_size, 6 ) = ZERO_INT64
    this%index_l( segment_size, 1 ) = .TRUE.
    this%index_l( segment_size, 2 ) = .TRUE.
    this%index_c( segment_size ) = ''
    
  END SUBROUTINE combine_segments

  SUBROUTINE maximize_segment( this, id  )
    !maximize this segment
    !payload len may be negative
    TYPE( segment ), INTENT( INOUT ) :: this
    INTEGER( INT32 ), INTENT( IN )   :: id

    this%index_i( id, 1 ) = this%index_i( id, 4 ) + pad_in_char
    this%index_i( id, 2 ) = this%index_i( id, 5 ) - ONE_INT64 - alignment_in_char
    this%index_i( id, 3 ) = this%index_i( id, 2 ) - this%index_i( id, 1 ) + ONE_INT64

  END SUBROUTINE maximize_segment
  
  SUBROUTINE remove_segment( this, remove_id, ierr )
    ! remove a segment from the indexes
    TYPE( segment ), INTENT( INOUT ) :: this
    INTEGER( INT32 ), INTENT( IN )   :: remove_id
    INTEGER( INT32 ), INTENT( OUT )  :: ierr
    
    TYPE( segment )                  :: this_backup
    INTEGER( INT32 )                 :: old_size
    INTEGER( INT32 )                 :: new_size

    IF( .NOT. ALLOCATED( this%index_i ) ) THEN
       ierr = -1
       WRITE( OUTPUT_UNIT, '(A)') "Cannot remove unallocated segments!"
       EXIT_ON_ERROR
    END IF

    old_size = SIZE( this%index_i, 1 )

    IF( old_size == 1 ) THEN
       ierr = -1
       WRITE( OUTPUT_UNIT, '(A)') "Cannot remove last segments!"
       EXIT_ON_ERROR
    END IF
    
    new_size = old_size - 1 !< decrease index by 1
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT ,'(A)') 'removing segment',this%index_c
#endif
    IF( old_size /= 0 ) THEN !< get back old indexes
       CALL backup_index( this, this_backup )
    END IF

    CALL allocate_index( this, new_size, ierr )

    IF( old_size /= 0 ) THEN !< get back old indexes
       CALL restore_index( this, this_backup, remove_id, ierr)
    END IF

  END SUBROUTINE remove_segment
  
  SUBROUTINE restore_index( this, this_backup, skip_id, ierr )
    TYPE( segment ), INTENT( INOUT ) :: this
    TYPE( segment ), INTENT( INOUT ) :: this_backup
    INTEGER( INT32 ), INTENT( IN )   :: skip_id
    INTEGER( INT32 ), INTENT( OUT )  :: ierr

    INTEGER( INT32 )                 :: new_size
    INTEGER( INT32 )                 :: old_size
    INTEGER( INT32 )                 :: offset
    INTEGER( INT32 )                 :: i
    INTEGER( INT32 )                 :: j

    old_size = SIZE( this_backup%index_i, 1 )
    new_size = SIZE( this%index_i, 1 )
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,2I4)') "old_size, new_size:", old_size, new_size
#endif
    IF( new_size > old_size ) THEN
       offset = +1
    ELSE
       offset = -1
    END IF
    
    DO j=1, 6
       DO i = 1, skip_id - 1
          this%index_i( i, j ) = this_backup%index_i( i, j )
       END DO
    END DO
    DO j=1, 6
       DO i = skip_id + 1, old_size
          this%index_i( i + offset, j) = this_backup%index_i( i, j )
       END DO
    END DO

    DO j=1, 2
       DO i = 1, skip_id - 1
          this%index_l( i, j ) = this_backup%index_l( i, j )
       END DO
    END DO
    DO j=1, 2
       DO i = skip_id + 1, old_size
          this%index_l( i + offset, j) = this_backup%index_l( i, j )
       END DO
    END DO
    DO i = 1, skip_id - 1
       this%index_c( i ) = this_backup%index_c( i )
    END DO

    DO i = 1, skip_id - 1
       this%index_c( i ) = this_backup%index_c( i )
    END DO
    DO i = skip_id + 1, old_size
       this%index_c( i + offset ) = this_backup%index_c( i )
    END DO

    DEALLOCATE( this_backup%index_i, STAT = ierr ) !< free old backup
    IF ( ierr /= 0 ) THEN
       WRITE( OUTPUT_UNIT, '(A)') "Could not initialize segment array"
       EXIT_ON_ERROR
    END IF
    DEALLOCATE( this_backup%index_l, STAT = ierr )
    IF ( ierr /= 0 ) THEN
       WRITE( OUTPUT_UNIT, '(A)') "Could not initialize segment array"
       EXIT_ON_ERROR
    END IF
    DEALLOCATE( this_backup%index_c, STAT = ierr )
    IF ( ierr /= 0 ) THEN
       WRITE( OUTPUT_UNIT, '(A)') "Could not initialize segment array"
       EXIT_ON_ERROR
    END IF
!    call move_alloc( this_backup%data, this%data)
  END SUBROUTINE restore_index
  
  SUBROUTINE allocate_index( this, new_size, ierr)
    TYPE( segment ), INTENT( INOUT ) :: this
    INTEGER( INT32 ), INTENT ( IN )  :: new_size
    INTEGER( INT32 ), INTENT ( OUT ) :: ierr

    ALLOCATE( this%index_i( new_size, 6 ), STAT = ierr)
    IF ( ierr /= 0 ) THEN
       WRITE( OUTPUT_UNIT, '(A)') "Could not initialize segment array"
       EXIT_ON_ERROR
    END IF
    ALLOCATE( this%index_l( new_size, 2 ), STAT = ierr)
    IF ( ierr /= 0 ) THEN
       WRITE( OUTPUT_UNIT, '(A)') "Could not initialize segment array"
       EXIT_ON_ERROR
    END IF
    ALLOCATE( this%index_c( new_size ), STAT = ierr)
    IF ( ierr /= 0 ) THEN
       WRITE( OUTPUT_UNIT, '(A)') "Could not initialize segment array"
       EXIT_ON_ERROR
    END IF
    
  END SUBROUTINE allocate_index
  
  SUBROUTINE deallocate_index( this, ierr)
    TYPE( segment ), INTENT( INOUT ) :: this
    INTEGER( INT32 ), INTENT ( OUT ) :: ierr

    DEALLOCATE( this%index_i, STAT = ierr)
    IF ( ierr /= 0 ) THEN
       WRITE( OUTPUT_UNIT, '(A)') "Could not deallocate segment array"
       EXIT_ON_ERROR
    END IF
    DEALLOCATE( this%index_l, STAT = ierr)
    IF ( ierr /= 0 ) THEN
       WRITE( OUTPUT_UNIT, '(A)') "Could not deallocate segment array"
       EXIT_ON_ERROR
    END IF
    DEALLOCATE( this%index_c, STAT = ierr)
    IF ( ierr /= 0 ) THEN
       WRITE( OUTPUT_UNIT,'(A)') "Could not deallocate segment array"
       EXIT_ON_ERROR
    END IF

  END SUBROUTINE deallocate_index

  SUBROUTINE backup_index( this, this_backup )
    TYPE( segment ), INTENT( INOUT ) :: this
    TYPE( segment ), INTENT( OUT )   :: this_backup

    CALL MOVE_ALLOC( this%index_i, this_backup%index_i )
    CALL MOVE_ALLOC( this%index_l, this_backup%index_l )
    CALL MOVE_ALLOC( this%index_c, this_backup%index_c )
    
  END SUBROUTINE backup_index

END MODULE segment_managment_utils


#if defined(_TEST_SEGMENT)
PROGRAM test_data_managment_utils
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_LOC
  use data_managment_utils
  USE segment_managment_utils
  IMPLICIT NONE

  TYPE( segment ) :: test_segment
  CHARACTER, POINTER, CONTIGUOUS :: test_ptr1( : ) => NULL()
  CHARACTER, POINTER, CONTIGUOUS :: test_ptr2( : ) => NULL()
  CHARACTER, POINTER, CONTIGUOUS :: test_ptr3( : ) => NULL()
  CHARACTER, POINTER, CONTIGUOUS :: test_ptr4( : ) => NULL()
  CHARACTER, POINTER, CONTIGUOUS :: test_ptr5( : ) => NULL()
  
  INTEGER( INT32 ) :: ierr
  INTEGER( INT64 ) :: len_req
  INTEGER( INT64 ) :: len_out
  CHARACTER( 256 ) :: tag1, tag2, tag3, tag4, tag5, tag6
  LOGICAL          :: free
  INTEGER( INT64 ) :: shift
  TYPE( C_PTR )    :: c_addr
  INTEGER( INT64 ) :: faddr

  
  len_req = INT( 1024, KIND = INT64 )
  len_req = len_req * len_req * len_req * INT( 10, KIND = INT64 ) + ONE_INT64

  CALL init_segment( test_segment, len_req, ierr )
  CALL is_free( test_segment, free )
  WRITE( OUTPUT_UNIT, '(A,L2)') "is free?", free
  tag1 = "test_ptr1"
  tag2 = "test_ptr2"
  tag3 = "test_ptr3"
  tag4 = "test_ptr4"
  tag5 = "test_ptr5"
  tag6 = "test_ptr6"
  
  CALL request_free_segment( test_ptr1, test_segment, len_req, tag1, ierr )
  WRITE( OUTPUT_UNIT, '(A,L2)') "test_ptr1 associated?", ASSOCIATED(test_ptr1)
  WRITE( OUTPUT_UNIT, '(A,I17)') "test_ptr1 sizes?", SIZE(test_ptr1, 1, KIND = INT64 )
  c_addr = C_LOC( test_ptr1( 1 ) )
  faddr = TRANSFER( c_addr , faddr )
  shift = ( alignment - MOD( faddr , alignment ) ) / sizeof_char
  WRITE( OUTPUT_UNIT, '(A,I17)' ) "shift needed to aling:", shift
  CALL is_free( test_segment, free )
  WRITE( OUTPUT_UNIT, '(A,L2)') "is free?", free
  CALL save_segment( test_ptr1, test_segment, len_req, tag1, ierr )
  test_ptr1 => NULL()
  CALL is_free( test_segment, free )
  WRITE( OUTPUT_UNIT, '(A,L2)') "is free?", free
  CALL request_saved_segment( test_ptr1, test_segment, len_req, tag1, ierr )
  WRITE( OUTPUT_UNIT, '(A,L2)') "test_ptr1 associated?", ASSOCIATED(test_ptr1)
  WRITE( OUTPUT_UNIT, '(A,I17)') "test_ptr1 sizes?", SIZE(test_ptr1, 1, KIND = INT64 )
  CALL is_free( test_segment, free )
  WRITE( OUTPUT_UNIT, '(A,L2)') "is free?", free
  CALL free_segment( test_ptr1, test_segment, len_req, tag1, ierr )
  test_ptr1 => NULL()
  CALL is_free( test_segment, free )
  WRITE( OUTPUT_UNIT, '(A,L2)') "is free?", free
  CALL request_free_segment( test_ptr1, test_segment, len_req / INT( 3, KIND = INT64 ) , tag1, ierr )
  CALL request_free_segment( test_ptr2, test_segment, len_req / INT( 14, KIND = INT64 ) , tag2, ierr )
  CALL request_free_segment( test_ptr3, test_segment, len_req / INT( 20, KIND = INT64 ) , tag3, ierr )
  CALL request_free_segment( test_ptr4, test_segment, len_req / INT( 15, KIND = INT64 ) , tag4, ierr )
  CALL request_free_segment( test_ptr5, test_segment, len_req / INT( 14, KIND = INT64 ) , tag5, ierr )
  CALL free_segment( test_ptr1, test_segment, len_req / INT( 3, KIND = INT64 )  , tag1, ierr )
  CALL free_segment( test_ptr2, test_segment, len_req / INT( 14, KIND = INT64 ) , tag2, ierr )
  CALL free_segment( test_ptr3, test_segment, len_req / INT( 20, KIND = INT64 ) , tag3, ierr )
  CALL free_segment( test_ptr4, test_segment, len_req / INT( 15, KIND = INT64 ) , tag4, ierr )
  CALL free_segment( test_ptr5, test_segment, len_req / INT( 14, KIND = INT64 ) , tag5, ierr )
  CALL finalize_segment( test_segment, len_out, ierr )
  WRITE( OUTPUT_UNIT, '(X,A,I17,X,A,I17,I4)') "len requested", len_req, "len allocated", len_out, ierr
  
END PROGRAM test_data_managment_utils
#endif
