#include "scratch_defines.h"
MODULE pool_managment_utils
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  USE segment_managment_utils

  IMPLICIT NONE
  PRIVATE
  
  TYPE, PRIVATE :: pool_ptr
     TYPE(segment), PRIVATE, ALLOCATABLE :: ptr !< needs to be an allocatable to allow growing and shrinking of nodes without copying data, if this would be missing, all arrays in scr would have to be pointers or targets and node would have to be a pointer
  END TYPE pool_ptr

  TYPE, PUBLIC :: memory_pool
     TYPE(pool_ptr), PRIVATE, ALLOCATABLE :: node(:) !< each memory pool can have multiple nodes, needs to be an allocatable array to safely distinguish between an initialized and uninitialized memory pool
  END TYPE memory_pool

  TYPE( memory_pool ), PUBLIC :: pool_int !< internal memory pool, only allocated if no external pool is provided
  
  PUBLIC :: reattach_usr_ptr
  PUBLIC :: detach_usr_ptr
  PUBLIC :: alloc_usr_ptr
  PUBLIC :: free_usr_ptr
  PUBLIC :: init_pool
  PUBLIC :: finalize_pool
  
CONTAINS

  SUBROUTINE reattach_usr_ptr( usr_ptr, this, requested_len, tag, ierr )
    !user requested a saved ptr again
    CHARACTER, POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr( : )
    TYPE( memory_pool ), INTENT( INOUT )          :: this
    INTEGER( INT64 ), INTENT( IN )                :: requested_len
    CHARACTER( 256 ), INTENT( IN )                :: tag
    INTEGER( INT32 ), INTENT( OUT )               :: ierr

    INTEGER( INT32)                               :: id

#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,I17,A)') "Trying to find requested user ptr", requested_len, tag
#endif

    DO id = 1, SIZE( this%node, 1 )
       CALL request_saved_segment( usr_ptr, this%node( id )%ptr, requested_len, tag,  ierr )
       !ierr < 0 nothing found, got to next node
       !ierr > 0 some fortran error appeared
       !ierr = 0 user_ptr is associated with some segment
       IF( ierr == 0 ) RETURN
       IF( ierr > 0 ) EXIT_ON_ERROR
    END DO
    
    ! never reach this point
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,I17,A)') "Cannot find requested user ptr", requested_len, tag
#endif
    EXIT_ON_ERROR
    
  END SUBROUTINE reattach_usr_ptr

  SUBROUTINE detach_usr_ptr( usr_ptr, this, requested_len, tag, ierr )
    CHARACTER, POINTER, CONTIGUOUS, INTENT( IN ) :: usr_ptr( : )
    TYPE( memory_pool ), INTENT( INOUT )         :: this
    INTEGER( INT64 ), INTENT( IN )               :: requested_len
    CHARACTER( 256 ), INTENT( IN )               :: tag
    INTEGER( INT32 ), INTENT( OUT )              :: ierr

    INTEGER( INT32)                              :: id

#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,I17,A)') "Trying to find requested user ptr", requested_len, tag
#endif

    DO id = 1, SIZE( this%node, 1 )
       CALL save_segment( usr_ptr, this%node( id )%ptr, requested_len, tag, ierr )
       !ierr < 0 nothing found, got to next node
       !ierr > 0 some fortran error appeared
       !ierr = 0 user_ptr is associated with some segment
       IF( ierr == 0 ) RETURN
       IF( ierr > 0 ) EXIT_ON_ERROR
    END DO

    !never reach this point
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,I17,A)') "Cannot find requested user ptr", requested_len, tag
#endif
    EXIT_ON_ERROR
    
  END SUBROUTINE detach_usr_ptr

  SUBROUTINE alloc_usr_ptr( usr_ptr, this, requested_len, tag, ierr )
    CHARACTER, POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr( : )
    TYPE( memory_pool ), INTENT( INOUT )          :: this
    INTEGER( INT64 ), INTENT( IN )                :: requested_len
    CHARACTER( 256 ), INTENT( IN )                :: tag
    INTEGER( INT32 ), INTENT( OUT )               :: ierr

    INTEGER( INT32 )                              :: id
    LOGICAL                                       :: free
    INTEGER( INT64 )                              :: temp

#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,I17,A)') "Trying to allocate requested user ptr", requested_len, tag
#endif

    DO id = 1, SIZE( this%node, 1 )
       CALL request_free_segment( usr_ptr, this%node( id )%ptr, requested_len, tag, ierr )
       !ierr < 0 nothing found, got to next node
       !ierr > 0 some fortran error appeared
       !ierr = 0 user_ptr is associated with some segment
       IF( ierr == 0 ) RETURN
       IF( ierr > 0 ) EXIT_ON_ERROR
    END DO

    DO id = SIZE( this%node, 1 ), 1, -1
       CALL is_free( this%node( id )%ptr, free )
       IF( free )THEN
          CALL remove_node_ptr( this%node, id, temp, ierr )
          IF( ierr /= 0 ) EXIT_ON_ERROR
       END IF
    END DO
    CALL add_node_ptr( this%node, requested_len, ierr )
    
    CALL request_free_segment( usr_ptr, this%node( SIZE( this%node, 1) )%ptr, requested_len, tag, ierr )
    IF( ierr == 0 ) RETURN
    IF( ierr > 0 ) EXIT_ON_ERROR
    IF( ierr < 0 ) EXIT_ON_ERROR

  END SUBROUTINE alloc_usr_ptr

  SUBROUTINE free_usr_ptr( usr_ptr, this, requested_len, tag, ierr )
    CHARACTER, POINTER, CONTIGUOUS, INTENT( IN ) :: usr_ptr( : )
    TYPE( memory_pool ), INTENT( INOUT )         :: this
    INTEGER( INT64 ), INTENT( IN )               :: requested_len
    CHARACTER( 256 ), INTENT( IN )               :: tag
    INTEGER( INT32 ), INTENT( OUT )              :: ierr

    INTEGER( INT32)                              :: id

#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,I17,A)') "Trying to find requested user ptr", requested_len, tag
#endif

    DO id = 1, SIZE( this%node, 1 )
       CALL free_segment( usr_ptr, this%node( id )%ptr, requested_len, tag, ierr )
       !ierr < 0 nothing found, got to next node
       !ierr > 0 some fortran error appeared
       !ierr = 0 user_ptr is associated with some segment
       IF( ierr == 0 ) EXIT
       IF( ierr > 0 ) EXIT_ON_ERROR
    END DO

    IF( ierr < 0 ) THEN
#if defined(_DEBUG)
       WRITE( OUTPUT_UNIT, '(A,I17,A)') "Cannot find requested user ptr", requested_len, tag
#endif
       EXIT_ON_ERROR
    END IF

    CALL cleanup_pool( this%node, ierr )

  END SUBROUTINE free_usr_ptr
  
  SUBROUTINE init_pool( this, requested_len, ierr )
    TYPE( memory_pool ), INTENT( OUT ) :: this
    INTEGER( INT64 ), INTENT( IN )     :: requested_len
    INTEGER( INT32 ), INTENT( OUT )    :: ierr

    CALL add_node_ptr( this%node, requested_len, ierr )
    
  END SUBROUTINE init_pool

  SUBROUTINE finalize_pool( this, pool_len, ierr )
    TYPE( memory_pool ), INTENT( INOUT ) :: this
    INTEGER( INT64 ), INTENT( OUT )      :: pool_len
    INTEGER( INT32 ), INTENT( OUT )      :: ierr

    INTEGER( INT32 )                     :: old_size
    INTEGER( INT32 )                     :: id
    INTEGER( INT64 )                     :: temp

    old_size = SIZE( this%node, 1 )
    pool_len = INT( 0, KIND = INT64 )
    
    DO id = old_size, 1, -1
       CALL remove_node_ptr( this%node, id, temp, ierr )
       IF( ierr /= 0 ) EXIT_ON_ERROR
       pool_len = pool_len + temp
    END DO

  END SUBROUTINE finalize_pool
  
  SUBROUTINE cleanup_pool( this, ierr )
    ! if more than one empty node is available, deallocate these nodes and reallocate a larger node
    ! since we call cleanup_pool every time we free a user ptr there can be at most two free nodes
    TYPE( pool_ptr ), ALLOCATABLE, INTENT( INOUT ) :: this( : )
    INTEGER( INT32 ), INTENT( OUT )                :: ierr

    LOGICAL                                        :: free
    INTEGER( INT64 )                               :: new_pool_len
    INTEGER( INT64 )                               :: temp
    INTEGER( INT32 )                               :: num_nodes
    INTEGER( INT32 )                               :: ids_free( 2 )
    INTEGER( INT32 )                               :: id
    INTEGER( INT32 )                               :: count_free

    num_nodes = SIZE( this, 1 )

    IF( num_nodes == 1 ) THEN
       ierr = 0
       RETURN
    END IF
    
    count_free = 0
    DO id = 1, num_nodes
       CALL is_free( this( id )%ptr, free )
       IF( free )THEN
          count_free = count_free +1
          ids_free( count_free ) = id
       END IF
    END DO
#if defined(_DEBUG)
    IF( count_free > 2 )THEN
       ierr = -1
       WRITE( OUTPUT_UNIT, '(A)') "This should never happen, more than two nodes are empty!"
       EXIT_ON_ERROR
    END IF
#endif
    IF( count_free == 2 )THEN
#if defined(_DEBUG)
       WRITE( OUTPUT_UNIT, '(A,2I4)' ) "Cleaning up the memory pool, concat id", ids_free
#endif
       new_pool_len = INT( 0, KIND = INT64 )
       !first remove the higher index
       DO id = 2, 1, -1
          CALL remove_node_ptr( this, ids_free( id ), temp, ierr )
          IF( ierr /= 0 ) EXIT_ON_ERROR
          new_pool_len = new_pool_len + temp
       END DO

       CALL add_node_ptr( this, new_pool_len, ierr )
       IF( ierr /= 0) EXIT_ON_ERROR
       
    END IF
    
  END SUBROUTINE cleanup_pool
  
  SUBROUTINE add_node_ptr( this, requested_len, ierr )
    ! adds a new memory node to the memory pool
    ! alings the data pointer
    ! initializes the index arrays
    TYPE( pool_ptr ), ALLOCATABLE, INTENT( INOUT ) :: this( : )
    INTEGER( INT64 ), INTENT( IN )                 :: requested_len
    INTEGER( INT32 ), INTENT( OUT )                :: ierr

    INTEGER( INT32 )                               :: num_nodes
    TYPE( pool_ptr ), ALLOCATABLE                  :: backup( : )


    num_nodes = 1
    IF( ALLOCATED( this ) )THEN
       num_nodes = SIZE( this, 1 ) + 1
       CALL backup_node_ptr( this, backup )
    END IF

    CALL allocate_node_ptr( this, num_nodes, ierr )
    IF ( ierr /= 0 ) EXIT_ON_ERROR

    CALL init_segment( this( num_nodes )%ptr, requested_len, ierr )
    IF ( ierr /= 0 ) EXIT_ON_ERROR
    
    IF( num_nodes > 1 )THEN
       CALL recover_node_ptr( this, backup, num_nodes, ierr )
       IF ( ierr /= 0 ) EXIT_ON_ERROR
    END IF

  END SUBROUTINE add_node_ptr

  SUBROUTINE remove_node_ptr( this, remove_id, removed_len, ierr )
    ! removes a memory node from the memory pool
    TYPE( pool_ptr ), ALLOCATABLE, INTENT( INOUT ) :: this( : )
    INTEGER( INT32 ), INTENT( IN )                 :: remove_id
    INTEGER( INT64 ), INTENT( OUT )                :: removed_len
    INTEGER( INT32 ), INTENT( OUT )                :: ierr

    INTEGER( INT32 )                               :: new_size
    TYPE( pool_ptr ), ALLOCATABLE                  :: backup( : )

    IF( .NOT. ALLOCATED( this ) ) THEN
       ierr = -1
       WRITE( OUTPUT_UNIT ,'(A)') "Cannot remove empty node_ptr"
       EXIT_ON_ERROR
    END IF

    CALL finalize_segment( this( remove_id )%ptr, removed_len, ierr )
    IF ( ierr /= 0 ) EXIT_ON_ERROR

    new_size = SIZE( this, 1 ) - 1
    IF( new_size > 0 ) THEN

       CALL backup_node_ptr( this, backup )

       ALLOCATE( this( new_size ), STAT = ierr )
       IF( ierr /= 0 )THEN
          WRITE( OUTPUT_UNIT, '(A1XI4)' ) "Cannot allocate node with size", new_size
       END IF
      
       CALL recover_node_ptr( this, backup, remove_id, ierr )
       IF ( ierr /= 0 ) EXIT_ON_ERROR
       
    ELSE

       DEALLOCATE( this( 1 )%ptr, STAT = ierr )
       IF( ierr /= 0 )THEN
          WRITE( OUTPUT_UNIT, '(A)' ) "Cannot deallocate last node_ptr"
          EXIT_ON_ERROR
       END IF

       DEALLOCATE( this, STAT = ierr )
       IF( ierr /= 0 )THEN
          WRITE( OUTPUT_UNIT, '(A)' ) "Cannot deallocate last node"
          EXIT_ON_ERROR
       END IF

    END IF

  END SUBROUTINE remove_node_ptr

  SUBROUTINE allocate_node_ptr( this, new_size, ierr )
    TYPE( pool_ptr ), ALLOCATABLE, INTENT( INOUT ) :: this( : )
    INTEGER( INT32 ), INTENT( IN )                 :: new_size
    INTEGER( INT32 ), INTENT( OUT )                :: ierr

    ALLOCATE( this( new_size ), STAT = ierr )
    IF( ierr /= 0 )THEN
       WRITE( OUTPUT_UNIT, '(A1XI4)' ) "Cannot allocate node with size", new_size
    END IF
    ALLOCATE( this( new_size )%ptr, STAT = ierr )
    IF( ierr /= 0 )THEN
       WRITE( OUTPUT_UNIT, '(A1XI4)' ) "Cannot allocate node_ptr with size", new_size
    END IF
       
  END SUBROUTINE allocate_node_ptr
  
  SUBROUTINE backup_node_ptr( this, this_backup )
    TYPE( pool_ptr ), ALLOCATABLE, INTENT( INOUT ) :: this(:)
    TYPE( pool_ptr ), ALLOCATABLE, INTENT( OUT )   :: this_backup(:)

    CALL MOVE_ALLOC( this, this_backup )
    
  END SUBROUTINE backup_node_ptr

  SUBROUTINE recover_node_ptr( this, this_backup, skip_id, ierr )
    TYPE( pool_ptr ), INTENT( INOUT )              :: this( : )
    TYPE( pool_ptr ), ALLOCATABLE, INTENT( INOUT ) :: this_backup( : )
    INTEGER( INT32 ), INTENT( IN )                 :: skip_id
    INTEGER( INT32 ), INTENT( OUT )                :: ierr

    INTEGER( INT32 )                               :: old_size
    INTEGER( INT32 )                               :: new_size
    INTEGER( INT32 )                               :: offset
    INTEGER( INT32 )                               :: i
    

    new_size = SIZE( this, 1 )
    old_size = SIZE( this_backup, 1 )
    DO i = 1, skip_id - 1
       CALL MOVE_ALLOC( this_backup( i )%ptr, this( i )%ptr )
    END DO
    
    IF( new_size > old_size )THEN
       offset = +1
    ELSE
       offset = -1
    END IF

    DO i = skip_id +1 , old_size
       CALL MOVE_ALLOC( this_backup( i )%ptr, this( i + offset )%ptr )
    END DO

#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, *) 'recover', skip_id, new_size, old_size
#endif
    IF( offset == -1 )THEN
       DEALLOCATE( this_backup( skip_id )%ptr, STAT = ierr )
       IF( ierr /= 0 )THEN
          WRITE( OUTPUT_UNIT, '(A)') "Could not deallocate memory pool"
          RETURN
       END IF
    END IF
       
    DEALLOCATE( this_backup, STAT = ierr )
    IF( ierr /= 0 )THEN
       WRITE( OUTPUT_UNIT, '(A)') "Could not deallocate memory pool"
       RETURN
    END IF

  END SUBROUTINE recover_node_ptr

END MODULE pool_managment_utils

#if defined(_TEST_POOL)
PROGRAM test_pool_managment_utils
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_LOC
  use data_managment_utils
  USE pool_managment_utils
  IMPLICIT NONE

  TYPE( memory_pool ) :: test_pool
  CHARACTER, POINTER, CONTIGUOUS :: test_ptr1( : ) => NULL()
  CHARACTER, POINTER, CONTIGUOUS :: test_ptr2( : ) => NULL()
  CHARACTER, POINTER, CONTIGUOUS :: test_ptr3( : ) => NULL()
  CHARACTER, POINTER, CONTIGUOUS :: test_ptr4( : ) => NULL()
  CHARACTER, POINTER, CONTIGUOUS :: test_ptr5( : ) => NULL()
  CHARACTER, POINTER, CONTIGUOUS :: test_ptr6( : ) => NULL()
  
  INTEGER( INT32 ) :: ierr
  INTEGER( INT64 ) :: len_req
  INTEGER( INT64 ) :: len_out
  CHARACTER( 256 ) :: tag1, tag2, tag3, tag4, tag5, tag6
  INTEGER( INT64 ) :: shift
  TYPE( C_PTR )                   :: c_addr
  INTEGER( INT64 )                :: faddr

  
  len_req = INT( 256, KIND = INT64 )
  len_req = len_req * len_req * len_req * INT( 10, KIND = INT64 ) + ONE_INT64

  CALL init_pool( test_pool, len_req, ierr )
  tag1 = "test_ptr1"
  tag2 = "test_ptr2"
  tag3 = "test_ptr3"
  tag4 = "test_ptr4"
  tag5 = "test_ptr5"
  tag6 = "test_ptr6"
  
  CALL alloc_usr_ptr( test_ptr1, test_pool, len_req, tag1, ierr )
  WRITE( OUTPUT_UNIT, '(A,L2)') "test_ptr1 associated?", ASSOCIATED(test_ptr1)
  WRITE( OUTPUT_UNIT, '(A,I17)') "test_ptr1 sizes?", SIZE(test_ptr1, 1, KIND = INT64 )
  c_addr = C_LOC( test_ptr1( 1 ) )
  faddr = TRANSFER( c_addr , faddr )
  shift = ( alignment - MOD( faddr , alignment ) ) / sizeof_char
  WRITE( OUTPUT_UNIT, '(A,I17)' ) "shift needed to aling:", shift

  CALL alloc_usr_ptr( test_ptr2, test_pool, len_req, tag2, ierr )
  WRITE( OUTPUT_UNIT, '(A,L2)') "test_ptr2 associated?", ASSOCIATED(test_ptr2)
  WRITE( OUTPUT_UNIT, '(A,I17)') "test_ptr2 sizes?", SIZE(test_ptr2, 1, KIND = INT64 )
  c_addr = C_LOC( test_ptr2( 1 ) )
  faddr = TRANSFER( c_addr , faddr )
  shift = ( alignment - MOD( faddr , alignment ) ) / sizeof_char
  WRITE( OUTPUT_UNIT, '(A,I17)' ) "shift needed to aling:", shift

  CALL detach_usr_ptr( test_ptr1, test_pool, len_req, tag1, ierr )
  test_ptr1 => NULL()
  
  CALL detach_usr_ptr( test_ptr2, test_pool, len_req, tag2, ierr )
  test_ptr2 => NULL()

  CALL reattach_usr_ptr( test_ptr1, test_pool, len_req, tag1, ierr )
  WRITE( OUTPUT_UNIT, '(A,L2)') "test_ptr1 associated?", ASSOCIATED(test_ptr1)
  WRITE( OUTPUT_UNIT, '(A,I17)') "test_ptr1 sizes?", SIZE(test_ptr1, 1, KIND = INT64 )
  
  CALL reattach_usr_ptr( test_ptr2, test_pool, len_req, tag2, ierr )
  WRITE( OUTPUT_UNIT, '(A,L2)') "test_ptr2 associated?", ASSOCIATED(test_ptr2)
  WRITE( OUTPUT_UNIT, '(A,I17)') "test_ptr2 sizes?", SIZE(test_ptr2, 1, KIND = INT64 )

  CALL free_usr_ptr( test_ptr1, test_pool, len_req, tag1, ierr )
  test_ptr1 => NULL()

  CALL free_usr_ptr( test_ptr2, test_pool, len_req, tag2, ierr )
  test_ptr2 => NULL()

  CALL alloc_usr_ptr( test_ptr1, test_pool, len_req, tag1, ierr )
  WRITE( OUTPUT_UNIT, '(A,L2)') "test_ptr1 associated?", ASSOCIATED(test_ptr1)
  WRITE( OUTPUT_UNIT, '(A,I17)') "test_ptr1 sizes?", SIZE(test_ptr1, 1, KIND = INT64 )
  c_addr = C_LOC( test_ptr1( 1 ) )
  faddr = TRANSFER( c_addr , faddr )
  shift = ( alignment - MOD( faddr , alignment ) ) / sizeof_char
  WRITE( OUTPUT_UNIT, '(A,I17)' ) "shift needed to aling:", shift

  CALL alloc_usr_ptr( test_ptr2, test_pool, len_req, tag2, ierr )
  WRITE( OUTPUT_UNIT, '(A,L2)') "test_ptr2 associated?", ASSOCIATED(test_ptr2)
  WRITE( OUTPUT_UNIT, '(A,I17)') "test_ptr2 sizes?", SIZE(test_ptr2, 1, KIND = INT64 )
  c_addr = C_LOC( test_ptr2( 1 ) )
  faddr = TRANSFER( c_addr , faddr )
  shift = ( alignment - MOD( faddr , alignment ) ) / sizeof_char
  WRITE( OUTPUT_UNIT, '(A,I17)' ) "shift needed to aling:", shift
  
  CALL free_usr_ptr( test_ptr1, test_pool, len_req, tag1, ierr )
  test_ptr1 => NULL()

  CALL free_usr_ptr( test_ptr2, test_pool, len_req, tag2, ierr )
  test_ptr2 => NULL()

  CALL alloc_usr_ptr( test_ptr1, test_pool, len_req / INT( 3, KIND = INT64 ) , tag1, ierr )
  CALL alloc_usr_ptr( test_ptr2, test_pool, len_req / INT( 14, KIND = INT64 ) , tag2, ierr )
  CALL alloc_usr_ptr( test_ptr3, test_pool, len_req / INT( 20, KIND = INT64 ) , tag3, ierr )
  CALL alloc_usr_ptr( test_ptr4, test_pool, len_req / INT( 15, KIND = INT64 ) , tag4, ierr )
  CALL alloc_usr_ptr( test_ptr5, test_pool, len_req / INT( 14, KIND = INT64 ) , tag5, ierr )
  CALL free_usr_ptr( test_ptr1, test_pool, len_req / INT( 3, KIND = INT64 )  , tag1, ierr )
  CALL free_usr_ptr( test_ptr4, test_pool, len_req / INT( 15, KIND = INT64 ) , tag4, ierr )
  CALL free_usr_ptr( test_ptr2, test_pool, len_req / INT( 14, KIND = INT64 ) , tag2, ierr )
  CALL free_usr_ptr( test_ptr5, test_pool, len_req / INT( 14, KIND = INT64 ) , tag5, ierr )
  CALL free_usr_ptr( test_ptr3, test_pool, len_req / INT( 20, KIND = INT64 ) , tag3, ierr )
  CALL finalize_pool( test_pool, len_out, ierr )
  WRITE( OUTPUT_UNIT, '(X,A,I17,X,A,I17,I4)') "len requested", len_req, "len allocated", len_out, ierr
  
END PROGRAM test_pool_managment_utils
#endif
