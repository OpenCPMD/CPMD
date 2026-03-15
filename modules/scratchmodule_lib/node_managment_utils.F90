SUBMODULE (scratch_data) data_managment_utils
  IMPLICIT NONE
  PRIVATE

  CHARACTER, PARAMETER :: guard_start = 'S'
  CHARACTER, PARAMETER :: guard_start = 'E'

CONTAINS
  
  SUBROUTINE allocate_node( this, requested_size, ierr )
    TYPE( scr ), POINTER, INTENT( OUT ) :: this
    INTEGER( INT64 ), INTENT( IN )      :: requested_size
    INTEGER( INT32 ), INTENT( OUT )     :: ierr


    ALLOCATE( this%array( requested_size ), STAT = ierr )
    IF ( ierr /= 0 ) THEN
       WRITE( OUTPUT_UNIT ,'A') "Could not allocate node"
       RETURN
    END IF

    this%array_size = requested_size

    CALL align_data( this )

    ierr = 0
    
  END SUBROUTINE allocate_node
  
  SUBROUTINE deallocate_node( this, ierr )
    TYPE( scr ), POINTER, INTENT( OUT ) :: this
    INTEGER( INT32 ), INTENT( OUT )     :: ierr

    IF( .NOT. ALLOCATED( this%array ) )THEN
       ierr = -1
       WRITE( OUTPUT_UNIT, 'A' ) "Can not deallocate unallocated array"
       RETURN
    END IF

    DEALLOCATE( this%array, STAT = ierr )
    IF ( ierr /= 0 ) THEN
       WRITE( OUTPUT_UNIT, 'A' ) "Can not deallocate node"
       RETURN
    END IF

    this%array_size = requested_size

    CALL align_data( this )

    ierr = 0
    
  END SUBROUTINE allocate_node

  SUBROUTINE align_array( this )
    !sets the pointer alinged_data to point to start at the first alingment boundary of data
    !and adjusts the lenght of the pointer to be a multiple of the alingment
    TYPE(scr), POINTER, INTENT(INOUT)     :: this

    INTEGER(INT64)                                   :: shift
    INTEGER(INT64)                                   :: start
    INTEGER(INT64)                                   :: shifted_start
    INTEGER(INT64)                                   :: shifted_end

    CALL calc_shift( this%array, shift )

    !substract the number of bytes to move the start to the next alignment boundary
    this%alinged_array_size = this%array_size - shift
    !substract the trailing block to get an alinged_size as a multiple of the alingment
    this%alinged_array_size = this%alinged_array_size - &
         MOD( this%alinged_array_size , INT( alignment_in_char, KIND = INT64 )
    start = 1
    shifted_start = start + shift
    shifted_end = alinged_size + shift

    this%alinged_array( start : alinged_size ) => this%array( shifted_start : shifted_end )

  END SUBROUTINE align_array

  SUBROUTINE calc_shift( arrayin, shift )
    ! calculates the padding needed to align the array to specified alignment
    ! in the worst case a padding of (alignment - 1) bytes is needed
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_LOC
    CHARACTER, TARGET, INTENT( IN ) :: arrayin(:)
    INTEGER(INT64), INTENT(OUT)        :: shift
    TYPE(C_PTR)                     :: c_addr
    INTEGER(INT64)                     :: faddr

    c_addr = C_LOC( arrayin( 1 ) )
    faddr = TRANSFER( c_addr , faddr )
    shift = ( alignment - MOD( faddr , alignment ) ) / sizeof_char

    IF( shift == alignment_in_char ) shift = 0

  END SUBROUTINE calc_shift

  SUBROUTINE set_guard( this, guards, ierr)
    TYPE( scr ), INTENT( INOUT )   :: this
    INTEGER( INT64 ), INTENT( IN ) :: guards(*)
    INTEGER( INT32 ), INTENT( OUT ) :: ierr

    INTEGER(INT64)               :: guard_start
    INTEGER(INT64)               :: guard_end
    INTEGER(INT64)               :: i

    IF( guards(4) > this%aligned_array_size )THEN
       ierr = -1
       WRITE( OUTPUT_UNIT, 'A') "Requested end of guard block is outside of aligned_array"
       RETURN
    END IF

    guard_start = guards( 1 )
    guard_end = guards( 2 )
    DO i = guard_start, guard_end
       this%array_array( i ) = guard_start
    END DO
    
    guard_start = guards( 3 )
    guard_end = guards( 4 )
    DO i = guard_start, guard_end
       this%aligned_array( i ) = guard_end
    END DO

    ierr = 0

  END SUBROUTINE set_guard
  
  SUBROUTINE check_guard( this, guards, ierr)
    TYPE( scr ), INTENT( INOUT )    :: this
    INTEGER( INT64 ), INTENT( IN )  :: guards(*)
    INTEGER( INT32 ), INTENT( OUT ) :: ierr

    INTEGER(INT64)                  :: guard_start
    INTEGER(INT64)                  :: guard_end
    INTEGER( INT64 )                :: i
    INTEGER( INT32 )                :: ierr1
    INTEGER( INT32 )                :: ierr2
    
    ierr1 = 0
    ierr2 = 0
    ierr = 0
    
    guard_start = guards( 1 )
    guard_end = guards( 2 )

    DO i = guard_start, guard_end
       IF( .NOT. this%aligned_arry( i ) == guard_start )THEN
          ierr1 = -1
       END IF
    END DO
    
    IF( ierr1 /= 0 ) THEN
       WRITE( OUTPUT_UNIT ,'A') "Warning memory corrupted in front of segment"
       ierr = -1
    END IF

    guard_start = guards( 3 )
    guard_end = guards( 4 )

    DO i = guard_start, guard_end
       IF( .NOT. this%aligned_array( i ) == guard_end ) THEN
          ierr2 = -1
       END IF
    END DO
    
    IF( ierr2 /= 0 ) THEN
       WRITE( OUTPUT_UNIT , 'A' ) "Warning memory corrupted after segment"
       ierr = -1
    END IF

  END SUBROUTINE check_guard

  SUBROUTINE associate_ptr( this, ptr_out, pos, ierr)
    TYPE( scr ), INTENT( INOUT )                  :: this
    CHARACTER, POINTER, CONTIGUOUS, INTENT( OUT ) :: ptr_out(:)
    INTEGER( INT64 ), INTENT( IN )                :: pos(*)
    INTEGER( INT32 ), INTENT( OUT )               :: ierr

    INTEGER( INT64 )                              :: start
    INTEGER( INT64 )                              :: end
    INTEGER( INT64 )                              :: len

    start = pos( 1 )
    end = pos( 2 )
    len = pos( 3 )

    IF( end > this%aligned_array_size ) THEN
       ierr = -1
       WRITE( OUTPUT_UNIT, 'A' ) "Cannot associate pointer, requested len is larger than available data"
       RETURN
    END IF

    ptr_out( 1 : len ) => this%alinged_array( start : end )
    CALL check_associate_ptr( this, ptr_out, start, end, ierr )
    
  END SUBROUTINE associate_ptr

  SUBROUTINE check_associate_ptr( this, ptr_in, pos, ierr)
    TYPE( scr ), INTENT( INOUT )                 :: this
    CHARACTER, POINTER, CONTIGUOUS, INTENT( IN ) :: ptr_in(:)
    INTEGER( INT64 ), INTENT( IN )               :: pos(*)
    INTEGER( INT32 ), INTENT( OUT )              :: ierr

    CHARACTER, POINTER, CONTIGUOUS               :: temp_ptr(:)
    INTEGER( INT64 )                              :: start
    INTEGER( INT64 )                              :: end
    INTEGER( INT64 )                              :: len

    start = pos( 1 )
    end = pos( 2 )
    len = pos( 3 )

    IF( end > this%aligned_array_size ) THEN
       ierr = -1
       WRITE( OUTPUT_UNIT, 'A' ) "Cannot check pointer association, requested len is larger than"
       WRITE( OUTPUT_UNIT, 'A' ) "available data"
       RETURN
    END IF

    temp_ptr( 1 : len ) => this%alinged_array( start : end )
    IF( ASSOCIATED( ptr_in, temp_ptr ) )THEN
       ierr = 0
    ELSE
       ierr = -1
       WRITE( OUTPUT_UNIT, 'A' ) "Tested pointer not associated with target segment"
    END IF
    
  END SUBROUTINE check_associate_ptr
  
END SUBMODULE data_managment_utils
