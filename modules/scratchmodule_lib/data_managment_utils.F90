#include "scratch_defines.h"

MODULE data_managment_utils
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  IMPLICIT NONE
  PRIVATE
  

  INTEGER( INT64 ), PUBLIC, PARAMETER :: alignment = INT( _ALIGNMENT, KIND = INT64 ) !< alignment in bytes
  INTEGER( INT64 ), PUBLIC, PARAMETER :: sizeof_char = INT( character_storage_size / 8, KIND = INT64 ) !< size of char in bytes
  INTEGER( INT64 ), PUBLIC, PARAMETER :: alignment_in_char = alignment / sizeof_char !< alignment in char
  INTEGER( INT64 ), PUBLIC, PARAMETER :: pad = INT( _PADDING, KIND = INT64 ) * alignment !< padding between segments in bytes
  INTEGER( INT64 ), PUBLIC, PARAMETER :: pad_in_char = pad / sizeof_char !<padding between segments in dp

  INTEGER( INT64 ), PUBLIC, PARAMETER :: One_INT64 = INT( 1, KIND = INT64 )
  INTEGER( INT64 ), PUBLIC, PARAMETER :: Zero_INT64 = INT( 0, KIND = INT64 )

  CHARACTER( LEN = 1 ), PARAMETER :: guard_start = 'S'
  CHARACTER( LEN = 1 ), PARAMETER :: guard_end = 'E'

  TYPE, PUBLIC :: scr
     CHARACTER, PRIVATE, ALLOCATABLE         :: array(:) !< actual allocated data array, character array => 1 element 1 byte
     CHARACTER, PRIVATE, POINTER, CONTIGUOUS :: aligned_array(:) => NULL() !< pointer to be used when accessing allocated data, points to the alignment boundary, length must be a multiple of alignment!
     INTEGER(INT64), PRIVATE                 :: array_size !< actual allocated node size
     INTEGER(INT64), PRIVATE                 :: aligned_array_size !< size remaining after shifting to alignment boundary
  END TYPE scr

  PUBLIC :: allocate_array
  PUBLIC :: deallocate_array
  PUBLIC :: associate_ptr
  PUBLIC :: check_associate_ptr
  
CONTAINS

  SUBROUTINE associate_ptr( ptr_out, this, pos, ierr )
    CHARACTER, POINTER, CONTIGUOUS, INTENT( OUT ) :: ptr_out( : )
    TYPE( scr ), INTENT( INOUT )                  :: this
    INTEGER( INT64 ), INTENT( IN )                :: pos( * )
    INTEGER( INT32 ), INTENT( OUT )               :: ierr

    INTEGER( INT64 )                              :: start
    INTEGER( INT64 )                              :: end
    INTEGER( INT64 )                              :: len
#if defined( _DEBUG )
    INTEGER( INT64 )                              :: shift
#endif
    
    start = pos( 1 )
    end = pos( 2 )
    len = pos( 3 )
    ierr = 0
    
    IF( end + alignment_in_char + alignment_in_char - &
         MOD( len, alignment_in_char) > this%aligned_array_size  ) THEN
       ierr = -1
       WRITE( OUTPUT_UNIT, '(A,3I17)' ) "Cannot associate pointer, requested len is larger than &
            available data", end, this%aligned_array_size, &
            this%aligned_array_size - alignment_in_char - alignment_in_char &
            + MOD( len, alignment_in_char )
       EXIT_ON_ERROR
    END IF

#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A)' ) "Size of array"
    WRITE( OUTPUT_UNIT, '(2I17)' ) SIZE( this%array, 1 , KIND = INT64 ), this%array_size

    WRITE( OUTPUT_UNIT, '(A)' ) "Size of alinged_array ptr"
    WRITE( OUTPUT_UNIT, '(2I17)' ) SIZE( this%aligned_array, 1 , KIND = INT64 ), this%aligned_array_size
#endif
    ptr_out( 1 : len ) => this%aligned_array( start : end )
    CALL set_guard( this%aligned_array, pos )
    
#if defined( _DEBUG)
    CALL calc_shift( ptr_out, shift )
    WRITE( OUTPUT_UNIT, '(A,4I17)') "associate, start, end, len, shift", start, end, len, shift
#endif
    
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A)' ) "Size of array"
    WRITE( OUTPUT_UNIT, '(2I17)' ) SIZE( this%array, 1 , KIND = INT64), this%array_size

    WRITE( OUTPUT_UNIT, '(A)' ) "Size of alinged_array ptr"
    WRITE( OUTPUT_UNIT, '(2I17)' ) SIZE( this%aligned_array, 1 , KIND = INT64), this%aligned_array_size
#endif

  END SUBROUTINE associate_ptr

  SUBROUTINE check_associate_ptr( ptr_in, this, pos, ierr )
    CHARACTER, POINTER, CONTIGUOUS, INTENT( IN ) :: ptr_in( : )
    TYPE( scr ), INTENT( INOUT )                 :: this
    INTEGER( INT64 ), INTENT( IN )               :: pos( * )
    INTEGER( INT32 ), INTENT( OUT )              :: ierr

    CHARACTER, POINTER, CONTIGUOUS               :: temp_ptr( : ) => NULL()
    INTEGER( INT64 )                             :: start
    INTEGER( INT64 )                             :: end
    INTEGER( INT64 )                             :: len

    start = pos( 1 )
    end = pos( 2 )
    len = pos( 3 )

    IF( end + alignment_in_char + alignment_in_char - &
         MOD( len, alignment_in_char) > this%aligned_array_size  ) THEN
       ierr = -1
       WRITE( OUTPUT_UNIT, '(A)' ) "Cannot check pointer association, requested len is larger than"
       WRITE( OUTPUT_UNIT, '(A,2I17)' ) "available data", end, this%aligned_array_size
       EXIT_ON_ERROR
    END IF

#if defined( _DEBUG)
    WRITE( OUTPUT_UNIT, '(A,3I17)') "check_associate, start, end, len", start, end, len
#endif

    temp_ptr( 1 : len ) => this%aligned_array( start : end )

    IF( ASSOCIATED( ptr_in, temp_ptr ) )THEN
       ierr = 0
    ELSE
       ierr = -100
       WRITE( OUTPUT_UNIT, '(A)' ) "Tested pointer not associated with target segment"
       EXIT_ON_ERROR
    END IF

    CALL check_guard( this%aligned_array, pos, ierr )
    
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A)' ) "Size of array"
    WRITE( OUTPUT_UNIT, '(2I17)' ) SIZE( this%array, 1 , KIND = INT64 ), this%array_size

    WRITE( OUTPUT_UNIT, '(A)' ) "Size of alinged_array ptr"
    WRITE( OUTPUT_UNIT, '(2I17)' ) SIZE( this%aligned_array, 1 , KIND = INT64 ), this%aligned_array_size
#endif

  END SUBROUTINE check_associate_ptr
  
  SUBROUTINE allocate_array( this, requested_size, ierr )
    TYPE( scr ), INTENT( OUT )      :: this
    INTEGER( INT64 ), INTENT( IN )  :: requested_size
    INTEGER( INT32 ), INTENT( OUT ) :: ierr

    INTEGER( INT64 )                :: new_size
    
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(6I4)' ) selected_int_kind( 15 ), INT64
#endif

    !in the worst case we need to add alignment_in_char bytes to be able to shift to the next
    !alignment boundary additionally we also want that every segment is a multiple of
    !alignment_in_char
    !we also add a guard block in the front and in the end of the array
    new_size = &
         !payload
         requested_size &
         !for shifting to the next alignment block
         + alignment_in_char &
         !payload extended to be a multiple of alignment
         + alignment_in_char - MOD( requested_size, alignment_in_char) &
         !segment will be padded in front and at the end
         + pad_in_char + alignment_in_char + alignment_in_char &
         !guard blocks for the aligned_array_pointer
         + pad_in_char + alignment_in_char + alignment_in_char
    
    ALLOCATE( this%array( new_size ), STAT = ierr )
    IF ( ierr /= 0 ) THEN
       WRITE( OUTPUT_UNIT ,'(A)') "Could not allocate node"
       EXIT_ON_ERROR
    END IF
#if defined(__HAS_OMP_OFFLOAD)
    !$omp target enter data map(alloc:this%array)
#endif
    this%array_size = new_size

    CALL align_array( this, ierr )
    IF( ierr /= 0 ) EXIT_ON_ERROR
    
    IF(( this%aligned_array_size  - pad_in_char - alignment_in_char - alignment_in_char ) &
         < requested_size )THEN
       ierr = -1
       WRITE( OUTPUT_UNIT, '(A)') "Something went wrong in the calculation of the arraylength"
    END IF
    
  END SUBROUTINE allocate_array
  
  SUBROUTINE deallocate_array( this, len, ierr )
    TYPE( scr ), INTENT( INOUT )    :: this
    INTEGER( INT64 ), INTENT( OUT ) :: len
    INTEGER( INT32 ), INTENT( OUT ) :: ierr

    INTEGER( INT64 )                :: pos( 3 )
    INTEGER( INT64 )                :: shift

    CALL calc_shift( this%array, shift )
    
    pos( 1 ) = shift + pad_in_char + ONE_INT64
    pos( 2 ) = this%aligned_array_size + shift + pad_in_char
    pos( 3 ) = this%aligned_array_size
    CALL check_guard( this%array, pos, ierr )
    IF( ierr /= 0 ) EXIT_ON_ERROR
    
    IF( .NOT. ALLOCATED( this%array ) )THEN
       ierr = -1
       WRITE( OUTPUT_UNIT, '(A)' ) "Can not deallocate unallocated array"
       EXIT_ON_ERROR
    END IF
#if defined(__HAS_OMP_OFFLOAD)
    !$omp target exit data map(delete:this%array)
#endif
    DEALLOCATE( this%array, STAT = ierr )
    IF ( ierr /= 0 ) THEN
       WRITE( OUTPUT_UNIT, '(A)' ) "Can not deallocate node"
       EXIT_ON_ERROR
    END IF

    len = this%array_size
    
    this%array_size = -1

  END SUBROUTINE deallocate_array
  
  SUBROUTINE align_array( this , ierr )
    !sets the pointer aligned_arry to point to start at the first aligment boundary of data
    !and adjusts the lenght of the pointer to be a multiple of the aligment
    TYPE( scr ), INTENT( INOUT )    :: this
    INTEGER( INT32 ), INTENT( OUT ) :: ierr
    
    INTEGER( INT64 )                :: shift
    INTEGER( INT64 )                :: pos_aligned(3)
    INTEGER( INT64 )                :: pos_temp(3)
    CHARACTER, POINTER, CONTIGUOUS  :: temp(:) => NULL()
    
    CALL calc_shift( this%array, shift )

    !substract the number of bytes to move the start to the next alignment boundary
    this%aligned_array_size = this%array_size - shift
    !add alignment - mod(size, alignment) to ensure length is a multiple of alignment_in_char
    this%aligned_array_size = this%aligned_array_size - alignment_in_char - alignment_in_char - &
         MOD( this%aligned_array_size , alignment_in_char ) - pad_in_char

    pos_aligned( 1 ) = shift + pad_in_char + ONE_INT64
    pos_aligned( 2 ) = this%aligned_array_size + shift + pad_in_char
    pos_aligned( 3 ) = this%aligned_array_size

#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,7I17)') 'Align array, pos_alinged', pos_aligned, pad_in_char, shift, &
         this%aligned_array_size, this%array_size
#endif

    IF( this%array_size < pos_aligned( 2 ) )THEN
       ierr = -1
       WRITE( OUTPUT_UNIT, '(A,2I17)' ) "Alinged array outside of allocated array", this%array_size, &
            pos_aligned( 2 )
       EXIT_ON_ERROR
    END IF
    
    CALL set_aligned_ptr( this%aligned_array, this%array, pos_aligned, ierr )
    IF( ierr /= 0 ) EXIT_ON_ERROR

    !set guard blocks for the aligned_array pointer
    CALL set_guard( this%array, pos_aligned )

    !set a dummy pointer
    pos_temp( 1 ) = pad_in_char + ONE_INT64
    pos_temp( 2 ) = this%aligned_array_size - alignment_in_char - alignment_in_char
    pos_temp( 3 ) = this%aligned_array_size - alignment_in_char - alignment_in_char - pad_in_char
    
    CALL associate_ptr( temp, this, pos_temp, ierr )
    IF( ierr /= 0 ) EXIT_ON_ERROR

    CALL check_associate_ptr( temp, this, pos_temp, ierr )
    IF( ierr /= 0 ) EXIT_ON_ERROR

    !set guard blocks for this dummy pointer
    CALL set_guard( this%aligned_array, pos_temp )
    !check guard blocks of the alinged_array pointer
    CALL check_guard( this%array, pos_aligned, ierr )
    IF( ierr /= 0 )EXIT_ON_ERROR
    NULLIFY( temp )
    
  END SUBROUTINE align_array

  SUBROUTINE set_aligned_ptr( ptr_out, this, pos, ierr )
    !special routine for setting the alinged_array pointer to point to the allocatable array
    CHARACTER, POINTER, CONTIGUOUS, INTENT( OUT ) :: ptr_out( : )
    CHARACTER, TARGET, INTENT( IN )               :: this( * )
    INTEGER( INT64 ), INTENT( IN )                :: pos( * )
    INTEGER( INT32 ), INTENT( OUT )               :: ierr
    
    INTEGER( INT64 )                              :: start
    INTEGER( INT64 )                              :: end
    INTEGER( INT64 )                              :: len

    start = pos( 1 )
    end = pos( 2 )
    len = pos( 2 ) - pos( 1 ) + ONE_INT64

    IF( len /= pos( 3 ) )THEN
       ierr = -1
       WRITE( OUTPUT_UNIT, '(A,4I17)' ) "Mismatch of pointer length requested", len, pos( 1:3 )
       EXIT_ON_ERROR
    END IF

    ptr_out( 1 : len ) => this( start : end )

    ierr = 0
    
  END SUBROUTINE set_aligned_ptr
  
  SUBROUTINE calc_shift( arrayin, shift )
    ! calculates the padding needed to align the array to specified alignment
    ! in the worst case a padding of (alignment - 1) bytes is needed
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_LOC
    CHARACTER, TARGET, INTENT( IN ) :: arrayin( : )
    INTEGER( INT64 ), INTENT( OUT ) :: shift
    TYPE( C_PTR )                   :: c_addr
    INTEGER( INT64 )                :: faddr

    c_addr = C_LOC( arrayin( 1 ) )
    faddr = TRANSFER( c_addr , faddr )
    shift = ( alignment - MOD( faddr , alignment ) ) / sizeof_char

    IF( shift == alignment_in_char ) shift = 0

#if defined( _DEBUG)
    WRITE( OUTPUT_UNIT, '(A,X,I16)') "shift", shift
#endif
  END SUBROUTINE calc_shift

  SUBROUTINE set_guard( this, pos )
    CHARACTER, INTENT( INOUT )      :: this( * )
    INTEGER( INT64 ), INTENT( IN )  :: pos( * )

    INTEGER( INT64 )                :: start
    INTEGER( INT64 )                :: end
    INTEGER( INT64 )                :: i

    !padding is used as guard block
    start = pos( 1 ) - pad_in_char
    end = pos( 1 ) - ONE_INT64
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,2I17)' ) "set start guard, start, end", start, end
#endif
    DO i = start, end
       this( i ) = guard_start
    END DO

    !alignment is used as guard block
    start = pos( 2 ) + ONE_INT64
    end =  pos( 2 ) + alignment_in_char + alignment_in_char - MOD( pos( 3 ), alignment_in_char )
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,2I17)' ) "set end guard, start, end", start, end
#endif
    DO i = start, end
       this( i ) = guard_end
    END DO

  END SUBROUTINE set_guard
  
  SUBROUTINE check_guard( this, pos, ierr )
    CHARACTER, INTENT( INOUT )      :: this( * )
    INTEGER( INT64 ), INTENT( IN )  :: pos( * )
    INTEGER( INT32 ), INTENT( OUT ) :: ierr

    INTEGER( INT64 )                :: start
    INTEGER( INT64 )                :: end
    INTEGER( INT64 )                :: i
    INTEGER( INT32 )                :: ierr1
    INTEGER( INT32 )                :: ierr2
    
    ierr1 = 0
    ierr2 = 0
    ierr = 0
    
    start = pos( 1 ) - pad_in_char
    end = pos( 1 ) - ONE_INT64
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,2I17)' ) "check start guard, start, end", start, end
#endif
    DO i = start, end
       IF( .NOT. this( i ) == guard_start )THEN
          ierr1 = -1
       END IF
    END DO

    IF( ierr1 /= 0 ) THEN
       WRITE( OUTPUT_UNIT ,'(A)') "Warning memory corrupted in front of segment"
       ierr = -99
       EXIT_ON_ERROR
    END IF

    start = pos( 2 ) + ONE_INT64
    end =  pos( 2 ) + alignment_in_char + alignment_in_char - MOD( pos( 3 ), alignment_in_char )
#if defined(_DEBUG)
    WRITE( OUTPUT_UNIT, '(A,3I17)' ) "check end guard, start, end", start, end, alignment_in_char
#endif
    DO i = start, end
       IF( .NOT. this( i ) == guard_end ) THEN
          ierr2 = -1
       END IF
    END DO
    
    IF( ierr2 /= 0 ) THEN
       WRITE( OUTPUT_UNIT , '(A)' ) "Warning memory corrupted after segment"
       ierr = -99
       EXIT_ON_ERROR
    END IF

  END SUBROUTINE check_guard

  
END MODULE data_managment_utils

#if defined(_TEST_DATA)
PROGRAM test_data_managment_utils
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  USE data_managment_utils
  IMPLICIT NONE

  TYPE( scr ) :: test_array
  INTEGER( INT32 ) :: ierr
  INTEGER( INT64 ) :: len_req
  INTEGER( INT64 ) :: len_out

  len_req = INT( 1024, KIND = INT64 ) * INT( 1024, KIND = INT64 ) * INT( 1024, KIND = INT64 ) * &
       INT( 10, KIND = INT64 )
  WRITE( OUTPUT_UNIT, '(A,3I17)' ) "allocated len =", len_req 
  CALL allocate_array( test_array, len_req, ierr )
  WRITE( OUTPUT_UNIT, '(A,3I17)' ) "allocated len =", len_req
  CALL deallocate_array( test_array, len_out, ierr )
  IF( len_req /= len_out ) WRITE( OUTPUT_UNIT, '(A,2I17)') "mismatch between requested and allocated &
       length", len_req, len_out

END PROGRAM test_data_managment_utils
#endif
