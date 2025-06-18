  SUBROUTINE request_internal_scratch_shape_dc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_dc1
  
  SUBROUTINE request_external_scratch_shape_dc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_dc1

  SUBROUTINE request_internal_scratch_bounds_dc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_dc1
  
  SUBROUTINE request_external_scratch_bounds_dc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_dc1

  SUBROUTINE request_saved_internal_scratch_shape_dc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_dc1
  
  SUBROUTINE request_saved_external_scratch_shape_dc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_dc1

  SUBROUTINE request_saved_internal_scratch_bounds_dc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_dc1

  SUBROUTINE request_saved_external_scratch_bounds_dc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_dc1

  SUBROUTINE free_internal_scratch_shape_dc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_dc1

  SUBROUTINE free_external_scratch_shape_dc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_dc1
  
  SUBROUTINE free_internal_scratch_bounds_dc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_dc1
  
  SUBROUTINE free_external_scratch_bounds_dc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_dc1

  SUBROUTINE save_internal_scratch_shape_dc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_dc1
  
  SUBROUTINE save_external_scratch_shape_dc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_dc1

  SUBROUTINE save_internal_scratch_bounds_dc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_dc1
  
  SUBROUTINE save_external_scratch_bounds_dc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_dc1
  SUBROUTINE request_internal_scratch_shape_dc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_dc2
  
  SUBROUTINE request_external_scratch_shape_dc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_dc2

  SUBROUTINE request_internal_scratch_bounds_dc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_dc2
  
  SUBROUTINE request_external_scratch_bounds_dc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_dc2

  SUBROUTINE request_saved_internal_scratch_shape_dc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_dc2
  
  SUBROUTINE request_saved_external_scratch_shape_dc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_dc2

  SUBROUTINE request_saved_internal_scratch_bounds_dc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_dc2

  SUBROUTINE request_saved_external_scratch_bounds_dc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_dc2

  SUBROUTINE free_internal_scratch_shape_dc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_dc2

  SUBROUTINE free_external_scratch_shape_dc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_dc2
  
  SUBROUTINE free_internal_scratch_bounds_dc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_dc2
  
  SUBROUTINE free_external_scratch_bounds_dc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_dc2

  SUBROUTINE save_internal_scratch_shape_dc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_dc2
  
  SUBROUTINE save_external_scratch_shape_dc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_dc2

  SUBROUTINE save_internal_scratch_bounds_dc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_dc2
  
  SUBROUTINE save_external_scratch_bounds_dc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_dc2
  SUBROUTINE request_internal_scratch_shape_dc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_dc3
  
  SUBROUTINE request_external_scratch_shape_dc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_dc3

  SUBROUTINE request_internal_scratch_bounds_dc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_dc3
  
  SUBROUTINE request_external_scratch_bounds_dc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_dc3

  SUBROUTINE request_saved_internal_scratch_shape_dc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_dc3
  
  SUBROUTINE request_saved_external_scratch_shape_dc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_dc3

  SUBROUTINE request_saved_internal_scratch_bounds_dc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_dc3

  SUBROUTINE request_saved_external_scratch_bounds_dc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_dc3

  SUBROUTINE free_internal_scratch_shape_dc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_dc3

  SUBROUTINE free_external_scratch_shape_dc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_dc3
  
  SUBROUTINE free_internal_scratch_bounds_dc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_dc3
  
  SUBROUTINE free_external_scratch_bounds_dc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_dc3

  SUBROUTINE save_internal_scratch_shape_dc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_dc3
  
  SUBROUTINE save_external_scratch_shape_dc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_dc3

  SUBROUTINE save_internal_scratch_bounds_dc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_dc3
  
  SUBROUTINE save_external_scratch_bounds_dc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_dc3
  SUBROUTINE request_internal_scratch_shape_dc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_dc4
  
  SUBROUTINE request_external_scratch_shape_dc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_dc4

  SUBROUTINE request_internal_scratch_bounds_dc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_dc4
  
  SUBROUTINE request_external_scratch_bounds_dc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_dc4

  SUBROUTINE request_saved_internal_scratch_shape_dc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_dc4
  
  SUBROUTINE request_saved_external_scratch_shape_dc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_dc4

  SUBROUTINE request_saved_internal_scratch_bounds_dc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_dc4

  SUBROUTINE request_saved_external_scratch_bounds_dc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_dc4

  SUBROUTINE free_internal_scratch_shape_dc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_dc4

  SUBROUTINE free_external_scratch_shape_dc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_dc4
  
  SUBROUTINE free_internal_scratch_bounds_dc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_dc4
  
  SUBROUTINE free_external_scratch_bounds_dc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_dc4

  SUBROUTINE save_internal_scratch_shape_dc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_dc4
  
  SUBROUTINE save_external_scratch_shape_dc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_dc4

  SUBROUTINE save_internal_scratch_bounds_dc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_dc4
  
  SUBROUTINE save_external_scratch_bounds_dc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_dc4
  SUBROUTINE request_internal_scratch_shape_dc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_dc5
  
  SUBROUTINE request_external_scratch_shape_dc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_dc5

  SUBROUTINE request_internal_scratch_bounds_dc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_dc5
  
  SUBROUTINE request_external_scratch_bounds_dc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_dc5

  SUBROUTINE request_saved_internal_scratch_shape_dc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_dc5
  
  SUBROUTINE request_saved_external_scratch_shape_dc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_dc5

  SUBROUTINE request_saved_internal_scratch_bounds_dc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_dc5

  SUBROUTINE request_saved_external_scratch_bounds_dc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_dc5

  SUBROUTINE free_internal_scratch_shape_dc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_dc5

  SUBROUTINE free_external_scratch_shape_dc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_dc5
  
  SUBROUTINE free_internal_scratch_bounds_dc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_dc5
  
  SUBROUTINE free_external_scratch_bounds_dc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_dc5

  SUBROUTINE save_internal_scratch_shape_dc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_dc5
  
  SUBROUTINE save_external_scratch_shape_dc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_dc5

  SUBROUTINE save_internal_scratch_bounds_dc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_dc5
  
  SUBROUTINE save_external_scratch_bounds_dc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_dc5
  SUBROUTINE request_internal_scratch_shape_dc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_dc6
  
  SUBROUTINE request_external_scratch_shape_dc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_dc6

  SUBROUTINE request_internal_scratch_bounds_dc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_dc6
  
  SUBROUTINE request_external_scratch_bounds_dc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_dc6

  SUBROUTINE request_saved_internal_scratch_shape_dc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_dc6
  
  SUBROUTINE request_saved_external_scratch_shape_dc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_dc6

  SUBROUTINE request_saved_internal_scratch_bounds_dc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_dc6

  SUBROUTINE request_saved_external_scratch_bounds_dc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_dc6

  SUBROUTINE free_internal_scratch_shape_dc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_dc6

  SUBROUTINE free_external_scratch_shape_dc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_dc6
  
  SUBROUTINE free_internal_scratch_bounds_dc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_dc6
  
  SUBROUTINE free_external_scratch_bounds_dc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_dc6

  SUBROUTINE save_internal_scratch_shape_dc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_dc6
  
  SUBROUTINE save_external_scratch_shape_dc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_dc6

  SUBROUTINE save_internal_scratch_bounds_dc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_dc6
  
  SUBROUTINE save_external_scratch_bounds_dc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    COMPLEX( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_dc6
  SUBROUTINE request_internal_scratch_shape_rc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_rc1
  
  SUBROUTINE request_external_scratch_shape_rc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_rc1

  SUBROUTINE request_internal_scratch_bounds_rc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_rc1
  
  SUBROUTINE request_external_scratch_bounds_rc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_rc1

  SUBROUTINE request_saved_internal_scratch_shape_rc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_rc1
  
  SUBROUTINE request_saved_external_scratch_shape_rc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_rc1

  SUBROUTINE request_saved_internal_scratch_bounds_rc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_rc1

  SUBROUTINE request_saved_external_scratch_bounds_rc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_rc1

  SUBROUTINE free_internal_scratch_shape_rc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_rc1

  SUBROUTINE free_external_scratch_shape_rc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_rc1
  
  SUBROUTINE free_internal_scratch_bounds_rc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_rc1
  
  SUBROUTINE free_external_scratch_bounds_rc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_rc1

  SUBROUTINE save_internal_scratch_shape_rc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_rc1
  
  SUBROUTINE save_external_scratch_shape_rc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_rc1

  SUBROUTINE save_internal_scratch_bounds_rc1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_rc1
  
  SUBROUTINE save_external_scratch_bounds_rc1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_rc1
  SUBROUTINE request_internal_scratch_shape_rc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_rc2
  
  SUBROUTINE request_external_scratch_shape_rc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_rc2

  SUBROUTINE request_internal_scratch_bounds_rc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_rc2
  
  SUBROUTINE request_external_scratch_bounds_rc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_rc2

  SUBROUTINE request_saved_internal_scratch_shape_rc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_rc2
  
  SUBROUTINE request_saved_external_scratch_shape_rc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_rc2

  SUBROUTINE request_saved_internal_scratch_bounds_rc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_rc2

  SUBROUTINE request_saved_external_scratch_bounds_rc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_rc2

  SUBROUTINE free_internal_scratch_shape_rc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_rc2

  SUBROUTINE free_external_scratch_shape_rc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_rc2
  
  SUBROUTINE free_internal_scratch_bounds_rc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_rc2
  
  SUBROUTINE free_external_scratch_bounds_rc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_rc2

  SUBROUTINE save_internal_scratch_shape_rc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_rc2
  
  SUBROUTINE save_external_scratch_shape_rc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_rc2

  SUBROUTINE save_internal_scratch_bounds_rc2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_rc2
  
  SUBROUTINE save_external_scratch_bounds_rc2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_rc2
  SUBROUTINE request_internal_scratch_shape_rc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_rc3
  
  SUBROUTINE request_external_scratch_shape_rc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_rc3

  SUBROUTINE request_internal_scratch_bounds_rc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_rc3
  
  SUBROUTINE request_external_scratch_bounds_rc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_rc3

  SUBROUTINE request_saved_internal_scratch_shape_rc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_rc3
  
  SUBROUTINE request_saved_external_scratch_shape_rc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_rc3

  SUBROUTINE request_saved_internal_scratch_bounds_rc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_rc3

  SUBROUTINE request_saved_external_scratch_bounds_rc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_rc3

  SUBROUTINE free_internal_scratch_shape_rc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_rc3

  SUBROUTINE free_external_scratch_shape_rc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_rc3
  
  SUBROUTINE free_internal_scratch_bounds_rc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_rc3
  
  SUBROUTINE free_external_scratch_bounds_rc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_rc3

  SUBROUTINE save_internal_scratch_shape_rc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_rc3
  
  SUBROUTINE save_external_scratch_shape_rc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_rc3

  SUBROUTINE save_internal_scratch_bounds_rc3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_rc3
  
  SUBROUTINE save_external_scratch_bounds_rc3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_rc3
  SUBROUTINE request_internal_scratch_shape_rc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_rc4
  
  SUBROUTINE request_external_scratch_shape_rc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_rc4

  SUBROUTINE request_internal_scratch_bounds_rc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_rc4
  
  SUBROUTINE request_external_scratch_bounds_rc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_rc4

  SUBROUTINE request_saved_internal_scratch_shape_rc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_rc4
  
  SUBROUTINE request_saved_external_scratch_shape_rc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_rc4

  SUBROUTINE request_saved_internal_scratch_bounds_rc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_rc4

  SUBROUTINE request_saved_external_scratch_bounds_rc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_rc4

  SUBROUTINE free_internal_scratch_shape_rc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_rc4

  SUBROUTINE free_external_scratch_shape_rc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_rc4
  
  SUBROUTINE free_internal_scratch_bounds_rc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_rc4
  
  SUBROUTINE free_external_scratch_bounds_rc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_rc4

  SUBROUTINE save_internal_scratch_shape_rc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_rc4
  
  SUBROUTINE save_external_scratch_shape_rc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_rc4

  SUBROUTINE save_internal_scratch_bounds_rc4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_rc4
  
  SUBROUTINE save_external_scratch_bounds_rc4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_rc4
  SUBROUTINE request_internal_scratch_shape_rc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_rc5
  
  SUBROUTINE request_external_scratch_shape_rc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_rc5

  SUBROUTINE request_internal_scratch_bounds_rc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_rc5
  
  SUBROUTINE request_external_scratch_bounds_rc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_rc5

  SUBROUTINE request_saved_internal_scratch_shape_rc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_rc5
  
  SUBROUTINE request_saved_external_scratch_shape_rc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_rc5

  SUBROUTINE request_saved_internal_scratch_bounds_rc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_rc5

  SUBROUTINE request_saved_external_scratch_bounds_rc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_rc5

  SUBROUTINE free_internal_scratch_shape_rc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_rc5

  SUBROUTINE free_external_scratch_shape_rc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_rc5
  
  SUBROUTINE free_internal_scratch_bounds_rc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_rc5
  
  SUBROUTINE free_external_scratch_bounds_rc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_rc5

  SUBROUTINE save_internal_scratch_shape_rc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_rc5
  
  SUBROUTINE save_external_scratch_shape_rc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_rc5

  SUBROUTINE save_internal_scratch_bounds_rc5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_rc5
  
  SUBROUTINE save_external_scratch_bounds_rc5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_rc5
  SUBROUTINE request_internal_scratch_shape_rc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_rc6
  
  SUBROUTINE request_external_scratch_shape_rc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_rc6

  SUBROUTINE request_internal_scratch_bounds_rc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_rc6
  
  SUBROUTINE request_external_scratch_bounds_rc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_rc6

  SUBROUTINE request_saved_internal_scratch_shape_rc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_rc6
  
  SUBROUTINE request_saved_external_scratch_shape_rc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_rc6

  SUBROUTINE request_saved_internal_scratch_bounds_rc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_rc6

  SUBROUTINE request_saved_external_scratch_bounds_rc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_rc6

  SUBROUTINE free_internal_scratch_shape_rc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_rc6

  SUBROUTINE free_external_scratch_shape_rc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_rc6
  
  SUBROUTINE free_internal_scratch_bounds_rc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_rc6
  
  SUBROUTINE free_external_scratch_bounds_rc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_rc6

  SUBROUTINE save_internal_scratch_shape_rc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_rc6
  
  SUBROUTINE save_external_scratch_shape_rc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_rc6

  SUBROUTINE save_internal_scratch_bounds_rc6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_rc6
  
  SUBROUTINE save_external_scratch_bounds_rc6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    COMPLEX( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_rc6
  SUBROUTINE request_internal_scratch_shape_d1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_d1
  
  SUBROUTINE request_external_scratch_shape_d1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_d1

  SUBROUTINE request_internal_scratch_bounds_d1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_d1
  
  SUBROUTINE request_external_scratch_bounds_d1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_d1

  SUBROUTINE request_saved_internal_scratch_shape_d1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_d1
  
  SUBROUTINE request_saved_external_scratch_shape_d1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_d1

  SUBROUTINE request_saved_internal_scratch_bounds_d1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_d1

  SUBROUTINE request_saved_external_scratch_bounds_d1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_d1

  SUBROUTINE free_internal_scratch_shape_d1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_d1

  SUBROUTINE free_external_scratch_shape_d1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_d1
  
  SUBROUTINE free_internal_scratch_bounds_d1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_d1
  
  SUBROUTINE free_external_scratch_bounds_d1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_d1

  SUBROUTINE save_internal_scratch_shape_d1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_d1
  
  SUBROUTINE save_external_scratch_shape_d1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_d1

  SUBROUTINE save_internal_scratch_bounds_d1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_d1
  
  SUBROUTINE save_external_scratch_bounds_d1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_d1
  SUBROUTINE request_internal_scratch_shape_d2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_d2
  
  SUBROUTINE request_external_scratch_shape_d2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_d2

  SUBROUTINE request_internal_scratch_bounds_d2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_d2
  
  SUBROUTINE request_external_scratch_bounds_d2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_d2

  SUBROUTINE request_saved_internal_scratch_shape_d2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_d2
  
  SUBROUTINE request_saved_external_scratch_shape_d2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_d2

  SUBROUTINE request_saved_internal_scratch_bounds_d2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_d2

  SUBROUTINE request_saved_external_scratch_bounds_d2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_d2

  SUBROUTINE free_internal_scratch_shape_d2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_d2

  SUBROUTINE free_external_scratch_shape_d2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_d2
  
  SUBROUTINE free_internal_scratch_bounds_d2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_d2
  
  SUBROUTINE free_external_scratch_bounds_d2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_d2

  SUBROUTINE save_internal_scratch_shape_d2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_d2
  
  SUBROUTINE save_external_scratch_shape_d2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_d2

  SUBROUTINE save_internal_scratch_bounds_d2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_d2
  
  SUBROUTINE save_external_scratch_bounds_d2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_d2
  SUBROUTINE request_internal_scratch_shape_d3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_d3
  
  SUBROUTINE request_external_scratch_shape_d3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_d3

  SUBROUTINE request_internal_scratch_bounds_d3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_d3
  
  SUBROUTINE request_external_scratch_bounds_d3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_d3

  SUBROUTINE request_saved_internal_scratch_shape_d3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_d3
  
  SUBROUTINE request_saved_external_scratch_shape_d3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_d3

  SUBROUTINE request_saved_internal_scratch_bounds_d3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_d3

  SUBROUTINE request_saved_external_scratch_bounds_d3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_d3

  SUBROUTINE free_internal_scratch_shape_d3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_d3

  SUBROUTINE free_external_scratch_shape_d3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_d3
  
  SUBROUTINE free_internal_scratch_bounds_d3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_d3
  
  SUBROUTINE free_external_scratch_bounds_d3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_d3

  SUBROUTINE save_internal_scratch_shape_d3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_d3
  
  SUBROUTINE save_external_scratch_shape_d3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_d3

  SUBROUTINE save_internal_scratch_bounds_d3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_d3
  
  SUBROUTINE save_external_scratch_bounds_d3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_d3
  SUBROUTINE request_internal_scratch_shape_d4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_d4
  
  SUBROUTINE request_external_scratch_shape_d4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_d4

  SUBROUTINE request_internal_scratch_bounds_d4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_d4
  
  SUBROUTINE request_external_scratch_bounds_d4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_d4

  SUBROUTINE request_saved_internal_scratch_shape_d4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_d4
  
  SUBROUTINE request_saved_external_scratch_shape_d4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_d4

  SUBROUTINE request_saved_internal_scratch_bounds_d4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_d4

  SUBROUTINE request_saved_external_scratch_bounds_d4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_d4

  SUBROUTINE free_internal_scratch_shape_d4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_d4

  SUBROUTINE free_external_scratch_shape_d4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_d4
  
  SUBROUTINE free_internal_scratch_bounds_d4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_d4
  
  SUBROUTINE free_external_scratch_bounds_d4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_d4

  SUBROUTINE save_internal_scratch_shape_d4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_d4
  
  SUBROUTINE save_external_scratch_shape_d4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_d4

  SUBROUTINE save_internal_scratch_bounds_d4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_d4
  
  SUBROUTINE save_external_scratch_bounds_d4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_d4
  SUBROUTINE request_internal_scratch_shape_d5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_d5
  
  SUBROUTINE request_external_scratch_shape_d5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_d5

  SUBROUTINE request_internal_scratch_bounds_d5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_d5
  
  SUBROUTINE request_external_scratch_bounds_d5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_d5

  SUBROUTINE request_saved_internal_scratch_shape_d5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_d5
  
  SUBROUTINE request_saved_external_scratch_shape_d5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_d5

  SUBROUTINE request_saved_internal_scratch_bounds_d5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_d5

  SUBROUTINE request_saved_external_scratch_bounds_d5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_d5

  SUBROUTINE free_internal_scratch_shape_d5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_d5

  SUBROUTINE free_external_scratch_shape_d5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_d5
  
  SUBROUTINE free_internal_scratch_bounds_d5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_d5
  
  SUBROUTINE free_external_scratch_bounds_d5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_d5

  SUBROUTINE save_internal_scratch_shape_d5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_d5
  
  SUBROUTINE save_external_scratch_shape_d5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_d5

  SUBROUTINE save_internal_scratch_bounds_d5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_d5
  
  SUBROUTINE save_external_scratch_bounds_d5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_d5
  SUBROUTINE request_internal_scratch_shape_d6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_d6
  
  SUBROUTINE request_external_scratch_shape_d6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_d6

  SUBROUTINE request_internal_scratch_bounds_d6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_d6
  
  SUBROUTINE request_external_scratch_bounds_d6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_d6

  SUBROUTINE request_saved_internal_scratch_shape_d6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_d6
  
  SUBROUTINE request_saved_external_scratch_shape_d6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_d6

  SUBROUTINE request_saved_internal_scratch_bounds_d6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_d6

  SUBROUTINE request_saved_external_scratch_bounds_d6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_d6

  SUBROUTINE free_internal_scratch_shape_d6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_d6

  SUBROUTINE free_external_scratch_shape_d6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_d6
  
  SUBROUTINE free_internal_scratch_bounds_d6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_d6
  
  SUBROUTINE free_external_scratch_bounds_d6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_d6

  SUBROUTINE save_internal_scratch_shape_d6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_d6
  
  SUBROUTINE save_external_scratch_shape_d6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_d6

  SUBROUTINE save_internal_scratch_bounds_d6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_d6
  
  SUBROUTINE save_external_scratch_bounds_d6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    REAL( REAL64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_d6
  SUBROUTINE request_internal_scratch_shape_r1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_r1
  
  SUBROUTINE request_external_scratch_shape_r1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_r1

  SUBROUTINE request_internal_scratch_bounds_r1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_r1
  
  SUBROUTINE request_external_scratch_bounds_r1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_r1

  SUBROUTINE request_saved_internal_scratch_shape_r1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_r1
  
  SUBROUTINE request_saved_external_scratch_shape_r1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_r1

  SUBROUTINE request_saved_internal_scratch_bounds_r1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_r1

  SUBROUTINE request_saved_external_scratch_bounds_r1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_r1

  SUBROUTINE free_internal_scratch_shape_r1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_r1

  SUBROUTINE free_external_scratch_shape_r1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_r1
  
  SUBROUTINE free_internal_scratch_bounds_r1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_r1
  
  SUBROUTINE free_external_scratch_bounds_r1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_r1

  SUBROUTINE save_internal_scratch_shape_r1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_r1
  
  SUBROUTINE save_external_scratch_shape_r1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_r1

  SUBROUTINE save_internal_scratch_bounds_r1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_r1
  
  SUBROUTINE save_external_scratch_bounds_r1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_r1
  SUBROUTINE request_internal_scratch_shape_r2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_r2
  
  SUBROUTINE request_external_scratch_shape_r2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_r2

  SUBROUTINE request_internal_scratch_bounds_r2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_r2
  
  SUBROUTINE request_external_scratch_bounds_r2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_r2

  SUBROUTINE request_saved_internal_scratch_shape_r2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_r2
  
  SUBROUTINE request_saved_external_scratch_shape_r2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_r2

  SUBROUTINE request_saved_internal_scratch_bounds_r2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_r2

  SUBROUTINE request_saved_external_scratch_bounds_r2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_r2

  SUBROUTINE free_internal_scratch_shape_r2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_r2

  SUBROUTINE free_external_scratch_shape_r2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_r2
  
  SUBROUTINE free_internal_scratch_bounds_r2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_r2
  
  SUBROUTINE free_external_scratch_bounds_r2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_r2

  SUBROUTINE save_internal_scratch_shape_r2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_r2
  
  SUBROUTINE save_external_scratch_shape_r2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_r2

  SUBROUTINE save_internal_scratch_bounds_r2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_r2
  
  SUBROUTINE save_external_scratch_bounds_r2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_r2
  SUBROUTINE request_internal_scratch_shape_r3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_r3
  
  SUBROUTINE request_external_scratch_shape_r3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_r3

  SUBROUTINE request_internal_scratch_bounds_r3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_r3
  
  SUBROUTINE request_external_scratch_bounds_r3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_r3

  SUBROUTINE request_saved_internal_scratch_shape_r3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_r3
  
  SUBROUTINE request_saved_external_scratch_shape_r3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_r3

  SUBROUTINE request_saved_internal_scratch_bounds_r3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_r3

  SUBROUTINE request_saved_external_scratch_bounds_r3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_r3

  SUBROUTINE free_internal_scratch_shape_r3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_r3

  SUBROUTINE free_external_scratch_shape_r3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_r3
  
  SUBROUTINE free_internal_scratch_bounds_r3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_r3
  
  SUBROUTINE free_external_scratch_bounds_r3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_r3

  SUBROUTINE save_internal_scratch_shape_r3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_r3
  
  SUBROUTINE save_external_scratch_shape_r3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_r3

  SUBROUTINE save_internal_scratch_bounds_r3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_r3
  
  SUBROUTINE save_external_scratch_bounds_r3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_r3
  SUBROUTINE request_internal_scratch_shape_r4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_r4
  
  SUBROUTINE request_external_scratch_shape_r4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_r4

  SUBROUTINE request_internal_scratch_bounds_r4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_r4
  
  SUBROUTINE request_external_scratch_bounds_r4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_r4

  SUBROUTINE request_saved_internal_scratch_shape_r4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_r4
  
  SUBROUTINE request_saved_external_scratch_shape_r4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_r4

  SUBROUTINE request_saved_internal_scratch_bounds_r4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_r4

  SUBROUTINE request_saved_external_scratch_bounds_r4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_r4

  SUBROUTINE free_internal_scratch_shape_r4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_r4

  SUBROUTINE free_external_scratch_shape_r4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_r4
  
  SUBROUTINE free_internal_scratch_bounds_r4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_r4
  
  SUBROUTINE free_external_scratch_bounds_r4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_r4

  SUBROUTINE save_internal_scratch_shape_r4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_r4
  
  SUBROUTINE save_external_scratch_shape_r4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_r4

  SUBROUTINE save_internal_scratch_bounds_r4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_r4
  
  SUBROUTINE save_external_scratch_bounds_r4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_r4
  SUBROUTINE request_internal_scratch_shape_r5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_r5
  
  SUBROUTINE request_external_scratch_shape_r5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_r5

  SUBROUTINE request_internal_scratch_bounds_r5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_r5
  
  SUBROUTINE request_external_scratch_bounds_r5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_r5

  SUBROUTINE request_saved_internal_scratch_shape_r5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_r5
  
  SUBROUTINE request_saved_external_scratch_shape_r5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_r5

  SUBROUTINE request_saved_internal_scratch_bounds_r5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_r5

  SUBROUTINE request_saved_external_scratch_bounds_r5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_r5

  SUBROUTINE free_internal_scratch_shape_r5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_r5

  SUBROUTINE free_external_scratch_shape_r5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_r5
  
  SUBROUTINE free_internal_scratch_bounds_r5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_r5
  
  SUBROUTINE free_external_scratch_bounds_r5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_r5

  SUBROUTINE save_internal_scratch_shape_r5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_r5
  
  SUBROUTINE save_external_scratch_shape_r5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_r5

  SUBROUTINE save_internal_scratch_bounds_r5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_r5
  
  SUBROUTINE save_external_scratch_bounds_r5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_r5
  SUBROUTINE request_internal_scratch_shape_r6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_r6
  
  SUBROUTINE request_external_scratch_shape_r6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_r6

  SUBROUTINE request_internal_scratch_bounds_r6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_r6
  
  SUBROUTINE request_external_scratch_bounds_r6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_r6

  SUBROUTINE request_saved_internal_scratch_shape_r6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_r6
  
  SUBROUTINE request_saved_external_scratch_shape_r6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_r6

  SUBROUTINE request_saved_internal_scratch_bounds_r6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_r6

  SUBROUTINE request_saved_external_scratch_bounds_r6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    REAL( REAL32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_r6

  SUBROUTINE free_internal_scratch_shape_r6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_r6

  SUBROUTINE free_external_scratch_shape_r6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_r6
  
  SUBROUTINE free_internal_scratch_bounds_r6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_r6
  
  SUBROUTINE free_external_scratch_bounds_r6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_r6

  SUBROUTINE save_internal_scratch_shape_r6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_r6
  
  SUBROUTINE save_external_scratch_shape_r6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_r6

  SUBROUTINE save_internal_scratch_bounds_r6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_r6
  
  SUBROUTINE save_external_scratch_bounds_r6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    REAL( REAL32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_r6
  SUBROUTINE request_internal_scratch_shape_il1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_il1
  
  SUBROUTINE request_external_scratch_shape_il1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_il1

  SUBROUTINE request_internal_scratch_bounds_il1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_il1
  
  SUBROUTINE request_external_scratch_bounds_il1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_il1

  SUBROUTINE request_saved_internal_scratch_shape_il1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_il1
  
  SUBROUTINE request_saved_external_scratch_shape_il1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_il1

  SUBROUTINE request_saved_internal_scratch_bounds_il1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_il1

  SUBROUTINE request_saved_external_scratch_bounds_il1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_il1

  SUBROUTINE free_internal_scratch_shape_il1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_il1

  SUBROUTINE free_external_scratch_shape_il1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_il1
  
  SUBROUTINE free_internal_scratch_bounds_il1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_il1
  
  SUBROUTINE free_external_scratch_bounds_il1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_il1

  SUBROUTINE save_internal_scratch_shape_il1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_il1
  
  SUBROUTINE save_external_scratch_shape_il1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_il1

  SUBROUTINE save_internal_scratch_bounds_il1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_il1
  
  SUBROUTINE save_external_scratch_bounds_il1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_il1
  SUBROUTINE request_internal_scratch_shape_il2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_il2
  
  SUBROUTINE request_external_scratch_shape_il2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_il2

  SUBROUTINE request_internal_scratch_bounds_il2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_il2
  
  SUBROUTINE request_external_scratch_bounds_il2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_il2

  SUBROUTINE request_saved_internal_scratch_shape_il2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_il2
  
  SUBROUTINE request_saved_external_scratch_shape_il2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_il2

  SUBROUTINE request_saved_internal_scratch_bounds_il2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_il2

  SUBROUTINE request_saved_external_scratch_bounds_il2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_il2

  SUBROUTINE free_internal_scratch_shape_il2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_il2

  SUBROUTINE free_external_scratch_shape_il2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_il2
  
  SUBROUTINE free_internal_scratch_bounds_il2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_il2
  
  SUBROUTINE free_external_scratch_bounds_il2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_il2

  SUBROUTINE save_internal_scratch_shape_il2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_il2
  
  SUBROUTINE save_external_scratch_shape_il2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_il2

  SUBROUTINE save_internal_scratch_bounds_il2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_il2
  
  SUBROUTINE save_external_scratch_bounds_il2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_il2
  SUBROUTINE request_internal_scratch_shape_il3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_il3
  
  SUBROUTINE request_external_scratch_shape_il3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_il3

  SUBROUTINE request_internal_scratch_bounds_il3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_il3
  
  SUBROUTINE request_external_scratch_bounds_il3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_il3

  SUBROUTINE request_saved_internal_scratch_shape_il3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_il3
  
  SUBROUTINE request_saved_external_scratch_shape_il3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_il3

  SUBROUTINE request_saved_internal_scratch_bounds_il3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_il3

  SUBROUTINE request_saved_external_scratch_bounds_il3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_il3

  SUBROUTINE free_internal_scratch_shape_il3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_il3

  SUBROUTINE free_external_scratch_shape_il3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_il3
  
  SUBROUTINE free_internal_scratch_bounds_il3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_il3
  
  SUBROUTINE free_external_scratch_bounds_il3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_il3

  SUBROUTINE save_internal_scratch_shape_il3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_il3
  
  SUBROUTINE save_external_scratch_shape_il3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_il3

  SUBROUTINE save_internal_scratch_bounds_il3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_il3
  
  SUBROUTINE save_external_scratch_bounds_il3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_il3
  SUBROUTINE request_internal_scratch_shape_il4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_il4
  
  SUBROUTINE request_external_scratch_shape_il4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_il4

  SUBROUTINE request_internal_scratch_bounds_il4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_il4
  
  SUBROUTINE request_external_scratch_bounds_il4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_il4

  SUBROUTINE request_saved_internal_scratch_shape_il4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_il4
  
  SUBROUTINE request_saved_external_scratch_shape_il4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_il4

  SUBROUTINE request_saved_internal_scratch_bounds_il4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_il4

  SUBROUTINE request_saved_external_scratch_bounds_il4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_il4

  SUBROUTINE free_internal_scratch_shape_il4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_il4

  SUBROUTINE free_external_scratch_shape_il4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_il4
  
  SUBROUTINE free_internal_scratch_bounds_il4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_il4
  
  SUBROUTINE free_external_scratch_bounds_il4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_il4

  SUBROUTINE save_internal_scratch_shape_il4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_il4
  
  SUBROUTINE save_external_scratch_shape_il4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_il4

  SUBROUTINE save_internal_scratch_bounds_il4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_il4
  
  SUBROUTINE save_external_scratch_bounds_il4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_il4
  SUBROUTINE request_internal_scratch_shape_il5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_il5
  
  SUBROUTINE request_external_scratch_shape_il5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_il5

  SUBROUTINE request_internal_scratch_bounds_il5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_il5
  
  SUBROUTINE request_external_scratch_bounds_il5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_il5

  SUBROUTINE request_saved_internal_scratch_shape_il5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_il5
  
  SUBROUTINE request_saved_external_scratch_shape_il5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_il5

  SUBROUTINE request_saved_internal_scratch_bounds_il5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_il5

  SUBROUTINE request_saved_external_scratch_bounds_il5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_il5

  SUBROUTINE free_internal_scratch_shape_il5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_il5

  SUBROUTINE free_external_scratch_shape_il5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_il5
  
  SUBROUTINE free_internal_scratch_bounds_il5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_il5
  
  SUBROUTINE free_external_scratch_bounds_il5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_il5

  SUBROUTINE save_internal_scratch_shape_il5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_il5
  
  SUBROUTINE save_external_scratch_shape_il5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_il5

  SUBROUTINE save_internal_scratch_bounds_il5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_il5
  
  SUBROUTINE save_external_scratch_bounds_il5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_il5
  SUBROUTINE request_internal_scratch_shape_il6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_il6
  
  SUBROUTINE request_external_scratch_shape_il6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_il6

  SUBROUTINE request_internal_scratch_bounds_il6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_il6
  
  SUBROUTINE request_external_scratch_bounds_il6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_il6

  SUBROUTINE request_saved_internal_scratch_shape_il6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_il6
  
  SUBROUTINE request_saved_external_scratch_shape_il6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_il6

  SUBROUTINE request_saved_internal_scratch_bounds_il6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_il6

  SUBROUTINE request_saved_external_scratch_bounds_il6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT64 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_il6

  SUBROUTINE free_internal_scratch_shape_il6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_il6

  SUBROUTINE free_external_scratch_shape_il6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_il6
  
  SUBROUTINE free_internal_scratch_bounds_il6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_il6
  
  SUBROUTINE free_external_scratch_bounds_il6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_il6

  SUBROUTINE save_internal_scratch_shape_il6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_il6
  
  SUBROUTINE save_external_scratch_shape_il6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_il6

  SUBROUTINE save_internal_scratch_bounds_il6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_il6
  
  SUBROUTINE save_external_scratch_bounds_il6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    INTEGER( INT64 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_il6
  SUBROUTINE request_internal_scratch_shape_i1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_i1
  
  SUBROUTINE request_external_scratch_shape_i1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_i1

  SUBROUTINE request_internal_scratch_bounds_i1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_i1
  
  SUBROUTINE request_external_scratch_bounds_i1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_i1

  SUBROUTINE request_saved_internal_scratch_shape_i1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_i1
  
  SUBROUTINE request_saved_external_scratch_shape_i1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_i1

  SUBROUTINE request_saved_internal_scratch_bounds_i1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_i1

  SUBROUTINE request_saved_external_scratch_bounds_i1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(1,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_i1

  SUBROUTINE free_internal_scratch_shape_i1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_i1

  SUBROUTINE free_external_scratch_shape_i1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_i1
  
  SUBROUTINE free_internal_scratch_bounds_i1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_i1
  
  SUBROUTINE free_external_scratch_bounds_i1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_i1

  SUBROUTINE save_internal_scratch_shape_i1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_i1
  
  SUBROUTINE save_external_scratch_shape_i1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_i1

  SUBROUTINE save_internal_scratch_bounds_i1(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_i1
  
  SUBROUTINE save_external_scratch_bounds_i1(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(1,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 1
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_i1
  SUBROUTINE request_internal_scratch_shape_i2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_i2
  
  SUBROUTINE request_external_scratch_shape_i2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_i2

  SUBROUTINE request_internal_scratch_bounds_i2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_i2
  
  SUBROUTINE request_external_scratch_bounds_i2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_i2

  SUBROUTINE request_saved_internal_scratch_shape_i2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_i2
  
  SUBROUTINE request_saved_external_scratch_shape_i2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_i2

  SUBROUTINE request_saved_internal_scratch_bounds_i2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_i2

  SUBROUTINE request_saved_external_scratch_bounds_i2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(2,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_i2

  SUBROUTINE free_internal_scratch_shape_i2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_i2

  SUBROUTINE free_external_scratch_shape_i2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_i2
  
  SUBROUTINE free_internal_scratch_bounds_i2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_i2
  
  SUBROUTINE free_external_scratch_bounds_i2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_i2

  SUBROUTINE save_internal_scratch_shape_i2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_i2
  
  SUBROUTINE save_external_scratch_shape_i2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_i2

  SUBROUTINE save_internal_scratch_bounds_i2(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_i2
  
  SUBROUTINE save_external_scratch_bounds_i2(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(2,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 2
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_i2
  SUBROUTINE request_internal_scratch_shape_i3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_i3
  
  SUBROUTINE request_external_scratch_shape_i3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_i3

  SUBROUTINE request_internal_scratch_bounds_i3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_i3
  
  SUBROUTINE request_external_scratch_bounds_i3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_i3

  SUBROUTINE request_saved_internal_scratch_shape_i3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_i3
  
  SUBROUTINE request_saved_external_scratch_shape_i3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_i3

  SUBROUTINE request_saved_internal_scratch_bounds_i3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_i3

  SUBROUTINE request_saved_external_scratch_bounds_i3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(3,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_i3

  SUBROUTINE free_internal_scratch_shape_i3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_i3

  SUBROUTINE free_external_scratch_shape_i3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_i3
  
  SUBROUTINE free_internal_scratch_bounds_i3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_i3
  
  SUBROUTINE free_external_scratch_bounds_i3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_i3

  SUBROUTINE save_internal_scratch_shape_i3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_i3
  
  SUBROUTINE save_external_scratch_shape_i3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_i3

  SUBROUTINE save_internal_scratch_bounds_i3(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_i3
  
  SUBROUTINE save_external_scratch_bounds_i3(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(3,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 3
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_i3
  SUBROUTINE request_internal_scratch_shape_i4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_i4
  
  SUBROUTINE request_external_scratch_shape_i4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_i4

  SUBROUTINE request_internal_scratch_bounds_i4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_i4
  
  SUBROUTINE request_external_scratch_bounds_i4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_i4

  SUBROUTINE request_saved_internal_scratch_shape_i4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_i4
  
  SUBROUTINE request_saved_external_scratch_shape_i4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_i4

  SUBROUTINE request_saved_internal_scratch_bounds_i4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_i4

  SUBROUTINE request_saved_external_scratch_bounds_i4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(4,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_i4

  SUBROUTINE free_internal_scratch_shape_i4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_i4

  SUBROUTINE free_external_scratch_shape_i4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_i4
  
  SUBROUTINE free_internal_scratch_bounds_i4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_i4
  
  SUBROUTINE free_external_scratch_bounds_i4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_i4

  SUBROUTINE save_internal_scratch_shape_i4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_i4
  
  SUBROUTINE save_external_scratch_shape_i4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_i4

  SUBROUTINE save_internal_scratch_bounds_i4(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_i4
  
  SUBROUTINE save_external_scratch_bounds_i4(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(4,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 4
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_i4
  SUBROUTINE request_internal_scratch_shape_i5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_i5
  
  SUBROUTINE request_external_scratch_shape_i5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_i5

  SUBROUTINE request_internal_scratch_bounds_i5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_i5
  
  SUBROUTINE request_external_scratch_bounds_i5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_i5

  SUBROUTINE request_saved_internal_scratch_shape_i5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_i5
  
  SUBROUTINE request_saved_external_scratch_shape_i5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_i5

  SUBROUTINE request_saved_internal_scratch_bounds_i5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_i5

  SUBROUTINE request_saved_external_scratch_bounds_i5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(5,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_i5

  SUBROUTINE free_internal_scratch_shape_i5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_i5

  SUBROUTINE free_external_scratch_shape_i5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_i5
  
  SUBROUTINE free_internal_scratch_bounds_i5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_i5
  
  SUBROUTINE free_external_scratch_bounds_i5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_i5

  SUBROUTINE save_internal_scratch_shape_i5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_i5
  
  SUBROUTINE save_external_scratch_shape_i5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_i5

  SUBROUTINE save_internal_scratch_bounds_i5(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_i5
  
  SUBROUTINE save_external_scratch_bounds_i5(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(5,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 5
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_i5
  SUBROUTINE request_internal_scratch_shape_i6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_internal_scratch_shape_i6
  
  SUBROUTINE request_external_scratch_shape_i6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_shape_i6

  SUBROUTINE request_internal_scratch_bounds_i6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_internal_scratch_bounds_i6
  
  SUBROUTINE request_external_scratch_bounds_i6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )           :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL alloc_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_external_scratch_bounds_i6

  SUBROUTINE request_saved_internal_scratch_shape_i6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char
    
    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )
    
  END SUBROUTINE request_saved_internal_scratch_shape_i6
  
  SUBROUTINE request_saved_external_scratch_shape_i6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:)
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:)
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( 1:s(1),1:s(2),1:s(3),1:s(4),1:s(5),1:s(6) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_shape_i6

  SUBROUTINE request_saved_internal_scratch_bounds_i6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_internal_scratch_bounds_i6

  SUBROUTINE request_saved_external_scratch_bounds_i6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )        :: pool_ext
    INTEGER( INT64 ), INTENT( IN )              :: s(6,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( OUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                  :: tagin
    INTEGER( INT32 ), INTENT( OUT )             :: ierr
    
    TYPE( C_PTR )                               :: c_addr
    CHARACTER( 256 )                            :: tag
    CHARACTER, POINTER, CONTIGUOUS              :: char_ptr(:) => NULL()
    INTEGER( INT32 ), POINTER, CONTIGUOUS                :: type_ptr(:) => NULL()
    INTEGER( INT64 )                            :: len_in_char
    INTEGER( INT64 )                            :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                 :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )

    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       ALLOCATE(usr_ptr(s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2)),STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    CALL reattach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    c_addr = C_LOC( char_ptr )
    CALL C_F_POINTER(c_addr, type_ptr, (/len_in_type/))
    usr_ptr( s(1,1):s(1,2),s(2,1):s(2,2),s(3,1):s(3,2),s(4,1):s(4,2),s(5,1):s(5,2),s(6,1):s(6,2) )&
         => type_ptr
    NULLIFY( type_ptr, char_ptr )

  END SUBROUTINE request_saved_external_scratch_bounds_i6

  SUBROUTINE free_internal_scratch_shape_i6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_internal_scratch_shape_i6

  SUBROUTINE free_external_scratch_shape_i6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
                                   STORAGE_SIZE( usr_ptr ) / character_storage_size

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))

    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )
    
  END SUBROUTINE free_external_scratch_shape_i6
  
  SUBROUTINE free_internal_scratch_bounds_i6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr
    
    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_internal_scratch_bounds_i6
  
  SUBROUTINE free_external_scratch_bounds_i6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL free_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE free_external_scratch_bounds_i6

  SUBROUTINE save_internal_scratch_shape_i6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    !quick return if len <=0
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_shape_i6
  
  SUBROUTINE save_external_scratch_shape_i6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )      :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:)
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = PRODUCT( s )
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_shape_i6

  SUBROUTINE save_internal_scratch_bounds_i6(s, usr_ptr, tagin, ierr)
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_int, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_internal_scratch_bounds_i6
  
  SUBROUTINE save_external_scratch_bounds_i6(pool_ext, s, usr_ptr, tagin, ierr)
    TYPE( memory_pool ), INTENT( INOUT )       :: pool_ext
    INTEGER( INT64 ), INTENT( IN )             :: s(6,2)
    INTEGER( INT32 ), POINTER, CONTIGUOUS, INTENT( INOUT ) :: usr_ptr(:,:,:,:,:,:)
    CHARACTER(*), INTENT( IN )                 :: tagin
    INTEGER( INT32 ), INTENT( OUT )            :: ierr

    TYPE( C_PTR )                              :: c_addr
    CHARACTER( 256 )                           :: tag
    CHARACTER, POINTER, CONTIGUOUS             :: char_ptr(:) => NULL()
    INTEGER( INT64 )                           :: len_in_char
    INTEGER( INT64 )                           :: len_in_type
    INTEGER( INT32 )                            :: i
    INTEGER( INT32 ), PARAMETER                :: type_size_in_char = &
         STORAGE_SIZE( usr_ptr ) / character_storage_size 

    tag = TRIM( ADJUSTL( tagin ) )
    len_in_type = INT( 1, KIND = INT64 )
    DO i = 1, 6
       len_in_type = len_in_type * &
            ( ABS( s( i, 2 )  -  s( i, 1 ) ) + INT( 1, KIND = INT64 ) )
    END DO
    IF( len_in_type .LE. 0 )THEN
       DEALLOCATE(usr_ptr,STAT=ierr)
       RETURN
    END IF
    len_in_char = len_in_type * type_size_in_char

    c_addr = C_LOC( usr_ptr )
    CALL C_F_POINTER(c_addr, char_ptr, (/len_in_char/))
    CALL detach_usr_ptr( char_ptr, pool_ext, len_in_char, tag, ierr )
    NULLIFY( char_ptr )

  END SUBROUTINE save_external_scratch_bounds_i6
