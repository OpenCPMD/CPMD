MODULE scratch_interface
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_LOC, C_F_POINTER
  USE pool_managment_utils
  IMPLICIT NONE
  PUBLIC
 
  !Module structure:
  !pool_managment_submod.f90:
  ! hosts all routines working on type(pool_ptr)
  ! calls routines from segment_managment_submod.f90
  !segment_managment_submod.f90:
  ! hosts all routines working on type(segment)
  ! calls routines from data_managment_submod.f90
  !data_managment_submod.f90
  ! hosts all routines working on type(scr)

  !Flow: request_scratch -> scrach_interface ( len -> len in char ) -> scratch_managment
  ! check allocation status of memory_pool
  ! if not allocated allocate a new pool, init pool
  ! loop over all nodes in this pool:
  ! -> search_free_segment
  ! a) found perfect ( payload len free == len request )
  ! b) found larger ( payload len free > len request )
  ! c) not found
  ! if a) -> update segment -> associate pointer -> return
  ! if b) -> split segment (insert segment / allocate bigger index )
  !       -> set index and remaining index -> a)
  ! if c) -> allocate new node -> init index -> a)
  !
  !Flow: free_scratch -> scratch_interface ( len -> len in char ) -> search_occupied_segment
  ! a) found perfect -> check guard blocks -> update index -> return
  ! b) return error code
  !
  !Flow: save_scratch 
  ! -> flow: free_scratch
  !
  !Flow: request_saved_scratch
  ! -> flow: free_scartch

  INTERFACE init_memory_pool
     MODULE PROCEDURE :: init_internal_pool
     MODULE PROCEDURE :: init_external_pool
  END INTERFACE init_memory_pool

  INTERFACE finalize_memory_pool
     MODULE PROCEDURE :: finalize_internal_pool
     MODULE PROCEDURE :: finalize_external_pool     
  END INTERFACE finalize_memory_pool
  
  INTERFACE request_scratch
     MODULE PROCEDURE :: request_external_scratch_shape_dc1
     MODULE PROCEDURE :: request_external_scratch_bounds_dc1
     MODULE PROCEDURE :: request_external_scratch_shape_dc2
     MODULE PROCEDURE :: request_external_scratch_bounds_dc2
     MODULE PROCEDURE :: request_external_scratch_shape_dc3
     MODULE PROCEDURE :: request_external_scratch_bounds_dc3
     MODULE PROCEDURE :: request_external_scratch_shape_dc4
     MODULE PROCEDURE :: request_external_scratch_bounds_dc4
     MODULE PROCEDURE :: request_external_scratch_shape_dc5
     MODULE PROCEDURE :: request_external_scratch_bounds_dc5
     MODULE PROCEDURE :: request_external_scratch_shape_dc6
     MODULE PROCEDURE :: request_external_scratch_bounds_dc6
     MODULE PROCEDURE :: request_internal_scratch_shape_dc1
     MODULE PROCEDURE :: request_internal_scratch_bounds_dc1
     MODULE PROCEDURE :: request_internal_scratch_shape_dc2
     MODULE PROCEDURE :: request_internal_scratch_bounds_dc2
     MODULE PROCEDURE :: request_internal_scratch_shape_dc3
     MODULE PROCEDURE :: request_internal_scratch_bounds_dc3
     MODULE PROCEDURE :: request_internal_scratch_shape_dc4
     MODULE PROCEDURE :: request_internal_scratch_bounds_dc4
     MODULE PROCEDURE :: request_internal_scratch_shape_dc5
     MODULE PROCEDURE :: request_internal_scratch_bounds_dc5
     MODULE PROCEDURE :: request_internal_scratch_shape_dc6
     MODULE PROCEDURE :: request_internal_scratch_bounds_dc6
     MODULE PROCEDURE :: request_external_scratch_shape_rc1
     MODULE PROCEDURE :: request_external_scratch_bounds_rc1
     MODULE PROCEDURE :: request_external_scratch_shape_rc2
     MODULE PROCEDURE :: request_external_scratch_bounds_rc2
     MODULE PROCEDURE :: request_external_scratch_shape_rc3
     MODULE PROCEDURE :: request_external_scratch_bounds_rc3
     MODULE PROCEDURE :: request_external_scratch_shape_rc4
     MODULE PROCEDURE :: request_external_scratch_bounds_rc4
     MODULE PROCEDURE :: request_external_scratch_shape_rc5
     MODULE PROCEDURE :: request_external_scratch_bounds_rc5
     MODULE PROCEDURE :: request_external_scratch_shape_rc6
     MODULE PROCEDURE :: request_external_scratch_bounds_rc6
     MODULE PROCEDURE :: request_internal_scratch_shape_rc1
     MODULE PROCEDURE :: request_internal_scratch_bounds_rc1
     MODULE PROCEDURE :: request_internal_scratch_shape_rc2
     MODULE PROCEDURE :: request_internal_scratch_bounds_rc2
     MODULE PROCEDURE :: request_internal_scratch_shape_rc3
     MODULE PROCEDURE :: request_internal_scratch_bounds_rc3
     MODULE PROCEDURE :: request_internal_scratch_shape_rc4
     MODULE PROCEDURE :: request_internal_scratch_bounds_rc4
     MODULE PROCEDURE :: request_internal_scratch_shape_rc5
     MODULE PROCEDURE :: request_internal_scratch_bounds_rc5
     MODULE PROCEDURE :: request_internal_scratch_shape_rc6
     MODULE PROCEDURE :: request_internal_scratch_bounds_rc6
     MODULE PROCEDURE :: request_external_scratch_shape_d1
     MODULE PROCEDURE :: request_external_scratch_bounds_d1
     MODULE PROCEDURE :: request_external_scratch_shape_d2
     MODULE PROCEDURE :: request_external_scratch_bounds_d2
     MODULE PROCEDURE :: request_external_scratch_shape_d3
     MODULE PROCEDURE :: request_external_scratch_bounds_d3
     MODULE PROCEDURE :: request_external_scratch_shape_d4
     MODULE PROCEDURE :: request_external_scratch_bounds_d4
     MODULE PROCEDURE :: request_external_scratch_shape_d5
     MODULE PROCEDURE :: request_external_scratch_bounds_d5
     MODULE PROCEDURE :: request_external_scratch_shape_d6
     MODULE PROCEDURE :: request_external_scratch_bounds_d6
     MODULE PROCEDURE :: request_internal_scratch_shape_d1
     MODULE PROCEDURE :: request_internal_scratch_bounds_d1
     MODULE PROCEDURE :: request_internal_scratch_shape_d2
     MODULE PROCEDURE :: request_internal_scratch_bounds_d2
     MODULE PROCEDURE :: request_internal_scratch_shape_d3
     MODULE PROCEDURE :: request_internal_scratch_bounds_d3
     MODULE PROCEDURE :: request_internal_scratch_shape_d4
     MODULE PROCEDURE :: request_internal_scratch_bounds_d4
     MODULE PROCEDURE :: request_internal_scratch_shape_d5
     MODULE PROCEDURE :: request_internal_scratch_bounds_d5
     MODULE PROCEDURE :: request_internal_scratch_shape_d6
     MODULE PROCEDURE :: request_internal_scratch_bounds_d6
     MODULE PROCEDURE :: request_external_scratch_shape_r1
     MODULE PROCEDURE :: request_external_scratch_bounds_r1
     MODULE PROCEDURE :: request_external_scratch_shape_r2
     MODULE PROCEDURE :: request_external_scratch_bounds_r2
     MODULE PROCEDURE :: request_external_scratch_shape_r3
     MODULE PROCEDURE :: request_external_scratch_bounds_r3
     MODULE PROCEDURE :: request_external_scratch_shape_r4
     MODULE PROCEDURE :: request_external_scratch_bounds_r4
     MODULE PROCEDURE :: request_external_scratch_shape_r5
     MODULE PROCEDURE :: request_external_scratch_bounds_r5
     MODULE PROCEDURE :: request_external_scratch_shape_r6
     MODULE PROCEDURE :: request_external_scratch_bounds_r6
     MODULE PROCEDURE :: request_internal_scratch_shape_r1
     MODULE PROCEDURE :: request_internal_scratch_bounds_r1
     MODULE PROCEDURE :: request_internal_scratch_shape_r2
     MODULE PROCEDURE :: request_internal_scratch_bounds_r2
     MODULE PROCEDURE :: request_internal_scratch_shape_r3
     MODULE PROCEDURE :: request_internal_scratch_bounds_r3
     MODULE PROCEDURE :: request_internal_scratch_shape_r4
     MODULE PROCEDURE :: request_internal_scratch_bounds_r4
     MODULE PROCEDURE :: request_internal_scratch_shape_r5
     MODULE PROCEDURE :: request_internal_scratch_bounds_r5
     MODULE PROCEDURE :: request_internal_scratch_shape_r6
     MODULE PROCEDURE :: request_internal_scratch_bounds_r6
     MODULE PROCEDURE :: request_external_scratch_shape_il1
     MODULE PROCEDURE :: request_external_scratch_bounds_il1
     MODULE PROCEDURE :: request_external_scratch_shape_il2
     MODULE PROCEDURE :: request_external_scratch_bounds_il2
     MODULE PROCEDURE :: request_external_scratch_shape_il3
     MODULE PROCEDURE :: request_external_scratch_bounds_il3
     MODULE PROCEDURE :: request_external_scratch_shape_il4
     MODULE PROCEDURE :: request_external_scratch_bounds_il4
     MODULE PROCEDURE :: request_external_scratch_shape_il5
     MODULE PROCEDURE :: request_external_scratch_bounds_il5
     MODULE PROCEDURE :: request_external_scratch_shape_il6
     MODULE PROCEDURE :: request_external_scratch_bounds_il6
     MODULE PROCEDURE :: request_internal_scratch_shape_il1
     MODULE PROCEDURE :: request_internal_scratch_bounds_il1
     MODULE PROCEDURE :: request_internal_scratch_shape_il2
     MODULE PROCEDURE :: request_internal_scratch_bounds_il2
     MODULE PROCEDURE :: request_internal_scratch_shape_il3
     MODULE PROCEDURE :: request_internal_scratch_bounds_il3
     MODULE PROCEDURE :: request_internal_scratch_shape_il4
     MODULE PROCEDURE :: request_internal_scratch_bounds_il4
     MODULE PROCEDURE :: request_internal_scratch_shape_il5
     MODULE PROCEDURE :: request_internal_scratch_bounds_il5
     MODULE PROCEDURE :: request_internal_scratch_shape_il6
     MODULE PROCEDURE :: request_internal_scratch_bounds_il6
     MODULE PROCEDURE :: request_external_scratch_shape_i1
     MODULE PROCEDURE :: request_external_scratch_bounds_i1
     MODULE PROCEDURE :: request_external_scratch_shape_i2
     MODULE PROCEDURE :: request_external_scratch_bounds_i2
     MODULE PROCEDURE :: request_external_scratch_shape_i3
     MODULE PROCEDURE :: request_external_scratch_bounds_i3
     MODULE PROCEDURE :: request_external_scratch_shape_i4
     MODULE PROCEDURE :: request_external_scratch_bounds_i4
     MODULE PROCEDURE :: request_external_scratch_shape_i5
     MODULE PROCEDURE :: request_external_scratch_bounds_i5
     MODULE PROCEDURE :: request_external_scratch_shape_i6
     MODULE PROCEDURE :: request_external_scratch_bounds_i6
     MODULE PROCEDURE :: request_internal_scratch_shape_i1
     MODULE PROCEDURE :: request_internal_scratch_bounds_i1
     MODULE PROCEDURE :: request_internal_scratch_shape_i2
     MODULE PROCEDURE :: request_internal_scratch_bounds_i2
     MODULE PROCEDURE :: request_internal_scratch_shape_i3
     MODULE PROCEDURE :: request_internal_scratch_bounds_i3
     MODULE PROCEDURE :: request_internal_scratch_shape_i4
     MODULE PROCEDURE :: request_internal_scratch_bounds_i4
     MODULE PROCEDURE :: request_internal_scratch_shape_i5
     MODULE PROCEDURE :: request_internal_scratch_bounds_i5
     MODULE PROCEDURE :: request_internal_scratch_shape_i6
     MODULE PROCEDURE :: request_internal_scratch_bounds_i6
  END INTERFACE request_scratch
  
  INTERFACE request_saved_scratch
     MODULE PROCEDURE :: request_saved_external_scratch_shape_dc1
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_dc1
     MODULE PROCEDURE :: request_saved_external_scratch_shape_dc2
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_dc2
     MODULE PROCEDURE :: request_saved_external_scratch_shape_dc3
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_dc3
     MODULE PROCEDURE :: request_saved_external_scratch_shape_dc4
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_dc4
     MODULE PROCEDURE :: request_saved_external_scratch_shape_dc5
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_dc5
     MODULE PROCEDURE :: request_saved_external_scratch_shape_dc6
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_dc6
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_dc1
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_dc1
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_dc2
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_dc2
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_dc3
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_dc3
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_dc4
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_dc4
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_dc5
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_dc5
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_dc6
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_dc6
     MODULE PROCEDURE :: request_saved_external_scratch_shape_rc1
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_rc1
     MODULE PROCEDURE :: request_saved_external_scratch_shape_rc2
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_rc2
     MODULE PROCEDURE :: request_saved_external_scratch_shape_rc3
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_rc3
     MODULE PROCEDURE :: request_saved_external_scratch_shape_rc4
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_rc4
     MODULE PROCEDURE :: request_saved_external_scratch_shape_rc5
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_rc5
     MODULE PROCEDURE :: request_saved_external_scratch_shape_rc6
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_rc6
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_rc1
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_rc1
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_rc2
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_rc2
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_rc3
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_rc3
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_rc4
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_rc4
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_rc5
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_rc5
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_rc6
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_rc6
     MODULE PROCEDURE :: request_saved_external_scratch_shape_d1
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_d1
     MODULE PROCEDURE :: request_saved_external_scratch_shape_d2
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_d2
     MODULE PROCEDURE :: request_saved_external_scratch_shape_d3
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_d3
     MODULE PROCEDURE :: request_saved_external_scratch_shape_d4
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_d4
     MODULE PROCEDURE :: request_saved_external_scratch_shape_d5
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_d5
     MODULE PROCEDURE :: request_saved_external_scratch_shape_d6
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_d6
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_d1
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_d1
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_d2
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_d2
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_d3
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_d3
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_d4
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_d4
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_d5
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_d5
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_d6
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_d6
     MODULE PROCEDURE :: request_saved_external_scratch_shape_r1
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_r1
     MODULE PROCEDURE :: request_saved_external_scratch_shape_r2
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_r2
     MODULE PROCEDURE :: request_saved_external_scratch_shape_r3
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_r3
     MODULE PROCEDURE :: request_saved_external_scratch_shape_r4
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_r4
     MODULE PROCEDURE :: request_saved_external_scratch_shape_r5
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_r5
     MODULE PROCEDURE :: request_saved_external_scratch_shape_r6
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_r6
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_r1
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_r1
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_r2
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_r2
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_r3
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_r3
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_r4
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_r4
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_r5
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_r5
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_r6
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_r6
     MODULE PROCEDURE :: request_saved_external_scratch_shape_il1
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_il1
     MODULE PROCEDURE :: request_saved_external_scratch_shape_il2
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_il2
     MODULE PROCEDURE :: request_saved_external_scratch_shape_il3
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_il3
     MODULE PROCEDURE :: request_saved_external_scratch_shape_il4
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_il4
     MODULE PROCEDURE :: request_saved_external_scratch_shape_il5
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_il5
     MODULE PROCEDURE :: request_saved_external_scratch_shape_il6
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_il6
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_il1
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_il1
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_il2
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_il2
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_il3
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_il3
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_il4
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_il4
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_il5
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_il5
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_il6
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_il6
     MODULE PROCEDURE :: request_saved_external_scratch_shape_i1
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_i1
     MODULE PROCEDURE :: request_saved_external_scratch_shape_i2
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_i2
     MODULE PROCEDURE :: request_saved_external_scratch_shape_i3
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_i3
     MODULE PROCEDURE :: request_saved_external_scratch_shape_i4
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_i4
     MODULE PROCEDURE :: request_saved_external_scratch_shape_i5
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_i5
     MODULE PROCEDURE :: request_saved_external_scratch_shape_i6
     MODULE PROCEDURE :: request_saved_external_scratch_bounds_i6
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_i1
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_i1
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_i2
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_i2
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_i3
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_i3
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_i4
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_i4
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_i5
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_i5
     MODULE PROCEDURE :: request_saved_internal_scratch_shape_i6
     MODULE PROCEDURE :: request_saved_internal_scratch_bounds_i6
  END INTERFACE request_saved_scratch

  INTERFACE free_scratch
     MODULE PROCEDURE :: free_external_scratch_shape_dc1
     MODULE PROCEDURE :: free_external_scratch_bounds_dc1
     MODULE PROCEDURE :: free_external_scratch_shape_dc2
     MODULE PROCEDURE :: free_external_scratch_bounds_dc2
     MODULE PROCEDURE :: free_external_scratch_shape_dc3
     MODULE PROCEDURE :: free_external_scratch_bounds_dc3
     MODULE PROCEDURE :: free_external_scratch_shape_dc4
     MODULE PROCEDURE :: free_external_scratch_bounds_dc4
     MODULE PROCEDURE :: free_external_scratch_shape_dc5
     MODULE PROCEDURE :: free_external_scratch_bounds_dc5
     MODULE PROCEDURE :: free_external_scratch_shape_dc6
     MODULE PROCEDURE :: free_external_scratch_bounds_dc6
     MODULE PROCEDURE :: free_internal_scratch_shape_dc1
     MODULE PROCEDURE :: free_internal_scratch_bounds_dc1
     MODULE PROCEDURE :: free_internal_scratch_shape_dc2
     MODULE PROCEDURE :: free_internal_scratch_bounds_dc2
     MODULE PROCEDURE :: free_internal_scratch_shape_dc3
     MODULE PROCEDURE :: free_internal_scratch_bounds_dc3
     MODULE PROCEDURE :: free_internal_scratch_shape_dc4
     MODULE PROCEDURE :: free_internal_scratch_bounds_dc4
     MODULE PROCEDURE :: free_internal_scratch_shape_dc5
     MODULE PROCEDURE :: free_internal_scratch_bounds_dc5
     MODULE PROCEDURE :: free_internal_scratch_shape_dc6
     MODULE PROCEDURE :: free_internal_scratch_bounds_dc6
     MODULE PROCEDURE :: free_external_scratch_shape_rc1
     MODULE PROCEDURE :: free_external_scratch_bounds_rc1
     MODULE PROCEDURE :: free_external_scratch_shape_rc2
     MODULE PROCEDURE :: free_external_scratch_bounds_rc2
     MODULE PROCEDURE :: free_external_scratch_shape_rc3
     MODULE PROCEDURE :: free_external_scratch_bounds_rc3
     MODULE PROCEDURE :: free_external_scratch_shape_rc4
     MODULE PROCEDURE :: free_external_scratch_bounds_rc4
     MODULE PROCEDURE :: free_external_scratch_shape_rc5
     MODULE PROCEDURE :: free_external_scratch_bounds_rc5
     MODULE PROCEDURE :: free_external_scratch_shape_rc6
     MODULE PROCEDURE :: free_external_scratch_bounds_rc6
     MODULE PROCEDURE :: free_internal_scratch_shape_rc1
     MODULE PROCEDURE :: free_internal_scratch_bounds_rc1
     MODULE PROCEDURE :: free_internal_scratch_shape_rc2
     MODULE PROCEDURE :: free_internal_scratch_bounds_rc2
     MODULE PROCEDURE :: free_internal_scratch_shape_rc3
     MODULE PROCEDURE :: free_internal_scratch_bounds_rc3
     MODULE PROCEDURE :: free_internal_scratch_shape_rc4
     MODULE PROCEDURE :: free_internal_scratch_bounds_rc4
     MODULE PROCEDURE :: free_internal_scratch_shape_rc5
     MODULE PROCEDURE :: free_internal_scratch_bounds_rc5
     MODULE PROCEDURE :: free_internal_scratch_shape_rc6
     MODULE PROCEDURE :: free_internal_scratch_bounds_rc6
     MODULE PROCEDURE :: free_external_scratch_shape_d1
     MODULE PROCEDURE :: free_external_scratch_bounds_d1
     MODULE PROCEDURE :: free_external_scratch_shape_d2
     MODULE PROCEDURE :: free_external_scratch_bounds_d2
     MODULE PROCEDURE :: free_external_scratch_shape_d3
     MODULE PROCEDURE :: free_external_scratch_bounds_d3
     MODULE PROCEDURE :: free_external_scratch_shape_d4
     MODULE PROCEDURE :: free_external_scratch_bounds_d4
     MODULE PROCEDURE :: free_external_scratch_shape_d5
     MODULE PROCEDURE :: free_external_scratch_bounds_d5
     MODULE PROCEDURE :: free_external_scratch_shape_d6
     MODULE PROCEDURE :: free_external_scratch_bounds_d6
     MODULE PROCEDURE :: free_internal_scratch_shape_d1
     MODULE PROCEDURE :: free_internal_scratch_bounds_d1
     MODULE PROCEDURE :: free_internal_scratch_shape_d2
     MODULE PROCEDURE :: free_internal_scratch_bounds_d2
     MODULE PROCEDURE :: free_internal_scratch_shape_d3
     MODULE PROCEDURE :: free_internal_scratch_bounds_d3
     MODULE PROCEDURE :: free_internal_scratch_shape_d4
     MODULE PROCEDURE :: free_internal_scratch_bounds_d4
     MODULE PROCEDURE :: free_internal_scratch_shape_d5
     MODULE PROCEDURE :: free_internal_scratch_bounds_d5
     MODULE PROCEDURE :: free_internal_scratch_shape_d6
     MODULE PROCEDURE :: free_internal_scratch_bounds_d6
     MODULE PROCEDURE :: free_external_scratch_shape_r1
     MODULE PROCEDURE :: free_external_scratch_bounds_r1
     MODULE PROCEDURE :: free_external_scratch_shape_r2
     MODULE PROCEDURE :: free_external_scratch_bounds_r2
     MODULE PROCEDURE :: free_external_scratch_shape_r3
     MODULE PROCEDURE :: free_external_scratch_bounds_r3
     MODULE PROCEDURE :: free_external_scratch_shape_r4
     MODULE PROCEDURE :: free_external_scratch_bounds_r4
     MODULE PROCEDURE :: free_external_scratch_shape_r5
     MODULE PROCEDURE :: free_external_scratch_bounds_r5
     MODULE PROCEDURE :: free_external_scratch_shape_r6
     MODULE PROCEDURE :: free_external_scratch_bounds_r6
     MODULE PROCEDURE :: free_internal_scratch_shape_r1
     MODULE PROCEDURE :: free_internal_scratch_bounds_r1
     MODULE PROCEDURE :: free_internal_scratch_shape_r2
     MODULE PROCEDURE :: free_internal_scratch_bounds_r2
     MODULE PROCEDURE :: free_internal_scratch_shape_r3
     MODULE PROCEDURE :: free_internal_scratch_bounds_r3
     MODULE PROCEDURE :: free_internal_scratch_shape_r4
     MODULE PROCEDURE :: free_internal_scratch_bounds_r4
     MODULE PROCEDURE :: free_internal_scratch_shape_r5
     MODULE PROCEDURE :: free_internal_scratch_bounds_r5
     MODULE PROCEDURE :: free_internal_scratch_shape_r6
     MODULE PROCEDURE :: free_internal_scratch_bounds_r6
     MODULE PROCEDURE :: free_external_scratch_shape_il1
     MODULE PROCEDURE :: free_external_scratch_bounds_il1
     MODULE PROCEDURE :: free_external_scratch_shape_il2
     MODULE PROCEDURE :: free_external_scratch_bounds_il2
     MODULE PROCEDURE :: free_external_scratch_shape_il3
     MODULE PROCEDURE :: free_external_scratch_bounds_il3
     MODULE PROCEDURE :: free_external_scratch_shape_il4
     MODULE PROCEDURE :: free_external_scratch_bounds_il4
     MODULE PROCEDURE :: free_external_scratch_shape_il5
     MODULE PROCEDURE :: free_external_scratch_bounds_il5
     MODULE PROCEDURE :: free_external_scratch_shape_il6
     MODULE PROCEDURE :: free_external_scratch_bounds_il6
     MODULE PROCEDURE :: free_internal_scratch_shape_il1
     MODULE PROCEDURE :: free_internal_scratch_bounds_il1
     MODULE PROCEDURE :: free_internal_scratch_shape_il2
     MODULE PROCEDURE :: free_internal_scratch_bounds_il2
     MODULE PROCEDURE :: free_internal_scratch_shape_il3
     MODULE PROCEDURE :: free_internal_scratch_bounds_il3
     MODULE PROCEDURE :: free_internal_scratch_shape_il4
     MODULE PROCEDURE :: free_internal_scratch_bounds_il4
     MODULE PROCEDURE :: free_internal_scratch_shape_il5
     MODULE PROCEDURE :: free_internal_scratch_bounds_il5
     MODULE PROCEDURE :: free_internal_scratch_shape_il6
     MODULE PROCEDURE :: free_internal_scratch_bounds_il6
     MODULE PROCEDURE :: free_external_scratch_shape_i1
     MODULE PROCEDURE :: free_external_scratch_bounds_i1
     MODULE PROCEDURE :: free_external_scratch_shape_i2
     MODULE PROCEDURE :: free_external_scratch_bounds_i2
     MODULE PROCEDURE :: free_external_scratch_shape_i3
     MODULE PROCEDURE :: free_external_scratch_bounds_i3
     MODULE PROCEDURE :: free_external_scratch_shape_i4
     MODULE PROCEDURE :: free_external_scratch_bounds_i4
     MODULE PROCEDURE :: free_external_scratch_shape_i5
     MODULE PROCEDURE :: free_external_scratch_bounds_i5
     MODULE PROCEDURE :: free_external_scratch_shape_i6
     MODULE PROCEDURE :: free_external_scratch_bounds_i6
     MODULE PROCEDURE :: free_internal_scratch_shape_i1
     MODULE PROCEDURE :: free_internal_scratch_bounds_i1
     MODULE PROCEDURE :: free_internal_scratch_shape_i2
     MODULE PROCEDURE :: free_internal_scratch_bounds_i2
     MODULE PROCEDURE :: free_internal_scratch_shape_i3
     MODULE PROCEDURE :: free_internal_scratch_bounds_i3
     MODULE PROCEDURE :: free_internal_scratch_shape_i4
     MODULE PROCEDURE :: free_internal_scratch_bounds_i4
     MODULE PROCEDURE :: free_internal_scratch_shape_i5
     MODULE PROCEDURE :: free_internal_scratch_bounds_i5
     MODULE PROCEDURE :: free_internal_scratch_shape_i6
     MODULE PROCEDURE :: free_internal_scratch_bounds_i6
  END INTERFACE free_scratch

  INTERFACE save_scratch
     MODULE PROCEDURE :: save_external_scratch_shape_dc1
     MODULE PROCEDURE :: save_external_scratch_bounds_dc1
     MODULE PROCEDURE :: save_external_scratch_shape_dc2
     MODULE PROCEDURE :: save_external_scratch_bounds_dc2
     MODULE PROCEDURE :: save_external_scratch_shape_dc3
     MODULE PROCEDURE :: save_external_scratch_bounds_dc3
     MODULE PROCEDURE :: save_external_scratch_shape_dc4
     MODULE PROCEDURE :: save_external_scratch_bounds_dc4
     MODULE PROCEDURE :: save_external_scratch_shape_dc5
     MODULE PROCEDURE :: save_external_scratch_bounds_dc5
     MODULE PROCEDURE :: save_external_scratch_shape_dc6
     MODULE PROCEDURE :: save_external_scratch_bounds_dc6
     MODULE PROCEDURE :: save_internal_scratch_shape_dc1
     MODULE PROCEDURE :: save_internal_scratch_bounds_dc1
     MODULE PROCEDURE :: save_internal_scratch_shape_dc2
     MODULE PROCEDURE :: save_internal_scratch_bounds_dc2
     MODULE PROCEDURE :: save_internal_scratch_shape_dc3
     MODULE PROCEDURE :: save_internal_scratch_bounds_dc3
     MODULE PROCEDURE :: save_internal_scratch_shape_dc4
     MODULE PROCEDURE :: save_internal_scratch_bounds_dc4
     MODULE PROCEDURE :: save_internal_scratch_shape_dc5
     MODULE PROCEDURE :: save_internal_scratch_bounds_dc5
     MODULE PROCEDURE :: save_internal_scratch_shape_dc6
     MODULE PROCEDURE :: save_internal_scratch_bounds_dc6
     MODULE PROCEDURE :: save_external_scratch_shape_rc1
     MODULE PROCEDURE :: save_external_scratch_bounds_rc1
     MODULE PROCEDURE :: save_external_scratch_shape_rc2
     MODULE PROCEDURE :: save_external_scratch_bounds_rc2
     MODULE PROCEDURE :: save_external_scratch_shape_rc3
     MODULE PROCEDURE :: save_external_scratch_bounds_rc3
     MODULE PROCEDURE :: save_external_scratch_shape_rc4
     MODULE PROCEDURE :: save_external_scratch_bounds_rc4
     MODULE PROCEDURE :: save_external_scratch_shape_rc5
     MODULE PROCEDURE :: save_external_scratch_bounds_rc5
     MODULE PROCEDURE :: save_external_scratch_shape_rc6
     MODULE PROCEDURE :: save_external_scratch_bounds_rc6
     MODULE PROCEDURE :: save_internal_scratch_shape_rc1
     MODULE PROCEDURE :: save_internal_scratch_bounds_rc1
     MODULE PROCEDURE :: save_internal_scratch_shape_rc2
     MODULE PROCEDURE :: save_internal_scratch_bounds_rc2
     MODULE PROCEDURE :: save_internal_scratch_shape_rc3
     MODULE PROCEDURE :: save_internal_scratch_bounds_rc3
     MODULE PROCEDURE :: save_internal_scratch_shape_rc4
     MODULE PROCEDURE :: save_internal_scratch_bounds_rc4
     MODULE PROCEDURE :: save_internal_scratch_shape_rc5
     MODULE PROCEDURE :: save_internal_scratch_bounds_rc5
     MODULE PROCEDURE :: save_internal_scratch_shape_rc6
     MODULE PROCEDURE :: save_internal_scratch_bounds_rc6
     MODULE PROCEDURE :: save_external_scratch_shape_d1
     MODULE PROCEDURE :: save_external_scratch_bounds_d1
     MODULE PROCEDURE :: save_external_scratch_shape_d2
     MODULE PROCEDURE :: save_external_scratch_bounds_d2
     MODULE PROCEDURE :: save_external_scratch_shape_d3
     MODULE PROCEDURE :: save_external_scratch_bounds_d3
     MODULE PROCEDURE :: save_external_scratch_shape_d4
     MODULE PROCEDURE :: save_external_scratch_bounds_d4
     MODULE PROCEDURE :: save_external_scratch_shape_d5
     MODULE PROCEDURE :: save_external_scratch_bounds_d5
     MODULE PROCEDURE :: save_external_scratch_shape_d6
     MODULE PROCEDURE :: save_external_scratch_bounds_d6
     MODULE PROCEDURE :: save_internal_scratch_shape_d1
     MODULE PROCEDURE :: save_internal_scratch_bounds_d1
     MODULE PROCEDURE :: save_internal_scratch_shape_d2
     MODULE PROCEDURE :: save_internal_scratch_bounds_d2
     MODULE PROCEDURE :: save_internal_scratch_shape_d3
     MODULE PROCEDURE :: save_internal_scratch_bounds_d3
     MODULE PROCEDURE :: save_internal_scratch_shape_d4
     MODULE PROCEDURE :: save_internal_scratch_bounds_d4
     MODULE PROCEDURE :: save_internal_scratch_shape_d5
     MODULE PROCEDURE :: save_internal_scratch_bounds_d5
     MODULE PROCEDURE :: save_internal_scratch_shape_d6
     MODULE PROCEDURE :: save_internal_scratch_bounds_d6
     MODULE PROCEDURE :: save_external_scratch_shape_r1
     MODULE PROCEDURE :: save_external_scratch_bounds_r1
     MODULE PROCEDURE :: save_external_scratch_shape_r2
     MODULE PROCEDURE :: save_external_scratch_bounds_r2
     MODULE PROCEDURE :: save_external_scratch_shape_r3
     MODULE PROCEDURE :: save_external_scratch_bounds_r3
     MODULE PROCEDURE :: save_external_scratch_shape_r4
     MODULE PROCEDURE :: save_external_scratch_bounds_r4
     MODULE PROCEDURE :: save_external_scratch_shape_r5
     MODULE PROCEDURE :: save_external_scratch_bounds_r5
     MODULE PROCEDURE :: save_external_scratch_shape_r6
     MODULE PROCEDURE :: save_external_scratch_bounds_r6
     MODULE PROCEDURE :: save_internal_scratch_shape_r1
     MODULE PROCEDURE :: save_internal_scratch_bounds_r1
     MODULE PROCEDURE :: save_internal_scratch_shape_r2
     MODULE PROCEDURE :: save_internal_scratch_bounds_r2
     MODULE PROCEDURE :: save_internal_scratch_shape_r3
     MODULE PROCEDURE :: save_internal_scratch_bounds_r3
     MODULE PROCEDURE :: save_internal_scratch_shape_r4
     MODULE PROCEDURE :: save_internal_scratch_bounds_r4
     MODULE PROCEDURE :: save_internal_scratch_shape_r5
     MODULE PROCEDURE :: save_internal_scratch_bounds_r5
     MODULE PROCEDURE :: save_internal_scratch_shape_r6
     MODULE PROCEDURE :: save_internal_scratch_bounds_r6
     MODULE PROCEDURE :: save_external_scratch_shape_il1
     MODULE PROCEDURE :: save_external_scratch_bounds_il1
     MODULE PROCEDURE :: save_external_scratch_shape_il2
     MODULE PROCEDURE :: save_external_scratch_bounds_il2
     MODULE PROCEDURE :: save_external_scratch_shape_il3
     MODULE PROCEDURE :: save_external_scratch_bounds_il3
     MODULE PROCEDURE :: save_external_scratch_shape_il4
     MODULE PROCEDURE :: save_external_scratch_bounds_il4
     MODULE PROCEDURE :: save_external_scratch_shape_il5
     MODULE PROCEDURE :: save_external_scratch_bounds_il5
     MODULE PROCEDURE :: save_external_scratch_shape_il6
     MODULE PROCEDURE :: save_external_scratch_bounds_il6
     MODULE PROCEDURE :: save_internal_scratch_shape_il1
     MODULE PROCEDURE :: save_internal_scratch_bounds_il1
     MODULE PROCEDURE :: save_internal_scratch_shape_il2
     MODULE PROCEDURE :: save_internal_scratch_bounds_il2
     MODULE PROCEDURE :: save_internal_scratch_shape_il3
     MODULE PROCEDURE :: save_internal_scratch_bounds_il3
     MODULE PROCEDURE :: save_internal_scratch_shape_il4
     MODULE PROCEDURE :: save_internal_scratch_bounds_il4
     MODULE PROCEDURE :: save_internal_scratch_shape_il5
     MODULE PROCEDURE :: save_internal_scratch_bounds_il5
     MODULE PROCEDURE :: save_internal_scratch_shape_il6
     MODULE PROCEDURE :: save_internal_scratch_bounds_il6
     MODULE PROCEDURE :: save_external_scratch_shape_i1
     MODULE PROCEDURE :: save_external_scratch_bounds_i1
     MODULE PROCEDURE :: save_external_scratch_shape_i2
     MODULE PROCEDURE :: save_external_scratch_bounds_i2
     MODULE PROCEDURE :: save_external_scratch_shape_i3
     MODULE PROCEDURE :: save_external_scratch_bounds_i3
     MODULE PROCEDURE :: save_external_scratch_shape_i4
     MODULE PROCEDURE :: save_external_scratch_bounds_i4
     MODULE PROCEDURE :: save_external_scratch_shape_i5
     MODULE PROCEDURE :: save_external_scratch_bounds_i5
     MODULE PROCEDURE :: save_external_scratch_shape_i6
     MODULE PROCEDURE :: save_external_scratch_bounds_i6
     MODULE PROCEDURE :: save_internal_scratch_shape_i1
     MODULE PROCEDURE :: save_internal_scratch_bounds_i1
     MODULE PROCEDURE :: save_internal_scratch_shape_i2
     MODULE PROCEDURE :: save_internal_scratch_bounds_i2
     MODULE PROCEDURE :: save_internal_scratch_shape_i3
     MODULE PROCEDURE :: save_internal_scratch_bounds_i3
     MODULE PROCEDURE :: save_internal_scratch_shape_i4
     MODULE PROCEDURE :: save_internal_scratch_bounds_i4
     MODULE PROCEDURE :: save_internal_scratch_shape_i5
     MODULE PROCEDURE :: save_internal_scratch_bounds_i5
     MODULE PROCEDURE :: save_internal_scratch_shape_i6
     MODULE PROCEDURE :: save_internal_scratch_bounds_i6
  END INTERFACE save_scratch

CONTAINS

  SUBROUTINE init_internal_pool( initial_len, ierr )
    INTEGER( INT64 ), INTENT( IN )  :: initial_len
    INTEGER( INT32 ), INTENT( OUT ) :: ierr

    CALL init_pool( pool_int, initial_len, ierr )
    
  END SUBROUTINE init_internal_pool

  SUBROUTINE init_external_pool( pool_ext, initial_len, ierr )
    TYPE( memory_pool ), INTENT( INOUT ) :: pool_ext
    INTEGER( INT64 ), INTENT( IN )       :: initial_len
    INTEGER( INT32 ), INTENT( OUT )      :: ierr

    CALL init_pool( pool_ext, initial_len, ierr )
    
  END SUBROUTINE init_external_pool

  SUBROUTINE finalize_internal_pool( finalized_pool_len, ierr )
    INTEGER( INT64 ), INTENT( OUT ) :: finalized_pool_len
    INTEGER( INT32 ), INTENT( OUT ) :: ierr

    CALL finalize_pool( pool_int, finalized_pool_len, ierr )
    
  END SUBROUTINE finalize_internal_pool

  SUBROUTINE finalize_external_pool( pool_ext, finalized_pool_len, ierr )
    TYPE( memory_pool ), INTENT( INOUT ) :: pool_ext
    INTEGER( INT64 ), INTENT( OUT )      :: finalized_pool_len
    INTEGER( INT32 ), INTENT( OUT )      :: ierr

    CALL finalize_pool( pool_ext, finalized_pool_len, ierr )
    
  END SUBROUTINE finalize_external_pool
  
#include "request_scratch.inc"
  
END MODULE scratch_interface
