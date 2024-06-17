module gpu
  IMPLICIT NONE

  LOGICAL :: update_first_to_gpu   =.TRUE.
  LOGICAL :: update_second_to_gpu  =.TRUE.
  LOGICAL :: update_third_to_gpu   =.TRUE.
  LOGICAL :: update_result_to_host =.TRUE.
  LOGICAL :: comm_buffers_on_host  =.TRUE.
end module gpu
