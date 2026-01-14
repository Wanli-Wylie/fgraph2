MODULE math_kinds
  IMPLICIT NONE
  INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(12)
  INTEGER, PARAMETER :: ik = SELECTED_INT_KIND(9)
END MODULE math_kinds

MODULE use_examples
  USE math_kinds, ONLY: rk, ik
  USE math_kinds, ONLY: real_kind => rk
  USE math_kinds, ONLY:
  IMPLICIT NONE
  REAL(rk) :: scale
CONTAINS
  SUBROUTINE init_scale(value)
    REAL(rk), INTENT(IN) :: value
    scale = value
  END SUBROUTINE init_scale
END MODULE use_examples

PROGRAM use_driver
  USE use_examples, ONLY: init_scale
  IMPLICIT NONE
  CALL init_scale(2.5)
END PROGRAM use_driver
