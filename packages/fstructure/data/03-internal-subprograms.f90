PROGRAM internal_demo
  IMPLICIT NONE
  INTEGER :: total

  total = add(3, 4)
  CALL report(total)

CONTAINS
  FUNCTION add(a, b) RESULT(sum)
    INTEGER, INTENT(IN) :: a, b
    INTEGER :: sum
    sum = a + b
  END FUNCTION add

  SUBROUTINE report(value)
    INTEGER, INTENT(IN) :: value
    PRINT *, "total =", value
  END SUBROUTINE report
END PROGRAM internal_demo
