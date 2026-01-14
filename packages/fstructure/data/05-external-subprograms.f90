PROGRAM external_driver
  IMPLICIT NONE
  INTEGER :: value
  value = double_it(7)
  CALL show_value(value)
END PROGRAM external_driver

INTEGER FUNCTION double_it(x)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: x
  double_it = 2 * x
END FUNCTION double_it

SUBROUTINE show_value(x)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: x
  PRINT *, "value =", x
CONTAINS
  SUBROUTINE show_detail(y)
    INTEGER, INTENT(IN) :: y
    PRINT *, "detail:", y
  END SUBROUTINE show_detail
END SUBROUTINE show_value
