MODULE geometry
  IMPLICIT NONE
  REAL, PARAMETER :: pi = 3.14159
CONTAINS
  SUBROUTINE circle_area(r, area)
    REAL, INTENT(IN) :: r
    REAL, INTENT(OUT) :: area
    area = pi * r * r
  END SUBROUTINE circle_area

  FUNCTION annulus_area(r1, r2) RESULT(area)
    REAL, INTENT(IN) :: r1, r2
    REAL :: area
    area = ring(r1, r2)
  CONTAINS
    FUNCTION ring(a, b) RESULT(value)
      REAL, INTENT(IN) :: a, b
      REAL :: value
      value = pi * (b*b - a*a)
    END FUNCTION ring
  END FUNCTION annulus_area
END MODULE geometry

PROGRAM module_driver
  USE geometry, ONLY: circle_area, annulus_area
  IMPLICIT NONE
  REAL :: a1, a2
  CALL circle_area(2.0, a1)
  a2 = annulus_area(1.0, 3.0)
  PRINT *, a1, a2
END PROGRAM module_driver
