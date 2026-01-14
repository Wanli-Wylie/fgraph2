MODULE mixed_mod
  IMPLICIT NONE
  INTEGER :: counter
CONTAINS
  SUBROUTINE bump()
    counter = counter + 1
  END SUBROUTINE bump
END MODULE mixed_mod

PROGRAM mixed_main
  USE mixed_mod
  IMPLICIT NONE
  counter = 0
  CALL bump()
  CALL bump()
  PRINT *, counter
CONTAINS
  SUBROUTINE reset()
    counter = 0
  END SUBROUTINE reset
END PROGRAM mixed_main
