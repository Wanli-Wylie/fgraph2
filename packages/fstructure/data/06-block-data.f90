MODULE common_block_mod
  IMPLICIT NONE
  INTEGER :: seed
  REAL :: weights(3)
  COMMON /shared_block/ seed, weights
END MODULE common_block_mod

BLOCK DATA shared_init
  IMPLICIT NONE
  INTEGER :: seed
  REAL :: weights(3)
  COMMON /shared_block/ seed, weights
  DATA seed /42/
  DATA weights /1.0, 2.0, 3.0/
END BLOCK DATA shared_init

PROGRAM block_data_driver
  USE common_block_mod
  IMPLICIT NONE
  PRINT *, seed, weights
END PROGRAM block_data_driver
