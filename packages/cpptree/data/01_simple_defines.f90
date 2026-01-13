! Example: Simple #define macros
! This file demonstrates basic object-like macro definitions

#define PI 3.141592653589793
#define E 2.718281828459045
#define MAX_ITERATIONS 1000
#define DEBUG_MODE 1
#define VERSION "1.0.0"

program simple_defines
    implicit none
    real(8) :: radius, area
    integer :: i
    
    radius = 5.0
    area = PI * radius * radius
    
    do i = 1, MAX_ITERATIONS
        ! Some computation
    end do
    
    print *, "Version: ", VERSION
end program simple_defines

