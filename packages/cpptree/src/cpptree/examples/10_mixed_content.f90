! Example: Mixed content with text blocks and directives
! This file demonstrates how preprocessor directives mix with regular code

program mixed_content
    implicit none
    
    ! Regular Fortran code
    integer :: i, j, k
    real(8) :: x, y, z
    
#define USE_OPTIMIZATION
#ifdef USE_OPTIMIZATION
    ! Optimized code path
    do i = 1, 100
        x = x + 1.0
    end do
#else
    ! Standard code path
    do i = 1, 100
        x = x + 1.0
        y = y + 0.5
    end do
#endif

    ! More regular code
    z = x + y
    
#define DEBUG_MODE
#ifdef DEBUG_MODE
    print *, "Debug: x = ", x, ", y = ", y, ", z = ", z
#endif

    ! Final code block
    print *, "Program completed"
end program mixed_content

