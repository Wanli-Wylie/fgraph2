! Example: #pragma directives
! This file demonstrates compiler-specific pragmas

#pragma once
#pragma pack(1)
#pragma GCC optimize("O3")
#pragma omp parallel
#pragma ivdep

program pragma_example
    implicit none
    integer :: i
    
    ! Some computation
    do i = 1, 100
        ! Loop body
    end do
end program pragma_example

