! Example: Multiple #elif branches
! This file demonstrates conditional compilation with many branches

#define COMPILER_TYPE 3

#if COMPILER_TYPE == 1
    ! GCC compiler
    subroutine gcc_specific()
        print *, "GCC compiler detected"
    end subroutine
#elif COMPILER_TYPE == 2
    ! Intel compiler
    subroutine intel_specific()
        print *, "Intel compiler detected"
    end subroutine
#elif COMPILER_TYPE == 3
    ! PGI compiler
    subroutine pgi_specific()
        print *, "PGI compiler detected"
    end subroutine
#elif COMPILER_TYPE == 4
    ! Cray compiler
    subroutine cray_specific()
        print *, "Cray compiler detected"
    end subroutine
#else
    ! Unknown compiler
    subroutine unknown_compiler()
        print *, "Unknown compiler"
    end subroutine
#endif

program multiple_elif
    implicit none
    
#if COMPILER_TYPE == 1
    call gcc_specific()
#elif COMPILER_TYPE == 2
    call intel_specific()
#elif COMPILER_TYPE == 3
    call pgi_specific()
#elif COMPILER_TYPE == 4
    call cray_specific()
#else
    call unknown_compiler()
#endif
end program multiple_elif

