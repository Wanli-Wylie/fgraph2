! Example: #if conditional compilation
! This file demonstrates #if, #elif, #else, #endif

#define PLATFORM 1

#if PLATFORM == 1
    ! Platform 1 specific code
    subroutine platform1_init()
        print *, "Initializing Platform 1"
    end subroutine
#elif PLATFORM == 2
    ! Platform 2 specific code
    subroutine platform2_init()
        print *, "Initializing Platform 2"
    end subroutine
#else
    ! Default platform code
    subroutine default_init()
        print *, "Initializing Default Platform"
    end subroutine
#endif

program conditional_if
    implicit none
#if PLATFORM == 1
    call platform1_init()
#elif PLATFORM == 2
    call platform2_init()
#else
    call default_init()
#endif
end program conditional_if

