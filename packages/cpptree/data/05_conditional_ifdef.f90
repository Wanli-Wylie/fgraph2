! Example: #ifdef and #ifndef conditional compilation
! This file demonstrates macro existence checking

#define DEBUG
#define FEATURE_X

#ifdef DEBUG
    ! Debug mode code
    subroutine debug_print(message)
        character(len=*), intent(in) :: message
        print *, "[DEBUG] ", message
    end subroutine
#endif

#ifndef RELEASE
    ! Non-release code
    subroutine development_code()
        print *, "Development code active"
    end subroutine
#endif

#ifdef FEATURE_X
    subroutine feature_x()
        print *, "Feature X enabled"
    end subroutine
#else
    subroutine feature_x()
        print *, "Feature X disabled"
    end subroutine
#endif

program conditional_ifdef
    implicit none
    
#ifdef DEBUG
    call debug_print("Program starting")
#endif

#ifndef RELEASE
    call development_code()
#endif

    call feature_x()
end program conditional_ifdef

