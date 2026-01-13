! Example: #undef directive
! This file demonstrates macro undefinition

#define TEMP_MACRO 42
#define OLD_API

program undef_example
    implicit none
    
    ! Use TEMP_MACRO here
    print *, "Using TEMP_MACRO"
    
    ! Undefine it
#undef TEMP_MACRO

    ! TEMP_MACRO is no longer defined here
    
    ! Undefine OLD_API
#undef OLD_API

#ifdef OLD_API
    ! This code won't be compiled
    print *, "Old API"
#endif
end program undef_example

