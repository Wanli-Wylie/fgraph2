! Example: Define/undefine cycles
! This file demonstrates redefining macros

#define MACRO_VALUE 10

program define_undef_cycle
    implicit none
    
    ! First use
    print *, "First MACRO_VALUE usage"
    
#undef MACRO_VALUE
#define MACRO_VALUE 20

    ! Second use with new value
    print *, "Second MACRO_VALUE usage"
    
#undef MACRO_VALUE
#define MACRO_VALUE 30

    ! Third use with another value
    print *, "Third MACRO_VALUE usage"
    
#undef MACRO_VALUE

    ! MACRO_VALUE is now undefined
#ifdef MACRO_VALUE
    print *, "This won't print"
#endif
end program define_undef_cycle

