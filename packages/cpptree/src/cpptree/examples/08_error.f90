! Example: #error directive
! This file demonstrates error generation

#define REQUIRED_VERSION 3

#if REQUIRED_VERSION < 2
#error "Version 2 or higher is required"
#endif

#ifndef COMPILER_SUPPORTED
#error "This compiler is not supported. Please use a supported compiler."
#endif

program error_example
    implicit none
    print *, "Program continues if no errors"
end program error_example

