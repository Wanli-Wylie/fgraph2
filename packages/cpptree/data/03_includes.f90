! Example: #include directives
! This file demonstrates various include patterns

#include "config.h"
#include "constants.f90"
#include <system_header.h>

program includes_example
    implicit none
    
    ! Code that uses definitions from included files
    print *, "Program using includes"
end program includes_example

