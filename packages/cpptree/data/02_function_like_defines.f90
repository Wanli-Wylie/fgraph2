! Example: Function-like #define macros
! This file demonstrates function-like macro definitions

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define SQUARE(x) ((x) * (x))
#define ABS(x) ((x) < 0 ? -(x) : (x))
#define SWAP(a, b) do { \
    typeof(a) _temp = (a); \
    (a) = (b); \
    (b) = _temp; \
} while(0)

program function_like_defines
    implicit none
    integer :: a, b, result
    
    a = 10
    b = 20
    
    result = MAX(a, b)
    result = MIN(a, b)
    result = SQUARE(5)
end program function_like_defines

