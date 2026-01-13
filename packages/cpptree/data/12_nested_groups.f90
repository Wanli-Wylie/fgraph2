! Example: Nested conditional groups
! This file demonstrates nested #if/#endif blocks

#define LEVEL1_ENABLED
#define LEVEL2_ENABLED
#define LEVEL3_ENABLED

program nested_groups
    implicit none
    
#ifdef LEVEL1_ENABLED
    print *, "Level 1 enabled"
    
    #ifdef LEVEL2_ENABLED
        print *, "Level 2 enabled"
        
        #ifdef LEVEL3_ENABLED
            print *, "Level 3 enabled"
        #else
            print *, "Level 3 disabled"
        #endif
        
    #else
        print *, "Level 2 disabled"
    #endif
    
#else
    print *, "Level 1 disabled"
#endif

end program nested_groups

