! Example: Complex conditional compilation
! This file demonstrates nested and complex conditional structures

#define ARCH_X86
#define OS_LINUX
#define FEATURE_A
#define FEATURE_B

#ifdef ARCH_X86
    #ifdef OS_LINUX
        subroutine linux_x86_init()
            print *, "Linux x86 initialization"
        end subroutine
    #elif defined(OS_WINDOWS)
        subroutine windows_x86_init()
            print *, "Windows x86 initialization"
        end subroutine
    #endif
#endif

#if defined(FEATURE_A) && defined(FEATURE_B)
    subroutine feature_ab()
        print *, "Both features A and B enabled"
    end subroutine
#elif defined(FEATURE_A)
    subroutine feature_a_only()
        print *, "Only feature A enabled"
    end subroutine
#elif defined(FEATURE_B)
    subroutine feature_b_only()
        print *, "Only feature B enabled"
    end subroutine
#else
    subroutine no_features()
        print *, "No features enabled"
    end subroutine
#endif

program complex_conditionals
    implicit none
    
#ifdef ARCH_X86
    #ifdef OS_LINUX
        call linux_x86_init()
    #endif
#endif

#if defined(FEATURE_A) && defined(FEATURE_B)
    call feature_ab()
#elif defined(FEATURE_A)
    call feature_a_only()
#elif defined(FEATURE_B)
    call feature_b_only()
#else
    call no_features()
#endif
end program complex_conditionals

