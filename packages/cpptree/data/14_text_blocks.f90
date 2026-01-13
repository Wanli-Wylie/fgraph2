! Example: Text blocks between directives
! This file demonstrates how text blocks appear in the AST

program text_blocks
    implicit none
    
    ! This is a text block before any directive
    integer :: var1, var2
    
#define DIRECTIVE_1
    ! Text block after directive 1
    var1 = 10
    
#define DIRECTIVE_2
    ! Another text block
    var2 = 20
    
    ! Text block with code
    print *, var1, var2
    
#undef DIRECTIVE_1
    ! Text block after undef
    var1 = var1 + 1
    
    ! Final text block
    print *, "Done"
end program text_blocks

