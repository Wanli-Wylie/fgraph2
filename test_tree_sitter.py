from tree_sitter import Language, Parser
import tree_sitter_fortran
import tree_sitter_c

# 初始化
FORTRAN_LANGUAGE = Language(tree_sitter_fortran.language())
parser = Parser(FORTRAN_LANGUAGE)

C_LANGUAGE = Language(tree_sitter_c.language())
parser = Parser(C_LANGUAGE)
# 你的混合代码
source_code = b"""
MODULE my_mod
#include "defs.h"
! #ifdef FAKE_NEWS
CONTAINS
  SUBROUTINE foo()
#ifdef MPI
    CALL mpi_init()
#else
    PRINT *, "Serial"
#endif
  END SUBROUTINE
END MODULE
"""

tree = parser.parse(source_code)
root = tree.root_node

# 遍历查看结构
def print_tree(node, depth=0):
    indent = "  " * depth
    # Tree-sitter 会把 #ifdef 识别为 preproc_if 之类的节点
    print(f"{indent}{node.type} [{node.start_point}-{node.end_point}]")
    for child in node.children:
        print_tree(child, depth + 1)

print_tree(root)