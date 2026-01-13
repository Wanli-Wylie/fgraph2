"""Examples module for cpptree.

This module provides access to example Fortran files with C preprocessor directives,
along with their tree-sitter parsed AST representations.
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional
from tree_sitter import Language, Parser, Tree
import tree_sitter_fortran
from fgraphutils import DomainDataLoader

# 初始化 tree-sitter-fortran parser（单例模式）
_parser: Optional[Parser] = None
_language: Optional[Language] = None


def _get_parser() -> Parser:
    """获取或初始化 tree-sitter Fortran 解析器（单例）。"""
    global _parser, _language
    
    if _parser is None:
        _language = Language(tree_sitter_fortran.language())
        _parser = Parser(_language)
    
    return _parser


class TreeSitterFortranDataLoader(DomainDataLoader[Tree]):
    """数据加载器，用于加载和解析 Fortran 示例文件。"""
    
    @property
    def file_extension(self) -> str:
        """Fortran 文件扩展名。"""
        return '.f90'
    
    @property
    def description(self) -> str:
        """描述这些测试数据的用途。"""
        return (
            "Fortran 文件示例，包含 C 预处理器指令，"
            "用于测试 cpptree 的解析和评估功能。"
        )
    
    def parse(self, content: str, case_id: str) -> Tree:
        """
        解析 Fortran 源代码并返回 AST 树。
        
        Args:
            content: 源代码内容
            case_id: 测试用例 ID（文件名，不含扩展名）
        
        Returns:
            tree-sitter 解析后的 Tree 对象
        """
        parser = _get_parser()
        tree = parser.parse(bytes(content, 'utf-8'))
        return tree


# 初始化数据加载器
# 使用相对于包根目录的 data 目录
_package_root = Path(__file__).parent.parent.parent.parent
_data_dir = _package_root / "data"
data_loader = TreeSitterFortranDataLoader(_data_dir)


__all__ = ['TreeSitterFortranDataLoader', 'data_loader']

