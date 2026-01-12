"""Examples module for cpptree.

This module provides access to example Fortran files with C preprocessor directives,
along with their tree-sitter parsed AST representations.
"""

from __future__ import annotations
from pathlib import Path
from typing import Optional
from tree_sitter import Language, Parser, Tree
import tree_sitter_fortran

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


def _get_examples_dir() -> Path:
    """获取 examples 目录的路径。"""
    return Path(__file__).parent


def list_examples() -> list[str]:
    """列出所有可用的示例文件名。
    
    Returns:
        示例文件名列表（不包括路径，只包含文件名）。
    """
    examples_dir = _get_examples_dir()
    # 查找所有 .f90 文件
    example_files = sorted(examples_dir.glob("*.f90"))
    return [f.name for f in example_files]


def get_example(name: str) -> tuple[str, Tree]:
    """获取指定示例文件的源代码和解析后的 AST。
    
    Args:
        name: 示例文件名（例如 "01_simple_defines.f90"）。如果只提供文件名
              不含扩展名，会自动添加 .f90。
    
    Returns:
        一个元组 (source_code, tree)，其中：
        - source_code: 文件的源代码内容（字符串）
        - tree: tree-sitter 解析后的 Tree 对象
    
    Raises:
        FileNotFoundError: 如果指定的示例文件不存在。
        ValueError: 如果文件名无效。
    """
    examples_dir = _get_examples_dir()
    
    # 确保文件名有 .f90 扩展名
    if not name.endswith('.f90'):
        name = f"{name}.f90"
    
    example_path = examples_dir / name
    
    if not example_path.exists():
        available = list_examples()
        raise FileNotFoundError(
            f"示例文件 '{name}' 不存在。\n"
            f"可用的示例文件: {', '.join(available)}"
        )
    
    # 读取源代码
    source_code = example_path.read_text(encoding='utf-8')
    
    # 使用 tree-sitter 解析
    parser = _get_parser()
    tree = parser.parse(bytes(source_code, 'utf-8'))
    
    return source_code, tree


def get_example_path(name: str) -> Path:
    """获取示例文件的路径。
    
    Args:
        name: 示例文件名（例如 "01_simple_defines.f90"）。
    
    Returns:
        示例文件的 Path 对象。
    
    Raises:
        FileNotFoundError: 如果指定的示例文件不存在。
    """
    examples_dir = _get_examples_dir()
    
    # 确保文件名有 .f90 扩展名
    if not name.endswith('.f90'):
        name = f"{name}.f90"
    
    example_path = examples_dir / name
    
    if not example_path.exists():
        available = list_examples()
        raise FileNotFoundError(
            f"示例文件 '{name}' 不存在。\n"
            f"可用的示例文件: {', '.join(available)}"
        )
    
    return example_path


__all__ = ['list_examples', 'get_example', 'get_example_path']

