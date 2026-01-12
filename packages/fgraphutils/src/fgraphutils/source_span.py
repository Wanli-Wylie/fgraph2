"""
Author: Wanli-Wylie Zhang
Date: 2026-01-12
Source span representation for Fortran source files.

This module provides a data structure to represent a contiguous region of source code
in a Fortran file, identified by file path and line/column positions.
"""

from pydantic import BaseModel

class SourceSpan(BaseModel):
    """
    Represents a span of source code in a Fortran source file.
    
    A source span defines a contiguous region of code using file path and
    start/end positions specified by line and column numbers. This is useful
    for tracking code locations, error reporting, and code analysis tools
    that need to reference specific portions of Fortran source files.
    
    Attributes:
        filepath: Path to the Fortran source file (relative or absolute).
        start_line: One-based line number where the span begins.
        start_column: One-based column number where the span begins.
        end_line: One-based line number where the span ends (inclusive).
        end_column: One-based column number where the span ends (inclusive).
    
    Example:
        >>> span = SourceSpan(
        ...     filepath="src/main.f90",
        ...     start_line=10,
        ...     start_column=5,
        ...     end_line=10,
        ...     end_column=20
        ... )
        >>> print(span)
        src/main.f90:10:5-10:20
    """
    filepath: str
    start_line: int
    start_column: int
    end_line: int
    end_column: int

    def __str__(self) -> str:
        return f"{self.filepath}:{self.start_line}:{self.start_column}-{self.end_line}:{self.end_column}"
    
    def __repr__(self) -> str:
        return f"SourceSpan(filepath={self.filepath}, start_line={self.start_line}, start_column={self.start_column}, end_line={self.end_line}, end_column={self.end_column})"