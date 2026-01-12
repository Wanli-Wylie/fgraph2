from __future__ import annotations
from typing import Literal, Union, Optional, Annotated
from pydantic import BaseModel, Field, model_validator, ConfigDict
from fgraphutils import SourceSpan

class TextBlock(BaseModel):
    """Represents a block of plain text in the C preprocessor AST.
    
    This node captures non-preprocessor-directive content, such as regular
    C code, comments, or whitespace that appears between preprocessor directives.
    """
    kind: Literal["text"] = "text"
    content: str
    
    span: SourceSpan
    model_config = ConfigDict(frozen=True)
    
class DefineNode(BaseModel):
    """Represents a #define preprocessor directive for simple macros.
    
    This node captures object-like macro definitions (e.g., #define PI 3.14159).
    For function-like macros, see FunctionDefineNode.
    """
    kind: Literal["define"] = "define"
    name: str
    value: str
    raw: str
    
    span: SourceSpan
    model_config = ConfigDict(frozen=True)
    
class FunctionDefineNode(BaseModel):
    """Represents a #define preprocessor directive for function-like macros.
    
    This node captures function-like macro definitions that take parameters
    (e.g., #define MAX(a, b) ((a) > (b) ? (a) : (b))).
    """
    kind: Literal["function_def"] = "function_def"
    name: str
    value: str
    parameters: list[str]
    raw: str
    
    span: SourceSpan
    model_config = ConfigDict(frozen=True)

class IncludeNode(BaseModel):
    """Represents an #include preprocessor directive.
    
    This node captures both angle-bracket includes (e.g., #include <stdio.h>)
    and quoted includes (e.g., #include "myheader.h").
    """
    kind: Literal["include"] = "include"
    target: str
    raw: str  # Original raw line
    
    span: SourceSpan
    model_config = ConfigDict(frozen=True)

class UnDefNode(BaseModel):
    """Represents an #undef preprocessor directive.
    
    This node captures macro undefinition directives that remove a previously
    defined macro (e.g., #undef DEBUG).
    """
    kind: Literal["undef"] = "undef"
    name: str
    raw: str  # Original raw line
    
    span: SourceSpan
    model_config = ConfigDict(frozen=True)

class PragmaNode(BaseModel):
    """Represents a #pragma preprocessor directive.
    
    This node captures compiler-specific pragma directives used to provide
    implementation-defined instructions to the compiler (e.g., #pragma once,
    #pragma pack(1)).
    """
    kind: Literal["pragma"] = "pragma"
    payload: str
    raw: str  # Original raw line
    
    model_config = ConfigDict(frozen=True)

class ErrorNode(BaseModel):
    """Represents a #error preprocessor directive.
    
    This node captures error directives that cause the preprocessor to emit
    a diagnostic message and stop processing (e.g., #error "Unsupported platform").
    """
    kind: Literal["error"] = "error"
    message: str
    raw: str  # Original raw line
    
    span: SourceSpan
    model_config = ConfigDict(frozen=True)

# forward-declare Node for type checking in pydantic models
Node = Annotated[Union["TextBlock", "DefineNode", "FunctionDefineNode", "IncludeNode", "UnDefNode", "PragmaNode", "ErrorNode", "ConditionalGroup"], Field(discriminator="kind")]

class ConditionalBranch(BaseModel):
    """Represents a branch of #if / #ifdef / #ifndef / #elif"""
    kind: Literal["if", "ifdef", "ifndef", "elif"]
    condition: str
    body: list[Node]
    raw: str  # '#if CONDITION'
    span: SourceSpan
    model_config = ConfigDict(frozen=True)

class ConditionalGroup(BaseModel):
    """Represents a complete conditional compilation group.
    
    This node captures the full structure of conditional compilation directives,
    including #if/#ifdef/#ifndef, optional #elif branches, optional #else,
    and the terminating #endif. The group represents a single logical unit
    of conditional compilation.
    """
    kind: Literal["conditional_group"] = "conditional_group"
    entry: ConditionalBranch
    elifs: Optional[list[ConditionalBranch]] = None
    else_body: Optional[list[Node]] = None
    else_raw: Optional[str] = None
    endif_raw: str = "#endif"

    span: SourceSpan
    model_config = ConfigDict(frozen=True)
    
class FileRoot(BaseModel):
    """Represents the root node of a parsed C source file.
    
    This node serves as the top-level container for all preprocessor directives
    and text blocks found in a C source file, maintaining the order and structure
    of the original file.
    """
    path: str
    items: list[Node]

    span: SourceSpan
    model_config = ConfigDict(frozen=True)
