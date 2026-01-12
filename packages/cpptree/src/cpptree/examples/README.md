# C Preprocessor Examples for Fortran

This directory contains Fortran source files demonstrating various C preprocessor directives that can be parsed by `cpptree`.

## File Descriptions

- **01_simple_defines.f90** - Basic object-like `#define` macros
- **02_function_like_defines.f90** - Function-like `#define` macros with parameters
- **03_includes.f90** - `#include` directives (both quoted and angle-bracket)
- **04_conditional_if.f90** - `#if`, `#elif`, `#else`, `#endif` conditional compilation
- **05_conditional_ifdef.f90** - `#ifdef` and `#ifndef` macro existence checks
- **06_undef.f90** - `#undef` macro undefinition
- **07_pragma.f90** - `#pragma` compiler-specific directives
- **08_error.f90** - `#error` error generation directives
- **09_complex_conditionals.f90** - Nested and complex conditional structures
- **10_mixed_content.f90** - Mixed preprocessor directives and regular Fortran code
- **11_multiple_elif.f90** - Conditional compilation with multiple `#elif` branches
- **12_nested_groups.f90** - Deeply nested conditional compilation groups
- **13_define_undef_cycle.f90** - Macro definition/undefinition cycles
- **14_text_blocks.f90** - Text blocks between preprocessor directives

## Preprocessor Directives Covered

### Basic Directives
- `#define` - Object-like and function-like macro definitions
- `#undef` - Macro undefinition
- `#include` - File inclusion

### Conditional Compilation
- `#if` - Conditional compilation based on expression
- `#ifdef` - Conditional compilation if macro is defined
- `#ifndef` - Conditional compilation if macro is not defined
- `#elif` - Alternative condition
- `#else` - Default branch
- `#endif` - End conditional block

### Special Directives
- `#pragma` - Compiler-specific pragmas
- `#error` - Error generation

### AST Node Types Demonstrated
- `TextBlock` - Plain text/code blocks
- `DefineNode` - Simple macro definitions
- `FunctionDefineNode` - Function-like macro definitions
- `IncludeNode` - Include directives
- `UnDefNode` - Undef directives
- `PragmaNode` - Pragma directives
- `ErrorNode` - Error directives
- `ConditionalBranch` - Individual conditional branches
- `ConditionalGroup` - Complete conditional compilation groups

These examples can be used to test and validate the `cpptree` parser's ability to correctly identify and structure all types of preprocessor directives in Fortran source files.

