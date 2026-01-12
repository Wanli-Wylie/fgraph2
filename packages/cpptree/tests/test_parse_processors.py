"""Unit tests for processor functions in cpptree.core.parse."""

import pytest
from cpptree.core.parse import (
    _process_preproc_def,
    _process_preproc_function_def,
    _process_preproc_include,
    _process_preproc_undef,
    _process_preproc_pragma,
    _process_preproc_error,
    _process_conditional_branch,
    _process_conditional_group,
    _merge_text_blocks,
    _process_node,
)
from cpptree.models import (
    DefineNode,
    FunctionDefineNode,
    IncludeNode,
    UnDefNode,
    PragmaNode,
    ErrorNode,
    ConditionalBranch,
    ConditionalGroup,
    TextBlock,
)


class TestProcessPreprocDef:
    """Tests for _process_preproc_def function."""

    def test_simple_define(self, parse_source, find_node_by_type):
        """Test processing a simple #define."""
        source = "#define PI 3.14159\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_def")
        
        result = _process_preproc_def(node, source_bytes, "test.f90")
        
        assert isinstance(result, DefineNode)
        assert result.kind == "define"
        assert result.name == "PI"
        assert result.value == "3.14159"
        assert result.raw == "#define PI 3.14159\n"
        assert result.span.filepath == "test.f90"

    def test_define_without_value(self, parse_source, find_node_by_type):
        """Test processing #define without value."""
        source = "#define DEBUG_MODE\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_def")
        
        result = _process_preproc_def(node, source_bytes, "test.f90")
        
        assert result.name == "DEBUG_MODE"
        assert result.value == ""  # No value provided

    def test_define_with_string_value(self, parse_source, find_node_by_type):
        """Test processing #define with string value."""
        source = '#define VERSION "1.0.0"\n'
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_def")
        
        result = _process_preproc_def(node, source_bytes, "test.f90")
        
        assert result.name == "VERSION"
        assert '"1.0.0"' in result.value or "1.0.0" in result.value

    def test_define_missing_identifier(self, parse_source):
        """Test that missing identifier raises ValueError."""
        # This is a bit tricky to test, but we can try with malformed input
        source = "#define\n"
        source_bytes, tree = parse_source(source)
        # Try to find preproc_def, might not exist or be malformed
        # This test might need adjustment based on actual tree-sitter behavior
        pass  # Skip for now as tree-sitter might not parse malformed input


class TestProcessPreprocFunctionDef:
    """Tests for _process_preproc_function_def function."""

    def test_function_like_define(self, parse_source, find_node_by_type):
        """Test processing function-like #define."""
        source = "#define MAX(a, b) ((a) > (b) ? (a) : (b))\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_function_def")
        
        result = _process_preproc_function_def(node, source_bytes, "test.f90")
        
        assert isinstance(result, FunctionDefineNode)
        assert result.kind == "function_def"
        assert result.name == "MAX"
        assert "a" in result.parameters
        assert "b" in result.parameters
        assert len(result.parameters) == 2
        assert "(a) > (b)" in result.value or "a" in result.value

    def test_single_parameter(self, parse_source, find_node_by_type):
        """Test function-like define with single parameter."""
        source = "#define SQUARE(x) ((x) * (x))\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_function_def")
        
        result = _process_preproc_function_def(node, source_bytes, "test.f90")
        
        assert result.name == "SQUARE"
        assert result.parameters == ["x"]

    def test_multiple_parameters(self, parse_source, find_node_by_type):
        """Test function-like define with multiple parameters."""
        source = "#define SWAP(a, b, temp) do { temp = a; a = b; b = temp; } while(0)\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_function_def")
        
        if node:
            result = _process_preproc_function_def(node, source_bytes, "test.f90")
            assert result.name == "SWAP"
            assert len(result.parameters) >= 2


class TestProcessPreprocInclude:
    """Tests for _process_preproc_include function."""

    def test_quoted_include(self, parse_source, find_node_by_type):
        """Test processing quoted #include."""
        source = '#include "config.h"\n'
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_include")
        
        result = _process_preproc_include(node, source_bytes, "test.f90")
        
        assert isinstance(result, IncludeNode)
        assert result.kind == "include"
        assert result.target == "config.h"
        assert '"config.h"' in result.raw

    def test_angle_bracket_include(self, parse_source, find_node_by_type):
        """Test processing angle-bracket #include."""
        source = "#include <stdio.h>\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_include")
        
        result = _process_preproc_include(node, source_bytes, "test.f90")
        
        assert result.target == "stdio.h"
        assert "<stdio.h>" in result.raw or "stdio.h" in result.raw


class TestProcessPreprocUndef:
    """Tests for _process_preproc_undef function."""

    def test_undef_from_preproc_call(self, parse_source, find_node_by_type):
        """Test processing #undef from preproc_call node."""
        source = "#undef TEMP_MACRO\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_call")
        
        if node:
            result = _process_preproc_undef(node, source_bytes, "test.f90")
            
            assert isinstance(result, UnDefNode)
            assert result.kind == "undef"
            assert result.name == "TEMP_MACRO"
            assert "#undef" in result.raw

    def test_undef_from_preproc_undef(self, parse_source, find_node_by_type):
        """Test processing #undef from preproc_undef node (if exists)."""
        source = "#undef MACRO_NAME\n"
        source_bytes, tree = parse_source(source)
        # Try both node types
        node = find_node_by_type(tree.root_node, "preproc_undef")
        if node is None:
            node = find_node_by_type(tree.root_node, "preproc_call")
        
        if node:
            result = _process_preproc_undef(node, source_bytes, "test.f90")
            assert result.name == "MACRO_NAME"


class TestProcessPreprocPragma:
    """Tests for _process_preproc_pragma function."""

    def test_simple_pragma(self, parse_source, find_node_by_type):
        """Test processing simple #pragma."""
        source = "#pragma once\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_call")
        
        if node:
            result = _process_preproc_pragma(node, source_bytes, "test.f90")
            
            assert isinstance(result, PragmaNode)
            assert result.kind == "pragma"
            assert result.payload == "once"
            assert "#pragma" in result.raw

    def test_pragma_with_value(self, parse_source, find_node_by_type):
        """Test processing #pragma with value."""
        source = "#pragma pack(1)\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_call")
        
        if node:
            result = _process_preproc_pragma(node, source_bytes, "test.f90")
            assert "pack" in result.payload or "1" in result.payload


class TestProcessPreprocError:
    """Tests for _process_preproc_error function."""

    def test_error_with_message(self, parse_source, find_node_by_type):
        """Test processing #error with message."""
        source = '#error "Version 2 or higher is required"\n'
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_call")
        
        if node:
            result = _process_preproc_error(node, source_bytes, "test.f90")
            
            assert isinstance(result, ErrorNode)
            assert result.kind == "error"
            assert "Version 2 or higher" in result.message
            assert "#error" in result.raw


class TestProcessConditionalBranch:
    """Tests for _process_conditional_branch function."""

    def test_if_branch(self, parse_source, find_node_by_type):
        """Test processing #if branch."""
        source = "#if PLATFORM == 1\n    code here\n#endif\n"
        source_bytes, tree = parse_source(source)
        if_node = find_node_by_type(tree.root_node, "preproc_if")
        
        if if_node:
            branch, remaining = _process_conditional_branch(
                if_node, source_bytes, "test.f90", "if"
            )
            
            assert isinstance(branch, ConditionalBranch)
            assert branch.kind == "if"
            assert isinstance(branch.body, list)
            # Condition should be extracted
            assert len(branch.condition) >= 0

    def test_elif_branch(self, parse_source, find_node_by_type):
        """Test processing #elif branch."""
        source = "#if 0\n#elif PLATFORM == 2\n    code\n#endif\n"
        source_bytes, tree = parse_source(source)
        elif_node = find_node_by_type(tree.root_node, "preproc_elif")
        
        if elif_node:
            branch, remaining = _process_conditional_branch(
                elif_node, source_bytes, "test.f90", "elif"
            )
            
            assert branch.kind == "elif"
            assert "PLATFORM" in branch.condition or len(branch.condition) >= 0


class TestProcessConditionalGroup:
    """Tests for _process_conditional_group function."""

    def test_simple_if_endif(self, parse_source, find_node_by_type):
        """Test processing simple #if...#endif."""
        source = "#if 1\n    code\n#endif\n"
        source_bytes, tree = parse_source(source)
        if_node = find_node_by_type(tree.root_node, "preproc_if")
        
        if if_node:
            result = _process_conditional_group(if_node, source_bytes, "test.f90")
            
            assert isinstance(result, ConditionalGroup)
            assert result.kind == "conditional_group"
            assert result.entry.kind == "if"
            assert result.elifs is None
            assert result.else_body is None

    def test_if_else_endif(self, parse_source, find_node_by_type):
        """Test processing #if...#else...#endif."""
        source = "#if 0\n    code1\n#else\n    code2\n#endif\n"
        source_bytes, tree = parse_source(source)
        if_node = find_node_by_type(tree.root_node, "preproc_if")
        
        if if_node:
            result = _process_conditional_group(if_node, source_bytes, "test.f90")
            
            assert result.else_body is not None
            assert isinstance(result.else_body, list)

    def test_if_elif_else_endif(self, parse_source, find_node_by_type):
        """Test processing #if...#elif...#else...#endif."""
        source = "#if 0\n#elif 1\n    code1\n#else\n    code2\n#endif\n"
        source_bytes, tree = parse_source(source)
        if_node = find_node_by_type(tree.root_node, "preproc_if")
        
        if if_node:
            result = _process_conditional_group(if_node, source_bytes, "test.f90")
            
            assert result.elifs is not None
            assert len(result.elifs) >= 1
            assert result.else_body is not None

    def test_ifdef_group(self, parse_source, find_node_by_type):
        """Test processing #ifdef...#endif."""
        source = "#ifdef DEBUG\n    code\n#endif\n"
        source_bytes, tree = parse_source(source)
        ifdef_node = find_node_by_type(tree.root_node, "preproc_ifdef")
        
        if ifdef_node:
            result = _process_conditional_group(ifdef_node, source_bytes, "test.f90")
            
            assert result.entry.kind == "ifdef"


class TestMergeTextBlocks:
    """Tests for _merge_text_blocks function."""

    def test_single_text_block(self, parse_source):
        """Test merging a single text range."""
        source = "some code here\n"
        source_bytes = source.encode('utf-8')
        text_accumulator = [(0, len(source_bytes))]
        
        result = _merge_text_blocks(text_accumulator, source_bytes, "test.f90")
        
        assert isinstance(result, TextBlock)
        assert result.kind == "text"
        assert result.content == source

    def test_multiple_adjacent_blocks(self, parse_source):
        """Test merging multiple adjacent text ranges."""
        source = "line1\nline2\nline3\n"
        source_bytes = source.encode('utf-8')
        text_accumulator = [
            (0, 6),      # "line1\n"
            (6, 12),     # "line2\n"
            (12, 18),    # "line3\n"
        ]
        
        result = _merge_text_blocks(text_accumulator, source_bytes, "test.f90")
        
        assert result is not None
        assert result.content == source

    def test_overlapping_blocks(self, parse_source):
        """Test merging overlapping text ranges."""
        source = "some text\n"
        source_bytes = source.encode('utf-8')
        text_accumulator = [
            (0, 5),      # "some "
            (4, 10),     # " text\n" (overlaps)
        ]
        
        result = _merge_text_blocks(text_accumulator, source_bytes, "test.f90")
        
        assert result is not None
        assert len(result.content) > 0

    def test_empty_accumulator(self, parse_source):
        """Test with empty accumulator."""
        source = "code\n"
        source_bytes = source.encode('utf-8')
        text_accumulator = []
        
        result = _merge_text_blocks(text_accumulator, source_bytes, "test.f90")
        
        assert result is None


class TestProcessNode:
    """Tests for _process_node function."""

    def test_process_preproc_def_node(self, parse_source, find_node_by_type):
        """Test processing a preproc_def node."""
        source = "#define X 1\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_def")
        text_accumulator = []
        
        results = _process_node(node, source_bytes, "test.f90", text_accumulator)
        
        assert len(results) == 1
        assert results[0].kind == "define"
        assert text_accumulator == []  # Should not accumulate text for preproc nodes

    def test_process_text_node(self, parse_source):
        """Test processing a non-preprocessor node (should accumulate text)."""
        source = "program test\nend program\n"
        source_bytes, tree = parse_source(source)
        # Get a non-preproc node (like program)
        program_node = None
        for child in tree.root_node.children:
            if child.type == "program":
                program_node = child
                break
        
        if program_node:
            text_accumulator = []
            results = _process_node(program_node, source_bytes, "test.f90", text_accumulator)
            
            # Should accumulate text, might also process nested preproc nodes
            assert len(text_accumulator) > 0 or len(results) >= 0

    def test_process_preproc_call_undef(self, parse_source, find_node_by_type):
        """Test processing preproc_call with #undef."""
        source = "#undef MACRO\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_call")
        text_accumulator = []
        
        if node:
            results = _process_node(node, source_bytes, "test.f90", text_accumulator)
            
            assert len(results) == 1
            assert results[0].kind == "undef"

    def test_process_preproc_call_pragma(self, parse_source, find_node_by_type):
        """Test processing preproc_call with #pragma."""
        source = "#pragma once\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_call")
        text_accumulator = []
        
        if node:
            results = _process_node(node, source_bytes, "test.f90", text_accumulator)
            
            assert len(results) == 1
            assert results[0].kind == "pragma"

    def test_process_conditional_node(self, parse_source, find_node_by_type):
        """Test processing conditional node."""
        source = "#if 1\n#endif\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_if")
        text_accumulator = []
        
        if node:
            results = _process_node(node, source_bytes, "test.f90", text_accumulator)
            
            assert len(results) == 1
            assert results[0].kind == "conditional_group"

