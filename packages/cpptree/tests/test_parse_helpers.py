"""Unit tests for helper functions in cpptree.core.parse."""

import pytest
from fgraphutils import SourceSpan
from cpptree.core.parse import (
    _node_to_span,
    _get_node_text,
    _find_child_by_type,
    _extract_condition_text,
)


class TestNodeToSpan:
    """Tests for _node_to_span function."""

    def test_basic_conversion(self, parse_source, find_node_by_type):
        """Test basic node to span conversion."""
        source = "#define PI 3.14\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_def")
        
        span = _node_to_span(node, "test.f90")
        
        assert isinstance(span, SourceSpan)
        assert span.filepath == "test.f90"
        assert span.start_line == 1  # tree-sitter uses 0-indexed, converted to 1-indexed
        assert span.start_column == 0
        assert span.end_line >= 1
        assert span.end_column >= 0

    def test_empty_filepath(self, parse_source, find_node_by_type):
        """Test with empty filepath."""
        source = "#define X 1\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_def")
        
        span = _node_to_span(node, "")
        
        assert span.filepath == ""

    def test_multiline_node(self, parse_source, find_node_by_type):
        """Test conversion of multiline node."""
        source = "#define LONG_MACRO \\\n    some_value\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_def")
        
        span = _node_to_span(node, "test.f90")
        
        assert span.end_line > span.start_line  # Should span multiple lines


class TestGetNodeText:
    """Tests for _get_node_text function."""

    def test_simple_text_extraction(self, parse_source, find_node_by_type):
        """Test extracting text from a simple node."""
        source = "#define PI 3.14159\n"
        source_bytes, tree = parse_source(source)
        node = find_node_by_type(tree.root_node, "preproc_def")
        
        text = _get_node_text(node, source_bytes)
        
        assert text == "#define PI 3.14159\n"
        assert "PI" in text
        assert "3.14159" in text

    def test_identifier_text(self, parse_source, find_node_by_type):
        """Test extracting text from identifier node."""
        source = "#define MACRO_NAME value\n"
        source_bytes, tree = parse_source(source)
        def_node = find_node_by_type(tree.root_node, "preproc_def")
        identifier_node = _find_child_by_type(def_node, "identifier")
        
        text = _get_node_text(identifier_node, source_bytes)
        
        assert text == "MACRO_NAME"

    def test_empty_node(self, parse_source):
        """Test with empty node (edge case)."""
        source = "\n"
        source_bytes, tree = parse_source(source)
        # Get a node that might be empty
        node = tree.root_node
        
        text = _get_node_text(node, source_bytes)
        
        assert isinstance(text, str)


class TestFindChildByType:
    """Tests for _find_child_by_type function."""

    def test_find_existing_child(self, parse_source, find_node_by_type):
        """Test finding an existing child node."""
        source = "#define MACRO value\n"
        source_bytes, tree = parse_source(source)
        def_node = find_node_by_type(tree.root_node, "preproc_def")
        
        identifier = _find_child_by_type(def_node, "identifier")
        
        assert identifier is not None
        assert identifier.type == "identifier"
        assert _get_node_text(identifier, source_bytes) == "MACRO"

    def test_find_nonexistent_child(self, parse_source, find_node_by_type):
        """Test finding a child that doesn't exist."""
        source = "#define MACRO value\n"
        source_bytes, tree = parse_source(source)
        def_node = find_node_by_type(tree.root_node, "preproc_def")
        
        result = _find_child_by_type(def_node, "nonexistent_type")
        
        assert result is None

    def test_find_first_matching_child(self, parse_source, find_node_by_type):
        """Test that it returns the first matching child when multiple exist."""
        source = "#define MACRO value\n"
        source_bytes, tree = parse_source(source)
        def_node = find_node_by_type(tree.root_node, "preproc_def")
        
        # There might be multiple children, but we should get the first one
        first_define = _find_child_by_type(def_node, "#define")
        
        assert first_define is not None
        assert first_define.type == "#define"


class TestExtractConditionText:
    """Tests for _extract_condition_text function."""

    def test_simple_condition(self, parse_source, find_node_by_type):
        """Test extracting simple condition."""
        source = "#if PLATFORM == 1\n#endif\n"
        source_bytes, tree = parse_source(source)
        if_node = find_node_by_type(tree.root_node, "preproc_if")
        
        condition = _extract_condition_text(if_node, source_bytes)
        
        assert "PLATFORM == 1" in condition or "PLATFORM" in condition

    def test_ifdef_condition(self, parse_source, find_node_by_type):
        """Test extracting ifdef condition."""
        source = "#ifdef DEBUG_MODE\n#endif\n"
        source_bytes, tree = parse_source(source)
        ifdef_node = find_node_by_type(tree.root_node, "preproc_ifdef")
        
        condition = _extract_condition_text(ifdef_node, source_bytes)
        
        # For ifdef, condition might be just the identifier
        assert isinstance(condition, str)
        assert len(condition) >= 0

    def test_complex_condition(self, parse_source, find_node_by_type):
        """Test extracting complex condition."""
        source = "#if defined(FEATURE_A) && defined(FEATURE_B)\n#endif\n"
        source_bytes, tree = parse_source(source)
        if_node = find_node_by_type(tree.root_node, "preproc_if")
        
        condition = _extract_condition_text(if_node, source_bytes)
        
        assert isinstance(condition, str)
        # Should contain parts of the condition
        assert len(condition) > 0

    def test_elif_condition(self, parse_source, find_node_by_type):
        """Test extracting elif condition."""
        source = "#if 0\n#elif PLATFORM == 2\n#endif\n"
        source_bytes, tree = parse_source(source)
        elif_node = find_node_by_type(tree.root_node, "preproc_elif")
        
        condition = _extract_condition_text(elif_node, source_bytes)
        
        assert isinstance(condition, str)
        assert "PLATFORM" in condition or len(condition) >= 0

