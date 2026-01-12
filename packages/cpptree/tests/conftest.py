"""Pytest configuration and fixtures for cpptree tests."""

import pytest
from tree_sitter import Language, Parser, Tree
import tree_sitter_fortran


@pytest.fixture
def parser():
    """Create a tree-sitter parser for Fortran."""
    language = Language(tree_sitter_fortran.language())
    return Parser(language)


@pytest.fixture
def parse_source(parser):
    """Parse source code and return (source_bytes, tree)."""
    def _parse(source: str) -> tuple[bytes, Tree]:
        source_bytes = source.encode('utf-8')
        tree = parser.parse(source_bytes)
        return source_bytes, tree
    return _parse


@pytest.fixture
def find_node_by_type():
    """Helper to find a node by type in the tree."""
    def _find(node, node_type: str):
        """Recursively find first node of given type."""
        if node.type == node_type:
            return node
        for child in node.children:
            result = _find(child, node_type)
            if result is not None:
                return result
        return None
    return _find

