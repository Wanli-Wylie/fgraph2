"""Parser for C preprocessor directives in Fortran files using tree-sitter.

Refactored to use Visitor pattern + Generator-based stream processing + Python pattern matching.
"""

from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Iterator, Optional, Generator
from tree_sitter import Node, Tree
from fgraphutils import SourceSpan

from cpptree.models import (
    FileRoot,
    PreprocessorNode as ModelNode,
    TextBlock,
    DefineNode,
    FunctionDefineNode,
    IncludeNode,
    UnDefNode,
    PragmaNode,
    ErrorNode,
    ConditionalGroup,
    ConditionalBranch,
)


# ============================================================================
# Helper Functions
# ============================================================================

def _node_to_span(node: Node, filepath: str = "") -> SourceSpan:
    """Convert a tree-sitter node to a SourceSpan.
    
    Args:
        node: The tree-sitter node
        filepath: The file path (defaults to empty string)
    
    Returns:
        A SourceSpan object with the node's location
    """
    start_point = node.start_point
    end_point = node.end_point
    
    return SourceSpan(
        filepath=filepath,
        start_line=start_point.row + 1,  # tree-sitter uses 0-indexed, SourceSpan uses 1-indexed
        start_column=start_point.column,
        end_line=end_point.row + 1,
        end_column=end_point.column,
    )


def _get_node_text(node: Node, source_code: bytes) -> str:
    """Extract text content from a tree-sitter node.
    
    Args:
        node: The tree-sitter node
        source_code: The source code as bytes
    
    Returns:
        The text content as a string
    """
    return source_code[node.start_byte:node.end_byte].decode('utf-8')


def _find_child_by_type(node: Node, node_type: str) -> Optional[Node]:
    """Find the first child node of a specific type.
    
    Args:
        node: The parent node
        node_type: The type of child to find
    
    Returns:
        The first matching child node, or None if not found
    """
    for child in node.children:
        if child.type == node_type:
            return child
    return None


def _find_children_by_type(node: Node, node_type: str) -> Iterator[Node]:
    """Find all child nodes of a specific type.
    
    Args:
        node: The parent node
        node_type: The type of children to find
    
    Yields:
        Matching child nodes
    """
    for child in node.children:
        if child.type == node_type:
            yield child


# ============================================================================
# Visitor Pattern Implementation
# ============================================================================

class NodeVisitor(ABC):
    """Base visitor class for traversing tree-sitter nodes."""
    
    def __init__(self, source_code: bytes, filepath: str):
        """Initialize the visitor.
        
        Args:
            source_code: The source code as bytes
            filepath: The file path
        """
        self.source_code = source_code
        self.filepath = filepath
    
    @abstractmethod
    def visit(self, node: Node) -> Iterator[ModelNode]:
        """Visit a node and yield processed model nodes.
        
        Args:
            node: The tree-sitter node to visit
        
        Yields:
            Processed model nodes
        """
        pass


class PreprocessorVisitor(NodeVisitor):
    """Visitor for processing preprocessor directives using pattern matching."""
    
    def visit(self, node: Node) -> Iterator[ModelNode]:
        """Visit a node and yield processed model nodes.
        
        Uses Python pattern matching to dispatch to appropriate handlers.
        
        Args:
            node: The tree-sitter node to visit
        
        Yields:
            Processed model nodes
        """
        match node.type:
            case "preproc_def":
                yield self._visit_preproc_def(node)
            case "preproc_function_def":
                yield self._visit_preproc_function_def(node)
            case "preproc_include":
                yield self._visit_preproc_include(node)
            case "preproc_undef":
                yield self._visit_preproc_undef(node)
            case "preproc_pragma":
                yield self._visit_preproc_pragma(node)
            case "preproc_error":
                yield self._visit_preproc_error(node)
            case "preproc_call":
                yield from self._visit_preproc_call(node)
            case "preproc_if" | "preproc_ifdef" | "preproc_ifndef":
                yield self._visit_conditional_group(node)
            case _:
                # Non-preprocessor node - process children recursively
                yield from self._visit_non_preproc_node(node)
    
    def _visit_preproc_def(self, node: Node) -> DefineNode:
        """Visit a preproc_def node (simple #define)."""
        raw = _get_node_text(node, self.source_code)
        
        identifier_node = _find_child_by_type(node, "identifier")
        if identifier_node is None:
            raise ValueError(f"preproc_def missing identifier: {raw}")
        
        name = _get_node_text(identifier_node, self.source_code)
        arg_node = _find_child_by_type(node, "preproc_arg")
        value = _get_node_text(arg_node, self.source_code) if arg_node else ""
        
        return DefineNode(
            name=name,
            value=value,
            raw=raw,
            span=_node_to_span(node, self.filepath),
        )
    
    def _visit_preproc_function_def(self, node: Node) -> FunctionDefineNode:
        """Visit a preproc_function_def node (function-like #define)."""
        raw = _get_node_text(node, self.source_code)
        
        identifier_node = _find_child_by_type(node, "identifier")
        if identifier_node is None:
            raise ValueError(f"preproc_function_def missing identifier: {raw}")
        
        name = _get_node_text(identifier_node, self.source_code)
        
        # Extract parameters from preproc_params
        parameters = []
        params_node = _find_child_by_type(node, "preproc_params")
        if params_node is not None:
            for param_child in params_node.children:
                if param_child.type == "identifier":
                    parameters.append(_get_node_text(param_child, self.source_code))
        
        arg_node = _find_child_by_type(node, "preproc_arg")
        value = _get_node_text(arg_node, self.source_code) if arg_node else ""
        
        return FunctionDefineNode(
            name=name,
            value=value,
            parameters=parameters,
            raw=raw,
            span=_node_to_span(node, self.filepath),
        )
    
    def _visit_preproc_include(self, node: Node) -> IncludeNode:
        """Visit a preproc_include node (#include)."""
        raw = _get_node_text(node, self.source_code)
        
        # Find the include target (string_literal, system_lib_string, or system_lib_name)
        target_node = (
            _find_child_by_type(node, "string_literal") or
            _find_child_by_type(node, "system_lib_string") or
            _find_child_by_type(node, "system_lib_name")
        )
        
        if target_node is None:
            raise ValueError(f"preproc_include missing target: {raw}")
        
        target = _get_node_text(target_node, self.source_code).strip('"<>')
        
        return IncludeNode(
            target=target,
            raw=raw,
            span=_node_to_span(node, self.filepath),
        )
    
    def _visit_preproc_undef(self, node: Node) -> UnDefNode:
        """Visit a preproc_undef node (#undef)."""
        raw = _get_node_text(node, self.source_code)
        
        arg_node = _find_child_by_type(node, "preproc_arg")
        if arg_node:
            name = _get_node_text(arg_node, self.source_code).strip()
        else:
            identifier_node = _find_child_by_type(node, "identifier")
            if identifier_node is None:
                raise ValueError(f"preproc_undef missing identifier: {raw}")

            name = _get_node_text(identifier_node, self.source_code)
        
        return UnDefNode(
            name=name,
            raw=raw,
            span=_node_to_span(node, self.filepath),
        )
    
    def _visit_preproc_pragma(self, node: Node) -> PragmaNode:
        """Visit a preproc_pragma node (#pragma)."""
        raw = _get_node_text(node, self.source_code)
        
        arg_node = _find_child_by_type(node, "preproc_arg")
        if arg_node:
            payload = _get_node_text(arg_node, self.source_code).strip()
        else:
            payload = raw.replace("#pragma", "").strip()
        
        return PragmaNode(
            payload=payload,
            raw=raw,
            span=_node_to_span(node, self.filepath),
        )
    
    def _visit_preproc_error(self, node: Node) -> ErrorNode:
        """Visit a preproc_error node (#error)."""
        raw = _get_node_text(node, self.source_code)
        
        arg_node = _find_child_by_type(node, "preproc_arg")
        if arg_node:
            message = _get_node_text(arg_node, self.source_code).strip('"')
        else:
            message = raw.replace("#error", "").strip().strip('"')
        
        return ErrorNode(
            message=message,
            raw=raw,
            span=_node_to_span(node, self.filepath),
        )
    
    def _visit_preproc_call(self, node: Node) -> Iterator[ModelNode]:
        """Visit a preproc_call node (can be #undef, #pragma, or #error)."""
        directive_node = _find_child_by_type(node, "preproc_directive")
        if directive_node:
            directive_text = _get_node_text(directive_node, self.source_code).strip()
            match directive_text:
                case "#undef":
                    yield self._visit_preproc_undef(node)
                case "#pragma":
                    yield self._visit_preproc_pragma(node)
                case "#error":
                    yield self._visit_preproc_error(node)
                case _:
                    # Unknown directive - treat as text
                    yield from self._yield_text_block(node)
        else:
            # No directive found - treat as text
            yield from self._yield_text_block(node)
    
    def _visit_conditional_group(self, node: Node) -> ConditionalGroup:
        """Visit a conditional group (preproc_if, preproc_ifdef, preproc_ifndef)."""
        # Determine the kind of conditional
        match node.type:
            case "preproc_ifdef":
                entry_kind = "ifdef"
            case "preproc_ifndef":
                entry_kind = "ifndef"
            case _:  # preproc_if
                entry_kind = "if"
        
        # Process the entry branch
        entry, remaining = self._process_conditional_branch(node, entry_kind)
        
        # Process elif branches
        elifs: list[ConditionalBranch] = []
        else_body: Optional[list[ModelNode]] = None
        else_raw: Optional[str] = None
        endif_raw = "#endif"
        
        # Process remaining children (elif, else, endif)
        i = 0
        while i < len(remaining):
            child = remaining[i]
            
            match child.type:
                case "preproc_elif":
                    elif_branch, next_remaining = self._process_conditional_branch(child, "elif")
                    elifs.append(elif_branch)
                    remaining = next_remaining
                    i = 0  # Reset to process new remaining children
                case "preproc_else":
                    else_raw = _get_node_text(child, self.source_code).split("\n")[0]
                    # Collect body until endif
                    else_body = []
                    i += 1
                    while i < len(remaining) and remaining[i].type != "#endif":
                        processed = list(self.visit(remaining[i]))
                        if processed:
                            else_body.extend(processed)
                        i += 1
                    break
                case "#endif":
                    endif_raw = _get_node_text(child, self.source_code).split("\n")[0]
                    break
                case _:
                    i += 1
        
        # Calculate span for the entire group
        start_span = entry.span
        end_node = node
        for descendant in node.children:
            if descendant.type == "#endif":
                end_node = descendant
                break
        
        group_span = SourceSpan(
            filepath=self.filepath,
            start_line=start_span.start_line,
            start_column=start_span.start_column,
            end_line=_node_to_span(end_node, self.filepath).end_line,
            end_column=_node_to_span(end_node, self.filepath).end_column,
        )
        
        return ConditionalGroup(
            entry=entry,
            elifs=elifs if elifs else None,
            else_body=else_body,
            else_raw=else_raw,
            endif_raw=endif_raw,
            span=group_span,
        )
    
    def _process_conditional_branch(
        self,
        node: Node,
        branch_kind: str,
    ) -> tuple[ConditionalBranch, list[Node]]:
        """Process a conditional branch (if/ifdef/ifndef/elif).
        
        Args:
            node: The conditional node
            branch_kind: The kind of branch ("if", "ifdef", "ifndef", or "elif")
        
        Returns:
            A tuple of (ConditionalBranch, remaining_children)
        """
        raw = _get_node_text(node, self.source_code).split("\n")[0]
        condition = self._extract_condition_text(node)
        
        # Collect body nodes (everything between the directive and next preproc directive)
        body: list[ModelNode] = []
        remaining_children: list[Node] = []
        
        in_body = False
        for child in node.children:
            match child.type:
                case "#if" | "#ifdef" | "#ifndef" | "#elif":
                    in_body = True
                case "preproc_elif" | "preproc_else" | "#endif":
                    remaining_children.append(child)
                    break
                case _:
                    if in_body:
                        # Process this child as part of the body
                        processed = list(self.visit(child))
                        if processed:
                            body.extend(processed)
        
        return ConditionalBranch(
            kind=branch_kind,
            condition=condition,
            body=body,
            raw=raw,
            span=_node_to_span(node, self.filepath),
        ), remaining_children
    
    def _extract_condition_text(self, node: Node) -> str:
        """Extract the condition text from a conditional node."""
        directive_types = ["#if", "#ifdef", "#ifndef", "#elif"]
        found_directive = False
        
        for child in node.children:
            if child.type in directive_types:
                found_directive = True
                continue
            
            if found_directive:
                # After the directive, we should find the condition
                if child.type in ["\n", "comment", "subroutine", "program", "preproc_elif", "preproc_else", "#endif"]:
                    break
                return _get_node_text(child, self.source_code).strip()
        
        # Fallback: try to extract from raw text
        raw = _get_node_text(node, self.source_code)
        for directive in ["#if", "#ifdef", "#ifndef", "#elif"]:
            if raw.startswith(directive):
                rest = raw[len(directive):].strip()
                if "\n" in rest:
                    return rest.split("\n")[0].strip()
                return rest.strip()
        
        return ""
    
    def _visit_non_preproc_node(self, node: Node) -> Iterator[ModelNode]:
        """Visit a non-preprocessor node, recursively processing children."""
        # Check if this node has preprocessor children
        has_preproc_children = any(
            child.type.startswith("preproc_") or
            (child.type == "preproc_call" and _find_child_by_type(child, "preproc_directive"))
            for child in node.children
        )
        
        if has_preproc_children:
            # Process children, yielding text blocks and preprocessor nodes
            text_ranges: list[tuple[int, int]] = []
            
            for child in node.children:
                if child.type.startswith("preproc_") or (
                    child.type == "preproc_call" and _find_child_by_type(child, "preproc_directive")
                ):
                    # Yield accumulated text before this preproc node
                    if text_ranges:
                        text_block = self._merge_text_ranges(text_ranges)
                        if text_block:
                            yield text_block
                        text_ranges = []
                    
                    # Process the preproc node
                    yield from self.visit(child)
                else:
                    # Non-preproc child - recursively process it
                    child_results = list(self.visit(child))
                    if child_results:
                        # If we have accumulated text, create a text block first
                        if text_ranges:
                            text_block = self._merge_text_ranges(text_ranges)
                            if text_block:
                                yield text_block
                            text_ranges = []
                        yield from child_results
                    else:
                        # No preproc nodes found - accumulate text
                        text_ranges.append((child.start_byte, child.end_byte))
            
            # Yield any remaining text
            if text_ranges:
                text_block = self._merge_text_ranges(text_ranges)
                if text_block:
                    yield text_block
        else:
            # No preproc children - yield as text block
            yield from self._yield_text_block(node)
    
    def _yield_text_block(self, node: Node) -> Iterator[TextBlock]:
        """Yield a text block for a node."""
        text_block = self._create_text_block(node.start_byte, node.end_byte)
        if text_block:
            yield text_block
    
    def _create_text_block(self, start_byte: int, end_byte: int) -> Optional[TextBlock]:
        """Create a text block from byte range."""
        if start_byte >= end_byte:
            return None
        
        content = self.source_code[start_byte:end_byte].decode('utf-8')
        
        # Calculate span
        text_before = self.source_code[:start_byte].decode('utf-8')
        lines_before = text_before.count('\n')
        col_before = len(text_before.split('\n')[-1])
        
        text_after = self.source_code[:end_byte].decode('utf-8')
        lines_after = text_after.count('\n')
        col_after = len(text_after.split('\n')[-1])
        
        return TextBlock(
            content=content,
            span=SourceSpan(
                filepath=self.filepath,
                start_line=lines_before + 1,
                start_column=col_before,
                end_line=lines_after + 1,
                end_column=col_after,
            ),
        )
    
    def _merge_text_ranges(
        self,
        text_ranges: list[tuple[int, int]],
    ) -> Optional[TextBlock]:
        """Merge accumulated text ranges into a single TextBlock."""
        if not text_ranges:
            return None
        
        # Sort and merge overlapping or adjacent ranges
        sorted_ranges = sorted(text_ranges)
        merged: list[tuple[int, int]] = []
        
        for start, end in sorted_ranges:
            if not merged:
                merged.append((start, end))
            else:
                last_start, last_end = merged[-1]
                if start <= last_end:
                    merged[-1] = (last_start, max(last_end, end))
                else:
                    merged.append((start, end))
        
        # Extract text from merged ranges
        text_parts = []
        for start, end in merged:
            text_parts.append(self.source_code[start:end].decode('utf-8'))
        
        content = "".join(text_parts)
        
        # Create span from first and last ranges
        first_start, _ = merged[0]
        _, last_end = merged[-1]
        
        return self._create_text_block(first_start, last_end)


# ============================================================================
# Stream Processing Functions
# ============================================================================

def process_tree_stream(
    root: Node,
    source_code: bytes,
    filepath: str = "",
) -> Generator[ModelNode, None, None]:
    """Process a tree-sitter tree using generator-based stream processing.
    
    Args:
        root: The root node of the tree-sitter tree
        source_code: The source code as bytes
        filepath: The file path
    
    Yields:
        Processed model nodes in order
    """
    visitor = PreprocessorVisitor(source_code, filepath)
    text_ranges: list[tuple[int, int]] = []
    
    for child in root.children:
        # Process child and collect results
        results = list(visitor.visit(child))
        
        if results:
            # Yield accumulated text before preprocessor nodes
            if text_ranges:
                text_block = visitor._merge_text_ranges(text_ranges)
                if text_block:
                    yield text_block
                text_ranges = []
            
            # Yield preprocessor nodes
            yield from results
        else:
            # No preprocessor nodes found - accumulate text
            text_ranges.append((child.start_byte, child.end_byte))
    
    # Yield any remaining text
    if text_ranges:
        text_block = visitor._merge_text_ranges(text_ranges)
        if text_block:
            yield text_block


# ============================================================================
# Legacy API Compatibility (for backward compatibility with tests)
# ============================================================================

def _process_preproc_def(node: Node, source_code: bytes, filepath: str) -> DefineNode:
    """Legacy function for backward compatibility."""
    visitor = PreprocessorVisitor(source_code, filepath)
    return visitor._visit_preproc_def(node)


def _process_preproc_function_def(node: Node, source_code: bytes, filepath: str) -> FunctionDefineNode:
    """Legacy function for backward compatibility."""
    visitor = PreprocessorVisitor(source_code, filepath)
    return visitor._visit_preproc_function_def(node)


def _process_preproc_include(node: Node, source_code: bytes, filepath: str) -> IncludeNode:
    """Legacy function for backward compatibility."""
    visitor = PreprocessorVisitor(source_code, filepath)
    return visitor._visit_preproc_include(node)


def _process_preproc_undef(node: Node, source_code: bytes, filepath: str) -> UnDefNode:
    """Legacy function for backward compatibility."""
    visitor = PreprocessorVisitor(source_code, filepath)
    return visitor._visit_preproc_undef(node)


def _process_preproc_pragma(node: Node, source_code: bytes, filepath: str) -> PragmaNode:
    """Legacy function for backward compatibility."""
    visitor = PreprocessorVisitor(source_code, filepath)
    return visitor._visit_preproc_pragma(node)


def _process_preproc_error(node: Node, source_code: bytes, filepath: str) -> ErrorNode:
    """Legacy function for backward compatibility."""
    visitor = PreprocessorVisitor(source_code, filepath)
    return visitor._visit_preproc_error(node)


def _process_conditional_branch(
    node: Node,
    source_code: bytes,
    filepath: str,
    branch_kind: str,
) -> tuple[ConditionalBranch, list[Node]]:
    """Legacy function for backward compatibility."""
    visitor = PreprocessorVisitor(source_code, filepath)
    return visitor._process_conditional_branch(node, branch_kind)


def _process_conditional_group(
    node: Node,
    source_code: bytes,
    filepath: str,
) -> ConditionalGroup:
    """Legacy function for backward compatibility."""
    visitor = PreprocessorVisitor(source_code, filepath)
    return visitor._visit_conditional_group(node)


def _process_node(
    node: Node,
    source_code: bytes,
    filepath: str,
    text_accumulator: list[tuple[int, int]],
) -> list[ModelNode]:
    """Legacy function for backward compatibility."""
    visitor = PreprocessorVisitor(source_code, filepath)
    results = list(visitor.visit(node))
    
    # Handle text accumulator for backward compatibility
    if not results and not any(
        child.type.startswith("preproc_") or
        (child.type == "preproc_call" and _find_child_by_type(child, "preproc_directive"))
        for child in node.children
    ):
        text_accumulator.append((node.start_byte, node.end_byte))
    
    return results


def _merge_text_blocks(
    text_accumulator: list[tuple[int, int]],
    source_code: bytes,
    filepath: str,
) -> Optional[TextBlock]:
    """Legacy function for backward compatibility."""
    visitor = PreprocessorVisitor(source_code, filepath)
    return visitor._merge_text_ranges(text_accumulator)


def _extract_condition_text(node: Node, source_code: bytes) -> str:
    """Legacy function for backward compatibility."""
    visitor = PreprocessorVisitor(source_code, "")
    return visitor._extract_condition_text(node)
