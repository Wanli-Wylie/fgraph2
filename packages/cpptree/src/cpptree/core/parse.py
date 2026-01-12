"""Parser for C preprocessor directives in Fortran files using tree-sitter."""

from __future__ import annotations
from typing import Optional
from tree_sitter import Node, Tree
from fgraphutils import SourceSpan

from cpptree.models import (
    FileRoot,
    Node as ModelNode,
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

def _process_preproc_def(node: Node, source_code: bytes, filepath: str) -> DefineNode:
    """Process a preproc_def node (simple #define).
    
    Args:
        node: The preproc_def node
        source_code: The source code as bytes
        filepath: The file path
    
    Returns:
        A DefineNode
    """
    raw = _get_node_text(node, source_code)
    
    # Find identifier (macro name)
    identifier_node = _find_child_by_type(node, "identifier")
    if identifier_node is None:
        raise ValueError(f"preproc_def missing identifier: {raw}")
    
    name = _get_node_text(identifier_node, source_code)
    
    # Find the value (preproc_arg)
    arg_node = _find_child_by_type(node, "preproc_arg")
    value = _get_node_text(arg_node, source_code) if arg_node else ""
    
    return DefineNode(
        name=name,
        value=value,
        raw=raw,
        span=_node_to_span(node, filepath),
    )


def _process_preproc_function_def(node: Node, source_code: bytes, filepath: str) -> FunctionDefineNode:
    """Process a preproc_function_def node (function-like #define).
    
    Args:
        node: The preproc_function_def node
        source_code: The source code as bytes
        filepath: The file path
    
    Returns:
        A FunctionDefineNode
    """
    raw = _get_node_text(node, source_code)
    
    # Find identifier (macro name)
    identifier_node = _find_child_by_type(node, "identifier")
    if identifier_node is None:
        raise ValueError(f"preproc_function_def missing identifier: {raw}")
    
    name = _get_node_text(identifier_node, source_code)
    
    # Extract parameters from preproc_params
    parameters = []
    params_node = _find_child_by_type(node, "preproc_params")
    if params_node is not None:
        for param_child in params_node.children:
            if param_child.type == "identifier":
                parameters.append(_get_node_text(param_child, source_code))
    
    # Find the value (preproc_arg)
    arg_node = _find_child_by_type(node, "preproc_arg")
    value = _get_node_text(arg_node, source_code) if arg_node else ""
    
    return FunctionDefineNode(
        name=name,
        value=value,
        parameters=parameters,
        raw=raw,
        span=_node_to_span(node, filepath),
    )


def _process_preproc_include(node: Node, source_code: bytes, filepath: str) -> IncludeNode:
    """Process a preproc_include node (#include).
    
    Args:
        node: The preproc_include node
        source_code: The source code as bytes
        filepath: The file path
    
    Returns:
        An IncludeNode
    """
    raw = _get_node_text(node, source_code)
    
    # Find the include target (string_literal, system_lib_string, or system_lib_name)
    target_node = _find_child_by_type(node, "string_literal")
    if target_node is None:
        target_node = _find_child_by_type(node, "system_lib_string")
    if target_node is None:
        target_node = _find_child_by_type(node, "system_lib_name")
    
    if target_node is None:
        raise ValueError(f"preproc_include missing target: {raw}")
    
    target = _get_node_text(target_node, source_code).strip('"<>')
    
    return IncludeNode(
        target=target,
        raw=raw,
        span=_node_to_span(node, filepath),
    )


def _process_preproc_undef(node: Node, source_code: bytes, filepath: str) -> UnDefNode:
    """Process a preproc_undef or preproc_call with #undef (#undef).
    
    Args:
        node: The preproc_undef or preproc_call node
        source_code: The source code as bytes
        filepath: The file path
    
    Returns:
        An UnDefNode
    """
    raw = _get_node_text(node, source_code)
    
    # For preproc_call, the name is in preproc_arg
    # For preproc_undef, it might be in identifier
    arg_node = _find_child_by_type(node, "preproc_arg")
    if arg_node:
        name = _get_node_text(arg_node, source_code).strip()
    else:
        identifier_node = _find_child_by_type(node, "identifier")
        if identifier_node is None:
            raise ValueError(f"preproc_undef missing identifier: {raw}")
        name = _get_node_text(identifier_node, source_code)
    
    return UnDefNode(
        name=name,
        raw=raw,
        span=_node_to_span(node, filepath),
    )


def _process_preproc_pragma(node: Node, source_code: bytes, filepath: str) -> PragmaNode:
    """Process a preproc_pragma or preproc_call with #pragma (#pragma).
    
    Args:
        node: The preproc_pragma or preproc_call node
        source_code: The source code as bytes
        filepath: The file path
    
    Returns:
        A PragmaNode
    """
    raw = _get_node_text(node, source_code)
    
    # Extract payload (everything after #pragma)
    # The payload is typically in preproc_arg
    arg_node = _find_child_by_type(node, "preproc_arg")
    if arg_node:
        payload = _get_node_text(arg_node, source_code).strip()
    else:
        # Try to extract from raw text
        payload = raw.replace("#pragma", "").strip()
    
    return PragmaNode(
        payload=payload,
        raw=raw,
        span=_node_to_span(node, filepath),
    )


def _process_preproc_error(node: Node, source_code: bytes, filepath: str) -> ErrorNode:
    """Process a preproc_error or preproc_call with #error (#error).
    
    Args:
        node: The preproc_error or preproc_call node
        source_code: The source code as bytes
        filepath: The file path
    
    Returns:
        An ErrorNode
    """
    raw = _get_node_text(node, source_code)
    
    # Extract error message
    arg_node = _find_child_by_type(node, "preproc_arg")
    if arg_node:
        message = _get_node_text(arg_node, source_code).strip('"')
    else:
        # Try to extract from raw text
        message = raw.replace("#error", "").strip().strip('"')
    
    return ErrorNode(
        message=message,
        raw=raw,
        span=_node_to_span(node, filepath),
    )


def _extract_condition_text(node: Node, source_code: bytes) -> str:
    """Extract the condition text from a conditional node.
    
    Args:
        node: The conditional node (preproc_if, preproc_elif, etc.)
        source_code: The source code as bytes
    
    Returns:
        The condition text
    """
    # Look for the condition expression (could be binary_expression, identifier, etc.)
    # Skip the directive token itself (#if, #ifdef, #ifndef, #elif)
    # The condition is typically the first non-directive, non-whitespace child
    # before any newline or body content
    
    # Find the directive token
    directive_types = ["#if", "#ifdef", "#ifndef", "#elif"]
    found_directive = False
    
    for child in node.children:
        if child.type in directive_types:
            found_directive = True
            continue
        
        if found_directive:
            # After the directive, we should find the condition
            # Stop at newline, body content, or other preproc directives
            if child.type in ["\n", "comment", "subroutine", "program", "preproc_elif", "preproc_else", "#endif"]:
                break
            # This is part of the condition
            return _get_node_text(child, source_code).strip()
    
    # Fallback: try to extract from raw text
    raw = _get_node_text(node, source_code)
    for directive in ["#if", "#ifdef", "#ifndef", "#elif"]:
        if raw.startswith(directive):
            # Extract everything after the directive until newline
            rest = raw[len(directive):].strip()
            if "\n" in rest:
                return rest.split("\n")[0].strip()
            return rest.strip()
    
    return ""


def _process_conditional_branch(
    node: Node,
    source_code: bytes,
    filepath: str,
    branch_kind: str,
) -> tuple[ConditionalBranch, list[Node]]:
    """Process a conditional branch (if/ifdef/ifndef/elif).
    
    Args:
        node: The conditional node
        source_code: The source code as bytes
        filepath: The file path
        branch_kind: The kind of branch ("if", "ifdef", "ifndef", or "elif")
    
    Returns:
        A tuple of (ConditionalBranch, remaining_children)
        where remaining_children are nodes that should be processed separately
    """
    raw = _get_node_text(node, source_code).split("\n")[0]  # Just the directive line
    condition = _extract_condition_text(node, source_code)
    
    # Collect body nodes (everything between the directive and next preproc directive)
    body: list[ModelNode] = []
    remaining_children: list[Node] = []
    
    in_body = False
    for child in node.children:
        if child.type in ["#if", "#ifdef", "#ifndef", "#elif"]:
            in_body = True
            continue
        elif child.type in ["preproc_elif", "preproc_else", "#endif"]:
            remaining_children.append(child)
            break
        elif in_body:
            # Process this child as part of the body
            processed = _process_node(child, source_code, filepath, [])
            if processed:
                body.extend(processed)
    
    return ConditionalBranch(
        kind=branch_kind,
        condition=condition,
        body=body,
        raw=raw,
        span=_node_to_span(node, filepath),
    ), remaining_children


def _process_conditional_group(
    node: Node,
    source_code: bytes,
    filepath: str,
) -> ConditionalGroup:
    """Process a conditional group (preproc_if, preproc_ifdef, preproc_ifndef).
    
    Args:
        node: The conditional node
        source_code: The source code as bytes
        filepath: The file path
    
    Returns:
        A ConditionalGroup
    """
    # Determine the kind of conditional
    if node.type == "preproc_ifdef":
        entry_kind = "ifdef"
    elif node.type == "preproc_ifndef":
        entry_kind = "ifndef"
    else:  # preproc_if
        entry_kind = "if"
    
    # Process the entry branch
    entry, remaining = _process_conditional_branch(node, source_code, filepath, entry_kind)
    
    # Process elif branches
    elifs: list[ConditionalBranch] = []
    else_body: Optional[list[ModelNode]] = None
    else_raw: Optional[str] = None
    endif_raw = "#endif"
    
    # Process remaining children (elif, else, endif)
    i = 0
    while i < len(remaining):
        child = remaining[i]
        
        if child.type == "preproc_elif":
            elif_branch, next_remaining = _process_conditional_branch(child, source_code, filepath, "elif")
            elifs.append(elif_branch)
            remaining = next_remaining
            i = 0  # Reset to process new remaining children
        elif child.type == "preproc_else":
            else_raw = _get_node_text(child, source_code).split("\n")[0]
            # Collect body until endif
            else_body = []
            i += 1
            while i < len(remaining) and remaining[i].type != "#endif":
                processed = _process_node(remaining[i], source_code, filepath, [])
                if processed:
                    else_body.extend(processed)
                i += 1
            break
        elif child.type == "#endif":
            endif_raw = _get_node_text(child, source_code).split("\n")[0]
            break
        else:
            i += 1
    
    # Calculate span for the entire group
    start_span = entry.span
    # Find the end by looking for #endif in the node's descendants
    end_node = node
    for descendant in node.children:
        if descendant.type == "#endif":
            end_node = descendant
            break
    
    group_span = SourceSpan(
        filepath=filepath,
        start_line=start_span.start_line,
        start_column=start_span.start_column,
        end_line=_node_to_span(end_node, filepath).end_line,
        end_column=_node_to_span(end_node, filepath).end_column,
    )
    
    return ConditionalGroup(
        entry=entry,
        elifs=elifs if elifs else None,
        else_body=else_body,
        else_raw=else_raw,
        endif_raw=endif_raw,
        span=group_span,
    )


def _process_node(
    node: Node,
    source_code: bytes,
    filepath: str,
    text_accumulator: list[tuple[int, int]],
) -> list[ModelNode]:
    """Process a single tree-sitter node.
    
    Args:
        node: The tree-sitter node to process
        source_code: The source code as bytes
        filepath: The file path
        text_accumulator: List of (start_byte, end_byte) tuples for text blocks
    
    Returns:
        A list of processed model nodes
    """
    results: list[ModelNode] = []
    
    if node.type.startswith("preproc_"):
        # Process preprocessor directives
        if node.type == "preproc_def":
            results.append(_process_preproc_def(node, source_code, filepath))
        elif node.type == "preproc_function_def":
            results.append(_process_preproc_function_def(node, source_code, filepath))
        elif node.type == "preproc_include":
            results.append(_process_preproc_include(node, source_code, filepath))
        elif node.type == "preproc_undef":
            results.append(_process_preproc_undef(node, source_code, filepath))
        elif node.type == "preproc_pragma":
            results.append(_process_preproc_pragma(node, source_code, filepath))
        elif node.type == "preproc_error":
            results.append(_process_preproc_error(node, source_code, filepath))
        elif node.type == "preproc_call":
            # preproc_call can be #undef, #pragma, or #error
            directive_node = _find_child_by_type(node, "preproc_directive")
            if directive_node:
                directive_text = _get_node_text(directive_node, source_code).strip()
                if directive_text == "#undef":
                    results.append(_process_preproc_undef(node, source_code, filepath))
                elif directive_text == "#pragma":
                    results.append(_process_preproc_pragma(node, source_code, filepath))
                elif directive_text == "#error":
                    results.append(_process_preproc_error(node, source_code, filepath))
                else:
                    # Unknown preproc_call - treat as text
                    text_accumulator.append((node.start_byte, node.end_byte))
            else:
                text_accumulator.append((node.start_byte, node.end_byte))
        elif node.type in ["preproc_if", "preproc_ifdef", "preproc_ifndef"]:
            results.append(_process_conditional_group(node, source_code, filepath))
        else:
            # Unknown preproc type - treat as text
            text_accumulator.append((node.start_byte, node.end_byte))
    else:
        # Non-preprocessor node - check if it contains preproc nodes
        # Recursively process children to find nested preproc directives
        has_preproc_children = False
        child_text_accumulator: list[tuple[int, int]] = []
        
        for child in node.children:
            if child.type.startswith("preproc_") or (
                child.type == "preproc_call" and _find_child_by_type(child, "preproc_directive")
            ):
                has_preproc_children = True
                # Process any accumulated text before this preproc node
                if child_text_accumulator:
                    text_block = _merge_text_blocks(child_text_accumulator, source_code, filepath)
                    if text_block:
                        results.append(text_block)
                    child_text_accumulator = []
                
                # Process the preproc node
                processed = _process_node(child, source_code, filepath, [])
                if processed:
                    results.extend(processed)
            else:
                # Non-preproc child - recursively process it
                child_processed = _process_node(child, source_code, filepath, child_text_accumulator)
                if child_processed:
                    # If we have accumulated text, create a text block first
                    if child_text_accumulator:
                        text_block = _merge_text_blocks(child_text_accumulator, source_code, filepath)
                        if text_block:
                            results.append(text_block)
                        child_text_accumulator = []
                    results.extend(child_processed)
        
        if has_preproc_children:
            # If we had preproc children, merge any remaining text
            if child_text_accumulator:
                text_block = _merge_text_blocks(child_text_accumulator, source_code, filepath)
                if text_block:
                    results.append(text_block)
        else:
            # No preproc children - accumulate entire node for text block
            text_accumulator.append((node.start_byte, node.end_byte))
    
    return results


def _merge_text_blocks(
    text_accumulator: list[tuple[int, int]],
    source_code: bytes,
    filepath: str,
) -> Optional[TextBlock]:
    """Merge accumulated text ranges into a single TextBlock.
    
    Args:
        text_accumulator: List of (start_byte, end_byte) tuples
        source_code: The source code as bytes
        filepath: The file path
    
    Returns:
        A TextBlock if there's any text, None otherwise
    """
    if not text_accumulator:
        return None
    
    # Sort and merge overlapping or adjacent ranges
    sorted_ranges = sorted(text_accumulator)
    merged: list[tuple[int, int]] = []
    
    for start, end in sorted_ranges:
        if not merged:
            merged.append((start, end))
        else:
            last_start, last_end = merged[-1]
            if start <= last_end:
                # Overlapping or adjacent - merge
                merged[-1] = (last_start, max(last_end, end))
            else:
                merged.append((start, end))
    
    # Extract text from merged ranges
    text_parts = []
    for start, end in merged:
        text_parts.append(source_code[start:end].decode('utf-8'))
    
    content = "".join(text_parts)
    
    # Create span from first and last ranges
    first_start, _ = merged[0]
    _, last_end = merged[-1]
    
    # Create dummy nodes for span calculation
    class DummyNode:
        def __init__(self, start_byte: int, end_byte: int):
            self.start_byte = start_byte
            self.end_byte = end_byte
            # Calculate approximate positions
            text_before = source_code[:start_byte].decode('utf-8')
            lines_before = text_before.count('\n')
            col_before = len(text_before.split('\n')[-1])
            
            text_after = source_code[:end_byte].decode('utf-8')
            lines_after = text_after.count('\n')
            col_after = len(text_after.split('\n')[-1])
            
            self.start_point = type('Point', (), {'row': lines_before, 'column': col_before})()
            self.end_point = type('Point', (), {'row': lines_after, 'column': col_after})()
    
    start_node = DummyNode(first_start, first_start)
    end_node = DummyNode(last_end, last_end)
    
    return TextBlock(
        content=content,
        span=SourceSpan(
            filepath=filepath,
            start_line=start_node.start_point.row + 1,
            start_column=start_node.start_point.column,
            end_line=end_node.end_point.row + 1,
            end_column=end_node.end_point.column,
        ),
    )


def parse(source_code: str, filepath: str = "") -> FileRoot:
    """Parse C preprocessor directives from source code.
    
    This function walks the tree-sitter generated syntax tree, processes
    nodes starting with preproc_*, and merges other nodes into TextBlocks.
    Conditional structures are handled with stack-based processing.
    
    Args:
        source_code: The source code to parse
        filepath: Optional file path for SourceSpan objects
    
    Returns:
        A FileRoot containing the parsed structure
    """
    from tree_sitter import Language, Parser
    import tree_sitter_fortran
    
    # Initialize parser
    language = Language(tree_sitter_fortran.language())
    parser = Parser(language)
    
    # Parse source code
    source_bytes = source_code.encode('utf-8')
    tree = parser.parse(source_bytes)
    
    # Process root node's children
    items: list[ModelNode] = []
    text_accumulator: list[tuple[int, int]] = []
    
    for child in tree.root_node.children:
        processed = _process_node(child, source_bytes, filepath, text_accumulator)
        if processed:
            # If we have accumulated text, create a text block first
            if text_accumulator:
                text_block = _merge_text_blocks(text_accumulator, source_bytes, filepath)
                if text_block:
                    items.append(text_block)
                text_accumulator = []
            
            items.extend(processed)
    
    # Handle any remaining text
    if text_accumulator:
        text_block = _merge_text_blocks(text_accumulator, source_bytes, filepath)
        if text_block:
            items.append(text_block)
    
    # Create FileRoot
    root_span = _node_to_span(tree.root_node, filepath)
    
    return FileRoot(
        path=filepath,
        items=items,
        span=root_span,
    )

