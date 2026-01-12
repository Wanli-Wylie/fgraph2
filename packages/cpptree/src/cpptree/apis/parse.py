from cpptree.core.parse import _process_node, _merge_text_blocks, _node_to_span
from cpptree.models import FileRoot, PreprocessorNode
from tree_sitter import Language, Parser
import tree_sitter_fortran

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
    
    
    # Initialize parser
    language = Language(tree_sitter_fortran.language())
    parser = Parser(language)
    
    # Parse source code
    source_bytes = source_code.encode('utf-8')
    tree = parser.parse(source_bytes)
    
    # Process root node's children
    items: list[PreprocessorNode] = []
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

