from cpptree.core.parse import process_tree_stream, _node_to_span
from cpptree.models import FileRoot
from tree_sitter import Language, Parser
import tree_sitter_fortran

def parse(source_code: str, filepath: str = "") -> FileRoot:
    """Parse C preprocessor directives from source code.
    
    This function uses a Visitor pattern with generator-based stream processing
    to walk the tree-sitter generated syntax tree, process nodes starting with
    preproc_*, and merge other nodes into TextBlocks.
    
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
    
    # Process tree using generator-based stream processing
    items = list(process_tree_stream(tree.root_node, source_bytes, filepath))
    
    # Create FileRoot
    root_span = _node_to_span(tree.root_node, filepath)
    
    return FileRoot(
        path=filepath,
        items=items,
        span=root_span,
    )

