# cpptree.examples Usage Guide

The `cpptree.examples` module provides functionality to access example Fortran files and their tree-sitter parsing results.

## Basic Usage

### List all available examples

```python
from cpptree.examples import list_examples

# Get all example file names
examples = list_examples()
print(examples)
# Output: ['01_simple_defines.f90', '02_function_like_defines.f90', ...]
```

### Get example source code and AST

```python
from cpptree.examples import get_example

# Get example file (automatically parsed)
source_code, tree = get_example("01_simple_defines.f90")

# Or without extension is also fine
source_code, tree = get_example("01_simple_defines")

# source_code is a string containing file content
print(f"File length: {len(source_code)} characters")

# tree is a tree-sitter Tree object
print(f"Root node type: {tree.root_node.type}")
print(f"Number of children: {len(tree.root_node.children)}")
```

### Get example file path

```python
from cpptree.examples import get_example_path

# Get file path (without parsing)
path = get_example_path("01_simple_defines.f90")
print(f"File path: {path}")
```

## Complete Example

```python
from cpptree.examples import list_examples, get_example

# 1. List all examples
print("Available example files:")
for name in list_examples():
    print(f"  - {name}")

# 2. Parse an example file
source_code, tree = get_example("04_conditional_if.f90")

# 3. View source code
print("\nSource code:")
print(source_code[:500])  # First 500 characters

# 4. Traverse AST
def print_ast(node, indent=0):
    """Recursively print AST nodes"""
    print("  " * indent + node.type)
    for child in node.children:
        print_ast(child, indent + 1)

print("\nAST structure:")
print_ast(tree.root_node)

# 5. Find nodes of specific type
def find_nodes(node, node_type, results=None):
    """Find all nodes of the specified type"""
    if results is None:
        results = []
    if node.type == node_type:
        results.append(node)
    for child in node.children:
        find_nodes(child, node_type, results)
    return results

# Find all comment nodes
comments = find_nodes(tree.root_node, "comment")
print(f"\nFound {len(comments)} comment nodes")
```

## API Reference

### `list_examples() -> list[str]`

Returns a list of all available example file names (excluding paths, only file names).

**Returns**: List of example file names

**Example**:
```python
examples = list_examples()
# ['01_simple_defines.f90', '02_function_like_defines.f90', ...]
```

### `get_example(name: str) -> tuple[str, Tree]`

Get the source code and parsed AST of the specified example file.

**Parameters**:
- `name`: Example file name (e.g., "01_simple_defines.f90"). If only the file name without extension is provided, .f90 will be automatically added.

**Returns**:
- A tuple `(source_code, tree)`, where:
  - `source_code`: The source code content of the file (string)
  - `tree`: The parsed Tree object from tree-sitter

**Exceptions**:
- `FileNotFoundError`: If the specified example file does not exist
- `ValueError`: If the file name is invalid

**Example**:
```python
source_code, tree = get_example("01_simple_defines.f90")
```

### `get_example_path(name: str) -> Path`

Get the path of the example file (without parsing).

**Parameters**:
- `name`: Example file name

**Returns**:
- Path object of the example file

**Exceptions**:
- `FileNotFoundError`: If the specified example file does not exist

**Example**:
```python
path = get_example_path("01_simple_defines.f90")
print(path)  # /path/to/cpptree/src/cpptree/examples/01_simple_defines.f90
```

## Notes

1. **Parser initialization**: `get_example()` automatically initializes the tree-sitter-fortran parser (using singleton pattern), and the first call may have a slight delay.

2. **File location**: Example files are located in the `cpptree/src/cpptree/examples/` directory.

3. **Encoding**: All files are read using UTF-8 encoding.

4. **Tree object**: The returned `Tree` object is a native tree-sitter object and can be manipulated using all tree-sitter APIs.
