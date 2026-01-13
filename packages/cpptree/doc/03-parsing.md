# 解析模块文档 (Parsing Module Documentation)

本文档描述了 `cpptree.core.parse` 模块的实现细节。该模块负责将 tree-sitter 生成的语法树转换为 C 预处理器抽象语法树（AST）的数据模型。

## 1. 模块概览 (Overview)

### 核心职责

`cpptree.core.parse` 模块是解析流程的核心实现层，它：

- **节点转换**: 将 tree-sitter 的 `Node` 对象转换为 `PreprocessorNode` 数据模型
- **递归处理**: 遍历语法树，识别并处理所有预处理指令
- **文本合并**: 将非预处理指令的内容合并为 `TextBlock` 节点
- **位置追踪**: 为每个节点生成 `SourceSpan`，记录其在源文件中的位置

### 架构设计

```
tree-sitter AST (Node)
    ↓
_process_node() [递归处理]
    ↓
特定处理器函数 (_process_preproc_*)
    ↓
数据模型 (PreprocessorNode)
```

### 依赖关系

- **输入**: tree-sitter 的 `Node` 对象和源代码字节流
- **输出**: `PreprocessorNode` 及其子类型的列表
- **依赖**: `cpptree.models` 中定义的数据模型

---

## 2. 辅助函数 (Helper Functions)

### `_node_to_span(node: Node, filepath: str = "") -> SourceSpan`

将 tree-sitter 节点转换为 `SourceSpan` 对象，用于记录源码位置信息。

**参数**:
- `node`: tree-sitter 节点对象
- `filepath`: 文件路径（默认为空字符串）

**返回**: `SourceSpan` 对象

**实现细节**:
- tree-sitter 使用 0 索引的行号，而 `SourceSpan` 使用 1 索引
- 自动转换行号索引：`start_line = start_point.row + 1`

**示例**:
```python
span = _node_to_span(node, "/path/to/file.c")
# span.start_line 为 1-indexed
```

### `_get_node_text(node: Node, source_code: bytes) -> str`

从 tree-sitter 节点中提取文本内容。

**参数**:
- `node`: tree-sitter 节点对象
- `source_code`: 源代码的字节流

**返回**: 节点对应的文本字符串（UTF-8 解码）

**实现细节**:
- 使用节点的 `start_byte` 和 `end_byte` 属性进行字节切片
- 自动进行 UTF-8 解码

### `_find_child_by_type(node: Node, node_type: str) -> Optional[Node]`

查找指定类型的第一个子节点。

**参数**:
- `node`: 父节点
- `node_type`: 要查找的子节点类型（字符串）

**返回**: 第一个匹配的子节点，如果未找到则返回 `None`

**用途**: 用于从预处理器节点中提取特定信息（如标识符、参数等）

---

## 3. 预处理器指令处理器 (Preprocessor Directive Processors)

### `_process_preproc_def(node, source_code, filepath) -> DefineNode`

处理简单的 `#define` 指令（对象式宏）。

**处理的节点类型**: `preproc_def`

**提取的信息**:
- `name`: 宏名称（从 `identifier` 子节点提取）
- `value`: 宏的值（从 `preproc_arg` 子节点提取，如果存在）
- `raw`: 完整的原始指令文本
- `span`: 源码位置信息

**错误处理**: 如果缺少 `identifier` 子节点，抛出 `ValueError`

**示例输入**:
```c
#define PI 3.14159
```

**示例输出**: `DefineNode(name="PI", value="3.14159", ...)`

### `_process_preproc_function_def(node, source_code, filepath) -> FunctionDefineNode`

处理函数式 `#define` 指令（函数式宏）。

**处理的节点类型**: `preproc_function_def`

**提取的信息**:
- `name`: 宏名称
- `parameters`: 参数列表（从 `preproc_params` 中提取所有 `identifier` 子节点）
- `value`: 宏的展开体（从 `preproc_arg` 提取）
- `raw`: 完整的原始指令文本
- `span`: 源码位置信息

**示例输入**:
```c
#define MAX(a, b) ((a) > (b) ? (a) : (b))
```

**示例输出**: `FunctionDefineNode(name="MAX", parameters=["a", "b"], value="((a) > (b) ? (a) : (b))", ...)`

### `_process_preproc_include(node, source_code, filepath) -> IncludeNode`

处理 `#include` 指令。

**处理的节点类型**: `preproc_include`

**提取的信息**:
- `target`: 包含的目标路径（去除引号或尖括号）
- `raw`: 完整的原始指令文本
- `span`: 源码位置信息

**实现细节**:
- 支持三种包含格式：
  - `string_literal`: `#include "header.h"`
  - `system_lib_string`: `#include <header.h>`
  - `system_lib_name`: 系统库名称
- 自动去除引号和尖括号：`.strip('"<>')`

**错误处理**: 如果找不到目标节点，抛出 `ValueError`

### `_process_preproc_undef(node, source_code, filepath) -> UnDefNode`

处理 `#undef` 指令。

**处理的节点类型**: `preproc_undef` 或 `preproc_call`（当指令为 `#undef` 时）

**提取的信息**:
- `name`: 要取消定义的宏名称
- `raw`: 完整的原始指令文本
- `span`: 源码位置信息

**实现细节**:
- 优先从 `preproc_arg` 子节点提取名称
- 如果不存在，则从 `identifier` 子节点提取
- 自动去除空白字符

### `_process_preproc_pragma(node, source_code, filepath) -> PragmaNode`

处理 `#pragma` 指令。

**处理的节点类型**: `preproc_pragma` 或 `preproc_call`（当指令为 `#pragma` 时）

**提取的信息**:
- `payload`: pragma 的具体内容（去除 `#pragma` 关键字后的文本）
- `raw`: 完整的原始指令文本
- `span`: 源码位置信息

**实现细节**:
- 优先从 `preproc_arg` 子节点提取
- 如果不存在，则从原始文本中手动提取（去除 `#pragma` 前缀）

### `_process_preproc_error(node, source_code, filepath) -> ErrorNode`

处理 `#error` 指令。

**处理的节点类型**: `preproc_error` 或 `preproc_call`（当指令为 `#error` 时）

**提取的信息**:
- `message`: 错误消息内容
- `raw`: 完整的原始指令文本
- `span`: 源码位置信息

**实现细节**:
- 优先从 `preproc_arg` 子节点提取消息
- 如果不存在，则从原始文本中手动提取（去除 `#error` 前缀和引号）

---

## 4. 条件编译处理 (Conditional Compilation Processing)

### `_extract_condition_text(node, source_code) -> str`

从条件编译节点中提取条件表达式文本。

**处理的节点类型**: `preproc_if`, `preproc_ifdef`, `preproc_ifndef`, `preproc_elif`

**实现策略**:
1. 遍历子节点，找到指令标记（`#if`, `#ifdef`, `#ifndef`, `#elif`）
2. 提取指令后的第一个非空白、非注释子节点作为条件
3. 如果上述方法失败，则从原始文本中手动提取（去除指令关键字，取第一行）

**返回**: 条件表达式字符串

### `_process_conditional_branch(node, source_code, filepath, branch_kind) -> tuple[ConditionalBranch, list[Node]]`

处理单个条件分支（`#if`, `#ifdef`, `#ifndef`, `#elif`）。

**参数**:
- `node`: 条件节点
- `source_code`: 源代码字节流
- `filepath`: 文件路径
- `branch_kind`: 分支类型（`"if"`, `"ifdef"`, `"ifndef"`, `"elif"`）

**返回**: `(ConditionalBranch, remaining_children)` 元组
- `ConditionalBranch`: 处理后的分支对象
- `remaining_children`: 剩余的未处理子节点（如 `#elif`, `#else`, `#endif`）

**处理逻辑**:
1. 提取条件表达式
2. 收集分支体内容（从指令后到下一个预处理器指令之间的所有节点）
3. 递归处理分支体内的所有节点
4. 返回分支对象和剩余节点

### `_process_conditional_group(node, source_code, filepath) -> ConditionalGroup`

处理完整的条件编译组（从 `#if` 到 `#endif`）。

**处理的节点类型**: `preproc_if`, `preproc_ifdef`, `preproc_ifndef`

**处理流程**:
1. **确定入口分支类型**: 根据节点类型确定是 `if`、`ifdef` 还是 `ifndef`
2. **处理入口分支**: 调用 `_process_conditional_branch` 处理第一个分支
3. **处理 `#elif` 分支**: 循环处理所有 `#elif` 分支
4. **处理 `#else` 分支**: 如果存在，收集 `#else` 到 `#endif` 之间的内容
5. **计算整体跨度**: 从入口分支到 `#endif` 的完整 `SourceSpan`

**返回**: `ConditionalGroup` 对象，包含：
- `entry`: 入口分支（`ConditionalBranch`）
- `elifs`: 可选的 `#elif` 分支列表
- `else_body`: 可选的 `#else` 分支内容
- `else_raw`: `#else` 指令的原始文本
- `endif_raw`: `#endif` 指令的原始文本
- `span`: 整个条件组的源码位置

**嵌套支持**: 条件组内部可以包含嵌套的条件组，通过递归调用 `_process_node` 实现。

---

## 5. 核心处理函数 (Core Processing Function)

### `_process_node(node, source_code, filepath, text_accumulator) -> list[ModelNode]`

核心递归处理函数，负责将 tree-sitter 节点转换为数据模型节点。

**参数**:
- `node`: 要处理的 tree-sitter 节点
- `source_code`: 源代码字节流
- `filepath`: 文件路径
- `text_accumulator`: 文本累积器，用于收集非预处理指令的文本范围 `(start_byte, end_byte)`

**返回**: `PreprocessorNode` 列表

**处理逻辑**:

#### 情况 1: 预处理器节点（`node.type.startswith("preproc_")`）

根据节点类型分发到相应的处理器：

| 节点类型 | 处理器函数 | 输出类型 |
|---------|-----------|---------|
| `preproc_def` | `_process_preproc_def` | `DefineNode` |
| `preproc_function_def` | `_process_preproc_function_def` | `FunctionDefineNode` |
| `preproc_include` | `_process_preproc_include` | `IncludeNode` |
| `preproc_undef` | `_process_preproc_undef` | `UnDefNode` |
| `preproc_pragma` | `_process_preproc_pragma` | `PragmaNode` |
| `preproc_error` | `_process_preproc_error` | `ErrorNode` |
| `preproc_call` | 根据指令类型分发 | `UnDefNode` / `PragmaNode` / `ErrorNode` |
| `preproc_if` / `preproc_ifdef` / `preproc_ifndef` | `_process_conditional_group` | `ConditionalGroup` |

**特殊处理**: `preproc_call` 节点需要检查其 `preproc_directive` 子节点来确定具体指令类型。

#### 情况 2: 非预处理器节点

递归处理子节点，寻找嵌套的预处理器指令：

1. **遍历所有子节点**
2. **识别预处理器子节点**: 如果子节点是预处理器节点，则：
   - 先处理之前累积的文本（创建 `TextBlock`）
   - 递归处理预处理器节点
   - 清空文本累积器
3. **处理非预处理器子节点**: 递归调用 `_process_node`，传递文本累积器
4. **文本累积**: 如果没有预处理器子节点，将整个节点范围添加到文本累积器

**文本块生成**: 当遇到预处理器节点时，会先处理累积的文本，确保文本块和预处理器节点按顺序输出。

---

## 6. 文本块合并 (Text Block Merging)

### `_merge_text_blocks(text_accumulator, source_code, filepath) -> Optional[TextBlock]`

将累积的文本范围合并为单个 `TextBlock` 节点。

**参数**:
- `text_accumulator`: `(start_byte, end_byte)` 元组列表
- `source_code`: 源代码字节流
- `filepath`: 文件路径

**返回**: `TextBlock` 对象，如果没有文本则返回 `None`

**合并算法**:
1. **排序**: 按 `start_byte` 对范围进行排序
2. **合并重叠或相邻的范围**: 
   - 如果当前范围的起始位置 ≤ 上一个范围的结束位置，则合并
   - 合并后的结束位置取两个范围的最大值
3. **提取文本**: 从合并后的范围中提取文本并拼接
4. **计算跨度**: 使用第一个和最后一个范围计算 `SourceSpan`

**实现细节**:
- 使用 `DummyNode` 类来模拟 tree-sitter 节点，用于计算行号和列号
- 通过计算换行符数量来确定行号
- 通过计算最后一行的字符数来确定列号

**用途**: 确保连续的文本内容被合并为单个 `TextBlock`，避免碎片化。

---

## 7. 解析流程示例 (Parsing Flow Example)

以下是一个完整的解析流程示例：

### 输入源代码
```c
#include <stdio.h>
#define PI 3.14159
int main() {
    return 0;
}
```

### 处理步骤

1. **tree-sitter 解析**: 生成语法树，根节点包含多个子节点
2. **遍历根节点的子节点**:
   - 子节点 1: `preproc_include` → `_process_preproc_include` → `IncludeNode`
   - 子节点 2: `preproc_def` → `_process_preproc_def` → `DefineNode`
   - 子节点 3: `function_definition` → 递归处理，发现无预处理器子节点 → 添加到文本累积器
3. **文本块生成**: 将累积的文本范围合并为 `TextBlock`
4. **输出**: `FileRoot(items=[IncludeNode(...), DefineNode(...), TextBlock(...)])`

### 条件编译示例

**输入**:
```c
#if DEBUG
    #define LOG(msg) printf("%s\n", msg)
#else
    #define LOG(msg)
#endif
```

**处理**:
1. 识别 `preproc_if` 节点
2. 调用 `_process_conditional_group`:
   - 处理入口分支（`#if DEBUG`）
   - 处理 `#else` 分支
   - 计算整体跨度
3. 输出: `ConditionalGroup` 包含两个 `ConditionalBranch`（entry 和 else_body）

---

## 8. 错误处理 (Error Handling)

### 常见错误情况

1. **缺少必需子节点**: 
   - 例如 `preproc_def` 缺少 `identifier` 子节点
   - 抛出 `ValueError` 异常，包含原始文本信息

2. **未知的预处理器类型**:
   - 遇到未实现的 `preproc_*` 类型
   - 将其作为文本处理，添加到文本累积器

3. **格式错误的指令**:
   - 例如 `#include` 缺少目标文件
   - 抛出 `ValueError` 异常

### 容错机制

- **文本累积**: 未知或格式错误的节点会被作为文本处理，不会中断整个解析过程
- **条件提取回退**: 如果无法从 AST 中提取条件表达式，会尝试从原始文本中提取

---

## 9. 性能考虑 (Performance Considerations)

### 优化策略

1. **字节级操作**: 使用 `start_byte` 和 `end_byte` 进行字节切片，避免字符串操作
2. **延迟文本提取**: 只在需要时提取文本内容，减少不必要的解码操作
3. **范围合并**: 合并相邻的文本范围，减少 `TextBlock` 节点数量

### 内存使用

- 源代码以字节流形式传递，避免多次编码/解码
- 文本累积器使用轻量级的 `(int, int)` 元组，而非字符串对象

---

## 10. 扩展性 (Extensibility)

### 添加新的预处理器指令支持

要添加对新指令的支持，需要：

1. **在 `models.py` 中定义新的节点类型**
2. **实现对应的 `_process_preproc_*` 函数**
3. **在 `_process_node` 中添加分发逻辑**

**示例**（伪代码）:
```python
def _process_preproc_warning(node, source_code, filepath) -> WarningNode:
    # 处理 #warning 指令
    ...

# 在 _process_node 中添加:
elif node.type == "preproc_warning":
    results.append(_process_preproc_warning(node, source_code, filepath))
```

---

## 11. 测试覆盖 (Test Coverage)

该模块的测试位于 `tests/` 目录：

- **`test_parse_helpers.py`**: 测试辅助函数（`_node_to_span`, `_get_node_text`, `_find_child_by_type`, `_extract_condition_text`）
- **`test_parse_processors.py`**: 测试各个处理器函数（`_process_preproc_*`）

详细的测试文档请参考 `tests/README.md`。

