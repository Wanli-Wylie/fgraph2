本文档描述了 C 预处理器（CPP）抽象语法树（AST）的数据模型定义。该模型基于 Python 的 pydantic 库构建，利用了强类型系统和不可变性（Immutability）来确保 AST 的完整性。

该模型旨在捕获 C 源代码中的预处理指令（如 `#include`, `#define`, `#if`）以及指令之间的普通文本块。

## 1. 模型概览 (Overview)

### 核心特性

- 不可变性 (Immutability): 所有模型配置为 `frozen=True`，确保 AST 构建后不会被意外修改。
- 区分联合 (Discriminated Union): 通过 `kind` 字段作为标识符（Discriminator），实现了多态类型的自动解析。
- 源码映射: 几乎所有节点都包含 `span: SourceSpan`，用于追踪节点在原始文件中的位置（行号/列号）。
- 递归结构: 条件编译组（Conditional Groups）支持递归嵌套。

### 类型层次图 (Class Diagram)

```
classDiagram
    class FileRoot {
        +path: str
        +items: list[PreprocessorNode]
    }

    class PreprocessorNode {
        <<Union>>
    }

    class ConditionalGroup {
        +kind: "conditional_group"
        +entry: ConditionalBranch
        +elifs: list[ConditionalBranch]
        +else_body: list[PreprocessorNode]
    }

    class ConditionalBranch {
        +kind: if/ifdef/ifndef/elif
        +condition: str
        +body: list[PreprocessorNode]
    }

    FileRoot --> PreprocessorNode : contains
    ConditionalBranch --> PreprocessorNode : contains (recursive)
    ConditionalGroup --> ConditionalBranch : consists of

    PreprocessorNode <|-- TextBlock
    PreprocessorNode <|-- DefineNode
    PreprocessorNode <|-- FunctionDefineNode
    PreprocessorNode <|-- IncludeNode
    PreprocessorNode <|-- UnDefNode
    PreprocessorNode <|-- PragmaNode
    PreprocessorNode <|-- ErrorNode
    PreprocessorNode <|-- ConditionalGroup
```

------

## 2. 基础内容节点 (Basic Content Nodes)

这些节点代表了源代码中最基础的组成部分：普通文本和简单的宏定义。

### `TextBlock`

表示非预处理指令的内容。这包括普通的 C 代码、注释或指令之间的空白。

| 字段  | 类型          | 描述         |
| --------- | ----------------- | ---------------- |
| `kind`    | `Literal["text"]` | 固定标识符。     |
| `content` | `str`             | 实际的文本内容。 |
| `span`    | `SourceSpan`      | 源码位置信息。   |

### `DefineNode`

表示对象式宏定义（Object-like Macro），即不带参数的宏。

- 示例: `#define PI 3.14159`

| 字段 | 类型            | 描述                 |
| -------- | ------------------- | ------------------------ |
| `kind`   | `Literal["define"]` | 固定标识符。             |
| `name`   | `str`               | 宏名称（如 `PI`）。      |
| `value`  | `str`               | 宏的值（如 `3.14159`）。 |
| `raw`    | `str`               | 原始指令行文本。         |

### `FunctionDefineNode`

表示函数式宏定义（Function-like Macro），即带参数的宏。

- 示例: `#define MIN(a, b) ((a) < (b) ? (a) : (b))`

| 字段     | 类型                  | 描述                      |
| ------------ | ------------------------- | ----------------------------- |
| `kind`       | `Literal["function_def"]` | 固定标识符。                  |
| `name`       | `str`                     | 宏名称。                      |
| `value`      | `str`                     | 宏的展开体。                  |
| `parameters` | `list[str]`               | 参数列表（如 `['a', 'b']`）。 |
| `raw`        | `str`                     | 原始指令行文本。              |

### `UnDefNode`

表示取消宏定义的指令。

- 示例: `#undef DEBUG`

| 字段 | 类型           | 描述             |
| -------- | ------------------ | -------------------- |
| `kind`   | `Literal["undef"]` | 固定标识符。         |
| `name`   | `str`              | 被取消定义的宏名称。 |

------

## 3. 指令与控制节点 (Directive Nodes)

### `IncludeNode`

表示文件包含指令。

| 字段 | 类型             | 描述                                                   |
| -------- | -------------------- | ---------------------------------------------------------- |
| `kind`   | `Literal["include"]` | 固定标识符。                                               |
| `target` | `str`                | 包含的目标路径或文件名（如 `<stdio.h>` 或 `"my_lib.h"`）。 |
| `raw`    | `str`                | 原始指令行。                                               |

### `PragmaNode`

表示编译器特定的 Pragma 指令。

- 注意：该节点不强制包含 `span` (根据当前定义)。

| 字段  | 类型            | 描述                                      |
| --------- | ------------------- | --------------------------------------------- |
| `kind`    | `Literal["pragma"]` | 固定标识符。                                  |
| `payload` | `str`               | Pragma 的具体内容（如 `once` 或 `pack(1)`）。 |

### `ErrorNode`

表示报错指令，通常用于在预处理阶段中断编译。

| 字段  | 类型           | 描述       |
| --------- | ------------------ | -------------- |
| `kind`    | `Literal["error"]` | 固定标识符。   |
| `message` | `str`              | 错误消息内容。 |

------

## 4. 条件编译结构 (Conditional Compilation)

条件编译的处理采用了嵌套结构，将 `#if`...`#endif` 块视为一个单一的逻辑单元 (`ConditionalGroup`)。

### 类型定义: `PreprocessorNode`

这是一个 Discriminated Union 类型，包含了所有可能的节点类型。用于在列表中存储异构节点。

### `ConditionalBranch`

这不是一个独立的顶层节点，而是 `ConditionalGroup` 的组成部分。它表示条件编译链中的一个分支（`if`, `ifdef`, `ifndef`, 或 `elif`）。

| 字段    | 类型                 | 描述                                                     |
| ----------- | ------------------------ | ------------------------------------------------------------ |
| `kind`      | `Literal`                | 分支类型: `"if"`, `"ifdef"`, `"ifndef"`, `"elif"`。          |
| `condition` | `str`                    | 条件表达式（如 `DEBUG_LEVEL > 0`）。                         |
| `body`      | `list[PreprocessorNode]` | 递归字段。该分支下的内容列表，可能包含嵌套的 `ConditionalGroup`。 |
| `raw`       | `str`                    | 分支开始的原始行文本。                                       |

### `ConditionalGroup`

表示完整的条件编译块。它封装了从 `#if` 开始到 `#endif` 结束的完整逻辑。

- 结构逻辑: `Entry Branch` -> `Optional Elifs` -> `Optional Else` -> `Endif`

| 字段    | 类型                            | 描述                                          |
| ----------- | ----------------------------------- | ------------------------------------------------- |
| `kind`      | `Literal["conditional_group"]`      | 固定标识符。                                      |
| `entry`     | `ConditionalBranch`                 | 入口分支（必须是 `#if`, `#ifdef` 或 `#ifndef`）。 |
| `elifs`     | `Optional[list[ConditionalBranch]]` | 零个或多个 `#elif` 分支。                         |
| `else_body` | `Optional[list[PreprocessorNode]]`  | `#else` 分支的内容列表（如果存在）。              |
| `else_raw`  | `Optional[str]`                     | `#else` 指令的原始文本。                          |
| `endif_raw` | `str`                               | `#endif` 指令的原始文本（默认为 `"#endif"`）。    |

------

## 5. 根节点 (Root Node)

### `FileRoot`

整个解析过程的顶层容器。代表一个完整的源文件。

| 字段 | 类型                 | 描述                 |
| -------- | ------------------------ | ------------------------ |
| `path`   | `str`                    | 文件路径。               |
| `items`  | `list[PreprocessorNode]` | 文件内的顶级节点列表。   |
| `span`   | `SourceSpan`             | 通常表示整个文件的跨度。 |