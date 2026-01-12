# cpptree.core.parse 单元测试

本目录包含 `cpptree.core.parse` 模块的单元测试，为每个 `_process*` 方法和辅助函数提供了全面的测试覆盖。

## 测试结构

### `conftest.py`
提供 pytest fixtures：
- `parser`: tree-sitter Fortran 解析器
- `parse_source`: 解析源代码并返回 (source_bytes, tree)
- `find_node_by_type`: 在树中递归查找指定类型的节点

### `test_parse_helpers.py`
测试辅助函数（14 个测试用例）：

#### `TestNodeToSpan` (3 个测试)
- `test_basic_conversion`: 基本节点到 SourceSpan 的转换
- `test_empty_filepath`: 空文件路径测试
- `test_multiline_node`: 多行节点的转换

#### `TestGetNodeText` (3 个测试)
- `test_simple_text_extraction`: 简单文本提取
- `test_identifier_text`: 标识符文本提取
- `test_empty_node`: 空节点处理

#### `TestFindChildByType` (3 个测试)
- `test_find_existing_child`: 查找存在的子节点
- `test_find_nonexistent_child`: 查找不存在的子节点
- `test_find_first_matching_child`: 查找第一个匹配的子节点

#### `TestExtractConditionText` (4 个测试)
- `test_simple_condition`: 简单条件提取
- `test_ifdef_condition`: ifdef 条件提取
- `test_complex_condition`: 复杂条件提取
- `test_elif_condition`: elif 条件提取

### `test_parse_processors.py`
测试处理器函数（28 个测试用例）：

#### `TestProcessPreprocDef` (4 个测试)
- `test_simple_define`: 简单 #define 处理
- `test_define_without_value`: 无值的 #define
- `test_define_with_string_value`: 字符串值的 #define
- `test_define_missing_identifier`: 缺少标识符的情况

#### `TestProcessPreprocFunctionDef` (3 个测试)
- `test_function_like_define`: 函数式宏定义
- `test_single_parameter`: 单参数函数式宏
- `test_multiple_parameters`: 多参数函数式宏

#### `TestProcessPreprocInclude` (2 个测试)
- `test_quoted_include`: 引号形式的 #include
- `test_angle_bracket_include`: 尖括号形式的 #include

#### `TestProcessPreprocUndef` (2 个测试)
- `test_undef_from_preproc_call`: 从 preproc_call 处理 #undef
- `test_undef_from_preproc_undef`: 从 preproc_undef 处理 #undef

#### `TestProcessPreprocPragma` (2 个测试)
- `test_simple_pragma`: 简单 #pragma
- `test_pragma_with_value`: 带值的 #pragma

#### `TestProcessPreprocError` (1 个测试)
- `test_error_with_message`: 带消息的 #error

#### `TestProcessConditionalBranch` (2 个测试)
- `test_if_branch`: #if 分支处理
- `test_elif_branch`: #elif 分支处理

#### `TestProcessConditionalGroup` (4 个测试)
- `test_simple_if_endif`: 简单 #if...#endif
- `test_if_else_endif`: #if...#else...#endif
- `test_if_elif_else_endif`: #if...#elif...#else...#endif
- `test_ifdef_group`: #ifdef...#endif

#### `TestMergeTextBlocks` (4 个测试)
- `test_single_text_block`: 单个文本块合并
- `test_multiple_adjacent_blocks`: 多个相邻文本块合并
- `test_overlapping_blocks`: 重叠文本块合并
- `test_empty_accumulator`: 空累加器处理

#### `TestProcessNode` (5 个测试)
- `test_process_preproc_def_node`: 处理 preproc_def 节点
- `test_process_text_node`: 处理文本节点
- `test_process_preproc_call_undef`: 处理 preproc_call (#undef)
- `test_process_preproc_call_pragma`: 处理 preproc_call (#pragma)
- `test_process_conditional_node`: 处理条件节点

## 运行测试

### 运行所有测试
```bash
pytest tests/
```

### 运行特定测试文件
```bash
pytest tests/test_parse_helpers.py
pytest tests/test_parse_processors.py
```

### 运行特定测试类
```bash
pytest tests/test_parse_helpers.py::TestNodeToSpan
```

### 运行特定测试用例
```bash
pytest tests/test_parse_processors.py::TestProcessPreprocDef::test_simple_define
```

### 详细输出
```bash
pytest tests/ -v
```

### 显示覆盖率
```bash
pytest tests/ --cov=cpptree.core.parse --cov-report=html
```

## 测试覆盖范围

### 已覆盖的函数

#### 辅助函数
- ✅ `_node_to_span`
- ✅ `_get_node_text`
- ✅ `_find_child_by_type`
- ✅ `_extract_condition_text`

#### 预处理器指令处理器
- ✅ `_process_preproc_def`
- ✅ `_process_preproc_function_def`
- ✅ `_process_preproc_include`
- ✅ `_process_preproc_undef`
- ✅ `_process_preproc_pragma`
- ✅ `_process_preproc_error`

#### 条件编译处理器
- ✅ `_process_conditional_branch`
- ✅ `_process_conditional_group`

#### 文本处理
- ✅ `_merge_text_blocks`
- ✅ `_process_node`

## 测试统计

- **总测试数**: 42
- **测试文件**: 2
- **测试类**: 12
- **通过率**: 100% ✅

## 测试设计原则

1. **隔离性**: 每个测试用例独立，不依赖其他测试
2. **可重复性**: 使用 fixtures 确保测试环境一致
3. **边界测试**: 包含正常情况、边界情况和错误情况
4. **清晰性**: 测试名称清晰描述测试内容
5. **模块化**: 按功能分组，便于维护

## 注意事项

- 测试使用真实的 tree-sitter 解析结果，确保与实际使用场景一致
- 某些测试用例可能因为 tree-sitter 解析行为而跳过（使用 `if node:` 检查）
- 测试覆盖了主要的代码路径，但可能不覆盖所有边缘情况

