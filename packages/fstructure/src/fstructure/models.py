from __future__ import annotations

from typing import Literal, Union, Optional
from pydantic import BaseModel, ConfigDict, Field
from fgraphutils import CodeSnippet

# ==========================================
# 1. 基础节点
# ==========================================

class FortranNode(BaseModel):
    """所有节点的基类"""
    kind: str
    snippet: CodeSnippet
    model_config = ConfigDict(frozen=True, extra="forbid")

class SpecCodeBlock(FortranNode):
    """
    规范部代码块
    """
    kind: Literal["spec_code_block"] = "spec_code_block"

class ExecCodeBlock(FortranNode):
    """
    可执行代码块
    """
    kind: Literal["exec_code_block"] = "exec_code_block"


# ==========================================
# 3. 递归容器定义 (关键：显式控制流)
# ==========================================

# 预声明，解决 Pydantic 递归引用
class ExecutionPart(BaseModel):
    """可执行部分的容器，代表一个语句序列"""
    content: list[ExecCodeBlock] = Field(default_factory=list)

# --- 3.1 控制流构造 (Constructs) ---
# 它们决定了代码的“形状”和“路径”，是承载作用域的骨架

# ==========================================
# 4. 聚合类型 (Union)
# ==========================================

SpecChild = Union[
    SpecCodeBlock,              # 声明语句
    "InterfaceBlockScope"     # 接口块容器
]

class RenameRule(BaseModel):
    local_name: str
    use_name: str

class UseStatement(FortranNode): # 建议继承 FortranNode 以保留位置信息
    kind: Literal["use_statement"] = "use_statement"
    module_name: str
    # nature: Literal["intrinsic", "non_intrinsic", "unspecified"] = "unspecified" # 可选：F03特性
    only_list: Optional[list[str]] = None  # None 表示导入所有；[] 表示不导入符号(仅为了依赖)
    renames: list[RenameRule] = Field(default_factory=list)

class ImportStatement(FortranNode):
    kind: Literal["import_statement"] = "import_statement"
    names: Optional[list[str]] = None # None 表示 IMPORT :: ALL (F2018), F03 通常有列表

# ==========================================
# 5. 顶层结构 (Subprograms / Modules)
# ==========================================

# --- 接口结构 ---
class InterfaceBodyScope(FortranNode):
    kind: Literal["interface_body_scope"] = "interface_body_scope"
    # NOTE: IMPORT statements
    import_statements: Optional[list[ImportStatement]] = Field(default_factory=list)
    spec_part: list[SpecChild] = Field(default_factory=list)

class InterfaceBlockScope(FortranNode):
    """接口块是容器，持有 Interface Bodies"""
    kind: Literal["interface_block_scope"] = "interface_block_scope"

    module_procedures: Optional[SpecCodeBlock] = None
    bodies: list[InterfaceBodyScope] = Field(default_factory=list)

# --- 程序单元结构 ---
class InternalSubprogramScope(FortranNode):
    kind: Literal["internal_subprogram_scope"] = "internal_subprogram_scope"
    # NOTE: Use statements
    use_statements: Optional[list[UseStatement]] = Field(default_factory=list)
    spec_part: list[SpecChild] = Field(default_factory=list)
    exec_part: ExecutionPart = Field(default_factory=ExecutionPart)

class ModuleSubprogramScope(FortranNode):
    kind: Literal["module_subprogram_scope"] = "module_subprogram_scope"
    # NOTE: Use statements
    use_statements: Optional[list[UseStatement]] = Field(default_factory=list)
    spec_part: list[SpecChild] = Field(default_factory=list)
    exec_part: ExecutionPart = Field(default_factory=ExecutionPart)
    internal_subprograms: list[InternalSubprogramScope] = Field(default_factory=list)

class ModuleScope(FortranNode):
    kind: Literal["module_scope"] = "module_scope"
    # NOTE: Use statements
    use_statements: Optional[list[UseStatement]] = Field(default_factory=list)
    spec_part: list[SpecChild] = Field(default_factory=list)
    module_subprograms: list[ModuleSubprogramScope] = Field(default_factory=list)

class MainProgramScope(FortranNode):
    kind: Literal["main_program_scope"] = "main_program_scope"
    # NOTE: Use statements
    use_statements: Optional[list[UseStatement]] = Field(default_factory=list)
    spec_part: list[SpecChild] = Field(default_factory=list)
    exec_part: ExecutionPart = Field(default_factory=ExecutionPart)
    internal_subprograms: list[InternalSubprogramScope] = Field(default_factory=list)

class ExternalSubprogramScope(FortranNode):
    kind: Literal["external_subprogram_scope"] = "external_subprogram_scope"
    # NOTE: Use statements
    use_statements: Optional[list[UseStatement]] = Field(default_factory=list)
    spec_part: list[SpecChild] = Field(default_factory=list)
    exec_part: ExecutionPart = Field(default_factory=ExecutionPart)
    internal_subprograms: list[InternalSubprogramScope] = Field(default_factory=list)
    
class BlockDataScope(FortranNode):
    kind: Literal["block_data_scope"] = "block_data_scope"
    spec_part: list[SpecChild] = Field(default_factory=list)

# 更新 Pydantic 的引用
ExecutionPart.model_rebuild()
InterfaceBodyScope.model_rebuild()
InterfaceBlockScope.model_rebuild()