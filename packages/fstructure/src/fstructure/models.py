from __future__ import annotations

from typing import Literal, Union, Optional
from pydantic import BaseModel, ConfigDict, Field
from fgraphutils import CodeSnippet

class FortranNode(BaseModel):
    """Base class for all nodes"""
    kind: str
    snippet: CodeSnippet
    model_config = ConfigDict(frozen=True, extra="forbid")

class SpecCodeBlock(FortranNode):
    """
    Specification part code block
    """
    kind: Literal["spec_code_block"] = "spec_code_block"

class ExecCodeBlock(FortranNode):
    """
    Executable code block
    """
    kind: Literal["exec_code_block"] = "exec_code_block"

class RenameRule(BaseModel):
    local_name: str
    remote_name: str
    
    model_config = ConfigDict(frozen=True, extra="forbid")

class UseStatement(FortranNode): # Recommended to inherit from FortranNode to retain position information
    kind: Literal["use_statement"] = "use_statement"
    module_name: str
    # nature: Literal["intrinsic", "non_intrinsic", "unspecified"] = "unspecified" # Optional: F03 feature
    only_list: Optional[tuple[str, ...]] = None  # None means import all; [] means import no symbols (only for dependency)
    renames: tuple[RenameRule, ...] = Field(default_factory=tuple)

# --- Program unit structures ---
class InternalSubprogramScope(FortranNode):
    name: str = None
    kind: Literal["internal_subprogram_scope"] = "internal_subprogram_scope"
    # NOTE: Use statements
    use_statements: Optional[tuple[UseStatement, ...]] = Field(default_factory=tuple)
    spec_part: tuple[SpecCodeBlock, ...] = Field(default_factory=tuple)
    exec_part: ExecCodeBlock
    
    @property
    def key(self) -> str:
        return f"internal_subprogram:{self.name}:{self.snippet.short_hash()}"

class ModuleSubprogramScope(FortranNode):
    name: str = None
    kind: Literal["module_subprogram_scope"] = "module_subprogram_scope"
    # NOTE: Use statements
    use_statements: Optional[tuple[UseStatement, ...]] = Field(default_factory=tuple)
    spec_part: SpecCodeBlock
    exec_part: ExecCodeBlock
    internal_subprograms: tuple[InternalSubprogramScope, ...] = Field(default_factory=list)
    
    @property
    def key(self) -> str:
        return f"module_subprogram:{self.name}:{self.snippet.short_hash()}"

class ModuleScope(FortranNode):
    name: str
    kind: Literal["module_scope"] = "module_scope"
    # NOTE: Use statements
    use_statements: Optional[tuple[UseStatement, ...]] = Field(default_factory=tuple)
    spec_part: SpecCodeBlock
    module_subprograms: tuple[ModuleSubprogramScope, ...] = Field(default_factory=tuple)
    
    @property
    def key(self) -> str:
        return f"module:{self.name}:{self.snippet.short_hash()}"

class MainProgramScope(FortranNode):
    name: Optional[str] = None
    kind: Literal["main_program_scope"] = "main_program_scope"
    # NOTE: Use statements
    use_statements: Optional[tuple[UseStatement, ...]] = Field(default_factory=tuple)
    spec_part: SpecCodeBlock
    exec_part: ExecCodeBlock
    internal_subprograms: tuple[InternalSubprogramScope, ...] = Field(default_factory=tuple)
    
    @property
    def key(self) -> str:
        return f"main_program:{self.name}:{self.snippet.short_hash()}"

class ExternalSubprogramScope(FortranNode):
    name: str = None
    kind: Literal["external_subprogram_scope"] = "external_subprogram_scope"
    # NOTE: Use statements
    use_statements: Optional[tuple[UseStatement, ...]] = Field(default_factory=tuple)
    spec_part: SpecCodeBlock
    exec_part: ExecCodeBlock
    internal_subprograms: tuple[InternalSubprogramScope, ...] = Field(default_factory=tuple)
    
    @property
    def key(self) -> str:
        return f"external_subprogram:{self.name}:{self.snippet.short_hash()}"

class BlockDataScope(FortranNode):
    kind: Literal["block_data_scope"] = "block_data_scope"
    spec_part: tuple[SpecCodeBlock, ...] = Field(default_factory=tuple)
    
class FileRoot(BaseModel):
    kind: Literal["file_root"] = "file_root"
    scopes: tuple[FortranNode, ...]