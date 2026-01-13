from __future__ import annotations

from typing import Literal
from pydantic import BaseModel, ConfigDict, Field


class FortranScope(BaseModel):
    """
    Base for all scope/visibility nodes.
    """
    kind: str
    model_config = ConfigDict(frozen=True, extra="forbid")


# -------------------------
# Leaf-ish expression scope
# -------------------------

class ImpliedDoScope(FortranScope):
    """implied-DO index var scope (list scope), can be nested."""
    kind: Literal["implied_do_scope"] = "implied_do_scope"
    children: list[ImpliedDoScope] = Field(default_factory=list)


class StatementFunctionScope(FortranScope):
    """statement function dummy-arg scope (one statement)"""
    kind: Literal["statement_function_scope"] = "statement_function_scope"
    children: list[ImpliedDoScope] = Field(default_factory=list)


# -------------------------
# Construct/statement scopes
# -------------------------

class ForallScope(FortranScope):
    """FORALL statement/construct scope"""
    kind: Literal["forall_scope"] = "forall_scope"
    # FORALL body can contain nested FORALL (and expressions => implied-DO).
    children: list[ForallScope | ImpliedDoScope] = Field(default_factory=list)


class AssociateScope(FortranScope):
    """ASSOCIATE ... END ASSOCIATE"""
    kind: Literal["associate_scope"] = "associate_scope"
    # Associate block is executable statements/constructs => may nest these constructs.
    children: list[
        AssociateScope
        | SelectTypeGuardScope
        | ForallScope
        | ImpliedDoScope
    ] = Field(default_factory=list)


class SelectTypeGuardScope(FortranScope):
    """each TYPE IS / CLASS IS / CLASS DEFAULT block"""
    kind: Literal["select_type_guard_scope"] = "select_type_guard_scope"
    # Each guard block is executable statements/constructs => may nest these constructs.
    children: list[
        AssociateScope
        | SelectTypeGuardScope
        | ForallScope
        | ImpliedDoScope
    ] = Field(default_factory=list)


# -------------------------
# F2003 scoping units
# -------------------------

class DerivedTypeDefinitionScope(FortranScope):
    """TYPE ... END TYPE (derived type definition scope)"""
    kind: Literal["derived_type_definition_scope"] = "derived_type_definition_scope"
    # No executable constructs; but initialization/bounds expressions may contain implied-DO.
    children: list[ImpliedDoScope] = Field(default_factory=list)


class InterfaceBodyScope(FortranScope):
    """interface body: SUBROUTINE/FUNCTION ... END within INTERFACE"""
    kind: Literal["interface_body_scope"] = "interface_body_scope"
    # Specification-only, but may still contain nested derived-type definitions / interface bodies,
    # and implied-DO inside specification expressions.
    children: list[
        DerivedTypeDefinitionScope
        | InterfaceBodyScope
        | ImpliedDoScope
    ] = Field(default_factory=list)


class ModuleSubprogramScope(FortranScope):
    """module procedure definition (after MODULE CONTAINS)"""
    kind: Literal["module_subprogram_scope"] = "module_subprogram_scope"
    children: list[
        DerivedTypeDefinitionScope
        | InterfaceBodyScope
        | StatementFunctionScope
        | InternalSubprogramScope
        | AssociateScope
        | SelectTypeGuardScope
        | ForallScope
        | ImpliedDoScope
    ] = Field(default_factory=list)


class InternalSubprogramScope(FortranScope):
    """internal procedure (after CONTAINS)"""
    kind: Literal["internal_subprogram_scope"] = "internal_subprogram_scope"
    # NOTE: internal subprograms cannot contain other internal subprograms.
    children: list[
        DerivedTypeDefinitionScope
        | InterfaceBodyScope
        | StatementFunctionScope
        | AssociateScope
        | SelectTypeGuardScope
        | ForallScope
        | ImpliedDoScope
    ] = Field(default_factory=list)


class ExternalSubprogramScope(FortranScope):
    """external SUBROUTINE/FUNCTION (top-level)"""
    kind: Literal["external_subprogram_scope"] = "external_subprogram_scope"
    children: list[
        DerivedTypeDefinitionScope
        | InterfaceBodyScope
        | StatementFunctionScope
        | InternalSubprogramScope
        | AssociateScope
        | SelectTypeGuardScope
        | ForallScope
        | ImpliedDoScope
    ] = Field(default_factory=list)


class MainProgramScope(FortranScope):
    """PROGRAM ... END PROGRAM"""
    kind: Literal["main_program_scope"] = "main_program_scope"
    children: list[
        DerivedTypeDefinitionScope
        | InterfaceBodyScope
        | StatementFunctionScope
        | InternalSubprogramScope
        | AssociateScope
        | SelectTypeGuardScope
        | ForallScope
        | ImpliedDoScope
    ] = Field(default_factory=list)


class ModuleScope(FortranScope):
    """MODULE ... END MODULE"""
    kind: Literal["module_scope"] = "module_scope"
    # No executable statements at module scope; no statement functions in module spec part.
    children: list[
        DerivedTypeDefinitionScope
        | InterfaceBodyScope
        | ModuleSubprogramScope
        | ImpliedDoScope
    ] = Field(default_factory=list)


class BlockDataScope(FortranScope):
    """BLOCK DATA ... END BLOCK DATA"""
    kind: Literal["block_data_scope"] = "block_data_scope"
    # No executable statements; interface blocks not allowed; derived-type definitions are allowed.
    children: list[
        DerivedTypeDefinitionScope
        | ImpliedDoScope
    ] = Field(default_factory=list)


# ---- resolve forward refs (important for recursive children) ----
for _m in [
    ImpliedDoScope,
    StatementFunctionScope,
    ForallScope,
    AssociateScope,
    SelectTypeGuardScope,
    DerivedTypeDefinitionScope,
    InterfaceBodyScope,
    ModuleSubprogramScope,
    InternalSubprogramScope,
    ExternalSubprogramScope,
    MainProgramScope,
    ModuleScope,
    BlockDataScope,
]:
    _m.model_rebuild()
