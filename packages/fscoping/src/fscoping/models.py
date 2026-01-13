from __future__ import annotations
from pydantic import BaseModel, ConfigDict
from typing import Literal

class FortranScope(BaseModel):
    """
    Abstract base class for all scoping / visibility units.
    """
    kind: Literal["scope"]

    model_config = ConfigDict(frozen=True)

# ---- F2003 standard scoping units ----

class ModuleScope(FortranScope):
    """MODULE ... END MODULE"""
    kind: Literal["module_scope"] = "module_scope"


class MainProgramScope(FortranScope):
    """PROGRAM ... END PROGRAM"""
    kind: Literal["main_program_scope"] = "main_program_scope"


class ExternalSubprogramScope(FortranScope):
    """External SUBROUTINE / FUNCTION (top-level)"""
    kind: Literal["external_subprogram_scope"] = "external_subprogram_scope"


class BlockDataScope(FortranScope):
    """BLOCK DATA ... END BLOCK DATA"""
    kind: Literal["block_data_scope"] = "block_data_scope"


class ModuleSubprogramScope(FortranScope):
    """Module procedure defined after MODULE CONTAINS"""
    kind: Literal["module_subprogram_scope"] = "module_subprogram_scope"


class InternalSubprogramScope(FortranScope):
    """Internal subprogram defined after CONTAINS"""
    kind: Literal["internal_subprogram_scope"] = "internal_subprogram_scope"


class InterfaceBodyScope(FortranScope):
    """SUBROUTINE / FUNCTION body inside INTERFACE"""
    kind: Literal["interface_body_scope"] = "interface_body_scope"


class DerivedTypeDefinitionScope(FortranScope):
    """TYPE ... END TYPE (derived type definition scope)"""
    kind: Literal["derived_type_definition_scope"] = "derived_type_definition_scope"
