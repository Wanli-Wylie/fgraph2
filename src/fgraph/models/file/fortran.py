from .base import BaseFile
from typing import Literal

class FortranFile(BaseFile):
    kind: Literal["f90", "f95", "f03", "f08", "f18", "f20", "f23"]
    