from abc import ABC, abstractmethod
from pathlib import Path
from typing import Generic, TypeVar, List, Union

# Define generic type T: represents the "in-memory object after loading"
# It can be a simple Dict, a complex AST Tree, or a custom DataClass
T = TypeVar('T')

class DomainDataLoader(ABC, Generic[T]):
    """
    Abstract base class for data loaders.
    
    Responsibilities: manage standardized input directories, define mappings from files
    to in-memory objects, and provide context for Agents.
    """

    def __init__(self, base_dir: Union[str, Path]):
        self.base_dir = Path(base_dir)
        if not self.base_dir.exists():
            raise FileNotFoundError(f"Directory not found: {self.base_dir}")

    @property
    @abstractmethod
    def file_extension(self) -> str:
        """Define the file extension for test case files managed by this module, e.g., '.f90' or '.json'"""
        pass

    @property
    @abstractmethod
    def description(self) -> str:
        """Describe the purpose of this data in natural language, to be used as part of the System Prompt."""
        pass

    def list_case_ids(self) -> List[str]:
        """
        List all available test case IDs (filenames without extension).
        Agents can use these IDs to request data.
        """
        files = sorted(self.base_dir.glob(f"*{self.file_extension}"))
        return [f.stem for f in files]

    def get_file_path(self, case_id: str) -> Path:
        """Get the physical file path by case ID"""
        return self.base_dir / f"{case_id}{self.file_extension}"

    def load_raw_content(self, case_id: str) -> str:
        """Read the raw text content"""
        path = self.get_file_path(case_id)
        if not path.exists():
            raise FileNotFoundError(f"Case {case_id} not found.")
        return path.read_text(encoding='utf-8')

    @abstractmethod
    def parse(self, content: str, case_id: str) -> T:
        """
        [Core Contract]
        Define how to convert a string to an in-memory object of type T.
        Must be deterministic, pure function logic.
        """
        pass

    def load(self, case_id: str) -> T:
        """[Template Method] Load and parse, the standard interface exposed to clients"""
        content = self.load_raw_content(case_id)
        return self.parse(content, case_id)
