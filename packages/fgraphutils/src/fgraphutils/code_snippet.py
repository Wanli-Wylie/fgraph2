from pydantic import BaseModel
import hashlib

class CodeSnippet(BaseModel):
    code: str

    def short_hash(self) -> str:
        # Truncated SHA256 hash to 16 characters
        return hashlib.sha256(self.code.encode()).hexdigest()[:16]