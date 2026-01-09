from pydantic import BaseModel

class BaseFile(BaseModel):
    path: str
    content: str
    hash: str