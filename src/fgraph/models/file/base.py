from pydantic import BaseModel

class BaseFile(BaseModel):
    filename: str
    content: str
    
class Span(BaseModel):
    offset: int
    length: int
    
class Location(BaseModel):
    file: str
    span: Span