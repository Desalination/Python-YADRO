import uuid
from pydantic import BaseModel
import uuid

class Molecule(BaseModel):
    id: str = str(uuid.uuid4())
    smile_notation: str
