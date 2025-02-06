from datab import Base

import json
from pydantic import BaseModel
from sqlalchemy import Column, String


class MoleculeDB(Base):
    __tablename__ = "Molecules"
    id = Column(String, primary_key=True, index=True)
    smile_notation = Column(String, index=True)


class MoleculeDBJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, MoleculeDB):
            return {"id": obj.id, "smile_notation": obj.smile_notation}
        return json.JSONEncoder.default(self, obj)


class MoleculeOut(BaseModel):
    id: str
    smile_notation: str

    class Config:
        from_attributes = True


class MoleculeIn(BaseModel):
    id: str
    smile_notation: str
