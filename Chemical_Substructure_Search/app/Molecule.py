from datab import Base

from pydantic import BaseModel
from sqlalchemy import Column, String


class MoleculeDB(Base):
    __tablename__ = "Molecules"
    id = Column(String, primary_key=True, index=True)
    smile_notation = Column(String, index=True)


class MoleculeOut(BaseModel):
    id: str
    smile_notation: str

    class Config:
        from_attributes = True


class MoleculeIn(BaseModel):
    id: str
    smile_notation: str
