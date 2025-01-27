from Molecule import MoleculeOut, MoleculeIn, MoleculeDB
from substructure_search import substructure_search
from datab import get_db

from fastapi import APIRouter, HTTPException, Depends
from rdkit import Chem
from sqlalchemy.orm import Session
from molecules_iterator import MoleculesIterator
import logging


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

router_mol = APIRouter()


@router_mol.get("/get/{id}", response_model=MoleculeOut)
async def get_mol(id: str, db: Session = Depends(get_db)):
    logger.info(f"Retrieving molecule with identifier: {id}.")
    db_molecule = db.query(MoleculeDB).filter(MoleculeDB.id == id).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Identifier not found.")
    logger.info(f"Retrieved molecule with identifier: {id}.")
    return db_molecule


@router_mol.get("/get_all", response_model=list[MoleculeOut])
async def get_list_mols(db: Session = Depends(get_db), limit: int = 100):
    logger.info(f"Retrieving list of 100 molecules.")
    mol_iter = MoleculesIterator(db, limit)
    mol_list = [it for it in mol_iter]
    logger.info(f"Retrieved list of {len(mol_list)} molecules.")
    return mol_list


@router_mol.get("/search", response_model=list[MoleculeOut])
async def get_substructure_search(substr: str, db: Session = Depends(get_db)) -> list:
    logger.info(f"Retrieving molecules with substructure: {substr}.")
    all_smiles = [sm.smile for sm in db.query(MoleculeDB).all()]
    logger.info(f"Retrieved molecules with substructure: {substr}.")
    return substructure_search(all_smiles, substr)


@router_mol.put("/update/{id}", response_model=MoleculeOut)
async def update_mol(id: str, upt_mol: MoleculeIn, db: Session = Depends(get_db)):
    logger.info(f"Updating molecule with id {id}.")
    if id != upt_mol.id:
        raise HTTPException(status_code=406, detail="Provided id doesn't match with provided molecule's id.")
    new_smile = Chem.MolFromSmiles(upt_mol.smile_notation, sanitize=False)
    if not new_smile:
        raise HTTPException(status_code=406, detail="Not valid smile.")
    smile = db.query(MoleculeDB).filter(MoleculeDB.id == id).first()
    if not smile:
        raise HTTPException(status_code=404, detail="Identifier not found.")
    smile.smile_notation = upt_mol.smile_notation
    db.commit()
    db.refresh(smile)
    logger.info(f"Molecule with id {id} updated.")
    return smile


@router_mol.delete("/delete/{id}", response_model=MoleculeOut)
async def delete_mol(id: str, db: Session = Depends(get_db)):
    logger.info(f"Deleting molecule with id {id}.")
    smile = db.query(MoleculeDB).filter(MoleculeDB.id == id).first()
    if not smile:
        raise HTTPException(status_code=404, detail="Identifier not found.")
    db.delete(smile)
    db.commit()
    logger.info(f"Molecule with id {id} deleted.")
    return smile


@router_mol.post("/add", response_model=MoleculeOut)
async def add_mol(mol: MoleculeIn, db: Session = Depends(get_db)):
    logger.info(f"Adding molecule {MoleculeIn}.")
    whole_mol = Chem.MolFromSmiles(mol.smile_notation)
    if not whole_mol:
        raise HTTPException(status_code=406, detail="Not valid smile")
    smile = db.query(MoleculeDB).filter(MoleculeDB.id == id).first()
    if smile:
        raise HTTPException(status_code=409, detail="Such identifier already exists.")
    new_smile = MoleculeDB(id=mol.id, smile_notaition=mol.smile_notation)
    db.add(new_smile)
    db.commit()
    db.refresh(new_smile)
    logger.info(f"Molecule {MoleculeIn} added.")
    return new_smile


