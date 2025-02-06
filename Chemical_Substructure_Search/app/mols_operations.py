from Molecule import MoleculeOut, MoleculeIn, MoleculeDB
from datab import get_db
from celery_worker import celery
from molecules_iterator import MoleculesIterator
from task import perform_substructure_search

from fastapi import APIRouter, HTTPException, Depends
from rdkit import Chem
from sqlalchemy.orm import Session
from utils import logger
from celery.result import AsyncResult

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


@router_mol.put("/update/{id}", response_model=MoleculeOut)
async def update_mol(id: str, upt_mol: MoleculeIn, db: Session = Depends(get_db)):
    logger.info(f"Updating molecule with id {id}.")
    if id != upt_mol.id:
        raise HTTPException(status_code=406, detail="Provided id doesn't match with provided molecule's id.")
    new_smile = Chem.MolFromSmiles(upt_mol.smile_notation, sanitize=False)
    if not new_smile:
        raise HTTPException(status_code=406, detail="Not valid smile.")
    db_molecule = db.query(MoleculeDB).filter(MoleculeDB.id == id).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Identifier not found.")
    db_molecule.smile_notation = upt_mol.smile_notation
    db.commit()
    db.refresh(db_molecule)
    logger.info(f"Molecule with id {id} updated.")
    return db_molecule


@router_mol.delete("/delete/{id}", response_model=MoleculeOut)
async def delete_mol(id: str, db: Session = Depends(get_db)):
    logger.info(f"Deleting molecule with id {id}.")
    db_molecule = db.query(MoleculeDB).filter(MoleculeDB.id == id).first()
    if not db_molecule:
        raise HTTPException(status_code=404, detail="Identifier not found.")
    db.delete(db_molecule)
    db.commit()
    logger.info(f"Molecule with id {id} deleted.")
    return db_molecule


@router_mol.post("/add", response_model=MoleculeOut)
async def add_mol(mol: MoleculeIn, db: Session = Depends(get_db)):
    logger.info(f"Adding molecule {MoleculeIn}.")
    whole_mol = Chem.MolFromSmiles(mol.smile_notation)
    if not whole_mol:
        raise HTTPException(status_code=406, detail="Not valid smile")
    db_molecule = db.query(MoleculeDB).filter(MoleculeDB.id == mol.id).first()
    if db_molecule:
        raise HTTPException(status_code=409, detail="Such identifier already exists.")
    new_smile = MoleculeDB(id=mol.id, smile_notation=mol.smile_notation)
    db.add(new_smile)
    db.commit()
    db.refresh(new_smile)
    logger.info(f"Molecule {MoleculeIn} added.")
    return new_smile


@router_mol.post("/search")
def search_substructure(substr: str):
    task = perform_substructure_search.apply_async(args=[substr])
    return {"task_id": task.id, "status": task.status}


@router_mol.get("/tasks/{task_id}")
def get_task_result(task_id: str):
    task_result = AsyncResult(task_id, app=celery)
    if task_result.state == 'PENDING':
        return {"task_id": task_id, "status": "Task is still processing"}
    elif task_result.state == 'SUCCESS':
        return {"task_id": task_id, "status": "Task completed", "result": task_result.result}
    else:
        return {"task_id": task_id, "status": task_result.state}
