from celery_worker import celery
from Molecule import MoleculeDB
from cache_funcs import get_cached_result, set_cache
from datab import SessionLocal
from utils import logger

from rdkit import Chem


@celery.task
def perform_substructure_search(substr: str):
    db = SessionLocal()
    sub_mol = Chem.MolFromSmiles(substr, sanitize=False)
    if not sub_mol:
        return {"error": "Invalid substructure SMILES string."}

    matching_molecules = get_cached_result(substr)

    if matching_molecules:
        logger.info(f"Substructure {substr} found in cache.")
        result = {
            "total_found": len(matching_molecules),
            "matching_molecules": [
                {"identifier": mol['id'], "smiles": mol['smile_notation']} for mol in matching_molecules
            ]
        }
    else:
        matching_molecules = db.query(MoleculeDB).filter(MoleculeDB.smile_notation.like(f"%{substr}%")).all()
        set_cache(substr, matching_molecules)
        result = {
            "total_found": len(matching_molecules),
            "matching_molecules": [
                {"identifier": mol.id, "smiles": mol.smile_notation} for mol in matching_molecules
            ]
        }

    return result
