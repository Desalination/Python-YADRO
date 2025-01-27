from fastapi import APIRouter, HTTPException, UploadFile, File
from rdkit import Chem

from models.Molecule import Molecule
from utils.substructure_search import substructure_search
# from src.models.Molecule import Molecule
# from src.utils.substructure_search import substructure_search

router_mol = APIRouter()

molecules_db = {
    "1": Molecule(id="1", smile_notation="CCO"),
    "2": Molecule(id="2", smile_notation="c1ccccc1"),
    "3": Molecule(id="3", smile_notation="CC(=O)O"),
    "4": Molecule(id="4", smile_notation="CC(=O)Oc1ccccc1C(=O)O")
}

@router_mol.get("/get/{id}")
def get_mol(id: str) -> Molecule:
    molecule = molecules_db.get(id)
    if not molecule:
        raise HTTPException(status_code=404, detail="Identifier not found.")
    return molecule

@router_mol.get("/get_all")
def get_list_mols():
    return list(molecules_db.values())

@router_mol.get("/search")
def get_substructure_search(substr: str) -> list:
    list_mols = list(molecules_db.values())
    return substructure_search(list_mols, substr)

@router_mol.put("/update/{id}")
def update_mol(id: str, upt_mol: Molecule):
    if id != upt_mol.id:
        raise HTTPException(status_code=406, detail="Provided id doesn't match with provided molecule's id.")
    whole_mol = Chem.MolFromSmiles(upt_mol.smile_notation)
    if not whole_mol:
        raise HTTPException(status_code=406, detail="Not valid smile.")
    molecule = molecules_db.get(id)
    if not molecule:
        raise HTTPException(status_code=404, detail="Identifier not found.")
    molecule.smile_notation = upt_mol.smile_notation
    return {"message": "Molecule updated", "molecule": molecule}

@router_mol.delete("/delete/{id}")
def delete_mol(id: str):
    if id in molecules_db:
        del molecules_db[id]
        return {"message": "Molecule deleted"}
    else:
        raise HTTPException(status_code=404, detail="Identifier not found.")

@router_mol.post("/add")
def add_mol(mol: Molecule):
    whole_mol = Chem.MolFromSmiles(mol.smile_notation)
    if not whole_mol:
        raise HTTPException(status_code=406, detail="Not valid smile")
    if mol.id in molecules_db:
        raise HTTPException(status_code=409, detail="Such identifier already exists.")
    molecules_db[mol.id] = mol
    return {"message": "New molecule added"}

# @router_mol.post("/upload")
# def upload_mols(file: UploadFile = File(...)):
#     try:
#         with open(file.filename, mode='rw', encoding="utf-8") as lines:
#             count_uploaded_mols = 0
#             for line in lines:
#                 try:
#                     smiles_id, smiles = line.split(',')
#                     smiles_id = smiles_id.strip()
#                     smiles = smiles.strip()
#                     mol = Molecule(id=smiles_id, smile_notation=smiles)
#                     if smiles_id not in molecules_db:
#                         molecules_db[mol.id] = mol
#                         count_uploaded_mols += 1
#                 except ValueError:
#                     continue
#                 return {"message": f"{count_uploaded_mols} molecules uploaded."}
#     except IOError:
#         return {"message": f"Couldn't open file."}
