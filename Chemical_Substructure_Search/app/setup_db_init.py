from Molecule import MoleculeDB
from datab import SessionLocal


def setup_initial_data(default_molls=None, to_init=True):
    if not default_molls and to_init:
        default_molls = [
            MoleculeDB(id="1", smile_notation="CCO"),
            MoleculeDB(id="2", smile_notation="c1ccccc1"),
            MoleculeDB(id="3", smile_notation="CC(=O)O"),
            MoleculeDB(id="4", smile_notation="CC(=O)Oc1ccccc1C(=O)O"),
        ]
    db = SessionLocal()
    try:
        for mol in default_molls:
            if not db.query(MoleculeDB).filter_by(id=mol.id).first():
                db.add(mol)
        db.commit()
    finally:
        db.close()
