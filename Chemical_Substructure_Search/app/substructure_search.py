from rdkit import Chem


def substructure_search(molecules, sub_smile):
    valid_molecules = []
    sub_mol = Chem.MolFromSmiles(sub_smile, sanitize=False)
    if not sub_mol or not molecules:
        return []
    for mol in molecules.values():
        mol_smile = mol.smile_notation
        mol_whole = Chem.MolFromSmiles(mol_smile)
        if mol_whole and mol_whole.HasSubstructMatch(sub_mol):
            valid_molecules.append(mol_smile)
    return valid_molecules