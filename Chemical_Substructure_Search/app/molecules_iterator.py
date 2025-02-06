from Molecule import MoleculeDB


class MoleculesIterator:
    def __init__(self, db, limit):
        self.result = list(db.query(MoleculeDB).offset(0).limit(limit))
        self.limit = limit
        self.index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.index >= self.limit or self.index >= len(self.result):
            raise StopIteration
        item = self.result[self.index]
        self.index += 1
        return item
