# from substructure_search import substructure_search
from datab import engine, Base, SessionLocal
from main import app

import pytest
from fastapi.testclient import TestClient
from celery.result import AsyncResult
import time
from utils import logger
# from cache_funcs import get_cached_result

client = TestClient(app)


@pytest.fixture()
def test_db():
    print("START TEST DB")
    db = SessionLocal()
    Base.metadata.drop_all(bind=engine)
    Base.metadata.create_all(bind=engine)
    yield
    Base.metadata.drop_all(bind=engine)
    db.close()
    print("END TEST DB")


def search_for_test(substr: str):
    s = time.perf_counter()
    substr = "O"
    response = client.post("/molecules/search", json={"substr": substr})

    task_id = response.json()['task_id']

    result = AsyncResult(task_id)

    while not result.ready():
        time.sleep(0.1)
    search_time_before = time.perf_counter() - s

    return result, search_time_before


def add_mols_for_test(mol, add_mol, id_start=0, count_of_mols=200):
    for i in range(id_start, count_of_mols + id_start):
        client.post("/molecules/add", json={"id": str(i), "smile_notation": mol})
        mol = add_mol + mol


def test_cache():
    db = SessionLocal()
    Base.metadata.drop_all(bind=engine)
    Base.metadata.create_all(bind=engine)

    mol = "CO"
    add_mol = "C"
    count_of_mols = 200
    add_mols_for_test(mol, add_mol, id_start=0, count_of_mols=200)
    substr = "O"
    result, search_time_before = search_for_test(substr)
    logger.info(f"Search time of NOT cached mols = {search_time_before}")

    assert result.status == 'SUCCESS', f"Task failed with status {result.status}"

    result_repeat, search_time_after = search_for_test(substr)
    logger.info(f"Search time of {count_of_mols} cached mols = {search_time_after}")
    logger.info(f"Search time decreased on: {search_time_before - search_time_after}")

    assert result_repeat.status == 'SUCCESS', f"Task failed with status {result_repeat.status}"

    db.close()

# def test_substructure_search_valid():
#     molecules = [
#         "CCO",
#         "CCCO",
#         "CCCCO",
#     ]
#     substructure = "O"
#
#     result = substructure_search(molecules, substructure)
#     assert result == ["CCO", "CCCO", "CCCCO"], "Substructure search failed"
#
#
# def test_substructure_search_invalid_substructure():
#     molecules = [
#         "CCO",
#         "CCCO",
#     ]
#     substructure = "invalid_smiles"
#
#     with pytest.raises(ValueError, match="Invalid substructure provided."):
#         substructure_search(molecules, substructure)
#
#
# def test_substructure_search_empty_molecules():
#     molecules = []
#     substructure = "O"
#
#     result = substructure_search(molecules, substructure)
#     assert result == [], "Substructure search should return empty list for empty molecules list"
#
#
# def test_substructure_search_some_match():
#     molecules = [
#         "CCO",
#         "CCCO",
#         "CCCCO",
#         "CCN",
#     ]
#     substructure = "O"
#
#     result = substructure_search(molecules, substructure)
#     assert result == ["CCO", "CCCO", "CCCCO"], "Substructure search failed to find matching molecules"
