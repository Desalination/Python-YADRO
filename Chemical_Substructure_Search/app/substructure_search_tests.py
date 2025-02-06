from substructure_search import substructure_search
from datab import engine, Base
from main import app

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient
from celery.result import AsyncResult
import time

client = TestClient(app)


@pytest.fixture()
def test_db():
    Base.metadata.create_all(bind=engine)
    yield
    Base.metadata.drop_all(bind=engine)


# def test_substructure_search_valid():
#     molecules = [
#         "CCO",
#         "CCCO",
#         "CCCCO",
#     ]
#     substructure = "O"

#     result = substructure_search(molecules, substructure)
#     assert result == ["CCO", "CCCO", "CCCCO"], "Substructure search failed"


# def test_substructure_search_invalid_substructure():
#     molecules = [
#         "CCO",
#         "CCCO",
#     ]
#     substructure = "invalid_smiles"

#     with pytest.raises(ValueError, match="Invalid substructure provided."):
#         substructure_search(molecules, substructure)


# def test_substructure_search_empty_molecules():
#     molecules = []
#     substructure = "O"

#     result = substructure_search(molecules, substructure)
#     assert result == [], "Substructure search should return empty list for empty molecules list"


# def test_substructure_search_some_match():
#     molecules = [
#         "CCO",
#         "CCCO",
#         "CCCCO",
#         "CCN",
#     ]
#     substructure = "O"

#     result = substructure_search(molecules, substructure)
#     assert result == ["CCO", "CCCO", "CCCCO"], "Substructure search failed to find matching molecules"


def test_cache(test_db):
    # Base.metadata.drop_all(bind=engine)
    mol = "CCO"
    for i in range(0, 100):
        client.post("/molecules/add", json={"id": str(i), "smile_notation": mol})
        mol = "C" + mol

    response = client.post("/molecules/search", json={"substr": "O"})
    # response = client.get('/molecules/' + str(task_id))
        
    data = response.json()
    print("DATA----------=", data)
    result = AsyncResult(data['task_id'])

    while not result.ready():
        time.sleep(1)

    assert result.status == 'SUCCESS', f"Task failes with status {result.status}"


