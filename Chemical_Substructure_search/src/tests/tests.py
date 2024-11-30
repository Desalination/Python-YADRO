# import httpx
#
# r = httpx.get("http://localhost:8000/molecules/get_all")
# print(r.json())


from fastapi.testclient import TestClient
from main import app
from models.Molecule import Molecule
# from src.main import app
# from src.models.Molecule import Molecule
client = TestClient(app)

# mol = Molecule(id="1", smile_notation="CC(=O)Oc1ccccc1C(=O)O")

# mol = {id: "1", smile_notation: "CC(=O)Oc1ccccc1C(=O)O"}

# response = client.get("/molecules/get_all")

response = client.get("/molecules/get/2")
print(response)
print(type(response))

# response = client.put("/molecules/update/1", json=mol)
# print(response.json())
#
# response = client.get("/molecules/search", params={"substr": "klm"})
# print(response.json())




