from os import getenv
from fastapi import FastAPI

from routers.molecules import router_mol
# from src.routers.molecules import router_mol

app = FastAPI()
app.include_router(router_mol, prefix="/molecules", tags=["Molecules"])

@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}