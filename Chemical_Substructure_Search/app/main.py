from setup_db_init import setup_initial_data
from mols_operations import router_mol, logger
from datab import engine, Base

import os
from fastapi import FastAPI
from fastapi.concurrency import asynccontextmanager
from fastapi.testclient import TestClient



@asynccontextmanager
async def lifespan(app: FastAPI):
    Base.metadata.create_all(bind=engine)
    setup_initial_data()
    yield


app = FastAPI(lifespan=lifespan)

app.include_router(router_mol, prefix="/molecules", tags=["Molecules"])

client = TestClient(app)
response_get2 = client.get("/molecules/get/2")
response_getAll = client.get("/molecules/get_all")
print(response_get2.content)
print(response_getAll.content)

@app.get("/")
async def get_server():
    return {"server_id": os.getenv("SERVER_ID", "1")}
