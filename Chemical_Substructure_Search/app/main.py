# from setup_db_init import setup_initial_data
from mols_operations import router_mol
from datab import engine, Base

import os
from fastapi import FastAPI

Base.metadata.create_all(bind=engine)
app = FastAPI()
app.include_router(router_mol, prefix="/molecules", tags=["Molecules"])


@app.get("/")
async def get_server():
    return {"server_id": os.getenv("SERVER_ID", "1")}
