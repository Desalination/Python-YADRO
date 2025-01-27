from db_config import settings

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import Session, sessionmaker

# DB_HOST: str = "localhost"
# DB_PORT: int = 5432
# DB_NAME: str = "db"
# DB_USERNAME: str = "user"
# DB_PASSWORD: str = "password"

DATABASE_URL = "postgresql://{username}:{password}@{host}:{port}/{db}".format(
    username=settings.DB_USERNAME,
    password=settings.DB_PASSWORD,
    host=settings.DB_HOST,
    port=settings.DB_PORT,
    db=settings.DB_NAME,
)

# DATABASE_URL = "postgresql://user:password@localhost/db"

engine = create_engine(DATABASE_URL, echo=True)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()

def get_db() -> Session:
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
