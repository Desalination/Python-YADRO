from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    model_config = SettingsConfigDict(
        case_sensitive=True,
        env_file=".env",
        validate_assignment=True,
        env_file_encoding="utf-8",
        extra="ignore",
    )

    # db:
    # image: postgres:15
    # container_name: postgres
    # environment:
    # POSTGRES_DRIVER: postgresql
    # POSTGRES_USER: user
    # POSTGRES_PASSWORD: password
    # POSTGRES_DB: db

    # "postgresql://postgres:postgres@db:5432/postgres"

    DB_HOST: str = "db"
    DB_PORT: int = 5432
    DB_NAME: str = "postgres"
    DB_USERNAME: str = "postgres"
    DB_PASSWORD: str = "postgres"


settings = Settings()
