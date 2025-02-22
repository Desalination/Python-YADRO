import redis
import json
from Molecule import MoleculeDBJSONEncoder

redis_client = redis.Redis(host='redis', port=6379, db=0)


def get_cached_result(key: str):
    result = redis_client.get(key)
    if result:
        return json.loads(result)
    return None


def set_cache(key: str, value: dict, expiration: int = 200):
    redis_client.setex(key, expiration, json.dumps(value, cls=MoleculeDBJSONEncoder))
