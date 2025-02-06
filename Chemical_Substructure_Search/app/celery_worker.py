from celery import Celery

celery = Celery(
    'task',
    broker=f'redis://{os.getenv("REDIS_HOST")}:6379/0',
    backend=f'redis://{os.getenv("REDIS_HOST")}:6379/0',
    include=['task']
)
celery.conf.update(task_track_started=True)
