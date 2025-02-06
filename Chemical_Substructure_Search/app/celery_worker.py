from celery import Celery

celery = Celery(
    'task',
    broker='redis://redis:6379/0',
    backend='redis://redis:6379/0',
    include=['task']
)
celery.conf.update(task_track_started=True)