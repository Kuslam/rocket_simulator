import os 
import datetime

def get_project_folder():
    # Currently in project/src/filepaths.py
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def get_logs_folder():
    return os.path.join(get_project_folder(), 'data')

def timestamp():
    return datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

def make_data_folder(name="run"):
    folder_name = f"{name}_{timestamp()}"
    log_folder = os.path.join(get_logs_folder(), folder_name)
    os.makedirs(log_folder)
    return log_folder