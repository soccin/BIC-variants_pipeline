import hashlib
import os.path
import date_helper

def get_md5(file_name):
    return hashlib.md5(open(file_name, 'rb').read()).hexdigest()

def get_file_timestamp(file_name):
    try:
        return date_helper.timestamp_to_str(os.path.getmtime(file_name))
    except OSError:
        return None
