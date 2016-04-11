import datetime

def now():
    """Date right now"""
    return datetime.datetime.now()

def pretty_today():
    """Current date in %Y-%m-%d format"""
    return datetime.datetime.now().strftime("%Y-%m-%d")

def current_year():
    """Current year"""
    return datetime.datetime.now().year

def pretty_now():
    """Current datetime in %Y-%m-%d %H:%M:%S format"""
    return str(now()).split(".")[0]

def timestamp_now():
    """Current datetime in %Y%m%d%H%M%S format"""
    return datetime.datetime.now().strftime("%Y%m%d%H%M%S") 

def datestamp_now():
    """Current date in %Y%m%d format"""
    return datetime.datetime.now().strftime("%Y%m%d") 

def str_to_date(date_str):
    """String should be in %Y-%m-%d %H:%M:%S format."""
    return datetime.datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S")

def fancy_str_to_date(date_str):
    """String should be in %b %d %Y format."""
    return datetime.datetime.strptime(date_str, "%b %d %Y")

def fancy_str_to_datetime(datetime_str):
    """String should be in %b %d %Y %H:%M format."""
    return datetime.datetime.strptime(datetime_str, "%b %d %Y %H:%M")

def date_to_str(date):
    """String will be in %Y-%m-%d %H:%M:%S format."""
    return date.strftime("%Y-%m-%d %H:%M:%S")

def timestamp_to_str(timestamp):
    return date_to_str(datetime.datetime.fromtimestamp(timestamp))

