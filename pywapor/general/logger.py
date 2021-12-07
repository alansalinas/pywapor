from python_log_indenter import IndentedLoggerAdapter
import logging
import os
from tqdm.contrib.logging import logging_redirect_tqdm

log_settings = logging.getLogger(__name__)
__formatter__ = logging.Formatter(fmt='%(message)s', datefmt = '%Y/%m/%d %H:%M:%S')
__handler__ = logging.StreamHandler()
__handler__.setFormatter(__formatter__)
log_settings.addHandler(__handler__)
log = IndentedLoggerAdapter(log_settings)
logging_redirect_tqdm(loggers = log_settings.handlers)

def adjust_logger(log_write, folder, log_level):
    if log_write and not any([isinstance(x, logging.FileHandler) for x in log_settings.handlers]):
        handler = logging.FileHandler(filename = os.path.join(folder, "log.txt"))
        # formatter = logging.Formatter(fmt='%(levelname)s %(asctime)s: %(message)s', datefmt = '%Y/%m/%d %H:%M:%S')
        formatter = logging.Formatter(fmt='%(levelname)s %(asctime)s: %(message)s')
        handler.setFormatter(formatter)
        log_settings.addHandler(handler)
        logging_redirect_tqdm(loggers = log_settings.handlers)
    log_settings.setLevel(log_level)