"""Importing `log` from this module starts a logging handler, which can be logged to 
with `log.info("log this info")` etc.
"""
from pywapor.general.log_indenter import IndentedLoggerAdapter
import logging
import os
from tqdm.contrib.logging import logging_redirect_tqdm
import warnings
import sys

class MyCoolFormatter(logging.Formatter):
    """
    See https://jupyterbook.org/en/stable/content/code-outputs.html#ansi-outputs
    for colors. E.g. below:
    30 gives black font (and 1;30 would give bold-black font)
    41 gives red background.
    """
    def format(self, record):
        if record.levelno > logging.INFO:
            self._style._fmt = "\x1b[30;41m%(message)s\x1b[0m"
        else:
            self._style._fmt = "%(message)s"
        return super().format(record)

log_settings = logging.getLogger(__name__)
__handler__ = logging.StreamHandler(stream=sys.stdout)
__handler__.setFormatter(MyCoolFormatter())
__handler__.setLevel("INFO")
log_settings.addHandler(__handler__)
log = IndentedLoggerAdapter(log_settings)

def adjust_logger(log_write, folder, log_level, testing = False):
    """Function to adjust the default logging settings that get initiated by
    importing this module.

    Parameters
    ----------
    log_write : bool
        Stop or start writing to log.
    folder : str
        Path to folder in which to store `"log.txt"`.
    log_level : str
        Set the log level.
    """
    if log_write and any([isinstance(x, logging.FileHandler) for x in log_settings.handlers]):
        handlers_to_remove = [x for x in log_settings.handlers if isinstance(x, logging.FileHandler)]
        for handler in handlers_to_remove:
            log_settings.removeHandler(handler)

    if log_write and not any([isinstance(x, logging.FileHandler) for x in log_settings.handlers]):
        if not os.path.isdir(folder):
            os.makedirs(folder)
        handler = logging.FileHandler(filename = os.path.join(folder, "log.txt"), encoding="utf-8")
        formatter = logging.Formatter('{asctime} {levelname:>9s}: {message}', style='{')
        handler.setFormatter(formatter)
        handler.setLevel("DEBUG")
        log_settings.addHandler(handler)
        logging_redirect_tqdm(loggers = log_settings.handlers)
        log_settings.propagate = False

    if testing:
        log_settings.propagate = True
        def fancy_warning(message):
            IndentedLoggerAdapter(log_settings).warning(message)
            warnings.warn(message, stacklevel=2)
        log.warning = fancy_warning

    log_settings.setLevel(log_level)

    if testing:
        from osgeo import gdal
        log.info(f"--> Using GDAL ({gdal.__version__}) from `{gdal.__file__}`.")

if __name__ == "__main__":

    adjust_logger(True, r"/Users/hmcoerver/Local/cog_test", "INFO")

    log.info("tst").add().info("tralala").sub()
    log.warning("balblblblasdas").add()
    log.info("hihihi").add()
    log.sub()
    log.warning("no").add()
    log.sub().sub().info("yes")