#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Set up logger
"""
from pathlib import Path
import sys 
import logging
from logging.handlers import RotatingFileHandler

module_name = 'tiskitpy'

class ShutdownHandler(logging.StreamHandler):
    def emit(self, record):
        super(ExitHandler, self).emit(record)
        if record.levelno >= logging.CRITICAL:
            sys.exit(1)
            
# Logging formatter supporting colorized output
class LogFormatter(logging.Formatter):

    COLOR_CODES = {
        logging.CRITICAL: "\033[1;35m", # bright/bold magenta
        logging.ERROR:    "\033[1;31m", # bright/bold red
        logging.WARNING:  "\033[1;33m", # bright/bold yellow
        logging.INFO:     "\033[0;30m", # black / dark gray
        logging.DEBUG:    "\033[1;30m"  # bright/bold black / dark gray
    }

    RESET_CODE = "\033[0m"

    def __init__(self, color, *args, **kwargs):
        super(LogFormatter, self).__init__(*args, **kwargs)
        self.color = color

    def format(self, record, *args, **kwargs):
        if (self.color == True and record.levelno in self.COLOR_CODES):
            record.color_on  = self.COLOR_CODES[record.levelno]
            record.color_off = self.RESET_CODE
        else:
            record.color_on  = ""
            record.color_off = ""
        return super(LogFormatter, self).format(record, *args, **kwargs)

def init_logger(file_level='DEBUG', console_level='INFO',
                 console_log_output="stdout"):
    """
    Create or open a rotating logging file and add it to ObsinfoConfiguration

    All arguments are only processed the first time this is called
    
    Args:
        file_level (str): level to start printing to file
        console_level (str): level to start printing to screen ()
        console_log_output (str): "stdout" or "stderr"
        
    valid levels are DEBUG, INFO, WARNING, ERROR, CRITICAL
    
    Returns: object of Logger class
    """
    logger = logging.getLogger(module_name)
    if logger.hasHandlers():  # don't recreate handlers
        return logger

    assert file_level.upper() in ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    assert console_level.upper() in ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
    assert console_log_output in ["stdout", "stderr"]

    logfile = Path.home().joinpath(f'.{module_name}', f'{module_name}log')
    if not logfile.parent.exists():
        logfile.parent.mkdir(parents=True, exist_ok=True)

    logger.setLevel(logging.DEBUG)
                             
    # File Handler
    fileHandler = RotatingFileHandler(logfile, maxBytes=200*1024, backupCount=3)
    fileHandler.setFormatter(logging.Formatter(
        fmt='{asctime} {levelname:<8s} {module}.{funcName}(): {message}',
        style='{', datefmt='%Y-%m-%d %H:%M:%S'))
    fileHandler.setLevel(file_level.upper())

    # Console Handler
    console_log_output = console_log_output.lower()
    if (console_log_output == "stdout"):
        console_log_output = sys.stdout
    elif (console_log_output == "stderr"):
        console_log_output = sys.stderr
    else:
        raise ValueError(f"Invalid console output: {console_log_output}")
    consoleHandler = logging.StreamHandler(console_log_output)
    consoleHandler.setFormatter(LogFormatter(
        # fmt="{color_on}[{threadName}] [{levelname:<8s}] {message}{color_off}",
        fmt="{color_on}[{levelname}] {message}{color_off}",
        # fmt="{color_on}[{module}][{funcName}()] [{levelname:<8s}] {message}{color_off}",
        color=True, style='{'))
    consoleHandler.setLevel(console_level.upper())
    
    logger.addHandler(fileHandler)
    logger.addHandler(consoleHandler)

    return logger
    
def change_level(logger, htype, level):
    """
    Change a handler's logging level
    
    Args:
        htype (str): handler type: 'console' or 'file'
        level (str): 'DEBUG', 'INFO', 'WARNING', 'ERROR', or 'CRITICAL'
    """
    assert htype.upper() in ['CONSOLE', 'FILE']
    assert level.upper() in ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']

    if htype.upper() == 'CONSOLE':
        handler_type = logging.StreamHandler
    else:
        handler_type = logging.RotatingFileHandler
    for handler in logger.handlers:
        if isinstance(handler, handler_type):
            handler.setLevel(level)

def change_console_level(logger, level):
    """
    Change consoleHandlers logging level
    
    Args:
        level (str): 'DEBUG', 'INFO', 'WARNING', 'ERROR', or 'CRITICAL'
    """
    change_level(logger, 'console', level)

def change_file_level(logger, level):
    """
    Change fileHandler's logging level
    
    Args:
        level (str): 'DEBUG', 'INFO', 'WARNING', 'ERROR', or 'CRITICAL'
    """
    change_level(logger, 'file', level)
