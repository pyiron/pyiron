# coding: utf-8
# Copyright (c) Max-Planck-Institut für Eisenforschung GmbH - Computational Materials Design (CM) Department
# Distributed under the terms of "New BSD License", see the LICENSE file.

import logging

__author__ = "Joerg Neugebauer"
__copyright__ = "Copyright 2019, Max-Planck-Institut für Eisenforschung GmbH - " \
                "Computational Materials Design (CM) Department"
__version__ = "1.0"
__maintainer__ = "Jan Janssen"
__email__ = "janssen@mpie.de"
__status__ = "production"
__date__ = "Sep 1, 2017"


"""
Set the logging level for pyiron
"""


def set_logging_level(level, channel=None):
    """
    Set level for logger

    Args:
        level (str): 'DEBUG, INFO, WARN'
        channel (int): 0: file_log, 1: stream, None: both
    """
    from pyiron.base.settings.generic import Settings
    s = Settings()
    if channel:
        s.logger.handlers[channel].setLevel(level)
    else:
        s.logger.handlers[0].setLevel(level)
        s.logger.handlers[1].setLevel(level)


def setup_logger():
    """
    Setup logger - logs are written to pyiron.log

    Returns:
        logging.getLogger: Logger
    """
    logger = logging.getLogger('pyiron_log')
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    ch = logging.StreamHandler()
    ch.setLevel(logging.WARN)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    try:
        fh = logging.FileHandler('pyiron.log')
    except PermissionError:
        pass
    else:
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger
