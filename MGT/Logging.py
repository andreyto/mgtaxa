### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Some support for logging"""

import logging

def assert_log(cond,msg,logger):
    """If condition is False, log error and raise AssertionError"""
    if not cond:
        logger.error(msg)
        raise AssertionError(msg)

def logging_config(detail="high",level=logging.DEBUG):
    """Some common logging configuration scenarious"""
    if detail == "high":
        logging.basicConfig(level=level, datefmt="%y-%m-%d %H:%M:%S",
            filemode="a",
            format=("%(asctime)s [%(levelname)5.5s] pid:%(process)-5s\t"
                + "thread:%(threadName)10.10s\t"
                + "%(module)20.20s.:%(funcName)-12.12s:%(lineno)-5s:\t"
                + "%(message)s"))
    else:
        raise ValueError("Unknown value detail={}".format(detail))

