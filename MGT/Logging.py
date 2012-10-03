### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Some support for logging"""

def assert_log(cond,msg,logger):
    """If condition is False, log error and raise AssertionError"""
    if not cond:
        logger.error(msg)
        raise AssertionError(msg)

