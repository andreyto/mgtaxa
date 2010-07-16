### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Some support for logging and debugging."""
__all__ = ["pdb","set_trace","Memuse","debugMsg"]

# This try block will ignore an error when importing
# this module from PyMol, which has a 'cmd' module that
# supecedes standard 'cmd' module and causes 'pdb'
# import to fail. Pdb will not be available in that case.
try:
    #This tries to use IPython shell support if available
    try:
        from IPython.Debugger import Pdb
        def set_trace():
            from IPython.Debugger import Pdb
            Pdb().set_trace()
        del Pdb
    except:
        import pdb
        def set_trace():
            import pdb
            pdb.set_trace()
        del pdb

    import pdb #standard pdb.set_trace() use
except AttributeError:
    pdb = None

def debugMsg(msg,level=1):
    """
    Conditionally print a message.
    Prints msg if the current debug level in the 
    global options instance is at least at level 
    supplied to this function."""
    if level <= options.debug:
        print msg

from Bits import Memuse

