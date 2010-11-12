"""Wrapper around GridWay DRMAA module to have it co-exist with other (SGE) DRMAA modules.
This module should be imported by user code instead of importing GridWay DRMAA module directly.
It imports all GridWay DRMAA attributes into its namespace.
Additionally, the GridWay DRMAA Python should be built with -W,-r$(GW_LOCATION)/lib
linker arguments (passed through LDFLAGS) as opposed to modifying LD_LIBRARY_PATH in
something like .bashrc at run-time, because both SGE and GridWay try to load shared C DRMAA 
libs with the same name libdrmaa.so.
The original GridWay module is also available as drmaa_base attribute in this module in
case you want to see its documentation with the help() command.
"""
import sys as _my_sys, os as _my_os, copy as _my_copy
_my_old_path=_my_copy.copy(_my_sys.path)
try:
    _my_sys.path.append(_my_os.path.join(_my_os.environ["GW_LOCATION"],"lib/python"))
    #Modyfing os.environ["LD_LIBRARY_PATH"] here would not affect loading of shared libs
    #into the current process, because the loader (itself a shared lib) gets loaded at
    #at start-up and reads the LD_LIBRARY_PATH only when.
    import DRMAA as drmaa_base
    from DRMAA import *
finally:
    _my_sys.path=_my_old_path

