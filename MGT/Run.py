### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Convenience methods to make it easier to run external programs"""

from subprocess import Popen, call, PIPE
import os,glob
import shlex

defineCalledProcessError = False
try:
    from subprocess import CalledProcessError
except ImportError:
    defineCalledProcessError = True

if defineCalledProcessError:

    class CalledProcessError(OSError):
        def __init__(self,returncode,cmd,*l,**kw):
            OSError.__init__(self,*l,**kw)
            self.cmd = cmd
            self.returncode = returncode

def run(*popenargs, **kwargs):
    kw = {}
    kw.update(kwargs)
    dryRun = kw.pop('dryRun',False)
    if dryRun:
        print popenargs
    else:
        # convert something like run("ls -l") into run("ls -l",shell=True)
        if isinstance(popenargs[0],str) and len(shlex.split(popenargs[0])) > 1:
            kw.setdefault("shell",True)
        try:
            _gdebug = options.debug
        except:
            _gdebug = 0
        if kw.pop('debug',_gdebug) > 0:
            print popenargs
        if kw.pop("supressAllOutput",False):
            stdnull = open(os.devnull,"w") #incompat with close_fds on Windows
            kw.setdefault("stdout",stdnull)
            kw.setdefault("stderr",stdnull)
        else:
            stdnull = None
        returncode = call(*popenargs,**kw)
        if stdnull:
            stdnull.close()
        if returncode != 0:
            raise CalledProcessError(returncode=returncode,cmd=str(popenargs))

def backsticks(*popenargs,**kwargs):
    """Similar to shell backsticks, e.g. a = `ls -1` <=> a = backsticks(['ls','-1']).
    If 'dryRun=True' is given as keyword argument, then 'dryRet' keyword must provide a value
    to return from this function."""
    kw = {}
    kw.update(kwargs)
    dryRun = kw.pop('dryRun',False)
    dryRet = kw.pop('dryRet',None)
    if dryRun:
        print popenargs
        return dryRet
    else:
        try:
            if options.debug > 0:
                print popenargs
        except:
            pass
        kw['stdout'] = PIPE
        p = Popen(*popenargs, **kw)
        retout = p.communicate()[0]
        if p.returncode != 0:
            raise CalledProcessError(returncode=p.returncode,cmd=str(popenargs))
        return retout

def makedir(path,dryRun=False):
    """Create one dir with pathname path or do nothing if it already exists.
    Same as Linux 'mkdir -p'."""
    run(["mkdir","-p",path],dryRun=dryRun)

def makedirs(paths,dryRun=False):
    """Create muiltiple dirs with the same semantics as makedir"""
    for path in paths:
        makedir(path=path,dryRun=dryRun)

def makeFilePath(fileName):
    """Assume that the argument is a file name and make all directories that are part of it"""
    dirName = os.path.dirname(fileName)
    if dirName not in ("","."):
        makedir(dirName)

#perhaps use shutil.rmtree instead?    
def rmdir(path,dryRun=False):
    run(["rm","-rf",path],dryRun=dryRun)

rmrf = rmdir

def rmf(path,dryRun=False):
    for f in glob.iglob(path):
        try:
            if os.path.exists(f):
                os.remove(f)
        except OSError:
            pass

def rmfMany(paths,dryRun=False):
    for f in paths:
        try:
            os.remove(f)
        except OSError:
            pass

def remakedir(path,dryRun=False):
    """Create an empty dir with a given path.
    If path already exists, it will be removed first."""
    rmrf(path,dryRun=dryRun)
    makedir(path,dryRun=dryRun)

def chmod(path,mode,opts='',dryRun=False):
    if isinstance(path,basestring):
        path = [path]
    else:
        path = list(path)
    run(["chmod"]+opts.split()+[mode]+path,dryRun=dryRun)

