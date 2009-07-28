"""Convenience methods to make it easier to run external programs"""

from subprocess import Popen, call, PIPE
import os

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
    dryRun = False
    if 'dryRun' in kw:
        dryRun = kw['dryRun']
        del kw['dryRun']
    if dryRun:
        print popenargs
    else:
        # convert something like run("ls -l") into run("ls -l",shell=True)
        if isinstance(popenargs[0],str) and len(popenargs[0].split()) > 1:
            kw.setdefault("shell",True)
        try:
            if options.debug > 0:
                print popenargs
        except:
            pass
        returncode = call(*popenargs,**kw)
        if returncode != 0:
            raise CalledProcessError(returncode=returncode,cmd=str(popenargs))

def backsticks(*popenargs,**kwargs):
    """Similar to shell backsticks, e.g. a = `ls -1` <=> a = backsticks(['ls','-1']).
    If 'dryRun=True' is given as keyword argument, then 'dryRet' keyword must provide a value
    to return from this function."""
    kw = {}
    kw.update(kwargs)
    dryRun = False
    if 'dryRun' in kw:
        dryRun = kw['dryRun']
        del kw['dryRun']
    dryRet = None
    if 'dryRet' in kw:
        dryRet = kw['dryRet']
        del kw['dryRet']
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
    run(["mkdir","-p",path],dryRun=dryRun)

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
    try:
        os.remove(path)
    except OSError:
        pass

def chmod(path,mode,opts='',dryRun=False):
    if isinstance(path,basestring):
        path = [path]
    else:
        path = list(path)
    run(["chmod"]+opts.split()+[mode]+path,dryRun=dryRun)

