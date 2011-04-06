### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Application interface to archiving and unarchiving of directories to/from flat files."""

__all__ = ["ArchiveApp"]

from MGT.App import *
from MGT.Common import *

import tarfile 

class ArchiveApp(App):
    
    @classmethod
    def makeOptionParserArgs(klass):
        from optparse import make_option
        option_list = [
            make_option("-m", "--mode",
            action="store", 
            type="choice",
            choices=("archive","extract"),
            dest="mode",
            help="Archive or extract"),
            
            make_option(None, "--path",
            action="append", 
            type="string",
            dest="path",
            help="In archiving mode, the path to be archived, multiple entries are allowed. "+\
                    "If it ends with '/', the content of the directory will be archived, "+\
                    "but not the directory itself (similar to how rsync treats ending '/', "+\
                    "and here it will create the 'tar bomb'). Also similar to rsync, "+\
                    "we do not store in archive the components of dirname(path). "+\
                    "In extraction mode, a single entry is allowed and must be a directory. "+\
                    "The archive content will be extracted into that directory."),
            
            make_option(None, "--archive",
            action="store", 
            type="string",
            dest="archive",
            help="Path of the archive"),

            make_option(None, "--compress",
            action="store", 
            type="choice",
            choices=("gz","bz2"),
            dest="compress",
            help="Optionally use compression when creating the archive: gz or bz2"),
            
            make_option(None, "--safe",
            action="store_true", 
            dest="safe",
            help="Consider archive content untrusted and make safety checks before extracting"),
        ]
        return Struct(usage = "Archive content of a given path or extract content of an archive.\n"+\
                "%prog [options]",option_list=option_list)

    def instanceOptionsPost(self,opt):
        """Set (in place) instance-specific options.
        This is called from __init__() and has access to the execution context (such as current dir)."""
        ## parseCmdLinePost will not modify options that are already defined, so we need to do it here
        if isinstance(opt.path,str):
            opt.path = [ opt.path ]
    
    def doWork(self,**kw):
        opt = self.opt
        if opt.mode == "archive":
            return self.archive(**kw)
        elif opt.mode == "extract":
            return self.extract(**kw)
        else:
            raise ValueError("Unknown opt.mode value: %s" % (opt.mode,))

    def archive(self,**kw):
        opt = self.opt
        flags="w"
        if opt.compress:
            if opt.compress in ("gz","bz2"):
                flags+=":%s" % (opt.compress,)
            else:
                raise ValueError("Invalid opt.compress value: %s" % (opt.compress,))
        tar = tarfile.open(opt.archive, flags)
        for name in opt.path:
            # here we both perform rsync-like treatment of trailing '/' and
            # only store the path components after the dirname(name)
            dirname,basename = os.path.split(name)
            if basename and basename != ".": # /dir/ -> ["dir",""], but /dir/. -> ["dir","."]
                tar.add(name,arcname=basename)
            else:
                for entry in os.listdir(name):
                    tar.add(os.path.join(name,entry),arcname=entry)
        tar.close()
        #TMP:
        #strToFile("EXIT_STATUS=0\n","stdout.wrapper","a")
    
    def extract(self,**kw):
        opt = self.opt
        assert len(opt.path)==1,"Only a single destination path is allowed when extracting"
        tar = tarfile.open(opt.archive, "r") #will auto-detect compression
        #DEBUG:
        #tar.list()
        if opt.safe:
            self.checkSafety(tar)
        tar.extractall(path=opt.path[0])
        tar.close()

    def iterMembers(self,fileNamePatt=None):
        """Return an iterator over members of an existing archive.
        This is meant to be called in-process when the mode is 'extract'.
        Use case: you need to know what is inside the archive before you 
        batch-submit the extract job.
        Parameters are taken from self.opt.
        @param archive Path of existing archive
        @param fileNamePatt File name pattern to filter archive members with 
        using fnmatch.fnmatch
        @return iterator over tarfile.TarInfo - like objects (specifically, always
        having attribute 'name')
        """
        opt = self.opt
        tar = tarfile.open(opt.archive, "r") #will auto-detect compression
        #DEBUG:
        #tar.list()
        for tarinfo in tar:
            if fileNamePatt is None:
                yield tarinfo
            else:
                if fnmatch.fnmatch(tarinfo.name,fileNamePatt):
                    yield tarinfo
        tar.close()
    
    @staticmethod
    def _tarinfoStr(tarinfo):
        return " ; ".join([ "%s : %s" % item for item in sorted(tarinfo.__dict__.items()) \
                if not item[0].startswith('_') and not item[0] == "buf" ])

    def checkSafety(self,tar):
        for tarinfo in tar:
            if os.path.isabs(tarinfo.name) or os.path.isabs(os.path.normpath(tarinfo.name)):
                raise ValueError("Archive failed safety check - absolute file name detected: %s" % \
                        (self._tarinfoStr(tarinfo),))
            elif ".." in tarinfo.name or ".." in os.path.normpath(tarinfo.name):
                raise ValueError("Archive failed safety check - upper directory reference is detected: %s" % \
                        (self._tarinfoStr(tarinfo),))
            elif not (tarinfo.isreg() or tarinfo.isdir()):
                #e.g. if archive was artificially manipulated to contain 
                #first A/B where B is a symlink to ../../something,
                #and then A/B/C, then C might be created as ../../something/C (my guess).
                raise ValueError("Archive failed safety check - non-regular files or dirs can lead to exploits: %s" % \
                        (self._tarinfoStr(tarinfo),))


if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(ArchiveApp)

