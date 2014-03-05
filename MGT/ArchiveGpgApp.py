### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


"""Application interface to archiving and unarchiving of directories to/from flat files."""

__all__ = ["ArchiveGpgApp"]

from MGT.App import *
from MGT.Common import *

import shlex

class ArchiveGpgApp(App):
    
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
            
            optParseMakeOption_Path(None, "--path",
            dest="path",
            help="In archiving mode, the path to be archived. "+\
                    "The directory itself is always archived if passed as 'path'; leading components are stripped; "+\
                    "if you want to extract it as a tar-bomb, pass to tar '--strip-components 1' when extracting. "+\
                    "In extraction mode, a single entry is allowed and must be a directory. "+\
                    "The archive content will be extracted into that directory."),
            
            optParseMakeOption_Path(None, "--archive",
            dest="archive",
            help="Path of the archive"),

            make_option(None, "--gpg-args",
            action="store", 
            type="string",
            dest="gpgArgs",
            help="Extra arguments to 'gpg' executable, to pass as a single string (e.g. '-r mgtaxa@jcvi.org'"),
            
            make_option(None, "--tar-args",
            action="store", 
            type="string",
            dest="tarArgs",
            help="Extra arguments to 'tar', to pass a single string (e.g. '--strip-components 1'"),
        ]
        return Struct(usage = "Use gpg-zip to archive and encrypt the content of a given path or extract content of an archive.\n"+\
                "%prog [options]",option_list=option_list)

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        opt = options
        globOpt = globals()["options"]
        assert opt.path
        assert opt.archive
        if not globOpt.isUndef("toolGpgKeyName"): 
            opt.setIfUndef("gpgArgs","-r %s" % (globOpt.toolGpgKeyName,))
    
    def initWork(self,**kw):
        opt = self.opt
        self.gpgZipExe = "gpg-zip"
    
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
        cmd = [self.gpgZipExe,"-e"]
        #gpg-zip is sensitive to options order, --xxx-args must be at the end
        cmd += ["--output",opt.archive]
        gpgArgs = "--sign"
        if opt.gpgArgs:
            gpgArgs += " "+opt.gpgArgs
        cmd += ["--gpg-args",opt.gpgArgs]
        path = opt.path
        assert len(shlex.split(path)) == 1,"I cannot work with output paths that contain whitespaces: %s" % (path,)
        dirPath,basePath = os.path.split(path)
        #tarArgs = "-over-write "
        #if opt.tarArgs:
        #    tarArgs += opt.tarArgs
        #cmd += ["--tar-args",tarArgs]
        if opt.tarArgs:
            cmd += ["--tar-args",opt.tarArgs]
        cmd.append(basePath)
        #gpg-zip will stop and ask if file exists
        if os.path.exists(opt.archive):
            os.remove(opt.archive)
        run(cmd,cwd=dirPath,debug=False,supressAllOutput=True)
    
    def extract(self,**kw):
        opt = self.opt
        outPath = opt.path
        cmd = [self.gpgZipExe,"-d"]
        #gpg-zip is sensitive to options order, --xxx-args must be at the end
        gpgArgs = "--verify"
        if opt.gpgArgs:
            gpgArgs += " "+opt.gpgArgs
        cmd += ["--gpg-args",opt.gpgArgs]
        #--tar-args are not passed through the shell, so no need to sanitize special symbols
        #but it will be split on whitespaces, so file name with whitespaces can potentially
        #inject arbitrary tar commands. Here we just refuse to work with such things.
        assert len(shlex.split(outPath)) == 1,"I cannot work with output paths that contain whitespaces: %s" % (outPath,)
        tarArgs = "-C " + outPath + " "
        if opt.tarArgs:
            tarArgs += opt.tarArgs
        cmd += ["--tar-args",tarArgs]
        cmd.append(opt.archive)
        makedir(outPath)
        run(cmd,debug=False,supressAllOutput=True)


    def iterMembers(self,fileNamePatt=None):
        """Return an iterator over members of an existing archive.
        This is meant to be called in-process when the mode is 'extract'.
        Use case: you need to know what is inside the archive before you 
        batch-submit the extract job.
        Parameters are taken from self.opt.
        @param archive Path of existing archive
        @param fileNamePatt File name pattern to filter archive members with 
        using fnmatch.fnmatch
        @return iterator over objects which have attribute 'name')
        """
        opt = self.opt
        #todo Use gpg-zip --list-archive
        raise NotImplementedError()

if __name__ == "__main__":
    #Allow to call this as script
    runAppAsScript(ArchiveGpgApp)

