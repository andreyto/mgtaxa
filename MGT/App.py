### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Common import *

__all__ = ["App","runAppAsScript"]

class App:
    """Interface to application object. 
    The application can be started from shell scripts or called in-process from Python code.
    In both modes, it provides a method to schedule a batch queue execution, possibly
    represented as a DAG of dependent jobs.
    If you need to run it as a script from shell, call runAppAsScript from module level,
    passing the proper derived class (see code example at the bottom of this source file).
    The uniform way of providing scripting interface is to place each App derived class in
    its own source file with the same name as the class, and include the mentioned
    code snippet at the bottom."""

    ## Derived classes should set this to a list of opt.mode values that can result in submision of new batch jobs.
    batchDepModes = tuple()

    def __init__(self,args=[],opt=Struct()):
        """Constructor.
        @param args optional command line arguments to parse -
        if executing as a module, should pass None or sys.argv[1:], else - [], which should result in 
        default values generated for all options not defined by opt.
        @param opt Struct instance - values defined here will override those parsed from args.
        Two basic use patterns: 
        if running as a shell script and parsing actual command line arguments:
            app = App(args=None)
            app.run()
        if calling from Python code:
            construct opt as a Struct() instance, specify only attributes with non-default values, then
            app = App(opt=opt) #that fills yet unset options with defaults
            app.run() #either runs in the same process, or submits itself or a set of other Apps to batch queue
        """
        optArgs, args = self.parseCmdLine(args=args,_explicitOpt=opt)
        optArgs.updateOtherMissing(opt)
        self._instanceOptionsPostBase(opt)
        self.instanceOptionsPost(opt)
        # this is now checked in parseCmdLine
        #if opt.optFile is not None:
        #    opt = loadObj(opt.optFile)
        self.opt = opt

    def run(self,**kw):
        """Run the application. 
        This is the method that is actually called by the user code. 
        Dispatches the work execution or batch submission depending on the options.
        @return list of sink BatchJob objects (if executed syncroniously, an empty list)"""
        opt = self.opt
        runMode = self._ajustRunMode(**kw)
        if runMode == "batch":
            return self.runBatch(**kw)
        elif runMode in ("inproc","batchDep"):
            curdir = os.getcwd()
            try:
                if "cwd" in opt:
                    os.chdir(opt["cwd"])
                #DBG:
                #time.sleep(5)
                self.initWork(**kw)
                ret = self.doWork(**kw)
            finally:
                if "cwd" in opt:
                    os.chdir(curdir)
            if ret is None:
                ret = []
            return ret
        else:
            raise ValueError(runMode)

    def _ajustRunMode(self,**kw):
        opt = self.opt
        runMode = opt.runMode
        runMode = kw.get("runMode",runMode) #keyword overrides just for this instance
        if options.app.runMode not in ("default","batch","batchDep"):
            #we get infinite loop with any "batch*"
            runMode = options.app.runMode
        depend = self._ajustDepend(**kw)
        if depend is not None and len(depend) > 0:
            #we can only run as "batch" or "batchDep" if we need to wait for dependency jobs
            #@todo make the current process to poll qstat at this location
            #if runMode is not "inproc".
            assert runMode in ("batch","batchDep")
        if runMode == "batchDep":
            if opt.mode not in self.batchDepModes:
                runMode = "batch"
        if hasattr(opt,"batchResume") and opt.batchResume:
            runMode = "inproc"
            opt.batchResume = False
        return runMode

    def _ajustDepend(self,**kw):
        """squash all None is 'depend' list. "None" can come from other run() calls when runMode=="inproc"."""
        depend = [ d for d in kw.get("depend",[None]) if d is not None ]
        return depend
        
    def initWork(self,**kw):
        """Perform common initialization right before doing the actual work in doWork().
        Must be redefined in the derived classes.
        Should not be called directly by the user except from initWork() in a derived class.
        This one can create large objects because they are not passed through the batch submission,
        but immediately used within the same process."""
        pass

    def doWork(self,**kw):
        """Do the actual work.
        Must be redefined in the derived classes.
        Should not be called directly by the user except from doWork() in a derived class.
        Should work with empty keyword dict, using only self.opt.
        If doing batch submision of other App instances, must return a list of sink (final) BatchJob objects."""
        pass

    def runBatch(self,**kw):
        """Submit this application to run in batch mode.
        Can be redefined in derived classes.
        self.getCmdOptFile() will provide support for constructing the command line for batch script.
        BatchRun.runBatch() should be used to submit the batch job.
        We provide a generic implementation here that just submits itself with current options.
        Derived classes can define more complex schemes, e.g. fan out N parallel jobs plus
        one that starts on their completion to collect the results. 
        Although batch jobs can re-use runBatch() to submit other jobs, the return value of the top-most
        runBatch() will not know about these sub-jobs and therefore could not define dependencies on them.
        The requirement to this method is that it must have a quick running top-level mode suitable
        to call on the login host.
        It can be called from both self.run() when 'batch' option is set, and directly from external code.
        @return list of BatchJob objects corresponding to the sinks (final vertices) in the DAG of submitted jobs.
        If kw["runMode"] or global options.app.runMode == "inproc" it will call self.run() instead of submitting a batch job."""
        opt = copy(self.opt)
        if opt.runMode == "batch":
            opt.runMode = "inproc" #avoid infinite loop
        opt.batchResume = True
        opt.optFile = None
        kw = kw.copy()
        kw.setdefault("scriptName",self.getAppName())
        dryRun = kw.pop("dryRun",False)
        cmd = self.getCmdOptFile(**kw)
        if not dryRun:
            dumpObj(opt,cmd.optFile)
        else:
            print "opt = \n", opt
        bkw = BatchSubmitter.defaultOptions()
        ## pull only relevant options
        bkw.updateFromOtherExisting(opt)
        ## kw options override all
        bkw.update(kw)
        bkw["depend"] = self._ajustDepend(**bkw.asDict())
        return [ runBatch(cmd.cmd,dryRun=dryRun,**bkw.asDict()) ]

    @classmethod
    def parseCmdLine(klass,args=None,_explicitOpt=None):
        """Obtain options from command line arguments.
        It is called from constructor or directly by the user.
        Derived classes must redefine makeOptionParserArgs() and optionally parseCmdLinePost() 
        if they need to parse additional arguments.
        The arguments defined here are needed by the App infrastructure: 
        --opt-file -> optFile (None); --batch -> batch (False).
        When called directly from the user code, the args is typically set to [] in order to obtain
        the default values for all options. Thus, it is important the the implementation provides
        reasonable defaults in a context independent way (e.g. without including the current directory info).
        @param args command line arguments (pass [] to get default values, pass None to use sys.argv[1:])"""
        from optparse import OptionParser, make_option
        option_list = [
            make_option(None, "--opt-file",
            action="store", type="string",dest="optFile",default=None,
            help="Load all program otions from this pickled file"),
            make_option(None, "--run-mode",
            action="store", type="choice",choices=("batch","inproc","batchDep"),dest="runMode",default="inproc",
            help="Choose to batch-run or in-process"),
        ]
        parseArgs = klass.makeOptionParserArgs()
        parseArgs.option_list.extend(option_list)
        parser = OptionParser(**parseArgs.asDict())
        (options, args) = parser.parse_args(args=args) #values=opt does not work
        options = Struct(options.__dict__)
        if _explicitOpt is not None:
            options.update(_explicitOpt)
        if options.optFile is not None:
            options = loadObj(options.optFile)
        else:
            klass.parseCmdLinePost(options=options,args=args,parser=parser)
        return options,args

    @classmethod
    def makeOptionParserArgs(klass):
        """Return a Struct with optparse.OptionParser constructor arguments specific to the application.
        The "option_list" attribute must be obtained with a sequence of calls to optparse.make_option.
        Other possible attributes can be e.g. "usage".
        This method will be called by parseCmdLine() and the returned "option_list" concatenated with the default
        one provided by the parseCmdLine().
        Must be redefined in the derived class only if there are any application specific command-line options."""
        return Struct(usage="%prog [options]",option_list=[])

    @classmethod
    def parseCmdLinePost(klass,options,args,parser):
        """Optionally modify options and args in-place.
        Called at the end of parseCmdLine to allow the derived classes customizing the option processing.
        @param options options returned by OptionParser and converted to Struct object
        @param args args returned by OptionParser
        @param parser OptionParser object used to parse the command line - needed here to call its error() method
        if necessary.
        options should be modified in place by this method"""
        pass

    @classmethod
    def defaultOptions(klass,_explicitOpt=None):
        return klass.parseCmdLine(args=[],_explicitOpt=_explicitOpt)

    @classmethod
    def fillWithDefaultOptions(klass,options):
        """Fill with default values those options that are not already set.
        @param options Existing options object - will be modified in place and also returned
        @return Modified options parameter"""
        #we use _explicitOpt to override any default options with those from 'options'
        #BEFORE parseCmdLinePost() is called
        optArgs,args = klass.defaultOptions(_explicitOpt=options)
        optArgs.updateOtherMissing(options)
        return options

    def _instanceOptionsPostBase(self,opt):
        """Set (in place) instance-specific options.
        This is called from __init__() and has access to the execution context (such as current dir).
        This method should not be redefined in derived classes. Redefine instanceOptionsPost, which
        is called right after this one."""
        opt.setdefault("cwd",os.getcwd())
        #set current data directory to cwd if not set already
        opt.setdefault("cdd",opt.cwd)
        
    
    def instanceOptionsPost(self,opt):
        """Set (in place) instance-specific options.
        This is called from __init__() and has access to the execution context (such as current dir)."""
        pass

    
    def getAppName(self):
        """Return mnemonic name for this application to use for example as a prefix of batch script name"""
        s = self.__class__.__name__
        return s[0].lower()+s[1:]

    def factory(self,**kw):
        """A factory function (class constructor by default) that creates a new instance of self.
        This is needed when the application needs to batch-submit a new instance of itself."""
        return self.__class__(**kw)

    def getCmd(self):
        """Return command line that starts this application from shell, w/o options.
        "python -c '...'" form is used if we can get the module name of self,
        and "python sys.argv[0]" otherwise."""
        import inspect
        modname = inspect.getmodule(self).__name__
        if modname == "__main__":
            #executed as script, module name is not available but we can execute the same way again
            return sys.executable + " " + sys.argv[0]
        else:
            #modname must have a full import path and we can use python -c '...'
            klassname = self.__class__.__name__
            return sys.executable + " -c 'import %s; %s.%s(args=None).run()'" % (modname,modname,klassname)

    def getCmdOptFile(self,cwd=os.getcwd(),**kw):
        """Generate unique file name for a new options pickle file and build full command line with it.
        @param cwd optional directory for the new file (current dir by default)
        @ret Struct(optFile,cmd) where optFile is file name, cmd is command line, 
        such as self.getCmd()+' --opt-file '+optFile."""
        out,optFile = makeTmpFile(suffix=".opt.pkl",prefix=self.getAppName()+'.',dir=cwd,withTime=True)
        out.close()
        return Struct(optFile=optFile,cmd=self.getCmd() + " --opt-file %s" % optFile)

    def getOpt(self):
        return self.opt

def runAppAsScript(klass):
    """Call this function from module level if running a module with App derived class as a script"""
    app = klass(args=None)
    return app.run()

## This is an example of the code that provides a direct scripting interface to an App derived class.
## It must be present in module with a class derived from App.
## App class name below must be replaced with a derived class name.
#
#if __name__ == "__main__":
#    runAppAsScript(App)
#
