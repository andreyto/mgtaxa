

class App:
    """Interface to application object. 
    The application can be started from shell scripts and called from Python code.
    In both modes, it provides a method to schedule a batch queue execution."""

    def __init__(self,**opt):
        """Constructor.
        @param **opt keyword dict - if empty, will be populated from command line arguments."""
        if len(opt) == 0:
            opt, args = self.parseCmdLine()
        if opt.optFile is not None:
            opt = loadObj(opt.optFile)
        self.opt = opt

    def run(self,**kw):
        """Run the application. 
        Should be redefined in derived classes. 
        Should execute 'return self.runBatch(**kw)' when the constructor option 'batch' is set.
        Should run correctly when called with empty keyword dict (using only constructor options)."""
        opt = self.opt
        if opt.batch:
            return self.runBatch(**kw)

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
        @ret list of BatchJob objects corresponding to the sinks (final vertices) in the DAG of submitted jobs"""
        kw = kw.copy()
        kw.setdefault("scriptName",self.getAppName())
        dryRun = kw.get("dryRun",False)
        opt = self.opt
        opt.batch = False #avoid infinite loop
        opt.optFile = None
        cmd = self.getCmdOptFile(**kw)
        if not dryRun:
            dumpObj(opt,cmd.optFile)
        else:
            print "opt = \n", opt
        return [ runBatch(cmd.cmd,dryRun=dryRun,**kw) ]

    def parseCmdLine(self):
        """Obtain options from command line arguments.
        It is called from constructor when its keyword parameter dict is empty.
        Derived classes must redefine this method and parse at least options 
        implemented here: --opt-file -> optFile (None); --batch -> batch (False)."""
        from optparse import OptionParser, make_option
        option_list = [
            make_option(None, "--opt-file",
            action="store", type="string",dest="optFile",default=None,
            help="Load all program otions from this pickled file"),
            make_option(None, "--batch",
            action="store_true", dest="batch",default=False,
            help="Submit itself to batch queue"),
        ]
        parser = OptionParser(usage = "%prog [options]",option_list=option_list)
        (options, args) = parser.parse_args()
        return options,args

    def getRunpyName(self):
        """Return the full import path for the module that, when executed, will create and run our App class.
        It will be executed as 'python -m runpy module_name args'.
        By default, we place App derived classes into modules with the same name as class name.
        So, we return the class name prepended with full import path (here we assume it is inside MGT namespace.
        Derived class can overide this method."""
        #@todo any way to automatically generate the full import path for the current module?
        return "MGT."+self.__class__.__name__

    def getAppName(self):
        """Return mnemonic name for this application to use for example as a prefix of batch script name"""
        s = self.__class__.__name__
        return s[0].lower()+s[1:]

    def getFactory(self):
        """Return a factory function (class object by default) that can be used to create a new instance of self.
        This is needed when the application needs to batch-submit a new instance of itself."""
        return self.__class__

    def getCmd(self):
        """Return command line that starts this application from shell, w/o options.
        python -m runpy mechanism is used (Python 2.5 is required)."""
        return "python -m runpy %s" % (self.getRunpyName(),)

    def getCmdOptFile(self,cwd=os.getcwd(),**kw):
        """Generate unique file name for a new options pickle file and build full command line with it.
        @param cwd optional directory for the new file (current dir by default)
        @ret Struct(optFile,cmd) where optFile is file name, cmd is command line, 
        such as self.getCmd()+' --opt-file '+optFile."""
        out,optFile = makeTmpFile(suffix=".opt.pkl",prefix=self.getAppName()+'.',dir=cwd)
        out.close()
        return Struct(optFile=optFile,cmd=self.getCmd() + " --opt-file %s" % optFile)


