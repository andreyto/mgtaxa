### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##


from MGT.Driver import *
import tempfile

def getProgOptions(args):
    from optparse import OptionParser, make_option
    option_list = [
        make_option("-i", "--in-seq",
        action="append", type="string",dest="inSeq",help="Input FASTA sequence file(s). Repeat this option to enter multiple files."),
        make_option("-o", "--out-dir",
        action="store", type="string",dest="outDir",default=".",help="Directory for output files"),
        make_option("-f", "--inp-format",
        action="store", type="choice",choices=("gos","ncbi","ca"),dest="inFormat",default="ca",
        help="Format of FASTA defline, used to extract unique sequence IDs. Known formats are: gos - JCVI GOS reads, ncbi, ca - Celera Accembler scaffolds file"),
        make_option("-P", "--sge-project",
        action="store", type="string",dest="PROJECT_CODE",help="JCVI SGE project code"),
    ]
    parser = OptionParser(usage = "usage: %prog [options]",option_list=option_list)
    (options, args) = parser.parse_args(args=args)

    return options,args


if __name__ == "__main__":
    opt,args = getProgOptions(args=None)
    assert opt.inSeq is not None and opt.PROJECT_CODE is not None, "Need at least --in-seq and --sge-project options"
    opt.inSeq = [ os.path.abspath(f) for f in opt.inSeq ]
    opt.outDir = os.path.abspath(opt.outDir)
    optSI = makeOpt()
    workDir = tempfile.mkdtemp(suffix='.tmp', prefix='pred-'+optSI.rank, dir=options.tmpDir)
    if opt.outDir is None:
        opt.outDir = os.getcwd()
    makedir(workDir)
    #os.chdir(workDir)
    optSI.runMode = "inproc" #"inproc" #"batchDep"
    optSI.cwd = workDir    
    optSI.inSeq = opt.inSeq
    optSI.minSampLen = 5000
    optSI.inFormat = opt.inFormat
    app = SeqImportApp(opt=optSI)
    jobs = app.run()
    optP = makeOptParamScanTest()
    optP.runMode = "inproc" #"inproc" #"inproc" #"batchDep"
    optP.cdd = pjoin(optP.cwd,optP.rank)
    optP.cwd = workDir
    optP.sampStorePred = optSI.cwd
    optP.prOpt.clOpt.thresh = [-2000]
    optP.prOpt.name = "t"
    optP.prOpt.outExportDir = opt.outDir
    modes = ["featPred","predict"]
    for mode in modes:
        optP.mode = mode
        app = ParamScanTestApp(opt=optP)
        jobs = app.run(depend=jobs)


