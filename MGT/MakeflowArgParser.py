"""Parser for Makeflow arguments to allow modifying them before passing on to Makeflow"""
from MGT.Common import *
import argparse

def add_makeflow_args(parser=None):
    
    if parser is None:
        
        parser = argparse.ArgumentParser(

        description="""Run Makeflow"""

        )

    parser.add_argument("-T","--batch-type", type=str, 
            help="""Makeflow: Batch system type [%(default)s].""",
            default="local")

    parser.add_argument("-B","--batch-options", type=str, 
            help="""Makeflow: Add these options to all batch submit files. 
            Quote in '' if spaces are present inside.""")

    parser.add_argument("-J","--max-remote", type=int, default=500,
            help="""Makeflow: maximum number of concurrent batch jobs [%(default)s]"""
            )
    parser.add_argument("-r","--retry-count", type=int, default=3,
            help="""Automatically retry failed batch jobs up to n times [%(default)s]"""
            )
    parser.add_argument("-S","--submission-timeout", type=int, default=30,
            help="""Time to retry failed batch job submission [%(default)s]"""
            )

    return parser

def parse_makeflow_args(args,parser=None):
    """Parse a subset of Makeflow arguments.
    @param args string with arguments
    @param parser If not None, should be already 
    populated ArgumentParser object. If None,
    new object will be created with @see add_makeflow_args
    @return (options,remaining arguments) as returned by
    parser.parse_known_args()
    """
    if parser is None:
        parser = add_makeflow_args()
    if is_string(args):
        args = shlex.split(args)
    return parser.parse_known_args(args)

def unparse_makeflow_args(opt,other_args=[]):
    """Take the Namespace object created by parse_makeflow_args() and return command line args.
    This should be aware of any possible options parsed into opt.
    There is no generic way to reverse argparse parsing because it is not necesserily
    a one-to-one relationship"""
    args = []
    opt_d = vars(opt) #will it sill work if Namespace is something weird?
    #for now parsing is simple, so we can do it automatically
    for key, val in sorted(opt_d.items()):
        if val is not None:
            args += ["--"+key.replace("_","-"),str(val)]
    return args + other_args

