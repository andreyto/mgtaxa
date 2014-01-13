"""Support for shallow parsing of command line interfaces (CLI) with rigid structure.
The use case is CLI of Galaxy tools glue code.
When writing tool XML, we have full control over the order of the generated
command line arguments. Instead of writing conditional expressions in the tool XML,
it is easier to always generate a command line with a fixed number and order of arguments,
using None when the value is missing, and parse this in a tool's Python wrapper.
For readability, it still makes sense to use named options, but their presence and
order is always the same."""

__all__ = ["rigid_cli_dispatch","rigid_cli_dict_to_argv"]

import sys

def rigid_cli_dispatch(cli_patterns,argv=None):
    """Dispatch command line from argv based on CLI patterns in cli_patterns.
    @param cli_patterns Dict that described the CLI. Example:
        cli_patterns = {
                "train" : {
                    "args" : "--inp-train-seq $param --inp-train-seq $param --db-imm-archive $param",
                    "func" : train_func
                    }
                }
    @param argv List of actual command line arguments such as sys.argv (which is the default).
    Example that matches the example above for cli_patterns:
        train --inp-train-seq some_value --inp-train-seq some_value2 --db-imm-archive other_value
    @post In the example above, argv line that starts with 'train' will result in a call
    to train_func(argv_dict,argv), where argv_dict is 
    {'--inp-train-seq': ['some_value','some_value2'], '--db-imm-archive': ['other_value']}
    and argv is ['--inp-train-seq', 'some_value', '--inp-train-seq', 'some_value2', '--db-imm-archive', 'other_value']
    """
    if argv is None:
        argv = sys.argv
    argv = sys.argv[1:]
    mode = argv[0]
    try:
        cli_patt = cli_patterns[mode]
    except KeyError:
        raise ValueError("Unknown mode: {} in command line {}".format(mode,argv))
    argv = argv[1:]
    cli_argv = cli_patt["args"].split()
    argv_dict = dict()
    for cli_w,w in zip(cli_argv,argv):
        if cli_w != "$param":
            assert cli_w == w,"Option name mismatch between CLI pattern and argument string: {} != {}".format(cli_w,w)
            key = cli_w
            argv_dict[key] = []
        else:
            argv_dict[key].append(w)
            
    cli_patt["func"](argv_dict=argv_dict,argv=argv)

def rigid_cli_dict_to_argv(argv_dict,remove_null=True,null="None"):
    """Convert dict with a structure as returned by rigid_cli_dispatch() into argv list.
    Use it to generate argv list for command line after manipulating the dict elements"""
    argv = []
    for key, vals in argv_dict.items():
        vals = [ v for v in vals if v != null ]
        if vals:
            for v in vals:
                argv += [key,v]
    return argv


