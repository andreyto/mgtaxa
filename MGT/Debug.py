import pdb
try:
    # this will drop into ipdb shell just like pdb.set_trace()
    # drops into pdb shell
    from IPython.Debugger import Tracer; debug_here = Tracer()
except ImportError:
    debug_here = pdb.set_trace()

