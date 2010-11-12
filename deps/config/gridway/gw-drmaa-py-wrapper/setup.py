from distutils.core import setup

def get_mod_doc_str(mod_file):
    """Return a doc-string from Python module.
    This assumes that doc-string is some string
    in triple quotes at the top of the file."""
    import re
    return re.split(re.compile(r'"""',re.MULTILINE),open(mod_file,"r").read())[1]

setup(
name='GridWay DRMAA wrapper',
version='0.1.0',
author='Andrey Tovchigrechko',
author_email='atovtchi@jcvi.org',
py_modules=["gw_drmaa"],
description='Wrapper module to be imported instead of GridWay DRMAA.py.',
long_description=get_mod_doc_str("gw_drmaa.py"),
)

