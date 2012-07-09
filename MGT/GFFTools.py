"""Graphical output for GFF3 files"""

## This module calls external programs from GT (Genome Tools) package

from MGT.Common import *

class GFF3Graphics:
    """Create graphical diagrams from GFF3 files"""
    
    def __init__(self,outFormat="pdf",exe=None,conf=None,width=800,noFileName=False):
        """Initialize the object setting up paths and common options"""
        if exe is None:
            exe = options.genomeTools.exe
        if conf is None:
            conf = options.genomeTools.sketchConf
        self.cmd = "%s sketch -force -style %s -format %s -width %s" % (exe,conf,outFormat,width)
        if noFileName:
            self.cmd += " -flattenfiles"

    def __call__(self,gffFile,picFile):
        """Conver GFF3 to graphics"""
        run(self.genCmd(picFile,gffFile))

    def genCmd(self,gffFile,picFile):
        return "%s %s %s" % (self.cmd,picFile,gffFile)

