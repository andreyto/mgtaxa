"Objects that describe predicted taxonomy"
from MGT.Common import *
from MGT.Svm import idDtype
from MGT.Hdf import *

class TaxaPred(object):
    """HDF5 data file that represents taxonomic predictions - one taxid per sample"""
    
    ## HDF5 group where data is stored
    group="taxaPred"
    
    def __init__(self,hdfFile,mode,**kw):
        if isinstance(hdfFile,str):
            self.hdfFile = pt.openFile(hdfFile, mode = mode,**kw)
        else:
            self.hdfFile = hdfFile

    def initForAppend(self,idSampDtype=None,idScoreDtype=None,expectedrows=None):
        if idSampDtype is None:
            idSampDtype = idDtype
        if idScoreDtype is None:
            idScoreDtype = idDtype
        f = self.hdfFile
        g = f.createGroup(f.root,self.group)
        g._v_attrs.kind = "TaxaPred"
        f.createEArray(g, 'idSamp', pt.Atom.from_dtype(n.dtype(idSampDtype)),(0,),"Sample ID",expectedrows=expectedrows)
        f.createEArray(g, 'predIdScore', pt.Atom.from_dtype(n.dtype(idScoreDtype)),(0,),
                "Predicted Score ID",expectedrows=expectedrows)
        f.createEArray(g, 'lenSamp', pt.Int64Atom(),(0,),"Sample ID",expectedrows=expectedrows)
        f.createEArray(g, 'predTaxid',pt.Int32Atom(),(0,),"Predicted Taxid",expectedrows=expectedrows)
        f.createEArray(g, 'predScore',pt.Float32Atom(),(0,),"Score of prediction",expectedrows=expectedrows)
        f.flush()
    
    def initFromSamp(self,idSamp,lenSamp,idScoreDtype=None):
        if idScoreDtype is None:
            idScoreDtype = idDtype
        assert len(idSamp) == len(lenSamp)
        f = self.hdfFile
        g = f.createGroup(f.root,self.group)
        g._v_attrs.kind = "TaxaPred"
        #@todo Maybe use CArray here and in ImmScores classes to avoid creating
        #full size temporaries just for use in array ctor. If reading of Imm
        #output is converted to streaming mode, classification will never have to load
        #all scores into RAM.
        #The line below will not work - PyTable can only create array
        #from numpy array or python sequence, but not from another tables.Array
        #f.createArray(g, 'idSamp', idSamp,"Sample ID")
        idSamp._f_copy(g)
        f.createArray(g, 'predIdScore',n.zeros(len(idSamp),dtype=n.dtype(idScoreDtype)), 
                "Predicted Score ID")
        lenSamp._f_copy(g)
        f.createArray(g, 'predTaxid',n.zeros(len(idSamp),dtype="i4"),"Predicted Taxid")
        f.createArray(g, 'predScore',n.zeros(len(idSamp),dtype="f4"),"Score of prediction")
        f.flush()
    
    def getData(self):
        return self.hdfFile.root._f_getChild(self.group)

    def close(self):
        self.hdfFile.close()

    def flush(self):
        self.hdfFile.flush()
