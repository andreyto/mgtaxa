"""Support for graphical display of various machine learning methods"""
from MGT.TaxaTree import *

class LabelMapper:
    def __init__(self,taxaTree):
        self.taxaTree = taxaTree
        self.topNodes = [ taxaTree.getNode(id) for id in (phageTailedTaxid,)+micVirTaxids ]
        
    def label(self,id):
        (idPref,idSuf) = id.rsplit('_',1)
        if idPref == 'NCBIVM':
            idTopNode = self.taxaTree.getNode(int(idSuf)).whichSupernode(self.topNodes).id
            if idTopNode == 2157:
                idTopNode = 2
            if idTopNode == 2:
                lab = "NCBI Microbial"
            elif idTopNode == 10239:
                lab = "NCBI Viral non-phage"
            elif idTopNode == 28883:
                lab = "NCBI Phage"
        elif idPref.startswith('GSIOVIR'):
            lab = 'GOS Viral Fraction'
        elif idPref in ('GSIOMIC_2157','GSIOMIC_2'):
            lab = 'GOS Microbial'
        elif idPref == 'GSIOMIC_10239':
            lab = 'GOS Viral Large Fraction'
        else:
            lab = idPref
        return lab

    def label2(self,id):
        (idPref,idSuf) = id.rsplit('_',1)
        if idPref == 'NCBIVM':
            idTopNode = self.taxaTree.getNode(int(idSuf)).whichSupernode(self.topNodes).id
            if idTopNode == 2157:
                idTopNode = 2
            lab = "%s_%s" % (idPref,idTopNode)
        elif idPref.startswith('GSIOVIR'):
            lab = 'GSIOVIR'
        elif idPref == 'GSIOMIC_10239':
            lab = 'GSIOVIR'
        elif idPref == 'GSIOMIC_2157':
            lab = 'NCBIVM_2'
        else:
            lab = idPref
        return lab

    def color2(self,lab,format="fullName"):
        if lab == 'NCBIVM_2':
            col = 'red'
        elif lab == 'NCBIVM_10239':
            col = 'blue'
        elif lab == 'NCBIVM_28883':
            col = 'blue'
        elif lab == 'GSIOVIR':
            col = 'blue'
        elif lab == 'GSIOMIC_2':
            col = 'red'
        else:
            raise ValueError("Unknown label for color assignment %s" % lab)
        if format == "oneLetterName":        
            col = col[0]
        return col

    def color3(self,lab,format="fullName"):
        if lab == 'NCBIVM_2':
            col = 'red'
        elif lab == 'NCBIVM_10239':
            col = 'blue'
        elif lab == 'NCBIVM_28883':
            col = 'cyan'
        elif lab == 'GSIOVIR':
            col = 'green'
        elif lab == 'GSIOMIC_2':
            col = 'yellow'
        else:
            raise ValueError("Unknown label for color assignment %s" % lab)
        if format == "oneLetterName":        
            col = col[0]
        return col
    
    def color(self,lab,format="fullName"):
        if lab == 'NCBI Microbial':
            col = 'red'
        elif lab == 'NCBI Viral non-phage':
            col = 'blue'
        elif lab == 'NCBI Phage':
            col = 'cyan'
        elif lab == 'GOS Viral Fraction':
            col = 'green'
        elif lab == 'GOS Microbial':
            col = 'yellow'
        elif lab == 'GOS Viral Large Fraction':
            col = 'violet'
        else:
            raise ValueError("Unknown label for color assignment %s" % lab)
        if format == "oneLetterName":        
            col = col[0]
        return col
