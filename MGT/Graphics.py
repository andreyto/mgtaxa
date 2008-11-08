"""Support for graphical display of various machine learning methods"""
from MGT.TaxaTree import *

class LabelMapper:
    def __init__(self,taxaTree):
        self.taxaTree = taxaTree
        self.topNodes = [ taxaTree.getNode(id) for id in (myovirTaxid,cyanobacTaxid,phageTailedTaxid,)+micVirTaxids ]
        
    def label(self,id):
        (idPref,idSuf) = id.rsplit('_',1)
        if idPref == 'NCBIVM':
            node = self.taxaTree.getNode(int(idSuf))
            idTopNode = node.whichSupernode(self.topNodes).id
            if idTopNode == archTaxid:
                idTopNode = bacTaxid
            if idTopNode == bacTaxid:
                lab = "NCBI Microbial"
            elif idTopNode == virTaxid:
                lab = "NCBI Viral non-phage"
            elif idTopNode == myovirTaxid:
                lab = "NCBI Myoviridae"
            elif idTopNode == phageTailedTaxid:
                print node.name
                if False and ("synechococcus" in node.name.lower() or "cyano" in node.name.lower()):
                    lab = "NCBI Cyanophage"
                else:
                    lab = "NCBI Phage"
            elif idTopNode == cyanobacTaxid:
                lab = "NCBI Cyanobac"
        elif idPref.startswith('GSIOVIR'):
            lab = 'GOS Viral Fraction'
        elif idPref in ('GSIOMIC_2157','GSIOMIC_2'):
            lab = 'GOS Microbial'
        elif idPref == 'GSIOMIC_10239':
            lab = 'GOS Viral Large Fraction'
        elif idPref.startswith('READ_GSIOVIR'):
            lab = 'READ Viral Fraction'
        elif idPref.startswith('READ_GSIOSM'):
            lab = 'READ Microbial'
        else:
            lab = 'Unknown'
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
    
    def color4(self,lab,format="fullName"):
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
        elif lab == 'READ Viral Fraction':
            col = 'black'
        elif lab == 'READ Microbial':
            col = 'brown'
        elif lab == 'NCBI Cyanobac':
            col = 'pink'
        else:
            raise ValueError("Unknown label for color assignment %s" % lab)
        if format == "oneLetterName":        
            col = col[0]
        return col
    
    def color(self,lab,format="fullName"):
        if lab in ('NCBI Microbial','GOS Microbial','READ Microbial'):
            col = 'yellow'
        elif lab in ('NCBI Viral non-phage','NCBI Phage','GOS Viral Fraction','GOS Viral Large Fraction','READ Viral Fraction'):
            col = 'blue'
        elif lab == 'NCBI Cyanobac':
            col = 'red'
        elif lab == 'NCBI Myoviridae':
            col = 'green'
        else:
            raise ValueError("Unknown label for color assignment %s" % lab)
        if format == "oneLetterName":        
            col = col[0]
        return col
