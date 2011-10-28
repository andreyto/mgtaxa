### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#   See COPYING file distributed along with the MGTAXA package for the
#   copyright and license terms.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
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
        elif lab == 'READ Viral Fraction':
            col = 'black'
        elif lab == 'READ Microbial':
            col = 'brown'
        elif lab == 'NCBI Cyanobac':
            col = 'pink'
        elif lab == 'NCBI Myoviridae':
            col = 'cyan'
        else:
            raise ValueError("Unknown label for color assignment %s" % lab)
        if format == "oneLetterName":        
            col = col[0]
        return col
    
    def color4(self,lab,format="fullName"):
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


def barHorizArea(data,xLabel,yLabel,outPrefix):
    """
    Make a horizontal bar histogram of assigned clade counts.
    The area of each bar is equal to the corresponding count.
    """
    #import here in case we want to set backend in the calling code,
    #which has to be done before importing pylab
    import pylab as pl
    labels = [ row[0] for row in data ]
    y = n.asarray([ row[1] for row in data ],dtype='f4')

    o = n.zeros(len(y),dtype='f4')
    area=y
    # Because we use horizontal bars, width and height roles are reversed
    width=n.sqrt(area)
    left=n.roll(n.add.accumulate(n.maximum(width,20)),1)
    left[0]=0
    left += n.arange(len(width))*10
    height=n.sqrt(area)
    bottom=o.copy()
    bottom[:]=0
    fig = pl.figure(figsize=(4, 5), dpi=600)
    #fig.subplots_adjust(right=0.5)
    ax = fig.add_subplot(111)
    recs = ax.bar(left=bottom,height=width,width=height,bottom=left,orientation='horizontal')
    txts = []
    min_bar_wid = height.min()
    txt_to_bar_dist = min_bar_wid*0.3
    for (rec,lab) in zip(recs,labels):
        txts.append(ax.text(rec.get_x()+rec.get_width()+txt_to_bar_dist,rec.get_y(),lab.upper(),ha='left', va='bottom',size=7,family="sans-serif",variant="small-caps"))
        #txts[-1].set_stretch("semi-expanded")
        #txts.append(ax.text(1.1*(rec.get_x()+rec.get_width()),rec.get_y(),lab.upper(),ha='left', va='bottom',size="x-small"))
    #hide vertical ticks and lines - they have no meaning
    pl.setp(ax.get_yticklabels(),visible=False)
    pl.setp(ax.get_yticklines(),visible=False)
    # run through all lines drawn for xticks
    for i, line in enumerate(ax.get_xticklines()):
        if i%2 == 1:   # odd indices
            line.set_visible(False)
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    state = dict(chLab=False)

    def on_draw(event):
        if not state["chLab"]:
            #make x labels to reflect the area of bars instead of their linear size
            ax.set_xticklabels([ float(x.get_text())**2 for x in ax.get_xticklabels() ],size="xx-small")
            state["chLab"] = True
        for (rec,txt,irec) in zip(recs,txts,xrange(len(recs))):
            bbox = txt.get_window_extent()
            # the figure transform goes from relative coords->pixels and we
            # want the inverse of that
            bboxi = bbox.inverse_transformed(fig.transFigure)

            #print str(bboxi)
            #t_x1,t_y1 = t.get_position()
            #x_lim = ax.get_xlim()
            if bboxi.xmax > 1 or irec == 1:
                txt.set_x(rec.get_x()+rec.get_width()-txt_to_bar_dist)
                txt.set_y(rec.get_y()+rec.get_height()*0.1)
                txt.set_ha("right")
                txt.set_color("white")
        #fig.canvas.draw()
        #return False


    #title('Clade counts', bbox={'facecolor':'0.8', 'pad':5})

    # We need to get the bounding box of the text element to move it if
    # it is outside the axes border. This is
    # only possible when the drawing backend is already active. Therefore,
    # that processing will be done inside the 'draw' event hadler.
    fig.canvas.mpl_connect('draw_event', on_draw)
    #pl.show()
    #this makes visible changes made by on_draw() event handler
    #we cannot call it from on_draw() because it would cause infinite recursion
    fig.canvas.draw()
    for format in ("svg","png"):
        fig.savefig(outPrefix+"."+format,format=format,dpi=600)

def barHorizArea2(data,xLabel,yLabel,outPrefix):
    import pylab as pl
    #data = [ (words[0],int(words[1])) for words in \
    #    [ line.strip().split('\t') for line in open(inpFile,'r') ][1:] ][:40]
    # data = [ row[:2] for row in data ]

    #labels = [ row[0] for row in data ]
    #y1 = n.asarray([ row[1] for row in data ],dtype='f4')
    #y2 = n.asarray([ row[2] for row in data ],dtype='f4')

    labels = []
    color = []
    group = []
    y = []
    for _group,(_lab,(_y1,_y2)) in enumerate(data):
        labels.append(_lab)
        labels.append("")
        color.append("y")
        color.append("c")
        group.append(_group)
        group.append(_group)
        y.append(_y1)
        y.append(_y2)
    y = n.asarray(y,dtype=float)
    group = n.asarray(group,dtype=int)

    o = n.zeros(len(y),dtype='f4')
    area=y
    # Because we use horizontal bars, width and height roles are reversed
    width=n.sqrt(area)
    #TMP:
    left=n.roll(n.add.accumulate(n.maximum(width,20)),1)
    #left=n.roll(n.add.accumulate(n.maximum(width,1)),1)
    left[0]=0
    left += 40
    #add empty padding between bars
    left += n.arange(len(width)) *10
    left += group*20
    height=n.sqrt(area)
    bottom=o.copy()
    bottom[:]=0
    fig = pl.figure(figsize=(4, 5), dpi=600) # 150
    #fig.subplots_adjust(right=0.5)
    ax = fig.add_subplot(111)
    recs = ax.bar(left=bottom,height=width,width=height,bottom=left,orientation='horizontal',color=color)
    txts = []
    min_bar_wid = height.min()
    #TMP:
    #txt_to_bar_dist = min_bar_wid*0.3
    txt_to_bar_dist = max(min_bar_wid*0.3,3)
    for group_start in xrange(0,len(recs),2):
        lab = labels[group_start]
        txt_x = max([ (rec.get_x()+rec.get_width()+txt_to_bar_dist) for rec in
            recs[group_start:group_start+2] ])
        if y[group_start] == 0:
            txt_y = recs[group_start+1].get_y()
        elif y[group_start+1] == 0:
            txt_y = recs[group_start].get_y()
        else:
            txt_y = (recs[group_start].get_y() + recs[group_start+1].get_y())/2.
        txt_y = (recs[group_start].get_y() + recs[group_start+1].get_y())/2.
        txts.append(ax.text(txt_x,txt_y,lab,ha='left', va='bottom',size=7,family="sans-serif",variant="small-caps"))
        txt_y = recs[group_start+1].get_y()
        txts.append(ax.text(txt_x,txt_y,"",ha='left', va='bottom',size=7,family="sans-serif",variant="small-caps"))
        #txts[-1].set_stretch("semi-expanded")
        #txts.append(ax.text(1.1*(rec.get_x()+rec.get_width()),rec.get_y(),lab.upper(),ha='left', va='bottom',size="x-small"))
    #hide vertical ticks - they have no meaning
    pl.setp(ax.get_yticklabels(), visible=False)
    # run through all lines drawn for xticks
    for i, line in enumerate(ax.get_xticklines()):
        if i%2 == 1:   # odd indices
            line.set_visible(False)
    pl.setp(ax.get_yticklines(),visible=False)
    ax.set_xlabel(xLabel)
    ax.set_ylabel(yLabel)
    state = dict(chLab=False)

    def on_draw(event):
        if not state["chLab"]:
            #make x labels to reflect the area of bars instead of their linear size
            ax.set_xticklabels([ int(float(x.get_text())**2) for x in ax.get_xticklabels() ],size="xx-small")
            state["chLab"] = True
        for (rec,txt,irec) in zip(recs,txts,xrange(len(recs))):
            bbox = txt.get_window_extent()
            # the figure transform goes from relative coords->pixels and we
            # want the inverse of that
            bboxi = bbox.inverse_transformed(fig.transFigure)

            #print str(bboxi)
            #t_x1,t_y1 = t.get_position()
            #x_lim = ax.get_xlim()
            if bboxi.xmax > 0.90:
                txt.set_x(txt.get_position()[0]/2)
                #txt.set_x(rec.get_x()+rec.get_width()-txt_to_bar_dist)
                #txt.set_y(rec.get_y()+rec.get_height()*0.1)
                #txt.set_ha("right")
                #txt.set_color("white")
        #fig.canvas.draw()
        #return False


    #title('Clade counts', bbox={'facecolor':'0.8', 'pad':5})

    # We need to get the bounding box of the text element to move it if
    # it is outside the axes border. This is
    # only possible when the drawing backend is already active. Therefore,
    # that processing will be done inside the 'draw' event hadler.
    fig.canvas.mpl_connect('draw_event', on_draw)
    pl.show()
    #this makes visible changes made by on_draw() event handler
    #we cannot call it from on_draw() because it would cause infinite recursion
    fig.canvas.draw()
    for format in ("png",):
        fig.savefig("%s.%s" % (outPrefix,format),format=format,dpi=600)

