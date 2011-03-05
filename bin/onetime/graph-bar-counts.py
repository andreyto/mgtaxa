"""
Make a horizontal bar histogram of assigned clade counts.
The area of each bar is equal to the corresponding count.
"""
from pylab import *
import numpy as n
import pdb
import os

inpFile = os.path.join(os.environ["MGT_HOME"],"example_data/graph-bar-counts.txt")

data = [ (words[0],int(words[1])) for words in \
    [ line.strip().split('\t') for line in open(inpFile,'r') ][1:] ][:40]

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
fig = figure(figsize=(4, 5), dpi=150)
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
#hide vertical ticks - they have no meaning
setp(ax.get_yticklabels(), visible=False)
# run through all lines drawn for xticks
for i, line in enumerate(ax.get_xticklines()):
    if i%2 == 1:   # odd indices
        line.set_visible(False)
setp(ax.get_yticklines(),visible=False)
ax.set_xlabel("# scaffolds")
ax.set_ylabel("genus")
state = dict(chLab=False)

def on_draw(event):
    if not state["chLab"]:
        #make x labels to reflect the area of bars instead of their linear size
        ax.set_xticklabels([ int(x.get_text())**2 for x in ax.get_xticklabels() ])
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
show()
#this makes visible changes made by on_draw() event handler
#we cannot call it from on_draw() because it would cause infinite recursion
fig.canvas.draw()
for format in ("svg","ps","png","eps"):
    fig.savefig("gos_scaf_counts."+format,format=format,dpi=600)

