"""
Make a horizontal bar histogram of assigned clade counts.
The area of each bar is equal to the corresponding count.
"""
from pylab import *
import numpy as n
import pdb
import os

from MGT.Sql import *


def prepDataAbundJoin(db,rank):
    """Create abundance table for ref and met sets joined by rank"""
    tbl_pref=rank if rank else "name"
    grp_fld_ref = "name_"+rank if rank else "name"
    grp_fld_met = "mgt_name_" + rank if rank else "mgt_name"
    #This works around SQLite inability to do a "full join" and
    #its buggy UNION
    db.createTableAs("%s_abund_ref" % (tbl_pref,),"""\
          SELECT %(grp_fld_ref)s AS name,
          count( * ) AS cnt
          FROM pred_annot_ref
          GROUP BY %(grp_fld_ref)s
          ORDER BY cnt DESC
          """ % dict(grp_fld_ref=grp_fld_ref))
    db.createTableAs("%s_abund_met" % (tbl_pref,),"""\
          SELECT %(grp_fld_met)s AS name,
          count( * ) AS cnt
          FROM pred_annot_met
          GROUP BY %(grp_fld_met)s
          ORDER BY cnt DESC
          """ % dict(grp_fld_met=grp_fld_met))
    db.createTableAs("%s_abund_ref_met" % (tbl_pref,),"""\
            SELECT a.name,a.cnt as cnt_ref,coalesce(b.cnt,0) as cnt_met
            FROM   %(tbl_pref)s_abund_ref a
            LEFT JOIN %(tbl_pref)s_abund_met b
            ON a.name = b.name""" % dict(tbl_pref=tbl_pref))
    db.createTableAs("%s_abund_met_ref" % (tbl_pref,),"""\
            SELECT b.name,coalesce(a.cnt,0) as cnt_ref,b.cnt as cnt_met
            FROM   %(tbl_pref)s_abund_met b
            LEFT JOIN %(tbl_pref)s_abund_ref a
            ON b.name = a.name""" % dict(tbl_pref=tbl_pref))
    db.createTableAs("%s_abund_full_join_1" % (tbl_pref,),"""\
            SELECT *
            FROM   %(tbl_pref)s_abund_ref_met""" % dict(tbl_pref=tbl_pref))
    db.ddl("""\
            INSERT INTO %(tbl_pref)s_abund_full_join_1
            SELECT * from %(tbl_pref)s_abund_met_ref""" % dict(tbl_pref=tbl_pref))
    res_tbl = "%s_abund_full_join" % (tbl_pref,)
    db.createTableAs(res_tbl,"""\
            SELECT DISTINCT *
            FROM   %(tbl_pref)s_abund_full_join_1
            ORDER BY max(cnt_ref,cnt_met) desc,(cnt_ref+cnt_met) desc""" % dict(tbl_pref=tbl_pref))
    return dict(res_tbl=res_tbl)

def prepDataAbundCross(db,rank):
    """Create abundance table for ref and met sets joined by contig ID.
    This shows the prevalent transitions of predictions for a given rank between
    ref and met sets"""
    tbl_pref=rank if rank else "name"
    grp_fld_ref = "name_"+rank if rank else "name"
    grp_fld_met = "mgt_name_" + rank if rank else "mgt_name"
    db.createTableAs("%s_pred_ref_met" % (tbl_pref,),"""\
          SELECT a.id,
          a.%(grp_fld_ref)s AS name_ref,
          b.%(grp_fld_met)s AS name_met
          FROM pred_annot_ref a, pred_annot_met b
          WHERE a.id=b.id
          """ % dict(grp_fld_ref=grp_fld_ref,grp_fld_met=grp_fld_met))
    res_tbl = "%s_pred_cross" % (tbl_pref,)
    db.createTableAs(res_tbl,"""\
          SELECT name_ref,name_met,count(*) as cnt
          FROM %(tbl_pref)s_pred_ref_met
          GROUP BY name_ref,name_met
          ORDER BY cnt DESC
          """ % dict(tbl_pref=tbl_pref))
    return dict(res_tbl=res_tbl)

    outCsv = openCompressed(opt.predOutStatsCsv,'w')
    sqlAsComment = True
    sqlpar = dict(lenMin = 5000)
    sql = """select name,count(*) as cnt from scaff_pred group by name order by cnt desc"""
    comment = "Count of assignments grouped by reference name"
    db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog="\n")
    sql = """select name_species,count(*) as cnt from scaff_pred group by name_species order by cnt desc"""
    comment = "Count of assignments grouped by reference species name"
    db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog="\n")
    sql = """select name_genus,count(*) as cnt from scaff_pred group by name_genus order by cnt desc"""
    comment = "Count of assignments grouped by reference genus name"
    db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment,epilog="\n")
    sql = """select name_family,count(*) as cnt from scaff_pred group by name_family order by cnt desc"""
    comment = "Count of assignments grouped by reference family name"
    db.exportAsCsv(sql,outCsv,comment=comment,sqlAsComment=sqlAsComment)
    outCsv.close()


def graph(data,xlabel,ylabel,outpref):
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
    for _group,(_lab,_y1,_y2) in enumerate(data):
        labels.append(_lab)
        labels.append("")
        color.append("g")
        color.append("b")
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
    #left=n.roll(n.add.accumulate(n.maximum(width,20)),1)
    left=n.roll(n.add.accumulate(n.maximum(width,1)),1)
    left[0]=0
    left += 1
    #add empty padding between bars
    #left += n.arange(len(width)) #*10
    left += group*2
    height=n.sqrt(area)
    bottom=o.copy()
    bottom[:]=0
    fig = figure(figsize=(4, 5), dpi=300) # 150
    #fig.subplots_adjust(right=0.5)
    ax = fig.add_subplot(111)
    recs = ax.bar(left=bottom,height=width,width=height,bottom=left,orientation='horizontal',color=color)
    txts = []
    min_bar_wid = height.min()
    #TMP:
    #txt_to_bar_dist = min_bar_wid*0.3
    txt_to_bar_dist = max(min_bar_wid*0.3,0.05)
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
        txts.append(ax.text(txt_x,txt_y,lab,ha='left', va='bottom',size=7,family="sans-serif",variant="small-caps"))
        txt_y = recs[group_start+1].get_y()
        txts.append(ax.text(txt_x,txt_y,"",ha='left', va='bottom',size=7,family="sans-serif",variant="small-caps"))
        #txts[-1].set_stretch("semi-expanded")
        #txts.append(ax.text(1.1*(rec.get_x()+rec.get_width()),rec.get_y(),lab.upper(),ha='left', va='bottom',size="x-small"))
    #hide vertical ticks - they have no meaning
    setp(ax.get_yticklabels(), visible=False)
    # run through all lines drawn for xticks
    for i, line in enumerate(ax.get_xticklines()):
        if i%2 == 1:   # odd indices
            line.set_visible(False)
    setp(ax.get_yticklines(),visible=False)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    state = dict(chLab=False)

    def on_draw(event):
        if not state["chLab"]:
            #make x labels to reflect the area of bars instead of their linear size
            ax.set_xticklabels([ int(float(x.get_text())**2) for x in ax.get_xticklabels() ])
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
    #show()
    #this makes visible changes made by on_draw() event handler
    #we cannot call it from on_draw() because it would cause infinite recursion
    fig.canvas.draw()
    for format in ("svg","ps","png","eps"):
        fig.savefig("ph-host-counts-%s.%s" % (outpref,format),format=format,dpi=600)


def main():
    
    #dbSqlite = os.path.join(os.environ["GOSII_WORK"],"ph-pred/asm_combined_454_large-gos-bac-comb/pred.db.sqlite")
    dbSqlite = "pred.db.sqlite"
    db = DbSqlLite(dbpath=dbSqlite)
    #for rank in ("","genus"):
    for rank in ("genus",""):
        res = prepDataAbundJoin(db=db,rank=rank)
        graph(data=db.selectAll("select * from %(res_tbl)s where not name ='null'" % res)[:20],
                xlabel="Count of assignments",
                ylabel=("Assigned %s" % (rank,)) if rank else "Lowest assignedclade",
                outpref=res["res_tbl"])
        prepDataAbundCross(db=db,rank=rank)

main()

