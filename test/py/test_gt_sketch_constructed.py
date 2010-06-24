#!/usr/bin/python
# -*- coding: utf-8 -*-

from gt.core import *
from gt.extended import *
from gt.annotationsketch import *
from gt.annotationsketch.custom_track import CustomTrack
from gt.core.gtrange import Range
import sys,os,pdb

if __name__ == "__main__":

    pngfile = "test_gt_sketch_constructed.png"
    seqid = "chromosome_21"
    nodes = []

  # construct a gene on the forward strand with two exons

    gene = FeatureNode.create_new(seqid, "gene", 100, 900, "+")
    exon = FeatureNode.create_new(seqid, "exon", 100, 200, "+")
    gene.add_child(exon)
    intron = FeatureNode.create_new(seqid, "intron", 201, 799, "+")
    gene.add_child(intron)
    exon = FeatureNode.create_new(seqid, "exon", 800, 900, "+")
    gene.add_child(exon)

  # construct a single-exon gene on the reverse strand
  # (within the intron of the forward strand gene)

    reverse_gene = FeatureNode.create_new(seqid, "gene", 400, 600, "-")
    reverse_exon = FeatureNode.create_new(seqid, "exon", 400, 600, "-")
    reverse_gene.add_child(reverse_exon)
  
  # construct a CRISPR

    crispr = FeatureNode.create_new(seqid, "CRISPR", 1000, 1300, ".")
    rep = FeatureNode.create_new(seqid, "repeat_unit", 1000, 1050, ".")
    crispr.add_child(rep)
    spa = FeatureNode.create_new(seqid, "spacer", 1051, 1100, ".")
    crispr.add_child(spa)
    rep = FeatureNode.create_new(seqid, "repeat_unit", 1101, 1150, ".")
    crispr.add_child(rep)
    spa = FeatureNode.create_new(seqid, "spacer", 1151, 1200, ".")
    crispr.add_child(spa)
    rep = FeatureNode.create_new(seqid, "repeat_unit", 1201, 1250, ".")
    crispr.add_child(rep)
    spa = FeatureNode.create_new(seqid, "spacer", 1251, 1300, ".")
    crispr.add_child(spa)

    style = Style()
    style.load_file(os.path.join(os.environ["MGT_HOME"],"etc/gt.sketch.default.style"))

    diagram = Diagram.from_array([gene, reverse_gene, crispr], Range(1, 1300),
                                 style)

    layout = Layout(diagram, 600, style)
    height = layout.get_height()
    canvas = CanvasCairoFile(style, 600, height)
    layout.sketch(canvas)

    canvas.to_file(pngfile)
    pdb.set_trace()

