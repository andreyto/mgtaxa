"""Rescoring of the raw ICM classifition scores into confidence scores."""

from MGT.Common import *
import math

__all__ = [ "ImmScoreRescaler" ]

# As of now, the logistic regression model that predicts accuracy based on
# the winning score, sample length and taxonomic level is built in R,
# and the coefficients are used here during prediction.

class ImmScoreRescaler(object):

    def __init__(self,predScoreRescaleModel=None):
        if predScoreRescaleModel is None:
            predScoreRescaleModel = options.icm.predScoreRescaleModel
        self.predScoreRescaleModel = predScoreRescaleModel
        self.len_samp_p = - 1/math.e

    def rescale(self,score,len_samp,i_lev_pre):
        ##apply the logistic regression model
        b = self.predScoreRescaleModel["lreg"]["coef"]
        score_b = score/len_samp
        len_samp_pre = len_samp**self.len_samp_p
        y =  b["(Intercept)"]
        y += b["score_b"]*score_b
        y += b["len_samp.pre"]*len_samp_pre
        y += b["i_lev_per.pre"]*n.select([i_lev_pre!=90],[i_lev_pre],115)
        y += b["score_b:len_samp.pre"]*score_b*len_samp_pre
        y *= -1
        return 1./(1.+n.exp(y))

