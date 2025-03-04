##@file SeqStats.py
#
# Utility class to store results (statistics) for a sequence once all coverage results have been processed for it

class SeqStats:
    def __init__(self,P_class, P_left, P_right, P_type, P_orient, ave_whole_cov, phage_plus_norm, phage_minus_norm, ArtcohesiveSeq,\
                 P_seqcoh, Redundant, Mu_like, added_whole_coverage, Permuted, termini_coverage_norm_close, picMaxPlus_norm_close, \
                 picMaxMinus_norm_close, gen_len, termini_coverage_close,ArtPackmode, termini, forward, reverse, ArtOrient, \
                 picMaxPlus_close, picMaxMinus_close, picOUT_norm_forw, picOUT_norm_rev, picOUT_forw, picOUT_rev, \
                 lost_perc, R1, R2, R3, picMaxPlus_host, picMaxMinus_host, drop_cov, added_paired_whole_coverage, P_concat):
        self.P_class=P_class # TODO: some information about the meaning of these fields would be welcome.
        self.P_left=P_left
        self.P_right=P_right
        self.P_type=P_type
        self.P_orient=P_orient
        self.ave_whole_cov=ave_whole_cov
        self.phage_plus_norm=phage_plus_norm
        self.phage_minus_norm=phage_minus_norm
        self.ArtcohesiveSeq=ArtcohesiveSeq
        self.P_seqcoh=P_seqcoh
        self.Redundant=Redundant
        self.Mu_like=Mu_like
        self.added_whole_coverage=added_whole_coverage
        self.Permuted=Permuted
        self.termini_coverage_norm_close=termini_coverage_norm_close
        self.picMaxPlus_norm_close=picMaxPlus_norm_close
        self.picMaxMinus_norm_close=picMaxMinus_norm_close
        self.gen_len=gen_len
        self.termini_coverage_close=termini_coverage_close
        self.ArtPackmode=ArtPackmode
        self.termini=termini
        self.forward=forward
        self.reverse=reverse
        self.ArtOrient=ArtOrient
        self.picMaxPlus_close=picMaxPlus_close
        self.picMaxMinus_close=picMaxMinus_close
        self.picOUT_norm_forw=picOUT_norm_forw
        self.picOUT_norm_rev=picOUT_norm_rev
        self.picOUT_forw=picOUT_forw
        self.picOUT_rev=picOUT_rev
        self.lost_perc=lost_perc
        self.R1=R1
        self.R2=R2
        self.R3=R3
        self.picMaxPlus_host=picMaxPlus_host
        self.picMaxMinus_host=picMaxMinus_host
        self.drop_cov=drop_cov
        self.added_paired_whole_coverage=added_paired_whole_coverage
        self.P_concat=P_concat

    def toFile(self,ficname): #TODO: implement me
        pass

# types of the elements of the class
# <type 'str'>
# <type 'numpy.int64'>
# <type 'numpy.int64'>
# <type 'str'>
# <type 'str'>
# <type 'float'>
# <class 'pandas.core.frame.DataFrame'>
# <class 'pandas.core.frame.DataFrame'>
# <type 'str'>
# <type 'str'>
# <type 'int'>
# <type 'int'>
# <type 'list'>
# <type 'str'>
# <type 'list'>
# <type 'list'>
# <type 'list'>
# <type 'int'>
# <type 'list'>
# <type 'str'>
# <type 'str'>
# <type 'str'>
# <type 'str'>
# <type 'str'>
# <type 'list'>
# <type 'list'>
# <type 'list'>
# <type 'list'>
# <type 'list'>
# <type 'list'>
# <type 'float'>
# <type 'float'>
# <type 'float'>
# <type 'float'>
# <type 'str'>
# <type 'str'>
# <type 'list'>
# <type 'list'>
# <type 'str'>
