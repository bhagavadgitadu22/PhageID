##@file seq_processing.py
#
# This file contains functions that are used when running phageterm on multiple machines on a calculation cluster.
# @param DR_Path directory path where to put DR content.
from __future__ import print_function

from time import gmtime, strftime
import os
import numpy as np
from utilities import checkReportTitle
from readsCoverage_res import loadRCRes
from common_readsCoverage_processing import processCovValuesForSeq
#from SeqStats import SeqStats
def sum_readsCoverage_for_seq(dir_cov_res,idx_refseq,nb_pieces,inDArgs,fParms,inRawDArgs,dir_seq_res,DR_path):
    if os.path.exists(DR_path):
        if not (os.path.isdir(DR_path)):
            raise RuntimeError("DR_path must point to a directory")
    else:
        os.mkdir(DR_path)
    DR = {"Headful (pac)": {}, "COS (5')": {}, "COS (3')": {}, "COS": {}, "DTR (short)": {}, "DTR (long)": {},
          "Mu-like": {}, "UNKNOWN": {}, "NEW": {}}
    print("going to load ",nb_pieces," files for sequence ",idx_refseq)
    print(strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))
    for i in range(0,nb_pieces):
        fic_name = os.path.join(dir_cov_res, "coverage" + str(idx_refseq) + "_" + str(i)+".npz")
        print("loading file: ",fic_name)
        print(strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))
        partial_res=loadRCRes(fic_name)
        #npzfile=np.load(fic_name)
        if i == 0:
            termini_coverage = partial_res.termini_coverage
            whole_coverage = partial_res.whole_coverage
            paired_whole_coverage = partial_res.paired_whole_coverage
            phage_hybrid_coverage = partial_res.phage_hybrid_coverage
            host_hybrid_coverage = partial_res.host_hybrid_coverage
            host_whole_coverage = partial_res.host_whole_coverage
            list_hybrid = partial_res.list_hybrid
            insert = partial_res.insert
            paired_missmatch = partial_res.paired_mismatch
            reads_tested = partial_res.reads_tested
        else:
            termini_coverage += partial_res.termini_coverage
            whole_coverage += partial_res.whole_coverage
            paired_whole_coverage += partial_res.paired_whole_coverage
            phage_hybrid_coverage += partial_res.phage_hybrid_coverage
            host_hybrid_coverage += partial_res.host_hybrid_coverage
            host_whole_coverage += partial_res.host_whole_coverage
            list_hybrid += partial_res.list_hybrid
            insert += partial_res.insert
            paired_missmatch += partial_res.paired_mismatch
            reads_tested += partial_res.reads_tested

    # fic_name = os.path.join(dir_seq_res, "coverage" + str(idx_refseq))
    # np.savez(fic_name, termini_coverage=termini_coverage, whole_coverage=whole_coverage,
    #          paired_whole_coverage=paired_whole_coverage, \
    #          phage_hybrid_coverage=phage_hybrid_coverage, host_hybrid_coverage=host_hybrid_coverage, \
    #          host_whole_coverage=host_whole_coverage, list_hybrid=list_hybrid, insert=insert,
    #          paired_missmatch=np.array(paired_missmatch))
    termini_coverage = termini_coverage.tolist()
    whole_coverage = whole_coverage.tolist()
    paired_whole_coverage = paired_whole_coverage.tolist()
    phage_hybrid_coverage = phage_hybrid_coverage.tolist()
    host_hybrid_coverage = host_hybrid_coverage.tolist()
    host_whole_coverage = host_whole_coverage.tolist()
    list_hybrid = list_hybrid.tolist()

    if sum(termini_coverage[0]) + sum(termini_coverage[1]) == 0:
        no_match_file="no_natch"+str(idx_refseq)
        f=open(os.path.join(dir_seq_res,no_match_file),"w")
        f.write((checkReportTitle(seq_name[idx_refseq])))
        f.close()

    print("finished sum, calling processCovValuesForSeq")
    print(strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))
    # TODO: having so many values in input and returned is ugly and bad for readibility and maintanability. Group those who are related in structures.
    refseq = inDArgs.refseq_liste[idx_refseq]
    S_stats=processCovValuesForSeq(refseq, inDArgs.hostseq, inDArgs.refseq_name, inDArgs.refseq_liste, fParms.seed,
                            inRawDArgs.analysis_name, inRawDArgs.tot_reads, \
                            idx_refseq, fParms.test_run, inRawDArgs.paired, fParms.edge, inRawDArgs.host,
                            fParms.test, fParms.surrounding, \
                            fParms.limit_preferred, fParms.limit_fixed, fParms.Mu_threshold, termini_coverage,
                            whole_coverage, \
                            paired_whole_coverage, phage_hybrid_coverage, host_hybrid_coverage,
                            host_whole_coverage, insert, list_hybrid, reads_tested, DR,DR_path)
    #fic_name = os.path.join(dir_seq_res, "seq_stats" + str(idx_refseq))
    # S_stats.toFile(fic_name) s_stats content is only used in the case where there is only 1 sequence. I'm not interested in this case here since sum_readsCoverage_for_seq will be used for viromes.
    # so, just drop s_stat and forget it...
    # Only writing DR content to file is usefuk in the case of a virome processing over several machines on a cluster.
    print("exit sum_readsCoverage_for_seq")
    print(strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))





