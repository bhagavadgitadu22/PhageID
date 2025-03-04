#! /usr/bin/env python
# -*- coding: utf-8 -*-
##@file phageterm.py
#
# main program
## PhageTerm software
#
#  Phageterm is a tool to determine phage termini and packaging strategy
#  and other useful informations using raw sequencing reads.
#  (This programs works with sequencing reads from a randomly
#  sheared DNA library preparations as Illumina TruSeq paired-end or similar)
#
#  ----------------------------------------------------------------------
#  Copyright (C) 2017 Julian Garneau
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or
#   (at your option) any later version.
#   <http://www.gnu.org/licenses/gpl-3.0.html>
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#  ----------------------------------------------------------------------
#
#  @author Julian Garneau <julian.garneau@usherbrooke.ca>
#  @author Marc Monot <marc.monot@pasteur.fr>
#  @author David Bikard <david.bikard@pasteur.fr>

 
### PYTHON Module
# Base
#import sys


from __future__ import print_function

# Multiprocessing
import multiprocessing
from multiprocessing import Manager


# Project

from utilities import checkReportTitle
from functions_PhageTerm import *
from common_readsCoverage_processing import processCovValuesForSeq,RemoveEdge
from main_utils import setOptions,checkOptArgsConsistency


### MAIN
def main():

    getopt=setOptions()
    inRawDArgs, fParms, tParms, inDArgs=checkOptArgsConsistency(getopt)

    # For each fasta in file
    DR = {"Headful (pac)":{}, "COS (5')":{}, "COS (3')":{}, "COS":{}, "DTR (short)":{}, "DTR (long)":{}, "Mu-like":{}, "UNKNOWN":{}, "NEW":{}}
    results_pos = 0
    no_match = []
    draw = 0 # used when one wants to draw some graphs.
    chk_handler = RCCheckpoint_handler(tParms.chk_freq, tParms.dir_chk, tParms.test_mode)

    if tParms.multi_machine:
        print("Running on cluster")
        print(tParms.dir_cov_mm, tParms.seq_id, tParms.dir_seq_mm, tParms.DR_path)
        if tParms.dir_cov_mm!=None and tParms.gpu_mapping_res_dir==None and tParms.dir_seq_mm==None: # perform mapping and readCoverage calculation and write results in file.
            # In that case we are processing data in an embarrassingly parallel way on a cluster.
            position = []
            read_indices = list(range(int(inRawDArgs.tot_reads)))
            part = chunks(read_indices, tParms.core)
            for i in range(tParms.core):
                position.append(next(part)[0])

            position = position + [int(inRawDArgs.tot_reads)]
            idx_refseq=chk_handler.getIdxSeq(tParms.core_id)
            print("starting processing at sequence: ",idx_refseq)
            for refseq in inDArgs.refseq_liste[idx_refseq:]:
                readsCoverage(inRawDArgs, refseq, inDArgs, fParms,None,tParms.core_id, position[tParms.core_id], position[tParms.core_id + 1],
                              tParms,chk_handler,idx_refseq)
                idx_refseq+=1
            if tParms.core_id==0:
                fname=os.path.join(tParms.dir_cov_mm,"nb_seq_processed.txt")
                f=open(fname,"w")
                f.write(str(idx_refseq))
                f.close()
            exit() # Consider that if we put results in files, it is because we are processing large datasets on a cluster.
        if tParms.dir_cov_mm!=None and tParms.seq_id!=None and tParms.dir_seq_mm!=None and tParms.DR_path!=None:
            from seq_processing import sum_readsCoverage_for_seq
            # in that case, we are processing all the results of readCoverage sequence by sequence in an embarrassingly parallel way on a cluster.
            sum_readsCoverage_for_seq(tParms.dir_cov_mm, tParms.seq_id, tParms.nb_pieces, inDArgs, fParms, inRawDArgs, tParms.dir_seq_mm,tParms.DR_path)
            exit()
        if tParms.dir_seq_mm!=None and tParms.dir_cov_mm==None and tParms.seq_id==None and tParms.DR_path!=None: # report generation
            from generate_report import loadDR,genReport
            loadDR(tParms.DR_path, DR)
            genReport(fParms, inDArgs, inRawDArgs, no_match, DR)
            exit()
    else: # mono machine original multi processing mode.
        ### COVERAGE
        print("\nCalculating coverage values, please wait (may take a while)...\n")
        start_run = time.time()

        if not fParms.test_run and tParms.core == 1:
            print("If your computer has more than 1 processor, you can use the -c or --core option to speed up the process.\n\n")


        for refseq in inDArgs.refseq_liste:
            jobs = []
            manager = Manager()
            return_dict = manager.dict()
            position = []

            read_indices = list(range(int(inRawDArgs.tot_reads)))
            part = chunks(read_indices, tParms.core)
            for i in range(tParms.core):
                position.append(next(part)[0])

            position = position + [int(inRawDArgs.tot_reads)]

            for i in range(0, tParms.core):
                tParms.core_id=i
                process = multiprocessing.Process(target=readsCoverage, args=(inRawDArgs, refseq, inDArgs, fParms,return_dict, i,position[i], position[i+1],
                                                                              tParms, chk_handler,results_pos))
                jobs.append(process)

            for j in jobs:
                j.start()

            for j in jobs:
                j.join()

            # merging results
            for core_id in range(tParms.core):
                if core_id == 0:
                    termini_coverage       = return_dict[core_id][0]
                    whole_coverage         = return_dict[core_id][1]
                    paired_whole_coverage  = return_dict[core_id][2]
                    phage_hybrid_coverage  = return_dict[core_id][3]
                    host_hybrid_coverage   = return_dict[core_id][4]
                    host_whole_coverage    = return_dict[core_id][5]
                    list_hybrid            = return_dict[core_id][6]
                    insert                 = return_dict[core_id][7].tolist()
                    paired_missmatch       = return_dict[core_id][8]
                    reads_tested           = return_dict[core_id][9]
                else:
                    termini_coverage      += return_dict[core_id][0]
                    whole_coverage        += return_dict[core_id][1]
                    paired_whole_coverage += return_dict[core_id][2]
                    phage_hybrid_coverage += return_dict[core_id][3]
                    host_hybrid_coverage  += return_dict[core_id][4]
                    host_whole_coverage   += return_dict[core_id][5]
                    list_hybrid           += return_dict[core_id][6]
                    insert                += return_dict[core_id][7].tolist()
                    paired_missmatch      += return_dict[core_id][8]
                    reads_tested          += return_dict[core_id][9]

            termini_coverage = termini_coverage.tolist()
            whole_coverage = whole_coverage.tolist()
            paired_whole_coverage = paired_whole_coverage.tolist()
            phage_hybrid_coverage = phage_hybrid_coverage.tolist()
            host_hybrid_coverage = host_hybrid_coverage.tolist()
            host_whole_coverage = host_whole_coverage.tolist()
            list_hybrid = list_hybrid.tolist()


                        # Estimate fParms.virome run time
            if fParms.virome:
                end_run = time.time()
                virome_run = int((end_run - start_run) * inDArgs.nbr_virome)
                print("\n\nThe fasta file tested contains: " + str(inDArgs.nbr_virome) + " contigs (mean length: " + str(
                    inDArgs.mean_virome) + ")")
                print("\nA complete run takes approximatively (" + str(tParms.core) + " core used) : " + EstimateTime(
                    virome_run) + "\n")
                exit()

            # Contigs without any match
            if sum(termini_coverage[0]) + sum(termini_coverage[1]) == 0:
                no_match.append((checkReportTitle(inDArgs.refseq_name[results_pos])))
                continue

            s_stats=processCovValuesForSeq(refseq,inDArgs.hostseq,inDArgs.refseq_name,inDArgs.refseq_liste,fParms.seed,inRawDArgs.analysis_name,inRawDArgs.tot_reads,\
                                       results_pos,fParms.test_run, inRawDArgs.paired,fParms.edge,inRawDArgs.host,fParms.test, fParms.surrounding,\
                                       fParms.limit_preferred,fParms.limit_fixed,fParms.Mu_threshold,termini_coverage,whole_coverage,\
                                       paired_whole_coverage,phage_hybrid_coverage,host_hybrid_coverage, host_whole_coverage,insert,list_hybrid,reads_tested,DR)


            results_pos += 1



        ### EXPORT Data
        if len(inDArgs.refseq_liste) == 1:
            # Test No Match
            if len(no_match) == 1:
                print("\n\nERROR: No reads match, please check your reference file.")
                exit()

            # Text report only
            if fParms.workflow:
                WorkflowReport(inRawDArgs.analysis_name, s_stats.P_class, s_stats.P_left, s_stats.P_right, s_stats.P_type, s_stats.P_orient, s_stats.ave_whole_cov)
            else:
                # Statistics
                ExportStatistics(inRawDArgs.analysis_name, whole_coverage, paired_whole_coverage, termini_coverage, s_stats.phage_plus_norm, s_stats.phage_minus_norm, inRawDArgs.paired, fParms.test_run)

            # Sequence
            ExportCohesiveSeq(inRawDArgs.analysis_name, s_stats.ArtcohesiveSeq, s_stats.P_seqcoh, fParms.test_run)
            ExportPhageSequence(inRawDArgs.analysis_name, s_stats.P_left, s_stats.P_right, RemoveEdge(refseq, fParms.edge), s_stats.P_orient, s_stats.Redundant, s_stats.Mu_like, \
                                s_stats.P_class, s_stats.P_seqcoh, fParms.test_run)

            # Report
            # TODO: just pass s_stat as argument; it will be cleaner.
            CreateReport(inRawDArgs.analysis_name, fParms.seed, s_stats.added_whole_coverage, draw, s_stats.Redundant, s_stats.P_left, s_stats.P_right, s_stats.Permuted, \
                         s_stats.P_orient, s_stats.termini_coverage_norm_close, \
                         s_stats.picMaxPlus_norm_close, s_stats.picMaxMinus_norm_close, s_stats.gen_len, inRawDArgs.tot_reads, s_stats.P_seqcoh, s_stats.phage_plus_norm, \
                         s_stats.phage_minus_norm, s_stats.ArtPackmode, s_stats.termini, s_stats.forward, s_stats.reverse, s_stats.ArtOrient, s_stats.ArtcohesiveSeq, \
                         s_stats.termini_coverage_close, s_stats.picMaxPlus_close, s_stats.picMaxMinus_close, \
                         s_stats.picOUT_norm_forw, s_stats.picOUT_norm_rev, s_stats.picOUT_forw, s_stats.picOUT_rev, s_stats.lost_perc, s_stats.ave_whole_cov, \
                         s_stats.R1, s_stats.R2, s_stats.R3, inRawDArgs.host, len(inDArgs.hostseq), host_whole_coverage, \
                         s_stats.picMaxPlus_host, s_stats.picMaxMinus_host, fParms.surrounding, s_stats.drop_cov, inRawDArgs.paired, insert, phage_hybrid_coverage,\
                         host_hybrid_coverage, s_stats.added_paired_whole_coverage, s_stats.Mu_like, fParms.test_run, s_stats.P_class, s_stats.P_type, s_stats.P_concat)

            if (inRawDArgs.nrt==True): # non regression tests, dump phage class name into file for later checking.
                fnrt=open("nrt.txt","w")
                fnrt.write(s_stats.P_class)
                fnrt.close()
        else:
            # Test No Match
            if len(no_match) == inDArgs.nbr_virome:
                print("\n\nERROR: No reads match, please check your reference file.")
                exit()

            # Report Resume
            multiReport     = SummaryReport(inRawDArgs.analysis_name, DR, no_match)
            multiCohSeq     = ""
            multiPhageSeq   = ""
            multiWorkflow   = "#analysis_name\tClass\tLeft\tPVal\tAdjPval\tRight\tPVal\tAdjPval\tType\tOrient\tCoverage\tComments\n"

            # No Match in workflow
            if fParms.workflow:
                for no_match_contig in no_match:
                    multiWorkflow += WorkflowReport(no_match_contig, "-", "-", "-", "-", "-", 0, 1)

            for DPC in DR:
                for DC in DR[DPC]:
                    stat_dict = DR[DPC][DC]  # splat this in everywhere
                    # Text report
                    if fParms.workflow:
                        multiWorkflow += WorkflowReport(phagename=DC, multi=1, **stat_dict)
                    # Sequence
                    idx_refseq=DR[DPC][DC]["idx_refseq_in_list"]
                    refseq=inDArgs.refseq_liste[idx_refseq]
                    multiCohSeq   += ExportCohesiveSeq(DC, stat_dict["ArtcohesiveSeq"], stat_dict["P_seqcoh"], fParms.test_run, 1)
                    multiPhageSeq += ExportPhageSequence(DC, stat_dict["P_left"], stat_dict["P_right"], RemoveEdge(refseq, fParms.edge), stat_dict["P_orient"], stat_dict["Redundant"], stat_dict["Mu_like"], stat_dict["P_class"], stat_dict["P_seqcoh"], fParms.test_run, 1)

                    # Report
                    multiReport = CreateReport(phagename=DC,
                                               draw=draw,
                                               multi=1,
                                               multiReport=multiReport,
                                               **stat_dict)

            # Workflow
            if not fParms.test:
                if fParms.workflow:
                    filoutWorkflow = open(inRawDArgs.analysis_name + "_workflow.txt", "w")
                    filoutWorkflow.write(multiWorkflow)
                    filoutWorkflow.close()

                # Concatene Sequences
                filoutCohSeq = open(inRawDArgs.analysis_name + "_cohesive-sequence.fasta", "w")
                filoutCohSeq.write(multiCohSeq)
                filoutCohSeq.close()

                filoutPhageSeq = open(inRawDArgs.analysis_name + "_sequence.fasta", "w")
                filoutPhageSeq.write(multiPhageSeq)
                filoutPhageSeq.close()

            # Concatene Report
            doc = SimpleDocTemplate("%s_PhageTerm_report.pdf" % inRawDArgs.analysis_name, pagesize=letter, rightMargin=10,leftMargin=10, topMargin=5, bottomMargin=10)
            doc.build(multiReport)


        # Real virome run time
        end_run = time.time()
        virome_run = int(end_run-start_run)
        print("\nThe fasta file tested contains: " + str(inDArgs.nbr_virome) + " contigs (mean length: " + str(inDArgs.mean_virome) + ")")
        print("The run has taken (" + str(tParms.core) + " core used) : " + EstimateTime(virome_run) + "\n")
        exit()



if __name__ == '__main__':
    main()









