from __future__ import print_function
import os
import pickle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
from reportlab.lib.pagesizes import letter, landscape
from functions_PhageTerm import SummaryReport,WorkflowReport,ExportCohesiveSeq,ExportPhageSequence,CreateReport
from common_readsCoverage_processing import RemoveEdge


def loadDR(DR_path,DR):
    for d in os.listdir(DR_path): # iterate over P_class subdirectories.
        if not os.path.isdir(os.path.join(DR_path,d)):
            err_str=DR_path+" should contain only directories."
            raise RuntimeError(err_str)
        for fic_name in os.listdir(os.path.join(DR_path,d)): # iterate over all files for a given P_class
            p=os.path.join(DR_path,d)
            fname=os.path.join(p,fic_name)
            with open(fname, 'rb') as f:
                loaded_items=pickle.load(f)
                # d is P_class name, fic_name is phagename.
                dict_tmp=dict()
                dict_tmp["phagename"]=loaded_items[0]
                dict_tmp["seed"]=loaded_items[1]
                dict_tmp["added_whole_coverage"]=loaded_items[2]
                dict_tmp["Redundant"]=loaded_items[3]
                dict_tmp["P_left"]=loaded_items[4]
                print("P_left=",dict_tmp["P_left"],type(dict_tmp["P_left"]))
                dict_tmp["P_right"] = loaded_items[5]
                print("P_right=",dict_tmp["P_right"],type(dict_tmp["P_right"]))
                dict_tmp["Permuted"]=loaded_items[6]
                dict_tmp["P_orient"] =loaded_items[7]
                dict_tmp["termini_coverage_norm_close"] =loaded_items[8]
                dict_tmp["picMaxPlus_norm_close"] =loaded_items[9]
                dict_tmp["picMaxMinus_norm_close"] =loaded_items[10]
                dict_tmp["gen_len"] =loaded_items[11]
                dict_tmp["tot_reads"] =loaded_items[12]
                dict_tmp["P_seqcoh"] =loaded_items[13]
                dict_tmp["phage_plus_norm"] =loaded_items[14]
                dict_tmp["phage_minus_norm"] =loaded_items[15]
                dict_tmp["ArtPackmode"] = loaded_items[16]
                dict_tmp["termini"] = loaded_items[17]
                dict_tmp["forward"] = loaded_items[18]
                dict_tmp["reverse"] = loaded_items[19]
                dict_tmp["ArtOrient"] = loaded_items[20]
                dict_tmp["ArtcohesiveSeq"] = loaded_items[21]
                dict_tmp["termini_coverage_close"] = loaded_items[22]
                dict_tmp["picMaxPlus_close"] = loaded_items[23]
                dict_tmp["picMaxMinus_close"] = loaded_items[24]
                dict_tmp["picOUT_norm_forw"] = loaded_items[25]
                dict_tmp["picOUT_norm_rev"] = loaded_items[26]
                dict_tmp["picOUT_forw"] = loaded_items[27]
                dict_tmp["picOUT_rev"] = loaded_items[28]
                dict_tmp["lost_perc"] = loaded_items[29]
                dict_tmp["ave_whole_cov"] = loaded_items[30]
                dict_tmp["R1"] = loaded_items[31]
                dict_tmp["R2"] = loaded_items[32]
                dict_tmp["R3"] = loaded_items[33]
                dict_tmp["host"] = loaded_items[34]
                dict_tmp["host_len"] = loaded_items[35]
                dict_tmp["host_whole_coverage"] = loaded_items[36]
                dict_tmp["picMaxPlus_host"] = loaded_items[37]
                dict_tmp["picMaxMinus_host"] = loaded_items[38]
                dict_tmp["surrounding"] = loaded_items[39]
                dict_tmp["drop_cov"] = loaded_items[40]
                dict_tmp["paired"] = loaded_items[41]
                dict_tmp["insert"] = loaded_items[42]
                dict_tmp["phage_hybrid_coverage"] = loaded_items[43]
                dict_tmp["host_hybrid_coverage"] = loaded_items[44]
                dict_tmp["added_paired_whole_coverage"] = loaded_items[45]
                dict_tmp["Mu_like"] = loaded_items[46]
                dict_tmp["test_run"] = loaded_items[47]
                dict_tmp["P_class"] = loaded_items[48]
                dict_tmp["P_type"] = loaded_items[49]
                dict_tmp["P_concat"] = loaded_items[50]
                dict_tmp["idx_refseq_in_list"] = loaded_items[51]
                DR [d][fic_name]=dict_tmp
            f.close()




def genReport(fParms,inDArgs,inRawDArgs,no_match,DR):
    # Test No Match
    if len(no_match) == inDArgs.nbr_virome:
        print("\n\nERROR: No reads match, please check your reference file.")
        exit()

    # Report Resume
    multiReport = SummaryReport(inRawDArgs.analysis_name, DR, no_match)
    multiCohSeq = ""
    multiPhageSeq = ""
    multiWorkflow = "#phagename\tClass\tLeft\tRight\tType\tOrient\tCoverage\tComments\n"

    # No Match in workflow
    if fParms.workflow:
        for no_match_contig in no_match:
            multiWorkflow += WorkflowReport(no_match_contig, "-", "-", "-", "-", "-", 0, 1)

    for DPC in DR:
        for DC in DR[DPC]:
            # Text report
            if fParms.workflow: # phagename, P_class, P_left, P_right, P_type, P_orient, ave_whole_cov, multi = 0
                multiWorkflow += WorkflowReport(DC, DR[DPC][DC]["P_class"], DR[DPC][DC]["P_left"],
                                                DR[DPC][DC]["P_right"],
                                                DR[DPC][DC]["P_type"], DR[DPC][DC]["P_orient"],
                                                DR[DPC][DC]["ave_whole_cov"], 1,DR[DPC][DC]["phage_plus_norm"],
                                                DR[DPC][DC]["phage_minus_norm"])

            # Sequence
            idx_refseq = DR[DPC][DC]["idx_refseq_in_list"]
            refseq = inDArgs.refseq_liste[idx_refseq]
            multiCohSeq += ExportCohesiveSeq(DC, DR[DPC][DC]["ArtcohesiveSeq"], DR[DPC][DC]["P_seqcoh"], fParms.test_run, 1)
            multiPhageSeq += ExportPhageSequence(DC, DR[DPC][DC]["P_left"], DR[DPC][DC]["P_right"], RemoveEdge(refseq, fParms.edge),
                                                 DR[DPC][DC]["P_orient"], DR[DPC][DC]["Redundant"], DR[DPC][DC]["Mu_like"],
                                                 DR[DPC][DC]["P_class"], DR[DPC][DC]["P_seqcoh"], fParms.test_run, 1)

            # Report
            draw=0 # TODO VL: ask what is the use of this parameter that is alwayes 0...
            multiReport = CreateReport(DC, DR[DPC][DC]["seed"], DR[DPC][DC]["added_whole_coverage"], draw,
                                       DR[DPC][DC]["Redundant"], DR[DPC][DC]["P_left"], DR[DPC][DC]["P_right"],
                                       DR[DPC][DC]["Permuted"], DR[DPC][DC]["P_orient"],
                                       DR[DPC][DC]["termini_coverage_norm_close"], DR[DPC][DC]["picMaxPlus_norm_close"],
                                       DR[DPC][DC]["picMaxMinus_norm_close"], DR[DPC][DC]["gen_len"],
                                       DR[DPC][DC]["tot_reads"], DR[DPC][DC]["P_seqcoh"], DR[DPC][DC]["phage_plus_norm"],
                                       DR[DPC][DC]["phage_minus_norm"], DR[DPC][DC]["ArtPackmode"], DR[DPC][DC]["termini"],
                                       DR[DPC][DC]["forward"], DR[DPC][DC]["reverse"], DR[DPC][DC]["ArtOrient"],
                                       DR[DPC][DC]["ArtcohesiveSeq"], DR[DPC][DC]["termini_coverage_close"],
                                       DR[DPC][DC]["picMaxPlus_close"], DR[DPC][DC]["picMaxMinus_close"],
                                       DR[DPC][DC]["picOUT_norm_forw"], DR[DPC][DC]["picOUT_norm_rev"],
                                       DR[DPC][DC]["picOUT_forw"], DR[DPC][DC]["picOUT_rev"], DR[DPC][DC]["lost_perc"],
                                       DR[DPC][DC]["ave_whole_cov"], DR[DPC][DC]["R1"], DR[DPC][DC]["R2"],
                                       DR[DPC][DC]["R3"], DR[DPC][DC]["host"], DR[DPC][DC]["host_len"],
                                       DR[DPC][DC]["host_whole_coverage"], DR[DPC][DC]["picMaxPlus_host"],
                                       DR[DPC][DC]["picMaxMinus_host"], DR[DPC][DC]["surrounding"], DR[DPC][DC]["drop_cov"],
                                       DR[DPC][DC]["paired"], DR[DPC][DC]["insert"], DR[DPC][DC]["phage_hybrid_coverage"],
                                       DR[DPC][DC]["host_hybrid_coverage"], DR[DPC][DC]["added_paired_whole_coverage"],
                                       DR[DPC][DC]["Mu_like"], fParms.test_run, DR[DPC][DC]["P_class"],
                                       DR[DPC][DC]["P_type"], DR[DPC][DC]["P_concat"], 1, multiReport)

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
    doc = SimpleDocTemplate("%s_PhageTerm_report.pdf" % inRawDArgs.analysis_name, pagesize=letter, rightMargin=10,
                            leftMargin=10, topMargin=5, bottomMargin=10)
    doc.build(multiReport)
