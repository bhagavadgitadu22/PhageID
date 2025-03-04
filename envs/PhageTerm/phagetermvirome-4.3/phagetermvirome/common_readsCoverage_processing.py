## @file common_readsCoverage_processing.py
#
# VL: here I gathered functions that are common to both GPU and mono/multi CPU versions.
# These functions are called after the mapping is done and all the counters are filled from mapping output results.
from __future__ import print_function

from time import gmtime, strftime
import heapq
import itertools

import numpy as np
import pandas as pd
# Statistics
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests
from sklearn.tree import DecisionTreeRegressor 

from utilities import checkReportTitle
from SeqStats import SeqStats

import os


k_no_match_for_contig=1

def wholeCov(whole_coverage, gen_len):
    """Calculate the coverage for whole read alignments and its average"""
    if gen_len == 0:
        return whole_coverage, 1
    total_cov = sum(whole_coverage[0]) + sum(whole_coverage[1])
    ave_whole_cov = float(total_cov) / (2 * float(gen_len))
    added_whole_coverage = [x + y for x, y in zip(whole_coverage[0], whole_coverage[1])]
    return added_whole_coverage, ave_whole_cov

def testwholeCov(added_whole_coverage, ave_whole_cov, test):
    """Return information about whole coverage."""
    if test:
        return ""
    if ave_whole_cov < 50:
        print("\nWARNING: average coverage is under the limit of the software (50)")
    elif ave_whole_cov < 200:
        print("\nWARNING: average coverage is low (<200), Li's method is presumably unreliable\n")
    drop_cov = []
    start_pos = last_pos = count_pos = 0
    for pos in range(len(added_whole_coverage)):
        if added_whole_coverage[pos] < (ave_whole_cov / 1.5):
            if pos == last_pos+1:
                count_pos += 1
                last_pos = pos
            else:
                if count_pos > 100:
                    drop_cov.append( (start_pos,last_pos+1) )
                last_pos = start_pos = pos
                count_pos = 0
            last_pos = pos
    return drop_cov

def maxPaired(paired_whole_coverage, whole_coverage):
    """Max paired coverage using whole coverage, counter edge effect with paired-ends."""
    pwc = paired_whole_coverage[:]
    wc = whole_coverage[:]
    for i in range(len(pwc)):
        for j in range(len(pwc[i])):
            if pwc[i][j] < wc[i][j]:
                pwc[i][j] = wc[i][j]
    return pwc

def replaceNormMean(norm_cov):
    """Replace the values not normalised due to covLimit by mean."""
    nc_sum = nc_count = 0
    for nc in norm_cov:
        if nc > 0:
            nc_sum += nc
            nc_count += 1
    if nc_count == 0:
        mean_nc = 0
    else:
        mean_nc = nc_sum / float(nc_count)
    for i in range(len(norm_cov)):
        if norm_cov[i] == 0:
            norm_cov[i] = mean_nc
    return norm_cov, mean_nc

def normCov(termini_coverage, whole_coverage, covLimit, edge):
    """Return the termini_coverage normalised by the whole coverage (% of coverage due to first base)."""
    normalised_coverage = [len(termini_coverage[0])*[0], len(termini_coverage[0])*[0]]
    termini_len = len(termini_coverage[0])
    mean_nc = [1,1]
    for i in range(len(termini_coverage)):
        for j in range(len(termini_coverage[i])):
            if j < edge or j > termini_len-edge:
                continue
            if whole_coverage[i][j] >= covLimit:
                if float(whole_coverage[i][j]) != 0:
                    normalised_coverage[i][j] = float(termini_coverage[i][j]) / float(whole_coverage[i][j])
                else:
                    normalised_coverage[i][j] = 0
            else:
                normalised_coverage[i][j] = 0
        normalised_coverage[i], mean_nc[i] = replaceNormMean(normalised_coverage[i])
    return normalised_coverage, mean_nc

def RemoveEdge(tableau, edge):
    return tableau[edge:-edge]

def usedReads(coverage, tot_reads):
    """Retrieve the number of reads after alignment and calculate the percentage of reads lost."""
    used_reads = sum(coverage[0]) + sum(coverage[1])
    lost_reads = tot_reads - used_reads
    lost_perc = (float(tot_reads) - float(used_reads))/float(tot_reads) * 100
    return used_reads, lost_reads, lost_perc

### PEAK functions
def picMax(coverage, nbr_pic):
    """COORDINATES (coverage value, position) of the nbr_pic largest coverage value."""
    if coverage == [[],[]] or coverage == []:
        return "", "", ""
    picMaxPlus = heapq.nlargest(nbr_pic, zip(coverage[0], itertools.count()))
    picMaxMinus = heapq.nlargest(nbr_pic, zip(coverage[1], itertools.count()))
    TopFreqH = max(max(np.array(list(zip(*picMaxPlus))[0])), max(np.array(list(zip(*picMaxMinus))[0])))
    return picMaxPlus, picMaxMinus, TopFreqH

def RemoveClosePicMax(picMax, gen_len, nbr_base):
    """Remove peaks that are too close of the maximum (nbr_base around)"""
    if nbr_base == 0:
        return picMax[1:], [picMax[0]]
    picMaxRC = picMax[:]
    posMax = picMaxRC[0][1]
    LimSup = posMax + nbr_base
    LimInf = posMax - nbr_base
    if LimSup < gen_len and LimInf >= 0:
        PosOut = list(range(LimInf,LimSup))
    elif LimSup >= gen_len:
        TurnSup = LimSup - gen_len
        PosOut = list(range(posMax,gen_len))+list(range(0,TurnSup)) + list(range(LimInf,posMax))
    elif LimInf < 0:
        TurnInf = gen_len + LimInf
        PosOut = list(range(0,posMax))+list(range(TurnInf,gen_len)) + list(range(posMax,LimSup))
    picMaxOK = []
    picOUT = []
    for peaks in picMaxRC:
        if peaks[1] not in PosOut:
            picMaxOK.append(peaks)
        else:
            picOUT.append(peaks)
    return picMaxOK, picOUT

def addClosePic(picList, picClose, norm = 0):
    """Add coverage value of close peaks to the top peak. Remove picClose in picList if exist."""
    if norm:
        if picClose[0][0] >= 0.5:
            return picList, [picClose[0]]
    picListOK = picList[:]
    cov_add = 0
    for cov in picClose:
        cov_add += cov[0]
        picListOK[cov[1]] = 0.01
    picListOK[picClose[0][1]] = cov_add
    return picListOK, picClose

def remove_pics(arr,n):
    '''Removes the n highest values from the array'''
    arr=np.array(arr)
    pic_pos=arr.argsort()[-n:][::-1]
    arr2=np.delete(arr,pic_pos)
    return arr2

def gamma(X):
    """Apply a gamma distribution."""
    X = np.array(X, dtype=np.int64)
    v = remove_pics(X, 3)

    dist_max = float(max(v))
    if dist_max == 0:
        return np.array([1.00] * len(X))

    actual = np.bincount(v)
    fit_alpha, fit_loc, fit_beta = stats.gamma.fit(v)
    expected = stats.gamma.pdf(np.arange(0, dist_max + 1, 1), fit_alpha, loc=fit_loc, scale=fit_beta) * sum(actual)

    return stats.gamma.pdf(X, fit_alpha, loc=fit_loc, scale=fit_beta)


# STATISTICS
def test_pics_decision_tree(whole_coverage, termini_coverage, termini_coverage_norm, termini_coverage_norm_close):
    """Fits a gamma distribution using a decision tree."""
    L = len(whole_coverage[0])
    res = pd.DataFrame({"Position": np.array(range(L)) + 1, "termini_plus": termini_coverage[0],
                        "SPC_norm_plus": termini_coverage_norm[0], "SPC_norm_minus": termini_coverage_norm[1],
                        "SPC_norm_plus_close": termini_coverage_norm_close[0],
                        "SPC_norm_minus_close": termini_coverage_norm_close[1], "termini_minus": termini_coverage[1],
                        "cov_plus": whole_coverage[0], "cov_minus": whole_coverage[1]})

    res["cov"] = res["cov_plus"].values + res["cov_minus"].values

    res["R_plus"] = list(map(float, termini_coverage[0])) // np.mean(termini_coverage[0])
    res["R_minus"] = list(map(float, termini_coverage[1])) // np.mean(termini_coverage[1])

    regr = DecisionTreeRegressor(max_depth=3, min_samples_leaf=100)
    X = np.arange(L)
    X = X[:, np.newaxis]
    y = res["cov"].values
    regr.fit(X, y)

    # Predict
    y_1 = regr.predict(X)
    res["covnode"] = y_1
    covnodes = np.unique(y_1)
    thres = np.mean(whole_coverage[0]) / 2
    covnodes = [n for n in covnodes if n > thres]

    for node in covnodes:
        X = res[res["covnode"] == node]["termini_plus"].values
        res.loc[res["covnode"] == node, "pval_plus"] = gamma(X)
        X = res[res["covnode"] == node]["termini_minus"].values
        res.loc[res["covnode"] == node, "pval_minus"] = gamma(X)

    res.loc[res.pval_plus > 1, 'pval_plus'] = 1.00
    res.loc[res.pval_minus > 1, 'pval_minus'] = 1.00
    res = res.fillna(1.00)

    res['pval_plus_adj'] = multipletests(res["pval_plus"].values, alpha=0.01, method="bonferroni")[1]
    res['pval_minus_adj'] = multipletests(res["pval_minus"].values, alpha=0.01, method="bonferroni")[1]

    res = res.fillna(1.00)

    res_plus = pd.DataFrame(
        {"Position": res['Position'], "SPC_std": res['SPC_norm_plus'] * 100, "SPC": res['SPC_norm_plus_close'] * 100,
         "pval_gamma": res['pval_plus'], "pval_gamma_adj": res['pval_plus_adj']})
    res_minus = pd.DataFrame(
        {"Position": res['Position'], "SPC_std": res['SPC_norm_minus'] * 100, "SPC": res['SPC_norm_minus_close'] * 100,
         "pval_gamma": res['pval_minus'], "pval_gamma_adj": res['pval_minus_adj']})

    res_plus.sort_values("SPC", ascending=False, inplace=True)
    res_minus.sort_values("SPC", ascending=False, inplace=True)

    res_plus.reset_index(drop=True, inplace=True)
    res_minus.reset_index(drop=True, inplace=True)

    return res, res_plus, res_minus

### SCORING functions
# Li's methodology
def ratioR1(TopFreqH, used_reads, gen_len):
    """Calculate the ratio H/A (R1) = highest frequency/average frequency. For Li's methodology."""
    AveFreq = (float(used_reads)/float(gen_len)/2)
    if AveFreq == 0:
        return 0, 0
    R1 = float(TopFreqH)/float(AveFreq)
    return R1, AveFreq

def ratioR(picMax):
    """Calculate the T1/T2 = Top 1st frequency/Second higher frequency. For Li's methodology."""
    T1 = picMax[0][0]
    T2 = max(1,picMax[1][0])
    R = float(T1)/float(T2)
    return round(R)


def packMode(R1, R2, R3):
    """Make the prognosis about the phage packaging mode and termini type. For Li's methodology."""
    packmode = "OTHER"
    termini = ""
    forward = ""
    reverse = ""

    if R1 < 30:
        termini = "Absence"
        if R2 < 3:
            forward = "No Obvious Termini"
        if R3 < 3:
            reverse = "No Obvious Termini"
    elif R1 > 100:
        termini = "Fixed"
        if R2 < 3:
            forward = "Multiple-Pref. Term."
        if R3 < 3:
            reverse = "Multiple-Pref. Term."
    else:
        termini = "Preferred"
        if R2 < 3:
            forward = "Multiple-Pref. Term."
        if R3 < 3:
            reverse = "Multiple-Pref. Term."

    if R2 >= 3:
        forward = "Obvious Termini"
    if R3 >= 3:
        reverse = "Obvious Termini"

    if R2 >= 3 and R3 >= 3:
        packmode = "COS"
    if R2 >= 3 and R3 < 3:
        packmode = "PAC"
    if R2 < 3 and R3 >= 3:
        packmode = "PAC"
    return packmode, termini, forward, reverse

### PHAGE Information
def orientation(picMaxPlus, picMaxMinus):
    """Return phage termini orientation."""
    if not picMaxPlus and not picMaxMinus:
        return "NA"
    if picMaxPlus and not picMaxMinus:
        return "Forward"
    if not picMaxPlus and picMaxMinus:
        return "Reverse"
    if picMaxPlus and picMaxMinus:
        if picMaxPlus[0][0] > picMaxMinus[0][0]:
            return "Forward"
        elif picMaxMinus[0][0] > picMaxPlus[0][0]:
            return "Reverse"
        elif picMaxMinus[0][0] == picMaxPlus[0][0]:
            return "NA"


def typeCOS(PosPlus, PosMinus, nbr_lim):
    """ Return type of COS sequence."""
    if PosPlus < PosMinus and abs(PosPlus-PosMinus) < nbr_lim:
        return "COS (5')", "Lambda"
    else:
        return "COS (3')", "HK97"

def sequenceCohesive(Packmode, refseq, picMaxPlus, picMaxMinus, nbr_lim):
    """Return cohesive sequence for COS phages."""
    if Packmode != 'COS':
        return '', Packmode
    PosPlus = picMaxPlus[0][1]
    PosMinus = picMaxMinus[0][1]

    SC_class, SC_type = typeCOS(PosPlus, PosMinus, nbr_lim)

    if SC_class == "COS (5')":
        if abs(PosMinus - PosPlus) < nbr_lim:
            seqcoh = refseq[min(PosPlus, PosMinus):max(PosPlus, PosMinus) + 1]
            return seqcoh, Packmode
        else:
            seqcoh = refseq[max(PosPlus, PosMinus) + 1:] + refseq[:min(PosPlus, PosMinus)]
            return seqcoh, Packmode

    elif SC_class == "COS (3')":
        if abs(PosMinus - PosPlus) < nbr_lim:
            seqcoh = refseq[min(PosPlus, PosMinus) + 1:max(PosPlus, PosMinus)]
            return seqcoh, Packmode
        else:
            seqcoh = refseq[max(PosPlus, PosMinus) + 1:] + refseq[:min(PosPlus, PosMinus)]
            return seqcoh, Packmode
    else:
        return '', Packmode

def selectSignificant(table, pvalue, limit):
    """Return significant peaks over a limit"""
    table_pvalue = table.loc[lambda df: df.pval_gamma_adj < pvalue, :]
    table_pvalue_limit = table_pvalue.loc[lambda df: df.SPC > limit, :]
    table_pvalue_limit.reset_index(drop=True, inplace=True)
    return table_pvalue_limit

def testMu(paired, list_hybrid, gen_len, used_reads, seed, insert, phage_hybrid_coverage, Mu_threshold, hostseq):
    """Return Mu if enough hybrid reads compared to theory."""
    if hostseq == "":
        return 0, -1, -1, ""
    if paired != "" and len(insert) != 0:
        insert_mean    = sum(insert) / len(insert)
    else:
        insert_mean    = max(100, seed+10)
    Mu_limit       = ((insert_mean - seed) / float(gen_len)) * used_reads/2
    test           = 0
    Mu_term_plus   = "Random"
    Mu_term_minus  = "Random"
    picMaxPlus_Mu, picMaxMinus_Mu, TopFreqH_phage_hybrid = picMax(phage_hybrid_coverage, 1)
    picMaxPlus_Mu  = picMaxPlus_Mu[0][1]
    picMaxMinus_Mu = picMaxMinus_Mu[0][1]

    # Orientation
    if list_hybrid[0] > list_hybrid[1]:
        P_orient = "Forward"
    elif list_hybrid[1] > list_hybrid[0]:
        P_orient = "Reverse"
    else:
        P_orient = ""

    # Termini
    if list_hybrid[0] > ( Mu_limit * Mu_threshold ):
        test = 1
        pos_to_check = range(picMaxPlus_Mu+1,gen_len) + range(0,100)
        for pos in pos_to_check:
            if phage_hybrid_coverage[0][pos] >= max(1,phage_hybrid_coverage[0][picMaxPlus_Mu]/4):
                Mu_term_plus = pos
                picMaxPlus_Mu = pos
            else:
                Mu_term_plus = pos
                break
    # Reverse
    if list_hybrid[1] > ( Mu_limit * Mu_threshold ):
        test = 1
        pos_to_check = range(0,picMaxMinus_Mu)[::-1] + range(gen_len-100,gen_len)[::-1]
        for pos in pos_to_check:
            if phage_hybrid_coverage[1][pos] >= max(1,phage_hybrid_coverage[1][picMaxMinus_Mu]/4):
                Mu_term_minus = pos
                picMaxMinus_Mu = pos
            else:
                Mu_term_minus = pos
                break
    return test, Mu_term_plus, Mu_term_minus, P_orient

### DECISION Process
def decisionProcess(plus_significant, minus_significant, limit_fixed, gen_len, paired, insert, R1, list_hybrid,
                    used_reads, seed, phage_hybrid_coverage, Mu_threshold, refseq, hostseq):
    """ ."""
    P_orient = "NA"
    P_seqcoh = ""
    P_concat = ""
    P_type = "-"
    Mu_like = 0
    P_left = "Random"
    P_right = "Random"
    # 2 peaks sig.
    if not plus_significant.empty and not minus_significant.empty:
        # Multiple
        if (len(plus_significant["SPC"]) > 1 or len(minus_significant["SPC"]) > 1):
            if not (plus_significant["SPC"][0] > limit_fixed or minus_significant["SPC"][0] > limit_fixed):
                Redundant = 1
                P_left = "Multiple"
                P_right = "Multiple"
                Permuted = "Yes"
                P_class = "-"
                P_type = "-"
                return Redundant, Permuted, P_class, P_type, P_seqcoh, P_concat, P_orient, P_left, P_right, Mu_like

        dist_peak = abs(plus_significant['Position'][0] - minus_significant['Position'][0])
        dist_peak_over = abs(abs(plus_significant['Position'][0] - minus_significant['Position'][0]) - gen_len)
        P_left = plus_significant['Position'][0]
        P_right = minus_significant['Position'][0]
        # COS
        if (dist_peak <= 2) or (dist_peak_over <= 2):
            Redundant = 0
            Permuted = "No"
            P_class = "COS"
            P_type = "-"
        elif (dist_peak < 20 and dist_peak > 2) or (dist_peak_over < 20 and dist_peak_over > 2):
            Redundant = 0
            Permuted = "No"
            P_class, P_type = typeCOS(plus_significant["Position"][0], minus_significant["Position"][0], gen_len / 2)
            P_seqcoh, packstat = sequenceCohesive("COS", refseq, [
                ((plus_significant["SPC"][0]), (plus_significant["Position"][0]) - 1)], [((minus_significant["SPC"][0]),
                                                                                          (
                                                                                          minus_significant["Position"][
                                                                                              0]) - 1)], gen_len / 2)
        # DTR
        elif (dist_peak <= 1000) or (dist_peak_over <= 1000):
            Redundant = 1
            Permuted = "No"
            P_class = "DTR (short)"
            P_type = "T7"
            P_seqcoh, packstat = sequenceCohesive("COS", refseq, [
                ((plus_significant["SPC"][0]), (plus_significant["Position"][0]) - 1)], [((minus_significant["SPC"][0]),
                                                                                          (
                                                                                          minus_significant["Position"][
                                                                                              0]) - 1)], gen_len / 2)
        elif (dist_peak <= 0.1 * gen_len) or (dist_peak_over <= 0.1 * gen_len):
            Redundant = 1
            Permuted = "No"
            P_class = "DTR (long)"
            P_type = "T5"
            P_seqcoh, packstat = sequenceCohesive("COS", refseq, [
                ((plus_significant["SPC"][0]), (plus_significant["Position"][0]) - 1)], [((minus_significant["SPC"][0]),
                                                                                          (
                                                                                          minus_significant["Position"][
                                                                                              0]) - 1)], gen_len / 2)
        else:
            Redundant = 1
            Permuted = "No"
            P_class = "-"
            P_type = "-"
    # 1 peak sig.
    elif not plus_significant.empty and minus_significant.empty or plus_significant.empty and not minus_significant.empty:
        Redundant = 1
        Permuted = "Yes"
        P_class = "Headful (pac)"
        P_type = "P1"
        if paired != "":
            if R1 == 0 or len(insert) == 0:
                P_concat = 1
            else:
                P_concat = round((sum(insert) / len(insert)) / R1) - 1
        if not plus_significant.empty:
            P_left = plus_significant['Position'][0]
            P_right = "Distributed"
            P_orient = "Forward"
        else:
            P_left = "Distributed"
            P_right = minus_significant['Position'][0]
            P_orient = "Reverse"
    # 0 peak sig.
    elif plus_significant.empty and minus_significant.empty:
        Mu_like, Mu_term_plus, Mu_term_minus, P_orient = testMu(paired, list_hybrid, gen_len, used_reads, seed, insert,
                                                                phage_hybrid_coverage, Mu_threshold, hostseq)
        if Mu_like:
            Redundant = 0
            Permuted = "No"
            P_class = "Mu-like"
            P_type = "Mu"
            P_left = Mu_term_plus
            P_right = Mu_term_minus
        else:
            Redundant = 1
            Permuted = "Yes"
            P_class = "-"
            P_type = "-"

    return Redundant, Permuted, P_class, P_type, P_seqcoh, P_concat, P_orient, P_left, P_right, Mu_like

# Processes coverage values for a sequence.
def processCovValuesForSeq(refseq,hostseq,refseq_name,refseq_liste,seed,analysis_name,tot_reads,results_pos,test_run, paired,edge,host,test, surrounding,limit_preferred,limit_fixed,Mu_threshold,\
                           termini_coverage,whole_coverage,paired_whole_coverage,phage_hybrid_coverage,host_hybrid_coverage, host_whole_coverage,insert,list_hybrid,reads_tested,DR,DR_path=None):

    print("\n\nFinished calculating coverage values, the remainder should be completed rapidly\n",
    strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))

    # WHOLE Coverage : Average, Maximum and Minimum
    added_whole_coverage, ave_whole_cov = wholeCov(whole_coverage, len(refseq))
    added_paired_whole_coverage, ave_paired_whole_cov = wholeCov(paired_whole_coverage, len(refseq))
    added_host_whole_coverage, ave_host_whole_cov = wholeCov(host_whole_coverage, len(hostseq))

    drop_cov = testwholeCov(added_whole_coverage, ave_whole_cov, test_run)

    # NORM pic by whole coverage (1 base)
    if paired != "":
        #paired_whole_coverage_test = maxPaired(paired_whole_coverage, whole_coverage)
        termini_coverage_norm, mean_nc = normCov(termini_coverage, paired_whole_coverage, max(10, ave_whole_cov / 1.5),
                                                 edge)
    else:
        termini_coverage_norm, mean_nc = normCov(termini_coverage, whole_coverage, max(10, ave_whole_cov / 1.5), edge)

    # REMOVE edge
    termini_coverage[0] = RemoveEdge(termini_coverage[0], edge)
    termini_coverage[1] = RemoveEdge(termini_coverage[1], edge)
    termini_coverage_norm[0] = RemoveEdge(termini_coverage_norm[0], edge)
    termini_coverage_norm[1] = RemoveEdge(termini_coverage_norm[1], edge)
    whole_coverage[0] = RemoveEdge(whole_coverage[0], edge)
    whole_coverage[1] = RemoveEdge(whole_coverage[1], edge)
    paired_whole_coverage[0] = RemoveEdge(paired_whole_coverage[0], edge)
    paired_whole_coverage[1] = RemoveEdge(paired_whole_coverage[1], edge)
    added_whole_coverage = RemoveEdge(added_whole_coverage, edge)
    added_paired_whole_coverage = RemoveEdge(added_paired_whole_coverage, edge)
    added_host_whole_coverage = RemoveEdge(added_host_whole_coverage, edge)
    phage_hybrid_coverage[0] = RemoveEdge(phage_hybrid_coverage[0], edge)
    phage_hybrid_coverage[1] = RemoveEdge(phage_hybrid_coverage[1], edge)
    host_whole_coverage[0] = RemoveEdge(host_whole_coverage[0], edge)
    host_whole_coverage[1] = RemoveEdge(host_whole_coverage[1], edge)
    host_hybrid_coverage[0] = RemoveEdge(host_hybrid_coverage[0], edge)
    host_hybrid_coverage[1] = RemoveEdge(host_hybrid_coverage[1], edge)
    refseq = RemoveEdge(refseq, edge)
    if host != "":
        hostseq = RemoveEdge(hostseq, edge)
    gen_len = len(refseq)
    host_len = len(hostseq)
    if test == "DL":
        gen_len = 100000

    # READS Total, Used and Lost
    used_reads, lost_reads, lost_perc = usedReads(termini_coverage, reads_tested)

    # PIC Max
    picMaxPlus, picMaxMinus, TopFreqH = picMax(termini_coverage, 5)
    picMaxPlus_norm, picMaxMinus_norm, TopFreqH_norm = picMax(termini_coverage_norm, 5)
    picMaxPlus_host, picMaxMinus_host, TopFreqH_host = picMax(host_whole_coverage, 5)

    ### ANALYSIS

    ## Close Peaks
    picMaxPlus, picOUT_forw = RemoveClosePicMax(picMaxPlus, gen_len, surrounding)
    picMaxMinus, picOUT_rev = RemoveClosePicMax(picMaxMinus, gen_len, surrounding)
    picMaxPlus_norm, picOUT_norm_forw = RemoveClosePicMax(picMaxPlus_norm, gen_len, surrounding)
    picMaxMinus_norm, picOUT_norm_rev = RemoveClosePicMax(picMaxMinus_norm, gen_len, surrounding)

    termini_coverage_close = termini_coverage[:]
    termini_coverage_close[0], picOUT_forw = addClosePic(termini_coverage[0], picOUT_forw)
    termini_coverage_close[1], picOUT_rev = addClosePic(termini_coverage[1], picOUT_rev)

    termini_coverage_norm_close = termini_coverage_norm[:]
    termini_coverage_norm_close[0], picOUT_norm_forw = addClosePic(termini_coverage_norm[0], picOUT_norm_forw, 1)
    termini_coverage_norm_close[1], picOUT_norm_rev = addClosePic(termini_coverage_norm[1], picOUT_norm_rev, 1)

    ## Statistical Analysis
    picMaxPlus_norm_close, picMaxMinus_norm_close, TopFreqH_norm = picMax(termini_coverage_norm_close, 5)
    phage_norm, phage_plus_norm, phage_minus_norm = test_pics_decision_tree(paired_whole_coverage, termini_coverage,
                                                                            termini_coverage_norm,
                                                                            termini_coverage_norm_close)
    # VL: comment that since the 2 different conditions lead to the execution of the same piece of code...
    # if paired != "":
    #     phage_norm, phage_plus_norm, phage_minus_norm = test_pics_decision_tree(paired_whole_coverage, termini_coverage,
    #                                                                             termini_coverage_norm,
    #                                                                             termini_coverage_norm_close)
    # else:
    #     phage_norm, phage_plus_norm, phage_minus_norm = test_pics_decision_tree(whole_coverage, termini_coverage,
    #                                                                             termini_coverage_norm,
    #                                                                             termini_coverage_norm_close)


    ## LI Analysis
    picMaxPlus_close, picMaxMinus_close, TopFreqH = picMax(termini_coverage_close, 5)

    R1, AveFreq = ratioR1(TopFreqH, used_reads, gen_len)
    R2 = ratioR(picMaxPlus_close)
    R3 = ratioR(picMaxMinus_close)

    ArtPackmode, termini, forward, reverse = packMode(R1, R2, R3)
    ArtOrient = orientation(picMaxPlus_close, picMaxMinus_close)
    ArtcohesiveSeq, ArtPackmode = sequenceCohesive(ArtPackmode, refseq, picMaxPlus_close, picMaxMinus_close,
                                                   gen_len / 2)

    ### DECISION Process

    # PEAKS Significativity
    plus_significant = selectSignificant(phage_plus_norm, 1.0 / gen_len, limit_preferred)
    minus_significant = selectSignificant(phage_minus_norm, 1.0 / gen_len, limit_preferred)

    # DECISION
    Redundant, Permuted, P_class, P_type, P_seqcoh, P_concat, P_orient, P_left, P_right, Mu_like = decisionProcess(
        plus_significant, minus_significant, limit_fixed, gen_len, paired, insert, R1, list_hybrid, used_reads,
        seed, phage_hybrid_coverage, Mu_threshold, refseq, hostseq)

    ### Results
    if len(refseq_liste) != 1:
        if P_class == "-":
            if P_left == "Random" and P_right == "Random":
                P_class = "UNKNOWN"
            else:
                P_class = "NEW"
        DR[P_class][checkReportTitle(refseq_name[results_pos])] = {"analysis_name": analysis_name, "seed": seed,
                                                                 "added_whole_coverage": added_whole_coverage,
                                                                 "Redundant": Redundant, "P_left": P_left,
                                                                 "P_right": P_right, "Permuted": Permuted,
                                                                 "P_orient": P_orient,
                                                                 "termini_coverage_norm_close": termini_coverage_norm_close,
                                                                 "picMaxPlus_norm_close": picMaxPlus_norm_close,
                                                                 "picMaxMinus_norm_close": picMaxMinus_norm_close,
                                                                 "gen_len": gen_len, "tot_reads": tot_reads,
                                                                 "P_seqcoh": P_seqcoh,
                                                                 "phage_plus_norm": phage_plus_norm,
                                                                 "phage_minus_norm": phage_minus_norm,
                                                                 "ArtPackmode": ArtPackmode, "termini": termini,
                                                                 "forward": forward, "reverse": reverse,
                                                                 "ArtOrient": ArtOrient,
                                                                 "ArtcohesiveSeq": ArtcohesiveSeq,
                                                                 "termini_coverage_close": termini_coverage_close,
                                                                 "picMaxPlus_close": picMaxPlus_close,
                                                                 "picMaxMinus_close": picMaxMinus_close,
                                                                 "picOUT_norm_forw": picOUT_norm_forw,
                                                                 "picOUT_norm_rev": picOUT_norm_rev,
                                                                 "picOUT_forw": picOUT_forw,
                                                                 "picOUT_rev": picOUT_rev, "lost_perc": lost_perc,
                                                                 "ave_whole_cov": ave_whole_cov, "R1": R1, "R2": R2,
                                                                 "R3": R3, "host": host, "host_len": host_len,
                                                                 "host_whole_coverage": host_whole_coverage,
                                                                 "picMaxPlus_host": picMaxPlus_host,
                                                                 "picMaxMinus_host": picMaxMinus_host,
                                                                 "surrounding": surrounding, "drop_cov": drop_cov,
                                                                 "paired": paired, "insert": insert,
                                                                 "phage_hybrid_coverage": phage_hybrid_coverage,
                                                                 "host_hybrid_coverage": host_hybrid_coverage,
                                                                 "added_paired_whole_coverage": added_paired_whole_coverage,
                                                                 "Mu_like": Mu_like, "test_run": test_run,
                                                                 "P_class": P_class, "P_type": P_type,
                                                                 "P_concat": P_concat,
                                                                 "idx_refseq_in_list": results_pos}

        if DR_path!=None: # multi machine cluster mode.
            strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
            P_class_dir=os.path.join(DR_path,P_class)
            if os.path.exists(P_class_dir):
                if not os.path.isdir(P_class_dir):
                    raise RuntimeError("P_class_dir is not a directory")
            else:
                os.mkdir(P_class_dir)
            import pickle
            fic_name=os.path.join(P_class_dir,checkReportTitle(refseq_name[results_pos]))
            items_to_save=(analysis_name,seed,added_whole_coverage,Redundant,P_left,P_right,Permuted, \
                           P_orient,termini_coverage_norm_close,picMaxPlus_norm_close,picMaxMinus_norm_close, \
                           gen_len,tot_reads,P_seqcoh,phage_plus_norm,phage_minus_norm,ArtPackmode,termini, \
                           forward,reverse,ArtOrient,ArtcohesiveSeq,termini_coverage_close,picMaxPlus_close, \
                           picMaxMinus_close,picOUT_norm_forw,picOUT_norm_rev,picOUT_forw,picOUT_rev, \
                           lost_perc,ave_whole_cov,R1,R2,R3,host,host_len,host_whole_coverage,picMaxPlus_host, \
                           picMaxMinus_host,surrounding,drop_cov,paired, insert,phage_hybrid_coverage, \
                           host_hybrid_coverage,added_paired_whole_coverage,Mu_like,test_run,P_class,P_type,\
                           P_concat,results_pos)
            with open(fic_name,'wb') as f:
                pickle.dump(items_to_save,f)
            f.close()

    return SeqStats(P_class, P_left, P_right, P_type, P_orient, ave_whole_cov, phage_plus_norm, phage_minus_norm, ArtcohesiveSeq, P_seqcoh, Redundant, Mu_like, \
        added_whole_coverage, Permuted, termini_coverage_norm_close, picMaxPlus_norm_close, picMaxMinus_norm_close, gen_len, termini_coverage_close, \
        ArtPackmode, termini, forward, reverse, ArtOrient, picMaxPlus_close, picMaxMinus_close, picOUT_norm_forw, picOUT_norm_rev, picOUT_forw, picOUT_rev, \
        lost_perc, R1, R2, R3, picMaxPlus_host, picMaxMinus_host, drop_cov, added_paired_whole_coverage, P_concat)
