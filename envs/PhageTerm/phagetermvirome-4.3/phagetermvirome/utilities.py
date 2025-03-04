## @file utilities.py
#
# Gather here utility methods for phageterm. Used in both CPU and GPU version.
#from string import maketrans
import re
import random
import sys

import numpy as np
import datetime

if sys.version_info < (3,):
    import string
    TRANSTAB = string.maketrans("ACGTN", "TGCAN")
else:
    TRANSTAB = str.maketrans("ACGTN", "TGCAN")

def checkReportTitle(report_title):
    """Normalise report title (take out any special char)"""
    default_title="Analysis_"
    right_now=datetime.datetime.now()
    default_title+=str(right_now.month)
    default_title+=str(right_now.day)
    default_title+="_"
    default_title+=str(right_now.hour)
    default_title+=str(right_now.minute)
    titleNorm = ""
    charok = list(range(48,58)) + list(range(65,91)) + list(range(97,123)) + [45,95]
    for char in report_title:
        if ord(char) in charok:
            titleNorm += char
    if len(titleNorm) > 1:
        return titleNorm[:20]
    else:
        return default

### SEQUENCE manipulation function
def changeCase(seq):
    """Change lower case to UPPER CASE for a sequence string."""
    return seq.upper()


def reverseComplement(seq, transtab=str.maketrans('ATGCN', 'TACGN')):
    """Reverse Complement a sequence."""
    return changeCase(seq).translate(transtab)[::-1]

def longest_common_substring(read, refseq):
    """Longest common substring between two strings."""
    m = [[0] * (1 + len(refseq)) for i in range(1 + len(read))]
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(read)):
        for y in range(1, 1 + len(refseq)):
            if read[x - 1] == refseq[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return read[x_longest - longest: x_longest]

def hybridCoverage(read, sequence, hybrid_coverage, start, end):
    """Return hybrid coverage."""
    aligned_part_only = longest_common_substring(read, sequence[start:end])
    for i in range(start, min(len(sequence),start+len(aligned_part_only))):
        hybrid_coverage[i]+=1
    return hybrid_coverage

## Determines if readPart maps against Sequence.
#
# @param readPart A part of a read (seed characters usually)
# @param sequence (a contig)
# It choses randomly a mapping position amongst all mappings found.
# It returns 2 numbers: the start and stop position of the chosen mapping location.
def applyCoverage(readPart, sequence):
    """Return a random match of a read onto the sequence. """
    position = []
    for pos in re.finditer(readPart,sequence):
        position.append(pos)
    if len(position) > 0:
        match = random.choice(position)
        return match.start(), match.end()
    else:
        return -1, -1

def correctEdge(coverage, edge):
    """Correction of the Edge coverage. """
    correctCov = np.array([len(coverage[0])*[0], len(coverage[0])*[0]])
    End = len(coverage[0])
    covSta = range(edge)
    covEnd = range(End-edge,End)
    for i in range(len(coverage)):
        for j in range(len(coverage[i])):
            correctCov[i][j] = coverage[i][j]
        for k in covSta:
            correctCov[i][k+edge] += coverage[i][k+End-edge]
        for l in covEnd:
            correctCov[i][l-edge] += coverage[i][l-End+edge]
    return correctCov

# utility class for storing results of decisionProcess function
class DecisionProcessOutput:
    def __init__(self, Redundant, Permuted, P_class, P_type, P_seqcoh, P_concat,
                 P_orient, P_left, P_right, Mu_like):
        pass

