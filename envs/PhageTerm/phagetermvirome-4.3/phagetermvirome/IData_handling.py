## @file IData_handling.py
#
# VL: Gather here the classes and functions useful for handling input data.
from __future__ import print_function

import gzip
from utilities import reverseComplement,changeCase
from time import gmtime, strftime
import datetime

try:
    import cPickle as pickle
except ImportError:  # python 3.x
    import pickle


## This class encapsulates the reference sequences, the host sequence if any and all useful information about the sequences.
#
# It is used both for searching the read extracts in the sequences and for exploiting the results
class refData:
    def __init__(self,refseq_list,seed,hostseq):
        self.refseq_list=refseq_list
        self.seed=seed
        self.hostseq=hostseq
        if hostseq!="":
            self.refseq_list.insert(0,hostseq)
        self.nb_sequences=len(refseq_list)

    def getIdxSeq(self,refseq):
        idx=-1
        found=False
        for s in self.refseq_list:
            idx += 1
            if s==refseq:
                found=True
                break
        if not found:
            raise RuntimeError("Couldn't find sequence in list of ref sequences.")
        return idx


    def IdxIsHostseq(self,idx_seq):
        if (((self.hostseq == "") and (idx_seq <= self.nb_sequences - 1)) or (
            (self.hostseq != "") and (idx_seq >0))):
            return False
        return True

    def getSeqSizesList(self):
        seq_sizes_list = []
        for seq in self.refseq_list:
            seq_sizes_list.append(len(seq))
        return seq_sizes_list


## Base class for handling read extracts.
#
# This class should not be used directly.
class ReadExtracts:
    def __init__(self,seed):
        self.seed = seed
        self.r_extracts_list = []
        self.nb_reads = 0
        self.nb_extracts=0

    ## Returns the list of read extracts from the loaded dataset, the number of reads and the total number of extracts
    def getRExtracts(self):
        return self.r_extracts_list,self.nb_reads,self.nb_extracts

## Class containing all the read extracts (PE reads) that must be mapped against a sequence.
class readExtractsPE(ReadExtracts):
    def __init__(self,seed):
        self.__init__(seed)


    def addRead(self, whole_PE1,whole_PE2):
        self.r_extracts_list.append(whole_PE1[:self.seed])
        self.r_extracts_list.append(whole_PE1[-self.seed:])
        self.r_extracts_list.append(whole_PE2[:self.seed])
        self.r_extracts_list.append(reverseComplement(whole_PE2)[:self.seed])
        self.r_extracts_list.append(reverseComplement(whole_PE2)[-self.seed:])
        self.nb_reads += 1
        self.nb_extracts += 5  # Number of extracts per read: 2 extracts for PE1 and 3 for PE2.



## Class containing all the read extracts (single reads) that must be mapped against a sequence.
class readsExtractsS(ReadExtracts):
    def __init__(self,seed):
        ReadExtracts.__init__(self,seed)

    ## Adds a read to the list of extracts
    #
    # @param whole_read The read as extracted from the fastq file
    # @param no_pair This paramenter s only used to make the distinction between Single and paired.
    # Note VL: I didn't use meta programming here because I thought that it would have a negative impact on performance.
    # TODO: test it when all the rest works.
    def addRead(self,whole_read,no_pair=""):
        read_part = whole_read[:self.seed]
        self.r_extracts_list.append(read_part)
        self.r_extracts_list.append(whole_read[-self.seed:])
        self.r_extracts_list.append(reverseComplement(whole_read)[:self.seed])
        self.r_extracts_list.append(reverseComplement(whole_read)[-self.seed:])
        self.nb_reads+=1
        self.nb_extracts += 4

## use objects of this class to store read offset (PE1 and PE2) in files.
class ReadInfo:
    def __init__(self, off_PE1,whole_read,seed,off_PE2=None):
        self.offset1=off_PE1
        self.offset2=off_PE2
        self.corlen = len(whole_read) - seed

## Gets the number of reads in the fastq file
# def getNbReads(fastq):
#     with open(fastq) as f:
#         for i, l in enumerate(f):
#             pass
#     nb_r=i+1
#     nb_r=nb_r/4
#     return nb_r



## loads a chunk of reads for mapping on GPU.
# Yields a ReadExtracts object plus a dictionnary of ReadInfo.
# keeps in memory the parsing state of the file.
# @param ch_size is in number of reads
# @reset_ids indicates whether or not we want read index to be reset to 0 at the beginning of each chunk.
def getChunk(fastq,seed,paired,ch_siz,reset_ids=True):
    new_chunk = False
    d_rinfo=dict()
    idx_read=0
    off2=None
    filin = open(fastq)
    line = filin.readline()
    read_paired=""
    if paired != "":
        RE=readExtractsPE(seed)
        filin_paired = open(paired)
        line_paired = filin_paired.readline()
    else:
        RE=readsExtractsS(seed)

    start = False
    num_line=0
    while line:
        # Read sequence
        read = line.split("\n")[0].split("\r")[0]
        if paired != "":
            read_paired = line_paired.split("\n")[0].split("\r")[0]
        if (read[0] == '@' and num_line%4 == 0): # make sure we don't take into account a quality score instead of a read.
            start = True
            off1=filin.tell()
            line = filin.readline()
            if paired != "":
                off2=filin_paired.tell()
                line_paired = filin_paired.readline()
            continue
        if (start == True):
            start = False
            readlen = len(read)
            if readlen < seed:
                line = filin.readline()
                if paired !="":
                    line_paired = filin_paired.readline() # also skip PE2 in that case
                continue
            RE.addRead(read,read_paired)
            d_rinfo[idx_read]=ReadInfo(off1,read,seed,off2)
            if (idx_read>0 and ((idx_read+1)%(ch_siz)==0)):
                yield RE,d_rinfo
                if (reset_ids):
                    idx_read=0
                    new_chunk=True
                if paired != "":
                    RE = readExtractsPE(seed)
                else:
                    RE = readsExtractsS(seed)
                d_rinfo = dict()
            if not new_chunk:
                idx_read+=1
            else:
                new_chunk=False

        line = filin.readline()
        if paired!="":
            line_paired = filin_paired.readline()
    filin.close()
    if paired !="":
        filin_paired.close()
    yield RE, d_rinfo

## dumps a dictionnary of ReadInfo objects indexed on read index.
#
# @param d_rinfo dictionnary to dump
# @param fic filename (incl. full path) where to dump
def dump_d_rinfo(d_rinfo,fic):
    with open(fic, 'wb') as fp:
        pickle.dump(d_rinfo, fp, protocol=pickle.HIGHEST_PROTOCOL)

## Loads a dictionnary of ReadInfo objects.
def load_d_rinfo(fic):
    with open(fic, 'rb') as fp:
        d_rinfo = pickle.load(fp)
    return d_rinfo


## loads all extracts of reads into a list for processing on GPU.
#
# returns 1 or 2 readExtracts objects plus a dictionnary of ReadInfo.
def getAllReads(fastq,seed,paired):
    d_rinfo=dict()
    idx_read=0
    off2=None
    filin = open(fastq)
    line = filin.readline()
    read_paired=""

    if paired != "":
        RE=readExtractsPE(seed)
        filin_paired = open(paired)
        line_paired = filin_paired.readline()
    else:
        RE=readsExtractsS(seed)

    start = False
    num_line=0
    while line:
        # Read sequence
        read = line.split("\n")[0].split("\r")[0]
        if paired != "":
            read_paired = line_paired.split("\n")[0].split("\r")[0]
        if (read[0] == '@' and num_line%4 == 0): # make sure we don't take into account a quality score instead of a read.
            start = True
            off1=filin.tell()
            line = filin.readline()
            if paired != "":
                off2=filin_paired.tell()
                line_paired = filin_paired.readline()
            continue
        if (start == True):
            start = False
            readlen = len(read)
            if readlen < seed:
                line = filin.readline()
                if paired !="":
                    line_paired = filin_paired.readline() # also skip PE2 in that case
                continue
            RE.addRead(read,read_paired)
            d_rinfo[idx_read]=ReadInfo(off1,read,seed,off2)
            idx_read+=1

        line = filin.readline()
        if paired!="":
            line_paired = filin_paired.readline()
    filin.close()
    if paired !="":
        filin_paired.close()
    return RE,d_rinfo

## use this class to retrieve reads from fastq file.
class ReadGetter:
    ## constructor
    #
    # @param fastq Name of the fastq file that contains the read
    # @param d_rinfo A dictionnary of ReadInfo objects that contains the offset indicating where the read starts in the file.
    # @param paired The name of the file containing the PE2 (defaults to "").
    def __init__(self,fastq,d_rinfo,paired=""):
        self.filin=open(fastq)
        self.d_rinfo=d_rinfo
        self.paired=paired
        if paired!="":
            self.filinp=open(fastq)

    def getOneRead(self,idx_read):
        read_paired=""
        self.filin.seek(self.d_rinfo[idx_read].offset1)
        read=self.filin.readline()
        if self.paired!="":
            self.filinp.seek(self.d_rinfo[idx_read].offset2)
            read_paired = self.filinp.readline()
        return read,read_paired


### READS Functions
def totReads(filin):
    """Verify and retrieve the number of reads in the fastq file before alignment"""
    if filin.endswith('.gz'):
        filein = gzip.open(filin, 'rb')
    else:
        filein = open(filin, 'r')

    line = 0
    while filein.readline():
        line += 1
    seq = float(round(line / 4))
    filein.close()
    return seq

### SEQUENCE parsing function
def genomeFastaRecovery(filin, limit_reference, edge, host_test = 0):
    """Get genome sequence from fasta file"""
    print("recovering genome from: ",filin)
    print(strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime()))
    if filin == "":
        return "", "", ""

    #infile = open(filin, 'r')
    infile = gzip.open(filin, "rt") if filin.endswith(".gz") else open(filin, 'r')
    name = []
    genome = []
    genome_line = ""
    genome_rejected = 0
    for line in infile:
        if line[0] == ">":
            if name != []:
                if len(genome_line) >= limit_reference:
                    genome.append(genome_line[-edge:] + genome_line + genome_line[:edge])
                else:
                    genome_rejected += 1
                    name = name[:-1]
                genome_line = ""
            name.append(line[1:].split('\r')[0].split('\n')[0])
        else:
            genome_line += changeCase(line).replace(' ', '').split('\r')[0].split('\n')[0]

    if len(genome_line) >= limit_reference:
        genome.append(genome_line[-edge:] + genome_line + genome_line[:edge])
        genome_line = ""
    else:
        genome_rejected += 1
        name = name[:-1]

    infile.close()

    if host_test:
        return "".join(genome)
    else:
        return genome, name, genome_rejected
    close(filin)

