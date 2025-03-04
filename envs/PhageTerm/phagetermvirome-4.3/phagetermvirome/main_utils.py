##@file main_utils.py
#
# Contains utility functions and classes for the main program.
# Aim is to make main simpler and smaller and thus improve testability (by allowing separate/independant testing of small program "subparts").

# Note about main program's options. This is to be discussed and subject to change.
# -g + --mapping_res_dir : Assume we are on a cluster. Perform mapping only and save results to files
# --mapping_res_dir+ --cov_res_dir : assume we are on a cluster. Process mapping results stored in files and puts the readsCoverage results in other files.
# Will use a job array in that case. each Phageterm will process 1 chunk for 1 sequence
from __future__ import print_function

from time import gmtime, strftime
import sys
import gzip
from optparse import OptionParser, OptionGroup
from utilities import checkReportTitle,changeCase
from IData_handling import totReads,genomeFastaRecovery

usage = """\n\nUsage: %prog -f reads.fastq -r phage_sequence.fasta [--report_title analysis_name -p reads_paired -s seed_lenght -d surrounding -t installation_test -c nbr_core -g host.fasta (warning increase process time) -l limit_multi-fasta -v virome_time]
[--mm --dir_cov_mm path_to_coverage_results -c nb_cores --core_id idx_core -p reads_paired -s seed_lenght -d surrounding -l limit_multi-fasta]
[--mm --dir_cov_mm path_to_coverage_results --dir_seq_mm path_to_sequence_results --DR_path path_to_results --seq_id index_of_sequence --nb_pieces nbr_of_read_chunks -p reads_paired -s seed_lenght -d surrounding -l limit_multi-fasta]
[--mm --DR_path path_to_results --dir_seq_mm path_to_sequence_results -p reads_paired -s seed_lenght -d surrounding -l limit_multi-fasta]

    Program: PhageTerm - Analyze phage termini and packaging mode using reads from high-throughput sequenced phage data
    Version: 4.1 (also py3_release_1)
    Contact: Julian Garneau <julian.garneau@usherbrooke.ca>
    Contact: David Bikard <david.bikard@pasteur.fr>
    Contact: Marc Monot <marc.monot@pasteur.fr>
    Contact: Veronique Legrand <vlegrand@pasteur.fr>

    You can perform a program test run upon installation using the "-t " option.
    Arguments for the -t option can be : C5, C3, DS, DL, M , H or V

    Example of test commands :
    PhageTerm.py -t C5       -> Test run for a 5\' cohesive end (e.g. Lambda)
    PhageTerm.py -t C3       -> Test run for a 3\' cohesive end (e.g. HK97)
    PhageTerm.py -t DS     -> Test run for a Direct Terminal Repeats end short (e.g. T7)
    PhageTerm.py -t DL     -> Test run for a Direct Terminal Repeats end long (e.g. T5)
    PhageTerm.py -t H       -> Test run for a Headful packaging (e.g. P1)
    PhageTerm.py -t M       -> Test run for a Mu-like packaging (e.g. Mu)
    PhageTerm.py -t V       -> Test run for a Virome data
    """


## checkFastaFile
#
#  Checking input Fasta file (file existence and format).
def checkFastaFile(filin):
    """Check sequence Fasta file given by user"""
    first_line = 1
    infil = gzip.open(filin, "rt") if filin.endswith(".gz") else open(filin, 'r')
    try:
        for line in infil:
            # Test '>'
            if first_line :
                if line[0] != '>':
                    return 1
                else:
                    first_line = 0
                    continue
            # Test 1st base per line : 'ATGCN>'
            base = changeCase(line[0])
            if base != 'A' and base != 'T' and base != 'C' and base != 'G' and base != 'N' and base != '\n' and base != '\r' and base != '>':
                infil.close()
                return 1
        infil.close()
        return 0
    except IOError:
        sys.exit('ERROR: No such file %s' % filin)

## setOptions
#
# Uses the OptionParser class. Defines all the options offered by phageterm and their default values if any.
# Also defines the usage message.
# Returns an optionParser object usable by the main program.
def setOptions():
    getopt = OptionParser(usage=usage)

    optreads = OptionGroup(getopt, 'Raw reads file in fastq format')
    optreads.add_option('-f', '--fastq', dest='fastq', metavar='FILE', help='Fastq reads from Illumina TruSeq')
    getopt.add_option_group(optreads)

    optref = OptionGroup(getopt, 'Phage genome in fasta format')
    optref.add_option('-r', '--ref', dest='reference', metavar='FILE',
                      help='Reference phage genome as contigs in fasta format')
    getopt.add_option_group(optref)

    optname = OptionGroup(getopt, 'Name of the phage being analyzed by the user')
    optname.add_option('--report_title', dest='analysis_name', metavar='STRING',
                       help='Manually enter the name of the analysis. Used as prefix for output file names. Default value is \"analysis_date_HHMM.')
    getopt.add_option_group(optname)

    optseed = OptionGroup(getopt, 'Lenght of the seed used for reads in the mapping process')
    optseed.add_option('-s', '--seed', dest='seed', metavar='INT', type="int",
                       help='Manually enter the lenght of the seed used for reads in the mapping process.')
    getopt.add_option_group(optseed)

    optsurround = OptionGroup(getopt, 'Lenght of the surrounding region considered for peak value cumulation')
    optsurround.add_option('-d', '--surrounding', dest='surround', type="int", metavar='INT',
                           help='Manually enter the lenght of the surrounding used to merge very close peaks in the analysis process.')
    getopt.add_option_group(optsurround)

    optcore = OptionGroup(getopt,
                          'GPU and multicore options. Default is 1 core and no GPU.')
    optcore.add_option('-c', '--core', dest='core', metavar='INT', type="int",
                       help='Manually enter the number of core you want to use.')
    getopt.add_option_group(optcore)
    #optcore.add_option('-u', '--gpu', dest='gpu', action="store_true", default=False, # VL: Keep that for later use maybe.
    #                   help='use this flag if you want to use GPU for read mapping')
    #optcore.add_option("--dir_mapping_res",dest='gpu_mapping_res_dir',metavar='STRING',default=None, help="directory where to put mapping results produced by GPU")
    # optcore.add_option("--idx_chunk",dest='idx_chunk',metavar='INT',default=None,help="index of the chunk for which we want to compute coverage")
    # optcore.add_option("--nb_chunks", dest='nb_chunks',metavar='INT', type="int",default=None,help="Indicate number of chunks wanted for GPU mapping. If None, phageTerm will automatically compute it")
    optmm=OptionGroup(getopt,"options for multi machine (or cluster mode)")
    optmm.add_option("--core_id",dest='core_id',metavar='INT',type="int",default=None,help="This option is used together with -c when running Pageterm on a cluster in parallel multimachine mode.")
    optmm.add_option("--mm",dest='multi_machine_mode',action='store_true',default=False,help="use this option to indicate that you want to use the cluster (or multi machine) mode.")
    optmm.add_option("--dir_cov_mm",dest='dir_cov_mm',metavar='STRING',default=None,help="directory where to put coverage results produced by Phageterm")
    optmm.add_option("--dir_seq_mm", dest='dir_seq_mm', metavar='STRING', default=None,
                       help="directory where to put per sequence results produced by Phageterm")
    optmm.add_option("--nb_pieces",dest='nb_pieces',metavar='INT',default=None,help="For per sequence processing after reads coverage has been done on the cluster")
    optmm.add_option("--DR_path",dest='DR_path',metavar='STRING',default=None,help="Directory where to put content of DR dictionnary (per sequence processing results)")
    optmm.add_option("--seq_id",dest='seq_id',metavar='INT',default=None,help="index of the sequence for which we want to compute coverage")
    getopt.add_option_group(optmm)

    optchk=OptionGroup(getopt,"options related to checkpoints.")
    optchk.add_option("--chk_freq",dest='chk_freq',metavar='INT',default=0,help="Frequency in minutes at which reads coverage (the longuest step in phageTerm) intermediate results must be saved ")
    optchk.add_option("--dir_chk",dest='dir_chk',metavar='STRING',default="",help="Directory where to put checkpoint files")
    getopt.add_option_group(optchk)

    opthost = OptionGroup(getopt, 'Host genome in fasta format')
    opthost.add_option('-g', '--host', dest='host', metavar='FILE',
                       help='Reference host genome as unique contig in fasta format')
    getopt.add_option_group(opthost)

    optpaired = OptionGroup(getopt, 'Use paired-end reads')
    optpaired.add_option('-p', '--paired', dest='paired', metavar='FILE',
                         help='Use paired-end reads to calculate real insert coverage')
    getopt.add_option_group(optpaired)

    optmean = OptionGroup(getopt, 'Defined phage mean coverage')
    optmean.add_option('-m', '--mean', dest='mean', metavar='INT', type="int", help='Defined phage mean coverage')
    getopt.add_option_group(optmean)

    optlimit = OptionGroup(getopt, 'Limit minimum fasta size (Default: 500)')
    optlimit.add_option('-l', '--limit', dest='limit', metavar='INT', type="int", help='Limit minimum fasta length')
    getopt.add_option_group(optlimit)

    optvirome = OptionGroup(getopt, 'Estimate execution time for a Virome')
    optvirome.add_option('-v', '--virome', dest='virome', metavar='INT', type="int",
                         help='Estimate execution time for a Virome')
    getopt.add_option_group(optvirome)

    opttest = OptionGroup(getopt, 'Perform a program test run upon installation')
    opttest.add_option('-t', '--test', dest='test', metavar='STRING',
                       help='Perform a program test run upon installation. If you want to perform a test run, use the "-t " option. Arguments for the -t option can be : C5, C3, DS, DL, H or M. C5 -> Test run for a 5\' cohesive end (e.g. Lambda); C3 -> Test run for a 3\' cohesive end (e.g. HK97); DS -> Test run for a short Direct Terminal Repeats end (e.g. T7); DL -> Test run for a long Direct Terminal Repeats end (e.g. T5); H -> Test run for a Headful packaging (e.g. P1); M -> Test run for a Mu-like packaging (e.g. Mu)')

    opttest.add_option('--nrt',dest='nrt',action='store_true',default=False,help='dump phage Class name to special file for non regression testing')
    getopt.add_option_group(opttest)

    return getopt

## User Raw data handling.
#
# This class provides encapsulation for raw data provided by the user as arguments to phageterm (input file names, testing mode if so, analysis_name, host and paired).
# It also performs checkings on the input files and computes the number of reads.
class inputRawDataArgs:
    def __init__(self,fastq,reference,host,analysis_name,paired,test,nrt):
        if test == "C5":
            print("\nPerforming a test run using test phage sequence with 5 prime cohesive overhang :")
            print("\npython PhageTerm.py -f test-data/COS-5.fastq -r test-data/COS-5.fasta --report_title TEST_cohesive_5_prime")
            fastq = "test-data/COS-5.fastq"
            reference = "test-data/COS-5.fasta"
            analysis_name = "Test-cohesive-5'"
        elif test == "C3":
            print("\nPerforming a test run using test phage sequence with 3 prime cohesive overhang:")
            print("\npython PhageTerm.py -f test-data/COS-3.fastq -r test-data/COS-3.fasta --report_title TEST_cohesive_3_prime")
            fastq = "test-data/COS-3.fastq"
            reference = "test-data/COS-3.fasta"
            analysis_name = "Test-cohesive-3'"
        elif test == "DS":
            print("\nPerforming a test run using test phage sequence with short direct terminal repeats (DTR-short) :")
            print("\npython PhageTerm.py -f test-data/DTR-short.fastq -r test-data/DTR-short.fasta --report_title TEST_short_direct_terminal_repeats")
            fastq = "test-data/DTR-short.fastq"
            reference = "test-data/DTR-short.fasta"
            analysis_name = "Test-short-direct-terminal-repeats"
        elif test == "DL":
            print("\nPerforming a test run using test phage sequence with long direct terminal repeats (DTR-long) :")
            print("\npython PhageTerm.py -f test-data/DTR-long.fastq -r test-data/DTR-long.fasta --report_title TEST_long_direct_terminal_repeats")
            fastq = "test-data/DTR-long.fastq"
            reference = "test-data/DTR-long.fasta"
            analysis_name = "Test-long-direct-terminal-repeats"
        elif test == "H":
            print("\nPerforming a test run using test phage sequence with headful packaging")
            print("\npython PhageTerm.py -f test-data/Headful.fastq -r test-data/Headful.fasta --report_title TEST_headful")
            fastq = "test-data/Headful.fastq"
            reference = "test-data/Headful.fasta"
            analysis_name = "Test-Headful"
        elif test == "M":
            print("\nPerforming a test run using test phage sequence with Mu-like packaging")
            print("\npython PhageTerm.py -f test-data/Mu-like_R1.fastq -p test-data/Mu-like_R2.fastq -r test-data/Mu-like.fasta --report_title TEST_Mu-like -g test-data/Mu-like_host.fasta")
            fastq = "test-data/Mu-like_R1.fastq"
            paired = "test-data/Mu-like_R2.fastq"
            reference = "test-data/Mu-like.fasta"
            host = "test-data/Mu-like_host.fasta"
            analysis_name = "Test-Mu-like"
        elif test == "V":
            print("\nPerforming a test run using virome data containing one example of each packaging mode")
            print("\npython PhageTerm.py -f test-data/Virome.fastq -r test-data/Virome.fasta --report_title TEST_Virome")
            fastq = "test-data/Virome.fastq"
            reference = "test-data/Virome.fasta"
            analysis_name = "Test-Virome"
        elif test==None:
            pass # Not a test, normal use.
        else:
            print("Unrecognized test run argument ('{}')!\nAllowed options are {}.".format(test, "C5, C3, DS, DL, H or M"))

        if host == None:
            host = ""
        if paired == None:
            paired = ""
        # CHECK inputs
        if analysis_name!=None:
            analysis_name = checkReportTitle(analysis_name)
            self.analysis_name = analysis_name
        else:
            self.analysis_name="NA"
        if checkFastaFile(reference):
            exit("ERROR in reference file")
        self.reference = reference
        if host != "":
            if checkFastaFile(host):
                exit("ERROR in reference file")
            self.host = host
        self.fastq=fastq
        self.paired=paired
        self.host=host
        self.nrt=nrt
        if (self.nrt==True):
            print("running nrt tests")

        # READS Number
        self.tot_reads = totReads(fastq)
        if paired != "":
            self.tot_reads_paired = totReads(paired)
            if (self.tot_reads != self.tot_reads_paired):
                print("\nWARNING: Number of reads between the two reads files differ, using single reads only\n")
                self.paired = ""


## User functional parameters handling
#
# Here gather user input parameters and global variable that define how the data will be processed from a functionnal point of view (ex: seed length...)
class functionalParms:
    def __init__(self,seed,surround,mean,limit,virome,test):
        if seed == None:
            seed = 20
        if seed < 15:
            seed = 15
        if surround == None:
            surround = 20
        self.seed=seed
        self.surrounding=surround

        if limit == None:
            limit = 500
        self.limit_reference=limit

        if virome == None:
            virome = 0
        if virome != 1:
            virome = 0
        self.virome=virome

        if mean == None:
            mean = 250
        self.mean=mean
        if test == None:
            self.test_run = 0
        else:
            self.test_run = 1
        self.test=test
        if test=="H" or test=="M" or test=="V":
            self.surrounding = 0
        if test=="V":
            self.workflow = 1
        # VARIABLE
        self.edge = 500
        self.insert_max = 1000
        self.limit_fixed = 35
        self.limit_preferred = 11
        self.Mu_threshold = 0.5
        self.draw = 0
        self.workflow = 0

## Derive other parameter from functional and raw parameters.
#
# Here, gather data derived from the rawInputData and updated according to the functionnal parameters.
# functionnal parameter workflow can also be updated.
class InputDerivedDataArgs:
    def __init__(self,inputRaw,fparms):
        # REFERENCE sequence recovery and edge adds
        self.refseq_liste, self.refseq_name, refseq_rejected = genomeFastaRecovery(inputRaw.reference, fparms.limit_reference, fparms.edge)
        self.nbr_virome = len(self.refseq_liste)
        if self.nbr_virome == 0:
            print("\nERROR: All the reference(s) sequence(s) are under the length limitation : " + str(
                fparms.limit_reference) + " (adapt your -l option)")
            exit()
        if self.nbr_virome > 1:
            fparms.workflow = 1
        length_virome = len("".join(self.refseq_liste))
        self.mean_virome = length_virome // self.nbr_virome
        if fparms.virome:
            self.refseq_liste, self.refseq_name, refseq_rejected = ["N" * int(self.mean_virome)], ["Test_virome"], 0
        if len(self.refseq_liste) == 1 and inputRaw.host != "":
            self.hostseq = genomeFastaRecovery(inputRaw.host, fparms.limit_reference, fparms.edge, 1)
            if len(self.hostseq[0]) != 0 and len(self.hostseq[0]) > len(self.refseq_liste[0]):
                print("\nHost length < Phage length : removing host sequence.")
                self.hostseq = ""
        else:
            self.hostseq = ""
            if len(self.refseq_liste) > 1:
                print("\nWARNING: Host analysis impossible with multiple fasta input\n")

## Handling of technical parameters given by the user
#
# Here gather user input parameters and former global variable that define how the data will be processed from a technical point of view (ex: multicore,gpu...)
# VL: here keep parameters related to gpu processing just in case GPU code would be needed one day for evolutions but they are not used.
class technicalParms:
    def __init__(self, core, gpu, mean, gpu_mapping_res_dir, nb_chunks, dir_cov_mm, seq_id, idx_chunk, \
                 core_id, dir_seq_mm, multi_machine_mode, DR_path, nb_pieces,chk_freq=0,dir_chk="",test_mode=False):
        self.chk_freq=chk_freq
        self.dir_chk=dir_chk
        self.multi_machine=multi_machine_mode
        self.core = core
        self.wanted_chunks = nb_chunks
        self.dir_cov_mm = dir_cov_mm
        self.DR_path=DR_path
        self.test_mode=test_mode # used for testing the checkpoint implementation.
        if nb_pieces!=None:
            self.nb_pieces=int(nb_pieces)
        else:
            self.nb_pieces =None
        if idx_chunk!=None:
            self.idx_chunk=int(idx_chunk)
        else:
            self.idx_chunk =None
        if seq_id!=None:
            self.seq_id=int(seq_id)
        else:
            self.seq_id=None
        self.core_id=core_id
        self.dir_seq_mm=dir_seq_mm
        if core == None:
            self.core = 1
        self.limit_coverage = max(50, mean * 2) / float(self.core)
        if gpu ==True and self.core > 1:
            print("Choose either multicore or gpu!")
            exit(1)
        self.gpu=gpu
        if gpu == None:
            self.gpu = False
        self.gpu_mapping_res_dir=gpu_mapping_res_dir
        if self.gpu==True and (self.dir_cov_mm != None or self.dir_seq_mm != None):
            print("when -g is used it is either to perform mapping only or whole process, --dir-cov_res/--dir_seq_res and -g are thus mutually exclusive")
            exit(1)
        if (self.gpu==True and self.core_id!=None):
            print("Inconsistency in options. -u/--gpu cannot be used with --core_id")
            exit(1)
        if self.chk_freq!=0 and self.dir_chk=="":
            print("Inconsistency in options: if frequency for checkpoints is not NULL (you activated checkpoints), you must also indicate in which directory to put them.")
            exit(1)
        if self.chk_freq==0 and self.dir_chk!="":
            print("Inconsistency in options: checkpoints are deactivated (frequency is 0) but you indicated directory for them!")
            exit(1)
        if self.multi_machine==True:
            if (self.dir_cov_mm==None and self.dir_seq_mm==None and self.DR_path==None):
                print("Please proivide path where to put results in multi machine mode")
                exit(1)
            elif  self.dir_cov_mm!=None and self.dir_seq_mm==None: # step 1: mapping+readsCoverage.
                self.checkOptConsistencyS1()
            elif  self.dir_cov_mm!=None and self.dir_seq_mm!=None: # step 2: per-sequence processing
                self.checkOptConsistencyS2()
            elif self.dir_cov_mm==None and self.dir_seq_mm!=None: # step 3: final report generation
                self.checkOptConsistencyS3()
            else:
                print("inconsistencies in options; please read documentation")
                print(usage)
                exit(1)
        else:
            if (self.dir_cov_mm!=None or self.dir_seq_mm!=None or self.DR_path!=None):
                print("Inconsistency in options: please use --mm if you intend to use multi machine mode")
                exit(1)
            if (self.chk_freq!=0 or self.dir_chk)!="":
                print("checkpoints can only be used in multi-machine mode")
                exit(0)
        ## GPU stuff, in case we need it one day
        # if (self.core>1 and self.core_id==None):
        #     if (self.gpu_mapping_res_dir!=None or self.dir_seq_res!=None or self.dir_cov_res!=None):
        #         print "Indicate core_id when processing mapping or coverage resuts on a cluster"
        #         exit(1)
        # if (self.core>1 and self.core_id!=None):
        #     if not((self.gpu_mapping_res_dir!=None and self.dir_cov_res!=None) or (self.dir_cov_res!=None and self.dir_seq_res==None)):
        #         print " Indicate both directory where to find intermediate results to process and directory where to put the results of this processing"
        #         exit(1)
        # if self.dir_cov_res!=None and (self.idx_seq!=None or self.idx_chunk!=None) and self.dir_seq_res==None and self.dir_mapping_res!=None:
        #     print "Please provide both index of sequence and chunk index. In case you have hostseq, it has index 0 by convention so --idx_seq must be >=1."
        #     exit(1)
        # if self.core<=1 and self.dir_cov_res!=None:
        #     print "Putting coverage results in files is usually used with multi-machine (cluster) mode"
        #     exit(1)

    def checkOptConsistencyS1(self):
        if self.core_id == None:
            print("Please indicate core_id when running mapping/coverage in multi machine mode")
            exit(1)
        if (self.core_id >= self.core):
            print("--core_id must be >=0 and <nb_cores")
            exit(1)
        if self.core == 1:
            print("Warning : running on only 1 core!")
        if self.DR_path != None:
            print("--DR_path is used at step 2 and step 3. It is incompatible with --dir_cov_res (step 1)")
            exit(1)
        if self.seq_id != None:
            print("--seq_id is only used at step 2. It is incompatible with --dir_cov_res (step 1)")
            exit(1)
        if self.nb_pieces != None:
            print("--nb_pieces is only used at step 2. It is incompatible with --dir_cov_res (step 1)")
            exit(1)

    def checkOptConsistencyS2(self):
        if self.DR_path == None:
            print("Please indicate DR_path when running per sequence processing in multi machine mode")
            exit(1)
        if self.seq_id == None:
            print("Please indicate index of sequence to process in multi machine mode.")
            exit(1)
        if self.nb_pieces == None:
            print(" Please indicate in how many number of packets the reads were mapped during step 1.")
            exit(1)
        if self.core_id != None:
            print("There is no need to specify --core_id doing step 2 in multi machine mode (per sequence processing of the results of step 1)")
            exit(1)
        if self.core != 1:
            print("There is no need to specify --core  doing step 2 in multi machine mode (per sequence processing of the results of step 1)")
            exit(1)

    def checkOptConsistencyS3(self):
        if self.DR_path == None:
            print("Please indicate DR_path for generating final report.")
            exit(1)
        if self.seq_id != None:
            print("--seq_id is incompatible with step 3 (report generation)")
            exit(1)
        if self.nb_pieces != None:
            print("--nb_pieces is incompatible with step 3 (report generation)")
            exit(1)
        if self.core_id != None:
            print("--core_id is incompatible with step 3 (report generation)")
            exit(1)
        if self.core != 1:
            print("--core_id is incompatible with step 3 (report generation)")
            exit(1)


## Checks options and arguments consistency and instantiates data structure for main.
#
# Consistency checkings and instantiation of technicalParms, inputDerivedDataArgs, functionalParms, inputRawDataArgs objects that are directly usable inside main.
def checkOptArgsConsistency(getopt):
    """

    :rtype:
    """
    options, arguments = getopt.parse_args()
    if options.fastq == None and options.test == None:
        getopt.error('\tNo reads file provided.\n\t\t\tUse -h or --help for more details\n')

    if options.reference == None and options.test == None:
        getopt.error('\tNo fasta reference file provided.\n\t\t\tUse -h or --help for more details\n')

    if options.analysis_name == None and options.test == None:
        analysis_name = "Analysis"

    inRawDArgs = inputRawDataArgs(options.fastq, options.reference, options.host, options.analysis_name, options.paired,
                                  options.test,options.nrt)
    fParms = functionalParms(options.seed, options.surround, options.mean, options.limit, options.virome, options.test)
    tParms = technicalParms(options.core, None, fParms.mean, None, None,
                            options.dir_cov_mm, options.seq_id, None, options.core_id,
                            options.dir_seq_mm, options.multi_machine_mode,
                            options.DR_path,options.nb_pieces,
                            float(options.chk_freq), options.dir_chk, False)
    inDArgs = InputDerivedDataArgs(inRawDArgs, fParms)
    return inRawDArgs, fParms, tParms, inDArgs  # TODO: make a version that returns only 1 structure gathering only the useful information.
