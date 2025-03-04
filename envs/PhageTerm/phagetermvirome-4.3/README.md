# PROGRAM
# =======

PhageTerm.py - run as command line in a shell


# VERSION
# =======

Version 4.3
Compatible with python 3.9 up to 3.13


# INTRODUCTION
# ============

PhageTermVirome software is a tool to determine phage genome termini and genome packaging mode on single phage or multiple contigs at once.
The software uses phage and virome sequencing reads obtained from libraries prepared with DNA fragmented randomly (e.g. Covaris fragmentation,
and library preparation using Illumina TruSeq). Phage or virome sequencing reads (fastq files) are aligned to the assembled phage genome or assembled
virome (fasta or multifasta files) in order to  calculate two types of coverage values (whole genome coverage and the Starting Position Coverage (SPC)). The starting position coverage is used to perform a detailed termini and packaging mode analysis. 

Mu-type phage analysis : can be done if user suspect the phage genome to be Mu-like type (Only for single phage genome analysis, not possible with multifasta file) :
User can also provide the host (bacterial) genome sequence. The Mu-type phage analysis will take the reads that does not match the phage
genome and align them on the bacterial genome using the same mapping function. The analysis to identify Mu-like phages is available only when providing a single phage genome (not possible if user provide a multi-fast file with multiple assembled phage contigs).


The previous PhageTerm program (single phage analysis only) is still available at https://sourceforge.net/projects/phageterm/ (for versions <3.0.0)


A Galaxy wrapper version is also available for the previous version at https://galaxy.pasteur.fr (only for the first version PhageTerm).
PhageTermVirome is not implemented on Galaxy yet).

Since version 4.1, PhageTerm can work in 2 modes:
- the usual mono machine mode (parallelization on several cores on the same machine). 
- a new multi machine mode (advanced users) with parallelization on several machines, using intermediate files for data exchange.
The default mode is mono machine.

Version 3.0.0 up to version 4.0 work with python 2.7

From version 4.0 to version 4.1, PhageTerm (now PhageTermVirome) works with python 3.7

version 4.1.1 works with python 3.9 up to python 3.13.


# PREREQUISITES
# =============

## For version 4.2
- python3.9
- poetry (https://python-poetry.org/docs/)
- An up to date list of dependencies can be found in the pyproject.toml file.

# INSTALLATION
# ============

- install python 3.9 (or higher; up to 1.13) and pip
- install poetry by typing

    pip install poetry

- to allow installation of deprecated package sklearn (use of sklearn will be replaced by scikit-learn in a next version), do

    export SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True

- download the .tar.gz file containing PhageTermVirome latest vesion from here : https://test.pypi.org/project/phagetermvirome/#files
and unzip it. This should create a "phagetermvirome-4.2" directory. Go into it and type:
    
    poetry install

- poetry just created a virtual environment for phageterm virome and installed all necessary dependencies in it.
Activate the virtual environment by typing:

    poetry shell

- type the following command (this is a workaround for a problem that is to be fixed in the next version)

    export PYTHONPATH=path/to/phagetermvirome

where "path/to/phagetermvirome" is the path to the phagetermvirome subdirectory in the directory where you unzipped the archive

- Type the following command to check that everything works well.

    phageterm --help



# COMMAND LINE USAGE
# ==================

Basic usage with mandatory options (PhageTermVirome needs at least one read file, but user can provide a second corresponding paired-end read file if available, using the -p option).

	phageterm -f reads.fastq -r phage_sequence(s).fasta

    
	Help:   
    
        phageterm -h
        phageterm --help


	After installation, we recommend users to perform a software run test, use any of the following:
    	-t TEST_VALUE, --test=TEST_VALUE
                    TEST_VALUE=C5   : Test run for a 5' cohesive end (e.g. Lambda)                        
               			TEST_VALUE=C3   : Test run for a 3' cohesive end (e.g. HK97)
               			TEST_VALUE=DS   : Test run for a short Direct Terminal Repeats end (e.g. T7)
               			TEST_VALUE=DL   : Test run for a long Direct Terminal Repeats end (e.g. T5)
               			TEST_VALUE=H    : Test run for a Headful packaging (e.g. P1)
               			TEST_VALUE=M    : Test run for a Mu-like packaging (e.g. Mu)


Non-mandatory options

[-p reads_paired -c nbr_core_threads --report_title name_to_write_on_report_outputs -s seed_lenght -d surrounding -g host.fasta -l contig_size_limit_multi-fasta -v virome_run_time_estimation]


Additional advanced options (only for multi-machine users)


[--mm --dir_cov_mm path_to_coverage_results -c nb_cores --core_id idx_core -p reads_paired -s seed_lenght -d surrounding -l limit_multi-fasta]
[--mm --dir_cov_mm path_to_coverage_results --dir_seq_mm path_to_sequence_results --DR_path path_to_results --seq_id index_of_sequence --nb_pieces nbr_of_read_chunks -p reads_paired -s seed_lenght -d surrounding -l limit_multi-fasta] [--mm --DR_path path_to_results --dir_seq_mm path_to_sequence_results -p reads_paired -s seed_lenght -d surrounding -l limit_multi-fasta]

    


   Detailed  ptions:


	Raw reads file in fastq format:
    -f INPUT_FILE, --fastq=INPUT_FILE
                        Fastq reads 
                        (NGS sequences from random fragmentation DNA only, 
                        e.g. Illumina TruSeq)
                        
	Phage genome(s) in fasta format:
    -r INPUT_FILE, --ref=INPUT_FILE
                        Reference phage genome(s) as unique contig in fasta format



    Other options common to both modes:

  Raw reads file in fastq format:
    -p INPUT_FILE, --paired=INPUT_FILE
                        Paired fastq reads
                        (NGS sequences from random fragmentation DNA only,
                        e.g. Illumina TruSeq)

	Analysis_name to write on output reports:
    --report_title USER_REPORT_NAME, --report_title=REPORT_NAME
                        Manually enter the name you want to have on your report outputs.
                        Used as prefix for output files.

	Lenght of the seed used for reads in the mapping process:
    -s SEED_LENGHT, --seed=SEED_LENGHT
                        Manually enter the lenght of the seed used for reads
                        in the mapping process (Default: 20).

	Number of nucleotides around the main peak to consider for merging adjacent significant peaks (set to 1 to discover secondary terminus but sites).  
    -d SUROUNDING_LENGHT, --surrounding=SUROUNDING_LENGHT
                        Manually enter the lenght of the surrounding used to
                        merge close peaks in the analysis process (Default: 20).

	Host genome in fasta format (option available only for analysis with a single phage genome):
    -g INPUT_FILE, --host=INPUT_FILE
                        Genome of reference host (bacterial genome) in fasta format
                        Warning: increase drastically process time
                        This option can be used only when analyzing a single phage genome (not available for virome contigs as multifasta)
                        
	Define phage mean coverage:
    -m MEAN_NBR, --mean=MEAN_NBR
                        Phage mean coverage to use (Default: 250).        

	Define phage mean coverage:
    -l LIMIT_FASTA, —limit=LIMIT_FASTA
                        Minimum phage fasta length (Default: 500).


    Options for mono machine (default) mode only
                
	Software run test:
    -t TEST_VALUE, --test=TEST_VALUE
                        TEST_VALUE=C5   : Test run for a 5' cohesive end (e.g. Lambda)                        
               			    TEST_VALUE=C3   : Test run for a 3' cohesive end (e.g. HK97)
               			    TEST_VALUE=DS   : Test run for a short Direct Terminal Repeats end (e.g. T7)
               			    TEST_VALUE=DL   : Test run for a long Direct Terminal Repeats end (e.g. T5)
               			    TEST_VALUE=H    : Test run for a Headful packaging (e.g. P1)
               			    TEST_VALUE=M    : Test run for a Mu-like packaging (e.g. Mu)

    Core processor number to use:
    -c CORE_NBR, --core=CORE_NBR
                        Number of core processor to use (Default: 1).



    Options for multi machine mode only

    Indicate that PhageTerm should run on several machines:
    --mm


    Options for step 1 of multi-machine mode (calculating reads coverage) on several machines

    Directory for coverage results:
    --dir_cov_mm=DIR_PATH/DIR_NAME
                        Directory where to put coverage results.
                        Note: it is up to the user to delete the files in this directory.

    Total number of cores to use
    -c CORE_NBR, --core=CORE_NBR
                        Total number used accross over all machines.

    Index of read chunk to process on current core
    --core_id=IDX
                A number between 0 and CORE_NBR-1

    Directory for checkpoint files:
    --dir_chk=DIR_PATH/DIR_NAME
                    Directory where phageTerm will put its ceckpoints.
                    Note: the directory must exist before launching phageTerm.
                    If the directory already contains a file, phageTerm will start from the results contained in this file.

    --chk_freq=FREQUENCY
                    The frequency in minutes at which checkpoints must be created.
                    Note: default value is 0 which means that no checkpoint is created.



    Options for step 2 of multi-machine mode (calculating per sequence statistics from reads coverage results) on several machines

    Directory for coverage results:
    --dir_cov_mm=DIR_PATH/DIR_NAME
                        Directory where to put coverage results.
                        Note: it is up to the user to delete the files in this directory.

    Directory for per sequence results
    --dir_seq_mm=DIR_PATH/DIR_NAME
                        Directory where to put the information if no match was found for one/several sequences.
                        Note: it is up to the user to delete the files in this directory.

    Directory for DR results
    --DR_path=DIR_PATH/DIR_NAME
                        Directory where to put the information necessary to step 3 (final report generation).
                        This information typically includes names of phage found and per sequence statistics.
                        Note: it is up to the user to delete the files in this directory.

    Sequence identifier
    --seq_id=IDX
            Index of the sequence to be processed by the current phageTerm process.
            Let N be the number of sequences given at the end of step 1.
            Then IDX is  number between 0 and N-1.

    Number of pieces
    --nb_pieces=NP
            Number of parts in which the reads were divided.
            Must be the same value as given via -c at step 1 (CORE_NBR).


    Options for step 3 of multi-machine mode (final report generation)

    Directory for DR results
    --DR_path=DIR_PATH/DIR_NAME
                        Directory where to read the information necessary to step 3 (final report generation).
                        This information typically includes names of phage found and per sequence statistics.
                        Note: it is up to the user to delete the files in this directory.

    Directory for per sequence results
    --dir_seq_mm=DIR_PATH/DIR_NAME
                        Directory where to get the information if no match was found for one/several sequences.
                        Note: it is up to the user to delete the files in this directory.




               
                        
# OUTPUT FILES
# ==========

	(i) Report (.pdf)
	
	(ii) Statistical table (.csv) 

	(iii) File containingg contains re-organized to stat at the predicted termini (.fasta)
	

# CONTACT
# =======

Julian Garneau <julian.garneau@usherbrooke.ca> 

Marc Monot <marc.monot@pasteur.fr>

David Bikard <david.bikard@pasteur.fr>

Véronique Legrand <vlegrand@pasteur.fr>
