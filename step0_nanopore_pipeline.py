"""
#### RUN AS ROOT W/ ENV (sudo -E) #### Maybe? I dunno. Running as user may work fine
step0_nanopore_pipeline.py
Marcus Viscardi     May 29, 2021

I have a ~prefix set of shell scripts that can go from raw ONT fast5 files
    all the way to called, mapped, & tail-called reads. This system works but
    is lacking in it's ability to handle multiple input formats and ease of
    use (I've had to reverse engineer these scripts every time I used them).

The goal here is to have a pipeline script similar to the pipelineWrappers we
    have for ribo-seq, including a input parsing and some more useful comments
    for myself (and others) in the future.
    
Side goal here is to try and update these scripts to use the newest base-
    callers from ONT, rather than the guppy_caller that was used during the
    design of nanopolish. This will also require making sure that nanopolish
    works with newer callers. I assume that it should as all it is really
    "reading" is the fast5 voltage trace, and using the basecalled reads to
    anchor where it's starting and finishing the polyA call.
    
ALSO, I wanna have this script check to see what parts of the pipeline have
    already been run, so I can avoid overwriting files unnecessarily. Even
    further, it would be nice to have this script be able to rerun specified
    steps as needed.
    
My recent nanopore run took about ~160hr of CPU time to get to 30%. That
    extrapolates out to >500 CPU hours to call a single (long) minIon run!!
    
New GPU based calling is WAY better. . . Why did I ever run with CPU?!
"""
# TODO: add regen tag whenever a upstream file is missing!
from os import path, listdir, mkdir
from argparse import ArgumentParser
from subprocess import Popen, CalledProcessError, PIPE
from typing import List
from glob import glob
from nanoporePipelineCommon import find_newest_matching_file, get_dt, minimap_bam_to_df

import pandas as pd
import numpy as np

pd.set_option("display.max_columns", None)


def live_cmd_call(command):
    with Popen(command, stdout=PIPE, shell=True,
               bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end="")
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)


def gene_names_to_gene_ids(tsv_path: str = "/data16/marcus/genomes/elegansRelease100"
                                           "/Caenorhabditis_elegans.WBcel235.100.gtf"
                                           ".tsv") -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t")[["gene_name", "gene_id"]].drop_duplicates(ignore_index=True)
    return df


#################################################################################
# Step0: Parse inputs
#################################################################################
def meshSetsAndArgs(skip_cli_dict: dict = None) -> dict:
    # First take inputs from the command line (this first so that the user can point
    #   to the settings file). Spit out a dictionary of called settings.
    def parseArgs() -> dict:
        parser = ArgumentParser(description="A pipeline to handle the outputs "
                                            "of nanopore runs!")
        # Required Arguments:
        parser.add_argument('settings', metavar='settings', type=str,
                            help="A .txt file with inputs all in arg|value format. "
                                 "Any of the below arguments will also function if put "
                                 "into this file.")
        # Arguments that can be included in the settings file
        parser.add_argument('--outputDir', metavar='outputDir', type=str, default=None,
                            help="Directory for the output of all resulting files.")
        parser.add_argument('--genomeDir', metavar='genomeDir', type=str, default=None,
                            help="Path to genome directory.")
        parser.add_argument('--dataDir', metavar='dataDir', type=str, default=None,
                            help="Path to sequencing data directory.")
        parser.add_argument('--threads', metavar='threads', type=int, default=None,
                            help="Number of threads to be used by nanopolish and minimap2. [20]")
        parser.add_argument('--guppyConfig', metavar='guppyConfig', type=str,
                            help="Configuration preset passed to the guppy_basecaller "
                                 "based on flowcell and kit used for run. Helpful "
                                 "table for picking a config @ https://denbi-nanopore-"
                                 "training-course.readthedocs.io/en/latest/basecalling/"
                                 "basecalling.html or just google 'guppy_basecalling'"
                                 "[rna_r9.4.1_70bps_hac.cfg]")
        parser.add_argument('--dropGeneWithHitsLessThan', metavar='dropGeneWithHitsLessThan',
                            type=int, default=None,
                            help="Number of threads to be used by guppy_basecaller, "
                                 "nanopolish and minimap2.")
        parser.add_argument('--altGenomeDirs', metavar='altGenomeDirs',
                            nargs='*', type=List[str], default=None,
                            help="Alternative genomes that can be used for filtering out reads "
                                 "that map to them")
        parser.add_argument('--stepsToRun', metavar='stepsToRun',
                            type=str, default=None,
                            help="Steps to run within the pipeline: (G)uppy basecalling, "
                                 "(M)inimap, (N)anopolish, (F)eature counts, (C)oncat files, "
                                 "merge with (P)andas, use f(L)air to ID transcripts "
                                 "(not default behavior), map pTRI nanopore (S)tandards, "
                                 "and/or random e(X)tra steps (plotting). [GMNFCPS]")
        parser.add_argument('--sampleID', metavar='sampleID',
                            type=int, default=None,
                            help="sampleID to pass to FLAIR [sample1]")
        parser.add_argument('--condition', metavar='condition',
                            type=int, default=None,
                            help="condition to pass to FLAIR [conditionA]")
        parser.add_argument('--minimapParam', metavar='minimapParam',
                            type=str, default=None,
                            help="Arguments to pass to minimap2. If used on the command line, "
                                 "be sure to surround in double-quotes! [\"-x splice -uf -k14\"]")
        # Flag Arguments
        parser.add_argument('-p', '--printArgs', action='store_true',
                            help="Boolean flag to show how arguments are overwritten/accepted")
        parser.add_argument('-n', '--nestedData', action='store_true',
                            help="Boolean flag that will account for data coming out of gridIONs "
                                 "as these will produce large nested fast5 dictionaries.")
        parser.add_argument('-r', '--regenerate', action='store_true',
                            help="Boolean flag to ignore previously produced files "
                                 "and generate all files anew")
        # parser.add_argument('-A', '--altGenomeProcessing', action='store_true',
        #                     help="Boolean flag to ")

        # Spit out namespace object from argParse
        args = parser.parse_args()
        # Quickly convert Namespace object to dictionary
        arg_dict = {arg: vars(args)[arg] for arg in vars(args)}

        if arg_dict['printArgs']:
            # Print arguments, specifies arguments which will not be passed in the arg_dict
            print("\nGiven Arguments (ArgParse):")
            for key, arg in arg_dict.items():
                if not arg:
                    print(f"\t{key} = {arg} -> (Will not be passed)")
                else:
                    print(f"\t{key} = {arg}")
        # Recreate dict without arguments that did not receive any input
        arg_dict = {k: v for k, v in arg_dict.items() if v is not None and v}
        return arg_dict

    # Second parse the settings file. This will also make a dictionary of all the
    #   called settings.
    def parseSettings(settings, printArgs=False, **other_kwargs) -> dict:
        """
        Will loop through and replace variables that are None
        """
        # first, parse the settings file to a dictionary called settingsDict
        settingsDict = {}
        with open(settings, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    line = line.strip()
                    if line != '':
                        line = line.split('|')
                        if len(line) == 2:
                            if line[0] == "altGenomeDirs":
                                settingsDict[line[0]] = line[1].split(",")
                            else:
                                settingsDict[line[0]] = line[1]
                        else:
                            print("\033[31;1m\nRemove pipes ('|') from settings "
                                  "file arguments (or rewrite parser)\n\033[0m")
                            exit()
        if printArgs:
            print(f"\nSettings Arguments (file: '{settings}')")
            for key, arg in settingsDict.items():
                if not arg:
                    print(f"\t{key} = {arg} -> (Will not be passed)")
                else:
                    print(f"\t{key} = {arg}")
        settingsDict = {k: v for k, v in settingsDict.items() if v is not None and v != ''}
        return settingsDict

    # Merge these dictionaries!
    absoluteDefDict = dict(printArgs=False,
                           nestedData=False,
                           regenerate=False,
                           altGenomeDirs=[],
                           threads=20,
                           stepsToRun="GMNFCPS",
                           sampleID="sample1",
                           condition="conditionA",
                           minimapParam="-x splice -uf -k14",
                           guppyConfig="rna_r9.4.1_70bps_hac.cfg")
    if skip_cli_dict:
        argDict = skip_cli_dict
    else:
        argDict = parseArgs()
    settingsDict = parseSettings(**argDict)
    finalArgDict = {}

    finalArgDict.update(absoluteDefDict)
    finalArgDict.update(settingsDict)
    finalArgDict.update(argDict)
    print("\033[1m\nPipeline Arguments:")

    # Absolute defaults overwritten by settings.txt then overwritten by CLI args
    for key, arg in finalArgDict.items():
        print(f"\t{key} = {arg}")
        if finalArgDict[key] == "True":
            finalArgDict[key] = True
        elif finalArgDict[key] == "False":
            finalArgDict[key] = False
        elif key == "altGenomeDirs":
            finalArgDict[key] = list(arg)
        else:
            try:
                finalArgDict[key] = int(arg)
            except ValueError or TypeError:
                finalArgDict[key] = str(arg)
    return finalArgDict


def buildOutputDirs(outputDir, stepsToRun, **kwargs) -> None:
    dirs_list = (("Z", outputDir),  # I am just going to use Z to mean always
                 ("G", "fastqs"),
                 ("M", "cat_files"),
                 ("N", "nanopolish"),  # TODO: split nanopolish from minimap2
                 ("F", "featureCounts"),
                 ("Z", "logs"),  # Z again
                 ("P", "merge_files"),
                 ("L", "flair"),
                 )
    print('\n')
    stepsToRun += "Z"
    for (steps_to_run_code, dir_to_make) in dirs_list:
        if steps_to_run_code in stepsToRun:
            if dir_to_make is outputDir:
                new_dir_path = outputDir
            else:
                new_dir_path = f"{outputDir}/{dir_to_make}"
            if not path.exists(new_dir_path):
                mkdir(new_dir_path)
            else:
                print(f"Directory @ {new_dir_path} already exists, skipping.")
    print("\n")


#################################################################################
# Step1: Trigger Docker Container which will do the actual guppy base-calling
#        (eventually I'll want to have this function outside of docker!)
#################################################################################
def guppy_basecall_w_docker(dataDir, outputDir, threads, guppyConfig, regenerate, **other_kwargs):
    prev_cat_fastq = path.exists(f"{outputDir}/cat_files/cat.fastq")
    if regenerate or not prev_cat_fastq:
        # TODO: do some math here instead of just splitting into 3
        callers = 3
        threads_per_caller = int(threads / callers)
        live_cmd_call(rf"""sudo docker run """
                      rf"""-v {dataDir}:/usr/src/prefix/data_dir """
                      rf"""-v {outputDir}:/usr/src/prefix/output_dir """
                      rf"""-it nanopore_empty """
                      rf"""guppy_basecaller --num_callers {callers} --cpu_threads_per_caller {threads_per_caller} """
                      rf"""-c {guppyConfig} -i /usr/src/prefix/data_dir/fast5 -s /usr/src/prefix/output_dir/fastqs """
                      rf"""2>&1 | tee {outputDir}/logs/{get_dt()}_guppy.log""")
        live_cmd_call(rf"cat {outputDir}/fastqs/*.fastq > {outputDir}/cat_files/cat.fastq")
    else:
        print(f"\n\nCalling already occurred. Based on file at:\n\t{outputDir}/cat_files/cat.fastq\n"
              f"Use the regenerate tag if you want to rerun calling.\n")


def guppy_basecall_w_gpu(dataDir, outputDir, threads, guppyConfig, regenerate, **other_kwargs):
    prev_cat_fastq = path.exists(f"{outputDir}/cat_files/cat.fastq")
    if regenerate or not prev_cat_fastq:
        guppy_log = f"{outputDir}/logs/{get_dt()}_guppy.log"
        call = rf"""guppy_basecaller -x "cuda:0" --num_callers 12 --gpu_runners_per_device 8 """ \
               rf"""-c {guppyConfig} -i {dataDir}/fast5 -s {outputDir}/fastqs """ \
               rf"""2>&1 | tee {guppy_log}"""
        print(f"Starting Guppy Basecalling with GPU @ {get_dt(for_print=True)}. . .\n"
              f"(The loading bar below will not change until the calling is completely finished!\n"
              f"if you want to watch progress run: \"watch less {guppy_log}\")\n\n"
              f"Running with call: {call}\n")
        # TODO: optimize the num_callers and gpu_runners_per_caller params!!
        live_cmd_call(call)
        live_cmd_call(rf"cat {outputDir}/fastqs/pass/*.fastq > {outputDir}/cat_files/cat.fastq")
        print(f"Finished Guppy Basecalling with GPU @ {get_dt(for_print=True)}. . .")
    else:
        print(f"\n\nCalling already occurred. Based on file at:\n\t{outputDir}/cat_files/cat.fastq\n"
              f"Use the regenerate tag if you want to rerun calling.\n")


#################################################################################
# Step2: Nanopolish index, map to genome with minimap2, samtools index/sort, &
#        Nanopolish calls of polyA tail lengths TODO: split minimap and nanopolish stuff
#################################################################################
def alternative_genome_filtering():
    pass


def minimap2_and_samtools(genomeDir, outputDir, threads, regenerate, minimapParam, **other_kwargs):
    minimap_flag = regenerate or not path.exists(f"{outputDir}/cat_files/cat.bam")
    if not minimap_flag:
        bam_length = path.getsize(f"{outputDir}/cat_files/cat.bam")
        if bam_length == 0:
            minimap_flag = True  # This is all to catch screwed up runs that have empty bam files!!!
    if minimap_flag:
        # TODO: Logging for the mapping call with minimap2
        genome_fa_file = glob(f"{genomeDir}/*allChrs.fa")
        if len(genome_fa_file) != 1:
            raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                      f"with one fa file that ends with 'allChrs.fa'")
        genome_bed_file = glob(f"{genomeDir}/*.bed")
        if len(genome_bed_file) != 1:
            raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                      f"with one bed file that ends with '.bed'")
        call = f"minimap2 -a {minimapParam} {genome_fa_file[0]} {outputDir}/cat_files/cat.fastq " \
               f"-t {threads} --junc-bed {genome_bed_file[0]} | samtools view -b - -o " \
               f"{outputDir}/cat_files/cat.bam"
        print(f"Starting minimap2 at {get_dt(for_print=True)}\nUsing call:\t{call}\n")
        live_cmd_call(call)
        print(f"\n\nFinished minimap2 at {get_dt(for_print=True)}\n")
    else:
        print(f"\n\nMiniMap2 already ran. Based on file at:"
              f"\n\t{outputDir}/cat_files/cat.bam\n"
              f"Use the regenerate tag if you want to rerun.\n")

    samtools_flag = regenerate or not path.exists(f"{outputDir}/cat_files/cat.sorted.bam")
    if samtools_flag:
        call = f"samtools sort -m 16G -T tmp -o {outputDir}/cat_files/cat.sorted.bam " \
               f"{outputDir}/cat_files/cat.bam && samtools index " \
               f"{outputDir}/cat_files/cat.sorted.bam"
        print(f"Starting samtools sort & index at {get_dt(for_print=True)}\nUsing call:\t{call}\n")
        live_cmd_call(call)
        print(f"\n\nFinished samtools sort & index at {get_dt(for_print=True)}")
    else:
        print(f"\n\nsamtools sort already ran. Based on file at:"
              f"\n\t{outputDir}/cat_files/cat.sorted.bam\n"
              f"Use the regenerate tag if you want to rerun.\n")


def nanopolish_index_and_polya(genomeDir, dataDir, outputDir, threads, regenerate, **other_kwargs):
    # First we have to index the fastq_files!
    nanopolish_index_flag = regenerate or not path.exists(f"{outputDir}/cat_files/cat.fastq.index.readdb")
    if nanopolish_index_flag:
        seq_sum_matches = [i for i in listdir(dataDir) if
                           path.isfile(path.join(dataDir, i)) and 'sequencing_summary' in i]
        if len(seq_sum_matches) != 1:
            raise IndexError(f"Too many/few matches for 'sequencing_summary' in sequencing directory!!\n"
                             f"{seq_sum_matches}")
        else:
            seq_sum_name = seq_sum_matches[0]
        call = f"nanopolish index --directory={dataDir}/fast5 " \
               f"--sequencing-summary={dataDir}/{seq_sum_name} " \
               f"{outputDir}/cat_files/cat.fastq"
        print(f"Starting nanopolish index at {get_dt(for_print=True)}\nUsing call:\t{call}\n"
              f"(There are limited outputs from this script, and it runs very slow)")
        live_cmd_call(call)
        print(f"Nanopolish index completed at {get_dt(for_print=True)}")
    else:
        print(f"\n\nNanopolish index already ran. Based on file at:"
              f"\n\t{outputDir}/cat_files/cat.fastq.index.readdb\n"
              f"Use the regenerate tag if you want to rerun.\n")

    # Now we are able to run the nanopolish polya script, this will throw an error if minimap2
    #   has not ran yet!!!!
    nanopolish_polya_flag = regenerate or not path.exists(f"{outputDir}/nanopolish/polya.tsv")
    if nanopolish_polya_flag:
        genome_fa_file = glob(f"{genomeDir}/*allChrs.fa")
        if len(genome_fa_file) != 1:
            raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                      f"with one fa files that ends with 'allChrs.fa'")
        call = f"nanopolish polya --threads={threads} --reads={outputDir}/cat_files/cat.fastq " \
               f"--bam={outputDir}/cat_files/cat.sorted.bam --genome={genome_fa_file[0]} " \
               f"> {outputDir}/nanopolish/polya.tsv"
        print(f"\nStarting nanopolish polyA at {get_dt(for_print=True)} (This takes a while)"
              f"\nUsing call:\t{call}\n")
        live_cmd_call(call)
        print(f"\n\nFinished nanopolish polyA at {get_dt(for_print=True)}")

        filter_call = f"head -n 1 {outputDir}/nanopolish/polya.tsv > {outputDir}/nanopolish/polya.passed.tsv; " \
                      f"grep PASS {outputDir}/nanopolish/polya.tsv >> {outputDir}/nanopolish/polya.passed.tsv"
        print(f"\nStarting filtering nanopolish polyA calls at {get_dt(for_print=True)}"
              f"\nUsing call:\t{filter_call}\n")
        live_cmd_call(filter_call)
        print(f"\n\nFinished filtering nanopolish polyA calls at {get_dt(for_print=True)}")
    else:
        print(f"\n\nNanopolish polyA already ran. Based on file at:"
              f"\n\t{outputDir}/nanopolish/polya.tsv\n"
              f"Use the regenerate tag if you want to rerun.\n")


#################################################################################
# Step4: Concatenate files (less important for single MinIon runs), and create
#        a single file that contains information from all the tools I am using
#################################################################################
def concat_files(outputDir, **other_kwargs):
    original_bam_file = f"{outputDir}/cat_files/cat.sorted.bam"
    bam_file_for_feature_counts = f"{outputDir}/cat_files/cat.sorted.mappedAndPrimary.bam"

    # Most of this is much less necessary as we are not getting the nested files that came out of Roach's gridION!
    calls = [f"samtools view  -b -F 0x904 {original_bam_file} > {bam_file_for_feature_counts}",
             # The above command will build a new bam file w/out reads w/ bit_flags:
             #    0x004, UNMAP           =   reads who's sequence didn't align to the genome
             #    0x100, SECONDARY       =   reads that are secondary alignments
             #    0x800, SUPPLEMENTARY   =   reads that are supplemental alignments
             f"samtools view {outputDir}/cat_files/cat.sorted.mappedAndPrimary.bam "
             f"> {outputDir}/cat_files/cat.sorted.mappedAndPrimary.sam",
             f"samtools view {outputDir}/cat_files/cat.sorted.bam "
             f"> {outputDir}/cat_files/cat.sorted.sam",
             ]
    print(f"Starting final cleanup at {get_dt(for_print=True)}\n")
    for num, call in enumerate(calls):
        print(f"\nStarting call ({num + 1} of {len(calls)}):\t{call}")
        live_cmd_call(call)
    print(f"\n\nFinished final cleanup at {get_dt(for_print=True)}")


#################################################################################
# Step3: featureCounts to identify the genes that reads map to, and how many hits
#        we have per gene
#################################################################################
def feature_counts(genomeDir, outputDir, regenerate, threads, **other_kwargs):
    feature_counts_flag = regenerate or \
                          not path.exists(f"{outputDir}/featureCounts/cat.sorted.mappedAndPrimary.bam.featureCounts")
    if feature_counts_flag:
        genome_gtf_file = glob(f"{genomeDir}/*.gtf")
        if len(genome_gtf_file) != 1:
            raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                      f"with one gtf files that ends with '.gtf'")
        call = f"featureCounts -L -T {threads} -R CORE -a {genome_gtf_file[0]} " \
               f"-o {outputDir}/featureCounts/{get_dt(for_file=True)} " \
               f"{outputDir}/cat_files/cat.sorted.mappedAndPrimary.bam " \
               f"2>&1 | tee -a {outputDir}/logs/{get_dt()}.featureCounts.log"
        print(f"Starting featureCounts at {get_dt(for_print=True)}\nUsing call:\t{call}\n")
        live_cmd_call(call)
        print(f"\n\nFinished featureCounts at {get_dt(for_print=True)}")

        filter_call = f"grep Assigned {outputDir}/featureCounts/cat.sorted.mappedAndPrimary.bam.featureCounts " \
                      f">> {outputDir}/featureCounts/cat.sorted.mappedAndPrimary.bam.Assigned.featureCounts"
        print(f"Filtering featureCounts calls at {get_dt(for_print=True)}\nUsing call:\t{filter_call}\n")
        live_cmd_call(filter_call)
        print(f"\n\nFinished filtering at {get_dt(for_print=True)}")
    else:
        print(f"\n\nPost calling processing already occurred. Based on file at:"
              f"\n\t{outputDir}/featureCounts/[...]\n"
              f"Use the regenerate tag if you want to rerun featureCounts.\n")


def merge_results(**other_kwargs):
    def create_merge_df(outputDir, keep_multimaps=False, print_info=False, **kwargs) -> pd.DataFrame:
        # First lets get the biggest one out of the way, importing the concatenated sam file:
        if not keep_multimaps:  # Pull the sam file that already had missed or secondaries dropped
            print(f"Starting to load SAM file from: {outputDir}/cat_files/cat.sorted.mappedAndPrimary.sam . . .",
                  end="")
            sam_df = pd.read_csv(f"{outputDir}/cat_files/cat.sorted.mappedAndPrimary.sam",
                                 sep="\t", names=range(22), low_memory=False, index_col=False)
        else:  # Otherwise: load the bam file that did not have those other reads dropped
            # Loading the bam with my function is slightly slower, but at least hides
            #   a good amount of the complicated bits!
            print(f"Starting to load BAM file from: {outputDir}/cat_files/cat.sorted.bam . . .",
                  end="")
            sam_df = minimap_bam_to_df(f"{outputDir}/cat_files/cat.sorted.bam",
                                       name_columns=False,
                                       drop_secondaries_and_unmapped=False).df
        print(f" Done!")
        # And lets rename columns!
        sam_header_names = ["read_id",  # string
                            "bit_flag",  # uint16
                            "chr_id",  # category
                            "chr_pos",  # uint32
                            "mapq",  # uint8
                            "cigar",  # string
                            "r_next",
                            "p_next",
                            "len",
                            "sequence",  # string
                            "phred_qual",  # string
                            ]
        extra_columns = ["num_mismatches",  # Keep, string
                         "best_dp_score",
                         "dp_score",
                         "num_ambiguous_bases",
                         "transcript_strand",  # Keep, string
                         "type_of_alignment",  # Keep, category (by the end)
                         "num_minimizes",
                         "chain_score",  # This column onwards is inconsistent,
                         "chain_score_top_secondary",  # b/c these are optional flags!
                         "gap_compressed_divergence",
                         "len_of_query_w_repeats"]
        sam_header_names += extra_columns
        sam_df = sam_df.rename(columns=dict(enumerate(sam_header_names)))

        for minimap_flag_column in ["type_of_alignment", "transcript_strand", "num_mismatches"]:
            sam_df[minimap_flag_column] = sam_df[minimap_flag_column].str.split(":").str[-1]

        # Make a list of columns to drop:
        extra_columns_to_drop = extra_columns + ["r_next", "p_next", "len"]
        # Remove the columns I want to keep from this "drop list"
        for col_to_keep in ["type_of_alignment", "transcript_strand", "num_mismatches"]:
            extra_columns_to_drop.remove(col_to_keep)
        # Drop the unsaved columns!
        sam_df = sam_df.drop(extra_columns_to_drop, axis=1)

        # Set dataframe column datatypes!
        # datatypes to use:
        o = "object"  # TODO: this could eventually be pd.StringDtype
        c = "category"
        ui8 = "uint8"
        ui16 = "uint16"
        ui32 = "uint32"

        df_dtypes = {"read_id": o,  # string
                     "bit_flag": ui16,  # uint16
                     "chr_id": c,  # category
                     "chr_pos": ui32,  # uint32
                     "mapq": ui8,  # uint8
                     "cigar": o,  # string
                     "sequence": o,  # string
                     "phred_qual": o,  # string
                     "num_mismatches": ui32,  # uint32, after parsing
                     "transcript_strand": c,  # category, after parsing. I don't really know what this column is...
                     #                              ALL of the values here are "+"?
                     "type_of_alignment": c,  # category, after parsing
                     }
        sam_df = sam_df.astype(df_dtypes)

        # Pull the 16 bit flag to get strand information (important for merge w/ featC later)
        sam_df["strand"] = (sam_df.bit_flag & 16).replace(to_replace={16: "-", 0: "+"})
        sam_df = sam_df.astype({"strand": c})

        if keep_multimaps:
            # Identify and drop reads that have the 4 bit flag: indicating they didn't map!
            sam_df = sam_df[(sam_df.bit_flag & 4) != 4]
        else:
            sam_df = sam_df[~sam_df.duplicated(subset="read_id", keep=False)]
        # Next lets pull in the featureCounts results (there is a population of duplicates here)
        featc_df = pd.read_csv(f"{outputDir}/featureCounts/cat.sorted.bam.Assigned.featureCounts",
                               sep="\t", names=["read_id", "qc_tag_featc", "qc_pass", "gene_id"])
        # Load up nanopolish polyA results (also a population of duplicates here!!)
        polya_df = pd.read_csv(f"{outputDir}/nanopolish/polya.passed.tsv", sep="\t")
        polya_df = polya_df.rename(columns={"readname": "read_id",
                                            "qc_tag": "qc_tag_polya"})  # b/c featC also has a qc_tag!
        names_df = gene_names_to_gene_ids()
        featc_df = featc_df.merge(names_df, on="gene_id")
        if print_info:
            print("#" * 100)
            print(f"\n\nSAM Dataframe info:")
            print(sam_df.info())
            print(f"\n\nfeatureCounts Dataframe info:")
            print(featc_df.info())
            print(f"\n\nPolyA Dataframe info:")
            print(polya_df.info())
        # LETS SMOOSH THEM ALL TOGETHER!!!
        # TODO: This is severely broken in terms of the featureCounts merge:
        #       B/c the featureCounts output only retains read_id, it doesn't
        #       have enough information to merge uniquely for reads that map
        #       more than once!! This means that read A that maps to gene X and Y
        #       is eventually producing 4 lines of data....
        # TODO: Revisiting on 10/26/2021: This is still broken. B/c the multiple
        #       hits (meaning multiple identical read_ids) in the bam/sam file are passed to
        #       feature counts, it propagates any multi-mappers. Currently I avoid this by
        #       dropping ANY read that hits more than once, meaning that propagation of
        #       multiple reads is avoided entirely. But it seems like I have a lot of good
        #       primary maps w/ trash secondaries! It would be really nice to retain those
        #       reads. . . Should I just switch to Josh's assignment method?
        #           OR: I could try to rename reads in the bam file w/ their map location,
        #           as this could help to uniquely identify primaries, and that info would
        #           get passed though featureCounts!
        # 8/24/2021: Sorta fixed this by just dropping all duplicate reads! (line 484)
        #            Back to square one... lol
        sam_featc_df = sam_df.merge(featc_df, how="left", on=["read_id"])
        merge_df = sam_featc_df.merge(polya_df, how="inner", left_on=["read_id"],
                                      right_on=["read_id"])
        merge_df.drop(columns=["contig", "position", "r_next", "p_next", "len"], inplace=True)
        merge_df = merge_df.drop_duplicates()
        merge_df = merge_df[merge_df["sequence"] != "*"]
        merge_df = merge_df[merge_df["mapq"] != 0]
        if print_info:
            print("\n\n")
            print("#" * 100)
            print(f"\n\nMerged Dataframe info:")
            print(merge_df.info())
        with open(f"{outputDir}/merge_files/{get_dt(for_file=True)}_mergedOnReads.tsv", "w") as merge_out_f:
            merge_df.to_csv(merge_out_f, sep="\t", index=False)
        return merge_df

    def compress_on_genes(merged_df, outputDir, dropGeneWithHitsLessThan=None,
                          output_to_file=True, **kwargs) -> pd.DataFrame:
        # This creates a pandas "groupby" object that can be used to extract info compressed on gene_ids
        print("\nMerging information from Minimap2, featureCounts and Nanopolish-PolyA:")
        merged_df["read_length"] = merged_df["sequence"].str.len()
        grouped_genes = merged_df.groupby("gene_id")

        gene_df = grouped_genes["read_id"].apply(len).to_frame(name="read_hits")
        gene_df["read_ids"] = grouped_genes["read_id"].apply(list).to_frame(name="read_ids")

        merged_df["read_length"] = merged_df["sequence"].str.len()
        gene_df["read_len_mean"] = grouped_genes["read_length"].apply(np.mean).to_frame(name="read_len_mean")
        gene_df["read_len_std"] = grouped_genes["read_length"].apply(np.std).to_frame(name="read_len_std")
        gene_df["read_lengths"] = grouped_genes["read_length"].apply(list).to_frame(name="read_lengths")

        gene_df["polya_lengths"] = grouped_genes["polya_length"].apply(list).to_frame(name="polya_lengths")
        gene_df["polya_mean"] = grouped_genes["polya_length"].apply(np.mean).to_frame(name="polya_mean")
        gene_df["polya_stdev"] = grouped_genes["polya_length"].apply(np.std).to_frame(name="polya_stdev")

        gene_df["genomic_starts"] = grouped_genes["chr_pos"].apply(list).to_frame(name="genomic_starts")
        gene_df["cigars"] = grouped_genes["cigar"].apply(list).to_frame(name="cigars")

        if dropGeneWithHitsLessThan:
            print(f"Dropping any genes with less than {dropGeneWithHitsLessThan} read hits")
            gene_df = gene_df[gene_df["read_hits"] >= dropGeneWithHitsLessThan]
        gene_df.sort_values("read_hits", ascending=False, inplace=True)
        print(f"Mean PolyA Length: {gene_df['polya_mean'].mean():.3f}")
        filename = f"{get_dt(for_file=True)}_compressedOnGenes"
        if output_to_file:
            gene_df.to_parquet(f"{outputDir}/merge_files/{filename}.parquet")
            gene_df.to_csv(f"{outputDir}/merge_files/{filename}.tsv", sep="\t")
            gene_df.drop(['read_ids',
                          'polya_lengths',
                          'read_lengths',
                          'genomic_starts',
                          'cigars'],
                         axis=1).to_csv(f"{outputDir}/merge_files/{filename}_simple.tsv", sep="\t")
        return gene_df

    print(f"Starting to merge all data at {get_dt(for_print=True)}\n")
    merge_df = create_merge_df(**other_kwargs)
    print(f"\n\nFinished merging all data at {get_dt(for_print=True)}")
    print(f"Starting to compress data on genes at {get_dt(for_print=True)}\n")
    gene_df = compress_on_genes(merge_df, **other_kwargs)
    print(f"\n\nFinished compressing data on genes at {get_dt(for_print=True)}")
    return merge_df, gene_df


def flair(outputDir, **other_kwargs):
    def run_and_load_flair(outputDir, genomeDir, threads, sampleID, condition, regenerate, **kwargs) -> pd.DataFrame:
        flair_map_path = f"{outputDir}/flair/flair.quantify.{sampleID}.{condition}.isoform.read.map.txt"
        if regenerate or not path.exists(flair_map_path):
            manifest_path = f"{outputDir}/flair/read_manifest.tsv"
            with open(manifest_path, "w") as f:
                f.write(f"{sampleID}\t{condition}\tbatch1\t{outputDir}/cat_files/cat.fastq")
            call = f"cd {outputDir}/flair ; python3 /data16/marcus/scripts/brooksLabUCSC_flair/flair.py " \
                   f"quantify -t {threads} --generate_map -r {manifest_path} " \
                   f"-i {genomeDir}/Caenorhabditis_elegans.WBcel235.cdna.all.fa ; cd -"
            live_cmd_call(call)
        fl_df = pd.read_csv(flair_map_path, sep="\t", header=None)
        fl_df[1] = fl_df[1].str.split(",")
        fl_df = fl_df.explode(1).reset_index(drop=True)
        fl_df.rename(columns={0: "transcript_id",
                              1: "read_id"},
                     inplace=True)
        print(f"Loaded map file from FLAIR. . .\n")
        return fl_df

    def merge_some_more(flair_df: pd.DataFrame, outputDir, dropGeneWithHitsLessThan: int = None,
                        output_to_file=True) -> pd.DataFrame:
        merge_tsv_path = f"{outputDir}/merge_files/{get_dt(for_file=True)}_mergedOnReads.tsv"
        if not path.exists(merge_tsv_path):
            older_merge_tsv = find_newest_matching_file(f"{outputDir}/merge_files/*_mergedOnReads.tsv")
            print(f"Could not find file @ {merge_tsv_path}\n"
                  f"Going to use the older file @ {older_merge_tsv}")
            merge_tsv_path = older_merge_tsv
        merge_df = pd.read_csv(merge_tsv_path, sep="\t")
        print(f"\nLoaded my merge file. . .\n")
        print(merge_df.info())
        super_df = merge_df.merge(flair_df, how="inner", on="read_id")
        super_df.to_csv(f"{outputDir}/merge_files/{get_dt(for_file=True)}_mergedWithTranscripts.tsv",
                        sep="\t", index=False)

        super_df["read_length"] = super_df["sequence"].str.len()
        grouped_transcripts = super_df.groupby("transcript_id")

        transcript_df = grouped_transcripts["transcript_id"].apply(len).to_frame(name="transcript_hits")

        transcript_df["read_len_mean"] = grouped_transcripts["read_length"].apply(np.mean).to_frame(
            name="read_len_mean")
        transcript_df["read_len_std"] = grouped_transcripts["read_length"].apply(np.std).to_frame(name="read_len_std")
        transcript_df["read_lengths"] = grouped_transcripts["read_length"].apply(list).to_frame(name="read_lengths")

        transcript_df["polya_lengths"] = grouped_transcripts["polya_length"]. \
            apply(list).to_frame(name="polya_lengths")
        transcript_df["polya_mean"] = grouped_transcripts["polya_length"].apply(np.mean).to_frame(name="polya_mean")
        transcript_df["polya_stdev"] = grouped_transcripts["polya_length"].apply(np.std).to_frame(name="polya_stdev")

        transcript_df["gene_ids"] = grouped_transcripts["gene_id"].apply(list).to_frame(name="gene_ids")
        transcript_df["gene_id"] = transcript_df["gene_ids"].apply(lambda x: pd.Series(x).dropna().mode()[0:1])
        transcript_df.drop(["gene_ids"], axis=1, inplace=True)

        transcript_df["genomic_starts"] = grouped_transcripts["chr_pos"]. \
            apply(list).to_frame(name="genomic_starts")
        transcript_df["cigars"] = grouped_transcripts["cigar"]. \
            apply(list).to_frame(name="cigars")

        if dropGeneWithHitsLessThan:
            print(f"Dropping any genes with less than {dropGeneWithHitsLessThan} transcript hits")
            transcript_df = transcript_df[transcript_df["transcript_hits"] >= dropGeneWithHitsLessThan]
        transcript_df.sort_values("transcript_hits", ascending=False, inplace=True)
        print(f"Mean PolyA Length: {transcript_df['polya_mean'].mean():.3f}")
        filename = f"{get_dt(for_file=True)}_compressedOnTranscripts"
        if output_to_file:
            with open(f"{outputDir}/merge_files/{filename}.tsv", "w") as out_f:
                transcript_df.to_csv(out_f, sep="\t")
            with open(f"{outputDir}/merge_files/{filename}_simple.tsv", "w") as out_f:
                transcript_df.drop(["genomic_starts",
                                    "cigars",
                                    "polya_lengths",
                                    "read_lengths"], axis=1).to_csv(out_f, sep="\t")
        return transcript_df

    flair_df = run_and_load_flair(outputDir, **other_kwargs)
    final_df = merge_some_more(flair_df, outputDir,
                               output_to_file=True)


def extra_steps(outputDir, df=None, **other_kwargs):
    import plotly.express as px
    if not isinstance(df, pd.DataFrame):
        path = find_newest_matching_file(f"{outputDir}/merge_files/*_mergedOnReads.tsv")
        print(f"Dataframe not passed, loading file from: {path}")
        df = pd.read_csv(path, sep="\t")
        print("Loaded.")
    sample_df = df.sample(n=10000)
    sample_df["dup"] = sample_df.duplicated(subset="read_id", keep=False).map({True: 1, False: 0})
    sample_df["multi"] = (sample_df.bit_flag & 256) == 256
    read_ids_of_multis = sample_df[sample_df["multi"]]["read_id"]
    sample_df = sample_df[~sample_df["read_id"].isin(read_ids_of_multis.to_list())]
    # sample_df = sample_df[sample_df.dup == 1]
    sample_df["len"] = sample_df["sequence"].apply(len)
    # sample_df["bit_flag"] = sample_df["bit_flag"].apply(str)
    # fig = px.box(sample_df, y="mapq", x="bit_flag", notched=True, points="all")
    fig = px.parallel_coordinates(sample_df[["dup", "bit_flag", "mapq", "len"]],
                                  )
    fig.show()
    print(sample_df.info())


def map_standards(outputDir, df: pd.DataFrame = None, **other_kwargs):
    from standardsAlignment.standardsAssignmentWithMinimap2 import align_standards, plot_value_counts
    if not isinstance(df, pd.DataFrame):
        merge_dir = f"{outputDir}/merge_files"
        try:
            merge_on_reads_path = find_newest_matching_file(f"{merge_dir}/*mergedOnReads.parquet")
            df = pd.read_parquet(merge_on_reads_path)
        except ValueError:
            print(f"Could not find a mergedOnReads parquet file in directory:\n\t{merge_dir}")
            merge_on_reads_path = find_newest_matching_file(f"{merge_dir}/*mergedOnReads.tsv")
            df = pd.read_csv(merge_on_reads_path, sep="\t", low_memory=False)
    df = df.merge(align_standards(compressed_df=df, keep_read_id=True, **other_kwargs), on="read_id")
    df.to_parquet(f"{outputDir}/merge_files/{get_dt(for_file=True)}_mergedOnReads.plusStandards.parquet")


def main(stepsToRun, **kwargs) -> (pd.DataFrame, pd.DataFrame) or None:
    return_value = None
    buildOutputDirs(**args)

    steps_dict = {"G": [guppy_basecall_w_gpu, "Guppy Basecalling"],
                  "A": [alternative_genome_filtering, "Filtering Alt. Genomes (not implemented)"],
                  # TODO: Add alternative genome filter mapping function here^^^
                  "M": [minimap2_and_samtools, "Minimap2 and SamTools"],
                  "N": [nanopolish_index_and_polya, "Nanopolish Index and polyA Calling"],
                  "F": [feature_counts, "FeatureCounts"],
                  "C": [concat_files, "Concatenate Files"],
                  "P": [merge_results, "Merging Results w/ Pandas"],
                  "S": [map_standards, "Mapping Standards (still experimental!!)"],
                  "L": [flair, "Calling Transcripts w/ Flair"],
                  "X": [extra_steps, "Running random eXtra steps"]
                  }
    for code, (step_func, step_name) in steps_dict.items():
        if code in stepsToRun:
            step_print = f"Starting step: \"{step_name}\" . . ."
            print("\n\n", step_print, f"\n{'#' * len(step_print)}\n", sep="")
            if code == "P":  # The pandas merge will both: save a file and return a tuple of dataframes
                return_value = step_func(**args)
            elif code == "X" or code == "S":
                if return_value:  # The extra steps function will accept a df, or load from the disk
                    step_func(df=return_value[0], **args)
                else:
                    step_func(**args)
            else:  # The rest of the functions work by loading and saving to the disk
                step_func(**args)
        else:
            step_print = f"Skipping step: \"{step_name}\" . . ."
            print("\n", step_print, f"\n{'#' * len(step_print)}\n", sep="")
    return return_value


if __name__ == '__main__':
    args = meshSetsAndArgs()
    main_output = main(**args)
