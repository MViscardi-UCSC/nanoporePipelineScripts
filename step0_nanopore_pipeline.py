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

Late 2021:
    Added features to incorporate TERA-Seq into this main pipeline!
"""
# TODO: add regen tag whenever a upstream file is missing!
from os import path, listdir, mkdir, environ
from argparse import ArgumentParser
from subprocess import Popen, CalledProcessError, PIPE
from typing import List
from glob import glob
import warnings

from tqdm import tqdm

from nanoporePipelineCommon import find_newest_matching_file, get_dt, \
    gene_names_to_gene_ids, SamOrBamFile, FastqFile, adjust_5_ends, assign_with_josh_method

import numpy as np
import pandas as pd

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)

# This is also in my .bashrc but that doesn't always seem to trigger,
#   ie. if I run the script from inside pycharm.
# If the HDF5 path isn't specified nanopolish freaks out, this solution is based on:
#   https://stackoverflow.com/questions/5971312/how-to-set-environment-variables-in-python
environ['HDF5_PLUGIN_PATH'] = '/usr/local/hdf5/lib/plugin'

pd.set_option("display.max_columns", None)


def live_cmd_call(command):
    with Popen(command, stdout=PIPE, shell=True,
               bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end="")
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)


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
                            help="Minimum number of reads per gene to have in the "
                                 "compressedOnGenes outputs.")
        parser.add_argument('--altGenomeDirs', metavar='altGenomeDirs',
                            nargs='*', type=List[str], default=None,
                            help="Alternative genomes that can be used for filtering out reads "
                                 "that map to them")
        parser.add_argument('--stepsToRun', metavar='stepsToRun',
                            type=str, default=None,
                            help="Steps to run within the pipeline: (G)uppy basecalling, "
                                 "(T)rim TERA-Seq adapters, "
                                 "(M)inimap, (N)anopolish, (F)eature counts, (C)oncat files, "
                                 "merge with (P)andas, use f(L)air to ID transcripts, "
                                 "map pTRI nanopore (S)tandards, and/or random e(X)tra "
                                 "steps (plotting). [GMNCFPS]")
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
        parser.add_argument('--tera3adapter', metavar='tera3adapter',
                            type=int, default=None,
                            help="Adapter to be trimmed for TERA3 [None]")
        parser.add_argument('--tera5adapter', metavar='tera5adapter',
                            type=int, default=None,
                            help="Adapter to be trimmed for 5TERA [None]")
        parser.add_argument('--extraGuppyOptions', metavar='extraGuppyOptions',
                            type=int, default=None,
                            help="String flags/options to be added to the gupy_basecaller call [None]")
        # Flag Arguments
        parser.add_argument('-p', '--printArgs', action='store_true',
                            help="Boolean flag to show how arguments are overwritten/accepted")
        parser.add_argument('-n', '--nestedData', action='store_true',
                            help="Boolean flag that will account for data coming out of gridIONs "
                                 "as these will produce large nested fast5 dictionaries.")
        parser.add_argument('-r', '--regenerate', action='store_true',
                            help="Boolean flag to ignore previously produced files "
                                 "and generate all files anew")
        parser.add_argument('-j', '--callWithJoshMethod', action='store_true',
                            help="Boolean flag to use Josh's read assignment method, rather"
                                 "than FeatureCounts")
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
                           stepsToRun="GMNCFPS",
                           sampleID="sample1",
                           condition="conditionA",
                           minimapParam="-x splice -uf -k14",
                           guppyConfig="rna_r9.4.1_70bps_hac.cfg",
                           tera3adapter=False,
                           tera5adapter=False,
                           extraGuppyOptions=False,
                           callWithJoshMethod=False)
    if skip_cli_dict:
        argDict = skip_cli_dict
    else:
        argDict = parseArgs()
    settingsDict = parseSettings(**argDict)
    finalArgDict = {}

    # Start with the absolute defaults
    finalArgDict.update(absoluteDefDict)
    # Overwrite any defaults with things found in the settings file
    finalArgDict.update(settingsDict)
    # Overwrite any defaults or settings file calls with the command line (or passed dict)
    finalArgDict.update(argDict)
    print("\033[1m\nPipeline Arguments:")

    # Absolute defaults overwritten by settings.txt then overwritten by CLI args
    for key, arg in finalArgDict.items():
        print(f"\t{key} = {arg}\t->", end="\t")
        if finalArgDict[key] == "True" or finalArgDict[key] is True:
            finalArgDict[key] = True
        elif finalArgDict[key] == "False" or finalArgDict[key] is False:
            finalArgDict[key] = False
        elif key == "altGenomeDirs":
            finalArgDict[key] = list(arg)
        else:
            try:
                finalArgDict[key] = int(arg)
            except ValueError or TypeError:
                finalArgDict[key] = str(arg)
        print(f"{key} = {arg}")
    return finalArgDict


def buildOutputDirs(stepsToRun, **kwargs) -> None:
    outputDir = kwargs["outputDir"]
    dirs_list = (("Z", outputDir),  # I am just going to use Z to mean always
                 ("G", "fastqs"),
                 ("M", "cat_files"),
                 ("N", "nanopolish"),
                 ("F", "featureCounts"),
                 ("Z", "logs"),
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
def guppy_basecall_w_gpu(dataDir, outputDir, threads, guppyConfig, regenerate, extraGuppyOptions, nestedData, **other_kwargs):
    # TODO: I may need to change the --trim_strategy for TERA3!! add an param here for tera3adapter,
    #       if that param is not None, than I'll probably want to add the '--trim_strategy none'!!
    prev_cat_fastq = path.exists(f"{outputDir}/cat_files/cat.fastq")
    if regenerate or not prev_cat_fastq:
        guppy_log = f"{outputDir}/logs/{get_dt()}_guppy.log"
        if isinstance(extraGuppyOptions, str):
            extra_guppy_options = extraGuppyOptions + " "
        else:
            extra_guppy_options = ""
        if nestedData:
            extra_guppy_options += "-r "
        call = rf"""guppy_basecaller -x "cuda:0" {extra_guppy_options}--num_callers 12 --gpu_runners_per_device 8 """ \
               rf"""-c {guppyConfig} -i {dataDir}/fast5 -s {outputDir}/fastqs """ \
               rf"""2>&1 | tee {guppy_log}"""
        print(f"Starting Guppy Basecalling with GPU @ {get_dt(for_print=True)}. . .\n"
              f"(The loading bar below will not change until the calling is completely finished!\n"
              f"if you want to watch progress run: \"watch less {guppy_log}\")\n\n"
              f"Running with call: {call}\n")
        # TODO: optimize the num_callers and gpu_runners_per_caller params!!
        live_cmd_call(call)
        print(f"Finished Guppy Basecalling with GPU @ {get_dt(for_print=True)}. . .\n")
        print(f"Starting to concatenate fastq files @ {get_dt(for_print=True)}. . .\n")
        live_cmd_call(rf"cat {outputDir}/fastqs/pass/*.fastq > {outputDir}/cat_files/cat.fastq")
        print(f"Finished to concatenating fastq files @ {get_dt(for_print=True)}. . .")
    else:
        print(f"\n\nCalling already occurred. Based on file at:\n\t{outputDir}/cat_files/cat.fastq\n"
              f"Use the regenerate tag if you want to rerun calling.\n")


#################################################################################
# Step2: Nanopolish index, map to genome with minimap2, samtools index/sort, &
#        Nanopolish calls of polyA tail lengths
#################################################################################
def alternative_genome_filtering(altGenomeDirs, outputDir, threads, minimapParam, regenerate, **other_kwargs):
    import filecmp
    altmap_flag = regenerate or not path.exists(f"{outputDir}/cat_files/cat.altGenome.sam")
    if not altmap_flag:
        bam_length = path.getsize(f"{outputDir}/cat_files/cat.altGenome.sam")
        if bam_length == 0:
            altmap_flag = True  # This is all to catch screwed up runs that have empty bam files!!!
    if altmap_flag:
        for alt_genome in altGenomeDirs:
            alt_genome_fa_file = glob(f"{alt_genome}/*.fsa")
            if len(alt_genome_fa_file) != 1:
                raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                          f"with one fsa file that ends with '.fsa'")
            alt_genome_bed_file = glob(f"{alt_genome}/*.bed")
            if len(alt_genome_bed_file) != 1:
                raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                          f"with one bed file that ends with '.bed'")
            call = f"minimap2 -a {minimapParam} {alt_genome_fa_file[0]} {outputDir}/cat_files/cat.fastq " \
                   f"-t {threads} --junc-bed {alt_genome_bed_file[0]} " \
                   f" --sam-hit-only " \
                   f"> {outputDir}/cat_files/cat.altGenome.sam"
            print(f"Starting alt_genome minimap2 at {get_dt(for_print=True)}\nUsing call:\t{call}\n")
            live_cmd_call(call)
            print(f"\n\nFinished alt_genome minimap2 at {get_dt(for_print=True)}\n")
    else:
        print(f"\n\nAlternative genome MiniMap2 already ran. Based on file at:"
              f"\n\t{outputDir}/cat_files/cat.altGenome.bam\n"
              f"Use the regenerate tag if you want to rerun.\n")

    fastq_backup_flag = regenerate or not path.exists(f"{outputDir}/cat_files/cat.pre_altGenome.fastq")
    if fastq_backup_flag:
        # First backup the fastq file that the real minimap2 call will need:
        call = f"cp {outputDir}/cat_files/cat.fastq {outputDir}/cat_files/cat.pre_altGenome.fastq"
        print(f"Starting fastq backup at {get_dt(for_print=True)}\nUsing call:\t{call}\n")
        live_cmd_call(call)
        print(f"\n\nFinished fastq backup at {get_dt(for_print=True)}")
    else:
        print(f"\n\nBacking up of cat.fastq file already happened. Based on file at:"
              f"\n\t{outputDir}/cat_files/cat.pre_altGenome.fastq\n"
              f"Use the regenerate tag if you want to rerun.\n")

    alt_filtering_flag = regenerate or filecmp.cmp(f"{outputDir}/cat_files/cat.pre_altGenome.fastq",
                                                   f"{outputDir}/cat_files/cat.fastq")
    if alt_filtering_flag:
        # Then we'll load the alt genome mapped reads from the sam file to a pd.Dataframe:
        #   The below call might break if minimap2 passed headers!
        print(f"Starting to load alt genome called sam file from: {outputDir}/cat_files/cat.altGenome.sam . . .")
        header_lines = 0
        with open(f"{outputDir}/cat_files/cat.altGenome.sam", 'r') as sam_file_quick:
            for line in sam_file_quick:
                if line.startswith('@'):
                    header_lines += 1
                else:
                    break
        alt_mapped_read_df = pd.read_csv(f"{outputDir}/cat_files/cat.altGenome.sam", sep="\t",
                                         usecols=[0], index_col=False, header=None,
                                         low_memory=False, skiprows=header_lines)
        alt_mapped_read_df.rename(columns={0: "read_id"}, inplace=True)
        print(f"Finished loading alt genome called sam file from: {outputDir}/cat_files/cat.altGenome.sam")
        # Next we'll build a dataframe of the fastq file!
        print(f"Starting to load fastq file from: {outputDir}/cat_files/cat.pre_altGenome.fastq . . .")
        fastq_file = FastqFile(f"{outputDir}/cat_files/cat.pre_altGenome.fastq")
        print(f"Finished loading fastq file from: {outputDir}/cat_files/cat.pre_altGenome.fastq . . .")
        fastq_file.filter_against(alt_mapped_read_df)
        fastq_file.save_to_fastq(f"{outputDir}/cat_files/cat.fastq")
    else:
        print(f"\n\nAlternative genome fastq filtering already ran. Based on the result that files:"
              f"\n\t{outputDir}/cat_files/cat.pre_altGenome.fastq and cat.fastq"
              f"\n\tAre 'the same'\n"
              f"Use the regenerate tag if you want to rerun.\n")
    # TODO: Add fastq file filtering based off what ends up in the cat.altGenome_hits.sam
    #       First copy the current cat.fastq to "cat.preAltGenome.fastq".
    #       Then make a new "cat.fastq" with only the reads that didn't map to the altGenome


def trim_tera_adapters(outputDir, threads, regenerate, tera3adapter, tera5adapter, **other_kwargs):
    # TODO: There is an error here that is leading to rerunning on already cutadapt'ed fastqs!!
    #  maybe solved w/ dropping the regen flag......

    # First step is to backup the original fastq file:
    fastq_backed_up = path.exists(f"{outputDir}/cat_files/cat.untrimmed.fastq")
    if not fastq_backed_up:
        # First backup the fastq file that the real minimap2 call will need:
        call = f"cp {outputDir}/cat_files/cat.fastq {outputDir}/cat_files/cat.untrimmed.fastq"
        print(f"Starting fastq backup at {get_dt(for_print=True)}\nUsing call:\t{call}\n")
        live_cmd_call(call)
        print(f"\n\nFinished fastq backup at {get_dt(for_print=True)}")
    else:
        print(f"\n\nBacking up of cat.fastq file already happened. Based on file at:"
              f"\n\t{outputDir}/cat_files/cat.untrimmed.fastq\n"
              f"Use the regenerate tag if you want to rerun.\n")

    # Then, we will want to check that adapter trimming hasn't already happened.
    #   This isn't quite as simple as file checking, but the easiest way I can
    #   think of will be to parse the first line of the fastq and check if
    #   'adapter' is in there:
    with open(f'{outputDir}/cat_files/cat.fastq', 'r') as fastq_file:
        first_line = fastq_file.readline()
        cutadapt_was_run = 'TERAADAPTER' in first_line
    if not cutadapt_was_run:
        cutadapt_call = None
        if isinstance(tera5adapter, str) and isinstance(tera3adapter, str):
            cutadapt_call = f"cutadapt --action=trim -j {threads} " \
                            f"-g TERA5={'X' + tera5adapter} --overlap 31 --error-rate 0.29 " \
                            "--rename '{id} TERAADAPTER5={adapter_name} {comment}' " \
                            f"{outputDir}/cat_files/cat.untrimmed.fastq | " \
                            f"cutadapt --action=trim -j {threads} " \
                            f"-g TERA3={tera3adapter} --overlap 16 --error-rate 0.18 " \
                            "--rename '{id} TERAADAPTER3={adapter_name} {comment}' " \
                            f"- > {outputDir}/cat_files/cat.fastq"
        elif isinstance(tera5adapter, str):
            cutadapt_call = f"cutadapt --action=trim -j {threads} " \
                            f"-g TERA5={'X' + tera5adapter} --overlap 31 --error-rate 0.29 " \
                            "--rename '{id} TERAADAPTER5={adapter_name} {comment}' " \
                            f"{outputDir}/cat_files/cat.untrimmed.fastq > {outputDir}/cat_files/cat.fastq"
        elif isinstance(tera3adapter, str):
            cutadapt_call = f"cutadapt --action=trim -j {threads} " \
                            f"-a TERA3={tera3adapter} --overlap 16 --error-rate 0.18 " \
                            "--rename '{id} TERAADAPTER3={adapter_name} {comment}' " \
                            f"{outputDir}/cat_files/cat.untrimmed.fastq > {outputDir}/cat_files/cat.fastq"
        else:
            warnings.warn(f"Please provide 5TERA and/or TERA3 adapters as strings!! "
                          f"You passed: {tera5adapter} and {tera3adapter}",
                          UserWarning)
            print(f"Skipping cutadapt, and moving backup file back to cat.fastq")
            call = f"mv {outputDir}/cat_files/cat.untrimmed.fastq {outputDir}/cat_files/cat.fastq"
            print(f"Starting undo of fastq backup at {get_dt(for_print=True)}\nUsing call:\t{call}\n")
            live_cmd_call(call)
            print(f"\n\nFinished undo of fastq backup at {get_dt(for_print=True)}")

        if isinstance(cutadapt_call, str):
            print(f"Starting cutadapt at {get_dt(for_print=True)}\nUsing call:\t{cutadapt_call}\n")
            live_cmd_call(cutadapt_call)
            print(f"\n\nFinished cutadapt at {get_dt(for_print=True)}")
    else:
        print(f"Skipping cutadapt for TERA-seq based on adapter comment(s) already being found in "
              f"cat.fastq file @ {outputDir}/cat_files/cat.fastq")


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


#################################################################################
# Step4: Concatenate files (less important for single MinIon runs), and create
#        a single file that contains information from all the tools I am using
#################################################################################
def __tera_adapter_tagging__(outputDir, tera3adapter, tera5adapter):
    import simplesam as ssam
    fastq_path = f'{outputDir}/cat_files/cat.fastq'
    
    # This open block is just to look at the first line of the fastq
    with open(fastq_path, 'r') as f:
        first_line = f.readline()
        # Check if each of the cutadapt comments were added, save for below.
        tera3_was_run = 'TERAADAPTER3' in first_line
        tera5_was_run = 'TERAADAPTER5' in first_line

        # If neither were run, just skip the rest of this method!
        if not tera3_was_run and not tera5_was_run:
            warnings.warn(f'Adapter comment not found in {fastq_path}, '
                          f'skipping adding TERA-seq tag to bam/sam files!')
            live_cmd_call(f"samtools view {outputDir}/cat_files/cat.sorted.bam "
                          f"> {outputDir}/cat_files/cat.sorted.sam")
            return None

    # If the adapter tag IS found in the cat.fastq, we'll add it to the bam/sam files!:
    tagged_fastq_df = FastqFile(fastq_path).df

    # For the two adapters, either extract the info if it's there, or default to false if not.
    if tera3_was_run:  # Parse out TERA3 adapter if it existed
        tagged_fastq_df['t3'] = tagged_fastq_df.comment.str.extract(r'TERAADAPTER3=(\S+)').replace({'no_adapter': '-',
                                                                                                    'TERA3': '+'})
    else:
        tagged_fastq_df['t3'] = '-'

    if tera5_was_run:  # Parse out TERA5 adapter if it existed
        tagged_fastq_df['t5'] = tagged_fastq_df.comment.str.extract(r'TERAADAPTER5=(\S+)').replace({'no_adapter': '-',
                                                                                                    'TERA5': '+'})
    else:
        tagged_fastq_df['t5'] = '-'

    # Finally we'll load and iterate through the bam file, creating a new sam file along the
    #   way and adding in the new tags!!
    tagged_fastq_df.set_index('read_id', inplace=True)
    input_bam = f'{outputDir}/cat_files/cat.sorted.bam'
    output_sam = f'{outputDir}/cat_files/cat.sorted.sam'
    print(f'Starting sam file tagging with TERA-seq adapter information @ {get_dt(for_print=True)}:')
    with ssam.Reader(open(input_bam, 'r')) as in_bam:
        with ssam.Writer(open(output_sam, 'w'), in_bam.header) as out_sam:
            row_iterator = tqdm(in_bam)
            for read in row_iterator:
                row_iterator.set_description(f"Tagging {read.qname}")
                read['t5'], read['t3'] = tagged_fastq_df.loc[read.qname, ['t5', 't3']].tolist()
                out_sam.write(read)
    print(f'Finished sam file tagging with TERA-seq adapter information @ {get_dt(for_print=True)}:')
    # Finally, we'll overwrite the old bam with the new, tagged sam file, and index it:
    call = f'samtools view -b {output_sam} -o {input_bam}'
    live_cmd_call(call)
    return None


def concat_files(outputDir, tera3adapter, tera5adapter, regenerate, **other_kwargs):
    original_bam_file = f"{outputDir}/cat_files/cat.sorted.bam"
    bam_file_with_only_mappedAndPrimary = f"{outputDir}/cat_files/cat.sorted.mappedAndPrimary.bam"
    if isinstance(tera3adapter, str) or isinstance(tera5adapter, str):
        __tera_adapter_tagging__(outputDir, tera3adapter, tera5adapter)
    else:
        # The below step has to happen in coordinance w/ re-tagging,
        #   so if re-tagging doesn't happen we still want to make the sam file!
        print(f"Unpacking bam file from:"
              f"\n\t{outputDir}/cat_files/cat.sorted.bam")
        live_cmd_call(f"samtools view {outputDir}/cat_files/cat.sorted.bam "
                      f"> {outputDir}/cat_files/cat.sorted.sam")
        print("Done.")

    concat_flag = regenerate or not path.exists(bam_file_with_only_mappedAndPrimary)
    if concat_flag:
        calls = [f"samtools view -b -F 0x904 {original_bam_file} > {bam_file_with_only_mappedAndPrimary}",
                 # The above command will build a new bam file w/out reads w/ bit_flags:
                 #    0x004, UNMAP           =   reads who's sequence didn't align to the genome
                 #    0x100, SECONDARY       =   reads that are secondary alignments
                 #    0x800, SUPPLEMENTARY   =   reads that are supplemental alignments
                 f"samtools index {bam_file_with_only_mappedAndPrimary}",
                 f"samtools view {outputDir}/cat_files/cat.sorted.mappedAndPrimary.bam "
                 f"> {outputDir}/cat_files/cat.sorted.mappedAndPrimary.sam",
                 ]
        print(f"Starting final cleanup at {get_dt(for_print=True)}\n")
        for num, call in enumerate(calls):
            print(f"\nStarting call ({num + 1} of {len(calls)}):\t{call}")
            live_cmd_call(call)
        print(f"\n\nFinished final cleanup at {get_dt(for_print=True)}")
    else:
        print(f"\n\nFile concatenation already ran. Based on file at:"
              f"\n\t{bam_file_with_only_mappedAndPrimary}\n"
              f"Use the regenerate tag if you want to rerun.\n")


def nanopolish_index_and_polya(genomeDir, dataDir, outputDir, threads, regenerate, **other_kwargs):
    # First we have to index the fastq_files!
    nanopolish_index_flag = regenerate or not path.exists(f"{outputDir}/cat_files/cat.fastq.index.readdb")
    if nanopolish_index_flag:
        seq_sum_matches = [i for i in listdir(dataDir) if
                           path.isfile(path.join(dataDir, i)) and 'sequencing_summary' in i]
        if len(seq_sum_matches) > 1:
            raise IndexError(f"Too many matches for 'sequencing_summary' in sequencing directory!!\n"
                             f"{seq_sum_matches}")
        if len(seq_sum_matches) == 0:
            warnings.warn(f"No matches for 'sequencing_summary' in sequencing directory!!\n"
                             f"{seq_sum_matches}\nGoing to continue w/out any seq summary! THIS WILL BE SLOW!")
            seq_sum_call = ""
        else:
            seq_sum_name = seq_sum_matches[0]
            seq_sum_call = f"--sequencing-summary={dataDir}/{seq_sum_name} "
        call = f"nanopolish index --directory={dataDir}/fast5 " \
               f"{seq_sum_call}" \
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
               f"--bam={outputDir}/cat_files/cat.sorted.mappedAndPrimary.bam --genome={genome_fa_file[0]} " \
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
# Step3: featureCounts to identify the genes that reads map to, and how many hits
#        we have per gene
#################################################################################
def feature_counts(genomeDir, outputDir, regenerate, threads, **other_kwargs):
    feature_counts_flag = regenerate or not path.exists(f"{outputDir}/featureCounts/"
                                                        f"cat.sorted.mappedAndPrimary.bam.featureCounts")
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
        print(f"\n\nfeatureCounts already occurred. Based on file at:"
              f"\n\t{outputDir}/featureCounts/[...]\n"
              f"Use the regenerate tag if you want to rerun featureCounts.\n")


def merge_results(**other_kwargs):
    def create_merge_df(outputDir, callWithJoshMethod, stepsToRun, genomeDir,
                        keep_multimaps=False, print_info=False,
                        **kwargs) -> pd.DataFrame:
        # First lets get the biggest one out of the way, importing the concatenated sam file:
        # 12/09/21: New SamOrBam class makes this wayyyy easier. Only downside is that
        #           it keeps all the flags I couldn't care less about!!
        sam_df = SamOrBamFile(f"{outputDir}/cat_files/cat.sorted.mappedAndPrimary.sam").df

        # Pull the 16 bit flag to get strand information (important for merge w/ featC later)
        sam_df["strand"] = (sam_df.bit_flag & 16).replace(to_replace={16: "-", 0: "+"})
        sam_df = sam_df.astype({'strand': 'category'})
        if keep_multimaps:
            # Identify and drop reads that have the 4 bit flag: indicating they didn't map!
            sam_df = sam_df[(sam_df.bit_flag & 4) != 4]
        else:
            # This shouldn't be dropping AS MANY reads now because I dumped secondary alignments with samtools
            sam_df = sam_df[~sam_df.duplicated(subset="read_id", keep=False)]

        if 'F' in stepsToRun:
            # Next lets pull in the featureCounts results
            featc_df = pd.read_csv(f"{outputDir}/featureCounts/cat.sorted.mappedAndPrimary.bam.Assigned.featureCounts",
                                   sep="\t", names=["read_id", "qc_tag_featc", "qc_pass_featc", "gene_id"])
            # Some reads will be ambiguously assigned to two different genes,
            #   for now I will drop these assignments (~1%):
            featc_df = featc_df.drop_duplicates()
            featc_df = featc_df[~featc_df.read_id.duplicated(keep=False)]  # The ~ inverts the filter :)
            # We can also load the gene_names dataframe and add it to the featureCounts one
            names_df = gene_names_to_gene_ids()
            featc_df = featc_df.merge(names_df, on="gene_id", how="left")

        # Load up nanopolish polyA results (also a population of duplicates here!!)
        polya_df = pd.read_csv(f"{outputDir}/nanopolish/polya.passed.tsv", sep="\t")
        # For some god-awful reason the chr_pos in polyA are -1 to those in the SAM file:
        polya_df["position"] = polya_df["position"] + 1
        polya_df = polya_df.rename(columns={"readname": "read_id",
                                            "qc_tag": "qc_tag_polya",
                                            "position": "chr_pos",
                                            "contig": "chr_id"})  # b/c featC also has a qc_tag!
        if print_info:
            print("#" * 100)
            print(f"\n\nSAM Dataframe info:")
            print(sam_df.info())
            if 'F' in stepsToRun:
                print(f"\n\nfeatureCounts Dataframe info:")
                print(featc_df.info())
            print(f"\n\nPolyA Dataframe info:")
            print(polya_df.info())
        # LETS SMOOSH THEM ALL TOGETHER!!!
        # This is severely broken in terms of the featureCounts merge:
        #       B/c the featureCounts output only retains read_id, it doesn't
        #       have enough information to merge uniquely for reads that map
        #       more than once!! This means that read A that maps to gene X and Y
        #       is eventually producing 4 lines of data....
        # 8/24/2021: Sorta fixed this by just dropping all duplicate reads! (line 484)
        #            Back to square one... lol
        # Revisiting on 10/26/2021: This is still broken. B/c the multiple
        #       hits (meaning multiple identical read_ids) in the bam/sam file are passed to
        #       feature counts, it propagates any multi-mappers. Currently I avoid this by
        #       dropping ANY read that hits more than once, meaning that propagation of
        #       multiple reads is avoided entirely. But it seems like I have a lot of good
        #       primary maps w/ trash secondaries! It would be really nice to retain those
        #       reads. . . Should I just switch to Josh's assignment method?
        #           OR: I could try to rename reads in the bam file w/ their map location,
        #           as this could help to uniquely identify primaries, and that info would
        #           get passed though featureCounts!
        # 10/28/2021: We'll call this provisionally solved, b/c I used some samtools features above
        #             to exclusively pass mapped, primary reads to featureCounts

        # These steps now retain reads w/out gene assignments or polyA tail calls!
        #   This filtering should be easy to do in later scripts.
        if "F" in stepsToRun:
            sam_featc_df = sam_df.merge(featc_df, how="left", on=["read_id"])
            merge_reads_df = sam_featc_df.merge(polya_df, how="left", on=["read_id", "chr_id", "chr_pos"])
        else:
            merge_reads_df = sam_df.merge(polya_df, how="left", on=["read_id", "chr_id", "chr_pos"])
        merge_reads_df = merge_reads_df.drop_duplicates()
        
        # Dropping unmapped or secondary reads, after the addition of the -F 0x904 with samtools this should do nothing:
        merge_reads_df = merge_reads_df[merge_reads_df["sequence"] != "*"]
        
        # Dropping terrible mapq scored reads, I don't think there are very many of these at all(?)
        merge_reads_df = merge_reads_df[merge_reads_df["mapq"] != 0]

        if callWithJoshMethod:
            # TODO: More standardized handling of this or the featureCounts called gene_ids.
            #       The current system is inconsistent and confuses downstream scripts!
            # TODO: Another option would be to just settle on one method for use!!
            merge_reads_df = assign_with_josh_method(merge_reads_df, genomeDir,
                                                     keepMultipleTranscriptInfo=False)
            if 'F' in stepsToRun:
                print(f"Reads that have matching assignments: "
                      f"{merge_reads_df[merge_reads_df.gene_id == merge_reads_df.gene_id_fromFeatureCounts].shape[0]}/"
                      f"{merge_reads_df.shape[0]}")
            print(f"Reads without Josh assignments: "
                  f"{merge_reads_df[merge_reads_df.gene_id.isna()].shape[0]}/"
                  f"{merge_reads_df.shape[0]}")
            print(f"For now (as of 12/13/2021) I'm keeping all of the above!")

        merge_reads_df['read_length'] = merge_reads_df['sequence'].apply(len)
        if print_info:
            print("\n\n")
            print("#" * 100)
            print(f"\n\nMerged Dataframe info:")
            print(merge_reads_df.info())
        if "gene_id" not in merge_reads_df.columns:
            raise NotImplementedError(f"Not genes_ids mapped! Please either run featureCounts,"
                                      f"or use the -j option to allow for using Josh's read assignment method!")
        merge_out_file = f"{outputDir}/merge_files/{get_dt(for_file=True)}_mergedOnReads"
        print(f"Saving compressed on reads files to:\n\t{merge_out_file} + .parquet/.tsv")
        merge_reads_df.to_csv(merge_out_file + ".tsv", sep="\t", index=False)
        # Added 10/28/2021: Parquet files are SOOOOOO much lighter and faster
        merge_reads_df.to_parquet(merge_out_file + ".parquet")
        return merge_reads_df

    def compress_on_genes(merged_df, outputDir, callWithJoshMethod,
                          dropGeneWithHitsLessThan=None, output_to_file=True,
                          **kwargs) -> pd.DataFrame:
        # This creates a pandas "groupby" object that can be used to extract info compressed on gene_ids
        print("\nMerging information from Minimap2, featureCounts and Nanopolish-PolyA:")
        if 'read_length' not in list(merged_df.columns):
            merged_df["read_length"] = merged_df["sequence"].str.len()
        for adapter_col in ['t5', 't3']:
            if adapter_col in merged_df.columns:
                merged_df[adapter_col].replace({'+': 1, '-': 0}, inplace=True)
        grouped_genes = merged_df.groupby(["gene_id", "gene_name"])

        # This next step now uses set rather than list, to ensure we count each read only once!!
        gene_df = grouped_genes["read_id"].apply(set).apply(len).to_frame(name="read_hits")
        gene_df["read_ids"] = grouped_genes["read_id"].apply(list).to_frame(name="read_ids")

        gene_df["read_len_mean"] = grouped_genes["read_length"].apply(np.mean).to_frame(name="read_len_mean")
        gene_df["read_len_std"] = grouped_genes["read_length"].apply(np.std).to_frame(name="read_len_std")
        gene_df["read_lengths"] = grouped_genes["read_length"].apply(list).to_frame(name="read_lengths")

        gene_df["polya_lengths"] = grouped_genes["polya_length"].apply(list).to_frame(name="polya_lengths")
        gene_df["polya_mean"] = grouped_genes["polya_length"].apply(np.mean).to_frame(name="polya_mean")
        gene_df["polya_stdev"] = grouped_genes["polya_length"].apply(np.std).to_frame(name="polya_stdev")

        gene_df["genomic_starts"] = grouped_genes["chr_pos"].apply(list).to_frame(name="genomic_starts")
        gene_df["cigars"] = grouped_genes["cigar"].apply(list).to_frame(name="cigars")

        for adapter_col in ['t5', 't3']:
            if adapter_col in merged_df.columns:
                gene_df[f"{adapter_col}_fraction"] = grouped_genes[adapter_col].sum() / grouped_genes[
                    adapter_col].apply(list).apply(len)

        if dropGeneWithHitsLessThan:
            print(f"Dropping any genes with less than {dropGeneWithHitsLessThan} read hits")
            gene_df = gene_df[gene_df["read_hits"] >= dropGeneWithHitsLessThan]
        gene_df.sort_values("read_hits", ascending=False, inplace=True)
        # print(f"Mean PolyA Length: {gene_df['polya_mean'].mean():.3f}")
        gene_df.reset_index(inplace=True)

        if output_to_file:
            filename = f"{get_dt(for_file=True)}_compressedOnGenes"
            output_file = f"{outputDir}/merge_files/{filename}"
            print(f"Saving compressed on genes files to:\n\t{output_file} + .parquet/.tsv")
            gene_df.to_parquet(f"{output_file}.parquet")
            light_gene_df = gene_df.drop(['read_ids',
                                          'polya_lengths',
                                          'read_lengths',
                                          'genomic_starts',
                                          'cigars'],
                                         axis=1)
            light_gene_df.to_csv(f"{output_file}_simple.tsv", sep="\t", index=False)
            light_gene_df.to_parquet(f"{output_file}_simple.parquet")
        return gene_df

    print(f"Starting to merge all data at {get_dt(for_print=True)}\n")
    merge_df = create_merge_df(**other_kwargs)
    print(f"\n\nFinished merging all data at {get_dt(for_print=True)}")
    print(f"Starting to compress data on genes at {get_dt(for_print=True)}\n")
    genes_df = compress_on_genes(merge_df, **other_kwargs)
    print(f"\n\nFinished compressing data on genes at {get_dt(for_print=True)}")
    return merge_df, genes_df


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
    stds_mini_df = align_standards(compressed_df=df, keep_read_id=True, **other_kwargs)
    if isinstance(stds_mini_df, pd.DataFrame):
        df = df.merge(stds_mini_df, on="read_id")
        out_file = f"{outputDir}/merge_files/{get_dt(for_file=True)}_mergedOnReads.plusStandards.parquet"
        print(f"Saving parquet file to:\n\t{out_file}")
        df.to_parquet(out_file)
        plot_value_counts(df)
    else:
        print(stds_mini_df)


def flair(outputDir, **other_kwargs):
    def run_and_load_flair(outputDir, genomeDir, threads, sampleID,
                           condition, regenerate, **kwargs) -> pd.DataFrame:
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
    final_flair_df = merge_some_more(flair_df, outputDir,
                                     output_to_file=True)
    return final_flair_df


def extra_steps(outputDir, genomeDir, minimapParam, threads, df=None, **other_kwargs):
    # adjust_5_ends(df, genomeDir, outputDir)
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
    print(f"MiniMap2 Call: {call}")


def main(stepsToRun, **kwargs) -> (pd.DataFrame, pd.DataFrame) or None:
    return_value = None
    buildOutputDirs(stepsToRun, **kwargs)
    kwargs['stepsToRun'] = stepsToRun

    # Current default is [GMNCFPS]
    steps_dict = {"G": [guppy_basecall_w_gpu, "Guppy Basecalling"],
                  "A": [alternative_genome_filtering, "Filtering Alt. Genomes (currently pretty rough and slow)"],
                  "T": [trim_tera_adapters, "Trimming TERA-Seq Adapters"],
                  "M": [minimap2_and_samtools, "Minimap2 and SamTools"],
                  "C": [concat_files, "Concatenate Files"],
                  "N": [nanopolish_index_and_polya, "Nanopolish Index and polyA Calling"],
                  "F": [feature_counts, "FeatureCounts"],
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
                return_value = step_func(**kwargs)
            elif code == "X" or code == "S":
                if return_value:  # The extra steps function will accept a df, or load from the disk
                    step_func(df=return_value[0], **kwargs)
                else:
                    step_func(**kwargs)
            else:  # The rest of the functions work by loading and saving to the disk
                step_func(**kwargs)
        else:
            step_print = f"Skipping step: \"{step_name}\" . . ."
            print("\n", step_print, f"\n{'#' * len(step_print)}\n", sep="")
    return return_value


if __name__ == '__main__':
    arguments = meshSetsAndArgs()
    main_output = main(**arguments)
