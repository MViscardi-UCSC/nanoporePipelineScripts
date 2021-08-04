"""
#### RUN AS ROOT W/ ENV (sudo -E) #### Maybe? I dunno. Running as user may work fine
step0_nanopore_pipeline.py
Marcus Viscardi     May 29, 2021

I have a ~working set of shell scripts that can go from raw ONT fast5 files
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
from datetime import datetime
from typing import List
from glob import glob

import pandas as pd
import numpy as np


def get_dt(for_print=False, for_output=False):
    now = datetime.now()
    if for_print:
        return str(now.strftime("%m/%d/%y @ %I:%M:%S %p"))
    elif for_output:
        return str(now.strftime("%y%m%d"))
    else:
        return str(now.strftime("%y%m%d_%I:%M:%S%p"))


def live_cmd_call(command):
    with Popen(command, stdout=PIPE, shell=True,
               bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end="")
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)


def find_newest_matching_file(path_str):
    # A tool that I'll want to use to grab the most recent file
    from glob import glob
    import os

    list_of_files = glob(path_str)
    try:
        latest_file = max(list_of_files, key=path.getctime)
        return latest_file
    except ValueError:
        raise ValueError(f"Failed to find any files matching \"{path_str}\"")


def gene_names_to_gene_ids(tsv_path: str = "/home/marcus/PycharmProjects/"
                                       "PersonalScripts/WBGene_to_geneName.tsv") -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t")
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
                            help="Number of threads to be used by nanopolish and minimap2.")
        parser.add_argument('--guppyConfig', metavar='guppy_config', type=str,
                            default='rna_r9.4.1_70bps_hac.cfg',
                            help="Configuration preset passed to the guppy_basecaller "
                                 "based on flowcell and kit used for run. Helpful "
                                 "table for picking a config @ https://denbi-nanopore-"
                                 "training-course.readthedocs.io/en/latest/basecalling/"
                                 "basecalling.html or just google 'guppy_basecalling'")
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
                            help="Steps to run within the pipeline (Normal Run = GMNFCP): "
                                 "(G)uppy basecalling, (M)inimap, (N)anopolish, "
                                 "(F)eature counts, (C)oncat files, merge with (P)andas, "
                                 "and/or use f(L)air to ID transcripts (not default behavior)")
        parser.add_argument('--sampleID', metavar='sampleID',
                            type=int, default=None,
                            help="sampleID to pass to FLAIR")
        parser.add_argument('--condition', metavar='condition',
                            type=int, default=None,
                            help="condition to pass to FLAIR")
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
    absoluteDefDict = {"printArgs": False,
                       "nestedData": False,
                       "regenerate": False,
                       "altGenomeDirs": [],
                       "threads": 20,
                       "stepsToRun": "GMNFCP",
                       "sampleID": "sample1",
                       "condition": "conditionA"}
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


def buildOutputDirs(outputDir, **kwargs) -> None:
    dirs_list = [outputDir,
                 "fastqs",
                 "cat_files",
                 "nanopolish",
                 "featureCounts",
                 "logs",
                 "merge_files",
                 "flair"
                 ]
    print('\n')
    for dir_to_make in dirs_list:
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
                      rf"""-v {dataDir}:/usr/src/working/data_dir """
                      rf"""-v {outputDir}:/usr/src/working/output_dir """
                      rf"""-it nanopore_empty """
                      rf"""guppy_basecaller --num_callers {callers} --cpu_threads_per_caller {threads_per_caller} """
                      rf"""-c {guppyConfig} -i /usr/src/working/data_dir/fast5 -s /usr/src/working/output_dir/fastqs """
                      rf"""2>&1 | tee {outputDir}/logs/{get_dt()}_guppy.log""")
        live_cmd_call(rf"cat {outputDir}/fastqs/*.fastq > {outputDir}/cat_files/cat.fastq")
    else:
        print(f"\n\nCalling already occurred. Based on file at:\n\t{outputDir}/cat_files/cat.fastq\n"
              f"Use the regenerate tag if you want to rerun calling.\n")


def guppy_basecall_w_gpu(dataDir, outputDir, threads, guppyConfig, regenerate, **other_kwargs):
    prev_cat_fastq = path.exists(f"{outputDir}/cat_files/cat.fastq")
    if regenerate or not prev_cat_fastq:
        # TODO: optimize the num_callers and gpu_runners_per_caller params!!
        live_cmd_call(rf"""guppy_basecaller -x "cuda:0" --num_callers 12 --gpu_runners_per_device 8 """
                      rf"""-c {guppyConfig} -i {dataDir}/fast5 -s {outputDir}/fastqs """
                      rf"""2>&1 | tee {outputDir}/logs/{get_dt()}_guppy.log""")
        live_cmd_call(rf"cat {outputDir}/fastqs/pass/*.fastq > {outputDir}/cat_files/cat.fastq")
    else:
        print(f"\n\nCalling already occurred. Based on file at:\n\t{outputDir}/cat_files/cat.fastq\n"
              f"Use the regenerate tag if you want to rerun calling.\n")


#################################################################################
# Step2: Nanopolish index, map to genome with minimap2, samtools index/sort, &
#        Nanopolish calls of polyA tail lengths TODO: split minimap and nanopolish stuff
#################################################################################
def alternative_genome_filtering():
    pass


def minimap_and_nanopolish(genomeDir, dataDir, outputDir, threads, regenerate, **other_kwargs):
    if regenerate or not path.exists(f"{outputDir}/cat_files/cat.fastq.index.readdb"):
        seq_sum_matches = [i for i in listdir(dataDir) if
                           path.isfile(path.join(dataDir, i)) and 'sequencing_summary' in i]
        if len(seq_sum_matches) != 1:
            raise IndexError(f"Too many/few matches for 'sequencing_summary' in sequencing directory!!\n"
                             f"{seq_sum_matches}")
        else:
            seq_sum_name = seq_sum_matches[0]
        call = f"nanopolish index --directory={dataDir}/fast5 " \
               f"--sequencing-summary={dataDir}/{seq_sum_name} " \
               f"{outputDir}/cat_files/cat.fastq"  # This still needs testing!
        print(f"Starting nanopolish index at {get_dt(for_print=True)}\nUsing call:\t{call}\n"
              f"(There are limited outputs from this script, and it runs very slow)")
        live_cmd_call(call)
        print(f"Nanopolish index completed at {get_dt(for_print=True)}")
    else:
        print(f"\n\nNanopolish index already ran. Based on file at:"
              f"\n\t{outputDir}/cat_files/cat.fastq.index.readdb\n"
              f"Use the regenerate tag if you want to rerun.\n")

    if regenerate or not path.exists(f"{outputDir}/cat_files/cat.bam"):
        # TODO: Logging for the mapping call with minimap2
        genome_fa_file = glob(f"{genomeDir}/*allChrs.fa")
        if len(genome_fa_file) != 1:
            raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                      f"with one fa file that ends with 'allChrs.fa'")
        genome_bed_file = glob(f"{genomeDir}/*.bed")
        if len(genome_bed_file) != 1:
            raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                      f"with one bed file that ends with '.bed'")
        call = f"minimap2 -a -x map-ont {genome_fa_file[0]} {outputDir}/cat_files/cat.fastq " \
               f"-t {threads} --junc-bed {genome_bed_file[0]} | samtools view -b - -o " \
               f"{outputDir}/cat_files/cat.bam"
        print(f"Starting minimap2 at {get_dt(for_print=True)}\nUsing call:\t{call}\n")
        live_cmd_call(call)
        print(f"\n\nFinished minimap2 at {get_dt(for_print=True)}\n")
    else:
        print(f"\n\nMiniMap2 already ran. Based on file at:"
              f"\n\t{outputDir}/cat_files/cat.bam\n"
              f"Use the regenerate tag if you want to rerun.\n")

    if regenerate or not path.exists(f"{outputDir}/cat_files/cat.sorted.bam"):
        call = f"samtools sort -T tmp -o {outputDir}/cat_files/cat.sorted.bam " \
               f"{outputDir}/cat_files/cat.bam && samtools index " \
               f"{outputDir}/cat_files/cat.sorted.bam"
        print(f"Starting samtools sort & index at {get_dt(for_print=True)}\nUsing call:\t{call}\n")
        live_cmd_call(call)
        print(f"\n\nFinished samtools sort & index at {get_dt(for_print=True)}")
    else:
        print(f"\n\nsamtools sort already ran. Based on file at:"
              f"\n\t{outputDir}/cat_files/cat.sorted.bam\n"
              f"Use the regenerate tag if you want to rerun.\n")

    if regenerate or not path.exists(f"{outputDir}/nanopolish/polya.tsv"):
        genome_fa_file = glob(f"{genomeDir}/*allChrs.fa")
        if len(genome_fa_file) != 1:
            raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                      f"with one fa files that ends with 'allChrs.fa'")
        call = f"nanopolish polya --threads={threads} --reads={outputDir}/cat_files/cat.fastq " \
               f"--bam={outputDir}/cat_files/cat.sorted.bam --genome={genome_fa_file[0]} " \
               f"> {outputDir}/nanopolish/polya.tsv"
        print(f"Starting nanopolish polyA at {get_dt(for_print=True)}\nUsing call:\t{call}\n")
        live_cmd_call(call)
        print(f"\n\nFinished nanopolish polyA at {get_dt(for_print=True)}")
    else:
        print(f"\n\nNanopolish polyA already ran. Based on file at:"
              f"\n\t{outputDir}/nanopolish/polya.tsv\n"
              f"Use the regenerate tag if you want to rerun.\n")


#################################################################################
# Step3: featureCounts to identify the genes that reads map to, and how many hits
#        we have per gene
#################################################################################
def feature_counts(genomeDir, outputDir, regenerate, threads, **other_kwargs):
    # feature_counts_file = os.path.exists(f"{outputDir}/featureCounts/")
    if regenerate or not path.exists(f"{outputDir}/featureCounts/cat.sorted.bam.featureCounts"):
        genome_gtf_file = glob(f"{genomeDir}/*.gtf")
        if len(genome_gtf_file) != 1:
            raise NotImplementedError(f"Currently this script only supports having genomeDirs "
                                      f"with one gtf files that ends with '.gtf'")
        call = f"featureCounts -L -T {threads} -R CORE -a {genome_gtf_file[0]} " \
               f"-o {outputDir}/featureCounts/{get_dt()} {outputDir}/cat_files/cat.sorted.bam " \
               f"2>&1 | tee -a {outputDir}/logs/{get_dt()}.featureCounts.log"
        print(f"Starting featureCounts at {get_dt(for_print=True)}\nUsing call:\t{call}\n")
        live_cmd_call(call)
        print(f"\n\nFinished featureCounts at {get_dt(for_print=True)}")
    else:
        print(f"\n\nPost calling processing already occurred. Based on file at:"
              f"\n\t{outputDir}/featureCounts/[...]\n"
              f"Use the regenerate tag if you want to rerun featureCounts.\n")


#################################################################################
# Step4: Concatenate files (less important for single MinIon runs), and create
#        a single file that contains information from all the tools I am using
#################################################################################
def final_touches(outputDir, **other_kwargs):
    # Most of this is much less necessary as we are not getting the nested files that came out of Roach's gridION!
    calls = [f"samtools view {outputDir}/cat_files/cat.sorted.bam > {outputDir}/cat_files/cat.sorted.sam",
             f"head -n 1 {outputDir}/nanopolish/polya.tsv > {outputDir}/nanopolish/polya.passed.tsv",
             f"grep PASS {outputDir}/nanopolish/polya.tsv >> {outputDir}/nanopolish/polya.passed.tsv",
             f"grep Assigned {outputDir}/featureCounts/cat.sorted.bam.featureCounts "
             f">> {outputDir}/featureCounts/cat.sorted.bam.Assigned.featureCounts",
             ]
    print(f"Starting final cleanup at {get_dt(for_print=True)}\n")
    for num, call in enumerate(calls):
        print(f"\nStarting call ({num + 1} of {len(calls)}):\t{call}")
        live_cmd_call(call)
    print(f"\n\nFinished final cleanup at {get_dt(for_print=True)}")


def merge_results(**other_kwargs):
    def create_merge_df(outputDir, print_info=False, **kwargs) -> pd.DataFrame:
        # First lets get the biggest one out of the way, importing the concatenated sam file:
        sam_df = pd.read_csv(f"{outputDir}/cat_files/cat.sorted.sam",
                             sep="\t", names=range(23))
        # Drop any read_ids that have come up multiple times:
        #   TODO: I am assuming these are multiply mapping? Is this wrong?
        sam_df = sam_df.drop_duplicates(subset=0, keep=False, ignore_index=True)
        # I have no idea what the additional tags from Minimap2 are, I'm dropping them for now:
        sam_df = sam_df[range(11)]
        # And lets rename columns while we are at it!
        sam_header_names = ["read_id",
                            "bit_flag",
                            "chr_id",
                            "chr_pos",
                            "mapq",
                            "cigar",
                            "r_next",
                            "p_next",
                            "len",
                            "sequence",
                            "phred_qual"]
        sam_df = sam_df.rename(columns=dict(enumerate(sam_header_names)))
        # Next lets pull in the featureCounts and nanopolish polyA results
        featc_df = pd.read_csv(f"{outputDir}/featureCounts/cat.sorted.bam.Assigned.featureCounts",
                               sep="\t", names=["read_id", "qc_tag", "qc_pass", "gene_id"])
        polya_df = pd.read_csv(f"{outputDir}/nanopolish/polya.passed.tsv", sep="\t")
        polya_df = polya_df.rename(columns={"readname": "read_id"})
        if print_info:
            print("#" * 100)
            print(f"\n\nSAM Dataframe info:")
            print(sam_df.info())
            print(f"\n\nfeatureCounts Dataframe info:")
            print(featc_df.info())
            print(f"\n\nPolyA Dataframe info:")
            print(polya_df.info())
        # LETS SMOOSH THEM ALL TOGETHER!!!
        sam_featc_df = sam_df.merge(featc_df, how="inner", on="read_id")
        merge_df = sam_featc_df.merge(polya_df, how="inner", on="read_id")
        merge_df.drop(["contig", "position", "r_next", "p_next", "len"], axis=1)
        merge_df = merge_df.drop_duplicates()
        if print_info:
            print("\n\n")
            print("#" * 100)
            print(f"\n\nMerged Dataframe info:")
            print(merge_df.info())
        with open(f"{outputDir}/merge_files/{get_dt(for_output=True)}_mergedOnReads.tsv", "w") as merge_out_f:
            merge_df.to_csv(merge_out_f, sep="\t", index=False)
        return merge_df

    def compress_on_genes(merged_df, outputDir, dropGeneWithHitsLessThan=None,  # dip_test=False,
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

        if dropGeneWithHitsLessThan:
            print(f"Dropping any genes with less than {dropGeneWithHitsLessThan} read hits")
            gene_df = gene_df[gene_df["read_hits"] >= dropGeneWithHitsLessThan]
        gene_df.sort_values("read_hits", ascending=False, inplace=True)
        # if dip_test:
        #     print(f"Starting dip test at {dt.now().strftime('%H:%M:%S')}")
        #     gene_df["dip_test"] = grouped_genes["polya_length"].apply(quick_dip).to_frame(name="dip_test")
        #     gene_df["modality"] = gene_df["dip_test"].apply(len)
        #     print(f"Finished dip test at {dt.now().strftime('%H:%M:%S')}")
        print(f"Mean PolyA Length: {gene_df['polya_mean'].mean():.3f}")
        filename = f"{get_dt(for_output=True)}_compressedOnGenes"
        if output_to_file:
            with open(f"{outputDir}/merge_files/{filename}.tsv", "w") as out_f:
                gene_df.to_csv(out_f, sep="\t")
            with open(f"{outputDir}/merge_files/{filename}_simple.tsv", "w") as out_f:
                gene_df.drop(["read_ids", "polya_lengths", 'read_lengths'], axis=1).to_csv(out_f, sep="\t")
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
        merge_tsv_path = f"{outputDir}/merge_files/{get_dt(for_output=True)}_mergedOnReads.tsv"
        if not path.exists(merge_tsv_path):
            older_merge_tsv = find_newest_matching_file(f"{outputDir}/merge_files/*_mergedOnReads.tsv")
            print(f"Could not find file @ {merge_tsv_path}\n"
                  f"Going to use the older file @ {older_merge_tsv}")
            merge_tsv_path = older_merge_tsv
        merge_df = pd.read_csv(merge_tsv_path, sep="\t")
        print(f"\nLoaded my merge file. . .\n")
        print(merge_df.info())
        super_df = merge_df.merge(flair_df, how="inner", on="read_id")
        super_df.to_csv(f"{outputDir}/merge_files/{get_dt(for_output=True)}_mergedWithTranscripts.tsv",
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
        filename = f"{get_dt(for_output=True)}_compressedOnTranscripts"
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


def main(stepsToRun, **kwargs) -> (pd.DataFrame, pd.DataFrame) or None:
    return_value = None
    buildOutputDirs(**args)

    steps_dict = {"G": [guppy_basecall_w_gpu, "Guppy Basecalling"],
                  "A": [alternative_genome_filtering, "Filtering Alt. Genomes (not implemented)"],
                  # TODO: Add alternative genome filter mapping function here^^^
                  "M": [minimap_and_nanopolish, "Minimap2 & Nanopolish"],  # Also N until I split them up!
                  "F": [feature_counts, "FeatureCounts"],
                  "C": [final_touches, "Concatinate Files"],
                  "P": [merge_results, "Merging Results w/ Pandas"],
                  "L": [flair, "Calling Transcripts w/ Flair"]
                  }
    for code, (step_func, step_name) in steps_dict.items():
        if code in stepsToRun:
            step_print = f"Starting step: \"{step_name}\" . . ."
            print("\n\n", step_print, f"\n{'#' * len(step_print)}\n", sep="")
            if code != "P":
                step_func(**args)
            else:
                return_value = step_func(**args)
        else:
            step_print = f"Skipping step: \"{step_name}\" . . ."
            print("\n", step_print, f"\n{'#' * len(step_print)}\n", sep="")
    return return_value


if __name__ == '__main__':
    args = meshSetsAndArgs()
    main_output = main(**args)
