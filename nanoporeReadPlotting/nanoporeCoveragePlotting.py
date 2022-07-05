"""
nanoporeCoveragePlotting.py
Marcus Viscardi,    July 05, 2022

Most of this is from numpySpeedTesting_andCoveragePlotting.ipynb,
    as this is where I got the really fast version of this working.
    This script is mostly to have it as a self-contained tool that
    can be given to other people and used easily. Doing so will likely
    require making the TERA part optional, and introducing an option to
    plot a region of a chromosome rather than a GTF specified gene.
"""
import numpy as np
import pandas as pd
import seaborn as sea
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

import re
from nanoporePipelineCommon import *

# TODO: Write this script...
#  Main efforts will be:
#   - adding CLI interface with argParse
#   - adding library loading step
#   - dropping 5TERA requirement
#   - introducing option to plot by locus, rather than by gene

def add_to_main_array_for_each_read(cigars_and_genomic_starts, chr_length, count_Ds_as_maps=False) -> np.array:
    coverage_array = np.zeros([chr_length], dtype=np.uint32)
    cigar_parsing_iterator = tqdm(cigars_and_genomic_starts, desc=f"Building coverage by adding to main array")

    gaps = ['N']  # introns
    maps = ['M']  # mapped segments
    if count_Ds_as_maps:
        maps.append('D')
    else:
        gaps.append('D')
    for cigar, genomic_start in cigar_parsing_iterator:
        genomic_pos = genomic_start
        parsed_cigar = re.findall(rf'(\d+)([MDNSIX])', cigar)
        parsed_cigar = [(int(length), code) for (length, code) in parsed_cigar]

        for length, code in parsed_cigar:
            if code in gaps:
                genomic_pos += length
            elif code in maps:
                coverage_array[genomic_pos:genomic_pos + length] += 1
                genomic_pos += length
    return coverage_array


def run_coverage_techiques(bam_df, chrs=("I", "II", "III", "IV", "V", "MtDNA"), count_Ds_as_maps=False):
    if "original_chr_pos" in bam_df.columns:
        gen_pos_col = "original_chr_pos"
    else:
        gen_pos_col = "chr_pos"

    # These are actually annotation ends:
    chr_max_length_dict = {'I': 15072426,
                           'II': 15279420,
                           'III': 13783459,
                           'IV': 17493829,
                           'V': 20922738,
                           'MtDNA': 13327}
    
    # Because of some python weirdness, we need to turn single chromosomes into lists!
    if not isinstance(chrs, (list, tuple)):
        chrs = [chrs]

    # First filter the chr dict, so we only use the ones that showed up in the method call:
    chr_max_length_dict = {chromosome: length for chromosome, length
                           in chr_max_length_dict.items()
                           if chromosome in chrs}
    array_dict = {}
    for chromosome, chr_length in chr_max_length_dict.items():
        print(f"Starting to build coverage array for chromosome: {chromosome}")
        chr_df = bam_df.query(f"chr_id == '{chromosome}'")
        cigars_and_genomic_start_positions = list(zip(chr_df.cigar.to_list(), chr_df[gen_pos_col].to_list()))
        coverage_array = add_to_main_array_for_each_read(cigars_and_genomic_start_positions, chr_length,
                                                         count_Ds_as_maps=count_Ds_as_maps)
        array_dict[chromosome] = coverage_array
    return array_dict


def get_gene_coordinates(gene_id=None, gene_name=None,
                         parsed_gtf_path="/data16/marcus/genomes/elegansRelease100/"
                                         "Caenorhabditis_elegans.WBcel235.100.gtf.parquet"
                         ) -> (str, int, int):
    # First make sure we got something to look up:
    gene_id_bool = isinstance(gene_id, str)
    gene_name_bool = isinstance(gene_name, str)
    if not gene_id_bool and not gene_name_bool:
        raise NotImplementedError(f"Please pass a gene_id or a gene_name!")
    # Load the parsed gtf_file
    try:
        gtf_df = pd.read_parquet(parsed_gtf_path)[["gene_id",
                                                   "gene_name",
                                                   "feature",
                                                   "chr",
                                                   "start",
                                                   "end"]].query("feature == 'gene'")
    except FileNotFoundError:
        raise FileNotFoundError(f"Please make sure there is a parsed gtf file at: {parsed_gtf_path}")

    # Get the gene of interest!
    try:
        if gene_id_bool:
            entry_of_interest = gtf_df.query(f"gene_id == '{gene_id}'").reset_index(drop=True).iloc[0].to_dict()
            gene_name = entry_of_interest["gene_name"]
        else:  # if gene_name_bool
            entry_of_interest = gtf_df.query(f"gene_name == '{gene_name}'").reset_index(drop=True).iloc[0].to_dict()
            gene_id = entry_of_interest["gene_id"]
    except IndexError:
        raise IndexError(f"Gene of interest (gene_id: {gene_id} / gene_name: {gene_name}) not found!")
    chromosome = entry_of_interest["chr"]
    start = entry_of_interest["start"]
    end = entry_of_interest["end"]
    print(f"Found entry for {gene_name} ({gene_id}) on chromosome {chromosome:>5} at ({start}, {end})")
    return chromosome, start, end


def coverage_plotting_5tera(bam_df_for_plot, gene_name, save_as=None, count_Ds_as_maps=False, rpm_normalize=False):
    gene_chr, start, end = get_gene_coordinates(gene_name=gene_name)
    chr_array_ad = run_coverage_techiques(bam_df_for_plot.query("t5 == '+'"),
                                          chrs=gene_chr,
                                          count_Ds_as_maps=count_Ds_as_maps)[gene_chr]
    chr_array_unad = run_coverage_techiques(bam_df_for_plot.query("t5 == '-'"),
                                            chrs=gene_chr,
                                            count_Ds_as_maps=count_Ds_as_maps)[gene_chr]
    locus_array_ad = chr_array_ad[start: end]
    locus_array_unad = chr_array_unad[start: end]
    if rpm_normalize:
        norm_factor = bam_df_for_plot.shape[0]
        # Turn the total number of read hits into the 'million of read hits'
        rpm_norm_factor = norm_factor / 1000000
        locus_array_ad = np.divide(locus_array_ad, rpm_norm_factor)
        locus_array_unad = np.divide(locus_array_unad, rpm_norm_factor)

    index_array = np.arange(len(locus_array_ad))
    zeros_array = np.zeros(len(locus_array_ad))

    fig, ax = plt.subplots(2, 1,
                           figsize=(8, 4),
                           sharex=True,
                           gridspec_kw={"height_ratios": [1, 4]}
                           )

    ax[1].fill_between(index_array, zeros_array, locus_array_unad, color='black')
    ax[0].fill_between(index_array, zeros_array, locus_array_ad, color='red')

    plt.tight_layout()
    plt.style.use("seaborn-paper")
    if isinstance(save_as, str):
        plt.savefig(save_as)
    plt.show()


