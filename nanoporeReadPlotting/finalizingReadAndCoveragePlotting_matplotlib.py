"""
finalizingReadAndCoveragePlotting_matplotlib.py
Marcus Viscardi,    January 16, 2023

This is really just taking the methods from finalizingReadPlotting_matplotlib.ipynb and
numpySpeedTesting_andCoveragePlotting.ipynb, then setting it up in a way to be import-able!
"""
import numpy as np
import pandas as pd
from tqdm import tqdm

import re

import seaborn as sea
import matplotlib.pyplot as plt
import sys

import warnings

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
sys.path.insert(0, '/data16/marcus/scripts/nanoporePipelineScripts')
from nanoporePipelineCommon import *
pd.set_option('display.width', 100)
pd.set_option('display.max_columns', None)


def _make_rectangle_patch(genome_start, length, y_center, thickness, color='gray'):
    return Rectangle((genome_start, y_center - (thickness / 2)), length, thickness,
                     facecolor=color,
                     edgecolor=color,
                     fill=True,
                     lw=0)


def _add_patches_from_cigars_and_gen_pos(axes, cigar, gen_start, y, strand, color='black', plot_introns=True,
                                         tail_length=None):
    # Parse the cigar string
    parsed_cigar = re.findall(rf'(\d+)([MDNSIX])', cigar)
    mdn_nums = [int(num) for num, char in parsed_cigar if char in "MDN"]
    gen_end = gen_start + sum(mdn_nums)
    genomic_read_length = gen_end - gen_start

    genome_loc = gen_start

    rectangle_patch_list = []
    first_n_length = 0
    for length, code in parsed_cigar:
        length = int(length)
        if code == 'S':
            pass
        elif code == 'M':
            rectangle_patch_list.append(_make_rectangle_patch(genome_loc, length, y, thickness=0.8))
            genome_loc += length
        elif code == 'D':
            if length < 50:
                rectangle_patch_list.append(_make_rectangle_patch(genome_loc, length, y, thickness=0.8))
            else:
                if plot_introns:
                    rectangle_patch_list.append(_make_rectangle_patch(genome_loc, length, y, thickness=0.001))
            genome_loc += length
        elif code == 'I':
            pass
        elif code == 'N':
            if plot_introns:
                rectangle_patch_list.append(_make_rectangle_patch(genome_loc, length, y, thickness=0.001))
            genome_loc += length
    axes.add_collection(PatchCollection(rectangle_patch_list, color=color))
    if isinstance(tail_length, float):
        if strand == "+":
            axes.add_patch(_make_rectangle_patch(genome_loc, tail_length, y, thickness=0.4, color='green'))
            genome_loc += tail_length
        else:
            axes.add_patch(_make_rectangle_patch(gen_start, -tail_length, y, thickness=0.4, color='green'))
    return genomic_read_length


def _row_apply_plot_cigar(row, axes, plot_introns=True):
    index = row.name
    cigar = row.cigar
    gen_start = row.original_chr_pos
    is_adapted = row.t5
    polya_length = row.polya_length
    strand = row.strand

    if is_adapted == '-':
        color = 'black'
    else:
        color = 'red'
    return _add_patches_from_cigars_and_gen_pos(axes, cigar, gen_start, index, strand,
                                                color=color,
                                                plot_introns=plot_introns,
                                                tail_length=polya_length)


def _get_gene_coordinates(gene_id=None, gene_name=None,
                          parsed_gtf_path="/data16/marcus/genomes/elegansRelease100"
                                          "/Caenorhabditis_elegans.WBcel235.100.gtf.parquet",
                          ) -> (str, str, int, int):
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
                                                   "end",
                                                   "strand"]].query("feature == 'gene'")
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
    strand = entry_of_interest["strand"]
    print(f"Found entry for {gene_name} ({gene_id}) on chromosome "
          f"{chromosome:>5} at ({start}, {end}) on the '{strand}' strand")
    return gene_name, chromosome, start, end, strand


def _add_to_main_array_for_each_read(cigars_and_genomic_starts, chr_length, count_Ds_as_maps=False) -> np.array:
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
            elif code in maps:  # TODO: the D "belongs" above, but not yet...
                coverage_array[genomic_pos:genomic_pos + length] += 1
                genomic_pos += length
    return coverage_array


def _run_coverage_calc(bam_df, chrs=("I", "II", "III", "IV", "V", "X", "MtDNA"), count_Ds_as_maps=False):
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
                           'X': 17718726,
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
        coverage_array = _add_to_main_array_for_each_read(cigars_and_genomic_start_positions, chr_length,
                                                          count_Ds_as_maps=count_Ds_as_maps)
        array_dict[chromosome] = coverage_array
    return array_dict


def coverage_plotting_5tera(bam_df_for_plot, gene_name=None, gene_id=None,
                            save_dir=None, save_suffix=None, count_Ds_as_maps=False,
                            rpm_normalize=False, provide_axes: (plt.Axes, plt.Axes) = (None, None)):
    if gene_name:
        _, gene_chr, start, end, strand = _get_gene_coordinates(gene_name=gene_name)
    elif gene_id:
        gene_name, gene_chr, start, end, strand = _get_gene_coordinates(gene_id=gene_id)
    else:
        raise LookupError(f"Please provide either the gene_name or gene_id for the target gene to plot!")
    chr_array_ad = _run_coverage_calc(bam_df_for_plot.query("t5 == '+'"),
                                      chrs=gene_chr, count_Ds_as_maps=count_Ds_as_maps)[gene_chr]
    chr_array_unad = _run_coverage_calc(bam_df_for_plot.query("t5 == '-'"),
                                        chrs=gene_chr, count_Ds_as_maps=count_Ds_as_maps)[gene_chr]
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

    if isinstance(provide_axes[0], plt.Axes) and isinstance(provide_axes[1], plt.Axes):
        axes: (plt.Axes, plt.Axes) = provide_axes
    else:
        fig, axes = plt.subplots(2, 1,
                                 figsize=(8, 4),
                                 sharex='all',
                                 gridspec_kw={"height_ratios": [1, 4]}
                                 )
    
    # Pycharm is being annoying and throwing warning about None not have Axes methods, maybe this will help?
    axes[0]: plt.Axes
    axes[1]: plt.Axes
    
    axes[0].fill_between(index_array, zeros_array, locus_array_ad, color='red')
    axes[1].fill_between(index_array, zeros_array, locus_array_unad, color='black')
    axes[0].set_xticks([])
    axes[1].set_xticks([])
    # axes[0].set_xlabel
    
    if isinstance(provide_axes[0], plt.Axes) and isinstance(provide_axes[1], plt.Axes):
        pass
    else:
        plt.tight_layout()
        plt.style.use("seaborn-paper")
        if isinstance(save_dir, str):
            save_path = f"{save_dir}/{get_dt(for_file=True)}_readCoveragePlotting_{gene_name}{save_suffix}"
            print(f"Saving plot to {save_path} + .png/.svg...")
            plt.savefig(save_path + ".svg")
            plt.savefig(save_path + ".png")
        plt.show()


def plot_reads(reads_df, gene_id_to_plot=None, gene_name_to_plot=None,
               save_dir=None, save_suffix="", plot_width_and_height=(25, 5),
               subsample_fraction=None, subsample_number=None,
               t5_pos_count=None, t5_neg_count=None,
               pad_x_axis_bounds_by=None, only_keep_reads_matched_to_gene=True,
               provided_axis=None):
    gene_name, chromosome, genomic_start, genomic_end, gene_strand = _get_gene_coordinates(gene_name=gene_name_to_plot,
                                                                                           gene_id=gene_id_to_plot)

    if isinstance(subsample_fraction, float):
        subsampled_reads_df = reads_df.sample(frac=subsample_fraction)
    elif isinstance(subsample_number, int):
        subsampled_reads_df = reads_df.sample(n=subsample_number)
    else:
        subsampled_reads_df = reads_df  # Just to have the same variable name!
    if only_keep_reads_matched_to_gene:
        all_gene_reads = subsampled_reads_df.query(f"gene_name == '{gene_name}'")
    else:
        raise NotImplementedError(f"This doesn't currently work...")
    gene_df_t5_pos = all_gene_reads.query("t5 == '+'")
    if isinstance(t5_pos_count, int):
        try:
            gene_df_t5_pos = gene_df_t5_pos.sample(t5_pos_count)
        except ValueError:
            print(f"Failed to subselect {t5_pos_count} adapted reads because only {gene_df_t5_pos.shape[0]} exist."
                  f"\nWe are just going to plot all that exist")
    gene_df_t5_neg = all_gene_reads.query("t5 == '-'")
    if isinstance(t5_neg_count, int):
        try:
            gene_df_t5_neg = gene_df_t5_neg.sample(t5_neg_count)
        except ValueError:
            print(f"Failed to subselect {t5_neg_count} adapted reads because only {gene_df_t5_neg.shape[0]} exist."
                  f"\nWe are just going to plot all that exist")
    gene_df = pd.concat([gene_df_t5_pos, gene_df_t5_neg])

    # Now for the actual plotting!!
    plt.style.use('default')
    if isinstance(provided_axis, plt.Axes):
        ax = provided_axis
    else:
        fig, ax = plt.subplots(figsize=plot_width_and_height)

    if gene_strand == "-":
        sort_order = ["t5", "chr_pos", "original_chr_pos", "read_length"]
        sort_order_ascending = [False, True, False, False]
    else:  # gene_strand == "+":
        sort_order = ["t5", "original_chr_pos", "chr_pos", "read_length"]
        sort_order_ascending = [False, False, False, False]
    tqdm.pandas(desc="Plotting Reads...")
    gene_df = gene_df.sort_values(sort_order, ascending=sort_order_ascending).reset_index(drop=True)
    gene_df.progress_apply(lambda row: _row_apply_plot_cigar(row, ax), axis=1)

    number_of_plotted_reads = gene_df.shape[0]
    ax.set_ylim(-1, number_of_plotted_reads + 1)

    if isinstance(pad_x_axis_bounds_by, int):
        ax.set_xlim(genomic_start - pad_x_axis_bounds_by,
                    genomic_end + pad_x_axis_bounds_by)
    else:
        ax.set_xlim(genomic_start, genomic_end)

    ax.set_xticks([])
    ax.set_yticks([])
    if isinstance(provided_axis, plt.Axes):
        pass
    else:
        if isinstance(save_dir, str):
            save_path = f"{save_dir}/{get_dt(for_file=True)}_readPlotting_{gene_name}{save_suffix}"
            print(f"Saving plot to {save_path} + .png/.svg...")
            plt.savefig(save_path + ".svg")
            plt.savefig(save_path + ".png")
    return gene_df


if __name__ == '__main__':
    read_df, compressed_df = load_and_merge_lib_parquets(["xrn-1-5tera", "xrn-1-5tera-smg-6"],
                                                          drop_sub_n=1, add_tail_groupings=False,
                                                          drop_failed_polya=False, group_by_t5=True)
    gene_to_plot = "xbp-1"

    plot_reads(read_df.query("lib == 'xrn-1-5tera-smg-6'"), gene_name_to_plot=gene_to_plot,
               # t5_pos_count=1, t5_neg_count=30,
               pad_x_axis_bounds_by=100, save_dir=f"./outputDir", save_suffix="_smg-6-KO_allReads")
    plot_reads(read_df.query("lib == 'xrn-1-5tera'"), gene_name_to_plot=gene_to_plot,
               # t5_pos_count=20, t5_neg_count=20,
               pad_x_axis_bounds_by=100, save_dir=f"./outputDir", save_suffix="_WT_allReads")

    coverage_plotting_5tera(reads_df.query("lib == 'xrn-1-5tera'"), gene_name=gene_to_plot,
                            count_Ds_as_maps=True, rpm_normalize=True)
    coverage_plotting_5tera(reads_df.query("lib == 'xrn-1-5tera-smg-6'"), gene_name=gene_to_plot,
                            count_Ds_as_maps=True, rpm_normalize=True)
