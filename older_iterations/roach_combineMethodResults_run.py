#!/usr/bin/python3
"""
roach_combineMethodResults_run.py
Marcus Viscardi     Jan 16, 2021

This script will combine the results that were spit out of:
    1. guppy_basecaller
    2. Minimap2
    3. nanopolish polya
    4. featureCounts
Which will all (hopefully) be combined on read_ids. The bam files from Minimap2 will be the "baseline" which all other
    results will be mapped to. On Josh's recommendation I will try to drop duplicate reads based on their MAPQ scores.
"""

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sea

from datetime import datetime as dt

from unidip import UniDip

pd.set_option("display.max_columns", None)

# Hardcode the full path to our parent directory:
OUT_DIR = r"/data16/marcus/working/210528_NanoporeRun_0639_L3s/output_dir"
DATE_FOR_FILES = dt.now().strftime('%y%m%d')


def create_merge_df(out_dir, sam_file, feature_counts_file,
                    polya_file, print_info=False) -> pd.DataFrame:
    # First lets get the biggest one out of the way, importing the concatenated sam file:
    sam_df = pd.read_csv(f"{out_dir}/{sam_file}",
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
    featc_df = pd.read_csv(f"{out_dir}/{feature_counts_file}",
                           sep="\t", names=["read_id", "qc_tag", "qc_pass", "gene_id"])
    polya_df = pd.read_csv(f"{out_dir}/{polya_file}", sep="\t")
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
    if print_info:
        print("\n\n")
        print("#" * 100)
        print(f"\n\nMerged Dataframe info:")
        print(merge_df.info())
    with open(f"{out_dir}/{DATE_FOR_FILES}_merged_on_reads.tsv", "w") as merge_out_f:
        merge_df.to_csv(merge_out_f, sep="\t")
    return merge_df


def compress_on_genes(merged_df, drop_sub=None, dip_test=False, output_to_file=True,
                      additional_file_name=None) -> pd.DataFrame:
    # This creates a pandas "groupby" object that can be used to extract info compressed on gene_ids
    print("\nMerging information from Minimap2, featureCounts and Nanopolish-PolyA:")
    grouped_genes = merged_df.groupby("gene_id")
    gene_df = grouped_genes["read_id"].apply(len).to_frame(name="read_hits")
    gene_df["read_ids"] = grouped_genes["read_id"].apply(list).to_frame(name="read_ids")
    gene_df["polya_lengths"] = grouped_genes["polya_length"].apply(np.msort).apply(list).to_frame(name="polya_lengths")
    gene_df["polya_mean"] = grouped_genes["polya_length"].apply(np.mean).to_frame(name="polya_mean")
    if drop_sub:
        print(f"Dropping any genes with less than {drop_sub} read hits")
        gene_df = gene_df[gene_df["read_hits"] >= drop_sub]
        gene_df.sort_values("read_hits", ascending=False, inplace=True)
    if dip_test:
        print(f"Starting dip test at {dt.now().strftime('%H:%M:%S')}")
        gene_df["dip_test"] = grouped_genes["polya_length"].apply(quick_dip).to_frame(name="dip_test")
        gene_df["modality"] = gene_df["dip_test"].apply(len)
        print(f"Finished dip test at {dt.now().strftime('%H:%M:%S')}")
    print(f"Mean PolyA Length: {gene_df['polya_mean'].mean():.3f}")
    if additional_file_name:
        filename = f"{DATE_FOR_FILES}_compressed_on_genes{additional_file_name}.tsv"
    else:
        filename = f"{DATE_FOR_FILES}_compressed_on_genes.tsv"
    if output_to_file:
        with open(f"{OUT_DIR}/{filename}", "w")\
                as out_f:
            gene_df.to_csv(out_f, sep="\t")
    return gene_df


def plt_histo(list_to_plot, bin_size=1, print_bins=False, title=None, x_label=None, density=True,
              y_label=None, cdf_not_pdf=True, output_filename=None, output_format='svg'):
    bins_list = np.arange(0, max(list_to_plot) + bin_size, bin_size)
    if cdf_not_pdf:
        _, bins, pats = plt.hist(list_to_plot, bins=bins_list, density=density,
                                 histtype='step', cumulative=True, lw=2, range=(0, max(list_to_plot)))
        pats[0].set_xy(pats[0].get_xy()[:-1])
        legend_loc = "lower right"
    else:
        _, bins, pats = plt.hist(list_to_plot, bins=bins_list, density=density, range=(0, max(list_to_plot)))
        legend_loc = "upper right"
    if print_bins:
        print(bins)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    legend_elements = [Patch(facecolor='w', edgecolor='w', label=f'bin size: {bin_size}'),
                       Patch(facecolor='w', edgecolor='w', label=f'n = {len(list_to_plot)}')]
    plt.legend(handles=legend_elements, loc=legend_loc)
    if output_filename:
        plt.savefig(f"{OUT_DIR}/{output_filename}.{output_format}", format=output_format)
    plt.show()


def quick_dip(data_in, p_cutoff=0.05, num_trials=100, merge_distance=1):
    data_list = list(data_in)
    sorted_data = np.msort(data_list)
    indexes = UniDip(sorted_data, alpha=p_cutoff, ntrials=num_trials, mrg_dst=merge_distance).run()
    polya_length_groups = []
    for a, b in indexes:
        if b >= len(data_list):
            b = -1
        x, y = sorted_data[a], sorted_data[b]
        polya_length_groups.append((x, y))
    return polya_length_groups


def ski_pelo_enriched_comp(gene_df, filepath):
    ski_pelo_genes = pd.read_csv(filepath, names=["gene_id"])
    ski_pelo_genes = ski_pelo_genes.drop_duplicates(subset="gene_id", keep=False, ignore_index=True)
    ski_pelo_merge_df = ski_pelo_genes.merge(gene_df, how='inner',  on=["gene_id"])
    print(ski_pelo_merge_df.info())
    return ski_pelo_merge_df


def plot_increasing_samples(gene_df):
    
    def apply_mean_subsample(list_to_subsample, sample_size_out_of_100=0):
        # This function is mapped to all the polya_length lists in the dataframe
        list_len = len(list_to_subsample)
        samples_to_take = list_len*sample_size_out_of_100//100
        samples = np.random.choice(list_to_subsample, samples_to_take, replace=False)
        mean = np.mean(samples)
        return mean

    prop_cycle = plt.rcParams['axes.prop_cycle']
    color_options = prop_cycle.by_key()['color']
    legend_elements = [Patch(facecolor='w', edgecolor='w', label=f'bin size: {1}'),
                       Patch(facecolor='w', edgecolor='w', label=f'total input genes: {gene_df.shape[0]}')]
    
    for i, sample_percent in enumerate(range(20, 120, 20)):
        subsampled_means = gene_df["polya_lengths"].apply(apply_mean_subsample,
                                                          sample_size_out_of_100=sample_percent).to_list()
        print(subsampled_means)
        bins_list = np.arange(0, max(subsampled_means) + 1, 1)
        _, _, pats = plt.hist(subsampled_means, bins=bins_list, density=True,
                              histtype='step', cumulative=True, lw=2, color=color_options[i])
        pats[0].set_xy(pats[0].get_xy()[:-1])
        legend_elements.append(Patch(facecolor=color_options[i], edgecolor='k', label=f'Subsampled {sample_percent}%'))
    plt.legend(handles=legend_elements, loc="lower right")
    plt.title("Subsampling PolyA Lengths to Calc Means (>20 hits/gene)")
    plt.ylabel("Fraction of genes with \"X\" polyA tail length")
    plt.xlabel("Mean tail length of gene (nts)")
    plt.xlim([0, 150])
    plt.show()


def plot_some_top_hits(gene_df: pd.DataFrame, how_many_plots=10, title=None,
                       x_limit=None, bin_size=1, series_named=True):
    figure_size = (6.4, 4.8*how_many_plots/10)
    fig, axs = plt.subplots(how_many_plots, 1, sharex="all", figsize=figure_size)
    max_x_value = 0
    for i, ax in enumerate(axs):
        gene_series = gene_df.iloc[i]
        if series_named:
            gene_name = gene_series.name
        else:
            gene_name = gene_series.gene_id
        polya_lengths = gene_series.polya_lengths
        bins_list = np.arange(0, int(max(polya_lengths)) + bin_size, bin_size)
        ax.hist(polya_lengths, bins=bins_list, density=True)
        local_max = max(polya_lengths)
        if local_max > max_x_value:
            max_x_value = local_max
        ax.set_ylabel(f"{gene_name}\n(n={len(polya_lengths)})",
                      rotation=45, ha="right",
                      labelpad=0, linespacing=0.75)
        ax.set_yticklabels([])
        ax.set_yticks([[0.025]])
        ax.get_xaxis().set_visible(False)
    for ax in axs:
        if x_limit:
            ax.set_xlim(0, x_limit)
        else:
            ax.set_xlim(0, max_x_value)
    axs[how_many_plots-1].get_xaxis().set_visible(True)
    axs[how_many_plots-1].set_xlabel(f"PolyA Tail Length (nts)")
    freq_text = plt.figtext(0.02, 0.98, "Frequencies:", ha="left", va="bottom",
                            fontsize=10, fontweight='semibold')
    plt.figtext(0.02, 0.97, "(tick @ 0.025)", ha="left", fontsize=10)
    if title:
        plt.suptitle(title)
    else:
        plt.suptitle(f"Top {how_many_plots} Genes' Tail Lengths:")
    plt.tight_layout()
    plt.subplots_adjust(left=None,  # 0.25,
                        bottom=None,  # 0.15,
                        right=None,
                        top=None,
                        wspace=None,
                        hspace=0.1)
    plt.show(bbox_inches="tight")


if __name__ == '__main__':
    all_merged_df = create_merge_df(OUT_DIR, "cat_allReads.sam",
                                    "cat_featureCounts_assigned.tsv",
                                    "cat_polya_passed.tsv")

    genes_df = compress_on_genes(all_merged_df,
                                 drop_sub=20)
    
    # ski_pelo_enriched_comp(genes_df,
    #                        "/data16/marcus/working/210119_SkiPeloTargets_fromStarDust/"
    #                        "170723_MSandM.wtAndSkiPelo_Bounds_-12_-14_S.DESeqgeneCts_diff"
    #                        "Expression_2.7319418642771283e-06Down.txt")
    
    # plt_histo(genes_df["read_hits"].to_list(), x_label="Number of reads per gene",
    #           y_label="Fraction of genes with x reads mapped", cdf_not_pdf=True,
    #           output_filename=f"{DATE_FOR_FILES}_readsPerGene", output_format='png',
    #           title="Number of Reads Mapped to Genes\n(genes with >20 reads)")
    # plt_histo(genes_df["polya_mean"].to_list(), x_label="Mean tail length of gene (nts)",
    #           y_label="Fraction of genes with x polyA tail length", cdf_not_pdf=True,
    #           bin_size=1, output_filename=f"{DATE_FOR_FILES}_meanTailLenForGenes_noDrop",
    #           title="Mean Tail Lengths for Genes\n(including genes w/low read counts)", output_format='png')
    
    # x = genes_df["polya_mean"].to_list()
    # y = genes_df["read_hits"].to_list()
    # g = sea.jointplot(data=genes_df, x="polya_mean", y="read_hits")
    # #g.plot_joint(sea.scatterplot, color='k', marker='+')
    # g.plot_marginals(sea.histplot, kde=True, zorder=1)
    # g.set_axis_labels(xlabel="Mean PolyA Tail Length", ylabel="Number of Read Hits")
    # g.fig.suptitle("Read Counts and Tail Length of Genes (genes with >20 reads)")
    # g.fig.tight_layout()
    # g.fig.subplots_adjust(top=0.95)
    # g.ax_joint.set_yscale('log')
    # plt.show()
    
    # print(genes_df.head(), genes_df.columns, sep='\n')
    # print(genes_df.loc["WBGene00006713", "polya_lengths"])
    # print("\n\n")
    # plot_increasing_samples(genes_df)
    
    # plot_some_top_hits(genes_df, x_limit=150)
    
    print("---Finished---")
