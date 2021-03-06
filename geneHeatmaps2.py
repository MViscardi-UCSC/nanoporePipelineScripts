"""
geneHeatmaps2.py
Marcus Viscardi,    August 14, 2021

So I managed to overcomplicate everything while trying
    different methods in geneHeatmaps.py, so now I am
    going to try and take a simpler approach.

The plan is to import the mergedOnReads.tsv rather than
    any of the compressedOnGenes files. B/c the majority
    of the early processing steps to get distances from
    stops are on a per read basis: I'll use the read/row
    format of the mergedOnReads file to be able to
    iterate through things more cleanly/obviously.

CIGAR parsing will still be necessary to "flip" negative
    strand reads to get their 5' ends.

The other introduced change will be to entirely use the
    Caenorhabditis_elegans.WBcel235.100.allChrs.txt file
    which is produced by Josh's prepareReadAssignmentFile3.py
This annotation file is in the format:
    CHR_chrpos  WBGene:strand   iso1:posRelStart:posRelStop|iso2:...
The only striking issue is that Josh's dropped any locations with
    overlapping genes (which makes heaps of sense for short read
    sequencing, but may weaken my nanopore stuff).

After getting stop distance for each read, I'll be able to compress
    these on genes OR transcripts, calculate CDFs and feed them into
    the same plot function from geneHeatmaps.py (or I can rewrite
    that bit too. . . ).

Hype!
"""
from typing import List

import numpy as np
import pandas as pd
import regex as re

from nanoporePipelineCommon import find_newest_matching_file, get_dt, assign_with_josh_method,\
    gene_names_to_gene_ids, parse_read_assignment_allChrs_txt


def load_reads_parquet(read_parquet, head=None) -> pd.DataFrame:
    print(f"Loading merged on reads file from: {read_parquet} ", end="")
    if isinstance(head, int):
        df = pd.read_parquet(read_parquet, nrows=head)
    else:
        df = pd.read_parquet(read_parquet)
    print(". ", end="")
    if "original_chr_pos" not in df.columns.to_list():
        print(f"\nMaking adjustments for 5' ends")
        df["original_chr_pos"] = df["chr_pos"]
        df["chr_pos"] = df.apply(lambda read: _flip_neg_strand_genes(read["original_chr_pos"],
                                                                     read["cigar"],
                                                                     read["strand"]),
                                 axis=1)
        df.to_parquet(read_parquet)
        print("(saved parquet w/ adjusted read_ends) ", end="")
    print(". ")
    return df


def _flip_neg_strand_genes(chr_position: int, cigar: str, strand: str) -> int:
    if strand == "+":
        read_end = chr_position
        return read_end
    else:
        numbers = list(map(int, re.findall(rf'(\d+)[MDNSI]', cigar)))
        cigar_chars = re.findall(rf'\d+([MDNSI])', cigar)
        mnd_nums, mnd_chars = [], []
        for i, cigar_char in enumerate(cigar_chars):
            if cigar_char in "MND":
                mnd_chars.append(cigar_char)
                mnd_nums.append(numbers[i])
        read_end = chr_position + sum(mnd_nums)
        return read_end


def plotly_cdf(cdf_values, x_values):
    import plotly.graph_objects as go

    fig = go.Figure()
    fig.add_scatter(x=x_values, y=cdf_values, line_shape='hv')
    fig.add_shape(type="line",
                  x0=0, x1=0, xref='x',
                  y0=0, y1=1, yref='paper',
                  line_color='darkred')
    fig.show()


def plotly_imshow_heatmap(output_name, y_labels, x_axis_cdfs, x_labels,
                          extra_annotation=None, title=None):
    import plotly.express as px
    fig = px.imshow(x_axis_cdfs,
                    x=x_labels,
                    y=y_labels,
                    aspect="free",
                    color_continuous_scale='RdBu_r')
    fig.add_shape(type="line",
                  x0=0, x1=0, xref='x',
                  y0=0, y1=1, yref='paper',
                  line_color='darkred')
    fig.add_annotation(x=0, y=-0.05, xref="x", yref="paper",
                       text="<b>Stop Codon</b>", showarrow=False,
                       font_size=15, )
    fig.add_annotation(x=0, y=-0.05, xref="paper", yref="paper",
                       text="<b>CDS</b>", showarrow=False,
                       font_size=15)
    fig.add_annotation(x=1, y=-0.05, xref="paper", yref="paper",
                       text="<b>3' UTR</b>", showarrow=False,
                       font_size=15)
    if isinstance(title, str):
        fig.update_layout(title=title)
    if isinstance(extra_annotation, str):
        fig.add_annotation(x=0, y=1.03, xref="paper", yref="paper",
                           text=f"<b>Library Working Directory:</b> {extra_annotation}", showarrow=False,
                           font_size=12)
    fig.update_traces(hovertemplate="<b>Gene Name: %{y}</b><br>Dist from Stop: %{x:+}nt"
                                    "<br>Cumulative Fraction of 5' Ends here: %{z:.3%}"
                                    "<extra></extra>", )
    fig.update_layout(coloraxis_colorbar=dict(tickformat=",.1%",
                                              title="CDF<br>of 5' Ends,<br>per Isoform",
                                              titleside="bottom"
                                              ),
                      hovermode="y")
    fig.write_html(f"./testOutputs/{get_dt(for_file=True)}_{output_name}.html")
    fig.show()


def manual_cdf(stop_distances: list, min_x: int, max_x: int) -> np.array:
    man_cdf = [0]
    stop_distances.sort(reverse=False)
    np_stop_dists = np.array(stop_distances)
    # Min and Max could also just be the bounds that get fed in here.
    #   Anything outside those bounds would just get rounded in?
    for x in range(min_x, max_x + 1):
        if x == min_x:  # If first number, sum all below
            hits_at_x = np.count_nonzero(np_stop_dists <= x)
        elif x == max_x:  # If last number, sum all above
            hits_at_x = np.count_nonzero(np_stop_dists >= x)
        else:  # Otherwise, store sum of values at that point
            hits_at_x = np.count_nonzero(np_stop_dists == x)
        sum_to_x = man_cdf[-1]
        sum_w_x = sum_to_x + hits_at_x
        man_cdf.append(sum_w_x)
    norm_sum = man_cdf[-1]
    normalized_cdf = np.array([x / norm_sum for x in man_cdf])
    return normalized_cdf


def plot_heatmap(smallish_df: pd.DataFrame, bounds: List[int] = (-300, 300),
                 title=None, extra_annotation=None, output_name="cdfHeatmap"):

    print("\n\n")
    # I hope I can use these max distances to act as max bounds for my heatmaps
    #   I'd introduce 0s anywhere I don't get hits
    smallish_df['max_dist'] = smallish_df["stop_distances"].apply(max)
    smallish_df['min_dist'] = smallish_df["stop_distances"].apply(min)

    # Because of some massive outliers, this min-max setup doesn't really work!!
    #   Since this I just changed to passing bounds!
    total_max = smallish_df['max_dist'].max()
    total_min = smallish_df['min_dist'].min()

    print("Max distance from stop codon:", total_max)
    print("Min distance from stop codon:", total_min)
    plotter_df = pd.DataFrame()
    plotter_df["gene_name"] = smallish_df["identifier"]
    plotter_df["norm_cdf"] = smallish_df.apply(lambda x: manual_cdf(x["stop_distances"],
                                                                    bounds[0],
                                                                    bounds[1]),
                                               axis=1)
    plotter_df["halfway_index"] = pd.DataFrame(plotter_df.apply(lambda x: np.where(x["norm_cdf"] >= 0.5)[0][0],
                                                                axis=1), index=plotter_df.index)
    plotter_df["quarter_index"] = pd.DataFrame(plotter_df.apply(lambda x: np.where(x["norm_cdf"] >= 0.25)[0][0],
                                                                axis=1), index=plotter_df.index)
    plotter_df["threequart_index"] = pd.DataFrame(plotter_df.apply(lambda x: np.where(x["norm_cdf"] >= 0.75)[0][0],
                                                                   axis=1), index=plotter_df.index)
    plotter_df.sort_values(by=[
        "quarter_index",
        "halfway_index",
        "threequart_index",
    ], inplace=True)
    plotter_df["norm_cdf_lists"] = pd.DataFrame(plotter_df.apply(lambda x: x["norm_cdf"].tolist(),
                                                                 axis=1), index=plotter_df.index)

    x_axis_cdfs = plotter_df["norm_cdf"].to_list()

    x_labels = list(range(bounds[0] - 1, bounds[1] + 1))

    y_labels = plotter_df["gene_name"].to_list()

    if plotter_df.shape[0] > 1:
        plotly_imshow_heatmap(output_name,
                              y_labels,
                              x_axis_cdfs,
                              x_labels,
                              title=title,
                              extra_annotation=extra_annotation)
    elif plotter_df.shape[0] == 1:
        plotly_cdf(x_axis_cdfs[0], x_labels)
    else:
        raise ValueError(f"The dataframe for plotting is empty!!")


def compress_on_transcripts(merged_df, drop_sub, save_as=None):
    transcript_groupby = merged_df.groupby(["gene_id_fromAssign", "transcript_id"])
    transcripts_df = transcript_groupby["to_stop"].apply(list).to_frame(name="stop_distances")
    transcripts_df["transcript_hits"] = transcript_groupby["transcript_id"].apply(len).to_frame(name="transcript_hits")
    transcripts_df = transcripts_df.reset_index().rename(columns={"gene_id_fromAssign": "gene_id"})
    transcripts_df = transcripts_df[transcripts_df["transcript_hits"] >= drop_sub]
    gene_df = gene_names_to_gene_ids()
    transcripts_df = transcripts_df.merge(gene_df, on="gene_id", how="inner")
    transcripts_df["identifier"] = transcripts_df["gene_name"].astype(str) + " (" + \
                                   transcripts_df["transcript_id"].astype(str) \
                                   + ")" + " [" + transcripts_df["transcript_hits"].astype(str) + "]"
    if isinstance(save_as, str):
        transcripts_df.to_parquet(save_as)
        print(f"Saved file to {save_as}")
    return transcripts_df


def parse_annies_deseq(csv_path: str) -> list:
    df = pd.read_csv(csv_path, names=["gene_id",
                                      "wt_1", "smg-1_1", "smg-6_1",
                                      "wt_2", "smg-1_2", "smg-6_2"],
                     skiprows=1)
    df[["gene_id", "strand"]] = df["gene_id"].str.split(":", expand=True)
    return df["gene_id"].tolist()


def main(library_strs, genome_dir="/data16/marcus/genomes/elegansRelease100", drop_sub=10,
         target_list=[], target_column="gene_name"):
    for lib in library_strs:
        merged_df, working_dir = load_tsv_and_assign_w_josh_method(lib, genome_dir=genome_dir)
        smaller_df = compress_on_transcripts(merged_df, drop_sub,
                                             save_as=f"{working_dir}/output_dir/merge_files/"
                                                     f"{get_dt(for_file=True)}_compressedOnTranscripts_"
                                                     f"fromJoshsSystem.parquet")
        if target_list:
            selected_df = smaller_df[smaller_df[target_column].isin(target_list)]
            if selected_df.shape[0] >= 1:
                smaller_df = selected_df
            else:
                print(f"{target_list} has zero hits in the {target_column} column... No filtering run")
        print(merged_df, smaller_df)
        plot_heatmap(smaller_df)


def load_tsv_and_assign_w_josh_method(library_str, genome_dir="/data16/marcus/genomes/elegansRelease100",
                                      nrows=None) -> (pd.DataFrame, str):
    """
    A method to load files and merge based on the chromosome position to find matching genes
    
    :param library_str: Must be one of the following: polyA, riboD, totalRNA, polyA2, totalRNA2, xrn-1
    :param genome_dir: Path to genome dir w/ processed parquet file
    :param nrows: Numder of rows to parse from the tsv file
    
    :return: A dataframe of passed genes & the path to the working directory
    """
    working_dir_dict = {"polyA": "210528_NanoporeRun_0639_L3s",  # Best (overkill) depth
                        "riboD": "210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3",  # Gross
                        "totalRNA": "210709_NanoporeRun_totalRNA_0639_L3",  # Low depth
                        "polyA2": "210719_nanoporeRun_polyA_0639_L3_replicate",  # Good depth
                        "totalRNA2": "210720_nanoporeRun_totalRNA_0639_L3_replicate",  # Good depth
                        "xrn-1": "210905_nanoporeRun_totalRNA_5108_xrn-1-KD",  # First XRN-1 Knockdown
                        }
    try:
        working_dir = f"/data16/marcus/working/{working_dir_dict[library_str]}"
    except KeyError:
        raise KeyError(f"Key: {library_str} isn't in the working directory dictionary keys : {working_dir_dict.keys()}")
    reads_df = load_reads_parquet(find_newest_matching_file(f"{working_dir}/output_dir/merge_files/"
                                                            f"*_mergedOnReads.parquet"),
                                  head=nrows)
    merged_df = assign_with_josh_method(reads_df, genome_dir, keepMultipleTranscriptInfo=True)
    return merged_df, working_dir


if __name__ == '__main__':
    annie_deseq_csv_path = "/data16/marcus/working/210204_smg-1and6_alteredGenes_fromAnnie/" \
                           "210209_GenesUpInSMG1AndSMG6.csv"
    # annies_deseq_genes = parse_annies_deseq(annie_deseq_csv_path)

    main(["polyA", "totalRNA", "polyA2", "totalRNA2"],
         drop_sub=10,
         # target_list=annies_deseq_genes,
         # target_column="gene_id",
         # target_list=['tbb-2',
         #              'tbb-1',
         #              'mec-7',
         #              'mec-12'],
         # target_column='gene_name'
         )
