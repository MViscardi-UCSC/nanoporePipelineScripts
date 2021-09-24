"""
geneMeta5Ends.py
Marcus Viscardi,    September 22, 2021

Going to make a script to plot 5' ends as a meta stop

Some additional details from talking with Josh:
1. Like in geneHeatmaps2.py, calculate the distance from read end to stop with josh's
2. I will want to normalize each position based on the number of genes at that position
    - This will help to avoid artificially creating a bell curve around the stop codon due
      to very small genes (rpl's and such)
3. I will also want to throw out genes with extremely short or unannotated 3' UTRs
    - This is similar to the issue with number 3, but is further hurt by bad annotations
"""
import pandas as pd
import numpy as np
from nanoporePipelineCommon import find_newest_matching_file, load_ski_pelo_targets

pd.set_option("display.max_columns", None)
pd.options.mode.chained_assignment = None

WORKING_DIR_DICT = {"polyA": "210528_NanoporeRun_0639_L3s",  # Best (overkill) depth
                    "riboD": "210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3",  # Gross
                    "totalRNA": "210709_NanoporeRun_totalRNA_0639_L3",  # Low depth
                    "polyA2": "210719_nanoporeRun_polyA_0639_L3_replicate",  # Good depth
                    "totalRNA2": "210720_nanoporeRun_totalRNA_0639_L3_replicate",  # Good depth
                    "xrn-1": "210905_nanoporeRun_totalRNA_5108_xrn-1-KD",  # First XRN-1 Knockdown
                    }


def load_reads_parquet(parquet_path, target_list=[], target_column=None) -> pd.DataFrame:
    df = pd.read_parquet(parquet_path)
    if target_list and target_column in list(df.columns):
        selected_df = df[df[target_column].isin(target_list)]
        if selected_df.shape[0] >= 1:
            df = selected_df
        else:
            print(f"{target_list} has zero hits in the {target_column} column... Moving forward without filtering.")
    elif target_column not in list(df.columns):
        print(f"Column {target_column} not found, please use one of these: {list(df.columns)}")
    return df


def manual_pdf(_stop_distances: list, min_x: int, max_x: int) -> np.array:
    man_pdf = []
    _stop_distances = np.flip(np.sort(_stop_distances))
    np_stop_dists = np.array(_stop_distances)
    for x in range(min_x, max_x + 1):
        hits_at_x = 0
        if x == min_x:
            hits_at_x = np.count_nonzero(np_stop_dists <= x)
        elif x == max_x:
            hits_at_x = np.count_nonzero(np_stop_dists >= x)
        else:
            hits_at_x = np.count_nonzero(np_stop_dists == x)
        man_pdf.append(hits_at_x)
    return man_pdf


def merge_tss_and_tes_to_cdf(stop_to_tes_pdf: list, stop_to_tss_pdf: list,
                             min_x: int, max_x: int, plot_it: bool = True) -> list:
    """
    So I am going to want to take transcript info from the list and use it to find
    each transcript's length relative to its stop
    
    :return: an array containing the distribution of transcripts spanning each point of the window
    """

    # Convert the pointwize pdf to a bell curve type thing(?),
    # Showing how many genes span each nucleotide along the window:
    spanning_transcript_count = []
    for i, (to_tss, to_tes) in enumerate(zip(stop_to_tss_pdf, stop_to_tes_pdf)):
        at_pos = to_tss - to_tes
        if i == 0:
            spanning_transcript_count.append(at_pos)
        else:
            sum_to_pos = spanning_transcript_count[-1]
            spanning_transcript_count.append(at_pos + sum_to_pos)
    print(list(zip(stop_to_tss_pdf, stop_to_tes_pdf)), spanning_transcript_count, sep="\n")
    if plot_it:
        import seaborn as sea
        import matplotlib.pyplot as plt
        quick_fig = sea.lineplot(x=range(min_x, max_x + 1), y=spanning_transcript_count)
        # quick_fig.set(xlim=(-10, 200))
        plt.show()
    return spanning_transcript_count


def annot_df_to_pdfs(annot_df: pd.DataFrame, max_x: int, min_x: int) -> [list, list]:
    """
    :param annot_df: annotation dataframe from the process_gtf method
    :param min_x: minimum range value of the window
    :param max_x: maximum range value of the window

    :return: 
    """
    stop_to_tss_pdf = manual_pdf([x * 1 for x in annot_df["stop_to_tss"].to_list()], min_x, max_x)
    stop_to_tes_pdf = manual_pdf([x * -1 for x in annot_df["stop_to_tes"].to_list()], min_x, max_x)
    return stop_to_tes_pdf, stop_to_tss_pdf


def process_gtf(transcripts: list, smallest_allowed_utr: int = None) -> pd.DataFrame:
    """
    Process a gtf file to look for the tss and tes distances from the stop codon for
    specified transcripts
    
    :param transcripts: list of all transcript_ids, to filter the GTF annotations by
    :param smallest_allowed_utr: int value of the shorter UTR that is accepted
    
    :return: a pandas dataframe with stop_to_tss and tes as columns
    """
    # First load the preprocessed (faster) GTF file:
    annot_df = pd.read_parquet("/data16/marcus/genomes/elegansRelease100/"
                               "Caenorhabditis_elegans.WBcel235.100.gtf.parquet")
    
    # Filter the GTF file for only transcripts that I found in my sequencing:
    annot_df = annot_df[annot_df["transcript_id"].isin(transcripts)]
    
    # Filter for only stop_codon and transcript feature rows:
    needed_columns = ["transcript_id", "strand", "start", "end"]
    stop_df = annot_df[annot_df.feature == "stop_codon"][needed_columns]
    transcript_df = annot_df[annot_df.feature == "transcript"][needed_columns]
    
    # Merge these two dataframes on the transcript ID
    annot_df = stop_df.merge(transcript_df, on=["transcript_id", "strand"],
                             suffixes=["_stop", "_transcript"])
    # Note: There seems to be a few transcript identities with 2 annotated stop codons.
    #       For now, the solution to this is just going to be dropping these weird cases,
    #       eventually I could use the annotated 3' UTR to pick the "real" stop.
    #       Example:    Gene: WBGene00018161;   Transcript: F38A5.2b.2
    
    print(f"Going to drop transcript_ids with multiple annotated stops or transcript regions"
          f"\nThere are {annot_df.duplicated(subset='transcript_id').shape[0]} transcript_ids"
          f" that match this case.")
    annot_df = annot_df.drop_duplicates(subset="transcript_id").sort_values("transcript_id", ignore_index=True)

    def calc_stop_to_tss_and_tes(row) -> [int, int]:
        # Split the row into variables:
        strand = row["strand"]
        start_stop = row["start_stop"]
        end_stop = row["end_stop"]
        start_trans = row["start_transcript"]
        end_trans = row["end_transcript"]

        # Calculate the distances from the stop codon to either end of the transcript
        if strand == "-":
            stop_to_tss = end_trans - end_stop
            stop_to_tes = start_stop - start_trans
        elif strand == "+":
            stop_to_tss = start_stop - start_trans
            stop_to_tes = end_trans - end_stop
        else:
            raise NotImplementedError(f"The strand for one gene was \"{strand}\"!?!?")
        return stop_to_tss, stop_to_tes

    annot_df[["stop_to_tss", "stop_to_tes"]] = pd.DataFrame(
        annot_df.apply(lambda x: calc_stop_to_tss_and_tes(row=x),
                       axis=1, result_type='expand'),
        index=annot_df.index)
    if isinstance(smallest_allowed_utr, int):
        annot_df = annot_df[annot_df["stop_to_tes"] >= smallest_allowed_utr]
    # Note: Transcripts in the stop_to_tes column with "0" as their value could be tossed,
    #       as they don't have annotated 3' UTRs!
    print(annot_df)
    return annot_df


def normalize_list_internally(un_normalized_list: list, norm_type="sum_to_max", normed_max=1.0):
    if norm_type == "sum_to_max":
        norm_factor = np.sum(un_normalized_list)
    elif norm_type == "max_value_is_max":
        norm_factor = np.max(np.abs(un_normalized_list))
    else:
        raise NotImplementedError(f"norm_type passed ({norm_type}) is not implemented at this time!")
    relative_norm_factor = normed_max / norm_factor
    normalized_array = [item * relative_norm_factor for item in un_normalized_list]
    return normalized_array


def process_to_pdf(compressed_reads_df: pd.DataFrame, min_max: (int, int) = (-300, 300),
                   smallest_allowed_utr: int = None) -> list:
    """
    This method takes the reads dataframe, and converts it to a normalized pdf over the min_max
    window
    
    :param compressed_reads_df: dataframe that was produced by the compress_on_transcripts method
           from geneHeatmaps2 (pass the save_as param!)
    :param min_max: the range over which the pdf will be produced
    :param smallest_allowed_utr: lower limit of annotated 3' UTR length that will be allowed through
    :return: a pdf across the min_max window of normalized number of 5'ends per position relative
             to the stop codon
    """
    print(compressed_reads_df.info())
    min_x, max_x = min_max
    # This will figure out how many potential genes there are at each location,
    #   which I can later use to normalize the stop_distances PDF (Hopefully. . . ?)
    annotation_df = process_gtf(compressed_reads_df.transcript_id.to_list(),
                                smallest_allowed_utr=smallest_allowed_utr)
    if isinstance(smallest_allowed_utr, int):
        passed_transcript_list = annotation_df["transcript_id"].to_list()
        compressed_reads_df = compressed_reads_df[compressed_reads_df.transcript_id.isin(passed_transcript_list)]

    # Turn each individual transcript into a pdf of "window size"
    compressed_reads_df["stop_distances_pdf"] = pd.DataFrame(
        compressed_reads_df.apply(lambda x: manual_pdf(x["stop_distances"],
                                                       min_x,
                                                       max_x),
                                  axis=1), index=compressed_reads_df.index)

    # Normalize each pdf so that their all only relative to their own transcript
    compressed_reads_df["normed_stops_pdf"] = pd.DataFrame(
        compressed_reads_df.apply(lambda x: normalize_list_internally(x["stop_distances_pdf"]),
                                  axis=1), index=compressed_reads_df.index)

    # Convert to a list of lists and compress all of these so that it's a single list of "window size"
    stop_distances = np.sum(compressed_reads_df.normed_stops_pdf.to_list(), axis=0)
    print(stop_distances)

    # Calculate normalization factor from annotations!
    stop_to_tes_pdf, stop_to_tss_pdf = annot_df_to_pdfs(annotation_df, max_x, min_x)
    normalization_factor_list = merge_tss_and_tes_to_cdf(stop_to_tss_pdf, stop_to_tes_pdf, min_x, max_x)

    # Normalize the stop_distances by the normalization factor list
    normalized_stop_distance_pdf = []
    for (pos, pos_txns) in zip(stop_distances, normalization_factor_list):
        if pos_txns != 0:
            normalized_stop_distance_pdf.append(pos / pos_txns)
            print(f"{pos:>6} / {pos_txns} = {pos / pos_txns}")
        else:
            normalized_stop_distance_pdf.append(pos / 1)
    return normalized_stop_distance_pdf


def plotly_pdf(stop_distances_pdf, window_min, window_max, mask_edges=True):
    import plotly.express as px
    print(stop_distances_pdf)
    if mask_edges:
        x = range(window_min, window_max + 1)[1:-1]
        y = stop_distances_pdf[1:-1]
    else:
        x = range(window_min, window_max + 1)
        y = stop_distances_pdf
    
    fig = px.line(x=x, y=y)
    
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
    title = "Meta Plot of 5' Ends of ONT dRNA-seq Reads"
    fig.update_layout(
        title=title,
        xaxis_title="Position along meta-transcript relative to STOP",
        yaxis_title="Relative fraction of reads",
        legend_title="Legend Title",
        font=dict(family="Courier New, monospace", size=18)
    )
    fig.show()


if __name__ == '__main__':
    working_dir = "/data16/marcus/working/" + WORKING_DIR_DICT["xrn-1"]
    range_to_plot = (-100, 100)
    stops_pdf = process_to_pdf(load_reads_parquet(find_newest_matching_file(f"{working_dir}/output_dir"
                                                                            f"/merge_files/*_compressedOnTranscripts_"
                                                                            f"fromJoshsSystem.parquet"
                                                                            ),
                                                  # target_list=["ets-4"],
                                                  # target_list=["E02C12.8"],
                                                  # target_list=load_ski_pelo_targets(as_df=False),
                                                  target_column="gene_id",
                                                  ),
                               min_max=range_to_plot, smallest_allowed_utr=50)
    plotly_pdf(stops_pdf, range_to_plot[0], range_to_plot[1])
