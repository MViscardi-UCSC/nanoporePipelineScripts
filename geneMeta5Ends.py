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
from nanoporePipelineCommon import find_newest_matching_file

WORKING_DIR_DICT = {"polyA": "210528_NanoporeRun_0639_L3s",  # Best (overkill) depth
                    "riboD": "210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3",  # Gross
                    "totalRNA": "210709_NanoporeRun_totalRNA_0639_L3",  # Low depth
                    "polyA2": "210719_nanoporeRun_polyA_0639_L3_replicate",  # Good depth
                    "totalRNA2": "210720_nanoporeRun_totalRNA_0639_L3_replicate",  # Good depth
                    "xrn-1": "210905_nanoporeRun_totalRNA_5108_xrn-1-KD",  # First XRN-1 Knockdown
                    }


def load_reads_parquet(parquet_path) -> pd.DataFrame:
    return pd.read_parquet(parquet_path)


def process_to_pdf(compressed_reads_df, min_max=(-300, 300)):
    # This is going to be a SUPER hacky first pass at getting this working,
    #   basically ignoring the stuff I note at the top of this file!
    def manual_pdf(_stop_distances: list, min_x: int, max_x: int) -> np.array:
        man_pdf = []
        _stop_distances.sort(reverse=False)
        np_stop_dists = np.array(_stop_distances)
        # Min and Max could also just be the bounds that get fed in here.
        #   Anything outside those bounds would just get rounded in?
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

    def transcripts_across_range(transcripts: list, min_x: int, max_x: int) -> np.array:
        # So I am going to want to take transcript info from the list and use it to find
        #   each transcript's length relative to its stop
        annot_df = pd.read_parquet("/data16/marcus/genomes/elegansRelease100/"
                                   "Caenorhabditis_elegans.WBcel235.100.gtf.parquet")
        annot_df = annot_df[annot_df["transcript_id"].isin(transcripts)]

        needed_columns = ["start", "end", "strand", "transcript_id"]

        stop_df = annot_df[annot_df.feature == "stop_codon"][needed_columns]
        transcript_df = annot_df[annot_df.feature == "transcript"][needed_columns]
        annot_df = stop_df.merge(transcript_df, on=["transcript_id", "strand"],
                                 suffixes=["_stop", "_transcript"])

        def calc_stop_to_tss(strand, start_stop, end_stop, start_trans, end_trans):
            if strand == "-":
                stop_to_tss = end_trans - end_stop
            elif strand == "+":
                stop_to_tss = start_trans - start_stop
            else:
                raise NotImplementedError(f"The strand for one gene was \"{strand}\"!?!?")
            return stop_to_tss

        annot_df["stop_to_tss"] = pd.DataFrame(annot_df.apply(lambda x: calc_stop_to_tss(x["strand"],
                                                                                         x["start_stop"],
                                                                                         x["end_stop"],
                                                                                         x["start_transcript"],
                                                                                         x["end_transcript"]), axis=1),
                                               index=annot_df.index)
        print(annot_df)
        stop_to_tss_pdf = manual_pdf(annot_df["stop_to_tss"].to_list(), min_x, max_x)  # This is for sure the wrong way to go about this
        return stop_to_tss_pdf

    print(compressed_reads_df.info())

    stop_distances = np.concatenate(compressed_reads_df.stop_distances.to_list()).ravel().tolist()
    stop_distances_pdf = manual_pdf(stop_distances, min_max[0], min_max[1])

    # This will figure out how many potential genes there are at each location,
    #   which I can use to normalize the stop_distances PDF (Hopefully. . . ?)
    normalization_factor_list = transcripts_across_range(compressed_reads_df.transcript_id.to_list(),
                                                         min_max[0], min_max[1])
    return stop_distances_pdf


def plotly_pdf(stop_distances_pdf, window_min, window_max):
    import plotly.express as px
    fig = px.line(x=range(window_min, window_max + 1)[1:-1], y=stop_distances_pdf[1:-1])
    fig.show()


if __name__ == '__main__':
    working_dir = "/data16/marcus/working/" + WORKING_DIR_DICT["xrn-1"]
    range_to_plot = (-300, 300)
    stops_pdf = process_to_pdf(load_reads_parquet(find_newest_matching_file(f"{working_dir}/output_dir"
                                                                            f"/merge_files/*_compressedOnTranscripts_"
                                                                            f"fromJoshsSystem.parquet")),
                               min_max=range_to_plot)
    plotly_pdf(stops_pdf, range_to_plot[0], range_to_plot[1])
