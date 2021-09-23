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

pd.set_option("display.max_columns", None)

WORKING_DIR_DICT = {"polyA": "210528_NanoporeRun_0639_L3s",  # Best (overkill) depth
                    "riboD": "210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3",  # Gross
                    "totalRNA": "210709_NanoporeRun_totalRNA_0639_L3",  # Low depth
                    "polyA2": "210719_nanoporeRun_polyA_0639_L3_replicate",  # Good depth
                    "totalRNA2": "210720_nanoporeRun_totalRNA_0639_L3_replicate",  # Good depth
                    "xrn-1": "210905_nanoporeRun_totalRNA_5108_xrn-1-KD",  # First XRN-1 Knockdown
                    }


def load_reads_parquet(parquet_path) -> pd.DataFrame:
    return pd.read_parquet(parquet_path)


def process_to_pdf(compressed_reads_df, min_max=(-300, 300), smallest_allowed_utr=None):
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

    def transcripts_across_range(transcripts: list, min_x: int, max_x: int,
                                 smallest_allowed_utr=None) -> np.array:
        """
        So I am going to want to take transcript info from the list and use it to find
        each transcript's length relative to its stop
        :param transcripts: list of all transcript_ids, to filter the GTF annotations by
        :param min_x: minimum range value of the window
        :param max_x: maximum range value of the window
        :return: a array containing the distribution of transcripts spanning each point of the window
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
              f"that match this case.")
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

        annot_df[["stop_to_tss", "stop_to_tes"]] = pd.DataFrame(annot_df.apply(lambda x: calc_stop_to_tss_and_tes(row=x),
                                                                               axis=1, result_type='expand'),
                                                                index=annot_df.index)
        if isinstance(smallest_allowed_utr, int):
            annot_df = annot_df[annot_df["stop_to_tes"] >= smallest_allowed_utr]
        # Note: Transcripts in the stop_to_tes column with "0" as their value could be tossed,
        #       as they don't have annotated 3' UTRs!
        print(annot_df)
        stop_to_tss_pdf = manual_pdf([x * -1 for x in annot_df["stop_to_tss"].to_list()], min_x, max_x)
        stop_to_tes_pdf = manual_pdf(annot_df["stop_to_tes"].to_list(), min_x, max_x)
        # Convert the pointwize pdf to a bell curve type thing(?):
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
        import seaborn as sea
        import matplotlib.pyplot as plt
        quick_fig = sea.lineplot(x=range(min_x, max_x+1), y=spanning_transcript_count)
        # quick_fig.set(xlim=(-10, 200))
        plt.show()
        return spanning_transcript_count  # TODO: Something is going wrong here, and it is due to the pdf function!!

    print(compressed_reads_df.info())

    stop_distances = np.concatenate(compressed_reads_df.stop_distances.to_list()).ravel().tolist()
    stop_distances_pdf = manual_pdf(stop_distances, min_max[0], min_max[1],
                                    smallest_allowed_utr=smallest_allowed_utr)

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
    range_to_plot = (-1000, 1000)
    stops_pdf = process_to_pdf(load_reads_parquet(find_newest_matching_file(f"{working_dir}/output_dir"
                                                                            f"/merge_files/*_compressedOnTranscripts_"
                                                                            f"fromJoshsSystem.parquet")),
                               min_max=range_to_plot,
                               smallest_allowed_utr=50)
    plotly_pdf(stops_pdf, range_to_plot[0], range_to_plot[1])
