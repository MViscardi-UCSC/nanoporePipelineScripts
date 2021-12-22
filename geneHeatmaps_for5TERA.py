"""
geneHeatmaps_for5TERA.py
Marcus Viscardi,    December 18, 2021

This is going to be an adaptation of geneHeatmaps2.py, mainly with the goal to
specifically plot 5TERA adapted reads (these should be specifically 5'monoP
containing mRNAs!)
"""
import time

from nanoporePipelineCommon import assign_with_josh_method, pick_libs_return_paths_dict, gene_names_to_gene_ids

from geneHeatmaps2 import manual_cdf

import numpy as np
import pandas as pd

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)


def np_cdf_from_hits(hit_list, bounds):
    # With tiny lists, the numpy function is slower than the manual method from geneHeatmaps2.py!
    # As list get longer, it seems like the numpy one is faster!
    #   The flip seems to be somewhere around 1000,
    #   but below that they are still really comparable:
    #       for 1,000 items:
    #           Numpy:  ~3.06E-4 sec
    #           Manual: ~2.50E-4 sec
    #       for 10,000 items:
    #           Numpy:  ~1.07E-3 sec
    #           Manual: ~2.78E-3 sec
    histo, bins = np.histogram(np.clip(hit_list, bounds[0], bounds[1]), bins=range(bounds[0], bounds[1] + 2),
                               density=True)
    return np.cumsum(histo)


def compress_on_transcripts(merged_df, drop_sub):
    transcript_groupby = merged_df.groupby(["gene_id", "transcript_id", "t5"])
    transcripts_df = transcript_groupby["to_stop"].apply(list).to_frame(name="stop_distances")
    transcripts_df["transcript_hits"] = transcript_groupby["transcript_id"].apply(len). \
        to_frame(name="transcript_hits")
    transcripts_df = transcripts_df.reset_index()
    transcripts_df = transcripts_df[transcripts_df["transcript_hits"] >= drop_sub]
    gene_df = gene_names_to_gene_ids()
    transcripts_df = transcripts_df.merge(gene_df, on="gene_id", how="inner")
    transcripts_df["identifier"] = transcripts_df["gene_name"].astype(str) + " (" + \
                                   transcripts_df["transcript_id"].astype(str) \
                                   + ")" + " [" + transcripts_df["transcript_hits"].astype(str) + "]"
    return transcripts_df


def make_dataframes_for_heatmaps(lib):
    _, lib_path = list(pick_libs_return_paths_dict([lib],
                                                   file_midfix="mergedOnReads",
                                                   file_suffix="parquet").items())[0]
    df = pd.read_parquet(lib_path)
    df: pd.DataFrame = df.drop(columns=["gene_id", "strand", "gene_name", "read_length",  # <- stuff from prev call w/JM
                                        'strand_fromFeatureCounts', 'qc_tag_featc', 'qc_pass_featc',  # <- featC stuff
                                        'gene_id_fromFeatureCounts', 'gene_name_fromFeatureCounts'])
    df = df.dropna(axis=0, subset=['leader_start',
                                   'adapter_start',
                                   'polya_start',
                                   'transcript_start',
                                   'read_rate',
                                   'polya_length',
                                   'qc_tag_polya'])
    df = assign_with_josh_method(df, "/data16/marcus/genomes/plus-pTRIxef_elegansRelease100",
                                 keepMultipleTranscriptInfo=True)
    # print(df)
    df_transcripts = compress_on_transcripts(df, 5)
    print(df_transcripts.query("t5 == '+'"))
    df_transcripts = df_transcripts.sort_values('transcriot_hits')
    return df, df_transcripts


def main(libs):
    df_dict = {}
    for lib in libs:
        df_dict[lib] = make_dataframes_for_heatmaps(lib)
    reads_df_dict = {lib: dfs[0] for lib, dfs in df_dict.items()}
    compressed_df_dict = {lib: dfs[1] for lib, dfs in df_dict.items()}
    print(compressed_df_dict)
    for lib, df in compressed_df_dict.items():
        # TODO: major error here due to stop_distances being lists of strings!! not numbers! WTF man...
        #       This is likely due to how the are currently stored in the damn readAssignment parquet!!!
        df['cdf'] = df.apply(lambda row: np_cdf_from_hits(row['stop_distances'], (-300, 300)), axis=1)
        print(df.cdf)


if __name__ == '__main__':
    main(["xrn-1-5tera", "xrn-1-5tera-smg-6"])
