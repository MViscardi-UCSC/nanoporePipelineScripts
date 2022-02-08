"""
geneHeatmaps_for5TERA.py
Marcus Viscardi,    December 18, 2021

This is going to be an adaptation of geneHeatmaps2.py, mainly with the goal to
specifically plot 5TERA adapted reads (these should be specifically 5'monoP
containing mRNAs!)

Jan 31, 2022:   I could probably just use the actual gene's identified by featureCounts,
                then find where those genes start and end (in the GTF file). Those
                distances would work just fine for heatmaps, and would be independent
                of transcript identity. Eventually I could reassess per transcript for
                the interesting genes.
                    The main issue with an approach like this is that it lacks the option
                    to localize my heatmap's view window to be around the relevant stops...
                    The second issue is that these would be intron-independent...
                    which is crappy.
                Conclusion: For now, lets stick with the joshAssign method, while it'll be
                            a bit more noisy and have an asterix, I'll at least be able to
                            get it working!
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
                                 keepMultipleTranscriptInfo=True, add_names=True)
    print(df)
    df_transcripts = compress_on_transcripts(df, 5)
    print(df_transcripts.query("t5 == '+'"))
    df_transcripts = df_transcripts.sort_values('transcript_hits')
    return df, df_transcripts


def main(libs):
    from tqdm import tqdm
    from dashForClickAndViewPolyATails import load_and_merge_lib_parquets
    from time import sleep
    tqdm.pandas()
    
    # This has been a great way to assign and handle several libraries together,
    #   I am very proud of this 'step'.
    reads_df, compressed_df = load_and_merge_lib_parquets(libs,
                                                          genomeDir=f"/data16/marcus/genomes/"
                                                                    f"plus-pTRIxef_elegansRelease100",
                                                          drop_failed_polya=False, drop_sub_n=1,
                                                          keep_transcript_info=True,
                                                          read_pos_in_groupby=True, group_by_t5=True,
                                                          subsample_each_lib=20000)
    
    bounds = (-300, 300)
    print(f"\nProcessing stop_distances to CDFs in the range of {bounds} around the stop codon:")
    compressed_df['stop_cdf'] = compressed_df.progress_apply(
        lambda row: np_cdf_from_hits(row['stop_distances'],
                                     bounds),
        axis=1)
    
    compressed_df = compressed_df.query("chr_id != 'MtDNA'")
    
    for lib in libs:
        plotter_df = compressed_df.query(f"lib == '{lib}'")
        plotter_df = plotter_df.query("t5 == '+'").query("transcript_hits >= 50")
        
        # Grabbed below form geneHeatmaps2.py, this is a rough way to sort the individual genes/transcripts
        plotter_df["halfway_index"] = pd.DataFrame(plotter_df.apply(lambda x: np.where(x["stop_cdf"] >= 0.5)[0][0],
                                                                    axis=1), index=plotter_df.index)
        plotter_df["quarter_index"] = pd.DataFrame(plotter_df.apply(lambda x: np.where(x["stop_cdf"] >= 0.25)[0][0],
                                                                    axis=1), index=plotter_df.index)
        plotter_df["threequart_index"] = pd.DataFrame(plotter_df.apply(lambda x: np.where(x["stop_cdf"] >= 0.75)[0][0],
                                                                       axis=1), index=plotter_df.index)
        plotter_df.sort_values(by=[
            "quarter_index",
            "halfway_index",
            "threequart_index",
        ], inplace=True)
        
        # This will create a 2D array, of width = nucleotide window and length = # of genes/transcripts
        x_axis_cdfs = plotter_df["stop_cdf"].to_list()
    
        x_labels = list(range(bounds[0], bounds[1]+1))
    
        y_labels = plotter_df["gene_name"].to_list()
        from geneHeatmaps2 import plotly_imshow_heatmap
        plotly_imshow_heatmap(f"t5_heatmap_{lib}",
                              y_labels,
                              x_axis_cdfs,
                              x_labels,
                              extra_annotation=lib)
        
        # Attempt to plot w/ clustering?
        x_labels_plus = [f"{i}nt" for i in range(-300, 300+1)]
        new_df = plotter_df.copy(deep=True)
        new_df[x_labels_plus] = pd.DataFrame(plotter_df.stop_cdf.to_list(), index=plotter_df.index)
        from sklearn.cluster import AgglomerativeClustering
        clusters = new_df.shape[0]
        clusters_max = 25
        clusters = min([clusters, clusters_max])
        cluster = AgglomerativeClustering(n_clusters=clusters, affinity='euclidean', linkage='ward')
        new_df['cluster_order'] = cluster.fit_predict(new_df[x_labels_plus])
        new_df = new_df.sort_values('cluster_order')
        x_axis_cdfs = new_df["stop_cdf"].to_list()
        x_labels = list(range(bounds[0], bounds[1]+1))
        y_labels = new_df["gene_name"].to_list()
        from geneHeatmaps2 import plotly_imshow_heatmap
        plotly_imshow_heatmap(f"t5_heatmap_{lib}_clusteredTo{clusters}groups",
                              y_labels,
                              x_axis_cdfs,
                              x_labels,
                              extra_annotation=f"{lib} (w/ clustering - {clusters} groups)")
        
    sleep(2)
    print(compressed_df)


if __name__ == '__main__':
    main(["xrn-1-5tera", "xrn-1-5tera-smg-6"])
