"""
General goal of this is to figure out what the mean read length for each gene is in libraries

This would be a nice metric to try and see if different mRNA selection methods are leading to
    an increase of degradation occurring
"""
import pandas as pd
import numpy as np


def load_merge_file(path):
    merge_df = pd.read_csv(path, sep="\t")
    merge_df["read_length"] = merge_df["sequence"].str.len()
    grouped_genes = merge_df.groupby("gene_id")
    gene_df = grouped_genes["read_id"].apply(len).to_frame(name="read_hits")
    gene_df["read_len_mean"] = grouped_genes["read_length"].apply(np.mean).to_frame(name="read_len_mean")
    gene_df["read_len_std"] = grouped_genes["read_length"].apply(np.std).to_frame(name="read_len_std")
    gene_df["read_lengths"] = grouped_genes["read_length"].apply(np.msort). \
        apply(list).to_frame(name="read_lengths")
    
    # gene_df["polya_lengths"] = grouped_genes["polya_length"].apply(np.msort). \
    #     apply(list).to_frame(name="polya_lengths")
    # gene_df["polya_mean"] = grouped_genes["polya_length"].apply(np.mean).to_frame(name="polya_mean")
    # gene_df["polya_stdev"] = grouped_genes["polya_length"].apply(np.std).to_frame(name="polya_stdev")
    
    gene_df = gene_df[gene_df["read_hits"] > 9]
    print(merge_df.columns)
    print(gene_df.sort_values(by="read_hits", ascending=False).head())


if __name__ == '__main__':
    print("Hi")
    pathdict = {
        "totalRNA": "/data16/marcus/working/210709_NanoporeRun_totalRNA_0639_L3/output_dir/merge_files/210709_01:38:16PM_mergedOnReads.tsv",
        "riboD": "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/output_dir/merge_files/210709_02:57:03PM__mergedOnReads.tsv",
        "polyA": "/data16/marcus/working/210528_NanoporeRun_0639_L3s/output_dir/merge_files/210601_05:05:27PM__mergedOnReads.tsv"}
    load_merge_file(pathdict["totalRNA"])
