"""
interestingGeneTailLengths.py
Marcus Viscardi, August 04, 2021

Script to plot some specific gene tail lengths!

Main goal is to see if my runs look like figure 2E:
    https://pubmed.ncbi.nlm.nih.gov/29106412/#&gid=article-figures&pid=figure-2-uid-1
"""

import pandas as pd
import seaborn as sea
import matplotlib.pyplot as plt
from older_iterations.roach_compareReplicates_run import load_comp_on_genes
from step0_nanopore_pipeline import find_newest_matching_file, gene_names_to_gene_ids


def make_split_violins2(df1, df2, gene_list):
    # df1, df2 = load_comp_on_genes(f"{OUTDIR1}/210124_compressed_on_genes.tsv"), \
    #            load_comp_on_genes("/data16/marcus/prefix/210709_NanoporeRun_totalRNA_0639_L3/output_dir/merge_files/210709_01:38:18PM_compressedOnGenes.tsv")
    #            # load_comp_on_genes(f"{OUTDIR2}/210128_compressed_on_genes.tsv")
    print(">> Dataframes loaded into memory.")
    new_dfs = []
    for i, df in enumerate((df1, df2)):
        i += 1
        # df["hits_rank"] = df['read_hits'].rank(ascending=False)
        # # TODO: drop the below line as it shrinks the dataset a bit to make it easier to manage:
        # df = df[df["hits_rank"] <= 15]
        df = df[df["gene_name"].isin(gene_list)]
        df_tails = df["polya_lengths"].apply(pd.Series).merge(df["gene_id"], right_index=True, left_index=True) \
            .melt(id_vars=["gene_id"], value_name="tail_length").dropna().drop("variable", axis=1)
        print(f">> Dataframe {i}'s tail lengths flattened")
        df_reads = df["read_ids"].apply(pd.Series).merge(df["gene_id"], right_index=True, left_index=True) \
            .melt(id_vars=["gene_id"], value_name="read_id").dropna().drop("variable", axis=1)
        print(f">> Dataframe {i}'s read ids flattened")
        df_merge = df_reads.merge(df_tails)
        df_merge = df_merge.merge(df[["gene_id", "polya_mean", "read_hits", "gene_name"]])
        new_dfs.append(df_merge)
        print(f">> Dataframe {i} finished")
    df1_long, df2_long = new_dfs
    df1_long["replicate"] = "Replicate 1"
    df2_long["replicate"] = "Replicate 2"
    super_df = pd.concat([df1_long, df2_long], ignore_index=True)
    print(super_df.info())

    # # PLOT IT:
    # top_2 = super_df[super_df["hits_rank"] <= 2]
    sea.set_theme(style="whitegrid")
    # ax = sea.swarmplot(x="tail_length", y="gene_id", hue="replicate",
    #                    data=super_df, dodge=True, alpha=.01, zorder=1)
    ax = sea.violinplot(x="gene_name", y="tail_length", hue="replicate",
                        data=super_df, split=True, bw=0.2, inner="quartiles")
    plt.show()
    print(" . . . Done!")


if __name__ == '__main__':
    prefix = "/data16/marcus/working"
    working_dir_dict = {"polyA": "210528_NanoporeRun_0639_L3s",  # Best (overkill) depth
                        "riboD": "210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3",  # Gross
                        "totalRNA": "210709_NanoporeRun_totalRNA_0639_L3",  # Low depth
                        "polyA2": "210719_nanoporeRun_polyA_0639_L3_replicate",  # Good depth
                        "totalRNA2": "210720_nanoporeRun_totalRNA_0639_L3_replicate",  # Good depth
                        }
    suffix = "output_dir/merge_files/*_compressedOnGenes.tsv"
    # Paths to compressedOnGenes.tsv:
    file1 = find_newest_matching_file(f"{prefix}/{working_dir_dict['totalRNA2']}/{suffix}")
    file2 = find_newest_matching_file(f"{prefix}/{working_dir_dict['polyA2']}/{suffix}")

    path_to_parsed_gtf = "/data16/marcus/scripts/nanoporePipelineScripts/" \
                         "Caenorhabditis_elegans.WBcel235.100.gtf.dataframe_parse.tsv"
    df = gene_names_to_gene_ids()
    df_list = [load_comp_on_genes(file1), load_comp_on_genes(file2)]
    df_list = [df_i.merge(df, on="gene_id") for df_i in df_list]
    make_split_violins2(df_list[0], df_list[1], ["rpl-21",
                                                 "daf-21",  # I think this has another name
                                                 "col-104",
                                                 "dpy-7",
                                                 "rde-12",
                                                 "egl-15",
                                                 "svh-1"])
