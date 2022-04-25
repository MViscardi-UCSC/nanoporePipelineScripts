#!/usr/bin/python3
"""
roach_compareReplicates_run.py
Marcus Viscardi     Jan 28, 2021

Trying to compare between the adult replicates from the Roach et al. paper
"""
import matplotlib.pyplot as plt
import seaborn as sea
import pandas as pd

pd.set_option("display.max_columns", None)

ABS_PATH = f"/data16/marcus/working/210106_nanopolish_wRoachData"
OUTDIR1 = f"{ABS_PATH}/output_dir_adultrep1"
OUTDIR2 = f"{ABS_PATH}/output_dir_adultrep2"


def load_comp_on_genes(filepath) -> pd.DataFrame:
    def fix_lists(item, is_float=False, is_int=False):
        item_list = item.strip("[]").split(", ")
        if is_float:
            items = [float(z) for z in item_list]
        elif is_int:
            items = [int(z) for z in item_list]
        else:
            item_list = [str(s.strip("'")) for s in item_list]
            items = item_list
        return items

    temp_df = pd.read_csv(filepath, sep="\t", converters={"polya_lengths": lambda x: fix_lists(x, is_float=True),
                                                          "read_ids": lambda x: fix_lists(x)})
    return temp_df


def make_joint_plot():
    df1 = load_comp_on_genes(f"{OUTDIR1}/210124_compressed_on_genes.tsv")
    df2 = load_comp_on_genes(f"{OUTDIR2}/210128_compressed_on_genes.tsv")
    for i, df in enumerate((df1, df2)):
        print(f"\n\nREPLICATE {i}:\n")
        print(df[["gene_id", "read_hits", "polya_mean"]].head(15))
    merge_df = df1.merge(df2, how="inner", on="gene_id", suffixes=["_rep1", "_rep2"])
    # merge_df["hit_diff"] = merge_df['read_hits_rep1'].sub(merge_df['read_hits_rep2'], axis=0).abs()
    merge_df["hits_rank_rep1"] = merge_df['read_hits_rep1'].rank(ascending=False)
    merge_df["hits_rank_rep2"] = merge_df['read_hits_rep2'].rank(ascending=False)
    merge_df["hit_rank_diff"] = merge_df['hits_rank_rep1'].sub(merge_df['hits_rank_rep2'], axis=0).abs()

    print(merge_df[["gene_id", "hit_rank_diff", "polya_mean_rep1", "polya_mean_rep2"]].head(15))

    merge_df[["gene_id", "read_hits_rep1", "hits_rank_rep1", "read_hits_rep2", "hits_rank_rep2",
              "hit_rank_diff", "polya_mean_rep1", "polya_mean_rep2"]] \
        .to_csv(f"{ABS_PATH}/comparing_adultreplicates/merged_compressed_on_genes.tsv", sep="\t", index=False)

    g = sea.jointplot(data=merge_df,
                      x="polya_mean_rep1", y="polya_mean_rep2",
                      #hue="hit_rank_diff", alpha=0.6,
                      )
    # g.plot_marginals(sea.histplot, kde=True, zorder=1)
    # g.plot_joint(plt.scatter, marker="+", zorder=1)
    g.set_axis_labels(xlabel="Mean PolyA Tail Length Rep1", ylabel="Mean PolyA Tail Length Rep2")
    g.fig.suptitle("PolyA Mean Lengths per Gene between Replicates")
    g.fig.tight_layout()
    g.fig.subplots_adjust(top=0.95)
    plt.show()


def make_reg_plot():
    def my_palplot(pal, size=1, ax=None):
        """Plot the values in a color palette as a horizontal array.
        Parameters
        ----------
        pal : sequence of matplotlib colors
            colors, i.e. as returned by seaborn.color_palette()
        size :
            scaling factor for size of plot
        ax :
            an existing axes to use
        """
        import numpy as np
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    
        n = len(pal)
        if ax is None:
            f, ax = plt.subplots(1, 1, figsize=(n * size, size))
        ax.imshow(np.arange(n).reshape(1, n),
                  cmap=mpl.colors.ListedColormap(list(pal)),
                  interpolation="nearest", aspect="auto")
        ax.set_xticks(np.arange(n) - .5)
        ax.set_yticks([-.5, .5])
        # Ensure nice border between colors
        ax.set_xticklabels(["" for _ in range(n)])
        # The proper way to set no ticks
        ax.yaxis.set_major_locator(ticker.NullLocator())
    df1 = load_comp_on_genes(f"{OUTDIR1}/210124_compressed_on_genes.tsv")
    df2 = load_comp_on_genes(f"{OUTDIR2}/210128_compressed_on_genes.tsv")
    for i, df in enumerate((df1, df2)):
        print(f"\n\nREPLICATE {i}:\n")
        print(df[["gene_id", "read_hits", "polya_mean"]].head(15))
    merge_df = df1.merge(df2, how="inner", on="gene_id", suffixes=["_rep1", "_rep2"])
    # merge_df["hit_diff"] = merge_df['read_hits_rep1'].sub(merge_df['read_hits_rep2'], axis=0).abs()
    merge_df["hits_rank_rep1"] = merge_df['read_hits_rep1'].rank(ascending=False)
    merge_df["hits_rank_rep2"] = merge_df['read_hits_rep2'].rank(ascending=False)
    merge_df["hit_rank_diff"] = merge_df['hits_rank_rep1'].sub(merge_df['hits_rank_rep2'], axis=0).abs()
    merge_df["total_hits"] = merge_df[['read_hits_rep1', 'read_hits_rep2']].sum(axis=1)

    print(merge_df[["gene_id", "hit_rank_diff", "polya_mean_rep1", "polya_mean_rep2"]].head(5))

    merge_df[["gene_id", "read_hits_rep1", "hits_rank_rep1", "read_hits_rep2", "hits_rank_rep2",
              "hit_rank_diff", "polya_mean_rep1", "polya_mean_rep2"]] \
        .to_csv(f"{ABS_PATH}/comparing_adultreplicates/merged_compressed_on_genes.tsv", sep="\t", index=False)

    sea.set_palette('flare')
    #fig = plt.figure(figsize=(7, 7))
    # axs = plt.subplot2grid((4, 1), (0, 0), rowspan=3), plt.subplot2grid((4, 1), (3, 0))
    fig, ax = plt.subplots()
    ax.set_xlim((38, 150))
    ax.set_ylim((38, 150))
    for cutoff in range(0, 1100, 100):
        # TODO: Maybe try to manually set the color each cycle, (https://stackoverflow.com/a/28430505/13316742),
        #       From there try to add a legend point with the color and the "cutoff" value!!
        g = sea.scatterplot(data=merge_df[merge_df.total_hits >= cutoff],
                            x="polya_mean_rep1", y="polya_mean_rep2",
                            alpha=0.5, ax=ax
                            )
    g = sea.scatterplot(data=merge_df, x="polya_mean_rep1", y="polya_mean_rep2", alpha=0.0, hue="total_hits")
    g = sea.regplot(data=merge_df, x="polya_mean_rep1", y="polya_mean_rep2",
                    scatter=False, color="k", ax=ax, truncate=False,
                    )
    g = sea.regplot(data=merge_df[merge_df.total_hits >= 1000], x="polya_mean_rep1", y="polya_mean_rep2",
                    scatter=False, color="b", ax=ax, truncate=False,
                    )
    g.set(xlabel="Mean PolyA Tail Length Rep1", ylabel="Mean PolyA Tail Length Rep2",
          xscale="log", yscale="log", title="PolyA Mean Lengths per Gene between Replicates")
    plt.tight_layout()
    plt.show()


def josh_plot():
    df1 = load_comp_on_genes(f"{OUTDIR1}/210124_compressed_on_genes.tsv")
    df2 = load_comp_on_genes(f"{OUTDIR2}/210128_compressed_on_genes.tsv")
    for i, df in enumerate((df1, df2)):
        print(f"\n\nREPLICATE {i}:\n")
        print(df[["gene_id", "read_hits", "polya_mean"]].head(5))
    merge_df = df1.merge(df2, how="inner", on="gene_id", suffixes=["_rep1", "_rep2"])
    # merge_df["hit_diff"] = merge_df['read_hits_rep1'].sub(merge_df['read_hits_rep2'], axis=0).abs()
    merge_df["hits_rank_rep1"] = merge_df['read_hits_rep1'].rank(ascending=False)
    merge_df["hits_rank_rep2"] = merge_df['read_hits_rep2'].rank(ascending=False)
    merge_df["hit_rank_diff"] = merge_df['hits_rank_rep1'].sub(merge_df['hits_rank_rep2'], axis=0).abs()
    merge_df["mean_hits"] = merge_df[['read_hits_rep1', 'read_hits_rep2']].mean(axis=1)
    merge_df["mean_hit_std"] = merge_df[['read_hits_rep1', 'read_hits_rep2']].std(axis=1)
    merge_df["hit_diff"] = merge_df['read_hits_rep1'].sub(merge_df['read_hits_rep2'], axis=0).abs()
    merge_df["hit_diff/mean_hits"] = merge_df['hit_diff'].div(merge_df['mean_hits'], axis=0).abs()
    g = sea.jointplot(data=merge_df, x="mean_hits", y="hit_diff/mean_hits")
    g.ax_joint.set(xscale="log")
    plt.show()
    print(merge_df[["gene_id", "read_hits_rep1",
                    "read_hits_rep2", "mean_hits",
                    "hit_diff/mean_hits"]].sort_values("hit_diff/mean_hits", 
                                                       ascending=False).head(5))


def make_split_violins():
    df1, df2 = load_comp_on_genes(f"{OUTDIR1}/210124_compressed_on_genes.tsv"),\
               load_comp_on_genes("/data16/marcus/working/210709_NanoporeRun_totalRNA_0639_L3/output_dir/merge_files/210709_01:38:18PM_compressedOnGenes.tsv")
               # load_comp_on_genes(f"{OUTDIR2}/210128_compressed_on_genes.tsv")
    print(">> Dataframes loaded into memory.")
    new_dfs = []
    for i, df in enumerate((df1, df2)):
        i += 1
        df["hits_rank"] = df['read_hits'].rank(ascending=False)
        # TODO: drop the below line as it shrinks the dataset a bit to make it easier to manage:
        #df = df[df["hits_rank"] <= 15]
        df_tails = df["polya_lengths"].apply(pd.Series).merge(df["gene_id"], right_index=True, left_index=True)\
            .melt(id_vars=["gene_id"], value_name="tail_length").dropna().drop("variable", axis=1)
        # >>>
        # The below is meant to limit the long tail above the violin plots!
        # print(df_tails.sort_values("tail_length", ascending=False).head())
        df_tails["tail_length"] = df_tails["tail_length"].apply(lambda x: x if x <= 200 else 200)
        # print(df_tails.sort_values("tail_length", ascending=False).head())
        # <<<
        print(f">> Dataframe {i}'s tail lengths flattened")
        df_reads = df["read_ids"].apply(pd.Series).merge(df["gene_id"], right_index=True, left_index=True) \
            .melt(id_vars=["gene_id"], value_name="read_id").dropna().drop("variable", axis=1)
        print(f">> Dataframe {i}'s read ids flattened")
        df_merge = df_reads.merge(df_tails)
        df_merge = df_merge.merge(df[["gene_id", "polya_mean", "read_hits", "hits_rank"]])
        new_dfs.append(df_merge)
        print(f">> Dataframe {i} finished")
    df1_long, df2_long = new_dfs
    df1_long["replicate"] = "Replicate 1"
    df2_long["replicate"] = "Replicate 2"
    super_df = pd.concat([df1_long, df2_long], ignore_index=True)
    super_df = super_df.sort_values("hits_rank")
    print(super_df.info)

    # PLOT IT:
    # I have to manually pick genes because of the way I set up the hits_rank column!
    top_genes_dict = {"WBGene00001168": "eef-1A",
                      "WBGene00006930": "vit-6",
                      "WBGene00004409": "rla-1",
                      "WBGene00004494": "rps-25",
                      "WBGene00006728": "ubq-2",
                      "WBGene00006789": "unc-54"}
    top = super_df[super_df["gene_id"].isin(list(top_genes_dict.keys()))]
    top = top.replace({"gene_id": top_genes_dict})
    sea.set_theme(style="whitegrid")
    ax = sea.violinplot(x="gene_id", y="tail_length", hue="replicate",
                        data=top, split=True, bw=0.2, inner="quartiles")
    # ax = sea.stripplot(x="gene_id", y="tail_length", hue="replicate", data=top_2, dodge=True)
    ax.set(ylim=(0, None))
    plt.show()


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
        df = df[df["gene_id"].isin(gene_list)]
        df_tails = df["polya_lengths"].apply(pd.Series).merge(df["gene_id"], right_index=True, left_index=True) \
            .melt(id_vars=["gene_id"], value_name="tail_length").dropna().drop("variable", axis=1)
        print(f">> Dataframe {i}'s tail lengths flattened")
        df_reads = df["read_ids"].apply(pd.Series).merge(df["gene_id"], right_index=True, left_index=True) \
            .melt(id_vars=["gene_id"], value_name="read_id").dropna().drop("variable", axis=1)
        print(f">> Dataframe {i}'s read ids flattened")
        df_merge = df_reads.merge(df_tails)
        df_merge = df_merge.merge(df[["gene_id", "polya_mean", "read_hits"]])
        new_dfs.append(df_merge)
        print(f">> Dataframe {i} finished")
    df1_long, df2_long = new_dfs
    df1_long["replicate"] = "Replicate 1"
    df2_long["replicate"] = "Replicate 2"
    super_df = pd.concat([df1_long, df2_long], ignore_index=True)
    print(super_df.info())

    # # PLOT IT:
    # top_2 = concat_df[concat_df["hits_rank"] <= 2]
    sea.set_theme(style="whitegrid")
    # ax = sea.swarmplot(x="tail_length", y="gene_id", hue="replicate",
    #                    data=concat_df, dodge=True, alpha=.01, zorder=1)
    ax = sea.violinplot(x="gene_id", y="tail_length", hue="replicate",
                        data=super_df, split=True, bw=0.2, inner="quartiles")
    plt.show()
    print(" . . . Done!")


if __name__ == '__main__':
    # make_reg_plot()
    #josh_plot()
    
    make_split_violins()
    
    
    # df2 = load_comp_on_genes(f"{OUTDIR2}/210128_compressed_on_genes.tsv")
    # df_tails = df2["polya_lengths"].apply(pd.Series).merge(df2["gene_id"], right_index=True, left_index=True) \
    #     .melt(id_vars=["gene_id"], value_name="tail_length").dropna().drop("variable", axis=1)
    # print(df_tails.shape[0])
    # for max_len in range(0, 25, 5):
    #     print(df_tails[df_tails.tail_length <= max_len].shape[0])
