"""
comparingMyRunToRoach.py
Marcus Viscardi     June 2, 2021

Based on the read counts it seems like I have a MUCH greater (~5X) read depth! To see if this is just an
    effect of my reads being shorter but having similar nucleotide depth leading to a skewed view, I want
    to quantify how many called and mapped bases each data set had.

Going to jump straing to the sorted sam files as these will be the most simple way to do this.
"""

import pandas as pd
from dashForClickAndView import load_and_merge_lib_parquets

pd.set_option("display.max_columns", None)


def sam_to_df(path_to_sam) -> pd.DataFrame:
    print(f"Starting import from file @ {path_to_sam}")
    sam_df = pd.read_csv(path_to_sam,
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
    return sam_df


def count_nts_mapped(df_from_sam: pd.DataFrame):
    print(f"\tCalculating read lengths for {df_from_sam.shape[0]} lines of .sam file")
    df_from_sam["read_length"] = df_from_sam["sequence"].str.len()
    print("\tCalculating sum, mean, median, max, and min")
    nt_sum = df_from_sam["read_length"].sum()
    mean = df_from_sam["read_length"].mean()
    median = df_from_sam["read_length"].median()
    min = df_from_sam["read_length"].min()
    max = df_from_sam["read_length"].max()
    return nt_sum, mean, median, min, max, df_from_sam.shape[0]


def old_main():
    path_to_roach = "/data16/marcus/prefix/210106_nanopolish_wRoachData/output_dir_youngadultrep1tech1/cat_allReads.sam"
    path_to_my = "/data16/marcus/prefix/210528_NanoporeRun_0639_L3s/output_dir/cat_files/cat.sorted.sam"
    dictionary = {"Roach data (young adult worms - replicate 1)": path_to_roach,
                  "My data (L3 worms)": path_to_my}
    results = []
    for name, path in dictionary.items():
        print(f"Working on {name}:")
        df = sam_to_df(path)
        results.append(count_nts_mapped(df))
    returns = results
    print(f"         Roach:\t Viscardi:\n"
          f"Reads:\t{returns[0][5]:10.0f}\t{returns[1][5]:10.0f}\n"
          f"Sum:  \t{returns[0][0]:10.0f}\t{returns[1][0]:10.0f}\n"
          f"Mean: \t{returns[0][1]:10.0f}\t{returns[1][1]:10.0f}\n"
          f"Median:\t{returns[0][2]:10.0f}\t{returns[1][2]:10.0f}\n"
          f"Min:  \t{returns[0][3]:10.0f}\t{returns[1][3]:10.0f}\n"
          f"Max:  \t{returns[0][4]:10.0f}\t{returns[1][4]:10.0f}")


def seaborn_tail_length_scatter(lib_list, min_hits, plot_column="mean_polya_length"):
    import seaborn as sea
    import matplotlib.pyplot as plt
    
    if len(lib_list) != 2:
        raise NotImplementedError(f"Please only provide 2 libraries, you passed: {lib_list}")
    
    reads_df, compressed_df = load_and_merge_lib_parquets(lib_list)
    x_lib, y_lib = lib_list

    min_hit_df = compressed_df[compressed_df['gene_hits'] >= min_hits]
    x_axis_df = min_hit_df[min_hit_df.lib == x_lib][["gene_id",
                                                             "gene_name",
                                                             "gene_hits",
                                                             "mean_polya_length"]]
    y_axis_df = min_hit_df[min_hit_df.lib == y_lib][["gene_id",
                                                             "gene_name",
                                                             "gene_hits",
                                                             "mean_polya_length"]]

    plot_df = pd.merge(x_axis_df, y_axis_df, on=["gene_id", "gene_name"],
                       suffixes=(f"_{x_lib}",
                                 f"_{y_lib}"))
    max_mean = min_hit_df[plot_column].max()
    max_mean += 10
    min_mean = min_hit_df[plot_column].min()
    min_mean -= 10
    
    sea.set_style("whitegrid")
    sea.set_context("notebook")
    plt.style.context("seaborn-whitegrid")
    plot_size = 5
    fig, ax = plt.subplots(figsize=(plot_size, plot_size))
    
    fig = sea.scatterplot(data=plot_df, ax=ax,
                          x=f"{plot_column}_{x_lib}",
                          y=f"{plot_column}_{y_lib}",
                          color=(0.2, 0.2, 0.2, 0.5))
    fig.set_xlim(round(min_mean/10)*10, round(max_mean/10)*10)
    fig.set_xticks(range(round(min_mean/10)*10, round(max_mean/10)*10, 15))
    fig.set_ylim(round(min_mean/10)*10, round(max_mean/10)*10)
    fig.set_yticks(range(round(min_mean/10)*10, round(max_mean/10)*10, 15))
    fig.set_title(f"Mean poly(A) Tail Lengths per Gene called by Nanopolish\n(hits cutoff of {min_hits}reads/gene)")
    fig.set_xlabel(f"Mean Tail Length Per Gene from: {x_lib} (nts)")
    fig.set_ylabel(f"Mean Tail Length Per Gene from: {y_lib} (nts)")
    sea.despine()
    plt.tight_layout()
    plt.show()
    print("done?")


if __name__ == '__main__':
    seaborn_tail_length_scatter(["totalRNA2", "polyA2"], 40)
