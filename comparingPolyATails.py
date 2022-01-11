"""
compariongPolyATails.py
Marcus Viscardi     June 2, 2021

Based on the read counts it seems like I have a MUCH greater (~5X) read depth! To see if this is just an
    effect of my reads being shorter but having similar nucleotide depth leading to a skewed view, I want
    to quantify how many called and mapped bases each data set had.

Going to jump straing to the sorted sam files as these will be the most simple way to do this.
"""

import pandas as pd
from dashForClickAndViewPolyATails import load_and_merge_lib_parquets

pd.set_option("display.max_columns", None)


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
