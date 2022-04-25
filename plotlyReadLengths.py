"""
plotlyReadLengths.py
Marcus Viscardi,    October 19, 2021

Plotting read lengths, specifically genes from MtDNA
"""
import pandas as pd
import numpy as np

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)

import plotly.express as px
import seaborn as sea
import matplotlib.pyplot as plt

from nanoporePipelineCommon import find_newest_matching_file, get_dt, load_read_assignments, pick_libs_return_paths_dict
from geneHeatmaps2 import load_tsv_and_assign_w_josh_method


def load_merged_on_reads(path_to_merged, lib_name: str = None, head=None):
    if lib_name:
        print(f"Starting to load library dataframe for: {lib_name} . . .", end="")
    if path_to_merged.endswith(".parquet"):
        merged_on_reads_df = pd.read_parquet(path_to_merged)
        if isinstance(head, int):
            merged_on_reads_df = merged_on_reads_df.head(head)
    elif path_to_merged.endswith(".tsv"):
        merged_on_reads_df = pd.read_csv(path_to_merged, sep="\t", nrows=head)
    else:
        raise NotImplementedError(f"Please provide file that ends w/ .tsv or .parquet, not:\n{path_to_merged}")
    merged_on_reads_df["read_len"] = merged_on_reads_df["sequence"].str.len()
    if lib_name:
        print(f"\rFinish loading library dataframe for: {lib_name}!")
    return merged_on_reads_df


def concat_reads_df(df_dict: dict):
    super_df = pd.DataFrame()
    for name, df in df_dict.items():
        df["lib"] = name
        super_df = pd.concat([super_df, df], ignore_index=True)
        print(f"Dataframe concatenated for library: {name}")
    return super_df


def plotly_lib_ecdfs(concatenated_df):
    """
    This creates a really beatiful plot that you can zoom into, but it is SOOOOO laggy
    I am going to try seaborn to make this a little easier to export
    """
    fig = px.ecdf(concatenated_df, x="read_len", color="lib", marginal="histogram")
    fig.update_xaxes(range=[0, 3000])
    fig.show()


def seaborn_lib_ecdfs(concatenated_df, x_limits=(1, 3500), save=True, save_as=False):
    sea.set_style("whitegrid")
    sea.set_context("notebook")
    plt.style.context("seaborn-whitegrid")
    fig, ax = plt.subplots(figsize=(6, 6))
    sea.ecdfplot(data=concatenated_df, x="read_length", hue="lib", ax=ax)
    ax.set(xscale="log")
    x_min, x_max = x_limits
    plt.xlim(x_min, x_max)
    # plt.legend(loc='lower right')
    if save and not save_as:
        plt.savefig(f"./testOutputs/{get_dt(for_file=True)}_readLen_ecdfs.svg")
    if isinstance(save_as, str):
        plt.savefig(save_as)
    plt.show()


def load_for_cds_based_plotting(path_dict, drop_unassigned=True, subset=None):
    read_assignment_df = load_read_assignments(f"/data16/marcus/genomes/elegansRelease100/"
                                               f"Caenorhabditis_elegans.WBcel235.100.allChrs.parquet")
    # Loop through each library name in the list and for each:
    #   1. Load the TSV
    #   2. Merge this w/ Josh's assign reads based on chr_pos
    #   3. Create a column to retain the library identity
    #   3. Concatenate these dataframe into one large dataframe
    #       NOTE: This structure can be seperated again based on
    #       the "lib" column added in the previous step
    df_dict = {}
    for library_name, file_path in path_dict.items():
        lib_df = load_merged_on_reads(file_path, lib_name=library_name, head=subset)
        df_dict[library_name] = lib_df

    # This is a cute way to quickly merge all of these dfs into one, while retaining lib info.
    #   B/c I am still a little scared of MultiIndexed dataframes, I used the reset_index steps
    #   to push the mutliindex back into columns. Maybe someday I'll use the multiindex!
    multi_df = pd.concat(df_dict.values(), keys=df_dict.keys())
    multi_df.index.set_names(("lib", "old_index"), inplace=True)
    super_df = multi_df.reset_index(level="lib").reset_index(drop=True)
    print(f"Starting merge . . .", end="")
    super_df = super_df.merge(read_assignment_df, on=["chr_id", "chr_pos"],
                              how="left", suffixes=["_fromReads",
                                                    "_fromAssign"])
    print(f"\rFinished merge!")
    if drop_unassigned:
        print(f"Read counts post gene assignment:  {super_df.shape[0]}")
        super_df = super_df[~super_df["gene_id_fromAssign"].isna()].reset_index(drop=True)
        print(f"Read counts post unassigned drop:  {super_df.shape[0]}")
        super_df = super_df[super_df["strand_fromReads"] == super_df["strand_fromAssign"]].reset_index(drop=True)
        print(f"Read counts post assignment check: {super_df.shape[0]}")
    return super_df


def handle_long_df_and_plot(super_df: pd.DataFrame, distance_from_start_cutoff=0, plotly_or_seaborn=None,
                            genes_or_txns="txns", filter_hits_less_than: int = 1):
    from tqdm import tqdm
    # First lets add a column to hold the past_start information
    super_df["to_start"] = super_df["to_start"].astype('int64')
    super_df["to_stop"] = super_df["to_stop"].astype('int64')
    tqdm.pandas(desc="Calculating dist past start codons")
    super_df["past_start"] = super_df.progress_apply(lambda row: _reads_past_start(row["to_start"],
                                                                                   cut_off=distance_from_start_cutoff),
                                                     axis=1)
    tqdm.pandas(desc="Calculating cds lengths")
    super_df["cds_len"] = super_df.progress_apply(lambda row: abs(row["to_start"] - row["to_stop"] + 4),
                                                  axis=1)
    print("Finished calculating read lengths")
    compress_list = ["lib", "gene_id_fromAssign"]
    if genes_or_txns == "txns":
        compress_list.append("transcript_id")
    group_by_txs = super_df.groupby(by=compress_list)
    grouped_df = pd.DataFrame(group_by_txs["past_start"].apply(np.mean))
    grouped_df["transcript_hits"] = group_by_txs["past_start"].apply(len)
    grouped_df = grouped_df[grouped_df["transcript_hits"] > filter_hits_less_than]
    if genes_or_txns == 'txns':
        tqdm.pandas(desc="Checking for only one cds per transcript")
    else:
        tqdm.pandas(desc="Getting average cds length for genes w/ > transcript")
    grouped_df["cds_len"] = group_by_txs["cds_len"].apply(lambda cds_lens:
                                                          _test_if_cds_len_consistent(cds_lens,
                                                                                      genes_or_txns=genes_or_txns,
                                                                                      zero_multi_txn_genes=True))
    grouped_df = grouped_df.query("cds_len != 0")

    # Binning data with the pandas cut function:
    cds_bins = [0, 250, 500, 750, 1000, 1250, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 10000]
    cds_bin_names = [f"{bin_start} to {cds_bins[i + 1]}" for i, bin_start in enumerate(cds_bins[:-1])]
    grouped_df["binned_cds_len"] = pd.cut(grouped_df["cds_len"], bins=cds_bins, labels=cds_bin_names)

    print("stopping point 1")

    # Plot it!
    if plotly_or_seaborn == "plotly":
        hover_col_list = ["gene_id_fromAssign", "transcript_hits"]
        if genes_or_txns == 'txns':
            hover_col_list.append('transcript_id')
        fig = px.scatter(grouped_df.reset_index().sort_values(['binned_cds_len', 'lib']),
                         x="cds_len", y="past_start",
                         color="lib", # points='outliers',
                         color_discrete_sequence=['#202020', '#A0A0A0', '#606060', '#C0C0C0'],
                         hover_data=hover_col_list,
                         category_orders=dict(zip(cds_bin_names, cds_bin_names)))
        fig.update_layout(template='plotly_white')
        fig.write_html("/home/marcus/Downloads/220302_deleteme.html")
        fig.show()
    elif plotly_or_seaborn == "seaborn":
        # first plot the strip plot, showing all the raw data
        fig = sea.stripplot(x="binned_cds_len",
                            y="past_start",
                            data=grouped_df,
                            size=4, color=".7",
                            order=cds_bin_names)

        # rotate the x labels so their easier to read!
        fig.set_xticklabels(fig.get_xticklabels(),
                            rotation=40, rotation_mode="anchor", ha='right')

        # plot the mean and median lines (mean in black, median in red!)
        sea.boxplot(showmeans=True,
                    meanline=True,
                    meanprops={'color': 'k',
                               'ls': '-',
                               'lw': 2},
                    medianprops={'visible': True,
                                 'color': 'r',
                                 'lw': 2},
                    whiskerprops={'visible': False},
                    zorder=10,
                    x="binned_cds_len",
                    y="past_start",
                    data=grouped_df,
                    showfliers=False,
                    showbox=False,
                    showcaps=False,
                    order=cds_bin_names,
                    ax=fig)
        # tight layout so everything is visible
        plt.tight_layout()
        plt.savefig(f"./testOutputs/{get_dt(for_file=True)}_readLen_binned.svg")
        plt.show()
    elif plotly_or_seaborn == "seaborn_cat":
        g = sea.catplot(data=grouped_df.reset_index(),
                        x="binned_cds_len",
                        y="past_start",
                        row="lib",
                        kind="strip",
                        order=cds_bin_names,
                        aspect=1.5, height=5,
                        color='gray',
                        jitter=0.25,
                        alpha=0.5,
                        )
        g.data = grouped_df.reset_index()
        g.map(sea.boxplot, x="binned_cds_len", y="past_start",
              showmeans=True,
              meanline=True,
              meanprops={'color': 'k',
                         'ls': '-',
                         'lw': 2},
              medianprops={'visible': True,
                           'color': 'r',
                           'lw': 2},
              whiskerprops={'visible': False},
              zorder=10,
              data=grouped_df,
              order=cds_bin_names,
              showfliers=False,
              showbox=False,
              showcaps=False)

        # Add a custom legend:
        from matplotlib.lines import Line2D
        custom_lines = [Line2D([0], [0], color='k', lw=2),
                        Line2D([0], [0], color='r', lw=2)]
        plt.legend(custom_lines, ["Mean", "Median"])

        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        save_fig_path = f"./testOutputs/{get_dt(for_file=True)}_readLen_binned_catplot.{genes_or_txns}.svg"
        plt.savefig(save_fig_path)
        plt.show()

    print("stopping point 2")


def _reads_past_start(to_start, cut_off: int = 0, is_past=100, is_not_past=0) -> int:
    if to_start <= cut_off:
        return_val = is_past
    elif to_start > cut_off:
        return_val = is_not_past
    else:
        return_val = is_not_past
    return return_val


def _test_if_cds_len_consistent(cds_lengths, genes_or_txns=None, zero_multi_txn_genes=False) -> int:
    cds_set = set(cds_lengths.to_list())
    cds_set_len = len(cds_set)
    if genes_or_txns == "txns" and cds_set_len > 1:
        raise NotImplementedError(f"Multiple CDS lengths?: {cds_set}")
    elif genes_or_txns == "genes" and cds_set_len > 1:
        if zero_multi_txn_genes:
            return 0
        else:
            return int(np.mean(list(cds_set)))
    else:  # Only one cds length:
        return int(next(iter(cds_set)))


def main_plot_ecdfs(path_dict: dict, plotly_or_seaborn: str):
    lib_df_dict = {name: load_merged_on_reads(path, lib_name=name) for name, path in path_dict.items()}
    long_df = concat_reads_df(lib_df_dict)
    if plotly_or_seaborn == "plotly":
        plotly_lib_ecdfs(long_df)
    elif plotly_or_seaborn == "seaborn":
        seaborn_lib_ecdfs(long_df)
    else:
        print(f"Please answer 'plotly_or_seaborn' param w/ 'plotly' or 'seaborn'... not '{plotly_or_seaborn}'!!")


def plot_ecdf_211109(libs):
    from dashForClickAndViewPolyATails import load_and_merge_lib_parquets
    reads_df, compressed_df = load_and_merge_lib_parquets(libs, drop_failed_polya=False,
                                                          drop_unassigned=False)
    seaborn_lib_ecdfs(reads_df)
    return reads_df, compressed_df


def plot_scatter_211109(libs, compressed_df, plot_column, min_hits=40, save=True):
    x_lib, y_lib = libs[:2]

    min_hit_df = compressed_df[compressed_df['gene_hits'] >= min_hits]
    x_axis_df = min_hit_df[min_hit_df.lib == x_lib][["gene_id",
                                                     "gene_name",
                                                     "gene_hits",
                                                     "mean_polya_length",
                                                     "mean_read_length"]]
    y_axis_df = min_hit_df[min_hit_df.lib == y_lib][["gene_id",
                                                     "gene_name",
                                                     "gene_hits",
                                                     "mean_polya_length",
                                                     "mean_read_length"]]

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
    fig.set_xlim(round(min_mean / 10) * 10, round(max_mean / 10) * 10)
    # fig.set_xticks(range(round(min_mean/10)*10, round(max_mean/10)*10, 500))
    fig.set_ylim(round(min_mean / 10) * 10, round(max_mean / 10) * 10)
    # fig.set_yticks(range(round(min_mean/10)*10, round(max_mean/10)*10, 500))
    fig.set_title(f"Mean Read Tail Lengths per Gene\n(hits cutoff of {min_hits}reads/gene)")
    fig.set_xlabel(f"Mean Read Length Per Gene from: {x_lib} (nts)")
    fig.set_ylabel(f"Mean Read Length Per Gene from: {y_lib} (nts)")
    fig.set(yscale='log', xscale='log')
    ticks = [100, 300, 1000, 3000, 10000]
    labels = [i for i in ticks]
    fig.set(yticks=ticks, yticklabels=labels,
            xticks=ticks, xticklabels=labels)
    sea.despine()
    plt.tight_layout()
    if save:
        plt.savefig(f"./testOutputs/{get_dt(for_file=True)}_readLen_scatter.svg")
    plt.show()


if __name__ == '__main__':
    run_with = ["polyA3", "totalRNA3", "polyA2", "totalRNA2"]
    force_long_df_build = False

    if not force_long_df_build:
        search_path = f"./testInputs/*_{'-'.join(run_with)}.longForReadLengths.parquet"
        try:
            longest_df_path = find_newest_matching_file(search_path)
            print(f"Found preprocessed at {longest_df_path}\nLoading now.", end='')
            longest_df = pd.read_parquet(longest_df_path)
            print(" Done.")
        except ValueError:
            print(f"Couldn't find pre-processed file at: {search_path}\nGoing to load from library files!")
            force_long_df_build = True
    if force_long_df_build:
        lib_path_dict = pick_libs_return_paths_dict(run_with)
        longest_df = load_for_cds_based_plotting(lib_path_dict, subset=None)
        save_path = f"./testInputs/{get_dt(for_file=True)}_{'-'.join(run_with)}.longForReadLengths.parquet"
        longest_df.to_parquet(save_path)
        print(f"Saved new long form file to: {save_path}")
    handle_long_df_and_plot(longest_df, distance_from_start_cutoff=50,
                            plotly_or_seaborn="plotly",
                            genes_or_txns="genes",
                            filter_hits_less_than=40)
    # main_plot_ecdfs(lib_path_dict, plotly_or_seaborn="seaborn")
    # main_plot_ecdfs(lib_path_dict, plotly_or_seaborn="seaborn")

    # reads_df, compressed_df = plot_ecdf_211109(run_with)
    # plot_scatter_211109(run_with, compressed_df, "mean_read_length")
    print("Done..")
