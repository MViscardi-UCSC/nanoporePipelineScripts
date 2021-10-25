"""
plotlyReadLengths.py
Marcus Viscardi,    October 19, 2021

Plotting read lengths, specifically genes from MtDNA
"""
import pandas as pd

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)

import plotly.express as px
from nanoporePipelineCommon import find_newest_matching_file
from geneHeatmaps2 import load_tsv_and_assign_w_josh_method, load_read_assignment_parquet


def load_merged_on_reads(path_to_merged, lib_name: str = None):
    if lib_name:
        print(f"Starting to load library dataframe for: {lib_name} . . .", end="")
    merged_on_reads_df = pd.read_csv(path_to_merged, sep="\t")
    merged_on_reads_df["read_len"] = merged_on_reads_df["sequence"].str.len()
    if lib_name:
        print(f"\bFinish loading library dataframe for: {lib_name}!")
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
    :param concatenated_df: 
    :return: 
    """
    fig = px.ecdf(concatenated_df, x="read_len", color="lib", marginal="histogram")
    fig.update_xaxes(range=[0, 3000])
    fig.show()


def seaborn_lib_ecdfs(concatenated_df):
    import seaborn as sea
    import matplotlib.pyplot as plt
    sea.ecdfplot(data=concatenated_df, x="read_len", hue="lib")
    plt.xlim(0, 3500)
    # plt.legend(loc='lower right')
    plt.show()


def load_for_cds_based_plotting(path_dict, drop_unassigned=True):
    read_assignment_df = load_read_assignment_parquet(f"/data16/marcus/genomes/elegansRelease100/"
                                                      f"Caenorhabditis_elegans.WBcel235.100.allChrs.parquet")
    # Loop through each library name in the list and for each:
    #   1. Load the TSV
    #   2. Merge this w/ Josh's assign reads based on chr_pos
    #   3. Create a column to retain the library identity
    #   3. Concatenate these dataframe into one large dataframe
    #       NOTE: This structure can be seperated again based on
    #       the "lib" column added in the previous step
    df_dict = {}
    for library_name, tsv_path in path_dict.items():
        lib_df = load_merged_on_reads(tsv_path, lib_name=library_name)
        # lib_df["lib"] = library_name
        # super_df = pd.concat([super_df, lib_df], ignore_index=True)
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
    print(f"\bFinished merge!")
    if drop_unassigned:
        super_df = super_df[~super_df["gene_id_fromAssign"].isna()].reset_index(drop=True)
        super_df = super_df[super_df["strand_fromReads"] == super_df["strand_fromAssign"]].reset_index(drop=True)
    return super_df


def main_plot_ecdfs(path_dict: dict, plotly_or_seaborn: str):
    lib_df_dict = {name: load_merged_on_reads(path, lib_name=name) for name, path in path_dict.items()}
    long_df = concat_reads_df(lib_df_dict)
    if plotly_or_seaborn == "plotly":
        plotly_lib_ecdfs(long_df)
    elif plotly_or_seaborn == "seaborn":
        seaborn_lib_ecdfs(long_df)
    else:
        print(f"Please answer 'plotly_or_seaborn' param w/ 'plotly' or 'seaborn'... not '{plotly_or_seaborn}'!!")


def pick_libs_return_paths_dict(lib_list: list):
    path_dict = {
        "riboD": "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/output_dir/"
                 "merge_files/*_mergedOnReads.tsv",
        "totalRNA": "/data16/marcus/working/210709_NanoporeRun_totalRNA_0639_L3/"
                    "output_dir/merge_files/*_mergedOnReads.tsv",
        "totalRNA2": "/data16/marcus/working/"
                     "210720_nanoporeRun_totalRNA_0639_L3_replicate/output_dir/"
                     "merge_files/*_mergedOnReads.tsv",
        "polyA": "/data16/marcus/working/210528_NanoporeRun_0639_L3s/output_dir/"
                 "merge_files/*_mergedOnReads.tsv",
        "polyA2": "/data16/marcus/working/210719_nanoporeRun_polyA_0639_L3_replicate/"
                  "output_dir/merge_files/*_mergedOnReads.tsv",
        "xrn-1": "/data16/marcus/working/210905_nanoporeRun_totalRNA_5108_xrn-1-KD/"
                 "output_dir/merge_files/*_mergedOnReads.tsv"
    }
    return_dict = {}
    for lib_key, tsv_path in path_dict.items():
        if lib_key in lib_list:
            return_dict[lib_key] = find_newest_matching_file(tsv_path)
    return return_dict


if __name__ == '__main__':
    run_with = ["polyA", "polyA2", "totalRNA2"]
    lib_path_dict = pick_libs_return_paths_dict(run_with)

    lonest_df = load_for_cds_based_plotting(lib_path_dict)
    # main_plot_ecdfs()
    print("Done..")
