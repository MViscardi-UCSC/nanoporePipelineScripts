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
from geneHeatmaps2 import load_tsv_and_assign_w_josh_method

def load_merged_on_reads(path_to_merged, lib_name: str = None):
    merged_on_reads_df = pd.read_csv(path_to_merged, sep="\t")
    merged_on_reads_df["read_len"] = merged_on_reads_df["sequence"].str.len()
    if lib_name:
        print(f"Finish loading library dataframe for: {lib_name}")
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


def load_for_cds_based_plotting(lib_list):
    
    # Loop through each library name in the list and for each:
    #   1. Load the TSV
    #   2. Merge this w/ Josh's assign reads based on chr_pos
    #   3. Create a column to retain the library identity
    #   3. Concatenate these dataframe into one large dataframe
    #       NOTE: This structure can be seperated again based on
    #       the "lib" column added in the previous step
    super_df = pd.DataFrame
    for library_name in lib_list:
        lib_df = load_tsv_and_assign_w_josh_method(library_name)
        lib_df["lib"] = library_name
        super_df = pd.concat([super_df, lib_df], ignore_index=True)
    return super_df


if __name__ == '__main__':
    pathdict = {
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
    
    run_with = ["polyA", "polyA2", "totalRNA2"]
    
    lonest_df = load_for_cds_based_plotting(run_with)
    # pathdict = {name: find_newest_matching_file(pathdict[name]) for name in run_with}
    # 
    # lib_df_dict = {name: load_merged_on_reads(path, lib_name=name) for name, path in pathdict.items()}
    # 
    # long_df = concat_reads_df(lib_df_dict)
    # # plotly_lib_ecdfs(long_df)
    # seaborn_lib_ecdfs(long_df)
    print("Done..")
