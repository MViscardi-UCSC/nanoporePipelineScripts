"""
plotyToCompare5TERAadapters.py
Marcus Viscardi,    December 12, 2021

Goal is to make a scatter between the two 5TERA libs
"""
import seaborn as sea
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

from nanoporePipelineCommon import pick_libs_return_paths_dict, gene_names_to_gene_ids

import pandas as pd

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)


def quick_test():
    merge_df = load_libraries()
    # fig = px.scatter(merge_df, x="t5_fraction_6", y="t5_fraction_wt", color='total_read_counts',
    #                  hover_name="gene_id", hover_data=["read_hits_6",
    #                                                    "read_hits_wt",
    #                                                    "t5_fraction_diff",
    #                                                    ],
    #                  trendline="ols")
    fig = px.scatter(merge_df, x="total_read_counts", y="t5_fraction_abs_diff",
                     color="t5_fraction_mean",
                     hover_name="gene_name", hover_data=["gene_id",
                                                         "read_hits_6",
                                                         "read_hits_wt",
                                                         "t5_fraction_diff",
                                                         "t5_fraction_wt",
                                                         "t5_fraction_6"],
                     log_x=True)
    fig.show()
    print("Done!")


def load_libraries(per_gene_cutoff=25):
    smg_6 = "xrn-1-5tera-smg-6"
    wt = "xrn-1-5tera"
    decode_dict = {smg_6: "smg-6",
                   wt: "wt"}
    df_dict = {}
    for lib, path in pick_libs_return_paths_dict([smg_6, wt], file_midfix='compressedOnGenes').items():
        lib = decode_dict[lib]
        df = pd.read_parquet(path)
        if "gene_id" not in df.columns:
            df.reset_index()
        if 'gene_name' not in df.columns:
            names_df = gene_names_to_gene_ids()
            df = df.merge(names_df, on='gene_id')
        df = df[["gene_id", "gene_name", "read_hits", "t5_fraction"]]
        df_dict[lib] = df
    print(df_dict)
    merge_df = pd.merge(df_dict["smg-6"], df_dict["wt"], on=["gene_id", "gene_name"], suffixes=("_6", "_wt"))
    merge_df['total_read_counts'] = merge_df['read_hits_6'] + merge_df['read_hits_wt']
    merge_df['t5_fraction_diff'] = (merge_df['t5_fraction_wt'] - merge_df['t5_fraction_6'])
    merge_df['t5_fraction_mean'] = (merge_df['t5_fraction_wt'] + merge_df['t5_fraction_6']) / 2
    merge_df['t5_fraction_abs_diff'] = merge_df['t5_fraction_diff'].abs()
    for suffix in ("_6", "_wt"):
        merge_df = merge_df[merge_df[f'read_hits{suffix}'] >= per_gene_cutoff]
    return merge_df


def with_dash_for_click_to_copy():
    merge_df = load_libraries()


if __name__ == '__main__':
    quick_test()
