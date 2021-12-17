"""
compareRunsWithPlotly.py
Marcus Viscardi, 8/3/21


"""
import os
import pandas as pd
import scipy.stats
from nanoporePipelineCommon import find_newest_matching_file, pick_libs_return_paths_dict, get_dt

import warnings
from pandas.core.common import SettingWithCopyWarning

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
pd.set_option("display.max_columns", None)


def load_df(path_to_file: str) -> pd.DataFrame:
    if path_to_file.endswith(".tsv"):
        return pd.read_csv(path_to_file, sep="\t")
    elif path_to_file.endswith(".parquet"):
        return pd.read_parquet(path_to_file)


def load_and_merge_from_dict(path_dict: dict, min_hit_cutoff: int = None,
                             add_names: bool = True, add_gc_frac: bool = True) -> (pd.DataFrame, list):
    test_path = list(path_dict.values())[0]
    if test_path.endswith(".tsv"):
        df_dict = {key: pd.read_csv(path, sep="\t") for (key, path) in path_dict.items()}
    elif test_path.endswith(".parquet"):
        df_dict = {key: pd.read_parquet(path) for (key, path) in path_dict.items()}
    else:
        raise NotImplementedError(f"Need file paths that end with parquet or tsv, not:\n{test_path}\n{path_dict}")
    df_key_list = list(df_dict.keys())
    print(f"\nCreating dataframes for: {df_key_list}...")
    for key, df in df_dict.items():
        print(f"{key:>10} dataframe length: {df.shape[0]:>6}")
    midmerge_df = df_dict[df_key_list[0]].merge(df_dict[df_key_list[1]],
                                                how="inner", on="gene_id",
                                                suffixes=["_" + df_key_list[0], "_" + df_key_list[1]])
    if len(df_key_list) >= 3:
        df_dict[df_key_list[2]] = df_dict[df_key_list[2]].add_suffix("_" + df_key_list[2])

        merge_df = midmerge_df.merge(df_dict[df_key_list[2]],
                                     how="inner", left_on="gene_id", right_on=f"gene_id_{df_key_list[2]}")
    else:
        merge_df = midmerge_df
    print(f"{'merge_df':>10} dataframe length: {merge_df.shape[0]:>6}")
    if add_names:
        path_to_names = "/data16/marcus/genomes/elegansRelease100/" \
                        "Caenorhabditis_elegans.WBcel235.100.gtf.parquet"
        names_df = pd.read_parquet(path_to_names)
        names_df = names_df[["gene_id", "gene_name", "chr"]].drop_duplicates()
        names_df["is_MtDNA"] = names_df["chr"] == "MtDNA"
        names_df["is_MtDNA"] = names_df["is_MtDNA"].replace([True, False], ["MtDNA", "Not_MtDNA"])

        merge_df = merge_df.merge(names_df, on="gene_id", how="left")
        merge_df['ident'] = merge_df["gene_name"] + " (" + merge_df["gene_id"] + ")"
    if add_gc_frac:
        path_to_gc = "/data16/marcus/genomes/elegansRelease100/" \
                        "Caenorhabditis_elegans.WBcel235.cdna.all.fa.GCcontent.parquet"
        gc_df = pd.read_parquet(path_to_gc).drop(columns=["chr_id"])

        merge_df = merge_df.merge(gc_df, on=["gene_id", "gene_name"], how="left")
    else:
        merge_df['ident'] = merge_df["gene_id"]
    if min_hit_cutoff:
        cutoff_merge_df = merge_df
        print(f"Applying cutoff of {min_hit_cutoff} reads/gene...")
        for key in df_key_list:
            cutoff_merge_df = cutoff_merge_df[cutoff_merge_df[f"read_hits_{key}"] >= min_hit_cutoff]
        for key in df_key_list:
            cutoff_merge_df[f"hits_rank_{key}"] = cutoff_merge_df[f'read_hits_{key}'].rank(ascending=False)
        print(f"{'cutoff_df':>10} dataframe length: {cutoff_merge_df.shape[0]:>6}\n")
        return cutoff_merge_df, df_key_list
    else:
        for key in df_key_list:
            merge_df[f"hits_rank_{key}"] = merge_df[f'read_hits_{key}'].rank(ascending=False)
        return merge_df, df_key_list


def plotly_from_triple_merge(merged_df, key_list, cutoff=None,
                             x=0, y=1, compare_column_prefix="polya_mean",
                             color_by: str = None, ols_trend: bool = True,
                             save_dir: str = None):
    import plotly.express as px
    import plotly.graph_objects as go
    if isinstance(x, int) and isinstance(y, int):
        x_key = key_list[x]
        y_key = key_list[y]
    elif isinstance(x, str) and isinstance(y, str):
        x_key = x
        y_key = y

    # Set title variable, based on what's being plotted
    if compare_column_prefix == "polya_mean":
        title = f"PolyA Mean Lengths per Gene between RNA Selection Methods"
    elif compare_column_prefix == "read_hits":
        title = f"Number of Reads per Gene between RNA Selection Methods ({x_key} & {y_key})"
    else:
        title = f"Comparing: {compare_column_prefix} between {x_key} and {y_key}"
    
    # Set trendline variable
    if ols_trend:
        ols_trend = "ols"
    else:
        ols_trend = None
    merged_df["tail_avg_stdev"] = merged_df[f'polya_stdev_{x_key}']. \
        add(merged_df[f'polya_stdev_{y_key}'], axis=0).abs().div(2)
    merged_df["read_len_std_mean"] = merged_df[f'read_len_std_{x_key}']. \
        add(merged_df[f'read_len_std_{y_key}'], axis=0).abs().div(2)
    merged_df["tail_length_diff"] = merged_df[f'polya_mean_{x_key}']. \
        sub(merged_df[f'polya_mean_{y_key}'], axis=0).abs()
    merged_df["tail_length_mean_mean"] = merged_df[f'polya_mean_{x_key}']. \
        add(merged_df[f'polya_mean_{y_key}'], axis=0).abs().div(2)
    merged_df["read_len_mean_mean"] = merged_df[f'read_len_mean_{x_key}']. \
        add(merged_df[f'read_len_mean_{y_key}'], axis=0).abs().div(2)
    if isinstance(color_by, str):
        color_column = color_by
    else:
        color_column = None
    labels_dict = {"gene_id": "WBGene ID",
                   "ident": "Identity",
                   "tail_length_diff": "Mean Tail Diff (nts)",
                   "read_rank_diff": "Difference in Read Rank",
                   "chr": "Chromosome",
                   "gene_gc": "GC Fraction from Annot's"}
    for key in [x_key, y_key]:
        labels_dict[f"polya_mean_{key}"] = f"{key} - Mean Tail Length (nts)"
        labels_dict[f"read_hits_{key}"] = f"{key} - Reads/Gene"
        labels_dict[f"hits_rank_{key}"] = f"{key} - Rank of Reads/Gene"
        labels_dict[f"read_len_mean_{key}"] = f"{key} - Mean Read Length/Gene"
    # Filter for only MtDNA genes?
    # merged_df = merged_df[merged_df["chr_id"] == "MtDNA"]
    fig = px.scatter(merged_df, x=f"{compare_column_prefix}_{x_key}", y=f"{compare_column_prefix}_{y_key}",
                     color=color_column,
                     opacity=.8,
                     hover_name='ident',
                     hover_data=[f"polya_mean_{x_key}", f"polya_mean_{y_key}",
                                 f'read_hits_{x_key}', f'read_hits_{y_key}',
                                 f"read_len_mean_{x_key}", f"read_len_mean_{y_key}",
                                 f"chr", f"gene_gc",
                                 "a_fraction", "t_fraction",
                                 "c_fraction", "g_fraction"],
                     title=title,
                     labels=labels_dict,
                     trendline=ols_trend,
                     color_continuous_scale='Bluered',
                     # marginal_x="violin", marginal_y="violin",
                     )
    if compare_column_prefix == "polya_mean":
        fig.add_trace(go.Scatter(x=[0, 300], y=[0, 300], line=dict(color="#4a4a4a", width=1, dash='dash'),
                                 showlegend=False, mode="lines"))
        min = 30
        max = 120
        fig.update_xaxes(range=[min, max])
        fig.update_yaxes(range=[min, max])
    elif compare_column_prefix == "read_hits":
        fig.update_xaxes(type="log")
        fig.update_yaxes(type="log")
        fig.add_trace(go.Scatter(x=[0, 10000], y=[0, 10000], line=dict(color="#4a4a4a", width=2, dash='dash'),
                                 showlegend=False, mode="lines"))
    elif compare_column_prefix == "hits_rank":
        fig.add_trace(go.Scatter(x=[0, 400], y=[0, 400], line=dict(color="#4a4a4a", width=1, dash='dash'),
                                 showlegend=False, mode="lines"))
    elif compare_column_prefix == "read_len_mean":
        fig.update_xaxes(type="log")
        fig.update_yaxes(type="log")
        fig.add_trace(go.Scatter(x=[0, 2000], y=[0, 2000], line=dict(color="#4a4a4a", width=1, dash='dash'),
                                 showlegend=False, mode="lines"))

    spearman = scipy.stats.spearmanr(merged_df[f'{compare_column_prefix}_{x_key}'].tolist(),
                                     merged_df[f'{compare_column_prefix}_{y_key}'].tolist())

    fig.add_annotation(dict(x=0, y=1.06, showarrow=False,
                            text=f"Dropped any genes w/ less than {cutoff} reads ({merged_df.shape[0]} genes passed)",
                            textangle=0, xref="paper", yref="paper"))
    fig.add_annotation(dict(x=0, y=1.03, showarrow=False,
                            text=f"{spearman}",
                            textangle=0, xref="paper", yref="paper"))
    fig.update_layout(template='plotly_white')

    fig.show()
    save_name = f"/{get_dt(for_file=True)}_assessing-{compare_column_prefix}" \
                f"{'-'+color_column if color_column else ''}_between-{x_key}-{y_key}"
    fig.write_image("./testOutputs" + save_name + ".svg",
                    width=600, height=600, scale=2)
    if isinstance(save_dir, str):
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        save_name = save_dir + save_name
        fig.write_image(save_name + ".png",
                        width=600, height=600, scale=2)


def plotter_helper(path_dict: dict, prefix: str, cut_off: int, one_to_drop: str = "",
                   color_by: str = None, drop_mt_dna: bool = False, drop_gc_less: bool = True,
                   ols_trend: bool = True, save_dir: str = None,):
    key_list = list(path_dict.keys())
    if one_to_drop:
        subset_list = [i for i in key_list if i != one_to_drop]
        subset_dict = {key: path_dict[key] for key in subset_list}
    else:
        subset_dict = path_dict
    merge_df, keys = load_and_merge_from_dict(subset_dict, min_hit_cutoff=cutoff)
    if drop_mt_dna:
        merge_df = merge_df[merge_df["chr"] != "MtDNA"]
    if drop_gc_less:
        merge_df = merge_df[~merge_df["gene_gc"].isna()]
    plotly_from_triple_merge(merge_df, keys,
                             cutoff=cut_off,
                             x=0,
                             y=1,
                             compare_column_prefix=prefix,
                             color_by=color_by,
                             ols_trend=ols_trend,
                             save_dir=save_dir,
                             )


if __name__ == '__main__':

    run_with = ["xrn-1-5tera-smg-6", "xrn-1-5tera"]

    pathdict = pick_libs_return_paths_dict(run_with, file_suffix="parquet", file_midfix="compressedOnGenes_simple")

    cutoff = 25

    prefix_dict = {1: "read_hits",
                   2: "hits_rank",
                   3: "read_len_mean",
                   4: "polya_mean"}
    
    prefix = prefix_dict[1]
    
    color_by_dict = {0: None,
                     1: "read_len_mean_mean",
                     2: "tail_length_diff",
                     3: "read_len_std_mean",
                     4: "tail_avg_stdev",
                     5: "chr",
                     6: "gene_gc",
                     7: "is_MtDNA",
                     8: "a_fraction",
                     9: "t_fraction",
                     10: "tail_length_mean_mean"
                     }
    
    color_by = color_by_dict[9]
    
    drop_mtDNA = False
    drop_gcLess = False
    trend_line = True
    save_directory = f"/home/marcus/Documents/{get_dt(for_file=True)}_plotlyFigures"
    
    if len(pathdict.keys()) == 3:
        for to_drop in pathdict.keys():
            plotter_helper(pathdict, prefix, cutoff, one_to_drop=to_drop,
                           color_by=color_by,
                           drop_mt_dna=drop_mtDNA,
                           drop_gc_less=drop_gcLess,
                           ols_trend=trend_line,
                           save_dir=save_directory,
                           )
    else:
        plotter_helper(pathdict, prefix, cutoff,
                       color_by=color_by,
                       drop_mt_dna=drop_mtDNA,
                       drop_gc_less=drop_gcLess,
                       ols_trend=trend_line,
                       save_dir=save_directory,
                       )
