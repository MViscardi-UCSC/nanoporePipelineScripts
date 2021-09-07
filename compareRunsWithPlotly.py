"""
compareRunsWithPlotly.py
Marcus Viscardi, 8/3/21


"""


import pandas as pd
import scipy.stats
from step0_nanopore_pipeline import find_newest_matching_file

import warnings
from pandas.core.common import SettingWithCopyWarning

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
pd.set_option("display.max_columns", None)


def load_df(path_to_file) -> pd.DataFrame:
    with open(path_to_file, 'r') as file:
        df = pd.read_csv(file, sep="\t")
    return df


def load_and_merge(file_paths: (str, str), min_hit_cutoff: int = None) -> pd.DataFrame:
    df1, df2 = [load_df(path) for path in file_paths]
    for df in (df1, df2):
        df["hits_rank"] = df['read_hits'].rank(ascending=False)
    print(f"df1 length: {df1.shape[0]}\ndf2 length: {df2.shape[0]}")
    merge_df = df1.merge(df2, how="inner", on="gene_id", suffixes=["_totalRNA", "_polyA"])
    merge_df["mean_tail_diff"] = merge_df['polya_mean_totalRNA'].sub(merge_df['polya_mean_polyA'], axis=0).abs()
    print(f"merge df length: {merge_df.shape[0]}")
    if min_hit_cutoff:
        cutoff_merge = merge_df[merge_df["read_hits_totalRNA"] >= min_hit_cutoff]
        cutoff_merge = cutoff_merge[cutoff_merge["read_hits_polyA"] >= min_hit_cutoff]
        print(f"after dropping genes w/ <{min_hit_cutoff} hits: {cutoff_merge.shape[0]}")
        return cutoff_merge
    else:
        return merge_df


def load_and_merge_from_dict(path_dict: dict, min_hit_cutoff: int = None,
                             add_names: bool = True) -> (pd.DataFrame, list):
    df_dict = {key: load_df(path) for (key, path) in path_dict.items()}
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
        path_to_names = "/home/marcus/PycharmProjects/PersonalScripts/WBGene_to_geneName.tsv"
        names_df = pd.read_csv(path_to_names, sep="\t")
        merge_df = merge_df.merge(names_df, on="gene_id", how="left")
        merge_df['ident'] = merge_df["gene_name"] + " (" + merge_df["gene_id"] + ")"
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
                             x=0, y=1, compare_column_prefix="polya_mean"):
    import plotly.express as px
    import plotly.graph_objects as go
    if isinstance(x, int) and isinstance(y, int):
        x_key = key_list[x]
        y_key = key_list[y]
    elif isinstance(x, str) and isinstance(y, str):
        x_key = x
        y_key = y

    if compare_column_prefix == "polya_mean":
        title = f"PolyA Mean Lengths per Gene between RNA Selection Methods"
    elif compare_column_prefix == "read_hits":
        title = f"Number of Reads per Gene between RNA Selection Methods"
    else:
        title = f"Comparing: {compare_column_prefix} between {x_key} and {y_key}"

    if compare_column_prefix == "polya_mean":
        # merged_df["read_rank_diff"] = merged_df[f'hits_rank_{x_key}']. \
        #     sub(merged_df[f'hits_rank_{y_key}'], axis=0).abs()
        # color_column = "read_rank_diff"
        merged_df["tail_avg_stdev"] = merged_df[f'polya_stdev_{x_key}']. \
            add(merged_df[f'polya_stdev_{y_key}'], axis=0).abs().div(2)
        color_column = "tail_avg_stdev"
    elif compare_column_prefix == "read_len_mean":
        merged_df["read_len_std_mean"] = merged_df[f'read_len_std_{x_key}']. \
            add(merged_df[f'read_len_std_{y_key}'], axis=0).abs().div(2)
        color_column = "read_len_std_mean"
    else:
        merged_df["tail_length_diff"] = merged_df[f'polya_mean_{x_key}']. \
            sub(merged_df[f'polya_mean_{y_key}'], axis=0).abs()
        color_column = "tail_length_diff"
    labels_dict = {"gene_id": "WBGene ID",
                   "ident": "Identity",
                   "tail_length_diff": "Mean Tail Length Diff (nts)",
                   "read_rank_diff": "Difference in Read Rank"}
    for key in [x_key, y_key]:
        labels_dict[f"polya_mean_{key}"] = f"{key} - Mean Tail Length (nts)"
        labels_dict[f"read_hits_{key}"] = f"{key} - Reads/Gene"
        labels_dict[f"hits_rank_{key}"] = f"{key} - Rank of Reads/Gene"
        labels_dict[f"read_len_mean_{key}"] = f"{key} - Mean Read Length/Gene"
    fig = px.scatter(merged_df, x=f"{compare_column_prefix}_{x_key}", y=f"{compare_column_prefix}_{y_key}",
                     color=color_column,
                     opacity=.8,
                     hover_name='ident',
                     hover_data=[f"polya_mean_{x_key}", f"polya_mean_{y_key}",
                                 f'read_hits_{x_key}', f'read_hits_{y_key}'],
                     title=title,
                     labels=labels_dict,
                     trendline="ols",
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

    fig.show()


def plotly_scatter(data, cutoff=None):
    import plotly.express as px

    fig = px.scatter(data, x="polya_mean_totalRNA", y="polya_mean_polyA",
                     color="hits_rank_totalRNA",
                     # color="mean_tail_diff",
                     hover_name='gene_id',
                     hover_data=['polya_mean_totalRNA', 'polya_mean_polyA', 'read_hits_totalRNA', 'read_hits_polyA'],
                     title=f"PolyA Mean Lengths per Gene between RNA Selection Methods\n(dropped genes w/ <{cutoff} hits)",
                     labels=dict(polya_mean_totalRNA="Total RNA - Mean Tail Length (nts)",
                                 polya_mean_polyA="PolyA - Mean Tail Length (nts)",
                                 read_hits_totalRNA="Total RNA - Reads/Gene",
                                 read_hits_polyA="PolyA - Reads/Gene",
                                 gene_id="WBGene ID",
                                 mean_tail_diff="Mean Tail Length Difference (nts)"),
                     trendline="ols")

    fig.update_layout(shapes=[{'type': 'line',
                               'yref': 'paper',
                               'xref': 'paper',
                               'y0': 0, 'y1': 1,
                               'x0': 0, 'x1': 1}])
    min = 30
    max = 120
    fig.update_xaxes(range=[min, max])
    fig.update_yaxes(range=[min, max])

    fig.show()


def plotly_w_slider(data, max_cutoff=50):
    """
    Grabbed from here: https://plotly.com/python/sliders/#sliders-in-plotly-express
    
    Not working?
    
    :param data: 
    :param max_cutoff: 
    :return: 
    """
    import plotly.graph_objects as go
    import plotly.express as px
    import numpy as np

    # Create figure
    fig = go.Figure()

    # Add traces, one for each slider step
    for step in np.arange(1, max_cutoff + 1, 1):
        step_data = data[data["read_hits_totalRNA"] >= step]
        fig.add_scatter(visible=False,
                        mode="markers",
                        name="cutoff = " + str(step),
                        hovertext=step_data["gene_id"],
                        x=step_data["polya_mean_totalRNA"], y=step_data["polya_mean_polyA"], )
        fig.update_layout(shapes=[{'type': 'line',
                                   'yref': 'paper',
                                   'xref': 'paper',
                                   'y0': 0, 'y1': 1,
                                   'x0': 0, 'x1': 1}])
    min = 30
    max = 120
    fig.update_xaxes(range=[min, max])
    fig.update_yaxes(range=[min, max])

    # Make 20th trace visible
    fig.data[20].visible = True

    # Create and add slider
    steps = []
    for i in range(len(fig.data)):
        step = dict(
            method="update",
            args=[{"visible": [False] * len(fig.data)},
                  {"title": f"Slider set to gene hit cutoff of: {i} reads"}],  # layout attribute
        )
        step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
        steps.append(step)

    sliders = [dict(
        active=20,
        currentvalue={"prefix": "Gene Hit Cutoff: "},
        pad={"t": 50},
        steps=steps
    )]

    fig.update_layout(
        sliders=sliders
    )

    fig.show()


def plotter_helper(path_dict: dict, prefix: str, cut_off: int, one_to_drop: str = ""):
    key_list = list(path_dict.keys())
    if one_to_drop:
        subset_list = [i for i in key_list if i != one_to_drop]
        subset_dict = {key: path_dict[key] for key in subset_list}
    else:
        subset_dict = path_dict
    merge_df, keys = load_and_merge_from_dict(subset_dict, min_hit_cutoff=cutoff)
    print(merge_df.columns)
    plotly_from_triple_merge(merge_df, keys,
                             cutoff=cut_off,
                             x=0,
                             y=1,
                             compare_column_prefix=prefix,
                             )


def plotly_3d(pathdict, prefix, cutoff):
    merge_df, key_list = load_and_merge_from_dict(pathdict, min_hit_cutoff=cutoff)
    import plotly.express as px
    merge_df["sum_tail_length"] = \
        merge_df[f'polya_mean_{key_list[0]}'].add(merge_df[f'polya_mean_{key_list[1]}'], axis=0).add(
            merge_df[f'polya_mean_{key_list[2]}'], axis=0)
    print(f"{prefix}_{list(pathdict.keys())[0]}")
    fig = px.scatter_3d(merge_df,
                        x=f"{prefix}_{key_list[0]}",
                        y=f"{prefix}_{key_list[1]}",
                        z=f"{prefix}_{key_list[2]}",
                        hover_name='gene_id', color="sum_tail_length",
                        hover_data=[f'polya_mean_{key_list[0]}', f'polya_mean_{key_list[1]}',
                                    f'polya_mean_{key_list[2]}',
                                    f'read_hits_{key_list[0]}', f'read_hits_{key_list[1]}',
                                    f'read_hits_{key_list[2]}'],
                        labels={f"polya_mean_{key_list[0]}": f"{key_list[0]} - Mean Tail Length (nts)",
                                f"polya_mean_{key_list[1]}": f"{key_list[1]} - Mean Tail Length (nts)",
                                f"polya_mean_{key_list[2]}": f"{key_list[2]} - Mean Tail Length (nts)",
                                f"read_hits_{key_list[0]}": f"{key_list[0]} - Reads/Gene",
                                f"read_hits_{key_list[1]}": f"{key_list[1]} - Reads/Gene",
                                f"read_hits_{key_list[2]}": f"{key_list[2]} - Reads/Gene",
                                f"hits_rank_{key_list[0]}": f"{key_list[0]} - Rank of Reads/Gene",
                                f"hits_rank_{key_list[1]}": f"{key_list[1]} - Rank of Reads/Gene",
                                f"hits_rank_{key_list[2]}": f"{key_list[2]} - Rank of Reads/Gene",
                                f"gene_id": "WBGene ID",
                                }, )
    fig.show()


def plotly_subplots_and_load(pathdict, prefix, cutoff):
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go
    import numpy as np

    merge_df, keys = load_and_merge_from_dict(pathdict, min_hit_cutoff=cutoff)

    fig = make_subplots(rows=1, cols=len(keys))

    label_dict = {f"gene_id": "WBGene ID"}
    if prefix == "polya_mean" or prefix == "hits_rank":
        scale_type = "linear"
    elif prefix == "read_hits":
        scale_type = "log"
    else:
        raise ValueError(f"The prefix {prefix} isn't coded into here... yet!")

    for key in keys:
        label_dict[f"polya_mean_{key}"] = f"{key} - Mean Tail Length (nts)"
        label_dict[f"read_hits_{key}"] = f"{key} - Reads/Gene"
        label_dict[f"hits_rank_{key}"] = f"{key} - Ranking of Reads/Gene"

    for plot_num in range(0, len(keys), 1):

        compare_plot = plot_num + 1
        if compare_plot >= len(keys):
            compare_plot = 0
        print(f"Grabbing plots from key numbers {plot_num} & {compare_plot}")

        plot_name = f"{prefix}_{keys[plot_num]}"
        comp_plot_name = f"{prefix}_{keys[compare_plot]}"
        row_num = 1
        col_num = plot_num + 1
        fig.add_trace(go.Scatter(x=merge_df[plot_name], y=merge_df[comp_plot_name],
                                 mode='markers', name=f"{keys[plot_num]} vs. {keys[compare_plot]}"),
                      row=row_num, col=col_num)
        fig.add_trace(go.Scatter(x=np.arange(0, 1000, 10), y=np.arange(0, 1000, 10),
                                 line=dict(color="#4a4a4a", width=1, dash='dash'),
                                 showlegend=False, mode="lines"),
                      row=row_num, col=col_num)
        fig.update_xaxes(title_text=label_dict[plot_name],
                         row=row_num, col=col_num, type=scale_type)
        fig.update_yaxes(title_text=label_dict[comp_plot_name],
                         row=row_num, col=col_num, type=scale_type)
    fig.update_layout(height=400, width=1000, showlegend=False,
                      title_text=f"Number of Reads per Gene between RNA Selection Methods")
    fig.show()


def plotly_click_to_clipboard():
    import plotly.graph_objects as go
    import numpy as np
    np.random.seed(1)

    x = np.random.rand(100)
    y = np.random.rand(100)

    f = go.FigureWidget([go.Scatter(x=x, y=y, mode='markers')])

    scatter = f.data[0]
    colors = ['#a3a7e4'] * 100
    scatter.marker.color = colors
    scatter.marker.size = [10] * 100
    f.layout.hovermode = 'closest'

    # create our callback function
    def update_point(trace, points, selector):
        c = list(scatter.marker.color)
        s = list(scatter.marker.size)
        for i in points.point_inds:
            c[i] = '#bae2be'
            s[i] = 20
            with f.batch_update():
                scatter.marker.color = c
                scatter.marker.size = s

    scatter.on_click(update_point)

    f.show()


if __name__ == '__main__':
    # path1, path2 = sys.argv

    pathdict = {
        # "riboD": "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/output_dir/"
        #          "merge_files/210713_compressedOnGenes_simple.tsv",
        # "totalRNA": find_newest_matching_file("/data16/marcus/working/210709_NanoporeRun_totalRNA_0639_L3/"
        #                                       "output_dir/merge_files/*_compressedOnGenes_simple.tsv"),
        "totalRNA2": find_newest_matching_file("/data16/marcus/working/"
                                               "210720_nanoporeRun_totalRNA_0639_L3_replicate/output_dir/"
                                               "merge_files/*_compressedOnGenes_simple.tsv"),
        # "polyA": find_newest_matching_file("/data16/marcus/working/210528_NanoporeRun_0639_L3s/output_dir/"
        #                                    "merge_files/*_compressedOnGenes_simple.tsv"),
        # "polyA2": find_newest_matching_file("/data16/marcus/working/210719_nanoporeRun_polyA_0639_L3_replicate/"
        #                                     "output_dir/merge_files/*_compressedOnGenes_simple.tsv"),
        "xrn-1": find_newest_matching_file("/data16/marcus/working/210905_nanoporeRun_totalRNA_5108_xrn-1-KD/"
                                           "output_dir/merge_files/*_compressedOnGenes_simple.tsv")
    }

    cutoff = 20
    prefix = "hits_rank"
    # prefix = "read_hits"
    # prefix = "polya_mean"
    # prefix = "read_len_mean"
    if len(pathdict.keys()) == 3:
        for to_drop in pathdict.keys():
            plotter_helper(pathdict, prefix, cutoff, one_to_drop=to_drop)
    else:
        plotter_helper(pathdict, prefix, cutoff)
    
    # plotly_click_to_clipboard()
