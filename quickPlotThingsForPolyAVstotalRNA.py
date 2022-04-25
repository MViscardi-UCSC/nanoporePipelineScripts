"""
quickPlotThingsForPolyAVstotalRNA.py
Marcus Viscardi,    February 09, 2022

Place to test out unique plots that I want to have for
the polyA vs totalRNA paper!
"""
import pandas as pd
import numpy as np
import plotly.graph_objs
from scipy import stats

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)

import seaborn as sea
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

from nanoporePipelineCommon import \
    library_reads_df_load_and_concat, \
    compress_concat_df_of_libs, \
    load_and_merge_lib_parquets, \
    find_newest_matching_file, \
    get_dt

OUTPUT_DIR = "/home/marcus/Insync/mviscard@ucsc.edu/" \
             "Google Drive/insync_folder/polyAPaperFigures/" \
             "figure_foldChangeAndTails/raw"


def compare_libraries_w_pA_over_total(force_compressed_df_build=False,
                                      cutoff=40,
                                      short_tail_cutoff=0.57,
                                      long_tail_cutoff=0.26,
                                      size_by_expression_quartile=True,
                                      save_files=True,
                                      locked_aspect=True,
                                      color_by_bins=False):
    libs_to_compare = ["polyA2", "totalRNA2",
                       "polyA3", "totalRNA3"]

    if not force_compressed_df_build:
        compressed_save_search = f"./testInputs/*_{'_'.join(libs_to_compare)}.compressed.parquet"
        reads_save_search = f"./testInputs/*_{'_'.join(libs_to_compare)}.reads.parquet"
        print(f"Looking for pre-processed file at: {compressed_save_search}...")
        try:
            compressed_parquet_path = find_newest_matching_file(compressed_save_search)
            compressed_df = pd.read_parquet(compressed_parquet_path)
            reads_save_path = find_newest_matching_file(reads_save_search)
            reads_df = pd.read_parquet(reads_save_path)
        except ValueError:
            print(f"Couldn't find pre-processed file at: {compressed_save_search}\nGoing to load from library files!")
            force_compressed_df_build = True
    if force_compressed_df_build:
        reads_df, compressed_df = load_and_merge_lib_parquets(libs_to_compare, drop_sub_n=1)
        compressed_save_path = f"./testInputs/{get_dt(for_file=True)}_{'_'.join(libs_to_compare)}.compressed.parquet"
        reads_save_path = f"./testInputs/{get_dt(for_file=True)}_{'_'.join(libs_to_compare)}.reads.parquet"
        compressed_df.to_parquet(compressed_save_path)
        reads_df.to_parquet(reads_save_path)
        print(f"Saved new compressed file to: {compressed_save_path}")
        print(f"Saved new reads file to:      {reads_save_path}")
    compressed_df = compressed_df.query(f'gene_hits >= {cutoff}')

    compressed_df.loc[:, 'tail_groupings_group'] = pd.cut(compressed_df[f'frac_sub_50_tails'],
                                                          bins=[0.0, long_tail_cutoff, short_tail_cutoff, 1.0],
                                                          labels=['long_tailed',
                                                                  'ungrouped',
                                                                  'short_tailed'],
                                                          include_lowest=True)

    compressed_df.loc[:, "lib_set"] = compressed_df.lib.str[-1].copy()
    compressed_df.loc[:, "lib_type"] = compressed_df.lib.str[:-1].copy()
    compressed_df.loc[:, 'mean_median_diff_polya_length'] = compressed_df["mean_polya_length"] - \
                                                            compressed_df["median_polya_length"]
    multi_df = compressed_df.set_index(["lib_set", "chr_id", "gene_id", "gene_name", "lib_type"]).sort_index().copy()
    group_2_df = multi_df.query("lib_set == '2'").copy()
    group_3_df = multi_df.query("lib_set == '3'").copy()
    retain_cols = ["chr_id",
                   "gene_id",
                   "gene_name",
                   "lib",
                   "gene_hits",
                   "gene_rpm",
                   "mean_polya_length",
                   "median_polya_length",
                   "mean_median_diff_polya_length",
                   "mean_read_length",
                   "tail_groupings_group"]
    dfs = {}
    for df in [group_2_df, group_3_df]:
        df1 = df.query("lib_type == 'polyA'").reset_index()[retain_cols].copy()
        df2 = df.query("lib_type == 'totalRNA'").reset_index()[retain_cols].copy()
        df_merge = df1.merge(df2, on=["chr_id", "gene_id", "gene_name"], suffixes=["_polyA", "_totalRNA"])
        lib_set = df.index[0][0]
        dfs[lib_set] = df_merge
    for lib_set, df in dfs.items():
        df[f"pA_over_total_rpm"] = df.gene_rpm_polyA / df.gene_rpm_totalRNA
        df[f"pA_over_total_mean_median_diff_polya_lengths"] = df["mean_median_diff_polya_length_polyA"] / \
                                                              df["mean_median_diff_polya_length_totalRNA"]
        df[f"avg_mean_median_diff_polya_lengths"] = (df["mean_median_diff_polya_length_polyA"] +
                                                     df["mean_median_diff_polya_length_totalRNA"]) / 2
        df[f"pA_minus_total_mean_diff_polya_lengths"] = df["mean_polya_length_polyA"] - \
                                                        df["mean_polya_length_totalRNA"]
        df[f"pA_minus_total_median_diff_polya_lengths"] = df["median_polya_length_polyA"] - \
                                                          df["median_polya_length_totalRNA"]
        df[f"gene_rpm_sum"] = df['gene_rpm_totalRNA'] + df['gene_rpm_polyA']
        df["gene_rpm_mean"] = df['gene_rpm_sum'] / 2
    dfs_list = list(dfs.values())
    final_merge_cols = ["chr_id",
                        "gene_id",
                        "gene_name",
                        "gene_rpm_polyA",
                        "gene_rpm_totalRNA",
                        "gene_rpm_sum",
                        "gene_rpm_mean",
                        "pA_over_total_rpm",
                        "tail_groupings_group_totalRNA",
                        "pA_over_total_mean_median_diff_polya_lengths",
                        "avg_mean_median_diff_polya_lengths",
                        "pA_minus_total_mean_diff_polya_lengths",
                        "pA_minus_total_median_diff_polya_lengths",
                        "mean_median_diff_polya_length_polyA",
                        "mean_median_diff_polya_length_totalRNA",
                        "mean_polya_length_polyA",
                        "mean_polya_length_totalRNA",
                        ]
    plot_df = pd.merge(dfs_list[0][final_merge_cols],
                       dfs_list[1][final_merge_cols],
                       on=["chr_id", "gene_name", "gene_id"],
                       suffixes=("_2", "_3"))
    plot_df['pA2_over_pA3_rpm'] = plot_df['gene_rpm_polyA_2'] / plot_df['gene_rpm_polyA_3']
    plot_df['log2_pA2_over_pA3_rpm'] = np.log2(plot_df['pA2_over_pA3_rpm'])
    plot_df['tail_group'] = plot_df[plot_df.tail_groupings_group_totalRNA_2 ==
                                    plot_df.tail_groupings_group_totalRNA_3].tail_groupings_group_totalRNA_2
    plot_df["tail_group"].fillna("ungrouped", inplace=True)
    for lib_set in [2, 3]:
        plot_df[f'log2_pA_over_total_rpm_{lib_set}'] = np.log2(plot_df[f'pA_over_total_rpm_{lib_set}'])
        # Plot mean_polyA - mean_tot for each lib against eachother:
        print(plot_df.columns)
        # Plot scatter of mean-median, colored by tail group
        # fig = px.scatter(plot_df,
        #                  x=f"mean_median_diff_polya_length_polyA_{lib_set}",
        #                  y=f"mean_median_diff_polya_length_totalRNA_{lib_set}",
        #                  color='tail_group',
        #                  hover_name="gene_name",
        #                  hover_data=["gene_rpm_polyA_2", "gene_rpm_polyA_3",
        #                              "gene_rpm_totalRNA_2", "gene_rpm_totalRNA_3",
        #                              ])
        # fig.add_trace(go.Scatter(x=[-20, 50],
        #                          y=[-20, 50],
        #                          mode='lines',
        #                          line=dict(color='black',
        #                                    dash='dash'),
        #                          showlegend=False))
        # fig.update_layout(template='plotly_white')
        # fig.show()
    # Plot mean-median average
    plot_df['overall_avg_mean_median_diff_polya_len'] = (plot_df[f"avg_mean_median_diff_polya_lengths_2"] +
                                                         plot_df[f"avg_mean_median_diff_polya_lengths_3"]) / 2
    if color_by_bins:
        if isinstance(color_by_bins, int):
            num_bins = color_by_bins
        else:
            num_bins = 10
        plot_df['binned_overall_avg_mean_median_diff_polya_len'] = pd.qcut(
            plot_df['overall_avg_mean_median_diff_polya_len'],
            num_bins,
            labels=False) + 1
        # fig = px.scatter(plot_df,
        #                  x=f"avg_mean_median_diff_polya_lengths_2",
        #                  y=f"avg_mean_median_diff_polya_lengths_3",
        #                  color='binned_overall_avg_mean_median_diff_polya_len',
        #                  color_continuous_scale=px.colors.sample_colorscale(px.colors.diverging.Portland,
        #                                                                     [(x / num_bins) for x in
        #                                                                      range(1, num_bins + 1)]),
        #                  hover_name="gene_name",
        #                  hover_data=["gene_rpm_polyA_2", "gene_rpm_polyA_3",
        #                              "gene_rpm_totalRNA_2", "gene_rpm_totalRNA_3",
        #                              ])
        # fig.add_trace(go.Scatter(x=[-20, 50],
        #                          y=[-20, 50],
        #                          mode='lines',
        #                          line=dict(color='black',
        #                                    dash='dash'),
        #                          showlegend=False))
        # fig.update_layout(template='plotly_white')
        # fig.show()
    plot_df["is_MtDNA"] = plot_df['chr_id'] == 'MtDNA'
    plot_df['gene_rpm_sum'] = plot_df['gene_rpm_sum_2'] + plot_df['gene_rpm_sum_3']
    plot_df['log_gene_rpm_sum'] = np.log2(plot_df['gene_rpm_sum'])
    plot_df['gene_rpm_mean'] = (plot_df['gene_rpm_mean_2'] + plot_df['gene_rpm_mean_3']) / 2
    plot_df['log_gene_rpm_mean'] = np.log2(plot_df['gene_rpm_mean'])

    plot_df['gene_rpm_binned'] = pd.qcut(plot_df['gene_rpm_mean'],
                                         4,
                                         labels=False) + 1

    # Labels dict used in plotly charts to make things a bit more readable
    labels_dict = {"log2_pA_over_total_rpm_3": "Log<sub>2</sub>(pA/total rpm); set 3",
                   "log2_pA_over_total_rpm_2": "Log<sub>2</sub>(pA/total rpm); set 2",
                   "gene_name": "Gene Name",
                   "is_MtDNA": "Mitochondrial RNA",
                   }
    for mean_or_median in [
        'mean',
        # 'median',
    ]:
        # I am creating these labels before even making the columns they refer to, little backwards,
        #   but I am trying to keep all the label making relatively close together!
        labels_dict[f"binned_diff_of_{mean_or_median}_tail_length"] = f"Binned diff of {mean_or_median}" \
                                                                      f"<br>tail lengths" \
                                                                      f"<br>between techs"
        labels_dict[f"avg_diff_of_{mean_or_median}_tail_length"] = f"Avg {mean_or_median} tail length" \
                                                                   f"<br>diff between techs"

        plot_df[f"avg_diff_of_{mean_or_median}_tail_length"] = (plot_df[
                                                                    f"pA_minus_total_{mean_or_median}"
                                                                    f"_diff_polya_lengths_2"] +
                                                                plot_df[
                                                                    f"pA_minus_total_{mean_or_median}"
                                                                    f"_diff_polya_lengths_3"]) / 2
        if color_by_bins:
            if isinstance(color_by_bins, int):
                num_bins = color_by_bins
            else:
                num_bins = 10
            plot_df[f'binned_diff_of_{mean_or_median}_tail_length'], bins = pd.qcut(
                plot_df[f'avg_diff_of_{mean_or_median}_tail_length'],
                num_bins,
                labels=False,
                retbins=True)
            plot_df[f'binned_diff_of_{mean_or_median}_tail_length'] += 1
            print(f"{mean_or_median.title()} Bins:\n", bins)
            # Plot mean_polyA - mean_tot for each lib against each other,
            #   colored by their place along y=x:
            color_col = None
            # Comment out below line for black points
            color_col = f'binned_diff_of_{mean_or_median}_tail_length'
            fig = px.scatter(plot_df,
                             x=f"pA_minus_total_{mean_or_median}_diff_polya_lengths_3",
                             y=f"pA_minus_total_{mean_or_median}_diff_polya_lengths_2",
                             color=color_col,
                             color_continuous_scale='Portland',
                             color_discrete_sequence=px.colors.sample_colorscale(px.colors.diverging.Portland,
                                                                                 [(x / num_bins) for x in
                                                                                  range(1, num_bins + 1)]),
                             hover_name="gene_name",
                             hover_data=["gene_rpm_polyA_2", "gene_rpm_polyA_3",
                                         "gene_rpm_totalRNA_2", "gene_rpm_totalRNA_3",
                                         ],
                             height=500,
                             width=500)
            fig.update_coloraxes(showscale=False)
            fig.update_xaxes(dtick=10)
            fig.update_yaxes(dtick=10)
            marker_dict = dict(line=dict(width=1, color='DarkSlateGrey'))
            marker_dict['size'] = 7
            if color_col is None:
                marker_dict['color'] = 'black'
            fig.update_traces(marker=marker_dict, selector=dict(mode='markers'))
            fig.add_trace(go.Scatter(x=[-10, 40],
                                     y=[-10, 40],
                                     mode='lines',
                                     line=dict(color='black',
                                               dash='dash'),
                                     showlegend=False))
            fig.update_layout(template='plotly_white')
            fig.show()
            if save_files:
                output_name = f"{OUTPUT_DIR}/{get_dt(for_file=True)}_"
                if color_by_bins:
                    output_name += "binned-tail-change-colors_"
                output_name += f"changeIn{mean_or_median.title()}TailScatter_{cutoff}rpg"
                fig.write_image(output_name + ".png")
                fig.write_image(output_name + ".svg")

    # Plot the 'fold change' in rpm between techniques between libs:
    for lib_set in [2, 3]:
        for lib_type in ['polyA', 'totalRNA']:
            labels_dict[f"gene_rpm_{lib_type}_{lib_set}"] = f"Gene RPM ({lib_type}{lib_set})"
        if color_by_bins:
            # Plot cdfs of rpm 'fold change' per technique, color by grouping!
            plot_df[f'binned_pA_minus_total_mean_diff_polya_lengths_{lib_set}'] = pd.qcut(
                plot_df[f"pA_minus_total_mean_diff_polya_lengths_{lib_set}"],
                num_bins,
                labels=False) + 1
            labels_dict[f'binned_pA_minus_total_mean_diff_polya_lengths_{lib_set}'] = f"Binned diff of mean" \
                                                                                      f"<br>tail lengths" \
                                                                                      f"<br>between techs (set{lib_set})"
        # fig = px.ecdf(plot_df.sort_values(f'binned_pA_minus_total_mean_diff_polya_lengths_{lib_set}'),
        #               x=f"log2_pA_over_total_rpm_{lib_set}",
        #               color=f'binned_pA_minus_total_mean_diff_polya_lengths_{lib_set}',
        #               color_discrete_sequence=px.colors.sample_colorscale(px.colors.diverging.Portland,
        #                                                                   [(x/num_bins) for x in range(1, num_bins+1)]),
        #               hover_name="gene_name",
        #               hover_data=["gene_rpm_polyA_2", "gene_rpm_polyA_3",
        #                           "gene_rpm_totalRNA_2", "gene_rpm_totalRNA_3",
        #                           ],
        #               labels=labels_dict)
        # fig.update_layout(template='plotly_white')
        # fig.show()
        # # Plot just the polyA duplicates:
        # fig = px.ecdf(plot_df.sort_values(f'binned_diff_of_mean_tail_length'),
        #               x=f"log2_pA2_over_pA3_rpm",
        #               color=f'binned_diff_of_mean_tail_length',
        #               color_discrete_sequence=px.colors.sample_colorscale(px.colors.diverging.Portland,
        #                                                                   [(x/num_bins) for x in range(1, num_bins+1)]),
        #               hover_name="gene_name",
        #               hover_data=["gene_rpm_polyA_2", "gene_rpm_polyA_3",
        #                           "gene_rpm_totalRNA_2", "gene_rpm_totalRNA_3",
        #                           ],
        #               labels=labels_dict)
        # fig.update_layout(template='plotly_white')
        # fig.show()

    if size_by_expression_quartile:
        size_col = "gene_rpm_binned"
    else:
        size_col = None

    fig = plot_rpm_fold_change(plot_df, labels_dict, cutoff, color_by_bins,
                               locked_aspect, size_col, size_by_expression_quartile,
                               mean_or_median="mean", show_color_legend=False, show_annotation=False)

    if save_files:
        output_name = f"{OUTPUT_DIR}/{get_dt(for_file=True)}_"
        if size_col:
            output_name += "mean-rpm-sized_"
        if color_by_bins:
            output_name += "binned-tail-change-colors_"
        output_name += f"techFoldChangeScatter_{cutoff}rpg"
        fig.write_html(output_name + ".html")
        fig.write_image(output_name + ".png", scale=10)
        fig.write_image(output_name + ".svg")

        # TODO: Better use of this:
        save_df = plot_df[
            ['chr_id', 'gene_id', 'gene_name', 'gene_rpm_polyA_2', 'mean_polya_length_polyA_2', 'gene_rpm_totalRNA_2',
             'mean_polya_length_totalRNA_2', 'pA_minus_total_mean_diff_polya_lengths_2', 'gene_rpm_polyA_3',
             'mean_polya_length_polyA_3', 'gene_rpm_totalRNA_3', 'mean_polya_length_totalRNA_3',
             'pA_minus_total_mean_diff_polya_lengths_3', 'avg_diff_of_mean_tail_length',
             'binned_diff_of_mean_tail_length']]
        save_df_column_dict = {
            # General Things (unchanged):
            'chr_id': 'chr_id',
            'gene_id': 'gene_id',
            'gene_name': 'gene_name',
            # Lib set 3 (converted to lib replicate 1 for paper):
            # Per lib things:
            "gene_rpm_polyA_3": "gene_rpm__selected1",
            "mean_polya_length_polyA_3": "mean_tail_length__selected1",
            "gene_rpm_totalRNA_3": "gene_rpm__unselected1",
            "mean_polya_length_totalRNA_3": "mean_tail_length__unselected1",
            # Per set things:
            "pA_minus_total_mean_diff_polya_lengths_3": "diff_in_mean_tail_lengths__selected1_minus_unselected1",
            # Lib set 2, stays as 2
            # Per lib things:
            "gene_rpm_polyA_2": "gene_rpm__selected2",
            "mean_polya_length_polyA_2": "mean_tail_length__selected2",
            "gene_rpm_totalRNA_2": "gene_rpm__unselected2",
            "mean_polya_length_totalRNA_2": "mean_tail_length__unselected2",
            # Per set things:
            "pA_minus_total_mean_diff_polya_lengths_2": "diff_in_mean_tail_lengths__selected2_minus_unselected2",
            # All Lib Things:
            #################
            "avg_diff_of_mean_tail_length": "avg_diff_in_mean_tail_lengths__set1_and_set2",
            "binned_diff_of_mean_tail_length": "deciles_of_avg_diff_in_mean_tail_lengths__set1_and_set2",
        }
        save_df_new_col_order = list(save_df_column_dict.values())
        save_df.rename(columns=save_df_column_dict, inplace=True)
        save_df = save_df[save_df_new_col_order]
        # Breakpoint below and change the save_csv param to True:
        save_csv = False
        print("Breakpoint.")
        if save_csv:
            save_df.to_csv("/home/marcus/Documents/figure4_supplementalTable.csv")

    # Plot plot gene rpm in polya against mean-median/mean-median
    # for lib_set in [2, 3]:
    #     fig = px.scatter(plot_df,
    #                      # x=f"gene_rpm_polyA_{lib_set}",
    #                      y=f'pA_over_total_mean_median_diff_polya_lengths_{lib_set}',
    #                      # y=f"gene_rpm_totalRNA_{lib_set}",
    #                      color=f'pA_over_total_mean_median_diff_polya_lengths_{lib_set}',
    #                      hover_name="gene_name",
    #                      hover_data={"gene_rpm_polyA_2", "gene_rpm_polyA_3",
    #                                  "gene_rpm_totalRNA_2", "gene_rpm_totalRNA_3",
    #                                  }
    #                      )
    #     fig.update_xaxes(type='log')
    #     fig.update_yaxes(type='log')
    #     fig.show()
    print(stats.spearmanr(plot_df.log2_pA_over_total_rpm_2, plot_df.log2_pA_over_total_rpm_3))

    print("breakpoint")
    return [compressed_df, reads_df]


def plot_rpm_fold_change(plot_df: pd.DataFrame, labels_dict: dict, cutoff: int,
                         color_by_bins: bool, locked_aspect: bool, size_col: str,
                         size_by_expression_quartile: bool, mean_or_median: str = "mean",
                         save_files=False, show_color_legend=True, show_annotation=True,
                         ) -> plotly.graph_objs.Figure:
    if locked_aspect:
        width_var, height_var = 500, 500
    else:
        width_var, height_var = None, None
    if color_by_bins:
        color_col = f'binned_diff_of_{mean_or_median}_tail_length'
    else:
        # color_col = f'avg_diff_of_{mean_or_median}_tail_length'
        color_col = None
    fig = px.scatter(plot_df,
                     x="log2_pA_over_total_rpm_3",
                     y="log2_pA_over_total_rpm_2",
                     color=color_col,
                     # color='binned_overall_avg_mean_median_diff_polya_len',
                     # color="tail_group",
                     symbol='is_MtDNA',
                     symbol_sequence=[0, 3],
                     size=size_col, size_max=6,
                     opacity=0.9,
                     color_continuous_scale=px.colors.diverging.Portland,
                     hover_name="gene_name",
                     hover_data=["gene_rpm_polyA_2", "gene_rpm_polyA_3",
                                 "gene_rpm_totalRNA_2", "gene_rpm_totalRNA_3",
                                 ],
                     labels=labels_dict,
                     width=width_var, height=height_var,
                     )
    marker_dict = dict(line=dict(width=1, color='DarkSlateGrey'))
    if size_col is None:
        marker_dict['size'] = 7
    if color_col is None:
        marker_dict['color'] = 'black'
    fig.update_traces(marker=marker_dict, selector=dict(mode='markers'))
    fig.add_trace(go.Scatter(x=[-20, 20],
                             y=[-20, 20],
                             mode='lines',
                             line=dict(color='black',
                                       dash='dash'),
                             showlegend=False))
    fig.add_trace(go.Scatter(x=[-20, 20],
                             y=[-0, 0],
                             mode='lines',
                             line=dict(color='grey'),
                             showlegend=False))
    fig.add_trace(go.Scatter(x=[-0, 0],
                             y=[-20, 20],
                             mode='lines',
                             line=dict(color='grey'),
                             showlegend=False))
    fig.update_layout(template='plotly_white',
                      legend=dict(
                          orientation="h",
                          yanchor="bottom",
                          y=1.02,
                          xanchor="right",
                          x=1),
                      )
    fig.update_xaxes(range=(-2, 1),
                     # domain=(0, 0.8),
                     )
    fig.update_yaxes(range=(-2, 1),
                     # domain=(0, 0.8),
                     scaleanchor='x',
                     scaleratio=1)

    if show_annotation:
        spearman_r, spearman_p = stats.spearmanr(plot_df["log2_pA_over_total_rpm_3"],
                                                 plot_df["log2_pA_over_total_rpm_2"])
        print_text = f"<b>Correlation:</b>"
        print_text += f"<br>  Spearman R = {spearman_r:.4f}" \
                      f"<br>  Spearman p-val = {spearman_p:.2E}"
        print_text += f"<br><b>Cutoff:</b>" \
                      f"<br>  Reads/Gene ≥ {cutoff}"
        fig.add_annotation(text=print_text,
                           x=0.99, xref='paper', xanchor='right',
                           y=0.01, yref='paper', yanchor='bottom',
                           align='right',
                           bordercolor="darkgray",
                           borderwidth=2,
                           borderpad=6,
                           bgcolor="lightgray",
                           font=dict(family="Courier New, monospace",
                                     size=16),
                           showarrow=False)
    if not show_color_legend:
        fig.update(layout_coloraxis_showscale=False)
    fig.show()
    return fig


def _plot_multi_violin(filtered_df):
    fig = px.violin(filtered_df, x="gene_name", y="polya_length",
                    color="lib", points="all", box=True,
                    color_discrete_sequence=['#202020', '#606060', '#A0A0A0', '#C0C0C0'])
    fig.update_layout(margin={'l': 0, 'b': 40, 't': 10, 'r': 40},
                      yaxis_title=f"Distribution of PolyA Tail Length Calls",
                      legend=dict(orientation="h",
                                  yanchor="bottom",
                                  y=1.02,
                                  xanchor="left",
                                  x=0),
                      template='plotly_white',
                      violingap=0.1, violingroupgap=0)
    fig.update_traces(meanline_visible=True,
                      points='all',  # show all points
                      side='positive',
                      spanmode='hard',
                      pointpos=-0.32,  # could maybe go back to both sides and zero this...
                      marker=dict(opacity=0.5),
                      # jitter=0.05,  # add some jitter on points for better visibility
                      scalemode='width',  # scale violin plot area with total count
                      )
    fig.show()
    return fig


if __name__ == '__main__':
    comp_df, read_df = compare_libraries_w_pA_over_total(cutoff=80, force_compressed_df_build=False,
                                                         long_tail_cutoff=0.28,
                                                         short_tail_cutoff=0.57,
                                                         size_by_expression_quartile=False,
                                                         save_files=True,
                                                         locked_aspect=True,
                                                         color_by_bins=10)
    column_list = ["lib",
                   "gene_id",
                   "gene_name",
                   "polya_length",
                   "read_id"]
    gene_name_lists = [
        # First set of interesting genes based on ... ?
        # ['rpl-20', 'rps-15', 'nduo-4'],
        # ['osm-11', 'act-5', 'vha-3'],
        # ['mxl-3', 'aco-2', 'ucr-1'],
        # ['clik-1', 'F56D3.1', 'lin-42'],
        # ['gldc-1', 'metr-1', 'mlt-11', ],
        # ['vha-12', 'sca-1', 'vgln-1'],
        # Second set of genes based on the technique RPM fold change plot:
        # ['rpl-2', 'rpl-4', 'rpl-9'],
        # ['asp-6', 'trap-2', 'ifb-2'],
        # ['eef-1A.1', 'rpl-20', 'rps-5'],
        # ['ola-1', 'abu-13', 'fat-1'],
        # ['ahcy-1']
        # An extreme from the mean-median plots for set 3:
        # ['C23H5.8', 'T05E12.6', 'F46G10.1'],
        # Fairly on diagonal gene, but strange distributions:
        # ['icl-1', ]
        # For RNA Club 4/10
        # ['mlt-11'],
        # ['gldc-1'],
        # ['cuc-1'],
        # ['C23G10.2'],
        # ['C23H5.8'],
        # ['F15B9.8'],
        # ['asp-13'],
        # For RNA Club, from tail length scatter:
        # ['gldc-1'],
        # ['C23G10.2'],
        # ['icl-1'],
        # ['mlt-11'],
        # For Josh 04/13/22:
        ['rpl-1', 'rpl-10', 'rpl-20', 'rps-10', 'rps-20']
    ]
    for gene_name_set in gene_name_lists:
        violin_df = read_df[read_df["gene_name"].isin(gene_name_set)][column_list]
        _plot_multi_violin(violin_df)
