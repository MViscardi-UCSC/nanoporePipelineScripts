"""
dashForClickAndViewPolyATails.py
Marcus Viscardi,    November 08, 2021

Trying to use dash to be able to quickly take genes that
are off diagonal in the read length of polyA comparison
scatter plots, pass those to another plot to see their
distributions! preferably on click or hover!

Mostly using:
https://dash.plotly.com/interactive-graphing
https://dash.plotly.com/dash-core-components/graph
https://community.plotly.com/t/generate-another-plot-from-data-obtained-from-hovering-over-the-one-graph-and-a-pandas-dataframe/51848/4
"""
import seaborn

from nanoporePipelineCommon import pick_libs_return_paths_dict, get_dt, \
    find_newest_matching_file, load_and_merge_lib_parquets
import pandas as pd
import numpy as np

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)


def my_first_attempt():
    import dash
    from dash import dcc
    from dash import html
    from dash.dependencies import Input, Output
    import plotly.express as px

    hardcoded_df = pd.read_parquet("/data16/marcus/working/211101_nanoporeSoftLinks/"
                                   "210709_NanoporeRun_totalRNA_0639_L3/"
                                   "output_dir/merge_files/"
                                   "211103_mergedOnReads.plusStandards.parquet")
    df = hardcoded_df
    # Columns:
    # ['read_id', 'bit_flag', 'chr_id', 'chr_pos', 'mapq', 'cigar', 'sequence',
    #    'phred_qual', 'num_mismatches', 'transcript_strand',
    #    'type_of_alignment', 'strand', 'qc_tag_featc', 'qc_pass_featc',
    #    'gene_id', 'gene_name', 'leader_start', 'adapter_start', 'polya_start',
    #    'transcript_start', 'read_rate', 'polya_length', 'qc_tag_polya',
    #    'read_length', 'adapter', 'stds_q_start', 'stds_q_end', 'stds_strand',
    #    'stds_chr_id', 'stds_chr_len', 'stds_r_start', 'stds_r_end',
    #    'stds_mlen', 'stds_blen', 'stds_mapq', 'stds_type_of_alignment',
    #    'stds_ts', 'stds_cigar']

    app = dash.Dash(__name__)

    app.layout = html.Div([
        html.Div([
            dcc.Dropdown(
                id='xaxis-column',
                options=[{'label': i, 'value': i} for i in df.columns],
                value='read_length',
            ),
            dcc.Dropdown(
                id='yaxis-column',
                options=[{'label': i, 'value': i} for i in df.columns],
                value='polya_length',
            )]),
        # Left block with one plot
        html.Div([
            dcc.Graph(
                id='primary-scatter',
                hoverData={'points': [{'customdata': 'WBGene00000001'}]}
            )
        ], style={'width': '49%',
                  'display': 'inline-block',
                  'padding': '0 20'}),

        # Right block with two plots
        html.Div([
            dcc.Graph(id='x-histogram'),
            dcc.Graph(id='y-histogram'),
        ], style={'width': '49%',
                  'display': 'inline-block'}),
    ])

    @app.callback(
        Output('primary-scatter', 'figure'),
        [Input('xaxis-column', 'value'),
         Input('yaxis-column', 'value')])
    def main_plot(xaxis_column_name, yaxis_columns_name):
        fig = px.scatter(df, x=xaxis_column_name, y=yaxis_columns_name,
                         custom_data=["gene_id", "gene_name"],
                         hover_name="gene_name")
        fig.update_layout(margin={'l': 40, 'b': 40, 't': 10, 'r': 0}, hovermode='closest')
        fig.update_layout(clickmode="event+select")
        fig.update_traces(marker_size=2)
        return fig

    def plot_histogram(dff, title, plot_col):
        fig = px.histogram(dff, x=plot_col)
        fig.update_xaxes(showgrid=False)
        fig.add_annotation(x=0, y=0.85, xanchor='left', yanchor='bottom',
                           xref='paper', yref='paper', showarrow=False, align='left',
                           text=title)
        fig.update_layout(height=225, margin={'l': 20, 'b': 30, 'r': 10, 't': 10})
        return fig

    @app.callback(
        Output("x-histogram", "figure"),
        [Input('primary-scatter', 'selectedData'),
         Input('xaxis-column', 'value')])
    def update_x_histo(hoverData, xaxis_column_name):
        gene_id = hoverData['points'][0]['customdata'][0]
        gene_name = hoverData['points'][0]['customdata'][1]
        dff = df[df["gene_id"] == gene_id]
        return plot_histogram(dff, gene_name, xaxis_column_name)

    @app.callback(
        Output("y-histogram", "figure"),
        [Input('primary-scatter', 'selectedData'),
         Input('yaxis-column', 'value')])
    def update_y_histo(hoverData, yaxis_column_name):
        gene_id = hoverData['points'][0]['customdata'][0]
        gene_name = hoverData['points'][0]['customdata'][1]
        dff = df[df["gene_id"] == gene_id]
        return plot_histogram(dff, gene_name, yaxis_column_name)

    if __name__ == '__main__':
        app.run_server(debug=True)


# def load_and_merge_lib_parquets(lib_list, genomeDir=f"/data16/marcus/genomes/elegansRelease100/",
#                                 drop_unassigned=True, drop_failed_polya=True,
#                                 drop_sub_n=5, keep_transcript_info=False,
#                                 subsample_each_lib=False, read_pos_in_groupby=False,
#                                 add_nucleotide_fractions=False,
#                                 group_by_t5=False) -> [pd.DataFrame, pd.DataFrame]:
#     read_assignment_path = find_newest_matching_file(f"{genomeDir}/*.allChrs.parquet")
#     read_assignment_df = pd.read_parquet(read_assignment_path)
#     # Loop through each library name in the list and for each:
#     #   1. Load the Parquet
#     #   2. Merge this w/ Josh's assign reads based on chr_pos
#     #   3. Create a column to retain the library identity
#     #   3. Concatenate these dataframe into one large dataframe
#     #       NOTE: This structure can be seperated again based on
#     #       the "lib" column added in the previous step
#     path_dict = pick_libs_return_paths_dict(lib_list, file_suffix="parquet")
#     df_dict = {}
#     for library_name, parquet_path in path_dict.items():
#         print(f"Loading parquet for {library_name} lib. . .")
#         lib_df = pd.read_parquet(parquet_path)
#         # if isinstance(subsample_each_lib, int):
#         #     lib_df - lib_df.sample(subsample_each_lib)
#         # elif subsample_each_lib:
#         #     lib_df = lib_df.sample(10000)
#         df_dict[library_name] = lib_df
# 
#     # This is a cute way to quickly merge all of these dfs into one, while retaining lib info.
#     #   B/c I am still a little scared of MultiIndexed dataframes, I used the reset_index steps
#     #   to push the mutliindex back into columns. Maybe someday I'll use the multiindex!
#     multi_df = pd.concat(df_dict.values(), keys=df_dict.keys())
#     multi_df.index.set_names(("lib", "old_index"), inplace=True)
#     super_df = multi_df.reset_index(level="lib").reset_index(drop=True)
#     print(f"Starting assignment merge . . .", end="")
#     
#     # Some major issues with ugly column names coming through, plan to clean them up:
#     read_assignment_cols = read_assignment_df.columns.to_list()
#     read_assignment_cols.remove('chr_id')
#     read_assignment_cols.remove('chr_pos')
#     read_assignment_cols_to_drop = [col for col in read_assignment_cols if col in super_df.columns]
#     super_df.drop(columns=read_assignment_cols_to_drop, inplace=True)
#     
#     super_df = super_df.merge(read_assignment_df, on=["chr_id", "chr_pos"],
#                               how="left", suffixes=["_originalOutput",
#                                                     ""])
#     keep_columns = ["lib",
#                     "read_id",
#                     "chr_id",
#                     "chr_pos",
#                     "gene_id",
#                     "gene_name",
#                     "cigar",
#                     "sequence",
#                     "polya_length",
#                     "strand",
#                     ]
#     if keep_transcript_info:
#         for col in ["transcript_id", "to_start", "to_stop"]:
#             keep_columns.append(col)
#     if group_by_t5:
#         keep_columns.append('t5')
#     super_df = super_df[keep_columns].drop_duplicates()
#     if drop_failed_polya:
#         print(f"\nRead counts pre-failed-polyA call drop:   {super_df.shape[0]}")
#         super_df = super_df[~super_df["polya_length"].isna()]
#         print(f"Read counts post-failed-polyA call drop:  {super_df.shape[0]}")
#     print(f"\rFinished assignment merge!")
#     
#     # Post-assignment cleanups:
#     if drop_unassigned:
#         print(f"Read counts post gene assignment:  {super_df.shape[0]}")
#         super_df = super_df[~super_df["gene_id"].isna()].reset_index(drop=True)
#         print(f"Read counts post unassigned drop:  {super_df.shape[0]}")
#         # 220201: I don't know why we were doing below...?
#         # concat_df = concat_df[concat_df["strand_originalOutput"] == concat_df["strand_forPlot"]].reset_index(drop=True)
#         # print(f"Read counts post consistent-assignment check: {concat_df.shape[0]}")
#         if keep_transcript_info:
#             super_df = super_df.astype({"to_stop": "int64",
#                                         "to_start": "int64"})
#     super_df["read_length"] = super_df["sequence"].apply(len)
# 
#     # Create the groupby dataframe:
#     print(f"Creating groupby dataframe merged on lib, chr_id, gene_id, and gene_name")
#     groupby_col_list = ["lib", "chr_id", "gene_id", "gene_name"]
#     if keep_transcript_info:
#         print(f"\t+transcript_id")
#         groupby_col_list.append("transcript_id")
#     if group_by_t5:
#         print(f"\t+t5 tag")
#         groupby_col_list.append("t5")
#     
#     # Holy crap, the observed=True helps to keep this from propagating out to 129,151,669,691,968 rows...
#     groupby_obj = super_df.groupby(groupby_col_list, observed=True)
#     if not keep_transcript_info:
#         compressed_prefix = "gene"
#     else:
#         compressed_prefix = "transcript"
#     compressed_df = groupby_obj["read_id"].apply(len).to_frame(name=f"{compressed_prefix}_hits")
#     compressed_df["mean_polya_length"] = groupby_obj["polya_length"].mean()
#     compressed_df["mean_read_length"] = groupby_obj["read_length"].mean()
#     if read_pos_in_groupby:
#         compressed_df['stop_distances'] = groupby_obj["to_stop"].apply(list).to_frame(name="stop_distances")
#         compressed_df['start_distances'] = groupby_obj["to_start"].apply(list).to_frame(name="stop_distances")
#     
#     # RPM and fractional hits calculations
#     # Need to first create columns of NA values, tobe overwritten
#     compressed_df[f"{compressed_prefix}_rpm"] = pd.NA
#     compressed_df[f"{compressed_prefix}_frac_hits"] = pd.NA
#     # Only look at one library at a time (so the normalization is per lib not whole df)
#     for lib in compressed_df.index.unique(level='lib').to_list():
#         # Create the 'norm_factor' which will be the total # of read hits in that lib
#         norm_factor = compressed_df.query(f"lib == '{lib}'")[f"{compressed_prefix}_hits"].sum()
#         # Turn the total number of read hits into the 'million of read hits'
#         rpm_norm_factor = norm_factor / 1000000
#         # For each library divide gene_hits by the rpm norm factor to get rpm
#         gene_rpm_series = compressed_df.query(f"lib == '{lib}'")[f"{compressed_prefix}_hits"] / rpm_norm_factor
#         # Use a series fill, so that we can fill that library's part of the DF without effecting others
#         compressed_df[f"{compressed_prefix}_rpm"] = compressed_df[f"{compressed_prefix}_rpm"].\
#             fillna(value=gene_rpm_series)
#         # Same as above, but with fraction of hits, rather than a rpm calc (practically same thing)
#         gene_frac_hits_series = compressed_df.query(f"lib == '{lib}'")[f"{compressed_prefix}_hits"] / norm_factor
#         compressed_df[f"{compressed_prefix}_frac_hits"] = compressed_df[f"{compressed_prefix}_frac_hits"].\
#             fillna(value=gene_frac_hits_series)
#     
#     # Requirement for min number of gene/transcript hits
#     if isinstance(drop_sub_n, int):
#         print(f"Gene counts pre sub-{drop_sub_n} {compressed_prefix}_hits drop:  {compressed_df.shape[0]}")
#         compressed_df = compressed_df[compressed_df[f"{compressed_prefix}_hits"] >= drop_sub_n]
#         print(f"Gene counts post sub-{drop_sub_n} {compressed_prefix}_hits drop:  {compressed_df.shape[0]}")
#     # Reset index at the end,
#     #   we didn't retain any info w/ the index, so it doesn't help much
#     compressed_df = compressed_df.reset_index()
#     
#     if add_nucleotide_fractions:
#         print(f"Adding nucleotide content information to genes!")
#         path_to_gc = "/data16/marcus/genomes/elegansRelease100/" \
#                      "Caenorhabditis_elegans.WBcel235.cdna.all.fa.GCcontent.parquet"
#         gc_df = pd.read_parquet(path_to_gc).drop(columns=["chr_id"])
#         compressed_df = compressed_df.merge(gc_df, on=["gene_id", "gene_name"], how="left")
#     
#     return super_df, compressed_df


def distributions_of_polya_tails(libs, force_compressed_df_build=False):
    import dash
    import json
    from dash import dcc, html, callback_context
    import dash_bootstrap_components as dbc
    from dash.dependencies import Input, Output
    import dash_daq as daq
    import plotly.express as px
    import plotly.graph_objects as go
    import plotly.io as pio

    if len(libs) < 2:
        raise ValueError(f"Please provide 2 or more libraries, only {len(libs)} given.")

    if not force_compressed_df_build:
        compressed_save_search = f"./testInputs/*_{'_'.join(libs)}.compressed.parquet"
        reads_save_search = f"./testInputs/*_{'_'.join(libs)}.reads.parquet"
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
        reads_df, compressed_df = load_and_merge_lib_parquets(libs, drop_sub_n=1)
        compressed_save_path = f"./testInputs/{get_dt(for_file=True)}_{'_'.join(libs)}.compressed.parquet"
        reads_save_path = f"./testInputs/{get_dt(for_file=True)}_{'_'.join(libs)}.reads.parquet"
        compressed_df.to_parquet(compressed_save_path)
        reads_df.to_parquet(reads_save_path)
        print(f"Saved new compressed file to: {compressed_save_path}")
        print(f"Saved new reads file to:      {reads_save_path}")
    lib_list = compressed_df.lib.unique().tolist()

    external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
    styles = {
        'pre': {
            'border': 'thin lightgrey solid',
            'overflowX': 'scroll',
            'overflowY': 'scroll',
        }
    }

    app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

    app.layout = html.Div(
        [
            # Option to change the libraries at the top of the plot
            html.Div(
                [
                    dcc.Markdown(
                        """**Select X and Y libraries**"""
                    ),
                    dcc.Dropdown(
                        id='xaxis-lib',
                        options=[{'label': i, 'value': i} for i in lib_list],
                        value=lib_list[0]
                    ),
                    dcc.Dropdown(
                        id='yaxis-lib',
                        options=[{'label': i, 'value': i} for i in lib_list],
                        value=lib_list[1]
                    ),
                    html.Div(
                        [
                            dcc.Markdown(
                                """**Select minimum reads/gene to allow**"""
                            ),
                            dcc.Slider(
                                id='min-hits-slider',
                                min=5, max=100,
                                value=40,
                                marks={str(n): str(n) for n in range(5, 105, 5)},
                                step=None
                            ),
                            daq.BooleanSwitch(
                                id='trendline-switch', on=False,
                                label="Trendline for scatter:"
                            ),
                        ]
                    ),
                ]
            ),
            # Plots row
            html.Div(
                [
                    html.Div(
                        [
                            # First plot with title and 8/12 width space
                            html.H3(
                                'Rocket Plot of Mean Tail Lengths'
                            ),
                            dcc.Graph(
                                id='primary-scatter',
                                hoverData={'points': [{'customdata': ['WBGene00010964', 'ctc-1', 'lots', 'lots']}]}
                            )
                        ],
                        className="six columns", style={'display': 'inline-block'}
                    ),

                    html.Div(
                        [
                            # Second plot with less for the width space
                            html.H3(
                                'Violin Plot of Selected Data'
                            ),
                            dcc.Graph(
                                id='violin-plot'
                            )
                        ],
                        className="six columns", style={'display': 'inline-block'}
                    ),
                ],
                className='row'
            ),
            # Info and Button Row
            html.Div(
                className='row',
                children=[
                    html.Div([
                        dcc.Markdown("""**Hover Data**
                    Mouse over values in the graph."""),
                        html.Pre(id='hover-data', style=styles['pre'])
                    ], className='four columns'),
                    html.Div([
                        dcc.Markdown("""
                **Selection Data**

                Click points to select data.
                Or to select multiple:
                Shift+click, or choose the lasso/rectangle select tool in the graph's menu
                bar and then select points in the graph.
            """),
                        html.Pre(id='selected-data', style=styles['pre']),
                    ], className='four columns'),
                    html.Div(
                        [
                            html.Div([
                                html.Button('Popout Scatter', id='btn-nclicks-1', n_clicks=0),
                                html.Button('Popout Violin', id='btn-nclicks-2', n_clicks=0)]),
                            html.Div([
                                html.Button('Save Scatter', id='btn-nclicks-3', n_clicks=0),
                                html.Button('Save Violin', id='btn-nclicks-4', n_clicks=0)]),
                            html.Div(id='container-button-timestamp'),
                        ]
                    )
                ]
            )
        ]
    )

    @app.callback(
        Output('primary-scatter', 'figure'),
        [Input('xaxis-lib', 'value'),
         Input('yaxis-lib', 'value'),
         Input('min-hits-slider', 'value'),
         Input('selected-data', 'children'),
         Input('trendline-switch', 'on')])
    def main_plot(xaxis_library, yaxis_library, min_hits,
                  selectedData, trendline_switch) -> go.Figure:
        if selectedData != "No points selected":
            selected_gene_ids = [hit["Gene ID"] for hit in json.loads(selectedData)]
            selected_gene_names = [hit["Gene Name"] for hit in json.loads(selectedData)]
        else:
            selected_gene_ids, selected_gene_names = [], []
        min_hit_df = compressed_df[compressed_df['gene_hits'] >= min_hits]
        max_mean_tail = min_hit_df.mean_polya_length.max()
        max_mean_tail += 10
        min_mean_tail = min_hit_df.mean_polya_length.min()
        min_mean_tail -= 10
        x_axis_df = min_hit_df[min_hit_df.lib == xaxis_library][["gene_id",
                                                                 "gene_name",
                                                                 "gene_hits",
                                                                 "mean_polya_length"]]
        y_axis_df = min_hit_df[min_hit_df.lib == yaxis_library][["gene_id",
                                                                 "gene_name",
                                                                 "gene_hits",
                                                                 "mean_polya_length"]]

        plot_df = pd.merge(x_axis_df, y_axis_df, on=["gene_id", "gene_name"],
                           suffixes=(f"_{xaxis_library}",
                                     f"_{yaxis_library}"))

        if trendline_switch:
            trend = "ols"
        else:
            trend = None
        fig = px.scatter(plot_df, x=f"mean_polya_length_{xaxis_library}", y=f"mean_polya_length_{yaxis_library}",
                         custom_data=["gene_id", "gene_name"],
                         hover_name="gene_name", hover_data=["gene_id",
                                                             f"gene_hits_{xaxis_library}",
                                                             f"gene_hits_{yaxis_library}"],
                         trendline=trend)
        fig.update_layout(margin={'l': 40, 'b': 40, 't': 10, 'r': 0},
                          hovermode='closest', clickmode="event+select",
                          xaxis_title=f"Mean polyA Tail Length (Lib: {xaxis_library})",
                          yaxis_title=f"Mean polyA Tail Length (Lib: {yaxis_library})",
                          template='plotly_white')
        fig.update_traces(marker=dict(size=7, color='darkgray', opacity=0.7))
        fig.add_trace(go.Scatter(x=[min_mean_tail, max_mean_tail],
                                 y=[min_mean_tail, max_mean_tail],
                                 mode='lines',
                                 line=dict(color='black',
                                           dash='dash'),
                                 showlegend=False))
        if selectedData != "No points selected":
            selected_df = plot_df[plot_df["gene_id"].isin(selected_gene_ids)]
            fig.add_trace(go.Scatter(mode='markers',
                                     x=selected_df[f"mean_polya_length_{xaxis_library}"],
                                     y=selected_df[f"mean_polya_length_{yaxis_library}"],
                                     customdata=list(zip(selected_df[f"gene_id"],
                                                         selected_df[f"gene_name"],
                                                         selected_df[f"gene_hits_{xaxis_library}"],
                                                         selected_df[f"gene_hits_{yaxis_library}"])),
                                     marker=dict(color='red',
                                                 size=7,
                                                 line=dict(color='black',
                                                           width=2)),
                                     name='Selected Genes',
                                     showlegend=False))
            fig.update_layout(legend=dict(orientation="h",
                                          yanchor="bottom",
                                          y=1.02,
                                          xanchor="left",
                                          x=0))
            # selected_data_points = px.scatter(plot_df[plot_df["gene_id"].isin(selected_gene_ids)],
            #                                   x=f"mean_polya_length_{xaxis_library}",
            #                                   y=f"mean_polya_length_{yaxis_library}",
            #                                   custom_data=["gene_id", "gene_name"],
            #                                   hover_name="gene_name", hover_data=["gene_id",
            #                                                                       f"gene_hits_{xaxis_library}",
            #                                                                       f"gene_hits_{yaxis_library}"],
            #                                   color=1)
            # fig.add_trace(selected_data_points.data[0])
        fig.update_xaxes(range=[min_mean_tail, max_mean_tail])
        fig.update_yaxes(range=[min_mean_tail, max_mean_tail])
        return fig

    def _plot_split_violin(filtered_df, x_lib, y_lib):
        fig = go.Figure()

        fig.add_trace(go.Violin(x=filtered_df['gene_name'][filtered_df['lib'] == x_lib],
                                y=filtered_df['polya_length'][filtered_df['lib'] == x_lib],
                                legendgroup=x_lib, scalegroup=x_lib, name=x_lib,
                                side='negative',
                                line_color='MediumPurple',
                                )
                      )
        fig.add_trace(go.Violin(x=filtered_df['gene_name'][filtered_df['lib'] == y_lib],
                                y=filtered_df['polya_length'][filtered_df['lib'] == y_lib],
                                legendgroup=y_lib, scalegroup=y_lib, name=y_lib,
                                side='positive',
                                line_color='LightSeaGreen',
                                )
                      )
        fig.update_traces(meanline_visible=True,
                          points='all',
                          jitter=0.05,
                          scalemode='count')
        fig.update_layout(violingap=0, violinmode='overlay',
                          margin={'l': 0, 'b': 40, 't': 10, 'r': 40},
                          yaxis_title=f"Distribution of PolyA Tail Length Calls",
                          legend=dict(orientation="h",
                                      yanchor="bottom",
                                      y=1.02,
                                      xanchor="left",
                                      x=0),
                          template='plotly_white')
        return fig

    def _plot_multi_violin(filtered_df):
        fig = px.violin(filtered_df, x="gene_name", y="polya_length",
                        color="lib", points="all")
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
                          pointpos=-0.1,  # could maybe go back to both sides and zero this...
                          marker=dict(opacity=0.5),
                          # jitter=0.05,  # add some jitter on points for better visibility
                          # scalemode='count',  # scale violin plot area with total count
                          )
        seaborn.set_palette('colorblind')
        return fig

    @app.callback(
        Output('violin-plot', 'figure'),
        [Input('primary-scatter', 'selectedData'),
         Input('xaxis-lib', 'value'),
         Input('yaxis-lib', 'value')])
    def handle_violin(selectedData, x_lib, y_lib):
        column_list = ["lib",
                       "gene_id",
                       "gene_name",
                       "polya_length",
                       "read_id"]
        gene_id_list = []
        if selectedData:
            for point_dict in selectedData['points']:
                gene_id_list.append(point_dict['customdata'][0])
            # gene_id, gene_name = selectedData['points'][0]['customdata'][:2]
            violin_df = reads_df[reads_df["gene_id"].isin(gene_id_list)][column_list]
        else:
            violin_df = pd.DataFrame(dict(zip(column_list, [[] for _ in range(len(column_list))])))
        # return _plot_split_violin(violin_df, x_lib, y_lib)
        return _plot_multi_violin(violin_df)

    @app.callback(
        Output('hover-data', 'children'),
        [Input('primary-scatter', 'hoverData'),
         Input('xaxis-lib', 'value'),
         Input('yaxis-lib', 'value')])
    def display_hover_data(hoverData, x_lib, y_lib):
        return_data = []
        if not hoverData:
            return None
        for point_dict in hoverData['points']:
            gene_id, gene_name, x_counts, y_counts = point_dict['customdata']
            return_data.append({'Gene Name': gene_name,
                                'Gene ID': gene_id,
                                f'{x_lib} hits': x_counts,
                                f'{y_lib} hits': y_counts})
        return json.dumps(return_data, indent=2)

    @app.callback(
        Output('selected-data', 'children'),
        [Input('primary-scatter', 'selectedData'),
         Input('xaxis-lib', 'value'),
         Input('yaxis-lib', 'value')])
    def display_selected_data(selectedData, x_lib, y_lib):
        if not selectedData:
            return "No points selected"
        return_data = []
        for point_dict in selectedData['points']:
            gene_id, gene_name, x_counts, y_counts = point_dict['customdata']
            return_data.append({'Gene Name': gene_name,
                                'Gene ID': gene_id,
                                f'{x_lib} hits': x_counts,
                                f'{y_lib} hits': y_counts})
        return json.dumps(return_data, indent=2)

    @app.callback(
        Output('container-button-timestamp', 'children'),
        [Input('btn-nclicks-1', 'n_clicks'),
         Input('btn-nclicks-2', 'n_clicks'),
         Input('btn-nclicks-3', 'n_clicks'),
         Input('btn-nclicks-4', 'n_clicks'),
         Input('primary-scatter', 'figure'),
         Input('violin-plot', 'figure')])
    def save_button_click(pop_scatter, pop_violin, save_scatter, save_violin, scatter_fig, violin_fig):
        changed_id = [p['prop_id'] for p in callback_context.triggered][0]
        if 'btn-nclicks-1' in changed_id:
            msg = 'Popout Scatter'
            # print(scatter_fig)
            pio.show(scatter_fig)
        elif 'btn-nclicks-2' in changed_id:
            msg = 'Popout Violin'
            # print(violin_fig)
            pio.show(violin_fig)
        elif 'btn-nclicks-3' in changed_id:
            msg = 'Save Scatter'
            save_file_loc = f"./testOutputs/{get_dt(for_file=True)}_dashScatter.svg"
            go.Figure(scatter_fig).write_image(save_file_loc,
                                               width=600, height=600, scale=2)
        elif 'btn-nclicks-4' in changed_id:
            msg = 'Save Violin'
            save_file_loc = f"./testOutputs/{get_dt(for_file=True)}_dashViolin.svg"
            go.Figure(violin_fig).write_image(save_file_loc,
                                              width=1000, height=500, scale=2)
        else:
            msg = '*None of the buttons*'
        return html.Div(f"{msg} was most recently clicked")

    app.run_server(debug=False, dev_tools_hot_reload=False, host='127.0.0.1', port=16000)


if __name__ == '__main__':
    from sys import argv

    libraries_to_run = argv[1:]
    print(f"Running w/ libraries: {libraries_to_run}")
    distributions_of_polya_tails(libraries_to_run,
                                 force_compressed_df_build=False,
                                 )
