"""
dashForClickAndViewAll.py
Marcus Viscardi,    January 11, 2022

So the goal is to mostly recapitulate what I was doing with the original
    dashForClickAndViewPolyATails.py without the violins, but with the
    addition of more options for what to plot.

One of the main reasons to do this is that I'll be able to "quickly" produce
    several of the figures that I want for the polyA-vs-total paper

A lot of this is going to be getting ripped from other scripts!

From dashForClickAndViewPolyATails.py:
    Mostly using:
    https://dash.plotly.com/interactive-graphing
    https://dash.plotly.com/dash-core-components/graph
    https://community.plotly.com/t/generate-another-plot-from-data-obtained-from-hovering-over-the-one-graph-and-a-pandas-dataframe/51848/4
"""
import os
from pprint import pprint

import numpy as np

from nanoporePipelineCommon import get_dt, find_newest_matching_file, load_and_merge_lib_parquets
import pandas as pd

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)


def main(libs, force_compressed_df_build=False, abs_min_cutoff=5):
    global COMPRESSED_DF, PLOT_DF, SELECTED_DF

    import dash
    from dash import dcc, html, callback_context, dash_table
    from dash.dependencies import Input, Output

    import plotly.express as px
    import plotly.graph_objects as go
    import plotly.io as pio

    from scipy import stats
    from math import log10

    if len(libs) < 2:
        raise ValueError(f"Please provide 2 or more libraries, only {len(libs)} given.")

    if not force_compressed_df_build:
        search_path = f"./testInputs/*_{'_'.join(libs)}.compressed.parquet"
        try:
            compressed_parquet_path = find_newest_matching_file(search_path)
            COMPRESSED_DF = pd.read_parquet(compressed_parquet_path)
        except ValueError:
            print(f"Couldn't find pre-processed file at: {search_path}\nGoing to load from library files!")
            force_compressed_df_build = True
    if force_compressed_df_build:
        reads_df, COMPRESSED_DF = load_and_merge_lib_parquets(libs,
                                                              drop_sub_n=abs_min_cutoff,
                                                              add_nucleotide_fractions=True,
                                                              add_tail_groupings=True)
        save_path = f"./testInputs/{get_dt(for_file=True)}_{'_'.join(libs)}.compressed.parquet"
        COMPRESSED_DF.to_parquet(save_path)
        print(f"Saved new compressed file to: {save_path}")
    lib_list = COMPRESSED_DF.lib.unique().tolist()
    COMPRESSED_DF['is_MtDNA'] = COMPRESSED_DF['chr_id'] == 'MtDNA'
    COMPRESSED_DF['median_mean_diff_tails'] = COMPRESSED_DF['median_polya_length'] - COMPRESSED_DF['mean_polya_length']
    plottable_columns = COMPRESSED_DF.select_dtypes(include='number').columns.to_list()
    color_by_columns = ["None", "is_MtDNA", "chr_id",
                        "t_fraction", "a_fraction",
                        "c_fraction", "g_fraction",
                        "tail_groupings_group",
                        "frac_sub_50_tails",
                        ]

    app = dash.Dash(__name__)

    app.layout = html.Div(
        [
            html.Div(  # Options Row
                [
                    html.Div(  # Left section of options row
                        [
                            html.Label('X Library:'),
                            dcc.RadioItems(
                                id='xaxis-lib',
                                options=[{'label': i, 'value': i} for i in lib_list],
                                value=lib_list[0],
                                persistence=True, persistence_type='local',
                                labelStyle={'display': 'inline-block', 'text-align': 'justify'},
                            ),

                            html.Br(),
                            html.Label('Y Library:'),
                            dcc.RadioItems(
                                id='yaxis-lib',
                                options=[{'label': i, 'value': i} for i in lib_list],
                                value=lib_list[1],
                                persistence=True, persistence_type='local',
                                labelStyle={'display': 'inline-block', 'text-align': 'justify'},
                            ),
                        ], style={'padding': 10, 'flex': 3}
                    ),
                    html.Div(  # Middle section of options row
                        [
                            html.Label('Options'),
                            dcc.Checklist(
                                id='options',
                                options=[
                                    {'label': 'Log X', 'value': 'log_x'},
                                    {'label': 'Log Y', 'value': 'log_y'},
                                    {'label': 'Trend Line', 'value': 'trendline'},
                                    {'label': 'Diagonal Line', 'value': 'diagonal'},
                                    {'label': 'Shape by MtDNA', 'value': 'MtDNA_shape'},
                                    # {'label': 'Color By Chr', 'value': 'color_chr'},
                                ],
                                value=['log_x', 'log_y'],
                                persistence=True, persistence_type='local',
                                labelStyle={'display': 'inline-block', 'text-align': 'justify'},
                            ),

                            html.Br(),
                            html.Label('Select by gene_name:'),
                            dcc.Input(
                                id='gene_name_select',
                                value='',
                                type='text',
                                debounce=True,
                            ),

                            html.Br(),
                            html.Label('Minimum number of hits/gene to allow:'),
                            dcc.Slider(
                                id='min-hits-slider',
                                min=5 * (abs_min_cutoff // 5), max=100,
                                value=40, persistence=True, persistence_type='local',
                                marks={str(n): str(n) for n in
                                       range(5 * (abs_min_cutoff // 5), 100 + abs_min_cutoff, 5)},
                                step=None
                            ),
                            html.Br(),
                            html.Label('Shorter tail group cutoff:', id="short-group-cutoff-label"),
                            dcc.Slider(
                                id='short-group-cutoff-slider',
                                min=0, max=1, step=0.01,
                                value=0.45,
                                persistence=True, persistence_type='local',
                                tooltip={"placement": "bottom", "always_visible": True},
                            ),
                            html.Br(),
                            html.Label('Longer tail group cutoff:', id="long-group-cutoff-label"),
                            dcc.Slider(
                                id='long-group-cutoff-slider',
                                min=0, max=1, step=0.01,
                                value=0.1,
                                persistence=True, persistence_type='local',
                                tooltip={"placement": "bottom", "always_visible": True},
                            ),
                        ], style={'padding': 10, 'flex': 6}
                    ),
                    html.Div(  # Right section of options row
                        [
                            html.Label('X Plot:'),
                            dcc.RadioItems(
                                id='xaxis-col',
                                options=[{'label': i, 'value': i} for i in plottable_columns],
                                value=plottable_columns[0],
                                persistence=True, persistence_type='local',
                                labelStyle={'display': 'inline-block', 'text-align': 'justify'},
                            ),

                            html.Br(),
                            html.Label('Y Plot:'),
                            dcc.RadioItems(
                                id='yaxis-col',
                                options=[{'label': i, 'value': i} for i in plottable_columns],
                                value=plottable_columns[0],
                                persistence=True, persistence_type='local',
                                labelStyle={'display': 'inline-block', 'text-align': 'justify'},
                            ),

                            html.Br(),
                            html.Label('Color By:'),
                            dcc.RadioItems(
                                id='color-col',
                                options=[{'label': i, 'value': i} for i in color_by_columns],
                                value=color_by_columns,
                                persistence=True, persistence_type='local',
                                labelStyle={'display': 'inline-block', 'text-align': 'justify'},
                            ),
                        ], style={'padding': 10, 'flex': 3}
                    ),
                ], style={'display': 'flex', 'flex-direction': 'row'},
            ),
            html.Div(  # Plotting Row
                [
                    html.Div(  # Left section of plotting row
                        [
                            dcc.Graph(
                                id='primary-scatter',
                            )
                        ], style={'padding': 10, 'flex': '0 0 600px'}
                    ),
                    html.Div(  # Right section of the plotting row, with buttons
                        [
                            html.Br(), html.Br(),
                            html.Div(
                                [
                                    html.Button('Popout Scatter', id='btn-nclicks-1', n_clicks=0),
                                    html.Button('Quick Save', id='btn-nclicks-2', n_clicks=0)
                                ]
                            ),
                            html.Div(id='container-button-timestamp'),
                            html.Br(), html.Br(),
                            html.Div(
                                [
                                    html.Button('Save Scatter As SVG', id='btn-nclicks-3', n_clicks=0),
                                    html.Button('Save Scatter As PNG', id='btn-nclicks-4', n_clicks=0),
                                    dcc.Download(id='download_image'),
                                ]
                            ),
                            html.Br(), html.Br(),
                            html.Div(
                                [
                                    html.Button('Reset Selected Data', id='btn-nclicks-5', n_clicks=0),
                                ]
                            )
                        ], style={'padding': 10, 'flex': 1}
                    )
                ], style={'display': 'flex', 'flex-direction': 'row'},
            ),
            html.Div(  # Selected data table
                id='selected-datatable'
            )
        ]
    )

    @app.callback(
        Output('primary-scatter', 'figure'),
        [Input('xaxis-lib', 'value'),
         Input('yaxis-lib', 'value'),
         Input('xaxis-col', 'value'),
         Input('yaxis-col', 'value'),
         Input('min-hits-slider', 'value'),
         Input('selected-datatable', 'children'),
         Input('options', 'value'),
         Input('color-col', 'value'),
         Input('short-group-cutoff-slider', 'value'),
         Input('long-group-cutoff-slider', 'value'),
         ])
    def main_plot(xaxis_library, yaxis_library,
                  xaxis_col, yaxis_col,
                  min_hits,
                  selected_data, options_list,
                  color_by_col, short_group_cutoff, long_group_cutoff) -> go.Figure:
        global COMPRESSED_DF, PLOT_DF, SELECTED_DF
        extended_xaxis_col = f"{xaxis_col}_{xaxis_library}"
        extended_yaxis_col = f"{yaxis_col}_{yaxis_library}"

        if color_by_col in ["tail_groupings_group", "frac_sub_50_tails"]:
            COMPRESSED_DF['tail_groupings_group'] = pd.cut(COMPRESSED_DF[f'frac_sub_50_tails'],
                                                           bins=[0.0, long_group_cutoff, short_group_cutoff, 1.0],
                                                           labels=['long_tailed',
                                                                   'ungrouped',
                                                                   'short_tailed'],
                                                           include_lowest=True)

        # Below tries to make some "smart" bounds around the plot!
        min_hit_df = COMPRESSED_DF[COMPRESSED_DF['gene_hits'] >= min_hits]

        max_plot_x = min_hit_df[min_hit_df['lib'] == xaxis_library][xaxis_col].max()
        min_plot_x = min_hit_df[min_hit_df['lib'] == xaxis_library][xaxis_col].min()

        max_plot_y = min_hit_df[min_hit_df['lib'] == yaxis_library][yaxis_col].max()
        min_plot_y = min_hit_df[min_hit_df['lib'] == yaxis_library][yaxis_col].min()

        # TODO: Something to handle if x and y libs are the same (just make one df and copy it?)
        x_axis_df = min_hit_df[min_hit_df.lib == xaxis_library].drop(columns='lib')
        y_axis_df = min_hit_df[min_hit_df.lib == yaxis_library].drop(columns='lib')
        if xaxis_library == yaxis_library:
            yaxis_library = 'y'
        merge_suffixes = (f"_{xaxis_library}",
                          f"_{yaxis_library}")
        PLOT_DF = pd.merge(x_axis_df, y_axis_df, on=["gene_id", "gene_name", "chr_id", "is_MtDNA",
                                                     "t_fraction", "a_fraction",
                                                     "c_fraction", "g_fraction"],
                           suffixes=merge_suffixes)

        if 'trendline' in options_list:
            trend = "ols"
        else:
            trend = None
        if 'MtDNA_shape' in options_list:
            shape_col = 'is_MtDNA'
        else:
            shape_col = None
        if color_by_col != "None":
            if color_by_col in PLOT_DF.columns.to_list():
                color_by = color_by_col
            elif color_by_col in ["tail_groupings_group", "frac_sub_50_tails"]:
                color_by = f"{color_by_col}_{yaxis_library}"
        else:
            color_by = None

        color_discrete_map = {'long_tailed': 'cornflowerblue',
                              'short_tailed': 'coral',
                              'ungrouped': 'slategrey',
                              False: 'black',
                              True: 'coral'}
        
        fig = px.scatter(PLOT_DF,
                         x=extended_xaxis_col,
                         y=extended_yaxis_col,
                         custom_data=["gene_id", "gene_name"],
                         symbol=shape_col,
                         symbol_sequence=[0, 3],
                         hover_name="gene_name", hover_data=["gene_id",
                                                             f"gene_hits_{xaxis_library}",
                                                             f"gene_hits_{yaxis_library}"],
                         trendline=trend, trendline_color_override="red",
                         color=color_by, color_continuous_scale='BlueRed', color_discrete_map=color_discrete_map,
                         labels={f"frac_sub_50_tails_{yaxis_library}": "% tails >50"})
        if 'trendline' in options_list:
            trend_results: pd.DataFrame = px.get_trendline_results(fig)
            for index, result in trend_results.iterrows():
                spearman_r, spearman_p = stats.spearmanr(PLOT_DF[extended_xaxis_col],
                                                         PLOT_DF[extended_yaxis_col])
                print_text = f"<b>Correlation & Trend:</b>"
                if not color_by:
                    print_text += f"<br>Trend R<sup>2</sup> = {result[0].rsquared:.4f}"
                print_text += f"<br>Spearman R = {spearman_r:.4f}" \
                              f"<br>Spearman p-val = {spearman_p:.2E}"
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
        fig.update_layout(margin={'l': 40, 'b': 40, 't': 10, 'r': 0},
                          width=600, height=600,
                          hovermode='closest', clickmode="event+select",
                          xaxis_title=f"{xaxis_col} (Lib: {xaxis_library})",
                          yaxis_title=f"{yaxis_col} (Lib: {yaxis_library})",
                          template='plotly_white')

        traces_dict = dict(size=4, color='black', opacity=0.8)
        if color_by_col != "None":
            if color_by_col in PLOT_DF.columns.to_list():
                traces_dict.pop('color')
            elif color_by_col in ["tail_groupings_group", "frac_sub_50_tails"]:
                traces_dict.pop('color')

        fig.update_traces(marker=traces_dict)
        if 'diagonal' in options_list:
            if xaxis_col == yaxis_col:
                min_xy = min(min_plot_x, min_plot_y)
                max_xy = max(max_plot_x, max_plot_y)
                fig.add_trace(go.Scatter(x=[min_xy, max_xy],
                                         y=[min_xy, max_xy],
                                         mode='lines',
                                         line=dict(color='black',
                                                   dash='dash'),
                                         showlegend=False))
            else:
                print(f"Diagonal doesn't really make sense unless you are comparing the same metric between libs!")
        if selected_data['type'] != "Pre":
            selected_df = SELECTED_DF
            fig.add_trace(go.Scatter(mode='markers',
                                     x=selected_df[extended_xaxis_col],
                                     y=selected_df[extended_yaxis_col],
                                     customdata=list(zip(selected_df[f"gene_id"],
                                                         selected_df[f"gene_name"],
                                                         selected_df[extended_xaxis_col],
                                                         selected_df[extended_yaxis_col])),
                                     marker=dict(color='red',
                                                 size=7,
                                                 line=dict(color='black',
                                                           width=2)),
                                     text=selected_df['gene_name'],
                                     hovertemplate="<b>%{text}</b><br>"
                                                   f"<br>{extended_xaxis_col}={'%{x}'}"
                                                   f"<br>{extended_yaxis_col}={'%{y}'}"
                                                   "<br>gene_id=%{customdata[0]}",
                                     name='Selected Genes',
                                     hoverlabel=dict(
                                         bgcolor="black",
                                     ),
                                     showlegend=False))

        fig.update_layout(legend=dict(  # orientation="h",
            yanchor="top", y=0.99,
            xanchor="left", x=0.01,
            itemsizing='constant',
            bordercolor="darkgray",
            borderwidth=2,
            bgcolor="lightgray",
            font=dict(family="Courier New, monospace",
                      size=14), ))
        if 'log_x' in options_list:
            fig.update_xaxes(type='log',
                             range=[log10(min_plot_x), log10(max_plot_x)],
                             )
        else:
            fig.update_xaxes(range=[min_plot_x, max_plot_x])
        if 'log_y' in options_list:
            fig.update_yaxes(type='log',
                             range=[log10(min_plot_y), log10(max_plot_y)],
                             )
        else:
            fig.update_yaxes(range=[min_plot_y, max_plot_y])

        # The below commented out code is trying to get the z-order right;
        #   Based on this link:
        #   https://stackoverflow.com/questions/58278444/how-do-i-control-which-trace-plots-on-top-with-plotly
        # fig.data = fig.data[::-1]

        return fig

    @app.callback(
        [
            Output('selected-datatable', 'children')
        ],
        [
            Input('primary-scatter', 'selectedData'),
            Input('gene_name_select', 'value'),
            Input('xaxis-lib', 'value'),
            Input('yaxis-lib', 'value'),
            Input('btn-nclicks-5', 'n_clicks'),
        ]
    )
    def display_selected_data(selectedData, select_by_name, x_lib, y_lib, reset_clicks):
        global COMPRESSED_DF, PLOT_DF, SELECTED_DF
        changed_id = [p['prop_id'] for p in callback_context.triggered][0]
        if 'btn-nclicks-5' in changed_id:
            SELECTED_DF = None
            return [html.Pre("No points selected")]
        if not selectedData and select_by_name == '':
            return [html.Pre("No points selected")]
        if select_by_name != '':
            try:
                name_selected_df = PLOT_DF.loc[PLOT_DF.gene_name == select_by_name]
            except KeyError:
                print(f"Could not find {select_by_name} in COMPRESSED_DF!!")
                name_selected_df = pd.DataFrame()
                if not selectedData:
                    return [html.Pre("No points selected")]
        else:
            name_selected_df = pd.DataFrame()
        if selectedData:
            selected_df = pd.DataFrame()
            for point_dict in selectedData['points']:
                gene_id, gene_name = point_dict['customdata'][:2]
                if gene_name != select_by_name:
                    selected_df = pd.concat([selected_df, PLOT_DF.loc[PLOT_DF.gene_name == gene_name]])
        else:
            selected_df = pd.DataFrame()
        SELECTED_DF = pd.concat([name_selected_df, selected_df])
        transposed_selected_df = SELECTED_DF.set_index('gene_name').T.reset_index()
        trans_df_columns = [{"name": i, "id": i} for i in transposed_selected_df.columns]
        return [dash_table.DataTable(data=transposed_selected_df.to_dict('records'), columns=trans_df_columns)]

    @app.callback(
        Output('container-button-timestamp', 'children'),
        [Input('btn-nclicks-1', 'n_clicks'),
         Input('btn-nclicks-2', 'n_clicks'),
         Input('primary-scatter', 'figure')])
    def quick_button_click(pop_scatter, save_scatter, scatter_fig):
        changed_id = [p['prop_id'] for p in callback_context.triggered][0]
        if 'btn-nclicks-1' in changed_id:
            msg = 'Popout Scatter'
            pio.show(scatter_fig)
        elif 'btn-nclicks-2' in changed_id:
            msg = 'Quick Save'
            save_file_loc = f"./testOutputs/{get_dt(for_file=True)}_dashScatter.svg"
            go.Figure(scatter_fig).write_image(save_file_loc,
                                               width=600, height=600, scale=2)
        else:
            msg = '*Neither of these buttons*'
        return html.Div(f"{msg} was most recently clicked")

    @app.callback(
        Output('download_image', 'data'),
        [Input('btn-nclicks-3', 'n_clicks'),
         Input('btn-nclicks-4', 'n_clicks'),
         Input('primary-scatter', 'figure')])
    def save_button_click(save_as_svg, save_as_png, scatter_fig):
        changed_id = [p['prop_id'] for p in callback_context.triggered][0]
        if 'btn-nclicks-3' in changed_id:
            temp_save_file_loc = f"./testOutputs/{get_dt(for_file=True)}_dashScatter.svg"
            go.Figure(scatter_fig).write_image(temp_save_file_loc,
                                               width=600, height=600, scale=2)
            temp_hold = dcc.send_file(temp_save_file_loc)
            os.remove(temp_save_file_loc)
            return temp_hold
        elif 'btn-nclicks-4' in changed_id:
            temp_save_file_loc = f"./testOutputs/{get_dt(for_file=True)}_dashScatter.png"
            go.Figure(scatter_fig).write_image(temp_save_file_loc,
                                               width=600, height=600, scale=2)
            temp_hold = dcc.send_file(temp_save_file_loc)
            os.remove(temp_save_file_loc)
            return temp_hold

    @app.callback(Output('short-group-cutoff-label', 'children'),
                  Output('long-group-cutoff-label', 'children'),
                  Input('yaxis-lib', 'value'),
                  Input('short-group-cutoff-slider', 'value'),
                  Input('long-group-cutoff-slider', 'value'),
                  Input('min-hits-slider', 'value'))
    def display_value(y_axis_lib, short_group_cutoff, long_group_cutoff, min_hits):
        COMPRESSED_DF['test_group'] = pd.cut(COMPRESSED_DF[f'frac_sub_50_tails'],
                                                       bins=[0.0, long_group_cutoff, short_group_cutoff, 1.0],
                                                       labels=['long_tailed',
                                                               'ungrouped',
                                                               'short_tailed'],
                                                       include_lowest=True)
        min_hit_df = COMPRESSED_DF.query(f"lib == '{y_axis_lib}'")
        min_hit_df = min_hit_df.query(f"gene_hits >= {min_hits}")
        counts = min_hit_df.value_counts('test_group')
        return (f"Shorter tail group cutoff: ({counts['short_tailed']}; "
                f"{counts['short_tailed'] / counts.sum()*100:0.1f}%)",
                f"Longer tail group cutoff: ({counts['long_tailed']}; "
                f"{counts['long_tailed'] / counts.sum()*100:0.1f}%)")
    
    app.run_server(debug=False, dev_tools_hot_reload=False)


if __name__ == '__main__':
    from sys import argv

    libraries_to_run = argv[1:]
    print(f"Running w/ libraries: {libraries_to_run}")
    main(libraries_to_run,
         force_compressed_df_build=False,
         abs_min_cutoff=1)
