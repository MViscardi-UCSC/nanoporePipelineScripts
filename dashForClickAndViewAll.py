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

from nanoporePipelineCommon import get_dt, find_newest_matching_file
from dashForClickAndViewPolyATails import load_and_merge_lib_parquets
import pandas as pd

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)


def main(libs, force_compressed_df_build=False, abs_min_cutoff=5):
    global COMPRESSED_DF, PLOT_DF
    import dash
    import json
    from dash import dcc, html, callback_context
    import dash_bootstrap_components as dbc
    from dash.dependencies import Input, Output
    import dash_daq as daq
    import plotly.express as px
    import plotly.graph_objects as go
    import plotly.io as pio

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
        reads_df, COMPRESSED_DF = load_and_merge_lib_parquets(libs, drop_sub_n=abs_min_cutoff)
        save_path = f"./testInputs/{get_dt(for_file=True)}_{'_'.join(libs)}.compressed.parquet"
        COMPRESSED_DF.to_parquet(save_path)
        print(f"Saved new compressed file to: {save_path}")
    lib_list = COMPRESSED_DF.lib.unique().tolist()
    plottable_columns = COMPRESSED_DF.select_dtypes(include='number').columns.to_list()

    PLOT_DF = pd.DataFrame()

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
                        ], style={'padding': 10, 'flex': 4}
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
                                    {'label': 'Color By Chr', 'value': 'color_chr'},
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
                                min=5*(abs_min_cutoff//5), max=100,
                                value=40, persistence=True, persistence_type='local',
                                marks={str(n): str(n) for n in range(5*(abs_min_cutoff//5), 100+abs_min_cutoff, 5)},
                                step=None
                            ),
                        ], style={'padding': 10, 'flex': 4}
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
                        ], style={'padding': 10, 'flex': 4}
                    ),
                ], style={'display': 'flex', 'flex-direction': 'row'},
            ),
            html.Div(  # Plotting Row
                [
                    html.Div(  # Left section of plotting row
                        [
                            # html.H3(
                            #     'Rocket Plot of Mean Tail Lengths'
                            # ),
                            dcc.Graph(
                                id='primary-scatter',
                            )
                        ], style={'padding': 10, 'flex': '0 0 600px'}
                    ),
                    html.Div(  # Right section of the plotting row, with selected-data
                        [
                            html.Pre(id='selected-data')
                        ], style={'padding': 10, 'flex': 1}
                    )
                ], style={'display': 'flex', 'flex-direction': 'row'},
            ),
            html.Div(  # Saving/Output Buttons Row
                [
                    html.Div(
                        [
                            html.Button('Popout Scatter', id='btn-nclicks-1', n_clicks=0),
                            html.Button('Quick Save', id='btn-nclicks-2', n_clicks=0)
                        ]
                    ),
                    html.Div(id='container-button-timestamp'),
                    html.Div(
                        [
                            html.Button('Save Scatter As SVG', id='btn-nclicks-3', n_clicks=0),
                            html.Button('Save Scatter As PNG', id='btn-nclicks-4', n_clicks=0),
                            dcc.Download(id='download_image'),
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
         Input('xaxis-col', 'value'),
         Input('yaxis-col', 'value'),
         Input('min-hits-slider', 'value'),
         Input('selected-data', 'children'),
         Input('options', 'value')])
    def main_plot(xaxis_library, yaxis_library,
                  xaxis_col, yaxis_col,
                  min_hits,
                  selectedData, options_list) -> go.Figure:
        extended_xaxis_col = f"{xaxis_col}_{xaxis_library}"
        extended_yaxis_col = f"{yaxis_col}_{yaxis_library}"
        if selectedData != "No points selected":
            selected_gene_ids = [hit["Gene ID"] for hit in json.loads(selectedData)]
            selected_gene_names = [hit["Gene Name"] for hit in json.loads(selectedData)]
        else:
            selected_gene_ids, selected_gene_names = [], []

        # Below tries to make some "smart" bounds around the plot!
        min_hit_df = COMPRESSED_DF[COMPRESSED_DF['gene_hits'] >= min_hits]

        max_plot_x = min_hit_df[min_hit_df['lib'] == xaxis_library][xaxis_col].max()
        min_plot_x = min_hit_df[min_hit_df['lib'] == xaxis_library][xaxis_col].min()

        max_plot_y = min_hit_df[min_hit_df['lib'] == yaxis_library][yaxis_col].max()
        min_plot_y = min_hit_df[min_hit_df['lib'] == yaxis_library][yaxis_col].min()

        # TODO: Something to handle if x and y libs are the same (just make one df and copy it?)
        x_axis_df = min_hit_df[min_hit_df.lib == xaxis_library].drop(columns='lib')
        y_axis_df = min_hit_df[min_hit_df.lib == yaxis_library].drop(columns='lib')

        PLOT_DF = pd.merge(x_axis_df, y_axis_df, on=["gene_id", "gene_name", "chr_id"],
                           suffixes=(f"_{xaxis_library}",
                                     f"_{yaxis_library}"))

        if 'trendline' in options_list:
            trend = "ols"
        else:
            trend = None
        if 'color_chr' in options_list:
            color_by = 'chr_id'
        else:
            color_by = None
        fig = px.scatter(PLOT_DF,
                         x=extended_xaxis_col,
                         y=extended_yaxis_col,
                         custom_data=["gene_id", "gene_name"],
                         hover_name="gene_name", hover_data=["gene_id",
                                                             f"gene_hits_{xaxis_library}",
                                                             f"gene_hits_{yaxis_library}"],
                         trendline=trend, trendline_color_override="red",
                         color=color_by)
        fig.update_layout(margin={'l': 40, 'b': 40, 't': 10, 'r': 0},
                          width=600, height=600,
                          hovermode='closest', clickmode="event+select",
                          xaxis_title=f"{xaxis_col} (Lib: {xaxis_library})",
                          yaxis_title=f"{yaxis_col} (Lib: {yaxis_library})",
                          template='plotly_white')

        traces_dict = dict(size=4, color='black', opacity=0.9)
        if 'color_chr' in options_list:
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
        if selectedData != "No points selected":
            selected_df = PLOT_DF[PLOT_DF["gene_id"].isin(selected_gene_ids)]
            fig.add_trace(go.Scatter(mode='markers',
                                     x=selected_df[extended_xaxis_col],
                                     y=selected_df[extended_yaxis_col],
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

        fig.update_layout(legend=dict(  # orientation="h",
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01))
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
        Output('selected-data', 'children'),
        [Input('primary-scatter', 'selectedData'),
         Input('gene_name_select', 'value'),
         Input('xaxis-lib', 'value'),
         Input('yaxis-lib', 'value')])
    def display_selected_data(selectedData, select_by_name, x_lib, y_lib):
        if not selectedData and select_by_name == '':
            return "No points selected"
        return_data = []
        if select_by_name != '':
            try:
                selected_df = COMPRESSED_DF.loc[COMPRESSED_DF.gene_name == select_by_name].reset_index().set_index(
                    'lib')
                return_data = [{'Gene Name': select_by_name,
                                'Gene ID': selected_df.loc[x_lib, 'gene_id'],
                                f'{x_lib} hits': str(selected_df.loc[x_lib, 'gene_hits']),
                                f'{y_lib} hits': str(selected_df.loc[y_lib, 'gene_hits'])}]
            except KeyError:
                print(f"Could not find {select_by_name} in COMPRESSED_DF!!")
                return_data = []
                if not selectedData:
                    return "No points selected"
        if selectedData:
            # pprint(selectedData)
            for point_dict in selectedData['points']:
                gene_id, gene_name, x_counts, y_counts = point_dict['customdata']
                if gene_name != select_by_name:
                    return_data.append({'Gene Name': gene_name,
                                        'Gene ID': gene_id,
                                        f'{x_lib} hits': x_counts,
                                        f'{y_lib} hits': y_counts})
        return json.dumps(return_data, indent=2)

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
            msg = '*None of the buttons*'
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
    
    app.run_server(debug=False, dev_tools_hot_reload=False)


if __name__ == '__main__':
    from sys import argv

    libraries_to_run = argv[1:]
    print(f"Running w/ libraries: {libraries_to_run}")
    main(libraries_to_run,
         #force_compressed_df_build=True,
         abs_min_cutoff=1)
