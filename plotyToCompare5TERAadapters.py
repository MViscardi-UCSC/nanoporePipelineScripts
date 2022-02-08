"""
plotyToCompare5TERAadapters.py
Marcus Viscardi,    December 12, 2021

Goal is to make a scatter between the two 5TERA libs
"""
import plotly.express as px
import scipy.stats

from nanoporePipelineCommon import pick_libs_return_paths_dict, gene_names_to_gene_ids, load_ski_pelo_targets

import pandas as pd
import numpy as np

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


def load_libraries(per_gene_cutoff=25, ski_pelo_col=True):
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
    print(df_dict['wt'].columns)
    merge_df = pd.merge(df_dict["smg-6"], df_dict["wt"], on=["gene_id", "gene_name"],
                        suffixes=("_6", "_wt"))
    merge_df['total_read_counts'] = merge_df['read_hits_6'] + merge_df['read_hits_wt']
    merge_df['t5_fraction_diff'] = (merge_df['t5_fraction_wt'] - merge_df['t5_fraction_6'])
    merge_df['t5_fraction_mean'] = (merge_df['t5_fraction_wt'] + merge_df['t5_fraction_6']) / 2
    merge_df['t5_fraction_abs_diff'] = merge_df['t5_fraction_diff'].abs()

    merge_df['t5_read_count_wt'] = merge_df['t5_fraction_wt'] * merge_df['read_hits_wt']
    merge_df['t5_read_count_6'] = merge_df['t5_fraction_6'] * merge_df['read_hits_6']
    merge_df['t5_total_reads'] = merge_df['t5_read_count_wt'] + merge_df['t5_read_count_6']
    merge_df['fraction_of_t5_from_wt'] = merge_df['t5_read_count_wt'] / merge_df['t5_total_reads']
    merge_df['per_gene_p_value'] = merge_df['read_hits_wt'] / merge_df['total_read_counts']

    t5_total_wt = merge_df['t5_read_count_wt'].sum()
    t5_total_6 = merge_df['t5_read_count_6'].sum()
    total_fraction_t5 = t5_total_wt / (t5_total_wt + t5_total_6)
    merge_df['bionomial_test'] = merge_df.apply(lambda row: scipy.stats.binom_test(
        x=row['t5_read_count_wt'],
        n=row['t5_total_reads'],
        p=total_fraction_t5),
                                                axis=1)
    merge_df['bionomial_test_per_gene_p'] = merge_df.apply(lambda row: scipy.stats.binom_test(
        x=row['t5_read_count_wt'],
        n=row['t5_total_reads'],
        p=row['per_gene_p_value']),
                                                           axis=1)
    merge_df['bionomial_test_nlog10'] = np.log10(merge_df['bionomial_test']) * -1
    merge_df['bionomial_test_nlog10_per_gene_p'] = np.log10(merge_df['bionomial_test_per_gene_p']) * -1
    print(t5_total_wt, t5_total_6, total_fraction_t5)

    if ski_pelo_col:
        ski_pelo_list = load_ski_pelo_targets()
        merge_df['ski_pelo_targ'] = merge_df.gene_id.isin(ski_pelo_list)

    for suffix in ("_6", "_wt"):
        merge_df = merge_df[merge_df[f'read_hits{suffix}'] >= per_gene_cutoff]
    print(merge_df.query("gene_id == 'WBGene00001340'"))
    return merge_df


def with_dash_for_click_to_copy():
    import dash
    import json
    from dash import dcc, html, callback_context
    import dash_bootstrap_components as dbc
    from dash.dependencies import Input, Output
    import pyperclip
    from webbrowser import open as open_webpage
    # pyperclip.copy('The text to be copied to the clipboard.')

    merge_df = load_libraries(per_gene_cutoff=1)

    app = dash.Dash(__name__)

    app.layout = dbc.Container(
        [
            dbc.Row(
                # Top row w/ buttons and slider
                [
                    dbc.Col(
                        [
                            dcc.Markdown(
                                """**Select minimum reads/gene to allow**""",
                            ),
                            dcc.Slider(
                                id='min-hits-slider',
                                min=5, max=100,
                                value=40,
                                marks={str(n): str(n) for n in range(5, 105, 5)},
                                step=None
                            ),
                            dcc.Dropdown(
                                id='x_col',
                                options=[{'label': i, 'value': i} for i in merge_df.columns],
                                value="total_read_counts"),
                            dcc.Dropdown(
                                id='y_col',
                                options=[{'label': i, 'value': i} for i in merge_df.columns],
                                value="fraction_of_t5_from_wt"),
                        ],
                        width=8,
                        style={'height': '100%'}
                    ),
                    dbc.Col(
                        [
                            html.Button(
                                'Copy selected GENE_NAME',
                                id='btn-nclicks-2',
                                n_clicks=0),
                            html.Button(
                                'Copy selected WBGENE_ID',
                                id='btn-nclicks-1',
                                n_clicks=0),
                            html.Button(
                                'Open WormBase',
                                id='btn-nclicks-3',
                                n_clicks=0),
                            html.Div(
                                id='container-button-timestamp'
                            ),
                            html.Div(
                                id='selected-data'
                            ),
                        ],
                        width=4,
                        style={'height': '100%'}
                    ),
                ],
                className="h-10",
            ),
            dbc.Row(
                # Row for graph!
                [
                    dcc.Graph(id='primary-scatter',
                              style=
                              {
                                  'width': '120vh',
                                  'height': '80vh',
                              }
                              ),
                ],
                className='h-60',
            ),
        ],
        style={'height': '100%', 'width': '100%'},
    )

    @app.callback(Output('primary-scatter', 'figure'),
                  [Input('min-hits-slider', 'value'),
                   Input('x_col', 'value'),
                   Input('y_col', 'value'), ])
    def main_plot(min_hits, x_col, y_col):
        plot_df = merge_df[merge_df["read_hits_wt"] >= min_hits].copy(deep=True)
        plot_df = plot_df[plot_df["read_hits_6"] >= min_hits]
        print(f"There are {plot_df.shape[0]} genes plotted.")
        fig = px.scatter(plot_df,
                         x=x_col,
                         y=y_col,
                         # x="total_read_counts",
                         # y='bionomial_test_nlog10',
                         # y="t5_fraction_diff",
                         # x='t5_total_reads', y='fraction_of_t5_from_6',
                         color="fraction_of_t5_from_wt",
                         color_continuous_scale='Bluered',
                         hover_name="gene_name", hover_data=["gene_id",
                                                             "read_hits_6",
                                                             "read_hits_wt",
                                                             "t5_fraction_diff",
                                                             "t5_fraction_wt",
                                                             "t5_fraction_6"],
                         custom_data=["gene_id", "gene_name"])
        fig.update_layout(margin={'l': 40, 'b': 40, 't': 10, 'r': 0},
                          hovermode='closest', clickmode="event+select",
                          template='plotly_white')
        fig.update_xaxes(type='log')
        # fig.update_yaxes(type='log')
        return fig

    @app.callback(
        Output('selected-data', 'children'),
        [Input('primary-scatter', 'selectedData')])
    def display_selected_data(selectedData):
        if not selectedData:
            return "No points selected"
        return_data = []
        for point_dict in selectedData['points']:
            gene_id, gene_name = point_dict['customdata'][0:2]
            # print(point_dict['customdata'])
            return_data.append({'Gene Name': gene_name,
                                'Gene ID': gene_id})
        return json.dumps(return_data, indent=2)

    @app.callback(
        Output('container-button-timestamp', 'children'),
        [Input('btn-nclicks-1', 'n_clicks'),
         Input('btn-nclicks-2', 'n_clicks'),
         Input('btn-nclicks-3', 'n_clicks'),
         Input('primary-scatter', 'figure'),
         Input('primary-scatter', 'selectedData')])
    def button_click(gene_id, gene_name, open_worm_base, scatter_fig, selectedData):
        changed_id = [p['prop_id'] for p in callback_context.triggered][0]
        if not selectedData:
            return_data = {'gene_name': [],
                           'gene_id': []}
        else:
            return_data = {'gene_name': [],
                           'gene_id': []}
            for point_dict in selectedData['points']:
                gene_id, gene_name = point_dict['customdata'][0:2]
                return_data["gene_name"].append(gene_name)
                return_data["gene_id"].append(gene_id)
        if 'btn-nclicks-2' in changed_id:
            msg = 'Copy gene_id'
            # print(scatter_fig)
            if len(return_data['gene_name']) == 1:
                pyperclip.copy(return_data['gene_name'][0])
            elif len(return_data['gene_name']) > 1:
                pyperclip.copy(', '.join(return_data['gene_name']))
        elif 'btn-nclicks-1' in changed_id:
            msg = 'Copy gene_name'
            if len(return_data['gene_id']) == 1:
                pyperclip.copy(return_data['gene_id'][0])
            elif len(return_data['gene_id']) > 1:
                pyperclip.copy(', '.join(return_data['gene_id']))
        elif 'btn-nclicks-3' in changed_id:
            msg = 'Open worm base'
            if len(return_data['gene_id']) == 1:
                wbgene_id = return_data['gene_id'][0]
                open_webpage(f"https://wormbase.org/species/c_elegans/gene/{wbgene_id}", new=1)
            elif len(return_data['gene_id']) > 1:
                for i, wbgene_id in enumerate(return_data['gene_id']):
                    if i == 0:
                        open_webpage(f"https://wormbase.org/species/c_elegans/gene/{wbgene_id}", new=1)
                    else:
                        open_webpage(f"https://wormbase.org/species/c_elegans/gene/{wbgene_id}", new=2)
        else:
            msg = '*None of the buttons*'
        return html.Div(f"{msg} was most recently clicked")

    app.run_server(debug=False, dev_tools_hot_reload=False, host='127.0.0.2', port=8000)


if __name__ == '__main__':
    # quick_test()
    with_dash_for_click_to_copy()
