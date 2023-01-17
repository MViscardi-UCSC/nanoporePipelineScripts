"""
dashPlotLibCounts.py
Marcus Viscardi,    December 27, 2022

Setting this up to hopefully be usable for other libraries and plot types.

Main going to plot a permutation matrix for X number of libraries that shows
all the different comparisons on one Dash/Plotly grid. Then, if you select a
gene on one plot, it gets highlighted in all other plots!
"""
import seaborn as sea
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

import numpy as np
import pandas as pd
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', None)

if __name__ == '__main__':
    target_prefix = "total_rpm"

    df = pd.read_parquet(f"221227_quad5TERA.counts.parquet")
    # Scatter matrix already allows us to select some points and have them selected in all plots.
    #   Is this all I needed...?!
    fig = px.scatter_matrix(df,
                            dimensions=[col for col in df.columns
                                        if col.startswith(target_prefix)],
                            labels={col: col.replace(f"{target_prefix}_", "") for col in df.columns},
                            )
    fig.write_html(f"221227_testing.html")
    fig.show()
