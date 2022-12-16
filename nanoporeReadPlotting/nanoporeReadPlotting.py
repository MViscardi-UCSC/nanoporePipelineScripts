"""
nanoporeReadPlotting.py
Marcus Viscardi,    July 05, 2022

Script to pull together the methods written in notebooks
    into a long term script, that will have a CLI interface!

Much of this will come from:
    - SimpleReadPlotting_cigars.py
    - finalizingReadPlotting_matplitlib.ipynb
    - testingReadPlotting_matplotlib.ipynb
    - testingReadPlotting_plotly.ipynb

Hopefully this will also integrate the plotting of the reference!
"""

from nanoporePipelineCommon import *

import seaborn as sea
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

import numpy as np
import pandas as pd
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', None)

