"""
SimpleReadPlotting_cigars.py
Marcus Viscardi,    May 23, 2022

Going to see if I can make read plots just with matplotlib

Maybe I'll also try with plotly? Actually upon looking,
plotly might be easier for the first pass

########

Some notes after making the notebook proof of concept.
1. Plotly has a hard time w/ that many objects being drawn? At least in the way that I am implimenting it!!
    Potential solutions to this:
     - Try to "dumb down" the cigar strings ahead of time, ie. compress short indels to just maps (like below).
       This should greatly reduce the number of objects that plotly has to render on screen!
            25M2D3M5I15M --> 50M
     - Load all of the objects into a list format.
       (based on second section below this: https://plotly.com/python/shapes/#shapedrawing-with-scatter-traces)
       This might make the number of objects far less, without needing to reduce complexity?
2. Try matplotlib... it should be better at this kind of thing!!
   The plotly API just looked easy to implement, and it was! but whatever was happening behind the scenes wasn't fast...

Some notes after doing basically the same thing in matplotlib:
1. WAYYYYY faster, and there were some tricks to drawing all objects in one go rather than iteratively, v cool.
2. Currently, it is FULLY functional for reads. But I really want to implement drawing the references from GTF

Now, I am going to bring over the code from  testingReadPlotting_matplotlib.ipynb!!
"""
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

from nanoporePipelineCommon import *

import re
import numpy as np
import pandas as pd
from tqdm import tqdm
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', None)

def make_rectangle_patch(genome_start, length, y_center, thickness, color='gray'):
    return Rectangle((genome_start, y_center-(thickness/2)), length, thickness,
                     facecolor=color,
                     edgecolor=color,
                     fill=True,
                     lw=0)


def plot_cigar_and_genome_start(axes, cigar, gen_start, y, color='black', tail_length=None):
    # Parse the cigar string
    parsed_cigar = re.findall(rf'(\d+)([MDNSIX])', cigar)
    mdn_nums = [int(num) for num, char in parsed_cigar if char in "MDN"]
    read_end = gen_start + sum(mdn_nums)

    genome_loc = gen_start

    rectangle_patch_list = []
    first_n_length = None
    for length, code in parsed_cigar:
        length = int(length)
        if code == 'S':
            if len(rectangle_patch_list) == 0:
                # genome_loc += length  # <- this was dumb and wrong...
                pass
        elif code == 'M':
            rectangle_patch_list.append(make_rectangle_patch(genome_loc, length, y, thickness=0.8))
            genome_loc += length
        elif code == 'D':
            if length < 50:
                rectangle_patch_list.append(make_rectangle_patch(genome_loc, length, y, thickness=0.8))
            else:
                rectangle_patch_list.append(make_rectangle_patch(genome_loc, length, y, thickness=0.001))
            genome_loc += length
        elif code == 'I':
            pass
        elif code == 'N':
            if not first_n_length:
                first_n_length = length
                # print(f"First N of length {length}, starting at {genome_loc}")
            rectangle_patch_list.append(make_rectangle_patch(genome_loc, length, y, thickness=0.001))
            genome_loc += length
    if read_end - gen_start < 1500:
        axes.add_collection(PatchCollection(rectangle_patch_list, color=color))
        if isinstance(tail_length, float):
            axes.add_patch(make_rectangle_patch(genome_loc, tail_length, y, thickness=0.4, color='green'))
            genome_loc += tail_length
        return read_end - gen_start, first_n_length
    else:
        return None


def row_apply_plot_cigar(row, axes):
    index = row.name
    cigar = row.cigar
    gen_start = row.chr_pos
    is_adapted = row.t5
    polya_length = row.polya_length

    if is_adapted == '-':
        color = 'black'
    else:
        color = 'red'
    return plot_cigar_and_genome_start(axes, cigar, gen_start, index, color=color, tail_length=polya_length)


def plot_gtf_annotation():
    pass