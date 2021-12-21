"""
geneHeatmaps_for5TERA.py
Marcus Viscardi,    December 18, 2021

This is going to be an adaptation of geneHeatmaps2.py, mainly with the goal to
specifically plot 5TERA adapted reads (these should be specifically 5'monoP
containing mRNAs!)
"""
import time

from nanoporePipelineCommon import assign_with_josh_method, pick_libs_return_paths_dict

from geneHeatmaps2 import manual_cdf
from scipy.stats import norm

import numpy as np
import pandas as pd

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)


def f(hit_list, bounds):
    # Somehow this is slower than the crap I wrote for geneHeatmaps2.py lol
    histo, bins = np.histogram(np.clip(hit_list, bounds[0], bounds[1]), bins=range(bounds[0], bounds[1] + 2),
                               density=True)
    return np.cumsum(histo)


def np_repeat_wrapper(hit_list, min_x, max_x, times=10000):
    for _ in range(times):
        f(hit_list, (min_x, max_x))

def man_repeat_wrapper(hit_list, min_x, max_x, times=10000):
    for _ in range(times):
        manual_cdf(hit_list, min_x=min_x, max_x=max_x)

if __name__ == '__main__':
    for length in [1000, 10000, 100000, 1000000, 10000000]:
        hits_list = list(np.random.randn(length))
        min_x = -3
        max_x = 3
        # With tiny lists, this function is slower than the manual method
        start1 = time.perf_counter()
        f(hits_list, (min_x, max_x))
        end1 = time.perf_counter()
        
        # As list get longer, it seems like the numpy one is faster!
        #   The flip seems to be somewhere around 1000,
        #   but below that they are still really comparable:
        #       for 1,000 items:
        #           Numpy:  ~3.06E-4 sec
        #           Manual: ~2.50E-4 sec
        #       for 10,000 items:
        #           Numpy:  ~1.07E-3 sec
        #           Manual: ~2.78E-3 sec
        start2 = time.perf_counter()
        manual_cdf(hits_list, min_x, max_x)
        end2 = time.perf_counter()
        
        print(f"List length of {len(hits_list)}:")
        print("\tnumpy: ", end1 - start1)
        print("\tmanual:", end2 - start2)
    
