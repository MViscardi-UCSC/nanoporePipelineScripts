"""
nanoporeDirectRNASeq_QC.py
Marcus Viscardi,    April 12, 2022

General goal is to actually make a QC script that will automatically
    (and consistently) spit out some of the things I have been putting
    into tables.
"""
import inspect

import seaborn as sea
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

from nanoporePipelineCommon import *

import numpy as np
import pandas as pd

pd.set_option('display.width', 200)
pd.set_option('display.max_columns', None)


def bases_called_and_mapped(dataframe=None, lib=None, bam_or_sam='bam'):
    if isinstance(dataframe, pd.DataFrame):
        bam_df = dataframe.copy()
    elif isinstance(lib, str):
        bam_path = pick_lib_return_path(lib,
                                        output_dir_folder='cat_files',
                                        file_midfix='cat.sorted.mappedAndPrimary',
                                        file_suffix=bam_or_sam)
        bam_obj = SamOrBamFile(bam_path)
        bam_df = bam_obj.df
    else:
        raise NotImplementedError(f"Please give a pandas dataframe or a "
                                  f"library name to this method "
                                  f"({inspect.currentframe().f_code.co_name})")
    num_bases = bam_df.sequence.apply(len).sum()
    if isinstance(lib, str):
        print('\n\n', lib, '\nNumber of bases sequenced and mapped: ', num_bases, '\n')
    return num_bases


if __name__ == '__main__':
    polyA3_obj = SamOrBamFile(pick_lib_return_path('polyA3',
                                                   output_dir_folder='cat_files',
                                                   file_midfix='cat.sorted.mappedAndPrimary',
                                                   file_suffix='bam'))
    totalRNA3_obj = SamOrBamFile(pick_lib_return_path('totalRNA3',
                                                      output_dir_folder='cat_files',
                                                      file_midfix='cat.sorted.mappedAndPrimary',
                                                      file_suffix='bam'))
    print(f"polyA3: number of bases called and mapped: {bases_called_and_mapped(dataframe=polyA3_obj.df)}")
    print(f"totalRNA3: number of bases called and mapped: {bases_called_and_mapped(dataframe=totalRNA3_obj.df)}")
    print('Done.')
