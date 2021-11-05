"""
geneMeta5end_CDSnormed.py
Marcus Viscardi,    November 04, 2021

This will plot something similar to geneMeta5Ends.py, but it will
    rather have the x-axis be % of CDS from the start codon to the stop,
    plus some distance on either side to accommodate UTRs
"""
import pandas as pd

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)
import seaborn as sea
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from nanoporePipelineCommon import find_newest_matching_file


def calc_percent_cds(reads_w_josh_assignment_df: pd.DataFrame) -> pd.DataFrame:
    print("hi")
    reads_w_josh_assignment_df["cds_length"] = reads_w_josh_assignment_df["to_start"] - \
                                               reads_w_josh_assignment_df["to_stop"]
    # negative to_stop values are in CDS, negative to_start values are in the 3'UTR
    reads_w_josh_assignment_df["percent_of_cds"] = - reads_w_josh_assignment_df["to_stop"] /\
                                                   reads_w_josh_assignment_df["cds_length"]
    return reads_w_josh_assignment_df


def load_parquet(outputDir):
    parquet_path = find_newest_matching_file(f"{outputDir}/merge_files/*mergedOnReadsPlus.parquet")
    return pd.read_parquet(parquet_path)


if __name__ == '__main__':
    df = calc_percent_cds(
        load_parquet(f"/data16/marcus/working/210720_nanoporeRun_totalRNA_0639_L3_replicate/output_dir"))
