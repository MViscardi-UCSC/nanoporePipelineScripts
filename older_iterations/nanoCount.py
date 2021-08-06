"""
nanoCount.py
Marcus Viscardi     July 21, 2021

Messing w/ nanocount to see if I can get it to work for the heatmap stuff
"""
from step0_nanopore_pipeline import minimap_and_nanopolish
from NanoCount.NanoCount import NanoCount
import pandas as pd

if __name__ == '__main__':
    # This doesn't work. Maybe because of the minimap2 settings I used??
    
    # An option would be to choose top transcripts/isoforms with a tool like
    #   NanoCount, then just running those top isoforms per gene for genes!
    path_to_bam = "/data16/marcus/prefix/210709_NanoporeRun_totalRNA_0639_L3/output_dir/cat_files/cat.sorted.bam"
    df = NanoCount(alignment_file=path_to_bam,
                   filter_bam_out="./aligned_reads_selected.bam").count_df
    print(df)
    print(df.head())
    print(df.info())
