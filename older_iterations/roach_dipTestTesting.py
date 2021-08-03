#!/usr/bin/python3
"""
roach_dipTestTesting.py
Marcus Viscardi     Jan 25, 2021

This script is meant to test different adjustments of the implementation of Hartigan's dip test in the UniDip package
"""
from roach_combineMethodResults_run import plot_some_top_hits
from roach_compareReplicates_run import load_comp_on_genes

ABS_PATH = f"/data16/marcus/working/210106_nanopolish_wRoachData"
OUTDIR1 = f"{ABS_PATH}/output_dir_adultrep1"
OUTDIR2 = f"{ABS_PATH}/output_dir_adultrep2"
genes_file = f"{OUTDIR1}/210124_compressed_on_genes.tsv"


if __name__ == '__main__':
    df = load_comp_on_genes(genes_file)
    genes_per_plot = 30
    how_many_plots = 5
    for i in range(0, genes_per_plot*how_many_plots, genes_per_plot):
        df_to_plot = df[["gene_id", "read_hits", "polya_mean", "polya_lengths"]].iloc[i:i+genes_per_plot]
        plot_some_top_hits(df_to_plot, how_many_plots=genes_per_plot, series_named=False,
                           x_limit=200, title=f"Hits[{i}:{i+genes_per_plot}]'s Tail Lengths")
