"""
geneAnnotationPlotting.py
Marcus Viscardi,    September 22, 2021

I want to plot gene annotation information so that I can make some
    choices regarding the metaSTOP for nanopore stuff.
"""

import pandas as pd
import seaborn as sea
import matplotlib.pyplot as plt
GENE_ANNOTATION_PATH = "/data16/marcus/genomes/elegansRelease100/" \
                       "Caenorhabditis_elegans.WBcel235.100.gtf.parquet"


def load_3_utrs():
    df = pd.read_parquet(GENE_ANNOTATION_PATH)
    df = df[df.gene_biotype == "protein_coding"]
    for thing in ["transcript", "stop_codon", "start_codon", "three_prime_utr"]:
        print(thing.title(), df[df.feature == thing].shape[0])
    df = df[df.feature == "three_prime_utr"]
    df["3' UTR length"] = df["end"] - df["start"]
    return df


def seaborn_cdf(df, columns_to_plot):
    for column_to_plot in columns_to_plot:
        fig = sea.ecdfplot(data=df, x=column_to_plot)
        fig.set(xlim=(-10, 1000))
        plt.show()


if __name__ == '__main__':
    seaborn_cdf(load_3_utrs(), ["3' UTR length"])
