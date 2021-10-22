"""
biases_wMtDNAgenes.py
Marcus Viscardi,    October 21, 2021

The goal here is to try and figure out what is going on with the MtDNA genes
that keep coming off diagonal between my polyA libraries and totalRNA.

The first step is going to be producing plots that cover the transcripts of
the MtDNA chromosome genes, and show the local T-fraction, to see if that
differs between the MtDNA genes that are coming off diagonal and those that
are not!

To do this I am going to turn the cDNA sequences of these genes into arrays
of nucleotides, then convert all non-T nucleotides to 0.0 and the T nucleotides
to be 1.0. By taking a rolling average across this array I should find local T
fractions.
"""
import pandas as pd
import numpy as np

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)

import seaborn as sea
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import mappy as mp
from tqdm import tqdm


def load_cdna_lib(cdna_path: str, filter_for_chr: str = "MtDNA",
                  subsample: int = None) -> pd.DataFrame:
    # Make an empty dataframe to store transcripts to:
    gc_df = pd.DataFrame()

    # B/C the mappy fastx reader doesn't have a built in __len__ value
    #   we have to feed that to tqdm separately so that it know the
    #   eventual total values it will iterate through:
    row_iterator = tqdm(mp.fastx_read(cdna_path, read_comment=True),
                        total=len(list(mp.fastx_read(cdna_path))))

    # Use the iterator we built to step through each cDNA in the
    #   fasta file:
    for i, (name, seq, _, comment) in enumerate(row_iterator):
        # Parse the comment line into a list:
        comment = comment.split(" ")

        # The chromosome id is in the second item of the comment list.
        #   Use this to decide if it's necessary to parse the rest!
        chr_id = comment[1].split(":")[2]
        if isinstance(filter_for_chr, str) and filter_for_chr == chr_id:
            gene_name = comment[5].split(":")[1]
            gene_id = comment[2].split(":")[1].split(".")[0]
            transcript_id = name
            main_dict = dict(chr_id=chr_id, gene_name=gene_name, gene_id=gene_id,
                             transcript_id=transcript_id, sequence=seq)
            gc_df = gc_df.append(main_dict,
                                 ignore_index=True)
        elif not isinstance(filter_for_chr, str):
            gene_name = comment[5].split(":")[1]
            gene_id = comment[2].split(":")[1].split(".")[0]
            transcript_id = name
            main_dict = dict(chr_id=chr_id, gene_name=gene_name, gene_id=gene_id,
                             transcript_id=transcript_id, sequence=seq)
            gc_df = gc_df.append(main_dict,
                                 ignore_index=True)
        # The "i" in the loop allows us to subsample reads if necessary
        if subsample and i > subsample:
            break

    # Manually add the two mt-rRNAs:
    MTCE33_seq = "TAAATTTTACTTGTTTAAATATTGGTATTGCATACTATTACAATAAAATTTCATGTTAATGAAAAATAGAAACAAAGGGTAGAGTAAATATTAGTTT" \
                 "TATTGTTTCATACTAAAAATTATATTTATTAGAGTTGATATGTCGACCTTTGTGATAACTGTTTTTATTTTTATATTAGAAAATTATATATTATATA" \
                 "ATTATTTTAGGAAATTTAAAATTTGAAGTGTTTTAAATTTATGTTTTACAACATTTTCCTAATTTTATTTAAGTTTAATTTTTAATTTAATAAAGTT" \
                 "TTATTAAATAAATAATTTGTAAATTAGTAAATTTTATAAATTTAATTTATTATTAAAATATAATTGAAGAACTTGAAGTCTTGATCAAATGTTTTTT" \
                 "AAAGACTTAGGCTTTATATTAAAGCTGGCTTCTGCCCTATGATATTTAAATGGCAGTCTTAGCGTGAGGACATTAAGGTAGCAAAATAATTTGTGCT" \
                 "TTTATTGAGTTCCAGTATGAATGAAGTTATTGGTTAGTTCTATTTATGTTTTATGTTTGAATTTAATTTTTATTTAAGAAAAAATAAATATATTTAT" \
                 "ACAAAGATAAGTCTTCGGAAATTCTGTTATTACACAATTAAATAATTGTGTAATAAATTTTCTAGGGCAGAATATTATATAATAGTATTTCACTATA" \
                 "TTTAATTTAAAGAATTACTCCGGAGTTAACAGAAAATCATACCTAATCTAGTACTTATAGTAAGGTAAGTTTTACATCGATGTTGTATTCAGATAAT" \
                 "CTAAGAGAGGAGAAGGCTTAGTAGTTTAGACTGTTCTTCTATTAATTAATCTGACGTGATATTAGTTTAATTCATTGTGAGATAGAATTGTTTATCT" \
                 "TGATAAATATTTATATTTAATACATTTAGTACGAAAGGAACATTGTAAAAGTTTTAAACTTTAAAGATTTTGAAATCTT"
    gc_df = gc_df.append(dict(chr_id="MtDNA", gene_name="MTCE.33", gene_id="WBGene00014472",
                              transcript_id="MTCE.33.1", sequence=MTCE33_seq),
                         ignore_index=True)
    MTCE7_seq = "TAAAGTTTTCTTTCAGGGAATTAAAATTTGATCATGGTTTAAGATGATTTAAAATGGTATTATCTAAATTTGATTTACAGAGTAGGCAATAAAAATTT" \
                "ACCTCGGCAATTTATCGCTTGTAAAATACTTGTTCCAGAATAATCGGCTAGACTTGTTAAAGCTTGTACTTTAATTGATGTTAATTATGAAATTATTA" \
                "TATTTTCTTTTAGATCTATGGTAGAATTTGGATTTATATTAGTGAATTTTCATAATTTTAAGATTTGTTGAACAAAGCAGATTAGTACCTGGTTAGAC" \
                "AAAAATTAAAAGAGCAGGAGTAAAGTTGTATTTAAACTGAAAAGATATTGGCAGACATTCTAAATTATCTTTGGAGGCTGAGTAGTAACTGAGAACCC" \
                "TCATTAACTACTTAATTTTTTGACTCGTGTATGATCGTTTATTTTATTCTTAAGGATTATAATAAAAAATTTTTAATTTATTAAAATAGATATATACC" \
                "CGGTTTATGATTTAAGAAACATTTGGCCTACAATATTTTATATTATGGATTTTAGTTTTAGTTAACTAAATGAAATTGTAAAAGACAGTAAAAAATTC" \
                "TTAATGTATTTTTGAAGATTATCTAGAAGTGGTACAAATCATCCATCAATTGCCCAAAGGGGAGTAAGTTGTAGTAAAGTAGATTTAGGGGAACCTGA" \
                "ATCTAGTAAT"
    gc_df = gc_df.append(dict(chr_id="MtDNA", gene_name="MTCE.7", gene_id="WBGene00014454",
                              transcript_id="MTCE.7.1", sequence=MTCE7_seq),
                         ignore_index=True)
    return gc_df


def explode_cdna_df_and_array_convert(cdna_df: pd.DataFrame, window_size=10) -> pd.DataFrame:
    # seq_df_dict = {}
    # row_iterator = tqdm(cdna_df.to_dict(orient='records'))
    # for row_dict in row_iterator:

    # Lets try converting the sequence column to a list then using explode
    #   rather than manually looping through stuff:
    cdna_df["sequence"] = cdna_df["sequence"].apply(list)
    nts_list = ['A', 'T', 'C', 'G']
    nt_keys = []
    for nt in nts_list:
        replacer_dict = {k: 0 for k in nts_list}
        replacer_dict[nt] = 1
        nt_key = f"{nt}s_in_seq"
        nt_keys.append(nt_key)
        cdna_df[nt_key] = cdna_df["sequence"].apply(lambda seq: list(map(replacer_dict.get, list(seq), list(seq))))
    cdna_df["seq_index"] = cdna_df["sequence"].apply(lambda seq: list(range(1, len(seq) + 1)))
    explode_column_list = ["sequence", "seq_index"] + nt_keys
    cdna_df = cdna_df.explode(explode_column_list, ignore_index=True)

    df_dict = {}
    for gene in cdna_df["gene_name"].unique():
        df_dict[gene] = cdna_df[cdna_df["gene_name"] == gene]
        for nt in nts_list:
            nt_key = f"{nt}_rolling_window"
            df_dict[gene][nt_key] = df_dict[gene][f"{nt}s_in_seq"].rolling(window_size).mean()
    concat_df = pd.concat(df_dict.values())
    plot_cdna_df(concat_df, window_size)
    return concat_df


def plot_cdna_df(exploded_cdna_df: pd.DataFrame, window_size):
    sea.set(font_scale=4)
    sea.set_style("white")
    sea.set_context("poster")
    plt.rcParams.update({"font.size": 60})
    g = sea.relplot(data=exploded_cdna_df, x="seq_index", y="T_rolling_window",
                    row="gene_name", kind="line", aspect=10,
                    facet_kws=dict(sharey=True, sharex=True),
                    linewidth=10)
    g.fig.suptitle(f"Rolling windows ({window_size}nts) of T fractions for MtDNA mRNAs", y=.99)
    g.fig.subplots_adjust(top=.98)
    g.refline(y=0.5)
    for ax in g.axes:
        ax = ax[0]
        ax.yaxis.set_ticks(np.arange(0.0, 1.0, 0.1))
        gene = ax.axes.get_title().split()[2]
        ax.set_ylabel(f"{gene}", fontsize=60)
        ax.set_xlabel('')
        ax.set_title('')
        for i, label in enumerate(ax.yaxis.get_ticklabels()):
            if i % 2 == 0:
                label.set_visible(False)
            else:
                label.set_visible(True)
                label.set_fontsize(45)
        ax.yaxis.grid(True)
        ax.xaxis.grid(False)
    from nanoporePipelineCommon import get_dt
    plt.savefig(f'./testOutputs/{get_dt(for_file=True)}_MtRNA_T-frac_{window_size}ntWindow.png')
    plt.show()


def cdna_df_to_t_stretches(cdna_df: pd.DataFrame, min_len: int = 4):
    import re
    from collections import Counter
    t_stretches_df = pd.DataFrame()
    for row_dict in cdna_df.to_dict(orient="records"):
        seq = row_dict["sequence"]
        matches = [(m.group(), m.start()) for m in re.finditer(r'([T])\1{3,}', seq) if len(m.group()) >= min_len]
        # print(matches)
        t_stretch_lengths = [len(ts) for (ts, loc) in matches]
        t_stretch_len_counts = Counter(t_stretch_lengths)
        t_stretch_len_counts = {f"{k}T_stretches": v*k for k, v in t_stretch_len_counts.items()}
        # print(t_stretch_lengths)
        append_dict = {**row_dict, **t_stretch_len_counts}
        t_stretches_df = t_stretches_df.append(append_dict, ignore_index=True)
    t_stretches_df.fillna(0.0, inplace=True)
    
    t_stretch_cols_and_dtypes = {#"4T_stretches": "int32",
                                 "5T_stretches": "int32",
                                 "6T_stretches": "int32",
                                 "7T_stretches": "int32",
                                 "8T_stretches": "int32",
                                 }
    t_stretch_cols = list(t_stretch_cols_and_dtypes.keys())
    t_stretches_df = t_stretches_df.astype(t_stretch_cols_and_dtypes)
    t_stretch_df = t_stretches_df[["gene_id", "gene_name"] + t_stretch_cols]
    t_stretch_df.drop(columns="gene_id").plot(x="gene_name", kind="bar")
    plt.tight_layout()
    plt.show()
    print(t_stretches_df)


if __name__ == '__main__':
    path_to_cdna = "/data16/marcus/genomes/elegansRelease100/Caenorhabditis_elegans.WBcel235.cdna.all.fa"
    df = load_cdna_lib(path_to_cdna)
    # for size in [10, 25, 50, 100]:
    #     explode_cdna_df_and_array_convert(df, window_size=size)
    cdna_df_to_t_stretches(df, min_len=5)
