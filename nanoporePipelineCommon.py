"""
nanoporePipelineCommon.py
Marcus Viscardi,    September 22, 2021

A common location for some often used methods.
"""
import pandas as pd


def get_dt(for_print=False, for_file=False):
    from datetime import datetime
    now = datetime.now()
    if for_print:
        return str(now.strftime("%m/%d/%y @ %I:%M:%S %p"))
    elif for_file:
        return str(now.strftime("%y%m%d"))
    else:
        return str(now.strftime("%y%m%d_%I:%M:%S%p"))


def tsv_to_parquet(tsv_path) -> str:
    parquet_path = tsv_path.rstrip("tsv") + "parquet"
    df = pd.read_csv(tsv_path, sep="\t")
    print(f"Saving new parquet to: {parquet_path}")
    df.to_parquet(parquet_path)
    return parquet_path


def find_newest_matching_file(path_str):
    # A tool that I'll want to use to grab the most recent file
    from os import path
    from glob import glob

    list_of_files = glob(path_str)
    try:
        latest_file = max(list_of_files, key=path.getctime)
        return latest_file
    except ValueError:
        raise ValueError(f"Failed to find any files matching \"{path_str}\"")


def load_ski_pelo_targets(as_df=False):
    df = pd.read_csv("/data16/marcus/working/210119_SkiPeloTargets_fromStarDust/"
                     "170723_MSandM.wtAndSkiPelo_Bounds_-12_-14_S.DESeqgeneCts_"
                     "diffExpression_2.7319418642771283e-06Down.txt", names=["gene_id"])
    if as_df:
        return df
    else:
        return df.gene_id.to_list()


def rev_compliment(seq: str, rna: bool = False) -> str:
    from Bio.Seq import Seq
    seq = Seq(seq)
    if rna:
        return seq.reverse_complement_rna()
    else:
        return seq.reverse_complement()


# Below class and methods are for loading BAM files directly to pandas dataframes
from typing import NamedTuple


class BamHeadersAndDf(NamedTuple):
    """
    This is a dataclass object to hold the a sam DF and
    the correct sam/bam headers
    """
    headers: str
    df: pd.DataFrame


def minimap_bam_to_df(bam_path, drop_secondaries_and_unmapped=True,
                      name_columns=True) -> BamHeadersAndDf:
    from subprocess import check_output
    from io import BytesIO

    if drop_secondaries_and_unmapped:
        drop_flag = '-F 0x904 '
    else:
        drop_flag = ''

    # First read the bam file into a tab-seperated string object:
    output = check_output(f"samtools view {drop_flag}{bam_path}", shell=True)

    # Use pandas to load this string object into a dataframe
    df = pd.read_csv(BytesIO(output),
                     encoding='utf8',
                     sep="\t",
                     names=range(22),
                     low_memory=False
                     )

    if name_columns:
        # Column names will make handling the dataframe easier,
        #   but they are not going to end up in the new sam/bam
        minimap_bam_header_names = ["read_id",
                                    "bit_flag",
                                    "chr_id",
                                    "chr_pos",
                                    "mapq",
                                    "cigar",
                                    "r_next",
                                    "p_next",
                                    "len",
                                    "sequence",
                                    "phred_qual",
                                    "num_mismatches",
                                    "best_dp_score",
                                    "dp_score",
                                    "num_ambiguous_bases",
                                    "transcript_strand",
                                    "type_of_alignment",
                                    "num_minimizers",
                                    18,
                                    19,
                                    20,
                                    21]
        df = df.rename(columns=dict(enumerate(minimap_bam_header_names)))

    header = check_output(f"samtools view -H {bam_path}", shell=True).decode("utf-8")
    output = BamHeadersAndDf(header, df)
    return output


def save_sorted_bam_obj(bam_obj: BamHeadersAndDf, output_path: str,
                        index: bool = False) -> None:
    from subprocess import run
    header, df = bam_obj
    buffer = header + df.to_csv(sep="\t",
                                header=False,
                                index=False)
    # subprocess.run accepts the input param to pass to the bash call!
    run(f"samtools view -S -b - | samtools sort -o {output_path}.sorted.bam",
        input=buffer.encode('utf-8'), shell=True)
    if index:
        run(f'samtools index {output_path}.sorted.bam', shell=True)


if __name__ == '__main__':
    # parquet = tsv_to_parquet("./Caenorhabditis_elegans.WBcel235.100.gtf.dataframe_parse.tsv")
    # print(pd.read_parquet(parquet).info())
    # print(load_ski_pelo_targets(as_df=True))
    print(rev_compliment("AAACCCGGGTTT"))
