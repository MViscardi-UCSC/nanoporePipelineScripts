"""
nanoporePipelineCommon.py
Marcus Viscardi,    September 22, 2021

A common location for some often used methods.
"""

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
    import pandas as pd
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
    import pandas as pd
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


if __name__ == '__main__':
    # parquet = tsv_to_parquet("./Caenorhabditis_elegans.WBcel235.100.gtf.dataframe_parse.tsv")
    # print(pd.read_parquet(parquet).info())
    # print(load_ski_pelo_targets(as_df=True))
    print(rev_compliment("AAACCCGGGTTT"))
