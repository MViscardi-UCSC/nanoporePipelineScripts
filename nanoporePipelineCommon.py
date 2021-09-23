"""
nanoporePipelineCommon.py
Marcus Viscardi,    September 22, 2021

A common location for some often used methods.
"""
import pandas as pd


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


if __name__ == '__main__':
    parquet = tsv_to_parquet("./Caenorhabditis_elegans.WBcel235.100.gtf.dataframe_parse.tsv")
    print(pd.read_parquet(parquet).info())
