"""
parquetViewer.py
Marcus Viscardi,    November 05, 2021

Quick script to load parquet files into pandasgui
"""
import pandas as pd
pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)
import pandasgui as pg
import sys

if __name__ == '__main__':
    try:
        parquet_path = sys.argv[1]
    except IndexError:
        print(f"Please provide a .parquet file as the first arg")
        sys.exit()
    if not parquet_path.endswith(".parquet"):
        print(f"Please provide a file that ends with '.parquet'")
        sys.exit()
    first_print = f"Starting to load file from: {parquet_path}"
    print(first_print, end="")
    parquet_df = pd.read_parquet(parquet_path)
    print('\b' * len(first_print), "Finished loading file!", sep="")
    pg.show(parquet_df)
