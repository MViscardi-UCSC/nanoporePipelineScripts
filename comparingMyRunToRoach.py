"""
comparingMyRunToRoach.py
Marcus Viscardi     June 2, 2021

Based on the read counts it seems like I have a MUCH greater (~5X) read depth! To see if this is just an
    effect of my reads being shorter but having similar nucleotide depth leading to a skewed view, I want
    to quantify how many called and mapped bases each data set had.

Going to jump straing to the sorted sam files as these will be the most simple way to do this.
"""

import pandas as pd

pd.set_option("display.max_columns", None)


def sam_to_df(path_to_sam) -> pd.DataFrame:
    print(f"Starting import from file @ {path_to_sam}")
    sam_df = pd.read_csv(path_to_sam,
                         sep="\t", names=range(23))
    # Drop any read_ids that have come up multiple times:
    #   TODO: I am assuming these are multiply mapping? Is this wrong?
    sam_df = sam_df.drop_duplicates(subset=0, keep=False, ignore_index=True)
    # I have no idea what the additional tags from Minimap2 are, I'm dropping them for now:
    sam_df = sam_df[range(11)]
    # And lets rename columns while we are at it!
    sam_header_names = ["read_id",
                        "bit_flag",
                        "chr_id",
                        "chr_pos",
                        "mapq",
                        "cigar",
                        "r_next",
                        "p_next",
                        "len",
                        "sequence",
                        "phred_qual"]
    sam_df = sam_df.rename(columns=dict(enumerate(sam_header_names)))
    return sam_df


def count_nts_mapped(df_from_sam: pd.DataFrame):
    print(f"\tCalculating read lengths for {df_from_sam.shape[0]} lines of .sam file")
    df_from_sam["read_length"] = df_from_sam["sequence"].str.len()
    print("\tCalculating sum, mean, median, max, and min")
    nt_sum = df_from_sam["read_length"].sum()
    mean = df_from_sam["read_length"].mean()
    median = df_from_sam["read_length"].median()
    min = df_from_sam["read_length"].min()
    max = df_from_sam["read_length"].max()
    return nt_sum, mean, median, min, max, df_from_sam.shape[0]


if __name__ == '__main__':
    path_to_roach = "/data16/marcus/working/210106_nanopolish_wRoachData/output_dir_youngadultrep1tech1/cat_allReads.sam"
    path_to_my = "/data16/marcus/working/210528_NanoporeRun_0639_L3s/output_dir/cat_files/cat.sorted.sam"
    dictionary = {"Roach data (young adult worms - replicate 1)": path_to_roach,
                  "My data (L3 worms)": path_to_my}
    results = []
    for name, path in dictionary.items():
        print(f"Working on {name}:")
        df = sam_to_df(path)
        results.append(count_nts_mapped(df))
    returns = results
    print(f"         Roach:\t Viscardi:\n"
          f"Reads:\t{returns[0][5]:10.0f}\t{returns[1][5]:10.0f}\n"
          f"Sum:  \t{returns[0][0]:10.0f}\t{returns[1][0]:10.0f}\n"
          f"Mean: \t{returns[0][1]:10.0f}\t{returns[1][1]:10.0f}\n"
          f"Median:\t{returns[0][2]:10.0f}\t{returns[1][2]:10.0f}\n"
          f"Min:  \t{returns[0][3]:10.0f}\t{returns[1][3]:10.0f}\n"
          f"Max:  \t{returns[0][4]:10.0f}\t{returns[1][4]:10.0f}")
