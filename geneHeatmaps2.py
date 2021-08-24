"""
geneHeatmaps2.py
Marcus Viscardi,    August 14, 2021

So I managed to overcomplicate everything while trying
    different methods in geneHeatmaps.py, so now I am
    going to try and take a simpler approach.

The plan is to import the mergedOnReads.tsv rather than
    any of the compressedOnGenes files. B/c the majority
    of the early processing steps to get distances from
    stops are on a per read basis: I'll use the read/row
    format of the mergedOnReads file to be able to
    iterate through things more cleanly/obviously.

CIGAR parsing will still be necessary to "flip" negative
    strand reads to get their 5' ends.

The other introduced change will be to entirely use the
    Caenorhabditis_elegans.WBcel235.100.allChrs.txt file
    which is produced by Josh's prepareReadAssignmentFile3.py
This annotation file is in the format:
    CHR_chrpos  WBGene:strand   iso1:posRelStart:posRelStop|iso2:...
The only striking issue is that Josh's dropped any locations with
    overlapping genes (which makes heaps of sense for short read
    sequencing, but may weaken my nanopore stuff).

After getting stop distance for each read, I'll be able to compress
    these on genes OR transcripts, calculate CDFs and feed them into
    the same plot function from geneHeatmaps.py (or I can rewrite
    that bit too. . . ).

Hype!
"""

import pandas as pd

from step0_nanopore_pipeline import find_newest_matching_file


def load_reads(read_tsv) -> pd.DataFrame:
    print(f"Loading merged on reads file from: {read_tsv} ", end="")
    df = pd.read_csv(read_tsv, sep="\t")
    print(". ")
    return df


def load_read_assignment_tsv(assignment_file_tsv, save_file=True) -> pd.DataFrame:
    print(f"Loading read assignment file from: {assignment_file_tsv} ", end="")
    df = pd.read_csv(assignment_file_tsv, sep="\t",
                     names=["chr_pos", "gene_id_strand", "transcript_ids"])
    print(". ", end="")
    df[["chr_id", "chr_pos"]] = df["chr_pos"].str.split("_", expand=True)
    print(". ", end="")
    df[["gene_id", "strand"]] = df["gene_id_strand"].str.split(":", expand=True)
    print(". ", end="")
    df.drop(columns="gene_id_strand", inplace=True)
    print(". ", end="")
    df["transcript_id"] = df["transcript_ids"].str.split("|")
    print(". ", end="")
    df = df.explode("transcript_id", ignore_index=True)
    print(". ", end="")
    df.drop(columns="transcript_ids", inplace=True)
    print(". ", end="")
    df[["transcript_id", "to_start", "to_stop"]] = df["transcript_id"].str.split(":", expand=True)
    print(". ", end="")
    df = df.astype({"chr_pos": "int32",
                    "to_start": "int32",
                    "to_stop": "int32",
                    "strand": "category",
                    "chr_id": "category",
                    })
    print(". ", end="")
    df = df[["chr_id", "chr_pos", "gene_id", "strand", "transcript_id", "to_start", "to_stop"]]
    print(". ", end="")
    if save_file:
        parquet_path = assignment_file_tsv.rstrip(".txt") + ".parquet"
        df.to_parquet(parquet_path)
        print("(saved parquet) ", end="")
    print(". ")
    return df


def load_read_assignment_parquet(assignment_file_parquet, save_file=True) -> pd.DataFrame:
    print(f"Loading read assignment file from: {assignment_file_parquet} ", end="")
    df = pd.read_parquet(assignment_file_parquet)
    print(". ")
    return df


def merge_on_chr_pos(read_assignment_df: pd.DataFrame, reads_df: pd.DataFrame,
                     subsample=None) -> pd.DataFrame:
    if isinstance(subsample, int):
        read_assignment_df = read_assignment_df.sample(subsample)
        reads_df = reads_df.sample(subsample)
    merge_df = reads_df.merge(read_assignment_df, on=["chr_id", "chr_pos"])
    return merge_df


def main(library_str, genome_dir="/data16/marcus/genomes/elegansRelease100"):
    working_dir_dict = {"polyA": "210528_NanoporeRun_0639_L3s",  # Best (overkill) depth
                        "riboD": "210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3",  # Gross
                        "totalRNA": "210709_NanoporeRun_totalRNA_0639_L3",  # Low depth
                        "polyA2": "210719_nanoporeRun_polyA_0639_L3_replicate",  # Good depth
                        "totalRNA2": "210720_nanoporeRun_totalRNA_0639_L3_replicate",  # Good depth
                        }
    try:
        working_dir = f"/data16/marcus/working/{working_dir_dict[library_str]}"
    except KeyError:
        raise KeyError(f"Key: {library_str} isn't in the working directory dictionary keys : {working_dir_dict.keys()}")
    
    reads_df = load_reads(find_newest_matching_file(f"{working_dir}/output_dir/merge_files/*_mergedOnReads.tsv"))
    read_assignment_df = load_read_assignment_parquet(f"{genome_dir}/"
                                                      f"Caenorhabditis_elegans.WBcel235.100.allChrs.parquet")
    print("Finished loading files!")
    for df in [reads_df, read_assignment_df]:
        print(df.info())
    merged_df = merge_on_chr_pos(read_assignment_df, reads_df, subsample=500000)
    breakpoint()


if __name__ == '__main__':
    main("totalRNA2")
