"""
nanoporePipelineCommon.py
Marcus Viscardi,    September 22, 2021

A common location for some often used methods.
"""
import pandas as pd


class FastqFile:
    def __init__(self, path, head=None):
        from tqdm import tqdm
        import mappy as mp
        self.path = path
        self.subset = isinstance(head, int)
        fastq_items_list = []  # List to hold fastq items as we iterate over
        # Below will allow us to track how long the fastq parse is taking:
        row_iterator = tqdm(mp.fastx_read(path, read_comment=True),
                            total=sum(1 for line in open(path)) // 4)
        print(f"Starting fastq iterative load @ {get_dt(for_print=True)}")
        for line, (read_id, sequence, quality, comment) in enumerate(row_iterator):
            fastq_items_list.append([read_id, sequence, "+", quality, comment])
            row_iterator.set_description(f"Processing {read_id}")
            if self.subset and line >= head:
                break
        # Convert the fastq items list into a pandas dataframe so it can be filtered by the alt_mapped_reads_df
        self.df = pd.DataFrame(fastq_items_list, columns=["read_id", "sequence", "plus", "quality", "comment"])

    def filter_against(self, df_w_read_id):
        # Cool way to only keep values that don't appear in alt_mapped_read_df:
        #   (From: https://tinyurl.com/22czvzua)
        # Also trying out query operator, some notes on this here:
        #   https://stackoverflow.com/questions/67341369/pandas-why-query-instead-of-bracket-operator
        self.df = pd.merge(self.df, df_w_read_id,
                           on="read_id",
                           indicator=True,
                           how="outer").query('_merge=="left_only"').drop('_merge', axis=1)

    def save_to_fastq(self, output_path):
        from csv import QUOTE_NONE
        self.df["read_id"] = "@" + self.df["read_id"] + " " + \
                             self.df["comment"]
        self.df.drop("comment", axis=1).to_csv(output_path,
                                               index=False,
                                               header=False,
                                               sep="\n",
                                               quoting=QUOTE_NONE)


def pick_libs_return_paths_dict(lib_list: list, file_suffix: str = "parquet"):
    path_dict = {
        "riboD": "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/output_dir/"
                 "merge_files/*_mergedOnReads.",
        "totalRNA": "/data16/marcus/working/210709_NanoporeRun_totalRNA_0639_L3/"
                    "output_dir/merge_files/*_mergedOnReads.",
        "totalRNA2": "/data16/marcus/working/"
                     "210720_nanoporeRun_totalRNA_0639_L3_replicate/output_dir/"
                     "merge_files/*_mergedOnReads.",
        "polyA": "/data16/marcus/working/210528_NanoporeRun_0639_L3s/output_dir/"
                 "merge_files/*_mergedOnReads.",
        "polyA2": "/data16/marcus/working/210719_nanoporeRun_polyA_0639_L3_replicate/"
                  "output_dir/merge_files/*_mergedOnReads.",
        "xrn-1": "/data16/marcus/working/210905_nanoporeRun_totalRNA_5108_xrn-1-KD/"
                 "output_dir/merge_files/*_mergedOnReads."
    }
    return_dict = {}
    for lib_key, file_prefix_path in path_dict.items():
        if lib_key in lib_list:
            return_dict[lib_key] = find_newest_matching_file(file_prefix_path + file_suffix)
    return return_dict


def load_read_assignments(assignment_file_parquet_path) -> pd.DataFrame:
    print(f"Loading read assignment file from: {assignment_file_parquet_path} ", end="")
    read_assignment_df = pd.read_parquet(assignment_file_parquet_path)
    print(". ")
    return read_assignment_df


def assign_w_josh_method(reads_df, genomeDir):
    
    def merge_on_chr_pos(read_assignment_df: pd.DataFrame, reads_df: pd.DataFrame) -> pd.DataFrame:
        print(f"Merging read assignments and reads at {get_dt(for_print=True)}")
        merge_df = reads_df.merge(read_assignment_df, on=["chr_id", "chr_pos"],
                                  how="left", suffixes=("_fromReads",
                                                        "_fromAssign"))
        # merge_df = merge_df[~(merge_df.gene_id_fromReads.isna() & merge_df.gene_id_fromAssign.isna())]
        # below call drops reads that don't get assigned by Josh's tool
        merge_df = merge_df[~merge_df.gene_id_fromAssign.isna()]
        merge_df = merge_df[merge_df.strand_fromReads == merge_df.strand_fromAssign]
        print(f"Done merging at {get_dt(for_print=True)}")
        return merge_df
    
    read_assignments_df = load_read_assignments(f"{genomeDir}/Caenorhabditis_elegans.WBcel235.100.allChrs.parquet")
    print("Finished loading files!")
    # for df in [reads_df, read_assignments_df]:
    #     print(df.info())
    merged_df = merge_on_chr_pos(read_assignments_df, reads_df)
    return merged_df


def gene_names_to_gene_ids(tsv_path: str = "/data16/marcus/genomes/elegansRelease100"
                                           "/Caenorhabditis_elegans.WBcel235.100.gtf"
                                           ".tsv") -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t")[["gene_name", "gene_id"]].drop_duplicates(ignore_index=True)
    return df


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
    fastq = FastqFile("./testOutputs/in.fastq")
    fastq.save_to_fastq("./testOutputs/out.fastq")
    print(fastq)
