"""
nanoporePipelineCommon.py
Marcus Viscardi,    September 22, 2021

A common location for some often used methods.
"""
import pandas as pd
import os

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)


class FastqFile:
    def __init__(self, path, head=None):
        from tqdm import tqdm
        import mappy as mp
        self.path = path
        self.subset = isinstance(head, int)
        fastq_items_list = []  # List to hold fastq items as we iterate over
        print(f"\nStarting fastq iterative load @ {get_dt(for_print=True)}", end="\n")
        # Below will allow us to track how long the fastq parse is taking:
        row_iterator = tqdm(mp.fastx_read(path, read_comment=True),
                            total=sum(1 for line in open(path)) // 4)
        for line, (read_id, sequence, quality, comment) in enumerate(row_iterator):
            fastq_items_list.append([read_id, sequence, "+", quality, comment])
            row_iterator.set_description(f"Processing {read_id}")
            if self.subset and line >= head:
                break
        # Convert the fastq items list into a pandas dataframe so it can be filtered by the alt_mapped_reads_df
        self.df = pd.DataFrame(fastq_items_list, columns=["read_id", "sequence", "plus", "quality", "comment"])

    def __str__(self):
        return str(self.df)

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


class SamOrBamFile:
    def __init__(self, path, header_source=None, subsample=None):
        from subprocess import check_output

        self.path = path
        if not os.path.exists(self.path):
            raise FileNotFoundError(f"Could not find file: {self.path}")
        self.max_cols = 30  # Probably shouldn't need to be changeable...

        self.file_type = self.path.split('.')[-1]
        if self.file_type not in ['bam', 'sam']:
            raise NotImplementedError(f"Please provide a file that ends with bam/sam, "
                                      f"your file ends in: {self.file_type}")
        self.sam_main_columns = ["read_id",
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
        # Anything not included below will either stay as a string/object,
        #   or whatever the simplesam decided to store it as (for tags):
        self.dtype_dict = {'bit_flag': 'uint16',
                           'chr_id': 'category',
                           'chr_pos': 'uint32',
                           'mapq': 'uint8',
                           'NM': 'uint16',
                           'ms': 'uint32',
                           'AS': 'uint32',
                           'nn': 'uint32',
                           'ts': 'category',
                           'tp': 'category',
                           'cm': 'uint32',
                           'rl': 'uint32',
                           't5': 'category',
                           't3': 'category'}
        self.tag_type_dict = {'NM': 'i',
                              'ms': 'i',
                              'AS': 'i',
                              'nn': 'i',
                              'ts': 'A',
                              'tp': 'A',
                              'cm': 'i',
                              's1': 'i',
                              's2': 'i',
                              'de': 'f',
                              'rl': 'i',
                              't5': 'A',
                              't3': 'A',
                              'SA': 'Z',
                              }

        if isinstance(header_source, str):
            self.header = check_output(f"samtools view -H {header_source}", shell=True).decode("utf-8")
            read_file_header = check_output(f"samtools view --no-PG -H {path}", shell=True).decode("utf-8")
            self.header_lines = len(read_file_header.split('\n')) - 1
        else:
            self.header = check_output(f"samtools view --no-PG -H {path}", shell=True).decode("utf-8")
            self.header_lines = len(self.header.split('\n')) - 1

        if isinstance(subsample, int):
            self.subsample = subsample
        else:
            self.subsample = None

        self.df = self.__build_df__()

    # def __build_df__(self):
    #     from subprocess import check_output
    #     from io import BytesIO
    #     initial_read_cols = list(range(self.max_cols))
    # 
    #     # These numbered headers are the variably placed tags! Parse 'em later!
    #     sam_header_names = self.sam_main_columns + list(range(11, self.max_cols))
    # 
    #     print(f"{get_dt(for_print=True)}: Starting to load sam file from: {self.path}")
    #     if self.path.endswith('.bam'):
    #         # First read the bam file into a tab-seperated string object:
    #         output = check_output(f"samtools view {self.path}", shell=True)
    # 
    #         # Use pandas to load this string object into a dataframe
    #         temp_df = pd.read_csv(BytesIO(output),
    #                               encoding='utf8',
    #                               sep="\t", names=initial_read_cols,
    #                               low_memory=False, index_col=False,
    #                               nrows=self.subsample,
    #                               quotechar='\0',
    #                               )
    #     else:
    #         temp_df = pd.read_csv(self.path, sep="\t", names=initial_read_cols,
    #                               low_memory=False, index_col=False,
    #                               skiprows=self.header_lines,
    #                               nrows=self.subsample,
    #                               quotechar='\0',
    #                               )
    #     print(f"SAM file loaded into a dataframe, starting tag parsing at {get_dt(for_print=True)}")
    #     temp_df = temp_df.rename(columns=dict(enumerate(sam_header_names)))
    #     temp_df = temp_df.fillna("*")
    #     temp_df['tags_list'] = temp_df[list(range(11, self.max_cols))].values.tolist()
    #     tags_df = temp_df['tags_list'].apply(lambda row: dict([(f"sam_tag|{tag_items[0]}:{tag_items[1]}", tag_items[2])
    #                                                            for tag_items in
    #                                                            [tag.split(':') for tag in row if tag != '*']])
    #                                          ).apply(pd.Series)
    #     final_df = pd.concat([temp_df.drop(['tags_list'], axis=1), tags_df],
    #                          axis=1).drop(list(range(11, self.max_cols)), axis=1)
    #     columns = list(final_df.columns)
    #     new_columns = columns[:-2] + [columns[-1]] + [columns[-2]]
    #     final_df = final_df.reindex(columns=new_columns)
    #     print(f"Completed tag parsing at {get_dt(for_print=True)}")
    #     return final_df

    def __build_df__(self):
        import simplesam as ssam
        from tqdm import tqdm
        from collections import ChainMap

        print(f"{get_dt(for_print=True)}: Starting to load {self.file_type} file from: {self.path}")
        sam_list_of_dicts = []
        with ssam.Reader(open(self.path, 'r')) as in_sam:
            to_subsample = isinstance(self.subsample, int)
            if to_subsample:
                max_to_load = self.subsample
            else:
                if self.file_type == 'bam':
                    max_to_load = len(in_sam)
                else:
                    max_to_load = sum(1 for line in open(self.path))
            sam_iterator = tqdm(enumerate(in_sam), total=max_to_load)
            for i, read in sam_iterator:
                sam_iterator.set_description(f"Processing {read.qname}")
                main_info_list = read.__str__().split('\t')[0:11]
                main_info = dict(zip(self.sam_main_columns, main_info_list))
                tag_info = read.tags
                all_info = dict(ChainMap(main_info, tag_info))
                sam_list_of_dicts.append(all_info)
                if to_subsample and i > self.subsample:
                    break
        final_df = pd.DataFrame(sam_list_of_dicts)
        self.tag_columns = [col for col in final_df if col not in self.sam_main_columns]
        final_df = final_df[self.sam_main_columns + self.tag_columns]
        
        # Filter the datatype dictionary, so we don't get errors from
        #   trying to retype columns that don't exist!
        self.dtype_dict = {k: v for k, v in self.dtype_dict.items() if k in list(final_df.columns)}
        # Adjust any datatypes we want:
        final_df = final_df.astype(self.dtype_dict)
        print(f"{get_dt(for_print=True)}: Finished loading {self.file_type} file!")
        return final_df

    def to_sam(self, output_path=None, escape_char="|", to_bam=False,
               infer_tag_types=False):
        import simplesam as ssam
        from subprocess import run
        from csv import QUOTE_NONE
        if escape_char not in "~}|":
            raise NotImplementedError("escape_char has to be |, }, or ~. Everything else is in the Phred Quals!")
        save_df = self.df.copy()
        if infer_tag_types:
            print(f"Inferring tag datatypes is sorta broken... Use at your own risk!")
            for col in self.tag_columns:
                save_df[col] = save_df[col].apply(lambda x: ssam.encode_tag(col, x))
        else:
            for col in self.tag_columns:
                try:
                    save_df[col] = col + ':' + self.tag_type_dict[col] + ':' + save_df[col].astype(str)
                except KeyError:
                    print(f"{col} not found in tag_type_dict and caused a KeyError??")
                    save_df[col] = col + ':Z:' + save_df[col].map(str)
        save_df["tags"] = save_df[self.tag_columns].values.tolist()
        save_df["tags"] = save_df["tags"].apply(lambda row: '\t'.join([x for x in row if not str(x).endswith('nan')]))
        save_df = save_df.drop(self.tag_columns, axis=1)

        temp_header = self.header + f"@CO\t{get_dt(for_print=True)}: Introduced edits with python . . . " \
                                    f"More info hopefully added to this later!\n"
        buffer = temp_header + save_df.to_csv(sep="\t",
                                              header=False,
                                              index=False,
                                              quoting=QUOTE_NONE,
                                              escapechar=escape_char,
                                              quotechar='\0').replace(escape_char, '')
        # subprocess.run accepts the input param to pass to the bash call!
        if output_path:
            if to_bam:
                run(f"samtools view -h --no-PG -b - > {output_path}",
                    input=buffer.encode('utf-8'), shell=True)
            else:
                run(f"samtools view -h --no-PG - > {output_path}",
                    input=buffer.encode('utf-8'), shell=True)
        else:
            print(save_df)
        return save_df

    # def to_sam(self, output_path, escape_char="|", to_bam=False):
    #     from subprocess import run
    #     from csv import QUOTE_NONE
    #     if escape_char not in "~}|":
    #         raise NotImplementedError("escape_char has to be |, }, or ~. Everything else is in the Phred Quals!")
    #     save_df = self.df.copy()
    #     sam_tag_cols = [col for col in list(save_df.columns) if col.startswith("sam_tag|")]
    #     for column in sam_tag_cols:
    #         save_df[column] = column.split("|")[1] + ":" + save_df[column]
    #     save_df["tags"] = save_df[sam_tag_cols].values.tolist()
    #     save_df["tags"] = save_df["tags"].apply(lambda row: '\t'.join([x for x in row if str(x) != 'nan']))
    #     save_df = save_df.drop(sam_tag_cols, axis=1)
    #     temp_header = self.header + f"@CO\t{get_dt(for_print=True)}: Introduced edits with python . . . " \
    #                                 f"More info hopefully added to this later!\n"
    #     buffer = temp_header + save_df.to_csv(sep="\t",
    #                                           header=False,
    #                                           index=False,
    #                                           quoting=QUOTE_NONE,
    #                                           escapechar=escape_char,
    #                                           quotechar='\0').replace(escape_char, '')
    #     # subprocess.run accepts the input param to pass to the bash call!
    #     if to_bam:
    #         run(f"samtools view -h --no-PG -b - > {output_path}",
    #             input=buffer.encode('utf-8'), shell=True)
    #     else:
    #         run(f"samtools view -h --no-PG - > {output_path}",
    #             input=buffer.encode('utf-8'), shell=True)

    def to_bam(self, output_path):
        self.to_sam(output_path=output_path, to_bam=True)


def pick_libs_return_paths_dict(lib_list: list, file_suffix: str = "parquet", file_midfix="mergedOnReads",
                                output_dir_folder="merge_files", return_all: bool = False) -> dict:
    output_dir_dict = {
        "riboD": "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/output_dir",
        "totalRNA": "/data16/marcus/working/210709_NanoporeRun_totalRNA_0639_L3/"
                    "output_dir",
        "totalRNA2": "/data16/marcus/working/"
                     "210720_nanoporeRun_totalRNA_0639_L3_replicate/output_dir",
        "polyA": "/data16/marcus/working/210528_NanoporeRun_0639_L3s/output_dir",
        "polyA2": "/data16/marcus/working/210719_nanoporeRun_polyA_0639_L3_replicate/output_dir",
        "xrn-1": "/data16/marcus/working/210905_nanoporeRun_totalRNA_5108_xrn-1-KD/output_dir",
        "xrn-1-5tera": "/data16/marcus/working/211118_nanoporeRun_totalRNA_5108_xrn-1-KD_5TERA/output_dir",
        "pTRI-stds": "/data16/marcus/working/211121_nanoporeRun_pTRIstds/output_dir",
    }
    if return_all:
        lib_list = output_dir_dict.keys()
    file_suffix = file_suffix.strip(".")
    return_dict = {}
    for lib_key, output_dir in output_dir_dict.items():
        if lib_key in lib_list:
            file_path = f"{output_dir}/{output_dir_folder}/*_{file_midfix}.{file_suffix}"
            print(f"Looking for file for {lib_key}, at {file_path}...", end=" ")
            return_dict[lib_key] = find_newest_matching_file(file_path)
            print(f"File Found.")
    if not return_dict:
        raise KeyError(f"No matching library keys found, please use keys from the following list: "
                       f"{list(output_dir_dict.keys())}. "
                       f"You passed the following keys: "
                       f"{lib_list}")
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
    test_sam_path = "./testInputs/pTRI_test.sorted.bam"
    # test_sam_path = "/data16/marcus/working/211101_nanoporeSoftLinks/" \
    #                 "211118_nanoporeRun_totalRNA_5108_xrn-1-KD_5TERA/" \
    #                 "output_dir/cat_files/cat.sorted.mappedAndPrimary.sam"
    sam = SamOrBamFile(test_sam_path,
                       # subsample=50000,
                       )
    # sam.to_sam(test_sam_path + '.resave.sam', escape_char="~")
    sam.to_sam(output_path=test_sam_path.rstrip('.sam') + '.resave.sam')
    print(sam)
    # fastq = FastqFile("./testOutputs/in.fastq")
    # fastq.save_to_fastq("./testOutputs/out.fastq")
    # print(fastq)
