"""
nanoporePipelineCommon.py
Marcus Viscardi,    September 22, 2021

A common location for some often used methods.
"""

import pandas as pd
import numpy as np
from tqdm import tqdm
import os
from subprocess import Popen, CalledProcessError, PIPE
from pathlib import Path
import datetime

import seaborn as sea
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

pd.set_option('display.width', 100)
pd.set_option('display.max_columns', None)

OUTPUT_DIR_DICT = {
    "riboD": "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/"
             "output_dir",
    "totalRNA": "/data16/marcus/working/210709_NanoporeRun_totalRNA_0639_L3/"
                "output_dir",
    "totalRNA1": "/data16/marcus/working/210709_NanoporeRun_totalRNA_0639_L3/"
                 "output_dir",
    "totalRNA2": "/data16/marcus/working/"
                 "210720_nanoporeRun_totalRNA_0639_L3_replicate/"
                 "output_dir",
    "polyA": "/data16/marcus/working/210528_NanoporeRun_0639_L3s/"
             "output_dir",
    "polyA1": "/data16/marcus/working/210528_NanoporeRun_0639_L3s/"
              "output_dir",
    "polyA2": "/data16/marcus/working/210719_nanoporeRun_polyA_0639_L3_replicate/"
              "output_dir",
    "xrn-1": "/data16/marcus/working/210905_nanoporeRun_totalRNA_5108_xrn-1-KD/"
             "output_dir",
    "xrn-1-5tera": "/data16/marcus/working/211118_nanoporeRun_totalRNA_5108_xrn-1-KD_5TERA/"
                   "output_dir",
    "pTRI-stds": "/data16/marcus/working/211121_nanoporeRun_pTRIstds/"
                 "output_dir",
    "xrn-1-5tera-smg-6": "/data16/marcus/working/211210_nanoporeRun_totalRNA_2102_xrn-1-KD_5TERA/"
                         "output_dir",
    "pTRI-stds-tera3": "/data16/marcus/working/211212_nanoporeRun_pTRIstds_TERA3/"
                       "output_dir",
    "polyA3": "/data16/marcus/working/220131_nanoporeRun_polyA_0639_L3_third/"
              "output_dir",
    "totalRNA3": "/data16/marcus/working/220131_nanoporeRun_totalRNA_0639_L3_third/"
                 "output_dir",
    "roach_L3_1": "/data16/marcus/working/220222_roach_analysis_revisit/L3_techRep1/"
                  "output_dir",
    "roach_L3_2": "/data16/marcus/working/220222_roach_analysis_revisit/L3_techRep2/"
                  "output_dir",
    "roach_L4_1": "/data16/marcus/working/220222_roach_analysis_revisit/L4_techRep1/"
                  "output_dir",
    "roach_L4_2": "/data16/marcus/working/220222_roach_analysis_revisit/L4_techRep2/"
                  "output_dir",
    "5tera_xrn-1-KD_wt": "/data16/marcus/working/221216_nanoporeRun_totalRNA_wt_xrn-1-KD_5TERA/"
                         "output_dir",
    "5tera_xrn-1-KD_smg-5": "/data16/marcus/working/221216_nanoporeRun_totalRNA_smg-5_xrn-1-KD_5TERA/"
                            "output_dir",
    "5tera_xrn-1-KD_smg-6": "/data16/marcus/working/221216_nanoporeRun_totalRNA_smg-6_xrn-1-KD_5TERA/"
                            "output_dir",
    "5tera_xrn-1-KD_smg-7": "/data16/marcus/working/221216_nanoporeRun_totalRNA_smg-7_xrn-1-KD_5TERA/"
                            "output_dir",
    "sPM58": "/data16/marcus/working/230207_nanoporeRun_totalRNA_sPM58_xrn-1-KD_5PTERA/"
             "output_dir",
    "sPM57": "/data16/marcus/working/230130_nanoporeRun_totalRNA_sPM57_xrn-1-KD_5PTERA/"
             "output_dir",
    "5tera_xrn-1-KD_wt_rerun": "/data16/marcus/working/230327_nanoporeRun_totalRNA_wt_xrn-1-KD_5TERA_rerun/"
                               "output_dir",
    "5tera_xrn-1-KD_smg-6_rerun": "/data16/marcus/working/230403_nanoporeRun_totalRNA_smg-6_xrn-1-KD_5TERA_rerun/"
                                  "output_dir",
    "5tera_xrn-1-KD_smg-5_rerun": "/data16/marcus/working/230410_nanoporeRun_totalRNA_smg-5_xrn-1-KD_5TERA_rerun/"
                                  "output_dir",
}
REV_OUTPUT_DIR_DICT = {v: k for k, v in OUTPUT_DIR_DICT.items()}

CONVERSION_DICT = {"xrn-1-5tera": "oldN2",
                   "xrn-1-5tera-smg-6": "oldS6",
                   "5tera_xrn-1-KD_wt": "newN2",
                   "5tera_xrn-1-KD_smg-5": "newS5",
                   "5tera_xrn-1-KD_smg-6": "newS6",
                   "5tera_xrn-1-KD_smg-7": "newS7",
                   "5tera_xrn-1-KD_wt_rerun": "newerN2",
                   "5tera_xrn-1-KD_smg-6_rerun": "newerS6",
                   "5tera_xrn-1-KD_smg-5_rerun": "newerS5",
                   "sPM57": "sPM57",
                   "sPM58": "sPM58",
                   }
REV_CONVERSION_DICT = {val: key for key, val in CONVERSION_DICT.items()}


class SamOrBamFile:
    def __init__(self, path, header_source=None, subsample=-1):
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
            # I don't remember what the plan was with storing this other header?:
            read_file_header = check_output(f"samtools view --no-PG -H {path}", shell=True).decode("utf-8")
            self.header_lines = len(self.header.split('\n')) - 1
        else:
            self.header = check_output(f"samtools view --no-PG -H {path}", shell=True).decode("utf-8")
            self.header_lines = len(self.header.split('\n')) - 1

        self.subsample = subsample

        self.df = self.__build_df__()

    def __build_df__(self):
        import simplesam as ssam
        from tqdm import tqdm
        from collections import ChainMap

        print(f"{get_dt(for_print=True)}: Starting to load {self.file_type} file from: {self.path}")
        sam_list_of_dicts = []
        with ssam.Reader(open(self.path, 'r')) as in_sam:
            to_subsample = self.subsample > 0
            if to_subsample:
                max_to_load = self.subsample
            else:
                if self.file_type == 'bam':
                    max_to_load = len(in_sam)  # TODO: Change this step so it ignores header lines!!
                else:
                    max_to_load = sum(1 for line in open(self.path) if not line.startswith(
                        "@"))  # TODO: This may cause an issue with nanopore reads who's names start with @!!
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

    def add_column_w_merge(self, df_to_merge_with: pd.DataFrame, cols_to_merge_on: list or tuple):
        self.df = self.df.merge(df_to_merge_with, on=cols_to_merge_on, how='left')
        print(f"Completed merge using columns: {cols_to_merge_on}")

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
                    print(f"For now, we are going to use simplesam to infer the correct type!")
                    save_df[col] = save_df[col].apply(lambda x: ssam.encode_tag(col, x))
                    # save_df[col] = col + ':Z:' + save_df[col].map(str)
        save_df["tags"] = save_df[self.tag_columns].values.tolist()
        save_df["tags"] = save_df["tags"].apply(lambda row: '\t'.join([x for x in row if not str(x).endswith('nan')]))
        save_df = save_df.drop(self.tag_columns, axis=1)

        temp_header = self.header + f"@CO\t{get_dt(for_print=True)}: Introduced edits with python . . . " \
                                    f"More info hopefully added to this later!\n"  # TODO: Add more info...
        buffer = temp_header + save_df.to_csv(sep="\t",
                                              header=False,
                                              index=False,
                                              quoting=QUOTE_NONE,
                                              escapechar=escape_char,
                                              quotechar='\0').replace(escape_char, '')
        # subprocess.run accepts the input param to pass to the bash call!
        if output_path:
            if to_bam:
                run(f"samtools view -h --no-PG -b - > {output_path} && samtools index {output_path}",
                    input=buffer.encode('utf-8'), shell=True)
            else:
                run(f"samtools view -h --no-PG - > {output_path}",
                    input=buffer.encode('utf-8'), shell=True)
        else:
            print(save_df)
        return save_df

    def to_bam(self, output_path):
        self.to_sam(output_path=output_path, to_bam=True)

    def __str__(self):
        return self.df.__str__()


class FastqFile:
    def __init__(self, path, subset_entries=-1):
        from tqdm import tqdm
        import mappy as mp
        self.path = path
        self.subset = isinstance(subset_entries, int) and subset_entries > 0
        fastq_items_list = []  # List to hold fastq items as we iterate over
        print(f"\nStarting fastq iterative load @ {get_dt(for_print=True)}", end="\n")
        # Below will allow us to track how long the fastq parse is taking:
        row_iterator = tqdm(mp.fastx_read(path, read_comment=True),
                            total=sum(1 for line in open(path)) // 4)
        for line, (read_id, sequence, quality, comment) in enumerate(row_iterator):
            fastq_items_list.append([read_id, sequence, "+", quality, comment])
            row_iterator.set_description(f"Processing {read_id}")
            if self.subset and line >= subset_entries:
                break
        # Convert the fastq items list into a pandas dataframe so it can be filtered by the alt_mapped_reads_df
        self.df = pd.DataFrame(fastq_items_list, columns=["read_id", "sequence", "plus", "quality", "comment"])

    def __str__(self):
        return self.save_to_fastq()

    def filter_against(self, filter_df_w_read_ids, overwrite_df=False):
        # Cool way to only keep values that don't appear in filter_df_w_read_ids:
        #   (From: https://tinyurl.com/22czvzua)
        # Also trying out query operator, some notes on this here:
        #   https://stackoverflow.com/questions/67341369/pandas-why-query-instead-of-bracket-operator
        filtered_df = pd.merge(self.df, filter_df_w_read_ids,
                               on="read_id",
                               indicator=True,
                               how="outer").query('_merge=="left_only"').drop('_merge', axis=1)
        if overwrite_df:
            self.df = filtered_df
        return filtered_df

    def filter_for(self, filter_df_w_read_ids, overwrite_df=False):
        # Only keep values that do appear in filter_df_w_read_ids:
        filtered_df = pd.merge(self.df, filter_df_w_read_ids, on="read_id", how="inner")
        if overwrite_df:
            self.df = filtered_df
        return filtered_df

    def save_to_fastq(self, output_path=None) -> str:
        from csv import QUOTE_NONE
        output_df = self.df.copy()
        output_df["read_id"] = ">@" + output_df["read_id"] + " " + output_df["comment"]
        if output_path:
            output_df.drop("comment", axis=1).to_csv(output_path,
                                                     index=False,
                                                     header=False,
                                                     sep="\n",
                                                     quoting=QUOTE_NONE)
            return f"Saved to {output_path}"
        else:
            return output_df.drop("comment", axis=1).to_csv(index=False,
                                                            header=False,
                                                            sep="\n",
                                                            quoting=QUOTE_NONE)


class NanoporeRun:
    def __init__(self, run_name: str = None, run_nickname: str = None,
                 run_dir_name: str = None, run_dir_abs_path: str = None):
        self.nanopore_run_soft_link_dir = f"/data16/marcus/nanoporeSoftLinks"
        self.run_name = run_name
        self.run_nickname = run_nickname
        self.run_dir_name = run_dir_name
        self.run_dir_abs_path = run_dir_abs_path

        if sum([bool(self.run_name), bool(self.run_nickname), bool(self.run_dir_name),
                bool(self.run_dir_abs_path)]) != 1:
            raise ValueError("Must provide one of run_name, run_nickname, run_dir_name, run_dir_abs_path")

        if self.run_name:
            self.run_nickname = CONVERSION_DICT[self.run_name]
            self.run_dir_abs_path = OUTPUT_DIR_DICT[self.run_name]
            self.run_dir_name = self.run_dir_abs_path.split("/")[-2]
        elif self.run_nickname:
            self.run_name = REV_CONVERSION_DICT[self.run_nickname]
            self.run_dir_abs_path = OUTPUT_DIR_DICT[self.run_name]
            self.run_dir_name = self.run_dir_abs_path.split("/")[-2]
        elif self.run_dir_name:
            self.run_dir_abs_path = f"{self.nanopore_run_soft_link_dir}/{self.run_dir_name}/output_dir"
            self.run_name = REV_OUTPUT_DIR_DICT[self.run_dir_abs_path]
            self.run_nickname = CONVERSION_DICT[self.run_name]
        elif self.run_dir_abs_path:
            self.run_dir_name = self.run_dir_abs_path.split("/")[-2]
            self.run_name = REV_OUTPUT_DIR_DICT[self.run_dir_abs_path]
            self.run_nickname = CONVERSION_DICT[self.run_name]

        self.run_dir = Path(self.run_dir_abs_path)
        if not self.run_dir.exists():
            raise ValueError(f"Run directory does not exist: {self.run_dir_abs_path}")

        self.cat_dir = self.run_dir / "cat_files"
        self.merge_dir = self.run_dir / "merge_files"
        self.fastq_dir = self.run_dir / "fastqs"
        self.nanopolish_dir = self.run_dir / "nanopolish"
        self.feature_counts_dir = self.run_dir / "featureCounts"
        self.logs_dir = self.run_dir / "logs"
        self.run_date = self.run_dir_name.split("_")[0]

        self.cat_files = list(self.cat_dir.glob("cat*"))
        # print([file.name for file in self.cat_files])

        self.had_standards = False
        self.merge_parquet_path_dict = {}
        for merge_file in self.merge_dir.glob("*.parquet"):
            short_name = merge_file.name[7:-8]
            if short_name in self.merge_parquet_path_dict.keys():
                self.merge_parquet_path_dict[short_name].append(merge_file)
            else:
                self.merge_parquet_path_dict[short_name] = [merge_file]
            if short_name == "mergedOnReads.plusStandards":
                self.had_standards = True
        self.merge_df_dict = {}
        self.mergedOnReads_df = None
        self.compressedOnGenes_df = None

        self.bam_path = self.cat_dir / "cat.sorted.mappedAndPrimary.bam"
        self.bam = None

        self.fastq_path = self.cat_dir / "cat.fastq"
        self.fastq = None

    def print_run_dirs(self) -> None:
        print(f"Run name: {self.run_name}\n"
              f"Run nickname: {self.run_nickname}\n"
              f"Run directory: {self.run_dir}\n")
        print(f"Run Dir contents:")
        for output_dir in self.run_dir.iterdir():
            print(f"--> {output_dir.name}")
            if len(list(output_dir.iterdir())) <= 15:
                for file in output_dir.iterdir():
                    print(f"    --> {file.name}")
            else:
                print(f"    --> ({len(list(output_dir.iterdir()))} files)")

    def _load_merge_parquet(self, file_midfix: str = "") -> pd.DataFrame:
        if file_midfix not in self.merge_parquet_path_dict.keys():
            raise ValueError(f"Short name {file_midfix} not found in merge_parquet_dict")
        if not self.had_standards and file_midfix == "mergedOnReads.plusStandards":
            raise ValueError(f"Run {self.run_name} did not have standards, try short_name='mergedOnReads'")

        if len(self.merge_parquet_path_dict[file_midfix]) == 1:
            merge_file = self.merge_parquet_path_dict[file_midfix][0]
        else:  # If there are multiple parquet files, we need to pick the most recent one
            merge_file = sorted(self.merge_parquet_path_dict[file_midfix], key=lambda x: x.stat().st_mtime)[-1]
        print(f"Loading {merge_file.name}...", end="")
        self.merge_df_dict[file_midfix] = pd.read_parquet(merge_file)
        print(f" Done. Loaded {self.merge_df_dict[file_midfix].shape[0]:,} rows.")
        if "mergedOnReads" in file_midfix:
            self.mergedOnReads_df = self.merge_df_dict[file_midfix]
        elif file_midfix == "compressedOnGenes":
            self.compressedOnGenes_df = self.merge_df_dict[file_midfix]
        return self.merge_df_dict[file_midfix]

    def load_mergedOnReads(self) -> pd.DataFrame:
        if self.had_standards:
            return self._load_merge_parquet(file_midfix="mergedOnReads.plusStandards")
        else:
            return self._load_merge_parquet(file_midfix="mergedOnReads")

    def load_compressedOnGenes(self) -> pd.DataFrame:
        return self._load_merge_parquet(file_midfix="compressedOnGenes")

    def load_bam(self, subsample=-1) -> SamOrBamFile:
        if not self.bam_path.exists():
            raise ValueError(f"BAM file does not exist: {self.bam_path}")
        if self.bam is None:
            self.bam = SamOrBamFile(str(self.bam_path), subsample=subsample)
        return self.bam

    def load_fastq(self, subsample=-1) -> FastqFile:
        if not self.fastq_path.exists():
            raise ValueError(f"FASTQ file does not exist: {self.fastq_path}")
        if self.fastq is None:
            self.fastq = FastqFile(str(self.fastq_path), subset_entries=subsample)
        return self.fastq

    def plot_standards(self, save=False, show=True, show_ambiguous=False,
                       show_non_standards=False, **kwargs) -> (plt.Figure, plt.Axes):
        if not self.had_standards:
            raise ValueError(f"Run {self.run_name} did not have standards")
        if self.mergedOnReads_df is None:
            self.load_mergedOnReads()
        if "figsize" not in kwargs.keys():
            kwargs["figsize"] = (12, 8)

        fig, ax = plt.subplots(**kwargs)
        assignments_df = self.mergedOnReads_df[['read_id', 'assignment', 'chr_id']]
        if not show_ambiguous:
            assignments_df = assignments_df[~assignments_df.assignment.str.contains("Ambiguous")]
        if not show_non_standards:
            assignments_df = assignments_df[~assignments_df.assignment.str.contains("NotAStandard")]
        assignments_df = assignments_df[['assignment']].value_counts(ascending=False).reset_index()
        sea.barplot(assignments_df, ax=ax,
                    x='assignment', y=0)
        ax.set_title(f"{self.run_nickname} Standards")
        ax.set_ylabel("Number of Reads")
        ax.set_xlabel("Assigned Standard")
        plt.tight_layout()
        if save:
            fig.savefig(f"{self.run_dir}/standards.png")
        if show:
            plt.show()
        return fig, ax

    def plot_read_length(self, save=False, show=True, **kwargs) -> (plt.Figure, plt.Axes):
        if self.mergedOnReads_df is None:
            self.load_mergedOnReads()
        if "figsize" not in kwargs.keys():
            kwargs["figsize"] = (12, 8)

        fig, ax = plt.subplots(**kwargs)
        sea.histplot(self.mergedOnReads_df['read_length'], ax=ax, kde=True)
        ax.set_title(f"{self.run_nickname} Read Lengths")
        ax.set_ylabel("Frequency")
        ax.set_xlabel("Read Length")
        plt.tight_layout()
        if save:
            fig.savefig(f"{self.run_dir}/read_lengths.png")
        if show:
            plt.show()
        return fig, ax

    def plot_tail_length_by_standard(self, save=False, show=True,
                                     y_lims=(0, 200), **kwargs) -> (plt.Figure, plt.Axes):
        if not self.had_standards:
            raise ValueError(f"Run {self.run_name} did not have standards")
        if self.mergedOnReads_df is None:
            self.load_mergedOnReads()
        if "figsize" not in kwargs.keys():
            kwargs["figsize"] = (12, 8)
        
        plot_df = self.mergedOnReads_df[['assignment', 'polya_length']]
        plot_df = plot_df[~plot_df.assignment.str.contains("Ambiguous")]
        
        fig, ax = plt.subplots(**kwargs)
        sea.boxplot(plot_df, ax=ax,
                    x='assignment', y='polya_length')
        ax.set_title(f"{self.run_nickname} Tail Lengths by Standard")
        ax.set_ylabel("Tail Length")
        ax.set_xlabel("Assigned Standard")
        ax.set_ylim(y_lims[0], y_lims[1])
        
        plt.tight_layout()
        if save:
            fig.savefig(f"{self.run_dir}/tail_lengths_by_standard.png")
        if show:
            plt.show()
        return fig, ax
    
    # TODO: Add ability to tag standards reads with "perfect" or "imperfect" based on presence of barcode


class NanoJAMDF(pd.DataFrame):
    pass


def load_and_merge_lib_parquets(lib_list, genomeDir=f"/data16/marcus/genomes/elegansRelease100/",
                                drop_unassigned=True, drop_failed_polya=True,
                                drop_sub_n=5, keep_transcript_info=False,
                                read_pos_in_groupby=False, add_nucleotide_fractions=False,
                                add_tail_groupings=True, tail_groupings_cutoff=50,
                                pass_list_columns=False,
                                use_josh_assignment=True,
                                group_by_t5=False) -> [pd.DataFrame, pd.DataFrame]:
    concat_df = library_reads_df_load_and_concat(lib_list, genomeDir=genomeDir,
                                                 drop_unassigned=drop_unassigned, drop_failed_polya=drop_failed_polya,
                                                 keep_transcript_info=keep_transcript_info, group_by_t5=group_by_t5,
                                                 use_josh_assignment=use_josh_assignment)

    compressed_df = compress_concat_df_of_libs(concat_df,
                                               drop_sub_n=drop_sub_n, read_pos_in_groupby=read_pos_in_groupby,
                                               add_nucleotide_fractions=add_nucleotide_fractions,
                                               add_tail_groupings=add_tail_groupings,
                                               tail_groupings_cutoff=tail_groupings_cutoff,
                                               pass_list_columns=pass_list_columns,
                                               keep_transcript_info=keep_transcript_info, group_by_t5=group_by_t5)

    return concat_df, compressed_df


def library_reads_df_load_and_concat(lib_list, genomeDir=f"/data16/marcus/genomes/elegansRelease100/",
                                     drop_unassigned=True, drop_failed_polya=True,
                                     keep_transcript_info=False,
                                     use_josh_assignment=True,
                                     group_by_t5=False):
    # TODO: Why am I not trusting the gene_assignment from featureCounts? I am using Josh's method here?!
    #       I really think that this is what is throwing off all the errors w/ Parissa's unc-54.
    #       Primarily b/c I can see the correct gene assignment in the *_mergedOnReads.parquet files!!
    if use_josh_assignment:
        # Initially load the read assignments file:
        read_assignment_path = find_newest_matching_file(f"{genomeDir}/*.allChrs.parquet")
        print(f"Loading readAssignments file from: {read_assignment_path}...", end=' ')
        read_assignment_df = pd.read_parquet(read_assignment_path)
        print("Done.")
    # Loop through each library name in the list and for each:
    #   1. Load the Parquet
    #   2. Merge this w/ Josh's assign reads based on chr_pos (if wanted!!)
    #   3. Create a column to retain the library identity
    #   3. Concatenate these dataframe into one large dataframe
    #       NOTE: This structure can be seperated again based on
    #       the "lib" column added in the previous step
    print(f"Looking for files for libraries: {lib_list}")
    path_dict = pick_libs_return_paths_dict(lib_list,
                                            output_dir_folder="merge_files",
                                            file_midfix="_mergedOnReads",
                                            file_suffix="parquet")
    # Dictionary to hold the library names and dataframes before the merge
    df_dict = {}
    for library_name, parquet_path in path_dict.items():
        # Load each individual library:
        print(f"Loading parquet for {library_name} lib...", end=' ')
        lib_df = pd.read_parquet(parquet_path)
        print("Done.")

        # Make sure that 5' end adjustments happened, do so if not!
        lib_df = adjust_5_ends(lib_df)
        # Store the library dataframe in the temporary dictionary
        df_dict[library_name] = lib_df
    # This is a cute way to quickly merge all of these dfs into one, while retaining lib info.
    #   B/c I am still a little scared of MultiIndexed dataframes, I used the reset_index steps
    #   to push the mutliindex back into columns. Maybe someday I'll use the multiindex!
    multi_df = pd.concat(df_dict.values(),
                         keys=df_dict.keys())
    multi_df.index.set_names(("lib", "old_index"), inplace=True)
    concat_df = multi_df.reset_index(level="lib").reset_index(drop=True)
    if use_josh_assignment:
        # Some major issues with ugly column names coming through, plan to clean them up:
        read_assignment_cols = read_assignment_df.columns.to_list()
        read_assignment_cols.remove('chr_id')
        read_assignment_cols.remove('chr_pos')
        # Only retain columns that don't show up in both the read assignment df and the super merge df:
        read_assignment_cols_to_drop = [col for col in read_assignment_cols if col in concat_df.columns]
        # Drop those 'unique' columns
        concat_df.drop(columns=read_assignment_cols_to_drop,
                       inplace=True)
        print(f"Starting assignment merge . . .", end="")
        # Add read assignments w/ josh's read_assignment dataframe
        concat_df = concat_df.merge(read_assignment_df, on=["chr_id", "chr_pos"],
                                    how="left", suffixes=("_originalOutput", ""))
        print(f"\rFinished assignment merge!          ")
    else:
        print(f"Skipping assignment with Josh method and relying on whatever assignment was made by the pipeline!")
    # To further clean up mixed columns, just retain the ones we care about!
    keep_columns = ["lib",
                    "read_id",
                    "chr_id",
                    "chr_pos",
                    "original_chr_pos",
                    "gene_id",
                    "gene_name",
                    "cigar",
                    "sequence",
                    "polya_length",
                    "strand",
                    ]
    if keep_transcript_info and use_josh_assignment:
        # Add transcript info to the columns we care about if requested
        for col in ["transcript_id", "to_start", "to_stop"]:
            keep_columns.append(col)
    else:
        print(f"Not keeping transcript information. . . (not using Josh assignment method will also force this!)")
    if group_by_t5:
        # Add 5TERA information to the columns we care about if requested
        keep_columns.append('t5')
    else:
        print(f"Not keeping 5TERA adapter information (if even present)...")
    # Drop unnecessary columns, reassess for duplicates.
    #   i.e. reads that mapped to two transcripts will
    #        otherwise be dups if we drop transcript info!
    print(f"Dropping duplicate columns. . .", end='')
    concat_df = concat_df[keep_columns].drop_duplicates()
    print(f"\rFinished dropping dup. columns.")
    # Only retain polyA passed reads if requested
    if drop_failed_polya:
        print(f"\nRead counts pre-failed-polyA call drop:   {concat_df.shape[0]}")
        concat_df = concat_df[~concat_df["polya_length"].isna()]
        print(f"Read counts post-failed-polyA call drop:  {concat_df.shape[0]}")
    # Post-assignment cleanups:
    if drop_unassigned:
        print(f"Read counts post gene assignment:  {concat_df.shape[0]}")
        concat_df = concat_df[~concat_df["gene_id"].isna()].reset_index(drop=True)
        print(f"Read counts post unassigned drop:  {concat_df.shape[0]}")
        if keep_transcript_info:
            concat_df = concat_df.astype({"to_stop": "int64",
                                          "to_start": "int64"})
    # Add a read length column!
    concat_df["read_length"] = concat_df["sequence"].apply(len)
    return concat_df


def compress_concat_df_of_libs(concat_df, drop_sub_n=5, read_pos_in_groupby=False, add_nucleotide_fractions=False,
                               keep_transcript_info=False, pass_list_columns=False,
                               add_tail_groupings=True, tail_groupings_cutoff=50,
                               group_by_t5=False):
    # Create the groupby dataframe:
    groupby_col_list = ["lib", "chr_id", "gene_id", "gene_name"]
    print(f"Creating groupby dataframe merged on: {groupby_col_list}")
    if keep_transcript_info:
        print(f"\t+ [transcript_id]")
        groupby_col_list.append("transcript_id")
    if group_by_t5:
        print(f"\t+ [t5] tag")
        groupby_col_list.append("t5")

    # Holy crap, the observed=True helps to keep this from propagating out to 129,151,669,691,968 rows...
    groupby_obj = concat_df.groupby(groupby_col_list, observed=True)

    # Change the compressed prefix so that I am count gene hits or transcript hits, depending on set up!
    if not keep_transcript_info:
        compressed_prefix = "gene"
    else:
        compressed_prefix = "transcript"

    tqdm.pandas(desc=f"Counting reads per {compressed_prefix}")
    compressed_df = groupby_obj["read_id"].progress_apply(len).to_frame(name=f"{compressed_prefix}_hits")

    compressed_df["mean_polya_length"] = groupby_obj["polya_length"].mean()
    compressed_df["median_polya_length"] = groupby_obj["polya_length"].median()

    compressed_df["mean_read_length"] = groupby_obj["read_length"].mean()
    compressed_df["median_read_length"] = groupby_obj["read_length"].median()

    if read_pos_in_groupby:
        tqdm.pandas(desc=f"Storing stop distances as lists")
        compressed_df['stop_distances'] = groupby_obj["to_stop"].progress_apply(list).to_frame(name="stop_distances")
        tqdm.pandas(desc=f"Storing start distances as lists")
        compressed_df['start_distances'] = groupby_obj["to_start"].progress_apply(list).to_frame(name="start_distances")
    if pass_list_columns:
        tqdm.pandas(desc=f"Storing polyA lengths as lists")
        compressed_df['polya_lengths'] = groupby_obj["polya_length"].progress_apply(list).to_frame(name="polya_lengths")
        tqdm.pandas(desc=f"Storing read lengths as lists")
        compressed_df['read_lengths'] = groupby_obj["read_length"].progress_apply(list).to_frame(name="read_lengths")
    if add_tail_groupings:
        print(f"Adding tail length groupings information, with a tail cutoff at {tail_groupings_cutoff} nts...",
              end=' ')
        tqdm.pandas(desc=f"Counting number of tails shorter than {tail_groupings_cutoff} for each gene")
        compressed_df[f"sub_{tail_groupings_cutoff}_tails"] = groupby_obj.polya_length. \
            progress_apply(lambda group_series: group_series[group_series <= tail_groupings_cutoff].count())
        compressed_df[f'frac_sub_{tail_groupings_cutoff}_tails'] = compressed_df[f"sub_{tail_groupings_cutoff}_tails"] \
                                                                   / compressed_df[f"{compressed_prefix}_hits"]
        compressed_df['tail_groupings_group'] = pd.cut(compressed_df[f'frac_sub_{tail_groupings_cutoff}_tails'],
                                                       bins=[0.0, 0.1, 0.45, 1.0], labels=['long_tailed',
                                                                                           'ungrouped',
                                                                                           'short_tailed'],
                                                       include_lowest=True)
        print(f"Done.")
    # RPM and fractional hits calculations
    # Need to first create columns of NA values, tobe overwritten
    compressed_df[f"{compressed_prefix}_rpm"] = pd.NA
    compressed_df[f"{compressed_prefix}_frac_hits"] = pd.NA
    if group_by_t5:
        compressed_df[f"{compressed_prefix}_t5group_rpm"] = pd.NA
    # Only look at one library at a time (so the normalization is per lib not whole df)
    for lib in compressed_df.index.unique(level='lib').to_list():
        # Create the 'norm_factor' which will be the total # of read hits in that lib
        norm_factor = compressed_df.query(f"lib == '{lib}'")[f"{compressed_prefix}_hits"].sum()
        # Turn the total number of read hits into the 'million of read hits'
        rpm_norm_factor = norm_factor / 1000000
        # For each library divide gene_hits by the rpm norm factor to get rpm
        rpm_series = compressed_df.query(f"lib == '{lib}'")[f"{compressed_prefix}_hits"] / rpm_norm_factor
        # Use a series fill, so that we can fill that library's part of the DF without effecting others.
        #       This step is reliant on the gene_rpm_series carrying over the indexes.
        #       That's how the 'fillna()' matches the right values!
        compressed_df[f"{compressed_prefix}_rpm"] = compressed_df[
            f"{compressed_prefix}_rpm"].fillna(value=rpm_series)
        # Same as above, but with fraction of hits, rather than a rpm calc (practically same thing)
        gene_frac_hits_series = compressed_df.query(f"lib == '{lib}'")[f"{compressed_prefix}_hits"] / norm_factor
        compressed_df[f"{compressed_prefix}_frac_hits"] = compressed_df[f"{compressed_prefix}_frac_hits"]. \
            fillna(value=gene_frac_hits_series)
        if group_by_t5:
            # We can also calculate an adapted-specific RPM:
            for adapted_or_not in ["+", "-"]:
                norm_factor = compressed_df.query(f"lib == '{lib}'") \
                    .query(f"t5 == '{adapted_or_not}'")[f"{compressed_prefix}_hits"].sum()
                rpm_norm_factor = norm_factor / 1_000_000
                rpm_series = compressed_df.query(f"lib == '{lib}'")[f"{compressed_prefix}_hits"] / rpm_norm_factor
                compressed_df[f"{compressed_prefix}_t5group_rpm"] = compressed_df[
                    f"{compressed_prefix}_rpm"].fillna(value=rpm_series)
    # Requirement for min number of gene/transcript hits
    if isinstance(drop_sub_n, int):
        print(f"Gene counts pre sub-{drop_sub_n} {compressed_prefix}_hits drop:  {compressed_df.shape[0]}")
        compressed_df = compressed_df[compressed_df[f"{compressed_prefix}_hits"] >= drop_sub_n]
        print(f"Gene counts post sub-{drop_sub_n} {compressed_prefix}_hits drop:  {compressed_df.shape[0]}")
    # Reset index at the end,
    #   we didn't retain any info w/ the index, so it doesn't help much
    compressed_df = compressed_df.reset_index()
    if add_nucleotide_fractions:
        if not keep_transcript_info:
            print(f"Adding nucleotide content information to genes!")
            path_to_gc = "/data16/marcus/genomes/elegansRelease100/" \
                         "Caenorhabditis_elegans.WBcel235.cdna.all.fa.GCcontent.genes.parquet"
            gc_df = pd.read_parquet(path_to_gc).drop(columns=["chr_id"])
            compressed_df = compressed_df.merge(gc_df, on=["gene_id", "gene_name"], how="left")
        else:
            print(f"Adding nucleotide content information to transcripts!")
            path_to_gc = "/data16/marcus/genomes/elegansRelease100/" \
                         "Caenorhabditis_elegans.WBcel235.cdna.all.fa.GCcontent.transcripts.parquet"
            gc_df = pd.read_parquet(path_to_gc).drop(columns=["chr_id"])
            compressed_df = compressed_df.merge(gc_df, on=["gene_id", "gene_name", "transcript_id"], how="left")
    return compressed_df


def live_cmd_call(command):
    with Popen(command, stdout=PIPE, shell=True,
               bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end="")
    if p.returncode != 0:
        raise CalledProcessError(p.returncode, p.args)
    else:
        return p.returncode


def sam_or_bam_class_testing():
    test_output_dir = "/data16/marcus/working/211101_nanoporeSoftLinks/220131_nanoporeRun_totalRNA_0639_L3_third/output_dir"
    test_sam_path = f"{test_output_dir}/cat_files/cat.sorted.mappedAndPrimary.bam"
    # test_sam_path = "/data16/marcus/working/211101_nanoporeSoftLinks/" \
    #                 "211118_nanoporeRun_totalRNA_5108_xrn-1-KD_5TERA/" \
    #                 "output_dir/cat_files/cat.sorted.mappedAndPrimary.sam"
    sam = SamOrBamFile(test_sam_path,
                       # subsample=50000,
                       )
    polya_df = pd.read_csv(f"{test_output_dir}/nanopolish/polya.passed.tsv", sep="\t")
    # For some god-awful reason the chr_pos in polyA are -1 to those in the SAM file:
    polya_df["position"] = polya_df["position"] + 1
    polya_df = polya_df.rename(columns={"readname": "read_id",
                                        "qc_tag": "qc_tag_polya",
                                        "position": "chr_pos",
                                        "contig": "chr_id",
                                        "polya_length": "pA"})
    polya_df = polya_df[['read_id', 'chr_id', 'chr_pos', 'pA']]
    sam.df = pd.merge(sam.df, polya_df, on=['read_id', 'chr_id', 'chr_pos'], how='left')
    sam.df.pA.fillna(0, inplace=True)
    sam.tag_columns.append('pA')
    # sam.to_sam(test_sam_path + '.resave.sam', escape_char="~")
    # sam.to_sam(output_path=test_sam_path.rstrip('.sam') + '.resave.sam')
    print(sam)


def boolDF_to_upsetPlot(input_df: pd.DataFrame,
                        show_percentages=False, show_counts=True,
                        file_name=None, min_subset_size=0, min_degree=0,
                        sort_by='cardinality', sort_categories_by='cardinality'):
    import upsetplot
    import matplotlib.pyplot as plt
    upset_format_data = upsetplot.from_indicators(lambda lambda_df: lambda_df.select_dtypes(bool), data=input_df)
    fig = plt.figure()
    upset = upsetplot.UpSet(upset_format_data, sort_by=sort_by, sort_categories_by=sort_categories_by,
                            min_subset_size=min_subset_size, min_degree=min_degree, show_percentages=show_percentages,
                            show_counts=show_counts)
    upset.plot(fig=fig)
    if isinstance(file_name, str):
        fig.savefig(file_name)


# "riboD", "totalRNA", "totalRNA2", "polyA", "polyA2",
# "xrn-1", "xrn-1-5tera", "pTRI-stds", "xrn-1-5tera-smg-6", "pTRI-stds-tera3"
def pick_libs_return_paths_dict(lib_list: list, output_dir_folder="merge_files",
                                file_midfix="_mergedOnReads", file_suffix="parquet",
                                return_all: bool = False, ignore_unmatched_keys: bool = False) -> dict:
    output_dir_dict = OUTPUT_DIR_DICT
    if not isinstance(lib_list, list) or not isinstance(lib_list, tuple):
        raise NotImplementedError(f"Please pass a list/tuple of library keys, "
                                  f"you passed a {type(lib_list)}. "
                                  f"If you only want one value, "
                                  f"please use the pick_lib_return_path() method.")
    if return_all:
        lib_list = output_dir_dict.keys()
        lib_list = [lib for lib in lib_list if lib not in ["polyA", "totalRNA"]]
    file_suffix = file_suffix.strip(".")
    return_dict = {}
    for lib_key in lib_list:
        if lib_key in output_dir_dict.keys():
            output_dir = output_dir_dict[lib_key]
            file_path = f"{output_dir}/{output_dir_folder}/*{file_midfix}.{file_suffix}"
            print(f"Looking for file for {lib_key}, at {file_path}...", end=" ")
            return_dict[lib_key] = find_newest_matching_file(file_path)
            print(f"File Found.")
        else:
            if ignore_unmatched_keys:
                continue
            else:
                raise FileNotFoundError(f"Could not find a matching library tag for {lib_key}\n"
                                        f"Your keys must be found in this list: {list(output_dir_dict.keys())}")
    if not return_dict:
        raise KeyError(f"No matching library keys found, please use keys from the following list: "
                       f"{list(output_dir_dict.keys())}. "
                       f"You passed the following keys: "
                       f"{lib_list}")
    return return_dict


def pick_lib_return_path(lib_key, output_dir_folder="merge_files",
                         file_midfix="_mergedOnReads", file_suffix="parquet", ) -> str:
    [(lib_key, lib_path)] = pick_libs_return_paths_dict([lib_key],
                                                        file_suffix=file_suffix,
                                                        file_midfix=file_midfix,
                                                        output_dir_folder=output_dir_folder).items()
    return lib_path


def load_read_assignments(assignment_file_parquet_path) -> pd.DataFrame:
    print(f"Loading read assignment file from: {assignment_file_parquet_path} ", end="")
    try:
        read_assignment_df = pd.read_parquet(assignment_file_parquet_path)
    except FileNotFoundError:
        raise FileNotFoundError(f"Couldn't find parquet file! Maybe try to run the method: "
                                f"parse_read_assignment_allChrs_txt in nanoporePipelineCommon.py or "
                                f"first: prepareReadAssignmentFily3.py in the ArribereLab git!!")
    return read_assignment_df


def old_assign_w_josh_method(reads_df, genomeDir):
    def merge_on_chr_pos(read_assignment_df: pd.DataFrame, reads_df: pd.DataFrame) -> pd.DataFrame:
        print(f"Merging read assignments and reads at {get_dt(for_print=True)}")
        merge_df = reads_df.merge(read_assignment_df, on=["chr_id", "chr_pos"],
                                  how="left", suffixes=("_fromFeatureCounts",
                                                        ""))
        # merge_df = merge_df[~(merge_df.gene_id_fromReads.isna() & merge_df.gene_id_fromAssign.isna())]
        # below call drops reads that don't get assigned by Josh's tool
        print(f"Reads unassigned by Josh method: {merge_df[~merge_df.gene_id.isna()].shape[0]} "
              f"/{merge_df.shape[0]}, [Dropping these!]")
        merge_df = merge_df[~merge_df.gene_id.isna()]

        print(f"Reads with different assignments between Josh method and FeatureCounts: "
              f"{merge_df[merge_df.strand_fromReads != merge_df.strand_fromAssign].shape[0]} "
              f"/{merge_df.shape[0]}, [Dropping these!]")
        print(f"Reads with the same assignments between Josh method and FeatureCounts: "
              f"{merge_df[merge_df.strand_fromReads == merge_df.strand_fromAssign].shape[0]} "
              f"/{merge_df.shape[0]}, [Dropping these!]")
        merge_df = merge_df[merge_df.strand_fromReads == merge_df.strand_fromAssign]
        print(f"Done merging at {get_dt(for_print=True)}")
        return merge_df

    read_assignments_path = find_newest_matching_file(f"{genomeDir}/*.allChrs.parquet")
    read_assignments_df = load_read_assignments(read_assignments_path)
    print("Finished loading files!")
    # for df in [reads_df, read_assignments_df]:
    #     print(df.info())
    merged_df = merge_on_chr_pos(read_assignments_df, reads_df)
    return merged_df


def assign_with_josh_method(merged_on_reads_df: pd.DataFrame, genomeDir: str,
                            keepMultipleTranscriptInfo=False, save_it=False,
                            add_names=False) -> pd.DataFrame:
    merged_on_reads_df = adjust_5_ends(merged_on_reads_df)
    print(f"Using Josh's read assignment method w/ 5'ends!")
    if keepMultipleTranscriptInfo:
        read_assignment_path = find_newest_matching_file(f"{genomeDir}/*.allChrs.parquet")
    else:
        read_assignment_path = find_newest_matching_file(f"{genomeDir}/*.allChrsLite.parquet")
    print(f"Loading readAssignment file from {read_assignment_path}")
    read_assignment_df = pd.read_parquet(read_assignment_path)
    if keepMultipleTranscriptInfo:
        read_assignment_df = read_assignment_df.astype({"to_stop": "int64",
                                                        "to_start": "int64"})
    print(f"Loaded  readAssignment file from {read_assignment_path}")

    print(f"Starting merge of dataframes to make assignments")
    merged_on_reads_and_assigned_df = merged_on_reads_df.merge(read_assignment_df, on=["chr_id", "chr_pos"],
                                                               how="left", suffixes=("",
                                                                                     "_fromJoshAssign"))
    print(f"Finished merge of dataframes to make assignments")
    return merged_on_reads_and_assigned_df


def adjust_5_ends(df: pd.DataFrame):
    import re
    from tqdm import tqdm
    tqdm.pandas()

    def _flip_neg_strand_genes(chr_position: int, cigar: str, strand: str) -> int:
        if strand == "+":
            read_end = chr_position
            return read_end
        else:  # if strand == "-":
            parsed_cigar = re.findall(rf'(\d+)([MDNSIX])', cigar)
            mdn_nums = [int(num) for num, char in parsed_cigar if char in "MDN"]
            read_end = chr_position + sum(mdn_nums)
            return read_end

    if "original_chr_pos" not in df.columns.to_list():
        print(f"\nMaking adjustments for 5' ends:")
        df["original_chr_pos"] = df["chr_pos"]
        df["chr_pos"] = df.progress_apply(lambda read: _flip_neg_strand_genes(read["original_chr_pos"],
                                                                              read["cigar"],
                                                                              read["strand"]),
                                          axis=1)
    else:
        print(f"'original_chr_pos' column already found in dataframe, skipping adjustment for 5'ends!")
    return df


def gene_names_to_gene_ids(
        parquet_path: str = "/data16/marcus/genomes/elegansRelease100/"
                            "Caenorhabditis_elegans.WBcel235.100.gtf.parquet") -> pd.DataFrame:
    df = pd.read_parquet(parquet_path)[["gene_name", "gene_id", "chr"]].drop_duplicates(ignore_index=True)
    return df


def get_dt(for_print=False, for_file=True):
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


def gtf_to_df(gtf_path: str) -> pd.DataFrame:
    gtf_header = ["chr",
                  "source",
                  "feature",
                  "start",
                  "end",
                  "score",
                  "strand",
                  "frame",
                  "attributes"]
    gtf_attributes = ["gene_name",
                      "gene_id",
                      "gene_version",
                      "gene_source",
                      "gene_biotype",
                      "transcript_id",
                      "transcript_source",
                      "transcript_biotype",
                      "exon_number",
                      "exon_id",
                      ]
    gtf_dtypes = {'chr': 'O',
                  'source': 'O',
                  'feature': 'O',
                  'start': 'Int64',
                  'end': 'Int64',
                  'score': 'O',
                  'strand': 'category',
                  'frame': 'O',
                  'gene_name': 'O',
                  'gene_id': 'O',
                  'gene_version': 'Int64',
                  'gene_source': 'O',
                  'gene_biotype': 'O',
                  'transcript_id': 'O',
                  'transcript_source': 'O',
                  'transcript_biotype': 'O',
                  'exon_number': 'float64',
                  'exon_id': 'O'}
    gtf_df = pd.read_csv(gtf_path, sep="\t", comment="#", names=gtf_header)
    for attribute in gtf_attributes:
        gtf_df[attribute] = gtf_df["attributes"].str.extract(rf"""{attribute} "(.*?)";""")
    gtf_df.drop(columns=["attributes"], inplace=True)
    gtf_df = gtf_df.astype(gtf_dtypes)
    return gtf_df


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


def parse_read_assignment_allChrs_txt(assignment_file_txt, multi_row_isoforms=False,
                                      save_file=True) -> pd.DataFrame:
    # Script to load josh's read assignment file from the *.allChrs.txt format,
    #   it can return a dataframe, but it will also (optionally) save a new parquet file!
    # Majorly rewritten on 12/14/2021, will now save a parquet with lists holding transcript/
    #   isoform information. Easily expanded with the df.explode() function!

    print(f"Loading read assignment file from: {assignment_file_txt}")
    df = pd.read_csv(assignment_file_txt, sep="\t",
                     names=["chr_pos", "gene_id_strand", "transcript_info"])
    print(f"\tLoaded file to dataframe.")

    df[["chr_id", "chr_pos"]] = df["chr_pos"].str.split("_", expand=True)
    print(f"\tParsed chr_id and chr_pos info.")

    df[["gene_id", "strand"]] = df["gene_id_strand"].str.split(":", expand=True)
    df.drop(columns="gene_id_strand", inplace=True)
    print(f"\tParsed gene_id and strand info.")

    # The below section is CRAZY, but allow me to 'quickly' refactor
    #   the transcript/start/stop information to lists: >>>
    df['transcript_info'] = df.transcript_info.str.split('|')  # First just split up multi-isoform spots
    print(f"\tSplit isoform info.")
    # Then split each isoform into the: trans_id, start, & stop information:
    df['transcript_info'] = df.transcript_info.apply(lambda tran_list:
                                                     np.array([tran.split(':') for tran in tran_list]).T)
    print(f"\tExtracted transcript, to_start and to_stop info.")
    # Then convert these nested listed into three columns of lists
    df[['transcript_id', 'to_start', 'to_stop']] = pd.DataFrame.from_records(df.transcript_info.to_numpy(),
                                                                             index=df.index,
                                                                             columns=['transcript_id',
                                                                                      'to_start',
                                                                                      'to_stop'])
    print(f"\tSaved transcript, to_start and to_stop info to new columns.")
    if multi_row_isoforms:
        # If we want each isoform to have it's own row, we can use explode!!
        df = df.explode(['transcript_id', 'to_start', 'to_stop'])
        additional_tag = ".isoformsAsRows"
        dtypes = {
            "chr_pos": "uint32",
            "to_start": "int32",
            "to_stop": "int32",
            "strand": "category",
            "chr_id": "category",
        }
        print(f"\tExpanded transcript, to_start and to_stop info to individual rows.")
    else:
        additional_tag = ""
        dtypes = {
            "chr_pos": "uint32",
            "to_start": "object",
            "to_stop": "object",
            "strand": "category",
            "chr_id": "category",
        }
        print(f"\tSkipping expansion of transcript, to_start and to_stop info to individual rows.")
    # <<< Maybe not that bad...? \s
    # Change datatypes:
    df = df.astype(dtypes)
    # Reorder columns:
    df = df[["chr_id", "chr_pos", "gene_id", "strand", "transcript_id", "to_start", "to_stop"]]
    print(". ", end="")
    if save_file:
        parquet_path = assignment_file_txt.rstrip(".txt") + additional_tag + ".parquet"
        df.to_parquet(parquet_path)
        print(f"Saved parquet to: {parquet_path}")
    return df


if __name__ == '__main__':
    run = NanoporeRun(run_nickname="newerS5")
    # run.print_run_dirs()
    run.plot_read_length()
    run.plot_tail_length_by_standard(y_lims=(0, 100))
