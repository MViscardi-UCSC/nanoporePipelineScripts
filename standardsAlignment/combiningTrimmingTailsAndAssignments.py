"""
testingCombinedTailsTrimmingAndAssignments.ipynb
Marcus Viscardi,    December 05, 2022

There are some questions regarding reads that get called as having zero tails by tailfindr.
It seems like ~80% of the 05A standards are called as having ZERO tails, could this be
because tailfindr needs some minimum length to recognize the homopolymer stretch?

With this in mind the hope would be that that minimum length is somewhere around the range
that the guppy_basecaller is able to pick out! If this is the case, we can use Porechop to
remove nano3P adapters, then count the number of actual basecalled A's at the 3' end of standards!!

Ideally, this code should work from base outputs (ie. Porechop trimmed fastq and
tailfindr CSVs), b/c this will afford us the most ability to recreate these results in
the future. It'll be a bit annoying b/c this will cause a bit of a slowdown due to
re-running the same crap over and over, but hopefully it's worthwhile!
- Maybe even start with raw fastqs from guppy, then do the Porechop step from there??
    - No, no, this would be a sh!tshow b/c I have to change files in the porechop directory
      in order to switch between dRNA and Nano3P
"""

import mappy as mp
import numpy as np
import pandas as pd

import seaborn as sea
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

import sys
import os

sys.path.insert(0, '/data16/marcus/scripts/nanoporePipelineScripts')
from nanoporePipelineCommon import *

sys.path.insert(0, '/data16/marcus/scripts/nanoporePipelineScripts/standardsAlignment')
import version2_mappingStandardsMethod as std_align

pd.set_option('display.width', 200)
pd.set_option('display.max_columns', None)

class StdLibrary:
    def __init__(self,
                 output_dir: str,
                 is_dRNA_not_Nano3P: bool,
                 ref_fasta: str = "/data16/marcus/scripts/nanoporePipelineScripts"
                                  "/standardsAlignment/220902_version2.0_releventSequences_wOutTails.fasta",
                 threads: int = 20,
                 ):
        if os.path.isfile(f"{output_dir}/cat_files/cat.trimmed.fastq"):
            self.cat_fastq_path = f"{output_dir}/cat_files/cat.trimmed.fastq"
            self.trimmed = True
        elif os.path.isfile(f"{output_dir}/cat_files/cat.fastq"):
            self.cat_fastq_path = f"{output_dir}/cat_files/cat.fastq"
            self.trimmed = False
        else:
            raise FileNotFoundError(f"Could not find a suitable fastq file in the output_dir!")
        
        self.tail_path = find_newest_matching_file(f"{output_dir}/tailfindr/*_tails.csv")
        self.tail_df = pd.read_csv(self.tail_path)
        
        self.ref_fasta = ref_fasta
        self.threads = threads
        
        self.is_dRNA = is_dRNA_not_Nano3P
        self.is_Nano3P = not is_dRNA_not_Nano3P
        
        if self.is_dRNA:
            lib_str = 'dRNA'
        else:
            lib_str = 'cDNA'
        
        # Preform standards assignment!
        self.raw_assignment_df = std_align.standards_aligner_v2(self.ref_fasta, self.cat_fastq_path,
                                                                testing_generated_data=False,
                                                                threads_for_aligner=threads,
                                                                library_type=lib_str)
        self.assignment_df = self.raw_assignment_df.copy(deep=True)
        self.assignment_df = self.assignment_df.reset_index(names=['read_id'])
        
        if self.is_dRNA:
            self.assignment_df.sequence = self.assignment_df.sequence.str.replace("U", "T")
        
        def mappy_cols_to_cols(mappy_obj_column):
            """
            The main reason to have this as its own script is to handle the errors
                that come from failed splits (which result from failed maps!)
            """
            mappy_attribute_names = ['q_st', 'q_en', 'strand', 'ctg', 'ctg_len',
                                     'r_st', 'r_en', 'mlen', 'blen', 'mapq',
                                     'tp', 'ts', 'cigar']
            if isinstance(mappy_obj_column, str):
                return_values = mappy_obj_column.split("\t")
            elif isinstance(mappy_obj_column, mp.Alignment):
                return_values = mappy_obj_column.__str__().split("\t")
            else:
                return_values = [None for _ in range(13)]
            if len(return_values) == 13:
                return return_values
            else:
                print(return_values)
        
        tqdm.pandas(desc="Extracting information from Mappy column")
        self.assignment_df[
            ['q_st', 'q_en', 'strand', 'ctg', 'ctg_len', 'r_st', 'r_en', 'mlen', 'blen', 'mapq', 'tp', 'ts',
             'cigar']] = self.assignment_df.progress_apply(lambda row: mappy_cols_to_cols(row['mappy_hit_obj']), axis=1,
                                                           result_type='expand')
        
        # TODO: There is some weirdness here regarding why not all reads appear in both dataframes!! Figure it out!
        self.merge_df = pd.merge(self.assignment_df, self.tail_df, on="read_id", how="inner")
        self.hit_df = self.merge_df[self.merge_df.assignment.isin(['00', '05', '10', '15', '30', '60'])]
        self.miss_df = self.merge_df[~self.merge_df.assignment.isin(['00', '05', '10', '15', '30', '60'])]
        self.hit_df['cigar'] = self.hit_df.cigar.str.rsplit(":", 1).str[1]
        self.hit_df = self.hit_df.astype({"r_st": int, "r_en": int, "q_st": int, "q_en": int})
        
        # TODO: It would be nice to have this not be hardcoded!
        barcode_dict = {'00': 'GGTGTTGTT',
                        '05': 'CGGCAATAA',
                        '10': 'TAATCGTCC',
                        '15': 'CCTTCTAGG',
                        '30': 'ACACACACC',
                        '60': 'AAGAGGAGG'}
        if self.is_Nano3P:
            # We'll flip (rev. comp.) the barcode dict if we are working with nano3P data (cDNA sequences)
            barcode_dict = {key: mp.revcomp(value) for key, value in barcode_dict.items()}
        
        def check_for_perfect_matches(_barcode_dict, **row):
            assignment = row['assignment']
            seq = row['sequence']
            barcode = _barcode_dict[assignment]
            perfect_barcode = barcode in seq
            return perfect_barcode
        
        self.hit_df['perfect_barcode'] = self.hit_df.apply(lambda row: check_for_perfect_matches(barcode_dict, **row),
                                                           axis=1)
        
        self.stds_ref_dict = {}
        for ref_id, sequence, _, comments in mp.fastx_read(self.ref_fasta, read_comment=True):
            cutdown_ref_id = ref_id[-6:-4]
            self.stds_ref_dict[cutdown_ref_id] = sequence
    
    def plot_perfect_barcodes(self, renderer=None, save_to=None):
        hit_df_groupby = self.hit_df.groupby("assignment")

        grouped_df = hit_df_groupby['assignment'].count().to_frame(name="total_count")
        grouped_df['perfect_match_count'] = hit_df_groupby['perfect_barcode'].sum().to_frame(name="perfect_match_count")
        grouped_df['mean_r_en'] = hit_df_groupby['r_en'].mean().to_frame(name="mean_r_en")

        grouped_df['imperfect_match_count'] = grouped_df.total_count - grouped_df.perfect_match_count
        grouped_df = grouped_df.reset_index()
        fig = px.bar(grouped_df,
                     x='assignment',
                     y=['perfect_match_count', 'imperfect_match_count'],
                     template="plotly_white")
        fig.update_layout(height=500, width=700)
        fig.update_layout(legend=dict(
            orientation='h',
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01
        ))
        if isinstance(renderer, str):
            fig.show(renderer=renderer)
        else:
            fig.show()
        if isinstance(save_to, str):
            if save_to.endswith("html"):
                fig.write_html(save_to)
            else:
                fig.write_image(save_to)
    
    def print_alignments(self,
                         number_to_print: int,
                         select_assignment=None,
                         constrict_to_width=None,
                         select_with_head=False,
                         additional_filtering_query=None):
        def print_mappy_hit_alignment_for_stds(read_id, cigar, r_st, r_en, q_st, q_en, strand, sequence, assignment,
                                               ref_dict=None, line_print_width=None,
                                               **other_row_items) -> None:
            import re
            ref_seq = ref_dict[assignment].upper()
            print(f"\nread_id={read_id}; assignment={assignment}")
            parsed_cigar = re.findall(rf'(\d+)([MDNSIX])', cigar)
            parsed_cigar = [(int(num), char) for num, char in parsed_cigar]
            ref_seq = ref_seq[r_st: r_en]
            ref_pos = 0
            sequence = sequence[q_st: q_en]
            read_pos = 0
            if strand == "-":
                sequence = mp.revcomp(sequence)
            top_line = ""
            middle_line = ""
            bottom_line = ""
            for length, code in parsed_cigar:
                if code == "M":  # Map (Read & Ref Match)
                    read_map_piece = sequence[read_pos:read_pos + length]
                    ref_map_piece = ref_seq[ref_pos:ref_pos + length]
                    perfect_matches = ""
                    for index, char in enumerate(read_map_piece):
                        try:
                            if char == ref_map_piece[index]:
                                perfect_matches += "|"
                            else:
                                perfect_matches += "•"
                        except IndexError:
                            perfect_matches += " "
                    top_line += read_map_piece
                    middle_line += perfect_matches
                    bottom_line += ref_map_piece
                    ref_pos += length
                    read_pos += length
                elif code == "I":  # Insert (Gap in Ref)
                    top_line += sequence[read_pos:read_pos + length]
                    middle_line += " " * length
                    bottom_line += " " * length
                    read_pos += length
                elif code == "D" or code == "N":  # Delete (Gap in Read)
                    top_line += " " * length
                    middle_line += " " * length
                    bottom_line += ref_seq[ref_pos:ref_pos + length]
                    ref_pos += length
            if isinstance(line_print_width, int):
                num_blocks = int(np.ceil(len(top_line) / line_print_width))
                print_blocks = []
                for block_index in range(num_blocks):
                    print_blocks.append([
                        top_line[block_index * line_print_width:(block_index + 1) * line_print_width],
                        middle_line[block_index * line_print_width:(block_index + 1) * line_print_width],
                        bottom_line[block_index * line_print_width:(block_index + 1) * line_print_width],
                    ])
                for top, mid, bot in print_blocks:
                    print()
                    print(f"Read: {top}")
                    print(f"      {mid}")
                    print(f"Ref:  {bot}")
            else:
                print(f"Read: {top_line}", f"      {middle_line}", f"Ref:  {bottom_line}", sep='\n')

        if isinstance(select_assignment, str):
            print_df = self.hit_df[self.hit_df.assignment == select_assignment]
        else:
            print_df = self.hit_df

        if isinstance(additional_filtering_query, str):
            print_df = print_df.query(additional_filtering_query)

        if select_with_head:
            print_df = print_df.head(number_to_print)
        else:
            print_df = print_df.sample(number_to_print)

        if isinstance(constrict_to_width, int):
            print_df.apply(lambda row: print_mappy_hit_alignment_for_stds(ref_dict=self.stds_ref_dict,
                                                                          line_print_width=constrict_to_width, **row),
                           axis=1)
        else:
            print_df.apply(
                lambda row: print_mappy_hit_alignment_for_stds(ref_dict=self.stds_ref_dict, line_print_width=175,
                                                               **row), axis=1)

if __name__ == '__main__':
    SY_nano3P = StdLibrary(f"/data16/marcus/nanoporeSoftLinks/"
                           f"230224_nanoporeRun_RNAStds_SY-TGIRT-50ng_Nano3P/output_dir", False)
    dRNA = StdLibrary(f"/data16/marcus/nanoporeSoftLinks/"
                      f"221112_nanoporeRun_ENO2RNAStds_dRNA/output_dir", True)
    # BHLT_nano3P = StdLibrary(f"/data16/marcus/working/230501_nanoporeRun_PureRNAStds_Nano3P/output_dir", False)
    print("done.")

