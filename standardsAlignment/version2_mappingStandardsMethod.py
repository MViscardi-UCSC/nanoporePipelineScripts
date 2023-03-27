"""
version2_mappingStandardsMethod.py
Marcus Viscardi,    September 21, 2022

Taking things from generatingFakeReadsForTesting.ipynb and pulling them into a single method.

Also trying to make lots of comments so that this will be easier to understand in the future.
    The version from 2021 was extremely hard to understand, and required a complete rewrite...
    Let's avoid that!

Generally, the new method consists of a few primary steps:
    1. Take reads and align them to a mini-genome that only has the 6 standards (w/out tails)
        which sole difference is their barcodes (8bp); because of this lack of difference reads
        that map to one will likely map to them all! Since we end up mapping to all 6 reference
        standards we need to be deliberate about how we pick up the "correct" alignment.
    2. We can first attempt to pick the correct barcode by looking for single standard that has
        the highest mlen score, which corresponds to the most matched bases. Usually when the
        barcode sequence has not been hit by any indels this will immediately pull out the
        correct standard. Sometimes multiple standards will have tied for the highest mlen,
        in these cases we'll continue down the process w/ step 3.
    3. If we have multiple standards tied for the number of matched bases with the read, it likely
        indicates that some indel has been introduced to the barcode making it harder for mappy to
        easily pick the right standard. In this case we can start to look at the NM score, which
        corresponds to the number of mismatches across the alignment. The alignment with highest
        mlen score and the lowest number of mismatches will correspond to an accurate assignment
        (almost always), because of the large levenshtein distances between barcodes.
    4. Finally, if there are multiple alignments that are tied for both the highest mlen score and
        the lowest NM score, we give up and call these ambiguous. The most common cause for these
        kinds of failures is a complete loss of the unique barcode region, resulting in a read that
        maps the exact same way to all 6 reference standards. The important thing about this system
        is that we can throw out these reads, if we allow an aligner like mappy to run with default
        settings, it would pick a (seemingly random) best standard and call that the primary
        alignment, than assign the others as secondaries. This would produce false positives that could
        tank our downstream analysis
"""
import mappy
import mappy as mp
from pprint import pprint
import pandas as pd
import numpy as np
from tqdm import tqdm
from typing import Union, List, Tuple, Any
import random
from nanoporePipelineCommon import find_newest_matching_file, pick_libs_return_paths_dict, get_dt
import textwrap

import seaborn as sea
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

pd.set_option("display.max_columns", None)


def print_mappy_hit_alignment(mappy_hit_obj: mappy.Alignment,
                              read_seq: str, ref_seq: str,
                              line_print_width=None,
                              do_not_print=False) -> Tuple[str, str, str]:
    import re
    parsed_cigar = re.findall(rf'(\d+)([MDNSIX])', mappy_hit_obj.cigar_str)
    parsed_cigar = [(int(num), char) for num, char in parsed_cigar]
    ref_seq = ref_seq[mappy_hit_obj.r_st: mappy_hit_obj.r_en].upper()
    ref_pos = 0
    read_seq = read_seq[mappy_hit_obj.q_st: mappy_hit_obj.q_en].upper()
    read_pos = 0
    if mappy_hit_obj.strand == -1:
        read_seq = mp.revcomp(read_seq)
    top_line = ""
    middle_line = ""
    bottom_line = ""
    for length, code in parsed_cigar:
        if code == "M":  # Map (Read & Ref Match)
            read_map_piece = read_seq[read_pos:read_pos + length]
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
            top_line += read_seq[read_pos:read_pos + length]
            middle_line += " " * length
            bottom_line += " " * length
            read_pos += length
        elif code == "D" or code == "N":  # Delete (Gap in Read)
            top_line += " " * length
            middle_line += " " * length
            bottom_line += ref_seq[ref_pos:ref_pos + length]
            ref_pos += length
    if do_not_print:
        return top_line, middle_line, bottom_line
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
            print("\n\n")
            print(f"Read: {top}")
            print(f"      {mid}")
            print(f"Ref:  {bot}")
        return top_line, middle_line, bottom_line
    else:
        print(top_line, middle_line, bottom_line, sep='\n')
        return top_line, middle_line, bottom_line


def apply_extract_mappy_hit_obj(**columns) -> mp.Alignment or None:
    # TODO: add the rest of the actual extracting... lol
    assigned_std = columns['assignment']
    correct_hit_obj = f"hit_obj_{assigned_std}"
    if assigned_std in ["00", "05", "10", "15", "30", "60"]:
        return columns[correct_hit_obj]
    else:
        return None


def per_read_mapping(read_id: str, sequence: str, aligner: mp.Aligner,
                     print_weird_results: bool = False) -> dict:
    per_read_dict = {  # This will be a dictionary to hold more general information about the read & mapping
        'sequence': sequence
    }

    # For each sequence, create empty dicts to store info about mapping hits
    hit_objs = {}
    # Loop through each of the mapping hits generated:
    for hit in aligner.map(sequence):
        hit_objs[hit.ctg] = hit  # and store each one

    if hit_objs:  # First check if there are any hits at all
        for aligned_std, hit_obj in hit_objs.items():
            aligned_std_int = aligned_std[-6:-4]
            per_read_dict[f"hit_obj_{aligned_std_int}"] = hit_obj
        mlen_dict = {key[-2:]: hit_obj.mlen for
                     key, hit_obj in per_read_dict.items()
                     if key.startswith("hit_obj")}
        mlen_list = list(mlen_dict.values())
        list_of_max_mlen = [std for std in mlen_dict if mlen_dict[std] == max(mlen_list)]
        if mlen_list.count(max(mlen_list)) == 1:  # There is only one std with a maximized map length
            assignment = list_of_max_mlen[0]
        elif mlen_list.count(max(mlen_list)) == 6:  # All stds are mapped, prob b/c the read is truncated
            assignment = "Ambiguous_All"
        else:  # >1 std mapped with equal length, we can check NM for more info
            NM_dict = {key[-2:]: hit_obj.NM for
                       key, hit_obj in per_read_dict.items()
                       if key.startswith("hit_obj") and key[-2:] in list_of_max_mlen}
            NM_list = list(NM_dict.values())
            if NM_list.count(min(NM_list)) == 1:  # There is one std that has the lowest mismatch count
                assignment = [std for std in NM_dict if NM_dict[std] == min(NM_list)][0]
            else:  # Not sure what these situations will look like! They should have been caught in elif above!
                str_of_tied_stds = "_".join(sorted(list(NM_dict.keys())))
                assignment = f"Ambiguous_subset_{str_of_tied_stds}"
                if print_weird_results:
                    print(f"\nWEIRD RESULTS for read: {read_id}; {assignment}")
                    for standard in sorted(list(NM_dict.keys())):
                        print(f"\nMap for: ENO2_finalStandard_{standard}Tail")
                        print_mappy_hit_alignment(per_read_dict[f'hit_obj_{standard}'],
                                                  sequence,
                                                  aligner.seq(f'ENO2_finalStandard_{standard}Tail'))
        per_read_dict['assignment'] = assignment
        return per_read_dict
    else:  # If we are here then nothing mapped!
        per_read_dict['assignment'] = 'Failed_To_Map'
        return per_read_dict


def align_stds_from_fastx(fastx: str, mappy_aligner: mp.Aligner,
                          print_weird_results: bool = False):
    if fastx.endswith("q"):
        lines_per_entry = 4
    elif fastx.endswith("a"):
        lines_per_entry = 2
    else:
        raise NotImplementedError(f"The input fastx file must end in fastq/fq or fasta/fa!!")
    read_iterator = tqdm(mp.fastx_read(fastx, read_comment=True),
                         # The Mappy library (from Minimap2) has a fastx file reader that is very fast
                         total=sum(1 for line in open(fastx)) // lines_per_entry,
                         # The mappy.fastx_read iterator doesn't have a __len__, so we
                         #    have to manually tell tqdm how long the fastq will be
                         desc=f"Reading fastx entries and assigning to standards"
                         )
    storage_dict = {}
    for read_id, sequence, _, comments in read_iterator:
        storage_dict[read_id] = per_read_mapping(read_id, sequence, mappy_aligner,
                                                 print_weird_results=print_weird_results)
    return storage_dict


def align_stds_from_merge_df(merge_df,
                             mappy_aligner: mp.Aligner,
                             print_weird_results: bool = False):
    """
    For this to work well in the pipeline, all I really need is to
    extract the reads that mapped initially to the ENO2 chromosome,
    then get the assignments from those. Finally, I'll return a very
    simple dataframe with just these columns:
        read_id
        standard_assignment
    This might even just be a series?!
    I can decide if we want more information later!
    
    :param merge_df: 
    :param mappy_aligner: 
    :param print_weird_results: 
    :return: 
    """
    # First lets just collapse the df into the columns we care about:
    columns_to_keep = ['read_id',
                       'chr_id',
                       'sequence']
    merge_df = merge_df[columns_to_keep]
    
    

def standards_aligner_v2(path_to_standards_ref_fasta: str,
                         fastx_file: str = None,
                         sam_df: pd.DataFrame = None,
                         mjv_compressed_df: pd.DataFrame = None,
                         threads_for_aligner: int = 10,
                         print_weird_results: bool = False,
                         testing_generated_data: bool = False,
                         library_type: str = "dRNA",
                         ) -> pd.DataFrame:
    if library_type == "dRNA":
        mappy_preset = "map-ont"
        extra_mappy_flag = 0x100000
    elif library_type == "cDNA":
        mappy_preset = "map-ont"
        extra_mappy_flag = 0x200000
    else:
        raise NotImplementedError(f"Please provide 'dRNA' or 'cDNA' for library_type, you provided: {library_type}")
    aligner = mp.Aligner(path_to_standards_ref_fasta,
                         preset=mappy_preset, k=14,  # These are the usual go-tos for mapping ONT dRNA-Seq
                         extra_flags=extra_mappy_flag,  # From: https://github.com/lh3/minimap2/blob/master/minimap.h
                         #           0x200000 = forces strand reverse strand alignment (b/c we have cDNA)
                         n_threads=threads_for_aligner,
                         )
    if not aligner:
        # This was in the tutorial, haven't seen this yet.
        raise Exception("ERROR: failed to load/build index")

    # Print the contigs (chromosomes) which in this case are the adapters in the fastq
    print(f"Mapping to adapters named: {aligner.seq_names}")

    if isinstance(fastx_file, str):
        # align_stds_from_fastx() will do the majority of the work here, see the method for details
        storage_dict = align_stds_from_fastx(fastx_file, aligner,
                                             print_weird_results=print_weird_results)

        # We convert the dictionary output from standards assignment to a dataframe,
        #   as this makes it a little easier to finesse.
        temp_df = pd.DataFrame.from_dict(storage_dict, orient='index')
        columns_to_keep = ['assignment',
                           'sequence']
        # At this point the read_ids are the indexes
        df = temp_df[columns_to_keep]

        # Here we can pull more information about the winning alignment to store in the dataframe if needed:
        tqdm.pandas(desc=f"Identifying assignments from Mappy objects")
        pd.options.mode.chained_assignment = None
        df['mappy_hit_obj'] = temp_df.progress_apply(lambda row: apply_extract_mappy_hit_obj(**row), axis=1)

        # This is a holdover from when I was testing with generated data
        if testing_generated_data:
            df['original_std'] = df.index.str[6:8]
            df['correct_hit'] = df.assignment == df.original_std
            pprint(df.correct_hit.value_counts())
            pd.options.mode.chained_assignment = 'warn'
    elif isinstance(mjv_compressed_df, pd.DataFrame):
        df = mjv_compressed_df
        # First check that the cerENO2 chr_id shows up in the compressed on reads dataframe:
        if "cerENO2" not in mjv_compressed_df['chr_id'].unique().tolist():
            return f"Chromosome cerENO2 not found in dataframe, please run minimap2 w/ genome that has cerENO2 added!!"
        # Todo: Write this section up. It would be nice to have this work with the compressed dataframes I'll actually
        #       be using for all this in the future!
    else:
        df = None
        raise NotImplementedError("Please provide either a fastq or a df")
    # pprint(df.value_counts('assignment'))
    if testing_generated_data:
        pprint(df[df['assignment'].isin(['00', '05', '10', '15', '30', '60'])]
               .groupby(['assignment',
                         'original_std'])
               .size().unstack(fill_value=0))
        print("\n")
        pprint(df
               .groupby(['assignment',
                         'original_std'])
               .size().unstack(fill_value=0))
    return df


def sequence_subindel_creater(seq,
                              insert_prob=0.01,
                              delete_prob=0.01,
                              sub_prob=0.05,
                              indel_size_factor=1.0,
                              upseq_truncate_prob=0.001,
                              downseq_truncate_prob=0.001) -> Tuple[Any, Tuple[Any, Any, Any]]:
    # TODO: Add ability to output multiple sequences
    seq = seq.upper()
    adjusted_seq = ""
    choices_per_nt = {"insert": insert_prob,
                      "delete": delete_prob,
                      "substitute": sub_prob,
                      "upseq_trunc": upseq_truncate_prob,
                      "downseq_trunc": downseq_truncate_prob,
                      "nothing": 1 - (insert_prob + delete_prob + sub_prob)}
    indel_sizes = [i + 1 for i in range(25)]
    indel_weights = [10 ** ((indel_size_factor * -i) + 1) for i in range(25)]
    nucleotide_choices = ['G', 'C', 'T', 'A']

    has_been_truncated = False
    trunc_pos = 0
    trunc_dir = "NA"

    index = 0
    while index < len(seq):
        nucl = seq[index]  # by putting this up here, we don't have to recall it below
        if nucl == "N":
            nucl = random.choice(nucleotide_choices)  # This will throw a random nucl instead of the N!
        choice = random.choices(list(choices_per_nt.keys()), weights=list(choices_per_nt.values()))[0]
        if choice == "insert":
            index += 1
            insert_size = random.choices(indel_sizes, weights=indel_weights)[0]
            for _ in range(insert_size):
                nucl += random.choice(nucleotide_choices)
        elif choice == "substitute":
            index += 1
            nucl = random.choice(nucleotide_choices)
        elif choice == 'delete':
            delete_size = random.choices(indel_sizes, weights=indel_weights)[0]
            index += delete_size
            nucl = ''
        elif choice == 'upseq_trunc':
            if not has_been_truncated:
                index += 1
                adjusted_seq = ''
                has_been_truncated = True
                trunc_pos = index
                trunc_dir = "up"
            else:  # This prevents multiple truncations from triggering
                index += 1
        elif choice == 'downseq_trunc':
            if not has_been_truncated:
                has_been_truncated = True
                trunc_pos = index
                trunc_dir = "down"
                break  # Break out of the loop, not adding any more nucleotides!
            else:  # This prevents multiple truncations from triggering
                index += 1
        else:
            index += 1
        adjusted_seq += nucl
    return adjusted_seq, (has_been_truncated, trunc_pos, trunc_dir)


def generate_test_std_reads(reference_fasta="/data16/marcus/scripts/nanoporePipelineScripts/standardsAlignment"
                                            "/220902_version2.0_releventSequences.fasta",
                            standards_to_use=(0, 5, 10, 15, 30, 60),
                            number_of_entries_in_fasta=1000,
                            insert_prob_quick=0.01,
                            delete_prob_quick=0.01,
                            sub_prob_quick=0.05,
                            upseq_trunc_prob_quick=0.001,
                            downseq_trunc_prob_quick=0.001,
                            indel_size_factor_quick=1.0,  # Inversely related to the rate of larger and larger indels!
                            output_file_suffix="generatedRNAStdsReads_withTruncations.fasta",
                            output_file_dir=".",
                            output_file_overwrite_full_path=None,
                            save_output=True,
                            ) -> str:
    if isinstance(output_file_overwrite_full_path, str):
        output_file_name = output_file_overwrite_full_path
    else:
        output_file_name = f"{output_file_dir}/{get_dt(for_file=True)}_{output_file_suffix}"

    finalStandards_dict = {}
    for read_id, sequence, quality, comments in mp.fastx_read(reference_fasta, read_comment=True):
        if "finalStandard" in read_id:
            finalStandards_dict[read_id] = (sequence, comments)
    print(f"Generating reads produced by: {list(finalStandards_dict.keys())}")
    print(f"Using parameters:\n"
          f"\tInsert rate={insert_prob_quick}\n"
          f"\tDelete rate={delete_prob_quick}\n"
          f"\tSubstitute rate={sub_prob_quick}\n"
          f"\tTruncation rates: up={upseq_trunc_prob_quick}; down={downseq_trunc_prob_quick}\n"
          f"\tIndel size weights:\n\t\t"
          f"{dict(zip([i + 1 for i in range(25)], [10 ** ((indel_size_factor_quick * -i) + 1) for i in range(25)]))}"
          f"\n\n")
    fasta_storage_dict = {}
    truncated_count = 0
    iterator = tqdm(range(number_of_entries_in_fasta))
    with open(output_file_name, 'w') as output_fasta:
        for i in iterator:
            template_standard_size = random.choice(standards_to_use)
            template_standard = f"ENO2_finalStandard_{template_standard_size:0>2}Tail"
            entry_num = f"{i + 1:0>{1 + int(np.ceil(np.log10(number_of_entries_in_fasta)))}}"
            generated_seq, trunc_info = sequence_subindel_creater(finalStandards_dict[template_standard][0],
                                                                  insert_prob=insert_prob_quick,
                                                                  delete_prob=delete_prob_quick,
                                                                  sub_prob=sub_prob_quick,
                                                                  indel_size_factor=indel_size_factor_quick,
                                                                  upseq_truncate_prob=upseq_trunc_prob_quick,
                                                                  downseq_truncate_prob=downseq_trunc_prob_quick)[:]
            truncated_count += trunc_info[0]
            seq_id = f"RNAStd{template_standard_size:0>2}.N{entry_num} template={template_standard} " \
                     f"entry_number={entry_num} truncated={trunc_info[0]} truncated_pos={trunc_info[1]} " \
                     f"truncated_dir={trunc_info[2]}"
            fasta_storage_dict[seq_id] = generated_seq
            file_line = f">{seq_id}\n{textwrap.fill(generated_seq, width=5000)}\n"
            if save_output:
                output_fasta.write(file_line)
            else:
                print(file_line)
    # print(truncated_count)
    return output_file_name


def plot_par_categ(df,
                   title=None):
    df['seq_len'] = df.sequence.str.len()
    df['original_std'] = df.original_std.astype(int)
    # df['assignment'] = df.assignment.astype(int)
    df['false_positive'] = df.assignment.isin(['00',
                                               '05',
                                               '10',
                                               '15',
                                               '30',
                                               '60']) & ~df.correct_hit
    fig = px.parallel_categories(df,
                                 template='plotly_white',
                                 color='original_std',
                                 dimensions=['assignment',
                                             'original_std',
                                             'correct_hit',
                                             # 'seq_len',
                                             'false_positive',
                                             ],
                                 title=title)
    fig.show(renderer='firefox')

def check_for_perfect_matches(rc_barcode_dict, **row):
    assignment = row['assignment']
    sequence = row['sequence']
    barcode = rc_barcode_dict[assignment]
    perfect_barcode = barcode in sequence
    return perfect_barcode


if __name__ == '__main__':
    title_text = "Testing presence of extra barcode in<br>generated sequences"
    # ref_fasta = "220902_version2.0_releventSequences.fasta"
    ref_fasta = "220902_version2.0_releventSequences_wOutTails.fasta"
    # ref_fasta = "/data16/marcus/genomes/plus_cerENO2_elegansRelease100/220923_allChrs_plus-cerENO2.allChrs.fa"
    gen_ref_fasta = "220902_version2.0_releventSequences_wOutTails_wExtraBarcode.fasta"
    # fasta = f"220920_generatedRNAStdsReads_withTruncations.fasta"
    
    # Nano3P Post Porechop trimming
    # fastq = "/data16/marcus/working/221028_nanoporeRun_ENO2RNAStds_Nano3P/output_dir/cat_files/cat.trimmed.fastq"
    
    # Seq Library
    lib_type = 'cDNA'  # or 'RNA'
    output_dir = f"/data16/marcus/working/230224_nanoporeRun_RNAStds_SY-TGIRT-50ng_Nano3P/output_dir"
    fastq = f"{output_dir}/cat_files/cat.fastq"
    assigned_df = standards_aligner_v2(ref_fasta, fastq,
                                       testing_generated_data=False,
                                       threads_for_aligner=20,
                                       library_type=lib_type)
    pprint(assigned_df)
    pprint(assigned_df.assignment.value_counts(normalize=True))
    fig = px.bar(assigned_df.assignment.value_counts(normalize=False))
    fig.show(renderer='firefox')
    assigned_df.to_csv(f"{get_dt(for_file=True)}_RNAStds_dRNA_Assignments.csv")

    barcode_dict = {
        '00': 'GGTGTTGTT',
        '05': 'CGGCAATAA',
        '10': 'TAATCGTCC',
        '15': 'CCTTCTAGG',
        '30': 'ACACACACC',
        '60': 'AAGAGGAGG',
    }
    
    rc_barcode_dict = {key: mp.revcomp(value) for key, value in barcode_dict.items()}

    hit_df = assigned_df[assigned_df.assignment.isin(['00', '05', '10', '15', '30', '60'])]
    
    if lib_type == 'cDNA':
        barcode_dict_to_check = rc_barcode_dict
    elif lib_type == 'RNA':
        barcode_dict_to_check = barcode_dict
    else:
        raise NotImplementedError
    # It's important to use the barcode_dict or the rc_barcode_dict depending on if we are using cDNA reads or RNA reads
    hit_df['perfect_barcode'] = hit_df.apply(lambda row: check_for_perfect_matches(barcode_dict_to_check,
                                                                                   **row), axis=1)

    hit_df_groupby = hit_df.groupby("assignment")

    grouped_df = hit_df_groupby['assignment'].count().to_frame(name="total_count")
    grouped_df['perfect_match_count'] = hit_df_groupby['perfect_barcode'].sum().to_frame(name="perfect_match_count")
    # grouped_df['mean_r_en'] = hit_df_groupby['r_en'].mean().to_frame(name="mean_r_en")

    grouped_df['imperfect_match_count'] = grouped_df.total_count - grouped_df.perfect_match_count
    grouped_df = grouped_df.reset_index()
    
    fig = px.bar(grouped_df, x='assignment',
                 y=['perfect_match_count',
                    'imperfect_match_count',
                    ])
    fig.update_layout(height=500, width=700)
    fig.update_layout(legend=dict(
        orientation='h',
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=0.01
    ))
    save_path = f"{output_dir}/{get_dt()}_barcodeHits"
    fig.write_html(f"{save_path}.html")
    fig.write_image(f"{save_path}.svg")
    fig.write_image(f"{save_path}.png")
    fig.show()
