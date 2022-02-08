"""
standardsAssignmentWithMinimap2.py
Marcus Viscardi,    October 07, 2021

So since cutadapt didn't seem to work very well, I am
    circling back to minimap2. I am hoping that I will
    be able to use the mapper with some stringent params
    in order to "map" and assign the correct adapters!

The general plan will be to take the reads that already
    mapped to the pTRI chromosome (the one that lacked
    any of the UMI information), then try to map those
    to versions of the same that include the UMIs.

Another potential option here that I spaced on before:
    I could use the mapping information from the first
    minimap2 pass to identify where in each read the
    UMI should be! There might be some wiggle room
    needed due to the prevalence of indels w/ ONT, but
    this seems like a feasible way to either strengthen
    my other tools, or just do this whole thing in a
    much simpler manner!
"""
# Mappy?! https://www.biostars.org/p/395092/
import mappy as mp
from pprint import pprint
import pandas as pd
from tqdm import tqdm
from typing import Union
from nanoporePipelineCommon import find_newest_matching_file, pick_libs_return_paths_dict, get_dt

pd.set_option("display.max_columns", None)


def align_standards(fastq_file=None, compressed_df=None, keep_read_id=False, bar_width=None,
                    keep_minimap_obj=False, threads=10, tail_less_reference=False,
                    **kwargs) -> Union[str, pd.DataFrame]:
    if tail_less_reference:  # This ref doesn't have the tails or HDV attached to the 3' ends
        path_to_genome = "/data16/marcus/genomes/plus-pTRIxef_elegansRelease100/" \
                         "intermediate_files/standards_with_indexes_without_tails_or_HDV.fa"
    else:  # This ref does...
        path_to_genome = "/data16/marcus/genomes/plus-pTRIxef_elegansRelease100/" \
                         "intermediate_files/standards_with_indexes.fa"
    aligner = mp.Aligner(path_to_genome,
                         preset="splice", k=14,  # These are the usual go-tos for mapping ONT dRNA-Seq
                         extra_flags=0x100000,  # From: https://github.com/lh3/minimap2/blob/master/minimap.h
                         #           0x100000 = forces strand matched alignment (b/c dRNA-seq retains strand info)
                         n_threads=threads,
                         )
    if not aligner:
        raise Exception("ERROR: failed to load/build index")  # This was in the tutorial, haven't seen this yet.

    # Print the contigs (chromosomes) which in this case are the adapters in the fastq
    print(f"Mapping to adapters named: {aligner.seq_names}")

    seq_assignments = []  # This list will hold all the alignments eventually

    # This method optionally takes a fastq file OR a compressed dataframe.
    #   The method will throw an error if neither are given.

    # If a fastq file is passed:
    if isinstance(fastq_file, str):
        # Create a tqdm iterator that is going to loop through the fastq,
        #   and provide a progress bar along the way!
        read_iterator = tqdm(mp.fastx_read(fastq_file),
                             # The Mappy library (from Minimap2) has a fastx file reader that is very fast
                             total=sum(1 for line in open(fastq_file)) // 4,
                             # The mappy.fastx_read iterator doesn't have a __len__, so we
                             #    have to manually tell tqdm how long the fastq will be
                             ncols=bar_width)
        for read_id, sequence, _ in read_iterator:
            seq_assignment = _loop_align_seq_to_adapters(aligner,
                                                         read_id,
                                                         sequence,
                                                         keep_read_id=keep_read_id,
                                                         was_fastq=True,
                                                         print_stuff=False,
                                                         print_per_seq=False,
                                                         print_per_seq_plus=False, )
            # Save the assignment dictionary from the aligner loop into the list:
            seq_assignments.append(seq_assignment)
            # Set the progress bar description, mostly cuz it's fun...
            read_iterator.set_description(f"Processing {read_id}")

    # if a compressed on reads dataframe is passed:
    elif isinstance(compressed_df, pd.DataFrame):
        # First check that the pTRI chr_id shows up in the compressed on reads dataframe:
        if "pTRI" not in compressed_df['chr_id'].unique().tolist():
            return f"Chromosome pTRI not found in dataframe, please run minimap2 w/ genome that has pTRI added!!"
        # Build a tqdm iterator to loop through the compressed on reads dataframe.
        #   We'll only give it the three columns that we care about, making the whole datastructure lighter
        row_iterator = tqdm(compressed_df[['read_id', 'sequence', 'chr_id']].to_dict(orient='records'),
                            ncols=bar_width)
        for row_dict in row_iterator:
            seq_assignment = _loop_align_seq_to_adapters(aligner,
                                                         keep_read_id=keep_read_id,
                                                         print_stuff=False,
                                                         print_per_seq=False,
                                                         print_per_seq_plus=False,
                                                         **row_dict,
                                                         # The **row_dict allows us to just pass the columns
                                                         #    of interest w/out needing to unpack them.
                                                         )
            # Save the assignment dictionary from the aligner loop into the list:
            seq_assignments.append(seq_assignment)
            # Set the progress bar description, mostly cuz it's fun...
            row_iterator.set_description(f"Chr: {row_dict['chr_id']:>5}; Read: {row_dict['read_id']}")
    else:
        raise NotImplementedError("Please provide either a fastq or a df")

    print("Finished Alignment")

    mappy_df = pd.DataFrame(seq_assignments)

    if not keep_minimap_obj:
        # These columns will be extracted from the mapping object:
        mappy_cols = ["stds_q_start",
                      "stds_q_end",
                      "stds_strand",
                      "stds_chr_id",
                      "stds_chr_len",
                      "stds_r_start",
                      "stds_r_end",
                      "stds_mlen",
                      "stds_blen",
                      "stds_mapq",
                      "stds_type_of_alignment",
                      "stds_ts",
                      "stds_cigar",
                      ]

        # Loop to extract and drop the mappy python object
        mappy_df[mappy_cols] = mappy_df.apply(lambda row: _loop_pull_mappy_obj_info(row["mapping_obj"],
                                                                                    columns=mappy_cols),
                                              axis=1, result_type="expand")
        mappy_df.drop(columns="mapping_obj", inplace=True)

        # Extract the last bit of information stored in the mappy flags
        for mappy_flag_col in ["stds_type_of_alignment", "stds_ts", "stds_cigar"]:
            mappy_df[mappy_flag_col] = mappy_df[mappy_flag_col].str.split(":").str[-1]
    return mappy_df


def _loop_align_seq_to_adapters(aligner, read_id, sequence, chr_id=None,
                                print_per_seq=False, print_per_seq_plus=False,
                                keep_read_id=False, was_fastq=False,
                                print_stuff=False, **kwargs):
    # This will quickly filter out those reads that didn't initially map to the pTRI chr
    if not was_fastq and chr_id != "pTRI":
        seq_assignment = {"adapter": "not_pTRI_standard", "mapping_obj": None}
        if keep_read_id:
            seq_assignment["read_id"] = read_id
        return seq_assignment

    if print_per_seq:
        print(f"Aligning read: {read_id}", end="\t")

    # For each sequence, create empty dicts to store info about mapping hits
    hit_objs = {}

    # Loop through each of the mapping hits generated:
    for hit in aligner.map(sequence):
        hit_objs[hit.ctg] = hit  # and store each one

    # If1: there are any mapping hits:
    if hit_objs:

        # If2: all the hits are not the same  M-length! (would be the same if they all
        #     had mapped to non-unique areas <- i.e. not the index region)
        hit_matches = {k: hit.mlen for k, hit in hit_objs.items()}  # count the number of Ms for each map
        if len(set(hit_matches.values())) > 1:
            # Find and store the longest map:
            max_match = max(hit_matches.values())
            # Make a list of hits that match the longest map:
            max_match_keys = [k for k, v in hit_matches.items() if v == max_match]

            # If3: there is only one longest mapping hit:
            if len(max_match_keys) == 1:
                # It wins, store it.
                seq_assignment = {"adapter": max_match_keys[0],
                                  "mapping_obj": hit_objs[max_match_keys[0]]}
                if keep_read_id:
                    seq_assignment["read_id"] = read_id
                if print_per_seq:
                    print(seq_assignment)
                    read_end = seq_assignment["mapping_obj"].q_en
                    print(sequence[read_end - 10:read_end + 2], read_end)

            # Else3: There is more than 1 longer mapping hit.
            #        i.e. Two indexes have the same number of Ms inside
            else:
                # Check if any of the hits map out to the index in chr space:
                if max([hit.r_en for k, hit in hit_objs.items()]) < 1452:
                    # If none due, I assume the read wasn't long enough (BAD ASSUMPTION)
                    seq_assignment = {"adapter": "ambiguous_subset_too_short",
                                      "mapping_obj": None}
                    if print_per_seq:
                        print("Ambiguous, subset of adapters too short")
                # Otherwise, it's likely the example case for If-Else loop #3
                else:
                    index_matches = {}

                    # For each map hit we'll make a string that matches in length w/ the ref chr:
                    for std, hit in hit_objs.items():

                        # Create variables to store info to for each map hit
                        map_seq = ''  # This will hold the gaped read sequence
                        temp_seq = sequence[hit.q_st:]  # Start the sequence at the first mapped base
                        ref = aligner.seq(hit.ctg)  # Start the reference at the first mapped base
                        alignment = ''  # This will hold pipes or spaces to match with M's from the CIGAR

                        # Loop through the cigar string, walking along the sequence and reference along the way:
                        for cig_count, cig_type in hit.cigar:

                            # IF the cigar tag is a map, walk same distance in both seq and ref
                            if cig_type == 0:  # (M)apped
                                map_seq += temp_seq[:cig_count]
                                temp_seq = temp_seq[cig_count:]
                                alignment += '|' * cig_count  # and add matches to th alignment
                            # IF the cigar tag is a deletion or intron, add spacers to the sequence
                            elif cig_type == 2 or cig_type == 3:  # (D)eletion or i(N)tron
                                map_seq += '.' * cig_count
                                alignment += ' ' * cig_count  # and add misses to the alignment
                            # IF the cigar tag is an insertion, only walk along the reference
                            elif cig_type == 1:  # (I)nsertion
                                temp_seq = temp_seq[cig_count:]
                            # Otherwise, catch any cases with other cigar tags (I should have any?)
                            else:
                                raise NotImplementedError(f"Unexpected code in cigar string: [{cig_count}, {cig_type}]")
                        if print_per_seq_plus:
                            print('\n' + hit.ctg + ':\t' + read_id + '\tMaps:' + str(hit.mlen))
                            print(f"  ref_start: {hit.r_st:>3};\t  ref_end: {hit.r_en:>3}\n"
                                  f"query_start: {hit.q_st:>3};\tquery_end: {hit.q_en:>3}")
                            print(ref[hit.r_st:], alignment, map_seq, sep='\n')

                        # I can't decide between keeping or ignoring the flanking CTCT's
                        # W/ CTCT flanking region
                        # umi_region_ref = ref[1452:1452+16]
                        # umi_region_alignment = alignment[1452-hit.r_st:1452-hit.r_st+16]
                        # umi_region_query = map_seq[1452-hit.r_st:1452-hit.r_st+16]

                        # W/out CTCT flanking region
                        umi_region_ref = ref[1456:1456 + 8]
                        umi_region_alignment = alignment[1456 - hit.r_st:1456 - hit.r_st + 8]
                        umi_region_query = map_seq[1456 - hit.r_st:1456 - hit.r_st + 8]

                        # Make a short list of 1's for matched nucleotides and 0's for mismatched:
                        umi_region_matches = [1 if x == y else 0 for x, y in zip(umi_region_ref, umi_region_query)]
                        # Count number of matches (1's):
                        umi_region_match_count = sum(umi_region_matches)
                        # Store that count in the original dictionary from way earlier:
                        index_matches[std] = umi_region_match_count

                        if print_stuff:
                            if len(umi_region_query) >= 2:
                                print(umi_region_ref,
                                      umi_region_alignment,
                                      umi_region_query,
                                      ''.join(map(str, umi_region_matches)),
                                      f"Match sum: {sum(umi_region_matches)}",
                                      sep='\n')
                                print('\n')
                            else:
                                print(umi_region_ref,
                                      f"*Too Short Match",
                                      sep='\n')
                    # Figure out the best UMI match we saw:
                    max_index_match = max(index_matches.values())
                    # Again, make a list of the one or multiple best match(es):
                    max_index_match_keys = [k for k, v in index_matches.items() if v == max_index_match]
                    # If more than one read equally best-matched the index region, we can't pick a winner
                    if len(max_index_match_keys) != 1:
                        seq_assignment = {"adapter": "ambiguous_subset_umis",
                                          "mapping_obj": None}
                    # Otherwise, we have one read with more index matches than the others! Wooooot.
                    else:
                        seq_assignment = {"adapter": max_index_match_keys[0],
                                          "mapping_obj": hit_objs[max_index_match_keys[0]]}
                if keep_read_id:
                    seq_assignment["read_id"] = read_id

        # Else2: all the hits were the same match length (these could be unambiguous,
        #       b/c one might have no errors!)
        else:
            # Check if any of the reads spanned out to the index region (starting at: 1452)
            if max({k: hit.r_en for k, hit in hit_objs.items()}.values()) < 1452:
                if len(sequence) >= 1452:
                    seq_assignment = {"adapter": "ambiguous_all_map_too_short",
                                      "mapping_obj": None}
                else:
                    seq_assignment = {"adapter": "ambiguous_all_read_too_short",
                                      "mapping_obj": None}
                if print_per_seq:
                    print("Ambiguous, all too short")
            else:  # TODO: Add the umi loop to this step... Maybe?
                seq_assignment = {"adapter": "ambiguous_all",
                                  "mapping_obj": None}
                if print_per_seq:
                    print("Ambiguous, all adapters")
            if keep_read_id:
                seq_assignment["read_id"] = read_id

    # Else1: there were not any mapping hits:
    else:
        seq_assignment = {"adapter": "no_match",
                          "mapping_obj": None}
        if keep_read_id:
            seq_assignment["read_id"] = read_id
        if print_per_seq:
            print("No match")
    return seq_assignment


def _loop_pull_mappy_obj_info(mappy_obj, columns=None, return_dict=False):
    # Method to take the mapping object created by mappy, extract the sam
    # format information, and put that into a format for a dataframe.
    if mappy_obj is None:
        if return_dict:
            return dict(zip(columns, [None for _ in columns]))
        else:
            return [None for _ in columns]
    else:
        mappy_list = str(mappy_obj).split("\t")
        if return_dict:
            return dict(zip(columns, mappy_list))
        else:
            return mappy_list


def plot_value_counts(standards_df: pd.DataFrame, x: str = "adapter",
                      plot_fractions=False, only_standards=False,
                      title=None) -> None:
    import seaborn as sea
    import matplotlib.pyplot as plt
    import matplotlib.transforms as transforms

    sea.set_palette("colorblind")

    if only_standards:
        standards_df = standards_df[standards_df.adapter.isin(['30A_standard',
                                                               '15A_standard',
                                                               '10A_standard',
                                                               '5A_standard',
                                                               ])]

    # Set x_labels so that we can force order and pass that to the
    # number of observations annotations
    x_labels = standards_df[x].value_counts().keys()

    # Build the plot! The seaborn "countplot" is basically a bar plot
    # with the dataframe_column.value_counts step built in!
    fig = sea.countplot(x=x, data=standards_df, order=x_labels)

    # Rotate x labels so they aren't overlapping
    fig.set_xticklabels(fig.get_xticklabels(),
                        rotation=40, rotation_mode="anchor", ha='right')

    # Log y-scale, as the number of unmatched is wildly more than matches
    fig.set(yscale='log')

    # Custom transform so that I can place number of observations w/
    # X based on data (so they're on the columns), and Y based on
    # the paper/figure!
    custom_transform = transforms.blended_transform_factory(fig.transData,
                                                            fig.transAxes)

    # Below adapted from: https://www.python-graph-gallery.com/38-show-number-of-observation-on-boxplot
    # Calculate number of obs per group, and throw them on the plot!
    n_string_list = standards_df[x].value_counts().values
    if plot_fractions:
        count_sum = sum(n_string_list)
        n_string_list = [f"{i * 100 / count_sum:.2f}%" for i in n_string_list.tolist()]
    else:
        n_string_list = [f"n: {i}" for i in n_string_list.tolist()]

    # Add it to the plot
    pos = range(len(n_string_list))
    for tick, x_label in zip(pos, fig.get_xticklabels()):
        fig.text(pos[tick],
                 0.01,
                 n_string_list[tick],
                 transform=custom_transform,
                 horizontalalignment='center',
                 size='small',
                 color='k')
    if isinstance(title, str):
        fig.set(title=title)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':

    # W/ Dataframes:
    # Find the newest matching files for libraries of interest:
    path_dict = pick_libs_return_paths_dict(
        [
            "pTRI-stds-tera3",
            "pTRI-stds",
        ],
        file_suffix="parquet",
        file_midfix="mergedOnReads")
    for lib_key, lib_path in path_dict.items():
        lib_key += '_fromDF'
        df = pd.read_parquet(lib_path)
        df['read_length'] = df['sequence'].str.len()
        min_length = 1350
        max_length = 1550
        df = df.query(f"{min_length} <= read_length <= {max_length}")
        stds_df = align_standards(compressed_df=df,
                                  tail_less_reference=True,
                                  keep_read_id=True,
                                  )
        extra_df = df.merge(stds_df[['read_id', 'adapter']], on='read_id')
        extra_df.to_csv(f"{get_dt(for_file=True)}_{lib_key}_compressed-plus_wLen{min_length}-{max_length}.tsv", sep='\t')
        extra_df.to_parquet(f"{get_dt(for_file=True)}_{lib_key}_compressed-plus_wLen{min_length}-{max_length}.parquet")
        for only_stds in (True, False):
            plot_value_counts(stds_df,
                              only_standards=only_stds,
                              title=f"{lib_key}: {min_length} <= read_len <= {max_length};"
                                    f" only_stds={only_stds}")

    # W/ Fastq:
    # path_dict = pick_libs_return_paths_dict(
    #     [
    #         "pTRI-stds-tera3",
    #         "pTRI-stds",
    #     ],
    #     file_suffix='fastq',
    #     file_midfix='cat',
    #     output_dir_folder='cat_files')
    # for lib_key, lib_path in path_dict.items():
    #     stds_df = align_standards(fastq_file=lib_path,
    #                               tail_less_reference=True,
    #                               keep_read_id=True,
    #                               )
    #     parquet_path = list(pick_libs_return_paths_dict([lib_key],
    #                                                     file_suffix="parquet",
    #                                                     file_midfix="mergedOnReads").values())[0]
    #     df = pd.read_parquet(parquet_path)
    #     extra_df = df.merge(stds_df[['read_id', 'adapter']], on='read_id', how='right')
    #     extra_df.to_csv(f"{get_dt(for_file=True)}_{lib_key}_compressed-plus.tsv", sep='\t')
    #     extra_df.to_parquet(f"{get_dt(for_file=True)}_{lib_key}_compressed-plus.parquet")
    # 
    #     lib_key += '_fromFastq'
    #     for plot_frac in (True, False):
    #         for only_stds in (True, False):
    #             plot_value_counts(stds_df,
    #                               plot_fractions=plot_frac,
    #                               only_standards=only_stds,
    #                               title=f"{lib_key}: plot_frac={plot_frac}; only_stds={only_stds}")
        print("Done.")
