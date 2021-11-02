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
from nanoporePipelineCommon import find_newest_matching_file

pd.set_option("display.max_columns", None)


def align_standards(fastq_file=None, compressed_df=None, keep_read_id=False, bar_width=None,
                    keep_minimap_obj=False, threads=10, **kwargs) -> Union[str, pd.DataFrame]:
    path_to_genome = "/data16/marcus/genomes/plus-pTRIxef_elegansRelease100/" \
                     "intermediate_files/standards_with_indexes.fa"
    aligner = mp.Aligner(path_to_genome,
                         preset="splice", k=14,
                         extra_flags=0x100000,  # From: https://github.com/lh3/minimap2/blob/master/minimap.h
                         n_threads=threads,
                         )
    if not aligner:
        raise Exception("ERROR: failed to load/build index")

    # Print the contigs which in this case are the adapters in the fastq
    print(f"Mapping to adapters named: {aligner.seq_names}")

    seq_assignments = []

    if isinstance(fastq_file, str):
        row_iterator = tqdm(mp.fastx_read(fastq_file),
                            total=sum(1 for line in open(fastq_file)) // 4,
                            ncols=bar_width)
        for read_id, sequence, _ in row_iterator:
            seq_assignment = _loop_align_seq_to_adapters(aligner, read_id,
                                                         sequence, keep_read_id=keep_read_id,
                                                         was_fastq=True)
            seq_assignments.append(seq_assignment)
            row_iterator.set_description(f"Processing {read_id}")
    elif isinstance(compressed_df, pd.DataFrame):
        if not "pTRI" in compressed_df["chr_id"].to_list():
            return f"Chromosome pTRI not found in dataframe, please run minimap2 w/ genome that has pTRI added!!"
        row_iterator = tqdm(compressed_df.to_dict(orient='records'),
                            ncols=bar_width)
        for row_dict in row_iterator:
            seq_assignment = _loop_align_seq_to_adapters(aligner, keep_read_id=keep_read_id, **row_dict)
            seq_assignments.append(seq_assignment)
            row_iterator.set_description(f"Chr: {row_dict['chr_id']:>5}; Read: {row_dict['read_id']}")
    else:
        raise NotImplementedError("Please provide either a fastq or a df")

    print("Finished Alignment")

    mappy_df = pd.DataFrame(seq_assignments)

    if not keep_minimap_obj:
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
                      "stds_cigar"]
        mappy_df[mappy_cols] = mappy_df.apply(lambda row: _loop_pull_mappy_obj_info(row["mapping_obj"],
                                                                                    columns=mappy_cols),
                                              axis=1, result_type="expand")
        mappy_df.drop(columns="mapping_obj", inplace=True)
        for mappy_flag_col in ["stds_type_of_alignment", "stds_ts", "stds_cigar"]:
            mappy_df[mappy_flag_col] = mappy_df[mappy_flag_col].str.split(":").str[-1]
    return mappy_df


def _loop_align_seq_to_adapters(aligner, read_id, sequence, chr_id=None,
                                print_per_seq=False, keep_read_id=False, was_fastq=False,
                                **kwargs):
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
        hit_objs[hit.ctg] = hit

    # If1: there are any mapping hits:
    if hit_objs:

        # If2: all the hits are not the same length! (would be the same if they all
        #     had mapped to non-unique areas <- ie. not the UMI)
        if len(set(hit_objs.values())) > 1:
            hit_matches = {k: hit.mlen for k, hit in hit_objs.items()}
            max_match = max(hit_matches.values())
            max_match_keys = [k for k, v in hit_matches.items() if v == max_match]

            # If3: there is only one longest mapping hit:
            if len(max_match_keys) == 1:
                seq_assignment = {"adapter": max_match_keys[0],
                                  "mapping_obj": hit_objs[max_match_keys[0]]}
                if keep_read_id:
                    seq_assignment["read_id"] = read_id
                if print_per_seq:
                    print(seq_assignment)
                    read_end = seq_assignment["mapping_obj"].q_en
                    print(sequence[read_end - 10:read_end + 2], read_end)

            # Else3: there is more than 1 longer mapping hit (pretty weird situation?)
            else:
                seq_assignment = {"adapter": "ambiguous_subset",
                                  "mapping_obj": None}
                if keep_read_id:
                    seq_assignment["read_id"] = read_id
                if print_per_seq:
                    print("Ambiguous, subset of adapters")

        # Else2: all the hits were the same match length (these could be unambiguous,
        #       b/c one might have no errors!)
        else:
            seq_assignment = {"adapter": "ambiguous_all",
                              "mapping_obj": None}
            if keep_read_id:
                seq_assignment["read_id"] = read_id
            if print_per_seq:
                print("Ambiguous, all adapters")
            # I think I am actually handling these situations b/c I am using the hit.mlen variable, which counts
            #   the number of Ms (matches) in the CIGAR!! Very cool.

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


def plot_value_counts(standards_df: pd.DataFrame, x: str = "adapter"):
    import seaborn as sea
    import matplotlib.pyplot as plt
    import matplotlib.transforms as transforms

    sea.set_palette("colorblind")

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
    n_string_list = [str(x) for x in n_string_list.tolist()]
    n_string_list = ["n: " + i for i in n_string_list]

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
    plt.tight_layout()
    plt.show()


def mini_pull_map_obj_info(mapping_obj, **df_row):
    return_dict = {"adapter_start": None,
                   "adapter_end": None,
                   "ref_start": None,
                   "ref_end": None,
                   }
    if mapping_obj:
        return_dict["adapter_start"] = int(mapping_obj.q_st)
        return_dict["adapter_end"] = int(mapping_obj.q_en)
        return_dict["ref_start"] = int(mapping_obj.r_st)
        return_dict["ref_end"] = int(mapping_obj.r_en)
    return return_dict


def print_alignments_wrapper(adapter_mapped_df,
                             adapters_fasta="/data16/marcus/genomes/plus-pTRIxef_elegansRelease100/"
                                            "intermediate_files/standards_with_indexes.fa"):
    def _print_adapter_alignments(adapter_dict, read_id, adapter,
                                  adapter_start, adapter_end, ref_start, ref_end,
                                  sequence, **row):
        if adapter in adapter_dict.keys():
            adapter_start, adapter_end, ref_start, ref_end = map(int, [adapter_start, adapter_end, ref_start, ref_end])
            print(f"\nRead: {read_id},\tMapped to: {adapter}")
            print(f"adapter [{ref_start:>4}, {ref_end:>4}]: {adapter_dict[adapter][int(ref_start):int(ref_end)]}")
            print(f"adapter [{adapter_start:>4}, {adapter_end:>4}]: {sequence[int(adapter_start):int(adapter_end)]}")

    adapter_seq_dict = {}
    for adapter_name, sequence, _ in mp.fastx_read(adapters_fasta):
        adapter_seq_dict[adapter_name] = sequence
    pprint(adapter_seq_dict)
    for row_dict in adapter_mapped_df.to_dict(orient='records'):
        _print_adapter_alignments(adapter_seq_dict, **row_dict)


if __name__ == '__main__':
    # path_to_fa = "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/" \
    #              "output_dir/cat_files/cat.fastq"
    # df = align_standards(fastq_file=path_to_fa)
    path_to_compressed = "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/" \
                         "output_dir/merge_files/*_mergedOnReads.parquet"
    path_to_compressed = find_newest_matching_file(path_to_compressed)
    df = pd.read_parquet(path_to_compressed)

    # I built this wrapper to do below work:
    df = df.merge(align_standards(compressed_df=df, keep_read_id=True, keep_minimap_obj=True), on="read_id")

    # This step pulls information from the mapping_obj objects:
    # TODO: This step should probably overwrite information from minimap2 regarding the standards.
    #       Or at least rename the gene_id and all that!
    df = pd.concat([df,
                    pd.DataFrame(df.apply(lambda row_dict:
                                          mini_pull_map_obj_info(**row_dict),
                                          axis=1).to_list(),
                                 index=df.index)],
                   axis=1)

    plot_value_counts(df)
    # Method 1 (current):            17.19 sec / 10,000 lines
    # Method 2 (using df.apply()):   18.62 sec / 10,000 lines
    # TODO: It would be good to have this merge back into the mergedOnReads dataframe!
    #       Probably need to save as a parquet thou!

    # Random Stuff:
    # print_alignments_wrapper(df)
    print("Done!")
