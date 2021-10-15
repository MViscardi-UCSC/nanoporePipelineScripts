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
import timeit

import mappy as mp
from pprint import pprint
import pandas as pd


def align_standards(fastq_file=None, compressed_df=None) -> pd.DataFrame:
    path_to_genome = "/data16/marcus/genomes/plus-pTRIxef_elegansRelease100/intermediate_files/5A_standard.fa"
    aligner = mp.Aligner(path_to_genome,
                         preset="splice", k=14,
                         extra_flags=0x100000,  # From: https://github.com/lh3/minimap2/blob/master/minimap.h
                         n_threads=10,
                         )
    if not aligner:
        raise Exception("ERROR: failed to load/build index")

    print(aligner.seq_names)  # Prints the contigs which in this case are the adapters in the fastq

    seq_assignments = []

    if isinstance(fastq_file, str):
        for read_id, sequence, _ in mp.fastx_read(fastq_file):
            seq_assignment = _loop_align_seq_to_adapters(aligner, read_id, sequence)
            seq_assignments.append(seq_assignment)
    elif isinstance(compressed_df, pd.DataFrame):
        # There is definitely a better way to do this than converting both to lists, but it's easy!
        for row_dict in compressed_df.to_dict(orient='records'):
            seq_assignment = _loop_align_seq_to_adapters(aligner, **row_dict)
            seq_assignments.append(seq_assignment)
    else:
        raise NotImplementedError("Please provide either a fastq or a df")

    print("Finished Alignment")
    return pd.DataFrame(seq_assignments)


def _loop_align_seq_to_adapters(aligner, read_id, sequence, print_per_seq=False, **kwargs):
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
                seq_assignment = {"read_id": read_id,
                                  "adapter": max_match_keys[0],
                                  "mapping_obj": hit_objs[max_match_keys[0]]}
                if print_per_seq:
                    print(seq_assignment)
                    read_end = seq_assignment["mapping_obj"].q_en
                    print(sequence[read_end - 10:read_end + 2], read_end)

            # Else3: there is more than 1 longer mapping hit (pretty weird situation?)
            else:
                seq_assignment = {"read_id": read_id,
                                  "adapter": "ambiguous_subset",
                                  "mapping_obj": None}
                if print_per_seq:
                    print("Ambiguous, subset of adapters")

        # Else2: all the hits were the same match length (these could be unambiguous,
        #       b/c one might have no errors!)
        else:
            seq_assignment = {"read_id": read_id,
                              "adapter": "ambiguous_all",
                              "mapping_obj": None}
            if print_per_seq:
                print("Ambiguous, all adapters")
            # I think I am actually handling these situations b/c I am using the hit.mlen variable, which counts
            #   the number of Ms (matches) in the CIGAR!! Very cool.

    # Else1: there were not any mapping hits:
    else:
        seq_assignment = {"read_id": read_id,
                          "adapter": "no_match",
                          "mapping_obj": None}
        if print_per_seq:
            print("No match")
    return seq_assignment


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


def mini_pull(map_obj):
    if map_obj:
        return map_obj.r_en
    else:
        return None


if __name__ == '__main__':
    # path_to_fa = "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/" \
    #              "output_dir/cat_files/cat.fastq"
    # df = align_standards(fastq_file=path_to_fa)
    # print(df)
    path_to_compressed = "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/output_dir/merge_files/210929_mergedOnReads.tsv"
    df = pd.read_csv(path_to_compressed, sep="\t", low_memory=False)
    df = df[df["chr_id"] == "pTRI"]
    df = df.head(10000)

    # Originally I built this wrapper to do below work:
    start1 = timeit.default_timer()
    out_df = align_standards(compressed_df=df)
    end1 = timeit.default_timer()

    # But I realized that I could achive the same end by doing the following! >>>
    start2= timeit.default_timer()
    path_to_genome = "/data16/marcus/genomes/plus-pTRIxef_elegansRelease100/intermediate_files/5A_standard.fa"
    aligner = mp.Aligner(path_to_genome,
                         preset="splice", k=14,
                         extra_flags=0x100000,  # From: https://github.com/lh3/minimap2/blob/master/minimap.h
                         n_threads=10)
    out_df = pd.DataFrame(df.apply(lambda row_dict:
                                   _loop_align_seq_to_adapters(aligner,
                                                               **row_dict),
                                   axis=1).to_list())
    end2 = timeit.default_timer()
    out_df["read_end"] = out_df.mapping_obj.apply(lambda map_obj: mini_pull(map_obj))
    # <<< So I think this is how I will add it to the main pipeline!

    plot_value_counts(out_df)
    print(f"Method 1 elapsed: {end1 - start1}")  # 17.19 sec / 10,000 lines
    print(f"Method 2 elapsed: {end2 - start2}")  # 18.62 sec / 10,000 lines
    # TODO: It would be good to have this merge back into the mergedOnReads dataframe!
