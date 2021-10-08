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
    
    seq_spans, seq_matches = {}, {}
    seq_assignments = []
    
    if isinstance(fastq_file, str):
        for name, seq, _ in mp.fastx_read(fastq_file):
            seq_assignment, hit_spans, hit_matches = _loop_align_seq_to_adapters(aligner, name, seq)
            seq_assignments.append(seq_assignment)
            seq_spans[name] = hit_spans
            seq_matches[name] = hit_matches
    elif isinstance(compressed_df, pd.DataFrame):
        # There is definitely a better way to do this than converting both to lists, but it's easy!
        for name, seq in zip(compressed_df["read_id"].to_list(), compressed_df["sequence"].to_list()):
            seq_assignment, hit_spans, hit_matches = _loop_align_seq_to_adapters(aligner, name, seq)
            seq_assignments.append(seq_assignment)
            seq_spans[name] = hit_spans
            seq_matches[name] = hit_matches
    else:
        raise NotImplementedError("Please provide either a fastq or a df")
    
    print("Finished Alignment")
    return pd.DataFrame(seq_assignments)


def _loop_align_seq_to_adapters(aligner, name, seq, print_per_seq=False):
    if print_per_seq:
        print(f"Aligning read: {name}", end="\t")
    hit_spans = {}
    hit_matches = {}
    hit_objs = {}
    for hit in aligner.map(seq):  # Loop through the mapping hits generated
        hit_spans[hit.ctg] = [hit.r_st, hit.r_en]
        hit_matches[hit.ctg] = hit.mlen
        hit_objs[hit.ctg] = hit
    if hit_matches:  # If there are any mapping hits:
        if len(set(hit_matches.values())) > 1:  # If all the hits are not the same length!
            max_match = max(hit_matches.values())
            max_match_keys = [k for k, v in hit_matches.items() if v == max_match]
            if len(max_match_keys) == 1:  # If there is only one longest mapping hit:
                seq_assignment = {"read_id": name,
                                  "adapter": max_match_keys[0],
                                  "mapping_obj": hit_objs[max_match_keys[0]]}
                if print_per_seq:
                    print(max_match_keys[0])
                    pprint(hit_matches)
                    pprint(hit_spans)
                    print()
            else:  # There is more than 1 longer mapping hit (pretty weird situation?)
                seq_assignment = {"read_id": name,
                                  "adapter": "ambiguous_subset",
                                  "mapping_obj": None}
                if print_per_seq:
                    print("Ambiguous, subset of adapters")
        else:  # All the hits were the same match length (these could be unambiguous, b/c one might have no errors!)
            seq_assignment = {"read_id": name,
                              "adapter": "ambiguous_all",
                              "mapping_obj": None}
            if print_per_seq:
                print("Ambiguous, all adapters")
            # I think I am actually handling these situations b/c I am using the hit.mlen variable, which counts
            #   the number of Ms (matches) in the CIGAR!! Very cool.
    else:  # There were not any mapping hits:
        seq_assignment = {"read_id": name,
                          "adapter": "no_match",
                          "mapping_obj": None}
        if print_per_seq:
            print("No match")
    return seq_assignment, hit_spans, hit_matches


def plot_value_counts(standards_df: pd.DataFrame, x: str = "adapter"):
    import seaborn as sea
    import matplotlib.pyplot as plt
    import matplotlib.transforms as transforms
    
    sea.set_palette("colorblind")
    # sea.set_style("darkgrid")
    obs_names = standards_df[x].value_counts().keys()

    fig = sea.countplot(x=x, data=standards_df, order=obs_names)
    fig.set_xticklabels(fig.get_xticklabels(), rotation=40, ha="right")
    fig.set(yscale="log")
    num_obs = standards_df[x].value_counts().values
    num_obs = [str(x) for x in num_obs.tolist()]
    num_obs = ["n: " + i for i in num_obs]

    # Custom transform so that I can set X based on data, and
    # Y based on the paper
    trans = transforms.blended_transform_factory(fig.transData,
                                                 fig.transAxes)
    # Add it to the plot
    pos = range(len(num_obs))
    for tick, label in zip(pos, fig.get_xticklabels()):
        fig.text(pos[tick],
                 0.01,
                 num_obs[tick],
                 transform=trans,
                 horizontalalignment='center',
                 size='x-small',
                 color='k')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    # path_to_fa = "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/" \
    #              "output_dir/cat_files/cat.fastq"
    # df = align_standards(fastq_file=path_to_fa)
    # print(df)
    path_to_compressed = "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/output_dir/merge_files/210929_mergedOnReads.tsv"
    df = pd.read_csv(path_to_compressed, sep="\t", low_memory=False)
    df = df[df["chr_id"] == "pTRI"]
    df = df.head(10000)
    out_df = align_standards(compressed_df=df)
    plot_value_counts(out_df)
    print(out_df)
    # TODO: It would be good to have this merge back into the mergedOnReads dataframe!
