"""
mappyTestingCIGARWalking.py
Marcus Viscardi,    October 13, 2021

I want to check that the CIGAR walking I have been doing with the
heatmap scripts is actually getting to the same "end" that minimap2
would have called with the original mapping!
"""

from geneHeatmaps2 import _flip_neg_strand_genes, load_reads_tsv
from nanoporePipelineCommon import find_newest_matching_file
from standardsAlignment.standardsAssignmentWithMinimap2 import _loop_align_seq_to_adapters
import mappy
from pprint import pprint


def mappy_aligner(path_to_genome):
    mp_aligner = mappy.Aligner(path_to_genome,
                               preset="splice", k=14,
                               extra_flags=0x100000,  # From: https://github.com/lh3/minimap2/blob/master/minimap.h
                               n_threads=10,
                               )
    if not mp_aligner:
        raise Exception("ERROR: failed to load/build index")

    print(mp_aligner.seq_names)  # Prints the contigs which in this case are the adapters in the fastq
    return mp_aligner


def _loop_align_seq(aligner, name, seq, print_per_seq=False):
    if print_per_seq:
        print(f"\n\nAligning read: {name}")
    # For each sequence, create empty dict to store info about mapping hits
    hit_objs = {}

    # Loop through each of the mapping hits generated:
    for hit in aligner.map(seq):
        hit_objs[hit.ctg] = hit

    # If1: there are any mapping hits:
    if hit_objs:
        print("\n", name)
        for chr, hit in hit_objs.items():
            if hit.trans_strand == 1:
                strand = "+"
            else:
                strand = "-"
            calced_other_end = _flip_neg_strand_genes(hit.r_st, hit.cigar_str, strand)
            if strand == "+":
                print(strand, hit.r_st, hit.r_en, calced_other_end, hit.r_st == calced_other_end, sep="\t")
            else:  # strand == "-"
                print(strand, hit.r_st, hit.r_en, calced_other_end, hit.r_en == calced_other_end, sep="\t")
    #         seq_assignment = {"read_id": name,
    #                           "adapter": max_match_keys[0],
    #                           "mapping_obj": hit_objs[max_match_keys[0]]}
    # return seq_assignment


if __name__ == '__main__':
    
    working_dir = "/data16/marcus/working/210905_nanoporeRun_totalRNA_5108_xrn-1-KD"
    df = load_reads_tsv(find_newest_matching_file(f"{working_dir}/output_dir/"
                                                  f"merge_files/*_mergedOnReads.tsv"), head=100)

    genome_path = "/data16/marcus/genomes/elegansRelease100/200430_allChrs.fa"
    aligner = mappy_aligner(genome_path)
    
    seq_assignments = []
    for name, seq in zip(df["read_id"].to_list(), df["sequence"].to_list()):
        seq_assignment = _loop_align_seq(aligner, name, seq, print_per_seq=False)
    #     seq_assignments.append(seq_assignment)
    # pprint(seq_assignments)
