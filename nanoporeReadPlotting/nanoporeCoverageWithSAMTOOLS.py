"""
nanoporeCoverageWithSAMTOOLS.py
Marcus Viscardi,    September 27, 2023

Going to try to rewrite a little of the nanoporeCoveragePlotting.py script to use samtools instead of cigar parsing!

It is FAST. I'll need to figure out how much faster it is than the cigar parsing method, but it is FAST.
"""
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm import tqdm

import nanoporePipelineCommon as npC

import pysam
import numpy


def simple_coverage_with_pysam(bam_path, target_chr=None):
    """
    This function will take a bam file path and return a numpy array of the coverage at each position in the chromosome.
    It is leveraging the pysam.AlignmentFile.pileup() method, which is very fast.
    """
    bam_obj = pysam.AlignmentFile(bam_path, "rb")  # Loads bam file as a pysam object
    chr_array = numpy.zeros([bam_obj.get_reference_length(target_chr)],
                            dtype=numpy.uint64)  # Makes a numpy array of zeros the length of the chr
    print(f"\nBuilding coverage array for chromosome: {target_chr}...", end=' ')
    for pileup_col in bam_obj.pileup(contig=target_chr):
        chr_array[pileup_col.reference_pos] = pileup_col.nsegments  # Store the coverage at each position
    print(f"Done.")
    bam_obj.close()
    return chr_array  # The returned array is just the entire chromosome, with coverage at each position


def run_coverage_with_pysam(bam_path_or_obj, target_chr=None, plot=False):
    """
    This function will take a bam file and return a numpy array of the coverage at each position in the chromosome.
    Added bells and whistles include:
        - A progress bar
        - The ability to take a pysam.AlignmentFile object instead of a path to a bam file
        - The ability to plot the coverage array
    """
    if isinstance(bam_path_or_obj, str) and Path(bam_path_or_obj).exists():
        pre_loaded_tag = False
        bam = pysam.AlignmentFile(bam_path_or_obj, "rb")
    elif isinstance(bam_path_or_obj, pysam.AlignmentFile):
        pre_loaded_tag = True
        bam = bam_path_or_obj
    else:
        raise TypeError(f"bam_path_or_obj must be a path to a bam file or a pysam.AlignmentFile object!")

    full_chr_len = bam.get_reference_length(target_chr)

    chr_array = numpy.zeros([full_chr_len], dtype=numpy.uint64)
    print(f"\nBuilding coverage array for chromosome: {target_chr}")
    with tqdm(desc=f"Iterating with pysam",
              total=full_chr_len) as progress_bar:
        for pileup_col in bam.pileup(contig=target_chr):
            progress_bar.update(pileup_col.reference_pos - progress_bar.n)
            chr_array[pileup_col.reference_pos] = pileup_col.nsegments
    print(f"Done.", end=" ")
    if not pre_loaded_tag:
        bam.close()

    if plot:
        print(f"Generating plot for chromosome: {target_chr}")
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.fill_between(x=numpy.arange(len(chr_array)),
                        y1=chr_array,
                        color='#303030')
        plt.show(dpi=300)
    else:
        print(f"Returning coverage array for chromosome: {target_chr}")

    return chr_array


def run_coverage_for_all_chrs(bam_path) -> dict:
    bam = pysam.AlignmentFile(bam_path, "rb")
    chr_dict = {}
    for chr_name in bam.references:
        chr_dict[chr_name] = run_coverage_with_pysam(bam, target_chr=chr_name, plot=False)
    bam.close()

    from pprint import pprint
    pprint(chr_dict)
    return chr_dict


def plot_all_chrs(chr_dict):
    fig, axes = plt.subplots(len(chr_dict), 1, figsize=(5, 3*len(chr_dict)))
    for i, (chr_name, chr_array) in enumerate(chr_dict.items()):
        print(f"Generating plot for chromosome: {chr_name}")
        axes[i].fill_between(x=numpy.arange(len(chr_array)),
                             y1=chr_array,
                             color='#303030')
        axes[i].set_title(chr_name)
    print(f"Drawing plot (slowly)...")
    plt.tight_layout()
    plt.show(dpi=96)


if __name__ == '__main__':
    bam_path = npC.pick_lib_return_path("newerN2", "cat_files", "cat.sorted.mappedAndPrimary", "bam")
    plot_all_chrs(run_coverage_for_all_chrs(bam_path))
