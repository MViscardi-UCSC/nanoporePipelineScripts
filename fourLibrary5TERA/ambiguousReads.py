"""
ambiguousReads.py
Marcus Viscardi,    July 03, 2023
Updated by Marcus Viscardi,    May 1, 2024

This script will iterate through a BAM file given from an NMD library, then it will identify reads that map to the
NMD-sensitive region of hardcoded genes. More specificially:
1. It will find the BAM file based on the library name by using the code in the nanoporePipelineCommon.py script.
2. It will load a hard-coded dictionary of genes that are known to be NMD-sensitive in C. elegans.
3. The NMD-sensitive genes ill be combined with a GTF reference annotations file to get the exact coordinates of the
   overall gene and other details.
4. Then for each NMD-sensitive gene, the script will iterate through the reads in the BAM file that overalap with the
   gene and determine if they overlap with the NMD-sensitive region of the gene.
5. The reads will be split into three categories based on their overlap with the NMD-sensitive region:
    a. Reads that overlap with the NMD-sensitive region (NMD_Targets)
    b. Reads that do not overlap with the NMD-sensitive region (Non_NMD_Targets)
    c. Reads that are ambiguous because they do not span the NMD-sensitive region (Ambiguous)
6. The script will output four BAM files:
    a. NMD_Targets.bam
    b. Non_NMD_Targets.bam
    c. Ambiguous.bam
    d. Merged.bam (all three of the above combined)
Importantly the outputted BAM files will have new tags added to them:
    - nC: This will be the category of the read (NMD_Targets, Non_NMD_Targets, Ambiguous)
    - gA: This will be the gene_id of the gene that the read overlapped with
The script will also output a CSV file with the details of the NMD-sensitive genes that were used in the analysis.
"""
import pysam
import pandas as pd
from pathlib import Path
import sys

sys.path.insert(0, '/data16/marcus/scripts/nanoporePipelineScripts')
from nanoporePipelineCommon import get_dt, pick_lib_return_path, CONVERSION_DICT

MIN_COVERAGE_OF_GENE = 50  # Min number of nucleotides mapped w/in the gene to consider it a hit
FRAC_COVERAGE_OF_REGION = 0.05  # Fraction of the region that must be covered to ID a isoform

NMD_ISO_REGIONS = {
    'rpl-30': {'start': 10_436_292,
               'end': 10_436_409,
               'region_is_target': True},
    'rps-22': {'start': 1_950_820,
               'end': 1_950_996,
               'region_is_target': True},
    'ubl-1': {'start': 3_068_573,
              'end': 3_068_581,
              'region_is_target': False},
    'rpl-7A': {'start': 4_389_805,
               'end': 4_389_886,
               'region_is_target': True},
    'rpl-3': {'start': 3_868_329,
              'end': 3_868_443,
              'region_is_target': True},
    'rpl-1': {'start': 2_876_019,
              'end': 2_876_144,
              'region_is_target': True},
    'rpl-12': {'start': 13_240_024,
               'end': 13_240_143,
               'region_is_target': True},
    'hel-1': {'start': 8_327_767,
              'end': 8_327_885,
              'region_is_target': True},
    'aly-3': {'start': 12_123_830,
              'end': 12_124_176,
              'region_is_target': True},
    'rsp-6': {'start': 7_790_498,
              'end': 7_790_745,
              'region_is_target': True},
    'K08D12.3': {'start': 1_710_371,
                 'end': 1_710_405,
                 'region_is_target': True},
    'R06C1.4': {'start': 11_931_061,
                'end': 11_931_173,
                'region_is_target': True},
    'C53H9.2': {'start': 1_833_384,
                'end': 1_833_408,
                'region_is_target': False},
    'rsp-5': {'start': 6_490_718,
              'end': 6_490_883,
              'region_is_target': True},
    # These following examples are a little shakier:
    'ZK228.4': {'start': 18_462_943,
                'end': 18_462_971,
                'region_is_target': True},  # Not sure about this one's NMD targeted region
    'rpl-26': {'start': 8_603_277,
               'end': 8_603_324,
               'region_is_target': True},  # rpl-26 is weird because it doesn't seem to have cleavage buildup/change
    'pqn-70': {'start': 11_227_418,
               'end': 11_227_453,
               'region_is_target': False},  # Somewhat unclear, but adding for posterity
    # 'examp': {'start': 3_068_573,
    #           'end': 3_068_581,
    #           'region_is_target': True},
}

CONVERSION_DICT = CONVERSION_DICT
REV_CONVERSION_DICT = {val: key for key, val in CONVERSION_DICT.items()}
LIB_NAMES = list(REV_CONVERSION_DICT.keys())

if __name__ == '__main__':
    # First additional command line argument will be the library name
    # The second will be the output directory
    if len(sys.argv) > 1:
        # Library name:
        if sys.argv[1] in LIB_NAMES:
            print(f"Using {sys.argv[1]} as lib_name")
            lib_name = sys.argv[1]
            bam_file = Path(pick_lib_return_path(REV_CONVERSION_DICT[lib_name], output_dir_folder='cat_files',
                                                 file_midfix='cat.sorted.mappedAndPrimary',
                                                 file_suffix='.bam'))
        elif sys.argv[1] in CONVERSION_DICT.keys():
            print(f"Using {REV_CONVERSION_DICT[sys.argv[1]]} as lib_name")
            lib_name = REV_CONVERSION_DICT[sys.argv[1]]
            bam_file = Path(pick_lib_return_path(REV_CONVERSION_DICT[lib_name], output_dir_folder='cat_files',
                                                 file_midfix='cat.sorted.mappedAndPrimary',
                                                 file_suffix='.bam'))
        else:
            print(f"Invalid library name: {sys.argv[1]}\nPlease use one of the following: {CONVERSION_DICT}")
            bam_file = None
            sys.exit(1)
        # Output directory:
        if sys.argv[2] == "input":
            long_output_names = False
            print(f"Saving all outputs in the input libraries output_dir under the directory: 'NMD_targets'")
            output_dir = bam_file.parent.parent / 'NMD_targets'
            output_dir.mkdir(exist_ok=True)
        else:
            long_output_names = True
            print(f"Saving all outputs in {sys.argv[2]}")
            output_dir = sys.argv[2]
            if not Path(output_dir).exists():
                Path(output_dir).mkdir(parents=True, exist_ok=True)
    else:
        print(f"Using default values for lib_name and output_dir")
        lib_name = "newerN2"
        output_dir = "/tmp"
        long_output_names = True
        bam_file = Path(pick_lib_return_path(REV_CONVERSION_DICT[lib_name], output_dir_folder='cat_files',
                                             file_midfix='cat.sorted.mappedAndPrimary', file_suffix='.bam'))

    print(f"Finished imports at: {get_dt(for_print=True)}")

    gtf_df = pd.read_parquet('/data16/marcus/genomes/'
                             'plus_cerENO2_elegansRelease100/'
                             '230327_allChrs_plus-cerENO2.gtf.parquet')
    targets_df = pd.DataFrame.from_dict(NMD_ISO_REGIONS,
                                        orient='index').reset_index().rename(columns={
        'index': 'gene_name',
        'start': 'nmd_region_start',
        'end': 'nmd_region_end',
        'region_is_target': 'nmd_region_is_target',
    })

    targets_df = targets_df.merge(
        gtf_df[gtf_df['gene_name'].isin(targets_df.gene_name) & (gtf_df['feature'] == 'gene')][
            ['gene_name', 'gene_id', 'chr', 'start', 'end', 'strand']])
    targets_df = targets_df[['gene_name',
                             'gene_id',
                             'chr',
                             'start',
                             'end',
                             'strand',
                             'nmd_region_start',
                             'nmd_region_end',
                             'nmd_region_is_target']]
    targets_df.to_csv(f"{output_dir}/identifiable_targets_df.csv", index=False)
    print(f"Saving targets to: {output_dir}/identifiable_targets_df.csv")

    bam = pysam.AlignmentFile(bam_file, 'rb')
    # Let's make temporary bam files for each of the below groups:
    # 1. Reads that hit the NMD-sensitive region
    # 2. Reads that did not hit the NMD-sensitive region
    # 3. Reads that were ambiguous b/c they don't span the putative NMD-sensitive region
    if long_output_names:
        general_name = f'{output_dir}/{bam_file.parent.parent.parent.stem}.all'
    else:
        general_name = f'{output_dir}/{get_dt()}_nmd_targets.all'
    
    # Output BAM paths:
    nmd_bam_path = Path(f'{general_name}.nmd.bam')
    non_nmd_bam_path = Path(f'{general_name}.non_nmd.bam')
    ambiguous_bam_path = Path(f'{general_name}.ambiguous.bam')
    
    # Output BAM pysam objects:
    nmd_bam = pysam.AlignmentFile(nmd_bam_path, 'wb', template=bam)
    non_nmd_bam = pysam.AlignmentFile(non_nmd_bam_path, 'wb', template=bam)
    ambiguous_bam = pysam.AlignmentFile(ambiguous_bam_path, 'wb', template=bam)

    for gene_id, gene_name in zip(targets_df['gene_id'], targets_df['gene_name']):
        chromosome, start, end, strand, nmd_start, nmd_end, region_is_nmd = targets_df.loc[targets_df['gene_id'] ==
                                                                                           gene_id,
                                                                                           [
                                                                                               'chr',
                                                                                               'start',
                                                                                               'end',
                                                                                               'strand',
                                                                                               'nmd_region_start',
                                                                                               'nmd_region_end',
                                                                                               'nmd_region_is_target',
                                                                                           ]].values[0]

        gene_is_reverse = (strand == '-')
        gene_is_forward = not gene_is_reverse
        last_edge = nmd_start if gene_is_reverse else nmd_end
        comparison_for_edge = f"<= {nmd_start}" if gene_is_reverse else f">= {nmd_end}"

        hit_count = 0
        non_hit_count = 0
        ambiguous_count = 0

        if region_is_nmd:
            hit = 'nmd_target'
            hit_bam = nmd_bam
            non_hit = 'non_nmd_target'
            non_hit_bam = non_nmd_bam
        else:
            hit = 'non_nmd_target'
            hit_bam = non_nmd_bam
            non_hit = 'nmd_target'
            non_hit_bam = nmd_bam

        min_nmd_overlap = int(FRAC_COVERAGE_OF_REGION * (nmd_end - nmd_start))
        if min_nmd_overlap <= 3:
            min_nmd_overlap = 1

        for i, read in enumerate(bam.fetch(chromosome, start, end)):
            if read.get_overlap(start, end) <= MIN_COVERAGE_OF_GENE:
                continue  # This will catch reads that are not in the gene! (50 is arbitrary)
                #           This is mainly because some reads "span" the gene, so they get fetched!
            if read.is_reverse != gene_is_reverse:
                continue  # This will catch reads that are in the wrong direction!
            if gene_is_reverse and eval(f"read.reference_end <= {nmd_start}"):
                ambiguous_count += 1
                read.set_tag('nC', 'ambiguous', value_type='Z')
                read.set_tag('gA', gene_id, value_type='Z')
                ambiguous_bam.write(read)
            elif not gene_is_reverse and eval(f"read.reference_start >= {nmd_end}"):
                ambiguous_count += 1
                read.set_tag('nC', 'ambiguous', value_type='Z')
                read.set_tag('gA', gene_id, value_type='Z')
                ambiguous_bam.write(read)
            elif read.get_overlap(nmd_start, nmd_end) >= min_nmd_overlap:
                hit_count += 1
                read.set_tag('nC', hit, value_type='Z')
                read.set_tag('gA', gene_id, value_type='Z')
                hit_bam.write(read)
            else:
                non_hit_count += 1
                read.set_tag('nC', non_hit, value_type='Z')
                read.set_tag('gA', gene_id, value_type='Z')
                non_hit_bam.write(read)

        if region_is_nmd:
            print(f"Gene: {gene_id} -- {gene_name} ({strand})\n\t"
                  f"nmd_count: {hit_count:>13,}\n\t"
                  f"non_nmd_count: {non_hit_count:>9,}\n\t"
                  f"ambig_count: {ambiguous_count:>11,}\n\t"
                  f"total_count: {hit_count + ambiguous_count + non_hit_count:>11,}")
        else:
            print(f"Gene: {gene_id} -- {gene_name} ({strand})\n\t"
                  f"nmd_count: {non_hit_count:>13,}\n\t"
                  f"non_nmd_count: {hit_count:>9,}\n\t"
                  f"ambig_count: {ambiguous_count:>11,}\n\t"
                  f"total_count: {hit_count + ambiguous_count + non_hit_count:>11,}")

    for bam in [nmd_bam, non_nmd_bam, ambiguous_bam]:
        file_name = bam.filename.decode()
        file_path = Path(file_name)
        print(f"Closing {file_path}")
        file_path_string = str(file_path)
        bam.close()
        pysam.sort(file_path_string, '-o', file_path_string)
        pysam.index(file_path_string)  # These samtools methods aren't in the pysam __init__ file, but they work!

    pysam.merge('-f', f'{general_name}.merge.bam',
                str(nmd_bam_path),
                str(non_nmd_bam_path),
                str(ambiguous_bam_path))
    pysam.sort(f'{general_name}.merge.bam', '-o', f'{general_name}.merge.bam')
    pysam.index(f'{general_name}.merge.bam')
    print(f"Merged all three categories into {general_name}.merge.bam")
    print(f"Done! Wrote all outputs to {output_dir}\n\n")
