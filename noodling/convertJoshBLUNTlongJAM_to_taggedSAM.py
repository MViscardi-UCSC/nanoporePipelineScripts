"""
convertJoshBLUNTlongJAM_to_taggedSAM.py
Marcus Viscardi,    August 09, 2023

Quick script to convert Josh's BLUNT JAM files to tagged SAM files
"""
from pathlib import Path
import pandas as pd
import sys
sys.path.insert(0, '/data16/marcus/scripts/nanoporePipelineScripts')
import nanoporePipelineCommon as npC
import pysam
from tqdm import tqdm
import argparse


def compress_blunt(blunt_string: str) -> str:
    compressed = ""
    count = 1

    for i in range(1, len(blunt_string)):
        if blunt_string[i] == blunt_string[i - 1]:
            count += 1
        else:
            compressed += str(count) + blunt_string[i - 1]
            count = 1

    # Add the last character and its count
    compressed += str(count) + blunt_string[-1]

    return compressed


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--input_longJAM", type=str, required=True,
                        help="Path to longJAM file from Josh's BLUNT pipeline")
    parser.add_argument("-s", "--input_SAM", type=str, required=True,
                        help="Path to SAM file from Josh's BLUNT pipeline")
    parser.add_argument("-o", "--output", type=str, required=False, default=None,
                        help="Path to output SAM file, defaults to **input_SAM_BLUNT**_tagged.sam")
    parser.add_argument("-d", "--drop_untagged", action="store_true",
                        help="Drop reads that are not in the longJAM file")
    args = parser.parse_args()
    
    longJAM_path = Path(args.input_longJAM)
    SAM_path = Path(args.input_SAM)
    if args.output is None:
        output_path = SAM_path.parent / f"{SAM_path.stem}_BLUNTtagged.sam"
    else:
        output_path = args.output
        if not output_path.endswith(".sam"):
            output_path += ".sam"
        output_path = Path(output_path)
    
    longJAM_df = pd.read_csv(longJAM_path, header=None, sep='\t', skiprows=14, names=list(range(24)))
    read_blunt_df = pd.concat(
        [
            longJAM_df.loc[::2, 0].reset_index(drop=True).rename('read_name'),
            longJAM_df.loc[1::2, 0].apply(lambda x: compress_blunt(x)).reset_index(drop=True).rename('blunt_string'),
        ], axis=1)
    
    sam = pysam.AlignmentFile(SAM_path, "rb")
    out_sam = pysam.AlignmentFile(output_path, "w", template=sam)
    for read in sam:
        try:
            read.set_tag("bL", read_blunt_df.loc[read.query_name, 'blunt_string'])
            out_sam.write(read)
        except KeyError:
            print(f"Could not find read: {read.query_name} in longJAM file,", end="")
            if args.drop_untagged:
                print(" dropping read.")
            else:
                print(" keeping read with NA for bL tag.")
                read.set_tag("bL", "NA")
                out_sam.write(read)
    print(f"Finished writing to {output_path}")
