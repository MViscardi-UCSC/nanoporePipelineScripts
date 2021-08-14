"""
geneHeatmaps2.py
Marcus Viscardi,    August 14, 2021

So I managed to overcomplicate everything while trying
    different methods in geneHeatmaps.py, so now I am
    going to try and take a simpler approach.

The plan is to import the mergedOnReads.tsv rather than
    any of the compressedOnGenes files. B/c the majority
    of the early processing steps to get distances from
    stops are on a per read basis: I'll use the read/row
    format of the mergedOnReads file to be able to
    iterate through things more cleanly/obviously.

CIGAR parsing will still be necessary to "flip" negative
    strand reads to get their 5' ends.

The other introduced change will be to entirely use the
    Caenorhabditis_elegans.WBcel235.100.allChrs.txt file
    which is produced by Josh's prepareReadAssignmentFile3.py
This annotation file is in the format:
    CHR_chrpos  WBGene:strand   iso1:posRelStart:posRelStop|iso2:...
The only striking issue is that Josh's dropped any locations with
    overlapping genes (which makes heaps of sense for short read
    sequencing, but may weaken my nanopore stuff).

After getting stop distance for each read, I'll be able to compress
    these on genes OR transcripts, calculate CDFs and feed them into
    the same plot function from geneHeatmaps.py (or I can rewrite
    that bit too. . . ).

Hype!
"""