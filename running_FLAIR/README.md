# running_FLAIR/README.md

Generally I have used the FLAIR quantify tool 
to try and get transcript assignments for my 
nanopore reads.
This is really just bending a small part of 
FLAIR to my usage.
In reality, it's a really powerful tool that 
I could be taking better advantage of!

What I want to do in this directory is see if 
I can really leverage the strengths of FLAIR 
(like using DRIMSeq or DESeq2) to my advantage.

***

## *flair_noodling.ipynb*

All I am really trying to do here is get a 
handle on running *FLAIR quantify* at scale 
and *FLAIR diffexp*/*FLAIR diff_iso_usage*.

My first goal with this is going to be automatically 
creating a reads_manifest for select libraries and 
getting the output file for downstream analysis.

***

### *FLAIR quantify*
This is what I have been using for all of my 
"isoform assignment" stuff in a very hacky way.
In reality, this is used to make an output 
counts file more than anything, but I haven't
really cared about all that...

#### Here are some notes from the [readthedocs](https://flair.readthedocs.io/en/latest/modules.html#flair-quantify):

    flair quantify -r reads_manifest.tsv -i isoforms.fa [options]

##### Required Arguments:

    --isoforms          Fasta of Flair collapsed isoforms
    --reads_manifest    Tab delimited file containing sample id, condition, batch,
                        reads.fq, where reads.fq is the path to the sample fastq file.

##### Reads manifest example:

    sample1      condition1      batch1  mydata/sample1.fq
    sample2      condition1      batch1  mydata/sample2.fq
    sample3      condition1      batch1  mydata/sample3.fq
    sample4      condition2      batch1  mydata/sample4.fq
    sample5      condition2      batch1  mydata/sample5.fq
    sample6      condition2      batch1  mydata/sample6.fq

*Note: Do not use underscores in the first three fields, see below for details.*

##### Optional arguments:

    --help              Show all options
    --output            Name base for output files (default: flair.quantify). You
                        can supply an output directory (e.g. output/flair_quantify).
    --threads           Number of processors to use (default 4).
    --temp_dir          Directory to put temporary files. use ./ to indicate current
                        directory (default: python tempfile directory).
    --sample_id_only    Only use sample id in output header instead of a concatenation
                        of id, condition, and batch.
    --quality           Minimum MAPQ of read assignment to an isoform (default 1).
    --trust_ends        Specify if reads are generated from a long read method with
                        minimal fragmentation.
    --generate_map      Create read-to-isoform assignment files for each sample.
    --isoform_bed       isoform .bed file, must be specified if --stringent or
                        --check-splice is specified.
    --stringent         Supporting reads must cover 80% of their isoform and extend
                        at least 25 nt into the first and last exons. If those exons
                        are themselves shorter than 25 nt, the requirement becomes
                        'must start within 4 nt from the start' or 'end within 4 nt
                        from the end'.
    --check_splice      Enforces coverage of 4 out of 6 bp around each splice site
                        and no insertions greater than 3 bp at the splice site.

***

### *FLAIR diffExp*

Again, here are their notes from [readthedocs](https://flair.readthedocs.io/en/latest/modules.html#flair-diffexp):

    flair diffExp -q counts_matrix.tsv --out_dir out_dir [options]

This module performs differential expression and differential usage analyses between exactly two conditions with 3 or more replicates. It does so by running these R packages:

* **DESeq2** on genes and isoforms. This tests for differential expression.

* **DRIMSeq** is used on isoforms only and tests for differential usage. This is done by testing if the ratio of isoforms changes between conditions.

If you do not have replicates you can use the diff_iso_usage standalone script.

***

### *FLAIR diff_iso_usage*

Again, here are their notes from [readthedocs](https://flair.readthedocs.io/en/latest/scripts.html#diff-iso-usage):

    diff_iso_usage counts_matrix colname1 colname2 diff_isos.txt

Requires four positional arguments to identify 
and calculate significance of alternative isoform 
usage between two samples using Fisher’s exact tests: 
1. counts_matrix.tsv from flair-quantify
2. the name of the column of the first sample
3. the name of the column of the second sample
4. txt output filename containing the p-value associated with differential isoform usage for each isoform. The more differentially used the isoforms are between the first and second condition, the lower the p-value.

Output file format columns are as follows:

* gene name

* isoform name

* p-value

* sample1 isoform count

* sample2 isoform count

* sample1 alternative isoforms for gene count

* sample2 alternative isoforms for gene count
