"""
samParsingFor5TERA.py
Marcus Viscardi,    November 19, 2021

Goal here is to test out parsing sam files and adding a new TAG that I can use
    for looking at 5TERA reads in IGV
    
11/21/21 - General Notes:
    So the main issue now is that I'll need to both run cutadapt early on to trim off
    adapters to allow for better mapping. And later integrate changed names from cutadapt.
    Importantly, I want to avoid changing these names for the in-between steps as they
    are currently working well and I'd rather not fuck with them.
    
    This sounds a lot like I will want to do 2 cutadapt calls... UGH.
    
    OR! We could add the cutadapt adapter identified thing as a comment in the fastqs.
    I then run with those trimmed and commented fastqs for the between steps (sadly the
    comments will be lost for a while between). Then use the fastq parsing class I have
    set up to extract the adapter information and reintegrate it to the final sam and bam
    files with simplesam.
    
    The above seems like the best way to go about this, so that's the plan:
        1. Make a copy of the fastq file called (cat.untrimmed.fastq)
        2. Run cutadapt with the rename part changing the comment rather than the read_id,
            and I'll want to have it trimming off the extra reads
        3. Run all the other stuff (ie. creating sam/bam files and so on)
        4. In a last step, parse the fastq comments and add those to the sam file as tags.
           -Probably use 't5' for 5TERA tags and 't3' for (eventual) TERA3 tags. On the sam
            format documentation (https://samtools.github.io/hts-specs/SAMv1.pdf) on pg 9
            they note that custom tags starting with X, Y, Z, or any lowercase letters
            are reserved for local use and should not have formally defined counterparts
            (this is to avoid overlaps!).
           -Additionally I think I'll do 'A' for if the (A)dapter is present or N for if the
           adapter is (N)ot present.
    
"""
import simplesam as ssam
from tqdm import tqdm

from step0_nanopore_pipeline import live_cmd_call
from nanoporePipelineCommon import FastqFile


def convert_read_id_tagged_sam(in_sam_path="/data16/marcus/working/"
                                           "211118_nanoporeRun_totalRNA_5108_xrn-1-KD_5TERA/"
                                           "output_dir/cat_files/deleteme.sam"):
    # The following loop will take stuff from a SAM file that was made by running minimap2 on
    #   a fastq that was run through cutadapt with the parameter: --rename '{id}:{adapter_name} {comment}'
    # This is a fine solution but is going to make an absolute MESS when it comes to the other things that
    #   read the fastq file (nanopolish and featureCounts, I think?)
    # For this reason, I'll need to very selectively keep or drop this other information as we go!
    #   Or maybe just make a copy of the fastq with the only change being the inclusion of the cutadapt
    #   tag to the name. Then I just feed the original to everything besides minimap2? then immediately
    #   after minimap2 I do the following to remove the cutadapt tag and make it a bam/sam tag?!
    with ssam.Reader(open(in_sam_path, 'r')) as in_sam:
        with ssam.Writer(open(in_sam_path.strip(".sam") + ".customComment.sam", 'w'),
                         in_sam.header) as out_sam:
            for read in in_sam:
                print(read)
                read["ZA"] = read.qname.split(":")[1]
                read.qname = read.qname.split(":")[0]
                print(read.tags["ZA"])
                out_sam.write(read)


if __name__ == '__main__':
    threads = 20
    output_dir_path = "/data16/marcus/working/211118_nanoporeRun_totalRNA_5108_xrn-1-KD_5TERA/output_dir"
    fastq_file_path = f"{output_dir_path}/cat_files/cat.fastq"
    output_fastq_path = f'{output_dir_path}/cat_files/cat.teraAdapters.fastq'
    # tera5_adapter = "AAUGAUACGGCGACCACCGAGAUCUACACUCUUUCCCUACACGACGCUCUUCCGANNN"  # Add X to front
    # cutadapt_call = f"cutadapt -g 5TERA={'X' + tera5_adapter} --action=trim -j {threads} " \
    #                f"--overlap 31 --error-rate 0.29 " \
    #                "--rename '{id} adapter={adapter_name} {comment}' " \
    #                f"{fastq_file_path} > {output_fastq_path}"
    # live_cmd_call(cutadapt_call)
    
    # Load the fastq with the additional adapter={} comment
    tagged_fastq = FastqFile(output_fastq_path)
    print(tagged_fastq)
    
    # \S+ matches any characters up to the first whitespace!
    tagged_fastq.df["adapter"] = tagged_fastq.df.comment.str.extract(r'adapter=(\S+)')
    tagged_fastq.df["tera5"] = tagged_fastq.df.adapter.replace({'no_adapter': '-',
                                                                'TERA3': '-',
                                                                '5TERA': '+'})
    tagged_fastq.df["tera3"] = tagged_fastq.df.adapter.replace({'no_adapter': '-',
                                                                'TERA3': '+',
                                                                '5TERA': '-'})
    tagged_fastq.df.set_index('read_id', inplace=True)
    with ssam.Reader(open(f"{output_dir_path}/cat_files/cat.sorted.mappedAndPrimary.bam", 'r')) as in_bam:
        with ssam.Writer(open(f"{output_dir_path}/cat_files/cat.sorted.mappedAndPrimary.tera.sam", 'w'),
                         in_bam.header) as out_sam:
            row_iterator = tqdm(in_bam)
            for read in row_iterator:
                row_iterator.set_description(f"Processing {read.qname}")
                read['t5'], read['t3'] = tagged_fastq.df.loc[read.qname, ['tera5', 'tera3']].tolist()
                out_sam.write(read)
