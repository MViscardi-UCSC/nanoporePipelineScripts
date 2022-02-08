"""
nanoporePipelineCommon.py
Marcus Viscardi,    September 22, 2021

A common location for some often used methods.
"""
import pandas as pd
import numpy as np
import os

pd.set_option('display.width', 400)
pd.set_option('display.max_columns', None)


class FastqFile:
    def __init__(self, path, head=None):
        from tqdm import tqdm
        import mappy as mp
        self.path = path
        self.subset = isinstance(head, int)
        fastq_items_list = []  # List to hold fastq items as we iterate over
        print(f"\nStarting fastq iterative load @ {get_dt(for_print=True)}", end="\n")
        # Below will allow us to track how long the fastq parse is taking:
        row_iterator = tqdm(mp.fastx_read(path, read_comment=True),
                            total=sum(1 for line in open(path)) // 4)
        for line, (read_id, sequence, quality, comment) in enumerate(row_iterator):
            fastq_items_list.append([read_id, sequence, "+", quality, comment])
            row_iterator.set_description(f"Processing {read_id}")
            if self.subset and line >= head:
                break
        # Convert the fastq items list into a pandas dataframe so it can be filtered by the alt_mapped_reads_df
        self.df = pd.DataFrame(fastq_items_list, columns=["read_id", "sequence", "plus", "quality", "comment"])

    def __str__(self):
        return str(self.df)

    def filter_against(self, df_w_read_id):
        # Cool way to only keep values that don't appear in alt_mapped_read_df:
        #   (From: https://tinyurl.com/22czvzua)
        # Also trying out query operator, some notes on this here:
        #   https://stackoverflow.com/questions/67341369/pandas-why-query-instead-of-bracket-operator
        self.df = pd.merge(self.df, df_w_read_id,
                           on="read_id",
                           indicator=True,
                           how="outer").query('_merge=="left_only"').drop('_merge', axis=1)

    def save_to_fastq(self, output_path):
        from csv import QUOTE_NONE
        self.df["read_id"] = "@" + self.df["read_id"] + " " + \
                             self.df["comment"]
        self.df.drop("comment", axis=1).to_csv(output_path,
                                               index=False,
                                               header=False,
                                               sep="\n",
                                               quoting=QUOTE_NONE)


class SamOrBamFile:
    def __init__(self, path, header_source=None, subsample=None):
        from subprocess import check_output

        self.path = path
        if not os.path.exists(self.path):
            raise FileNotFoundError(f"Could not find file: {self.path}")
        self.max_cols = 30  # Probably shouldn't need to be changeable...

        self.file_type = self.path.split('.')[-1]
        if self.file_type not in ['bam', 'sam']:
            raise NotImplementedError(f"Please provide a file that ends with bam/sam, "
                                      f"your file ends in: {self.file_type}")
        self.sam_main_columns = ["read_id",
                                 "bit_flag",
                                 "chr_id",
                                 "chr_pos",
                                 "mapq",
                                 "cigar",
                                 "r_next",
                                 "p_next",
                                 "len",
                                 "sequence",
                                 "phred_qual"]
        # Anything not included below will either stay as a string/object,
        #   or whatever the simplesam decided to store it as (for tags):
        self.dtype_dict = {'bit_flag': 'uint16',
                           'chr_id': 'category',
                           'chr_pos': 'uint32',
                           'mapq': 'uint8',
                           'NM': 'uint16',
                           'ms': 'uint32',
                           'AS': 'uint32',
                           'nn': 'uint32',
                           'ts': 'category',
                           'tp': 'category',
                           'cm': 'uint32',
                           'rl': 'uint32',
                           't5': 'category',
                           't3': 'category'}
        self.tag_type_dict = {'NM': 'i',
                              'ms': 'i',
                              'AS': 'i',
                              'nn': 'i',
                              'ts': 'A',
                              'tp': 'A',
                              'cm': 'i',
                              's1': 'i',
                              's2': 'i',
                              'de': 'f',
                              'rl': 'i',
                              't5': 'A',
                              't3': 'A',
                              'SA': 'Z',
                              }

        if isinstance(header_source, str):
            self.header = check_output(f"samtools view -H {header_source}", shell=True).decode("utf-8")
            read_file_header = check_output(f"samtools view --no-PG -H {path}", shell=True).decode("utf-8")
            self.header_lines = len(read_file_header.split('\n')) - 1
        else:
            self.header = check_output(f"samtools view --no-PG -H {path}", shell=True).decode("utf-8")
            self.header_lines = len(self.header.split('\n')) - 1

        if isinstance(subsample, int):
            self.subsample = subsample
        else:
            self.subsample = None

        self.df = self.__build_df__()

    def __build_df__(self):
        import simplesam as ssam
        from tqdm import tqdm
        from collections import ChainMap

        print(f"{get_dt(for_print=True)}: Starting to load {self.file_type} file from: {self.path}")
        sam_list_of_dicts = []
        with ssam.Reader(open(self.path, 'r')) as in_sam:
            to_subsample = isinstance(self.subsample, int)
            if to_subsample:
                max_to_load = self.subsample
            else:
                if self.file_type == 'bam':
                    max_to_load = len(in_sam)
                else:
                    max_to_load = sum(1 for line in open(self.path))
            sam_iterator = tqdm(enumerate(in_sam), total=max_to_load)
            for i, read in sam_iterator:
                sam_iterator.set_description(f"Processing {read.qname}")
                main_info_list = read.__str__().split('\t')[0:11]
                main_info = dict(zip(self.sam_main_columns, main_info_list))
                tag_info = read.tags
                all_info = dict(ChainMap(main_info, tag_info))
                sam_list_of_dicts.append(all_info)
                if to_subsample and i > self.subsample:
                    break
        final_df = pd.DataFrame(sam_list_of_dicts)
        self.tag_columns = [col for col in final_df if col not in self.sam_main_columns]
        final_df = final_df[self.sam_main_columns + self.tag_columns]

        # Filter the datatype dictionary, so we don't get errors from
        #   trying to retype columns that don't exist!
        self.dtype_dict = {k: v for k, v in self.dtype_dict.items() if k in list(final_df.columns)}
        # Adjust any datatypes we want:
        final_df = final_df.astype(self.dtype_dict)
        print(f"{get_dt(for_print=True)}: Finished loading {self.file_type} file!")
        return final_df

    def to_sam(self, output_path=None, escape_char="|", to_bam=False,
               infer_tag_types=False):
        import simplesam as ssam
        from subprocess import run
        from csv import QUOTE_NONE
        if escape_char not in "~}|":
            raise NotImplementedError("escape_char has to be |, }, or ~. Everything else is in the Phred Quals!")
        save_df = self.df.copy()
        if infer_tag_types:
            print(f"Inferring tag datatypes is sorta broken... Use at your own risk!")
            for col in self.tag_columns:
                save_df[col] = save_df[col].apply(lambda x: ssam.encode_tag(col, x))
        else:
            for col in self.tag_columns:
                try:
                    save_df[col] = col + ':' + self.tag_type_dict[col] + ':' + save_df[col].astype(str)
                except KeyError:
                    print(f"{col} not found in tag_type_dict and caused a KeyError??")
                    print(f"For now, we are going to use simplesam to infer the correct type!")
                    save_df[col] = save_df[col].apply(lambda x: ssam.encode_tag(col, x))
                    # save_df[col] = col + ':Z:' + save_df[col].map(str)
        save_df["tags"] = save_df[self.tag_columns].values.tolist()
        save_df["tags"] = save_df["tags"].apply(lambda row: '\t'.join([x for x in row if not str(x).endswith('nan')]))
        save_df = save_df.drop(self.tag_columns, axis=1)

        temp_header = self.header + f"@CO\t{get_dt(for_print=True)}: Introduced edits with python . . . " \
                                    f"More info hopefully added to this later!\n"
        buffer = temp_header + save_df.to_csv(sep="\t",
                                              header=False,
                                              index=False,
                                              quoting=QUOTE_NONE,
                                              escapechar=escape_char,
                                              quotechar='\0').replace(escape_char, '')
        # subprocess.run accepts the input param to pass to the bash call!
        if output_path:
            if to_bam:
                run(f"samtools view -h --no-PG -b - > {output_path} && samtools index {output_path}",
                    input=buffer.encode('utf-8'), shell=True)
            else:
                run(f"samtools view -h --no-PG - > {output_path}",
                    input=buffer.encode('utf-8'), shell=True)
        else:
            print(save_df)
        return save_df

    def to_bam(self, output_path):
        self.to_sam(output_path=output_path, to_bam=True)


def load_and_merge_lib_parquets(lib_list, genomeDir=f"/data16/marcus/genomes/elegansRelease100/",
                                drop_unassigned=True, drop_failed_polya=True,
                                drop_sub_n=5, keep_transcript_info=False,
                                read_pos_in_groupby=False, add_nucleotide_fractions=False,
                                group_by_t5=False) -> [pd.DataFrame, pd.DataFrame]:
    concat_df = library_reads_df_load_and_concat(lib_list, genomeDir, drop_unassigned, drop_failed_polya,
                                                keep_transcript_info, group_by_t5)

    compressed_df = compress_concat_df_of_libs(concat_df, drop_sub_n, read_pos_in_groupby, add_nucleotide_fractions,
                                               keep_transcript_info, group_by_t5)

    return concat_df, compressed_df


def library_reads_df_load_and_concat(lib_list, genomeDir, drop_unassigned, drop_failed_polya, keep_transcript_info,
                                     group_by_t5):
    # Initially load the read assignments file:
    read_assignment_path = find_newest_matching_file(f"{genomeDir}/*.allChrs.parquet")
    print(f"Loading readAssignments file from: {read_assignment_path}...", end=' ')
    read_assignment_df = pd.read_parquet(read_assignment_path)
    print("Done.")
    # Loop through each library name in the list and for each:
    #   1. Load the Parquet
    #   2. Merge this w/ Josh's assign reads based on chr_pos
    #   3. Create a column to retain the library identity
    #   3. Concatenate these dataframe into one large dataframe
    #       NOTE: This structure can be seperated again based on
    #       the "lib" column added in the previous step
    print(f"Looking for files for libraries: {lib_list}")
    path_dict = pick_libs_return_paths_dict(lib_list,
                                            output_dir_folder="merge_files",
                                            file_midfix="_mergedOnReads",
                                            file_suffix="parquet")
    # Dictionary to hold the library names and dataframes before the merge
    df_dict = {}
    for library_name, parquet_path in path_dict.items():
        # Load each individual library:
        print(f"Loading parquet for {library_name} lib...", end=' ')
        lib_df = pd.read_parquet(parquet_path)
        print("Done.")

        # Make sure that 5' end adjustments happened, do so if not!
        lib_df = adjust_5_ends(lib_df)
        # Store the library dataframe in the temporary dictionary
        df_dict[library_name] = lib_df
    # This is a cute way to quickly merge all of these dfs into one, while retaining lib info.
    #   B/c I am still a little scared of MultiIndexed dataframes, I used the reset_index steps
    #   to push the mutliindex back into columns. Maybe someday I'll use the multiindex!
    multi_df = pd.concat(df_dict.values(),
                         keys=df_dict.keys())
    multi_df.index.set_names(("lib", "old_index"), inplace=True)
    concat_df = multi_df.reset_index(level="lib").reset_index(drop=True)
    # Some major issues with ugly column names coming through, plan to clean them up:
    read_assignment_cols = read_assignment_df.columns.to_list()
    read_assignment_cols.remove('chr_id')
    read_assignment_cols.remove('chr_pos')
    # Only retain columns that don't show up in both the read assignment df and the super merge df:
    read_assignment_cols_to_drop = [col for col in read_assignment_cols if col in concat_df.columns]
    # Drop those 'unique' columns
    concat_df.drop(columns=read_assignment_cols_to_drop,
                  inplace=True)
    print(f"Starting assignment merge . . .", end="")
    # Add read assignments w/ josh's read_assignment dataframe
    concat_df = concat_df.merge(read_assignment_df, on=["chr_id", "chr_pos"],
                              how="left", suffixes=["_originalOutput",
                                                    ""])
    print(f"\rFinished assignment merge!")
    # To further clean up mixed columns, just retain the ones we care about!
    keep_columns = ["lib",
                    "read_id",
                    "chr_id",
                    "chr_pos",
                    "gene_id",
                    "gene_name",
                    "cigar",
                    "sequence",
                    "polya_length",
                    "strand",
                    ]
    if keep_transcript_info:
        # Add transcript info to the columns we care about if requested
        for col in ["transcript_id", "to_start", "to_stop"]:
            keep_columns.append(col)
    if group_by_t5:
        # Add 5TERA information to the columns we care about if requested
        keep_columns.append('t5')
    # Drop unnecessary columns, reassess for duplicates.
    #   i.e. reads that mapped to two transcripts will
    #        otherwise be dups if we drop transcript info!
    concat_df = concat_df[keep_columns].drop_duplicates()
    # Only retain polyA passed reads if requested
    if drop_failed_polya:
        print(f"\nRead counts pre-failed-polyA call drop:   {concat_df.shape[0]}")
        concat_df = concat_df[~concat_df["polya_length"].isna()]
        print(f"Read counts post-failed-polyA call drop:  {concat_df.shape[0]}")
    # Post-assignment cleanups:
    if drop_unassigned:
        print(f"Read counts post gene assignment:  {concat_df.shape[0]}")
        concat_df = concat_df[~concat_df["gene_id"].isna()].reset_index(drop=True)
        print(f"Read counts post unassigned drop:  {concat_df.shape[0]}")
        # 220201: I don't know why we were doing below...?
        # concat_df = concat_df[concat_df["strand_originalOutput"] == concat_df["strand_forPlot"]].reset_index(drop=True)
        # print(f"Read counts post consistent-assignment check: {concat_df.shape[0]}")
        if keep_transcript_info:
            concat_df = concat_df.astype({"to_stop": "int64",
                                        "to_start": "int64"})
    # Add a read length column!
    concat_df["read_length"] = concat_df["sequence"].apply(len)
    return concat_df


def compress_concat_df_of_libs(super_df, drop_sub_n, read_pos_in_groupby, add_nucleotide_fractions, keep_transcript_info,
                               group_by_t5):
    # Create the groupby dataframe:
    groupby_col_list = ["lib", "chr_id", "gene_id", "gene_name"]
    print(f"Creating groupby dataframe merged on: {groupby_col_list}")
    if keep_transcript_info:
        print(f"\t+ [transcript_id]")
        groupby_col_list.append("transcript_id")
    if group_by_t5:
        print(f"\t+ [t5] tag")
        groupby_col_list.append("t5")
    # Holy crap, the observed=True helps to keep this from propagating out to 129,151,669,691,968 rows...
    groupby_obj = super_df.groupby(groupby_col_list, observed=True)
    if not keep_transcript_info:
        compressed_prefix = "gene"
    else:
        compressed_prefix = "transcript"
    compressed_df = groupby_obj["read_id"].apply(len).to_frame(name=f"{compressed_prefix}_hits")
    compressed_df["mean_polya_length"] = groupby_obj["polya_length"].mean()
    compressed_df["mean_read_length"] = groupby_obj["read_length"].mean()
    if read_pos_in_groupby:
        compressed_df['stop_distances'] = groupby_obj["to_stop"].apply(list).to_frame(name="stop_distances")
        compressed_df['start_distances'] = groupby_obj["to_start"].apply(list).to_frame(name="stop_distances")
    # RPM and fractional hits calculations
    # Need to first create columns of NA values, tobe overwritten
    compressed_df[f"{compressed_prefix}_rpm"] = pd.NA
    compressed_df[f"{compressed_prefix}_frac_hits"] = pd.NA
    # Only look at one library at a time (so the normalization is per lib not whole df)
    for lib in compressed_df.index.unique(level='lib').to_list():
        # Create the 'norm_factor' which will be the total # of read hits in that lib
        norm_factor = compressed_df.query(f"lib == '{lib}'")[f"{compressed_prefix}_hits"].sum()
        # Turn the total number of read hits into the 'million of read hits'
        rpm_norm_factor = norm_factor / 1000000
        # For each library divide gene_hits by the rpm norm factor to get rpm
        gene_rpm_series = compressed_df.query(f"lib == '{lib}'")[f"{compressed_prefix}_hits"] / rpm_norm_factor
        # Use a series fill, so that we can fill that library's part of the DF without effecting others
        compressed_df[f"{compressed_prefix}_rpm"] = compressed_df[f"{compressed_prefix}_rpm"]. \
            fillna(value=gene_rpm_series)
        # Same as above, but with fraction of hits, rather than a rpm calc (practically same thing)
        gene_frac_hits_series = compressed_df.query(f"lib == '{lib}'")[f"{compressed_prefix}_hits"] / norm_factor
        compressed_df[f"{compressed_prefix}_frac_hits"] = compressed_df[f"{compressed_prefix}_frac_hits"]. \
            fillna(value=gene_frac_hits_series)
    # Requirement for min number of gene/transcript hits
    if isinstance(drop_sub_n, int):
        print(f"Gene counts pre sub-{drop_sub_n} {compressed_prefix}_hits drop:  {compressed_df.shape[0]}")
        compressed_df = compressed_df[compressed_df[f"{compressed_prefix}_hits"] >= drop_sub_n]
        print(f"Gene counts post sub-{drop_sub_n} {compressed_prefix}_hits drop:  {compressed_df.shape[0]}")
    # Reset index at the end,
    #   we didn't retain any info w/ the index, so it doesn't help much
    compressed_df = compressed_df.reset_index()
    if add_nucleotide_fractions:
        if not keep_transcript_info:
            print(f"Adding nucleotide content information to genes!")
            path_to_gc = "/data16/marcus/genomes/elegansRelease100/" \
                         "Caenorhabditis_elegans.WBcel235.cdna.all.fa.GCcontent.genes.parquet"
            gc_df = pd.read_parquet(path_to_gc).drop(columns=["chr_id"])
            compressed_df = compressed_df.merge(gc_df, on=["gene_id", "gene_name"], how="left")
        else:
            print(f"Adding nucleotide content information to transcripts!")
            path_to_gc = "/data16/marcus/genomes/elegansRelease100/" \
                         "Caenorhabditis_elegans.WBcel235.cdna.all.fa.GCcontent.transcripts.parquet"
            gc_df = pd.read_parquet(path_to_gc).drop(columns=["chr_id"])
            compressed_df = compressed_df.merge(gc_df, on=["gene_id", "gene_name", "transcript_id"], how="left")
    return compressed_df


def sam_or_bam_class_testing():
    test_output_dir = "/data16/marcus/working/211101_nanoporeSoftLinks/220131_nanoporeRun_totalRNA_0639_L3_third/output_dir"
    test_sam_path = f"{test_output_dir}/cat_files/cat.sorted.mappedAndPrimary.bam"
    # test_sam_path = "/data16/marcus/working/211101_nanoporeSoftLinks/" \
    #                 "211118_nanoporeRun_totalRNA_5108_xrn-1-KD_5TERA/" \
    #                 "output_dir/cat_files/cat.sorted.mappedAndPrimary.sam"
    sam = SamOrBamFile(test_sam_path,
                       # subsample=50000,
                       )
    polya_df = pd.read_csv(f"{test_output_dir}/nanopolish/polya.passed.tsv", sep="\t")
    # For some god-awful reason the chr_pos in polyA are -1 to those in the SAM file:
    polya_df["position"] = polya_df["position"] + 1
    polya_df = polya_df.rename(columns={"readname": "read_id",
                                        "qc_tag": "qc_tag_polya",
                                        "position": "chr_pos",
                                        "contig": "chr_id",
                                        "polya_length": "pA"})
    polya_df = polya_df[['read_id', 'chr_id', 'chr_pos', 'pA']]
    sam.df = pd.merge(sam.df, polya_df, on=['read_id', 'chr_id', 'chr_pos'], how='left')
    sam.df.pA.fillna(0, inplace=True)
    sam.tag_columns.append('pA')
    # sam.to_sam(test_sam_path + '.resave.sam', escape_char="~")
    # sam.to_sam(output_path=test_sam_path.rstrip('.sam') + '.resave.sam')
    print(sam)


# "riboD", "totalRNA", "totalRNA2", "polyA", "polyA2",
# "xrn-1", "xrn-1-5tera", "pTRI-stds", "xrn-1-5tera-smg-6", "pTRI-stds-tera3"
def pick_libs_return_paths_dict(lib_list: list, file_suffix: str = "parquet", file_midfix="_mergedOnReads",
                                output_dir_folder="merge_files", return_all: bool = False) -> dict:
    output_dir_dict = {
        "riboD": "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/output_dir",
        "totalRNA": "/data16/marcus/working/210709_NanoporeRun_totalRNA_0639_L3/"
                    "output_dir",
        "totalRNA2": "/data16/marcus/working/"
                     "210720_nanoporeRun_totalRNA_0639_L3_replicate/output_dir",
        "polyA": "/data16/marcus/working/210528_NanoporeRun_0639_L3s/output_dir",
        "polyA2": "/data16/marcus/working/210719_nanoporeRun_polyA_0639_L3_replicate/output_dir",
        "xrn-1": "/data16/marcus/working/210905_nanoporeRun_totalRNA_5108_xrn-1-KD/output_dir",
        "xrn-1-5tera": "/data16/marcus/working/211118_nanoporeRun_totalRNA_5108_xrn-1-KD_5TERA/output_dir",
        "pTRI-stds": "/data16/marcus/working/211121_nanoporeRun_pTRIstds/output_dir",
        "xrn-1-5tera-smg-6": "/data16/marcus/working/211210_nanoporeRun_totalRNA_2102_xrn-1-KD_5TERA/output_dir",
        "pTRI-stds-tera3": "/data16/marcus/working/211212_nanoporeRun_pTRIstds_TERA3/output_dir",
        "polyA3": "/data16/marcus/working/220131_nanoporeRun_polyA_0639_L3_third/output_dir",
        "totalRNA3": "/data16/marcus/working/220131_nanoporeRun_totalRNA_0639_L3_third/output_dir",
    }
    if return_all:
        lib_list = output_dir_dict.keys()
    file_suffix = file_suffix.strip(".")
    return_dict = {}
    for lib_key, output_dir in output_dir_dict.items():
        if lib_key in lib_list:
            file_path = f"{output_dir}/{output_dir_folder}/*{file_midfix}.{file_suffix}"
            print(f"Looking for file for {lib_key}, at {file_path}...", end=" ")
            return_dict[lib_key] = find_newest_matching_file(file_path)
            print(f"File Found.")
    if not return_dict:
        raise KeyError(f"No matching library keys found, please use keys from the following list: "
                       f"{list(output_dir_dict.keys())}. "
                       f"You passed the following keys: "
                       f"{lib_list}")
    return return_dict


def load_read_assignments(assignment_file_parquet_path) -> pd.DataFrame:
    print(f"Loading read assignment file from: {assignment_file_parquet_path} ", end="")
    try:
        read_assignment_df = pd.read_parquet(assignment_file_parquet_path)
    except FileNotFoundError:
        raise FileNotFoundError(f"Couldn't find parquet file! Maybe try to run the method: "
                                f"parse_read_assignment_allChrs_txt in nanoporePipelineCommon.py or "
                                f"first: prepareReadAssignmentFily3.py in the ArribereLab git!!")
    return read_assignment_df


def old_assign_w_josh_method(reads_df, genomeDir):
    def merge_on_chr_pos(read_assignment_df: pd.DataFrame, reads_df: pd.DataFrame) -> pd.DataFrame:
        print(f"Merging read assignments and reads at {get_dt(for_print=True)}")
        merge_df = reads_df.merge(read_assignment_df, on=["chr_id", "chr_pos"],
                                  how="left", suffixes=("_fromFeatureCounts",
                                                        ""))
        # merge_df = merge_df[~(merge_df.gene_id_fromReads.isna() & merge_df.gene_id_fromAssign.isna())]
        # below call drops reads that don't get assigned by Josh's tool
        print(f"Reads unassigned by Josh method: {merge_df[~merge_df.gene_id.isna()].shape[0]} "
              f"/{merge_df.shape[0]}, [Dropping these!]")
        merge_df = merge_df[~merge_df.gene_id.isna()]
        
        print(f"Reads with different assignments between Josh method and FeatureCounts: "
              f"{merge_df[merge_df.strand_fromReads != merge_df.strand_fromAssign].shape[0]} "
              f"/{merge_df.shape[0]}, [Dropping these!]")
        print(f"Reads with the same assignments between Josh method and FeatureCounts: "
              f"{merge_df[merge_df.strand_fromReads == merge_df.strand_fromAssign].shape[0]} "
              f"/{merge_df.shape[0]}, [Dropping these!]")
        merge_df = merge_df[merge_df.strand_fromReads == merge_df.strand_fromAssign]
        print(f"Done merging at {get_dt(for_print=True)}")
        return merge_df

    read_assignments_path = find_newest_matching_file(f"{genomeDir}/*.allChrs.parquet")
    read_assignments_df = load_read_assignments(read_assignments_path)
    print("Finished loading files!")
    # for df in [reads_df, read_assignments_df]:
    #     print(df.info())
    merged_df = merge_on_chr_pos(read_assignments_df, reads_df)
    return merged_df


def assign_with_josh_method(merged_on_reads_df: pd.DataFrame, genomeDir: str,
                            keepMultipleTranscriptInfo=False, save_it=False,
                            add_names=False) -> pd.DataFrame:
    merged_on_reads_df = adjust_5_ends(merged_on_reads_df)
    print(f"Using Josh's read assignment method w/ 5'ends!")
    if keepMultipleTranscriptInfo:
        read_assignment_path = find_newest_matching_file(f"{genomeDir}/*.allChrs.parquet")
    else:
        read_assignment_path = find_newest_matching_file(f"{genomeDir}/*.allChrsLite.parquet")
    print(f"Loading readAssignment file from {read_assignment_path}")
    read_assignment_df = pd.read_parquet(read_assignment_path)
    if keepMultipleTranscriptInfo:
        read_assignment_df = read_assignment_df.astype({"to_stop": "int64",
                                                        "to_start": "int64"})
    print(f"Loaded  readAssignment file from {read_assignment_path}")
    
    print(f"Starting merge of dataframes to make assignments")
    merged_on_reads_and_assigned_df = merged_on_reads_df.merge(read_assignment_df, on=["chr_id", "chr_pos"],
                                                               how="left", suffixes=("",
                                                                                     "_fromJoshAssign"))
    print(f"Finished merge of dataframes to make assignments")
    return merged_on_reads_and_assigned_df


def adjust_5_ends(df: pd.DataFrame):
    import re
    from tqdm import tqdm
    tqdm.pandas()

    def _flip_neg_strand_genes(chr_position: int, cigar: str, strand: str) -> int:
        if strand == "+":
            read_end = chr_position
            return read_end
        else:  # if strand == "-":
            parsed_cigar = re.findall(rf'(\d+)([MDNSIX])', cigar)
            mdn_nums = [int(num) for num, char in parsed_cigar if char in "MDN"]
            read_end = chr_position + sum(mdn_nums)
            return read_end

    if "original_chr_pos" not in df.columns.to_list():
        print(f"\nMaking adjustments for 5' ends:")
        df["original_chr_pos"] = df["chr_pos"]
        df["chr_pos"] = df.progress_apply(lambda read: _flip_neg_strand_genes(read["original_chr_pos"],
                                                                              read["cigar"],
                                                                              read["strand"]),
                                          axis=1)
    else:
        print(f"'original_chr_pos' column already found in dataframe, skipping adjustment for 5'ends!")
    return df


def gene_names_to_gene_ids(parquet_path: str = "/data16/marcus/genomes/elegansRelease100"
                                               "/Caenorhabditis_elegans.WBcel235.100.gtf"
                                               ".parquet") -> pd.DataFrame:
    df = pd.read_parquet(parquet_path)[["gene_name", "gene_id", "chr"]].drop_duplicates(ignore_index=True)
    return df


def get_dt(for_print=False, for_file=False):
    from datetime import datetime
    now = datetime.now()
    if for_print:
        return str(now.strftime("%m/%d/%y @ %I:%M:%S %p"))
    elif for_file:
        return str(now.strftime("%y%m%d"))
    else:
        return str(now.strftime("%y%m%d_%I:%M:%S%p"))


def tsv_to_parquet(tsv_path) -> str:
    parquet_path = tsv_path.rstrip("tsv") + "parquet"
    df = pd.read_csv(tsv_path, sep="\t")
    print(f"Saving new parquet to: {parquet_path}")
    df.to_parquet(parquet_path)
    return parquet_path


def find_newest_matching_file(path_str):
    # A tool that I'll want to use to grab the most recent file
    from os import path
    from glob import glob

    list_of_files = glob(path_str)
    try:
        latest_file = max(list_of_files, key=path.getctime)
        return latest_file
    except ValueError:
        raise ValueError(f"Failed to find any files matching \"{path_str}\"")


def load_ski_pelo_targets(as_df=False):
    df = pd.read_csv("/data16/marcus/working/210119_SkiPeloTargets_fromStarDust/"
                     "170723_MSandM.wtAndSkiPelo_Bounds_-12_-14_S.DESeqgeneCts_"
                     "diffExpression_2.7319418642771283e-06Down.txt", names=["gene_id"])
    if as_df:
        return df
    else:
        return df.gene_id.to_list()


def rev_compliment(seq: str, rna: bool = False) -> str:
    from Bio.Seq import Seq
    seq = Seq(seq)
    if rna:
        return seq.reverse_complement_rna()
    else:
        return seq.reverse_complement()


def parse_read_assignment_allChrs_txt(assignment_file_txt, multi_row_isoforms=False,
                                      save_file=True) -> pd.DataFrame:
    # Script to load josh's read assignment file from the *.allChrs.txt format,
    #   it can return a dataframe, but it will also (optionally) save a new parquet file!
    # Majorly rewritten on 12/14/2021, will now save a parquet with lists holding transcript/
    #   isoform information. Easily expanded with the df.explode() function!
    
    print(f"Loading read assignment file from: {assignment_file_txt}")
    df = pd.read_csv(assignment_file_txt, sep="\t",
                     names=["chr_pos", "gene_id_strand", "transcript_info"])
    print(f"\tLoaded file to dataframe.")
    
    df[["chr_id", "chr_pos"]] = df["chr_pos"].str.split("_", expand=True)
    print(f"\tParsed chr_id and chr_pos info.")
    
    df[["gene_id", "strand"]] = df["gene_id_strand"].str.split(":", expand=True)
    df.drop(columns="gene_id_strand", inplace=True)
    print(f"\tParsed gene_id and strand info.")

    # The below section is CRAZY, but allow me to 'quickly' refactor
    #   the transcript/start/stop information to lists: >>>
    df['transcript_info'] = df.transcript_info.str.split('|')  # First just split up multi-isoform spots
    print(f"\tSplit isoform info.")
    # Then split each isoform into the: trans_id, start, & stop information:
    df['transcript_info'] = df.transcript_info.apply(lambda tran_list:
                                                     np.array([tran.split(':') for tran in tran_list]).T)
    print(f"\tExtracted transcript, to_start and to_stop info.")
    # Then convert these nested listed into three columns of lists
    df[['transcript_id', 'to_start', 'to_stop']] = pd.DataFrame.from_records(df.transcript_info.to_numpy(),
                                                                             index=df.index,
                                                                             columns=['transcript_id',
                                                                                      'to_start',
                                                                                      'to_stop'])
    print(f"\tSaved transcript, to_start and to_stop info to new columns.")
    if multi_row_isoforms:
        # If we want each isoform to have it's own row, we can use explode!!
        df = df.explode(['transcript_id', 'to_start', 'to_stop'])
        additional_tag = ".isoformsAsRows"
        dtypes = {
            "chr_pos": "uint32",
            "to_start": "int32",
            "to_stop": "int32",
            "strand": "category",
            "chr_id": "category",
        }
        print(f"\tExpanded transcript, to_start and to_stop info to individual rows.")
    else:
        additional_tag = ""
        dtypes = {
            "chr_pos": "uint32",
            "to_start": "object",
            "to_stop": "object",
            "strand": "category",
            "chr_id": "category",
        }
        print(f"\tSkipping expansion of transcript, to_start and to_stop info to individual rows.")
    # <<< Maybe not that bad...? \s
    # Change datatypes:
    df = df.astype(dtypes)
    # Reorder columns:
    df = df[["chr_id", "chr_pos", "gene_id", "strand", "transcript_id", "to_start", "to_stop"]]
    print(". ", end="")
    if save_file:
        parquet_path = assignment_file_txt.rstrip(".txt") + additional_tag + ".parquet"
        df.to_parquet(parquet_path)
        print(f"Saved parquet to: {parquet_path}")
    return df


if __name__ == '__main__':
    # original_genome = '/data16/marcus/genomes/elegansRelease100/' \
    #                   'Caenorhabditis_elegans.WBcel235.100.allChrs.txt'
    # pTRI_genome = '/data16/marcus/genomes/plus-pTRIxef_elegansRelease100/' \
    #               '210928_allChrs_plus-pTRI.allChrs.txt'
    # test_df = parse_read_assignment_allChrs_txt(original_genome,
    #                                             save_file=True)
    # test_df2 = parse_read_assignment_allChrs_txt(pTRI_genome,
    #                                              save_file=True)
    # print(test_df)
    # print(test_df2)

    # assign_with_josh_method("", '/data16/marcus/genomes/elegansRelease100/',
    #                         keepMultipleTranscriptInfo=True, save_it=True, add_names=True)
    # assign_with_josh_method("", '/data16/marcus/genomes/elegansRelease100/',
    #                         keepMultipleTranscriptInfo=False, save_it=True)
    # assign_with_josh_method("", '/data16/marcus/genomes/plus-pTRIxef_elegansRelease100/',
    #                         keepMultipleTranscriptInfo=True, save_it=True, add_names=True)
    # assign_with_josh_method("", '/data16/marcus/genomes/plus-pTRIxef_elegansRelease100/',
    #                         keepMultipleTranscriptInfo=False, save_it=True)

    # sam_or_bam_class_testing()
    df = load_and_merge_lib_parquets(['polyA2', 'polyA3'],
                                     keep_transcript_info=True,
                                     add_nucleotide_fractions=True)
