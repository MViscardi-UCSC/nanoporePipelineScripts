"""
geneHeatmaps.py
Marcus Viscardi,    July 20,2021

Goal here is to make a CDF for 5'-ends of nanopore reads anchored to the
    polyA site for each gene. From there we can take the CDF info for each
    gene and make a heatmap of that, allowing the visualization of global
    5'-end trends.
Still not sure that the polyA site is the best anchoring spot...
I ended up going with stop codons as the anchoring site (much better
annotated)

Biggest issue here that wasn't a problem for short riboseq reads is that
    the long reads from nanopore will span introns, meaning I will need to
    make sure I am counting distance in mRNA space (exons) not genomic
    space (exons+introns)!!

August 4, 2021:
Realized that flair is losing a lot of the reads that I am interested that
are short and only map to the 3'UTR. This is because, for genes with more
than 1 isoform, these short reads don't have enough information to identify
them as one isoform ID or the other, so flair tosses them.
Because of this issue I am going to try to use CIGAR parsing to "walk" along
transcript & genomic space to find the distance to the stop codon.
    CIGAR parsing:
        Transcript walking -> add up Ms from CIGAR
        Genomic walking -> add up Ms and Ds/Ns from CIGAR
"""
import pandas as pd
import numpy as np
from typing import List
import regex as re

from step0_nanopore_pipeline import find_newest_matching_file, get_dt

pd.set_option("display.max_columns", None)


def gtf_to_df(gtf_path: str) -> pd.DataFrame:
    gtf_header = ["chr",
                  "source",
                  "feature",
                  "start",
                  "end",
                  "score",
                  "strand",
                  "frame",
                  "attributes"]
    gtf_attributes = ["gene_name",
                      "gene_id",
                      "gene_version",
                      "gene_source",
                      "gene_biotype",
                      "transcript_id",
                      "transcript_source",
                      "transcript_biotype",
                      "exon_number",
                      "exon_id",
                      ]
    gtf_df = pd.read_csv(gtf_path, sep="\t", comment="#", names=gtf_header)
    for attribute in gtf_attributes:
        gtf_df[attribute] = gtf_df["attributes"].str.extract(rf"""{attribute} "(.*?)";""")
    gtf_df.drop(columns=["attributes"], inplace=True)
    return gtf_df


def parse_parsed_gtf(parsed_gtf_path: str, for_isoforms=True) -> pd.DataFrame:
    def drop_dup_stops(stop_starts: List[int], stop_ends: List[int],
                       ends_not_starts) -> List[int]:
        stacked = np.c_[stop_starts, stop_ends]
        u_stacked = np.unique(stacked, axis=0)
        new_stop_starts = [x[0] for x in u_stacked]
        new_stop_ends = [x[1] for x in u_stacked]
        if not ends_not_starts:
            return new_stop_starts
        else:
            return new_stop_ends
    print("Loading parsed gtf. . .")
    df = pd.read_csv(parsed_gtf_path, sep="\t")
    if for_isoforms:
        print("Starting to process exon info from gtf. . .")
        exon_df = df[df["feature"] == "exon"].reset_index(drop=True)[["transcript_id",
                                                                      "gene_id",
                                                                      "start", "end",
                                                                      "exon_number",
                                                                      "strand",
                                                                      "gene_biotype",
                                                                      "gene_name"]]
        exon_group = exon_df.groupby("transcript_id")
    
        crush_exon_df = exon_group["start"].apply(list).to_frame(name="exon_starts")
        crush_exon_df["exon_ends"] = exon_group["end"].apply(list).to_frame(name="exon_ends")
        crush_exon_df["exon_nums"] = exon_group["exon_number"].apply(list).to_frame(name="exon_nums")
    
        for single_column in ["gene_id", "gene_biotype", "gene_name", "strand"]:
            crush_exon_df[f"{single_column}s"] = exon_group[f"{single_column}"].apply(list).to_frame(
                name=f"{single_column}s")
            crush_exon_df[f"{single_column}"] = crush_exon_df[f"{single_column}s"].apply(
                lambda x: pd.Series(x).dropna().mode()[0:1])
            crush_exon_df.drop([f"{single_column}s"], axis=1, inplace=True)
        print("Finished parsing exon information. . .")
    ###############################################
    print("Starting on stop information. . .")
    stop_df = df[df["feature"] == "stop_codon"].reset_index(drop=True)[["transcript_id",
                                                                        "gene_id",
                                                                        "gene_name",
                                                                        "start", "end",
                                                                        "exon_number",
                                                                        "strand"]]
    stop_group = stop_df.groupby(["gene_id", "gene_name", "strand"])
    crush_stop_df = stop_group["start"].apply(list).to_frame(name="stop_starts")
    crush_stop_df["stop_ends"] = stop_group["end"].apply(list).to_frame(name="stop_ends")
    print("Finished with stop information. . .")
    if for_isoforms:
        print("Merging stop info and exon info. . .")
        large_df = crush_exon_df.merge(crush_stop_df, on="transcript_id", how="inner")
        large_df = large_df.reset_index()[["transcript_id",
                                           "gene_id", "gene_name", "gene_biotype",
                                           "exon_starts", "exon_ends", "exon_nums",
                                           "stop_start", "stop_end",
                                           "strand"]]
    else:
        large_df = crush_stop_df.reset_index()
        # TODO: Drop duplicate stop locations!!!
        print("Dropping duplicate stop locations. . . ")
        large_df["new_stop_starts"] = pd.DataFrame(large_df.apply(lambda x: drop_dup_stops(x["stop_starts"],
                                                                                           x["stop_ends"],
                                                                                           ends_not_starts=False),
                                                                  axis=1),
                                                   index=large_df.index)
        large_df["new_stop_ends"] = pd.DataFrame(large_df.apply(lambda x: drop_dup_stops(x["stop_starts"],
                                                                                           x["stop_ends"],
                                                                                           ends_not_starts=True),
                                                                  axis=1),
                                                 index=large_df.index)
        large_df.drop(["stop_starts", "stop_ends"], axis=1, inplace=True)
        large_df.rename(columns={"new_stop_ends": "stop_ends",
                                 "new_stop_starts": "stop_starts"},
                        inplace=True)
        large_df["stop_count"] = large_df["stop_ends"].apply(len)
        print("Finished with duplicate stop locations. . . ")
        # print(large_df.head())
    # print(large_df.info())
    return large_df


def load_merge_df(merge_tsv):
    merge_df = pd.read_csv(merge_tsv, sep="\t")
    return merge_df


def load_flair_map(flair_map_path: str) -> pd.DataFrame:
    fl_df = pd.read_csv(flair_map_path, sep="\t", header=None)
    fl_df[1] = fl_df[1].str.split(",")
    fl_df = fl_df.explode(1).reset_index(drop=True)
    fl_df.rename(columns={0: "transcript_id",
                          1: "read_id"},
                 inplace=True)
    return fl_df


def load_compressed_on_transcripts(comp_path: str) -> pd.DataFrame:
    from ast import literal_eval
    df = pd.read_csv(comp_path, sep="\t", converters={"read_lengths": literal_eval,
                                                      "polya_lengths": literal_eval,
                                                      "genomic_starts": literal_eval,
                                                      "cigars": literal_eval})
    return df


def merge_gtf_and_transcripts_plus(path_to_parsed_gtf, path_to_compressedOnTranscripts,
                                   output_path=None) -> pd.DataFrame:
    compressed_trans_df = load_compressed_on_transcripts(path_to_compressedOnTranscripts)
    gtf_df = parse_parsed_gtf(path_to_parsed_gtf)
    mega_df = compressed_trans_df.merge(gtf_df, on=["transcript_id", "gene_id"])
    if output_path:
        print(f"Writing Parquet file to {output_path}")
        mega_df.to_parquet(output_path)
    return mega_df


def get_stop_distances_transcripts(parquet_path=None, path_to_parsed_gtf=None,
                                   path_to_compressedOnTranscripts=None,
                                   minify=False, dropSubHits=None,
                                   stranded=None) -> [pd.DataFrame, pd.DataFrame]:
    def cigar_parse_to_genomic_dist(cigar) -> int:
        from regex import findall
        # TODO: Add something to remove soft clipped reads from length!!
        for code in ["M", "D"]:
            reg_ex = rf'(\d+){code}'
            ints = list(map(int, findall(reg_ex, cigar)))
            summed = sum(ints)
            try:
                maxed = max(ints)
            except ValueError:
                print(f"Weird. . . no {code} in this CIGAR: {cigar}")
            if code == "M":
                m_sum = summed
                m_max = maxed
            elif code == "D":
                d_sum = summed
                d_max = maxed
        length_sum = m_sum + d_sum
        return length_sum

    def fix_read_locations(strand: str, genomic_starts: List[int],
                           cigars: List[str]) -> List[int]:
        if strand == "+":
            return genomic_starts
        elif strand == "-":
            gen_lengths = map(cigar_parse_to_genomic_dist, cigars)
            new_genomic_starts = [length + genomic_starts[i] for i, length in enumerate(gen_lengths)]
            return new_genomic_starts
        else:
            raise NotImplementedError(f"We need strand information as '+' or '-'!! Not: {strand}")

    def pick_spanning_exons(exon_nums: list, exon_starts: list,
                            exon_ends: list, genomic_starts: List[int],
                            strand: str) -> List[int]:
        """
        :return: A list containing the exon_num that each reads sits within
        """
        # TODO: There is currently a major error due to stranded-ness of exons
        #       For (-) strand transcripts, the exons are in reverse order
        #       but correct orientation!:
        #           {+: {starts:[ 1, 3, 5], ends:[ 2, 4, 6]}}
        #           {-: {starts:[10, 8, 6], ends:[11, 9, 7]}}
        #       In order to solve this I should have a check BEFORE EVERYTHING,
        #           which will assess what direction I should be prefix through
        #           exons in (ie. backwards for - strand genes/transcripts)
        # TODO: I have mostly solved the above, but the reads from "-" strand genes
        #       are still running into the issue of having their "starts" treated
        #       as their "ends". . . I think I have to parse CIGARs to get around
        #       this damn issue. . . . :'(
        read_exons = []
        for read_loc in genomic_starts:
            return_exon = 0
            for i, exon_num in enumerate(exon_nums):
                # This captures read_locs that are within exons
                start_pos = exon_starts[i]
                end_pos = exon_ends[i]
                exon_l_to_r = start_pos <= end_pos
                internal_to_exon = start_pos <= read_loc <= end_pos
                if internal_to_exon:
                    # This data structure will hold exon map, read loc, and the distance to that exon (if missed)
                    return_exon = (exon_num, read_loc, 0)
                elif end_pos <= read_loc <= start_pos:
                    raise NotImplementedError("I dunno whats happening here. . . "
                                              "But it seems like the annotation is backwards?")
                elif not exon_l_to_r:
                    raise NotImplementedError("I dunno whats happening here. . . "
                                              "But it seems like the annotation is backwards?")

            if return_exon != 0:  # Unsure why I don't have this up above?
                read_exons.append(return_exon)
            else:  # Now this part handles if the read end is outside of annotated exons
                #                                              (even if only slightly!!)
                #     {+: {starts:[ 1, 3, 5], ends:[ 2, 4, 6]}}
                #     {-: {starts:[10, 8, 6], ends:[11, 9, 7]}}
                if strand == "+":
                    start = 0
                    end = -1
                elif strand == "-":
                    start = -1
                    end = 0
                else:
                    raise NotImplementedError(f"We need strand information as '+' or '-'!! Not: {strand}")

                if read_loc <= exon_starts[start]:  # The first item in exon starts is the lowest for +
                    #                                  The last item in exon starts is the lowest for -
                    mis_annotated = (exon_nums[start], read_loc, exon_starts[start] - read_loc)
                elif read_loc >= exon_ends[end]:  # The last item in exon ends is the highest for +
                    #                              The first item in exon ends is the highest for -
                    mis_annotated = (exon_nums[end], read_loc, read_loc - exon_ends[end])
                else:  # These will be intronic reads that don't map to outer edges of gene or exons

                    # First I can compute the distance to each exon start
                    distances_to_starts = {e: abs(exon_starts[i] - read_loc) for (i, e) in enumerate(exon_nums)}
                    # And return the minimum distance and the hit exon
                    nearest_distance, nearest_exon = min(distances_to_starts), \
                                                     min(distances_to_starts, key=distances_to_starts.get)
                    mis_annotated = (nearest_exon, read_loc, nearest_distance)
                read_exons.append(mis_annotated)
        return read_exons

    def pick_stop_exon(exon_nums: list, exon_starts: list,
                       exon_ends: list, strand: str, stop_loc: list) -> List[int]:
        stop_loc = stop_loc[0]
        for i, exon_num in enumerate(exon_nums):
            return_exon = None
            start_pos = exon_starts[i]
            end_pos = exon_ends[i]
            internal_to_exon = start_pos <= stop_loc <= end_pos
            if internal_to_exon:
                return_exon = exon_num
                break  # The break will ensure we don't overwrite the value
                #         once we manage to get it!!
        if not return_exon:
            # This doesn't seem to be happening
            raise NotImplementedError(f"Stop codon not inside any exons?")
        stop_list = (return_exon, stop_loc)
        return stop_list

    def calc_exon_sizes(exon_nums: list, exon_starts: list,
                        exon_ends: list, stop_list: list) -> List[int]:
        stop_exon, stop_loc = stop_list
        exon_sizes = []
        for i, exon_num in enumerate(exon_nums):
            if stop_exon != exon_num:
                exon_size = exon_ends[i] - exon_starts[i]
            elif stop_exon == exon_num:
                exon_size = stop_loc - exon_starts[i]
            exon_sizes.append(exon_size)
        return exon_sizes

    def calc_end_stop_dist(exon_nums: list, exon_starts: list, read_end_exons: list,
                           exon_ends: list, exon_sizes: list, stop_list: list,
                           strand: str) -> List[int]:
        distance_list = []
        stop_exon, stop_pos = stop_list
        for (exon_num, read_loc, added_dist) in read_end_exons:
            if exon_num == stop_exon:
                to_subtract = stop_pos
            else:
                exon_index = int(exon_num - 1)
                to_subtract = exon_ends[exon_index]
            exon_index = int(exon_num - 1)
            match_exon_size = to_subtract - read_loc
            remaining_exon_sizes = exon_sizes[exon_index + 1:]
            distance = match_exon_size + sum(remaining_exon_sizes)
            # if strand == "-":  # THIS WAS NOT THE FIX. lol
            #     distance = distance * -1
            distance_list.append(distance)
        return distance_list

    def calc_cdf_per_row(stop_distances):
        N = len(stop_distances)
        x = np.arange(1, N + 1)
        cdf = np.cumsum(stop_distances)
        return x, cdf

    if parquet_path:
        mega_df = pd.read_parquet(parquet_path)
    elif path_to_compressedOnTranscripts and path_to_parsed_gtf:
        mega_df = merge_gtf_and_transcripts_plus(path_to_parsed_gtf, path_to_compressedOnTranscripts)
    else:
        raise ImportError(f"Please provide either a parquet path or the path to GTF and compressedOnTranscripts")

    if minify:
        if isinstance(minify, int):
            mega_df = mega_df.sample(minify, axis=0)
        else:
            mega_df = mega_df.sample(10, axis=0)
    if dropSubHits and isinstance(dropSubHits, int):
        mega_df = mega_df[mega_df["transcript_hits"] >= dropSubHits]
    if stranded and isinstance(stranded, str):
        mega_df = mega_df[mega_df["strand"] == stranded]
    mega_df["new_genomic_starts"] = pd.DataFrame(mega_df.apply(lambda x: fix_read_locations(x['strand'],
                                                                                            x['genomic_starts'],
                                                                                            x['cigars']),
                                                               axis=1),
                                                 index=mega_df.index)
    mega_df["read_end_exons"] = pd.DataFrame(mega_df.apply(lambda x: pick_spanning_exons(x['exon_nums'],
                                                                                         x['exon_starts'],
                                                                                         x['exon_ends'],
                                                                                         x['new_genomic_starts'],
                                                                                         x['strand']),
                                                           axis=1),
                                             index=mega_df.index)
    mega_df["stop_exons"] = pd.DataFrame(mega_df.apply(lambda x: pick_stop_exon(x['exon_nums'],
                                                                                x['exon_starts'],
                                                                                x['exon_ends'],
                                                                                x['strand'],
                                                                                x['stop_start']),
                                                       axis=1),
                                         index=mega_df.index)
    mega_df["exon_sizes"] = pd.DataFrame(mega_df.apply(lambda x: calc_exon_sizes(x['exon_nums'],
                                                                                 x['exon_starts'],
                                                                                 x['exon_ends'],
                                                                                 x['stop_exons']),
                                                       axis=1),
                                         index=mega_df.index)
    mega_df["stop_distances"] = pd.DataFrame(mega_df.apply(lambda x: calc_end_stop_dist(x['exon_nums'],
                                                                                        x['exon_starts'],
                                                                                        x['read_end_exons'],
                                                                                        x['exon_ends'],
                                                                                        x['exon_sizes'],
                                                                                        x['stop_exons'],
                                                                                        x["strand"]),
                                                           axis=1),
                                             index=mega_df.index)
    # mega_df[["cdf_x", "stop_cdf"]] = pd.DataFrame(
    #     mega_df.apply(lambda x: calc_cdf_per_row(x['stop_distances']),
    #                   axis=1).to_list(),
    #     index=mega_df.index)
    return mega_df, mega_df[["transcript_id",
                             "gene_id",
                             "gene_name",
                             "transcript_hits",
                             "stop_distances",
                             "strand",
                             # "stop_cdf", "cdf_x",
                             ]]


def get_stop_distances_genes(parquet_path=None, path_to_parsed_gtf=None,
                             path_to_compressedOnGenes=None,
                             minify=False, dropSubHits=None,
                             stranded=None, save_parquet_path=None) -> [pd.DataFrame, pd.DataFrame]:
    def find_dist(gen_start, cigar, stops) -> List[int]:
        distance_list = []
        numbers = list(map(int, re.findall(rf'(\d+)[MDNSI]', cigar)))
        cigar_chars = re.findall(rf'\d+([MDNSI])', cigar)
        MND_nums, MND_chars = [], []
        for i, cigar_char in enumerate(cigar_chars):
            if cigar_char in "MND":
                MND_chars.append(cigar_char)
                MND_nums.append(numbers[i])
                assert len(MND_chars) == len(MND_nums)

        # Account for softclips pushing the left edge of reads:
        if cigar_chars[0] == "S":
            # So this will add the softclipped nucleotides, moving the read
            # "start" to the right in chromosome space.
            gen_start += numbers[0]
        # Start handling negitive strand genes
        if strand == "-":
            # First change our start position to the other end of the read
            # by "walking" along MNS
            gen_start += sum(MND_chars)
        return distance_list
    
    def calc_stop_dist_cigar(genomic_starts, cigars, strand, stop_starts, stop_ends) -> List[int]:
        """
        Pseudocode:
            If +: find distance from genomic_start to stop codon by parsing out lengths of Ms & Ds/Ns
            (later on I should probably prefer just Ns as to avoid inaccuracy of nanopore). I'll be
            ignoring Is from the CIGAR for now (always?)
        
        :param genomic_starts: list of interagers that are the genomic position of the left of the read
        :param cigars: list of strings that are the CIGARs
        :param strand: (hopefully) as + or a -
        :param stop_starts: locations of the stop codons
        :param stop_ends: 
        :return: 
        """
        # TODO!!!: This
        #       Focus on taking the negative or positive strand genes into account
        #       Also try to figure out if you need to drop softclips or ignore them!
        dist_list = []
        if strand == "+":
            stops = stop_starts
        elif strand == "-":
            stops = stop_ends
        else:
            raise NotImplementedError(f"Please provide strand information as \"+\" or \"-\", not: \"{strand}\"")
        for i, gen_start in enumerate(genomic_starts):
            dist_list.append(find_dist(gen_start, cigars[i], stops))
        return dist_list
    if parquet_path:
        print(f"Loading parquet from: {parquet_path}")
        mega_df = pd.read_parquet(parquet_path)
    elif path_to_compressedOnGenes and path_to_parsed_gtf:
        from ast import literal_eval
        gene_df = pd.read_parquet(path_to_compressedOnGenes)
        gtf_df = parse_parsed_gtf(path_to_parsed_gtf, for_isoforms=False)
        mega_df = gene_df.merge(gtf_df, on=["gene_id"])
        if isinstance(save_parquet_path, str):
            mega_df.to_parquet(save_parquet_path)
    else:
        raise ImportError(f"Please provide either a parquet path or the path to GTF and compressedOnTranscripts")

    if minify:
        if isinstance(minify, int):
            mega_df = mega_df.sample(minify, axis=0)
        else:
            mega_df = mega_df.head()
    if dropSubHits and isinstance(dropSubHits, int):
        mega_df = mega_df[mega_df["gene_hits"] >= dropSubHits]
    if stranded and isinstance(stranded, str):
        mega_df = mega_df[mega_df["strand"] == stranded]
    mega_df.sort_values(by="stop_count", ascending=False, inplace=True)
    # print(mega_df[["gene_id", "stop_count", "stop_starts"]])
    # print(mega_df.head())
    
    # TODO: Below call:
    # mega_df["stop_distances"] = pd.DataFrame(mega_df.apply(lambda row: calc_stop_dist_cigar(row["genomic_starts"],
    #                                                                                         row["cigars"],
    #                                                                                         row["strand"],
    #                                                                                         row["stop_starts"],
    #                                                                                         row["stop_ends"]), axis=1),
    #                                          index=mega_df.index)
    
    return mega_df, mega_df[["gene_id",
                             "gene_name",
                             "read_hits",
                             # "stop_distances",  # TODO: This is the big thing to implement!!
                             "strand",
                             ]]


def plot_heatmap(smallish_df: pd.DataFrame, bounds: List[int] = (-100, 500),
                 title=None, extra_annotation=None):
    def manual_cdf(stop_distances: list, min_x: int, max_x: int) -> np.array:
        man_cdf = [0]
        stop_distances.sort(reverse=False)
        np_stop_dists = np.array(stop_distances)
        # Min and Max could also just be the bounds that get fed in here.
        #   Anything outside those bounds would just get rounded in?
        for x in range(min_x, max_x + 1):
            if x == min_x:
                hits_at_x = np.count_nonzero(np_stop_dists <= x)
            elif x == max_x:
                hits_at_x = np.count_nonzero(np_stop_dists >= x)
            else:
                hits_at_x = np.count_nonzero(np_stop_dists == x)
            sum_to_x = man_cdf[-1]
            sum_w_x = sum_to_x + hits_at_x
            man_cdf.append(sum_w_x)
        norm_sum = man_cdf[-1]
        normalized_cdf = np.array([x / norm_sum for x in man_cdf])
        # rank_index = np.where(normalized_cdf >= 0.5)[0][0]  # try to find where each hits halfway? Maybe hacky?
        # return_tuple = [normalized_cdf.tolist(), rank_index]
        # if not rank_index:
        #     breakpoint()
        # return return_tuple
        return normalized_cdf

    print("\n\n")
    # I hope I can use these max distances to act as max bounds for my heatmaps
    #   I'd introduce 0s anywhere I don't get hits
    smallish_df['max_dist'] = smallish_df["stop_distances"].apply(max)
    smallish_df['min_dist'] = smallish_df["stop_distances"].apply(min)
    total_max = smallish_df['max_dist'].max()
    total_min = smallish_df['min_dist'].min()

    print("Max distance from stop codon:", total_max)
    print("Min distance from stop codon:", total_min)
    plotter_df = pd.DataFrame()
    plotter_df["gene_name"] = smallish_df["gene_name"]
    # Because of some massive outliers, this min-max setup doesn't really work!!
    plotter_df["norm_cdf"] = smallish_df.apply(lambda x: manual_cdf(x["stop_distances"],
                                                                    bounds[0],
                                                                    bounds[1]),
                                               axis=1)
    plotter_df["halfway_index"] = pd.DataFrame(plotter_df.apply(lambda x: np.where(x["norm_cdf"] >= 0.5)[0][0],
                                                                axis=1), index=plotter_df.index)
    plotter_df["quarter_index"] = pd.DataFrame(plotter_df.apply(lambda x: np.where(x["norm_cdf"] >= 0.25)[0][0],
                                                                axis=1), index=plotter_df.index)
    plotter_df["threequart_index"] = pd.DataFrame(plotter_df.apply(lambda x: np.where(x["norm_cdf"] >= 0.75)[0][0],
                                                                   axis=1), index=plotter_df.index)
    plotter_df.sort_values(by=[
        "threequart_index",
        "quarter_index",
        "halfway_index",
    ], inplace=True)
    plotter_df["norm_cdf_lists"] = pd.DataFrame(plotter_df.apply(lambda x: x["norm_cdf"].tolist(),
                                                                 axis=1), index=plotter_df.index)
    x_axis_cdfs = plotter_df["norm_cdf"].to_list()

    x_labels = list(range(bounds[0] - 1, bounds[1] + 1))
    import plotly.express as px

    y_labels = plotter_df["gene_name"].to_list()

    fig = px.imshow(x_axis_cdfs,
                    x=x_labels,
                    y=y_labels,
                    aspect="free",
                    color_continuous_scale='RdBu_r')
    fig.add_shape(type="line",
                  x0=0, x1=0, xref='x',
                  y0=0, y1=1, yref='paper',
                  line_color='darkred')
    fig.add_annotation(x=0, y=-0.05, xref="x", yref="paper",
                       text="<b>Stop Codon</b>", showarrow=False,
                       font_size=15, )
    fig.add_annotation(x=0, y=-0.05, xref="paper", yref="paper",
                       text="<b>3' UTR</b>", showarrow=False,
                       font_size=15)
    fig.add_annotation(x=1, y=-0.05, xref="paper", yref="paper",
                       text="<b>CDS</b>", showarrow=False,
                       font_size=15)
    if isinstance(title, str):
        fig.update_layout(title=title)
    if isinstance(extra_annotation, str):
        fig.add_annotation(x=0, y=1.03, xref="paper", yref="paper",
                           text=f"<b>Library Working Directory:</b> {extra_annotation}", showarrow=False,
                           font_size=12)
    fig.update_traces(hovertemplate="<b>Gene Name: %{y}</b><br>Dist from Stop: %{x:+}nt"
                                    "<br>Cumulative Fraction of 5' Ends here: %{z:.3%}"
                                    "<extra></extra>", )
    fig.update_layout(coloraxis_colorbar=dict(tickformat=",.1%",
                                              title="CDF<br>of 5' Ends,<br>per Isoform",
                                              titleside="bottom"
                                              ),
                      hovermode="y")
    fig.show()


def flair_main(working_dir, path_to_merge, parsed_gtf):
    try:
        found_parquet_path = find_newest_matching_file(f"{path_to_merge}/*_superMerge.parquet")
    except ValueError:
        try:
            path_to_new_merge_tsv = find_newest_matching_file(f"{path_to_merge}/*_compressedOnTranscripts.tsv")
        except ValueError:
            raise NotImplementedError(f"Run the pipeline with '--stepsToRun L' on library @ {working_dir}")
        parquet_path = f"{path_to_merge}/{get_dt(for_output=True)}_superMerge.parquet"
        merge_gtf_and_transcripts_plus(parsed_gtf, path_to_new_merge_tsv,
                                       output_path=parquet_path)
        found_parquet_path = parquet_path
    minNumHits = 5  # or None
    strand_to_plot = None
    # or "-" or None
    df, smaller_df = get_stop_distances_transcripts(parquet_path=found_parquet_path,
                                                    # minify=500,
                                                    dropSubHits=minNumHits,
                                                    stranded=strand_to_plot,
                                                    )
    print("Total number of reads that passed all cutoffs:",
          smaller_df.transcript_hits.sum(),
          "\nMapping to a total number of", smaller_df.shape[0], "transcripts")
    if minNumHits and strand_to_plot:
        plot_title = f"5' End CDF Heatmap Relative to Stop Codon, transcripts on \"{strand_to_plot}\" strand " \
                     f"and each having >= {minNumHits} hits"
    elif minNumHits:
        plot_title = f"5' End CDF Heatmap Relative to Stop Codon, transcripts each having >= {minNumHits} hits"
    elif strand_to_plot:
        plot_title = f"5' End CDF Heatmap Relative to Stop Codon, transcripts on \"{strand_to_plot}\" strand"
    else:
        plot_title = None
    # Add step to filter out an only show ski/pelo targets: 
    skipelo_path = "/data16/marcus/working/210119_SkiPeloTargets_fromStarDust/" \
                   "170723_MSandM.wtAndSkiPelo_Bounds_-12_-14_S." \
                   "DESeqgeneCts_diffExpression_2.7319418642771283e-06Down.txt"
    skipelo_df = pd.read_csv(skipelo_path, names=["gene_id"])
    smaller_df = smaller_df.merge(skipelo_df, on="gene_id", how="inner")
    plot_heatmap(smaller_df.loc[:], bounds=[-750, 2000],
                 title=plot_title, extra_annotation=working_dir)
    print(smaller_df[["gene_name", "strand", "transcript_id"]])


def cigar_main(working_dir, path_to_merge, parsed_gtf):
    try:
        found_parquet_path = find_newest_matching_file(f"{path_to_merge}/*_superMerge-genes.parquet")
        df, smaller_df = get_stop_distances_genes(parquet_path=found_parquet_path)
    except ValueError:
        try:
            path_to_new_merge_tsv = find_newest_matching_file(f"{path_to_merge}/*_compressedOnGenes.parquet")
        except ValueError:
            raise NotImplementedError(f"Rerun the pipeline with '--stepsToRun P' on library @ {working_dir}")
        parquet_path = f"{path_to_merge}/{get_dt(for_output=True)}_superMerge-genes.parquet"
        print(f"Saving new 'supermerge-genes.parquet' to: {parquet_path}")
        df, smaller_df = get_stop_distances_genes(path_to_parsed_gtf=parsed_gtf,
                                                  path_to_compressedOnGenes=path_to_new_merge_tsv,
                                                  save_parquet_path=parquet_path)
    print(smaller_df)


if __name__ == '__main__':
    working_dir_dict = {"polyA": "210528_NanoporeRun_0639_L3s",  # Best (overkill) depth
                        "riboD": "210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3",  # Gross
                        "totalRNA": "210709_NanoporeRun_totalRNA_0639_L3",  # Low depth
                        "polyA2": "210719_nanoporeRun_polyA_0639_L3_replicate",  # Good depth
                        "totalRNA2": "210720_nanoporeRun_totalRNA_0639_L3_replicate",  # Good depth
                        }

    working_dir_name = working_dir_dict["totalRNA2"]

    path_to_merge_dir = f"/data16/marcus/working/{working_dir_name}/output_dir/merge_files"
    path_to_parsed_gtf = "/data16/marcus/scripts/nanoporePipelineScripts/" \
                         "Caenorhabditis_elegans.WBcel235.100.gtf.dataframe_parse.tsv"
    # flair_main(working_dir_name, path_to_merge_dir, path_to_parsed_gtf)
    cigar_main(working_dir_name, path_to_merge_dir, path_to_parsed_gtf)
