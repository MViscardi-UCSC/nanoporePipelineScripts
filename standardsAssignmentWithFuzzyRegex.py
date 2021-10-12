"""
standardsAssignmentWithFuzzyRegex.py
Marcus Viscardi,    September 29, 2021

Right now this script uses fuzzy matching built into the
regex library to try and pick the best match for each 
"""
import regex

import pandas as pd

pd.set_option("display.max_columns", None)

TEST_PATH = "/data16/marcus/working/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/" \
            "output_dir/merge_files/210929_mergedOnReads.tsv"


def load_standards_from_merged_on_reads(merged_on_reads_path: str, filter=True) -> pd.DataFrame:
    df = pd.read_csv(merged_on_reads_path, sep="\t", low_memory=False)
    if filter:
        df = df[df["chr_id"] == "pTRI"]
    return df


def parse_test_UMIs(path_to_csv="/data16/marcus/scripts/nanoporePipelineScripts/standardsUMIs.csv",
                    keep_rev_comp=False) -> pd.DataFrame:
    df = pd.read_csv(path_to_csv)
    if keep_rev_comp:
        return df
    else:
        return df[["identity", "umi"]]


def fuzzy_matching_UMIs(search_item, search_inside):
    from thefuzz import fuzz, process
    return fuzz.partial_ratio(search_item, search_inside)


def grep_substitution_matching(search_item, search_inside, subs: int = 2, error_type: str = "s"):
    search_item_plus = "(" + search_item + ")" + "{" + error_type + "<=" + str(subs) + "}"
    # print(search_item)
    results = regex.findall(search_item_plus, search_inside)
    if results:
        print(f"Search: {search_item}")
        for i, result in enumerate(results):
            print(f" Hit{i + 1:>2}: {result}")
        if len(results) == 1:
            return 1
        else:
            return -1
    else:
        print("---No hit---")
        return 0


def regex_pandas_for_umi(read_id, read_seq, umis: dict, errors: int = 3,
                         print_stuff=False):
    hit_dict = {}
    error_pattern = "{" + "s" + "<=" + str(errors) + "}"
    # Get results for each UMI option:
    for umi_key, umi in umis.items():
        pattern = "(" + umi + ")" + error_pattern
        result = regex.search(pattern, read_seq)
        if result:
            hit_dict[umi_key] = result
    if not hit_dict:
        # Return fail if not UMIs match at all
        return "FAIL", "No Match", None, None, []

    span_start = {k: v.span()[0] for k, v in hit_dict.items()}
    span_end = {k: v.span()[1] for k, v in hit_dict.items()}
    fuzzy_changes = {k: v.fuzzy_changes[0] for k, v in hit_dict.items()}

    if len(hit_dict) == 1:
        # Return pass if only one UMI matched
        success_key = list(hit_dict.keys())[0]
        return "PASS", success_key, span_start[success_key], span_end[success_key], fuzzy_changes[success_key]

    # Unpack the successful maps that had multiple hits (used to figure
    #   out the best match!)
    fuzzy_subs = {k: v.fuzzy_counts[0] for k, v in hit_dict.items()}
    fuzzy_sums = {k: sum(v.fuzzy_counts) for k, v in hit_dict.items()}

    if print_stuff:
        print(f"Read ID: {read_id}",
              f"Span Start:    {span_start}",
              f"Span End:      {span_end}",
              f"Substitutions: {fuzzy_subs}",
              "", sep="\n")
    else:
        min_fuzz = min(fuzzy_sums.values())
        min_fuzz_keys = [k for k, v in fuzzy_sums.items() if v == min_fuzz]
        if len(min_fuzz_keys) == 1:
            success_key = min_fuzz_keys[0]
            return "PASS", success_key, span_start[success_key], span_end[success_key], fuzzy_changes[success_key]
        else:
            return "FAIL", "ambiguous", None, None, []


def multi_violin_seaborn(data, x="umi", y="polya_length", title: str = None):
    import seaborn as sea
    import matplotlib.pyplot as plt
    import matplotlib.transforms as transforms
    
    # Force order so that it matches w/ number of observations
    obs_names = data[x].value_counts().keys()
    
    ax = sea.violinplot(data=data, x=x, y=y,
                        order=obs_names)
    
    if isinstance(title, str):
        ax.set_title(title)

    # Custom transform so that I can place number of observations w/
    # X based on data (so they're on the columns), and Y based on
    # the paper/figure!
    custom_transform = transforms.blended_transform_factory(ax.transData,
                                                            ax.transAxes)
    
    # Below adapted from: https://www.python-graph-gallery.com/38-show-number-of-observation-on-boxplot
    # Calculate number of obs per group
    n_string_list = data[x].value_counts().values
    n_string_list = [str(x) for x in n_string_list.tolist()]
    n_string_list = ["n: " + i for i in n_string_list]
    
    # Add it to the plot
    pos = range(len(n_string_list))
    for tick, label in zip(pos, ax.get_xticklabels()):
        ax.text(pos[tick],
                0.01,
                n_string_list[tick],
                transform=custom_transform,
                horizontalalignment='center',
                size='x-small',
                color='k')
    plt.show()


def print_alignments(read_id, seq, span_start, span_end, fuzzy_changes, umi, umi_dict):
    span_start, span_end = int(span_start), int(span_end)
    match = seq[span_start:span_end]
    print(f"Read ID: {read_id}")
    print(f"UMI {umi:>3}: {umi_dict[umi]}")
    print(f"Match:   {match}")
    if fuzzy_changes:
        print("         ", end="")
        for i in range(len(match)):
            if i+span_start in fuzzy_changes:
                print("^", end="")
            else:
                print(" ", end="")
    print("\n")


def main(longer_match=False, errors: int = 3, head: int = None,
         title: str = None, print_matches=False):
    reads_df = load_standards_from_merged_on_reads(TEST_PATH)
    if isinstance(head, int):
        reads_df = reads_df.head(head)
    if longer_match:
        umi_df = parse_test_UMIs(path_to_csv="/data16/marcus/scripts/"
                                             "nanoporePipelineScripts/"
                                             "standardsUMIs_plusUpstreamSeq.csv")
    else:
        umi_df = parse_test_UMIs()
    # umi_df = umi_df[umi_df["identity"] == "30A"]
    # hit_counter = 0
    # for read_seq in reads_df["sequence"].to_list():
    #     for umi in umi_df["umi"].to_list():
    #         # hit_counter += grep_substitution_matching(umi, read_seq, subs=2, error_type="e")
    #         pattern = "(" + umi + ")" + "{s<=2}"
    #         result = re.search(pattern, read_seq)
    #         if result:
    #             print(result.fuzzy_counts)
    # print("\n\n\n", hit_counter)
    # print(reads_df.columns)
    umi_dict = dict((zip(umi_df.identity, umi_df.umi)))
    for key, val in umi_dict.items():
        print(key, val)
    reads_df[["umi_qc", "umi", "match_start", "match_end", "fuzzy_changes"]] = reads_df.apply(
        lambda x: regex_pandas_for_umi(x["read_id"],
                                       x["sequence"],
                                       umi_dict,
                                       errors=errors),
        axis=1, result_type='expand')
    multi_violin_seaborn(reads_df, x="umi", y="polya_length", title=title)
    if print_matches:
        print("")
        mapped_reads_df = reads_df[reads_df["umi_qc"] == "PASS"]
        mapped_reads_df.reset_index(drop=True, inplace=True)
        mapped_reads_df.apply(lambda x: print_alignments(x["read_id"],
                                                         x["sequence"],
                                                         x["match_start"],
                                                         x["match_end"],
                                                         x["fuzzy_changes"],
                                                         x["umi"],
                                                         umi_dict), axis=1)
        print("done")


if __name__ == '__main__':
    errors_allowed = 18
    head_subset = None  # 1000
    if isinstance(head_subset, int):
        main(longer_match=True, errors=errors_allowed, head=head_subset,
             title=f"Allowing {errors_allowed} of 's' type errors (subset={head_subset}):",
             print_matches=True)
    else:
        main(longer_match=True, errors=errors_allowed,
             title=f"Allowing {errors_allowed} of 's' type errors",
             print_matches=True)
