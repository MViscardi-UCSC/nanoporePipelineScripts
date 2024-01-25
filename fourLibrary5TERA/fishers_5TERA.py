"""
fishers_5TERA.py
Marcus Viscardi,    August 28, 2023

The plan is to take the code from chi2_fishers_testing.ipynb
and make it a little more streamlined and readable for posterity.

FIRST: I'll just copy everything in... then we can worry about making it better.
"""
import sys
import warnings
sys.path.insert(0, '/data16/marcus/scripts/nanoporePipelineScripts')
import nanoporePipelineCommon as npCommon

from tqdm.notebook import tqdm

import seaborn as sea
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default = "browser"

import numpy as np
import pandas as pd
import statistics as stats
pd.set_option('display.width', 200)
pd.set_option('display.max_columns', None)

CONVERSION_DICT = {"xrn-1-5tera": "oldN2",
                   "xrn-1-5tera-smg-6": "oldS6",
                   "5tera_xrn-1-KD_wt": "newN2",
                   "5tera_xrn-1-KD_smg-5": "newS5",
                   "5tera_xrn-1-KD_smg-6": "newS6",
                   "5tera_xrn-1-KD_smg-7": "newS7",
                   "5tera_xrn-1-KD_wt_rerun": "newerN2",
                   "5tera_xrn-1-KD_smg-6_rerun": "newerS6",
                   "5tera_xrn-1-KD_smg-5_rerun": "newerS5",
                   "sPM57": "sPM57",
                   "sPM58": "sPM58",
                   }
REV_CONVERSION_DICT = {val: key for key, val in CONVERSION_DICT.items()}

print(f"Imports done at {npCommon.get_dt(for_print=True)}")

regenerate = False
libs_to_load = sorted({
    'oldN2',
    'newN2',
    'newerN2',
    'oldS6',
    # 'newS6',
    'newerS6',
    # 'newS5',
    'newerS5',
    'newS7',
})

try:
    if regenerate:
        raise FileNotFoundError

    reads_df_raw_path = npCommon.find_newest_matching_file(
        f"./output_files/mega_merge_parquets/*_{'-'.join(libs_to_load)}_merged5TERA.reads_df.parquet")
    compressed_df_genes_raw_path = npCommon.find_newest_matching_file(
        f"./output_files/mega_merge_parquets/*_{'-'.join(libs_to_load)}_merged5TERA.compressed_df.parquet")
    print(f"Found preprocessed files at:\n\t{reads_df_raw_path}\nand:\n\t{compressed_df_genes_raw_path}")

    reads_df_genes_raw = pd.read_parquet(reads_df_raw_path)
    compressed_df_genes_raw = pd.read_parquet(compressed_df_genes_raw_path)
except FileNotFoundError:
    print(
        f"Could not find preprocessed files matching these libs: {'/'.join(libs_to_load)}\n"
        f"Going to create new ones from scratch! This will take longer.")
    reads_df_genes_raw, compressed_df_genes_raw = npCommon.load_and_merge_lib_parquets(
        [REV_CONVERSION_DICT[lib] for lib in libs_to_load],
        drop_sub_n=1,
        add_tail_groupings=False,
        drop_failed_polya=False,
        group_by_t5=True,
        use_josh_assignment=False)
    print(f"Saving new parquets to speed up future runs.")
    reads_df_genes_raw.to_parquet(
        f"./output_files/mega_merge_parquets/{npCommon.get_dt()}_{'-'.join(libs_to_load)}"
        f"_merged5TERA.reads_df.parquet")
    compressed_df_genes_raw.to_parquet(
        f"./output_files/mega_merge_parquets/{npCommon.get_dt()}_{'-'.join(libs_to_load)}"
        f"_merged5TERA.compressed_df.parquet")
print(f"Lib load done @ {npCommon.get_dt(for_print=True)}")

compressed_df_genes_short = compressed_df_genes_raw.copy()[
    ["lib", "chr_id", "gene_id", "gene_name", "t5", "gene_hits", "gene_rpm"]]
compressed_df_genes_short.query("gene_name == 'rpl-12'")

conversion_dict = CONVERSION_DICT
ans = [y for x, y in compressed_df_genes_short.groupby(['lib', 't5'], as_index=False)]
df_dict = {}
for i, df in enumerate(ans):
    lib = df.lib.unique()[0]
    t5 = df.t5.unique()[0]
    df = df[["chr_id", "gene_id", "gene_name", "gene_hits", "gene_rpm"]]
    df = df.rename(columns={col: f'{col}_{conversion_dict[lib]}_t5{t5}' for col in df.columns
                            if col not in ["chr_id", "gene_id", "gene_name"]})
    df_dict[(conversion_dict[lib], t5)] = df.set_index(["chr_id", "gene_id", "gene_name"])
    # print((conversion_dict[lib], t5))
    # print(df_dict[(conversion_dict[lib], t5)].query("gene_name == 'rpl-12'"))

super_df = pd.concat(df_dict.values(), axis=1, join='outer').fillna(0)

# This step will calculate total hits and the fraction adapted for each gene

filter_df_raw = pd.DataFrame()
for lib in libs_to_load:
    for rpm_or_hits in ["rpm", "hits"]:
        super_df[f"total_gene_{rpm_or_hits}_{lib}"] = super_df[f"gene_{rpm_or_hits}_{lib}_t5+"] + super_df[
            f"gene_{rpm_or_hits}_{lib}_t5-"]
    super_df[f"fraction_adapted_{lib}"] = super_df[f"gene_hits_{lib}_t5+"] / super_df[f"total_gene_hits_{lib}"]

    cols_to_carry_over = [col for col in super_df.columns if lib in col]
    filter_df_raw[cols_to_carry_over] = super_df[cols_to_carry_over]

from scipy.stats import chi2_contingency, chisquare, fisher_exact, boschloo_exact, barnard_exact
def row_chi2(row, target_lib_1, target_lib_2):
    array = np.array([[row[f"gene_hits_{target_lib_1}_t5-"], row[f"gene_hits_{target_lib_2}_t5-"]],
                      [row[f"gene_hits_{target_lib_1}_t5+"], row[f"gene_hits_{target_lib_2}_t5+"]]])
    try:
        chi2, p, deg_of_free, expected = chi2_contingency(array)
        return chi2, p
    except ValueError:
        return None, None

def row_fishers_exact(row, target_lib_1, target_lib_2, hits_or_rpm='hits', alternative='two-sided'):
    array = np.array([[row[f"gene_{hits_or_rpm}_{target_lib_1}_t5-"], row[f"gene_{hits_or_rpm}_{target_lib_2}_t5-"]],
                      [row[f"gene_{hits_or_rpm}_{target_lib_1}_t5+"], row[f"gene_{hits_or_rpm}_{target_lib_2}_t5+"]]])
    if alternative not in ['two-sided', 'greater', 'less']:
        raise KeyError(f"Please use 'two-sided', 'greater', or 'less' "
                       f"for the alternative hypothesis input for fisher's exact test!!")
    odds_ratio, p_value = fisher_exact(array, alternative=alternative)
    return odds_ratio, p_value

def row_boschloo_exact(row, target_lib_1, target_lib_2, hits_or_rpm='hits', alternative='two-sided', sampling_points=32):
    array = np.array([[row[f"gene_{hits_or_rpm}_{target_lib_1}_t5-"], row[f"gene_{hits_or_rpm}_{target_lib_2}_t5-"]],
                      [row[f"gene_{hits_or_rpm}_{target_lib_1}_t5+"], row[f"gene_{hits_or_rpm}_{target_lib_2}_t5+"]]])
    if alternative not in ['two-sided', 'greater', 'less']:
        raise KeyError(f"Please use 'two-sided', 'greater', or 'less' "
                       f"for the alternative hypothesis input for Boschloo's exact test!!")
    boschloo_result = boschloo_exact(array, alternative=alternative, n=sampling_points)
    return boschloo_result.statistic, boschloo_result.pvalue

def row_barnard_exact(row, target_lib_1, target_lib_2, hits_or_rpm='hits', alternative='two-sided', sampling_points=32):
    array = np.array([[row[f"gene_{hits_or_rpm}_{target_lib_1}_t5-"], row[f"gene_{hits_or_rpm}_{target_lib_2}_t5-"]],
                      [row[f"gene_{hits_or_rpm}_{target_lib_1}_t5+"], row[f"gene_{hits_or_rpm}_{target_lib_2}_t5+"]]])
    if alternative not in ['two-sided', 'greater', 'less']:
        raise KeyError(f"Please use 'two-sided', 'greater', or 'less' "
                       f"for the alternative hypothesis input for Barnard's exact test!!")
    barnard_result = barnard_exact(array, alternative=alternative, n=sampling_points)
    return barnard_result.statistic, barnard_result.pvalue


from pathlib import Path

filter_df = filter_df_raw.copy()

stat_test_record_file = Path(
    f"/home/marcus/Insync/mviscard@ucsc.edu/Google Drive/insync_folder/NMD_cleavage_and_deadenylation_paper/"
    f"raw_figures_from_python/{npCommon.get_dt()}_stat_test_record_file.txt")
# If we are subsampling I want to add _subsampled to the record file:

stat_test_record_file.write_text(f"{npCommon.get_dt(for_print=True)}\n")


def print_and_write(string):
    print(string)
    with open(stat_test_record_file, 'a') as f:
        f.write(string + "\n")


# p-value setpoint and the applied cutoffs will be used for a Bonferroni correction
#   Currently no genes will be dropped based on these filters/cutoffs
base_sig_cutoff = 0.05

cumulative_min_read_cutoff = 100
filter_with_rpm_or_hits = "hits"
filter_col_target = "total"  # "total" or "adapted" or "unadapted"
tests_to_run = [
    "fishers",
    # "boschloo",
    # "barnard",
    # "chi2",
]

p_value_cutoff_dict = {}

lib_comparisons_to_test = []
# lib_comparisons_to_test += list(zip(["newN2"]*3, ["newS5", "newS6", "newS7"]))
# lib_comparisons_to_test += [("newN2", "newS6")]
# lib_comparisons_to_test += [("oldN2", "oldS6")]
lib_comparisons_to_test += [("newerN2", "newerS6")]
lib_comparisons_to_test += [("newerN2", "newerS5")]
lib_comparisons_to_test += [("newN2", "newS7")]
# lib_comparisons_to_test += [("newerN2", "oldN2")]
# lib_comparisons_to_test += [("sPM57", "sPM58")]


gene_lists = {}

for libs in lib_comparisons_to_test:
    first_lib, second_lib = libs
    # Run the tests:
    with warnings.catch_warnings():
        if "chi2" in tests_to_run:
            tqdm.pandas(desc=f"Calculating Chi Squared for {first_lib} and {second_lib}")
            filter_df[[f"{first_lib}_v_{second_lib}_chi2_test_result",
                       f"{first_lib}_v_{second_lib}_chi2_p_value"]] = filter_df.progress_apply(
                lambda row: row_chi2(row, first_lib, second_lib), axis=1, result_type="expand")
        if "fishers" in tests_to_run:
            tqdm.pandas(desc=f"Calculating Fisher's exact for {first_lib} and {second_lib}")
            filter_df[[f"{first_lib}_v_{second_lib}_fishers_test_result",
                       f"{first_lib}_v_{second_lib}_fishers_p_value"]] = filter_df.progress_apply(
                lambda row: row_fishers_exact(row, first_lib, second_lib, hits_or_rpm='hits', alternative='less'),
                axis=1, result_type="expand")

        # Barnard and Boschloo take so so long to run and don't provide much of a power increase...
        if "barnard" in tests_to_run:
            tqdm.pandas(desc=f"Calculating Barnard's exact for {first_lib} and {second_lib}")
            filter_df[[f"{first_lib}_v_{second_lib}_barnard_test_result",
                       f"{first_lib}_v_{second_lib}_barnard_p_value"]] = filter_df.progress_apply(
                lambda row: row_barnard_exact(row, first_lib, second_lib, hits_or_rpm='hits', alternative='less',
                                              sampling_points=4), axis=1, result_type="expand")
        if "boschloo" in tests_to_run:
            tqdm.pandas(desc=f"Calculating Boschloo exact for {first_lib} and {second_lib}")
            filter_df[[f"{first_lib}_v_{second_lib}_boschloo_test_result",
                       f"{first_lib}_v_{second_lib}_boschloo_p_value"]] = filter_df.progress_apply(
                lambda row: row_boschloo_exact(row, first_lib, second_lib, hits_or_rpm='hits', alternative='less',
                                               sampling_points=4), axis=1, result_type="expand")
    # Make adjustments for multiple testing:
    for stat_test in tests_to_run:
        lib_cols_for_correction = []
        for lib in libs:
            filter_col_converter = {"total": f"total_gene_{filter_with_rpm_or_hits}_{lib}",
                                    "adapted": f"gene_{filter_with_rpm_or_hits}_{lib}_t5+",
                                    "unadapted": f"gene_{filter_with_rpm_or_hits}_{lib}_t5-"}
            lib_cols_for_correction.append(filter_col_converter[filter_col_target])
        cumulative_col_name = f"cumulative_{filter_col_target}_{filter_with_rpm_or_hits}_{first_lib}_{second_lib}"
        filter_df[cumulative_col_name] = filter_df[lib_cols_for_correction[0]] + filter_df[lib_cols_for_correction[1]]
        number_of_genes_passing_cutoff = filter_df[filter_df[cumulative_col_name] >= cumulative_min_read_cutoff].shape[
            0]
        passed_cutoff_df = filter_df[filter_df[cumulative_col_name] >= cumulative_min_read_cutoff]
        gene_lists[(libs, stat_test)] = passed_cutoff_df.index
        adjusted_sig_cutoff = base_sig_cutoff / number_of_genes_passing_cutoff
        p_value_cutoff_dict[(libs, stat_test)] = (adjusted_sig_cutoff, number_of_genes_passing_cutoff)

        print_and_write(
            f"There were {number_of_genes_passing_cutoff} genes that passed the cutoff of having "
            f">={cumulative_min_read_cutoff} cumulative {filter_col_target} {filter_with_rpm_or_hits} "
            f"between {first_lib} and {second_lib}\n\tA Bonferroni correction with this in mind will "
            f"expect a p value of {adjusted_sig_cutoff:.3g} for a significant {stat_test} test result.")
        print_and_write(
            f"***But I think this isn't a great way to filter, I am going to recalculate the cutoff based on the "
            f"number of genes that pass the filter for ALL LIBRARIES TESTED... This does mean the pool of libs "
            f"tested matters a lot more!***")
# Now we need to recalculate the cutoff based on the number of genes that pass the filter for ALL LIBRARIES TESTED:
# First we'll create the tested_genes_df and remove all duplicates, then we'll calculate the adjusted cutoff:
stat_test_sig_cutoffs = {}
stat_tested_genes_dfs = {}
for stat_test in tests_to_run:
    tested_genes_df = pd.DataFrame(columns=['gene_id', 'gene_name'])
    for (libs, list_stat_test), gene_list in gene_lists.items():
        if stat_test != list_stat_test:
            continue
        lib1, lib2 = libs
        gene_list_better = gene_list.to_frame().reset_index(drop=True)[["gene_id", "gene_name"]]
        gene_list_better[f'tested_for_{stat_test}_{lib1}-v-{lib2}'] = True
        tested_genes_df = pd.merge(tested_genes_df, gene_list_better, how='outer', on=['gene_id', 'gene_name'])
        tested_genes_df = tested_genes_df.fillna(False)
        tested_genes_df = tested_genes_df.drop_duplicates()
    tested_genes_df[f'tested_{stat_test}_for_all'] = tested_genes_df.all(axis='columns')
    total_genes_tested_count = tested_genes_df.shape[0]
    universally_tested_genes_count = tested_genes_df[tested_genes_df[f'tested_{stat_test}_for_all']].shape[0]
    adjusted_sig_cutoff = base_sig_cutoff / universally_tested_genes_count
    stat_test_sig_cutoffs[stat_test] = adjusted_sig_cutoff
    stat_tested_genes_dfs[stat_test] = tested_genes_df
    print_and_write(
        f"Of the {total_genes_tested_count} genes tested for {stat_test}, {universally_tested_genes_count} "
        f"were tested among all libraries.")
    # Add a column to the filter_df, that stores a Bool
    # depending on if the gene was in the tested_genes_df[f'tested_{stat_test}_for_all'] column:
    temp_filter_df = filter_df.copy().reset_index()
    temp_filter_df = pd.merge(temp_filter_df, tested_genes_df[['gene_id', 'gene_name', f'tested_{stat_test}_for_all']],
                              how='left', on=['gene_id', 'gene_name'])
    temp_filter_df[f'tested_{stat_test}_for_all'] = temp_filter_df[f'tested_{stat_test}_for_all'].fillna(False)
    filter_df = temp_filter_df.set_index(['chr_id', 'gene_id', 'gene_name'])

for stat_test in tests_to_run:
    for (first_lib, second_lib) in lib_comparisons_to_test:
        try:
            filter_df.sort_values(f"{first_lib}_v_{second_lib}_{stat_test}_p_value")
            filter_df[f"neg_log10_{first_lib}_v_{second_lib}_{stat_test}_p_value"] = -np.log10(
                filter_df[f"{first_lib}_v_{second_lib}_{stat_test}_p_value"])
            filter_df[f"{first_lib}_v_{second_lib}_{stat_test}_significant"] = \
                filter_df[f"{first_lib}_v_{second_lib}_{stat_test}_p_value"] <= stat_test_sig_cutoffs[stat_test]
            # Some nice printing!
            #####################
            print_and_write(f"For {first_lib} v {second_lib} {stat_test} test:")
            print_and_write(
                f"    {filter_df[filter_df[f'{first_lib}_v_{second_lib}_{stat_test}_significant']].shape[0]} "
                f"genes were significant at a p value of {stat_test_sig_cutoffs[stat_test]:.3g} "
                f"(Bonferroni corrected for all genes tested)")
            sig_genes = filter_df[filter_df[f'{first_lib}_v_{second_lib}_{stat_test}_significant']][
                [f'{first_lib}_v_{second_lib}_{stat_test}_p_value']].reset_index()[
                ['gene_id', 'gene_name', f'{first_lib}_v_{second_lib}_{stat_test}_p_value']].sort_values(
                f'{first_lib}_v_{second_lib}_{stat_test}_p_value')
            sig_genes_list = sig_genes.values.tolist()
            sig_genes_string = '\n'.join(
                [f"        {gene_id}    {gene_name:<10}    p-val: {p_value:.3g}" for gene_id, gene_name, p_value in
                 sig_genes_list])
            print_and_write(f"    Genes:\n{sig_genes_string}")
            #####################

        except KeyError:
            print(f"Couldn't find columns corresponding to '{stat_test}'!! Be sure spelling is correct!")
    print(
        f"Total Count of genes tested: {filter_df.shape[0]}, "
        f"with {filter_df[filter_df[f'tested_{stat_test}_for_all']].shape[0]} genes tested among all libraries.")
# Let's just roll back to a simple save, writing the entire df to a csv
df_output_path = stat_test_record_file.parent / f"{npCommon.get_dt()}_statTests_{'-'.join(tests_to_run)}_largeDF.csv"
print_and_write(f"Saving the entire df to {df_output_path}...")
filter_df.to_csv(df_output_path)
print_and_write('done.')
