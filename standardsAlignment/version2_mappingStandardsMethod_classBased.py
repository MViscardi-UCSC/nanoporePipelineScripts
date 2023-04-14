"""
version2_mappingStandardsMethod_classBased.py
Marcus Viscardi,    March 27, 2023

This is really just all the code from version2_mappingStandardsMethod.py,
but I am spending the time to turn this into a class. I really think this
will help the deeper levels of variables being carried from method to method.

Plus it shouldn't take THAT long, and it will be a good trial for writing my
class_based_nanopore_pipeline.py code!

ALSO, I'm going to trim a lot of the fat from version2_mappingStandardsMethod.py,
as much of that came from the initial setup and building of the algorithm.
"""
import mappy as mp
import pandas as pd
from tqdm import tqdm
from pathlib import Path


class StandardsAlignerENO2:
    def __init__(self, path_to_standards_ref_fasta: str = None, fastx_file: str = None,
                 sam_path: str = None, mjv_compressed_df: pd.DataFrame = None,
                 threads_for_aligner: int = 10, library_type: str = "dRNA",
                 try_to_align_everything=False):
        if isinstance(path_to_standards_ref_fasta, str) and Path(path_to_standards_ref_fasta).exists():
            self.stds_ref_fasta = path_to_standards_ref_fasta
        else:
            self.stds_ref_fasta = "/data16/marcus/scripts/nanoporePipelineScripts/standardsAlignment/" \
                                  "220902_version2.0_releventSequences_wOutTails.fasta"
            print(f"Using default path to reference standards' barcodes:\n{self.stds_ref_fasta}")

        self.library_type = library_type
        # From: https://github.com/lh3/minimap2/blob/master/minimap.h
        #           0x100000 = forces stranded alignment (for dRNA)
        #           0x200000 = forces reverse strand alignment (for cDNA)
        if self.library_type == "dRNA":
            self.mappy_preset = "map-ont"
            self.extra_mappy_flag = 0x100000
        elif self.library_type == "cDNA":
            self.mappy_preset = "map-ont"
            self.extra_mappy_flag = 0x200000
        else:
            raise NotImplementedError(f"Please provide 'dRNA' or 'cDNA' for library_type, "
                                      f"you provided: {self.library_type}")

        self.aligner = mp.Aligner(self.stds_ref_fasta,
                                  preset=self.mappy_preset, k=14,
                                  extra_flags=self.extra_mappy_flag,
                                  n_threads=threads_for_aligner)
        if not self.aligner:
            # This was in the tutorial, haven't seen this yet.
            print(self.stds_ref_fasta)
            raise Exception("ERROR: failed to load/build index")
        else:
            print(f"Aligner indexes build for adapters named: {self.aligner.seq_names}")

        if isinstance(fastx_file, str):
            self.input_type = "fastx"
            self.input_path = fastx_file
        elif isinstance(sam_path, str):
            self.input_type = "sam"
            self.input_path = sam_path
        elif isinstance(mjv_compressed_df, pd.DataFrame):
            self.input_type = "merge_df"
            self.input_df = mjv_compressed_df
        else:
            raise NotImplementedError("Please provide either a path to a fastx or a sam file, or a merge_df DataFrame.")
        self.align_everything = try_to_align_everything
        self.output_df = None
        self.alignments_run = False

    def run_alignments(self, report_alignments=False):
        if self.input_type == "merge_df":
            if "cerENO2" not in self.input_df['chr_id'].unique().tolist():
                if not self.align_everything:
                    print(f"\nChromosome cerENO2 not found in dataframe, "
                          f"please run minimap2 w/ genome that has cerENO2 added!!\n\n"
                          f"No reads will end up being assigned!\n")
                    return self.input_df
                else:
                    print(f"\nChromosome cerENO2 not found in dataframe, "
                          f"please run minimap2 w/ genome that has cerENO2 added!!\n\n"
                          f"Because of the align_everything tag, we'll try to align things!\n")
            assignment_df = self._align_standards_from_merge_df()
            self.output_df = self.input_df.merge(assignment_df, on='read_id', how='left')
            self.alignments_run = True
            if report_alignments:
                print(self.output_df.assignment.value_counts().head(10))
            return self.output_df

    def _align_standards_from_merge_df(self):
        # First lets just collapse the df into the columns we care about:
        columns_to_keep = ['read_id',
                           'chr_id',
                           'sequence']
        cutdown_input_df = self.input_df[columns_to_keep].copy()
        tqdm.pandas(desc="Assigning ENO2 reads to barcodes/standards")
        cutdown_input_df['assignment'] = cutdown_input_df.progress_apply(lambda row:
                                                                         self.__per_read_mapping(row['read_id'],
                                                                                                 row['sequence'],
                                                                                                 'return_assignment',
                                                                                                 chr_id=row['chr_id']),
                                                                         axis=1)
        return cutdown_input_df[['read_id', 'assignment']]

    def __per_read_mapping(self, read_id: str, sequence: str, dict_or_assignment: str, chr_id: str = 'cerENO2'):
        per_read_dict = {  # This will be a dictionary to hold more general information about the read & mapping
            'read_id': read_id,
            'sequence': sequence,
        }
        if chr_id != 'cerENO2' and not self.align_everything:
            per_read_dict['assignment'] = 'NotAStandard'
            if dict_or_assignment == 'return_dict':
                return per_read_dict
            else:
                return per_read_dict['assignment']

        # For each sequence, create empty dicts to store info about mapping hits
        hit_objs = {}
        # Loop through each of the mapping hits generated:
        for hit in self.aligner.map(sequence):
            hit_objs[hit.ctg] = hit  # and store each one

        if hit_objs:  # First check if there are any hits at all
            for aligned_std, hit_obj in hit_objs.items():
                aligned_std_int = aligned_std[-6:-4]
                per_read_dict[f"hit_obj_{aligned_std_int}"] = hit_obj
            mlen_dict = {key[-2:]: hit_obj.mlen for
                         key, hit_obj in per_read_dict.items()
                         if key.startswith("hit_obj")}
            mlen_list = list(mlen_dict.values())
            list_of_max_mlen = [std for std in mlen_dict if mlen_dict[std] == max(mlen_list)]
            if mlen_list.count(max(mlen_list)) == 1:  # There is only one std with a maximized map length
                assignment = list_of_max_mlen[0]
            elif mlen_list.count(max(mlen_list)) == 6:  # All stds are mapped, prob b/c the read is truncated
                assignment = "Ambiguous_All"
            else:  # >1 std mapped with equal length, we can check NM for more info
                NM_dict = {key[-2:]: hit_obj.NM for
                           key, hit_obj in per_read_dict.items()
                           if key.startswith("hit_obj") and key[-2:] in list_of_max_mlen}
                NM_list = list(NM_dict.values())
                if NM_list.count(min(NM_list)) == 1:  # There is one std that has the lowest mismatch count
                    assignment = [std for std in NM_dict if NM_dict[std] == min(NM_list)][0]
                else:  # Not sure what these situations will look like! They should have been caught in elif above!
                    str_of_tied_stds = "_".join(sorted(list(NM_dict.keys())))
                    assignment = f"Ambiguous_subset_{str_of_tied_stds}"
            per_read_dict['assignment'] = assignment
            if dict_or_assignment == 'return_dict':
                return per_read_dict
            else:
                return per_read_dict['assignment']
        else:  # If we are here then nothing mapped!
            per_read_dict['assignment'] = 'Failed_To_Map'
            if dict_or_assignment == 'return_dict':
                return per_read_dict
            else:
                return per_read_dict['assignment']

    def plot_alignments(self,
                        save_dir_path=""):
        import matplotlib.pyplot as plt
        if not self.alignments_run:
            print(f"Please run alignments first!")
            return None
        self.output_df.assignment.value_counts().head(7).plot(kind='pie',
                                                              autopct=lambda x:
                                                              '{:,.0f}'.format(
                                                                  x * (self.output_df['assignment'].count()) / 100
                                                              ))
        plt.show()
        ax = self.output_df.assignment.value_counts().head(9).plot(kind='bar', logy=False,
                                                                   figsize=(7, 4))
        for p in ax.patches:
            ax.annotate(f"{p.get_height():,g}", (p.get_x() * 1.000, p.get_height() * 1.005))
        plt.tight_layout()
        if save_dir_path != "":
            save_dir_path_obj = Path(save_dir_path)
            if not save_dir_path_obj.exists() or not save_dir_path_obj.is_dir():
                print(f"Provide save directory: {save_dir_path}\neither doesn't exist or isn't a directory!")
                plt.show()
                return None
            save_path_obj = save_dir_path_obj / f"{get_dt()}_assignmentBarPlot"
            for file_type in [".svg", ".png"]:
                plt.savefig(save_path_obj.with_suffix(file_type))
        plt.show()
        return None


if __name__ == '__main__':
    from nanoporePipelineCommon import get_dt, pick_lib_return_path
    merge_df_path = pick_lib_return_path("5tera_xrn-1-KD_smg-6_rerun")
    merge_df = pd.read_parquet(merge_df_path)
    aligner_obj = StandardsAlignerENO2(mjv_compressed_df=merge_df,
                                       try_to_align_everything=False,
                                       )
    print(merge_df.chr_id.unique())
    aligner_obj.run_alignments(report_alignments=True)
    aligner_obj.plot_alignments(save_dir_path="/home/marcus/Desktop/")
