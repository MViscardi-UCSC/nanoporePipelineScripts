"""
nanopore_pipeline_piper.py
Marcus Viscardi,    October 31, 2021

Goal of this script is just to feed the nanopore pipeline all of my scripts!
"""
from datetime import datetime
from pprint import pprint
from step0_nanopore_pipeline import meshSetsAndArgs, main, live_cmd_call
from nanoporePipelineCommon import get_dt


def per_library_run(lib_name, path_to_settings):
    print(f"##Starting run for {lib_name} @ {get_dt(for_print=True)}##")
    start = datetime.now()
    settings_dict = {"settings": path_to_settings,
                     # Steps to run within the pipeline: *=not default
                     #     (G)uppy basecalling,
                     #     *(T)rim TERA-Seq adapters,
                     #     (M)inimap,
                     #     (N)anopolish,
                     #     (F)eature counts,
                     #     (C)oncat files,
                     #     merge with (P)andas,
                     #     *use f(L)air to ID transcripts,
                     #     *map pTRI nanopore (S)tandards,
                     #     *and/or random e(X)tra steps (plotting)
                     "stepsToRun": "L",
                     "genomeDir": "/data16/marcus/genomes/elegansRelease100",
                     "callWithJoshMethod": False,
                     "regenerate": True,
                     }
    # For roach data:
    if lib_name.startswith('roach'):
        settings_dict["nestedData"] = True
        settings_dict["outputDir"] = f"/data16/marcus/working/220222_roach_analysis_revisit/{lib_name}/output_dir"
        settings_dict["dataDir"] = f"/data16/marcus/sequencing/220222_roach_data_revisit/{lib_name}"
    args = meshSetsAndArgs(skip_cli_dict=settings_dict)
    main_out = main(**args)
    end = datetime.now()
    print(f"##Finished run for {lib_name} @ {get_dt(for_print=True)}##")
    print(f"## This took a time of: {str(end - start)}")
    return str(end - start)


if __name__ == '__main__':
    roach_settings_dict = "/220222_roach_analysis_revisit/220222_roachRevisit_L3andL4_settingsFile.txt"
    lib_dict = dict(
        polyA1="/210528_NanoporeRun_0639_L3s/"
               "211028_polyA-updated_settingsFile.txt",
        # riboD="/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/"
        #       "210928and211029_riboD-and-yeastCarrier_settingsFile.txt",
        # totalRNA1="/210709_NanoporeRun_totalRNA_0639_L3/"
        #           "211029_totalRNA_0639_L3_settingsFile.txt",
        polyA2="/210719_nanoporeRun_polyA_0639_L3_replicate/"
               "211031_polyA_0639_L3_settingsFile.txt",
        totalRNA2="/210720_nanoporeRun_totalRNA_0639_L3_replicate/"
                  "211031_totalRNA_0639_L3_settingsFile.txt",
        # xrn1="/210905_nanoporeRun_totalRNA_5108_xrn-1-KD/"
        #      "211031_totalRNA_5108_xrn-1-KD_settingsFile.txt",
        # wt_5tera="/211118_nanoporeRun_totalRNA_5108_xrn-1-KD_5TERA/"
        #          "211118_totalRNA_5108_xrn-1-KD_settingsFile.txt",
        # smg6_5tera="/211210_nanoporeRun_totalRNA_2102_xrn-1-KD_5TERA/"
        #            "211210_totalRNA_2102_xrn-1-KD_settingsFile.txt",
        # totalRNA3="/220131_nanoporeRun_totalRNA_0639_L3_third/"
        #           "220131_totalRNA_0639_L3_third_settingsFile.txt",
        # polyA3="/220131_nanoporeRun_polyA_0639_L3_third/"
        #        "220131_polyA_0639_L3_third_settingsFile.txt",
        # L3_techRep1=roach_settings_dict,
        # L3_techRep2=roach_settings_dict,
        # L4_techRep1=roach_settings_dict,
        # L4_techRep2=roach_settings_dict,
    )
    working_path = "/data16/marcus/working"

    time_dict = {}
    file_to_save_times_to = f"/data16/marcus/scripts/nanoporePipelineScripts/" \
                            f"testOutputs/{get_dt(for_file=True)}_piperSuperRun.txt"
    with open(file_to_save_times_to, "w") as f:
        f.write(f"Run of libs on {get_dt(for_print=True)}\n")
        f.write(f"  Library:\tTime To Run:\n")
    for library, settings_file in lib_dict.items():
        settings_file = working_path + settings_file
        print("###########################################\n" * 10,
              library, settings_file)
        time_dict[library] = per_library_run(lib_name=library,
                                             path_to_settings=settings_file)
        with open(file_to_save_times_to, "a") as f:
            f.write(f"{library:>10}\t{str(time_dict[library])}\n")
    pprint(time_dict)
