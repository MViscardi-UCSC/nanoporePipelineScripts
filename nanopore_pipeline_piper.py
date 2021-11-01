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
                     # I could add other over-ride settings here if wanted!
                     }
    args = meshSetsAndArgs(skip_cli_dict=settings_dict)
    main_out = main(**args)
    # pipeline_path = "/data16/marcus/scripts/nanoporePipelineScripts/step0_nanopore_pipeline.py"
    # live_cmd_call(f"python3 {pipeline_path} {path_to_settings}")
    end = datetime.now()
    print(f"##Finished run for {lib_name} @ {get_dt(for_print=True)}##")
    print(f"## This took a time of: {str(end - start)}")
    return str(end - start)


if __name__ == '__main__':
    lib_dict = dict(polyA1="/210528_NanoporeRun_0639_L3s/"
                           "211028_polyA-updated_settingsFile.txt",
                    riboD="/210706_NanoporeRun_riboD-and-yeastCarrier_0639_L3/"
                          "210928and211029_riboD-and-yeastCarrier_settingsFile.txt",
                    totalRNA1="/210709_NanoporeRun_totalRNA_0639_L3/"
                              "211029_totalRNA_0639_L3_settingsFile.txt",
                    polyA2="/210719_nanoporeRun_polyA_0639_L3_replicate/"
                           "211031_polyA_0639_L3_settingsFile.txt",
                    totalRNA2="/210720_nanoporeRun_totalRNA_0639_L3_replicate/"
                              "211031_totalRNA_0639_L3_settingsFile.txt",
                    xrn1="/210905_nanoporeRun_totalRNA_5108_xrn-1-KD/"
                         "211031_totalRNA_5108_xrn-1-KD_settingsFile.txt")
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
