# nano3pTailFindR.R
# Marcus Viscardi,  Oct 19, 2022

library(tailfindr)
library(arrow)

# Holy F*&#ing S%#t, this helps:
# So the vbz plugin from ONT is extreamly important for being
# able to read the compressed nanopore hdf5 files (fast5s)
Sys.setenv(HDF5_PLUGIN_PATH = "/usr/local/hdf5/lib/plugin")

# For My Nano3P Stds Run:
Nano3P_output_dir_path <- '/data16/marcus/nanoporeSoftLinks/221028_nanoporeRun_ENO2RNAStds_Nano3P/output_dir/'
# Nano3P_df <- find_tails(fast5_dir = paste0(Nano3P_output_dir_path, 'fastqs/workspace'),
#                         save_directory = paste0(Nano3P_output_dir_path, 'tailfindr'),
#                         csv_filename = paste0(format(Sys.Date(), "%y%m%d"), '_tails.csv'),
#                         num_cores = 10,
#                         )
# write_parquet(df, paste0(Nano3P_output_dir_path, 'tailfindr/', format(Sys.Date(), "%y%m%d"), '_tails.parquet'))

# For Samira's Nano3P Stds Run:
Nano3P_output_dir_path <- '/data16/marcus/nanoporeSoftLinks/230224_nanoporeRun_RNAStds_SY-TGIRT-50ng_Nano3P/output_dir/'
Nano3P_df <- find_tails(fast5_dir = paste0(Nano3P_output_dir_path, 'fastqs/workspace'),
                        save_dir = paste0(Nano3P_output_dir_path, 'tailfindr'),
                        csv_filename = paste0(format(Sys.Date(), "%y%m%d"), '_tails.csv'),
                        num_cores = 10,
                        )
write_parquet(df, paste0(Nano3P_output_dir_path, 'tailfindr/', format(Sys.Date(), "%y%m%d"), '_tails.parquet'))

# For dRNA Stds Run:
dRNA_output_dir_path <- '/data16/marcus/nanoporeSoftLinks/221112_nanoporeRun_ENO2RNAStds_dRNA/221206_output_dir/'
# dRNA_df <- find_tails(fast5_dir = paste0(dRNA_output_dir_path, 'fastqs/workspace'),
#                       save_directory = paste0(dRNA_output_dir_path, 'tailfindr'),
#                       csv_filename = paste0(format(Sys.Date(), "%y%m%d"), '_tails.csv'),
#                       num_cores = 20,
#                       )
# write_parquet(df, paste0(dRNA_output_dir_path, 'tailfindr/', format(Sys.Date(), "%y%m%d"), '_tails.parquet'))

# With custom-cdna, oMV-38 as end_primer, and '' as front_primer:     Nothing in columns
# With no settings for custom-cdna it seems to work but everything is tagged "invalid"
#   This was an issue with the install not the run params, see below

# I uninstalled and reinstalled everything and that did SOMETHING, so the test 1000 reads ~worked:
#   They had a fairly low success rate, but any success at all is new... So I am going to run with
#   it for all of the reads!

# Confirmed "working" for Nano3P, trying to get dRNA going now!!
