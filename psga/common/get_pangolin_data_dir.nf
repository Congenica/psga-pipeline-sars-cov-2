/*
 * Obtain the full path to the pangolin-data directory from the Reference Data CSV
 * file.
 * pangolin-data is distributed as a Python package and pangolin will use it by
 * default when the data directory is not specified in the call (i.e. no "--datadir").
 * We want to control the version of pangolin-data independently of the version of
 * the pipeline, so we pass it as Reference Data.
 */
process get_pangolin_data_dir {
  input:
    path reference_data_csv

  output:
    stdout

  shell:
  '''
  python ${PSGA_ROOT_PATH}/scripts/common/reference_data.py get-location\
    "!{reference_data_csv}" \
    "pangolin-data"
  '''
}
