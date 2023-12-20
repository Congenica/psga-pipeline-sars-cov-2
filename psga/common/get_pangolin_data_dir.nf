/*
 * Obtain the full path to the pangolin-data directory from the Reference Data CSV
 * file.
 * pangolin-data is distributed as a Python package and pangolin will use it by
 * default when the data directory is not specified in the call (i.e. no "--datadir").
 * We want to control the version of pangolin-data independently of the version of
 * the pipeline, so we pass it as Reference Data.
 */
def get_pangolin_data_dir(reference_data_csv) {
  getLocationCLI = """
    python ${PSGA_ROOT_PATH}/scripts/common/reference_data.py \
    get-location \
    ${reference_data_csv} \
    pangolin-data
  """
  process = getLocationCLI.execute()
  pangolin_data_dir = process.in.text
  if (!pangolin_data_dir) {
    // The call to `reference_data.py` did not produce anything in STDOUT
    throw new Exception("Error getting pangolin-data data directory: ${process.err.text}")
  }
  return pangolin_data_dir
}
