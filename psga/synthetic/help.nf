def printPathogenConfig() {

    log.info"""
        Synthetic pathogen
        ===================
        Environment variables:
        * <list pathogen-specific env vars here>

        Parameters:
        * <list pathogen-specific parameters here>
    """.stripIndent()
}

def printPathogenHelp() {
    log.info"""
        Synthetic pathogen
          Description:
            <describe the pathogen analysis here>

          Usage:
            nextflow run . --run [analysis_run] --metadata [csv] [parameters]

          Mandatory environment variables:
            <any mandatory env vars here>
                                    <description>

          Mandatory parameters:
            <any mandatory parameters here>    <description>

          Optional parameters:
            <any optional parameters here>    <description>
    """.stripIndent()
}
