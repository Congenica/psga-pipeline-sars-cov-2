-- Deploy covid-pipeline:02-apptables to pg
-- requires: 01-appschema

BEGIN;

  -- Create tables, enums, et al for the SARS COV2 Application

  SET LOCAL search_path = sars_cov_2;

  CREATE TABLE IF NOT EXISTS "analysis_run" (
    "analysis_run_id" SERIAL PRIMARY KEY
   ,"analysis_run_name" VARCHAR
   ,"pipeline_version" VARCHAR
   ,"pangolearn_version" VARCHAR
   ,"pangolin_version" VARCHAR
   ,"pango_version" VARCHAR
  );

  COMMENT ON TABLE "analysis_run" IS
    'Analysis run of samples';
    COMMENT ON COLUMN "analysis_run"."analysis_run_id" IS
      'Primary key, serial sequence generated in the database';
    COMMENT ON COLUMN "analysis_run"."analysis_run_name" IS
      'The dir of the analysis run';
    COMMENT ON COLUMN "analysis_run"."pipeline_version" IS
      'The version of this pipeline';
    COMMENT ON COLUMN "analysis_run"."pangolearn_version"
        IS 'The version of pangoLEARN used by Pangolin';
    COMMENT ON COLUMN "analysis_run"."pangolin_version"
        IS 'The version of Pangolin';
    COMMENT ON COLUMN "analysis_run"."pango_version"
        IS 'The version of pango used by Pangolin';

  CREATE TYPE "pangolin_status" AS ENUM (
    'unknown',
    'fail'
   ,'passed_qc'
  );

  CREATE TABLE IF NOT EXISTS "sample" (
    "sample_id" SERIAL PRIMARY KEY
   ,"analysis_run_id" INTEGER
      REFERENCES "analysis_run" ("analysis_run_id")
      ON DELETE NO ACTION ON UPDATE NO ACTION DEFERRABLE INITIALLY DEFERRED
   ,"sample_name" VARCHAR
   ,"date_collected" TIMESTAMP WITHOUT TIME ZONE
   ,"pangolin_lineage" VARCHAR
   ,"gisaid_id" VARCHAR
   ,"genome_length" INTEGER
   ,"pangolin_status" "pangolin_status"
      NOT NULL DEFAULT 'unknown'
   ,"metadata_loaded" BOOLEAN
   ,"genbank_submit_id" VARCHAR
   ,"pangolin_conflict" DOUBLE PRECISION
   ,"pangolin_ambiguity_score" DOUBLE PRECISION
   ,"scorpio_call" VARCHAR
   ,"scorpio_support" DOUBLE PRECISION
   ,"scorpio_conflict" DOUBLE PRECISION
   ,"note" VARCHAR
  );

  COMMENT ON TABLE "sample" IS
    'Sample data';
    COMMENT ON COLUMN "sample"."sample_id" IS
      'Primary key, serial sequence generated in the database';
    COMMENT ON COLUMN "sample"."analysis_run_id" IS
      'Primary key, and foreign key to analysis_run table';
    COMMENT ON COLUMN "sample"."sample_name" IS
      'Lab sample identifier';
    COMMENT ON COLUMN "sample"."date_collected" IS
      'timestamp of sample collection';
    COMMENT ON COLUMN "sample"."pangolin_lineage" IS
      'Viral phylogenetic lineage';
    COMMENT ON COLUMN "sample"."gisaid_id" IS
      'GISAID identifier (GISAID is a global science initiative and primary source that provides open-access to genomic data of the novel coronavirus responsible for COVID-19)';
    COMMENT ON COLUMN "sample"."genome_length" IS
      'Number of base pair in the virus genome';
    COMMENT ON COLUMN "sample"."pangolin_status"
        IS 'Reported pangolin lineage status';
    COMMENT ON COLUMN "sample"."metadata_loaded"
        IS 'Flag indicating that the metadata was loaded';
    COMMENT ON COLUMN "sample"."genbank_submit_id"
        IS 'Unique identifier of GenBank submission session, which was used to submit a sample.';
    COMMENT ON COLUMN "sample"."pangolin_conflict"
        IS 'In the pangoLEARN decision tree model, a given sequence gets assigned to the most likely category based on known diversity. If a sequence can fit into more than one category, the conflict score will be greater than 0 and reflect the number of categories the sequence could fit into. If the conflict score is 0, this means that within the current decision tree there is only one category that the sequence could be assigned to.';
    COMMENT ON COLUMN "sample"."pangolin_ambiguity_score"
        IS 'This score is a function of the quantity of missing data in a sequence. It represents the proportion of relevant sites in a sequnece which were imputed to the reference values. A score of 1 indicates that no sites were imputed, while a score of 0 indicates that more sites were imputed than were not imputed. This score only includes sites which are used by the decision tree to classify a sequence.';
    COMMENT ON COLUMN "sample"."scorpio_call"
        IS 'If a query is assigned a constellation by scorpio this call is output in this column. The full set of constellations searched by default can be found at the constellations repository';
    COMMENT ON COLUMN "sample"."scorpio_support"
        IS 'The support score is the proportion of defining variants which have the alternative allele in the sequence.';
    COMMENT ON COLUMN "sample"."scorpio_conflict"
        IS 'The conflict score is the proportion of defining variants which have the reference allele in the sequence. Ambiguous/other non-ref/alt bases at each of the variant positions contribute only to the denominators of these scores';
    COMMENT ON COLUMN "sample"."note"
        IS 'If any conflicts from the decision tree, this field will output the alternative assignments. If the sequence failed QC this field will describe why. If the sequence met the SNP thresholds for scorpio to call a constellation, it’ll describe the exact SNP counts of Alt, Ref and Amb (Alternative, reference and ambiguous) alleles for that call';


  CREATE TABLE IF NOT EXISTS "sample_qc" (
    "sample_id" INTEGER
      REFERENCES "sample" ("sample_id")
      ON DELETE NO ACTION ON UPDATE NO ACTION DEFERRABLE INITIALLY DEFERRED
   ,"pct_n_bases" DOUBLE PRECISION
   ,"pct_covered_bases"	DOUBLE PRECISION
   ,"longest_no_n_run" INTEGER
   ,"num_aligned_reads"	DOUBLE PRECISION
   ,"qc_pass" BOOLEAN
   ,"qc_plot" BYTEA
  );

  COMMENT ON TABLE "sample_qc" IS
    'Sample Quality Control, QC data for a sample';
    COMMENT ON COLUMN "sample_qc"."sample_id" IS
      'Primary key, and foreign key to sample table';
    COMMENT ON COLUMN "sample_qc"."pct_n_bases" IS
      'percentage of N bases in sequenced genome';
    COMMENT ON COLUMN "sample_qc"."pct_covered_bases" IS
      'percentage of covered bases in sequenced genome';
    COMMENT ON COLUMN "sample_qc"."longest_no_n_run" IS
      'longest sequence run without N base';
    COMMENT ON COLUMN "sample_qc"."num_aligned_reads" IS
      'Number of reads succesfully aligned and mapped';
    COMMENT ON COLUMN "sample_qc"."qc_pass" IS
      'Has sample passed QC';
    COMMENT ON COLUMN "sample_qc"."qc_plot" IS
      'QC plot image';

COMMIT;
